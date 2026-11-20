//! `andes cohort`: plan the analysis of several projects' SDRFs.

use std::error::Error;
use std::path::{Path, PathBuf};

use crate::cli::{CohortCommand, CohortPlanArgs, CohortRunArgs, CohortXlFdrArgs};
use crate::cohort::{self, Plan};
use crate::config::RunConfig;
use output::stats::thousands;

pub(crate) fn run(cmd: CohortCommand) -> Result<(), Box<dyn Error>> {
    match cmd {
        CohortCommand::Plan(args) => plan(args),
        CohortCommand::XlFdr(args) => xl_fdr(args),
        CohortCommand::Run(args) => run_cohort(args),
    }
}

/// The cohort's plan, or every problem logged and an error.
fn load_plan(config: &Path) -> Result<Plan, Box<dyn Error>> {
    let cfg = RunConfig::load(config)?;
    let base = config.parent().unwrap_or_else(|| Path::new("."));
    match cohort::plan(&cfg, base) {
        Ok(p) => {
            for w in &p.warnings {
                log::warn!(target: "cohort", "{w}");
            }
            Ok(p)
        }
        Err(problems) => {
            for w in &problems.warnings {
                log::warn!(target: "cohort", "{w}");
            }
            for e in &problems.errors {
                log::error!(target: "cohort", "{e}");
            }
            let n = problems.errors.len();
            Err(format!(
                "{n} error{} in {}; no plan written",
                if n == 1 { "" } else { "s" },
                config.display()
            )
            .into())
        }
    }
}

fn write_plan(plan: &Plan, out: &Path) -> Result<(), Box<dyn Error>> {
    if let Some(dir) = out.parent().filter(|d| !d.as_os_str().is_empty()) {
        std::fs::create_dir_all(dir)
            .map_err(|e| format!("cannot create {}: {e}", dir.display()))?;
    }
    let json = serde_json::to_string_pretty(plan)?;
    std::fs::write(out, json + "\n").map_err(|e| format!("cannot write {}: {e}", out.display()))?;
    log::info!(target: "cohort", "wrote {}", out.display());
    Ok(())
}

/// Give the spectrum file its run's name inside the run folder, so that the PIN `SpecId`
/// prefix is the run name. A link left by an earlier plan that points elsewhere is replaced.
fn stage(file: &Path, at: &Path) -> std::io::Result<()> {
    if at.symlink_metadata().is_ok() {
        if std::fs::read_link(at).is_ok_and(|t| t == file) {
            return Ok(());
        }
        std::fs::remove_file(at)?;
    }
    #[cfg(unix)]
    {
        std::os::unix::fs::symlink(file, at)
    }
    #[cfg(not(unix))]
    {
        std::fs::hard_link(file, at)
    }
}

/// Marker a run's folder gets once its search has finished, so that a rerun skips it; it
/// holds the run's [`fingerprint`].
const DONE: &str = "search.done";

/// A file's size and modification time (ns since the epoch).
fn identity(meta: &std::fs::Metadata) -> (u64, u128) {
    let mtime = meta
        .modified()
        .ok()
        .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
        .map_or(0, |d| d.as_nanos());
    (meta.len(), mtime)
}

/// Every string of a config value, depth first.
fn strings<'a>(v: &'a serde_yaml::Value, out: &mut Vec<&'a str>) {
    match v {
        serde_yaml::Value::String(s) => out.push(s),
        serde_yaml::Value::Sequence(seq) => seq.iter().for_each(|x| strings(x, out)),
        serde_yaml::Value::Mapping(m) => m.values().for_each(|x| strings(x, out)),
        _ => {}
    }
}

/// What a run's search depends on: the spectrum file, the config's sections other than
/// `cohort`, the run's search flags, and every file those name (the database, a model),
/// each file by path, size and modification time. A finished search is kept only while its
/// fingerprint is unchanged, so editing the cohort section to add a project does not search
/// the other projects again, while replacing the database at the same path does.
fn fingerprint(r: &cohort::PlannedRun, config: &Path) -> Result<String, Box<dyn Error>> {
    let meta = std::fs::metadata(&r.file)
        .map_err(|e| format!("{}: cannot read {}: {e}", r.run, r.file.display()))?;
    let (len, mtime) = identity(&meta);
    let text = std::fs::read_to_string(config)
        .map_err(|e| format!("cannot read {}: {e}", config.display()))?;
    let mut yaml: serde_yaml::Value = serde_yaml::from_str(&text)?;
    if let Some(m) = yaml.as_mapping_mut() {
        m.remove("cohort");
    }
    let base = config.parent().unwrap_or_else(|| Path::new("."));
    let mut named = Vec::new();
    strings(&yaml, &mut named);
    named.extend(r.args.iter().map(String::as_str));
    let mut inputs = std::collections::BTreeSet::new();
    for s in named {
        let p = Path::new(s);
        let p = if p.is_absolute() {
            p.to_path_buf()
        } else {
            base.join(p)
        };
        if let Some(m) = std::fs::metadata(&p).ok().filter(|m| m.is_file()) {
            let (l, t) = identity(&m);
            inputs.insert(format!("input\t{}\t{l}\t{t}", p.display()));
        }
    }
    Ok(format!(
        "file\t{}\t{len}\t{mtime}\nargs\t{}\n{}\nconfig\t{}",
        r.file.display(),
        r.args.join("\u{1f}"),
        inputs.into_iter().collect::<Vec<_>>().join("\n"),
        serde_yaml::to_string(&yaml)?
    ))
}

/// The outputs the plan's groups read from a run's folder.
fn planned_outputs<'a>(plan: &'a Plan, r: &cohort::PlannedRun) -> Vec<&'a PathBuf> {
    plan.groups
        .iter()
        .flat_map(|g| g.inputs.iter())
        .filter(|p| p.starts_with(&r.dir))
        .collect()
}

/// Each run's project, by run name.
type Projects = std::collections::BTreeMap<String, String>;

/// The run -> project map of a plan written by `andes cohort plan`.
fn projects_of_plan(path: &Path) -> Result<Projects, Box<dyn Error>> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    let v: serde_json::Value =
        serde_json::from_str(&text).map_err(|e| format!("{}: not a plan: {e}", path.display()))?;
    let runs = v
        .get("runs")
        .and_then(|r| r.as_array())
        .ok_or_else(|| format!("{}: no `runs` list; not a plan", path.display()))?;
    let mut out = Projects::new();
    for r in runs {
        match (
            r.get("run").and_then(|x| x.as_str()),
            r.get("project").and_then(|x| x.as_str()),
        ) {
            (Some(run), Some(project)) => {
                out.insert(run.to_string(), project.to_string());
            }
            _ => {
                return Err(format!("{}: a run without `run` or `project`", path.display()).into())
            }
        }
    }
    Ok(out)
}

fn run_cohort(args: CohortRunArgs) -> Result<(), Box<dyn Error>> {
    let config = std::fs::canonicalize(&args.config)
        .map_err(|e| format!("cannot read {}: {e}", args.config.display()))?;
    let plan = load_plan(&config)?;
    report(&plan);
    write_plan(&plan, &plan.output_dir.join("plan.json"))?;
    let exe = std::env::current_exe()?;
    let n = plan.runs.len();
    let mut skipped = 0;
    let mut redone = 0;
    for (i, r) in plan.runs.iter().enumerate() {
        let print = fingerprint(r, &config)?;
        let missing = planned_outputs(&plan, r).into_iter().find(|p| !p.exists());
        match std::fs::read_to_string(r.dir.join(DONE)) {
            Ok(old) if old == print && missing.is_none() => {
                skipped += 1;
                continue;
            }
            Ok(old) if old == print => {
                redone += 1;
                log::info!(
                    target: "cohort",
                    "{}: {} is missing; searching again",
                    r.run,
                    missing.map_or(String::new(), |p| p.display().to_string())
                );
                std::fs::remove_file(r.dir.join(DONE))?;
            }
            Ok(_) => {
                redone += 1;
                log::info!(
                    target: "cohort",
                    "{}: its file, config or search flags changed; searching again",
                    r.run
                );
                std::fs::remove_file(r.dir.join(DONE))?;
            }
            Err(_) => {}
        }
        std::fs::create_dir_all(&r.dir)
            .map_err(|e| format!("cannot create {}: {e}", r.dir.display()))?;
        let staged = r.dir.join(&r.stage_as);
        stage(&r.file, &staged).map_err(|e| {
            format!(
                "cannot stage {} as {}: {e}",
                r.file.display(),
                staged.display()
            )
        })?;
        let log_path = r.dir.join("andes.log");
        let log = std::fs::File::create(&log_path)
            .map_err(|e| format!("cannot write {}: {e}", log_path.display()))?;
        let mut cmd = std::process::Command::new(&exe);
        cmd.arg("--config")
            .arg(&config)
            .arg("--spectrum")
            .arg(&staged)
            .arg("--output-pin")
            .arg(r.dir.join(format!("{}.pin", r.run)));
        if let Some(t) = args.threads {
            cmd.arg("--threads").arg(t.to_string());
        }
        cmd.args(&r.args).stdout(log.try_clone()?).stderr(log);
        let start = std::time::Instant::now();
        log::info!(target: "cohort", "search {}/{n}: {}", i + 1, r.run);
        let status = cmd
            .status()
            .map_err(|e| format!("{}: cannot start {}: {e}", r.run, exe.display()))?;
        if !status.success() {
            return Err(format!(
                "{}: the search failed ({status}); see {}",
                r.run,
                log_path.display()
            )
            .into());
        }
        std::fs::write(r.dir.join(DONE), &print)?;
        log::info!(
            target: "cohort",
            "search {}/{n}: {} done in {:.0} s",
            i + 1,
            r.run,
            start.elapsed().as_secs_f64()
        );
    }
    if skipped > 0 {
        log::info!(
            target: "cohort",
            "{} of {} runs already searched with the same settings; kept",
            thousands(skipped),
            thousands(n)
        );
    }
    if redone > 0 {
        log::info!(
            target: "cohort",
            "{} of {} runs searched again after a change",
            thousands(redone),
            thousands(n)
        );
    }
    let projects: Projects = plan
        .runs
        .iter()
        .map(|r| (r.run.clone(), r.project.clone()))
        .collect();
    let mut not_run: Vec<&str> = Vec::new();
    for g in &plan.groups {
        for stage in &g.stages {
            match *stage {
                "xl-fdr" => {
                    std::fs::create_dir_all(&g.dir)
                        .map_err(|e| format!("cannot create {}: {e}", g.dir.display()))?;
                    log::info!(target: "cohort", "group {}: crosslink FDR over {} runs", g.name, g.inputs.len());
                    xl_fdr_tables(&g.inputs, &g.dir.join(&g.name), &projects)?;
                }
                other => {
                    if !not_run.contains(&other) {
                        not_run.push(other);
                    }
                }
            }
        }
    }
    if !not_run.is_empty() {
        log::warn!(
            target: "cohort",
            "`andes cohort run` does not yet run the {} stage(s); those groups' runs are searched \
             and their outputs are in {}",
            not_run.join(", "),
            plan.output_dir.join("runs").display()
        );
    }
    Ok(())
}

/// One CSM read from a run's table.
struct CsmRow {
    run: String,
    spectrum: String,
    class: xl::decoy::PairClass,
    link: xl::rollup::LinkClass,
    score: f64,
    residue_pair: String,
    protein_pair: String,
}

fn pair_class(s: &str) -> Option<xl::decoy::PairClass> {
    use xl::decoy::PairClass::*;
    [TT, TD, DD].into_iter().find(|c| c.name() == s)
}

fn link_class(s: &str) -> Option<xl::rollup::LinkClass> {
    xl::rollup::LinkClass::ALL
        .into_iter()
        .find(|c| c.name() == s)
}

/// Read one run's CSM table into `out`. Rows with the score `NA` (spectra without a CSM)
/// are counted, not read; any other score that is not a number, or a scored row without
/// its residue pair or protein pair, is an error. Returns (scored, unscored).
fn read_csms(path: &Path, out: &mut Vec<CsmRow>) -> Result<(usize, usize), String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let run = name.strip_suffix(".xl.tsv").unwrap_or(&name).to_string();
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
    let col = |c: &str| {
        header.iter().position(|h| *h == c).ok_or_else(|| {
            format!(
                "{}: no `{c}` column; the table must come from an andes --xl search of this \
                 version (with the residue_pair, protein_pair and score_source columns)",
                path.display()
            )
        })
    };
    let (sp, cl, li, sc) = (col("spectrum")?, col("class")?, col("link")?, col("score")?);
    let (rp, pp) = (col("residue_pair")?, col("protein_pair")?);
    let src = col("score_source")?;
    let (mut n, mut unscored) = (0, 0);
    for (k, line) in lines.enumerate() {
        let f: Vec<&str> = line.split('\t').collect();
        let get = |i: usize| f.get(i).copied().unwrap_or("");
        let bad = |what: &str| format!("{}: row {}: bad {what}", path.display(), k + 2);
        if get(sc) == "NA" {
            unscored += 1;
            continue;
        }
        if get(src) != "pair" {
            return Err(format!(
                "{}: row {}: score_source `{}`; the run was rescored on its own \
                 (--xl-rescore-intra/--xl-rescore-inter), so its scores are not on the scale of \
                 the other runs. Search the runs without those flags to pool them",
                path.display(),
                k + 2,
                get(src)
            ));
        }
        let score = get(sc)
            .parse::<f64>()
            .ok()
            .filter(|s| s.is_finite())
            .ok_or_else(|| bad("score"))?;
        let pair = |i: usize, what: &str| {
            Some(get(i))
                .filter(|v| !v.is_empty() && *v != "NA")
                .map(str::to_string)
                .ok_or_else(|| bad(what))
        };
        out.push(CsmRow {
            run: run.clone(),
            spectrum: get(sp).to_string(),
            class: pair_class(get(cl)).ok_or_else(|| bad("class"))?,
            link: link_class(get(li)).ok_or_else(|| bad("link"))?,
            score,
            residue_pair: pair(rp, "residue_pair")?,
            protein_pair: pair(pp, "protein_pair")?,
        });
        n += 1;
    }
    Ok((n, unscored))
}

fn xl_fdr(args: CohortXlFdrArgs) -> Result<(), Box<dyn Error>> {
    let projects = match &args.plan {
        Some(p) => projects_of_plan(p)?,
        None => Projects::new(),
    };
    xl_fdr_tables(&args.tables, &args.output, &projects)
}

/// The q cuts of the per-project error report.
const PROJECT_REPORT_CUTS: [f64; 2] = [0.01, 0.05];

/// The project of a run: from the plan when one is given, else the part of the run name
/// before its first `.` (the default run name is `<project>.<file stem>`).
fn project_of<'a>(run: &'a str, projects: &'a Projects) -> &'a str {
    projects
        .get(run)
        .map(String::as_str)
        .unwrap_or_else(|| run.split_once('.').map_or(run, |(p, _)| p))
}

/// The decoy estimate of a share's false target-target CSMs, as a fraction of its TT:
/// (TD - DD) / TT, and (TD + DD) / TT for homotypic links, whose decoys are only TD and DD.
fn share_fdr(c: &xl::fdr::ClassCounts, link: xl::rollup::LinkClass) -> f64 {
    if c.tt == 0 {
        return 0.0;
    }
    let false_tt = if link == xl::rollup::LinkClass::Homotypic {
        (c.td + c.dd) as f64
    } else {
        (c.td as f64 - c.dd as f64).max(0.0)
    };
    false_tt / c.tt as f64
}

/// Every project's share of the CSMs each pooled cut of `PROJECT_REPORT_CUTS` accepts, per link
/// class, for the projects with a scored CSM of that class.
fn project_error(
    csms: &[CsmRow],
    q: &[Option<f64>],
    projects: &Projects,
) -> Vec<output::xl_projects::XlProjectRow> {
    use std::collections::BTreeMap;
    let mut counts: BTreeMap<(&str, usize), [xl::fdr::ClassCounts; PROJECT_REPORT_CUTS.len()]> =
        BTreeMap::new();
    for (c, q) in csms.iter().zip(q) {
        let li = xl::rollup::LinkClass::ALL
            .iter()
            .position(|&l| l == c.link)
            .unwrap_or(0);
        let e = counts
            .entry((project_of(&c.run, projects), li))
            .or_default();
        for (k, cut) in PROJECT_REPORT_CUTS.iter().enumerate() {
            if q.is_some_and(|v| v <= *cut) {
                e[k].add(c.class);
            }
        }
    }
    let mut rows = Vec::new();
    for ((project, li), at) in counts {
        let link = xl::rollup::LinkClass::ALL[li];
        for (c, cut) in at.iter().zip(PROJECT_REPORT_CUTS) {
            rows.push(output::xl_projects::XlProjectRow {
                project: project.to_string(),
                link: link.name().to_string(),
                q_cut: cut,
                tt: c.tt,
                td: c.td,
                dd: c.dd,
                decoy_fdr: share_fdr(c, link),
            });
        }
    }
    rows
}

fn report_project_error(rows: &[output::xl_projects::XlProjectRow]) {
    let shown: Vec<&output::xl_projects::XlProjectRow> = rows.iter().filter(|r| r.tt > 0).collect();
    if shown.is_empty() {
        return;
    }
    let width = shown
        .iter()
        .map(|r| r.project.len())
        .max()
        .unwrap_or(7)
        .max(7);
    log::info!(
        target: "cohort",
        "each project's share of the pooled list (decoy estimate = (TD - DD) / TT; \
         homotypic (TD + DD) / TT):"
    );
    log::info!(
        target: "cohort",
        "{:<width$}  {:<9}  {:>5}  {:>7}  {:>5}  {:>5}  {:>8}",
        "project",
        "link",
        "q",
        "TT",
        "TD",
        "DD",
        "estimate"
    );
    for r in shown {
        log::info!(
            target: "cohort",
            "{:<width$}  {:<9}  {:>4}%  {:>7}  {:>5}  {:>5}  {:>7.1}%",
            r.project,
            r.link,
            r.q_cut * 100.0,
            thousands(r.tt),
            thousands(r.td),
            thousands(r.dd),
            r.decoy_fdr * 100.0
        );
    }
}

/// Crosslink q-values over every CSM of `tables`, written to `<output>.xl.*`.
fn xl_fdr_tables(
    tables: &[PathBuf],
    output: &Path,
    projects: &Projects,
) -> Result<(), Box<dyn Error>> {
    let mut csms = Vec::new();
    for t in tables {
        let (n, unscored) = read_csms(t, &mut csms)?;
        log::info!(
            target: "cohort",
            "{}: {} scored CSMs, {} spectra without a CSM",
            t.display(),
            thousands(n),
            thousands(unscored)
        );
    }
    let items: Vec<xl::rollup::Item> = csms
        .iter()
        .map(|c| xl::rollup::Item {
            key: c.residue_pair.clone(),
            score: c.score,
            class: c.class,
            link: c.link,
            n_csms: 1,
        })
        .collect();
    let ppi: Vec<String> = csms.iter().map(|c| c.protein_pair.clone()).collect();
    let rows: Vec<usize> = (0..csms.len()).collect();
    let levels = crate::xl_fdr::levels_of(items, ppi, rows);
    crate::xl_fdr::report_levels(&levels);
    let with = |ext: &str| {
        let mut s = output.to_path_buf().into_os_string();
        s.push(ext);
        PathBuf::from(s)
    };
    let csm_path = with(".xl.tsv");
    let mut w = std::io::BufWriter::new(
        std::fs::File::create(&csm_path)
            .map_err(|e| format!("cannot write {}: {e}", csm_path.display()))?,
    );
    use std::io::Write;
    writeln!(
        w,
        "run\tspectrum\tclass\tlink\tscore\tq_value\tresidue_pair\tprotein_pair"
    )?;
    let mut q = vec![None; csms.len()];
    for (&row, v) in levels.csm_rows.iter().zip(&levels.csm_q) {
        q[row] = *v;
    }
    for (c, q) in csms.iter().zip(&q) {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}",
            c.run,
            c.spectrum,
            c.class.name(),
            c.link.name(),
            c.score,
            q.map_or("NA".to_string(), |v| format!("{v:.5}")),
            c.residue_pair,
            c.protein_pair
        )?;
    }
    w.flush()?;
    crate::xl_fdr::write_levels(
        &with(".xl.residue_pairs.tsv"),
        &with(".xl.ppi.tsv"),
        &levels,
        &|row| format!("{}:{}", csms[row].run, csms[row].spectrum),
    )?;
    let shares = project_error(&csms, &q, projects);
    report_project_error(&shares);
    output::xl_projects::write_xl_projects(&with(".xl.projects.parquet"), &shares)?;
    log::info!(
        target: "cohort",
        "wrote {} ({} CSMs from {} runs), the residue-pair and protein-pair tables and the \
         per-project error table",
        csm_path.display(),
        thousands(csms.len()),
        thousands(tables.len())
    );
    Ok(())
}

fn plan(args: CohortPlanArgs) -> Result<(), Box<dyn Error>> {
    let plan = load_plan(&args.config)?;
    report(&plan);
    let out = args
        .output
        .unwrap_or_else(|| plan.output_dir.join("plan.json"));
    write_plan(&plan, &out)
}

fn report(plan: &Plan) {
    let projects: std::collections::BTreeSet<&str> =
        plan.runs.iter().map(|r| r.project.as_str()).collect();
    log::info!(
        target: "summary",
        "cohort {}: {} projects, {} runs, {} group{}",
        plan.cohort,
        projects.len(),
        plan.runs.len(),
        plan.groups.len(),
        if plan.groups.len() == 1 { "" } else { "s" }
    );
    let width = plan
        .groups
        .iter()
        .map(|g| g.name.len())
        .max()
        .unwrap_or(5)
        .max(5);
    log::info!(target: "summary", "{:<width$}  kind   runs  stages  projects", "group");
    for g in &plan.groups {
        log::info!(
            target: "summary",
            "{:<width$}  {:<5}  {:>4}  {}  {}",
            g.name,
            g.kind.name(),
            g.runs.len(),
            g.stages.join(" > "),
            g.projects.join(" ")
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn csm_rows_are_read_by_header_and_unscored_rows_skipped() {
        let t = tempfile::tempdir().unwrap();
        let p = t.path().join("A.r1.xl.tsv");
        std::fs::write(
            &p,
            "spectrum\tscore_source\tclass\tlink\tq_value\tscore\tresidue_pair\tprotein_pair\n\
             A.r1/s1\tpair\tTT\tinter\t0.01\t12.5\tP1:5 -- P2:9\tP1 -- P2\n\
             A.r1/s2\tpair\tTD\tintra\tNA\tNA\tP1:5 -- P1:30\tP1 -- P1\n\
             A.r1/s3\tpair\tDD\thomotypic\t0.5\t3.0\tP3:7 -- P3:7\tP3 -- P3\n",
        )
        .unwrap();
        let mut rows = Vec::new();
        assert_eq!(read_csms(&p, &mut rows).unwrap(), (2, 1));
        assert_eq!(rows[0].run, "A.r1");
        assert_eq!(rows[0].class, xl::decoy::PairClass::TT);
        assert_eq!(rows[1].link, xl::rollup::LinkClass::Homotypic);
        let old = t.path().join("B.xl.tsv");
        std::fs::write(&old, "spectrum\tclass\tlink\tscore\n").unwrap();
        assert!(read_csms(&old, &mut rows)
            .unwrap_err()
            .contains("residue_pair"));
    }

    #[test]
    fn a_malformed_score_or_missing_pair_is_an_error() {
        let t = tempfile::tempdir().unwrap();
        let head = "spectrum\tscore_source\tclass\tlink\tscore\tresidue_pair\tprotein_pair\n";
        for (name, row, what) in [
            (
                "a.xl.tsv",
                "s1\tpair\tTT\tinter\t12x\tP1:5 -- P2:9\tP1 -- P2",
                "row 2: bad score",
            ),
            (
                "b.xl.tsv",
                "s1\tpair\tTT\tinter\tnan\tP1:5 -- P2:9\tP1 -- P2",
                "bad score",
            ),
            (
                "c.xl.tsv",
                "s1\tpair\tTT\tinter\t12.5\t\tP1 -- P2",
                "bad residue_pair",
            ),
            (
                "e.xl.tsv",
                "s1\tpercolator\tTT\tinter\t2.1\tP1:5 -- P2:9\tP1 -- P2",
                "rescored on its own",
            ),
            (
                "d.xl.tsv",
                "s1\tpair\tTT\tinter\t12.5\tP1:5 -- P2:9",
                "bad protein_pair",
            ),
        ] {
            let p = t.path().join(name);
            std::fs::write(&p, format!("{head}{row}\n")).unwrap();
            let e = read_csms(&p, &mut Vec::new()).unwrap_err();
            assert!(e.contains(what), "{name}: {e}");
        }
    }

    #[test]
    fn staging_replaces_a_link_to_another_file() {
        let t = tempfile::tempdir().unwrap();
        let (a, b, at) = (
            t.path().join("a.raw"),
            t.path().join("b.raw"),
            t.path().join("r.raw"),
        );
        std::fs::write(&a, "a").unwrap();
        std::fs::write(&b, "b").unwrap();
        stage(&a, &at).unwrap();
        stage(&a, &at).unwrap();
        stage(&b, &at).unwrap();
        assert_eq!(std::fs::read_to_string(&at).unwrap(), "b");
    }

    #[test]
    fn the_fingerprint_follows_the_search_settings_but_not_the_cohort_section() {
        let t = tempfile::tempdir().unwrap();
        let file = t.path().join("r1.raw");
        std::fs::write(&file, "spectra").unwrap();
        let config = t.path().join("c.yaml");
        let run = |args: &[&str]| cohort::PlannedRun {
            run: "A.r1".into(),
            project: "A".into(),
            sdrf: t.path().join("a.sdrf.tsv"),
            file: file.clone(),
            dir: t.path().join("runs/A.r1"),
            stage_as: "A.r1.raw".into(),
            kind: cohort::Kind::Xl,
            group: "xl".into(),
            args: args.iter().map(|s| s.to_string()).collect(),
        };
        std::fs::write(&config, "search: {threads: 4}\ncohort: {name: one}\n").unwrap();
        let base = fingerprint(&run(&["--xl", "DSSO"]), &config).unwrap();
        std::fs::write(&config, "search: {threads: 4}\ncohort: {name: two}\n").unwrap();
        assert_eq!(fingerprint(&run(&["--xl", "DSSO"]), &config).unwrap(), base);
        assert_ne!(fingerprint(&run(&["--xl", "DSBU"]), &config).unwrap(), base);
        std::fs::write(&config, "search: {threads: 8}\ncohort: {name: two}\n").unwrap();
        assert_ne!(fingerprint(&run(&["--xl", "DSSO"]), &config).unwrap(), base);
        std::fs::write(&config, "search: {threads: 4}\ncohort: {name: two}\n").unwrap();
        std::fs::write(&file, "other spectra").unwrap();
        assert_ne!(fingerprint(&run(&["--xl", "DSSO"]), &config).unwrap(), base);

        std::fs::write(t.path().join("db.fasta"), ">P1\nPEPTIDEK\n").unwrap();
        std::fs::write(&config, "io: {database: db.fasta}\ncohort: {name: two}\n").unwrap();
        let with_db = fingerprint(&run(&["--xl", "DSSO"]), &config).unwrap();
        assert!(with_db.contains("db.fasta"), "{with_db}");
        std::fs::write(t.path().join("db.fasta"), ">P1\nPEPTIDEKR\n").unwrap();
        assert_ne!(
            fingerprint(&run(&["--xl", "DSSO"]), &config).unwrap(),
            with_db,
            "a database replaced at the same path"
        );
    }

    #[test]
    fn projects_come_from_the_plan_when_one_is_given() {
        let t = tempfile::tempdir().unwrap();
        let p = t.path().join("plan.json");
        std::fs::write(
            &p,
            r#"{"runs": [{"run": "A.v1.r1", "project": "A.v1"}, {"run": "r2", "project": "B"}]}"#,
        )
        .unwrap();
        let projects = projects_of_plan(&p).unwrap();
        assert_eq!(project_of("A.v1.r1", &projects), "A.v1");
        assert_eq!(project_of("r2", &projects), "B");
        assert_eq!(project_of("C.r3", &projects), "C");
        std::fs::write(&p, r#"{"runs": [{"run": "x"}]}"#).unwrap();
        assert!(projects_of_plan(&p).is_err());
    }

    #[test]
    fn each_project_share_is_counted_at_the_pooled_cuts() {
        use xl::decoy::PairClass::{DD, TD, TT};
        use xl::rollup::LinkClass::{Homotypic, Intra};
        let row = |run: &str, class, link| CsmRow {
            run: run.to_string(),
            spectrum: String::new(),
            class,
            link,
            score: 0.0,
            residue_pair: String::new(),
            protein_pair: String::new(),
        };
        let csms = vec![
            row("A.r1", TT, Intra),
            row("A.r1", TT, Intra),
            row("A.r2", TD, Intra),
            row("A.r2", TT, Intra),
            row("B.r1", TT, Intra),
            row("B.r1", DD, Intra),
            row("B.r1", TD, Homotypic),
            row("B.r1", TT, Homotypic),
        ];
        let q = vec![
            Some(0.001),
            Some(0.02),
            Some(0.03),
            None,
            Some(0.004),
            Some(0.04),
            Some(0.01),
            Some(0.01),
        ];
        let rows = project_error(&csms, &q, &Projects::new());
        let at = |p: &str, link: &str, cut: f64| {
            rows.iter()
                .find(|r| r.project == p && r.link == link && r.q_cut == cut)
                .map(|r| (r.tt, r.td, r.dd, r.decoy_fdr))
                .unwrap()
        };
        assert_eq!(at("A", "intra", 0.01), (1, 0, 0, 0.0));
        assert_eq!(at("A", "intra", 0.05), (2, 1, 0, 0.5));
        assert_eq!(at("B", "intra", 0.05), (1, 0, 1, 0.0));
        assert_eq!(at("B", "homotypic", 0.01), (1, 1, 0, 1.0));
        assert_eq!(rows.len(), 6);
        assert_eq!(project_of("PXD1.run.with.dots", &Projects::new()), "PXD1");
    }
}
