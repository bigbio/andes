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
fn fingerprint(
    r: &cohort::PlannedRun,
    args: &[String],
    config: &Path,
) -> Result<String, Box<dyn Error>> {
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
    named.extend(args.iter().map(String::as_str));
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
        args.join("\u{1f}"),
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

/// The group's run outputs by project, in plan order.
fn inputs_by_project(plan: &Plan, inputs: &[PathBuf]) -> Vec<(String, Vec<PathBuf>)> {
    let mut out: Vec<(String, Vec<PathBuf>)> = Vec::new();
    for p in inputs {
        let project = plan
            .runs
            .iter()
            .find(|r| p.starts_with(&r.dir))
            .map_or_else(String::new, |r| r.project.clone());
        match out.iter_mut().find(|(q, _)| *q == project) {
            Some((_, v)) => v.push(p.clone()),
            None => out.push((project, vec![p.clone()])),
        }
    }
    out
}

/// The database a run searches: its `--database` flag, else the config's.
fn database_of(r: &cohort::PlannedRun, cfg: &RunConfig, config: &Path) -> Option<PathBuf> {
    let flag = r
        .args
        .iter()
        .position(|a| a == "--database")
        .and_then(|i| r.args.get(i + 1))
        .map(PathBuf::from);
    let p = flag.or_else(|| cfg.io.database.clone())?;
    Some(if p.is_absolute() {
        p
    } else {
        config.parent().unwrap_or_else(|| Path::new(".")).join(p)
    })
}

fn run_cohort(args: CohortRunArgs) -> Result<(), Box<dyn Error>> {
    let config = std::fs::canonicalize(&args.config)
        .map_err(|e| format!("cannot read {}: {e}", args.config.display()))?;
    let plan = load_plan(&config)?;
    report(&plan);
    write_plan(&plan, &plan.output_dir.join("plan.json"))?;
    let exe = std::env::current_exe()?;
    let cfg = RunConfig::load(&config)?;
    let n = plan.runs.len();
    let mut skipped = 0;
    let mut redone = 0;
    for (i, r) in plan.runs.iter().enumerate() {
        let run_args = r.args.clone();
        let print = fingerprint(r, &run_args, &config)?;
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
        cmd.current_dir(config.parent().unwrap_or(&config))
            .arg("--config")
            .arg(&config)
            .arg("--spectrum")
            .arg(&staged)
            .arg("--output-pin")
            .arg(r.dir.join(format!("{}.pin", r.run)));
        if let Some(t) = args.threads {
            cmd.arg("--threads").arg(t.to_string());
        }
        cmd.args(&run_args).stdout(log.try_clone()?).stderr(log);
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
                    for (project, inputs) in inputs_by_project(&plan, &g.inputs) {
                        log::info!(
                            target: "cohort",
                            "group {}, project {project}: crosslink FDR over {} runs",
                            g.name,
                            inputs.len()
                        );
                        let db = plan
                            .runs
                            .iter()
                            .find(|r| r.project == project)
                            .and_then(|r| database_of(r, &cfg, &config));
                        xl_fdr_tables(
                            &inputs,
                            &g.dir.join(format!("{}.{project}", g.name)),
                            &projects,
                            db.as_deref(),
                        )?;
                    }
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
    /// The index of its table among the tables read, and its `spec_id` in that table (the
    /// `SpecId` of its row in the table's sibling PIN): the key that joins rescored PIN rows
    /// back to it. Scan numbers are not unique (one search may read several spectrum files)
    /// and not always present (MGF without scans), so they are not used for the join.
    table: usize,
    spec_id: String,
    spectrum: String,
    class: xl::decoy::PairClass,
    link: xl::rollup::LinkClass,
    score: f64,
    score_source: &'static str,
    residue_pair: String,
    protein_pair: String,
    /// The two peptides (bare, I read as L) and the proteins each maps to, for the
    /// near-copy check of inter-links; empty when the table has no such columns.
    peptides: [Vec<u8>; 2],
    proteins: [Vec<String>; 2],
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
fn read_csms(path: &Path, table: usize, out: &mut Vec<CsmRow>) -> Result<(usize, usize), String> {
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
    let spec_id = col("spec_id")?;
    let opt = |c: &str| header.iter().position(|h| *h == c);
    let sides = [
        (opt("alpha_peptide"), opt("alpha_proteins")),
        (opt("beta_peptide"), opt("beta_proteins")),
    ];
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
            table,
            spec_id: get(spec_id).to_string(),
            spectrum: get(sp).to_string(),
            class: pair_class(get(cl)).ok_or_else(|| bad("class"))?,
            link: link_class(get(li)).ok_or_else(|| bad("link"))?,
            score,
            score_source: "pair",
            residue_pair: pair(rp, "residue_pair")?,
            protein_pair: pair(pp, "protein_pair")?,
            peptides: sides
                .map(|(p, _)| p.map_or_else(Vec::new, |i| crate::xl_ambiguity::folded(get(i)))),
            proteins: sides.map(|(_, a)| {
                a.map_or_else(Vec::new, |i| {
                    get(i)
                        .split(';')
                        .filter(|x| !x.is_empty())
                        .map(str::to_string)
                        .collect()
                })
            }),
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
    xl_fdr_tables(
        &args.tables,
        &args.output,
        &projects,
        args.database.as_deref(),
    )
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

/// The sibling PIN of a run's CSM table: `<base>.xl.pin` beside `<base>.xl.tsv`.
fn pin_of(table: &Path) -> PathBuf {
    let name = table
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    table.with_file_name(format!(
        "{}.xl.pin",
        name.strip_suffix(".xl.tsv").unwrap_or(&name)
    ))
}

/// Replace the pair score of the CSMs of `classes` with the score of one Percolator model
/// trained on all CSMs of `tables` (their sibling PINs pooled, the `dropped` features left
/// out), so that their q-values use the evidence of the whole group. Does nothing, with a
/// warning, when a table has no PIN beside it or Percolator cannot train a model on the
/// pooled PIN (too few confident targets, as in a small purified-complex project).
fn rescore(
    tables: &[PathBuf],
    csms: &mut [CsmRow],
    classes: &[xl::rollup::LinkClass],
    dropped: &[&str],
    what: &str,
) -> Result<(), Box<dyn Error>> {
    use std::collections::HashMap;
    use std::io::{BufRead, BufReader, Write};
    if !csms.iter().any(|c| classes.contains(&c.link)) {
        return Ok(());
    }
    let pins: Vec<PathBuf> = tables.iter().map(|t| pin_of(t)).collect();
    if let Some(missing) = pins.iter().find(|p| !p.exists()) {
        log::warn!(
            target: "cohort",
            "{} not found; {what} keep the pair score (rescoring pools every run's PIN)",
            missing.display()
        );
        return Ok(());
    }
    let dir = tempfile::Builder::new()
        .prefix("andes-xl-rescore")
        .tempdir()?;
    let pooled = dir.path().join("rescore.pin");
    let mut out = std::io::BufWriter::new(std::fs::File::create(&pooled)?);
    let mut key_of: HashMap<String, (usize, String)> = HashMap::new();
    let mut head: Option<Vec<String>> = None;
    let mut keep: Vec<usize> = Vec::new();
    let (mut si, mut sc) = (0, 0);
    for (k, pin) in pins.iter().enumerate() {
        let mut lines = BufReader::new(
            std::fs::File::open(pin).map_err(|e| format!("cannot read {}: {e}", pin.display()))?,
        )
        .lines();
        let h: Vec<String> = lines
            .next()
            .transpose()?
            .unwrap_or_default()
            .split('\t')
            .map(str::to_string)
            .collect();
        match &head {
            None => {
                let pos = |c: &str| {
                    h.iter()
                        .position(|x| x == c)
                        .ok_or_else(|| format!("{}: no `{c}` column", pin.display()))
                };
                si = pos("SpecId")?;
                sc = pos("ScanNr")?;
                keep = (0..h.len())
                    .filter(|&i| {
                        let c = h[i].as_str();
                        !["FileName", "ExpMass", "CalcMass"].contains(&c) && !dropped.contains(&c)
                    })
                    .collect();
                let cols: Vec<&str> = keep.iter().map(|&i| h[i].as_str()).collect();
                writeln!(out, "{}", cols.join("\t"))?;
                head = Some(h);
            }
            Some(first) if *first != h => {
                return Err(format!(
                    "{}: its PIN columns differ from {}; search all runs with the same andes version",
                    pin.display(),
                    pins[0].display()
                )
                .into());
            }
            Some(_) => {}
        }
        for line in lines {
            let line = line?;
            let mut f: Vec<String> = line.split('\t').map(str::to_string).collect();
            let scan = f.get(sc).cloned().unwrap_or_default();
            let n: u64 = scan
                .parse()
                .map_err(|_| format!("{}: ScanNr `{scan}` is not a number", pin.display()))?;
            let id = format!("{k}|{}", f[si]);
            let spec_id = std::mem::replace(&mut f[si], id.clone());
            f[sc] = (k as u64 * 1_000_000_000 + n).to_string();
            key_of.insert(id, (k, spec_id));
            let cols: Vec<&str> = keep
                .iter()
                .filter_map(|&i| f.get(i).map(String::as_str))
                .collect();
            writeln!(out, "{}", cols.join("\t"))?;
        }
    }
    out.flush()?;
    drop(out);
    let backend = output::PercolatorBackend::Native {
        memory: crate::memlimit::rescore_memory(None, dir.path().to_path_buf()),
    };
    let t = std::time::Instant::now();
    let mut scored = match output::run_percolator_pooled(&backend, &pooled, &[]) {
        Ok(s) => s,
        Err(e) => {
            log::warn!(
                target: "cohort",
                "{what} keep the pair score: Percolator could not rescore them ({e})"
            );
            return Ok(());
        }
    };
    let decoy_text = std::fs::read_to_string(output::decoy_results_path(&pooled))?;
    scored.extend(output::parse_psm_results(&decoy_text).map_err(|e| format!("{e}"))?);
    let mut score: HashMap<(usize, &str), f64> = HashMap::new();
    for (id, psm) in &scored {
        if let (Some((k, spec_id)), Some(s)) = (key_of.get(id), psm.score) {
            score.insert((*k, spec_id.as_str()), s);
        }
    }
    let mut unscored = 0usize;
    let mut n = 0usize;
    for c in csms.iter_mut().filter(|c| classes.contains(&c.link)) {
        match score.get(&(c.table, c.spec_id.as_str())) {
            Some(&s) => {
                c.score = s;
                c.score_source = "percolator";
                n += 1;
            }
            None => {
                c.score = f64::NEG_INFINITY;
                c.score_source = "percolator";
                unscored += 1;
            }
        }
    }
    log::info!(
        target: "cohort",
        "{what} rescored by one Percolator model over {} runs: {} CSMs ({} without a \
         score, ranked last) [{:.1}s]",
        thousands(tables.len()),
        thousands(n),
        thousands(unscored),
        t.elapsed().as_secs_f64()
    );
    Ok(())
}

/// Proteins a peptide's near-copies are collected for, at most; a CSM lists fewer.
const NEAR_COPY_CAP: usize = 64;

/// Move to `LinkClass::InterAmbiguous` every inter-link CSM with a peptide that has a
/// near-identical copy (`xl_ambiguity`) in a protein of `db` the CSM does not list for it.
/// Returns the number moved.
fn flag_ambiguous(csms: &mut [CsmRow], db: &Path) -> Result<usize, Box<dyn Error>> {
    use std::collections::HashMap;
    let t = std::time::Instant::now();
    let mut ids: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut peptides: Vec<Vec<u8>> = Vec::new();
    for c in csms
        .iter()
        .filter(|c| c.link == xl::rollup::LinkClass::Inter)
    {
        for p in &c.peptides {
            if !p.is_empty() && !ids.contains_key(p) {
                ids.insert(p.clone(), peptides.len());
                peptides.push(p.clone());
            }
        }
    }
    let copies = crate::xl_ambiguity::near_copies(db, &peptides, NEAR_COPY_CAP)?;
    let mut n = 0;
    for c in csms
        .iter_mut()
        .filter(|c| c.link == xl::rollup::LinkClass::Inter)
    {
        let ambiguous = (0..2).any(|side| {
            ids.get(&c.peptides[side])
                .is_some_and(|&i| copies[i].iter().any(|a| !c.proteins[side].contains(a)))
        });
        if ambiguous {
            c.link = xl::rollup::LinkClass::InterAmbiguous;
            n += 1;
        }
    }
    log::info!(
        target: "cohort",
        "{} of the inter-link CSMs have a peptide with a near-identical copy in another protein \
         of {}; their FDR is estimated apart (inter_ambiguous) [{:.1}s]",
        thousands(n),
        db.display(),
        t.elapsed().as_secs_f64()
    );
    Ok(n)
}

/// Crosslink q-values over every CSM of `tables`, written to `<output>.xl.*`. With the
/// database the runs searched, inter-links whose peptide has a near-copy in another protein
/// are estimated apart and every inter-link is rescored; without it inter-links keep the
/// pair score.
fn xl_fdr_tables(
    tables: &[PathBuf],
    output: &Path,
    projects: &Projects,
    database: Option<&Path>,
) -> Result<(), Box<dyn Error>> {
    let mut csms = Vec::new();
    for (k, t) in tables.iter().enumerate() {
        let (n, unscored) = read_csms(t, k, &mut csms)?;
        log::info!(
            target: "cohort",
            "{}: {} scored CSMs, {} spectra without a CSM",
            t.display(),
            thousands(n),
            thousands(unscored)
        );
    }
    rescore(
        tables,
        &mut csms,
        &[xl::rollup::LinkClass::Intra],
        &[],
        "intra-protein links",
    )?;
    match database {
        Some(db) => {
            flag_ambiguous(&mut csms, db)?;
            rescore(
                tables,
                &mut csms,
                &[
                    xl::rollup::LinkClass::Inter,
                    xl::rollup::LinkClass::InterAmbiguous,
                ],
                &crate::xl_fdr::GAP_FEATURES,
                "inter-links",
            )?;
        }
        None => log::warn!(
            target: "cohort",
            "no database given; inter-links keep the pair score and are not checked for \
             near-copies"
        ),
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
        "run\tspectrum\tclass\tlink\tscore\tscore_source\tq_value\tresidue_pair\tprotein_pair"
    )?;
    let mut q = vec![None; csms.len()];
    for (&row, v) in levels.csm_rows.iter().zip(&levels.csm_q) {
        q[row] = *v;
    }
    for (c, q) in csms.iter().zip(&q) {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}",
            c.run,
            c.spectrum,
            c.class.name(),
            c.link.name(),
            c.score,
            c.score_source,
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
#[path = "cohort_cmd/tests/cohort_cmd.rs"]
mod tests;
