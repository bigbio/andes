//! `andes cohort`: plan the analysis of several projects' SDRFs.

use std::error::Error;
use std::path::{Path, PathBuf};

use crate::cli::{CohortCommand, CohortPlanArgs, CohortRunArgs, CohortXlFdrArgs};
use crate::cohort::{self, Plan};
use crate::config::RunConfig;

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
/// prefix is the run name.
fn stage(file: &Path, at: &Path) -> std::io::Result<()> {
    if at.symlink_metadata().is_ok() {
        return Ok(());
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

/// Marker a run's folder gets once its search has finished, so that a rerun skips it.
const DONE: &str = "search.done";

fn run_cohort(args: CohortRunArgs) -> Result<(), Box<dyn Error>> {
    let config = std::fs::canonicalize(&args.config)
        .map_err(|e| format!("cannot read {}: {e}", args.config.display()))?;
    let plan = load_plan(&config)?;
    report(&plan);
    write_plan(&plan, &plan.output_dir.join("plan.json"))?;
    let exe = std::env::current_exe()?;
    let n = plan.runs.len();
    let mut skipped = 0;
    for (i, r) in plan.runs.iter().enumerate() {
        if r.dir.join(DONE).exists() {
            skipped += 1;
            continue;
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
        std::fs::write(r.dir.join(DONE), "")?;
        log::info!(
            target: "cohort",
            "search {}/{n}: {} done in {:.0} s",
            i + 1,
            r.run,
            start.elapsed().as_secs_f64()
        );
    }
    if skipped > 0 {
        log::info!(target: "cohort", "{skipped} of {n} runs already searched; kept");
    }
    let mut not_run: Vec<&str> = Vec::new();
    for g in &plan.groups {
        for stage in &g.stages {
            match *stage {
                "xl-fdr" => {
                    std::fs::create_dir_all(&g.dir)
                        .map_err(|e| format!("cannot create {}: {e}", g.dir.display()))?;
                    log::info!(target: "cohort", "group {}: crosslink FDR over {} runs", g.name, g.inputs.len());
                    xl_fdr_tables(&g.inputs, &g.dir.join(&g.name))?;
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

fn read_csms(path: &Path, out: &mut Vec<CsmRow>) -> Result<usize, String> {
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
                "{}: no `{c}` column; the table must come from an andes --xl search that writes \
                 residue_pair and protein_pair",
                path.display()
            )
        })
    };
    let (sp, cl, li, sc) = (col("spectrum")?, col("class")?, col("link")?, col("score")?);
    let (rp, pp) = (col("residue_pair")?, col("protein_pair")?);
    let mut n = 0;
    for (k, line) in lines.enumerate() {
        let f: Vec<&str> = line.split('\t').collect();
        let get = |i: usize| f.get(i).copied().unwrap_or("");
        let Ok(score) = get(sc).parse::<f64>() else {
            continue;
        };
        let bad = |what: &str| format!("{}: row {}: bad {what}", path.display(), k + 2);
        out.push(CsmRow {
            run: run.clone(),
            spectrum: get(sp).to_string(),
            class: pair_class(get(cl)).ok_or_else(|| bad("class"))?,
            link: link_class(get(li)).ok_or_else(|| bad("link"))?,
            score,
            residue_pair: get(rp).to_string(),
            protein_pair: get(pp).to_string(),
        });
        n += 1;
    }
    Ok(n)
}

fn xl_fdr(args: CohortXlFdrArgs) -> Result<(), Box<dyn Error>> {
    xl_fdr_tables(&args.tables, &args.output)
}

/// Crosslink q-values over every CSM of `tables`, written to `<output>.xl.*`.
fn xl_fdr_tables(tables: &[PathBuf], output: &Path) -> Result<(), Box<dyn Error>> {
    let mut csms = Vec::new();
    for t in tables {
        let n = read_csms(t, &mut csms)?;
        log::info!(target: "cohort", "{}: {n} scored CSMs", t.display());
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
    log::info!(
        target: "cohort",
        "wrote {} ({} CSMs from {} runs) and the residue-pair and protein-pair tables",
        csm_path.display(),
        csms.len(),
        tables.len()
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
            "spectrum\tclass\tlink\tq_value\tscore\tresidue_pair\tprotein_pair\n\
             A.r1/s1\tTT\tinter\t0.01\t12.5\tP1:5 -- P2:9\tP1 -- P2\n\
             A.r1/s2\tTD\tintra\tNA\tNA\tP1:5 -- P1:30\tP1 -- P1\n\
             A.r1/s3\tDD\thomotypic\t0.5\t3.0\tP3:7 -- P3:7\tP3 -- P3\n",
        )
        .unwrap();
        let mut rows = Vec::new();
        assert_eq!(read_csms(&p, &mut rows).unwrap(), 2);
        assert_eq!(rows[0].run, "A.r1");
        assert_eq!(rows[0].class, xl::decoy::PairClass::TT);
        assert_eq!(rows[1].link, xl::rollup::LinkClass::Homotypic);
        let old = t.path().join("B.xl.tsv");
        std::fs::write(&old, "spectrum\tclass\tlink\tscore\n").unwrap();
        assert!(read_csms(&old, &mut rows)
            .unwrap_err()
            .contains("residue_pair"));
    }
}
