//! `andes cohort`: plan the analysis of several projects' SDRFs.

use std::error::Error;
use std::path::{Path, PathBuf};

use crate::cli::{CohortCommand, CohortPlanArgs, CohortXlFdrArgs};
use crate::cohort::{self, Plan};
use crate::config::RunConfig;

pub(crate) fn run(cmd: CohortCommand) -> Result<(), Box<dyn Error>> {
    match cmd {
        CohortCommand::Plan(args) => plan(args),
        CohortCommand::XlFdr(args) => xl_fdr(args),
    }
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
    let mut csms = Vec::new();
    for t in &args.tables {
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
        let mut s = args.output.clone().into_os_string();
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
        args.tables.len()
    );
    Ok(())
}

fn plan(args: CohortPlanArgs) -> Result<(), Box<dyn Error>> {
    let cfg = RunConfig::load(&args.config)?;
    let base = args.config.parent().unwrap_or_else(|| Path::new("."));
    let plan = match cohort::plan(&cfg, base) {
        Ok(p) => p,
        Err(problems) => {
            for w in &problems.warnings {
                log::warn!(target: "cohort", "{w}");
            }
            for e in &problems.errors {
                log::error!(target: "cohort", "{e}");
            }
            let n = problems.errors.len();
            return Err(format!(
                "{n} error{} in {}; no plan written",
                if n == 1 { "" } else { "s" },
                args.config.display()
            )
            .into());
        }
    };
    for w in &plan.warnings {
        log::warn!(target: "cohort", "{w}");
    }
    report(&plan);
    let out = args
        .output
        .unwrap_or_else(|| PathBuf::from(format!("{}.plan.json", plan.cohort)));
    let json = serde_json::to_string_pretty(&plan)?;
    std::fs::write(&out, json + "\n")
        .map_err(|e| format!("cannot write {}: {e}", out.display()))?;
    log::info!(target: "cohort", "wrote {}", out.display());
    Ok(())
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
