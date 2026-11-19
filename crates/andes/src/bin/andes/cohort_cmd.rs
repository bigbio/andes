//! `andes cohort`: plan the analysis of several projects' SDRFs.

use std::error::Error;
use std::path::{Path, PathBuf};

use crate::cli::{CohortCommand, CohortPlanArgs};
use crate::cohort::{self, Plan};
use crate::config::RunConfig;

pub(crate) fn run(cmd: CohortCommand) -> Result<(), Box<dyn Error>> {
    match cmd {
        CohortCommand::Plan(args) => plan(args),
    }
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
