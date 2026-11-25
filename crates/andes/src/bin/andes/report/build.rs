//! Builds the run report of a search from what the run gathered (`Trace`) and its results: the
//! parameter rows from `ParamSources`, the counts at each FDR level from the rescored PSMs and the
//! protein groups, the Methods paragraph and its references.

use std::collections::{HashMap, HashSet};

use model::tolerance::{PrecursorTolerance, Tolerance};
use output::{CountRow, RunReport};

use super::methods::{methods_paragraph, references, MethodsFacts};
use super::{InputHasher, StageClock};
use crate::cli::Cli;
use crate::params::sources::ParamSources;

/// The q-value thresholds the report counts PSMs, peptides and protein groups at, besides the
/// run's own `--psm-fdr`.
const REPORT_FDR_LEVELS: [f64; 3] = [0.001, 0.01, 0.05];

/// Tolerances are printed to this many parts of their unit (three decimals).
const TOLERANCE_ROUNDING: f64 = 1000.0;

/// What a search gathers for its run report while it runs.
pub(crate) struct Trace {
    pub(crate) clock: StageClock,
    pub(crate) hasher: InputHasher,
    pub(crate) sources: ParamSources,
    pub(crate) command_line: String,
    /// `(role, model id)` of each model the run used.
    pub(crate) models: Vec<(String, String)>,
    pub(crate) decoy_strategy: search::decoy::DecoyStrategy,
    pub(crate) decoy_key: Option<String>,
}

/// The results of a search the report counts and describes.
pub(crate) struct Results<'a> {
    pub(crate) cli: &'a Cli,
    pub(crate) params: &'a search::SearchParams,
    pub(crate) fragment_tolerance: Tolerance,
    pub(crate) spectra: usize,
    pub(crate) rank1_target_psms: usize,
    pub(crate) rescored: Option<&'a HashMap<String, output::PercolatorPsm>>,
    pub(crate) protein_groups: Option<&'a output::Inference>,
}

fn tolerance_text(t: &Tolerance) -> String {
    let round = |v: f64| (v * TOLERANCE_ROUNDING).round() / TOLERANCE_ROUNDING;
    match t {
        Tolerance::Ppm(v) => format!("{} ppm", round(*v)),
        Tolerance::Da(v) => format!("{} Da", round(*v)),
    }
}

fn precursor_text(params: &search::SearchParams) -> String {
    let PrecursorTolerance { left, right } = params.precursor_tolerance;
    match &params.open_search {
        Some(o) => format!("an open-search window of {:+} to {:+} Da", o.delta_min, o.delta_max),
        None if left == right => tolerance_text(&left),
        None => format!("-{} / +{}", tolerance_text(&left), tolerance_text(&right)),
    }
}

/// How the decoys were made, as the Methods sentence says it.
pub(crate) fn decoys_text(strategy: search::decoy::DecoyStrategy, seed: u64) -> String {
    use search::decoy::DecoyStrategy as D;
    match strategy {
        D::Shuffle => format!(
            "shuffling each target peptide with its cleavage residues in place (seed {seed})"
        ),
        D::Reverse => "reversing each target peptide with its cleavage residues in place".into(),
        D::Sequon => format!(
            "shuffling each target peptide with its cleavage residues and N-glycosylation \
             sequons in place (seed {seed})"
        ),
        D::Mutate => format!("mutating each target peptide (seed {seed})"),
        D::None => "the database itself (andes generated no decoys)".into(),
    }
}

/// Spectra and rank-1 target PSMs, then for each FDR level the target PSMs, their distinct
/// peptidoforms and the target protein groups at q at most that level; the levels only when the
/// run was rescored.
fn counts(r: &Results) -> Vec<CountRow> {
    let mut rows = vec![
        CountRow { what: "spectra".into(), fdr: None, n: r.spectra },
        CountRow { what: "rank-1 target PSMs".into(), fdr: None, n: r.rank1_target_psms },
    ];
    let Some(map) = r.rescored else {
        return rows;
    };
    let mut levels: Vec<f64> = REPORT_FDR_LEVELS.to_vec();
    let own = r.cli.psm_fdr();
    if !levels.contains(&own) {
        levels.push(own);
        levels.sort_by(f64::total_cmp);
    }
    for q in levels {
        let passing: Vec<&output::PercolatorPsm> =
            map.values().filter(|p| p.q_value <= q).collect();
        let peptides: HashSet<&str> = passing.iter().map(|p| p.peptide.as_str()).collect();
        rows.push(CountRow { what: "PSMs".into(), fdr: Some(q), n: passing.len() });
        rows.push(CountRow { what: "peptides".into(), fdr: Some(q), n: peptides.len() });
        if let Some(inf) = r.protein_groups {
            let n = inf.groups.iter().filter(|g| g.q_value <= q).count();
            rows.push(CountRow { what: "protein groups".into(), fdr: Some(q), n });
        }
    }
    rows
}

fn mod_names(params: &search::SearchParams, fixed: bool) -> Vec<String> {
    params
        .aa_set
        .distinct_mods()
        .iter()
        .filter(|m| m.fixed == fixed && m.parts.is_empty())
        .map(|m| model::modification::spec::describe(m))
        .collect()
}

/// The run report of a finished search. Waits for the input hashes; the wait is the
/// `hash inputs` stage.
pub(crate) fn run_report(mut trace: Trace, r: &Results, parquet_name: &str) -> RunReport {
    trace.clock.stage("hash inputs");
    let inputs = trace.hasher.join();
    let stages = trace.clock.finish();
    let cli = r.cli;
    let version = env!("CARGO_PKG_VERSION");
    let database = cli
        .database
        .as_deref()
        .and_then(|p| p.file_name())
        .map_or_else(|| "the protein database".to_string(), |n| n.to_string_lossy().into_owned());
    let decoys = decoys_text(trace.decoy_strategy, cli.decoy_seed);
    let enzymes = std::iter::once(&r.params.enzyme)
        .chain(&r.params.extra_enzymes)
        .map(|e| e.name())
        .collect::<Vec<_>>()
        .join(" + ");
    let fixed = mod_names(r.params, true);
    let variable = mod_names(r.params, false);
    let (methods, keys) = methods_paragraph(&MethodsFacts {
        version,
        database: &database,
        enzyme: &enzymes,
        missed_cleavages: r.params.max_missed_cleavages,
        precursor_tolerance: &precursor_text(r.params),
        fragment_tolerance: &tolerance_text(&r.fragment_tolerance),
        fixed_mods: &fixed,
        variable_mods: &variable,
        decoys: &decoys,
        rescored: r.rescored.is_some(),
        psm_fdr: r.rescored.map(|_| cli.psm_fdr()),
        protein_fdr: cli.protein_fdr.filter(|_| r.protein_groups.is_some()),
    });
    RunReport {
        andes_version: version.to_string(),
        command_line: trace.command_line,
        inputs,
        params: crate::params::render::report_rows(&trace.sources),
        stages,
        counts: counts(r),
        models: trace.models,
        decoys: format!("{decoys}, prefix {}", cli.db.decoy_prefix),
        decoy_key: trace.decoy_key.unwrap_or_else(|| "none (no decoys)".to_string()),
        methods,
        references: references(&keys),
        parquet_name: parquet_name.to_string(),
    }
}
