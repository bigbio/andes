//! `andes inference`: protein groups over the rescored PSMs of one or more
//! runs, written as the QPX protein-group view.
//!
//! The database is loaded with the search's options, so protein families,
//! knockoff entries and decoy pairs are the ones the PSMs were searched
//! against. The class and knockoff q-values are applied as in a search, then
//! [`output::infer_protein_groups`] runs once over every run's PSMs: one
//! parsimony and one picked competition, so the q-values are dataset-level.

use std::collections::HashMap;
use std::path::Path;

use search::SearchIndex;

use crate::cli::InferenceArgs;
use crate::model_select::parse_enzymes;

/// A run name from a Percolator result file name.
fn run_name(path: &Path) -> String {
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    for suffix in [
        ".percolator.target.psms.txt",
        ".percolator.decoy.psms.txt",
        ".target.psms.txt",
        ".psms.txt",
        ".tsv",
        ".txt",
    ] {
        if let Some(stem) = name.strip_suffix(suffix) {
            return stem.to_string();
        }
    }
    name
}

/// Every file's rows, keyed `<file index>:<PSMId>` so SpecIds repeated across
/// per-run files stay distinct.
fn read_all(
    paths: &[std::path::PathBuf],
) -> Result<HashMap<String, output::PercolatorPsm>, String> {
    let mut all = HashMap::new();
    for (i, p) in paths.iter().enumerate() {
        let text = std::fs::read_to_string(p).map_err(|e| format!("{}: {e}", p.display()))?;
        let rows = output::parse_psm_results(&text).map_err(|e| format!("{}: {e}", p.display()))?;
        for (id, mut psm) in rows {
            let key = if paths.len() > 1 {
                format!("{i}:{id}")
            } else {
                id
            };
            psm.psm_id = key.clone();
            all.insert(key, psm);
        }
    }
    Ok(all)
}

pub(crate) fn run_inference(a: InferenceArgs) -> Result<(), Box<dyn std::error::Error>> {
    if a.psms.len() != a.decoy_psms.len() {
        return Err(format!(
            "--psms has {} files and --decoy-psms {}: give the decoy results of each target file",
            a.psms.len(),
            a.decoy_psms.len()
        )
        .into());
    }
    let runs: Vec<String> = if a.runs.is_empty() {
        a.psms.iter().map(|p| run_name(p)).collect()
    } else {
        a.runs.clone()
    };
    let (enzyme, _) = parse_enzymes(&a.enzyme)?;
    let opts = input::UniProtOptions {
        variants: a.uniprot_variants,
        isoforms: a.uniprot_isoforms,
        chains: a.uniprot_chains,
        evidence: a.uniprot_evidence,
        knockoffs: Some(input::KnockoffOptions {
            per_variant: a.knockoffs,
            enzyme,
        }),
    };
    let loaded = input::load_protein_database(&a.database, &opts)?;
    let mut db = loaded.db;
    if !a.annotations.is_empty() {
        input::load_annotations(&a.annotations, &mut db, opts.knockoffs)?;
    }
    let mut idx = SearchIndex::from_target_db(&db, &a.decoy_prefix);
    idx.decoy_suffix = a.decoy_suffix.clone();

    let mut targets = read_all(&a.psms)?;
    let decoys = read_all(&a.decoy_psms)?;
    let proteoform_fdr = a.proteoform_fdr.unwrap_or(a.fdr.max(0.05));
    if idx.has_alternative_entries {
        output::share_isobaric_references(&mut targets, &idx, a.fdr.max(proteoform_fdr));
    }
    output::apply_class_fdr(&mut targets, &decoys, &idx);
    let knockoffs =
        output::apply_knockoff_fdr(&mut targets, &idx, proteoform_fdr, a.knockoffs > 0)
            .knockoff_rows;
    let cutoffs = output::Cutoffs::new(&targets, &idx, a.fdr, proteoform_fdr);
    let inf = output::infer_protein_groups(&targets, &decoys, &idx, |id| {
        targets[id].q_value <= cutoffs.of(id) && !knockoffs.contains(id)
    });

    let path = {
        let mut p = a.output.clone().into_os_string();
        p.push(".pg.parquet");
        std::path::PathBuf::from(p)
    };
    if let Some(dir) = path.parent().filter(|d| !d.as_os_str().is_empty()) {
        std::fs::create_dir_all(dir)?;
    }
    output::write_pg_parquet(&path, &inf, &runs, runs.len() == 1)?;
    let at = |q: f64| inf.groups.iter().filter(|g| g.q_value <= q).count();
    eprintln!(
        "Wrote {} ({} target groups over {} run(s): {} at protein q<=0.01, {} at q<=0.05; {} decoy groups)",
        path.display(),
        inf.groups.len(),
        runs.len(),
        at(0.01),
        at(0.05),
        inf.decoy_groups.len()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn run_names_come_from_percolator_file_names() {
        assert_eq!(
            run_name(Path::new("out/a01.percolator.target.psms.txt")),
            "a01"
        );
        assert_eq!(run_name(Path::new("b.tsv")), "b");
    }
}
