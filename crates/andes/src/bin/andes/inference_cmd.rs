use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use search::SearchIndex;

use crate::cli::{DatabaseArgs, InferenceArgs, ProteinArgs};
use crate::model_select::parse_enzymes;
use crate::search::load_claim_database;

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

/// The PSMs of every file, and the runs they come from: the distinct values of a file's run
/// column in first-seen order (a pooled Percolator result names each PSM's spectrum file), or
/// the file's own name when it has no such column.
type PsmFiles = (HashMap<String, output::PercolatorPsm>, Vec<String>);

fn read_all(paths: &[PathBuf]) -> Result<PsmFiles, String> {
    let mut all = HashMap::new();
    let mut runs: Vec<String> = Vec::new();
    let mut known_runs: HashSet<String> = HashSet::new();
    for p in paths {
        let err = |e: &dyn std::fmt::Display| format!("{}: {e}", p.display());
        let file = std::fs::File::open(p).map_err(|e| err(&e))?;
        let mut reader =
            output::PsmResultReader::new(std::io::BufReader::with_capacity(1 << 20, file))
                .map_err(|e| err(&e))?;
        let mut seen: HashSet<u128> = HashSet::new();
        let mut named_runs = false;
        while let Some(row) = reader.next_row().map_err(|e| err(&e))? {
            if let Some(run) = row.run.filter(|r| !r.is_empty()) {
                named_runs = true;
                if !known_runs.contains(run) {
                    known_runs.insert(run.to_string());
                    runs.push(run.to_string());
                }
            }
            if !seen.insert(id_hash(row.psm_id)) {
                return Err(err(&format!(
                    "Percolator result has duplicate PSMId `{}` (join key must be unique)",
                    row.psm_id
                )));
            }
            let psm = output::PercolatorPsm {
                pooled_q_value: None,
                psm_id: String::new(),
                q_value: row.q_value,
                pep: row.pep,
                score: row.score,
                peptide: row.peptide.to_string(),
                proteins: row.proteins(),
            };
            all.insert(all.len().to_string(), psm);
        }
        if !named_runs {
            let run = run_name(p);
            if known_runs.insert(run.clone()) {
                runs.push(run);
            }
        }
    }
    all.shrink_to_fit();
    Ok((all, runs))
}

fn id_hash(id: &str) -> u128 {
    use std::hash::{BuildHasher, BuildHasherDefault, DefaultHasher, Hash, Hasher};
    let half = |salt: u8| {
        let mut h = BuildHasherDefault::<DefaultHasher>::default().build_hasher();
        salt.hash(&mut h);
        id.hash(&mut h);
        h.finish()
    };
    (u128::from(half(0)) << 64) | u128::from(half(1))
}

type SiteMods = (Vec<f64>, Vec<(model::modification::ModLocation, f64)>);

fn annotated_site_masses(db: &DatabaseArgs) -> Result<SiteMods, Box<dyn std::error::Error>> {
    let aa = crate::model_select::search_aa_set(db, None)
        .map_err(|e| format!("building the amino-acid set: {e}"))?;
    Ok((
        aa.distinct_mods()
            .iter()
            .filter_map(|m| m.annotated_site_mass)
            .collect(),
        aa.fixed_terminal_deltas(),
    ))
}

/// Protein inference over Percolator result files, with the runs the target PSMs come from.
pub(crate) fn infer_from_files(
    database: &Path,
    a: &ProteinArgs,
    psms: &[PathBuf],
    decoy_psms: &[PathBuf],
) -> Result<(output::Inference, Vec<String>), Box<dyn std::error::Error>> {
    if psms.len() != decoy_psms.len() {
        return Err(format!(
            "--psms has {} files and --decoy-psms {}: give the decoy results of each target file",
            psms.len(),
            decoy_psms.len()
        )
        .into());
    }
    let (targets, runs) = read_all(psms)?;
    let (decoys, _) = read_all(decoy_psms)?;
    Ok((infer_from_maps(database, a, targets, &decoys)?, runs))
}

/// Protein inference over Percolator results already in memory (target and decoy PSMs).
pub(crate) fn infer_from_maps(
    database: &Path,
    a: &ProteinArgs,
    mut targets: HashMap<String, output::PercolatorPsm>,
    decoys: &HashMap<String, output::PercolatorPsm>,
) -> Result<output::Inference, Box<dyn std::error::Error>> {
    let (enzyme, _) = parse_enzymes(&a.db.enzyme)?;
    let (restricted, fixed_terminal) = annotated_site_masses(&a.db)?;
    let opts = input::UniProtOptions {
        variants: a.db.uniprot_variants,
        isoforms: a.db.uniprot_isoforms,
        chains: a.db.uniprot_chains,
        evidence: a.db.uniprot_evidence,
        knockoffs: Some(input::KnockoffOptions {
            per_variant: a.db.knockoffs,
            enzyme,
        }),
    };
    let db = load_claim_database(database, &a.db.annotations, &opts)?
        .loaded
        .db;
    let mut idx = SearchIndex::from_target_db(&db, &a.db.decoy_prefix);
    idx.decoy_suffix = a.db.decoy_suffix.clone();

    let proteoform_fdr = a.proteoform_fdr();
    if idx.has_alternative_entries {
        output::share_isobaric_references(&mut targets, &idx, a.psm_fdr.max(proteoform_fdr));
    }
    output::apply_class_fdr(&mut targets, decoys, &idx);
    let mut knockoffs =
        output::apply_knockoff_fdr(&mut targets, &idx, proteoform_fdr, a.db.knockoffs > 0)
            .knockoff_rows;
    if !restricted.is_empty() {
        knockoffs.extend(
            output::site_claims(
                &mut targets,
                &idx,
                a.psm_fdr,
                proteoform_fdr,
                &restricted,
                None,
                &fixed_terminal,
            )
                .knockoff_rows,
        );
    }
    let cutoffs = output::Cutoffs::new(&targets, &idx, a.psm_fdr, proteoform_fdr);
    Ok(output::infer_protein_groups(
        &targets,
        decoys,
        &idx,
        |id| targets[id].q_value <= cutoffs.of(id) && !knockoffs.contains(id),
    ))
}

pub(crate) fn run_inference(a: InferenceArgs) -> Result<(), Box<dyn std::error::Error>> {
    let (inf, file_runs) = infer_from_files(&a.database, &a.proteins, &a.psms, &a.decoy_psms)?;
    let runs: Vec<String> = if a.runs.is_empty() {
        file_runs
    } else {
        a.runs.clone()
    };

    let path = {
        let mut p = a.output.clone().into_os_string();
        p.push(".pg.parquet");
        std::path::PathBuf::from(p)
    };
    if let Some(dir) = path.parent().filter(|d| !d.as_os_str().is_empty()) {
        std::fs::create_dir_all(dir)?;
    }
    output::write_pg_parquet(
        &path,
        &inf,
        &runs,
        runs.len() == 1,
        a.proteins.protein_fdr,
        None,
        None,
    )?;
    let at = |q: f64| inf.groups.iter().filter(|g| g.q_value <= q).count();
    log::info!(
        target: "inference",
        "wrote {} ({} of {} target groups at protein q<={} over {} run(s): {} at q<=0.01, {} at q<=0.05)",
        path.display(),
        at(a.proteins.protein_fdr),
        inf.groups.len(),
        a.proteins.protein_fdr,
        runs.len(),
        at(0.01),
        at(0.05),
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

    const HEADER: &str = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n";

    fn write(dir: &Path, name: &str, rows: &[&str]) -> PathBuf {
        let path = dir.join(name);
        std::fs::write(&path, format!("{HEADER}{}", rows.join("\n"))).unwrap();
        path
    }

    #[test]
    fn read_all_keeps_rows_of_equal_ids_in_different_files() {
        let dir = tempfile::tempdir().unwrap();
        let row = "a_1_1\t2.5\t0.001\t0.0005\tR.PEPTIDEK.A\tsp|P1\tsp|P2\n";
        let a = write(dir.path(), "a.txt", &[row]);
        let b = write(
            dir.path(),
            "b.txt",
            &[row, "b_1_1\t1.0\t0.02\t0.03\tK.OTHERK.R\tsp|P3"],
        );
        let (all, runs) = read_all(&[a, b]).unwrap();
        assert_eq!(all.len(), 3);
        assert_eq!(runs, ["a", "b"]);
        let first = all.values().find(|p| p.peptide == "R.PEPTIDEK.A").unwrap();
        assert_eq!(first.proteins, "sp|P1;sp|P2");
        assert_eq!(first.score, Some(2.5));
    }

    #[test]
    fn a_pooled_result_names_the_runs_of_its_run_column() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("pooled.target.psms.txt");
        std::fs::write(
            &path,
            "PSMId\tfilename\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n\
             a_1_1\tFr02\t2.5\t0.001\t0.0005\tR.PEPTIDEK.A\tsp|P1\n\
             a_2_1\tFr01\t2.0\t0.002\t0.001\tR.OTHERK.A\tsp|P2\n\
             a_3_1\tFr02\t1.0\t0.003\t0.002\tR.THIRDK.A\tsp|P3\n",
        )
        .unwrap();
        let plain = write(
            dir.path(),
            "single.txt",
            &["b_1_1\t1.0\t0.02\t0.03\tK.OTHERK.R\tsp|P3"],
        );
        let (all, runs) = read_all(&[path, plain]).unwrap();
        assert_eq!(all.len(), 4);
        assert_eq!(runs, ["Fr02", "Fr01", "single"]);
    }

    #[test]
    fn read_all_rejects_a_repeated_id_within_a_file() {
        let dir = tempfile::tempdir().unwrap();
        let row = "a_1_1\t2.5\t0.001\t0.0005\tR.PEPTIDEK.A\tsp|P1";
        let a = write(dir.path(), "a.txt", &[row, row]);
        let err = read_all(&[a]).unwrap_err();
        assert!(err.contains("duplicate PSMId `a_1_1`"), "{err}");
    }

    #[derive(clap::Parser)]
    struct ProteinWrap {
        #[command(flatten)]
        proteins: ProteinArgs,
    }

    fn counted_with(args: &[&str], config: Option<&str>) -> (usize, usize) {
        use clap::{CommandFactory, FromArgMatches};
        let dir = tempfile::tempdir().unwrap();
        let db = dir.path().join("db.peff");
        std::fs::write(
            &db,
            "# PEFF 1.0\n# //\n# DbName=test\n# //\n\
             >sp:P1 \\DbUniqueId=P1 \\PName=Site test \\Length=40 \
             \\ModResUnimod=(15|UNIMOD:21|Phospho)\n\
             MGGGSPGGGKAGGGSPGGGKLGGGSPGGGKVGGGSPGGGK\n",
        )
        .unwrap();
        let targets = write(
            dir.path(),
            "t.txt",
            &[
                "s_1_1\t5.0\t0.001\t0.001\tK.AGGGS+79.96633PGGGK.L\tsp:P1",
                "s_2_1\t4.0\t0.001\t0.001\tK.LGGGS+79.96633PGGGK.V\tsp:P1",
                "s_3_1\t3.0\t0.001\t0.001\tK.VGGGSPGGGK.-\tsp:P1",
            ],
        );
        let decoys = write(
            dir.path(),
            "d.txt",
            &["s_4_1\t-1.0\t0.5\t0.9\tK.GGPGSGGAGK.L\tXXX_sp:P1"],
        );
        let mut argv: Vec<String> = ["andes"].iter().chain(args).map(|s| s.to_string()).collect();
        if let Some(content) = config {
            let path = dir.path().join("run.yaml");
            std::fs::write(&path, content).unwrap();
            argv.extend(["--config".to_string(), path.display().to_string()]);
        }
        let m = ProteinWrap::command().try_get_matches_from(argv).unwrap();
        let mut a = ProteinWrap::from_arg_matches(&m).unwrap().proteins;
        crate::config::apply_protein_config(&mut a, &m).unwrap();
        let (inf, _) = infer_from_files(&db, &a, &[targets], &[decoys]).unwrap();
        let group = inf
            .groups
            .iter()
            .find(|g| g.accessions.iter().any(|acc| acc == "sp:P1"))
            .expect("a group for sp:P1");
        (group.peptides, group.psms)
    }

    #[test]
    fn annotated_site_mods_on_knockoff_sites_count_for_no_group() {
        let config = "search:\n  variable_mods:\n    - {name: Phospho, sites: S, restrict_to: annotated}\n";
        assert_eq!(counted_with(&[], Some(config)), (2, 2));
        assert_eq!(
            counted_with(&["--variable-mods", "Phospho (S)"], Some(config)),
            (3, 3)
        );
    }

    #[test]
    fn knockoff_site_placements_count_without_an_annotated_site_mod() {
        assert_eq!(counted_with(&["--variable-mods", "Phospho (S)"], None), (3, 3));
        assert_eq!(counted_with(&[], None), (3, 3));
    }
}
