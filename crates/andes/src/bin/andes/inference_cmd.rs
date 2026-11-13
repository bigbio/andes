//! `andes inference`: protein groups over the rescored PSMs of one or more
//! runs, written as the QPX protein-group view.
//!
//! The database is loaded with the search's options, so protein families,
//! knockoff entries and decoy pairs are the ones the PSMs were searched
//! against. The class and knockoff q-values are applied as in a search, then
//! [`output::infer_protein_groups`] runs once over every run's PSMs: one
//! parsimony and one picked competition, so the q-values are dataset-level.

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use search::SearchIndex;

use crate::cli::{InferenceArgs, ProteinArgs};
use crate::model_select::parse_enzymes;
use crate::search::load_claim_database;

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

/// Every file's rows, each under a key of its own: the row's number among all
/// rows read, in decimal.
///
/// Nothing downstream reads a key's text: the FDR stages and the inference
/// use keys only to find and to tell rows apart, and none of their results
/// depends on the keys (each sorts by rank and breaks ties as the map
/// iterates, which is random either way). A dataset-level PSM file has tens
/// of millions of rows, so keying by the `PSMId` itself (and keeping it again
/// in `psm_id`, and the whole file as one string while parsing) would cost
/// tens of gigabytes. Files are read one line at a time; `psm_id` is left
/// empty. A `PSMId` repeated within a file is still an error, as in
/// [`output::parse_psm_results`], found by a 128-bit hash of each id.
fn read_all(paths: &[PathBuf]) -> Result<HashMap<String, output::PercolatorPsm>, String> {
    let mut all = HashMap::new();
    for p in paths {
        let err = |e: &dyn std::fmt::Display| format!("{}: {e}", p.display());
        let file = std::fs::File::open(p).map_err(|e| err(&e))?;
        let mut reader =
            output::PsmResultReader::new(std::io::BufReader::with_capacity(1 << 20, file))
                .map_err(|e| err(&e))?;
        let mut seen: HashSet<u128> = HashSet::new();
        while let Some(row) = reader.next_row().map_err(|e| err(&e))? {
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
    }
    all.shrink_to_fit();
    Ok(all)
}

/// A 128-bit hash of a PSM id: two SipHash runs with different salts, so two
/// distinct ids of one file collide with odds near 2^-128 per pair.
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

/// Masses of the `sites=annotated` modifications in the search's `--mods`
/// file; none without a file, as the built-in modifications have none.
/// The `sites=annotated` masses and the fixed terminal mods of a mods file.
type SiteMods = (Vec<f64>, Vec<(model::modification::ModLocation, f64)>);

fn annotated_site_masses(mods: Option<&Path>) -> Result<SiteMods, Box<dyn std::error::Error>> {
    let Some(path) = mods else {
        return Ok((Vec::new(), Vec::new()));
    };
    let aa = model::AminoAcidSetBuilder::new_standard()
        .add_mods_from_file(path)
        .map_err(|e| format!("loading mods from {}: {e}", path.display()))?
        .build()
        .map_err(|e| format!("building amino-acid set from {}: {e}", path.display()))?;
    Ok((
        aa.distinct_mods()
            .iter()
            .filter_map(|m| m.annotated_site_mass)
            .collect(),
        aa.fixed_terminal_deltas(),
    ))
}

/// Protein groups over the rescored PSMs in `psms` and `decoy_psms` (as
/// `andes inference`), from `database` read with the search's options.
pub(crate) fn infer_from_files(
    database: &Path,
    a: &ProteinArgs,
    psms: &[PathBuf],
    decoy_psms: &[PathBuf],
) -> Result<output::Inference, Box<dyn std::error::Error>> {
    if psms.len() != decoy_psms.len() {
        return Err(format!(
            "--psms has {} files and --decoy-psms {}: give the decoy results of each target file",
            psms.len(),
            decoy_psms.len()
        )
        .into());
    }
    let (enzyme, _) = parse_enzymes(&a.enzyme)?;
    let (restricted, fixed_terminal) = annotated_site_masses(a.mods.as_deref())?;
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
    let db = load_claim_database(database, &a.annotations, &opts)?
        .loaded
        .db;
    let mut idx = SearchIndex::from_target_db(&db, &a.decoy_prefix);
    idx.decoy_suffix = a.decoy_suffix.clone();

    let mut targets = read_all(psms)?;
    let decoys = read_all(decoy_psms)?;
    let proteoform_fdr = a.proteoform_fdr();
    if idx.has_alternative_entries {
        output::share_isobaric_references(&mut targets, &idx, a.psm_fdr.max(proteoform_fdr));
    }
    output::apply_class_fdr(&mut targets, &decoys, &idx);
    let mut knockoffs =
        output::apply_knockoff_fdr(&mut targets, &idx, proteoform_fdr, a.knockoffs > 0)
            .knockoff_rows;
    // Placements of a `sites=annotated` modification on a knockoff site exist
    // only for the site q-values. Only those modifications make such rows.
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
        &decoys,
        &idx,
        |id| targets[id].q_value <= cutoffs.of(id) && !knockoffs.contains(id),
    ))
}

pub(crate) fn run_inference(a: InferenceArgs) -> Result<(), Box<dyn std::error::Error>> {
    let runs: Vec<String> = if a.runs.is_empty() {
        a.psms.iter().map(|p| run_name(p)).collect()
    } else {
        a.runs.clone()
    };
    let inf = infer_from_files(&a.database, &a.proteins, &a.psms, &a.decoy_psms)?;

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
    )?;
    let at = |q: f64| inf.groups.iter().filter(|g| g.q_value <= q).count();
    eprintln!(
        "Wrote {} ({} of {} target groups at protein q<={} over {} run(s): {} at q<=0.01, {} at q<=0.05)",
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
        let all = read_all(&[a, b]).unwrap();
        assert_eq!(all.len(), 3);
        let first = all.values().find(|p| p.peptide == "R.PEPTIDEK.A").unwrap();
        assert_eq!(first.proteins, "sp|P1;sp|P2");
        assert_eq!(first.score, Some(2.5));
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

    /// Peptides and PSMs counted for the group of a protein with one
    /// annotated phosphosite (S15) and three knockoff sites (S5, S25, S35),
    /// from three PSMs: phospho on the real site, phospho on a knockoff site,
    /// and an unmodified peptide. `mods` is the content of the `--mods` file.
    fn counted_with(mods: Option<&str>) -> (usize, usize) {
        use clap::Parser;
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
        let mut argv = vec!["andes".to_string()];
        if let Some(content) = mods {
            let path = dir.path().join("mods.txt");
            std::fs::write(&path, content).unwrap();
            argv.extend(["--mods".to_string(), path.display().to_string()]);
        }
        let a = ProteinWrap::try_parse_from(argv).unwrap().proteins;
        let inf = infer_from_files(&db, &a, &[targets], &[decoys]).unwrap();
        let group = inf
            .groups
            .iter()
            .find(|g| g.accessions.iter().any(|acc| acc == "sp:P1"))
            .expect("a group for sp:P1");
        (group.peptides, group.psms)
    }

    #[test]
    fn annotated_site_mods_on_knockoff_sites_count_for_no_group() {
        assert_eq!(
            counted_with(Some("79.966331,S,opt,any,Phospho,sites=annotated\n")),
            (2, 2)
        );
    }

    #[test]
    fn knockoff_site_placements_count_without_an_annotated_site_mod() {
        // The same modification allowed on any S: a placement on a knockoff
        // site is an ordinary identification.
        assert_eq!(counted_with(Some("79.966331,S,opt,any,Phospho\n")), (3, 3));
        assert_eq!(counted_with(None), (3, 3));
    }
}
