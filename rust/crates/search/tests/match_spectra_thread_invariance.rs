//! Thread-count invariance: match_spectra must produce bit-identical output
//! regardless of the Rayon thread count, because each spectrum's full pipeline
//! (scoring + GF + spec_e_value assignment) runs entirely on one Rayon worker
//! — there is no FP-accumulation non-determinism across thread counts, only
//! wall time changes.

mod common;
use common::*;

use std::fs::File;
use std::io::BufReader;

use input::{FastaReader, MgfReader};
use model::{Enzyme, Tolerance};
use model::tolerance::PrecursorTolerance;
use search::{match_spectra, SearchIndex, SearchParams, TopNQueue};

fn run_search(thread_count: usize) -> (Vec<TopNQueue>, Vec<search::candidate_gen::Candidate>) {
    // Use a scoped pool via `install` (NOT `build_global`) so the test does
    // not conflict with any global pool initialization done elsewhere.
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build()
        .expect("build pool");

    let target = FastaReader::load_all(BufReader::new(
        File::open(fixture("rust/test-fixtures/BSA.fasta")).unwrap(),
    ))
    .unwrap();
    let idx = SearchIndex::from_target_db(&target, "XXX_");
    let aa = aa_set();
    let scorer = rank_scorer();

    let mut params = SearchParams::default_tryptic(aa.clone());
    params.enzyme = Enzyme::Trypsin;
    params.precursor_tolerance = PrecursorTolerance::symmetric(Tolerance::Ppm(20.0));
    params.charge_range = 2..=3;
    params.isotope_error_range = -1..=2;

    let mgf_file = File::open(fixture("rust/test-fixtures/test.mgf")).unwrap();
    let spectra: Vec<_> = MgfReader::new(BufReader::new(mgf_file))
        .filter_map(|r| r.ok())
        .collect();

    pool.install(|| match_spectra(&spectra, &idx, &params, &scorer, 0.5, "XXX_"))
}

#[test]
fn match_spectra_output_invariant_across_thread_counts() {
    let (q1, cands_a) = run_search(1);
    let (q4, cands_b) = run_search(4);

    assert_eq!(q1.len(), q4.len(), "queue count differs");

    let mut spectra_with_psms = 0;
    for (i, (qa, qb)) in q1.iter().zip(q4.iter()).enumerate() {
        let psms_a = qa.clone().into_sorted_vec();
        let psms_b = qb.clone().into_sorted_vec();
        assert_eq!(
            psms_a.len(),
            psms_b.len(),
            "spectrum {}: PSM count differs ({} vs {})",
            i,
            psms_a.len(),
            psms_b.len()
        );
        if !psms_a.is_empty() {
            spectra_with_psms += 1;
            for (j, (a, b)) in psms_a.iter().zip(psms_b.iter()).enumerate() {
                let pep_a: String = cands_a[a.primary_candidate_idx() as usize]
                    .peptide
                    .residues
                    .iter()
                    .map(|aa| aa.residue as char)
                    .collect();
                let pep_b: String = cands_b[b.primary_candidate_idx() as usize]
                    .peptide
                    .residues
                    .iter()
                    .map(|aa| aa.residue as char)
                    .collect();
                assert_eq!(
                    pep_a, pep_b,
                    "spectrum {} PSM rank {}: peptide differs ({} vs {})",
                    i, j, pep_a, pep_b
                );
                assert_eq!(
                    a.charge_used, b.charge_used,
                    "spectrum {} PSM rank {}: charge differs",
                    i, j
                );
                assert_eq!(
                    a.score.to_bits(),
                    b.score.to_bits(),
                    "spectrum {} PSM rank {}: score differs ({} vs {})",
                    i, j, a.score, b.score
                );
                assert_eq!(
                    a.spec_e_value.to_bits(),
                    b.spec_e_value.to_bits(),
                    "spectrum {} PSM rank {}: spec_e_value differs ({} vs {})",
                    i, j, a.spec_e_value, b.spec_e_value
                );
            }
        }
    }
    assert!(
        spectra_with_psms > 0,
        "no spectra produced PSMs — fixture problem"
    );
    eprintln!(
        "Verified bit-identical output across thread counts on {} spectra with PSMs",
        spectra_with_psms
    );
}
