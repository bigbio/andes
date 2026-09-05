//! End-to-end integration test for own-derived geometry.
//!
//! Proves the full path: `derive_geometry` (corpus → partition/segment geometry)
//! → `RankScorer` → `StatsAccumulator` → `Estimator` yields a scorable model
//! whose geometry comes from the corpus, not from the seed template.

use std::path::Path;

use model::amino_acid::AminoAcid;
use model::peptide::Peptide;
use model::spectrum::Spectrum;
use scoring_crate::param_model::Param;
use scoring_crate::scoring::rank_scorer::RankScorer;

use model_train::accumulate::StatsAccumulator;
use model_train::counts::CountStats;
use model_train::estimate::{Estimator, EstimatorConfig};
use model_train::geometry::{derive_geometry, GeometryConfig};
use model_train::store::ModelStore;

const PROTON: f64 = 1.007_276_49;
const H2O: f64 = 18.010_564_68;

fn make_peptide(seq: &[u8]) -> Peptide {
    let residues = seq
        .iter()
        .map(|&r| AminoAcid::standard(r).unwrap())
        .collect();
    Peptide::new(residues, b'_', b'-')
}

/// Bundled hcd_qexactive_tryp — used only for its NON-geometry metadata
/// (tolerance / activation / deconvolution); its geometry is what we replace.
fn load_base() -> Param {
    let bundled = Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../resources/models"
    ));
    ModelStore::open(bundled)
        .expect("open bundled model store")
        .load_param("hcd_qexactive_tryp")
        .expect("load hcd_qexactive_tryp")
}

#[test]
fn derived_geometry_trains_a_scorable_model_end_to_end() {
    let base = load_base();

    // Own corpus of (charge, parent_mass): the partition tiers are derived from
    // THIS, not inherited from the seed.
    let corpus: Vec<(i32, f32)> = vec![
        (2, 700.0),
        (2, 800.0),
        (2, 900.0),
        (2, 1000.0),
        (3, 1200.0),
        (3, 1500.0),
    ];
    // occupancy=1, max=2 → each charge gets min(n_masses, 2) = 2 tiers.
    let cfg = GeometryConfig {
        num_segments: 2,
        max_rank: 150,
        mass_tier_occupancy: 1,
        max_mass_tiers: 2,
        max_fragment_charge: 2,
    };
    let derived = derive_geometry(&corpus, &base, &cfg);

    // Geometry is the derived one, not the seed's.
    assert_eq!(derived.num_segments, 2);
    assert_eq!(derived.max_rank, 150);
    assert_eq!(derived.min_charge, 2);
    assert_eq!(derived.max_charge, 3);
    // 2 charges × 2 tiers × 2 segments = 8 partitions, each with a Noise entry.
    assert_eq!(derived.partitions.len(), 8);
    for part in &derived.partitions {
        let frags = derived
            .frag_off_table
            .get(part)
            .expect("frag entry per partition");
        assert!(
            frags.iter().any(|f| f.ion_type.is_noise()),
            "Noise per partition"
        );
    }

    // (1) RankScorer accepts the geometry-only template (empty tables → no panic).
    let scorer = RankScorer::new(&derived);

    // (2) Accumulate a spectrum for a charge-2 peptide into the DERIVED partitions.
    let seq = b"PEPTIDE";
    let peptide = make_peptide(seq);
    let charge = 2u8;
    let precursor_mz = (peptide.mass() + charge as f64 * PROTON) / charge as f64;

    // Place peaks AT the peptide's own b/y ion m/z (charge 1) so the model's
    // enumerated fragments actually match and record real ranks.
    let res: Vec<f64> = seq
        .iter()
        .map(|&r| AminoAcid::standard(r).unwrap().mass)
        .collect();
    let mut peaks: Vec<(f64, f32)> = Vec::new();
    let mut cum = 0.0;
    for (i, &m) in res.iter().take(res.len() - 1).enumerate() {
        cum += m;
        peaks.push((cum + PROTON, (2000 - i as i32 * 120).max(100) as f32)); // b(i+1)+
    }
    let mut cum_y = 0.0;
    for (i, &m) in res.iter().rev().take(res.len() - 1).enumerate() {
        cum_y += m;
        peaks.push((
            cum_y + H2O + PROTON,
            (1900 - i as i32 * 120).max(100) as f32,
        )); // y(i+1)+
    }
    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let spectrum = Spectrum {
        title: "synthetic".into(),
        precursor_mz,
        precursor_intensity: None,
        precursor_charge: Some(charge as i32),
        rt_seconds: None,
        scan: None,
        peaks,
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
    };
    let mut stats = CountStats::default();
    StatsAccumulator::new(&scorer).accumulate(&mut stats, &spectrum, &peptide, charge);

    // ★ THE BUG + FIX: a GEOMETRY-ONLY template (empty rank_dist_table) must still
    // enumerate ion MEMBERSHIP from frag_off_table and record matched-fragment
    // ranks. Pre-fix, the accumulator enumerated ions from the (empty) trained
    // LLR cache → ZERO facts → uniform/flat tables → degenerate, slow scoring.
    // Post-fix it enumerates membership → real ranks land in counts.rank.
    assert!(
        !stats.rank.is_empty(),
        "geometry-only accumulate recorded NO rank facts — trained tables would be flat (the bug)"
    );

    // (3) Estimate fills the learned tables from the accumulated counts.
    let estimated = Estimator::new(EstimatorConfig::default()).estimate(&stats, &derived);

    // (4) The trained own-geometry model is scorable: non-empty rank tables, every
    //     slot strictly positive (Laplace floor), and RankScorer accepts it.
    assert!(
        !estimated.rank_dist_table.is_empty(),
        "trained model has rank tables"
    );
    for ion_tab in estimated.rank_dist_table.values() {
        for v in ion_tab.values() {
            assert!(
                v.iter().all(|&x| x > 0.0),
                "every rank-dist slot strictly positive"
            );
        }
    }
    let _trained_scorer = RankScorer::new(&estimated);
    // The trained model kept the derived geometry.
    assert_eq!(estimated.num_segments, 2);
    assert_eq!(estimated.max_rank, 150);
}
