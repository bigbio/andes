//! Regression guard for the per-spectrum activation-routing fix
//! (merge commit `bc8cff6` on `rust-implement`).
//!
//! Asserts that scoring scan=28787 of PXD001819's `UPS1_5000amol_R1.mzML`
//! with the **auto-detected** `CID_HighRes_Tryp.param` produces a
//! `RawScore` within tolerance of the Java baseline. Without the fix —
//! when the search defaulted to `HCD_HighRes_Tryp_TMT.param` for every
//! spectrum — the same PSM scored ≈ 108 instead of ≈ 235, dropping a
//! valid identification and inflating the Rust↔Java SpecEValue gap.
//!
//! The two load-bearing assertions are:
//!   1. The mzML parser sets `spec.activation_method == ActivationMethod::CID`
//!      from the `<activation>` cvParam `MS:1000133`. This is what triggers
//!      auto-routing in `bin/msgf-rust` — losing the cvParam in extraction
//!      or in the parser breaks the fix silently.
//!   2. The resulting score is within ±15 of Java's 225 baseline.
//!
//! **Scope**: only scan=28787 is locked in here. Sister scans (28825, 33606,
//! 32395) referenced in the original fix plan need fresh Java baselines —
//! their published numbers were captured under the wrong-param config —
//! so they're deferred until those baselines are re-verified.

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use input::MzMLReader;
use model::activation::ActivationMethod;
use model::amino_acid::AminoAcid;
use model::peptide::Peptide;
use scoring::scoring::score_psm;
use scoring::{Param, RankScorer, ScoredSpectrum};

/// Java MS-GF+ reports RawScore=225 for this PSM when run with the
/// canonical CID/HighRes/Tryp parameters on the unmodified mzML.
const EXPECTED_JAVA_BASELINE_RAWSCORE: i32 = 225;

/// Tolerance covers float-precision and prefix-mass rounding drift
/// between Java and Rust. The verified Rust value at the fix's merge
/// (`bc8cff6`) was 235. Do **not** widen this to make a regressed test
/// pass — investigate the divergence first.
const TOLERANCE: i32 = 15;

/// Fragment tolerance Da used by the production CID search path (see
/// `bin/msgf-rust.rs` and `match_engine.rs` — both use 0.5 Da for CID).
const FRAGMENT_TOLERANCE_DA: f64 = 0.5;

/// Repo-relative path: `astral-speed/rust/crates/scoring` → workspace root.
fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..") // crates/scoring → crates → rust → astral-speed
        .canonicalize()
        .expect("canonicalize workspace root")
}

fn fixture_path() -> PathBuf {
    workspace_root().join("src/test/resources/benchmark/PXD001819/scan_28787.mzML")
}

fn param_path() -> PathBuf {
    workspace_root().join("src/main/resources/ionstat/CID_HighRes_Tryp.param")
}

fn build_peptide_ivneefdqleedtpvyk() -> Peptide {
    // K.IVNEEFDQLEEDTPVYK.L
    // pre='K' (preceding residue in the source protein), post='L'.
    let residues: Vec<AminoAcid> = b"IVNEEFDQLEEDTPVYK"
        .iter()
        .map(|&r| {
            AminoAcid::standard(r).unwrap_or_else(|| panic!("standard AA lookup failed for {:?}", r as char))
        })
        .collect();
    Peptide::new(residues, b'K', b'L')
}

#[test]
fn score_psm_scan_28787_ivneefdqleedtpvyk_matches_java_baseline() {
    // ── 1. Load fixture ────────────────────────────────────────────────
    let fixture = fixture_path();
    assert!(
        fixture.exists(),
        "missing fixture: {fixture:?} — extract scan=28787 from PXD001819 \
         UPS1_5000amol_R1.mzML and place it at this path"
    );
    let file = File::open(&fixture).expect("open fixture mzML");
    let reader = MzMLReader::new(BufReader::new(file));

    let spec = reader
        .filter_map(|r| r.ok())
        .find(|s| s.scan == Some(28787))
        .expect("scan=28787 not found in fixture");

    // ── 2. Activation routing — the load-bearing path ──────────────────
    // Without this cvParam (MS:1000133), the binary would default to HCD
    // and load the wrong `.param` file, regressing the fix silently.
    assert_eq!(
        spec.activation_method,
        Some(ActivationMethod::CID),
        "fixture spectrum lost its <activation> cvParam — auto-routing \
         would fall back to HCD and the score would regress"
    );

    // ── 3. Build scorer with the param Java would pick ─────────────────
    let param_path = param_path();
    let param = Param::load_from_file(&param_path)
        .unwrap_or_else(|e| panic!("load {param_path:?}: {e}"));
    let scorer = RankScorer::new(&param);

    // ── 4. Build the peptide and ScoredSpectrum ────────────────────────
    let peptide = build_peptide_ivneefdqleedtpvyk();
    // Charge 2+ matches the PSM's reported charge in Java's output and
    // the `<cvParam … MS:1000041 … 2>` in the fixture's selectedIon.
    let charge: u8 = 2;
    let scored_spec = ScoredSpectrum::new(&spec, &scorer, charge);

    // ── 5. Score and assert ────────────────────────────────────────────
    let raw_score = score_psm(&scored_spec, &peptide, &scorer, charge, FRAGMENT_TOLERANCE_DA);
    let raw_score_i32 = raw_score as i32;

    let lo = EXPECTED_JAVA_BASELINE_RAWSCORE - TOLERANCE;
    let hi = EXPECTED_JAVA_BASELINE_RAWSCORE + TOLERANCE;
    assert!(
        (lo..=hi).contains(&raw_score_i32),
        "RawScore={raw_score_i32} outside Java baseline window {lo}..={hi}. \
         Verified Rust value at fix merge bc8cff6 was 235. If this assertion \
         fires, investigate the score divergence — DO NOT widen TOLERANCE \
         without root-causing the change."
    );
}
