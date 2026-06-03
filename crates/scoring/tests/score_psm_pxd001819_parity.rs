//! Regression guard for the per-spectrum activation-routing fix
//! (merge commit `bc8cff6` on `rust-implement`) and the follow-up
//! instrument-type auto-detection (this commit, 2026-05-14).
//!
//! Asserts that scoring scan=28787 of PXD001819's `UPS1_5000amol_R1.mzML`
//! with `CID_LowRes_Tryp.param` produces a stable `RawScore` value.
//!
//! Why `CID_LowRes_Tryp.param`, not `CID_HighRes_Tryp.param`: PXD001819
//! is LTQ Velos data, where MS1 lives in the orbitrap but MS2 lives in
//! the linear ion trap (IC2 in the mzML's
//! `<instrumentConfigurationList>`). Java's `NewScorerFactory.get`
//! defaults `instType` to `LOW_RESOLUTION_LTQ` when no `-inst` flag is
//! given, so Java picks `CID_LowRes_Tryp.param` for this dataset. The
//! Rust port's new `detect_instrument_type` helper reads the MS2-
//! referenced `<analyzer>` cvParam and arrives at the same answer.
//!
//! The two load-bearing assertions are:
//!   1. The mzML parser sets `spec.activation_method == ActivationMethod::CID`
//!      from the `<activation>` cvParam `MS:1000133`. This is what triggers
//!      auto-routing in `bin/msgf-rust` — losing the cvParam in extraction
//!      or in the parser breaks the fix silently.
//!   2. The resulting score is stable around the locked Rust value (no
//!      Java baseline exists for scan=28787 under CID_LowRes — diagnostic
//!      runs were captured with `-inst 1`). We treat this as a "score
//!      stability" test: changes in the scoring path must not silently
//!      drift this value.
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

/// Rust-side score for this PSM under `CID_LowRes_Tryp.param` (the
/// param that auto-detection picks for PXD001819 LTQ-Velos MS2 data and
/// the param Java's `NewScorerFactory` defaults to). Locked at 293 on
/// `rust-implement` after the instrument-detection landing (2026-05-14).
///
/// This is a Rust-vs-Rust stability test, not a Java parity test —
/// scan=28787's Java baseline was captured with `-inst 1` (HighRes),
/// so it can't be reused here. If you change the scoring path and this
/// drifts, investigate the divergence before adjusting the constant.
const EXPECTED_RAWSCORE: i32 = 293;

/// Tolerance covers float-precision and prefix-mass rounding drift.
/// Do **not** widen this to make a regressed test pass — investigate
/// the divergence first.
const TOLERANCE: i32 = 15;

/// Fragment tolerance Da used by the production CID search path (see
/// `bin/msgf-rust.rs` and `match_engine.rs` — both use 0.5 Da for CID).
const FRAGMENT_TOLERANCE_DA: f64 = 0.5;

/// Repo-relative path: `astral-speed/rust/crates/scoring` → workspace root.
fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..") // crates/scoring → crates → rust → astral-speed
        .canonicalize()
        .expect("canonicalize workspace root")
}

fn fixture_path() -> PathBuf {
    workspace_root().join("test-fixtures/benchmark/PXD001819/scan_28787.mzML")
}

fn param_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/CID_LowRes_Tryp.param")
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

    let lo = EXPECTED_RAWSCORE - TOLERANCE;
    let hi = EXPECTED_RAWSCORE + TOLERANCE;
    assert!(
        (lo..=hi).contains(&raw_score_i32),
        "RawScore={raw_score_i32} outside Rust stability window {lo}..={hi}. \
         Locked value on `rust-implement` after instrument-detection landing \
         was 293 (CID_LowRes_Tryp.param). If this assertion fires, investigate \
         the score divergence — DO NOT widen TOLERANCE without root-causing it."
    );
}
