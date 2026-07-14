//! Per-annotated-ion ASSIGNMENT-AWARE features for the decoy-aware rich-ion-LLR.
//! Single source of truth for train + serve (mirrors frag_features.rs). The
//! chemistry block is reused verbatim from frag_features; the appended block
//! encodes complement coherence + this ion's own match quality.
use crate::frag_features::{extract_frag_features, N_FRAG_FEATURES};
use crate::scoring::fragment_ions::{predict_by_ions, IonKind};
use crate::scoring::ScoredSpectrum;
use model::peptide::Peptide;

pub const FEAT_ION_RANK: usize = N_FRAG_FEATURES;
pub const FEAT_MASS_ERR_PPM: usize = N_FRAG_FEATURES + 1;
pub const FEAT_COMPLEMENT_PRESENT: usize = N_FRAG_FEATURES + 2;
pub const FEAT_COMPLEMENT_PPM_DELTA: usize = N_FRAG_FEATURES + 3;
pub const N_ION_FEATURES: usize = N_FRAG_FEATURES + 4;
pub const COMPLEMENT_PPM_SENTINEL: f32 = 1000.0;

/// Signed ppm of an observed peak m/z `obs` relative to a theoretical m/z `theo`.
#[inline]
fn ppm(obs: f64, theo: f64) -> f64 {
    ((obs - theo) / theo) * 1e6
}

/// Convert a feature tolerance (ppm or Da) into an absolute Da window at `mz`.
#[inline]
fn tol_da(mz: f64, feature_tol: f64, feature_tol_is_ppm: bool) -> f64 {
    if feature_tol_is_ppm {
        mz * feature_tol / 1e6
    } else {
        feature_tol
    }
}

/// The complementary ion of `(kind, position)` in a peptide of length `n`:
/// (B, i) ↔ (Y, n−i). Returns the complement's `(kind, position)`, or `None`
/// when the complementary position is degenerate/out of range
/// (`n − position` must satisfy `1 ≤ n − position < n`).
#[inline]
fn complement(kind: IonKind, position: u32, n: usize) -> Option<(IonKind, u32)> {
    let n = n as u32;
    if position >= n {
        return None;
    }
    let comp_pos = n - position;
    if comp_pos < 1 || comp_pos >= n {
        return None;
    }
    let comp_kind = match kind {
        IonKind::B => IonKind::Y,
        IonKind::Y => IonKind::B,
        // ETD complements: c ↔ z• (N/C-terminal partners across the same bond).
        IonKind::C => IonKind::Z,
        IonKind::Z => IonKind::C,
    };
    Some((comp_kind, comp_pos))
}

/// Per-annotated-ion feature vector: the `frag_features` chemistry block
/// (indices `0..N_FRAG_FEATURES`) followed by four assignment-aware features
/// describing THIS ion's observed match and its complement's coherence.
///
/// `position` is 1-based; `precursor_charge` is the PSM precursor charge and
/// `frag_charge` the fragment charge of this ion. `feature_tol` /
/// `feature_tol_is_ppm` define the peak-match window (ppm → `mz·tol/1e6` Da;
/// else `tol` Da).
///
/// The appended features:
/// - [`FEAT_ION_RANK`]: 1-based intensity rank of the matched peak (0.0 if the
///   ion does not match any peak within tolerance).
/// - [`FEAT_MASS_ERR_PPM`]: signed ppm `((obs − theo)/theo)·1e6` of the match
///   (0.0 when unmatched).
/// - [`FEAT_COMPLEMENT_PRESENT`]: 1.0 when the complementary ion (same frag
///   charge) also matches a peak, else 0.0.
/// - [`FEAT_COMPLEMENT_PPM_DELTA`]: `|ppm_self − ppm_complement|` when both
///   match, else [`COMPLEMENT_PPM_SENTINEL`].
#[allow(clippy::too_many_arguments)]
pub fn extract_ion_features(
    peptide: &Peptide,
    scored_spec: &ScoredSpectrum<'_>,
    kind: IonKind,
    position: u32,
    precursor_charge: u8,
    frag_charge: u8,
    feature_tol: f64,
    feature_tol_is_ppm: bool,
) -> [f32; N_ION_FEATURES] {
    let mut f = [0.0f32; N_ION_FEATURES];

    // Chemistry block (verbatim from frag_features; nce = 0.0).
    let chem = extract_frag_features(peptide, kind, position, precursor_charge, frag_charge, 0.0);
    f[..N_FRAG_FEATURES].copy_from_slice(&chem);

    let n = peptide.residues.len();

    // Enumerate predicted 1+/2+ b/y ions ONCE and resolve self + complement m/z.
    let predicted = predict_by_ions(peptide, 1..=2);
    let comp = complement(kind, position, n);

    let mut self_mz: Option<f64> = None;
    let mut comp_mz: Option<f64> = None;
    for ion in &predicted {
        if ion.charge != frag_charge {
            continue;
        }
        if ion.kind == kind && ion.position == position {
            self_mz = Some(ion.mz);
        }
        if let Some((ck, cp)) = comp {
            if ion.kind == ck && ion.position == cp {
                comp_mz = Some(ion.mz);
            }
        }
    }

    // This ion's own match: rank + signed ppm.
    let self_ppm: Option<f64> = if let Some(mz) = self_mz {
        let td = tol_da(mz, feature_tol, feature_tol_is_ppm);
        // nearest_peak_full returns (rank, intensity, peak_mz).
        if let Some((rank, _intensity, obs_mz)) = scored_spec.nearest_peak_full(mz, td) {
            let p = ppm(obs_mz, mz);
            f[FEAT_ION_RANK] = rank as f32;
            f[FEAT_MASS_ERR_PPM] = p as f32;
            Some(p)
        } else {
            None
        }
    } else {
        None
    };

    // Complement coherence: present + ppm-delta.
    f[FEAT_COMPLEMENT_PPM_DELTA] = COMPLEMENT_PPM_SENTINEL;
    if let (Some(sp), Some(cmz)) = (self_ppm, comp_mz) {
        let td = tol_da(cmz, feature_tol, feature_tol_is_ppm);
        if let Some((_rank, _intensity, obs_mz)) = scored_spec.nearest_peak_full(cmz, td) {
            let comp_ppm = ppm(obs_mz, cmz);
            f[FEAT_COMPLEMENT_PRESENT] = 1.0;
            f[FEAT_COMPLEMENT_PPM_DELTA] = (sp - comp_ppm).abs() as f32;
        }
    }

    f
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::fragment_ions::predict_by_ions;
    use crate::scoring::ScoredSpectrum;
    use crate::testutil::tiny_param_with_ions;
    use crate::RankScorer;
    use model::amino_acid::AminoAcid;
    use model::peptide::Peptide;
    use model::spectrum::Spectrum;

    fn pep(seq: &str) -> Peptide {
        Peptide::new(
            seq.bytes().map(|b| AminoAcid::standard(b).unwrap()).collect(),
            b'K',
            b'R',
        )
    }

    fn spectrum(peaks: Vec<(f64, f32)>, mz: f64, z: i32) -> Spectrum {
        let mut peaks = peaks;
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        Spectrum {
            title: "test".to_string(),
            precursor_mz: mz,
            precursor_intensity: None,
            precursor_charge: Some(z),
            rt_seconds: None,
            scan: None,
            peaks,
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
        }
    }

    /// Find a predicted ion's m/z by (kind, position, charge).
    fn ion_mz(p: &Peptide, kind: IonKind, position: u32, charge: u8) -> f64 {
        predict_by_ions(p, 1..=2)
            .into_iter()
            .find(|i| i.kind == kind && i.position == position && i.charge == charge)
            .unwrap_or_else(|| panic!("ion {kind:?} {position} z{charge} not predicted"))
            .mz
    }

    #[test]
    fn complement_present_set_when_both_halves_match() {
        // PEPTIDEK, n=8. b2's complement is y6 (n - 2 = 6), same frag charge 1.
        let p = pep("PEPTIDEK");
        let b2 = ion_mz(&p, IonKind::B, 2, 1);
        let y6 = ion_mz(&p, IonKind::Y, 6, 1);
        // Place exact-m/z peaks at b2 and its complement y6 (+ a junk peak).
        let spec = spectrum(vec![(b2, 1000.0), (y6, 800.0), (50.0, 5.0)], 500.0, 2);
        let scorer = RankScorer::new(&tiny_param_with_ions());
        let ss = ScoredSpectrum::new(&spec, &scorer, 2);

        let f = extract_ion_features(&p, &ss, IonKind::B, 2, 2, 1, 20.0, true);

        assert_eq!(f.len(), N_ION_FEATURES);
        assert_eq!(f[FEAT_COMPLEMENT_PRESENT], 1.0);
        assert!(
            f[FEAT_COMPLEMENT_PPM_DELTA] < 50.0,
            "complement ppm delta {} should be small (both exact)",
            f[FEAT_COMPLEMENT_PPM_DELTA]
        );
        // self matched → rank set, well below the sentinel.
        assert!(f[FEAT_ION_RANK] >= 1.0, "self rank should be set");
    }

    #[test]
    fn complement_absent_uses_sentinel() {
        // Only b2 present; y6 (its complement) absent.
        let p = pep("PEPTIDEK");
        let b2 = ion_mz(&p, IonKind::B, 2, 1);
        let spec = spectrum(vec![(b2, 1000.0), (50.0, 5.0)], 500.0, 2);
        let scorer = RankScorer::new(&tiny_param_with_ions());
        let ss = ScoredSpectrum::new(&spec, &scorer, 2);

        let f = extract_ion_features(&p, &ss, IonKind::B, 2, 2, 1, 20.0, true);

        assert_eq!(f.len(), N_ION_FEATURES);
        assert_eq!(f[FEAT_COMPLEMENT_PRESENT], 0.0);
        assert_eq!(f[FEAT_COMPLEMENT_PPM_DELTA], COMPLEMENT_PPM_SENTINEL);
    }
}
