//! Site-determining-ion / mod-localization features for the PTM-refinement
//! cascade. ORTHOGONAL to RankScore: where RankScore measures overall b/y-rank
//! evidence, these features measure whether the spectrum LOCALIZES a variable
//! modification — i.e. whether the b/y fragments that carry the mod's mass
//! shift (and bracket the modified residue) are actually observed.
//!
//! A correctly-localized true modified peptide has its site-determining ions
//! matched; a false or mislocalized candidate of the same mass does not.
//!
//! All features are 0.0 for an unmodified peptide (no variable-mod site), so
//! emitting them is purely additive — the default (non-modified) path is
//! unchanged.
//!
//! ## Span / shift conventions (1-based positions)
//!
//! For a peptide of length `n` with a variable mod at 1-based residue position
//! `p`:
//! - `b_i` spans residues `1..=i` ⇒ it carries the mod (is "shifted") iff
//!   `p <= i`.
//! - `y_j` spans residues `(n-j+1)..=n` ⇒ it is shifted iff `p > n - j`.
//!
//! ## Bracketing (localization)
//!
//! Position `p` is localized when a matched ion ends just *before* `p` and a
//! matched ion just *includes* `p`, so the mass jump pins the mod to `p`:
//! - N-side: `b_{p-1}` (unshifted) + `b_p` (shifted) both matched.
//! - C-side: `y_{n-p}` (unshifted) + `y_{n-p+1}` (shifted) both matched.

use crate::scoring::fragment_ions::{predict_by_ions, IonKind};
use crate::scoring::ScoredSpectrum;
use model::peptide::Peptide;

/// Index of [`mod_site_features`] output: count of matched SHIFTED ions
/// (fragments that carry the mod's mass shift). Confirms the mod is present.
pub const FEAT_SITE_SHIFTED_MATCHED: usize = 0;
/// Matched shifted ions / total theoretical shifted ions.
pub const FEAT_SITE_SHIFTED_FRAC: usize = 1;
/// Summed intensity of matched shifted ions / summed intensity of ALL matched
/// ions (0.0 if nothing matched).
pub const FEAT_SITE_INTENS_FRAC: usize = 2;
/// 1.0 if ANY mod site has a matched bracketing (site-determining) ion pair
/// that localizes the mod, else 0.0.
pub const FEAT_SITE_LOCALIZED: usize = 3;
/// Number of bracketing (site-determining) matched ions summed over mod sites.
pub const FEAT_SITE_DET_COUNT: usize = 4;
/// Number of features returned by [`mod_site_features`].
pub const N_MOD_SITE_FEATURES: usize = 5;

/// Convert a feature tolerance (ppm or Da) into an absolute Da window at `mz`.
/// Mirrors `ion_features::tol_da` exactly for train/serve parity.
#[inline]
fn tol_da(mz: f64, feature_tol: f64, feature_tol_is_ppm: bool) -> f64 {
    if feature_tol_is_ppm {
        mz * feature_tol / 1e6
    } else {
        feature_tol
    }
}

/// Is the charge-1 b/y ion `(kind, position)` shifted (mod-bearing) for a
/// variable mod at 1-based residue position `mod_pos` in a peptide of length
/// `n`? `b_i` carries the mod iff `mod_pos <= i`; `y_j` iff `mod_pos > n - j`.
#[inline]
fn ion_is_shifted(kind: IonKind, position: u32, mod_pos: u32, n: u32) -> bool {
    match kind {
        // c-ions carry the mod like b (prefix), z•-ions like y (suffix).
        IonKind::B | IonKind::C => mod_pos <= position,
        IonKind::Y | IonKind::Z => mod_pos > n - position,
    }
}

/// Site-determining-ion / mod-localization features for a (possibly modified)
/// peptide vs a spectrum. All 0.0 for an unmodified peptide. Orthogonal to
/// RankScore — measures whether the spectrum LOCALIZES the variable mod(s).
///
/// Uses `feature_tol` / `feature_tol_is_ppm` exactly as `extract_ion_features`
/// (ppm → `mz·tol/1e6` Da; else `tol` Da) so train and serve agree.
///
/// Variable-mod sites are residues carrying a non-`fixed` `Modification`
/// (`aa.mod_` Some with `mod_.fixed == false`). Fixed mods (e.g.
/// Carbamidomethyl on every C) are not localization-relevant and are ignored.
pub fn mod_site_features(
    peptide: &Peptide,
    scored_spec: &ScoredSpectrum<'_>,
    feature_tol: f64,
    feature_tol_is_ppm: bool,
) -> [f32; N_MOD_SITE_FEATURES] {
    let mut f = [0.0f32; N_MOD_SITE_FEATURES];

    let n = peptide.residues.len();
    if n < 2 {
        return f;
    }

    // 1-based positions of VARIABLE-modified residues. A residue is a mod site
    // iff it carries a Modification that is NOT fixed.
    let mod_sites: Vec<u32> = peptide
        .residues
        .iter()
        .enumerate()
        .filter_map(|(i, aa)| {
            aa.mod_
                .as_ref()
                .filter(|m| !m.fixed)
                .map(|_| (i + 1) as u32)
        })
        .collect();

    if mod_sites.is_empty() {
        return f;
    }

    let n_u = n as u32;

    // Enumerate charge-1 b/y ions once; resolve each to a peak match.
    // `matched[(kind, position)] -> intensity` for matched ions only.
    let predicted = predict_by_ions(peptide, 1..=1);

    let mut shifted_matched: u32 = 0;
    let mut shifted_total: u32 = 0;
    let mut shifted_intens_sum: f64 = 0.0;
    let mut all_matched_intens_sum: f64 = 0.0;

    // Track which (kind, position) ions matched, for the bracketing test.
    // Positions are 1..=n-1 for each series; store as two boolean rows.
    let mut b_matched = vec![false; n]; // index by position (0 unused)
    let mut y_matched = vec![false; n];

    for ion in &predicted {
        // mz for charge-1 ions.
        let td = tol_da(ion.mz, feature_tol, feature_tol_is_ppm);
        let matched = scored_spec.nearest_peak_full(ion.mz, td);

        // An ion is "shifted" if its span includes ANY mod site.
        let is_shifted = mod_sites
            .iter()
            .any(|&p| ion_is_shifted(ion.kind, ion.position, p, n_u));
        if is_shifted {
            shifted_total += 1;
        }

        if let Some((_rank, intensity, _peak_mz)) = matched {
            all_matched_intens_sum += intensity as f64;
            if is_shifted {
                shifted_matched += 1;
                shifted_intens_sum += intensity as f64;
            }
            let pos = ion.position as usize;
            match ion.kind {
                IonKind::B | IonKind::C => {
                    if pos < b_matched.len() {
                        b_matched[pos] = true;
                    }
                }
                IonKind::Y | IonKind::Z => {
                    if pos < y_matched.len() {
                        y_matched[pos] = true;
                    }
                }
            }
        }
    }

    f[FEAT_SITE_SHIFTED_MATCHED] = shifted_matched as f32;
    f[FEAT_SITE_SHIFTED_FRAC] = if shifted_total > 0 {
        shifted_matched as f32 / shifted_total as f32
    } else {
        0.0
    };
    f[FEAT_SITE_INTENS_FRAC] = if all_matched_intens_sum > 0.0 {
        (shifted_intens_sum / all_matched_intens_sum) as f32
    } else {
        0.0
    };

    // Bracketing / localization. For each mod site p:
    //   N-side pair: b_{p-1} (unshifted) + b_p (shifted) both matched.
    //   C-side pair: y_{n-p} (unshifted) + y_{n-p+1} (shifted) both matched.
    // FEAT_SITE_LOCALIZED = 1 if ANY site has a matched bracketing pair.
    // FEAT_SITE_DET_COUNT = number of matched site-determining ions summed
    //   over mod sites (each member of a confirmed pair counts).
    let mut any_localized = false;
    let mut det_count: u32 = 0;
    for &p in &mod_sites {
        let pu = p as usize;

        // N-side: positions p-1 (before) and p (includes). p-1 == 0 means the
        // mod is on residue 1, so there is no b_{p-1} prefix to anchor against.
        if pu >= 2 {
            let before = b_matched[pu - 1]; // b_{p-1}
            let includes = pu < n && b_matched[pu]; // b_p (b_n doesn't exist)
            if before && includes {
                any_localized = true;
                det_count += 2;
            }
        }

        // C-side: y_{n-p} (unshifted, before from C-side) and y_{n-p+1}
        // (shifted, includes p). n-p == 0 means the mod is on the last
        // residue, so there is no y_{n-p} suffix to anchor against.
        let np = n_u - p; // n - p
        if np >= 1 {
            let before = y_matched[np as usize]; // y_{n-p}
            let inc_pos = (np + 1) as usize; // y_{n-p+1}
            let includes = inc_pos < n && y_matched[inc_pos];
            if before && includes {
                any_localized = true;
                det_count += 2;
            }
        }
    }

    f[FEAT_SITE_LOCALIZED] = if any_localized { 1.0 } else { 0.0 };
    f[FEAT_SITE_DET_COUNT] = det_count as f32;

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
    use model::modification::{ModLocation, Modification, ResidueSpec};
    use model::peptide::Peptide;
    use model::spectrum::Spectrum;
    use std::sync::Arc;

    /// Build an unmodified peptide from a residue string.
    fn pep(seq: &str) -> Peptide {
        Peptide::new(
            seq.bytes()
                .map(|b| AminoAcid::standard(b).unwrap())
                .collect(),
            b'K',
            b'R',
        )
    }

    /// Build a peptide from `seq`, attaching a VARIABLE mod of `delta` Da to the
    /// residue at 0-based index `mod_idx`. Mirrors how the model tests construct
    /// a modified residue (`AminoAcid::with_mod` + a `Modification { fixed:false }`).
    fn pep_with_var_mod(seq: &str, mod_idx: usize, delta: f64) -> Peptide {
        let bytes: Vec<u8> = seq.bytes().collect();
        let m = Modification {
            name: "Oxidation".to_string(),
            mass_delta: delta,
            residue: ResidueSpec::Specific(bytes[mod_idx]),
            location: ModLocation::Anywhere,
            fixed: false, // VARIABLE — this marks a localization-relevant site.
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
        };
        let arc = Arc::new(m);
        let residues: Vec<AminoAcid> = bytes
            .iter()
            .enumerate()
            .map(|(i, &b)| {
                let aa = AminoAcid::standard(b).unwrap();
                if i == mod_idx {
                    aa.with_mod(arc.clone())
                } else {
                    aa
                }
            })
            .collect();
        Peptide::new(residues, b'K', b'R')
    }

    /// Build a peptide from `seq` with a FIXED mod on the residue at `mod_idx`.
    fn pep_with_fixed_mod(seq: &str, mod_idx: usize, delta: f64) -> Peptide {
        let bytes: Vec<u8> = seq.bytes().collect();
        let m = Modification {
            name: "Carbamidomethyl".to_string(),
            mass_delta: delta,
            residue: ResidueSpec::Specific(bytes[mod_idx]),
            location: ModLocation::Anywhere,
            fixed: true,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
        };
        let arc = Arc::new(m);
        let residues: Vec<AminoAcid> = bytes
            .iter()
            .enumerate()
            .map(|(i, &b)| {
                let aa = AminoAcid::standard(b).unwrap();
                if i == mod_idx {
                    aa.with_mod(arc.clone())
                } else {
                    aa
                }
            })
            .collect();
        Peptide::new(residues, b'K', b'R')
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

    /// All charge-1 b/y ion m/z values of `p`, as (kind, position, mz).
    fn ions_of(p: &Peptide) -> Vec<(IonKind, u32, f64)> {
        predict_by_ions(p, 1..=1)
            .into_iter()
            .map(|i| (i.kind, i.position, i.mz))
            .collect()
    }

    #[test]
    fn unmodified_peptide_returns_all_zero() {
        let p = pep("PEPTIDEK");
        // Put peaks at every b/y ion so matching is not the issue — the mod-site
        // detector must short-circuit to all-zero because there is no var mod.
        let peaks: Vec<(f64, f32)> = ions_of(&p)
            .into_iter()
            .map(|(_, _, mz)| (mz, 1000.0))
            .collect();
        let spec = spectrum(peaks, 500.0, 2);
        let scorer = RankScorer::new(&tiny_param_with_ions());
        let ss = ScoredSpectrum::new(&spec, &scorer, 2);

        let f = mod_site_features(&p, &ss, 20.0, true);
        assert_eq!(f, [0.0; N_MOD_SITE_FEATURES]);
    }

    #[test]
    fn fixed_mod_is_not_a_localization_site() {
        // A FIXED mod (e.g. Carbamidomethyl) is not localization-relevant: the
        // detector must treat the peptide as having no variable-mod site.
        let p = pep_with_fixed_mod("PEPCIDEK", 3, 57.02146);
        let peaks: Vec<(f64, f32)> = ions_of(&p)
            .into_iter()
            .map(|(_, _, mz)| (mz, 1000.0))
            .collect();
        let spec = spectrum(peaks, 500.0, 2);
        let scorer = RankScorer::new(&tiny_param_with_ions());
        let ss = ScoredSpectrum::new(&spec, &scorer, 2);

        let f = mod_site_features(&p, &ss, 20.0, true);
        assert_eq!(f, [0.0; N_MOD_SITE_FEATURES]);
    }

    #[test]
    fn modified_peptide_with_shifted_peaks_is_localized() {
        // PEPM[+15.99]TIDEK: variable Oxidation on M (0-based index 3, 1-based
        // position 4). Synthesize peaks at ALL b/y ions so the shifted ions and
        // the bracketing pair are present.
        let p = pep_with_var_mod("PEPMTIDEK", 3, 15.99491);
        let peaks: Vec<(f64, f32)> = ions_of(&p)
            .into_iter()
            .map(|(_, _, mz)| (mz, 1000.0))
            .collect();
        let spec = spectrum(peaks, 500.0, 2);
        let scorer = RankScorer::new(&tiny_param_with_ions());
        let ss = ScoredSpectrum::new(&spec, &scorer, 2);

        let f = mod_site_features(&p, &ss, 20.0, true);
        assert!(
            f[FEAT_SITE_SHIFTED_MATCHED] > 0.0,
            "expected matched shifted ions, got {}",
            f[FEAT_SITE_SHIFTED_MATCHED]
        );
        assert_eq!(
            f[FEAT_SITE_LOCALIZED], 1.0,
            "mod should be localized when b_{{p-1}}/b_p (and y pair) all match"
        );
        assert!(
            f[FEAT_SITE_DET_COUNT] >= 2.0,
            "at least one bracketing pair (2 ions) should be counted, got {}",
            f[FEAT_SITE_DET_COUNT]
        );
        // Shifted fraction is in (0, 1]; intensity fraction in (0, 1].
        assert!(f[FEAT_SITE_SHIFTED_FRAC] > 0.0 && f[FEAT_SITE_SHIFTED_FRAC] <= 1.0);
        assert!(f[FEAT_SITE_INTENS_FRAC] > 0.0 && f[FEAT_SITE_INTENS_FRAC] <= 1.0);
    }

    #[test]
    fn modified_peptide_without_shifted_peaks_is_not_localized() {
        // Same modified peptide, but only the UNSHIFTED prefix ions are present
        // (b1..b3 span residues before the M at position 4; they never carry the
        // mod). No shifted ion matches ⇒ shifted_matched == 0 and localized == 0.
        let p = pep_with_var_mod("PEPMTIDEK", 3, 15.99491);
        let n = p.length() as u32;
        let unshifted_b: Vec<(f64, f32)> = predict_by_ions(&p, 1..=1)
            .into_iter()
            .filter(|i| {
                matches!(i.kind, IonKind::B) && !ion_is_shifted(i.kind, i.position, 4, n)
                // mod at pos 4
            })
            .map(|i| (i.mz, 1000.0))
            .collect();
        assert!(
            !unshifted_b.is_empty(),
            "fixture must have some unshifted b ions"
        );
        let spec = spectrum(unshifted_b, 500.0, 2);
        let scorer = RankScorer::new(&tiny_param_with_ions());
        let ss = ScoredSpectrum::new(&spec, &scorer, 2);

        let f = mod_site_features(&p, &ss, 20.0, true);
        assert_eq!(
            f[FEAT_SITE_SHIFTED_MATCHED], 0.0,
            "no shifted ion present ⇒ shifted_matched must be 0"
        );
        assert_eq!(
            f[FEAT_SITE_LOCALIZED], 0.0,
            "no shifted ion ⇒ the mod cannot be localized"
        );
        assert_eq!(f[FEAT_SITE_DET_COUNT], 0.0);
        // Unshifted ions DID match, so intensity fraction (shifted/all) is 0.
        assert_eq!(f[FEAT_SITE_INTENS_FRAC], 0.0);
    }

    #[test]
    fn ion_is_shifted_span_logic() {
        // n = 9, mod at position 4 (the M in PEPMTIDEK).
        // b_i shifted iff 4 <= i.
        assert!(!ion_is_shifted(IonKind::B, 3, 4, 9));
        assert!(ion_is_shifted(IonKind::B, 4, 4, 9));
        assert!(ion_is_shifted(IonKind::B, 8, 4, 9));
        // y_j shifted iff 4 > 9 - j  ⇔  j > 5.
        assert!(!ion_is_shifted(IonKind::Y, 5, 4, 9)); // 4 > 4 is false
        assert!(ion_is_shifted(IonKind::Y, 6, 4, 9)); // 4 > 3 is true
        assert!(ion_is_shifted(IonKind::Y, 8, 4, 9));
    }
}
