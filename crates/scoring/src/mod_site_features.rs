use crate::scoring::fragment_ions::{predict_by_ions, IonKind};
use crate::scoring::ScoredSpectrum;
use model::peptide::Peptide;

pub const FEAT_SITE_SHIFTED_MATCHED: usize = 0;
pub const FEAT_SITE_SHIFTED_FRAC: usize = 1;
pub const FEAT_SITE_INTENS_FRAC: usize = 2;
pub const FEAT_SITE_LOCALIZED: usize = 3;
pub const FEAT_SITE_DET_COUNT: usize = 4;
pub const N_MOD_SITE_FEATURES: usize = 5;

#[inline]
fn tol_da(mz: f64, feature_tol: f64, feature_tol_is_ppm: bool) -> f64 {
    if feature_tol_is_ppm {
        mz * feature_tol / 1e6
    } else {
        feature_tol
    }
}

#[inline]
fn ion_is_shifted(kind: IonKind, position: u32, mod_pos: u32, n: u32) -> bool {
    match kind {
        IonKind::B | IonKind::C => mod_pos <= position,
        IonKind::Y | IonKind::Z => mod_pos > n - position,
    }
}

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

    let predicted = predict_by_ions(peptide, 1..=1);

    let mut shifted_matched: u32 = 0;
    let mut shifted_total: u32 = 0;
    let mut shifted_intens_sum: f64 = 0.0;
    let mut all_matched_intens_sum: f64 = 0.0;

    let mut b_matched = vec![false; n];
    let mut y_matched = vec![false; n];

    for ion in &predicted {
        let td = tol_da(ion.mz, feature_tol, feature_tol_is_ppm);
        let matched = scored_spec.nearest_peak_full(ion.mz, td);

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

    let mut any_localized = false;
    let mut det_count: u32 = 0;
    for &p in &mod_sites {
        let pu = p as usize;

        if pu >= 2 {
            let before = b_matched[pu - 1];
            let includes = pu < n && b_matched[pu];
            if before && includes {
                any_localized = true;
                det_count += 2;
            }
        }

        let np = n_u - p;
        if np >= 1 {
            let before = y_matched[np as usize];
            let inc_pos = (np + 1) as usize;
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

    fn pep(seq: &str) -> Peptide {
        Peptide::new(
            seq.bytes()
                .map(|b| AminoAcid::standard(b).unwrap())
                .collect(),
            b'K',
            b'R',
        )
    }

    fn pep_with_var_mod(seq: &str, mod_idx: usize, delta: f64) -> Peptide {
        let bytes: Vec<u8> = seq.bytes().collect();
        let m = Modification {
            name: "Oxidation".to_string(),
            mass_delta: delta,
            residue: ResidueSpec::Specific(bytes[mod_idx]),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        let arc = model::modification::leak_mod(m);
        let residues: Vec<AminoAcid> = bytes
            .iter()
            .enumerate()
            .map(|(i, &b)| {
                let aa = AminoAcid::standard(b).unwrap();
                if i == mod_idx {
                    aa.with_mod(arc)
                } else {
                    aa
                }
            })
            .collect();
        Peptide::new(residues, b'K', b'R')
    }

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
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        let arc = model::modification::leak_mod(m);
        let residues: Vec<AminoAcid> = bytes
            .iter()
            .enumerate()
            .map(|(i, &b)| {
                let aa = AminoAcid::standard(b).unwrap();
                if i == mod_idx {
                    aa.with_mod(arc)
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
            isolation_target_mz: None,
            mobility_window: None,
        }
    }

    fn ions_of(p: &Peptide) -> Vec<(IonKind, u32, f64)> {
        predict_by_ions(p, 1..=1)
            .into_iter()
            .map(|i| (i.kind, i.position, i.mz))
            .collect()
    }

    #[test]
    fn unmodified_peptide_returns_all_zero() {
        let p = pep("PEPTIDEK");
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
        assert!(f[FEAT_SITE_SHIFTED_FRAC] > 0.0 && f[FEAT_SITE_SHIFTED_FRAC] <= 1.0);
        assert!(f[FEAT_SITE_INTENS_FRAC] > 0.0 && f[FEAT_SITE_INTENS_FRAC] <= 1.0);
    }

    #[test]
    fn modified_peptide_without_shifted_peaks_is_not_localized() {
        let p = pep_with_var_mod("PEPMTIDEK", 3, 15.99491);
        let n = p.length() as u32;
        let unshifted_b: Vec<(f64, f32)> = predict_by_ions(&p, 1..=1)
            .into_iter()
            .filter(|i| {
                matches!(i.kind, IonKind::B) && !ion_is_shifted(i.kind, i.position, 4, n)
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
        assert_eq!(f[FEAT_SITE_INTENS_FRAC], 0.0);
    }

    #[test]
    fn ion_is_shifted_span_logic() {
        assert!(!ion_is_shifted(IonKind::B, 3, 4, 9));
        assert!(ion_is_shifted(IonKind::B, 4, 4, 9));
        assert!(ion_is_shifted(IonKind::B, 8, 4, 9));
        assert!(!ion_is_shifted(IonKind::Y, 5, 4, 9));
        assert!(ion_is_shifted(IonKind::Y, 6, 4, 9));
        assert!(ion_is_shifted(IonKind::Y, 8, 4, 9));
    }
}
