//! Peptide-channel backbone chance LLR for glycopeptides (redesign step 6).
//!
//! `ChanceLLR_masked(p) = Σ_matched w_I · (−ln p_chance)` evaluated on the
//! PEPTIDE-CHANNEL spectrum (oxonium and this backbone's Y rungs masked before ranks,
//! base peak and local density are computed — `ScoredSpectrum::new_with_excluded_mz`).
//!
//! Two things the bare b/y scorer cannot do, both from the design doc §3.2:
//!
//! 1. **HexNAc-stub max-over-forms.** Under HCD an ion that spans the glycosite is
//!    observed bare, with one HexNAc, or with two (chitobiose stub) — which one depends
//!    on collision energy, not on the peptide. Decorating every spanning ion with the
//!    FULL glycan replaced the bare masses and regressed 41%; here a spanning ion takes
//!    the MAX over the stub forms so the bare form is never lost and a stub hit is
//!    never invisible. Non-spanning ions are bare only.
//! 2. **Isotope-gated fragment charges to z−1.** Glycopeptides are z3–z6 and the bare
//!    backbone is often 1.5–2.5 kDa, so z2/z3 fragments are real; but a z≥3 fragment
//!    match at 20 ppm is a coin toss in a dense spectrum, so it only counts when the
//!    M+1 isotope peak is also present at `mz + 1.00335/z`.
//!
//! Per-ion chance probability is `ρ(mz)·2·tol`, ρ = local peak density of the masked
//! spectrum, exactly the form `cz_structure_features` uses on the c/z channel — the
//! only c/z variant that ever separated. Both outputs are per-predicted-ion means so
//! long decoys do not win by ion count.
//!
//! Decoy symmetry: the spectrum mask depends on backbone mass only; the ion set is
//! the peptide's own, as for every backbone score. Nothing here reads the composition.

use crate::scoring::fragment_ions::{predict_by_ions, IonKind};
use crate::scoring::scored_spectrum::ScoredSpectrum;
use model::peptide::Peptide;

/// HexNAc residue mass (C8H13NO5), the stub left on a glycosite-spanning HCD fragment.
pub const HEXNAC: f64 = 203.079_373;
/// C13−C12 spacing used for the isotope confirmation of z≥3 fragments.
const ISOTOPE: f64 = 1.003_354_83;
/// Half-width (Da) of the local-density window; same as the c/z channel.
const DENSITY_HW: f64 = 50.0;
/// Fragment charge at or above which a match must be isotope-confirmed.
const ISOTOPE_GATE_MIN_CHARGE: u8 = 3;
/// Hard cap on the fragment charge sweep (doc: `1..min(z−1, 5)`).
const MAX_FRAG_CHARGE: u8 = 5;

/// Neutral-mass shifts a glycosite-spanning ion is scored under: bare, +HexNAc,
/// +2 HexNAc. The bare form is always first so `stub_forms = &DEFAULT_STUB_FORMS[..1]`
/// degrades to a plain (mask-only) chance LLR for A/B.
pub const DEFAULT_STUB_FORMS: [f64; 3] = [0.0, HEXNAC, 2.0 * HEXNAC];

/// Output of [`masked_backbone_llr`]. All zero for `n < 2` or an empty spectrum.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct MaskedBackboneLlr {
    /// Mean over predicted (ion, charge) pairs of `w_I · (−ln p_chance)`, `w_I` the
    /// matched peak's base-peak fraction on the masked spectrum. The PIN column
    /// `ChanceLlrMasked`.
    pub chance_llr: f32,
    /// Fraction of predicted (ion, charge) pairs that matched (any admitted form).
    /// The PIN column `ExplainedMasked`.
    pub explained: f32,
    /// Predicted (ion, charge) pairs after the charge sweep.
    pub n_predicted: u32,
    /// Matched pairs whose best form carried a HexNAc stub (shift > 0). A nonzero
    /// count is the wiring proof that stub forms are consulted at all.
    pub n_stub_matched: u32,
}

/// Fragment charges scored for a precursor of charge `precursor_charge`:
/// `1..=min(z−1, 5)`, never empty.
#[inline]
pub fn fragment_charge_max(precursor_charge: u8) -> u8 {
    precursor_charge.saturating_sub(1).clamp(1, MAX_FRAG_CHARGE)
}

/// True when a b/y ion at `position` spans `glycosite` (0-based residue index) in a
/// peptide of length `n`, i.e. the fragment carries the glycan (or its stub).
#[inline]
pub fn spans_glycosite(kind: IonKind, position: usize, n: usize, glycosite: usize) -> bool {
    match kind {
        // b_k = prefix residues [0, k): spans when glycosite < k.
        IonKind::B | IonKind::C => glycosite < position,
        // y_k = suffix residues [n−k, n): spans when glycosite >= n−k.
        IonKind::Y | IonKind::Z => glycosite >= n.saturating_sub(position),
    }
}

/// Peptide-channel backbone chance LLR of `peptide` on the (already masked)
/// `scored_spec`. `glycosite` is a 0-based residue index; a value `>= n` means no
/// ion spans it, which reduces this to a plain mask-only chance LLR.
pub fn masked_backbone_llr(
    scored_spec: &ScoredSpectrum,
    peptide: &Peptide,
    glycosite: usize,
    precursor_charge: u8,
    tol_ppm: f64,
    stub_forms: &[f64],
) -> MaskedBackboneLlr {
    let n = peptide.length();
    if n < 2 || stub_forms.is_empty() {
        return MaskedBackboneLlr::default();
    }
    let (peaks, _) = scored_spec.active_peaks_and_ranks();
    if peaks.is_empty() {
        return MaskedBackboneLlr::default();
    }
    let base = peaks.iter().map(|&(_, i)| i).fold(0.0f32, f32::max).max(1e-9);

    let zmax = fragment_charge_max(precursor_charge);
    let ions = predict_by_ions(peptide, 1..=zmax);
    let n_predicted = ions.len() as u32;
    if n_predicted == 0 {
        return MaskedBackboneLlr::default();
    }

    let bare_only = &stub_forms[..1];
    let mut sum_llr = 0.0f64;
    let mut n_matched = 0u32;
    let mut n_stub_matched = 0u32;
    for ion in &ions {
        let z = ion.charge as f64;
        let forms = if spans_glycosite(ion.kind, ion.position as usize, n, glycosite) {
            stub_forms
        } else {
            bare_only
        };
        // Best form = the one with the largest intensity-weighted surprise, so a
        // bright stub peak beats a dim bare one and vice versa; a form that is not
        // observed contributes nothing (no penalty: the form choice is instrument-,
        // not peptide-, determined).
        let mut best: Option<(f64, f64)> = None; // (term, shift)
        for &shift in forms {
            let mz = ion.mz + shift / z;
            let tol_da = (mz * tol_ppm * 1e-6).max(0.01);
            let Some((_, intensity, _)) = scored_spec.nearest_peak_full(mz, tol_da) else {
                continue;
            };
            if ion.charge >= ISOTOPE_GATE_MIN_CHARGE
                && scored_spec
                    .nearest_peak_full(mz + ISOTOPE / z, tol_da)
                    .is_none()
            {
                continue;
            }
            let obs_w = ((intensity / base) as f64).clamp(0.0, 1.0);
            let rho = scored_spec.local_peak_density(mz, DENSITY_HW);
            let p_chance = (rho * 2.0 * tol_da).clamp(1e-12, 1.0);
            let term = obs_w * (-p_chance.ln()).max(0.0);
            if best.is_none_or(|(t, _)| term > t) {
                best = Some((term, shift));
            }
        }
        if let Some((term, shift)) = best {
            sum_llr += term;
            n_matched += 1;
            if shift > 0.0 {
                n_stub_matched += 1;
            }
        }
    }
    MaskedBackboneLlr {
        chance_llr: (sum_llr / n_predicted as f64) as f32,
        explained: n_matched as f32 / n_predicted as f32,
        n_predicted,
        n_stub_matched,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::rank_scorer::RankScorer;
    use crate::testutil::tiny_param;
    use model::amino_acid::AminoAcid;
    use model::spectrum::Spectrum;

    fn pep(seq: &[u8]) -> Peptide {
        let residues: Vec<AminoAcid> =
            seq.iter().map(|&r| AminoAcid::standard(r).unwrap()).collect();
        Peptide::new(residues, b'_', b'-')
    }

    fn spectrum(peaks: Vec<(f64, f32)>, z: u8) -> Spectrum {
        Spectrum {
            title: "t".into(),
            precursor_mz: 900.0,
            precursor_intensity: None,
            precursor_charge: Some(z as i32),
            rt_seconds: None,
            scan: None,
            peaks,
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
        }
    }

    fn ion_mz(p: &Peptide, kind: IonKind, position: u32, charge: u8) -> f64 {
        predict_by_ions(p, 1..=charge)
            .into_iter()
            .find(|i| i.kind == kind && i.position == position && i.charge == charge)
            .map(|i| i.mz)
            .expect("ion exists")
    }

    /// P E P T N Y S K: sequon N at index 4. y3 (Y,S,K) is bare; b6 (P,E,P,T,N,Y)
    /// spans the site.
    const SEQ: &[u8] = b"PEPTNYSK";
    const SITE: usize = 4;

    #[test]
    fn empty_spectrum_and_short_peptide_are_zero() {
        let scorer = RankScorer::new(&tiny_param());
        let s = spectrum(vec![], 2);
        let ss = ScoredSpectrum::new(&s, &scorer, 2);
        assert_eq!(
            masked_backbone_llr(&ss, &pep(SEQ), SITE, 2, 20.0, &DEFAULT_STUB_FORMS),
            MaskedBackboneLlr::default()
        );
        let s = spectrum(vec![(300.0, 10.0)], 2);
        let ss = ScoredSpectrum::new(&s, &scorer, 2);
        assert_eq!(
            masked_backbone_llr(&ss, &pep(b"K"), 0, 2, 20.0, &DEFAULT_STUB_FORMS),
            MaskedBackboneLlr::default()
        );
    }

    #[test]
    fn bare_non_spanning_ion_scores() {
        let p = pep(SEQ);
        let y3 = ion_mz(&p, IonKind::Y, 3, 1);
        let scorer = RankScorer::new(&tiny_param());
        let s = spectrum(vec![(y3, 100.0), (1500.0, 5.0)], 2);
        let ss = ScoredSpectrum::new(&s, &scorer, 2);
        let r = masked_backbone_llr(&ss, &p, SITE, 2, 20.0, &DEFAULT_STUB_FORMS);
        assert_eq!(r.n_predicted, 14, "2·(8−1) ions at charge 1");
        assert!(r.chance_llr > 0.0, "{r:?}");
        assert!((r.explained - 1.0 / 14.0).abs() < 1e-6, "{r:?}");
        assert_eq!(r.n_stub_matched, 0);
    }

    #[test]
    fn spanning_ion_matches_on_hexnac_stub_form_only_when_admitted() {
        let p = pep(SEQ);
        let b6_stub = ion_mz(&p, IonKind::B, 6, 1) + HEXNAC;
        let scorer = RankScorer::new(&tiny_param());
        let s = spectrum(vec![(b6_stub, 100.0), (1500.0, 5.0)], 2);
        let ss = ScoredSpectrum::new(&s, &scorer, 2);

        let with = masked_backbone_llr(&ss, &p, SITE, 2, 20.0, &DEFAULT_STUB_FORMS);
        assert_eq!(with.n_stub_matched, 1, "{with:?}");
        assert!(with.chance_llr > 0.0);

        let bare = masked_backbone_llr(&ss, &p, SITE, 2, 20.0, &DEFAULT_STUB_FORMS[..1]);
        assert_eq!(bare, MaskedBackboneLlr { n_predicted: 14, ..Default::default() });

        // A non-spanning ion never takes a stub form: y2 (S,K) + HexNAc is not an ion.
        let y2_stub = ion_mz(&p, IonKind::Y, 2, 1) + HEXNAC;
        let s = spectrum(vec![(y2_stub, 100.0), (1500.0, 5.0)], 2);
        let ss = ScoredSpectrum::new(&s, &scorer, 2);
        let r = masked_backbone_llr(&ss, &p, SITE, 2, 20.0, &DEFAULT_STUB_FORMS);
        assert_eq!(r.n_stub_matched, 0, "{r:?}");
        assert_eq!(r.explained, 0.0);
    }

    #[test]
    fn bare_form_never_lost_when_stub_forms_admitted() {
        // Spanning b6 observed BARE: max-over-forms must still credit it (the
        // regression the full-glycan decoration caused).
        let p = pep(SEQ);
        let b6 = ion_mz(&p, IonKind::B, 6, 1);
        let scorer = RankScorer::new(&tiny_param());
        let s = spectrum(vec![(b6, 100.0), (1500.0, 5.0)], 2);
        let ss = ScoredSpectrum::new(&s, &scorer, 2);
        let r = masked_backbone_llr(&ss, &p, SITE, 2, 20.0, &DEFAULT_STUB_FORMS);
        assert!(r.chance_llr > 0.0, "{r:?}");
        assert_eq!(r.n_stub_matched, 0);
    }

    #[test]
    fn glycosite_past_end_means_no_spanning_ions() {
        let p = pep(SEQ);
        let b6_stub = ion_mz(&p, IonKind::B, 6, 1) + HEXNAC;
        let scorer = RankScorer::new(&tiny_param());
        let s = spectrum(vec![(b6_stub, 100.0), (1500.0, 5.0)], 2);
        let ss = ScoredSpectrum::new(&s, &scorer, 2);
        let r = masked_backbone_llr(&ss, &p, usize::MAX, 2, 20.0, &DEFAULT_STUB_FORMS);
        assert_eq!(r.explained, 0.0, "{r:?}");
    }

    #[test]
    fn fragment_charge_sweep_is_z_minus_one_capped_at_five() {
        assert_eq!(fragment_charge_max(1), 1);
        assert_eq!(fragment_charge_max(2), 1);
        assert_eq!(fragment_charge_max(4), 3);
        assert_eq!(fragment_charge_max(9), 5);
        let p = pep(SEQ);
        let scorer = RankScorer::new(&tiny_param());
        let s = spectrum(vec![(1500.0, 5.0)], 4);
        let ss = ScoredSpectrum::new(&s, &scorer, 4);
        let r = masked_backbone_llr(&ss, &p, SITE, 4, 20.0, &DEFAULT_STUB_FORMS);
        assert_eq!(r.n_predicted, 14 * 3);
    }

    #[test]
    fn z3_fragment_needs_isotope_confirmation() {
        let p = pep(SEQ);
        let y5_z3 = ion_mz(&p, IonKind::Y, 5, 3);
        let scorer = RankScorer::new(&tiny_param());
        // Without the M+1 peak: not counted.
        let s = spectrum(vec![(y5_z3, 100.0), (1500.0, 5.0)], 4);
        let ss = ScoredSpectrum::new(&s, &scorer, 4);
        let r = masked_backbone_llr(&ss, &p, SITE, 4, 20.0, &DEFAULT_STUB_FORMS);
        assert_eq!(r.explained, 0.0, "{r:?}");
        // With it: counted.
        let s = spectrum(
            vec![(y5_z3, 100.0), (y5_z3 + ISOTOPE / 3.0, 40.0), (1500.0, 5.0)],
            4,
        );
        let ss = ScoredSpectrum::new(&s, &scorer, 4);
        let r = masked_backbone_llr(&ss, &p, SITE, 4, 20.0, &DEFAULT_STUB_FORMS);
        assert!(r.explained > 0.0 && r.chance_llr > 0.0, "{r:?}");
        // A z1 ion needs no confirmation.
        let y3 = ion_mz(&p, IonKind::Y, 3, 1);
        let s = spectrum(vec![(y3, 100.0), (1500.0, 5.0)], 4);
        let ss = ScoredSpectrum::new(&s, &scorer, 4);
        let r = masked_backbone_llr(&ss, &p, SITE, 4, 20.0, &DEFAULT_STUB_FORMS);
        assert!(r.explained > 0.0, "{r:?}");
    }

    #[test]
    fn intensity_weights_the_surprise_and_output_is_finite() {
        let p = pep(SEQ);
        let y3 = ion_mz(&p, IonKind::Y, 3, 1);
        let scorer = RankScorer::new(&tiny_param());
        let bright = spectrum(vec![(y3, 100.0), (1500.0, 100.0)], 2);
        let dim = spectrum(vec![(y3, 10.0), (1500.0, 100.0)], 2);
        let a = masked_backbone_llr(
            &ScoredSpectrum::new(&bright, &scorer, 2), &p, SITE, 2, 20.0, &DEFAULT_STUB_FORMS,
        );
        let b = masked_backbone_llr(
            &ScoredSpectrum::new(&dim, &scorer, 2), &p, SITE, 2, 20.0, &DEFAULT_STUB_FORMS,
        );
        assert!(a.chance_llr > b.chance_llr, "{a:?} vs {b:?}");
        assert!(a.chance_llr.is_finite() && b.chance_llr.is_finite());
    }
}
