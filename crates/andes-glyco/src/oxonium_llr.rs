//! Composition-conditioned oxonium log-likelihood ratio.
//!
//! WHY this exists separately from [`crate::oxonium`]: the gate in that module is
//! SPECTRUM-level — it answers "is this scan glycosylated at all?" and returns the
//! same number for every candidate of a scan, so it has zero per-candidate
//! resolving power. Its `sialic_consistency` companion is per-candidate but only
//! flips the SIGN of an observed intensity, which means a composition that CLAIMS
//! a monosaccharide whose diagnostic ion is ABSENT pays nothing: `-0.0 == 0.0`.
//! Absence of a diagnostic ion is therefore never evidence against a composition.
//!
//! The field does the opposite. pGlyco3 (Zeng et al., Nat Commun 2021) rejects a
//! glycan whose monosaccharide-diagnostic B ions are missing, and PTM-Shepherd's
//! diagnostic-ion module (Polasky et al., Mol Cell Proteomics 2022,
//! doi:10.1016/j.mcpro.2022.100205) scores oxonium per monosaccharide CLASS with
//! intensity-weighted hit ratios above 1 and MISS ratios below 1 — i.e. a signed,
//! two-sided likelihood ratio, not a one-sided reward. Ablating the miss side in
//! that work cost few identifications overall but lost control of sialic-acid and
//! fucose false positives, which is exactly the failure mode seen here.
//!
//! Split of work: [`oxonium_profile`] is computed ONCE per spectrum (the cheap
//! spectrum-level part); [`oxonium_composition_llr`] is the per-candidate part and
//! is the only place a diagnostic ion gets to influence WHICH composition wins.
//!
//! Nothing here is wired into the selector or the PIN feature vector.

use crate::backbone::SpectrumStats;
use crate::glycan_db::GlycanComp;
use crate::oxonium::{NEUAC_OXONIUM_MZ, NEUGC_OXONIUM_MZ};

/// Monosaccharide classes that have a usable diagnostic oxonium signature.
///
/// Grouped by CLASS rather than by individual ion because the individual ions of
/// one monosaccharide are strongly correlated (they are successive water/ring
/// losses of the same residue), so treating them as independent evidence would
/// multiply-count one observation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MonoClass {
    HexNAc,
    HexHexNAc,
    Fuc,
    NeuAc,
    NeuGc,
}

/// Iteration order for every per-class array in this module.
pub const MONO_CLASSES: [MonoClass; 5] = [
    MonoClass::HexNAc,
    MonoClass::HexHexNAc,
    MonoClass::Fuc,
    MonoClass::NeuAc,
    MonoClass::NeuGc,
];

/// HexNAc oxonium series (HexNAc+H and its ring/water losses).
pub const HEXNAC_OXONIUM_MZ: [f64; 4] = [138.05496, 168.06552, 186.07608, 204.08665];
/// HexHexNAc oxonium (366.1395) plus the HexNAc2 ion (407.1664), which is the
/// chitobiose-core B ion reported by PTM-Shepherd for the HexNAc class extension.
pub const HEXHEXNAC_OXONIUM_MZ: [f64; 2] = [366.13947, 407.16637];
/// Fucose-containing B ions. Free Fuc oxonium (147.065) is deliberately EXCLUDED:
/// it is isobaric-adjacent to abundant peptide immonium/side-chain ions and is the
/// classical false-positive source. The retained ions all carry Fuc bound to
/// HexNAc, which is what a core- or antenna-fucosylated N-glycan actually yields:
/// 350.1446 = HexNAc+Fuc, 512.1974 = HexNAc+Hex+Fuc.
/// 803.2925 (NeuAc+Hex+HexNAc+Fuc+H, the sialyl-Lewis B ion; NOT HexNAc2+Hex+Fuc,
/// which is 715.277) is deliberately absent: it needs NeuAc as well as Fuc, so
/// crediting it to the Fuc class alone would reward a Fuc-only composition with an
/// ion it cannot produce. It belongs to a joint NeuAc∧Fuc class if ever wired.
/// (Diagnostic-ion set as tabulated by Polasky et al. 2022, doi:10.1016/j.mcpro.2022.100205.)
pub const FUC_OXONIUM_MZ: [f64; 2] = [350.14456, 512.19738];
/// NeuAc class = the two Neu5Ac oxonium ions already used elsewhere in the crate
/// plus 657.2349 (NeuAc-Hex-HexNAc), the sialylated-antenna B ion. 657 is grouped
/// with NeuAc rather than given its own class because it is only produced when a
/// sialic acid is present, so it is evidence for the same latent variable.
pub const NEUAC_CLASS_OXONIUM_MZ: [f64; 3] = [NEUAC_OXONIUM_MZ[0], NEUAC_OXONIUM_MZ[1], 657.23488];

impl MonoClass {
    /// Diagnostic singly-charged m/z values for this class.
    pub fn diagnostic_mz(self) -> &'static [f64] {
        match self {
            MonoClass::HexNAc => &HEXNAC_OXONIUM_MZ,
            MonoClass::HexHexNAc => &HEXHEXNAC_OXONIUM_MZ,
            MonoClass::Fuc => &FUC_OXONIUM_MZ,
            MonoClass::NeuAc => &NEUAC_CLASS_OXONIUM_MZ,
            MonoClass::NeuGc => &NEUGC_OXONIUM_MZ,
        }
    }

    /// Index into the per-class arrays of [`OxoniumProfile`] / [`OxoniumLlrParams`].
    #[inline]
    pub fn idx(self) -> usize {
        match self {
            MonoClass::HexNAc => 0,
            MonoClass::HexHexNAc => 1,
            MonoClass::Fuc => 2,
            MonoClass::NeuAc => 3,
            MonoClass::NeuGc => 4,
        }
    }

    /// Does `comp` CLAIM this monosaccharide class?
    ///
    /// HexHexNAc is claimed only when both residues are present, because the ion
    /// is a disaccharide B ion, not a sum of two independent residues.
    #[inline]
    pub fn claimed_by(self, comp: &GlycanComp) -> bool {
        match self {
            MonoClass::HexNAc => comp.hexnac > 0,
            MonoClass::HexHexNAc => comp.hexnac > 0 && comp.hex > 0,
            MonoClass::Fuc => comp.fuc > 0,
            MonoClass::NeuAc => comp.neuac > 0,
            MonoClass::NeuGc => comp.neugc > 0,
        }
    }
}

/// Observed base-peak-normalised intensity per monosaccharide class.
///
/// One value per class = the MAX normalised intensity over that class's diagnostic
/// ions (max, not sum, for the correlation reason given on [`MonoClass`]).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct OxoniumProfile {
    /// Per-class normalised intensities, indexed by [`MonoClass::idx`].
    pub obs: [f32; 5],
}

impl OxoniumProfile {
    #[inline]
    pub fn get(&self, class: MonoClass) -> f32 {
        self.obs[class.idx()]
    }
}

/// Spectrum-level pass: one profile per scan, reused by every candidate.
///
/// `stats.base` is the base-peak intensity; normalising by it makes the profile
/// comparable across scans of wildly different absolute intensity, which is a
/// precondition for shared (fitted) likelihood constants.
pub fn oxonium_profile(
    peaks: &[(f64, f32)],
    stats: &SpectrumStats,
    tol_ppm: f64,
) -> OxoniumProfile {
    let base = (stats.base as f32).max(1e-9);
    let mut obs = [0.0f32; 5];
    for class in MONO_CLASSES {
        let mut best_class = 0.0f32;
        for &mz in class.diagnostic_mz() {
            // Absolute floor of 0.01 Th mirrors the tolerance handling of the
            // existing gate so the two agree on which peak is "the" oxonium.
            let tol = (mz * tol_ppm / 1e6).max(0.01);
            let mut best = 0.0f32;
            for &(pmz, pi) in peaks {
                if (pmz - mz).abs() <= tol && pi > best {
                    best = pi;
                }
            }
            if best > best_class {
                best_class = best;
            }
        }
        obs[class.idx()] = best_class / base;
    }
    OxoniumProfile { obs }
}

/// Per-class two-component likelihood constants.
///
/// The likelihood ratio is modelled as a logistic in `ln(o + eps)`:
///
/// ```text
/// P(present | o) = sigmoid( slope * (ln(o + eps) - loc) )
/// ln( P(o|present) / P(o|absent) ) = slope * (ln(o + eps) - loc)
/// ```
///
/// so the log-ratio is the LOGIT itself — no exponentials, no division, finite at
/// `o == 0` (it evaluates to `slope * (ln(eps) - loc)`, a large negative number),
/// which is the numerically awkward case a naive `ln(p/q)` would blow up on.
/// `loc` is the normalised intensity at which presence and absence are equally
/// likely; `slope` is how sharply the class discriminates.
///
/// ⚠ EVERY constant below is a PLACEHOLDER chosen from published diagnostic-ion
/// behaviour, NOT fitted here. They are meant to be replaced by a fit on own
/// entrapment-labelled data (empirical `P(o | class present)` vs
/// `P(o | class absent)` per class). They live in this struct precisely so the fit
/// can be swapped in without touching a single call site.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OxoniumLlrParams {
    /// Equal-likelihood normalised intensity, stored as `ln(o)`. Index = class idx.
    pub loc: [f32; 5],
    /// Logistic slope per class.
    pub slope: [f32; 5],
    /// Weight applied when the composition CLAIMS the class (hit side).
    pub w_present: [f32; 5],
    /// Weight applied when the composition does NOT claim the class (miss side).
    pub w_absent: [f32; 5],
    /// Regulariser inside the log, so `o == 0` stays finite.
    pub eps: f32,
    /// Symmetric clamp on the per-class logit before weighting, so one saturated
    /// or one missing ion cannot dominate the whole score.
    pub clamp: f32,
}

impl Default for OxoniumLlrParams {
    fn default() -> Self {
        // ln of the placeholder half-way normalised intensities.
        const LN_1PCT: f32 = -4.60517; // ln(0.01)
        const LN_2PCT: f32 = -3.912_023; // ln(0.02)
        Self {
            loc: [
                0.0,     // HexNAc: unused, weight 0
                0.0,     // HexHexNAc: unused, weight 0
                LN_1PCT, // Fuc
                LN_2PCT, // NeuAc
                LN_2PCT, // NeuGc
            ],
            slope: [0.0, 0.0, 0.9, 1.1, 1.0],
            // WHY HexNAc and HexHexNAc carry ZERO weight: 204.0867 and 366.1395 are
            // produced by essentially every N-glycan, so their presence is constant
            // across the candidate set of a scan and carries no per-candidate
            // information — exactly the spectrum-level degeneracy this module was
            // written to avoid re-introducing. They stay in the profile because the
            // gate and future fits want them, not because they score.
            w_present: [0.0, 0.0, 1.0, 1.0, 0.8],
            // WHY the Fuc MISS weight is much smaller than the NeuAc miss weight:
            // the evidence is asymmetric. Fuc-containing B ions are produced
            // erratically — fucose is labile, migrates during HCD, and core-fucose
            // often survives on the Y ions instead of yielding 350/512 — so a
            // fucosylated glycan frequently shows NO fucose oxonium. Their PRESENCE
            // is still strong evidence (nothing else makes 350.1446), so the hit
            // side keeps full weight while the miss side is discounted. Sialic acid
            // is the opposite: Neu5Ac is the most reliably produced oxonium in the
            // whole series, so its absence is genuinely informative and the miss
            // side keeps full weight.
            w_absent: [0.0, 0.0, 0.25, 1.0, 0.8],
            eps: 1e-4,
            clamp: 6.0,
        }
    }
}

impl OxoniumLlrParams {
    /// Per-class log-likelihood ratio `ln(P(o|present)/P(o|absent))`, clamped.
    #[inline]
    fn class_logit(&self, class: MonoClass, o: f32) -> f32 {
        let i = class.idx();
        let x = (o.max(0.0) + self.eps).ln();
        let s = self.slope[i] * (x - self.loc[i]);
        s.clamp(-self.clamp, self.clamp)
    }
}

/// Per-candidate composition-conditioned oxonium LLR.
///
/// ```text
/// llr = sum over classes m:
///     if comp claims m:  w_present[m] * ln( P(o_m|present) / P(o_m|absent) )
///     else:              w_absent[m]  * ln( P(o_m|absent)  / P(o_m|present) )
/// ```
///
/// The second branch is the term the existing `sialic_consistency` cannot express:
/// it is POSITIVE when a composition correctly declines a monosaccharide the
/// spectrum does not show, and NEGATIVE when the spectrum shows one the
/// composition ignored. Together with the first branch's negative value at
/// `o ≈ 0`, a composition can now be PENALISED for claiming a residue with no
/// diagnostic support — the pGlyco3 / PTM-Shepherd behaviour.
pub fn oxonium_composition_llr(profile: &OxoniumProfile, comp: &GlycanComp) -> f32 {
    oxonium_composition_llr_with(profile, comp, &OxoniumLlrParams::default())
}

/// [`oxonium_composition_llr`] with explicit (e.g. fitted, or cofrag-adjusted)
/// constants.
pub fn oxonium_composition_llr_with(
    profile: &OxoniumProfile,
    comp: &GlycanComp,
    params: &OxoniumLlrParams,
) -> f32 {
    let mut acc = 0.0f32;
    for class in MONO_CLASSES {
        let logit = params.class_logit(class, profile.get(class));
        let i = class.idx();
        acc += if class.claimed_by(comp) {
            params.w_present[i] * logit
        } else {
            params.w_absent[i] * -logit
        };
    }
    acc
}

/// Run-level prevalence of each diagnostic class: the fraction of ALL spectra in
/// the run in which that class's oxonium is observed. Index = [`MonoClass::idx`].
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ClassPrevalence {
    pub frac: [f32; 5],
}

/// Baseline prevalence above which a class is treated as cofragmentation-corrupted.
pub const COFRAG_PREVALENCE_BASELINE: f32 = 0.5;

/// Scale down the weights of classes whose diagnostic ion fires in an implausibly
/// large fraction of the run's spectra.
///
/// WHY: this is PTM-Shepherd's own named residual failure mode. In brain tissue
/// data the Neu5Ac oxonium appears in >90% of scans while only ~40% of matches are
/// sialylated — co-fragmenting sialylated precursors deposit the ion into the
/// spectra of unsialylated ones. When that happens the ion stops being conditioned
/// on the candidate's composition, so BOTH sides of the ratio are corrupted: the
/// hit side rewards every sialylated candidate indiscriminately and the miss side
/// punishes every unsialylated one. Both weights are therefore scaled by the same
/// factor, driving the class smoothly toward the zero-weight (uninformative)
/// treatment the core HexNAc ions already get.
///
/// Factor is `1` up to [`COFRAG_PREVALENCE_BASELINE`], then falls linearly to `0`
/// at prevalence `1`. The run-level prevalence is NOT computed here — it needs a
/// pass over the whole run and belongs to the caller that owns the spectrum list.
pub fn cofrag_downweight(
    params: &OxoniumLlrParams,
    prevalence: &ClassPrevalence,
) -> OxoniumLlrParams {
    let mut out = *params;
    for class in MONO_CLASSES {
        let i = class.idx();
        let p = prevalence.frac[i].clamp(0.0, 1.0);
        let factor = if p <= COFRAG_PREVALENCE_BASELINE {
            1.0
        } else {
            (1.0 - p) / (1.0 - COFRAG_PREVALENCE_BASELINE)
        };
        out.w_present[i] *= factor;
        out.w_absent[i] *= factor;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn comp(hexnac: u8, hex: u8, fuc: u8, neuac: u8, neugc: u8) -> GlycanComp {
        GlycanComp {
            hexnac,
            hex,
            fuc,
            neuac,
            neugc,
            mass: 0.0,
        }
    }

    /// Core-only N-glycan spectrum: strong HexNAc/HexHexNAc, nothing else.
    fn core_peaks() -> Vec<(f64, f32)> {
        vec![
            (204.08665, 1000.0),
            (186.07608, 400.0),
            (366.13947, 600.0),
            (500.0, 50.0),
        ]
    }

    fn with_ion(mut peaks: Vec<(f64, f32)>, mz: f64, inten: f32) -> Vec<(f64, f32)> {
        peaks.push((mz, inten));
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        peaks
    }

    fn profile_of(peaks: &[(f64, f32)]) -> OxoniumProfile {
        let stats = SpectrumStats::new(peaks);
        oxonium_profile(peaks, &stats, 20.0)
    }

    #[test]
    fn profile_normalises_by_base_peak_and_maxes_within_class() {
        let peaks = core_peaks();
        let p = profile_of(&peaks);
        // base peak is 204.08665 at 1000 -> HexNAc class saturates at 1.0
        assert!((p.get(MonoClass::HexNAc) - 1.0).abs() < 1e-6);
        assert!((p.get(MonoClass::HexHexNAc) - 0.6).abs() < 1e-6);
        assert_eq!(p.get(MonoClass::NeuAc), 0.0);
        assert_eq!(p.get(MonoClass::NeuGc), 0.0);
        assert_eq!(p.get(MonoClass::Fuc), 0.0);
    }

    #[test]
    fn neuac_657_counts_toward_the_neuac_class() {
        let peaks = with_ion(core_peaks(), 657.23488, 300.0);
        let p = profile_of(&peaks);
        assert!((p.get(MonoClass::NeuAc) - 0.3).abs() < 1e-6);
    }

    #[test]
    fn claiming_neuac_with_no_neuac_oxonium_is_penalised() {
        // THE missing penalty: `sialic_consistency` returns -0.0 == 0.0 here.
        let p = profile_of(&core_peaks());
        assert_eq!(p.get(MonoClass::NeuAc), 0.0);
        let sialylated = comp(4, 5, 0, 2, 0);
        let bare = comp(4, 5, 0, 0, 0);
        let params = OxoniumLlrParams::default();
        let i = MonoClass::NeuAc.idx();
        let logit = params.class_logit(MonoClass::NeuAc, 0.0);
        assert!(logit < 0.0, "absent ion must give a negative logit");
        let claimed = params.w_present[i] * logit;
        assert!(claimed < 0.0, "claiming NeuAc without evidence must cost");
        // and the whole-score ordering follows
        assert!(oxonium_composition_llr(&p, &sialylated) < oxonium_composition_llr(&p, &bare));
    }

    #[test]
    fn not_claiming_neuac_when_oxonium_is_strong_is_penalised() {
        let peaks = with_ion(core_peaks(), 274.09211, 800.0);
        let p = profile_of(&peaks);
        assert!(p.get(MonoClass::NeuAc) > 0.5);
        let params = OxoniumLlrParams::default();
        let i = MonoClass::NeuAc.idx();
        let logit = params.class_logit(MonoClass::NeuAc, p.get(MonoClass::NeuAc));
        let unclaimed = params.w_absent[i] * -logit;
        assert!(unclaimed < 0.0, "ignoring a strong NeuAc ion must cost");
        assert!(
            oxonium_composition_llr(&p, &comp(4, 5, 0, 0, 0))
                < oxonium_composition_llr(&p, &comp(4, 5, 0, 2, 0))
        );
    }

    #[test]
    fn matching_cases_are_positive() {
        // claims NeuAc, NeuAc ion present -> positive
        let peaks = with_ion(core_peaks(), 292.10267, 700.0);
        let p = profile_of(&peaks);
        assert!(oxonium_composition_llr(&p, &comp(4, 5, 0, 1, 0)) > 0.0);
        // claims no NeuAc/NeuGc/Fuc, none present -> positive
        let p0 = profile_of(&core_peaks());
        assert!(oxonium_composition_llr(&p0, &comp(4, 5, 0, 0, 0)) > 0.0);
    }

    #[test]
    fn fuc_absence_is_weaker_evidence_than_neuac_absence() {
        // Asymmetry rationale is documented on OxoniumLlrParams::default.
        let params = OxoniumLlrParams::default();
        let f = MonoClass::Fuc.idx();
        let a = MonoClass::NeuAc.idx();
        assert!(params.w_absent[f] < params.w_absent[a]);
        // presence side is NOT discounted for fucose
        assert!(params.w_present[f] >= params.w_absent[a]);

        // and the behaviour follows: claiming Fuc without its ion costs strictly
        // less than claiming NeuAc without its ion, at equal (zero) observation.
        let p = profile_of(&core_peaks());
        let base = oxonium_composition_llr(&p, &comp(4, 5, 0, 0, 0));
        let fuc_cost = base - oxonium_composition_llr(&p, &comp(4, 5, 1, 0, 0));
        let sia_cost = base - oxonium_composition_llr(&p, &comp(4, 5, 0, 1, 0));
        assert!(fuc_cost > 0.0 && sia_cost > 0.0);
        assert!(fuc_cost < sia_cost);
    }

    #[test]
    fn fuc_presence_is_strong_evidence() {
        let peaks = with_ion(core_peaks(), 512.19738, 500.0);
        let p = profile_of(&peaks);
        assert!(
            oxonium_composition_llr(&p, &comp(4, 5, 1, 0, 0))
                > oxonium_composition_llr(&p, &comp(4, 5, 0, 0, 0))
        );
    }

    #[test]
    fn core_classes_carry_no_discriminative_weight() {
        let params = OxoniumLlrParams::default();
        for class in [MonoClass::HexNAc, MonoClass::HexHexNAc] {
            let i = class.idx();
            assert_eq!(params.w_present[i], 0.0);
            assert_eq!(params.w_absent[i], 0.0);
        }
        // Behavioural check: wiping the core ions out of the spectrum entirely must
        // not move the score by a single bit.
        let with_core = profile_of(&core_peaks());
        let without_core = profile_of(&[(500.0, 50.0)]);
        let c = comp(4, 5, 0, 0, 0);
        assert_eq!(
            oxonium_composition_llr(&with_core, &c),
            oxonium_composition_llr(&without_core, &c)
        );
        // ...and neither can a composition with vs without Hex, which only changes
        // the HexHexNAc claim.
        assert_eq!(
            oxonium_composition_llr(&with_core, &comp(4, 5, 0, 0, 0)),
            oxonium_composition_llr(&with_core, &comp(4, 0, 0, 0, 0))
        );
    }

    #[test]
    fn neugc_is_scored_independently_of_neuac() {
        let peaks = with_ion(core_peaks(), 308.09759, 600.0);
        let p = profile_of(&peaks);
        assert!(p.get(MonoClass::NeuGc) > 0.5);
        assert_eq!(p.get(MonoClass::NeuAc), 0.0);
        assert!(
            oxonium_composition_llr(&p, &comp(4, 5, 0, 0, 1))
                > oxonium_composition_llr(&p, &comp(4, 5, 0, 1, 0))
        );
    }

    #[test]
    fn zero_observation_is_finite_and_clamped() {
        let params = OxoniumLlrParams::default();
        let empty = OxoniumProfile::default();
        for class in MONO_CLASSES {
            let l = params.class_logit(class, 0.0);
            assert!(l.is_finite(), "{class:?} logit not finite at o == 0");
            assert!(l.abs() <= params.clamp + 1e-6);
        }
        for c in [
            comp(4, 5, 0, 0, 0),
            comp(6, 7, 1, 3, 1),
            comp(2, 3, 0, 0, 0),
        ] {
            assert!(oxonium_composition_llr(&empty, &c).is_finite());
        }
        // negative intensity (should never happen, but must not produce NaN)
        let mut weird = OxoniumProfile::default();
        weird.obs[MonoClass::NeuAc.idx()] = -1.0;
        assert!(oxonium_composition_llr(&weird, &comp(4, 5, 0, 1, 0)).is_finite());
    }

    #[test]
    fn saturated_and_absent_scores_are_bounded_by_the_clamp() {
        let params = OxoniumLlrParams::default();
        let bound: f32 = params.clamp
            * MONO_CLASSES
                .iter()
                .map(|c| params.w_present[c.idx()].max(params.w_absent[c.idx()]))
                .sum::<f32>();
        let full = OxoniumProfile { obs: [1.0; 5] };
        let empty = OxoniumProfile::default();
        for p in [full, empty] {
            for c in [comp(4, 5, 1, 2, 1), comp(2, 0, 0, 0, 0)] {
                assert!(oxonium_composition_llr(&p, &c).abs() <= bound + 1e-4);
            }
        }
    }

    #[test]
    fn deterministic_across_repeated_and_reordered_input() {
        let peaks = with_ion(with_ion(core_peaks(), 274.09211, 300.0), 350.14456, 200.0);
        let c = comp(4, 5, 1, 1, 0);
        let a = oxonium_composition_llr(&profile_of(&peaks), &c);
        let b = oxonium_composition_llr(&profile_of(&peaks), &c);
        assert_eq!(a, b);
        let mut rev = peaks.clone();
        rev.reverse();
        assert_eq!(a, oxonium_composition_llr(&profile_of(&rev), &c));
    }

    #[test]
    fn cofrag_downweight_shrinks_only_prevalent_classes() {
        let params = OxoniumLlrParams::default();
        let mut prev = ClassPrevalence::default();
        prev.frac[MonoClass::NeuAc.idx()] = 0.9; // the brain-tissue case
        prev.frac[MonoClass::Fuc.idx()] = 0.3; // unremarkable
        let adj = cofrag_downweight(&params, &prev);

        let a = MonoClass::NeuAc.idx();
        let f = MonoClass::Fuc.idx();
        assert!((adj.w_present[a] - params.w_present[a] * 0.2).abs() < 1e-6);
        assert!((adj.w_absent[a] - params.w_absent[a] * 0.2).abs() < 1e-6);
        assert_eq!(adj.w_present[f], params.w_present[f]);
        assert_eq!(adj.w_absent[f], params.w_absent[f]);

        // the corrupted class now moves the score less than it used to
        let peaks = with_ion(core_peaks(), 274.09211, 900.0);
        let p = profile_of(&peaks);
        let c = comp(4, 5, 0, 2, 0);
        assert!(
            oxonium_composition_llr_with(&p, &c, &adj).abs()
                < oxonium_composition_llr_with(&p, &c, &params).abs()
        );
    }

    #[test]
    fn cofrag_downweight_at_full_prevalence_zeroes_the_class() {
        let params = OxoniumLlrParams::default();
        let prev = ClassPrevalence { frac: [1.0; 5] };
        let adj = cofrag_downweight(&params, &prev);
        let peaks = with_ion(core_peaks(), 274.09211, 900.0);
        let p = profile_of(&peaks);
        assert_eq!(
            oxonium_composition_llr_with(&p, &comp(4, 5, 1, 2, 1), &adj),
            0.0
        );
    }

    #[test]
    fn llr_is_monotone_in_the_observed_intensity_for_a_claiming_composition() {
        let params = OxoniumLlrParams::default();
        let c = comp(4, 5, 0, 1, 0);
        let mut last = f32::NEG_INFINITY;
        for o in [0.0f32, 0.005, 0.02, 0.1, 0.5, 1.0] {
            let mut p = OxoniumProfile::default();
            p.obs[MonoClass::NeuAc.idx()] = o;
            let v = oxonium_composition_llr_with(&p, &c, &params);
            assert!(v >= last, "not monotone at o={o}");
            last = v;
        }
    }
}
