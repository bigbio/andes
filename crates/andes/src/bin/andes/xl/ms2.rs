//! Fragment-intensity evidence of the two peptides of a crosslink spectrum match: how well the
//! intensities of each peptide's fragments that do not carry the link agree with the bundled
//! andes-ms2 prediction of that peptide as a linear one.
//!
//! A fragment without the link site has the mass of the linear peptide's fragment, so the
//! linear prediction applies to it; a fragment with the link site carries the partner and is
//! left out. A random partner that matches a few fragments by chance matches them without the
//! peptide's intensity pattern.
//!
//! The features apply to MS-cleavable crosslinkers only. For a non-cleavable linker, a partner
//! that is in the spectrum but not linked to the other peptide has a genuine intensity pattern
//! too, so the features reward it, and the decoys, which never have one, cannot model that.

use std::borrow::Cow;

use model::{Peptide, Spectrum};
use onnx_infer::predictors::form::parse_peptide;
use onnx_infer::predictors::spectral::{ms2_features_where, FragmentTolerance};
use onnx_infer::predictors::{NativeMs2, Predictor};
use rayon::prelude::*;
use search::crosslinks::{Csm, PeptideHit, XlCandidates, XlSpectrumResult};
use xl::crosslinker::Crosslinker;

use crate::acquisition::Acquisition;
use crate::cli::{Cli, PredictFeature};
use crate::xl::fdr::XlRun;

/// The PIN columns `--predict ms2` adds to a crosslink PIN.
pub(crate) const FEATURES: [&str; 4] = [
    "alpha_ms2_angle",
    "beta_ms2_angle",
    "min_ms2_angle",
    "ms2_missing",
];

/// Index of the spectral angle in `ms2_features_where`.
const ANGLE: usize = 6;

pub(crate) type Ms2Row = [f64; FEATURES.len()];

/// The features of a CSM without a prediction (a non-cleavable linker's CSM, or a peptide the
/// model cannot encode): zero, marked missing.
pub(crate) const MISSING: Ms2Row = [0.0, 0.0, 0.0, 1.0];

/// Whether the features apply to a search with these linkers: one of them is MS-cleavable.
pub(crate) fn applies(linkers: &[Crosslinker]) -> bool {
    linkers.iter().any(|l| l.cleavage.is_some())
}

pub(crate) struct XlMs2 {
    model: NativeMs2,
    tol: FragmentTolerance,
}

impl XlMs2 {
    /// The model when `--predict` asks for `ms2` and one of `linkers` is MS-cleavable, else
    /// `None`.
    pub(crate) fn from_cli(
        cli: &Cli,
        acquisition: &Acquisition,
        linkers: &[Crosslinker],
    ) -> Result<Option<Self>, String> {
        if !cli.predict.contains(&PredictFeature::Ms2) {
            return Ok(None);
        }
        if !applies(linkers) {
            log::info!(
                target: "xl",
                "--predict ms2: the fragment-intensity features apply to MS-cleavable \
                 crosslinkers only; none of this search's linkers is, so the crosslink PIN has none"
            );
            return Ok(None);
        }
        Ok(Some(XlMs2 {
            model: crate::dia::rt::andes_ms2(cli, acquisition)?,
            tol: crate::dia::rt::feature_tolerance(cli, acquisition.high_res),
        }))
    }

    /// One feature row per result with a CSM of an MS-cleavable linker, aligned with `results`.
    pub(crate) fn rows(
        &self,
        spectra: &[Spectrum],
        results: &[XlSpectrumResult],
        runs: &[XlRun],
    ) -> Vec<Option<Ms2Row>> {
        results
            .par_iter()
            .map(|r| {
                let c = r.csm.as_ref()?;
                let run = &runs[r.linker];
                run.linker.cleavage.as_ref()?;
                Some(self.row(&spectra[r.spec_idx], c, &run.table))
            })
            .collect()
    }

    fn row(&self, spec: &Spectrum, c: &Csm, table: &XlCandidates) -> Ms2Row {
        let peaks = sorted(&spec.peaks);
        let z = predictor_charge(c.charge);
        let side = |h: &PeptideHit| {
            self.side(
                &table.peptides()[h.pep as usize].peptide,
                h.site.pos,
                z,
                &peaks,
            )
        };
        let (a, b) = (side(&c.alpha), side(&c.beta));
        let missing = u8::from(a.is_none() || b.is_none()) as f64;
        let get = |s: &Option<[f64; 8]>, k: usize| s.map_or(0.0, |v| v[k]);
        let (aa, ba) = (get(&a, ANGLE), get(&b, ANGLE));
        [aa, ba, aa.min(ba), missing]
    }

    fn side(&self, p: &Peptide, site: usize, z: u8, peaks: &[(f64, f32)]) -> Option<[f64; 8]> {
        let (form, mass) = parse_peptide(&p.to_string())?;
        if !self.model.supports_ms2(&form, z) {
            return None;
        }
        let pred = self.model.ms2(&mut [], &[(&form, z)]).ok()?.pop()?;
        Some(ms2_features_where(
            &pred,
            self.model.ions(),
            &mass,
            z,
            peaks,
            self.tol,
            &|i, b| unlinked(i, b, site),
        ))
    }
}

/// Whether the fragment at cleavage index `i` (the b ion holds residues `0..=i`) leaves out
/// the link site at residue index `site`.
fn unlinked(i: usize, b_ion: bool, site: usize) -> bool {
    if b_ion {
        i < site
    } else {
        i >= site
    }
}

/// The precursor charge the linear model is asked for. A crosslinked precursor splits its
/// charge over two peptides, so each is predicted at 2, or at 3 from a charge of 5 up.
fn predictor_charge(z: u8) -> u8 {
    if z >= 5 {
        3
    } else {
        2
    }
}

fn sorted(peaks: &[(f64, f32)]) -> Cow<'_, [(f64, f32)]> {
    if peaks.windows(2).all(|w| w[0].0 <= w[1].0) {
        Cow::Borrowed(peaks)
    } else {
        let mut p = peaks.to_vec();
        p.sort_by(|a, b| a.0.total_cmp(&b.0));
        Cow::Owned(p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_fragment_holding_the_link_site_is_left_out() {
        assert!(unlinked(1, true, 2));
        assert!(!unlinked(2, true, 2));
        assert!(unlinked(2, false, 2));
        assert!(!unlinked(1, false, 2));
        assert!(!unlinked(0, true, 0));
        assert!(unlinked(0, false, 0));
    }

    #[test]
    fn the_features_apply_only_to_a_search_with_an_ms_cleavable_linker() {
        let dsso = Crosslinker::dsso();
        let mut dss = dsso.clone();
        dss.cleavage = None;
        assert!(applies(&[dsso.clone()]));
        assert!(!applies(&[dss.clone()]));
        assert!(applies(&[dss.clone(), dsso]));
        assert!(!applies(&[]));
        assert_eq!(MISSING.len(), FEATURES.len());
        assert_eq!(MISSING[FEATURES.len() - 1], 1.0);
    }

    #[test]
    fn each_peptide_of_a_crosslink_is_predicted_at_a_partial_charge() {
        assert_eq!(predictor_charge(3), 2);
        assert_eq!(predictor_charge(4), 2);
        assert_eq!(predictor_charge(5), 3);
        assert_eq!(predictor_charge(7), 3);
    }
}
