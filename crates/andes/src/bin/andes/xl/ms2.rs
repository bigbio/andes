//! Fragment-intensity evidence of the two peptides of a crosslink spectrum match: how well the
//! intensities of each peptide's fragments that do not carry the link agree with the bundled
//! andes-ms2 prediction of that peptide as a linear one.
//!
//! A fragment without the link site has the mass of the linear peptide's fragment, so the
//! linear prediction applies to it; a fragment with the link site carries the partner and is
//! left out. A random partner that matches a few fragments by chance matches them without the
//! peptide's intensity pattern.

use std::borrow::Cow;

use model::{Peptide, Spectrum};
use onnx_infer::predictors::form::parse_peptide;
use onnx_infer::predictors::spectral::{ms2_features_where, FragmentTolerance};
use onnx_infer::predictors::{NativeMs2, Predictor};
use rayon::prelude::*;
use search::crosslinks::{Csm, PeptideHit, XlCandidates, XlSpectrumResult};

use crate::acquisition::Acquisition;
use crate::cli::{Cli, PredictFeature};
use crate::xl::fdr::XlRun;

/// The PIN columns `--predict ms2` adds to a crosslink PIN.
pub(crate) const FEATURES: [&str; 7] = [
    "alpha_ms2_angle",
    "beta_ms2_angle",
    "min_ms2_angle",
    "alpha_ms2_pearson",
    "beta_ms2_pearson",
    "min_ms2_pearson",
    "ms2_missing",
];

/// Indices of the spectral angle and the Pearson correlation in `ms2_features_where`.
const ANGLE: usize = 6;
const PEARSON: usize = 0;

pub(crate) type Ms2Row = [f64; FEATURES.len()];

pub(crate) struct XlMs2 {
    model: NativeMs2,
    tol: FragmentTolerance,
}

impl XlMs2 {
    /// The model when `--predict` asks for `ms2`, else `None`.
    pub(crate) fn from_cli(cli: &Cli, acquisition: &Acquisition) -> Result<Option<Self>, String> {
        if !cli.predict.contains(&PredictFeature::Ms2) {
            return Ok(None);
        }
        Ok(Some(XlMs2 {
            model: crate::dia::rt::andes_ms2(cli, acquisition)?,
            tol: crate::dia::rt::feature_tolerance(cli, acquisition.high_res),
        }))
    }

    /// One feature row per result with a CSM, aligned with `results`.
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
                Some(self.row(&spectra[r.spec_idx], c, &runs[r.linker].table))
            })
            .collect()
    }

    fn row(&self, spec: &Spectrum, c: &Csm, table: &XlCandidates) -> Ms2Row {
        let peaks = sorted(&spec.peaks);
        let z = predictor_charge(c.charge);
        let side = |h: &PeptideHit| {
            self.side(&table.peptides()[h.pep as usize].peptide, h.site.pos, z, &peaks)
        };
        let (a, b) = (side(&c.alpha), side(&c.beta));
        let missing = u8::from(a.is_none() || b.is_none()) as f64;
        let get = |s: &Option<[f64; 8]>, k: usize| s.map_or(0.0, |v| v[k]);
        let (aa, ba) = (get(&a, ANGLE), get(&b, ANGLE));
        let (ap, bp) = (get(&a, PEARSON), get(&b, PEARSON));
        [aa, ba, aa.min(ba), ap, bp, ap.min(bp), missing]
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
    fn each_peptide_of_a_crosslink_is_predicted_at_a_partial_charge() {
        assert_eq!(predictor_charge(3), 2);
        assert_eq!(predictor_charge(4), 2);
        assert_eq!(predictor_charge(5), 3);
        assert_eq!(predictor_charge(7), 3);
    }
}
