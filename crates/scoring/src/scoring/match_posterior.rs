//! Posterior probability that a matched peak was produced by the ion it was assigned to.
//!
//! A peak at `x` assigned to ion `a` is explained by `a` or by any other theoretical ion of the
//! peptide within the fragment tolerance, or by chance. With ion presence priors `pi`, a Gaussian
//! mass error of width `sigma` (a quarter of the run's fragment tolerance, the relation the
//! tolerance is learned with) and the local peak density `rho` per Da as the chance term:
//!
//! `P(a | x) = pi_a N(x - m_a; sigma) / (sum_k pi_k N(x - m_k; sigma) + rho(x))`.
//!
//! The ion presence priors come from the scoring model: one minus the probability mass its rank
//! distribution puts on "no peak within the rank limit" for that ion type in the partition of
//! the precursor charge, parent mass and fragment segment.

use crate::param_model::{IonType, Param};
use crate::scoring::fragment_error::TOL_SIGMAS;
use crate::scoring::fragment_ions::IonKind;
use crate::scoring::strong_score::{FragmentIon, FragmentLoss, MatchedIon, DENSITY_HW};
use crate::scoring::ScoredSpectrum;
use model::mass::{H2O, NH3, PROTON};

/// Largest difference in Da between a model ion offset and the offset expected for a series and
/// loss class, for the two to count as the same ion type.
const OFFSET_MATCH_DA: f64 = 0.01;

/// Prior used for every ion when the model carries no rank distributions for the partition.
const UNIFORM_PRIOR: f64 = 1.0;

/// Smallest ion presence prior, so an ion type the model never saw still competes.
const MIN_PRIOR: f64 = 1e-6;

/// Posteriors are kept at least this far from 0 and 1 before the log-odds.
const POSTERIOR_EPSILON: f64 = 1e-6;

/// Ion presence priors of one scored spectrum, read from the scoring model.
pub struct IonPriors<'a> {
    param: &'a Param,
    charge: u8,
    parent_mass: f64,
}

impl<'a> IonPriors<'a> {
    pub fn new(param: &'a Param, charge: u8, parent_mass: f64) -> Self {
        Self {
            param,
            charge,
            parent_mass,
        }
    }

    /// Probability that the theoretical ion `ion` at `mz` has a peak in the spectrum.
    pub fn prior(&self, ion: FragmentIon, mz: f64) -> f64 {
        let seg = self.param.segment_num(mz, self.parent_mass);
        let part = self.param.partition_for(self.charge, self.parent_mass, seg);
        let Some(table) = self.param.rank_dist_table.get(&part) else {
            return UNIFORM_PRIOR;
        };
        let presence = |dist: &Vec<f32>| {
            (1.0 - f64::from(dist.last().copied().unwrap_or(1.0))).clamp(0.0, 1.0)
        };
        let noise = table.get(&IonType::Noise).map(presence).unwrap_or(0.0);
        let want_prefix = matches!(ion.kind, IonKind::B | IonKind::C);
        let expected = expected_offset(ion);
        table
            .iter()
            .filter(|(t, _)| {
                !t.is_noise()
                    && t.charge() == Some(1)
                    && t.is_prefix() == want_prefix
                    && match ion.loss {
                        FragmentLoss::Modification(c) => t.loss_class() == c,
                        _ => t
                            .offset()
                            .is_some_and(|o| (f64::from(o) - expected).abs() < OFFSET_MATCH_DA),
                    }
            })
            .map(|(_, d)| presence(d))
            .fold(None, |best: Option<f64>, p| {
                Some(best.map_or(p, |b| b.max(p)))
            })
            .unwrap_or(noise)
            .max(MIN_PRIOR)
    }
}

fn expected_offset(ion: FragmentIon) -> f64 {
    let base = match ion.kind {
        IonKind::B | IonKind::C => PROTON,
        IonKind::Y | IonKind::Z => PROTON + H2O,
    };
    match ion.loss {
        FragmentLoss::Water => base - H2O,
        FragmentLoss::Ammonia => base - NH3,
        FragmentLoss::None | FragmentLoss::Modification(_) => base,
    }
}

/// Summary of the posteriors of the matched peaks of one PSM.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct MatchPosteriors {
    pub mean: f32,
    pub log_odds_sum: f32,
}

/// Posterior of every matched peak and its summary (mean posterior, sum of log-odds).
pub fn match_posteriors(
    scored_spec: &ScoredSpectrum<'_>,
    matched_ions: &[MatchedIon],
    theoretical: &[(f64, FragmentIon)],
    priors: &IonPriors<'_>,
    feature_tol: f64,
    feature_tol_is_ppm: bool,
) -> MatchPosteriors {
    if matched_ions.is_empty() {
        return MatchPosteriors::default();
    }
    let mut sum_p = 0.0_f64;
    let mut sum_logit = 0.0_f64;
    for &(_, obs, pred, _, ion) in matched_ions {
        let tol_da = if feature_tol_is_ppm {
            obs * feature_tol / 1e6
        } else {
            feature_tol
        };
        let sigma = tol_da / TOL_SIGMAS;
        let density = |theo: f64| {
            let z = (obs - theo) / sigma;
            (-0.5 * z * z).exp() / (sigma * (2.0 * std::f64::consts::PI).sqrt())
        };
        let own = priors.prior(ion, pred) * density(pred);
        let others: f64 = theoretical
            .iter()
            .filter(|&&(theo, other)| (theo - obs).abs() <= tol_da && other != ion)
            .map(|&(theo, other)| priors.prior(other, theo) * density(theo))
            .sum();
        let chance = scored_spec.local_peak_density(obs, DENSITY_HW);
        let denom = own + others + chance;
        let p = if denom > 0.0 { own / denom } else { 0.0 };
        let p = p.clamp(POSTERIOR_EPSILON, 1.0 - POSTERIOR_EPSILON);
        sum_p += p;
        sum_logit += (p / (1.0 - p)).ln();
    }
    MatchPosteriors {
        mean: (sum_p / matched_ions.len() as f64) as f32,
        log_odds_sum: sum_logit as f32,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::fragment_ions::predict_by_ions;
    use crate::scoring::strong_score::theoretical_fragment_ions;
    use model::amino_acid::AminoAcid;
    use model::peptide::Peptide;
    use model::spectrum::Spectrum;

    fn pep(seq: &[u8]) -> Peptide {
        let residues: Vec<AminoAcid> = seq
            .iter()
            .map(|&r| AminoAcid::standard(r).unwrap())
            .collect();
        Peptide::new(residues, b'K', b'R')
    }

    fn param_without_tables() -> Param {
        let mut param = crate::testutil::tiny_param();
        param.rank_dist_table.clear();
        param
    }

    fn param_with_b_and_noise() -> Param {
        use rustc_hash::FxHashMap;
        let mut param = crate::testutil::tiny_param();
        let part = param.partition_for(2, 1500.0, 0);
        let mut table: FxHashMap<IonType, Vec<f32>> = FxHashMap::default();
        table.insert(
            IonType::Prefix {
                charge: 1,
                offset_bits: (PROTON as f32).to_bits(),
                loss_class: 0,
            },
            vec![0.5, 0.3, 0.0, 0.2],
        );
        table.insert(IonType::Noise, vec![0.01, 0.01, 0.03, 0.95]);
        param.rank_dist_table.insert(part, table);
        param
    }

    #[test]
    fn priors_follow_the_model_and_unmodelled_ions_fall_back_to_noise() {
        let param = param_with_b_and_noise();
        let priors = IonPriors::new(&param, 2, 1500.0);
        let b = FragmentIon {
            kind: IonKind::B,
            position: 3,
            loss: FragmentLoss::None,
        };
        let y = FragmentIon {
            kind: IonKind::Y,
            position: 3,
            loss: FragmentLoss::None,
        };
        let b_water = FragmentIon {
            loss: FragmentLoss::Water,
            ..b
        };
        assert!((priors.prior(b, 400.0) - 0.8).abs() < 1e-6);
        assert!((priors.prior(y, 400.0) - 0.05).abs() < 1e-6);
        assert!((priors.prior(b_water, 400.0) - 0.05).abs() < 1e-6);
    }

    #[test]
    fn a_coincident_ion_halves_the_posterior_and_an_unambiguous_ion_keeps_it() {
        let param = param_without_tables();
        let priors = IonPriors::new(&param, 2, 1500.0);
        let peptide = pep(b"GAKAG");
        let predicted = predict_by_ions(&peptide, 1..=1);
        let theoretical = theoretical_fragment_ions(&predicted);
        let b2 = predicted
            .iter()
            .find(|p| p.kind == IonKind::B && p.position == 2)
            .unwrap();
        let y3 = predicted
            .iter()
            .find(|p| p.kind == IonKind::Y && p.position == 3)
            .unwrap();
        let spectrum = Spectrum {
            peaks: vec![(b2.mz, 1000.0), (y3.mz, 1000.0)],
            ..Default::default()
        };
        let ss = ScoredSpectrum::new_without_filtering(&spectrum);
        let coincident = [(
            1000.0_f32,
            b2.mz,
            b2.mz,
            true,
            FragmentIon::of_predicted(b2),
        )];
        let alone = [(
            1000.0_f32,
            y3.mz,
            y3.mz,
            false,
            FragmentIon::of_predicted(y3),
        )];
        let pc = match_posteriors(&ss, &coincident, &theoretical, &priors, 20.0, true);
        let pa = match_posteriors(&ss, &alone, &theoretical, &priors, 20.0, true);
        assert!((f64::from(pc.mean) - 0.5).abs() < 0.01, "{}", pc.mean);
        assert!(pa.mean > 0.99, "{}", pa.mean);
        assert!(pa.log_odds_sum > pc.log_odds_sum);
    }

    #[test]
    fn a_higher_prior_competitor_lowers_the_posterior() {
        let param = param_with_b_and_noise();
        let priors = IonPriors::new(&param, 2, 1500.0);
        let peptide = pep(b"GAKAG");
        let predicted = predict_by_ions(&peptide, 1..=1);
        let theoretical = theoretical_fragment_ions(&predicted);
        let b2 = predicted
            .iter()
            .find(|p| p.kind == IonKind::B && p.position == 2)
            .unwrap();
        let y2 = predicted
            .iter()
            .find(|p| p.kind == IonKind::Y && p.position == 2)
            .unwrap();
        let spectrum = Spectrum {
            peaks: vec![(b2.mz, 1000.0)],
            ..Default::default()
        };
        let ss = ScoredSpectrum::new_without_filtering(&spectrum);
        let as_b = [(
            1000.0_f32,
            b2.mz,
            b2.mz,
            true,
            FragmentIon::of_predicted(b2),
        )];
        let water_loss = FragmentIon {
            loss: FragmentLoss::Water,
            ..FragmentIon::of_predicted(y2)
        };
        let as_y_loss = [(1000.0_f32, b2.mz, b2.mz, false, water_loss)];
        let pb = match_posteriors(&ss, &as_b, &theoretical, &priors, 20.0, true);
        let py = match_posteriors(&ss, &as_y_loss, &theoretical, &priors, 20.0, true);
        assert!(
            (f64::from(pb.mean) - 0.8 / 0.85).abs() < 0.01,
            "{}",
            pb.mean
        );
        assert!(
            (f64::from(py.mean) - 0.05 / 0.85).abs() < 0.01,
            "{}",
            py.mean
        );
    }
}
