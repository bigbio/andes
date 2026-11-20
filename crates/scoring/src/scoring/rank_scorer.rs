use std::collections::HashMap;

use crate::param_model::{IonType, Param, Partition};

#[derive(Debug, Clone)]
pub struct RankScorer {
    param: Param,
    pub(crate) log_table: HashMap<(Partition, IonType), Vec<f32>>,
    pub(crate) partition_ion_logs: HashMap<Partition, Vec<(IonType, Vec<f32>)>>,
    pub(crate) partition_loss_ion_logs: HashMap<Partition, Vec<(IonType, Vec<f32>)>>,
    max_rank: u32,
    fragment_tol_override: Option<model::tolerance::Tolerance>,
    fragment_error: crate::scoring::fragment_error::SharedFragmentError,
    acquisition_high_res: Option<bool>,
    recentre_default_tolerance: bool,
}

const FREQ_FLOOR: f32 = 1e-9;

#[inline]
fn finite_llr(ion_freq: f32, noise_freq: f32, norm: f32) -> f32 {
    let num = if ion_freq.is_finite() && ion_freq > 0.0 {
        ion_freq
    } else {
        FREQ_FLOOR
    };
    let den_raw = noise_freq * norm;
    let den = if den_raw.is_finite() && den_raw > 0.0 {
        den_raw
    } else {
        FREQ_FLOOR
    };
    let llr = (num / den).ln();
    debug_assert!(
        llr.is_finite(),
        "rank LLR must be finite (ion={ion_freq}, noise={noise_freq}, norm={norm})"
    );
    llr
}

impl RankScorer {
    pub fn new(param: &Param) -> Self {
        let mut log_table: HashMap<(Partition, IonType), Vec<f32>> = HashMap::new();

        for (partition, ion_table) in &param.rank_dist_table {
            let Some(noise_freqs) = ion_table.get(&IonType::Noise) else {
                continue;
            };

            for (ion_type, ion_freqs) in ion_table {
                let ion_charge = match ion_type {
                    IonType::Prefix { charge, .. } | IonType::Suffix { charge, .. } => *charge,
                    IonType::Noise => continue,
                };

                let charge_or_seg = (ion_charge as u32).min(param.num_segments as u32) as f32;

                let rank_count = ion_freqs.len().min(noise_freqs.len());
                let table: Vec<f32> = (0..rank_count)
                    .map(|r| finite_llr(ion_freqs[r], noise_freqs[r], charge_or_seg))
                    .collect();

                log_table.insert((*partition, *ion_type), table);
            }
        }

        let mut partition_ion_logs: HashMap<Partition, Vec<(IonType, Vec<f32>)>> = HashMap::new();
        for (&partition, ions) in &param.partition_ion_types_cache {
            let paired: Vec<(IonType, Vec<f32>)> = ions
                .iter()
                .filter(|ion| !ion.is_loss())
                .filter_map(|&ion| {
                    log_table
                        .get(&(partition, ion))
                        .map(|table| (ion, table.clone()))
                })
                .collect();
            partition_ion_logs.insert(partition, paired);
        }

        let mut partition_loss_ion_logs: HashMap<Partition, Vec<(IonType, Vec<f32>)>> =
            HashMap::new();
        for (&(partition, ion), table) in &log_table {
            if ion.is_loss() {
                partition_loss_ion_logs
                    .entry(partition)
                    .or_default()
                    .push((ion, table.clone()));
            }
        }
        for logs in partition_loss_ion_logs.values_mut() {
            logs.sort_by_key(|(ion, _)| match *ion {
                IonType::Prefix {
                    charge,
                    offset_bits,
                    loss_class,
                } => (0u8, charge, offset_bits, loss_class),
                IonType::Suffix {
                    charge,
                    offset_bits,
                    loss_class,
                } => (1u8, charge, offset_bits, loss_class),
                IonType::Noise => (2u8, 0, 0, 0),
            });
        }

        Self {
            param: param.clone(),
            log_table,
            partition_ion_logs,
            partition_loss_ion_logs,
            max_rank: param.max_rank as u32,
            fragment_tol_override: None,
            acquisition_high_res: None,
            fragment_error: crate::scoring::fragment_error::SharedFragmentError::new(),
            recentre_default_tolerance: false,
        }
    }

    pub fn partition_ion_logs(&self, partition: &Partition) -> &[(IonType, Vec<f32>)] {
        self.partition_ion_logs
            .get(partition)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    pub fn partition_loss_ion_logs(&self, partition: &Partition) -> &[(IonType, Vec<f32>)] {
        self.partition_loss_ion_logs
            .get(partition)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    pub fn has_loss_tables(&self) -> bool {
        !self.partition_loss_ion_logs.is_empty()
    }

    pub fn max_rank(&self) -> u32 {
        self.max_rank
    }

    pub fn param(&self) -> &Param {
        &self.param
    }

    pub fn set_fragment_tol_override(&mut self, tol: Option<model::tolerance::Tolerance>) {
        self.fragment_tol_override = tol;
    }

    /// Recentres the fragment tolerance override on the run's measured offset when `on`; the
    /// override then stands for a `--dia` default, not for the user's `--fragment-tol-*`.
    pub fn set_recentre_default_tolerance(&mut self, on: bool) {
        self.recentre_default_tolerance = on;
    }

    pub fn recentres_default_tolerance(&self) -> bool {
        self.recentre_default_tolerance
    }

    pub fn set_acquisition_high_res(&mut self, high_res: Option<bool>) {
        self.acquisition_high_res = high_res;
    }

    pub fn acquisition_is_high_res(&self) -> bool {
        self.acquisition_high_res
            .unwrap_or_else(|| self.param.data_type.instrument.is_high_resolution())
    }

    pub fn feature_match_tolerance(&self) -> model::tolerance::Tolerance {
        use model::tolerance::Tolerance;
        self.fragment_tol_override.unwrap_or_else(|| {
            if self.acquisition_is_high_res() {
                Tolerance::Ppm(20.0)
            } else {
                Tolerance::Da(0.5)
            }
        })
    }

    pub fn set_fragment_error(&self, m: Option<crate::scoring::FragmentErrorModel>) {
        self.fragment_error.set(m);
    }

    pub fn fragment_error(&self) -> Option<crate::scoring::FragmentErrorModel> {
        self.fragment_error.get()
    }

    pub fn pin_match_tolerance(&self) -> model::tolerance::Tolerance {
        if let Some(t) = self.fragment_tol_override {
            return t;
        }
        match self.fragment_error.get() {
            Some(m) if self.acquisition_is_high_res() => {
                model::tolerance::Tolerance::Ppm(m.feature_tolerance_ppm())
            }
            _ => self.feature_match_tolerance(),
        }
    }

    pub fn node_score_window(&self) -> model::tolerance::Tolerance {
        use model::tolerance::Tolerance;
        let stored = self.param.mme;
        if self.fragment_tol_override.is_some() || !self.acquisition_is_high_res() {
            return stored;
        }
        match (
            stored,
            self.fragment_error.get().and_then(|m| m.node_window_da()),
        ) {
            (Tolerance::Da(v), Some(w)) => Tolerance::Da(w.min(v)),
            _ => stored,
        }
    }

    pub fn fragment_recentre_ppm(&self) -> f64 {
        match self.fragment_error.get() {
            Some(m)
                if (self.fragment_tol_override.is_none() || self.recentre_default_tolerance)
                    && self.acquisition_is_high_res() =>
            {
                m.offset_ppm
            }
            _ => 0.0,
        }
    }

    pub fn ppm_gaussian(&self) -> (f64, f64) {
        match self.fragment_error.get() {
            Some(m) if self.fragment_tol_override.is_none() && self.acquisition_is_high_res() => {
                (0.0, m.gaussian_sigma_ppm())
            }
            _ => (
                0.0,
                crate::scoring::fragment_error::DEFAULT_GAUSSIAN_SIGMA_PPM,
            ),
        }
    }

    pub fn node_score(&self, partition: Partition, ion_type: IonType, rank: u32) -> f32 {
        let Some(table) = self.log_table.get(&(partition, ion_type)) else {
            return 0.0;
        };
        if self.max_rank == 0 {
            return 0.0;
        }
        let observed_rank = rank.clamp(1, self.max_rank);
        table
            .get((observed_rank - 1) as usize)
            .copied()
            .unwrap_or(0.0)
    }

    pub fn missing_ion_score(&self, partition: Partition, ion_type: IonType) -> f32 {
        let Some(table) = self.log_table.get(&(partition, ion_type)) else {
            return 0.0;
        };
        table.get(self.max_rank as usize).copied().unwrap_or(0.0)
    }

    pub fn ion_existence_score(&self, partition: Partition, index: usize, prob_peak: f32) -> f32 {
        let Some(table) = self.param.ion_existence_table.get(&partition) else {
            return 0.0;
        };
        if index >= table.len() {
            return 0.0;
        }
        if !prob_peak.is_finite() || !(0.0..=1.0).contains(&prob_peak) {
            return 0.0;
        }
        let noise_baseline = match index {
            0 => (1.0 - prob_peak) * (1.0 - prob_peak),
            3 => prob_peak * prob_peak,
            _ => prob_peak * (1.0 - prob_peak),
        };
        let ion_prob = if table[index] == 0.0 {
            0.01
        } else {
            table[index]
        };
        (ion_prob / noise_baseline).ln()
    }

    pub fn error_score(&self, partition: Partition, error: f32) -> f32 {
        let esf = self.param.error_scaling_factor;
        if esf == 0 {
            return 0.0;
        }
        let mut err_index = (error * esf as f32).round() as i32;
        if err_index > esf {
            err_index = esf;
        } else if err_index < -esf {
            err_index = -esf;
        }
        err_index += esf;
        let idx = err_index as usize;

        let ion_err = match self.param.ion_err_dist_table.get(&partition) {
            Some(v) => v,
            None => return 0.0,
        };
        let noise_err = match self.param.noise_err_dist_table.get(&partition) {
            Some(v) => v,
            None => return 0.0,
        };
        if idx >= ion_err.len() || idx >= noise_err.len() {
            return 0.0;
        }
        let ion_f = ion_err[idx];
        let noise_f = noise_err[idx];
        if ion_f <= 0.0 || noise_f <= 0.0 {
            return 0.0;
        }
        (ion_f / noise_f).ln()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::tiny_param;

    fn low_res_param() -> Param {
        use model::instrument::InstrumentType;
        let mut p = tiny_param();
        p.data_type.instrument = InstrumentType::LowRes;
        p
    }

    #[test]
    fn feature_match_tolerance_defaults_by_resolution() {
        use model::tolerance::Tolerance;
        let hi = RankScorer::new(&tiny_param());
        assert_eq!(hi.feature_match_tolerance(), Tolerance::Ppm(20.0));
        let lo = RankScorer::new(&low_res_param());
        assert_eq!(lo.feature_match_tolerance(), Tolerance::Da(0.5));
    }

    #[test]
    fn feature_match_tolerance_follows_the_acquisition() {
        use model::tolerance::Tolerance;
        let mut lo_model = RankScorer::new(&low_res_param());
        lo_model.set_acquisition_high_res(Some(true));
        assert!(lo_model.acquisition_is_high_res());
        assert_eq!(lo_model.feature_match_tolerance(), Tolerance::Ppm(20.0));
        let mut hi_model = RankScorer::new(&tiny_param());
        hi_model.set_acquisition_high_res(Some(false));
        assert_eq!(hi_model.feature_match_tolerance(), Tolerance::Da(0.5));
        hi_model.set_acquisition_high_res(None);
        assert_eq!(hi_model.feature_match_tolerance(), Tolerance::Ppm(20.0));
        lo_model.set_fragment_tol_override(Some(Tolerance::Ppm(7.0)));
        assert_eq!(lo_model.feature_match_tolerance(), Tolerance::Ppm(7.0));
    }

    #[test]
    fn measured_error_sets_pin_tolerance_on_high_res_only() {
        use crate::scoring::FragmentErrorModel;
        use model::tolerance::Tolerance;
        let hi = RankScorer::new(&tiny_param());
        assert_eq!(hi.pin_match_tolerance(), Tolerance::Ppm(20.0));
        assert_eq!(hi.ppm_gaussian(), (0.0, 7.0));
        hi.set_fragment_error(Some(FragmentErrorModel::new(-2.5, 1.6, 100, 100)));
        let Tolerance::Ppm(t) = hi.pin_match_tolerance() else {
            panic!("ppm expected");
        };
        assert!((t - 6.4).abs() < 1e-9, "{t}");
        assert_eq!(hi.feature_match_tolerance(), Tolerance::Ppm(20.0));
        let (c, s) = hi.ppm_gaussian();
        assert!(c == 0.0 && (s - 2.4).abs() < 1e-9);
        assert_eq!(hi.fragment_recentre_ppm(), -2.5);
        let lo = RankScorer::new(&low_res_param());
        lo.set_fragment_error(Some(FragmentErrorModel::new(-2.5, 1.6, 100, 100)));
        assert_eq!(lo.pin_match_tolerance(), Tolerance::Da(0.5));
        assert_eq!(lo.ppm_gaussian(), (0.0, 7.0));
        assert_eq!(lo.fragment_recentre_ppm(), 0.0);
        let mut ov = RankScorer::new(&tiny_param());
        ov.set_fragment_error(Some(FragmentErrorModel::new(-2.5, 1.6, 100, 100)));
        ov.set_fragment_tol_override(Some(Tolerance::Ppm(15.0)));
        assert_eq!(ov.pin_match_tolerance(), Tolerance::Ppm(15.0));
        assert_eq!(ov.ppm_gaussian(), (0.0, 7.0));
        assert_eq!(ov.fragment_recentre_ppm(), 0.0);
        hi.set_fragment_error(None);
        assert_eq!(hi.fragment_recentre_ppm(), 0.0);
        assert_eq!(hi.pin_match_tolerance(), Tolerance::Ppm(20.0));
    }

    #[test]
    fn the_default_override_is_recentred_only_when_asked() {
        use crate::scoring::FragmentErrorModel;
        use model::tolerance::Tolerance;
        let m = FragmentErrorModel::new(-3.3, 1.3, 1000, 200);
        let mut s = RankScorer::new(&tiny_param());
        s.set_fragment_tol_override(Some(Tolerance::Ppm(7.0)));
        s.set_fragment_error(Some(m));
        assert_eq!(s.fragment_recentre_ppm(), 0.0, "off: the override wins");
        s.set_recentre_default_tolerance(true);
        assert_eq!(s.fragment_recentre_ppm(), -3.3);
        assert_eq!(s.feature_match_tolerance(), Tolerance::Ppm(7.0));
        assert_eq!(s.pin_match_tolerance(), Tolerance::Ppm(7.0));
        assert_eq!(s.ppm_gaussian(), (0.0, 7.0));
        s.set_fragment_error(None);
        assert_eq!(s.fragment_recentre_ppm(), 0.0, "unmeasured: the default");
        let mut lo = RankScorer::new(&low_res_param());
        lo.set_fragment_tol_override(Some(Tolerance::Da(0.5)));
        lo.set_recentre_default_tolerance(true);
        lo.set_fragment_error(Some(m));
        assert_eq!(lo.fragment_recentre_ppm(), 0.0);
        assert_eq!(lo.feature_match_tolerance(), Tolerance::Da(0.5));
    }

    #[test]
    fn node_score_window_is_learned_on_high_res_only() {
        use crate::scoring::FragmentErrorModel;
        use model::tolerance::Tolerance;
        let mut p = tiny_param();
        p.mme = Tolerance::Da(0.5);
        let hi = RankScorer::new(&p);
        assert_eq!(hi.node_score_window(), Tolerance::Da(0.5));
        hi.set_fragment_error(Some(FragmentErrorModel::new(-2.5, 1.6, 100, 100)));
        assert_eq!(
            hi.node_score_window(),
            Tolerance::Da(0.5),
            "no nominal measurement"
        );
        let m = FragmentErrorModel::new(-2.5, 1.6, 100, 100).with_node_error_da(0.01, 0.04);
        hi.set_fragment_error(Some(m));
        let Tolerance::Da(w) = hi.node_score_window() else {
            panic!("Da expected");
        };
        assert!((w - 0.25).abs() < 1e-12, "{w}");
        hi.set_fragment_error(Some(
            FragmentErrorModel::new(0.0, 1.0, 100, 100).with_node_error_da(0.0, 0.3),
        ));
        assert_eq!(
            hi.node_score_window(),
            Tolerance::Da(0.5),
            "capped at the stored window"
        );
        let mut lo_p = low_res_param();
        lo_p.mme = Tolerance::Da(0.5);
        let lo = RankScorer::new(&lo_p);
        lo.set_fragment_error(Some(m));
        assert_eq!(lo.node_score_window(), Tolerance::Da(0.5));
        let mut ov = RankScorer::new(&p);
        ov.set_fragment_error(Some(m));
        ov.set_fragment_tol_override(Some(Tolerance::Ppm(15.0)));
        assert_eq!(ov.node_score_window(), Tolerance::Da(0.5));
    }

    #[test]
    fn feature_match_tolerance_honors_override() {
        use model::tolerance::Tolerance;
        let mut s = RankScorer::new(&low_res_param());
        s.set_fragment_tol_override(Some(Tolerance::Da(0.6)));
        assert_eq!(s.feature_match_tolerance(), Tolerance::Da(0.6));
        s.set_fragment_tol_override(None);
        assert_eq!(s.feature_match_tolerance(), Tolerance::Da(0.5));
    }

    #[test]
    fn finite_llr_matches_plain_ln_for_valid_inputs() {
        let got = super::finite_llr(0.6, 0.1, 1.0);
        let want = (0.6_f32 / (0.1_f32 * 1.0)).ln();
        assert_eq!(got, want);
    }

    #[test]
    fn finite_llr_stays_finite_on_corrupt_inputs() {
        assert!(
            super::finite_llr(0.0, 0.1, 1.0).is_finite(),
            "zero ion freq"
        );
        assert!(
            super::finite_llr(0.6, 0.0, 1.0).is_finite(),
            "zero noise freq"
        );
        assert!(super::finite_llr(0.6, 0.1, 0.0).is_finite(), "zero norm");
        assert!(
            super::finite_llr(f32::NAN, 0.1, 1.0).is_finite(),
            "NaN ion freq"
        );
        assert!(
            super::finite_llr(0.6, f32::NAN, 1.0).is_finite(),
            "NaN noise freq"
        );
        assert!(
            super::finite_llr(-1.0, 0.1, 1.0).is_finite(),
            "negative ion freq"
        );
        assert!(
            super::finite_llr(0.6, -0.1, 1.0).is_finite(),
            "negative noise freq"
        );
    }

    #[test]
    fn rank_scorer_new_tables_are_finite_even_with_zero_freqs() {
        let mut param = tiny_param();
        let part = Partition {
            charge: 2,
            parent_mass: 1500.0,
            seg_num: 0,
        };
        let zeroed = IonType::Prefix {
            charge: 1,
            offset_bits: 0.0_f32.to_bits(),
            loss_class: 0,
        };
        param
            .rank_dist_table
            .get_mut(&part)
            .unwrap()
            .insert(zeroed, vec![0.0, 0.0, 0.0, 0.0]);
        let scorer = RankScorer::new(&param);
        for table in scorer.log_table.values() {
            assert!(
                table.iter().all(|v| v.is_finite()),
                "all LLR entries must be finite"
            );
        }
    }

    #[test]
    fn node_score_log_formula() {
        let param = tiny_param();
        let scorer = RankScorer::new(&param);
        let part = Partition {
            charge: 2,
            parent_mass: 1500.0,
            seg_num: 0,
        };
        let ion = IonType::Prefix {
            charge: 1,
            offset_bits: 0.0_f32.to_bits(),
            loss_class: 0,
        };

        let s1 = scorer.node_score(part, ion, 1);
        assert!(
            (s1 - 6.0_f32.ln()).abs() < 1e-5,
            "rank1: got {s1}, expected {}",
            6.0_f32.ln()
        );

        let s2 = scorer.node_score(part, ion, 2);
        assert!((s2 - 1.5_f32.ln()).abs() < 1e-5);
    }

    #[test]
    fn rank_above_max_clamps() {
        let param = tiny_param();
        let scorer = RankScorer::new(&param);
        let part = Partition {
            charge: 2,
            parent_mass: 1500.0,
            seg_num: 0,
        };
        let ion = IonType::Prefix {
            charge: 1,
            offset_bits: 0.0_f32.to_bits(),
            loss_class: 0,
        };

        let s5 = scorer.node_score(part, ion, 5);
        let expected = (0.05_f32 / 0.3_f32).ln();
        assert!((s5 - expected).abs() < 1e-5);
    }

    #[test]
    fn missing_ion_score_uses_last_slot() {
        let param = tiny_param();
        let scorer = RankScorer::new(&param);
        let part = Partition {
            charge: 2,
            parent_mass: 1500.0,
            seg_num: 0,
        };
        let ion = IonType::Prefix {
            charge: 1,
            offset_bits: 0.0_f32.to_bits(),
            loss_class: 0,
        };

        let s_missing = scorer.missing_ion_score(part, ion);
        let expected = (0.001_f32 / 0.4_f32).ln();
        assert!((s_missing - expected).abs() < 1e-5);
    }

    #[test]
    fn chargeorseg_uses_min_of_ion_charge_and_num_segments() {
        let mut param = tiny_param();
        let part = Partition {
            charge: 2,
            parent_mass: 1500.0,
            seg_num: 0,
        };
        let ion3 = IonType::Prefix {
            charge: 3,
            offset_bits: 0.0_f32.to_bits(),
            loss_class: 0,
        };
        let ion_freqs = vec![0.6_f32, 0.3, 0.05, 0.001];
        param
            .rank_dist_table
            .get_mut(&part)
            .unwrap()
            .insert(ion3, ion_freqs);

        let scorer = RankScorer::new(&param);
        let s1 = scorer.node_score(part, ion3, 1);
        assert!((s1 - 6.0_f32.ln()).abs() < 1e-5);
    }

    #[test]
    fn unknown_partition_returns_zero() {
        let param = tiny_param();
        let scorer = RankScorer::new(&param);
        let unknown = Partition {
            charge: 99,
            parent_mass: 0.0,
            seg_num: 0,
        };
        let ion = IonType::Prefix {
            charge: 1,
            offset_bits: 0.0_f32.to_bits(),
            loss_class: 0,
        };
        assert_eq!(scorer.node_score(unknown, ion, 1), 0.0);
        assert_eq!(scorer.missing_ion_score(unknown, ion), 0.0);
    }

    #[test]
    fn loss_ions_routed_to_loss_logs_not_intact_path() {
        let mut param = tiny_param();
        let part = Partition {
            charge: 2,
            parent_mass: 1500.0,
            seg_num: 0,
        };
        let intact = IonType::Prefix {
            charge: 1,
            offset_bits: 0.0_f32.to_bits(),
            loss_class: 0,
        };
        let loss = IonType::Prefix {
            charge: 1,
            offset_bits: 0.0_f32.to_bits(),
            loss_class: 1,
        };
        param
            .rank_dist_table
            .get_mut(&part)
            .unwrap()
            .insert(loss, vec![0.6_f32, 0.3, 0.05, 0.001]);
        param
            .partition_ion_types_cache
            .insert(part, vec![intact, loss]);

        let scorer = RankScorer::new(&param);

        assert!(scorer
            .partition_ion_logs(&part)
            .iter()
            .all(|(ion, _)| !ion.is_loss()));
        assert!(scorer
            .partition_ion_logs(&part)
            .iter()
            .any(|(ion, _)| *ion == intact));
        assert!(scorer.has_loss_tables());
        assert!(scorer
            .partition_loss_ion_logs(&part)
            .iter()
            .any(|(ion, _)| *ion == loss));
        assert!(!RankScorer::new(&tiny_param()).has_loss_tables());
    }

    #[test]
    fn ion_existence_score_out_of_domain_prob_peak_is_bounded_zero() {
        let mut param = tiny_param();
        let part = Partition {
            charge: 2,
            parent_mass: 1500.0,
            seg_num: 0,
        };
        param
            .ion_existence_table
            .insert(part, vec![0.7, 0.1, 0.1, 0.1]);
        let scorer = RankScorer::new(&param);

        for &index in &[0usize, 1, 2, 3] {
            let s = scorer.ion_existence_score(part, index, 1.5);
            assert!(
                s.is_finite(),
                "out-of-domain prob_peak must give finite score"
            );
            assert_eq!(s, 0.0, "out-of-domain edge contributes neutral 0");
            let s_nan = scorer.ion_existence_score(part, index, f32::NAN);
            assert_eq!(s_nan, 0.0);
        }

        let s_valid = scorer.ion_existence_score(part, 3, 0.3);
        assert!(s_valid.is_finite());
        assert_ne!(
            s_valid, 0.0,
            "valid prob_peak must score on the normal path"
        );
    }

    #[test]
    fn unknown_ion_returns_zero() {
        let param = tiny_param();
        let scorer = RankScorer::new(&param);
        let part = Partition {
            charge: 2,
            parent_mass: 1500.0,
            seg_num: 0,
        };
        let unknown_ion = IonType::Suffix {
            charge: 1,
            offset_bits: 0.0_f32.to_bits(),
            loss_class: 0,
        };
        assert_eq!(scorer.node_score(part, unknown_ion, 1), 0.0);
    }
}
