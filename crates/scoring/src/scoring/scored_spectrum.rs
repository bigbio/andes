use std::sync::OnceLock;

use crate::param_model::{IonType, Param, Partition, PrecursorOffsetFrequency};
use crate::scoring::rank_scorer::RankScorer;
use model::mass::nominal_from;
use model::peptide::Peptide;
use model::protocol::Protocol;
use model::spectrum::Spectrum;

const PROTON: f64 = model::mass::PROTON;

pub(crate) type SegmentPartitionCache = Vec<(Partition, Vec<(IonType, Vec<f32>)>)>;
pub(crate) type SegmentPartitionSlice<'a> = &'a [(Partition, Vec<(IonType, Vec<f32>)>)];
type DeconvResult = (Option<Vec<(f64, f32)>>, Option<Vec<u32>>);

struct RankKeptCtx {
    prob_peak: f32,
    main_ion: IonType,
    parent_mass: f64,
    charge: u8,
    segment_partition_cache: SegmentPartitionCache,
    prefix_score_cache: Vec<f32>,
    suffix_score_cache: Vec<f32>,
}

#[inline]
fn soft_blend(matched: f32, missing: f32, dm: f32, sigma: f32) -> f32 {
    let sigma = sigma.max(1e-6);
    let z = dm / sigma;
    let w = (-0.5 * z * z).exp();
    w * matched + (1.0 - w) * missing
}

#[derive(Debug, Clone, Copy)]
pub struct IonMatchFact {
    pub partition: Partition,
    pub ion_type: IonType,
    pub rank: Option<u32>,
    pub error_bin: Option<u32>,
}

#[derive(Debug, Clone)]
pub struct ScoredSpectrum<'a> {
    spec: &'a Spectrum,
    ranks: Vec<u32>,
    deconv_peaks: Option<Vec<(f64, f32)>>,
    deconv_ranks: Option<Vec<u32>>,
    kept_count: usize,
    total_intensity: f64,
    max_active_intensity: f64,
    pub(crate) prob_peak: f32,
    pub(crate) main_ion: IonType,
    pub(crate) parent_mass: f64,
    pub(crate) charge: u8,
    segment_partition_cache: SegmentPartitionCache,
    prefix_score_cache: Vec<f32>,
    suffix_score_cache: Vec<f32>,
    observed_mass_cache: std::cell::RefCell<Vec<f64>>,
    gbdt_logit_by_rank: Vec<f32>,
}

enum PeakFilterEnv {
    Disabled,
    Override(f64, usize),
    Unset,
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ScoringSettings {
    pub peak_filter: Option<(f64, usize)>,
}

static SCORING_SETTINGS: OnceLock<ScoringSettings> = OnceLock::new();

pub fn init_scoring_settings(settings: ScoringSettings) {
    let _ = SCORING_SETTINGS.set(settings);
}

#[inline]
pub(crate) fn scoring_settings() -> ScoringSettings {
    SCORING_SETTINGS.get().copied().unwrap_or_default()
}

fn peak_filter_env() -> &'static PeakFilterEnv {
    static CACHE: OnceLock<PeakFilterEnv> = OnceLock::new();
    CACHE.get_or_init(|| {
        let (w, k) = match scoring_settings().peak_filter {
            Some((w, k)) => (Some(w), Some(k)),
            None => (None, None),
        };
        match (w, k) {
            (Some(w), _) if w <= 0.0 => PeakFilterEnv::Disabled,
            (Some(w), Some(k)) => PeakFilterEnv::Override(w, k),
            _ => PeakFilterEnv::Unset,
        }
    })
}

impl<'a> ScoredSpectrum<'a> {
    pub fn new(spec: &'a Spectrum, scorer: &RankScorer, charge: u8) -> Self {
        let param = scorer.param();
        let n = spec.peaks.len();

        let filter_entries: &[PrecursorOffsetFrequency] = param
            .precursor_off_map
            .get(&(charge as i32))
            .map(Vec::as_slice)
            .or_else(|| {
                if param.precursor_off_map.is_empty() {
                    return None;
                }
                param
                    .precursor_off_map
                    .iter()
                    .min_by_key(|(k, _)| ((**k - charge as i32).abs(), **k))
                    .map(|(_, v)| v.as_slice())
            })
            .unwrap_or(&[]);

        let neutral_mass = (spec.precursor_mz - PROTON) * (charge as f64);
        let filter_mzs: Vec<(f64, f64)> = filter_entries
            .iter()
            .filter_map(|pof| {
                let c = (charge as i32 - pof.reduced_charge) as f64;
                if c <= 0.0 {
                    return None;
                }
                let filter_mz = (neutral_mass + c * PROTON) / c + (pof.offset as f64);
                let tol_da = pof.tolerance.as_da(filter_mz);
                Some((filter_mz, tol_da))
            })
            .collect();

        let mut ranks = vec![u32::MAX; n];
        let mut kept: Vec<(usize, f32, f64)> = Vec::with_capacity(n);
        for (i, &(mz, intensity)) in spec.peaks.iter().enumerate() {
            let filtered = filter_mzs.iter().any(|&(fmz, tol)| (mz - fmz).abs() <= tol);
            if !filtered {
                kept.push((i, intensity, mz));
            }
        }

        let window_kk: Option<(f64, usize)> = match peak_filter_env() {
            PeakFilterEnv::Disabled => None,
            PeakFilterEnv::Override(w, k) => Some((*w, *k)),
            PeakFilterEnv::Unset => match param.data_type.protocol {
                Protocol::TMT | Protocol::ITRAQ | Protocol::ITRAQPhospho => Some((100.0, 20)),
                _ => None,
            },
        };
        if let Some((w, k)) = window_kk {
            if w > 0.0 && k > 0 && !kept.is_empty() {
                kept.sort_by(|a, b| {
                    let wa = (a.2 / w).floor() as i64;
                    let wb = (b.2 / w).floor() as i64;
                    wa.cmp(&wb)
                        .then_with(|| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal))
                });
                let mut out: Vec<(usize, f32, f64)> = Vec::with_capacity(kept.len());
                let mut cur_win = i64::MIN;
                let mut in_win = 0usize;
                for &p in kept.iter() {
                    let win = (p.2 / w).floor() as i64;
                    if win != cur_win {
                        cur_win = win;
                        in_win = 0;
                    }
                    if in_win < k {
                        out.push(p);
                        in_win += 1;
                    }
                }
                kept = out;
            }
        }

        let kept_count = kept.len();

        let total_intensity: f64 = kept.iter().map(|&(_, intensity, _)| intensity as f64).sum();

        kept.sort_by(|a, b| {
            b.1.partial_cmp(&a.1)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal))
        });
        for (rank_minus_one, &(orig_idx, _, _)) in kept.iter().enumerate() {
            ranks[orig_idx] = (rank_minus_one + 1) as u32;
        }

        let parent_mass = neutral_mass;

        let (deconv_peaks, deconv_ranks): DeconvResult = if param.apply_deconvolution {
            let tol = param.deconvolution_error_tolerance as f64;
            let (dp, dr) = deconvolute_spectrum(&spec.peaks, &ranks, charge, tol);
            (Some(dp), Some(dr))
        } else {
            (None, None)
        };

        let mme_raw = scorer.node_score_window().raw_value();
        let approx_num_bins = if mme_raw > 0.0 {
            parent_mass / (mme_raw * 2.0)
        } else {
            1.0
        };
        let active_count = match &deconv_peaks {
            Some(dp) => dp.len(),
            None => kept_count,
        };
        let peak_count = if active_count == 0 { 1 } else { active_count } as f64;
        let prob_peak = (peak_count / approx_num_bins.max(1.0)) as f32;

        let last_seg = (param.num_segments - 1).max(0) as usize;
        let part = param.partition_for(charge, parent_mass, last_seg);
        let main_ion = main_ion_from_param(param, part);

        let num_segs = param.num_segments.max(0) as usize;
        let segment_partition_cache: SegmentPartitionCache = (0..num_segs)
            .map(|seg| {
                let p = param.partition_for(charge, parent_mass, seg);
                let logs = scorer.partition_ion_logs(&p).to_vec();
                (p, logs)
            })
            .collect();

        let cache_len = (nominal_from(parent_mass).max(0) as usize) + 1;
        let mut prefix_score_cache = vec![0.0; cache_len];
        let mut suffix_score_cache = vec![0.0; cache_len];
        let (cache_peaks, cache_ranks): (&[(f64, f32)], &[u32]) =
            match (&deconv_peaks, &deconv_ranks) {
                (Some(dp), Some(dr)) => (dp.as_slice(), dr.as_slice()),
                _ => (spec.peaks.as_slice(), ranks.as_slice()),
            };

        let gbdt_logit_by_rank: Vec<f32> = match param.gbdt_peak_model.as_ref() {
            None => Vec::new(),
            Some(model) => {
                use crate::peak_features::{extract_peak_features, PeakFeatureCtx};
                let max_r = cache_ranks
                    .iter()
                    .copied()
                    .filter(|&r| r != u32::MAX)
                    .max()
                    .unwrap_or(0);
                let mut by_rank = vec![0.0_f32; (max_r as usize) + 1];
                let ctx = PeakFeatureCtx::for_spectrum(
                    spec.precursor_mz,
                    charge,
                    parent_mass,
                    cache_peaks,
                    &scorer.node_score_window(),
                );
                let feats = extract_peak_features(cache_peaks, cache_ranks, &ctx);
                let live: Vec<usize> = cache_ranks
                    .iter()
                    .enumerate()
                    .filter(|&(_, &r)| r != u32::MAX && (r as usize) < by_rank.len())
                    .map(|(i, _)| i)
                    .collect();
                if !live.is_empty() {
                    let rows: Vec<&[f32]> = live.iter().map(|&i| feats[i].as_slice()).collect();
                    let mut out = vec![0.0_f32; rows.len()];
                    model.predict_logit_batch(&rows, &mut out);
                    for (&i, &v) in live.iter().zip(out.iter()) {
                        by_rank[cache_ranks[i] as usize] = v;
                    }
                }
                by_rank
            }
        };

        let (fill_peaks, fill_ranks): (Vec<(f64, f32)>, Vec<u32>) = cache_peaks
            .iter()
            .zip(cache_ranks)
            .filter(|&(_, &r)| r != u32::MAX)
            .map(|(&p, &r)| (p, r))
            .unzip();
        for nominal_mass in 1..cache_len {
            let node_nominal = nominal_mass as f64;
            prefix_score_cache[nominal_mass] = Self::directional_node_score_inner(
                &fill_peaks,
                &fill_ranks,
                &segment_partition_cache,
                scorer,
                node_nominal,
                true,
                charge,
                parent_mass,
                &gbdt_logit_by_rank,
            );
            suffix_score_cache[nominal_mass] = Self::directional_node_score_inner(
                &fill_peaks,
                &fill_ranks,
                &segment_partition_cache,
                scorer,
                node_nominal,
                false,
                charge,
                parent_mass,
                &gbdt_logit_by_rank,
            );
        }

        let parent_nominal = nominal_from(parent_mass).max(0) as usize;
        let observed_mass_cache =
            std::cell::RefCell::new(vec![f64::NEG_INFINITY; parent_nominal + 1]);

        let max_active_intensity = match (&deconv_peaks, &deconv_ranks) {
            (Some(peaks), Some(dranks)) => max_active_intensity_of(peaks, dranks),
            _ => max_active_intensity_of(&spec.peaks, &ranks),
        };

        Self {
            spec,
            ranks,
            kept_count,
            total_intensity,
            max_active_intensity,
            prob_peak,
            main_ion,
            parent_mass,
            charge,
            segment_partition_cache,
            prefix_score_cache,
            suffix_score_cache,
            deconv_peaks,
            deconv_ranks,
            observed_mass_cache,
            gbdt_logit_by_rank,
        }
    }

    pub fn new_without_filtering(spec: &'a Spectrum) -> Self {
        let n = spec.peaks.len();
        let kept: Vec<(usize, f32, f64)> = spec
            .peaks
            .iter()
            .enumerate()
            .map(|(i, &(mz, intensity))| (i, intensity, mz))
            .collect();
        let kept_count = kept.len();
        let ranks = vec![u32::MAX; n];
        let prob_peak = 1.0_f32;
        let main_ion = IonType::Prefix {
            charge: 1,
            offset_bits: 0.0_f32.to_bits(),
            loss_class: 0,
        };
        let charge = spec.precursor_charge.map(|z| z.max(1) as u8).unwrap_or(2);
        let parent_mass = (spec.precursor_mz - PROTON) * (charge as f64);
        Self::rank_kept(
            spec,
            kept,
            kept_count,
            ranks,
            RankKeptCtx {
                prob_peak,
                main_ion,
                parent_mass,
                charge,
                segment_partition_cache: Vec::new(),
                prefix_score_cache: Vec::new(),
                suffix_score_cache: Vec::new(),
            },
        )
    }

    fn rank_kept(
        spec: &'a Spectrum,
        mut kept: Vec<(usize, f32, f64)>,
        kept_count: usize,
        mut ranks: Vec<u32>,
        ctx: RankKeptCtx,
    ) -> Self {
        let total_intensity: f64 = kept.iter().map(|&(_, intensity, _)| intensity as f64).sum();
        kept.sort_by(|a, b| {
            b.1.partial_cmp(&a.1)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal))
        });
        for (rank_minus_one, &(orig_idx, _, _)) in kept.iter().enumerate() {
            ranks[orig_idx] = (rank_minus_one + 1) as u32;
        }
        let max_active_intensity = max_active_intensity_of(&spec.peaks, &ranks);
        Self {
            spec,
            ranks,
            kept_count,
            total_intensity,
            max_active_intensity,
            prob_peak: ctx.prob_peak,
            main_ion: ctx.main_ion,
            parent_mass: ctx.parent_mass,
            charge: ctx.charge,
            segment_partition_cache: ctx.segment_partition_cache,
            prefix_score_cache: ctx.prefix_score_cache,
            suffix_score_cache: ctx.suffix_score_cache,
            deconv_peaks: None,
            deconv_ranks: None,
            observed_mass_cache: std::cell::RefCell::new(Vec::new()),
            gbdt_logit_by_rank: Vec::new(),
        }
    }

    pub fn main_ion_direction(&self) -> bool {
        self.main_ion.is_prefix()
    }

    #[inline]
    pub fn active_peaks_and_ranks(&self) -> (&[(f64, f32)], &[u32]) {
        match (&self.deconv_peaks, &self.deconv_ranks) {
            (Some(peaks), Some(ranks)) => (peaks.as_slice(), ranks.as_slice()),
            _ => (self.spec.peaks.as_slice(), self.ranks.as_slice()),
        }
    }

    #[inline]
    pub fn is_deconvoluted(&self) -> bool {
        self.deconv_peaks.is_some()
    }

    pub fn parent_mass(&self) -> f64 {
        self.parent_mass
    }

    pub fn dump_active_peaks(&self) -> Vec<(u32, f64, f32)> {
        let (peaks, ranks) = self.active_peaks_and_ranks();
        let mut out: Vec<(u32, f64, f32)> = peaks
            .iter()
            .zip(ranks.iter())
            .filter(|(_, &rank)| rank != u32::MAX)
            .map(|(&(mz, intensity), &rank)| (rank, mz, intensity))
            .collect();
        out.sort_by(|a, b| a.0.cmp(&b.0));
        out
    }

    pub fn cached_split_score(&self, prefix_nominal: i32, suffix_nominal: i32) -> Option<i32> {
        if prefix_nominal < 0 || suffix_nominal < 0 {
            return None;
        }
        let pref = *self.prefix_score_cache.get(prefix_nominal as usize)?;
        let suff = *self.suffix_score_cache.get(suffix_nominal as usize)?;
        Some((pref + suff).round() as i32)
    }

    pub fn cached_split_score_f32(&self, prefix_nominal: i32, suffix_nominal: i32) -> Option<f32> {
        if prefix_nominal < 0 || suffix_nominal < 0 {
            return None;
        }
        let pref = *self.prefix_score_cache.get(prefix_nominal as usize)?;
        let suff = *self.suffix_score_cache.get(suffix_nominal as usize)?;
        Some(pref + suff)
    }

    pub fn cached_prefix_score(&self, prefix_nominal: i32) -> Option<f32> {
        if prefix_nominal < 0 {
            return None;
        }
        self.prefix_score_cache
            .get(prefix_nominal as usize)
            .copied()
    }

    pub fn cached_suffix_score(&self, suffix_nominal: i32) -> Option<f32> {
        if suffix_nominal < 0 {
            return None;
        }
        self.suffix_score_cache
            .get(suffix_nominal as usize)
            .copied()
    }

    pub fn charge(&self) -> u8 {
        self.charge
    }

    pub fn set_main_ion_for_test(&mut self, ion: IonType) {
        self.main_ion = ion;
    }

    pub fn peak_count(&self) -> usize {
        self.spec.peaks.len()
    }

    pub fn peak_count_after_filtering(&self) -> usize {
        self.kept_count
    }

    pub fn local_peak_density(&self, mz: f64, hw: f64) -> f64 {
        if hw <= 0.0 {
            return 0.0;
        }
        let peaks = self.active_peaks_and_ranks().0;
        let lo = peaks.partition_point(|&(m, _)| m < mz - hw);
        let hi = peaks.partition_point(|&(m, _)| m <= mz + hw);
        (hi - lo) as f64 / (2.0 * hw)
    }

    pub fn total_intensity(&self) -> f64 {
        self.total_intensity
    }

    pub fn max_active_intensity(&self) -> f64 {
        self.max_active_intensity
    }

    pub fn nearest_peak_full(&self, target_mz: f64, tolerance_da: f64) -> Option<(u32, f32, f64)> {
        let (peaks, ranks) = self.active_peaks_and_ranks();
        if peaks.is_empty() {
            return None;
        }
        let lo_mz = target_mz - tolerance_da;
        let hi_mz = target_mz + tolerance_da;
        let start = peaks.partition_point(|&(mz, _)| mz < lo_mz);
        let mut best: Option<(usize, f32)> = None;
        for i in start..peaks.len() {
            let (mz, intensity) = peaks[i];
            if mz > hi_mz {
                break;
            }
            if ranks[i] == u32::MAX {
                continue;
            }
            if best
                .as_ref()
                .is_none_or(|(_, best_int)| intensity > *best_int)
            {
                best = Some((i, intensity));
            }
        }
        best.map(|(i, _)| {
            let (peak_mz, intensity) = peaks[i];
            (ranks[i], intensity, peak_mz)
        })
    }

    pub fn nearest_peak_rank(&self, target_mz: f64, tolerance_da: f64) -> Option<u32> {
        let (peaks, ranks) = self.active_peaks_and_ranks();
        if peaks.is_empty() {
            return None;
        }
        let lo_mz = target_mz - tolerance_da;
        let hi_mz = target_mz + tolerance_da;
        let start = peaks.partition_point(|&(mz, _)| mz < lo_mz);
        let mut best: Option<(usize, f32)> = None;
        for i in start..peaks.len() {
            let (mz, intensity) = peaks[i];
            if mz > hi_mz {
                break;
            }
            if ranks[i] == u32::MAX {
                continue;
            }
            if best
                .as_ref()
                .is_none_or(|(_, best_int)| intensity > *best_int)
            {
                best = Some((i, intensity));
            }
        }
        best.map(|(i, _)| ranks[i])
    }

    #[cfg(test)]
    pub(crate) fn peak_rank_at(&self, idx: usize) -> Option<u32> {
        let r = *self.ranks.get(idx)?;
        if r == u32::MAX {
            None
        } else {
            Some(r)
        }
    }

    pub fn node_score(
        &self,
        prefix_nominal: f64,
        suffix_nominal: f64,
        scorer: &RankScorer,
        charge: u8,
        parent_mass: f64,
        fragment_tolerance_da: f64,
    ) -> i32 {
        let pref = self.directional_node_score(
            prefix_nominal,
            true,
            scorer,
            charge,
            parent_mass,
            fragment_tolerance_da,
        );
        let suff = self.directional_node_score(
            suffix_nominal,
            false,
            scorer,
            charge,
            parent_mass,
            fragment_tolerance_da,
        );
        (pref + suff).round() as i32
    }

    pub fn node_score_f32(
        &self,
        prefix_nominal: f64,
        suffix_nominal: f64,
        scorer: &RankScorer,
        charge: u8,
        parent_mass: f64,
        fragment_tolerance_da: f64,
    ) -> f32 {
        let pref = self.directional_node_score(
            prefix_nominal,
            true,
            scorer,
            charge,
            parent_mass,
            fragment_tolerance_da,
        );
        let suff = self.directional_node_score(
            suffix_nominal,
            false,
            scorer,
            charge,
            parent_mass,
            fragment_tolerance_da,
        );
        pref + suff
    }

    fn directional_node_score(
        &self,
        nominal_mass: f64,
        is_prefix: bool,
        scorer: &RankScorer,
        charge: u8,
        parent_mass: f64,
        _fragment_tolerance_da: f64,
    ) -> f32 {
        let (peaks, ranks) = self.active_peaks_and_ranks();
        Self::directional_node_score_inner(
            peaks,
            ranks,
            &self.segment_partition_cache,
            scorer,
            nominal_mass,
            is_prefix,
            charge,
            parent_mass,
            &self.gbdt_logit_by_rank,
        )
    }

    #[allow(
        clippy::too_many_arguments,
        reason = "private inner driver tightly coupled to the scoring loop; all args are distinct"
    )]
    fn directional_node_score_inner(
        peaks: &[(f64, f32)],
        ranks: &[u32],
        segment_partition_cache: SegmentPartitionSlice<'_>,
        scorer: &RankScorer,
        nominal_mass: f64,
        is_prefix: bool,
        charge: u8,
        parent_mass: f64,
        gbdt_logit_by_rank: &[f32],
    ) -> f32 {
        let max_rank = scorer.max_rank();
        let max_rank_idx = max_rank as usize;
        let mut total = 0.0_f32;
        visit_directional_node_ion_matches(
            peaks,
            ranks,
            segment_partition_cache,
            scorer,
            nominal_mass,
            None,
            is_prefix,
            charge,
            parent_mass,
            |_, _, matched_peak, logs, theo_mz, tol_da| {
                let missing = if max_rank_idx < logs.len() {
                    logs[max_rank_idx]
                } else {
                    0.0
                };
                let score = match matched_peak {
                    Some((rank, peak_mz)) => {
                        let idx = rank.min(max_rank).max(1) as usize - 1;
                        let base = if idx < logs.len() { logs[idx] } else { 0.0 };
                        let gbdt = gbdt_logit_by_rank
                            .get(rank as usize)
                            .copied()
                            .unwrap_or(0.0);
                        let matched = base + gbdt;
                        let dm = (peak_mz - theo_mz) as f32;
                        soft_blend(matched, missing, dm, tol_da as f32)
                    }
                    None => missing,
                };
                total += score;
            },
        );
        total
    }

    #[allow(
        clippy::too_many_arguments,
        reason = "mirrors node_score's split-scoring signature"
    )]
    pub fn loss_node_score(
        &self,
        prefix_nominal: f64,
        suffix_nominal: f64,
        scorer: &RankScorer,
        charge: u8,
        parent_mass: f64,
        prefix_losses: &[(f64, u8)],
        suffix_losses: &[(f64, u8)],
    ) -> i32 {
        if prefix_losses.is_empty() && suffix_losses.is_empty() {
            return 0;
        }
        let pref = self.directional_loss_node_score(
            prefix_nominal,
            true,
            scorer,
            charge,
            parent_mass,
            prefix_losses,
        );
        let suff = self.directional_loss_node_score(
            suffix_nominal,
            false,
            scorer,
            charge,
            parent_mass,
            suffix_losses,
        );
        (pref + suff).round() as i32
    }

    fn directional_loss_node_score(
        &self,
        nominal_mass: f64,
        is_prefix: bool,
        scorer: &RankScorer,
        charge: u8,
        parent_mass: f64,
        active_losses: &[(f64, u8)],
    ) -> f32 {
        if active_losses.is_empty() {
            return 0.0;
        }
        let (peaks, ranks) = self.active_peaks_and_ranks();
        let max_rank = scorer.max_rank();
        let max_rank_idx = max_rank as usize;
        let mut total = 0.0_f32;
        visit_directional_loss_ion_matches(
            peaks,
            ranks,
            &self.segment_partition_cache,
            scorer,
            nominal_mass,
            is_prefix,
            charge,
            parent_mass,
            active_losses,
            |_, _, rank, logs, _, _| {
                let score = match rank {
                    Some(rank) => {
                        let idx = rank.min(max_rank).max(1) as usize - 1;
                        if idx < logs.len() {
                            logs[idx]
                        } else {
                            0.0
                        }
                    }
                    None => {
                        if max_rank_idx < logs.len() {
                            logs[max_rank_idx]
                        } else {
                            0.0
                        }
                    }
                };
                total += score;
            },
        );
        total
    }

    pub fn observed_node_mass(
        &self,
        node_nominal: i32,
        scorer: &RankScorer,
        charge: u8,
        _parent_mass: f64,
    ) -> Option<f64> {
        let _ = charge;
        if node_nominal == 0 {
            return Some(0.0);
        }

        let idx = node_nominal as usize;
        {
            let cache = self.observed_mass_cache.borrow();
            if idx < cache.len() {
                let cached = cache[idx];
                if cached == f64::INFINITY {
                    return None;
                }
                if cached.is_finite() {
                    return Some(cached);
                }
            }
        }

        let theo_mz = self.main_ion.mz(node_nominal as f64);
        let tol_da = scorer.node_score_window().as_da(theo_mz);
        let (peaks, ranks) = self.active_peaks_and_ranks();
        let lo_mz = theo_mz - tol_da;
        let hi_mz = theo_mz + tol_da;
        let start = peaks.partition_point(|&(mz, _)| mz < lo_mz);
        let mut best_peak_mz: Option<(f64, f32)> = None;
        for i in start..peaks.len() {
            let (mz, intensity) = peaks[i];
            if mz > hi_mz {
                break;
            }
            if ranks[i] == u32::MAX {
                continue;
            }
            if best_peak_mz
                .as_ref()
                .is_none_or(|&(_, best_int)| intensity > best_int)
            {
                best_peak_mz = Some((mz, intensity));
            }
        }
        let result = best_peak_mz.map(|(peak_mz, _)| self.main_ion.mass_from_mz(peak_mz));

        {
            let mut cache = self.observed_mass_cache.borrow_mut();
            if idx < cache.len() {
                cache[idx] = match result {
                    Some(m) => m,
                    None => f64::INFINITY,
                };
            }
        }

        result
    }

    pub(crate) fn edge_existence_facts(
        &self,
        cur_nominal: i32,
        prev_nominal: i32,
        scorer: &RankScorer,
        charge: u8,
        parent_mass: f64,
    ) -> (Partition, usize, Option<f64>, Option<f64>) {
        let cur_mass = self.observed_node_mass(cur_nominal, scorer, charge, parent_mass);
        let prev_mass = self.observed_node_mass(prev_nominal, scorer, charge, parent_mass);

        let mut idx = 0usize;
        if cur_mass.is_some() {
            idx += 1;
        }
        if prev_mass.is_some() {
            idx += 2;
        }

        let last_seg = (scorer.param().num_segments - 1).max(0) as usize;
        let part = match self.segment_partition_cache.get(last_seg) {
            Some((p, _)) => *p,
            None => scorer.param().partition_for(charge, parent_mass, last_seg),
        };

        (part, idx, cur_mass, prev_mass)
    }

    pub fn edge_score(
        &self,
        cur_nominal: i32,
        prev_nominal: i32,
        theo_aa_mass: f64,
        scorer: &RankScorer,
        charge: u8,
        parent_mass: f64,
    ) -> i32 {
        if scorer.param().error_scaling_factor == 0 {
            return 0;
        }
        if scorer.param().ion_existence_table.is_empty() {
            return 0;
        }

        let (part, idx, cur_mass, prev_mass) =
            self.edge_existence_facts(cur_nominal, prev_nominal, scorer, charge, parent_mass);

        let mut s = scorer.ion_existence_score(part, idx, self.prob_peak);

        if idx == 3 {
            let delta = cur_mass.unwrap() - prev_mass.unwrap() - theo_aa_mass;
            s += scorer.error_score(part, delta as f32);
        }

        s.round() as i32
    }

    pub fn ion_match_facts(&self, peptide: &Peptide, scorer: &RankScorer) -> Vec<IonMatchFact> {
        let param = scorer.param();
        let max_rank = scorer.max_rank();
        let esf = param.error_scaling_factor;

        let (peaks, ranks) = self.active_peaks_and_ranks();

        let n = peptide.length();
        if n < 2 {
            return Vec::new();
        }

        let peptide_nominal = peptide.nominal_residue_mass();
        let mut prefix_nominal_arr: Vec<i32> = Vec::with_capacity(n + 1);
        prefix_nominal_arr.push(0);
        let mut prefix_acc = 0.0_f64;
        for s in 1..n {
            let aa = &peptide.residues[s - 1];
            prefix_acc += aa.mass + aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta);
            prefix_nominal_arr.push(nominal_from(prefix_acc));
        }
        let last_aa = &peptide.residues[n - 1];
        prefix_acc += last_aa.mass + last_aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta);
        prefix_nominal_arr.push(nominal_from(prefix_acc));

        let mut out: Vec<IonMatchFact> = Vec::new();

        let emit_losses = scorer
            .param()
            .data_type
            .activation
            .predicts_neutral_losses()
            && peptide.residues.iter().any(|aa| {
                aa.mod_
                    .as_ref()
                    .is_some_and(|m| m.loss_class != 0 && !m.neutral_losses.is_empty())
            });

        #[allow(clippy::needless_range_loop)]
        for split in 1..n {
            let prefix_nom = prefix_nominal_arr[split] as f64;
            let suffix_nom = (peptide_nominal - prefix_nominal_arr[split]) as f64;

            for &(is_prefix, nominal_mass) in &[(true, prefix_nom), (false, suffix_nom)] {
                visit_directional_node_ion_matches(
                    peaks,
                    ranks,
                    &self.segment_partition_cache,
                    scorer,
                    nominal_mass,
                    None,
                    is_prefix,
                    self.charge,
                    self.parent_mass,
                    |partition, ion, matched_peak, _logs, theo_mz, _tol_da| {
                        let (rank, error_bin) = ion_fact_rank_and_error_from_match(
                            matched_peak,
                            theo_mz,
                            max_rank,
                            esf,
                        );
                        out.push(IonMatchFact {
                            partition,
                            ion_type: ion,
                            rank,
                            error_bin,
                        });
                    },
                );
            }

            if emit_losses {
                let prefix_losses = span_losses(peptide, 0..split);
                let suffix_losses = span_losses(peptide, split..n);
                for &(is_prefix, nominal_mass, losses) in &[
                    (true, prefix_nom, &prefix_losses),
                    (false, suffix_nom, &suffix_losses),
                ] {
                    visit_directional_loss_ion_facts(
                        peaks,
                        ranks,
                        &self.segment_partition_cache,
                        scorer,
                        nominal_mass,
                        None,
                        is_prefix,
                        self.charge,
                        self.parent_mass,
                        losses,
                        |partition, loss_ion, rank, theo_mz, tol_da| {
                            let (rank, error_bin) = ion_fact_rank_and_error(
                                rank, peaks, ranks, theo_mz, tol_da, max_rank, esf,
                            );
                            out.push(IonMatchFact {
                                partition,
                                ion_type: loss_ion,
                                rank,
                                error_bin,
                            });
                        },
                    );
                }
            }
        }
        out
    }

    pub fn noise_match_facts(
        &self,
        peptide: &Peptide,
        scorer: &RankScorer,
    ) -> Vec<(Partition, Option<u32>, Option<u32>)> {
        if peptide.length() < 2 {
            return Vec::new();
        }
        let mut rev = peptide.residues.clone();
        rev.reverse();
        let decoy = Peptide::new(rev, peptide.pre, peptide.post);
        self.ion_match_facts(&decoy, scorer)
            .into_iter()
            .map(|f| (f.partition, f.rank, f.error_bin))
            .collect()
    }

    pub fn dense_noise_facts(
        &self,
        peptide: &Peptide,
        scorer: &RankScorer,
        n_samples: usize,
    ) -> Vec<(Partition, Option<u32>, Option<u32>)> {
        let max_rank = scorer.max_rank();
        let esf = scorer.param().error_scaling_factor;
        if n_samples == 0 {
            return Vec::new();
        }
        let (peaks, ranks) = self.active_peaks_and_ranks();
        let peptide_nominal = peptide.nominal_residue_mass();
        let lo: i64 = 57;
        let hi: i64 = (peptide_nominal as i64 - 57).max(lo + 1);
        let mut out: Vec<(Partition, Option<u32>, Option<u32>)> = Vec::with_capacity(n_samples * 2);
        for i in 0..n_samples {
            let nominal = (lo + (hi - lo) * i as i64 / n_samples as i64) as f64;
            for &is_prefix in &[true, false] {
                visit_directional_node_ion_matches(
                    peaks,
                    ranks,
                    &self.segment_partition_cache,
                    scorer,
                    nominal,
                    None,
                    is_prefix,
                    self.charge,
                    self.parent_mass,
                    |partition, _ion, matched_peak, _logs, theo_mz, _tol_da| {
                        let (r, ebin) = ion_fact_rank_and_error_from_match(
                            matched_peak,
                            theo_mz,
                            max_rank,
                            esf,
                        );
                        out.push((partition, r, ebin));
                    },
                );
            }
        }
        out
    }
}

#[allow(
    clippy::too_many_arguments,
    reason = "mirrors the scoring loop's argument bundle"
)]
#[allow(
    clippy::too_many_arguments,
    reason = "private inner matcher for the node-scoring hot path; every argument is an orthogonal, already-resolved input and bundling them would add a struct build per call"
)]
fn visit_directional_node_ion_matches<F>(
    peaks: &[(f64, f32)],
    ranks: &[u32],
    segment_partition_cache: SegmentPartitionSlice<'_>,
    scorer: &RankScorer,
    nominal_mass: f64,
    exact_mass: Option<f64>,
    is_prefix: bool,
    charge: u8,
    parent_mass: f64,
    mut visit: F,
) where
    F: FnMut(Partition, IonType, Option<(u32, f64)>, &[f32], f64, f64),
{
    use crate::param_model::IonType;

    let param = scorer.param();
    let node_window = scorer.node_score_window();
    let mme = &node_window;
    let num_segs = param.num_segments as usize;
    let use_cache = !segment_partition_cache.is_empty();
    #[allow(clippy::needless_range_loop)]
    for seg in 0..num_segs {
        let (partition, ion_logs_slice): (Partition, &[(IonType, Vec<f32>)]) = if use_cache {
            (
                segment_partition_cache[seg].0,
                segment_partition_cache[seg].1.as_slice(),
            )
        } else {
            let p = param.partition_for(charge, parent_mass, seg);
            let logs = scorer.partition_ion_logs(&p);
            (p, logs)
        };
        if !ion_logs_slice.is_empty() {
            for (ion, logs) in ion_logs_slice {
                let theo_mz = match (is_prefix, *ion) {
                    (true, IonType::Prefix { .. }) => match exact_mass {
                        Some(m) => ion.mz_exact(m),
                        None => ion.mz(nominal_mass),
                    },
                    (false, IonType::Suffix { .. }) => match exact_mass {
                        Some(m) => ion.mz_exact(m),
                        None => ion.mz(nominal_mass),
                    },
                    _ => continue,
                };
                if param.segment_num(theo_mz, parent_mass) != seg {
                    continue;
                }
                let tol_da = mme.as_da(theo_mz);
                let matched = nearest_peak_rank_and_mz_in(peaks, ranks, theo_mz, tol_da);
                visit(partition, *ion, matched, logs, theo_mz, tol_da);
            }
        } else {
            for &ion in param.ion_types_for_partition_slice(charge, parent_mass, seg) {
                let theo_mz = match (is_prefix, ion) {
                    (true, IonType::Prefix { .. }) => match exact_mass {
                        Some(m) => ion.mz_exact(m),
                        None => ion.mz(nominal_mass),
                    },
                    (false, IonType::Suffix { .. }) => match exact_mass {
                        Some(m) => ion.mz_exact(m),
                        None => ion.mz(nominal_mass),
                    },
                    _ => continue,
                };
                if param.segment_num(theo_mz, parent_mass) != seg {
                    continue;
                }
                let tol_da = mme.as_da(theo_mz);
                let matched = nearest_peak_rank_and_mz_in(peaks, ranks, theo_mz, tol_da);
                visit(partition, ion, matched, &[], theo_mz, tol_da);
            }
        }
    }
}

#[allow(
    clippy::too_many_arguments,
    reason = "mirrors visit_directional_node_ion_matches plus the active-loss slice"
)]
fn visit_directional_loss_ion_matches<F>(
    peaks: &[(f64, f32)],
    ranks: &[u32],
    segment_partition_cache: SegmentPartitionSlice<'_>,
    scorer: &RankScorer,
    nominal_mass: f64,
    is_prefix: bool,
    charge: u8,
    parent_mass: f64,
    active_losses: &[(f64, u8)],
    mut visit: F,
) where
    F: FnMut(Partition, IonType, Option<u32>, &[f32], f64, f64),
{
    use crate::param_model::IonType;
    if active_losses.is_empty() {
        return;
    }
    let param = scorer.param();
    let node_window = scorer.node_score_window();
    let mme = &node_window;
    let num_segs = param.num_segments as usize;
    let use_cache = !segment_partition_cache.is_empty();
    #[allow(clippy::needless_range_loop)]
    for seg in 0..num_segs {
        let partition = if use_cache {
            segment_partition_cache[seg].0
        } else {
            param.partition_for(charge, parent_mass, seg)
        };
        let loss_logs = scorer.partition_loss_ion_logs(&partition);
        if loss_logs.is_empty() {
            continue;
        }
        for (ion, logs) in loss_logs {
            match (is_prefix, ion) {
                (true, IonType::Prefix { .. }) | (false, IonType::Suffix { .. }) => {}
                _ => continue,
            }
            let cls = ion.loss_class();
            let ion_charge = ion.charge().unwrap_or(1).max(1) as f64;
            let base_mz = ion.mz(nominal_mass);
            for &(loss, lcls) in active_losses {
                if lcls != cls {
                    continue;
                }
                let theo_mz = base_mz - loss / ion_charge;
                if theo_mz <= 0.0 || param.segment_num(theo_mz, parent_mass) != seg {
                    continue;
                }
                let tol_da = mme.as_da(theo_mz);
                let rank = nearest_peak_rank_in(peaks, ranks, theo_mz, tol_da);
                visit(partition, *ion, rank, logs, theo_mz, tol_da);
            }
        }
    }
}

pub(crate) fn span_losses(peptide: &Peptide, range: std::ops::Range<usize>) -> Vec<(f64, u8)> {
    let mut v: Vec<(f64, u8)> = Vec::new();
    for i in range {
        if let Some(m) = peptide.residues[i].mod_.as_ref() {
            if m.loss_class != 0 {
                for &l in &m.neutral_losses {
                    if !v
                        .iter()
                        .any(|&(el, ec)| ec == m.loss_class && (el - l).abs() < 1e-9)
                    {
                        v.push((l, m.loss_class));
                    }
                }
            }
        }
    }
    v
}

fn ion_fact_rank_and_error(
    rank: Option<u32>,
    peaks: &[(f64, f32)],
    ranks: &[u32],
    theo_mz: f64,
    tol_da: f64,
    max_rank: u32,
    esf: i32,
) -> (Option<u32>, Option<u32>) {
    let matched = rank.map(|r| {
        let mz = matched_peak_mz(peaks, ranks, theo_mz, tol_da).unwrap_or(theo_mz);
        (r, mz)
    });
    ion_fact_rank_and_error_from_match(matched, theo_mz, max_rank, esf)
}

fn ion_fact_rank_and_error_from_match(
    matched: Option<(u32, f64)>,
    theo_mz: f64,
    max_rank: u32,
    esf: i32,
) -> (Option<u32>, Option<u32>) {
    match matched {
        Some((r, peak_mz)) => {
            let clamped = r.min(max_rank).max(1);
            let ebin = if esf > 0 {
                let error_da = (peak_mz - theo_mz) as f32;
                let mut idx = (error_da * esf as f32).round() as i32;
                if idx > esf {
                    idx = esf;
                } else if idx < -esf {
                    idx = -esf;
                }
                idx += esf;
                Some(idx as u32)
            } else {
                None
            };
            (Some(clamped), ebin)
        }
        None => (None, None),
    }
}

#[allow(
    clippy::too_many_arguments,
    reason = "mirrors visit_directional_node_ion_matches plus the active-loss slice"
)]
fn visit_directional_loss_ion_facts<F>(
    peaks: &[(f64, f32)],
    ranks: &[u32],
    segment_partition_cache: SegmentPartitionSlice<'_>,
    scorer: &RankScorer,
    nominal_mass: f64,
    exact_mass: Option<f64>,
    is_prefix: bool,
    charge: u8,
    parent_mass: f64,
    active_losses: &[(f64, u8)],
    mut visit: F,
) where
    F: FnMut(Partition, IonType, Option<u32>, f64, f64),
{
    use crate::param_model::IonType;
    if active_losses.is_empty() {
        return;
    }
    let param = scorer.param();
    let node_window = scorer.node_score_window();
    let mme = &node_window;
    let num_segs = param.num_segments as usize;
    let use_cache = !segment_partition_cache.is_empty();
    #[allow(clippy::needless_range_loop)]
    for seg in 0..num_segs {
        let (partition, intact_logs): (Partition, &[(IonType, Vec<f32>)]) = if use_cache {
            (
                segment_partition_cache[seg].0,
                segment_partition_cache[seg].1.as_slice(),
            )
        } else {
            let p = param.partition_for(charge, parent_mass, seg);
            (p, scorer.partition_ion_logs(&p))
        };
        for (ion, _logs) in intact_logs {
            let (icharge, ioff) = match (is_prefix, ion) {
                (
                    true,
                    IonType::Prefix {
                        charge,
                        offset_bits,
                        ..
                    },
                ) => (*charge, *offset_bits),
                (
                    false,
                    IonType::Suffix {
                        charge,
                        offset_bits,
                        ..
                    },
                ) => (*charge, *offset_bits),
                _ => continue,
            };
            let ion_charge = icharge.max(1) as f64;
            let base_mz = match exact_mass {
                Some(m) => ion.mz_exact(m),
                None => ion.mz(nominal_mass),
            };
            for &(loss, cls) in active_losses {
                let theo_mz = base_mz - loss / ion_charge;
                if theo_mz <= 0.0 || param.segment_num(theo_mz, parent_mass) != seg {
                    continue;
                }
                let tol_da = mme.as_da(theo_mz);
                let rank = nearest_peak_rank_in(peaks, ranks, theo_mz, tol_da);
                let loss_ion = if is_prefix {
                    IonType::Prefix {
                        charge: icharge,
                        offset_bits: ioff,
                        loss_class: cls,
                    }
                } else {
                    IonType::Suffix {
                        charge: icharge,
                        offset_bits: ioff,
                        loss_class: cls,
                    }
                };
                visit(partition, loss_ion, rank, theo_mz, tol_da);
            }
        }
    }
}

fn matched_peak_mz(
    peaks: &[(f64, f32)],
    ranks: &[u32],
    target_mz: f64,
    tolerance_da: f64,
) -> Option<f64> {
    if peaks.is_empty() {
        return None;
    }
    let lo_mz = target_mz - tolerance_da;
    let hi_mz = target_mz + tolerance_da;
    let start = peaks.partition_point(|&(mz, _)| mz < lo_mz);
    let mut best: Option<(f64, f32)> = None;
    for i in start..peaks.len() {
        let (mz, intensity) = peaks[i];
        if mz > hi_mz {
            break;
        }
        if ranks[i] == u32::MAX {
            continue;
        }
        if best
            .as_ref()
            .is_none_or(|(_, best_int)| intensity > *best_int)
        {
            best = Some((mz, intensity));
        }
    }
    best.map(|(mz, _)| mz)
}

fn nearest_peak_rank_in(
    peaks: &[(f64, f32)],
    ranks: &[u32],
    target_mz: f64,
    tolerance_da: f64,
) -> Option<u32> {
    nearest_peak_rank_and_mz_in(peaks, ranks, target_mz, tolerance_da).map(|(rank, _mz)| rank)
}

fn nearest_peak_rank_and_mz_in(
    peaks: &[(f64, f32)],
    ranks: &[u32],
    target_mz: f64,
    tolerance_da: f64,
) -> Option<(u32, f64)> {
    if peaks.is_empty() {
        return None;
    }
    let lo_mz = target_mz - tolerance_da;
    let hi_mz = target_mz + tolerance_da;
    let start = peaks.partition_point(|&(mz, _)| mz < lo_mz);
    let mut best: Option<(usize, f32)> = None;
    for i in start..peaks.len() {
        let (mz, intensity) = peaks[i];
        if mz > hi_mz {
            break;
        }
        if ranks[i] == u32::MAX {
            continue;
        }
        if best
            .as_ref()
            .is_none_or(|(_, best_int)| intensity > *best_int)
        {
            best = Some((i, intensity));
        }
    }
    best.map(|(i, _)| (ranks[i], peaks[i].0))
}

fn max_active_intensity_of(peaks: &[(f64, f32)], ranks: &[u32]) -> f64 {
    peaks
        .iter()
        .zip(ranks.iter())
        .filter(|(_, &rank)| rank != u32::MAX)
        .map(|(&(_, intensity), _)| intensity as f64)
        .fold(0.0_f64, f64::max)
}

#[inline]
pub fn scored_for_charge<'s, 'a>(
    scored_per_charge: &'s [(u8, ScoredSpectrum<'a>)],
    z: u8,
) -> Option<&'s ScoredSpectrum<'a>> {
    scored_per_charge
        .iter()
        .find(|(charge, _)| *charge == z)
        .map(|(_, scored)| scored)
}

fn deconvolute_spectrum(
    peaks: &[(f64, f32)],
    ranks: &[u32],
    precursor_charge: u8,
    tol: f64,
) -> (Vec<(f64, f32)>, Vec<u32>) {
    const ISOTOPE: f64 = model::mass::ISOTOPE;

    let n = peaks.len();
    if n == 0 {
        return (Vec::new(), Vec::new());
    }
    let mut ignore = vec![false; n];
    let mut out: Vec<(f64, f32, u32)> = Vec::with_capacity(n);
    let charge_i32 = precursor_charge as i32;

    for i in 0..n {
        if ignore[i] {
            continue;
        }
        let (mut p_mz, p_int) = peaks[i];
        let p_rank = ranks[i];

        for ion_charge_i in 2..charge_i32.min(4) {
            let ion_charge = ion_charge_i as f64;
            let expected_diff = ISOTOPE / ion_charge;
            let mut is_deconvoluted = false;
            for j in (i + 1)..n {
                let (p2_mz, p2_int) = peaks[j];
                let diff = p2_mz - p_mz - expected_diff;
                if diff > -tol && diff < tol {
                    ignore[j] = true;
                    let p_new_mz = ion_charge * p_mz - (ion_charge - 1.0) * PROTON;
                    let p2_new_mz = ion_charge * p2_mz - (ion_charge - 1.0) * PROTON;
                    p_mz = p_new_mz;
                    is_deconvoluted = true;

                    let neutral_mass = ion_charge * (p2_mz - PROTON);
                    let p3_diff_expected =
                        model::isotope::averagine_m1_m2_spacing(neutral_mass) / ion_charge;
                    for k in (j + 1)..n {
                        let (p3_mz, p3_int) = peaks[k];
                        let diff2 = p3_mz - p2_mz - p3_diff_expected;
                        if diff2 > -tol && diff2 < tol {
                            ignore[k] = true;
                            let p3_new_mz = ion_charge * p3_mz - (ion_charge - 1.0) * PROTON;
                            out.push((p3_new_mz, p3_int, ranks[k]));
                            break;
                        } else if diff2 > tol {
                            break;
                        }
                    }
                    out.push((p2_new_mz, p2_int, ranks[j]));
                    break;
                } else if diff > tol {
                    break;
                }
            }
            if is_deconvoluted {
                break;
            }
        }
        out.push((p_mz, p_int, p_rank));
    }

    out.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let mut out_peaks: Vec<(f64, f32)> = Vec::with_capacity(out.len());
    let mut out_ranks: Vec<u32> = Vec::with_capacity(out.len());
    for (mz, intensity, rank) in out {
        out_peaks.push((mz, intensity));
        out_ranks.push(rank);
    }
    (out_peaks, out_ranks)
}

fn main_ion_from_param(param: &Param, partition: crate::param_model::Partition) -> IonType {
    let fallback = IonType::Prefix {
        charge: 1,
        offset_bits: 0.0_f32.to_bits(),
        loss_class: 0,
    };
    let num_segments = param.num_segments.max(1) as usize;
    let mut ion_freq: std::collections::HashMap<IonType, f32> = std::collections::HashMap::new();
    for seg in 0..num_segments {
        let part = crate::param_model::Partition {
            charge: partition.charge,
            parent_mass: partition.parent_mass,
            seg_num: seg as i32,
        };
        if let Some(frag_list) = param.frag_off_table.get(&part) {
            for f in frag_list {
                if matches!(f.ion_type, IonType::Noise) {
                    continue;
                }
                *ion_freq.entry(f.ion_type).or_insert(0.0) += f.frequency;
            }
        }
    }
    let tiebreak = |ion: IonType| -> (u8, i32, u32, u8) {
        match ion {
            IonType::Suffix {
                charge,
                offset_bits,
                loss_class,
            } => (0, charge, offset_bits, loss_class),
            IonType::Prefix {
                charge,
                offset_bits,
                loss_class,
            } => (1, charge, offset_bits, loss_class),
            IonType::Noise => (2, 0, 0, 0),
        }
    };
    let mut best_ion: Option<IonType> = None;
    let mut best_freq = f32::NEG_INFINITY;
    let mut best_key = (u8::MAX, i32::MAX, u32::MAX, u8::MAX);
    for (&ion, &freq) in &ion_freq {
        let key = tiebreak(ion);
        if freq > best_freq || (freq == best_freq && key < best_key) {
            best_freq = freq;
            best_key = key;
            best_ion = Some(ion);
        }
    }
    best_ion.unwrap_or(fallback)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn main_ion_deterministic_with_placeholder_zero_frequencies() {
        use crate::param_model::{FragmentOffsetFrequency, IonType, Partition};
        let part = Partition {
            charge: 2,
            parent_mass: 1500.0,
            seg_num: 0,
        };
        let zero = |ion| FragmentOffsetFrequency {
            ion_type: ion,
            frequency: 0.0,
        };
        let b1 = IonType::Prefix {
            charge: 1,
            offset_bits: 1.0_f32.to_bits(),
            loss_class: 0,
        };
        let y1 = IonType::Suffix {
            charge: 1,
            offset_bits: 19.0_f32.to_bits(),
            loss_class: 0,
        };
        let b2 = IonType::Prefix {
            charge: 2,
            offset_bits: 1.0_f32.to_bits(),
            loss_class: 0,
        };
        let y2 = IonType::Suffix {
            charge: 2,
            offset_bits: 19.0_f32.to_bits(),
            loss_class: 0,
        };

        let mut p = crate::testutil::tiny_param();
        p.num_segments = 1;
        p.frag_off_table.clear();
        p.frag_off_table.insert(
            part,
            vec![zero(b1), zero(y1), zero(b2), zero(y2), zero(IonType::Noise)],
        );
        p.rebuild_cache();

        let a = main_ion_from_param(&p, part);
        let b = main_ion_from_param(&p, part);
        assert_eq!(a, b, "main ion must be deterministic across calls");
        assert_eq!(
            a, y1,
            "tied-zero frequencies must resolve to y1 (Suffix, charge 1), not a hash-order ion"
        );
    }

    /// The scorer blends with σ equal to the full fragment tolerance (`soft_blend(.., tol_da)`),
    /// so a peak at the window edge keeps exp(-1/2), about 61 %, of the matched credit.
    const TOL_DA: f32 = 0.5;

    #[test]
    fn soft_match_centre_gives_full_credit() {
        let matched = 4.0_f32;
        let missing = -2.0_f32;
        let s = soft_blend(matched, missing, 0.0, TOL_DA);
        assert!(
            (s - matched).abs() < 1e-6,
            "centre peak must keep full matched credit, got {s}"
        );
    }

    #[test]
    fn soft_match_edge_keeps_most_of_the_matched_credit() {
        let matched = 4.0_f32;
        let missing = -2.0_f32;
        let s = soft_blend(matched, missing, TOL_DA, TOL_DA);
        assert!(
            s > 0.5 * (matched + missing),
            "with σ = tolerance an edge peak leans toward matched, got {s}"
        );
        let w = (-0.5f32).exp();
        let expected = w * matched + (1.0 - w) * missing;
        assert!(
            (s - expected).abs() < 1e-5,
            "edge blend math mismatch: {s} vs {expected}"
        );
    }

    #[test]
    fn soft_match_monotonic_in_abs_dm() {
        let matched = 5.0_f32;
        let missing = -1.0_f32;
        let mut prev = f32::INFINITY;
        for k in 0..=25 {
            let dm = k as f32 * 0.02;
            let s = soft_blend(matched, missing, dm, TOL_DA);
            assert!(
                s <= prev + 1e-6,
                "blend must be non-increasing in |Δm|: dm={dm} s={s} prev={prev}"
            );
            prev = s;
        }
        let plus = soft_blend(matched, missing, 0.1, TOL_DA);
        let minus = soft_blend(matched, missing, -0.1, TOL_DA);
        assert!(
            (plus - minus).abs() < 1e-6,
            "blend must be symmetric in Δm sign"
        );
    }

    #[test]
    fn soft_match_large_factor_approaches_hard() {
        let matched = 3.0_f32;
        let missing = -2.0_f32;
        let huge_sigma = 1000.0_f32;
        let s = soft_blend(matched, missing, 0.4, huge_sigma);
        assert!(
            (s - matched).abs() < 1e-3,
            "σ ≫ tol must degenerate to hard matching, got {s}"
        );
    }

    #[test]
    fn nearest_peak_rank_and_mz_single_scan_matches_split_helpers() {
        let peaks = vec![(99.95, 10.0_f32), (100.05, 50.0_f32), (100.5, 5.0_f32)];
        let ranks = vec![2_u32, 1, 3];
        let tol = 0.1;
        let combined = nearest_peak_rank_and_mz_in(&peaks, &ranks, 100.0, tol);
        let rank = nearest_peak_rank_in(&peaks, &ranks, 100.0, tol);
        let mz = matched_peak_mz(&peaks, &ranks, 100.0, tol);
        assert_eq!(
            combined.map(|(r, _)| r),
            rank,
            "single-scan rank must match"
        );
        assert_eq!(combined.map(|(_, m)| m), mz, "single-scan m/z must match");
        assert_eq!(combined, Some((1, 100.05)));
        assert_eq!(
            nearest_peak_rank_and_mz_in(&peaks, &ranks, 200.0, tol),
            None
        );
    }

    use crate::param_model::{IonType, Partition, SpecDataType};
    use crate::scoring::rank_scorer::RankScorer;
    use crate::testutil::tiny_param_with_ions;

    fn spec(peaks: &[(f64, f32)]) -> Spectrum {
        Spectrum {
            title: "test".into(),
            precursor_mz: 500.0,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks: peaks.to_vec(),
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
            mobility_window: None,
        }
    }

    #[test]
    fn isobaric_windowed_peak_filter_caps_dense_window() {
        use model::protocol::Protocol;
        let peaks: Vec<(f64, f32)> = (0..30)
            .map(|i| (200.0 + i as f64 * 0.9, (100 - i) as f32))
            .collect();
        let s = spec(&peaks);

        let mut p = tiny_param_with_ions();
        p.data_type.protocol = Protocol::Automatic;
        let scorer = RankScorer::new(&p);
        assert_eq!(
            ScoredSpectrum::new(&s, &scorer, 2)
                .dump_active_peaks()
                .len(),
            30,
            "non-isobaric must keep all peaks"
        );

        let mut pt = tiny_param_with_ions();
        pt.data_type.protocol = Protocol::TMT;
        let scorer_t = RankScorer::new(&pt);
        let kept = ScoredSpectrum::new(&s, &scorer_t, 2).dump_active_peaks();
        assert_eq!(kept.len(), 20, "TMT must cap to top-20 per 100-Da window");
        let min_kept = kept
            .iter()
            .map(|&(_, _, it)| it)
            .fold(f32::INFINITY, f32::min);
        assert!(
            min_kept >= 81.0,
            "kept peaks must be the most intense (got min {min_kept})"
        );
    }

    #[test]
    fn loss_node_score_scores_loss_shifted_peak_via_trained_table() {
        let part = Partition {
            charge: 2,
            parent_mass: 1000.0,
            seg_num: 0,
        };
        let loss_ion = IonType::Prefix {
            charge: 1,
            offset_bits: 0.0_f32.to_bits(),
            loss_class: 1,
        };

        let mut p = tiny_param_with_ions();
        p.rank_dist_table
            .get_mut(&part)
            .unwrap()
            .insert(loss_ion, vec![0.6_f32, 0.3, 0.05, 0.001]);

        let nominal = 500.0_f64;
        let loss_mass = 162.0528_f64;
        let theo_mz = loss_ion.mz(nominal) - loss_mass;

        let s = spec(&[(theo_mz, 100.0)]);
        let scorer = RankScorer::new(&p);
        let ss = ScoredSpectrum::new(&s, &scorer, 2);

        let active = [(loss_mass, 1u8)];
        let with_loss = ss.loss_node_score(nominal, 0.0, &scorer, 2, 1000.0, &active, &[]);
        assert!(
            with_loss > 0,
            "loss ion at rank 1 should score positive, got {with_loss}"
        );

        assert_eq!(
            ss.loss_node_score(nominal, 0.0, &scorer, 2, 1000.0, &[], &[]),
            0
        );

        assert_eq!(
            ss.loss_node_score(nominal, 0.0, &scorer, 2, 1000.0, &[(loss_mass, 2u8)], &[]),
            0
        );

        let plain = RankScorer::new(&tiny_param_with_ions());
        let ss2 = ScoredSpectrum::new(&s, &plain, 2);
        assert_eq!(
            ss2.loss_node_score(nominal, 0.0, &plain, 2, 1000.0, &active, &[]),
            0
        );
    }

    #[test]
    fn ion_match_facts_emits_loss_facts_for_loss_peptide() {
        use model::amino_acid::AminoAcid;
        use model::modification::{ModLocation, Modification, ResidueSpec};
        use model::peptide::Peptide;

        let p = tiny_param_with_ions();
        let scorer = RankScorer::new(&p);

        let m = Modification {
            name: "L".into(),
            mass_delta: 0.0,
            residue: ResidueSpec::Specific(b'T'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: vec![162.0528],
            loss_class: 1,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        let t = AminoAcid::standard(b'T').unwrap().with_mod(m);
        let w1 = AminoAcid::standard(b'W').unwrap();
        let w2 = AminoAcid::standard(b'W').unwrap();
        let pep = Peptide::new(vec![t, w1, w2], b'_', b'-');

        let intact = IonType::Prefix {
            charge: 1,
            offset_bits: 0.0_f32.to_bits(),
            loss_class: 0,
        };
        let prefix_nominal = nominal_from(pep.residues[0].mass + pep.residues[1].mass);
        let theo = intact.mz(prefix_nominal as f64) - 162.0528;

        let s = Spectrum {
            title: "loss".into(),
            precursor_mz: 600.0,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks: vec![(theo, 100.0)],
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
            mobility_window: None,
        };
        let ss = ScoredSpectrum::new(&s, &scorer, 2);

        let facts = ss.ion_match_facts(&pep, &scorer);
        assert!(
            facts
                .iter()
                .any(|f| f.ion_type.loss_class() == 1 && f.rank.is_some()),
            "expected a matched glyco loss fact (rank Some, class 1)"
        );

        let plain = Peptide::new(
            vec![
                AminoAcid::standard(b'W').unwrap(),
                AminoAcid::standard(b'W').unwrap(),
                AminoAcid::standard(b'A').unwrap(),
            ],
            b'_',
            b'-',
        );
        assert!(
            ss.ion_match_facts(&plain, &scorer)
                .iter()
                .all(|f| !f.ion_type.is_loss()),
            "plain peptide must emit no loss facts"
        );
    }

    #[test]
    fn prob_peak_uses_raw_mme_value_not_da_converted() {
        use crate::param_model::SpecDataType;
        use model::activation::ActivationMethod;
        use model::instrument::InstrumentType;
        use model::protocol::Protocol;
        use model::tolerance::Tolerance;
        use rustc_hash::FxHashMap;

        let precursor_mz = 501.007_276_49_f64;
        let s = Spectrum {
            title: "prob_peak_test".into(),
            precursor_mz,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks: vec![(100.0, 1.0), (200.0, 2.0), (300.0, 3.0)],
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
            mobility_window: None,
        };

        let param = Param {
            version: 10001,
            data_type: SpecDataType {
                activation: ActivationMethod::HCD,
                instrument: InstrumentType::QExactive,
                enzyme: None,
                protocol: Protocol::Automatic,
            },
            mme: Tolerance::Ppm(20.0),
            apply_deconvolution: false,
            deconvolution_error_tolerance: 0.0,
            charge_hist: vec![(2, 100)],
            min_charge: 2,
            max_charge: 2,
            num_segments: 1,
            partitions: vec![],
            num_precursor_off: 0,
            precursor_off_map: FxHashMap::default(),
            frag_off_table: FxHashMap::default(),
            max_rank: 3,
            rank_dist_table: FxHashMap::default(),
            error_scaling_factor: 0,
            ion_err_dist_table: FxHashMap::default(),
            noise_err_dist_table: FxHashMap::default(),
            ion_existence_table: FxHashMap::default(),
            partition_ion_types_cache: FxHashMap::default(),
            gbdt_peak_model: None,
            frag_intensity_model: None,
            rich_ion_model: None,
        };

        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new(&s, &scorer, 2);

        let parent_mass = (precursor_mz - PROTON) * 2.0;
        let raw_mme = 20.0_f64;
        let approx_num_bins = parent_mass / (raw_mme * 2.0);
        let expected_prob_peak = (3.0_f64 / approx_num_bins.max(1.0)) as f32;

        let wrong_approx_num_bins = parent_mass / (parent_mass * 20e-6 * 2.0);
        let wrong_prob_peak = (3.0_f64 / wrong_approx_num_bins.max(1.0)) as f32;

        assert!(
            (expected_prob_peak - wrong_prob_peak).abs() > 0.001,
            "test precondition failed: Ppm raw vs Da-converted did not produce different prob_peak values"
        );

        assert!(
            (ss.prob_peak - expected_prob_peak).abs() < 1e-5,
            "prob_peak={} but expected={} (raw-mme formula). Wrong Da-converted value would be {}",
            ss.prob_peak,
            expected_prob_peak,
            wrong_prob_peak
        );
    }

    fn deconv_param(apply: bool) -> Param {
        use model::activation::ActivationMethod;
        use model::instrument::InstrumentType;
        use model::protocol::Protocol;
        use model::tolerance::Tolerance;
        use rustc_hash::FxHashMap;
        Param {
            version: 10001,
            data_type: SpecDataType {
                activation: ActivationMethod::HCD,
                instrument: InstrumentType::QExactive,
                enzyme: None,
                protocol: Protocol::Automatic,
            },
            mme: Tolerance::Ppm(20.0),
            apply_deconvolution: apply,
            deconvolution_error_tolerance: 0.05,
            charge_hist: vec![(2, 100)],
            min_charge: 2,
            max_charge: 4,
            num_segments: 1,
            partitions: vec![],
            num_precursor_off: 0,
            precursor_off_map: FxHashMap::default(),
            frag_off_table: FxHashMap::default(),
            max_rank: 3,
            rank_dist_table: FxHashMap::default(),
            error_scaling_factor: 0,
            ion_err_dist_table: FxHashMap::default(),
            noise_err_dist_table: FxHashMap::default(),
            ion_existence_table: FxHashMap::default(),
            partition_ion_types_cache: FxHashMap::default(),
            gbdt_peak_model: None,
            frag_intensity_model: None,
            rich_ion_model: None,
        }
    }

    #[test]
    fn deconv_active_for_charge_2_produces_input_equivalent_peaks() {
        let s = Spectrum {
            title: "deconv_test".into(),
            precursor_mz: 501.007_276_49_f64,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks: vec![(100.0, 1.0), (200.0, 2.0), (300.0, 3.0)],
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
            mobility_window: None,
        };
        let param = deconv_param(true);
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new(&s, &scorer, 2);

        let parent_mass = (s.precursor_mz - PROTON) * 2.0;
        let approx = parent_mass / (20.0_f64 * 2.0);
        let expected = (3.0_f64 / approx.max(1.0)) as f32;
        assert!(
            (ss.prob_peak - expected).abs() < 1e-5,
            "charge=2 deconv-active spectrum: prob_peak={} expected={} (active_count=3)",
            ss.prob_peak,
            expected
        );
    }

    #[test]
    fn deconv_active_for_charge_3_uses_post_deconv_peak_count_for_prob_peak() {
        const ISOTOPE: f64 = 1.003_354_83;
        let p1 = 100.0;
        let p2 = p1 + ISOTOPE / 2.0;
        let p3 = 500.0;
        let s = Spectrum {
            title: "deconv_charge3".into(),
            precursor_mz: 401.0,
            precursor_intensity: None,
            precursor_charge: Some(3),
            rt_seconds: None,
            scan: None,
            peaks: vec![(p1, 10.0), (p2, 5.0), (p3, 1.0)],
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
            mobility_window: None,
        };
        let param = deconv_param(true);
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new(&s, &scorer, 3);

        let active_count = ss.deconv_peaks.as_ref().map(|p| p.len()).unwrap_or(0);
        assert!(
            active_count >= 1,
            "deconv_peaks should be populated for charge=3 + apply_deconvolution=true"
        );
        let parent_mass = (401.0 - PROTON) * 3.0;
        let approx = parent_mass / (20.0_f64 * 2.0);
        let expected = (active_count as f64 / approx.max(1.0)) as f32;
        assert!(
            (ss.prob_peak - expected).abs() < 1e-5,
            "charge=3 deconv-active spectrum: prob_peak={} expected={} (post-deconv count={})",
            ss.prob_peak,
            expected,
            active_count
        );
    }

    #[test]
    fn deconv_off_uses_kept_count_for_prob_peak() {
        let s = Spectrum {
            title: "no_deconv".into(),
            precursor_mz: 501.007_276_49_f64,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks: vec![(100.0, 1.0), (200.0, 2.0), (300.0, 3.0), (400.0, 4.0)],
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
            mobility_window: None,
        };
        let param = deconv_param(false);
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new(&s, &scorer, 2);

        let parent_mass = (s.precursor_mz - PROTON) * 2.0;
        let approx = parent_mass / (20.0_f64 * 2.0);
        let expected = (4.0_f64 / approx.max(1.0)) as f32;
        assert!(
            (ss.prob_peak - expected).abs() < 1e-5,
            "deconv-off: prob_peak={} expected={} (kept_count=4)",
            ss.prob_peak,
            expected
        );
        assert!(
            ss.deconv_peaks.is_none(),
            "deconv_peaks must be None when apply_deconvolution=false"
        );
    }

    #[test]
    fn observed_node_mass_picks_highest_intensity_peak_in_window() {
        use model::mass::INTEGER_MASS_SCALER;
        let node_nominal = 100_i32;
        let theo_mz = node_nominal as f64 / INTEGER_MASS_SCALER as f64;
        let closer_mz = theo_mz + 0.09;
        let farther_mz = theo_mz + 0.39;
        let s = spec(&[(closer_mz, 1.0), (farther_mz, 100.0)]);
        let param = tiny_param_with_ions();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let result = ss.observed_node_mass(node_nominal, &scorer, 2, 1000.0);
        let result_mass = result.expect("should find a peak in the window");
        let expected_mass = farther_mz;
        let wrong_mass = closer_mz;
        assert!(
            (result_mass - expected_mass).abs() < 1e-6,
            "expected highest-intensity (farther) peak mass {expected_mass:.6}, \
             got {result_mass:.6} (closest/wrong would be {wrong_mass:.6})"
        );
    }

    #[test]
    fn node_score_does_not_panic_on_empty_spectrum() {
        let s = spec(&[]);
        let param = tiny_param_with_ions();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let n = ss.node_score(100.0, 900.0, &scorer, 2, 1000.0, 0.5);
        assert!(
            n <= 0,
            "missing-ion score on empty spectrum should be non-positive, got {n}"
        );
    }

    #[test]
    fn node_score_nonzero_when_peak_matches_prefix_ion() {
        use model::mass::INTEGER_MASS_SCALER;
        let nominal = 100.0_f64;
        let b1_mz = nominal / INTEGER_MASS_SCALER as f64;
        let s = spec(&[(50.0, 1.0), (b1_mz, 100.0), (200.0, 2.0)]);
        let param = tiny_param_with_ions();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let n = ss.node_score(nominal, 900.0, &scorer, 2, 1000.0, 0.5);
        assert!(
            n > 0,
            "expected positive node_score when b-ion peak present, got {n}"
        );
    }

    #[test]
    fn node_score_prefix_only_match() {
        use model::mass::INTEGER_MASS_SCALER;
        let nominal = 57.0_f64;
        let mz = nominal / INTEGER_MASS_SCALER as f64;
        let s = spec(&[(mz, 50.0), (300.0, 1.0)]);
        let param = tiny_param_with_ions();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let n = ss.node_score(nominal, 900.0, &scorer, 2, 1000.0, 0.5);
        assert!(n > 0, "prefix-only match: expected positive score, got {n}");
    }

    #[test]
    fn node_score_no_matching_ions_returns_negative_or_zero() {
        let s = spec(&[(5000.0, 100.0)]);
        let param = tiny_param_with_ions();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let n = ss.node_score(100.0, 900.0, &scorer, 2, 1000.0, 0.5);
        assert!(
            n <= 0,
            "missing ion should produce non-positive score, got {n}"
        );
    }

    #[test]
    fn node_score_nominal_mass_zero_prefix_returns_zero() {
        let s = spec(&[]);
        let param = tiny_param_with_ions();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let n = ss.node_score(0.0, 1000.0, &scorer, 2, 1000.0, 0.5);
        assert!(
            n <= 0,
            "source-node score with empty spectrum should be non-positive, got {n}"
        );
    }

    #[test]
    fn edge_score_returns_zero_when_table_empty() {
        let s = spec(&[(100.0, 1.0)]);
        let mut param = tiny_param_with_ions();
        param.ion_existence_table.clear();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let e = ss.edge_score(150, 100, 50.0, &scorer, 2, 1000.0);
        assert_eq!(e, 0);
    }

    #[test]
    fn edge_score_returns_zero_when_error_scaling_factor_zero() {
        let s = spec(&[(100.0, 1.0)]);
        let param = tiny_param_with_ions();
        assert_eq!(param.error_scaling_factor, 0);
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let e = ss.edge_score(150, 100, 50.0, &scorer, 2, 1000.0);
        assert_eq!(e, 0);
    }

    #[test]
    fn edge_score_nonzero_with_existence_table() {
        use crate::param_model::{FragmentOffsetFrequency, SpecDataType};
        use model::activation::ActivationMethod;
        use model::instrument::InstrumentType;
        use model::protocol::Protocol;
        use model::tolerance::Tolerance;
        use rustc_hash::FxHashMap;

        let part = Partition {
            charge: 2,
            parent_mass: 1000.0,
            seg_num: 0,
        };
        let prefix1 = IonType::Prefix {
            charge: 1,
            offset_bits: 0.0_f32.to_bits(),
            loss_class: 0,
        };
        let noise = IonType::Noise;

        let ion_freqs = vec![0.6_f32, 0.3, 0.05, 0.001];
        let noise_freqs = vec![0.1_f32, 0.2, 0.3, 0.4];

        let mut ion_table: FxHashMap<IonType, Vec<f32>> = FxHashMap::default();
        ion_table.insert(prefix1, ion_freqs);
        ion_table.insert(noise, noise_freqs);

        let mut rank_dist_table: FxHashMap<Partition, FxHashMap<IonType, Vec<f32>>> =
            FxHashMap::default();
        rank_dist_table.insert(part, ion_table);

        let mut frag_off_table = FxHashMap::default();
        frag_off_table.insert(
            part,
            vec![FragmentOffsetFrequency {
                ion_type: prefix1,
                frequency: 0.7,
            }],
        );

        let error_scaling_factor = 2_i32;
        let dist_len = (error_scaling_factor as usize) * 2 + 1;

        let mut ion_err_dist_table: FxHashMap<Partition, Vec<f32>> = FxHashMap::default();
        ion_err_dist_table.insert(part, vec![0.1_f32, 0.2, 0.4, 0.2, 0.1]);

        let mut noise_err_dist_table: FxHashMap<Partition, Vec<f32>> = FxHashMap::default();
        noise_err_dist_table.insert(part, vec![0.05_f32, 0.1, 0.7, 0.1, 0.05]);

        let mut ion_existence_table: FxHashMap<Partition, Vec<f32>> = FxHashMap::default();
        ion_existence_table.insert(part, vec![0.1_f32, 0.3, 0.3, 0.5]);

        let _ = dist_len;

        let mut param = Param {
            version: 10001,
            data_type: SpecDataType {
                activation: ActivationMethod::HCD,
                instrument: InstrumentType::QExactive,
                enzyme: None,
                protocol: Protocol::Automatic,
            },
            mme: Tolerance::Da(0.5),
            apply_deconvolution: false,
            deconvolution_error_tolerance: 0.0,
            charge_hist: vec![(2, 100)],
            min_charge: 2,
            max_charge: 2,
            num_segments: 1,
            partitions: vec![part],
            num_precursor_off: 0,
            precursor_off_map: FxHashMap::default(),
            frag_off_table,
            max_rank: 3,
            rank_dist_table,
            error_scaling_factor,
            ion_err_dist_table,
            noise_err_dist_table,
            ion_existence_table,
            partition_ion_types_cache: FxHashMap::default(),
            gbdt_peak_model: None,
            frag_intensity_model: None,
            rich_ion_model: None,
        };
        param.rebuild_cache();

        let s = spec(&[]);
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let e = ss.edge_score(150, 100, 50.0, &scorer, 2, 1000.0);
        assert_ne!(
            e, 0,
            "edge_score should be nonzero with populated existence table"
        );
    }

    #[test]
    fn empty_spectrum_yields_no_ranks() {
        let s = spec(&[]);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        assert_eq!(ss.peak_count(), 0);
        assert!(ss.nearest_peak_rank(500.0, 0.1).is_none());
    }

    #[test]
    fn highest_intensity_gets_rank_1() {
        let s = spec(&[(100.0, 1.0), (200.0, 5.0), (300.0, 3.0)]);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        assert_eq!(ss.peak_count(), 3);
        assert_eq!(ss.nearest_peak_rank(200.0, 0.1), Some(1));
        assert_eq!(ss.nearest_peak_rank(300.0, 0.1), Some(2));
        assert_eq!(ss.nearest_peak_rank(100.0, 0.1), Some(3));
    }

    #[test]
    fn nearest_peak_within_tolerance() {
        let s = spec(&[(100.0, 1.0), (200.5, 5.0), (300.0, 3.0)]);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        assert_eq!(ss.nearest_peak_rank(200.4, 0.2), Some(1));
        assert_eq!(ss.nearest_peak_rank(200.5, 0.001), Some(1));
        assert_eq!(ss.nearest_peak_rank(200.4, 0.05), None);
    }

    #[test]
    fn ties_broken_deterministically() {
        let s = spec(&[(100.0, 5.0), (200.0, 5.0)]);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let r1 = ss.nearest_peak_rank(100.0, 0.1).unwrap();
        let r2 = ss.nearest_peak_rank(200.0, 0.1).unwrap();
        assert_ne!(r1, r2);
        assert!(r1 == 1 || r2 == 1);
        assert!(r1 == 2 || r2 == 2);
    }

    #[test]
    fn nearest_peak_full_picks_max_intensity_within_tolerance() {
        let s = spec(&[(99.6, 5.0), (100.0, 9.0), (100.5, 1.0)]);
        let ss = ScoredSpectrum::new_without_filtering(&s);

        let (rank, intensity, mz) = ss.nearest_peak_full(100.45, 0.6).unwrap();
        assert_eq!(
            rank, 1,
            "max-intensity peak (100.0) must win, not nearest-mz (100.5)"
        );
        assert_eq!(mz, 100.0);
        assert_eq!(intensity, 9.0);

        assert_eq!(ss.nearest_peak_rank(100.45, 0.6), Some(1));
    }

    #[test]
    fn nearest_peak_rank_matches_linear_scan_on_many_peaks() {
        let mut peaks: Vec<(f64, f32)> = (0..100)
            .map(|i| (i as f64 * 10.0 + 0.5, (100 - i) as f32))
            .collect();
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let s = Spectrum {
            title: "many".into(),
            precursor_mz: 500.0,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks,
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
            mobility_window: None,
        };
        let ss = ScoredSpectrum::new_without_filtering(&s);

        for target in [50.5, 100.0, 250.0, 333.7, 500.5, 750.5, 999.5] {
            let tol = 5.0_f64;
            let bs_result = ss.nearest_peak_rank(target, tol);
            let bf_result = {
                let mut best: Option<(usize, f64)> = None;
                for (i, &(mz, _)) in s.peaks.iter().enumerate() {
                    if (mz - target).abs() <= tol
                        && best.as_ref().is_none_or(|(_, d)| (mz - target).abs() < *d)
                    {
                        best = Some((i, (mz - target).abs()));
                    }
                }
                best.map(|(i, _)| ss.peak_rank_at(i).unwrap_or(u32::MAX))
            };
            assert_eq!(
                bs_result, bf_result,
                "binary search and linear scan differ at target {target}"
            );
        }
    }

    #[test]
    fn gbdt_term_raises_node_score_for_signal_peak() {
        use crate::gbdt_eval::{GbdtPeakModel, Tree};
        use model::mass::INTEGER_MASS_SCALER;

        let gbdt = GbdtPeakModel {
            packed: Default::default(),
            n_features: 1,
            feature_names: crate::gbdt_eval::indexed_feature_names(1),
            apply_sigmoid: true,
            trees: vec![Tree {
                feature: vec![-1],
                threshold: vec![0.0],
                left: vec![-1],
                right: vec![-1],
                value: vec![4.0],
                default_left: vec![1],
            }],
            iso_x: vec![],
            iso_y: vec![],
        };

        const PROTON_F64: f64 = 1.007_276_49;
        let nominal = 100.0_f64;
        let b1_mz = nominal / INTEGER_MASS_SCALER as f64;
        let precursor_mz = (1000.0 + 2.0 * PROTON_F64) / 2.0;

        let peaks_vec = vec![(b1_mz, 1000.0_f32)];
        let spec_plain = Spectrum {
            title: "gbdt_test".into(),
            precursor_mz,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks: peaks_vec.clone(),
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
            mobility_window: None,
        };

        let param_plain = tiny_param_with_ions();
        let scorer_plain = RankScorer::new(&param_plain);
        let s_plain = ScoredSpectrum::new(&spec_plain, &scorer_plain, 2);
        let parent_mass_plain = s_plain.parent_mass;
        let base =
            s_plain.directional_node_score(nominal, true, &scorer_plain, 2, parent_mass_plain, 0.5);

        let mut param_gbdt = tiny_param_with_ions();
        param_gbdt.gbdt_peak_model = Some(gbdt);
        let scorer_gbdt = RankScorer::new(&param_gbdt);
        let spec_gbdt = Spectrum {
            title: "gbdt_test".into(),
            precursor_mz,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks: peaks_vec,
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
            mobility_window: None,
        };
        let s_gbdt = ScoredSpectrum::new(&spec_gbdt, &scorer_gbdt, 2);
        let parent_mass_gbdt = s_gbdt.parent_mass;
        let with =
            s_gbdt.directional_node_score(nominal, true, &scorer_gbdt, 2, parent_mass_gbdt, 0.5);

        assert!(
            with > base + 3.0,
            "expected GBDT bump ≥ 3.0: base={base}, with={with}"
        );
    }
}

#[cfg(test)]
mod precursor_filter_tests {
    use super::*;
    use crate::param_model::{Param, PrecursorOffsetFrequency, SpecDataType};
    use model::activation::ActivationMethod;
    use model::instrument::InstrumentType;
    use model::protocol::Protocol;
    use model::tolerance::Tolerance;
    use rustc_hash::FxHashMap;

    fn param_with_precursor_filter() -> Param {
        let mut precursor_off_map: FxHashMap<i32, Vec<PrecursorOffsetFrequency>> =
            FxHashMap::default();
        precursor_off_map.insert(
            2,
            vec![PrecursorOffsetFrequency {
                reduced_charge: 2,
                offset: 0.0,
                tolerance: Tolerance::Da(0.5),
                frequency: 1.0,
            }],
        );

        Param {
            version: 10001,
            data_type: SpecDataType {
                activation: ActivationMethod::HCD,
                instrument: InstrumentType::QExactive,
                enzyme: None,
                protocol: Protocol::Automatic,
            },
            mme: Tolerance::Ppm(20.0),
            apply_deconvolution: false,
            deconvolution_error_tolerance: 0.0,
            charge_hist: vec![(2, 100)],
            min_charge: 2,
            max_charge: 2,
            num_segments: 1,
            partitions: vec![],
            num_precursor_off: 1,
            precursor_off_map,
            frag_off_table: FxHashMap::default(),
            max_rank: 3,
            rank_dist_table: FxHashMap::default(),
            error_scaling_factor: 0,
            ion_err_dist_table: FxHashMap::default(),
            noise_err_dist_table: FxHashMap::default(),
            ion_existence_table: FxHashMap::default(),
            partition_ion_types_cache: FxHashMap::default(),
            gbdt_peak_model: None,
            frag_intensity_model: None,
            rich_ion_model: None,
        }
    }

    fn make_spec(precursor_mz: f64, peaks: &[(f64, f32)]) -> Spectrum {
        Spectrum {
            title: "test".into(),
            precursor_mz,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks: peaks.to_vec(),
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
            mobility_window: None,
        }
    }

    fn param_with_precursor_filter_rc0() -> Param {
        let mut precursor_off_map: FxHashMap<i32, Vec<PrecursorOffsetFrequency>> =
            FxHashMap::default();
        precursor_off_map.insert(
            2,
            vec![PrecursorOffsetFrequency {
                reduced_charge: 0,
                offset: 0.0,
                tolerance: Tolerance::Da(0.5),
                frequency: 1.0,
            }],
        );

        Param {
            version: 10001,
            data_type: SpecDataType {
                activation: ActivationMethod::HCD,
                instrument: InstrumentType::QExactive,
                enzyme: None,
                protocol: Protocol::Automatic,
            },
            mme: Tolerance::Ppm(20.0),
            apply_deconvolution: false,
            deconvolution_error_tolerance: 0.0,
            charge_hist: vec![(2, 100)],
            min_charge: 2,
            max_charge: 2,
            num_segments: 1,
            partitions: vec![],
            num_precursor_off: 1,
            precursor_off_map,
            frag_off_table: FxHashMap::default(),
            max_rank: 3,
            rank_dist_table: FxHashMap::default(),
            error_scaling_factor: 0,
            ion_err_dist_table: FxHashMap::default(),
            noise_err_dist_table: FxHashMap::default(),
            ion_existence_table: FxHashMap::default(),
            partition_ion_types_cache: FxHashMap::default(),
            gbdt_peak_model: None,
            frag_intensity_model: None,
            rich_ion_model: None,
        }
    }

    #[test]
    fn precursor_peak_is_filtered_out() {
        let s = make_spec(500.0, &[(100.0, 1.0), (300.0, 5.0), (500.0, 100.0)]);
        let param = param_with_precursor_filter_rc0();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new(&s, &scorer, 2);

        assert!(
            ss.nearest_peak_rank(500.0, 0.1).is_none(),
            "precursor peak at 500.0 should be filtered, but a peak at that m/z was found"
        );

        assert_eq!(ss.nearest_peak_rank(300.0, 0.1), Some(1));
        assert_eq!(ss.nearest_peak_rank(100.0, 0.1), Some(2));
    }

    #[test]
    fn non_precursor_peaks_kept() {
        let s = make_spec(500.0, &[(100.0, 1.0), (200.0, 50.0), (300.0, 5.0)]);
        let param = param_with_precursor_filter_rc0();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new(&s, &scorer, 2);

        assert_eq!(ss.peak_count_after_filtering(), 3);
        assert_eq!(ss.nearest_peak_rank(200.0, 0.1), Some(1));
    }

    #[test]
    fn missing_precursor_off_map_falls_back_to_unfiltered() {
        let mut param = param_with_precursor_filter_rc0();
        param.precursor_off_map.clear();
        let s = make_spec(500.0, &[(100.0, 1.0), (500.0, 100.0)]);
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new(&s, &scorer, 2);
        assert_eq!(ss.peak_count_after_filtering(), 2);
    }

    #[test]
    fn invalid_reduced_charge_skipped() {
        let param = param_with_precursor_filter();
        let s = make_spec(500.0, &[(100.0, 1.0), (500.0, 100.0)]);
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new(&s, &scorer, 2);
        assert_eq!(ss.peak_count_after_filtering(), 2);
    }
}
