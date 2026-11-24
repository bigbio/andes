//! A spectrum prepared for scoring: ranked peaks, node and edge scores per peptide, and the
//! ion-match facts the rescoring features read.

mod deconvolve;
mod matches;
mod settings;

use model::mass::nominal_from;
use model::peptide::Peptide;
use model::protocol::Protocol;
use model::spectrum::Spectrum;

use crate::params::{IonType, Partition, PrecursorOffsetFrequency};
use crate::psm::rank::RankScorer;
use deconvolve::{deconvolute_spectrum, main_ion_from_param, DeconvResult};
pub(crate) use matches::span_losses;
use matches::{
    ion_fact_rank_and_error, ion_fact_rank_and_error_from_match, max_active_intensity_of,
    visit_directional_loss_ion_facts, visit_directional_loss_ion_matches,
    visit_directional_node_ion_matches,
};
pub use settings::{init_scoring_settings, ScoringSettings};
use settings::{peak_filter_env, PeakFilterEnv};

const PROTON: f64 = model::mass::PROTON;

pub(crate) type SegmentPartitionCache = Vec<(Partition, Vec<(IonType, Vec<f32>)>)>;
pub(crate) type SegmentPartitionSlice<'a> = &'a [(Partition, Vec<(IonType, Vec<f32>)>)];

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
    prefix_score_cache: LazyNodeScores,
    suffix_score_cache: LazyNodeScores,
    node_peaks: Vec<(f64, f32)>,
    node_ranks: Vec<u32>,
    node_scorer: Option<&'a RankScorer>,
    observed_mass_cache: std::cell::RefCell<Vec<f64>>,
    gbdt_logit_by_rank: Vec<f32>,
}

/// Prefix or suffix node scores by nominal mass, each computed on its first lookup.
#[derive(Debug, Clone, Default)]
struct LazyNodeScores {
    values: Vec<std::cell::Cell<f32>>,
    known: Vec<std::cell::Cell<bool>>,
}

impl LazyNodeScores {
    /// `len` scores, none computed yet except nominal mass 0 (score 0).
    fn pending(len: usize) -> Self {
        let mut known = vec![std::cell::Cell::new(false); len];
        if let Some(k) = known.first_mut() {
            *k = std::cell::Cell::new(true);
        }
        Self {
            values: vec![std::cell::Cell::new(0.0); len],
            known,
        }
    }

    /// Scores already computed.
    fn known(values: Vec<f32>) -> Self {
        Self {
            known: vec![std::cell::Cell::new(true); values.len()],
            values: values.into_iter().map(std::cell::Cell::new).collect(),
        }
    }
}

impl<'a> ScoredSpectrum<'a> {
    pub fn new(spec: &'a Spectrum, scorer: &'a RankScorer, charge: u8) -> Self {
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
        let (cache_peaks, cache_ranks): (&[(f64, f32)], &[u32]) =
            match (&deconv_peaks, &deconv_ranks) {
                (Some(dp), Some(dr)) => (dp.as_slice(), dr.as_slice()),
                _ => (spec.peaks.as_slice(), ranks.as_slice()),
            };

        let gbdt_logit_by_rank: Vec<f32> = match param.gbdt_peak_model.as_ref() {
            None => Vec::new(),
            Some(model) => {
                use crate::features::peaks::{extract_peak_features, PeakFeatureCtx};
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
            prefix_score_cache: LazyNodeScores::pending(cache_len),
            suffix_score_cache: LazyNodeScores::pending(cache_len),
            node_peaks: fill_peaks,
            node_ranks: fill_ranks,
            node_scorer: Some(scorer),
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
            prefix_score_cache: LazyNodeScores::known(ctx.prefix_score_cache),
            suffix_score_cache: LazyNodeScores::known(ctx.suffix_score_cache),
            node_peaks: Vec::new(),
            node_ranks: Vec::new(),
            node_scorer: None,
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
        let pref = self.cached_node_score(prefix_nominal as usize, true)?;
        let suff = self.cached_node_score(suffix_nominal as usize, false)?;
        Some((pref + suff).round() as i32)
    }

    pub fn cached_split_score_f32(&self, prefix_nominal: i32, suffix_nominal: i32) -> Option<f32> {
        if prefix_nominal < 0 || suffix_nominal < 0 {
            return None;
        }
        let pref = self.cached_node_score(prefix_nominal as usize, true)?;
        let suff = self.cached_node_score(suffix_nominal as usize, false)?;
        Some(pref + suff)
    }

    pub fn cached_prefix_score(&self, prefix_nominal: i32) -> Option<f32> {
        if prefix_nominal < 0 {
            return None;
        }
        self.cached_node_score(prefix_nominal as usize, true)
    }

    pub fn cached_suffix_score(&self, suffix_nominal: i32) -> Option<f32> {
        if suffix_nominal < 0 {
            return None;
        }
        self.cached_node_score(suffix_nominal as usize, false)
    }

    pub fn charge(&self) -> u8 {
        self.charge
    }

    /// The prefix (or suffix) node score at nominal mass `nominal`, computed on first use;
    /// `None` beyond the parent mass.
    fn cached_node_score(&self, nominal: usize, is_prefix: bool) -> Option<f32> {
        let cache = if is_prefix {
            &self.prefix_score_cache
        } else {
            &self.suffix_score_cache
        };
        let value = cache.values.get(nominal)?;
        if !cache.known[nominal].get() {
            let scorer = self
                .node_scorer
                .expect("a node score left to compute has its scorer");
            value.set(Self::directional_node_score_inner(
                &self.node_peaks,
                &self.node_ranks,
                &self.segment_partition_cache,
                scorer,
                nominal as f64,
                is_prefix,
                self.charge,
                self.parent_mass,
                &self.gbdt_logit_by_rank,
            ));
            cache.known[nominal].set(true);
        }
        Some(value.get())
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

#[cfg(test)]
#[path = "spectrum/tests/spectrum.rs"]
mod tests;

#[cfg(test)]
#[path = "spectrum/tests/precursor.rs"]
mod precursor_filter_tests;
