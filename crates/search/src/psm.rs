use std::cmp::Reverse;
use std::collections::BinaryHeap;

#[derive(Debug, Clone, Default)]
pub struct PsmFeatures {
    pub num_matched_main_ions: u32,
    pub longest_b: u32,
    pub longest_y: u32,
    pub longest_y_pct: f32,
    pub matched_ion_ratio: f32,

    pub explained_ion_current_ratio: f32,
    pub n_term_ion_current_ratio: f32,
    pub c_term_ion_current_ratio: f32,
    pub ms2_ion_current: f32,
    pub isolation_window_efficiency: f32,

    pub mean_error_top7: f32,
    pub stdev_error_top7: f32,
    pub mean_rel_error_top7: f32,
    pub stdev_rel_error_top7: f32,

    pub edge_score: i32,

    pub precursor_isotope_kl: f32,
    pub precursor_snr: f32,

    pub delta_raw_score: f32,

    pub tailor_score: f32,

    pub rank_score_float: f32,

    pub ppm_gaussian_score: f32,
    pub longest_complementary_ladder: u32,
    pub complementary_ion_balance: f32,
    pub neutral_loss_ion_count: u32,
    pub mean_matched_intensity_rank: f32,
    pub doubly_charged_matched_ion_count: u32,

    pub chance_match_surprise: f32,
    pub unique_match_fraction: f32,

    pub intensity_signal: f32,

    pub frag_pred_explained: f32,
    pub frag_pred_chance_llr: f32,
    pub frag_topk_observed: f32,

    pub rich_ion_llr: f32,

    pub is_refinement: u32,
    pub num_mods: u32,
    pub refine_mod_class: u32,

    pub mod_site_shifted_matched: f32,
    pub mod_site_shifted_frac: f32,
    pub mod_site_intens_frac: f32,
    pub mod_site_localized: f32,
    pub mod_site_det_count: f32,

    pub mass_competition_evidence: f32,
    pub match_posterior_mean: f32,
    pub match_log_odds: f32,
    pub candidate_rank_entropy: f32,
    pub listwise_score_gap: f32,

    pub strong_score: f32,

    pub strong_score_cal: f32,

    pub delta_rt: f32,
    pub abs_delta_rt: f32,
    pub delta_rt_norm: f32,
    pub ref_competitor: u32,
    pub ref_best_raw: f32,
    pub is_variant: u32,
    pub is_isoform: u32,
    pub is_chain: u32,
    pub delta_ref_score: f32,
    pub isotope_mimic: u32,
    pub predicted_rt_min: f32,
    pub isobaric_rt_margin: f32,
    pub open_site: u16,
    pub open_loc_score: f32,
    pub open_loc_margin: f32,
    pub open_delta: f32,
    pub open_used_shift: bool,
}

pub const TAILOR_MIN_CANDIDATES: u32 = 100;

pub const TAILOR_QUANTILE: f64 = 0.01;

pub fn tailor_denominator<S: std::hash::BuildHasher>(
    hist: &std::collections::HashMap<i32, u32, S>,
    total: u32,
) -> f64 {
    if total < TAILOR_MIN_CANDIDATES {
        return 1.0;
    }
    let top_k = (TAILOR_QUANTILE * total as f64).ceil() as u64;
    let top_k = top_k.max(1);

    let mut scores: Vec<i32> = hist.keys().copied().collect();
    scores.sort_unstable_by(|a, b| b.cmp(a));
    let mut cum: u64 = 0;
    let mut quantile_score: Option<i32> = None;
    for s in scores {
        cum += hist[&s] as u64;
        if cum >= top_k {
            quantile_score = Some(s);
            break;
        }
    }
    match quantile_score {
        Some(s) if s > 0 => s as f64,
        _ => 1.0,
    }
}

#[derive(Debug, Clone)]
pub struct PsmMatch {
    pub spectrum_idx: usize,
    pub candidate_idxs: Vec<u32>,
    pub charge_used: u8,
    pub mass_error_ppm: f64,
    pub score: f32,
    pub rank_score: f32,
    pub edge_score: i32,
    pub activation_method: Option<model::activation::ActivationMethod>,
    pub features: PsmFeatures,
    pub isotope_offset: i8,
    pub precursor_mz_override: Option<f64>,
    pub site_evidence: Option<Box<crate::site_evidence::SiteEvidence>>,
    pub tie_key: u64,
}

pub fn tie_key_of(cand: &crate::candidate_gen::Candidate) -> u64 {
    (crate::match_engine::peptidoform_key(&cand.peptide) << 1) | u64::from(cand.is_decoy)
}

#[inline]
fn score_key(x: f32) -> f32 {
    if x.is_nan() {
        f32::NEG_INFINITY
    } else {
        x + 0.0
    }
}

impl PsmMatch {
    pub fn rank_cmp(&self, other: &Self) -> std::cmp::Ordering {
        score_key(self.rank_score).total_cmp(&score_key(other.rank_score))
    }

    pub fn primary_candidate_idx(&self) -> u32 {
        self.candidate_idxs[0]
    }
}

impl PartialEq for PsmMatch {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == std::cmp::Ordering::Equal
    }
}

impl Eq for PsmMatch {}

impl PartialOrd for PsmMatch {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for PsmMatch {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.rank_cmp(other)
            .then_with(|| score_key(self.score).total_cmp(&score_key(other.score)))
            .then_with(|| other.tie_key.cmp(&self.tie_key))
            .then_with(|| other.charge_used.cmp(&self.charge_used))
            .then_with(|| {
                (other.isotope_offset.unsigned_abs(), other.isotope_offset)
                    .cmp(&(self.isotope_offset.unsigned_abs(), self.isotope_offset))
            })
            .then_with(|| {
                other
                    .mass_error_ppm
                    .abs()
                    .total_cmp(&self.mass_error_ppm.abs())
                    .then_with(|| other.mass_error_ppm.total_cmp(&self.mass_error_ppm))
            })
            .then_with(|| other.spectrum_idx.cmp(&self.spectrum_idx))
            .then_with(|| {
                let mz = |p: &Self| p.precursor_mz_override.map(f64::to_bits);
                mz(other).cmp(&mz(self))
            })
            .then_with(|| other.candidate_idxs.cmp(&self.candidate_idxs))
    }
}

fn rank_score_desc(a: &PsmMatch, b: &PsmMatch) -> std::cmp::Ordering {
    b.cmp(a)
}

#[derive(Debug, Clone)]
pub struct TopNQueue {
    capacity: u32,
    heap: BinaryHeap<Reverse<PsmMatch>>,
}

impl TopNQueue {
    pub fn new(capacity: u32) -> Self {
        Self {
            capacity,
            heap: BinaryHeap::with_capacity(capacity as usize),
        }
    }

    pub fn push(&mut self, m: PsmMatch) {
        if self.heap.len() < self.capacity as usize {
            self.heap.push(Reverse(m));
        } else if let Some(Reverse(top)) = self.heap.peek() {
            match m.rank_cmp(top) {
                std::cmp::Ordering::Greater => {
                    self.heap.pop();
                    self.heap.push(Reverse(m));
                }
                std::cmp::Ordering::Equal => {
                    self.heap.push(Reverse(m));
                }
                std::cmp::Ordering::Less => {
                }
            }
        }
    }

    pub fn force_push(&mut self, m: PsmMatch) {
        self.heap.push(Reverse(m));
    }

    pub fn len(&self) -> usize {
        self.heap.len()
    }
    pub fn is_empty(&self) -> bool {
        self.heap.is_empty()
    }

    pub fn worst_rank_score(&self) -> Option<f32> {
        self.heap.peek().map(|std::cmp::Reverse(m)| m.rank_score)
    }

    pub fn capacity(&self) -> u32 {
        self.capacity
    }

    pub fn iter_psms(&self) -> impl Iterator<Item = &PsmMatch> {
        self.heap.iter().map(|Reverse(m)| m)
    }

    pub fn drain_into_vec(&mut self) -> Vec<PsmMatch> {
        self.heap.drain().map(|Reverse(m)| m).collect()
    }

    pub fn update_in_place<F: FnMut(&mut PsmMatch)>(&mut self, mut f: F) {
        let mut v = std::mem::take(&mut self.heap).into_vec();
        for Reverse(m) in &mut v {
            f(m);
        }
        self.heap = BinaryHeap::from(v);
    }

    pub fn fill_post_topn<F: FnMut(&mut PsmMatch)>(&mut self, mut f: F) {
        let mut psms: Vec<PsmMatch> = self.heap.drain().map(|Reverse(m)| m).collect();
        for psm in &mut psms {
            f(psm);
        }
        for psm in psms {
            self.heap.push(Reverse(psm));
        }
    }

    pub fn remap_candidate_idxs<F: FnMut(u32) -> u32>(&mut self, mut remap: F) {
        let mut psms: Vec<PsmMatch> = self.heap.drain().map(|Reverse(m)| m).collect();
        for psm in &mut psms {
            for idx in &mut psm.candidate_idxs {
                *idx = remap(*idx);
            }
        }
        for psm in psms {
            self.heap.push(Reverse(psm));
        }
    }

    pub fn reorder_by_strong_score(&mut self) {
        let cap = self.capacity as usize;
        let mut psms: Vec<PsmMatch> = self.heap.drain().map(|Reverse(m)| m).collect();
        for psm in &mut psms {
            psm.rank_score = psm.features.strong_score;
        }
        psms.sort_by(rank_score_desc);
        Self::retain_top_with_ties(&mut psms, cap);
        for psm in psms {
            self.heap.push(Reverse(psm));
        }
    }

    fn retain_top_with_ties(psms: &mut Vec<PsmMatch>, cap: usize) {
        if psms.len() <= cap {
            return;
        }
        if cap == 0 {
            psms.clear();
            return;
        }
        let cutoff = psms[cap - 1].clone();
        psms.retain(|p| p.rank_cmp(&cutoff) != std::cmp::Ordering::Less);
    }

    pub fn trim_to_capacity(&mut self, cap: u32) {
        let cap = cap as usize;
        if self.heap.len() <= cap {
            return;
        }
        let mut psms: Vec<PsmMatch> = self.heap.drain().map(|Reverse(m)| m).collect();
        psms.sort_by(rank_score_desc);
        Self::retain_top_with_ties(&mut psms, cap);
        for psm in psms {
            self.heap.push(Reverse(psm));
        }
    }

    pub fn peek_top(&self) -> Option<&PsmMatch> {
        self.heap.iter().map(|Reverse(m)| m).max_by(|a, b| a.cmp(b))
    }

    pub fn into_sorted_vec(self) -> Vec<PsmMatch> {
        let mut v: Vec<PsmMatch> = self.heap.into_iter().map(|Reverse(m)| m).collect();
        v.sort_by(rank_score_desc);
        v
    }

    pub fn into_rank_sorted_vec(self) -> Vec<PsmMatch> {
        self.into_sorted_vec()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_match(spectrum_idx: usize, score: f32) -> PsmMatch {
        PsmMatch {
            spectrum_idx,
            candidate_idxs: vec![0],
            charge_used: 2,
            mass_error_ppm: 0.0,
            score,
            rank_score: score,
            edge_score: 0,
            activation_method: None,
            features: PsmFeatures::default(),
            isotope_offset: 0,
            precursor_mz_override: None,
            site_evidence: None,
            tie_key: 0,
        }
    }

    #[test]
    fn update_in_place_keeps_the_output_order_of_ties() {
        let mut q = TopNQueue::new(32);
        for i in 0..20usize {
            let mut m = make_match(i, (i % 3) as f32);
            m.candidate_idxs = vec![i as u32];
            q.force_push(m);
        }
        let order = |q: &TopNQueue| -> Vec<u32> {
            q.clone()
                .into_rank_sorted_vec()
                .iter()
                .map(|m| m.candidate_idxs[0])
                .collect()
        };
        let before = order(&q);
        q.update_in_place(|m| m.features.strong_score += 1.0);
        assert_eq!(order(&q), before);
        assert!(q.iter_psms().all(|m| m.features.strong_score == 1.0));
    }

    fn make_match_with_rank(spectrum_idx: usize, score: f32, rank_score: f32) -> PsmMatch {
        let mut m = make_match(spectrum_idx, score);
        m.rank_score = rank_score;
        m
    }

    #[test]
    fn empty_queue() {
        let q = TopNQueue::new(5);
        assert!(q.is_empty());
        assert_eq!(q.len(), 0);
    }

    #[test]
    fn queue_below_capacity_keeps_everything() {
        let mut q = TopNQueue::new(5);
        for s in [1.0, 2.0, 3.0] {
            q.push(make_match(0, s));
        }
        assert_eq!(q.len(), 3);
        let sorted = q.into_sorted_vec();
        assert_eq!(
            sorted.iter().map(|m| m.score).collect::<Vec<_>>(),
            vec![3.0, 2.0, 1.0]
        );
    }

    #[test]
    fn queue_at_capacity_keeps_top_n_by_score() {
        let mut q = TopNQueue::new(3);
        for s in [1.0, 5.0, 2.0, 4.0, 3.0] {
            q.push(make_match(0, s));
        }
        assert_eq!(q.len(), 3);
        let sorted = q.into_sorted_vec();
        assert_eq!(
            sorted.iter().map(|m| m.score).collect::<Vec<_>>(),
            vec![5.0, 4.0, 3.0]
        );
    }

    #[test]
    fn lower_score_dropped_when_full() {
        let mut q = TopNQueue::new(2);
        q.push(make_match(0, 5.0));
        q.push(make_match(0, 3.0));
        assert_eq!(q.len(), 2);
        q.push(make_match(0, 1.0));
        let sorted = q.into_sorted_vec();
        assert_eq!(
            sorted.iter().map(|m| m.score).collect::<Vec<_>>(),
            vec![5.0, 3.0]
        );
    }

    #[test]
    fn topn_queue_keeps_ties_at_capacity() {
        let mut q = TopNQueue::new(1);
        q.push(make_match(0, 100.0));
        q.push(make_match(0, 100.0));
        q.push(make_match(0, 100.0));
        assert_eq!(
            q.len(),
            3,
            "all three tied PSMs should be retained at capacity=1 (tie retention, R-1)"
        );
    }

    #[test]
    fn force_push_bypasses_capacity_eviction() {
        let mut q = TopNQueue::new(1);
        q.push(make_match(0, 100.0));
        let mut secondary = make_match(0, 50.0);
        secondary.candidate_idxs = vec![1];
        q.force_push(secondary);
        let drained = q.drain_into_vec();
        assert_eq!(
            drained.len(),
            2,
            "force_push must retain BOTH primary and secondary despite capacity=1"
        );
    }

    #[test]
    fn dedup_pepseq_score_aggregates_candidate_idxs() {
        let mut q = TopNQueue::new(10);
        let mut a = make_match(0, 50.0);
        a.candidate_idxs = vec![10];
        let mut b = make_match(0, 50.0);
        b.candidate_idxs = vec![20];
        let mut c = make_match(0, 40.0);
        c.candidate_idxs = vec![30];

        q.push(a);
        q.push(b);
        q.push(c);
        assert_eq!(q.len(), 3, "all three PSMs initially retained");

        let drained = q.drain_into_vec();
        assert_eq!(drained.len(), 3);

        let deduped = simple_dedup_by_score_for_test(drained);

        assert_eq!(deduped.len(), 2, "should collapse to 2 unique-score groups");

        let mut score_50 = deduped
            .iter()
            .find(|p| (p.score as i32) == 50)
            .unwrap()
            .candidate_idxs
            .clone();
        score_50.sort();
        assert_eq!(
            score_50,
            vec![10, 20],
            "score=50 should aggregate both idxs"
        );

        let score_40 = &deduped
            .iter()
            .find(|p| (p.score as i32) == 40)
            .unwrap()
            .candidate_idxs;
        assert_eq!(*score_40, vec![30]);
    }

    fn simple_dedup_by_score_for_test(psms: Vec<PsmMatch>) -> Vec<PsmMatch> {
        use std::collections::HashMap;
        let mut groups: HashMap<i32, PsmMatch> = HashMap::new();
        for psm in psms {
            let key = psm.score as i32;
            groups
                .entry(key)
                .and_modify(|existing| {
                    existing
                        .candidate_idxs
                        .extend(psm.candidate_idxs.iter().copied())
                })
                .or_insert(psm);
        }
        groups.into_values().collect()
    }

    #[test]
    fn psm_match_clones_correctly() {
        let m = make_match(7, 4.2);
        let cloned = m.clone();
        assert_eq!(cloned.spectrum_idx, 7);
        assert_eq!(cloned.score, 4.2);
        assert_eq!(cloned.rank_score, 4.2);
    }

    #[test]
    fn psm_match_orders_by_rank_score_descending() {
        let better = make_match_with_rank(0, 5.0, 10.0);
        let worse = make_match_with_rank(0, 5.0, 1.0);
        assert!(
            better > worse,
            "PSM with higher rank_score should be Ord-greater (better in the min-heap)"
        );
    }

    #[test]
    fn queue_keeps_best_rank_score_psms_when_full() {
        let mut q = TopNQueue::new(2);
        q.push(make_match_with_rank(0, 5.0, 1.0));
        q.push(make_match_with_rank(0, 5.0, 100.0));
        assert_eq!(q.len(), 2);
        q.push(make_match_with_rank(0, 5.0, 50.0));
        assert_eq!(q.len(), 2);
        let sorted = q.into_sorted_vec();
        let ranks: Vec<f32> = sorted.iter().map(|m| m.rank_score).collect();
        assert!(ranks.contains(&100.0), "best rank 100.0 should be retained");
        assert!(ranks.contains(&50.0), "medium rank 50.0 should be retained");
        assert!(!ranks.contains(&1.0), "worst rank 1.0 should be evicted");
    }

    #[test]
    fn iter_psms_yields_all_psms() {
        let mut q = TopNQueue::new(5);
        for s in [1.0_f32, 2.0, 3.0] {
            q.push(make_match(0, s));
        }
        let scores: Vec<f32> = {
            let mut v: Vec<f32> = q.iter_psms().map(|m| m.score).collect();
            v.sort_by(|a, b| b.partial_cmp(a).unwrap());
            v
        };
        assert_eq!(scores, vec![3.0, 2.0, 1.0]);
    }

    #[test]
    fn psm_match_default_isotope_offset_is_zero() {
        let m = make_match(0, 1.0);
        assert_eq!(
            m.isotope_offset, 0,
            "isotope_offset sentinel should be 0 before match_engine populates it"
        );
    }

    #[test]
    fn psm_features_default_is_zero() {
        let f = PsmFeatures::default();
        assert_eq!(f.num_matched_main_ions, 0);
        assert_eq!(f.longest_b, 0);
        assert_eq!(f.longest_y, 0);
        assert_eq!(f.longest_y_pct, 0.0);
        assert_eq!(f.matched_ion_ratio, 0.0);
        assert_eq!(f.explained_ion_current_ratio, 0.0);
        assert_eq!(f.n_term_ion_current_ratio, 0.0);
        assert_eq!(f.c_term_ion_current_ratio, 0.0);
        assert_eq!(f.ms2_ion_current, 0.0);
        assert_eq!(f.isolation_window_efficiency, 0.0);
        assert_eq!(f.mean_error_top7, 0.0);
        assert_eq!(f.stdev_error_top7, 0.0);
        assert_eq!(f.mean_rel_error_top7, 0.0);
        assert_eq!(f.stdev_rel_error_top7, 0.0);
        assert_eq!(f.strong_score, 0.0);
        assert_eq!(f.strong_score_cal, 0.0);
    }

    #[test]
    fn trim_to_zero_capacity_empties_the_queue() {
        let mut q = TopNQueue::new(10);
        for s in [5.0, 4.0, 3.0] {
            q.push(make_match(0, s));
        }
        q.trim_to_capacity(0);
        assert!(q.is_empty());
    }

    #[test]
    fn trim_to_capacity_keeps_ties_at_cutoff() {
        let mut q = TopNQueue::new(10);
        for s in [5.0, 4.0, 3.0, 3.0, 2.0] {
            q.push(make_match(0, s));
        }
        q.trim_to_capacity(3);
        let scores: Vec<f32> = q.into_sorted_vec().iter().map(|m| m.score).collect();
        assert_eq!(scores, vec![5.0, 4.0, 3.0, 3.0]);
    }

    #[test]
    fn reorder_by_strong_score_updates_ranking_fields() {
        let mut q = TopNQueue::new(2);
        let low = PsmMatch {
            spectrum_idx: 0,
            candidate_idxs: vec![0],
            charge_used: 2,
            mass_error_ppm: 0.0,
            score: 50.0,
            rank_score: 50.0,
            edge_score: 0,
            activation_method: None,
            features: PsmFeatures {
                strong_score: 1.0,
                ..Default::default()
            },
            isotope_offset: 0,
            precursor_mz_override: None,
            site_evidence: None,
            tie_key: 0,
        };
        let mut high = low.clone();
        high.candidate_idxs = vec![1];
        high.score = 100.0;
        high.rank_score = 100.0;
        high.features.strong_score = 5.0;
        q.push(low);
        q.push(high);
        q.reorder_by_strong_score();
        let top = q.peek_top().unwrap();
        assert!((top.features.strong_score - 5.0).abs() < f32::EPSILON);
        assert!((top.score - 100.0).abs() < f32::EPSILON);
        assert!((top.rank_score - 5.0).abs() < f32::EPSILON);
    }

    #[test]
    fn psm_match_default_features_is_zeroed() {
        let m = make_match(0, 1.0);
        assert_eq!(
            m.features.num_matched_main_ions, 0,
            "features.num_matched_main_ions should default to 0"
        );
        assert_eq!(
            m.features.longest_b, 0,
            "features.longest_b should default to 0"
        );
        assert_eq!(
            m.features.longest_y, 0,
            "features.longest_y should default to 0"
        );
        assert_eq!(
            m.features.longest_y_pct, 0.0,
            "features.longest_y_pct should default to 0.0"
        );
        assert_eq!(
            m.features.matched_ion_ratio, 0.0,
            "features.matched_ion_ratio should default to 0.0"
        );
        assert_eq!(
            m.features.explained_ion_current_ratio, 0.0,
            "explained_ion_current_ratio should default to 0.0"
        );
        assert_eq!(
            m.features.n_term_ion_current_ratio, 0.0,
            "n_term_ion_current_ratio should default to 0.0"
        );
        assert_eq!(
            m.features.c_term_ion_current_ratio, 0.0,
            "c_term_ion_current_ratio should default to 0.0"
        );
        assert_eq!(
            m.features.ms2_ion_current, 0.0,
            "ms2_ion_current should default to 0.0"
        );
        assert_eq!(
            m.features.isolation_window_efficiency, 0.0,
            "isolation_window_efficiency should default to 0.0"
        );
        assert_eq!(
            m.features.mean_error_top7, 0.0,
            "mean_error_top7 should default to 0.0"
        );
        assert_eq!(
            m.features.stdev_error_top7, 0.0,
            "stdev_error_top7 should default to 0.0"
        );
        assert_eq!(
            m.features.mean_rel_error_top7, 0.0,
            "mean_rel_error_top7 should default to 0.0"
        );
        assert_eq!(
            m.features.stdev_rel_error_top7, 0.0,
            "stdev_rel_error_top7 should default to 0.0"
        );
    }

    #[test]
    fn psm_match_with_nan_rank_score_orders_as_worst() {
        let nan_rank = make_match_with_rank(0, 5.0, f32::NAN);
        let finite = make_match_with_rank(0, 0.0, 1.0);
        assert_eq!(
            nan_rank.cmp(&finite),
            std::cmp::Ordering::Less,
            "NaN rank_score should sort as worse (Less) than a finite value"
        );
    }

    use std::collections::HashMap;

    fn hist_from(pairs: &[(i32, u32)]) -> (HashMap<i32, u32>, u32) {
        let mut h = HashMap::new();
        let mut total = 0u32;
        for &(s, c) in pairs {
            *h.entry(s).or_insert(0) += c;
            total += c;
        }
        (h, total)
    }

    #[test]
    fn tailor_denom_top1pct_quantile_basic() {
        let (h, total) = hist_from(&[(5, 990), (50, 10)]);
        assert_eq!(total, 1000);
        assert_eq!(tailor_denominator(&h, total), 50.0);
        assert!(((50.0f64 / tailor_denominator(&h, total)) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn tailor_denom_ceil_rounds_up() {
        let (h, total) = hist_from(&[(100, 1), (90, 1), (10, 148)]);
        assert_eq!(total, 150);
        assert_eq!(tailor_denominator(&h, total), 90.0);
    }

    #[test]
    fn tailor_denom_small_n_falls_back_to_one() {
        let (h, total) = hist_from(&[(5, 50), (50, 5)]);
        assert_eq!(total, 55);
        assert!(total < TAILOR_MIN_CANDIDATES);
        assert_eq!(tailor_denominator(&h, total), 1.0);
    }

    #[test]
    fn tailor_denom_nonpositive_quantile_falls_back_to_one() {
        let (h, total) = hist_from(&[(-10, 990), (-1, 10)]);
        assert_eq!(total, 1000);
        assert_eq!(tailor_denominator(&h, total), 1.0);
    }

    #[test]
    fn tailor_denom_exactly_min_candidates_calibrates() {
        let (h, total) = hist_from(&[(10, 99), (40, 1)]);
        assert_eq!(total, 100);
        assert_eq!(tailor_denominator(&h, total), 40.0);
    }

    #[test]
    fn psm_match_two_nan_rank_scores_compare_equal() {
        let a = make_match_with_rank(0, 5.0, f32::NAN);
        let b = make_match_with_rank(0, 5.0, f32::NAN);
        assert_eq!(
            a.cmp(&b),
            std::cmp::Ordering::Equal,
            "Two PSMs with NaN rank_score should compare Equal"
        );
    }

    fn tied(tie_key: u64, cand: u32) -> PsmMatch {
        let mut m = make_match(0, 10.0);
        m.tie_key = tie_key;
        m.candidate_idxs = vec![cand];
        m
    }

    fn keys(v: &[PsmMatch]) -> Vec<u64> {
        v.iter().map(|m| m.tie_key).collect()
    }

    #[test]
    fn exact_ties_order_by_content_not_by_push_order() {
        let psms: Vec<PsmMatch> = [(30, 0), (10, 2), (20, 1), (40, 3)]
            .into_iter()
            .map(|(k, c)| tied(k, c))
            .collect();
        let mut outputs = Vec::new();
        for rotation in 0..psms.len() {
            let mut q = TopNQueue::new(10);
            for m in psms.iter().cycle().skip(rotation).take(psms.len()) {
                q.push(m.clone());
            }
            assert_eq!(q.peek_top().unwrap().tie_key, 10);
            outputs.push(keys(&q.into_rank_sorted_vec()));
        }
        assert!(
            outputs.iter().all(|o| *o == vec![10, 20, 30, 40]),
            "{outputs:?}"
        );
    }

    #[test]
    fn rank_ties_are_still_retained_past_capacity() {
        let mut q = TopNQueue::new(1);
        q.push(tied(2, 0));
        q.push(tied(1, 1));
        assert_eq!(q.len(), 2, "a rank_score tie at capacity is kept");
        let mut better = tied(3, 2);
        better.rank_score = 11.0;
        q.push(better);
        assert_eq!(keys(&q.into_rank_sorted_vec()), vec![3, 1]);
    }

    #[test]
    fn total_order_breaks_ties_in_the_documented_order() {
        use std::cmp::Ordering::*;
        let base = tied(5, 0);
        let mut higher_score = base.clone();
        higher_score.score = 11.0;
        assert_eq!(higher_score.cmp(&base), Greater, "score breaks a rank tie");
        let mut lower_charge = base.clone();
        lower_charge.charge_used = 1;
        assert_eq!(lower_charge.cmp(&base), Greater);
        let mut iso = base.clone();
        iso.isotope_offset = 1;
        assert_eq!(base.cmp(&iso), Greater, "isotope offset 0 beats 1");
        let mut err = base.clone();
        err.mass_error_ppm = 3.0;
        assert_eq!(base.cmp(&err), Greater, "smaller |mass error| wins");
        assert_eq!(
            tied(4, 9).cmp(&base),
            Greater,
            "smaller tie_key wins over candidate index"
        );
        let mut neg_zero = base.clone();
        neg_zero.rank_score = -0.0;
        let mut pos_zero = base.clone();
        pos_zero.rank_score = 0.0;
        assert_eq!(neg_zero.rank_cmp(&pos_zero), Equal);
        let mut nan = base.clone();
        nan.rank_score = f32::NAN;
        assert_eq!(nan.rank_cmp(&neg_zero), Less, "NaN ranks worst");
    }
}
