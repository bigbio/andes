//! PSM (peptide-spectrum match) data + top-N ranking queue.

use std::cmp::Reverse;
use std::collections::BinaryHeap;


/// Per-PSM fragment-ion feature columns computed from the scoring machinery
/// and emitted into the Percolator `.pin` file.
///
/// Filled by `compute_psm_features` in `match_engine.rs` after `score_psm`.
/// Fields use `Default` (all zero) as the safe sentinel before computation.
#[derive(Debug, Clone, Default)]
pub struct PsmFeatures {
    /// Number of unique fragment positions where a b- or y-ion at charge 1
    /// matched a peak within the fragment tolerance. Each position counts
    /// at most once per ion series, but can contribute 1 from b AND 1 from y.
    pub num_matched_main_ions: u32,
    /// Length of the longest contiguous run of matched b-ions
    /// (b1, b2, … must all match to form the run).
    pub longest_b: u32,
    /// Length of the longest contiguous run of matched y-ions.
    pub longest_y: u32,
    /// `longest_y as f32 / peptide.length() as f32` — fraction in 0.0..=1.0.
    pub longest_y_pct: f32,
    /// `num_matched_main_ions as f32 / peptide.length() as f32` — fraction
    /// of peptide positions covered by matched b/y ions.
    pub matched_ion_ratio: f32,

    // ── Ion-current ratios ─────────────────────────────────────────────────

    /// `n_term_ion_current_ratio + c_term_ion_current_ratio`.
    pub explained_ion_current_ratio: f32,
    /// Sum of matched b-ion intensities divided by total MS2 ion current.
    pub n_term_ion_current_ratio: f32,
    /// Sum of matched y-ion intensities divided by total MS2 ion current.
    pub c_term_ion_current_ratio: f32,
    /// Raw sum of all peak intensities in the MS2 spectrum (no log10).
    pub ms2_ion_current: f32,
    /// Isolation-window efficiency. Not available from the Spectrum object;
    /// always emitted as 0.0.
    pub isolation_window_efficiency: f32,

    // ── Top-7 mass-error statistics ────────────────────────────────────────

    /// Mean of absolute Da errors for the top-7 most-intense matched ions.
    pub mean_error_top7: f32,
    /// Population standard deviation of absolute Da errors for top-7 ions
    /// (formula: `sqrt(E[x²] - mean²)`).
    pub stdev_error_top7: f32,
    /// Mean of signed relative errors (ppm) for the top-7 most-intense matched ions.
    pub mean_rel_error_top7: f32,
    /// Population standard deviation of signed relative errors (ppm) for top-7 ions.
    pub stdev_rel_error_top7: f32,
}

#[derive(Debug, Clone)]
pub struct PsmMatch {
    pub spectrum_idx: usize,
    /// Index into the `&[Candidate]` slice owned by `PreparedSearch.candidates`.
    /// Replaces the inlined `Candidate` clone: previously each push to the queue
    /// cloned the full `Candidate` (including its `Peptide.residues: Vec<...>`),
    /// allocating millions of times per large-fasta search. Now the queue stores
    /// only a 4-byte index and consumers (writers, feature extraction, GF) look
    /// up the `Candidate` by index when needed.
    ///
    /// Every real PSM points at a valid index into `PreparedSearch.candidates`.
    /// There is no "synthetic / no backing Candidate" sentinel — test fixtures
    /// that don't need to resolve back use `0` as a placeholder and avoid
    /// touching the candidates slice from inside the test.
    pub candidate_idx: u32,
    pub charge_used: u8,
    /// Signed: positive when peptide mass exceeds spectrum's implied mass.
    pub mass_error_ppm: f64,
    /// Higher is better. Real spectral-similarity score.
    pub score: f32,
    /// SpecEValue: lower is better. Default 1.0 = "not yet computed"
    /// / "no signal". Set by `compute_spec_e_values_for_spectrum` after the
    /// per-candidate scoring loop.
    pub spec_e_value: f64,
    /// De-novo score: `gf_group.max_score() - 1` for the GF that scored
    /// this peptide. Set during `compute_spec_e_values_for_spectrum`.
    /// Sentinel: `i32::MIN` if not yet computed.
    pub de_novo_score: i32,
    /// Activation method captured from `param.data_type.activation` at scoring
    /// time. `None` if unknown or not yet set.
    pub activation_method: Option<model::activation::ActivationMethod>,
    /// `spec_e_value * num_distinct_peptides_at_length`. Sentinel: `1.0`.
    /// Approximate: uses the candidate-set size filtered by the same length as
    /// a proxy for `num_distinct_peptides` when no suffix-array helper exists.
    pub e_value: f64,
    /// Fragment-ion feature columns computed after `score_psm`.
    /// Defaults to all-zero until `compute_psm_features` runs.
    pub features: PsmFeatures,
    /// The isotope offset that produced the precursor match: 0 = monoisotopic,
    /// +N = spectrum precursor was N C13 peaks above the true monoisotopic.
    /// Default range −1..=2. Threaded from `MassError::isotope_offset`
    /// (precursor_matching.rs) via match_engine.rs. Written as the PIN
    /// `isotope_error` column.
    pub isotope_offset: i8,
}

impl PartialEq for PsmMatch {
    fn eq(&self, other: &Self) -> bool {
        self.spec_e_value == other.spec_e_value && self.score == other.score
    }
}

impl Eq for PsmMatch {}

impl PartialOrd for PsmMatch {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// Primary: `spec_e_value` ascending (lower = better).
/// Secondary: `score` descending (higher = better).
///
/// This ordering is used by `TopNQueue`'s min-heap (via `Reverse<PsmMatch>`):
/// the heap's "minimum" element is the one with the *largest* spec_e_value
/// (worst), so `push` evicts it when over capacity.
impl Ord for PsmMatch {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        use std::cmp::Ordering;
        // "Better" PSM = smaller spec_e_value, then larger score.
        // NaN spec_e_value or score is treated as worst (sorts last / loses to finite).
        // Map NaN spec_e_value → +infinity (worst, since smaller is better).
        // Map NaN score        → -infinity (worst, since larger is better).
        let self_sev  = if self.spec_e_value.is_nan()  { f64::INFINITY }      else { self.spec_e_value };
        let other_sev = if other.spec_e_value.is_nan() { f64::INFINITY }      else { other.spec_e_value };
        match other_sev.partial_cmp(&self_sev).unwrap_or(Ordering::Equal) {
            Ordering::Equal => {
                let self_score  = if self.score.is_nan()  { f32::NEG_INFINITY } else { self.score };
                let other_score = if other.score.is_nan() { f32::NEG_INFINITY } else { other.score };
                self_score.partial_cmp(&other_score).unwrap_or(Ordering::Equal)
            }
            ord => ord,
        }
    }
}

#[derive(Debug, Clone)]
pub struct TopNQueue {
    capacity: u32,
    /// Min-heap (via Reverse): smallest score sits at top, easy to pop
    /// when over capacity.
    heap: BinaryHeap<Reverse<PsmMatch>>,
}

impl TopNQueue {
    pub fn new(capacity: u32) -> Self {
        Self { capacity, heap: BinaryHeap::with_capacity(capacity as usize) }
    }

    /// Insert a PSM. The queue keeps **at least** `capacity` of the *best*
    /// PSMs, plus any additional PSMs tied with the current worst.
    ///
    /// "Best" = smallest `spec_e_value` first (then largest `score` for ties).
    /// The min-heap (via `Reverse<PsmMatch>`) puts the *worst* PSM at the top
    /// so it can be evicted when a strictly-better PSM arrives.
    ///
    /// Before `compute_spec_e_values_for_spectrum` runs, all PSMs have
    /// `spec_e_value = 1.0` and the secondary `score` key governs eviction.
    ///
    /// **Tie handling (R-1, 2026-05-18):** when the queue is at capacity and
    /// a new PSM is `Equal` (in `Ord` terms) to the worst retained PSM, the
    /// new PSM is inserted WITHOUT evicting the tied one. This matches
    /// Java's `DBScanner.java:540` (`size < n OR score == worst → add`).
    /// As a result, the queue can grow beyond `capacity` when ties exist;
    /// `capacity` becomes a *minimum* top-N, not a hard cap.
    pub fn push(&mut self, m: PsmMatch) {
        if self.heap.len() < self.capacity as usize {
            self.heap.push(Reverse(m));
        } else if let Some(Reverse(top)) = self.heap.peek() {
            match m.cmp(top) {
                std::cmp::Ordering::Greater => {
                    // m is strictly better than the worst retained PSM: evict
                    // the worst, insert m.
                    self.heap.pop();
                    self.heap.push(Reverse(m));
                }
                std::cmp::Ordering::Equal => {
                    // R-1 (2026-05-18): Java's DBScanner.java:540 keeps tied
                    // PSMs at capacity (and DBScanner.java:745 keeps SpecE
                    // ties on the per-spectrum merge). Rust now matches.
                    // The queue may exceed `capacity` when ties exist —
                    // `capacity` becomes a *minimum* top-N, not a hard cap.
                    // Spec:
                    // docs/parity-analysis/specs/2026-05-18-r1-tie-retention-test-design.md
                    self.heap.push(Reverse(m));
                }
                std::cmp::Ordering::Less => {
                    // m is strictly worse than the worst retained PSM: drop.
                }
            }
        }
    }

    pub fn len(&self) -> usize { self.heap.len() }
    pub fn is_empty(&self) -> bool { self.heap.is_empty() }

    /// Iterate over all PSMs in the queue (order not guaranteed).
    pub fn iter_psms(&self) -> impl Iterator<Item = &PsmMatch> {
        self.heap.iter().map(|Reverse(m)| m)
    }

    /// Apply `f` to each retained PSM in-place. Used for filling in
    /// post-finalization fields (e.g. `features`) that are NOT part of
    /// `PsmMatch::cmp` and therefore do not affect heap ordering.
    ///
    /// Implementation drains the heap, applies `f`, and re-pushes — this is
    /// O(N log N) on a small `N` (top-N, typically 1-10) and avoids the
    /// std-library restriction that `BinaryHeap::iter_mut()` is not exposed
    /// (it would let callers break the heap invariant). Since features do
    /// not participate in ordering, the re-push is logically a no-op for
    /// retention.
    ///
    /// This is distinct from `update_psm_enrichment` only in intent
    /// (post-top-N feature fill vs Phase-7 score/e-value enrichment) — the
    /// mechanism is identical.
    pub fn fill_post_topn<F: FnMut(&mut PsmMatch)>(&mut self, mut f: F) {
        let mut psms: Vec<PsmMatch> = self.heap.drain().map(|Reverse(m)| m).collect();
        for psm in &mut psms {
            f(psm);
        }
        for psm in psms {
            self.heap.push(Reverse(psm));
        }
    }

    /// Return the best PSM (smallest `spec_e_value`, then largest `score`)
    /// without removing it. Returns `None` if the queue is empty.
    ///
    /// The heap is a min-heap on `Reverse<PsmMatch>` so the *worst* entry sits
    /// at the top (for cheap eviction). To find the *best* entry we iterate
    /// all elements and take the max in natural `PsmMatch` ordering.
    /// Cost is O(N) — acceptable for the small top-N queues used in practice.
    pub fn peek_top(&self) -> Option<&PsmMatch> {
        self.heap.iter().map(|Reverse(m)| m).max_by(|a, b| a.cmp(b))
    }

    /// Apply `f` to each PSM to compute its `spec_e_value`, then rebuild
    /// the heap so the ordering invariant holds.
    ///
    /// Draining + re-inserting is O(N log N) — cheap for small N (top-10).
    pub fn update_spec_e_values<F: Fn(&PsmMatch) -> f64>(&mut self, f: F) {
        let mut psms: Vec<PsmMatch> = self.heap.drain().map(|Reverse(m)| m).collect();
        for psm in &mut psms {
            psm.spec_e_value = f(psm);
        }
        for psm in psms {
            self.heap.push(Reverse(psm));
        }
    }

    /// Apply `f` to each PSM in-place (mutable borrow), then rebuild the heap.
    ///
    /// Used by enrichment to set `de_novo_score`, `e_value`, and other
    /// fields that don't affect ordering. The heap is rebuilt after all mutations
    /// (O(N) heapify) to maintain the invariant.
    pub fn update_psm_enrichment<F: FnMut(&mut PsmMatch)>(&mut self, mut f: F) {
        let mut psms: Vec<PsmMatch> = self.heap.drain().map(|Reverse(m)| m).collect();
        for psm in &mut psms {
            f(psm);
        }
        for psm in psms {
            self.heap.push(Reverse(psm));
        }
    }

    /// Drain into a Vec sorted best-first (smallest spec_e_value, then largest score).
    pub fn into_sorted_vec(self) -> Vec<PsmMatch> {
        let mut v: Vec<PsmMatch> = self.heap.into_iter().map(|Reverse(m)| m).collect();
        v.sort_by(|a, b| b.cmp(a));
        v
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_match(spectrum_idx: usize, score: f32) -> PsmMatch {
        // Test-only PSM: candidate_idx = 0 is a sentinel for queue-ordering tests
        // that never resolve back to a real Candidate. Tests that need to read
        // peptide / protein metadata must build their own &[Candidate] alongside.
        PsmMatch {
            spectrum_idx,
            candidate_idx: 0,
            charge_used: 2,
            mass_error_ppm: 0.0,
            score,
            spec_e_value: 1.0,  // default sentinel: "not yet computed"
            de_novo_score: i32::MIN,  // sentinel: not yet computed
            activation_method: None,
            e_value: 1.0,  // sentinel: not yet computed
            features: PsmFeatures::default(),
            isotope_offset: 0,
        }
    }

    fn make_match_with_evalue(spectrum_idx: usize, score: f32, spec_e_value: f64) -> PsmMatch {
        let mut m = make_match(spectrum_idx, score);
        m.spec_e_value = spec_e_value;
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
        for s in [1.0, 2.0, 3.0] { q.push(make_match(0, s)); }
        assert_eq!(q.len(), 3);
        let sorted = q.into_sorted_vec();
        // All spec_e_value = 1.0 (default) → secondary sort by score descending.
        assert_eq!(sorted.iter().map(|m| m.score).collect::<Vec<_>>(),
                   vec![3.0, 2.0, 1.0]);
    }

    #[test]
    fn queue_at_capacity_keeps_top_n_by_score() {
        let mut q = TopNQueue::new(3);
        for s in [1.0, 5.0, 2.0, 4.0, 3.0] { q.push(make_match(0, s)); }
        assert_eq!(q.len(), 3);
        let sorted = q.into_sorted_vec();
        // All spec_e_value = 1.0 → secondary score keeps top-3 by score.
        assert_eq!(sorted.iter().map(|m| m.score).collect::<Vec<_>>(),
                   vec![5.0, 4.0, 3.0]);
    }

    #[test]
    fn lower_score_dropped_when_full() {
        let mut q = TopNQueue::new(2);
        q.push(make_match(0, 5.0));
        q.push(make_match(0, 3.0));
        assert_eq!(q.len(), 2);
        q.push(make_match(0, 1.0));
        let sorted = q.into_sorted_vec();
        assert_eq!(sorted.iter().map(|m| m.score).collect::<Vec<_>>(),
                   vec![5.0, 3.0]);
    }

    #[test]
    fn topn_queue_keeps_ties_at_capacity() {
        // R-1 fix: Java's DBScanner keeps tied PSMs at capacity
        // (DBScanner.java:540 raw-score retention; DBScanner.java:745 SpecE
        // merge). Rust's TopNQueue must mirror this — strict-greater eviction
        // was dropping ties Java keeps, plausibly causing the Astral 14K raw-
        // target gap. See
        // docs/parity-analysis/notes/2026-05-18-piecewise-fixes-dont-work.md
        // (Open: retention layer, R-1).
        let mut q = TopNQueue::new(1);
        q.push(make_match(0, 100.0));
        q.push(make_match(0, 100.0));
        q.push(make_match(0, 100.0));
        assert_eq!(
            q.len(),
            3,
            "all three tied PSMs should be retained at capacity=1 (Java parity, R-1)"
        );
    }

    #[test]
    fn psm_match_clones_correctly() {
        let m = make_match(7, 4.2);
        let cloned = m.clone();
        assert_eq!(cloned.spectrum_idx, 7);
        assert_eq!(cloned.score, 4.2);
        assert_eq!(cloned.spec_e_value, 1.0);
    }

    // -----------------------------------------------------------------------
    // SpecEValue ordering tests
    // -----------------------------------------------------------------------

    #[test]
    fn psm_match_orders_by_spec_e_value_ascending_then_score_descending() {
        // Lower spec_e_value means "better" → should sort before (greater in
        // natural Ord so the min-heap can evict the worst).
        let better = make_match_with_evalue(0, 5.0, 0.001);
        let worse  = make_match_with_evalue(0, 5.0, 0.5);
        // "better" is greater in natural order (because lower e-value wins).
        assert!(better > worse,
            "PSM with lower spec_e_value should be Ord-greater (better in the min-heap)");

        // Tie-break by score descending.
        let high_score = make_match_with_evalue(0, 10.0, 0.01);
        let low_score  = make_match_with_evalue(0, 3.0,  0.01);
        assert!(high_score > low_score,
            "when spec_e_value equal, higher score should be Ord-greater");
    }

    #[test]
    fn queue_keeps_best_spec_e_value_psms_when_full() {
        // Three PSMs with same score but different spec_e_values; capacity = 2.
        let mut q = TopNQueue::new(2);
        q.push(make_match_with_evalue(0, 5.0, 0.5));   // worst
        q.push(make_match_with_evalue(0, 5.0, 0.001)); // best
        assert_eq!(q.len(), 2);
        // Push a medium one; it should evict the worst (0.5).
        q.push(make_match_with_evalue(0, 5.0, 0.1));
        assert_eq!(q.len(), 2);
        let sorted = q.into_sorted_vec();
        // Should keep 0.001 and 0.1 (best two).
        let evalues: Vec<f64> = sorted.iter().map(|m| m.spec_e_value).collect();
        assert!(evalues.contains(&0.001), "best e-value 0.001 should be retained");
        assert!(evalues.contains(&0.1),   "medium e-value 0.1 should be retained");
        assert!(!evalues.contains(&0.5),  "worst e-value 0.5 should be evicted");
    }

    #[test]
    fn update_spec_e_values_applies_to_all_psms() {
        let mut q = TopNQueue::new(5);
        for s in [1.0_f32, 2.0, 3.0] {
            q.push(make_match(0, s));
        }
        // Set spec_e_value = 1.0 / score for each PSM.
        q.update_spec_e_values(|psm| 1.0 / psm.score as f64);
        let sorted = q.into_sorted_vec();
        // After update: score 3.0 → e=0.333, score 2.0 → e=0.5, score 1.0 → e=1.0.
        // Best e-value first.
        assert!((sorted[0].spec_e_value - 1.0 / 3.0).abs() < 1e-9);
        assert!((sorted[1].spec_e_value - 0.5).abs() < 1e-9);
        assert!((sorted[2].spec_e_value - 1.0).abs() < 1e-9);
    }

    #[test]
    fn iter_psms_yields_all_psms() {
        let mut q = TopNQueue::new(5);
        for s in [1.0_f32, 2.0, 3.0] { q.push(make_match(0, s)); }
        let scores: Vec<f32> = {
            let mut v: Vec<f32> = q.iter_psms().map(|m| m.score).collect();
            v.sort_by(|a, b| b.partial_cmp(a).unwrap());
            v
        };
        assert_eq!(scores, vec![3.0, 2.0, 1.0]);
    }

    // -----------------------------------------------------------------------
    // isotope_offset field
    // -----------------------------------------------------------------------

    #[test]
    fn psm_match_default_isotope_offset_is_zero() {
        let m = make_match(0, 1.0);
        assert_eq!(m.isotope_offset, 0,
            "isotope_offset sentinel should be 0 before match_engine populates it");
    }

    // -----------------------------------------------------------------------
    // Enrichment field sentinel defaults
    // -----------------------------------------------------------------------

    #[test]
    fn psm_match_default_de_novo_score_is_min() {
        let m = make_match(0, 1.0);
        assert_eq!(m.de_novo_score, i32::MIN,
            "de_novo_score sentinel should be i32::MIN before enrichment");
    }

    #[test]
    fn psm_match_default_e_value_is_one() {
        let m = make_match(0, 1.0);
        assert_eq!(m.e_value, 1.0,
            "e_value sentinel should be 1.0 before enrichment");
    }

    // -----------------------------------------------------------------------
    // PsmFeatures struct and default initialization
    // -----------------------------------------------------------------------

    #[test]
    fn psm_features_default_is_zero() {
        let f = PsmFeatures::default();
        assert_eq!(f.num_matched_main_ions, 0);
        assert_eq!(f.longest_b, 0);
        assert_eq!(f.longest_y, 0);
        assert_eq!(f.longest_y_pct, 0.0);
        assert_eq!(f.matched_ion_ratio, 0.0);
        // Ion-current + error-stat columns (9 fields)
        assert_eq!(f.explained_ion_current_ratio, 0.0);
        assert_eq!(f.n_term_ion_current_ratio, 0.0);
        assert_eq!(f.c_term_ion_current_ratio, 0.0);
        assert_eq!(f.ms2_ion_current, 0.0);
        assert_eq!(f.isolation_window_efficiency, 0.0);
        assert_eq!(f.mean_error_top7, 0.0);
        assert_eq!(f.stdev_error_top7, 0.0);
        assert_eq!(f.mean_rel_error_top7, 0.0);
        assert_eq!(f.stdev_rel_error_top7, 0.0);
    }

    #[test]
    fn psm_match_default_features_is_zeroed() {
        let m = make_match(0, 1.0);
        assert_eq!(m.features.num_matched_main_ions, 0,
            "features.num_matched_main_ions should default to 0");
        assert_eq!(m.features.longest_b, 0,
            "features.longest_b should default to 0");
        assert_eq!(m.features.longest_y, 0,
            "features.longest_y should default to 0");
        assert_eq!(m.features.longest_y_pct, 0.0,
            "features.longest_y_pct should default to 0.0");
        assert_eq!(m.features.matched_ion_ratio, 0.0,
            "features.matched_ion_ratio should default to 0.0");
        // Ion-current + error-stat columns (9 fields)
        assert_eq!(m.features.explained_ion_current_ratio, 0.0,
            "explained_ion_current_ratio should default to 0.0");
        assert_eq!(m.features.n_term_ion_current_ratio, 0.0,
            "n_term_ion_current_ratio should default to 0.0");
        assert_eq!(m.features.c_term_ion_current_ratio, 0.0,
            "c_term_ion_current_ratio should default to 0.0");
        assert_eq!(m.features.ms2_ion_current, 0.0,
            "ms2_ion_current should default to 0.0");
        assert_eq!(m.features.isolation_window_efficiency, 0.0,
            "isolation_window_efficiency should default to 0.0");
        assert_eq!(m.features.mean_error_top7, 0.0,
            "mean_error_top7 should default to 0.0");
        assert_eq!(m.features.stdev_error_top7, 0.0,
            "stdev_error_top7 should default to 0.0");
        assert_eq!(m.features.mean_rel_error_top7, 0.0,
            "mean_rel_error_top7 should default to 0.0");
        assert_eq!(m.features.stdev_rel_error_top7, 0.0,
            "stdev_rel_error_top7 should default to 0.0");
    }

    // -----------------------------------------------------------------------
    // Issue 8: NaN-safe Ord impl
    // -----------------------------------------------------------------------

    #[test]
    fn psm_match_with_nan_spec_evalue_orders_as_worst() {
        // NaN spec_e_value should sort as WORSE than any finite value.
        // "Better" = greater in natural Ord (used by the min-heap via Reverse).
        let nan_sev = make_match_with_evalue(0, 5.0, f64::NAN);
        let finite  = make_match_with_evalue(0, 0.0, 1.0);
        assert_eq!(
            nan_sev.cmp(&finite),
            std::cmp::Ordering::Less,
            "NaN spec_e_value should sort as worse (Less) than a finite value"
        );
    }

    #[test]
    fn psm_match_with_nan_score_orders_as_worst() {
        // When spec_e_value ties, NaN score should sort as worse than any finite score.
        let nan_score     = make_match_with_evalue(0, f32::NAN, 0.01);
        let finite_score  = make_match_with_evalue(0, 0.0,      0.01);
        assert_eq!(
            nan_score.cmp(&finite_score),
            std::cmp::Ordering::Less,
            "NaN score should sort as worse (Less) than a finite score at equal spec_e_value"
        );
    }

    #[test]
    fn psm_match_two_nan_spec_evalues_compare_equal() {
        // Two PSMs both with NaN spec_e_value and same score → Equal.
        let a = make_match_with_evalue(0, 5.0, f64::NAN);
        let b = make_match_with_evalue(0, 5.0, f64::NAN);
        assert_eq!(
            a.cmp(&b),
            std::cmp::Ordering::Equal,
            "Two PSMs with NaN spec_e_value and equal score should compare Equal"
        );
    }
}
