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

    // ── Additive Java-parity features ──────────────────────────────────────
    /// Per-bond edge score sum, mirroring Java's `DBScanScorer.getScore`
    /// edge loop (IES + error_score per bond). Emitted as a NEW `EdgeScore`
    /// PIN column alongside the unchanged `RawScore`, so Percolator can
    /// learn weights without disrupting the existing RawScore distribution
    /// (which destroyed discrimination in iter17/iter18 when blended into
    /// RawScore directly). Computed via `psm_edge_score` in `score_psm.rs`.
    pub edge_score: i32,

    // ── Chimeric MS1 precursor-envelope features (Task 3) ──────────────────
    /// KL-divergence between the observed precursor isotope envelope (from the
    /// linked MS1) and the averagine theoretical envelope. Higher = poorer match
    /// = likely spurious co-isolation. 0.0 when MS1/feature unavailable.
    pub precursor_isotope_kl: f32,
    /// Observed monoisotopic precursor intensity / median MS1 intensity (SNR
    /// proxy). 0.0 when unavailable.
    pub precursor_snr: f32,

    // ── Chimeric Phase 3 shared-fragment competition features ──────────────
    /// Number of this peptide's matched charge-1 b/y peaks that were NOT already
    /// claimed by a more-confident peptide on the same scan (greedy, rank-1
    /// first). For rank-1 (and the entire non-chimeric path) this equals the
    /// full matched-peak count. Emitted as the additive `UniqueMatchedIons`
    /// PIN column. See `shared_fragment.rs`.
    pub unique_matched_ions: u32,
    /// Σ intensity(unique matched peaks) / Σ intensity(all matched peaks of this
    /// PSM), in 0.0..=1.0. 1.0 for rank-1 / non-chimeric. Additive
    /// `UniqueExplainedFraction` PIN column.
    pub unique_explained_fraction: f32,
    /// Fraction of this peptide's matched peaks already claimed by a
    /// more-confident peptide: |matched ∩ claimed| / |matched|, in 0.0..=1.0.
    /// 0.0 for rank-1 / non-chimeric. Additive `SharedFracClaimed` PIN column.
    pub shared_frac_claimed: f32,
}

#[derive(Debug, Clone)]
pub struct PsmMatch {
    pub spectrum_idx: usize,
    /// Indices into the `&[Candidate]` slice owned by `PreparedSearch.candidates`.
    /// Length is always ≥ 1. The first index (`candidate_idxs[0]`) is the
    /// "primary" candidate — used by callers that need a single Candidate
    /// (most do; see `primary_candidate_idx()`). Multiple indices accumulate
    /// when the R-2 pepSeq+score dedup pass merges multiple Candidates that
    /// share the same peptide sequence and rounded score (typically the same
    /// peptide matched against multiple proteins, e.g. shared tryptic
    /// peptides in target+decoy concat). The PIN writer iterates this Vec to
    /// emit one tab-separated `Proteins` column per row, matching Java parity
    /// for the Proteins column in PIN output.
    ///
    /// Every real PSM has length ≥ 1 with valid indices into
    /// `PreparedSearch.candidates`. Test fixtures that don't need to resolve
    /// back use `vec![0]` as a placeholder and avoid touching the candidates
    /// slice from inside the test.
    pub candidate_idxs: Vec<u32>,
    pub charge_used: u8,
    /// Signed: positive when peptide mass exceeds spectrum's implied mass.
    pub mass_error_ppm: f64,
    /// Pin RawScore = `node_score + cleavage_credit`. Higher is better.
    /// This is what gets emitted in the `RawScore` PIN column (unchanged
    /// from iter19's design). Used by Percolator as one of many features.
    pub score: f32,
    /// iter33: queue-ordering score = `node + cleavage + edge`. Java's
    /// `DBScanScorer.getScore` returns `node + edge` and Java parity adds
    /// cleavage, so Java's `match.score` (used by its `PriorityQueue`
    /// ordering) is `node + cleavage + edge`. Rust's pin RawScore stays at
    /// `node + cleavage` for Percolator distribution stability (iter19); the
    /// SEPARATE `EdgeScore` PIN column carries the `+edge` contribution.
    /// `rank_score` mirrors Java's queue-ordering key without changing the
    /// pin RawScore distribution.
    ///
    /// **No automatic default**: PsmMatch does not implement `Default`, and
    /// callers MUST set `rank_score` explicitly. Test fixtures that build
    /// PsmMatch literals should set `rank_score = score` for pre-iter33
    /// behavior (no edge contribution to ranking). The `match_engine.rs`
    /// candidate loop computes `rank_score = score + edge_score as f32`.
    pub rank_score: f32,
    /// Per-PSM edge_score = `psm_edge_score(...)` for this candidate.
    /// Computed at queue-insertion time in `match_engine.rs` and reused by
    /// `compute_psm_features` to populate the iter19 `EdgeScore` PIN column
    /// (avoids the recompute). Default 0 — features extraction will compute
    /// it on the fly if it remains 0 (e.g. for test fixtures).
    pub edge_score: i32,
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
    /// `spec_e_value * num_distinct_peptides_at_length`. Set in
    /// `compute_spec_e_values_for_spectrum` using
    /// `SearchIndex::num_distinct_peptides_at_length` (counts distinct bare
    /// residue sequences at that length over the enumerated candidate set).
    /// Sentinel before enrichment: `1.0`.
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
    /// Chimeric cascade Pass-2 secondaries only: the co-isolated precursor m/z
    /// this PSM was actually scored at (the SECONDARY peptide's precursor, not
    /// the scan's selected/primary m/z). `None` for every ordinary PSM, in which
    /// case the PIN writer uses the spectrum's `precursor_mz`. `Some(mz)` makes
    /// the PIN writer compute ExpMass / dm / absdm from the secondary's own
    /// precursor mass instead of the primary's selected m/z.
    pub precursor_mz_override: Option<f64>,
}

impl PsmMatch {
    /// Returns the first (primary) candidate index. Callers that need to
    /// resolve back to a single Candidate use this; PIN writer iterates
    /// `candidate_idxs` directly to emit the multi-protein `Proteins` column.
    pub fn primary_candidate_idx(&self) -> u32 {
        self.candidate_idxs[0]
    }
}

impl PartialEq for PsmMatch {
    fn eq(&self, other: &Self) -> bool {
        // iter37 HIGH-2: PartialEq MUST agree with `Ord::cmp` (Rust contract
        // a == b ⇒ a.cmp(b) == Equal). Ord uses (spec_e_value, rank_score)
        // post-iter33, so PartialEq must compare the same fields. Pre-iter37
        // this compared `score` (= node + cleavage), violating the contract
        // for any pair of PSMs with equal `score` but different `rank_score`
        // (= `score + edge`). BinaryHeap behavior was technically undefined
        // for those pairs.
        self.spec_e_value == other.spec_e_value && self.rank_score == other.rank_score
    }
}

impl Eq for PsmMatch {}

impl PartialOrd for PsmMatch {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// Primary: `spec_e_value` ascending (lower = better).
/// Secondary: `rank_score` descending (higher = better).
///
/// iter33: `rank_score` is the Java-aligned queue-ordering key `node +
/// cleavage + edge`. Pre-iter33 the secondary key was just `score`
/// (= node + cleavage); post-iter33 it's `rank_score` (= node + cleavage +
/// edge) so the queue selects Java-equivalent top-1 PSMs even though the
/// PIN RawScore distribution (iter19) stays unchanged at `node + cleavage`.
///
/// For pre-iter33 callers / test fixtures that never set `rank_score`, the
/// default of 0.0 means an unset `rank_score` would lose to a set one. The
/// `match_engine` candidate loop always sets both `score` and `rank_score`;
/// fixtures that build PsmMatch manually should set `rank_score = score`
/// to preserve old behavior.
///
/// This ordering is used by `TopNQueue`'s min-heap (via `Reverse<PsmMatch>`):
/// the heap's "minimum" element is the one with the *largest* spec_e_value
/// (worst), so `push` evicts it when over capacity.
impl Ord for PsmMatch {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        use std::cmp::Ordering;
        // "Better" PSM = smaller spec_e_value, then larger rank_score.
        // NaN values are treated as worst (sort last / lose to finite).
        let self_sev  = if self.spec_e_value.is_nan()  { f64::INFINITY }      else { self.spec_e_value };
        let other_sev = if other.spec_e_value.is_nan() { f64::INFINITY }      else { other.spec_e_value };
        match other_sev.partial_cmp(&self_sev).unwrap_or(Ordering::Equal) {
            Ordering::Equal => {
                let self_rank  = if self.rank_score.is_nan()  { f32::NEG_INFINITY } else { self.rank_score };
                let other_rank = if other.rank_score.is_nan() { f32::NEG_INFINITY } else { other.rank_score };
                self_rank.partial_cmp(&other_rank).unwrap_or(Ordering::Equal)
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
    /// Java parity: `size < n OR score == worst → add`.
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
                    // R-1 (2026-05-18): Java parity keeps tied
                    // PSMs at capacity (and SpecE ties on the per-spectrum
                    // merge). Rust now matches.
                    // The queue may exceed `capacity` when ties exist —
                    // `capacity` becomes a *minimum* top-N, not a hard cap.
                    self.heap.push(Reverse(m));
                }
                std::cmp::Ordering::Less => {
                    // m is strictly worse than the worst retained PSM: drop.
                }
            }
        }
    }

    /// Add a PSM unconditionally, bypassing the capacity/eviction logic of `push`.
    /// Used for chimeric Pass-2 secondary peptides, which are legitimate extra
    /// emissions (a distinct co-isolated peptide on the same scan), NOT competitors
    /// for the primary's top-N slot.
    pub fn force_push(&mut self, m: PsmMatch) {
        self.heap.push(Reverse(m));
    }

    pub fn len(&self) -> usize { self.heap.len() }
    pub fn is_empty(&self) -> bool { self.heap.is_empty() }

    /// Return the `rank_score` of the queue's WORST retained PSM in O(1).
    ///
    /// The min-heap stores `Reverse<PsmMatch>` so `heap.peek()` returns the
    /// PSM with the LOWEST `Ord` value — the candidate that would be
    /// evicted first if a strictly better PSM arrived. Returns `None` if
    /// the queue is empty.
    ///
    /// iter34: used by the per-candidate two-stage gating in
    /// `match_engine.rs` — candidates whose `pin_score + max_edge_bonus`
    /// cannot exceed the worst retained `rank_score` skip the expensive
    /// `psm_edge_score` computation entirely.
    pub fn worst_rank_score(&self) -> Option<f32> {
        self.heap.peek().map(|std::cmp::Reverse(m)| m.rank_score)
    }

    /// Queue capacity (the top-N target). Used by callers that need to
    /// distinguish "queue has spare capacity, accept everything" from
    /// "queue at capacity, must beat worst".
    pub fn capacity(&self) -> u32 { self.capacity }

    /// Iterate over all PSMs in the queue (order not guaranteed).
    pub fn iter_psms(&self) -> impl Iterator<Item = &PsmMatch> {
        self.heap.iter().map(|Reverse(m)| m)
    }

    /// Drain all PSMs from the queue, returning them in an unordered Vec.
    /// Leaves the queue empty after the call. The returned Vec preserves no
    /// particular order — callers that need ordering should sort the result.
    ///
    /// Cost: O(N) drain + Vec collection. Cheap for small N (top-N typically ≤ 10).
    pub fn drain_into_vec(&mut self) -> Vec<PsmMatch> {
        self.heap.drain().map(|Reverse(m)| m).collect()
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
        // Test-only PSM: candidate_idxs[0] = 0 is a sentinel for queue-ordering tests
        // that never resolve back to a real Candidate. Tests that need to read
        // peptide / protein metadata must build their own &[Candidate] alongside.
        PsmMatch {
            spectrum_idx,
            candidate_idxs: vec![0],
            charge_used: 2,
            mass_error_ppm: 0.0,
            score,
            rank_score: score,  // iter33 fixture default: rank_score = score
            edge_score: 0,
            spec_e_value: 1.0,  // default sentinel: "not yet computed"
            de_novo_score: i32::MIN,  // sentinel: not yet computed
            activation_method: None,
            e_value: 1.0,  // sentinel: not yet computed
            features: PsmFeatures::default(),
            isotope_offset: 0,
            precursor_mz_override: None,
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
        // R-1 fix: Java parity keeps tied PSMs at capacity (raw-score
        // retention and SpecE merge). Rust's TopNQueue must mirror this —
        // strict-greater eviction
        // was dropping ties Java keeps, plausibly causing the Astral 14K raw-
        // target gap that R-1 + R-2 closed.
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
    fn force_push_bypasses_capacity_eviction() {
        // Chimeric Pass-2: a capacity-1 queue holds the primary (Pass 1, top-1).
        // A Pass-2 secondary is a distinct co-isolated peptide, NOT a competitor
        // for the primary's slot — force_push must ADD it without eviction.
        let mut q = TopNQueue::new(1);
        // Primary: best score. Secondary: strictly worse — `push` would drop it.
        q.push(make_match(0, 100.0));
        let mut secondary = make_match(0, 50.0);
        secondary.candidate_idxs = vec![1]; // distinct peptide
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
        // R-2.2 (2026-05-18): synthetic test for pepSeq+score dedup. Two PSMs
        // with the same (peptide_residue, score) key should collapse to one
        // PsmMatch with both candidate_idxs aggregated into the surviving Vec.
        //
        // We use drain_into_vec to extract PSMs, then assert the dedup helper
        // collapses them correctly.

        let mut q = TopNQueue::new(10);
        // Three PSMs: two share (peptide=0, score=50), one is distinct (peptide=1, score=40)
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

        // Caller (match_engine) provides the key function. Here we use
        // a synthetic key based on score only (test scaffolding — real
        // dedup uses peptide_residue + rounded_score from candidates).
        let deduped = simple_dedup_by_score_for_test(drained);

        // Expect: 2 groups — score=50 with idxs [10,20], score=40 with [30]
        assert_eq!(deduped.len(), 2, "should collapse to 2 unique-score groups");

        let mut score_50 = deduped.iter().find(|p| (p.score as i32) == 50).unwrap().candidate_idxs.clone();
        score_50.sort();
        assert_eq!(score_50, vec![10, 20], "score=50 should aggregate both idxs");

        let score_40 = &deduped.iter().find(|p| (p.score as i32) == 40).unwrap().candidate_idxs;
        assert_eq!(*score_40, vec![30]);
    }

    /// Test-only dedup that groups by score alone (real production
    /// dedup_pepseq_score in match_engine.rs uses peptide_residue + score).
    fn simple_dedup_by_score_for_test(psms: Vec<PsmMatch>) -> Vec<PsmMatch> {
        use std::collections::HashMap;
        let mut groups: HashMap<i32, PsmMatch> = HashMap::new();
        for psm in psms {
            let key = psm.score as i32;
            groups
                .entry(key)
                .and_modify(|existing| existing.candidate_idxs.extend(psm.candidate_idxs.iter().copied()))
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
