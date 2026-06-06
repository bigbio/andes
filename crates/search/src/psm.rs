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
    /// learn weights without disrupting the existing RawScore distribution.
    /// Kept separate from `RawScore` rather than blended into it.
    /// Computed via `psm_edge_score` in `score_psm.rs`.
    pub edge_score: i32,

    // ── Chimeric MS1 precursor-envelope features ───────────────────────────
    /// KL-divergence between the observed precursor isotope envelope (from the
    /// linked MS1) and the averagine theoretical envelope. Higher = poorer match
    /// = likely spurious co-isolation. 0.0 when MS1/feature unavailable.
    pub precursor_isotope_kl: f32,
    /// Observed monoisotopic precursor intensity / median MS1 intensity (SNR
    /// proxy). 0.0 when unavailable.
    pub precursor_snr: f32,

    /// Top-1 RawScore dominance: `RawScore(best) − RawScore(2nd-best distinct
    /// peptide)` for this spectrum, a positive "lead over the runner-up"
    /// confidence signal. Emitted as the additive `DeltaRawScore` PIN column on
    /// the rank-1 row only (0.0 elsewhere). Captured during candidate scoring
    /// (NOT from the retained queue), so it is populated even at `top_n = 1`
    /// where the runner-up is otherwise evicted — and without perturbing the
    /// GF `min_score` / SpecEValue of any emitted PSM (purely additive).
    pub delta_raw_score: f32,

    /// Tailor per-spectrum score calibration (Yang et al., JPR 2020):
    /// `RawScore / denom`, where `denom` is this spectrum's RawScore at the
    /// top-1% quantile of its candidate-score distribution (the score at rank
    /// `ceil(0.01 * N)` from the top). Dividing each PSM's RawScore by its own
    /// spectrum's high-quantile score makes RawScores comparable across spectra
    /// — the role the removed generating function used to play. Emitted as the
    /// additive `TailorScore` PIN column. Falls back to `denom = 1.0`
    /// (TailorScore = RawScore) when the spectrum has too few candidates
    /// (< 100) or the quantile score is `<= 0`, so it never divides by zero or
    /// amplifies noise on tiny candidate sets. Purely additive: the histogram
    /// it is computed from never alters which candidates are scored or kept.
    pub tailor_score: f32,

    // ── Strong-score Stage-1 additive bolt-ons (new PIN columns) ───────────
    /// `Σ exp(-½ (ppmᵢ/σ)²)` over all matched b/y ions (σ = 7 ppm). A
    /// Gaussian-kernel "tight-match evidence" sum: a fragment matched at a few
    /// ppm contributes ~1, one scattered near the tolerance edge contributes
    /// ~0. Turns fragment mass accuracy into evidence (the rank model discards
    /// it), so a real high-res PSM whose ions cluster at low ppm scores far
    /// above a coincidental match whose ions scatter. ~0 for all low-res PSMs
    /// (every ppm is huge), so Percolator simply down-weights it there.
    pub ppm_gaussian_score: f32,
    /// Number of backbone cleavage sites where BOTH the b-ion and its
    /// complementary y-ion (bond `i` ↔ `b_i` + `y_{n-i}`) are observed.
    /// Both halves of one cleavage matching by chance is rare, so this is a
    /// strong "real peptide" signal that `num_matched_main_ions` (b OR y per
    /// position) does not capture.
    pub complementary_ion_count: u32,
    /// Longest consecutive run of complementary cleavage sites (both `b_i`
    /// and `y_{n-i}` matched). A contiguous ladder is stronger evidence than
    /// scattered complementary pairs.
    pub longest_complementary_ladder: u32,
    /// Fraction of the spectrum's top-20 most-intense peaks that are NOT
    /// matched by any predicted b/y ion (within feature tolerance). Lower =
    /// better — a real PSM explains most big peaks; a coincidental match
    /// leaves intense unexplained signal.
    pub unexplained_top_intensity_fraction: f32,
    /// Count of matched b/y ions that also have a neutral-loss partner peak at
    /// −H2O (−18.0106 Da) or −NH3 (−17.0265 Da) within feature tolerance.
    /// Strong CID/TMT signal that the intensity-rank model underuses.
    pub neutral_loss_ion_count: u32,
    /// Mean intensity-rank of matched b/y ions (rank 1 = most intense peak).
    /// Real PSMs match dominant peaks; coincidental matches hit weak peaks.
    pub mean_matched_intensity_rank: f32,

    // ── Strong-score Stage-2: competition/null denominator (the moat) ──────
    /// `Σ max(0, -ln(ρᵢ·Δᵢ))` over matched ions — the per-peak Poisson
    /// chance-match "surprise". `ρᵢ` = local peak density (peaks/Da) around the
    /// matched peak, `Δᵢ` = the match window (Da). A match in a sparse region at
    /// tight tolerance is improbable by chance (large surprise = strong
    /// evidence); a match in a crowded region within a wide window is nearly
    /// free (≈0). This is the deterministic null no ML rescorer computes — it
    /// weighs each match by how *non-coincidental* it is.
    pub chance_match_surprise: f32,
}

/// Number of candidates below which Tailor calibration is skipped (denom = 1.0).
/// On tiny candidate sets the top-1% quantile is noisy / undefined, so we fall
/// back to the uncalibrated RawScore rather than amplify noise.
pub const TAILOR_MIN_CANDIDATES: u32 = 100;

/// Tailor top-quantile fraction (q): the denominator is the RawScore at the
/// top `q` of the candidate-score distribution (Yang et al. use q = 0.01).
pub const TAILOR_QUANTILE: f64 = 0.01;

/// Compute the Tailor calibration denominator for one spectrum from a histogram
/// of its candidate RawScores.
///
/// `hist` maps a (rounded, integer) RawScore to the number of candidates that
/// achieved it; `total` is the total candidate count (== sum of histogram
/// counts). Returns the RawScore `s` at the top-`TAILOR_QUANTILE` of the
/// distribution: the smallest `s` such that
/// `count(candidates with RawScore >= s) <= ceil(q * total)` — i.e. the score
/// at rank `ceil(q * total)` counting from the highest score down.
///
/// Returns the fallback `1.0` (no calibration) when `total < TAILOR_MIN_CANDIDATES`
/// or when the resolved quantile score is `<= 0`, so `RawScore / denom` never
/// divides by zero and never amplifies noise on tiny candidate sets.
pub fn tailor_denominator<S: std::hash::BuildHasher>(
    hist: &std::collections::HashMap<i32, u32, S>,
    total: u32,
) -> f64 {
    if total < TAILOR_MIN_CANDIDATES {
        return 1.0;
    }
    // Number of candidates that define the "top" tail. ceil(q * total), at least 1.
    let top_k = (TAILOR_QUANTILE * total as f64).ceil() as u64;
    let top_k = top_k.max(1);

    // Walk scores from highest to lowest, accumulating counts. The quantile
    // score is the lowest score still inside the top-`top_k` candidates: the
    // first score at which the running count (candidates with RawScore >= s)
    // reaches/exceeds `top_k`.
    let mut scores: Vec<i32> = hist.keys().copied().collect();
    scores.sort_unstable_by(|a, b| b.cmp(a)); // descending
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
        _ => 1.0, // <= 0 quantile score → fall back (avoids div-by-zero / sign flip)
    }
}

#[derive(Debug, Clone)]
pub struct PsmMatch {
    pub spectrum_idx: usize,
    /// Indices into the `&[Candidate]` slice owned by `PreparedSearch.candidates`.
    /// Length is always ≥ 1. The first index (`candidate_idxs[0]`) is the
    /// "primary" candidate — used by callers that need a single Candidate
    /// (most do; see `primary_candidate_idx()`). Multiple indices accumulate
    /// when the pepSeq+score dedup pass merges multiple Candidates that
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
    /// This is what gets emitted in the `RawScore` PIN column. Used by
    /// Percolator as one of many features.
    pub score: f32,
    /// Queue-ordering score = `node + cleavage + edge`. Java's
    /// `DBScanScorer.getScore` returns `node + edge` and Java parity adds
    /// cleavage, so Java's `match.score` (used by its `PriorityQueue`
    /// ordering) is `node + cleavage + edge`. Rust's pin RawScore stays at
    /// `node + cleavage` for Percolator distribution stability; the
    /// SEPARATE `EdgeScore` PIN column carries the `+edge` contribution.
    /// `rank_score` mirrors Java's queue-ordering key without changing the
    /// pin RawScore distribution.
    ///
    /// **No automatic default**: PsmMatch does not implement `Default`, and
    /// callers MUST set `rank_score` explicitly. Test fixtures that build
    /// PsmMatch literals should set `rank_score = score` for behavior with
    /// no edge contribution to ranking. The `match_engine.rs`
    /// candidate loop computes `rank_score = score + edge_score as f32`.
    pub rank_score: f32,
    /// Per-PSM edge_score = `psm_edge_score(...)` for this candidate.
    /// Computed at queue-insertion time in `match_engine.rs` and reused by
    /// `compute_psm_features` to populate the `EdgeScore` PIN column
    /// (avoids the recompute). Default 0 — features extraction will compute
    /// it on the fly if it remains 0 (e.g. for test fixtures).
    pub edge_score: i32,
    /// Activation method captured from `param.data_type.activation` at scoring
    /// time. `None` if unknown or not yet set.
    pub activation_method: Option<model::activation::ActivationMethod>,
    /// Fragment-ion feature columns computed after `score_psm`.
    /// Defaults to all-zero until `compute_psm_features` runs.
    pub features: PsmFeatures,
    /// The isotope offset that produced the precursor match: 0 = monoisotopic,
    /// +N = spectrum precursor was N C13 peaks above the true monoisotopic.
    /// Default range −1..=2. Threaded from `MassError::isotope_offset`
    /// (precursor_matching.rs) via match_engine.rs. Written as the PIN
    /// `isotope_error` column.
    pub isotope_offset: i8,
    /// Chimeric Pass-2 secondaries only: the co-isolated precursor m/z this PSM was
    /// scored at. `Some(mz)` makes the PIN/TSV writers compute ExpMass / dm / absdm
    /// from the secondary's own precursor mass; `None` (every ordinary PSM) falls
    /// back to the spectrum's `precursor_mz`.
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
        // PartialEq MUST agree with `Ord::cmp` (Rust contract
        // a == b ⇒ a.cmp(b) == Equal). Ord ranks by `rank_score` alone,
        // so PartialEq compares the same field. Comparing `score`
        // (= node + cleavage) would violate the contract for any pair of
        // PSMs with equal `score` but different `rank_score`
        // (= `score + edge`), leaving BinaryHeap behavior undefined for
        // those pairs. Delegate to `cmp` (rather than a raw `==`) so the two
        // agree on NaN: `cmp` maps NaN `rank_score` to NEG_INFINITY, making
        // two NaN-ranked PSMs compare Equal — a raw float `==` would report
        // them unequal and break the Eq/Ord contract.
        self.cmp(other) == std::cmp::Ordering::Equal
    }
}

impl Eq for PsmMatch {}

impl PartialOrd for PsmMatch {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// `rank_score` descending (higher = better).
///
/// `rank_score` is the Java-aligned queue-ordering key `node +
/// cleavage + edge`, so the queue selects Java-equivalent top-1 PSMs even
/// though the PIN RawScore distribution stays unchanged at `node + cleavage`.
/// With the generating function removed, `rank_score` is the sole ranking
/// signal for both queue retention and output ordering.
///
/// For callers / test fixtures that never set `rank_score`, the
/// default of 0.0 means an unset `rank_score` would lose to a set one. The
/// `match_engine` candidate loop always sets both `score` and `rank_score`;
/// fixtures that build PsmMatch manually should set `rank_score = score`
/// to preserve old behavior.
///
/// This ordering is used by `TopNQueue`'s min-heap (via `Reverse<PsmMatch>`):
/// the heap's "minimum" element is the one with the *smallest* rank_score
/// (worst), so `push` evicts it when over capacity.
impl Ord for PsmMatch {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        use std::cmp::Ordering;
        // "Better" PSM = larger rank_score.
        // NaN values are treated as worst (sort last / lose to finite).
        let self_rank  = if self.rank_score.is_nan()  { f32::NEG_INFINITY } else { self.rank_score };
        let other_rank = if other.rank_score.is_nan() { f32::NEG_INFINITY } else { other.rank_score };
        self_rank.partial_cmp(&other_rank).unwrap_or(Ordering::Equal)
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
    /// "Best" = largest `rank_score` first. The min-heap (via
    /// `Reverse<PsmMatch>`) puts the *worst* PSM at the top so it can be
    /// evicted when a strictly-better PSM arrives.
    ///
    /// **Tie handling:** when the queue is at capacity and
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
                    // Java parity keeps tied PSMs at capacity (and SpecE
                    // ties on the per-spectrum merge).
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
    /// Used by the per-candidate two-stage gating in
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
    pub fn fill_post_topn<F: FnMut(&mut PsmMatch)>(&mut self, mut f: F) {
        let mut psms: Vec<PsmMatch> = self.heap.drain().map(|Reverse(m)| m).collect();
        for psm in &mut psms {
            f(psm);
        }
        for psm in psms {
            self.heap.push(Reverse(psm));
        }
    }

    /// Return the best PSM (largest `rank_score`) without removing it.
    /// Returns `None` if the queue is empty.
    ///
    /// The heap is a min-heap on `Reverse<PsmMatch>` so the *worst* entry sits
    /// at the top (for cheap eviction). To find the *best* entry we iterate
    /// all elements and take the max in natural `PsmMatch` ordering.
    /// Cost is O(N) — acceptable for the small top-N queues used in practice.
    pub fn peek_top(&self) -> Option<&PsmMatch> {
        self.heap.iter().map(|Reverse(m)| m).max_by(|a, b| a.cmp(b))
    }

    /// Drain into a Vec sorted best-first (largest `rank_score`).
    pub fn into_sorted_vec(self) -> Vec<PsmMatch> {
        let mut v: Vec<PsmMatch> = self.heap.into_iter().map(|Reverse(m)| m).collect();
        v.sort_by(|a, b| b.cmp(a));
        v
    }

    /// Drain into a Vec sorted best-first by `rank_score` DESCENDING (the
    /// larger RawScore = better), then by `score` descending as a stable
    /// tiebreak. This is the canonical output ordering: `rank_score` is the
    /// sole ranking signal now that the generating function is removed.
    pub fn into_rank_sorted_vec(self) -> Vec<PsmMatch> {
        let mut v: Vec<PsmMatch> = self.heap.into_iter().map(|Reverse(m)| m).collect();
        v.sort_by(|a, b| {
            let ar = if a.rank_score.is_nan() { f32::NEG_INFINITY } else { a.rank_score };
            let br = if b.rank_score.is_nan() { f32::NEG_INFINITY } else { b.rank_score };
            br.partial_cmp(&ar)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| {
                    let asc = if a.score.is_nan() { f32::NEG_INFINITY } else { a.score };
                    let bsc = if b.score.is_nan() { f32::NEG_INFINITY } else { b.score };
                    bsc.partial_cmp(&asc).unwrap_or(std::cmp::Ordering::Equal)
                })
        });
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
            rank_score: score,  // fixture default: rank_score = score
            edge_score: 0,
            activation_method: None,
            features: PsmFeatures::default(),
            isotope_offset: 0,
            precursor_mz_override: None,
        }
    }

    /// Build a fixture PSM whose `rank_score` is set independently of `score`.
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
        for s in [1.0, 2.0, 3.0] { q.push(make_match(0, s)); }
        assert_eq!(q.len(), 3);
        let sorted = q.into_sorted_vec();
        // rank_score = score (fixture default) → sort by score descending.
        assert_eq!(sorted.iter().map(|m| m.score).collect::<Vec<_>>(),
                   vec![3.0, 2.0, 1.0]);
    }

    #[test]
    fn queue_at_capacity_keeps_top_n_by_score() {
        let mut q = TopNQueue::new(3);
        for s in [1.0, 5.0, 2.0, 4.0, 3.0] { q.push(make_match(0, s)); }
        assert_eq!(q.len(), 3);
        let sorted = q.into_sorted_vec();
        // rank_score = score → keeps top-3 by score.
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
        // Java parity keeps tied PSMs at capacity (raw-score retention and
        // SpecE merge). Rust's TopNQueue must mirror this — strict-greater
        // eviction would drop ties Java keeps.
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
        // Synthetic test for pepSeq+score dedup. Two PSMs
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
        assert_eq!(cloned.rank_score, 4.2);
    }

    // -----------------------------------------------------------------------
    // rank_score ordering tests
    // -----------------------------------------------------------------------

    #[test]
    fn psm_match_orders_by_rank_score_descending() {
        // Higher rank_score means "better" → Ord-greater so the min-heap can
        // evict the worst.
        let better = make_match_with_rank(0, 5.0, 10.0);
        let worse  = make_match_with_rank(0, 5.0, 1.0);
        assert!(better > worse,
            "PSM with higher rank_score should be Ord-greater (better in the min-heap)");
    }

    #[test]
    fn queue_keeps_best_rank_score_psms_when_full() {
        // Three PSMs with same score but different rank_scores; capacity = 2.
        let mut q = TopNQueue::new(2);
        q.push(make_match_with_rank(0, 5.0, 1.0));   // worst
        q.push(make_match_with_rank(0, 5.0, 100.0)); // best
        assert_eq!(q.len(), 2);
        // Push a medium one; it should evict the worst (rank 1.0).
        q.push(make_match_with_rank(0, 5.0, 50.0));
        assert_eq!(q.len(), 2);
        let sorted = q.into_sorted_vec();
        let ranks: Vec<f32> = sorted.iter().map(|m| m.rank_score).collect();
        assert!(ranks.contains(&100.0), "best rank 100.0 should be retained");
        assert!(ranks.contains(&50.0),  "medium rank 50.0 should be retained");
        assert!(!ranks.contains(&1.0),  "worst rank 1.0 should be evicted");
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
    fn psm_match_with_nan_rank_score_orders_as_worst() {
        // NaN rank_score should sort as WORSE than any finite value.
        // "Better" = greater in natural Ord (used by the min-heap via Reverse).
        let nan_rank = make_match_with_rank(0, 5.0, f32::NAN);
        let finite   = make_match_with_rank(0, 0.0, 1.0);
        assert_eq!(
            nan_rank.cmp(&finite),
            std::cmp::Ordering::Less,
            "NaN rank_score should sort as worse (Less) than a finite value"
        );
    }

    // -----------------------------------------------------------------------
    // Tailor per-spectrum calibration
    // -----------------------------------------------------------------------

    use std::collections::HashMap;

    /// Build a histogram with `score` repeated `count` times, summing total.
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
        // 1000 candidates: 990 at score 5, 10 at score 50 (the top 1%).
        // top_k = ceil(0.01 * 1000) = 10. Walking from highest: 10 candidates
        // at score 50 reach top_k exactly → quantile score = 50.
        let (h, total) = hist_from(&[(5, 990), (50, 10)]);
        assert_eq!(total, 1000);
        assert_eq!(tailor_denominator(&h, total), 50.0);
        // TailorScore for a RawScore-50 PSM = 50 / 50 = 1.0.
        assert!(((50.0f64 / tailor_denominator(&h, total)) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn tailor_denom_ceil_rounds_up() {
        // 150 candidates → ceil(0.01 * 150) = ceil(1.5) = 2 candidates in the
        // top tail. 1 candidate at 100, 1 at 90, rest at 10. Walking down:
        // 100 (cum 1) < 2, then 90 (cum 2) >= 2 → quantile score = 90.
        let (h, total) = hist_from(&[(100, 1), (90, 1), (10, 148)]);
        assert_eq!(total, 150);
        assert_eq!(tailor_denominator(&h, total), 90.0);
    }

    #[test]
    fn tailor_denom_small_n_falls_back_to_one() {
        // < TAILOR_MIN_CANDIDATES (100): no calibration → denom 1.0, so
        // TailorScore == RawScore.
        let (h, total) = hist_from(&[(5, 50), (50, 5)]);
        assert_eq!(total, 55);
        assert!(total < TAILOR_MIN_CANDIDATES);
        assert_eq!(tailor_denominator(&h, total), 1.0);
    }

    #[test]
    fn tailor_denom_nonpositive_quantile_falls_back_to_one() {
        // Enough candidates but the top-1% quantile score is <= 0 (all negative
        // RawScores). Must fall back to 1.0 to avoid div-by-zero / sign flip.
        let (h, total) = hist_from(&[(-10, 990), (-1, 10)]);
        assert_eq!(total, 1000);
        // top tail is the 10 candidates at -1, which is <= 0 → fallback.
        assert_eq!(tailor_denominator(&h, total), 1.0);
    }

    #[test]
    fn tailor_denom_exactly_min_candidates_calibrates() {
        // Exactly TAILOR_MIN_CANDIDATES → calibrates (>= threshold).
        // 100 candidates, top_k = ceil(0.01*100) = 1. Highest score = 40.
        let (h, total) = hist_from(&[(10, 99), (40, 1)]);
        assert_eq!(total, 100);
        assert_eq!(tailor_denominator(&h, total), 40.0);
    }

    #[test]
    fn psm_match_two_nan_rank_scores_compare_equal() {
        // Two PSMs both with NaN rank_score → Equal.
        let a = make_match_with_rank(0, 5.0, f32::NAN);
        let b = make_match_with_rank(0, 5.0, f32::NAN);
        assert_eq!(
            a.cmp(&b),
            std::cmp::Ordering::Equal,
            "Two PSMs with NaN rank_score should compare Equal"
        );
    }
}
