//! Convert accumulated [`CountStats`] into a scoring [`Param`] with Laplace
//! smoothing and partition backoff.
//!
//! # `rank_dist_table` representation
//!
//! [`RankScorer::new`] expects **raw frequency (probability) vectors**, NOT
//! log-scores.  Each vector has length `max_rank + 1`:
//!
//! - Indices `0 .. max_rank-1` → observed ranks 1 .. max_rank.
//! - Index `max_rank` → the "missing ion" sentinel slot.
//!
//! Both `IonType::Noise` and each non-noise ion get an entry of this length
//! for every partition.  `RankScorer::new` computes
//!
//! ```text
//! log_score[i] = ln(ion_freq[i] / (noise_freq[i] * charge_or_seg))
//! ```
//!
//! so the frequencies must be strictly positive everywhere.  Laplace smoothing
//! guarantees that.
//!
//! # Partition backoff hierarchy
//!
//! When a partition's total rank-count `n < cfg.min_count` we blend the
//! empirical distribution with a coarser "parent":
//!
//! 1. **Segment collapse** – collapse the `seg_num` dimension: for the same
//!    `(charge, parent_mass.to_bits())`, sum counts across all segments.
//! 2. **Global pool** – sum all counts across all partitions for the ion type.
//!
//! The blended value is `(n*emp + w*parent) / (n+w)` where `w = backoff_weight`.
//!
//! # High-res vs low-res `error_scaling_factor`
//!
//! The template's `error_scaling_factor` is used verbatim unless overridden via
//! `EstimatorConfig::error_scaling_factor_override`.  The template already
//! encodes the correct ESF for its instrument class (high-res instruments such
//! as QExactive / OrbitrapAstral / TimsTOF use ESF ≈ 40–100; LowRes uses a
//! smaller value or 0).

use rustc_hash::FxHashMap;

use scoring_crate::param_model::{IonType, Param, Partition};

use crate::counts::CountStats;

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

/// Hyper-parameters for the [`Estimator`].
#[derive(Debug, Clone)]
pub struct EstimatorConfig {
    /// Laplace pseudo-count added to **every** signal (prefix/suffix) rank/error
    /// bin before normalising.  Must be > 0.  Default: `1.0`.
    pub pseudo: f32,
    /// Laplace pseudo-count for the **noise** rank distribution.  Noise is
    /// abundantly sampled and strongly concentrated (random theoretical
    /// positions are almost always *absent*), so the signal `pseudo=1.0` spread
    /// over ~`max_rank` slots over-flattens it — inflating `noise_freq[r]` at
    /// the ranks where real ions land and compressing `ln(ion/noise)` node
    /// scores (the diagnosed estimator dilution). A much smaller value keeps the
    /// noise model sharp while still flooring every slot > 0 for log-safety.
    /// Must be > 0.  Default: `0.05`.
    pub noise_pseudo: f32,
    /// Partition total below which backoff blending is applied.  Default: `50`.
    pub min_count: u64,
    /// Prior weight `w` in `(n*emp + w*parent) / (n+w)`.  Default: `20.0`.
    pub backoff_weight: f32,
    /// If `Some(esf)`, override the template's `error_scaling_factor`; if
    /// `None`, copy from the template.  Default: `None`.
    pub error_scaling_factor_override: Option<i32>,
}

impl Default for EstimatorConfig {
    fn default() -> Self {
        Self {
            pseudo: 1.0,
            noise_pseudo: 0.05,
            min_count: 50,
            backoff_weight: 20.0,
            error_scaling_factor_override: None,
        }
    }
}

// ---------------------------------------------------------------------------
// Estimator
// ---------------------------------------------------------------------------

/// Turns [`CountStats`] into a scoring [`Param`] with Laplace smoothing and
/// partition backoff.
pub struct Estimator {
    cfg: EstimatorConfig,
}

impl Estimator {
    /// Create a new `Estimator` with the given configuration.
    pub fn new(cfg: EstimatorConfig) -> Self {
        Self { cfg }
    }

    /// Build a [`Param`] from accumulated counts.
    ///
    /// `template` supplies scalar metadata that is **not** learned from counts:
    /// `data_type` / activation / instrument / enzyme / protocol, `mme`
    /// tolerance, `apply_deconvolution` / `deconvolution_error_tolerance`,
    /// `num_segments`, and the ion-type / frag-offset layout.
    ///
    /// The learned tables — `rank_dist_table`, `ion_err_dist_table`,
    /// `noise_err_dist_table`, `ion_existence_table`, `charge_hist`,
    /// `partitions` — are built from `counts`.
    pub fn estimate(&self, counts: &CountStats, template: &Param) -> Param {
        self.estimate_with_prior(counts, template, None)
    }

    /// Like [`estimate`], but also accepts an optional independent prior model.
    ///
    /// When `prior` is `Some(p)`, sparse partitions (total rank-count below
    /// `cfg.min_count`) shrink toward the prior's rank distribution rather than
    /// (only) the corpus-internal pool.  The prior's distribution for a given
    /// `(partition, ion)` is used as **Level 0** of the backoff hierarchy — it
    /// is consulted before segment-collapse (Level 1) and the global pool
    /// (Level 2).  Passing `None` is identical to calling [`estimate`].
    ///
    /// **Invariant:** every distribution vector in `prior.rank_dist_table` must
    /// be strictly positive in every slot (the same requirement
    /// [`RankScorer::new`] places on the trained model). A prior produced by
    /// [`estimate`] always satisfies this; a hand-crafted prior with zero
    /// entries can yield `-inf`/`NaN` node scores when a partition has zero
    /// observed counts (the blend collapses to the prior verbatim).
    pub fn estimate_with_prior(
        &self,
        counts: &CountStats,
        template: &Param,
        prior: Option<&Param>,
    ) -> Param {
        // `pseudo` must be > 0: it is the Laplace smoothing mass that guarantees
        // every distribution slot is strictly positive. With `pseudo == 0` an
        // all-zero raw count would normalize to all-zero probabilities, and
        // `RankScorer::new` would then compute `ln(0) = -inf` and break scoring.
        assert!(
            self.cfg.pseudo > 0.0,
            "EstimatorConfig.pseudo must be > 0 (got {})",
            self.cfg.pseudo
        );
        assert!(
            self.cfg.noise_pseudo > 0.0,
            "EstimatorConfig.noise_pseudo must be > 0 (got {})",
            self.cfg.noise_pseudo
        );
        let max_rank = template.max_rank;
        let esf = self.cfg.error_scaling_factor_override
            .unwrap_or(template.error_scaling_factor);

        let rank_dist_table = self.build_rank_dist_table(counts, template, max_rank, prior);
        let (ion_err_dist_table, noise_err_dist_table) =
            self.build_error_tables(counts, template, esf, prior);
        let ion_existence_table = self.build_existence_table(counts, template, prior);
        let charge_hist = build_charge_hist(&counts.charge);
        let (min_charge, max_charge) = charge_range(&charge_hist, template);

        // Sorted partitions (same invariant as the binary loader).
        let mut partitions: Vec<Partition> = rank_dist_table.keys().copied().collect();
        partitions.sort();

        let mut param = Param {
            version: template.version,
            data_type: template.data_type.clone(),
            mme: template.mme,
            apply_deconvolution: template.apply_deconvolution,
            deconvolution_error_tolerance: template.deconvolution_error_tolerance,
            charge_hist,
            min_charge,
            max_charge,
            num_segments: template.num_segments,
            partitions,
            num_precursor_off: template.num_precursor_off,
            precursor_off_map: template.precursor_off_map.clone(),
            frag_off_table: template.frag_off_table.clone(),
            max_rank,
            rank_dist_table,
            error_scaling_factor: esf,
            ion_err_dist_table,
            noise_err_dist_table,
            ion_existence_table,
            partition_ion_types_cache: FxHashMap::default(),
        };
        // Rebuild the per-partition ion-type cache required by RankScorer::new
        // and ion_types_for_partition_slice.
        param.rebuild_cache();
        param
    }

    // -----------------------------------------------------------------------
    // rank_dist_table
    // -----------------------------------------------------------------------

    fn build_rank_dist_table(
        &self,
        counts: &CountStats,
        template: &Param,
        max_rank: i32,
        prior: Option<&Param>,
    ) -> FxHashMap<Partition, FxHashMap<IonType, Vec<f32>>> {
        // Array length: max_rank observed-rank slots + 1 missing-ion slot.
        let n_slots = (max_rank + 1) as usize;
        let pseudo = self.cfg.pseudo;
        let noise_pseudo = self.cfg.noise_pseudo;
        let min_count = self.cfg.min_count;
        let w = self.cfg.backoff_weight;

        // Universe of partitions: union of template's frag_off_table and
        // partitions present in counts.
        let mut all_partitions: std::collections::HashSet<Partition> =
            template.frag_off_table.keys().copied().collect();
        for &(part, _) in counts.rank.keys() {
            all_partitions.insert(part);
        }

        // Per-partition ion lists (from the template's fragment-offset layout).
        let ion_lists: FxHashMap<Partition, Vec<IonType>> = template
            .frag_off_table
            .iter()
            .map(|(&part, frags)| {
                let ions: Vec<IonType> = frags.iter().map(|f| f.ion_type).collect();
                (part, ions)
            })
            .collect();

        // Global pool: per-ion, sum across all partitions.
        let global_pool = build_global_pool(&counts.rank, n_slots);

        // Segment-collapsed pool: per-(charge, parent_mass.to_bits()), sum
        // across all seg_num values.
        let seg_collapsed = build_seg_collapsed(&counts.rank, n_slots);

        let mut out: FxHashMap<Partition, FxHashMap<IonType, Vec<f32>>> =
            FxHashMap::default();

        for &part in &all_partitions {
            // Only emit partitions that have an ion list in the template;
            // RankScorer requires a Noise entry for every populated partition.
            let ions = match ion_lists.get(&part) {
                Some(v) if !v.is_empty() => v.clone(),
                _ => continue,
            };

            let mut ion_table: FxHashMap<IonType, Vec<f32>> = FxHashMap::default();

            // Precompute the segment-collapsed parent map for this partition's
            // (charge, parent_mass.to_bits()) key (used by all ions in this partition).
            let seg_key = (part.charge, part.parent_mass.to_bits());
            let seg_parent: Option<&FxHashMap<IonType, Vec<u64>>> =
                seg_collapsed.get(&seg_key);

            // Helper: compute a normalised parent vector for `ion` using the
            // given pseudo-count (signal vs noise differ — see `noise_pseudo`).
            let parent_vec = |ion: IonType, ps: f32| -> Vec<f32> {
                // Level 0: independent prior model (own-data broad prior).
                if let Some(pr) = prior {
                    if let Some(dist) = pr.rank_dist_table.get(&part).and_then(|m| m.get(&ion)) {
                        if dist.len() == n_slots {
                            // Returned as-is: the prior is already normalized;
                            // `ps` only applies to Levels 1/2 where raw counts
                            // still need Laplace smoothing.
                            return dist.clone();
                        }
                    }
                }
                // Level 1: segment-collapse.
                if let Some(seg_map) = seg_parent {
                    if let Some(raw) = seg_map.get(&ion) {
                        let n: u64 = raw.iter().sum();
                        if n >= min_count {
                            return normalize_with_pseudo(raw, n_slots, ps);
                        }
                    }
                }
                // Level 2: global pool.
                let graw = global_pool.get(&ion)
                    .map(|v| v.as_slice())
                    .unwrap_or(&[]);
                normalize_with_pseudo(graw, n_slots, ps)
            };

            for &ion in &ions {
                let raw = counts.rank.get(&(part, ion))
                    .map(|v| v.as_slice())
                    .unwrap_or(&[]);
                let n: u64 = raw.iter().sum();
                let emp = normalize_with_pseudo(raw, n_slots, pseudo);
                let blended = if n < min_count {
                    blend(&emp, &parent_vec(ion, pseudo), n as f32, w)
                } else {
                    emp
                };
                ion_table.insert(ion, blended);
            }

            // Noise is required by RankScorer::new. It uses a much smaller
            // pseudo-count (`noise_pseudo`) so its sharp, concentrated shape is
            // preserved rather than flattened across all `max_rank` slots.
            let noise_raw = counts.rank.get(&(part, IonType::Noise))
                .map(|v| v.as_slice())
                .unwrap_or(&[]);
            let noise_n: u64 = noise_raw.iter().sum();
            let noise_emp = normalize_with_pseudo(noise_raw, n_slots, noise_pseudo);
            let noise_dist = if noise_n < min_count {
                blend(&noise_emp, &parent_vec(IonType::Noise, noise_pseudo), noise_n as f32, w)
            } else {
                noise_emp
            };
            ion_table.insert(IonType::Noise, noise_dist);

            out.insert(part, ion_table);
        }

        out
    }

    // -----------------------------------------------------------------------
    // error tables
    // -----------------------------------------------------------------------

    fn build_error_tables(
        &self,
        counts: &CountStats,
        template: &Param,
        esf: i32,
        prior: Option<&Param>,
    ) -> (FxHashMap<Partition, Vec<f32>>, FxHashMap<Partition, Vec<f32>>) {
        if esf <= 0 {
            return (FxHashMap::default(), FxHashMap::default());
        }
        let dist_len = (esf as usize) * 2 + 1;
        let pseudo = self.cfg.pseudo;
        let min_count = self.cfg.min_count;
        let w = self.cfg.backoff_weight;

        // Global pool for ion and noise error distributions.
        let mut global_ion_raw = vec![0u64; dist_len];
        let mut global_noise_raw = vec![0u64; dist_len];
        for v in counts.error.values() {
            for (i, &c) in v.iter().enumerate() {
                if i < dist_len {
                    global_ion_raw[i] = global_ion_raw[i].saturating_add(c);
                }
            }
        }
        for v in counts.noise_error.values() {
            for (i, &c) in v.iter().enumerate() {
                if i < dist_len {
                    global_noise_raw[i] = global_noise_raw[i].saturating_add(c);
                }
            }
        }
        let global_ion_norm = normalize_with_pseudo(&global_ion_raw, dist_len, pseudo);
        let global_noise_norm = normalize_with_pseudo(&global_noise_raw, dist_len, pseudo);

        let mut ion_out: FxHashMap<Partition, Vec<f32>> = FxHashMap::default();
        let mut noise_out: FxHashMap<Partition, Vec<f32>> = FxHashMap::default();

        for &part in &template.partitions {
            let ion_raw = counts.error.get(&part).map(|v| v.as_slice()).unwrap_or(&[]);
            let ion_n: u64 = ion_raw.iter().sum();
            let ion_emp = normalize_with_pseudo(ion_raw, dist_len, pseudo);
            let ion_parent: Vec<f32> = prior
                .and_then(|p| p.ion_err_dist_table.get(&part))
                .filter(|d| d.len() == dist_len)
                .cloned()
                .unwrap_or_else(|| global_ion_norm.clone());
            let ion_dist = if ion_n < min_count {
                blend(&ion_emp, &ion_parent, ion_n as f32, w)
            } else {
                ion_emp
            };
            ion_out.insert(part, ion_dist);

            let noise_raw = counts.noise_error.get(&part).map(|v| v.as_slice()).unwrap_or(&[]);
            let noise_n: u64 = noise_raw.iter().sum();
            let noise_emp = normalize_with_pseudo(noise_raw, dist_len, pseudo);
            let noise_parent: Vec<f32> = prior
                .and_then(|p| p.noise_err_dist_table.get(&part))
                .filter(|d| d.len() == dist_len)
                .cloned()
                .unwrap_or_else(|| global_noise_norm.clone());
            let noise_dist = if noise_n < min_count {
                blend(&noise_emp, &noise_parent, noise_n as f32, w)
            } else {
                noise_emp
            };
            noise_out.insert(part, noise_dist);
        }

        (ion_out, noise_out)
    }

    // -----------------------------------------------------------------------
    // ion existence table
    // -----------------------------------------------------------------------

    fn build_existence_table(
        &self,
        counts: &CountStats,
        template: &Param,
        prior: Option<&Param>,
    ) -> FxHashMap<Partition, Vec<f32>> {
        const N_EX: usize = 4;
        let pseudo = self.cfg.pseudo;
        let min_count = self.cfg.min_count;
        let w = self.cfg.backoff_weight;

        // Global pool.
        let mut global_raw = [0u64; N_EX];
        for (&(_part, idx), &c) in &counts.existence {
            if (idx as usize) < N_EX {
                global_raw[idx as usize] = global_raw[idx as usize].saturating_add(c);
            }
        }
        let global_norm = normalize_with_pseudo(&global_raw, N_EX, pseudo);

        let mut out: FxHashMap<Partition, Vec<f32>> = FxHashMap::default();
        for &part in &template.partitions {
            let raw: Vec<u64> = (0..N_EX as u32)
                .map(|idx| counts.existence.get(&(part, idx)).copied().unwrap_or(0))
                .collect();
            let n: u64 = raw.iter().sum();
            let emp = normalize_with_pseudo(&raw, N_EX, pseudo);
            let parent: Vec<f32> = prior
                .and_then(|p| p.ion_existence_table.get(&part))
                .filter(|d| d.len() == N_EX)
                .cloned()
                .unwrap_or_else(|| global_norm.clone());
            let dist = if n < min_count {
                blend(&emp, &parent, n as f32, w)
            } else {
                emp
            };
            out.insert(part, dist);
        }
        out
    }
}

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

/// MS-GF+-style widening rank-window smoothing of a rank-distribution vector.
///
/// `dist` has length `max_rank + 1`: indices `0..max_rank` are observed ranks
/// 1..max_rank, index `max_rank` is the missing-ion sentinel. The half-width of
/// the moving average widens with rank — the discriminative head (ranks 1-3) and
/// the missing-ion sentinel are left untouched; the noisy tail is averaged with
/// neighbors. The result is renormalized to sum 1.0.
///
/// Window schedule (published method, smoothingRanks {3,5,10,20,50}):
/// rank<3 -> hw 0 (no smoothing), <5 -> 1, <10 -> 2, <20 -> 3, <50 -> 4, else -> 5.
pub fn smooth_rank_window(dist: &[f32], max_rank: usize) -> Vec<f32> {
    let n = dist.len();
    let last = max_rank.min(n.saturating_sub(1));
    let mut out = dist.to_vec();
    let halfwidth = |r: usize| -> usize {
        if r < 3 { 0 } else if r < 5 { 1 } else if r < 10 { 2 }
        else if r < 20 { 3 } else if r < 50 { 4 } else { 5 }
    };
    for (i, slot) in out.iter_mut().enumerate().take(last) {
        let hw = halfwidth(i);
        if hw == 0 { continue; }
        let lo = i.saturating_sub(hw);
        let hi = (i + hw + 1).min(last);
        let mut s = 0.0f32;
        let mut c = 0usize;
        for v in dist.iter().take(hi).skip(lo) { s += *v; c += 1; }
        if c > 0 { *slot = s / c as f32; }
    }
    let tot: f32 = out.iter().sum();
    if tot > 0.0 {
        for x in &mut out { *x /= tot; }
    }
    out
}

/// Normalise a raw-count slice into a probability vector of length `n_slots`
/// with Laplace smoothing (`pseudo` added to every bin before dividing).
///
/// Short slices are zero-padded to `n_slots`.  The result sums to 1.0 (within
/// floating-point error) and every entry is strictly positive (>= pseudo / total).
fn normalize_with_pseudo(raw: &[u64], n_slots: usize, pseudo: f32) -> Vec<f32> {
    let mut v: Vec<f32> = Vec::with_capacity(n_slots);
    for i in 0..n_slots {
        let c = raw.get(i).copied().unwrap_or(0) as f32;
        v.push(c + pseudo);
    }
    let total: f32 = v.iter().sum();
    if total > 0.0 {
        for x in &mut v {
            *x /= total;
        }
    }
    v
}

/// Bayesian blend: `(n * empirical + w * prior) / (n + w)` element-wise.
///
/// Both slices have the same length (guaranteed by the callers above).
fn blend(emp: &[f32], prior: &[f32], n: f32, w: f32) -> Vec<f32> {
    let denom = n + w;
    if denom <= 0.0 {
        return prior.to_vec();
    }
    emp.iter()
        .zip(prior.iter())
        .map(|(&e, &p)| (n * e + w * p) / denom)
        .collect()
}

/// Global pool: per-`IonType`, sum rank counts across all partitions.
fn build_global_pool(
    rank: &FxHashMap<(Partition, IonType), Vec<u64>>,
    n_slots: usize,
) -> FxHashMap<IonType, Vec<u64>> {
    let mut pool: FxHashMap<IonType, Vec<u64>> = FxHashMap::default();
    for (&(_part, ion), v) in rank {
        let entry = pool.entry(ion).or_insert_with(|| vec![0u64; n_slots]);
        if entry.len() < n_slots {
            entry.resize(n_slots, 0);
        }
        for (i, &c) in v.iter().enumerate() {
            if i < n_slots {
                entry[i] = entry[i].saturating_add(c);
            }
        }
    }
    pool
}

/// Segment-collapsed pool: per `(charge, parent_mass.to_bits())`, sum rank
/// counts across all `seg_num` values.
fn build_seg_collapsed(
    rank: &FxHashMap<(Partition, IonType), Vec<u64>>,
    n_slots: usize,
) -> FxHashMap<(i32, u32), FxHashMap<IonType, Vec<u64>>> {
    let mut out: FxHashMap<(i32, u32), FxHashMap<IonType, Vec<u64>>> = FxHashMap::default();
    for (&(part, ion), v) in rank {
        let key = (part.charge, part.parent_mass.to_bits());
        let ion_map = out.entry(key).or_default();
        let entry = ion_map.entry(ion).or_insert_with(|| vec![0u64; n_slots]);
        if entry.len() < n_slots {
            entry.resize(n_slots, 0);
        }
        for (i, &c) in v.iter().enumerate() {
            if i < n_slots {
                entry[i] = entry[i].saturating_add(c);
            }
        }
    }
    out
}

/// Build `Vec<(charge, num_specs)>` from the raw charge-count map.
fn build_charge_hist(charge: &FxHashMap<i32, u64>) -> Vec<(i32, i32)> {
    let mut hist: Vec<(i32, i32)> = charge
        .iter()
        .map(|(&ch, &n)| (ch, n.min(i32::MAX as u64) as i32))
        .collect();
    hist.sort_by_key(|(ch, _)| *ch);
    hist
}

/// Derive (min_charge, max_charge) from the histogram, falling back to the
/// template when the histogram is empty.
fn charge_range(hist: &[(i32, i32)], template: &Param) -> (i32, i32) {
    if hist.is_empty() {
        return (template.min_charge, template.max_charge);
    }
    let min = hist.iter().map(|(c, _)| *c).min().unwrap_or(template.min_charge);
    let max = hist.iter().map(|(c, _)| *c).max().unwrap_or(template.max_charge);
    (min, max)
}
