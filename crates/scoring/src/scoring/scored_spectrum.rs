//! Per-spectrum precomputed state for scoring.
//!
//! Provides peak ranking by intensity + nearest-peak-by-mz lookup, plus
//! precursor-peak filtering before ranking.
//!
//! ## Precursor-peak filtering formula
//!
//! For each `(reduced_charge, offset, tolerance)` entry in
//! `precursor_off_map[charge]`:
//!
//! ```text
//! neutral_mass = (precursor_mz - PROTON) * charge
//! c            = charge - reduced_charge
//! filter_mz    = (neutral_mass + c * PROTON) / c + offset
//! ```
//!
//! Any peak whose m/z is within `tolerance` Da of `filter_mz` is excluded
//! from ranking. `offset` is in m/z space, added after dividing by `c`.
//!
//! Also exposes `prob_peak`, `main_ion`, `node_score`, `edge_score`,
//! and `observed_node_mass` for the GF DP graph traversal.

use std::sync::OnceLock;

use crate::param_model::{IonType, Param, Partition, PrecursorOffsetFrequency};
use crate::scoring::rank_scorer::RankScorer;
use model::mass::nominal_from;
use model::spectrum::Spectrum;

const PROTON: f64 = 1.007_276_49;

/// iter31 P-2: cache the (MSGF_TRACE_IONS && MSGF_TRACE_PEP) env-var probe
/// once instead of calling `env::var_os` twice per `directional_node_score_inner`
/// invocation. The inner loop fires for every (spectrum × split × segment)
/// triple in the score_psm cache build.
fn trace_ions_enabled() -> bool {
    static CELL: OnceLock<bool> = OnceLock::new();
    *CELL.get_or_init(|| {
        std::env::var_os("MSGF_TRACE_IONS").is_some()
            && std::env::var_os("MSGF_TRACE_PEP").is_some()
    })
}

#[derive(Debug, Clone)]
pub struct ScoredSpectrum<'a> {
    spec: &'a Spectrum,
    /// Per-peak rank (1 = highest intensity), aligned with `spec.peaks`
    /// indices. `ranks[i]` is the rank of the peak at index `i` in the
    /// original `spec.peaks` array. Ties broken by ascending m/z.
    /// Peaks filtered out by precursor-peak filtering receive rank `u32::MAX`.
    ///
    /// When deconvolution is applied (see `deconv_peaks`), the active
    /// rank list is `deconv_ranks`, NOT this field. This field is
    /// retained for `nearest_peak_full` / `nearest_peak_rank` which by
    /// design operate on the original spectrum peaks.
    ranks: Vec<u32>,
    /// Deconvoluted peak list when `param.apply_deconvolution = true`.
    /// Each entry is `(mz, intensity)` after charge-reducing multi-charge
    /// isotope clusters to charge-1 mass (`new_mz = ionCharge * mz - (ionCharge - 1) * PROTON`).
    /// Sorted ascending by m/z so binary search lookups stay O(log n).
    /// Mirrors Java's `Spectrum.getDeconvolutedSpectrum` and is consumed
    /// by `directional_node_score_inner` and `observed_node_mass`.
    /// `None` when deconvolution is not applied — callers fall back to
    /// `spec.peaks` / `ranks` (the original spectrum).
    deconv_peaks: Option<Vec<(f64, f32)>>,
    /// Ranks aligned with `deconv_peaks`. Each original peak's rank is
    /// preserved on the deconvoluted peak (Java's
    /// `setRanksOfPeaks` runs BEFORE `getDeconvolutedSpectrum`).
    /// `None` exactly when `deconv_peaks` is `None`.
    deconv_ranks: Option<Vec<u32>>,
    /// Number of peaks that survived precursor-peak filtering (used for
    /// `peak_count_after_filtering`).
    kept_count: usize,
    /// Raw sum of all peak intensities in the original spectrum.
    total_intensity: f64,
    /// Probability that a random m/z bin contains a peak.
    /// `prob_peak = peak_count / max(approxNumBins, 1)` where
    /// `approxNumBins = parentMass / (mme.getValue() * 2)`.
    ///
    /// For `new_without_filtering` (tests / unit use) this is set to a
    /// sentinel value of `1.0` — callers relying on `edge_score` accuracy
    /// should use the `new` constructor with a full `Param`.
    pub(crate) prob_peak: f32,
    /// The "main ion" for this spectrum's precursor partition. Used by
    /// `observed_node_mass` to look up the observed peak closest to a
    /// theoretical node mass. Set to a Prefix(charge=1, offset=0) fallback
    /// when `new_without_filtering` is used, or derived from the scorer's
    /// table when `new` is used.
    pub(crate) main_ion: IonType,
    /// Spectrum-level parent mass (= `(precursor_mz - PROTON) * charge`),
    /// the OBSERVED neutral mass. Used by `score_psm` / `node_score` for
    /// partition + segment selection so that all candidates at this
    /// spectrum see the same partition (a per-spectrum parent_mass,
    /// regardless of any candidate's nominal/iso-offset mass).
    pub(crate) parent_mass: f64,
    /// The charge state used to construct this ScoredSpectrum.
    pub(crate) charge: u8,
    /// Per-segment (partition, paired (ion, log_table)) cache. Precomputed at
    /// ScoredSpectrum construction (constant for this spectrum's
    /// (charge, parent_mass)). Replaces per-call `partition_for` binary
    /// search + `partition_ion_logs` HashMap lookup in
    /// `directional_node_score`.
    ///
    /// Indexed by segment number `[0..num_segments)`. For the test-fixture
    /// constructor `new_without_filtering` (no Param / RankScorer in scope)
    /// the cache is empty; the hot path tolerates length 0 by simply
    /// iterating no segments and returning 0.0.
    segment_partition_cache: Vec<(Partition, Vec<(IonType, Vec<f32>)>)>,
    /// FastScorer-style directional node-score tables indexed by nominal
    /// residue mass. Populated for production `new()` so candidate scoring
    /// can do array lookups instead of recomputing per-split node scores.
    /// Left empty in `new_without_filtering`, where callers fall back to the
    /// exact uncached path.
    prefix_score_cache: Vec<f32>,
    suffix_score_cache: Vec<f32>,
    /// iter36: spectrum-wide cache for `observed_node_mass(node_nominal)`.
    /// Indexed by `node_nominal` (i32 → usize). Each cell uses an f64 sentinel
    /// encoding:
    ///
    ///   - `f64::NEG_INFINITY` → uncached (not yet computed)
    ///   - `f64::INFINITY`     → cached / no peak in tolerance window
    ///   - any finite value    → cached / observed peak mass
    ///
    /// `RefCell` for interior mutability — ScoredSpectrum is constructed and
    /// consumed within a single Rayon worker thread; no cross-thread sharing,
    /// so single-threaded interior mutability is safe. Note: this REMOVES the
    /// `Sync` auto-derived bound on ScoredSpectrum, which is acceptable
    /// because callers only hand out `&ScoredSpectrum` within one thread.
    ///
    /// Without this cache, `observed_node_mass` was 11.56% of Astral wall
    /// (per iter35 perf profile) — each call did a binary_search over peaks
    /// + linear scan. iter33's per-candidate `psm_edge_score` calls it twice
    /// per edge × 9 edges × 16M candidates ≈ 290M times per Astral spectrum,
    /// repeatedly for the same `node_nominal` values.
    observed_mass_cache: std::cell::RefCell<Vec<f64>>,
}

impl<'a> ScoredSpectrum<'a> {
    /// Construct, filtering precursor peaks at offsets from
    /// `param.precursor_off_map[charge]` before ranking. Also computes
    /// `prob_peak` and selects `main_ion` from the scorer.
    ///
    /// `charge` is the precursor charge of `spec`; if `spec.precursor_charge`
    /// is `Some(z)`, callers typically pass `z`; if `None`, pass the charge
    /// being tried by the search loop.
    ///
    /// Any peak whose m/z is within the tolerance of a precursor filter m/z
    /// gets rank `u32::MAX` and is effectively invisible to `nearest_peak_rank`.
    pub fn new(spec: &'a Spectrum, scorer: &RankScorer, charge: u8) -> Self {
        let param = scorer.param();
        let n = spec.peaks.len();

        // Collect filter m/z values from param.precursor_off_map for this charge.
        let filter_entries: &[PrecursorOffsetFrequency] = param
            .precursor_off_map
            .get(&(charge as i32))
            .map(Vec::as_slice)
            .unwrap_or(&[]);

        // Compute each filter m/z:
        // neutral_mass = (precursor_mz - PROTON) * charge
        // c = charge - reduced_charge
        // filter_mz = (neutral_mass + c * PROTON) / c + offset
        let neutral_mass = (spec.precursor_mz - PROTON) * (charge as f64);
        let filter_mzs: Vec<(f64, f64)> = filter_entries
            .iter()
            .filter_map(|pof| {
                let c = (charge as i32 - pof.reduced_charge) as f64;
                if c <= 0.0 {
                    // Would produce division by zero or negative charge; skip.
                    return None;
                }
                let filter_mz = (neutral_mass + c * PROTON) / c + (pof.offset as f64);
                let tol_da = pof.tolerance.as_da(filter_mz);
                Some((filter_mz, tol_da))
            })
            .collect();

        // Determine which peaks survive filtering.
        let mut ranks = vec![u32::MAX; n];
        let mut kept: Vec<(usize, f32, f64)> = Vec::with_capacity(n);
        for (i, &(mz, intensity)) in spec.peaks.iter().enumerate() {
            let filtered = filter_mzs
                .iter()
                .any(|&(fmz, tol)| (mz - fmz).abs() <= tol);
            if !filtered {
                kept.push((i, intensity, mz));
            }
        }

        let kept_count = kept.len();

        // MS2IonCurrent / ion-current-ratio denominator: Java zeroes precursor
        // peak intensities via `Spectrum.filterPrecursorPeaks` BEFORE
        // PSMFeatureFinder.computeSumIonCurrent iterates the spec
        // (NewScoredSpectrum.java:44-45). Those zeroed peaks then contribute
        // 0 to MS2IonCurrent. Rust filters precursor peaks for rank
        // assignment but the original `spec.peaks` is unmodified, so summing
        // it directly OVER-COUNTS by the precursor-peak intensity. Use the
        // kept set (post-precursor-filter) for the running sum, matching
        // Java's effective denominator. (2026-05-19 PIN diff harness flagged
        // MS2IonCurrent as ~1.6x over Java; this is the source.)
        let total_intensity: f64 = kept.iter().map(|&(_, intensity, _)| intensity as f64).sum();

        // Ranks must be computed BEFORE the FastScorer cache below reads them.
        // The cache calls `directional_node_score_inner(&ranks, ...)` which
        // feeds into `nearest_peak_rank_in` to determine which rank-slot's
        // log score to use. If ranks were all u32::MAX at that point every
        // matched ion would pick the LAST rank slot, producing systematically
        // wrong scores (negative RawScores, near-zero Percolator @ 1% FDR).
        kept.sort_by(|a, b| {
            b.1.partial_cmp(&a.1)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal))
        });
        for (rank_minus_one, &(orig_idx, _, _)) in kept.iter().enumerate() {
            ranks[orig_idx] = (rank_minus_one + 1) as u32;
        }

        let parent_mass = neutral_mass; // = (precursor_mz - PROTON) * charge

        // iter30 C-1: apply Java-parity isotope-cluster deconvolution FIRST,
        // BEFORE prob_peak is computed (Java's `NewScoredSpectrum.java:76-88`
        // does deconv first, then probPeak from the post-deconv spectrum).
        //
        // No `charge > 2` guard — Java's `applyDeconvolution` is unconditional;
        // `deconvolute_spectrum` is a no-op for charge ≤ 2 because its inner
        // loop `for ion_charge_i in 2..charge.min(4)` runs zero iterations.
        // The guard previously skipped deconvolution for Astral charge-2 HCD
        // spectra (a large fraction of the data), introducing a per-spectrum
        // divergence in both `prob_peak` and the prefix/suffix node-score
        // cache.
        let (deconv_peaks, deconv_ranks): (Option<Vec<(f64, f32)>>, Option<Vec<u32>>) =
            if param.apply_deconvolution {
                let tol = param.deconvolution_error_tolerance as f64;
                let (dp, dr) = deconvolute_spectrum(&spec.peaks, &ranks, charge, tol);
                (Some(dp), Some(dr))
            } else {
                (None, None)
            };

        // iter30 C-2: compute prob_peak from the ACTIVE peak list (post-deconv
        // if applied; else kept_count). Java: `probPeak = spec.size() /
        // max(approxNumBins, 1)` where `spec` is the post-deconv spectrum
        // (`NewScoredSpectrum.java:83-88`).
        //
        // parent_mass    = (precursor_mz - PROTON) * charge
        // approxNumBins  = parent_mass / (mme.raw_value() * 2)
        // prob_peak      = max(active_count, 1) / max(approxNumBins, 1)
        let mme_raw = param.mme.raw_value();
        let approx_num_bins = if mme_raw > 0.0 { parent_mass / (mme_raw * 2.0) } else { 1.0 };
        let active_count = match &deconv_peaks {
            Some(dp) => dp.len(),
            None => kept_count,
        };
        let peak_count = if active_count == 0 { 1 } else { active_count } as f64;
        let prob_peak = (peak_count / approx_num_bins.max(1.0)) as f32;

        // Select main_ion: per-partition main ion for (charge, parent_mass, last_seg).
        let last_seg = (param.num_segments - 1).max(0) as usize;
        let part = param.partition_for(charge, parent_mass, last_seg);
        let main_ion = main_ion_from_param(param, part);

        // Precompute the (partition, paired (ion, log_table)) for every
        // segment. This is constant for this spectrum's (charge,
        // parent_mass), so caching here removes a `partition_for` binary
        // search + `partition_ion_logs` HashMap lookup from every call to
        // `directional_node_score`. `partition_ion_logs` returns a
        // borrowed slice; `.to_vec()` clones it to owned so the cache can
        // outlive the borrow on `scorer`.
        let num_segs = param.num_segments.max(0) as usize;
        let segment_partition_cache: Vec<(Partition, Vec<(IonType, Vec<f32>)>)> = (0..num_segs)
            .map(|seg| {
                let p = param.partition_for(charge, parent_mass, seg);
                let logs = scorer.partition_ion_logs(&p).to_vec();
                (p, logs)
            })
            .collect();

        let cache_len = (nominal_from(parent_mass).max(0) as usize) + 1;
        let mut prefix_score_cache = vec![0.0; cache_len];
        let mut suffix_score_cache = vec![0.0; cache_len];
        // Choose the active peak list / rank list ONCE, then reuse for the
        // whole cache fill. When deconvolution was applied, the cache is
        // built against the charge-reduced spectrum (matching Java).
        let (cache_peaks, cache_ranks): (&[(f64, f32)], &[u32]) =
            match (&deconv_peaks, &deconv_ranks) {
                (Some(dp), Some(dr)) => (dp.as_slice(), dr.as_slice()),
                _ => (spec.peaks.as_slice(), ranks.as_slice()),
            };
        for nominal_mass in 1..cache_len {
            let node_nominal = nominal_mass as f64;
            prefix_score_cache[nominal_mass] = Self::directional_node_score_inner(
                cache_peaks,
                cache_ranks,
                &segment_partition_cache,
                scorer,
                node_nominal,
                true,
                charge,
                parent_mass,
            );
            suffix_score_cache[nominal_mass] = Self::directional_node_score_inner(
                cache_peaks,
                cache_ranks,
                &segment_partition_cache,
                scorer,
                node_nominal,
                false,
                charge,
                parent_mass,
            );
        }

        // iter36: spectrum-wide observed_node_mass cache.
        // Size = (parent_nominal + 1) so node_nominal in [0, parent_nominal]
        // is directly indexable. Cap at parent_mass (in Da) → nominal mass
        // ≈ parent_mass × INTEGER_MASS_SCALER. Add small margin for isotope
        // tolerance + rounding.
        let parent_nominal = nominal_from(parent_mass).max(0) as usize;
        let observed_mass_cache = std::cell::RefCell::new(vec![f64::NEG_INFINITY; parent_nominal + 1]);

        Self {
            spec,
            ranks,
            kept_count,
            total_intensity,
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
        }
    }

    /// Constructor that skips precursor-peak filtering. Convenient for
    /// tests; preserves the simpler unfiltered API.
    ///
    /// Sets `prob_peak = 1.0` and `main_ion = Prefix(charge=1, offset=0)`
    /// as sentinels. For accurate `edge_score` computations, use `new`.
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
        let main_ion = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
        // Sentinel: derive parent_mass from spec.precursor_mz with charge defaulted to
        // spec.precursor_charge or 2. Tests using this constructor are typically
        // not sensitive to partition selection.
        let charge = spec.precursor_charge.map(|z| z.max(1) as u8).unwrap_or(2);
        let parent_mass = (spec.precursor_mz - PROTON) * (charge as f64);
        // No Param / RankScorer in scope; segment_partition_cache is left
        // empty. `directional_node_score` tolerates an empty cache: the
        // outer loop iterates zero times and the function returns 0.0.
        // The test-fixture path doesn't need the per-segment optimization.
        let segment_partition_cache: Vec<(Partition, Vec<(IonType, Vec<f32>)>)> = Vec::new();
        let prefix_score_cache: Vec<f32> = Vec::new();
        let suffix_score_cache: Vec<f32> = Vec::new();
        Self::rank_kept(
            spec, kept, kept_count, ranks, prob_peak, main_ion, parent_mass, charge,
            segment_partition_cache,
            prefix_score_cache,
            suffix_score_cache,
        )
    }

    /// Shared ranking logic: sort `kept` by intensity DESC / mz ASC and
    /// write ranks back into the `ranks` vec. Returns the finished
    /// `ScoredSpectrum`.
    fn rank_kept(
        spec: &'a Spectrum,
        mut kept: Vec<(usize, f32, f64)>,
        kept_count: usize,
        mut ranks: Vec<u32>,
        prob_peak: f32,
        main_ion: IonType,
        parent_mass: f64,
        charge: u8,
        segment_partition_cache: Vec<(Partition, Vec<(IonType, Vec<f32>)>)>,
        prefix_score_cache: Vec<f32>,
        suffix_score_cache: Vec<f32>,
    ) -> Self {
        let total_intensity: f64 = kept.iter().map(|&(_, intensity, _)| intensity as f64).sum();
        kept.sort_by(|a, b| {
            // Higher intensity first; if equal, lower m/z first.
            b.1.partial_cmp(&a.1)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal))
        });
        for (rank_minus_one, &(orig_idx, _, _)) in kept.iter().enumerate() {
            ranks[orig_idx] = (rank_minus_one + 1) as u32;
        }
        Self {
            spec,
            ranks,
            kept_count,
            total_intensity,
            prob_peak,
            main_ion,
            parent_mass,
            charge,
            segment_partition_cache,
            prefix_score_cache,
            suffix_score_cache,
            deconv_peaks: None,
            deconv_ranks: None,
            // iter36: empty cache for test fixtures (rank_kept path). All
            // observed_node_mass queries fall through to compute on every call.
            observed_mass_cache: std::cell::RefCell::new(Vec::new()),
        }
    }

    /// Returns `true` if the main ion is a prefix ion (b-ion direction),
    /// `false` if it is a suffix ion (y-ion direction). Used by
    /// `PrimitiveAaGraph` to decide which end is the graph source.
    pub fn main_ion_direction(&self) -> bool {
        self.main_ion.is_prefix()
    }

    /// Return the active peak list and aligned rank vector for the per-node
    /// scoring path. When deconvolution is applied (HCD/CID-HighRes/ETD/QExactive
    /// params with `apply_deconvolution=true`), this returns the
    /// charge-reduced peak list. Otherwise it returns the original spectrum's
    /// peaks and their ranks.
    #[inline]
    fn active_peaks_and_ranks(&self) -> (&[(f64, f32)], &[u32]) {
        match (&self.deconv_peaks, &self.deconv_ranks) {
            (Some(peaks), Some(ranks)) => (peaks.as_slice(), ranks.as_slice()),
            _ => (self.spec.peaks.as_slice(), self.ranks.as_slice()),
        }
    }

    /// Spectrum-level parent mass (= `(precursor_mz - PROTON) * charge`).
    /// This is the OBSERVED neutral mass of the spectrum at the charge
    /// state used to construct this `ScoredSpectrum`, NOT the candidate
    /// peptide's mass.
    pub fn parent_mass(&self) -> f64 {
        self.parent_mass
    }

    /// Return a cached `round(prefix_score + suffix_score)` split score when
    /// both nominal masses are in-bounds for this spectrum's FastScorer-style
    /// tables. Returns `None` when the cache is unavailable or either index is
    /// out of range, allowing callers to fall back to the exact node-score path.
    pub fn cached_split_score(&self, prefix_nominal: i32, suffix_nominal: i32) -> Option<i32> {
        if prefix_nominal < 0 || suffix_nominal < 0 {
            return None;
        }
        let pref = *self.prefix_score_cache.get(prefix_nominal as usize)?;
        let suff = *self.suffix_score_cache.get(suffix_nominal as usize)?;
        Some((pref + suff).round() as i32)
    }

    /// Trace-only accessor: raw `prefix_score_cache[prefix_nominal]` if in
    /// range, mirroring Java's `FastScorer.prefixScore[prefixMass]`. Returns
    /// `None` for an out-of-range index or an empty cache (the
    /// `new_without_filtering` test path leaves the cache empty). This is
    /// consumed by `score_psm`'s trace branch only; the hot scoring path
    /// continues to read through `cached_split_score`.
    pub fn cached_prefix_score(&self, prefix_nominal: i32) -> Option<f32> {
        if prefix_nominal < 0 {
            return None;
        }
        self.prefix_score_cache.get(prefix_nominal as usize).copied()
    }

    /// Trace-only accessor companion to [`cached_prefix_score`]. Mirrors
    /// Java's `FastScorer.suffixScore[suffixMass]`.
    pub fn cached_suffix_score(&self, suffix_nominal: i32) -> Option<f32> {
        if suffix_nominal < 0 {
            return None;
        }
        self.suffix_score_cache.get(suffix_nominal as usize).copied()
    }

    /// Charge state used when this `ScoredSpectrum` was constructed.
    pub fn charge(&self) -> u8 {
        self.charge
    }

    /// For tests only: mutate the main_ion to a different ion type.
    /// Allows test code to exercise both prefix and suffix direction paths.
    /// Not gated by `#[cfg(test)]` so that integration tests in `tests/`
    /// can call it (integration test binaries compile the crate without
    /// the test cfg).
    pub fn set_main_ion_for_test(&mut self, ion: IonType) {
        self.main_ion = ion;
    }

    /// Total number of peaks in the original spectrum (before any filtering).
    pub fn peak_count(&self) -> usize {
        self.spec.peaks.len()
    }

    /// Number of peaks that survived precursor-peak filtering (and were ranked).
    pub fn peak_count_after_filtering(&self) -> usize {
        self.kept_count
    }

    /// Total intensity of all peaks in the original spectrum (before any
    /// filtering). This is the raw MS2 ion current used by the Java
    /// `PSMFeatureFinder.computeSumIonCurrent()` method.
    ///
    /// Returns 0.0 for an empty spectrum.
    pub fn total_intensity(&self) -> f64 {
        self.total_intensity
    }

    /// Find the **highest-intensity** peak within `tolerance_da` of
    /// `target_mz` and return `(rank, intensity, peak_mz)`, or `None` if
    /// no peak falls within the window. Filtered-out peaks
    /// (rank == `u32::MAX`) are never returned.
    ///
    /// Intensity-max selection (same semantics as `nearest_peak_rank`).
    /// Used by `compute_psm_features` for ion-current ratio and
    /// error-stat columns. Closest-by-m/z selection would disagree with
    /// the intensity-comparator selection and affect PIN feature columns
    /// even when the rank lookup matches.
    pub fn nearest_peak_full(&self, target_mz: f64, tolerance_da: f64) -> Option<(u32, f32, f64)> {
        if self.spec.peaks.is_empty() {
            return None;
        }
        let lo_mz = target_mz - tolerance_da;
        let hi_mz = target_mz + tolerance_da;
        let start = self.spec.peaks.partition_point(|&(mz, _)| mz < lo_mz);
        let mut best: Option<(usize, f32)> = None; // (peak_index, intensity)
        for i in start..self.spec.peaks.len() {
            let (mz, intensity) = self.spec.peaks[i];
            if mz > hi_mz {
                break;
            }
            if self.ranks[i] == u32::MAX {
                continue;
            }
            if best.as_ref().map_or(true, |(_, best_int)| intensity > *best_int) {
                best = Some((i, intensity));
            }
        }
        best.map(|(i, _)| {
            let (peak_mz, intensity) = self.spec.peaks[i];
            (self.ranks[i], intensity, peak_mz)
        })
    }

    /// Find the **highest-intensity** peak within `tolerance_da` of `target_mz`,
    /// and return its rank. Returns `None` if no peak falls within the window.
    ///
    /// Returns the most-intense peak in the window (intensity-max
    /// selection); the caller then reads the peak's rank. For LowRes CID
    /// with mme = 0.5 Da, windows frequently contain multiple peaks;
    /// selecting the most-intense matches rank-based scoring exactly.
    /// Closest-by-m/z selection yields systematically higher (worse) rank
    /// numbers and is a dominant cause of top-1 flips.
    ///
    /// Filtered-out peaks (rank == `u32::MAX`) are never returned.
    ///
    /// `spec.peaks` is sorted ascending by m/z (the MGF reader guarantees
    /// this). Binary search (`partition_point`) locates the first
    /// peak with `mz >= target_mz - tolerance_da`; the forward scan then
    /// stops as soon as `mz > target_mz + tolerance_da`, so only the O(k)
    /// peaks in the window are visited.
    pub fn nearest_peak_rank(&self, target_mz: f64, tolerance_da: f64) -> Option<u32> {
        if self.spec.peaks.is_empty() {
            return None;
        }
        let lo_mz = target_mz - tolerance_da;
        let hi_mz = target_mz + tolerance_da;
        // Find first peak with mz >= lo_mz via binary search.
        let start = self.spec.peaks.partition_point(|&(mz, _)| mz < lo_mz);
        // Track (peak_index, intensity); pick max intensity (intensity-comparator selection).
        let mut best: Option<(usize, f32)> = None;
        for i in start..self.spec.peaks.len() {
            let (mz, intensity) = self.spec.peaks[i];
            if mz > hi_mz {
                break;
            }
            // Skip filtered-out peaks.
            if self.ranks[i] == u32::MAX {
                continue;
            }
            if best.as_ref().map_or(true, |(_, best_int)| intensity > *best_int) {
                best = Some((i, intensity));
            }
        }
        best.map(|(i, _)| self.ranks[i])
    }

    /// Return the rank of the peak at index `idx`, or `None` if the peak has
    /// been filtered out (rank == `u32::MAX`) or `idx` is out of bounds.
    ///
    /// Primarily used by tests to compare binary-search results against
    /// brute-force linear scans.
    #[cfg(test)]
    pub(crate) fn peak_rank_at(&self, idx: usize) -> Option<u32> {
        let r = *self.ranks.get(idx)?;
        if r == u32::MAX { None } else { Some(r) }
    }

    // -----------------------------------------------------------------------
    // GF DP scoring methods
    // -----------------------------------------------------------------------

    /// Combined node score for a peptide split position:
    /// `round(prefix_score(prefix_nominal) + suffix_score(suffix_nominal))`.
    ///
    /// `prefix_nominal` and `suffix_nominal` are the float node masses in Da
    /// (not integer nominal-mass indices). `parent_mass` is the precursor
    /// neutral mass. `fragment_tolerance_da` is the m/z window for peak lookup.
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
            prefix_nominal, true, scorer, charge, parent_mass, fragment_tolerance_da,
        );
        let suff = self.directional_node_score(
            suffix_nominal, false, scorer, charge, parent_mass, fragment_tolerance_da,
        );
        (pref + suff).round() as i32
    }

    /// Score for a single directional (prefix or suffix) node at `nominal_mass`.
    ///
    /// **Fragment tolerance:** the per-ion peak-lookup window comes from
    /// `scorer.param().mme.as_da(theo_mz)`. The `fragment_tolerance_da`
    /// argument is retained for backward compat but **ignored** for ion
    /// matching — the param's `mme` is the source of truth here, not a
    /// global search-level fragment tolerance. A hardcoded 0.5 Da happens
    /// to match LowRes CID's mme but is wrong for any other
    /// instrument/protocol.
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
        )
    }

    fn directional_node_score_inner(
        peaks: &[(f64, f32)],
        ranks: &[u32],
        segment_partition_cache: &[(Partition, Vec<(IonType, Vec<f32>)>)],
        scorer: &RankScorer,
        nominal_mass: f64,
        is_prefix: bool,
        charge: u8,
        parent_mass: f64,
    ) -> f32 {
        use crate::param_model::IonType;
        let param = scorer.param();
        let mme = &param.mme;
        let max_rank = scorer.max_rank();
        let max_rank_idx = max_rank as usize;
        let num_segs = param.num_segments as usize;
        let mut total = 0.0_f32;
        let use_cache = !segment_partition_cache.is_empty();
        // Trace gating: only fire when explicitly enabled AND a peptide-trace
        // env var is set (matches `score_psm`'s gating). iter31 P-2: cache the
        // env probe in a OnceLock — was firing `env::var_os` twice per call,
        // which on Astral runs is ~hundreds of millions of acquisitions of the
        // global env lock.
        let trace_ions = trace_ions_enabled();
        for seg in 0..num_segs {
            let ion_logs_slice: &[(IonType, Vec<f32>)] = if use_cache {
                segment_partition_cache[seg].1.as_slice()
            } else {
                let p = param.partition_for(charge, parent_mass, seg);
                scorer.partition_ion_logs(&p)
            };
            if trace_ions {
                eprintln!(
                    "TRACE_RUST_IONS\tnominal={:.3}\tis_prefix={}\tseg={}\tnum_ions={}",
                    nominal_mass, is_prefix, seg, ion_logs_slice.len()
                );
            }
            for (ion, logs) in ion_logs_slice {
                let theo_mz = match (is_prefix, *ion) {
                    (true, IonType::Prefix { .. }) => ion.mz(nominal_mass),
                    (false, IonType::Suffix { .. }) => ion.mz(nominal_mass),
                    _ => continue,
                };
                if param.segment_num(theo_mz, parent_mass) != seg {
                    continue;
                }
                let tol_da = mme.as_da(theo_mz);
                let score = match nearest_peak_rank_in(peaks, ranks, theo_mz, tol_da) {
                    Some(rank) => {
                        let idx = rank.min(max_rank).max(1) as usize - 1;
                        if idx < logs.len() { logs[idx] } else { 0.0 }
                    }
                    None => {
                        if max_rank_idx < logs.len() { logs[max_rank_idx] } else { 0.0 }
                    }
                };
                total += score;
            }
        }
        total
    }

    /// Return the observed node mass for `node_nominal`, or `None` if no
    /// peak is near the theoretical m/z of the main ion.
    ///
    /// Computes `theo_mz = main_ion.mz(node_mass)`, then returns
    /// `main_ion.mass_from_mz(peak_mz)` for the highest-intensity peak
    /// within `mme.as_da(theo_mz)` of `theo_mz`. Returns `Some(0.0)`
    /// at the source node by convention.
    pub fn observed_node_mass(
        &self,
        node_nominal: i32,
        scorer: &RankScorer,
        charge: u8,
        _parent_mass: f64,
    ) -> Option<f64> {
        let _ = charge; // not needed in formula; kept for API symmetry
        if node_nominal == 0 {
            // Source node mass is exactly 0 by convention.
            return Some(0.0);
        }

        // iter36: check spectrum-wide cache first.
        //
        // Sentinel encoding in self.observed_mass_cache:
        //   NEG_INFINITY → uncached, compute now
        //   INFINITY     → cached / no peak found in tolerance window
        //   finite       → cached observed peak mass
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
                // NEG_INFINITY → fall through to compute.
            }
        }

        let theo_mz = self.main_ion.mz(node_nominal as f64);
        let tol_da = scorer.param().mme.as_da(theo_mz);
        // Select the highest-intensity peak within [theo_mz - tol_da, theo_mz + tol_da].
        // Intensity-comparator selection: pick the maximum-intensity peak in the window.
        // Skip filtered peaks (ranks[i] == u32::MAX).
        // Uses the deconvoluted peak list when `param.apply_deconvolution = true` —
        // edge scoring lives downstream of node scoring and must see the same peaks.
        let (peaks, ranks) = self.active_peaks_and_ranks();
        let lo_mz = theo_mz - tol_da;
        let hi_mz = theo_mz + tol_da;
        let start = peaks.partition_point(|&(mz, _)| mz < lo_mz);
        let mut best_peak_mz: Option<(f64, f32)> = None; // (mz, intensity)
        for i in start..peaks.len() {
            let (mz, intensity) = peaks[i];
            if mz > hi_mz {
                break;
            }
            if ranks[i] == u32::MAX {
                continue;
            }
            if best_peak_mz.as_ref().map_or(true, |&(_, best_int)| intensity > best_int) {
                best_peak_mz = Some((mz, intensity));
            }
        }
        let result = best_peak_mz.map(|(peak_mz, _)| self.main_ion.mass_from_mz(peak_mz));

        // iter36: store result in the spectrum-wide cache. Only if idx fits.
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

    /// Edge score for the GF DP.
    ///
    /// If `param.ion_existence_table` is empty (edge scoring not supported),
    /// returns 0. Otherwise:
    ///   1. Look up observed node masses for `cur_nominal` and `prev_nominal`.
    ///   2. `ion_existence_index` = (cur observed?) + 2*(prev observed?).
    ///   3. `score = ion_existence_score(part, idx, prob_peak)`.
    ///   4. If `idx == 3` (both observed), also add `error_score(cur_mass - prev_mass - theo_aa_mass)`.
    ///   5. Return `round(score) as i32`.
    pub fn edge_score(
        &self,
        cur_nominal: i32,
        prev_nominal: i32,
        theo_aa_mass: f64,
        scorer: &RankScorer,
        charge: u8,
        parent_mass: f64,
    ) -> i32 {
        // supportEdgeScores() ↔ errorScalingFactor != 0.
        if scorer.param().error_scaling_factor == 0 {
            return 0;
        }
        if scorer.param().ion_existence_table.is_empty() {
            return 0;
        }

        // 1. Observed masses for cur and prev nodes.
        let cur_mass = self.observed_node_mass(cur_nominal, scorer, charge, parent_mass);
        let prev_mass = self.observed_node_mass(prev_nominal, scorer, charge, parent_mass);

        // 2. ion_existence_index: 1 if cur observed, +2 if prev observed.
        let mut idx = 0usize;
        if cur_mass.is_some() { idx += 1; }
        if prev_mass.is_some() { idx += 2; }

        // 3. Partition for this spectrum — Java uses the "last segment" partition
        //    stored at construction time.
        //
        // iter38 P-9b: per-edge `param.partition_for(charge, parent_mass, last_seg)`
        // was 3.26% of Astral wall (~144M calls under iter33's per-candidate
        // edge scoring). The partition is constant for this ScoredSpectrum's
        // `(charge, parent_mass)` and is already cached in
        // `segment_partition_cache`. Use that instead of re-running the binary
        // search per edge.
        let last_seg = (scorer.param().num_segments - 1).max(0) as usize;
        let part = match self.segment_partition_cache.get(last_seg) {
            Some((p, _)) => *p,
            None => scorer.param().partition_for(charge, parent_mass, last_seg),
        };

        // 4. Ion existence score.
        let mut s = scorer.ion_existence_score(part, idx, self.prob_peak);

        // 5. If both observed, add error score.
        if idx == 3 {
            let delta = cur_mass.unwrap() - prev_mass.unwrap() - theo_aa_mass;
            s += scorer.error_score(part, delta as f32);
        }

        s.round() as i32
    }
}

fn nearest_peak_rank_in(peaks: &[(f64, f32)], ranks: &[u32], target_mz: f64, tolerance_da: f64) -> Option<u32> {
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
        if best.as_ref().map_or(true, |(_, best_int)| intensity > *best_int) {
            best = Some((i, intensity));
        }
    }
    best.map(|(i, _)| ranks[i])
}

/// Java-parity isotope-cluster deconvolution.
///
/// Mirrors `Spectrum.getDeconvolutedSpectrum(toleranceBetweenIsotopes)` in
/// `astral-speed/src/main/java/edu/ucsd/msjava/msutil/Spectrum.java`.
///
/// Input is the spectrum's peak list (sorted ascending by m/z) plus the
/// rank vector aligned with it (rank 1 = highest intensity; `u32::MAX`
/// for filtered peaks). Returns `(peaks, ranks)` of the deconvoluted
/// spectrum, sorted ascending by m/z.
///
/// Algorithm: for each peak `p[i]` (not already consumed), look for a
/// matching +1/ionCharge isotope `p[j]`. If found at `ionCharge ∈ {2, 3}`
/// (and `ionCharge < precursor_charge`), charge-reduce all clustered
/// peaks (`new_mz = ionCharge * mz - (ionCharge - 1) * PROTON`) and look
/// forward for a +2/ionCharge third isotope. Ranks are preserved
/// per-peak because Java's `setRanksOfPeaks` runs BEFORE deconvolution.
///
/// `precursor_charge` is the spectrum's precursor charge (matches Java's
/// `this.getCharge()`). For `precursor_charge <= 2`, no charge-reduction
/// candidates exist (loop `2 < charge` is empty), so the output equals
/// the input modulo a mass-sort.
fn deconvolute_spectrum(
    peaks: &[(f64, f32)],
    ranks: &[u32],
    precursor_charge: u8,
    tol: f64,
) -> (Vec<(f64, f32)>, Vec<u32>) {
    // Java: Composition.ISOTOPE = C13 - C ≈ 1.00335483.
    const ISOTOPE: f64 = 1.003_354_83;
    // Java: (Composition.C14 - Composition.C13) ≈ 0.999_886_17.
    const C14_MINUS_C13: f64 = 0.999_886_17;

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

        // Java's inner loop: `for (ionCharge = 2; ionCharge < charge && ionCharge < 4; ionCharge++)`
        for ion_charge_i in 2..charge_i32.min(4) {
            let ion_charge = ion_charge_i as f64;
            let expected_diff = ISOTOPE / ion_charge;
            let mut is_deconvoluted = false;
            // Look forward for p2 = p1's +1 isotope.
            for j in (i + 1)..n {
                let (p2_mz, p2_int) = peaks[j];
                let diff = p2_mz - p_mz - expected_diff;
                if diff > -tol && diff < tol {
                    // Match: charge-reduce p1 (mutate locally for output) and p2.
                    ignore[j] = true;
                    let p_new_mz = ion_charge * p_mz - (ion_charge - 1.0) * PROTON;
                    let p2_new_mz = ion_charge * p2_mz - (ion_charge - 1.0) * PROTON;
                    // Save p1's new mass; we'll push it after the inner loop
                    // (Java does `deconvSpec.add(p)` at the end of the outer loop).
                    p_mz = p_new_mz;
                    is_deconvoluted = true;

                    // Look for p3 = p2's +1 isotope (uses C14_MINUS_C13 / ion_charge).
                    let p3_diff_expected = C14_MINUS_C13 / ion_charge;
                    for k in (j + 1)..n {
                        let (p3_mz, p3_int) = peaks[k];
                        let diff2 = p3_mz - p2_mz - p3_diff_expected;
                        if diff2 > -tol && diff2 < tol {
                            ignore[k] = true;
                            let p3_new_mz =
                                ion_charge * p3_mz - (ion_charge - 1.0) * PROTON;
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
        // Add p1 (possibly mutated) to output.
        out.push((p_mz, p_int, p_rank));
    }

    // Sort by m/z ascending, ties broken by rank (stable on ties is fine).
    out.sort_by(|a, b| {
        a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut out_peaks: Vec<(f64, f32)> = Vec::with_capacity(out.len());
    let mut out_ranks: Vec<u32> = Vec::with_capacity(out.len());
    for (mz, intensity, rank) in out {
        out_peaks.push((mz, intensity));
        out_ranks.push(rank);
    }
    (out_peaks, out_ranks)
}

/// Select the main ion for `partition` from `param.rank_dist_table`.
///
/// Picks the Prefix ion with the highest freq at rank-1 index (index 0).
/// Falls back to `Prefix { charge: 1, offset_bits: 0 }` if the table is empty.
///
/// Note: selection currently uses per-partition rank-1 prefix-ion frequency
/// from `rank_dist_table`. A fuller selection would aggregate `frag_off_table`
/// across segments and consider all ion types; for HCD these agree, for
/// ETD/ECD they may diverge.
fn main_ion_from_param(param: &Param, partition: crate::param_model::Partition) -> IonType {
    // Mirrors Java's `NewRankScorer.determineIonTypes` (lines 611-640).
    // Aggregates `frag_off_table` frequencies ACROSS ALL SEGMENTS for the same
    // `(charge, parent_mass)` partition and picks the overall highest-frequency
    // ion — regardless of prefix/suffix type. For HCD/QExactive this typically
    // selects a y-ion (suffix), giving `main_ion_direction() = false`.
    //
    // Previous Rust behavior filtered to `is_prefix()` only, forcing direction
    // always true. That mismatched Java's `getMainIonType` and produced wrong
    // EdgeScore values for HCD spectra (iter28 trace: scan 47106 EdgeScore
    // Rust -18 vs Java +8). See
    // `docs/parity-analysis/notes/2026-05-21-iter27-pin-diff.md`.
    let fallback = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
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
    let mut best_ion: Option<IonType> = None;
    let mut best_freq = f32::NEG_INFINITY;
    for (&ion, &freq) in &ion_freq {
        if freq > best_freq {
            best_freq = freq;
            best_ion = Some(ion);
        }
    }
    best_ion.unwrap_or(fallback)
}

#[cfg(test)]
mod tests {
    use super::*;
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
        }
    }

    // --- prob_peak uses raw mme value ---

    /// Verify that `prob_peak` is computed using the raw stored mme value,
    /// not the Da-converted form. For `Tolerance::Ppm(20.0)`:
    ///   Expected: approxNumBins = parent_mass / (mme.raw_value() * 2)
    ///                           = parent_mass / (20.0 * 2)
    ///   NOT:      parent_mass / (as_da(parent_mass) * 2)
    ///                           = parent_mass / (parent_mass * 20e-6 * 2)
    #[test]
    fn prob_peak_uses_raw_mme_value_not_da_converted() {
        use model::activation::ActivationMethod;
        use model::instrument::InstrumentType;
        use crate::param_model::SpecDataType;
        use model::protocol::Protocol;
        use model::tolerance::Tolerance;
        use std::collections::HashMap;

        // Spectrum: precursor_mz=501.00727649 → neutral_mass≈(501.007-PROTON)*2≈1000.0 Da,
        // charge=2.
        let precursor_mz = 501.007_276_49_f64; // ≈ (1000/2) + PROTON
        let s = Spectrum {
            title: "parity_test".into(),
            precursor_mz,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks: vec![(100.0, 1.0), (200.0, 2.0), (300.0, 3.0)],
            activation_method: None,
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
            precursor_off_map: HashMap::new(),
            frag_off_table: HashMap::new(),
            max_rank: 3,
            rank_dist_table: HashMap::new(),
            error_scaling_factor: 0,
            ion_err_dist_table: HashMap::new(),
            noise_err_dist_table: HashMap::new(),
            ion_existence_table: HashMap::new(),
            partition_ion_types_cache: HashMap::new(),
        };

        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new(&s, &scorer, 2);

        // Expected: raw_value = 20.0, parent_mass ≈ (501.007276 - PROTON) * 2.
        let parent_mass = (precursor_mz - PROTON) * 2.0;
        let raw_mme = 20.0_f64;
        let approx_num_bins = parent_mass / (raw_mme * 2.0);
        let expected_prob_peak = (3.0_f64 / approx_num_bins.max(1.0)) as f32;

        // The Da-converted form would be: parent_mass / (parent_mass * 20e-6 * 2) ≈ 25_000.0,
        // giving prob_peak ≈ 3/25000 = 0.00012, not the raw-value result ≈ 3/100 = 0.06.
        let wrong_approx_num_bins = parent_mass / (parent_mass * 20e-6 * 2.0);
        let wrong_prob_peak = (3.0_f64 / wrong_approx_num_bins.max(1.0)) as f32;

        // Sanity: raw and Da results must differ significantly for this to be a meaningful test.
        assert!(
            (expected_prob_peak - wrong_prob_peak).abs() > 0.001,
            "test precondition failed: Ppm raw vs Da-converted did not produce different prob_peak values"
        );

        assert!(
            (ss.prob_peak - expected_prob_peak).abs() < 1e-5,
            "prob_peak={} but expected={} (raw-mme formula). Wrong Da-converted value would be {}",
            ss.prob_peak, expected_prob_peak, wrong_prob_peak
        );
    }

    // --- iter30 C-1 + C-2 deconvolution tests ---

    /// Helper: build a minimal Param with apply_deconvolution toggleable.
    fn deconv_param(apply: bool) -> Param {
        use model::activation::ActivationMethod;
        use model::instrument::InstrumentType;
        use model::protocol::Protocol;
        use model::tolerance::Tolerance;
        use std::collections::HashMap;
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
            precursor_off_map: HashMap::new(),
            frag_off_table: HashMap::new(),
            max_rank: 3,
            rank_dist_table: HashMap::new(),
            error_scaling_factor: 0,
            ion_err_dist_table: HashMap::new(),
            noise_err_dist_table: HashMap::new(),
            ion_existence_table: HashMap::new(),
            partition_ion_types_cache: HashMap::new(),
        }
    }

    /// T-1: For charge-2 spectra with `apply_deconvolution=true`, the deconv
    /// path must be exercised (no early guard) and the output must equal the
    /// input mathematically — because `deconvolute_spectrum`'s inner loop is
    /// `for ion_charge_i in 2..charge.min(4)` which produces an empty range
    /// for charge=2. Iter30 C-1 dropped the `charge > 2` guard so this case
    /// follows Java's unconditional `applyDeconvolution()` branch.
    #[test]
    fn deconv_active_for_charge_2_produces_input_equivalent_peaks() {
        let s = Spectrum {
            title: "deconv_test".into(),
            precursor_mz: 501.007_276_49_f64,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            // Three peaks; none of them is at the deconvolution-tolerance
            // window for charge ≥ 2 since the inner loop is empty for charge=2.
            peaks: vec![(100.0, 1.0), (200.0, 2.0), (300.0, 3.0)],
            activation_method: None,
        };
        let param = deconv_param(true);
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new(&s, &scorer, 2);

        // prob_peak should be derived from the same 3 peaks (deconv is a
        // no-op for charge=2). Active peak count = 3.
        let parent_mass = (s.precursor_mz - PROTON) * 2.0;
        let approx = parent_mass / (20.0_f64 * 2.0);
        let expected = (3.0_f64 / approx.max(1.0)) as f32;
        assert!(
            (ss.prob_peak - expected).abs() < 1e-5,
            "charge=2 deconv-active spectrum: prob_peak={} expected={} (active_count=3)",
            ss.prob_peak, expected
        );
    }

    /// T-2: For charge-3 spectra with `apply_deconvolution=true`, `prob_peak`
    /// MUST be computed from the post-deconvolution peak count, not the
    /// pre-deconvolution kept_count. Java's `NewScoredSpectrum.java:83-88`
    /// derives `probPeak` from `spec.size()` AFTER `spec` is replaced by the
    /// deconvoluted spectrum. Iter30 C-2 enforces this ordering.
    #[test]
    fn deconv_active_for_charge_3_uses_post_deconv_peak_count_for_prob_peak() {
        // Pick a charge=3 spectrum whose peaks include an isotope cluster
        // that the deconvolution algorithm will merge.
        //
        // Construct two peaks at charge=2 m/z separation: ISOTOPE/2 ≈ 0.5017 Da apart
        // and a third for the inner-inner loop. The deconvolution will recognize
        // these as a +2 isotope cluster and reduce them to charge-1 m/z. The
        // OUTPUT peak count differs from the input peak count.
        //
        // For two peaks (the "two-pattern" case), Java's algorithm KEEPS the
        // first, RE-EMITS the second (charge-reduced). So output count == input
        // count when no +3 peak follows. Add a peak FAR from the cluster so it
        // also survives unchanged. The point: even if count is preserved here,
        // the m/z values change → prob_peak's bin model is unaffected since
        // approx_num_bins is parent_mass-derived; what matters is that the
        // value is computed from the active list.
        const ISOTOPE: f64 = 1.003_354_83;
        let p1 = 100.0;
        let p2 = p1 + ISOTOPE / 2.0; // ≈ 100.5017
        let p3 = 500.0; // unrelated peak
        let s = Spectrum {
            title: "deconv_charge3".into(),
            precursor_mz: 401.0,
            precursor_intensity: None,
            precursor_charge: Some(3),
            rt_seconds: None,
            scan: None,
            peaks: vec![(p1, 10.0), (p2, 5.0), (p3, 1.0)],
            activation_method: None,
        };
        let param = deconv_param(true);
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new(&s, &scorer, 3);

        // Whatever the deconvoluted peak count is, prob_peak should match it.
        let active_count = ss.deconv_peaks.as_ref().map(|p| p.len()).unwrap_or(0);
        assert!(active_count >= 1, "deconv_peaks should be populated for charge=3 + apply_deconvolution=true");
        let parent_mass = (401.0 - PROTON) * 3.0;
        let approx = parent_mass / (20.0_f64 * 2.0);
        let expected = (active_count as f64 / approx.max(1.0)) as f32;
        assert!(
            (ss.prob_peak - expected).abs() < 1e-5,
            "charge=3 deconv-active spectrum: prob_peak={} expected={} (post-deconv count={})",
            ss.prob_peak, expected, active_count
        );
    }

    /// T-2b: When `apply_deconvolution=false`, prob_peak follows the pre-deconv
    /// kept count (existing behavior). Sanity check to ensure C-2 doesn't
    /// flip the deconv-off path.
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
        };
        let param = deconv_param(false);
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new(&s, &scorer, 2);

        // No deconv path → active = kept = 4.
        let parent_mass = (s.precursor_mz - PROTON) * 2.0;
        let approx = parent_mass / (20.0_f64 * 2.0);
        let expected = (4.0_f64 / approx.max(1.0)) as f32;
        assert!(
            (ss.prob_peak - expected).abs() < 1e-5,
            "deconv-off: prob_peak={} expected={} (kept_count=4)",
            ss.prob_peak, expected
        );
        assert!(ss.deconv_peaks.is_none(), "deconv_peaks must be None when apply_deconvolution=false");
    }

    // --- observed_node_mass picks highest-intensity ---

    #[test]
    fn observed_node_mass_picks_highest_intensity_peak_in_window() {
        // Two peaks within the MME window of theo_mz; the higher-intensity one wins.
        // tiny_param_with_ions uses Tolerance::Da(0.5) → window ±0.5 Da.
        // main_ion = Prefix { charge: 1, offset_bits: 0 }
        //
        // theo_mz = (node_nominal / INTEGER_MASS_SCALER) / charge + offset
        //         = (100 / 0.999497) / 1 + 0.0 ≈ 100.05028
        //
        // Place two peaks both within ±0.5 of theo_mz ≈ 100.050:
        //   peak A at 100.14 (delta ≈ 0.09, low intensity 1.0) — CLOSER
        //   peak B at 100.44 (delta ≈ 0.39, high intensity 100.0) — FARTHER but HIGHER intensity
        // Highest-intensity wins → peak B.
        use model::mass::INTEGER_MASS_SCALER;
        let node_nominal = 100_i32;
        // theo_mz with offset=0: real_mass / 1 + 0 = nominal / INTEGER_MASS_SCALER
        let theo_mz = node_nominal as f64 / INTEGER_MASS_SCALER as f64;
        let closer_mz = theo_mz + 0.09; // delta 0.09 < 0.39
        let farther_mz = theo_mz + 0.39; // still within ±0.5
        let s = spec(&[(closer_mz, 1.0), (farther_mz, 100.0)]);
        let param = tiny_param_with_ions(); // mme = Da(0.5)
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let result = ss.observed_node_mass(node_nominal, &scorer, 2, 1000.0);
        let result_mass = result.expect("should find a peak in the window");
        // main_ion.mass_from_mz(peak_mz) with offset=0, charge=1: (mz - 0) * 1 = mz
        let expected_mass = farther_mz;
        let wrong_mass = closer_mz;
        assert!(
            (result_mass - expected_mass).abs() < 1e-6,
            "expected highest-intensity (farther) peak mass {expected_mass:.6}, \
             got {result_mass:.6} (closest/wrong would be {wrong_mass:.6})"
        );
    }

    // --- node_score and edge_score ---

    #[test]
    fn node_score_does_not_panic_on_empty_spectrum() {
        // Spectrum with no peaks; every ion is missing → all contributions
        // come from missing_ion_score. With no matching peaks the missing
        // score for Prefix(charge=1) is log(0.001/0.4) < 0, but we also
        // include the suffix side which has no ions. Sum rounds to a small
        // negative.
        let s = spec(&[]);
        let param = tiny_param_with_ions();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let n = ss.node_score(100.0, 900.0, &scorer, 2, 1000.0, 0.5);
        // With empty ion_types_for_segment the suffix side contributes 0,
        // and no suffix ions are in the table → suffix score is 0.
        // The prefix missing-ion score is negative → total rounds negative or 0.
        assert!(n <= 0, "missing-ion score on empty spectrum should be non-positive, got {n}");
    }

    #[test]
    fn node_score_nonzero_when_peak_matches_prefix_ion() {
        // Place a high-intensity peak at the predicted b1 m/z for a node of
        // nominal mass = 100. Prefix ion: Prefix(charge=1, offset=0).
        // theo_mz = (nominal / INTEGER_MASS_SCALER) / 1 + 0
        //         = 100 / 0.999497 ≈ 100.0503
        use model::mass::INTEGER_MASS_SCALER;
        let nominal = 100.0_f64;
        let b1_mz = nominal / INTEGER_MASS_SCALER as f64; // charge=1, offset=0
        let s = spec(&[(50.0, 1.0), (b1_mz, 100.0), (200.0, 2.0)]);
        let param = tiny_param_with_ions();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        // prefix_nominal = 100, suffix_nominal = 900 (doesn't matter, no suffix ions in table).
        let n = ss.node_score(nominal, 900.0, &scorer, 2, 1000.0, 0.5);
        // Peak at b1_mz gets rank 1 (highest intensity = 100.0).
        // node_score(rank=1, Prefix) = log(0.6 / (0.1 * 1)) = log(6) > 0.
        // Total suffix = 0. Round(log(6)) = round(1.79) = 2.
        assert!(n > 0, "expected positive node_score when b-ion peak present, got {n}");
    }

    #[test]
    fn node_score_prefix_only_match() {
        // Only prefix ions in table; suffix side always contributes 0.
        // theo_mz = (nominal / INTEGER_MASS_SCALER) / 1 + 0
        use model::mass::INTEGER_MASS_SCALER;
        let nominal = 57.0_f64; // roughly glycine residue mass
        let mz = nominal / INTEGER_MASS_SCALER as f64; // charge=1, offset=0
        let s = spec(&[(mz, 50.0), (300.0, 1.0)]);
        let param = tiny_param_with_ions();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let n = ss.node_score(nominal, 900.0, &scorer, 2, 1000.0, 0.5);
        // Peak at mz is rank 1. score = log(0.6 / 0.1) = log(6) ≈ 1.79 → rounds to 2.
        assert!(n > 0, "prefix-only match: expected positive score, got {n}");
    }

    #[test]
    fn node_score_no_matching_ions_returns_negative_or_zero() {
        // With a peak far from any ion, all ions are missing → negative score.
        let s = spec(&[(5000.0, 100.0)]); // peak far from any fragment ion
        let param = tiny_param_with_ions();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let n = ss.node_score(100.0, 900.0, &scorer, 2, 1000.0, 0.5);
        // missing_ion_score for Prefix(1) = log(0.001/0.4) < 0 → n <= 0.
        assert!(n <= 0, "missing ion should produce non-positive score, got {n}");
    }

    #[test]
    fn node_score_nominal_mass_zero_prefix_returns_zero() {
        // nominal_mass = 0 is the source node. This impl evaluates
        // ions_for_node(0.0, …) directly. With prefix_nominal=0 and
        // suffix_nominal=1000 (parent mass), and no peaks in the spectrum,
        // the missing-ion score for the Prefix ion governs. The suffix
        // nominal = 1000 > parent_mass → ions_for_node produces no suffix
        // ions for that degenerate case. Net result: non-positive score.
        let s = spec(&[]);
        let param = tiny_param_with_ions();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let n = ss.node_score(0.0, 1000.0, &scorer, 2, 1000.0, 0.5);
        // Score is non-positive (missing-ion penalty applies).
        assert!(n <= 0, "source-node score with empty spectrum should be non-positive, got {n}");
    }

    #[test]
    fn edge_score_returns_zero_when_table_empty() {
        // No ion_existence_table → edge_score returns 0.
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
        // error_scaling_factor == 0 ↔ supportEdgeScores() == false → returns 0.
        let s = spec(&[(100.0, 1.0)]);
        let param = tiny_param_with_ions(); // error_scaling_factor defaults to 0
        assert_eq!(param.error_scaling_factor, 0);
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let e = ss.edge_score(150, 100, 50.0, &scorer, 2, 1000.0);
        assert_eq!(e, 0);
    }

    #[test]
    fn edge_score_nonzero_with_existence_table() {
        // Build a param with error_scaling_factor > 0 and a populated
        // ion_existence_table. Check that edge_score is computed (non-zero).
        use model::activation::ActivationMethod;
        use model::instrument::InstrumentType;
        use crate::param_model::{FragmentOffsetFrequency, SpecDataType};
        use model::protocol::Protocol;
        use model::tolerance::Tolerance;
        use std::collections::HashMap;

        let part = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
        let prefix1 = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
        let noise = IonType::Noise;

        let ion_freqs = vec![0.6_f32, 0.3, 0.05, 0.001];
        let noise_freqs = vec![0.1_f32, 0.2, 0.3, 0.4];

        let mut ion_table: HashMap<IonType, Vec<f32>> = HashMap::new();
        ion_table.insert(prefix1, ion_freqs);
        ion_table.insert(noise, noise_freqs);

        let mut rank_dist_table: HashMap<Partition, HashMap<IonType, Vec<f32>>> = HashMap::new();
        rank_dist_table.insert(part, ion_table);

        let mut frag_off_table = HashMap::new();
        frag_off_table.insert(part, vec![FragmentOffsetFrequency {
            ion_type: prefix1,
            frequency: 0.7,
        }]);

        // error_scaling_factor = 2 → dist_len = 5; ion_existence = 4 entries
        let error_scaling_factor = 2_i32;
        let dist_len = (error_scaling_factor as usize) * 2 + 1;

        let mut ion_err_dist_table: HashMap<Partition, Vec<f32>> = HashMap::new();
        ion_err_dist_table.insert(part, vec![0.1_f32, 0.2, 0.4, 0.2, 0.1]);

        let mut noise_err_dist_table: HashMap<Partition, Vec<f32>> = HashMap::new();
        noise_err_dist_table.insert(part, vec![0.05_f32, 0.1, 0.7, 0.1, 0.05]);

        let mut ion_existence_table: HashMap<Partition, Vec<f32>> = HashMap::new();
        // [nn, ?, ?, yy] = [0.1, 0.3, 0.3, 0.5]
        ion_existence_table.insert(part, vec![0.1_f32, 0.3, 0.3, 0.5]);

        let _ = dist_len; // used for documentation

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
            precursor_off_map: HashMap::new(),
            frag_off_table,
            max_rank: 3,
            rank_dist_table,
            error_scaling_factor,
            ion_err_dist_table,
            noise_err_dist_table,
            ion_existence_table,
            partition_ion_types_cache: HashMap::new(),
        };
        param.rebuild_cache();

        // No peaks in spectrum → cur_mass = None, prev_mass = None → idx = 0 (nn).
        let s = spec(&[]);
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let e = ss.edge_score(150, 100, 50.0, &scorer, 2, 1000.0);
        // ion_existence_score(part, 0, prob_peak): ionExistenceProb[0]=0.1,
        // noiseExistenceProb = (1-p)^2. With many bins prob_peak ≈ 0.
        // log(0.1 / ~1.0) = ~log(0.1) ≈ -2.3 → rounds to -2.
        // Confirm the table is used (non-zero result).
        assert_ne!(e, 0, "edge_score should be nonzero with populated existence table");
    }

    #[test]
    fn directional_node_score_segment_cache_sanity() {
        use crate::param_model::Param;
        use std::path::PathBuf;
        let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        path.push("../../resources/ionstat/CID_LowRes_Tryp.param");
        let param = Param::load_from_file(&path).expect("param loads");
        let scorer = RankScorer::new(&param);
        let peaks: Vec<(f64, f32)> = (0..100).map(|i| (50.0 + i as f64 * 19.5, 100.0 - i as f32)).collect();
        let spec = Spectrum {
            title: "parity".into(), precursor_mz: 800.0, precursor_intensity: None,
            precursor_charge: Some(2), rt_seconds: None, scan: None, peaks,
            activation_method: None,
        };
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        let mut state: u64 = 0xCAFEBABEDEADBEEF;
        let mut next = || { state ^= state << 13; state ^= state >> 7; state ^= state << 17; state };
        for _ in 0..200 {
            let nominal_mass = 100.0 + (next() % 2400) as f64;
            let is_prefix = (next() & 1) == 0;
            let charge = 2 + (next() % 3) as u8;
            let parent_mass = 600.0 + (next() % 2400) as f64;
            let val = ss.directional_node_score(nominal_mass, is_prefix, &scorer, charge, parent_mass, 0.0);
            assert!(val.is_finite() || val == 0.0,
                "non-finite directional_node_score at nominal={nominal_mass} prefix={is_prefix} charge={charge} parent_mass={parent_mass}: {val}");
        }
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
        // Peaks sorted ascending by m/z (the MGF reader guarantees this).
        let s = spec(&[(100.0, 1.0), (200.0, 5.0), (300.0, 3.0)]);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        assert_eq!(ss.peak_count(), 3);
        // Peak at m/z 200 has the highest intensity (5.0) → rank 1.
        // The lookup window of 0.1 should find it.
        assert_eq!(ss.nearest_peak_rank(200.0, 0.1), Some(1));
        // Peak at m/z 300 has intensity 3.0 → rank 2.
        assert_eq!(ss.nearest_peak_rank(300.0, 0.1), Some(2));
        // Peak at m/z 100 has intensity 1.0 → rank 3 (lowest).
        assert_eq!(ss.nearest_peak_rank(100.0, 0.1), Some(3));
    }

    #[test]
    fn nearest_peak_within_tolerance() {
        let s = spec(&[(100.0, 1.0), (200.5, 5.0), (300.0, 3.0)]);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        // Target 200.4 with tol 0.2 → finds peak at 200.5 (within 0.1).
        assert_eq!(ss.nearest_peak_rank(200.4, 0.2), Some(1));
        // Target 200.5 with tol 0.001 → exact match.
        assert_eq!(ss.nearest_peak_rank(200.5, 0.001), Some(1));
        // Target 200.4 with tol 0.05 → outside window, no match.
        assert_eq!(ss.nearest_peak_rank(200.4, 0.05), None);
    }

    #[test]
    fn ties_broken_deterministically() {
        // Two peaks with identical intensity — the lower m/z gets rank 1
        // (matching Java's behavior of sort stability + ties going to
        // earlier-indexed peaks).
        let s = spec(&[(100.0, 5.0), (200.0, 5.0)]);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        // Both peaks should have a defined rank; the test asserts the
        // ranking is total (no two peaks share a rank).
        let r1 = ss.nearest_peak_rank(100.0, 0.1).unwrap();
        let r2 = ss.nearest_peak_rank(200.0, 0.1).unwrap();
        assert_ne!(r1, r2);
        assert!(r1 == 1 || r2 == 1);
        assert!(r1 == 2 || r2 == 2);
    }

    #[test]
    fn closest_among_multiple_in_tolerance() {
        // Multiple peaks within the tolerance window; the closest wins.
        let s = spec(&[(99.5, 1.0), (100.0, 5.0), (100.5, 2.0)]);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        // Target 100.1 with tol 0.6: all three are within. Closest is 100.0 → rank 1.
        assert_eq!(ss.nearest_peak_rank(100.1, 0.6), Some(1));
    }

    #[test]
    fn nearest_peak_rank_matches_linear_scan_on_many_peaks() {
        // Build a spectrum with 100 peaks across 0.0 - 1000.0 m/z, varying intensities.
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
        };
        let ss = ScoredSpectrum::new_without_filtering(&s);

        // For several target m/z values, the binary-search result must match
        // what a brute-force linear scan produces.
        for target in [50.5, 100.0, 250.0, 333.7, 500.5, 750.5, 999.5] {
            let tol = 5.0_f64; // wide window
            let bs_result = ss.nearest_peak_rank(target, tol);
            // Brute force: scan all peaks, pick closest within tolerance.
            let bf_result = {
                let mut best: Option<(usize, f64)> = None;
                for (i, &(mz, _)) in s.peaks.iter().enumerate() {
                    if (mz - target).abs() <= tol
                        && best.as_ref().map_or(true, |(_, d)| (mz - target).abs() < *d)
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
}

#[cfg(test)]
mod precursor_filter_tests {
    use super::*;
    use model::activation::ActivationMethod;
    use model::instrument::InstrumentType;
    use crate::param_model::{Param, PrecursorOffsetFrequency, SpecDataType};
    use model::protocol::Protocol;
    use model::tolerance::Tolerance;
    use std::collections::HashMap;

    /// Build a Param with a single precursor offset entry: charge 2,
    /// reduced_charge 2, offset 0.0 Da (the precursor itself), tolerance 0.5 Da.
    fn param_with_precursor_filter() -> Param {
        let mut precursor_off_map: HashMap<i32, Vec<PrecursorOffsetFrequency>> = HashMap::new();
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
            frag_off_table: HashMap::new(),
            max_rank: 3,
            rank_dist_table: HashMap::new(),
            error_scaling_factor: 0,
            ion_err_dist_table: HashMap::new(),
            noise_err_dist_table: HashMap::new(),
            ion_existence_table: HashMap::new(),
            partition_ion_types_cache: HashMap::new(),
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
        }
    }

    /// Verify the filter_mz formula for reduced_charge=2, offset=0:
    /// neutral_mass = (500.0 - PROTON) * 2 = 997.985450...
    /// c = 2 - 2 = 0 → filtered (c <= 0), so no filtering happens.
    ///
    /// Re-check: the task says "charge 2, reduced_charge 2" for the
    /// precursor itself. With c = charge - reduced_charge = 0, that
    /// would be division by zero. Real param files use reduced_charge < charge.
    ///
    /// Let's use reduced_charge=0 for the precursor filter test:
    /// c = 2 - 0 = 2; filter_mz = (neutral + 2*PROTON) / 2 + 0 = precursor_mz.
    fn param_with_precursor_filter_rc0() -> Param {
        let mut precursor_off_map: HashMap<i32, Vec<PrecursorOffsetFrequency>> = HashMap::new();
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
            frag_off_table: HashMap::new(),
            max_rank: 3,
            rank_dist_table: HashMap::new(),
            error_scaling_factor: 0,
            ion_err_dist_table: HashMap::new(),
            noise_err_dist_table: HashMap::new(),
            ion_existence_table: HashMap::new(),
            partition_ion_types_cache: HashMap::new(),
        }
    }

    #[test]
    fn precursor_peak_is_filtered_out() {
        // precursor m/z = 500.0, charge 2, reduced_charge=0:
        // c = 2 - 0 = 2
        // neutral_mass = (500.0 - PROTON) * 2 ≈ 997.9855 Da
        // filter_mz = (997.9855 + 2 * PROTON) / 2 + 0.0 = 500.0 (the precursor m/z)
        //
        // A peak AT 500.0 (the precursor m/z itself, very high intensity) should be filtered.
        // Peaks must be sorted ascending by m/z (MGF reader invariant).
        let s = make_spec(500.0, &[(100.0, 1.0), (300.0, 5.0), (500.0, 100.0)]);
        let param = param_with_precursor_filter_rc0();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new(&s, &scorer, 2);

        // The precursor peak (500.0) should be filtered out (rank u32::MAX, not returned).
        assert!(
            ss.nearest_peak_rank(500.0, 0.1).is_none(),
            "precursor peak at 500.0 should be filtered, but a peak at that m/z was found"
        );

        // The other peaks should still be present and ranked.
        // (300.0, 5.0) is now rank 1 (highest among non-filtered);
        // (100.0, 1.0) is rank 2.
        assert_eq!(ss.nearest_peak_rank(300.0, 0.1), Some(1));
        assert_eq!(ss.nearest_peak_rank(100.0, 0.1), Some(2));
    }

    #[test]
    fn non_precursor_peaks_kept() {
        // Without filtering hitting any peak, all peaks should be present.
        // The filter is at precursor m/z = 500.0 ± 0.5, no peak in this set is there.
        let s = make_spec(500.0, &[(100.0, 1.0), (200.0, 50.0), (300.0, 5.0)]);
        let param = param_with_precursor_filter_rc0();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new(&s, &scorer, 2);

        assert_eq!(ss.peak_count_after_filtering(), 3);
        assert_eq!(ss.nearest_peak_rank(200.0, 0.1), Some(1));
    }

    #[test]
    fn missing_precursor_off_map_falls_back_to_unfiltered() {
        // If param has no precursor offsets for this charge, all peaks
        // are kept and ranked normally.
        let mut param = param_with_precursor_filter_rc0();
        param.precursor_off_map.clear();
        let s = make_spec(500.0, &[(100.0, 1.0), (500.0, 100.0)]);
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new(&s, &scorer, 2);
        assert_eq!(ss.peak_count_after_filtering(), 2);
    }

    #[test]
    fn invalid_reduced_charge_skipped() {
        // reduced_charge >= charge → c = 0 → skip (no div-by-zero).
        // Using param_with_precursor_filter which has reduced_charge=2, charge=2.
        let param = param_with_precursor_filter();
        let s = make_spec(500.0, &[(100.0, 1.0), (500.0, 100.0)]);
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new(&s, &scorer, 2);
        // No filtering occurred (c <= 0 was skipped) → both peaks kept.
        assert_eq!(ss.peak_count_after_filtering(), 2);
    }
}
