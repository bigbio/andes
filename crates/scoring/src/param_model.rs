//! The `Param` scoring model: the in-memory partition geometry, ion-type
//! layout, and learned scoring tables that the engine scores against.
//! Loaded from the canonical Parquet model store (`model-train::ModelStore`).

use std::cmp::Ordering;
use std::hash::{Hash, Hasher};
use std::sync::Arc;

use rustc_hash::FxHashMap;

use crate::gbdt_eval::GbdtPeakModel;
use model::activation::ActivationMethod;
use model::enzyme::Enzyme;
use model::instrument::InstrumentType;
use model::protocol::Protocol;
use model::tolerance::Tolerance;

#[derive(Debug, Clone, PartialEq)]
pub struct Param {
    pub version: i32,
    pub data_type: SpecDataType,
    pub mme: Tolerance,
    pub apply_deconvolution: bool,
    pub deconvolution_error_tolerance: f32,
    pub charge_hist: Vec<(i32, i32)>,
    pub min_charge: i32,
    pub max_charge: i32,
    pub num_segments: i32,
    pub partitions: Vec<Partition>,
    pub num_precursor_off: i32,
    pub precursor_off_map: FxHashMap<i32, Vec<PrecursorOffsetFrequency>>,
    pub frag_off_table: FxHashMap<Partition, Vec<FragmentOffsetFrequency>>,
    pub max_rank: i32,
    pub rank_dist_table: FxHashMap<Partition, FxHashMap<IonType, Vec<f32>>>,
    pub error_scaling_factor: i32,
    pub ion_err_dist_table: FxHashMap<Partition, Vec<f32>>,
    pub noise_err_dist_table: FxHashMap<Partition, Vec<f32>>,
    pub ion_existence_table: FxHashMap<Partition, Vec<f32>>,
    /// Pre-filtered ion-type list per partition (Noise excluded), populated
    /// at load time. Used by `ion_types_for_partition_slice` to avoid
    /// per-call Vec allocation in the node-scoring DP hot path.
    /// Call `rebuild_cache()` after manually constructing a `Param` in tests
    /// or any context where the cache was not populated at construction time.
    pub partition_ion_types_cache: FxHashMap<Partition, Vec<IonType>>,
    /// Optional peptide-agnostic GBDT per-peak signal/noise model. Populated by
    /// the store reader from the manifest row's `gbdt_model_bytes` blob; `None`
    /// for stores / slugs without a trained GBDT (scoring is
    /// then byte-identical to the pre-GBDT engine).
    pub gbdt_peak_model: Option<GbdtPeakModel>,
    /// Optional GBDT fragment-intensity model (regressor, raw `predict_value`
    /// output). Populated from the `frag_intensity_model_bytes` manifest column;
    /// `None` for any store written before this column existed or for slugs
    /// without a trained intensity regressor.
    ///
    /// Wrapped in `Arc` so callers can share one model across parallel search
    /// threads without cloning the tree arrays.
    pub frag_intensity_model: Option<Arc<GbdtPeakModel>>,
    /// Optional GBDT rich-ion LLR classifier (logistic, raw `predict_value`
    /// output). Populated from the `rich_ion_model_bytes` manifest column;
    /// `None` for any store written before this column existed or for slugs
    /// without a trained rich-ion classifier.
    ///
    /// Wrapped in `Arc` so callers can share one model across parallel search
    /// threads without cloning the tree arrays.
    pub rich_ion_model: Option<Arc<GbdtPeakModel>>,
}

/// Build the per-partition ion-type cache (Noise excluded). Single source of
/// truth for cache construction (the Parquet store reader and
/// `Param::rebuild_cache`).
fn build_partition_ion_types_cache(
    frag_off_table: &FxHashMap<Partition, Vec<FragmentOffsetFrequency>>,
) -> FxHashMap<Partition, Vec<IonType>> {
    let mut cache: FxHashMap<Partition, Vec<IonType>> =
        FxHashMap::with_capacity_and_hasher(frag_off_table.len(), Default::default());
    for (&part, frag_list) in frag_off_table {
        let mut ions: Vec<IonType> = Vec::with_capacity(frag_list.len());
        for fof in frag_list {
            if !matches!(fof.ion_type, IonType::Noise) {
                ions.push(fof.ion_type);
            }
        }
        cache.insert(part, ions);
    }
    cache
}

impl Param {
    /// Find the partition matching `(charge, parent_mass, seg_num)` via a
    /// floor lookup (the largest partition ≤ target by lex order on
    /// `(charge, parent_mass.to_bits(), seg_num)`).
    ///
    /// Falls back gracefully:
    /// - If no partition matches the requested charge: use the smallest
    ///   charge available with the requested mass + segment.
    /// - If charge > all available: use the largest available charge.
    pub fn find_partition(&self, charge: i32, parent_mass: f32, seg_num: i32) -> Option<Partition> {
        if self.partitions.is_empty() {
            return None;
        }

        // Build the target partition for the floor lookup.
        let target = Partition { charge, parent_mass, seg_num };

        // partitions is already sorted (loader invariant). Find the largest
        // partition <= target via binary search.
        let pos = self.partitions.partition_point(|p| p <= &target);
        if pos > 0 {
            // partitions[pos - 1] is the largest <= target.
            let candidate = self.partitions[pos - 1];
            if candidate.charge == charge {
                return Some(candidate);
            }
            // Floor returned a partition with smaller charge: if no
            // exact-charge match, find smallest available charge, then floor
            // on (smallest_charge, parent_mass, seg_num).
        }

        // Fall back: find smallest charge in partitions, retry.
        let min_charge = self.partitions.iter().map(|p| p.charge).min()?;
        let max_charge = self.partitions.iter().map(|p| p.charge).max()?;
        let fallback_charge = if charge < min_charge {
            min_charge
        } else if charge > max_charge {
            max_charge
        } else {
            // charge is in range but had no exact partition (a gap between
            // available charges). Prefer the floor candidate — the nearest
            // partition <= target — over jumping to `partitions.last()` (max
            // charge / mass / segment), which would score against an unrelated
            // partition. If the target sorts below every partition (pos == 0,
            // e.g. parent_mass below the smallest trained mass) there is no
            // floor, so fall through to the clamp-and-retry path with min_charge.
            let pos = self.partitions.partition_point(|p| p <= &target);
            if pos > 0 {
                return self.partitions.get(pos - 1).copied();
            }
            min_charge
        };
        let fallback_target = Partition { charge: fallback_charge, parent_mass, seg_num };
        let fallback_pos = self.partitions.partition_point(|p| p <= &fallback_target);
        if fallback_pos > 0 {
            let candidate = self.partitions[fallback_pos - 1];
            if candidate.charge == fallback_charge {
                return Some(candidate);
            }
        }
        // Last resort: just return any partition with the fallback charge.
        self.partitions.iter().find(|p| p.charge == fallback_charge).copied()
    }

    /// Compute the segment number for a peak m/z relative to the peptide's
    /// parent mass.
    pub fn segment_num_for(&self, peak_mz: f64, parent_mass: f64) -> i32 {
        if parent_mass <= 0.0 || self.num_segments <= 0 {
            return 0;
        }
        let seg = (peak_mz / parent_mass * self.num_segments as f64) as i32;
        seg.min(self.num_segments - 1).max(0)
    }

    /// Alias for `segment_num_for` matching the name used by the node-scoring DP code
    /// (`param.segment_num(theo_mz, parent_mass)`).
    #[inline]
    pub fn segment_num(&self, peak_mz: f64, parent_mass: f64) -> usize {
        self.segment_num_for(peak_mz, parent_mass) as usize
    }

    /// Collect the unique ion types (Prefix and Suffix, not Noise) whose
    /// partition has `seg_num == seg`. Derived from `frag_off_table` keys
    /// (ion-type membership lives in `frag_off_table`, not `rank_dist_table`).
    ///
    /// Returned in stable insertion order; duplicates suppressed.
    pub fn ion_types_for_segment(&self, seg: usize) -> Vec<IonType> {
        let mut seen: std::collections::HashSet<IonType> = std::collections::HashSet::new();
        let mut out: Vec<IonType> = Vec::new();
        for (partition, frag_list) in &self.frag_off_table {
            if partition.seg_num as usize != seg {
                continue;
            }
            for fof in frag_list {
                let ion = fof.ion_type;
                if matches!(ion, IonType::Noise) {
                    continue;
                }
                if seen.insert(ion) {
                    out.push(ion);
                }
            }
        }
        out
    }

    /// Find the partition for `(charge, parent_mass, seg_num)` using the
    /// floor-lookup semantics of `find_partition`. Returns a synthetic
    /// partition if none is found (so callers don't need to unwrap).
    pub fn partition_for(&self, charge: u8, parent_mass: f64, seg_num: usize) -> Partition {
        self.find_partition(charge as i32, parent_mass as f32, seg_num as i32)
            .unwrap_or(Partition {
                charge: charge as i32,
                parent_mass: parent_mass as f32,
                seg_num: seg_num as i32,
            })
    }

    /// Ion types for the SPECIFIC partition `(charge, parent_mass, seg)`.
    ///
    /// Selects the partition's ion list from `frag_off_table` rather than
    /// the segment-wide union returned by `ion_types_for_segment`. Used
    /// in the per-node scoring path.
    pub fn ion_types_for_partition(&self, charge: u8, parent_mass: f64, seg: usize) -> Vec<IonType> {
        // Compat shim — callers in hot paths should use
        // `ion_types_for_partition_slice` to avoid the allocation.
        self.ion_types_for_partition_slice(charge, parent_mass, seg).to_vec()
    }

    /// Slice-borrowing version of `ion_types_for_partition`. Reads from the
    /// pre-filtered `partition_ion_types_cache` populated at param-load time.
    /// Zero allocations per call. Used by the node-scoring DP hot path.
    pub fn ion_types_for_partition_slice(&self, charge: u8, parent_mass: f64, seg: usize) -> &[IonType] {
        let part = self.partition_for(charge, parent_mass, seg);
        self.partition_ion_types_cache
            .get(&part)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }


    /// Rebuild the `partition_ion_types_cache` from `frag_off_table`.
    /// Call this after manually constructing a `Param` in tests or any
    /// context where the cache was not populated at construction time.
    /// The Parquet store reader builds the cache automatically.
    pub fn rebuild_cache(&mut self) {
        self.partition_ion_types_cache = build_partition_ion_types_cache(&self.frag_off_table);
    }
}


#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct SpecDataType {
    pub activation: ActivationMethod,
    pub instrument: InstrumentType,
    pub enzyme: Option<Enzyme>,
    pub protocol: Protocol,
}

#[derive(Debug, Clone, Copy)]
pub struct Partition {
    pub charge: i32,
    pub parent_mass: f32,
    pub seg_num: i32,
}

impl PartialEq for Partition {
    fn eq(&self, other: &Self) -> bool {
        self.charge == other.charge
            && self.parent_mass.to_bits() == other.parent_mass.to_bits()
            && self.seg_num == other.seg_num
    }
}

impl Eq for Partition {}

impl Hash for Partition {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.charge.hash(state);
        self.parent_mass.to_bits().hash(state);
        self.seg_num.hash(state);
    }
}

impl Ord for Partition {
    fn cmp(&self, other: &Self) -> Ordering {
        // Lex order: charge → seg_num → parent_mass.
        // The order is load-bearing: a charge → parent_mass → seg_num order
        // produces wrong floor-lookup results for `find_partition` (seg=0
        // queries would return a seg=1 partition with the same parent_mass
        // tier, resolving to the wrong rank distribution table).
        self.charge.cmp(&other.charge)
            .then_with(|| self.seg_num.cmp(&other.seg_num))
            .then_with(|| self.parent_mass.to_bits().cmp(&other.parent_mass.to_bits()))
    }
}

impl PartialOrd for Partition {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum IonType {
    /// `offset_bits` is `f32::to_bits` so the type can derive Eq/Hash;
    /// recover the float via `offset()`.
    /// `loss_class`: 0 = intact (no neutral loss); 1.. = per-mod-class loss pool.
    Prefix { charge: i32, offset_bits: u32, loss_class: u8 },
    Suffix { charge: i32, offset_bits: u32, loss_class: u8 },
    Noise,
}

impl IonType {
    pub fn offset(&self) -> Option<f32> {
        match self {
            IonType::Prefix { offset_bits, .. } | IonType::Suffix { offset_bits, .. } => {
                Some(f32::from_bits(*offset_bits))
            }
            IonType::Noise => None,
        }
    }

    pub fn charge(&self) -> Option<i32> {
        match self {
            IonType::Prefix { charge, .. } | IonType::Suffix { charge, .. } => Some(*charge),
            IonType::Noise => None,
        }
    }

    pub fn is_prefix(&self) -> bool { matches!(self, IonType::Prefix { .. }) }
    pub fn is_suffix(&self) -> bool { matches!(self, IonType::Suffix { .. }) }
    pub fn is_noise(&self) -> bool { matches!(self, IonType::Noise) }

    /// Loss-class id: 0 = intact; 1.. = a per-mod-class neutral-loss pool.
    pub fn loss_class(&self) -> u8 {
        match self {
            IonType::Prefix { loss_class, .. } | IonType::Suffix { loss_class, .. } => *loss_class,
            IonType::Noise => 0,
        }
    }
    /// True if this is a neutral-loss-shifted fragment ion (any loss class).
    pub fn is_loss(&self) -> bool { self.loss_class() != 0 }

    /// Compute the predicted m/z for this ion type given a **nominal** node mass.
    ///
    /// Formula:
    ///   `real_mass = node_nominal / INTEGER_MASS_SCALER`
    ///   `mz = real_mass / charge + offset`
    ///
    /// The `offset` field already includes the proton mass contribution
    /// (for b-ions: `offset = PROTON ≈ 1.00728`; for y-ions: `offset = H2O + PROTON ≈ 19.018`).
    /// The `INTEGER_MASS_SCALER` division converts integer nominal mass back to real
    /// monoisotopic mass before dividing by charge.
    ///
    /// For `Noise`, returns 0.0.
    pub fn mz(&self, node_nominal: f64) -> f64 {
        match self {
            IonType::Prefix { charge, offset_bits, .. } | IonType::Suffix { charge, offset_bits, .. } => {
                let offset = f32::from_bits(*offset_bits) as f64;
                let c = *charge as f64;
                // real_mass = node_nominal / INTEGER_MASS_SCALER
                // mz        = real_mass / charge + offset
                let real_mass = node_nominal / model::mass::INTEGER_MASS_SCALER as f64;
                real_mass / c + offset
            }
            IonType::Noise => 0.0,
        }
    }

    /// Inverse of `mz`: given an observed peak m/z, recover the real node mass (in Da).
    ///
    /// Formula: `real_mass = (mz - offset) * charge`
    ///
    /// Returns the real monoisotopic node mass (Da), NOT nominal mass.
    /// For `Noise`: returns 0.0.
    pub fn mass_from_mz(&self, mz: f64) -> f64 {
        match self {
            IonType::Prefix { charge, offset_bits, .. } | IonType::Suffix { charge, offset_bits, .. } => {
                let offset = f32::from_bits(*offset_bits) as f64;
                let c = *charge as f64;
                (mz - offset) * c
            }
            IonType::Noise => 0.0,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrecursorOffsetFrequency {
    pub reduced_charge: i32,
    pub offset: f32,
    pub tolerance: Tolerance,
    pub frequency: f32,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FragmentOffsetFrequency {
    pub ion_type: IonType,
    pub frequency: f32,
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loss_class_is_distinct_key_from_intact() {
        use std::collections::HashMap;
        let intact  = IonType::Prefix { charge: 1, offset_bits: 1.0f32.to_bits(), loss_class: 0 };
        let glyco   = IonType::Prefix { charge: 1, offset_bits: 1.0f32.to_bits(), loss_class: 1 };
        let phospho = IonType::Prefix { charge: 1, offset_bits: 1.0f32.to_bits(), loss_class: 2 };
        let mut m = HashMap::new();
        m.insert(intact, "i"); m.insert(glyco, "g"); m.insert(phospho, "p");
        assert_eq!(m.len(), 3);
        assert!(!intact.is_loss());
        assert!(glyco.is_loss() && phospho.is_loss());
        assert_eq!(glyco.loss_class(), 1);
        assert_eq!(intact.loss_class(), 0);
    }

    #[test]
    fn partition_eq_via_to_bits() {
        let a = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
        let b = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
        assert_eq!(a, b);
        let c = Partition { charge: 2, parent_mass: 1000.0001, seg_num: 0 };
        assert_ne!(a, c);
    }

    #[test]
    fn partition_ord_lex_order() {
        let a = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
        let b = Partition { charge: 2, parent_mass: 1000.0, seg_num: 1 };
        let c = Partition { charge: 3, parent_mass: 500.0,  seg_num: 0 };
        assert!(a < b);
        assert!(b < c);
    }

    #[test]
    fn partition_hash_consistent_with_eq() {
        use std::collections::HashSet;
        let a = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
        let b = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
        let set: HashSet<_> = [a, b].into_iter().collect();
        assert_eq!(set.len(), 1);
    }

    #[test]
    fn ion_type_helpers() {
        let p = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits(), loss_class: 0 };
        let s = IonType::Suffix { charge: 1, offset_bits: 0.0_f32.to_bits(), loss_class: 0 };
        let n = IonType::Noise;
        assert!(p.is_prefix());  assert!(!p.is_suffix()); assert!(!p.is_noise());
        assert!(!s.is_prefix()); assert!(s.is_suffix());  assert!(!s.is_noise());
        assert!(!n.is_prefix()); assert!(!n.is_suffix()); assert!(n.is_noise());
        assert_eq!(p.charge(), Some(1));
        assert_eq!(n.charge(), None);
    }

    #[test]
    fn ion_type_offset_round_trip() {
        let i = IonType::Prefix { charge: 2, offset_bits: 1.5_f32.to_bits(), loss_class: 0 };
        assert_eq!(i.offset(), Some(1.5));
    }


    fn make_param() -> Param {
        use model::activation::ActivationMethod;
        use model::instrument::InstrumentType;
        use model::protocol::Protocol;
        use model::tolerance::Tolerance;

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
            charge_hist: vec![],
            min_charge: 2,
            max_charge: 3,
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
    fn find_partition_exact_charge_match() {
        let mut param = make_param();
        param.partitions = vec![
            Partition { charge: 2, parent_mass: 500.0, seg_num: 0 },
            Partition { charge: 2, parent_mass: 500.0, seg_num: 1 },
            Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 },
            Partition { charge: 3, parent_mass: 500.0, seg_num: 0 },
        ];
        // Sort matches the loader invariant.
        param.partitions.sort();

        // Partition Ord: charge → seg_num → parent_mass.
        // Sorted order: (2,seg0,500), (2,seg0,1000), (2,seg1,500), (3,seg0,500).
        // Target (2, 800.0, seg0): floor is (2,seg0,500) — same charge, same seg,
        // and 500.0 < 800.0. The next candidate (2,seg0,1000) is above 800.0.
        // seg1 partitions are NOT considered because seg_num 1 > 0 = target seg.
        let p = param.find_partition(2, 800.0, 0).expect("find");
        assert_eq!(p.charge, 2);
        assert_eq!(p.parent_mass, 500.0);
        assert_eq!(p.seg_num, 0);
    }

    #[test]
    fn find_partition_low_charge_fallback() {
        let mut param = make_param();
        param.partitions = vec![
            Partition { charge: 2, parent_mass: 500.0, seg_num: 0 },
            Partition { charge: 3, parent_mass: 500.0, seg_num: 0 },
        ];
        param.partitions.sort();

        // Target charge 1 (below all): falls back to smallest charge = 2.
        let p = param.find_partition(1, 500.0, 0).expect("find with fallback");
        assert_eq!(p.charge, 2);
    }

    #[test]
    fn find_partition_high_charge_fallback() {
        let mut param = make_param();
        param.partitions = vec![
            Partition { charge: 2, parent_mass: 500.0, seg_num: 0 },
            Partition { charge: 3, parent_mass: 500.0, seg_num: 0 },
        ];
        param.partitions.sort();

        // Target charge 5 (above all): falls back to largest = 3.
        let p = param.find_partition(5, 500.0, 0).expect("find with fallback");
        assert_eq!(p.charge, 3);
    }

    #[test]
    fn segment_num_clamps_to_max() {
        let mut param = make_param();
        param.num_segments = 3;
        // peak_mz / parent_mass × num_segments = floor calculation
        assert_eq!(param.segment_num_for(50.0, 100.0), 1);
        assert_eq!(param.segment_num_for(99.0, 100.0), 2);
        assert_eq!(param.segment_num_for(100.0, 100.0), 2);  // clamped
        assert_eq!(param.segment_num_for(120.0, 100.0), 2);  // clamped
    }

    #[test]
    fn ion_type_mz_prefix_charge1_offset0() {
        // mz = (node_nominal / INTEGER_MASS_SCALER) / charge + offset
        // For Prefix(charge=1, offset=0): mz = (node_nominal / 0.999497) / 1 + 0
        use model::mass::INTEGER_MASS_SCALER;
        let ion = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits(), loss_class: 0 };
        let node_nominal = 100.0_f64;
        let expected = (node_nominal / INTEGER_MASS_SCALER as f64) / 1.0;
        assert!((ion.mz(node_nominal) - expected).abs() < 1e-9);
    }

    #[test]
    fn ion_type_mz_prefix_charge2() {
        // mz = (node_nominal / INTEGER_MASS_SCALER) / charge + offset
        // For Prefix(charge=2, offset=0): mz = (node_nominal / 0.999497) / 2
        use model::mass::INTEGER_MASS_SCALER;
        let ion = IonType::Prefix { charge: 2, offset_bits: 0.0_f32.to_bits(), loss_class: 0 };
        let node_nominal = 200.0_f64;
        let expected = (node_nominal / INTEGER_MASS_SCALER as f64) / 2.0;
        assert!((ion.mz(node_nominal) - expected).abs() < 1e-9);
    }

    #[test]
    fn ion_type_mz_prefix_with_b_ion_offset() {
        // Realistic b-ion case: offset = PROTON (≈1.00728).
        // mz = (node_nominal / INTEGER_MASS_SCALER) / charge + PROTON
        use model::mass::{PROTON, INTEGER_MASS_SCALER};
        let b_ion = IonType::Prefix { charge: 1, offset_bits: (PROTON as f32).to_bits(), loss_class: 0 };
        let node_nominal = 100.0_f64;
        let expected = (node_nominal / INTEGER_MASS_SCALER as f64) / 1.0 + PROTON;
        assert!((b_ion.mz(node_nominal) - expected).abs() < 1e-4);
    }

    #[test]
    fn ion_type_mz_suffix_same_formula_as_prefix() {
        // Suffix uses the same mz formula as prefix.
        let offset = 18.01_f32;
        let prefix = IonType::Prefix { charge: 1, offset_bits: offset.to_bits(), loss_class: 0 };
        let suffix = IonType::Suffix { charge: 1, offset_bits: offset.to_bits(), loss_class: 0 };
        let node_nominal = 150.0_f64;
        assert!((prefix.mz(node_nominal) - suffix.mz(node_nominal)).abs() < 1e-9);
    }

    #[test]
    fn ion_type_mz_noise_returns_zero() {
        assert_eq!(IonType::Noise.mz(100.0), 0.0);
    }

    #[test]
    fn ion_type_mass_from_mz_roundtrip() {
        // mass_from_mz(mz) = (mz - offset) * charge
        // Returns the REAL monoisotopic mass (Da), not nominal mass.
        // Round-trip: mz(nominal) → mass_from_mz(mz) = (nominal/scaler/c+offset - offset)*c
        //           = (nominal / scaler) = real_mass  (NOT the original nominal input).
        use model::mass::INTEGER_MASS_SCALER;
        let offset = 1.00782_f32; // realistic b-ion offset
        let ion = IonType::Prefix { charge: 1, offset_bits: offset.to_bits(), loss_class: 0 };
        let node_nominal = 100.0_f64;
        let mz = ion.mz(node_nominal);
        let recovered_real_mass = ion.mass_from_mz(mz);
        // Recovered mass should equal node_nominal / INTEGER_MASS_SCALER (real mass)
        let expected_real_mass = node_nominal / INTEGER_MASS_SCALER as f64;
        assert!((recovered_real_mass - expected_real_mass).abs() < 1e-4,
            "mass_from_mz returned {recovered_real_mass}, expected real mass {expected_real_mass}");
    }

    #[test]
    fn param_defaults_gbdt_model_to_none() {
        let p = crate::testutil::tiny_param();
        assert!(p.gbdt_peak_model.is_none(), "fresh param must carry no GBDT model");
    }

    #[test]
    fn param_defaults_frag_intensity_model_to_none() {
        let p = crate::testutil::tiny_param();
        assert!(
            p.frag_intensity_model.is_none(),
            "fresh param must carry no frag_intensity_model"
        );
    }

    #[test]
    fn frag_intensity_model_can_be_set_and_read() {
        use std::sync::Arc;
        use crate::gbdt_eval::{GbdtPeakModel, Tree};
        let model = Arc::new(GbdtPeakModel {
            n_features: 1,
            apply_sigmoid: false,
            trees: vec![Tree {
                feature: vec![-1],
                threshold: vec![0.0],
                left: vec![-1],
                right: vec![-1],
                value: vec![3.5],
                default_left: vec![1],
            }],
            iso_x: vec![],
            iso_y: vec![],
        });
        let mut p = crate::testutil::tiny_param();
        p.frag_intensity_model = Some(Arc::clone(&model));
        let v = p.frag_intensity_model.as_ref().unwrap().predict_value(&[]);
        assert!((v - 3.5).abs() < 1e-5, "expected 3.5, got {v}");
    }

    #[test]
    fn ion_types_for_segment_returns_unique() {
        use model::activation::ActivationMethod;
        use model::instrument::InstrumentType;
        use model::protocol::Protocol;
        use model::tolerance::Tolerance;

        let part = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
        let prefix = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits(), loss_class: 0 };
        let suffix = IonType::Suffix { charge: 1, offset_bits: 0.0_f32.to_bits(), loss_class: 0 };

        // Populate frag_off_table (the source of truth for ion_types_for_segment).
        let mut frag_off_table: FxHashMap<Partition, Vec<FragmentOffsetFrequency>> = FxHashMap::default();
        frag_off_table.insert(part, vec![
            FragmentOffsetFrequency { ion_type: prefix, frequency: 0.7 },
            FragmentOffsetFrequency { ion_type: suffix, frequency: 0.6 },
        ]);

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
            charge_hist: vec![],
            min_charge: 2,
            max_charge: 2,
            num_segments: 1,
            partitions: vec![part],
            num_precursor_off: 0,
            precursor_off_map: FxHashMap::default(),
            frag_off_table,
            max_rank: 2,
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
        param.rebuild_cache();

        let seg0 = param.ion_types_for_segment(0);
        // Should return prefix and suffix (not noise), no duplicates.
        assert_eq!(seg0.len(), 2);
        assert!(seg0.iter().all(|i| !i.is_noise()));
        assert!(seg0.iter().any(|i| i.is_prefix()));
        assert!(seg0.iter().any(|i| i.is_suffix()));

        // Segment 1 has no partitions → empty.
        let seg1 = param.ion_types_for_segment(1);
        assert!(seg1.is_empty());
    }
}
