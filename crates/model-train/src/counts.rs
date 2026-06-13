//! Sufficient statistics (raw histogram counts) accumulated from confident PSMs.
//!
//! [`CountStats`] accumulates all the raw counts needed to estimate a scoring
//! model ([`scoring_crate::param_model::Param`]).  It supports exact `add`,
//! `sub` (saturating), and `scale` so that per-source statistics can be
//! combined or subtracted incrementally.

use rustc_hash::FxHashMap;
use scoring_crate::param_model::{IonType, Partition};

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Grow `vec` in-place so `vec[idx]` exists, padding with zeros.
#[inline]
fn ensure_len(vec: &mut Vec<u64>, idx: usize) {
    if idx >= vec.len() {
        vec.resize(idx + 1, 0);
    }
}

/// Add `src` into `dst` element-wise, growing `dst` as needed.
fn vec_add(dst: &mut Vec<u64>, src: &[u64]) {
    if src.len() > dst.len() {
        dst.resize(src.len(), 0);
    }
    for (d, &s) in dst.iter_mut().zip(src.iter()) {
        *d = d.saturating_add(s);
    }
}

/// Saturating-subtract `src` from `dst` element-wise, growing `dst` as needed.
fn vec_sub(dst: &mut Vec<u64>, src: &[u64]) {
    if src.len() > dst.len() {
        dst.resize(src.len(), 0);
    }
    for (d, &s) in dst.iter_mut().zip(src.iter()) {
        *d = d.saturating_sub(s);
    }
}

/// Scale every element of `src` by `factor` (round half-away from zero) into
/// a new `Vec<u64>`.
fn vec_scale(src: &[u64], factor: f32) -> Vec<u64> {
    src.iter().map(|&c| (c as f32 * factor).round() as u64).collect()
}

/// Trim trailing zeros from `vec` so that `derive(PartialEq)` works correctly.
fn trim_trailing_zeros(vec: &mut Vec<u64>) {
    while vec.last() == Some(&0) {
        vec.pop();
    }
}

// ---------------------------------------------------------------------------
// Public type
// ---------------------------------------------------------------------------

/// Raw histogram counts (sufficient statistics) accumulated from confident PSMs.
///
/// All five maps are updated together by `add`/`sub`/`scaled`.  Two
/// `CountStats` are equal if and only if every non-zero entry in each map
/// matches; trailing zeros in Vec histograms do not affect equality because
/// every mutating method trims them.
#[derive(Clone, Default, Debug)]
pub struct CountStats {
    /// Per-rank histogram: `rank[( Partition, IonType )][rank_index]` = count.
    pub rank: FxHashMap<(Partition, IonType), Vec<u64>>,
    /// Mass-error histogram: `error[Partition][bin_index]` = count.
    pub error: FxHashMap<Partition, Vec<u64>>,
    /// Noise mass-error histogram (same structure as `error`).
    pub noise_error: FxHashMap<Partition, Vec<u64>>,
    /// Ion-existence counts: `existence[( Partition, idx )]` = count.
    pub existence: FxHashMap<(Partition, u32), u64>,
    /// Charge histogram: `charge[charge_value]` = count.
    pub charge: FxHashMap<i32, u64>,
}

// ---------------------------------------------------------------------------
// Manual PartialEq — ignores trailing zeros in all Vec histograms
// ---------------------------------------------------------------------------

/// Compare two `FxHashMap<K, Vec<u64>>` treating trailing zeros as
/// insignificant (i.e. `[1, 0]` == `[1]`).
fn map_vec_eq<K>(a: &FxHashMap<K, Vec<u64>>, b: &FxHashMap<K, Vec<u64>>) -> bool
where
    K: Eq + std::hash::Hash,
{
    // All keys present in `a` must match in `b`.
    for (k, av) in a {
        let av = logical_slice(av);
        match b.get(k) {
            Some(bv) => {
                if av != logical_slice(bv) {
                    return false;
                }
            }
            None => {
                if !av.is_empty() {
                    return false;
                }
            }
        }
    }
    // Keys present in `b` but not `a` must be logically zero.
    for (k, bv) in b {
        if !a.contains_key(k) && !logical_slice(bv).is_empty() {
            return false;
        }
    }
    true
}

/// Slice of `vec` without trailing zeros.
#[inline]
fn logical_slice(vec: &[u64]) -> &[u64] {
    let end = vec.iter().rposition(|&x| x != 0).map_or(0, |i| i + 1);
    &vec[..end]
}

impl PartialEq for CountStats {
    fn eq(&self, other: &Self) -> bool {
        map_vec_eq(&self.rank, &other.rank)
            && map_vec_eq(&self.error, &other.error)
            && map_vec_eq(&self.noise_error, &other.noise_error)
            && self.existence == other.existence
            && self.charge == other.charge
    }
}

impl Eq for CountStats {}

// ---------------------------------------------------------------------------
// Core API
// ---------------------------------------------------------------------------

impl CountStats {
    /// Create an empty `CountStats` (equivalent to `Default::default()`).
    pub fn new() -> Self {
        Self::default()
    }

    // --- bump helpers -------------------------------------------------------

    /// Increment the count for `(partition, ion, rank)` by 1.
    pub fn bump_rank(&mut self, p: Partition, ion: IonType, rank: u32) {
        let vec = self.rank.entry((p, ion)).or_default();
        let idx = rank as usize;
        ensure_len(vec, idx);
        vec[idx] += 1;
    }

    /// Increment the mass-error histogram count for `(partition, bin)` by 1.
    pub fn bump_error(&mut self, p: Partition, bin: u32) {
        let vec = self.error.entry(p).or_default();
        let idx = bin as usize;
        ensure_len(vec, idx);
        vec[idx] += 1;
    }

    /// Increment the noise mass-error histogram count for `(partition, bin)` by 1.
    pub fn bump_noise_error(&mut self, p: Partition, bin: u32) {
        let vec = self.noise_error.entry(p).or_default();
        let idx = bin as usize;
        ensure_len(vec, idx);
        vec[idx] += 1;
    }

    /// Increment the ion-existence count for `(partition, idx)` by 1.
    pub fn bump_existence(&mut self, p: Partition, idx: u32) {
        *self.existence.entry((p, idx)).or_default() += 1;
    }

    /// Increment the charge-histogram count for `charge` by 1.
    pub fn bump_charge(&mut self, charge: i32) {
        *self.charge.entry(charge).or_default() += 1;
    }

    // --- read accessor ------------------------------------------------------

    /// Return the count at `(partition, ion, rank)`, or 0 if absent.
    pub fn rank_count(&self, p: &Partition, ion: IonType, rank: u32) -> u64 {
        self.rank
            .get(&(*p, ion))
            .and_then(|v| v.get(rank as usize).copied())
            .unwrap_or(0)
    }

    // --- arithmetic ---------------------------------------------------------

    /// Add all counts from `other` into `self` element-wise.
    pub fn add(&mut self, other: &Self) {
        // rank
        for (k, src) in &other.rank {
            vec_add(self.rank.entry(*k).or_default(), src);
        }
        // error
        for (k, src) in &other.error {
            vec_add(self.error.entry(*k).or_default(), src);
        }
        // noise_error
        for (k, src) in &other.noise_error {
            vec_add(self.noise_error.entry(*k).or_default(), src);
        }
        // existence
        for (&k, &v) in &other.existence {
            *self.existence.entry(k).or_default() = self
                .existence
                .get(&k)
                .copied()
                .unwrap_or(0)
                .saturating_add(v);
        }
        // charge
        for (&k, &v) in &other.charge {
            *self.charge.entry(k).or_default() = self
                .charge
                .get(&k)
                .copied()
                .unwrap_or(0)
                .saturating_add(v);
        }
    }

    /// Subtract all counts of `other` from `self`, saturating at zero.
    pub fn sub(&mut self, other: &Self) {
        // rank
        for (k, src) in &other.rank {
            if let Some(dst) = self.rank.get_mut(k) {
                vec_sub(dst, src);
                trim_trailing_zeros(dst);
            }
            // If the key is absent in self, nothing to subtract.
        }
        // error
        for (k, src) in &other.error {
            if let Some(dst) = self.error.get_mut(k) {
                vec_sub(dst, src);
                trim_trailing_zeros(dst);
            }
        }
        // noise_error
        for (k, src) in &other.noise_error {
            if let Some(dst) = self.noise_error.get_mut(k) {
                vec_sub(dst, src);
                trim_trailing_zeros(dst);
            }
        }
        // existence
        for (&k, &v) in &other.existence {
            let entry = self.existence.entry(k).or_default();
            *entry = entry.saturating_sub(v);
        }
        // charge
        for (&k, &v) in &other.charge {
            let entry = self.charge.entry(k).or_default();
            *entry = entry.saturating_sub(v);
        }
    }

    /// Return a new `CountStats` with every count multiplied by `factor`
    /// (rounded to nearest `u64`).
    pub fn scaled(&self, factor: f32) -> Self {
        let rank = self
            .rank
            .iter()
            .map(|(&k, v)| (k, vec_scale(v, factor)))
            .collect();
        let error = self
            .error
            .iter()
            .map(|(&k, v)| (k, vec_scale(v, factor)))
            .collect();
        let noise_error = self
            .noise_error
            .iter()
            .map(|(&k, v)| (k, vec_scale(v, factor)))
            .collect();
        let existence = self
            .existence
            .iter()
            .map(|(&k, &v)| (k, (v as f32 * factor).round() as u64))
            .collect();
        let charge = self
            .charge
            .iter()
            .map(|(&k, &v)| (k, (v as f32 * factor).round() as u64))
            .collect();
        Self { rank, error, noise_error, existence, charge }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn test_partition() -> Partition {
        Partition { charge: 2, parent_mass: 1000.0_f32, seg_num: 0 }
    }

    fn test_prefix_ion() -> IonType {
        IonType::Prefix { charge: 1, offset_bits: 1.007_f32.to_bits(), loss_class: 0 }
    }

    #[test]
    fn counts_add_then_sub_is_identity() {
        let p = test_partition();
        let ion = test_prefix_ion();
        let mut a = CountStats::new();
        a.bump_rank(p, ion, 3);
        let b = a.clone();
        a.add(&b);   // count = 2
        a.sub(&b);   // count = 1 again
        assert_eq!(a, b);
    }

    #[test]
    fn scale_halves_counts() {
        let p = test_partition();
        let ion = test_prefix_ion();
        let mut a = CountStats::new();
        a.bump_rank(p, ion, 3);
        a.bump_rank(p, ion, 3); // count = 2
        let scaled = a.scaled(0.5);
        assert_eq!(scaled.rank_count(&p, ion, 3), 1); // 2 * 0.5 = 1
    }

    #[test]
    fn sub_saturates_at_zero() {
        let p = test_partition();
        let ion = test_prefix_ion();
        let mut a = CountStats::new();
        a.bump_rank(p, ion, 3); // 1

        let mut big = CountStats::new();
        big.bump_rank(p, ion, 3);
        big.bump_rank(p, ion, 3); // 2

        a.sub(&big); // 1 - 2 saturates to 0
        assert_eq!(a.rank_count(&p, ion, 3), 0);
    }

    // --- additional coverage ------------------------------------------------

    #[test]
    fn bump_error_and_add() {
        let p = test_partition();
        let mut a = CountStats::new();
        a.bump_error(p, 5);
        let b = a.clone();
        a.add(&b);
        let v = a.error.get(&p).unwrap();
        assert_eq!(v[5], 2);
    }

    #[test]
    fn bump_noise_error_and_sub_saturates() {
        let p = test_partition();
        let mut a = CountStats::new();
        a.bump_noise_error(p, 2);
        let mut big = CountStats::new();
        big.bump_noise_error(p, 2);
        big.bump_noise_error(p, 2);
        a.sub(&big);
        let v = a.noise_error.get(&p).unwrap();
        assert_eq!(v.get(2).copied().unwrap_or(0), 0);
    }

    #[test]
    fn bump_existence_and_scale() {
        let p = test_partition();
        let mut a = CountStats::new();
        a.bump_existence(p, 0);
        a.bump_existence(p, 0);
        let scaled = a.scaled(0.5);
        assert_eq!(*scaled.existence.get(&(p, 0)).unwrap(), 1);
    }

    #[test]
    fn bump_charge_and_add() {
        let mut a = CountStats::new();
        a.bump_charge(2);
        a.bump_charge(2);
        let b = a.clone();
        a.add(&b);
        assert_eq!(*a.charge.get(&2).unwrap(), 4);
    }

    #[test]
    fn default_equals_new() {
        let a: CountStats = Default::default();
        let b = CountStats::new();
        assert_eq!(a, b);
    }

    #[test]
    fn rank_count_absent_returns_zero() {
        let a = CountStats::new();
        let p = test_partition();
        let ion = test_prefix_ion();
        assert_eq!(a.rank_count(&p, ion, 99), 0);
    }

    #[test]
    fn trailing_zeros_do_not_break_equality() {
        let p = test_partition();
        let ion = test_prefix_ion();
        let mut a = CountStats::new();
        a.bump_rank(p, ion, 0); // vec = [1]

        // Manually create a vec with trailing zeros.
        let mut b = CountStats::new();
        b.rank.insert((p, ion), vec![1, 0, 0]);

        assert_eq!(a, b, "trailing zeros must not affect equality");
    }

    #[test]
    fn scaled_by_zero_gives_all_zeros() {
        let p = test_partition();
        let ion = test_prefix_ion();
        let mut a = CountStats::new();
        a.bump_rank(p, ion, 1);
        a.bump_error(p, 3);
        a.bump_existence(p, 0);
        a.bump_charge(3);
        let z = a.scaled(0.0);
        assert_eq!(z.rank_count(&p, ion, 1), 0);
        assert_eq!(z.error.get(&p).map(|v| v.get(3).copied().unwrap_or(0)).unwrap_or(0), 0);
        assert_eq!(z.existence.get(&(p, 0)).copied().unwrap_or(0), 0);
        assert_eq!(z.charge.get(&3).copied().unwrap_or(0), 0);
    }

    #[test]
    fn add_is_commutative_for_disjoint_keys() {
        let p1 = Partition { charge: 2, parent_mass: 500.0, seg_num: 0 };
        let p2 = Partition { charge: 3, parent_mass: 800.0, seg_num: 1 };
        let ion = test_prefix_ion();

        let mut a = CountStats::new();
        a.bump_rank(p1, ion, 0);

        let mut b = CountStats::new();
        b.bump_rank(p2, ion, 0);

        let mut ab = a.clone();
        ab.add(&b);

        let mut ba = b.clone();
        ba.add(&a);

        assert_eq!(ab, ba);
    }
}
