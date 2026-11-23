use rustc_hash::FxHashMap;
use scoring_crate::params::{IonType, Partition};

#[inline]
fn ensure_len(vec: &mut Vec<u64>, idx: usize) {
    if idx >= vec.len() {
        vec.resize(idx + 1, 0);
    }
}

fn vec_add(dst: &mut Vec<u64>, src: &[u64]) {
    if src.len() > dst.len() {
        dst.resize(src.len(), 0);
    }
    for (d, &s) in dst.iter_mut().zip(src.iter()) {
        *d = d.saturating_add(s);
    }
}

fn vec_sub(dst: &mut Vec<u64>, src: &[u64]) {
    if src.len() > dst.len() {
        dst.resize(src.len(), 0);
    }
    for (d, &s) in dst.iter_mut().zip(src.iter()) {
        *d = d.saturating_sub(s);
    }
}

#[inline]
fn scale_count(c: u64, factor: f32) -> u64 {
    (c as f64 * f64::from(factor)).round() as u64
}

fn vec_scale(src: &[u64], factor: f32) -> Vec<u64> {
    src.iter().map(|&c| scale_count(c, factor)).collect()
}

fn trim_trailing_zeros(vec: &mut Vec<u64>) {
    while vec.last() == Some(&0) {
        vec.pop();
    }
}

#[derive(Clone, Default, Debug)]
pub struct CountStats {
    pub rank: FxHashMap<(Partition, IonType), Vec<u64>>,
    pub error: FxHashMap<Partition, Vec<u64>>,
    pub noise_error: FxHashMap<Partition, Vec<u64>>,
    pub existence: FxHashMap<(Partition, u32), u64>,
    pub charge: FxHashMap<i32, u64>,
}

fn map_vec_eq<K>(a: &FxHashMap<K, Vec<u64>>, b: &FxHashMap<K, Vec<u64>>) -> bool
where
    K: Eq + std::hash::Hash,
{
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
    for (k, bv) in b {
        if !a.contains_key(k) && !logical_slice(bv).is_empty() {
            return false;
        }
    }
    true
}

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

impl CountStats {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn bump_rank(&mut self, p: Partition, ion: IonType, rank: u32) {
        let vec = self.rank.entry((p, ion)).or_default();
        let idx = rank as usize;
        ensure_len(vec, idx);
        vec[idx] += 1;
    }

    pub fn bump_error(&mut self, p: Partition, bin: u32) {
        let vec = self.error.entry(p).or_default();
        let idx = bin as usize;
        ensure_len(vec, idx);
        vec[idx] += 1;
    }

    pub fn bump_noise_error(&mut self, p: Partition, bin: u32) {
        let vec = self.noise_error.entry(p).or_default();
        let idx = bin as usize;
        ensure_len(vec, idx);
        vec[idx] += 1;
    }

    pub fn bump_existence(&mut self, p: Partition, idx: u32) {
        *self.existence.entry((p, idx)).or_default() += 1;
    }

    pub fn bump_charge(&mut self, charge: i32) {
        *self.charge.entry(charge).or_default() += 1;
    }

    pub fn rank_count(&self, p: &Partition, ion: IonType, rank: u32) -> u64 {
        self.rank
            .get(&(*p, ion))
            .and_then(|v| v.get(rank as usize).copied())
            .unwrap_or(0)
    }

    pub fn add(&mut self, other: &Self) {
        for (k, src) in &other.rank {
            vec_add(self.rank.entry(*k).or_default(), src);
        }
        for (k, src) in &other.error {
            vec_add(self.error.entry(*k).or_default(), src);
        }
        for (k, src) in &other.noise_error {
            vec_add(self.noise_error.entry(*k).or_default(), src);
        }
        for (&k, &v) in &other.existence {
            *self.existence.entry(k).or_default() = self
                .existence
                .get(&k)
                .copied()
                .unwrap_or(0)
                .saturating_add(v);
        }
        for (&k, &v) in &other.charge {
            *self.charge.entry(k).or_default() =
                self.charge.get(&k).copied().unwrap_or(0).saturating_add(v);
        }
    }

    pub fn sub(&mut self, other: &Self) {
        for (k, src) in &other.rank {
            if let Some(dst) = self.rank.get_mut(k) {
                vec_sub(dst, src);
                trim_trailing_zeros(dst);
            }
        }
        for (k, src) in &other.error {
            if let Some(dst) = self.error.get_mut(k) {
                vec_sub(dst, src);
                trim_trailing_zeros(dst);
            }
        }
        for (k, src) in &other.noise_error {
            if let Some(dst) = self.noise_error.get_mut(k) {
                vec_sub(dst, src);
                trim_trailing_zeros(dst);
            }
        }
        for (&k, &v) in &other.existence {
            let entry = self.existence.entry(k).or_default();
            *entry = entry.saturating_sub(v);
        }
        for (&k, &v) in &other.charge {
            let entry = self.charge.entry(k).or_default();
            *entry = entry.saturating_sub(v);
        }
    }

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
            .map(|(&k, &v)| (k, scale_count(v, factor)))
            .collect();
        let charge = self
            .charge
            .iter()
            .map(|(&k, &v)| (k, scale_count(v, factor)))
            .collect();
        Self {
            rank,
            error,
            noise_error,
            existence,
            charge,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_partition() -> Partition {
        Partition {
            charge: 2,
            parent_mass: 1000.0_f32,
            seg_num: 0,
        }
    }

    fn test_prefix_ion() -> IonType {
        IonType::Prefix {
            charge: 1,
            offset_bits: 1.007_f32.to_bits(),
            loss_class: 0,
        }
    }

    #[test]
    fn counts_add_then_sub_is_identity() {
        let p = test_partition();
        let ion = test_prefix_ion();
        let mut a = CountStats::new();
        a.bump_rank(p, ion, 3);
        let b = a.clone();
        a.add(&b);
        a.sub(&b);
        assert_eq!(a, b);
    }

    #[test]
    fn scale_halves_counts() {
        let p = test_partition();
        let ion = test_prefix_ion();
        let mut a = CountStats::new();
        a.bump_rank(p, ion, 3);
        a.bump_rank(p, ion, 3);
        let scaled = a.scaled(0.5);
        assert_eq!(scaled.rank_count(&p, ion, 3), 1);
    }

    #[test]
    fn sub_saturates_at_zero() {
        let p = test_partition();
        let ion = test_prefix_ion();
        let mut a = CountStats::new();
        a.bump_rank(p, ion, 3);

        let mut big = CountStats::new();
        big.bump_rank(p, ion, 3);
        big.bump_rank(p, ion, 3);

        a.sub(&big);
        assert_eq!(a.rank_count(&p, ion, 3), 0);
    }

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
        a.bump_rank(p, ion, 0);

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
        assert_eq!(
            z.error
                .get(&p)
                .map(|v| v.get(3).copied().unwrap_or(0))
                .unwrap_or(0),
            0
        );
        assert_eq!(z.existence.get(&(p, 0)).copied().unwrap_or(0), 0);
        assert_eq!(z.charge.get(&3).copied().unwrap_or(0), 0);
    }

    #[test]
    fn add_is_commutative_for_disjoint_keys() {
        let p1 = Partition {
            charge: 2,
            parent_mass: 500.0,
            seg_num: 0,
        };
        let p2 = Partition {
            charge: 3,
            parent_mass: 800.0,
            seg_num: 1,
        };
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

    #[test]
    fn scale_by_one_is_exact_above_f32_mantissa() {
        let big = (1u64 << 24) + 1;
        let p = test_partition();
        let ion = test_prefix_ion();
        let mut a = CountStats::new();
        a.rank.insert((p, ion), vec![0, big]);
        a.charge.insert(2, big);
        let s = a.scaled(1.0);
        assert_eq!(s, a);
        assert_eq!(s.rank_count(&p, ion, 1), big);
        assert_eq!(s.charge[&2], big);
        assert_eq!(a.scaled(0.5).charge[&2], big / 2 + 1);
    }
}
