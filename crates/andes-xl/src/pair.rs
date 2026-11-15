//! Candidate lookup by mass and pair enumeration by summed mass.
//!
//! Candidates are referred to by `u32` index into the caller's candidate
//! table; a pair is `(u32, u32)` with the first index not larger than the
//! second, so a homodimer `(i, i)` is one pair.

/// Candidate indices sorted by neutral mass (ties by index).
#[derive(Debug, Clone, Default)]
pub struct MassIndex {
    masses: Vec<f64>,
    order: Vec<u32>,
    sorted_masses: Vec<f64>,
}

impl MassIndex {
    /// Index the candidates whose neutral masses are `masses`.
    pub fn new(masses: Vec<f64>) -> Self {
        let mut order: Vec<u32> = (0..masses.len() as u32).collect();
        order.sort_by(|&a, &b| {
            masses[a as usize]
                .total_cmp(&masses[b as usize])
                .then(a.cmp(&b))
        });
        let sorted_masses = order.iter().map(|&i| masses[i as usize]).collect();
        MassIndex {
            masses,
            order,
            sorted_masses,
        }
    }

    /// Number of indexed candidates.
    pub fn len(&self) -> usize {
        self.masses.len()
    }

    /// Whether the index is empty.
    pub fn is_empty(&self) -> bool {
        self.masses.is_empty()
    }

    /// Neutral mass of candidate `i`.
    pub fn mass(&self, i: u32) -> f64 {
        self.masses[i as usize]
    }

    /// Candidates with mass in `[lo, hi]`, in mass order.
    pub fn in_range(&self, lo: f64, hi: f64) -> &[u32] {
        let a = self.sorted_masses.partition_point(|&m| m < lo);
        let b = self.sorted_masses.partition_point(|&m| m <= hi);
        &self.order[a..b.max(a)]
    }

    /// Candidates within `tol` of `mass`.
    pub fn within(&self, mass: f64, tol: f64) -> &[u32] {
        self.in_range(mass - tol, mass + tol)
    }
}

/// Every pair `(i, j)`, `i <= j`, with `mass(i) + mass(j)` in
/// `[target - tol, target + tol]`, sorted.
///
/// Walks the mass-sorted candidates as the lighter partner and looks the
/// heavier one up by binary search, so the cost is `O(n log n + pairs)`.
pub fn pairs_by_sum(index: &MassIndex, target: f64, tol: f64) -> Vec<(u32, u32)> {
    let mut out = Vec::new();
    let n = index.sorted_masses.len();
    for a in 0..n {
        let ma = index.sorted_masses[a];
        // The lighter member is at most half the sum.
        if ma > (target + tol) / 2.0 + 1e-9 {
            break;
        }
        // The window is widened by a rounding slack; the exact test below
        // decides membership the same way as the brute-force reference.
        let lo = target - tol - ma - 1e-9;
        let hi = target + tol - ma + 1e-9;
        let b0 = index.sorted_masses.partition_point(|&m| m < lo).max(a);
        for b in b0..n {
            let mb = index.sorted_masses[b];
            if mb > hi {
                break;
            }
            if (ma + mb - target).abs() > tol {
                continue;
            }
            let (i, j) = (index.order[a], index.order[b]);
            out.push((i.min(j), i.max(j)));
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

/// Brute-force reference for [`pairs_by_sum`]: every `i <= j` checked.
pub fn pairs_by_sum_brute(masses: &[f64], target: f64, tol: f64) -> Vec<(u32, u32)> {
    let mut out = Vec::new();
    for i in 0..masses.len() {
        for j in i..masses.len() {
            if (masses[i] + masses[j] - target).abs() <= tol {
                out.push((i as u32, j as u32));
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic pseudo-random masses (xorshift).
    fn masses(n: usize, seed: u64) -> Vec<f64> {
        let mut s = seed;
        (0..n)
            .map(|_| {
                s ^= s << 13;
                s ^= s >> 7;
                s ^= s << 17;
                500.0 + (s % 2_000_000) as f64 / 1000.0
            })
            .collect()
    }

    #[test]
    fn pair_enumeration_matches_brute_force() {
        for seed in 1..20u64 {
            let m = masses(300, seed * 7919);
            let idx = MassIndex::new(m.clone());
            for &(target, tol) in &[(2500.0, 0.5), (1500.0, 2.0), (3000.0, 0.05), (1000.0, 0.0)] {
                assert_eq!(
                    pairs_by_sum(&idx, target, tol),
                    pairs_by_sum_brute(&m, target, tol),
                    "seed {seed} target {target} tol {tol}"
                );
            }
        }
    }

    #[test]
    fn homodimers_and_duplicate_masses_are_enumerated_once() {
        let m = vec![1000.0, 1000.0, 1500.0];
        let idx = MassIndex::new(m.clone());
        let got = pairs_by_sum(&idx, 2000.0, 1e-6);
        assert_eq!(got, vec![(0, 0), (0, 1), (1, 1)]);
        assert_eq!(got, pairs_by_sum_brute(&m, 2000.0, 1e-6));
    }

    #[test]
    fn range_lookup_is_inclusive_and_mass_ordered() {
        let idx = MassIndex::new(vec![3.0, 1.0, 2.0, 2.0]);
        assert_eq!(idx.in_range(2.0, 3.0), &[2, 3, 0]);
        assert_eq!(idx.within(1.0, 0.0), &[1]);
        assert!(idx.within(5.0, 0.1).is_empty());
    }
}
