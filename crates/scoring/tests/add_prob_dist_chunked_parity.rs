//! Verify chunked `add_prob_dist` is bit-identical to scalar across 10 random inputs.
//!
//! Background: Task 5 of the suffix-array refactor splits the `add_prob_dist`
//! inner loop into 4-wide chunks so LLVM can auto-vectorize on AVX2 / NEON.
//! Each destination index is unique (no cross-lane sum), so the chunked form
//! must produce IDENTICAL float bits to the scalar form. This test asserts
//! that property across 10 randomized fixtures spanning the parameter shapes
//! that appear in the production DP (various sizes, score_diff offsets,
//! aa_prob values, and pre-existing `self` contents).
//!
//! If you remove the scalar variant from the production crate, port its
//! reference body into this test file — the test is the only consumer.

use scoring::gf::score_dist::ScoreDist;

/// Reference scalar implementation, frozen here so the parity test outlives
/// the deletion of the scalar variant from the production crate. Mirrors the
/// pre-Task-5 body of `ScoreDist::add_prob_dist`.
fn add_prob_dist_scalar(
    dst: &mut ScoreDist,
    src: &ScoreDist,
    score_diff: i32,
    aa_prob: f64,
) {
    let other_min = src.min_score();
    let other_max = src.max_score();
    let self_min = dst.min_score();
    let self_max = dst.max_score();
    let t_start = other_min.max(self_min - score_diff);
    let t_end = other_max.min(self_max - score_diff);
    for t in t_start..t_end {
        let src_idx = (t - other_min) as usize;
        let dst_idx = (t + score_diff - self_min) as usize;
        let cur = dst.get_probability((t + score_diff));
        dst.set_prob((t + score_diff), cur + src_p(src, src_idx) * aa_prob);
        let _ = dst_idx; // silence
    }
}

fn src_p(d: &ScoreDist, idx: usize) -> f64 {
    // The only way to read by raw idx without exposing internals is via
    // get_probability(min + idx).
    d.get_probability(d.min_score() + idx as i32)
}

#[test]
fn chunked_matches_scalar_bit_for_bit() {
    // xorshift64* — deterministic; 10 iterations is plenty given each
    // covers an independent (size, offset, prob, contents) sample.
    let mut state: u64 = 0x1234_5678_90AB_CDEF;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };

    for iter in 0..10 {
        // Random sizes: pick self/other ranges in [4, 200) to exercise both
        // sub-chunk and multi-chunk paths (4-wide split: chunks + remainder).
        let self_len = 4 + (next() % 200) as i32;
        let other_len = 4 + (next() % 200) as i32;
        // Random min anchors in [-50, 50) so score_diff sweeps both signs.
        let self_min = -50 + (next() % 100) as i32;
        let other_min = -50 + (next() % 100) as i32;
        let self_max = self_min + self_len;
        let other_max = other_min + other_len;
        // score_diff: any int in [-150, 150) — sometimes makes t_start > t_end
        // (no-op), sometimes makes overlap partial, sometimes full.
        let score_diff = -150 + (next() % 300) as i32;
        // aa_prob: a non-trivial multiplier in [0, 1).
        let aa_prob = (next() as f64 / u64::MAX as f64).clamp(0.0, 1.0);

        // Two identical self distributions: scalar baseline + chunked target.
        let mut self_a = ScoreDist::new(self_min, self_max, false, true);
        let mut self_b = ScoreDist::new(self_min, self_max, false, true);
        // Pre-fill self with random contents so we test += (not just =).
        for i in 0..self_len {
            let v = (next() as f64 / u64::MAX as f64) * 1e-3;
            self_a.set_prob(self_min + i, v);
            self_b.set_prob(self_min + i, v);
        }
        // src: random contents.
        let mut src = ScoreDist::new(other_min, other_max, false, true);
        for i in 0..other_len {
            let v = (next() as f64 / u64::MAX as f64) * 1e-3;
            src.set_prob(other_min + i, v);
        }

        // Apply scalar reference.
        add_prob_dist_scalar(&mut self_a, &src, score_diff, aa_prob);
        // Apply production (chunked) variant.
        self_b.add_prob_dist(&src, score_diff, aa_prob);

        // Bit-identity check across the full self range.
        for i in 0..self_len {
            let s = self_min + i;
            let a = self_a.get_probability(s);
            let b = self_b.get_probability(s);
            assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "iter {} idx {}: scalar={:?} chunked={:?} \
                 (self_len={}, other_len={}, self_min={}, other_min={}, \
                 score_diff={}, aa_prob={})",
                iter, i, a, b,
                self_len, other_len, self_min, other_min, score_diff, aa_prob,
            );
        }
    }
}
