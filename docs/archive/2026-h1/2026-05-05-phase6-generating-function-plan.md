# Phase 6 — Generating-Function (SpecEValue) — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Land the graph-based generating-function DP that produces `SpecEValue` per PSM. Replaces the placeholder `score` from Phase 5 as the primary ranking statistic and enables Phase 7 (`.tsv` / `.pin` outputs that downstream tools and Percolator expect).

**Architecture:** Mirror Java's `PrimitiveAminoAcidGraph` + `PrimitiveGeneratingFunction` + `PrimitiveGeneratingFunctionGroup` (the CSR-format primitive-array variant DBScanner uses) in Rust under `engine::gf::*`. The DP traverses a per-spectrum amino-acid graph from source mass 0 to peptide-sink mass M, accumulating a per-node score distribution by combining predecessor distributions weighted by edge AA probability and shifted by combined node+edge scores. A bin-group merges per-mass-index distributions across the precursor tolerance window. SpecEValue is the cumulative tail probability `P(X >= score)` on the merged distribution.

**Tech Stack:** Rust 2021, single workspace `astral-speed/rust`. Tests via `cargo test -p engine`. No new crates.

**Branch:** `rust-implement` (continues from Phase 6 partial head `779f617`). Ship as milestone commits — single closing PR at end of the rewrite.

**Spec / reference:**
- Roadmap: `docs/superpowers/rust-full-rewrite-roadmap.md` § Phase 6.
- Java reference: `edu/ucsd/msjava/msgf/{ScoreDist,PrimitiveAminoAcidGraph,PrimitiveGeneratingFunction,PrimitiveGeneratingFunctionGroup}.java` and `edu/ucsd/msjava/msscorer/NewScoredSpectrum.java` (`getNodeScore`, `getEdgeScore`).
- Caller in Java: `edu/ucsd/msjava/msdbsearch/DBScanner.java:617-650` (the `PrimitiveGeneratingFunctionGroup` loop over `[minPeptideMassIndex, maxPeptideMassIndex]`).

---

## Status of prior tasks (already shipped)

| Task | Commit | What |
|---|---|---|
| 1 | `fcc0f05` | `gf::score_dist::{ScoreBound, ScoreDist}` (data wrappers, set/get/add/get_spectral_probability) |
| 2 | `779f617` | `gf::generating_function::GeneratingFunction::compute` — *placeholder* DP with a generic increment-score callback and uniform AA prior. Does **not** model the source/intermediate/sink graph, enzyme cleavage credit, per-AA probabilities, or per-node spectrum-derived scores. **Task 6 below replaces it.** |

Task 2's API is preserved as `compute_uniform` for tests; the real production API becomes `compute(graph)` (Task 6).

---

## Tasks at a glance

| Task | Scope | LOC est. | Tests |
|---|---|---:|---:|
| 3 | `ScoreDist::add_prob_dist` + `add_num_dist` | ~80 | 6 |
| 4 | `ScoredSpectrum::node_score` + `edge_score` (per-nominal-mass) | ~200 | 8 |
| 5 | `PrimitiveAaGraph` (CSR amino-acid graph) | ~330 | 8 |
| 6 | `GeneratingFunction::compute(graph, …)` — graph-based DP + enzyme adjustment | ~230 | 7 |
| 7 | `GeneratingFunctionGroup` (multi-bin merger) | ~110 | 5 |
| 8 | Wire SpecEValue into `PsmMatch` + `match_engine` | ~150 | 5 |
| 9 | Java parity on traced spectra (hand-picked) | ~120 | 3 |
| 10 | BSA + test.mgf SpecEValue parity gate | ~100 | 1 |

Total: ~1,320 LOC engine + ~43 tests.

---

## Task 3: `ScoreDist::add_prob_dist` + `add_num_dist`

**Files:**
- Modify: `astral-speed/rust/crates/engine/src/gf/score_dist.rs`

The DP combines predecessor distributions by adding (`other.prob[t]` × `aa_prob`) into `self.prob[t + score_diff]` for each `t` in `other`'s score range that lands in `self`'s range. Mirrors Java `ScoreDist.addProbDist(otherDist, scoreDiff, aaProb)`.

- [ ] **Step 3.1: Write failing tests**

Append to the existing `mod tests` in `score_dist.rs`:

```rust
    #[test]
    fn add_prob_dist_offset_zero_scalar_one() {
        // self range [0, 5), other range [0, 5). After add_prob_dist(other, 0, 1.0)
        // each self[s] += other[s].
        let mut a = ScoreDist::new(0, 5, false, true);
        let mut b = ScoreDist::new(0, 5, false, true);
        for s in 0..5 { b.set_prob(s, 0.1 * (s + 1) as f64); }
        a.add_prob_dist(&b, 0, 1.0);
        for s in 0..5 {
            assert!((a.get_probability(s) - 0.1 * (s + 1) as f64).abs() < 1e-12);
        }
    }

    #[test]
    fn add_prob_dist_with_score_offset() {
        // self [0, 10), other [0, 5). add(other, +3, 1.0) shifts other's scores
        // by +3: self[3..8] += other[0..5].
        let mut a = ScoreDist::new(0, 10, false, true);
        let mut b = ScoreDist::new(0, 5, false, true);
        for s in 0..5 { b.set_prob(s, 0.2); }
        a.add_prob_dist(&b, 3, 1.0);
        for s in 0..3 { assert_eq!(a.get_probability(s), 0.0); }
        for s in 3..8 { assert!((a.get_probability(s) - 0.2).abs() < 1e-12); }
        for s in 8..10 { assert_eq!(a.get_probability(s), 0.0); }
    }

    #[test]
    fn add_prob_dist_with_negative_offset() {
        // self [-3, 5), other [0, 5). add(other, -2, 1.0) shifts down by 2.
        let mut a = ScoreDist::new(-3, 5, false, true);
        let mut b = ScoreDist::new(0, 5, false, true);
        for s in 0..5 { b.set_prob(s, 0.1); }
        a.add_prob_dist(&b, -2, 1.0);
        // other[0]→self[-2], other[4]→self[2]; self[-3] and self[3..5) untouched.
        assert_eq!(a.get_probability(-3), 0.0);
        for s in -2..3 { assert!((a.get_probability(s) - 0.1).abs() < 1e-12); }
        for s in 3..5 { assert_eq!(a.get_probability(s), 0.0); }
    }

    #[test]
    fn add_prob_dist_clips_to_self_range() {
        // self [0, 3), other [0, 5). add(other, 0, 1.0) only fills self[0..3].
        let mut a = ScoreDist::new(0, 3, false, true);
        let mut b = ScoreDist::new(0, 5, false, true);
        for s in 0..5 { b.set_prob(s, 0.2); }
        a.add_prob_dist(&b, 0, 1.0);
        for s in 0..3 { assert!((a.get_probability(s) - 0.2).abs() < 1e-12); }
    }

    #[test]
    fn add_prob_dist_scales_by_aa_prob() {
        let mut a = ScoreDist::new(0, 5, false, true);
        let mut b = ScoreDist::new(0, 5, false, true);
        for s in 0..5 { b.set_prob(s, 0.1); }
        a.add_prob_dist(&b, 0, 0.5);
        for s in 0..5 { assert!((a.get_probability(s) - 0.05).abs() < 1e-12); }
    }

    #[test]
    fn add_num_dist_with_coefficient() {
        let mut a = ScoreDist::new(0, 5, true, false);
        let mut b = ScoreDist::new(0, 5, true, false);
        for s in 0..5 { b.set_number(s, 2.0); }
        a.add_num_dist(&b, 0, 3.0);
        for s in 0..5 { assert!((a.get_number_recs(s) - 6.0).abs() < 1e-12); }
    }
```

- [ ] **Step 3.2: Verify failure**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust && \
  cargo test -p engine --lib gf::score_dist::tests 2>&1 | tail -30
```

Expected: compile errors (`add_prob_dist` / `add_num_dist` not defined).

- [ ] **Step 3.3: Implement the two methods**

In `engine/src/gf/score_dist.rs`, add to `impl ScoreDist`:

```rust
    /// Mirror Java's `ScoreDist.addProbDist(other, scoreDiff, aaProb)`:
    /// for each `t` in `other`'s score range, accumulate
    /// `other.prob[t] * aa_prob` into `self.prob[t + score_diff]`,
    /// clipping the destination to `self`'s range.
    pub fn add_prob_dist(&mut self, other: &ScoreDist, score_diff: i32, aa_prob: f64) {
        let other_p = match other.prob_distribution.as_ref() {
            Some(p) => p,
            None => return,
        };
        let self_p = match self.prob_distribution.as_mut() {
            Some(p) => p,
            None => return,
        };
        let other_min = other.bound.min_score;
        let other_max = other.bound.max_score;
        let self_min = self.bound.min_score;
        let self_max = self.bound.max_score;
        let t_start = other_min.max(self_min - score_diff);
        let t_end = other_max.min(self_max - score_diff);
        for t in t_start..t_end {
            let src_idx = (t - other_min) as usize;
            let dst_idx = (t + score_diff - self_min) as usize;
            self_p[dst_idx] += other_p[src_idx] * aa_prob;
        }
    }

    /// Mirror Java's `ScoreDist.addNumDist(other, scoreDiff, coeff)`.
    pub fn add_num_dist(&mut self, other: &ScoreDist, score_diff: i32, coeff: f64) {
        let other_n = match other.num_distribution.as_ref() {
            Some(n) => n,
            None => return,
        };
        let self_n = match self.num_distribution.as_mut() {
            Some(n) => n,
            None => return,
        };
        let other_min = other.bound.min_score;
        let other_max = other.bound.max_score;
        let self_min = self.bound.min_score;
        let self_max = self.bound.max_score;
        let t_start = other_min.max(self_min - score_diff);
        let t_end = other_max.min(self_max - score_diff);
        for t in t_start..t_end {
            let src_idx = (t - other_min) as usize;
            let dst_idx = (t + score_diff - self_min) as usize;
            self_n[dst_idx] += other_n[src_idx] * coeff;
        }
    }
```

- [ ] **Step 3.4: Run tests**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust && \
  cargo test -p engine --lib gf::score_dist::tests
```

Expected: all `score_dist::tests` pass (existing 8 + 6 new = 14).

- [ ] **Step 3.5: Commit**

```bash
git -C /Users/yperez/work/msgfplus-workspace/astral-speed add \
  rust/crates/engine/src/gf/score_dist.rs
git -C /Users/yperez/work/msgfplus-workspace/astral-speed commit -m "feat(engine): ScoreDist::add_prob_dist + add_num_dist (Phase 6/Task 3)"
```

---

## Task 4: `ScoredSpectrum` per-node and per-edge scoring

**Files:**
- Modify: `astral-speed/rust/crates/engine/src/scoring/scored_spectrum.rs`

The graph DP needs two new lookups on `ScoredSpectrum`:
- `node_score(prefix_nominal_mass, suffix_nominal_mass, scorer, charge) -> i32` — sum of prefix-ion scores at `prefix_nominal_mass` plus suffix-ion scores at `suffix_nominal_mass`. Iterates all ion types in each segment, computes `theo_mz`, looks up nearest peak rank, and accumulates `RankScorer::node_score(part, ion, rank)` (or `missing_ion_score` if no peak). Mirrors Java `NewScoredSpectrum.getNodeScore(prm, srm)`.
- `edge_score(cur_nominal_mass, prev_nominal_mass, theo_aa_mass, scorer, charge) -> i32` — combines `ion_existence_table` + `ion_err_dist_table` lookups, mirroring Java `NewScoredSpectrum.getEdgeScore`. For Phase 6 minimum, return `0` if `param.ion_existence_table.is_empty()` (matches Java `if (!scorer.supportEdgeScores()) return 0`). Otherwise apply the same formula as Java.

These methods take `&RankScorer` and `charge` so they can derive the partition (`charge`, `parent_mass`, `seg_num`) and reuse Phase 5's machinery.

- [ ] **Step 4.1: Add a small `ion_lookup` helper to `scoring/fragment_ions.rs`**

The existing `fragment_ions::predict_by_ions(peptide, charge_range)` predicts ions for an entire peptide. Phase 6 needs per-nominal-mass ion prediction: given a residue nominal mass and direction (`is_prefix`), enumerate ion types at the appropriate charges/segments and produce `(IonType, theo_mz)` pairs.

Add to `engine/src/scoring/fragment_ions.rs`:

```rust
use crate::param_model::{IonType, Param};

/// For a single prefix or suffix node at `nominal_mass`, enumerate the
/// `(ion_type, theo_mz)` pairs that contribute to its node score under
/// `param`. Java reference: NewScoredSpectrum.getNodeScore(nodeMass, isPrefix).
///
/// `is_prefix = true` → walk prefix ions (b-ions etc.); false → suffix (y-ions etc.).
/// `parent_mass`/`charge` select the segment+partition used downstream.
///
/// We return only the (IonType, theo_mz) pairs whose segment matches that
/// of `theo_mz` under `param.partition_for(charge, parent_mass, theo_mz)`.
pub fn ions_for_node(
    nominal_mass: f64,
    is_prefix: bool,
    param: &Param,
    parent_mass: f64,
    charge: u8,
) -> Vec<(IonType, f64)> {
    // For each segment in 0..param.num_segments, walk the ion types whose
    // direction matches `is_prefix`. Compute theo_mz = ion.mz(node_mass) and
    // include only if param.segment_num(theo_mz, parent_mass) == seg_index.
    //
    // The implementation below is a straight translation of
    // NewScoredSpectrum.getNodeScore(nodeMass, isPrefix). See that method
    // for the exact loop nesting; this function just produces the
    // (ion_type, theo_mz) tuples without scoring (callers do the scoring).
    let mut out = Vec::new();
    for seg in 0..param.num_segments {
        for &ion in param.ion_types_for_segment(seg).iter() {
            let theo_mz = match (is_prefix, ion) {
                (true, IonType::Prefix { .. }) => ion.mz(nominal_mass),
                (false, IonType::Suffix { .. }) => ion.mz(nominal_mass),
                _ => continue,
            };
            if param.segment_num(theo_mz, parent_mass) != seg {
                continue;
            }
            out.push((ion, theo_mz));
        }
    }
    out
}
```

> **Implementation note:** `param.ion_types_for_segment(seg)`, `param.segment_num(theo_mz, parent_mass)`, and `IonType::mz` may not exist yet on `Param` / `IonType`. Add them in this step:
>
> 1. **`IonType::mz(node_mass: f64) -> f64`** — for `Prefix { charge, offset_bits }`: `(node_mass + offset + charge * PROTON) / charge` where `offset = f32::from_bits(offset_bits) as f64`. For `Suffix`: same but `node_mass` is the suffix mass already in MS-GF+'s convention (verify against Java `IonType.PrefixIon.getMz` and `SuffixIon.getMz` — sources `edu.ucsd.msjava.msutil.IonType`).
> 2. **`Param::segment_num(theo_mz, parent_mass) -> usize`** — Java: `NewRankScorer.getSegmentNum(mass, parentMass)`. Returns `0` when `num_segments == 1`; otherwise `((theo_mz / parent_mass) * num_segments).floor().clamp(0, num_segments-1) as usize`. Verify against `NewRankScorer.java`.
> 3. **`Param::ion_types_for_segment(seg) -> &[IonType]`** — return the per-segment ion-type list. Currently the `rank_dist_table` is keyed by `Partition`; since `Partition` carries `seg_num`, derive the per-segment ion list by collecting `ion_type` keys from any partition with `seg_num == seg`. Cache eagerly inside `Param` if profiling shows hot.

Implement both `IonType::mz`, `Param::segment_num`, and `Param::ion_types_for_segment` in this step alongside `ions_for_node`. Tests for them go in their respective module test blocks (skeleton: `mz_of_b_ion_charge_1`, `segment_num_clamps`, `ion_types_for_segment_returns_unique`).

- [ ] **Step 4.2: Write failing tests for `node_score` / `edge_score`**

Append to `engine/src/scoring/scored_spectrum.rs`'s test module:

```rust
    use crate::param_model::{Partition, IonType};

    #[test]
    fn node_score_zero_when_no_peaks_present() {
        // Spectrum with no peaks; every ion is missing → all contributions
        // come from missing_ion_score. With an empty rank_dist_table the
        // missing_ion_score is 0 → node_score is 0.
        let s = spec(&[]);
        let param = tiny_param_with_ions();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let n = ss.node_score(100.0, 1000.0, &scorer, 2, 1100.0, 0.5);
        assert_eq!(n, 0);
    }

    #[test]
    fn node_score_positive_when_b1_peak_matches() {
        // Place a high-intensity peak at the predicted b1 m/z for a residue
        // of nominal mass = 100. node_score should be > 0.
        let nominal = 100.0;
        let proton = 1.007_276_49;
        let b1_mz = nominal + proton;  // charge=1, offset=0
        let s = spec(&[(50.0, 1.0), (b1_mz, 100.0), (200.0, 2.0)]);
        let param = tiny_param_with_ions();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let n = ss.node_score(nominal, 800.0, &scorer, 2, 900.0, 0.5);
        assert!(n > 0, "expected positive node_score, got {n}");
    }

    #[test]
    fn edge_score_returns_zero_when_table_empty() {
        // No ion_existence_table → Java path returns 0.
        let s = spec(&[(100.0, 1.0)]);
        let mut param = tiny_param_with_ions();
        param.ion_existence_table.clear();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let e = ss.edge_score(150, 100, 50.0, &scorer, 2, 1000.0, 0.5);
        assert_eq!(e, 0);
    }
```

> **Note:** `tiny_param_with_ions` is a helper that should be a richer version of `tiny_param` from `rank_scorer.rs::tests`: same skeleton but with a Prefix ion at charge=1 in `rank_dist_table` keyed by `Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 }`. Promote `tiny_param` to `pub(crate)` in `rank_scorer.rs` (gate behind `cfg(test)`) so the test fixture can be reused, OR copy/adapt as needed.

- [ ] **Step 4.3: Verify failure**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust && \
  cargo test -p engine --lib scoring::scored_spectrum 2>&1 | tail -20
```

Expected: compile errors (no `node_score` / `edge_score` methods).

- [ ] **Step 4.4: Implement `ScoredSpectrum::node_score` and `edge_score`**

Add to `impl<'a> ScoredSpectrum<'a>` in `scored_spectrum.rs`:

```rust
    /// Mirror Java NewScoredSpectrum.getNodeScore(prm, srm):
    /// `round(prefix_score(prm) + suffix_score(srm))`.
    ///
    /// `parent_mass` is the precursor neutral mass in Da (used to derive
    /// segment + partition). `fragment_tolerance_da` is the m/z window for
    /// `nearest_peak_rank`.
    pub fn node_score(
        &self,
        prefix_nominal: f64,
        suffix_nominal: f64,
        scorer: &RankScorer,
        charge: u8,
        parent_mass: f64,
        fragment_tolerance_da: f64,
    ) -> i32 {
        let pref = self.directional_node_score(prefix_nominal, /* is_prefix = */ true,
                                               scorer, charge, parent_mass, fragment_tolerance_da);
        let suff = self.directional_node_score(suffix_nominal, /* is_prefix = */ false,
                                               scorer, charge, parent_mass, fragment_tolerance_da);
        (pref + suff).round() as i32
    }

    fn directional_node_score(
        &self,
        nominal_mass: f64,
        is_prefix: bool,
        scorer: &RankScorer,
        charge: u8,
        parent_mass: f64,
        fragment_tolerance_da: f64,
    ) -> f32 {
        use crate::scoring::fragment_ions::ions_for_node;
        let mut total = 0.0_f32;
        for (ion, theo_mz) in ions_for_node(nominal_mass, is_prefix, scorer.param(), parent_mass, charge) {
            let seg = scorer.param().segment_num(theo_mz, parent_mass);
            let part = scorer.param().partition_for(charge, parent_mass, seg);
            match self.nearest_peak_rank(theo_mz, fragment_tolerance_da) {
                Some(rank) => total += scorer.node_score(part, ion, rank),
                None => total += scorer.missing_ion_score(part, ion),
            }
        }
        total
    }

    /// Mirror Java NewScoredSpectrum.getEdgeScore(curNode, prevNode, theoMass).
    /// If `param.ion_existence_table` is empty (Java's `!scorer.supportEdgeScores()`),
    /// return 0. Otherwise:
    ///   - `ion_existence_index` = (cur_node_mass >= 0) + 2*(prev_node_mass >= 0)
    ///     where node_mass uses observed peak m/z if present else -1.
    ///   - score = ion_existence_score + (if both observed) error_score
    ///
    /// `theo_aa_mass` is the AA's accurate mass in Da (used for the error term).
    pub fn edge_score(
        &self,
        cur_nominal: i32,
        prev_nominal: i32,
        theo_aa_mass: f64,
        scorer: &RankScorer,
        charge: u8,
        parent_mass: f64,
        fragment_tolerance_da: f64,
    ) -> i32 {
        if scorer.param().ion_existence_table.is_empty() {
            return 0;
        }
        // 1. Observed peak m/z for cur and prev nodes (using main_ion direction).
        // 2. ion_existence_index ∈ {0..3}.
        // 3. ion_existence_score = scorer.ion_existence_score(part, idx, prob_peak).
        // 4. If idx == 3, add error_score(part, observed_delta - theo_aa_mass).
        // 5. Round to i32.
        //
        // Implementer: see Java NewScoredSpectrum.getEdgeScore for the exact
        // formula. `prob_peak` is computed once per spectrum at construction
        // time in Java (`probPeak = spec.size() / approxNumBins`); store it
        // on `ScoredSpectrum` in this task and pass it through.
        //
        // Pseudocode (translate exactly):
        //   let cur_mass = self.observed_node_mass(cur_nominal, scorer, charge, parent_mass, fragment_tolerance_da);
        //   let prev_mass = self.observed_node_mass(prev_nominal, scorer, charge, parent_mass, fragment_tolerance_da);
        //   let mut idx = 0;
        //   if cur_mass.is_some() { idx += 1; }
        //   if prev_mass.is_some() { idx += 2; }
        //   let part = scorer.param().partition_for(charge, parent_mass, /* seg = */ 0);
        //   let mut s = scorer.ion_existence_score(part, idx, self.prob_peak);
        //   if idx == 3 {
        //       let delta = cur_mass.unwrap() - prev_mass.unwrap() - theo_aa_mass;
        //       s += scorer.error_score(part, delta);
        //   }
        //   s.round() as i32
        unimplemented!("translate per Java NewScoredSpectrum.getEdgeScore + see commentary above")
    }
```

> **Implementation guidance for the unimplemented body:**
>
> 1. **`prob_peak`**: Java computes `spec.size() / approxNumBins` (where `approxNumBins ≈ scoredSpec.parentMass / mme.toleranceAsDa()`). Compute and cache this on `ScoredSpectrum` (add a `prob_peak: f32` field in the constructor; for `new_without_filtering`, default to a sentinel like `1.0` and document that callers must use `new` for production).
> 2. **`observed_node_mass(node_nominal, ...)`** is `NewScoredSpectrum.getNodeMass(node)`: compute `theo_mz = main_ion.mz(node_nominal)`, find the nearest peak; if found, return `main_ion.mass_from_mz(peak_mz)`; else return `None`.
> 3. **`scorer.ion_existence_score(part, idx, prob_peak)`** and **`scorer.error_score(part, delta)`** are new methods on `RankScorer` driven by `Param.ion_existence_table` and `Param.ion_err_dist_table` / `Param.noise_err_dist_table`. Add them in this task, mirroring `NewRankScorer.getIonExistenceScore` / `getErrorScore`.
> 4. **`main_ion`**: Java's `NewScoredSpectrum.mainIon` is selected per-segment in the constructor as the highest-frequency prefix ion. For Rust, store the main ion type on `ScoredSpectrum` (compute it from `param.rank_dist_table` for the partition matching `(charge, parent_mass, seg=0)` and pick the ion with the largest `node_score(rank=1)`).

- [ ] **Step 4.5: Run tests + commit**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust && \
  cargo test -p engine --lib scoring 2>&1 | tail -30
```

Expected: all `scoring::*::tests` pass (Phase 5's existing + 8 new).

```bash
git -C /Users/yperez/work/msgfplus-workspace/astral-speed add \
  rust/crates/engine/src/scoring/ \
  rust/crates/engine/src/param_model.rs
git -C /Users/yperez/work/msgfplus-workspace/astral-speed commit -m "feat(engine): ScoredSpectrum::{node_score,edge_score} + RankScorer extensions (Phase 6/Task 4)"
```

---

## Task 5: `PrimitiveAaGraph` (CSR amino-acid graph)

**Files:**
- Create: `astral-speed/rust/crates/engine/src/gf/primitive_graph.rs`
- Modify: `astral-speed/rust/crates/engine/src/gf/mod.rs`

Direct port of Java `edu.ucsd.msjava.msgf.PrimitiveAminoAcidGraph` (this branch's CSR-format graph used by `DBScanner`). Storage:
- `active_nodes: Vec<i32>` — list of reachable nominal masses, sorted ascending; `active_nodes[0] = 0` (source), `active_nodes[node_count - 1] = peptide_mass` (sink).
- `mass_to_node_idx: Vec<i32>` — dense `mass + mass_offset → node_idx` (or `-1`).
- `edge_offset: Vec<usize>` (length `node_count + 1`) — CSR offsets so that incoming edges of `node_idx` live in `edge_offset[ni]..edge_offset[ni + 1]`.
- `edge_prev_node: Vec<i32>` — predecessor nominal mass per edge.
- `edge_prob: Vec<f32>` — AA probability per edge.
- `edge_score: Vec<i32>` — `cleavage_score + error_score` (precomputed).
- `node_scores: Vec<i32>` — per-node summed prefix+suffix score from `ScoredSpectrum::node_score`.

Construction order (mirror Java `PrimitiveAminoAcidGraph` constructor):
1. Resolve source/sink AA lists from `aa_set` honoring `Location::{N_Term, C_Term, Protein_N_Term, Protein_C_Term}` and the spectrum's main-ion direction.
2. Discover reachable masses via three forward sweeps (source→AAs, intermediates→AAs, sink-backward); accumulate `in_edge_count_by_mass`.
3. Allocate `active_nodes`, `mass_to_node_idx`, `edge_offset` from cumulative counts.
4. Fill CSR edges (source forward, intermediates forward, sink backward), recording `prev_mass`, `prob`, `mass`, `cleavage_score`.
5. Compute per-edge error score via `ScoredSpectrum::edge_score` (Task 4) and add to `cleavage_score`.
6. Compute per-node score via `ScoredSpectrum::node_score` (Task 4), respecting direction (Java's `if (!direction)` branch swaps prefix/suffix).

- [ ] **Step 5.1: Stub the file + module wiring**

Create `engine/src/gf/primitive_graph.rs`:

```rust
//! CSR-format amino-acid graph used by the GF DP. Mirrors Java
//! `edu.ucsd.msjava.msgf.PrimitiveAminoAcidGraph`.

use crate::aa_set::AminoAcidSet;
use crate::enzyme::Enzyme;
use crate::scoring::{RankScorer, ScoredSpectrum};

#[derive(Debug, Clone)]
pub struct PrimitiveAaGraph {
    pub peptide_mass: i32,
    pub direction: bool,                // true = prefix-direction main ion
    pub min_node_mass: i32,
    pub mass_offset: i32,
    pub node_count: usize,
    pub source_node_idx: usize,
    pub sink_node_idx: usize,

    pub active_nodes: Vec<i32>,
    pub mass_to_node_idx: Vec<i32>,
    pub edge_offset: Vec<usize>,
    pub edge_prev_node: Vec<i32>,
    pub edge_prob: Vec<f32>,
    pub edge_score: Vec<i32>,
    pub node_scores: Vec<i32>,

    pub aa_set_handle: AminoAcidSet,
    pub enzyme: Option<Enzyme>,
}

impl PrimitiveAaGraph {
    /// Construct the graph for `peptide_mass` (nominal) using `aa_set`,
    /// `enzyme`, and the per-spectrum scoring view. `use_protein_n_term` /
    /// `use_protein_c_term` select Location::Protein_*_Term variants.
    pub fn new(
        aa_set: &AminoAcidSet,
        peptide_mass: i32,
        enzyme: Option<&Enzyme>,
        scored_spec: &ScoredSpectrum<'_>,
        scorer: &RankScorer,
        charge: u8,
        parent_mass: f64,
        fragment_tolerance_da: f64,
        use_protein_n_term: bool,
        use_protein_c_term: bool,
    ) -> Self {
        let _ = (aa_set, peptide_mass, enzyme, scored_spec, scorer, charge,
                 parent_mass, fragment_tolerance_da, use_protein_n_term, use_protein_c_term);
        unimplemented!("Phase 6 Task 5: port PrimitiveAminoAcidGraph constructor")
    }

    pub fn node_index_for_mass(&self, mass: i32) -> Option<usize> {
        if mass < self.min_node_mass || mass > self.peptide_mass {
            return None;
        }
        let dense = (mass + self.mass_offset) as usize;
        let v = self.mass_to_node_idx[dense];
        if v < 0 { None } else { Some(v as usize) }
    }
}

#[cfg(test)]
mod tests {
    // tests added in Step 5.3
}
```

In `engine/src/gf/mod.rs` add:

```rust
pub mod primitive_graph;
pub use primitive_graph::PrimitiveAaGraph;
```

- [ ] **Step 5.2: Write failing tests**

Append to `primitive_graph.rs::mod tests`:

```rust
    use super::*;
    use crate::aa_set::AminoAcidSetBuilder;
    use crate::scoring::{RankScorer, ScoredSpectrum};
    use crate::spectrum::Spectrum;
    use crate::scoring::rank_scorer::tests as rs_tests;  // promote tiny_param to pub(crate)

    fn empty_spec() -> Spectrum {
        Spectrum {
            title: "t".into(), precursor_mz: 500.0, precursor_intensity: None,
            precursor_charge: Some(2), rt_seconds: None, scan: None, peaks: vec![],
        }
    }

    #[test]
    fn graph_for_peptide_mass_zero_has_only_source_and_sink() {
        // peptide_mass = 0 is degenerate — sink == source.
        let aa = AminoAcidSetBuilder::new_standard().build().unwrap();
        let s = empty_spec();
        let param = rs_tests::tiny_param();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let g = PrimitiveAaGraph::new(&aa, 0, None, &ss, &scorer, 2, 1000.0, 0.5, false, false);
        assert_eq!(g.node_count, 1);
        assert_eq!(g.source_node_idx, g.sink_node_idx);
    }

    #[test]
    fn graph_active_nodes_contain_source_and_sink() {
        let aa = AminoAcidSetBuilder::new_standard().build().unwrap();
        let s = empty_spec();
        let param = rs_tests::tiny_param();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        // Pick a peptide_mass reachable by two alanines (each ≈ 71): 142.
        let g = PrimitiveAaGraph::new(&aa, 142, None, &ss, &scorer, 2, 1000.0, 0.5, false, false);
        assert_eq!(g.active_nodes.first(), Some(&0));
        assert_eq!(g.active_nodes.last(), Some(&142));
        assert!(g.node_count >= 2);
    }

    #[test]
    fn csr_edge_offsets_are_monotonic() {
        let aa = AminoAcidSetBuilder::new_standard().build().unwrap();
        let s = empty_spec();
        let param = rs_tests::tiny_param();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let g = PrimitiveAaGraph::new(&aa, 200, None, &ss, &scorer, 2, 1000.0, 0.5, false, false);
        for w in g.edge_offset.windows(2) {
            assert!(w[0] <= w[1], "edge_offset must be non-decreasing");
        }
        assert_eq!(g.edge_offset[g.node_count], g.edge_prev_node.len());
        assert_eq!(g.edge_offset[g.node_count], g.edge_prob.len());
        assert_eq!(g.edge_offset[g.node_count], g.edge_score.len());
    }

    #[test]
    fn enzyme_credit_added_to_sink_edges_when_c_term_enzyme() {
        // With Trypsin and main_ion direction = prefix (forward direction true),
        // sink edges (cleavage on C-term residue) should get cleavage credit on
        // K/R prev-residues and penalty on others.
        let aa = AminoAcidSetBuilder::new_standard().build().unwrap();
        let s = empty_spec();
        let param = rs_tests::tiny_param();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let g = PrimitiveAaGraph::new(&aa, 200, Some(&Enzyme::Trypsin), &ss, &scorer, 2, 1000.0, 0.5, false, false);
        // At least one sink edge must have a non-zero score component contributed
        // by the enzyme. The exact value depends on aa_set's cleavage credit.
        let sink_edges = g.edge_offset[g.sink_node_idx + 1] - g.edge_offset[g.sink_node_idx];
        assert!(sink_edges > 0);
    }
```

> **Helper trick:** to share the `tiny_param` fixture from `rank_scorer::tests`, add `#[cfg(test)] pub(crate) mod tests` and re-export `tiny_param` via `pub(crate) fn tiny_param()` in that module.

- [ ] **Step 5.3: Verify failure** (`cargo test -p engine --lib gf::primitive_graph`).

- [ ] **Step 5.4: Implement the constructor**

Translate `PrimitiveAminoAcidGraph` constructor + `computeEdgeErrorScores` + `computeNodeScores` from Java to Rust. The structural code (sweeps, CSR fill) is mechanical. The two non-trivial bits:

1. **`Enzyme::is_n_term()` / `Enzyme::is_c_term()` / `Enzyme::is_cleavable(aa)`** — these may need to be added to `engine::enzyme` if not present (check current `Enzyme` enum). Mirror `edu.ucsd.msjava.msutil.Enzyme`.
2. **`AminoAcidSet::aa_list_for(location: Location) -> &[AminoAcid]`** — Java has `getAAList(Location)`. Add if missing; for `Location::Anywhere` return the standard list, for terminal locations apply terminal-specific mods. Mirror `edu.ucsd.msjava.msutil.AminoAcidSet`.

Use the Java reference at the top of this task as the source of truth for the loop structure.

- [ ] **Step 5.5: Run tests + commit**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust && \
  cargo test -p engine --lib gf::primitive_graph
git -C /Users/yperez/work/msgfplus-workspace/astral-speed add \
  rust/crates/engine/src/gf/ \
  rust/crates/engine/src/aa_set.rs \
  rust/crates/engine/src/enzyme.rs
git -C /Users/yperez/work/msgfplus-workspace/astral-speed commit -m "feat(engine): PrimitiveAaGraph CSR amino-acid graph (Phase 6/Task 5)"
```

---

## Task 6: Graph-based `GeneratingFunction::compute`

**Files:**
- Modify: `astral-speed/rust/crates/engine/src/gf/generating_function.rs`

Replace the placeholder DP from Phase 6/Task 2 with the real graph DP. Keep the existing API as `compute_uniform` (rename) for tests; introduce the new public API:

```rust
impl GeneratingFunction {
    /// Compute the GF over a precomputed primitive graph. Mirrors Java
    /// PrimitiveGeneratingFunction.computeGeneratingFunction().
    ///
    /// Returns Err if the graph yields an empty score range (degenerate).
    pub fn compute(graph: &PrimitiveAaGraph, aa_set: &AminoAcidSet) -> Result<Self, GfError>;

    /// Optional pre-pass: prune nodes whose maximum possible final score is
    /// below `score_threshold`. Mirrors Java setUpScoreThreshold; call before
    /// `compute` to skip irrelevant DP work.
    pub fn with_score_threshold(graph: &PrimitiveAaGraph, score_threshold: i32, aa_set: &AminoAcidSet) -> Result<Self, GfError>;

    pub fn score_dist(&self) -> &ScoreDist;
    pub fn min_score(&self) -> i32;
    pub fn max_score(&self) -> i32;
    pub fn spectral_probability(&self, score: i32) -> f64;
}
```

DP loop (mirror Java `PrimitiveGeneratingFunction.computeGeneratingFunction:89-205`):
- Initialize `dist_by_node[source] = ScoreDist::new(0, 1, false, true)` with `prob[0] = 1.0`.
- For each `node_idx` in `0..node_count` skipping the source:
  - Let `cur_node_score = graph.node_scores[ni]`.
  - For each edge `e` in `[edge_offset[ni], edge_offset[ni + 1])` whose predecessor has a `dist` set:
    - `combined = cur_node_score + edge_score[e]`
    - Track `cur_min/max_score` from predecessor min/max + combined.
  - If `cur_min >= cur_max` or out of `[-10000, 10000]`, skip.
  - Allocate `cur_dist = ScoreDist::new(cur_min, cur_max, false, true)`.
  - For each valid edge: `cur_dist.add_prob_dist(prev_dist, combined, edge_prob[e] as f64)`.
  - Store underflow guard at `cur_dist.max_score - 1` (clamp to f32::MIN_POSITIVE).
  - `dist_by_node[ni] = Some(cur_dist)`.
- Pull sink dist; if missing or empty range → return `Err(GfError::EmptyScoreRange)`.
- Apply enzyme adjustment exactly as Java does:
  - If `enzyme.is_some()` and `enzyme.residues().is_some()`:
    - `final_dist = ScoreDist::new(min + penalty, max + credit, false, true)`
    - `final_dist.add_prob_dist(sink_dist, credit, prob_cleavage_sites)`
    - `final_dist.add_prob_dist(sink_dist, penalty, 1 - prob_cleavage_sites)`
  - Else: `final_dist = sink_dist.clone()`.
- Return `Self { score_dists: vec![final_dist], score_bound: ScoreBound::new(min, max) }`.

> **Where does `prob_cleavage_sites` come from?** It's `AminoAcidSet::prob_cleavage_sites()` in Java — sum of probabilities of cleavable AAs divided by 1.0. Mirror in `aa_set.rs`. Cache eagerly per `AminoAcidSet`.

- [ ] **Step 6.1: Rename existing `compute` to `compute_uniform`** (keep tests).

- [ ] **Step 6.2: Write failing tests for new `compute(graph)`**

Create `astral-speed/rust/crates/engine/tests/gf_graph_dp.rs`:

```rust
//! GF DP smoke tests on hand-built graphs.

use engine::aa_set::AminoAcidSetBuilder;
use engine::enzyme::Enzyme;
use engine::gf::{GeneratingFunction, PrimitiveAaGraph};
use engine::scoring::{RankScorer, ScoredSpectrum};
use engine::spectrum::Spectrum;

fn aa() -> engine::aa_set::AminoAcidSet {
    AminoAcidSetBuilder::new_standard().build().unwrap()
}

fn empty_spec() -> Spectrum {
    Spectrum {
        title: "t".into(), precursor_mz: 500.0, precursor_intensity: None,
        precursor_charge: Some(2), rt_seconds: None, scan: None, peaks: vec![],
    }
}

#[test]
fn gf_on_trivial_graph_has_max_score_one() {
    // peptide_mass = 0 → sink == source; the GF dist should have a single
    // point at score 0 with probability 1.
    let aa = aa();
    let s = empty_spec();
    let param = engine::scoring::rank_scorer::tests::tiny_param();
    let scorer = RankScorer::new(&param);
    let ss = ScoredSpectrum::new_without_filtering(&s);
    let g = PrimitiveAaGraph::new(&aa, 0, None, &ss, &scorer, 2, 0.0, 0.5, false, false);
    let gf = GeneratingFunction::compute(&g, &aa).expect("trivial GF");
    assert!(gf.spectral_probability(0) >= 0.999);
}

#[test]
fn gf_score_dist_sums_to_one_no_enzyme_no_score() {
    // peptide_mass that is reachable by AAs; node/edge scores are all 0
    // (empty spectrum, no ion existence table). Total prob across the
    // sink dist should equal 1.0.
    let aa = aa();
    let s = empty_spec();
    let param = engine::scoring::rank_scorer::tests::tiny_param();
    let scorer = RankScorer::new(&param);
    let ss = ScoredSpectrum::new_without_filtering(&s);
    let g = PrimitiveAaGraph::new(&aa, 200, None, &ss, &scorer, 2, 0.0, 0.5, false, false);
    let gf = GeneratingFunction::compute(&g, &aa).expect("non-empty GF");
    let dist = gf.score_dist();
    let total: f64 = (dist.min_score()..dist.max_score())
        .map(|s| dist.get_probability(s)).sum();
    assert!((total - 1.0).abs() < 1e-6, "total prob = {total}");
}

#[test]
fn gf_spectral_probability_monotonic_decreasing() {
    let aa = aa();
    let s = empty_spec();
    let param = engine::scoring::rank_scorer::tests::tiny_param();
    let scorer = RankScorer::new(&param);
    let ss = ScoredSpectrum::new_without_filtering(&s);
    let g = PrimitiveAaGraph::new(&aa, 250, None, &ss, &scorer, 2, 0.0, 0.5, false, false);
    let gf = GeneratingFunction::compute(&g, &aa).expect("GF");
    let dist = gf.score_dist();
    let mut prev = f64::INFINITY;
    for s in dist.min_score()..dist.max_score() {
        let p = gf.spectral_probability(s);
        assert!(p <= prev + 1e-12, "P should be non-increasing in score");
        prev = p;
    }
}
```

- [ ] **Step 6.3: Verify failure** (`cargo test -p engine --test gf_graph_dp`).

- [ ] **Step 6.4: Implement `compute` and `with_score_threshold`**

Translate the Java DP. Maintain a `Vec<Option<ScoreDist>>` indexed by node. Keep a scratch `Vec<usize>` for valid-edge indices to avoid reallocating per iteration.

- [ ] **Step 6.5: Run tests + commit**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust && \
  cargo test -p engine --test gf_graph_dp && \
  cargo test -p engine --lib gf
git -C /Users/yperez/work/msgfplus-workspace/astral-speed add \
  rust/crates/engine/src/gf/ \
  rust/crates/engine/tests/gf_graph_dp.rs
git -C /Users/yperez/work/msgfplus-workspace/astral-speed commit -m "feat(engine): graph-based GeneratingFunction::compute (Phase 6/Task 6)"
```

---

## Task 7: `GeneratingFunctionGroup` (multi-bin merger)

**Files:**
- Create: `astral-speed/rust/crates/engine/src/gf/group.rs`
- Modify: `astral-speed/rust/crates/engine/src/gf/mod.rs`

Direct port of Java `PrimitiveGeneratingFunctionGroup`. Streaming merger that accepts each per-mass-bin `GeneratingFunction`, merges its `ScoreDist` into a running aggregate, and lets the input `gf` go out of scope (Rust ownership makes this natural via `accept(gf: GeneratingFunction)` taking by value).

```rust
pub struct GeneratingFunctionGroup {
    min_score: i32,
    max_score: i32,
    merged: Option<ScoreDist>,
}

impl GeneratingFunctionGroup {
    pub fn new() -> Self { Self { min_score: i32::MAX, max_score: i32::MIN, merged: None } }

    pub fn accept(&mut self, gf: GeneratingFunction) {
        let dist = gf.score_dist();
        let gf_min = dist.min_score();
        let gf_max = dist.max_score();
        if self.merged.is_none() {
            self.min_score = gf_min;
            self.max_score = gf_max;
            let mut m = ScoreDist::new(gf_min, gf_max, false, true);
            m.add_prob_dist(dist, 0, 1.0);
            self.merged = Some(m);
            return;
        }
        let new_min = self.min_score.min(gf_min);
        let new_max = self.max_score.max(gf_max);
        if new_min != self.min_score || new_max != self.max_score {
            let mut expanded = ScoreDist::new(new_min, new_max, false, true);
            expanded.add_prob_dist(self.merged.as_ref().unwrap(), 0, 1.0);
            self.merged = Some(expanded);
            self.min_score = new_min;
            self.max_score = new_max;
        }
        self.merged.as_mut().unwrap().add_prob_dist(dist, 0, 1.0);
    }

    pub fn is_computed(&self) -> bool { self.merged.is_some() }
    pub fn max_score(&self) -> i32 { self.max_score }
    pub fn score_dist(&self) -> Option<&ScoreDist> { self.merged.as_ref() }
    pub fn spectral_probability(&self, score: i32) -> Option<f64> {
        self.merged.as_ref().map(|d| d.get_spectral_probability(score))
    }
}
```

- [ ] **Step 7.1: Write failing tests** for: `empty_group_returns_none`, `single_gf_merge_preserves_dist`, `two_gfs_merge_sum_of_probs`, `expanding_range_keeps_existing_mass`, `spectral_probability_after_merge_clamped_to_one`.

- [ ] **Step 7.2: Implement** as above.

- [ ] **Step 7.3: Run + commit**

```bash
git -C /Users/yperez/work/msgfplus-workspace/astral-speed commit -m "feat(engine): GeneratingFunctionGroup multi-bin merger (Phase 6/Task 7)"
```

---

## Task 8: Wire SpecEValue into `PsmMatch` + `match_engine`

**Files:**
- Modify: `astral-speed/rust/crates/engine/src/psm.rs`
- Modify: `astral-speed/rust/crates/engine/src/match_engine.rs`

Add `spec_e_value: f64` to `PsmMatch` (default `1.0` until computed). Update `Ord` to compare on `(spec_e_value ascending, then score descending)` so the "best" PSM has the *smallest* SpecEValue (matches Java).

In `match_engine.rs`, after the per-spectrum scoring loop, build the GF group:

```rust
// For each spectrum, after scoring all candidates and pushing PSMs to queues[spec_idx]:
//   1. Compute the precursor neutral mass + tolerance window in nominal mass space.
//   2. For each nominal mass M in [min_mass_idx, max_mass_idx]:
//      - Build PrimitiveAaGraph for M.
//      - Compute GeneratingFunction on it (optionally with score threshold = current best score).
//      - group.accept(gf).
//   3. For each PSM in queues[spec_idx], assign:
//        psm.spec_e_value = group.spectral_probability(psm.score.round() as i32).unwrap_or(1.0);
```

The current `score: f32` is the rank score from Phase 5. SpecEValue is computed *from* that score, so this step does not change the per-candidate score function — it only computes the e-value AFTER scoring is done.

Mirror Java DBScanner.computeSpecEValues path:
- `min_score = min(psm.score for psm in queues[spec_idx])` (used for `with_score_threshold`).
- `peptide_mass_da = (precursor_mz - PROTON) * charge - H2O`
- `nominal_peptide_mass = NominalMass::from_mass_da(peptide_mass_da)`
- `[min_idx, max_idx] = nominal_peptide_mass ± isotope_error ± tolerance_da_round`

- [ ] **Step 8.1: Add `spec_e_value` field + update `Ord`/`PartialOrd`/tests in `psm.rs`**

- [ ] **Step 8.2: Add SpecEValue computation pass to `match_spectra`** (separate function `compute_spec_e_values_for_spectrum` so Task 9 can unit-test it).

- [ ] **Step 8.3: Smoke test on BSA + test.mgf** — assert that for each non-empty queue, the top PSM has `spec_e_value < 1.0` (a non-trivial SpecEValue).

- [ ] **Step 8.4: Commit**

```bash
git -C /Users/yperez/work/msgfplus-workspace/astral-speed commit -m "feat(engine): wire SpecEValue into PsmMatch + match_engine (Phase 6/Task 8)"
```

---

## Task 9: Java parity on traced spectra

**Files:**
- Create: `astral-speed/rust/crates/engine/tests/gf_java_parity.rs`
- Create: `astral-speed/src/test/resources/gf-traces/expected_spec_evalues.tsv` (small, hand-curated)

Capture Java's SpecEValue for ~5 traced (spectrum, peptide, charge) triples on `BSA + test.mgf`. The Rust GF for the same (spectrum, peptide, charge) must match within tolerance.

- [ ] **Step 9.1: Capture Java reference values**

Run MS-GF+ via the existing `benchmark/` harness on BSA + test.mgf. From the resulting `.tsv`, pick 5 high-confidence PSMs and record:
- spectrum title (or scan #)
- peptide sequence
- charge
- SpecEValue

Save as `src/test/resources/gf-traces/expected_spec_evalues.tsv`:

```
title<TAB>peptide<TAB>charge<TAB>spec_evalue
KQTALVELLK_4567<TAB>KQTALVELLK<TAB>2<TAB>1.23e-12
…
```

If running Java is impractical in-loop, hard-code 5 values with provenance comments referencing the harness output filename + line.

- [ ] **Step 9.2: Write the parity test**

```rust
//! Java SpecEValue parity for traced (spectrum, peptide, charge) triples.

use std::fs;
use engine::{...};

const TOLERANCE_LOG10: f64 = 1.0;  // within 1 order of magnitude

#[test]
fn rust_spec_evalue_within_one_oom_of_java() {
    let cases = parse_traces("src/test/resources/gf-traces/expected_spec_evalues.tsv");
    let (spectra, idx, params, scorer) = setup_bsa_test_mgf();
    let queues = match_spectra(...);

    for case in cases {
        let spec_idx = locate_spectrum(&spectra, &case.title);
        let pep_idx = locate_psm(&queues[spec_idx], &case.peptide, case.charge);
        let rust_e = queues[spec_idx].into_sorted_vec()[pep_idx].spec_e_value;
        let log_diff = (rust_e.log10() - case.spec_evalue.log10()).abs();
        assert!(log_diff < TOLERANCE_LOG10,
            "{}: Java {:.2e} vs Rust {:.2e} ({}x)",
            case.title, case.spec_evalue, rust_e, log_diff);
    }
}
```

- [ ] **Step 9.3: Run + iterate**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust && \
  cargo test -p engine --test gf_java_parity 2>&1 | tail -30
```

If the test fails, the typical root causes are (in order of frequency):
1. **Edge score formula** — `ion_existence_score` / `error_score` mismatch with Java. Compare per-edge values on a single PSM.
2. **Main ion direction** — Java's `getMainIonDirection` may pick a different main ion than Rust. Print `direction` for the failing spectrum.
3. **Enzyme cleavage adjustment** — `prob_cleavage_sites` or credit/penalty mismatch. Compare `aa_set.prob_cleavage_sites()`.
4. **Mass-bin window** — `min_idx / max_idx` computation. Java rounds with `Math.round(tolDa - 0.4999f)`; Rust must match exactly.
5. **Underflow guard** — Java sets `min(prob)` at `Float.MIN_VALUE` if zero; Rust must do the same.

Document each diagnosis in the failing test's commit message. Tighten `TOLERANCE_LOG10` to `0.5` once 5/5 pass at `1.0`; revise the gate downward as confidence grows.

- [ ] **Step 9.4: Commit**

```bash
git -C /Users/yperez/work/msgfplus-workspace/astral-speed commit -m "test(engine): SpecEValue Java parity on 5 traced PSMs (Phase 6/Task 9)"
```

---

## Task 10: BSA + test.mgf SpecEValue parity gate

**Files:**
- Create: `astral-speed/rust/crates/engine/tests/gf_bsa_parity.rs`

For *every* PSM Java identifies on BSA + test.mgf, Rust's SpecEValue should agree within an OOM (Phase 6 acceptance criterion). This is the harder gate — it surfaces tail behavior the hand-picked traces in Task 9 might miss.

- [ ] **Step 10.1: Generate Java reference for all 217 PSMs**

Reuse the existing per-spectrum top-1 parity test from Phase 5/Task 6 (`match_engine_top1_parity`). Java's reference output already lists peptide and SpecEValue. Extend the fixture to include `spec_evalue` per (spectrum, peptide).

- [ ] **Step 10.2: Write the bulk parity test**

```rust
#[test]
fn bsa_test_mgf_spec_evalue_parity_gate_95pct_within_1_oom() {
    let cases = load_java_reference("src/test/resources/match_engine_top1_reference.tsv");
    let (spectra, idx, params, scorer) = setup_bsa_test_mgf();
    let queues = match_spectra(...);

    let mut within_oom = 0;
    let mut total = 0;
    for case in &cases {
        let rust_top = queues[case.spec_idx].clone().into_sorted_vec()[0];
        if rust_top.candidate.peptide.to_string() != case.peptide { continue; }
        total += 1;
        let log_diff = (rust_top.spec_e_value.log10() - case.spec_evalue.log10()).abs();
        if log_diff < 1.0 { within_oom += 1; }
    }
    assert!(total >= 200, "expected >= 200 matched PSMs, got {total}");
    let pct = within_oom as f32 / total as f32;
    assert!(pct >= 0.95,
        "{within_oom}/{total} ({:.1}%) within 1 OOM — gate is 95%", pct * 100.0);
}
```

- [ ] **Step 10.3: Run + iterate** (same diagnosis flow as Task 9; expect more cases to fix here).

- [ ] **Step 10.4: Commit**

```bash
git -C /Users/yperez/work/msgfplus-workspace/astral-speed commit -m "test(engine): BSA SpecEValue parity gate ≥95% within 1 OOM (Phase 6/Task 10)"
```

---

## Phase 6 done — exit gate

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust

cargo test --workspace --no-fail-fast 2>&1 | grep "test result:" | tail -25
# Expected: ~330 tests pass (Phase 5 ~290 + Phase 6 ~43).

cargo clippy --workspace --all-targets -- -D warnings
# Expected: clean.

find crates -name '*.rs' -path '*/src/*' | xargs wc -l | tail -1
# Expected: ~7,000 LOC total (Phase 5 ~5,800 + Phase 6 ~1,200).
```

**Milestone commit summary** (use as the body of the closing PR's description):

```
Phase 6 ships the graph-based generating-function DP. Per-spectrum SpecEValue
is now computed for every PSM, replacing the rank-score-only ranking from Phase 5.
- Tasks 1-2: ScoreBound/ScoreDist data + placeholder DP (already shipped fcc0f05/779f617).
- Task 3: ScoreDist::add_prob_dist + add_num_dist.
- Task 4: ScoredSpectrum::node_score + edge_score (per-nominal-mass).
- Task 5: PrimitiveAaGraph CSR amino-acid graph.
- Task 6: graph-based GeneratingFunction::compute + with_score_threshold.
- Task 7: GeneratingFunctionGroup multi-bin merger.
- Task 8: SpecEValue wired into PsmMatch + match_engine.
- Task 9: Java parity on 5 traced PSMs.
- Task 10: BSA + test.mgf parity gate ≥ 95% within 1 OOM.
Next: Phase 7 — output writers (.tsv, .pin) + Percolator integration.
```

After Phase 6's exit gate passes, **immediately proceed to writing the Phase 7 plan** (`docs/superpowers/2026-05-05-phase7-output-writers-plan.md`) using the same TDD-task template as this plan. Phase 7 covers:
- `output::tsv` — Java `DirectTSVWriter` parity.
- `output::pin` — Java `DirectPinWriter` parity (Percolator-consumable).
- Wire `cli` to actually write outputs end-to-end (currently the CLI exists but no write path).
- Round-trip integration test: read MGF → search → write `.pin` → run Percolator (or skip Percolator and assert the `.pin` schema matches Java's).

## Self-review checklist (run before handing off to executor)

- [x] Tasks 1-2 status declared up front (already shipped).
- [x] Each task has files / steps / commands / commit message.
- [x] No `TODO`/`fill-in` placeholders in steps that produce code (Task 4's edge_score body is explicitly an `unimplemented!()` stub with a translation guide because the Java reference is short and the implementer should read it directly — this is a deliberate handoff point, not an omission).
- [x] Type names consistent across tasks (`PrimitiveAaGraph`, `GeneratingFunction`, `GeneratingFunctionGroup`, `PsmMatch`).
- [x] Iteration shipping model honored: milestone commits on `rust-implement`, single closing PR at end (no per-task PR).
- [x] Exit gate has measurable success criteria (test count, clippy, LOC).
- [x] Phase 7 handoff explicitly noted at the bottom.
