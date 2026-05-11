# Known Java↔Rust Divergences (Parity Follow-Ups)

_Last updated: 2026-05-10_

This is the canonical register of **explicit follow-ups, not regressions**, that
distinguish msgf-rust from Java MS-GF+ on the discovery-style pinned BSA +
test.mgf and PXD001819 fixtures. They are the gaps to track when comparing
the two engines numerically.

For the **closed** parity work, see
[`reports/2026-05-10-parity-report.md`](reports/2026-05-10-parity-report.md)
and the project memory entry `pxd001819-fdr-parity-achieved`.

---

## 1. SpecEValue / GF tails (still not bit-identical to Java)

`gf_java_parity.rs` documents that gates were loosened by ~4 orders of
magnitude because several traced PSMs disagree with Java in **both**
directions — consistent with cumulative differences in the GF DP (node/edge
scores, `add_prob_dist` merges, score thresholds, rounding) rather than one
missing term.

Rust special-cases `score >= max_score` for spectral probability to mirror
Java's tail behavior (`match_engine.rs:419-436`). Identity-level parity
(top-1 PSM agreement) can pass while `ln(SpecProb)` parity does not. The
remaining ~7.8% PXD001819 class flips (per `2026-05-10-parity-report.md`)
are dominated by SpecE-only swaps (72.8% of remaining flips), confirming
this is the load-bearing gap.

**Where:**
- `rust/crates/scoring/src/gf/generating_function.rs` — `compute_inner`,
  `setup_score_threshold`, `f32::from_bits(1)` underflow guard
- `rust/crates/scoring/src/gf/score_dist.rs` — `add_prob_dist`,
  `spectral_probability`
- `rust/crates/scoring/tests/gf_java_parity.rs` — loosened gates and
  hypothesis notes

**To close:** trace one flipping scan end-to-end through `compute_inner`,
diff the per-node `ScoreDist` vs Java's `GeneratingFunctionGroup` output,
identify whether divergence is in (a) per-edge probability accumulation,
(b) score-threshold pruning, or (c) `f32` vs `f64` precision.

---

## 2. E-value (Rust uses an MVP proxy, not Java's suffix-array count) — ITER 1 LANDED 2026-05-10

The queue-derived proxy was replaced by
`SearchIndex.num_distinct_peptides_at_length` (commits `f5f6884`, `a547c39`,
`95fa9bc`, `3e416a3`; report:
[`reports/2026-05-10-evalue-iter1-report.md`](reports/2026-05-10-evalue-iter1-report.md)).

**Status:** Percolator @ 1% FDR rose from 14,798 → 14,850 (+52). EValue
column is now strictly better than the proxy (finite-EValue PSM count
15,929 → 27,779; median value is more conservative).

**Residual gap:** median ratio Rust/Java = **0.0368**, 0% of full-match
PSMs within ±5%. Two open follow-ups:

1. **Mod-aware distinct counting.** Java's `PeptideEnumerator.getNumDistinctPeptides`
   likely counts `PEPTIDE` and `PEPTID+15.99M` as different distinct
   peptides. Rust currently dedupes by bare residue bytes. Fix: extend the
   seen-set key to include the modified peptide form. To investigate:
   `grep -n "getNumDistinctPeptides" src/main/java/edu/ucsd/msjava/msdbsearch/PeptideEnumerator.java`.

2. **Wall regression: 5-6m → 9m17s.** Suspect HashSet<Vec<u8>> allocator
   pressure (~5-10M Vec<u8> per enumerate_candidates pass). Mitigation:
   hash residue bytes to `u64` (FxHash / xxHash) and dedupe by `u64`
   instead of `Vec<u8>`. Saves ~360 MB transient allocation + reduces
   ~10M heap allocations.

---

## 3. Spectrum deconvolution

Java applies `spec.getDeconvolutedSpectrum(...)` when
`scorer.applyDeconvolution()` is true (`NewScoredSpectrum.java`). The Rust
`Param` carries `apply_deconvolution`, but `ScoredSpectrum::new` does **not**
branch on it — so for instrument/protocol combinations where Java turns
deconvolution on, Rust remains on the raw (filtered, ranked) spectrum path.

**Where:**
- `rust/crates/scoring/src/scoring/scored_spectrum.rs` — `ScoredSpectrum::new`
- `rust/crates/scoring/src/param_model.rs` — `Param.apply_deconvolution` (carried but unused)
- Java reference: `NewScoredSpectrum.java:54-55`
  (`if (scorer.applyDeconvolution()) spec = spec.getDeconvolutedSpectrum(...)`)

**To close:** port `Spectrum.getDeconvolutedSpectrum` (charge-state
deconvolution) into the model crate, branch in `ScoredSpectrum::new` when
`param.apply_deconvolution` is true. PXD001819 (LowRes CID, deconv off) does
not exercise this; Astral and HCD HighRes runs typically do.

---

## 4. Top-N semantics after SpecEValue

Rust's `TopNQueue` orders by best `spec_e_value`, then `score`
(`rust/crates/search/src/psm.rs`). Java's hot path orders by `score` into
`PriorityQueue<DatabaseMatch>` and only later reorganizes by spectral
probability / SpecEValue when merging
(`generateSpecIndexDBMatchMap`).

For `top_n > 1`, the **multiset of retained PSMs can differ** even when
top-1 agrees: Rust may evict a higher-`score` candidate that Java would have
kept (because Rust's eviction key is SpecEValue, set later in the pipeline).

PXD001819 currently runs with `top_n = 1`, so this gap is dormant in the
benchmark. Becomes load-bearing the moment we want top-N PSMs for downstream
analysis.

**Where:**
- `rust/crates/search/src/psm.rs` — `PsmMatch::cmp` (line ~120) and
  `TopNQueue::push`
- Java reference: `DBScanner.scan` insertion path + `generateSpecIndexDBMatchMap`

**To close:** for `top_n > 1`, switch `TopNQueue` ordering to `score` first
during insertion, run SpecEValue computation across the retained set, then
re-rank by SpecEValue. Adds a re-sort after Phase 6 but matches Java's
semantics exactly.

---

## 5. Candidate generation architecture

Rust enumerates peptides by walking proteins (+ Met-cleaved branch) with
`expand_mod_combinations`; Java's hot path walks the suffix array via
peptide mass indexing. Comments in `candidate_gen.rs` claim equivalent
candidate sets modulo iteration order, but anything that depends on
**search order, early cutoffs, or memory limits** could still diverge
outside the parity fixtures.

Specifically known: `expand_mod_combinations` only emits
`ModLocation::Anywhere` variants; N-term, C-term, ProtNTerm, ProtCTerm
mod locations are **not expanded** (also tracked as Phase 6 root-cause #2
in `project_phase6_parity_root_causes.md`).

**Where:**
- `rust/crates/search/src/candidate_gen.rs` — `enumerate_candidates`,
  `expand_mod_combinations`
- Java reference: `CandidatePeptideGrid.java`,
  `CandidatePeptideGridConsideringMetCleavage.java`

**To close:** extend `expand_mod_combinations` to read all four
`ModLocation` contexts. Validation: parity test on a fixture that includes
N-term modifications (e.g., Acetyl-ProtN-term in `mods.txt`).

---

## 6. Stale top-of-file comment in `candidate_gen.rs` — FIXED 2026-05-10

The file's module doc said "NO variable-mod expansion (Task 4)" but
`expand_mod_combinations` has shipped. Updated to describe the actual
state (variable-mod expansion implemented for `ModLocation::Anywhere`,
with terminal contexts tracked as item #5 above).

---

## Code-quality notes (Rust)

**Good:**
- Heavy cross-links to Java line-level behavior make maintenance easier
  (`match_engine.rs`, `scored_spectrum.rs`, `gf_java_parity.rs`).

**Watch:**
- `eprintln!` on every multi-spectrum run (yield/progress lines in
  `match_engine.rs`) may be noisy for library/CLI UX. Consider a
  logger or a feature flag (`--quiet` already exists?) when promoting
  msgf-rust beyond benchmark use.
- Per-spectrum `par_iter()` is a reasonable analogue of Java worker
  threads but can change floating-point accumulation order wherever
  order-sensitive ops appear. Most scoring here is per-candidate /
  per-spectrum, so the risk is modest — but `ln(SpecProb)` tail
  divergence (item #1 above) could be partly attributable to
  reduction-order effects from Rayon. Worth an explicit single-thread
  vs multi-thread parity test if item #1 starts moving.

## Verification

A full `cargo test -p search` run is heavy (BSA matching tests take
minutes). The parity fixtures (`gf_java_parity.rs`, `match_engine_bsa.rs`,
`match_engine_java_parity.rs`) are the right gates to consult when any
of the items above are touched. For release-mode parity claims,
run `cargo test -p search --release` so floating-point behavior matches
the binary that ships.

## Bottom line

For **discovery-style parity** against Java on the pinned BSA + test.mgf
+ PXD001819 setup, the Rust code is thoughtfully aligned with
`DBScanner` + `NewScoredSpectrum` + `NewRankScorer`, including the
2026-05-09 fixes (cleavage credit, charge selection, nominal buckets,
GF window, terminal flags) and the 2026-05-10 theo_mz formula fix.

For **numerical parity** (SpecProb tails, E-value, deconvolution-on
protocols, n > 1 ranking), treat Java as reference and keep the gaps
above as **explicit follow-ups, not regressions**, in the parity
checklist.
