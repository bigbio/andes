# Iter18: atomic Java-mirror (R-3 + C-5b + units + edge) FAILS

_2026-05-21. The compensating-mistakes hypothesis from [[2026-05-20-iter17-edge-score-regresses]] does not hold. Applying all 4 reverted audit-tier fixes simultaneously regresses Astral 1% FDR by 31% (18,157 vs iter16's 26,432). Branch iter18-atomic-mirror discarded; rust-implement HEAD remains `a4e01f8b`._

## Bench summary

| Iter | Fixes applied | T/D | 1% FDR | 5% FDR |
|---|---|---:|---:|---:|
| iter16 (baseline) | C-4 + HIGH-2 + MS2IonCurrent only | 92,977 / 56,570 | **26,432** | ~28K |
| iter17 (single fix) | + edge-scoring (off-by-one) | 82,614 / 58,358 | 18,494 | 20,981 |
| iter17b (off-by-one fixed) | + edge-scoring (correct) | 82,493 / 58,322 | 18,449 | 20,927 |
| **iter18 (atomic)** | + R-3 + C-5b + units + edge | 82,331 / 58,168 | **18,157** | 20,699 |

Delta iter17 → iter18: only -292 PSMs from adding R-3 (-3,093 solo), C-5b (-602 solo), units (-479 solo) on top of edge-scoring. The hypothesized compensation never materialized — the three reverts barely moved anything once edge-scoring was already in.

## Why the hypothesis failed

The "compensating mistakes" theory held that Percolator's weights are miscalibrated for Rust's current distributions, and applying ALL Java-faithful fixes simultaneously would let Percolator re-learn coherent weights. The data refute this:

- T/D ratio in iter18 (1.42) ≈ T/D in iter17 (1.42), both far from iter16's 1.65. The label-flip behaviour of edge scoring is the dominant determinant; the other features can't rescue lost targets.
- The 4 reverted fixes are not compensating each other — each fix moves a feature toward Java's value INDEPENDENTLY, and they remain locally negative for Percolator together too.

## n=8 audit pattern (updated)

| Fix | Type | Astral 1% FDR delta solo | Atomic with edge? |
|---|---|---:|---|
| C-4 (enzN/enzC/enzInt) | ADD new dimensions | +1,718 | (shipped) |
| HIGH-2 (E-value denom) | MODIFY existing | +498 | (shipped) |
| MS2IonCurrent denom | MODIFY existing | flat (-/+60) | (shipped) |
| GF SinkUnreachable retry | MODIFY existing | flat | (shipped) |
| R-3 (minDeNovoScore filter) | MODIFY existing | -3,093 | tested in iter18 |
| C-5b (longest_y_pct denom) | MODIFY existing | -602 | tested in iter18 |
| Units fix (Da→ppm) | MODIFY existing | -479 | tested in iter18 |
| edge-scoring (score_psm) | MODIFY existing | -8,000 | tested in iter18 |

Sum of solo deltas for the 4 reverted fixes: **-12,174 PSMs** (worst case if independent and negative).
iter18 actual delta vs iter16: **-8,275 PSMs**.

The actual atomic delta (-8,275) is LESS negative than the worst-case sum (-12,174), so there IS some compensation. But not enough to make the atomic path positive. The dominant edge-scoring regression overshadows everything.

Rule reinforced: **NEVER modify existing PIN feature distributions piecewise OR atomically without simultaneously fixing the upstream score_psm AND something else that compensates the lost discriminator information.**

## What this implies for closing the 26% gap

Adding edge scoring to `score_psm` makes Rust's RawScore "more Java-like" but Rust's downstream features (RawScore distribution shape, DeNovoScore range, lnSpecEValue range) end up MORE harmful for Percolator's discrimination. The shift to Java-like scores DESTROYS Percolator's ability to distinguish targets from decoys on Rust's pin.

Two viable paths remain:

1. **Additive EdgeScore PIN column** (haven't tried): emit `EdgeScore = sum(per-bond edge_score)` as a NEW PIN column alongside the unchanged RawScore. Percolator learns weights on the augmented feature space. Per the n=8 pattern, additive fixes (like C-4) are the only ones that have shipped wins.

2. **Native MS-GF+ SpecEValue-only FDR** (clean reset): skip Percolator entirely, do FDR on raw `spec_e_value` like the original MS-GF+ tool. Removes the Percolator-recalibration trap. Rust's score_psm would benefit from edges without breaking anything downstream.

3. **(Speculative) Joint Percolator training**: train Percolator on Java's pin then re-use weights for Rust's pin. Not how Percolator works natively; would need custom infrastructure.

4. **(Deferred) Atomic-mirror with FAR more changes**: would need to also fix MeanRelErrorTop7/StdevRelErrorTop7 (large divergence in diff harness), NumMatchedMainIons (+3 Rust vs Java), longest_b (+2), ExplainedIonCurrentRatio (compressed). Multi-week effort with high risk.

## Recommended next-step

Option 1 (Additive EdgeScore column) is lowest-risk, additive, and directly testable in one bench. Should be tried before Option 2 (which is a bigger pipeline change).

## Branch discarded

`iter18-atomic-mirror` (5 commits) is kept locally as a record but not pushed. Rust-implement HEAD remains `a4e01f8b` (iter14 stable baseline 26,461 PSMs).
