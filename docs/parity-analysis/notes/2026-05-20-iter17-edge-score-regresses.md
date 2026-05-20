# Iter17 + iter17b: edge scoring in score_psm REGRESSED Astral by -30% (REVERTED)

_2026-05-20. The edge-score port from [[2026-05-20-edge-score-fixed]] regressed Astral 1% FDR from 26,432 (iter16) → 18,449. The off-by-one fix didn't recover it. Reverted; n=7 confirmation of the "piecewise modifies distribution → Percolator regression" pattern._

## The hypothesis

The audit in [[2026-05-20-edge-score-fixed]] established:
1. Java's `DBScanScorer.getScore` is node + cleavage + edge.
2. Rust's `score_psm` was node + cleavage only (missing the per-bond IES + error_score that the GF DP graph already adds).
3. Per-bucket avg edge scores match Java HCD (idx0 -4.00, idx1 -1.00, idx2 -1.00, idx3 ~+1.00).

So adding the DBScanScorer edge loop to `score_psm` should bring RawScore in line with Java's, and Percolator @ 1% FDR should approach Java's 35,818 (we were at 26,432).

## What actually happened

| Iter | Description | T/D | 1% FDR | 5% FDR |
|---|---|---:|---:|---:|
| iter12 | C-4 baseline (rust-implement HEAD pre-fix) | 92,909 / 56,442 | 26,401 | ~28K |
| iter13–16 | small fixes (units, MS2IonCurrent, GF retry) | 92,946 / 56,543 | 26,401–26,461 | flat |
| iter17  | edge-scoring fix 2d63ff84 (off-by-one BUG in reverse loop) | 82,614 / 58,358 | **18,494** | 20,981 |
| iter17b | + off-by-one fix 683e8796 (reverse loop now `(1..n).rev()`) | 82,493 / 58,322 | **18,449** | 20,927 |

The off-by-one fix barely moved the needle (-45 PSMs). The regression is structural.

## Why edge scoring breaks Percolator @ 1% FDR

Diagnostics:

- **Targets dropped 92.9K → 82.6K (-11%)**; decoys gained 56.4K → 58.3K (+3%). The top-1 candidate per spectrum CHANGED — fewer real targets, more decoy noise. Edge scoring re-ranks candidates differently than Percolator's learned weights expect.
- The pattern matches **R-3** (charge-ambiguous spectra), **C-5b** (NumMatchedMainIons), **units fix** (Da→ppm in MeanErrorTop7/StdevErrorTop7), and now **edge scoring**: all four are "Java-faithful piecewise fixes that MODIFY existing column distributions" → all four regress.

Speculative mechanism: when only one engine layer is changed to be "more Java-like", the downstream Percolator weights (which were trained on the pre-fix distribution) no longer discriminate as well — even if the change moves individual PSMs toward Java's score values. The RawScore distribution's shape matters as much as its absolute values.

## n=7 audit pattern (vs the n=6 we had)

ADDITIVE fixes (Rust didn't compute the feature; new dimension in PIN):
- C-4 enzN/enzC/enzInt: +1,718 PSMs ✓

DISTRIBUTION-MODIFYING fixes (Rust computed feature; values change):
- R-3 (charge-ambiguous): -3,093 PSMs ❌
- C-5b (NumMatchedMainIons): -602 PSMs ❌
- HIGH-2 (E-value denominator): +498 PSMs ✓ (exception — kept)
- units fix (Da→ppm): -479 PSMs ❌
- MS2IonCurrent denom: ~flat (-/+60) ≈
- GF SinkUnreachable retry: ~flat (within noise) ≈
- **edge-score (this iter): -8,000 PSMs** ❌

Rule reinforced: **piecewise Java-mirror fixes that modify existing per-PSM feature values regress Percolator. Only ADDITIVE (new dimensions) or noise-level (no Percolator-visible change) fixes are safe one-at-a-time.**

## The fundamental tension

Java uses node + cleavage + edge scoring everywhere; it gets 35,818 PSMs @ 1% FDR. Rust is missing edges in score_psm; it gets 26,432. Adding edges to Rust's score_psm to match Java drops Rust to 18K. So the "Java-faithful" fix is locally regressing while we know Java's behavior produces the higher count.

This implies the gap isn't a single-fix away — multiple things need to change simultaneously. A piecewise approach won't work because each fix breaks the local Percolator calibration.

## What to do instead

1. **Don't ship piecewise edge-scoring**. Filed under [[piecewise-alignment-doesnt-work]].
2. To close the 26% gap, options ordered by likelihood:
   - **Re-train Percolator weights** by feeding it both Java's PIN and Rust's PIN simultaneously; learn the right scale. Not how Percolator actually works, but worth thinking about whether feature reweighting could close the gap.
   - **Mirror Java's full pipeline (score_psm AND PIN features AND GF DP) atomically** in one branch, ship as a single change. High risk, high investment. Would need its own bisect harness.
   - **Switch to a non-Percolator FDR scheme** (e.g., MS-GF+'s native SpecEValue-only FDR). Cleaner separation; Rust's score_psm would match Java better one-shot. Risk: loses the ML uplift Percolator provides.
   - **Find an ADDITIVE feature** that captures the edge-score information without modifying RawScore. E.g., add `EdgeScore` as a NEW PIN column alongside RawScore. Percolator can learn weights for it on top of existing features. Low risk per the n=7 pattern.

Option 4 (additive `EdgeScore` column) is the highest expected value next step.

## Reverted commits

- `2d63ff84` — fix(scoring): add DBScanScorer-style edge scoring to score_psm
- `683e8796` — fix(scoring): off-by-one in score_psm edge loop — start at i=1, not i=0

Reverts in `1494cd6f` (main fix) and `addc2927` (off-by-one).

## Preserved (not reverted)

- `a25ba894` — docs(parity): edge-score divergence localized
- `aa93abe0` — docs(parity): audit findings for edge-score re-fix
- `14de1971` — diag(trace): selected-partition + dump-all flag for msgf-trace
- `18d2de48` — docs(parity): score_psm divergence — Rust scores Java target peptides 20+ pts lower

These remain useful empirical artifacts.
