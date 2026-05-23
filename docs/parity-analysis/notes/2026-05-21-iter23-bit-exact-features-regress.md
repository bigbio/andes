# Iter23: bit-exact feature parity REGRESSES Percolator (REVERTED)

_2026-05-21. The cleanest Java feature alignment in the entire parity-work history caused a -1,404 PSM regression at 1% FDR. Decisive evidence that **per-feature Java parity ≠ Percolator FDR parity**._

## What was done

iter23 extended iter22b's partition-ion-list per-bond loop with per-bond charge-1 max-intensity picks (`best_b_charge1` / `best_y_charge1`), mirroring Java's `getMassErrorWithIntensity` selection exactly. The `matched_ions` buffer was cleared and re-populated from these picks, so:
- NumMatchedMainIons = count of (bond, direction) tuples with ANY charge-1 match
- Top-7 error stats computed from the max-intensity matched peak per (bond, direction)

Java's `PSMFeatureFinder.getMassErrorWithIntensity` (NewScoredSpectrum.java:284) does exactly this.

## Result

**Bit-exact feature parity. -1,404 PSMs @ 1% FDR.**

| Feature | iter22b median Δ | **iter23 median Δ** |
|---|---:|---:|
| NumMatchedMainIons | -1 | **0** ✓ |
| MeanErrorTop7 | -1.47 | **+0.013** ✓ |
| MeanRelErrorTop7 | +1.46 | **-0.016** ✓ |
| StdevErrorTop7 | -0.34 | **+0.003** ✓ |
| StdevRelErrorTop7 | -0.42 | **+0.005** ✓ |
| longest_b/y | 0 | 0 |
| Intensity ratios | bit-exact | bit-exact |

| Iter | 1% FDR | Δ vs iter22b |
|---|---:|---:|
| iter22b baseline | 31,006 | — |
| iter23 (bit-exact) | **29,602** | **-1,404** |

T/D ratio: 1.647 (preserved). Top-1 candidate selection unchanged. The regression is purely from feature distribution changes Percolator's weights can't accommodate.

## What this teaches

Per-feature parity to Java is NOT a viable path to closing the FDR gap. Percolator's discriminative weights are trained on the SHAPE of Rust's feature distributions, not Java's. Bringing features to bit-exact-Java values disrupts that calibration even when individually correct.

The only feature-level fix that has shipped a win is iter20's tolerance fix, and that worked because Rust's pre-fix feature values were carrying **active misinformation** (noise peaks matched in a too-wide window). Removing the noise was a net gain regardless of Percolator calibration.

iter23 — and earlier iter17/iter18 (edge-scoring) — confirm that "Java is right, Rust is wrong" framing is the wrong frame. The correct frame is: **which feature changes give Percolator MORE signal vs which give it equivalent or noisier signal**.

## Implications for closing the remaining 13.4% gap

The 25% top-1 label-flip rate (Java/Rust pick different peptides) is the dominant FDR driver. Closing those flips requires **score_psm-level** convergence — but iter17/iter18 proved that adding edge scoring to score_psm regresses too, because it changes RawScore distribution AND top-1 selection in ways Percolator can't accommodate.

This leaves a narrow window of options:
1. **Native MS-GF+ SpecEValue-only FDR** — skip Percolator entirely; FDR computed from raw spec_e_value. score_psm fixes (edge scoring) would help here because spec_e_value uses GF DP which already includes edges. Cleanest reset.
2. **Re-train Percolator on Java's PIN as a starting calibration** — speculative; needs custom infrastructure.
3. **Find feature changes that ADD signal without disrupting existing distributions** — only iter20 fit this. EdgeScore (iter19) was technically additive but flat because correlated with existing features.

Path 1 (SpecEValue-only FDR) is the most concrete remaining lever. It would change the FDR scheme but unlock the score_psm edge-scoring fix's actual benefit.

## Commits

- `a1eb10bd` fix(features): NumMatchedMainIons + error stats use partition charge-1 max-intensity pick (iter23) — REVERTED
- `46775de4` Revert "fix(features): ..." (iter23 revert)

Stable branch HEAD: `b6b88341` (iter22b — 31,006 PSMs @ 1% FDR).
