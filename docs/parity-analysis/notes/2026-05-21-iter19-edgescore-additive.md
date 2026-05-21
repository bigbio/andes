# Iter19: additive EdgeScore PIN column — SAFE but FLAT

_2026-05-21. Added `EdgeScore` as a NEW PIN column (per-bond DBScanScorer edge sum, RawScore unchanged). Astral 1% FDR: 26,333 (vs iter16 baseline 26,432, Δ -99 within noise). Confirms the n=8 audit rule for additive features: safe, no regression — but doesn't close the 26% gap to Java._

## Hypothesis & implementation

Per n=8, the only safe pattern is ADDITIVE features. The DBScanScorer edge sum is computed per-bond exactly as in Java (mirrors `getEdgeScoreInt` via existing `ScoredSpectrum::edge_score`), but emitted as a SEPARATE PIN column called `EdgeScore` instead of being blended into RawScore (iter17/18 showed that blending regresses -30%).

Implementation:
- `psm_edge_score()` in `scoring/scoring/psm_score.rs`: mirrors Java's reverse/forward edge loop (fromIndex=1, toIndex=n+1).
- `PsmFeatures::edge_score: i32` field.
- `compute_psm_features` now takes `charge: u8`; calls `psm_edge_score`.
- PIN writer emits `EdgeScore` column between `matchedIonRatio` and `Peptide`.
- Schema parity tests updated to accept "Java header + 1 (EdgeScore)".

## Bench results

| Metric | iter16 (baseline) | iter19 (+EdgeScore) | Δ |
|---|---:|---:|---:|
| Targets | 92,977 | 92,767 | -210 (flat) |
| Decoys | 56,570 | 56,499 | -71 (flat) |
| **1% FDR** | **26,432** | **26,333** | **-99 (noise)** |
| 5% FDR | 31,783 | 31,697 | -86 (noise) |
| T/D ratio | 1.643 | 1.642 | preserved ✓ |

EdgeScore distribution across 149,266 PSM rows:
- mean: +61
- min: -50
- max: +919

Wide range with Java-like positive mean (matches HCD/Trypsin expectation per the BSA test fixture audit).

## Why it's flat

Percolator's existing 14 feature columns (RawScore, DeNovoScore, lnSpecEValue, NumMatchedMainIons, longest_b, longest_y, longest_y_pct, ExplainedIonCurrentRatio, NTermIonCurrentRatio, CTermIonCurrentRatio, MS2IonCurrent, MeanErrorTop7, StdevErrorTop7, lnDeltaSpecEValue, matchedIonRatio) already carry enough information to extract the discriminative signal that EdgeScore could add. EdgeScore is correlated with NumMatchedMainIons + lnSpecEValue (both shift when ion existence matches the partition's IES tables).

Same explanation as iter14's MS2IonCurrent fix: "coherent multi-feature fixes at the source are SAFE but not POWERFUL when Percolator is already working around the distortion via other features."

## What this means for closing the 26% gap

EdgeScore was the last additive Java-parity feature I could think of. Adding it doesn't help, which means the remaining gap is NOT at the PIN-feature level. It's at one of:

1. **Top-1 candidate selection** (25% label-flip rate per pin-diff-harness). Different top-1 PSM → different feature row → can't be fixed by adding more features. Need to fix the candidate scoring itself, but iter17/18 proved adding edges to RawScore regresses.

2. **`num_distinct_peptides_at_length`** denominator in lnSpecEValue. Per [[rust-java-known-divergences]]. If Rust counts differently from Java, spec_e_value range shifts.

3. **Candidate enumeration**: mod expansion, missed cleavages, charge handling edge cases.

4. **Peak preprocessing**: deconvolution, ranking, filtering before scoring.

## Decision

- **Ship iter19 as a PR.** No regression, adds Java-comparable per-PSM diagnostic information that's useful for future trace investigations, and the additive feature pattern is the only one with shipped wins.
- **Next gap-close target**: investigate the 25% label-flip rate (option 2 from earlier — diagnose Astral RawScore distribution histograms vs Java).

## Commit

- `d8a8e66f` feat(scoring): additive EdgeScore PIN column (iter19, n=8 audit safe pattern)
