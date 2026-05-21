# Iter20: feature-counting tolerance fix — +4,650 PSMs @ 1% FDR (+17.7%)

_2026-05-21. Localized that Rust's `compute_psm_features` was using `scorer.param().mme.as_da(p.mz)` = 0.5 Da for HCD_QExactive_Tryp.param — ~50× wider than Java's hardcoded 20 ppm for high-resolution instruments. Fixed; Astral 1% FDR jumped from 26,333 (iter19) to **30,983**, closing the gap to Java from 26% to 13.5% in a single change._

## The bug

Java's `PSMFeatureFinder.java:51-54`:
```java
if (scorer.getSpecDataType().getInstrumentType().isHighResolution())
    mme = new Tolerance(20f, true);    // 20 ppm
else
    mme = new Tolerance(0.5f, false);   // 0.5 Da
```

This hardcoded tolerance is used for feature-counting fragment matches (b/y at charge 1). The param's `mme` value is the coarser binning tolerance used by the rank-distribution scoring tables (appropriate for node-score table lookup), NOT the precise fragment-matching tolerance.

Rust's `match_engine.rs:741` was using `scorer.param().mme.as_da(p.mz)` for the same feature-counting code path. For HCD_QExactive_Tryp.param this gives `mme = Da(0.5)`. At m/z 500, this is 0.5 Da window vs Java's 20 ppm = 0.01 Da window — a **50× wider** matching window.

The wider window in Rust matched noise peaks Java skipped, which:
1. Inflated `NumMatchedMainIons` (+3 vs Java)
2. Lengthened `longest_b` (+2 vs Java)
3. Compressed all intensity ratios (more low-intensity noise peaks → matched-ion intensity sum had lots of low values → ratio smaller)

## The fix (cf287c4a)

- Added `InstrumentType::is_high_resolution()` mirroring Java's `InstrumentType.isHighResolution()` (true for HighRes, TOF, QExactive; false for LowRes).
- `compute_psm_features` now uses hardcoded 20 ppm (high-res) / 0.5 Da (low-res) instead of `param.mme`.
- Edge scoring and node scoring continue to use `param.mme` (Java does the same; bit-exact match preserved).
- Unit test offset: 0.01 Da → 0.0005 Da to fit within the new 20 ppm window.

## Bench results

| Iter | Description | T/D | 1% FDR | 5% FDR |
|---|---|---:|---:|---:|
| iter16 baseline | C-4 + HIGH-2 + MS2IonCurrent | 92,977 / 56,570 | 26,432 | 31,783 |
| iter17 | + edge-scoring in RawScore | 82,614 / 58,358 | 18,494 (regress) | 20,981 |
| iter18 | + R-3 + C-5b + units + edge | 82,331 / 58,168 | 18,157 (regress) | 20,699 |
| iter19 | + additive EdgeScore PIN col | 92,767 / 56,499 | 26,333 (flat) | 31,697 |
| **iter20** | **+ feature-tolerance fix** | **92,994 / 56,467** | **30,983** | **34,307** |
| Java reference | | | 35,818 | — |

**Δ vs iter19: +4,650 PSMs @ 1% FDR (+17.7%), +2,610 @ 5%.**
**Gap to Java: 26% → 13.5%.**

## Per-feature alignment (iter20 vs Java diff harness)

| Feature | iter16 median Δ | iter20 median Δ | Notes |
|---|---:|---:|---|
| NumMatchedMainIons | +3 | -1 | converged; slight overcorrection (Rust matches 1.3 fewer) |
| longest_b | +2 | 0 | bit-exact ✓ |
| longest_y | 0 | 0 | unchanged (already good) |
| ExplainedIonCurrentRatio | -0.017 | -0.026 | slight regression — Rust ratio still under Java |
| NTermIonCurrentRatio | -0.001 | -0.005 | slight regression |
| CTermIonCurrentRatio | -0.015 | -0.018 | flat |
| RawScore | -2 | -2 | unchanged (score_psm uses param.mme = 0.5Da, same as Java) |
| DeNovoScore | -10 | -10 | unchanged |
| lnSpecEValue | -0.72 | -0.72 | unchanged |
| MS2IonCurrent | 0 | 0 | bit-exact (iter14 fix preserved) |

Net: feature counts converged dramatically, ratios moderately worse — but Percolator extracts +4,650 PSMs of signal that was previously buried in noise.

## n=8 audit pattern — REFINED to n=9

Prior rule: "MODIFYING-EXISTING-DISTRIBUTION fixes regress; only ADDITIVE fixes ship."

Iter20 is a clear counter-example: it modifies existing feature values (NumMatchedMainIons, longest_b, intensity ratios) and gains +4,650 PSMs.

**Refined rule:**
- Fixes that change **top-1 candidate selection** (= RawScore distribution) regress, because Percolator's weights are calibrated against the current top-1 set. Edge-scoring blended into RawScore (iter17/iter18) fits this.
- Fixes that **clean up noise in feature values without affecting RawScore-based ranking** can improve. The tolerance fix is this category — Rust was finding spurious "matches" that were actively misleading Percolator.

Empirical test: T/D ratio.
- Iter17/iter18 T/D = 1.42 (top-1 selection diverged catastrophically)
- Iter20 T/D = 1.647 (≈ iter16's 1.643; top-1 selection unchanged)
- The T/D ratio's deviation from baseline = early indicator of whether a fix is in the "selection-changing" (dangerous) or "noise-cleanup" (safe) category.

## What's still off in iter20

The residual 13.5% gap appears to come from features still divergent:
- **MeanRelErrorTop7 / StdevRelErrorTop7** — large spread (units fix reverted; may now be safe to re-apply on top of iter20).
- **MeanErrorTop7 / StdevErrorTop7** — units fix again.
- **RawScore / DeNovoScore / lnSpecEValue** — the structural score_psm divergence; edge-scoring would close it but regresses.
- **ExplainedIonCurrentRatio / NTermIonCurrentRatio / CTermIonCurrentRatio** — slight residual divergence; may need a separate audit of Java's intensity-sum computation.

## Next steps

1. **Re-test units fix on top of iter20** (was -479 solo on iter16 baseline; may now be net-positive given the cleaner feature distribution).
2. **Ship iter20 as a PR** — single largest Astral 1% FDR improvement since C-4.
3. **Re-audit intensity ratios** — small residual divergence may have a localized fix.
4. **Defer edge-scoring** — still regresses on top-1 selection; needs a different approach (e.g., additive EdgeScore column from iter19 is on the table).

## Commits

- `d8a8e66f` feat(scoring): additive EdgeScore PIN column (iter19, flat)
- `cf287c4a` fix(features): hardcoded 20ppm/0.5Da feature tolerance like Java PSMFeatureFinder (iter20, +4,650 PSMs)

Branch: `iter19-additive-edge` (now iter20 effectively — should rename to `iter20-feature-tolerance-fix` on PR).
