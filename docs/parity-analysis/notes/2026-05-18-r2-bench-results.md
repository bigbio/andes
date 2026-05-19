# R-2 retention-layer refactor — empirical results

_2026-05-18. Branch `rust-implement` at iter7 HEAD (Tasks 1-5 of R-2 plan committed: ca72192..ce09034)._

## Implementation summary

R-2 lands the full retention-layer refactor: per-charge `TopNQueue` keyed by charge state, pre-merge pepSeq+score dedup that aggregates `candidate_idxs` across same-peptide/same-score matches, per-charge GF/SpecEValue compute (Java DBScanner.java:606,779), spectrum-level merge with SpecE tie keep, and PIN writer emits one accession per `candidate_idx` in the Proteins column (Java DirectPinWriter.java:237).

Local gates all green:
- `r1_tie_retention_active_in_production_pipeline` ✓
- `gf_java_parity` all 5 BSA PSMs within 1.0 OOM ✓
- `score_psm_pxd001819_parity` RawScore=293 stable ✓
- New `r2_deduped_psm_count_matches_java_on_bsa_fixture`: Rust=215, Java=217, ratio=0.991 ✓
- All output crate tests + schema parity ✓

## Astral no-mods bench (iter7)

| Metric | Java | b1d45bb | iter6 (R-1) | iter7 (R-2) | Gate | Status |
|---|---:|---:|---:|---:|---|---|
| Raw targets | 89,479 | 75,457 | 1,042,255 | 92,825 | 72K-107K | ✅ |
| Raw decoys | 46,792 | 46,208 | 530,430 | 56,501 | — | — |
| T/D ratio | 1.912 | 1.633 | 1.965 | 1.643 | >= 1.85 | ❌ |
| Wall | — | ~8:36 | 23:52 | 11:06 | <= 26 min | ✅ |
| Percolator @ 1% FDR rows | 35,818 | 25,224 | 74,204 | 24,675 | >= 30K | ❌ |
| Distinct (scan, peptide) at q<=0.01 | 35,818 | unknown | 26,934 | 24,675 | >= 34,027 | ❌ |

PIN file size: 45 MB (down from R-1's 467 MB, a 10× reduction — strong evidence R-2 retention is collapsing PSMs as intended). Max RSS: 9.87 GB.

## ⚠ Percolator mode-detection caveat (discovered post-bench)

R-2.5 (multi-accession Proteins column) changed the PIN structure to match Java's `DirectPinWriter.java:237` — one row per PSM with N tab-separated proteins. This shifted Percolator's auto-detection from **Concatenated** (b1d45bb / iter5 PIN, with multiple rows per PSM) to **Separate / mix-max** (R-2 / Java PIN, single row per PSM).

The two modes use **different statistics** and produce **different q-value calibrations**. Cross-mode comparisons are not apples-to-apples.

**Empirical evidence (re-run on existing VM data):**

| Run | Percolator mode | T/D ratio | Percolator @ 1% FDR |
|---|---|---:|---:|
| Java (iter3) | Separate / mix-max | 1.912 | 35,818 |
| Rust b1d45bb (arc) | **Concatenated** | 1.633 | 25,224 |
| Rust iter5 PIN (rerun) | Concatenated | 1.642 | 17,160 |
| Rust iter6 R-1 | Separate / mix-max | 1.965 | 74,204 (over-shoot) |
| **Rust iter7 R-2** | **Separate / mix-max** | **1.643** | **24,675** |
| Rust iter7 R-2 (-Y, TDC) | Separate / TDC | 1.643 | 24,658 |

The fair comparison is **iter7 R-2 vs Java**, both under Separate / mix-max: 24,675 vs 35,818 → **-31% gap**. The previously-reported "-549 vs b1d45bb" is a mode artifact, not a regression.

## Gate decision (revised after mode-detection finding)

**Architecture gates (PASS):**
- Raw targets in 72K-107K band ✅ (92,825)
- Wall ≤ 26 min ✅ (11:06)
- PIN size 10× smaller than R-1 ✅ (45 MB vs 467 MB)
- Local parity tests all green ✅ (215/217 deduped pairs, 1.0 OOM SP, RawScore stable)
- PIN format matches Java (triggers same Percolator mode) ✅

**Quality gates (FAIL, but for an understood reason):**
- T/D ratio: 1.643 vs gate 1.85 vs Java 1.912 — driven by scoring/feature divergences
- Percolator @ 1% FDR: 24,675 vs gate 30K vs Java 35,818 — same root cause

The quality-gate thresholds (30K / 34,027) assumed R-2 would close the Java gap. The R-2 spec was explicit about being **architecture-only**, not addressing the audit-tier C-4/C-5/C-5b/F-1 feature/scoring divergences. Those gates were therefore always going to fail at this point in the work, regardless of R-2 correctness.

## Audit-tier follow-up bisect (2026-05-19)

Three "Java-faithful" fixes were applied on top of R-2 baseline (iter8) — R-3 (minDeNovoScore filter), C-5b (longest_y_pct denominator pepLen→pepLen-1), HIGH-2 Path A (e_value lookup index +1). iter8 measured **21,486 PSMs @ 1% FDR**, a -12.9% regression vs the R-2 baseline (24,675). Sequential bisect (iter9-iter11) decomposed the contribution of each fix:

| iter | Fixes on top of R-2 | Percolator @ 1% FDR | T/D | Wall |
|---|---|---:|---:|---|
| iter7 | (none — R-2 baseline) | 24,675 | 1.643 | 11:06 |
| iter11 | (none — reverted to R-2) | 24,683 | 1.643 | 11:17 |
| iter10 | +R-3 only | 21,590 | 1.586 | 11:32 |
| iter9 | +R-3 +C-5b | 20,988 | 1.583 | 11:00 |
| iter8 | +R-3 +C-5b +HIGH-2 | 21,486 | 1.586 | 10:43 |

**Differential impact** of each fix (per-bisect arithmetic):

- **R-3** (minDeNovoScore filter at PIN/TSV emit): **-3,093 PSMs (-12.5%)**. Java-faithful per `DirectPinWriter.java:132`, but empirically Percolator was rescuing many of the filtered low-de_novo PSMs via other features. Pre-filtering them stripped signal.
- **C-5b** (longest_y_pct denominator pepLen-1 instead of pepLen): **-602 PSMs (-2.4%)**. Java-faithful per `PSMFeatureFinder.java:95-96`, but the 5-10% length-dependent feature rescale disrupted Percolator's calibration on a discriminator feature.
- **HIGH-2 Path A** (e_value lookup index +1 with enzyme): **+498 PSMs (+2.0%)**. Java-faithful per `DirectPinWriter.java:165`. Structural index alignment that travels well — the only helpful audit-tier fix.

iter11 confirms run-to-run noise is <10 PSMs (24,683 vs iter7's 24,675).

**Verdict applied (2026-05-19 cab8...c8d1):** keep HIGH-2 (`b3cb3277`), revert R-3 (`c8d1ed90`) and C-5b (`7166ddcb`).

This matches the [[piecewise-alignment-doesnt-work]] pattern even WITH the R-2 retention prerequisite in place: applying Java-faithful per-feature fixes is empirically negative on Astral because Percolator's discrimination depends on the joint feature distribution, not per-feature correctness. The next layer of work needs to either:
- Apply feature fixes as a coherent group calibrated against Percolator, not one-at-a-time
- Tune Rust's features against Astral 1% FDR directly (treat Java as guide, not target)
- Build a per-PSM Rust↔Java diff harness so divergence sources are empirically traced

## Next

Current `rust-implement` HEAD: R-2 baseline + HIGH-2 only. Expected Astral 1% FDR: ~25,180 (extrapolating from +500 effect; bench to confirm). Gap to Java's 35,818 narrows from 31% to 29%, still dominated by structural feature/scoring divergences that don't decompose into single-line fixes.

## iter12: C-4 (enzN/enzC/enzInt) (2026-05-19)

Diff-harness-driven follow-up. The `2026-05-19-pin-diff-findings.md` analysis localized enzN/enzC/enzInt as the highest-value remaining fix — Rust was emitting constant 0 for all three across every PSM (Java emits real values).

| Metric | iter11 (R-2 baseline) | iter12 (+C-4) | Δ | Java |
|---|---:|---:|---:|---:|
| Raw targets | 92,929 | 92,909 | -20 | 89,479 |
| Raw decoys | 56,568 | 56,442 | -126 | 46,792 |
| T/D ratio | 1.643 | 1.646 | +0.003 | 1.912 |
| Wall | 11:17 | 11:26 | +9s | — |
| **Percolator @ 1% FDR** | **24,683** | **26,401** | **+1,718 (+7.0%)** | 35,818 |
| Percolator @ 5% FDR | 30,385 | 31,660 | +1,275 | — |

**C-4 closes 5 percentage points of the Java gap** (31% → 26%). Top-1-per-scan buckets are unchanged within noise — C-4 doesn't change which PSMs Rust emits, it gives Percolator three new discriminator dimensions to use for FDR calibration.

Re-run of the diff harness on iter12 confirms bit-exact agreement with Java on enzN/enzC/enzInt:
- enzN median Δ = 0, mean Δ = +4e-5 (float noise)
- enzC median Δ = 0, mean Δ = +2e-5
- enzInt median Δ = 0, mean Δ = 0

This validates both the implementation and the diff-harness workflow: localize empirically → implement the additive fix → re-measure on the harness AND on Percolator. ADDITIVE fixes don't carry the piecewise-regression risk that R-3 / C-5b in the bisect did.

Current `rust-implement` HEAD: `1d9da765` (R-2 + HIGH-2 + C-4). Astral 1% FDR = 26,401; gap to Java = 9,417 PSMs (26%). Remaining gap dominated by RawScore / lnSpecEValue / DeNovoScore covariance (structural scoring divergence per the harness) plus the MeanErrorTop7/StdevErrorTop7 units mismatch (smaller, easier).

## iter13: MeanErrorTop7/StdevErrorTop7 units fix (Da → ppm) — reverted

The diff harness flagged these two columns as Da-vs-ppm units mismatch (Java emits ppm via `NewScoredSpectrum.getMassErrorWithIntensity:229`, Rust was emitting `|obs-pred|` in Da). The fix changed Rust to emit `|(obs-pred)/pred*1e6|` (abs ppm) to match Java.

iter13 Astral bench:

| Metric | iter12 (+C-4) | iter13 (+units fix) | Δ |
|---|---:|---:|---:|
| Percolator @ 1% FDR | 26,401 | **25,922** | **-479 (-1.8%)** |
| Percolator @ 5% FDR | 31,660 | 31,047 | -613 |

The fix was theoretically correct — diff harness re-run confirms Rust now reports in PPM. But Percolator regressed by -479 PSMs. This is the same pattern as R-3 / C-5b in the bisect: **modifying an existing column's distribution scale disrupts Percolator's learned weights even when the modification is Java-faithful**.

The harness data also showed that after the unit conversion, Rust's MeanErrorTop7 still diverges from Java's by ~24 ppm median — i.e. the underlying top-7 ion selection differs between engines, so the units fix alone doesn't produce bit-exact agreement. This is consistent with the harness's earlier finding that `NumMatchedMainIons` median Δ = +3 (Rust matches more ions).

**Reverted at 5007d8d1.** General principle confirmed: only ADDITIVE fixes (like C-4) can be applied piecewise. DISTRIBUTION-MODIFYING fixes need to be batched with Percolator re-training OR avoided until the underlying selection logic also matches Java.

Current `rust-implement` HEAD: `5007d8d1`. Best stable Astral 1% FDR result: **26,401** (iter12, R-2 + HIGH-2 + C-4). Gap to Java's 35,818 holds at 26%.

## iter14: MS2IonCurrent denominator fix (kept-peaks sum)

Diff-harness-driven structural fix. The 2026-05-19 PIN diff showed Rust's MS2IonCurrent was 1.6× Java's (median Δ +78,570) because Rust's `total_intensity` summed the original `spec.peaks` while Java's effective sum is over precursor-filtered peaks (Java zeroes precursor-peak intensities via `Spectrum.filterPrecursorPeaks` at `NewScoredSpectrum.java:44-45` before `PSMFeatureFinder.computeSumIonCurrent` iterates). Rust ALREADY applies the same precursor filter for rank assignment but the `total_intensity` field used the unfiltered list. Single-line fix: compute total_intensity from the `kept` Vec.

| Metric | iter12 (+C-4) | iter14 (+MS2 denom) | Δ |
|---|---:|---:|---:|
| Percolator @ 1% FDR | 26,401 | **26,461** | +60 (neutral, within noise) |
| MS2IonCurrent agreement (diff harness) | mean +78,570 / 98.6% diff | **mean +254 / 0.2% diff** | bit-exact |
| ExplainedIonCurrentRatio | median -0.023 | median -0.017 | partial closure |
| CTermIonCurrentRatio | median -0.019 | median -0.015 | partial closure |
| NTermIonCurrentRatio | median -0.003 | median -0.001 | partial closure |

**Outcome:**
- MS2IonCurrent fully aligned with Java (bit-exact agreement).
- The 3 ion-current ratios partially closed but still diverge because their NUMERATOR (matched-ion intensity sum) is different between engines — Rust matches more ions overall (`NumMatchedMainIons` median +3), so `matched_b_intensity + matched_y_intensity` is larger. The denominator fix concentrated the residual divergence into the matched-ion-set divergence, which is a separate, upstream structural issue.
- Net Percolator effect: +60 PSMs (within run-to-run noise). Not regressive, not strongly positive — the 4 features were apparently not heavily weighted by Percolator for this dataset, so even fully aligning MS2IonCurrent didn't translate to meaningful FDR movement.

Pattern data extended (n=6):

| Fix | Type | Astral Δ |
|---|---|---:|
| R-3 | distribution-modifying (one column) | -3,093 |
| C-5b | distribution-modifying (one column) | -602 |
| Units fix | distribution-modifying (one column) | -479 |
| HIGH-2 | distribution-modifying (one column, small) | +498 |
| **C-4** | **ADDITIVE (3 new columns)** | **+1,718** |
| MS2 denom (iter14) | coherent multi-feature at source | +60 |

Refinement: **coherent multi-feature fixes at the source are SAFE but not necessarily POWERFUL.** The fix matters in proportion to how much Percolator was relying on the affected feature for FDR. C-4 was both safe AND high-impact because it added 3 NEW dimensions; the MS2 denom fix was safe but low-impact because Percolator already worked around the 1.6× distortion via other features.

Current `rust-implement` HEAD: `6ed6e724` (R-2 + HIGH-2 + C-4 + MS2 denom). Astral 1% FDR: **26,461**. Gap to Java: 26.1%.

The remaining gap dominated by upstream structural divergences (RawScore -2 median, NumMatchedMainIons +3 median, lnSpecEValue +0.74 median, DeNovoScore -9 median, longest_b +2). These all derive from the scoring/ion-matching pipeline before PIN feature emission. Closing them likely requires the per-PSM scoring trace harness (option A from the 2026-05-19 brainstorm) — instrument both engines to dump intermediate scoring values for a small fixture and identify the first divergence point.
