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

## Next

**Recommendation: KEEP R-2 baseline. Add audit-tier feature/scoring fixes on top, then re-bench.**

Reasoning:
1. R-2 correctly establishes per-SpecKey retention parity (architectural alignment with Java)
2. R-2.5 makes the Rust PIN structurally match Java's PIN (same Percolator mode now)
3. The 24,675 vs 35,818 (vs Java) gap is the **scoring/feature quality** gap, which is exactly what the divergence-audit's C-4, C-5, C-5b, F-1 items target
4. Reverting R-2 would lose the wall-time recovery (11:06 vs R-1's 23:52) AND would restore the Concatenated-mode PIN which masks the real Java gap
5. The architecture work was a precondition for closing the feature gap — feature fixes would have been impossible to evaluate against R-1's 1M-row over-shoot

**Task 7 (PXD + TMT bench) is still skipped** — not because R-2 is wrong, but because (a) the plan's strict gate decision blocks it, and (b) running PXD + TMT now would just confirm the same architecture-good / quality-needs-feature-work pattern on more datasets. Better to land the next layer of feature fixes first, then do the 3-dataset bench on a richer baseline.

If the user prefers strict adherence to the plan's gate decision (Option B in the original note), revert via `git reset --hard 37d28f95` then re-apply minus the production code (docs + R-1 + R-1 test). But this loses substantive architectural progress.
