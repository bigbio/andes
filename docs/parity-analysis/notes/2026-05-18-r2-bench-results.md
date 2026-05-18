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

## Gate decision

**2 of 4 main gates pass; 3 of 5 metric gates fail.** Per plan's gate decision: do NOT proceed to Task 7 (PXD + TMT) without addressing the failures.

**What R-2 fixed (positive signals):**
- R-1's 11.6× raw-target over-shoot is gone (1,042,255 → 92,825, within healthy 72K-107K band)
- Wall time recovered to 11:06 (vs R-1's 23:52 — 2× speedup, beats b1d45bb)
- PIN size dropped 10× (467 MB → 45 MB) — retention is correctly capping per-spectrum candidates

**What R-2 did NOT fix (real signals):**
- Percolator @ 1% FDR: 24,675 — small regression from b1d45bb's 25,224 (-2.2%); well below Java's 35,818 (-31%) and the 30K gate
- T/D ratio: 1.643 vs Java's 1.912 — Percolator can't separate targets from decoys as well as Java
- These point to **PSM quality/discriminability**, not retention. Architecture parity is established; scoring + feature divergences (audit-tier C-4, C-5, C-5b, F-1) still need closing.

## Next

User decision required. Three options:

**Option A — Keep R-2 baseline, continue with feature fixes (recommended).** R-2 correctly establishes per-SpecKey retention parity (Tasks 1-5 are architecturally aligned with Java). The remaining gap is in PSM scoring/feature quality, which is what the divergence-audit C/F items target. Build feature fixes on top of R-2; re-bench when feature work lands.

**Option B — Revert R-2 to b1d45bb baseline.** Strict reading of the plan's "if any gate fails, decide on revert." Loses the architecture work and the wall-time recovery, but matches the original gate's literal interpretation. The 24,675 vs 25,224 Percolator regression (-2.2%) becomes the proximate trigger.

**Option C — Investigate the small regression before deciding.** The 549-PSM drop from b1d45bb (24,675 vs 25,224) is small enough that it could be one specific divergence (e.g. per-charge GF SpecE calibration shifting some borderline PSMs out of 1% FDR). A targeted diff between b1d45bb's and iter7's `astral_iter7.target.psms.txt` at the 1%/1.1% q-value boundary would reveal whether the drop is structural or noise.

Task 7 (PXD + TMT bench) skipped pending decision.
