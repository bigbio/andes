# Gate run — chimeric NO_RESCORE vs Java: real PSM gains, blocked on speed + TMT

**Date:** 2026-05-29
**Branch:** `feat/chimeric-dda-plus`
**Config:** Rust chimeric, `MSGF_CHIMERIC_NO_RESCORE=1` (Phase-1 emission; the
entrapment work proved it FDR-honest and the residual rescore unnecessary), normal
target+decoy DB, cal=auto. Java = upstream v2024.03.26.

## Result (@1% FDR Percolator PSMs; wall)

| Dataset | Java PSMs | Rust-chim PSMs | ΔPSM | Java wall | Rust-chim wall | speed |
|---|---:|---:|---:|---:|---:|---:|
| PXD001819 | 14,989 | **18,234** | **+21.6%** | 1:20.9 | 1:33.5 | 1.16× slower |
| Astral | 35,818 | **77,287** | **+115.8%** | 6:17.7 | 17:04.1 | 2.71× slower |
| TMT | 10,194 | 9,686 | −5.0% | 3:09.9 | 4:20.5 | 1.37× slower |

Rust-chim PSM counts are **entrapment-validated FDR-honest**
(`2026-05-29-entrapment-fdp-reversal.md`: PXD FDP 0.34%, Astral 0.42%, fair 50/50
noise floor). So the gains are *real* IDs, not inflation.

## Gate status: NOT cleared (2 crisp blockers)

The merge gate = beat Java on **PSMs AND speed, all 3 datasets**.
- **PSMs:** win 2/3 (PXD +21.6%, Astral +116%), **lose TMT (−5%)**.
- **Speed:** **lose 3/3** — wide-window candidate explosion (Astral 17 min vs 6 min).

This is NOT a refutation. The PSM gains are real and large. Two engineering blockers
remain:

### Blocker 1 — Speed (wide-window candidate explosion)
The full-isolation-window search multiplies candidates/spectrum; the protein-walk +
bucket scan slows ~1.2–2.7×. **Known solution: fragment-ion index candidate
generator** (the project's planned "speed-v2"; the exact mechanism MSFragger uses to
do wide-window search fast). The PSMs are already there — they just need to be found
faster. This is the highest-leverage next build for the gate.

### Blocker 2 — TMT PSMs (−5%)
Chimeric barely moves TMT (9,686 vs narrow Rust 9,605; Java 10,194). TMT co-isolation
isn't yielding gains. TMT remains its own problem — the GF SpecEValue-shape divergence
(Lever-2a), orthogonal to chimeric.

## Implications / path

1. **Chimeric is reinstated as a real PSM lever** (PXD+Astral), trust-proven. The
   four prior "refutations" were broken-ruler (count-vs-narrow) artifacts.
2. **The gate now hinges on SPEED** for PXD+Astral (fragment index) and on **TMT**
   (Lever-2a / GF-shape) for PSMs. Neither is a chimeric-FDR problem.
3. The rescore (Phase 3) is dead weight — NO_RESCORE is faster and equally
   trustworthy; keep `--chimeric` on the NO_RESCORE path.
4. **Fairness caveat (still open):** Rust searches a wider window than narrow Java.
   The honest "beats Java" claim should be cross-checked by running the SAME
   entrapment-FDP ruler on Java, and at peptide level. Pending.

## Recommendation

The PSM win is real; the binding constraint is **speed**. Pursue the **fragment-index
candidate generator** (speed-v2) as the enabler that could let chimeric clear the
gate on PXD+Astral, and handle TMT separately via Lever-2a. Until speed lands,
chimeric does not clear the gate; nothing ships ([[merge-gate-beat-java]]).

## Artifacts
- VM: `/srv/data/msgf-bench/gate-chimeric-norescore/`
- Java baseline: `/srv/data/msgf-bench/bench-java-3ds-percolator/`
