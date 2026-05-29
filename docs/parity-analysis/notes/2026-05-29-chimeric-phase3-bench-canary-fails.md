# Chimeric Phase 3 bench — residual SpecEValue re-score does NOT deflate the Astral canary

**Date:** 2026-05-29
**Branch:** `feat/chimeric-dda-plus` (impl commit `3456d386`)
**Binary:** chimeric-build rebuilt with Phase 3 (greedy shared-fragment competition
+ residual SpecEValue re-score + additive unique-evidence PIN columns), cal=auto.

## Result (Percolator @1% FDR, `run_bench_chimeric.sh`, off vs on)

| Dataset | off @1% | on @1% | Δ | decoy fraction off→on |
|---|---:|---:|---:|---|
| PXD001819 | 14,808 | 18,306 | +3,498 (+23.6%) | 0.40 → 0.83 |
| **Astral (canary)** | **36,715** | **77,444** | **+40,729 (+111%)** | 0.52 → 0.83 |
| TMT | 9,605 | 9,362 | −243 (−2.5%) | 0.56 → 0.90 |

Java references: PXD 14,974, Astral ~36,271, TMT ~10,115.

PIN rows: PXD 39,171→365,336; Astral 137,641→**1,229,207**; TMT 46,332→373,241.
Wall (on): PXD 7:03, Astral **40:09**, TMT 12:58 (residual re-score rebuilds a
ScoredSpectrum + GF group per rank≥2 PSM — ~5× the off wall on Astral).

## Verdict: the canary FIRES — Phase 3 does NOT make chimeric trustworthy

Astral has narrow isolation windows → minimal real co-isolation → a *correct*
method must leave the count ~flat at ~36.7k. Instead on=**77,444 (+111%)**, even
slightly *above* the pre-Phase-3 chimeric (~72k, Phase-2 bench). The residual
SpecEValue re-score **did not deflate the inflation.** Decoy fraction stays ~0.83
(near 1:1) on all three datasets — the structural-inflation signature is intact.

PXD's "18,306 > Java 14,974" is **not a real win** — the Astral canary proves the
method inflates (same conclusion as Phase 2). TMT is net-negative (−2.5%),
consistent with chimeric being a no-op/loss on narrow windows.

## Root cause: a per-PSM score change cannot fix multiple-testing FDR inflation

The residual re-score deflates a spurious **target** that stole/coincidentally
matched peaks — but it deflates the spurious **decoy** at the same mass
*symmetrically*. PSM-level target-decoy FDR (the q-value curve) is driven by the
*relative* target-vs-decoy score distribution; shifting both sides down together
leaves the q-value threshold essentially where it was. So even though individual
theft rows demonstrably deflate (BSA smoke: 37 fully-stolen rows → negative
residual RawScore, lnSpecE≈0), the *aggregate* 1% FDR count does not move.

This is the empirical confirmation of the Phase-2 post-mortem's requirement #2:
the broken part is the **FDR model** (PSM-level TDC over ~5 PSMs/scan), not the
per-peptide score. Fragment competition is real (overlap diagnostic confirmed
theft), and the residual re-score is mechanically correct, but neither is
*sufficient* — you cannot restore a credible FDR by rescoring within an inflated
multi-PSM-per-scan set. The two inflation drivers (Phase-2 addendum) — coincidental
multiple-testing (~28% of overlap scans keep their own peaks) and the structural
T/D edge of wide-window multi-emission — are both untouched by a symmetric
per-PSM rescore.

## Implications

- **Phase 3 (residual SpecEValue re-score) is REFUTED as a chimeric gate-clearer.**
  Do not pursue the wall optimization (node-score-only fallback) — the approach
  doesn't deflate the canary even at full GF fidelity, so making it faster is moot.
- **Chimeric still does NOT pass [[merge-gate-beat-java]]** — untrustworthy on the
  Astral control; net-negative on TMT.
- The implementation is **kept on the branch** (committed `3456d386`) as a
  validated negative result + reusable machinery (greedy competition, residual
  re-score, unique-evidence columns). `--chimeric off` is byte-identical; nothing
  ships to `dev`.

## What would actually be required (unchanged from Phase 2, now doubly confirmed)

A **per-scan / peptide-level FDR** (or rank-1-vs-rest competition that models the
multi-PSM-per-scan structure, à la MSFragger-DDA+ + Philosopher) — a structural
change to how FDR is computed, NOT another per-PSM score or feature. Fragment
competition would then operate *within* that credible model. Until that lands,
chimeric stays shelved.

## References

- `2026-05-29-chimeric-fragment-overlap-diagnostic.md` — theft confirmed (bimodal).
- `2026-05-28-chimeric-phase2-bench.md` — soft feature insufficient; the 3
  requirements (the per-scan-FDR requirement is the one this bench re-confirms).
- `superpowers/specs/2026-05-29-chimeric-phase3-shared-fragment-design.md` — design.
