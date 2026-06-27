# Own-geometry bundle — quality A/B vs the prior shipped bundle (2026-06-26)

The own-geometry bundle (commit `d4357def`) re-derives every model's partition
geometry from andes's own corpus instead of copying the MS-GF+ seed template —
closing the last structural carryover. Because that is **result-changing**, it
was validated against the prior shipped (seed-geometry) own-only bundle before
adoption.

## Method
Per-model **pinned** A/B (isolates the geometry+retrain change; routing
correctness is covered by unit tests). Same andes binary, same target-only
entrapment FASTA per dataset, `--candidate-index ram`, default scoring. Metric =
PSMs at **1% true entrapment-FDP** (1:1 entrapment → FDP = 2·ENT/total; deepest
score-sorted prefix with FDP ≤ 1%) — mode-independent, the canonical comparison.

| dataset | model | OLD (seed-geometry) | NEW (own-geometry) | Δ |
|---|---|---|---|---|
| Astral (high-res) | `hcd_astral_tryp` | 40,643 (1.00%) | 36,730 (1.00%) | **−9.6%** |
| UPS1 (low-res LFQ) | `cid_lowres_tryp` | 15,823 (1.00%) | 14,919 (0.99%) | **−5.7%** |
| TMT (low-res CID) | `cid_lowres_tryp_tmt` | 10,582 (0.98%) | 11,215 (1.00%) | **+6.0%** |

## Reading
- Own-geometry **regresses Astral and UPS1, improves TMT**, at matched FDP — a
  real quality change, not an FDP artifact.
- vs the field: NEW Astral 36,730 still leads Comet/MSFragger (~28–29k); NEW UPS1
  14,919 falls **below Java MS-GF+ (~15.9k)** — the prior bundle was at parity.
- **Confound:** OLD and NEW differ in geometry **and** training corpus. The
  own-geometry `hcd_astral_tryp` trained on 218,738 PSMs vs the prior wave's
  larger (~481k) corpus, so part of the Astral loss may be corpus size, not
  geometry. Not cleanly attributable without a same-corpus re-run.

## Decision (user, 2026-06-26)
**Ship own-geometry as-is** — accept the regression in exchange for full geometry
independence (no MS-GF+ code, constants, or geometry) and the Apache flip. Open
item before the *public* benchmark claims are republished: re-run the headline
multi-engine q≤1% table with this bundle (the prior table — incl. "ties Java on
UPS1" — was measured on the seed-geometry bundle).
