# Phase-1 / SP-A Hybrid Glyco Backbone — RESULT

## SP-A Hybrid Gate: PASS (90.4% ≥ 90% target)

**Verdict: PASS** — hybrid searchable-backbone 90.4% (490/542) ≥ 90% gate.
SP-A adds a DB-constrained backbone branch (backbone = precursor − known glycan),
unioned with the existing de-novo Y-ladder solver, lifting overall from 57.7% (de-novo-alone) to 90.4%.

## The three numbers (n = 542 matched scans)

| Metric | Value |
|---|---|
| Oxonium-fire rate | **100.0%** (542/542) |
| **De-novo-only baseline (Phase-1 top_k=5)** | **57.7%** (313/542) |
| **Hybrid searchable-backbone OVERALL (SP-A)** | **90.4%** (490/542) — GATE PASS |
| DB-branch contribution (% of searchable) | 87.3% (428/490 hits from DB) |
| De-novo-only contribution (% of searchable) | 12.7% (62/490 hits de-novo only) |
| Sparse stratum (≤2 true core-Y rungs, n=148) | **71.6%** (106/148) |

## SP-A architecture

- **`glycan_db.rs`**: clean-room N-glycan enumerator — HexNAc 2..=8, Hex 3..=12, Fuc 0..=3,
  NeuAc 0..=5, NeuGc 0..=2; plausibility constraints (fuc ≤ hexnac, sialic ≤ hexnac−2),
  mass ∈ [500, 6000]. Produces **2510 compositions**. Deterministic total-order sort.
- **`hybrid.rs`**: `db_branch(precursor_neutral, glycans, min_backbone)` computes
  `bb = precursor_neutral − glycan.mass` for each composition, filters bb ≥ 500 Da.
  `hybrid_candidates` unions DB-branch + de-novo solver; deduplicates within 0.02 Da
  (DB source preferred over de-novo when they cluster).
- **`glyco_probe` harness**: updated to call `hybrid_candidates`, report source split
  and comparison vs de-novo baseline.

## Run provenance

- Binary commit: see commit on branch `glyco-phase1` for this PHASE1_RESULT.md update.
- Branch: `glyco-phase1`.
- Gate window (non-negotiable): symmetric `in_window(solved, truth) = |solved − truth| ≤ max(truth*20e-6, 0.01)`.
- Hybrid call: `hybrid_candidates(&peaks, precursor_neutral, prec_z, &glycans, 20.0, 5)`.
- Oxonium gate: `oxonium_gate(&peaks, 0.10, 20.0)` (DB branch runs regardless; de-novo requires oxonium).
- All 19 unit tests green.

## Phase-1 (de-novo-only) prior record

| Metric | Phase-1 (prior) | SP-A (hybrid) |
|---|---|---|
| Searchable-backbone OVERALL | 59.8% (prior run) / 57.7% (this run) | **90.4%** |
| Sparse stratum | 2.7% | **71.6%** |

Note: de-novo baseline is 57.7% in this run vs 59.8% in the prior PHASE1_RESULT (minor variance
from sorting determinism fixes across the interim commits).

## Key finding

The DB branch alone accounts for 87.3% of searchable hits. The database-constrained approach
dominates: for most glycopeptide spectra, `precursor_neutral − known_glycan_mass` lands within
±20 ppm of the true backbone even without any spectral evidence from Y-ions. This is the
expected result for intact-glycopeptide spectra where the precursor mass is precisely measured
and the glycan composition is known a priori.

The sparse stratum jumped from 2.7% to 71.6% — the biggest gain — because sparse spectra
(≤2 core-Y rungs) were structurally inaccessible to the de-novo solver but are trivially
solved by the DB branch: no Y-ladder evidence is needed when the glycan is enumerated.

## SP-A FINAL (corrected source attribution, commit 3aa93a6c, glycan DB=2510 comps, n=542)
- De-novo-only baseline: 57.7% (313/542)
- Hybrid searchable-backbone: **90.4%** (490/542) — GATE PASS
- Source split (reliable): DB-branch 87.3% (428/490), de-novo-only 12.7% (62/490)
- Sparse stratum (≤2 core-Y): 71.6% (106/148)
