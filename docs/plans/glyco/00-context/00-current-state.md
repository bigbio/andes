# andes glyco — current state, failure diagnosis, what we have vs pending

*Authoritative status as of 2026-07-02 (branch `glyco-phase1`, HEAD ~`35d31bb9`).*

## 1. The bottleneck chain (measured, cross-validated)

Truth set: PXD025455, file `HCC_pool_Late_Fc3_r1`, **523** the reference engine-nglycan truth
scans (our own re-search; the dataset was originally a commercial glyco engine).

```
~90% backbone-findable (generation near-ceiling; PHASE1_RESULT.md)
   → 154/523 recovered @1% FDR      (29.4% — the "stuck <30%")
   → 83/154 top-1 correct           (ranking loses 71)   ← SP-B / G2  (~70% headroom)
   → 66 true @1% FDR                (111 false pass)      ← 2D-FDR / G3
```

Two independent analyses (an external AI review the user relayed, plus our own
`SPA2_RESULT.md`) and our measurements all converge:
**generation is solved; identification (ranking + FDR) is the failure.**

## 2. What we measured this session

- **G1 (glycan-Y-first candidate selection) — VERIFIED.** Fast harness
  `ANDES_GLYCO_SCANS=<file>` (commit `90159852`) subsets the driver to the 523
  truth scans (~8 min/arm vs ~100 min). Clean A/B, one variable
  (`ANDES_GLYCO_YINDEX`): backbone-findability **59.3% → 69.8% @0.05 Da**, robust
  to tolerance (**+7.2 pts even at 0.005 Da ≈ exact**) → genuine coverage, not
  chance. Generation-side, modest (generation was already ~90%).
- **G3 (glycan-axis decoy) — SCORER SOUND, UNIFIED SCHEME REFUTED.** The decoy
  Y-ladder scorer `glycan_y_intensity_decoy` (commit `b3b298cd`, TDD: target >
  decoy) is correct. But feeding the ~352K glycan-decoy rows into Percolator under
  **one `Label`** *crashed* recovery **29.4% → 4.4%** (they differ from targets
  only in YLadder → dominate the −1 pile → Percolator over-weights YLadder and
  under-ranks genuine targets). **Lesson: 2D-FDR must be the a glyco search engine separate-axis
  post-process, not a unified Percolator pile.**

## 3. What we have (implemented on `glyco-phase1`)

- Oxonium gate; de-novo Y-ladder backbone solver + complement scoring.
- Clean-room glycan DB (≈600 common / 2510 full compositions).
- Y-first cascade + DB union/fallback (`hybrid.rs`).
- Peptide-first fragment index (`ANDES_GLYCO_PEPTIDE_FIRST=1`, default on).
- Composition-specific `YLadderScore` (was a dead 0.0 feature).
- Glycan-Y index + two-axis retention (`ANDES_GLYCO_YINDEX`, opt-in) — **G1**.
- Glycan-axis decoy ladder + paired PIN rows (`ANDES_GLYCO_DECOY`, opt-in) — **G3**.
- Cross-spectrum transfer scaffold (`ANDES_GLYCO_CROSSSPECTRUM`, opt-in, NULL v1 —
  needs RT gating).
- Fast dev harness (`ANDES_GLYCO_SCANS`); many correctness fixes (isotope, H2O
  convention, CalcMass, determinism, decoy-gating perf).

## 4. What is pending (the real work)

| Tag | Item | Why |
|---|---|---|
| **P1.0 / SP-B** | Learned peptide-axis glyco scoring (regime-matched strong model + Y0/Y1 anchor feature) | fixes ranking (loses 71/154) — **#1 lever** |
| **G4** | RT-gated cross-spectrum transfer (a cross-spectrum glyco engine) | recovers sparse-b/y stratum single-spectrum scoring can't rank |
| **G3′** | Separate-axis 2D-FDR post-process (not unified pile) | valid glyco FDR; fixes true-FDP |
| G0 | DET-1 (`hybrid.rs:318` → `total_cmp`), P0.3 (Y0/Y1-only quorum, *measure first*), P0.4 (probe isotope fidelity) | correctness hygiene |

## 5. The SP-B / G2 seam (from the reuse audit)

- The glyco driver scores backbone b/y via `score_psm(ss, …, scorer, …)`
  (`glyco_search.rs:577`) using the **single** `prepared.scorer` selected by
  `SelectionKey{activation, instrument, enzyme, experiment_class}`. **No
  glyco-specific scorer path** — auto-selection is purely a matter of the
  `SelectionKey` built at `build_selection_key` (`andes.rs:4809`).
- The model store is **Hive-partitioned by `protocol=<P>`**
  (`resources/models/protocol=TMT`, `=Phosphorylation`, …). `andes train` ingests
  externally-labeled PSMs from a **flat parquet** (every row a label, no bootstrap
  search); `split_store_by_protocol` writes a new partition generically.
- **~90% reusable.** Load-bearing new work:
  1. `Protocol::NGlyco` variant (model enum + CLI) + a `build_selection_key` arm →
     `experiment_class={"glyco"}` (the `glyco` catalog slug already exists).
  2. **Glycan-stripped backbone training rows** — the glycan mass must NOT shift
     backbone b/y; rows carry the bare peptide. *This is the real modeling decision.*
  3. Truth-TSV → flat-parquet converter (a script).

## 6. Why retraining alone is likely insufficient (the honest premise)

Stepped-HCD glyco b/y is *physically* sparse. a cross-spectrum glyco engine (Nat Commun 2022):
direct b/y IDs the backbone in only **~11%** of spectra; its +48% PSM gain came
from **cross-spectrum transfer**, not scoring. A better-calibrated model reshapes a
likelihood over the same sparse peaks — it cannot manufacture b/y ions the
fragmentation never deposited. So:
- **Kill-gate on decoy-separated ranking** (true peptide's b/y-score vs a
  same-backbone decoy's), not find-rate.
- **Pair retraining with the Y0/Y1 peptide-mass-anchor feature** — a
  peptide-mass-conditioned feature *does* discriminate competing peptides (unlike
  backbone-level oxonium/YLadder), and the anchor is high-intensity even when b/y is
  dead.
- **If the gate fails (likely), the fix is cross-spectrum transfer (G4).**

See `50-roadmap/` for the phased plan and `00-context/01-spb-brainstorm-synthesis.md`
for the full design dialogue.
