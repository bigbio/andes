# andes Second-Pass Redesign — Unified-FDR Anchored Reinterpretation (+ confident-subset fragment index)

**Date:** 2026-06-18 · **Branch:** feat/ptm-refinement-cascade (supersedes the cascade in `2026-06-17-andes-ptm-refinement-cascade-design.md`)

**Status:** design, approved-direction. Replaces the disjoint-union cascade with a per-scan unified-FDR competition (Phase A) and a confident-subset fragment-ion open search (Phase B).

---

## 1. Why (the evidence that forced the redesign)

The shipped cascade (anchored Pass-2, reversed-base-peptide decoys, separate "refine" PIN, disjoint-union FDR) was validated to *recover* ~45,411 modified PSMs@nominal-1% on PXD001468 — but the entrapment-FDP cross-check (2026-06-18, VM) showed the nominal 1% is really **5.75% true combined entrapment-FDP** (≈5.7× optimistic; honest 1% = ~31,500 PSMs). Four controlled probes on `b1931.raw` (1:1 ENT_ entrapment DB, RichIonLLR-active model) isolated *why*:

| Probe | Result | Lesson for the design |
|---|---|---|
| **EXP-1** — can a score gate separate real-mod from FORWARD entrapment? | Best ROC-AUC **0.88** (RichIonLLR + NumMatchedMainIons); RankScore useless (0.53), and adding it *dilutes* | A score gate **trims but cannot fix** the leak — some forward false positives are genuinely score-overlapping. RankScore is the wrong selection score for the modified regime. |
| **EXP-2** — precursor-mass open search headroom | **Coincidence** — 6 random non-PTM offsets match **80%** of unidentified scans; the +15.995 ox bin does not rise above noise | Precursor-mass-keyed open search is **noise**. Real open search MUST be **MS2 fragment-ion constrained**, not a precursor window. |
| **EXP-3** — double-assignment channel | **68.3%** of refine-accepted scans are *also* Pass-1-identified (unmodified) | Two-thirds of "modified" IDs are redundant reinterpretations admitted on a *separate* FDR axis — pure artifact of the disjoint union. |
| **EXP-4** — decoy hardness (reverse vs shuffle) | Shuffle **6.04%** vs reverse 5.75% true-FDP — *worse* | The 5.7× optimism is **not** a decoy-quality problem. The false positives are *forward* (real-looking tryptic) peptides no decoy null models; the construction itself is optimistic. |

**Root cause:** the second pass runs as a *separate search on a separate FDR axis* ("Pass-1 PIN ⊎ refine PIN"). That single architectural choice produces both the 68.3% double-counting (EXP-3) and the 5.7× FDR optimism (EXP-4) — because the "unidentified" gate is the weak internal `rank_score` TDC while the report is Percolator, so a scan can win in *both* lists, and the modified subgroup is never forced to compete against the unmodified explanation in one calibrated null.

**Design principle (what the best engines do):** modified forms compete **per-scan, in one unified target-decoy list** (MSFragger / MaxQuant main search). andes keeps its memory-feasible **anchored** scoping (never expand the whole proteome — global semi/PTM OOMs at 31 GB) but drops the disjoint union.

## 2. Objectives & constraints (the standing andes targets)

- **More PSMs than the majority of engines** at *honest* FDR — the deliverable is net-new + correctly-upgraded scans at **true** ≤1% entrapment-FDP, not nominal-1% TDC.
- **Similar speed** — the current refinement phase is ~9 s; the redesign must stay in the seconds regime (Phase A ≲2× refine cost; Phase B's index is built over the *confident subset only* to stay fast).
- **Own data / no patent** — reuse the own-trained RankScorer + RichIonLLR; no MS2PIP/NN; public datasets only (PXD001468 is the canonical PTM-rich testbed, already staged).
- **Additive-PIN discipline** — never modify existing score columns (proven to regress Percolator); new signals enter as additive columns. The *FDR construction* changes; the per-candidate scoring math does not.

## 3. Architecture

### Phase A — Unified-FDR anchored reinterpretation (cheap; the correctness fix; ships the honest number)

Replace the disjoint union with one per-scan competition + one FDR:

1. **Pass-1** unchanged: full-DB search, all spectra → per-scan `TopNQueue` of unmodified candidates; confident base peptides via the existing TDC walk.
2. **Anchored modified candidates over ALL scans (not just "unidentified").** Build the anchored modified candidate pool ({confident base peptides × refine-tier mods} + paired decoys) and search it precursor-gated against **every** scan. Because candidate generation is precursor-windowed, a scan only scores the modified candidates whose modified precursor mass matches it — so "all scans" costs ≈1.7× the current unidentified-only refine, still seconds. This **eliminates the internal-TDC-vs-Percolator gap** (no reliance on a separate "unidentified" definition): the modified hypothesis competes head-to-head with the unmodified one on every scan.
3. **Unified candidate index space (the core engineering obstacle).** The shipped code keeps two PINs because a merged `PsmMatch::candidate_idxs` would alias across the Pass-1 vs Pass-2 candidate slices (different `protein_index` spaces → wrong accession/OOB). Phase A resolves this: the Pass-2 anchored candidates are appended to a **single combined candidate list** with a consistent index space (or carry an explicit `source + offset` the PIN/TSV writer resolves through), so unmodified and modified PSMs for the same scan live in **one queue** and **one PIN**.
4. **Per-scan best-PSM competition.** Each scan's emitted PSMs include its unmodified winner *and* its modified winner(s); the single downstream Percolator/mokapot run picks the best PSM per scan. A modified PSM is reported **only if it beats its own scan's unmodified interpretation** → the 68.3% redundant channel collapses to legitimate *upgrades* (marginal unmodified ID → confident modified ID = a real correctness/PSM gain) plus genuinely net-new modified scans.
5. **One calibrated FDR.** A single target-decoy list; group-FDR by mod-class via mokapot `--group-column` (subgroup validity, Bogdanow/Selbach); modified and unmodified rows share the decoy null so the q-values are mutually calibrated. `is_refinement`, `num_mods`, `refine_mod_class` remain additive PIN columns; RichIonLLR + NumMatchedMainIons (the AUC-0.88 separators) are the discriminative features the unified rescorer leans on for the modified group.

### Phase B — Confident-subset fragment-ion index, open/delta-mass search (the novel, high-PSM ceiling)

EXP-2 says coverage expansion beyond the 5-mod tier requires an MS2 index. The novel andes move: build the index over **only the confident-peptide subset**, not the whole proteome.

1. **Fragment-ion index** over the confident base-peptide backbones (theoretical b/y m/z bucketed). Small (thousands of peptides) → fast to build, memory-trivial — the andes-specific innovation that keeps speed where MSFragger/MaxQuant index the whole DB.
2. **Open per-scan match:** observed peaks → index → candidate base peptides sharing many **unshifted** fragments; precursor delta = scan precursor − base mass over an open range (e.g. [−100, +300] Da) → **localize** the delta by which shifted-ion series explains the residual peaks (MSFragger localization-aware: shifted + regular ions, Yu et al. 2020). Discovers **any** modification, not just the tier.
3. **Delta-aware RichIonLLR:** extend the per-ion LLR to the shifted ladder so the score *selects* on the modified evidence (closing EXP-1's "RankScore is the wrong selection score" gap).
4. **Two-axis FDR into the same unified list:** reversed/shuffled anchored decoys (peptide axis) + **delta-decoys** (random non-chemistry mass shifts) for the mass-shift FDR axis. Feeds Phase A's per-scan competition + group-FDR.
5. **Phase-B headroom gate (run before heavy build):** an empirical probe — does the confident-subset fragment index recover documented PXD001468 mods at honest entrapment-FDP that the 5-mod tier misses? Only commit the full build if the MS2-constrained headroom is real (EXP-2 showed precursor-mass headroom is not).

## 4. Data flow

```
spectra ─▶ Pass-1 search (full DB) ─▶ per-scan unmod queues ─┐
                                      confident base peptides ─┤
                                                              ▼
        anchored modified candidates (tier; Phase A) ───▶ precursor-gated search over ALL scans
        confident-subset fragment index (open; Phase B) ─┘            │
                                                                       ▼
                              MERGE per-scan: {unmod ∪ mod} in ONE candidate index space
                                                                       ▼
                              ONE PIN (additive is_refinement/num_mods/refine_mod_class + RichIonLLR…)
                                                                       ▼
              Percolator/mokapot: single run, --group-column, best-PSM-per-scan ─▶ report @ entrapment-calibrated q
```

## 5. FDR design (non-negotiable)

- **Single unified target-decoy list.** No disjoint union; no separate refine FDR axis.
- **Per-scan best-PSM competition** before FDR → modified must beat unmodified on its own scan → kills EXP-3's 68.3%.
- **Group-FDR by mod-class** (mokapot `--group-column`) for subgroup validity; sparse classes (<k≈20 decoys) fold into "other".
- **Symmetric expansion** of targets and decoys (same mods, same caps).
- **Phase B adds a delta-decoy axis** so the mass-shift assignment has its own null.
- **Mandatory entrapment-FDP gate:** every milestone is measured with the existing peptide-level entrapment harness (`build_entrapment_db.py` + `pep_entrap_fdp.py`/`pep_entrap_curve.py`); the merged true combined-FDP must be ≤1%. This is the ship-gate — nominal Percolator q is not trusted (EXP-4).

## 6. Code structure (deltas to the current tree)

- `crates/search/src/refinement.rs`: keep `confident_base_peptides`/`refinement_aa_set`/`mod_count_and_class`; **replace** `run_refinement`'s disjoint-union output with a merge that returns modified candidates in a **combined index space** + per-scan tags; refine over all scans (precursor-gated) not only `unidentified_spectrum_indices` (the latter is retired as the FDR-defining gate, optionally kept as a cheap prefilter).
- `crates/search/src/match_engine.rs`: candidate-index unification so Pass-2 PSMs merge into the Pass-1 queues without aliasing; per-scan competition retains both interpretations.
- `crates/output/src/pin.rs` + `crates/search/src/psm.rs`: ONE PIN with the additive group columns (already present); ensure modified + unmodified rows for a scan share `SpecId` scan so Percolator does best-per-scan.
- `crates/andes/src/bin/andes.rs`: `--refine` now emits one PIN; drop the `.refine.pin` second file (or keep behind a debug flag). Report wording updated (no "disjoint union").
- **Phase B (new):** `crates/search/src/fragment_index.rs` (confident-subset b/y index), delta-aware extension in `crates/scoring` (RichIonLLR over shifted ions), delta-decoy generation. Gated behind the Phase-B headroom experiment.

## 7. Validation gates (each milestone)

1. **Entrapment-FDP** on PXD001468 (PXD001468 b1931 fraction): merged true combined-FDP ≤1% for the unmodified AND modified groups, and per-mod-class. Baseline to beat: disjoint-union 5.75%@nominal-1%.
2. **Net-new + upgraded scans** at honest 1% vs the disjoint-union baseline — must not lose honest PSMs; target a gain from upgrades + removed double-counts being reallocated.
3. **Speed:** refine phase wall-clock within ~2× of current (seconds); Phase B index build amortized.
4. **Determinism A/B:** stripping the new/additive columns reproduces the Pass-1 baseline byte-identically (the proven safety check).
5. **Cross-engine:** PSMs@honest-1% vs MSFragger/MaxQuant open search on the same fraction (the "beat the field" check).

## 8. Phasing

- **Phase A (this spec's MVP):** unified per-scan competition + single calibrated FDR over the existing bounded tier. Data-mandated correctness fix; no new search infra. The reported modified count will be *lower than* the disjoint-union's honest ~31.5k (the 68.3% double-counts are removed) but *honest and cleaner* — net contribution = genuinely net-new modified scans + correctness upgrades of marginal unmodified IDs, all at true ≤1% entrapment-FDP. The honest baseline to beat is whatever the disjoint-union nets after deduplication, not its inflated 45,411.
- **Phase B (same infra, follow-on):** confident-subset fragment index → open/delta-mass + delta-aware RichIonLLR + delta-decoy FDR. The coverage ceiling ("more mods than the field"), gated by its own headroom experiment.

## 9. Non-goals (deferred)

- Whole-proteome open search / fragment index (the thing that OOMs / that EXP-2 proves needs MS2 constraint anyway) — andes stays confident-subset-anchored.
- In-engine FDR — andes remains a PIN emitter; group-FDR via mokapot downstream.
- Precursor-mass-only open search — empirically refuted (EXP-2, 80% coincidence).

## 10. References

- Disjoint-union cascade it replaces: `internal-docs/specs/2026-06-17-andes-ptm-refinement-cascade-design.md`.
- Empirical basis: this session's EXP-1..4 (PXD001468 `b1931.raw`, 1:1 ENT_ entrapment, `pep_entrap_*.py`).
- Yu et al., "Identification of modified peptides using localization-aware open search," *Nat Commun* 2020, 11:4065 (shifted + regular fragment ions; fast mass calibration) — Phase B blueprint.
- Bogdanow/Zauber/Selbach, *MCP* 2016 (subgroup FDR); Fondrie & Noble, *JPR* 2021 (mokapot grouped confidence) — group-FDR.
- Kertesz-Farkas/Keich/Noble, *JPR* 2015 (cascade-search per-tier decoys); Tyanova/Temu/Cox, *Nat Protoc* 2016 (MaxQuant dependent-peptide / main-search per-scan competition).
