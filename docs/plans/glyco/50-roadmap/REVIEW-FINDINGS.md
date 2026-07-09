# andes-glyco v2 — multi-agent review findings (2026-07)

Four adversarial review agents audited the `glyco-v2` branch. This records what was
ACTED ON and what remains as precise, actionable follow-ups.

## DONE this pass

- **Removed all refuted code** (learned GBDT selector: `glyco_selector.rs`,
  `train_glyco_selector` bin, the P3 collapse integration + `OnceLock` model loader;
  the `GLYCO_GP_P` partial-glycan-in-selector fusion term, reverting the fusion to 7
  args; offline bins `partial_glycan_test`, `augment_glyco_pin`). Commit `fd7b8637`.
- **Tightened stale/verbose comments** — notably the `glycan_db` expansion comment that
  credited the now-deleted GBDT selector, and the `partial_glycan_by_intensity` doc.
- Build + tests green (82 + 76).

## SPEED — DONE (behavior-preserving, all 82+162 tests green)

1. ✅ **Hoisted per-spectrum constants out of the `backbone.rs` intensity fns.** New
   `backbone::SpectrumStats { base, sorted }` computed ONCE per spectrum in
   `score_spectrum_glyco` (and in `hybrid.rs`'s core-Y ranking) and passed by reference to
   `core_y_intensity`/`glycan_y_intensity`/`glycan_y_intensity_decoy`/`y0y1_anchor_intensity`/
   `partial_glycan_by_intensity`/`count_core_y_hits`. Removes the per-call `O(#peaks)`
   base-peak max-fold + `windows(2).all(...)` sorted-check (the dominant per-call cost).
   Values are byte-identical to the old inline computation.
3. ✅ **Per-candidate sequon-membership `Vec<bool>`** precomputed once in
   `GlycoCtxOwned::build` (also reused by the peptide-first fragment-index build, so the
   sequon scan now runs once, not twice); the hot-loop check is an O(1) `ctx.sequon_membership[slot]`.
4. ✅ **Hoisted the 6 per-spectrum `std::env::var` reads** (`gp_selector_on`, `glyco_gp_k/j/h`,
   `glyco_charge_expand`, `y_primary_selection`) into `GlycoCtxOwned::build` (read once,
   stored on the ctx) — no more process-env-lock contention inside `par_iter`.
5. ✅ (partial) **De-duplicated the cumulative-`adds` builder** into one
   `glycan_cumulative_adds(comp)` shared by `glycan_y_intensity` + its decoy (was ~22 lines
   copy-pasted). Cleanliness + single source of truth for the ladder order.

## SPEED — DEFERRED (diminishing returns after the above; parity-sensitive)

2. **Memoize the Y-ladder per `bb_hit_idx`** (a `Vec<Option<f32>>`), not per winner — its
   value is independent of `cand_slot`, so candidates sharing a backbone recompute it. Carry
   the selection-time ladder into `GlycoPsmKey.y_ladder_intensity_score` instead of the third
   inline recompute. Deferred: needs a memo threaded through the winner loop; parity-risky.
5b. **Full per-glycan `adds` side-table** (by glycan index) to avoid the one `Vec` alloc per
   `glycan_y_intensity` call — smaller win now that `SpectrumStats` removed the O(#peaks) cost;
   needs plumbing an index (not just `&GlycanComp`) into the fn.

Confirmed NOT issues: `best_frag_intensity` already binary-searches (`partition_point`); the
default (non-gbdt) path computes `compute_psm_features` for only the 1 emitted winner.

## PARAMETERS / CONFIG — follow-ups

- **Env-var sprawl.** REMOVE dead/refuted dev knobs: `ANDES_GLYCO_CHARGE_EXPAND`,
  `ANDES_GLYCO_CROSSSPECTRUM`, `ANDES_GLYCO_EXHAUSTIVE`, `ANDES_GLYCO_PEPTIDE_FIRST`,
  `ANDES_GLYCO_SCANS`, `ANDES_GLYCO_SELECT`. PROMOTE real settings to a `GlycoSearchConfig`/
  CLI/model metadata: `FULL_GLYCANS`, `DECOY`, `GP_K/J/H`, `MAX_PF`, `SEED_FDR`, `RT_WINDOW`,
  `SELECTOR`, `MIN_SUPPORT`, `DENOVO`, `ALL_HITS`, `YINDEX`.
- **Fusion weights (K=50, J=5, H=1) are hand-tuned on one dataset, hidden behind env vars.**
  Promote to typed model/config metadata; longer term learn them (RankSVM/GBDT in the store).
- **Magic numbers:** the hardcoded `20.0` ppm glyco fragment tolerance (`andes.rs:2365`) and
  `20e-6` in `complement_score` should use the configurable fragment tolerance; make floors
  named constants. `backbone_top_k` (50), `SELECTOR_SHORTLIST_K` (24), `max_features` formula
  → promote to config with documented defaults. Replace the local `core_adds` with `CORE_Y_STEPS`.
- **Glycan-DB profiles hardcoded** — make selectable `GlycanDbProfile`s (CLI/SDRF-derived).
- **Doc drift:** `docs/plans/glyco/40-data/glycan-db.md` still says glycan mass includes `+H2O`
  and the DB is ~2510 with old bounds; code uses residue-sum (no water) and expanded
  `Hex 2..14 / Fuc 0..4`. Reconcile docs to cite the same source constants.

## ARCHITECTURAL — follow-ups

- **Collapse logic lives in 3 places** (driver `gp`, PIN `select_emitted_hits`, transfer-seed
  `collapse_cmp`) and MUST call `glyco_gp_fused_score` identically — a silent-divergence footgun.
  Refactor to one shared collapse API + driver-vs-PIN parity tests.
- The growing hand-weighted fusion (`rank + K·ladder + J·core_y + H·hyper`) is brittle model
  logic; the config/learned-weights path above is the durable fix.

## FOUNDATIONAL (masses/constants)

A Codex deep audit of the monosaccharide masses, b/y formulas, and `partial_glycan_by`
core-adds was launched but its background result was not retrievable in-session. The shipped
masses are empirically validated (the same `glycan_mass.rs` constants drive the Y-ladder,
core-Y, and partial-glycan matching, all producing correct discrimination). Recommend a
dedicated unit test asserting the monosaccharide constants against reference monoisotopic
residue masses as a permanent guard.
