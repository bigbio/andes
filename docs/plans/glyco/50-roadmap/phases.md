# andes N-glyco — per-phase implementation specs (G0 → G4 + G3′)

> ⚠️ **CORRECTION (2026-07-03) — read [`../LESSONS.md`](../LESSONS.md) first.** Any
> gate below citing the **523-scan truth** or **`154/523 → 83 → 66`** is **VOID**
> (multi-row-PIN artifact, fixed in commit 7c269aab; correct answer is 0 glyco-PSMs
> @1% FDR). Re-express all gates as **unique glyco scans/glycopeptides @ FDR on
> top-1-collapsed PINs, vs a valid FDR-controlled reference, entrapment-validated.**


*Engineer-executable checklist. Companion to `50-roadmap/roadmap.md` (sequencing +
thesis). Branch `glyco-phase1`, HEAD ~`35d31bb9`. File:line anchors from
`00-context/02-code-inventory.md`; theory from `20-theory/`; masses/notation from
`30-standards/`; data policy from `40-data/` + roadmap §5.*

**Read the roadmap first.** This file is the *how*; the roadmap is the *why/when*.
Phases map onto the prompt's labels as: **G0** = correctness hygiene · **P0** =
Y0/Y1 anchor + **SP-B kill-gate** (roadmap "G2") · **P1** = harvest + regime-matched
retrain / SP-B model (roadmap "G3") · **G3′** = 2D-FDR post-process · **G4** =
cross-spectrum transfer. G1 (glycan-Y-first generation) is already VERIFIED this
session and only needs the promotion step folded into G0's exit.

## Hard constraints (apply to EVERY task below — non-negotiable)

- **FDR = Percolator ONLY, never Mokapot.** 2D-FDR (G3′) is a *thin post-process*
  over two vanilla Percolator runs + one inclusion-exclusion merge. No andes-internal
  FDR engine, no finite-mixture FDR. Grep the Percolator PIN mode (Concatenated vs
  Separate) — cross-mode counts are not comparable
  (`feedback_percolator_mode_detection_caveat`).
- **Clean-room, published papers only.** Reference-OK (algorithms, not code):
  **a glyco search engine2/a glyco search engine** (PMC5585273 / PMC8493561; repo Apache-2.0 but *runtime
  license-gated at i.pfind.net* — paper only, never vendor the binary),
  **a cross-spectrum glyco engine** (github.com/DICP-1809/a cross-spectrum glyco engine, **Apache-2.0** — the G4
  reference), **GlycReSoft** (PMC11263600, Apache-2.0), **O-Pair/an open-source glyco engine**
  (**GPL-3.0** — cite PMC7606753 / PMC8933705, do **NOT** copy code into Apache
  andes). **FORBIDDEN code:** **a commercial glyco engine** (commercial, Protein Metrics/Dotmatics) and
  **the reference glyco engine/FragPipe** (UM-proprietary academic-only) — usable for
  labels/notation cross-checks only, never algorithms.
- **Additive-only PIN features.** Y0/Y1 anchor and all glyco evidence are *additive*
  PIN columns; never fuse glycan evidence into the peptide ranking score, never
  modify an existing PIN feature (`feedback_piecewise_alignment_doesnt_work`).
- **Kill-gate on decoy-separated top-1 ranking, not find-rate.** Generation is
  solved; the honest metric for every scoring change is target/decoy *separation* on
  the peptide axis.
- **Differentiate, don't clone:** andes = glycan-Y-first generation + own learned
  regime-matched model + in-process RT-gated cross-spectrum + Percolator-native
  2D-FDR. Not a re-implementation of the reference engine/a comparison search engine fragment-index open search.

## The measurement harness (used by every gate — read once)

- **Fast dev harness:** `ANDES_GLYCO_SCANS=<file>` subsets the driver to the truth
  scan numbers (one per line), ~8 min/arm vs ~100 min full
  (`glyco_search.rs:216`, use `:308`). This is the A/B rig for **every** gate below.
- **Truth set:** PXD025455 `HCC_pool_Late_Fc3_r1`, **523** N-glyco truth scans
  (own the reference engine-nglycan re-search). Frozen EVAL holdout — never train on any
  PXD025455 file (roadmap §5). Baseline chain to beat:
  `154/523 recovered @1% → 83/154 top-1 correct → 66 true (111 false pass)`
  (`00-current-state.md §1`).
- **One variable per A/B; one-commit=one-truth; verify provenance** (binary commit /
  model SHA / truth row count) before trusting any number
  (`feedback_experiment_hygiene`). Flat / byte-identical A/B = RED FLAG (toggle not
  wired), investigate before reporting.
- **Percolator is the FDR engine** for recovery numbers; the fast harness alone gives
  find-rate + top-1, which is enough for G0/P0 gates. Recovery@1%FDR gates (P1/G3′/G4)
  need a Percolator run on the emitted `.glyco.pin`.

---

## Phase G0 — Correctness hygiene (do first; no gate on IDs, gate on invariants)

**Objective.** Remove known defects so every later A/B measures signal, not bugs.
Pin the mass/water convention. Fold in the G1 promotion decision. Nothing here should
*change* IDs except by removing nondeterminism; if an item moves find-rate, that is a
latent bug being fixed and must be measured.

**Tasks (checklist).**

- [ ] **DET-1 — nondeterministic truncation sort.** `hybrid.rs:318` sorts backbone
  hits by `core_y_intensity` with `partial_cmp(&a.0).unwrap_or(Ordering::Equal)` then
  a `to_bits()` mass tiebreak. The intensity `partial_cmp().unwrap_or(Equal)` silently
  treats any NaN as a tie → truncation order can jitter. Change the primary key to
  `b.0.total_cmp(&a.0)` (the mass tiebreak at `:319-320` is already a total order;
  keep it). Note: the *selection* axes in `glyco_search.rs:621-656` already use
  `total_cmp` — this is the one remaining `partial_cmp` on a hot sort path.
- [ ] **P0.4 — probe isotope fidelity.** Audit `hybrid_candidates_with_isotope`
  (`hybrid.rs:159`) + the isotope sweep in the driver (`glyco_search.rs:333-450`):
  confirm each isotope offset uses `proton 1.007276` (not neutral-H 1.00783) and the
  C13 step (`1.003355`) — a wrong constant is ~0.5 mDa/charge drift that flips
  near-isobaric compositions (`30-standards/masses.md`).
- [ ] **P0.3 — Y0/Y1-only quorum retention (MEASURE BEFORE SHIPPING).** The de-novo
  solver quorum (`backbone.rs:109` `solve_backbone_min`, rescue at `hybrid.rs:200-203`)
  may drop backbones supported only by Y0/Y1. Y0/Y1 are the *one* peptide-specific
  anchor (P0/SP-B lever) — dropping them pre-scoring wastes the signal. Add Y0/Y1 to
  the retention quorum, then **A/B find-rate on the 523**: ship only if find-rate is
  unchanged-or-up. This is a de-risk for P0, not an independent win.
- [ ] **H2O convention — confirm pinned, do not re-derive.** Glycopeptide neutral =
  `M_peptide (Σresidues + one H2O + mods) + glycan_residue_sum`. Glycan is a *single
  delta on the sequon residue*; **never add +18.010565 to an attached glycan** (that
  water is only for free/released glycans). Verify `write_glyco_psm_row`
  `CalcMass = peptide.mass() + glycan_mass` (`glyco_pin.rs:164`) and the backbone
  scorers add `+ H2O` to the *peptide* backbone once and only once (e.g.
  `hybrid.rs` `backbone_mass + H2O` at the `core_y_intensity` call). Residue masses
  are the 6-decimal table in roadmap §4 — assert against it in a test, don't retype.
- [ ] **G1 promotion decision (fold in here).** G1 (`ANDES_GLYCO_YINDEX`) is VERIFIED
  (find-rate 59.3%→69.8% @0.05 Da, +7.2 pts @0.005 Da). Decide the default flip *after*
  DET-1 lands (so the A/B is on a deterministic build). Keep peptide-first
  (`ANDES_GLYCO_PEPTIDE_FIRST`, `glyco_search.rs:200`) as a **fallback branch only**,
  never the spine (roadmap §7 anti-goal). Promotion is a one-line default change +
  a re-run of the 523 A/B; do not promote if it regresses eval.

**GATE (invariants, not IDs).**
1. All 58 existing glyco tests green (`crates/andes-glyco/src/`, `glyco_search.rs`;
   count via `grep -rn "#\[test\]"`).
2. Find-rate A/B on the 523 truth scans (`ANDES_GLYCO_SCANS`) **unchanged or up**
   for DET-1 and P0.4; **measured and non-regressive** for P0.3 before it ships.
3. **Bit-reproducible:** two full invocations on the same input produce
   byte-identical `.glyco.pin` (the DET-1 acceptance test). Diff sorted PIN.

**Risks / kill-criteria.**
- P0.3 *reduces* find-rate (Y0/Y1 quorum admits noise backbones that crowd
  truncation) → **do not ship P0.3**; leave the quorum as-is and rely on P0's anchor
  feature at scoring time instead. This is a measured no-op, not a failure.
- DET-1 changes IDs (should not) → a NaN was previously masking a real ordering bug;
  investigate before proceeding (`superpowers:systematic-debugging`).

**TDD test ideas.**
- `det1_truncation_is_total_order`: build two `Vec<BackboneHit>` differing only in
  push order with equal `core_y_intensity`; assert the truncated set is identical
  (deterministic) across both orders.
- `attached_glycan_no_extra_water`: assert
  `glycopeptide_neutral(peptide, glycan) == peptide.mass() + glycan_residue_sum`
  to ≤1 mDa for a known GPSM (e.g. `Hex(5)HexNAc(4)` on a fixture peptide);
  regression-guards the −18.0106 double-count.
- `residue_masses_match_standard`: table-assert the 12 residue constants
  (`glycan_mass.rs:5-17` + roadmap §4) to 6 decimals.
- `full_run_is_bit_reproducible`: run the fast harness twice on a 5-scan fixture,
  assert `.glyco.pin` bytes equal.

---

## Phase P0 — Y0/Y1 anchor feature + SP-B kill-gate  ★ cheap decision point

**Objective.** Add the one peptide-specific PIN feature and **decide whether a
learned model can rank the peptide axis at all** *before* spending the harvest/retrain
(P1). The honest premise (`00-current-state.md §6`): retraining alone is likely
insufficient because stepped-HCD b/y is physically ~11% sparse. Y0/Y1 is the only
recoverable peptide-mass-conditioned signal (`why-andes-fails-and-succeed.md §3`).

**Tasks (checklist).**

- [ ] **Add `Y0Y1AnchorScore` — an additive-only PIN column.** Compute at phase-2
  feature extraction (`glyco_search.rs:691-778`, alongside the `GlycoPsmKey` build at
  `:752`). Definition (`30-standards/masses.md` Y-ion convention):
  - `Y0` neutral = `M_peptide`; observed 1+ m/z = `M_peptide + 1.007276`.
  - `Y1` neutral = `M_peptide + 203.079373` (one innermost HexNAc); 1+ adds a proton.
  - Feature = peptide-mass-conditioned matched intensity of Y0 and Y1 (and their 2+
    charge states) at the model fragment tolerance. Because it is a function of the
    *peptide* mass, it discriminates competing peptides at one backbone window — which
    oxonium/YLadder/CoreYHits/GlycanMass cannot (`why-andes-fails-and-succeed.md §2`).
  - **Additive only:** new column in `GlycoPsmKey` (`glyco_psm.rs:36-59`) + header in
    `write_glyco_header` (`glyco_pin.rs:107-112`). **Never** fold it into
    `score_psm`/`psm_edge_score` at `glyco_search.rs:577-578` (the ranking score
    stays the intact backbone b/y model). Modifying an existing feature regresses
    Percolator (`feedback_parity_tuning_lessons`).
- [ ] **Build the decoy-separated kill-gate harness.** For each of the 523 truth
  scans, at the true backbone mass window, score TWO candidates with the *same*
  backbone-mass hypothesis: (a) the true peptide, (b) a **reversed-peptide decoy**
  (same residue multiset, sequon fixed by N-X-S/T). Emit
  `(backbone_b/y_score + Y0Y1AnchorScore)` for both. This reuses the search's own
  reversed-decoy generator; the gate lives in the fast harness, not a new binary.
- [ ] **Measure separation, not find-rate.** Report: (1) top-1 correct count on the
  523 with the anchor feature added to the PIN and a Percolator run (target: 83 → up);
  (2) target-vs-decoy separation on the kill-gate — e.g. AUROC / Mann-Whitney of
  `score_true` vs `score_reversed`, and the fraction of scans where true > decoy.

**GATE (go/no-go — this is the cheap decision).**
- **GO** if the anchor feature **both** lifts top-1 (83 → measurably higher on the
  523 via Percolator) **and** gives measurable target/decoy separation on the
  kill-gate (true > reversed materially above 0.5 fraction / AUROC materially > 0.5).
  → proceed to **P1** (harvest + retrain as the calibration layer the anchor rides on).
- **NO-GO** if separation stays flat (true ≈ reversed) — the honest-premise outcome.
  → **skip P1's retrain-for-ranking rationale and jump to G4** (cross-spectrum
  transfer), because no single-spectrum model can rank a signal that isn't there. P1's
  wiring (protocol variant) may still be worth doing as infra, but not as *the ranking
  fix*.

**Risks / kill-criteria.**
- Anchor fires but does not separate (Y0/Y1 present for *both* competing peptides at
  near-isobaric backbone mass) → separation flat → **NO-GO**, this is the designed
  kill outcome, not a bug. Report it and pivot to G4.
- Anchor collinear with an existing feature → Percolator flat @1% despite raised raw
  separation (`recent-ideas-impl-divergence`) → check the PIN feature correlation;
  still additive-safe, but the *gate* is top-1 + kill-gate separation, not Δrecovery.
- Reversed-peptide decoy accidentally re-generates a real tryptic peptide (palindrome
  / near-palindrome) → contaminates the gate; dedup decoys against the target DB.

**TDD test ideas.**
- `y0y1_anchor_targets`: for a fixture GPSM, assert Y0 m/z = `M_pep + 1.007276` and
  Y1 m/z = `M_pep + 203.079373 + 1.007276` to ≤1 mDa; assert the feature is 0 when
  neither ion is present and >0 when a synthetic Y0 peak is injected.
- `anchor_is_additive_only`: golden-diff the *ranking* score (`score_psm+edge`) with
  and without the anchor feature enabled — must be byte-identical (proves no fusion
  into ranking).
- `anchor_discriminates_competing_peptides`: two peptides at equal backbone mass,
  distinct sequence; inject Y0/Y1 for peptide A only; assert
  `anchor(A) > anchor(B)` (the peptide-axis discrimination the SPA2 features lack).
- `reversed_decoy_preserves_backbone_mass`: assert reversed-peptide decoy mass ==
  target mass (same multiset) and still contains an N-X-S/T sequon.

---

## Phase P1 — Harvest + regime-matched retrain (SP-B model; the calibration layer)

**Objective.** Give the peptide axis a learned, stepped-HCD-matched model under
`protocol=NGlyco` — the calibration the P0 anchor feature rides on. Runs **only if P0
= GO**. Reuse the existing `protocol=` model-store infra (~90% reusable,
`00-current-state.md §5`).

**Tasks (checklist).**

- [ ] **(a) `Protocol::NGlyco` wiring.**
  - Add `NGlyco` to `model::protocol::Protocol` (`crates/model/src/protocol.rs:4`)
    with `name()`/`from_name()` arms (`:14`, `:26`). Keep N first (andes is N-glyco;
    `sequon.rs:13` is N-X-S/T only). O-glyco = a later second variant.
  - Add `#[clap(name = "N-glyco")] NGlyco` to the CLI `Protocol` enum
    (`andes.rs:69`) and map it through the CLI→model `Protocol` conversion.
    `--glyco` (`andes.rs:391`) stays the mode flag and forces `protocol=NGlyco`
    internally (cleaner minimal change than a new `--protocol N-glyco`).
  - Add the `build_selection_key` arm (`andes.rs:4809`): `Protocol::NGlyco =>
    "NGlyco"` as `protocol_for_store`.
  - Add the `protocol_to_experiment_class` arm (`store/read.rs:255`):
    `"NGlyco" => parse_experiment_class("glyco")` — single opaque slug, mirroring the
    `iTRAQPhospho` precedent at `:262`. The `"glyco"` catalog slug already exists
    (`catalog.rs:95`, `inference: None`). `select_nearest` (`select.rs:261`) already
    WARN-degrades to the standard base when `protocol=NGlyco` is absent → incremental
    rollout is safe (variant added before model trained = today's behaviour).
- [ ] **(b) Glycan-stripped backbone training rows — THE real modeling decision.**
  Training accumulates from labeled `(peptide, spectrum, charge)` rows via
  `accumulate` (`accumulate.rs:62` → `ScoredSpectrum::new` + `ion_match_facts` at
  `:73,:80`), a bare-peptide matcher. For glyco this is correct *iff* the peptide is
  the **glycan-stripped backbone** (Asn carries NO glycan mass) and the spectrum is
  the raw stepped-HCD glyco MS2. **The glycan mass must not shift backbone b/y.**
  - Sub-task 1: emit labeled rows `(backbone_peptide, glyco_spectrum, charge)` with
    the bare backbone.
  - Sub-task 2 (optional, higher-value): register glyco Y-ion offsets as
    `loss_class=1` (schema already reserves `1=glyco` at `store/schema.rs:211`,
    columns `ion_loss_class :213` + `frag_off_loss_classes :222`; read/write paths
    exist at `store/read.rs:491-558` / `store/write.rs:393-470`) so the model *learns*
    the trimannosyl-core ladder instead of penalizing it as missing. **No producer
    exists today** — this is new `train`-side code, gated by the existing slug.
- [ ] **(c) Truth-TSV → flat-parquet converter (new script).** Turn
  `(scan, backbone_peptide, glycan, charge)` truth rows into the labeled flat-parquet
  `accumulate` consumes (every row a label, no bootstrap search). No converter exists
  in-tree (`store/write.rs` writes the *trained model*, not labeled input). Join key =
  `(raw_file_stem, native scan number)`; parse `scan=N` from the spectrum id (mzML
  index ≠ scan; fold I→L for the join key only, keep original in the report).
- [ ] **(d) Harvest on Codon + tiered labels + retrain.** Use the `codon-cluster`
  skill. Corpus (roadmap §5, all PRIDE-verified, **mixed species = anti-leakage**):
  PXD005411 (mouse brain, a glyco search engine2 ref, sample already at
  `40-data/collection/psms_pxd005411.tsv`), PXD016175 (human IgG plasma), PXD030670
  (human saliva QE — closest instrument to eval), PXD020254 (Lumos stepped-HCD).
  Demote PXD011239 (EThcD, wrong fragmentation) to train-only; **drop** PXD057219
  (venom, no deposited glyco results). **Tiered labels, NOT pure consensus:**
  Tier A = multi-engine consensus (peptide + glycan-composition agreement) = gold;
  Tier B = single-engine Y-ion/oxonium-confirmed **hard cases** (consensus alone
  excludes the sparse-b/y cases the model exists to rank). Target ≥100–300k
  glyco-PSMs. Retrain the own `strong` spectral model; write
  `resources/models/protocol=NGlyco` via `split_store_by_protocol`.

**GATE.**
- Retrained `protocol=NGlyco` model + P0 anchor **beats the field-default model on
  eval top-1** on the 523 (via Percolator on `.glyco.pin`) **AND** clears the P0
  kill-gate on **held-out** data (a PXD025455 file never in any harness A/B — leakage
  guard).
- **Entrapment-FDP honest** (roadmap §6, a glyco search engine2 validator PMC5585273): pad the search
  DB with a foreign glycome+proteome; any GPSM with a foreign-only glycan OR
  foreign-only peptide is a de-facto false positive → the decoy-estimated FDR must
  not under-cover the entrapment FDP.

**Risks / kill-criteria.**
- Retrained model ≈ field-default on eval top-1 → confirms P0's honest premise
  (calibration can't manufacture absent b/y) → **stop treating P1 as the ranking fix;
  the win is P0 anchor + G4**. Keep the wiring (it's needed for G4/G3′ routing) but
  don't chase model tweaks.
- **Do NOT prefer GBDT for SP-B.** GlycReSoft tested gradient-boosted trees, saw
  substantial overfitting, and chose a regularized multinomial/log-linear model
  (PMC11263600). Regime-match / regularize, not GBDT (roadmap §7).
- **Do NOT copy static hand-tuned weights** (a glyco search engine `w≈0.35`, a commercial glyco engine
  `w_pep ≫ w_Y ≫ w_oxo`, a cross-spectrum glyco engine `freq^0.3`, O-Pair Table S1). Borrow the score
  *decomposition* (intensity × quartic-mass-error × ion-ratio × core-ratio); **learn**
  the combiner regime-matched.
- Label leakage: any PXD025455 file in the train corpus → invalidates the eval.
  Assert the harvest manifest excludes all PXD025455 files (roadmap §5).
- NeuGc/NeuAc confusion in mouse data (PXD005411): NeuGc 307.0903 ≠ NeuAc 291.0954
  (16.00 Da). Ensure the converter maps species-correct sialic acids or sialylated
  masses silently mis-join.

**TDD test ideas.**
- `nglyco_protocol_roundtrips`: `Protocol::from_name("N-glyco").name() == "N-glyco"`;
  `build_selection_key(NGlyco).protocol_for_store == "NGlyco"`;
  `protocol_to_experiment_class("NGlyco")` contains `"glyco"`.
- `select_nearest_degrades_when_glyco_absent`: with no `protocol=NGlyco` partition,
  selecting NGlyco WARN-falls-back to the standard base (asserts safe rollout).
- `training_rows_are_glycan_stripped`: for a truth GPSM, the emitted training-row
  peptide mass == bare backbone mass (Asn carries no glycan); assert no +glycan delta
  leaked into the labeled peptide.
- `truth_tsv_join_key`: `parse_scan("controllerType=0 controllerNumber=1 scan=12345")
  == 12345`; I→L folding applied to the join key only, original preserved.
- `entrapment_foreign_only_is_fp`: a GPSM whose glycan is in the foreign padding only
  is flagged FP by the validator.

---

## Phase G3′ — 2D-FDR post-process (runs WITH P1/G4, not after)

**Objective.** Valid glycopeptide FDR without a unified decoy pile and without a
second FDR engine. Two vanilla Percolator runs + one inclusion-exclusion merge
(`20-theory/glyco-fdr.md`). This is the fix for the **111 false passers** and the
refuted-unified-pile crash (29.4%→4.4%, `00-current-state.md §2`).

**Tasks (checklist).**

- [ ] **Peptide-axis Percolator run.** Emit ONE PIN with **reversed-peptide decoys
  only**; `Label` = peptide target/decoy. Features = backbone b/y
  (`RankScore`/learned `ScoreP`) **plus the P0 `Y0Y1AnchorScore`**. Run Percolator →
  per-PSM `score_P`, `q_P`.
- [ ] **Glycan-axis Percolator run.** Emit a SEPARATE PIN whose decoys are the
  **Y-rung-shifted glycan decoys** — shift all peptide+Y ions **EXCEPT Y0/Y1** (keep
  the peptide-mass anchors intact so the decoy competes on the same backbone), or
  per-fragment random 1–30 Da (`20-theory/glyco-fdr.md §3`, a glyco search engine2 recipe
  PMC5585273 / a cross-spectrum glyco engine PMC8990002). `Label` = glycan target/decoy; features =
  glycan-structure evidence only (`YLadderScore`, `OxoniumScore`, `CoreYHits`). The
  glycan-decoy Y-ladder scorer (`backbone.rs:498` `glycan_y_intensity_decoy`, splitmix64
  seed `:479`) and paired decoy PIN rows (`glyco_pin.rs:320-330`, `glycandecoy_`
  accession `:247`) already exist and are TDD-sound (target > decoy). Run Percolator →
  `score_G`, `q_G`.
- [ ] **Inclusion-exclusion merge (the only new logic).** On the joint ranking (e.g.
  `w·score_G + (1−w)·score_P`, `w` *learned*, not the a glyco search engine 0.35 constant), at each
  cut compute `FDR_P`, `FDR_G`, and `FDR_{P∩G}` (GPSMs flagged decoy on *both* axes),
  then `FDR = FDR_P + FDR_G − FDR_{P∩G}`; accept GPSMs with `FDR ≤ 0.01`. Thin
  post-process reading the two Percolator outputs — no andes-internal FDR.

**GATE.**
- Combined 2D-q ≤ 1% **holds against the entrapment-FDP** (roadmap §6) — decoy FDR
  must not under-cover the foreign-padding FDP.
- **Recovery does NOT collapse** — contrast the refuted unified pile (29.4%→4.4%).
  Recovery on the 523 at valid 2D-q ≤ 1% should exceed the 66-true baseline and track
  the top-1 correct count (target: kill the 111 false passers without dropping true
  IDs).

**Risks / kill-criteria.**
- **Do NOT feed glycan-decoy rows into one unified `Label` pile** (roadmap §7). They
  differ only in YLadder → flood the −1 pile → Percolator over-weights YLadder →
  recovery crash. Two separate axes, always. Re-check the mode (Concatenated vs
  Separate) with a PIN grep.
- Glycan decoy accidentally shifts Y0/Y1 → destroys the peptide anchor on the glycan
  axis and double-penalizes → assert Y0/Y1 are held fixed in the glycan decoy.
- Combined q under-covers entrapment FDP → the merge or a decoy recipe is
  mis-specified → **do not ship**; the entrapment gate is the honesty check.

**TDD test ideas.**
- `glycan_decoy_preserves_y0_y1`: assert the Y-rung-shift leaves Y0 and Y1 m/z
  unchanged and shifts ≥1 interior Y-rung (the recipe invariant).
- `glycan_decoy_scores_below_target`: regression-lock `glycan_y_intensity_decoy <
  glycan_y_intensity` for a fixture (already covered — extend to the shifted-decoy PIN
  row).
- `inclusion_exclusion_bounds_union`: with synthetic `(FDR_P, FDR_G, FDR_{P∩G})`,
  assert `FDR = FDR_P + FDR_G − FDR_{P∩G}` and monotone non-decreasing along the cut.
- `two_axis_pins_are_separate`: assert the peptide PIN carries no glycan-only feature
  columns and the glycan PIN carries no peptide-anchor decoy rows (no cross-contam).

---

## Phase G4 — In-process cross-spectrum transfer  ★ ceiling-breaker

**Objective.** Break past the ~11% single-spectrum b/y ceiling — andes's unique
differentiator. Activate the existing scaffold; add **RT gating** + cosine-weighted
backbone-evidence transfer from confident donor PSMs to unassigned same-backbone
spectra (a cross-spectrum glyco engine mechanism, Apache-2.0 reference, +33.5–178.5% PSMs from
*transfer, not scoring*, PMC8990002). Pays only once P0/P1 lift donor-glycoform top-1
confidence.

**Tasks (checklist).**

- [ ] **Activate + RT-gate the scaffold.** The two-pass design is COMPLETE and gated
  OFF (`ANDES_GLYCO_CROSSSPECTRUM`, `glyco_search.rs:203`; early return `:797-799`
  skips pass 2; whitelist `crossspectrum.rs:25` `GlycoformWhitelist`, `transfer :58`,
  `nearest_glycan :85`; 3 unit tests `:119-155`). Pass 1 (`:791`) = normal generation;
  build the whitelist from pass-1 PSMs with `core_y_hits ≥ CONF_MIN_CORE_Y=3`
  (`:804-814`); pass 2 (`:831`) injects transferred backbones into non-confident
  oxonium-positive spectra via `whitelist.transfer` (`:856`) through the same
  dedup/score path (`process_one(..., &transfer)` at `:876`); pass-2 supersedes pass-1
  (`:880-882`). **New work:** add an **RT window** so transfer only fires between
  co-eluting donor/acceptor spectra (avoid transferring across co-eluting
  confounders, `why-andes-fails-and-succeed.md §4`); weight the transferred
  backbone-fragment-frequency prior by donor-acceptor **cosine similarity**.
- [ ] **Truth A/B on the 523** (never before P0/P1 lift donor confidence — transfer
  from a wrong donor propagates error). Measure net glyco-PSM gain at 2D-q ≤ 1%.

**GATE.**
- **Net glyco-PSM gain** on the 523 eval (target: recover a meaningful fraction of the
  sialylated / short-peptide sparse-b/y stratum that single-spectrum scoring cannot
  rank) **with G3′ 2D-FDR held at 1%**. Gain must be from transfer (ablate: transfer
  ON vs OFF, one variable).

**Risks / kill-criteria.**
- Transfer from a low-confidence / wrong donor → error propagation → recovery drops or
  FDP inflates → **tighten the donor confidence floor** (`CONF_MIN_CORE_Y`) and/or the
  RT window; if gain stays negative at 2D-q≤1%, **kill G4 for this dataset** (the
  scaffold stays opt-in).
- No RT gate → transfer across co-eluting confounders inflates false PSMs (this is why
  the NULL-v1 scaffold is gated OFF). RT gating is mandatory, not optional.
- Gain not entrapment-validated → run the §6 entrapment FDP with transfer ON.

**TDD test ideas.**
- `transfer_respects_rt_window`: two donor/acceptor spectra outside the RT window →
  `whitelist.transfer` yields no transferred backbone; inside the window → it does.
- `transfer_requires_confident_donor`: a pass-1 PSM with `core_y_hits < CONF_MIN_CORE_Y`
  does not enter the whitelist.
- `transfer_cosine_weighting_monotone`: higher donor-acceptor cosine → higher
  transferred-prior weight (monotone).
- `pass2_supersedes_pass1`: an acceptor spectrum re-ID'd in pass 2 replaces its pass-1
  (non-)assignment deterministically.

---

## Cross-phase dependency + gate summary

| Phase | Depends on | Gate (one line) | Kill → pivot |
|---|---|---|---|
| **G0** | none | 58 tests green · find-rate ≥ · bit-reproducible | P0.3 regresses → drop P0.3 |
| **P0** (SP-B kill-gate) | G0 | top-1 83→up **AND** decoy-separated separation on 523 | flat separation → skip P1-as-fix, jump G4 |
| **P1** (SP-B model) | P0 = GO | NGlyco model+anchor beats default top-1 on held-out **AND** entrapment-honest | ≈default → keep wiring, win is P0+G4 |
| **G3′** (2D-FDR) | P0 (anchor) | 2D-q≤1% holds vs entrapment · recovery no-collapse | unified pile / Y0Y1-shift → don't ship |
| **G4** (cross-spectrum) | P1 (donor confidence) | net PSM gain at 2D-q≤1%, transfer-attributable | negative gain → keep opt-in, kill for dataset |

## Cited sources

- Fang et al., **a cross-spectrum glyco engine**, *Nat Commun* 2022 — **PMC8990002** (cross-spectrum
  transfer; ~11% direct-b/y ceiling; +33.5–178.5%). **Apache-2.0**,
  github.com/DICP-1809/a cross-spectrum glyco engine — the G4 reference.
- Liu et al., **a glyco search engine 2.0**, *Nat Commun* 8:438, 2017 — **PMC5585273** (2D
  inclusion-exclusion FDR; Y-rung-shift glycan decoy; entrapment validator). PXD005411.
- Zeng et al., **a glyco search engine**, *Nat Methods* 2021 — **PMC8493561** (glycan-Y-first ion
  index; biosynthetic DB). Repo Apache-2.0, **runtime license-gated (i.pfind.net)**.
- Lu et al., **O-Pair**, *Nat Methods* 2020 — **PMC7606753**; Klein & Zaia,
  multi-attribute glycan FDR, *MCP* 2022 — **PMC8933705**. an open-source glyco engine **GPL-3.0**
  (paper-cite only, do NOT copy code).
- **GlycReSoft** — **PMC11263600** (Apache-2.0; separated FDR axes; GBDT-overfit
  finding → use regularized log-linear, not GBDT).
- ppmFixer, *Glycobiology* 2024 (6-decimal mass necessity).
- Repo internal: `00-context/00-current-state.md`, `00-context/02-code-inventory.md`,
  `20-theory/{glyco-fdr,why-andes-fails-and-succeed,fragmentation-why-hard}.md`,
  `30-standards/{masses,notations}.md`, `40-data/*`, `PHASE1_RESULT.md`,
  `SPA2_RESULT.md`.
