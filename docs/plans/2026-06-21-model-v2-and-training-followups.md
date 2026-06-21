# Plan — Model v2 decision + training/independence follow-ups

**Date:** 2026-06-21
**Status:** proposed (awaiting go)
**Scope:** resolve the v2 retrain ship/no-ship decision, and the three engineering tracks (silent train failures, PROTON→CODATA, partition-structure re-derivation).

---

## 0. Context & baseline

The v2 retrain was assembled (`own_models_v2.parquet`, 39 models, version 2, zero MS-GF+ `7061` rows) and benchmarked on the VM vs the **shipped v1 bundle** (which is *itself* already own-trained — no `7061`), at 1% true entrapment-FDP:

| Dataset (model exercised) | v2 | v1 (shipped) | verdict |
|---|--:|--:|---|
| TMT a05058 (`cid_lowres_tryp`) | 11,341 | 11,335 | ~par |
| UPS1 5000amol (`cid_lowres_tryp`) | **16,655** | 16,507 | **+148 win** |
| Astral LFQ (`hcd_qexactive_tryp`) | 21,775 | **32,077** | **−32% regression** |

The Astral regression is real (lower FDP for v2: 0.74% vs 0.97%; more candidate rows but fewer confident PSMs → the retrained `hcd_qexactive_tryp` tables are noisier). **Conclusion: do not ship v2 wholesale.** The `cid_lowres_tryp` retrain is good; the `hcd_qexactive_tryp` retrain is bad.

**★ VERIFIED 2026-06-21 — 2 of the 39 shipped models are still byte-identical MS-GF+ seed (version-masked).** A stamp-independent table-value hash of `resources/models.parquet` vs the seed (`7061`) shows `hcd_qexactive_tryp_phosphorylation` and `hcd_highres_tryp_phosphorylation` (hash `0180be4c`, 371 tables) are the **original MS-GF+ phospho tables, restamped version 7061→1** (silent-train-failure → seed fallback). The other 37/39 have own table values. → **The version stamp is NOT proof of independence**; this elevates E1 (esp. the fail-loud failsafe) to required, and adds a hard gate: **retrain or drop these 2 phospho models before any "own models" claim.** (Table-value independence = 37/39; geometry independence = 0/39 until E3.)

**Two code facts established this session (grounding the plan):**
- `model-train/src/estimate.rs:109-112` — the trained `Param` **copies** `num_segments` + ion-type/frag-offset **layout** from the seed template; only the table *values* are learned. → partition **geometry** is still MS-GF+-shaped in every "own" model.
- `model-train/src/accumulate.rs` + `estimate.rs` — `accumulate()` runs only on accepted PSMs; `estimate()` falls back to the template when counts are empty. → 0 confident bootstrap PSMs ⇒ **silent** seed copy (the train-failure mechanism).

---

## Part 1 — Models: one sequenced path (resolves options a / b / c)

Don't pick one of a/b/c in isolation — sequence them by payoff/risk. (a) keep-v1 is the always-available safety net; (b) v3 banks the only real win now; (c) is the deeper fix that can later beat v1 on Astral too.

### Phase M0 — Ship v3 hybrid *(option b; immediate, ~1–2h)*
Bank the `cid_lowres_tryp` win with **zero** Astral risk.
- **Assemble `own_models_v3.parquet`** = new v2 `cid_lowres_tryp` rows (TMT par + UPS1 +148) **+ v1 `hcd_qexactive_tryp`** rows (Astral 32,077) **+ v1 for all other bases**.
  - Mechanism: in the assembler, source the `hcd_qexactive_tryp` base (and the other untouched bases) from the **shipped v1 `resources/models.parquet`**, and only `cid_lowres_tryp` from the v2 store. (Per-base store-selection map; the existing `MAP` already routes slugs→bases.)
- **Validate**: assembler invariants (39 models, no `7061`) + VM re-benchmark TMT/UPS1/Astral at 1% entrapment-FDP.
- **Gate to ship**: Astral ≥ v1 (32,077) **and** UPS1 > v1 (16,507) **and** TMT ≥ v1 (11,335). If any fails unexpectedly → fall back to **keep v1** (option a), no harm.
- **Output**: a shippable bundle that is strictly ≥ v1 on every dataset.

### Phase M1 — Diagnose the Astral regression *(option c; parallel deeper track, ~1–2 days)*
Systematic-debugging (root cause **before** any fix). Symptom: more candidates, fewer confident PSMs → noisier rank/intensity tables for the v2 `hcd_qexactive_tryp`.

Hypotheses, tested one variable at a time (re-train the single slug, benchmark Astral each time):
- **H1 — corpus composition.** The v2 corpus added lower-quality / off-regime PSMs that diluted the tables. → Diff the v1 vs v2 `hcd_qexactive_tryp` corpus manifest (PXDs, #PSMs, instrument). Re-train on **v1's PXD subset only**; if Astral recovers → corpus dilution is the cause.
- **H2 — over-smoothing / backoff.** More PSMs but flatter tables (Laplace smoothing + segment-collapse). → Compare `rank_dist_table` sharpness (entropy) v1 vs v2; try tighter `min_count` / smoothing.
- **H3 — noise calibration.** `ANDES_DENSE_NOISE` differed → missing-ion penalty miscalibrated (the accumulate.rs noise comment: wrong noise shape ⇒ "0 PSMs at 1% FDR" failure mode). → Re-train with matched noise setting.
- **H4 — labeling.** Bootstrap `train_fdr` admitted more borderline PSMs → noisier facts. → Tighten train_fdr; re-train.
- **Deliverable**: the variable that caused −32% + a **v4 `hcd_qexactive_tryp`** that beats v1 on Astral.

### Phase M2 — Ship the improved Astral model *(gated on M1)*
If M1 yields a `hcd_qexactive_tryp` that beats v1, fold it into the bundle (v4) and re-benchmark all 3.

---

## Part 2 — Engineering follow-ups

### Track E1 — The 2 silent train failures *(cid_lowres_tryp_tmt, hcd_highres_nocleavage_phosphorylation)*
**Symptom**: `train-from-msnet` accumulates 0 PSMs → seed store, while logging "trained OK".
**Grounded hypothesis** (from code): the bootstrap seed search (`labeled.rs::bootstrap_labels`) accepts 0 confident target PSMs at `train_fdr` for these slugs — most likely because the MSnet flat's modifications (TMT6plex; phospho on STY) aren't applied when building the search peptides → peptides searched unmodified → mass mismatch → no matches → empty counts → `estimate()` returns the template (seed).

**Debugging plan (systematic, root-cause first):**
1. **Reproduce**: re-run train-from-msnet for `cid_lowres_tryp_tmt` with verbose logging; confirm **0 accepted PSMs** (vs 0 *facts* — different layer).
2. **Instrument the boundaries**: log #flat PSMs read → #with parseable mods → #peptides built *with* the mod → #confident at `train_fdr`. Find the layer that collapses to 0.
3. **Confirm the mod-mapping hypothesis**: check that the flat's mod tokens (e.g. `TMT6plex`@K/n-term, `Phospho`@STY) map to the slug's fixed/variable mods used in the bootstrap search.
4. **Fix at the source**: apply the flat's mods when constructing the labeled peptides, *or* seed the slug's search params with the correct mods (TMT6plex fixed on K + n-term; Phospho variable on STY).
5. **★ Durable failsafe (most important)**: make `train-from-msnet` **error loudly** when accepted PSMs == 0 (or below a floor), instead of silently copying the seed. This prevents *any* future silent regression and would have caught this immediately.
- **Gate**: both slugs accumulate >0 PSMs and produce `rank_dist`/`ion_err` tables that **differ from the seed** (verify byte-difference vs the template); failsafe unit-tested.

### Track E2 — PROTON → CODATA *(independence #6; ~1–2h, low risk)*
- Change `PROTON 1.00727649 → 1.007276466879` (CODATA) in the constants.
- **Impact**: ~2e-8 Da — below every tolerance → **no PSM flips**, but it shifts exact masses in pinned goldens (PIN ExpMass/CalcMass) and breaks the pinned `PROTON.to_bits()` test.
- **Steps**: change the constant → update the pinned-bits test → regenerate the affected parity goldens → run the full suite → **A/B a sample search before/after to confirm PSM counts are byte-identical** (proves result-neutrality).
- **Gate**: tests green + identical PSM counts on the sample search. Removes the last "MS-GF+ truncation tell" on a non-copyrightable physical constant.

### Track E3 — Partition-structure re-derivation *(Proposal B; the last data-axis independence item; biggest)*
- **Confirmed**: `estimate()` copies `num_segments` + partition/frag-offset **layout** from the seed template → every model still carries MS-GF+ geometry (structural lineage of version `7061`), even with retrained tables.
- **Goal**: derive the geometry from the **own corpus** instead of the seed.
  1. Analyze corpus parent-mass/charge distribution → derive mass-tier boundaries (data quantiles) + empirical charge range.
  2. Derive `num_segments` from data (held-out likelihood: 1 vs 2 vs 3) rather than hardcoded 2.
  3. New `Param` constructor that builds structure from derived geometry (no `template` dependency for layout).
  4. Re-train **all** slugs on the new geometry; re-benchmark all 3 datasets.
- **Risk**: the seed geometry is tuned; changing it can regress. **Gate**: per-dataset entrapment-FDP ≥ current on TMT/UPS1/Astral.
- **Framing**: this is an **independence** deliverable (removes the last seed-copied structure → models become fully own in *both* tables and geometry), not a performance play. Pairs with the retrained tables to let us truthfully claim full model independence.
- **Sequencing**: **last** — needs E1 (a reliable, fail-loud training pipeline) and a clean benchmark baseline first.

---

## 3. Recommended sequencing & dependencies

| Order | Task | Why here | Effort | Risk |
|---|---|---|---|---|
| 1 | **E2** PROTON→CODATA | quick, independent, closes independence #6 | ~1–2h | low |
| 1 | **M0** ship v3 hybrid | banks UPS1 +148, zero Astral risk | ~1–2h | low |
| 2 | **E1** train-failure fix + **failsafe** | unblocks reliable retraining; failsafe is highest-value durable fix | ~0.5–1d | low/med |
| 3 | **M1** Astral regression diagnosis | needs E1's reliable pipeline; can then beat v1 on Astral (→M2) | ~1–2d | med |
| 4 | **E3** partition-structure re-derivation | deep independence item; needs E1 + clean baseline; full retrain+bench | ~3–5d | med/high |

**Parallelism**: E2 + M0 can run together immediately (different systems: E2 local, M0 Codon/VM). E1 should precede M1 and E3.

## 4. Success criteria (definition of done)
- **M0**: v3 bundle benchmarks ≥ v1 on all 3 datasets at ≤ v1 FDP → shipped as `resources/models.parquet`.
- **E1**: both slugs train non-trivially (tables ≠ seed); train-from-msnet errors on 0-PSM (unit-tested).
- **E2**: PROTON is CODATA-sourced; suite green; sample search PSM-count identical.
- **M1/M2**: root cause documented; a `hcd_qexactive_tryp` that ≥ v1 on Astral (or a documented "v1 is the ceiling for this corpus").
- **E3**: geometry derived from corpus; all-3-dataset FDP ≥ current → models fully own (tables + structure).

## 5. Notes / non-goals
- No PR/relicense claim of "full model independence" until **E3** lands (geometry) — until then the honest claim is "own-trained tables on MS-GF+-shaped geometry".
- The Astral regression and the train failures may share a corpus-handling root cause (both are training-pipeline robustness) — investigate E1 and M1 with that in mind.
- andes-full packaging, peak-filter knob, Rescore P5/P6 are out of scope for this plan (separate tracks).
