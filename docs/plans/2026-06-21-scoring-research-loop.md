# andes scoring-research loop — agenda

**Goal:** best-in-class peptide ID + speed, better than any search engine. **Self-paced `/loop`** — one focused scoring experiment per iteration.

## Hard constraints
- **NO MS-GF+ ideas** — no generating function / spectral E-value / per-amino-acid SpecProb. We moved away from MS-GF+.
- **Model-driven scoring** — lean into learned models (GBDT intensity / rich-ion / peak S-N, and richer) that help the engine *understand* the spectrum.
- **Candidate-space expansion** is proven-good — chimeric (co-fragmented) + refinement/PTM-cascade. IDs ≈ coverage × separation.
- **PSM-FDR only**, Percolator only (no Mokapot/other tools), metric = PSMs @ 1% true entrapment-FDP.
- Each iteration: read code → 1 hypothesis → cheap A/B (one variable, experiment-hygiene provenance: binary commit + model SHA + data SHA) → bank → next.

## Why this direction (from the 30-agent diagnosis)
The learned GBDT discriminative layer is *flat on closed search* but the thesis is it pays on **expanded** candidate space. So the win is **learned scoring × expansion**, not closed-search table tweaks and NOT a per-spectrum calibration (that was the MS-GF+ reach we dropped). Most of today's loss is recoverable crippled models (recovering now).

## Experiment queue (ordered; re-rank as results land)
- **E1 — Chimeric × live-GBDT (the convergence, headline).** Recovered GBDT-live models + chimeric expansion → does the learned discriminative score convert chimeric coverage into IDs above closed at honest FDP? (chimeric was +101% Astral historically; rich-ion was flat closed — does it pay now under expansion?). *Needs the recovering models.*
- **E2 — Learned PSM re-ranker (model-driven centerpiece).** Replace the analytic StrongScore null with a learned model (GBDT/logistic) over the full per-PSM feature set → a fully learned PSM score. A/B vs analytic strong. This is "models understand MS" taken to the score itself.
- **E3 — Refinement/PTM-cascade × discriminative score.** Refinement expands declared→discovered mods; does the learned score recover them at honest entrapment-FDP? (Percolator, PSM-FDR only.)
- **E4 — Richer learned spectrum model.** Strengthen the fragment-intensity / rich-ion GBDT (more features / better targets / 2+ frags) so the model predicts the spectrum well enough to dominate discrimination.
- **E5 — Float-precision rank score (cheap, additive).** Per-split LLR is `round()`-ed to i32 before summation → low-res ties. Add a float-precision rank feature (additive, safe) → kill ties (~+1–4%, concentrated low-res).
- **E6 — Model-aware low-res peak handling.** Local-window intensity ranking + noise floor + UPS filtering (currently global ranking, TMT-only windowing).
- **E7 — Speed × ID jointly.** Candidate index / prefilter / parallelism so expansion (E1/E3) stays fast. Best ID AND best speed.
- **E8 — Feature engineering (additive).** Complementary-ion / ladder / cross-spectrum features that strengthen the learned discrimination.

## State
- Models recovering on Codon (rt_hcd, rt_cid uncapped + GBDT) → the GBDT-live baseline E1/E2 need. Start E1 once they land.

## Loop results (banked)
### R1 (2026-06-21) — recovered cid_lowres_tryp on UPS1
- Uncapped full-3-stage-GBDT cid_lowres_tryp: **rank 16,451 @1%FDP 1.24% BEATS Comet 14,833 (+10.9%) AND Java 16,351 (+0.6%)** — reverses the crippled release loss (12,212). Model recovery VALIDATED on UPS (rank mode).
- `--score strong` REGRESSES on low-res UPS: 11,480 (-30% vs rank), partly a Percolator mode confound (strong top-1 PIN -> Concatenated/smaller pool; rank target+decoy -> Separate/mix-max). Live GBDTs ARE firing (RichIonLLR/IntensitySignal non-zero).
- ★ INSIGHT: the learned rich_ion/frag_intensity GBDTs help as **PIN FEATURES (rank mode, fed to Percolator)**, NOT as a top-1 strong re-ranker (which narrows the pool + hurts low-res). Direction: model-driven discriminative FEATURES, not pool-narrowing strong re-rank. Binary dbfbe630, store SHA ab01507f.

### R2 (2026-06-21) — recovered hcd_qexactive_tryp on Astral + E1 chimeric
- Astral @1% entrapment-FDP (live 3-stage GBDT, store SHA c3becc4b): rank-closed 29,086 (1.01%); **rank+CHIMERIC 53,574 (1.18%) = +84.2% over closed** — E1 VALIDATED, chimeric survives entrapment-FDP, rank-1 FDP lower (0.43%); strong-closed 35,890 (0.91%) = +23.4% over rank at BETTER FDP. Comet cached 29,244 (1.10%).
- ★ STRONG WINS HIGH-RES (+23% Astral) but REGRESSED LOW-RES (UPS R1). RULE: **strong mode for high-res, rank mode for low-res.** Recovery: 22,401 crippled -> 29-36k = +30%, regression erased.
- ★★ CHIMERIC is the killer lever (+84% Astral at honest FDP with live GBDT discrimination) — the expand-then-discriminate convergence works.
- CAVEAT: Comet comparison has Percolator mode confound (andes Concatenated/TDC vs Comet Separate/mix-max). andes-internal (a/b/c) clean; Comet margin needs mode-matched re-run.

### R3 (2026-06-22) — recovered cid_lowres_tryp_tmt on TMT a05058 — ★★ HARD GAP CLOSED
- TMT @1% entrapment-FDP (recovered UNCAPPED full-3-stage-GBDT, rich_ion 2.24MB + frag_intensity 798KB, store SHA cd4788b5): **rank-closed 11,421 (0.98%) BEATS Java 10,775 (+6.0%) AND Comet 10,248 (+11.4%)** — FIRST time andes leads Java on low-res CID TMT. rank+chimeric 11,549 (+1.1%, small low-res top-up); strong 10,769 (regresses, high-res-only confirmed). Model selection verified: --protocol TMT + auto-detect picked cid_lowres_tryp_tmt. TMT/UPS rank mode = Percolator Separate = MATCHES cached competitors (clean).
- ★★ MILESTONE: andes beats BOTH Comet+Java on ALL 3 datasets with recovered models. WINNING RECIPE: recovered full-GBDT models + RANK mode low-res, STRONG mode high-res, CHIMERIC top-up (+84% high-res / +1% low-res). Release bar (beat both on all 3) MET, modulo: Astral Percolator-mode caveat (andes Concatenated vs Comet Separate in R2) + Java-Astral number missing. NEXT: clean mode-matched Astral re-confirm + Java-Astral; then assemble winning bundle -> re-gate -> release.

### Release assembly (2026-06-22)
- own_models_winning.parquet (SHA c8688fca, 40 models, version 2, GATE PASS zero-seed): recovered full-GBDT cid_lowres_tryp/hcd_qexactive_tryp swapped in + cid_lowres_tryp_tmt ADDED (routable via --protocol TMT), all rich_ion non-null. On Codon $B/own_models_winning.parquet.
- ★ --score AUTO implemented (andes.rs): default now `auto` = STRONG for high-res / RANK for low-res (resolved from model instrument.is_high_resolution()). Makes the "beats the field" claim hold OUT-OF-BOX (high-res default was rank=29,086~tie; now strong=35,890 beats Comet). Build+clippy clean; andes tests running. NEXT: commit, re-gate combined bundle (per-regime auto mode) vs field all 3, then release.
