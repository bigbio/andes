# andes Glyco Engine — Deep-Review Synthesis (8-agent workflow, 2026-07-07)

> Full synthesis of 7 parallel deep reviews (code model inventory, GBDT/training, scoring internals, broad field glyco survey, RT+intensity+rescoring survey, learning-to-rank survey, data-grounded the reference engine-vs-andes gap) + synthesis. Source of the corrected generation-first roadmap.

This is a synthesis task. I have all seven reviews inline and need to produce a rigorous four-part synthesis. No tool calls are needed to gather more data — the task is to integrate what's given. Let me write the synthesis directly as my final response.

# andes Glyco Engine — Synthesis of 7 Deep Reviews

## PART 1 — Complete andes Model / Capability Inventory

### 1A. What andes HAS (models & learned artifacts)

| # | Model / capability | What it does | Trained on | In glyco SELECTOR? | Source |
|---|---|---|---|---|---|
| a | **Rank model** (`Param`+`RankScorer`) → **RawScore** | Per-fragment intensity-RANK LLR `ln(ion_freq[rank]/(noise_freq[rank]·norm))`, summed over cleavage sites → integer RawScore that RANKS candidates | Public PRIDE PSMs (rank_dist/ion_err/noise_err/ion_existence histograms) | **YES — this IS the selector.** `glyco_search.rs:801` `score_psm(...)` | R1(a), R2, R3 §3 |
| — | **Edge score** (`psm_edge_score`) | Complementary per-bond credit along dominant ion direction; additive | derived from rank tables | **YES** — `rk = sc + ei` at `glyco_search.rs:802-803` | R3 §1, §3 |
| b1 | **Peak GBDT** (peptide-AGNOSTIC signal/noise) | `predict_logit` per peak added to node score by rank; 18 spectrum-only features | LightGBM + isotonic on public data | **YES (indirectly)** — folded into `score_psm` node scores | R1(b1), R2-D |
| b2 | **Frag-intensity GBDT** (peptide-CONDITIONED) | Predicts log-relative b/y intensity → drives IntensitySignal cosine + FragPredExplained/ChanceLLR/TopKObserved; 19 features | LightGBM regressor on public data | **NO — winner-only.** `compute_psm_features` at `glyco_search.rs:976` runs AFTER selection | R1(b2), R2-D, R3 §2-3 |
| b3 | **Rich-ion GBDT** (decoy-aware per-ion LLR) | Σ`predict_logit` over matched b/y → RichIonLLR PIN feature; 23 features | LightGBM logistic, target-vs-decoy-matched ions | **NO — winner-only** (same `:976`) | R1(b3), R2-D, R3 |
| b4 | **Coarse IntensityModel** (fallback table) | Mean/spread log-rel intensity per context key; fallback when b2 absent | Aggregated public PSMs | **NO — winner-only** | R1(b4) |
| c | **Mass calibration** (precursor recal + offset table) | Precursor m/z offset correction before windowing | Sampled pre-pass + trained offset histogram | **DELIBERATELY OFF for glyco** (`glyco_pin.rs:716` `PrecursorCalMode::Off`) | R1(c) |
| d | **Glyco Y-ladder scoring** (`glycan_y_intensity`, `core_y_intensity`) | Hand-written summed matched-intensity Y-ladder / base_peak; **no learned weights** | Heuristic (not trained) | **YES** — Y-ladder is the PRIMARY sort key (`y_primary=true`, `glyco_psm.rs:53`) | R2-D, R3 §3, R4 |
| e | **Oxonium gate** (`oxonium.rs:48`) | Fixed fractional threshold on diagnostic oxonium ions | Heuristic (not trained) | YES (gate, pre-scoring) | R2-D |
| f | **y0y1_anchor_score** | Y0/Y1 core anchor score | Heuristic | **NO — winner-only** (`glyco_search.rs:1047`) | R3 §3 |
| g | **RT transfer** (`transfer_rt_delta`) | Raw RT-delta heuristic across glycoforms | Heuristic, NOT a learned RT model | Not a selector score | R2-D |
| h | **Rescorer GBDT** (`rescore.rs:163`) | Training-time-only Percolator-style CV GBDT; andes does NOT compute production FDR | On-the-fly per run | N/A (Percolator is prod FDR) | R1 |

**Feature-count contracts (load-bearing):** Peak=18, Frag=19, Rich-ion=23 (=19+4), hard-validated at load (`gbdt_eval.rs:198`). [R1, R2]

### 1B. What andes is MISSING

| Missing capability | What it would do | Does andes have pieces? | Source |
|---|---|---|---|
| **RT-prediction model** | Predict backbone+glycan RT → ΔRT rescoring feature | **NO model.** `predicted_rt` is an always-null placeholder (`qpx.rs:472` appends null unconditionally). Zero predictor code repo-wide. Has GBDT engine + `rt_seconds` carried | R1(d) [strong, multi-grep], R2 |
| **Glycan-Y intensity model** (learned) | Learn Y0/Y1/+HexNAc/+Hex ladder intensities per monosaccharide/rung/charge | **NO.** `loss_class` column scaffolded in store schema (0/1/2) but trainers never emit glyco rows | R2-D, R4 #6 |
| **Oxonium intensity/pattern model** | Learn oxonium intensities vs fixed threshold | NO | R2-D, R4 #2 |
| **Intensity model in SELECTOR** | Use b2/b3 to pick backbone, not just rescore winner | Pieces EXIST (`ctx.intensity_model`, `frag_intensity_model`, `rich_ion_model` all on-hand at selection site) — only wiring missing | R3 §4 |
| **Separate learned glycan-vs-peptide fused score** | a glyco search engine-style RankSVM fusing ScoreG+ScoreP | NO learned fusion; Y-ladder+rank are fixed-order tiebreak | R4, R6 §1 |
| **Learning-to-rank selector** (RankSVM/LambdaMART) | Listwise/pairwise per-spectrum ranking of competing backbones | NO — additive fixed-weight only. Has pure-Rust GBDT engine (LambdaMART-capable) | R6 §4-5 |
| **Super-linear coverage term** (hyperscore-style `n_b!·n_y!`) | Reward CONTIGUOUS backbone coverage super-linearly | NO — RawScore is purely additive/independent | R6 §1,4 |
| **2D / multi-dim FDR** (peptide×glycan) | Joint + independent glycan/peptide FDR | NO (stays Percolator-boundary; glycan dim = separate thin decoy) | R4 #9 |
| **Glycan-level decoy** | Random-mass / Y-shift glycan decoy for glycan-axis FDR | Partial: `glycan_y_intensity_decoy` (a glyco search engine-style Y-shift) exists heuristically | R4 #1, R2-D |
| **Two-file glyco/non-glyco FDR separation** | O-Pair-style separate q-value pools | NO | R4 #10 |
| **Shared-backbone pooling across glycoform ladder** | Borrow backbone evidence from co-eluting siblings | NO | R4 #3 |
| **Multi-site glyco localization DAG** | O-Pair graph over S/T/N with glycan-accumulation edges | NO | R4 #12 |
| **Activation-aware ETD/EThcD c/z ions** | Glycan-retaining ions for localization | NO — HCD-centric b/y | R4 #15 |
| **De novo glycan sequencing fall-through** | Sequence glycan tree when DB fails | NO | R4 #14 |
| **Entrapment-glycome validation** (glycan axis) | Foreign-glycome cross-hit FDP check | NO for glycan axis | R4 #11 |

---

## PART 2 — Field-vs-andes Gap Matrix

| Foundational capability | Field does | andes does | Delta | Pieces on hand? |
|---|---|---|---|---|
| **Separate glycan/peptide scores** | a glyco search engine 2.0: distinct ScoreG+ScoreP, **learned RankSVM fusion** (α/β/γ/w) [R6] | Y-ladder + RawScore as **fixed-order tiebreak** (`collapse_cmp`) [R3] | No *learned* fusion of the two axes | Both sub-scores computed; fusion weights not learned |
| **Naked-backbone intensity in selector** | the reference glyco engine: hyperscore over glyco-aware ion set (coverage×intensity) at selection [R6] | Selector = rank-LLR ONLY; intensity model runs winner-only [R3 §3] | Intensity discrimination absent at selection | **YES** — models on-hand at `glyco_search.rs:782-828`; caveat: b2 EXTRAPOLATES on glycan-sized mod-delta, over-predicts backbone b/y for glyco [R2] |
| **Y-ladder-primary ranking** | a glyco search engine glycan-first indexing; StrucGP Y-ladder modules [R4] | Y-ladder IS primary (`y_primary=true`) [R3] | **andes MATCHES field here** — but Y-ladder is heuristic, not learned; note R7: Y-ladder never over-ranks truth, so it's not the scoring failure | Y-ladder present; learned version missing |
| **Peptide-conditioned intensity** | Prosit/DeepLC/AlphaPeptDeep; glyco: DeepGP/Deea glyco search engine (GNN/tree-LSTM) [R5] | Frag-intensity GBDT (b2) — but backbone-only, survivor-biased, glycan-blind [R2] | No glycan-structure encoder; b2 extrapolates past training support | Backbone GBDT yes; glycan encoder no |
| **2D FDR (peptide×glycan)** | a glyco search engine, GproDIA, GlycReSoft SVM [R4] | 1D PSM FDR via Percolator only [R4] | Glycan-axis FDR absent | Percolator boundary must hold; glycan dim = separate thin layer |
| **RT-prediction rescoring** | DeepLC/Prosit-iRT; glyco: DeepGP + Klein relative-RT (**+9.7% alone, +15.7% combined**, resolves isobaric glycoforms) [R5] | **NONE** — `predicted_rt` always null [R1] | Entire orthogonal axis missing | GBDT engine + rt_seconds present; per-monosaccharide model is ~5 coefficients [R5] |
| **Learning-to-rank selection** | a glyco search engine RankSVM; LambdaMART (top proteomics rank-learner); Casanovo-DB learned score (+31–102% pre-FDR) [R6] | Fixed-weight additive argmax [R6] | Selector sets the ceiling Percolator can't lift | Pure-Rust GBDT engine (LambdaMART-capable); truth.tsv corpus exists |

**Cross-cutting structural insight [R6 §3]:** Percolator only re-ranks candidates the SELECTOR surfaced; it cannot resurrect a demoted true backbone. The selector — not the rescorer — sets the hard ceiling. This is exactly the 80%-gen / 17%-peptide signature.

---

## PART 3 — Why andes Misses What the reference engine Gets (data-grounded, R7)

andes reproduces only **22.2%** (116/523) of the reference engine's glyco IDs. The loss breakdown:

| Failure mode | Count | % | Mechanism |
|---|---|---|---|
| **`truth_absent` — GENERATION loss** | 301 | **57.6%** | Correct backbone never enumerated |
| `truth_outranked` — SCORING loss | 106 | 20.3% | Correct candidate present, not top-1 |
| `top1_correct` — andes wins | 116 | 22.2% | — |
| `no_candidates` | 0 | 0% | — |

**The loss is dominated by GENERATION, ~3:1 over scoring.** Fixing the scorer alone caps recovery at ~20 of 78 missing points.

**The systematic bias is BACKBONE MASS (+ charge/precursor-mass correlates), NOT the glycan:**
- Backbone mass: WON median **1348 Da** vs MISSED **1753 Da**. Win rate: bb<1000→68%, 1800–2200→11%, **>2200 Da→0% (84% absent)**.
- Charge: z2→56%, z3→23%, z4→13%, **z5→0% win (100% absent)**.
- **Glycan mass is identical** in won vs missed (median 2205 vs 2206) — glycan is NOT the discriminator.
- No glycan-class blind spot (sialylated dominates both won 83% / outranked 93%).

**Generation-loss root cause (mass mis-partitioning):** In 89% of absent scans, andes' winner picks a backbone **>50 Da smaller** than truth (median −688 Da) and **compensates with an oversized glycan** (2409 vs 2206). **andes assigns too much precursor mass to glycan, too little to peptide** → enumerates short decoy backbones + oversized glycans instead of the true long backbone.

**Scoring-loss root cause:** 0/106 outranked cases have truth with a higher Y-ladder than the winner. Losses come through the **rank/RawScore tiebreaker** (rank_gap median 7 in winner's favor); 53.8% `yladder_tie_loses_rank` + 46.2% `truth_loses_y_ladder`. **The scoring lever is the rank-score backbone-fragment term, NOT the Y-ladder.**

**Bottom line:** andes misses **large-backbone, high-charge (z≥4), high-precursor-mass N-glycopeptides**, primarily at candidate GENERATION via precursor-mass mis-partitioning — a backbone/mass problem, not a glycan-class or Y-ion problem.

---

## PART 4 — Refined Roadmap (ranked by leverage/cost)

The single most important correction from R7: **the old Phase-1 premise (wire intensity/y0y1 into the selector) attacks the 20% scoring loss while the 58% generation loss goes unaddressed.** Generation must lead.

### Tier 0 — GENERATION (the 58% loss; nothing else can recover it) — **NEW, highest leverage**

**P0. Fix precursor-mass partitioning for large-backbone / high-charge glycopeptides.** [FOUNDATION]
- *Why:* Directly attacks the 301/523 generation losses. R7 proves andes systematically under-sizes the backbone and over-sizes the glycan; z≥4 and bb>2200 Da are ~0% win. No downstream scoring/RT/LTR work can recover a backbone never enumerated [R6 §3, R7].
- *Cost:* Medium — candidate-generation logic, not ML. Investigate charge-state handling (z5 = 100% absent suggests a charge ceiling/blind spot — cf. memory "charge-1-only blind spot"), backbone mass-window upper bound, and glycan/backbone mass-split enumeration order.
- *Pieces on hand:* Yes — generation cascade exists; this is a bound/partition fix.
- *Metric:* clean-truth backbone find-rate (NOT 1%-FDR count) [R6 guard].
- *Flag:* R7 is single-source (one 523-scan the reference engine-truth analysis) — validate the mass-partitioning mechanism on a second dataset before large refactor.

### Tier 1 — SELECTOR (the 20% scoring loss; sets the Percolator ceiling)

**P1a. Add super-linear contiguity/coverage term to the selector.** [FOUNDATION, additive-safe]
- *Why:* R6 §1,4 — RawScore is purely additive; treats 3 scattered ≈ 3 contiguous matches. R7 shows scoring loss runs through the rank-score tiebreak. A hyperscore-style `n_b!·n_y!` contiguity term over glyco-offset-corrected b/y directly strengthens that tiebreak. Purely additive → low regression risk (matches "additive/top-1-preserving only" lesson).
- *Cost:* Low. *Pieces:* Yes.

**P1b. Wire assignment-aware intensity LLRs into the selector (not the cosine).** [FOUNDATION]
- *Why:* R3 §3-4 — selector uses rank-LLR only; b2/b3 run winner-only. Inject `frag_llr_battery`/`rich_ion_llr` (NOT `intensity_signal` cosine, which R3 notes gives no lift) into the Phase-1 loop.
- *Cost:* Medium — main cost is GBDT eval × candidate fan-out; mitigate by gating to top-K-by-rank shortlist (re-rank, not full re-score) [R3 §4].
- *Caveat (R2):* b2 EXTRAPOLATES on glycan-sized mod-delta and over-predicts backbone b/y for glyco spectra — validate it actually helps glyco selection before shipping; may need glyco-aware retraining first.
- *Pieces:* Yes — all models on-hand at selection site.

### Tier 2 — ORTHOGONAL RESCORING AXES (compose with Percolator)

**P2. RT-prediction rescoring feature (per-monosaccharide relative-RT).** [FOUNDATION, promoted UP from unlisted] — **best cost/leverage new lever**
- *Why:* R5 is emphatic — for glyco, RT is *disproportionately* valuable because the dominant error is glycan-composition ambiguity among co-eluting isobaric glycoforms, which MS2 can't resolve but RT can. Klein 2024: **RT alone +9.7% high-confidence, −44% low-confidence; resolves quasi-isobaric misassignments (56% of one class reassigned).** The glycan RT shift is **additive per monosaccharide** (Klein/Zaia 2019: sialic +1.94, fucose/high-mannose −%ACN).
- *Cost:* **LOW** — the whole trick is ~5 fitted per-monosaccharide coefficients on top of a lightweight backbone-RT GBDT, **self-calibrated per run** from confident targets (own-data, no external model, no patent — clean on all 4 andes objectives). Emit `abs_delta_rt`, `signed_delta_rt`, `dRT_rank_among_glycoforms` (the last kills the isobaric-glycoform error mode).
- *Pieces:* Yes — GBDT engine + `rt_seconds` carried; `predicted_rt` schema field already reserved (just needs filling, R1).
- *Note:* Fixes the always-null `predicted_rt` gap [R1(d)].

### Tier 3 — LEARNED FUSION & 2D-FDR (structural, higher cost)

**P3a. a glyco search engine-style RankSVM learned fusion of backbone-LLR + glyco-evidence.** [FOUNDATION]
- *Why:* R6 §5(B), R4 — the one proven-for-glyco learned selector. Replaces the fixed Y-ladder/rank tiebreak with trained weights on the truth.tsv corpus. Linear, interpretable, fast.
- *Cost:* Medium. *Pieces:* Yes — sub-scores computed; truth corpus exists.

**P3b. LambdaMART listwise selector on the existing GBDT engine.** [FOUNDATION, higher ceiling]
- *Why:* R6 §5(C) — per-spectrum listwise objective (query=spectrum, label=truth) optimizing rank-of-true; LambdaMART was top proteomics rank-learner. Surfaces the correct backbone before top-k truncation.
- *Cost:* Medium-high. *Pieces:* Yes — pure-Rust GBDT engine (`feat/gbdt-stronger-models`).

**P3c. Glycan-Y intensity model + 2D/glycan-level FDR.** [FOUNDATION]
- *Why:* R2, R4 #1,#6,#9 — learn Y-ladder intensities (`loss_class` scaffold already in store schema); add glycan-level decoy (partial `glycan_y_intensity_decoy` exists) + thin 2D-FDR layer (stays Percolator-boundary; glycan dim = separate estimator).
- *Cost:* High (new trainer path + benchmarking). *Pieces:* Partial — schema scaffold + heuristic Y-decoy.

### Tier 4 — RESEARCH-TIER (defer)

- **Glyco-aware GNN/tree-LSTM spectrum prediction** (DeepGP/Deea glyco search engine-style) [R4 #6, R5] — highest ceiling for isomer discrimination but heaviest, hardest to keep own-data/patent-clean. [TRICK→FOUNDATION long-term]
- **Shared-backbone pooling across glycoform ladder** (a cross-spectrum glyco engine) [R4 #3] — cheap non-ML win to rescue weak-backbone spectra; worth a spike after Tier 0-1.
- **ETD/EThcD c/z ion support** [R4 #15], **multi-site DAG localization** [R4 #12], **de novo glycan fall-through** [R4 #14] — acquisition/capability expansions, out of current scope.

### Ranked summary (leverage ÷ cost)

1. **P0 — precursor-mass partitioning fix** (58% of loss; medium cost) — *do first, nothing else recovers generation loss*
2. **P2 — per-monosaccharide RT rescoring** (proven +9.7% glyco, ~5 coefficients, own-data) — *best new cost/leverage*
3. **P1a — additive contiguity term** (low cost, additive-safe)
4. **P1b — assignment-aware intensity LLR in selector** (medium; needs glyco-validation of b2 extrapolation)
5. **P3a — RankSVM fusion** (proven-for-glyco learned selector)
6. **P3b — LambdaMART listwise selector** (higher ceiling)
7. **P3c — glycan-Y intensity + 2D-FDR** (structural, high cost)

**Single-source / confidence flags:**
- R7's generation-vs-scoring split and mass-partitioning mechanism = **one 523-scan analysis** — reconfirm on a second dataset before the P0 refactor.
- R5's glyco RT gains (+9.7%/+15.7%) = **Klein 2024 (single group), a glyco search engine-relative** — directionally strong and mechanistically sound (additive per-monosaccharide is independently supported by Klein/Zaia 2019), but the exact % won't transfer to andes; treat as "high-value, validate magnitude."
- R2's "b2 over-predicts backbone b/y for glyco" = mechanistic inference from training-support extrapolation, not a measured andes benchmark — verify empirically before relying on b2 in the glyco selector (P1b).
- R6's Casanovo-DB pre-FDR gains (+31–102%) are for **non-glyco tryptic** search — cited as motivation for learned selectors, not a glyco promise.
---

## ADDENDUM (2026-07-07) — P0/P0b find-rate A/B RESOLVES generation-vs-scoring

Two conclusions from earlier need correction, now data-grounded on the 523 truth scans (all-hits, driver charge):

**1. P0 charge-expand (ANDES_GLYCO_CHARGE_EXPAND) REFUTED.** truth_absent 301→309 (worse). Charge under-calling was the wrong hypothesis; expanding charges just adds competing candidates.

**2. The "58% generation loss" is ~half ACCEPTANCE-TRUNCATION, and retention DOESN'T CONVERT — scoring is the true binding constraint.** Re-diagnosis (base truth_absent=301, top1_correct=115):
| config | truth_absent | top1_correct | outranked |
|---|---|---|---|
| base | 301 | 115 | 107 |
| YINDEX=1 | 239 | 119 | 165 |
| top-k 500 | 200 | 115 | 208 |
| YINDEX+top-k500 | **159** | **114** | 250 |
YINDEX+top-k recovers 142 truth backbones into the pool (301→159) but top1_correct stays FLAT (~115) — all recovered land in truth_OUTRANKED. So: (a) ~142 = acceptance-truncation (b/y-rank top-k gate drops large/weak-b/y true backbones), recoverable by Y-aware retention; (b) retention alone is USELESS + would hurt @1% (outranked noise, cf. isolation exp); (c) top1_correct is SCORER-capped at ~115 regardless of retention → **SCORING is the ultimate binding lever** (vindicates the foundational L1/L2: peptide-conditioned scoring in the selector, which the generation-first correction had demoted); (d) ~159 = true generation-loss residual (deeper charge/mass/DB, secondary).

**REVISED PRIORITY:** the SELECTOR-scoring fixes (P1a contiguity + P1b intensity-in-selector + P3a a glyco search engine separate glycan/peptide scores) are THE lever, PAIRED with Y-aware retention (default AXIS-2 YINDEX on) so the true candidate survives to be scored. Validate on top1_correct + @1% TOGETHER (retention without scoring hurts @1%). P0 generation residual (~159) is secondary. P2 RT (orthogonal rescoring) still valuable.
