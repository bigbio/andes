# Glyco Scoring — Foundational Roadmap (to beat the reference engine/a glyco search engine)

**Status:** ROADMAP for review — no implementation committed. Produced 2026-07-06
from a 5-agent consensus study (a glyco search engine2.0 + the reference glyco engine + fragmentation
physics + andes code audit + SOTA intensity models) **independently confirmed**
by a Codex adversarial review of the scoring code. Sources triangulate.

## 0. The question this answers  — CORRECTED 2026-07-07 by the deep-review

andes beats an open-source glyco engine (253/97 @1% FDR, deterministic) but reproduces only
**22.2% (116/523)** of the reference engine's curated truth. **An 8-agent deep review with a
data-grounded gap analysis CORRECTED the earlier "scoring is the binding
constraint" conclusion.** The true breakdown of the 523 truth scans:

| Failure mode | Count | % | Fix tier |
|---|---|---|---|
| **Generation loss** — true backbone never enumerated/emitted | **301** | **57.6%** | **Tier 0 (do first)** |
| Scoring loss — present but out-ranked | 106 | 20.3% | Tier 1 |
| andes wins | 116 | 22.2% | — |

**Generation dominates ~3:1.** Root cause (R7, data-grounded): **precursor-mass
mis-partitioning** — andes systematically assigns too much mass to the glycan and
too little to the peptide. In 89% of missed scans its winner's backbone is >50 Da
*smaller* than truth (median −688 Da) with an oversized glycan. The bias is
**backbone mass + charge, not glycan**: backbone >2200 Da → **0% win**; charge z5
→ **0% win (100% absent)** = a charge blind spot; glycan mass is *identical* in
won vs missed. This also explains why enlarging the glycan list *lowered* @1% (the
isolation experiment) — more glycans = more wrong short-backbone/big-glycan
alternatives. **Fixing the scorer alone caps recovery at ~20 of the 78 missing
points; generation must lead.**

Full model inventory, field-vs-andes matrix, and the gap data:
`50-roadmap/deep-review-synthesis.md`. §1–§3 below (physics, foundations, scoring
gaps) remain valid for the Tier-1 scoring work; the roadmap in §4 is superseded
by the corrected tiers in §4b.

## 1. The physics (why it's hard) — the foundation everything rests on

An intact N-glycopeptide is two molecules with very different bond strengths.
Glycosidic bonds dissociate at far lower energy than peptide amides, so under
HCD the glycan comes apart first and acts as an energy sink. Three orthogonal
ion populations result:

| Ion | Physically encodes | Correct scoring role |
|---|---|---|
| **Oxonium/B** (204, 366, 274/292) | glycan **class only** (detached sugar cations) | **GATE, never a discriminative score** — composition-degenerate |
| **Y-ions** (Y0=bare peptide, Y1..core ladder) | peptide **MASS** + glycan composition (not sequence) | primary, abundant **backbone-MASS** discriminator |
| **b/y** (c/z under ETD) | backbone **SEQUENCE** — *which* peptide | scarce, information-rich, physically **suppressed** in HCD |

**Consequence (the hard truth):** glycan mass and peptide *mass* are
well-determined from one spectrum; peptide backbone *identity* is the
bottleneck and is often mass-degenerate. The correct backbone is pinned in two
stages — (i) **Y-ions rank the backbone mass** and kill most wrong candidates
cheaply; (ii) **b/y sequence ions (or their intensity pattern) resolve which
peptide** among survivors. Bare b/y match-count alone is the *noisiest* signal
on exactly these spectra.

## 2. The five foundations SOTA shares (a glyco search engine2/3, the reference engine, Deea glyco search engine)

Distilled from the consensus; these are **foundations, not tricks**:

- **F1 — Two independent scores on orthogonal axes** (glycan-Y score + peptide-b/y
  score), fused with a learned weight. Right-peptide/wrong-glycan is a *distinct*
  error from wrong-peptide/right-glycan; only two axes separate them.
- **F2 — Glycan = precursor mass-offset; peptide scored on the NAKED backbone.**
  Aligns the theoretical ion set with what HCD physically produces (the reference engine:
  21 matched ions / >50% ion current vs 8 for a variable-mod treatment).
- **F3 — Y-ladder (trimannosyl core) as the primary backbone-MASS discriminator
  and a cheap high-precision gate** (a glyco search engine's ≥2 core-Y filter *before* peptide
  scoring).
- **F4 — Peptide-conditioned intensity prediction** turns match-count (~1 bit/peak)
  into graded spectral similarity and makes *expected-but-absent* peaks
  informative — the only lever that raises separation without new ions.
- **F5 — 2D target-decoy FDR** on separate glycan and peptide axes
  (`FDR = FDR_G + FDR_P − FDR_{G∩P}`). A glycoPSM is false if *either* moiety is
  wrong; 1D FDR can't tell which.

**Tricks (real but incidental, not foundations):** the exact hyperscore factorial
algebra; fragment-ion indexing (speed, not discrimination); the specific fusion
constants/optimizer; specific NN architectures.

## 3. What andes structurally lacks (code-grounded; consensus + Codex agree)

- **The default selector ranks backbones by bare b/y rank-LLR** (`rk = score_psm +
  psm_edge_score`, `glyco_search.rs:801-803`); the glycan contributes **nothing**
  to the rank — it only sets the mass window. This is the physically-suppressed
  signal (§1) used as the *primary* discriminator. [Codex HIGH; consensus Gap-1]
- **The reliable Y-evidence is only a tiebreaker.** AXIS-1 truncation
  (`glyco_search.rs:845-859`) sorts by the noisy b/y rank with `core_y_hits` as
  tiebreak; the AXIS-2 Y-rescue is **off by default** (`ANDES_GLYCO_YINDEX`). A
  correct weak-b/y backbone is dropped before the Y-primary collapse ever sees
  it. Truncation is b/y-primary but the final collapse is Y-primary — an internal
  inconsistency. [Codex HIGH; consensus Gap-2]
- **The one peptide-conditioned glyco feature is inert.** `y0y1_anchor_score`
  (Y0/Y1 depend on peptide mass) is computed **after** the winner is collapsed and
  emitted only as a PIN column — it never enters the ranker or `collapse_cmp`, so
  it cannot rescue the correct backbone, and Percolator never sees the losing
  alternatives. Every other glyco feature is backbone/spectrum-level
  (non-discriminative for *which* peptide). [Codex HIGH; consensus Gap-3]
- **No 2D FDR.** The only glyco decoy is the glycan-axis shifted-composition decoy;
  there is no reversed-peptide decoy scored on the glyco fragment set, and
  everything is one 1D Percolator pile (`glyco_pin.rs:426-434`). [Codex MEDIUM;
  consensus Gap-4/F5]

**Net:** a selection-architecture failure. The true backbone is generated but
ranked by noise, while the informative signals are tiebreakers, off, or inert.

## 4. The roadmap (phased; each phase gated by validation)

Ordered by leverage-per-cost. Foundations first; tricks as cheap companions.
**Every phase is validated on the generated-but-outranked truth scans** (the
`glyco_outrank_audit.py` harness) — not just total @1%, to avoid the
candidate-flood trap (isolation lesson) and FDP-blind gains (refine-no-op lesson).

### Phase 1 — Fix the selector to use the signals andes ALREADY computes (FOUNDATION, low cost)
The highest-leverage, cheapest change; both sources rank it first.
- **L1 (foundation):** move the **peptide-conditioned** signal into the scalar that
  drives *truncation and collapse* — combine `y0y1_anchor_score` (and a
  normalized peptide-conditioned Y-ladder term) with b/y into `rk`/`collapse_cmp`,
  instead of ranking by bare b/y and emitting the anchor as an inert column.
- **L4 (high-value trick):** make AXIS-1 truncation Y-aware (rank retained
  backbones on a Y-inclusive scalar) and/or default the AXIS-2 Y-rescue on, so
  weak-b/y correct backbones survive to the selector.
- **Gate:** truth-outranked recovery ↑ AND total @1% not worse AND decoys@1%
  controlled. Add Codex's regression: *assert generated true backbones survive
  pre-feature retention.*
- **Cost:** medium (wiring exists; signals already computed). **Do this first.**

### Phase 2 — Peptide-conditioned intensity similarity (FOUNDATION) — LARGELY ALREADY BUILT
**Key finding (2026-07-07, code-verified):** andes already has the intensity model
and it already scores glyco PSMs. `--score strong` = a GBDT fragment-intensity
model trained on PRIDE data (high-res; Astral +23%). In the glyco path,
`compute_psm_features(..., ctx.intensity_model)` (`glyco_search.rs:976-982`) already
computes the intensity-derived features (`IntensitySignal`, `FragPredExplained`,
`FragPredChanceLLR`, `FragTopKObserved`) and emits them as glyco PIN columns — **but
only for the winner already chosen by rank-LLR.** So the "intensity foundation" is
NOT a model-building project; it is the SAME wiring fix as L1.
- **L2a (do WITH L1, low cost):** use the existing `intensity_model` /
  strong-score peptide-b/y similarity as (part of) the backbone **selector**
  scalar, not just a post-collapse PIN column. The model predicts peptide b/y
  intensities, which transfer to the NAKED glyco backbone (F2). High-res only
  (Fc3 = Q-Exactive HF ✓); do NOT use on low-res ion-trap glyco.
- **L2b (later, medium/high):** the model does NOT cover glycan-**Y**-ion
  intensities. Add a glycan-Y intensity prior (rule-based core-Y ladder first;
  learned tree-LSTM/GNN glycan branch only if L2a proves out) to get the glycan
  axis (F1's ScoreG analogue).
- **Gate:** additive/selector changes measured on outranked scans; parity-safe;
  entrapment/2D-validate (Phase 3) before trusting.

### Phase 3 — Two-axis scoring + 2D FDR (FOUNDATION + the validation mechanism)
- **L3:** expose two separable axes (peptide-similarity vs glycan-similarity) and a
  **peptide-axis glyco decoy** (reversed peptide scored on the glyco fragment set)
  alongside the existing glycan-axis decoy. Combine as a **thin post-process of
  Percolator** (peptide-q × glycan-q with the union rule) — never an
  andes-internal FDR (respects the Percolator-only boundary).
- **Why here:** L3 is *also the validation lever* for L1/L2 — it distinguishes
  right-peptide/wrong-glycan from wrong-peptide/right-glycan and catches
  FDP-blind gains. **Cost:** medium; depends on L2's separable axes.

### Companion tricks (cheap, low-risk, anytime)
- **L5:** additive combinatorial-coverage (factorial-style) PIN feature — rich
  Y/b-y match sets give super-linear separation vs additive rank-LLR. Additive,
  not a `score_psm` rewrite (parity lesson).
- **L6:** enforce oxonium strictly as a gate/QC, never a discriminative score
  contributor (all studies agree it's composition-degenerate).

## 4b. CORRECTED ROADMAP (deep-review, generation-first; supersedes §4)

Ranked by leverage ÷ cost. "Pieces on hand?" = does andes already have the code.

### Tier 0 — GENERATION (the 58% loss; nothing downstream can recover it)
**P0 — Fix precursor-mass partitioning for large-backbone / high-charge glycopeptides.**
[FOUNDATION, medium cost, NOT ML]
- Attacks the 301/523 generation losses. andes under-sizes the backbone, over-sizes
  the glycan; z≥4 and bb>2200 Da are ~0% win. Investigate: the **charge-state
  ceiling** (z5 = 100% absent — cf. memory "charge-1-only blind spot"), the
  **backbone mass-window upper bound**, and the **glycan/backbone mass-split
  enumeration order** (stop enumerating short-backbone + oversized-glycan decoys
  ahead of the true long backbone). Pieces: yes (generation cascade exists).
- **Metric:** clean-truth backbone find-rate (NOT 1%-FDR count).
- **Flag:** R7 is one 523-scan analysis — reconfirm the mass-partitioning mechanism
  on a 2nd dataset before a large refactor.

### Tier 1 — SELECTOR (the 20% scoring loss; sets the Percolator ceiling)
- **P1a — additive super-linear contiguity/coverage term** (hyperscore-style
  `n_b!·n_y!` over glyco-offset-corrected b/y). [FOUNDATION, low cost, additive-safe]
  RawScore is purely additive (3 scattered ≈ 3 contiguous); scoring loss runs
  through the rank-score tiebreak, not the Y-ladder (R7: 0/106 lose on Y-ladder).
- **P1b — wire assignment-aware intensity LLRs into the selector** (the
  `frag_llr`/`rich_ion_llr`, NOT the `intensity_signal` cosine which gives no lift).
  [FOUNDATION, medium] Gate to a top-K-by-rank shortlist (re-rank, not full
  re-score). **Caveat:** the frag-intensity GBDT *extrapolates* on glycan-sized
  mod-deltas and may over-predict backbone b/y for glyco — validate it helps glyco
  selection (may need glyco-aware retraining) before shipping. Pieces: yes.

### Tier 2 — ORTHOGONAL RESCORING (composes with Percolator)
**P2 — RT-prediction rescoring (per-monosaccharide relative-RT).** [FOUNDATION,
**LOW cost — best new cost/leverage**]
- For glyco, RT is disproportionately valuable: the dominant residual error is
  ambiguity among co-eluting isobaric glycoforms, which MS2 can't resolve but RT
  can. Published: RT alone +9.7% high-confidence, resolves isobaric misassignments;
  glycan RT shift is **additive per monosaccharide** (sialic +, fucose/high-mannose −).
- The whole thing is ~5 fitted per-monosaccharide coefficients on a lightweight
  backbone-RT model, **self-calibrated per run** from confident targets (own-data,
  no external model, no patent — clean on all 4 andes objectives). Emit
  `abs_delta_rt`, `signed_delta_rt`, `dRT_rank_among_glycoforms` as PIN features.
- **Fills the always-null `predicted_rt` gap** (`qpx.rs:472`). Pieces: yes (GBDT
  engine + `rt_seconds` carried; schema field reserved).
- **Flag:** the +9.7%/+15.7% figures are single-group (Klein 2024); mechanism is
  sound, exact magnitude for andes must be measured.
- **Concrete engine-wide design (SOTA-reviewed):** see
  `50-roadmap/rt-prediction-design.md` — GBDT RT-index (reuse `GbdtPeakModel`) +
  per-run linear calibration + `DeltaRT`/`AbsDeltaRT`/`DeltaRTNorm` in
  `psm_feature_values` (Commit 1, benefits regular search); per-monosaccharide
  glycan offset (NeuAc own sign) + `DeltaRTRank`-among-glycoforms (Commit 2).
  GBDT≈SSRCalc tier (R≈0.96–0.98), sufficient as a rescoring feature.

### Tier 3 — LEARNED FUSION & 2D-FDR (structural, higher cost)
- **P3a — a glyco search engine-style RankSVM learned fusion** of backbone-LLR + glyco-evidence
  (replaces the fixed Y-ladder/rank tiebreak; trained on truth.tsv). [FOUNDATION, medium]
- **P3b — LambdaMART listwise selector** on the existing pure-Rust GBDT engine
  (query=spectrum, label=truth; surfaces the true backbone before top-k truncation).
  [FOUNDATION, medium-high]
- **P3c — glycan-Y intensity model + glycan-level decoy + thin 2D-FDR** (the
  `loss_class` store scaffold + heuristic `glycan_y_intensity_decoy` already exist;
  glycan-axis FDR stays a thin post-process, Percolator-boundary intact). [FOUNDATION, high]

### Tier 4 — RESEARCH (defer)
Glyco-aware GNN/tree-LSTM spectrum prediction (DeepGP/Deea glyco search engine) · shared-backbone
pooling across the glycoform ladder (a cross-spectrum glyco engine — cheap non-ML spike worth
trying after Tier 0-1) · ETD/EThcD c/z ions · multi-site DAG localization · de-novo
glycan fall-through.

### Ranked summary
1. **P0** precursor-mass partitioning (58% of loss) — *do first*
2. **P2** per-monosaccharide RT rescoring (~5 coeffs, own-data, proven +9.7%) — *best new lever*
3. **P1a** additive contiguity term (low cost, additive-safe)
4. **P1b** intensity LLR in selector (needs glyco-validation)
5. **P3a** RankSVM fusion → **P3b** LambdaMART → **P3c** glycan-Y intensity + 2D-FDR

## 4c. GOVERNING PRINCIPLE — model changes are ENGINE-WIDE, not glyco-only (user, 2026-07-07)
Any change to a learned model or the GBDT engine must be a **general andes
capability**, not a glyco-only bolt-on. Concretely:
- **RT prediction (P2)** is a general model: a backbone/peptide RT predictor that
  benefits *regular* peptide search too (a general ΔRT rescoring feature), with the
  per-monosaccharide glycan offset as a glyco-specific *addition* on top. Train,
  store (parquet ModelStore), and expose it engine-wide; glyco consumes it.
- **Any GBDT change** (new features, retraining, LambdaMART, glyco-aware intensity)
  flows through the shared `model-train` + `scoring` crates and the model store, so
  the rank/intensity models improve for all search modes, not just glyco.
- Rationale: avoids a divergent glyco-only model fork, keeps one training/eval
  path, and lifts the whole engine. Glyco-specific pieces (glycan-Y intensity,
  per-monosaccharide RT offset) are *extensions* of the general models.

## 5. Cross-cutting discipline (do not skip)
- **Validate on outranked truth scans + 2D/entrapment FDR**, not total @1% alone.
- **High-res only** — SOTA intensity models assume Orbitrap; do not expect
  transfer to low-res ion-trap glyco.
- **Additive PIN features only** for scoring changes; modifying existing features
  regresses Percolator (andes parity lesson).
- **Confirm single-study numbers** before hard-coding: the 50 kcal/mol charge
  thermodynamics, exact Deea glyco search engine percentages.
- **a glyco search engine 2.0 constants VERIFIED first-hand** (PMC5585273, 2026-07-07) — usable
  starting values for L1/L2:
  - `ScoreG = Σ log(intensᵢ)·(1 − merritolᵢ/4)·ratioion^0.56·ratiocore^0.42`
    (Y-ions + diagnostic glycan ions; `ratiocore` = matched/theoretical
    trimannosyl-core Y-ions — first-class factor).
  - `ScoreP = Σ log(intensᵢ)·(1 − merritolᵢ/4)·ratioion^0.94` (b/y; note:
    intensity-weighted, NOT match-count — andes's bare rank-LLR match is *below*
    even this 2017 baseline).
  - `ScoreGP = 0.35·ScoreG + 0.65·ScoreP` (fusion weight learned by Ranking-SVM).
  - Comprehensive FDR `= FDR_G + FDR_P − FDR_{G∩P}`; glycan decoy = finite-mixture
    method; ≥2 trimannosyl-core-Y pre-gate; top-300 peak filtration (vs 50).

## 6. Recommendation (CORRECTED)
The deep-review overturns the scoring-first plan: **58% of the gap is GENERATION,
not scoring.** So:
1. **P0 first — precursor-mass partitioning / charge blind spot.** No scoring, RT,
   or learned-selector work can recover a backbone that is never enumerated. Gate
   on clean-truth backbone find-rate; reconfirm the mass-partitioning mechanism on
   a 2nd dataset before a large refactor.
2. **P2 in parallel — per-monosaccharide RT rescoring.** Best new cost/leverage
   (~5 self-calibrated coefficients, own-data, fills the null `predicted_rt`),
   orthogonal to P0 so it composes cleanly.
3. **Then Tier 1 (P1a/P1b)** for the 20% scoring loss, then Tier 3 learned
   fusion/2D-FDR.

Everything remains gated behind validation on the truth scans + 2D/entrapment FDR;
no implementation is committed until you approve the sequence.
