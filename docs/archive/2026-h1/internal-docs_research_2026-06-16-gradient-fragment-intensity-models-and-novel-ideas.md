# Gradient/tree-based (and adjacent) MS/MS fragment-intensity models — survey + novel ideas for andes

Date: 2026-06-16
Author: research agent (for the andes lead engineer)
Audience: andes scoring/model-train owners. Grounded in the v3 design
(`internal-docs/specs/2026-06-15-andes-own-fragmentation-score-v3-design.md`) and the actual code:
`crates/scoring/src/frag_features.rs`, `crates/model-train/src/gbdt/train.rs`,
`crates/model-train/src/gbdt/frag_dataset.rs`, `crates/scoring/src/scoring/strong_score.rs`.

**Hard constraints repeated up front (every recommendation respects these):** OWN models; pure-Rust GBDT;
NO neural nets; NO MS2PIP/Prosit wrapping or re-implementation; the model is a RESCORING feature (top-N,
additive PIN columns; the rank-LLR `RankScore` stays byte-identical); andes already has a noise/chance-match
denominator (`strong_score.rs`), which is the structural differentiator most intensity models lack. Targets
are per-regime: HCD/Astral high-res, low-res CID, TMT.

---

## 0. TL;DR for the lead engineer

The field has converged on one idea: **predict per-fragment relative intensity from local sequence context,
then use `similarity(predicted, observed)` as a rescoring feature.** The tree-based lineage
(MS2PBPI 2014 → MS2PIP RF 2013 → MS2PIP XGBoost 2019) is exactly what andes's v3 regressor reproduces (box 1).
The NN lineage (Prosit/pDeep/DeepMass/Predfull/AlphaPeptDeep) adds (a) sequence+charge+NCE conditioning,
(b) a **spectral-angle training loss** instead of per-fragment MSE, (c) **transfer learning** across
instruments/PTMs, and (d) full-vocabulary prediction (internal ions). We can borrow (a)–(d) as *concepts*
without any NN.

The single most under-exploited lever in the whole field — and the one andes is uniquely positioned to pull —
is that **everyone optimizes intensity R²/spectral-angle, but the metric that actually matters is
target/decoy separation at 1% FDR.** andes's GBDT trainer and strong-score fusion are *the right substrate*
to optimize the discriminative objective directly. That, plus andes's existing noise denominator, is where
we beat the field rather than tie it.

Top recommendations (full roadmap in Part C): **(1)** keep the v3 regressor MVP, but **(2)** add a
listwise spectral-angle objective (not MSE) — strong literature support; **(3)** formalize the
predicted-intensity ÷ chance-match per-fragment LLR (extends our existing denominator — genuine novelty);
**(4)** add relational complement/neighbor features (cheap, novel for a tree model); **(5)** quantile /
uncertainty-weighted cosine; **(6)** isotonic/Platt FDR-calibration of the final fused feature. Defer
multi-task and conformal until the above are measured.

---

# PART A — How the field's models work

For each model: approach · features · target · training algorithm · key result · **what andes should borrow / why it may not apply.**

## A.1 MS2PBPI — Zhou, Han, Yao, et al., *Anal. Chem.* 2014, 86(15):7446-7454 (the closest ancestor of v3)
- **Approach.** Partition the matched-fragment population into "regions" by *fragmentation pathway* (ion type,
  charge, mobile-proton class, residue context) using **binary trees** that split the bulk data into tens to
  >1000 regions; fit **one stochastic-gradient-boosting regression-tree (SGBT) model per region**. Hundreds of
  small models = a hierarchical mixture-of-experts.
- **Features.** Sequence/fragmentation-pathway descriptors: residues flanking the cleavage, proton-mobility
  category, ion type/charge, position; designed to mirror the chemistry (mobile/partially-mobile/non-mobile
  proton regimes).
- **Target.** Relative fragment-ion intensity (regression).
- **Training.** Stochastic gradient boosting (Friedman) of regression trees, per region.
- **Key result.** Predicts unmodified + modified peptide spectra with good consistency across **ion-trap**
  instruments (low- and high-res), **outperforming MassAnalyzer and PeptideART**. Known weakness: degrades
  on instrument/fragmentation regimes far from training (a TOF/beam-type cliff vs the ion-trap training).
- **andes:** This *is* andes-v3 box 1, with a key refinement available to us: MS2PBPI's "regions" map onto
  andes's **per-slug models** + per-feature splits. **Borrow:** the proton-mobility region idea as an explicit
  GBDT feature (cheap; see B.4/B.6). **Why its exact design doesn't apply:** hundreds of disjoint per-region
  models fragment the (already thin) corpus; a single conditioned GBDT with mobility *as a feature* + per-slug
  specialization is the modern, data-efficient equivalent. The instrument cliff is the explicit warning behind
  recommendation B.8 (base+residual transfer over many thin models).
  Source: https://pubs.acs.org/doi/10.1021/ac501094m · https://pubmed.ncbi.nlm.nih.gov/25032905/

## A.2 MS2PIP — the primary comparator (CompOmics)
Two eras; both are pure tree-based, so this is the closest production analogue to andes-v3.

**Era 1 — Random Forest (Degroeve & Martens, *Bioinformatics* 2013, 29(24):3199-3203; web server 2015).**
- **Approach.** Per-(ion-type, charge, length-bucket) **random-forest regression**: a separate model
  `D_clf` per fragment ion type `f`, partitioned further by charge (+2/+3) and peptide length (8–28).
- **Features (the canonical feature set worth copying verbatim).** Positional one-hots
  `seq_<pos>_<amino>` and modification `seq_<pos>_<mod>`; per-position **chemical properties**
  `seq_<pos>_<chem>` = hydrophobicity, basicity, helicity, pI; **global averages** `avg_<chem>`;
  mass features `pep_mz`, `ion_mz`, `ion_mz_other` (the *complement* ion's m/z — note: MS2PIP already passes
  the sibling ion's mass as a feature, a weak form of relational context); **composition counts** `I_<amino>`.
- **Target.** Peak intensities normalized to the spectrum's total ion current, **log2-transformed**.
  (andes uses `ln(obs/base_peak)` — same family, base-peak vs TIC normalization.)
- **Training data.** ~3.97M Orbitrap PSMs → merged to ~62k (+2) / ~11k (+3) representative peptides.
- **Key result.** Higher Pearson r than PeptideART across datasets; fast, lightweight.

**Era 2 — Gradient Tree Boosting / XGBoost (Gabriels, Martens, Degroeve, *NAR* 2019, 47(W1):W295; "MS²PIPc").**
- **Change.** Random forest → **XGBoost gradient-boosted trees**, more + more-diverse training data,
  per-fragmentation-method (CID/HCD) and per-instrument/label models (TripleTOF, Orbitrap-LTQ, Q-Exactive;
  later TMT, iTRAQ, phospho, immonium). Predicts **b, y, b²⁺, y²⁺**.
- **Key result.** Big accuracy gain especially for higher charges (the +2/+3 gap shrinks); the per-ion-type +
  per-method model matrix is the deployment pattern.
- **Rescoring consumption (this is the integration andes mirrors).** **MS²Rescore** (Declercq et al.,
  *MCP* 2022/2024) feeds MS2PIP-predicted vs observed spectra into Percolator as similarity features computed
  **separately for b-ions, y-ions, or both** — **Pearson correlation, cosine similarity, spectral angle, and
  related metrics** — alongside DeepLC ΔRT features. Reported **+46% PSMs / +36% peptides at 1% FDR** vs plain
  Percolator on immunopeptidomics. **MSBooster** (Yang et al., *Nat. Commun.* 2023) does the same against
  MSFragger output; its default features are **`unweighted_spectral_entropy`** (entropy-based MS2 similarity),
  **`delta_RT_loess`**, and **`pred_RT_real_units`**, plus a top-20 predicted/observed peak **intersection**
  count — combined with MSFragger's `hyperscore` for Percolator.
- **andes:** This is the blueprint. **Borrow now:** (i) the per-position **chemical-property** features
  (basicity/hydrophobicity/helicity/pI) — andes currently passes only residue *index*, which forces the trees
  to relearn chemistry; adding 4 cheap per-flank chem scalars is the single highest-ROI feature add (B.6).
  (ii) The MS²Rescore practice of emitting **separate b-only / y-only similarity** features, not just a pooled
  cosine. (iii) MSBooster's **spectral-entropy similarity** as a second additive PIN column alongside cosine
  (different geometry, cheap, empirically strong). **Why the rest doesn't apply:** we will *not* run MS2PIP
  (independence). The "intensity-only, per-fragment-independent, no noise model" gaps (v3 design §"Why a
  native model can out-inform") are exactly what Part B exploits.
  Sources: https://academic.oup.com/bioinformatics/article/29/24/3199/193560 ·
  https://academic.oup.com/nar/article/47/W1/W295/5480903 ·
  https://www.mcponline.org/article/S1535-9476(24)00088-4/fulltext ·
  https://www.nature.com/articles/s41467-023-40129-9 · https://github.com/Nesvilab/MSBooster

## A.3 PeptideART (feed-forward ANN) and MassAnalyzer (physics/kinetic) — the two contrasts MS2PBPI/MS2PIP beat
- **PeptideART** (Arnold et al.; "machine-learning approach to explore spectra intensity patterns").
  Ensemble of **feed-forward neural networks**, one multi-output net per important fragment, modelling
  normalized peak intensities directly; trained on ~41k PSMs. Higher accuracy than MassAnalyzer but
  data-hungry and superseded by tree models then by deep recurrent nets.
  Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2529326/
- **MassAnalyzer** (Zhang, *Anal. Chem.* 2004/2005). A **kinetic / physical-chemistry simulation** of peptide
  fragmentation built on the **mobile-proton hypothesis**: proton mobility + charge-remote cleavage rates
  generate the spectrum from first principles (no learning). Strong interpretability; weaker raw accuracy than
  data-driven models, and brittle to instrument specifics.
- **andes:** We are firmly in the data-driven (MS2PBPI/MS2PIP) camp. **Borrow from MassAnalyzer the *priors*,
  not the simulator:** proton-mobility class, enhanced cleavage N-terminal to Pro and C-terminal to D/E,
  suppressed cleavage around basic residues — encode these as **monotonic GBDT constraints / features**
  (B.4) so a thin corpus generalizes the way physics says it should. **Why the simulator doesn't apply:** a
  full kinetic model is neither pure-Rust-cheap nor more accurate than a conditioned GBDT.

## A.4 PepNovo / BoostRank — Frank, *J. Proteome Res.* 2009, 8(5):2226-2240 ("Predicting Intensity Ranks of Peptide Fragment Ions") — **the most directly relevant to andes's rank core**
- **Approach.** Predicts **intensity RANKS, not values**, with **RankBoost** (a boosting learning-to-rank
  algorithm). Rationale (quoted intent): peptide fragmentation involves competing chemical pathways that make
  *generative* probabilistic models hard, so use a **discriminative ranking model** instead.
- **Features.** Simple sequence-based features: terminal/adjacent amino acids, amino-acid composition, relative
  peak location, fragment-type-specific (b/y) features.
- **Target.** Ordinal rank of each fragment within the spectrum (e.g. predicting which peaks land at observed
  ranks 1,3,5,10 out of ~36).
- **Training.** RankBoost for ~300,000 rounds (pairwise rank loss), with growing feature counts.
- **Key result.** Accurate rank prediction; **used inside PepNovo+ scoring and to design MRM transitions** and
  improve identification. The conceptual landmark: **ranks are more instrument/method-robust than absolute
  intensities and sidestep generative modelling.**
- **andes:** This is the strongest external validation of andes's whole philosophy (v2 design §2 non-goals
  already argues "we model ranks, instrument-robust and cheap"). **Borrow:** (i) the **pairwise/listwise rank
  objective** — for the v3 *fragment* regressor, a rank/spectral-angle loss is more faithful to the
  scale-invariant cosine we actually emit than per-fragment MSE (B.2). (ii) Validates that the existing
  **rank-LLR core need not be touched** — the intensity model is the *additive* layer, exactly as planned.
  **Why a pure rank model doesn't fully apply:** the cosine numerator wants a (relative) magnitude, not just
  an order; quantile/value regression with a rank-aware loss is the sweet spot, not pure RankBoost.
  Sources: https://pubs.acs.org/doi/abs/10.1021/pr800677f · https://pubmed.ncbi.nlm.nih.gov/19256476/

## A.5 NN models — CONCEPTS ONLY (we will not use an NN; extract transferable ideas)
All of these are RNN/Transformer/CNN intensity predictors. **The architecture is explicitly out of scope
(independence + pure-Rust + no-NN).** We mine them for *modelling ideas* that survive in a GBDT.

- **Prosit** (Gessulat et al., *Nat. Methods* 2019, 16:509-518). BiGRU encoder/decoder; **inputs = sequence +
  precursor charge + NCE** via a second "meta" encoder modulating the latent (one model covers all
  charges/energies); **6 ion predictions per position** for up to 29 positions; **trained with normalized
  spectral angle (SA) as the loss**, chosen for sensitivity. Trained on ProteomeTools (~21M synthetic
  spectra). Database-search integration gave more IDs at >10× lower FDR; generalizes to other proteases,
  DIA libraries, metaproteomes. **Transferable ideas: (1) NCE/charge as continuous conditioning covariates
  (andes has `FEAT_NCE` but currently trains it at 0.0 — wasted); (2) spectral-angle objective; (3) one
  conditioned model > many disjoint models.**
  Source: https://www.nature.com/articles/s41592-019-0426-7
- **pDeep / pDeep2 / pDeep3** (Zhou/Tan et al., 2017–2019). BiLSTM; **pDeep2 adds transfer learning for PTMs**
  (fine-tune a base model on small modified-peptide data → >80% of PTM spectra reach Pearson > 0.9);
  pDeep3 BiLSTM, multi-activation (HCD/ETD/EThcD). **Transferable idea: transfer learning = train a big base
  model, fine-tune a small correction per regime/PTM** — the data-efficient answer to thin per-slug corpora
  (B.8). Sources: https://pubmed.ncbi.nlm.nih.gov/31283184/ ·
  https://pmc.ncbi.nlm.nih.gov/articles/PMC11165591/
- **DeepMass:Prism** (Tiwary et al., *Nat. Methods* 2019). RNN intensity predictor; productionized as a
  Google service. Same conceptual family as Prosit. Source: https://www.nature.com/articles/s41592-019-0428-5
- **Predfull** (Liu et al., 2020). CNN predicting the **full spectrum on an m/z grid** rather than a fixed b/y
  list — captures **internal ions, neutral losses, immonium** without an explicit ion vocabulary.
  **Transferable idea: don't restrict to b/y — model a/neutral-loss/internal ions** (v3 box 4). For a GBDT we
  keep an explicit (extended) ion list rather than a grid, but the lesson is the *vocabulary breadth*.
- **AlphaPeptDeep / peptdeep** (Zeng et al., *Nat. Commun.* 2022). Modular framework; **MS2 model only ~4M
  params (vs 64M Prosit-Transformer), ~40× faster inference**, generic chemical-element embedding for
  arbitrary PTMs. **Transferable idea: small models are enough — accuracy isn't bottlenecked by capacity but
  by training signal and conditioning**, which argues a GBDT can be competitive *as a rescoring feature* even
  if it loses raw-SA benchmarks. Source: https://www.nature.com/articles/s41467-022-34904-3
- **timsTOF-Prosit fine-tune** (Wilhelm/Wassef et al., *Nat. Commun.* 2024). Fine-tuning the intensity model
  on 302k synthetic **non-tryptic** peptides gave **up to 2.8× more HLA-I IDs** on timsTOF. **Transferable
  idea: the biggest ID gains come from matching the model to the regime (instrument + cleavage chemistry)** —
  reinforces andes's per-slug / per-regime plan and the enzyme-aware corpus.
  Source: https://www.nature.com/articles/s41467-024-48322-0

**Why the NN itself stays out of scope (state plainly):** a trained NN is not pure-Rust-inferable without
a heavy runtime, would re-introduce an external-weights dependency, and breaks the independence goal. The
*ideas* above (conditioning, SA loss, transfer/residual, vocabulary breadth) all transplant into the GBDT.

## A.6 The rescoring-integration pattern (MS²Rescore / MSBooster / Percolator) — andes's exact pattern
The universal recipe: **predicted spectrum → similarity metric(s) vs observed → additive feature(s) → semi-
supervised target/decoy classifier (Percolator) → FDR.** Similarity metrics used in the wild: cosine /
spectral angle, Pearson, **spectral entropy** (MSBooster default), top-N peak intersection, separate
b/y/both variants. **This is precisely `intensity_signal` (cosine) + the `--score strong` fusion in
`strong_score.rs`, emitted as PIN columns.** Two cheap wins fall straight out: emit **b-only/y-only** cosines
(MS²Rescore) and an **entropy-based similarity** (MSBooster) as extra additive columns — Percolator picks the
useful ones; additive features are the proven-safe pattern in our memory (parity-tuning lessons).
Sources: https://www.mcponline.org/article/S1535-9476(24)00088-4/fulltext ·
https://github.com/compomics/ms2rescore · https://www.nature.com/articles/s41467-023-40129-9

---

# PART B — Novel concepts to push andes PAST these models

For each: idea · why it can beat point-intensity regression · how to implement against our files · value/effort · risk.
"Strong evidence" = literature-backed; "andes-first" = promising but unproven.

## B.1 Discriminative-objective training (train for target/decoy separation, not intensity R²) — **andes-first, highest ceiling**
- **Idea.** Every model in Part A minimizes intensity error (MSE/SA). The metric that actually moves PSMs is
  **separation between target and decoy PSMs at 1% FDP.** Train a thin component (or reweight the regressor's
  objective) so the *emitted feature* maximizes target/decoy AUC. Concretely: keep the v3 intensity regressor
  as a feature generator, but add a second pass — a small GBDT **classifier** whose label is target(1)/decoy(0)
  PSM and whose features are the v3 similarity components (cosine, b-only, y-only, entropy, the LLR of B.7).
  This is a **learning-to-rank over PSMs** (the exact thing Percolator does, but we can pre-shape the feature).
- **Why it beats point regression.** Intensity R² and PSM discrimination are correlated but not identical;
  a fragment whose intensity is hard to predict but *diagnostic* (e.g. a proline-effect peak) should be
  up-weighted for discrimination even if it hurts global R². BoostRank (A.4) already showed discriminative
  rank training beats generative intensity modelling for *identification*.
- **Implementation.** We already have BOTH trainers: `train_gbdt` (logistic classifier, with AUC gate at
  `train.rs:100`/`:396`) and `train_gbdt_regression`. Build a PSM-level dataset (positives = gold target PSMs,
  negatives = top decoy PSMs) with the v3 similarity components as features; reuse `train_gbdt`. **Caveat
  (critical):** this is a *rescoring* feature that then goes into Percolator, so we must avoid double-dipping —
  train the discriminative head on a **separate fold / decoy split** from the one Percolator sees, or it
  inflates FDR. The memory's entrapment-FDP gate is the honest check.
- **Value/effort.** Value: **high** (directly optimizes the objective in `feedback_andes_objective_function`).
  Effort: medium (both trainers exist; the work is dataset plumbing + fold discipline).
- **Risk.** medium-high: overfitting to decoys / FDR inflation if folds leak; must gate on true entrapment FDP.

## B.2 Listwise spectral-angle / cosine loss (optimize the feature we actually emit) — **strong evidence (Prosit)**
- **Idea.** `intensity_signal` emits a per-spectrum **cosine**; the regressor is trained with **per-fragment
  MSE** on `ln(obs/base)` (`frag_dataset.rs` + `train_gbdt_regression`). These objectives disagree: MSE on
  log-intensity over-weights tiny/!matched fragments and is not scale-invariant, whereas cosine is. Prosit
  switched to **normalized spectral angle precisely for this reason** and got big gains.
- **Why it beats MSE.** Optimizing the deployment metric removes the train/serve objective mismatch; the model
  spends capacity making the *direction* of the predicted intensity vector right, which is all cosine cares about.
- **Implementation (pure-Rust, no NN).** Two tractable routes inside the existing boosting loop:
  (1) **Per-spectrum normalization in the target/gradient.** Group rows by `groups` (already run+peptide) and,
  each round, normalize predictions to unit L2 per group before computing the residual — i.e. a **listwise
  gradient**: `grad_i = (p̂_i/‖p̂‖ − o_i/‖o‖) · ∂(p̂_i/‖p̂‖)/∂p̂_i`. This is a modest change to the gradient
  computation at `train.rs:614` (the OLS gradient loop) — keep trees/leaves the same.
  (2) Cheaper proxy: train on **per-spectrum-normalized targets** (`o_i/‖o‖` instead of raw `ln(obs/base)`) and
  keep MSE — most of the benefit, near-zero code. Start with (2), graduate to (1) if it pays.
- **Value/effort.** Value: **medium-high**. Effort: low (proxy) / medium (true listwise gradient).
- **Risk.** low. The gate metric (`pearson_r2` at `train.rs:441`) should be supplemented with a held-out
  spectral-angle metric so we measure what we optimize.

## B.3 Distributional / quantile intensity regression → uncertainty-weighted cosine — **strong evidence (GBDT quantile is standard) + andes-first weighting**
- **Idea.** Predict intensity **quantiles** (e.g. q10/q50/q90) with **pinball (quantile) loss** instead of a
  point estimate. The interval width q90−q50−q10 is a **per-fragment uncertainty**. Then compute an
  **uncertainty-weighted cosine** in `intensity_signal`: down-weight high-variance fragments
  (`w_i = 1/(1+width_i)`), so confidently-predicted ions dominate the similarity.
- **Why it beats point regression.** Fragment intensities are heteroscedastic (proline/charge-remote sites are
  predictable; mobile-proton interior cleavages are noisy). A point cosine treats a wild-card fragment the same
  as a reliable one; an uncertainty-weighted cosine is a **principled, per-fragment confidence** — strictly
  more information than MS2PIP's flat cosine.
- **Implementation.** `train_gbdt_regression` currently uses OLS gradient (`grad = pred − y`, `hess = 1`,
  `train.rs:614`). Pinball loss for quantile τ has gradient `grad_i = (pred_i > y_i ? τ−1 : τ)` with constant
  hessian — a ~3-line change; train 2–3 models (one per quantile) into the model store. In `strong_score.rs`,
  add the `w_i` weights to `spectral_cosine_similarity` (a weighted-cosine variant) and feed the median as the
  predicted intensity.
- **Value/effort.** Value: **medium-high** (novel feature most rescorers lack). Effort: medium (3 models +
  weighted cosine + store plumbing).
- **Risk.** medium: 3× training cost; the weighting scheme needs tuning; gate on FDP.

## B.4 Monotonic / physics constraints in the GBDT — **strong evidence (XGBoost/LightGBM support it; thin-corpus benefit documented)**
- **Idea.** Constrain the regressor to be **monotone** in physically-grounded features so it generalizes on
  thin per-slug corpora: e.g. y-ion intensity monotone **increasing** with C-terminal basicity (R/K), b/y
  enhancement monotone with the proline-effect flag, intensity monotone **decreasing** with local proton
  competition. MassAnalyzer's physics (A.3) supplies the priors; monotone GBDTs are a documented regularizer
  that "reduce overfitting… particularly valuable when you don't have a lot of data."
- **Why it beats unconstrained regression.** Per-slug corpora are small (the memory's recurring "corpus is the
  substrate" problem). Monotone constraints inject domain knowledge as a prior → better extrapolation to
  residue contexts unseen in a thin slug.
- **Implementation.** This is the most invasive: `fit_tree` (`crates/model-train/src/gbdt/tree.rs`, not shown)
  must reject splits that violate a per-feature monotone direction (track min/max bounds per node, the standard
  XGBoost algorithm). Add a `monotone: Vec<i8>` to `TreeParams`. Requires care to keep determinism
  (`deterministic_for_same_seed` test must still pass).
- **Value/effort.** Value: medium (mostly a thin-corpus / generalization win, not a fat-corpus win).
  Effort: **high** (core tree-builder change).
- **Risk.** medium: wrong monotone direction *hurts*; needs ablation per feature. Defer until B.1/B.2/B.6 land.

## B.5 Multi-task GBDT: jointly predict intensity + presence — **andes-first**
- **Idea.** The refuted v1 was a **peak-presence classifier** (`train_gbdt`, AUC 0.8949 gate); v3 is an
  **intensity regressor**. Unify them: predict **P(fragment observed)** *and* **E(intensity | observed)**.
  The cosine numerator should use `P(present) · E(intensity)` so that a fragment the model expects to be
  *absent* contributes ~0 even if a noise peak happens to sit at its m/z — a built-in chance-match guard at
  the fragment level.
- **Why it beats intensity-only.** MS2PIP/Prosit predict intensity unconditionally; a predicted-bright ion that
  is simply *missing* in a true spectrum (common: incomplete fragmentation) is penalized, while a decoy that
  coincidentally has a peak there is rewarded. Separating presence from magnitude fixes both. Alpha-Frag
  (presence-only) showed presence prediction alone improves IDs — combining is strictly richer.
- **Implementation.** We already have both trainers. Train the v1 classifier (presence) on ALL enumerated ions
  (matched=1/unmatched=0) and the v3 regressor (intensity) on matched ions only — *exactly the two datasets we
  already build*. In `intensity_signal`, predicted intensity `:= sigmoid(presence_logit) · exp(reg_value)`.
- **Value/effort.** Value: **medium-high** (resurrects sunk v1 work; principled). Effort: low-medium (both
  models exist; only the fusion in `strong_score.rs` is new).
- **Risk.** low-medium. Calibrate the presence model (PAVA isotonic already in `train_gbdt`).

## B.6 Relational / joint structure without NN (the v3 box-2 differentiator) — **andes-first, cheap**
- **Idea.** Add features that capture fragment-ladder *correlation*, which no tree intensity model uses well:
  - **Complement coupling:** for b_i, a feature = predicted/observed presence+intensity of its complement
    y_{n−i} (same cleavage). Real ladders have correlated complements; decoy chance-matches don't.
  - **Neighbor/ladder context:** b_i conditioned on b_{i±1} presence (running-ladder evidence).
  - **Residue-pair / motif features:** the cleavage-site dipeptide (already have N/C flank *indices*; add the
    interaction explicitly), proline/Asp/Glu motif flags, and **per-flank chemical scalars**
    (basicity/hydrophobicity/helicity/pI — copy MS2PIP's `seq_<pos>_<chem>`).
- **Why it beats per-fragment-independent regression.** MS2PIP passes only the sibling's *m/z* (`ion_mz_other`),
  not its *observed evidence*; andes can pass the **observed complement/neighbor evidence at inference**,
  turning the similarity into a joint-ladder likelihood. This is information MS2PIP/Prosit structurally lack.
- **Implementation.** Extend `frag_features.rs` (currently 10 features). Add: `FEAT_NFLANK_CHEM_*` (4),
  `FEAT_COMPLEMENT_OBS` (observed relative intensity of y_{n−i} — requires the spectrum at feature-extraction
  time, so this variant is computed in `strong_score.rs`/`compute_psm_features` rather than the pure
  peptide-only `extract_frag_features`), `FEAT_NEIGHBOR_OBS`. **Important:** features that depend on the
  observed spectrum can't go in the training oracle `extract_frag_features` as-is (it's peptide-only by
  design); split into "peptide-only" (trainable) vs "observed-context" (inference-time) feature groups, and
  the observed-context ones become **additional PIN columns** rather than regressor inputs. Pure-sequence ones
  (chem scalars, residue-pair, motif) go straight into `extract_frag_features` + `frag_dataset.rs`.
- **Value/effort.** Value: **high** for the pure-sequence chem features (low effort), **medium-high** for the
  observed-complement PIN feature. Effort: low (chem) / medium (complement/neighbor plumbing).
- **Risk.** low. Additive features are the proven-safe pattern.

## B.7 Noise/chance-match likelihood ratio — formalize predicted-intensity ÷ chance-match — **andes-first, the genuine differentiator**
- **Idea.** andes is the only model here with a **noise/chance-match denominator** (`ChanceMatchSurprise`,
  `mass_competition_evidence`, `local_peak_density` in `strong_score.rs`). Promote it from a spectrum-level
  subtraction to a **per-fragment log-likelihood ratio**:
  `LLR_i = ln P(observed peak at m/z_i | real fragment of predicted intensity) − ln P(observed peak | chance)`
  where the chance term uses the local peak density ρ (peaks/Da) already computed, and the real term uses the
  predicted intensity (and its uncertainty from B.3, and presence from B.5). Sum LLR_i over matched ions →
  a **fragmentation log-likelihood-ratio feature**, the natural Bayesian generalization of `fuse_strong_score`.
- **Why it beats cosine.** Cosine rewards agreement but is blind to *how surprising* the agreement is. A match
  in a sparse region is far stronger evidence than the same match in a crowded region; LLR encodes that.
  This is the core "more information than MS2PIP/MS-GF+" claim of the v3 design, made rigorous.
- **Implementation.** `strong_score.rs` already has all the pieces: predicted intensity (regressor),
  `local_peak_density(obs, DENSITY_HW)`, ambiguity count. Add a `frag_llr(pred_int, obs_int, rho, tol)`
  returning the per-ion LLR; sum it as a new additive PIN column next to `IntensitySignal`. Subtract the
  base-rate / prior odds (the design's open calibration TODO at v3 §"Prerequisites").
- **Value/effort.** Value: **high** (unique, directly extends our differentiator). Effort: medium (the math +
  base-rate calibration). 
- **Risk.** medium: the per-fragment chance model must be calibrated on entrapment data (B.10) or the LLR is
  mis-scaled; but Percolator is robust to monotone rescalings, so even an approximate LLR usually helps.

## B.8 Base + per-regime residual / transfer adaptation — **strong evidence (pDeep2 transfer, timsTOF-Prosit, MS2PBPI cliff is the warning)**
- **Idea.** Instead of N independent thin per-slug models (MS2PBPI's failure mode — its TOF cliff), train **one
  big base regressor** on all regimes, then a **small per-slug residual/correction model** fit on that slug's
  (thin) data. Prediction = `base(x) + residual_slug(x)`. This is gradient-boosting-native: continue boosting
  the base model for a few extra rounds on the slug's data (warm-start), which is *literally what our boosting
  loop already does* — just initialize `raw_train`/`raw_val` from the base model's predictions instead of 0.
- **Why it beats per-slug-from-scratch.** Thin slugs (LysC, TMT-CID, ETD) can't support a full model; transfer
  learning is exactly how pDeep2 (PTMs) and timsTOF-Prosit (non-tryptic, 2.8× gain) solved this. Borrows the
  shared fragmentation chemistry and only learns the slug-specific delta.
- **Implementation.** In `train_gbdt_regression`, add an optional `base_model: Option<&GbdtPeakModel>` whose
  per-row predictions seed `raw_train`/`raw_val` (instead of `vec![0.0…]` at `train.rs:596`). Store both the
  base trees and the residual trees (concatenate the `trees` vec — the SoA walker sums them all anyway).
- **Value/effort.** Value: **high** for thin slugs (TMT/LysC/ETD are the andes weak spots per memory). Effort:
  low-medium (warm-start is a tiny change to the init).
- **Risk.** low. Worst case = base-only (residual contributes nothing). Fits the memory's per-regime story.

## B.9 Conformal prediction for calibrated per-fragment intervals — **strong method, andes-first in proteomics**
- **Idea.** Wrap the regressor (or B.3's quantiles) in **conformalized quantile regression (CQR)** to get
  *distribution-free, finite-sample-valid* per-fragment intervals, then use the interval width in the weighted
  cosine (B.3). CQR calibrates the quantile model on a held-out split so coverage is guaranteed.
- **Why it beats raw quantiles.** GBDT quantiles are not calibrated (the τ-quantile rarely has exactly τ
  coverage); CQR fixes coverage with a single held-out conformity quantile — cheap and theoretically clean.
- **Implementation.** We already hold out a group-disjoint validation set (`val_rows`). Compute conformity
  scores `s_i = max(q_lo − y, y − q_hi)` on val, take the (1−α) quantile, widen intervals by it. ~15 lines on
  top of B.3.
- **Value/effort.** Value: medium (refinement of B.3, not standalone). Effort: low (given B.3). 
- **Risk.** low. **Only pursue after B.3 shows the weighted cosine helps** — otherwise it's polishing an
  unused knob.

## B.10 FDR-calibration of the final fused feature (isotonic/Platt on entrapment data) — **strong evidence, near-free**
- **Idea.** Map the final fused `strong` feature (and the B.7 LLR) through an **isotonic** (PAVA — we already
  have `pava` in `model-train/src/gbdt/isotonic.rs`) or Platt calibration fit on **entrapment** PSMs, so the
  feature reads as a monotone, well-scaled evidence value before Percolator.
- **Why.** Percolator is invariant to monotone transforms of a single feature, so the gain here is mostly
  cross-feature comparability and stabilizing the per-spectrum z-scoring (`strong_score_zscore`). Most useful
  when the LLR (B.7) is mis-scaled.
- **Implementation.** Reuse `pava`; fit on entrapment splits per the memory's true-FDP harness.
- **Value/effort.** Value: low-medium. Effort: low (PAVA exists).
- **Risk.** low. Don't over-invest; Percolator already does most of the work.

## B.11 Other genuinely-novel directions (briefly)
- **Hashed k-mer sequence kernels for context without embeddings (andes-first).** NN models get context from
  learned embeddings; a GBDT can get cheap context via **feature-hashing of flanking k-mers** (e.g. hash the
  ±3 window dipeptides into a fixed-width bucket count vector). Captures motif effects no per-residue feature
  does, with zero embedding machinery. Effort: low-medium; risk: hashing collisions add noise — gate it.
- **Collision energy as a learned continuous covariate (strong evidence; near-free).** Prosit conditions on
  NCE. andes has `FEAT_NCE` but `frag_dataset.rs`/`intensity_signal` **hardcode `nce=0.0`** (see
  `frag_dataset.rs:69` and `strong_score.rs:102`). **Plumb the real parsed NCE through** so the one feature
  earns its slot — immediate, near-zero-cost win for HCD/Astral where NCE varies. This is arguably a *bug fix*,
  not a research idea, and should ship first.
- **Co-fragmentation / chimera-aware intensity (andes-first).** andes already has `--chimeric` two-pass. For
  chimeric spectra, the observed intensity at a fragment m/z is a *mixture* of co-isolated precursors;
  modelling the expected interference (predicted intensities of the co-isolated peptide) as a per-fragment
  *competition* term would denoise the cosine for the primary PSM. Effort: high; high novelty; defer.
- **Mass-error-aware matching folded into the likelihood (strong evidence; cheap).** v3 box 5. The LLR (B.7)
  should weight a match by how well the **observed m/z error** fits the per-slug error model (andes already
  models mass error / has `param.mme`). A 0.2-ppm match is stronger evidence than a 15-ppm match at the
  tolerance edge; fold a Gaussian-in-ppm term into `frag_llr`. Effort: low; pairs naturally with B.7.

---

# PART C — Recommendation: prioritized roadmap beyond the v3 MVP

Ranked by (expected PSM-lift-at-honest-FDP ÷ effort), respecting no-NN / own-model. The MVP = the bare
per-fragment intensity regressor feeding the existing cosine (`intensity_signal`). Everything below is *on top*.

### Tier 0 — fix/ship first (hours–days; not really "research")
0a. **Plumb real NCE** through `frag_dataset.rs` and `intensity_signal` (currently hardcoded `0.0`). The
   feature already exists; it's dead. *Strong evidence (Prosit). Near-zero effort.*
0b. **Per-spectrum-normalized regression target** (B.2 route 2): train on `o_i/‖o‖`, keep MSE. Aligns the
   train objective with the deployed cosine. *Strong evidence (Prosit SA). Low effort.*
0c. **Add a held-out spectral-angle gate metric** alongside `pearson_r2` so we measure what we deploy.

### Tier 1 — high lift / modest effort (the core of v3.1)
1. **Per-flank chemical-property features** (B.6: basicity/hydrophobicity/helicity/pI) into
   `extract_frag_features`. Copies MS2PIP's most valuable features; lets the trees stop relearning chemistry.
   *Strong evidence. Low effort.*
2. **Per-fragment LLR feature** (B.7): predicted-intensity ÷ chance-match, summed, as a new additive PIN
   column. This is andes's genuine differentiator and the v3 design's central thesis made concrete. Fold in
   the **mass-error-in-ppm** term (B.11). *andes-first. Medium effort, high ceiling.*
3. **Base + per-regime residual / warm-start transfer** (B.8) for thin slugs (TMT/LysC/ETD). Tiny change to
   the boosting init; directly attacks the memory's documented per-regime weak spots. *Strong evidence
   (pDeep2/timsTOF-Prosit). Low-medium effort.*
4. **b-only / y-only cosines + spectral-entropy similarity** as extra additive PIN columns (A.2/A.6). Free
   ride on the MS²Rescore/MSBooster recipe; Percolator selects. *Strong evidence. Low effort.*

### Tier 2 — high ceiling, more work / more risk (v3.2)
5. **Multi-task presence × intensity** (B.5): resurrect the v1 classifier and fuse `P(present)·E(intensity)`.
   *andes-first. Low-medium effort, principled.*
6. **Listwise spectral-angle gradient** (B.2 route 1) if route 2 paid off. *Strong evidence. Medium effort.*
7. **Quantile regression → uncertainty-weighted cosine** (B.3). *Standard method, andes-first weighting.
   Medium effort.*
8. **Discriminative head** (B.1): a GBDT trained on target/decoy separation over the similarity components.
   Highest ceiling but the most FDR-discipline risk — **must** train on folds disjoint from Percolator's and
   gate on true entrapment FDP. *andes-first. Medium-high effort/risk.*

### Tier 3 — defer until Tier 1–2 measured
9. **Monotonic physics constraints** (B.4) — only if thin-corpus generalization is the bottleneck after
   transfer learning (B.8); it's a core tree-builder change. *Strong method, high effort.*
10. **Conformal calibration** (B.9) — only if B.3's weighted cosine helps.
11. **Hashed k-mer context** (B.11), **chimera-aware competition** (B.11) — exploratory.
12. **FDR isotonic calibration** (B.10) — fold in whenever the LLR (item 2) lands.

### Strong-evidence vs andes-first (honesty ledger)
- **Strong literature evidence (someone has shown the lift):** NCE conditioning, spectral-angle objective,
  chemical-property features, transfer/residual adaptation, the rescoring-similarity pattern incl.
  entropy/b-only/y-only, quantile GBDT as a method, monotone constraints as a regularizer.
- **andes-would-be-first (promising, unproven — measure, don't assume):** the per-fragment predicted÷chance
  LLR feature, multi-task presence×intensity fusion, discriminative-objective head, uncertainty-weighted
  cosine, hashed-kmer context, chimera-aware competition. These are exactly where the v3 design's "more
  information than MS2PIP/MS-GF+" claim lives — and where, per `feedback_andes_objective_function`, the
  payoff (beat the field on PSMs) is largest. **All must clear the honest 1% true-entrapment-FDP gate before
  any public-release PR** (v3 design "The gate"); none of them are assumed to work.

---

## Appendix — sources
- MS2PBPI: Zhou et al., *Anal. Chem.* 2014. https://pubs.acs.org/doi/10.1021/ac501094m · https://pubmed.ncbi.nlm.nih.gov/25032905/
- MS2PIP (RF): Degroeve & Martens, *Bioinformatics* 2013. https://academic.oup.com/bioinformatics/article/29/24/3199/193560
- MS2PIP web server 2015: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4489309/
- MS2PIP (XGBoost / MS2PIPc): Gabriels et al., *NAR* 2019. https://academic.oup.com/nar/article/47/W1/W295/5480903
- MS²Rescore: Declercq et al., *MCP* 2022/2024. https://www.mcponline.org/article/S1535-9476(24)00088-4/fulltext · https://github.com/compomics/ms2rescore
- MSBooster: Yang et al., *Nat. Commun.* 2023. https://www.nature.com/articles/s41467-023-40129-9 · https://github.com/Nesvilab/MSBooster
- PeptideART: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2529326/
- BoostRank / PepNovo: Frank, *J. Proteome Res.* 2009. https://pubs.acs.org/doi/abs/10.1021/pr800677f · https://pubmed.ncbi.nlm.nih.gov/19256476/
- Prosit: Gessulat et al., *Nat. Methods* 2019. https://www.nature.com/articles/s41592-019-0426-7
- pDeep2 (transfer learning): Zeng et al., 2019. https://pubmed.ncbi.nlm.nih.gov/31283184/
- DeepMass:Prism: Tiwary et al., *Nat. Methods* 2019. https://www.nature.com/articles/s41592-019-0428-5
- AlphaPeptDeep/peptdeep: Zeng et al., *Nat. Commun.* 2022. https://www.nature.com/articles/s41467-022-34904-3
- timsTOF Prosit fine-tune (non-tryptic, 2.8× HLA): *Nat. Commun.* 2024. https://www.nature.com/articles/s41467-024-48322-0
- Systematic DL fragment-intensity assessment (Prosit/DeepMass/pDeep3/AlphaPeptDeep/Predfull): *J. Proteome Res.* 2024. https://pmc.ncbi.nlm.nih.gov/articles/PMC11165591/
- GBDT quantile/pinball + conformal (CQR): https://arxiv.org/pdf/2304.11732 · https://scikit-learn.org/dev/auto_examples/ensemble/plot_gradient_boosting_quantile.html
- Monotonic constraints in GBDT: https://xgboost.readthedocs.io/en/latest/tutorials/monotonic.html
