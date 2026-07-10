# RT Prediction — SOTA review + engine-wide andes design (P2)

**Status:** review + design (2026-07-07), NO implementation. From a deep multi-agent
SOTA review. Governing principle: **engine-wide** — a general peptide/backbone RT
predictor (benefits regular search), with the glycan RT offset as a glyco extension.

## SOTA summary (what to reuse, what to skip)
- **Deep nets** (Prosit-iRT R≈1.00, DeepLC CNN-on-atomic-composition, AutoRT
  MAE 0.57 min, Chronologer harmonized-corpus RCNN, DeepRTplus CapsNet R²≈0.994)
  are the accuracy ceiling but need GPUs/large corpora and aren't clean-room-cheap.
- **Classic/GBDT tier** (SSRCalc R²≈0.95–0.98; ELUDE SVR-on-~60-features) reaches
  **R≈0.96–0.98** — parity with SSRCalc, ~30–50% higher RT error than deep nets.
  **That gap is small and largely closed by per-run calibration** when RT is just
  one PIN feature among ~50. Right cost/accuracy point for andes (own-trained,
  patent-free, pure-Rust GBDT, fast).
- RT is a **secondary** rescoring lever (fragment intensity is stronger): adds a
  **complementary ~10–17%** ID lift, decisive where MS2 is ambiguous — HLA,
  low-input, and **isobaric glycoform disambiguation**.

## Glyco RT — the additive-offset foundation (Krokhin/Klein & Zaia 2019, R²=0.995)
`RT_glyco = RT_backbone + Σ (per-monosaccharide coefficient)`. **CRITICAL sign
caveat:** on reversed-phase, glycans elute the glycopeptide EARLIER (hydrophilic)
— *except sialic acid*, which retains LATER (carboxylate ion-pairing): asialo
−1.37% → mono −0.47% → di +0.61% → tri +1.94% ACN. So **NeuAc gets its own
(positive) coefficient; do NOT lump it with neutral hexoses (negative).**
Coefficients are RP-specific (HILIC flips sign) and single-group — **re-fit on
andes's own PRIDE data, do not hardcode Krokhin's numbers.**

## andes design (reuses existing machinery — file:line confirmed)
- **Model:** GBDT regressor to a dimensionless RT INDEX. Reuse `GbdtPeakModel`
  (`scoring/gbdt_eval.rs`, `predict_value`, `apply_sigmoid=false`). Store as a keyed
  blob in the parquet ModelStore (`gbdt_model_bytes` column already exists).
- **Features (~30, engineered from `Peptide`):** 20 AA composition counts; length;
  mass; summed + mean Kyte-Doolittle hydropathy (reuse `frag_features::hydropathy`);
  **N-/C-terminal residue identity + hydropathy** (SSRCalc's biggest lever);
  strong-retention (W,F,L,I,Y,M) vs polar/charged counts; net-charge/pI proxy;
  **mod count + summed `mod_.mass_delta`** + optional Δatom-composition (DeepLC-style,
  for partial unseen-PTM generalization).
- **Per-run self-calibration (iRT mechanism):** fit `RT_obs = a·index_pred + b` on
  the run's confident 1%-FDR targets (robust linear; upgrade to piecewise only if
  residuals show gradient non-linearity). Decouples transferable chemistry from
  per-run LC geometry.
- **Glycan offset = glyco-only additive extension:** `index_glyco = index_backbone
  + Σ c_mono·count` from `GlycanComp{hexnac,hex,fuc,neuac,neugc}`; NeuAc own sign.
  Zero for regular peptides (same code path serves both).
- **PIN features (append to the single source of truth `output::pin::psm_feature_values`;
  the PIN currently has NO RT feature):**
  - `DeltaRT` (signed, minutes), `AbsDeltaRT` (tools disagree on sign — emit both),
    `DeltaRTNorm` (÷ gradient span, transfers across runs).
  - **glyco-only `DeltaRTRank`** = rank of this candidate's `AbsDeltaRT` among
    competing glycoforms sharing the backbone for that scan — the GlycReSoft
    disambiguation trick, the **highest-value glyco RT feature** (neutral when <2
    co-identified glycoforms). Null/0 when `rt_seconds` is None.
- **Units:** `Spectrum.rt_seconds` is seconds → convert to minutes; keep f64; emit
  Double. Fills the always-null `predicted_rt` (`qpx.rs:472`).

## Sequencing + validation
- **Commit 1 (engine-wide):** GBDT backbone RT index + per-run linear calibration +
  `DeltaRT`/`AbsDeltaRT`/`DeltaRTNorm` in `psm_feature_values`. Validate on Astral
  (regular high-res) — additive/top-1-preserving; gate on non-regression.
- **Commit 2 (glyco extension):** per-monosaccharide offset + `DeltaRTRank`. Validate
  on glyco backbone find-rate + @1%.
- **Confidence flags:** GBDT-vs-deep gap is *inferred* (no published GBDT-vs-DeepLC
  head-to-head); glyco ID gains (+9.7%/+15.7%) are single-group/mouse-tissue —
  mechanism sound, magnitude must be measured on andes data; RP per-monosaccharide
  coefficients single-group — re-fit, don't hardcode.
