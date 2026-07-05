# Glyco fragment-intensity model — design (2026-07-05)

## Why this, and why now
Every "re-rank the existing candidates" lever is refuted with clean, review-verified
controls (Y-ladder-primary, glyco b/y rank model, 2-pass Percolator re-collapse,
2D-FDR ×2 — see `spb-design.md`). The reason is structural: the correct-vs-wrong
BACKBONE signal is in neither the sparse peptide b/y ions nor the target/decoy
labels. The one lever the refutations do NOT touch is a **generative spectral-fit
score**: predict the fragment-intensity pattern a (peptide, glycan) hypothesis
should produce, and score each hypothesis by how well it explains the OBSERVED
spectrum. A wrong backbone predicts the wrong Y-ladder → poor fit; the correct one
fits. That is a different mechanism from ranking or TD competition.

Crucially this is NOT a new subsystem — andes already has the machinery.

## What andes already has (reuse verbatim)
- `crates/scoring/src/intensity_model.rs` — `IntensityModel`: context-conditioned
  relative-intensity lookup with sparse-key backoff. Key = `IntensityKey {
  ion_type, flank_n, flank_c, pos_bin, charge, nce_bin }`. Trained by
  `andes train-intensity` → `intensity_model.parquet`.
- `crates/scoring/src/frag_features.rs` + the GBDT frag regressor (`GbdtPeakModel`,
  trained by `andes train-intensity-gbdt`) — the richer v3 model the STRONG score
  uses (`strong_score.rs`, S1 = context-intensity spectral similarity).
- The glyco pipeline ALREADY consumes the strong score (`IntensitySignal`,
  `FragPredExplained` are glyco-PIN columns).
- Glyco ion m/z are already defined: `glycan_mass::CORE_OXONIUM_MZ[5]`,
  `CORE_Y_STEPS[5]`, `oxonium::{NEUAC,NEUGC}_OXONIUM_MZ`.

## THE two limits to fix
1. **Vocabulary.** `IntensityIonType` is `{B, Y}` only ("plain b/y only in Phase T").
   In glyco spectra b/y are the SPARSE, suppressed ions (the reason the rank model
   failed); the model has no vocabulary for the RICH ions — the Y-ladder (Y0, Y1,
   core-Y rungs) and oxonium — which carry the backbone/glycan signal.
2. **Regime.** It is trained on regular DDA, so even its b/y predictions are
   mismatched to stepped-HCD glyco.

## Design

### 1. Ion vocabulary (extend `IntensityIonType`)
Add glyco ion types:
- `Y0` (bare peptide backbone), `Y1` (peptide + 1 HexNAc) — the two most diagnostic
  backbone-mass anchors.
- `YCore(k)` for k in 0..len(CORE_Y_STEPS) — peptide + cumulative core-glycan
  (trimannosyl) rungs.
- `Oxonium(i)` — the core oxonium set (204/138/168 HexNAc, 366 HexHexNAc, 512, 657)
  + sialic (274/292 NeuAc, 290/308 NeuGc).
Keep `B`, `Y` for the (sparse but present) peptide backbone ions.

### 2. Context key (glyco ions condition on COMPOSITION, not peptide flanks)
b/y intensity is peptide-local (flank residues). Glyco-ion intensity is
GLYCAN-conditioned + energy-conditioned. Split the key by ion class:
- b/y ions: existing `IntensityKey` unchanged.
- glyco ions: a NEW key `GlycoIntensityKey { glyco_ion_type, glycan_class,
  precursor_charge, nce_bin }`, where `glycan_class` is a coarse composition bucket
  (high-mannose / complex / hybrid / sialylated / fucosylated) so the model stays
  low-cardinality with backoff (drop glycan_class → drop nce_bin → drop charge).
  Rationale: oxonium/Y-ladder intensities depend on sialic/fuc content and NCE, not
  on the peptide sequence.

### 3. Training (extend the existing accumulator, not a new engine)
- Reuse `andes train-intensity` (or `-gbdt`) with a glyco-aware fragment enumerator:
  for a labeled (peptide, glycan) PSM, enumerate b/y + Y0/Y1/YCore + oxonium at
  their theoretical m/z (functions already exist in `andes-glyco`), match to peaks,
  record base-peak-normalised relative intensities keyed as above. Accumulate
  mean/var per key (lookup) or emit rows for the GBDT.
- Labels: the multi-dataset harvest (`40-data/pride-datasets.md`) — PXD005411
  (a glyco search engine2 mouse-brain), PXD016175 (a glyco search engine2 IgG), PXD030670 (a commercial glyco engine saliva) — PLUS
  the an open-source glyco engine (O-Pair) oracle now stood up (`40-data/multitool-truth-validation.md`)
  to label additional PXD025455 Fc-runs. TRAIN excludes the Fc3_r1 eval run.

### 4. Scoring (the strong score already spans the ions once they exist)
`strong_score.rs` computes predicted-vs-observed similarity over the fragment set.
Once the glyco ions are in the vocabulary + predicted, the strong score's
`intensity_signal` naturally spans {b, y, Y-ladder, oxonium} → a (peptide, glycan)
hypothesis is scored by how well it explains the WHOLE observed pattern. This is
the generative-fit backbone/glycan-correctness score.

### 5. Integration — attack the COLLAPSE (the bottleneck)
The per-scan collapse currently selects by b/y `rank_score` (the ranking bottleneck:
51/196 backbone-correct). Add a `collapse_cmp` mode `ANDES_GLYCO_SELECT=strong` that
selects by the glyco STRONG score (generative fit) instead of / blended with b/y
rank. The strong score also enters Percolator as a feature (already partly there).
Because the strong score fits the RICH Y-ladder (not sparse b/y), it can pick the
correct backbone where rank_score cannot.

## S0 FEASIBILITY GATE (2026-07-05) — PASSED (signal is learnable + discriminative)
Before any harvest/code, checked whether glyco-ion relative intensities are
predictable, using the an open-source glyco engine matched-ion intensities already on the VM (222
confident N-glyco PSMs, `feasibility_intensity.py`). Result:
- LEARNABLE: oxonium D-ion relative intensities are highly CONSISTENT within a
  glycan composition — CV ≈ 0.05–0.2 for D204/D138/D168/D366/D274/D292, present in
  ~100% of PSMs (per composition: H5N4A2F1 n=24, H6N5A3 n=20, H5N4A1 n=17,
  H6N5A3F1 n=16). Robust statistics — the OPPOSITE of the sparse-b/y that sank the
  rank model.
- DISCRIMINATIVE: NeuAc oxonium (D274+D292) = 12.1% mean, present 190/190 (100%) in
  sialylated (A>0) vs 0.00%, present 0/32 (0%) in non-sialylated. Perfect class
  separation.
⇒ The design's core bet holds: glyco ions (oxonium + Y-ladder) are abundant,
low-variance, class-discriminative → a class-conditioned intensity model WILL have
signal. Proceed to S1/S2 with confidence. (The Y-ladder / (M0−loss) ions carry the
BACKBONE-mass signal similarly; oxonium checked here as the composition-diagnostic
half.)

## Staged plan + gates
- **S1 — Data.** Harvest PXD005411/PXD016175/PXD030670 result tables; run the
  an open-source glyco engine oracle on 2–3 non-Fc3_r1 PXD025455 runs. Assemble a labeled
  (peptide, glycan, spectrum) corpus, Fc3_r1 held out. Gate: ≥ a few k confident
  glyco PSMs spanning glycan classes + NCE.
- **S2 — Vocabulary + accumulator.** Extend `IntensityIonType` + the glyco fragment
  enumerator + the training accumulator. TDD the enumerator (known m/z). Gate:
  round-trip a labeled PSM → expected ion intensities.
- **S3 — Train.** Produce the glyco intensity model (lookup first — simplest;
  GBDT if it pays). Gate: predicted intensities correlate with held-out observed
  (per-ion-class calibration).
- **S4 — Score + integrate.** Wire the glyco strong score; add
  `ANDES_GLYCO_SELECT=strong` collapse mode.
- **S5 — Validate.** On Fc3_r1 vs the 196 consensus: **success = backbone-correct
  > 51/196** (and total @1% FDR ≥ 261), geometry/labels controlled as in the SP-B
  protocol. Cross-dataset check on a harvested run.

## Honest risks
- **Why this could beat the rank model:** the Y-ladder + oxonium are ABUNDANT (not
  sparse), so their per-class intensity statistics are robust — the opposite of the
  sparse-b/y problem that sank the rank model.
- **Why it could still fail:** if glyco-ion intensity is dominated by NCE/instrument
  variation the model can't condition on, the fit score is noisy. Mitigation:
  stepped-HCD-specific nce_bin + coarse glycan_class backoff; validate calibration
  (S3 gate) BEFORE integrating.
- **Data sufficiency:** composition × NCE × charge is high-cardinality; rely on the
  coarse `glycan_class` bucket + backoff, and the multi-dataset harvest for volume.
- **FDR boundary unchanged:** the strong score is a scoring feature; FDR stays
  Percolator on the 1-PSM/scan PIN. This does not touch the FDR model.
- **Don't repeat the confounds:** geometry not applicable (intensity model, not
  rank geometry), but KEEP the held-out eval + 2-tool consensus metric + one-variable
  A/B, and control instrument/NCE between train and test.

## Non-goals
No new FDR machinery. No borrowing a glyco search engine/a commercial glyco engine/O-Pair CODE (papers + running them
as oracles only). Not a full deep-learning predictor in v1 — start with the existing
lookup/GBDT context model extended to glyco ions; escalate to a sequence model only
if the context model's calibration ceiling demands it.
