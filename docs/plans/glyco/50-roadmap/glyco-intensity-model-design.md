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

## S1 + S3 RESULTS (2026-07-05)
S1 (clean data): stood up the oracle pipeline for a non-eval run — download
`HCC_pool_Late_Fc5_r2.raw` → ThermoRawFileParser convert to mzML (the bundled
Thermo license prompt crashes an open-source glyco engine non-interactively; conversion sidesteps
it) → an open-source glyco engine N-glyco. Result: **154 confident N-glyco PSMs on Fc5_r2**
(held-out-from-Fc3_r1 training data). Pipeline is now reusable (`run_s1_clean.sh`).

S3 (cross-run calibration, oxonium): train per-glycan-class oxonium lookup on
Fc5_r2, predict held-out Fc3_r1 → **cosine median 0.989 / mean 0.983 (n=222)**.
STRONG cross-run generalization ⇒ a trained intensity model transfers.
NUANCE: the class-AGNOSTIC baseline is 0.986 — nearly as high. The core oxonium
profile is NEAR-UNIVERSAL across N-glyco spectra, so oxonium GENERALIZES but is
only weakly class-DISCRIMINATIVE in cosine (the S0 discrimination is the binary
presence/absence of the sialic markers, not a graded profile). ⇒ Oxonium confirms
"glyco spectrum of ~this class"; it does NOT pick the backbone. The
backbone-discriminative signal for the RANKING bottleneck is the Y-LADDER (peptide
+ glycan-loss ions) — must validate THAT next (S3b): does the predicted Y-ladder
separate the correct backbone from wrong ones? That is the crux, not oxonium.

S3b (Y-ladder learnability): the glycan-Y-ladder is present in 222/222 PSMs. The
DOMINANT ion (full-glycan-loss → the Y1/backbone anchor) is reliable — present
100%, CV ≈ 0.2–0.37 per composition. But the INTERMEDIATE ladder rungs are NOISY —
CV ≈ 0.4–1.0+, present 40–100%. So: the backbone-ANCHOR ion (Y1) is learnable and
reliable (good — it pins the backbone mass), but the full Y-ladder INTENSITY
PATTERN is much noisier than oxonium.

## HONEST ASSESSMENT (2026-07-05) — feasibility PASSES, net-new ranking value UNPROVEN
The feasibility gates pass: glyco ions are learnable (S0), the model generalizes
cross-run (S3, cosine 0.989), the clean-data pipeline works (S1, 154 PSMs). BUT the
nuances stack into a real caveat for the RANKING bottleneck specifically:
- Oxonium is NEAR-UNIVERSAL → confirms glyco-class, does NOT pick the backbone.
- The Y-ladder's RELIABLE part (the Y1/backbone anchor) is ALREADY used by andes
  (`glycan_y_intensity` / `YLadderScore` / the collapse tiebreak).
- The Y-ladder PATTERN beyond the anchor is noisy → hard to model precisely.
⇒ The signals are real and learnable, but much of the RELIABLE signal is ALREADY in
andes's features. Whether a calibrated generative-fit BEATS the existing YLadderScore
for backbone RANKING is UNPROVEN — it is not guaranteed by the feasibility results.

THE decisive test (do BEFORE building the full model): for truth PSMs where andes
generates the correct backbone but ranks it below a wrong one, does the Y-ladder
PATTERN-FIT (spectral angle vs a learned per-class template) score the correct
backbone ABOVE the wrong one MORE OFTEN than the current YLadderScore (sum) does?
That is a candidate-level analysis on andes's own output (emit accepted backbones +
their Y-ladder pattern-fit vs the truth). If yes → build S2–S5. If the pattern-fit
does not out-separate the existing sum → the intensity model won't move 51/196, and
the honest conclusion is that andes is at its ceiling for THIS data/fragmentation
(HCD), and the real lever is orthogonal fragmentation (EThcD, idea C) or better
GENERATION (idea A), not a richer HCD intensity model.

## ⚠ CONCLUSION REOPENED (2026-07-05, Codex adversarial review)
The "decisive test" below is METHODOLOGICALLY FLAWED and does NOT close the lever:
it compares the correct backbone to a SYNTHETIC mass-shifted decoy, but the design
called for testing against andes's REAL alternative-backbone competitors — the
specific wrong (backbone, glycan) splits of the same precursor that actually beat
the correct backbone in the collapse (the 61 consensus losses). A random mass shift
does not model those competitors. So `0.679 vs 0.682` does NOT prove the intensity
model can't help ranking. TODO before trusting "single-spectrum exhausted": rerun
pattern-vs-sum on the ACTUAL per-scan candidate sets (andes ALL_HITS emits them),
specifically the 61 wrong-backbone consensus losses — does pattern-fit reorder
those real competitors? Until then this conclusion is PROVISIONAL, not closed.

## DECISIVE TEST (2026-07-05) — [SUPERSEDED: flawed decoy model, see above]
### intensity model will NOT move the ranking bottleneck
Ran the candidate-level measurement (`pattern_vs_sum.py`, real Fc3_r1 mzML peaks):
for 222 correct PSMs, extracted the core-Y ladder for the CORRECT backbone vs a
mass-shifted DECOY backbone, and compared how well the SUM (what andes uses) vs the
PATTERN-fit (cosine to the learned template) separates them.

| metric | correct>decoy rate | rescues SUM's failures |
|---|---|---|
| SUM (andes today) | 0.679 | — |
| PATTERN-fit | 0.682 | 7/428 = 1.6% |

IDENTICAL discrimination; the pattern rescues only 1.6% of the sum's failures (it
fails in the SAME cases). Template = [Y0 0.055, Y1 0.128, Y2 0.014, Y3 0.009, Y4
0.011, Y5 0.006] — the core-Y ladder is essentially Y0+Y1 (two ions); the higher
rungs are ~0.01 (near-noise). There is NO rich intensity pattern to exploit.
And the other ions don't help BACKBONE ranking: oxonium is backbone-INDEPENDENT (a
mass-shifted backbone has identical oxonium), b/y are sparse.

⇒ MEASURED CONCLUSION: a richer HCD glyco intensity model does NOT beat andes's
existing YLadderScore for backbone ranking. The feasibility gates passed (ions
learnable, model generalizes) but the discriminative signal is Y0+Y1, already
captured. The intensity-model lever is CLOSED for the ranking bottleneck. andes is
at its HCD ceiling (261 @1% FDR > an open-source glyco engine 222; 51/196 backbone-correct).

The genuine remaining levers are NOT HCD scoring refinements:
- **Orthogonal fragmentation (EThcD/ETD, idea C)** — c/z ions pin the backbone that
  HCD's suppressed b/y + two-ion core-Y cannot. The real information gap.
- **Better GENERATION (idea A)** — fewer wrong-backbone candidates per scan, so the
  moderate (~68%) core-Y discrimination has fewer competitors to lose to.

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
