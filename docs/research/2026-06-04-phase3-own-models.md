# Phase 3 — Own models from public/MSnet data (plan)

**Goal:** replace the 39 bundled models (`cid_lowres_tryp`, `hcd_qexactive_tryp`, …) —
which are **MS-GF+'s trained `.param` tables** (the last derived *data*) — with
SIMAS-native models trained from **public** data. This removes the final derivation →
enables relicense, **and** is where the **low-res 6%** gets recovered.

## The core challenge (why this is research, not "press train")

Tonight proved the bootstrap estimator **dilutes**: a 15-file, cross-dataset-validated
PXD016999 model scored **−4.3%** vs the curated `cid_lowres` — *worse* than even the
prior one-file model. More data made it worse. So the dilution is an **estimator
property, not a data-volume problem**. Naive retraining on MSnet will not match curated
quality. We must fix the estimator first.

## Plan

### Step 1 — Diagnose the dilution *(cheap, on data we have; start here)*
Compare a bootstrap-trained model's tables to `cid_lowres`'s, term by term:
- rank-distribution tables (`rank_dist_table`): are the trained distributions **flatter /
  softer** than curated? By how much, and where (which ranks / partitions)?
- ion-existence + error tables: same comparison.
- Hypotheses to test: (a) Laplace/add-k smoothing over-softens sparse partitions;
  (b) the label set (TDC ≤1% PSMs) is biased toward what the seed already scores well
  (confirmation bias) → narrower observed-ion statistics; (c) curated tables were
  hand-tuned/sharpened beyond raw MLE.
- Output: a precise characterization of *what* differs → tells us the fix.

### Step 2 — Fix the estimator
Candidate fixes (pick by Step-1 evidence):
- **Shrinkage toward a non-derived prior** (an analytic fragmentation prior, or a
  cross-partition pooled estimate) instead of flat smoothing — sharpens sparse partitions
  without overfitting.
- **Sharpening / temperature** on the estimated rank distributions to counter the
  softening.
- **Better labels** — iterate labels with the improving model (careful: EM degraded
  before), or richer label sources (MSnet's confident reanalysis PSMs).
- Validate each candidate on held-out data (cross-dataset, not train=test).

### Step 3 — Source public/MSnet training data (**quality-prefiltered**)
Harvest confident PSMs across **activation × instrument × enzyme × label** from PRIDE/MSnet
reanalyses (Parquet results) + raw where needed. Map coverage to the 39 model slugs;
prioritize the ones that matter (CID/HCD × tryp × {none, TMT, phospho}). Note: MSnet is
mostly label-free + 2 TMT projects — supplement TMT/iTRAQ from ProteomeTools/PRIDE.

**PREFILTER to only the *really good* data — quality over quantity** (tonight proved more
data *diluted*). Two layers:
- **Dataset-level:** keep only well-characterized, high-quality reanalyses (high ID rate,
  clean instrument/acquisition, trusted submitters); drop noisy/low-yield projects.
- **PSM-level:** strict confidence (well below 1% q, high score margin / unambiguous
  rank-1), clean spectra (good fragment coverage, low co-isolation). Only these confident
  PSMs feed the estimator — a small, *sharp* training set beats a large, noisy one for
  this estimator (the dilution lesson). Define explicit selection thresholds in Step 2's
  validation.

### Step 4 — Train + validate the full model set
Train SIMAS-native models with the fixed estimator; write a fresh `models.parquet` with
**zero MS-GF+-derived rows**; retire the `.param` fixtures + (optionally) the `.param`
loader. **Gate:** yield ≥ curated on Astral + a05058 + a cross-dataset hold-out
(entrapment-validated). Recovering the low-res a05058 toward 11,128 is the headline test.

## Done = independence
GF-free (Phase 1) + own scoring code (Phase 2) + own models (Phase 3) ⇒ no MS-GF+-derived
substance ⇒ Phase 4 relicense (Apache-2.0, drop NOTICE requirement).

## First move
**Step 1 — the dilution diagnosis.** Cheap, uses the trained model from tonight
(`tmt_store.parquet`) vs `cid_lowres`, and it determines the whole estimator fix.
