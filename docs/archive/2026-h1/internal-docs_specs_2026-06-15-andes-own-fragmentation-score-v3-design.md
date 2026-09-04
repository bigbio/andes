# andes own fragmentation-aware scoring model (v3) — design

Date: 2026-06-15 (updated 2026-06-16)
Status: APPROVED to build (user 2026-06-16: "merge RawScore + GBDT into one" = v3). Supersedes the A2
peak-model and GBDT-v2 designs (subsumed here). Location: internal (never in the public repo).

## DECISION + STATUS (2026-06-16)
The merge = box 1 below: a peptide-CONDITIONED GBDT that predicts per-fragment expected relative intensity
(REGRESSION) REPLACES the coarse `IntensityModel` lookup table inside `RawScore`'s `intensity_signal`
cosine. One trained fragment model powers RawScore; the v1 peptide-AGNOSTIC peak classifier is subsumed.
- Prereqs (§ below): **DONE** — mod-aware labels + per-peak `param.mme` tolerance (commit 5002df10, 18
  gbdt tests). 
- v1 substrate **GATE PASSED**: held-out AUC 0.8949 (commit 01133f51) ⇒ the GBDT machinery + corrected
  labels discriminate; safe to invest in the regressor.
- Score naming (post 6f268a1a + 64710d17): **RankScore** = rank-LLR (always), **RawScore** = fused
  (= `intensity_signal − null`, the column v3 improves). v3 changes ONLY what feeds `intensity_signal`.
- Baseline to beat: Astral `--score strong` = 38,909 (RE-VERIFY after the 64710d17 clobber fix, queued
  behind the running Java-Astral).
Integration stays the proven safe pattern: additive PIN feature(s) in rescoring; RankScore byte-identical.

## Goal

andes ships its **own** fragmentation-aware scoring model that is (1) **fast**, (2) **beats Java
MS-GF+**, (3) **beats MS2PIP**, and (4) **captures more information than either** — with **no external
tool, no neural net, pure-Rust inference, own-trained models** (the campaign's independence goal). It is
NOT an MS2PIP wrapper and NOT an MS2PIP re-implementation: MS2PIP's signal (sequence-conditioned fragment
intensity) is *one component* of a richer discriminative likelihood.

## Why a native model can out-inform MS-GF+ and MS2PIP

- **MS-GF+** = a generating-function probability over the rank/match structure. No fragment-intensity
  model, no fragmentation chemistry. (andes's in-transition bundle is MS-GF+-derived.)
- **MS2PIP** = a GBDT predicting *per-fragment relative intensity* from *local sequence context*. Powerful,
  but **intensity-only, per-fragment-independent, and has no noise/chance-match model**. Its discriminative
  feature is `spectral_similarity(predicted, observed)`.
- **andes-v3** subsumes MS2PIP's intensity prediction as one term and adds, as a single fused likelihood:
  1. **Sequence-conditioned fragment expectation** (flanks ±k, residue pairs, proline flag, position,
     length, charge, mods) — *matches MS2PIP*.
  2. **Relational / joint structure** — model an ion conditioned on its **complement** (b_i ↔ y_{n−i},
     same cleavage) and **neighbors** (b_i ↔ b_{i±1}); real ladders are correlated, decoy chance-matches
     are not. *MS2PIP does not model this.*
  3. **Noise / chance-match likelihood** — andes's per-rank/peak-density noise model → the LLR
     *denominator*. Turns "predicted intensity" into "real-fragmentation-vs-noise likelihood ratio."
     *MS2PIP has no noise model.*
  4. **Full fragment vocabulary** — b/y + a + neutral losses + multiple charges.
  5. **Fragment mass-error / accuracy** — andes already models mass error; fold it in.
  6. **Per-regime specialization** — per-slug models (activation/instrument/enzyme/protocol).
  7. **Discriminative fusion** — the `--score strong` numerator/denominator: similarity numerator ÷
     (chance-match surprise + mass competition + rank entropy + listwise gap).

  Boxes 2–7 are information MS2PIP and MS-GF+ do not use. The "more information" claim is structural; the
  "better" claim is **to be measured**, not assumed (see Gate + Risks).

## Why it's fast

- Pure-Rust GBDT (the existing zero-dep SoA walker). No Python, no NN.
- Evaluated in **rescoring**: for the **top-N candidates per spectrum** (post-search), never in the
  16M-candidate DP hot loop. This is the proven, Percolator-safe pattern; the per-PSM rescoring hook
  (`compute_psm_features` in `match_engine.rs`) already exists and already has the peptide + spectrum +
  charge in hand. Additive PIN feature(s) — `RawScore` stays byte-identical.

## Architecture (leverages existing scaffolding — big head start)

Already built (Phase V) and reused:
- `IntensityModel` (flank-conditioned mean log-rel-intensity table) — the **coarse v0** of box 1; we
  *enrich* it, not replace the plumbing.
- `intensity_signal` / the `IntensitySignal` PIN column — predicted-vs-observed cosine similarity,
  per-candidate, additive. The rescoring feature hook.
- `--score strong` numerator/denominator framework (`strong_score.rs`): `strong = intensity_signal − null`.
- Rust GBDT trainer (`model-train/src/gbdt/`) + SoA walker (`gbdt_eval.rs`) + the model-store blob column.
- Fragment enumeration (`fragment_ions.rs`), observed matching (`nearest_peak_full`), `spectral_cosine_similarity`.

New (the v3 model):
- **A GBDT fragment-expectation predictor**: input = per-(annotated)-ion feature vector
  {ion-type, charge, flanks ±2, residue pair, proline flag, position-frac, length, mods, **complement
  presence/rank, neighbor presence/rank**}; output = expected relative intensity (regression) — replacing
  the coarse `IntensityModel` lookup with a richer GBDT fit. (Regression target → may need a regression
  mode in the trainer, or binned-intensity classification.)
- **A noise/chance-match likelihood** combined with the predicted intensity → a per-ion LLR (the v2 idea,
  folded in) feeding the similarity/strong-score numerator.
- The similarity (and component features) emitted as additive PIN columns; `--score strong` fuses them.

## Prerequisites / fixes (blocking)
- **Mod-aware labels** (the v1 bug): ✅ DONE (5002df10) — oracle uses modified residue masses from the
  parsed `Peptide`; no string-parse misalignment.
- **Per-peak fragment tolerance**: ✅ DONE (5002df10) — labeling uses `param.mme` (ppm-aware on high-res),
  not a hardcoded 0.5 Da scalar.
- **Calibration / base-rate** correctness for any LLR term (subtract prior odds; isotonic on a
  representative split): still required for the regressor/LLR term.

## The gate (how we know it's "better")
A/B at honest 1% true entrapment-FDP, uniform Percolator, on Astral / TMT / UPS / LysC, with wall-time:
- **Beat Java MS-GF+** (already close: own rank core beats it on Astral) — must hold with v3 on.
- **Beat the open-source field** (Sage/Comet/ProSE) on PSMs@1%.
- **Match-or-beat MS2PIP-level lift** — measured against the field leaders that *use* MS2PIP/Prosit
  rescoring (MSFragger+MSBooster is the external reference; we do not run MS2PIP in the product). If a
  one-off MS2PIP ceiling probe is ever wanted it is a *measurement*, not a dependency — but the user has
  opted out; we target beating the field directly.
- **Fast**: rescoring adds bounded per-spectrum cost (top-N only); search wall stays competitive (the
  Phase-V gate already checks "strong ≤ ~110% of rank" wall).
- No public-release PR until the gate passes.

## Honest risks
- **"Beat MS2PIP" is an empirical bar, not a given.** MS2PIP is trained on far more data; our edge has to
  come from boxes 2–7 (relational + noise + fusion) + a good own corpus. We measure at each step; if the
  fragmentation lever doesn't lift past the field, the lever is corpus (proven +8.3%) and we say so.
- **Regression-target GBDT** (intensity) is a new trainer mode vs the current classification trainer.
- **Relational features** add modeling + training-data complexity (complement/neighbor must be computed
  consistently train and infer).
- **Calibration / fusion** (strong-score numerator÷denominator) must be tuned; Phase-V built it but its
  gate result must be re-measured (the prior output dir is gone).
- Corpus quality remains the substrate — a richer model on a thin corpus still underperforms.

## Build sequence (high level; to be turned into a task-by-task plan)
1. Fix mod-aware labels + per-peak tolerance (prereqs).
2. Re-measure the EXISTING coarse `IntensityModel` + `--score strong` on Astral (baseline + the
   intensity-lever's current value), since the Phase-V result is lost.
3. Enrich the model: GBDT fragment-expectation predictor with sequence + **relational** features
   (regression/binned), replacing the coarse lookup; keep the noise-LLR + strong-score fusion.
4. Train per-slug on own gold PSMs (Codon); benchmark each gate dataset on the VM; iterate.
5. Ship the own-only store once the gate passes (Phase B PR).
