# andes Rich-Ion-LLR — decoy-aware per-annotated-ion discriminative score (design)

**Date:** 2026-06-16 · **Branch:** feat/gbdt-stronger-models

**Goal:** For each matched b/y ion, score `P(real fragment of THIS peptide | assignment-aware features)` with a pure-Rust GBDT, convert to a per-ion likelihood ratio `LLR = logit(P) − logit(π)`, and sum over matched ions into one discriminative score. This fuses the assignment-consistency properties (complement coherence + cleavage chemistry first; explain-away later) **inside one GBDT** rather than as separate PIN columns. **RankScore (the stable default) is never modified;** this lives in the RawScore/strong path and as an additive PIN column.

**Tech stack:** Rust (`scoring`, `model-train`, `search`, `output`); the existing pure-Rust histogram-GBDT classifier (`model-train/src/gbdt/train.rs` logistic + `gbdt_eval.rs`); gold-PSM flats; the (now-correct) decoy generation (commit 1bbb6ff3); Percolator 3.7.1 + entrapment-FDP for validation.

---

## Why (the evidence)
- **Intensity accuracy is a dead lever** (cosine ceiling, confirmed 3× on Astral incl. de-confounded + charge-extended): a per-peak/per-fragment *quality* or *intensity* score is **isotropic** — it lifts the target and competing decoys equally, so it doesn't separate them.
- **Discrimination lives in assignment CONSISTENCY** — does *this peptide's claimed set* of peaks cohere? That is peptide-aware (non-isotropic). The cheapest, highest-conviction such signals are **complement-pair coherence** and **cleavage-pattern chemistry**, both of which a decoy structurally cannot fake.

## Non-negotiables (or this repeats the cosine null)
1. **Assignment-aware features** — every feature depends on the peptide-peak assignment (complement-of-this-peptide, cleavage-of-this-bond), not peptide-agnostic peak quality.
2. **Decoy-aware training** — positives = target-matched ions; negatives = **decoy-matched** ions (generated via the corrected decoy path), NOT signal-vs-junk.
3. **Calibrated LLR** — per-ion output = `logit(P(signal|x)) − logit(π)` (prior-subtract), so the summed score is additive (fixes the latent logit≠LLR bias).
4. **Validate on net near-threshold crossings + entrapment-FDP** — count net target PSMs crossing q=1% (minus decoys crossing), NOT held-out peak AUC (which is dominated by easy peaks and lies).

## Features — per matched annotated ion (Phase 1)
Assignment-aware, computable from `(spectrum, peptide)` alone:
- **Cleavage chemistry** (reuse `frag_features.rs`): ion type, n/c-flank residues, proline-flank, proton-mobility, basic-flank, hydropathy, position-frac, peptide length, fragment charge. *The GBDT learns cleavage propensity from these internally — no separate P(cleave) table.*
- **Complement coherence:** complement-present (is `b_i`'s `y_{n−i}` matched?), complement-ppm-delta `|ppm_b − ppm_y|`, complement-rank-agreement.
- **Local reliability (used in combination):** intensity rank, local peak density, signed mass-error ppm.

**Deferred to a later layer:** explain-away `claim_count` (needs candidate-set training data, not in static flats), isotope/neutral-loss partner support, sequential-ladder neighbor.

## Model + integration
- Reuse the pure-Rust histogram-GBDT classifier (`train_gbdt`, logistic, held-out AUC) + `gbdt_eval` SoA walker.
- Per-ion `P(signal|features)` → `LLR = logit(P) − logit(π)`; sum over matched ions → `RichIonLLR`.
- **Deploy as an additive PIN column first** (RankScore untouched; the parity-proven-safe integration). Later option: make the summed LLR the RawScore engine (replacing the dead cosine) if it wins.
- Stored as a new GBDT blob in the model store (per regime; back-compat: absent → feature emits 0.0).

## Training data (decoy-aware)
- **Positives:** matched b/y ions of gold PSMs (target peptide) from the flats.
- **Negatives:** decoy-matched ions — generate the decoy peptide the SAME way the search does (reverse/shuffle the protein → digest, per fix 1bbb6ff3), match its theoretical ions to the spectrum → those matched ions are the negatives. (MVP simplification allowed: reverse the gold peptide directly; flagged as a proxy to revisit.)
- **Splits:** group-disjoint by `(seq,charge)`; 3-way train/early-stop/**test** (per the 2026-06-16 audit — the current 2-way over-reports); report the gate metric on the untouched test fold; multi-seed for variance.

## Validation gates
- **Phase 0:** decoy-vs-entrapment audit (does the reverse-decoy score distribution match the true entrapment-false distribution?).
- **Held-out (offline):** target-matched-vs-decoy-matched AUC (the boundary that matters) + PSM-level `RichIonLLR`-sum target/decoy separation. NOT signal-vs-junk AUC.
- **Decisive (online):** A/B on Astral (high-res) + TMT/UPS (low-res) — net PSMs@1% threshold crossings + entrapment-FDP, **one variable at a time**, byte-identical-RawScore check before each, `--score rank` (additive) on low-res / `--score strong` on high-res.
- **Ship condition:** beats the RankScore baseline at honest FDP on at least one regime, neutral elsewhere; additive feature never regresses at fixed FDP.

## File structure
- `crates/scoring/src/ion_features.rs` **(new):** `extract_annotated_ion_features(peptide, scored_spec, ion, ctx) -> [f32; N]` — assignment-aware (complement + chemistry); the single source for train + serve (mirrors `frag_features.rs`).
- `crates/model-train/src/gbdt/ion_dataset.rs` **(new):** decoy-aware training-data builder (positives from gold flats, negatives from decoy peptides).
- `crates/scoring/src/scoring/strong_score.rs`: `rich_ion_llr(...)` summation + additive emission (0.0 when no model).
- `crates/search/src/match_engine.rs` + `crates/search/src/psm.rs` + `crates/output/src/pin.rs`: emit the `RichIonLLR` additive PIN column.
- `crates/model-train/src/store/*` + `crates/scoring/src/param_model.rs`: new GBDT blob column (nullable, back-compat).

## Phasing
- **P1 (this spec):** complement + chemistry features → decoy-aware train → held-out target/decoy AUC + crossings A/B.
- **P2:** explain-away (candidate-set training data) + isotope/loss/neighbor support; PSM-level aggregates (complement intensity correlation, runner-up collapse) as PIN columns.
- **P3:** make the summed LLR the RawScore engine if P1/P2 beat RankScore; per-regime deployment (rank+additive low-res, strong high-res).
