# Andes model-improvement roadmap (5-agent deep research synthesis)

**Date:** 2026-06-21
**Source:** 5 parallel read-only research agents — harvesting, cross-species, non-tryptic/HLA, training-methodology, model-architecture.
**Goal:** the final fully-own andes model set that BEATS the reference engine/a comparison search engine on confident PSMs at 1% true entrapment-FDP across regimes.

---

## ★★ THE HEADLINE FINDING (architecture agent)

**The shipped `resources/models.parquet` is the PRE-GBDT schema — the entire discriminative GBDT layer is DORMANT, emitting zeros in production.** There is no `frag_intensity_model_bytes` / `rich_ion_model_bytes` / `gbdt_model_bytes` column, so on every shipped model: `intensity_signal` falls back to a coarse table or 0, `RichIonLLR = 0`, the frag-LLR battery = (0,0,0). **`RawScore`/`StrongScore` and the exact features engineered to escape the "cosine ceiling" are emitting zeros.** The discriminative architecture is fully implemented and train/serve-correct — but no shipped model activates it. This is almost certainly the single biggest PSM lever, AND it gates honest evaluation of everything else (you can't A/B levers on zeros).

Two prior leads are **already fixed** (don't re-spend): charge-1 fragment blind spot (now `1..=2`, commit caa595e9); chimeric `max_n` hardcode (now data-driven default 4).

---

## Cross-cutting convergences (where independent agents agreed)

1. **Label quality / FDR** — *harvesting + training.* Flats are labeled by the reference engine `expect ≤ 0.01` **per-PSM, not a dataset-wide TDC q-value** → a false-label tail taught to the estimator as truth. Training agent independently found train_fdr=0.01 over heterogeneous corpora injects flat label-noise into abundant partitions. **Same root issue.**
2. **Noise model / over-training (the Astral −32%)** — *training + harvesting.* The default reversed-peptide noise is sampled at *signal/peak density*, so pooling heterogeneous spectra inflates `noise_freq` at the low ranks where real ions live → compresses every `ln(ion/noise)`. `ANDES_DENSE_NOISE` (density-decoupled) exists but is **off by default**. This is THE regression lever and the **precondition** for safely adding more data.
3. **Group / subset FDR** — *HLA + architecture.* Near-zero-risk lever; the PIN already carries charge one-hots, `peplen`, enzN/enzC, mod-class groups. Converts discriminative separation into recovered IDs in low-coverage subsets (the "Separation" half of IDs ≈ Coverage × Separation).
4. **Beat-the-teacher ceiling** — *harvesting (strategic).* Training on the reference engine's labels caps the ceiling at "what the reference engine finds" + bakes in its blind spots → at odds with *beating* it, and an independence concern. Long-term fix: self-label with andes once it leads (the in-repo `bootstrap_labels` TDC path).
5. **Species-agnostic = free coverage** — *cross-species.* The model keys only on fragmentation physics (zero organism/taxon in code); non-human spectra are valid training signal within a matched (activation,instrument,mod) slug — and already in production (E. coli/Candida/yeast in hcd_qexactive_tryp/uvpd/phospho).
6. **Geometry is seed-derived** — *architecture + the E3 independence item.* All 39 models = `num_segments=2`, `max_rank=150` (inherited, never re-derived). Independence + resolution lever.

---

## Unified prioritized roadmap

### Tier 0 — Foundation (do WITH the final-model retrain; they define the recipe)
- **F1. Dataset-wide TDC q≤0.01 labels in the flats** (harvesting #1). One edit to `mzml_pepxml_to_flat.py`; retroactively cleans all 18 slugs. Attacks the false-label tail + FDP honesty.
- **F2. Density-decoupled noise default for high-res training** (training #B1). Flag flip (`ANDES_DENSE_NOISE` / fixed-density reversed sampler). Fixes the Astral −32% and is the precondition for adding data.
- **F3. Per-slug corpus cap / subsampling** (harvesting #3, training #B5). Also fixes the `cid_lowres_tryp_tmt` GBDT OOM (498k PSMs). Learning-curve knee likely ~100–250k.

### Tier 1 — ★ Activate the discriminative layer (the big lever)
- **A1. Retrain + repack the models WITH the GBDT blob columns** (architecture #1). Turns RawScore/RichIonLLR/frag-LLR battery from zeros to live. Medium effort (Codon retrain + repack, no engine change). Validate: non-zero PIN columns + PSMs@1% entrapment-FDP vs shipped. **This is the first thing to measure.**

### Tier 2 — Near-zero-risk wins (parallel, cheap)
- **W1. Group/subset FDR** (HLA #3 + arch #2) — per charge / ntt-class / chimeric-secondary.
- **W2. Enable 1+ precursors + an `--hla` preset** (HLA #1+#5) — `charge_min=2` silently drops 1+ HLA peptides; pure recall.
- **W3. Cheap additive PIN features** (arch #3) — `ln(numCandidates)`, `log(MS2IonCurrent)`, all-ladder ppm-stdev (one-at-a-time, parity discipline).

### Tier 3 — Coverage expansion (SAFE only after F2 + A1)
- **C1. Cross-species pooling** (cross-species #1) — E. coli (PXD018176, 242k) + yeast flats *already on Codon* into hcd_qexactive_tryp; weight-swept; A/B Astral. Fills high-charge/high-mass tails.
- **C2. timsTOF curation** (harvesting #2) — the one data-rich un-curated gap; add cid_tof/hcd_tof slugs + a **timsTOF-HLA** slug. (Do NOT chase ETD/UVPD bottom-up — genuinely scarce in PRIDE.)
- **C3. Semi-tryptic (`--ntt semi`) + group-FDR** (arch #4) — biggest free coverage expansion, must pair with W1.
- **C4. Length-normalized scoring for HLA** (HLA #2) — the additive rank LLR is length-biased toward longer peptides; fight it with `RawScore/peplen` + length prior.

### Tier 4 — Bigger / longer bets
- **B1. Dedicated HLA model + corpus** (HLA #6) — 1+/13–25mer/class-II corpus (Sarkizova/Abelin); highest HLA ceiling.
- **B2. Re-derive partition geometry from corpus** (arch #5 = E3) — own geometry (independence) + resolution; 2 vs 3 vs 4 segments sweep.
- **B3. Held-out likelihood + early-stop + seed averaging** (training #B4) — makes regressions detectable at train time.
- **B4. Self-labeling with andes** (harvesting #7) — removes the competitor ceiling once andes leads a regime.
- **B5. RT features + learned re-scorer** (arch #3/#6) — RT is the biggest missing feature category; learned re-scorer only after A1.

### Audits (cheap, do early)
- Verify `uvpd_qexactive_tryp` provenance (295k PSMs implausible for UVPD-bottom-up — likely mislabeled HCD).
- The assembler still relabels cid/hcd tables onto 14+2 slugs (not independently trained) — intersects the independence campaign.

---

## How this shapes the FINAL-MODEL retrain recipe

The final fully-own model retrain should bake in Tier 0 + Tier 1 from the start:
**F1 (q-value labels) + F2 (density-decoupled noise) + F3 (corpus cap) → train each data-ready slug → A1 (GBDT blobs packed in) → own geometry (B2/E3) → assemble best-own-per-slug → independence gate → benchmark.**
These aren't separate from the campaign — they ARE the recipe that makes the final model both *independent* and *field-beating*.

## Validation discipline (all agents agreed)
Every change A/B'd one variable at a time, on the existing VM benchmark (TMT a05058 / UPS1 / Astral) at **1% true entrapment-FDP** (not reported FDR), with experiment-hygiene provenance (binary commit + model SHA + flat SHAs). Gate keep/drop on Astral; confirm no UPS/TMT regression before banking.

## Single highest-leverage first move
**A1 — activate the dormant GBDT layer** (repack the 3 benchmark models with the blob columns, trained with F1+F2). It is the cheapest path from "architecture built but inert" to "live," and the prerequisite for honestly measuring every other lever.
