---
name: open-data-models
description: Rules for training, shipping and documenting andes models (scoring tables, GBDTs, andes-rt, andes-ms2, MBR weights) so they are built only from open data and stay independent and free of benchmark leakage. Use when training, porting or shipping any model or weights, editing resources/models or resources/predictors, or writing NOTICE/TRAIN.md text about models.
---

# Models from open data

Every model andes ships is trained by the andes authors on openly licensed public data,
with a reimplemented, published method. No third-party weights are redistributed.

## Training data

- Use only data that is public and licence-compatible with redistributing the derived
  weights under Apache-2.0 (public ProteomeXchange/PRIDE/MassIVE reanalyses, open
  benchmark sets). No data under embargo, no private or collaborator data, no data whose
  licence forbids derivative models.
- Never train on a dataset used to benchmark andes. Check the training corpus against the
  datasets under docs/benchmarks/ before any training run.
- Keep held-out or entrapment validation data separate from training data (TRAIN.md §7).

## Independence

- A shipped model must not be numerically identical to a third-party model (for example
  the MS-GF+ seed tables). Run `scripts/check_models_independence.py` against the seed
  before shipping a model store; a failed training run can silently keep seed tables.
- Third-party predictors (AlphaPeptDeep, DeepLC) are loaded from user-supplied files and
  never committed (`crates/onnx-infer/MODELS.md`).
- Reimplement methods from the paper's description; do not copy code or weights from
  repositories whose licence is not Apache-2.0 compatible.

## Recording provenance

- NOTICE says only that a model was trained by the andes authors on open data, its
  licence, and the paper whose method it follows. It does not list datasets.
- How it was trained (corpus description, exclusions, command lines, seeds, held-out
  result) goes in the commit message and PR description.
- The model store keeps its source ledger (TRAIN.md §5) so a source can be removed later.

## Shipping gate

Before a model replaces a bundled one: held-out yield at 1% FDR is at least the current
model's, entrapment FDP is at or below nominal (see `fdr-validation`), and the independence
check passes.
