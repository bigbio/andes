---
name: fdr-validation
description: Statistical rules for claiming that an andes change improves identifications or quantification - target-decoy, entrapment, seeds, held-out data. Use when a change alters scores, features, decoys, rescoring, MBR or quantification, or when reporting PSM/peptide/protein counts or DE results.
---

# FDR and statistical validation

More PSMs at a nominal 1% FDR is not an improvement unless the true error rate stays at or
below 1%. A feature that leaks decoy/target identity inflates counts while raising true FDP.

## Before claiming a gain

1. Search with an entrapment database (peptide-level shuffle: `docs/benchmarks/build_shuffled_entrap.py`)
   and report the true FDP next to the count (`docs/benchmarks/methods.md`).
2. Rescore with Percolator over at least three seeds; report the range, not the best seed.
3. Use at least two datasets of different kinds (for example high-res HCD and low-res CID)
   so a gain is not specific to one instrument.
4. For training or tuning, evaluate on data not used to fit or tune anything.

## Leakage checks for a new PIN feature

- Can the feature differ between targets and decoys for a reason other than match quality
  (decoy generation, protein-level information, peptide length bias, a database-derived
  prior)? If yes, it leaks.
- Compute the feature's target-vs-decoy distribution on known-wrong matches (entrapment
  hits): they should look like decoys.
- A feature that only helps after Percolator and raises entrapment FDP is rejected.

## Quantification

- Use designs with a known answer (spike-in ratios, mixed proteomes) and the shared
  evaluation in `docs/benchmarks/quant-eval/`.
- Report accuracy and precision together with the number of quantified proteins; report
  DE counts with their observed false-discovery proportion.
- Missing-value imputation and MBR transfers must be counted and reported separately.

## Reporting

Counts at the threshold, the true FDP range, seeds, datasets and commit. If another tool is
better on a metric, say so (the benchmark pages do).
