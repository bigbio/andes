# Heritage

andes began as an effort to improve **MS-GF+** (Kim & Pevzner, *Nat Commun* 2014;5:5277)
and has since become an independent project:

- **Scoring engine (code)** — clean-room reimplemented in Rust. The MS-GF+
  generating-function scoring was removed; ranking uses an empirical
  rank/intensity log-likelihood calibration.
- **Statistical models** — trained on our own curated mass-spectrometry data
  (`resources/ionstat/models.parquet`).
- **Implementation** — an independent Rust codebase.

We gratefully acknowledge MS-GF+ (originally developed at the University of
California, San Diego) as the origin and inspiration for this work, and we cite:

> Kim, S. and Pevzner, P.A. (2014). "MS-GF+ makes progress towards a universal
> database search tool for proteomics." *Nature Communications* 5:5277.

This repository's history begins at a clean root. The earlier
MS-GF+-derived development history is retained privately for provenance and is
not redistributed under this license.
