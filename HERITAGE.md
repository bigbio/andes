# Heritage

andes began as an effort to improve **MS-GF+** (Kim & Pevzner, *Nat Commun* 2014;5:5277)
and has since become an independent project:

- **Scoring engine (code)** — clean-room reimplemented in Rust. The MS-GF+
  generating-function scoring was removed; ranking uses an empirical
  rank/intensity log-likelihood calibration.
- **Statistical models** — trained on our own curated mass-spectrometry data
  (`resources/models.parquet`).
- **Implementation** — an independent Rust codebase.

We gratefully acknowledge MS-GF+ (originally developed at the University of
California, San Diego) as the origin and inspiration for this work.


