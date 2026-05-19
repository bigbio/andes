# Per-PSM Rust↔Java PIN diff harness

Stand-alone Python script that consumes a Java MS-GF+ `.pin` and an
`msgf-rust` `.pin` produced for the same dataset and emits an empirical
localization report:

  - Top-1-per-scan disagreement buckets (agreement / ranking flip / label
    flip / one-sided)
  - Per-feature distribution diff on the agreement bucket (where both
    engines pick the same target peptide for the same scan)
  - CSV of per-PSM diffs for further drill-down

The harness was introduced 2026-05-19 to break the
"piecewise-alignment-doesn't-work" trap: rather than guess which Java
behaviour to mirror next, the report shows empirically where Java and
Rust differ, ranked by magnitude.

## Usage

```bash
python3 benchmark/parity/analyze_rust_java_pin_diff.py \
    --java path/to/astral-java.pin \
    --rust path/to/astral-rust.pin \
    --out-dir path/to/output-dir
```

Outputs:

  - `<out-dir>/report.md` — human-readable bucket counts + per-feature
    diff table
  - `<out-dir>/per_psm_diff.csv` — one row per (scan, peptide) in the
    agreement bucket, with one column per numeric PIN feature carrying
    the Rust-minus-Java delta

## Conventions

- Both PIN files must be from the same dataset / FASTA / CLI flags.
- The harness reads PIN files BEFORE Percolator (raw scoring output);
  it does not consume the post-Percolator `*.target.psms.txt`.
- The peptide-residue key strips Percolator flanking and mod-mass tokens
  (mirroring `crates/search/tests/common/mod.rs::strip_flanking_and_mods`).
- Sentinel filter: rows where either side carries a `< -1e8` DeNovoScore
  placeholder (Rust's i32::MIN for GF-uncomputed PSMs) are excluded from
  per-feature stats to avoid skewing the mean.

## Dependencies

Python stdlib only. No pandas / numpy / scipy. Matches the project's
existing local-benchmark convention.
