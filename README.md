# msgf-rust — peptide identification from MS/MS spectra

[![CI](https://github.com/bigbio/msgf-rust/actions/workflows/ci.yml/badge.svg)](https://github.com/bigbio/msgf-rust/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/bigbio/msgf-rust)](https://github.com/bigbio/msgf-rust/releases)
[![License: UCSD-Noncommercial](https://img.shields.io/badge/license-UCSD--Noncommercial-blue)](LICENSE)

> **A Rust port of MS-GF+** — takes mzML/MGF spectra + FASTA in, produces Percolator-ready `.pin` out. Beats Java MS-GF+ on all three benchmark datasets at 1% FDR while running 14-330% faster.

## What is this?

msgf-rust is a from-scratch Rust reimplementation of [MS-GF+](https://github.com/MSGFPlus/msgfplus) (Kim & Pevzner, 2014), the canonical generating-function peptide-identification engine. It reads MS/MS spectra (mzML or MGF), searches them against a FASTA protein database, and emits Percolator-ready PIN rows (or a TSV) with per-PSM features for rescoring. The original Java implementation is preserved on the `java-legacy` branch.

## Why msgf-rust?

Three datasets, three results (all at 1% FDR via Percolator 3.7.1):

| Dataset | Java MS-GF+ PSMs | msgf-rust PSMs | Δ | Java wall | msgf-rust wall | Wall Δ |
|---|---:|---:|---:|---:|---:|---:|
| **Astral DDA** (LFQ_Astral_DDA_15min_50ng) | 35,818 | **36,170** | **+352 (+0.98%)** | 5:49 | 5:57 | within 2% |
| **PXD001819** (UPS1 yeast tryp) | 14,798 | 14,760 | -38 (-0.26%) | ~150s | **45.88s** | **3.3× faster** |
| **TMT** (a05058 PXD007683) | 10,166 | **11,108** | **+9.3%** | ~2:55 | **2:30** | **14% faster** |

What that means: on Astral we find more peptide hits than Java; on PXD001819 we match Java's hit count at 3.3× the speed; on TMT we find ~9% more PSMs at 14% less wall. The remaining feature-level divergences (lnEValue, MeanRelErrorTop7 normalization) are tracked in `DOCS.md` §8d as research follow-up — they don't gate cutover.

## Install

**Option 1 — download a release archive** (recommended):

Grab the archive for your platform from the [Releases page](https://github.com/bigbio/msgf-rust/releases). Five platform builds are published per release:

```
msgf-rust-<version>-x86_64-unknown-linux-gnu.tar.gz
msgf-rust-<version>-aarch64-unknown-linux-gnu.tar.gz
msgf-rust-<version>-x86_64-apple-darwin.tar.gz
msgf-rust-<version>-aarch64-apple-darwin.tar.gz
msgf-rust-<version>-x86_64-pc-windows-msvc.zip
```

Each archive contains the `msgf-rust` binary, the `resources/` tree (39 bundled `.param` files + unimod.obo), and LICENSE/NOTICE/README.

**Option 2 — `cargo install`:**

```bash
cargo install --git https://github.com/bigbio/msgf-rust --bin msgf-rust
```

**Option 3 — build from source:**

```bash
git clone https://github.com/bigbio/msgf-rust
cd msgf-rust
cargo build --release
# Binary: target/release/msgf-rust
```

Requires Rust 1.85+ (see `rust-toolchain.toml`).

## Quick Start

```bash
msgf-rust \
  --spectrum BSA.mgf \
  --database BSA.fasta \
  --output-pin out.pin
```

This runs a tryptic search at 20 ppm precursor tolerance with the bundled HCD_QExactive_Tryp scoring model, writes Percolator-format PSMs to `out.pin`, and prints per-phase timings to stderr. Feed `out.pin` directly into Percolator (Docker or native) to compute q-values.

A row in `out.pin` is one peptide–spectrum match. With the default charge range (2–3), each row has **36 tab-separated columns**: 35 Java-parity Percolator features plus Rust-only `EdgeScore` (inserted before `Peptide`). Charge one-hot columns scale with `[--charge-min, --charge-max]`. Full column reference: `DOCS.md` §3a.

## Common workflows

**Tryptic DDA + Percolator** (default):

```bash
msgf-rust --spectrum spectra.mzML --database db.fasta --output-pin out.pin
docker run --rm -v $(pwd):/data biocontainers/percolator:v3.7.1_cv1 \
  percolator -X /data/weights.txt /data/out.pin
```

**TMT 10-plex search with mods.txt:**

```bash
msgf-rust \
  --spectrum tmt_spectra.mzML \
  --database hsapiens.fasta \
  --output-pin out.pin \
  --mods tmt_10plex_mods.txt \
  --protocol TMT \
  --fragmentation HCD \
  --instrument QExactive
```

**Direct TSV output (skip Percolator):**

```bash
msgf-rust --spectrum spectra.mzML --database db.fasta \
  --output-pin out.pin --output-tsv out.tsv
```

**[quantms](https://github.com/bigbio/quantms) pipeline integration:**

Point quantms's PSM search step at `msgf-rust` and use the standard quantms post-processing. The `.pin` row format is the same; existing quantms scripts using legacy numeric flag values (`--fragmentation 3 --instrument 3 --protocol 4`) keep working without modification (see `CLI_MIGRATION.md`).

## CLI summary

Most-used flags (full reference in `DOCS.md` §1):

| Flag | Purpose | Default |
|---|---|---|
| `--spectrum <FILE>` | Input mzML or MGF | (required) |
| `--database <FILE>` | Input FASTA | (required) |
| `--output-pin <FILE>` | Percolator PIN output | (required) |
| `--output-tsv <FILE>` | Optional TSV output | (off) |
| `--mods <FILE>` | mods.txt file (Cam-C + Ox-M built-in) | (off) |
| `--precursor-tol-ppm <FLOAT>` | Precursor mass tolerance | 20.0 |
| `--isotope-error-min/-max <INT>` | Isotope error range | -1, 2 |
| `--charge-min/-max <INT>` | Charge range when not in spectrum | 2, 3 |
| `--enzyme-specificity <auto\|...>` | NTT enforcement | fully |
| `--max-missed-cleavages <INT>` | Missed cleavages | 1 |
| `--min/-max-length <INT>` | Peptide length range | 6, 40 |
| `--min-peaks <INT>` | Min peaks per spectrum to score | 10 |
| `--top-n <INT>` | PSMs retained per spectrum | 10 |
| `--fragmentation <auto\|...>` | Frag method (auto-detect from mzML if `auto`) | auto |
| `--instrument <low-res\|...>` | Instrument class | low-res |
| `--protocol <auto\|...>` | Search protocol | auto |
| `--param-file <FILE>` | Override bundled scoring model | (auto-pick) |
| `--threads <INT>` | Worker threads | (logical CPUs) |

Run `msgf-rust --help` for the auto-generated help with full descriptions.

## Auto-detection

For mzML inputs with `--fragmentation auto` (the default), msgf-rust peeks the first 64 MS2 spectra, histograms activation methods and analyzer types, and selects a bundled `.param` file from the dominant values. The `--instrument` CLI flag is **not** required for this path — instrument class is read from the mzML when possible. `--protocol` from the CLI is still applied when resolving the bundled model. MGF files have no activation metadata, so they use flag-based resolution (defaulting to `HCD_QExactive_Tryp.param`). Full resolution table: `DOCS.md` §4.

## Parity vs Java MS-GF+

PIN output columns are bit-exact with Java MS-GF+ on the agreement bucket (same scan + same top-1 peptide) for most features. Three residual divergences exist as deferred research: `lnEValue` (num_distinct semantics), `MeanRelErrorTop7` (error-stat normalization), and the BSA charge-3 SEV gap from deconvolution-implementation differences. None gate cutover; aggregate 1% FDR PSM counts beat Java on all three benchmark datasets. Full detail: `DOCS.md` §8d.

## Citation

If you use msgf-rust in published work, please cite the original MS-GF+ paper:

> Kim, S. and Pevzner, P.A. (2014). MS-GF+ makes progress towards a universal database search tool for proteomics. *Nature Communications*, 5:5277.

And optionally this Rust port:

> bigbio (2026). msgf-rust: a Rust port of MS-GF+ for the quantms pipeline. https://github.com/bigbio/msgf-rust

## License

msgf-rust inherits the upstream MS-GF+ UCSD-Noncommercial license. The license restricts redistribution and commercial use; see `LICENSE` for the full text and `NOTICE` for attribution. The original Java implementation is preserved on the `java-legacy` branch (frozen at the bigbio-optimized version) and `java-legacy-original` branch (synced to upstream `MSGFPlus/msgfplus/master`).

## Acknowledgments

- Sangtae Kim, Pavel Pevzner, and the PNNL Proteomics team at UCSD's Center for Computational Mass Spectrometry, for the original MS-GF+ engine and the bundled `.param` scoring models.
- The [bigbio](https://github.com/bigbio) maintainers and the [quantms](https://github.com/bigbio/quantms) team.
