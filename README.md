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

Point quantms's PSM search step at `msgf-rust` and use the standard quantms post-processing. The `.pin` row format is the same; existing quantms scripts using legacy numeric flag values (`--fragmentation 3 --instrument 3 --protocol 4`) keep working without modification (see [`docs/CLI_MIGRATION.md`](docs/CLI_MIGRATION.md)).

## CLI summary

Most-used flags (full reference in `DOCS.md` §1):

Required:

| Flag | Purpose |
|---|---|
| `--spectrum <FILE>` | Input mzML, MGF, or Thermo `.raw` (auto-detected by extension; `.raw` needs the `thermo` feature + .NET 8) |
| `--database <FILE>` | Input FASTA (targets only; decoys generated) |
| `--output-pin <FILE>` | Percolator PIN output |

Optional (default in **bold**):

| Flag | Purpose | Default |
|---|---|---|
| `--output-tsv <FILE>` | Also write a TSV | **none** |
| `--mods <FILE>` | mods.txt file | **Cam-C fixed + Ox-M variable** |
| `--precursor-tol-ppm <FLOAT>` | Precursor mass tolerance (ppm) | **20.0** |
| `--precursor-cal <off\|auto\|on>` | Learn + apply a precursor ppm shift | **off** |
| `--isotope-error-min/-max <INT>` | Isotope-error range | **-1, 2** |
| `--charge-min/-max <INT>` | Charge range when absent in the spectrum | **2, 3** |
| `--enzyme-specificity <fully\|semi\|non-specific>` | Tolerable termini (NTT) | **fully** |
| `--max-missed-cleavages <INT>` | Missed cleavages | **1** |
| `--min-length/-max-length <INT>` | Peptide length range | **6, 40** |
| `--min-peaks <INT>` | Min peaks per spectrum to score | **10** |
| `--top-n <INT>` | PSMs retained per spectrum | **10** |
| `--fragmentation <auto\|CID\|ETD\|HCD\|UVPD>` | Fragmentation (auto-detected from mzML) | **auto** |
| `--instrument <low-res\|high-res\|TOF\|QExactive>` | Instrument class | **low-res** |
| `--protocol <auto\|phospho\|iTRAQ\|iTRAQ-phospho\|TMT\|standard>` | Search protocol | **auto** |
| `--param-file <FILE>` | Override the bundled scoring model | **auto-pick** |
| `--decoy-prefix <STR>` | Prefix for generated decoys | **XXX_** |
| `--ms-level <INT>` | MS level to search; MS1/MS3+ (e.g. TMT SPS-MS3) filtered out (mzML or `.raw`) | **2** |
| `--threads <INT>` | Worker threads | **logical CPUs** |
| `--chimeric` | Two-pass co-isolated-peptide cascade (mzML or Thermo `.raw`) | **off** — see below |

Run `msgf-rust --help` for the auto-generated help with full descriptions and the legacy numeric flag aliases.

## Chimeric / co-isolated peptides (`--chimeric`, experimental)

DDA scans frequently co-isolate more than one precursor, and the second peptide is normally lost. With `--chimeric` (mzML or Thermo `.raw`), msgf-rust runs a **two-pass cascade**: Pass 1 is the normal top-1 search; Pass 2 then detects co-isolated precursors in each scan's MS1 isolation window (averagine envelope match) and runs a targeted search for the second peptide on the *residual* spectrum (the primary's matched peaks removed), emitting it as an extra PSM. This recovers co-isolated identifications without the FDR inflation of a blind wide-window search — gains are entrapment-FDP validated. It is **opt-in and off by default**; the default engine is unchanged.

## Reading Thermo `.raw` files

msgf-rust reads native Thermo `.raw` directly — pass `--spectrum sample.raw`, no other flags; the format is auto-detected by extension just like mzML/MGF, and `--chimeric` works on `.raw` too. Output is parity-identical to searching the equivalent mzML (validated scan-for-scan on a 2.4 GB Orbitrap Astral run).

There are two ways to use it:

- **Pre-built release archives (recommended) — nothing to install.** The macOS (x64/arm64), Windows (x64), and Linux (x64) archives bundle a self-contained .NET 8 runtime next to the binary, so `.raw` reading works out of the box.
- **Building from source** with `--features thermo`. Then `.raw` reading needs the **.NET 8 runtime** installed (the build itself does not need the .NET SDK — the RawFileReader assemblies are vendored):
  - Linux: `sudo dnf install dotnet-runtime-8.0` (RHEL/Fedora) or `apt-get install dotnet-runtime-8.0` (Debian/Ubuntu), or `curl -sSL https://dot.net/v1/dotnet-install.sh | bash -s -- --channel 8.0 --runtime dotnet`
  - macOS: `brew install dotnet@8`
  - Windows: the [.NET 8 Desktop/Runtime installer](https://dotnet.microsoft.com/download/dotnet/8.0)
  - Build needs rustc ≥ 1.88: `RUSTUP_TOOLCHAIN=stable cargo build --release -p msgf-rust --features thermo`

The runtime is auto-discovered: a bundled `dotnet/` next to the binary is used automatically; otherwise an existing `DOTNET_ROOT` or a system install is used. mzML/MGF reading never loads .NET. RawFileReader is under Thermo's license — see `crates/input/THERMO_LICENSE.txt`.

**Containers:** base on a .NET 8 runtime image (or add the runtime), e.g.

```dockerfile
FROM mcr.microsoft.com/dotnet/runtime:8.0
COPY msgf-rust /usr/local/bin/msgf-rust   # built with --features thermo
ENTRYPOINT ["msgf-rust"]
```

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
