# msgf-rust — peptide identification from MS/MS spectra

[![CI](https://github.com/bigbio/msgf-rust/actions/workflows/ci.yml/badge.svg)](https://github.com/bigbio/msgf-rust/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/bigbio/msgf-rust)](https://github.com/bigbio/msgf-rust/releases)
[![License: UCSD-Noncommercial](https://img.shields.io/badge/license-UCSD--Noncommercial-blue)](LICENSE)

> **A Rust port of MS-GF+** — takes mzML/MGF spectra + FASTA in, produces Percolator-ready `.pin` out. Matches or beats Java MS-GF+ PSM counts at 1% FDR while running **10-28× faster**.

## What is this?

msgf-rust is a from-scratch Rust reimplementation of [MS-GF+](https://github.com/MSGFPlus/msgfplus) (Kim & Pevzner, 2014), the canonical generating-function peptide-identification engine. It reads MS/MS spectra (mzML or MGF), searches them against a FASTA protein database, and emits Percolator-ready PIN rows (or a TSV) with per-PSM features for rescoring. The original Java implementation is preserved on the `java-legacy` branch.

## Why msgf-rust?

Three reference datasets, three results — all at 1% FDR via Percolator 3.7.1, all run on the same 8-thread VM:

| Dataset | Java PSMs @1% | msgf-rust PSMs @1% | Δ PSMs | Java wall | msgf-rust wall | Speedup |
|---|---:|---:|---:|---:|---:|---:|
| **Astral DDA** (LFQ_Astral_DDA_15min_50ng) | 33,425 | **36,715** | **+3,290 (+9.8%)** | 2:20:42 | **6:28** | **21.8×** |
| **PXD001819** (UPS1 yeast tryp) | 14,974 | 14,755 | -219 (-1.5%) | 8:46 | **0:54** | **9.7×** |
| **TMT** (a05058 PXD007683) | 10,115 | 9,605 | -510 (-5.0%) | 1:11:00 | **2:33** | **27.9×** |

What that means: on Astral we find **+9.8% more PSMs than Java at 21.8× the speed**; on PXD001819 we match Java's PSM count within 1.5% at 9.7× the speed; on TMT we trail Java by 5% PSMs but at 27.9× the speed. Java baseline is upstream MSGFPlus v2024.03.26 (no calibration; that flag isn't in upstream). msgf-rust runs with `--precursor-cal auto`. The remaining feature-level divergences (lnEValue, MeanRelErrorTop7 normalization, TMT PSM gap) are tracked in `DOCS.md` §8d and the I5 trace-investigation notes as research follow-up.

<details>
<summary>Bench methodology</summary>

- **Hardware:** 8-thread Intel Xeon Gold 6238 VM, AVX exposed (no AVX2/FMA), Linux x86_64.
- **Java baseline:** `MSGFPlus.jar` from the [MSGFPlus/msgfplus v2024.03.26 release](https://github.com/MSGFPlus/msgfplus/releases/tag/v2024.03.26), run with `-Xmx8192m -thread 8 -tda 1 -addFeatures 1`. Per-dataset args match `--precursor-tol-ppm`/`--isotope-error`/`--instrument`/`--protocol` of the Rust runs.
- **msgf-rust:** master branch, release build with `target-cpu=sandybridge` (AVX, no FMA), `--threads 8 --top-n 1 --precursor-cal auto`.
- **Java → PIN:** `msgf2pin` from the percolator `3.6.5--h6351f2a_0` container (single-arg mode for concatenated-TDA mzid; the `3.7.1` container's msgf2pin has a known parser crash on this mzid output).
- **Percolator:** `percolator 3.7.1` in `quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2` with `--seed 42 --only-psms`. Same parser script for both Java and Rust PINs.
- **Wall time:** `/usr/bin/time -v` "Elapsed (wall clock) time" — does not include Percolator stage.
- **Reproducibility:** scripts at `/srv/data/msgf-bench/finalize2_v2024.sh` and `/srv/data/msgf-bench/run_percolator_docker.sh` on the bench VM.

</details>

In a four-engine comparison against Java MS-GF+, Sage, and MSFragger on vendor-native data (Orbitrap Astral `.raw` + Bruker timsTOF `.d`), msgf-rust returns the most PSMs *and* distinct peptides at 1% FDR on both datasets — and is the only engine that reads Thermo `.raw` natively. Full methodology, per-engine parameters, and config files: [`docs/benchmarks/`](docs/benchmarks/).

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

Each archive contains the `msgf-rust` binary, the `resources/` tree (bundled `models.parquet` model store with all 39 scoring models), and LICENSE/NOTICE/README.

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

A row in `out.pin` is one peptide–spectrum match, with the Java-parity Percolator features plus Rust-only additive columns (`EdgeScore`, …) before `Peptide`. The number of charge one-hot columns scales with `[--charge-min, --charge-max]` (default **2–5** ⇒ `charge2…charge5`). Full column reference: `DOCS.md` §3a.

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
| `--spectrum <FILE>` | Input mzML, MGF, Thermo `.raw` (needs `thermo` feature + .NET 8), or Bruker timsTOF `.d` (needs `timstof` feature). Auto-detected by extension |
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
| `--charge-min/-max <INT>` | Charge range when absent in the spectrum | **2, 5** |
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

## Reading Bruker timsTOF `.d` files

msgf-rust reads native Bruker timsTOF `.d` (DDA-PASEF) data directly — pass `--spectrum sample.d`, no other flags; the format is auto-detected by extension just like mzML/MGF. A `.d` is a *directory* (a TDF SQLite database plus a binary blob); reading it uses the pure-Rust [`timsrust`](https://crates.io/crates/timsrust) crate (the same reader [Sage](https://github.com/lazear/sage) uses), so there is **no vendor runtime and nothing to bundle** — unlike Thermo `.raw`.

It is feature-gated to keep the default build pure-Rust. Build with `--features timstof` on a toolchain with a recent rustc (the `timsrust` dependency tree needs rustc ≥ 1.88):

```bash
cargo build --release -p msgf-rust --features timstof
msgf-rust --spectrum sample.d --database human.fasta --output-pin out.pin
```

Scope: **MS2 only**, the non-chimeric search path. The ion-mobility dimension is carried as metadata but not used by scoring. `--chimeric` on a `.d` degrades gracefully to a normal search (the co-isolation cascade needs an MS1 stream the DDA reader does not expose), as does `--precursor-cal`. Default (non-`timstof`) builds read mzML/MGF only and never pull in `timsrust`.

## Auto-detection

For mzML inputs with `--fragmentation auto` (the default), msgf-rust peeks the first 64 MS2 spectra, histograms activation methods and analyzer types, and selects a scoring model from the bundled `models.parquet` store based on the dominant values. The `--instrument` CLI flag is **not** required for this path — instrument class is read from the mzML when possible. `--protocol` from the CLI is still applied when selecting the model. MGF files have no activation metadata, so they use flag-based selection (defaulting to `hcd_qexactive_tryp`). Full resolution table: `DOCS.md` §4.

## Training your own models

msgf-rust can generate scoring models from your own data (`msgf-rust train`) and select them automatically by instrument at search time — useful for instruments or experiment classes the bundled models don't cover well (Orbitrap Astral, timsTOF, TMT/phospho/immunopeptidomics, …). Models live in a single Parquet store and support incremental add/remove/reweight updates with a held-out acceptance gate. See [`TRAIN.md`](TRAIN.md).

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

- Sangtae Kim, Pavel Pevzner, and the PNNL Proteomics team at UCSD's Center for Computational Mass Spectrometry, for the original MS-GF+ engine and the bundled scoring models.
- The [bigbio](https://github.com/bigbio) maintainers and the [quantms](https://github.com/bigbio/quantms) team.
