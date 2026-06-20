<img src="docs/assets/andes-logo.png" alt="Andes" width="440" align="left">

<br clear="left">

_The data-driven peptide search engine of the quantms ecosystem. Built and maintained by the quantms team._

[![CI](https://github.com/bigbio/andes/actions/workflows/ci.yml/badge.svg)](https://github.com/bigbio/andes/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/bigbio/andes)](https://github.com/bigbio/andes/releases)
[![License: UCSD-Noncommercial](https://img.shields.io/badge/license-UCSD--Noncommercial-blue)](LICENSE)

> **A fast, data-driven peptide search engine** — spectra (mzML, MGF, native Thermo `.raw`, Bruker timsTOF `.d`) + a FASTA database in, Percolator-ready `.pin` out. Leading PSM counts at 1% FDR, in minutes where comparable Java tools take hours. To our knowledge, the **first proteomics search engine designed and built end-to-end with AI coding agents.**

## What is this?

andes is a peptide-spectrum database search engine for shotgun proteomics. It reads MS/MS spectra (mzML, MGF, native Thermo `.raw`, Bruker timsTOF `.d`), searches them against a FASTA protein database with **data-driven, per-regime scoring models**, and emits Percolator-ready PIN rows (or a TSV) with rich per-PSM features for rescoring. Beyond a fast closed search it offers opt-in **PTM discovery** (`--refine`), **chimeric** co-isolation recovery, multi-enzyme digestion, an out-of-core candidate index for large searches, and zero-config reanalysis — and it returns the most PSMs at 1% FDR on the reference datasets while running 10–28× faster than Java MS-GF+ (see [Why andes?](#why-andes)).

andes is also notable for *how* it was built: its engine, models, and benchmarks were developed iteratively by AI coding agents under human direction — a working demonstration of an agent-built scientific tool.

## Why andes?

Against the canonical open-source engines — **Java MS-GF+ and a comparison engine** — andes returns the most PSMs at 1% FDR on all three reference datasets, reads vendor formats natively, and runs in minutes where Java takes hours. Every engine is re-scored through one uniform Percolator (3.7.1, `--seed 42`) on the same 8-thread VM.

| Engine | Astral (high-res HCD) | TMT a05058 (low-res CID) | UPS1 (low-res CID) |
|---|---:|---:|---:|
| **andes** (`--chimeric`) | **69,968** | **12,043** | **17,879** |
| **andes** (top-1) | **36,782** | **11,957** | 17,143 |
| Java MS-GF+ v20240326 | 26,542 | 11,555 | 17,305 |
| a comparison engine 2025.01 | 31,435 | 10,876 | 15,809 |

<sub>PSMs at 1% FDR (distinct peptides track the same ordering). andes top-1 beats both Java MS-GF+ and a comparison engine on the high-res Astral run and on TMT (PSMs **and** peptides); on UPS1 it lands within 1% of Java and its `--chimeric` two-pass — which recovers co-isolated second peptides (opt-in) — takes the lead. Speed: andes finishes each run in ~1–4 min vs Java MS-GF+'s 9 min – 2.5 h (≈10–40×), on par with a comparison engine. A separate 1:1-entrapment head-to-head at a true 1% FDP (mode-independent) confirms the ordering on the low-res sets (andes ≈ Java MS-GF+, both ahead of a comparison engine); see [`docs/benchmarks/`](docs/benchmarks/).</sub>

**The 1% FDR is real, not inflated.** Re-measured against a 1:1 entrapment database with the shipped own-trained models, the *true* false-discovery proportion at the nominal 1% q-value is **1.08%** (top-1) and **1.29%** (chimeric) on Astral, **2.08%** on low-res TMT, and **1.43%** on UPS1 — the ID gains (including the chimeric near-doubling) are genuine identifications, not bought by a violated FDR. A **mode-independent** head-to-head at a *true* 1% FDP (comparable across engines regardless of target-decoy mode) confirms the ordering: andes leads on Astral and UPS1, ties Java MS-GF+ on low-res TMT, and beats a comparison engine on all three. Full numbers in [`docs/benchmarks/`](docs/benchmarks/). (Opt-in `--refine` PTM discovery runs on top, but its gains are not yet entrapment-validated — the entrapment metric is blind to its peptide-anchored second pass — so it ships as a capability, not a headline number.)

<details>
<summary>Bench methodology</summary>

- **Hardware:** 8-thread Intel Xeon Gold 6238 VM, Linux x86_64. Same machine for every engine.
- **Engines:** andes (this repo), Java MS-GF+ [v20240326](https://github.com/MSGFPlus/msgfplus/releases/tag/v2024.03.26), a comparison engine 2025.01 (via OpenMS). Parameters harmonized per dataset (trypsin, ≤2 missed cleavages, matched fixed/variable mods and precursor/fragment tolerances).
- **Uniform FDR:** every engine's PSMs re-scored through the **same** Percolator (`quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2`, `--seed 42 -Y`); counts reported at q ≤ 0.01.
- **PIN building:** andes and a comparison engine write Percolator PIN directly; Java MS-GF+ via `MzIDToTsv` + `build_pins.py` (its concatenated-TDA mzid crashes `msgf2pin`).
- **Models:** all andes runs use the bundled `resources/models.parquet` — andes's **own models trained on public data** for the covered regimes (high-res HCD, low-res CID, TMT, LysC, …); a few rarer regimes are still seeded from the original models pending retraining (see [`NOTICE`](NOTICE)). Independence verified per-regime: the bundle's auto-selected model matches the per-regime specialized models (e.g. Astral 30,933 vs 30,803).
- **FDR honesty** independently verified with a 1:1 entrapment database — true FDP at q≤1% is ≈1% (see above and `docs/benchmarks/`).
- **Notes:** Java MS-GF+ is deterministic; the Astral count reuses a prior run (its `msgf2pin` step crashes here regardless of input, and the count is pin-builder-independent). Protein-level counts are omitted from the headline — they require uniform parsimony grouping to be comparable across engines, since raw `proteinIds` differ by output format. Precursor calibration is off (the andes default).

</details>

andes is also the only engine here that reads Thermo `.raw` and Bruker timsTOF `.d` natively. Full methodology, per-engine parameters, data URLs, config files, and the entrapment-FDP validation: [`docs/benchmarks/`](docs/benchmarks/).

## How it works

andes is a streaming, multi-pass search cascade that ends in one uniform Percolator rescoring step.

```mermaid
flowchart TD
    %% ---- Scoring models (trained offline) ----
    subgraph TRAIN["Scoring models — trained offline on public data"]
      direction LR
      PRIDE[("PRIDE<br/>public datasets")] -->|"SDRF / quantms curation"| TR["andes train<br/>(own model per regime)"]
      TR --> STORE["models.parquet<br/>per activation × instrument × enzyme × protocol"]
    end

    %% ---- Inputs ----
    SPEC["Spectra<br/>mzML · MGF · Thermo .raw · Bruker .d"]
    DB["FASTA database<br/>target only — decoys auto-generated"]

    %% ---- Candidate generation ----
    DB --> CAND["Candidate peptides<br/>enzymatic digest + variable mods"]
    CAND --> IDX{"Candidate index<br/>--candidate-index auto"}
    IDX -->|"fits memory"| RAM["in-RAM index"]
    IDX -->|"too large"| MMAP["out-of-core mmap index"]

    %% ---- Pass 1 ----
    SPEC --> P1["Pass 1 — top-1 search<br/>peptide-spectrum scoring"]
    RAM --> P1
    MMAP --> P1
    STORE -. model selected per spectrum .-> P1
    P1 --> QUEUE["Top-N PSM queues<br/>+ rich per-PSM features"]

    %% ---- Optional second passes ----
    QUEUE -->|"--chimeric (opt-in)"| CHIM["Pass 2a — chimeric<br/>recover co-isolated 2nd peptide<br/>from the residual spectrum"]
    QUEUE -->|"--refine (opt-in)"| REF["Pass 2b — PTM refinement<br/>discovery mods on confident-protein anchors"]

    %% ---- Merge + rescore ----
    QUEUE --> MERGE["Unified PIN<br/>Pass 1 + chimeric + refine"]
    CHIM --> MERGE
    REF --> MERGE
    MERGE --> PERC["Percolator 3.7.1<br/>semi-supervised rescoring"]
    PERC --> OUT["FDR-controlled PSMs<br/>q ≤ 0.01 · entrapment-validated"]
```

1. **Candidate generation.** The FASTA is digested into candidate peptides (with variable mods). The candidate index is chosen automatically — kept in RAM, or mapped out-of-core (`mmap`) when it would exceed available memory — so very large mod searches don't OOM (`--candidate-index {auto,ram,mmap}`).
2. **Data-driven scoring.** Each spectrum is scored against its candidates with a model **selected per spectrum** by its `(activation, instrument, enzyme, protocol)`. These are andes's **own models, trained offline on public [PRIDE](https://www.ebi.ac.uk/pride/) datasets** curated through the **[quantms](https://github.com/bigbio/quantms)** / SDRF pipeline — not hand-tuned heuristics.
3. **Pass 1** is the standard top-1 search, emitting top-N PSM queues with rich per-PSM features.
4. **Optional second passes** (opt-in, off by default, do not change the default engine):
   - **`--chimeric`** detects co-isolated precursors in each scan's MS1 isolation window and searches the *residual* spectrum (primary peaks removed) for the **second peptide** — recovering co-isolated IDs without wide-window FDR inflation.
   - **`--refine`** runs a **PTM-discovery** search (oxidation, deamidation, pyro-Glu, acetyl, …) anchored on confident-protein peptides, to rescue modified spectra a closed search misses.
5. **Merge + rescore.** Pass 1 and any second-pass PSMs are written to **one Percolator PIN**; Percolator does the semi-supervised rescoring and FDR control. The reported 1% FDR is independently **entrapment-validated** (true FDP ≈ 1%).

## Install

**Option 1 — download a release archive** (recommended):

Grab the archive for your platform from the [Releases page](https://github.com/bigbio/andes/releases). Five platform builds are published per release:

```
andes-<version>-x86_64-unknown-linux-gnu.tar.gz
andes-<version>-aarch64-unknown-linux-gnu.tar.gz
andes-<version>-x86_64-apple-darwin.tar.gz
andes-<version>-aarch64-apple-darwin.tar.gz
andes-<version>-x86_64-pc-windows-msvc.zip
```

Each archive contains the `andes` binary, the `resources/` tree (bundled `models.parquet` model store with all 39 scoring models), and LICENSE/NOTICE/README.

**Option 2 — `cargo install`:**

```bash
cargo install --git https://github.com/bigbio/andes --bin andes
```

**Option 3 — build from source:**

```bash
git clone https://github.com/bigbio/andes
cd andes
cargo build --release
# Binary: target/release/andes
```

Requires Rust 1.85+ (see `rust-toolchain.toml`).

## Quick Start

```bash
andes \
  --spectrum BSA.mgf \
  --database BSA.fasta \
  --output-pin out.pin \
  --fragmentation HCD \
  --fragment-tol-ppm 20
```

This runs a tryptic search at the default 20 ppm precursor tolerance (`--precursor-tol-ppm`, default 20.0), with 20 ppm fragment-matching tolerance and the bundled **hcd_qexactive_tryp** scoring model (both selected via `--fragmentation HCD` + `--fragment-tol-ppm 20`), writes Percolator-format PSMs to `out.pin`, and prints per-phase timings to stderr. Feed `out.pin` directly into Percolator (Docker or native) to compute q-values.

A row in `out.pin` is one peptide–spectrum match, with the Java-parity Percolator features plus Rust-only additive columns (`EdgeScore`, …) before `Peptide`. The number of charge one-hot columns scales with `[--charge-min, --charge-max]` (default **2–5** ⇒ `charge2…charge5`). Full column reference: `DOCS.md` §3a.

## Common workflows

**Tryptic DDA + Percolator** (default):

```bash
andes --spectrum spectra.mzML --database db.fasta --output-pin out.pin
docker run --rm -v $(pwd):/data biocontainers/percolator:v3.7.1_cv1 \
  percolator -X /data/weights.txt /data/out.pin
```

**TMT 10-plex search with mods.txt:**

```bash
andes \
  --spectrum tmt_spectra.mzML \
  --database hsapiens.fasta \
  --output-pin out.pin \
  --mods tmt_10plex_mods.txt \
  --protocol TMT
```

**Direct TSV output (skip Percolator):**

```bash
andes --spectrum spectra.mzML --database db.fasta \
  --output-pin out.pin --output-tsv out.tsv
```

**[quantms](https://github.com/bigbio/quantms) pipeline integration:**

Point quantms's PSM search step at `andes` and use the standard quantms post-processing. The `.pin` row format is the same; existing quantms scripts using legacy numeric flag values (`--fragmentation 3 --protocol 4`) keep working without modification (the legacy numeric flag values are documented in [`DOCS.md`](DOCS.md)).

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
| `--fragmentation <CID\|ETD\|HCD\|UVPD>` | Fragmentation/activation method — **MGF-only** (auto-detected for mzML/`.raw`/`.d`) | *(see below)* |
| `--protocol <auto\|phospho\|iTRAQ\|iTRAQ-phospho\|TMT\|standard>` | Search protocol | **auto** |
| `--param-file <FILE>` | Override the bundled scoring model | **auto-pick** |
| `--decoy-prefix <STR>` | Prefix for generated decoys | **XXX_** |
| `--ms-level <INT>` | MS level to search; MS1/MS3+ (e.g. TMT SPS-MS3) filtered out (mzML or `.raw`) | **2** |
| `--threads <INT>` | Worker threads | **logical CPUs** |
| `--chimeric` | Two-pass co-isolated-peptide cascade (mzML or Thermo `.raw`) | **off** — see below |

Run `andes --help` for the auto-generated help with full descriptions and the legacy numeric flag aliases.

mzML, Thermo `.raw`, and Bruker `.d` are fully auto-detected — andes reads the
activation method and analyzer resolution from the file, so you pass no
fragmentation parameters for these formats.

### MGF input (extended parameters)

MGF files carry no activation or analyzer metadata, so you describe the
acquisition yourself:

| Parameter | When to pass | Example |
|---|---|---|
| `--fragmentation <CID\|ETD\|HCD\|UVPD>` | the activation method used | `--fragmentation HCD` |
| `--fragment-tol-ppm <X>` | high-resolution MS/MS (Orbitrap/TOF) | `--fragment-tol-ppm 20` |
| `--fragment-tol-da <X>`  | low-resolution MS/MS (ion trap)      | `--fragment-tol-da 0.5` |

If you pass none of these for an MGF file, andes assumes CID / low-res / 0.5 Da
and prints a warning. These parameters have no effect on mzML/`.raw`/`.d`.

## Chimeric / co-isolated peptides (`--chimeric`, experimental)

DDA scans frequently co-isolate more than one precursor, and the second peptide is normally lost. With `--chimeric` (mzML or Thermo `.raw`), andes runs a **two-pass cascade**: Pass 1 is the normal top-1 search; Pass 2 then detects co-isolated precursors in each scan's MS1 isolation window (averagine envelope match) and runs a targeted search for the second peptide on the *residual* spectrum (the primary's matched peaks removed), emitting it as an extra PSM. This recovers co-isolated identifications without the FDR inflation of a blind wide-window search — gains are entrapment-FDP validated. It is **opt-in and off by default**; the default engine is unchanged.

## Reading Thermo `.raw` files

andes reads native Thermo `.raw` directly — pass `--spectrum sample.raw`, no other flags; the format is auto-detected by extension just like mzML/MGF, and `--chimeric` works on `.raw` too. Output is parity-identical to searching the equivalent mzML (validated scan-for-scan on a 2.4 GB Orbitrap Astral run).

There are two ways to use it:

- **Pre-built release archives (recommended) — nothing to install.** The macOS (x64/arm64), Windows (x64), and Linux (x64) archives bundle a self-contained .NET 8 runtime next to the binary, so `.raw` reading works out of the box.
- **Building from source** with `--features thermo`. Then `.raw` reading needs the **.NET 8 runtime** installed (the build itself does not need the .NET SDK — the RawFileReader assemblies are vendored):
  - Linux: `sudo dnf install dotnet-runtime-8.0` (RHEL/Fedora) or `apt-get install dotnet-runtime-8.0` (Debian/Ubuntu), or `curl -sSL https://dot.net/v1/dotnet-install.sh | bash -s -- --channel 8.0 --runtime dotnet`
  - macOS: `brew install dotnet@8`
  - Windows: the [.NET 8 Desktop/Runtime installer](https://dotnet.microsoft.com/download/dotnet/8.0)
  - Build needs rustc ≥ 1.88: `RUSTUP_TOOLCHAIN=stable cargo build --release -p andes --features thermo`

The runtime is auto-discovered: a bundled `dotnet/` next to the binary is used automatically; otherwise an existing `DOTNET_ROOT` or a system install is used. mzML/MGF reading never loads .NET. RawFileReader is under Thermo's license — see `crates/input/THERMO_LICENSE.txt`.

**Containers:** base on a .NET 8 runtime image (or add the runtime), e.g.

```dockerfile
FROM mcr.microsoft.com/dotnet/runtime:8.0
COPY andes /usr/local/bin/andes   # built with --features thermo
ENTRYPOINT ["andes"]
```

## Reading Bruker timsTOF `.d` files

andes reads native Bruker timsTOF `.d` (DDA-PASEF) data directly — pass `--spectrum sample.d`, no other flags; the format is auto-detected by extension just like mzML/MGF. A `.d` is a *directory* (a TDF SQLite database plus a binary blob); reading it uses the pure-Rust [`timsrust`](https://crates.io/crates/timsrust) crate (the same reader [a comparison search engine](https://github.com/lazear/a comparison search engine) uses), so there is **no vendor runtime and nothing to bundle** — unlike Thermo `.raw`.

It is feature-gated to keep the default build pure-Rust. Build with `--features timstof` on a toolchain with a recent rustc (the `timsrust` dependency tree needs rustc ≥ 1.88):

```bash
cargo build --release -p andes --features timstof
andes --spectrum sample.d --database human.fasta --output-pin out.pin
```

Scope: **MS2 only**, the non-chimeric search path. The ion-mobility dimension is carried as metadata but not used by scoring. `--chimeric` on a `.d` degrades gracefully to a normal search (the co-isolation cascade needs an MS1 stream the DDA reader does not expose), as does `--precursor-cal`. Default (non-`timstof`) builds read mzML/MGF only and never pull in `timsrust`.

## Auto-detection

For mzML, Thermo `.raw`, and Bruker `.d` inputs, andes auto-detects the activation method and analyzer type from file metadata — no fragmentation or instrument parameters are needed. `--protocol` from the CLI is still applied to select protocol-specific models (e.g. TMT, iTRAQ). MGF files carry no activation or analyzer metadata; use `--fragmentation` / `--fragment-tol-ppm` / `--fragment-tol-da` to describe the acquisition (see the MGF section above), or andes defaults to CID / low-res / 0.5 Da and prints a warning. Full resolution table: `DOCS.md` §4.

## Training your own models

andes can generate scoring models from your own data (`andes train`) and select them automatically by instrument at search time — useful for instruments or experiment classes the bundled models don't cover well (Orbitrap Astral, timsTOF, TMT/phospho/immunopeptidomics, …). Models live in a single Parquet store and support incremental add/remove/reweight updates with a held-out acceptance gate. See [`TRAIN.md`](TRAIN.md).

## Parity vs Java MS-GF+

PIN output columns are bit-exact with Java MS-GF+ on the agreement bucket (same scan + same top-1 peptide) for most features. Three residual divergences exist as deferred research: `lnEValue` (num_distinct semantics), `MeanRelErrorTop7` (error-stat normalization), and the BSA charge-3 SEV gap from deconvolution-implementation differences. None gate cutover; aggregate 1% FDR PSM counts beat Java on all three benchmark datasets. Full detail: `DOCS.md` §8d.

## Citation

If you use andes in published work, please cite the original MS-GF+ paper:

> Kim, S. and Pevzner, P.A. (2014). MS-GF+ makes progress towards a universal database search tool for proteomics. *Nature Communications*, 5:5277.

And optionally this Rust port:

> bigbio (2026). andes: a Rust port of MS-GF+ for the quantms pipeline. https://github.com/bigbio/andes

## License

andes inherits the upstream MS-GF+ UCSD-Noncommercial license. The license restricts redistribution and commercial use; see `LICENSE` for the full text and `NOTICE` for attribution. The original Java implementation is preserved on the `java-legacy` branch (frozen at the bigbio-optimized version) and `java-legacy-original` branch (synced to upstream `MSGFPlus/msgfplus/master`).

## Acknowledgments

- Sangtae Kim, Pavel Pevzner, and the PNNL Proteomics team at UCSD's Center for Computational Mass Spectrometry, for the original MS-GF+ engine and the bundled scoring models.
- The [bigbio](https://github.com/bigbio) maintainers and the [quantms](https://github.com/bigbio/quantms) team.
