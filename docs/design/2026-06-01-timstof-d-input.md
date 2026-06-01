# Native Bruker timsTOF `.d` input — design

**Branch:** `feat/timstof-d-input` (off `dev`)
**Status:** implemented (draft PR)

## Goal

`msgf-rust --spectrum sample.d ...` reads native Bruker timsTOF `.d`
(DDA-PASEF) data directly — no other flags. The format is auto-detected by the
`.d` path extension, exactly like mzML/MGF/Thermo `.raw` today. mzML/MGF keep
working untouched, and the default (non-`timstof`) build never pulls in the new
dependency.

## Approach

A `crates/input/src/timstof.rs` backend built on the pure-Rust
[`timsrust`](https://crates.io/crates/timsrust) crate (v0.4.2) — **the same
reader Sage uses** for timsTOF. A Bruker `.d` is a DIRECTORY containing a TDF
SQLite database (`analysis.tdf`) plus a binary blob (`analysis.tdf_bin`);
`timsrust` reads both natively.

### Why `timsrust`

- **Pure Rust, no vendor runtime, no bundling.** Unlike the Thermo `.raw`
  reader (which needs the hosted .NET 8 runtime), `timsrust` compiles into the
  binary and reads a `.d` with zero external dependencies at runtime. This makes
  the feature strictly simpler than Thermo: no release-packaging or
  runtime-discovery work.
- **Sage-proven.** Sage reads timsTOF with this exact crate
  (`SpectrumReader::build().with_path(...).finalize()` → `get`/`get_all`), so
  the DDA grouping + precursor extraction are battle-tested on real data.

### Dependency isolation (the key invariant)

- The `timsrust` dep sits behind a cargo feature **`timstof`**, off by default
  at both the `input` crate level and the `msgf-rust` binary level.
- All `.d` code is behind `#[cfg(feature = "timstof")]`, so the default build is
  pure-Rust and reads mzML/MGF only. `timsrust` and its transitive tree (arrow,
  parquet, zstd, …) only compile under `--features timstof`.
- The lockfile change from declaring the optional dep is purely **additive** —
  no shared/default dependency is upgraded or downgraded — so the default
  `cargo check` MSRV is unchanged.

## Components

1. **`crates/input/Cargo.toml`** — optional dep `timsrust = { version = "0.4.2",
   optional = true }`; feature `timstof = ["dep:timsrust"]`.
2. **`crates/msgf-rust/Cargo.toml`** — feature `timstof = ["input/timstof"]`
   (non-default).
3. **`crates/input/src/timstof.rs`** — `TimsTofReader`:
   - Opens a `.d` directory by path via `SpectrumReader::new(path)`.
   - Iterates by index (`len`/`get`) → the existing `model::Spectrum`: 1-based
     scan number, precursor m/z + charge + intensity, retention time (carried in
     **seconds**, as `timsrust` reports it), isolation window (center + width →
     symmetric lower/upper offsets), centroided peaks `(m/z, intensity)`
     ascending by m/z.
   - `Option<Spectrum>` skip guard: a spectrum with no DDA precursor or a
     non-positive precursor m/z is skipped (mirrors the mzML/Thermo readers).
   - Implements `Iterator<Item = Result<Spectrum, TimsTofParseError>>`, exactly
     what the binary's `send_chunks` streaming consumer expects.
4. **`crates/input/src/lib.rs`** — `#[cfg(feature = "timstof")] pub mod timstof;`
   + re-export `TimsTofReader` / `TimsTofParseError`.
5. **`crates/msgf-rust/src/bin/msgf-rust.rs`** — extension dispatch adds
   `is_d` (`.d`) alongside `is_mzml`/`is_mgf`; `.d` is routed to `TimsTofReader`
   in the non-chimeric streaming reader, feature-gated, with a clear error when
   built without `--features timstof`. `is_mgf` now means "neither mzML nor
   `.d`". Precursor calibration is skipped for `.d` (a follow-up).

## Data-model parity

The produced `Spectrum` is the same struct the mzML/MGF/Thermo readers emit, so
the entire search/scoring path is format-agnostic. Mapping specifics:

- **Scan number / title:** `timsrust`'s 0-based `index` → 1-based `scan`; title
  is `scan=N` (scan-keyed like mzML/Thermo, so the PIN `SpecID`/`ScanNum`
  columns line up and the TSV writer uses its non-MGF, scan-based path).
- **RT:** `timsrust` reports precursor RT in **seconds**; the model stores
  seconds, so it is carried as-is. (Sage divides by 60 to get minutes for its
  own internal representation — we do not.)
- **Charge:** `Option<usize>` → `Option<i32>`; `Some(0)` is treated as unknown
  so the engine sweeps the configured charge range.
- **Isolation window:** `isolation_mz` (center) + `isolation_width` (total) →
  symmetric `isolation_lower_offset`/`isolation_upper_offset` = width / 2.
- **Activation:** left `None` — DDA-PASEF is beam-type CID/HCD-style but the
  `.d` records no discrete activation cvParam the param resolver keys on, so the
  resolver uses its default (or an explicit `--fragmentation`/`--instrument`).

## Scope

- **In scope:** DDA-PASEF, **MS2 only**, the non-chimeric search path.
  `timsrust`'s `SpectrumReader` already groups TIMS frames into centroided MS2
  fragment spectra each carrying their DDA precursor — exactly the search unit.
- **Out of scope (this PR):**
  - **Ion mobility.** The mobility dimension is extra metadata MS-GF+ scoring
    does not use; we carry RT/precursor as usual and ignore mobility for the
    base search. (Future idea: expose 1/K0 as an additive Percolator feature —
    not implemented here.)
  - **Chimeric / MS1-link on `.d`.** The co-isolation cascade needs an MS1
    stream, which the DDA `SpectrumReader` does not expose. `--chimeric` on `.d`
    degrades gracefully to a normal search (like MGF), with a warning.
  - **Precursor calibration on `.d`** (skipped + warned; a follow-up).

## Testing

- **Unit (`cfg(feature = "timstof")`, in-module):** `convert` precursor/charge/
  intensity/RT/isolation mapping, the precursor-skip and zero-charge guards, and
  `extract_peaks` ascending-sort — all against small synthetic `timsrust`
  records (no `.d` on disk needed).
- **Integration (`crates/input/tests/timstof_d_loads.rs`,
  `cfg(feature = "timstof")`):** opens a real `.d` and asserts MS2 invariants
  (positive precursor m/z, ascending peaks, scan numbers). A **no-op unless**
  `MSGF_TEST_D` points at a `.d` directory, so CI without a `.d` stays green.

## Benchmark dataset

**PXD072598** — HeLa, DDA-PASEF, timsTOF Pro 2, with FragPipe/MSFragger
reference results (`combined_peptide.tsv`, `fragger.params`). Smallest `.d`:
`HeLa_IAA_F51_1.d.zip` (1.11 GB) at
`ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/03/PXD072598/`. Human HeLa →
search against a human UniProt FASTA. (`.d` files are 1–3.5 GB; none is
downloaded on the dev machine.)

Converting one `.d` to mzML (e.g. via the Bruker→mzML path in ProteoWizard /
`tdf2mzml`) is an OPTIONAL aid to cross-check the native read against the mzML
reader on the same acquisition; the PRIMARY goal is reading the `.d` natively.

## Build / validation note

The dev machine has Homebrew rustc 1.87 (no rustup) and the repo pins `1.87.0`.
`timsrust`'s transitive tree likely needs rustc ≥ 1.88, so the `timstof` feature
**cannot be built locally**; this is expected. What was verified locally: the
default (non-`timstof`) `cargo check -p msgf-rust` stays green, and the lockfile
change is additive-only. The `timstof` build and a live `.d` read need
validation on a machine with rustc ≥ 1.88 and a `.d` file.

## Out of scope (separate follow-up branches)

- timsTOF **chimeric** (needs an MS1/frame-level accessor `timsrust` does expose
  via `FrameReader` — a larger effort than the DDA `SpectrumReader`).
- Ion-mobility as a Percolator feature.
- `.d` precursor calibration.

## Licensing

`timsrust` is Apache-2.0 (OSI). Nothing vendor-licensed is bundled.
