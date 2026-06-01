# Native Thermo `.raw` input — design

**Branch:** `feat/thermo-raw-input` (off `dev`)
**Status:** design (awaiting go-ahead to implement)

## Goal

`msgf-rust --spectrum sample.raw ...` reads Thermo `.raw` files natively on
Linux / macOS / Windows (and containers). **The only thing the user installs is
the .NET 8 runtime.** No new CLI flags — format is auto-detected by extension,
exactly like mzML/MGF today. mzML/MGF keep working with zero .NET involvement.

## Approach

A hand-rolled `crates/input/src/thermo.rs` backend built on the
[`thermorawfilereader`](https://crates.io/crates/thermorawfilereader) crate,
which wraps Thermo's **official RawFileReader** assemblies via a hosted .NET 8
runtime and exchanges spectra as FlatBuffers. This is the same engine `mzdata`
uses; we wrap it ourselves rather than adopt the whole `mzdata` umbrella.

Why this over a self-contained C# sidecar: far less code to own, uses the
official library exactly as the ecosystem already does, and meets the "install
.NET 8 and it just works" requirement. The trade-off — it needs the .NET 8
*runtime* present to open a `.raw` — is precisely the prerequisite the user
accepted.

### Dependency isolation (the key invariant)

- The `thermorawfilereader` dep sits behind a cargo feature **`thermo`**.
- **Build needs no .NET SDK**: `dotnetrawfilereader-sys` vendors the prebuilt C#
  assemblies + Thermo DLLs (~5 MB) and only compiles `netcorehost` glue.
- **The .NET runtime is loaded lazily** — only when a `ThermoRawReader` is
  constructed (i.e. only when a `.raw` is actually opened). So a binary built
  *with* `thermo`, run on a machine *without* .NET, still reads mzML/MGF fine;
  only `.raw` inputs surface a clear, actionable error.
- `thermo` is a **default feature** (released binaries include it), with
  `--no-default-features` available for a pure-Rust, no-vendor build.

Runtime discovery is automatic via hostfxr (`DOTNET_ROOT` / global install) —
no parameter. `DOTNET_RAWFILEREADER_BUNDLE_PATH` can override if ever needed.

## Components

1. **`crates/input/Cargo.toml`** — optional dep `thermorawfilereader`; feature
   `thermo = ["dep:thermorawfilereader"]`; `thermo` added to the workspace/binary
   default features.
2. **`crates/input/src/thermo.rs`** — `ThermoRawReader`:
   - Opens a `.raw` by filesystem path (RawFileReader only reads paths, not
     streams — acceptable; `--spectrum` is always a path).
   - Iterates scans → the existing `Spectrum` model: scan number, `ms_level`,
     retention time, precursor m/z + charge + isolation window, centroided
     peaks (m/z, intensity). Uses the lib's centroid/label stream for FT scans.
   - Mirrors the mzML reader's public surface: a plain spectrum iterator **and**
     an MS1-linked variant (`Ms1Link` + `read_with_ms1`-style) so `--chimeric`
     works on `.raw` (Milestone 2).
3. **`crates/input/src/lib.rs`** — re-export `ThermoRawReader` under
   `#[cfg(feature = "thermo")]`.
4. **`crates/msgf-rust/src/bin/msgf-rust.rs`** — extend extension dispatch
   (currently `is_mzml` vs MGF): add `.raw`/`.RAW` → `ThermoRawReader`.
   - feature on, .NET present → reads natively.
   - feature on, .NET missing → clear error: ".raw needs the .NET 8 runtime;
     install it from https://dotnet.microsoft.com/download. mzML/MGF are
     unaffected."
   - feature off (`--no-default-features`) → ".raw support not built in."

## Data-model parity

The produced `Spectrum` must be indistinguishable from the mzML reader's output
so the entire search/scoring path is unchanged. Correctness is proven by a
parity test (below): same `.raw`, converted to mzML via ThermoRawFileParser,
must yield matching scans/precursors/peaks.

## Testing

- **Unit:** FlatBuffers-spectrum → `Spectrum` mapping (precursor, isolation
  window, charge, peaks) against a small synthetic record.
- **Integration (`cfg(feature = "thermo")`):** open a tiny real `.raw` fixture
  and assert scan count / first MS2 precursor / peak count. Skips with a printed
  message when the .NET 8 runtime is absent (so CI without .NET still passes).
- **Parity:** assert `ThermoRawReader(sample.raw)` ≈ `MzMLReader(sample.mzML)`
  (same file converted) — scan-for-scan precursor + peak agreement within
  tolerance. The strongest correctness gate.

## Build / container / docs

- Released binary: `cargo build --release` (thermo default) — no .NET SDK.
- Runtime: install the **.NET 8 runtime**. Document the one-liner per OS
  (`apt install dotnet-runtime-8.0`, Homebrew `dotnet`, Windows installer).
- Containers: base on / add `mcr.microsoft.com/dotnet/runtime:8.0` (or apt the
  runtime). Document a minimal Dockerfile.
- README: a "Reading Thermo .raw" section + the Thermo RawFileReader license
  note (the DLLs ship inside the crate under Thermo's license; we reference it).

## Milestones

1. **MS2 read + dispatch + parity.** Wire crate/feature/dispatch; stream
   MS2 + precursor; search a `.raw` end-to-end; pass the mzML-parity test.
2. **MS1 streaming.** Expose the MS1-linked read so `--chimeric` works on `.raw`.
3. **Tests + fixture + container/CI docs.**

## Out of scope (separate follow-up branches)

- Bruker **timsTOF `.d`** (via `timsrust`, adds the ion-mobility dimension).
- **Sciex `.wiff`** (no open Rust reader; needs ProteoWizard / a converter).
- Profile-mode peak-picking beyond what RawFileReader's centroid stream gives.

## Licensing

Thermo's RawFileReader is under Thermo's license (not OSI). The DLLs are vendored
*inside* `dotnetrawfilereader-sys`, not committed to this repo; we add the
license reference to our docs/NOTICE.
