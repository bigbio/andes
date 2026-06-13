# andes documentation

This is the full reference for the `andes` binary and its outputs. For a quick start and benchmark summary, see [`README.md`](README.md).

Run `andes --help` for auto-generated help derived from the same `Cli` struct documented below.

---

## Contents

1. [CLI reference](#1-cli-reference)
2. [Mods.txt format](#2-modstxt-format)
3. [Output formats](#3-output-formats)
4. [Auto-detection](#4-auto-detection)
5. [Building from source](#5-building-from-source)
6. [Training new scoring models](#6-training-new-scoring-models)
7. [Isobaric labeling](#7-isobaric-labeling)
8. [Legacy numeric values & behavior notes](#8-legacy-numeric-values--behavior-notes)
9. [License and citation](#9-license-and-citation)

---

## 1. CLI reference

All flags use kebab-case long options (`--flag-name`). Several flags also accept legacy Java MS-GF+ numeric values (see §8). The CLI is implemented in `crates/andes/src/bin/andes.rs`.

### Input formats

`--spectrum` auto-detects the reader from the file extension — there is **no format flag** to set:

| Extension | Reader | Build requirement | Runtime requirement | Notes |
|---|---|---|---|---|
| `.mzML` / `.mzml` | mzML (streaming) | always built | none | Full activation + instrument auto-detection (§4). |
| `.raw` | Thermo RawFileReader | `--features thermo` (release archives ship it) | .NET 8 runtime — **bundled in the release archives** (nothing to install); from source, install .NET 8 | Native Thermo; parity-identical to the equivalent mzML. Supports `--chimeric`. Activation/instrument read from vendor metadata (§4). |
| `.d` | Bruker timsTOF (`timsrust`) | `--features timstof` | none (pure Rust) | DDA-PASEF, **MS2 only**; auto-routed to the `cid_tof_tryp` model. A `.d` is a *directory*. `--chimeric` / `--precursor-cal` degrade to a normal search. |
| any other (e.g. `.mgf`) | MGF | always built | none | No MS-level/activation metadata; treated as MS2 with flag-based model resolution. |

Native `.raw`/`.d` search **MS2 (identification) scans only** — MS1 and MS3+ scans (e.g. TMT SPS-MS3 reporter-quant) are filtered at load so `--ms-level 3` cannot accidentally search reporter scans. Default builds (no extra features) read mzML/MGF only; see [`README.md`](README.md) for `.raw`/`.d` install details and container recipes.

### Required

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--spectrum` | path | *(required)* | Input spectrum file. Reader auto-selected by extension — mzML, MGF, Thermo `.raw`, or Bruker timsTOF `.d` (see *Input formats* above). | Java `-s <FILE>` |
| `--database` | path | *(required)* | Target FASTA database. Decoys are generated automatically by reversing target sequences (see `--decoy-prefix`). | Java `-d <FILE>` |
| `--output-pin` | path | *(required)* | Output Percolator `.pin` file path. Always written unless the process exits with an error before the write phase. | Java `-o <FILE>` (when `-outputFormat pin`) |

### Search parameters

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--precursor-tol-ppm` | f64 | `20.0` | Symmetric precursor mass tolerance in parts per million. | Java `-t 20ppm` |
| `--charge-min` | u8 | `2` | Minimum precursor charge to try when the spectrum record does not specify charge. Must be ≤ `--charge-max` (inverted ranges are rejected at startup). | *(no direct Java flag; set via param file in Java)* |
| `--charge-max` | u8 | `5` | Maximum precursor charge to try when charge is missing from the spectrum. Must be ≥ `--charge-min`. The default range is **2–5**. | *(same)* |
| `--enzyme-specificity` | enum | `fully` | Enzymatic cleavage enforcement at peptide termini (Number of Tolerable Termini). `fully`: both termini must be cleavage sites (Java `-ntt 2`). `semi`: at least one terminus (Java `-ntt 1`). `non-specific`: neither required (Java `-ntt 0`). | `--ntt` alias; numeric `0`/`1`/`2` |
| `--max-missed-cleavages` | u32 | `1` | Maximum missed enzymatic cleavages allowed per candidate peptide. | Java `-maxMissedCleavages 1` |
| `--min-length` | u32 | `6` | Minimum peptide length in residues (excluding flanking context). | Java `-minLength 6` |
| `--max-length` | u32 | `40` | Maximum peptide length in residues. | Java `-maxLength 40` |
| `--top-n` | u32 | `10` | Maximum PSMs retained per spectrum (ranked by `RawScore`, best-first). | Java `-n 10` |
| `--isotope-error-min` | i8 | `-1` | Minimum isotope error offset to evaluate during precursor matching. Must be ≤ `--isotope-error-max`. | Java `-ti -1,2` (first value) |
| `--isotope-error-max` | i8 | `2` | Maximum isotope error offset to evaluate. Must be ≥ `--isotope-error-min`. | Java `-ti -1,2` (second value) |
| `--min-peaks` | u32 | `10` | Minimum number of MS2 peaks required to score a spectrum; spectra below this threshold are skipped. | Java `-minNumPeaks 10` |

### Modifications

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--mods` | path | *(off)* | Path to a Java-format `mods.txt` file describing fixed and variable modifications. When omitted, built-in defaults apply: Carbamidomethyl on C (fixed) and Oxidation on M (variable, max 3 per peptide). Composition strings (e.g. `C2H3N1O1`) are **not** supported — use numeric Da masses. | Java `-mod <FILE>` |
| | | | Hidden alias: `--mod` (singular). | |

### Scoring

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--fragmentation` | enum | *(none)* | **MGF-only.** Activation method for bundled model resolution: `CID`, `ETD`, `HCD`, `UVPD`. Auto-detected from file metadata for mzML/`.raw`/`.d`; has no effect on those formats. | Java `-m`; numeric `1`=CID, `2`=ETD, `3`=HCD, `4`=UVPD |
| `--fragment-tol-ppm` | float | *(none)* | **MGF-only.** Fragment-matching tolerance in ppm for high-resolution MS/MS (Orbitrap/TOF). No effect on mzML/`.raw`/`.d`. | *(no Java equivalent)* |
| `--fragment-tol-da` | float | *(none)* | **MGF-only.** Fragment-matching tolerance in Da for low-resolution MS/MS (ion trap). No effect on mzML/`.raw`/`.d`. | *(no Java equivalent)* |
| `--protocol` | enum | `auto` | Search protocol suffix for bundled model resolution. Named: `auto`, `phospho`, `iTRAQ`, `iTRAQ-phospho`, `TMT`, `standard`. | Java `-protocol`; numeric `0`=auto, `1`=phospho, `2`=iTRAQ, `3`=iTRAQ-phospho, `4`=TMT, `5`=standard |
| `--param-file` | path | *(auto)* | Explicit path to a `.param` scoring model file. When set, overrides all auto-detection and bundled resolution. Required when running a release binary outside the source tree if bundled resources are not present. | Java `-conf` / model path |
| `--model-store` | path | *(bundled)* | Path to a Parquet model store to use instead of the bundled `resources/ionstat/models.parquet`. Model selection reads from this store when set. | *(no Java equivalent)* |
| `--model` | string | *(auto-select)* | Exact model ID to load from the model store, skipping automatic selection by `(--fragmentation, --protocol)`. Useful for searching with a freshly-trained model (see `andes train`). | *(no Java equivalent)* |

**Bundled default:** mzML/`.raw`/`.d` inputs auto-detect activation and analyzer from file metadata (§4). For MGF inputs with no `--fragmentation` / `--fragment-tol-*` flags, andes defaults to `cid_lowres_tryp` (CID / low-res / 0.5 Da) and prints a warning.

**Model selection** (when `--param-file` is not set, resolved from `resources/ionstat/models.parquet`):

1. Build a selection key: `{Frag}_{Inst}_Trypsin` with optional protocol experiment class (e.g. `tmt`).
2. Exact match on the key → use that model.
3. If protocol-specific model absent, retry without the protocol class.
4. Final fallback: `cid_tof_tryp` (HCD + TOF/HighRes), `etd_lowres_tryp` (ETD), or `cid_lowres_tryp` (everything else).

**Normalisation rules** (mirrors Java `NewScorerFactory`):

- MGF with no `--fragmentation` flag → treated as `CID` / `low-res` for model selection.
- HCD + low-res analyzer (auto-detected or defaulted) → upgraded to `QExactive`.

Only tryptic enzyme models are in the store; other enzymes require `--param-file` with a binary `.param` file.

**Protocol values** (`--protocol`): selects a protocol-specific scoring model and, for isobaric labels, engages the dense-spectrum windowed peak filter. It does **not** add modifications — supply the label's fixed/variable mods in `mods.txt` (§2, §7).

| Value | Numeric | Model experiment class | Effect | Companion mods (in `mods.txt`) |
|---|---|---|---|---|
| `auto` *(default)* | `0` | Automatic | No protocol suffix; selects the generic model. Identical model selection to `standard`. | — |
| `phospho` | `1` | Phosphorylation | Prefers a phospho-specific model (e.g. `hcd_qexactive_tryp_phosphorylation`) when present; else falls back to the generic model. | Phospho variable mod (`+79.96633` on S/T/Y) |
| `iTRAQ` | `2` | iTRAQ | Prefers an iTRAQ-specific model; engages the dense-spectrum peak filter for crowded reporter-ion spectra. | iTRAQ fixed mods on K and peptide N-term |
| `iTRAQ-phospho` | `3` | iTRAQPhospho | Prefers a combined iTRAQ+phospho model; engages the dense-spectrum peak filter. | iTRAQ fixed mods + phospho variable mod |
| `TMT` | `4` | TMT | Prefers a TMT-specific model (e.g. `hcd_qexactive_tryp_tmt`); engages the dense-spectrum peak filter. | TMT fixed mods on K and peptide N-term |
| `standard` | `5` | Automatic | Explicit "no special protocol" — identical model selection to `auto`. | — |

When no protocol-specific model exists for the resolved `(fragmentation, instrument)` key, selection falls back to the generic model but the requested isobaric protocol is still stamped onto it so the dense-spectrum filter engages (e.g. `--protocol TMT` with no CID-TMT model resolves to `cid_lowres_tryp` with TMT filtering applied).

### Calibration

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--precursor-cal` | enum | `off` | Precursor-mass calibration: `off`, `auto`, or `on`. `auto`/`on` run a pre-pass that learns a systematic ppm shift from confident PSMs, then tighten the precursor tolerance for the main search; `auto` skips the correction when the sample is too small to be reliable. Opt-in only (default `off`). No effect on native `.raw` or `.d` input — calibration is not yet supported for those formats, so it is skipped (with a warning) and the search proceeds uncalibrated. | Java `-precursorCal auto\|on\|off` |

### Chimeric cascade

Opt-in two-pass search for co-isolated (co-fragmented) peptides. Requires an MS1 stream, so it runs on **mzML or Thermo `.raw`** only; on MGF/`.d` it warns and falls back to a normal search.

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--chimeric` | flag | *(off)* | Enable the two-pass chimeric cascade. Pass 1 is the normal top-1 search; Pass 2 detects co-isolated precursors in each scan's MS1 isolation window (averagine envelope match) and runs a targeted search for the second peptide on the *residual* spectrum (the primary's matched peaks removed), emitting it as an extra PSM. Forces top-1 per pass and always searches MS2 (`--ms-level` is ignored). Gains are entrapment-FDP validated. Experimental. | *(no Java equivalent)* |
| `--isolation-halfwidth` | f64 | `1.5` | Fallback isolation-window half-width in Da, used only when the mzML/`.raw` omits the per-scan isolation-window offsets. | *(no Java equivalent)* |

### Runtime

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--threads` | usize | logical CPU count | Rayon worker threads for the search loop. Pool is initialised once per process. | Java `-thread N` |
| `--ms-level` | u8 | `2` | MS level to search. Defaults to MS2 (identification); MS1 and MS3+ scans (e.g. TMT SPS-MS3 reporter-quant) are filtered at load so they never enter the search loop. Applies to mzML. Native `.raw`/`.d` always search MS2 regardless of this flag (a warning is printed if overridden), as does the chimeric cascade. MGF has no MS-level metadata and is always MS2. | *(no Java equivalent)* |
| `--max-spectra` | usize | `0` | Bench mode: process only the first N MS2 spectra. `0` = full input. When > 0, TSV output is skipped (PIN is still written). | *(no Java equivalent)* |
| `--decoy-prefix` | string | `XXX_` | Prefix prepended to reversed decoy protein accessions during index construction. | Java decoy tag in `-tda` workflows |

### Output

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--output-tsv` | path | *(off)* | Optional tab-separated PSM report (§3b). Skipped in bench mode (`--max-spectra > 0`). | Java `-outputFormat 1` with output path |

**Environment variable:** set `MSGF_RSS_PROBE=1` on Linux to print `VmRSS` checkpoints to stderr during long runs (debugging memory use). The legacy name `MSGFRUST_RSS_PROBE=1` is still accepted with a one-line deprecation warning and will be removed in the next quality cleanup.

---

## 2. Mods.txt format

andes reads the same modification file format as Java MS-GF+. The parser lives in `crates/model/src/modification.rs` and `crates/model/src/aa_set.rs`.

### Grammar

Each non-comment line is five comma-separated fields:

```text
<mass>,<aa>,<fix|opt>,<location>,<name>
```

| Field | Rule |
|---|---|
| `<mass>` | Numeric monoisotopic mass delta in Da. Composition strings (`C2H3N1O1`) are **not** supported in andes. |
| `<aa>` | Single uppercase ASCII letter, or `*` (wildcard). Multi-residue strings like `STY` are **not** supported — declare one line per residue. |
| `<fix\|opt>` | `fix` = fixed (static) modification; `opt` = variable modification. Case-insensitive. |
| `<location>` | One of `any`, `N-term`, `C-term`, `Prot-N-term`, `Prot-C-term` (case-insensitive; hyphens optional). |
| `<name>` | Human-readable modification name (used in logs; not written to mzIdentML — that format is not supported). |

**Special directive:** a line `NumMods=N` sets the maximum number of variable modifications per peptide. Parsed separately and applied to `SearchParams.max_variable_mods_per_peptide`. Default when absent: `3`.

**Comments:** lines whose first non-whitespace character is `#` are ignored. Inline `# ...` comments are stripped from the end of a line (Java `stripComment` semantics). Blank lines are ignored.

**Conflicts:** a fixed and variable mod targeting the same `(residue, location)` slot is rejected at build time.

### Example (a) — Carbamidomethyl C + Oxidation M

```text
NumMods=3
57.02146,C,fix,any,Carbamidomethyl
15.99491,M,opt,any,Oxidation
```

When `--mods` is omitted, andes uses these two modifications as built-in defaults.

### Example (b) — TMT 10-plex on K and peptide N-term

```text
NumMods=2
57.02146,C,fix,any,Carbamidomethyl
229.162932,K,fix,any,TMT10plex
229.162932,*,fix,N-term,TMT10plex
```

Pair with `--protocol TMT` to select the `hcd_qexactive_tryp_tmt` model from the store (activation and analyzer are auto-detected from mzML/`.raw`; see §4, §7).

### Example (c) — Phosphorylation on S, T, Y

```text
NumMods=3
57.02146,C,fix,any,Carbamidomethyl
79.966331,S,opt,any,Phospho
79.966331,T,opt,any,Phospho
79.966331,Y,opt,any,Phospho
```

Pair with `--protocol phospho` to prefer a phosphorylation-specific model (e.g. `hcd_qexactive_tryp_phosphorylation`) from the store when one is available.

---

## 3. Output formats

andes writes Percolator `.pin` (always) and optionally `.tsv`. Implementation: `crates/output/src/pin.rs`, `crates/output/src/tsv.rs`.

### 3a. PIN columns

Tab-separated, one header row, one row per PSM. Rows are sorted best-first within each spectrum by `RawScore` (the GF-free rank score) — the generating function and all of its derived score columns have been removed. With the default charge range `--charge-min 2 --charge-max 5`, the header has **39 columns** (listed below in order). The `chargeN` one-hots track the `--charge-min`…`--charge-max` range: one column per charge state, so narrowing/widening the range removes/adds one `chargeN` column each (e.g. a 2–3 range yields just `charge2 charge3`).

`--chimeric` does **not** change the column set; the header is the same 39 columns. What it changes: it populates the `PrecursorIsotopeKL` / `PrecursorSNR` features (`0.0` otherwise) from a linked MS1, and — because a single scan can then emit several rows — the `SpecId` of multi-row scans gains a trailing per-row index (see `SpecId` below). Multi-row scans also occur without `--chimeric` whenever rank-1 candidates tie on `RawScore`.

Columns in emission order (39 with the default 2–5 charge range):

| # | Column | Type | Description |
|---|---|---|---|
| 1 | `SpecId` | string | `{specID}_{scan}_{rank}` PSM identifier. On scans that emit more than one row (under `--chimeric`, or when candidates tie on `RawScore`), a per-row index is appended → `{specID}_{scan}_{rank}_{rowIdx}` to keep SpecIds unique. |
| 2 | `Label` | int | `+1` target, `-1` decoy (by **source protein**, not peptide sequence). |
| 3 | `ScanNr` | int | MS2 scan number from the input file. |
| 4 | `ExpMass` | float | Experimental neutral precursor mass (Da): `precursor_mz × charge − charge × proton`. |
| 5 | `CalcMass` | float | Theoretical neutral peptide mass (includes H₂O). |
| 6 | `mass` | float | Duplicate of `ExpMass` (OpenMS PercolatorAdapter convention). |
| 7 | `RawScore` | int | Rounded MS-GF+ score — the sole score/rank column. |
| 8 | `isotope_error` | int | Winning isotope offset (−1…2 by default). |
| 9 | `peplen` | int | Peptide residue count **+ 2** (includes flanking pre/post residues). |
| 10 | `dm` | float | Precursor mass error (Da) after isotope correction. |
| 11 | `absdm` | float | Absolute value of `dm`. |
| 12… | `charge2` … `chargeK` | 0/1 | One-hot encoding of assigned precursor charge; one column per state in the `--charge-min`…`--charge-max` range. |
| | `enzN` | 0/1 | N-terminal boundary consistent with enzyme rules. |
| | `enzC` | 0/1 | C-terminal boundary consistent with enzyme rules. |
| | `enzInt` | int | Count of internal enzymatic cleavage positions in the peptide. |
| | `NumMatchedMainIons` | int | Matched charge-1 b/y fragment positions. |
| | `longest_b` | int | Longest contiguous matched b-ion run. |
| | `longest_y` | int | Longest contiguous matched y-ion run. |
| | `longest_y_pct` | float | `longest_y / peptide.length()` (6 decimal places). |
| | `ExplainedIonCurrentRatio` | float | Matched b+y intensity / total MS2 intensity. |
| | `NTermIonCurrentRatio` | float | Matched b-ion intensity / total MS2 intensity. |
| | `CTermIonCurrentRatio` | float | Matched y-ion intensity / total MS2 intensity. |
| | `MS2IonCurrent` | float | Sum of all MS2 peak intensities (not log-scaled). |
| | `IsolationWindowEfficiency` | float | Always `0.0` (not available from parsed spectra). |
| | `MeanErrorTop7` | float | Mean absolute Da error of top-7 most-intense matched ions. |
| | `StdevErrorTop7` | float | Population stdev of absolute Da errors (top-7). |
| | `MeanRelErrorTop7` | float | Mean signed ppm error of top-7 ions. |
| | `StdevRelErrorTop7` | float | Population stdev of signed ppm errors (top-7). |
| | `matchedIonRatio` | float | `NumMatchedMainIons / peptide.length()`. |
| | `EdgeScore` | int | Per-bond DBScanScorer edge sum (IES + error score). Additive feature — not present in Java MS-GF+ PIN output. |
| | `PrecursorIsotopeKL` | float | KL divergence of the precursor isotope envelope vs the theoretical model (MS1). `0.0` unless `--chimeric` populates it from a linked MS1 scan. |
| | `PrecursorSNR` | float | Precursor signal-to-noise ratio from the MS1 envelope. `0.0` unless `--chimeric` populates it. |
| | `DeltaRawScore` | float | `RawScore(best) − RawScore(2nd-best distinct peptide)` for the spectrum; emitted on the rank-1 row only (`0.0` on lower-ranked rows). |
| | `TailorScore` | float | Tailor per-spectrum calibration (Yang et al., JPR 2020): `RawScore` ÷ the spectrum's top-1% RawScore quantile — makes RawScores comparable across spectra. |
| 38 | `Peptide` | string | `pre.SEQUENCE.post` with `+mass` mod annotations. |
| 39 | `Proteins` | string | Protein accession(s); decoy accessions carry `--decoy-prefix`. Multiple accessions tab-separated when one peptide maps to several proteins. |

### 3b. TSV columns

Tab-separated human-readable report. The `Title` column appears **only for MGF** inputs (Java parity).

**MGF header** (`is_mgf = true`):

| Column | Type | Description |
|---|---|---|
| `#SpecFile` | string | Bare filename of the input spectrum file. |
| `SpecID` | string | Spectrum identifier (MGF title, or `scan=N`). |
| `ScanNum` | int | Scan number. |
| `Title` | string | MGF `TITLE=` field. |
| `FragMethod` | string | Activation method name (`HCD`, `CID`, …) or `UNKNOWN`. |
| `Precursor` | float | Precursor m/z (4 decimal places). |
| `IsotopeError` | int | Winning isotope offset (same value as PIN `isotope_error`). |
| `PrecursorError(ppm)` | float | Mass error in ppm when tolerance is ppm mode; column named `PrecursorError(Da)` in Da mode. |
| `Charge` | int | Assigned precursor charge. |
| `Peptide` | string | Annotated peptide sequence with modifications. |
| `Protein` | string | Single protein accession (primary candidate). |
| `MSGFScore` | int | Rounded `RawScore` — the sole score column (the generating function and its derived score columns have been removed). |

**mzML header** — same as above **without** the `Title` column (11 columns total).

Decoy PSMs are included in TSV output; downstream tools label them via Percolator or manual filtering.

### 3c. PIN vs TSV — which to use

Use **PIN** when the goal is FDR calibration or rescoring: Percolator, MS²Rescore, Mokapot, and quantms-style pipelines consume `.pin` directly and learn feature weights from the full Percolator feature set (including `EdgeScore`). Use **TSV** for spreadsheet inspection, custom reporting, or tools that expect Java MS-GF+'s flat PSM table. You can emit both in one run with `--output-pin` and `--output-tsv`. For production quantms workflows, PIN is the standard path; TSV is optional diagnostics.

---

## 4. Auto-detection

For **mzML** inputs when `--fragmentation auto` (the default), andes peeks the input file before loading the full dataset:

1. **Activation method** — histogram of `<activation>` cvParams across the first 64 MS2 spectra; dominant method wins. Mixed methods trigger an stderr warning but the dominant method is still used file-wide.
2. **Instrument class** — scans `<instrumentConfiguration>` / analyzer cvParams via `input::detect_instrument_type`; dominant analyzer among MS2 spectra wins. `None` → `low-res` (Java `LOW_RESOLUTION_LTQ` default).

mzML/`.raw`/`.d` inputs are **zero-config**: andes reads the activation method and analyzer resolution directly from file metadata. `--fragmentation` and `--fragment-tol-*` are **MGF-only** parameters and have no effect on these formats. `--protocol` from the CLI is always applied to pick protocol-specific models from the parquet store (e.g. the `tmt` experiment-class entry).

MGF files carry no activation or analyzer metadata. Supply `--fragmentation <method>` (e.g. `HCD`, `CID`) and, for high-resolution MS/MS, `--fragment-tol-ppm <X>`, or for low-resolution MS/MS, `--fragment-tol-da <X>`. If none of these flags are provided for an MGF file, andes defaults to CID / low-res / 0.5 Da and prints a warning to stderr.

### Native Thermo `.raw`

A `.raw` file carries the activation method and analyzer in vendor metadata, so andes reads them directly (no mzML peek) and routes through the same parquet-store selection as mzML — e.g. beam-type CID (HCD) on an Orbitrap → `hcd_qexactive_tryp`. `--protocol` from the CLI still selects protocol-specific models (`tmt`, `itraq`). No fragmentation or instrument parameters are required or have any effect.

### Native Bruker timsTOF `.d`

timsTOF DDA-PASEF is beam-type CID on a TOF analyzer, so `.d` input auto-routes to the **`cid_tof_tryp`** model in the parquet store. `--protocol` still applies. Searched **MS2 only**; the ion-mobility dimension is carried as metadata but not used by scoring.

### Activation CV mapping (mzML `<activation>` cvParam accession → method)

| CV accession | Name (PSI-MS) | andes method | Notes |
|---|---|---|---|
| `MS:1000133` | collision-induced dissociation | CID | |
| `MS:1000422` | beam-type collision-induced dissociation (HCD) | HCD | |
| `MS:1000598` | electron transfer dissociation | ETD | |
| `MS:1000599` | pulsed Q dissociation | CID | Java collapses PQD → CID (`NewScorerFactory`) |
| `MS:1000435` | photodissociation | UVPD | Java UVPD mapping |
| `MS:1000250` | electron capture dissociation | ETD | Mapped to ETD (no dedicated ECD variant) |

### Instrument detection (analyzer cvParam → class)

| Analyzer family | Examples | Instrument class |
|---|---|---|
| Ion trap / linear ion trap | `MS:1000264`, Velos, LTQ | `low-res` |
| Orbitrap / Fusion | `MS:1000480`, Fusion Lumos | `QExactive` |
| FT-ICR | `MS:1000480` (FT) | `high-res` |
| TOF | `MS:1000128` | `TOF` |

### Bundled model store (`resources/ionstat/models.parquet`)

All 39 scoring models ship with the binary as a single Parquet model store
(`resources/ionstat/models.parquet`). The store covers the same
fragmentation × instrument × protocol matrix that the original 39 binary `.param`
files covered (CID/ETD/HCD/UVPD × LowRes/HighRes/TOF/QExactive × Trypsin, with
protocol variants for Phospho, TMT, iTRAQ, iTRAQPhospho). The individual binary
`.param` files are no longer shipped on disk; git history preserves them if the
store needs to be regenerated via `cargo run --example gen_bundled_store`.

**When auto-detection fails** (missing activation block, unknown CV term, or running outside the source tree without bundled resources): andes falls back to the `hcd_qexactive_tryp` model for default-flag runs, or to the resolution ladder in §1 for explicit flags. If no model resolves in the store, the process exits with an error instructing you to pass `--param-file <PATH>` with an external binary `.param` file.

---

## 5. Building from source

**Requirements:** Rust **1.85+** (workspace pins **1.87.0** in `rust-toolchain.toml` because transitive dependencies use `edition = "2024"`).

```bash
git clone https://github.com/bigbio/andes
cd andes
cargo build --release
# Binary: target/release/andes   (mzML + MGF; pure Rust)
```

**Native vendor formats** are feature-gated (the default build stays pure-Rust):

```bash
# Thermo .raw — needs rustc >= 1.88 and, at run time, the .NET 8 runtime
RUSTUP_TOOLCHAIN=stable cargo build --release -p andes --features thermo

# Bruker timsTOF .d — pure Rust, no vendor runtime
cargo build --release -p andes --features timstof

# Both at once (what the release archives ship for desktop/server targets)
RUSTUP_TOOLCHAIN=stable cargo build --release -p andes --features "thermo timstof"
```

See [`README.md`](README.md) (§Reading Thermo `.raw` / §Reading Bruker timsTOF `.d`) for the .NET 8 install, the bundled-runtime release archives, and container recipes.

Run the full workspace test suite:

```bash
cargo test --release --workspace
```

**CI-skipped tests:** GitHub Actions (`.github/workflows/ci.yml`) skips seven tests that fail on a clean checkout or are tracked as follow-up work. The release binary is unaffected.

| Skipped test | Reason |
|---|---|
| `charge_missing_spectrum_uses_per_charge_scored_spec` | `min_peaks` filter regression (pre-iter32 baseline) |
| `spectrum_without_charge_tries_charge_range` | same category |
| `known_peptide_appears_in_top_n` | same category |
| `read_bsa_canno_text_format` | Maven fixture under `target/test-classes/` not generated in CI |
| `read_tryp_pig_bov_revcat_csarr_cnlcp` | same |
| `tryp_pig_bov_revcat_full_set_loads` | same |
| `match_spectra_output_invariant_across_thread_counts` | Rayon tie-breaking nondeterminism when scores tie |

Reproduce the CI test invocation:

```bash
cargo test --release --workspace -- \
  --skip charge_missing_spectrum_uses_per_charge_scored_spec \
  --skip spectrum_without_charge_tries_charge_range \
  --skip known_peptide_appears_in_top_n \
  --skip read_bsa_canno_text_format \
  --skip read_tryp_pig_bov_revcat_csarr_cnlcp \
  --skip tryp_pig_bov_revcat_full_set_loads \
  --skip match_spectra_output_invariant_across_thread_counts
```

Release archives bundle the binary, the `models.parquet` model store (all 39 scoring models), and `unimod.obo` under `resources/` — see [`README.md`](README.md) §Install.

---

## 6. Training new scoring models

andes includes a native Rust training engine — **`andes train`** — that generates scoring models from your own data and writes them into the same Parquet model store the bundled models live in. No Java `ScoringParamGen` round-trip is needed.

Training is **bootstrap-supervised**: andes searches your spectra with a seed model, keeps the confident PSMs (target-decoy q ≤ `--train-fdr`), and re-estimates the per-partition rank and mass-error distributions from them. Trained models are auto-selected by instrument/protocol at search time, and the store supports incremental add / remove / reweight / decay updates with a held-out acceptance gate.

```bash
andes train \
  --spectra mydata.mzML \
  --database mydata.fasta \
  --seed-model hcd_qexactive_tryp \
  --out-store models.parquet \
  --model-id astral_tryp \
  --train-fdr 0.01
```

Then search with it:

```bash
andes --spectrum more.mzML --database mydata.fasta --output-pin out.pin \
  --model-store models.parquet --model astral_tryp
```

See **[`TRAIN.md`](TRAIN.md)** for the full guide: where to get training data, the experiment-class catalog, incremental training (`--update --add` / `--remove-source` / `--reweight` / `--decay`), and how to evaluate a candidate model on held-out data before committing it.

andes also still loads any Java MS-GF+ binary `.param` file directly via `--param-file` — the binary reader is retained for custom/external models and for the migration tooling. The 39 original models are bundled in `resources/ionstat/models.parquet`; the individual `.param` files are no longer shipped on disk (git history preserves them, and `cargo run -p model-train --example gen_bundled_store` regenerates the store from them).

---

## 7. Isobaric labeling

TMT and iTRAQ searches require **both** protocol-aware scoring models **and** correct fixed modifications in `mods.txt`. Set `--protocol TMT` or `--protocol iTRAQ` (or legacy `--protocol 4` / `--protocol 2`) so the model selector prefers protocol-specific models such as `hcd_qexactive_tryp_tmt` or `hcd_qexactive_tryp_itraq` from the bundled store.

### TMT (10-plex example)

**Mod masses:** TMT10plex = **229.162932 Da** on lysine and peptide N-terminus (Unimod). Carbamidomethyl on C is standard.

**mods.txt:**

```text
NumMods=2
57.02146,C,fix,any,Carbamidomethyl
229.162932,K,fix,any,TMT10plex
229.162932,*,fix,N-term,TMT10plex
```

**Command:**

```bash
andes \
  --spectrum tmt_spectra.mzML \
  --database hsapiens.fasta \
  --output-pin out.pin \
  --mods tmt_10plex_mods.txt \
  --protocol TMT
```

### iTRAQ (8-plex example)

**Mod masses:** iTRAQ8plex = **304.20536 Da** on K and peptide N-terminus.

**mods.txt:**

```text
NumMods=2
57.02146,C,fix,any,Carbamidomethyl
304.20536,K,fix,any,iTRAQ8plex
304.20536,*,fix,N-term,iTRAQ8plex
```

**Command:**

```bash
andes \
  --spectrum itraq_spectra.mzML \
  --database hsapiens.fasta \
  --output-pin out.pin \
  --mods itraq_8plex_mods.txt \
  --protocol iTRAQ
```

For phospho-enriched isobaric data use `--protocol iTRAQ-phospho` (legacy `--protocol 3`) and include phospho variable mods in `mods.txt` (§2 example c).

---

## 8. Legacy numeric values & behavior notes

For backward compatibility, the routing flags accept legacy 0…N numeric values in
addition to their canonical named values; clap parses named values
case-insensitively (`--fragmentation hcd` ≡ `HCD`).

| Flag | Numeric | Named | Notes |
|---|---|---|---|
| `--fragmentation` *(MGF-only)* | `1` | `CID` | |
| `--fragmentation` *(MGF-only)* | `2` | `ETD` | |
| `--fragmentation` *(MGF-only)* | `3` | `HCD` | |
| `--fragmentation` *(MGF-only)* | `4` | `UVPD` | |
| `--instrument` *(MGF-only, hidden)* | `0` | `low-res` | Legacy alias for `--fragment-tol-da 0.5` model key |
| `--instrument` *(MGF-only, hidden)* | `1` | `high-res` | Legacy alias for `--fragment-tol-ppm` model key |
| `--instrument` *(MGF-only, hidden)* | `2` | `TOF` | Legacy alias for TOF model key |
| `--instrument` *(MGF-only, hidden)* | `3` | `QExactive` | Legacy alias for QExactive model key |
| `--protocol` | `0` | `auto` | |
| `--protocol` | `1` | `phospho` | |
| `--protocol` | `2` | `iTRAQ` | |
| `--protocol` | `3` | `iTRAQ-phospho` | |
| `--protocol` | `4` | `TMT` | |
| `--protocol` | `5` | `standard` | |
| `--enzyme-specificity` (alias `--ntt`) | `0` | `non-specific` | |
| `--enzyme-specificity` (alias `--ntt`) | `1` | `semi` | |
| `--enzyme-specificity` (alias `--ntt`) | `2` | `fully` | |

### Behavior notes

- **Spectrum inputs:** mzML, MGF, native Thermo `.raw` (`thermo` feature), and native
  Bruker timsTOF `.d` (`timstof` feature) — see §1 *Input formats*.
- **Identification output:** Percolator PIN (always) plus an optional TSV; no mzIdentML.
- **Decoys:** always auto-generated by reversing target sequences at search time
  (prefix configurable via `--decoy-prefix`, default `XXX_`).
- **Enzyme:** Trypsin in the bundled models; other enzymes require a custom
  `--param-file`.
- **Modifications:** numeric Da masses only (composition strings are not parsed).
- **Memory:** spectra are processed in chunked streaming (5000/chunk), so large mzML
  files do not load fully into memory.

---

## 9. License and citation

andes is distributed under the **UCSD Noncommercial License** — the same terms as upstream MS-GF+. The license permits copying, modification, and distribution for **educational, research, and non-profit** purposes without fee, provided the copyright notice and liability paragraphs are retained. **Commercial use requires written permission** from the UCSD Technology Transfer Office (see `LICENSE` for contact details).

The software is provided **"as is"** without warranty. See [`LICENSE`](LICENSE) for the full upstream text and [`NOTICE`](NOTICE) for port attribution.

### Citation

If you use andes in published work, cite the original MS-GF+ paper:

> Kim, S. and Pevzner, P.A. (2014). MS-GF+ makes progress towards a universal database search tool for proteomics. *Nature Communications*, 5:5277.

And optionally this Rust port:

> bigbio (2026). andes: a Rust port of MS-GF+ for the quantms pipeline. https://github.com/bigbio/andes

The original Java implementation is preserved on the `java-legacy` branch; upstream MS-GF+ lives at https://github.com/MSGFPlus/msgfplus.
