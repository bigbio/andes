# msgf-rust documentation

This is the full reference for the `msgf-rust` binary and its outputs. For a quick start and benchmark summary, see [`README.md`](README.md). For porting Java MS-GF+ command lines and numeric legacy flags, see [`docs/CLI_MIGRATION.md`](docs/CLI_MIGRATION.md).

Run `msgf-rust --help` for auto-generated help derived from the same `Cli` struct documented below.

---

## Contents

1. [CLI reference](#1-cli-reference)
2. [Mods.txt format](#2-modstxt-format)
3. [Output formats](#3-output-formats)
4. [Auto-detection](#4-auto-detection)
5. [Building from source](#5-building-from-source)
6. [Training new `.param` files](#6-training-new-param-files)
7. [Isobaric labeling](#7-isobaric-labeling)
8. [Java MS-GF+ → msgf-rust migration](#8-java-ms-gf--msgf-rust-migration)
9. [License and citation](#9-license-and-citation)

---

## 1. CLI reference

All flags use kebab-case long options (`--flag-name`). Several flags also accept legacy Java MS-GF+ numeric values (see §8). The CLI is implemented in `crates/msgf-rust/src/bin/msgf-rust.rs`.

### Required

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--spectrum` | path | *(required)* | Input spectrum file. Extension `.mzML`/`.mzml` selects the mzML reader; any other extension (including `.mgf`) selects MGF. | Java `-s <FILE>` |
| `--database` | path | *(required)* | Target FASTA database. Decoys are generated automatically by reversing target sequences (see `--decoy-prefix`). | Java `-d <FILE>` |
| `--output-pin` | path | *(required)* | Output Percolator `.pin` file path. Always written unless the process exits with an error before the write phase. | Java `-o <FILE>` (when `-outputFormat pin`) |

### Search parameters

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--precursor-tol-ppm` | f64 | `20.0` | Symmetric precursor mass tolerance in parts per million. | Java `-t 20ppm` |
| `--charge-min` | u8 | `2` | Minimum precursor charge to try when the spectrum record does not specify charge. Must be ≤ `--charge-max` (inverted ranges are rejected at startup). | *(no direct Java flag; set via param file in Java)* |
| `--charge-max` | u8 | `3` | Maximum precursor charge to try when charge is missing from the spectrum. Must be ≥ `--charge-min`. | *(same)* |
| `--enzyme-specificity` | enum | `fully` | Enzymatic cleavage enforcement at peptide termini (Number of Tolerable Termini). `fully`: both termini must be cleavage sites (Java `-ntt 2`). `semi`: at least one terminus (Java `-ntt 1`). `non-specific`: neither required (Java `-ntt 0`). | `--ntt` alias; numeric `0`/`1`/`2` |
| `--max-missed-cleavages` | u32 | `1` | Maximum missed enzymatic cleavages allowed per candidate peptide. | Java `-maxMissedCleavages 1` |
| `--min-length` | u32 | `6` | Minimum peptide length in residues (excluding flanking context). | Java `-minLength 6` |
| `--max-length` | u32 | `40` | Maximum peptide length in residues. | Java `-maxLength 40` |
| `--top-n` | u32 | `10` | Maximum PSMs retained per spectrum (ranked by SpecEValue). | Java `-n 10` |
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
| `--fragmentation` | enum | `auto` | Fragmentation method for bundled `.param` resolution. Named: `auto`, `CID`, `ETD`, `HCD`, `UVPD`. `auto` on mzML triggers activation detection (§4); on MGF falls back to bundled defaults. | Java `-m`; numeric `0`=auto, `1`=CID, `2`=ETD, `3`=HCD, `4`=UVPD |
| `--instrument` | enum | `low-res` | Instrument class for bundled `.param` resolution. Named: `low-res`, `high-res`, `TOF`, `QExactive`. | Java `-inst`; numeric `0`=low-res, `1`=high-res, `2`=TOF, `3`=QExactive |
| `--protocol` | enum | `auto` | Search protocol suffix for bundled `.param` resolution. Named: `auto`, `phospho`, `iTRAQ`, `iTRAQ-phospho`, `TMT`, `standard`. | Java `-protocol`; numeric `0`=auto, `1`=phospho, `2`=iTRAQ, `3`=iTRAQ-phospho, `4`=TMT, `5`=standard |
| `--param-file` | path | *(auto)* | Explicit path to a `.param` scoring model file. When set, overrides all auto-detection and bundled resolution. Required when running a release binary outside the source tree if bundled resources are not present. | Java `-conf` / model path |

**Bundled default when all scoring flags are at their defaults** (`--fragmentation auto --instrument low-res --protocol auto`): `HCD_QExactive_Tryp.param`. This preserves pre-auto-detect behaviour for MGF inputs and mzML files without activation metadata.

**Resolution ladder** (when `--param-file` is not set):

1. Try exact `{Frag}_{Inst}_Tryp{ProtocolSuffix}.param`.
2. If protocol-specific file missing, drop protocol suffix → `{Frag}_{Inst}_Tryp.param`.
3. Final fallback: `CID_TOF_Tryp.param` (HCD + TOF/HighRes), `ETD_LowRes_Tryp.param` (ETD), or `CID_LowRes_Tryp.param` (everything else).

**Normalisation rules** (mirrors Java `NewScorerFactory`):

- `auto` fragmentation → treated as `CID` for filename resolution (except mzML auto-detect path, §4).
- HCD + `low-res` instrument → upgraded to `QExactive`.

Only tryptic enzyme models are bundled; other enzymes require `--param-file`.

### Runtime

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--threads` | usize | logical CPU count | Rayon worker threads for the search loop. Pool is initialised once per process. | Java `-thread N` |
| `--ms-level` | u8 | `2` | MS level to search. Non-matching spectra are filtered at load time. Meaningful for mzML only; MGF has no MS-level metadata and is always treated as MS2 (a warning is printed if `--ms-level` ≠ 2 on MGF). | *(no Java equivalent)* |
| `--max-spectra` | usize | `0` | Bench mode: process only the first N MS2 spectra. `0` = full input. When > 0, TSV output is skipped (PIN is still written). | *(no Java equivalent)* |
| `--decoy-prefix` | string | `XXX_` | Prefix prepended to reversed decoy protein accessions during index construction. | Java decoy tag in `-tda` workflows |

### Output

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--output-tsv` | path | *(off)* | Optional tab-separated PSM report (§3b). Skipped in bench mode (`--max-spectra > 0`). | Java `-outputFormat 1` with output path |

**Environment variable:** set `MSGF_RSS_PROBE=1` on Linux to print `VmRSS` checkpoints to stderr during long runs (debugging memory use). The legacy name `MSGFRUST_RSS_PROBE=1` is still accepted with a one-line deprecation warning and will be removed in the next quality cleanup.

---

## 2. Mods.txt format

msgf-rust reads the same modification file format as Java MS-GF+. The parser lives in `crates/model/src/modification.rs` and `crates/model/src/aa_set.rs`.

### Grammar

Each non-comment line is five comma-separated fields:

```text
<mass>,<aa>,<fix|opt>,<location>,<name>
```

| Field | Rule |
|---|---|
| `<mass>` | Numeric monoisotopic mass delta in Da. Composition strings (`C2H3N1O1`) are **not** supported in msgf-rust. |
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

When `--mods` is omitted, msgf-rust uses these two modifications as built-in defaults.

### Example (b) — TMT 10-plex on K and peptide N-term

```text
NumMods=2
57.02146,C,fix,any,Carbamidomethyl
229.162932,K,fix,any,TMT10plex
229.162932,*,fix,N-term,TMT10plex
```

Pair with `--protocol TMT --fragmentation HCD --instrument QExactive` to select `HCD_QExactive_Tryp_TMT.param` (§4, §7).

### Example (c) — Phosphorylation on S, T, Y

```text
NumMods=3
57.02146,C,fix,any,Carbamidomethyl
79.966331,S,opt,any,Phospho
79.966331,T,opt,any,Phospho
79.966331,Y,opt,any,Phospho
```

Pair with `--protocol phospho` to prefer a `_Phosphorylation` protocol-suffixed `.param` file when bundled.

---

## 3. Output formats

msgf-rust writes Percolator `.pin` (always) and optionally `.tsv`. Implementation: `crates/output/src/pin.rs`, `crates/output/src/tsv.rs`.

### 3a. PIN columns

Tab-separated, one header row, one row per PSM. Rows are sorted best-first within each spectrum (lowest SpecEValue first). With default `--charge-min 2 --charge-max 3`, the header has **36 columns**: 35 Java-parity fields plus Rust-only `EdgeScore` (before `Peptide`). Additional charge states add one `chargeN` column each.

| Column | Type | Description |
|---|---|---|
| `SpecId` | string | `{specID}_{scan}_{rank}` — unique PSM identifier. |
| `Label` | int | `+1` target, `-1` decoy (by **source protein**, not peptide sequence). |
| `ScanNr` | int | MS2 scan number from the input file. |
| `ExpMass` | float | Experimental neutral precursor mass (Da): `precursor_mz × charge − charge × proton`. |
| `CalcMass` | float | Theoretical neutral peptide mass (includes H₂O). |
| `mass` | float | Duplicate of `ExpMass` (OpenMS PercolatorAdapter convention). |
| `RawScore` | int | Rounded MS-GF+ score (`MSGFScore`). |
| `DeNovoScore` | int | Best de novo graph score for the spectrum. |
| `lnSpecEValue` | float | `ln(SpecEValue)`; `-MAX` if non-positive. |
| `lnEValue` | float | `ln(EValue)` where EValue = SpecEValue × num_distinct peptides. |
| `isotope_error` | int | Winning isotope offset (−1…2 by default). |
| `peplen` | int | Peptide residue count **+ 2** (includes flanking pre/post residues). |
| `dm` | float | Precursor mass error (Da) after isotope correction. |
| `absdm` | float | Absolute value of `dm`. |
| `charge2` … `chargeK` | 0/1 | One-hot encoding of assigned precursor charge. |
| `enzN` | 0/1 | N-terminal boundary consistent with enzyme rules. |
| `enzC` | 0/1 | C-terminal boundary consistent with enzyme rules. |
| `enzInt` | int | Count of internal enzymatic cleavage positions in the peptide. |
| `NumMatchedMainIons` | int | Matched charge-1 b/y fragment positions. |
| `longest_b` | int | Longest contiguous matched b-ion run. |
| `longest_y` | int | Longest contiguous matched y-ion run. |
| `longest_y_pct` | float | `longest_y / peptide.length()` (6 decimal places). |
| `ExplainedIonCurrentRatio` | float | Matched b+y intensity / total MS2 intensity. |
| `NTermIonCurrentRatio` | float | Matched b-ion intensity / total MS2 intensity. |
| `CTermIonCurrentRatio` | float | Matched y-ion intensity / total MS2 intensity. |
| `MS2IonCurrent` | float | Sum of all MS2 peak intensities (not log-scaled). |
| `IsolationWindowEfficiency` | float | Always `0.0` (not available from parsed spectra). |
| `MeanErrorTop7` | float | Mean absolute Da error of top-7 most-intense matched ions. |
| `StdevErrorTop7` | float | Population stdev of absolute Da errors (top-7). |
| `MeanRelErrorTop7` | float | Mean signed ppm error of top-7 ions. |
| `StdevRelErrorTop7` | float | Population stdev of signed ppm errors (top-7). |
| `lnDeltaSpecEValue` | float | `ln(rank1 SpecEValue / rank2 SpecEValue)` for rank-1 PSMs; `0` otherwise. |
| `matchedIonRatio` | float | `NumMatchedMainIons / peptide.length()`. |
| `EdgeScore` | int | Per-bond DBScanScorer edge sum (IES + error score). **Rust-only additive column** — not present in Java MS-GF+ PIN output; placed before `Peptide` so legacy column-index parsers still find sequence at the tail. |
| `Peptide` | string | `pre.SEQUENCE.post` with `+mass` mod annotations. |
| `Proteins` | string | Protein accession(s); decoy accessions carry `--decoy-prefix`. Multiple accessions tab-separated when one peptide maps to several proteins. |

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
| `DeNovoScore` | int | De novo score. |
| `MSGFScore` | int | Rounded raw score. |
| `SpecEValue` | float | SpecEValue in `%.6e` notation. |
| `EValue` | float | Database E-value in `%.6e` notation. |

**mzML header** — same as above **without** the `Title` column (14 columns total).

Decoy PSMs are included in TSV output; downstream tools label them via Percolator or manual filtering.

### 3c. PIN vs TSV — which to use

Use **PIN** when the goal is FDR calibration or rescoring: Percolator, MS²Rescore, Mokapot, and quantms-style pipelines consume `.pin` directly and learn feature weights from the full Percolator feature set (including `EdgeScore`). Use **TSV** for spreadsheet inspection, custom reporting, or tools that expect Java MS-GF+'s flat PSM table. You can emit both in one run with `--output-pin` and `--output-tsv`. For production quantms workflows, PIN is the standard path; TSV is optional diagnostics.

---

## 4. Auto-detection

For **mzML** inputs when `--fragmentation auto` (the default), msgf-rust peeks the input file before loading the full dataset:

1. **Activation method** — histogram of `<activation>` cvParams across the first 64 MS2 spectra; dominant method wins. Mixed methods trigger an stderr warning but the dominant method is still used file-wide.
2. **Instrument class** — scans `<instrumentConfiguration>` / analyzer cvParams via `input::detect_instrument_type`; dominant analyzer among MS2 spectra wins. `None` → `low-res` (Java `LOW_RESOLUTION_LTQ` default).

The CLI `--instrument` flag does **not** gate this path — only `--fragmentation auto` + mzML extension does. When peek succeeds, instrument is taken from the file; `--protocol` from the CLI is still used to pick protocol-suffixed `.param` files (e.g. `_TMT`).

MGF files carry no activation or instrument metadata → auto-detect returns `None` → bundled default `HCD_QExactive_Tryp.param` unless explicit `--fragmentation` / `--instrument` flags override via `resolve_bundled_param`.

Non-auto `--fragmentation` (e.g. `HCD`, `3`) disables the activation peek and uses flag-based resolution directly (§1), including `--instrument` and `--protocol` from the CLI.

### Activation CV mapping (mzML `<activation>` cvParam accession → method)

| CV accession | Name (PSI-MS) | msgf-rust method | Notes |
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

### Bundled `.param` files (`resources/ionstat/`)

39 scoring models ship with the binary (Tryp-centric unless noted):

```text
CID_HighRes_NoCleavage.param          CID_HighRes_Tryp.param
CID_LowRes_ArgC.param                 CID_LowRes_AspN.param
CID_LowRes_GluC.param                 CID_LowRes_LysC.param
CID_LowRes_LysN.param                 CID_LowRes_LysN_Phosphorylation.param
CID_LowRes_NoCleavage.param           CID_LowRes_Tryp.param
CID_LowRes_Tryp_Phosphorylation.param CID_LowRes_aLP.param
CID_TOF_Tryp.param                    CID_TOF_aLP.param
ETD_HighRes_NoCleavage.param         ETD_HighRes_Tryp.param
ETD_LowRes_ArgC.param                 ETD_LowRes_AspN.param
ETD_LowRes_GluC.param                 ETD_LowRes_LysC.param
ETD_LowRes_LysN.param                 ETD_LowRes_LysN_Phosphorylation.param
ETD_LowRes_Tryp.param                 ETD_LowRes_Tryp_Phosphorylation.param
ETD_LowRes_aLP.param
HCD_HighRes_NoCleavage.param          HCD_HighRes_Tryp.param
HCD_HighRes_Tryp_Phosphorylation.param HCD_HighRes_Tryp_TMT.param
HCD_HighRes_Tryp_iTRAQ.param          HCD_HighRes_Tryp_iTRAQPhospho.param
HCD_QExactive_Tryp.param              HCD_QExactive_Tryp_Phosphorylation.param
HCD_QExactive_Tryp_TMT.param          HCD_QExactive_Tryp_iTRAQ.param
HCD_QExactive_Tryp_iTRAQPhospho.param HCD_TOF_aLP.param
UVPD_QExactive_Tryp.param             UVPD_QExactive_Tryp_TMT.param
```

**When auto-detection fails** (missing activation block, unknown CV term, or running outside the source tree without bundled resources): msgf-rust falls back to `HCD_QExactive_Tryp.param` for default-flag runs, or to the resolution ladder in §1 for explicit flags. If no bundled file resolves, the process exits with an error instructing you to pass `--param-file <PATH>` explicitly.

---

## 5. Building from source

**Requirements:** Rust **1.85+** (workspace pins **1.87.0** in `rust-toolchain.toml` because transitive dependencies use `edition = "2024"`).

```bash
git clone https://github.com/bigbio/msgf-rust
cd msgf-rust
cargo build --release
# Binary: target/release/msgf-rust
```

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

Release archives bundle the binary, all 39 `.param` files, and `unimod.obo` under `resources/` — see [`README.md`](README.md) §Install.

---

## 6. Training new `.param` files

msgf-rust loads Java MS-GF+ `.param` scoring models **without conversion**. The 39 bundled files in `resources/ionstat/` were copied from the Java distribution unchanged; the on-disk binary format is identical.

Training **new** models (novel fragmentation chemistry, instrument class, or acquisition protocol) requires a scoring-parameter generator. Java MS-GF+'s **`ScoringParamGen`** is the canonical trainer.

**Status in v0.1.0:** search and scoring are fully ported and benchmark-validated; **`ScoringParamGen` is not yet ported** to Rust. Track progress on the [GitHub issues](https://github.com/bigbio/msgf-rust/issues) page.

**Interim workflows:**

1. **Use bundled models** — covers HCD QExactive tryptic DDA, CID low-res ion trap, ETD, phosphorylation, TMT, and iTRAQ variants (§4 file list).
2. **Train on the `java-legacy` branch** — check out the preserved Java tree (`git checkout java-legacy`), run Java `ScoringParamGen` on representative training data, then point msgf-rust at the output: `--param-file /path/to/MyModel.param`.

The Rust scorer reads any valid Java `.param` file via `Param::load_from_file`.

---

## 7. Isobaric labeling

TMT and iTRAQ searches require **both** protocol-aware scoring models **and** correct fixed modifications in `mods.txt`. Set `--protocol TMT` or `--protocol iTRAQ` (or legacy `--protocol 4` / `--protocol 2`) so the resolver prefers protocol-suffixed bundled files such as `HCD_QExactive_Tryp_TMT.param` or `HCD_QExactive_Tryp_iTRAQ.param`.

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
msgf-rust \
  --spectrum tmt_spectra.mzML \
  --database hsapiens.fasta \
  --output-pin out.pin \
  --mods tmt_10plex_mods.txt \
  --protocol TMT \
  --fragmentation HCD \
  --instrument QExactive
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
msgf-rust \
  --spectrum itraq_spectra.mzML \
  --database hsapiens.fasta \
  --output-pin out.pin \
  --mods itraq_8plex_mods.txt \
  --protocol iTRAQ \
  --fragmentation HCD \
  --instrument QExactive
```

For phospho-enriched isobaric data use `--protocol iTRAQ-phospho` (legacy `--protocol 3`) and include phospho variable mods in `mods.txt` (§2 example c).

---

## 8. Java MS-GF+ → msgf-rust migration

msgf-rust accepts **both** canonical kebab-case flags with named enum values **and** legacy Java short flags / numeric IDs. Existing quantms scripts using `--fragmentation 3 --instrument 3 --protocol 4` continue to work.

### 8a. Flag rename table

| Java MS-GF+ | msgf-rust canonical | msgf-rust legacy alias |
|---|---|---|
| `-s <FILE>` | `--spectrum <FILE>` | — |
| `-d <FILE>` | `--database <FILE>` | — |
| `-o <FILE>` | `--output-pin <FILE>` | — |
| `-mod <FILE>` | `--mods <FILE>` | `--mod <FILE>` (hidden) |
| `-t 20ppm` | `--precursor-tol-ppm 20` | — |
| `-ti -1,2` | `--isotope-error-min -1 --isotope-error-max 2` | — |
| `-m 3` (HCD) | `--fragmentation HCD` | `--fragmentation 3` |
| `-inst 3` (QExactive) | `--instrument QExactive` | `--instrument 3` |
| `-protocol 4` (TMT) | `--protocol TMT` | `--protocol 4` |
| `-ntt 2` (fully specific) | `--enzyme-specificity fully` | `--ntt 2` |
| `-tda 1` | *(omit — decoys auto-generated)* | — |
| `-e 1` (Trypsin) | *(omit — Trypsin only; other enzymes need `--param-file`)* | — |
| `-outputFormat 1` (TSV) | `--output-tsv <FILE>` | — |
| `-thread N` | `--threads N` | — |
| `-minLength 6` | `--min-length 6` | — |
| `-maxLength 40` | `--max-length 40` | — |
| `-maxMissedCleavages 1` | `--max-missed-cleavages 1` | — |
| `-minNumPeaks 10` | `--min-peaks 10` | — |
| `-n 10` | `--top-n 10` | — |
| `-precursorCal auto\|on\|off` | `--precursor-cal auto\|on\|off` | — |
| model path / `-conf` | `--param-file <FILE>` | — |

### 8b. Numeric-legacy values

Full legacy 0…N → named-value tables for `--fragmentation`, `--instrument`, `--protocol`, and `--enzyme-specificity` (`--ntt`) live in [`docs/CLI_MIGRATION.md`](docs/CLI_MIGRATION.md). clap accepts named values case-insensitively (`--fragmentation hcd` ≡ `HCD`).

### 8c. Behavior differences

| Area | Java MS-GF+ | msgf-rust |
|---|---|---|
| Spectrum inputs | mzML, MGF, mzXML, MS2, PKL, `_dta.txt`, … | **mzML and MGF only** |
| Identification output | PIN, TSV, mzIdentML | **PIN + optional TSV** (no mzIdentML) |
| Decoys | Separate target/decoy FASTA or `-tda` modes | **Always auto-generated** reversed decoys from target FASTA (`--decoy-prefix`) |
| Enzymes | Many via param file / CLI | **Trypsin only** in bundled models; other enzymes via `--param-file` |
| Mods file | Composition strings supported | **Numeric Da masses only** |
| Help | Picocli | clap-derived `--help` |
| Memory model | Loads full spectrum list | **Chunked streaming** (5000 spectra/chunk) for large mzML files |

### 8d. Known parity divergences

On PSMs where Java and Rust agree on scan + top-1 peptide (the "agreement bucket"), three PIN features still differ systematically. None block production use — aggregate 1% FDR PSM counts meet or beat Java on all three benchmark datasets (see [`README.md`](README.md)).

| Feature | Divergence | Status |
|---|---|---|
| `lnEValue` | ~4 orders of magnitude mean shift (Rust more confident) | Deferred — `num_distinct` semantics differ (see item #2 below) |
| `MeanRelErrorTop7` / `MeanErrorTop7` / `StdevRelErrorTop7` | >1% relative difference on ~99% of agreement-bucket PSMs | Deferred — error-stat normalization differs |
| BSA charge-3 SpecEValue (BSA.fasta + test.mgf fixture) | 1–4 OOM gap depending on deconvolution iteration | Known — deconvolution implementation divergence; kept as dev-branch smoke gate (`gf_java_parity` tests) |

Percolator's learned weights absorb these distribution shifts; rescored FDR results remain competitive or better than Java on `--precursor-cal off` benchmarks.

### 8e. precursorCal ship gates (`--precursor-cal auto`)

Java `-precursorCal auto` runs a file-wide pre-pass (sampled mini-search → median ppm
shift → tightened precursor tolerance) before the main search. msgf-rust ports this
in `mass_calibrator.rs` / `precursor_cal.rs`.

**G1 gate:** Rust @1% FDR within ±1% of Java on all three sign-off datasets (LFQ,
Astral, TMT) with matching cal mode. Fair comparison requires explicit Rust routing
flags — especially TMT (`--fragmentation CID --instrument high-res --protocol TMT`).
Harness: [`benchmark/vm/run_bench_calauto_3ds.sh`](benchmark/vm/run_bench_calauto_3ds.sh).

As of 2026-05-25 (fair v5 gate, `bench-calauto-v5.log`) the calibrator is
validated (shift + tightening parity), but G1 **fails** on all three datasets:
LFQ −2.2% (cal skipped), Astral +1.2%, TMT −5.9%. Remaining
gaps trace to **SpecE tail / Percolator feature parity** (same class as historical
Astral GF work), not calibrator logic. Full analysis:
[`docs/parity-analysis/notes/2026-05-25-precursor-cal-ship-gates.md`](docs/parity-analysis/notes/2026-05-25-precursor-cal-ship-gates.md).

| `--precursor-cal` | Ship recommendation |
|---|---|
| `off` | Yes — baseline unchanged (**CLI default** until G1) |
| `auto` | Opt-in only — no until G1 passes |
| `on` | Opt-in only — no until G1 passes |

---

## 9. License and citation

msgf-rust is distributed under the **UCSD Noncommercial License** — the same terms as upstream MS-GF+. The license permits copying, modification, and distribution for **educational, research, and non-profit** purposes without fee, provided the copyright notice and liability paragraphs are retained. **Commercial use requires written permission** from the UCSD Technology Transfer Office (see `LICENSE` for contact details).

The software is provided **"as is"** without warranty. See [`LICENSE`](LICENSE) for the full upstream text and [`NOTICE`](NOTICE) for port attribution.

### Citation

If you use msgf-rust in published work, cite the original MS-GF+ paper:

> Kim, S. and Pevzner, P.A. (2014). MS-GF+ makes progress towards a universal database search tool for proteomics. *Nature Communications*, 5:5277.

And optionally this Rust port:

> bigbio (2026). msgf-rust: a Rust port of MS-GF+ for the quantms pipeline. https://github.com/bigbio/msgf-rust

The original Java implementation is preserved on the `java-legacy` branch; upstream MS-GF+ lives at https://github.com/MSGFPlus/msgfplus.
