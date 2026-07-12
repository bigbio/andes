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
| `.raw` | Thermo RawFileReader | `--features thermo` (release archives ship it) | .NET 8 runtime — **bundled in the release archives** (nothing to install); from source, install .NET 8 | Native Thermo; results are identical to searching the equivalent mzML. Supports `--chimeric`. Activation/instrument read from vendor metadata (§4). |
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
| `--fragmentation` | enum | `auto` | Fragmentation method for bundled model resolution. Named: `auto`, `CID`, `ETD`, `HCD`, `UVPD`. `auto` on mzML triggers activation detection (§4); on MGF falls back to bundled defaults. | Java `-m`; numeric `0`=auto, `1`=CID, `2`=ETD, `3`=HCD, `4`=UVPD |
| `--instrument` | enum | `low-res` | Instrument class for bundled model resolution. Named: `low-res`, `high-res`, `TOF`, `QExactive`. | Java `-inst`; numeric `0`=low-res, `1`=high-res, `2`=TOF, `3`=QExactive |
| `--protocol` | enum | `auto` | Search protocol suffix for bundled model resolution. Named: `auto`, `phospho`, `iTRAQ`, `iTRAQ-phospho`, `TMT`, `standard`. | Java `-protocol`; numeric `0`=auto, `1`=phospho, `2`=iTRAQ, `3`=iTRAQ-phospho, `4`=TMT, `5`=standard |
| `--param-file` | path | *(auto)* | Explicit path to a `.param` scoring model file. When set, overrides all auto-detection and bundled resolution. Required when running a release binary outside the source tree if bundled resources are not present. | Java `-conf` / model path |
| `--model-store` | path | *(bundled)* | Path to a Parquet model store to use instead of the bundled `resources/models.parquet`. Model selection reads from this store when set. | *(no Java equivalent)* |
| `--model` | string | *(auto-select)* | Exact model ID to load from the model store, skipping automatic selection by `(--fragmentation, --instrument, --protocol)`. Useful for searching with a freshly-trained model (see `andes train`). | *(no Java equivalent)* |

**Bundled default when all scoring flags are at their defaults** (`--fragmentation auto --instrument low-res --protocol auto`): `hcd_qexactive_tryp` (from the parquet model store). This preserves pre-auto-detect behaviour for MGF inputs and mzML files without activation metadata.

**Model selection** (when `--param-file` is not set, resolved from `resources/models.parquet`):

1. Build a selection key: `{Frag}_{Inst}_Trypsin` with optional protocol experiment class (e.g. `tmt`).
2. Exact match on the key → use that model.
3. If protocol-specific model absent, retry without the protocol class.
4. Final fallback: `cid_tof_tryp` (HCD + TOF/HighRes), `etd_lowres_tryp` (ETD), or `cid_lowres_tryp` (everything else).

**Normalisation rules:**

- `auto` fragmentation → treated as `CID` for model selection (except mzML auto-detect path, §4).
- HCD + `low-res` instrument → upgraded to `QExactive`.

Only tryptic enzyme models are in the store; other enzymes require `--param-file` with a binary `.param` file.

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

**Environment variable:** set `ANDES_RSS_PROBE=1` on Linux to print `VmRSS` checkpoints to stderr during long runs (debugging memory use). See §9 for the full list of internal environment variables.

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

**Comments:** lines whose first non-whitespace character is `#` are ignored. Inline `# ...` comments are stripped from the end of a line. Blank lines are ignored.

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

Pair with `--protocol TMT --fragmentation HCD --instrument QExactive` to select the `hcd_qexactive_tryp_tmt` model from the store (§4, §7).

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

Tab-separated, one header row, one row per PSM. Rows are sorted best-first within each spectrum by `RankScore` (the GF-free rank-LLR score) — the generating function and all of its derived score columns have been removed. The `chargeN` one-hots track the `--charge-min`…`--charge-max` range: one column per charge state, so narrowing/widening the range removes/adds one `chargeN` column each (e.g. a 2–3 range yields just `charge2 charge3`). With the default 2–5 range the full column set is the 66 columns listed below in emission order.

There are **two score columns**, easy to confuse:

* **`RankScore`** (col 7) — the rank-LLR score; the **ranking** signal that orders candidates within a spectrum (this was historically called `RawScore`).
* **`RawScore`** (col 62) — the fused strong-score `signal − null`; the **headline discriminative** feature Percolator weights most (historically `StrongScore`). With `--score strong` it also becomes the ranking signal.

Most of the columns after `matchedIonRatio` are **additive** features: extra evidence Percolator can learn weights for without perturbing the core score distribution. Several are **0.0 unless a flag/model is active** — see the note after the table.

`--chimeric` does **not** change the column set. It populates `PrecursorIsotopeKL` / `PrecursorSNR` (`0.0` otherwise) from a linked MS1, and — because a scan can then emit several rows — appends a per-row index to multi-row `SpecId`s (see below). Multi-row scans also occur without `--chimeric` whenever rank-1 candidates tie.

| # | Column | Type | Range | Description |
|---|---|---|---|---|
| 1 | `SpecId` | string | — | `{specID}_{scan}_{rank}` PSM id; multi-row scans get a `_{rowIdx}` suffix to stay unique. |
| 2 | `Label` | int | {−1, +1} | `+1` target, `−1` decoy (by **source protein**, TDC convention). |
| 3 | `ScanNr` | int | ≥0 | MS2 scan number. |
| 4 | `ExpMass` | float | >0 | Experimental neutral precursor mass (Da): `mz×z − z×proton`. |
| 5 | `CalcMass` | float | >0 | Theoretical neutral peptide mass (Da, incl. H₂O). |
| 6 | `mass` | float | >0 | Duplicate of `ExpMass` (PercolatorAdapter convention). |
| 7 | `RankScore` | int | unbounded | **Rank-LLR ranking score** (orders candidates within a spectrum). |
| 8 | `isotope_error` | int | [−1, 2] | Winning ¹³C isotope offset. |
| 9 | `peplen` | int | ≥6 | Residue count **+ 2** (includes flanking pre/post). |
| 10 | `dm` | float | signed | Precursor mass error (Da) after isotope correction. |
| 11 | `absdm` | float | ≥0 | `\|dm\|`. |
| 12–15 | `charge2`…`charge5` | 0/1 | one-hot | One-hot precursor charge; one column per state in `--charge-min`…`--charge-max`. |
| 16 | `enzN` | 0/1 | one-hot | N-terminal boundary consistent with the enzyme rule. |
| 17 | `enzC` | 0/1 | one-hot | C-terminal boundary consistent with the enzyme rule. |
| 18 | `enzInt` | int | ≥0 | Count of internal positions matching the enzyme rule. |
| 19 | `NumMatchedMainIons` | int | [0, peplen−1] | Matched charge-1 b/y fragment positions. |
| 20 | `longest_b` | int | [0, peplen−1] | Longest contiguous matched b-ion run. |
| 21 | `longest_y` | int | [0, peplen−1] | Longest contiguous matched y-ion run. |
| 22 | `longest_y_pct` | float | [0, 1] | `longest_y / peplen`. |
| 23 | `ExplainedIonCurrentRatio` | float | [0, 1] | Matched b+y intensity / total MS2 ion current. |
| 24 | `NTermIonCurrentRatio` | float | [0, 1] | Matched b-ion intensity / total MS2 ion current. |
| 25 | `CTermIonCurrentRatio` | float | [0, 1] | Matched y-ion intensity / total MS2 ion current. |
| 26 | `MS2IonCurrent` | float | ≥0 | Sum of all MS2 peak intensities (not log-scaled). |
| 27 | `IsolationWindowEfficiency` | float | 0.0 | Always `0.0` (not available from parsed spectra). |
| 28 | `MeanErrorTop7` | float | ≥0 | Mean absolute ppm error of the top-7 most-intense matched ions. |
| 29 | `StdevErrorTop7` | float | ≥0 | Population stdev of absolute ppm errors (top-7). |
| 30 | `MeanRelErrorTop7` | float | signed | Mean signed ppm error (top-7). |
| 31 | `StdevRelErrorTop7` | float | ≥0 | Population stdev of signed ppm errors (top-7). |
| 32 | `matchedIonRatio` | float | [0, 1] | `NumMatchedMainIons / peplen`. |
| 33 | `EdgeScore` | int | unbounded | Per-bond edge-score sum (ion-existence + error); additive (Kim et al. 2014). |
| 34 | `PrecursorIsotopeKL` | float | ≥0 | KL divergence of precursor isotope envelope vs averagine. **0.0 unless `--chimeric`.** |
| 35 | `PrecursorSNR` | float | ≥0 | Precursor SNR from the MS1 envelope. **0.0 unless `--chimeric`.** |
| 36 | `DeltaRankScore` | float | ≥0 | `RankScore(best) − RankScore(2nd-best distinct peptide)`; rank-1 row only, else 0.0. |
| 37 | `TailorScore` | float | ≥0 | `RankScore ÷` spectrum's top-1% quantile (Yang et al. 2020); cross-spectrum comparability. |
| 38 | `PpmGaussianScore` | float | ≥0 | `Σ exp(−½(ppm/7)²)` over matched ions — mass-accuracy evidence the rank score discards. |
| 39 | `NeutralLossIonCount` | int | ≥0 | Matched b/y ions with −H₂O/−NH₃ partner peaks. |
| 40 | `LongestComplementaryLadder` | int | [0, peplen−1] | Longest run of bonds where both bᵢ and y₍ₙ₋ᵢ₎ matched. |
| 41 | `ComplementaryIonBalance` | float | ≥0 | `Σ 1/(1+\|rankᵦ−rankᵧ\|)` over complementary bonds. |
| 42 | `MeanMatchedIntensityRank` | float | ≥1 | Mean intensity-rank of matched ions (1 = most intense; lower is better). |
| 43 | `DoublyChargedMatchedIonCount` | int | ≥0 | Matched charge-2 b/y ions. |
| 44 | `UniqueMatchFraction` | float | [0, 1] | Within-peptide peak-explanation uniqueness. |
| 45 | `ChanceMatchSurprise` | float | ≥0 | `Σ max(0, −ln(ρ·Δ))` — how improbable the matches are by chance (null moat). |
| 46 | `IntensitySignal` | float | [0, 1] | Cosine sim. of predicted vs observed intensities. **0.0 without an intensity model.** |
| 47 | `FragPredExplained` | float | [0, 1] | `Σ(matched·pred)/Σpred`. **0.0 without a frag-intensity model.** |
| 48 | `FragPredChanceLLR` | float | ≥0 | `Σ matched·pred·max(0,−ln p_chance)`. **0.0 without a frag-intensity model.** |
| 49 | `FragTopKObserved` | float | [0, 1] | Top-K predicted-most-intense ions observed. **0.0 without a frag-intensity model.** |
| 50 | `RichIonLLR` | float | unbounded | Decoy-aware per-annotated-ion LLR sum. **0.0 without a rich-ion model.** |
| 51 | `IsRefinement` | 0/1 | one-hot | 1 if the PSM came from the Pass-2 refinement search. **0 without `--refine`.** |
| 52 | `NumMods` | int | ≥0 | Variable-modification count on the matched peptide. |
| 53 | `RefinementModClass` | int | [0, 99] | Mod-class id for subgroup-FDR grouping. **0 without `--refine`.** |
| 54 | `ModSiteShiftedMatched` | int | ≥0 | Matched mod-bearing (mass-shifted) b/y ions. **0 for unmodified peptides.** |
| 55 | `ModSiteShiftedFrac` | float | [0, 1] | Matched shifted ÷ total shifted ions. |
| 56 | `ModSiteIntensFrac` | float | [0, 1] | Shifted-ion intensity ÷ all matched-ion intensity. |
| 57 | `ModSiteLocalized` | 0/1 | one-hot | 1 if a bracketing ion pair localizes the mod. |
| 58 | `ModSiteDetCount` | int | ≥0 | Count of site-determining (bracketing) ions over all mod sites. |
| 59 | `MassCompetitionEvidence` | float | ≥0 | `Σ 1/(1+ambiguity+ρ)` — alternative-mass competition null term. |
| 60 | `CandidateRankEntropy` | float | ≥0 | Softmax entropy over the retained top-K candidate scores (spectrum-level). |
| 61 | `ListwiseScoreGap` | float | signed | Top-1 − top-2 `RankScore` in the retained queue. |
| 62 | `RawScore` | float | unbounded | **Headline fused strong-score** `signal − null` — the primary discriminative feature. |
| 63 | `RawScoreCal` | float | signed | Per-spectrum z-scored `RawScore` (significance calibration). |
| 64 | `RankScoreFloat` | float | unbounded | Unrounded `RankScore` (continuous split-sum) — finer-grained ranking feature than the integer `RankScore`. |
| 65 | `Peptide` | string | — | `pre.SEQUENCE.post` with `+mass` mod annotations. |
| 66 | `Proteins` | string | — | Protein accession(s), tab-separated for shared peptides; decoys carry `--decoy-prefix`. |

**Conditional columns** (always present in the header, but `0.0`/`0` unless their condition holds):

* `PrecursorIsotopeKL`, `PrecursorSNR` — need `--chimeric` + a linked MS1.
* `IntensitySignal`, `FragPredExplained`, `FragPredChanceLLR`, `FragTopKObserved` — need a trained intensity / frag-intensity model.
* `RichIonLLR` — needs a trained rich-ion model.
* `IsRefinement`, `NumMods`, `RefinementModClass`, `ModSite*` — populated by `--refine` (and the `ModSite*` block only on modified peptides).
* `DeltaRankScore` — emitted on the rank-1 row only.

### 3b. TSV columns

Tab-separated human-readable report. The `Title` column appears **only for MGF** inputs.

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
| `RawScore` | int | Rounded raw score — the sole score column (the generating function and its derived score columns have been removed). |

**mzML header** — same as above **without** the `Title` column (11 columns total).

Decoy PSMs are included in TSV output; downstream tools label them via Percolator or manual filtering.

### 3c. PIN vs TSV — which to use

Use **PIN** when the goal is FDR calibration or rescoring: Percolator, MS²Rescore, Mokapot, and quantms-style pipelines consume `.pin` directly and learn feature weights from the full Percolator feature set (including `EdgeScore`). Use **TSV** for spreadsheet inspection, custom reporting, or tools that expect a flat PSM table. You can emit both in one run with `--output-pin` and `--output-tsv`. For production quantms workflows, PIN is the standard path; TSV is optional diagnostics.

### 3d. Run summary (`statistics.log`)

andes auto-resolves the scoring model and the precursor/fragment tolerances from the input metadata, so the parameters a search **ends** with are not necessarily the CLI inputs: precursor calibration tightens the window, and a high-res model carries (e.g.) a 20 ppm fragment tolerance even when the input named none. To make a run's true parameters recoverable, andes prints a summary to stderr at the end of every search **and** writes a `statistics.log` next to the PIN (in the PIN's parent directory). Implementation: `crates/output/src/stats.rs`.

The summary records the **final** precursor tolerance (+ calibration mode), the **final** fragment tolerance (the resolved model's `mme`), the number of spectra with a match, the pre-FDR rank-1 target/decoy PSM split, and a **per-modification PSM tally** — for each modification (fixed like Carbamidomethyl and variable like Oxidation/Acetyl), how many rank-1 target PSMs carry it, plus an `(unmodified)` count.

```text
──────── andes run summary ────────
  Final precursor tolerance : Symmetric(10.0 ppm) (calibration: Auto)
  Final fragment tolerance  : 0.5 Da
  Spectra with a match      : 48210
  Rank-1 PSMs (pre-FDR)     : 31204 target, 17006 decoy
  PTM report (rank-1 target PSMs carrying each modification):
    Carbamidomethyl : 28933
    Oxidation       :  6120
    Acetyl          :   341
    (unmodified)    :  2150
  ───────────────────────────────────
```

Counts are **pre-FDR**, taken over each spectrum's best (rank-1) candidate; final FDR control happens downstream in Percolator. The tally is most useful with `--refine`, where it shows exactly which discovered PTMs were identified and at what volume. (`statistics.log` matches the gitignore `*.log*` pattern — it is a per-run output artifact, not a tracked file.)

### 3e. QPX `.idparquet` bundle (`--output-parquet`)

`--output-parquet <DIR>` writes an **OpenMS-compatible QPX 1.0** Parquet bundle — a directory (conventionally ending in `.idparquet`) containing `psms.parquet`, `proteins.parquet`, and `search_params.parquet`. The schema (column names, Arrow types, nested `list<element: …>` structures, and the per-file metadata keys `qpx_version`/`file_type`/`uuid`/`creation_date`/`software_provider`/`creator`) matches what OpenMS's `QPXFile` writer emits byte-for-byte, so the files are interchangeable with OpenMS / [quantms](https://github.com/bigbio/quantms) tooling. Implementation: `crates/output/src/qpx.rs`. Reuses the workspace's existing `arrow`/`parquet` stack — no new heavy dependency.

`psms.parquet` carries one row per PSM with `sequence`, `peptidoform`, `modifications` (name + Unimod accession + positions), `precursor_charge`, `calculated_mz`/`observed_mz`, `is_decoy`, `scan`/`rt`, `protein_accessions` (with flanks + offsets), the spectrum `mz_array`/`intensity_array`, the headline `score` (`andes:RawScore`), and an `additional_scores` list carrying the other andes features (`RankScore`, `TailorScore`, `DeltaRankScore`, `EdgeScore`, `RichIonLLR`, …). `search_params.parquet` records the resolved engine/tolerances/enzyme/modifications.

Fields andes does **not** compute pre-rescoring are written null: `posterior_error_probability` and the q-value are Percolator's job (downstream), and `predicted_rt`/`ion_mobility`/per-peak `charge_array`/`ion_type_array` are not produced. `proteins.parquet` lists the distinct accessions seen in PSMs (andes does no protein inference). Emit it alongside `--output-pin`/`--output-tsv`:

```bash
andes --spectrum spectra.mzML --database db.fasta \
  --output-pin out.pin --output-parquet out.idparquet
```

---

## 4. Auto-detection

For **mzML** inputs when `--fragmentation auto` (the default), andes peeks the input file before loading the full dataset:

1. **Activation method** — histogram of `<activation>` cvParams across the first 64 MS2 spectra; dominant method wins. Mixed methods trigger an stderr warning but the dominant method is still used file-wide.
2. **Instrument class** — scans `<instrumentConfiguration>` / analyzer cvParams via `input::detect_instrument_type`; dominant analyzer among MS2 spectra wins. `None` → `low-res` (the low-resolution ion-trap default).

Precedence: whether auto-detection runs is gated **only** by `--fragmentation auto` (the default) on an mzML/`.raw`/`.d` input — *not* by `--instrument`. When it runs and the peek succeeds, the **detected** instrument is used and any `--instrument` value on the command line is **ignored** for model selection; to force an instrument, set an explicit `--fragmentation` (e.g. `HCD`) so the auto path is disabled and the flags drive resolution (§1). `--protocol` from the CLI is always applied to pick protocol-specific models from the parquet store (e.g. the `tmt` experiment-class entry).

MGF files carry no activation or instrument metadata → auto-detect returns `None` → bundled default `hcd_qexactive_tryp` model (from the parquet store) unless explicit `--fragmentation` / `--instrument` flags override the store selection key.

Non-auto `--fragmentation` (e.g. `HCD`, `3`) disables the activation peek and uses flag-based resolution directly (§1), including `--instrument` and `--protocol` from the CLI.

### Native Thermo `.raw`

A `.raw` file carries the activation method and analyzer in vendor metadata, so andes reads them directly (no mzML peek) and routes through the same parquet-store selection as mzML — e.g. beam-type CID (HCD) on an Orbitrap → `hcd_qexactive_tryp`. `--protocol` from the CLI still selects protocol-specific models (`tmt`, `itraq`); explicit `--fragmentation`/`--instrument` are not required.

### Native Bruker timsTOF `.d`

timsTOF DDA-PASEF is beam-type CID on a TOF analyzer, so `.d` input auto-routes to the **`cid_tof_tryp`** model in the parquet store. `--protocol` still applies. Searched **MS2 only**; the ion-mobility dimension is carried as metadata but not used by scoring.

### Activation CV mapping (mzML `<activation>` cvParam accession → method)

| CV accession | Name (PSI-MS) | andes method | Notes |
|---|---|---|---|
| `MS:1000133` | collision-induced dissociation | CID | |
| `MS:1000422` | beam-type collision-induced dissociation (HCD) | HCD | |
| `MS:1000598` | electron transfer dissociation | ETD | |
| `MS:1000599` | pulsed Q dissociation | CID | PQD is scored as CID |
| `MS:1000435` | photodissociation | UVPD | |
| `MS:1000250` | electron capture dissociation | ETD | Mapped to ETD (no dedicated ECD variant) |

### Instrument detection (analyzer cvParam → class)

| Analyzer family | Examples | Instrument class |
|---|---|---|
| Ion trap / linear ion trap | `MS:1000264`, Velos, LTQ | `low-res` |
| Orbitrap / Fusion | `MS:1000480`, Fusion Lumos | `QExactive` |
| FT-ICR | `MS:1000480` (FT) | `high-res` |
| TOF | `MS:1000128` | `TOF` |

### Bundled model store (`resources/models.parquet`)

All 39 scoring models ship with the binary as a single Parquet model store
(`resources/models.parquet`). The store covers the full
fragmentation × instrument × protocol matrix (CID/ETD/HCD/UVPD ×
LowRes/HighRes/TOF/QExactive × Trypsin, with protocol variants for Phospho, TMT,
iTRAQ, iTRAQPhospho).

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

andes includes a native Rust training engine — **`andes train`** — that generates scoring models from your own data and writes them into the same Parquet model store the bundled models live in.

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

andes ships its own model store at `resources/models.parquet`, containing all 39 bundled scoring models. The `--param-file` flag can additionally load an external binary model file directly for custom or externally supplied models.

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
andes \
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

## 8. Legacy numeric values & behavior notes

For backward compatibility, the routing flags accept legacy 0…N numeric values in
addition to their canonical named values; clap parses named values
case-insensitively (`--fragmentation hcd` ≡ `HCD`).

| Flag | Numeric | Named |
|---|---|---|
| `--fragmentation` | `0` | `auto` |
| `--fragmentation` | `1` | `CID` |
| `--fragmentation` | `2` | `ETD` |
| `--fragmentation` | `3` | `HCD` |
| `--fragmentation` | `4` | `UVPD` |
| `--instrument` | `0` | `low-res` |
| `--instrument` | `1` | `high-res` |
| `--instrument` | `2` | `TOF` |
| `--instrument` | `3` | `QExactive` |
| `--protocol` | `0` | `auto` |
| `--protocol` | `1` | `phospho` |
| `--protocol` | `2` | `iTRAQ` |
| `--protocol` | `3` | `iTRAQ-phospho` |
| `--protocol` | `4` | `TMT` |
| `--protocol` | `5` | `standard` |
| `--enzyme-specificity` (alias `--ntt`) | `0` | `non-specific` |
| `--enzyme-specificity` (alias `--ntt`) | `1` | `semi` |
| `--enzyme-specificity` (alias `--ntt`) | `2` | `fully` |

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

## 9. Glycopeptide search (experimental) & advanced knobs

Enable N-glycopeptide search with `--glyco` (requires a `.glyco.pin` output; the
backbone model is the N-X-S/T sequon). Cross-spectrum backbone transfer is an
opt-in second pass via `--glyco-transfer`. All glyco tuning is exposed as **hidden
CLI flags** (advanced; the shipped defaults are validated and rarely need changing):

| Flag | Default | Purpose |
|---|---|---|
| `--glyco-backbone-top-k` | 50 | Max backbone candidates per spectrum (set large to approximate an exhaustive ceiling). |
| `--glyco-gp-k` / `-gp-j` / `-gp-h` | 50 / 5 / 1 | Weights of the `gp` fused per-scan selector `rank + K·ladder + J·core_y + H·hyper` (the shipped, default selector). |
| `--glyco-pf-charge` | 2 | Charge states the peptide-first fragment index covers (b/y at 1..=N, clamped 1..=3); targets high-charge glycopeptides. |
| `--glyco-max-pf` | 1024 | Max peptide-first candidates per spectrum. |
| `--glyco-decoy` | off | Emit paired glycan-axis decoy rows for experimental 2D (peptide × glycan) FDR. |
| `--glyco-transfer` | off | Enable cross-spectrum backbone transfer (two-pass). |
| `--glyco-transfer-seed-fdr` | 0.05 | q-value threshold for confident donor seeds. |
| `--glyco-rt-window` | 1800 | RT co-elution window (seconds) for transfer. |
| `--glyco-transfer-ungated` | off | Skip the RT co-elution gate (unsafe research opt-in). |
| `--glyco-transfer-min-support` | 1 | Minimum independent-donor graph support to inject a transfer. |
| `--glyco-transfer-core-y` | 3 | Acceptor-side core-Y quorum (incl. mandatory Y1) to accept a transfer. |
| `--debug-glyco` | off | **Diagnostic only:** emit all candidate rows per scan (incl. de-novo) + transfer diagnostics. The resulting PIN must NEVER be fed to an FDR tool. |

FDR is always computed externally (Percolator) — andes never computes FDR itself.

### Internal environment variables

A few remaining `ANDES_*` environment variables are **internal / advanced** and are
NOT part of the supported search interface. They do not change default search output.

- **Model training** (only read by the hidden `train*` subcommands): `ANDES_GEO_SEGMENTS`, `ANDES_GEO_MAX_RANK`, `ANDES_GEO_OCCUPANCY`, `ANDES_GEO_MAX_TIERS`, `ANDES_GEO_MAX_FRAG_CHARGE` (partition-geometry derivation), `ANDES_SEED_GEOMETRY` (reuse seed geometry), `ANDES_DENSE_NOISE` (noise sampler), `ANDES_V1_STORE` / `ANDES_V1_OUT`, `ANDES_TRAIN_BENCH`.
- **Advanced scoring**: `ANDES_PEAK_WINDOW` / `ANDES_PEAK_PER_WINDOW` (windowed peak filtering; unset = model-tolerance default).
- **Read-only developer instrumentation** (no effect on output): `ANDES_RSS_PROBE` (memory logging), `ANDES_CHIMERIC_OVERLAP` (fragment-overlap diagnostic), `Andes_TRACE_IONS` / `Andes_TRACE_PEP` (ion/peptide trace logging).
- **Test harness only**: `ANDES_TEST_D`, `ANDES_TEST_RAW`, `ANDES_TEST_PERCOLATOR_BIN`.

---

## 10. License and citation

andes is licensed under the **Apache License 2.0**. See [`LICENSE`](LICENSE) for the full text and [`NOTICE`](NOTICE) for attribution and the project's origin in MS-GF+.

The software is provided **"as is"** without warranty.

### Citation

If you use andes in published work, please cite both andes and the foundational MS-GF+ paper:

> bigbio (2026). andes: a data-driven peptide search engine for the quantms ecosystem. https://github.com/bigbio/andes

> Kim, S. and Pevzner, P.A. (2014). MS-GF+ makes progress towards a universal database search tool for proteomics. *Nature Communications*, 5:5277.

andes originated from MS-GF+ (https://github.com/MSGFPlus/msgfplus); see [`NOTICE`](NOTICE).
