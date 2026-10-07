# andes documentation

The full reference for the `andes` binary and its outputs; install and quick start are in
[`README.md`](README.md). `andes --help` shows the flags a normal run needs. The rest are
**advanced**: hidden from `--help` because the engine derives them and reports what it chose.
Use them to reproduce a measurement or override a derivation.

Advanced (hidden) flags include: `--candidate-index`, `--fragment-index`, `--gbdt-max-trees`,
`--peak-filter`, `--ethcd-activation`, `--isotope-error`, and the glyco tuning group
(`--glyco-tol-ppm`, `--glyco-retrieval-tol-ppm`, `--glyco-max-peaks`, `--glyco-y-max-charge`,
`--glyco-cz-max-charge`, `--glyco-scans`).

---

## 1. CLI reference

Flags are kebab-case long options, defined in `crates/andes/src/bin/andes/cli.rs`.

### Input formats

`--spectrum` picks the reader from the file extension; there is no format flag.

| Extension | Reader | Build requirement | Runtime requirement | Notes |
|---|---|---|---|---|
| `.mzML` / `.mzml` | mzML (streaming) | always built | none | Full activation + instrument auto-detection (§4). |
| `.raw` | Thermo RawFileReader | `--features thermo` | .NET 8 (bundled in release archives) | Same results as the equivalent mzML. Supports `--chimeric`. |
| `.d` | Bruker timsTOF (`timsrust`) | `--features timstof` | none | DDA-PASEF, MS2 only; routed to `cid_tof_tryp`. `--chimeric` / `--precursor-cal` degrade to a normal search. |
| any other (e.g. `.mgf`) | MGF | always built | none | No metadata; MS2, model from flags. |

Native `.raw`/`.d` search MS2 only (MS3 reporter scans are filtered at load).

### Required

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--spectrum` | path | *(required)* | Input spectrum file; reader chosen by extension (see above). | Java `-s <FILE>` |
| `--database` | path | *(required)* | Target FASTA. Decoys are generated (see `--decoy-strategy`, `--decoy-prefix`). | Java `-d <FILE>` |
| `--output-pin` | path | *(required)* | Output Percolator `.pin`. | Java `-o <FILE>` (when `-outputFormat pin`) |

### Search parameters

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--precursor-tol` | string | `20ppm` | Symmetric precursor tolerance, e.g. `20ppm` or `0.02da`. | Java `-t 20ppm` |
| `--enzyme` | enum | `trypsin` | `trypsin`, `chymotrypsin`, `lysc`, `aspn`, `gluc`, `lysn`, `argc`, `alphalp`, `nocleavage`, `nonspecific` (alias `elastase`). A comma list (`gluc,trypsin`) uses every enzyme listed. | Java `-e` |
| `--charge` | `MIN..MAX` | `2..5` | Charges tried when the spectrum has none. | *(no direct Java flag; set via param file in Java)* |
| `--enzyme-specificity` | enum | `fully` | Tolerable termini: `fully` (Java `-ntt 2`), `semi` (`-ntt 1`), `non-specific` (`-ntt 0`). | Java `-ntt` |
| `--max-missed-cleavages` | u32 | `1` | Missed cleavages per peptide. | Java `-maxMissedCleavages 1` |
| `--min-length` | u32 | `6` | Minimum peptide length. | Java `-minLength 6` |
| `--max-length` | u32 | `50` | Maximum peptide length. | Java `-maxLength 40` |
| `--top-n` | u32 | `10` | PSMs kept per spectrum. | Java `-n 10` |
| `--isotope-error` | `MIN..MAX` | `-1..2` | Isotope-error offsets tried. | Java `-ti -1,2` |
| `--min-peaks` | u32 | `10` | Spectra with fewer MS2 peaks are skipped. | Java `-minNumPeaks 10` |

### Modifications

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--mods` | path | *(off)* | Java-format `mods.txt` (§2). Without it: Carbamidomethyl C (fixed), Oxidation M (variable, max 3 per peptide). Numeric Da masses only. | Java `-mod <FILE>` |

### Scoring

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--fragmentation` | enum | `auto` | `auto`, `CID`, `ETD`, `HCD`, `UVPD`. `auto` reads the file (§4); on MGF it falls back to CID and warns. | Java `-m` |
| `--protocol` | enum | `auto` | `auto`, `phospho`, `iTRAQ`, `iTRAQ-phospho`, `TMT`, `standard`. An explicit value selects the protocol model. `auto` keeps the model; if it finds TMT/iTRAQ reporter ions it turns on the isobaric peak filter and (with no `--mods`) adds the tag as a fixed mod. | Java `-protocol` |
| `--score` | enum | `auto` | What ranks candidates and fills `RawScore`: `rank` (low-res), `strong` (high-res), or `auto` (by the model's instrument class). | — |
| `--gbdt-max-trees` | u32 | `100` | Trees per GBDT ensemble, `0` = all. 100 is 33–41% faster and identification-neutral (2026-09). `--glyco` uses all trees unless set. | — |
| `--peak-filter` | `WINDOW_DA:PEAKS` | protocol default | Keep the `PEAKS` most intense peaks per `WINDOW_DA`. Unset = `100:20` for isobaric data, else off; window `0` forces off. | — |
| `--ethcd-activation` | enum | `hcd` | EThcD/ETciD routing: `hcd` (no EThcD model exists) or `etd` (c/z path). | — |
| `--model-store` | path | *(bundled)* | Model store instead of `resources/models/` (directory or one `models.parquet`). | — |
| `--model` | string | *(auto-select)* | Load this model id, skipping selection. | — |

**Model selection** tries an exact (activation, instrument class, enzyme, protocol) match, then
drops the protocol, then takes the closest instrument class. The instrument class comes from
the file; use `--model-store` plus `--model` for an unbundled model.

### Calibration

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--precursor-cal` | enum | `off` | `off`, `auto`, `on`. A pre-pass learns a ppm shift from confident PSMs and tightens the precursor tolerance; `auto` skips it on small samples. Off by default: on the benchmark sets the pre-pass took 9–17% of the run and never fired, and where it fired (UPS1, larger sample) the tightened window lost 7.6% of PSMs. Skipped (with a warning) on `.raw` and `.d`. | Java `-precursorCal auto\|on\|off` |

**Candidate retrieval is automatic.** On high-res fragment matching, andes builds the candidate
index out-of-core (cached as `andes-candidx-<hash>.bin` in the system temp directory, ~150–250 MB
on the benchmark databases, reused by later runs on the same database and settings) and retrieves
the 100 candidates whose singly-charged b/y ions best match each spectrum (at least 3 matched).
On Astral this gave +22.2% PSMs over full-window enumeration at an equal true FDP of 1.00%; on
phospho, 164 s against 5,954 s. Low-res data stays on in-RAM enumeration (forcing the index took
TMT from 12,281 to 3,613 PSMs at 1% and UPS1 from 15,838 to 10,312), as do `--refine` and
`--glyco`; `--chimeric` runs on either path. The engine prints its choice; `--candidate-index ram` or
`--fragment-index off` restore enumeration
([measurements](docs/benchmarks/README.md#choosing-the-retrieval-strategy)).

### Runtime

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--threads` | usize | logical CPU count | Worker threads. | Java `-thread N` |
| `--ms-level` | u8 | `2` | MS level to search in mzML. `.raw`/`.d`, MGF and `--chimeric` always search MS2. | — |
| `--max-spectra` | usize | `0` | Bench mode: first N MS2 spectra only (`0` = all); skips TSV output. | — |

### Output

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--output-tsv` | path | *(off)* | Tab-separated PSM report (§3b). | Java `-outputFormat 1` with output path |
| `--output-parquet` | dir | *(off)* | QPX `.idparquet/` bundle (§3e). | — |

No environment variable changes a search result; the binary reads none. The test-harness
variables are listed in [`docs/ENV_VARS.md`](docs/ENV_VARS.md).

---

## 1a. Workflow parameters (grouped by experimental design)

Opt-in modes, each switched on by one parent flag.

### Decoys & FDR strategy

For a pre-built target+decoy FASTA, use `--decoy-strategy none` with `--decoy-prefix` or
`--decoy-suffix`.

| Flag | Type | Default | Description |
|---|---|---|---|
| `--decoy-strategy` | enum | `reverse` | `reverse`, `shuffle`, `sequon-reverse`, or `none`. **Use `sequon-reverse` with `--glyco`**: plain reversal turns N-X-S/T into S/T-X-N and makes glyco q-values anti-conservative. |
| `--decoy-prefix` | string | `XXX_` | Accession prefix marking a decoy. |
| `--decoy-suffix` | string | *(off)* | Decoy accession suffix (OpenMS `_rev` convention). |
| `--decoy-seed` | u64 | fixed | *(advanced)* RNG seed for `shuffle` decoys. |

### Chimeric cascade

Needs MS1: **mzML or Thermo `.raw`** only (MGF/`.d` warn and run a normal search).

| Flag | Type | Default | Description |
|---|---|---|---|
| `--chimeric` | flag | *(off)* | Pass 2 finds co-isolated precursors in the MS1 isolation window (averagine match) and searches the residual spectrum for a second peptide. Forces top-1 per pass and MS2. Experimental. |
| `--chimeric-max-coisolated` | u32 | `4` | *(advanced)* Max co-isolated precursors per scan. |
| `--chimeric-max-kl` | f64 | `0.3` | *(advanced)* Max isotope-envelope KL divergence to accept a co-isolated precursor. |

### Refine — secondary chemistry cascade

A second pass over confident proteins, on the spectra pass 1 missed, with a **fixed tier**:
oxidation on M/P/K, deamidation on N/Q, the two pyro-Glu losses and protein N-terminal acetyl,
at most two per peptide. It does **not** discover modifications; andes has no open search.

| Flag | Type | Default | Description |
|---|---|---|---|
| `--refine` | flag | *(off)* | Enable the PTM-refinement cascade (Pass-2 over confident proteins). |
| `--refine-config` | path | *(tier default)* | *(advanced)* YAML tier config: the mod set, max variable mods, the high-res-only gate and `entrapment: true` (one shuffled entrapment anchor per Pass-2 anchor, to count Pass-2 false discoveries; diagnostic only). |
| `--refine-select-psm-fdr` | fraction | `0.01` | *(advanced)* PSM-FDR of the confident set that seeds Pass-2. |

**FDR for refined PSMs.** Pass-2 rows (`IsRefinement = 1`) share the PIN. Under one merged 1%
threshold they measured **~2.8% FDP** on Astral; thresholded on their own, 0.8–1.1%. Threshold
them separately when you need 1% on modified PSMs.

### Rescoring & FDR filtering

andes **does not compute FDR itself**; these flags run a rescorer in-process. Leave them off
when a pipeline (e.g. quantms) owns rescoring. `--fdr`/`--pep` alone are ignored with a warning.

| Flag | Type | Default | Description |
|---|---|---|---|
| `--rescore` | flag | *(off)* | Run Percolator and write rescored, filtered output. |
| `--rescore-native` | flag | *(off)* | *(advanced)* Built-in GBDT rescorer (3-fold CV folded by spectrum); non-production fallback. |
| `--fdr` | fraction | `0.01` | q-value threshold applied by a rescoring run. |
| `--pep` | fraction | *(off)* | *(advanced)* PEP threshold, ANDed with `--fdr`. |
| `--percolator-bin` / `--percolator-docker` / `--percolator-image` / `--percolator-args` | — | auto | *(advanced)* Backend: `--percolator-bin`, else `$PATH`, else Docker; `--percolator-args` passes flags through. |
| `--keep-pin` | bool | `true` | *(advanced)* Keep the intermediate PIN after rescoring. |

### Glycopeptide search

See [§9](#9-glycopeptide-search-experimental--advanced-knobs).

### Isobaric labeling (TMT / iTRAQ)

See [§7](#7-isobaric-labeling).

---

## 1b. Configuration file (`--config`)

`andes --config run.yaml` sets any parameter from YAML, in sections `io`, `search`, `scoring`,
`decoys`, `chimeric`, `refine`, `rescoring`, `glyco` (template:
[`config.example.yaml`](config.example.yaml)). Keys are optional; **`CLI flag > --config >
built-in default`**; unknown keys are a hard error. Values use the CLI strings
(`precursor_tol: 20ppm`, `charge: "2..5"`, `isotope_error: "-1..2"`, `enzyme: gluc,trypsin`).

```yaml
# run.yaml (minimal)
io:
  spectrum: [sample.mzML]
  database: human.fasta
  output_pin: out.pin
search:
  precursor_tol: 20ppm
  enzyme: trypsin
glyco:
  enabled: true
```

---

## 2. Mods.txt format

The Java MS-GF+ format (`crates/model/src/modification.rs`). Each line has five fields:

```text
<mass>,<aa>,<fix|opt>,<location>,<name>
```

| Field | Rule |
|---|---|
| `<mass>` | Monoisotopic mass delta in Da. Composition strings (`C2H3N1O1`) are **not** supported. |
| `<aa>` | One uppercase letter or `*`. `STY` is **not** supported; use one line per residue. |
| `<fix\|opt>` | `fix` = fixed, `opt` = variable. Case-insensitive. |
| `<location>` | `any`, `N-term`, `C-term`, `Prot-N-term`, `Prot-C-term` (case-insensitive; hyphens optional). |
| `<name>` | Name used in logs. |

`NumMods=N` sets the maximum variable mods per peptide (default `3`). `#` starts a comment. A
fixed and a variable mod on the same `(residue, location)` is rejected.

### Example (a) — Carbamidomethyl C + Oxidation M

```text
NumMods=3
57.02146,C,fix,any,Carbamidomethyl
15.99491,M,opt,any,Oxidation
```

These two are the built-in defaults when `--mods` is omitted.

### Example (b) — TMT 10-plex on K and peptide N-term

```text
NumMods=2
57.02146,C,fix,any,Carbamidomethyl
229.162932,K,fix,any,TMT10plex
229.162932,*,fix,N-term,TMT10plex
```

Add `--protocol TMT` to select `hcd_qexactive_tryp_tmt` (§7).

### Example (c) — Phosphorylation on S, T, Y

```text
NumMods=3
57.02146,C,fix,any,Carbamidomethyl
79.966331,S,opt,any,Phospho
79.966331,T,opt,any,Phospho
79.966331,Y,opt,any,Phospho
```

Add `--protocol phospho` to prefer a phospho model (e.g. `hcd_qexactive_tryp_phosphorylation`).

---

## 3. Output formats

Percolator `.pin` (always), plus optional `.tsv` and QPX parquet. Implementation:
`crates/output/src/pin.rs`, `crates/output/src/tsv.rs`.

### 3a. PIN columns

Tab-separated, one row per PSM, best-first within each spectrum by `RankScore`. The `chargeN`
one-hots follow `--charge`; the default 2–5 gives the 66 columns below. Do not confuse
**`RankScore`** (col 7, orders candidates; historically `RawScore`) with **`RawScore`** (col 62,
the fused strong score Percolator weights most; historically `StrongScore`; it also ranks under
`--score strong`).

| # | Column | Type | Range | Description |
|---|---|---|---|---|
| 1 | `SpecId` | string | — | `{specID}_{scan}_{rank}`; multi-row scans add `_{rowIdx}`. |
| 2 | `Label` | int | {−1, +1} | `+1` target, `−1` decoy (by source protein). |
| 3 | `ScanNr` | int | ≥0 | MS2 scan number. |
| 4 | `ExpMass` | float | >0 | Experimental neutral precursor mass (Da). |
| 5 | `CalcMass` | float | >0 | Theoretical neutral peptide mass (Da). |
| 6 | `mass` | float | >0 | Duplicate of `ExpMass`. |
| 7 | `RankScore` | int | unbounded | **Rank-LLR ranking score**. |
| 8 | `isotope_error` | int | [−1, 2] | Winning ¹³C isotope offset. |
| 9 | `peplen` | int | ≥6 | Residue count **+ 2** (flanks). |
| 10 | `dm` | float | signed | Precursor mass error (Da) after isotope correction. |
| 11 | `absdm` | float | ≥0 | `\|dm\|`. |
| 12–15 | `charge2`…`charge5` | 0/1 | one-hot | One-hot precursor charge. |
| 16 | `enzN` | 0/1 | one-hot | N-terminus fits the enzyme rule. |
| 17 | `enzC` | 0/1 | one-hot | C-terminus fits the enzyme rule. |
| 18 | `enzInt` | int | ≥0 | Internal cleavage sites. |
| 19 | `NumMatchedMainIons` | int | [0, peplen−1] | Matched charge-1 b/y fragment positions. |
| 20 | `longest_b` | int | [0, peplen−1] | Longest contiguous matched b-ion run. |
| 21 | `longest_y` | int | [0, peplen−1] | Longest contiguous matched y-ion run. |
| 22 | `longest_y_pct` | float | [0, 1] | `longest_y / peplen`. |
| 23 | `ExplainedIonCurrentRatio` | float | [0, 1] | Matched b+y intensity / total MS2 ion current. |
| 24 | `NTermIonCurrentRatio` | float | [0, 1] | Matched b-ion intensity / total MS2 ion current. |
| 25 | `CTermIonCurrentRatio` | float | [0, 1] | Matched y-ion intensity / total MS2 ion current. |
| 26 | `MS2IonCurrent` | float | ≥0 | Sum of all MS2 peak intensities (not log-scaled). |
| 27 | `IsolationWindowEfficiency` | float | 0.0 | Always `0.0`. |
| 28 | `MeanErrorTop7` | float | ≥0 | Mean absolute ppm error, top-7 matched ions. |
| 29 | `StdevErrorTop7` | float | ≥0 | Population stdev of absolute ppm errors (top-7). |
| 30 | `MeanRelErrorTop7` | float | signed | Mean signed ppm error (top-7). |
| 31 | `StdevRelErrorTop7` | float | ≥0 | Population stdev of signed ppm errors (top-7). |
| 32 | `matchedIonRatio` | float | [0, 1] | `NumMatchedMainIons / peplen`. |
| 33 | `EdgeScore` | int | unbounded | Per-bond edge-score sum (ion-existence + error); additive. |
| 34 | `PrecursorIsotopeKL` | float | ≥0 | Precursor envelope KL vs averagine. **Always 0.0** (kept for positional consumers). |
| 35 | `PrecursorSNR` | float | ≥0 | Precursor SNR from the MS1 envelope. **0.0 unless `--chimeric`.** |
| 36 | `DeltaRankScore` | float | ≥0 | `RankScore(best) − RankScore(2nd-best distinct peptide)`; rank-1 row only, else 0.0. |
| 37 | `TailorScore` | float | ≥0 | `RankScore ÷` spectrum's top-1% quantile; cross-spectrum comparability. |
| 38 | `PpmGaussianScore` | float | ≥0 | `Σ exp(−½(ppm/7)²)` over matched ions. |
| 39 | `NeutralLossIonCount` | int | ≥0 | Matched b/y ions with −H₂O/−NH₃ partner peaks. |
| 40 | `LongestComplementaryLadder` | int | [0, peplen−1] | Longest run of bonds where both bᵢ and y₍ₙ₋ᵢ₎ matched. |
| 41 | `ComplementaryIonBalance` | float | ≥0 | `Σ 1/(1+\|rankᵦ−rankᵧ\|)` over complementary bonds. |
| 42 | `MeanMatchedIntensityRank` | float | ≥1 | Mean intensity rank of matched ions (1 = most intense). |
| 43 | `DoublyChargedMatchedIonCount` | int | ≥0 | Matched charge-2 b/y ions. |
| 44 | `UniqueMatchFraction` | float | [0, 1] | Within-peptide peak-explanation uniqueness. |
| 45 | `ChanceMatchSurprise` | float | ≥0 | `Σ max(0, −ln(ρ·Δ))`: improbability of chance matches. |
| 46 | `IntensitySignal` | float | [0, 1] | Predicted vs observed intensity cosine. **0.0 without an intensity model.** |
| 47 | `FragPredExplained` | float | [0, 1] | `Σ(matched·pred)/Σpred`. **0.0 without a frag-intensity model.** |
| 48 | `FragPredChanceLLR` | float | ≥0 | `Σ matched·pred·max(0,−ln p_chance)`. **0.0 without a frag-intensity model.** |
| 49 | `FragTopKObserved` | float | [0, 1] | Top-K predicted-most-intense ions observed. **0.0 without a frag-intensity model.** |
| 50 | `RichIonLLR` | float | unbounded | Decoy-aware per-annotated-ion LLR sum. **0.0 without a rich-ion model.** |
| 51 | `IsRefinement` | 0/1 | one-hot | 1 for a Pass-2 refinement PSM. **0 without `--refine`.** |
| 52 | `NumMods` | int | ≥0 | Variable modifications on the peptide. |
| 53 | `RefinementModClass` | int | [0, 99] | Mod-class id for subgroup-FDR grouping. **0 without `--refine`.** |
| 54 | `ModSiteShiftedMatched` | int | ≥0 | Matched mass-shifted b/y ions. **0 for unmodified peptides.** |
| 55 | `ModSiteShiftedFrac` | float | [0, 1] | Matched shifted ÷ total shifted ions. |
| 56 | `ModSiteIntensFrac` | float | [0, 1] | Shifted-ion intensity ÷ all matched-ion intensity. |
| 57 | `ModSiteLocalized` | 0/1 | one-hot | 1 if a bracketing ion pair localizes the mod. |
| 58 | `ModSiteDetCount` | int | ≥0 | Site-determining (bracketing) ions over all sites. |
| 59 | `MassCompetitionEvidence` | float | ≥0 | `Σ 1/(1+ambiguity+ρ)`, alternative-mass competition. |
| 60 | `CandidateRankEntropy` | float | ≥0 | Softmax entropy of the retained top-K scores. |
| 61 | `ListwiseScoreGap` | float | signed | Top-1 − top-2 `RankScore` in the retained queue. |
| 62 | `RawScore` | float | unbounded | **Fused strong score** `signal − null`; the primary feature. |
| 63 | `RawScoreCal` | float | signed | Per-spectrum z-scored `RawScore`. |
| 64 | `RankScoreFloat` | float | unbounded | Unrounded `RankScore` (continuous split-sum). |
| 65 | `Peptide` | string | — | `pre.SEQUENCE.post` with `+mass` mod annotations. |
| 66 | `Proteins` | string | — | Protein accession(s), tab-separated; decoys carry `--decoy-prefix`. |

Columns marked **0.0 unless/without** stay in the header and are zero when their condition
does not hold.

### 3b. TSV columns

**MGF header** (decoys included):

| Column | Type | Description |
|---|---|---|
| `#SpecFile` | string | Input file name. |
| `SpecID` | string | MGF title or `scan=N`. |
| `ScanNum` | int | Scan number. |
| `Title` | string | MGF `TITLE=` field. |
| `FragMethod` | string | Activation (`HCD`, `CID`, …) or `UNKNOWN`. |
| `Precursor` | float | Precursor m/z (4 decimal places). |
| `IsotopeError` | int | Same as PIN `isotope_error`. |
| `PrecursorError(ppm)` | float | Mass error; `PrecursorError(Da)` in Da mode. |
| `Charge` | int | Assigned precursor charge. |
| `Peptide` | string | Peptide with modifications. |
| `Protein` | string | Primary protein accession. |
| `RawScore` | int | Rounded raw score (the only score column). |

**mzML header**: the same without `Title` (11 columns).

### 3c. PIN vs TSV — which to use

**PIN** for FDR and rescoring (Percolator, MS²Rescore, Mokapot, quantms); **TSV** for inspection.

### 3d. Run summary (`statistics.log`)

Calibration and model choice can change the tolerances mid-run, so every search ends by
printing the **final** tolerances, the pre-FDR rank-1 target/decoy split and a
per-modification tally to stderr and to `statistics.log` next to the PIN.

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

### 3e. QPX `.idparquet` bundle (`--output-parquet`)

`--output-parquet <DIR>` writes an **OpenMS-compatible QPX 1.0** bundle (`psms.parquet`,
`proteins.parquet`, `search_params.parquet`) matching OpenMS's `QPXFile` schema, for OpenMS and
[quantms](https://github.com/bigbio/quantms). `psms.parquet` carries the peptidoform,
modifications, charge, m/z, `is_decoy`, scan/rt, proteins, the spectrum arrays, `score`
(`andes:RawScore`) and the other features in `additional_scores`. PEP and q-value are null
until rescoring; andes does no protein inference.

```bash
andes --spectrum spectra.mzML --database db.fasta \
  --output-pin out.pin --output-parquet out.idparquet
```

---

## 4. Auto-detection

With `--fragmentation auto`, andes takes the dominant `<activation>` cvParam over the first 64
MS2 spectra of an mzML (mixed methods warn) and the dominant analyser (none → `low-res`).
`.raw` carries both in vendor metadata (HCD on an Orbitrap → `hcd_qexactive_tryp`); `.d` routes
to `cid_tof_tryp`. **MGF** has none: andes assumes CID / low-res / 0.5 Da (`cid_lowres_tryp`)
and warns; `--fragment-tol-ppm` implies a high-res instrument, `--fragment-tol-da` a low-res
one. `--protocol` applies on top.

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

### Bundled model store (`resources/models/`)

17 models in `resources/models/protocol=<Automatic|TMT|Phosphorylation|iTRAQ>/models.parquet`
(9 / 3 / 4 / 1 models), listed in [`README.md`](README.md#supported-models). If detection
fails, andes falls back to `hcd_qexactive_tryp` or the closest regime and names it in the run
summary.

---

## 5. Building from source

Rust **1.85+** (`rust-toolchain.toml` pins **1.87.0**).

```bash
git clone https://github.com/bigbio/andes
cd andes
cargo build --release
# Binary: target/release/andes   (mzML + MGF; pure Rust)
```

**Native vendor formats** are feature-gated:

```bash
# Thermo .raw — needs rustc >= 1.88 and, at run time, the .NET 8 runtime (point DOTNET_ROOT
# at an install that contains shared/Microsoft.NETCore.App, or "One of the dependent
# libraries is missing" is the error you get).
RUSTUP_TOOLCHAIN=stable cargo build --release -p andes --features thermo

# Bruker timsTOF .d — pure Rust, no vendor runtime
cargo build --release -p andes --features timstof

# Both at once (what the release archives ship for desktop/server targets)
RUSTUP_TOOLCHAIN=stable cargo build --release -p andes --features "thermo timstof"
```

Tests:

```bash
cargo test --release --workspace
```

CI skips seven tests: three hit a `min_peaks` filter regression, three need Maven fixtures
under `target/test-classes/`, and one hits Rayon tie-breaking nondeterminism:

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

---

## 6. Training new scoring models

`andes train-from-search` searches your data with a seed model, keeps PSMs at q ≤
`--train-fdr`, and writes a model into a Parquet store; search with it via
`--model-store <path> --model <id>`. Incremental updates (`--update --add` /
`--remove-source` / `--reweight` / `--decay`) pass a held-out acceptance gate. `andes train
--in <parquet>` trains from externally labelled PSMs instead. See **[`TRAIN.md`](TRAIN.md)**.

---

## 7. Isobaric labeling

Set `--protocol TMT` or `--protocol iTRAQ` (selects `hcd_qexactive_tryp_tmt` /
`hcd_qexactive_tryp_itraq`) and declare the label as a fixed mod.

### TMT (10-plex example)

TMT10plex = **229.162932 Da** on K and peptide N-terminus (Unimod); mods.txt as in §2 example (b).

```bash
andes \
  --spectrum tmt_spectra.mzML \
  --database hsapiens.fasta \
  --output-pin out.pin \
  --mods tmt_10plex_mods.txt \
  --protocol TMT
```

### iTRAQ (8-plex example)

iTRAQ8plex = **304.20536 Da** on K and peptide N-terminus.

```text
NumMods=2
57.02146,C,fix,any,Carbamidomethyl
304.20536,K,fix,any,iTRAQ8plex
304.20536,*,fix,N-term,iTRAQ8plex
```

```bash
andes \
  --spectrum itraq_spectra.mzML \
  --database hsapiens.fasta \
  --output-pin out.pin \
  --mods itraq_8plex_mods.txt \
  --protocol iTRAQ
```

For phospho-enriched iTRAQ use `--protocol iTRAQ-phospho` plus phospho mods (§2 example c).

---

## 8. Legacy numeric values & behavior notes

Legacy MS-GF+ numeric values (e.g. `--fragmentation 3`, `--protocol 4`) are **no longer
accepted**; use the names (case-insensitive, `--fragmentation hcd` ≡ `HCD`). The MS-GF+
`-ntt` setting is `--enzyme-specificity fully|semi|non-specific`.

### Behavior notes

- No mzIdentML output.
- Bundled models are trypsin-trained except the three low-res LysC/ArgC/GluC models.
- Spectra stream in chunks of 5000, so large mzML files are not loaded whole.

---

## 9. Glycopeptide search (experimental) & advanced knobs

`--glyco` searches intact N-glycopeptides (N-X-S/T sequon backbones) and writes a
`.glyco.pin`. The defaults are validated; most knobs are hidden.

**Main flags:**

| Flag | Default | Purpose |
|---|---|---|
| `--glyco-tol-ppm` | 20 | Fragment tolerance for glyco matching (oxonium, core-Y, backbone mass, c/z). **Raise it on ion-trap MS2** or the oxonium gate never fires. |
| `--glyco-glycan-gdb <FILE>` | — | pGlyco-style `.gdb` glycan database; structure (core- vs antenna-fucose) is kept. Wins over `--glyco-species`. |
| `--glyco-species <NAME>` | — | Bundled database: `human`, `human-multi`, `mouse`, `mouse-large`, `high-mannose`. One of the two is required. |
| `--precursor-mono` | `auto` | Correct each precursor to the monoisotope its MS1 envelope supports (mzML or `.raw`): tests "recorded = M+k", k = 0..6, and moves down k−1 isotopes when k > 0 clearly wins. On pGlyco2 mouse liver 88 of 3,824 reference scans were recorded exactly +4 high (issue #64). Fitted spectra are searched with a `0..1` isotope window instead of `0..2`. Adds `MonoShift`/`MonoFit`/`MonoFitGain`/`MonoSNR`. An explicit `--isotope-error` is honoured; without MS1 output is byte-identical to `off`. Under `--glyco` the isotope-error window defaults to 0..=2 (dropping −1 measured +81 backbone-correct @1%). |
| `--glyco-max-peaks` | 0 (no cap) | Peaks the **generation** stage considers; 300–500 helps very dense scans. |
| `--glyco-retrieval-tol-ppm` | tol-ppm on high-res, the model's Da window on low-res | Retrieval window only; 20 ppm on high-res measured 6.9x faster at no identification cost. |
| `--glyco-y-max-charge` | 3 | Maximum glycan-Y fragment charge. |
| `--glyco-cz-max-charge` | derived | Maximum c/z charge on ETD. |
| `--glyco-hcd-pair` | **on** | ETD, single-file runs: backbones from the paired HCD scan, c/z scored on the ETD scan (+153 backbone-correct @1%). Disabled with a warning for multi-file runs. |
| `--glyco-min-core-y` | 0 | Require N trimannosyl-core Y ions before reporting. |
| `--glyco-scans <FILE>` / `--debug-glyco` | off | Diagnostics. Never feed a `--debug-glyco` PIN to an FDR tool. |

**Tuning knobs** (`hide = true`):

| Flag | Default | Purpose |
|---|---|---|
| `--glyco-index-sequon-only` | off | Index only N-X-S/T peptides: the mouse entrapment recipe needs 3.4 GB instead of ~27 GB; 16 of 7,113 rows (0.2%) differ in `RawScore`/`CandidateRankEntropy` only. |
| `--precursor-mono-dump <FILE>` | off | Write every `--precursor-mono` envelope fit as TSV. |
| `--glyco-peptide-first` | off | Peptide-first b/y fragment-index fallback instead of the default mass-driven full-glycan-list branch: ~2% more glycoPSMs on a mouse-liver fraction at ~13x the run time (79 vs 6 min) and 29.7 vs 19.9 GB peak memory. |

**Fixed settings** (validated; no flag): selector `rank + K·ladder + J·core_y + H·hyper` with
K/J/H = 10/5/1, ETD c/z weight 15 and isotope-offset penalty 1; 150 backbone candidates per
spectrum after the DB/de-novo union; ETD c/z evidence can rescue a backbone from truncation;
on ETD the rank/edge/hyperscore path scores the backbone with its intact glycan (+33
backbone-correct @1%); the best enumerated candidate is promoted over a de-novo winner;
peptide-first index at fragment charges 1..=2 with at most 1,024 candidates per spectrum;
`--precursor-mono` uses max shift 6, minimum fit 0.90, minimum gain 0.15, minimum SNR 3.0,
10 ppm MS1 tolerance and holds back 1 isotope.

Deleted after A/B tests (2026-09-05): `--glyco-split-election`, `--glyco-gp-g`,
`--glyco-gp-m`, `--glyco-isobar-rep`, `--glyco-y-index`, `--glyco-decorated-features`,
`--glyco-cz-intensity`, `--glyco-per-spectrum-model`, `--tight-highres-scoring`,
`--glyco-y-tree`, `--glyco-oxonium-llr`, `--glyco-rank-masked`, `--glyco-chance-llr-masked`,
`--glyco-transfer` (and its five knobs) and `--glyco-decoy`
([*Refuted*](docs/benchmarks/README.md#refuted--do-not-re-try-without-new-evidence)).
Removed 2026-10 (defaults fixed as above, or diagnostics/opt-ins that were never adopted):
`--glyco-gp-k/-j/-h/-cz/-iso`, `--glyco-backbone-top-k`, `--glyco-pf-charge`,
`--glyco-max-pf`, `--glyco-retrieval-tol-da`, `--glyco-isotope-error`,
`--glyco-elect-top-k`, `--glyco-cz-gate`, `--glyco-cz-multisite`, `--glyco-etd-rank-glycan`,
`--glyco-enum-fallback`, `--glyco-pair-y-on-gen`, `--glyco-etd-require-oxonium`,
`--glyco-min-raw-score(-quantile)`, `--glyco-min-matched-ions`,
`--glyco-sialic-oxonium-min-frac`, `--glyco-pin-curated`, `--glyco-diag-splits` and the
`--precursor-mono-*` tuning flags.

**Run Percolator on a glyco PIN with `--trainFDR 0.05`.** A glyco run has too few positives
at the default 1% training threshold for a stable fit. Pooled human-plasma glyco PIN, same
rows, five seeds each:

| `--trainFDR` | PSMs @1% | agreeing with the reference identifications |
|---|---:|---:|
| 0.01 (Percolator's default) | 210.6 ± 50.9 | 171.0 ± 37.6 |
| **0.05** | **384.6 ± 19.9** | **301.4 ± 11.1** |
| 0.10 | 348.8 ± 9.6 | 283.2 ± 6.1 |

The default is worse and erratic (112 vs 256 PSMs across seeds on identical input). Do **not**
filter low-scoring rows out of the PIN instead: they are about half decoys and Percolator
needs them to place a threshold.

---

## 10. License and citation

andes is licensed under the **Apache License 2.0**. See [`LICENSE`](LICENSE) for the full text and [`NOTICE`](NOTICE) for attribution and the project's origin in MS-GF+. The software is provided **"as is"** without warranty.

### Citation

If you use andes in published work, please cite both andes and the foundational MS-GF+ paper:

> bigbio (2026). andes: a data-driven peptide search engine for the quantms ecosystem. https://github.com/bigbio/andes

> Kim, S. and Pevzner, P.A. (2014). MS-GF+ makes progress towards a universal database search tool for proteomics. *Nature Communications*, 5:5277.

andes originated from MS-GF+ (https://github.com/MSGFPlus/msgfplus); see [`NOTICE`](NOTICE).
