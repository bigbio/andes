# MS-GF+ Output Formats

[MS-GF+ Documentation home](readme.md) · [ChangeLog](changelog.md)

MS-GF+ writes two output formats: Percolator `.pin` (default) and a tab-separated `.tsv`. The `.mzid` format has been removed; downstream tools should rescore the `.pin` with Percolator or consume the `.tsv` directly.

Select the format with `-outputFormat`:

| Flag | Format | Extension | Typical use |
|---|---|---|---|
| `-outputFormat pin` (default) | Percolator `.pin` | `.pin` | Feed to Percolator / MS²Rescore / Mokapot for FDR-calibrated rescoring |
| `-outputFormat tsv` | Tab-separated values | `.tsv` | Direct inspection / downstream tools that consume TSV |

`-outputFormat` accepts the named values `pin` and `tsv` (case-insensitive). Numeric forms (`0`, `1`) accepted by older releases are no longer recognised — pass the named value instead.

The output path (`-o`) must use the matching extension. If `-o` is omitted, MS-GF+ writes `<SpectrumFileName>.pin` (or `.tsv`) in the spectrum file's directory.

---

## The `.pin` file (default)

A Percolator input file — tab-separated, one header row followed by one row per peptide-spectrum match (PSM). The column schema matches OpenMS `PercolatorAdapter` so the same pin can be consumed by Percolator, MS²Rescore, or Mokapot without translation.

### Row format

```
SpecId  Label  ScanNr  ExpMass  CalcMass  mass  RawScore  DeNovoScore  lnSpecEValue  lnEValue  isotope_error  peplen  dm  absdm  charge2 … chargeK  enzN  enzC  enzInt  NumMatchedMainIons  longest_b  longest_y  longest_y_pct  ExplainedIonCurrentRatio  NTermIonCurrentRatio  CTermIonCurrentRatio  MS2IonCurrent  IsolationWindowEfficiency  MeanErrorTop7  StdevErrorTop7  MeanRelErrorTop7  StdevRelErrorTop7  lnDeltaSpecEValue  matchedIonRatio  Peptide  Proteins
```

Rows are written in scoring order (rank-1 PSM per spectrum, unless `-n > 1`).

### Column reference

**Identity + label (positions 1–3)**

| Column | Meaning |
|---|---|
| `SpecId` | Unique PSM id: `<specID>_<scanNr>_<rank>`. SpecID comes from the spectrum file (nativeID for mzML). |
| `Label` | `1` for PSMs whose protein list contains at least one target, `-1` when every protein is a decoy. Percolator uses this to train the target/decoy classifier. |
| `ScanNr` | MS2 scan number from the spectrum file. |

**Mass + charge features (positions 4–16)**

| Column | Meaning |
|---|---|
| `ExpMass` | Experimental precursor mass (Da) = `precursorMz * charge`. |
| `CalcMass` | Theoretical peptide mass (Da) from the matched sequence. |
| `mass` | Duplicate of `ExpMass` — kept for OpenMS `PercolatorAdapter` layout parity. |
| `RawScore` | MS-GF+ `MSGFScore` (integer) — the sum-of-node-scores from the peptide's backbone graph. |
| `DeNovoScore` | MS-GF+ `DeNovoScore` — the best-path score through the spectrum graph ignoring peptide constraints. |
| `lnSpecEValue` | `log(SpecEValue)`. SpecEValue is MS-GF+'s exact E-value from the generating function. `-Double.MAX_VALUE` if non-positive. |
| `lnEValue` | `log(EValue)`. EValue is SpecEValue × number of distinct candidate peptides. |
| `isotope_error` | Rounded `(ExpMass − CalcMass) / 1.00335` — which isotope peak was fragmented (−2, −1, 0, 1, 2). |
| `peplen` | Peptide length *plus 2* (includes the flanking residues recorded in the suffix array). |
| `dm` | Mass delta `adjustedExpMz − theoMz` in Da (after isotope-error correction). |
| `absdm` | `abs(dm)`. |
| `charge2 … chargeK` | One-hot charge encoding over `params.getMinCharge()..getMaxCharge()`. |

**Enzymatic-boundary features (positions 17–19)**

| Column | Meaning |
|---|---|
| `enzN` | `1` if the peptide's N-terminal boundary is enzymatically valid (or a protein boundary); `0` otherwise. Rules mirror OpenMS `PercolatorInfile::isEnz_`. |
| `enzC` | `1` if the peptide's C-terminal boundary is enzymatically valid; `0` otherwise. |
| `enzInt` | Number of enzymatically cleavable positions *internal* to the peptide (≤ `peplen − 1`). For trypsin, typically `0` unless the PSM spans missed cleavages. |

**Ion-structure features (positions 20–24)**

| Column | Meaning |
|---|---|
| `NumMatchedMainIons` | Number of predicted main ions (b/y, any charge) matched to observed peaks. |
| `longest_b` | Longest consecutive run of matched b-ions along the backbone. |
| `longest_y` | Longest consecutive run of matched y-ions. Tends to be higher-weighted by Percolator than `longest_b` on tryptic data. |
| `longest_y_pct` | `longest_y / (peptide.size() − 1)` — fraction of peptide-bond positions covered by the longest y-ion run. |

**Ion-current features (positions 25–29)**

| Column | Meaning |
|---|---|
| `ExplainedIonCurrentRatio` | `NTermIonCurrentRatio + CTermIonCurrentRatio`. |
| `NTermIonCurrentRatio` | Summed intensity of matched b-ions, divided by total MS2 ion current for the spectrum. |
| `CTermIonCurrentRatio` | Same for y-ions. |
| `MS2IonCurrent` | Summed intensity of all observed peaks in the spectrum. |
| `IsolationWindowEfficiency` | Reserved (emitted as 0 in the current release). |

**Fragment-mass-error features (positions 30–33)**

Computed across the seven highest-intensity matched peaks in the spectrum.

| Column | Meaning |
|---|---|
| `MeanErrorTop7` | Mean ppm error. |
| `StdevErrorTop7` | Stdev of ppm error. |
| `MeanRelErrorTop7` | Mean relative error (intensity-weighted). |
| `StdevRelErrorTop7` | Stdev of relative error. |

`NaN` / `Infinity` values (produced when fewer than two ions are matched) are automatically rewritten to `0` so Percolator does not terminate on non-finite feature values.

**Generating-function / rank features (positions 34–35)**

| Column | Meaning |
|---|---|
| `lnDeltaSpecEValue` | `log(rank1_SpecEValue / rank2_SpecEValue)` for the rank-1 PSM of a spectrum; `0` otherwise or if the rank-2 candidate is missing. Larger (more negative) = more separation from the runner-up. |
| `matchedIonRatio` | `NumMatchedMainIons / peplen` — peptide-length-normalized ion-match density. |

**Peptide + protein (positions 36+)**

| Column | Meaning |
|---|---|
| `Peptide` | Percolator-style annotation: `flanking.PEPTIDE_WITH_MODS.flanking`. Modifications are written as `+mass` after the modified residue, or as `[+mass]-PEPTIDE` at peptide N-term. Fixed modifications ARE included. |
| `Proteins` | One column per matching protein accession (tab-separated). Prefixed with the decoy tag (default `XXX_`) when the match is a decoy. |

### Example file

A trimmed sample `.pin` from a PXD001819 search (20 target PSMs + 10 decoys, curated for peptide-sequence diversity so every feature column has data) is at [`examples/pxd001819_example.pin`](examples/pxd001819_example.pin). Open it in any TSV viewer to inspect the schema; the full MS-GF+ output on the same dataset has ~39 K rows.

### Using the `.pin` with Percolator

```bash
percolator --seed 42 \
  --results-psms out.target.psms.txt \
  --decoy-results-psms out.decoy.psms.txt \
  --weights out.weights.txt \
  --only-psms \
  run.pin
```

Percolator output is also tab-separated: `PSMId`, `score`, `q-value`, `posterior_error_prob`, `peptide`, `proteinIds`. Count PSMs at 1% FDR with `awk -F'\t' 'NR>1 && $3+0<=0.01 {c++} END {print c}' out.target.psms.txt`.

### Peptide modification notation

| Source | Example |
|---|---|
| Variable (on residue) | `K+229.163` (TMT6plex on K) |
| Variable (N-term) | `[+42.011]-PEPTIDE` (N-term acetyl) |
| Variable (C-term) | `PEPTIDE-[+0.984]` |
| Fixed (auto-applied to every occurrence) | rendered same as variable, e.g. `+57.021` after every C for carbamidomethyl |
| Stacked mods on one residue | `C+57.021+14.016` |

---

## The `.tsv` file

Tab-separated one-row-per-PSM table. Contains the same scoring columns as the pin plus an "#SpecFile" leading column identifying the source spectrum file. Useful when you want to inspect PSMs without running Percolator, or feed a downstream tool that accepts a flat PSM table.

Column order:

```
#SpecFile  SpecID  ScanNum  [Title]  FragMethod  Precursor  IsotopeError  PrecursorError(ppm)  Charge  Peptide  Protein  DeNovoScore  MSGFScore  SpecEValue  EValue  [QValue]  [PepQValue]  [additional features…]
```

- `Title` column is only present for MGF input.
- `QValue` / `PepQValue` only present when `-tda 1` (target-decoy search).
- Additional feature columns (same set as pin's extra features) are appended when `-addFeatures 1` is supplied.

---

## Migrating from `.mzid`

Previous MS-GF+ releases wrote `.mzid` (mzIdentML) by default. That output has been fully removed in the next release. To get equivalent FDR-calibrated PSM tables:

1. Run MS-GF+ with default settings (produces `.pin`).
2. Rescore the `.pin` with Percolator (or MS²Rescore / Mokapot).
3. Load the resulting `out.target.psms.txt` with your downstream tool.

For legacy `.mzid` post-processing, use MS-GF+ v2026.03.25 or earlier, or an external mzid toolchain.
