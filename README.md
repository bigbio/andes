# MS-GF+ (bigbio fork for quantms)

[![CI](https://github.com/bigbio/msgfplus/actions/workflows/ci.yml/badge.svg)](https://github.com/bigbio/msgfplus/actions/workflows/ci.yml)

> **This is a lightweight fork of [MS-GF+](https://github.com/MSGFPlus/msgfplus) maintained by [bigbio](https://github.com/bigbio) for use in the [quantms](https://github.com/bigbio/quantms) pipeline.** It contains targeted performance improvements (streaming mzML parsing, reduced memory footprint) and CI/release automation. The primary maintained input formats in this fork are mzML and MGF; legacy text-based readers (including mzXML) remain available for compatibility.
>
> **For the full-featured, officially maintained version of MS-GF+** and the latest upstream features, please use the original repository:
>
> **[https://github.com/MSGFPlus/msgfplus](https://github.com/MSGFPlus/msgfplus)**

## What is MS-GF+?

MS-GF+ (aka MSGF+ or MSGFPlus) performs peptide identification by scoring
MS/MS spectra against peptides derived from a protein sequence database.
It supports the HUPO PSI standard input file (mzML) and additional legacy spectrum inputs, and saves results in
the mzIdentML format, though results can easily be transformed to TSV.
ProteomeXchange supports Complete data submissions using MS-GF+ search results.

MS-GF+ is developed by Sangtae Kim and the PNNL Proteomics team at the
Center for Computational Mass Spectrometry, University of California, San Diego.

## What is different in this fork?

- **Streaming mzML parser** -- replaces the in-memory preload with a single-pass StAX parser, significantly reducing memory usage for large files
- **Primary maintained formats: mzML and MGF** -- legacy formats (including mzXML) are still available, but not the main optimization target in this fork
- **Java 17 minimum** -- updated from Java 8
- **CI/CD** -- GitHub Actions for automated testing and releases
- **Direct TSV output** -- optional TSV output alongside mzIdentML

## Requirements

- Java Runtime 17 or higher (use 64-bit Java)
- At least 2 GB of memory (4 GB+ recommended); larger FASTA files require more memory

## Installation

Download the latest release from the [Releases page](https://github.com/bigbio/msgfplus/releases). The zip contains `MSGFPlus.jar` with all dependencies bundled.

## Quick Start

```bash
# Basic search
java -Xmx4G -jar MSGFPlus.jar \
  -s spectra.mzML \
  -d database.fasta \
  -o results.mzid

# TMT search with target-decoy analysis
java -Xmx8G -jar MSGFPlus.jar \
  -s spectra.mzML \
  -d database.fasta \
  -tda 1 \
  -t 20ppm \
  -ti -1,2 \
  -inst 1 \
  -e 1 \
  -protocol 4 \
  -mod mods.txt \
  -o results.mzid

# Convert mzid output to TSV
java -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv \
  -i results.mzid \
  -o results.tsv
```

## Parameters

### Required

| Flag | Name | Description |
|------|------|-------------|
| `-s` | SpectrumFile | Input spectrum file (`*.mzML`, `*.mzXML`, `*.mgf`, `*.ms2`, `*.pkl`, `*_dta.txt`). Spectra should be centroided. |
| `-d` | DatabaseFile | Protein sequence database (`*.fasta`, `*.fa`, `*.faa`). |

### Core Search Parameters

| Flag | Name | Default | Description |
|------|------|---------|-------------|
| `-o` | OutputFile | `[input].mzid` | Output file path (`.mzid` format). |
| `-conf` | ConfigurationFile | — | Configuration file; command-line options override config file settings. |
| `-t` | PrecursorMassTolerance | `20ppm` | Precursor mass tolerance (e.g., `2.5Da`, `20ppm`, or `0.5Da,2.5Da` for asymmetric). |
| `-ti` | IsotopeErrorRange | `0,1` | Range of allowed isotope peak errors (e.g., `-1,2`). |
| `-tda` | TDA | `0` | Target-decoy analysis: `0` = don't search decoy database, `1` = search decoy database. |
| `-decoy` | DecoyPrefix | `XXX` | Prefix for decoy protein names in the FASTA file. |

### Fragmentation and Instrument

| Flag | Name | Default | Description |
|------|------|---------|-------------|
| `-m` | FragmentationMethodID | `0` | `0` = As written in spectrum or CID if no info, `1` = CID, `2` = ETD, `3` = HCD, `4` = UVPD. |
| `-inst` | InstrumentID | `0` | `0` = Low-res LCQ/LTQ, `1` = Orbitrap/FTICR/Lumos (default for HCD), `2` = TOF, `3` = Q-Exactive. |

### Enzyme and Digestion

| Flag | Name | Default | Description |
|------|------|---------|-------------|
| `-e` | EnzymeID | `1` | `0` = Unspecific, `1` = Trypsin, `2` = Chymotrypsin, `3` = Lys-C, `4` = Lys-N, `5` = Glu-C, `6` = Arg-C, `7` = Asp-N, `8` = alphaLP, `9` = No cleavage, `10` = TrypPlusC. |
| `-ntt` | NTT | `2` | Number of tolerable termini: `0` = non-specific, `1` = semi-specific, `2` = fully specific. |
| `-maxMissedCleavages` | MaxMissedCleavages | `-1` | Maximum missed cleavages (`-1` = no limit). |

### Peptide Filtering

| Flag | Name | Default | Description |
|------|------|---------|-------------|
| `-minLength` | MinPepLength | `6` | Minimum peptide length to consider. |
| `-maxLength` | MaxPepLength | `40` | Maximum peptide length to consider. |
| `-minCharge` | MinCharge | `2` | Minimum precursor charge (if not in spectrum file). |
| `-maxCharge` | MaxCharge | `3` | Maximum precursor charge (if not in spectrum file). |
| `-msLevel` | MSLevel | `2` | MS level(s) to search (e.g., `2` or `2,3` for MS2+MS3). |

### Modifications and Protocol

| Flag | Name | Default | Description |
|------|------|---------|-------------|
| `-mod` | ModificationFileName | — | Modification file path. If not specified, uses standard amino acids with fixed Carbamidomethyl C. |
| `-numMods` | NumMods | `3` | Maximum number of dynamic (variable) modifications per peptide. |
| `-protocol` | ProtocolID | `0` | `0` = Automatic, `1` = Phosphorylation, `2` = iTRAQ, `3` = iTRAQPhospho, `4` = TMT, `5` = Standard. |

### Output and Performance

| Flag | Name | Default | Description |
|------|------|---------|-------------|
| `-n` | NumMatchesPerSpec | `1` | Number of matches per spectrum to report. Values >1 may skew FDR. |
| `-addFeatures` | AddFeatures | `0` | `0` = basic scores, `1` = additional features (enable for Percolator). |
| `-thread` | NumThreads | All cores | Number of concurrent threads. |
| `-tasks` | NumTasks | `0` | Override task count: `0` = auto, `>0` = exact count, `<0` = multiplier of threads. |
| `-verbose` | Verbose | `0` | `0` = total progress only, `1` = per-thread progress. |
| `-ccm` | ChargeCarrierMass | `1.00727649` | Mass of charge carrier (proton). |

### Advanced Parameters

| Flag | Name | Default | Description |
|------|------|---------|-------------|
| `-minNumPeaks` | MinNumPeaksPerSpectrum | `10` | Minimum number of peaks per spectrum. |
| `-iso` | NumIsoforms | `128` | Number of isoforms to consider per peptide. |
| `-ignoreMetCleavage` | IgnoreMetCleavage | `0` | `0` = consider N-term Met cleavage, `1` = ignore. |
| `-allowDenseCentroidedPeaks` | AllowDenseCentroidedPeaks | `0` | `0` = skip spectra failing density check, `1` = allow dense centroided spectra. |

## Modification File Format

Modifications are specified in a text file passed via `-mod`. Each line defines a static or dynamic modification:

```
# Format: Mass_or_Composition, Residues, ModType, Position, Name

# Static modifications
StaticMod=C2H3N1O1,  C,   fix, any,       Carbamidomethyl   # Fixed alkylation
StaticMod=229.1629,   *,   fix, N-term,    TMT6plex
StaticMod=229.1629,   K,   fix, any,       TMT6plex

# Dynamic modifications
DynamicMod=O1,        M,   opt, any,       Oxidation         # Oxidized methionine
DynamicMod=HO3P,      STY, opt, any,       Phospho           # Phosphorylation
DynamicMod=H-1N-1O1,  NQ,  opt, any,       Deamidated        # Deamidation

# Position options: any, N-term, C-term, Prot-N-term, Prot-C-term
```

See [`docs/examples/MSGFPlus_Params.txt`](docs/examples/MSGFPlus_Params.txt) for a complete example configuration file, and [`docs/examples/README.md`](docs/examples/README.md) for what else lives in that folder. Long-form usage topics (MzID→TSV, BuildSA, changelog, and so on) live under [`docs/README.md`](docs/README.md).

## Building from Source

```bash
# Requires Java 17+ and Maven (same as CI)
mvn -B verify

# The shaded JAR is produced at target/MSGFPlus.jar
```

## Publications

Kim S. and Pevzner P.A.,
"MS-GF+ makes progress towards a universal database search tool for proteomics,"
*Nat Commun.* 2014 Oct 31; 5:5277.
[doi: 10.1038/ncomms6277](https://doi.org/10.1038/ncomms6277)

Kim S., Gupta N., and Pevzner P.A.,
"Spectral Probabilities and Generating Functions of Tandem Mass Spectra: A Strike against Decoy Databases,"
*J Proteome Res.* 2008 Aug; 7(8):3354-63.
[doi: 10.1021/pr8001244](https://doi.org/10.1021/pr8001244)

## Contact

For the official MS-GF+ tool: [MSGFPlus/msgfplus](https://github.com/MSGFPlus/msgfplus)

PNNL Proteomics: proteomics@pnnl.gov
Sangtae Kim: sangtae.kim (at) gmail.com

For this fork (quantms integration): [bigbio/msgfplus](https://github.com/bigbio/msgfplus)
