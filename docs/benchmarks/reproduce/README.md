# Reproducing the benchmarks

Three scripts with no hardcoded paths. You need `andes` on `PATH` (or `$ANDES`), Docker for
Percolator, `curl` and `python3`.

```bash
DATA=~/andes-bench                     # anywhere with ~200 GB free for the raw files

./build_databases.sh "$DATA"           # ~40 MB from UniProt, about a minute
./fetch_spectra.sh   "$DATA" tmt ups1  # spectra from PRIDE
./run.sh             "$DATA" tmt ups1  # search + Percolator + results table
```

Add `glyco-mouse` (12.6 GB; the quick tier uses only the first fraction) to
`fetch_spectra.sh` for the glyco dataset, and run it through [`../glyco/`](../glyco/) rather
than `run.sh`: the deep tier pools the five fractions before Percolator, the quick tier runs
`MouseLiver-Z-T-1` alone.

## What each script does

- **`build_databases.sh`** builds the three search databases from UniProt's human, yeast and
  E. coli reviewed proteomes, including the UPS1 entrapment database (yeast targets +
  `ENTRAP_`-tagged E. coli), and **prints the measured `T/E` factor** you need for a true FDP.
- **`fetch_spectra.sh`** resolves download URLs through the PRIDE API, skips files already
  present, and prints a sha256 for each.
- **`run.sh`** searches, rescores through the pinned Percolator container, and prints wall
  time, PSMs at `q ≤ 0.01` and entrapment hits, with a provenance line (binary, platform,
  threads, date).

## Read `.raw` natively — do not convert

andes reads Thermo `.raw` directly (`--features thermo`, plus the .NET 8 runtime). Native
reading uses Thermo's own RawFileReader, so it is the reference for what a file contains:

| `MouseLiver-Z-T-1.raw` | MS2 | glyco rows |
|---|---:|---:|
| **native (reference)** | **45,905** | **41,929** |
| TRFP 1.4.3 | 45,905 | 41,929 |
| TRFP 2.0.0 | 33,892 | 31,279 |

TRFP 2.0.0 dropped 26% of the MS2 scans on this file, and 30% of the identifications with them.
On another file (a human-plasma glyco file, since retired) both versions agreed exactly, so the
effect is file-dependent and no converter can be assumed safe. If you must convert, use
**1.4.3** and state the version with any number you publish.

## What will not reproduce byte-for-byte

**Databases: UniProt is versioned.** A build on 2026-09-04 gave:

| database | this build | ours | note |
|---|---:|---:|---|
| `tmt_db.fasta` | 26,483 | 26,483 | exact match |
| `hye.fasta` | 30,886 | 31,889 | ours was ProteoBench's own file (sha256 `d9ac434d…`); every URL it was served from now returns 404, so this is a Human/Yeast/E.coli reconstruction |
| `yeast_entrap.fasta` | 10,470 | 11,264 | different UniProt release |

Counts will differ slightly; quote the sha256 and sequence count each script prints alongside
any number you report.

**Spectra: all datasets fetch.** The PRIDE API serves at most 100 files per page, and the
script pages through the whole listing (PXD070049 has 2,173 files). Every file resolved on
2026-09-05:

| dataset | accession | files | size |
|---|---|---|---:|
| astral | PXD070049 | `LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.raw` | 2.58 GB |
| tmt | PXD007683 | `a05058.raw` | 0.54 GB |
| ups1 | PXD001819 | `UPS1_5000amol_R1.raw` | 1.70 GB |
| glyco-mouse | PXD005553 | `MouseLiver-Z-T-{1..5}.raw` (pGlyco2 liver) | 12.6 GB |

## Reading the output

`PSMs@q0.01` is Percolator's *claim*, comparable across engines run through this protocol. A
measured error rate needs an entrapment database:

```
FDP = (entrapment hits / total accepted) x (1 + T/E)
```

with `T/E` from `build_databases.sh`; **never assume 1:1**. Of the standard databases only
UPS1 has an entrapment component; the glyco mouse database is 1:1 by construction (factor
exactly 2).
