# Reproducing the benchmarks

Three scripts, no hardcoded paths, nothing specific to our machines. You need `andes` on
`PATH` (or `$ANDES`), Docker for Percolator, `curl` and `python3`.

```bash
DATA=~/andes-bench                     # anywhere with ~200 GB free for the raw files

./build_databases.sh "$DATA"           # ~40 MB from UniProt, about a minute
./fetch_spectra.sh   "$DATA" tmt ups1  # spectra from PRIDE
./run.sh             "$DATA" tmt ups1  # search + Percolator + results table
```

Add `glyco-mouse` (12.6 GB; the quick tier uses only the first fraction) to
`fetch_spectra.sh` for the glyco dataset; it is run through the harness in
[`../glyco/`](../glyco/) rather than `run.sh`, because fractions must be pooled before Percolator.

## What each script does

**`build_databases.sh`** pulls the human, yeast and E. coli reviewed proteomes from
UniProt's REST API and assembles the three search databases, including the UPS1 entrapment
database (yeast targets + `ENTRAP_`-tagged E. coli). It then **measures and prints the
entrapment scaling factor `T/E` for the database it just built** — you need that number to
compute a true FDP, and it is not 1.

**`fetch_spectra.sh`** resolves real download URLs through the PRIDE API rather than
hardcoding archive paths, so it survives layout changes. Files already present are skipped,
and a sha256 is printed for each so you can confirm you have the same inputs we did.

**`run.sh`** searches each dataset, rescores through the pinned Percolator container, and
prints wall time, PSMs at `q ≤ 0.01`, and entrapment hits — with a provenance line naming
the binary, platform, thread count and date.

## ⚠ Read `.raw` natively — do not convert

andes reads Thermo `.raw` directly (`--features thermo`, plus the .NET 8 runtime at search
time). **Use it.** Conversion is an extra step that can silently change your results, and
native reading costs nothing:

| input (a human-plasma glyco file, since retired from the benchmark) | andes spectra | PIN rows | wall |
|---|---:|---:|---:|
| `.raw` native | 24,857 | 7,257 | 152 s |
| TRFP 1.4.3 mzML | 24,857 | 7,257 | 149 s |

Identical output, same speed, one less moving part.

### Why this matters: a converter version silently changed results by 30%

Verified against native reading as the reference — native uses Thermo's own RawFileReader,
so it is the ground truth for what a file contains:

| `MouseLiver-Z-T-1.raw` | MS2 | glyco rows |
|---|---:|---:|
| **native (reference)** | **45,905** | **41,929** |
| TRFP 1.4.3 | 45,905 | 41,929 |
| TRFP 2.0.0 | 33,892 | 31,279 |

**TRFP 1.4.3 is correct; 2.0.0 dropped 26% of the MS2 scans on this file**, and with them
30% of the identifications. This surfaced only because the same file was searched on two
machines and the counts disagreed.

**It is file-dependent, not a blanket property of 2.0.0.** On the plasma file above, both
converter versions produced identical output. So you cannot assume a given converter is
safe for your data — which is the argument for skipping conversion altogether.

If you must convert, use **1.4.3**, and state the converter version with any number you
publish. All figures here were produced with native reading or 1.4.3, which agree exactly.

## What will not reproduce byte-for-byte, stated up front

**Databases: UniProt is versioned.** A build on 2026-09-04 gave:

| database | this build | ours | note |
|---|---:|---:|---|
| `tmt_db.fasta` | 26,483 | 26,483 | exact match |
| `hye.fasta` | 30,886 | 31,889 | ours was ProteoBench's own file (sha256 `d9ac434d…`); every URL it was served from now returns 404, so this is a Human/Yeast/E.coli reconstruction |
| `yeast_entrap.fasta` | 10,470 | 11,264 | different UniProt release |

Counts will therefore differ slightly from the published table. That is expected and is why
every script prints the sha256 and sequence count of what it actually built — quote those
alongside any number you report.

**Spectra: all five datasets fetch.** An earlier version of this document said the Astral
file could not be found in PXD070049. That was this script's own bug: the PRIDE API serves
at most 100 files per page, the project has 2,173, and the script read one page. It now
pages through the whole listing, and every file the benchmarks use resolved on 2026-09-05:

| dataset | accession | files | size |
|---|---|---|---:|
| astral | PXD070049 | `LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.raw` | 2.58 GB |
| tmt | PXD007683 | `a05058.raw` | 0.54 GB |
| ups1 | PXD001819 | `UPS1_5000amol_R1.raw` | 1.70 GB |
| glyco-mouse | PXD005553 | `MouseLiver-Z-T-{1..5}.raw` (pGlyco2 liver) | 12.6 GB |

## Reading the output

`PSMs@q0.01` is Percolator's *claim*, comparable across engines run through this same
protocol. A measured error rate needs an entrapment database:

```
FDP = (entrapment hits / total accepted) x (1 + T/E)
```

with `T/E` from `build_databases.sh`. **Never assume 1:1 (factor 2)** — doing so has
understated the true error in this project twice. Of the standard databases only UPS1 has an
entrapment component, so only it yields a true FDP; the glyco mouse database is 1:1 by
construction (factor exactly 2).
