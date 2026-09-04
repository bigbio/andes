# Reproducing the benchmarks

Three scripts, no hardcoded paths, nothing specific to our machines. You need `andes` on
`PATH` (or `$ANDES`), Docker for Percolator, `curl` and `python3`.

```bash
DATA=~/andes-bench                     # anywhere with ~200 GB free for the raw files

./build_databases.sh "$DATA"           # ~40 MB from UniProt, about a minute
./fetch_spectra.sh   "$DATA" tmt ups1  # spectra from PRIDE
./run.sh             "$DATA" tmt ups1  # search + Percolator + results table
```

Add `glyco-mouse` or `glyco-plasma` to `fetch_spectra.sh` for the glyco datasets; those are
run through the harness in [`../glyco/`](../glyco/)
rather than `run.sh`, because they must be pooled before Percolator.

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

## ⚠ PIN YOUR RAW CONVERTER — it changes the answer by 30%

Measured 2026-09-04 on one file (`MouseLiver-Z-T-1.raw`, 72,558 MS scans in the raw):

| ThermoRawFileParser | spectra written | MS2 | glycoPSMs @1% |
|---|---:|---:|---:|
| **1.4.3** | ~56,000 | **45,905** | **6,554** |
| 2.0.0 | 42,774 | 33,893 | 4,605 |

Same raw file, same search, same database, same flags (`-f 2`). **v2.0.0 dropped roughly
30,000 scans**, 26% of the MS2, and with them 30% of the identifications. This was found
while chasing what looked like a code regression between two andes builds; it was the
converter all along.

So: **pin the converter version and report it with any number**. The figures published here
were produced with **ThermoRawFileParser 1.4.3**. A benchmark run with a different converter
is not comparable to them, and the difference is far larger than any engine change we have
measured this year.

## Two things that will not reproduce exactly, stated up front

**The Astral row cannot be fetched at all.** Our records name
`LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.raw` in PXD070049, but that project does not
contain that file — verified against its complete 100-file listing on 2026-09-04. Either
the accession or the filename in our records is wrong. `fetch_spectra.sh` refuses that
dataset and says so rather than downloading nothing quietly. The other datasets are
unaffected.

**UniProt is versioned, so databases are not byte-identical to ours.** A build on
2026-09-04 gave:

| database | this build | ours | note |
|---|---:|---:|---|
| `tmt_db.fasta` | 26,483 | 26,483 | exact match |
| `hye.fasta` | 30,886 | 31,889 | ProteoBench's own file is not fetchable (no listing); this is a Human/Yeast/E.coli reconstruction |
| `yeast_entrap.fasta` | 10,470 | 11,264 | different UniProt release |

Counts will therefore differ slightly from the published table. That is expected and is why
every script prints the sha256 and sequence count of what it actually built — quote those
alongside any number you report.

## Reading the output

`PSMs@q0.01` is Percolator's *claim*, comparable across engines run through this same
protocol. A measured error rate needs an entrapment database:

```
FDP = (entrapment hits / total accepted) x (1 + T/E)
```

with `T/E` from `build_databases.sh`. **Never assume 1:1 (factor 2)** — doing so has
understated the true error in this project twice. Only the UPS1 database here has an
entrapment component, so only it yields a true FDP.
