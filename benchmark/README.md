# MS-GF+ Benchmarks

Only the **CI benchmark scaffold** is committed under this directory; heavy
local-only harnesses and artifacts (`data/`, `results/`, prebuilt JARs, etc.)
are intentionally gitignored and not distributed with the fork.

## Datasets

| Dataset | PXD | Instrument | Type | Spectra Source | FASTA / SDRF |
|---------|-----|------------|------|----------------|--------------|
| LFQ (Yeast+UPS1) | PXD001819 | LTQ Orbitrap Velos, CID | LFQ | [PRIDE FTP](https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/benchmarks/lfq/LTQOrbitrapVelos/PXD001819/) | [quantms-test-datasets](https://github.com/bigbio/quantms-test-datasets) |

## CI benchmark scaffold

```
benchmark/
  README.md                      # this file
  ci/
    README.md                    # how to run/compare CI benchmark(s)
    PXD001819/
      run_ci.sh                  # downloads public data, runs MS-GF+, emits metrics
      compare_metrics.py         # compares metrics against baseline ranges
      baseline.tsv               # acceptable metric ranges
```

Per-dataset inputs (mod config, search-flag documentation) live under
`src/test/resources/benchmark/<PXD>/` so each dataset's configuration travels
with the repo while bulky downloads remain outside it.

## PXD001819 CI benchmark (public data)

End-to-end search on public PRIDE mzML + `quantms-test-datasets` FASTA, with
metrics compared to a checked-in baseline TSV. The workflow runs on a fixed
self-hosted runner profile (`self-hosted,linux,msgf-benchmark`) to keep
CPU/memory comparable between runs.

See [`benchmark/ci/README.md`](ci/README.md) for commands.
