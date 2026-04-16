# MS-GF+ Benchmark

Benchmarks comparing the bigbio fork against the official MS-GF+ baseline across
multiple datasets covering different instruments, fragmentation methods, and labeling strategies.

## Data Sources

| Dataset | PXD | Instrument | Type | Spectra Source | FASTA / SDRF |
|---------|-----|------------|------|----------------|--------------|
| TMT (Human) | internal | HCD | TMT | local | local |
| LFQ (Yeast+UPS1) | PXD001819 | LTQ Orbitrap Velos, CID | LFQ | [PRIDE FTP](https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/benchmarks/lfq/LTQOrbitrapVelos/PXD001819/) | [quantms-test-datasets](https://github.com/bigbio/quantms-test-datasets) |
| LFQ (Mixed-species) | PXD028735 | Q Exactive HF, HCD | LFQ | [PRIDE FTP](https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/benchmarks/) | [quantms-test-datasets](https://github.com/bigbio/quantms-test-datasets) |

## Directory Structure

```
benchmark/
  baseline/MSGFPlus.jar    # Official MS-GF+ JAR (built from dev branch)
  new/MSGFPlus.jar         # bigbio fork JAR (built from feature branch)
  data/                    # Downloaded benchmark data (mzML, FASTA, mods)
  results/                 # Benchmark output files and metrics
  download_data.sh         # Downloads PXD001819 benchmark data
  run_benchmark_multi.sh   # Multi-dataset benchmark (for remote server)
  run_local_benchmark.sh   # Quick local benchmark using repo test data
```

## Quick Start (Remote Server)

```bash
# 1. Download benchmark data
bash download_data.sh

# 2. Place JARs
cp /path/to/official/MSGFPlus.jar baseline/
cp /path/to/new/MSGFPlus.jar new/

# 3. Run all datasets
bash run_benchmark_multi.sh

# 4. Run single dataset
bash run_benchmark_multi.sh yeast
```

## Quick Start (Local)

```bash
# Uses test.mgf from src/test/resources (slower, ~6 min on laptop)
bash run_local_benchmark.sh
```

## PXD001819 CI benchmark (public data)

End-to-end search on public PRIDE mzML + `quantms-test-datasets` FASTA, with metrics compared to a checked-in baseline TSV. Config lives under `src/test/resources/benchmark/PXD001819/`; CI scaffold lives under `benchmark/ci/PXD001819/` so future datasets can add sibling folders.

The workflow runs on a fixed self-hosted runner profile (`self-hosted,linux,msgf-benchmark`) to keep CPU/memory comparable between runs.

See [`.claude/plans/pxd001819-ci-benchmark.md`](../.claude/plans/pxd001819-ci-benchmark.md) and [`benchmark/ci/README.md`](ci/README.md).
