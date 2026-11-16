# CI benchmarks

Behind the [CI / benchmarks page](https://andes.quantms.org/ci/). `.github/workflows/benchmark.yml`
runs this on every push to `main` and appends one JSON line per dataset to `results.jsonl` on the
`benchmarks` branch; the page reads that file in the browser, so new results need no site rebuild.

| file | what |
|---|---|
| `datasets.json` | the datasets: spectra (a `fetch_spectra.sh` set), database, extra andes flags, and `entrapment` / `quant` evaluation |
| `run_ci.py` | runs andes + Percolator per dataset, records IDs, cost, true FDP and quant ratios |
| `test_run_ci.py` | known-answer check of the quantification metrics |
| `example.jsonl` | two real records (andes `0f630204`, TMT and UPS1, one laptop): the record format, and a file to preview the page with |

Local run (same as CI):

```bash
DATA=~/andes-bench
../reproduce/build_databases.sh "$DATA"
../reproduce/fetch_spectra.sh "$DATA" tmt
python3 run_ci.py --data "$DATA" --andes target/release/andes --out results.jsonl PXD007683
```

Preview the page against any results file: `https://andes.quantms.org/ci/?src=<url-of-results.jsonl>`,
for example the raw URL of `example.jsonl`.

**Adding a dataset:** add its files to `../reproduce/fetch_spectra.sh` (accession and filename
pattern), its database to `build_databases.sh` if new, then one entry in `datasets.json`. A
label-free design with known ratios adds `quant` with a `design.tsv` (run, condition) and
`amounts.tsv` (species, condition, amount) as in `../quant-eval/examples/`; the last condition
is the reference. Run it with the workflow's `datasets` input, then add it to the repo variable
`BENCH_DATASETS`.

**Runner:** set the repo variable `BENCH_RUNNER` (e.g. `["self-hosted","benchmark"]`) and
`BENCH_DATA` (a persistent data directory) to use a dedicated machine; runtime and memory are
only comparable on one machine. The default hosted runner (4 vCPU, 16 GB) runs the single-run
sets on each push (`BENCH_DATASETS`), and `proteobench-astral` (6 runs, 15.5 GB) once a week
(`BENCH_WEEKLY_DATASETS`), downloaded each time. On a self-hosted runner, add
`proteobench-astral` to `BENCH_DATASETS` to run it on every push.
