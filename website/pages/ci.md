# CI / benchmarks

Every merge to `main` builds andes and searches public PRIDE datasets end to end (search,
Percolator, protein groups). Larger sets, including match-between-runs quantification on a
design with known ratios, run once a week. Each run records the identifications at `q ≤ 0.01`, the measured error rate where
the database carries an entrapment, quantification accuracy, and the cost, so a change in any of
them is visible against the commit that caused it. Identifications are portable; wall time, CPU
time and memory are those of one run on the benchmark runner and vary between runs by 5–15%.

The harness is [`docs/benchmarks/ci/`](https://github.com/bigbio/andes/tree/main/docs/benchmarks/ci);
a dataset is added there with one entry. Results are kept on the
[`benchmarks` branch](https://github.com/bigbio/andes/blob/benchmarks/results.jsonl).

<div id="andes-ci" data-src="https://raw.githubusercontent.com/bigbio/andes/benchmarks/results.jsonl">
<p>Loading benchmark results…</p>
</div>
