# Rust↔Java PIN diff report

- Java:  `/srv/data/msgf-bench/bench-merged-results/astral-java.pin` (136271 rows, 121654 scans)
- Rust:  `/srv/data/msgf-bench/bench-iter27-results/astral-rust-iter27.pin` (149685 rows, 121677 scans)
- Total unique scans: 121681

## Top-1-per-scan buckets

| Bucket | Count | % of total |
|---|---:|---:|
| both_target_same_peptide | 45,870 | 37.70% |
| both_target_diff_peptide | 18,892 | 15.53% |
| java_target_rust_decoy | 16,623 | 13.66% |
| rust_target_java_decoy | 13,632 | 11.20% |
| java_only_target | 3 | 0.00% |
| rust_only_target | 14 | 0.01% |
| both_decoy | 26,633 | 21.89% |
| java_only_decoy | 1 | 0.00% |
| rust_only_decoy | 13 | 0.01% |
| both_missing | 0 | 0.00% |
| **total** | **121,681** | 100.00% |

### Sample disagreements

**both_target_diff_peptide** (first 5):

- `(11, 'QKYFYR', 'QSELDRR')`
- `(13, 'TPEEGEK', 'MTYTTR')`
- `(14, 'MLQTPESR', 'MSNRNNNK')`
- `(15, 'PEWTQR', 'SSRSSHR')`
- `(29, 'MPPAGGPR', 'MAAQGAPR')`

**rust_target_java_decoy** (first 5):

- `(20, 'EQCDTK', 'QMNDEK')`
- `(23, 'LHNYEDMLEKNK', 'FQARCCPLQNQK')`
- `(25, 'RYMEYR', 'PYSYNHR')`
- `(33, 'QMAASPAK', 'SSENPLR')`
- `(37, 'WQDFTK', 'YCAERK')`

**java_target_rust_decoy** (first 5):

- `(21, 'NEEQSR', 'TEAPCGK')`
- `(45, 'RCSEDK', 'EKDCDK')`
- `(47, 'VWFSQIEYIVLR', 'FQLLEKYEPLNR')`
- `(48, 'VEEQEK', 'TAQEWK')`
- `(58, 'ELYETR', 'FMELQK')`

**java_only_target** (first 3):

- `(13201, 'MCYGYGCGCGSFCR')`
- `(14402, 'MDPNCSCATGGSCSCASSCKCK')`
- `(30441, 'MDLSCSCATGGSCTCASSCK')`

**rust_only_target** (first 5):

- `(106652, 'DCPSCK')`
- `(113243, 'LYYALK')`
- `(113285, 'MKRGQR')`
- `(114265, 'LIVLKEK')`
- `(114482, 'YNEGCR')`

## Per-feature diff (agreement bucket: same scan + peptide, both target)

_50,466 PSMs in agreement bucket._

Sorted by mean |Δ| (Rust - Java), descending:

| Feature | n | mean Δ | median Δ | stdev | p5 | p95 | mean \|Δ\| | mean rel Δ | %frac \|relΔ\|>1% |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| MS2IonCurrent | 50,466 | +252.7 | +0 | 7909 | -4 | +4.1 | 254.1 | +0.000258 | 0.2% |
| DeNovoScore | 50,466 | -13.09 | -13 | 12.28 | -32 | +4 | 14.89 | -0.1788 | 98.1% |
| RawScore | 50,466 | -2.891 | -2 | 10.48 | -21 | +12 | 8.562 | +0.1294 | 96.4% |
| lnEValue | 50,466 | -6.758 | -6.881 | 4.087 | -13.12 | +0.1029 | 7.073 | -11.17 | 99.6% |
| lnSpecEValue | 50,466 | -2.189 | -2.325 | 4.115 | -8.584 | +4.792 | 3.67 | -0.1291 | 96.4% |
| MeanRelErrorTop7 | 50,466 | +1.338 | +1.462 | 2.744 | -3.254 | +5.475 | 2.269 | -6.267 | 99.5% |
| NumMatchedMainIons | 50,466 | -1.296 | -1 | 2.36 | -5 | +2 | 1.961 | -0.147 | 80.2% |
| MeanErrorTop7 | 50,466 | -1.695 | -1.466 | 1.89 | -5.094 | +0.7979 | 1.941 | -0.2931 | 99.1% |
| StdevRelErrorTop7 | 50,466 | -1.009 | -0.418 | 2.676 | -6.548 | +2.398 | 1.722 | -0.08367 | 98.1% |
| StdevErrorTop7 | 50,466 | -0.451 | -0.3355 | 1.703 | -3.885 | +2.145 | 1.154 | -0.07068 | 98.5% |
| longest_y | 50,466 | -0.1274 | +0 | 0.848 | -1 | +0 | 0.22 | -0.009756 | 9.5% |
| matchedIonRatio | 50,466 | +0.04442 | +0.02381 | 0.2151 | -0.2848 | +0.425 | 0.1716 | +0.01366 | 96.9% |
| lnDeltaSpecEValue | 50,466 | -0.1313 | +0 | 1.255 | +0 | +0 | 0.1313 | +nan | — |
| longest_y_pct | 50,466 | -0.06723 | -0.05494 | 0.06946 | -0.1667 | -0.005 | 0.07206 | -0.1015 | 100.0% |
| longest_b | 50,466 | +0.03836 | +0 | 0.3911 | +0 | +0 | 0.06432 | +0.02031 | 3.7% |
| ExplainedIonCurrentRatio | 50,466 | -0.004182 | -1e-08 | 0.01876 | -0.02749 | +4.7e-07 | 0.004671 | -0.03985 | 19.5% |
| CTermIonCurrentRatio | 50,466 | -0.00432 | -5.85e-09 | 0.01888 | -0.02808 | +3.7e-07 | 0.004644 | -0.05176 | 18.6% |
| NTermIonCurrentRatio | 50,466 | +0.0001374 | +0 | 0.002249 | -1e-07 | +0.0005453 | 0.0003215 | +0.02227 | 10.2% |
| enzC | 50,466 | +0 | +0 | 0.008903 | +0 | +0 | 7.926e-05 | -4.061e-05 | 0.0% |
| dm | 50,466 | +1.85e-06 | +1.125e-07 | 2.801e-05 | -3.922e-05 | +4.957e-05 | 2.229e-05 | +0.0008993 | 60.7% |
| absdm | 50,466 | +1.519e-06 | +8.8e-07 | 2.801e-05 | -4.196e-05 | +4.836e-05 | 2.228e-05 | +0.0097 | 60.7% |
| enzN | 50,466 | +1.982e-05 | +0 | 0.004451 | +0 | +0 | 1.982e-05 | +0 | 0.0% |
| isotope_error | 50,466 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| peplen | 50,466 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| charge2 | 50,466 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| charge3 | 50,466 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| charge4 | 50,466 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| enzInt | 50,466 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| IsolationWindowEfficiency | 50,466 | +0 | +0 | 0 | +0 | +0 | 0 | +nan | — |

## Notes

- Δ = (Rust value) - (Java value).
- `mean rel Δ` = mean of (Δ / |java|) over PSMs with |java| > 1e-12.
- `%frac |relΔ|>1%` = fraction of PSMs where the relative diff exceeds 1%.
- Agreement bucket restricts to scans + peptides present as target on BOTH sides; this strips ranking-flip + retention-only effects so the table measures FEATURE divergence specifically.
