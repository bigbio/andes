# Rust↔Java PIN diff report

- Java:  `/srv/data/msgf-bench/bench-merged-results/astral-java.pin` (136271 rows, 121654 scans)
- Rust:  `/srv/data/msgf-bench/bench-iter33-results/astral-rust.pin` (124471 rows, 121677 scans)
- Total unique scans: 121681

## Top-1-per-scan buckets

| Bucket | Count | % of total |
|---|---:|---:|
| both_target_same_peptide | 69,264 | 56.92% |
| both_target_diff_peptide | 6,548 | 5.38% |
| java_target_rust_decoy | 5,573 | 4.58% |
| rust_target_java_decoy | 5,201 | 4.27% |
| java_only_target | 3 | 0.00% |
| rust_only_target | 12 | 0.01% |
| both_decoy | 35,064 | 28.82% |
| java_only_decoy | 1 | 0.00% |
| rust_only_decoy | 15 | 0.01% |
| both_missing | 0 | 0.00% |
| **total** | **121,681** | 100.00% |

### Sample disagreements

**both_target_diff_peptide** (first 5):

- `(13, 'TPEEGEK', 'EAQTEGR')`
- `(14, 'MLQTPESR', 'FLYYDEK')`
- `(15, 'PEWTQR', 'RDLSER')`
- `(32, 'NEYNDR', 'MQEEEK')`
- `(72, 'PQYQKNTR', 'YYAHLGRR')`

**rust_target_java_decoy** (first 5):

- `(25, 'RYMEYR', 'WWHSYR')`
- `(33, 'QMAASPAK', 'AESWLR')`
- `(42, 'KLPAESR', 'LQRSPR')`
- `(50, 'AYPDSKGR', 'IYWADAR')`
- `(92, 'SQDVQR', 'QEMPPR')`

**java_target_rust_decoy** (first 5):

- `(36, 'SEQETR', 'CIQEEK')`
- `(47, 'VWFSQIEYIVLR', 'FQLLEKYEPLNR')`
- `(58, 'ELYETR', 'FMELQK')`
- `(65, 'ACPYYRSR', 'YYMIQGQR')`
- `(82, 'NLVQPFICTYCDK', 'DLRCEFSEPVYSR')`

**java_only_target** (first 3):

- `(13201, 'MCYGYGCGCGSFCR')`
- `(14402, 'MDPNCSCATGGSCSCASSCKCK')`
- `(30441, 'MDLSCSCATGGSCTCASSCK')`

**rust_only_target** (first 5):

- `(106652, 'DECGSR')`
- `(109707, 'AFLPLLK')`
- `(113285, 'MKRGQR')`
- `(114265, 'IITLQKK')`
- `(114482, 'CCIMGR')`

## Per-feature diff (agreement bucket: same scan + peptide, both target)

_73,568 PSMs in agreement bucket._

Sorted by mean |Δ| (Rust - Java), descending:

| Feature | n | mean Δ | median Δ | stdev | p5 | p95 | mean \|Δ\| | mean rel Δ | %frac \|relΔ\|>1% |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| MS2IonCurrent | 73,568 | +227.6 | +0 | 7554 | -4 | +4.1 | 229.1 | +0.0002286 | 0.2% |
| RawScore | 73,568 | -2.75 | -2 | 9.837 | -20 | +11 | 7.731 | +0.08281 | 95.6% |
| lnEValue | 73,568 | -2.966 | -3.349 | 3.431 | -7.86 | +3.359 | 3.851 | -6.138 | 99.2% |
| lnSpecEValue | 73,568 | +1.583 | +1.177 | 3.471 | -3.322 | +8.046 | 2.772 | +0.05723 | 95.4% |
| MeanRelErrorTop7 | 73,568 | +1.319 | +1.488 | 3.041 | -3.872 | +5.988 | 2.467 | -5.202 | 99.5% |
| MeanErrorTop7 | 73,568 | -1.784 | -1.521 | 2.046 | -5.449 | +1.039 | 2.079 | -0.2867 | 99.0% |
| StdevRelErrorTop7 | 73,568 | -1.173 | -0.47 | 2.893 | -7.047 | +2.678 | 1.963 | -0.08092 | 98.0% |
| NumMatchedMainIons | 73,568 | -1.397 | -1 | 2.247 | -5 | +2 | 1.942 | -0.1708 | 80.1% |
| StdevErrorTop7 | 73,568 | -0.4602 | -0.3216 | 1.837 | -4.111 | +2.341 | 1.281 | -0.04469 | 98.4% |
| DeNovoScore | 73,568 | +0.05596 | +0 | 3.373 | -3 | +4 | 1.04 | +0.002041 | 28.6% |
| longest_y | 73,568 | -0.09439 | +0 | 0.7469 | -1 | +0 | 0.179 | -0.006405 | 8.2% |
| matchedIonRatio | 73,568 | +0.02564 | +0.00556 | 0.2109 | -0.3 | +0.3968 | 0.1656 | -0.008698 | 97.1% |
| longest_y_pct | 73,568 | -0.06381 | -0.05147 | 0.06619 | -0.1667 | -0.004762 | 0.0684 | -0.1018 | 100.0% |
| longest_b | 73,568 | +0.03425 | +0 | 0.3606 | +0 | +0 | 0.05856 | +0.01994 | 3.7% |
| lnDeltaSpecEValue | 73,568 | -0.006047 | +0 | 0.1908 | +0 | +0 | 0.006047 | +nan | — |
| ExplainedIonCurrentRatio | 73,568 | -0.003079 | -6e-09 | 0.0159 | -0.01876 | +4.8e-07 | 0.003515 | -0.0317 | 18.6% |
| CTermIonCurrentRatio | 73,568 | -0.00319 | -4e-09 | 0.016 | -0.01903 | +3.4e-07 | 0.003479 | -0.04253 | 17.6% |
| NTermIonCurrentRatio | 73,568 | +0.0001113 | +0 | 0.001956 | -6.765e-08 | +0.0003862 | 0.0002653 | +0.02863 | 9.7% |
| enzC | 73,568 | +1.359e-05 | +0 | 0.009754 | +0 | +0 | 9.515e-05 | -4.198e-05 | 0.0% |
| enzN | 73,568 | +4.078e-05 | +0 | 0.006386 | +0 | +0 | 4.078e-05 | +0 | 0.0% |
| dm | 73,568 | +1.182e-06 | -1.061e-06 | 2.757e-05 | -3.855e-05 | +4.895e-05 | 2.19e-05 | -0.002893 | 56.7% |
| absdm | 73,568 | +8.839e-07 | +1.58e-07 | 2.757e-05 | -4.186e-05 | +4.681e-05 | 2.189e-05 | +0.01116 | 56.7% |
| isotope_error | 73,568 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| peplen | 73,568 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| charge2 | 73,568 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| charge3 | 73,568 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| charge4 | 73,568 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| enzInt | 73,568 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| IsolationWindowEfficiency | 73,568 | +0 | +0 | 0 | +0 | +0 | 0 | +nan | — |

## Notes

- Δ = (Rust value) - (Java value).
- `mean rel Δ` = mean of (Δ / |java|) over PSMs with |java| > 1e-12.
- `%frac |relΔ|>1%` = fraction of PSMs where the relative diff exceeds 1%.
- Agreement bucket restricts to scans + peptides present as target on BOTH sides; this strips ranking-flip + retention-only effects so the table measures FEATURE divergence specifically.
