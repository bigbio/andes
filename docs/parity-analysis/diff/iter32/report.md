# Rust↔Java PIN diff report

- Java:  `/srv/data/msgf-bench/bench-merged-results/astral-java.pin` (136271 rows, 121654 scans)
- Rust:  `/srv/data/msgf-bench/bench-iter32-results/astral-rust.pin` (149616 rows, 121677 scans)
- Total unique scans: 121681

## Top-1-per-scan buckets

| Bucket | Count | % of total |
|---|---:|---:|
| both_target_same_peptide | 45,881 | 37.71% |
| both_target_diff_peptide | 18,864 | 15.50% |
| java_target_rust_decoy | 16,640 | 13.68% |
| rust_target_java_decoy | 13,637 | 11.21% |
| java_only_target | 3 | 0.00% |
| rust_only_target | 13 | 0.01% |
| both_decoy | 26,628 | 21.88% |
| java_only_decoy | 1 | 0.00% |
| rust_only_decoy | 14 | 0.01% |
| both_missing | 0 | 0.00% |
| **total** | **121,681** | 100.00% |

### Sample disagreements

**java_target_rust_decoy** (first 5):

- `(9, 'HAAENPGK', 'DEYIGPK')`
- `(21, 'NEEQSR', 'TEAPCGK')`
- `(29, 'MPPAGGPR', 'HARAGCK')`
- `(45, 'RCSEDK', 'EKDCDK')`
- `(47, 'VWFSQIEYIVLR', 'FQLLEKYEPLNR')`

**both_target_diff_peptide** (first 5):

- `(11, 'QKYFYR', 'QSELDRR')`
- `(13, 'TPEEGEK', 'MQHETK')`
- `(14, 'MLQTPESR', 'FLYYDEK')`
- `(15, 'PEWTQR', 'SSRSSHR')`
- `(32, 'NEYNDR', 'MQEEEK')`

**rust_target_java_decoy** (first 5):

- `(23, 'LHNYEDMLEKNK', 'FQARCCPLQNQK')`
- `(25, 'RYMEYR', 'IQQDSGCK')`
- `(33, 'QMAASPAK', 'SSENPLR')`
- `(37, 'WQDFTK', 'YCAERK')`
- `(42, 'KLPAESR', 'LQRSPR')`

**java_only_target** (first 3):

- `(13201, 'MCYGYGCGCGSFCR')`
- `(14402, 'MDPNCSCATGGSCSCASSCKCK')`
- `(30441, 'MDLSCSCATGGSCTCASSCK')`

**rust_only_target** (first 5):

- `(106652, 'DCPSCK')`
- `(113285, 'MKRGQR')`
- `(114265, 'LIVLKEK')`
- `(114482, 'YNEGCR')`
- `(117256, 'MTPHVMK')`

## Per-feature diff (agreement bucket: same scan + peptide, both target)

_50,480 PSMs in agreement bucket._

Sorted by mean |Δ| (Rust - Java), descending:

| Feature | n | mean Δ | median Δ | stdev | p5 | p95 | mean \|Δ\| | mean rel Δ | %frac \|relΔ\|>1% |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| MS2IonCurrent | 50,480 | +249.7 | +0 | 7895 | -4 | +4.06 | 251.2 | +0.0002466 | 0.2% |
| RawScore | 50,480 | -2.895 | -2 | 10.48 | -21 | +12 | 8.563 | +0.1294 | 96.4% |
| lnEValue | 50,480 | -2.907 | -3.396 | 3.832 | -8.346 | +4.141 | 4.058 | -7.002 | 98.9% |
| lnSpecEValue | 50,480 | +1.662 | +1.151 | 3.875 | -3.788 | +8.845 | 3.101 | +0.04392 | 95.5% |
| MeanRelErrorTop7 | 50,480 | +1.341 | +1.463 | 2.746 | -3.238 | +5.477 | 2.268 | -6.268 | 99.5% |
| NumMatchedMainIons | 50,480 | -1.295 | -1 | 2.36 | -5 | +2 | 1.961 | -0.147 | 80.2% |
| MeanErrorTop7 | 50,480 | -1.694 | -1.466 | 1.889 | -5.092 | +0.7995 | 1.941 | -0.2928 | 99.1% |
| StdevRelErrorTop7 | 50,480 | -1.008 | -0.4182 | 2.675 | -6.544 | +2.41 | 1.721 | -0.08505 | 98.1% |
| StdevErrorTop7 | 50,480 | -0.4512 | -0.336 | 1.703 | -3.885 | +2.149 | 1.154 | -0.07209 | 98.5% |
| DeNovoScore | 50,480 | +0.142 | +0 | 3.617 | -3 | +4 | 1.057 | +0.002925 | 27.7% |
| longest_y | 50,480 | -0.1274 | +0 | 0.8479 | -1 | +0 | 0.2199 | -0.009781 | 9.5% |
| matchedIonRatio | 50,480 | +0.04448 | +0.02381 | 0.2151 | -0.2847 | +0.425 | 0.1716 | +0.01371 | 96.9% |
| lnDeltaSpecEValue | 50,480 | -0.09642 | +0 | 0.8924 | +0 | +0 | 0.09642 | +nan | — |
| longest_y_pct | 50,480 | -0.06726 | -0.05494 | 0.06947 | -0.1667 | -0.004762 | 0.07208 | -0.1016 | 100.0% |
| longest_b | 50,480 | +0.03851 | +0 | 0.391 | +0 | +0 | 0.06426 | +0.02041 | 3.7% |
| ExplainedIonCurrentRatio | 50,480 | -0.004179 | -1e-08 | 0.01875 | -0.02749 | +4.7e-07 | 0.004668 | -0.03982 | 19.5% |
| CTermIonCurrentRatio | 50,480 | -0.004317 | -6e-09 | 0.01888 | -0.02808 | +3.7e-07 | 0.004641 | -0.05172 | 18.6% |
| NTermIonCurrentRatio | 50,480 | +0.0001377 | +0 | 0.002248 | -1e-07 | +0.0005455 | 0.0003208 | +0.02222 | 10.1% |
| enzC | 50,480 | +0 | +0 | 0.008902 | +0 | +0 | 7.924e-05 | -4.06e-05 | 0.0% |
| dm | 50,480 | +1.834e-06 | +9e-08 | 2.802e-05 | -3.929e-05 | +4.956e-05 | 2.229e-05 | +0.000823 | 60.7% |
| absdm | 50,480 | +1.536e-06 | +9.1e-07 | 2.802e-05 | -4.196e-05 | +4.836e-05 | 2.228e-05 | +0.009508 | 60.7% |
| enzN | 50,480 | +1.981e-05 | +0 | 0.004451 | +0 | +0 | 1.981e-05 | +0 | 0.0% |
| isotope_error | 50,480 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| peplen | 50,480 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| charge2 | 50,480 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| charge3 | 50,480 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| charge4 | 50,480 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| enzInt | 50,480 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| IsolationWindowEfficiency | 50,480 | +0 | +0 | 0 | +0 | +0 | 0 | +nan | — |

## Notes

- Δ = (Rust value) - (Java value).
- `mean rel Δ` = mean of (Δ / |java|) over PSMs with |java| > 1e-12.
- `%frac |relΔ|>1%` = fraction of PSMs where the relative diff exceeds 1%.
- Agreement bucket restricts to scans + peptides present as target on BOTH sides; this strips ranking-flip + retention-only effects so the table measures FEATURE divergence specifically.
