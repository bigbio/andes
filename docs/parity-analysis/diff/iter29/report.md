# Rust↔Java PIN diff report

- Java:  `/srv/data/msgf-bench/bench-merged-results/astral-java.pin` (136271 rows, 121654 scans)
- Rust:  `/srv/data/msgf-bench/bench-iter29-results/astral-rust-iter29.pin` (149577 rows, 121677 scans)
- Total unique scans: 121681

## Top-1-per-scan buckets

| Bucket | Count | % of total |
|---|---:|---:|
| both_target_same_peptide | 45,856 | 37.69% |
| both_target_diff_peptide | 18,893 | 15.53% |
| java_target_rust_decoy | 16,636 | 13.67% |
| rust_target_java_decoy | 13,633 | 11.20% |
| java_only_target | 3 | 0.00% |
| rust_only_target | 14 | 0.01% |
| both_decoy | 26,632 | 21.89% |
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
- `(32, 'NEYNDR', 'MQEEEK')`

**rust_target_java_decoy** (first 5):

- `(20, 'EQCDTK', 'QMNDEK')`
- `(23, 'LHNYEDMLEKNK', 'FQARCCPLQNQK')`
- `(25, 'RYMEYR', 'PYSYNHR')`
- `(33, 'QMAASPAK', 'SSENPLR')`
- `(37, 'WQDFTK', 'YCAERK')`

**java_target_rust_decoy** (first 5):

- `(21, 'NEEQSR', 'TEAPCGK')`
- `(29, 'MPPAGGPR', 'HARAGCK')`
- `(45, 'RCSEDK', 'EKDCDK')`
- `(47, 'VWFSQIEYIVLR', 'FQLLEKYEPLNR')`
- `(48, 'VEEQEK', 'TAQEWK')`

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

_50,450 PSMs in agreement bucket._

Sorted by mean |Δ| (Rust - Java), descending:

| Feature | n | mean Δ | median Δ | stdev | p5 | p95 | mean \|Δ\| | mean rel Δ | %frac \|relΔ\|>1% |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| MS2IonCurrent | 50,450 | +251.6 | +0 | 7906 | -4 | +4.082 | 253 | +0.0002487 | 0.2% |
| RawScore | 50,450 | -2.898 | -2 | 10.48 | -21 | +12 | 8.564 | +0.1287 | 96.4% |
| lnEValue | 50,450 | -2.857 | -3.348 | 3.852 | -8.314 | +4.252 | 4.041 | -6.94 | 98.9% |
| lnSpecEValue | 50,450 | +1.712 | +1.199 | 3.895 | -3.763 | +8.952 | 3.127 | +0.04634 | 95.5% |
| MeanRelErrorTop7 | 50,450 | +1.344 | +1.463 | 2.737 | -3.233 | +5.477 | 2.267 | -6.28 | 99.5% |
| NumMatchedMainIons | 50,450 | -1.295 | -1 | 2.36 | -5 | +2 | 1.961 | -0.1469 | 80.2% |
| MeanErrorTop7 | 50,450 | -1.695 | -1.466 | 1.887 | -5.089 | +0.7928 | 1.941 | -0.2929 | 99.1% |
| StdevRelErrorTop7 | 50,450 | -1.009 | -0.4194 | 2.67 | -6.535 | +2.392 | 1.719 | -0.08391 | 98.1% |
| DeNovoScore | 50,450 | +0.3781 | +0 | 4.477 | -3 | +6 | 1.417 | +0.007059 | 34.3% |
| StdevErrorTop7 | 50,450 | -0.4525 | -0.3365 | 1.7 | -3.882 | +2.141 | 1.152 | -0.07113 | 98.5% |
| longest_y | 50,450 | -0.1275 | +0 | 0.848 | -1 | +0 | 0.2199 | -0.009755 | 9.5% |
| matchedIonRatio | 50,450 | +0.04452 | +0.02381 | 0.2151 | -0.2848 | +0.425 | 0.1716 | +0.01379 | 96.9% |
| lnDeltaSpecEValue | 50,450 | -0.09411 | +0 | 0.8988 | +0 | +0 | 0.09411 | +nan | — |
| longest_y_pct | 50,450 | -0.06726 | -0.05494 | 0.06945 | -0.1667 | -0.004926 | 0.07208 | -0.1015 | 100.0% |
| longest_b | 50,450 | +0.03845 | +0 | 0.3913 | +0 | +0 | 0.06438 | +0.02039 | 3.7% |
| ExplainedIonCurrentRatio | 50,450 | -0.004183 | -1e-08 | 0.01876 | -0.0275 | +4.7e-07 | 0.004672 | -0.03985 | 19.5% |
| CTermIonCurrentRatio | 50,450 | -0.004321 | -6e-09 | 0.01888 | -0.02808 | +3.7e-07 | 0.004646 | -0.05177 | 18.6% |
| NTermIonCurrentRatio | 50,450 | +0.0001379 | +0 | 0.002249 | -1e-07 | +0.0005492 | 0.0003216 | +0.02266 | 10.2% |
| enzC | 50,450 | +0 | +0 | 0.008904 | +0 | +0 | 7.929e-05 | -4.062e-05 | 0.0% |
| dm | 50,450 | +1.836e-06 | +8e-08 | 2.801e-05 | -3.928e-05 | +4.959e-05 | 2.23e-05 | +0.0007948 | 60.7% |
| absdm | 50,450 | +1.541e-06 | +9.185e-07 | 2.802e-05 | -4.196e-05 | +4.836e-05 | 2.229e-05 | +0.00948 | 60.7% |
| enzN | 50,450 | +1.982e-05 | +0 | 0.004452 | +0 | +0 | 1.982e-05 | +0 | 0.0% |
| isotope_error | 50,450 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| peplen | 50,450 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| charge2 | 50,450 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| charge3 | 50,450 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| charge4 | 50,450 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| enzInt | 50,450 | +0 | +0 | 0 | +0 | +0 | 0 | +0 | 0.0% |
| IsolationWindowEfficiency | 50,450 | +0 | +0 | 0 | +0 | +0 | 0 | +nan | — |

## Notes

- Δ = (Rust value) - (Java value).
- `mean rel Δ` = mean of (Δ / |java|) over PSMs with |java| > 1e-12.
- `%frac |relΔ|>1%` = fraction of PSMs where the relative diff exceeds 1%.
- Agreement bucket restricts to scans + peptides present as target on BOTH sides; this strips ranking-flip + retention-only effects so the table measures FEATURE divergence specifically.
