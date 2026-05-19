# Rust GF DP fails on ~5% of Astral PSMs where Java succeeds

_2026-05-19. Discovered via the iter12 PIN diff harness: 4.7% of all Rust PSM rows on Astral carry `DeNovoScore = -2147483648` (Rust's `i32::MIN` sentinel) and `lnSpecEValue = 0` (which is `ln(spec_e_value=1.0)`, the sentinel default at `match_engine.rs:302`)._

The sentinel value means the GF DP failed to compute for that spectrum and the per-PSM enrichment (`update_psm_enrichment` at `match_engine.rs:611-617`) never ran. Percolator is then fed broken lnSpec=0 values for ~5K-7K target PSMs per Astral run, which the FDR procedure can't discriminate.

## Empirical characterization (iter12 Astral PIN)

| Statistic | Count |
|---|---:|
| Total Rust rows | 149,351 |
| Sentinel rows (DeNovoScore < -1e9) | **6,982 (4.7%)** |
| Sentinel TARGET rows | 5,656 |
| Sentinel charge-2 rows | 6,819 (97.7% of sentinels) |
| Sentinel charge-3 rows | 155 |
| Sentinel charge-4 rows | 8 |

**Per-charge failure rate:**

| Charge | All rows | Sentinel rows | Failure rate |
|---|---:|---:|---:|
| 2 | 108,422 | 6,819 | **6.3%** |
| 3 | 38,186 | 155 | 0.4% |
| 4 | 2,464 | 8 | 0.3% |

The failure is **dramatically charge-2 specific** — 16× more likely for charge-2 spectra than higher charges. Sentinel population also skews toward short peptides (residues 6-8, peplen 8-10 in PIN convention).

## Smoking gun: same scan, same peptide, different outcome

Scan 100619, peptide `AFLASPEYVNLPINGNGK` (18 residues, charge 2):

| Field | Java | Rust |
|---|---:|---:|
| Label | 1 | 1 |
| RawScore | 119 | 97 |
| DeNovoScore | 119 | **-2147483648 (sentinel)** |
| lnSpecEValue | -48.1303 | **0 (sentinel)** |
| NumMatchedMainIons | 18 | 25 |
| Peptide | `K.AFLASPEYVNLPINGNGK.Q` | `K.AFLASPEYVNLPINGNGK.Q` |

Same scan, same peptide. Java's GF DP succeeded; Rust's failed. Java's DeNovoScore = 119 means `gf.getMaxScore() - 1 = 119`, so `gf.getMaxScore() = 120` — a finite computed distribution. Rust returned `i32::MIN`, meaning `group.is_computed()` returned false (no bin in the precursor-mass window produced a valid GF distribution).

## Why GF fails in Rust

`compute_spec_e_values_for_spectrum` (match_engine.rs:548-577) iterates the precursor nominal-mass window (~4 bins for Astral with iso_error -1..+2 and 10 ppm tolerance), builds a `PrimitiveAaGraph` per bin, and calls `GeneratingFunction::with_score_threshold(graph, min_score, aa_set)`. The GF can fail in two modes (`generating_function.rs:22-29`):

- `EmptyScoreRange { min, max }` when `max_score <= min_score` of the sink distribution (line 611)
- `SinkUnreachable` when `sink_dist` is `None` (line 606)

If ALL 4 bins fail, `group.is_computed()` returns false and the spectrum's PSMs keep their sentinel `de_novo_score = i32::MIN` and `spec_e_value = 1.0`.

The 6.3% charge-2 failure rate plus skew toward short peptides suggests the failure is in the AA graph for narrow precursor-mass windows where the path set through the graph is too constrained. Java handles this case (per `DBScanner.java:644` it asserts `specProb > 0` AFTER computing) — Java's GF apparently returns a valid (possibly trivial) distribution where Rust's returns `Err`.

## Next steps (require pride-linux-vm access)

1. **Identify the dominant failure mode**: add a thread-local counter for `EmptyScoreRange` vs `SinkUnreachable` returns in `compute_spec_e_values_for_spectrum`; emit aggregate counts in the yield-accounting summary. Run on Astral to see which mode dominates.

2. **Run `msgf-trace` on a known-failing scan**: scan 100619 with `--java_top1 K.AFLASPEYVNLPINGNGK.Q` will surface the per-bin node-score breakdown and pinpoint the first bin that fails plus the failure mode.

3. **Compare Java's GF outcome** for the same (scan, peptide) on the same fixture to confirm whether Java's GF produced a trivial-distribution success where Rust returned `Err`, OR whether Java's bin set is wider / its graph construction differs.

## Hypothesis for the fix

If the dominant failure is `EmptyScoreRange` (max == min in sink dist), the fix is to treat this as a valid (trivial) distribution rather than `Err`. The downstream `spectral_probability(score)` would then return 1.0 for `score >= the single value` and 0.0 otherwise — sensible behaviour. Caller code at `match_engine.rs:581-602` already has the `unwrap_or(1.0)` fallback so it would degrade gracefully.

If the dominant failure is `SinkUnreachable`, the AA graph construction itself is too narrow for these spectra and needs investigation. This is a bigger fix.

## Expected impact if fixed

5,656 target Astral PSMs would get valid `lnSpecEValue` + `DeNovoScore` instead of broken sentinel values. Percolator could then use these for FDR discrimination instead of treating them as outliers. Estimated upper bound: closing ~5K of the 9,400 PSM gap to Java would put Rust at ~31,500 PSMs (within ~12% of Java's 35,818). Lower bound depends on how Percolator currently weights the broken rows — could be smaller if the sentinel rows are already being effectively ignored.

This is the single highest-leverage remaining lead per the diff harness data.
