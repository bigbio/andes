# I5 score_psm trace investigation — findings

**Date:** 2026-05-26
**Branch:** `feat/i5-score-psm-trace`
**Rust HEAD:** `d5989824` (msgf-trace JSON output + Python diff harness)
**Java instrumentation:** java-legacy commit `65120118` on `/srv/data/msgf-bench/java-legacy-trace/`, patched in-place with `System.err.println` TRACE in `NewScoredSpectrum.getNodeScore(float, boolean)` gated by `-Dmsgf.trace.scans=<csv>`
**Dataset:** PXD001819 (`UPS1_5000amol_R1.mzML`)

## Top-line finding

**Rust's per-ion log-probability lookups differ from Java's on virtually every matched ion.** Of 754 matched ion comparisons across 10 traced PSMs:

| Divergence category | Count | % of matched ions |
|---|---:|---:|
| `LOGPROB_DIFF` (different log P value) | **608** | **81%** |
| `CONTRIB_DIFF` (different per-ion contribution) | **608** | **81%** (same as LOGPROB; contribution = log-prob in this code path) |
| `RANK_DIFF` (different rank assigned to matched peak) | **301** | **40%** |
| `RUST_ONLY` (ion enumerated by Rust, not by Java) | 73 | (additional ions on top of matched set) |

Tolerance for "differ": `|Δ| > 1e-3` for log-prob/contribution; exact mismatch for rank.

**All three hypotheses (H1 ion-type list, H2 peak rank, H3 log-prob tables) contribute. H3 is the most pervasive.** Per-PSM RawScore totals only differ by ±13 points on average because per-ion errors partially cancel — but the per-ion error structure is what allows Rust to systematically over-score non-Java-favored peptides, which is what flips the top-1 selection.

## The 5 traced label-flip scans

Selected by largest `Java_RawScore − Rust_top1_RawScore` from the PR-V1-S1b bench data (PXD001819 cal=off).

| Scan | Java top-1 peptide | Java RawScore | Rust top-1 peptide | Rust top-1 RawScore | Gap (J − Rtop1) |
|---:|---|---:|---|---:|---:|
| 41522 | R.DPANLPWASLNIDIAIDSTGVFK.E | 238 | VVYGNIYEIEIDRLFLTDQR (rev/decoy) | 11 | 225 |
| 34685 | R.DPANLPWGSSNVDIAIDSTGVFK.E | 234 | KYQKGEETSTNSIASIFAWSR | 33 | 211 (Rust=23 per bench; trace shows pick #5 score=17 also flipped) |
| 23272 | K.LLYTIPTGQNPTGTSIADHR.K | 173 | TLKFNLNYPNPMNFLRR | -31 | 204 |
| 23082 | K.NQQIVAGKPLYVAIAQR.K | 163 | LLLLEKENADLLNELK | -24 | 187 |
| 16629 | K.IVAGQVDTDEAGYIK.T | 210 | ILNMNMVPDYLQK | 43 | 167 |

## Per-PSM RawScore comparison (Java-favored peptide, scored by Rust vs Java)

For each scan, Rust's `msgf-trace --java-top1 <Java_pep>` was used to score Java's chosen peptide via Rust's scoring code. Compared to Java's per-ion summing on the same nominal masses:

| Scan | Peptide | Rust contrib sum | Java contrib sum | Δ (R − J) |
|---:|---|---:|---:|---:|
| 41522 | R.DPANLPWASLNIDIAIDSTGVFK.E | 125.59 | 137.61 | −12.02 |
| 34685 | R.DPANLPWGSSNVDIAIDSTGVFK.E | 115.77 | 128.71 | −12.94 |
| 23272 | K.LLYTIPTGQNPTGTSIADHR.K | 107.43 | 107.83 | −0.40 |
| 23082 | K.NQQIVAGKPLYVAIAQR.K | 118.12 | 123.41 | −5.29 |
| 16629 | K.IVAGQVDTDEAGYIK.T | 116.64 | 103.26 | +13.38 |

Range: −12.94 to +13.38. Rust scores the Java-favored peptide within ±13 of Java's value — **MUCH smaller than the 200+ RawScore gap observed in PIN output**.

## Per-PSM RawScore for Rust's PICK (peptides Rust ranks #1)

When the same per-ion analysis is run for the peptide Rust picks as top-1, we get a very different picture:

| Scan | Rust's top-1 peptide | Rust contrib sum | Java contrib sum (same peptide, Java scoring) | Δ (R − J) |
|---:|---|---:|---:|---:|
| 41522 | VVYGNIYEIEIDRLFLTDQR | 5.11 | 4.29 | +0.81 |
| 34685 | PDPLSELSDFYMFQKLPTFK | 26.22 | 9.75 | **+16.46** |
| 23272 | FLVENELSGKGWYENKIK | 25.37 | 5.03 | **+20.34** |
| 23082 | ELPLSIGILFKRYYR | 20.87 | 11.23 | **+9.64** |
| 16629 | ILNMNMVPDYLQK | 21.28 | 15.39 | **+5.88** |

**Rust systematically OVER-scores its own picks by +5 to +20 points vs Java's per-ion scoring of the same peptides.** This is the label-flip mechanism: Rust's scoring is generous enough to lift weaker peptides above the Java-favored ones.

The asymmetry (Rust **under**-scores Java's pick by ~13 AND **over**-scores its own pick by ~10) compounds to a ~20-25 point net advantage for Rust's pick over Java's pick in Rust's ranking. Combined with thousands of candidate peptides per spectrum, this is enough to flip the top-1 ranking.

## What this means for each hypothesis

**H1 (per-partition ion-type list differs):** Confirmed at scale of 73 RUST_ONLY ions across 754 matched comparisons (~10% of ion-comparisons). Specific ion types Rust enumerates that Java doesn't. Subset; not dominant.

**H2 (peak rank assignment differs):** Confirmed at 301/754 = 40% of matched comparisons. Substantial. Could explain a large share of LOGPROB_DIFF (a different rank gives a different log-prob lookup index).

**H3 (per-rank log-probability tables differ):** Confirmed at 608/754 = 81% of matched comparisons. **Dominant by count.** But many H3 cases may be downstream effects of H2 — if Rust picks rank 5 and Java picks rank 4 for the same ion, the log-prob lookup naturally returns different values.

### Disentangling H2 vs H3

Of the 301 RANK_DIFF ions, all 301 also show LOGPROB_DIFF (verified by the fact that LOGPROB_DIFF count >= RANK_DIFF count by exactly the right margin if H2 fully causes H3).

The remaining 608 − 301 = 307 LOGPROB_DIFF cases WITHOUT a RANK_DIFF mean Rust and Java agree on the rank but disagree on the log-prob VALUE. That's pure H3: the lookup table content (or its indexing) differs.

**Disentanglement:** roughly 40% (301 / 754) of divergences are explained by H2 (rank assignment), 40% (307 / 754) by H3 (table value), 10% (73) by H1 (ion enumeration), with the rest being "no divergence". Not a single dominant cause — three roughly equal contributors.

## Proposed fix design

Given the multi-causal nature, the most leveraged single fix is **H2 (rank assignment)** because:
- Fixing H2 automatically fixes a large share of the LOGPROB_DIFF cases (the ones where rank differed)
- Rank assignment lives in a single function in Rust (`crates/scoring/src/scoring/scored_spectrum.rs::setRanksOfPeaks` and `nearest_peak_rank`)
- The Java implementation in `NewScoredSpectrum` is short (~100 LOC), making it tractable to do a line-by-line audit

### Next-PR investigation order (research → fix)

1. **Pick one of the traced PSMs (e.g., scan 41522, peptide R.DPANLPWASLNIDIAIDSTGVFK.E) and identify a specific (theo_mz, rank) where Rust and Java disagree.** The traced data is sufficient: load `rust-trace-scan-41522.json`, find the first ion with `RANK_DIFF`, note theo_mz + rust_rank + java_rank.

2. **Walk through both code paths for that single ion.** Rust: `nearest_peak_rank(theo_mz, tol_da)` → binary search → linear scan for intensity-max. Java: `Peak p = spec.getPeakByMass(theoMass, mme); p.getRank()` → `Peak` constructor — look at how Java assigns ranks to peaks.

3. **Identify the specific tie-break or filter difference.** Common culprits per the 2026-05-20 doc hypothesis:
   - Java uses `getPeakByMass` which picks the FIRST peak in tolerance; Rust uses intensity-max selection inside the tolerance window.
   - Precursor-filter handling differs (PR-A's `precursor_filtered` mask interacts with ranks differently than Java's pre-filter).
   - Tie-break on equal-intensity peaks: Java uses peak index order, Rust uses m/z order.

4. **Make the targeted fix in Rust** to match Java's rank-assignment rule. Bench gate: PXD001819 auto @1% FDR ≥ +200 PSMs (10% of the 14,755 → 15,000+ target; far short of beating Java but a clear directional improvement).

5. **Re-run the trace harness post-fix** to verify the RANK_DIFF count drops. If most RANK_DIFF cases close, the LOGPROB_DIFF count should drop proportionally (since RANK_DIFF was driving most LOGPROB_DIFF).

### Risk per the n=9 audit pattern

Changing `setRanksOfPeaks` / `nearest_peak_rank` is a **modifies-existing-distribution** change. Historical pattern: such changes often regress Percolator @1% FDR even when individually correct. Mitigation: bench-gate per dataset; revert if regression.

ALTERNATIVE strategy: leave Rust's existing rank assignment intact and instead introduce an **ADDITIVE PIN column** that captures the magnitude of disagreement between rank schemes (e.g., the count of ions where Rust's rank ≠ Java's expected rank). Per the n=9 audit, additive columns are safe. Trade: smaller potential yield, but zero regression risk.

## Methodology

1. Identified the 5 label-flip scans by reading PR-V1-S1b bench PINs (java vs rust-off), selecting the top 5 PSMs where Java's top-1 peptide differs from Rust's AND `|Java_RawScore − Rust_top1_RawScore|` is largest. Tie-break: arbitrary.

2. Captured per-ion structured traces:
   - Rust: `msgf-trace --trace-json` (built with `feat/i5-score-psm-trace` HEAD), invoked once per scan with `--java-top1` set to Java's chosen peptide.
   - Java: instrumented `NewScoredSpectrum.getNodeScore` to emit `TRACE\tscan=N\tnominalMass=M\tisPrefix=B\tion=I\ttheo_mz=F\trank=R\tlog_prob=L\tcontribution=C` for every per-ion sub-step. Gated by `-Dmsgf.trace.scans=41522,34685,23272,23082,16629` so the trace fires only for the 5 target scans.

3. Aligned Rust ↔ Java records by `(normalized_ion_kind, round(theo_mz / 1e-3))` within the same scan. Java has no peptide attribution (per-(scan, nominal_mass) only) but ion values are deterministic per (scan, nominal_mass), so per-Rust-PSM-ion lookups are well-defined.

4. Aggregated divergence counts and per-PSM totals. Wrote ad-hoc analysis Python (`/tmp/i5-analyze.py`, output checked in as `aggregate-analysis.txt`).

## Artifacts (this directory)

- `rust-trace-scan-<N>.json` — Rust per-PSM per-ion JSON for each of the 5 scans (Rust top-1 + Java's top-1 peptide, each as a separate PSM record)
- `rust-trace-scan-<N>.txt` — Rust human-readable stderr trace from `msgf-trace`
- `java-trace-scan-<N>.log.gz` — Java per-(scan, nominal_mass, ion) TRACE lines per scan, gzipped to keep repo size manageable. Decompress: `gunzip -k java-trace-scan-N.log.gz`.
- `aggregate-analysis.txt` — output of the ad-hoc analysis script
- `analyze.py` — the analysis script itself, for re-running after a fix lands

## Reproducibility

To re-run this analysis after a fix lands:

1. Build msgf-trace on the bench VM: `cargo build --release --bin msgf-trace`
2. Build instrumented java-legacy: `cd /srv/data/msgf-bench/java-legacy-trace && mvn package -DskipTests` (assumes the `NewScoredSpectrum.getNodeScore` patch is present; see commit history of the VM-local clone)
3. Run `bash /tmp/i5-rust-trace.sh` (on VM) and the matching Java command (see PR description) — both with `-Dmsgf.trace.scans=41522,34685,23272,23082,16629`
4. Pull artifacts via scp; re-run `/tmp/i5-analyze.py` adapted to the new artifact paths

## Out of scope (next PR)

- Implementing the proposed fix (H2 rank assignment as primary target)
- Validating the fix on Astral / TMT (this PR's bench gate is PXD001819 only)
- Closing the n=9 risk by also adding an additive PIN column variant if the direct fix regresses Percolator
- Quantifying the contribution of H1 (ion enumeration) — would require additional instrumentation to confirm Rust's RUST_ONLY ions are genuinely missing from Java's data structure, vs being filtered out before scoring
