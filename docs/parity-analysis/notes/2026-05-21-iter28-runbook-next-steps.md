# iter28 wrap-up — concrete runbook for closing the DeNovoScore -13 floor

_2026-05-21. The iter27/28 audit decomposed the gap into Layer 1 (score_psm) + Layer 2 (GF max headroom). Both need real instrumentation; this is the runbook so a future session can pick it up cold._

## State

- Astral 1% FDR after iter27: **31,298** (vs Java 35,818, gap 12.6%).
- iter27 pin-diff agreement bucket: DeNovoScore Δ median = -13, RawScore Δ median = -2; per-peptide-length spread is roughly constant -10 to -16, length-independent → constant offset NOT per-AA.
- All three obvious Layer 2 causes ruled out in iter28: cleavage credit constants, Acetyl-Prot-N-term variant inclusion, AA prior probabilities.

## Correlation in pin-diff agreement bucket

`pearson(DeNovoScore Δ, ·)` on 50,466 PSMs:

| feature   | r |
|---|---:|
| lnEValue  | +0.37 |
| lnSpecEValue | +0.37 |
| RawScore  | +0.34 |
| longest_y_pct | -0.08 |
| longest_b | +0.06 |
| ...all others | \|r\| < 0.05 |

DeNovoScore covaries strictly with the other GF-derived scores (RawScore, lnSpecEValue, lnEValue) and with NOTHING else. This is consistent with a single per-edge node/error scoring bug affecting both `score_psm` and `compute_inner` (which share the same per-edge scoring helpers).

## Step 1 — Layer 1 (score_psm -7 RawScore) instrumentation runbook

**Goal:** Find which of the three open hypotheses is right:
1. per-partition ion-type list differs
2. peak rank assignment differs  
3. per-rank log-probability tables differ

**Test case:** Pick scan **42510** peptide **CSACNVWR** from `docs/parity-analysis/diff/iter27/per_psm_diff.csv` (RawScore Δ = -7, DeNovoScore Δ = -13, fully tryptic). Same scan on both engines should produce same matched peptide.

**Java side:** Add a per-split log inside `FastScorer.getScore`. Suggested diff:

```java
// at the top of FastScorer.getScore
boolean trace = scorer.getScanNumArr() != null
    && scorer.getScanNumArr().length > 0
    && scorer.getScanNumArr()[0] == 42510
    && pepSeq != null
    && pepSeq.contains("CSACNVWR");

for (int i = startIdx; i < endIdx; i++) {
    int s = ...;  // per-split score
    if (trace) {
        System.err.printf("TRACE_JAVA\tpep=%s\tsplit=%d\tprefMass=%d\tsuffMass=%d\tcontribution=%d%n",
                          pepSeq, i, prefixGrid[i], suffixGrid[i], s);
    }
    total += s;
}
```

**Rust side:** Already exists (gated by `MSGF_TRACE_PEP` env var, see `crates/scoring/src/scoring/psm_score.rs:146-208`). Run with:

```bash
MSGF_TRACE_PEP=CSACNVWR ./target/release/msgf-rust ... > rust.tsv 2>&1
```

**Diff format:** Pair the per-split lines by `split` value; compute Δ per split. Outliers reveal the structural divergence (whether by-rank, by-ion-type, or by-peak-mass).

## Step 2 — Layer 2 (GF max-headroom -6) per-mass-bin runbook

**Goal:** For scan 32227 (YDCSFCGK, where Rust RawScore is HIGHER than Java by 4 but DeNovoScore is LOWER by 13), dump `gf.getMaxScore()` for each peptide_mass index in the tolerance window.

**Java side:** Add a one-line log to `DBScanner.computeSpecEValue`:

```java
for (int peptideMassIndex = minPeptideMassIndex; peptideMassIndex <= maxPeptideMassIndex; peptideMassIndex++) {
    ...
    PrimitiveGeneratingFunction gfi = new PrimitiveGeneratingFunction(graph);
    gfi.setUpScoreThreshold(minScore);
    gf.accept(gfi);
    if (specIndex == 32227) {  // tracing flag
        System.err.printf("GFMAX_JAVA\tscan=%d\tmass_idx=%d\tmax=%d%n",
                          32227, peptideMassIndex, gfi.getMaxScore());
    }
}
```

**Rust side:** Similar log in `crates/search/src/match_engine.rs` around line 584:

```rust
for nominal_mass_idx in min_peptide_mass_idx..=max_peptide_mass_idx {
    ...
    let gf = GeneratingFunction::with_score_threshold(&graph, min_score, aa_set)?;
    eprintln!("GFMAX_RUST\tscan={}\tmass_idx={}\tmax={}",
              spec_idx, nominal_mass_idx, gf.max_score());
    group.accept(gf);
}
```

**Hypothesis:** Either Rust enumerates the same bins with lower max (per-bin DP divergence — same root cause as Layer 1) or Rust enumerates FEWER bins (peptide_mass tolerance window logic divergence).

## Step 3 — num_distinct semantic fix (lnEValue divergence, separate)

Lower priority; documented as item #2 in `known-divergences.md`. Java uses `CompactSuffixArray.computeNumDistinctPeptides` (counts all DB substrings via SA-LCP). Rust uses `SearchIndex.ensure_distinct_peptide_counts` (counts enzyme-filtered candidates). Ratio explains `lnEValue Δ ≈ -4.56` exactly.

**Fix candidate:** Replace Rust's `enumerate_candidates` walk with an SA-LCP substring count, matching Java's `computeNumDistinctPeptides`. Significant work because Rust uses an `FxHashSet<u64>` fingerprint per length; the SA-LCP approach is structurally different.

Per the n=9 audit pattern this is "modifying-existing-distribution", so expect Percolator FDR within noise. Ship as low-priority Java-alignment commit.

## What NOT to try

- Trypsin efficiency 0.95 → 0.99999 alignment: only changes penalty (-3 → -11); doesn't affect max_score; expected no impact on Astral PSMs.
- Re-add edge-scoring to score_psm: rejected at iter17 (regressed -8K PSMs) and iter18 atomic-mirror (regressed -8K). Top-1-changing per n=9.
- Bit-exact Top7 error stats: rejected at iter23 (regressed -1,404 PSMs). Top-1-changing per n=9.
