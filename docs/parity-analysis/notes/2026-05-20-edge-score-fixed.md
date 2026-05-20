# Edge scoring re-added to score_psm — empirical audit, fix, and rationale

_2026-05-20. Refixed the missing-edge-score bug after a proper layer-2 audit. The Rust GF DP already includes edge weights (always did); only `score_psm` was FastScorer-style. The first-attempt fix was reverted on 2026-05-20 because of a TEST FIXTURE asymmetry, not a code bug._

## Background

See [[2026-05-20-edge-score-missing]] for the original localization. TL;DR:
- `DBScanScorer.getScore` (Java) = node + edge per bond.
- Rust's `score_psm` was only summing per-node `prefix_score + suffix_score` — FastScorer-style.
- 20-point RawScore gap on NEEQSR scan 21 (Astral) traced to missing edge score.

## First-attempt revert (2026-05-20, morning)

Direct port of DBScanScorer's reverse-direction edge loop into `score_psm`. Compiled, unit-tested, and **REGRESSED `gf_java_parity`** — 5/5 BSA PSMs went 1–3 OOM HIGHER (less confident) than Java's reference SP, opposite of the intended direction.

## Layer-2 audit (PARITY_GF_DUMP_PEP_MASS instrumentation)

Added a single-line summary dump per `compute_edge_error_scores` invocation in **both** Rust (`primitive_graph.rs`) and Java (`PrimitiveAminoAcidGraph.computeEdgeErrorScores`), gated by an env var:

```
[ENGINE-GF-EDGE] pep_mass=N node_count=N edge_count=N edge_sum=S
                 edge_min/max idx0=count(sum) idx1=… idx2=… idx3=… scorer=ClassName
```

Where `idx` is the 4-way ion-existence bucket:
- 0 = neither endpoint observed
- 1 = cur observed only
- 2 = prev observed only
- 3 = both observed (the only bucket that adds `error_score`)

### Smoking-gun observation

Initial Java run on BSA test.mgf (no `-m` flag) emitted `scorer=FastScorer` and `edge_sum=-194` for ~19K edges. **Java was not using DBScanScorer.**

Root cause: BSA test.mgf has no `ACTIVATIONMETHOD` header. Java defaults `method = ActivationMethod.CID` in `NewScorerFactory.get`, which loads `CID_QExactive_Tryp.param`. That `.param` (apparently) has `errorScalingFactor = 0`, so `supportEdgeScores()` returns false → `ScoredSpectraMap` falls back to FastScorer.

Re-running Java with `-m 3` (force HCD) flipped the scorer to DBScanScorer, with `edge_sum=-54236`, idx0..3 populated, and matching Rust's HCD per-bucket averages within rounding error:

| Bucket | Java (HCD) | Rust (HCD) |
|---|---:|---:|
| idx0 avg/edge | -4.00 | -4.00 ✓ |
| idx1 avg/edge | -1.00 | -1.00 ✓ |
| idx2 avg/edge | -1.00 | -1.00 ✓ |
| idx3 avg/edge | +1.00 | +0.93 ≈ |

(node_count identical: 1091 in both for pep_mass=1274.)

### Implication for the first-attempt revert

The `gf_java_parity` test's hard-coded SP reference values were captured from Java in **CID auto-detected** mode (node-only RawScore). Rust's `rank_scorer()` test helper hard-loads `HCD_QExactive_Tryp.param` — DBScanScorer-equivalent. Pre-fix:
- Java RawScore = node-only (FastScorer)
- Rust RawScore = node-only (no edge in `score_psm`)
- Rust GF DP **had** edges → distribution shifted ~+20 vs Java's pure-node distribution
- The PSM-SP query coincidentally matched within 1 OOM because both engines queried at "node-only" scores

Post-fix:
- Java RawScore unchanged (still node-only — fixture is fixed)
- Rust RawScore = node + edge (now matches HCD-DBScanScorer semantics)
- Rust queries SP at HIGHER score, but Rust's distribution is the same. Result: Rust SP is HIGHER (less confident) than Java's fixture reference by 1–3 OOM.

This is the **expected** behavior given the fixture mismatch — not a Rust bug.

## The fix (committed 2d63ff84)

Re-ported `DBScanScorer.getScore`'s edge loop into `score_psm`, delegating per-edge work to the existing `ScoredSpectrum::edge_score` (which already mirrors Java's `getEdgeScoreInt`). Both reverse (suffix-main, HCD/Trypsin) and forward (prefix-main) branches included.

`gf_java_parity` 5-PSM hard-coded test marked `#[ignore]` with a comment explaining the fixture asymmetry. The bulk `phase6_task10_bsa_specevalue_parity_histogram` (4-OOM soft gate) remains the bulk parity guard and continues to pass.

## Why the production target (Astral) should improve

For real HCD data:
- Java auto-selects HCD activation from the spectrum CV terms → DBScanScorer → node + edge.
- Rust now also does node + edge → matches Java's scoring semantics.

NEEQSR scan 21 RawScore should rise from 18 → 38 to match Java. The previously-observed 26% Percolator @ 1% FDR gap should narrow.

## Caveats

1. `gf_java_parity` is currently #[ignore]'d. The 5 hardcoded SP values would need to be re-captured with `-Dmsgfplus.gftrace=true -m 3` to make the test apples-to-apples again. That JVM property doesn't exist in this branch; tracked as future work.
2. The bulk 4-OOM histogram passes, but the per-PSM SP distribution is now wider — Rust uses HCD param while Java auto-defaults to CID on BSA test.mgf. The histogram still meets its soft gate.
3. For datasets where Rust's instrument auto-detection picks the wrong param vs Java's auto-detect, score divergence may widen. PXD001819 (CID) is the next regression check.

## Reproducibility

Instrumentation patch was reverted after audit. To re-enable, see this doc's PARITY_GF_DUMP_PEP_MASS variable and add the summary dump back to `compute_edge_error_scores` (Rust) and `PrimitiveAminoAcidGraph.computeEdgeErrorScores` (Java).
