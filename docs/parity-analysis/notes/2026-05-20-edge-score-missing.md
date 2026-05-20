# Root cause of Rust↔Java score divergence: Rust's `score_psm` is missing edge scoring

_2026-05-20. Per-PSM trace harness localized the 20-point RawScore gap to **missing edge scores** in Rust's `score_psm`._

## Reproduction (scan 21 of Astral LFQ Condition A REP1, peptide `R.NEEQSR.D`, charge 2)

With **both engines using HCD_QExactive_Tryp.param** (Java with `-inst 3`, Rust auto-detects QExactive):

| Step | Java | Rust |
|---|---:|---:|
| Partition selection | (c=2, pm=761.42, seg=0) | (c=2, pm=761.42, seg=0) ✓ |
| Per-partition ion list (seg=0) | 5 ions: S_1_19, P_1_1, P_1_-27, S_1_1, S_1_2 | 5 ions: same ✓ |
| Per-node score (bond 1: pref(114)+suf(629)) | 0 | 0 ✓ |
| Per-node score (bond 2: pref(243)+suf(500)) | -1 | -1 ✓ |
| Per-node score (bond 3: pref(372)+suf(371)) | 11 | 11 ✓ |
| Per-node score (bond 4: pref(500)+suf(243)) | 1 | 1 ✓ |
| Per-node score (bond 5: pref(587)+suf(156)) | 3 | 3 ✓ |
| `getScore()` (node-score sum) | **14** | **14** ✓ |
| `cleavageScore` (N-term + C-term cleavage credits) | **4** | **4** ✓ |
| `edgeScore` (per-bond ion-existence + error) | **+20** | **0** ❌ |
| **Total RawScore** | **38** | **18** ❌ |

Per-node scoring is bit-exact. Cleavage credits match. The 20-point gap is **edge scoring** that Rust's `score_psm` does not call.

## Java's two-tier scorer architecture

Java picks the scorer at `ScoredSpectraMap.java:266-270`:

```java
if (scorer.supportEdgeScores()) {
    specKeyScorerMap.put(specKey, new DBScanScorer(scoredSpec, maxNominalPeptideMass));
} else {
    specKeyScorerMap.put(specKey, new FastScorer(scoredSpec, maxNominalPeptideMass));
}
```

- `FastScorer`: returns `sum_over_bonds(prefixScore[i] + suffixScore[i])` — node only.
- `DBScanScorer extends FastScorer`: `getScore = super.getScore() + edgeScore` where:
  ```java
  for each bond:
      edgeScore += getIonExistenceScore(partition, ionExistenceIndex, probPeak)
                 + (both ions found: getErrorScore(partition, curMass - prevMass - theoMass))
  ```

For HCD_QExactive_Tryp, `supportEdgeScores()` is true → Java uses **DBScanScorer**. Rust always uses FastScorer-style scoring.

## Rust's existing edge-scoring infrastructure (unused by `score_psm`)

Rust DOES parse the ion-existence table from .param and use it for the GF DP:

- `Param.ion_existence_table` (param_model.rs:37) — parsed from .param section 9 (line 385-399)
- `RankScorer::ion_existence_score(partition, idx, prob_peak)` (rank_scorer.rs:158) — lookup function
- `compute_edge_error_scores(...)` (primitive_graph.rs:633) — used during GF graph build

But `score_psm` (psm_score.rs:29-129) only computes:
```rust
for s in 1..n {
    let contribution = scored_spec.node_score(prefix_nominal, suffix_nominal, ...);
    total += contribution;
}
```

No call to `ion_existence_score`, no per-bond edge computation. The scorer used for the per-PSM RawScore on the production search path is structurally **FastScorer-style**, not DBScanScorer-style.

## Why the per-PSM trace harness missed this before

The earlier trace investigation (2026-05-20-score-psm-divergence.md) ran Java with `-inst 1` (HighRes), loading `HCD_HighRes_Tryp.param` (429 KB), while Rust auto-detects `QExactive` and loads `HCD_QExactive_Tryp.param` (741 KB). Different .param files, different partition counts (Java 92 vs Rust 140), different ion lists, different scores. That made everything look divergent.

After re-running Java with `-inst 3` (QExactive), the .param parsing aligns perfectly (both engines: 140 partitions, identical (c=2, seg=0) partition selection, identical 5-ion list). The remaining divergence is then localized cleanly to edge scoring.

## The fix

Implement DBScanScorer-style edge scoring inside `score_psm`. Per Java's `DBScanScorer.getScore`:

```java
@Override
public int getScore(...) {
    int nodeScore = super.getScore(...);  // = FastScorer.getScore = sum of (prefix+suffix) per bond
    int edgeScore = 0;
    if (!isNodeMassPRM) {  // reverse direction (suffix-main, typical HCD)
        int nominalPeptideMass = nominalPrefixMassArr[toIndex - 1];
        for (int i = toIndex - 2; i >= fromIndex; i--)
            edgeScore += getEdgeScoreInt(
                nominalPeptideMass - nominalPrefixMassArr[i],
                nominalPeptideMass - nominalPrefixMassArr[i + 1],
                (float)(prefixMassArr[i + 1] - prefixMassArr[i]));
    } else {  // forward direction (prefix-main)
        for (int i = fromIndex; i <= toIndex - 2; i++)
            edgeScore += getEdgeScoreInt(
                nominalPrefixMassArr[i],
                nominalPrefixMassArr[i - 1],
                (float)(prefixMassArr[i] - prefixMassArr[i - 1]));
    }
    return nodeScore + edgeScore;
}
```

`getEdgeScoreInt(curNominalMass, prevNominalMass, theoMass)`:
- `nodeMass[i] = scoredSpec.getNodeMass(NominalMass(i))` — the OBSERVED main-ion m/z if a peak matches, else -1
- `ionExistenceIndex` = `(curMass >= 0 ? 1 : 0) + (prevMass >= 0 ? 2 : 0)` → 0/1/2/3
- `edgeScore = ion_existence_score(partition, ionExistenceIndex, prob_peak)`
- If both ions exist: `+ error_score(partition, curMass - prevMass - theoMass)`

Implementing this in Rust requires:
1. A per-spectrum `node_mass[]` array (= main-ion observed m/z if peak exists, -1 otherwise) — analogous to `DBScanScorer.nodeMass`. Rust already has `ScoredSpectrum::observed_node_mass` which returns `Option<f64>`. Need a cached `Vec<Option<f64>>` indexed by nominal mass.
2. In `score_psm`, iterate bonds and accumulate edge score per Java's reverse-direction loop (since HCD/Trypsin → suffix-main → reverse).
3. Lookup `ion_existence_score` and `error_score` from the cached partition.

The fix is ~50-80 LOC. Mirrors Rust's existing `compute_edge_error_scores` in primitive_graph.rs.

## Expected impact

For NEEQSR scan 21: Java's edgeScore = +20, RawScore goes 14 → 38. The same pattern likely applies to most agreement-bucket PSMs.

Across the 49,538 agreement-bucket PSMs (diff harness), the per-PSM RawScore divergence is median ~+22 (Rust under Java). Implementing edge scoring should close most of that. If Percolator's discrimination improves with the higher (correctly-scaled) RawScores, the 26% gap to Java may close significantly.

This is conceptually an ADDITIVE fix at the SCORE level — Rust currently doesn't compute edge scoring at all. Unlike R-3 or C-5b (which modified an existing feature distribution), this ADDS a missing per-PSM scoring component. By the n=6 audit pattern, additive fixes don't trigger Percolator-recalibration regressions; they introduce new (correct) signal.

## Reproducibility

Trace re-run on scan 21 with correct `-inst 3`:
```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm.ebi.ac.uk '
  MSGF_TRACE_PARTITIONS=1 MSGF_TRACE_SCAN=21 \
    MSGF_TRACE_NODE_MASSES=114,156,243,371,372,500,587,629 \
    java -Xmx32g -jar /srv/data/msgf-bench/MSGFPlus-traced.jar \
      -s ... -d ... -inst 3 -e 1 -protocol 0 ... 2> trace.log
'
```

Rust trace at the same scan:
```bash
$TRACE --spectrum ... --database ... --param HCD_QExactive_Tryp.param \
       --scan 21 --java-top1 "R.NEEQSR.D"
```
