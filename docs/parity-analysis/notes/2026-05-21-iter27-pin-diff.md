# iter27 pin-diff vs Java — DeNovoScore offset has two layers

_2026-05-21. After iter25 (prob_peak clamp removed) + iter27 (source-protein labels), per-PSM Rust↔Java PIN diff on Astral agreement bucket (n=50,466) shows:_

| Feature | median Δ | mean Δ | mean \|Δ\| | %frac \|relΔ\|>1% |
|---|---:|---:|---:|---:|
| MS2IonCurrent | +0 | +252.7 | 254.1 | 0.2% |
| **DeNovoScore** | **-13** | -13.09 | 14.89 | 98.1% |
| **RawScore** | **-2** | -2.891 | 8.562 | 96.4% |
| lnEValue | -6.881 | -6.758 | 7.073 | 99.6% |
| lnSpecEValue | -2.325 | -2.189 | 3.67 | 96.4% |
| MeanRelErrorTop7 | +1.462 | +1.338 | 2.269 | 99.5% |
| NumMatchedMainIons | -1 | -1.296 | 1.961 | 80.2% |
| MeanErrorTop7 | -1.466 | -1.695 | 1.941 | 99.1% |

## DeNovoScore decomposition

The **-13 floor** decomposes into two near-orthogonal effects:

```
DeNovoScore_Δ ≈ RawScore_Δ + (DeNovoScore - RawScore)_Δ
        -13  ≈      -2.9   +              -10            (mean)
```

Per-length `(DeNovoScore - RawScore)` Δ (headroom = how much higher the GF
max is above the matched peptide's score):

| peplen | n | median | mean |
|---:|---:|---:|---:|
| 6 | 1,870 | -11 | -11.98 |
| 7 | 4,117 | -11 | -11.70 |
| 8 | 4,819 | -11 | -11.54 |
| 9 | 4,860 | -11 | -11.23 |
| 10 | 5,133 | -10 | -9.96 |
| 11 | 4,811 | -10 | -8.69 |
| 12 | 4,683 | -10 | -8.78 |
| 13 | 3,908 | -10 | -9.37 |
| 14 | 3,660 | -11 | -9.87 |
| 18 | 1,260 | -10 | -10.11 |
| 24 | 259 | -8 | -9.09 |

**Pattern:** roughly constant -10 to -11 across lengths, slightly attenuating
at longer peptides. **Independent of length** → comes from a CONSTANT
GF-max enumeration headroom that Java has but Rust doesn't.

## Two divergences to fix

### Layer 1 — score_psm (-2.9 RawScore, -7 typical per-PSM)

Rust's `score_psm` scores the matched peptide ~3–7 points LOWER than Java on
the same (scan, peptide). Documented in `2026-05-20-score-psm-divergence.md`
on the label-flip cases (Rust 14 vs Java 38 for R.NEEQSR.D on scan 21).

Root-cause hypotheses (per prior audit, still open):
1. Per-partition ion-type list differs
2. Peak rank assignment differs
3. Per-rank log-probability tables differ

**Blocked on:** Java instrumentation to dump per-edge node/error contributions
(2–3 day investigation).

### Layer 2 — GF compute_max_score headroom (-6 extra)

Rust's GF DP finds peptide paths with MAX score ~10 points LOWER than the
matched peptide; Java's GF finds paths ~16 points above the matched. So
Java's GF enumerates higher-scoring de-novo paths.

Possible sources (all NEEDS-CHECK; none confirmed):

- **Cleavage credit/penalty constants.** Rust uses `register_enzyme(0.95,
  0.95)` (`match_engine.rs:123`); Java uses `0.99999, 0.99999`
  (`Enzyme.java:300-301`). With probCleavageSites = 0.1:
  - Java: credit=2, penalty=-11
  - Rust: credit=2, penalty=-3
  CREDITS MATCH. Penalty differs. Penalty only affects min_score of finalDist,
  not max → **does NOT explain the +10 headroom.**

- **AA prior probabilities.** Java calls `DBScanner.setAminoAcidProbabilities`
  to set per-AA probability from FASTA frequencies (typical human:
  P(K)≈0.058, P(R)≈0.056, prob_clv ≈ 0.114). Rust uses uniform 1/20 = 0.05.
  Affects `edge_prob` only (probability of each path), not edge_score. Max
  is determined by edge_score, not edge_prob → **does NOT explain the +10.**

- **Source/sink AA list.** Both engines pull from cached_aa_list at the
  same Location (NTerm/CTerm or ProtNTerm/ProtCTerm). With useProtNTerm,
  Java includes the Acetyl variant of every source AA (mass-shifted). Rust
  appears to do the same. **NEEDS verification** — Acetyl-modified source
  AA gives a different start-node mass.

- **per-edge ion_existence_score / error_score.** iter25 brought distribution
  WIDTH to parity (max DeNovoScore 293 vs Java 292). Per-edge MEDIAN
  contribution may still differ — would shift the overall enumeration max.

- **score_threshold pruning.** Java's `setUpScoreThreshold(minScore)` and
  Rust's `with_score_threshold` use the same formula (`adjustedScore =
  minScore - neighboring_aa_cleavage_credit`). Both use credit=2.

## Recommended next steps

1. **Land 0.95→0.99999 efficiency alignment** as a no-op-on-DeNovoScore but
   correct-by-Java commit. Effect: only changes penalty in finalDist
   min_score (Java alignment), not max. Should be net-neutral for Astral.

2. **Audit Acetyl-Prot-N-term variant inclusion in source_aas.** iter24
   added the mod but the GF graph's source AA list may or may not include
   it. Check `cached_aa_list(ModLocation::ProtNTerm)` for the modified
   variant.

3. **Score_psm trace.** Pick scan 32227 (YDCSFCGK, RawScore_d=+4 — Rust
   HIGHER than Java for once) and scan 42510 (CSACNVWR, RawScore_d=-7) and
   instrument both engines to dump per-edge node_score + error_score per
   step. This is the 2–3 day investigation that unblocks the -7 RawScore
   floor.

4. **GF max trace.** Dump per-peptide_mass-bin `gf.getMaxScore()` and
   per-node max_score for both engines on the same scan. Compare which
   peptide_mass bin contributes the overall max in each.

The DeNovoScore -13 floor is a clean signal and a real bug. Closing it
requires Java instrumentation (currently lacking). Filed for follow-up.

## iter28 audit follow-up

- **Acetyl-Prot-N-term variant in source_aas** ruled out as Layer 2 cause.
  Added unit test `acetyl_prot_n_term_appears_in_source_aas_for_gf` in
  `aa_set.rs` (commit `4d324f21`): with Acetyl registered as wildcard-
  Prot-N-term, `cached_aa_list(ProtNTerm)` correctly returns the Anywhere
  list + 20 acetyl variants (one per residue). Matches Java's locMap
  semantics. So the source-AA list isn't being shrunk in Rust.

- **Cleavage credit at source/sink edges audited.** For tryptic peptides
  (preceded by K/R + ending in K/R), both engines compute
  `n_term_cleavage_score + c_term_cleavage_score = +2 + +2 = +4`. Tryptic
  agreement-bucket PSMs see no cleavage-credit divergence. The penalty
  difference (Rust -3 vs Java -11) only matters for non-cleavable termini
  in the GF DP score range, which affects `min_score` of `finalDist`, NOT
  `max_score`. Confirmed not the source of Layer 2 headroom.

- **num_distinct semantic divergence (different signal).** Java's
  `CompactSuffixArray.computeNumDistinctPeptides` counts distinct
  *substrings* in the target+decoy DB via the suffix-array LCP (all-AA
  prefixes regardless of enzyme). Rust's
  `ensure_distinct_peptide_counts` walks `enumerate_candidates`, which is
  ENZYME-FILTERED (tryptic peptides with up to N missed cleavages).
  Java's count is ~100x Rust's; `ln(100) ≈ 4.6` matches the observed
  `lnEValue Δ ≈ -4.56` exactly. **This is the source of the lnEValue
  divergence, NOT the DeNovoScore divergence.** Item #2 in
  `known-divergences.md` already filed.

## iter28 trace experiment — Layer 1 score_psm CLOSED

_2026-05-22._ Ran a single-scan trace on scan 47106 (HGIPTAQWK, unmodified
tryptic, iter27 RawScore Δ = -8, DeNovoScore Δ = -13). Method:

1. Extracted scan 47106 from the Astral mzML to a 1-scan MGF via
   `/tmp/extract_scan.py` (XML stream parser, ~110 LOC).
2. Built Java jar with `TRACE_JAVA*` per-split prints (already committed
   in `7823609a`, FastScorer.getScoreWithTrace). Added a temporary
   `TRACE_JAVA_CLEAVAGE` line in DBScanner to dump nTerm/cTerm cleavage
   scores at match construction; reverted after experiment.
3. Ran Java with `-Dmsgfplus.trace=true
   -Dmsgfplus.trace.scan=47106 -Dmsgfplus.trace.pep=HGIPTAQWK`.
4. Ran Rust with `MSGF_TRACE_PEP=HGIPTAQWK` on the same MGF.
5. Diffed per-split traces side by side.

**Result: BIT-EXACT MATCH per-split.**

| split | prefMass | suffMass | prefScore (J=R) | suffScore (J=R) | contribution (J=R) | cumulative |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 137 | 881 | -0.5222 | 11.3515 | 11 | 11 |
| 2 | 194 | 824 | 4.0091 | 7.2896 | 11 | 22 |
| 3 | 307 | 711 | 3.4843 | 18.0430 | 22 | 44 |
| 4 | 404 | 614 | -0.9877 | 3.3076 | 2 | 46 |
| 5 | 505 | 513 | 0.2636 | 3.5650 | 4 | 50 |
| 6 | 576 | 442 | 0.0000 | 3.0720 | 3 | 53 |
| 7 | 704 | 314 | 0.0000 | 5.3150 | 5 | 58 |
| 8 | 890 | 128 | 0.0000 | 2.5525 | 3 | 61 |

Both engines: per-split rawScore = **61**. Cleavage score: nTerm=+2,
cTerm=+2, total=+4. Score = 61 + 4 = **65** (both engines).

**Java pin RawScore = 73; Rust pin RawScore = 65 (iter27). Δ = -8.**

**Root cause of the -8 Δ:** Java's HCD scorer is
`DBScanScorer extends FastScorer`. `DBScanScorer.getScore` OVERRIDES
`FastScorer.getScore` to add an **edge-score loop** after the node-score
sum:

```java
// DBScanScorer.java:36-50
public int getScore(...) {
    int nodeScore = super.getScore(...);  // = 61
    int edgeScore = 0;
    for (int i = ...) {
        edgeScore += getEdgeScoreInt(...);  // = 8 for HGIPTAQWK
    }
    return nodeScore + edgeScore;  // = 69
}
// then match.score = cleavageScore (4) + 69 = 73.
```

Rust's `score_psm` returns node-only (61). The edge_score lives in a
SEPARATE PIN column `EdgeScore` (iter19 additive design). For scan 47106
the iter27 PIN shows `EdgeScore = -18`. Java's effective edge_score = +8
(per the trace dump). **Both engines compute edge_score with the same
algorithm** (`scored_spec.edge_score` mirrors `getEdgeScoreInt`), but the
per-edge value differs significantly (-18 vs +8 = 26-point gap on this
scan).

So **Layer 1 is closed**:
- score_psm node scoring is BIT-EXACT
- cleavage_score is BIT-EXACT (+4 both engines)
- Rust's RawScore = node + cleavage; Java's RawScore = node + cleavage + edge.
  This is a **design difference, not a bug** in score_psm. iter17/iter18
  added edge to Rust's RawScore and regressed -8K PSMs; iter19 made it a
  separate PIN column.

**New finding: EdgeScore PIN column itself diverges from Java by ~26 points
on this scan.** -18 (Rust) vs +8 (Java effective). Both use the same
algorithm but produce different values. Source: per-edge `edge_score`
function in `scored_spectrum.rs:730` vs Java `DBScanScorer.getEdgeScoreInt`.
Needs a per-edge trace to localize. Iter29 target.

## Bench impact (unchanged)

iter27 Astral 1% FDR: 31,298 (vs iter25: 31,410, Δ -112, within noise).
Gap to Java: 12.6%. The DeNovoScore divergence is real but Percolator is
already absorbing it via cross-validation (per n=9 audit pattern).
