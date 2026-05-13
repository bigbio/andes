# `score_psm` under-scoring on PXD001819/Astral — forensic finding

**Date:** 2026-05-13
**Branch state at write time:** `rust-implement @ beb6912` (local; NOT pushed)
**Status:** identified, not yet fixed. Investigation/fix is the next iteration.

## TL;DR

After the candidate-enumeration fixes that closed the TMT recall gap, Percolator
@ 1% FDR comparison surfaced a second, deeper bug:

| Dataset | Java @ 1% FDR | Rust @ 1% FDR | Δ | Δ % |
|---|---:|---:|---:|---:|
| **PXD001819** | 14,989 | 11,623 | **−3,366** | **−22%** |
| **Astral** | 35,818 | 24,828 | **−10,990** | **−31%** |
| TMT | 10,194 | 10,548 | +354 | +3.5% |

On PXD001819 and Astral, Rust's `score_psm` produces RawScore values roughly
**1/3 of Java's** for **identical (peptide, scan, charge)** inputs. The downstream
cascade kills SpecEValue computation (early return → `i32::MIN` /
`lnSpecEValue=0` sentinels) and shreds Percolator discrimination. Bug predates
this session's commits — past memory shows Rust at 14,839 @ 1% FDR on
PXD001819 (matching Java) on 2026-05-10; current Rust is at 11,623. Most
likely culprit: commit `0af1a37` (FastScorer prefix/suffix score cache) or
`be50dab` (compute_psm_features post-top-N hoist), both 2026-05-12.

## Smoking-gun example

PXD001819, scan=28787, both Java and Rust picked the same peptide:

| Column | Java | Rust |
|---|---:|---:|
| ScanNr | 28787 | 28787 |
| Peptide | `K.IVNEEFDQLEEDTPVYK.L` | `K.IVNEEFDQLEEDTPVYK.L` |
| charge | 2 | 2 |
| isotope_error | 0 | 0 |
| **RawScore** | **297** | **108** |
| **DeNovoScore** | 306 | **−2,147,483,648** (`i32::MIN`) |
| **lnSpecEValue** | −44.16 | **0** (= log(1.0) sentinel) |
| **lnEValue** | −28.60 | 0 |
| NumMatchedMainIons | 28 | 28 |
| longest_b | 14 | 14 |
| longest_y | 15 | 14 |
| ExplainedIonCurrentRatio | 0.61 | 0.31 |

**Same peptide. Same scan. Same charge. RawScore differs by 189 points.**

`NumMatchedMainIons` *matches* (both 28) — Rust IS finding the same b/y ions
Java does — but the per-ion contribution to RawScore is roughly 1/3. The
score-summing path is the bug, not the ion-matching path.

## Statistical confirmation across PXD001819

On 19,726 same-peptide scans (mod-mass normalized to 2 decimal places) on the
post-fix Rust pin (`bench-fixmod2-results/pxd001819-rust.pin`) vs Java
(`bench-merged-results/pxd001819-java.pin`):

| Statistic | Value |
|---|---:|
| Mean Java RawScore − Rust RawScore | +33.6 |
| Median \|gap\| | 31.0 |
| Max \|gap\| | 189.0 |
| Same RawScore (\|gap\|≤1) | 2.4% |
| Within 5 points (\|gap\|≤5) | 8.9% |
| Mean lnSpecEValue gap (Java − Rust) | −6.48 (Rust assigns much better probabilities) |
| Same lnSpecEValue (\|gap\|≤0.1) | 0.7% |

Top 10 worst gaps all involve unmodified Trypsin-tryptic peptides at lengths
17-20 — peptides without any modification at all. So **the bug doesn't depend
on mod handling** (the candidate-gen fix was orthogonal to this).

```
scan   Java   Rust   gap   peptide
28787   297    108   189   K.IVNEEFDQLEEDTPVYK.L
28825   305    116   189   K.IVNEEFDQLEEDTPVYK.L
28699   298    116   182   K.IVNEEFDQLEEDTPVYK.L  (same peptide, different scan)
33606   318    136   182   R.LESYVASIEQTVTDPVLSSK.L
32395   329    153   176   R.AVGSLTFDENYNLLDTSGVAK.V
28729   292    118   174   K.IVNEEFDQLEEDTPVYK.L
19338   284    112   172   K.EAC+57.021DWYAHSLNYNTPGGK.L  (CAM-C, no var mods)
30774   297    128   169   K.NINSETTDEQFQELFAK.F
21785   297    129   168   R.AEQLYEGPADDANC+57.021IAIK.N
27086   305    137   168   K.APEGELGDSLQTAFDEGK.D
```

The score gap is **systematic** (Rust always lower, never higher) and
**roughly proportional** to Java's RawScore (~1/3 ratio across the board).

## Cascade: how a RawScore bug becomes a Percolator collapse

1. Rust's `score_psm` undercounts by ~3× for some PSMs.
2. The undercounted PSM enters `compute_spec_e_values_for_spectrum` with
   `min_score = round(raw_score)`.
3. The GF DP builds graphs for each mass bin in the precursor window and
   asks each graph: "can your max_score reach `min_score`?"
4. Because `min_score` is derived from the under-scoring `score_psm` (which
   should be FAR higher), but the GF DP itself ALSO uses the same
   per-split prefix/suffix scores (which ALSO undercount), the GF graphs'
   max_score values reach the requested min_score for the easy cases but
   fail to construct a valid GF for many cases.
5. When `GeneratingFunction::with_score_threshold` returns `Err` for every
   mass bin in the window, `group.is_computed()` returns false, and
   `compute_spec_e_values_for_spectrum` early-returns without calling
   `update_psm_enrichment`.
6. The PSM is written to the PIN with `de_novo_score = i32::MIN` and
   `spec_e_value = 1.0` (the default sentinels from match_engine.rs:303,305).
7. Percolator gets a feature vector with bogus `lnSpecEValue=0` and a wildly
   negative `DeNovoScore = i32::MIN` for many target PSMs but not decoys,
   destroying its target/decoy discrimination.

The visible symptom (sentinel leak in PSMs) is real and is the immediate
write-path bug, but it's a SECONDARY consequence of the RawScore
undercounting. Fixing the sentinel write WITHOUT fixing RawScore would
hide the bug, not solve it.

## What is NOT the bug (ruled out)

- **Candidate enumeration:** confirmed by same-peptide diagnostic. Rust IS
  enumerating IVNEEFDQLEEDTPVYK. The bug is in scoring after enumeration.
- **Precursor m/z handling:** the Thermo Trailer fix (commit `dd2e4c8`) was
  a no-op — Rust was already reading the corrected monoisotopic m/z from
  selectedIon.
- **Fixed-mod slot counting / wildcard N-term variant doubling:** the two
  candidate-gen fixes (`a54e7e9`, `1ad272b`) closed the TMT recall gap and
  affect doubly-modded peptides; the score_psm bug surfaces on UNMODIFIED
  peptides and so cannot be downstream of those fixes.
- **Charge / isotope_error:** the diagnostic strictly filters to same
  isotope_error agreement. Both engines pick iso_error=0 on the same
  affected PSMs.
- **Peptide-string normalization artifact:** the diagnostic uses
  mod-mass-rounded matching, so `+229.163` vs `+229.16293` resolve to the
  same peptide. Same residue sequence + same flanking is verified.
- **Java RawScore field semantics:** confirmed equivalent to Rust's
  `psm.score`. Java's `getScore()` returns `cleavageScore + rawScore` from
  `FastScorer.getScore` (which sums prefix+suffix per split, matching
  Rust's `score_psm` algorithm in principle).

## Working hypotheses, ranked

### Hypothesis A (highest probability): `FastScorer` cache mis-sized or mis-populated

Rust's `ScoredSpectrum::new` precomputes
`prefix_score_cache: Vec<f32>` and `suffix_score_cache: Vec<f32>` indexed by
integer nominal mass:

```rust
let cache_len = (nominal_from(parent_mass).max(0) as usize) + 1;
for nominal_mass in 1..cache_len {
    prefix_score_cache[nominal_mass] =
        directional_node_score_inner(.., node_nominal as f64, true, ..);
    suffix_score_cache[nominal_mass] =
        directional_node_score_inner(.., node_nominal as f64, false, ..);
}
```

`score_psm` calls `cached_split_score` which sums `prefix_cache[i] + suffix_cache[i]`. Bug surface:
- `cache_len` derived from `parent_mass`. If `parent_mass` for ScoredSpectrum
  is different from what `peptide_nominal` resolves against, some splits
  fall outside the cache and silently return `None`, falling through to
  `node_score` ... which would NOT undercount; so this is probably not it.
- More likely: the loop populates `prefix_score_cache[nominal_mass]` with the
  result of `directional_node_score_inner(node_nominal, true, ..)` — but the
  initial `let ranks = vec![u32::MAX; n]` (scored_spectrum.rs:140) might be
  too few ranks if the cache loop runs BEFORE the rank table is populated.
- Past commit `0af1a37` had a known similar bug at the FastScorer integration:
  "fix the issue of Percolator" — earlier in this session the user asked
  about this. We landed the prefix/suffix cache fix, but the GF cache and
  the ranks-table population order may still have an edge case.

### Hypothesis B (medium probability): `directional_node_score_inner` undercounts

The cached and on-demand paths both call `directional_node_score_inner`.
The function iterates ion types from `segment_partition_cache`, looks up
each ion's expected peak by m/z in `spec.peaks`, applies a rank-based
log-prob score. If the iteration is missing ion types (e.g., only b/y, not
b2+/y2+/a/y-H2O) that Java's scorer considers, the per-split score is
systematically lower.

Java's `FastScorer.precomputeLogScoreTables` (FastScorer.java around line
40-60) populates prefix/suffix scores by iterating `getNodeScore` over
prefix/suffix node positions. `getNodeScore` enumerates ALL ion types in
the partition. If Rust's `directional_node_score_inner` enumerates only a
subset of partition ion types, it undercounts proportionally.

### Hypothesis C (lower probability): `parent_mass` / segment selection off

If Rust's `parent_mass` resolves to a different segment than Java's, the
ion-frequency lookup tables come from a different partition, producing
different scores. Less likely because partition lookups use the same
`partition_for(charge, parent_mass, last_seg)` formula on both sides.

## Reproducibility recipe

For the next investigator:

```bash
# 1. Reproduce the divergence with these exact inputs:
#    Spectrum:  /srv/data/msgf-bench/data/UPS1_5000amol_R1.mzML  (scan 28787)
#    FASTA:     /srv/data/msgf-bench/data/PXD001819_uniprot_yeast_ups.fasta
#    Peptide:   K.IVNEEFDQLEEDTPVYK.L
#    Charge:    2
#    Param:     HCD_QExactive_Tryp.param (Rust default)
#    Expected Java RawScore = 297
#    Actual Rust RawScore   = 108

# 2. Use msgf-trace (binary already exists) with --print-score-dist to get
#    Rust's per-split breakdown:
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm \
  /srv/data/msgf-bench/fixmod2-build/rust/target/release/msgf-trace \
  --spectrum /srv/data/msgf-bench/data/UPS1_5000amol_R1.mzML \
  --database /srv/data/msgf-bench/data/PXD001819_uniprot_yeast_ups.fasta \
  --param /srv/data/msgf-bench/fixmod2-build/src/main/resources/ionstat/HCD_QExactive_Tryp.param \
  --scan 28787 \
  --java-top1 K.IVNEEFDQLEEDTPVYK.L \
  --peptide IVNEEFDQLEEDTPVYK \
  --print-score-dist

# 3. Run Java MS-GF+ with the trace instrumentation enabled
#    (TRACE/-Dmsgfplus.trace=true in FastScorer.java) on the same scan +
#    peptide; compare per-split scores side-by-side.

# 4. The first per-split where Rust < Java is the bug entry point.
```

## What's safe to ship

Nothing from this session is shipped externally. Local branch state:
- `rust-implement @ beb6912` carries all this session's work (5 merges)
- Earlier release-quality work (the TMT / candidate-gen fixes) is correct
  and provides the per-dataset wins observed for TMT
- PXD001819 and Astral show the score_psm bug at FDR-controlled level
- Do NOT release `rust-implement` until the score_psm bug is rooted and
  fixed; the local-pin Percolator gap is the gate.

## Related commits (this session, in order)

```
beb6912 Merge fix/thermo-monoisotopic-precursor: candidate-gen + mzML trailer
1ad272b fix(search): drop Anywhere variants when fixed terminal mod is mandatory
a54e7e9 fix(search): fixed mods should not count against max_variable_mods_per_peptide
dd2e4c8 fix(input): read Thermo Trailer Extra Monoisotopic M/Z for precursor
076e1d4 Merge feat/param-fallback: Java NewScorerFactory ladder + coderabbit cleanup
e9edcb8 chore: cleanup from coderabbit review
1cfa402 feat(cli): Java NewScorerFactory fallback for bundled .param resolution
3bd9fc9 Merge feat/psm-candidate-handle: candidate_idx replaces clone in PsmMatch
dfbf4f9 perf(search): replace PsmMatch.candidate clone with candidate_idx handle
9cf9549 Merge feat/match-engine-hot-path: hoist edge_score inner-loop constants
120ae36 perf(scoring): hoist spectrum-constants out of edge_score inner loop
```

## Suspected root commits (pre-session)

```
0af1a37 (2026-05-12) perf(scoring): Track A — FastScorer prefix/suffix score cache per spectrum
be50dab (2026-05-12) perf(search): hoist compute_psm_features to post-top-N finalization
```

Past memory recorded Rust @ 1% FDR ≈ 14,839 on PXD001819 on 2026-05-10 —
matching Java. Current state is 11,623. The regression was introduced
between then and now.
