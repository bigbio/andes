# Per-PSM scoring divergence: Rust scores Java target peptides 20+ points lower

_2026-05-20. msgf-trace on a representative label-flip scan (scan 21 from Astral LFQ Condition A REP1) localizes the root cause of the 25% top-1 label-flip rate to per-split node-score divergence in `score_psm`._

## Pattern in label-flip scans (iter16 PIN diff)

5 sample `java_target_rust_decoy` scans:

| scan | Java top-1 (target) | Java RawScore | Rust top-1 | Rust RawScore |
|---:|---|---:|---|---:|
| 21 | R.NEEQSR.D | 38 | K.TEAPC+57.02146GK.P (decoy) | 32 |
| 29 | -.M+15.995PPAGGPR.A | 13 | K.HARAGC+57.02146K.F (decoy) | 7 |
| 45 | K.RC+57.021SEDK.P | 17 | K.EKDC+57.02146DK.R (decoy) | 8 |
| 47 | K.VWFSQIEYIVLR.N | 19 | K.FQLLEKYEPLNR.Y (decoy) | 9 |
| 48 | K.VEEQEK.V | 23 | K.EIEGEGK.L (decoy) | 19 |

**Java's chosen target always scores higher than Rust's chosen decoy in Java's PIN, by 6-22 points.** This means Java IS finding the target with high confidence; Rust either doesn't enumerate it or scores it lower than the decoy.

## msgf-trace on scan 21 (representative case)

`msgf-trace --scan 21 --java-top1 "R.NEEQSR.D"` confirms Rust DOES enumerate `NEEQSR` (4 candidates from 2 protein entries). It appears at **#5 in Rust's top-10 queue**, NOT as top-1:

```
Rust top-10:
  #1: NVCGADK    score=27.00  iso_off=0  decoy=true  (Java says decoy too — different decoy)
  #5: NEEQSR     score=18.00  iso_off=1  decoy=false (= Java's top-1)
```

Per-split node-score breakdown for `R.NEEQSR.D` from Rust's `score_psm`:

```
spectrum_parent_mass=762.3336, peptide_mass=761.3304, peptide_nominal=743
split=1 aa[0]=N  pref_nom=114  suf_nom=629  score=0   (matched 2 sum +2.13, missing 5 sum -2.05)
split=2 aa[1]=E  pref_nom=243  suf_nom=500  score=-1  (matched 4 sum +0.77, missing 3 sum -1.41)
split=3 aa[2]=E  pref_nom=372  suf_nom=371  score=11  (matched 6 sum +13.17, missing 1 sum -1.68)
split=4 aa[3]=Q  pref_nom=500  suf_nom=243  score=1   (matched 2 sum +1.59, missing 1 sum -0.75)
split=5 aa[4]=S  pref_nom=587  suf_nom=156  score=3   (matched 3 sum +2.77, missing 0 sum 0)
breakdown_total = 14
PSM.score (from queue) = 18  (= breakdown 14 + cleavage credit +4)
```

**Java's RawScore = 38. Rust's breakdown_total = 14. The 24-point gap is in per-split node scoring.**

For Java's 38 (= 34 + cleavage credit 4), each split must average ~7 vs Rust's ~3. Splits 1, 2, and 4 are where Rust is weakest (scores 0, -1, +1).

## Hypothesis space

Three places `score_psm` could diverge from Java:

1. **Per-partition ion-type list**: Rust's trace shows 5 ions for seg=0 + 5 ions for seg=1 = 7 unique ions tried per split. If Java loads a different ion-mask from the same .param file, Java would match more (or score differently weighted) ions per split.

2. **Peak rank assignment**: For each matched ion, the score is read from a log-probability table indexed by peak rank. If Rust's ranks differ from Java's (because filtered peak counts differ, or because ranks are computed differently), the per-ion scores diverge.

3. **Per-rank log-probability tables**: Both engines load the same `frag_off_table` per partition from the .param. If Rust's table parsing produces different values than Java's, the score lookup itself diverges. Existing tests already cover this for BSA at high precision, so it's the least likely.

## What msgf-trace doesn't tell us alone

Without a side-by-side Java trace (Java currently has no per-split breakdown dump), we know Rust's number but not which specific ion-by-ion lookup diverges. To confirm hypothesis (1) vs (2) vs (3):

- **(1) ion-type list parity**: emit Rust's per-partition (charge, parent_mass_seg) ion list at startup; compare to a Java dump of the same .param. Need to instrument Java to dump.
- **(2) rank divergence**: dump Rust's ranks for the top-50 peaks of scan 21; compare to Java's ranks for the same peaks. Java's `NewScoredSpectrum` exposes ranks indirectly via `getRankOfPeak`. Need Java instrumentation.
- **(3) log-probability table parity**: parse the .param both engine ways; assert per-(partition, rank) score tables match byte-for-byte.

## Empirical scoring divergence quantification

The 5 label-flip samples above show Java RawScore - Rust RawScore (for Rust's chosen decoy) of 6 to 22 points. For Rust's score of Java's target peptide (not top-1 in Rust), the gap is larger — 20-24 points for scan 21.

Extrapolating: if Rust scores ALL agreement-bucket PSMs ~20 points lower than Java does (on average), Rust's score distribution is COMPRESSED relative to Java's. This compression makes targets and decoys harder to separate at the top, producing the 25% label-flip rate.

## Action items

1. **Compare per-partition ion lists between Java and Rust** on the HCD_QExactive_Tryp param. Instrument Java to print the ion list per partition; compare to Rust's existing diagnostic output (already prints per-partition ion-list sizes at startup).

2. **Spot-check 3-5 more label-flip scans** with msgf-trace — confirm the 20+ point gap pattern is universal across diverse PSMs, not specific to scan 21.

3. **Investigate hypothesis (2) — peak rank**: extend msgf-trace to dump Rust's ranks for the matched ions of a known PSM; instrument Java's `getRankOfPeak` to dump the same; diff.

This is a 2-3 day investigation requiring Java instrumentation. Until then, the structural scoring divergence remains the dominant unaddressed source of the 26% Astral gap.

## Reproducibility

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm.ebi.ac.uk '
  TRACE=/srv/data/msgf-bench/track-iter16-build/rust/target/release/msgf-trace
  PARAM=/srv/data/msgf-bench/track-iter16-build/src/main/resources/ionstat/HCD_QExactive_Tryp.param
  SPEC=/srv/data/msgf-bench/astral-data/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML
  DB=/srv/data/msgf-bench/astral-data/ProteoBenchFASTA_MixedSpecies_HYE.fasta
  $TRACE --spectrum $SPEC --database $DB --param $PARAM --scan 21 --java-top1 "R.NEEQSR.D"
'
```
