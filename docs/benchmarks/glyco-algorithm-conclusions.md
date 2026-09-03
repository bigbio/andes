# Glyco algorithm — measured conclusions

Running record of what has been MEASURED on the plasma benchmark (PXD030622) with
entrapment FDP. Negative results are kept, not buried. A result that contradicts an
earlier claim REPLACES it here.

Ground rules used throughout:

- Yield alone never justifies a change; it is always paired with entrapment FDP.
- FDP falling toward 0.00% means **over-conservative**, not correct.
- Percolator's `q_min = 1/T_top` makes IDs@1% a step function, so fractions are POOLED.

## 2026-08-27 — 5-seed replication refutes the gate/selector "wins"

Re-scored the same pooled plasma PINs under 5 Percolator seeds. Searches unchanged; only
the FDR step is seeded, which is where the step-function instability lives.

| arm | mean PSMs | sd | range | vs base | effect/SE | mean FDP |
|---|---|---|---|---|---|---|
| pooled (baseline) | 223.8 | 56.4 | 128-267 | - | - | 0.81% |
| `--glyco-gp-m 10` | 241.2 | 81.0 | 126-354 | +17.4 | +0.39 | 1.46% |
| `--glyco-min-matched-ions 2` | 215.0 | 57.6 | 118-272 | -8.8 | -0.24 | 0.79% |
| combo | 150.0 | 63.1 | 111-262 | -73.8 | -1.95 | 5.15% |

**This REPLACES the earlier single-replicate claims.** The `gp_m 10` "+51% (383 PSMs)"
and the `min-matched-ions 2` "+24% glycopeptides (380 PSMs)" results were single draws
from distributions whose within-arm spread is 130-230 PSMs. Both effects are inside noise
(|effect/SE| < 0.4).

The only arm with a real signal is `combo` (`min-core-y 1` + `min-matched-ions 4` +
`gp-m 1`), and it is negative on **both** axes: fewer PSMs *and* ~6x higher entrapment FDP
(8.42/8.22/9.13% on three of five seeds). Over-gating does not trade yield for
calibration here — it loses both. A same-day single-replicate sweep agrees in direction:
`mi4` 124 PSMs and `combo` 115 against a 365 baseline.

**Methodological consequence:** pooling fractions is necessary but **not sufficient**. At
glyco counts the FDP estimate is seed-unstable too (0.00% on most seeds, 4.04/4.24/8.42%
on others), so a single-replicate "0.00% FDP" is not evidence of calibration. Replicate
across seeds before believing any glyco delta.

## Standing conclusions (unchanged by the above)

- **Generation is not the bottleneck.** Three measured expansions (wider glycan box,
  two-axis Y retention, isobar resolver) all moved yield DOWN and FDP DOWN.
- **The NeuGc isobar was a real, species-wrong-list defect** and its removal is the one
  entrapment-validated win: 268 -> 365 glycoPSMs at FDP 1.87% -> 0.55%, i.e. more IDs at
  lower true error.
- The headroom is **separation**, not coverage: the fused selector's heaviest terms are
  per-backbone and cannot discriminate between peptides competing for one scan.

## 2026-08-27 (later) — NeuGc arm replicated, and the design's detection floor

Same protocol applied to the campaign's one "validated" win. Same binary, same entrapment
FASTA, three fractions pooled; the only variable is `--glyco-taxon mammal` (NeuGc kept)
against the shipped `auto` (which excludes NeuGc on this data).

| arm | mean PSMs | sd | range | mean FDP |
|---|---|---|---|---|
| excl — NeuGc dropped (shipped default) | 223.8 | 56.4 | 128-267 | 0.81% |
| incl — NeuGc kept | 200.6 | 75.4 | 113-264 | 1.63% |

excl − incl = **+23.2 PSMs (+11.6%), effect/SE +0.55**. The single-replicate figure was
268 → 365 (+36%); five seeds give 201 → 224.

### The detection floor

Pooled within-arm sd ≈ 66 PSMs, so at 5 seeds/arm (80% power, α=0.05) the smallest
detectable effect is **~117 PSMs, i.e. ~58% relative change**. Detecting the observed
23-PSM difference would need ~129 seeds/arm; a 50-PSM difference still needs ~27.

This **corrects the framing in the entry above**: `gp_m 10` and `min-matched-ions 2` are
*not demonstrable at this power* rather than *refuted*. The design cannot separate a real
+10-20% effect from zero. Only `combo` (−73.8) approached the floor.

**Consequence for how this campaign is run:** stop spending seeds on sub-50-PSM glyco
effects — the instrument cannot resolve them at any sane cost. Use measures without
Percolator's `q_min = 1/T_top` step function: external agreement (Byonic mass-agreement
moved 74% → 91% under NeuGc exclusion — orthogonal, and far more decisive than any yield
delta), entrapment FDP at matched yield, or fixed-score-threshold counts.

**NeuGc exclusion stays the default.** The +36% yield claim is retired, but the change is
justified independently of Percolator: humans lack functional CMAH, every mainstream
human glycan list ships zero NeuGc, Byonic agreement improves sharply, and FDP direction
favours it. Right call, wrong headline number.

## 2026-08-29 — the primary defect was junk emission, not the selector

A four-way diagnostic (per-scan join against MSFragger on identical plasma data; data
forensics; a line-verified pipeline audit; the MSFragger mechanism from its paper and
docs) converged on one fact: **90.5% of andes's emitted glyco rows sit on scans that
contain no glycopeptide at all.** andes emitted a best guess for every scan clearing the
oxonium gate, with no evidence floor of any kind. Those rows are target/decoy coin flips
(measured target fraction 0.558), and Percolator learns its scoring model from them:
over the full PIN every feature's target/decoy AUC collapses to ~0.5 — RawScore INVERTS
to 0.471 — while on real glyco scans the same features separate fine (RawScore 0.663).
This one mechanism explains the immovable ~43% decoy-winner rate and why sixteen
selector-reweighting ablation arms all measured null.

Two corrections to earlier numbers, found during the same forensics: agreement with
MSFragger is 130/605 scans (21.5%), not 81 — the old comparator counted every
modified peptide as a disagreement — and 33% of MSFragger's rank-1 "glyco scans" are
its own decoys. On scans MSFragger is confident about (hyperscore ≥ 20), andes already
agrees ~82%; the engine was never catastrophically wrong on well-determined spectra.

### The emission floor, and its dose-response (pooled 3 fractions, 5 seeds, entrapment)

| floor (RawScore) | pooled rows | mean glycoPSMs @1% | sd | range | mean FDP |
|---|---|---|---|---|---|
| off | 22,592 | 244.4 | 93.7 | 101–363 | 1.72% |
| **3** | **5,072** | **256.8** | **16.5** | **236–277** | 1.76% |
| 6 | 2,884 | 192.0 | 59.0 | 119–246 | 0.88% |
| 10 | 1,714 | 153.0 | ~103 | 0–249 | erratic |

**The response is non-monotonic: the gate must clean without starving.** No floor
drowns Percolator in junk (seed sd 93.7); a moderate floor removes the junk bulk while
leaving ~5,000 pooled rows to train on (sd collapses to 16.5 — the demonstrated effect;
the +12 PSM yield delta is inside noise); harsher floors starve the training set, the
instability returns, and at floor 10 one seed returned zero (the `q_min = 1/T_top` step
function reappears at small row counts). At the R1 level, every one of the 130
externally-agreed correct answers survives even floor 10 — correct winners live far
above the decoy score distribution.

### The shippable form is run-adaptive, not a constant

An absolute score floor tuned on one dataset does not transfer (the July `min-core-y 2`
"plasma fix" cost mouse 161 of 707 IDs the same way). `--glyco-min-raw-score-quantile`
derives the floor from the run's own decoy winners — the run's null — and applies it
identically to target and decoy scans. Validated on R1: q=0.90 → derived floor 6.96,
q=0.95 → 12.23, q=0.99 → 20.76, holding 130/130, 129, and 116 agreements respectively.
The measured yield optimum (absolute floor ≈ 3) corresponds to **q = 0.775**; the
recommended setting is **0.75**, leaning toward retention because the starvation cliff
is steeper than the junk cost. Both flags default OFF pending the cross-dataset (mouse)
check; the starvation boundary depends on absolute pooled row count, so larger datasets
likely tolerate higher quantiles.

### Confirmation of the adaptive operating point (2026-08-29, pooled, 5 seeds)

`--glyco-min-raw-score-quantile 0.775` independently derived floors of 3.081 / 3.159 /
3.089 on the three fractions (2.5% spread — the per-run calibration is stable) and
reproduced the absolute floor-3 arm: 5,022 pooled rows, **243.0 mean glycoPSMs @1%
(sd 16.3, range 229–264)** against the absolute arm's 256.8 (sd 16.5) — within noise.
Seed stability, the demonstrated effect of the gate, is fully retained.

Caveat, recorded rather than glossed: the mean entrapment FDP point estimate is higher
in the adaptive arm (4.36% vs 1.76%), but at ~240 accepted PSMs one entrapment hit moves
FDP by ~2%, the per-seed values span 0–8.35% in both arms' lineages, and the difference
is ~1.5σ — statistically indistinguishable. FDP calibration at this yield scale remains
seed-unstable and unresolved; it is the main reason both flags stay OFF by default until
the mouse cross-check.

## 2026-08-29 (later) — the fit was starved and confounded: +50% from model-side fixes

With the emission gate in place, an autopsy of the remaining losses on an external
engine's confident scans found the 90 largest-bucket rows ALREADY carry the correct
peptide and dominate the threshold-setting decoys on every feature (RawScore 67 vs 7) —
the failure was Percolator's fit, not the evidence. Two mechanisms, both measured on the
same pooled PIN with no new searches:

1. **Training starvation.** Default `--trainFDR 0.01` leaves ~250 positives per
   iteration; `-F 0.05` trains on ~400 and the cross-validated fit stabilises.
   Single-feature TDC exposes the tail: RawScore alone yields 0 at 1% (a few decoys
   ride high in the tail) while the chance-calibrated ChanceMatchSurprise yields 201.
2. **Column confounds.** 19 columns were structurally constant or byte-identical
   duplicates (RawScoreCal ≡ RawScore under the collapse), and several per-scan
   spectrum-level columns (MS2IonCurrent, CandidateRankEntropy, ListwiseScoreGap,
   DeltaRankScore) let the small-sample SVM fit scan-quality confounds.

| arm (same PIN, 5 seeds) | mean @1% | sd | entrapment FDP |
|---|---|---|---|
| gate baseline, full PIN, `-F 0.01` | 256.8 | 16.5 | 1.76% |
| dead columns dropped + `-F 0.05` | 348.4 | 47 | 0.00 ×5 |
| **curated 52 columns + `-F 0.05`** | **384.6** | **23** | **0.00 ×5** |

**+128 PSMs (+50%) at ~10σ; 65% of the external engine (384.6/589), from 43% at the
start of the day.** Confident-scan recount: 131 → 210 accepted with the identical
peptide; wrong-peptide accepts stayed ≈0 (1); the ranking bucket (19 decoy winners) was
untouched, exactly as the mechanism predicts.

Shipped as: the glyco PIN always drops the 14 structurally-dead columns, and
`--glyco-pin-curated` (default off, HCD-oriented) emits the validated 52-column set.
Recommended pairing: `percolator --trainFDR 0.05` — rescoring remains pipeline-side;
this is a documented recommendation, not an andes feature. Guarded by
`crates/output/tests/glyco_pin_no_dead_columns.rs`.

Also closed today, measured post-gate: `--glyco-decorated-features` confirmed negative
(−29%, −2.3σ; under HCD backbone fragments LOSE the glycan, so bare-mass features are
the correct model), and `--glyco-decoy` through single-axis TDC breaks Percolator
outright (decoy twins near-tie their targets; needs a separate glycan-axis FDR consumer).
