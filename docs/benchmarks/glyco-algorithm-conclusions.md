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
