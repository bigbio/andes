# Deep audit of the remaining 12.4% Astral gap (post-iter24)

_2026-05-21. After iter24's acetyl-mod fix closed half the gap (26% → 12.4%), a deep audit of label-flip scans + per-feature distributions identified the dominant remaining divergence: **Rust's GF DP produces a score distribution ~8× wider than Java's** (max DeNovoScore 2,350 vs Java's 292). This is a structural scoring calibration issue, not a candidate-enumeration or feature-extraction bug._

## Audit summary

### jtrd label flips (Java target, Rust decoy)

16,437 jtrd flips remain after iter24. Modification analysis of Java's wins:
- 65.4% have NO mods (mods are NOT the dominant remaining issue)
- 24.2% have 1 mod (mostly Cam-C, Ox)
- 5.4% have Acetyl-Prot-N-term
- 90.7% have K./R. flanking (regular tryptic peptides)

The Rust picks (decoys that beat Java's targets) have ESSENTIALLY THE SAME mod distribution. So the flips are NOT a mod-handling bug — they're scoring divergences on standard peptides.

### Decoy/target set divergence

Top-1 PIN PSM sets:
- Java targets: 70,416 distinct bare sequences
- Rust targets: 75,368 distinct bare sequences (~5K more)
- Target overlap: 39,557 (37.2% of union)
- Java decoys: 43,513
- Rust decoys: 51,307 (~8K more)
- Decoy overlap: 17,887 (23.3% of union)

**Decoy ALGORITHM is identical** — both engines do full naive sequence reversal (confirmed by reading Java `ReverseDB.java:82-86` vs Rust `decoy.rs:13-18`). So the 76.7% non-overlap is a DOWNSTREAM EFFECT of scoring divergence, not an algorithm bug.

### Scan-level trace of jtrd flip (scan 472)

Java: K.HSWTPAR.F (RawScore=11), Rust: K.YYTSVPK.G (RawScore=20).

Rust enumerates HSWTPAR at top-N #3 with RawScore=11. Java's pick IS in Rust's candidate set, just ranked lower because YYTSVPK scored higher on Rust's score_psm.

This confirms the pattern from earlier scans (21 NEEQSR, 9 HAAENPGK): enumeration is fine, **scoring is the divergence**. Rust's score_psm + GF DP produce different rankings than Java's for the same peptide pool.

### TDC FDR (no Percolator) — biggest finding

| | TDC @ 1% FDR | Percolator @ 1% FDR | Percolator boost |
|---|---:|---:|---:|
| Java | 24,561 | 35,818 | +11,257 (45%) |
| Rust iter22b | 8,302 | 31,006 | +22,704 (273%) |
| Rust iter24 (+acetyl) | 8,506 | 31,390 | +22,884 (269%) |

**Rust's raw scoring is 3× worse than Java's** (8K vs 24K TDC). Percolator masks this via the 14-feature discrimination (Rust gets 270%+ boost vs Java's 45%). SpecEValue-only FDR is NOT a viable lever — it would tank Rust to 8.5K PSMs.

### DeNovoScore distribution audit

The cleanest evidence of the scoring calibration mismatch:

| Stat | Java | Rust iter24 |
|---|---:|---:|
| min | 0 | -13 |
| p5 | 35 | 21 |
| median | 68 | 71 |
| **p95** | **118** | **1,278** |
| **max** | **292** | **2,350** |

**42,399 Rust PSMs (28%) have DeNovoScore > Java's MAXIMUM**. This is Rust's GF DP computing per-mass `max_score` distributions that extend 8× wider than Java's.

DeNovoScore = `gf_max_score - 1` (the highest score achievable on ANY path through the AA graph at the candidate peptide's mass). The 28% target / 28% decoy split for high values shows it doesn't directly discriminate — it's NOISE in Percolator's feature space.

Capping DeNovoScore at 292 (Java's max) gave only +52 PSMs (within noise) — Percolator was already absorbing the noise via cross-validation. So the high DeNovoScore values don't directly hurt FDR, but they reflect a SYSTEMIC scoring calibration mismatch.

## Root-cause synthesis

Per-edge / per-node SCORES match Java within rounding (audit 2026-05-20). Decoy ALGORITHM matches Java exactly. Candidate ENUMERATION matches Java (with iter24's mod fix).

**The remaining divergence is in the GF DP's score-distribution accumulation**: Rust's AA graph evidently allows more high-scoring paths through the mass-bound search space than Java's, inflating `gf_max_score` 8×. This propagates to:
- Wider DeNovoScore distribution (cosmetic noise in Percolator)
- More confident lnSpecEValue at any given RawScore (Rust's spec_e_value median is 0.72 LOWER than Java's = MORE confident) — yet Rust's RawScore is 2 lower than Java's per identical PSM
- Per-PSM ranking divergence when comparing peptides at different mass bins (the label flips)

## Why per-feature fixes can't close this

Per the n=9 audit:
- iter17/18 added edge scoring to score_psm → -8K PSMs (changed top-1 selection)
- iter22b made intensity ratios bit-exact → flat
- iter23 made NumMatchedMainIons + error stats bit-exact → -1,404 PSMs

Fixing the underlying GF DP score-distribution shape would change the spec_e_value AND DeNovoScore AND RawScore feature distributions simultaneously. Percolator's calibration would shift; based on the n=9 pattern, this would likely regress unless Percolator's weights are also retrained.

## What COULD close the gap (not attempted)

1. **Audit Rust's GF DP score accumulation algorithm** in detail vs Java's `PrimitiveGeneratingFunction`. Find where Rust's max_score grows 8× faster. Fix the per-mass-bin distribution accumulation. ~1-2 weeks of work.

2. **Retrain Percolator on Java's PIN as starting weights** — feed Java's pin first, save weights, apply to Rust's pin. Not how Percolator natively works but possibly achievable via Percolator's `--init-weights` flag. ~1 day to validate.

3. **Cap DeNovoScore at Java's range AND adjust spec_e_value scaling** — try various calibration adjustments to mimic Java's narrower score distribution. Ad-hoc but might help marginally.

4. **Add new features that compensate for the calibration mismatch** — e.g., a `score_normalized = RawScore / DeNovoScore` ratio, or `score_z = (RawScore - dist_mean) / dist_stddev`. Additive per n=9; flat by default but could help.

## Bench numbers

| Iter | 1% FDR | Δ vs iter16 | Gap to Java |
|---|---:|---:|---:|
| iter16 baseline | 26,432 | — | 26% |
| iter20 tolerance fix | 30,983 | +4,551 | 13.5% |
| iter22b partition-ion ratios | 31,006 | +4,574 | 13.4% |
| iter24 + Acetyl mod | 31,390 | +4,958 | 12.4% |
| (Java reference) | 35,818 | | — |

**Total Astral gain: 18.8% over iter16 baseline (+4,958 PSMs).** Path to closing the remaining 12.4% requires GF DP score-distribution algorithm work, not feature-level audits.
