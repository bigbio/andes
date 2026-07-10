# andes-glyco v2 campaign — code + experiments summary (2026-07)

Benchmark: PXD025455 Fc3_r1 (HCC_pool_Late_Fc3_r1, stepped-HCD N-glycopeptides), 523
the reference engine-confident backbones as truth; metric = backbone-correct @1% FDR (Percolator),
validated NUMERICally (CalcMass − GlycanMass vs truth backbone, tol 0.05 Da).

## Headline result

**andes: ~323/523 backbone-correct @1%, FDR-honest (2 decoys) — decisively beats
an open-source glyco engine (~222).** The gap to the reference engine's raw 523 is an EVIDENCE WALL, not a scoring
gap (see Conclusion).

## Code shipped (validated, on by default)

| change | file(s) | effect |
|---|---|---|
| `gp` fused selector `rank + K·ladder + J·core_y_hits + H·hyperscore` | glyco_psm.rs, glyco_search.rs, glyco_pin.rs | 221→301 @1% (the foundation) |
| count-rewarding hyperscore `ln(N_matched!)` | scoring/psm_score.rs | +7 (294→301) |
| glycan-list expansion (Fuc 3→4, Hex 3..12→2..14) | glycan_db.rs | surfaces more true backbones (part of the 323) |
| per-spectrum calibration features (Tailor/RankScoreFloat/strong_score/listwise) | glyco_search.rs | fixes dead-0 glyco features (correctness; @1%-neutral) |
| partial-glycan b/y evidence FEATURE (`PartialGlycanBY`) | backbone.rs, glyco_psm.rs, glyco_pin.rs, glyco_search.rs | +6 (317→323); sequence-specific evidence |

## Experiments (what moved @1%, what didn't)

- **gp selector ladder** 221 → gp 287 → gp2 (+core-Y) 294 → hyperscore 301. ✓
- **glycan expansion** (full 2510 list) 301 → 319 @1% — expansion PAYS once the peptide
  axis is strong (earlier "expansion crashes" was a weak-pipeline artifact). ✓
- **partial-glycan b/y feature** → 323. ✓ (offline: +30% true−decoy gap on fail-FDR set.)
- **REFUTED (@1%, held-out):** learned GBDT selector (+11 top-1 but −16 @1%: top-1 ≠ FDR
  yield), EdgeScore term, partial-glycan in the *selector* (P>0), retention (top-k 500,
  neutral), 2D glycan-decoy FDR (glycan-decoys are peptide-twins → break Percolator),
  intermediate-Y feature (poisoned Percolator), calibration features (neutral),
  glyco-regime base retrain (flat — a better model can't extract ions that aren't there).

## Diagnostics that pin the ceiling

- **FDR sweep:** andes maxes at **337 @100% FDR** vs 523 → NOT an FDR problem; the missing
  186 are never scored as a winning candidate.
- **The 186 missing are large / high-charge glycopeptides** (recovery 63%/54%/25% at
  z3/z4/z5; missing median backbone 1771 vs recovered 1606 Da). Biggest bucket = 119
  "outranked" (true backbone present, a wrong one outscores it on thin b/y).
- **Glycan-gap:** ≤7 recoverable with sulfate/phospho glycans; ~10 non-standard
  (likely the reference engine artifacts andes is correct to miss).

## The one law learned

Every change to winner **selection** (learned selector, EdgeScore, partial-glycan-P) HURT
@1%, because Percolator scores the *emitted winner's* features and the newly-selected
weak-evidence winners fail FDR while displacing gp's strong ones. Only **additive
features** helped, and modestly. **gp is FDR-optimal for selection — leave it.**

## Conclusion

~323 is andes' honest @1% ceiling on this stepped-HCD dataset, bounded by intrinsically
sparse backbone b/y in large/high-charge glycopeptides. Beating the reference engine's 523 needs a
different fragmentation modality (**ETD/EThcD** → backbone c/z ions the peptide keeps
under glycosylation) or an **orthogonal-truth benchmark** (synthetic/entrapment) to test
whether andes misses real IDs or the reference engine over-calls — not more scoring/selection/feature
work. Design for the next modality: `STRONGER-SPECTRAL-MODEL.md`.
