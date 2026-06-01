# Chimeric Phase 1 bench — FDR inflation finding

**Date:** 2026-05-28
**Branch:** `feat/chimeric-dda-plus` (Tasks 1-5; full-window search + multi-PSM emission, no MS1/shared-fragment refinement yet)
**Binary:** chimeric-build (sandybridge/LTO), cal=auto, top-N bumped 1→5 under --chimeric.

## Result (PSMs @1% FDR, Percolator --only-psms)

| Dataset | OFF | ON | Δ | OFF→ON wall | OFF tgt:dcy | ON tgt:dcy |
|---|---:|---:|---:|---:|---:|---:|
| PXD001819 | 14,755 | 17,015 | +15.3% | 56s → 87s | 2.5:1 | 1.2:1 |
| Astral    | 36,715 | 71,347 | +94%   | 6:16 → ~20min | 1.9:1 | 1.2:1 |
| TMT       | 9,605  | 9,608  | flat   | 2:34 → 3:55 | 1.8:1 | 1.1:1 |

## Conclusion: Phase 1 alone inflates FDR — not a real gain

- **Astral +94% is impossible as real sensitivity:** Astral has narrow isolation
  windows, so genuine co-isolation is minimal. A 2× gain there is an artifact.
- **Decoy fraction collapses toward 1:1 in every on-mode:** emitting ~5 PSMs/scan
  across the full window floods the set with secondary matches that are mostly
  random (target≈decoy). PSM-level Percolator cannot reject them because the
  decoy model doesn't encode the "few real peptides per scan" constraint.
- **Large wall cost** from the wide-window candidate explosion (Astral 6min→20min).

This confirms the design's flagged multi-PSM-FDR risk and the MSFragger-DDA+
insight: the **MS1 targeted-XIC isotope filter (Phase 2)** and **greedy
shared-fragment rescoring (Phase 3)** are LOAD-BEARING, not optional — they are
what removes the spurious co-IDs that Phase 1 alone over-counts.

## Status of Phase 1 code

- Mechanically correct: multi-PSM emission, unique SpecIds, isolation-window
  gating all work; `--chimeric` OFF is bit-identical (local golden + off-mode
  bench). Safe as default-off infrastructure for Phases 2-3.
- NOT shippable as a user-facing sensitivity feature on its own (counts are
  FDR-inflated). Either land behind a default-off experimental flag with this
  caveat, or hold the PR until Phase 2 (MS1 XIC) controls the FDR.

## Next step

Phase 2: MS1 targeted-XIC + isotope KL-divergence as a filter / additive
Percolator feature — the minimum needed to make chimeric counts trustworthy.
