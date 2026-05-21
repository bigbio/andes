# Iter21/22/22b — Java feature parity cleanups on top of iter20

_2026-05-21. Three feature-level fixes on top of iter20's tolerance win. iter21 (units) flat. iter22 broken (nominal-mass theo m/z drift). iter22b (accurate-mass theo m/z) is bit-exact for intensity ratios. Percolator FDR stays at iter20's level (~30,983–31,006). Confirms n=9 pattern: feature-level convergence doesn't translate to FDR gains once Percolator has the signal via correlated features._

## Bench summary

| Iter | Description | T/D | 1% FDR | 5% FDR |
|---|---|---:|---:|---:|
| iter16 baseline | C-4 + HIGH-2 + MS2IonCurrent | 92,977 / 56,570 | 26,432 | 31,783 |
| **iter20** | + feature-tolerance fix (20ppm) | 92,994 / 56,467 | **30,983** | 34,307 |
| iter21 | + units fix (Da → ppm Mean/StdevErrorTop7) | 92,828 / 56,507 | 30,888 (-95, noise) | 34,254 |
| iter22 | + partition-ion-list intensity sums (NOMINAL bug) | 92,877 / 56,444 | 30,500 (-388 vs iter21) | 34,333 |
| **iter22b** | partition-ion-list w/ ACCURATE mass | 92,900 / 56,488 | **31,006 (+23 vs iter20)** | 34,465 |

## Feature alignment trajectory (median Δ vs Java)

| Feature | iter16 | iter20 | iter22 | **iter22b** |
|---|---:|---:|---:|---:|
| NumMatchedMainIons | +3 | -1 | -1 | -1 (residual) |
| longest_b | +2 | 0 | -1 | **0** ✓ |
| longest_y | 0 | 0 | -3 | **0** ✓ |
| longest_y_pct | -0.038 | -0.106 | -0.356 | -0.056 |
| ExplainedIonCurrentRatio | -0.017 | -0.026 | -0.050 | **-1e-08** ✓ |
| CTermIonCurrentRatio | -0.015 | -0.018 | -0.036 | **-6e-09** ✓ |
| NTermIonCurrentRatio | -0.001 | -0.005 | -0.011 | **+0.00014** ✓ |

iter22b achieves bit-exact agreement with Java on the three intensity-ratio columns. The residual -1 on NumMatchedMainIons is from Java's per-bond max-intensity selection across multiple charge-1 ion types vs Rust's b/y-only count.

## iter22 root cause + fix

iter22 introduced partition-ion-list iteration for intensity sums but
called `IonType::mz(nominal)` which internally does
`real_mass = nominal / INTEGER_MASS_SCALER (0.999497)`. The recovered
"real_mass" drifts ~0.014 Da/residue from the true accurate residue mass
(NEEQSR's N: nominal 114 → 114.057 vs accurate 114.043), well outside
the 20 ppm window for high-resolution feature counting.

Net effect: iter22 found FEWER matches than iter20, REGRESSING intensity
ratios further (Explained -0.026 → -0.050).

iter22b fix: compute theo_mz directly from accurate residue mass:
```
theo_mz = prm_accurate / charge + offset
```
where `prm_accurate` accumulates `aa.mass + mod_.mass_delta` per residue
— matching Java's `peptide.get(i).getAccurateMass()` flow. Bypasses
IonType::mz entirely for the feature-counting path without affecting
the GF DP / score_psm paths (which use nominal-mass-indexed lookups
correctly via INTEGER_MASS_SCALER).

## iter21 result (units fix on iter20)

The units fix (Da → ppm for MeanErrorTop7/StdevErrorTop7) was -479 PSMs
solo on the iter16 baseline. On iter20's baseline: **flat (-95, noise)**.

The previous regression was caused by interaction with the wrong-tolerance
feature counting in iter16: Mean/StdevErrorTop7 computed from noise peaks
matched in the 0.5 Da window were polluting the Percolator features.
Once iter20 fixed the tolerance, the ppm units fix becomes neutral —
Percolator already extracts the signal via other features.

## n=9 audit pattern reinforced

iter21/iter22b confirm the refined n=9 pattern:
- **Top-1 selection preserved** → safe (T/D ratio stays around 1.65).
- **Feature alignment improves** → parity-pure but FDR-flat unless the
  feature was previously providing MISLEADING signal that Percolator
  was misled by. iter20 was the only feature fix where the previous
  feature value was actively misleading (noise-polluted) — that's why
  iter20 won big and iter21/22b are flat.

## Branch state

`iter19-additive-edge` HEAD has 7 commits stacked on `rust-implement`:
```
10d7874a fix(features): use accurate residue mass for partition-ion theo m/z (iter22b)
9d7cb846 fix(features): partition ion list intensity sums (iter22 — BROKEN, fixed by 10d7874a)
6b20edaa fix(search): MeanErrorTop7/StdevErrorTop7 units (iter21 — flat on iter20)
cf287c4a fix(features): hardcoded 20ppm/0.5Da feature tolerance (iter20 — THE WIN)
834b9ac0 docs(parity): iter20 results
c6820594 docs(parity): iter19 EdgeScore SAFE-FLAT
d8a8e66f feat(scoring): additive EdgeScore PIN column (iter19 — flat)
```

Ship recommendation: rebase to squash iter22 + iter22b into one (since
iter22 alone is broken). Then PR all 4 stable changes + docs.
