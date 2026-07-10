# Glyco FDR theory — peptide-axis vs glycan-axis decoys, 2D-FDR, and a Percolator-only realization

*Theory note for andes glyco (branch `glyco-phase1`). Clean-room: all algorithms below are
from published papers; no a commercial glyco engine/the reference engine code is used. See license notes at the end.*

## 1. Why a glycopeptide needs TWO FDR axes

A glycopeptide-spectrum match (GPSM) is a pair `(peptide P, glycan G)`. A wrong ID can be
wrong in the peptide, in the glycan, or in both. A single decoy pile (one `Label`) cannot
separate these error modes: the peptide backbone b/y evidence and the glycan Y-ladder /
oxonium evidence live on different, weakly-correlated feature subspaces. This is exactly the
failure andes measured — feeding ~352 K glycan-decoy rows into Percolator under one `Label`
collapsed recovery **29.4% → 4.4%**: the decoys differ from targets only in the Y-ladder, so
they flood the −1 pile, Percolator over-weights `YLadder`, and genuine targets are under-ranked
(SPA2/PHASE1 notes). The fix is the a glyco search engine **separate-axis (2-dimensional) FDR**.

## 2. a glyco search engine 2.0 two-dimensional FDR (exact)

a glyco search engine 2.0 controls FDR on the peptide moiety and the glycan moiety independently, then
combines by inclusion–exclusion ([Liu 2017, *Nat Commun* 8:438, PMC5585273](https://pmc.ncbi.nlm.nih.gov/articles/PMC5585273/)):

```
FDR^x = FDR^Px + FDR^Gx − FDR^{P∩G}x
```

- `FDR^Px` — peptide-axis FDR at the GPSM score cut `x`, from **reversed-sequence peptide
  decoys** (standard target-decoy competition on the backbone score `ScoreP`).
- `FDR^Gx` — glycan-axis FDR at cut `x`, from **glycan decoys** (recipe §3).
- `FDR^{P∩G}x` — the intersection: fraction of GPSMs that are **both** peptide-decoy **and**
  glycan-decoy at cut `x`. Subtracting it removes the double-counted joint-false term so
  `FDR^x` is a proper upper bound on the union of the two error modes.

Each axis is `FDR = N_decoy / N_target` (cumulative, sorted by descending GPSM score).
a glyco search engine scores each GPSM as a linear mixture `ScoreGP = w·ScoreG + (1−w)·ScoreP` (`w = 0.35`
in the paper); the 2D-FDR is applied on the combined ranking.

## 3. Glycan decoy construction (clean-room recipes)

**a glyco search engine Y-rung shift** (the canonical recipe): apply a random mass shift to **all
peptide+Y ions EXCEPT the two anchors** — `Y0` (bare peptide, no monosaccharide) and `Y1`
(peptide + 1 innermost HexNAc). Y0/Y1 are the peptide-mass anchors and are left intact so the
glycan decoy still competes on the *same backbone*; only the glycan structure evidence is
scrambled. A GPSM matching <2 trimannosyl-core ions is filtered ([PMC5585273](https://pmc.ncbi.nlm.nih.gov/articles/PMC5585273/)).
This is the direct analogue of what andes needs: shift Y-rungs above Y1, keep Y0/Y1.

**a cross-spectrum glyco engine / random-shift recipe**: generate a decoy spectrum per GPSM by shifting the
*m/z* of each MS2 fragment by a random **1–30 m/z**; peptide and glycan FDR each use
`FDR = 2·N_decoy / (N_decoy + N_target)` (concatenated TDC form)
([Fang 2022, *Nat Commun* 13:1900, PMC8990002](https://pmc.ncbi.nlm.nih.gov/articles/PMC8990002/)).

**Composition/monosaccharide shuffle** (multi-attribute variant): build a decoy glycan with
the *same number of fragment ions of each type* as its target but **randomly shifted masses
(1–20 Da/ion)**, or shift the intact glycan mass by a random value inside the mass tolerance
with a random isotope error — decoys keep the target's fragment *cardinality* but not its
masses ([Klein & Zaia 2022, *Mol Cell Proteomics* 21:100201, PMC8933705](https://pmc.ncbi.nlm.nih.gov/articles/PMC8933705/)).

## 4. Entrapment FDR (ground-truth validator)

Independent of decoys, a glyco search engine 2.0 validated calibration with **entrapment**: search a target
organism (yeast) against a database padded with a foreign glycome+proteome (mouse); any GPSM
whose glycan is mouse-only **or** whose peptide is mouse-only is a de-facto false positive.
The entrapment FDP then checks whether the decoy-estimated FDR is honest ([PMC5585273](https://pmc.ncbi.nlm.nih.gov/articles/PMC5585273/)).
For andes this is the correct kill-gate for the 2D-FDR — decoy FDR that under-covers the
entrapment FDP is not shippable.

## 5. Why the UNIFIED Percolator pile fails, and the separate-axis realization

Percolator learns **one** hyperplane over the PIN feature vector and one target/decoy split.
Spectrum-level glyco features (`OxoniumScore`, `YLadderScore`, `CoreYHits`, `GlycanMass`) are
*identical for the target peptide and a competing decoy peptide at the same backbone mass*
(SPA2 note) — they give zero peptide-axis discrimination, yet a unified glycan-decoy pile
makes Percolator weight them, sinking real targets. So the 2D-FDR **cannot** be a single
Percolator run; it must be a thin **post-process over Percolator PSM scores**.

**Percolator-only 2D-FDR procedure** (respects the FDR = Percolator-only constraint):

1. **Peptide axis (native Percolator).** Emit ONE PIN with reversed-peptide decoys only;
   `Label` = peptide target/decoy. Features = backbone b/y (`RankScore`/learned glyco `ScoreP`)
   **plus** the Y0/Y1 peptide-mass-anchor feature (discriminates competing peptides even when
   b/y is dead). Run Percolator → per-PSM `score_P`, `q_P`.
2. **Glycan axis (second native Percolator run).** Emit a SEPARATE PIN whose decoys are the
   §3 Y-rung-shifted glycan decoys (shift Y>Y1, keep Y0/Y1); `Label` = glycan target/decoy;
   features = glycan-structure evidence only (`YLadder`, `Oxonium`, `CoreYHits`). Run
   Percolator → `score_G`, `q_G`. Keeping the piles in two runs is what avoids the 29.4→4.4
   collapse.
3. **Combine per GPSM** on the sorted joint ranking (e.g. by `w·score_G + (1−w)·score_P`):
   at each cut compute `FDR^P`, `FDR^G`, and `FDR^{P∩G}` (GPSMs flagged decoy on *both*
   axes), then `FDR = FDR^P + FDR^G − FDR^{P∩G}`; accept GPSMs with `FDR ≤ 0.01`.
4. **Validate** with the §4 entrapment FDP before trusting the cut.

This keeps Percolator as the sole FDR engine (two vanilla runs + one inclusion–exclusion
merge), and differentiates andes — glycan-Y-first candidate selection, own learned `ScoreP`,
in-process cross-spectrum transfer — rather than cloning the reference engine/a comparison search engine.

## 6. Licenses (clean-room provenance)

- **a glyco search engine 2.0 / a glyco search engine** — algorithm published (PMC5585273); the **a glyco search engine binary is
  license-gated** (per-user application at i.pfind.net, *not* an OSI license). Use the
  **papers only**, never binary/decompiled logic.
- **a cross-spectrum glyco engine** — [github.com/DICP-1809/a cross-spectrum glyco engine](https://github.com/DICP-1809/a cross-spectrum glyco engine),
  **Apache-2.0** → clean-room reference OK (cross-spectrum transfer, random-shift decoys).
- **O-Pair / an open-source glyco engine** — [github.com/smith-chem-wisc/an open-source glyco engine](https://github.com/smith-chem-wisc/an open-source glyco engine),
  **MIT** → clean-room reference OK (graph localization, concatenated TDC q-values).

Sources: [a glyco search engine 2.0 PMC5585273](https://pmc.ncbi.nlm.nih.gov/articles/PMC5585273/) ·
[a cross-spectrum glyco engine PMC8990002](https://pmc.ncbi.nlm.nih.gov/articles/PMC8990002/) ·
[Multi-attribute glycan FDR PMC8933705](https://pmc.ncbi.nlm.nih.gov/articles/PMC8933705/) ·
[O-Pair *Nat Methods* PMC7606753](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7606753/).
