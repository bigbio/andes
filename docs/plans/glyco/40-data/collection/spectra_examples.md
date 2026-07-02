# Annotated N-glycopeptide example spectra (REAL, HCD, Q Exactive Plus)

*Written 2026-07-02 for the andes glyco campaign. Purpose: give andes concrete,
ground-truthed stepped-HCD glycopeptide spectra to reason about — the fragmentation
regime the identifier (ranking + FDR) must survive. See
`../../00-context/00-current-state.md` for why generation is solved and
identification is the failure.*

## Provenance (this is REAL data, not synthetic)

| Field | Value |
|---|---|
| Dataset | **PXD030670** — "HILIC contributes to integrated N-glycoproteomics analysis of human saliva for lung cancer", Liu et al., *J Proteome Res* 2022;21(7):1589-1602 (doi:10.1021/acs.jproteome.1c00701) |
| Raw file | `HILIC-cancer-Intact_glycopeptides-3.raw` (885 MB, PRIDE FTP) |
| Instrument | **Thermo Q Exactive Plus Orbitrap**, HCD (from raw metadata) |
| Enrichment | ZIC-HILIC intact-glycopeptide enrichment |
| Ground truth | a commercial glyco engine 3.10.10 search table `HILIC-Intact_glycopeptides.xlsx` (309 mammalian N-glycans; 1% FDR, a commercial glyco engine score > 100, manually validated by the authors) |
| Conversion | ThermoRawFileParser 1.4.5 → centroided MGF (25,719 scans) |
| Matching | oxonium-gated MS2 (204.087 HexNAc > 2% base peak) → precursor neutral mass matched to a commercial glyco engine peptide+glycan within **≤15 ppm** (all 5 below ≤1.1 ppm) |

**Selection logic.** The a commercial glyco engine quant table gives peptide + glycan composition +
protein (the lowercase `n` marks the sequon Asn) but **no scan numbers**. I recomputed
each glycopeptide's neutral mass (bare-peptide monoisotopic mass, Cys carbamidomethyl
fixed, + glycan residue masses) and matched it to the nearest oxonium-positive MS2
precursor. 165 of 299 file-present glycopeptides matched an oxonium spectrum at ≤15 ppm;
the 5 below were chosen for **fragmentation diversity** (high-mannose, core-fucosylated
biantennary, sialylated, short peptide) and sub-ppm precursor accuracy.

### Reference masses used (monoisotopic, clean-room — from Unimod / textbook, no vendor code)

- Proton 1.007276; H2O 18.010565; Cys(+CAM) 160.03065.
- Glycan residue (dehydrated) masses: **HexNAc 203.079373, Hex 162.052824,
  Fuc 146.057909, NeuAc 291.095417, NeuGc 307.090331**.
- Oxonium reference ions (singly charged m/z): 126.055 / 138.055 / 144.065 / 168.066 /
  186.076 (HexNAc fragment series), **204.0867 HexNAc**, 163.0601 Hex, 274.0921 NeuAc−H2O,
  292.1027 NeuAc, **366.1395 HexNAc+Hex**, 512.1974 HexNAc+Hex+Fuc.
- Core-Y ladder = peptide(+H)_z + stepwise core: +HexNAc (Y1), +2HexNAc (Y2),
  then +Hex ×1..3 up the trimannosyl chitobiose core (Y0…Y core). Core-fucose adds
  146.0579 to Y1 (the diagnostic **peptide+HexNAc+Fuc** ion).

---

## Example 1 — IgG1 Fc, core-fucosylated (the canonical N297 glycopeptide)

**`EEQYNSTYR` + HexNAc(4)Hex(3)Fuc(1)** — IGHG1_HUMAN, sequon N-S-T.
scan 1633, z=2, precursor m/z 1317.5251, neutral 2633.039 (0.7 ppm), RT 8.2 min, 113 peaks.
Peptide backbone 1188.505; glycan 1444.534. This is the FA2G0 / "G0F"-type IgG Fc glycoform.

```
OXONIUM (relative to base peak = 204.087):
   204.0865  100.0%  HexNAc              <- base peak
   138.0548   74.0%  HexNAc fragment
   126.0550   38.7%  HexNAc fragment
   168.0654   28.5%  HexNAc fragment
   186.0759   24.2%  HexNAc-H2O
   366.1393   14.2%  HexNAc+Hex
   144.0654   12.2%  HexNAc fragment
   163.0598    0.9%  Hex

CORE-Y LADDER (peptide + core, singly charged unless noted):
   1392.5917   19.3%  Y1  pep+HexNAc              <- 2nd most intense glyco-Y
   1538.6498   12.7%  Y1F pep+HexNAc+Fuc          <- DIAGNOSTIC core-fucose ion
   1595.6826    3.6%  Y2  pep+HexNAc2
   1757.7147    2.7%  Y   pep+HexNAc2Hex
   1919.7667    2.7%  Y   pep+HexNAc2Hex2
   1189.5040    2.3%  Y0  peptide (bare)

BACKBONE b/y (bare peptide) — SPARSE:
   175.1187  y1;  338.1837 y2;  439.2303 y3;  526.2663 y4;  387.1521 b3   (all <4%)
   -> 5 backbone ions, none above 4% base peak.
```

**Annotation logic.** Base peak is 204.087 (HexNAc). The peptide+HexNAc (Y1, 1392.59)
and peptide+HexNAc+Fuc (Y1F, 1538.65) together pin the backbone mass AND prove the fucose
is on the core (Y1F present). The bare-peptide b/y ions exist but are all <4% — Percolator
cannot rank the peptide from these alone. This spectrum is why **backbone b/y is not a
discriminating feature** on glyco (SPA2_RESULT.md, Problem 2).

## Example 2 — IgG4 Fc, core-fucosylated (fuller Y-ladder)

**`EEQFNSTYR` + HexNAc(4)Hex(3)Fuc(1)** — IGHG4_HUMAN. scan 2554, z=3,
precursor m/z 873.3548, neutral 2617.044 (0.4 ppm), RT 12.6 min, 661 peaks.

```
OXONIUM:  204.087 100% (base); 138.055 83%; 126.055 60%; 168.066 28%; 186.076 24%;
          366.139 8.8%; 512.197 0.2% (HexNAc+Hex+Fuc, weak but present)
CORE-Y (this spectrum resolves the WHOLE core ladder, unusually complete):
   1376.598 16.9% Y1;  1522.653 8.1% Y1F(+Fuc);  1579.677 1.4% Y2;
   1741.731 1.2% +Hex; 1903.790 1.7% +Hex2; 2065.853 0.7% +Hex3 (trimannosyl core)
   (also the 2+ series: 688.80, 790.33, 871.37, 952.39, 1033.40)
BACKBONE b/y: b1 130.05, b2 259.09, b3 387.15, b4 534.22, b5 648.27;
              y1 175.12, y2 338.18, y3 439.24, y4 526.26, y7 915.43  -> 10 ions, all <2.1%
```

**Annotation logic.** At z=3 with 661 peaks this is a "good" glyco spectrum: the full
core-Y ladder (Y0→Y-core) is walkable at both 1+ and 2+, which is exactly what the andes
de-novo Y-solver + DB branch exploit (PHASE1_RESULT.md). Yet backbone b/y is still ≤2% —
richer peaks do NOT rescue backbone-based ranking.

## Example 3 — IgA1, high-mannose Man5 (backbone-rich outlier)

**`LAGKPTHVNVSVVMAEVDGTCY` + HexNAc(2)Hex(5)** — IGHA1_HUMAN. scan 13752, z=3,
precursor m/z 1188.5263, neutral 3562.558 (0.2 ppm), RT 61.4 min, 739 peaks.
Peptide backbone 2346.135 (long peptide); glycan 1216.423 (Man5 = HexNAc2Hex5).

```
OXONIUM:  204.087 100%; 138.055 99%; 168.066 41%; 366.139 24%; 163.060 27% (Hex
          elevated — high-mannose signature); 186.076 22%
CORE-Y (2+ ladder, strong):
   1275.623 40.4% Y1(2+);  1377.163 7.3% Y2;  1458.191 6.2% +Hex;
   1539.222 4.5% +Hex2;    1620.249 1.7% +Hex3 (core)
BACKBONE b/y — UNUSUALLY RICH: 25 ions matched
   b2 185.13, b4 370.25, b7 705.40, b8 804.47, b10 1017.58, b12 1203.69, b13 1302.75,
   b14 1433.79, b16 1633.87, b18 1847.96 ...
   y2 342.11 (18.6%!), y5 615.21, y6 714.28, y7 843.32, y9 1045.39, y10 1144.46, y12 1330.58
```

**Annotation logic.** The elevated 163.060 (Hex) and 366.139 vs a modest sialic signal
is the high-mannose fingerprint. This is the **backbone-rich tail of the distribution** —
a long, basic-residue-rich peptide fragments well even under HCD, giving 25 b/y ions and
a 18.6% y2. andes CAN rank spectra like this from backbone alone; the problem is they are
the minority. Keep this example as the "easy" contrast case.

## Example 4 — IgJ, sialylated biantennary (the HARDEST stratum)

**`ENISDPTSPLR` + HexNAc(4)Hex(5)NeuAc(2)** — IGJ_HUMAN. scan 2660, z=3,
precursor m/z 1145.1332, neutral 3432.382 (1.2 ppm), RT 13.1 min, 138 peaks.
Glycan 2204.772 (disialylated biantennary, A2G2S2).

```
OXONIUM:  204.087 100%; 138.055 98%; 126.055 43%; 168.066 39%; 366.139 22%; 186.076 26%;
          274.091 0.7% (NeuAc-H2O, weak) ; 292.103 essentially absent
CORE-Y LADDER:  *** NONE detected at >2% ***  (sialic-acid labile loss strips the antennae;
          the intact peptide+core Y ions did not survive)
BACKBONE b/y:  y1 175.12 ONLY (1.8%)  -> 1 ion total
```

**Annotation logic.** This is the worst case and the reason G4 (cross-spectrum transfer)
exists. Disialylated glycans lose NeuAc so readily under HCD that neither the Y-ladder
NOR the backbone survives: 1 backbone ion, no walkable Y-ladder. Single-spectrum scoring
**cannot** identify this — it is a cross-spectrum glyco engine's ~89% of spectra where direct b/y fails
(Nat Commun 2022). The precursor mass + oxonium prove it is a glycopeptide; only
composition-DB matching (SP-A DB branch) + cross-spectrum transfer can name it. Keep this
as the kill-case for any "just retrain the b/y model" proposal.

## Example 5 — short peptide, bisected/fucosylated

**`NLTATK` + HexNAc(5)Hex(3)Fuc(1)** — PCX1_HUMAN. scan 1073, z=2,
precursor m/z 1147.9959, neutral 2293.978 (0.4 ppm), RT 5.5 min, 152 peaks.
Peptide backbone 646.365 (only 6 residues); glycan 1647.613 (bisected + fucosylated).

```
OXONIUM:  204.087 100%; 138.055 62%; 126.055 40%; 168.066 25%; 186.076 24%; 366.139 7.3%
CORE-Y:   850.451 8.3% Y1;  996.504 2.7% Y1F(+Fuc);  647.372 7.3% Y0;  1053.535 3.0% Y2
BACKBONE b/y:  b2 228.13, b4 400.22, b5 501.27; y1 147.11, y2 248.16  -> 5 ions, all <2%
```

**Annotation logic.** Short peptide → few possible b/y ions in the first place, and they
are weak. Y0 (bare peptide, 647.37) is actually visible here because the peptide is light
enough to fall in a well-populated m/z region. The Y1/Y1F pair again localizes the core
fucose. Illustrates that peptide length bounds backbone evidence independently of the
fragmentation regime.

---

## Fragmentation pattern — the three takeaways for andes

1. **Oxonium-dominated, universally.** 204.087 (HexNAc) is the base peak in all 5, with the
   138/126/168/186 HexNAc-fragment satellites and 366.139 (HexNAc+Hex) always present.
   Oxonium fires in 16,048 / 25,719 MS2 in this file — a reliable glyco GATE but
   **spectrum-level, so it cannot separate competing peptides at one backbone mass**
   (the exact null-discrimination trap of SPA2_RESULT.md Problem 2 and the refuted unified
   Percolator pile of the current-state G3 note).

2. **The core-Y ladder is the real backbone anchor, and its completeness is glycan-dependent.**
   Neutral/paucimannose and high-mannose (Ex 1-3) give a partial-to-full peptide+HexNAc…core
   ladder including the **diagnostic peptide+HexNAc+Fuc** ion that localizes core fucose;
   **sialylated glycans (Ex 4) destroy it** via labile NeuAc loss. This is why SP-B should be
   built as a *peptide-mass-anchored Y0/Y1 feature* (current-state §6), not a b/y model:
   Y1 is high-intensity even when b/y is dead.

3. **Backbone b/y is sparse and never intense — exactly as predicted.** 4/5 examples give
   ≤10 backbone ions, all <4% of base peak; only the long basic peptide (Ex 3, 25 ions) is
   backbone-rich, and the disialylated case (Ex 4) gives **one** ion. Backbone-only RankScore
   therefore cannot rank glyco PSMs (0 @1% FDR baseline). The lever is not a better b/y model
   over the same dead peaks but **(a)** a peptide-anchored Y-feature (SP-B) and **(b)**
   cross-spectrum transfer for the sialylated/short-peptide strata (G4).

## Clean-room / license note (respect for all downstream reuse)

- Oxonium m/z, glycan residue masses, and the core-Y logic here are from **published
  literature and Unimod**, not from any engine's source. No a commercial glyco engine (commercial, Protein
  Metrics) or the reference engine (UM-proprietary) code was read or copied — only a commercial glyco engine's *output
  table* was used as ground-truth labels, which is fair use of a public PRIDE result file.
- Clean-room algorithmic references for the identification work these examples motivate:
  **a glyco search engine / a cross-spectrum glyco engine** (Zeng et al., *Nat Commun* 2021/2022; Apache-2.0,
  github.com/pFindStudio/a glyco search engine, github.com/pFindStudio/a cross-spectrum glyco engine) for the separate
  glycan/peptide 2D-FDR and cross-spectrum transfer; **O-Pair / an open-source glyco engine** (Lu et al.,
  *Nat Methods* 2020; MIT-licensed, github.com/smith-chem-wisc/an open-source glyco engine) for graph
  localization. FDR remains **Percolator-only** with 2D-FDR as a thin separate-axis
  post-process (never a unified pile, never Mokapot) per the andes feedback constraints.

## Reproduce

```
# 1. truth table (52 KB)
curl -s 'ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2023/02/PXD030670/HILIC-Intact_glycopeptides.xlsx' -o truth.xlsx
# 2. raw (885 MB) + convert
curl -s 'ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2023/02/PXD030670/HILIC-cancer-Intact_glycopeptides-3.raw' -o cancer3.raw
ThermoRawFileParser -i cancer3.raw -o . -f 0 -m 0     # -> cancer3.mgf (centroided MGF)
# 3. annotate scans 1633, 2554, 13752, 2660, 1073 with the masses above.
```

Annotation scripts used this session (not committed; in the session scratchpad):
`annotate.py` (oxonium-gate + ppm precursor match), `fullannot.py` (per-scan oxonium /
core-Y / b/y annotation).
