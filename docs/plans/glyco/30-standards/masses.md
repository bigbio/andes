# Glyco standardized mass reference (andes)

*Authoritative monoisotopic mass conventions for N-/O-glycopeptide search. All values
cross-checked against Unimod, ProForma 2.0, a glyco search engine, FragPipe/the reference glyco engine, and
GlycoMod. Computed from CODATA/AME atomic masses (H 1.00782503, C 12 exact, N 14.00307401,
O 15.99491462, P 30.97376151, S 31.97207069); proton = 1.00727647.*

## 1. Constants

| Constant | Value (Da) | Note |
|---|---|---|
| Proton (charge carrier) | **1.007276** | mass of H⁺ = H − electron; use for m/z, NOT 1.00783 |
| Water (H₂O) | **18.010565** | |
| Electron | 0.000549 | rarely needed explicitly |

## 2. Monosaccharide **residue** masses (in-glycan = free monosaccharide − H₂O)

Glycans polymerize by glycosidic condensation, so each residue in a chain contributes
its free mass minus one water. These are the **residue** masses andes must sum.

| Monosaccharide | Symbol | Residue formula | Residue mass (Da) | Free mass (+H₂O) |
|---|---|---|---:|---:|
| Hexose (Man/Gal/Glc) | Hex | C₆H₁₀O₅ | **162.052824** | 180.063388 |
| N-acetylhexosamine (GlcNAc/GalNAc) | HexNAc | C₈H₁₃NO₅ | **203.079373** | 221.089938 |
| Deoxyhexose (Fucose) | dHex / Fuc | C₆H₁₀O₄ | **146.057909** | 164.058473 |
| N-acetylneuraminic acid (sialic) | NeuAc / Sia | C₁₁H₁₇NO₈ | **291.095417** | 309.105982 |
| N-glycolylneuraminic acid | NeuGc | C₁₁H₁₇NO₉ | **307.090331** | 325.100896 |
| Pentose (Xyl/Ara) | Pent | C₅H₈O₄ | **132.042259** | 150.052823 |
| Hexuronic acid (GlcA/IdoA) | HexA | C₆H₈O₆ | **176.032088** | 194.042653 |
| Ketodeoxynonulosonic acid | Kdn | C₉H₁₄O₈ | **250.068867** | 268.079432 |
| Phosphorylation (on glycan) | Phospho | +HPO₃ | **+79.966331** | additive delta |
| Sulfation (on glycan) | Sulfo | +SO₃ | **+79.956815** | additive delta |

These match Unimod/ProForma 2.0 (Hex 162.0528, HexNAc 203.0793, dHex 146.0579,
NeuAc 291.0954, Pent 132.0422) and a glyco search engine's `glycan.ini` values
(203.07937, 162.05282, 146.05791, 291.09542, 307.09033) to 5 decimals.
FragPipe/the reference glyco engine explicitly uses **203.07937** as the HexNAc "remainder" residue mass.

## 3. Glycan composition mass

A glycan composition **G** = {nᵢ × residueᵢ} contributes, when **attached to Asn/Ser/Thr**,
the pure residue sum with **no extra water** (the glycosidic bond to the peptide is itself a
condensation; the reducing-end anomeric OH is consumed in the N-/O-glycosidic linkage):

```
glycan_residue_sum(G) = Σ nᵢ · residue_massᵢ          # attached form — USE THIS
```

Example, biantennary complex **HexNAc(4)Hex(5)NeuAc(2)Fuc(1)** ("A2G2S2F"):
`4·203.079373 + 5·162.052824 + 2·291.095417 + 1·146.057909 = 2350.831693 Da`.

A **free/released** glycan (reducing end unlinked, e.g. GlycoMod, PNGase-F released) adds one
water: `free_glycan = glycan_residue_sum + 18.010565`. GlycoMod/permethylation tools work in
the free-reducing-end convention; **andes searching intact glycopeptides must NOT add this water.**

## 4. Glycopeptide precursor (neutral monoisotopic) mass

```
M_glycopeptide = M_peptide + glycan_residue_sum(G)
```
where `M_peptide` = Σ residue masses + H₂O + fixed/variable peptide mods (Cam-C, etc.).
The glycan is a single delta mass on the sequon residue; **no second water is added** — this
is the andes/the reference engine "mass-offset" convention. Precursor m/z at charge z:
`(M_glycopeptide + z·1.007276) / z`.

## 5. Y-ion (glycan-retaining backbone) conventions

Y ions = intact peptide backbone + partial glycan (peptide fragmentation absent). Neutral
Y-species mass = `M_peptide + (partial glycan residue sum)`; observed singly-charged m/z adds
one proton. Canonical N-glycan core ladder (a glyco search engine/a commercial glyco engine/O-Pair agree):

| Y label | Composition on peptide | Neutral mass |
|---|---|---|
| **Y0** | peptide only | `M_peptide` |
| Y1 | + HexNAc | `M_peptide + 203.079373` |
| Y2 | + 2 HexNAc | `M_peptide + 406.158746` |
| Y3 | + 2 HexNAc + Hex | `+ 568.211570` |
| Y4/Y5 | + HexNAc₂Hex₂ / HexNAc₂Hex₃ | trimannosyl-core rungs |

**Y0 observed (1+) = M_peptide + 1.007276**; Y1 (1+) = M_peptide + 203.079373 + 1.007276.
Y0/Y1 are the high-intensity peptide-mass anchors (the SP-B anchor feature) — present even
when b/y is dead. Multiply-charged Y: `(neutral + z·1.007276)/z`.

## 6. Divergence flags for andes' clean-room DB

- **Water double-count (highest risk).** If the glycan enumerator stores *free* monosaccharide
  masses or a *free-glycan* total, `precursor − glycan` will be off by −18.0106 Da per glycan.
  Enforce **residue masses, one water total, on the peptide only.** (PHASE1 notes already flag a
  "H2O convention" fix — pin it to this doc.)
- **Proton vs H mass.** Use 1.007276 for all charge carriers/m/z. Using 1.00783 (neutral H)
  shifts every Y/precursor m/z by ~0.5 mDa/charge — silent low-grade calibration drift.
- **HexNAc rounding.** 203.0793 vs 203.07937: 7 mDa — within 20 ppm at backbone mass but can
  flip near-isobaric HexNAc₄ ≈ Hex₅ (+2 Da) calls; store full 6-decimal `203.079373`
  (see ppmFixer, *Glycobiology* 2024, on a glyco search engine near-isobaric mismatches).
- **NeuGc/Kdn/HexA/sulfo/phospho** must exist in the DB for non-human, sulfated, or
  phosphomannose glycans; absence silently drops those glycopeptides (not an FDR issue —
  a coverage hole).
- **Sequon-only attachment.** The −H₂O residue convention is only valid for the *attached*
  form; any released-glycan tooling andes reuses (GlycoMod-style) is in the +H₂O frame —
  never mix the two in one code path.

## Sources
- Unimod / ProForma 2.0 (Hex 162.0528, HexNAc 203.0793, dHex 146.0579, NeuAc 291.0954, Pent 132.0422): arXiv:2109.11352 (PSI ProForma 2.0).
- a glyco search engine (glycan-first, clean-room Apache-family reference): PMC4853738 (a glyco search engine), a glyco search engine (*Nat Methods* 2021); ppmFixer near-isobaric note: *Glycobiology* 34(4) cwae006 (2024).
- FragPipe/the reference glyco engine residue definitions (HexNAc remainder 203.07937): fragpipe.nesvilab.org/docs/tutorial_glyco.html (UM-proprietary — **spec/values only, do not copy code**).
- a commercial glyco engine N-linked conventions (commercial — reference only, do not copy): support.proteinmetrics.com.
- O-Pair / an open-source glyco engine (permissive, clean-room OK): Y-ion + core-ladder conventions.
- GlycoMod (free-reducing-end +H₂O frame contrast): web.expasy.org/glycomod/glycomod-doc.html.

*Licenses: a glyco search engine (academic/open), O-Pair-an open-source glyco engine (permissive) = clean-room-safe to
reimplement from paper. a commercial glyco engine (commercial) and the reference engine (UM-proprietary) = read published
values only; never transcribe source. All masses here are physical constants / published
values, not code.*
