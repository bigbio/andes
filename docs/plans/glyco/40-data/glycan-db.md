# Standardized N-Glycan Database for andes (clean-room)

*Reference for the andes N-glycan composition DB (`andes-glyco/src/glycan_db.rs`, currently ≈600 common / 2510 full). Composition-only is sufficient for N-glyco search — no branching topology is needed, because scoring uses the intact-glycan mass, oxonium ions, and the Y-ladder (peptide + successive monosaccharide losses), all of which are composition-determined.*

## 1. Canonical N-glycan core & sequon

Every N-glycan is attached at the **N-X-S/T sequon** (X ≠ Pro; occasionally N-X-C) and shares the **common trimannosyl core (pentasaccharide) Man₃GlcNAc₂ = Hex₃HexNAc₂ ("M3")** ([Essentials of Glycobiology, NBK1917](https://www.ncbi.nlm.nih.gov/books/NBK1917/)). In H/N/F/S/G notation the core is **H3N2**. Three maturation classes, all built on this core:

- **High-mannose (oligomannose):** core + α-Man extensions only → **H5–H9 N2**, no Fuc/NeuAc/NeuGc (Man5…Man9).
- **Hybrid:** one arm mannose-terminated, one arm GlcNAc-initiated (has ≥1 antenna HexNAc but retained mannoses).
- **Complex:** both core mannoses substituted with GlcNAc antennae (bi-/tri-/tetra-/penta-antennary), elaborated with Gal (Hex), Fuc, and sialic acid (NeuAc/NeuGc).
- **Core fucosylation:** α1,6-Fuc on the innermost GlcNAc (adds one F); **sialylation:** NeuAc (S) caps galactose (human); NeuGc (G) is largely non-human/contaminant.

## 2. Composition alphabet & masses (validation targets)

Five monosaccharide building blocks; validate andes residue masses (monoisotopic, residue = monomer − H₂O) to <1 mDa:

| Symbol | Residue | Monoisotopic residue mass (Da) |
|---|---|---|
| H (Hex: Man/Gal/Glc) | C₆H₁₀O₅ | 162.052824 |
| N (HexNAc: GlcNAc/GalNAc) | C₈H₁₃NO₅ | 203.079373 |
| F (dHex: Fuc) | C₆H₁₀O₄ | 146.057909 |
| S (NeuAc / sialic acid) | C₁₁H₁₇NO₈ | 291.095417 |
| G (NeuGc) | C₁₁H₁₇NO₉ | 307.090331 |

Composition mass = Σ residues + H₂O (18.010565). E.g. core H3N2 = 3·162.0528 + 2·203.0794 + 18.0106 = **892.317** Da. **Validate the andes DB by checking every composition mass against a glyco search engine's canonical strings and GlyConnect composition masses** (both publish composition→mass), and by confirming the core M3 and the high-mannose series (Man5=H5N2 1216.42, Man9=H9N2 1864.63) land exactly.

## 3. Biosynthetic plausibility rules (clean-room, from a glyco search engine)

a glyco search engine builds its list by a **mammalian N-glycan biosynthetic simulation** rather than raw combinatorics, then simplifies structures to compositions ([Zeng et al., Nat Methods 2021, PMC8648562](https://pmc.ncbi.nlm.nih.gov/articles/PMC8648562/)). Reported constraints to reproduce (algorithm, not code — clean-room safe):

- **N (HexNAc): 2 ≤ N ≤ 8**; **H (Hex): 3 ≤ H ≤ 9** (up to 12 if allowing large paucimannose/hyper-mannose); **F ≤ 3** (F > 4 non-canonical); **S ≤ 4–5**; **G ≤ 2**.
- **Core floor:** H ≥ 3 and N ≥ 2 always (the trimannosyl core is obligatory).
- **Antenna coupling:** sialic acid ≤ available galactose ≤ antenna GlcNAc, i.e. roughly **S+G ≤ N−2** and Fuc ≤ N (andes already encodes `sialic ≤ hexnac−2`, `fuc ≤ hexnac`).
- **Oligomannose restriction:** for H > 5 with N = 2, retain only **N2H5–N2H9, no F/S/G** (high-mannose has no antennae, so no sialylation/fucosylation).
- Mass window ≈ **500–6000 Da** (andes uses [500,6000]).

andes' current rules (`glycan_db.rs`: N 2–8, H 3–12, F 0–3, S 0–5, G 0–2, fuc≤hexnac, sialic≤hexnac−2) are **consistent with a glyco search engine's** and produce 2510 compositions — the right order of magnitude.

## 4. Practical search-list sizes & sources

- **~180–200 compositions** — the reference glyco engine / FragPipe default N-glycan mass-offset list (mouse-derived, widely reused for human) ([Polasky et al., Nat Methods 2020](https://www.nature.com/articles/s41592-020-0967-9)). Fast, high-abundance-biased.
- **~600 common** — andes' curated "common human" tier; matches the abundant high-mannose + core-fuc + bi/tri-antennary sialylated space actually seen in HCD glyco datasets.
- **~1,000–1,234** — **GlyConnect** curated human compositions (1,041 compositions in the May-2020 release, [Compozitor, PMC8014996](https://pmc.ncbi.nlm.nih.gov/articles/PMC8014996/)) and **a glyco search engine's built-in list (1,234 compositions / 6,662 structures)**. These are the **gold curated references** — validate andes' common tier ⊆ GlyConnect ∪ a glyco search engine.
- **~2,500–8,000** — full biosynthetic-simulation space (andes 2510; larger if S/G/F ceilings relaxed). Higher coverage, larger false-match space.

**Sources for validation:** [GlyConnect / Compozitor (Expasy)](https://glyconnect.expasy.org/) and [GlyGen](https://www.glygen.org/) (aggregates GlyConnect/GlyTouCan; provides composition + monoisotopic mass + GlyTouCan accession per glycan); a glyco search engine's shipped `.gdb` canonical-string lists (Apache-2.0); UniCorn theoretical N-glycan database ([Akune et al.](https://www.sciencedirect.com/science/article/abs/pii/S0008621516301823)).

## 5. Recommended validation procedure for andes' clean-room DB

1. **Composition set:** assert andes' ~600 common tier is a subset of (GlyConnect human ∪ a glyco search engine built-in); flag any composition andes emits that neither curated source contains (likely biosynthetically implausible → tighten rules).
2. **Coverage:** assert every truth-set glycan in PXD025455 (523 scans) is present in andes' full 2510 list (a miss here is a hard generation gap, not a scoring gap).
3. **Mass parity:** for each shared composition, |mass_andes − mass_GlyConnect| < 1 mDa.
4. **Biosynthetic sanity:** no composition violates §3 (e.g. S>0 with N=2; F>N; H<3).

## Licensing (clean-room)

Databases/rules above are **descriptive facts and published algorithms**, freely usable. Reference implementations that are license-safe to read: **a glyco search engine (Apache-2.0)**, **a cross-spectrum glyco engine (Apache-2.0)**, **O-Pair / an open-source glyco engine (permissive)**. **Do NOT** copy from **a commercial glyco engine (commercial)** or **the reference glyco engine (UM-proprietary)** — their default glycan *lists* are just composition tables (facts, reusable), but their code is not. andes must generate its list from the biosynthetic rules, not import theirs.
