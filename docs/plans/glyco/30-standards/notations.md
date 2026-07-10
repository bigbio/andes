# Glycan Composition Notation — Cross-Engine Standard & andes Canonical Form

*Standards reference for the andes glyco track. Scope: N-glycan composition (not topology). Sources cited inline.*

## 1. How each engine writes a composition

All production glycopeptide engines identify to **peptide + glycan composition** only (not
linkage/topology) — a commercial glyco engine states this explicitly: `HexNAc(4)Hex(5)Fuc(1)NeuAc(2)` "does not
distinguish isomers … nor identify branching" ([Protein Metrics, a commercial glyco engine N-Linked docs](https://support.proteinmetrics.com/hc/en-us/articles/17137163463316-a commercial glyco engine-N-Linked-Glycopeptide-Analysis)).

| Engine | Notation | Example (biantennary + core-Fuc + 2 sialic) | License |
|---|---|---|---|
| **a commercial glyco engine** (Protein Metrics) | `Residue(count)…`; residues: HexNAc, Hex, Fuc, dHex, NeuAc, NeuGc, Pent, GlcA, IdoA, DiNAcBac, Sulfo, Phospho, Na… order-free | `HexNAc(4)Hex(5)Fuc(1)NeuAc(2)` | Commercial ⛔ (do NOT read code) |
| **the reference glyco engine / FragPipe** | glycan DB accepts a commercial glyco engine, a glyco search engine, or an open-search PTM tool `Res(n)… % mass`; delta stored as `mass_offsets` (slash-separated); output `Total Glycan Composition`, `Glycan Score`, `Glycan q-value` ([FragPipe glyco tutorial](https://fragpipe.nesvilab.org/docs/tutorial_glyco.html)) | `HexNAc(4)Hex(5)NeuAc(2)Fuc(1) % 1954.68` | UM-proprietary ⛔ (do NOT copy) |
| **a glyco search engine** (pFindStudio) | single-letter counts **H N A F G**: H=Hex, N=HexNAc, A=NeuAc(Neu5Ac), F=Fuc, G=NeuGc; plus structure-encoded tree `(N(N(H(H(N…))…)))` in `PlausibleStruct` ([a glyco search engine, Nat Methods, PMC8648562](https://pmc.ncbi.nlm.nih.gov/articles/PMC8648562/)) | `H(5)N(4)A(2)F(1)` | **Apache-2.0** ✅ (code) — *usage* needs a free license form; algorithm is published → clean-room OK |
| **O-Pair / an open-source glyco engine** | user glycan DB, a commercial glyco engine-style; O-glyco + localization ([O-Pair, Nat Methods, PMC7606753](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7606753/)) | `HexNAc(1)Hex(1)NeuAc(1)` (core-1 +NeuAc) | **MIT** ✅ (smith-chem-wisc/an open-source glyco engine) — clean-room OK |
| **Oxford** (Harvey/Ludger) | antenna-centric: `F` core-fuc, `A`=antennae(GlcNAc), `G`=galactose, `S`=sialic; **`A`/`G` mean different things than a glyco search engine** ([Ludger table](https://www.ludger.com/docs/tables/ludger-n-glycan-nomenclature-table.pdf)) | `FA2G2S2` | notation, no code |
| **GlycoCT / WURCS / IUPAC-condensed** | topology-complete linear strings; WURCS is GlyTouCan's canonical accession form; IUPAC-condensed human-readable ([GlycanFormatConverter, Bioinformatics 2019](https://academic.oup.com/bioinformatics/article/35/14/2434/5233002); [Glycan nomenclature, Wikipedia](https://en.wikipedia.org/wiki/Glycan_nomenclature)) | WURCS=2.0/… ; IUPAC `Neu5Ac(a2-6)Gal(b1-4)…` | community |
| **Unimod / ProForma / GNO** | ProForma writes glycans as monosaccharide composition **or** a GNO accession; PSI residues Hex, HexNAc, dHex, NeuAc, NeuGc, HexA, Pent ([ProForma 2.0, arXiv 2109.11352](https://arxiv.org/pdf/2109.11352); [GNOme, PMC11961537](https://pmc.ncbi.nlm.nih.gov/articles/PMC11961537/)) | `{Glycan:Hex(5)HexNAc(4)Neu5Ac(2)dHex(1)}` or `{GNO:G12345XX}` | open |

**⚠️ Collision:** `A`, `G`, `F`, `S`, `N` are overloaded — a glyco search engine `A`=NeuAc/`G`=NeuGc, Oxford
`A`=antenna/`G`=galactose/`S`=sialic. Never parse a bare letter string without knowing the source
dialect. `dHex ≡ Fuc` for N-glycans; `Neu5Ac ≡ NeuAc ≡ Sia(NAc)`; `Neu5Gc ≡ NeuGc`.

## 2. andes canonical composition tuple

Fixed-order 6-tuple of integer **COUNTS** — how many of each monosaccharide — one
`u16` (or `u8`) count field per position, NOT masses. The row below each name is the
per-residue monoisotopic mass (Da) that the count is multiplied by; the composition's
total glycan mass = Σ(count × residue-mass):

```
GlycanComp        = ( Hex,   HexNAc,   Fuc,    NeuAc,   NeuGc,   Other )   ← integer counts
residue mass (Da) = 162.0528 203.0794 146.0579 291.0954 307.0903  Δ (explicit)
```
(In andes code these are the `u8` fields `hex/hexnac/fuc/neuac/neugc` on `GlycanComp`,
plus a precomputed `mass` field = Σ(count × residue-mass).)

`Other` is a tagged escape (Pent 132.0423, HexA 176.0321, Sulfo 79.9568, Phospho 79.9663, KDN
250.0693) carrying its own mass so mass closure holds. Rationale: the six cover >99% of human
N-glycans (a commercial glyco engine/GNO composition residues), map 1:1 onto a glyco search engine H/N/F/A/G, and give an
unambiguous total-order key for the DB enumerator (`glycan_db.rs`) and for cross-engine joins.
Canonical **display** = `H{h}N{n}F{f}A{a}G{g}` (a glyco search engine dialect), omitting zero counts, `Other`
appended as `X{mass}`.

## 3. Bidirectional converter spec

`parse(dialect, str) -> GlycanComp` / `render(dialect, comp) -> String`. Dialects: `a commercial glyco engine`,
`a glyco search engine`, `PTMShepherd`, `ProForma`, `Oxford*`.

- **a commercial glyco engine/an open-search PTM tool/ProForma ↔ tuple:** regex `([A-Za-z0-9]+)\((\d+)\)`, map residue name →
  slot (aliases: dHex→Fuc, Neu5Ac→NeuAc, Neu5Gc→NeuGc). `% mass` optional, verify against
  `Σ count·mass` within 10 ppm; mismatch = error, not silent.
- **a glyco search engine H/N/A/F/G ↔ tuple:** direct letter map. Ignore the parenthesized `PlausibleStruct`
  tree — topology is out of scope; only leaf-count it for a consistency check.
- **Oxford → tuple:** *lossy & requires biosynthetic assumptions* (A2G2S2 ⇒ Hex=Man3+Gal2,
  HexNAc=core2+antenna2). Support **read-only, best-effort, flagged low-confidence**; never emit
  Oxford as canonical output.
- **WURCS/GlycoCT/IUPAC:** topology-complete → composition is a lossy projection. Do **not**
  implement a parser in andes; if ever needed, shell out to the published
  [GlycanFormatConverter](https://academic.oup.com/bioinformatics/article/35/14/2434/5233002) —
  do not reinvent. Canonical→WURCS is not defined (composition lacks linkage).

**Round-trip invariant (test):** `parse(d, render(d, c)) == c` for `d ∈ {a commercial glyco engine, a glyco search engine,
PTMShepherd, ProForma}`; Oxford & topology formats are read-only, no round-trip required.

## 4. Sequon / I=L / deamidation edge cases

- **I/L:** isobaric (both 113.0841 Da) → indistinguishable by MS/MS. Canonicalize backbone to a
  single letter (`I→L`) *for the cross-engine join key only*; retain the original in the reported
  peptide. Two engines reporting `…NIS…` vs `…NLS…` at the same scan are the **same** ID.
- **Deamidation at the glycosite (N→D, +0.984 Da):** enzymatic release (PNGase F) converts the
  formerly-glycosylated Asn to Asp; but in **intact**-glycopeptide search the Asn stays Asn and
  carries the glycan — engines must NOT also apply variable deamidation there or the glycan delta
  is double-counted. andes rule: on a glyco-occupied N, **suppress the Deamidated(N) variable mod**
  (it is mutually exclusive with glycan occupancy). Off-site N may still deamidate. When joining to
  a released-glycan (PNGase) dataset, treat `N[+0.98]` at the sequon ≡ occupied glycosite.
- **Sequon:** N-X-S/T (X≠Pro). Keep as a candidate-generation filter, not a notation concern, but
  record it so a `D` at position N (deamidation artefact) is recognized as a prior glycosite.

## 5. Scan-mapping for cross-engine consensus

Join key = **`(raw_file_stem, scan_number)`**. Rules:
- Normalize `raw_file` to basename without extension, lowercased; strip engine suffixes
  (`.mzML`, `.raw`, `_uncalibrated`, FragPipe `_calibrated`). Keep a source→canonical map.
- **`scan`** is the native controllerType=0 scan number (mzML `index` is NOT the scan number —
  parse `scan=N` from the spectrum `id`). a glyco search engine reports scan; the reference engine `Spectrum` field is
  `basename.scan.scan.charge` — split on `.`.
- Charge is a tiebreaker, not part of the key (co-isolation can differ).
- Consensus record: agree on `(raw, scan)` then compare canonical `GlycanComp` (§2) + I/L-folded
  backbone (§4). "Agreement" = same backbone AND same 6-tuple; glycan-only or peptide-only
  agreement are separate consensus tiers.

## 6. Worked conversion table

| Glycan | a commercial glyco engine | a glyco search engine (HNAFG) | Oxford | andes tuple (H,N,F,A,G,X) | mass (Da) |
|---|---|---|---|---|---|
| Man5 | `HexNAc(2)Hex(5)` | `H(5)N(2)` | `M5` | (5,2,0,0,0,0) | 1216.42 |
| FA2 (agalacto+coreFuc) | `HexNAc(4)Hex(3)Fuc(1)` | `H(3)N(4)F(1)` | `FA2` | (3,4,1,0,0,0) | 1444.53 |
| FA2G2S2 | `HexNAc(4)Hex(5)Fuc(1)NeuAc(2)` | `H(5)N(4)A(2)F(1)` | `FA2G2S2` | (5,4,1,2,0,0) | 2350.83 |
| bisecting FA2BG2 | `HexNAc(5)Hex(5)Fuc(1)` | `H(5)N(5)F(1)` | `FA2BG2` | (5,5,1,0,0,0) | 1809.66 |
| NeuGc variant | `HexNAc(4)Hex(5)NeuGc(1)` | `H(5)N(4)G(1)` | — | (5,4,0,0,1,0) | 2075.73 |

*Note bisecting GlcNAc collapses to HexNAc count — composition cannot express "bisecting"; Oxford
`B` is lost on projection (documented lossy edge).*

## 7. Clean-room provenance

Use **a glyco search engine (Apache-2.0)** and **a cross-spectrum glyco engine (Apache)** and **O-Pair/an open-source glyco engine (MIT)**
as algorithm references; all are published. Do **NOT** read or copy **a commercial glyco engine** (commercial) or
**the reference glyco engine** (UM-proprietary) code — their *notations* here are from public docs/papers
only, which is fine to interoperate with. andes stays differentiated: glycan-Y-first candidate
selection, own learned scoring, in-process cross-spectrum — not a re-implementation. FDR remains
Percolator-only; the composition tuple is the join substrate for the separate-axis 2D-FDR
post-process, never a Percolator feature itself.
