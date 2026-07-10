# andes glyco — differentiation strategy vs the field

*2026-07-02. Scope: intact N-glycopeptide search. Read `00-context/00-current-state.md`
first: andes already does ~90% backbone **generation**; it fails at **identification**
(ranking + Percolator-native FDR). This doc positions andes against the field and states
what NOT to build.*

## Comparison matrix

| Engine | Generation axis | Glycan handling | Peptide scoring | FDR | Cross-spectrum | License | Speed |
|---|---|---|---|---|---|---|---|
| **a glyco search engine** | Glycan-first (glycan ion-index) | DB (built-in glycan DB) | Fixed b/y + glycan match; sceHCD+ETxxD | Separate 2-D glycan/peptide FDR + a glyco search engineSite localization | No | **Apache-2.0** (repo) BUT runtime needs email-issued license code | 5–40× faster than peers (glycan ion-index) |
| **the reference glyco engine** | Open/mass-offset (peptide-first, fragment-index) | Mass-offset list / glyco DB; labile mode | Hyperscore on backbone b/y (glycan as delta) | Target-decoy; Philosopher/glyco-FDR | No | UM **academic / non-commercial**, commercial license required | Ultrafast (fragment-index open search) |
| **StrucGP** | DB, structure-first | Modular B/Y branch decomposition (core/subtype/branch) → structures | B/Y-pattern module matching | Target-decoy + structure scoring | No | Free binary, non-standard | Moderate |
| **a cross-spectrum glyco engine** | Glycan DB-independent peptide matching | Glycan-ladder stepping; unbiased glycan discovery | Backbone match + **spectrum expansion** | Own FDR | **Yes — same-backbone cross-spectrum transfer** (the +33.5–178.5% lever) | Free, **email-issued license code** (not OSI) | Moderate |
| **O-Pair (an open-source glyco engine)** | Ion-indexed open search (HCD) | Total glycan mass, then **graph localization** on EThcD | an open-source glyco engine scoring | Target-decoy (an open-source glyco engine) | No | **MIT** (open source) | >2000× faster localization vs prior O-glyco tools |
| **a commercial glyco engine** | DB + wildcard | Glycan DB + glycan wildcard (unanticipated glycans) | Proprietary intensity/fragment scoring | a commercial glyco engine 2-D score / PEP | No | **Commercial** (Protein Metrics) — closed | Moderate |
| **andes** | **Glycan-Y-first cascade + DB union** (~90% backbone) | Clean-room enumerator (2510 comps) + Y-ladder de-novo | **Own learned per-fragment models** (planned SP-B) | **Percolator-native** (thin 2-D post-process) | **Yes — in-process RT-gated transfer** (planned G4) | **Apache-2.0, pure Rust** | Fast (Rust fragment-index) |

## andes' defensible unique position

No competitor combines all four of andes' axes; each is individually precedented but the
**stack is unoccupied**:

1. **Glycan-Y-first candidate selection** — like a glyco search engine's glycan-first, but andes indexes on
   the intact-peptide Y-ladder (Y0/Y1 anchors), which is high-intensity in stepped-HCD even
   when backbone b/y is dead. Measured: backbone-findability 59.3 → 69.8% @0.05 Da (G1,
   `00-current-state.md`).
2. **Own learned per-fragment models** — the field uses fixed/hyperscore backbone matching
   (a glyco search engine, the reference engine, StrucGP) or proprietary (a commercial glyco engine). andes reuses its `andes train`
   partitioned model store to fit a **regime-matched glyco fragment model** (SP-B). This is
   the #1 identification lever: the SP-A2 baseline showed backbone-only scoring gives
   **0 IDs @1% FDR** because spectrum-level oxonium/YLadder features don't discriminate
   competing peptides at one backbone mass (`SPA2_RESULT.md`).
3. **Percolator-native 2-D FDR** — a glyco search engine does separate glycan/peptide FDR *inside* its own
   engine; andes gets the same statistics as a **thin Percolator post-process** (glycan axis
   and peptide axis rescored separately, then intersected). Hard constraint: Percolator only,
   never Mokapot; 2-D FDR is post-process, not a unified PIN pile (the unified pile crashed
   recovery 29.4 → 4.4%, G3).
4. **In-process RT-gated cross-spectrum transfer** — a cross-spectrum glyco engine's key result: direct b/y
   IDs the backbone in only a minority of spectra; its large gains came from **transferring
   backbone evidence across spectra sharing a peptide** (+33.5–178.5% glyco-PSMs vs a commercial glyco engine /
   the reference glyco engine / StrucGP / a glyco search engine). andes does this **in-process** (no separate tool,
   no re-search) and **RT-gated** to suppress false transfer (G4). This is the fix for the
   sparse-b/y stratum that no amount of single-spectrum scoring can rank.

**License moat:** andes is the only **Apache-2.0, pure-Rust, single-binary** engine in the set.
a glyco search engine's repo is Apache but the runtime still requires an emailed license code; a cross-spectrum glyco engine
is email-issued; the reference engine is UM non-commercial; a commercial glyco engine is closed commercial. Only O-Pair
(MIT, but C#/.NET, O-glyco-focused) is comparably open. andes is uniquely
redistributable/embeddable (e.g. quantms Nextflow) with no license gate.

## Clean-room provenance (hard constraint)

- **Reference OK (algorithms from published papers):** a glyco search engine (glycan-first, 2-D FDR),
  O-Pair (ion-index + graph localization — MIT, source-readable), a cross-spectrum glyco engine
  (spectrum-expansion transfer). Cite papers, re-derive from method text, do not vendor code.
- **Do NOT copy code:** a commercial glyco engine (commercial, closed) and the reference engine (UM-proprietary,
  non-commercial). Use their *published* behavior only as a benchmark target.
- **Note:** memory previously assumed a glyco search engine/a cross-spectrum glyco engine were freely-licensed source.
  Correction — both gate their **runtime** behind an issued license code; a glyco search engine's *repo* is
  Apache-2.0 (algorithms citable) but a cross-spectrum glyco engine's is an email-only license. Treat both as
  **paper-only** clean-room references, not code to import.

## What NOT to build

- **No O-Pair-style localization graph for single-sequon N-glyco.** O-Pair's graph theory
  solves *O*-glyco site ambiguity (many S/T, unknown occupancy). N-glyco here is single-sequon
  (N-X-S/T) — the site is determined by the sequon; a localization graph adds cost and buys
  nothing. (If multi-sequon N-glyco arises later, revisit — not now.)
- **No brute-force open / wildcard mass search.** the reference engine open-search and a commercial glyco engine wildcard
  already occupy that niche and it inflates the candidate space (SP-A2: ~20 (peptide,glycan)
  candidates/scan already; open search makes the needle-in-haystack worse, not better).
  andes' edge is the *constrained* glycan-Y-first cascade, not an unconstrained delta-mass sweep.
- **No re-implementation of the reference engine/a comparison search engine fragment-index open search.** andes differentiates
  by Y-first generation + learned scoring + cross-spectrum, not by being a faster open-search clone.

## Sources

- a glyco search engine — Nat Methods 2021, [10.1038/s41592-021-01306-0](https://www.nature.com/articles/s41592-021-01306-0);
  [PMC8648562](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8648562/); repo
  [github.com/pFindStudio/a glyco search engine](https://github.com/pFindStudio/a glyco search engine) (Apache-2.0 repo, license-code runtime).
- the reference glyco engine — Nat Methods 2020, [10.1038/s41592-020-0967-9](https://www.nature.com/articles/s41592-020-0967-9);
  [PMC7606558](https://pmc.ncbi.nlm.nih.gov/articles/PMC7606558/); UM
  [academic license](https://available-inventions.umich.edu/product/the reference engine-ultrafast-and-comprehensive-identification-of-peptides-from-tandem-mass-spectra).
- StrucGP — Nat Methods 2021, [10.1038/s41592-021-01209-0](https://www.nature.com/articles/s41592-021-01209-0).
- a cross-spectrum glyco engine — Nat Commun 2022, [10.1038/s41467-022-29530-y](https://www.nature.com/articles/s41467-022-29530-y);
  [PMC8990002](https://pmc.ncbi.nlm.nih.gov/articles/PMC8990002/); repo
  [github.com/DICP-1809/a cross-spectrum glyco engine](https://github.com/DICP-1809/a cross-spectrum glyco engine) (email-issued license).
- O-Pair / an open-source glyco engine — Nat Methods 2020, [10.1038/s41592-020-00985-5](https://www.nature.com/articles/s41592-020-00985-5);
  repo [github.com/smith-chem-wisc/an open-source glyco engine](https://github.com/smith-chem-wisc/an open-source glyco engine) (MIT).
- a commercial glyco engine — [Protein Metrics](https://www.proteinmetrics.com/products/a commercial glyco engine) (commercial);
  wildcard glyco [PMC8724605](https://pmc.ncbi.nlm.nih.gov/articles/PMC8724605/).
