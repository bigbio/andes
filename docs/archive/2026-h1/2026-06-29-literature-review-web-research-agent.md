# Web-Research Agent Brief — DDA Scoring Literature (≥55 DOIs)

**Parent:** [`2026-06-29-literature-review-brief.md`](2026-06-29-literature-review-brief.md)  
**Mission:** Adversarially verify patents, licenses, and low-res claims. Acquire missing PDFs. Do **not** re-derive math already in parent doc.

**Success criterion for the engine:** PSMs at **1% paired entrapment-FDP** (Wen & Noble 2025), not reported TDC FDR.

---

## Search protocol (mandatory)

For each DOI below:
1. Fetch abstract + methods equations (not just title).
2. Record: **low-res tested?** (Y/N/uncertain), **patent?**, **code license**, **weights license**.
3. Cross-check patent numbers on Google Patents — never assert freedom without citation.
4. Flag `[ACQUIRE]` if paywalled and not in [`internal-docs/papers/`](internal-docs/papers/REFERENCES.md).

**Off-limits (do not recommend implementing):** US **8,639,447** generating-function spectral probability (active to **2030-07-25**).

---

## Tier A — Must read (calibration + low-res gap)

| # | DOI | Why search |
|---|-----|------------|
| A1 | https://doi.org/10.1021/acs.jproteome.9b00736 | **Tailor** — per-spectrum Q99 calibration; low-res benchmark vs exact p |
| A2 | https://doi.org/10.1074/mcp.O113.036327 | **Exact XCorr p-value** via DP — patent-adjacent pattern |
| A3 | https://doi.org/10.1021/pr5010983 | **Keich & Noble** — Monte Carlo 10K decoy DBs per spectrum for empirical p |
| A4 | https://doi.org/10.1021/acs.jproteome.8b00206 | **Res-ev + combined p-value** — high-res; contrast with low-res XCorr arm |
| A5 | https://doi.org/10.1021/pr8001244 | MS-GF **spectral probability** (patented implementation) — understand only |
| A6 | https://doi.org/10.1038/ncomms6277 | MS-GF+ rank score + GF E-value |
| A7 | https://doi.org/10.1021/pr0499491 | **OMSSA** Poisson E-value — ion-trap motivation |
| A8 | https://doi.org/10.1021/pr101065j | **Andromeda** binomial peak-depth — local 100 Th windows |
| A9 | https://doi.org/10.1038/s41592-025-02719-x | **Entrapment FDP** + paired estimator (r=1) |
| A10 | https://doi.org/10.1021/pr8011107 | Statistical XCorr calibration (parametric precursor) |

---

## Tier B — Primary scoring functions

| # | DOI | Topic |
|---|-----|-------|
| B1 | https://doi.org/10.1016/1044-0305(94)80016-2 | SEQUEST XCorr original |
| B2 | https://doi.org/10.1021/pr800420s | Fast XCorr / background subtraction |
| B3 | https://doi.org/10.1093/bioinformatics/bth023 | X!Tandem hyperscore |
| B4 | https://doi.org/10.1021/ac025676e | Hyperscore E-value / survival function |
| B5 | https://doi.org/10.1038/nmeth.4256 | MSFragger hyperscore + index (license!) |
| B6 | https://doi.org/10.1021/acs.jproteome.3c00486 | Sage (MIT) |
| B7 | https://doi.org/10.1021/pr101196n | Tide / fast SEQUEST reimplementation |
| B8 | https://doi.org/10.1186/1471-2105-8-327 | MyriMatch hypergeometric |
| B9 | https://doi.org/10.1021/pr8007374 | Frank **rank prediction** |
| B10 | https://doi.org/10.1021/pr8006788 | Frank **rank-based PSM score** |
| B11 | https://doi.org/10.1021/acs.jproteome.6b00290 | DRIP DBN alignment |
| B12 | https://doi.org/10.1093/bioinformatics/btn189 | Klammer DBN fragmentation model |
| B13 | https://doi.org/10.1074/mcp.M700022-MCP200 | Morpheus (dot-product, high-res) |
| B14 | https://doi.org/10.1021/pr300631t | De-Noise post-processor (ion-trap LCQ/LTQ) |
| B15 | https://doi.org/10.1021/pr401026y | Empirical multidimensional PSM scoring (Zubarev) |

---

## Tier C — FDR, validation, inference

| # | DOI | Topic |
|---|-----|-------|
| C1 | https://doi.org/10.1021/pr700600n | Target-decoy (Käll) |
| C2 | https://doi.org/10.1038/nmeth1113 | Percolator |
| C3 | https://doi.org/10.1074/mcp.T900012-MCP200 | PeptideProphet |
| C4 | https://doi.org/10.1021/ac0341261 | ProteinProphet |
| C5 | https://doi.org/10.1074/mcp.M900317-MCP200 | MAYU protein FDR |
| C6 | https://doi.org/10.1002/pmic.201500431 | Protein-level FDR semantics |
| C7 | https://doi.org/10.1016/j.jprot.2010.08.009 | Nesvizhskii survey (error rates) |
| C8 | https://doi.org/10.1002/rcm.4417 | Empirical FDR estimation (Goloborodko) |
| C9 | https://doi.org/10.1101/2024.06.01.596967 | Entrapment preprint (same as A9) |
| C10 | https://doi.org/10.1038/s41467-025-58728-z | MSFragger-DDA+ entrapment validation |

---

## Tier D — Candidate generation & search space

| # | DOI | Topic |
|---|-----|-------|
| D1 | https://doi.org/10.1038/nmeth.1889 | InsPecT tags |
| D2 | https://doi.org/10.1021/pr0500111 | Peptide sequence tags (Frank) |
| D3 | https://doi.org/10.1021/ac048788h | PepNovo |
| D4 | https://doi.org/10.1038/s41467-024-49731-x | Casanovo (Apache weights) |
| D5 | https://doi.org/10.1038/s41587-024-01382-9 | InstaNovo (NC weights) |
| D6 | https://doi.org/10.1021/pr101196n | Tide indexing |
| D7 | PMC13232765 | Comet fragment-ion index (2025) — https://pmc.ncbi.nlm.nih.gov/articles/PMC13232765/ |
| D8 | https://doi.org/10.1021/pr8001244 | Spectral dictionaries / tag prefilters |
| D9 | https://doi.org/10.1089/cmb.2014.0165 | GF on spectral networks (related patent family) |
| D10 | https://doi.org/10.1074/mcp.M110.003731 | Crux cascaded search |

---

## Tier E — Spectrum processing & libraries

| # | DOI | Topic |
|---|-----|-------|
| E1 | https://doi.org/10.1002/pmic.200600625 | SpectraST library cosine |
| E2 | https://doi.org/10.1038/nmeth.1240 | SpectraST consensus libraries |
| E3 | https://doi.org/10.1021/pr900473s | SpectraST decoy libraries |
| E4 | https://doi.org/10.1371/journal.pcbi.1008724 | Spec2Vec embeddings |
| E5 | https://doi.org/10.1021/pr0700693 | Peptide-centric analysis (Ting) |
| E6 | https://doi.org/10.1038/s41467-020-18138-8 | Quandenser quant-first |
| E7 | https://doi.org/10.1021/pr0709777 | Spectral counting / intensity transforms |
| E8 | https://doi.org/10.1021/pr0601336 | Dynamic exclusion / peak picking context |

---

## Tier F — Learned prediction & rescoring

| # | DOI | Topic |
|---|-----|-------|
| F1 | https://doi.org/10.1038/s41592-019-0426-7 | Prosit |
| F2 | https://doi.org/10.1093/nar/gkz299 | MS2PIP |
| F3 | https://doi.org/10.1038/nbt.4313 | pDeep |
| F4 | https://doi.org/10.1038/s41467-023-40129-9 | MSBooster |
| F5 | https://doi.org/10.1016/j.mcpro.2022.100266 | MS²Rescore immunopeptidomics |
| F6 | https://doi.org/10.1021/acs.jproteome.3c00785 | MS²Rescore 3.0 modular |
| F7 | https://doi.org/10.1074/mcp.M900222-MCP200 | Elude RT (Percolator family) |
| F8 | https://doi.org/10.1021/pr9009233 | mProphet / RT features |
| F9 | https://doi.org/10.1038/nmeth.1711 | PeptideAtlas / training data context |
| F10 | https://doi.org/10.1021/acs.jproteome.0c01013 | AlphaPeptDeep (verify license) |

---

## Tier G — Patents, commercial, licenses (verify aggressively)

| # | Source | Verify |
|---|--------|--------|
| G1 | https://patents.google.com/patent/US8639447B2/en | GF method — expiry **2030-07-25**, UCSD |
| G2 | https://doi.org/10.1074/mcp.T600050-MCP200 | Paragon — commercial probability engine |
| G3 | https://github.com/MSGFPlus/msgfplus/blob/master/LICENSE.txt | MS-GF+ UC non-profit |
| G4 | https://github.com/Nesvilab/MSFragger | Academic-only |
| G5 | https://github.com/percolator/percolator/blob/master/license.txt | Apache-2.0 |
| G6 | https://github.com/UWPR/Comet | Apache-2.0 |
| G7 | https://github.com/lazear/sage | MIT |
| G8 | https://instadeepai.github.io/InstaNovo/license/ | Code Apache / weights NC |
| G9 | https://github.com/Noble-Lab/casanovo | Apache-2.0 |
| G10 | https://github.com/melodi-lab/dripToolkit | OSL-3.0 vs paper "Apache" — reconcile |

---

## Tier H — Recent exact-calibration extensions (2023–2025)

| # | DOI / source | Topic |
|---|--------------|-------|
| H1 | https://doi.org/10.1002/pmic.202300145 | Faster XPV / HR-XPV generalization (Bhimani et al.) |
| H2 | https://doi.org/10.1021/acs.jproteome.3c00224 | Crux toolkit update (Tailor integration) — PMC10284583 |
| H3 | https://doi.org/10.1101/831776v1 | Tailor bioRxiv (preprint details) |
| H4 | https://doi.org/10.1021/pr0706698 | RAId aPS (GF-related — check patents) |
| H5 | https://doi.org/10.1186/1752-0509-4-154 | FDR post-processing survey |

---

## Verification checklist (deliver back to parent agent)

```markdown
## Patent verification log
- [ ] US8639447: claims read; list score-enumeration claims verbatim
- [ ] Paragon US patents (search assignee AB Sciex)
- [ ] RAId / aPS patent status
- [ ] Any Tailor patent? (expect none)

## License matrix (code | weights | commercial OK?)
- [ ] Each tool in Tier G

## Low-res evidence table
| Method | Ion-trap/CID 0.5Da tested? | Citation figure/table |
|--------|---------------------------|------------------------|

## Acquire list
- [ ] PMC2689316 (MS-GF 2010)
- [ ] PMC6342018 (res-ev)
- [ ] PMC4185971 (DRIP)
- [ ] PMC2738854 (Frank ranks)
```

---

## Original hypotheses for research agent to test

These are **engine-team ideas**, not literature facts — your job is to find supporting or contradicting evidence.

| ID | Hypothesis | Search for |
|----|------------|------------|
| H-A | Tailor on rank-LLR beats Tailor on hyperscore at 0.5 Da when N_candidates ≥ 100 | Any rank-calibration paper; simulate from A1 |
| H-B | `ChanceMatchSurprise` (local ρ·Δ null) approximates OMSSA Poisson at wide tolerance | Compare Geer 2004 λ to density-based surprise |
| H-C | Integer rounding in `score_psm` costs ≥3% PSMs low-res vs float path | Frank 2009 rank discretization; engine ablation |
| H-D | Percolator cannot recover spectrum-heterogeneity without calibrated input feature | Keich 2014; Percolator + exact p synergy in A2 |
| H-E | GF patent claims may not cover **renewal/saddlepoint on site-visit CGF** (different object) | Claim text US8639447 vs RS³ design doc |
| H-F | Entrapment paired FDP penalizes miscalibrated scores more than TDC | Wen 2025 supplementary simulations |

---

## Local PDFs already collected

See [`internal-docs/papers/REFERENCES.md`](internal-docs/papers/REFERENCES.md):
- Tailor, OMSSA, Wen-Noble entrapment, MSFragger-DDA+, Comet fast-XCorr, Quandenser, peptide-centric, etc.

**Acquire-list (priority):** PMC2689316, PMC6342018, PMC4185971, PMC2738854, PMC2206012, PMC10374903, PMC6602496.

---

## Output format (return to user)

1. Completed verification checklist  
2. Annotated bibliography (55+ entries, 1-line low-res verdict each)  
3. "Surprises" — anything that contradicts parent brief  
4. 3–5 highest-confidence **new** references not in this list
