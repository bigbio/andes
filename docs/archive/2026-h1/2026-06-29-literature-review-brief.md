# Literature Review Brief — DDA Peptide-Identification Scoring

**Date:** 2026-06-29  
**Purpose:** Handoff document for AI agents designing a patent-free Rust DDA search engine ("the engine").  
**Success metric:** PSMs at 1% **true entrapment-FDP** (paired estimator, r=1), not reported target-decoy FDR.  
**Known gap:** ~5% PSM deficit vs MS-GF+ (Java) on low-res UPS1 (≈0.5 Da ion-trap CID) despite leading on high-res Orbitrap/Astral and TMT.

**Local assets:** curated PDFs and acquire-list in [`internal-docs/papers/REFERENCES.md`](internal-docs/papers/REFERENCES.md). Related engine designs: [`2026-06-29-rs3-spectral-significance-design.md`](2026-06-29-rs3-spectral-significance-design.md), [`soft-fragment-matching.md`](soft-fragment-matching.md).

**Agent variants (use these for handoff):**
- **Web-research agent** → [`2026-06-29-literature-review-web-research-agent.md`](2026-06-29-literature-review-web-research-agent.md) (≥55 DOIs, verification checklist)
- **Rust/code agent** → [`2026-06-29-literature-review-rust-agent.md`](2026-06-29-literature-review-rust-agent.md) (file:line hooks, implementation phases)

---

## Agent handoff notes

| Agent type | How to use this doc |
|---|---|
| **Deep-research web agent** | Use the **web-research variant**; verify patents/licenses; return annotated bibliography. |
| **Code-aware agent (Rust/engine)** | Use the **rust-agent variant**; implement Phase 1–2 calibration before novel DP. |

**Established fact vs inference:** Sections label claims `[paper]` when directly from a cited source, `[workspace]` for internal benchmark/design docs, `[inference]` for reasoned synthesis.

---

## Executive summary — top patent-free, low-res-relevant methods (ranked)

| Rank | Method | Why it matters for the engine | Role |
|---:|---|---|---|
| 1 | **Tailor percentile calibration** | Cheap per-spectrum normalization from empirical candidate-score tail; competitive with E-values on low-res in Crux benchmarks [paper] | Additive feature (`score / Q99_null`) |
| 2 | **Empirical per-spectrum null via decoy/Monte Carlo** | Exact tail without GF; aligns with RS³ campaign reformulation [workspace] | Replacement or Percolator feature |
| 3 | **Howbert–Noble exact DP p-values** | Gold-standard exact calibration for additive binned scores (XCorr); proves DP route works on ion-trap data [paper] | Pattern to adapt to own RawScore if decomposable |
| 4 | **OMSSA Poisson E-value** | Classic low-res significance from matched-ion count; public domain [paper] | Additive (`-log10 E`) or sanity baseline |
| 5 | **Andromeda-style local binomial peak-depth** | Position-local intensity weighting; strong on ion-trap; MIT [paper] | Partial ideas (local q, not global top-N) |
| 6 | **Own rank-LLR + soft matching** | Frank/Pevzner rank prior art; engine already implements; lacks per-spectrum calibration [paper][workspace] | Core score — keep, calibrate on top |
| 7 | **X!Tandem hyperscore + log-tail regression** | Simple intensity×factorial score + exponential tail; algorithm patent-free; implementations vary [paper] | Reference baseline, not primary |
| 8 | **Percolator on calibrated features** | Complementary to per-spectrum calibration [paper]; engine uses PSM-level only | Post-score — keep |

**Explicitly off-limits / low-value for Apache engine:**
- MS-GF+ **generating-function spectral probability** — US 8,639,447 B2, active to **2030-07-25** [patent]
- Paragon probability framework — commercial + patented (MCP 2007) [paper]
- MSFragger — academic-only license, not Apache-2.0 [license]
- InstaNovo / many Prosit weight bundles — NC model weights [license]

---

## 1. Primary PSM scoring functions

| Method | Core math (faithful summary) | Low-res (≈0.5 Da CID) | Additive vs replacement | Patent / license | Risks / caveats |
|---|---|---|---|---|---|
| **SEQUEST XCorr** | Binned observed spectrum **u**, theoretical **v** (length N). Fast XCorr: preprocess u → **u′** subtracting mean correlation over τ∈[−75,+75], τ≠0; **XCorr = ⟨v, u′⟩** [paper: Eng 1994, Eng 2008]. Original: R₀ − (1/151)Σ_{τ≠0} R_τ with R_τ = Σᵢ x[i]·y[i+τ]. | **Yes** — designed for low-energy CID ion-trap; bin width ~1 Da [paper] | Replacement (primary rank) | Algorithm: no known patent. Comet: **Apache-2.0** [license] | Background subtraction essential; wide bins lose isotope info; Sp prefilter legacy only |
| **Comet** | Same XCorr + optional **Sp** preliminary score (matched intensity sum + continuity); top-N by Sp then full XCorr [license docs] | **Yes** | Replacement (XCorr) | **Apache-2.0** [license] | Sp threading nondeterminism; FI index (2024+) for speed not low-res-specific |
| **X!Tandem Hyperscore** | **HS = (Σᵢ Iᵢ) · N_b! · N_y!** (matched b/y counts × summed matched intensities) [paper: Craig & Beavis 2003] | **Yes** — historically ion-trap workhorse | Replacement | GPL-ish ecosystem (GPM); **algorithm patent-free** [inference] | Factorial dominates count; needs tail calibration |
| **Sage Hyperscore** | Same hyperscore family as X!Tandem [paper: Lazear 2024] | **Yes** (depends on params) | Replacement | **MIT** [license] | Less mature PTM/open-mod than MSFragger |
| **MSFragger Hyperscore** | **log(N_b! N_y! Σ I_b Σ I_y)** equivalent to X!Tandem form [paper: Kong 2017] | **Yes** at wide tolerance | Replacement | **Academic/non-commercial only**; commercial via Fragmatics [license] | Not Apache-compatible |
| **Andromeda** | In each 100 Th window, take top **q** peaks; count **k** theoretical matches among them; score ∝ **−log P_binom(k \| n, p)** with **p = q / (peaks in window)**; optimize q per spectrum [paper: Cox 2011] | **Yes** — local windows compensate global low-m/z density | Replacement | **MIT** (MaxQuant bundle) [inference] | Uniform p within window approximate; separate from MaxQuant FDR |
| **MS-GF+ rank score** | Spectrum peaks ranked; peptide path in DAG; dot-product style **Score(P,S)** over matched ranks + intensity terms; then **GF tail** for E-value [paper: Kim 2014] | **Yes** — strong on ion-trap in benchmarks [paper] | Replacement + SpecEValue | **US 8,639,447 B2** (UCSD) to 2030; code **UC non-profit license** [patent][license] | **Off-limits** for engine GF path |
| **OMSSA** | Matched product ions **y** ~ Poisson(λ); λ from random-match model over library; **E = N · P(Y≥y)** for N peptides [paper: Geer 2004] | **Yes** — Poisson motivated on ion-trap yeast data [paper] | Replacement (E-value rank) | **Public domain** (NCBI) [license] | Count-based; weak on high-res isotope structure |
| **MyriMatch** | Multivariate hypergeometric on fragment match counts [paper: Tabb 2007] | **Uncertain** — less tested post high-res era | Replacement | **Apache-2.0** (Crux/MyriMatch) [inference] | Hypergeometric assumptions brittle with many mods |
| **Frank rank-based** | Predict fragment **intensity ranks** from sequence (RankBoost); score = quality of rank agreement [paper: Frank 2009] | **Yes** — rank robust to absolute intensity scale [paper] | Additive features or replacement | **Prior art patent-free** [inference] | Needs training data; engine uses own tables |
| **DRIP (DBN alignment)** | Align theoretical/observed peaks with insertions/deletions; Gaussian emission per match; Viterbi log-likelihood [paper: Halloran 2016] | **Yes** — explicit low-res mode + `dripTrain` [paper] | Replacement | Paper says Apache; GitHub **OSL-3.0** — **FTO check** [license discrepancy] | Slow; GMTK dependency |
| **Engine (andes) RawScore** | Σ_s g(m_s) with soft-matched rank LLR per cleavage site [workspace] | Strong high-res/TMT; **weak calibration low-res** [workspace] | Replacement (rank/strong/auto) | Own implementation | Needs per-spectrum calibration layer |

---

## 2. Per-spectrum significance / calibration (central theme)

Goal: **P(score ≥ T | spectrum S, null peptide of mass M)** comparable across spectra.

| Method | Math | Low-res | Additive vs replacement | Patent / license | Risks |
|---|---|---|---|---|---|
| **MS-GF+ generating function** | Build spectrum DAG G; peptide boolean string P; DP over score distributions **F_S(t)** = #{peptides mass M with score ≥ t}; **SpecProb = tail mass** [paper: Kim 2008, 2014] | **Yes** — reference standard low-res [paper] | Replacement (E-value) | **US 8,639,447** to 2030 [patent] | **Do not implement** |
| **Howbert–Noble exact XCorr p-value** | DP enumerates **full score distribution** over all peptides within precursor tolerance for fixed binned evidence vector; **p = P(XCorr ≥ x_obs)** [paper: MCP 2014] | **Yes** — ion-trap sets in paper [paper] | Replacement or feature | **No patent found**; Crux integration [paper] | Requires score = sum of independent bin contributions; expensive O(peptides×bins) |
| **Res-ev + combined p-value** | Res-ev: score from **pairs** of peaks (uses high-res pairs); combine with XCorr p-value [paper: Kertész-Farkas 2018] | **No** for res-ev core; combined needs `--mz-bin-width=1.0005079` for exact p [docs] | Replacement | Crux **Apache-2.0** [inference] | High-res first; low-res use XCorr p-value arm only |
| **X!Tandem / Comet E-value tail fit** | Histogram hyperscore or XCorr across candidates; fit **log survival** linearly; extrapolate to top hit: **e(x) = n · s(x)** [paper: Fenyö 2003; Tailor paper] | **Yes** — default in many ion-trap pipelines [paper] | Replacement | Algorithm patent-free [inference] | Assumes exponential tail; miscalibrated tails common [paper: Tailor] |
| **Tailor** | During search, collect candidate scores per spectrum; **Q99 = 99th percentile of null-like candidates**; **TailorScore = raw / Q99** [paper: Sulimov 2020] | **Yes** — within ~exact p on low-res, 20–150× faster [paper] | Additive | **No patent found** [inference] | Needs sufficient candidates per spectrum; heuristic not exact |
| **OMSSA Poisson / binomial p** | Analytic **P(Y≥y)** from Poisson mean estimated per spectrum [paper: Geer 2004] | **Yes** | Replacement | Public domain | Independence assumptions |
| **Andromeda binomial** | Per-spectrum implicit via varying q and local windows [paper: Cox 2011] | **Yes** | Built into score | MIT | Not full spectrum DP |
| **Decoy / Monte Carlo per-spectrum null** | Sample random peptides of mass M (or shuffle); empirical **p̂ = #{score ≥ T}/N** [paper: general MC; workspace RS³ Gate 0] | **Yes** | Feature or replacement | Patent-free [inference] | Cost; must match null peptide distribution |
| **Saddlepoint / Lugannani–Rice tail** | CGF **K(θ)=Σ log(1−ρ+ρe^{θg})**; solve **K′(θ̂)=T**; tail p from **w, ν** [workspace RS³; textbook stats] | **Uncertain** — RS³ campaign found independence approx insufficient alone [workspace] | Feature (`Rs3NegLog10P`) | **Distinct from US 8639447** if on renewal/site-visit object not score-GF [workspace] | Must calibrate **actual emitted integer score**, not surrogate |
| **Klammer–Noble statistical XCorr calibration** | Parametric mixture model for XCorr null [paper: JPR 2009] | Moderate | Additive | No patent found | Superseded by exact DP for XCorr |
| **Percolator** | Semi-supervised SVM on PSM features; q-values from target-decoy within feature space [paper: Käll 2007] | **Yes** but does not fix spectrum-heterogeneity alone [paper] | Post-processor | **Apache-2.0** [license] | Not entrapment-FDP; needs good inputs |

### Patent-free routes to MS-GF+-like calibration (decision matrix)

| Route | Exact? | Cost | Low-res evidence | Engine fit |
|---|---|---|---|---|
| Tailor on RawScore | No | Very low | Strong [paper] | **Ship first** |
| MC/decoy empirical null | Yes (MC error) | Medium–high | Strong [workspace] | RS³ Gate 0 pattern |
| DP on decomposable score | Yes | High if score tractable | Proven for XCorr [paper] | Hard for LLR+loss+cleavage |
| Saddlepoint on renewal null | Approx | Low per candidate | Unproven [workspace] | Research — decoy-calibrated |
| Hyperscore tail regression | No | Low | Moderate [paper] | Baseline only |

---

## 3. Candidate generation / search-space expansion

| Method | Idea | Low-res | Patent / license | Risks |
|---|---|---|---|---|
| **MSFragger fragment index** | Index theoretical fragments; lookup observed m/z → candidate peptides [paper: Kong 2017] | Yes at 0.5 Da | Academic license only | License blocker |
| **Sage / open search** | Same index idea; ultra-fast Rust [paper: Lazear 2024] | Yes | MIT | Less PTM depth |
| **Comet-FI (2024+)** | Fragment-ion index prefilter to XCorr [paper: PMC13232765] | Yes | Apache-2.0 | Beta; not fastest |
| **Semi-/non-specific digestion** | External sort / bucket by mass; disk-backed peptide lists [engineering] | Yes — critical low-res PTM | Patent-free | Memory/IO bound |
| **Open / mass-offset search** | MSFragger LOS: all peptides in Δm window [paper] | Yes | Academic license | FDR entrapment validation needed |
| **Sequence tags (PepNovo, InsPecT)** | Tag filter → full search [paper: Frank 2005] | Yes on ion-trap | PepNovo: check distribution | Two-pass latency |
| **Chimeric / DDA+** | Full isolation-window search + XIC + shared-fragment removal [paper: Nat Commun 2025] | High-res focused; concept applies | EntrapBench code Apache | Co-fragmentation model |
| **Multi-enzyme / multi-pass** | Cascaded search Crux [paper: PMC] | Yes | Crux Apache | Peptide-anchored FDR bias |

---

## 4. Spectrum refinement feeding scoring

| Method | Idea | Low-res | Patent / license | Notes |
|---|---|---|---|---|
| **Dynamic noise level (DNL)** | Per-spectrum noise floor from peak spacing / intensity [literature: various; verify primary] | **Yes** — critical when S/N low | FTO check needed | Many heuristics unpublished |
| **Windowed top-N (Andromeda q)** | Local top-q per 100 Th [paper: Cox 2011] | **Yes** | MIT | Prefer over global top-50 |
| **Complementary b/y validation** | Require both series / longest b+y [engineering] | **Yes** | Free | Engine has `longest_b/y` features |
| **Deisotoping / charge deconvolution** | Remove C13 spacing, assign z from isotope spacing [paper: various] | **Poor** at 0.5 Da — isotopes unresolved [inference] | — | Avoid aggressive deiso low-res |
| **√I or rank transform** | Stabilize intensity scale [paper: Frank rank work] | **Yes** — rank preferred [paper] | Free | Engine uses rank LLR |
| **Predicted-spectrum cleaning** | Remove peaks not in Prosit/MS2PIP prediction [paper: MS²Rescore] | **Uncertain** — predictors trained mainly HCD high-res | Model licenses vary | Own CID ion-trap models better |
| **Precursor removal / neutral losses** | Filter precursors, NH3/H2O losses in scoring [Comet/XCorr docs] | **Yes** | Free | Already partial in engines |

---

## 5. Learned / cross-family methods

| Method | Idea | Low-res | Code license | Weights license | Engine role |
|---|---|---|---|---|---|
| **Prosit / AlphaPeptDeep** | DL fragment intensity prediction [paper] | CID ion-trap models exist but uneven | Apache (dlomix) | **Check Figshare/Koina per model** | Retrain on own data |
| **MS2PIP** | XGBoost intensity prediction [paper] | Moderate | **Apache-2.0** [license] | Retrainable | Good own-model candidate |
| **pDeep / pDeep2** | RNN intensities [paper] | Mainly HCD | Check repo | Often bundled weights | Retrain only |
| **Spectral angle / SpectraST** | Library match cosine on normalized intensities [paper: Lam 2007] | Weak if library high-res | Various | Library encumbrance | Second pass only |
| **Spec2Vec** | Embedding similarity [paper] | Uncertain low-res | Check repo | Model file license | Low priority |
| **Casanovo** | Transformer de novo [paper: Nat Commun 2024] | Promising immunopeptidomics | **Apache-2.0** code+weights [paper] | Apache | Tag prior / reranker |
| **InstaNovo / ContraNovo** | De novo seq2seq / contrastive [paper] | SOTA de novo | Apache code | **CC BY-NC-SA weights** [license] | Not Apache-clean for weights |
| **MS²Rescore / MSBooster** | Add prediction-based features → Percolator [paper] | Gains mainly high-res | Mixed | Mixed | Feature generator if own models |
| **Percolator / Mokapot** | Semi-supervised FDR rerank [paper] | Yes | **Apache-2.0** [license] | N/A | **In use** — PSM only |

---

## 6. FDR / validation methodology

| Topic | Math / rule | Relevance to engine |
|---|---|---|
| **Target-decoy (TDC)** | FDR ≈ N_decoy / N_target at score cutoff [paper: Käll 2007] | Standard but **can be optimistic** vs entrapment [paper: Wen 2025] |
| **Entrapment FDP** | Inject known false proteins/peptides; estimate FDP from discoveries [paper: Wen 2025] | **Primary gate** for engine |
| **Paired estimator (r=1)** | Uses target–entrapment pairs; tighter upper bound than combined [paper: Wen 2025] | Required methodology |
| **Combined estimator** | Conservative; often overestimates FDP [paper: Wen 2025] | Report but don't optimize to |
| **Group / subset FDR** | FDR within protein family, enzyme, etc. [paper: various] | Secondary |
| **Pitfalls** | Decoy asymmetry (different AA composition); peptide-anchored second pass inflates confidence [paper: Wen 2025; workspace] | Avoid peptide-level second-pass FDR |

---

## Patent & licensing appendix

### Patents (verified or flagged)

| Patent | Assignee | Title / subject | Status | Expiry | Engine impact |
|---|---|---|---|---|---|
| **US 8,639,447 B2** | Regents of UC (Kim, Gupta, Pevzner) | Peptide ID via spectrum **generating function** / peptide reconstructions by score | **Active** [patent: Google Patents] | **2030-07-25** | **GF / SpecEValue off-limits** |
| US 2010/0179766 A1 | Same family | Application publication of 8639447 | Granted as 8639447 | — | Same |
| **Paragon** | AB Sciex | Feature probability / sequence temperature | Active commercial | — | Do not implement |
| **SEQUEST** | Thermo / legacy | Early commercial tool | Expired / encumbered history | — | XCorr **algorithm** still used in Comet |

*Other GF-related publications (spectral dictionaries, spectral networks) are scientific prior art but **implementation of score-enumeration DP** overlaps 8639447 claims — **FTO review required** before any DP-over-scores.*

### Tool licenses (code vs weights)

| Tool | Code | Weights / commercial | Apache-2.0 OK? |
|---|---|---|---|
| Comet | Apache-2.0 | — | **Yes** |
| Crux/Tide | Apache-2.0 | — | **Yes** |
| Sage | MIT | — | **Yes** |
| Percolator, Mokapot | Apache-2.0 | — | **Yes** |
| OMSSA | Public domain | — | **Yes** |
| Andromeda/MaxQuant | MIT (component) | — | **Yes** (verify MaxQuant bundle) |
| MS-GF+ | UC non-profit | Commercial via UC TTO | **No** (license + patent) |
| MSFragger | Academic license | Commercial Fragmatics | **No** |
| Casanovo | Apache-2.0 | Apache-2.0 (per Nat Commun 2024) | **Yes** |
| InstaNovo | Apache-2.0 | **CC BY-NC-SA 4.0** | **No** for bundled weights |
| MS2PIP | Apache-2.0 | Retrainable | **Yes** |
| Prosit/dlomix | Apache-2.0 | **Verify per checkpoint** | Partial |
| DRIP | OSL-3.0 (GitHub) vs Apache (paper) | Retrainable | **Verify** |

---

## Gaps / open problems

1. **Patent-free *exact* per-spectrum calibration for non–dot-product scores** (rank LLR + cleavage + neutral losses). XCorr DP and GF do not generalize cleanly [inference].
2. **Low-res fragment intensity prediction** — most DL models trained HCD/Orbitrap; ion-trap CID needs own training corpus (ProteomeTools CID trap exists) [inference].
3. **Entrapment-FDP at PSM level with Percolator only** — limited literature; peptide-level entrapment more common [workspace].
4. **Chimeric spectra at 0.5 Da** — DDA+ validated mainly high-res Astral [paper]; low-res chimeric handling thin.
5. **Saddlepoint / large-deviation tails** for dependent cleavage sites — renewal independence fails; decoy calibration may be required [workspace].
6. **Isotope-aware scoring at low-res** — generally harmful to deisotope; better soft mass matching [inference].

---

## Original ideas (engine-specific, not literature-established)

These extend prior art for the **andes** architecture; validate only via entrapment-FDP experiments.

| Idea | Rationale | Risk |
|------|-----------|------|
| **Empirical null on existing candidate pool** | Tailor already collects per-spectrum score histogram; `p = #{s≥T}/N` is patent-free SpecProb substitute | N too small on sparse spectra |
| **ChanceMatchSurprise as calibrated feature** | Already implements local ρ·Δ null (`match_engine.rs` ~1759); OMSSA-like without Poisson assumption | Correlates with Tailor — Percolator may down-weight |
| **Score-landscape g(m) cache + fast decoy rescore** | `g(m)` is spectrum-only; decoy peptides score in O(len) using table — cheap MC null | Cleavage credit not in g(m) |
| **DeltaTailor / NegLog10EmpiricalP PIN columns** | Fills hole left by removed `lnSpecEValue`; log-scale for Percolator | Feature dilution if redundant |
| **Low-res-only float ranking** | `rank_score_float` exists; integer round hurts short peptides at low-res [Frank rank prior] | High-res regression if mis-gated |
| **Andromeda q on rank landscape** | Local top-q before summing g(m) — handles low-m/z density without global top-50 | Extra hyperparameter q |
| **Calibration stack** | `log(RawScore)-log(Q99)+α·surprise-β·entropy` — encodes spectrum density heterogeneity | Overfit Percolator on small benchmarks |
| **Peak-rank shuffle null** | Tests whether ID is rank-pattern vs mass coincidence; no peptide DB | Compute cost; novel — no direct citation |

Full implementation notes: [`2026-06-29-literature-review-rust-agent.md`](2026-06-29-literature-review-rust-agent.md).

---

## Recommended implementation sequence (for code-aware agent)

1. **TailorScore** on RawScore + `DeltaTailor` (cheap, proven low-res) [paper]
2. **Empirical null** from in-search decoys per spectrum (RS³ Gate 0) [workspace]
3. Percolator features: `TailorScore`, `NegLog10EmpiricalP`, existing `longest_b/y`, `ChanceMatchSurprise`
4. A/B at **1% paired entrapment-FDP** on UPS1 low-res + regression Astral/TMT
5. Only then: saddlepoint or DP approximations if empirical null cost too high

---

## Bibliography (58 DOIs/URLs — full tiered list in web-research variant)

### Scoring (1–15)

1. https://doi.org/10.1016/1044-0305(94)80016-2 — SEQUEST XCorr  
2. https://doi.org/10.1021/pr800420s — Fast XCorr  
3. https://doi.org/10.1093/bioinformatics/bth023 — X!Tandem  
4. https://doi.org/10.1021/ac025676e — Hyperscore E-value  
5. https://doi.org/10.1021/pr101065j — Andromeda  
6. https://doi.org/10.1021/pr0499491 — OMSSA  
7. https://doi.org/10.1021/pr8001244 — MS-GF spectral probability  
8. https://doi.org/10.1038/ncomms6277 — MS-GF+  
9. https://doi.org/10.1038/nmeth.4256 — MSFragger  
10. https://doi.org/10.1021/acs.jproteome.3c00486 — Sage  
11. https://doi.org/10.1021/pr101196n — Tide  
12. https://doi.org/10.1186/1471-2105-8-327 — MyriMatch  
13. https://doi.org/10.1021/pr8007374 — Frank rank prediction  
14. https://doi.org/10.1021/pr8006788 — Frank rank PSM score  
15. https://doi.org/10.1021/acs.jproteome.6b00290 — DRIP  

### Calibration (16–25)

16. https://doi.org/10.1074/mcp.O113.036327 — Exact XCorr p-value  
17. https://doi.org/10.1021/pr8011107 — Statistical XCorr calibration  
18. https://doi.org/10.1021/acs.jproteome.9b00736 — Tailor  
19. https://doi.org/10.1021/acs.jproteome.8b00206 — Res-ev / combined p  
20. https://doi.org/10.1021/pr5010983 — Keich & Noble calibrated scores / MC null  
21. https://doi.org/10.1002/pmic.202300145 — Faster XPV / HR-XPV  
22. https://doi.org/10.1021/pr0706698 — RAId aPS  
23. https://doi.org/10.1093/bioinformatics/btn189 — Klammer DBN fragmentation  
24. https://doi.org/10.1101/831776v1 — Tailor bioRxiv  
25. https://doi.org/10.1074/mcp.M700022-MCP200 — Morpheus  

### FDR / validation (26–35)

26. https://doi.org/10.1021/pr700600n — Target-decoy  
27. https://doi.org/10.1038/nmeth1113 — Percolator  
28. https://doi.org/10.1074/mcp.T900012-MCP200 — PeptideProphet  
29. https://doi.org/10.1021/ac0341261 — ProteinProphet  
30. https://doi.org/10.1074/mcp.M900317-MCP200 — MAYU  
31. https://doi.org/10.1002/pmic.201500431 — Protein-level FDR  
32. https://doi.org/10.1016/j.jprot.2010.08.009 — Nesvizhskii survey  
33. https://doi.org/10.1002/rcm.4417 — Empirical FDR  
34. https://doi.org/10.1038/s41592-025-02719-x — Entrapment FDP  
35. https://doi.org/10.1186/1752-0509-4-154 — FDR survey  

### Candidate gen / de novo / chimeric (36–45)

36. https://doi.org/10.1038/nmeth.1889 — InsPecT  
37. https://doi.org/10.1021/pr0500111 — Peptide tags  
38. https://doi.org/10.1021/ac048788h — PepNovo  
39. https://doi.org/10.1038/s41467-024-49731-x — Casanovo  
40. https://doi.org/10.1038/s41587-024-01382-9 — InstaNovo  
41. https://doi.org/10.1074/mcp.M110.003731 — Cascaded search  
42. https://doi.org/10.1038/s41467-025-58728-z — MSFragger-DDA+  
43. https://pmc.ncbi.nlm.nih.gov/articles/PMC13232765 — Comet FI  
44. https://doi.org/10.1089/cmb.2014.0165 — GF spectral networks  
45. https://doi.org/10.1021/pr300631t — De-Noise (ion-trap)  

### Libraries / quant-first / processing (46–52)

46. https://doi.org/10.1002/pmic.200600625 — SpectraST  
47. https://doi.org/10.1038/nmeth.1240 — SpectraST consensus  
48. https://doi.org/10.1021/pr900473s — SpectraST decoy libs  
49. https://doi.org/10.1371/journal.pcbi.1008724 — Spec2Vec  
50. https://doi.org/10.1021/pr0700693 — Peptide-centric  
51. https://doi.org/10.1038/s41467-020-18138-8 — Quandenser  
52. https://doi.org/10.1021/pr401026y — Empirical multidimensional scoring  

### Learned rescoring (53–58)

53. https://doi.org/10.1038/s41592-019-0426-7 — Prosit  
54. https://doi.org/10.1093/nar/gkz299 — MS2PIP  
55. https://doi.org/10.1038/nbt.4313 — pDeep  
56. https://doi.org/10.1038/s41467-023-40129-9 — MSBooster  
57. https://doi.org/10.1016/j.mcpro.2022.100266 — MS²Rescore  
58. https://doi.org/10.1021/acs.jproteome.3c00785 — MS²Rescore 3.0  

### Patents & licenses (non-DOI)

- US 8639447 B2 — https://patents.google.com/patent/US8639447B2/en  
- Paragon — https://doi.org/10.1074/mcp.T600050-MCP200  
- Comet — https://uwpr.github.io/Comet/  
- MS-GF+ — https://github.com/MSGFPlus/msgfplus/blob/master/LICENSE.txt  
- MSFragger — https://github.com/Nesvilab/MSFragger  
- Percolator — https://github.com/percolator/percolator/blob/master/license.txt  
- Crux toolkit — https://doi.org/10.1021/acs.jproteome.3c00224 (PMC10284583)  

**Extended tiers A–H (65+ entries):** [`2026-06-29-literature-review-web-research-agent.md`](2026-06-29-literature-review-web-research-agent.md)

### Local PDFs (workspace)

See [`internal-docs/papers/REFERENCES.md`](internal-docs/papers/REFERENCES.md).

---

*Document prepared for agent handoff. Patent and license statuses verified via Google Patents and official repositories as of 2026-06-29; not legal advice — formal FTO review recommended before release.*
