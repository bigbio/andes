# Andes competitive search plan: Comet speed, MSFragger glyco depth

**Date:** 2026-09-03  
**Branch reviewed:** `feat/glyco-pin-curation`  
**Scope:** standard DDA speed and PSM yield; N-glycopeptide speed, yield, and error control  
**Decision standard:** correct identifications at matched empirical error, not nominal counts alone

## Executive decision

Andes already produces more standard-search PSMs than the **classic Comet** configuration used in the existing benchmark: +22.3% on Astral and +17.2% on low-resolution TMT. It is not yet faster: the current post-optimization measurements are 555 s versus 217 s on Astral and 143 s versus 80 s on TMT. Those results do **not** establish competitiveness with current Comet-FI; the comparator must be upgraded to Comet 2026.02.2 and run in fragment-index mode.

The updated high-resolution glyco retrieval default is the largest successful speed change so far. Moving retrieval from the rank model's 0.5 Da window to 20 ppm reduced six-fraction mouse wall time from 144 to 21 minutes (6.9x) while correct identifications moved from 3183 to 3198. Plasma moved from 380 +/- 8 to 399 +/- 23 glycoPSMs and approximately 7x faster search. This should be retained after the candidate-recall diagnostic and two implementation defects below are resolved.

The glyco identification gap remains real: the best comparable plasma result is about 399 versus 587 MSFragger-Glyco glycoPSMs, so parity requires approximately **47% more** Andes identifications at matched peptide and glycan error. The failed split election must not ship: on mouse it promoted more nominal identifications but fewer correct identifications and materially higher true error. The next ID work is retrieval calibration plus separate peptide/glycan evidence and FDR, not another hand-weighted global election.

## What the current code actually does

### Standard DDA path

The standard hot loop first scores every precursor-matched candidate with the MS-GF-style node score, retains a small finalist set (normally 25 under strong scoring), then computes expensive learned features for those finalists. The two-stage edge gate is effective and skips most edge-score calls after the queue fills.

The measured profile says approximately 89% of standard-search time is tree inference and only about 2% is candidate work. That agrees with the code:

- `Tree::eval` performs branch-dependent pointer chasing through six parallel arrays for every tree and row.
- Fragment-intensity inference is batched tree-outer/row-inner, which removed one duplicate ensemble walk, but still allocates ion, feature, row, raw-output, and final-output vectors for every finalist.
- The spectrum peak model is still evaluated one peak at a time with `predict_logit`, despite the batch evaluator already existing.
- Final feature construction calls `score_psm` again to reconstruct `RankScoreFloat`, after the rounded score was already computed in the candidate loop.

Conclusion: adding a standard peptide fragment index is not the first speed move. It attacks the measured 2% while leaving the 89% intact.

**Native profile of the current binary (2026-09-03, Codon, `perf record -g --call-graph=dwarf`, 5,760-spectrum fixture vs E. coli, 8 threads, 4.8 s total, 3.4 s in the search):** `GbdtPeakModel::predict_value_batch` is **77.0% self time**; 82% of wall sits under the finalist feature fill (`fill_post_topn` → `compute_psm_features`), of which 52.6% is the fragment-intensity ensemble and the remainder the rich-ion ensemble; the candidate loop closure is 4.7% self, node scoring 3.4%, nearest-peak lookup 1.5%, allocation under 2% in total. The premise above is confirmed on the current code: standard search is ensemble-evaluation-bound. That is also why every item in A1/A2 measured neutral — none of them is inside the 77%. The only standard-path levers left are the ensemble evaluation itself (compiled/bitvector inference; every shipped tree has ≤ 64 leaves, measured) and the model's size (A3).

### Glyco path

The glyco path has the opposite shape. It is candidate-generation and repeated exact-scoring bound:

- The peptide fragment index is `FxHashMap<bin, Vec<posting>>`.
- Every query allocates hash sets for consumed candidate/ion pairs and per-peak candidate membership, plus a hash map of candidate counts.
- The query votes over the full sequon peptide index before peptide mass is conditioned by precursor-minus-glycan hypotheses.
- The fragment index is an additive rescue route; the exhaustive backbone/glycan lattice is still generated and scored.
- Phase one materializes a new `Vec<usize>` for the peptide mass bucket of every backbone, linearly searches charge-scored spectra, then calls both `score_psm` and `psm_edge_score` for every surviving peptide.
- Each retention axis allocates an index vector and fully sorts it, even though only top K is required.
- Backbone and glycan-Y solvers also use per-spectrum hash aggregation.

This architecture explains why the retrieval-tolerance change had such large leverage and why selector-weight experiments did not.

**Native profile at the new 20 ppm default (2026-09-03, mouse Frac1, 149 s glyco phase, 8 threads):** the fragment index no longer appears in the profile at all (`FragmentIndex::build` 0.26%, query below the noise); 64% of self time is the inlined body of `score_spectrum_glyco` itself, 12% is `hashbrown` rehashing (per-spectrum maps growing without capacity), `score_psm` 3.5%, `ion_match_facts` 3.2%, `nearest_peak_full` 2.5%, `best_frag_intensity` 1.9%, allocator ~2.5%. **Item 5 (CSR postings) is therefore withdrawn as a speed item** — it targets a cost the ppm window already removed. The inlined-callee profile attributes that 64% as ~55% hash-table work (membership 18%, entry lookup 21%, insert 14%) plus 11% rehashing. **Converting those maps to FxHash with preallocation measured NEUTRAL** (mouse Frac1, byte-identical PIN, 179.3 s vs 178.3 s), so it was kept only for a determinism fix it happened to carry - the election's split map is summed by `log_mean_exp`, and float addition is order-dependent, so std's per-process random iteration order varied that score run to run. **Methodological caveat for the next profile:** `perf --children` credits an inlined frame with its caller's entire subtree, so an inlined-frame share is an upper bound on cost, not a recoverable one; only self time and an A/B measure a lever. The glyco path's remaining cost is spread across the driver's per-candidate work (`score_psm` 4.3%, `ion_match_facts` 4.0%, `cz_hyperscore_psm` 4.6%, `psm_edge_score` 3.5% inclusive) with no single dominant callee - i.e. no cheap win left there; the next real glyco lever is architectural (items 6-7: mass-conditioned retrieval and the two-route cascade).

## Review of the updated high-resolution retrieval default

The change is directionally correct and its low-resolution behavior is preserved, but four gates remain before calling it production-ready.

1. **Resolution is inferred from the selected model, not necessarily the acquisition.** The default checks `param.data_type.instrument.is_high_resolution()`. Model routing intentionally maps some high-resolution CID and ETD cases to `LowRes`; those acquisitions can therefore retain the 0.5 Da retrieval path. Resolve retrieval resolution from detected MS2 metadata, then treat model selection as an independent scoring decision.
2. **The ppm bin-neighborhood guarantee ends at m/z 5000.** The bin width is sized from `PPM_BIN_MAX_MZ = 5000`, while queries inspect only the adjacent three bins. A valid within-ppm match above that m/z can be more than one bin away. Either derive the maximum indexed fragment m/z during construction or query `ceil(window/bin_width)` neighboring bins.
3. **There is no explicit fixed-Da compatibility control.** The ppm CLI accepts only positive values, so after the automatic high-resolution default there is no clean command-line A/B for the old fixed 0.5 Da retrieval. Add a mutually exclusive `--glyco-retrieval-tol-da`, or an explicit retrieval mode.
4. **The promised recall accounting is still missing.** Before promotion, record postings visited, candidates above threshold, truth candidate rank before and after the cap, correct winners, and true FDP. Require at least 99% recall of independently correct baseline candidates before final scoring.

The changed high-resolution golden is useful, but stable winning peptides in 120 scans is not a substitute for the production candidate-recall gate. The competition columns changed because the candidate population changed; that is expected and means the change is not byte-identical by design.

## The competitor baseline must be reset

### Comet

The existing Andes benchmark used classic Comet. The current target must be **Comet 2026.02.2 in FI mode**, because recent Comet work materially changes the speed baseline:

- fragment-ion indexing performs a cheap candidate screen before exact XCorr;
- the published Comet-FI work reports 400-800x fewer exact XCorr evaluations in typical human searches;
- current releases use per-worker arenas/scratch storage, asynchronous spectrum readahead, contiguous CSR-style structures, cached masses, and reduced hash/allocation traffic;
- 2026.02.2 also fixed an FI modification-scoring defect that changed phosphopeptide PSM yield by roughly 4-4.6%, so older FI results are not a valid identification baseline.

Report both cold and warm index timing. Comet's index-build cost is amortized over multiple files, so a single blended number can make either engine look artificially favorable.

### MSFragger-Glyco

The target should be MSFragger 4.4.1 plus the current FragPipe glyco workflow, not only the search JAR. MSFragger 4.3 specifically accelerated HCD/CID glyco and labile searches. The full identification claim also includes PTM-Shepherd's multiattribute glycan assignment and separate glycan FDR.

Report two comparisons:

1. search-only glycoPSMs, with identical spectra, FASTA, glycan list, enzyme, modifications, tolerances, isotope hypotheses, and PSM-level FDR policy;
2. end-to-end intact glycopeptides after peptide FDR **and** glycan-composition FDR.

Match-between-glycans or match-between-runs must be a separate third result. It can legitimately expand biological coverage, but it is not a spectrum-search PSM and must not be mixed into the base-search count.

## Workstream A: make standard DDA faster than Comet without surrendering PSMs

### A0. Freeze a current benchmark contract

Run Andes and Comet 2026.02.2 classic/FI on Astral and low-resolution TMT, plus one modification-heavy workload. Pin CPU, threads, compiler target, input format, FASTA hash, modifications, tolerances, output scope, and rescoring. Record:

- cold wall time, warm-search wall time, CPU time, peak RSS, and index size;
- candidates visited and exact-score calls per spectrum;
- PSMs and unique peptides at matched picked/competition FDR;
- overlap, Andes-only, Comet-only, and entrapment-supported correctness.

**Gate:** no public “faster than Comet” claim until Andes beats current Comet-FI on at least two named workloads. Classic-Comet parity is only an intermediate milestone.

### A1. Remove avoidable inference overhead, byte-identically

Implement as independent changes, each with a benchmark and golden output check:

1. Batch the per-spectrum peak GBDT with `predict_logit_batch`.
2. Introduce per-worker reusable scratch storage for peak features, fragment features, row references, raw predictions, matched ions, and histograms.
3. Preserve the rounded node score or cleavage credit in the retained PSM so `RankScoreFloat` does not call `score_psm` again.
4. Fuse fragment generation, feature extraction, inference, and conversion into caller-provided buffers.
5. Profile an array-of-compact-nodes tree representation against the current six-vector SoA layout. Keep it only if exact predictions are bit-identical and native profiles improve.

**Gate:** every item is output-identical and independently reduces production wall time. Do not merge a bundle that hides regressions.

### A2. Eliminate repeated peptide-only model work

Fragment-intensity features are a function of peptide, fragment identity, precursor charge, and a currently constant NCE value of 0. That makes their raw predictions reusable across spectra.

Measure three alternatives:

- a bounded per-worker lazy cache keyed by peptide and precursor charge;
- precomputed raw `f32` predictions stored beside candidates for common charges;
- precomputed fragment features with inference still performed in batches.

The right choice depends on hit reuse and memory. Do not compress predictions until a non-bit-identical accuracy study is explicitly approved.

### A3. Only then consider model simplification

The existing 300-to-100-tree option cuts wall time about 41%, not 67%, and showed a small TMT identification loss. Treat tree pruning or distillation as a scoring-model experiment:

- distill to 64/100/150 trees using ranking loss, not raw regression error alone;
- validate PSM and peptide yield on held-out high- and low-resolution data;
- require noninferiority at matched empirical FDP.

### A4. Preserve and extend the PSM advantage

Andes already beats the old classic-Comet PSM count on two datasets. Generalize that claim rather than optimizing one fixture:

- stratify disagreements by charge, length, missed cleavage, modification, precursor interference, and score margin;
- route high-resolution scoring models from acquisition metadata correctly;
- benchmark the existing co-isolation cascade against MSFragger DDA+ and Comet on wide-isolation Astral data;
- train/rescore on held-out datasets and report picked peptide/protein FDR as well as PSM FDR.

The success metric is additional **correct** PSMs. A larger target list at a higher entrapment FDP is a failure.

## Workstream B: exceed MSFragger-Glyco

### B0. Build the oracle diagnostic before another algorithm change

For every truth-supported or high-confidence comparator scan, classify the failure at four boundaries:

1. true peptide absent from generated database candidates;
2. present but missed or capped by retrieval;
3. retained but loses the Andes election;
4. wins search but is rejected by Percolator/glycan FDR.

Report candidate rank distributions (median, p90, p99), not just recall. This tells us whether to work on retrieval, scoring, or calibration. Previous evidence that Andes agrees with MSFragger on roughly 82% of its high-confidence scans means the remaining gap is likely concentrated in weak spectra and final separation.

### B1. Make the fragment index contiguous and allocation-free per spectrum

Replace `FxHashMap<bin, Vec<posting>>` with CSR-style contiguous postings:

- `bin_offsets: Vec<u32>` plus a flat posting array;
- candidate counts in generation-stamped arrays indexed by candidate ID;
- one-peak/one-candidate and one-ion/one-candidate state encoded with generation stamps or sorted local postings rather than hash sets;
- stored document frequency per bin for rarity weighting;
- exact per-ion tolerance validation retained.

**Gate:** identical retrieved candidate/count pairs under differential randomized tests, bounded memory on the largest glycan database, and a production timing win.

### B2. Apply precursor/glycan mass constraints before fragment votes

This is the next high-leverage architectural change. Sort sequon peptides by mass. For each precursor charge/isotope hypothesis, form the union of peptide-mass intervals implied by the allowed glycan masses. Use that union as an eligibility mask while visiting fragment postings.

The current order—vote over the whole sequon database, then test glycan subtraction—does needless work and admits more coincidental competitors. Mass conditioning should improve both speed and the score distribution presented to Percolator.

**Gate:** at candidate budgets 25/50/100/200, retain at least 99% of independently correct baseline peptides and do not increase true FDP.

### B3. Promote a calibrated two-route retrieval cascade

Use complementary evidence routes:

```text
spectrum
  |-- peptide fragments --> supported sequon peptides --> glycan by mass subtraction
  `-- complementary Y --> supported glycans/backbones --> peptide mass slice
                                      |
                                      v
                         bounded label-blind union
                                      |
                                      v
                              exact scoring once
```

The exhaustive lattice becomes a confidence-triggered fallback when both routes are weak. Shadow mode must first record every baseline winner that each route would retain.

Candidate confidence may use only frozen, label-blind observables: independent-ion count, rarity-weighted vote, score margin, core-Y quorum, precursor mass fit, and oxonium compatibility. Target/decoy handling and caps must be symmetric.

### B4. Improve retrieval evidence without repeating the failed election

Evaluate offline, in this order:

1. 20 ppm distinct b/y count;
2. mask confidently assigned oxonium and core-Y peaks from the peptide channel;
3. weight matches by inverse posting frequency;
4. add bounded peak-intensity evidence;
5. add glycopeptide-specific retrieval channels: glycosite-retaining b/y+HexNAc where physically appropriate, and candidate-specific Y0/Y1/Y2 complements.

Do not simply add shared oxonium/Y evidence to every candidate. The failed election demonstrated that shared per-backbone terms can widen separation on wrong winners. Candidate comparisons must use **exclusive explained evidence**: peaks uniquely explained by one peptide/glycan hypothesis versus its nearest competitor.

### B5. Separate peptide identity, glycan identity, and localization confidence

MSFragger's modern workflow closes an important statistical gap after search: PTM-Shepherd assigns glycans using Y ions, oxonium ions, mass error, and isotope evidence, then controls glycan FDR separately.

Andes should produce three calibrated quantities:

- peptide/backbone score and peptide q-value;
- glycan-composition score and glycan q-value;
- site-localization probability or localization q-value where the data support it.

Build glycan target/decoy or entrapment evaluation into the benchmark. Optimize the joint accepted set only after each component is calibrated. The former split election, placeholder glycan weight of 1.0, and shared Y terms are not acceptable substitutes.

### B6. Recover weak repeated glycoforms as a separate post-search layer

The new FragPipe match-between-glycans work shows that learned RT/IM shifts between related glycans can recover low-abundance glycoforms with precursor-level FDR. Andes already has transfer and RT machinery, so a conservative analogous layer is plausible:

- seed only from high-confidence peptide+glycan assignments;
- require isotope-envelope, mass, RT/IM-shift, and optional Y0/Y1 support;
- lock peptide and target/decoy label through transfer;
- estimate transfer FDR independently;
- report transferred glycoforms separately from MS2-derived glycoPSMs.

This can increase glycopeptide coverage, but it must not be used to claim more searched PSMs than MSFragger.

## Ordered 30-day execution plan

| Order | Deliverable | Expected leverage | Promotion gate |
|---|---|---:|---|
| 1 | ~~Fix acquisition-resolution default, >5000 m/z ppm lookup, explicit fixed-Da A/B, resolved-tolerance logging~~ **Done 2026-09-03**: window resolved from detected analyzer metadata (metadata-less input falls back to the `--fragment-tol-*` unit); query walks `ceil(window/bin_width)` bins each side, with a forced multi-bin unit test; `--glyco-retrieval-tol-da` added, mutually exclusive with the ppm flag; one `glyco retrieval window:` line per run. Goldens unchanged. | correctness | all glyco goldens + randomized boundary tests |
| 2 | Produce retrieval oracle and candidate-recall/rank dump on mouse and pooled plasma. **Mouse Frac1 measured 2026-09-03**: true-candidate recall 97.0% at 20-30 ppm, 97.7% at 40, 98.5% at 60 (the four lost sit >60 ppm off, none ever won; pool size unchanged at 473/scan; shipped top-1 111-112 at every width). The 99% line is not reached at any ppm width; the default was promoted on the end-to-end result and the deviation is recorded in the glyco roadmap. Still owed: plasma, rank distributions, before/after the cap. | decision-enabling | >=99% truth-candidate recall for default |
| 3 | Rebenchmark against Comet 2026.02.2 FI and MSFragger 4.4.1/FragPipe | establishes real gap | reproducible cold/warm harness |
| 4 | Batch peak GBDT and add per-worker scratch in standard search. **Measured 2026-09-03** (Codon, 5,760-spectrum standard fixture vs E. coli, 8 threads, 3 interleaved reps, identical output checksum on every run): batching the per-peak model (`cc3743d7`) and recovering the cleavage credit instead of re-walking `score_psm` (`8be751de`) are byte-identical and **speed-neutral** (5.1-5.25 s throughout); the per-worker fragment-prediction cache (A2 lazy variant) was also neutral (5.08-5.48 s vs 5.11-6.10 s) and was reverted. The per-spectrum peak model and the second score walk were not material; the remaining scratch-buffer and compact-node items should be profiled on a production-sized workload before implementation. | high standard speed | byte-identical output, repeated timing win |
| 5 | ~~CSR glyco postings plus generation-stamped counters~~ **Withdrawn 2026-09-03**: at the 20 ppm default the index is below profile noise (see the glyco-path profile note above); the cost moved to the driver body (64%) and map rehashing (12%). | high glyco speed | identical retrieval, bounded RSS |
| 6 | Mass-conditioned fragment voting in shadow mode | very high glyco speed/precision | recall and true-FDP gate |
| 7 | Activate bounded two-route cascade with exhaustive fallback | transformative | no correct-ID loss; large exact-score reduction |
| 8 | Offline exclusive-evidence and separate glycan-score study | high glyco IDs | improves correct IDs at matched peptide+glycan FDP |
| 9 | Fit/freeze joint calibration; only then reconsider election | high glyco IDs | held-out mouse and plasma gains |
| 10 | Evaluate conservative RT/IM glycoform transfer separately | coverage | independent precursor-level FDR |

## Quantitative milestone ladder

### Standard search

- **Milestone 1:** preserve the existing +17-22% PSM lead and reach classic-Comet wall-time parity.
- **Milestone 2:** reach current Comet-FI warm-search parity on two named workloads.
- **Milestone 3:** beat Comet-FI by at least 10% wall time with no lower correct-PSM count and no higher empirical FDP.

### Glyco search

- **Milestone 1:** promote the 20 ppm default after recall/correctness fixes; preserve approximately 399 plasma glycoPSMs and the 6.9-7x speed gain.
- **Milestone 2:** reduce exact candidate scoring by at least 10x via mass-conditioned two-route retrieval while retaining >=99% truth candidates.
- **Milestone 3:** reach at least 500 plasma glycoPSMs at matched peptide and glycan FDP.
- **Milestone 4:** exceed the current 587 comparator result on the identical end-to-end workflow and hold out mouse, with confidence intervals excluding parity.

## Stop conditions

- Do not build a standard-search fragment index unless a new profile shows candidate work has become material.
- Do not hand-write SIMD/bitvector tree traversal before the exact compiled ensemble benchmark shows it beats batching, scratch reuse, and compact nodes.
  **Measured 2026-09-03** (Codon, AMD EPYC 9555, one thread, the shipped `hcd_astral_tryp` fragment-intensity ensemble: 300 trees, 37,614 nodes, 19 features, 36,864 rows sampled between each feature's split-threshold range; every evaluator within 8e-6 of the Rust reference): our `predict_value_batch` 8.7 µs/row at PSM-sized batches (48 rows) and 7.4 µs/row whole-set; treelite/tl2cgen compiled C 9.8 / 8.6 µs/row; LightGBM's predictor 11.9 / 11.4 µs/row. tl2cgen does not beat the walker, but **lleaves (LLVM-compiled trees, `llvmlite<0.44`) does: 5.96 µs/row at 48-row batches and 3.27 µs/row whole-set** — the gap between those two is Python call overhead, so the in-kernel advantage over our walker is ≈2.6×, exact to 8e-6. lleaves compiles each tree into straight-line branch code with thresholds and leaf values as immediates (no threshold/child array loads), which is what a data-driven walker cannot do. Bringing that into andes means either ahead-of-time code generation for the bundled models at build time or a JIT (Cranelift) at model-load time; with the 77% profile share, a 2.6× kernel is worth ≈1.9× end-to-end on the standard path — the largest standard-path lever now measured, and a decision about build complexity rather than about algorithms. The standard-path lever is therefore the QuickScorer-family layout (all shipped trees ≤ 64 leaves) or a smaller model (A3), not compilation. Scaffolding: `scripts/gbdt_bench/`, `crates/scoring/tests/gbdt_bench.rs`.
- Do not ship the split election in its current form.
- Do not activate glyco pruning until shadow-mode truth recall is measured.
- Do not count nominal Percolator gains when entrapment FDP worsens.
- Do not mix MS1 transfer/MBG-like identifications with MS2 search PSMs in competitive claims.

## Evidence ledger and sources

### Repository evidence

- Standard benchmark: [`docs/benchmarks/2026-08-23-andes-vs-comet-refresh.md`](../benchmarks/2026-08-23-andes-vs-comet-refresh.md)
- Glyco conclusions: [`docs/benchmarks/glyco-algorithm-conclusions.md`](../benchmarks/glyco-algorithm-conclusions.md)
- Current glyco speed roadmap and measured 20 ppm result: [`docs/plans/glyco/50-roadmap/2026-09-02-redesign-benchmark-and-speed-plan.md`](glyco/50-roadmap/2026-09-02-redesign-benchmark-and-speed-plan.md)
- Failed election design/results: [`docs/glyco-scoring-redesign-2026-09.md`](../glyco-scoring-redesign-2026-09.md)
- GBDT evaluator: [`crates/scoring/src/gbdt_eval.rs`](../../crates/scoring/src/gbdt_eval.rs)
- Standard hot loop: [`crates/search/src/match_engine.rs`](../../crates/search/src/match_engine.rs)
- Glyco fragment index: [`crates/search/src/glyco_fragment_index.rs`](../../crates/search/src/glyco_fragment_index.rs)
- Glyco search: [`crates/search/src/glyco_search.rs`](../../crates/search/src/glyco_search.rs)

### Primary and official external sources

- Comet current releases and 2026.02 line: [Comet release index](https://uwpr.github.io/Comet/releases/)
- Comet-FI implementation, timing, memory, and amortization notes: [Comet fragment-ion indexing notes](https://uwpr.github.io/Comet/notes/20241001_FI.html)
- Comet-FI paper: [Comet fragment-ion indexing for enhanced peptide sequencing](https://pmc.ncbi.nlm.nih.gov/articles/PMC13232765/)
- MSFragger algorithm: [MSFragger: ultrafast and comprehensive peptide identification](https://pmc.ncbi.nlm.nih.gov/articles/PMC5409104/)
- MSFragger-Glyco algorithm: [Fast and comprehensive N- and O-glycoproteomics analysis with MSFragger-Glyco](https://pmc.ncbi.nlm.nih.gov/articles/PMC7606558/)
- Current MSFragger version history: [MSFragger changelog](https://github.com/Nesvilab/MSFragger/blob/master/CHANGELOG.md)
- Multiattribute glycan assignment and glycan FDR: [Multiattribute glycan identification and FDR control](https://pmc.ncbi.nlm.nih.gov/articles/PMC8933705/)
- Glycan-first complementary-ion indexing: [pGlyco3](https://pmc.ncbi.nlm.nih.gov/articles/PMC8648562/)
- Post-search coverage expansion with explicit precursor-level FDR: [Match-between-glycans](https://pmc.ncbi.nlm.nih.gov/articles/PMC12934668/)

## Confidence statement

High confidence: the two paths require different optimization strategies; current standard search is inference-bound; current glyco search is retrieval/candidate-bound; the 20 ppm retrieval change is a genuine large speed gain; the present split election is harmful on mouse; the old Comet comparator is obsolete.

Medium confidence: CSR postings, generation-stamped counters, mass-conditioned voting, and the two-route cascade are the highest-leverage next glyco changes. They are strongly supported by the code and successful external architectures, but their exact speed and recall effects must be measured in Andes.

Hypothesis requiring experiment: separate exclusive glycan evidence plus joint peptide/glycan calibration can close the remaining approximately 47% plasma gap. The literature supports the mechanism, but the existing Andes experiments do not yet prove the gain.
