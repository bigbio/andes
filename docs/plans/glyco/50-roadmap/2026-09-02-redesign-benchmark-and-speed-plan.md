# Andes glyco — redesign disposition, benchmark repair, and speed plan

**Date:** 2026-09-02
**Branch:** `feat/glyco-pin-curation`
**Status:** execution plan; v2 benchmark COMPLETE (36/36 searches, 20 evaluations) — Stage 1 gate met, Stage 2 resolved to the refutation branch (see 1.1)
**Inputs:** the implemented scoring redesign and two review passes, today's mouse/plasma benchmark results, the entrapment-evaluator fix, and the glyco performance review.

## 1. Decisions already supported by evidence

### 1.1 Split election: refuted in its current form

On the clean mouse benchmark, across all five Percolator seeds, the split election produced:

- **+554 nominal identifications** at 1% reported FDR;
- **−58 correct identifications**;
- approximately **4.5× the true error rate**.

This is a structural failure, not a small calibration miss in the reported count. The election promotes higher-scoring wrong winners, and Percolator then separates those wrong winners more confidently. A larger nominal accepted set is therefore not evidence of improvement.

**v2 confirmation (binary with all five second-review fixes, five seeds, same data):**

| arm | glycoPSMs @1% | correct @1% | true FDP |
|---|---|---|---|
| def0 baseline | 4833 | 3137 | 1.04% |
| defC election (v1) | 5387 | 3079 | 4.64% |
| defC election (v2) | 5200 | 3049 | 3.71% |
| curC election on curated policy (v2) | 4361 | 2916 | 0.97% (cur0: 4409 / 2960 / 1.18%) |

The fixes reduced the damage without changing its direction: v2 still reports 367 more nominal identifications while getting 88 fewer right, at 3.6× the true error. Under the curated policy the election is simply worse. The additive-column arms are neutral in both binaries (defA v2: 4833 / 3137 / 1.16%; defB v2: 4814 / 3130 / 1.02%).

**Disposition:** `--glyco-split-election` is **refuted in its current form** (Stage 2 gate, refutation branch). It stays default-off and experimental. Do not tune its placeholder weight; a second attempt goes through Section 3.

### 1.2 Additive redesign columns: neutral so far

The two additive-column arms are identification-neutral on the valid benchmark. Keep the columns available behind their flags because they are useful experimental observables, but do not promote them to defaults or claim a yield benefit.

### 1.3 Plasma benchmark: not interpretable in the current setup

The plasma arm ranges from 0 to 335 identifications within one arm. Seed 3 reaches the known target-decoy step-function floor in every arm, and the current entrapment database roughly doubles the sequon-bearing search space relative to the database behind the 384.6 reference baseline.

**Disposition:** do not compare these plasma results with 384.6, do not pool them into a redesign conclusion, and do not spend more scoring compute on this setup. Plasma becomes actionable only after the benchmark is made like-for-like.

### 1.4 Review fix: accepted

The pair-on-generation flag now requires an actual paired scan, matching its documentation. Glyco goldens remain byte-identical and all test suites pass. This is an independent correctness fix and is retained regardless of the election result.

## 2. Ordered execution plan

Each stage has a gate. Do not start a dependent stage merely because its code is ready.

### Stage 1 — Finish and read the v2 benchmark

The v2 binary contains today's five review fixes. Complete the existing 36-search matrix and its chained evaluations before launching another redesign experiment.

Record, for every arm and seed:

- nominal PSMs at 1% FDR;
- correct identifications against the independent truth set;
- true FDP;
- target and decoy counts above the threshold;
- winner changes relative to the shipped selector;
- results stratified by paired/unpaired scan and fragmentation regime.

**Gate:** all 36 searches and evaluations complete, with summary files checked for completeness. At 11/36, no final v2 conclusion is permitted.

### Stage 2 — Dispose of the current election

Use mouse as the decision dataset because its five-seed result is readable. Treat plasma as unavailable until Stage 3.

- If v2 again increases nominal IDs while reducing correct IDs or materially increasing true FDP, mark the present split election **refuted** and keep it default-off.
- If v2 reverses the result, investigate which review fix caused the reversal before considering shipment; require replication rather than reading a single favourable seed.
- Keep the additive columns behind their existing flags unless they independently improve a valid endpoint.

**Gate for shipment:** the election must improve correct identifications without worsening true FDP, across the predefined seeds and both usable fragmentation regimes. Nominal 1% identifications alone cannot pass the gate.

### Stage 3 — Make plasma measurable

Rebuild the plasma comparison before using it for any design decision:

1. Obtain or reconstruct the like-for-like yeast entrapment database used by the 384.6 baseline. Note: the cluster file named as a human+yeast entrapment set holds paired-shuffle HUMAN sequences as its `ENTRAP_` entries (equal counts), so a real *S. cerevisiae* reviewed proteome must be fetched (`fetch_uniprot.py`) and appended with `build_entrap.py`.
2. Confirm digestion, modifications, glycan list, isotope range, precursor/fragment tolerances, decoy construction, and Percolator settings match the reference run.
3. Pool more than three sceHCD files if necessary. Three files currently sit at the convergence boundary where Percolator can return either hundreds of rows or none.
4. Run a baseline-only seed sweep first. Do not run experimental arms until the baseline distribution is stable and non-zero.
5. Store database hashes, input file list, engine commit, command/config, Percolator version, seed, and evaluation output with every summary.

**Gate:** the baseline must be reproducible and comparable to the 384.6 setup, with no seed collapsing to zero because of the step-function floor. If that cannot be achieved, plasma remains a descriptive dataset rather than a decision benchmark.

**Result (2026-09-02, commit `b9f25527`, yeast entrapment database, routing fixed, three sceHCD files pooled, five seeds):** curated baseline 392 / 378 / 376 / 381 / 371 glycoPSMs (mean 380 ± 8, FDP 0-1.6%), no seed at the floor — this reproduces the 384.6 ± 23 reference. Default-policy baseline 374 / 352 / 284 / 342 / 364 (mean 343 ± 35): the curated policy is the better one on plasma, the reverse of mouse. Gate passed; plasma is a decision benchmark again. Comparability note: the reference used E. coli entrapment (10.77× sequon correction), this database uses yeast (3.19×).

### Stage 4 — Fix model fallback independently

Every current glyco run in both regimes falls back to the Astral model because the model store has no matching high-resolution HCD tryptic model. This is pre-existing and equal across benchmark arms, so it does not explain the arm-to-arm election result; it can still miscalibrate the entire benchmark.

0. First ask why the store's existing `hcd_qexactive_tryp` (HCD, Q Exactive, Trypsin) does not match a key of `activation=HCD, instrument=HighRes, enzyme=Trypsin`: the auto-detected generic `HighRes` instrument is not a store instrument name, so this is most likely a routing bug in the selection key, not a missing model.
1. Add an explicit high-resolution HCD tryptic fallback policy, with a loud resolved-model report.
2. Prefer a correctly trained high-resolution HCD tryptic model if suitable training data exist; otherwise define and validate the least-wrong deterministic fallback.
3. Test model routing separately from the election and additive columns.
4. Re-establish the shipped baseline after the routing change before rerunning redesign arms.

**Gate:** correct model selection is regression-tested, and the model-only A/B reports correct IDs and true FDP. Do not attribute a model-routing change to the scoring redesign.

**Result (2026-09-02, commit `b9f25527`, mouse, five seeds):** the cause was routing, not a missing model — the generic `HighRes` key fell through to the alphabetically first same-activation model. With the sibling-resolution step (`HighRes` ↔ `QExactive`) the baseline resolves to `hcd_qexactive_tryp`: correct identifications 3137 → 3183, true FDP 1.04% → 0.92%, all five seeds above the previous maximum. Gate passed; regression test in `select.rs`.

### Stage 5 — Preserve and commit the verified work

The current tree contains the redesign, two review passes, the actual-pairing guard, and the entrapment evaluator fix. Before committing:

1. Wait for the v2 result so the documentation records the final disposition accurately.
2. Run the workspace tests, glyco goldens, flag-on guard, evaluator tests, and clippy. The verification log must include the two harness fixes that made the benchmark readable: `eval_entrap.py` now sizes the entrapment space by the same substring test it counts hits with, and a Percolator result under 50 rows is reported as `PERCOLATOR_FAILED` rather than as a data point.
3. Confirm default-off redesign flags leave the shipped default output byte-identical.
4. Review the diff for generated benchmark artifacts and exclude transient files.
5. Commit correctness fixes and experimental redesign code in separable commits where practical; do not present the election as validated.

**Gate:** clean verification log and documentation that distinguishes retained correctness fixes, neutral additive features, and the refuted experimental selector.

## 3. A second election attempt, only after calibration prerequisites

Do not iterate on the current election merely by changing its weight. A new attempt requires the two missing elements already identified by the redesign document:

1. **Exclusive Y0/Y1/Y2 split anchor.** A peak cannot support several competing backbone splits simultaneously. Implement explicit peak assignment and an absence penalty where oxonium evidence predicts that anchor ions should exist.
2. **Fitted glycan contribution.** Replace `--glyco-gp-g = 1.0` and the placeholder Y-tree/oxonium priors with parameters fitted on training data that are separate from the final benchmark.

Before an engine A/B:

- measure split-level calibration and reliability, not only top-1 accuracy;
- report correct- and wrong-winner score distributions;
- require improvement in correct top-1 selections at a fixed false-winner rate;
- freeze the fitted parameters before evaluating mouse or repaired plasma holdouts.

**Stop condition:** if the exclusive anchor plus calibrated glycan score cannot improve correct split selection offline, do not implement another end-to-end election.

## 4. Independent speed track

Speed work must not be mixed with election experiments. Candidate pruning changes the competition distribution, so every non-byte-identical speed change needs the same correctness and FDP gates as a scoring change.

### Stage S0 — Define the claim

First name the PATH. Two different speed problems exist and must not share a plan: the standard DDA search (3.9×/2.0× behind Comet on Astral/TMT; measured profile 89% tree inference, 2% candidate work — see the 2026-09-02 speed research memo) and the glyco path (≈16× behind MSFragger-Glyco; candidate-generation-bound). The stages below are for the GLYCO path. The Stage S1 gate (no inference kernel before a profile) is already satisfied for the standard path and must not block that work.

Decide whether the public target is:

- **2× faster than the current Andes glyco path**, or
- **2× faster than Comet's fragment-index mode** on a named, like-for-like workload.

The existing benchmark describes Andes glyco as approximately 16× slower than its comparison engine. If that engine remains the target, “2× faster than it” requires approximately a **32× Andes speedup**. Do not use “2× faster” without naming the competitor, version, mode, dataset, hardware, thread count, search space, and whether index construction and rescoring are included.

### Stage S1 — Measure before redesigning inference

Spend one bounded profiling day on the exact compiled ensemble and production binary before implementing a custom bitvector/tree traversal.

Instrument and report:

- fragment-index construction and per-spectrum query time;
- backbone solving;
- DB-lattice candidate generation;
- phase-one `score_psm` and edge-score calls;
- unique versus repeated `(peptide, charge, scoring-view)` evaluations;
- selector evidence computation;
- winner feature extraction;
- candidates and postings visited per scan;
- thread scaling and peak memory.

Profile native binaries. The current local macOS release is x86_64/Sandy Bridge running through Rosetta on Apple Silicon, so it is suitable for before/after checks but not for a hardware-potential or competitor claim.

**Gate:** an attributed CPU and allocation profile on at least one HCD and one paired ETD/EThcD dataset. No handwritten inference kernel before this profile shows compiled tree traversal is material.

### Stage S1A — Resolve the retrieval/scoring tolerance coupling first

The code audit found a concrete high-resolution retrieval problem that must be measured before broader optimization:

- `andes.rs` hard-codes `fragment_tol_da = 0.5` when constructing `PreparedSearch`;
- `GlycoCtxOwned::build` passes that value into `FragmentIndex::build`;
- `FragmentIndex::query` then admits peptide b/y matches at `abs(observed - theoretical) <= 0.5 Da`;
- this happens even when the glyco/Y/c-z evidence paths use a high-resolution ppm tolerance.

The 0.5 Da value belongs to the currently selected rank-scoring model. It should not automatically define the tolerance of a high-resolution candidate-retrieval index. At m/z 500, 20 ppm is 0.01 Da, so the present retrieval window is about fifty times wider on each exact comparison. The practical effect may be many more posting visits, chance b/y matches, and wrong peptide-first candidates. This is a code-derived hypothesis; quantify it on the mouse spectra before claiming it explains the observed wrong winners.

Implement the experiment as a clean decoupling, not by changing the rank model:

1. Add a separately named peptide-fragment **retrieval** tolerance, resolved from instrument/high-resolution settings and reported in the run metadata.
2. Keep the legacy rank scorer at its validated tolerance and model while this retrieval-only A/B runs. The previous harmful `--tight-highres-scoring` result does not test this hypothesis.
3. Make the index support a per-ion ppm acceptance window, with a small absolute floor, while keeping a sufficiently fine fixed binning scheme for lookup.
4. Compare 0.5 Da retrieval with 20 ppm retrieval on identical candidates and spectra.
5. Record postings visited, peptides above the six-match threshold, valid peptide/glycan hypotheses, true peptide candidate rank, correct winners, nominal IDs, and true FDP.
6. Report the interaction with the existing per-spectrum peptide-first candidate cap: a narrower window changes WHICH candidates survive the cap, not only how many, so record recall both before and after the cap.

**Gate:** adopt the decoupled tolerance only if it preserves at least 99% of independently correct baseline candidates before final scoring and does not reduce correct IDs at matched true FDP. If recall falls, inspect the lost spectra by fragment charge, peptide length, and fragmentation regime before widening the window.

**Result on plasma (yeast entrapment database, three files pooled, five seeds):** curated policy 380 ± 8 → 399 ± 23 with the 20 ppm window, default policy 343 ± 35 → 365 ± 25; wall time per file 6.6-8.5 min → 1.1 min (**7×**). Identifications neutral to slightly up on both policies; entrapment hits are 0-5 per seed, so FDP is not resolvable at this scale. The additive-column arms are neutral on plasma as well (374 ± 50 and 373 ± 36 against 343 ± 35).

**Result (2026-09-02, commit `b9f25527`, mouse, five seeds, `--glyco-retrieval-tol-ppm 20` vs the same binary at 0.5 Da):** correct identifications 3198 vs 3183 (overlapping seeds, neutral to slightly positive), true FDP 0.97% vs 0.92%, and **wall time 21 min vs 144 min over the six fractions (6.9× faster)**. The identification half of the gate is met; the candidate-recall diagnostic (step 5) is still to be produced before the ppm window becomes the high-resolution default.

### Stage S2 — Quick semantics-preserving reductions (time-boxed)

Implement and benchmark separately:

1. Materialise one candidate-evidence record containing rank, edge, hyperscore, matched-ion count, ladder, paired rank, and c/z evidence. Reuse it in election, fallback, evidence gates, and output.
2. Make `hyper` and `matched_ions` consume the same `hyperscore_psm_with_matches` result instead of traversing ions twice.
3. Compute the best c/z glycosite once and derive all c/z columns from that localization.
4. Reuse target Y-tree topology and node matches when producing the mass-shifted decoy twin.
5. Reserve local maps from measured counts.

**Gate:** byte-identical PIN output on goldens and representative runs, plus a measured wall-time improvement. Keep only changes that survive repeated timing.

### Stage S3 — Replace hash-heavy kernels

1. Replace the backbone solver's two `HashMap` aggregation passes with a reusable flat vote buffer followed by sort/reduce, or a sparse dense-bin workspace with generation stamps.
2. Store the peptide fragment index as contiguous CSR-style postings rather than one heap-allocated vector per hash bin.
3. Replace per-spectrum candidate `HashMap`/`HashSet` bookkeeping with generation-stamped arrays indexed by candidate id.
4. Batch theoretical Y-node m/z values and merge them against sorted peaks rather than binary-searching per node and charge.

**Gate:** identical candidates and scores under differential tests, bounded memory on the largest benchmark database, and improvement on production-sized spectra rather than only the 120-scan fixture.

### Stage S4 — Specify the fragment-first primary cascade

The present peptide fragment index is an additive rescue path; the full precursor-minus-glycan lattice is still generated and phase-one scored. Lowering `glyco-backbone-top-k` therefore does not remove the dominant work.

Specify and instrument a bounded two-route cascade, but do not activate its pruning until the Stage S5 retrieval-recall gate passes:

```text
spectrum
   |-- fragment index --> supported peptides --> glycan by subtraction
   `-- strong Y ions ---> supported backbone masses --> peptide mass lookup
                                  |
                                  v
                         bounded candidate union
                                  |
                                  v
                           full scoring once
```

Use the exhaustive DB lattice only as a measured fallback for scans where both evidence routes are inadequate. Candidate caps and tie-breaking must remain target/decoy symmetric.

Measure at each threshold:

- fraction of baseline winning targets and decoys retained;
- independent-truth candidate recall;
- correct IDs and true FDP after Percolator;
- candidate evaluations and wall time.

**Gate:** a complete shadow-mode accounting of which baseline candidates each route and the fallback would retain. This architecture is the only identified route with enough leverage for a roughly 32× target, but its pruning remains disabled until Stage S5 demonstrates no material loss of correct glycoPSMs or glycopeptides at matched true FDP.

### Stage S5 — Make retrieval selective enough to improve IDs as well as speed

The current peptide-first score is an unweighted count of matched b/y ions with a fixed minimum of six. It queries every generated peak, including peaks likely to be glycan-derived, and uses per-spectrum hash sets/maps to enforce one-peak/one-ion counting. This is both allocation-heavy and vulnerable to common-ion and glycan-peak coincidences. A smaller candidate set is useful only if it retains the true peptide more often than the wrong high-scoring competitors.

Evaluate the following retrieval scores offline, in this order:

1. **High-resolution count:** the existing distinct-ion count with the decoupled ppm tolerance.
2. **Channel-masked count:** exclude only peaks confidently assigned to known oxonium ions or core-Y evidence; do not apply a blanket low-m/z cutoff.
3. **Rarity-weighted vote:** weight a fragment match by inverse posting frequency, for example `log((N + 1) / (df_bin + 1))`, so a rare sequence-specific ion contributes more than a common fragment bin.
4. **Rarity plus intensity:** add a bounded intensity term so one dominant peak cannot overwhelm several independent sequence ions.
5. **Mass-conditioned vote:** before visiting/scoring a posting, require its peptide mass to fall in the union of precursor-minus-known-glycan windows for the scan's charge/isotope hypotheses.

The mass condition is the high-leverage change. Build a peptide-mass-sorted index, derive allowed peptide intervals from each precursor and glycan mass, and expose those candidates as a bitmap/generation-stamped array to the fragment query. This fuses two constraints that are currently applied in opposite order: the fragment index first votes over the entire sequon database, and glycan-by-subtraction filters only afterward.

For every score, produce an oracle table on truth-labelled scans:

| retrieval diagnostic | required measurement |
|---|---|
| true peptide recall | fraction retained at each candidate budget |
| true peptide rank | median, p90, p99 and lost-tail spectra |
| selectivity | postings visited and candidates retained per scan |
| competition | number of target and decoy candidates entering exact scoring |
| final accuracy | correct IDs and true FDP after the unchanged scorer/Percolator |
| final speed | phase-one calls, wall time, CPU and peak RSS |

Use a confidence cascade rather than a universal hard cap:

- high-confidence mass-conditioned peptide retrieval becomes the primary route;
- strong core-Y evidence contributes a complementary glycan-first route;
- their bounded, label-blind union receives exact scoring once;
- low-confidence scans fall back to the exhaustive lattice.

Confidence must be based on frozen, target/decoy-symmetric observables such as absolute rare-ion vote, number of independent ions, and the score margin. Do not use the unsafe rule “skip the full path if any glycan evidence exists”: the current core-Y prefilter already demonstrates that a wrong evidence-bearing backbone can suppress a weak true one.

This structure follows the successful indexing ideas in the primary literature without copying another engine's score: MSFragger uses precursor-mass ordering inside its fragment index to score candidates simultaneously, while pGlyco 2/3 use coarse glycan/core-Y filtering before expensive fine scoring. The Andes-specific contribution should be a calibrated two-route union that preserves weak-b/y and weak-Y spectra through different routes.

### Stage S6 — Remove exact-scoring duplication after the candidate pool is fixed

The audit also found repeat work that is safe to optimize but unlikely to deliver the full speed target alone:

1. `score_psm` walks peptide residues to accumulate prefix masses; `psm_edge_score` immediately allocates two prefix arrays and walks the same residues again. Add a shared prefix-mass/nominal sidecar or a combined scoring entry point, preserving the two output features exactly.
2. Phase one creates a `Vec<usize>` of mass-bucket candidates for every backbone before scoring it. Replace this with a borrowed range/iterator or a reusable scratch buffer.
3. Charge-indexed scored spectra are repeatedly found with linear `.iter().find(...)` calls inside candidate loops. Resolve them once into a small charge-indexed table.
4. If the same bare `(candidate, charge, scoring-view)` survives through several equivalent backbone hypotheses, cache its phase-one score. Decorated ETD variants require glycan mass and site in the key and should be measured separately.

**Gate:** each item needs a byte-identical PIN golden and a production-sized timing result. Merge these independently so a semantic regression cannot hide inside the architectural candidate change.

### Stage S7 — Revisit winner quality only after retrieval is calibrated

The mouse result says the present split election amplifies wrong winners. Do not combine its retry with the retrieval work. First freeze a candidate generator that has high truth-set recall and materially fewer competitors. Then measure whether candidate multiplicity alone reduces the wrong-winner tail under the shipped selector.

If an election is still needed, compare candidates using **differential evidence**: shared oxonium and shared Y peaks should cancel, while exclusively explained Y0/Y1/Y2 and peptide fragments should decide the pair. Fit any glycan weight on training data, freeze it, and evaluate on the mouse holdout at matched true FDP. This is the earliest point at which the redesign election should be reopened.

### Stage S8 — Proposed implementation slices

Keep the work reviewable and bisectable:

1. **Instrumentation only:** timings, posting counts, candidate recall/rank dump; no search-result changes.
2. **Retrieval-tolerance decoupling:** ppm-aware index query behind a flag.
3. **Mass-conditioned index:** allowed-candidate mask, still using the existing unweighted count.
4. **Retrieval-score study:** masking and rarity/intensity weights, offline first.
5. **Primary cascade:** bounded two-route union plus explicit exhaustive fallback.
6. **Scoring-kernel cleanup:** shared prefixes and materialised finalist evidence.
7. **Data-structure kernel:** CSR postings and generation-stamped counters, after semantics are frozen.

Each slice must report both speed and independent-truth accuracy. The promotion criterion is multiplicative: retain correct identifications at matched true FDP **and** reduce total wall time. A fast arm that changes the target/decoy competition without passing the FDP gate is a failed scoring experiment, not a speed optimization.

## 5. Reporting template

Every completed experiment should add one row to a durable summary:

| field | required value |
|---|---|
| hypothesis | one sentence, defined before the run |
| code | commit and dirty/clean state |
| data | dataset, files, database hash, glycan list |
| configuration | full config/command and model selected |
| repetitions | seeds and whether files were pooled |
| accuracy | nominal IDs, correct IDs, true FDP, unique glycopeptides |
| performance | index, search, rescoring, total wall, CPU, RSS, threads |
| verdict | supported, neutral, refuted, or unreadable |
| next action | ship, retain experimental, revise prerequisite, or stop |

## 6. Immediate queue

1. ~~Let the current v2 matrix and chained evaluations finish; read it at 36/36.~~ Done (2026-09-02 evening).
2. ~~Write the v2 result into this plan and the redesign document's implementation-status section.~~ Done.
3. ~~Close or reopen the election strictly by the Stage 2 gate.~~ Closed as refuted.
4. ~~Reconstruct the like-for-like plasma entrapment benchmark.~~ Done: baseline 380 ± 8 reproduces 384.6. Database built (`fasta/plasma_entrap_yeast.fasta`: deposited human + 6733 reviewed yeast as `ENTRAP_`; sequon correction 3.19×); baseline-only seed sweep queued.
5. ~~Fix and independently benchmark high-resolution HCD model routing.~~ Done: +46 correct, see Stage 4.
6. Verify and commit the current tree with an honest experimental-status message.
7. Define the named speed target, then add Stage S1 instrumentation before changing inference or candidate generation.
8. ~~Run the Stage S1A retrieval-only tolerance A/B on mouse; keep the rank model unchanged.~~ Done: 6.9× faster, IDs neutral; recall diagnostic still owed.
9. Produce the Stage S5 oracle table before implementing a primary fragment-first cascade.

## 7. Primary algorithm references

- Kong et al., [MSFragger: ultrafast and comprehensive peptide identification in shotgun proteomics](https://pmc.ncbi.nlm.nih.gov/articles/PMC5409104/) — theoretical fragment indexing with precursor-mass ordering and simultaneous candidate scoring.
- Polasky et al., [Fast and Comprehensive N- and O-glycoproteomics analysis with MSFragger-Glyco](https://pmc.ncbi.nlm.nih.gov/articles/PMC7606558/) — mass-offset glyco search made practical by fragment-ion indexing.
- Liu et al., [pGlyco 2.0 enables precision N-glycoproteomics with comprehensive quality control](https://pmc.ncbi.nlm.nih.gov/articles/PMC5585273/) — coarse glycan/core-Y candidate filtering followed by fine scoring.
- Zeng et al., [Precise, fast and comprehensive analysis of intact glycopeptides and modified glycans with pGlyco3](https://pmc.ncbi.nlm.nih.gov/articles/PMC8648562/) — glycan-first complementary-Y indexing and explicit core-ion handling.
