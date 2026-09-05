# Glyco scoring redesign — factored per-candidate likelihood with a two-stage election

> **Status 2026-09-05.** Items 4 and 5 (the split election and its `--glyco-gp-g` weight)
> were measured as losing (table below) and have been **removed from the code**, together
> with `--glyco-gp-m`. The additive-column items (Y-tree, oxonium LLR, masked rank/chance
> LLR) remain as hidden research flags, measured identification-neutral. This document is
> kept as the design record; the flags it names in those sections no longer exist.

**Date:** 2026-09-02. **Branch audited:** `feat/glyco-pin-curation` @ `69fbf354`.
**Inputs:** six independent read-only audits (candidate generation / mass split; fused selector;
glycan-evidence scoring; backbone scoring; PIN + decoy interface; competitor literature at formula
level) cross-checked against the measured record in the session memory.
**Target:** intact N-glycoPSMs and glycopeptides @1% PSM FDR (Percolator) above MSFragger-Glyco on
the like-for-like plasma benchmark (andes 254 → curated-PIN 385 vs MSFragger 587) and on PXD011533.

---

## 1. Where the IDs are lost (measured, consolidated)

| Stage | Finding | Status |
|---|---|---|
| Generation | 97% of truth glycopeptides are generated; every generation-side widening tried moved yield DOWN. | MEASURED |
| Emission | 90.5% of emitted rows sat on scans with no glycopeptide; a RawScore>3 floor removes 83% of them. Curated PIN + floor: 256.8 → 384.6 @1%. | MEASURED |
| Selection | On contested scans 42.6% of winners are decoys; 96.9% of those sit at a DIFFERENT backbone mass than truth. Split election is 98.7% correct on MSFragger-confident scans; failures concentrate on MSFragger's own low-margin scans. | MEASURED |
| Selector | `S = rank + 10·ladder + 5·core_y + 1·hyper (+15·cz)`. `ladder`/`core_y` are per-split constants (5–38 distinct values per 600 candidates); the only per-candidate discriminator with real power (RawScore, truth median rank 2) is computed AFTER the argmax, for the winner only. | MEASURED |
| Percolator | One row per scan, no runner-up. Glycan-decoy twin copies 45 of ~48 features verbatim → Percolator finds no direction (0 PSMs, 5 seeds). | MEASURED |
| Composition | No composition-diagnostic oxonium ion influences which composition wins (gate only, `--glyco-sialic-oxonium-min-frac` defaults 0). Isobars resolved by float sort order. 2.5 composition strings per glycan mass vs Byonic 1.0. | MEASURED |
| Y-ladder | Single linear chain; Fuc appended after all antennae, so Y1+Fuc / Y2+Fuc are never theorised; no absence penalty; base-peak normalised (base is usually the 204 oxonium). | code-verified |
| Pairing | Under `--glyco-hcd-pair`, `core_y` reads the HCD partner while `ladder`, `YHitFrac`, `SialicConsistency` read the ETD scan. One argmax, two spectra. | code-verified |
| Backbone | Bare backbone; hyperscore discards intensity; no oxonium/Y masking before ranking (rank model trained on plain tryptic peptides sees 20–40 intense glycan peaks above the backbone ions); full-glycan decoration REPLACED bare masses and regressed (−41%); HexNAc-stub inside the LLR never tried. | code-verified / MEASURED |

The unifying picture: **the split (peptide mass vs glycan mass) is never chosen; it falls out of a
global argmax in which each of ~17 mass hypotheses independently rolls a weak sequence-score die,
padded by per-split constants.** A decoy in a peptide-dense window wins the max-order statistic.
Percolator then receives the single winner with no margin information on either axis.

## 2. What every top engine does that andes does not

From the formula-level review (Polasky 2020/2022, Zeng 2021, Liu 2017, Lu 2020, Fang 2022):

1. A **glycan score separate from the backbone score, consulted at selection time.** MSFragger-Glyco's
   2020 search mixed Y ions into one hyperscore and its own Methods name that as a sensitivity
   hazard; the 2022 PTM-Shepherd glycan assignment is what fixed it.
2. **Composition-diagnostic oxonium with an absence penalty** (PTM-Shepherd: hit ratios >1, miss
   ratios <1, per NeuAc/NeuGc/Fuc class; pGlyco3 gates compositions on diagnostic ions).
3. **Y evidence square-root normalised** so large glycans do not win by rung count; only ions
   UNIQUE to a candidate score in pairwise comparison.
4. **Per-peak mass-error weighting** (pGlyco `(1−merr/tol)^4`; PTM-Shepherd `log|Δm/σ|`).
5. **Backbone and glycan combined with a fitted weight** (pGlyco `Score_GP = w·G + (1−w)·P`,
   w = 0.35 fitted by ranking SVM), with `ratio_core` as a heavily weighted term.
6. **Mass-shifted-Y glycan decoys and two-axis error control** — expressed in andes as PIN features,
   never a second FDR engine.

## 3. The proposed scoring system

Everything below is additive to the existing PIN (RankScore and the ~50 closed-search columns are
untouched). The SELECTOR is replaced, not reweighted — reweighting a linear fusion of per-split
constants is measured exhausted.

### 3.1 Stage 0 — spectrum preparation (once per scan)

- Oxonium profile `O = {o_m}` for monosaccharide classes m ∈ {HexNAc, HexNAc-Hex, NeuAc (274/292),
  NeuGc (290/308), Fuc (512.197, 350.145), HexNAc2 (407.166), NeuAc-Hex-HexNAc (657.235)} as
  base-peak fractions. Used as a GATE and as an input to per-candidate terms; never as a summand.
- Peptide-channel view of the spectrum: a peak list with oxonium m/z, the charge-reduced precursor,
  and the Y-rung m/z of the candidate's OWN backbone mass excluded, re-ranked by intensity, with its
  own base peak. The mask depends on backbone mass only, so it is identical for a target and its
  reversed decoy (decoy-symmetric).
- Under HCD/ETD pairing: the HCD partner supplies the glycan channel, the ETD scan the c/z channel.
  Every term declares which spectrum it reads; no term mixes them.

### 3.2 Stage 1 — per-candidate log-likelihood ratios (for ALL accepted candidates, before any argmax)

For candidate c = (peptide p, composition g, backbone mass m):

```
P(c)  backbone LLR     = RankScore(p)                               [existing, untouched]
                        + λ_I · ChanceLLR_masked(p)                  [new: Σ_matched p_ion·w_I·(−ln p_chance)
                                                                      on the peptide-channel spectrum;
                                                                      glycosite-spanning ions = max over
                                                                      {bare, +HexNAc, +2HexNAc} forms;
                                                                      fragment charge 1..min(z−1,5) with
                                                                      isotope confirmation for z≥3;
                                                                      length-normalised as cz already is]
                        + ETD: CzChanceLlr(p,g)                       [existing form; the only c/z variant
                                                                      that ever separated]

G(c)  glycan LLR       = Σ_t 2·(a_hit^t·√U_t + a_miss^t·√V_t)        [Y-tree: U/V = matched/missed nodes of
                                                                      g's OWN composition-specific Y set
                                                                      (core, core+Fuc, antennae, sialic-loss
                                                                      shadows); classes t = {Hex/HexNAc-only,
                                                                      Fuc-containing}; intensity-weighted hits]
                        + Σ_m [ g∋m ? log P(o_m|m present)/P(o_m|absent)
                                    : log P(o_m|absent)/P(o_m|present) ] [oxonium composition consistency,
                                                                      presence AND absence, per class]
                        + β · log(σ_unmod / |Δm_ppm|)                [per-candidate mass-error term]
                        + c · log(ratio_core)                        [pGlyco core term]

A(c)  split anchor     = Y0/Y1/Y2 evidence at backbone mass m under EXCLUSIVE peak assignment:
                         a peak already claimed as Y_k of a competing split with higher intensity-
                         weighted likelihood cannot also count for m. Absence of Y1 when the scan
                         has strong oxonium is penalised.
```

All three are log-odds, so they sit on one scale and need no K/J/H weights. `a_hit`, `a_miss`, the
oxonium ratios and `σ_unmod` are fitted on OWN data (logistic regression of ion presence on
target-vs-entrapment labels), not copied from PTM-Shepherd's supplement.

### 3.3 Stage 2 — split election (multiplicity-controlled)

```
S_split(m) = logsumexp_{c: bb(c)=m} [P(c) + A(c) + G(c)]  −  log n_peptides_in_window(m)
m* = argmax_m S_split(m);   DeltaBackbone = S_split(m*) − second-best split
```

The `−log n` term removes the max-of-N lottery (F6): a split with 40 sequon peptides in its window
no longer gets 40 draws. Keep the top-2 splits for the next stage. A de-novo (composition-less)
candidate may enter only with an explicit `log P(novel)` prior; the enum-fallback that today
promotes a losing enumerated runner-up on ~22% of scans is deleted.

### 3.4 Stage 3 — glycoform and peptide election inside m*

```
c* = argmax_{c: bb(c)=m*}  P(c) + w·G(c)        w fitted by ranking SVM on labelled scans (pGlyco: 0.35)
DeltaGlycan  = G(c*) − best G of a different composition at m*
DeltaPeptide = P(c*) − best P of a different peptide at m*
DeltaTD      = P(c*) − best decoy peptide at m*
```

Composition ties (NeuAc/NeuGc, Hex+Fuc isobars) are now decided by G, i.e. by oxonium consistency and
the composition-specific Y set, instead of by float sort order.

### 3.5 Stage 4 — emission and PIN (one row per scan, Percolator PSM-only)

Keep the emission floor (measured). New ADDITIVE columns on the single winner row:

| column | definition | axis |
|---|---|---|
| `BackboneLLR`, `ChanceLLRMasked`, `RankScoreMasked` | P(c*) and its new parts | peptide |
| `GlycanLLR`, `GlycanYLLR`, `OxoniumCompLLR`, `YHighPriorMissing`, `MassErrLLR` | G(c*) and its parts | glycan |
| `DeltaBackbone`, `NSplitsInWindow`, `CompositionDensity` (log #compositions fitting the residual) | split | both |
| `DeltaGlycan`, `NCompAtBackbone` | glycoform | glycan |
| `DeltaPeptide`, `DeltaTD` | sequence | peptide |
| `GlycanDecoyGap` = `GlycanLLR − GlycanLLR(mass-shifted-Y twin, RECOMPUTED with the same tree)` | glycan-axis null as a FEATURE, not a row | glycan |
| `DeltaMassLogFreq` | empirical log-frequency of the winner's glycan mass across the run (additive analogue of PeptideProphet's extended mass model) | glycan |

Retire once the above is live (after one A/B, via the dead-column guard): spectrum-level
`OxoniumScore`/`NCoreOxoniumIons`, `TailorScore` (≡ RankScore), `IsGlycanDb` (constant 1),
the axis-mixed `DeltaRankScore`, and raw `GlycanMass` (size-only, rewards big decoy-fitting glycans).

## 4. Why this should exceed MSFragger, and where it will not

- MSFragger's search score still mixes Y ions into one hyperscore and relies on a post-hoc mass prior
  and a post-hoc glycan assignment; here the factored per-candidate likelihood decides the election
  itself and hands Percolator margins on all three axes. That is a strictly richer decision than
  either engine's search stage.
- Where the gap is not scoring: MSFragger spends more error budget (FDP 0.85% vs 0.39%, few events)
  and is ~16x faster. Speed is a separate track; do not trade it here.
- Historical refutations do NOT apply as-is: promoting `partial_glycan_by`/`y0y1` into the OLD linear
  fusion regressed because they were added on top of per-split constants; the design above removes
  those constants and replaces the argmax, which has never been tested.

## 5. Falsifiable build order (one variable per step, pooled fractions, 5 seeds, entrapment FDP)

Seed floor is ~117 PSMs; the gap to MSFragger is ~200 from the curated baseline, so the target is
measurable but individual steps may not be. Use the offline re-ranking gate to size steps first.

0. **Offline re-rank gate (zero engine changes, hours).** On the existing `--debug-glyco` candidate
   dumps (plasma 334 scans × ~600 candidates with MSFragger-confident labels; mouse frac2) compute
   truth-top-1 rate and margin for: current selector; `rank + RawScore`; `+ exclusive Y0/Y1/Y2
   anchor`; `+ oxonium consistency`; `− log n_window`. Any step that does not raise top-1 on the
   low-margin stratum is dropped before it is built.
1. **Wrong-spectrum and dead-code fixes** (cheap, correctness): `SialicConsistency` on `gen_peaks`;
   ladder and `YHitFrac` read the HCD partner under pairing; delete enum-fallback; ETD scans no
   longer bypass the oxonium gate into a full 600-split enumeration. Measure junk fraction and IDs.
2. **`RankScoreMasked` PIN column** (peptide-channel spectrum, `score_psm` unchanged). Tests the
   rank-poisoning hypothesis directly; if Percolator weights it, stages 3.1/3.2 build on the mask.
3. **Composition-specific Y-tree + oxonium composition LLR as PIN columns** (`GlycanYLLR`,
   `OxoniumCompLLR`, `YHighPriorMissing`, `GlycanDecoyGap` with the recomputed twin). Fit the class
   ratios on own labels. Retire `OxoniumScore`/`SialicConsistency` after the A/B.
4. **Split election (3.3) with `DeltaBackbone`/`CompositionDensity`** — the stage where 96.9% of
   the wrong winners live. Measure decoy-winner fraction on contested scans (223/524 today) before
   IDs.
5. **Glycoform/peptide election (3.4) with fitted w** and the delta columns.
6. **HexNAc-stub max-over-forms inside `ChanceLLRMasked`** and isotope-gated fragment charges to
   z−1, PIN-only first, then into P(c).
7. Only then: own-data glyco intensity / c/z model as an additional LLR term.

A/B on BOTH regimes every time (mouse PXD011533 6-frac pooled and plasma PXD030622 pooled); a
plasma-only win that empties mouse of decoys has happened before (664 → 0).

## 6. Provenance of claims

Code locations behind each finding are in the six audit reports (session 2026-09-02); the headline
ones: selector `glyco_psm.rs:75-116` and `glyco_search.rs:2150-2165`; post-argmax feature block
`glyco_search.rs:2245-2510`; linear ladder `backbone.rs:599-612`; oxonium gate `oxonium.rs:48-69`
and the unused sialic gate `oxonium.rs:298-306`; pairing spectrum mix `glyco_search.rs:1203-1209,
1516`; glycan-decoy twin copy `glyco_pin.rs:168-345`; no peak masking `scored_spectrum.rs:423-430`;
generation split union `hybrid.rs:441-556`.

---

## 7. Implementation status (2026-09-02)

Everything below is in the tree, builds clean, and passes the workspace suite plus
clippy. **Every new flag defaults OFF and the default PIN is byte-identical**: both glyco
goldens were regenerated and verified column-by-column to show all 77 pre-existing
columns unchanged on every row, with the 11 added columns all zero.

| Step | State | Flag |
|---|---|---|
| 0 Offline re-rank gate | `scripts/glyco_rerank_gate.py` | — |
| 1 Wrong-spectrum / gate defects | done | `--glyco-pair-y-on-gen`, `--glyco-enum-fallback`, `--glyco-etd-require-oxonium` |
| 1 `SialicConsistency` on the gate's spectrum | done, unconditional | — (byte-identical unpaired) |
| 2 Peptide-channel rank | done | `--glyco-rank-masked` |
| 3 Composition-specific Y-tree | done | `--glyco-y-tree` |
| 3 Oxonium-composition LLR | done | `--glyco-oxonium-llr` |
| 4 Split election with multiplicity control | done | `--glyco-split-election` |
| 5 Glycoform/peptide election and margins | done | same flag, `--glyco-gp-g` |
| 6 HexNAc-stub max-over-forms in the backbone LLR, isotope-gated charges to z−1 | done, PIN-only (`ChanceLlrMasked`, `ExplainedMasked`) | `--glyco-chance-llr-masked` |
| 7 Own-data glyco intensity / c/z model | NOT built | — |

### Review pass (2026-09-02, three independent reviewers on the uncommitted tree)

Fixed in the tree after the review, all verified by the regenerated goldens (90
pre-existing columns byte-identical on every row of both) and by a new flag-ON
guard (`crates/andes/tests/glyco_flags_on.rs`) that runs the fixture with every
redesign flag set and requires each new column to take ≥2 distinct values:

- **Y-tree antennae ignored core fucose.** Every antenna node of a Fuc-containing
  composition was emitted bare, predicting up to ~24 ions a core-fucosylated
  glycopeptide can never show, so the tree penalised fucosylation by that miss count.
  Antenna nodes now carry the +Fuc form when the composition has a fucose.
- **Noise-level hits cancelled a full miss penalty.** A hit now needs ≥1% of the base
  peak (`MIN_HIT_FRAC`); below that the node is a miss.
- **Decoy twin could land on a target node.** Shifts are now symmetric ±(1–30) Da and
  re-drawn when within 0.5 Da of any target node.
- **Multiplicity control counted rows, not peptides.** `−log n` now uses the number of
  DISTINCT peptides in the split, with the logsumexp over each peptide's best score.
- **Split key was a fixed grid.** Two masses 5 ppm apart could straddle a bucket edge;
  splits are now single-linkage clusters of the deduped backbone masses
  (`glyco_psm::split_ids_by_clustering`, swept test over 2000 masses).
- **De-novo candidates entered the election with G = 0**, a systematic bonus over
  every enumerated candidate that misses a node. They are excluded from the election
  until a `log P(novel)` prior is fitted.
- **Election mixed spectra under pairing.** The Y-tree half of G read the ETD scan
  while the oxonium half read the HCD partner; both now read the HCD partner.
- **Oxonium half of G was gated on a different flag** (`--glyco-oxonium-llr`); the
  profile is now built whenever the election is on.
- **Mask missed Y isotopologues.** M+1/M+2 of every masked Y rung are now masked too.
- **`YTreeDecoyGap` was a label leak** with `--glyco-decoy` (0.0 on every decoy row);
  it is now negated on the decoy row, as `SialicConsistency` is.
- **803.29 oxonium was mislabelled** as HexNAc2+Hex+Fuc (that is 715.28); it is the
  sialyl-Lewis ion and needs NeuAc, so it is dropped from the Fuc class.
- **The offline gate's multiplicity rule was dead** (keyed on the intact mass, which
  every candidate shares). Rewritten: backbone-mass clustering, distinct-peptide `n`,
  a low-margin stratum, (file, scan) keys for pooled input, `shipped` = emitted row
  order, and two circular-truth guards.

### Second review pass (external, same day) — fixed

- **The offline gate's "shipped" order was not the shipped order.** Under
  `--debug-glyco` the rows were sorted by the fused selector, then passed through a
  HashMap and re-sorted by `RankScore`. The fused position is now remembered and
  restored for the diagnostic dump (the collapse path keeps its rank-DESC order).
- **Masked peaks still reached the isotope deconvolution** (rank `u32::MAX` made them
  unmatchable but they could seed or be consumed by clusters and shift a neighbour's
  m/z). The masked path now deconvolutes the surviving peaks only; `new` is unchanged.
- **Composition-only fucose was treated as known core fucose.** `GlycanComp` carries
  counts, not topology, so an antenna-fucosylated glycopeptide was penalised by every
  +Fuc antenna node it cannot show. The Y-tree is now a max over the two placements
  (`y_node_topologies`), target and decoy twin alike.
- **Glycan-decoy rows copied the target's `YTreeHitFrac` and `YTreeHighPriorMissing`.**
  The twin's own values are stored and emitted.
- **`DeltaGlycan`/`DeltaPeptide` were margins of the fused score.** They are now
  axis-specific: `DeltaPeptide` is a difference of P (rank + hyper + cz) over different
  peptides in the split, `DeltaGlycan` a difference of G over different compositions.
- **Calibration** (finding 6) stands as a limitation, not a fix: the class priors, the
  5% chance rate and the oxonium logistic parameters are placeholders, so `--glyco-gp-g`
  = 1.0 is a starting point, not a fitted weight. Fit on target/entrapment labels before
  reading the fused score probabilistically.

Deviations from §3 that stand, on purpose, and must be read before an A/B:

- The election score is `rank + H·hyper + CZ·cz + G·(Y-tree + oxonium)`; the anchor
  term `A(c)`, the per-candidate mass-error term and the `ratio_core` term of §3.2 are
  NOT built. `H` keeps its shipped default 1.0 (set `--glyco-gp-h 0` to test the pure
  log-odds form).
- Stage 3 runs on the top-1 split only; the runner-up split contributes `DeltaBackbone`.
- `DeltaGlycan`/`DeltaPeptide` are margins of the fused election score, not of `G`
  and `P` separately, and one weight `--glyco-gp-g` serves both stages.
- `DeltaTD`, `CompositionDensity`, `NCompAtBackbone` and `DeltaMassLogFreq` of §3.5
  are not implemented (the election is label-blind, so `DeltaTD` would have to be a
  writer-side feature). §3.5's `GlycanYLLR`/`YHighPriorMissing`/`GlycanDecoyGap`/
  `NSplitsInWindow` are the code's `YTreeLLR`/`YTreeHighPriorMissing`/`YTreeDecoyGap`/
  `NSplitsConsidered`.
- The Y-tree uses six structural classes (Y0, core, core+Fuc, antenna, sialylated,
  intact) with placeholder priors, not the two classes of §3.2; the class ratios are
  still to be fitted on own labels.
- `--glyco-split-election` has no effect under `--debug-glyco` (a multi-row dump has
  no single winner to attach margins to); the engine now warns.

New modules: `andes_glyco::glycan_y_tree` (composition-specific Y node set with
sqrt-normalised hit/miss LLR and a mass-shifted-Y twin, 7 tests),
`andes_glyco::oxonium_llr` (per-class diagnostic-ion LLR with the absence penalty and a
cofragmentation down-weight hook, 15 tests), `ScoredSpectrum::new_with_excluded_mz` in
the scoring crate (masks m/z windows before ranks, base peak and noise density are
computed; `new` proved unchanged field-by-field, 7 tests).

New PIN columns, default-policy only (the curated set stays pinned to its measured
configuration, and a guard test enforces both halves of that): `YTreeLLR`,
`YTreeHitFrac`, `YTreeHighPriorMissing`, `YTreeDecoyGap`, `OxoniumCompLLR`,
`RankScoreMasked`, `MaskedPeakCount`, `ChanceLlrMasked`, `ExplainedMasked`,
`DeltaBackbone`, `DeltaGlycan`, `DeltaPeptide`, `NSplitsConsidered`.

### Two defects found by building it

- **The split key was a constant.** Bucketing a backbone mass by dividing it by a
  mass-PROPORTIONAL tolerance evaluates to `1/ppm` at every mass, so the whole lattice
  collapsed into one split and the election had nothing to elect between. Caught only
  because the emitted split count was constant 1 on a fixture holding 123 splits. Now
  `glyco_psm::split_bucket`, in log-mass space, with four unit tests.
- **The PIN header and the row writer drifted apart.** The new block was declared after
  `SialicConsistency` in the header and written before `CoreYHits` in the writer,
  shifting five unrelated columns on every row. Caught by diffing the regenerated golden
  column-by-column rather than accepting it. The golden's own regeneration recipe was
  also missing `--glyco-taxon human`, which the test passes; following the recipe
  produced 25 wrong rows. Both the recipe and a warning are now in the test's header.

### MEASURED (2026-09-02, Codon, five Percolator seeds, pooled fractions, sequon-corrected entrapment)

Mouse PXD011533 rep1 Frac1-6 against the 5088-PSM deposited truth set. `def0` = default
policy + RawScore floor; arms add one flag set each. v1 = the binary before the second
review pass, v2 = after it.

| arm | glycoPSMs @1% | correct @1% | true FDP |
|---|---|---|---|
| def0 baseline | 4833 | 3137 | 1.04% |
| + `--glyco-rank-masked --glyco-chance-llr-masked` (v2) | 4833 | 3137 | 1.16% |
| + `--glyco-y-tree --glyco-oxonium-llr` (v2) | 4814 | 3130 | 1.02% |
| + `--glyco-split-election` (v1) | 5387 | 3079 | 4.64% |
| + `--glyco-split-election` (v2) | 5200 | 3049 | 3.71% |

**The split election (§3.3-3.4) is refuted as built.** It raises the nominal count while
lowering the number of correct identifications, at 3.6-4.5x the true error, on every seed
in both binaries; the review fixes reduced the damage without changing the direction.
Under the curated policy it is worse with no yield gain (4361 / 2916 / 0.97% vs 4409 /
2960 / 1.18%). The additive columns (§3.5) are identification-neutral. Plasma PXD030622
(sceHCD R1-R3) was at the noise floor in this setup (0-335 IDs within one arm; the
entrapment database was not the 384.6 one) and decided nothing. Full record:
`docs/plans/glyco/50-roadmap/2026-09-02-redesign-benchmark-and-speed-plan.md`.

### What is NOT measured

The election's yield and error are measured above. Still unmeasured: the additive
columns on a plasma benchmark that can resolve them, and everything under the model
fallback caveat in the roadmap (every run used `hcd_astral_tryp` by fallback).

One caution from the flag-ON guard run on the 120-scan low-res fixture (synthetic
E. coli background, NOT a benchmark): with every redesign flag set the election
changed the winner on 23 of 120 scans and the number of decoy-labelled winners went
from 3 to 9. On this fixture the direction is unfavourable. It is not a measurement,
but it is the first thing the pooled A/B must look at (decoy-winner fraction on
contested scans, step 4 of §5) before any yield number is read. The
smoke test of the re-rank gate used this repo's own collapsed winners as "truth", which
is circular and is not evidence. Run step 0 against real reference labels first, then
the pooled five-seed end-to-end A/B with entrapment FDP, on BOTH regimes.
