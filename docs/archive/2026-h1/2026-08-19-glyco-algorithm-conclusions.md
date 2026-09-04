# Glyco algorithm: what is measured, what is refuted, where the headroom is

**2026-08-19.** Every number below is measured, not argued. Where a mechanism was traced
in code but never measured, it says so. Several confident mechanisms in this document were
*wrong* — those are kept, because the pattern of which arguments failed is the most useful
thing here.

Benchmark substrate throughout: **PXD030622 human plasma sceHCD R1-R3**, the depositors'
own `uniprot-human20410.fasta`, pooled across fractions, one Percolator run at seed 42,
scored against their deposited Byonic results and against a 24,942-entry human + E. coli
`ENTRAP_` database.

---

## 1. The headline: generation is NOT the glyco bottleneck

**Six** levers were measured on plasma with entrapment, and **not one raised yield.** Four
expanded the candidate space, one changed the feature model, one gated compositions on
oxonium evidence (§2b — its apparent +34% was a partial-gate artifact and is retracted).
**The only change in the whole campaign that ADDED identifications was a REMOVAL:** dropping
240 biologically impossible NeuGc compositions, +36% at 3.4x lower true error. That
asymmetry is the central result — expansion is paid for in threshold tightening, and only
removal of what should never have been enumerated pays. Baseline is `--glyco-no-neugc`: 365 glycoPSMs @ 0.55%.

| lever | what it expands | glycoPSMs | entrapment FDP |
|---|---|---|---|
| baseline | — | **365** | 0.55% |
| `--glyco-glycan-list reference-human` (312 -> 1,229 comps) | glycan space | **228 (-37%)** | 0.00% |
| `--glyco-y-index` (AXIS 2 two-axis retention) | surviving backbones | **340 (-6.8%)** | 0.29% |
| `--glyco-isobar-rep` | (selection) | 226 vs 247 on the plain DB | — |
| `--glyco-decorated-features` | (feature model) | **215 (−41%)** | 0.00% |

**Mechanism.** A larger surviving candidate pool gives decoys more places to fit,
Percolator's threshold tightens, and real identifications are left behind. FDP heading
toward 0.00% is the tell that a search has become over-conservative, not more correct.
This reproduces pGlyco's "the glycan database is not the larger the better"
(182-comp DB -> 0.8% glycan FDR vs 1,234 -> 4.0%) from the yield side.

**This contradicts the recorded decomposition** of the plasma gap as ~12.5% generation /
~14% selection. On plasma the generation share measures at zero or negative.

**Corollary that matters more than the finding.** In a target-decoy + Percolator pipeline,
candidate-space expansion is **not free** — it is paid for in threshold tightening.
Expansion only pays when the added candidates are *enriched for truth*. The proof is the
mirror image: the one change that ADDED identifications was a **removal**.

### 1a. AXIS 2 is nearly redundant with AXIS 1
Emitted glyco rows per fraction, baseline vs `--glyco-y-index`: 7,205 -> 7,223,
7,700 -> 7,717, 7,631 -> 7,656. **+60 rows on ~22,500 = +0.27%.** Glycan-Y retention keeps
almost exactly the backbones peptide-b/y retention already kept. An audit had rated this
the "highest expected plasma yield per line changed"; it was not.

---

## 2. The one big glyco win: NeuGc is a species error, not a tuning knob

`n_glycan_list_common()` enumerated `neugc 0..=1` unconditionally.

**The arithmetic is exact, not approximate.** Fuc = Hex − O and NeuGc = NeuAc + O, so
Hex1NeuAc1 and Fuc1NeuGc1 are the SAME elemental formula (C17H27NO13):
`(Hex+NeuAc) − (Fuc+NeuGc) = 0.000000000000 Da`. **No mass accuracy at any resolving power
separates them.** (An earlier claim of a "~1 µDa gap" in this campaign was an artifact of
the code's own 6-dp rounded constants.)

Measured on the shipped list: **600 compositions over only 460 distinct masses; 140 masses
(30%) carry more than one composition; 100% of those collisions involve NeuGc.** Excluding
NeuGc gives 360 compositions over 360 masses — **zero collisions, by construction**.

**Biology.** Humans lack a functional CMAH gene and cannot synthesise NeuGc (Chou et al.,
PNAS 1998; Alu-mediated, ~2.8 Mya). Dietary NeuGc is incorporated but sits ~10,000x below
NeuAc in serum (3.3 pg/µL, Seo et al. 2021) and is "essentially undetectable on human plasma
proteins" (Muchmore et al. 1998). **Mice have functional CMAH** — which is why the
mouse-developed glyco benchmarks never surfaced this.

**Measured, entrapment-validated:**

| | default | `--glyco-no-neugc` |
|---|---|---|
| glycoPSMs @1% | 268 | **365 (+36%)** |
| entrapment hits | 5 | **2** |
| **entrapment FDP** | **1.87% — OPTIMISTIC** | **0.55% — CONSERVATIVE** |

**More identifications at 3.4x lower true error, from one change.** Both axes improving
together is rare and is the strongest available evidence that this is a real fix.

**The mouse counter-arm closes the argument (2026-08-21).** A smaller search space is not
automatically better, so the human result alone does not prove the *species* claim. Run the
same flag where NeuGc is genuine and it must HURT. It does. Mouse liver, 3 fractions pooled,
entrapment DB:

| arm | flags | glycoPSMs @1% | glycopeptides | entrapment FDP |
|---|---|---|---|---|
| A | `--glyco-taxon mammal` (NeuGc kept) | **20,515** | 4,404 | 0.19% CONSERVATIVE |
| B | `--glyco-no-neugc` | **19,112 (−6.8%)** | 3,985 (−9.5%) | 0.20% CONSERVATIVE |
| **C** | **`--glyco-taxon auto` — the SHIPPED DEFAULT** | **20,515** | **4,404** | **0.19% CONSERVATIVE** |

**Arm C is identical to Arm A digit for digit**, so the shipped default costs mouse users
nothing. Auto got there from two independent signals that both said keep: 99.75% of sialylated
spectra carry a NeuGc oxonium at >=10% of NeuAc, and the FASTA is 17,267/17,267 CMAH-competent.
The same logic narrows the list on human plasma. The species gate is therefore safe to ship ON.

FDP is flat (0.19% -> 0.20%), so this is pure yield loss, not the usual tighter-gate trade.
+36% on human plasma and −6.8% on mouse liver is precisely the asymmetry CMAH predicts and
that a generic "fewer candidates" story cannot produce. **Consequence for the release:
`--glyco-no-neugc` is a species assertion, never a global default.**

**It also dissolves the recorded "plasma FDR floor."** That was attributed to glycan SIZE
creating a dense composition space. It was not: the space was dense because half the list
was a NeuGc shadow of the other half, and those shadows were exactly where the decoys were
fitting. Remove them and the FDR estimate becomes honest.

**Field check.** Every mainstream human-labelled glycan list ships ZERO NeuGc — verified by
direct file inspection: FragPipe Human-253/708 and its shipped glyco-N workflows,
GlycReSoft, MetaMorpheus, GPQuest, Protein Prospector "Human", and Byonic's 160-entry list
as deposited in PXD030622. NeuGc appears only in lists labelled *mammalian* or *mouse*.
andes was the only engine including it by default. Klein et al. (Nat Commun 2024) warn
verbatim that applying a mammalian glycome to a human sample "could skew results towards
incorrect compositions".

**Keep it opt-in, never inferred from the FASTA.** Dietary NeuGc is real and can be enriched
in carcinomas, and recombinant human protein from murine/CHO hosts genuinely carries it —
GlyCounter measured 30% FEWER IDs when excluding NeuGc from an NS0-expressed sample.

### 2b. Sialic oxonium gating fixes CALIBRATION, not yield (and a retraction)

`--glyco-sialic-oxonium-min-frac` keeps NeuGc in the list but admits a NeuAc/NeuGc claim
only when the matching oxonium (274/292 or 290/308) clears a fraction of base peak.

⚠ **A first measurement of 358 glycoPSMs @ 0.84% FDP, reported as +34%, is RETRACTED.** A
five-lens review found the gate was applied only inside `db_branch`, while three other
routes push `Source::Db` ungated — the peptide-first union (ON BY DEFAULT), the
glycan-Y-first block, and cross-spectrum transfer — and that it filtered AFTER the argmin,
dropping the annotation instead of falling through to the non-sialic isobaric twin. Pruning
one generator of four biased WHICH GENERATOR WON, and that rebalancing produced the apparent
gain, not the sialic evidence.

Re-measured with the gate covering all four generators and filtering before the argmin:

| arm | glycoPSMs | entrapment FDP | verdict |
|---|---|---|---|
| default (NeuGc in, ungated) | 268 | **1.87%** | **OPTIMISTIC** |
| `--glyco-no-neugc` (species exclusion) | **365** | 0.55% | CONSERVATIVE |
| sialic gate 2%, full | 267 | **0.00%** | CONSERVATIVE |
| sialic gate 5%, full | 241 | **0.00%** | CONSERVATIVE |

**What it actually does:** converts an OPTIMISTIC 1.87% into a CONSERVATIVE 0.00% at no
yield cost versus default. The default's q-values genuinely understate the real error and
the gate repairs that — a correctness gain worth having. But it buys no identifications,
and FDP pinned at 0.00% at both thresholds is this document's own tell for
over-conservatism. Monotonic in the wrong direction: 2% -> 267, 5% -> 241. The untested and
informative direction is LOOSER (0.005, 0.01), not stricter.

**`--glyco-no-neugc` remains the best plasma lever (365 @ 0.55%). The gate stays
default-off.**

**★ The sharpest methodological lesson of the campaign lives here:** a partially-applied
feature produced a large apparent win by biasing which subsystem dominated, with nothing to
do with its stated mechanism. The measurement was honest about the code as it stood — the
code was not doing what its name said. **Verify a feature's SCOPE before trusting its
effect size.**

### 2a. You cannot resolve from fragments what should never have been enumerated
`--glyco-isobar-rep` resolves the same collisions on Y-ladder evidence. Measured: **−21
PSMs**, and compositions-per-mass got *worse* (2.47 -> 2.70), because the resolution is
per-spectrum — two spectra of the same glycan still pick different members of the pair.
Removing the impossible compositions works; adjudicating between them does not.

---

## 3. Where the headroom actually is: separation, not coverage

At **0.29–0.55% measured true error against a nominal 1%**, the search is leaving room on
the table. That is a separation problem.

The audit located it precisely: the fused selector's two heaviest terms — `ladder`
(weight 10) and `core_y_hits` (weight 5) — are **pure functions of `bb_hit_idx`**. They are
mathematically identical across every peptide competing at the same backbone mass, so they
contribute **zero peptide-level discrimination**. On HCD `gp_cz` is 0, which leaves
`rank + 1.0·hyperscore` doing all the real work against per-backbone terms of magnitude
30–60.

Meanwhile the genuinely discriminative quantities — `RawScore`, `strong_score`,
`IntensitySignal` — are computed **after** the collapse has already committed.

**Lever 2 was measured and is REFUTED — and the chemistry explains why.**

Decorating the ~40-column PIN feature vector with the glycan (so glycosite-spanning
fragments sit at peptide+glycan mass) was rated a top structural defect by audit. Measured:

| | baseline | + decorated features |
|---|---|---|
| glycoPSMs @1% | **365** | **215 (−41%)** |
| unique glycopeptides | 161 | 115 |
| entrapment FDP | 0.55% | 0.00% |

**The premise was wrong.** Under HCD a glycopeptide fragments at the GLYCOSIDIC bonds
first: the dominant products are oxonium ions, Y-ions (peptide + partial glycan), and b/y
ions of the backbone AFTER the glycan is lost. **The bare backbone is therefore the correct
theoretical ladder for b/y ions.** Decorating moves half the predicted ladder to
peptide+glycan masses where there are few or no peaks. This also finally explains the
pre-existing `let decorate = false` on the scoring path, which carried a measured −16 with
no rationale attached.

**That leaves ONE untested selection-side lever:** get a per-candidate discriminator into
the collapse. `RawScore` / `strong_score` / `IntensitySignal` are computed AFTER the
collapse has committed, so the best evidence never influences which candidate is emitted.
This is the only structural fix in the campaign that has not yet been measured.

---

## 4. Non-glyco scoring defects found by the same campaign

**`DeltaRankScore` was measuring noise.** It tracked "best vs second-best distinct peptide"
but keyed distinctness on `nominal_residue_mass()` — and the candidate window IS a
nominal-mass window, so every genuine sequence competitor took the same-mass branch and
never set `second_raw`. The emitted delta was the lead over a candidate 1–2 Da away, i.e. a
different isotope hypothesis. Re-keyed to peptidoform identity:
**UPS1 +1.36% entrapment-validated** (16,034.8 -> 16,253.2 PSMs, 5/5 seeds, arms
non-overlapping; of 218 added PSMs only 8.4 are entrapment hits, so **96.2% are real**).

Note the direction: the old column was not *inflated*, it was **noise**. An incorrect feature
is not a biased feature, it is an uninformative one — replacing it ADDS discrimination.

**High-res models were trained through a ~50x-too-tight window.** `IonType::mz` rebuilds
theoretical m/z from the INTEGER nominal node mass, displacing it by a **median 52 ppm**
(90th pct 133 ppm; measured over 300,333 simulated b-ion positions). Training then matched
that displaced position inside **20 ppm**, so only ~20% of theoretical positions were
matchable and the rest were recorded "missing" even when the peak was present.

Measured in the shipped store (rank_dist missing-slot frequency, prefix/suffix):

| model | prefix | suffix |
|---|---|---|
| cid_lowres_tryp (trained @ 0.5 Da) | 0.452 | 0.421 |
| **hcd_qexactive_tryp** (trained @ 20 ppm) | **0.928** | **0.828** |
| hcd_astral_tryp | 0.854 | 0.773 |
| etd_highres_tryp | 0.888 | 0.836 |

A clean split on `is_high_resolution()`. The learned absent-ion penalty collapses to
**~−0.04 nats** (vs ~−0.8 low-res), so a candidate can miss nearly every predicted fragment
almost for free and the score degenerates to a sum of positive matched evidence.
**This is a plausible mechanism for the long-standing "closed-search scoring is at ceiling"
result, which was measured on exactly these high-res models.** Fixed on the training side
(exact m/z + match at the model's own `mme`), verified serve-neutral; **the payoff requires
a retrain**, gated on the missing-slot frequency falling 0.93 -> ~0.45.

**GATE PASSED, with a control (2026-08-21).** A fresh corpus was harvested from PXD009875
(8 flats, 145,385 PSMs) and `hcd_qexactive_tryp` was trained from it TWICE, changing nothing
but the fragment-match window:

| model (identical corpus) | mean missing | median |
|---|---|---|
| control, `--fragment-tol-ppm 20` (old behaviour) | 0.8953 | **0.9268** |
| fixed, seed `mme` = Da(0.5) | **0.5732** | **0.6491** |
| shipped `hcd_qexactive_tryp` (defective) | 0.9063 | 0.9284 |
| shipped `cid_lowres_tryp` (healthy reference) | 0.5213 | 0.5837 |

**Replicated on a 3.6x larger corpus (526,234 PSMs / 40 flats):** control 0.9234, fixed
**0.5563** (mean 0.4946). The decisive detail is that the CONTROL DID NOT MOVE — both arms got
the identical data increase, the control shifted 0.0034 while the fixed arm shifted 0.0928. So
corpus size cannot explain the fix; only the window can. The fixed model is now below the
healthy low-res reference (0.5563 vs 0.5837).

The control reproduces the shipped defect to within 0.002 on data the shipped model never saw,
which rules out the corpus as the cause and pins it on the window. The fix lands the statistic
in the same regime as the healthy low-res model. (The "~0.45" target was an approximation;
low-res itself is 0.52/0.58.) **The training-window commit stays.**

**This is a training statistic, not an identification gain.** Nothing has been searched with the
retrained model yet. The absent-ion penalty should now carry real information, but whether that
converts into PSMs at 1% FDR is UNMEASURED, and on this campaign's record (5 of 7 levers
negative) that is not a safe assumption. The next step is to serve the retrained store on a
high-res benchmark, multi-seed, against the shipped model.

Also corrected: `ANDES_TIGHT_HIGHRES` tightened SERVING to 20 ppm (Astral 36,719 -> 28,894)
as a "train/serve mismatch" fix. Serving was never the broken side — its 0.5 Da window
absorbs a 0.036 Da displacement.

---

## 5. The recurring defect class: computed, then not consulted

Six instances found in one campaign. All the same shape — correct behaviour implemented,
then made unreachable by a literal, with a comment describing behaviour the code no longer
has:

1. `let yindex_on = false` — AXIS 2 retention, comment promised "opt-in for a clean A/B"
   with no flag existing.
2. `let isobar_rep = false` — evidence-based isobar resolution, comment claimed it was
   "removed rather than left as an opt-in switch"; it was neither.
3. `let decorate = false` — decorated scoring peptide.
4. **Sialic oxonium ions**: `NEUAC_OXONIUM_MZ` / `NEUGC_OXONIUM_MZ` are defined,
   `neuac_obs`/`neugc_obs` are computed — and **never consulted for composition
   assignment**. This is the field's actual mechanism for breaking the NeuGc/NeuAc
   degeneracy (pGlyco: require 274/292 for NeuAc, 290/308 for NeuGc).
5. `edge_score` returning 0 on every low-res model via an `error_scaling_factor` guard.
6. `--refine` PSMs emitted with synthetic `BASEPEP_<n>` accessions — `merge_into_pass1`
   offsets the protein index into the concatenated DB but never maps back, so **every
   modified PSM from a refinement run was unattributable to a protein**. Fixed.

**⚠ Oxonium gating needs an intensity threshold, not a presence test.** Chalkley & Baker
(MCP 2025) measured 40,466 mouse-liver spectra carrying the m/z 290 NeuGc oxonium among
glycopeptides containing NO NeuGc — ~70% of spectra with a NeuGc oxonium had none, from
co-isolation. PTM-Shepherd's calibration is the model to copy: hit/miss probability ratios
NeuAc 2/0.05, NeuGc 2/0.05, **dHex 2/0.5** — fucose absence weighted 10x weaker than sialic
absence.

---

## 6. Proteoform inference: what is and is not buildable

Asked whether andes could report a probability for the most prominent proteoform from
bottom-up data with a one-gene-one-protein DB.

**Not buildable honestly: a joint proteoform probability composed from independent
per-site marginals.** Histone H4 has >10^10 theoretical proteoforms from 58 PTMs at 17
sites; seven labs measuring intact H4 observed **75 above 0.01% abundance**
(Aebersold et al. 2018). Eight orders of magnitude — a product-of-marginals model is not
miscalibrated there, it is the wrong functional form. The field has never settled whether
that gap is biology or detection limits, so "nature is sparse, the argmax is probably right"
rests on the open question.

**The state of the art is a clean dichotomy.** Every method that NAMES a proteoform gets
its evidence from intact mass / top-down (Proteoform Suite). Every purely bottom-up method
refuses to name it, outputting anonymous "proteoform groups" of co-varying peptides (COPF,
PeCorA, ProteoForge, AlphaQuant, BP-Quant). COPF is the only calibrated one and it
calibrates "do these peptides split into >=2 populations", not "is proteoform X present".

**What IS solid.** Bottom-up and middle-down **agree well on marginal PTM abundances**
(Sidoli et al. 2015). The failure is specifically in the **joint**. So per-site occupancy is
defensible; the joint is not.

**The genuine opening.** The measured joint exists and has never been transferred: interplay
scores `I = log(f12/(f1·f2))` from middle-down histone work, CrossTalkDB, Apache-licensed
implementations. Transferability **r = 0.47–0.68 across cell lines** vs ~0 for randomised
PTM codes — real but moderate, arguing for **shrinkage toward independence, not a hard
prior**. Conditionals **flip sign**: ΔI(K36me2|K27un) = +0.11 vs ΔI(K36me2|K27me2) = −0.37.
No published method learns co-occurrence from top-down/middle-down and transfers it into
bottom-up inference. (Caveat: "not found by a thorough search", not exhaustively excluded.)

**Benchmark substrate already exists:** ABRF **iPRG-2016**, "Inferring Proteoforms from
Bottom-up Proteomics Data" — 8 submissions, task "perceived as difficult", and **the
ground-truth dataset is public**.

**Terminomics is the one axis where the phasing objection does not apply** (termini are
observed, not inferred across a broken linkage) — but its dominant failure mode is NOT FDR.
In-source fragmentation accounts for ~22% of semi-tryptic peptides; sample prep swings
partial-tryptic IDs 28.4% -> 2.8%; trypsin storage buffer swings non-specific activity
20% -> 1%. Those are real peptides, correctly identified, **wrongly interpreted** —
invisible to target-decoy. And Nt-acetylation does NOT prove a mature terminus: downstream
initiation sites are *more* acetylated (87% vs 72%).

Also: **FDR cannot be inherited across aggregation levels.** A PrSM-level FDR under-reported
the true protein-level FDR by **24-fold** (LeDuc et al. 2019). A proteoform-level error rate
must be estimated at the level the claim is made.

---

## 7. Methodological lessons (the part most likely to save time)

1. **A traced mechanism is a hypothesis, not a result.** In this campaign, mechanisms that
   were confirmed in code and refuted by measurement: the glyco default decoy strategy
   (recommended `sequon-reverse`; entrapment showed the DEFAULTS win, 1066 @ 0.47% vs 974 @
   0.72%), three candidate-space expansions, and the isobar resolver. Mechanisms confirmed
   AND validated: NeuGc removal, DeltaRankScore.
2. **Yield alone will ship a bad change, and so will FDP alone.** Read them together. FDP
   falling toward 0.00% means over-conservatism, not correctness.
3. **A q-value gain is not an identification gain** until entrapment says so. One glyco
   result in this campaign looked like +28% on q-values and had to be re-checked before it
   held up (it did, at +36% and 3.4x lower true error).
4. **Predictions were wrong more often than right, and in a consistent direction** — "more
   candidates -> more IDs" failed three times. In a target-decoy pipeline the intuition is
   backwards.
5. **Verify agent findings before acting.** Several were correct and load-bearing (the
   training-window defect, the neutral-loss path gap found by CodeRabbit); one rated the
   highest-value lever (AXIS 2) which measured negative.
6. **Check the harness before believing the result.** A `SERVE DIFFERS` verdict in this
   campaign came from a deleted baseline binary and `cmp` failing on missing files; a golden
   "regenerated" per its own documented recipe differed because the recipe named `test.mgf`
   while the test runs `test.mgf.gz`.

---

## 8. Next actions, ranked by evidence

1. **Retrain the high-res models** with the training-window fix. Gate on missing-slot
   frequency 0.93 -> ~0.45 *before* anything downstream. Then re-run one previously-refuted
   discriminator (RS³ / RichIonLLR): if a shelved discriminator comes back positive, the
   "scoring is at ceiling" conclusion was an artifact.
2. **Intensity-thresholded sialic oxonium scoring.** The ions are already computed; the
   calibration to copy is published; it is the evidence-based version of the species gate
   and works where NeuGc genuinely belongs.
3. **Selection-side, not generation-side:** decorated PIN features, and a per-candidate
   discriminator inside the collapse.
4. **Do NOT invest further in glycan-axis decoys.** pGlyco1 built one, pGlyco3 abandoned it
   for hard evidence gates, and andes's `--glyco-decoy` is independently broken (emits a
   second row per scan, violating the top-1-per-scan rule its own writer documents as
   required for a valid per-scan FDR).

## Definitive baseline vs MSFragger (2026-08-26)

The first glyco comparison run at **matched error control**. Everything before this compared an
FDR-controlled andes count against a score-thresholded Byonic count, which is not a comparison.

| | andes | MSFragger 4.2 | Byonic |
|---|---|---|---|
| glycoPSMs @1% FDR | 254 | **587** | 539 (score threshold, no FDR column) |
| unique glyco peptides | 64 | **91** | 79 |
| entrapment FDP | 0.39% | 0.85% | not measurable |
| wall, 3 files | ~4,200 s | **260 s** | - |

**andes reaches 43% of MSFragger's glycoPSMs and 70% of its peptides, ~16x slower.** Same
database, tolerances, decoys, Percolator invocation, entrapment proteins and extraction code;
MSFragger searched the same composition space (1000 N-glycan offsets, NeuGc excluded).

**This closes a question the earlier analysis left open.** andes matches or beats Byonic on
distinct entities (51 vs 43 compositions, 81% of peptides), which raised the possibility that the
recall "gap" was an artifact of comparing against an unknown threshold. Against a properly
FDR-controlled competitor that possibility is eliminated: the gap is real.

**It is also not a tuning problem.** Nine measured levers are negative or neutral - six
generation-side expansions, two selector re-weightings (down-weighting the heavy terms costs IDs:
254 -> 110 and 254 -> 217), and the decoy strategy (sequon-reverse 224 vs reverse 254,
replicating the earlier mouse result on human). Knobs cannot reach 587.

**andes is more conservative, not less accurate**: FDP 0.39% vs 0.85%. Some of MSFragger's lead is
spending more of the error budget, but not most of it.

WARNING - **every FDP here rests on very few entrapment hits.** andes's 0.39% is 1 hit in 253; the
Poisson 95% CI is 0.01%-2.20%. The earlier claim of "2.5x unused error budget" is NOT supported.

WARNING - **MSFragger's raw accepted count is 3,126, of which 2,539 are UNMODIFIED** (mass offset
0.0 is in the list so plain peptides can win). Only 587 carry a glycan. Comparing the raw number
would overstate the gap ~5x. The PIN has no delta-mass column; glyco PSMs are recovered as
`ExpMass - theoretical_peptide_mass > 800`.

The remaining direction is architectural, not parametric: the ~11 constant PIN columns that give
Percolator nothing to separate on, a graded glycan score, and log-domain intensity terms that a
linear model can actually use.
