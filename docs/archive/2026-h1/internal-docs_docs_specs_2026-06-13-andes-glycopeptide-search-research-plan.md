# Intact-glycopeptide search mode in andes — Research Plan

**Date:** 2026-06-13
**Status:** Research / feasibility plan — not yet scoped for implementation
**Author basis:** feasibility study of the andes codebase (post neutral-loss milestone) + intact-glycopeptide literature/engine survey

> This is a **research roadmap**, not an implementation spec. It records the
> landscape, a feasibility verdict against andes's actual architecture, what
> reuses the just-shipped neutral-loss primitive, a phased approach, and what to
> prototype first. Implementation of any phase requires its own
> brainstorm → spec → plan cycle. Companion to the crosslinking research plan
> ([2026-06-13-andes-crosslinking-research-plan.md](2026-06-13-andes-crosslinking-research-plan.md)).

## 0. Verdict up front

An intact-glycopeptide search **mode** is **viable on andes and a much better fit
than crosslinking** — because, unlike a crosslink, a glycopeptide is **still one
peptide**. The glycan is a single labile precursor adduct on a normal backbone, so
andes's core invariant `suffix_mass = peptide_nominal − prefix_mass`
([psm_score.rs:234-311](../../crates/scoring/src/scoring/psm_score.rs)) **holds at
the peptide-backbone mass**. There is no two-peptide co-scoring rewrite. The
problem is not the scoring core; it is **(a) determining the backbone mass when the
precursor = peptide + unknown glycan**, **(b) a glycan composition search space**,
and **(c) a glycan-level FDR layer** — all of which sit *around* the existing
engine rather than rewriting it.

- A **Phase-0 "glyco-lite" mode already ships today** as a stopgap: declare a few
  glycan compositions as variable mods + `loss=` (the neutral-loss milestone on
  `feat/enzyme-support`). It finds only the *declared* compositions, searches them
  combinatorially (precursor window per composition), and has **no glycan-level
  FDR**. Useful for a known single glycoform; not a glyco search engine.
- A **real N-glyco, stepped-HCD, oxonium-triggered MVP** is a **multi-week
  feature** that reuses the linear engine almost entirely. The non-obvious new
  parts are the **glycan composition list + candidate expansion**, **Y0/Y1
  backbone-mass determination**, and **2D (peptide × glycan) FDR**. It warrants a
  **`--glyco` mode** (and likely an `andes-glyco` crate for the glycan DB + glyco
  FDR), but **not a separate engine** the way non-cleavable crosslinking did.

**Recommendation:** pursue **Phase 0 (already shipped) → Phase 1 prototype
(oxonium gate + Y-ion backbone-mass determination, *measurement only*, no scoring
change)**. Gate the full MVP (Phases 2–4) on the Phase-1 prototype showing andes
can recover backbone masses on a real stepped-HCD glyco file (target dataset
**PXD025455**, Q Exactive HF). This is the make-or-break, exactly as the doublet
detector was for crosslinking.

## 1. Motivation

Intact glycopeptide analysis identifies, at one shot, *which peptide*, *which
site*, and *which glycan composition* — the highest-value, hardest problem in
glycoproteomics. andes today cannot do it: it has no glycan search space, no
backbone-mass-from-Y-ion logic, and no glycan FDR. The neutral-loss milestone
gave andes the *fragment-scoring* primitive for labile mods but explicitly **not**
a glyco search engine (design doc §"Out of scope":
[2026-06-13-andes-neutral-loss-glyco-design.md](2026-06-13-andes-neutral-loss-glyco-design.md)).
This plan answers: *what would a real glyco mode take, what does it reuse, and is
it the right shape for andes?*

## 2. Background — intact glycopeptides for an engine designer

- **The species reaching the MS** is one peptide bearing one (N-glyco) or several
  (O-glyco) glycans. Precursor mass = `peptide backbone + glycan(s)`. The backbone
  is a normal tryptic peptide; the glycan is a large (often 1–3 kDa), **labile**
  adduct that fragments *off* the backbone under collisional activation.
- **N-glyco** attaches at the **N-X-S/T sequon** (X ≠ P), one glycan per peptide —
  the tractable case. **O-glyco** attaches at S/T with **no sequon and multiple
  glycans per peptide** — the localization-combinatorics case (O-Pair's domain;
  out of MVP scope).
- **Three fragment families in one stepped-HCD scan** (the key enabler):
  1. **Oxonium ions** — low-m/z glycan markers (HexNAc **204.087**, **138.055**,
     Hex-HexNAc **366.140**, NeuAc **292.103 / 274.092**). They *flag* a spectrum
     as glyco but carry no peptide info.
  2. **Y-ion ladder** — intact peptide + progressively trimmed glycan. **Y0 = bare
     peptide backbone**, **Y1 = peptide + one HexNAc**. The Y-ladder is how you
     read the **backbone mass** and the glycan composition.
  3. **Backbone b/y ions** — the peptide's own sequence ions, mostly
     **glycan-stripped** under HCD (so they appear at the *deglycosylated* backbone
     mass). These are the ions andes already scores — at the backbone mass, not the
     precursor mass.
- **Stepped-HCD (sceHCD, e.g. 20-30-40%)** deliberately produces all three families
  in **one MS2 scan** — no MS3, no paired ETD required. This is why the MVP targets
  it: it maps onto andes's single-MS2 model with no acquisition-pairing logic.
- **The combinatorial trap:** the precursor mass alone does not tell you the
  backbone mass, because `precursor = backbone + glycan` is one equation with two
  unknowns. You resolve it **either** by reading Y0/Y1 from the spectrum (measure
  the backbone) **or** by trying every glycan composition (subtract each from the
  precursor and search the residual mass). This is the central design fork (§4).

## 3. State of the art (engines, algorithms, licensing)

The five production engines (pGlyco3, MSFragger-Glyco, Byonic, StrucGP, O-Pair)
all instantiate the **same five-stage skeleton**, differing mainly in glycan-DB
philosophy and FDR rigor:

1. **Oxonium trigger / spectrum gating.** Universally anchored on **HexNAc
   204.087** plus {138.055, 366.140, 274.092}. Two threshold patterns:
   - **Intensity gate:** summed oxonium intensity ≥ **10% of base peak**
     (pGlyco3, MSFragger; lowered to 2.5% for AI-ETD). ([MSFragger-Glyco, Nat
     Methods 2020](https://www.nature.com/articles/s41592-020-0967-9))
   - **k-of-n peak rule:** Byonic ≥ 2 of {204.087, 274.092, 366.140} at ±0.01 Da;
     StrucGP ≥ 2 of {138.055, 204.087} among the top-20 peaks. ([Byonic glyco,
     MCP 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8724605); [StrucGP, Nat
     Methods 2021](https://www.nature.com/articles/s41592-021-01209-0))
2. **Y-ion backbone-mass determination.** Either read the **Y0/Y1 ladder**
   directly (Byonic, StrucGP), do **glycan-first Y-complementary indexing then
   `peptide = precursor − glycan`** (pGlyco3's distinctive trick — index
   `precursor − peak = glycan − glycan_Y`, require ≥ 2 trimannosyl-core Y ions for
   N-glyco, pass the top-100 glycan candidates to the peptide search), **or** an
   **open / mass-offset search** where the glycan rides as a precursor delta-mass
   and the **glycan-stripped backbone b/y ions** drive the peptide ID
   (MSFragger-Glyco, O-Pair). ([pGlyco3, Nat Methods 2021](https://www.nature.com/articles/s41592-021-01306-0))
3. **Peptide search at the backbone mass.** A conventional b/y database search at
   the recovered backbone mass, **sequon-constrained for N-glyco**.
4. **Glycan composition/structure assignment.** From a **curated composition DB**
   (pGlyco3 **182–1,234** compositions; MSFragger **182 N / 110 O**; Byonic 9–309 +
   a wildcard ±300 Da; O-Pair 12–32 combinatorial) **or** a **modular/de-novo
   reconstruction** (StrucGP: 4 cores × 17 branches; pGlycoNovo). The composition
   lists are small enough to enumerate.
5. **Combined score + multi-level FDR.** A **peptide-backbone score** (b/y
   matching) combined with a **glycan score** (Y-ladder + oxonium matching):
   - pGlyco: explicit weighted sum `ScoreGP = w·ScoreG + (1−w)·ScoreP`, **w ≈ 0.35**
     (≈ 65% weight on the peptide).
   - MSFragger / Byonic: folded into one hyperscore / Byonic-score, glycan
     refined post-hoc.
   - **FDR is the real differentiator:**
     - **pGlyco** — true **3-level** (glycan / peptide / glycopeptide) with
       spectrum-based glycan decoys + finite-mixture model, combined by
       inclusion-exclusion `FDR_GP = FDR_G + FDR_P − FDR_(G∩P)`.
     - **MSFragger/FragPipe, StrucGP, O-Pair** — a **separate glycan-level
       target-decoy** stage on top of standard reverse-sequence peptide decoys.
       FragPipe's decoy glycan keeps the nominal composition but **shifts the
       intact mass within tolerance + a random isotope error and scrambles each
       Y/oxonium fragment by a unique random 1–20 Da** (no monosaccharide
       swapping). This is the **simplest defensible glycan-decoy method** and the
       recommended andes target. ([FragPipe glycan FDR, MCP 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC8933705))
     - **Byonic** — 2D target-decoy (PSM × protein) with ML-PEP; no dedicated
       glycan-composition decoy.

**Acquisition fit:** all are MS2-based (no MS3). pGlyco3/StrucGP recommend
stepped-HCD; O-Pair needs paired HCD+EThcD; MSFragger is activation-aware. **The
MVP target (stepped-HCD, single MS2) is the most andes-friendly.**

**Licensing (relevant to andes's Apache-2.0 goal):**

| Engine | License | Open? | Safe to read clean-room? |
|---|---|---|---|
| **O-Pair / MetaMorpheus** | **MIT** | **Yes** (C#) | **Yes** — best reference; study method, don't paste code |
| **Glyco-Decipher** | **Apache-2.0** | Yes | **Yes** — same target license; ideal reference (verify repo has source not just a JAR) |
| pGlyco3 | repo says Apache-2.0 but ships **binary-only + activation-gated** | No | **Read papers only** — treat like pLink2 (license-gated) |
| MSFragger-Glyco | free academic, **commercial license** (Fragmatics/U-M); closed | No | **Read papers only** |
| Byonic | **commercial/proprietary** (Bruker) | No | **Avoid entirely** |
| StrucGP | **no license** (MAC-locked activation) | No | **Read papers only** |
| GlycoPAT | custom non-SPDX License.pdf | source avail. | **Caution** — read the PDF before any reuse |

**Patent posture (clean-room conclusion).** The **general intact-glycopeptide
method appears NOT patent-encumbered**: oxonium-triggering (HCD-pd-ETD / product-
dependent triggering) is an *acquisition* method published independently by many
groups since ~2010; Y-ion backbone determination, glycan-composition DB matching,
and glyco 2D/3-level FDR are published across **many independent groups**
(pGlyco, Byonic, O-Pair, StrucGP, Glyco-Decipher) — the hallmark of a clean-room-
implementable method. The **only located patents** are Protein Metrics' (M. Bern):
**US20100124785A1** (wildcard-modification search) and **US10546736B2 /
US20210335589A1 / US20170236697A1** (interactive MS-data analysis GUI). These are
narrow. **Guardrails:** (1) do **not** implement an open/**wildcard-modification**
search in the manner US20100124785A1 claims (the MVP's curated-composition list is
clear of it); (2) do **not** clone Byonic's interactive peak-labeling GUI. As with
crosslinking, a **patent review of the specific chosen methods is required before
implementation**, but no blocking patent on the general method was found.

## 4. Feasibility against andes's architecture — component map

For each essential component: **(a) REUSABLE as-is**, **(b) EXTENDABLE**, or
**(c) NET-NEW**, with the andes file/function.

| Component | andes today | Verdict |
|---|---|---|
| **Backbone b/y GF scoring** | `score_psm` walks one peptide, `suffix = peptide_nominal − prefix`; per-(partition,IonType) rank LLR; mass-indexed node DP ([psm_score.rs:234-311](../../crates/scoring/src/scoring/psm_score.rs), [scored_spectrum.rs](../../crates/scoring/src/scoring/scored_spectrum.rs)) | **(a) REUSABLE as-is** at the **backbone mass**. A glycopeptide's glycan-stripped b/y ions ARE the backbone peptide's b/y ions; the single-peptide invariant holds. This is the whole reason glyco fits and crosslinking does not. |
| **Candidate generation** | per-protein enzymatic walk + variable-mod combinatorics; precursor-mass-windowed bucket index `nominal(peptide.mass − H2O) → cand idxs` ([candidate_gen.rs:40-203](../../crates/search/src/candidate_gen.rs), [match_engine.rs:127-136,87-110](../../crates/search/src/match_engine.rs)) | **(b) EXTENDABLE** — the engine and the bucket index are reusable, but the **mass window must be opened at the *backbone* mass, not the precursor mass** (see §5, the hard problem). Two routes: search at the Y-determined backbone mass (one window) or at `precursor − each glycan composition` (N windows). |
| **N-glyco sequon constraint** | none — andes is sequence-agnostic in candidate gen | **(c) NET-NEW** (small): restrict glyco candidates to peptides containing **N-X-S/T (X≠P)**. A candidate filter, not an engine change. |
| **Glycan adduct representation** | `Modification` = single-residue scalar mass delta + `neutral_losses` + `loss_class` ([modification.rs:28-46](../../crates/model/src/modification.rs)) | **(b) EXTENDABLE.** A *single declared* glycan fits as a mod with `loss=`/`class=glyco` TODAY (Phase 0). A **glycan composition with its Y-ladder masses** needs a richer first-class **`Glycan` concept** (composition → monosaccharide counts → ordered Y-ion masses → total delta). Unlike crosslinking's bivalent linker, this is still a *monovalent* adduct — no relational/two-peptide concept needed. |
| **Oxonium trigger** | none — andes has no spectrum-level diagnostic-ion gating | **(c) NET-NEW** (cheap): a pre-scoring spectrum classifier — sum/count oxonium peaks vs a threshold to flag a scan as glyco and route it to the glyco path. Pure peak lookup; no scoring-core change. |
| **Y0/Y1 backbone-mass determination** | none | **(c) NET-NEW** (the make-or-break): given the precursor and the spectrum, find the Y-ion ladder and deduce the bare-peptide mass. This is the component with no andes analogue and the Phase-1 prototype target. |
| **Glycan composition DB** | none (catalog has a `glyco` *model* slug [catalog.rs:95](../../crates/model-train/src/catalog.rs), but no composition list) | **(c) NET-NEW**: a curated **N-glycan composition list** (start from MSFragger's/pGlyco's **182**) → monosaccharide counts → total mass + Y-ladder masses. Small, static, shippable as a resource file. |
| **Neutral-loss fragment scoring** | per-class pooled loss IonType (`loss_class:u8`) + activation-gated loss-shifted b/y prediction + peptide-aware additive loss pass in `score_psm` ([fragment_ions.rs:146-238](../../crates/scoring/src/scoring/fragment_ions.rs), [psm_score.rs:286-340](../../crates/scoring/src/scoring/psm_score.rs)) | **(a) REUSABLE / (b) EXTENDABLE** — see §6. This is the **glycan-score fragment primitive**: the Y-ladder rungs and glycan-retaining backbone ions are exactly "loss-shifted b/y ions tagged with a glyco loss class," already predicted, scored against a pooled per-class table, and activation-gated (correctly suppressed for ETD). The glyco mode *consumes* this, it does not rebuild it. |
| **Combined peptide+glycan score** | single integer `RawScore`; additive PIN features (DeltaRawScore, Tailor, edge score) feed Percolator | **(b) EXTENDABLE** — emit the **glycan score as additive PIN feature column(s)** (oxonium count, Y-ladder match count/intensity-ratio, glycan mass error) next to `RawScore`. Percolator/the FDR layer combines them — no need to hand-tune pGlyco's `w ≈ 0.35`; let the rescorer learn the weight. Respects the "additive-only PIN features" rule. |
| **Output schema** | one `Peptide`/`Proteins`/mass per PIN row ([pin.rs:11-20](../../crates/output/src/pin.rs)); accession in TSV-only `Modifications` column | **(b) EXTENDABLE** — one peptide per row still holds (unlike crosslinking). Add **glycan composition + site** to the **TSV** (the `Modifications` column already emits `pos:CURIE`; extend to carry the glycan composition string). PIN stays additive-feature-only (its `Proteins` column is rest-of-line — same constraint the neutral-loss work documented). |
| **FDR** | whole-protein-reversal decoys ([decoy.rs:11-25](../../crates/search/src/decoy.rs)); Percolator/TDC at the PSM level | **(c) NET-NEW layer** — peptide FDR reuses reverse decoys + Percolator, but a **separate glycan-level target-decoy** (FragPipe-style: random intact-mass shift + 1–20 Da fragment scramble per glycan) and a **2D combination** (`FDR_GP`) is a new post-search layer. The single biggest net-new piece after backbone-mass determination. |
| **Precursor calibration / isotope error** | `adjusted_observed_neutral_mass`, isotope-error window ([match_engine.rs:97-109](../../crates/search/src/match_engine.rs)) | **(a) REUSABLE** — glyco precursors are heavier and isotope-error-prone, but the existing machinery applies unchanged once the backbone mass is the search target. |

**Summary:** the **scoring core, candidate engine, decoy/Percolator peptide-FDR,
PIN feature plumbing, and the neutral-loss fragment primitive are all reusable or
extendable.** The **net-new** pieces are a *bounded* set: oxonium gate (cheap),
**Y0/Y1 backbone determination (hard)**, glycan composition DB (static resource),
sequon filter (small), and **glycan-level 2D FDR (a new post-search layer)**.
There is **no structural blocker** of the kind crosslinking hit (no two-peptide
co-scoring, no n² pairing, no bivalent-linker concept).

## 5. The hard problem specific to andes — backbone mass when `precursor = peptide + unknown glycan`

andes's candidate generation is **precursor-mass-windowed**: it opens a bucket
window centered on the observed neutral mass and scores only candidates whose
`nominal(peptide.mass − H2O)` falls in it ([match_engine.rs:87-110,127-136](../../crates/search/src/match_engine.rs)).
For a glycopeptide the observed precursor = `backbone + glycan`, so a window at the
precursor mass contains **no correct backbone peptide** (it's off by the glycan
mass, 1–3 kDa). You must search at the **backbone** mass. Two routes:

### Route (i) — Y-ion-determine the backbone mass first, then search at that mass

Read **Y0 (bare peptide)** / **Y1 (peptide + HexNAc)** from the stepped-HCD
spectrum, deduce the backbone neutral mass, open **one** candidate window there,
and run the ordinary search. The glycan = `precursor − backbone`, assigned to a
composition afterward.

- **Pro:** *one* search window per spectrum → **same cost as a normal search**;
  reuses `candidate_nominal_bounds` + the bucket index unchanged; the backbone b/y
  ions are scored by `score_psm` as-is.
- **Con:** depends entirely on **detecting Y0/Y1**, which can be weak/absent at low
  collision energy; a wrong backbone mass = a missed ID. This is the make-or-break
  prototype (§7). pGlyco3 mitigates with its glycan-first Y-complementary index +
  core-ion minimums; the MVP can start with a simpler "find the dominant Y-ladder
  spacing" heuristic.

### Route (ii) — search peptide + each composition (combinatorial)

For each glycan composition `g` in the DB, search a window at `precursor − mass(g)`.

- **Pro:** does not depend on detecting Y0/Y1; robust to weak glycan fragmentation;
  this is essentially what **Phase-0 glyco-lite does today** (each declared
  composition is a variable mod → its own precursor match).
- **Con:** **N windows per spectrum** (N = DB size, **~182**). With the
  GF mass-indexed DP this is **182× the candidate-scoring work per glyco scan** —
  and crucially, **the DP node cache is mass-indexed / peptide-agnostic**, so each
  backbone-mass window is a *different* set of nodes; the cache does not amortize
  across compositions. The neutral-loss work already hit this exact constraint: it
  found the DP "is mass-indexed/peptide-agnostic and the loss shift is
  peptide-specific, so loss scoring CANNOT flow through the cache — it is a
  **peptide-aware additive pass** instead"
  ([handoff doc](../plans/2026-06-13-andes-glyco-training-benchmark-handoff.md);
  [plan Task 7](../plans/2026-06-13-andes-neutral-loss-glyco-plan.md)). The same
  lesson applies: the glycan delta is per-composition, so the per-composition
  windows cannot share node DP state. Route (ii) is correct but **O(N) slower**.

### The GF mass-indexed-DP constraint, and why it favors Route (i)

The node-score DP keys on **nominal fragment mass**, independent of which peptide
produced it. That is what makes andes fast (the cache amortizes across candidates
at the *same* mass) and is exactly why:

- **Route (i)** is cheap: one backbone mass → one window → the DP cache works
  normally, and the glyco-specific evidence (Y-ladder, oxonium) is scored as an
  **additive pass** (the neutral-loss primitive's pattern, §6).
- **Route (ii)** is expensive: N backbone masses → N windows, no cache sharing.

**Recommended design:** **Route (i) primary** (Y-determine backbone → one window),
with **Route (ii) as a bounded fallback** for high-oxonium scans where Y0/Y1
detection fails — restrict the fallback to a small high-confidence composition
subset to cap the blow-up. This mirrors pGlyco3 (glycan-first candidate gen, but
capped at top-100) and MSFragger (mass-offset list) rather than a naive full cross
product. **Stepped-HCD is the enabler:** because oxonium + Y-ladder + backbone b/y
all appear in one scan, Route (i) has the Y-ladder it needs *and* the backbone b/y
for `score_psm` *in the same spectrum* — no MS3, no paired-scan logic.

## 6. What reuses the neutral-loss primitive (the key leverage)

The just-shipped neutral-loss feature (`feat/enzyme-support`) is the **fragment-
level glycan-scoring primitive** the glyco mode builds on. Concretely:

- **Y-ladder rungs are loss-shifted b/y ions.** A glycopeptide's Y-ion ladder
  (peptide + progressively trimmed glycan) and its glycan-retaining backbone ions
  are exactly what `predict_by_ions_with_losses` emits: for a residue carrying a
  glyco mod with `neutral_losses = [Hex 162.0528, Hex2 324.1056, …]`, it predicts
  one loss-shifted ion per declared loss at `mz_intact − L/z`, tagged
  `loss_class = glyco (1)` ([fragment_ions.rs:146-273](../../crates/scoring/src/scoring/fragment_ions.rs)).
- **They're scored against a trained, pooled, glyco-class table.** `score_psm`
  adds a **peptide-aware additive loss pass** that probes the model's pooled
  per-class glyco loss table at the shifted m/z, gated on activation +
  `has_loss_tables` + a declared loss mod ([psm_score.rs:286-340](../../crates/scoring/src/scoring/psm_score.rs)).
  This is **the additive-pass pattern** the glyco mode needs (§5): glyco evidence
  scored *additively* on top of the backbone b/y RawScore, not threaded through the
  mass-indexed DP cache.
- **Activation gating is already correct.** `predicts_neutral_losses()` returns
  `false` for ETD/electron-based methods ([activation.rs](../../crates/model/src/activation.rs)) —
  matching the glyco literature (HCD strips glycan → loss ions; ETD/EThcD preserves
  it → glycan-retaining c/z ions instead). The glyco MVP inherits this for free.
- **The model store + training already carry it.** The parquet store serializes
  `loss_class` (back-compat), the `glyco` model slug exists
  ([catalog.rs:95](../../crates/model-train/src/catalog.rs)), and `andes train` can
  learn the glyco loss-class table from a stepped-HCD corpus (handoff doc).

**What the neutral-loss primitive does NOT provide** (the glyco mode's net-new
work, restated): the **glycan search space** (composition DB), **Y0/Y1 backbone-
mass determination** (it scores Y-rungs given a *known* peptide; it does not
*solve* for the unknown backbone), **joint peptide/glycan mass solving**, the
**oxonium spectrum-level trigger**, and **glycan-level FDR**. The primitive is the
fragment scorer; the glyco mode is the search/identification/FDR layer around it.

## 7. Phased roadmap

- **Phase 0 — glyco-lite (SHIPPED, document the limits).** A single (or few)
  declared glycoform searched as variable mod(s) + `loss=`/`class=glyco`
  (neutral-loss milestone). Finds only the declared compositions; searches them
  combinatorially (Route ii at tiny N); **no glycan FDR; no Y-determination**.
  *Action:* document it as "known-glycoform mode," not a glyco search engine.
- **Phase 1 — oxonium gate + Y0/Y1 backbone-mass prototype (the make-or-break;
  *measurement only, no scoring change*).** (1) **oxonium classifier** — flag a
  scan as glyco from {204.087, 138.055, 366.140, 292.103} vs an intensity/k-of-n
  threshold; (2) **Y-ladder reader** — detect Y0/Y1 and deduce the backbone neutral
  mass. Validate on **PXD025455** (stepped-HCD, Q Exactive HF): does andes recover
  backbone masses that, when searched, hit the known glycopeptides? **Gate Phases
  2–4 on this.** Pure analysis tool; touches no scoring core.
- **Phase 2 — glycan composition DB + sequon-constrained backbone search
  (`--glyco` mode; medium).** (1) ship a curated **N-glycan composition list**
  (≈182) → masses + Y-ladder; (2) **Route (i)** backbone-mass search (one window at
  the Y-determined mass) reusing `candidate_gen` + `score_psm`, with an **N-X-S/T
  sequon filter**; (3) assign `glycan = precursor − backbone` to the nearest
  composition; (4) **Route (ii) bounded fallback** for Y-detection failures. Reuses
  the linear engine end-to-end.
- **Phase 3 — glycan scoring + combined evidence (reuses the neutral-loss
  primitive).** Score the **Y-ladder + glycan-retaining ions** via the existing
  per-class glyco loss table (§6); emit **oxonium count, Y-ladder match, glycan
  mass error** as **additive PIN feature columns**; let Percolator weight peptide
  vs glycan evidence (no hand-tuned `w`). Train the `glyco` model loss table on the
  PXD025455 corpus (handoff doc) — gated on the VM benchmark.
- **Phase 4 — glycan-level 2D FDR.** Reverse-sequence peptide decoys (existing) +
  a **FragPipe-style glycan decoy** (random intact-mass shift within tol + a random
  isotope error + each Y/oxonium fragment scrambled 1–20 Da); combine to a
  **glycopeptide-level FDR** (`FDR_GP`). New post-search layer; the largest net-new
  piece after Phase 1.
- **Phase 5 (reconsider whether to do at all) — O-glyco / multi-site
  localization.** No sequon, multiple glycans per peptide → O-Pair's graph-based
  site-localization combinatorics. A different, harder sub-problem; only after N-
  glyco MVP proves out. Likely its own effort.

## 8. Prototype & benchmark first (make-or-break, before committing to Phase 2)

1. **Oxonium classifier + Y0/Y1 backbone-mass reader** on **PXD025455**
   (stepped-HCD, Q Exactive HF — the target N-glyco dataset). Metric: of scans an
   external engine (pGlyco3/MSFragger published IDs) calls glycopeptides, what
   fraction does andes flag as glyco (oxonium sensitivity) and recover the correct
   backbone mass for (Y-determination accuracy)? **This is the make-or-break of the
   whole mode**, exactly as the doublet detector was for crosslinking.
2. **Route-(i) backbone search + composition assignment** on the same file, scored
   against pGlyco3/MSFragger IDs: at a controlled FDP, does andes's *existing*
   `score_psm` + a one-window backbone search + nearest-composition assignment
   recover known N-glycopeptides? Tells you whether the linear core suffices before
   investing in glycan scoring (Phase 3) and 2D FDR (Phase 4).

## 9. Hardest parts (ranked)

1. **Y0/Y1 backbone-mass determination** (Phase 1) — the component with no andes
   analogue; everything downstream depends on getting the backbone mass right.
2. **Glycan-level 2D FDR** (Phase 4) — a new post-search layer (glycan decoys +
   `FDR_GP` combination) beyond Percolator's peptide FDR.
3. **Glycan composition DB + Route-(ii) cost control** (Phase 2) — bounding the
   combinatorial fallback given the mass-indexed DP cannot amortize across
   compositions (§5).
4. **Robust oxonium gating** across collision energies / instruments (Phase 1).

## 10. Independence & licensing

Glyco fits the Apache-2.0 independence goal **if** andes uses its own algorithms
and permissively-licensed references only. **Read freely:** **O-Pair/MetaMorpheus
(MIT)** and **Glyco-Decipher (Apache-2.0)** — the latter matches andes's target
license exactly. **Papers only (do not touch internals):** pGlyco3,
MSFragger-Glyco, Byonic, StrucGP — binary-only, license-gated, or proprietary
(same posture as pLink2 for crosslinking). The general method (oxonium trigger +
Y-ion backbone mass + glycan composition DB + glyco FDR) is multiply, independently
published → clean-room implementable. **Guardrails:** avoid Protein Metrics'
**wildcard-modification-search** patent (US20100124785A1) — the MVP's *curated*
composition list steers clear of it — and Byonic's interactive-visualization GUI
patents. **Required:** a patent review of the specific chosen Y-determination /
glycan-decoy / 2D-FDR methods before implementation; no blocking patent on the
general method was found, but no exhaustive IP search was done.

## 11. Open questions / decisions before any implementation

- Is intact glyco strategically in scope for andes (alongside linear ID +
  independence + the neutral-loss milestone)?
- N-glyco only (MVP) acceptable as the product, or is O-glyco required (much
  larger)?
- `--glyco` mode in the main binary vs a separate `andes-glyco` crate for the
  glycan DB + glyco FDR (recommend: mode in the binary; DB + FDR in a crate)?
- Route (i) Y-determination as primary with bounded Route (ii) fallback, or
  composition-first (Route ii) for robustness at the cost of speed?
- Is **PXD025455** sufficient for the Phase-1 prototype *and* the Phase-3 training
  corpus, or is a second stepped-HCD dataset needed?
- Glycan composition list: adopt MSFragger's/pGlyco's published 182 N-glycan list
  (clean-room from the paper), or curate from a public glycan DB (GlyGen)?

## 12. References

**Engines / algorithms**
- pGlyco3 (Y-complementary index, 3-level FDR, sceHCD): https://www.nature.com/articles/s41592-021-01306-0 · pGlyco 2.0 (scoring/FDR foundation): https://www.nature.com/articles/s41467-017-00535-2
- MSFragger-Glyco (mass-offset open search, oxonium gate): https://www.nature.com/articles/s41592-020-0967-9 · MSFragger-Labile: https://www.mcponline.org/article/S1535-9476(23)00048-8/fulltext
- FragPipe glycan FDR (decoy-glycan method — the recommended target): https://pmc.ncbi.nlm.nih.gov/articles/PMC8933705
- Byonic glycoproteomics (peak filtering, Y0/Y1, wildcard): https://pmc.ncbi.nlm.nih.gov/articles/PMC8724605
- StrucGP (modular structural sequencing, dual-energy HCD): https://www.nature.com/articles/s41592-021-01209-0
- O-Pair / MetaMorpheus (MIT; graph localization, open search): https://www.nature.com/articles/s41592-020-00985-5 · https://github.com/smith-chem-wisc/MetaMorpheus
- Glyco-Decipher (Apache-2.0 reference): https://pmc.ncbi.nlm.nih.gov/articles/PMC8990002/ · https://github.com/DICP-1809/Glyco-Decipher

**Patents (guardrails)**
- Protein Metrics wildcard-modification search: https://patents.google.com/patent/US20100124785
- Protein Metrics interactive MS analysis GUI: https://patents.google.com/patent/US10546736B2/en

### andes code references
- single-peptide invariant `suffix = peptide_nominal − prefix` (reusable at backbone mass): `crates/scoring/src/scoring/psm_score.rs:234-311`
- neutral-loss primitive — prediction: `crates/scoring/src/scoring/fragment_ions.rs:146-273`; additive scoring pass: `crates/scoring/src/scoring/psm_score.rs:286-340`
- activation gate (ETD suppresses loss ions): `crates/model/src/activation.rs`
- precursor-mass-windowed candidate gen + bucket index (extend to backbone mass): `crates/search/src/match_engine.rs:87-110,127-136`, `crates/search/src/candidate_gen.rs:40-203`
- `Modification` scalar + `neutral_losses`/`loss_class` (extend to a `Glycan` concept): `crates/model/src/modification.rs:28-46`
- reverse-sequence decoys + Percolator peptide FDR (reuse; add glycan-level layer): `crates/search/src/decoy.rs:11-25`, `crates/output/src/pin.rs:11-20`
- `glyco` model slug (training target): `crates/model-train/src/catalog.rs:95`
- neutral-loss training/benchmark handoff (Phase-3 corpus): `docs/plans/2026-06-13-andes-glyco-training-benchmark-handoff.md`
