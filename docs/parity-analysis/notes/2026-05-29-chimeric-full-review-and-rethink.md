# Chimeric search — full review, re-think, and decomposed plan

**Date:** 2026-05-29
**Branch:** `feat/chimeric-dda-plus`
**Purpose:** Step all the way back. Trace everything we tried, compare against how
the established open-source engines actually do it (deep-research survey, cited
below), and re-decompose into the smallest viable steps.

---

## Part 1 — Internal trace: what we tried, and the one root cause

We pursued chimeric (wide-isolation-window, multi-peptide-per-scan) search to
recover co-fragmented IDs. Four levers, all bench-refuted against
reversed-decoy TDC + Percolator:

| Lever | What | Result | Note |
|---|---|---|---|
| (a) Phase 1 | wide-window multi-emission | Astral +94–97% inflation | `chimeric-phase1-bench` |
| (b) Phase 2 | MS1 isotope-KL **soft** Percolator feature | insufficient (+89% after hard-filter test) | `chimeric-phase2-bench` |
| (c) Phase 3 | greedy fragment competition + **residual SpecEValue re-score** | Astral canary +111% | `chimeric-phase3-bench-canary-fails` |
| (d) | **rank-stratified separate FDR** | rank-1 alone +40% (Astral 51,579) | `rank-stratified-fdr-bench` |

**Root cause we isolated (d):** wide-window search inflates IDs *even at rank-1*
(one best peptide per scan) with a **healthy target-decoy ratio** (Astral rank-1
decoy_frac 0.39, *better* than off's 0.52) yet an implausible count (51k on a
15-min run, truth ≈36k). That is **coincidental real-DB-sequence targets winning
by chance** — a *target* with no decoy twin, so reversed-decoy TDC scores it ~99%
correct. The inflation is in the wide-window *selection*, and TDC is structurally
blind to it.

We concluded "no post-hoc FDR re-stratification can fix this." **The survey below
shows that conclusion is right — and that we were also measuring with the wrong
ruler.**

---

## Part 2 — External survey: what the established engines actually do

Deep-research (21 primary sources, 25 claims adversarially verified, 21 confirmed).
Full report: workflow `wf_4fb9718a-59c`. Key, verified findings:

### 2.1 MSFragger-DDA+ (Nat Commun 2025) — the flagship — uses OUR refuted levers
- **Full-isolation-window search against ALL peptides in the window**, fragment-ion
  index + hyperscore = our lever (a), verbatim. [3-0]
- **MS1 XIC + Kullback-Leibler isotope filter is POST-search and SOFT** ("performs
  a full isolation window search and then detects and utilizes MS1 precursor
  peaks"; "PSMs with low-quality XICs are discarded"). **No hard MS1 admission
  gate** = our lever (b), verbatim. [3-0]
- **Greedy shared-fragment competition** (remove top PSM's peaks, drop it,
  recompute scores, re-rank; Yu et al. 2023) = our lever (c), verbatim — only the
  score function differs (hyperscore vs RawScore+GF SpecE). [3-0]
- **FDR = the STANDARD FragPipe chain** (MSBooster → Percolator → ProteinProphet →
  Philosopher), **no chimeric-specific FDR**, rank≥2 rides the same postprocessor
  as rank-1 — i.e. it does NOT do our lever (d) at all. [3-0]

**So the flagship engine does (a)+(b)+(c) and explicitly NOT (d) — and never claims
those mechanisms fix coincidental-target inflation.** Our levers were mainstream;
our error was expecting them to control FDR by themselves.

### 2.2 How DDA+ actually earns trust: ENTRAPMENT, not the search
DDA+ proves trustworthiness **empirically with an entrapment database** (each
target protein shuffled into one entrapment protein, C-termini fixed) measuring
**false-discovery proportion (FDP)** — reporting ~57% more peptide sequences at
similarly-low/lower FDP than peers (Fig 2d). [3-0; the "genuine not inflation"
sub-claim was 2-1 — developers' own benchmark, co-fragmented peptides have lower
quality.]

### 2.3 Reversed-decoy TDC is known to be blind here — entrapment is the field check
Wen/Freestone/Riffle/Noble, **Nat Methods 2025**: "many analysis pipelines
implement variants of the procedure that potentially fail to control the FDR";
entrapment "expand[s] the tool's input … so that its search space includes
verifiably false entrapment discoveries … hidden from the tool itself." DIA/
protein-level inflation up to ~48% observed. **FDRBench (Noble-Lab, github)** builds
shuffled/foreign-species entrapment DBs and estimates FDP via lower-bound,
combined, and **paired** (per-peptide) methods. [3-0] **This is exactly the tool
for our symptom** — it exposes coincidental targets TDC cannot.

### 2.4 The bigger reframe: the postprocessor may be the real culprit
Freestone/Noble/Keich, **JPR 2024**: the wide/open-window FDR-control failures
"uncovered a potential general problem … in the machine learning postprocessors
Percolator and PeptideProphet," not in TDC itself (Percolator cross-validation can
indirectly peek at labels via near-identical PSMs split across folds). [2-1, hedged]
**This dovetails with our own memory note**: TDC FDR (no Percolator) is Java 24,561
vs Rust 8,506, but Percolator boosts Rust +270% via correlated features. **Meaning:
we have been measuring chimeric trustworthiness with reversed-decoy-TDC + Percolator
— precisely the instrument the literature says is blind to coincidental targets AND
prone to postprocessor artifacts. We refuted four levers with a broken ruler.**

### 2.5 Peer engines (Sage, etc.)
- **Sage** (Rust): `wide_window` and `chimera` are two independent flags (both
  default off); chimera emits multiple PSMs/scan; **no chimeric-specific FDR**.
  Uses **picked-PEPTIDE FDR** (reverses tryptic peptides, not proteins) + hierarchical
  spectrum/peptide/protein q-values. [3-0] Picked-peptide FDR is a worthwhile
  refinement but is **not** evidence Sage solves coincidental-target inflation.
- **The actual discriminator (inferred, not isolated by any single source):**
  deep-learning **spectral-prediction rescoring** — MSBooster/Prosit features
  (predicted-vs-observed fragment-intensity correlation + delta-RT) — is what gives
  the postprocessor real power over coincidental targets in practice. Caveat: no
  source isolated rescoring as THE causal FDR-trust factor; it's a strong synthesis.
- **"ProSE"**: not found as a real search/rescoring tool in any verified source —
  likely a confusion with **Prosit** (the DL fragment-intensity predictor). Treat as
  Prosit unless you meant something specific.

---

## Part 3 — The re-think (what this changes)

1. **We can't see truth with our current ruler.** Reversed-decoy TDC + Percolator
   is blind to coincidental targets and artifact-prone. Every chimeric verdict (and
   arguably the whole Java-vs-Rust PSM-gain comparison) needs an **entrapment-FDP**
   ground truth to be trusted. This is the single highest-leverage missing piece —
   and it is **independent of chimeric**: it would also tell us the true FDP of our
   baseline Rust-vs-Java numbers (which our own memory flags as Percolator-inflated).
2. **The trust mechanism is downstream, not in-search.** The field does NOT add a
   chimeric-specific FDR. It (i) emits cheaply, (ii) discriminates with **DL spectral
   rescoring**, (iii) **validates with entrapment**. Our missing piece vs DDA+ is the
   DL rescoring (MSBooster-style), not a cleverer in-search FDR.
3. **Lever (d) was correctly refuted** and is confirmed a dead-end: the flagship
   doesn't do rank-stratified FDR; coincidental-target inflation is real and not
   fixable by re-stratification.
4. **Gate reality unchanged.** Even fully built, chimeric is no-op on Astral /
   net-negative on TMT → does not clear [[merge-gate-beat-java]]. BUT the two
   building blocks below have value **beyond chimeric**.

---

## Part 4 — Decomposed plan (smallest viable steps, ranked by evidence ÷ effort)

### Step 1 — Entrapment-FDP harness (FDRBench). *Do this first; project-wide value.*
Build/adopt an entrapment-DB + FDP measurement (FDRBench: shuffled or
foreign-species entrapment; paired per-peptide method). Run it on what we already
have: baseline narrow Rust, Java, and the chimeric PINs.
- **Evidence:** highest — Nat Methods 2025 + FDRBench are the field standard; directly
  targets our exact symptom.
- **Effort:** low–moderate (tooling, no engine change; entrapment DB + re-search +
  FDP script).
- **Payoff:** a *trustworthy ruler*. Confirms (or overturns) the coincidental-target
  story, and — crucially — re-grounds the **whole** Rust-vs-Java comparison, not just
  chimeric. Likely reveals whether our Percolator counts are real.
- **Smallest version:** entrapment DB + re-run + FDP on the existing narrow baseline
  first (chimeric off) — answers "is our normal FDR even honest?" before touching
  chimeric.

### Step 2 — DL spectral-prediction rescoring features (MSBooster/Prosit-style). *The real discriminator.*
Add predicted-fragment-intensity correlation + delta-RT as **additive** PIN features
(predict via an external model; we already have a lean-external-handoff pattern from
Phase-2). This is what DDA+ relies on to make wide-window trustworthy.
- **Evidence:** high (mainstream; the inferred causal discriminator).
- **Effort:** moderate–high (model integration / external handoff), but **additive**
  and **gate-relevant beyond chimeric** — could lift PXD/TMT narrow-search IDs too.
- **Gate by Step 1:** only trust its gains via entrapment FDP.

### Step 3 — Picked-peptide hierarchical FDR (Sage-style). *Cheap refinement.*
Reverse tryptic peptides (not proteins) for picked-peptide FDR + hierarchical
spectrum/peptide q-values.
- **Evidence:** medium (refinement, not a coincidental-target fix).
- **Effort:** low–moderate. Worth folding in once Step 1 exists to measure it.

### Explicitly NOT doing
- Another in-search chimeric FDR scheme (lever d class) — refuted internally and
  absent from the flagship.
- A hard MS1 admission gate — the flagship deliberately doesn't; MS1-KL stays a
  soft post-filter (our lever b already matches state of the art).

---

## Recommendation

**Pivot the immediate work to Step 1 (entrapment-FDP harness), decoupled from
chimeric**, because it re-grounds every PSM number we report — including the
gate-blocking PXD/TMT comparison — with a ruler the field trusts and ours isn't.
Chimeric itself stays shelved (all four levers refuted; gate-irrelevant). If Step 1
shows our baseline FDR is honest and we still want chimeric sensitivity later,
Step 2 (DL rescoring) is the evidence-backed way to make it trustworthy — and it
helps the narrow search too.

## Caveats (from the survey)
- DDA+ entrapment results are the developers' own (no independent replication);
  co-fragmented peptides have lower quality (the "genuine IDs" claim was 2-1).
- "Postprocessor is the culprit, not TDC" is a strong working hypothesis, not
  settled (2-1, hedged in-source).
- "DL rescoring is THE discriminator" is synthesis, not a directly isolated result —
  Step 1 + a rescoring-on/off entrapment-FDP test would settle it.
- "ProSE" unconfirmed as a tool (likely Prosit).

## Sources (verified, primary)
- MSFragger-DDA+: nature.com/articles/s41467-025-58728-z · pmc PMC11978857 · biorxiv 2024.10.12.618041
- Entrapment / FDR validity: Nat Methods 2025 s41592-025-02719-x · pmc PMC12240826 · github Noble-Lab/FDRBench · JPR 2024 10.1021/acs.jproteome.3c00902
- Sage: github lazear/sage/blob/master/DOCS.md
