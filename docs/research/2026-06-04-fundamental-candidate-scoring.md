# Fundamental deterministic candidate-scoring improvements for SIMAS

**2026-06-04** — synthesis of the internal idea docs (`internal-docs/2026-06-04-simas-
relational-candidate-rescoring-architecture.md`, `…cluster-latent-psm-math.md`,
`…tmt-lfq-transformational-algorithms-research.md`, `…high-resolution-generating-
function-investigation.md`), the papers (`internal-docs/papers/`), and the
parity-analysis findings.

**Frame (the SIMAS identity):** SIMAS is a *deterministic candidate-generation-and-
scoring* engine for **one MS run**. Its scores should deterministically capture three
axes: **(1) spectrum explanation**, **(2) the candidate's relation to other candidates'
ranks** (competitive/listwise), **(3) experiment-level context**. ONNX/Prosit rescoring
and Percolator are *separate later deciders*. Everything below is **deterministic and
in-engine** unless flagged "needs data/models" or "needs multi-run."

---

## The validation guardrail (read first — it gates everything)

Reversed-decoy TDC + Percolator is **blind to coincidental real-DB targets** and is
postprocessor-artifact-prone: it refuted 4 chimeric levers, and rank-1 wide-window
*inflated +40% on Astral while keeping a HEALTHY decoy fraction* (`2026-05-29-rank-
stratified-fdr-bench.md`). **So no scoring change is "real" until it passes an
entrapment-FDP check.** Building the **entrapment-FDP ruler (FDRBench)** is the single
highest-leverage missing piece — it re-grounds *every* number, not just chimeric. **Do
this first; gate all of the below on it.**

---

## Axis 1 — Spectrum explanation

### 1A. GF-free per-spectrum significance (OMSSA Poisson tail) — **fundamental**
The GF is gone; the dominant ranking-flip axis was *significance*, not RawScore
(72.8% `spec_e_swap_only`, `2026-05-25-spece-tail-exploration.md`). The closed-form,
model-free replacement: OMSSA's Poisson tail. With fragment tol `t`, theoretical
fragments `h`, observed peaks `v`, precursor mass `m`: `μ = 2thv/m`; per-spectrum
significance of `y` matched ions `= Σ_{x<y} P(x,μ)`; E-value `= N·[1−(…)ᴺ]`. **Naturally
per-spectrum-calibrated** (noisier spectrum → larger μ → higher bar). Pure-deterministic,
no model, no decoy DB. This is *the* deterministic significance that fills the low-res
hole the GF left.

### 1B. Fold high-res accuracy INTO the score (not only Percolator features)
The GF forced ~0.5 Da matching (integer mass axis) and wasted the instrument's 20 ppm —
accuracy currently lives only in Percolator (`MeanErrorTop7`/`StdevErrorTop7`). GF-free,
SIMAS can match at true ppm and add a **mass-error-consistency** term + an **averagine
mass-defect-band noise filter** (a peak whose fractional mass is outside the expected band
for its nominal mass is unlikely real → down-weight). Deterministic, in-engine.

### 1C. Residual-spectrum + precursor-existence scoring
Score what's **unexplained** after the candidate (clean residual = single peptide;
structured = chimeric). Plus the **targeted-XIC + KL-divergence isotope gate** (does this
precursor exist with the right isotope envelope — DDA+). SIMAS reads `.raw`/`.d` natively,
so the XICs are free. Deterministic.

### 1D. Cheap local-null calibrators (zero-cost, add regardless)
±75-bin **shifted-background subtraction** (Tide XCorr: correct alignment must beat random
offsets) and the **consecutive-ion-run bonus** `R` (contiguous b/y ladders are strong).

---

## Axis 2 — Relation to other candidates (the SIMAS differentiator)

### 2A. DeltaRawScore — **ship it now** *(already a proven win)*
`RawScore(best) − RawScore(2nd-best distinct)`. Benched **+129 PXD / +12 TMT / +104 Astral
@1% FDR, zero wall cost, no regression** (`feat/delta-raw-score` `bea5d697`), shelved only
by the old simultaneous-gate. Pure additive competitive feature; the realized orthogonal-
dimension lever. **Lowest-risk fundamental win available.**

### 2B. Top-K listwise candidate-set features (Prototype A — the cheap bridge)
Emit top-K per scan internally; compute: delta(cal), **rank entropy** / softmax over
candidates, **ambiguity group size** (#within ε of best), **top-1 vs top-2 shared-fragment
fraction**, **unique matched-ion count**, candidate-set p-value. Makes "how much better than
the alternatives" *native* rather than bolted-on. Deterministic.

### 2C. Tailor refinements + Šidák candidate-count correction
We have a basic Tailor. The missing pieces (Tailor paper): **rank-3 position floor**
(`Q₁₀₀ = s_{min(3,⌊N/100⌋)}` — never use rank-1/2 as the reference), **N≥30 null-padding**
from outside the precursor window (null-only, not assignable), **positive-shift guard** for
any negative sub-score. Then **`p = 1 − (1−p′)ᶜ`** (cascaded search) to make per-spectrum
significance comparable across scans with different candidate counts. Deterministic.

### 2D. Greedy residual rescoring for co-isolated candidates (the chimeric win, formalized)
Fragment "theft" is **real** (Astral mean 0.37 overlap, 38% of scans ≥0.5; bimodal,
`2026-05-28-chimeric-fragment-overlap-diagnostic.md`). DDA+ pattern: rank candidates, top
claims its peaks, **mask explained peaks, rescore the rest on the residual**. Turns
independent scores into **unique-evidence** competitive scores — the deterministic core of
the validated +80–101% Astral chimeric gain. **Caveat:** the ~28% coincidental-overlap mode
means this alone does NOT fix FDR inflation → must be gated by 2C + the entrapment ruler.

---

## Axis 3 — Experiment-level (single-run-applicable, deterministic)

### 3A. Cluster-latent evidence pooling — **the most distinctive contribution**
Your own `cluster-latent-psm-math.md`: repeated/similar spectra are noisy observations of
ONE latent peptide `z`. Deterministic kernel:
`cluster_score(p) = Σ_i calibrated_score(s_i,p) + α·coverage − β·ambiguity − γ·impurity`,
gated by `transfer_weight = P(cluster pure)`, with **leave-one-out back-transfer** (a
spectrum can't be in its own consensus) and **cluster-level TDC** (one cluster = one
discovery, not N). "Weak + weak + same-latent constraint = strong." Needs only SIMAS's
calibrated scores + cheap cosine clustering + candidate unions — **no NN, no extra data**.
This is the "experiment graph" the architecture doc centers on, in its leanest deterministic
form, and it's the most novel-vs-MSFragger/MS-GF+/Percolator idea here.

### 3B. Stratum/relevance-aware calibration (cascaded search)
Calibrate a candidate against **its own relevance stratum's null** (unmod vs open-mod,
tryptic vs semi/non-tryptic), `k≈20` guard, "search relevant-first, commit, remove." A
1-mod candidate shouldn't compete in a flat null. Deterministic, structural.

### 3C. MS1-precursor-dominance feature (Quandenser, single-run part)
Is the candidate's precursor the **dominant co-isolated MS1 species** in its isolation
window (`s ≥ max/2`)? A within-run, deterministic competitive feature for chimeric scans.

### 3D. Experiment-conditioned fragmentation tables *(needs the two-pass + model-train)*
Search → learn THIS run's per-charge/segment rank distributions + neutral-loss + reporter
behavior (count tables w/ smoothing, **not** a NN) → rescore. Deterministic; reuses the
`model-train` substrate; ties directly into Phase 3 (own models). Pillars 1+3.

---

## The deep root cause (the low-res 6%) — honest placement
Both the Astral score-psm trace and the TMT diagnosis independently localize it: **Rust
under-scores the true peptide's RawScore on CID so badly it leaves the top-10 on 95% of
flips** — the per-rank log-prob **table application (H3) + ion-type enumeration (H1)** on
CID fragmentation (H2/peak-rank is proven NULL, `2026-05-28-phase2-peak-rank-parity.md`).
This is a **node-scoring / model** problem → **Phase 3 (own models)** is the real fix, not a
new feature. **Caution (Rule-2):** historically only *top-1-restoring* node-score changes
gained PSMs; modifying the node score regresses Percolator (n=12). So touch the node score
only via retrained tables (Phase 3), gated on entrapment.

---

## Recommended order (deterministic, low→high risk)
0. **Entrapment-FDP ruler (FDRBench)** — the guardrail; do first.
1. **DeltaRawScore** — ship the proven win (2A).
2. **Tailor refinements + Šidák count correction** (2C) and **OMSSA Poisson significance**
   (1A) — the GF-free deterministic significance the engine now lacks.
3. **Top-K listwise features** (2B) — the competitive axis, cheap.
4. **Cluster-latent pooling** (3A) — the distinctive experiment-axis contribution; validate
   hard with LOO + cluster-TDC + entrapment.
5. **Residual / unique-evidence + MS1-dominance** (1C/2D/3C) — the chimeric win made
   principled.
6. **Experiment-conditioned tables** (3D) — folds into Phase 3.

NOT SIMAS's job (keep in the separate rescoring layer): transformer/contrastive rerankers,
Prosit/MS2PIP predicted-intensity, MSBooster/MS2Rescore, GLEAMS embeddings, de-novo engines,
all quant/MBR. The cross-run recurrence term and full peptide-centric detection need
multi-run data and are out of single-run scope until a cohort mode exists.
