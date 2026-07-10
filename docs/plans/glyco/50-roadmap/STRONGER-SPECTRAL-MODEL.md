# andes-glyco stronger fragment/spectral model — design (2026-07-09)

> Goal: make the WEAK backbones FDR-passable. The @1% ceiling (~319 of 523 the reference engine
> truth) is bounded by EVIDENCE STRENGTH, confirmed flat/negative across six levers
> (learned selector, GBDT-feature, retention, expansion, list-widening, calibration).
> No scoring/selection/feature trick moves it — the true backbones that fail FDR lack
> the discriminative fragment signal to separate from decoys. This model attacks THAT.

## Two non-negotiable design constraints (learned the hard way)

1. **Optimize target/decoy SEPARATION at @1%, not top-1 correctness.** The learned
   selector recovered +11 top-1 but LOST @1% (120 vs 136) because it picked
   weak-evidence winners that fail FDR. Every component here is validated by held-out /
   cross-dataset **@1% yield**, never top-1.
2. **The output must be EVIDENCE-STRONG, FDR-aligned.** It augments/replaces the base
   b/y score and is consumed by Percolator; it does not merely re-select winners.

## Why the current evidence is thin here (root causes)

- **Regime mismatch.** The b/y rank + intensity models are trained on STANDARD
  peptides. Stepped-HCD GLYCOPEPTIDE backbones fragment differently: the glycan carries
  charge/energy, so bare-backbone b/y are sparse and their intensity distribution
  differs. andes scores these spectra with the wrong prior.
- **Untapped evidence (the key gap).** andes scores only the BARE backbone b/y
  (deglycosylated) plus the glycan Y-ladder (Y0..Yn on the intact backbone). It does NOT
  score the **peptide-fragment + partial-glycan** ions — b_i/y_i still bearing the
  innermost glycan (b_i+HexNAc, y_j+HexNAc, +HexNAc2, +core-fucose). In stepped-HCD the
  glycan sheds STEPWISE, so these are abundant — a whole second ladder of
  peptide-SEQUENCE evidence currently thrown away.
- **Linear FDR model.** Percolator (linear SVM) can't capture nonlinear spectral-match
  structure; a learned score can pre-compute it.

## The model — four components, cheapest/highest-ROI first

### A. Glyco-regime base model (retrain, not rebuild) — do first
Retrain the b/y rank + fragment-intensity models on GLYCOPEPTIDE-BACKBONE spectra
(stepped-HCD, deglycosylated b/y), labels = the reference engine/a glyco search engine confident backbone IDs.
andes ALREADY has the entry point: `andes train-from-search --labels <scan,peptide,charge>`
(SP-B glyco training path). Regime-matched intensities → the SAME RankScore/strong_score
features become sharper on exactly the weak spectra. Cheap, existing infra, low risk.

### B. Partial-glycan b/y ladder (the novel evidence multiplier) — the bet
For each backbone cleavage i, score not just b_i/y_i (bare) but the retained-glycan
series: b_i+HexNAc, b_i+HexNAc2, b_i+HexNAc2Hex, y_j+HexNAc, … (peptide fragment +
innermost core glycan). Predict their intensities (regime model) and match. This
MULTIPLIES peptide-sequence coverage precisely where bare b/y is sparse — directly
attacking the thin-evidence wall. Mechanism is physical (stepwise glycan shedding in
stepped-HCD) and currently unexploited in andes' backbone score. Emitted as new
evidence features (matched partial-glycan b/y count/intensity/ladder-length) AND folded
into the base score.

### C. Full-spectrum coherence + complementary evidence
Add complementary b/y pairs (b_i + y_{n-i} = backbone+H2O), internal fragments,
immonium ions, and water/ammonia neutral losses; jointly score observed vs the FULL
predicted glycopeptide fragmentation (peptide b/y ⊕ partial-glycan b/y ⊕ Y-ladder ⊕
oxonium). Every extra matched ion is FDR-discriminating evidence for sparse spectra.

### D. FDR-aware deep spectral scorer (ambitious tier) — only if A–C plateau
A small neural spectral-match model: (peaks, backbone sequence, glycan) → one match
score, trained CONTRASTIVELY (true backbone vs reversed-decoy backbone vs same-scan
wrong-competitor) to MAXIMIZE target/decoy margin. Captures patterns the GBDT + hand
features miss. High cost, uncertain; last resort. Keep native-Rust-inferable (own the
stack; the GBDT engine is the fallback if a NN is too heavy).

## Training data & labels

- Positives: the reference engine/a glyco search engine confident backbone IDs, **cross-dataset** (Fc3_r1 +
  Fc5_r2 + others) so the model generalizes, not memorizes one run.
- Decoys: reversed-sequence backbones (for the contrastive/FDR objective).
- **Circularity caveat:** truth = the reference engine's calls, so the metric caps at 523. To
  truly EXCEED the reference engine, validate on ORTHOGONAL truth (synthetic/entrapment glyco
  standards), not the reference engine overlap.

## Integration (FDR-aligned)

- The model sharpens the BASE b/y score (regime model) and adds evidence features
  (partial-glycan ladder, complementary, coherence) — all EVIDENCE-STRONG so they
  correlate with FDR survival. gp fusion still selects (it's FDR-aligned); this makes
  the evidence it selects on stronger.
- Emit new features as ADDITIVE PIN columns for Percolator (which now actually gets its
  calibration features too — the dead-feature fix landed this session).
- **Validate every step at @1% target/decoy separation (held-out + cross-dataset).**

## Build order & honest expectation

1. **A (regime retrain)** — cheapest, existing CLI; test @1%. Expect a modest but real
   lift if regime mismatch is a factor (it was for low-res; likely here too).
2. **B (partial-glycan ladder)** — HIGHEST upside, offline-testable first: on the weak
   (present-but-fail-FDR) backbones, does adding partial-glycan ion matches widen the
   true-vs-decoy score gap? If yes, this is the breakthrough. Build the scorer, wire the
   features, validate @1%.
3. **C (complementary/internal)** — incremental hardening.
4. **D (neural)** — only if A–C stall; biggest effort, own-stack constraint.

**Ceiling reality:** A–C can plausibly move weak-but-present backbones (≈399–423 present
→ more passing FDR) toward the low-400s @1%. The **45 non-decomposable glycan-gap**
backbones remain unreachable without modified-glycan DB work (sulfate/phospho) — or they
are the reference engine artifacts andes is correct to miss. Beating 523 outright needs BOTH: this
model (evidence) AND the glycan-DB (coverage) AND orthogonal truth (to prove it).

**One-line thesis:** the peptide evidence isn't missing from the spectrum — it's on the
**partial-glycan b/y ions andes never scores**. Model those (B) on a glyco-regime base
(A), score the whole spectrum coherently (C), and the weak backbones gain the FDR
margin that six scoring/selection tricks could not manufacture.
