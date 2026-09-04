# andes Unique-Scoring Campaign Plan (pre-implementation)

**Date:** 2026-06-29
**Status:** synthesized design campaign — pending adversarial review (Codex) + independent judge, then implementation plan.
**Inputs:** the RS³ design spec ([2026-06-29-rs3-spectral-significance-design.md](2026-06-29-rs3-spectral-significance-design.md)) + three parallel explorer agents (candidate generation; spectrum-refinement scoring; de novo / cross-family inspiration) and their literature sub-agents.

> Objective (on record): own data / no patents / beat ALL tools on PSMs at 1% **true entrapment-FDP** / max speed. FDR is **Percolator-only**. Safe-change rule: **additive features help; modifying the matched-peak set / existing score regresses Percolator**. Thesis: **IDs ≈ coverage × separation** — discriminators only pay when candidate space is expanded.

---

## 0. The organizing idea

RS³ computes a per-spectrum significance of andes's own score over a **score landscape** `g(m) = prefix_node(m) + suffix_node(M−m)` evaluated densely over the mass axis. Two consequences shape the whole campaign:

1. **Anything additive into the node/matched-ion LLR is inherited by `g(m)` for free, symmetrically for targets and decoys, without mutating the peak set** — the safe corner of this codebase. All the recommended *scoring* refinements are of this form.
2. **`g(m)` is a node-product significance, so deleting/moving a real peak is catastrophic** (it zeros a node). This rules *out* the destructive refinements (top-K deletion, iterative peak cleaning, aggressive deisotoping) that other engines use.
3. RS³ is a **separation** lever. Separation converts to IDs only when **coverage** expands. So the campaign pairs RS³ + additive node features (separation) with a candidate-space expansion (coverage), or the scoring gains stay flat (the documented "discriminator flat on closed search" result).

The campaign is therefore three pillars: **Separation**, **Coverage**, **Discipline**.

---

## Pillar A — Separation: RS³ + additive node features (the core novelty)

All additive, own-data, patent-free, feeding `g(m)` without touching the peak set.

| # | Idea | Source | Why it's net-positive | Regime | Risk |
|---|---|---|---|---|---|
| A0 | **RS³** renewal-saddlepoint spectral significance | RS³ spec | Patent-free per-spectrum calibration; the centerpiece. Closes the low-res null-heterogeneity gap. | both (low-res workhorse) | independence approx — gate 0 prototype |
| A1 | **Complementary-ion node reinforcement** — bonus to `g(m)` when prefix(m) AND suffix(M−m) both present | refine #1 | Coincident prefix+suffix collapses the low-res noise floor of `g(m)` (two independent 0.5 Da chance matches multiply). Reuses existing `has_complement`. **Most RS³-synergistic, lowest effort.** | low-res ★ | low |
| A2 | **Predicted-intensity node weighting** from andes's OWN frag-intensity GBDT (asymmetric `explained`/`chance_llr` form, NOT cosine) | refine #2 / denovo #1 | The literature's #1 rescoring lever, realized as RS³'s flagship feature; Apache-clean (no Prosit weights). | both (esp. low quality) | collinearity w/ existing strong battery |
| A3 | **Additive local noise-floor S/N feature** (per-peak, into the peak GBDT) | refine #3 | The canonical low-res denoiser, reframed *additively* so it generalizes the TMT-only windowed top-K to UPS1/LFQ low-res **without deleting peaks**. | low-res | low |
| A4 | **`delta_to_second` competition margin** PIN feature | denovo #5 | Near-free separation signal (SpectraST/Comet "score gap"); complements RS³ significance with separation. | both | very low |
| A5 | **Own-trained spectral-angle / Spec2Vec features** | denovo #1/#4 | Orthogonal "forward-model physics fit" axis; train own predictor (MS²PIP is Apache + has an ion-trap CID model as reference). | high-res ★ | collinearity; train cost |
| A6 | **Casanovo-style de novo log-prob** as a re-ranking PIN feature / tag prefilter | denovo #3 | A third orthogonal "learned sequence prior" axis. Casanovo is Apache-2.0. Opt-in (transformer cost vs max-speed). | high-res | speed; opt-in only |

**Three orthogonal evidence axes for Percolator** (the design thesis): RS³ = *statistical significance* (low-res) · spectral-angle/Spec2Vec = *forward-model physics* (high-res) · de novo log-prob = *learned sequence prior*. Orthogonality is what makes the scoring uniquely strong; additivity is what keeps Percolator from regressing.

**Explicitly dropped** (RS³ fragility / FDR-leak / redundant): destructive iterative peak-cleaning; aggressive deisotoping beyond the existing bounded z-collapse; raw sqrt/rank intensity transforms (already captured by the rank model). Charge-deconvolution beyond z≤3 is **deferred** behind a hard entrapment gate (it mutates the peak set + RS³'s null density).

---

## Pillar B — Coverage: expand the space so separation pays

| # | Idea | Source | Why | Patent | Risk |
|---|---|---|---|---|---|
| B1 | **Semi-tryptic via streaming external-sort index build** (revive the abandoned lever) | candgen #1 | The ntt=1 enumeration **already exists & is test-guarded**; only `build_base_peptide_index`'s in-RAM collect+sort OOMs. Replace with spill+merge external sort; the serve path is already out-of-core. **Biggest free coverage multiplier for RS³** (10–40× candidates/spectrum → needs the calibrated null). | clean | medium (external sort + group-FDR) |
| B2 | **Fragment-ion inverted index** (Comet-FI Apache / Sage MIT designs, NOT MSFragger code) | candgen #2 | Speed substrate that makes B1 + mass-offset + deeper chimeric affordable under expansion. No direct ID gain; enables the rest. | clean (build from open) | high effort |
| B3 | **Chimeric depth sweep** (existing two-pass; `ANDES_CHIMERIC_OVERLAP` diagnostic) | candgen #3 | Near-zero code cost; high payoff on dense high-res (Astral). RS³ + chimeric pair naturally (explained vs residual peak set). | caution (vendor IP near deconv) | FDR (fragment theft) |
| B4 | **Mass-offset / restricted-open** (curated Δ list, not unrestricted open) | candgen #4 | PTM-rich / immunopeptidomics coverage; needs B2 first. | clean (published) | FDR per-offset group |

**Dropped:** de novo sequence-tag prefilter (dominated by B2 for andes's objective; needs a whole de novo engine).

---

## Pillar C — Discipline (non-negotiable, every item)

- **Gate 0 for RS³**: numerical prototype vs a brute-force decoy-peptide null before any engine wiring.
- **One variable at a time**, provenance-stamped (binary commit + model SHA + data SHA).
- **1% true entrapment-FDP**, not reported FDR. **Percolator only.**
- **Astral gate first** (keep only if ≥ flat where andes already leads), then prove the gain on **low-res UPS1** (the target), confirm no TMT regression.
- **Expanded-space items (B1/B4) require group/subset FDR** (emit the group column; Percolator post-process) — the expansion inflates the composite null otherwise.
- **Additive-only** for all Pillar-A scoring features; peak-set mutations (B-side excluded) are gated hard.

---

## Sequencing (prove cheap → scale)

1. **A0 RS³ (gate 0 prototype → additive PIN features)** + **A1 complementary-ion node bonus** + **A4 delta_to_second**. All additive, low effort, directly attack low-res. Astral-gate then UPS1.
2. **A3 local noise-floor feature** + **A2 predicted-intensity node weighting** (own GBDT). Retrain peak/frag models as needed.
3. **B1 semi-tryptic external-sort** — the coverage multiplier; re-A/B the Pillar-A features on the expanded space (this is where separation finally pays).
4. **B2 fragment-ion index** — speed substrate; then **B3 chimeric sweep**, **B4 mass-offset**.
5. **A5 spectral-angle / A6 de novo log-prob** — the bigger learned axes, last, opt-in for A6.

Rationale: Pillar A is the unique IP (a calibrated, own-data, patent-free scorer with three orthogonal axes). Pillar B is what makes Pillar A convert to PSMs. The first milestone (step 1) is cheap and falsifiable; if RS³ + complementary-ion don't move low-res UPS1 at honest FDP, we learn that before investing in B.

---

## Licensing bottom line

Every *method* here is patent-free (MS-GF+ generating function avoided by construction; XCorr patent expired 2014; complementary-ion physics, external sort, fragment indexing, spectral angle all unencumbered). The only traps are specific *artifacts*: **do not ship/call Prosit weights (CC-BY-NC) or Koina; pDeep has no license**. Train andes's own predictors (Apache). Reference designs: Comet (Apache), Sage (MIT), MS²PIP/Casanovo/AlphaPeptDeep/Spec2Vec (Apache).

---

---

## CAMPAIGN RULING (adversarial review + independent judge, 2026-06-29) — supersedes §Sequencing

A Codex adversarial review + an independent judge (every file:line re-verified against `crates/`)
ruled the campaign. Outcome: the organizing thesis (additive separation + expanded coverage) is
sound, but the first milestone as written was the part most likely to fail, and several Pillar-A
features are already present, redundant, or unsafe score mutations. RS³ is the genuine novelty and
is salvageable — but only after a reformulation.

### RS³ reformulation (fixes the two CRITICAL findings together)
Do **not** implement the analytic independent-bin CGF on the `g(m)` surrogate. Instead:
1. **Calibrate the REAL emitted score** (`pin_score = score_psm() + cleavage_credit`, integer node
   sum + loss term), reusing the production `score_psm` path — not the float `g(m)`.
2. **Replace the analytic null with a per-spectrum DECOY-CALIBRATED EMPIRICAL null**: score a fixed
   budget of mass-feasible decoy peptides (real AA random walks conditioned to hit M — these obey
   cleavage-site exclusion by construction, so the independence error evaporates) against the
   spectrum; take the tail. This promotes the spec's Gate-0 brute-force from a test oracle to the
   method itself. The renewal `ρ(m)`/`u(m)` may stay as an importance-sampling proposal for cheap
   tail-focused draws; the p-value comes from the empirical tail (optionally a Lugannani–Rice fit).
3. **Patent posture:** an empirical Monte-Carlo tail is *further* from US 8,639,447 (it builds no
   score-generating object). Still requires counsel FTO before any release claim — downgrade
   "patent-free by construction" to "FTO pending."
4. **Kill condition:** if Gate 0 (empirical vs analytic agreement on ≥5 varied low-res spectra)
   fails within a bounded prototype, KILL RS³ and fall back to the existing additive calibration
   features (`TailorScore`, `ChanceMatchSurprise`, `RawScoreCal`).

### Go / no-go (authoritative)
| Item | Verdict | First falsifiable test |
|---|---|---|
| **C0** remove dead `IsolationWindowEfficiency` (alone) | **GO** | PIN byte-diff: column gone, all else identical; Percolator counts unchanged |
| **RS³ reformulated** (empirical decoy null of real score) | **GATE** | Gate 0 agreement ≥5 varied low-res spectra; fail → KILL |
| RS³ original (independent-bin CGF on `g(m)`) | **CUT** | superseded by reformulation |
| **A0** RS³ additive PIN features (`Rs3NegLog10P`/`Rs3StdScore`/`Rs3Delta`) | **GATE** (after Gate 0) | Astral keep-iff-≥flat → UPS1 must gain @1% entrapment-FDP |
| **A1** complementary-ion node bonus | **CUT** | redundant w/ `LongestComplementaryLadder`/`ComplementaryIonBalance` (pin.rs:181-182) + score mutation |
| **A2** predicted-intensity node weighting | **CUT** | collinear w/ `IntensitySignal`/`FragPred*` (pin.rs:187-188) + score mutation |
| **A3** additive local S/N (noise-floor) feature | **GATE** | one PIN column; UPS1/LFQ low-res A/B; flat = cut |
| **A4** delta_to_second margin | **CUT (exists)** | already emitted as `DeltaRankScore` (pin.rs:177) |
| **A5** own-trained spectral-angle | **GATE** | one PIN column; ablation vs frag battery on Astral; flat = cut |
| **A6** de novo log-prob (Casanovo, own-trained) | **GATE** (opt-in, last) | Astral opt-in PIN feature; speed budget holds; flat = cut |
| **B2** fragment-ion index (Comet-FI/Sage designs) | **GO** (substrate) | identical PSM output + measured speedup |
| **B3** chimeric depth sweep | **GO** (cheapest expansion; the separation substrate) | `ANDES_CHIMERIC_OVERLAP` sweep, Astral PSM gain at flat true-FDP |
| **B1** semi-tryptic (external-sort + serve scaling) | **GATE** (after B2) | full-DB throughput + PSM-level entrapment-FDP holds |
| **B4** mass-offset / restricted-open | **GATE** (after B2) | per-offset PSM-level entrapment-FDP |
| **group-FDR** | **CUT until built** | not implemented (Percolator `--only-psms`); build as own PR or validate B1/B4 at PSM-FDP |

### Corrected sequencing
1. **C0** (alone). 2. **RS³ Gate-0 prototype** (reformulated; research spike, decides go/no-go).
3. **B3 chimeric expansion** (cheap) — provides the expanded-space substrate so separation features
   are tested where they can actually pay (resolves the "flat on closed search" sequencing flaw).
4. **RS³ additive PIN features**, one at a time, on the expanded path. 5. **A3 / A5** behind
   ablation gates. 6. **B2 → B1 → B4** coverage scaling. Every step one variable, entrapment-FDP,
   Astral-gate-then-UPS1.

### Cut from the campaign
A1, A2, A4 (redundant/exists/mutation), the `--score rs3` ranking hook (defer indefinitely — ranking
on an under-validated score is the highest-risk lever), and any group-FDR dependency until built.

---

## Open questions for review/judge (now answered by the ruling above)

1. Is RS³ + complementary-ion (A0+A1) enough to move low-res UPS1 *before* B1 expansion, or is the closed-search "discriminator flat" result going to bite (i.e., must B1 come first)?
2. A2/A5 collinearity with the existing `strong` battery — real risk of flat-at-1%. Which single predicted-intensity feature pays, and does it duplicate `frag_llr_battery`?
3. Is B1's group-FDR a thin Percolator post-process (per the FDR-boundary rule) or does it need more?
4. Sequencing: should B2 (fragment index) precede B1, since B1 at full-DB scale may still be slow without it?
