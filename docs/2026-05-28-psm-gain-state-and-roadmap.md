# msgf-rust — PSM-gain: state of play & roadmap

**Date:** 2026-05-28
**Objective function:** maximize **PSMs @1% FDR (Percolator)**, hold **speed ≈ current** (no wall regression beyond ~3%).
**Authoritative gate:** merge nothing to `dev` until Rust beats Java on **both PSMs and speed on all 3 datasets**.

This document consolidates every investigation note, spec, plan, and review under `docs/`
into one place: where we are, what we've proven, the ranked actions worth doing next,
and the realistic outcome. It supersedes the scattered notes as the entry point;
the source notes are cited inline and remain for detail.

---

## 1. Where we are

Speed is **won everywhere** — Rust is at or faster than Java on all three datasets
(V1 speed milestone, PR #36; iter31/32 perf cluster). The gate is blocked purely on
**PSM counts on PXD001819 and TMT**.

| Dataset | Rust PSMs @1% | Java PSMs @1% | PSM Δ | Speed | Gate |
|---|---:|---:|---:|:--:|:--:|
| Astral DDA | ~36,715 | ~36,271 (≤33,425 narrow-mod cfg) | **Rust ahead** | Rust ≥ Java | ✅ |
| PXD001819 (UPS1 yeast tryptic) | 14,808 | 14,974 | **−1.1%** | Rust ≥ Java | ❌ PSMs |
| TMT (PXD007683, a05058) | 9,605 | ~10,115–10,212 | **−5.0%** | Rust ≥ Java | ❌ PSMs |

**Blocker = PXD −1.1% and TMT −5.0% on PSMs.** TMT is the larger gap.

**Shipped:** V1 speed (PR #36) · README 3-dataset bench table (PR #39) · precursor-cal
opt-in `--precursor-cal {auto,on,off}` + Phase-1b calibration nudge (PR #40, +53 PXD).

**Parked branches (NOT merged):**
- `feat/chimeric-dda-plus` — shelved, see §4 dead-ends.
- `feat/id-rate-pxd001819-tmt` — scoring-path lever exhausted, see §2/§3.

---

## 2. What we've learned (the rules — do not relitigate)

Twelve-plus bench-gated iterations produced a hard-won empirical pattern. These are
the rules that constrain every future change:

**Rule 1 — The scoring path is at/near Java parity.** There is no remaining "scoring
bug" big enough to fix for a large gain:
- Peak-rank assignment (H2) is **bit-identical** to Java — debunked the I5 doc's
  central hypothesis (`parity-analysis/notes/2026-05-28-phase2-peak-rank-parity.md`).
- BSA fixture: 217/217 top-1 parity.
- Instrument / fragment-tolerance / calibration resolution already correct.

**Rule 2 — The Percolator audit dichotomy (n≈9+):**
- **Additive PIN columns** are *safe* but usually *flat* — Percolator already extracts
  the signal via correlated columns (EdgeScore iter19 flat; isotope-KL Phase 2 flat).
  **Exception:** a column adding a genuinely *orthogonal* dimension can win big
  (C-4 enzN/enzC/enzInt: **+1,718** Astral).
- **Modifying an existing column's distribution** *regresses* Percolator @1% FDR even
  when individually "more correct" — **UNLESS it moves Rust's top-1 selection toward
  Java's.** Every big win was top-1-restoring: iter20 high-res tolerance **+4,650**,
  iter33 edge-in-ranking **+3,705**, iter29 main-ion fix **+379**.
- **Bit-exact per-feature parity ≠ FDR parity** (iter23: bit-exact features regressed
  **−1,404**). Percolator trains on the *shape* of Rust's distributions, not on
  Java-identical values.

**Rule 3 — +10% via parity/scoring tweaks is off the table.** The scoring path is
essentially correct; a *large* gain requires an **algorithmic** change, not a parity
fix (`2026-05-28-phase2-peak-rank-parity.md`).

**Rule 4 — Multi-PSM-per-scan emission inflates FDR.** PSM-level target-decoy does not
encode the "few real peptides per scan" constraint (chimeric Phase 1: Astral +94%,
decoy fraction → 1.2:1).

---

## 3. The remaining gap, decomposed

Where are PXD's −1.1% and TMT's −5.0% actually lost? The PIN-diff + parity harnesses
(`benchmark/parity/`) localize it precisely:

- **72.8% of PXD ranking flips are `spec_e_swap_only`** — RawScore agrees with Java,
  but Rust's **`lnSpecEValue`** ranks a different top-1
  (`2026-05-25-spece-tail-exploration.md`).
- **Rust's lnSpecE tail is ~2.4× heavier than Java's**: 11.5% of agreement-bucket PSMs
  have `lnSpecE > −10` (a bad tail) vs Java's 1.8%; mean Δ +0.87.
- **TMT loses on ranking, not on rows:** Rust emits **more** PIN rows (29,765 vs
  28,201) but Percolator passes **fewer** at 1% FDR (9,605 vs 10,212) — pure ranking
  loss via lnSpecE / ion-feature drift (`2026-05-25-precursor-cal-ship-gates.md`).
- **Root cause:** the **GF score-distribution *shape*** (deconvolution / ion_existence
  → node scores → merged `ScoreDist`) diverges. Boundary/threshold fixes
  (`java_match_score`, `score≥max→0`, sink-retry removal, int subtraction) did **not**
  close it — same class as the BSA charge-3 `gf_java_parity` deconvolution gap.

**Conclusion:** the gate-blocking gap is concentrated in the **GF SpecEValue
distribution on low-quality / tail spectra** — the one place Rust still measurably
diverges *and* it directly drives the FDR-determining feature.

**Critical corollary (selection vs. rescoring):** the 72.8% `spec_e_swap_only` flips
are lost at **within-scan top-1 selection time** — Rust's `TopNQueue` orders by
`spec_e_value` first (known-divergences §4), so a wrong SpecE picks the wrong peptide
*before* Percolator ever sees the scan. **A Percolator feature can only re-rank PSMs
that were emitted; it cannot recover a top-1 that SpecE never selected.** Therefore
new features (Lever 3) improve discrimination among emitted PSMs but **cannot, by
themselves, recover the flips** — those need either a SpecE-*shape* fix (Lever 2a) or a
change to *what we select top-1 by* (Lever 2b).

**Related known divergence — `lnEValue`:** Rust's `lnEValue = ln(spec_e_value ×
num_distinct_peptides_at_length)` is **structurally ~27× off** from Java (median ratio
0.0368), *not* explained by mod-aware counting (known-divergences §2; iter33 diff: the
single most-divergent emitted feature, mean Δ −2.97, 99.2% of PSMs differ). Since
`lnEValue ≈ lnSpecE + f(peplen)` and we already emit `lnSpecE` and `peplen` separately,
this miscalibrated column is largely redundant — a candidate to drop or fix (Lever 3b).

---

## 4. Action plan (ranked by expected value: PSM-yield ÷ risk ÷ cost)

### Lever 1 — Search-space / mod-set parity audit on PXD + TMT  *(do first: cheap, safe, historically high-yield)*

**Why:** iter24 proved a search-space completeness fix (a missing Acetyl-Prot-N-term
mod) was a **large, pure-sensitivity gain with zero Percolator-distribution risk**
(+384 Astral). We have *not* audited whether PXD and TMT search the same effective
mod set + `.param` + enzyme config as Java.

**Experiment (measurement, then enable):**
- Diff the effective fixed+variable mod set, `.param`, enzyme, missed-cleavages, and
  isotope-error range, Rust vs Java, for **PXD** (UPS1 yeast tryptic) and **TMT**.
- For TMT specifically: confirm it runs `CID_HighRes_Tryp.param` (a prior bench bug ran
  it as `CID_LowRes`, `2026-05-25-precursor-cal-ship-gates.md`) **and** the full TMT
  mod set — TMT fixed on K + peptide-N-term, M-ox, and check TMT-on-S/T/H, pyro-Glu,
  N-term acetyl, deamidation against what Java enumerates.
- Any mod/param Java searches and Rust doesn't → enable it.

**Gate:** pure sensitivity gain (no existing-distribution change → no Rule-2 risk).
**Speed:** more mods → more candidates → some wall; spend the existing speed surplus,
measure, keep within ~3%.
**Expected:** could be large (iter24 precedent); at worst zero. **Best risk-adjusted EV.**

### Lever 2 — Recover the SpecE-driven top-1 flips  *(highest upside for the gate)*

Two routes to the same prize (the 72.8% `spec_e_swap_only` flips), per §3's corollary.

#### Lever 2a — GF SpecE-distribution-shape parity  *(research + targeted fix)*

**Why:** §3 — the gap lives here. 72.8% spec_e swaps + TMT ranking loss + 2.4× tail.
This is the only remaining *measurable* scoring divergence and it drives the FDR
feature directly. A fix that makes Rust's SpecE distribution *match* Java's is
**top-1-restoring at the SpecE level** → the Rule-2 *safe* class of modify-distribution
change (the 72.8% spec_e swaps would converge toward Java's picks).

**Experiment (localize before fixing):**
- GF trace on the top-20 agreement-bucket scans with `|Δ lnSpecE| > 3` on PXD + TMT;
  diff Java vs Rust GF score distributions with
  `benchmark/parity/diff_gf_distribution.py`.
- Attribute the tail to one of: **(a)** node-score shape from deconvolution (BSA
  charge-3 class), **(b)** per-bin merge widening (`group.rs::add_prob_dist` summing
  isotope/tolerance bins), or **(c)** `min_score` / threshold interaction.

**Fix:** targeted to whichever of (a)/(b)/(c) dominates. Bench-gate per dataset; revert
in place if Percolator regresses (Rule 2 risk is real for non-top-1-restoring shape
changes).
**Leading indicator:** drive the lnSpecE bad-tail % (Rust 11.5% → toward Java 1.8%)
*before* trusting the Percolator @1% delta.
**Speed:** if the fix needs the unthresholded GF or finer binning, spend the speed
surplus; measure wall.
**Expected:** the principled route to closing TMT −5% and PXD −1.1%.

#### Lever 2b — Change the within-scan top-1 selection criterion  *(cheaper, but Rule-2/Rule-4 risk)*

**Why:** if Rust's SpecE is unreliable, stop letting it pick the top-1. Two variants:
- **Re-rank selection by RawScore (or a RawScore+delta composite).** RawScore is
  parity-grade (peak ranks bit-identical; agreement-bucket RawScore Δ ≈ 0.35). **Risk:**
  RawScore reintroduces the peptide-length bias SpecE was *designed* to remove, and
  diverges from Java's SpecE-based selection — the Rule-2 regression zone (moving *away*
  from Java's top-1). Bench-gate hard.
- **Emit top-2 distinct peptides per scan (narrow mass window, NOT the chimeric wide
  window)** so Percolator sees Java's peptide and decides via the full feature model
  (incl. the Lever-3a delta). **Risk:** this is multi-PSM-per-scan emission → Rule 4
  FDR inflation; far smaller than chimeric's wide-window top-10, but must be hard-gated
  on the target:decoy ratio and the Astral control.

**Decision:** 2a is the principled fix; 2b is the cheap probe. Try 2b's top-2-narrow
variant first as a *measurement* (does emitting Java's peptide + a good delta let
Percolator recover the flips without inflating decoys?), and fall back to 2a if it
inflates.

### Lever 3 — Orthogonal additive Percolator features  *(safety net: zero regression risk, modest yield)*

Cheap, safe, do alongside Lever 1. Add **new dimensions uncorrelated** with existing
columns (the C-4 win pattern, not the EdgeScore/isotope-KL flat pattern). Note these
sharpen discrimination among *emitted* PSMs but cannot recover the §3 selection flips.

- **Lever 3a — `DeltaRawScore` = RawScore(best) − RawScore(2nd-best distinct peptide).**
  ✅ **BENCHED — clean win, keep** (branch `feat/delta-raw-score`, commit `bea5d697`;
  bench `2026-05-28-delta-raw-score-bench.md`). Captured as a per-spectrum scalar during
  candidate scoring (so it's populated even at `top_n=1` without feeding the GF
  `min_score` → no existing column changes). **+129 PXD / +12 TMT / +104 Astral @1% FDR,
  zero wall cost, no regression**, decoy structure unchanged. Closes the PXD gap to Java
  from −1.1% → **−0.25%**. Not mergeable alone (gate needs PXD AND TMT). The C-4 pattern
  realized: a new orthogonal dimension that wins, unlike the flat EdgeScore/isotope-KL.
- **Lever 3b — drop `lnEValue`.** ❌ **BENCHED — noise, discarded.** −8 PXD / +9 TMT /
  +18 Astral (sub-0.2%, within Percolator noise) and it *costs* PXD. The fix variant is
  not clean either: Java's `getEValue = specProb × numPeptides` and Rust's
  formula+length-arg already match (HIGH-2, 2026-05-18); the residual ~27× divergence is
  a `num_distinct` count difference + the SpecE tail, and `lnEValue` is redundant with
  `lnSpecE + peplen`. **lnEValue question closed; keep the column (schema stays
  Java-faithful).**
- Further candidates if needed: explained-intensity-fraction-at-top-1; spectrum-quality
  score. Bench-gate; expected flat-to-small, never negative.

### Dead-ends — do NOT re-attempt without a genuinely new idea

| Path | Why dead | Evidence |
|---|---|---|
| Chimeric DDA+ (Phase 1+2) | Inflates FDR; back-end GF is single-precursor-centered; even a correct precursor-hypothesis + shared-fragment rebuild helps only wide-window data, is no-op on Astral, **net-negative on TMT** → never clears the gate | `2026-05-28-chimeric-phase1/phase2-bench.md` |
| Per-feature bit-exact Java parity | Regresses Percolator (−1,404) | iter23 |
| Modify-distribution scoring fixes that don't restore top-1 (R-3, C-5b, units, edge-in-score) | n=8 regressions | memory audit log |
| Peak-rank assignment fix (I5 H2) | Bug doesn't exist — already bit-identical | `2026-05-28-phase2-peak-rank-parity.md` |
| Extending MassCalibrator | Gaps are downstream of the calibrator | `2026-05-25-precursor-cal-ship-gates.md` |

---

## 5. Sequencing & realistic outcome

**Progress (2026-05-28):** Lever 3a benched and kept (+129/+12/+104, free); Lever 3b
benched and discarded (noise). Net: PXD now −0.25% to Java, Astral still ahead, **TMT
still −4.9% — the remaining gate blocker.** 3a/3b can't move TMT (confidence/redundant
features don't recover SpecE selection flips).

**Order (remaining):** the gate now hinges on **TMT**. Next: **Lever 1 TMT mod/param
config audit** (cheap first probe — does Rust search the same effective mod set as Java
on TMT? the iter24-class check), then **Lever 2** for the SpecE-driven flips — the 2b
top-2-narrow *probe* (cheap), falling back to the 2a GF-shape fix (~1–2 weeks,
research-grade) if 2b inflates decoys.

**Realistic future outcome:**
- **Lever 1** finding a mod/param gap (plausible — iter24 precedent) is a direct
  PXD+TMT sensitivity gain that could, on its own, flip the gate. Lowest-risk, do it first.
- **Lever 2** is the best-evidence route to closing the SpecE tail → plausibly closes
  TMT −5% and PXD −1.1% and beats Java by a few %, **clearing the gate**.
- **+10% (the original stretch) is not reachable by these** — it needs a new
  algorithmic capability (chimeric done right, or a new candidate-generation/scoring
  model), a separate research project that *also* would not help TMT.
- **Speed:** every lever is framed to spend the existing speed surplus; hold wall ≈
  current; revert any lever that breaks the ~3% wall budget.

**Definition of done (gate cleared):** PXD ≥ Java **and** TMT ≥ Java on PSMs, Astral
still ≥ Java, speed ≥ Java on all three.

---

## Source notes consolidated here

- `parity-analysis/notes/2026-05-28-phase2-peak-rank-parity.md` — scoring at parity; H2 null.
- `parity-analysis/notes/2026-05-26-score-psm-trace-findings.md` — I5 label-flip trace (H1/H2/H3).
- `parity-analysis/notes/2026-05-25-spece-tail-exploration.md` — lnSpecE tail + GF-shape root cause.
- `parity-analysis/notes/2026-05-25-precursor-cal-ship-gates.md` — G1 gates, TMT ranking loss.
- `parity-analysis/notes/2026-05-28-chimeric-phase1-bench.md` / `…-phase2-bench.md` — chimeric refutation.
- `parity-analysis/diff/iter33/report.md` — Astral PIN-diff agreement-bucket feature deltas.
- `superpowers/specs|plans/2026-05-28-id-rate-pxd001819-tmt-*.md` — the ID-rate phased plan (Phase 1b shipped; Phase 2 superseded by the H2-null finding).
- `superpowers/specs/2026-05-28-chimeric-dda-plus-integration-design.md` — chimeric design (shelved).
</content>
</invoke>
