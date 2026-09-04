# RS³ — Renewal-Saddlepoint Spectral Significance

**Date:** 2026-06-29
**Status:** ⚠️ REVISED BY CAMPAIGN RULING (2026-06-29). The original construction below (analytic
saddlepoint over an *independent-bin* CGF, calibrating the float surrogate `g(m)`) was found
to (a) drop site-visit covariance → mis-calibrated null, and (b) calibrate a surrogate, not
andes's real emitted integer/loss/cleavage-credited score. **Superseded by the decoy-calibrated
empirical-null reformulation** — see "Campaign ruling" in
[2026-06-29-unique-scoring-campaign-plan.md](2026-06-29-unique-scoring-campaign-plan.md).
The §2 math is retained as background; do NOT implement it as written. Gate 0 is now a hard go/no-go.
**Scope:** a novel, own, patent-free **per-spectrum significance calibration** of andes's
RawScore, emitted as a Percolator feature and (per user decision, "full RS³") allowed to
influence ranking if the math warrants it.

> **Prior-rule override (explicit).** [docs/plans/2026-06-21-scoring-research-loop.md](../../plans/2026-06-21-scoring-research-loop.md)
> recorded a deliberate rule: *"NOT a per-spectrum calibration — that was the MS-GF+ reach
> we dropped."* The user has explicitly overridden that rule for RS³ on the grounds that RS³
> is (a) patent-free (it is **not** the MS-GF+ generating function), (b) own-data/own-model,
> and (c) computed by a mathematically distinct route. This document supersedes that rule for
> the spectral-significance direction. The rule's underlying *evidence* (per-spectrum
> calibration was deprioritized because the GBDT discriminator was flat on closed search) is
> respected by the validation plan: RS³ must prove a net PSM gain at honest entrapment-FDP or
> it does not ship.

---

## 1. Problem & objective

andes leads the field on high-res (Astral) and TMT, but **trails MS-GF+ (Java) ~5% PSMs on
low-res UPS1** (post-soft-matching, v0.3.0). The investigated root cause: andes's summed
RawScore is **not calibrated per spectrum**. At 0.5 Da fragment tolerance the random-match
rate is high and *heterogeneous across spectra* (peak density, precursor mass, charge,
peptide length all shift the null), so a single raw-score scale is poorly comparable — the
exact failure mode MS-GF+'s spectral probability fixes.

The two known ways to get per-spectrum calibration are both unavailable or weak here:
- **MS-GF+ generating function** — exact, but **patented** (US 8,639,447 B2, UC, active to
  2030-07-25). Off-limits for an Apache project. This is why andes stripped its `lnSpecEValue`
  features.
- **X!Tandem / MSFragger / Andromeda approximations** (exponential tail-fit, uniform
  binomial) — patent-clear, but andes already ships the structural pieces they bolt on
  (`longest_b/y`, `ChanceMatchSurprise`, `TailorScore`, `RawScoreCal`), so they add little.

**Objective:** a per-spectrum significance for andes's *own* RawScore, computed from andes's
*own* trained tables, by a route that is neither the patented generating function nor a crude
uniform binomial — and prove it converts the low-res calibration gap into PSMs at 1% true
entrapment-FDP without regressing Astral/TMT.

---

## 2. The mathematics (RS³)

### 2.1 The score as a renewal functional

andes's RawScore is a sum over cleavage sites
([psm_score.rs:234-368](../../../crates/scoring/src/scoring/psm_score.rs)):

```
Score = Σ_s g(m_s),     g(m) ≡ prefix_node(m) + suffix_node(M − m)
```

where `m_s` are the candidate's prefix masses, `M` is the precursor neutral mass, and `g(m)`
is the **spectrum-determined score landscape** — exactly what the existing node kernel
([scored_spectrum.rs:909-987](../../../crates/scoring/src/scoring/scored_spectrum.rs))
returns if a fragment lands at prefix mass `m` (soft-matched rank LLR, both ion sides).

Under a **random peptide of mass M**, the cleavage sites are the partial sums of an i.i.d.
amino-acid-mass random walk conditioned to hit `M`. So the *null* score is an **additive
functional of a renewal process**:

```
Score = Σ_m V(m) · g(m),     V(m) ∈ {0,1} = "a random peptide of mass M cleaves at mass m"
```

The MS-GF+ generating function enumerates *peptides by score* via DP (the patented object).
RS³ instead characterizes the *site-visit process* `V` and convolves *score contributions by
mass position* — a different object with a classical, tail-accurate solution.

### 2.2 Renewal visit density (spectrum-independent, precomputed)

By renewal theory, the probability a random peptide of mass `M` cleaves at mass `m` is the
**renewal bridge** density

```
ρ(m) = u(m) · u(M − m) / u(M),     u(m) = Σ_a π(a) · u(m − mass(a)),  u(0)=1
```

where `π(a)` is the null amino-acid distribution (uniform, or natural abundance) and `u` is
the classical renewal function (probability a renewal occurs at `m`). In the bulk
`u(m) → 1/μ`, `μ ≈ 111.1 Da`. **`u` and `ρ` depend only on the AA mass alphabet — not on the
spectrum and not on any score** → precompute once at startup. This is *not* the patented
score-generating function (it never counts reconstructions per score).

### 2.3 Exact null via cumulant generating function + saddlepoint

Treating site visits as independent (controlled approximation; renewal correlations are
short-range, ~one AA mass, and the bridge term absorbs the first-order structure):

```
K(θ) = Σ_m log( 1 − ρ(m) + ρ(m) · e^{θ g(m)} )       (cumulant generating function)
```

For an observed top score `T`, the spectral p-value `P(Score ≥ T)` is the **Lugannani–Rice
saddlepoint tail**:

```
solve  K'(θ̂) = T
w = sgn(θ̂) · sqrt( 2(θ̂ T − K(θ̂)) )
ν = θ̂ · sqrt( K''(θ̂) )
p ≈ 1 − Φ(w) + φ(w) · (1/ν − 1/w)
```

`K, K', K''` are sums over mass bins; `θ̂` is a 1-D Newton solve (a handful of `K` evals).
Saddlepoint is exact to high order in the tail — the regime p-values live in.

### 2.4 What RS³ emits

- `Rs3NegLog10P` = `−log10(p̂)` — the primary calibrated significance.
- `Rs3StdScore` = `w` (signed, standardized) — a calibrated z-like score, robust when
  `p̂` underflows.

Both are computed for every scored candidate (or at least rank-1 + the delta to rank-2).

### 2.5 Why this is novel and patent-distinct

RS³ **subsumes the field as special cases**: Andromeda's binomial is `ρ(m)≡p`,
`g(m)∈{c,0}`; X!Tandem/MSFragger's tail-fit is an *empirical estimate* of this same tail.
RS³ is the principled, trained, analytic version — position-resolved `ρ(m)`, trained `g(m)`,
exact-in-the-tail. It never enumerates peptide reconstructions by score (the US 8,639,447
claim); it computes a renewal functional's tail by saddlepoint, a textbook statistical
technique on a different mathematical object. No engine is known to do this.

---

## 3. Architecture & integration

RS³ is a self-contained scoring module that consumes existing per-spectrum data and emits new
features. It does **not** modify the existing RawScore/strong-score math (the safe,
proven-additive pattern), except for the optional ranking hook in §3.4.

### 3.1 New module: `crates/scoring/src/scoring/rs3.rs`

Responsibilities, each independently testable:

| Unit | Input | Output | Depends on |
|---|---|---|---|
| `RenewalTable::build(π, max_mass, bin)` | AA mass alphabet, null dist | `u(m)`, cached | nothing (startup-once, `OnceLock`) |
| `score_landscape(scored_spec, M, bin) -> Vec<f32>` | a prepared spectrum + precursor mass | `g(m)` over bins | existing node kernel (`scored_spectrum.rs`) |
| `visit_density(&RenewalTable, M) -> Vec<f32>` | renewal table, M | `ρ(m)` over bins | RenewalTable |
| `cgf(ρ, g) -> Cgf` then `saddlepoint_tail(&Cgf, T)` | ρ, g, observed score | `(neg_log10_p, w)` | nothing |

`bin` (mass-axis resolution) defaults to the model tolerance (`mme`) so RS³ self-scales to
instrument resolution — parameter-free, consistent with the soft-matching σ design.

### 3.2 Wiring into the search/scoring path

- The score landscape `g(m)` is evaluated **once per spectrum** (not per candidate) — a
  single pass over the mass axis reusing the node kernel. Per-candidate cost is then one
  table lookup of `T` plus the (shared) saddlepoint on the spectrum's `(ρ, g)` — `θ̂` depends
  on `T` so the Newton solve is per-candidate but cheap.
- Hook point: where per-PSM features are assembled in
  [match_engine.rs:1449-1948](../../../crates/search/src/match_engine.rs) (the
  `PsmFeatures` population, [psm.rs:12-234](../../../crates/search/src/psm.rs)).

### 3.3 PIN features (additive — the safe lever)

Add to the `PsmFeatures` struct and the PIN writer
([pin.rs](../../../crates/output/src/pin.rs)):

- `Rs3NegLog10P` (rank-1 and any retained candidate)
- `Rs3StdScore` (`w`)
- `Rs3DeltaNegLog10P` = rank-1 `−log10 p` minus rank-2 `−log10 p` (a calibrated deltaCn
  analog; rank-1 row only, else 0 — mirrors the existing `DeltaRankScore` convention)

These are additive columns; no existing feature changes. **Also: remove the dead
`IsolationWindowEfficiency` column** (always 0.0 — confirmed in the PIN inventory) while
touching `pin.rs`/`psm.rs`, since it is pure noise to Percolator.

### 3.4 Optional ranking hook ("full RS³")

Per the user's "full RS³" decision, expose RS³ as a selectable ranking score in addition to
the PIN feature: extend the existing `--score {rank,strong,auto}` selector (in
[andes.rs](../../../crates/andes/src/bin/andes.rs)) with `rs3`, which ranks candidates by
`w` (or `−log10 p`) instead of the raw summed LLR. **Default stays `auto`** — `rs3` ranking
ships only if the §5 A/B shows it beats `auto` at honest FDP. The additive PIN features ship
independently of the ranking hook (they are the low-risk path).

---

## 4. Validation gate 0 — numerical prototype (before any engine wiring)

A throwaway prototype (Rust test or a small script) that, on a handful of real low-res
spectra:

1. Builds `ρ(m)` and `g(m)` for a known top PSM.
2. Computes RS³ `p̂` via saddlepoint.
3. Computes a **brute-force empirical null**: score N=10⁵–10⁶ random mass-feasible decoy
   peptides against the same spectrum, take the tail fraction ≥ `T`.
4. **Gate:** RS³ `p̂` agrees with the brute-force tail within ~0.2 log10 units across several
   spectra of varying peak density and charge. If it does not, diagnose (independence
   approximation? bridge density? binning?) and add the §6 Edgeworth correction before
   proceeding. **No engine wiring until gate 0 passes.**

This gate is the cheapest possible test that the novel math is correct, and it isolates the
one real risk (the independence assumption) before any expensive benchmark.

---

## 5. Validation gate 1 — benchmark A/B (experiment-hygiene)

Per [feedback_experiment_hygiene] and [feedback_astral-gates-milestone-commits]: one variable
at a time, provenance-stamped (binary commit + model SHA + data SHA), at **1% true
entrapment-FDP** (not reported FDR), Percolator only.

1. **Astral gate first.** Add the RS³ PIN features (ranking unchanged), retrain nothing —
   the features are model-derived but require no retrain. A/B Astral with vs without the RS³
   columns. **Keep only if ≥ flat** (a calibration should be ~neutral where the null is
   already homogeneous). A regression here kills it.
2. **UPS1 low-res** (the target). Same binary, same features. **This is where the gain must
   appear** — expect the per-spectrum calibration to recover PSMs lost to null heterogeneity.
3. **TMT** (low-res CID, regression canary). Confirm no regression.
4. **Optional ranking hook** (`--score rs3`): only if 1–3 are net-positive, A/B `rs3` ranking
   vs `auto` on all three. Ship the ranking default flip only on a clean win.

Bank/drop rule: net PSM gain on UPS1 **and** no Astral/TMT regression at honest FDP → ship as
additive features (milestone commit on a feature branch). Ranking flip is a separate gate.

---

## 6. Risks & mitigations

| Risk | Mitigation |
|---|---|
| **Site-visit independence is an approximation** (adjacent sites are correlated) | The renewal-bridge `ρ(m)=u(m)u(M−m)/u(M)` absorbs first-order structure; residual is short-range. Gate 0 measures the error directly. If too large: add a second-order **Edgeworth correction** using the renewal pair-covariance, or a small empirical recalibration of `w`. |
| **Collinearity with existing calibration features** (`TailorScore`, `RawScoreCal`, `ChanceMatchSurprise`) | Percolator gains from collinear features are flat ([parity-tuning-lessons]). Add RS³ features one at a time; if flat, that is an honest negative result — RS³ does not ship. |
| **Saddlepoint instability on sparse spectra** (few candidates / near-degenerate `K''`) | Guard `θ̂` range and `K'' > 0`; fall back to `w` (or a normal approx) and flag low-confidence fits, degrading gracefully. |
| **Null AA distribution choice** (`π` uniform vs natural-abundance) | Treat as an internal constant, not a knob; test both in gate 0, pick the one matching the brute-force null; document. |
| **Reopens a deprioritized direction** | Honored by the hard FDP gate (§5) — RS³ ships only on a measured win, exactly addressing why per-spectrum calibration was dropped before. |

---

## 7. Independence / patent posture

- RS³ uses **andes's own trained tables** (`rank_dist`, soft-matched node scores) for `g(m)`
  and the **AA mass alphabet** for `ρ(m)` — no MS-GF+ values, no MS-GF+ geometry.
- It is **not** the generating function: it never determines "the number of peptide
  reconstructions at each score" (the US 8,639,447 independent-claim language). It computes a
  renewal functional's tail by saddlepoint. This is a *different mathematical object* and a
  classical statistical technique.
- Recommend a brief counsel FTO confirmation before any release claim, consistent with the
  existing [independence-license-status] track — but the design deliberately avoids the one
  known blocking patent.

---

## 8. Out of scope

- New training data / unsupported-model gaps (Astral+TMT, timsTOF+HLA) — separate workstream
  (WS-2), tracked from the model-coverage gap analysis.
- Partition-geometry independence — shipped (own-geometry bundle).
- Any change to the FDR boundary — Percolator remains the only production FDR tool
  ([feedback_andes_fdr_boundary]); RS³ is an andes *scoring/feature* change only.
- Comet-XCorr / Andromeda-binomial reimplementations — explicitly *not* pursued; RS³ is the
  own, novel alternative.

---

## 9. Definition of done

1. Gate 0: RS³ saddlepoint `p̂` matches the brute-force decoy null within tolerance on
   varied spectra (unit-tested).
2. RS³ module + additive PIN features (`Rs3NegLog10P`, `Rs3StdScore`, `Rs3DeltaNegLog10P`),
   dead `IsolationWindowEfficiency` removed; parity tests stay green (additive only).
3. Gate 1: net PSM gain on UPS1 at 1% true entrapment-FDP, no Astral/TMT regression.
4. (Optional) `--score rs3` ranking ships only on a clean A/B win vs `auto`.
5. One milestone commit on a feature branch; single closing PR ([iteration-shipping-model]).
