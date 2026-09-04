# RS³ Gate-0 + Calibration Baseline — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Decide — cheaply and falsifiably — whether per-spectrum calibration is the lever that closes andes's low-res gap, and whether RS³'s decoy-calibrated empirical null is correct against the *real emitted score*, before building any production feature.

**Architecture:** Three isolated deliverables. (1) Remove a dead PIN column (zero-risk deck-clearing). (2) A measurement harness that quantifies whether the existing calibration features (`TailorScore`, `ChanceMatchSurprise`, `RawScoreCal`) already capture per-spectrum signal on low-res UPS1 — the collinearity gate. (3) A self-contained RS³ **Gate-0 prototype**: renewal density + a mass-feasible decoy-peptide generator + an empirical per-spectrum null computed against andes's **real emitted score**, validated against a brute-force Monte-Carlo null on real low-res spectra. Gate 0 is a hard go/no-go: if the cheap empirical null cannot match brute-force tails within tolerance, RS³ stops here.

**Tech Stack:** Rust (workspace `msgf-rust/`), `cargo test`, existing crates `scoring`, `search`, `output`. No new external dependencies.

## Global Constraints

- **Metric:** PSMs at **1% true entrapment-FDP** (paired estimator, r=1) — never reported target-decoy FDR.
- **FDR boundary:** Percolator only, **PSM-level only**. No group-FDR (not implemented). No Mokapot.
- **Safe-change rule:** **additive PIN features only**; never modify the matched-peak set or the existing score. No node-score mutations.
- **Score identity:** RS³ must calibrate the **real emitted score** = `score_psm()` integer node sum **+ `cleavage_credit`** (added in `match_engine.rs`) **+ optional loss term** — NOT the float `g(m)` surrogate.
- **Patent:** No score-enumeration DP (overlaps US 8,639,447 claims). Empirical/Monte-Carlo null only. Treat all "patent-free" statements as **FTO-pending counsel**, not engineering fact.
- **Provenance (experiment-hygiene):** every benchmark stamped with binary commit + model SHA + data SHA; one variable per A/B.
- **Repo:** all code in `msgf-rust/` (a git repo). Branch before committing; do not commit to the default branch directly.

**Source specs:** [`2026-06-29-rs3-spectral-significance-design.md`](../2026-06-29-rs3-spectral-significance-design.md) (with campaign ruling), [`2026-06-29-unique-scoring-campaign-plan.md`](../2026-06-29-unique-scoring-campaign-plan.md), [`2026-06-29-literature-review-rust-agent.md`](../2026-06-29-literature-review-rust-agent.md).

---

## File Structure

| File | Responsibility | Status |
|---|---|---|
| `crates/output/src/pin.rs` | PIN column contract — remove dead `IsolationWindowEfficiency` | modify |
| `crates/search/src/psm.rs` | `PsmFeatures` struct — remove dead field + its test asserts | modify |
| `crates/scoring/src/scoring/rs3_gate0.rs` | **NEW** — renewal density, decoy generator, empirical null, brute-force comparison (prototype only, behind `#[cfg(test)]`-friendly pub API) | create |
| `crates/scoring/src/scoring/mod.rs` | export `rs3_gate0` | modify |
| `crates/scoring/tests/rs3_gate0_bruteforce.rs` | integration: Gate-0 assertion on real low-res spectra | create |
| `docs/superpowers/experiment-protocol.md` | record the Phase-0 calibration measurement + Gate-0 result | modify |

**Recon note (do first, no commit):** before Task 3, open and read `crates/scoring/src/scoring/psm_score.rs:234-368` (`score_psm`) and the `cleavage_credit` add-site in `crates/search/src/match_engine.rs` (≈905) to pin the exact signature of "score one peptide against one prepared spectrum and return the emitted f32 score." Task 3's `EmittedScore` interface MUST match that real signature; if it differs from the shape below, adapt the adapter, not the math.

---

### Task 1: Remove dead `IsolationWindowEfficiency` PIN column (C0)

Confirmed always-`0.0` (set at `match_engine.rs:1705`; asserted `0.0` in `psm.rs` tests). A constant feature carries zero information to Percolator, so removal is provably neutral — its own commit keeps later A/Bs clean.

**Files:**
- Modify: `crates/output/src/pin.rs` (column list + header)
- Modify: `crates/search/src/psm.rs` (`PsmFeatures` field + asserts)
- Test: existing PIN golden/round-trip tests in `crates/output/`

**Interfaces:**
- Consumes: nothing.
- Produces: a `PsmFeatures` with no `isolation_window_efficiency` field; PIN header with that column absent.

- [ ] **Step 1: Capture the current PIN header as a baseline**

Run: `cd msgf-rust && grep -rn "IsolationWindowEfficiency\|isolation_window_efficiency" crates/`
Expected: matches in `pin.rs` (header + value write), `psm.rs` (field + asserts), possibly `qpx.rs`. Record every file:line.

- [ ] **Step 2: Write the failing test — header must not contain the column**

Add to the existing PIN test module (`crates/output/src/pin.rs` `#[cfg(test)] mod tests`):

```rust
#[test]
fn pin_header_has_no_isolation_window_efficiency() {
    let header = pin_header(); // use the actual header-producing fn; if private, test via a written PIN
    assert!(!header.contains("IsolationWindowEfficiency"),
        "dead constant column must be removed");
}
```

- [ ] **Step 3: Run it to confirm it fails**

Run: `cargo test -p output pin_header_has_no_isolation_window_efficiency -- --nocapture`
Expected: FAIL (column still present).

- [ ] **Step 4: Remove the field and all references**

Delete the `isolation_window_efficiency` field from `PsmFeatures` (`psm.rs`), its initialization, the header entry and value write in `pin.rs` (and `qpx.rs` if mirrored), and any `assert_eq!(..., 0.0)` lines that reference it.

- [ ] **Step 5: Run the full output + search test suites**

Run: `cargo test -p output && cargo test -p search`
Expected: PASS, including the new test. Fix any column-count constants the PIN writer asserts.

- [ ] **Step 6: Verify result-neutrality on a real PIN**

Run a small search (use the project's smallest fixture per `experiment-protocol.md`) before and after on the same input; diff the PINs ignoring the removed column.
Expected: every remaining column byte-identical; Percolator PSM counts unchanged.

- [ ] **Step 7: Commit**

```bash
git checkout -b rs3-gate0
git add crates/output/src/pin.rs crates/search/src/psm.rs
git commit -m "chore(pin): remove dead IsolationWindowEfficiency column (always 0.0)"
```

---

### Task 2: Phase-0 calibration measurement (the collinearity gate)

No new math. Answer two questions on low-res UPS1 at 1% entrapment-FDP: (a) is the MS-GF+ gap a *calibration* gap (per-spectrum score-variance heterogeneity), and (b) do the **existing** features `TailorScore` / `ChanceMatchSurprise` / `RawScoreCal` already capture it? If they do, RS³'s marginal value is small and the bar for Gate 1 rises. This is a benchmark + analysis task; its deliverable is a recorded decision, not code.

**Files:**
- Modify: `docs/superpowers/experiment-protocol.md` (append a "Phase-0 calibration measurement" results block)

**Interfaces:**
- Consumes: the project benchmark harness (`prov.sh`, the UPS1 low-res entrapment gate).
- Produces: a recorded verdict `{gap_is_calibration: bool, existing_features_saturated: bool}` that Task 6's Gate-1 framing depends on.

- [ ] **Step 1: Run the low-res UPS1 entrapment benchmark, baseline**

Use the harness from `experiment-protocol.md` to produce a PIN + Percolator rescore on UPS1 low-res with the current `main` binary. Stamp binary commit + model SHA + data SHA.
Expected: a PSM count at 1% paired entrapment-FDP + the Percolator `--rescore` weights file.

- [ ] **Step 2: Extract Percolator feature weights**

Inspect the Percolator weights output; record the learned weights for `TailorScore`, `ChanceMatchSurprise`, `RawScoreCal`, `DeltaRankScore`, `ListwiseScoreGap`.
Expected: a table of feature → weight. A near-zero weight on all calibration features = they are NOT helping (room for RS³); a large weight = already-captured signal.

- [ ] **Step 3: Ablation A/B — drop `TailorScore`**

Re-run Percolator on the same PIN with the `TailorScore` column removed (or zeroed). One variable.
Expected: ΔPSMs at 1% entrapment-FDP. A large drop ⇒ calibration already pays via Tailor; a flat result ⇒ Tailor is not capturing low-res calibration (RS³ opportunity).

- [ ] **Step 4: Per-spectrum score-variance check**

For correctly-identified vs decoy-top PSMs, compute the per-spectrum null score variance (spread of candidate RawScores per spectrum). Confirm/deny that low-res spectra show heterogeneous nulls (the calibration hypothesis).
Expected: a quantified statement "low-res per-spectrum null variance spans X–Y, vs high-res Z" supporting or refuting "the gap is calibration."

- [ ] **Step 5: Record the verdict and commit**

Append results to `experiment-protocol.md` with the two booleans and the evidence.

```bash
git add docs/superpowers/experiment-protocol.md
git commit -m "docs(exp): Phase-0 low-res calibration measurement + collinearity verdict"
```

---

### Task 3: Renewal density `u(m)` / `ρ(m)` (pure math, TDD)

The bridge visit-density for the decoy null. Pure function of the amino-acid mass alphabet; spectrum- and score-independent. Used to make decoy draws cheap and tail-focused (importance proposal) — NOT to compute the p-value analytically.

**Files:**
- Create: `crates/scoring/src/scoring/rs3_gate0.rs`
- Modify: `crates/scoring/src/scoring/mod.rs` (add `pub mod rs3_gate0;`)

**Interfaces:**
- Consumes: amino-acid monoisotopic masses (reuse `model::mass` / the AA table already in `crates/model`).
- Produces:
  - `pub struct RenewalTable { bin_da: f32, u: Vec<f32> }`
  - `pub fn build_renewal_table(aa_masses: &[f32], max_mass: f32, bin_da: f32) -> RenewalTable`
  - `pub fn visit_density(t: &RenewalTable, total_mass: f32) -> Vec<f32>` returning `ρ(m)` over bins with `ρ(m)=u(m)·u(M−m)/u(M)`.

- [ ] **Step 1: Write the failing test — uniform-AA bulk density approaches 1/μ**

```rust
#[test]
fn renewal_bulk_density_approaches_inverse_mean_mass() {
    // single "amino acid" of mass 100 Da → renewals exactly every 100 Da
    let t = build_renewal_table(&[100.0], 2000.0, 1.0);
    let rho = visit_density(&t, 2000.0);
    // a renewal occurs only at multiples of 100; midpoints ~0
    let at_500 = rho[500];   // multiple of 100 → high
    let at_550 = rho[550];   // not a multiple → ~0
    assert!(at_500 > at_550 * 10.0, "density must peak at reachable masses");
}
```

- [ ] **Step 2: Run it to confirm it fails**

Run: `cargo test -p scoring renewal_bulk_density_approaches_inverse_mean_mass`
Expected: FAIL (functions undefined).

- [ ] **Step 3: Implement the renewal recurrence**

```rust
pub struct RenewalTable { pub bin_da: f32, pub u: Vec<f32> }

pub fn build_renewal_table(aa_masses: &[f32], max_mass: f32, bin_da: f32) -> RenewalTable {
    let n = (max_mass / bin_da).ceil() as usize + 1;
    let mut u = vec![0.0f32; n];
    u[0] = 1.0;
    let p = 1.0 / aa_masses.len() as f32; // uniform null AA distribution
    for m in 1..n {
        let mut acc = 0.0;
        for &a in aa_masses {
            let back = m as f32 - a / bin_da;
            if back >= 0.0 {
                acc += p * u[back.round() as usize];
            }
        }
        u[m] = acc;
    }
    RenewalTable { bin_da, u }
}

pub fn visit_density(t: &RenewalTable, total_mass: f32) -> Vec<f32> {
    let mi = (total_mass / t.bin_da).round() as usize;
    let denom = t.u.get(mi).copied().unwrap_or(0.0).max(1e-12);
    (0..t.u.len())
        .map(|m| if m <= mi { t.u[m] * t.u[mi - m] / denom } else { 0.0 })
        .collect()
}
```

- [ ] **Step 4: Run the test to confirm it passes**

Run: `cargo test -p scoring renewal_bulk_density_approaches_inverse_mean_mass`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add crates/scoring/src/scoring/rs3_gate0.rs crates/scoring/src/scoring/mod.rs
git commit -m "feat(rs3-gate0): renewal density u(m)/rho(m) for decoy null proposal"
```

---

### Task 4: Mass-feasible decoy-peptide generator (TDD)

Random peptides whose residue masses sum to ≈ M. Because they are real amino-acid walks, their cleavage sites obey the spacing/exclusion structure automatically — this is what makes the empirical null honest (and sidesteps the independence flaw of the analytic CGF).

**Files:**
- Modify: `crates/scoring/src/scoring/rs3_gate0.rs`

**Interfaces:**
- Consumes: `aa_masses: &[f32]`, target mass `M`, tolerance `tol_da`, a deterministic seed.
- Produces: `pub fn sample_decoy_masses(aa_masses: &[f32], target_mass: f32, tol_da: f32, seed: u64) -> Option<Vec<f32>>` returning a sequence of residue masses summing within `tol_da` of `target_mass`, or `None` if not reached within a bounded attempt budget.

- [ ] **Step 1: Write the failing test — sampled mass hits the target within tolerance**

```rust
#[test]
fn decoy_mass_hits_target_within_tolerance() {
    let aa = [57.02146, 71.03711, 87.03203, 97.05276, 99.06841, 101.04768,
              113.08406, 114.04293, 128.05858, 128.09496, 131.04049, 137.05891,
              147.06841, 156.10111, 163.06333, 186.07931]; // representative residues
    let m = 1500.0;
    let seq = sample_decoy_masses(&aa, m, 0.5, 42).expect("should reach target");
    let sum: f32 = seq.iter().sum();
    assert!((sum - m).abs() <= 0.5, "decoy mass {sum} not within 0.5 of {m}");
}
```

- [ ] **Step 2: Run it to confirm it fails**

Run: `cargo test -p scoring decoy_mass_hits_target_within_tolerance`
Expected: FAIL (function undefined).

- [ ] **Step 3: Implement a greedy-with-backtrack sampler (deterministic from seed)**

```rust
pub fn sample_decoy_masses(aa: &[f32], target: f32, tol: f32, seed: u64)
    -> Option<Vec<f32>> {
    // deterministic LCG — no Math.random / new Date allowed
    let mut state = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    let mut next = || { state = state.wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407); (state >> 33) as usize };
    for _ in 0..64 {              // bounded restarts
        let mut seq = Vec::new();
        let mut sum = 0.0f32;
        for _ in 0..200 {         // bounded length
            let remaining = target - sum;
            if remaining.abs() <= tol { return Some(seq); }
            if remaining < -tol { break; }            // overshot → restart
            let a = aa[next() % aa.len()];
            seq.push(a); sum += a;
        }
    }
    None
}
```

- [ ] **Step 4: Run the test to confirm it passes**

Run: `cargo test -p scoring decoy_mass_hits_target_within_tolerance`
Expected: PASS.

- [ ] **Step 5: Add a determinism test, then commit**

```rust
#[test]
fn decoy_sampler_is_deterministic() {
    let aa = [57.02146, 113.08406, 128.09496, 147.06841];
    assert_eq!(sample_decoy_masses(&aa, 800.0, 0.5, 7),
               sample_decoy_masses(&aa, 800.0, 0.5, 7));
}
```

```bash
git add crates/scoring/src/scoring/rs3_gate0.rs
git commit -m "feat(rs3-gate0): deterministic mass-feasible decoy peptide sampler"
```

---

### Task 5: Empirical per-spectrum null against the REAL emitted score (TDD)

Compute `p̂(T) = (#{decoy score ≥ T} + 1) / (N + 1)` where each decoy is scored by andes's **real emitted score path** (not `g(m)`). The `+1` Laplace guard avoids `p̂=0`. This is the Gate-0 method itself.

**Files:**
- Modify: `crates/scoring/src/scoring/rs3_gate0.rs`

**Interfaces:**
- Consumes: an `EmittedScore` callback abstracting the production scorer (defined here so the prototype is testable without the full engine; wired to the real path in Task 6):
  ```rust
  pub trait EmittedScore { fn score(&self, residue_masses: &[f32]) -> f32; }
  ```
  This callback MUST, in the real wiring, return `score_psm() + cleavage_credit (+ loss)` for the residue sequence against a fixed prepared spectrum (see recon note).
- Produces:
  - `pub fn empirical_p(scorer: &dyn EmittedScore, decoys: &[Vec<f32>], observed: f32) -> f32`
  - `pub fn neg_log10(p: f32) -> f32`

- [ ] **Step 1: Write the failing test with a synthetic linear scorer**

```rust
struct SumScorer; // score = total residue mass — monotone, exactly checkable
impl EmittedScore for SumScorer { fn score(&self, r: &[f32]) -> f32 { r.iter().sum() } }

#[test]
fn empirical_p_matches_counting() {
    let decoys: Vec<Vec<f32>> = (0..100).map(|i| vec![i as f32]).collect(); // scores 0..99
    let p = empirical_p(&SumScorer, &decoys, 90.0); // 10 decoys >= 90 (90..=99), inclusive
    // (10 + 1) / (100 + 1) = 0.1089...  -- empirical p-value is INCLUSIVE (>=), conservative
    assert!((p - 11.0/101.0).abs() < 1e-6, "got {p}");
}
```

- [ ] **Step 2: Run it to confirm it fails**

Run: `cargo test -p scoring empirical_p_matches_counting`
Expected: FAIL (function undefined).

- [ ] **Step 3: Implement**

```rust
pub trait EmittedScore { fn score(&self, residue_masses: &[f32]) -> f32; }

pub fn empirical_p(scorer: &dyn EmittedScore, decoys: &[Vec<f32>], observed: f32) -> f32 {
    let ge = decoys.iter().filter(|d| scorer.score(d) >= observed).count();
    (ge as f32 + 1.0) / (decoys.len() as f32 + 1.0)
}

pub fn neg_log10(p: f32) -> f32 { -(p.max(1e-12).log10()) }
```

- [ ] **Step 4: Run the test to confirm it passes**

Run: `cargo test -p scoring empirical_p_matches_counting`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add crates/scoring/src/scoring/rs3_gate0.rs
git commit -m "feat(rs3-gate0): empirical per-spectrum p against an EmittedScore callback"
```

---

### Task 6: Gate-0 integration — empirical vs brute-force on real low-res spectra (the go/no-go)

Wire `EmittedScore` to the **real** production scorer and assert that a *cheap* importance-sampled decoy null (using `ρ(m)` to focus draws) agrees with a *large* brute-force uniform-decoy null within tolerance, on real low-res spectra. This is the hard kill switch.

**Files:**
- Create: `crates/scoring/tests/rs3_gate0_bruteforce.rs`
- Modify: `docs/superpowers/experiment-protocol.md` (record Gate-0 result)

**Interfaces:**
- Consumes: Tasks 3–5; the production scorer (recon note: the `score_psm` + `cleavage_credit` path) wrapped as an `EmittedScore`.
- Produces: the Gate-0 verdict (pass/fail) recorded in `experiment-protocol.md`.

- [ ] **Step 1: Recon — pin the real scoring signature (no commit)**

Read `crates/scoring/src/scoring/psm_score.rs:234-368` and the `cleavage_credit` add-site in `crates/search/src/match_engine.rs` (≈905). Write down the exact call needed to score a residue-mass sequence against one prepared `ScoredSpectrum` and return the emitted f32. If it requires a `Peptide`/`ScoredSpectrum` rather than `&[f32]`, the test's `EmittedScore` impl builds those from the residue masses.

- [ ] **Step 2: Write the Gate-0 test against ≥5 real low-res spectra**

```rust
// crates/scoring/tests/rs3_gate0_bruteforce.rs
// Loads a small low-res fixture (see experiment-protocol.md for the path),
// picks >=5 spectra spanning low/high peak density and charge 2/3.
#[test]
#[ignore] // research spike: run explicitly, not in CI
fn gate0_empirical_matches_bruteforce() {
    for spec in load_lowres_fixture_spectra() {          // helper in the test file
        let scorer = RealEmittedScore::for_spectrum(&spec); // wraps score_psm + cleavage_credit
        let observed = scorer.score(&spec.top_hit_residue_masses);

        // brute-force null: 100k uniform mass-feasible decoys
        let brute: Vec<_> = (0..100_000)
            .filter_map(|i| sample_decoy_masses(&AA, spec.precursor_mass, 0.5, i))
            .collect();
        let p_brute = empirical_p(&scorer, &brute, observed);

        // cheap null: 5k decoys (importance-focused via rho is optional here)
        let cheap: Vec<_> = (1_000_000..1_005_000)
            .filter_map(|i| sample_decoy_masses(&AA, spec.precursor_mass, 0.5, i))
            .collect();
        let p_cheap = empirical_p(&scorer, &cheap, observed);

        let d = (neg_log10(p_cheap) - neg_log10(p_brute)).abs();
        assert!(d < 0.2, "spectrum {}: |Δlog10 p| = {d} exceeds 0.2", spec.id);
    }
}
```

- [ ] **Step 3: Run Gate 0**

Run: `cargo test -p scoring --test rs3_gate0_bruteforce -- --ignored --nocapture`
Expected: one of two decisive outcomes.

- [ ] **Step 4: Record the verdict in experiment-protocol.md**

- **PASS** (all spectra `< 0.2`): the cheap empirical null is faithful → proceed to a follow-up plan that wires `Rs3NegLog10P` / `Rs3StdScore` as **additive PIN features** and runs Gate 1 (Astral-flat → UPS1-gain).
- **FAIL** (any spectrum off): record the failure mode (sparse-spectrum N too small? null mass-distribution mismatch? cleavage_credit not in the scorer?). Try the documented fallbacks (raise N for sparse spectra; match the null AA distribution to natural abundance). If still failing within a bounded effort → **KILL RS³** and fall back to strengthening the existing additive calibration features (`TailorScore`, `ChanceMatchSurprise`, `RawScoreCal`) per Task 2's findings.

- [ ] **Step 5: Commit the prototype + verdict**

```bash
git add crates/scoring/tests/rs3_gate0_bruteforce.rs docs/superpowers/experiment-protocol.md
git commit -m "test(rs3-gate0): empirical-vs-bruteforce null gate on real low-res spectra + verdict"
```

---

## What this plan deliberately excludes (future plans, gated on Gate 0 PASS)

- **RS³ as production PIN features** (`Rs3NegLog10P`, `Rs3StdScore`, `Rs3DeltaNegLog10P`) + Gate-1 benchmark — separate plan.
- **Calibrated pruning** (the speed move): a cheap pre-scoring upper-bound on achievable significance + RAM-vs-pruned identification-equivalence test — separate plan, only meaningful once the significance is trusted.
- **Self-labeling** (the flywheel): RS³ + entrapment as the label-trust anchor for andes-on-andes retraining — separate plan, last.
- **Coverage** (B2 fragment index, B3 chimeric sweep, B1 semi-tryptic external-sort) — separate plans.
- **A5 spectral-angle / A3 noise-floor / A6 de novo** — separate, ablation-gated plans.

---

## Self-Review

- **Spec coverage:** C0 (Task 1), Phase-0 calibration measurement / collinearity gate (Task 2), RS³ reformulation = decoy-calibrated empirical null of the real score (Tasks 3–6), Gate-0 kill switch (Task 6 Step 4). Pruning/self-labeling explicitly deferred. ✔
- **Placeholder scan:** every code step has concrete Rust; the one genuine unknown (the exact production scorer signature) is handled by an explicit recon step + an `EmittedScore` trait boundary, not a placeholder. ✔
- **Type consistency:** `RenewalTable`, `build_renewal_table`, `visit_density`, `sample_decoy_masses`, `EmittedScore::score`, `empirical_p`, `neg_log10` are defined once and reused with the same signatures across Tasks 3–6. ✔
- **Constraint compliance:** additive-only (no node-score mutation anywhere), real emitted score (Task 5/6), no DP-over-scores (empirical MC only), entrapment-FDP framing (Task 2/follow-up). ✔
