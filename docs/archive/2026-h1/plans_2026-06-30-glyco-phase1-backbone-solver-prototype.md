# Glyco Phase 1 — Oxonium Gate + Backbone-Mass Solver Prototype (Implementation Plan)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a measurement-only prototype that, on real intact-N-glycopeptide stepped-HCD data, detects glyco spectra (oxonium gate) and solves the **peptide backbone mass** de novo from the Y-ion ladder — then measures whether that solved backbone falls inside andes's existing candidate-mass window. This single number (the "searchable-backbone rate") is the **go/no-go gate** for building the full `--glyco` search mode.

**Architecture:** A new, self-contained `andes-glyco` crate with three pure modules (glycan masses, oxonium gate, backbone solver) + a probe binary. It **reuses andes's spectrum reader and candidate-mass-window function read-only**, touches **no scoring code**, emits **no PSMs/PIN**, and changes **nothing** in the shipped engine. Methodologically it follows the reconciled v2 design (Y-complementary voting + core-Y quorum + bounded top-K) and is informed by — but copies nothing from — the published StrucGP method (Shen et al., Nat Commun 2021).

**Tech Stack:** Rust (workspace `msgf-rust/`), `cargo test`, existing crates `input` (mzML reader) and `search`/`model` (candidate-mass bounds). No new external dependencies.

## Global Constraints

- **Phase 1 is MEASUREMENT-ONLY.** No change to scoring, ranking, PIN/TSV output, FDR, or any shipped behavior. No PSMs claimed. The deliverable is a number + a report.
- **Clean-room.** All glyco math (monosaccharide masses, oxonium m/z, trimannosyl-core Y pattern) is implemented from **public/published values** (Unimod, the StrucGP/pGlyco3 papers' *described methods*) — never transcribed from the StrucGP binary at `/Users/yperez/Downloads/main.exe_extracted`. Do not copy its constant arrays, `structureCoding` alphabet, branch tables, or scoring exponents.
- **Determinism.** No `rand`, no system time. Any sampling uses a seeded LCG (the codebase forbids nondeterminism).
- **Target regime (MVP):** N-linked glycopeptides, sequon **N-X-S/T (X≠P)**, **stepped-HCD / HCD** (collisional). No O-glyco, no de-novo structure in Phase 1.
- **Dataset:** **PXD025455** (Lubman, NASH-HCC serum; intact N-glycopeptides; stepped-HCD; Q Exactive HF; Thermo .raw; Byonic-validated with `.pepXML` ground truth). FTP: `https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/05/PXD025455/`.
- **The go/no-go gate (v2 §6.3):** on ground-truth glyco spectra, the solved backbone mass must fall inside `candidate_nominal_bounds`' opened window for **≥ 0.70** of spectra, **AND** the rate must hold (not collapse) in the **sparse stratum** (spectra with ≤1 core-Y ion). Below that → the de-novo-backbone architecture is refuted for this engine; stop before building the mode.
- **Repo:** all code in `msgf-rust/` (git). Branch before committing; never commit to the default branch directly.

**Source design:** [2026-06-13-andes-glyco-mode-design-v2-improvements.md](../internal-docs/docs/specs/2026-06-13-andes-glyco-mode-design-v2-improvements.md) (authoritative), [2026-06-13-andes-glycopeptide-search-research-plan.md](../internal-docs/docs/specs/2026-06-13-andes-glycopeptide-search-research-plan.md), [2026-06-20-glyco-neutral-loss-and-maxsbm.md](../2026-06-20-glyco-neutral-loss-and-maxsbm.md).

---

## Roadmap (context — only Phase 1 is in this plan)

- **Phase 0 (shipped):** additive neutral-loss/glyco scoring primitive for a *known* glycoform (inert until a glyco model is trained). Not in scope here.
- **Phase 1 (THIS PLAN):** oxonium gate + backbone-mass solver, measurement-only → the searchable-backbone gate.
- **Phase 2:** `--glyco` MVP — composition DB + sequon filter + Route-i backbone search reusing `candidate_gen`/`score_psm`.
- **Phase 3:** glyco scoring — retained-core ladder via the primitive + Y-core/oxonium matcher → additive PIN features.
- **Phase 4:** 2D glycan FDR — decoy (FragPipe fragment-scramble **or** StrucGP-style random precursor-mass-shift) + sequential TDC `max(q_pep,q_gly)`; the 1%×1% true-entrapment-FDP gate vs MSFragger-Glyco + pGlyco3 lives here.

---

## File Structure

| File | Responsibility | Status |
|---|---|---|
| `crates/andes-glyco/Cargo.toml` | New crate manifest; deps: `input`, `model` (read-only) | create |
| `crates/andes-glyco/src/lib.rs` | crate root, module exports | create |
| `crates/andes-glyco/src/glycan_mass.rs` | monosaccharide + oxonium + core-Y mass constants (clean-room) | create |
| `crates/andes-glyco/src/oxonium.rs` | the two-part oxonium gate (pure fn on a spectrum) | create |
| `crates/andes-glyco/src/backbone.rs` | Y-complementary mass-difference-ladder backbone solver | create |
| `crates/andes-glyco/src/bin/glyco_probe.rs` | the validation harness binary (real-data measurement) | create |
| `Cargo.toml` (workspace) | add `crates/andes-glyco` to members | modify |
| `crates/andes-glyco/tests/data/` | tiny synthetic spectra fixtures for unit tests | create |

**Recon note (do before Task 5, no commit):** read `crates/input/src/mzml.rs` for the exact spectrum/peak API (how to get an MS2 scan's `(m/z, intensity)` peaks + precursor m/z + charge), and find andes's precursor candidate-mass window function (grep `candidate_nominal_bounds` / `mass_window` in `crates/search/src/`). Task 5's interface must match the real signatures; if they differ from the shapes below, adapt the harness, not the math.

---

### Task 1: Glyco mass constants (clean-room, from public values)

Monosaccharide residue masses, core oxonium m/z, and the trimannosyl-core Y-ion offsets — all public/published numbers. These drive the gate and the solver.

**Files:**
- Create: `crates/andes-glyco/Cargo.toml`, `crates/andes-glyco/src/lib.rs`, `crates/andes-glyco/src/glycan_mass.rs`
- Modify: `Cargo.toml` (workspace members)

**Interfaces:**
- Produces:
  - `pub const HEXNAC: f64 = 203.07937;` `pub const HEX: f64 = 162.05282;` `pub const FUC: f64 = 146.05791;` `pub const NEUAC: f64 = 291.09542;` `pub const NEUGC: f64 = 307.09033;` (monoisotopic residue masses, water already subtracted)
  - `pub const PROTON: f64 = 1.0072765;`
  - `pub const CORE_OXONIUM_MZ: [f64; 5] = [138.05496, 168.06552, 186.07608, 204.08665, 366.13947];` (HexNAc-derived + HexHexNAc; singly-charged m/z)
  - `pub const CORE_Y_STEPS: [f64; 5] = [203.07937, 406.15874, 568.21156, 730.26438, 892.31720];` (Y1..Y5 offsets above the peptide backbone: +HexNAc, +2HexNAc, +2HexNAc+Hex, +2HexNAc+2Hex, +2HexNAc+3Hex — the trimannosyl core ladder)
  - `pub const MONO_STEPS: [f64; 5] = [HEXNAC, HEX, FUC, 365.13219, 324.10565];` (single-monosaccharide + common combo steps the Y-ladder walks)

- [ ] **Step 1: Write the failing test (constants are self-consistent)**

```rust
// crates/andes-glyco/src/glycan_mass.rs  (#[cfg(test)] mod tests)
#[test]
fn core_y_steps_are_cumulative_core() {
    // Y2 = Y1 + HexNAc; Y3 = Y2 + Hex; Y4 = Y3 + Hex; Y5 = Y4 + Hex
    assert!((CORE_Y_STEPS[1] - (CORE_Y_STEPS[0] + HEXNAC)).abs() < 1e-4);
    assert!((CORE_Y_STEPS[2] - (CORE_Y_STEPS[1] + HEX)).abs() < 1e-4);
    assert!((CORE_Y_STEPS[3] - (CORE_Y_STEPS[2] + HEX)).abs() < 1e-4);
    assert!((CORE_Y_STEPS[4] - (CORE_Y_STEPS[3] + HEX)).abs() < 1e-4);
}
```

- [ ] **Step 2: Create the crate + run the test to see it fail**

Create `crates/andes-glyco/Cargo.toml`:
```toml
[package]
name = "andes-glyco"
version = "0.3.0"
edition = "2021"

[dependencies]
input = { path = "../input" }
model = { path = "../model" }
```
Add `"crates/andes-glyco"` to the `members` list in the workspace `Cargo.toml`. Create `src/lib.rs` with `pub mod glycan_mass;`.
Run: `cargo test -p andes-glyco core_y_steps_are_cumulative_core`
Expected: FAIL (constants not yet defined).

- [ ] **Step 3: Add the constants (verbatim values above) to `glycan_mass.rs`**

- [ ] **Step 4: Run the test to confirm it passes**

Run: `cargo test -p andes-glyco core_y_steps_are_cumulative_core`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git checkout -b glyco-phase1
git add crates/andes-glyco/ Cargo.toml
git commit -m "feat(andes-glyco): new crate + clean-room glyco mass constants"
```

---

### Task 2: Two-part oxonium gate (TDD)

Decide whether an MS2 scan is glyco. v2 §1: fire iff (A) summed core-oxonium intensity ≥ `min_frac` × base peak **AND** (B) ≥2 distinct core oxonium ions present above an absolute floor (≥1% base peak).

**Files:**
- Create: `crates/andes-glyco/src/oxonium.rs` (+ `pub mod oxonium;` in lib.rs)

**Interfaces:**
- Consumes: `glycan_mass::{CORE_OXONIUM_MZ}`.
- Produces:
  - `pub struct OxoniumEvidence { pub fired: bool, pub summed_frac: f32, pub n_core_ions: u8 }`
  - `pub fn oxonium_gate(peaks: &[(f64, f32)], min_frac: f32, tol_ppm: f64) -> OxoniumEvidence` — `peaks` = (m/z, intensity), sorted or not; `tol_ppm` match window (default 20.0) with a 0.01 Th floor.

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn oxonium_gate_fires_on_glyco_spectrum() {
    // base peak intensity 100; two core oxonium ions (204.087, 138.055) at 15 each = 30% summed
    let peaks = vec![(500.0, 100.0), (204.0867, 15.0), (138.055, 15.0), (700.0, 5.0)];
    let e = oxonium_gate(&peaks, 0.10, 20.0);
    assert!(e.fired); assert_eq!(e.n_core_ions, 2); assert!(e.summed_frac >= 0.29);
}
#[test]
fn oxonium_gate_silent_on_nonglyco() {
    let peaks = vec![(500.0, 100.0), (700.0, 5.0), (204.5, 30.0)]; // 204.5 not within tol of 204.0867
    assert!(!oxonium_gate(&peaks, 0.10, 20.0).fired);
}
```

- [ ] **Step 2: Run to confirm both fail**

Run: `cargo test -p andes-glyco oxonium_gate`
Expected: FAIL (fn undefined).

- [ ] **Step 3: Implement**

```rust
use crate::glycan_mass::CORE_OXONIUM_MZ;
pub struct OxoniumEvidence { pub fired: bool, pub summed_frac: f32, pub n_core_ions: u8 }

pub fn oxonium_gate(peaks: &[(f64, f32)], min_frac: f32, tol_ppm: f64) -> OxoniumEvidence {
    let base = peaks.iter().map(|p| p.1).fold(0.0f32, f32::max).max(1e-9);
    let floor = 0.01 * base;
    let mut summed = 0.0f32; let mut n = 0u8;
    for &mz in CORE_OXONIUM_MZ.iter() {
        let tol = (mz * tol_ppm / 1e6).max(0.01);
        // best matching peak for this oxonium m/z
        let mut best = 0.0f32;
        for &(pmz, pi) in peaks { if (pmz - mz).abs() <= tol && pi > best { best = pi; } }
        if best >= floor { summed += best; n += 1; }
    }
    let frac = summed / base;
    OxoniumEvidence { fired: frac >= min_frac && n >= 2, summed_frac: frac, n_core_ions: n }
}
```

- [ ] **Step 4: Run to confirm both pass**

Run: `cargo test -p andes-glyco oxonium_gate`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add crates/andes-glyco/src/
git commit -m "feat(andes-glyco): two-part oxonium gate (v2 design)"
```

---

### Task 3: Y-complementary backbone-mass solver (TDD)

The make-or-break math. From the (low-energy) MS2 peaks + precursor neutral mass, recover the **peptide backbone mass** by finding a chain of peaks spaced by monosaccharide masses that matches the trimannosyl-core Y-ladder, then subtracting the core. v2: complementary-mass voting + core-Y quorum + top-K. This is a clean-room implementation of the published Y-ion-fingerprint idea (StrucGP/pGlyco3), not a copy.

**Files:**
- Modify: `crates/andes-glyco/src/backbone.rs` (+ `pub mod backbone;`)

**Interfaces:**
- Consumes: `glycan_mass::{CORE_Y_STEPS, PROTON}`.
- Produces:
  - `pub struct BackboneCandidate { pub backbone_mass: f64, pub core_y_hits: u8, pub votes: u32 }`
  - `pub fn solve_backbone(peaks: &[(f64, f32)], precursor_neutral: f64, precursor_z: u8, tol_ppm: f64, top_k: usize) -> Vec<BackboneCandidate>` — returns up to `top_k` candidate backbone masses, sorted by (core_y_hits desc, votes desc). Empty if no core-Y quorum (≥2 hits) is reached.

- [ ] **Step 1: Write the failing test (synthetic Y-ladder)**

```rust
// peptide backbone neutral mass = 1500.0; build singly-charged Y ions = backbone + core step + PROTON
#[test]
fn solve_backbone_recovers_known_mass_from_core_ladder() {
    let bb = 1500.0;
    let mut peaks: Vec<(f64,f32)> = crate::glycan_mass::CORE_Y_STEPS.iter()
        .map(|&s| (bb + s + crate::glycan_mass::PROTON, 50.0)).collect();
    peaks.push((bb + crate::glycan_mass::PROTON, 40.0)); // Y0 = bare backbone+H
    peaks.push((999.9, 5.0)); // noise
    let precursor = bb + 1444.53; // backbone + a HexNAc2Hex5 glycan (~1444.53)
    let out = solve_backbone(&peaks, precursor, 2, 20.0, 5);
    assert!(!out.is_empty());
    assert!((out[0].backbone_mass - bb).abs() < 0.02, "got {}", out[0].backbone_mass);
    assert!(out[0].core_y_hits >= 2);
}
```

- [ ] **Step 2: Run to confirm it fails**

Run: `cargo test -p andes-glyco solve_backbone_recovers_known_mass_from_core_ladder`
Expected: FAIL (fn undefined).

- [ ] **Step 3: Implement the voting solver**

```rust
use crate::glycan_mass::{CORE_Y_STEPS, PROTON};
pub struct BackboneCandidate { pub backbone_mass: f64, pub core_y_hits: u8, pub votes: u32 }

pub fn solve_backbone(peaks: &[(f64, f32)], precursor_neutral: f64, precursor_z: u8,
                      tol_ppm: f64, top_k: usize) -> Vec<BackboneCandidate> {
    use std::collections::HashMap;
    // bucket votes for a backbone mass at 0.01-Da resolution
    let mut votes: HashMap<i64, (u32, [bool; 6])> = HashMap::new();
    let neutral = |mz: f64, z: f64| (mz - PROTON) * z; // peak neutral mass at charge z
    for &(pmz, _pi) in peaks {
        for z in 1..=precursor_z.max(1) {
            let pn = neutral(pmz, z as f64);
            if pn <= 0.0 || pn > precursor_neutral { continue; }
            // this peak could be Y0 (bare backbone) or Y_r (backbone + core step r)
            // candidate backbone = pn - core_step (Y0: step 0)
            for (ri, step) in std::iter::once(0.0).chain(CORE_Y_STEPS.iter().copied()).enumerate() {
                let bb = pn - step;
                if bb <= 0.0 { continue; }
                let key = (bb * 100.0).round() as i64;
                let tol_key = ((bb * tol_ppm / 1e6).max(0.01) * 100.0).round() as i64;
                // accumulate into nearby keys within tolerance
                for k in (key - tol_key)..=(key + tol_key) {
                    let e = votes.entry(k).or_insert((0, [false; 6]));
                    e.0 += 1;
                    if ri < 6 { e.1[ri] = true; }
                }
            }
        }
    }
    let mut cands: Vec<BackboneCandidate> = votes.into_iter().map(|(k, (v, hits))| {
        BackboneCandidate { backbone_mass: k as f64 / 100.0,
            core_y_hits: hits.iter().filter(|&&h| h).count() as u8, votes: v }
    }).filter(|c| c.core_y_hits >= 2).collect();   // core-Y quorum
    cands.sort_by(|a, b| b.core_y_hits.cmp(&a.core_y_hits).then(b.votes.cmp(&a.votes)));
    cands.dedup_by(|a, b| (a.backbone_mass - b.backbone_mass).abs() < 0.05);
    cands.truncate(top_k);
    cands
}
```

- [ ] **Step 4: Run to confirm it passes**

Run: `cargo test -p andes-glyco solve_backbone_recovers_known_mass_from_core_ladder`
Expected: PASS.

- [ ] **Step 5: Add a negative test (no quorum → empty), then commit**

```rust
#[test]
fn solve_backbone_empty_without_core_quorum() {
    let peaks = vec![(700.0, 50.0), (1234.5, 50.0)]; // no core-Y ladder
    assert!(solve_backbone(&peaks, 2500.0, 2, 20.0, 5).is_empty());
}
```
```bash
git add crates/andes-glyco/src/
git commit -m "feat(andes-glyco): Y-complementary backbone-mass voting solver (v2 design, clean-room)"
```

---

### Task 4: Stage PXD025455 data + ground truth

Get the real spectra + per-scan truth needed for the validation. Measurement task; deliverable is staged files + a manifest, not code.

**Files:**
- Create: `crates/andes-glyco/tests/data/PXD025455.manifest.md` (records the staged paths + how truth was obtained)

**Interfaces:**
- Produces: on disk — one converted `.mzML` (a single PXD025455 fraction) + a per-scan ground-truth table `truth.tsv` with columns `scan, backbone_peptide_mass, precursor_mz, precursor_z, glycan_composition`.

- [ ] **Step 1: Download one fraction + the Byonic truth**

From `https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/05/PXD025455/`: fetch one `.raw` (e.g. an HCC serum fraction) and its Byonic `.pepXML` if present. Convert `.raw → .mzML` with `ThermoRawFileParser` (the project uses it; available via the `thermorawfileparser.sif` container or a local install). Record exact filenames in the manifest.

- [ ] **Step 2: Build the per-scan truth table**

If a Byonic `.pepXML` exists: parse it to `truth.tsv` (scan → identified peptide → compute its monoisotopic backbone mass; glycan composition string; precursor m/z + z from the mzML). If no Byonic file is downloadable, generate truth by running **pGlyco3 or MSFragger-Glyco** on the fraction (per the v2 §6.1 fallback) and parse its glycoPSM output. Document which path was used.

- [ ] **Step 3: Sanity-check + commit the manifest (not the data)**

Confirm `truth.tsv` has ≥500 confident glyco scans and the mzML loads. Commit only the manifest (data files are large/external — `.gitignore` the data dir).
```bash
git add crates/andes-glyco/tests/data/PXD025455.manifest.md crates/andes-glyco/tests/data/.gitignore
git commit -m "docs(andes-glyco): PXD025455 staging manifest for the backbone-solver gate"
```

---

### Task 5: Validation harness — the searchable-backbone measurement

The probe binary: for each ground-truth glyco scan, oxonium-gate it, solve the backbone, and check whether any solved candidate's mass falls inside andes's precursor candidate window for the *true* backbone. Report the rate, **stratified by core-Y richness**.

**Files:**
- Create: `crates/andes-glyco/src/bin/glyco_probe.rs`

**Interfaces:**
- Consumes: `oxonium::oxonium_gate`, `backbone::solve_backbone`, the `input` crate mzML reader, and andes's candidate-mass-window predicate (recon: `candidate_nominal_bounds`/`mass_window` in `crates/search/src/` — wrap it as `fn in_candidate_window(solved: f64, truth: f64, tol_ppm: f64) -> bool`, the simplest faithful form being `|solved-truth| <= max(truth*tol_ppm/1e6, 0.01)` if the engine's window is a symmetric precursor tolerance; confirm against the real fn in recon and use it).
- Produces: stdout report + `crates/andes-glyco/tests/data/PHASE1_RESULT.md`.

- [ ] **Step 1: Recon the spectrum + candidate-window APIs (no commit)**

Read `crates/input/src/mzml.rs` (MS2 peaks + precursor m/z/z accessors) and grep `crates/search/src` for the candidate-mass window. Write the exact calls into the harness below; if the window is not a plain tolerance, replicate its bound check faithfully.

- [ ] **Step 2: Write the harness**

```rust
// crates/andes-glyco/src/bin/glyco_probe.rs
// usage: glyco_probe <fraction.mzML> <truth.tsv>
// For each truth glyco scan: oxonium_gate -> solve_backbone -> is any candidate within
// the candidate window of the TRUE backbone mass? Stratify by core_y_hits (<=1 vs >=2).
fn main() {
    let args: Vec<String> = std::env::args().collect();
    let (mzml, truth) = (&args[1], &args[2]);
    // load truth: scan -> (backbone_mass, prec_mz, prec_z)
    // load mzML MS2 scans (input crate)
    // counters: total, fired, searchable; and the same within the sparse (<=1 core-Y) stratum
    // for each truth scan:
    //   peaks = scan.ms2_peaks(); ev = oxonium_gate(&peaks, 0.10, 20.0);
    //   if !ev.fired { continue; }                      // soft route: only count gated scans
    //   prec_neutral = (prec_mz - PROTON)*prec_z;
    //   cands = solve_backbone(&peaks, prec_neutral, prec_z, 20.0, 5);
    //   searchable = cands.iter().any(|c| in_candidate_window(c.backbone_mass, truth_bb, 20.0));
    //   tally overall + by stratum (use cands.first().core_y_hits or the scan's core-Y count)
    // print: oxonium-fire rate, overall searchable rate, sparse-stratum searchable rate
}
```
(Fill the body using the exact `input`/`search` calls from recon; keep `in_candidate_window` faithful to the engine's precursor window.)

- [ ] **Step 3: Build + run on the staged fraction**

Run: `cargo run -p andes-glyco --bin glyco_probe -- crates/andes-glyco/tests/data/<fraction>.mzML crates/andes-glyco/tests/data/truth.tsv`
Expected: prints `oxonium-fire`, `searchable-overall`, `searchable-sparse` rates over ≥500 scans, no panic.

- [ ] **Step 4: Record the result**

Write `PHASE1_RESULT.md` with the three rates, the dataset/fraction, the ground-truth source (Byonic vs pGlyco3/MSFragger), and the binary commit. Then evaluate the gate:
- **PASS** (searchable-overall ≥ 0.70 AND sparse-stratum does not collapse, e.g. ≥ ~0.55): the de-novo backbone architecture works → proceed to a Phase-2 plan (`--glyco` MVP).
- **FAIL**: record the failure mode (oxonium gate too strict/loose? voting picks the wrong mass? sparse stratum collapses?). Try the documented knobs (oxonium `min_frac`, solver `top_k`, charge range) one at a time; if still < 0.70 within a bounded effort → **STOP**: the backbone-first architecture is refuted for this engine, and the full glyco mode should not be built on it.

- [ ] **Step 5: Commit**

```bash
git add crates/andes-glyco/src/bin/glyco_probe.rs crates/andes-glyco/tests/data/PHASE1_RESULT.md
git commit -m "feat(andes-glyco): searchable-backbone validation harness + Phase-1 gate result"
```

---

## What this plan deliberately excludes (future, gated on Phase-1 PASS)

- The `--glyco` search mode, glycan composition DB, sequon candidate filter, glyco scoring features, glycan decoy + 2D FDR, glyco model training, and any PIN/TSV change. All deferred to Phase 2–4 plans, each its own writing-plans cycle.
- Any change to shipped scoring/FDR/output. Phase 1 is inert by construction.

## Self-Review

- **Spec coverage:** oxonium gate (Task 2 ← v2 §1), backbone solver (Task 3 ← v2 §2 / StrucGP Y-ladder), the searchable-backbone gate ≥0.70 stratified (Task 5 ← v2 §6.3), measurement-only/no-scoring-change (Global Constraints), clean-room (Global Constraints + Task 1/3 notes). ✔
- **Placeholder scan:** the one genuine unknown — the exact `input`/`search` API signatures — is handled by an explicit recon step + a faithful `in_candidate_window` wrapper, not a placeholder. The harness body references real calls to be filled from recon (a measurement binary, not a unit under TDD). ✔
- **Type consistency:** `OxoniumEvidence`, `oxonium_gate`, `BackboneCandidate`, `solve_backbone`, the mass constants — defined once, reused with the same signatures across Tasks 2/3/5. ✔
- **Constraint compliance:** no scoring/PIN/FDR touched; clean-room constants; deterministic; PXD025455; the gate is the sole deliverable. ✔
