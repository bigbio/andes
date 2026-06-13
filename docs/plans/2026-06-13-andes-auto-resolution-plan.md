# Auto-detect instrument; remove `--instrument`; MGF-only frag/tol params — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove `--instrument` from the `andes` search CLI and make `--fragmentation` + a new `--fragment-tol-ppm/-da` hidden MGF-only extended params; metadata formats (mzML/.raw/.d) stay zero-config.

**Architecture:** The fragment-matching tolerance (today a hardcoded `20 ppm` high-res / `0.5 Da` low-res branch at 3 sites) is centralized into a single `RankScorer::feature_match_tolerance()` method with an optional override. The CLI sets the override and the resolution class only for metadata-less input (MGF); metadata-bearing input auto-detects as today. No data-driven inference pre-pass.

**Tech Stack:** Rust (cargo workspace), clap CLI, existing model store. Tests are `#[cfg(test)]` modules + `crates/*/tests/*.rs`.

**Spec:** `docs/specs/2026-06-13-andes-auto-resolution-design.md`

**Branch:** `feat/enzyme-support` (current). Build gate: `cargo build --workspace` + `cargo clippy --workspace -- -D warnings` stay clean. Pre-existing failing tests on this branch (NOT ours): `train_from_msnet::fragment_tolerance_override_changes_model` + 3 `match_engine_smoke` min_peaks tests — confirm any new failure is ours via `git stash`.

---

## Task 1: Centralize the feature-match tolerance on `RankScorer`

**Files:**
- Modify: `crates/scoring/src/scoring/rank_scorer.rs` (struct ~L28-50, `new()` ~L55, return literal)
- Test: `crates/scoring/src/scoring/rank_scorer.rs` (`#[cfg(test)]` module at bottom)

- [ ] **Step 1: Write the failing test**

Add to the `#[cfg(test)] mod tests` block in `rank_scorer.rs` (use the existing `testutil` helpers — `crate::testutil::high_res_param()` / `low_res_param()` exist with `mme: Ppm(20)` / `Da(0.5)`; if names differ, use the two builders at `crates/scoring/src/testutil.rs:46` and `:111`):

```rust
#[test]
fn feature_match_tolerance_defaults_by_resolution() {
    use model::tolerance::Tolerance;
    let hi = RankScorer::new(&crate::testutil::high_res_param());
    assert_eq!(hi.feature_match_tolerance(), Tolerance::Ppm(20.0));

    let lo = RankScorer::new(&crate::testutil::low_res_param());
    assert_eq!(lo.feature_match_tolerance(), Tolerance::Da(0.5));
}

#[test]
fn feature_match_tolerance_honors_override() {
    use model::tolerance::Tolerance;
    let mut s = RankScorer::new(&crate::testutil::low_res_param());
    s.set_fragment_tol_override(Some(Tolerance::Da(0.6)));
    assert_eq!(s.feature_match_tolerance(), Tolerance::Da(0.6));
    s.set_fragment_tol_override(None);
    assert_eq!(s.feature_match_tolerance(), Tolerance::Da(0.5));
}
```

If the testutil builder names differ, open `crates/scoring/src/testutil.rs` and use the actual `pub fn` names; both return a `Param` with the stated `mme`/instrument.

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p scoring feature_match_tolerance 2>&1 | tail -20`
Expected: FAIL — `no method named feature_match_tolerance` / `set_fragment_tol_override`.

- [ ] **Step 3: Implement the field + methods**

In the `RankScorer` struct (after the `max_rank: u32,` field), add:

```rust
    /// Optional CLI override for the fragment-matching tolerance (MGF input).
    /// `None` ⇒ derive from the model's instrument resolution class (the
    /// historical 20 ppm high-res / 0.5 Da low-res default).
    fragment_tol_override: Option<model::tolerance::Tolerance>,
```

In `new()`, at the final `Self { ... }` return literal, add the initializer:

```rust
            fragment_tol_override: None,
```

Add these methods inside `impl RankScorer` (near `param()` at ~L133):

```rust
    /// Set (or clear) the CLI fragment-tolerance override. Used by the search
    /// binary for metadata-less (MGF) input; metadata-bearing input leaves it
    /// `None` so behavior is byte-identical to auto-detection.
    pub fn set_fragment_tol_override(&mut self, tol: Option<model::tolerance::Tolerance>) {
        self.fragment_tol_override = tol;
    }

    /// Effective fragment-matching tolerance for Percolator-feature and
    /// peak-match counting. Returns the CLI override when set, else the
    /// instrument-derived default: 20 ppm for high-resolution analyzers
    /// (Kim et al., Nat Commun 5:5277, 2014), 0.5 Da for ion-trap low-res.
    pub fn feature_match_tolerance(&self) -> model::tolerance::Tolerance {
        use model::tolerance::Tolerance;
        self.fragment_tol_override.unwrap_or_else(|| {
            if self.param.data_type.instrument.is_high_resolution() {
                Tolerance::Ppm(20.0)
            } else {
                Tolerance::Da(0.5)
            }
        })
    }
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p scoring feature_match_tolerance 2>&1 | tail -20`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add crates/scoring/src/scoring/rank_scorer.rs
git commit -m "feat(scoring): RankScorer::feature_match_tolerance() + override (no behavior change)"
```

---

## Task 2: Route the 3 feature/matching sites through the method

**Files:**
- Modify: `crates/search/src/match_engine.rs` (matched_peak_keys ~L996-997; compute_psm_features ~L1069-1074)
- Modify: `crates/search/src/coisolation.rs` (~L123-124)

These changes are byte-identical when no override is set (the regression guard). `Tolerance::Ppm(20.0).as_da(mz)` == `mz * 20.0 / 1e6` and `Tolerance::Da(0.5).as_da(mz)` == `0.5`.

- [ ] **Step 1: Verify the regression suite passes BEFORE the change (baseline)**

Run: `cargo test -p search -p scoring -p output 2>&1 | tail -25`
Expected: the parity/golden tests pass (note any pre-existing failures listed in the header so you can distinguish them).

- [ ] **Step 2: Edit `match_engine.rs` matched_peak_keys (~L996-997)**

Replace:

```rust
    let tol_is_ppm = scorer.param().data_type.instrument.is_high_resolution();
    let tol = if tol_is_ppm { 20.0_f64 } else { 0.5_f64 };
    for p in &predicted {
        let tol_da = if tol_is_ppm { p.mz * tol / 1e6 } else { tol };
```

with:

```rust
    let feat_tol = scorer.feature_match_tolerance();
    for p in &predicted {
        let tol_da = feat_tol.as_da(p.mz);
```

(`use model::tolerance::Tolerance;` is not needed here — `as_da` is a method.)

- [ ] **Step 3: Edit `match_engine.rs` compute_psm_features (~L1069-1074)**

Replace:

```rust
    let feature_tol = if scorer.param().data_type.instrument.is_high_resolution() {
        20.0_f64 // ppm
    } else {
        0.5_f64 // Da
    };
    let feature_tol_is_ppm = scorer.param().data_type.instrument.is_high_resolution();
```

with:

```rust
    let feat_tol = scorer.feature_match_tolerance();
    let feature_tol_is_ppm = matches!(feat_tol, model::tolerance::Tolerance::Ppm(_));
    let feature_tol = feat_tol.raw_value(); // numeric value in the unit (ppm or Da)
```

Verify downstream uses of `feature_tol` in this function still compute the per-ion window the same way (search for `feature_tol` below L1074 — it is used as `theo_mz * feature_tol / 1e6` for ppm and as a constant for Da, matching `raw_value()`). If any site uses it differently, replace that arithmetic with `feat_tol.as_da(theo_mz)`.

- [ ] **Step 4: Edit `coisolation.rs` (~L123-124)**

Replace:

```rust
    let tol_is_ppm = scorer.param().data_type.instrument.is_high_resolution();
    let tol = if tol_is_ppm { 20.0_f64 } else { 0.5_f64 };
    for p in &predicted {
        let tol_da = if tol_is_ppm { p.mz * tol / 1e6 } else { tol };
```

with:

```rust
    let feat_tol = scorer.feature_match_tolerance();
    for p in &predicted {
        let tol_da = feat_tol.as_da(p.mz);
```

- [ ] **Step 5: Run the regression suite — must be unchanged**

Run: `cargo test -p search -p scoring -p output 2>&1 | tail -25`
Expected: same pass/fail set as Step 1 (no NEW failures). Confirm the PIN-schema and `score_psm_*_parity` golden tests still pass.

- [ ] **Step 6: Clippy + commit**

```bash
cargo clippy -p search -p scoring -- -D warnings 2>&1 | tail -5
git add crates/search/src/match_engine.rs crates/search/src/coisolation.rs
git commit -m "refactor(search): route feature-match tol through RankScorer::feature_match_tolerance() (byte-identical)"
```

---

## Task 3: Add `--fragment-tol-ppm/-da` to the search CLI and wire the override

**Files:**
- Modify: `crates/andes/src/bin/andes.rs` (search `Cli` struct ~L185-270; scorer construction ~L1129)
- Test: `crates/andes/tests/cli_smoke.rs`

- [ ] **Step 1: Add the CLI fields**

In the search `Cli` struct (next to `protocol`, ~L227), add:

```rust
    /// Fragment-matching tolerance in ppm for **MGF input only** (high-resolution
    /// MS/MS). Has no effect on mzML/.raw/.d (analyzer auto-detected). Mutually
    /// exclusive with `--fragment-tol-da`.
    #[arg(long = "fragment-tol-ppm", hide = true, conflicts_with = "fragment_tol_da")]
    fragment_tol_ppm: Option<f64>,

    /// Fragment-matching tolerance in Da for **MGF input only** (low-resolution
    /// ion-trap MS/MS). Has no effect on mzML/.raw/.d. Mutually exclusive with
    /// `--fragment-tol-ppm`.
    #[arg(long = "fragment-tol-da", hide = true)]
    fragment_tol_da: Option<f64>,
```

- [ ] **Step 2: Add a helper to resolve the override `Tolerance`**

Near the other free helpers (e.g. before `cli_flags_to_activation_instrument` ~L2953), add:

```rust
/// Resolve the CLI fragment-tolerance override (MGF only) into a `Tolerance`.
/// `--fragment-tol-ppm` ⇒ `Ppm`; `--fragment-tol-da` ⇒ `Da`; none ⇒ `None`.
fn cli_fragment_tol_override(
    fragment_tol_ppm: Option<f64>,
    fragment_tol_da: Option<f64>,
) -> Option<model::tolerance::Tolerance> {
    use model::tolerance::Tolerance;
    fragment_tol_ppm
        .map(Tolerance::Ppm)
        .or_else(|| fragment_tol_da.map(Tolerance::Da))
}
```

- [ ] **Step 3: Wire the override onto the scorer (metadata-less only)**

At the scorer construction (~L1129 `let scorer = RankScorer::new(&param);`), make it `mut` and set the override only when the input has no analyzer metadata (MGF, or detection returned no instrument). Use `is_mgf` (already computed ~L1038) and the detected-instrument `Option`:

```rust
    let mut scorer = RankScorer::new(&param);
    // Fragment-tol override applies to metadata-less input only. For
    // mzML/.raw/.d the analyzer is auto-detected, so the override is ignored
    // (warn once if the user passed it anyway).
    let frag_tol_override = cli_fragment_tol_override(cli.fragment_tol_ppm, cli.fragment_tol_da);
    let instrument_was_detected = detected_activation_instrument
        .map(|(_, inst)| inst.is_some())
        .unwrap_or(false);
    if frag_tol_override.is_some() {
        if instrument_was_detected {
            eprintln!(
                "WARN: --fragment-tol-* ignored — instrument auto-detected from metadata"
            );
        } else {
            scorer.set_fragment_tol_override(frag_tol_override);
        }
    }
```

(Placement note: this must come AFTER `detected_activation_instrument` is bound (~L1046) and after `param`/`scorer` are created. Keep `scorer` `mut`.)

- [ ] **Step 4: Add a CLI smoke test**

In `crates/andes/tests/cli_smoke.rs`, add a test that the flags parse and conflict (follow the file's existing `assert_cmd`/`Command` pattern; if it shells out to the built binary, use `--help` parse-only style already present). Example using the crate's existing harness:

```rust
#[test]
fn fragment_tol_flags_are_mutually_exclusive() {
    // Mirror the existing cli_smoke harness for invoking the binary.
    let out = run_andes(&["--spectrum", "x.mgf", "--database", "x.fasta",
                          "--output-pin", "x.pin",
                          "--fragment-tol-ppm", "20", "--fragment-tol-da", "0.5"]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("cannot be used with"));
}
```

If `cli_smoke.rs` has no `run_andes` helper, copy the invocation style already used by the nearest existing test in that file.

- [ ] **Step 5: Build + run**

Run: `cargo build -p andes 2>&1 | tail -5 && cargo test -p andes fragment_tol 2>&1 | tail -15`
Expected: build OK; the mutual-exclusion test PASSES.

- [ ] **Step 6: Commit**

```bash
git add crates/andes/src/bin/andes.rs crates/andes/tests/cli_smoke.rs
git commit -m "feat(cli): hidden MGF-only --fragment-tol-ppm/-da; override wired to scorer"
```

---

## Task 4: Remove `--instrument`, hide `--fragmentation`, rework model-selection routing

**Files:**
- Modify: `crates/andes/src/bin/andes.rs` (Cli `instrument` field ~L223-224; `--fragmentation` arg ~L222-223; routing ~L1079-1118; `cli_flags_to_activation_instrument` ~L2953; remove `cli_instrument_to_instrument_type` ~L2936, `parse_instrument` ~L3430, `Instrument` enum ~L66 if unused; `resolve_bundled_param` ~L2997 reference impl)
- Modify: tests that pass search `--instrument`: `crates/andes/tests/cli_smoke.rs`, `crates/andes/tests/store_selection_equivalence.rs`, `crates/search/tests/mass_calibrator_integration.rs`, `crates/search/tests/match_engine_smoke.rs`
- NOTE: the `train`/`train-from-msnet` subcommand `--instrument` (a `String` model-tag field, ~L431) is SEPARATE — do NOT remove it. `crates/model-train/tests/*` and `train_from_msnet.rs` reference the train flag and stay as-is.

- [ ] **Step 1: Write the failing behavior tests (MGF routing)**

Add to `crates/andes/tests/cli_smoke.rs` (or a new `mgf_routing.rs` integration test). These assert via stderr the selected model id (`eprintln!("Param model: {model_id} (from store)")` already logs it). Use a tiny MGF + FASTA fixture (reuse an existing one under `crates/input/tests/` / `test-fixtures/`; pick the smallest `.mgf`):

```rust
#[test]
fn mgf_no_flags_defaults_to_cid_lowres_with_warning() {
    let out = run_andes(&["--spectrum", FIXTURE_MGF, "--database", FIXTURE_FASTA,
                          "--output-pin", &tmp("a.pin")]);
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("cid_lowres"), "stderr: {err}");
    assert!(err.to_lowercase().contains("assuming"), "expected default warning: {err}");
}

#[test]
fn mgf_fragment_tol_ppm_selects_high_res() {
    let out = run_andes(&["--spectrum", FIXTURE_MGF, "--database", FIXTURE_FASTA,
                          "--output-pin", &tmp("b.pin"), "--fragment-tol-ppm", "20"]);
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("qexactive"), "stderr: {err}");
}
```

Define `FIXTURE_MGF`/`FIXTURE_FASTA` to existing fixtures used elsewhere in the andes tests (grep `cli_smoke.rs` for the paths it already uses).

- [ ] **Step 2: Run to verify failure**

Run: `cargo test -p andes mgf_ 2>&1 | tail -20`
Expected: FAIL (today MGF-no-flags resolves to `hcd_qexactive`, not `cid_lowres`).

- [ ] **Step 3: Remove the `instrument` field and hide `--fragmentation`**

Delete the search `Cli` `instrument: Instrument` field (~L223-224). On the `fragmentation` arg (~L222-223) add `hide = true` and update the doc comment to say it sets the activation method for MGF input only:

```rust
    /// Fragmentation/activation method for **MGF input only**. mzML/.raw/.d
    /// auto-detect this. Named: auto, CID, ETD, HCD, UVPD.
    #[arg(long, hide = true, default_value = "auto", value_parser = parse_fragmentation)]
    fragmentation: Fragmentation,
```

- [ ] **Step 4: Rework `cli_flags_to_activation_instrument` → metadata-less resolver**

Replace `cli_flags_to_activation_instrument` (~L2953) with a resolver that takes the detected activation (if any) and the CLI flags, and applies decision E. Drop the `instrument: Instrument` parameter entirely:

```rust
/// Resolve (activation, instrument) for model selection on metadata-less input
/// (MGF, or mzML/.raw with no analyzer metadata). Resolution class comes from
/// the `--fragment-tol-*` unit; activation from detected method, else
/// `--fragmentation`, else the class-consistent default. When nothing
/// disambiguates, decision E: CID / LowRes (→ `cid_lowres_tryp`) + a warning.
fn resolve_metadataless_selection(
    detected_activation: Option<ActivationMethod>,
    fragmentation: Fragmentation,
    fragment_tol_ppm: Option<f64>,
    fragment_tol_da: Option<f64>,
) -> (ActivationMethod, Option<InstrumentType>) {
    // Resolution class from the fragment-tol unit.
    let instrument: Option<InstrumentType> = if fragment_tol_ppm.is_some() {
        Some(InstrumentType::QExactive) // high-res family
    } else if fragment_tol_da.is_some() {
        Some(InstrumentType::LowRes)
    } else {
        None // unknown → build_selection_key defaults to LowRes
    };

    // Activation: detected > explicit --fragmentation > class-consistent default.
    let explicit = cli_fragmentation_to_activation_opt(fragmentation); // None when Auto
    let activation = detected_activation.or(explicit).unwrap_or_else(|| {
        match instrument {
            Some(InstrumentType::QExactive)
            | Some(InstrumentType::HighRes)
            | Some(InstrumentType::TOF) => ActivationMethod::HCD,
            _ => ActivationMethod::CID, // decision E default
        }
    });

    // Decision-E warning: nothing told us activation OR resolution.
    if detected_activation.is_none() && explicit.is_none() && instrument.is_none() {
        eprintln!(
            "WARN: MGF input with no --fragmentation/--fragment-tol; assuming \
             CID / low-res / 0.5 Da. Pass --fragmentation and --fragment-tol-ppm/-da \
             to override."
        );
    }
    (activation, instrument)
}
```

Add the `cli_fragmentation_to_activation_opt` helper (returns `None` for `Auto`, else the mapped method) next to the existing `cli_fragmentation_to_activation` (search for it; it maps Auto→CID — the new `_opt` version returns `None` for Auto so "auto" doesn't masquerade as an explicit choice):

```rust
fn cli_fragmentation_to_activation_opt(f: Fragmentation) -> Option<ActivationMethod> {
    match f {
        Fragmentation::Auto => None,
        Fragmentation::Cid  => Some(ActivationMethod::CID),
        Fragmentation::Etd  => Some(ActivationMethod::ETD),
        Fragmentation::Hcd  => Some(ActivationMethod::HCD),
        Fragmentation::Uvpd => Some(ActivationMethod::UVPD),
    }
}
```

- [ ] **Step 5: Rework the routing block in `main` (~L1079-1118)**

Replace the `(activation, instrument_opt)` resolution with metadata-first precedence:

```rust
        let (activation, instrument_opt): (ActivationMethod, Option<InstrumentType>) =
            match detected_activation_instrument {
                // Metadata gave us a concrete instrument → it wins (hot path).
                Some((method, Some(inst))) => {
                    eprintln!(
                        "Param resolver: auto-detected activation = {} (instrument = {}) from {}",
                        method.name(), inst.name(), spectrum_path.display()
                    );
                    (method, Some(inst))
                }
                // Detected activation but no instrument → metadata-less resolver,
                // seeding the detected activation.
                Some((method, None)) => resolve_metadataless_selection(
                    Some(method), cli.fragmentation, cli.fragment_tol_ppm, cli.fragment_tol_da,
                ),
                // No metadata at all (MGF, or detection failed).
                None => resolve_metadataless_selection(
                    None, cli.fragmentation, cli.fragment_tol_ppm, cli.fragment_tol_da,
                ),
            };
```

Remove the now-unused `auto_route_eligible` gating if it becomes dead (clippy will flag it). Keep `detected_activation_instrument` computation, but note it currently only runs when `cli.fragmentation == Auto` (~L1045) — change that gate to run detection regardless of `--fragmentation` (since `--fragmentation` is now MGF-only and must NOT disable detection):

```rust
    let auto_route_eligible = is_mzml || is_raw || is_d;
```

- [ ] **Step 6: Delete dead helpers**

Remove `cli_instrument_to_instrument_type` (~L2936) and `parse_instrument` (~L3430) and its tests (`parse_instrument_rejects_out_of_range_numeric` ~L3508). Remove the `Instrument` enum (~L66) IF `cargo build` reports it unused (the train subcommand uses a `String`, not this enum). Update the `--param-file` doc comment (~L193) to drop the `--instrument` reference. For `resolve_bundled_param` (~L2997, `#[allow(dead_code)]` reference impl) — it takes an `Instrument`; if `store_selection_equivalence` still uses it, leave it but feed it a fixed `Instrument::LowRes` from the test's own local enum (the test defines its own `enum Instrument` at L31, independent of the CLI one).

- [ ] **Step 7: Fix the broken tests**

`git stash`-free: just build and fix compile errors. Update:
- `crates/andes/tests/cli_smoke.rs`, `store_selection_equivalence.rs`, `crates/search/tests/mass_calibrator_integration.rs`, `match_engine_smoke.rs`: remove any `--instrument <x>` search args / `instrument:` CLI fields they set.
- `store_selection_equivalence.rs`: its reference ladder expects the OLD all-defaults default (`hcd_qexactive`). Update the no-flags expectation to `cid_lowres` (decision E), OR scope the equivalence test to explicit-flag cases only and add a separate assertion for the new default. Pick whichever keeps the test meaningful; document the change in the test comment.

Run: `cargo build --workspace 2>&1 | tail -20` — fix every error.

- [ ] **Step 8: Run the MGF tests + full suite**

Run: `cargo test -p andes mgf_ 2>&1 | tail -20` → PASS.
Run: `cargo test --workspace 2>&1 | tail -30` → only the pre-existing known failures remain (confirm via the header list; `git stash` + re-run if unsure a failure is pre-existing).

- [ ] **Step 9: Clippy + commit**

```bash
cargo clippy --workspace -- -D warnings 2>&1 | tail -10
git add -A
git commit -m "feat(cli): remove --instrument; --fragmentation MGF-only; metadata-first model selection"
```

---

## Task 5: Documentation — README separation + DOCS.md + workflow examples

**Files:**
- Modify: `README.md` (CLI summary table ~L139-158; Java-baseline note L33; TMT example L108-110; legacy-flags note L122)
- Modify: `DOCS.md` (§1 Scoring table L75-77; §4 Auto-detection L277-294; any `--instrument` mention)
- Modify: `benchmark/ci/run_bench_calauto_3ds.sh` IF it passes `--instrument` (grep first)

- [ ] **Step 1: README — remove `--instrument`, mark frag/tol MGF-only**

In the optional-args table (~L151-153): delete the `--instrument` row (L152). Change the `--fragmentation` row to note it is MGF-only. Then insert, immediately after the CLI summary tables, the MGF section from the spec §4:

```markdown
mzML, Thermo `.raw`, and Bruker `.d` are fully auto-detected — andes reads the
activation method and analyzer resolution from the file, so you pass no
fragmentation/instrument parameter for these formats.

### MGF input (extended parameters)

MGF files carry no activation or analyzer metadata, so you describe the
acquisition yourself:

| Parameter | When to pass | Example |
|---|---|---|
| `--fragmentation <CID\|ETD\|HCD\|UVPD>` | the activation method used | `--fragmentation HCD` |
| `--fragment-tol-ppm <X>` | high-resolution MS/MS (Orbitrap/TOF) | `--fragment-tol-ppm 20` |
| `--fragment-tol-da <X>`  | low-resolution MS/MS (ion trap)      | `--fragment-tol-da 0.5` |

If you pass none of these for an MGF file, andes assumes CID / low-res / 0.5 Da
and prints a warning. These parameters have no effect on mzML/`.raw`/`.d`.
```

- [ ] **Step 2: README — fix the example and notes**

- TMT example (~L108-110): remove `--instrument QExactive` (and the trailing `\`); the `.raw`/mzML input auto-detects. If the example uses MGF, replace `--instrument QExactive` with `--fragment-tol-ppm 20`.
- Java-baseline note (L33): drop `--instrument` from the listed matched args.
- Legacy-flags note (L122): remove `--instrument 3` from the `--fragmentation 3 --instrument 3 --protocol 4` example.

- [ ] **Step 3: DOCS.md updates**

- §1 table (L76): delete the `--instrument` row. Edit the `--fragmentation` row (L75) and add `--fragment-tol-ppm`/`--fragment-tol-da` rows, all marked **MGF-only**.
- §1 "Bundled default" note (L82): the all-defaults default is now described by format — metadata formats auto-detect; metadata-less defaults to `cid_lowres_tryp` with a warning.
- §4 (L277-294): state metadata formats are zero-config and the extended params are MGF-only; remove "`--instrument` value on the command line is ignored" language (the flag no longer exists); update the MGF paragraph (L284) to describe the extended-param contract + the CID/low-res default.
- Grep `DOCS.md` for any remaining `--instrument` and fix.

- [ ] **Step 4: CI benchmark script**

Run: `grep -n "instrument" benchmark/ci/run_bench_calauto_3ds.sh`
If it passes a search `--instrument`, remove it (auto-detected) or replace with `--fragment-tol-*` if the input is MGF. If no hit, no change.

- [ ] **Step 5: Verify no stale references remain**

Run: `grep -rn "\-\-instrument" README.md DOCS.md TRAIN.md benchmark/ crates/andes/src crates/andes/tests | grep -v "train\|model-train"`
Expected: no search-command `--instrument` references remain (train-subcommand mentions are fine).

- [ ] **Step 6: Commit**

```bash
git add README.md DOCS.md benchmark/ci/run_bench_calauto_3ds.sh
git commit -m "docs: remove --instrument; document MGF-only --fragmentation/--fragment-tol params"
```

---

## Final verification (whole-feature gate)

- [ ] **Step 1: Workspace build + clippy**

Run: `cargo build --workspace 2>&1 | tail -5 && cargo clippy --workspace -- -D warnings 2>&1 | tail -10`
Expected: clean.

- [ ] **Step 2: Full test suite**

Run: `cargo test --workspace 2>&1 | tail -30`
Expected: only the pre-existing known failures from the plan header; no new failures.

- [ ] **Step 3: Byte-identical regression confirmation (no override path)**

Confirm the PIN/score golden parity tests pass (these prove metadata-bearing selection + features are unchanged): `cargo test --workspace parity 2>&1 | tail -15` and `cargo test --workspace golden 2>&1 | tail -15`.

- [ ] **Step 4: Hand back for the VM benchmark gate**

The project requires Astral + UPS1 + a05058 TMT through the same Percolator for any feature change. The refactor is byte-identical with no override, so model selection on those (all metadata-bearing) datasets should be unchanged — but state this explicitly to the user and let them run the VM benchmark as the final merge gate. Do NOT claim a benchmark pass without the VM run.

---

## Self-review notes (author)

- **Spec coverage:** Decisions A (remove `--instrument`, Task 4), B (`--fragmentation` hidden MGF-only, Task 4), C (`--fragment-tol-*`, Task 3), D (no inference — nothing to build), E (CID/LowRes default + warning, Task 4 Step 4), F (README/DOCS, Task 5). §3 centralization = Tasks 1-2. All covered.
- **Byte-identical guard:** Tasks 1-2 keep override `None` ⇒ identical tolerance; Task 2 Step 5 and Final Step 3 enforce it via golden/parity tests.
- **Type consistency:** `feature_match_tolerance()`/`set_fragment_tol_override()` (Task 1) used verbatim in Tasks 2-3. `Tolerance::{Ppm,Da,as_da,raw_value}` exist (`crates/model/src/tolerance.rs`). `resolve_metadataless_selection`/`cli_fragmentation_to_activation_opt`/`cli_fragment_tol_override` all defined before use.
- **Known risk:** exact line numbers drift; every step gives the surrounding code to match. The `compute_psm_features` `feature_tol` downstream arithmetic (Task 2 Step 3) must be eyeballed — flagged in-step.
