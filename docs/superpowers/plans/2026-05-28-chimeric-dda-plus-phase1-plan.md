# Chimeric DDA+ — Phase 1 Implementation Plan (full-window search + multi-PSM emission)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.
>
> Phase 1 of the design `docs/superpowers/specs/2026-05-28-chimeric-dda-plus-integration-design.md`. Option **A** (build on the existing `bucket_index` scan; fragment-index deferred). Phases 2–3 (MS1 XIC refinement, shared-fragment rescoring) are separate plans.

**Goal:** Add a `--chimeric` mode that searches each MS2 against the full isolation window and emits the top-N distinct-peptide PSMs per scan. Default off ⇒ output bit-identical to current.

**Architecture:** Reuse the existing per-spectrum candidate **mass-range scan** (`bucket_index.range`) and the **per-charge `TopNQueue`** machinery. The only new data is the isolation-window width (parsed from mzML); the only behavioral switch is "candidate window = isolation window" + "retain top-N distinct peptides" when `--chimeric` is set.

**Tech Stack:** Rust workspace (`model`, `input`, `search`, `output`, `msgf-rust`); quick-xml mzML parser; Percolator PIN output.

---

## Conventions

- Build: `cargo build --release -p msgf-rust`
- Targeted tests: `cargo test --release -p <crate> <test_name>`
- **Bit-identical gate (the load-bearing safety net):** `cargo test --release -p msgf-rust precursor_cal_off_pin_tsv_match_golden_after_sort` → `ok. 1 passed`. Must stay green for EVERY task (chimeric defaults off).
- Clippy: `cargo clippy --workspace --all-targets -- -D warnings` → clean.
- Commit after each task.

---

### Task 1: Isolation-window offsets on the `Spectrum` model

**Files:**
- Modify: `crates/model/src/spectrum.rs` (`Spectrum` struct + constructors)

- [ ] **Step 1: Failing test** — add to `spectrum.rs` tests:
```rust
#[test]
fn spectrum_isolation_offsets_default_none() {
    let s = Spectrum::empty();
    assert_eq!(s.isolation_lower_offset, None);
    assert_eq!(s.isolation_upper_offset, None);
}
```
Run: `cargo test --release -p model spectrum_isolation_offsets_default_none` → FAIL (no field).

- [ ] **Step 2: Add fields.** In `struct Spectrum`, after `precursor_charge`:
```rust
/// Isolation-window lower offset in Da (selected m/z − lower = window start).
/// `None` when the mzML omits `<isolationWindow>`. Used only by `--chimeric`.
pub isolation_lower_offset: Option<f64>,
/// Isolation-window upper offset in Da (selected m/z + upper = window end).
pub isolation_upper_offset: Option<f64>,
```
Set both to `None` in every `Spectrum` constructor/literal in this file (the `empty()`/test builders). Use `grep -n "precursor_charge:" crates/model/src/spectrum.rs` to find them all.

- [ ] **Step 3:** `cargo test --release -p model spectrum_isolation_offsets_default_none` → PASS. Then `cargo build --release -p msgf-rust` and fix any other `Spectrum { .. }` literals across crates (`grep -rn "Spectrum {" crates/ --include=*.rs`) by adding `isolation_lower_offset: None, isolation_upper_offset: None,` (or `..Default::default()` if the struct derives Default).

- [ ] **Step 4: Commit** `feat(model): isolation-window offset fields on Spectrum (chimeric prep)`

---

### Task 2: Parse `<isolationWindow>` in the mzML reader

**Files:**
- Modify: `crates/input/src/mzml.rs` (State enum, SpectrumBuilder, cvParam handling, builder→Spectrum)
- Test: `crates/input/src/mzml.rs` inline tests (there are existing mzML-snippet parse tests — follow their pattern)

CV accessions: `MS:1000827` isolation window target m/z, `MS:1000828` isolation window lower offset, `MS:1000829` isolation window upper offset.

- [ ] **Step 1: Failing test** — add an mzML snippet test mirroring the existing parser tests, with a `<precursor><isolationWindow>` block:
```rust
#[test]
fn parses_isolation_window_offsets() {
    let xml = r#"<spectrum id="scan=1" defaultArrayLength="0">
      <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2"/>
      <precursorList count="1"><precursor>
        <isolationWindow>
          <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="500.0"/>
          <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.5"/>
          <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="1.5"/>
        </isolationWindow>
        <selectedIonList count="1"><selectedIon>
          <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="500.5"/>
        </selectedIon></selectedIonList>
      </precursor></precursorList>
      <binaryDataArrayList count="0"></binaryDataArrayList>
    </spectrum>"#;
    let spec = parse_single_spectrum_for_test(xml); // use the existing test helper; if none, parse via MzMLReader over a cursor
    assert_eq!(spec.isolation_lower_offset, Some(1.5));
    assert_eq!(spec.isolation_upper_offset, Some(1.5));
}
```
(If no single-spectrum test helper exists, follow whatever harness the existing `mzml.rs` tests use to feed a string through the parser.) Run → FAIL.

- [ ] **Step 2: Add parse state + capture.**
  - Add `IsolationWindow` to the parser `State` enum.
  - On `Event::Start`/`Empty` for `b"isolationWindow"` while `state == State::Spectrum`, set `state = State::IsolationWindow`; on its end, pop back to `State::Spectrum`. (Mirror how `selectedIon`/`activation` are handled around mzml.rs:463-473.)
  - Add fields to `SpectrumBuilder`: `isolation_lower_offset: Option<f64>`, `isolation_upper_offset: Option<f64>`.
  - In `apply_cv_param` (where existing cvParams are matched), when `state == State::IsolationWindow`, match `MS:1000828` → `isolation_lower_offset`, `MS:1000829` → `isolation_upper_offset` (parse value as f64).
  - In the builder→`Spectrum` finalize, copy both fields through.

- [ ] **Step 3:** Run the test → PASS. Run the full `input` crate tests (`cargo test --release -p input`) to confirm no parse regression.

- [ ] **Step 4: Commit** `feat(mzml): capture isolation-window lower/upper offsets`

---

### Task 3: `--chimeric` flag, isolation-width fallback, chimeric top-N

**Files:**
- Modify: `crates/msgf-rust/src/bin/msgf-rust.rs` (CLI), `crates/search/src/search_params.rs` (`SearchParams`)

- [ ] **Step 1:** In `SearchParams`, add:
```rust
/// Full-isolation-window chimeric search (MSFragger-DDA+ style). Default false.
pub chimeric: bool,
/// Fallback isolation half-width (Da) used when the mzML lacks
/// `<isolationWindow>` offsets. Only consulted when `chimeric` is true.
pub chimeric_isolation_halfwidth_da: f64,
```
Set defaults in `default_tryptic` (and any other constructor): `chimeric: false, chimeric_isolation_halfwidth_da: 1.5`.

- [ ] **Step 2:** In `msgf-rust.rs` `Cli`, add:
```rust
/// Search the full isolation window per MS2 and emit multiple distinct-peptide
/// PSMs per scan (chimeric / co-fragmented peptides; MSFragger-DDA+ style).
#[arg(long, default_value = "false")]
chimeric: bool,
/// Fallback isolation half-width in Da when the mzML omits isolation offsets.
#[arg(long, default_value = "1.5")]
isolation_halfwidth: f64,
```
Wire `params.chimeric = cli.chimeric; params.chimeric_isolation_halfwidth_da = cli.isolation_halfwidth;` where `SearchParams` is assembled. When `--chimeric` and `--top-n` is still its default of 1, bump the effective top-N (e.g. `if cli.chimeric && top_n == 1 { top_n = 5 }`) and log it.

- [ ] **Step 3:** `cargo build --release -p msgf-rust`; `./target/release/msgf-rust --help | grep -E "chimeric|isolation-halfwidth"` shows both. Bit-identical gate green (flag off ⇒ unchanged).

- [ ] **Step 4: Commit** `feat(cli): --chimeric + --isolation-halfwidth flags (wiring only)`

---

### Task 4: Widen candidate enumeration to the isolation window when chimeric

**Files:**
- Modify: `crates/search/src/match_engine.rs` (`run_chunk_inner`, the per-spectrum candidate-window derivation ~L243-266)

Today the candidate window is `nominal_center − iso_max − widen_right .. nominal_center − iso_min + widen_left`, derived from `precursor_mz` ± precursor tolerance, per charge. For chimeric, the window must instead span the **isolation window** in m/z, converted to neutral nominal mass per charge.

- [ ] **Step 1: Failing test** — add a `match_engine` test: a synthetic spectrum with `isolation_lower_offset = Some(2.0)`, `isolation_upper_offset = Some(2.0)`, two candidate peptides whose neutral masses correspond to precursors at `selected_mz` and `selected_mz − 1.2` (inside the window but OUTSIDE ±tol). With `params.chimeric = true` and a charge, assert BOTH candidate indices appear in the enumerated window set; with `chimeric = false`, assert only the on-precursor one does. (Factor the window-derivation into a testable helper if needed.)
Run → FAIL.

- [ ] **Step 2: Implement.** In the per-charge window loop, when `params.chimeric`:
  - Compute the isolation window m/z bounds:
    `lo_mz = selected_mz − isolation_lower_offset.unwrap_or(halfwidth)`,
    `hi_mz = selected_mz + isolation_upper_offset.unwrap_or(halfwidth)`.
  - Convert each bound to neutral nominal mass at charge `z`:
    `neutral = (mz − PROTON) * z − H2O`; `nominal = nominal_from(adjusted_observed_neutral_mass(neutral, shift_ppm))`.
  - Set `min_nominal = nominal_from(lo) − iso_max`, `max_nominal = nominal_from(hi) − iso_min` (keep the isotope-error widening; drop the ±precursor-tol widening since the window already spans wider).
  - Feed these bounds to the existing `bucket_index.range(min_nominal..=max_nominal)` extend.
  - When `!params.chimeric`: unchanged code path (bit-identical).

- [ ] **Step 3:** Test → PASS. Bit-identical gate green (chimeric off). `cargo test --release -p search` green.

- [ ] **Step 4: Commit** `feat(search): widen candidate window to isolation window under --chimeric`

---

### Task 5: Retain + emit top-N distinct peptides per scan

**Files:**
- Modify: `crates/search/src/match_engine.rs` (the per-charge-queue merge + `matches_precursor` gating), `crates/output/src/pin.rs` (SpecId uniqueness)

The widened window (Task 4) now feeds many candidates into the per-charge `TopNQueue`s. `dedup_pepseq_score` already collapses one-peptide-many-proteins; `TopNQueue` keeps top-N by SpecE. With effective top-N>1 (Task 3) the spectrum queue already holds multiple distinct peptides. Two gaps to close:

- **5a — `matches_precursor` gating.** In the candidate loop, scoring still calls `matches_precursor(spec, peptide, z, offset, tol, shift)`. Under chimeric, a co-isolated peptide's mass is offset from the *selected* precursor by up to the window half-width — it would fail the ±tol precursor check and be dropped. When `params.chimeric`, replace the precursor-match gate with an **isolation-window membership** check: accept the candidate if its theoretical neutral mass (for some isotope offset) falls within `[lo_mz, hi_mz]` converted to neutral mass at charge `z` (± the fragment-irrelevant precursor tol for the *matched* peptide's own m/z). Compute `mass_error_ppm` against the candidate's own nearest in-window precursor m/z (for the PIN `dm`/`isotope_error` columns), not the selected m/z.

- [ ] **Step 1: Failing test** — synthetic chimeric spectrum where two distinct peptides both have fragment support; with `chimeric=true, top_n=5`, assert the final per-spectrum queue contains BOTH peptides (distinct sequences). With `chimeric=false`, assert only the on-precursor peptide. Run → FAIL.

- [ ] **Step 2: Implement 5a** (the gating swap above), guarded by `if params.chimeric { … } else { existing matches_precursor }`.

- [ ] **Step 3: SpecId uniqueness (5b).** In `pin.rs::write_psm_row`, when more than one PSM is emitted for a scan, ensure `SpecId` is unique per row. Inspect the current SpecId construction; if it is `<basename>_<scan>_<scan>_<charge>` (or similar) it may already collide across distinct peptides at the same scan/charge. Append a per-scan emission rank (`_<k>`) so rows are unique. Add/extend a `pin.rs` test asserting unique SpecIds for two PSMs sharing a scan.

- [ ] **Step 4:** Tests → PASS. Bit-identical gate green (chimeric off path untouched). `cargo test --release -p search -p output` green. Clippy clean.

- [ ] **Step 5: Commit** `feat(search/output): emit top-N distinct-peptide PSMs per scan under --chimeric`

---

### Task 6: Workspace gate + VM bench (the decision point)

**Files:** none (validation)

- [ ] **Step 1:** Full workspace tests under the CI skip list (see PR #40 conventions) → all green. Clippy clean. Bit-identical gate green.

- [ ] **Step 2: Bit-identical off-mode bench (must hold).** On the VM, build and run all 3 datasets WITHOUT `--chimeric`; sorted-row PIN diff vs the current master baseline → identical. (Proves the feature is a true no-op when off.)

- [ ] **Step 3: Chimeric-on bench (the measurement).** Re-run PXD001819 + TMT WITH `--chimeric` (+ Astral as a regression guard). Capture: PSMs @1% FDR (Percolator, `--only-psms`), wall time, target/decoy balance, mean PSMs-per-scan. Compare PXD/TMT @1% vs the current baseline (PXD 14,808, TMT 9,605).

- [ ] **Step 4: Decide (per the spec's open question).** If PXD/TMT gain meaningfully and wall is acceptable → ship Phase 1, proceed to Phase 2. If wall is unacceptable → the fragment-index enabler becomes a prerequisite (spec Cross-cutting decision B); pause Phase 1 ship and plan the index. If no PSM gain → investigate (window too wide/narrow? emission collapsing?) before shipping.

- [ ] **Step 5: Open the Phase 1 PR** to `dev` with the bench table (off=bit-identical, on=PXD/TMT/Astral PSM + wall), and the ship/defer decision.

---

## Self-review notes

- **Spec coverage:** spec Phase 1a→Task 2; 1b→Task 1; 1c→Task 4; 1d→Task 5; 1e→Task 5b; gate→Task 6. Covered.
- **No-op-when-off is the invariant:** every task keeps the `chimeric=false` path on the existing code, gated by `if params.chimeric`, and the bit-identical golden test guards it at each step.
- **Deliberately not pre-written:** Task 2's exact cvParam-match site and Task 5's `matches_precursor`/SpecId edits depend on current signatures the implementer must read (`apply_cv_param`, the candidate loop, `write_psm_row`); each task pins the file + function + behavior + test rather than a fabricated diff, since those internals weren't fully read at plan time.
- **Type consistency:** `isolation_lower_offset`/`isolation_upper_offset` (`Option<f64>`), `chimeric` (`bool`), `chimeric_isolation_halfwidth_da` (`f64`) used consistently across Tasks 1–5.
