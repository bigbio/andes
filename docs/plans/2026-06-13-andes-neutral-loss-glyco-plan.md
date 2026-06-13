# Neutral-Loss-Aware Scoring (glyco / Unimod 393) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Let andes score labile-mod glycopeptides by modeling user-declared neutral losses (stepwise b/y loss ions) via a pooled loss-flag `IonType`, plus capture/emit the Unimod accession for quantms interop.

**Architecture:** A modification can declare neutral losses (mods.txt `loss=` attribute). Fragments spanning a loss-bearing residue additionally predict loss-shifted b/y ions tagged with a `loss: true` `IonType` discriminant; the GF model carries a *pooled* loss-ion distribution family (one per partition/series/charge). Inert and byte-identical when no mod declares a loss. The `accession=` attribute is captured and emitted in a new additive PIN/TSV column.

**Tech Stack:** Rust workspace (`model`, `scoring`, `model-train`, `output`, `andes`), parquet model store, clap CLI. TDD + per-task commits.

**Spec:** `docs/specs/2026-06-13-andes-neutral-loss-glyco-design.md`

**Branch:** decide at execution (recommend `feat/glyco-neutral-loss` off current HEAD so it inherits recent scoring work). Build gate: `cargo build --workspace` + `cargo clippy --workspace -- -D warnings` clean.

**Pre-existing failing tests (NOT ours):** `train_from_msnet::fragment_tolerance_override_changes_model`, 3 `match_engine_smoke` min_peaks tests, `chemistry_constants::integer_mass_scaler_matches_residue_table_mean`. Confirm any new failure is ours via `git stash`.

**Scope of THIS plan:** SP1 (representation + config + accession) and SP2 (model format + scoring). SP3 (training) and SP4 (benchmark validation) are outlined at the end — they are gated on the glyco dataset (user-provided) + the VM and become a follow-on plan.

---

## File structure

| File | Responsibility | Phase |
|---|---|---|
| `crates/model/src/modification.rs` | `neutral_losses` field; `key=value` grammar (`loss=`, `accession=`); errors | SP1 |
| `crates/scoring/src/param_model.rs` | `IonType` loss discriminant + accessors | SP1 |
| `crates/scoring/src/scoring/fragment_ions.rs` | predict loss-shifted b/y ions | SP1 |
| `crates/output/src/pin.rs`, `tsv.rs` | additive `Modifications` (accession) column | SP1 |
| `DOCS.md` §2 + new glyco template under `resources/` | config docs + example | SP1 |
| `crates/model-train/src/store/{write,read}.rs` | serialize/deserialize the loss bit (back-compat) | SP2 |
| `crates/scoring/src/scoring/rank_scorer.rs` | score loss ions via pooled table; absent → missing | SP2 |

---

## Task 1: `Modification.neutral_losses` + `key=value` grammar

**Files:** Modify `crates/model/src/modification.rs` (struct ~L29-37; `ModParseError` ~L55-68; `from_mods_txt_line` ~L74-118; tests at bottom).

- [ ] **Step 1: Write failing tests** (add to `#[cfg(test)] mod tests`):

```rust
#[test]
fn parses_loss_and_accession_attributes() {
    let m = Modification::from_mods_txt_line(
        "340.100562,K,opt,any,Glucosylgalactosyl,loss=162.0528;324.1056,accession=UNIMOD:393"
    ).unwrap();
    assert_eq!(m.residue, ResidueSpec::Specific(b'K'));
    assert!(!m.fixed);
    assert_eq!(m.neutral_losses, vec![162.0528, 324.1056]);
    assert_eq!(m.accession.as_deref(), Some("UNIMOD:393"));
}

#[test]
fn five_field_line_has_no_losses_or_accession() {
    let m = Modification::from_mods_txt_line("57.02146,C,fix,any,Carbamidomethyl").unwrap();
    assert!(m.neutral_losses.is_empty());
    assert_eq!(m.accession, None);
}

#[test]
fn rejects_unknown_attr_and_bad_loss() {
    assert!(matches!(
        Modification::from_mods_txt_line("1.0,K,opt,any,X,frobnicate=7"),
        Err(ModParseError::UnknownModAttr { .. })));
    assert!(matches!(
        Modification::from_mods_txt_line("1.0,K,opt,any,X,loss=abc"),
        Err(ModParseError::BadNeutralLoss { .. })));
    assert!(matches!(
        Modification::from_mods_txt_line("1.0,K,opt,any,X,nokey"),
        Err(ModParseError::BadModAttr { .. })));
}
```

- [ ] **Step 2: Run to verify fail**

Run: `cargo test -p model parses_loss_and_accession 2>&1 | tail -15`
Expected: FAIL (`neutral_losses` field missing).

- [ ] **Step 3: Add the field** to `Modification` (after `accession`):

```rust
    /// User-declared neutral-loss masses (Da) for this mod's fragment ions.
    /// Empty ⇒ no loss ions predicted (default; byte-identical to pre-feature).
    pub neutral_losses: Vec<f64>,
```

- [ ] **Step 4: Add error variants** to `ModParseError`:

```rust
    #[error("unknown mod attribute key {key:?} (expected loss|accession)")]
    UnknownModAttr { key: String },
    #[error("malformed mod attribute {field:?} (expected key=value)")]
    BadModAttr { field: String },
    #[error("invalid neutral-loss value {value:?} (expected positive number < 2000)")]
    BadNeutralLoss { value: String },
```

- [ ] **Step 5: Rewrite the parser tail.** Replace the `splitn(5)` + `len()!=5` guard and the `Ok(Modification { … })` construction so the first 5 fields are positional and any remaining comma fields are `key=value`:

```rust
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 5 {
            return Err(ModParseError::WrongFieldCount { got: parts.len() });
        }
        let (mass_s, residues_s, fixity_s, location_s, name_s) =
            (parts[0].trim(), parts[1].trim(), parts[2].trim(), parts[3].trim(), parts[4].trim());
        // (mass/residue/fixity/location parsing unchanged — keep existing blocks)

        let mut neutral_losses: Vec<f64> = Vec::new();
        let mut accession: Option<String> = None;
        for attr in &parts[5..] {
            let attr = attr.trim();
            if attr.is_empty() { continue; }
            let (key, value) = attr.split_once('=')
                .ok_or_else(|| ModParseError::BadModAttr { field: attr.to_string() })?;
            match key.trim().to_ascii_lowercase().as_str() {
                "loss" => {
                    for tok in value.split(';') {
                        let tok = tok.trim();
                        if tok.is_empty() { continue; }
                        let v: f64 = tok.parse()
                            .map_err(|_| ModParseError::BadNeutralLoss { value: tok.to_string() })?;
                        if !(v > 0.0 && v < 2000.0) {
                            return Err(ModParseError::BadNeutralLoss { value: tok.to_string() });
                        }
                        neutral_losses.push(v);
                    }
                }
                "accession" => accession = Some(value.trim().to_string()),
                other => return Err(ModParseError::UnknownModAttr { key: other.to_string() }),
            }
        }

        Ok(Modification {
            name: name_s.to_string(),
            mass_delta, residue, location, fixed,
            accession,
            neutral_losses,
        })
```

(Keep the existing mass/residue/fixity/location parsing blocks verbatim between the field binding and the attribute loop.)

- [ ] **Step 6: Fix other `Modification { … }` constructors** the compiler flags (test helpers in this file, `aa_set.rs:199`, `amino_acid.rs` tests, and any others): add `neutral_losses: Vec::new(),`. Run `cargo build -p model 2>&1 | tail -20` and fix each.

- [ ] **Step 7: Run tests + commit**

Run: `cargo test -p model 2>&1 | tail -15` → new tests PASS, existing pass.
```bash
git add crates/model/src/modification.rs crates/model/src/aa_set.rs crates/model/src/amino_acid.rs
git commit -m "feat(model): mods.txt key=value tail — loss= + accession= attributes"
```

---

## Task 2: `IonType` loss discriminant (pooled key)

**Files:** Modify `crates/scoring/src/param_model.rs` (`IonType` enum ~L494-499; accessors `offset`/`charge`/`is_prefix` ~L500-525). Ripple: ~58 construction sites across the workspace (compiler-driven).

- [ ] **Step 1: Write failing test** (in param_model.rs tests):

```rust
#[test]
fn loss_ion_is_distinct_key_from_intact() {
    use std::collections::HashMap;
    let intact = IonType::Prefix { charge: 1, offset_bits: 1.0f32.to_bits(), loss: false };
    let lost   = IonType::Prefix { charge: 1, offset_bits: 1.0f32.to_bits(), loss: true };
    let mut m = HashMap::new();
    m.insert(intact, "a"); m.insert(lost, "b");
    assert_eq!(m.len(), 2);
    assert!(!intact.is_loss());
    assert!(lost.is_loss());
}
```

- [ ] **Step 2: Run to verify fail** — `cargo test -p scoring loss_ion_is_distinct 2>&1 | tail` → FAIL (no `loss` field / `is_loss`).

- [ ] **Step 3: Add the field + accessor.** Change the variants:

```rust
    Prefix { charge: i32, offset_bits: u32, loss: bool },
    Suffix { charge: i32, offset_bits: u32, loss: bool },
    Noise,
```

Add accessor in `impl IonType`:

```rust
    /// True if this is a neutral-loss-shifted fragment ion (pooled loss key).
    pub fn is_loss(&self) -> bool {
        matches!(self, IonType::Prefix { loss: true, .. } | IonType::Suffix { loss: true, .. })
    }
```

Update `offset()`/`charge()`/`is_prefix()`/`is_suffix()` patterns to add `, ..` if they don't already use it (the `{ offset_bits, .. }` form is forward-compatible).

- [ ] **Step 4: Compiler-driven ripple.** Run `cargo build --workspace 2>&1 | tee /tmp/build.txt | tail -30`. For EVERY flagged `IonType::Prefix { … }` / `IonType::Suffix { … }` **construction** literal (≈58 across `model-train/src/{counts,store/read}.rs`, `scoring/src/{param_model,testutil,scoring/*}.rs`, `search/src/{match_engine,coisolation}.rs`, tests, examples), add `loss: false`. **Do NOT** touch `{ .. }` match patterns (forward-compatible). Loss-`true` is produced ONLY in Task 3. Re-run build until clean.

- [ ] **Step 5: Run + commit**

Run: `cargo test -p scoring 2>&1 | tail -15` (existing pass — intact keys unchanged, `loss:false` ≡ old).
```bash
git add -A
git commit -m "feat(scoring): IonType pooled loss discriminant (loss:false ≡ prior behavior)"
```

---

## Task 3: Predict loss-shifted b/y ions

**Files:** Modify `crates/scoring/src/scoring/fragment_ions.rs` (`predict_by_ions` ~L98-144; module header L1-5). Inert when no residue carries losses.

- [ ] **Step 1: Write failing test:**

```rust
#[test]
fn emits_loss_ions_for_loss_bearing_residue() {
    use model::modification::Modification;
    // Build a 3-residue peptide with a +340 mod carrying losses [162.0528, 324.1056] on residue 2.
    let pep = peptide_with_loss_mod(vec![162.0528, 324.1056]); // helper: K at pos 1 modified
    let ions = predict_by_ions(&pep, 1..=1);
    // Intact b/y still present; plus loss-shifted ions tagged loss:true.
    assert!(ions.iter().any(|p| !p.ion_type.is_loss()));
    let loss_ions: Vec<_> = ions.iter().filter(|p| p.ion_type.is_loss()).collect();
    assert!(!loss_ions.is_empty(), "expected loss-shifted ions");
    // A loss ion sits exactly 162.0528 (or 324.1056) below some intact ion of the same series/charge.
    assert!(loss_ions.iter().any(|l| ions.iter().any(|i|
        !i.ion_type.is_loss() && (i.mz - l.mz - 162.0528).abs() < 1e-4)));
}

#[test]
fn no_loss_ions_when_no_loss_mod() {
    let pep = plain_peptide(); // existing helper / construct unmodified
    let ions = predict_by_ions(&pep, 1..=1);
    assert!(ions.iter().all(|p| !p.ion_type.is_loss()));
}
```

(Add `peptide_with_loss_mod` near the existing test helpers; it builds a `Peptide` whose modified residue's `mod_.neutral_losses` is the passed vec.)

- [ ] **Step 2: Run to verify fail** — `cargo test -p scoring emits_loss_ions 2>&1 | tail` → FAIL.

- [ ] **Step 3: Implement.** In `predict_by_ions`, while walking fragment positions, track whether the current prefix/suffix span includes a residue whose `mod_` has non-empty `neutral_losses`; collect the union of those loss masses for the span. After pushing the intact ion `(ion_type, mz)` at charge `z`, for each loss `L` in the span's loss set, also push a loss ion:

```rust
        // After the intact ion is pushed for this (series, position, charge z):
        if !span_losses.is_empty() {
            let base_offset_bits = /* same offset_bits as the intact ion */;
            for &loss in &span_losses {
                let loss_mz = mz - loss / z as f64;
                let loss_ion_type = match intact_ion_type {
                    IonType::Prefix { charge, offset_bits, .. } =>
                        IonType::Prefix { charge, offset_bits, loss: true },
                    IonType::Suffix { charge, offset_bits, .. } =>
                        IonType::Suffix { charge, offset_bits, loss: true },
                    IonType::Noise => continue,
                };
                out.push(PredictedIon { ion_type: loss_ion_type, mz: loss_mz /* + other fields copied from intact */ });
            }
        }
```

`span_losses` is built incrementally: for a prefix ion ending at residue `i`, it is the union of `neutral_losses` over residues `0..=i` that carry a mod with losses (v1: union, no cross-products). For suffix ions, the mirror span. Update the module header L3 from "no neutral losses" to "b/y plus user-declared neutral-loss ions (pooled loss IonType)."

- [ ] **Step 4: Run to verify pass** — `cargo test -p scoring fragment 2>&1 | tail` → PASS, existing fragment tests unchanged.

- [ ] **Step 5: Commit**
```bash
git add crates/scoring/src/scoring/fragment_ions.rs
git commit -m "feat(scoring): predict neutral-loss-shifted b/y ions for loss-bearing residues (inert when none)"
```

---

## Task 4: Additive `Modifications` (accession) output column

**Files:** Modify `crates/output/src/pin.rs` (header column list ~L158-256; row writer ~L514) and `crates/output/src/tsv.rs` (analogous). Additive — the existing `Peptide` column is untouched.

- [ ] **Step 1: Write failing test** (pin.rs tests): assert the header contains `"Modifications"` (appended, after existing columns) and that a peptide whose modified residue has `accession = Some("UNIMOD:393")` emits `6:UNIMOD:393` (1-based residue position) in that column; a peptide with no accession emits an empty field.

```rust
#[test]
fn pin_emits_additive_modifications_accession_column() {
    let header = pin_header(); // existing header-builder used by tests
    assert_eq!(header.last().map(String::as_str), Some("Modifications"));
    // build a PSM whose residue 6 carries accession UNIMOD:393 → row's last col == "6:UNIMOD:393"
    // (mirror the existing pin row test harness)
}
```

- [ ] **Step 2: Run to verify fail** — `cargo test -p output pin_emits_additive_modifications 2>&1 | tail` → FAIL.

- [ ] **Step 3: Implement.** Append `"Modifications".to_string()` to the PIN header `cols` vector (LAST, after every existing column — additive). In the row writer, after the existing columns, write a tab + a `;`-joined list of `"{pos}:{curie}"` for each modified residue whose `mod_.accession` is `Some`, 1-based position; empty string if none. Do the same in `tsv.rs`. Pull the accession from `cand.peptide.residues[i].mod_.as_ref().and_then(|m| m.accession.clone())`.

- [ ] **Step 4: Run + confirm byte-identical existing columns** — `cargo test -p output 2>&1 | tail -15`. The PIN-schema parity test will now expect the extra column; update its expected column count (the column is additive/appended — confirm it's the LAST column). Verify the `Peptide` column value is unchanged.

- [ ] **Step 5: Commit**
```bash
git add crates/output/src/pin.rs crates/output/src/tsv.rs
git commit -m "feat(output): additive Modifications column emitting Unimod accessions (Peptide column unchanged)"
```

---

## Task 5: Docs + glyco template

**Files:** Modify `DOCS.md` §2 (Mods.txt format). Create `resources/mods/glyco_example.txt` (or alongside existing example mods).

- [ ] **Step 1:** In DOCS §2, document the `key=value` tail: `loss=<m1;m2;…>` and `accession=<CURIE>`; the 5-field-unchanged rule; the comma-in-name caveat; the validation errors. Add a **common neutral-loss table**: Hex `162.0528`, HexNAc `203.0794`, NeuAc `291.0954`, phospho (H₃PO₄) `97.9769`, sulfo `79.9568`. Add the copy-paste line:
  `340.100562,K,opt,any,Glucosylgalactosyl,loss=162.0528;324.1056,accession=UNIMOD:393`
  and a one-paragraph "how losses are scored" note (pooled loss IonType; needs a glyco model — SP3). State attributes are orthogonal to `fix|opt`/`location`/`NumMods`.
- [ ] **Step 2:** Create `resources/mods/glyco_example.txt` with the Carbamidomethyl-C fixed line + the Glucosylgalactosyl variable line above + a `NumMods=2` header, as a ready template.
- [ ] **Step 3: Commit**
```bash
git add DOCS.md resources/mods/glyco_example.txt
git commit -m "docs: mods.txt loss=/accession= attributes, common-loss table, glyco template"
```

---

## Task 6: Model-store serialization of the loss bit (back-compat)

**Files:** Modify `crates/model-train/src/store/write.rs` (ion-type encoding) and `crates/model-train/src/store/read.rs` (~L350-365 decode).

- [ ] **Step 1: Write failing round-trip test** (model-train tests): build a `Param` containing a `loss:true` IonType entry in `rank_dist_table`, write it to an in-memory/temp parquet store, read it back, assert the `loss:true` entry survives; and assert a store written WITHOUT the loss column still reads (loss defaults `false`).

- [ ] **Step 2: Run to verify fail.**

- [ ] **Step 3: Implement.** The flat ion-type encoding currently stores `(is_prefix_f, ion_charge, offset_bits_f, …)` and the reader reconstructs at [read.rs:352-362](../../crates/model-train/src/store/read.rs). Add a loss flag to the per-ion-type encoding (an extra trailing value, e.g. `loss_f` 1.0/0.0). Reader: when the loss value is present read it; when absent (older files) default `false`. Construct `IonType::Prefix { charge, offset_bits, loss: loss_f > 0.5 }` (and Suffix). Keep `is_prefix_f` semantics (>0.5 Prefix, <-0.5 Noise, else Suffix).

- [ ] **Step 4: Run round-trip + the bundled-store load test** — `cargo test -p model-train 2>&1 | tail -15` and `cargo test --workspace param_loads_all_bundled 2>&1 | tail` (the existing 39 bundled models must still load — they have no loss entries → all intact, unchanged).

- [ ] **Step 5: Commit**
```bash
git add crates/model-train/src/store/write.rs crates/model-train/src/store/read.rs
git commit -m "feat(model-train): serialize IonType loss bit in parquet store (back-compat default false)"
```

---

## Task 7: Score loss ions via the pooled table

**Files:** Modify `crates/scoring/src/scoring/rank_scorer.rs` (table build ~L55-90; node-score lookup ~L169-195). Byte-identical for peptides with no loss ions / models with no loss tables.

- [ ] **Step 1: Write failing tests:** (a) a `RankScorer` built from a `Param` that HAS a `loss:true` rank table returns that table's LLR for a loss ion at a given rank; (b) a scorer built from a Param WITHOUT loss tables returns the missing-ion score for a loss ion (no panic); (c) for a no-loss peptide, the per-ion scores are byte-identical to a pre-change baseline (snapshot a few node scores).

- [ ] **Step 2: Run to verify fail.**

- [ ] **Step 3: Implement.** `RankScorer::new` already builds `log_table`/`partition_ion_logs` from `param.rank_dist_table` keyed by `IonType` — a `loss:true` key flows through automatically (no special-casing needed for table build). In the node-score lookup, a `loss:true` ion whose `(partition, ion_type)` is absent must fall to the missing-ion/absent path exactly as an unknown intact ion does today (verify the existing lookup already returns the absent score for a missing key; if it panics/unwraps, guard it). Ensure `partition_ion_types_cache`/`ion_types_for_partition_slice` include `loss:true` types when present so the DP iterates them.

- [ ] **Step 4: Run + byte-identical guard** — `cargo test -p scoring 2>&1 | tail` and `cargo test --workspace score_psm 2>&1 | tail` (no-loss scoring unchanged).

- [ ] **Step 5: Commit**
```bash
git add crates/scoring/src/scoring/rank_scorer.rs
git commit -m "feat(scoring): score pooled loss ions; absent loss table → missing-ion (no-loss byte-identical)"
```

---

## Final SP1+SP2 verification

- [ ] `cargo build --workspace` + `cargo clippy --workspace -- -D warnings` clean.
- [ ] `cargo test --workspace` — only the known pre-existing failures; confirm `score_psm_*_parity`, `output_pin_schema_parity` (column-count updated for the additive column), `param_loads_all_bundled`, `migration_parity` all pass.
- [ ] **Byte-identical guard:** a standard search (no loss mods, bundled non-glyco model) predicts no loss ions and the existing PIN columns are unchanged (the only diff is the appended `Modifications` column, empty for non-accession mods). State this to the user.
- [ ] At this point loss ions are predicted, representable, serializable, and *scored when a model has loss tables* — but no bundled model has them yet. That is SP3.

---

## SP3 + SP4 — Outline (follow-on plan, gated on the glyco dataset + VM)

**Do not start until the user provides the glyco training/benchmark dataset.** Then write a dedicated plan:

- **SP3 — Training:** in `crates/model-train/src/estimate.rs` / `counts.rs`, the accumulator already keys facts by `IonType`; feed it confident glyco PSMs (whose mods.txt declares losses) so predicted `loss:true` ions accumulate into the pooled loss rank/error/existence tables. Add a `glyco` model build (catalog slug exists, `catalog.rs:95`) producing a model WITH loss tables; wire into the store. TDD: synthetic glyco PSM → estimator produces a non-empty pooled loss distribution.
- **SP4 — Validation (benchmark-gated, user runs VM):** search the glyco dataset with the glyco model + `resources/mods/glyco_example.txt`, FDP/entrapment-controlled → demonstrate PSM gain vs the no-loss baseline. Confirm **byte-identical** on Astral/UPS1/a05058 TMT (no loss mods). Deliverable handed to the user for the VM run; do NOT claim a benchmark pass without it.

---

## Self-review notes (author)

- **Spec coverage:** SP1 decisions A/C/D/F/G → Tasks 1 (grammar+accession capture), 2 (pooled IonType), 3 (prediction), 4 (accession output), 5 (docs). SP2 → Tasks 6 (store), 7 (scoring). SP3/SP4 → outline. Decision B (data available) is an SP3 input.
- **Byte-identical guards** at Tasks 2, 4, 7 + final verification (no-loss ⇒ unchanged; additive column only).
- **Type consistency:** `neutral_losses: Vec<f64>` (Task 1) used in Task 3; `IonType{…, loss: bool}` + `is_loss()` (Task 2) used in Tasks 3/6/7; accession `Option<String>` (Task 1) used in Task 4.
- **Known risk / honesty:** Task 2 touches ~58 construction sites — mechanical and compiler-driven (`loss: false`), but the single largest churn; `{ .. }` match patterns are unaffected. Line numbers drift; every task gives surrounding-code anchors.
