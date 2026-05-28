# PXD001819 + TMT ID-rate Improvement Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
>
> **This is a bench-gated investigation, not a deterministic build.** The ultimate
> pass/fail gate for each change is the Percolator 1%-FDR PSM count on the VM, not a
> unit test. Phase 0 is a hard prerequisite: its measurements decide whether Phase 1
> has any yield. Phase 2's exact code edit is determined by a trace whose result is
> not known in advance — that task is written as an investigation protocol with a
> concrete fix location, candidate rules, and verification, NOT a pre-written diff.

**Goal:** Lift PXD001819 (14,755) and TMT (9,605) PSM counts @1% FDR toward +10% over current Rust, without regressing wall >3% or regressing Astral (36,715).

**Architecture:** One feature branch (`feat/id-rate-pxd001819-tmt`), one commit per change, each bench-gated on all three datasets via the existing VM harness, with in-place revert when the gate fails. Highest-leverage / lowest-risk first.

**Tech Stack:** Rust (workspace crates: `scoring`, `search`, `input`, `output`, `msgf-rust`); bench VM `pride-linux-vm` via SSH control socket `/tmp/msgfplus-bench.sock`; Percolator 3.7.1 in Docker; Python stdlib for the I5 trace harness.

---

## Conventions used throughout

**Build (local):** `cargo build --release -p msgf-rust` (the committed `.cargo/config.toml` sets `target-cpu=sandybridge`).

**Bit-identical regression gate (local):**
```bash
cargo test --release -p msgf-rust precursor_cal_off_pin_tsv_match_golden_after_sort
```
Expected: `test result: ok. 1 passed`. (Phase 2 changes top-1 selection, so this gate will legitimately change — see Task 2a Step 6 for how to regenerate goldens. Phases 0/1-diagnostic/3 must keep it green.)

**Workspace tests (local), CI skip list:**
```bash
cargo test --release --workspace -- \
  --skip charge_missing_spectrum_uses_per_charge_scored_spec \
  --skip spectrum_without_charge_tries_charge_range \
  --skip known_peptide_appears_in_top_n \
  --skip read_bsa_canno_text_format \
  --skip read_tryp_pig_bov_revcat_csarr_cnlcp \
  --skip tryp_pig_bov_revcat_full_set_loads \
  --skip match_spectra_output_invariant_across_thread_counts
```

**Clippy gate (local):** `cargo clippy --workspace --all-targets -- -D warnings` → clean.

**VM bench (per change):** ship source, build, run 3 datasets cal=auto, run Percolator. The canonical baseline numbers to beat/hold:

| Dataset | Rust @1% FDR (baseline) | Wall (baseline) |
|---|---:|---:|
| PXD001819 | 14,755 | ~0:54 |
| TMT | 9,605 | ~2:33 |
| Astral | 36,715 | ~6:28 |

**Bench ship gate (per change):** gains PXD or TMT @1% FDR; regresses none of PXD/TMT/Astral beyond ~±0.3% (≈ ±45 PXD, ±29 TMT, ±110 Astral); wall within ~3% on all three. Otherwise `git revert` the commit in place.

**VM socket precondition:** every VM step needs `ssh -S /tmp/msgfplus-bench.sock -O check root@pride-linux-vm` to succeed. If it fails, STOP and ask the human to run `ssh -M -S /tmp/msgfplus-bench.sock -fN root@pride-linux-vm`.

---

## Phase 0 — Diagnostic (prerequisite; no code change)

### Task 0: Measure what PXD001819 and TMT actually resolve to

**Files:** none (measurement only). Produces `docs/parity-analysis/notes/2026-05-28-id-rate-phase0-diagnostic.md` (gitignored notes dir; local record).

- [ ] **Step 1: Verify the VM socket is alive**

Run:
```bash
ssh -S /tmp/msgfplus-bench.sock -O check root@pride-linux-vm
```
Expected: `Master running (pid=NNNN)`. If not, STOP and request the human re-establish it.

- [ ] **Step 2: Capture each Rust run's resolver + calibration stderr**

The previous bench logs are at `/srv/data/msgf-bench/bench-v2024-results/{pxd001819,tmt,astral}-rust-auto.log`. Extract the resolver + calibration lines:
```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm \
  "for d in pxd001819 tmt astral; do echo \"== \$d ==\"; grep -iE 'Param resolver|instrument =|precursor.?cal|MassCalibrator|calibrat|confident|high.res|tolerance' /srv/data/msgf-bench/bench-v2024-results/\${d}-rust-auto.log; done"
```
Expected: lines like `Param resolver: auto-detected dominant activation method = HCD (instrument = ...)`. Record the resolved instrument string per dataset.

- [ ] **Step 3: Read the true instrument analyzer from each mzML**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm \
  "for f in /srv/data/msgf-bench/data/UPS1_5000amol_R1.mzML /srv/data/msgf-bench/tmt-data/a05058.mzML /srv/data/msgf-bench/astral-data/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML; do echo \"== \$f ==\"; grep -aoE 'MS:100(0079|0083|0084|0264|0484|0624|1906)|FTMS|ITMS|orbitrap|Q Exactive|Astral|ion trap' \$f | sort | uniq -c | head; done"
```
Map cvParams: `MS:1000484` orbitrap, `MS:1000079` FTICR, `MS:1000264` ion trap, `MS:1000083` radial ejection ion trap. This tells the *true* MS2 analyzer.

- [ ] **Step 4: Determine the tolerance branch actually taken**

For each dataset, the scoring/feature tolerance is high-res (20 ppm) iff the resolved `InstrumentType::is_high_resolution()` is true. Find which `InstrumentType` enum value `is_high_resolution()` returns true for:
```bash
grep -n -A8 "fn is_high_resolution" crates/model/src/*.rs crates/scoring/src/*.rs 2>/dev/null
```
Cross-reference the Step-2 resolved instrument against this. Record per dataset: `tolerance = 20ppm | 0.5Da`.

- [ ] **Step 5: Write the diagnostic table**

Create `docs/parity-analysis/notes/2026-05-28-id-rate-phase0-diagnostic.md` with the filled table:
```markdown
| Dataset | True MS2 analyzer | Resolved instrument | is_high_resolution | Feature/scoring tol | Calibration fired? |
|---|---|---|---|---|---|
| PXD001819 | ... | ... | ... | ... | ... |
| TMT | ... | ... | ... | ... | ... |
| Astral | (orbitrap/Astral) | (high-res) | true | 20ppm | yes (reference) |
```

- [ ] **Step 6: Decide Phase 1 scope (record in the same note)**

Decision rules:
- If PXD or TMT shows `is_high_resolution = false` / `tol = 0.5Da` but the true analyzer is orbitrap/FT → **Phase 1a is live** (iter20-style win available).
- If calibration shows "skipped" on PXD or TMT → **Phase 1b is live**.
- If both are already high-res + calibrated → skip Phase 1, go to Phase 2.

Write the decision explicitly. No commit (notes dir is gitignored); this task gates the rest.

---

## Phase 1 — Config levers (only the sub-tasks Phase 0 marked "live")

### Task 1a: Instrument-resolution / tolerance fix (CONDITIONAL on Phase 0 Step 6)

**Files:**
- Modify: `crates/msgf-rust/src/bin/msgf-rust.rs` (instrument resolution at ~L585–631) — exact edit depends on the Phase-0 root cause (see Step 2).
- Possibly: `crates/model/src/<instrument module>.rs` (`is_high_resolution` mapping) or `crates/search/src/match_engine.rs::compute_psm_features` (tolerance branch).

- [ ] **Step 1: Reproduce the mis-resolution locally**

If PXD001819's mzML is available locally, run; otherwise reason from the Phase-0 evidence. Confirm the code path: does `detect_instrument_type` return the orbitrap type, and does that type's `is_high_resolution()` return true? If detection returns `None` and the code falls back to low-res (`crates/msgf-rust/src/bin/msgf-rust.rs:605` comment "None ⇒ resolver picks LowRes"), the bug is in detection, not the tolerance branch.

- [ ] **Step 2: Make the targeted fix**

Two likely shapes (pick per Step 1 evidence):
- **Detection gap:** `detect_instrument_type` fails to recognize the analyzer cvParam present in the mzML. Add the missing cvParam mapping in `crates/input/src/mzml.rs` (the `CV_ANALYZER_*` constants near L20–40 and the match at ~L722–735).
- **Tolerance-branch gap:** detection is correct but `compute_psm_features` / scoring tolerance still uses 0.5 Da. Route the resolved high-res `InstrumentType` into the `is_high_resolution()` check.

Show the diff in the commit; do not guess it here.

- [ ] **Step 3: Local gates**

Run the bit-identical gate, workspace tests, clippy (commands in Conventions). The bit-identical golden is cal=off low-res, so a high-res-detection fix should NOT change it (different code path) — expected still green. If it changes, investigate before proceeding.

- [ ] **Step 4: Commit**
```bash
git add -A
git commit -m "fix(resolve): detect high-res instrument so 20ppm tolerance engages on PXD/TMT"
```

- [ ] **Step 5: VM bench + gate**

Ship, build, bench all 3 datasets (use the bench recipe in Task 2a Step 5). Apply the bench ship gate. If it regresses, `git revert HEAD` and record the result in the Phase-0 note.

### Task 1b: Calibration engagement (CONDITIONAL on Phase 0 Step 6)

**Files:**
- Modify: the MassCalibrator confident-PSM guard (find it):
```bash
grep -rn "confident\|min.*confident\|200\|calibrat" crates/scoring/src crates/search/src --include=*.rs | grep -i cal | head
```

- [ ] **Step 1: Locate the skip guard and log the actual count**

Find where calibration decides to skip (the "<200 confident PSMs" guard from memory). Confirm via the Phase-0 stderr what count PXD/TMT hit.

- [ ] **Step 2: Decide the change**

If the guard is skipping with, say, 150 confident PSMs, evaluate lowering the threshold OR using a wider confident-PSM definition. This is a parameter change to a single constant — show it in the commit.

- [ ] **Step 3: Local gates + commit**
```bash
cargo test --release -p msgf-rust precursor_cal_off_pin_tsv_match_golden_after_sort
git add -A && git commit -m "tune(cal): lower confident-PSM guard so calibration engages on PXD/TMT"
```

- [ ] **Step 4: VM bench + gate** (same recipe as Task 2a Step 5; revert if regresses).

---

## Phase 2 — Label-flip gap (hot-path; investigation + bench gate)

### Task 2a: H2 peak-rank assignment — match Java's `getPeakByMass` rule

**Files:**
- Modify: `crates/scoring/src/scoring/scored_spectrum.rs` — `nearest_peak_rank_in` (~L897–918) and/or the rank-assignment sort in `ScoredSpectrum::new` (~L199–237).
- Test: `crates/scoring/src/scoring/scored_spectrum.rs` (inline `#[cfg(test)]` module).

- [ ] **Step 1: Verify the VM socket; rebuild msgf-trace if needed**

```bash
ssh -S /tmp/msgfplus-bench.sock -O check root@pride-linux-vm
```
The I5 artifacts are committed at `docs/parity-analysis/notes/score-psm-trace-artifacts/`. The 5 traced scans: 41522, 34685, 23272, 23082, 16629.

- [ ] **Step 2: Identify one concrete RANK_DIFF ion**

From the committed artifacts, find the first RANK_DIFF ion on scan 41522:
```bash
cd docs/parity-analysis/notes/score-psm-trace-artifacts
python3 analyze.py 2>/dev/null | grep -i rank_diff | head
# Then inspect the specific ion in rust-trace-scan-41522.json:
python3 -c "import json; d=json.load(open('rust-trace-scan-41522.json')); [print(p['peptide'], [i for i in p['ions'] if 'rank' in i][:3]) for p in d]" 2>/dev/null | head
```
Record: `(theo_mz, rust_rank, java_rank)` for the first divergent ion. Decompress the Java side if needed: `gunzip -k java-trace-scan-41522.log.gz`.

- [ ] **Step 3: Read Java's `getPeakByMass` rule**

The Java reference is on the bench VM clone `/srv/data/msgf-bench/java-legacy-trace/`. Find the peak-selection rule:
```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm \
  "grep -rn -A20 'getPeakByMass' /srv/data/msgf-bench/java-legacy-trace/src/main/java/ | head -40"
```
Determine Java's tie-break when multiple peaks fall in tolerance: first-in-list? closest-m/z? highest-intensity? Rust currently uses **highest-intensity, strict `>` (first wins on tie)**. Record the exact divergence.

- [ ] **Step 4: Write a failing unit test that encodes Java's rule**

In `scored_spectrum.rs` tests, construct a tiny spectrum with two peaks inside one tolerance window where Java's rule and Rust's current rule disagree (use the real (theo_mz, ranks) from Step 2 if it reproduces, else a synthetic case matching the rule found in Step 3). Assert `nearest_peak_rank_in` returns Java's choice.
```rust
#[test]
fn nearest_peak_rank_matches_java_getpeakbymass_tiebreak() {
    // two peaks within tol of target; Java picks <rule from Step 3>
    // assert the rank returned equals Java's pick
}
```
Run: `cargo test --release -p scoring nearest_peak_rank_matches_java -- --nocapture` → Expected: FAIL.

- [ ] **Step 5: Implement the rule change; verify the unit test passes**

Change the selection in `nearest_peak_rank_in` (and the rank-assignment sort if the divergence is there) to match Java. Run the unit test → Expected: PASS. Then run the full `scoring` crate tests.

- [ ] **Step 6: Regenerate the bit-identical golden (top-1 legitimately changes)**

This change alters top-1 selection, so `precursor_cal_off.pin/.tsv` goldens shift. Regenerate them per the repo's golden-update procedure:
```bash
grep -rn "precursor_cal_off_pin_tsv_match_golden" crates/msgf-rust/tests/ | head
# Follow the test's documented regeneration path (typically an env var or a
# regen helper); inspect the test to confirm before regenerating.
```
Commit the regenerated goldens together with the code change so the gate stays meaningful.

- [ ] **Step 7: Commit**
```bash
git add -A
git commit -m "fix(scoring): match Java getPeakByMass tie-break in peak-rank assignment (H2)"
```

- [ ] **Step 8: Re-run the I5 trace harness on the VM; confirm RANK_DIFF drops**

Rebuild msgf-trace on the VM, re-run the 5-scan trace, re-run `analyze.py`. Expected: RANK_DIFF count drops from 301; LOGPROB_DIFF drops proportionally. If RANK_DIFF does NOT drop, the rule change was wrong — revisit Step 3.

- [ ] **Step 9: VM bench + gate (the decisive test)**

Bench recipe:
```bash
# ship source (rsync/scp changed crates), then on VM:
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm \
  "cd /srv/data/msgf-bench/iter5-build && PATH=\$HOME/.cargo/bin:\$PATH cargo build --release -p msgf-rust"
# run 3 datasets cal=auto (reuse run_bench_*.sh pattern), then percolator:
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm \
  "for l in pxd001819-rust-auto astral-rust-auto tmt-rust-auto; do bash /srv/data/msgf-bench/run_percolator_docker.sh /srv/data/msgf-bench/bench-<TAG>-results/\$l.pin /srv/data/msgf-bench/bench-<TAG>-percolator \$l; done"
```
Apply the bench ship gate. **If Percolator regresses on any dataset beyond noise, `git revert HEAD` (and the goldens commit) — this is the n=9+ audit risk materializing.** Record the result.

### Task 2b: H3 log-prob indexing (CONDITIONAL — only if 2a passed its gate)

**Files:**
- Modify: `crates/scoring/src/scoring/rank_scorer.rs` (table indexing / clamping).
- Test: inline `#[cfg(test)]` in `rank_scorer.rs`.

- [ ] **Step 1: Re-run analyze.py; isolate a pure-H3 ion**

After 2a, find an ion with LOGPROB_DIFF but NO RANK_DIFF (same rank, different value) — there were 307 such cases pre-fix. Record `(partition, ion_type, rank, rust_logprob, java_logprob)`.

- [ ] **Step 2: Compare Rust's table lookup to Java's for that ion**

Trace `rank_scorer.rs` indexing (`idx = min(rank, max_rank).max(1) - 1`; missing uses `max_rank` slot) against Java's `getNodeScore` table lookup for the same (partition, ion, rank). Identify the off-by-one / clamping / table-content difference.

- [ ] **Step 3: Failing unit test encoding Java's lookup**

Construct a `RankScorer` from a known partition table and assert the log-prob at a specific rank matches Java's value for the recorded ion. Run → Expected FAIL.

- [ ] **Step 4: Fix indexing; unit test passes; regenerate goldens; commit**
```bash
cargo test --release -p scoring <test_name>
git add -A && git commit -m "fix(scoring): align per-rank log-prob table indexing with Java (H3)"
```

- [ ] **Step 5: VM bench + gate** (same recipe as 2a Step 9; revert if regresses).

---

## Phase 3 — Additive PIN features (safety net; only if Phases 1–2 < +10%)

### Task 3a: Add `MeanMatchedRank` PIN column

**Files:**
- Modify: `crates/search/src/psm.rs` (`PsmFeatures` — add field).
- Modify: `crates/search/src/match_engine.rs::compute_psm_features` (compute it).
- Modify: `crates/output/src/pin.rs` (`write_header` + `write_psm_row` — emit it, inserted before `Peptide` next to `EdgeScore`).
- Test: `crates/output/src/pin.rs` tests + `crates/output/tests/output_pin_schema_parity.rs`.

- [ ] **Step 1: Failing test for the new column in the header**

In `crates/output/tests/output_pin_schema_parity.rs`, add an assertion that the header contains `MeanMatchedRank` at the expected position. Run → Expected FAIL.

- [ ] **Step 2: Add the field to `PsmFeatures`**
```rust
// crates/search/src/psm.rs, in PsmFeatures
pub mean_matched_rank: f32,
```

- [ ] **Step 3: Compute it in `compute_psm_features`**

After the matched-ion loop, average `nearest_peak_rank` over matched b/y ions (0.0 if none matched):
```rust
// crates/search/src/match_engine.rs, in compute_psm_features
let mean_matched_rank = if matched_rank_count > 0 {
    matched_rank_sum as f32 / matched_rank_count as f32
} else { 0.0 };
// ...set features.mean_matched_rank = mean_matched_rank;
```
(Accumulate `matched_rank_sum`/`matched_rank_count` inside the existing matched-ion loop where `nearest_peak_rank` is already called.)

- [ ] **Step 4: Emit the column**

In `crates/output/src/pin.rs::write_header`, insert `MeanMatchedRank` immediately before `EdgeScore`; in `write_psm_row`, write `features.mean_matched_rank` in the same position. Keep `EdgeScore` and `Peptide`/`Proteins` last.

- [ ] **Step 5: Update schema-parity test + run**
```bash
cargo test --release -p output
cargo test --release --workspace -- <CI skip list>
```
Expected: PASS (column count +1 everywhere it's asserted).

- [ ] **Step 6: Commit**
```bash
git add -A && git commit -m "feat(pin): add additive MeanMatchedRank feature column"
```

- [ ] **Step 7: VM bench + gate**

Additive columns never change top-1, so PXD/TMT/Astral PSM *search* output is byte-identical pre-Percolator except the new column; the only delta is whether Percolator extracts signal from it. Bench all 3; keep if it gains, revert if it regresses (it should at worst be flat).

### Task 3b: Add `ScoreFractionTop1Split` PIN column (only if 3a flat and still < target)

**Files:**
- Modify: `crates/scoring/src/scoring/psm_score.rs::score_psm` — return the max single-split contribution alongside the sum (change return type to `(f32, f32)` or add a sibling fn to avoid disturbing callers).
- Modify: `crates/search/src/psm.rs`, `match_engine.rs`, `crates/output/src/pin.rs` as in 3a.

- [ ] **Step 1: Failing header test** (as 3a Step 1, for `ScoreFractionTop1Split`).
- [ ] **Step 2: Track max split in `score_psm`**

Inside the split loop that accumulates `total += contribution`, also track `max_contrib = max_contrib.max(contribution)`. Expose `max_contrib / total` (guard total==0 → 0.0). To avoid perturbing the hot `score_psm` signature used in ranking, add a separate `score_psm_with_split_max(...) -> (f32, f32)` used only in `compute_psm_features`, OR compute the fraction in `compute_psm_features` by re-deriving splits (prefer the former to avoid double work).

- [ ] **Step 3–6:** wire field → compute → emit → schema test → bench gate, exactly as 3a Steps 2–7.

---

## Final: cumulative bench + branch close-out

### Task 4: Cumulative bench, PR, ship/revert reconciliation

- [ ] **Step 1: Full 3-dataset bench at branch HEAD**

Confirm the cumulative PXD/TMT/Astral @1% FDR and walls vs baseline. Build the result table.

- [ ] **Step 2: Reconcile against the stretch target**

Record cumulative deltas. If +10% reached, note it. If not, the shipped subset is whatever net-gained under the hard gate (the +10% is a direction, not a revert-all gate per the spec).

- [ ] **Step 3: Open the PR**
```bash
git push -u origin feat/id-rate-pxd001819-tmt
gh pr create --base dev --title "perf(id-rate): close PXD001819 + TMT label-flip gap vs Java" --body "<results table + per-phase outcomes + reverts>"
```

- [ ] **Step 4: Watch CI to green** (Lint + 3 OS tests + CodeRabbit), fix any failures.

---

## Self-review notes

- **Spec coverage:** Phase 0 → Task 0; Phase 1a → Task 1a; Phase 1b → Task 1b; Phase 2 H2 → Task 2a; Phase 2 H3 → Task 2b; Phase 3 → Tasks 3a/3b; success-criteria reconciliation → Task 4. All spec sections covered.
- **Conditionality is intentional:** Tasks 1a/1b/2b/3a/3b are explicitly gated on prior-task outcomes — this matches the spec's "outcome of each phase decides the next." An executor must honor the CONDITIONAL markers, not run every task blindly.
- **Phase 2 code is deliberately not pre-written:** the exact tie-break/indexing edit is unknowable until Steps 2–3 of each task read the real trace + Java source. The tasks pin the file:line, the candidate rules, the unit-test shape, and the bench gate — which is the maximum honest specificity here.
- **Golden regeneration** is called out explicitly in 2a/2b because those tasks change top-1 selection and will move the committed PIN/TSV goldens.
