# R-1 Retention-Layer Empirical Test — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` (recommended) or `superpowers:executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Empirically test whether `TopNQueue::push`'s strict-greater eviction is the dominant cause of the 14K-PSM raw-target gap on the Astral no-mods bench (Java 89,479 vs Rust 75,457). One file changed, one function modified, one unit test added, one bench cycle.

**Architecture:** Replace `Ordering::Greater` check in `TopNQueue::push` with a `match` that keeps tied PSMs at capacity (Java `DBScanner.java:540`-equivalent). Queue may exceed `capacity` when ties exist — this is intentional and matches Java's unbounded-tied-PSM behavior. Validate via a new unit test, then bench Astral on `pride-linux-vm` to measure raw target/decoy delta vs Java baseline.

**Tech Stack:** Rust (cargo), Bash (bench harness via SSH ControlMaster), awk (PSM counting from `.pin` Label column).

**Spec:** [`docs/parity-analysis/specs/2026-05-18-r1-tie-retention-test-design.md`](../specs/2026-05-18-r1-tie-retention-test-design.md)

**Branch:** `rust-implement`, starting from HEAD `01e7062` (spec commit on top of `b1d45bb`).

---

## File Structure

**Modified:**
- `rust/crates/search/src/psm.rs` — modify `TopNQueue::push` (around line 163-173), add one unit test in `#[cfg(test)] mod tests` (after line 332)

**No new files created.**

**Bench-side files** (on pride-linux-vm, not in repo):
- `/srv/data/msgf-bench/msgf-rust-iter6.tgz` — source tarball
- `/srv/data/msgf-bench/track-iter6-build/` — extracted source + target build
- `/srv/data/msgf-bench/bench-iter6-results/astral-rust-r1.{pin,log}` — bench output

**Doc updates after bench (Task 5):**
- `docs/parity-analysis/notes/2026-05-18-r1-bench-results.md` — new file documenting the empirical result

---

## Task 1: Write the failing unit test for tie retention

**Files:**
- Modify: `rust/crates/search/src/psm.rs` (add test in `#[cfg(test)] mod tests` block, after line 332)

The existing `make_match(spectrum_idx, score)` helper builds a `PsmMatch` with default `spec_e_value = 1.0` (sentinel). Three calls with the same score produce three identical PSMs in `Ord` terms — perfect for the tie test.

- [ ] **Step 1: Add the failing test**

Append to `rust/crates/search/src/psm.rs` immediately after the existing `lower_score_dropped_when_full` test (around line 323), inside the `#[cfg(test)] mod tests` block:

```rust
    #[test]
    fn topn_queue_keeps_ties_at_capacity() {
        // R-1 fix: Java's DBScanner keeps tied PSMs at capacity
        // (DBScanner.java:540 raw-score retention; DBScanner.java:745 SpecE
        // merge). Rust's TopNQueue must mirror this — strict-greater eviction
        // was dropping ties Java keeps, plausibly causing the Astral 14K raw-
        // target gap. See
        // docs/parity-analysis/notes/2026-05-18-piecewise-fixes-dont-work.md
        // (Open: retention layer, R-1).
        let mut q = TopNQueue::new(1);
        q.push(make_match(0, 100.0));
        q.push(make_match(0, 100.0));
        q.push(make_match(0, 100.0));
        assert_eq!(
            q.len(),
            3,
            "all three tied PSMs should be retained at capacity=1 (Java parity, R-1)"
        );
    }
```

- [ ] **Step 2: Run the test, verify it FAILS**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust
cargo test -p search --lib psm::tests::topn_queue_keeps_ties_at_capacity -- --nocapture
```

Expected: FAIL with assertion `left: 1, right: 3` (current Rust drops two of the three tied PSMs).

- [ ] **Step 3: Commit the failing test**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed
git add rust/crates/search/src/psm.rs
git commit -m "test(search): failing test exposes TopNQueue ties dropped at capacity (R-1)"
```

---

## Task 2: Implement the R-1 fix

**Files:**
- Modify: `rust/crates/search/src/psm.rs` lines 163-173 (the `push` function body)

- [ ] **Step 1: Replace `push` body to keep ties at capacity**

Find the existing function at `rust/crates/search/src/psm.rs:163-173`:

```rust
    pub fn push(&mut self, m: PsmMatch) {
        if self.heap.len() < self.capacity as usize {
            self.heap.push(Reverse(m));
        } else if let Some(Reverse(top)) = self.heap.peek() {
            // `m > top` in natural ordering means m is better.
            if m.cmp(top) == std::cmp::Ordering::Greater {
                self.heap.pop();
                self.heap.push(Reverse(m));
            }
        }
    }
```

Replace with:

```rust
    pub fn push(&mut self, m: PsmMatch) {
        if self.heap.len() < self.capacity as usize {
            self.heap.push(Reverse(m));
        } else if let Some(Reverse(top)) = self.heap.peek() {
            match m.cmp(top) {
                std::cmp::Ordering::Greater => {
                    // m is strictly better than the worst retained PSM: evict
                    // the worst, insert m.
                    self.heap.pop();
                    self.heap.push(Reverse(m));
                }
                std::cmp::Ordering::Equal => {
                    // R-1 (2026-05-18): Java's DBScanner.java:540 keeps tied
                    // PSMs at capacity (and DBScanner.java:745 keeps SpecE
                    // ties on the per-spectrum merge). Rust now matches.
                    // The queue may exceed `capacity` when ties exist —
                    // `capacity` becomes a *minimum* top-N, not a hard cap.
                    // Spec:
                    // docs/parity-analysis/specs/2026-05-18-r1-tie-retention-test-design.md
                    self.heap.push(Reverse(m));
                }
                std::cmp::Ordering::Less => {
                    // m is strictly worse than the worst retained PSM: drop.
                }
            }
        }
    }
```

Also update the docstring directly above the function (around line 155-162). The current docstring says:

```rust
    /// Insert a PSM. The queue keeps the `capacity` *best* PSMs.
    ///
    /// "Best" = smallest `spec_e_value` first (then largest `score` for ties).
    /// The min-heap (via `Reverse<PsmMatch>`) puts the *worst* PSM at the top
    /// so it can be evicted when over capacity.
    ///
    /// Before `compute_spec_e_values_for_spectrum` runs, all PSMs have
    /// `spec_e_value = 1.0` and the secondary `score` key governs eviction.
```

Replace with:

```rust
    /// Insert a PSM. The queue keeps **at least** `capacity` of the *best*
    /// PSMs, plus any additional PSMs tied with the current worst.
    ///
    /// "Best" = smallest `spec_e_value` first (then largest `score` for ties).
    /// The min-heap (via `Reverse<PsmMatch>`) puts the *worst* PSM at the top
    /// so it can be evicted when a strictly-better PSM arrives.
    ///
    /// Before `compute_spec_e_values_for_spectrum` runs, all PSMs have
    /// `spec_e_value = 1.0` and the secondary `score` key governs eviction.
    ///
    /// **Tie handling (R-1, 2026-05-18):** when the queue is at capacity and
    /// a new PSM is `Equal` (in `Ord` terms) to the worst retained PSM, the
    /// new PSM is inserted WITHOUT evicting the tied one. This matches
    /// Java's `DBScanner.java:540` (`size < n OR score == worst → add`).
    /// As a result, the queue can grow beyond `capacity` when ties exist;
    /// `capacity` becomes a *minimum* top-N, not a hard cap.
```

- [ ] **Step 2: Run the unit test, verify it PASSES**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust
cargo test -p search --lib psm::tests::topn_queue_keeps_ties_at_capacity -- --nocapture
```

Expected: PASS with `test result: ok. 1 passed`.

- [ ] **Step 3: Run the full `psm` test module, verify no regressions in adjacent tests**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust
cargo test -p search --lib psm:: -- --nocapture 2>&1 | tail -20
```

Expected: ALL tests pass (`empty_queue`, `queue_below_capacity_keeps_everything`, `queue_at_capacity_keeps_top_n_by_score`, `lower_score_dropped_when_full`, `psm_match_clones_correctly`, `topn_queue_keeps_ties_at_capacity`, and any spec-e-value ordering tests). Several of these tests use `make_match(0, s)` with all-tied `spec_e_value=1.0` — the new tie-keeping behavior may cause them to retain *more* PSMs than before. Check each that asserts a queue size or sorted-vec length carefully.

Specifically, **`queue_at_capacity_keeps_top_n_by_score` is the test most likely affected**: it pushes `[1.0, 5.0, 2.0, 4.0, 3.0]` into capacity=3 and expects `[5.0, 4.0, 3.0]`. All values are unique (no ties), so the test should still pass. Run it explicitly:

```bash
cargo test -p search --lib psm::tests::queue_at_capacity_keeps_top_n_by_score -- --nocapture
```

Expected: PASS.

- [ ] **Step 4: Run the broader scoring + search test suites, verify no regressions**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust
cargo test -p scoring -p search 2>&1 | grep -E "^test result|FAIL" | head -30
```

Expected: All `ok` except the **pre-existing** 3 failures in `match_engine_smoke` (`charge_missing_spectrum_uses_per_charge_scored_spec`, `known_peptide_appears_in_top_n`, `spectrum_without_charge_tries_charge_range`) that fail at baseline `b1d45bb` too (verified in earlier sessions). No NEW failures.

If any other test fails (e.g., `gf_java_parity`, `score_psm_pxd001819_parity`, `score_psm_vs_gf_dp_edge_parity`, any other `match_engine_*` test, any `output_pin_*` test): investigate before continuing. The R-1 change should NOT affect anything except queue retention; if a downstream consumer depends on `queue.len() <= capacity` as an invariant, it will break here.

- [ ] **Step 5: Commit the R-1 fix**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed
git add rust/crates/search/src/psm.rs
git commit -m "fix(search): TopNQueue keeps tied PSMs at capacity (R-1)

Java's DBScanner.java:540 keeps tied PSMs at capacity during raw-score
retention; DBScanner.java:745 does the same at the per-spectrum SpecE
merge. Rust's TopNQueue::push used strict-greater eviction, dropping
ties Java keeps. On the Astral no-mods bench this plausibly caused the
14K-PSM raw-target gap (Java 89,479 / Rust 75,457).

The queue may now grow past 'capacity' when ties exist — 'capacity'
becomes a minimum top-N, not a hard cap. Java has unbounded ties on
score (or SpecE) at the worst-retained slot; Rust now matches.

This is an empirical hypothesis test, not a parity push. Spec at
docs/parity-analysis/specs/2026-05-18-r1-tie-retention-test-design.md.
Production-bench validation follows."
```

---

## Task 3: Sync to pride-linux-vm and build release binary

**Files:**
- Create on VM: `/srv/data/msgf-bench/msgf-rust-iter6.tgz` (source tarball)
- Create on VM: `/srv/data/msgf-bench/track-iter6-build/` (extracted source + cargo target)

Prerequisites:
- SSH ControlMaster socket alive at `/tmp/msgfplus-bench.sock`. If not, the user must re-establish via `ssh -M -S /tmp/msgfplus-bench.sock -fN root@pride-linux-vm.ebi.ac.uk` (clear stale host key first with `ssh-keygen -R pride-linux-vm.ebi.ac.uk` if needed).
- Bundled `.param` files already on VM at `/srv/data/msgf-bench/track-iter3-build/src/main/resources/ionstat/` (left over from earlier iterations; reused via `cp -r`).

- [ ] **Step 1: Tar the local source (excluding target/)**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed
tar --exclude='rust/target' --exclude='rust/**/target' -czf /tmp/msgf-rust-iter6.tgz rust
ls -lh /tmp/msgf-rust-iter6.tgz
```

Expected: tarball ~260K, output line like `-rw-r--r--  1 yperez  staff  265K ... /tmp/msgf-rust-iter6.tgz`.

- [ ] **Step 2: Upload tarball + extract on VM**

```bash
scp -o ControlPath=/tmp/msgfplus-bench.sock /tmp/msgf-rust-iter6.tgz \
    root@pride-linux-vm.ebi.ac.uk:/srv/data/msgf-bench/

ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm.ebi.ac.uk '
    rm -rf /srv/data/msgf-bench/track-iter6-build &&
    mkdir -p /srv/data/msgf-bench/track-iter6-build &&
    cd /srv/data/msgf-bench/track-iter6-build &&
    tar -xzf /srv/data/msgf-bench/msgf-rust-iter6.tgz 2>&1 | grep -v "LIBARCHIVE" | head -3 &&
    cp -r /srv/data/msgf-bench/track-iter3-build/src ./ &&
    cd rust &&
    printf "[toolchain]\nchannel = \"stable\"\ncomponents = [\"rustfmt\", \"clippy\"]\n" > rust-toolchain.toml &&
    grep -q "Equal => {" crates/search/src/psm.rs && echo "R-1 fix present" || echo "MISSING R-1 fix"
'
```

Expected output: `R-1 fix present` (sentinel check that the `Equal` branch made it through the tarball).

- [ ] **Step 3: Build release binary in background**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm.ebi.ac.uk '
    cd /srv/data/msgf-bench/track-iter6-build/rust &&
    nohup bash -c "cargo build --release -p msgf-rust > /tmp/cargo-build-iter6.log 2>&1" >/dev/null 2>&1 &
    echo "build PID: $!"
'
```

Expected: prints `build PID: <pid>`. Build kicks off in background.

- [ ] **Step 4: Wait for the build, verify binary exists**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm.ebi.ac.uk '
    until ! pgrep -x cargo >/dev/null && ! pgrep -x rustc >/dev/null; do
        sleep 30
    done
    echo "build done"
    tail -3 /tmp/cargo-build-iter6.log
    ls -la /srv/data/msgf-bench/track-iter6-build/rust/target/release/msgf-rust
'
```

Expected output: `build done`, last build log lines like `Finished 'release' profile [optimized] target(s) in N s`, and binary listing showing the executable (~2 MB).

If build FAILS: check `/tmp/cargo-build-iter6.log` on the VM. The most common failure is the toolchain — the prior `rust-toolchain.toml` was overridden to `stable` (Rust 1.95.0 on the VM); if a transitive dep needs an even newer Rust, this step will surface it. Do NOT proceed to Task 4 until the binary exists.

---

## Task 4: Run the Astral no-mods bench + count results

**Files:**
- Create on VM: `/srv/data/msgf-bench/bench-iter6-results/astral-rust-r1.{pin,log}`

- [ ] **Step 1: Launch the bench in nohup background**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm.ebi.ac.uk '
    RUST=/srv/data/msgf-bench/track-iter6-build/rust/target/release/msgf-rust
    OUT=/srv/data/msgf-bench/bench-iter6-results
    mkdir -p "$OUT"
    nohup bash -c "/usr/bin/time -v $RUST \
        --spectrum /srv/data/msgf-bench/astral-data/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML \
        --database /srv/data/msgf-bench/astral-data/ProteoBenchFASTA_MixedSpecies_HYE.fasta \
        --output-pin $OUT/astral-rust-r1.pin \
        --precursor-tol-ppm 10 --isotope-error-min=-1 --isotope-error-max=2 \
        --ntt 2 --max-missed-cleavages 2 --min-peaks 10 \
        --min-length 6 --max-length 40 --charge-min 2 --charge-max 4 \
        --top-n 1 --threads 8 > $OUT/astral-rust-r1.log 2>&1" >/dev/null 2>&1 &
    echo "iter6 astral PID: $!"
'
```

Expected: prints `iter6 astral PID: <pid>`. **Important:** use `nohup` so the bench survives SSH disconnections (6 dropped in earlier sessions).

- [ ] **Step 2: Poll for bench completion (15-30 min wall, allow up to 26 min per spec)**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm.ebi.ac.uk '
    until ! pgrep -f "astral-rust-r1.pin" >/dev/null; do
        sleep 60
    done
    echo "bench done"
    grep -E "Elapsed|Maximum resident|Exit status" /srv/data/msgf-bench/bench-iter6-results/astral-rust-r1.log | head -3
    ls -la /srv/data/msgf-bench/bench-iter6-results/astral-rust-r1.pin
'
```

Expected: `bench done`, wall time line like `Elapsed (wall clock) time (h:mm:ss): 8:36.69`, max RSS ~9.7 GB, exit status 0, pin file ~36 MB.

If wall exceeds 26 min: queue is growing pathologically. Kill the process and follow the rollback plan (Task 5 / failure branch).

If exit status is non-zero: check the log for the error. R-1 should not cause runtime errors, but if `queue.len()` is now used as an indexed bound somewhere in the codebase that assumed `<= capacity`, this is where it'd surface.

- [ ] **Step 3: Count raw target/decoy PSMs from the .pin**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm.ebi.ac.uk '
    PIN=/srv/data/msgf-bench/bench-iter6-results/astral-rust-r1.pin
    T=$(awk -F"\t" "NR>1 && \$2==1 {c++} END {print c+0}" "$PIN")
    D=$(awk -F"\t" "NR>1 && \$2==-1 {c++} END {print c+0}" "$PIN")
    RATIO=$(awk -v t=$T -v d=$D "BEGIN {if (d>0) printf \"%.3f\", t/d; else print \"inf\"}")
    echo "iter6 astral-r1 targets=$T decoys=$D T/D=$RATIO"
    echo "vs Java baseline: targets=89479 decoys=46792 T/D=1.912"
    echo "vs b1d45bb (pre-R-1): targets=75457 decoys=46208 T/D=1.633"
    echo "success bar: targets >= 86674, T/D >= 1.47, wall <= 26min"
'
```

Expected output: prints the iter6 counts + baselines for direct comparison.

- [ ] **Step 4: Capture the iter6 result snapshot for the report**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm.ebi.ac.uk '
    PIN=/srv/data/msgf-bench/bench-iter6-results/astral-rust-r1.pin
    LOG=/srv/data/msgf-bench/bench-iter6-results/astral-rust-r1.log
    T=$(awk -F"\t" "NR>1 && \$2==1 {c++} END {print c+0}" "$PIN")
    D=$(awk -F"\t" "NR>1 && \$2==-1 {c++} END {print c+0}" "$PIN")
    WALL=$(grep "Elapsed (wall clock)" "$LOG" | awk -F": " "{print \$NF}")
    RSS=$(grep "Maximum resident" "$LOG" | awk "{print \$NF}")
    EXIT=$(grep "Exit status" "$LOG" | awk "{print \$NF}")
    cat <<EOF
iter6_astral_r1_results:
  targets: $T
  decoys: $D
  wall_time: $WALL
  max_rss_kb: $RSS
  exit_status: $EXIT
EOF
'
```

Copy the output into the local clipboard / scratch — Task 5 reads it.

---

## Task 5: Analyze, document, decide (commit / commit-and-plan-next / revert)

**Files:**
- Create: `docs/parity-analysis/notes/2026-05-18-r1-bench-results.md`

- [ ] **Step 1: Create the results note**

Use the actual values from Task 4 Step 4. Template:

```markdown
# R-1 retention-layer test — empirical results

_2026-05-18. Branch `rust-implement` at HEAD (iter6 = b1d45bb + R-1)._

## Spec

[`docs/parity-analysis/specs/2026-05-18-r1-tie-retention-test-design.md`](../specs/2026-05-18-r1-tie-retention-test-design.md)

## Bench config

Astral no-mods, threads=8, `pride-linux-vm`. Matches Java's reference bench config exactly
(no `--mod`, `--precursor-tol-ppm 10 --isotope-error-min=-1 --isotope-error-max=2`).

## Results

| Metric | Java baseline | b1d45bb (pre-R-1) | iter6 (with R-1) | vs Java | vs b1d45bb |
|---|---:|---:|---:|---:|---:|
| Raw targets | 89,479 | 75,457 | **<INSERT>** | <INSERT> | <INSERT> |
| Raw decoys | 46,792 | 46,208 | **<INSERT>** | <INSERT> | <INSERT> |
| T/D ratio | 1.912 | 1.633 | **<INSERT>** | — | — |
| Wall time | — | ~8:36 | **<INSERT>** | — | — |
| Max RSS | — | ~9.7 GB | **<INSERT>** | — | — |

## Decision per spec outcome table

<INSERT: one of:>

- ≥86,674 targets, T/D ≥1.47: **Hypothesis validated.** R-1 is the dominant cause of the
  14K-PSM gap. Plan R-2 (per-charge GF) for the next iteration. The R-1 commit stays.
- 78,000-86,673 targets: R-1 is real but partial. R-1 + R-2 coupled; next iteration runs
  both together. The R-1 commit stays.
- 75,458-77,999 targets: R-1 effect is small. Commit stays IFF T/D ratio + tests OK; pivot
  investigation to other layers (R-3, R-4, F-1, or audit-tier scoring fixes with retention
  now correct).
- ≤75,457 targets OR T/D < 1.47: **Hypothesis falsified.** Revert R-1. Document why ties
  weren't the gap. Pivot.

## Next iteration

<INSERT: based on the decision above>
```

Fill in `<INSERT>` with the actual values from Task 4 Step 4. Compute the deltas:
- `vs Java`: `(iter6 - 89479)` (negative if iter6 < Java)
- `vs b1d45bb`: `(iter6 - 75457)` (positive if R-1 closed the gap)

- [ ] **Step 2: Commit the results note**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed
git add docs/parity-analysis/notes/2026-05-18-r1-bench-results.md
git commit -m "docs(parity): R-1 retention-layer test empirical results

iter6 = b1d45bb + R-1 fix (TopNQueue keeps tied PSMs at capacity).

Astral no-mods bench (threads=8):
  targets: <INSERT> (Java 89,479 / b1d45bb 75,457)
  decoys:  <INSERT> (Java 46,792 / b1d45bb 46,208)
  T/D:     <INSERT> (Java 1.912 / b1d45bb 1.633)
  wall:    <INSERT> (b1d45bb ~8:36, spec cap 26min)

Decision: <INSERT one-line outcome verdict>
Next:     <INSERT plan or revert>"
```

(Replace `<INSERT>` with the actual values.)

- [ ] **Step 3: If hypothesis FALSIFIED, revert R-1**

Only do this step if Task 5 Step 1's decision branch is "Revert R-1":

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed
git revert HEAD~1   # reverts the R-1 fix commit (Task 2's commit)
# Note: the test commit (Task 1) stays — it documents the divergence even
# if the fix is wrong. The test now passes a posteriori vs the original
# Rust behavior with capacity=1, so it would FAIL again after revert.
# Mark it #[ignore] before pushing the revert:
```

Edit `rust/crates/search/src/psm.rs` to add `#[ignore = "R-1 hypothesis falsified 2026-05-18; see docs/parity-analysis/notes/2026-05-18-r1-bench-results.md"]` above `#[test] fn topn_queue_keeps_ties_at_capacity`. Then:

```bash
git add rust/crates/search/src/psm.rs
git commit -m "test(search): #[ignore] R-1 tie test pending hypothesis re-examination"
```

- [ ] **Step 4: If hypothesis VALIDATED or PARTIAL, leave R-1 committed**

The Task 2 commit stays. No further action this iteration.

- [ ] **Step 5: Update memory entries with the empirical result**

If the result is significant either way (validated, partial, or falsified — all are interesting findings), update:

- `/Users/yperez/.claude/projects/-Users-yperez-work-msgfplus-workspace/memory/feedback_piecewise_alignment_doesnt_work.md` — add a paragraph in the "Why" section recording the iter6 result and what it means for the retention-layer hypothesis.

- `/Users/yperez/.claude/projects/-Users-yperez-work-msgfplus-workspace/memory/MEMORY.md` — update the line for `piecewise-alignment-doesnt-work` if the conclusion changes substantially.

For a clean validation (≥86,674 targets): the memory should say "R-1 confirmed: tie retention was the dominant cause of the 14K-target gap; fix-order is now R-2 → R-3 → R-4 → F-1 → audit-tier features."

For a falsification: the memory should say "R-1 falsified: tie retention is not the dominant cause; the 14K-target gap has another source. Audit + R-1 together do not explain Astral. Need a fundamentally different investigation strategy."

For a partial result: "R-1 closes X PSMs but R-2 is needed for the rest; couple them in the next iteration."

---

## Self-Review (writing-plans skill checklist)

**1. Spec coverage:**

| Spec section | Plan task |
|---|---|
| Problem | covered in plan intro |
| Goal | covered in plan intro |
| Success criteria (≥86,674 targets, T/D ≥1.47, wall ≤26min, tests pass) | Task 4 Step 3 + Task 5 Step 1 outcome table |
| Architecture (`Ordering::Equal` branch in push) | Task 2 Step 1 (exact diff) |
| Testing (one unit test) | Task 1 (test code) + Task 2 Step 2 (verify pass) |
| Bench protocol (nohup, no-mods, threads=8) | Task 3 + Task 4 (commands match spec verbatim) |
| Rollback plan (revert + ignore test) | Task 5 Step 3 |
| Out of scope (R-2, R-3, R-4, F-1) | not in plan, explicitly noted in plan intro |
| Risks (queue grows unbounded, downstream perf) | Task 4 Step 2 (wall-time check) |
| Estimated effort (1-1.5h) | Tasks 1+2 ~45min, Task 3 ~20min, Task 4 ~25min, Task 5 ~15min — totals ~1h45m which is in the spec's 1-1.5h ballpark (slightly over but within tolerance) |

All spec requirements have a corresponding task. No gaps.

**2. Placeholder scan:**

- `<INSERT>` in Task 5 Step 1 + Step 2: these are FILL-IN-FROM-BENCH-RESULTS markers, not unimplemented work — explicitly marked with `<INSERT: ...>` so an agentic worker fills them with the empirical values from Task 4. This is appropriate because the plan can't know the answer in advance.
- No "TBD", no "TODO", no "implement later", no "add appropriate error handling", no "write tests for the above", no "similar to Task N".

**3. Type consistency:**

- `PsmMatch` (struct) — consistent.
- `TopNQueue::push(&mut self, m: PsmMatch)` — consistent signature throughout.
- `make_match(spectrum_idx, score)` — the existing test helper (line 259 of `psm.rs`) — used consistently in Task 1.
- `std::cmp::Ordering::Greater | Equal | Less` — consistent enum names.

No type / name mismatches.

---

## Execution Handoff

Plan complete and saved to `docs/parity-analysis/plans/2026-05-18-r1-tie-retention-test-plan.md`.
