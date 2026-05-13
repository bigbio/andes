# `score_psm` under-scoring fix — design

**Date:** 2026-05-13
**Author:** brainstorm session captured to spec
**Status:** approved design; ready for implementation plan
**Prerequisite reading:** `docs/parity-analysis/reports/2026-05-13-score-psm-undercount-finding.md` (the forensic finding)

## Goal

Fix the `score_psm` under-scoring bug that causes Rust to compute RawScore
values roughly 1/3 of Java's for identical (peptide, scan, charge) inputs on
PXD001819 and Astral. The bug yields **−22%** and **−31%** Percolator @ 1%
FDR vs Java respectively. The TMT win recovered earlier in the iteration
(this session) must not regress.

Acceptance criteria for "fixed":

1. A new `score_psm_pxd001819_scan_28787_matches_java` unit test asserts
   `RawScore == 297` (Java's value) for `IVNEEFDQLEEDTPVYK` at charge 2 on
   PXD001819 scan 28787. Plus 2–3 sister unit tests on other worst-gap scans
   (28825 / 33606 / 32395) for bug-class detection.
2. All existing 318 workspace lib tests still pass.
3. `gf_java_parity` (release, 1.0 OOM tolerance) still passes.
4. **VM bench Percolator @ 1% FDR:** PXD001819 ≥ 14,800, Astral ≥ 33,000,
   TMT ≥ 10,500 (no regression).
5. Any perf wins introduced by the suspected commits (most likely `0af1a37`
   FastScorer prefix/suffix cache) are PRESERVED — fix must be surgical, not
   a full revert.

## Architecture

Six sequential phases on a dedicated worktree `astral-speed-score-fix`
on a new branch `fix/score-psm-undercount` based on `rust-implement`.
Phase boundaries are firm: each phase's output is the next phase's input.
If a phase blocks (e.g., bisect identifies a commit whose diff doesn't
obviously explain a 3× score drop), the process ESCALATES rather than
guessing forward.

The investigation strategy is **bisect-driven** (find the exact bad commit
via git bisect with a deterministic RawScore oracle), the fix is
**surgical** (single commit on top of the bad-commit's intent, preserving
any perf), and validation is **two-tier** (fast unit tests for CI,
Percolator gate run on the VM bench).

## Phase 1 — Bisect infrastructure

**Deliverable:** `scripts/bisect-score-psm.sh` committed to
`fix/score-psm-undercount`.

The script is invoked once per commit by `git bisect run`. Per call:

1. `cargo build --release --bin msgf-rust` — on failure, `exit 125`
   (git-bisect "skip this commit").
2. `/usr/bin/time -l <binary> --spectrum <PXD001819 mzML> --database
   <fasta> --output-pin /tmp/bisect.pin --precursor-tol-ppm 5
   --isotope-error-min=0 --isotope-error-max=1 --threads 1 --top-n 1
   --max-spectra 30000`
   Single-threaded for determinism (no rayon thread-count nondeterminism);
   `--max-spectra 30000` bounds wall time per step (scan=28787 is well
   inside the first 30k spectra).
3. `awk -F'\t' 'NR>1 && $3==28787 {print $7}' /tmp/bisect.pin >
   /tmp/bisect.score` — RawScore lives in column 7.
4. Append `<commit-sha>,<rawscore>` to `/tmp/bisect-trace.csv` (the
   per-step diagnostic log retained across bisect steps).
5. Exit code:
   - `0` (good) if RawScore ≥ 290
   - `1` (bad) if RawScore < 200
   - `125` (skip) if the pin is missing the scan or build failed

The two-band threshold (200 / 290) gives a safety margin: known-good is
~297, known-bad is ~108, so 200 cleanly separates them.

**Speed estimate:** ~7 min per step on Mac with `--threads 1` and
`--max-spectra 30000`. Bisect range of ~30 commits → 5-6 steps → ~45 min
total wall time.

## Phase 2 — Bisect execution

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed-score-fix
git bisect start
git bisect bad 68de6a7                              # current rust-implement HEAD
git bisect good <commit-from-May-10-parity-region>  # found by inspecting git log
git bisect run scripts/bisect-score-psm.sh
```

Finding the "good" commit: search `git log --until="2026-05-10 23:59"
--first-parent rust-implement` and find a commit whose state matches the
2026-05-10 memory (Rust 14,839 @ 1% FDR on PXD001819). Most likely
`5d912fc` ("docs: GF tails iter 2 closed - SP-vs-SP parity at 1.0 OOM")
or its immediate parent. Verify by running the bisect script ONCE on
that commit and confirming exit code 0 before starting `git bisect run`.

Save the CSV trace into the fix branch as
`docs/parity-analysis/notes/2026-05-13-score-bisect-trace.csv` for audit.

## Phase 3 — Diagnose the bad commit

For the identified commit:

1. `git show <bad-sha>` — read the full diff
2. Re-read the affected hunks against the forensic doc's three hypotheses
   (FastScorer cache, `directional_node_score_inner` ion-type subset,
   `parent_mass`/segment selection). Identify which hypothesis matches.
3. Write a one-sentence FALSIFIABLE hypothesis: *"The bug is that
   [specific line in commit X] introduces [specific behavior] which causes
   `score_psm` to undercount by N for splits where [property]."*
4. Validate by reading surrounding code: trace what the function did
   pre-commit vs post-commit for the `IVNEEFDQLEEDTPVYK` case.
5. If the diff doesn't obviously explain the 3× score drop, **ESCALATE**.
   Two escalation paths:
   - **(a)** Add per-split `eprintln!` instrumentation inside `score_psm`
     to surface where in the sum the gap appears.
   - **(b)** Run Java with `-Dmsgfplus.trace=true
     -Dmsgfplus.trace.scan=28787
     -Dmsgfplus.trace.pep=IVNEEFDQLEEDTPVYK` to get per-split Java scores;
     diff against Rust's per-split.

If the bisect CSV reveals a STAIRSTEP pattern (multiple contributing
commits, each adding some part of the gap), repeat Phase 3 for each
contributor.

**Acceptance to exit Phase 3:** a written one-paragraph diagnosis pinning
the bug to specific source lines + a falsifiable test distinguishing
"fixed" from "not fixed" (more specific than the bisect oracle).

## Phase 4 — Surgical fix

Single commit on `fix/score-psm-undercount` containing only the bug fix
(no test changes — those are Phase 5). Constraints:

- Touches the SAME file(s) as the bad commit
- Preserves the bad commit's INTENT (e.g., if bad commit added a perf
  optimization, the fix keeps the perf path but corrects the score
  arithmetic)
- Commit message names the bad commit SHA and quotes the Java reference
- "Surgical" means scoped to fixing the identified bug, not necessarily
  one-liner. A fix that adds a missing ion-type iteration or restructures
  cache population to match on-demand semantics still counts as surgical
  if it's the minimum change needed to restore correctness.

**Acceptable shapes:**

- Replace a constant (e.g., wrong array bound, wrong loop limit)
- Add a missing ion-type iteration to `directional_node_score_inner`
- Fix the cache population path to match on-demand semantics

**Not acceptable:**

- Refactoring beyond the fix
- "While we're at it" cleanups in adjacent code
- Reverting the bad commit wholesale (sacrifices perf wins)

## Phase 5 — Tests + validation

**New tests** committed in a SEPARATE commit immediately after the Phase-4
fix commit (so bisecting future regressions is unambiguous — the fix and
the regression guards are distinct in git history):

1. **`score_psm_pxd001819_scan_28787_matches_java`** in
   `rust/crates/scoring/src/scoring/psm_score.rs` test module. Loads a
   small fixture (commit a scan=28787-only mzML extracted from PXD001819
   into `benchmark/fixtures/`, plus the existing PXD001819 fasta is
   available via test resources) or uses the full PXD001819 mzML if test
   speed allows. Runs `score_psm` for `IVNEEFDQLEEDTPVYK` at charge 2.
   Asserts `RawScore == 297` (Java's value at the time the test is
   written).
2. **Sister tests** on 2-3 other worst-gap scans (28825, 33606, 32395) to
   detect the bug class, not just the one example. If all four are fixed
   by the single Phase-4 change, confidence is high.

The fixture format question (full mzML vs scan-extracted slice) is
decided in writing-plans based on test-runtime budget; either works for
correctness.

**Validation gates** (must all pass before merging into `rust-implement`):

- `cargo build --workspace --release` clean
- `cargo test --workspace --lib` all tests passing (318 current + the
  new regression tests from this iteration — count grows with how many
  sister-scan assertions land)
- `cargo test -p search --test gf_java_parity --release` passes at 1.0
  OOM tolerance
- VM bench: PXD001819 full search + Percolator @ 1% FDR ≥ **14,800**
  (~99% of Java's 14,989)
- VM bench: Astral full search + Percolator @ 1% FDR ≥ **33,000** (~92%
  of Java's 35,818)
- VM bench: TMT full search + Percolator @ 1% FDR ≥ **10,500** (no
  regression from current 10,548)

If Astral gate fails but PXD001819 passes, the bug had a dataset-specific
second cause. Document the residual gap and decide as a separate
iteration — but the per-scan diagnostic suggests the same `score_psm`
path runs on both, so a single fix is expected to cover both datasets.

## Phase 6 — Cleanup

1. Update the forensic doc
   `docs/parity-analysis/reports/2026-05-13-score-psm-undercount-finding.md`
   with a new `## Resolution` section: bad commit SHA, root cause
   one-paragraph, fix commit SHA, post-fix Percolator numbers.
2. Final review pass: remove any `eprintln!` left from Phase-3
   instrumentation; remove commented-out exploratory code; verify the
   bisect CSV is committed as an audit artifact.
3. Merge `fix/score-psm-undercount` into `rust-implement` via
   `git merge --no-ff` with a merge message summarizing the bench gates
   and PSM recovery numbers.

## Data flow

```
rust-implement HEAD (68de6a7, broken on PXD001819/Astral)
       │
       ├──> worktree astral-speed-score-fix, branch fix/score-psm-undercount
       │
       ├──> Phase 1: scripts/bisect-score-psm.sh committed
       │
       ├──> Phase 2: git bisect run -> bad commit SHA identified;
       │            /tmp/bisect-trace.csv -> docs/parity-analysis/notes/
       │
       ├──> Phase 3: diagnosis paragraph + falsifiable test specification
       │
       ├──> Phase 4: surgical fix commit
       │
       ├──> Phase 5: regression-test commit(s); VM bench validation
       │
       ├──> Phase 6: forensic doc update + merge --no-ff
       │
       v
rust-implement HEAD (fixed; Percolator gates satisfied)
```

## Error handling / escalation

| Situation | Escalation |
|---|---|
| Build fails at multiple intermediate commits | Verify `cargo --version` consistent with `rust-toolchain.toml`; skip via exit 125 |
| Bisect lands on a "merge" commit | Re-bisect on linear ancestry (`git bisect start --first-parent`) |
| Bisect lands on a commit whose diff is opaque | Escalate to Phase-3 path (a) or (b): per-split instrumentation or Java trace |
| Per-step RawScore CSV shows stairstep (multiple contributors) | Run Phase 3 for each contributor; each may need its own surgical fix commit |
| Fix passes scan=28787 unit test but Astral Percolator still fails | Re-run diagnostic on a worst-gap Astral scan; treat as second-iteration scope |
| Fix passes all but TMT regresses | Revert the fix; the change interacts with the candidate-gen fixes from this session — escalate to a wider design |

## Testing

Two tiers, gating different things:

**Tier 1 — CI-fast (runs in `cargo test --workspace --lib`):**
- `score_psm_pxd001819_scan_28787_matches_java` (primary regression gate)
- 2-3 sister tests on other worst-gap scans

These run in seconds; they're the permanent guard against this bug class
recurring.

**Tier 2 — VM bench (runs manually, not CI):**
- `PXD001819 / Astral / TMT` end-to-end → Percolator @ 1% FDR thresholds

This tier validates the whole pipeline including the score → SpecEValue →
Percolator-feature cascade. Run before merging the fix.

## Scope boundaries

**In scope:**
- The `score_psm` under-scoring bug
- Minimal supporting changes inside the bad commit's file(s)
- One or more regression tests
- The forensic doc update
- The fix-branch merge

**Explicitly out of scope:**
- Other column-divergence bugs flagged in the forensic doc
  (`MS2IonCurrent` scale, `ExplainedIonCurrentRatio`, `MeanErrorTop7`
  scale, etc.). They contribute to Percolator's gap but are independent
  bugs. Note their existence; do not fix them here.
- The `lnDeltaSpecEValue` zero-stub (related to `top-N=1` design choice;
  separate iteration).
- The `i32::MIN` sentinel write-path leak. It's a SECONDARY symptom of
  the `score_psm` under-scoring (when score is wrong, SpecE early-returns,
  sentinels leak). Once `score_psm` is fixed, the SpecE early-return
  stops firing and the sentinels disappear naturally. If after Phase 5
  any PSMs still leak sentinels, treat as a separate finding.
- Any perf optimization of the fix itself.

## Risk register

| Risk | Mitigation |
|---|---|
| Bisect determinism breaks (different threads = different scores) | `--threads 1` + verify with two consecutive runs at HEAD |
| Bad commit is a merge with many sub-commits | `git bisect start --first-parent`; re-bisect within the offending sub-branch if needed |
| Surgical fix breaks the bad commit's perf intent | Re-run a quick perf check (hot-path timing on PXD001819, 8 threads) before merging |
| Astral fails Percolator gate even with PXD001819 fixed | Run focused diagnostic on the worst-gap Astral scan; defer to a follow-up iteration if the cause is independent |
| Bisect CSV shows stairstep pattern | Address each contributor in its own commit; not all may be salvageable surgically |
| Fix introduces a different correctness regression | The `gf_java_parity` 1.0 OOM gate + 318 lib tests provide a wide safety net |

## Out-of-band notes

- This work happens on `fix/score-psm-undercount`. The branch will be
  merged into `rust-implement` after all Phase 5 gates pass. **Nothing on
  this branch ships to any external release.** External shipping is gated
  on the broader "rust-implement is release-quality" decision the user
  has not yet made.
- The forensic doc's `## What's safe to ship` section is the standing
  guidance: nothing from this session is shipped externally until
  `score_psm` is rooted and fixed AND the Percolator gates close.
