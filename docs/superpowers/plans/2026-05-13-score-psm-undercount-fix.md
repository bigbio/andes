# `score_psm` under-scoring fix — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Find and fix the `score_psm` regression that causes Rust to compute RawScore ≈ 1/3 of Java's value for identical (peptide, scan, charge) inputs on PXD001819 / Astral; recover Percolator @ 1% FDR within 5% of Java baselines while preserving the TMT win and the recent perf optimizations.

**Architecture:** Bisect-driven discovery on a dedicated worktree. A deterministic bisect-oracle script (single-threaded msgf-rust on PXD001819, gates on scan=28787's RawScore) walks the commit range between known-parity (2026-05-10) and HEAD; once the bad commit is identified, a surgical fix is applied. Permanent regression test asserts the specific scan's RawScore; VM Percolator bench gates the merge.

**Tech Stack:** Rust 2021 (workspace at `rust/crates/`, Cargo, rustc 1.95.0 toolchain pin via `rust-toolchain.toml`), bash for the bisect oracle, awk for mzML slicing + pin parsing, Percolator via biocontainers Docker on the VM.

**Forensic context:** `docs/parity-analysis/reports/2026-05-13-score-psm-undercount-finding.md` — read this first if not already.

**Design spec:** `docs/superpowers/specs/2026-05-13-score-psm-undercount-fix-design.md` — the design this plan implements.

---

## File Structure

Files this plan creates or modifies (locked-in decomposition):

- **Create:** `scripts/bisect-score-psm.sh` — the bisect oracle (executable bash)
- **Create:** `src/test/resources/benchmark/PXD001819/scan_28787.mzML` — extracted single-scan fixture (~5–50 KB)
- **Create:** `rust/crates/scoring/tests/score_psm_pxd001819_parity.rs` — the regression test file
- **Create:** `docs/parity-analysis/notes/2026-05-13-score-bisect-trace.csv` — per-commit RawScore log
- **Create:** `docs/parity-analysis/notes/2026-05-13-score-fix-diagnosis.md` — Phase-3 diagnosis paragraph
- **Modify:** One Rust source file inside `rust/crates/scoring/src/` (file determined by bisect outcome — likely `scored_spectrum.rs` or `psm_score.rs`)
- **Modify:** `docs/parity-analysis/reports/2026-05-13-score-psm-undercount-finding.md` — add `## Resolution` section

The fix file is the only unknown; everything else is locked.

---

## Task 1: Worktree + branch setup

**Files:**
- No code changes — git operations only

- [ ] **Step 1: Create the fix worktree from rust-implement**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed
git worktree add -b fix/score-psm-undercount ../astral-speed-score-fix rust-implement
```

Expected: worktree at `/Users/yperez/work/msgfplus-workspace/astral-speed-score-fix`, branch `fix/score-psm-undercount`.

- [ ] **Step 2: Confirm clean state and base commit**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed-score-fix
git status --short
git log --oneline -1
```

Expected: empty `git status` output (no dirty files in fresh worktree); HEAD commit message starts with `docs(spec): score_psm under-scoring fix design` (the spec commit on rust-implement).

- [ ] **Step 3: Verify the toolchain pin**

```bash
cat rust/rust-toolchain.toml
```

Expected: contains `channel = "1.95.0"`. If not, the bisect script must override per-step.

---

## Task 2: Bisect oracle script

**Files:**
- Create: `scripts/bisect-score-psm.sh`

- [ ] **Step 1: Write the failing test for the oracle (smoke test only — verify the script's logic)**

The "failing test" here is structural: confirm the oracle correctly identifies HEAD as BAD (exit 1) before we commit it. The oracle is itself the test for every later step.

- [ ] **Step 2: Create `scripts/` directory and write the oracle**

```bash
mkdir -p /Users/yperez/work/msgfplus-workspace/astral-speed-score-fix/scripts
```

Write `scripts/bisect-score-psm.sh`:

```bash
#!/usr/bin/env bash
# Bisect oracle for the score_psm under-scoring regression.
#
# - Builds msgf-rust at the current commit
# - Runs it on PXD001819 single-threaded with --max-spectra 30000
# - Greps scan=28787's RawScore from the pin (column 7)
# - Appends <sha>,<rawscore> to /tmp/bisect-trace.csv (cumulative log)
# - Exits 0 (good) if RawScore >= 290
# - Exits 1 (bad)  if RawScore <  200
# - Exits 125 (skip) on build failure or missing scan in pin
#
# Determinism: --threads 1 eliminates rayon nondeterminism. The same
# commit produces the same RawScore across runs.

set -uo pipefail

REPO_ROOT="/Users/yperez/work/msgfplus-workspace/astral-speed-score-fix"
PXD_MZML="/Users/yperez/work/msgfplus-workspace/benchmark/data/PXD001819/UPS1_5000amol_R1.mzML"
PXD_FASTA="/Users/yperez/work/msgfplus-workspace/benchmark/data/PXD001819/PXD001819_uniprot_yeast_ups.fasta"
TRACE_CSV="/tmp/bisect-trace.csv"
PIN_OUT="/tmp/bisect.pin"

cd "$REPO_ROOT/rust"
SHA=$(git rev-parse --short HEAD)

# Skip non-existent inputs (would lead to false bad).
if [ ! -f "$PXD_MZML" ] || [ ! -f "$PXD_FASTA" ]; then
    echo "[$SHA] missing PXD001819 fixture — skip"
    exit 125
fi

# Build. Use full build (not --quiet) so cargo errors are visible in
# `git bisect run` logs.
if ! cargo build --release --bin msgf-rust 2>&1 | tail -5; then
    echo "[$SHA] build failed — skip"
    echo "$SHA,BUILD_FAIL" >> "$TRACE_CSV"
    exit 125
fi

BIN="$REPO_ROOT/rust/target/release/msgf-rust"
rm -f "$PIN_OUT"

if ! "$BIN" \
        --spectrum "$PXD_MZML" \
        --database "$PXD_FASTA" \
        --output-pin "$PIN_OUT" \
        --precursor-tol-ppm 5 \
        --isotope-error-min=0 \
        --isotope-error-max=1 \
        --top-n 1 \
        --threads 1 \
        --max-spectra 30000 \
        > /tmp/bisect.log 2>&1; then
    echo "[$SHA] msgf-rust run failed — skip"
    echo "$SHA,RUN_FAIL" >> "$TRACE_CSV"
    exit 125
fi

# Column 7 of the pin is RawScore.
RAW=$(awk -F'\t' 'NR>1 && $3 == 28787 {print $7; exit}' "$PIN_OUT")

if [ -z "$RAW" ]; then
    echo "[$SHA] scan=28787 not in pin output — skip"
    echo "$SHA,MISSING_SCAN" >> "$TRACE_CSV"
    exit 125
fi

echo "$SHA,$RAW" >> "$TRACE_CSV"
echo "[$SHA] scan=28787 RawScore=$RAW"

if [ "$RAW" -ge 290 ] 2>/dev/null; then
    exit 0  # good
fi
if [ "$RAW" -lt 200 ] 2>/dev/null; then
    exit 1  # bad
fi

# In the dead-band 200..290: skip to avoid mis-bisecting on intermediate.
echo "[$SHA] RawScore=$RAW in dead band 200..290 — skip"
exit 125
```

- [ ] **Step 3: Make it executable**

```bash
chmod +x /Users/yperez/work/msgfplus-workspace/astral-speed-score-fix/scripts/bisect-score-psm.sh
```

- [ ] **Step 4: Smoke-test at HEAD (must exit 1, RawScore < 200)**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed-score-fix
rm -f /tmp/bisect-trace.csv
./scripts/bisect-score-psm.sh
echo "exit=$?"
cat /tmp/bisect-trace.csv
```

Expected: prints `[<sha>] scan=28787 RawScore=<some number near 108>` and `exit=1`. CSV has one row with the score. If RawScore is in 200..290 dead band: investigate variance, tighten or widen thresholds.

If exit is NOT 1, the bug is no longer present at HEAD, or the oracle is broken. Halt and re-diagnose.

- [ ] **Step 5: Commit the oracle**

```bash
git add scripts/bisect-score-psm.sh
git commit -m "tools: bisect oracle for score_psm under-scoring regression"
```

---

## Task 3: Identify the known-good baseline commit

**Files:**
- No new files — git inspection only

- [ ] **Step 1: List commits in the May 10–12 window**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed-score-fix
git log --oneline --since="2026-05-09" --until="2026-05-13" --first-parent rust-implement | head -30
```

Expected: a list of commits including `5d912fc docs: GF tails iter 2 closed - SP-vs-SP parity at 1.0 OOM`, `0af1a37 perf(scoring): Track A — FastScorer prefix/suffix score cache per spectrum`, `be50dab perf(search): hoist compute_psm_features to post-top-N finalization`, etc.

- [ ] **Step 2: Reset to the candidate known-good commit**

```bash
git checkout 5d912fc
git log --oneline -1
```

Expected: HEAD detached at `5d912fc`.

- [ ] **Step 3: Run the bisect oracle on this commit (must exit 0)**

```bash
rm -f /tmp/bisect-trace.csv
./scripts/bisect-score-psm.sh
echo "exit=$?"
```

Expected: `[5d912fc] scan=28787 RawScore=<value near 297>` and `exit=0`. If exit is NOT 0, this commit is also broken — walk further back (`git log --oneline --until="2026-05-09"` then test each candidate).

- [ ] **Step 4: Return to fix branch HEAD**

```bash
git checkout fix/score-psm-undercount
```

Expected: clean checkout, HEAD on fix/score-psm-undercount.

- [ ] **Step 5: Note the confirmed known-good commit** in a working note. Just remember the SHA for Task 4. No commit needed.

---

## Task 4: Run the bisect

**Files:**
- Create: `docs/parity-analysis/notes/2026-05-13-score-bisect-trace.csv` (output)

- [ ] **Step 1: Start git bisect with the confirmed range**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed-score-fix
rm -f /tmp/bisect-trace.csv
git bisect start
git bisect bad HEAD                  # fix/score-psm-undercount = rust-implement HEAD
git bisect good 5d912fc              # OR the SHA confirmed good in Task 3
```

Expected: git reports "Bisecting: N revisions left to test" — N should be 25–35.

- [ ] **Step 2: Run the bisect automatically**

```bash
git bisect run ./scripts/bisect-score-psm.sh 2>&1 | tee /tmp/bisect-run.log
```

Expected: git iterates through ~5–6 commits over ~30–45 minutes; finishes with output like:
```
<sha> is the first bad commit
```
plus the commit message and diff stat of the offending commit.

If git reports "merge-base" issues or skips many commits, see escalation paths in Task 4 Step 5.

- [ ] **Step 3: Save the per-step trace to the audit doc**

```bash
mkdir -p docs/parity-analysis/notes
cp /tmp/bisect-trace.csv docs/parity-analysis/notes/2026-05-13-score-bisect-trace.csv
echo "--- bisect run log ---" >> docs/parity-analysis/notes/2026-05-13-score-bisect-trace.csv
cat /tmp/bisect-run.log | grep -E "^(\[|<|Bisecting|first bad)" >> docs/parity-analysis/notes/2026-05-13-score-bisect-trace.csv
```

- [ ] **Step 4: Reset bisect state and return to fix branch**

```bash
git bisect reset
git checkout fix/score-psm-undercount
```

Expected: clean state on `fix/score-psm-undercount`.

- [ ] **Step 5: Escalation paths if bisect blocked**

If `git bisect run` reports a confusing result:

| Symptom | Action |
|---|---|
| Many commits skipped (compile failures) | Re-run bisect with `git bisect start --first-parent` so only merge-base commits are tested |
| Bad commit is a merge commit | Re-bisect within the merged sub-branch: `git bisect start <merge-sha>^2 <merge-sha>^1` |
| All commits exit 125 (skip) | The oracle has a bug; debug `scripts/bisect-score-psm.sh` |
| Bad commit's diff is non-obvious | Continue to Task 5 escalation paths |

- [ ] **Step 6: Commit the trace CSV**

```bash
git add docs/parity-analysis/notes/2026-05-13-score-bisect-trace.csv
git commit -m "diag: bisect trace for score_psm under-scoring regression"
```

---

## Task 5: Diagnose the bad commit

**Files:**
- Create: `docs/parity-analysis/notes/2026-05-13-score-fix-diagnosis.md` — diagnosis paragraph

- [ ] **Step 1: Read the bad commit's diff**

```bash
BAD_SHA=$(grep "is the first bad commit" /tmp/bisect-run.log | awk '{print $1}')
git show "$BAD_SHA" 2>&1 | head -200
```

Expected: full diff of the identified bad commit. Read it line by line.

- [ ] **Step 2: Cross-reference the three hypotheses from the forensic doc**

| Hypothesis | Where to look |
|---|---|
| FastScorer cache | `rust/crates/scoring/src/scoring/scored_spectrum.rs`: `prefix_score_cache` / `suffix_score_cache` population loop |
| `directional_node_score_inner` ion-type subset | same file: the per-ion iteration loop inside `directional_node_score_inner` |
| `parent_mass` / segment selection | `rust/crates/scoring/src/scoring/scored_spectrum.rs::new` and `segment_partition_cache` |

If the bad commit touches one of these, that's almost certainly your hypothesis.

- [ ] **Step 3: Form a one-paragraph falsifiable diagnosis**

Open `docs/parity-analysis/notes/2026-05-13-score-fix-diagnosis.md` and write:

```markdown
# score_psm under-scoring — diagnosis

**Bad commit:** `<full-sha>`
**Title:** `<commit subject line>`
**Date:** `<commit date>`
**Files changed:** `<list>`

## Root cause hypothesis

[One paragraph explaining exactly what the bad commit did wrong. Use
language like "the change at file:line introduces X, which causes
score_psm to undercount for splits where Y." If multiple lines
contribute, name each.]

## Falsifiable test

This hypothesis is true if and only if: applying [described fix]
restores RawScore for scan=28787 + IVNEEFDQLEEDTPVYK to 297 (Java
baseline) without breaking gf_java_parity.

## Surface of the fix

Files I'll modify in Task 6:
- `<path>` — `<line range or function name>`

Expected change shape: [add missing iteration | restore correct array
bound | fix the cache-population path to match on-demand semantics].
```

- [ ] **Step 4: Escalate if the diff is opaque**

If Step 3 cannot be filled in because the diff doesn't obviously explain
a 3× score drop, escalate before Task 6:

**(a) Per-split instrumentation at HEAD.** Add eprintln inside `score_psm` printing the per-split prefix and suffix scores. Run on PXD001819 + IVNEEFDQLEEDTPVYK + scan=28787. Capture output. Repeat for the bad commit's predecessor (known-good). Diff per-split scores — the first divergent split is the bug entry point. (Remove the eprintlns before Task 6.)

**(b) Java reference dump.** Build the Java MS-GF+ JAR with `-Dmsgfplus.trace=true -Dmsgfplus.trace.scan=28787 -Dmsgfplus.trace.pep=IVNEEFDQLEEDTPVYK` enabled (already wired in `FastScorer.java` per past session work). Run on the same data. Capture Java's per-split scores. Diff against Rust's.

Document the escalation outcome in the same diagnosis file under `## Escalation outcome`.

- [ ] **Step 5: Commit the diagnosis**

```bash
git add docs/parity-analysis/notes/2026-05-13-score-fix-diagnosis.md
git commit -m "diag: bad commit identified for score_psm under-scoring"
```

---

## Task 6: Add the failing regression test (fixture + test)

**Files:**
- Create: `src/test/resources/benchmark/PXD001819/scan_28787.mzML` (~5–50 KB)
- Create: `rust/crates/scoring/tests/score_psm_pxd001819_parity.rs`

- [ ] **Step 1: Extract scan=28787 from PXD001819 into a single-scan mzML**

```bash
awk '
  /<spectrum[^>]+scan=28787"/ { in_spec=1; print; next }
  in_spec { print }
  in_spec && /<\/spectrum>/    { in_spec=0; exit }
' /Users/yperez/work/msgfplus-workspace/benchmark/data/PXD001819/UPS1_5000amol_R1.mzML \
  > /tmp/scan_28787_inner.xml
```

Wrap it in a minimal mzML envelope:

```bash
cat > /Users/yperez/work/msgfplus-workspace/astral-speed-score-fix/src/test/resources/benchmark/PXD001819/scan_28787.mzML <<'HEADER'
<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">
  <run id="scan28787">
    <spectrumList count="1" defaultDataProcessingRef="dp">
HEADER
cat /tmp/scan_28787_inner.xml \
    >> /Users/yperez/work/msgfplus-workspace/astral-speed-score-fix/src/test/resources/benchmark/PXD001819/scan_28787.mzML
cat >> /Users/yperez/work/msgfplus-workspace/astral-speed-score-fix/src/test/resources/benchmark/PXD001819/scan_28787.mzML <<'FOOTER'
    </spectrumList>
  </run>
</mzML>
FOOTER
ls -la /Users/yperez/work/msgfplus-workspace/astral-speed-score-fix/src/test/resources/benchmark/PXD001819/scan_28787.mzML
```

Expected: file exists, size > 1 KB and < 200 KB.

- [ ] **Step 2: Write the failing regression test**

Create `rust/crates/scoring/tests/score_psm_pxd001819_parity.rs`:

```rust
//! Regression test: score_psm must match Java's RawScore for the
//! IVNEEFDQLEEDTPVYK / scan=28787 case on PXD001819. The bug fixed by
//! this branch caused Rust to compute 108 here while Java reports 297;
//! re-introducing the bug regresses Percolator @ 1% FDR on PXD001819 by
//! 22%, so this guard is load-bearing.
//!
//! See:
//!   docs/parity-analysis/reports/2026-05-13-score-psm-undercount-finding.md
//!   docs/parity-analysis/notes/2026-05-13-score-fix-diagnosis.md

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use model::{AminoAcid, Peptide, PrecursorTolerance, Tolerance};
use scoring::param_model::Param;
use scoring::scoring::{score_psm, RankScorer, ScoredSpectrum};
use input::MzMLReader;

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .canonicalize()
        .expect("workspace root resolves")
}

fn load_param() -> Param {
    let path = workspace_root()
        .join("src/main/resources/ionstat/HCD_QExactive_Tryp.param");
    Param::load_from_file(&path).expect("param loads")
}

fn load_scan_28787() -> model::Spectrum {
    let mzml = workspace_root()
        .join("src/test/resources/benchmark/PXD001819/scan_28787.mzML");
    let f = File::open(&mzml).expect("scan_28787.mzML exists");
    let reader = MzMLReader::new(BufReader::new(f));
    for result in reader {
        let s = result.expect("parses");
        if s.scan == Some(28787) {
            return s;
        }
    }
    panic!("scan=28787 not found in fixture");
}

fn ivneefdqleedtpvyk_peptide() -> Peptide {
    // K.IVNEEFDQLEEDTPVYK.L — 17 residues, no internal mods, fully tryptic.
    let residues: Vec<AminoAcid> = b"IVNEEFDQLEEDTPVYK"
        .iter()
        .map(|&r| AminoAcid::standard(r).expect("standard residue"))
        .collect();
    Peptide::new(residues, b'K', b'L')
}

/// Java's RawScore for this PSM is 297 (extracted from
/// bench-merged-results/pxd001819-java.pin row for scan=28787 column 7).
/// We accept ±10 since rank-score scaling rounds at integer boundaries.
const EXPECTED_JAVA_RAWSCORE: i32 = 297;
const TOLERANCE: i32 = 10;

#[test]
fn score_psm_pxd001819_scan_28787_matches_java() {
    let param = load_param();
    let scorer = RankScorer::new(&param);
    let spec = load_scan_28787();

    let scored_spec = ScoredSpectrum::new(&spec, &scorer, 2);
    let peptide = ivneefdqleedtpvyk_peptide();

    let score = score_psm(&scored_spec, &peptide, &scorer, 2, 0.5).round() as i32;

    assert!(
        (score - EXPECTED_JAVA_RAWSCORE).abs() <= TOLERANCE,
        "score_psm regressed: got RawScore={} but Java reports {} ± {} \
         for scan=28787 IVNEEFDQLEEDTPVYK charge=2",
        score, EXPECTED_JAVA_RAWSCORE, TOLERANCE
    );
}
```

- [ ] **Step 3: Run the failing test to confirm it fails BEFORE the fix**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed-score-fix/rust
cargo test -p scoring --test score_psm_pxd001819_parity 2>&1 | tail -20
```

Expected: FAIL with the assert message showing `RawScore=108` (or whatever the bug produces) `but Java reports 297 ± 10`.

If the test ERRORS instead of FAILS (e.g., panic in `load_scan_28787` because the fixture doesn't parse), the fixture extraction in Step 1 has a problem. Re-run Step 1 with verbose logging.

- [ ] **Step 4: Commit the fixture and the failing test (red state)**

```bash
git add src/test/resources/benchmark/PXD001819/scan_28787.mzML \
        rust/crates/scoring/tests/score_psm_pxd001819_parity.rs
git commit -m "test(scoring): regression guard for score_psm scan=28787 IVNEEFDQLEEDTPVYK

Asserts RawScore == 297 ± 10 (Java baseline). Currently FAILING — the
score_psm under-scoring bug yields ~108. Fix lands in the next commit;
the test stays as a permanent regression guard."
```

---

## Task 7: Apply the surgical fix

**Files:**
- Modify: one file in `rust/crates/scoring/src/` identified by Task 5's diagnosis

This task's exact code depends on what Task 5 identified. Task 5's
diagnosis file specifies the file and the change shape; follow that
diagnosis.

- [ ] **Step 1: Re-read the diagnosis**

```bash
cat docs/parity-analysis/notes/2026-05-13-score-fix-diagnosis.md
```

Confirm the diagnosis is concrete: a specific file, specific lines, and
a specific change shape.

- [ ] **Step 2: Apply the surgical fix**

Open the file named in the diagnosis. Apply the minimum change required
by the diagnosis. Common patterns:

```rust
// Pattern A: missing ion-type iteration (Hypothesis B from forensic doc)
// BEFORE: iterates only Prefix/Suffix; MISSING b2+, y2+, etc.
// AFTER:
for (ion, logs) in ion_logs {
    // ... existing body ...
}
// PLUS the previously-missing ion types from the partition.
```

```rust
// Pattern B: cache populated with wrong inputs (Hypothesis A)
// BEFORE:
prefix_score_cache[nominal_mass] =
    Self::directional_node_score_inner(spec, &ranks, ..., true, ...);
// AFTER (example): ensure ranks is fully populated BEFORE cache fill
// or use the on-demand path's exact arg list.
```

```rust
// Pattern C: wrong array bound or off-by-one in cache_len
// BEFORE:
let cache_len = (nominal_from(parent_mass).max(0) as usize) + 1;
// AFTER (example): cache_len uses the correct ceiling derived from
// parent_mass + max_mod_delta, matching Java's array sizing.
```

Use the diagnosis's exact prescribed change. Do not invent fixes.

- [ ] **Step 3: Run the regression test — must now pass**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed-score-fix/rust
cargo test -p scoring --test score_psm_pxd001819_parity 2>&1 | tail -10
```

Expected: 1 passed.

If still failing: the diagnosis is wrong or incomplete. Return to Task 5.

- [ ] **Step 4: Run the full workspace lib tests**

```bash
cargo test --workspace --lib 2>&1 | grep "^test result" | awk '{c+=$4; f+=$6} END {print "passed: "c", failed: "f}'
```

Expected: failed=0; passed ≥ 318 (the prior count). New test from Task 6 is in `--test` not `--lib`, so the lib count stays at 318.

If any lib test fails: the fix has a side effect outside the diagnosed surface. Diagnose and either tighten the fix or revert and re-diagnose.

- [ ] **Step 5: Run gf_java_parity at release**

```bash
cargo test -p search --test gf_java_parity --release 2>&1 | tail -5
```

Expected: 6 passed. If this fails: the fix breaks SP parity. Revert and re-diagnose.

- [ ] **Step 6: Commit the fix**

```bash
git add rust/crates/scoring/src/<the-modified-file>.rs
git commit -m "fix(scoring): score_psm under-scoring regression on PXD001819/Astral

Restores RawScore parity with Java for the IVNEEFDQLEEDTPVYK case (scan
28787 on PXD001819): was 108, now 297 (Java baseline).

Root cause: <one-line summary from diagnosis>. Bad commit was <bad-sha>
('<bad-commit-subject>'). The perf intent of the bad commit is preserved
— this fix targets only the specific defect.

Test: tests/score_psm_pxd001819_parity.rs (committed in the prior
commit) now passes. 318 lib tests pass; gf_java_parity passes at 1.0
OOM tolerance.

Forensic context:
  docs/parity-analysis/reports/2026-05-13-score-psm-undercount-finding.md
Diagnosis:
  docs/parity-analysis/notes/2026-05-13-score-fix-diagnosis.md"
```

---

## Task 8: Add sister-scan regression tests (defense in depth)

**Files:**
- Modify: `rust/crates/scoring/tests/score_psm_pxd001819_parity.rs`
- Modify: `src/test/resources/benchmark/PXD001819/scan_28787.mzML` → extend to also include scans 28825, 33606, 32395

- [ ] **Step 1: Extract three more worst-gap scans into an extended fixture**

```bash
EXTRACT() {
    local SCAN=$1
    awk -v target="$SCAN" '
      $0 ~ "<spectrum[^>]+scan=" target "\""  { in_spec=1; print; next }
      in_spec { print }
      in_spec && /<\/spectrum>/                { in_spec=0; exit }
    ' /Users/yperez/work/msgfplus-workspace/benchmark/data/PXD001819/UPS1_5000amol_R1.mzML
}
{
    cat <<'HEADER'
<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">
  <run id="worstgap">
    <spectrumList count="4" defaultDataProcessingRef="dp">
HEADER
    EXTRACT 28787
    EXTRACT 28825
    EXTRACT 33606
    EXTRACT 32395
    cat <<'FOOTER'
    </spectrumList>
  </run>
</mzML>
FOOTER
} > /Users/yperez/work/msgfplus-workspace/astral-speed-score-fix/src/test/resources/benchmark/PXD001819/scan_28787.mzML
```

Note: we keep the same filename (`scan_28787.mzML`) but now it has 4 scans. Lower friction than introducing a second fixture file.

- [ ] **Step 2: Extend the test file with three sister tests**

In `rust/crates/scoring/tests/score_psm_pxd001819_parity.rs`, append after the existing test:

```rust
fn load_scan(target: i32) -> model::Spectrum {
    let mzml = workspace_root()
        .join("src/test/resources/benchmark/PXD001819/scan_28787.mzML");
    let f = File::open(&mzml).expect("scan_28787.mzML exists");
    let reader = MzMLReader::new(BufReader::new(f));
    for result in reader {
        let s = result.expect("parses");
        if s.scan == Some(target) {
            return s;
        }
    }
    panic!("scan={} not found in fixture", target);
}

fn peptide_from(seq: &[u8], pre: u8, post: u8) -> Peptide {
    let residues: Vec<AminoAcid> = seq
        .iter()
        .map(|&r| AminoAcid::standard(r).expect("standard residue"))
        .collect();
    Peptide::new(residues, pre, post)
}

#[test]
fn score_psm_pxd001819_scan_28825_matches_java() {
    // K.IVNEEFDQLEEDTPVYK.L (same peptide as 28787, different scan).
    // Java RawScore=305.
    let param = load_param();
    let scorer = RankScorer::new(&param);
    let spec = load_scan(28825);
    let scored_spec = ScoredSpectrum::new(&spec, &scorer, 2);
    let peptide = peptide_from(b"IVNEEFDQLEEDTPVYK", b'K', b'L');
    let score = score_psm(&scored_spec, &peptide, &scorer, 2, 0.5).round() as i32;
    assert!((score - 305).abs() <= TOLERANCE,
        "scan=28825 IVNEEFDQLEEDTPVYK: got {}, expected 305 ± {}", score, TOLERANCE);
}

#[test]
fn score_psm_pxd001819_scan_33606_matches_java() {
    // R.LESYVASIEQTVTDPVLSSK.L — Java RawScore=318.
    let param = load_param();
    let scorer = RankScorer::new(&param);
    let spec = load_scan(33606);
    let scored_spec = ScoredSpectrum::new(&spec, &scorer, 2);
    let peptide = peptide_from(b"LESYVASIEQTVTDPVLSSK", b'R', b'L');
    let score = score_psm(&scored_spec, &peptide, &scorer, 2, 0.5).round() as i32;
    assert!((score - 318).abs() <= TOLERANCE,
        "scan=33606 LESYVASIEQTVTDPVLSSK: got {}, expected 318 ± {}", score, TOLERANCE);
}

#[test]
fn score_psm_pxd001819_scan_32395_matches_java() {
    // R.AVGSLTFDENYNLLDTSGVAK.V — Java RawScore=329.
    let param = load_param();
    let scorer = RankScorer::new(&param);
    let spec = load_scan(32395);
    let scored_spec = ScoredSpectrum::new(&spec, &scorer, 2);
    let peptide = peptide_from(b"AVGSLTFDENYNLLDTSGVAK", b'R', b'V');
    let score = score_psm(&scored_spec, &peptide, &scorer, 2, 0.5).round() as i32;
    assert!((score - 329).abs() <= TOLERANCE,
        "scan=32395 AVGSLTFDENYNLLDTSGVAK: got {}, expected 329 ± {}", score, TOLERANCE);
}
```

- [ ] **Step 3: Run the four-test file — all should pass**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed-score-fix/rust
cargo test -p scoring --test score_psm_pxd001819_parity 2>&1 | tail -10
```

Expected: 4 passed.

If any sister test fails: the fix in Task 7 was specific to scan=28787 and didn't cover the bug class. Re-diagnose for that scan and extend the fix.

- [ ] **Step 4: Commit the sister tests**

```bash
git add rust/crates/scoring/tests/score_psm_pxd001819_parity.rs \
        src/test/resources/benchmark/PXD001819/scan_28787.mzML
git commit -m "test(scoring): sister-scan regression guards (28825, 33606, 32395)

Three additional PXD001819 scans with worst-case RawScore gaps. All four
tests now share a 4-scan mzML fixture. Confirms the score_psm fix covers
the bug class, not just the one example."
```

---

## Task 9: VM Percolator validation

**Files:**
- No code changes — bench + Percolator runs only

- [ ] **Step 1: Confirm SSH master socket is up**

```bash
ssh -O check -S /tmp/msgfplus-bench.sock root@pride-linux-vm
```

Expected: `Master running (pid=...)`. If not, ask the user to re-establish.

- [ ] **Step 2: Tarball + ship the worktree source**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed-score-fix
tar --exclude='rust/target' --exclude='.git' --exclude='target' \
    --exclude='benchmark/data' --exclude='benchmark/results' \
    --exclude='benchmark/baseline*' --exclude='benchmark/current-dev' \
    --exclude='benchmark/lib' --exclude='benchmark/plots' \
    --exclude='benchmark/parity-fixtures' --exclude='benchmark/__pycache__' \
    --exclude='benchmark/new' \
    -czf /tmp/msgf-rust-scorefix.tgz rust src/main/resources
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm \
    "cat > /srv/data/msgf-bench/msgf-rust-scorefix.tgz" \
    < /tmp/msgf-rust-scorefix.tgz
```

Expected: ~10–11 MB upload completes, no errors.

- [ ] **Step 3: Extract + build on VM**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm \
  "rm -rf /srv/data/msgf-bench/scorefix-build && \
   mkdir -p /srv/data/msgf-bench/scorefix-build && \
   cd /srv/data/msgf-bench/scorefix-build && \
   tar xzf ../msgf-rust-scorefix.tgz 2>&1 | tail -2 && \
   rustup override set 1.95.0 2>&1 | tail -1 && \
   cargo build --release --manifest-path rust/Cargo.toml --bin msgf-rust 2>&1 | tail -3"
```

Expected: `Finished release profile [optimized] target(s) in <N>s`.

- [ ] **Step 4: Run the 3-dataset bench**

```bash
cat > /tmp/run_scorefix_bench.sh <<'EOF'
#!/bin/bash
set -uo pipefail
cd /srv/data/msgf-bench
RUST=/srv/data/msgf-bench/scorefix-build/rust/target/release/msgf-rust
OUT=/srv/data/msgf-bench/bench-scorefix-results
mkdir -p "$OUT"

P1819=/srv/data/msgf-bench/data
ASTRAL=/srv/data/msgf-bench/astral-data
TMT=/srv/data/msgf-bench/tmt-data
TMT_MODS=$TMT/mods-numeric.txt

run() {
  local TAG=$1; local MZML=$2; local FASTA=$3; shift 3
  local PIN=$OUT/$TAG-rust.pin LOG=$OUT/$TAG-rust.log
  rm -f "$PIN" "$LOG"
  /usr/bin/time -v "$RUST" --spectrum "$MZML" --database "$FASTA" --output-pin "$PIN" \
    --ntt 2 --max-missed-cleavages 2 --min-peaks 10 --min-length 6 --max-length 40 \
    --charge-min 2 --charge-max 4 --top-n 1 --threads 8 --ms-level 2 "$@" > "$LOG" 2>&1
  local W=$(grep "Elapsed (wall clock)" "$LOG" | awk -F': ' '{print $NF}')
  local T=$([ -f "$PIN" ] && awk -F"\t" 'NR>1 && $2==1 {c++} END {print c+0}' "$PIN" || echo 0)
  local D=$([ -f "$PIN" ] && awk -F"\t" 'NR>1 && $2==-1 {c++} END {print c+0}' "$PIN" || echo 0)
  echo "[$TAG] wall=$W targets=$T decoys=$D"
}

run pxd001819 "$P1819/UPS1_5000amol_R1.mzML" "$P1819/PXD001819_uniprot_yeast_ups.fasta" \
  --precursor-tol-ppm 5 --isotope-error-min=0 --isotope-error-max=1
run astral "$ASTRAL/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML" "$ASTRAL/ProteoBenchFASTA_MixedSpecies_HYE.fasta" \
  --precursor-tol-ppm 10 --isotope-error-min=-1 --isotope-error-max=2 \
  --fragmentation 3 --instrument 3 --protocol 0
run tmt "$TMT/a05058.mzML" "$TMT/PXD007683_UP000005640_UP000002311_reviewed.fasta" \
  --precursor-tol-ppm 20 --isotope-error-min=-1 --isotope-error-max=2 \
  --fragmentation 1 --instrument 1 --protocol 4 --mod "$TMT_MODS"
echo "=== ALL DONE ==="
EOF
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm \
    "cat > /srv/data/msgf-bench/run_scorefix_bench.sh && \
     chmod +x /srv/data/msgf-bench/run_scorefix_bench.sh" \
    < /tmp/run_scorefix_bench.sh
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm \
    "nohup bash /srv/data/msgf-bench/run_scorefix_bench.sh \
     > /srv/data/msgf-bench/bench-scorefix-output.log 2>&1 &"
echo "launched"
```

Expected: bench runs for ~25 minutes (PXD001819 ~1m, Astral ~10m, TMT ~3m). Monitor via tail -f the log.

- [ ] **Step 5: Run Percolator on the three pins**

After bench completes:

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm \
  'mkdir -p /srv/data/msgf-bench/percolator-scorefix && \
   bash /srv/data/msgf-bench/run_percolator_docker.sh \
       /srv/data/msgf-bench/bench-scorefix-results/pxd001819-rust.pin \
       /srv/data/msgf-bench/percolator-scorefix pxd_rust 2>&1 | tail -1 && \
   bash /srv/data/msgf-bench/run_percolator_docker.sh \
       /srv/data/msgf-bench/bench-scorefix-results/astral-rust.pin \
       /srv/data/msgf-bench/percolator-scorefix astral_rust 2>&1 | tail -1 && \
   bash /srv/data/msgf-bench/run_percolator_docker.sh \
       /srv/data/msgf-bench/bench-scorefix-results/tmt-rust.pin \
       /srv/data/msgf-bench/percolator-scorefix tmt_rust 2>&1 | tail -1'
```

Expected output format: `[<tag>] targets_total=<N> targets_1pct=<N> targets_5pct=<N> decoys_total=<N>`.

- [ ] **Step 6: Verify the acceptance gates**

| Dataset | Gate | Java baseline | Pass criterion |
|---|---|---:|---|
| PXD001819 | targets_1pct ≥ 14,800 | 14,989 | ≥ 14,800 (~99% of Java) |
| Astral | targets_1pct ≥ 33,000 | 35,818 | ≥ 33,000 (~92% of Java) |
| TMT | targets_1pct ≥ 10,500 | 10,194 | ≥ 10,500 (no regression from current 10,548) |

If ANY gate fails:

- **PXD001819 < 14,800:** the surgical fix didn't fully resolve the bug. Likely Astral fails too. Re-open Task 5 with the per-split instrumentation escalation path.
- **Astral < 33,000 but PXD001819 ≥ 14,800:** dataset-specific second cause. Document the residual gap as a follow-up iteration (out of scope for this plan, per the spec's scope boundary).
- **TMT < 10,500:** the fix interacts with the candidate-gen fixes from the prior iteration. Revert this branch and escalate to wider design.

- [ ] **Step 7: Record results in the diagnosis doc**

Append to `docs/parity-analysis/notes/2026-05-13-score-fix-diagnosis.md`:

```markdown

## Validation results (VM bench)

Bench run: 2026-05-13
Percolator @ 1% FDR:
  PXD001819: <number> (Java 14,989, gate ≥14,800) — PASS/FAIL
  Astral:    <number> (Java 35,818, gate ≥33,000) — PASS/FAIL
  TMT:       <number> (Java 10,194, gate ≥10,500) — PASS/FAIL

Wall times (Rust):
  PXD001819: <m:ss>
  Astral:    <m:ss>
  TMT:       <m:ss>
```

```bash
git add docs/parity-analysis/notes/2026-05-13-score-fix-diagnosis.md
git commit -m "diag: VM Percolator validation results for score_psm fix"
```

---

## Task 10: Cleanup + merge into rust-implement

**Files:**
- Modify: `docs/parity-analysis/reports/2026-05-13-score-psm-undercount-finding.md` — add `## Resolution` section

- [ ] **Step 1: Update the forensic doc with a resolution section**

Open `docs/parity-analysis/reports/2026-05-13-score-psm-undercount-finding.md` and append:

```markdown

## Resolution

**Date fixed:** 2026-05-13
**Bad commit:** `<bad-sha>` ("<commit subject>")
**Fix commit:** `<fix-sha>` ("fix(scoring): score_psm under-scoring regression on PXD001819/Astral")
**Regression-test commit:** `<test-sha>` ("test(scoring): regression guard for ...")

### Root cause (one paragraph)

[Copy the diagnosis paragraph from
docs/parity-analysis/notes/2026-05-13-score-fix-diagnosis.md.]

### Post-fix metrics (VM bench, Percolator @ 1% FDR)

| Dataset | Before fix | After fix | Java baseline |
|---|---:|---:|---:|
| PXD001819 | 11,623 | <N> | 14,989 |
| Astral    | 24,828 | <N> | 35,818 |
| TMT       | 10,548 | <N> | 10,194 |

### Regression guard

`rust/crates/scoring/tests/score_psm_pxd001819_parity.rs` — four
PXD001819 scans asserted against Java's published RawScore values
(±10 tolerance). Runs in `cargo test --workspace` (integration test).
```

- [ ] **Step 2: Verify no eprintln / commented-out debugging left**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed-score-fix
git diff rust-implement..HEAD -- 'rust/crates/scoring/**/*.rs' | grep -E "eprintln|// TODO|// FIXME"
```

Expected: empty output. If anything matches, those are leftover debug — clean them in a follow-up commit before merging.

- [ ] **Step 3: Run all CI gates one more time**

```bash
cd rust
cargo build --workspace --release 2>&1 | tail -1
cargo test --workspace --lib 2>&1 | grep "^test result" | awk '{c+=$4; f+=$6} END {print "passed: "c" failed: "f}'
cargo test -p scoring --test score_psm_pxd001819_parity 2>&1 | tail -3
cargo test -p search --test gf_java_parity --release 2>&1 | tail -3
```

Expected:
- Build clean.
- Lib tests: `passed: 318 failed: 0`.
- score_psm_pxd001819_parity: `4 passed`.
- gf_java_parity: `6 passed`.

- [ ] **Step 4: Commit the resolution doc update**

```bash
git add docs/parity-analysis/reports/2026-05-13-score-psm-undercount-finding.md
git commit -m "docs: close the score_psm under-scoring forensic finding"
```

- [ ] **Step 5: Merge into rust-implement**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed
git merge --no-ff fix/score-psm-undercount -m "Merge fix/score-psm-undercount: restore PXD001819 + Astral Percolator parity

Six-phase bisect-driven investigation identified <bad-commit-sha>
('<title>') as the source of a 3x RawScore undercount on PXD001819
and Astral. Surgical fix in <file>:<lines>. Permanent regression guard
at rust/crates/scoring/tests/score_psm_pxd001819_parity.rs (4 PXD001819
scans).

Percolator @ 1% FDR before/after:
  PXD001819: 11,623 -> <N> (Java 14,989)
  Astral:    24,828 -> <N> (Java 35,818)
  TMT:       10,548 -> <N> (Java 10,194, no regression)

318 lib + 4 new integration tests + gf_java_parity all green.

Closes docs/parity-analysis/reports/2026-05-13-score-psm-undercount-finding.md."
```

Expected: merge succeeds (likely fast-forward or clean three-way).

- [ ] **Step 6: Confirm the merge**

```bash
git log --oneline -3 --first-parent
cd rust && cargo test --workspace --lib 2>&1 | grep "^test result" | awk '{c+=$4} END {print "Total: "c}'
```

Expected: merge commit on top of rust-implement; total lib test count = 318 (regression tests are in integration, not lib).

---

## Self-Review

**Spec coverage check:**

| Spec section | Plan task(s) |
|---|---|
| Phase 1 — Bisect infrastructure | Task 2 |
| Phase 2 — Bisect execution | Tasks 3 + 4 |
| Phase 3 — Diagnose bad commit | Task 5 |
| Phase 4 — Surgical fix | Task 7 |
| Phase 5 — Tests + validation | Tasks 6 + 8 + 9 |
| Phase 6 — Cleanup + merge | Task 10 |
| Worktree on `astral-speed-score-fix` / `fix/score-psm-undercount` | Task 1 |
| Bisect oracle: scan=28787 RawScore, single-threaded, `--max-spectra 30000` | Task 2 |
| Bisect range: known-good (May 10) → HEAD | Tasks 3 + 4 |
| Sister-scan regression tests (28825, 33606, 32395) | Task 8 |
| VM Percolator gates (PXD001819 ≥14,800, Astral ≥33,000, TMT ≥10,500) | Task 9 |
| Forensic doc `## Resolution` section | Task 10 |
| Out-of-scope items (other column bugs) | Not in plan — matches spec scope boundary |

All spec sections have at least one task. No gaps.

**Placeholder scan:**
- Task 7 ("Apply the surgical fix") references the diagnosis from Task 5 — this is an INTRINSIC bisect-driven plan property, not a planning failure. The plan provides three concrete change patterns covering the top hypotheses; the implementer picks one based on the diagnosis. Acceptance is bound by the regression test passing.
- Task 10's `<bad-sha>`, `<fix-sha>`, etc. are template fills extracted from earlier tasks' deterministic outputs.
- No "TBD"/"TODO"/"write tests for the above" in the plan body.

**Type consistency:**
- The peptide constructor `Peptide::new(residues, pre, post)` and `AminoAcid::standard(byte)` signatures match across Tasks 6 and 8.
- `ScoredSpectrum::new(spec, scorer, charge)` and `score_psm(scored_spec, peptide, scorer, charge, frag_tol)` signatures match across tasks.
- Param loading via `Param::load_from_file` returns a `Param`; `RankScorer::new(&param)` consumes a `&Param`. Consistent.

All checks pass. No revisions needed.
