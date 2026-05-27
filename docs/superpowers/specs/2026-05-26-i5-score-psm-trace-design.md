# Design — I5 score_psm trace investigation (research-only PR)

**Date:** 2026-05-26
**Branch:** `feat/i5-score-psm-trace` (from `origin/dev @ 42a6d54f`)
**Status:** Spec for review

## Problem

PR-V1 shipped a 10–15% wall reduction (FxHashMap on hot scoring tables). Wall is no longer the bottleneck for the +5%/dataset PSM goal — the bottleneck is now per-PSM scoring divergence between Rust and Java.

A prior diagnostic session (2026-05-20, captured in project auto-memory) ran `msgf-trace` on 5 label-flip PSMs from PXD001819 and found:

> "Rust scores the Java-favored target peptide R.NEEQSR.D at 14 (per-split breakdown) vs Java's RawScore 38. 20-24 point gap on the SAME (spectrum, peptide). Rust DOES enumerate the peptide (it's at #5 in Rust's top-10 queue), so candidate enumeration is fine — the divergence is in per-split node scoring inside score_psm. Pattern is universal across 5 label-flip samples (Java RawScore 13-38 vs Rust top-1 7-32, 6-22 point gap)."

Three hypotheses:
- **H1** — per-partition ion-type list differs (Rust's `partition_ion_logs` enumerates a different IonType set than Java's per-partition table)
- **H2** — peak rank assignment differs (Rust's `setRanksOfPeaks` (after precursor-filter) yields different ranks per peak)
- **H3** — per-rank log-probability tables differ (the `rank_dist_table[partition][ion_type][rank]` lookup returns different values)

That session ended with "Closing this requires Java instrumentation to dump ranks/ions for diff comparison — 2-3 day investigation." This is that investigation.

## Goal

Identify the dominant root cause (one of H1/H2/H3 or a compound) of the per-PSM scoring divergence. Output: written analysis with side-by-side evidence on the same 5 label-flip PSMs + a proposed fix design for the next PR.

**No production code changes** in this PR. Diagnostic-binary extensions (`msgf-trace`) and a Python diff harness are the only Rust code.

## Non-goals

- Implementing the fix (next PR)
- Any change to `crates/*/src/` other than `crates/msgf-rust/src/bin/msgf-trace.rs`
- Datasets other than PXD001819 (per the brainstorm; pattern is reportedly universal)
- Java repo changes committed to msgf-rust (instrumented Java patch lives in a separate java-legacy worktree on the bench VM)
- Rebasing on top of PR-V1 (this branch is off dev; PR-V1's perf changes are orthogonal to scoring correctness)

## Architecture — 4 components

### Component 1 — Rust trace extensions

File: `crates/msgf-rust/src/bin/msgf-trace.rs` (already 729 LOC, used for the 2026-05-20 finding).

Extend with structured JSON output for per-PSM per-ion diagnostics:

```json
{
  "scan": 21,
  "peptide": "R.NEEQSR.D",
  "charge": 2,
  "rust_top_rank_score": 14,
  "ions": [
    {
      "ion_type": "Prefix(c=1, off=0.0)",
      "theo_mz": 130.0498,
      "observed_peak_mz": 130.0501,
      "matched": true,
      "rank_assigned": 7,
      "max_rank_in_partition": 150,
      "log_prob_at_rank": -0.43,
      "score_contribution": -0.43
    },
    ...
  ],
  "partition": {
    "charge": 2,
    "parent_mass_tier": 1500.0,
    "seg_num": 0,
    "ion_types_count": 24,
    "ion_types": ["Prefix(c=1, off=0)", "Suffix(c=1, off=0)", ...]
  }
}
```

Output file: `--trace-json <PATH>`. Existing human-readable stderr trace stays; the JSON is additive.

Implementation: capture the per-ion data inside the existing per-split-breakdown loop; serialize with `serde_json` (already in the workspace).

### Component 2 — Java instrumentation (out-of-repo)

On the bench VM (`pride-linux-vm`):

1. Verify JDK 17 + Maven installed (`java -version; mvn -version`)
2. Clone java-legacy into a new dir: `git clone <local> /srv/data/msgf-bench/java-legacy-trace && git checkout 65120118`
3. Add `System.err.println` traces in:
   - `src/main/java/edu/ucsd/msjava/msdbsearch/DBScanScorer.java::score(...)` — log per-ion score contribution + ion type + rank
   - `src/main/java/edu/ucsd/msjava/msutil/NewScoredSpectrum.java::setRanksOfPeaks()` — log final rank assignment per peak
   - `src/main/java/edu/ucsd/msjava/msscorer/NewRankScorer.java::errorScore(...)` and the rank-lookup method — log per-rank table value
4. Each `eprintln` outputs a structured line: `TRACE\t<scan>\t<peptide>\t<field>=<value>`
5. `mvn package -DskipTests` → `target/MSGFPlus-trace.jar`
6. Run on the same 5 label-flip scans, redirect stderr to JSON-ish log

The Java patch + build artifacts live in `/srv/data/msgf-bench/java-legacy-trace/` ONLY. The instrumented JAR is NOT committed to msgf-rust. The analysis doc cites the patch's commit SHA on the java-legacy clone for reproducibility.

### Component 3 — Python diff harness

File: `benchmark/ci/diff_score_psm_traces.py` (the `benchmark/ci/` dir is the existing carve-out for committed bench tooling).

Behavior:
- Inputs: Rust trace JSON (one JSON object per scan) + Java trace log (TRACE lines, parsed into a JSON-equivalent dict)
- For each (scan, peptide) pair, align records by (ion_type_key, theoretical_mz) within a small tolerance
- Output: stdout table per (scan, peptide), columns: `IonType | Theo_mz | Rust rank | Java rank | Rust log-prob | Java log-prob | Rust contrib | Java contrib | DIVERGE?`
- Summary footer: total Rust score, total Java score, divergence count by category (rank mismatch, log-prob mismatch, ion-type-list mismatch)

Uses only stdlib (`json`, `argparse`, `collections`). No new deps.

### Component 4 — Analysis doc

File: `docs/parity-analysis/notes/2026-05-26-score-psm-trace-findings.md` — needs `.gitignore` allowlist entry alongside the existing `2026-05-25-precursor-cal-ship-gates.md`-style allowlist.

Contents:
1. Methodology (which scans, which Java commit, which Rust HEAD)
2. Five side-by-side example PSMs (diff-harness output per PSM)
3. Aggregated divergence counts by category (H1/H2/H3)
4. Code-level root cause: Rust file:line + Java file:line for the divergent path; one paragraph explaining the divergence
5. **Proposed fix design** (no code; high-level):
   - What code path to change
   - What direction (e.g., "Rust's setRanksOfPeaks needs to apply the same tie-break rule as Java")
   - Expected PSM-count impact, rough order of magnitude
   - Risk class per the n=9 audit pattern (additive vs. modifying-existing-distribution)

### Verification / success criteria

- 5+ PSMs traced with full side-by-side data
- Function-level localization: "Rust's `X::y` at file:line produces value A where Java's `Z.w` at file:line produces value B; root cause is C"
- Proposed fix design exists with the above structure
- Trace artifacts (Rust JSON + Java log + diff outputs) committed to `docs/parity-analysis/notes/score-psm-trace-artifacts/` (allowlist-relevant), small enough to commit (5 PSMs × ~kB each = tens of kB)

If after 3 days the investigation has not produced a single function-level localization but HAS produced data: ship the data + a "pending" finding doc and pause for human triage.

## Out-of-scope safety net

- **No production code change.** The `msgf-trace` binary is diagnostic — extending its JSON output cannot affect production `msgf-rust` behavior. CI bit-identical regression gate still passes trivially.
- **No Java production change.** Instrumented JAR is local-to-bench-VM; production benches still use the canonical `MSGFPlus.jar`.

## Risks & mitigations

| Risk | Mitigation |
|---|---|
| Bench VM lacks JDK 17 / Maven | Check first; install via conda or `dnf install java-17-openjdk-devel maven` |
| `java-legacy @ 65120118` doesn't build cleanly on VM | Bisect to a nearby buildable commit; document the SHA used |
| 5 PSMs produce 5 different "dominant" hypotheses | Doc reports each independently; next PR addresses them in priority order |
| Instrumented JAR's PSM counts diverge from canonical (the trace itself broke things) | Add an integrity check: run instrumented JAR vs canonical on a 100-spectrum subset; PSM counts should match within rayon-noise ±5 |
| Trace data explodes in volume (5 PSMs × dozens of ions × multiple ranks) | Cap output: matched ions only; rank list ≤ partition max_rank; per-PSM JSON ≤ 10 kB |
| Python harness misaligns Rust ↔ Java ions due to mod-name differences | Align by (theoretical_mz, ion_kind) with mz tolerance ≤ 0.001 Da; emit warnings for unmatched on either side |
| Investigation reveals divergence is in MULTIPLE places, no single root cause | OK — doc reports the full picture; fix PR can address them sequentially or pick the highest-impact first |

## Sequencing (single PR, ~3 commits)

```
feat/i5-score-psm-trace (off origin/dev @ 42a6d54f)
  ↓
Commit 1: extend msgf-trace with --trace-json output + per-ion structured fields
  ↓
Commit 2: add benchmark/ci/diff_score_psm_traces.py harness
  ↓
[out-of-repo, bench VM] Java instrumentation; build; run on 5 PSMs
  ↓
Commit 3: trace artifacts + analysis doc; gitignore allowlist entry
  ↓
PR open with the analysis doc as the PR description summary
```

## Time estimate

2-3 working days:
- Day 1 morning: extend `msgf-trace` with JSON output (commit 1)
- Day 1 afternoon: write diff harness (commit 2); verify bench VM Java toolchain
- Day 2 morning: instrument Java on VM, build, run on 5 PSMs
- Day 2 afternoon: run Rust traces; diff; preliminary findings
- Day 3 morning: write analysis doc (commit 3)
- Day 3 afternoon: iterate if needed; spec self-review; push + open PR

## Open questions

None — all design points resolved in brainstorming.

## Related documents

- Project memory: 2026-05-20 score_psm divergence finding (local-only at `docs/parity-analysis/notes/2026-05-20-score-psm-divergence.md` on a prior worktree, not in repo)
- `docs/parity-analysis/reports/2026-05-13-score-psm-undercount-finding.md` — earlier under-scoring investigation (different bug, since fixed)
- PR-V1 (`feat/quality-perf-id-rate`, in review at PR #36) — speed PR; orthogonal to this scoring-correctness work
- `docs/parity-analysis/notes/2026-05-25-spece-tail-exploration.md` — SpecE-tail context; the per-PSM scoring divergence is upstream of the lnSpecE distribution drift documented there
