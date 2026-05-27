# I5 score_psm trace investigation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Identify the dominant root cause of the Rust↔Java per-PSM scoring divergence (Rust ~14 vs Java ~38 RawScore on the same spectrum+peptide) for 5 known label-flip PSMs on PXD001819, by capturing structured per-ion traces on both sides and diffing them. Output: written analysis + proposed fix design for the next PR.

**Architecture:** Three small artifacts: (a) extend `msgf-trace` with `--trace-json` for per-PSM per-ion JSON output, (b) instrument java-legacy on the bench VM with `System.err.println` traces, (c) Python diff harness that aligns the two outputs and emits side-by-side rows. No production code changes; CI bit-identical regression gate passes trivially.

**Tech Stack:** Rust 2024 edition pinned to 1.87.0; JSON output written manually via `write!` (no new serde dep); Java instrumentation against `java-legacy @ 65120118` built with Maven on bench VM (`pride-linux-vm`); Python 3 stdlib for the diff harness.

**Spec:** `docs/superpowers/specs/2026-05-26-i5-score-psm-trace-design.md`

---

## File map

**Created in this PR:**
- `crates/msgf-rust/src/bin/msgf-trace.rs` — extended (existing 729 LOC; add `--trace-json` flag + per-ion JSON output writer)
- `benchmark/ci/diff_score_psm_traces.py` — Python diff harness
- `docs/parity-analysis/notes/2026-05-26-score-psm-trace-findings.md` — analysis doc (allowlisted in `.gitignore`)
- `docs/parity-analysis/notes/score-psm-trace-artifacts/` — directory with the 5-PSM Rust JSON traces + Java trace logs + diff outputs (small, ~tens of kB)
- `.gitignore` — allowlist entries for the new note + artifacts dir

**Out-of-repo (bench VM only):**
- `/srv/data/msgf-bench/java-legacy-trace/` — fresh clone of `java-legacy` branch with instrumentation patch
- `/srv/data/msgf-bench/java-legacy-trace/target/MSGFPlus-trace.jar` — built instrumented JAR

---

## The 5 label-flip PSMs (from 2026-05-20 finding)

Per project memory, the 2026-05-20 investigation found 5 scans on PXD001819 where Rust and Java disagree on top-1 peptide. The flagship example is **scan 21** where Rust scores Java-favored peptide `R.NEEQSR.D` at 14 vs Java's RawScore 38.

The exact 5 scan IDs are documented in the 2026-05-20 doc (local-only at the time, may need re-derivation):

```bash
# To re-derive on bench VM if the original list is unavailable:
ssh root@pride-linux-vm 'cd /srv/data/msgf-bench/bench-pr-v1-s1b-results && \
  python3 /srv/data/msgf-bench/diff_top1.py \
    pxd001819-java.pin pxd001819-rust-off.pin | head -20'
```

A small re-derivation script (5 scans of the largest |Java RawScore − Rust top-1 RawScore| where both agree on the peptide candidate enumeration) can be added if the 2026-05-20 list is missing. For this plan, assume the scans are available; document the actual scan IDs in the analysis doc.

---

## Pre-flight (run before Task 1)

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed
git branch --show-current
# Expect: feat/i5-score-psm-trace

git log origin/dev..HEAD --oneline | wc -l
# Expect: 1 (the spec commit f943aa7e)

git status --short
# Expect: empty (clean tree)

cargo build --release -p msgf-rust --bin msgf-trace 2>&1 | tail -3
# Expect: Finished release profile

cargo test --release --workspace -- \
  --skip charge_missing_spectrum_uses_per_charge_scored_spec \
  --skip spectrum_without_charge_tries_charge_range \
  --skip known_peptide_appears_in_top_n \
  --skip read_bsa_canno_text_format \
  --skip read_tryp_pig_bov_revcat_csarr_cnlcp \
  --skip tryp_pig_bov_revcat_full_set_loads \
  --skip match_spectra_output_invariant_across_thread_counts 2>&1 | grep -E "^test result" | grep -vE "0 passed.*0 failed.*0 ignored" | tail -5
# Expect: all 0 failed.
```

If pre-flight fails, STOP and investigate.

---

## Task 1: Extend `msgf-trace` with `--trace-json` output

**Goal:** Add a flag that, when set, writes per-PSM per-ion structured JSON to a file alongside the existing human-readable stderr trace.

**Files:**
- Modify: `crates/msgf-rust/src/bin/msgf-trace.rs`

- [ ] **Step 1: Add the CLI flag**

Open `crates/msgf-rust/src/bin/msgf-trace.rs`. Find the `struct Cli` definition (around line 30). After the existing `--java-top1` field, add:

```rust
    /// Output structured per-PSM per-ion JSON to this path (additive; the
    /// existing human-readable stderr trace is unaffected).
    #[arg(long)]
    trace_json: Option<PathBuf>,
```

- [ ] **Step 2: Locate the per-split breakdown loop**

In the same file, find where the per-split / per-ion breakdown is computed for the top-1 PSM (and the optional `--java-top1` peptide). Look for the loop that calls `directional_node_score_inner` or `partition_ion_logs` or `nearest_peak_rank` — that's the data source for the JSON.

```bash
grep -nE "partition_ion_logs|nearest_peak_rank|directional_node_score|partition_for" crates/msgf-rust/src/bin/msgf-trace.rs | head -20
```

Identify the line ranges where the per-ion data is produced.

- [ ] **Step 3: Add a JSON-writer module to msgf-trace.rs**

Near the top of the file (after imports, before the `Cli` struct), add:

```rust
// ─── Per-PSM JSON trace output (additive; no new deps) ─────────────────────
//
// Hand-written JSON via `write!` macros: small output (~5-10 KB per PSM),
// no serde dependency, and the diff harness parses on the Python side
// where stdlib json is sufficient.

use std::io::Write as _;

struct TraceJson<W: std::io::Write> {
    out: W,
    first_psm: bool,
}

impl<W: std::io::Write> TraceJson<W> {
    fn new(mut out: W) -> std::io::Result<Self> {
        out.write_all(b"[\n")?;
        Ok(Self { out, first_psm: true })
    }

    fn begin_psm(
        &mut self,
        scan: i32,
        peptide: &str,
        charge: u8,
        rust_rank_score: i32,
    ) -> std::io::Result<()> {
        if !self.first_psm {
            self.out.write_all(b",\n")?;
        }
        self.first_psm = false;
        write!(
            self.out,
            "  {{\n    \"scan\": {},\n    \"peptide\": \"{}\",\n    \"charge\": {},\n    \"rust_rank_score\": {},\n    \"ions\": [",
            scan, escape_json(peptide), charge, rust_rank_score
        )
    }

    fn end_psm(&mut self) -> std::io::Result<()> {
        self.out.write_all(b"\n    ]\n  }")
    }

    fn ion(
        &mut self,
        first_ion: bool,
        ion_type: &str,
        theo_mz: f64,
        rank_assigned: Option<u32>,
        max_rank: u32,
        log_prob: f32,
        contribution: f32,
    ) -> std::io::Result<()> {
        if !first_ion {
            self.out.write_all(b",")?;
        }
        let rank_str = rank_assigned
            .map(|r| r.to_string())
            .unwrap_or_else(|| "null".to_string());
        write!(
            self.out,
            "\n      {{\"ion_type\": \"{}\", \"theo_mz\": {:.6}, \"rank\": {}, \"max_rank\": {}, \"log_prob\": {:.6}, \"contribution\": {:.6}}}",
            escape_json(ion_type), theo_mz, rank_str, max_rank, log_prob, contribution
        )
    }

    fn finish(mut self) -> std::io::Result<()> {
        self.out.write_all(b"\n]\n")
    }
}

fn escape_json(s: &str) -> String {
    s.replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
        .replace('\t', "\\t")
}
```

- [ ] **Step 4: Wire the JSON writer into the per-split breakdown loop**

In `fn main()`, after parsing the CLI, before the per-split-breakdown loop, add:

```rust
    let mut trace_json: Option<TraceJson<_>> = match cli.trace_json {
        Some(ref path) => {
            let file = File::create(path).map_err(|e| {
                eprintln!("Failed to create --trace-json output {}: {}", path.display(), e);
                e
            })?;
            Some(TraceJson::new(std::io::BufWriter::new(file))?)
        }
        None => None,
    };
```

Then INSIDE the per-PSM per-split-breakdown loop where the human-readable stderr is already being emitted, add parallel JSON emissions:

```rust
    // Inside the loop where you iterate over `(rust top-1, optional java_top1)`:
    if let Some(ref mut tj) = trace_json {
        tj.begin_psm(cli.scan, &peptide_label, charge, rust_rank_score as i32)?;
        let mut first_ion = true;
        for seg in 0..num_segs {
            let partition = param.partition_for(charge, parent_mass, seg);
            let ion_logs = scorer.partition_ion_logs(&partition);
            for (ion, logs) in ion_logs {
                let theo_mz = ion.mz(nominal_mass);  // adjust to whatever drives the inner loop
                let tol_da = param.mme.as_da(theo_mz);
                let rank = ss.nearest_peak_rank(theo_mz, tol_da);
                let max_rank = scorer.max_rank();
                let (log_prob, contribution) = match rank {
                    Some(r) => {
                        let idx = (r.min(max_rank).max(1) as usize) - 1;
                        let lp = if idx < logs.len() { logs[idx] } else { 0.0 };
                        (lp, lp)
                    }
                    None => {
                        // No peak: missed-ion slot is logs[max_rank as usize] if present.
                        let lp = logs.get(max_rank as usize).copied().unwrap_or(0.0);
                        (lp, lp)
                    }
                };
                tj.ion(
                    first_ion,
                    &format!("{:?}", ion),
                    theo_mz,
                    rank,
                    max_rank,
                    log_prob,
                    contribution,
                )?;
                first_ion = false;
            }
        }
        tj.end_psm()?;
    }
```

The exact details of where this slots into the existing 729-line file depend on the current structure. **Step 4a:** before writing the loop body, READ the existing `main()` function and figure out:
- Where is `peptide_label` available (the peptide being scored)?
- Where is `parent_mass` computed?
- Where is `num_segs` (`param.num_segments`)?
- Where is `nominal_mass` derived per inner iteration?

Use those bindings in your insertion. If the existing code uses different field names, adapt.

- [ ] **Step 5: Close the JSON document at end of main**

At the bottom of `main()`, just before the final `ExitCode::SUCCESS` return:

```rust
    if let Some(tj) = trace_json {
        tj.finish()?;
    }
```

- [ ] **Step 6: Build + smoke test**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed
cargo build --release -p msgf-rust --bin msgf-trace 2>&1 | tail -3
# Expect: Finished

./target/release/msgf-trace --help 2>&1 | grep -A 1 "trace-json"
# Expect: --trace-json line with description
```

- [ ] **Step 7: Functional smoke test (local fixture)**

```bash
# Use a small in-tree fixture so we don't depend on bench VM data.
./target/release/msgf-trace \
  --spectrum test-fixtures/test.mgf \
  --database test-fixtures/BSA.fasta \
  --param resources/ionstat/HCD_QExactive_Tryp.param \
  --scan 1 \
  --trace-json /tmp/smoke-trace.json 2>&1 | tail -5

# Validate JSON parses:
python3 -c "import json; j=json.load(open('/tmp/smoke-trace.json')); print(f'PSMs: {len(j)}, first ions: {len(j[0][\"ions\"])}' if j else 'empty')"
# Expect: at least one PSM with at least one ion record, JSON parses cleanly.
```

- [ ] **Step 8: Workspace tests + clippy**

```bash
cargo test --release --workspace -- \
  --skip charge_missing_spectrum_uses_per_charge_scored_spec \
  --skip spectrum_without_charge_tries_charge_range \
  --skip known_peptide_appears_in_top_n \
  --skip read_bsa_canno_text_format \
  --skip read_tryp_pig_bov_revcat_csarr_cnlcp \
  --skip tryp_pig_bov_revcat_full_set_loads \
  --skip match_spectra_output_invariant_across_thread_counts 2>&1 | grep -E "^test result" | grep -vE "0 passed.*0 failed.*0 ignored" | tail -5

cargo clippy --workspace --all-targets 2>&1 | tail -3
```

Both must pass. `msgf-trace` is a diagnostic binary so any new code there doesn't affect production correctness.

- [ ] **Step 9: Commit**

```bash
git add crates/msgf-rust/src/bin/msgf-trace.rs
git commit -m "$(cat <<'COMMIT_EOF'
feat(msgf-trace): per-PSM per-ion JSON output via --trace-json

Adds a structured output mode to the diagnostic trace binary so its
per-split breakdown can be diffed against Java's instrumentation
output. JSON is written by hand (no new serde dep) since the volume
is small (~5-10 KB per PSM). The existing human-readable stderr
output is unaffected.

No production code change; msgf-trace is a separate binary from
msgf-rust.
COMMIT_EOF
)"
```

---

## Task 2: Python diff harness

**Goal:** Take a Rust trace JSON file + a Java trace log file, produce a side-by-side per-ion comparison.

**Files:**
- Create: `benchmark/ci/diff_score_psm_traces.py`

- [ ] **Step 1: Create the script**

```bash
mkdir -p benchmark/ci
```

Create `benchmark/ci/diff_score_psm_traces.py` with:

```python
#!/usr/bin/env python3
"""
Diff per-PSM per-ion trace outputs from Rust (msgf-trace --trace-json) and
Java (instrumented java-legacy stderr). For each (scan, peptide) PSM, align
records by (ion_kind, theoretical mz tolerance 1e-3 Da) and emit a side-by-side
table.

Usage:
    diff_score_psm_traces.py --rust rust-trace.json --java java-trace.log \\
        [--mz-tol 1e-3] [--scan SCAN] [--peptide PEP]

Outputs to stdout. Exit code 0 = success.

Rust JSON shape (per PSM):
    {
      "scan": int,
      "peptide": str,
      "charge": int,
      "rust_rank_score": int,
      "ions": [
        {"ion_type": str, "theo_mz": float, "rank": int|null,
         "max_rank": int, "log_prob": float, "contribution": float},
        ...
      ]
    }

Java log shape (one line per ion):
    TRACE\\tscan=<int>\\tpeptide=<str>\\tion=<str>\\ttheo_mz=<float>\\trank=<int>\\tlog_prob=<float>\\tcontribution=<float>
"""

import argparse
import collections
import json
import sys


def parse_java_log(path: str) -> dict:
    """Returns {(scan, peptide): [{ion fields}, ...]}."""
    out = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.startswith("TRACE\t"):
                continue
            fields = {}
            for part in line.split("\t")[1:]:
                if "=" not in part:
                    continue
                k, v = part.split("=", 1)
                fields[k] = v
            try:
                scan = int(fields["scan"])
                peptide = fields["peptide"]
                ion = {
                    "ion_type": fields.get("ion", "?"),
                    "theo_mz": float(fields.get("theo_mz", "nan")),
                    "rank": int(fields["rank"]) if fields.get("rank", "") not in ("", "-1", "null") else None,
                    "log_prob": float(fields.get("log_prob", "nan")),
                    "contribution": float(fields.get("contribution", "nan")),
                }
            except (KeyError, ValueError) as e:
                print(f"WARN: skipping malformed Java TRACE line: {line[:80]}... ({e})", file=sys.stderr)
                continue
            out[(scan, peptide)].append(ion)
    return out


def parse_rust_json(path: str) -> dict:
    """Returns {(scan, peptide): [{ion fields}, ...]}."""
    out = {}
    with open(path) as fh:
        data = json.load(fh)
    for psm in data:
        key = (psm["scan"], psm["peptide"])
        out[key] = psm["ions"]
    return out


def normalize_ion_kind(s: str) -> str:
    """Map both Rust and Java ion-type representations to a normalized key.

    Rust format example:    `Prefix { charge: 1, offset_bits: 0 }`
    Java format example:    `b/1+ off=0.0`  (or whatever Java's TRACE emits)
    Normalize to:           `b/1+0.0` or `y/1+0.0` or `Noise`.
    """
    s = s.strip()
    if "Noise" in s:
        return "Noise"
    # Rust: `Prefix { charge: <c>, offset_bits: <bits-as-int> }`
    if s.startswith("Prefix"):
        # extract charge and offset_bits, reconstruct as `b/<c>+<off_f32>`
        import re
        m = re.search(r"charge:\s*(\d+).*offset_bits:\s*(\d+)", s)
        if m:
            charge = int(m.group(1))
            off_bits = int(m.group(2))
            # Decode f32::from_bits(u32) — use struct to avoid float imports
            import struct
            off = struct.unpack(">f", struct.pack(">I", off_bits))[0]
            return f"b/{charge}+{off:.5f}"
    if s.startswith("Suffix"):
        import re, struct
        m = re.search(r"charge:\s*(\d+).*offset_bits:\s*(\d+)", s)
        if m:
            charge = int(m.group(1))
            off_bits = int(m.group(2))
            off = struct.unpack(">f", struct.pack(">I", off_bits))[0]
            return f"y/{charge}+{off:.5f}"
    # Java format (placeholder; tighten when actual Java TRACE format is known)
    return s


def align_and_diff(rust_ions: list, java_ions: list, mz_tol: float = 1e-3):
    """Yields rows: (key, rust, java, diverge_flags) per matched/unmatched ion."""
    java_by_key = collections.defaultdict(list)
    for ion in java_ions:
        key = (normalize_ion_kind(ion["ion_type"]), round(ion["theo_mz"] / mz_tol))
        java_by_key[key].append(ion)

    matched_java = set()
    for rust_ion in rust_ions:
        rust_key = (
            normalize_ion_kind(rust_ion["ion_type"]),
            round(rust_ion["theo_mz"] / mz_tol),
        )
        candidates = java_by_key.get(rust_key, [])
        java_ion = candidates.pop(0) if candidates else None
        if java_ion is not None:
            matched_java.add(id(java_ion))
        flags = []
        if java_ion is None:
            flags.append("RUST_ONLY")
        else:
            if rust_ion["rank"] != java_ion["rank"]:
                flags.append("RANK_DIFF")
            if abs(rust_ion["log_prob"] - java_ion["log_prob"]) > 1e-4:
                flags.append("LOGPROB_DIFF")
            if abs(rust_ion["contribution"] - java_ion["contribution"]) > 1e-4:
                flags.append("CONTRIB_DIFF")
        yield (rust_key, rust_ion, java_ion, flags)

    # Any remaining Java ions not matched in Rust:
    for ion in java_ions:
        if id(ion) in matched_java:
            continue
        key = (normalize_ion_kind(ion["ion_type"]), round(ion["theo_mz"] / mz_tol))
        yield (key, None, ion, ["JAVA_ONLY"])


def format_row(rust_key, rust_ion, java_ion, flags):
    def fmt(v, w):
        if v is None:
            return "-" * w
        if isinstance(v, float):
            return f"{v:>{w}.4f}"
        return f"{str(v):>{w}}"
    return "  ".join([
        fmt(rust_key[0], 22),
        fmt((rust_ion or java_ion)["theo_mz"], 10),
        fmt(rust_ion["rank"] if rust_ion else None, 5),
        fmt(java_ion["rank"] if java_ion else None, 5),
        fmt(rust_ion["log_prob"] if rust_ion else None, 9),
        fmt(java_ion["log_prob"] if java_ion else None, 9),
        fmt(rust_ion["contribution"] if rust_ion else None, 9),
        fmt(java_ion["contribution"] if java_ion else None, 9),
        ",".join(flags) if flags else "",
    ])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rust", required=True, help="Rust trace JSON from msgf-trace --trace-json")
    ap.add_argument("--java", required=True, help="Java instrumented trace log (TRACE lines)")
    ap.add_argument("--mz-tol", type=float, default=1e-3, help="m/z alignment tolerance (Da)")
    ap.add_argument("--scan", type=int, default=None, help="Restrict to one scan")
    ap.add_argument("--peptide", default=None, help="Restrict to one peptide")
    args = ap.parse_args()

    rust = parse_rust_json(args.rust)
    java = parse_java_log(args.java)

    all_keys = sorted(set(rust.keys()) | set(java.keys()))
    for key in all_keys:
        scan, pep = key
        if args.scan is not None and scan != args.scan:
            continue
        if args.peptide is not None and pep != args.peptide:
            continue
        print(f"\n=== scan={scan} peptide={pep} ===")
        rust_ions = rust.get(key, [])
        java_ions = java.get(key, [])
        if not rust_ions and not java_ions:
            print("  (no data on either side)")
            continue
        print("  ion_type                theo_mz     R_rk   J_rk    R_logP   J_logP   R_ctrb    J_ctrb   flags")
        rust_total = 0.0
        java_total = 0.0
        category_counts = collections.Counter()
        for row in align_and_diff(rust_ions, java_ions, args.mz_tol):
            print("  " + format_row(*row))
            if row[1] is not None:
                rust_total += row[1]["contribution"]
            if row[2] is not None:
                java_total += row[2]["contribution"]
            for f in row[3]:
                category_counts[f] += 1
        print(f"  TOTAL contribution: rust={rust_total:.4f}  java={java_total:.4f}  delta={rust_total - java_total:+.4f}")
        if category_counts:
            print(f"  DIVERGENCES: {dict(category_counts)}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Make executable + smoke test**

```bash
chmod +x benchmark/ci/diff_score_psm_traces.py

# Synthetic test: create tiny rust + java trace inputs and run
cat > /tmp/rust-smoke.json <<'EOF'
[
  {"scan": 1, "peptide": "K.PEPTIDE.D", "charge": 2, "rust_rank_score": 10,
   "ions": [
     {"ion_type": "Prefix { charge: 1, offset_bits: 0 }", "theo_mz": 100.05, "rank": 5, "max_rank": 150, "log_prob": -0.4, "contribution": -0.4},
     {"ion_type": "Suffix { charge: 1, offset_bits: 0 }", "theo_mz": 200.10, "rank": null, "max_rank": 150, "log_prob": -2.1, "contribution": -2.1}
   ]}
]
EOF

cat > /tmp/java-smoke.log <<'EOF'
TRACE	scan=1	peptide=K.PEPTIDE.D	ion=b/1+0.00000	theo_mz=100.05	rank=4	log_prob=-0.35	contribution=-0.35
TRACE	scan=1	peptide=K.PEPTIDE.D	ion=y/1+0.00000	theo_mz=200.10	rank=-1	log_prob=-2.05	contribution=-2.05
EOF

python3 benchmark/ci/diff_score_psm_traces.py --rust /tmp/rust-smoke.json --java /tmp/java-smoke.log
# Expect: a table showing rust=5 vs java=4 (RANK_DIFF) + LOGPROB_DIFF + CONTRIB_DIFF
# Total delta: rust=-2.5, java=-2.4, delta=-0.1.
```

- [ ] **Step 3: Commit**

```bash
git add benchmark/ci/diff_score_psm_traces.py
git commit -m "$(cat <<'COMMIT_EOF'
feat(diff-harness): Python diff for Rust vs Java per-PSM ion traces

Aligns msgf-trace JSON output against java-legacy instrumented TRACE
lines by (ion_kind, theo_mz). Emits side-by-side per-ion rows with
RANK_DIFF / LOGPROB_DIFF / CONTRIB_DIFF flags + per-PSM totals.
stdlib-only; runs on any Python 3 install.
COMMIT_EOF
)"
```

---

## Task 3: Bench VM Java instrumentation

**Goal:** Build an instrumented `MSGFPlus-trace.jar` on the bench VM and capture the 5-PSM trace log.

**Files:** none in this repo (all changes live on the bench VM under `/srv/data/msgf-bench/java-legacy-trace/`).

- [ ] **Step 1: Verify VM Java toolchain + reactivate VM socket if needed**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm 'java -version 2>&1 | head -3; mvn -version 2>&1 | head -3'
```

Expected: Java 17 (or 11) and Maven 3.x. If missing, install:

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm 'dnf install -y java-17-openjdk-devel maven 2>&1 | tail -5'
```

- [ ] **Step 2: Clone java-legacy on VM**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm 'cd /srv/data/msgf-bench && \
  rm -rf java-legacy-trace && \
  git clone https://github.com/bigbio/msgf-rust.git java-legacy-trace && \
  cd java-legacy-trace && \
  git checkout 65120118 && \
  git log -1 --format="%h %s"'
```

If the commit `65120118` isn't reachable (e.g., the java-legacy branch was removed), bisect from the most recent commit on the `java-legacy` or `java-legacy-original` branch.

- [ ] **Step 3: Apply instrumentation patch on the VM**

```bash
# Edit DBScanScorer.java to add TRACE prints in the score path.
# Pattern: in the score-summing inner loop, before adding ion contribution to total:
#   System.err.println("TRACE\tscan=" + scanNum + "\tpeptide=" + peptideStr + "\tion=" + ionType + "\ttheo_mz=" + theoMz + "\trank=" + rank + "\tlog_prob=" + logProb + "\tcontribution=" + contribution);
```

Use `sed` or paste a patch via stdin from the controller side. The exact insertion line depends on java-legacy's code structure. Reference patch shape (the actual lines to add, given by the agent on demand):

```java
// In DBScanScorer.java, score(...) method, inside the per-ion loop:
double contribution = /* existing per-ion score */;
System.err.println(
    "TRACE\tscan=" + scanNum +
    "\tpeptide=" + peptideStr +
    "\tion=" + ionType.toString() +
    "\ttheo_mz=" + theoMz +
    "\trank=" + rank +
    "\tlog_prob=" + logProb +
    "\tcontribution=" + contribution
);
totalScore += contribution;
```

Apply via heredoc/scp; commit on the VM-side clone (not pushed):

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm 'cd /srv/data/msgf-bench/java-legacy-trace && \
  # patch applied via Edit on VM-side files; commit:
  git add -A && \
  git commit -m "diag: TRACE per-ion prints for I5 investigation" && \
  git log -1 --format="%h %s"'
```

Note the SHA — cite it in the analysis doc.

- [ ] **Step 4: Build instrumented JAR**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm 'cd /srv/data/msgf-bench/java-legacy-trace && \
  mvn package -DskipTests 2>&1 | tail -10'
# Expect: BUILD SUCCESS; target/MSGFPlus-*.jar exists.
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm 'ls -la /srv/data/msgf-bench/java-legacy-trace/target/*.jar | head'
```

If build fails, capture the error, downgrade to a nearby buildable commit on java-legacy, document the actual SHA used.

- [ ] **Step 5: Identify the 5 label-flip scans**

If the 2026-05-20 doc is unavailable, derive from current PR-V1 bench data:

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm 'python3 <<EOF
# Read java + rust PINs from PR-V1 bench output; find scans where top-1 peptide diverges
import csv, collections
java_top = {}
rust_top = {}
with open("/srv/data/msgf-bench/bench-pr-v1-s1b-results/pxd001819-java.pin") as fh:
    rows = csv.DictReader(fh, delimiter="\t")
    for r in rows:
        scan = r.get("ScanNr") or r.get("scan")
        if scan and scan not in java_top:
            java_top[scan] = (r.get("Peptide"), r.get("RawScore"))
with open("/srv/data/msgf-bench/bench-pr-v1-s1b-results/pxd001819-rust-off.pin") as fh:
    rows = csv.DictReader(fh, delimiter="\t")
    for r in rows:
        scan = r.get("ScanNr") or r.get("scan")
        if scan and scan not in rust_top:
            rust_top[scan] = (r.get("Peptide"), r.get("RawScore"))
divergent = []
for s, (jp, jr) in java_top.items():
    rp, rr = rust_top.get(s, (None, None))
    if rp is None or jp == rp:
        continue
    try:
        jr_i = int(float(jr))
        rr_i = int(float(rr))
        divergent.append((jr_i - rr_i, s, jp, rp))
    except (TypeError, ValueError):
        continue
divergent.sort(reverse=True)
for d, s, jp, rp in divergent[:10]:
    print(f"scan={s}\tjava_RawScore-rust_top1_RawScore={d}\tjava_peptide={jp}\trust_peptide={rp}")
EOF'
```

This produces a ranked list of label-flip scans with the largest scoring gap. Pick the top 5; record them in the analysis doc.

- [ ] **Step 6: Run instrumented Java on the 5 scans**

For each of the 5 scans, run the instrumented JAR against PXD001819 spectra, redirect stderr (where TRACE lines go) to a per-scan log:

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm 'mkdir -p /srv/data/msgf-bench/i5-trace-out && \
  cd /srv/data/msgf-bench && \
  for SCAN in <5-scan-ids-here>; do \
    java -Xmx8192m -jar java-legacy-trace/target/MSGFPlus-*.jar \
      -s data/UPS1_5000amol_R1.mzML \
      -d data/PXD001819_uniprot_yeast_ups.fasta \
      -mod mods.txt \
      -o /tmp/java-trace-$SCAN.mzid \
      -tda 1 -t 5ppm -ti 0,1 -m 0 -inst 0 -e 1 -protocol 0 -ntt 2 \
      -minLength 6 -maxLength 40 -minNumPeaks 10 \
      -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1 \
      -msLevel 2 -thread 8 \
      2>/srv/data/msgf-bench/i5-trace-out/java-trace-scan-$SCAN.log; \
  done'
```

Note: the instrumented JAR will produce TRACE lines for ALL scans it processes, not just the 5 we care about. The Python diff harness will filter by `--scan`. Alternative: add a scan filter inside the Java instrumentation (e.g., `if (scanNum != TARGET_SCAN) return;`) to keep log volume manageable.

If log size is unmanageable (>1 GB), add a runtime filter in Java code (a `Set<Integer>` of target scans, only print TRACE when contained).

- [ ] **Step 7: Run msgf-rust trace on the same 5 scans**

```bash
# Make sure msgf-rust binary is up to date with Task 1's commit
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm 'cd /srv/data/msgf-bench/pr-v1-s1b-build && /root/.cargo/bin/cargo build --release --bin msgf-trace 2>&1 | tail -3'

# Or: scp updated source from local, rebuild
# (skip if VM build is fresh)

# Run msgf-trace on each scan with --trace-json
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm 'cd /srv/data/msgf-bench && \
  for SCAN in <5-scan-ids-here>; do \
    pr-v1-s1b-build/target/release/msgf-trace \
      --spectrum data/UPS1_5000amol_R1.mzML \
      --database data/PXD001819_uniprot_yeast_ups.fasta \
      --param resources/ionstat/HCD_QExactive_Tryp.param \
      --scan $SCAN \
      --java-top1 "<java-top1-peptide-for-this-scan>" \
      --trace-json /srv/data/msgf-bench/i5-trace-out/rust-trace-scan-$SCAN.json \
      > /srv/data/msgf-bench/i5-trace-out/rust-trace-scan-$SCAN.txt 2>&1; \
  done'
```

- [ ] **Step 8: Run the diff harness for each scan**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm 'cd /srv/data/msgf-bench && \
  for SCAN in <5-scan-ids-here>; do \
    echo "=== scan $SCAN diff ==="; \
    python3 /srv/data/msgf-bench/diff_score_psm_traces.py \
      --rust /srv/data/msgf-bench/i5-trace-out/rust-trace-scan-$SCAN.json \
      --java /srv/data/msgf-bench/i5-trace-out/java-trace-scan-$SCAN.log \
      --scan $SCAN > /srv/data/msgf-bench/i5-trace-out/diff-scan-$SCAN.txt; \
    tail -5 /srv/data/msgf-bench/i5-trace-out/diff-scan-$SCAN.txt; \
  done'
```

(Make sure to scp `benchmark/ci/diff_score_psm_traces.py` to the VM as `/srv/data/msgf-bench/diff_score_psm_traces.py` first, or run from a clone of this branch on the VM.)

- [ ] **Step 9: Pull artifacts to local**

```bash
mkdir -p docs/parity-analysis/notes/score-psm-trace-artifacts
scp -o ControlPath=/tmp/msgfplus-bench.sock \
  'root@pride-linux-vm:/srv/data/msgf-bench/i5-trace-out/*' \
  docs/parity-analysis/notes/score-psm-trace-artifacts/
ls -la docs/parity-analysis/notes/score-psm-trace-artifacts/
# Expect: ~15 files (5 rust json + 5 java log + 5 diff txt). Total ~50-500 KB.
```

Note: the Java log files may be large. If any exceed 1 MB, filter them down to TRACE lines for the 5 target scans only:

```bash
for f in docs/parity-analysis/notes/score-psm-trace-artifacts/java-trace-scan-*.log; do
  scan=$(basename "$f" .log | sed 's/java-trace-scan-//')
  grep "TRACE.*scan=${scan}\b" "$f" > "${f}.filtered" && mv "${f}.filtered" "$f"
done
```

- [ ] **Step 10: No commit yet** (artifacts staged in Task 4 alongside the analysis doc).

---

## Task 4: Write the analysis doc + .gitignore allowlist

**Goal:** Read the diff outputs from Task 3 Step 8, identify the dominant root cause, write the analysis doc with side-by-side evidence and a proposed fix design.

**Files:**
- Create: `docs/parity-analysis/notes/2026-05-26-score-psm-trace-findings.md`
- Modify: `.gitignore` (allowlist the new note + artifacts dir)

- [ ] **Step 1: Read the 5 diff outputs**

```bash
for s in <5-scan-ids-here>; do
  echo "=== scan $s ==="
  cat docs/parity-analysis/notes/score-psm-trace-artifacts/diff-scan-${s}.txt
done
```

For each scan, identify:
- Are there RANK_DIFF flags? If yes, how many ions show rank mismatch?
- Are there LOGPROB_DIFF flags? Where do they cluster?
- Are there CONTRIB_DIFF flags driven by rank or by log-prob?
- Are there RUST_ONLY / JAVA_ONLY ions (ion-type-list mismatch)?

Tally divergence categories across all 5 scans. The category with the most ion-level divergences AND the largest score-delta contribution is the dominant root cause.

- [ ] **Step 2: Localize to code**

Once a dominant category is identified:

- **H1 dominant** (ion-type-list mismatch): inspect Rust's `crates/scoring/src/scoring/rank_scorer.rs::partition_ion_logs` vs Java's `NewRankScorer.getIonProbabilities(Partition)` or equivalent. Capture the file:line on both sides where the ion-type set is constructed.
- **H2 dominant** (rank mismatch): inspect Rust's `crates/scoring/src/scoring/scored_spectrum.rs::nearest_peak_rank` + `setRanksOfPeaks`-equivalent vs Java's `NewScoredSpectrum.setRanksOfPeaks`. Particularly check the precursor-filter handling and rank tie-break behavior.
- **H3 dominant** (log-prob mismatch): inspect Rust's `crates/scoring/src/param_model.rs::partition_for` + the rank index calculation (`r.min(max_rank).max(1) as usize - 1`) vs Java's analogous lookup.

Document the divergence with code citations.

- [ ] **Step 3: Write the analysis doc**

Create `docs/parity-analysis/notes/2026-05-26-score-psm-trace-findings.md`:

```markdown
# I5 score_psm trace investigation — findings

**Date:** 2026-05-26
**Branch:** feat/i5-score-psm-trace
**Java instrumentation:** java-legacy @ <commit-sha-on-vm-clone> (out-of-repo)
**Dataset:** PXD001819 (UPS1_5000amol_R1.mzML)

## Five label-flip PSMs traced

| Scan | Java top-1 peptide | Java RawScore | Rust top-1 peptide | Rust RawScore | Δ |
|---:|---|---:|---|---:|---:|
| <scan1> | ... | ... | ... | ... | ... |
| <scan2> | ... | ... | ... | ... | ... |
| <scan3> | ... | ... | ... | ... | ... |
| <scan4> | ... | ... | ... | ... | ... |
| <scan5> | ... | ... | ... | ... | ... |

Trace artifacts: `score-psm-trace-artifacts/{rust-trace-scan-N.json, java-trace-scan-N.log, diff-scan-N.txt}`.

## Aggregate divergence counts (5 PSMs combined)

| Category | Count | % of total divergences |
|---|---:|---:|
| RANK_DIFF | <N> | <P>% |
| LOGPROB_DIFF | <N> | <P>% |
| CONTRIB_DIFF | <N> | <P>% |
| RUST_ONLY | <N> | <P>% |
| JAVA_ONLY | <N> | <P>% |

## Dominant root cause

<one-paragraph description with code citations>

**Rust:** `crates/<path>:<line>`
**Java:** `<path>:<line>` (in java-legacy clone)

The divergence arises because <explanation>.

## Proposed fix design

**Code path to change:** <Rust file:line>
**Direction:** <what to change; do NOT include actual code>
**Expected PSM impact:** estimated +<N>% on PXD001819 (~+<NN> PSMs at 1% FDR). On Astral and TMT, likely <similar / different> based on <reasoning>.
**Risk class:** <additive / modifies-existing-distribution> per the n=9 audit pattern.
**Bench gate for the fix PR:** PXD001819 auto @1% FDR ≥ +<NN> PSMs; no regression on Astral / TMT.

## Methodology

1. Identified 5 label-flip PSMs from PR-V1 bench (largest |Java RawScore − Rust top-1 RawScore| where peptide differs).
2. Captured per-ion structured traces:
   - Rust: `msgf-trace --trace-json` (commit <task1-sha>)
   - Java: java-legacy with `System.err.println` patches in `DBScanScorer.score()` (java-legacy clone commit <vm-sha>)
3. Aligned Rust ↔ Java records by (ion_kind, theo_mz) tolerance 1e-3 Da.
4. Diff harness: `benchmark/ci/diff_score_psm_traces.py` (commit <task2-sha>).

## Out of scope (next PR)

- Implementing the fix
- Validating the fix on Astral / TMT (the bench gate is PXD001819 only, but Astral / TMT should be monitored for regressions)
```

Replace all `<...>` placeholders with actual values from your investigation.

- [ ] **Step 4: Update .gitignore allowlist**

Open `.gitignore`. Find the existing parity-analysis allowlist:

```gitignore
docs/parity-analysis/*
!docs/parity-analysis/notes/
!docs/parity-analysis/notes/2026-05-25-precursor-cal-ship-gates.md
!docs/parity-analysis/notes/2026-05-25-spece-tail-exploration.md
```

Add:

```gitignore
!docs/parity-analysis/notes/2026-05-26-score-psm-trace-findings.md
!docs/parity-analysis/notes/score-psm-trace-artifacts/
!docs/parity-analysis/notes/score-psm-trace-artifacts/*
```

- [ ] **Step 5: Confirm files are tracked**

```bash
git check-ignore docs/parity-analysis/notes/2026-05-26-score-psm-trace-findings.md && echo "STILL_IGNORED" || echo "TRACKED"
# Expect: TRACKED

git check-ignore docs/parity-analysis/notes/score-psm-trace-artifacts/diff-scan-21.txt && echo "STILL_IGNORED" || echo "TRACKED"
# Expect: TRACKED
```

(Adjust the example scan-id to one of the 5 actual scans.)

- [ ] **Step 6: Stage and commit**

```bash
# Stage allowlist + analysis doc + artifacts
git add .gitignore
git add docs/parity-analysis/notes/2026-05-26-score-psm-trace-findings.md
git add docs/parity-analysis/notes/score-psm-trace-artifacts/

git status --short
# Expect: 4 new entries (gitignore + note + artifacts dir + diff harness already-committed).

git commit -m "$(cat <<'COMMIT_EOF'
docs(i5): per-PSM trace findings + 5-PSM artifacts (PXD001819)

Identifies the dominant root cause of the Rust vs Java per-PSM scoring
divergence on PXD001819 label-flip PSMs. Methodology + artifacts +
proposed fix design (no code in this PR; fix lands separately).

Dominant cause: <H1|H2|H3> — Rust's <code path> diverges from Java's
<code path>.

Trace artifacts (Rust JSON + Java TRACE log + diff outputs for 5
PSMs) committed under docs/parity-analysis/notes/score-psm-trace-artifacts/
for reproducibility.

Out of scope: fix implementation; next PR after this.
COMMIT_EOF
)"
```

Replace the placeholder `<H1|H2|H3> — Rust's <code path> diverges from Java's <code path>` in the message with the actual finding before running the commit.

---

## Task 5: Push + open PR

- [ ] **Step 1: Final workspace check**

```bash
cargo build --release --workspace 2>&1 | tail -3
# Expect: Finished

cargo test --release --workspace -- \
  --skip charge_missing_spectrum_uses_per_charge_scored_spec \
  --skip spectrum_without_charge_tries_charge_range \
  --skip known_peptide_appears_in_top_n \
  --skip read_bsa_canno_text_format \
  --skip read_tryp_pig_bov_revcat_csarr_cnlcp \
  --skip tryp_pig_bov_revcat_full_set_loads \
  --skip match_spectra_output_invariant_across_thread_counts 2>&1 | grep -E "^test result" | grep -vE "0 passed.*0 failed.*0 ignored" | tail -5
# Expect: all 0 failed.
```

- [ ] **Step 2: Confirm commit ladder**

```bash
git log origin/dev..HEAD --oneline
# Expect:
#   <task4 sha> docs(i5): per-PSM trace findings ...
#   <task2 sha> feat(diff-harness): ...
#   <task1 sha> feat(msgf-trace): per-PSM per-ion JSON output ...
#   f943aa7e   docs(spec): I5 score_psm trace investigation design
```

- [ ] **Step 3: Push**

```bash
git push -u origin feat/i5-score-psm-trace 2>&1 | tail -3
```

- [ ] **Step 4: Open PR**

```bash
gh pr create --base dev --head feat/i5-score-psm-trace \
  --title "diag(i5): score_psm trace findings + diff harness (no production code change)" \
  --body "$(cat <<'PR_BODY'
## Summary

Research-only PR. Identifies the dominant root cause of the Rust vs
Java per-PSM scoring divergence (Rust ~14 vs Java ~38 RawScore on the
same spectrum+peptide). The actual fix is a separate PR after this.

## Finding

<one-sentence summary; e.g. "Rust's nearest_peak_rank applies a
different tie-break rule than Java's setRanksOfPeaks, causing
systematic rank inflation on PSMs with multiple peaks in the same
window">

Full analysis with side-by-side evidence on 5 label-flip PSMs from
PXD001819: `docs/parity-analysis/notes/2026-05-26-score-psm-trace-findings.md`.

## What this PR contains

- `crates/msgf-rust/src/bin/msgf-trace.rs` — extended with `--trace-json`
  for per-PSM per-ion structured output (no production code change;
  diagnostic binary)
- `benchmark/ci/diff_score_psm_traces.py` — Python diff harness
- `docs/parity-analysis/notes/2026-05-26-score-psm-trace-findings.md` — analysis
- `docs/parity-analysis/notes/score-psm-trace-artifacts/` — Rust + Java
  traces + diff outputs for 5 PSMs (reproducibility)

## What this PR does NOT contain

- The fix itself (next PR)
- Production code changes (`msgf-trace` is a separate binary)
- Java repo changes (java-legacy instrumentation lives on bench VM)
- Datasets other than PXD001819

## Verification

- [x] `cargo clippy --workspace --all-targets` clean
- [x] Workspace tests green under existing CI skip list
- [x] `precursor_cal_bit_identical` regression gate green (no
  production code change → trivially passes)
- [ ] CodeRabbit review pass
- [ ] CI matrix green

## Next PR

The proposed fix from the analysis doc, bench-gated on PXD001819
@1% FDR.
PR_BODY
)"
```

Replace the `<one-sentence summary>` placeholder with the actual finding from Task 4.

- [ ] **Step 5: Confirm PR open**

```bash
gh pr view --json number,title,state,statusCheckRollup --jq '{number, state, checks: [.statusCheckRollup[]? | {name, status}]}'
```

---

## Self-review

I checked the plan against the spec section-by-section:

**1. Spec coverage:**
- Component 1 (Rust trace extensions) → Task 1 ✓
- Component 2 (Java instrumentation, out-of-repo) → Task 3 ✓
- Component 3 (Python diff harness) → Task 2 ✓
- Component 4 (analysis doc + artifacts) → Task 4 ✓
- Verification / success criteria (5+ PSMs, function-level localization, fix design) → Task 4 ✓
- Out-of-scope safety net (no production code change) → Task 1 (msgf-trace is diagnostic) + Task 3 (Java patch out-of-repo) ✓

**2. Placeholder scan:** The plan contains `<5-scan-ids-here>` and `<one-sentence summary>` style placeholders intentionally — they are inputs the implementer fills in from the live investigation. Each is documented as such. No "TBD" or "implement later" instructions for things that should be specified upfront.

**3. Type consistency:** The JSON field names (`ion_type`, `theo_mz`, `rank`, `max_rank`, `log_prob`, `contribution`) are used identically across Task 1 (writer), Task 2 (parser), and Task 4 (analysis). The Java TRACE format (tab-separated `key=value`) is used identically in Task 2's parser and Task 3's emitter.

**Known soft spots:**
- The exact Java instrumentation patch lines depend on the actual java-legacy source structure at SHA `65120118`. Task 3 Step 3 provides the pattern; the agent fills in line-specific edits.
- The 5 scan IDs depend on either the 2026-05-20 doc (local-only) OR a re-derivation script (Task 3 Step 5). If re-derivation produces a different set, that's acceptable; document the actual scans used.
- If the diff harness reveals that NONE of H1/H2/H3 dominates and the cause is more subtle (e.g., a numeric-precision issue in a different code path), the analysis doc reports that honestly and the next PR has a wider scope.
