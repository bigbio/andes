# R-1 retention-layer empirical test — design spec

_2026-05-18. Branch `rust-implement` at HEAD `588a630`._

## Problem

`msgf-rust` produces 75,457 target PSMs on the Astral no-mods bench; Java MS-GF+ produces
89,479. The 14,022-PSM raw-target gap exists **before** Percolator even sees the data.
A 2026-05-18 code review identified that Rust's `TopNQueue::push`
([`psm.rs:163-173`](../../../rust/crates/search/src/psm.rs#L163-L173)) uses
strict-greater eviction, dropping tied PSMs that Java explicitly keeps at three
separate retention points:
- [`DBScanner.java:540`](../../../src/main/java/edu/ucsd/msjava/msdbsearch/DBScanner.java#L540)
  — per-SpecKey raw-score retention keeps ties at capacity
- [`DBScanner.java:719-733`](../../../src/main/java/edu/ucsd/msjava/msdbsearch/DBScanner.java#L719-L733)
  — pre-merge dedup keyed by `pepSeq + score`
- [`DBScanner.java:745`](../../../src/main/java/edu/ucsd/msjava/msdbsearch/DBScanner.java#L745)
  — per-spectrum SpecE merge keeps ties at capacity

The audit + iter3-5 cycle never tested the retention layer. Per the empirical investigation
([`docs/parity-analysis/notes/2026-05-18-piecewise-fixes-dont-work.md`](../../parity-analysis/notes/2026-05-18-piecewise-fixes-dont-work.md)),
this is the most-likely upstream prerequisite for any scoring/feature fix to land
cleanly.

## Goal

Empirically test the hypothesis: **R-1 (tied-PSM eviction) is the dominant cause of the
14K-PSM raw-target gap on Astral.**

This is a focused experiment, not a parity push. The smallest possible code change that
isolates the variable.

## Success criteria

After the R-1 change, the Astral bench (no mods, threads=8, same fixture and CLI flags
Java's bench uses) must produce:

- **Raw target count ≥ 86,674** (closes ≥80% of the 75,457 → 89,479 gap)
- **AND** target-to-decoy ratio does not worsen by ≥10% (currently 1.63; must stay ≥1.47)
- **AND** wall time does not exceed 3× the b1d45bb baseline (~8:36 on Astral threads=8, so cap at ~26 min)
- **AND** existing tests still pass: `gf_java_parity`, `score_psm_vs_gf_dp_edge_parity`,
  `score_psm_pxd001819_parity`, `match_engine_java_parity`

Outcome decision tree (raw target count after R-1, before Percolator):

| Range | Interpretation | Next |
|---|---|---|
| ≥ 86,674 | **Hypothesis validated.** R-1 was the dominant cause. | Commit R-1. Plan R-2 next iteration. |
| 78,000 – 86,673 | R-1 real but partial. R-2 likely coupled. | Commit R-1. Next iteration: R-1+R-2 together. |
| 75,458 – 77,999 | R-1 effect small. Gap is elsewhere. | Commit R-1 if T/D ratio + tests OK; pivot investigation. |
| ≤ 75,457 | R-1 had no effect or lost PSMs. | Revert. Hypothesis falsified. |
| Any count with T/D ratio < 1.47 | R-1 retains noise more than signal. | Revert. Hypothesis falsified. |

## Architecture

One file changed. One function modified.

### `rust/crates/search/src/psm.rs`

Current `TopNQueue::push` (line 163-173) uses strict-greater eviction. Replace with a
match on `Ordering` that keeps ties at capacity (matching Java `DBScanner.java:540`
behavior):

```rust
pub fn push(&mut self, m: PsmMatch) {
    if self.heap.len() < self.capacity as usize {
        self.heap.push(Reverse(m));
    } else if let Some(Reverse(top)) = self.heap.peek() {
        match m.cmp(top) {
            std::cmp::Ordering::Greater => {
                self.heap.pop();
                self.heap.push(Reverse(m));
            }
            std::cmp::Ordering::Equal => {
                // R-1 fix (2026-05-18): Java keeps ties at capacity
                // (DBScanner.java:540, :745). Queue grows past capacity
                // when ties exist, matching Java's unbounded-tied-PSM
                // behavior. See docs/parity-analysis/notes/
                // 2026-05-18-piecewise-fixes-dont-work.md and
                // docs/superpowers/specs/2026-05-18-r1-tie-retention-test-design.md
                self.heap.push(Reverse(m));
            }
            std::cmp::Ordering::Less => { /* drop */ }
        }
    }
}
```

Update the docstring on `TopNQueue::push` to reflect that the queue may exceed
`capacity` when tied PSMs exist (this is intentional — the capacity becomes a *minimum*
top-N, not a hard cap).

### What this does NOT change

- `PsmMatch::cmp` ordering (preserved as-is)
- Queue iteration / drain order (preserved — `into_sorted_vec` etc work the same)
- Downstream consumers (`compute_psm_features`, `compute_spec_e_values_for_spectrum`,
  `pin.rs`) all iterate `queue.iter_psms()` — they will see more PSMs but no API change
- The capacity field on `TopNQueue` (kept; just no longer strictly enforced)

### Architectural note for future iterations (out of scope)

The Rust port has ONE `TopNQueue` per spectrum. Java has ONE `PriorityQueue<DatabaseMatch>`
per `SpecKey` (i.e., per `(spectrum, charge)`), and merges across charges per-spectrum
later. This architectural difference means Rust's tie semantics, even after R-1, won't be
perfectly faithful to Java's per-charge tie keeping. For the raw-target-count hypothesis
test, this is acceptable. The full per-SpecKey refactor is "Approach C" and deferred to
a future iteration.

## Testing

Add one targeted unit test in `psm.rs` `#[cfg(test)] mod tests`:

```rust
#[test]
fn topn_queue_keeps_ties_at_capacity() {
    // Three tied PSMs into capacity=1 queue. Java keeps all three;
    // Rust pre-R-1 kept only one. Asserts R-1 fix is live.
    let mut q = TopNQueue::new(1);
    // ... construct three PsmMatch with identical (spec_e_value, score) ...
    q.push(psm_a);
    q.push(psm_b);
    q.push(psm_c);
    assert_eq!(q.len(), 3, "all three tied PSMs should be retained at capacity=1");
}
```

Use the same `PsmMatch` constructor pattern as the existing tests in `psm.rs`. The
test must not depend on any spectrum/scoring path.

## Bench protocol

1. **Build iter6 = b1d45bb + R-1.** Sync to pride-linux-vm via tarball:
   ```
   tar --exclude='rust/target' -czf /tmp/msgf-rust-iter6.tgz rust
   scp -o ControlPath=/tmp/msgfplus-bench.sock /tmp/msgf-rust-iter6.tgz \
       root@pride-linux-vm.ebi.ac.uk:/srv/data/msgf-bench/
   ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm.ebi.ac.uk '
       rm -rf /srv/data/msgf-bench/track-iter6-build &&
       mkdir -p /srv/data/msgf-bench/track-iter6-build &&
       cd /srv/data/msgf-bench/track-iter6-build &&
       tar -xzf /srv/data/msgf-bench/msgf-rust-iter6.tgz &&
       cp -r /srv/data/msgf-bench/track-iter3-build/src ./ &&
       cd rust &&
       printf "[toolchain]\nchannel = \"stable\"\ncomponents = [\"rustfmt\", \"clippy\"]\n" > rust-toolchain.toml &&
       cargo build --release -p msgf-rust
   '
   ```

2. **Run no-mods Astral bench** (matches Java's `-precursorCal off -tda 1 -t 10ppm` config).
   **Important:** run via `nohup` because the SSH ControlMaster has dropped 6 times during
   prior sessions; bench must survive transient SSH disconnections:
   ```
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
   ```
   Then poll for completion: `until ! pgrep -f astral-rust-r1.pin >/dev/null; do sleep 60; done`

3. **Measure:**
   - Raw target count: `awk -F"\t" 'NR>1 && $2==1 {c++} END {print c+0}' astral-rust-r1.pin`
   - Raw decoy count: `awk -F"\t" 'NR>1 && $2==-1 {c++} END {print c+0}' astral-rust-r1.pin`
   - T/D ratio: target / decoy
   - Wall time: from `time` output

4. **Compare against:**
   - b1d45bb baseline: 75,457 targets / 46,208 decoys / T/D = 1.63
   - Java baseline: 89,479 targets / 46,792 decoys / T/D = 1.91
   - Success bar: 86,674 targets, T/D ≥ 1.47

5. **(Optional) Run Percolator** if raw count moves significantly. The hypothesis test
   does NOT require Percolator — raw counts are sufficient. But Percolator @ 1% FDR is
   informative.

## Rollback plan

If any of the rollback criteria fire (T/D worsens, wall blows up, test regresses, or
raw count moves the wrong way):

1. Revert the change in `psm.rs` to its pre-R-1 state.
2. Add a follow-up note to
   [`docs/parity-analysis/notes/2026-05-18-piecewise-fixes-dont-work.md`](../../parity-analysis/notes/2026-05-18-piecewise-fixes-dont-work.md)
   recording the empirical falsification of the R-1 hypothesis. The audit catalog stays;
   the next iteration's strategy changes.
3. No production impact — b1d45bb baseline remains the production HEAD.

## Out of scope (NOT this iteration)

- R-2 (per-charge GF computation) — deferred. If R-1 succeeds, this is next.
- R-3 (`minDeNovoScore` PIN filter)
- R-4 (`lnEValue` denominator length-indexing)
- F-1 (`matched_ion_ratio` denominator)
- All 17 audit-catalog scoring + feature fixes (C-1..C-6, B-1, B-2)
- Strengthening `match_engine_java_parity.rs` (acknowledged-too-weak, but blocking nothing
  in THIS experiment since we're measuring raw .pin counts directly, not via the test)
- Per-SpecKey architectural refactor (Approach C)
- TMT + PXD001819 benches — only Astral matters for the R-1 hypothesis
- Mods configuration changes — explicitly no-mods to match Java's reference bench

## Risks

| Risk | Mitigation |
|---|---|
| Queue grows unboundedly on pathological spectra (e.g., 100K tied PSMs) | Watch wall time + bench output. If queue size goes crazy, revert + add a soft cap (Approach B). |
| `compute_psm_features` / `compute_spec_e_values_for_spectrum` perf scales linearly with queue size; getting slower | Sanity-check wall time at each step. Acceptable up to 3× baseline. |
| Tied PSMs introduce different protein/peptide variants, breaking downstream parsers | Existing tests cover `match_engine_java_parity` (peptide identity for top-1). If it regresses, investigate. |
| R-1 changes interact with `into_sorted_vec` ordering (which is consumed by `pin.rs::write_spectrum_rows` for rank assignment) | Tied PSMs get equal `cmp` → stable sort keeps insertion order. Should not affect rank, but spot-check `lnDeltaSpecEValue` doesn't change semantics. |

## Estimated effort

- Implementation: 30 min (one diff + unit test)
- Sync + build on VM: 15 min
- Astral bench: 8-13 min
- Analysis: 15 min
- **Total: 1-1.5 hours wall time**

## Deliverables

1. One commit on `rust-implement` branch: `fix(search): TopNQueue keeps tied PSMs at capacity (R-1)`
2. New unit test in `psm.rs`
3. Bench results documented in the iteration note (extend
   [`docs/parity-analysis/notes/2026-05-18-piecewise-fixes-dont-work.md`](../../parity-analysis/notes/2026-05-18-piecewise-fixes-dont-work.md)
   or new file under `docs/parity-analysis/notes/2026-05-18-r1-bench-results.md`)
4. Decision per the outcome table — commit + plan next iteration, or revert + document.
