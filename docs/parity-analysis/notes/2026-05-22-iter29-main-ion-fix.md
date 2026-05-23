# iter29: main_ion_direction fix — DeNovoScore now at parity, +379 PSMs

_2026-05-22. First top-1-changing fix to gain Astral PSMs (n=10 audit pattern broken). DeNovoScore agreement-bucket median Δ collapsed from -13 to **0**._

## The fix (1-block)

In `rust/crates/scoring/src/scoring/scored_spectrum.rs:915`,
`main_ion_from_param`:

```rust
// Before (BUG):
for (ion, freqs) in table {
    if !ion.is_prefix() {          // ← filters to prefix ions ONLY
        continue;
    }
    let freq_at_rank1 = freqs.first().copied().unwrap_or(0.0);
    ...
}
```

was unconditionally selecting a prefix ion as `main_ion`, forcing
`main_ion_direction() = true` for every spectrum. For HCD/QExactive
Astral data, Java's `NewRankScorer.determineIonTypes` picks a **y-ion
(suffix)** as the most-frequent ion → direction = false → reverse loop.

iter28 single-scan trace (`docs/parity-analysis/notes/2026-05-21-iter27-pin-diff.md`)
localized this on scan 47106:

| trace | direction | edge_score sum |
|---|---|---:|
| Java | reverse (y-main) | +8 |
| Rust iter27 | forward (b-main) | -18 |

The fix mirrors Java exactly: aggregate `frag_off_table` frequencies
across all segments for the same `(charge, parent_mass)` partition, then
pick the overall highest-frequency ion regardless of type.

## Verification on scan 47106

| | RawScore | DeNovoScore | EdgeScore |
|---|---:|---:|---:|
| Java | 73 | 73 | +8 (effective) |
| Rust iter27 | 65 | 60 | -18 |
| Rust iter29 | 65 | **73** ✓ | **+8** ✓ |

DeNovoScore matches Java bit-for-bit. EdgeScore matches Java's effective
contribution bit-for-bit. RawScore stays at 65 because Rust's pin RawScore
column is node+cleavage-only by design (edge lives in the separate
iter19 `EdgeScore` PIN column).

## Astral bench

| Metric | iter27 | **iter29** | Java | gap |
|---|---:|---:|---:|---:|
| 1% FDR (PSMs) | 31,298 | **31,677** | 35,818 | 11.6% |
| 5% FDR (PSMs) | n/a | 34,713 | — | — |
| Targets total | 92,764 | 92,781 | n/a | — |
| Decoys total | 56,921 | 56,796 | n/a | — |
| T/D ratio | 1.630 | 1.634 | 2.02 | — |
| Wall (8 threads) | n/a | 7:32 | — | — |

**Δ +379 PSMs / +1.2%.** Gap to Java closed from 12.6% to 11.6% (1 point).

Cumulative since iter16 baseline (26,432): **+5,245 PSMs / +19.8%**.

## Pin-diff: DeNovoScore parity achieved

Agreement-bucket (50,450 PSMs where both engines pick the same target
peptide):

| Feature | iter27 median Δ | iter29 median Δ | iter27 mean \|Δ\| | iter29 mean \|Δ\| |
|---|---:|---:|---:|---:|
| **DeNovoScore** | -13 | **0** | 14.89 | **1.42** |
| RawScore | -2 | -2 | 8.56 | 8.56 |
| lnSpecEValue | -2.33 | +1.20 | 3.67 | 3.13 |
| NumMatchedMainIons | -1 | -1 | 1.96 | 1.96 |

DeNovoScore is now at near-bit-exact parity with Java. The first iter
to truly close the GF-max gap.

RawScore stays at -2 (unchanged) because Rust's pin RawScore column
deliberately excludes edge_score (it lives in a separate `EdgeScore` PIN
column per iter19). Java's pin RawScore includes edge by virtue of
`DBScanScorer.getScore` overriding `FastScorer.getScore`.

lnSpecEValue went from -2.33 → +1.20: Rust's spec_e_value is now
slightly tighter (more confident) than Java's on average, reflecting the
GF DP's improved per-edge scoring.

## Why this fix doesn't regress (n=10 audit pattern)

Prior top-1-changing fixes (iter3, iter17, iter18, iter23, units fix)
regressed Astral 1% FDR. This one helped. Why?

**The prior fixes shifted Rust's scoring distribution AWAY from Java's
intended behavior.** Percolator had learned weights against Rust's
current (buggy but consistent) distribution; modifying that distribution
broke the learned weights.

**This fix RESTORES Java's intended direction.** Percolator's weights
were always tuned for Java-direction scoring; Rust was just emitting
features in a flipped frame. Restoring direction lets Percolator's
existing feature weights extract signal that was previously masked by
the direction mismatch.

The audit pattern refines from "top-1-changing regresses" to "top-1-
changing that DIVERGES from Java's intended behavior regresses; top-1-
changing that RESTORES Java's intended behavior gains".

## Commits

- `994cf1a0` fix(scoring): `main_ion_from_param` picks overall most-
  frequent ion, not prefix-only (iter29)
- `c7566100` docs+trace: per-edge trace localizes the bug
- `7e4b3d50` docs: iter28 trace closes Layer 1

Trace artifacts on EBI VM at `/tmp/trace-{rust,java}-47106-*.{err,stderr}`.
