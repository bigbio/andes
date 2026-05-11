# E-value SearchIndex iteration 1 — PXD001819 results

## Context

Task 5 of the e-value SearchIndex plan validates that auto-populating
`SearchIndex.num_distinct_peptides_at_length` (T2-1/2/3 + T2-4) and consuming
it in the Phase 7 e-value block produces EValue column values closer to
Java's, without regressing Percolator targets @ 1% FDR.

- Branch: `rust-implement`
- Prior commits gating this iteration: f5f6884, a547c39, 95fa9bc
- Pre-evalue snapshot: `benchmark/results/PXD001819-parity/rust.pin.pre-evalue`
- Post-evalue artifact: `benchmark/results/PXD001819-parity/rust.pin`

## Rust arm wall time

`match_spectra` wall: 242.51s. End-to-end `real`: 9m17.871s with 12 worker
threads. PSM yield: 37,112 / 37,918 spectra with PSMs (1 spectrum skipped by
min_peaks). Yield is unchanged from the pre-evalue snapshot, as expected
(scoring is unaltered; only the EValue column is affected).

## Ratio histogram (Rust EValue / Java EValue, full-match PSMs)

### Post-evalue (this iteration)

```
N=27779 full-match PSMs with finite EValues
Median ratio: 0.0368
P10 / P90:    0.0045 / 0.0997
Within +/-5%: 4 (0.0%)
Min / Max:    0.0000 / 32.9899
```

### Pre-evalue baseline (for direction-of-change reference)

```
[pre-evalue baseline] N=15929
Median ratio: 0.0000
P10 / P90:    0.0000 / 0.0000
Within +/-5%: 0 (0.0%)
Min / Max:    0.0000 / 0.0001
```

Two observations from the comparison:

1. The count of full-match PSMs with finite (non-zero, positive) Rust
   EValues grew from 15,929 to 27,779 — more PSMs now have a usable
   EValue column at all.
2. The median ratio moved from effectively 0 (Rust EValues many orders of
   magnitude below Java) to 0.0368 — still ~27x lower than Java in the
   median, but moving in the correct direction by roughly four orders of
   magnitude. P90 is now within ~10x of Java.

The acceptance gate ([0.95, 1.05] median, >= 95% within +/-5%) is not met
this iteration. Residual discrepancy is consistent with one or more of the
remaining algorithmic divergences tracked in
`project_phase6_parity_root_causes.md` (notably the RawScore scale mismatch
and survival-function inputs) still feeding into the e-value path.

## Percolator @ FDR

```
[rust] targets_total=27199 targets_1pct=14850 targets_5pct=16617 decoys_total=9913 q_col=3
```

## Decision

- **Median ratio:** 0.0368 (target: [0.95, 1.05]) — FAIL
- **% within +/-5%:** 0.0% (target: >= 95%) — FAIL
- **Percolator @ 1% FDR:** 14,850 targets (baseline: 14,798; acceptance: >= 14,798) — PASS
- **Direction vs proxy:** strictly better (median ratio moved from ~0 to
  0.0368; finite-EValue PSM count grew from 15,929 to 27,779)
- **Status:** DONE_WITH_CONCERNS

The EValue ratio gate fails this iteration, but the change is strictly
better than the pre-evalue proxy and Percolator @ 1% FDR exceeds the
baseline (14,850 >= 14,798). Per the plan's residual-discrepancy clause,
ship and document. Further iterations should target the upstream parity
gaps (RawScore scale, survival-function inputs) before re-running this
acceptance gate.
