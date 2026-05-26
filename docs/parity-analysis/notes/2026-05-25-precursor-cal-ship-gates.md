# precursorCal ship gates and follow-ups (2026-05-25)

Branch: `feat/precursor-cal`. Benchmark harness: `benchmark/vm/run_bench_calauto_3ds.sh`.

## G1 gate (blocks `--precursor-cal auto` default)

Rust @1% FDR (Percolator target PSMs, **q-value column 3**) must be within **±1%**
of Java with matching calibration mode (`--precursor-cal auto` / `-precursorCal auto`).

| Dataset | Java @1% | Rust v6 fair gate (2026-05-25, SpecE fixes) | Δ | G1 |
|---------|----------|-----------------------------------------------|---|-----|
| LFQ (PXD001819) | 15,088 | 14,755 (cal skipped, 193/200 SpecE) | −2.2% | fail |
| Astral | 36,271 | 36,715 | +1.2% | fail |
| TMT | 10,212 | 9,605 | −5.9% | fail |

Prior v5 and v6 are **identical** — SpecE boundary fixes (`java_match_score`,
score≥max→0, sink-retry removal, threshold int subtract) did not move Percolator @1%.

Prior v4 **smoke** TMT (+7.2%, 10,952) used a pre-truncation binary and tighter
window (20→6.94 ppm). The v5 fair gate includes `java_rank_score_int` (Java `(int)`
score lookup) and measured tightening 20→8.14 ppm — SpecE/Percolator ranking shifts
are expected until the SpecE tail follow-up lands.

Re-run after each PR change: `benchmark/vm/run_bench_calauto_3ds.sh` on the bench VM.

## Calibrator status — done for this PR

The MassCalibrator pre-pass is **validated** on Astral and TMT:

- Residue-mass bug fixed (`Peptide::residue_mass()` vs Java `(mz−H)×z−H₂O`)
- Learned shifts match Java (~±0.03 ppm Astral, ~±0.05 ppm TMT)
- Post-cal tolerance tightening implemented (`robust_sigma_ppm` → narrowed ppm window)

LFQ: Rust finds **193/200** confident SpecE PSMs in the pre-pass (Java learns
+0.101 ppm from 200). See optional fix below — not a calibrator logic bug.

**Do not extend calibrator work** to close G1; remaining gaps are downstream.

## Bench validity lesson (v3 invalid TMT)

The ad-hoc `bench-calauto-rust-v3` nohup script omitted explicit routing flags.
TMT Rust ran with auto-detected `CID_LowRes_Tryp.param` instead of
`CID_HighRes_Tryp.param`, inflating the apparent +9.2% gap. The committed VM
script requires `--fragmentation CID --instrument high-res --protocol TMT`.

## Follow-up: SpecE tail / Percolator feature parity (same class as Astral GF gap)

**Owner:** scoring / PIN / Percolator — **not** MassCalibrator.

Evidence from fair TMT v5 gate (Java vs Rust, both after cal + tightening, correct
`CID_HighRes_Tryp.param`):

| Metric | Java | Rust v5 |
|--------|------|---------|
| PIN target rows | 28,201 | 29,765 |
| @1% FDR | 10,212 | 9,605 (−5.9%) |
| Cal shift / tighten | −0.754 ppm, 20→6.67 ppm | −0.681 ppm, 20→8.14 ppm |

Rust emits more PIN rows but Percolator rescored **fewer** targets at 1% FDR than
Java — lnSpecE / ion-feature distribution drift (same failure mode as historical
Astral GF / SpecE tail work, `DOCS.md` §8d, `gf_java_parity` smoke gate). Prior
v4 smoke (+7.2%) used a pre-`java_rank_score_int` binary; score truncation is
Java-correct but shifts Percolator ranking until the tail follow-up lands.

Suggested investigation (defer to follow-up PR):

See also [`2026-05-25-spece-tail-exploration.md`](2026-05-25-spece-tail-exploration.md)
for local PXD001819 harness results (72.8% `spec_e_swap_only` flips).

1. Agreement-bucket lnSpecE / RawScore / ion_existence histograms (reuse
   `benchmark/parity/analyze_rust_java_pin_diff.py`)
2. GF score-int truncation vs rounding at threshold bins (partial fix in this PR
   for LFQ calibrator edge cases)
3. Percolator weight stability when lnSpecE tail is corrected

## Optional (lower priority): LFQ 193→200 SpecE

Hypothesis: Rust used `rank_score.round() as i32` for GF `spectral_probability`
lookup while Java uses `(int) score` truncation. Near the SpecE < 1e−6 filter,
rounding can flip a handful of calibrator PSMs.

**Change in this PR:** `java_rank_score_int()` (`as i32`) in
`compute_spec_e_values_for_spectrum` for `min_score` and per-PSM lookup.

**v5 smoke result:** still **193/200** SpecE passes (`305 failed SpecE` out of 498
sampled PSMs). Truncation alone does not close the LFQ calibrator gap; treat as
deferred alongside the broader SpecE tail follow-up.

## Ship recommendation

| Mode | Ship? |
|------|-------|
| `--precursor-cal off` | Yes — baseline parity unchanged (Rust CLI default until G1) |
| `--precursor-cal auto` | Opt-in only — **no** until G1 passes |

Merge the port with docs + bench harness. CLI defaults to `off`; switch default to
`auto` (Java parity) only after G1 green.
