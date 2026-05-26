# SpecE tail exploration (2026-05-25)

Follow-up to `2026-05-25-precursor-cal-ship-gates.md`. **Not a ship item** — investigation
notes for the GF / Percolator parity PR.

## What `lnSpecEValue` actually is

| Layer | Field | Meaning |
|-------|-------|---------|
| GF lookup | `psm.spec_e_value` | Raw tail `P(X ≥ score_int)` from merged bin distributions |
| PIN | `lnSpecEValue` | `ln(spec_e_value)` — **not** multiplied by `num_distinct` |
| PIN | `lnEValue` | `ln(spec_e × num_distinct)` — separate known divergence (§8d) |

GF input score is **`rank_score`** (node + cleavage + edge), cast with Java `(int)` truncation
via `java_rank_score_int()` in `match_engine.rs`. PIN **`RawScore`** is node + cleavage only.

## Code path (hot spots) — updated 2026-05-25

```
compute_spec_e_values_for_spectrum
  min_score = min(java_match_score(psm))     ← round(pin)+edge, Java getScore()
  per bin: GeneratingFunction::with_score_threshold(min_score)
           on SinkUnreachable → skip bin (Java accept early-return; no unthresholded retry)
  group.accept → linear sum of per-bin ScoreDist (group.rs)
  spec_e = group.spectral_probability(java_match_score(psm))
           score >= max_score → 0 (Java empty tail; score_dist.rs)
```

`setup_score_threshold` now uses plain `int` subtraction like Java (not `saturating_sub` at 0).

Historical **SinkUnreachable retry** (2026-05-20) was removed — it merged unthresholded
distributions Java never merges (~6954 bin skips on PXD001819).

## Fix iteration (2026-05-25) — applied, measured

| Change | Files | PXD001819 agreement lnSpecE mean Δ | Bad tail (lnSpecE > −10, target rows) |
|--------|-------|-------------------------------------|----------------------------------------|
| Baseline `rust.pin` | — | **+0.867** | Rust **17.0%** vs Java **4.8%** |
| + java_match_score, score≥max→0, skip sink retry | `match_engine.rs`, `score_dist.rs` | **+0.867** (unchanged) | **17.7%** (slightly worse) |
| + setup_score_threshold int subtract | `generating_function.rs` | **+0.867** (unchanged) | **17.7%** |

BSA histogram unchanged: 99.5% within 2 OOM; charge-3 deconv outliers remain (max 2.879 OOM).

**Conclusion:** Lookup/threshold boundary fixes are necessary Java parity hygiene but do **not**
close the lnSpecE tail gap. Dominant signal remains **GF score-distribution shape**
(node/edge scores → different merged `ScoreDist`), not score-int truncation or sink-retry policy alone.

Regenerated PIN: `benchmark/results/PXD001819-parity/rust.pin.spece-fix2` (58s local, fair flags).

**VM v6 gate (2026-05-25, SpecE fixes on bench):** Percolator @1% unchanged vs v5
(LFQ −2.2%, Astral +1.2%, TMT −5.9%). Log: `/srv/data/msgf-bench/bench-calauto-v6.log`.

Calibrator proxy (top-1 target, spec_e ≤ 1e−6): Java 17243/27011 vs Rust fix 16251/27220 — still ~1000 short, consistent with tail gap.

## Prior code path (superseded)

```
  min_score = min(java_rank_score_int(rank_score))     ← replaced by java_match_score
  on SinkUnreachable → retry GeneratingFunction::compute (unthresholded)  ← removed
  if score_int >= max_score → lookup at max_score - 1  ← removed
```

Tail math: `score_dist.rs:get_spectral_probability` — cumulative sum from `score` to `max`,
clamped to 1.0.

Historical **SinkUnreachable** fix (2026-05-20) removed ~5% sentinel rows on Astral; LFQ
parity fixture now has **0** sentinel DeNovoScore rows.

## Local evidence: PXD001819 parity PINs (no precursorCal)

Harness: `benchmark/parity/analyze_rust_java_pin_diff.py` +
`analyze_parity.py` on `benchmark/results/PXD001819-parity/{java,rust}.pin`.

### Agreement bucket (same scan + top-1 peptide, both target)

| Metric | Java | Rust |
|--------|------|------|
| lnSpecE p50 | −20.05 | −19.36 |
| lnSpecE p90 | −11.53 | −9.78 |
| lnSpecE > −10 (bad tail) | **1.8%** | **11.5%** |
| mean lnSpecE Δ (rust−java) | — | **+0.87** (~2.4× spec_e) |
| median lnSpecE Δ | — | **+1.00** |
| rust worse (Δ > 1) | — | 50.1% of agreement PSMs |

RawScore mean |Δ| ≈ 0.35 on agreement bucket — **much smaller** than lnSpecE.

### Ranking flips (9301 scans where top-1 peptide differs)

| Mode | Count | % |
|------|------:|---:|
| **spec_e_swap_only** | 6774 | **72.8%** |
| both_swap | 1817 | 19.5% |
| raw_swap | 710 | 7.6% |

`analyze_parity.py` conclusion: **SpecE / GF distribution**, not RawScore math.

### lnEValue (separate from SpecE tail)

Agreement-bucket mean lnEValue Δ ≈ **−3.68** (~exp(3.68) ≈ 40×) — dominated by
`num_distinct` semantics, tracked in DOCS §8d. Does not explain lnSpecE tail alone.

## TMT v5 gate interaction (VM, fair flags + precursorCal)

| Observation | Implication |
|-------------|-------------|
| Rust PIN rows 29,765 vs Java 28,201 | More candidates retained |
| @1% FDR 9,605 vs 10,212 (**−5.9%**) | Worse **Percolator ranking**, not fewer rows |
| `java_rank_score_int` vs prior `round()` | Truncation **raises** spec_e near integer boundaries → worse lnSpecE |
| Cal shift/tighten OK | Not a calibrator problem |

Prior v4 smoke (+7.2%) mixed pre-truncation binary + different tightening; not comparable.

## LFQ calibrator 193/200

`java_rank_score_int` did **not** reach 200 confident SpecE PSMs. Failures are broad
(305/498 fail SpecE filter), not a 7-PSM rounding cliff. Defer to SpecE tail work.

## Hypothesized mechanisms (priority order)

1. **GF score distribution shape** — deconvolution / ion_existence → different node scores
   → different `ScoreDist` after DP + enzyme neighboring-AA merge (`generating_function.rs:615-631`).
   Same class as BSA charge-3 `gf_java_parity` gap (deconv implementation divergence).

2. **Per-bin merge widening tail** — summing distributions across isotope/tolerance bins
   (`group.rs:add_prob_dist`) inflates `P(X ≥ s)` vs Java if Java merges differently.

3. **`min_score` threshold interaction** — queue min of truncated rank scores sets
   `with_score_threshold`; interacts with SinkUnreachable retry path. Lower truncated
   min_score → different pruned vs unpruned mix.

4. **`score_int >= max_score` guard** — returns `spectral_probability(max_score − 1)`;
   worth counting how often this fires on agreement-bucket outliers.

5. **Score-int boundary** — truncate vs round at `.0` boundaries; Java-correct but moves
   mass in the tail; insufficient alone for LFQ 193→200.

## Recommended next steps (when VM socket is up)

1. **TMT agreement-bucket replay** — same scripts as PXD001819 on
   `bench-calauto-results/tmt-{java,rust}.pin`; expect spec_e_swap dominance.

2. **GF trace on worst lnSpecE outliers** — top-20 agreement scans with |Δ lnSpecE| > 3;
   `-Dmsgfplus.gftrace=true` / `msgf-trace`; diff with `benchmark/parity/diff_gf_distribution.py`.

3. **Instrument `score_int >= max_score` + GF diagnostics** — per-run counts already
   in stderr atomics; add score-boundary counter for calibrator-near PSMs.

4. **Percolator feature dump** — compare learned weights Java vs Rust PIN on TMT when
   lnSpecE tail is shifted; confirms coupling vs independent bug.

5. **Do not** extend MassCalibrator — shift/tighten validated; G1 gaps are downstream.

## Local commands

```bash
# Agreement-bucket feature diff
python3 benchmark/parity/analyze_rust_java_pin_diff.py \
  --java benchmark/results/PXD001819-parity/java.pin \
  --rust benchmark/results/PXD001819-parity/rust.pin \
  --out-dir /tmp/pin-diff

# Ranking-mode classification
python3 benchmark/parity/analyze_parity.py \
  --java-pin benchmark/results/PXD001819-parity/java.pin \
  --rust-pin benchmark/results/PXD001819-parity/rust.pin \
  --output /tmp/parity-report.md

# Coarse GF smoke
cargo test -p search --test gf_java_parity
```
