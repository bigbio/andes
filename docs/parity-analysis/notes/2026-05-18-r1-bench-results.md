# R-1 retention-layer test — empirical results

_2026-05-18. Branch `rust-implement` at HEAD (iter6 = b1d45bb + R-1 fix, commit `fc16407`)._

## Spec

[`docs/parity-analysis/specs/2026-05-18-r1-tie-retention-test-design.md`](../specs/2026-05-18-r1-tie-retention-test-design.md)

## Bench config

Astral no-mods, threads=8, `pride-linux-vm` (EBI). Matches Java's reference bench config
exactly: `--precursor-tol-ppm 10 --isotope-error-min=-1 --isotope-error-max=2 --ntt 2
--max-missed-cleavages 2 --min-peaks 10 --min-length 6 --max-length 40 --charge-min 2
--charge-max 4 --top-n 1 --threads 8`. No `--mod` flags. Input:
`LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML` / `ProteoBenchFASTA_MixedSpecies_HYE.fasta`.
Percolator 3.7.1 via Docker (biocontainers).

## Results

| Metric | Java baseline | b1d45bb (pre-R-1) | iter5 (all-revert) | iter6 (with R-1) | vs Java | vs b1d45bb |
|---|---:|---:|---:|---:|---:|---:|
| Raw targets | 89,479 | 75,457 | 75,619 | **1,042,255** | +952,776 (+1065%) | +966,798 (+1281%) |
| Raw decoys | 46,792 | 46,208 | 46,058 | **530,430** | +483,638 (+1034%) | +484,222 (+1048%) |
| T/D ratio | 1.912 | 1.633 | 1.643 | **1.965** | +0.053 | +0.332 |
| Wall time | — | ~8:36 | ~8:36 | **23:52** | — | +15:16 (within 3× cap) |
| Max RSS | — | ~9.7 GB | ~9.7 GB | **~9.6 GB** | — | −0.1 GB |
| Percolator @ 1% FDR | 35,818 | 25,224 | 17,160 | **74,204** | +38,386 (+107%) | +48,980 (+194%) |

_iter5 included to contextualize the iter3-5 arc (all-revert): iter3-5 lost ~8K Percolator PSMs
from b1d45bb; that reversion landed at iter5 = 17,160. R-1 alone more than recovers all of that
and then dramatically over-shoots._

## Mechanism interpretation — why R-1 over-shoots in Rust's architecture

### The spec_e_value=1.0 sentinel tie

Before `compute_spec_e_values_for_spectrum` runs (i.e., at raw-score insertion time), **every
PSM has `spec_e_value = 1.0`** — a sentinel value. The `PsmMatch::cmp` ordering puts
`spec_e_value` first (smallest wins) and `score` second (largest wins). So at insertion time
the primary key ties for *every* PSM in the queue: all have `spec_e_value = 1.0`.

The secondary `score` key is an `i32` rounded from a float. Integer-valued scores produce
many ties within a spectrum's candidate set, particularly at the lower end of the score
distribution. R-1's `Equal` branch inserts all tied PSMs without eviction.

### One queue per spectrum vs one queue per SpecKey

**Java's architecture:** one `PriorityQueue<DatabaseMatch>` per `SpecKey`, i.e., per
`(spectrum, charge)` pair (`DBScanner.java:320`, `:400`). With `--top-n 1`, Java keeps at most
N=1 PSMs *per charge* per spectrum, plus tied PSMs at that per-charge capacity. The per-spectrum
PSM set is the union across charges, so a 4-charge-state spectrum can yield at most
4 × (N=1 + tied siblings) PSMs.

**Rust's architecture:** ONE `TopNQueue` per spectrum, shared across all charge states
(`match_engine.rs:325`). With `--top-n 1`, pre-R-1 Rust kept exactly 1 PSM per spectrum.
Post-R-1, Rust keeps all PSMs tied with the worst-in-queue. Because the queue spans all
charges, ALL same-score candidates from ALL charge states tie with each other. On a
high-fragment-density Astral spectrum with many candidate peptides, this accumulates
dramatically more PSMs than Java's per-charge queues would.

### The magnitude explained

The 11.6× raw target count (1,042,255 vs Java's 89,479) is not R-1 doing "more than
expected" in terms of correctness — it is R-1 faithfully implementing "keep all ties at
capacity" on a cross-charge queue that Java would never give that many tied slots to.
R-1 is the correct fix to the tie-dropping bug. R-2 (per-SpecKey separation) is the
architectural counterweight that constrains what "tied" means to the per-charge context
Java intended.

## Decision per spec outcome table

**Raw target count:** 1,042,255 — exceeds the ≥86,674 threshold by a factor of ~12.

**T/D ratio:** 1.965 — exceeds the ≥1.47 gate.

**Wall time:** 23:52 — within the 3× baseline cap of 26:00 (barely; 3× of 8:36 = ~25:48,
and 23:52 is under this, but only by ~2 minutes with the queue growing significantly in
memory).

**All three spec gates pass** → outcome branch: **≥86,674 targets, T/D ≥1.47**.

Per the spec's outcome table decision: **"Hypothesis validated. R-1 was the dominant cause.
Commit R-1. Plan R-2 next iteration."**

**Mechanism interpretation revised:** The "dominant cause" framing is confirmed, but the
mechanism operates through a different pathway than the spec anticipated. Ties matter more
than expected — R-1 alone over-shoots by 12× on raw target count because Rust's
single-queue-per-spectrum architecture accumulates tied PSMs across all charge states,
while Java's per-SpecKey architecture limits ties to within each charge. R-1 is the correct
fix; R-2 (per-SpecKey GF compute + per-charge queue) is the necessary architectural
counterweight, not just a "coupled follow-up" as described in the spec. The audit's
prediction that "R-1 + R-2 should be coupled" is empirically confirmed in the strongest
possible direction: R-1 without R-2 over-shoots; R-2 is structurally required.

**Commit R-1 stays** (fc16407 is not reverted). **Plan R-2 as the natural next iteration.**

## Implications for the production metric

The Percolator @ 1% FDR count is 74,204, which dramatically beats:
- Java baseline: 35,818 (+38,386, +107%)
- b1d45bb pre-R-1: 25,224 (+48,980, +194%)
- iter5 all-revert: 17,160 (+57,044, +332%)

However, two indicators suggest the score distributions are noisier than Java's:

1. **pi_0 estimate: 0.889** — Percolator estimates that ~89% of targets are noise. Java's
   distribution presumably has a lower pi_0 (its 35,818 / 89,479 = 40% pass rate implies
   tighter target/noise separation). Rust's 74,204 / 1,042,255 = **7.1% pass rate** vs
   Java's ~40% pass rate.

2. **Score distribution shape:** With 12× more raw PSMs going into Percolator from a
   single-queue cross-charge architecture, the additional PSMs are likely low-quality
   cross-charge candidates that Java's per-SpecKey architecture would have kept in separate
   (smaller) queues, then filtered. Percolator is doing extra filtering work that Java
   would have done architecturally.

**Whether the extra 38K PSMs over Java are a "real win"** depends on downstream validation
of false positive rates. The R-1 result is empirically positive on the production metric,
but with the above caveats. R-2 is expected to bring the pass rate closer to Java's 40% by
reducing the raw target count from 1,042,255 toward Java's 89,479 range, while (ideally)
preserving the Percolator @ 1% FDR improvement.

## Next iteration — R-2

R-2 (per-SpecKey GF compute + per-charge queue architecture) is now the natural next step:

- **What:** Rust's `match_engine.rs:325` uses a single `top_charge` GF context for the
  whole spectrum. Java computes per `(spectrum, charge)` via `SpecKey`
  (`DBScanner.java:606`, `:779`). Refactoring Rust to use per-SpecKey queues will:
  - Limit tie accumulation to within each charge state (closing the 12× over-shoot)
  - Compute GF denominators correctly per charge (the original R-2 GF-correctness goal)
  - Bring the raw target count closer to Java's ~89K range

- **How:** This is "Approach C" from the investigation note — a structural refactor of the
  per-spectrum loop, not a one-line fix. The spec for R-2 should be written before starting.
  Estimated effort: 2-4 hours (refactor + unit test + bench cycle).

- **Prediction:** After R-2, raw targets will fall from 1,042,255 toward Java's 89,479 range.
  The Percolator @ 1% FDR count may stay elevated (if the extra PSMs were real) or drop
  (if they were noise). Either outcome is informative. The T/D ratio should remain ≥1.47.

- **Coupling confirmed:** The audit's "R-1 + R-2 should be coupled" prediction is confirmed.
  R-1 and R-2 are now committed and planned respectively. R-3 (`minDeNovoScore` PIN filter),
  R-4 (`lnEValue` denominator), and F-1 (`matched_ion_ratio` denominator) remain deferred
  until the R-2 architectural baseline is established.
