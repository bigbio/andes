# Tree-count default: measured speedup, and an audit of the entrapment metric (2026-09-04)

Two results. The first is a clean speed measurement of the default changed in #47. The
second is an audit finding: **the entrapment-FDP metric quoted in the README overstates
what two of its three columns actually measure**, and one of its scaling assumptions is
wrong. Both were produced on the same host on the same day.

Provenance for everything below: commit `1b8520f8` (merged `main`), benchmark VM, 8
threads, Percolator 3.7.1 biocontainer `--seed 42 -Y`, counts at `q ≤ 0.01`.

## 1. `--gbdt-max-trees` default: 1.64x faster, identification-neutral

Same binary, same host, same day, same data — the tree-count default is the only
variable, so this is a controlled A/B rather than a comparison across runs:

| Astral arm | trees | wall | PSMs @ q≤0.01 |
|---|---:|---:|---:|
| previous default | all (300) | 400 s | 38,402 |
| **new default** | **100** | **244 s** | **38,394** |

**1.64x faster for −8 PSMs (−0.02%)** — inside run-to-run noise. An independent run of the
new default earlier the same day gave 262 s / 38,394, so treat the wall time as ~250 s
± 8%, not as a single-second-precision figure.

UPS1, same binary and session: **50 s, 15,838 PSMs**.

### What this does NOT establish

The previously published head-to-head — "450 s vs 217 s on Astral" (2026-08-28) — predates
this default **and** a different binary. It has not been re-run, and **Comet is no longer
installed on the benchmark host**, so no refreshed andes-vs-Comet wall-time ratio is
claimed here. The 1.64x above is andes against its own previous default, nothing more.

## 2. Entrapment-metric audit — read before quoting a "1% true FDP"

The README describes its headline table as "PSMs at 1% true entrapment-FDP (1:1 entrapment
database; FDP = 2·ENT/total)". Measured against the databases actually in use, that
description does not hold for two of the three columns.

### Astral is not entrapment-validated

The Astral run against `hye.fasta` returns **0 entrapment hits, because the database
contains no entrapment sequences at all** (31,889 sequences, zero `ENTRAP_`). An entrapment
FDP cannot be computed from it, at any threshold. `2026-08-28-config-matrix.md` already
stated this ("not entrapment-validated and should not be quoted as though it were"); this
run confirms it directly.

A genuine 1:1 entrapment cross-check on Astral *was* performed in `2026-06-15`, using a
purpose-built database (real HYE proteins plus an equal set of foreign sequences). That
measurement stands on its own terms; it is simply not the configuration benchmarked since.

### UPS1's entrapment database is not 1:1, so the factor of 2 is wrong

`yeast_entrap.fasta` is **6,733 target (yeast) against 4,531 entrapment (E. coli)**. The
correct estimator is `FDP = ENT/total × (1 + T/E)`, not `2 × ENT/total`:

| basis for T/E | ratio | factor | measured FDP at nominal 1% |
|---|---:|---:|---:|
| protein count | 1.49 | 2.49 | **2.61%** |
| searchable tryptic peptides (734,280 vs 303,537) | 2.42 | 3.42 | **3.58%** |
| *(the assumed 1:1)* | *1.00* | *2.00* | *2.10%* |

Measured here: 166 entrapment hits on 15,838 PSMs. The protein-count figure (2.61%)
reproduces the 2.57% independently derived in `2026-08-28-config-matrix.md` from a separate
run, so the two agree.

**Peptide space is the better basis** — it is what the search actually samples, and E. coli
proteins are shorter than yeast ones, which is why the two bases differ by 37%. On either
basis the nominal 1% on UPS1 is really **2.6–3.6% true FDP**, not 1%.

### Consequence

Cross-arm and cross-engine comparisons under one methodology are unaffected: every arm
carries identical scaling, which is what those tables are for. **Absolute "1% true FDP"
claims are not supported** for Astral (no entrapment component) or UPS1 (mis-scaled, and
~2.6–3.6% in truth). Fixing this properly means rebuilding the UPS1 entrapment database
near 1:1 and adding an entrapment component to the Astral database — filed as follow-up.

## 3. Provenance gap in the README headline table

The README's headline counts — 36,873 / 11,163 / 15,061 for andes and 28,401 for Comet —
**do not appear in any document in `docs/benchmarks/`**. The nearest documented run
(`2026-06-26-owngeometry-ab.md`) reports 36,730 / 11,215 / 14,919, close but different.
Those numbers should be re-derived from a recorded run or replaced with measured ones
before being quoted further; they are left in place here, flagged, rather than silently
restated or silently deleted.
