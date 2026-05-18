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

## Headline correction — the 74,204 PSM "win" is artifact, not signal

Initial framing of this result emphasized the Percolator @ 1% FDR count (74,204 vs Java
35,818 = +107%) as a production improvement. **That framing is superseded** by the deduped
peptide-level comparison below. The +38K Percolator-row "win" is row-count inflation from
tied retention without Java's per-SpecKey + dedup + protein-index aggregation pipeline.
By the metrics that actually matter (unique identifications, spectrum coverage, score
confidence), **Rust + R-1 alone is worse than Java**, not better.

## Honest comparison — deduped peptide and scan coverage at q≤0.01

| Metric | Java | Rust iter6 | Δ (Rust − Java) |
|---|---:|---:|---:|
| PSMs (raw rows in `.target.psms.txt`) | 35,818 | 74,204 | +38,386 |
| **Distinct (scan, peptide) pairs** | **35,818** | **26,934** | **−8,884** |
| **Unique scans identified** | **35,411** | **25,972** | **−9,439** |
| **Unique peptides identified** | **22,925** | **18,345** | **−4,580** |
| Avg PSMs per identified scan | 1.01 | 2.86 | — |

The +38K "extra" PSMs are tied duplicates per scan in Rust. After deduping to distinct
(scan, peptide) assignments, **Java identifies 8,884 more PSMs**, **9,439 more spectra**,
and **4,580 more peptides** than Rust. R-1 alone produces redundancy, not coverage.

### Set overlap on (scan, peptide) at q≤0.01

| Bucket | Count |
|---|---:|
| Both engines agree (same scan, same peptide) | 24,622 |
| **Java-only (Rust missed)** | **11,196** |
| Rust-only (extra unique to Rust) | 2,312 |

Of the 35,818 distinct Java PSMs at 1% FDR, Rust agrees on 24,622 (69%) and **uniquely
loses 11,196**. Rust gains 2,312 PSMs Java doesn't find. Net: −8,884.

### Percolator score comparison on the 24,622 overlapping PSMs

| | Java | Rust |
|---|---:|---:|
| Mean Percolator score | 1.129 | 0.215 |
| Stdev | 0.683 | 0.190 |
| Pearson r (Java vs Rust) | — | 0.786 |

Rust's score distribution is **compressed** (mean ~5× smaller, stdev ~3.6× smaller).
99.2% of overlapping PSMs have a *lower* Percolator score in Rust than in Java (median
delta = −0.84). Score-delta histogram across 24,622 overlapping PSMs:

| Rust − Java | Count | % |
|---|---:|---:|
| < −2.0 | 962 | 3.9% |
| −2 to −1 | 8,566 | 34.8% |
| −1 to 0 | 14,890 | 60.5% |
| 0 to +1 | 204 | 0.8% |
| > +1 | 0 | 0.0% |

The compressed-discriminator pattern can't be fixed by R-2 alone — see "Remaining
divergences" below.

## Mechanism interpretation — why R-1 over-shoots in Rust's architecture

### The spec_e_value=1.0 sentinel tie

Before `compute_spec_e_values_for_spectrum` runs (i.e., at raw-score insertion time), **every
PSM has `spec_e_value = 1.0`** — a sentinel value. The `PsmMatch::cmp` ordering puts
`spec_e_value` first (smallest wins) and `score` second (largest wins). So at insertion time
the primary key ties for *every* PSM in the queue: all have `spec_e_value = 1.0`.

The secondary `score` key is an `i32` rounded from a float. Integer-valued scores produce
many ties within a spectrum's candidate set, particularly at the lower end of the score
distribution. R-1's `Equal` branch inserts all tied PSMs without eviction.

### One queue per spectrum vs Java's per-SpecKey pipeline

Java's pipeline has **four layers** that together limit how many tied PSMs reach the .pin:

1. **Per-SpecKey raw-score retention with tie keep** (`DBScanner.java:540`) — one queue per
   `(spectrum, charge)`. Ties at the worst score are kept, but only within that charge.
2. **Pre-merge pepSeq + score dedup** (`DBScanner.java:719`) — collapses degenerate matches
   (same peptide, same score) BEFORE the spectrum-level merge.
3. **GF / SpecEValue computed per SpecKey** (`DBScanner.java:606`, `:779`) — each
   (spectrum, charge) gets its own GF distribution; SpecE is calibrated per-charge.
4. **Spectrum-level merge by SpecEValue** with tie keep (`DBScanner.java:745`) — when
   merging the per-charge per-SpecKey queues into a per-spectrum result, ties on SpecE are
   kept (mirrors layer 1 at the merge stage).

Rust has **only layer 1** (and only post-R-1; pre-R-1 it had a *broken* layer 1 too):

- **Layer 1 (Rust):** one `TopNQueue` per spectrum, shared across all charges
  (`match_engine.rs:325`). R-1 made the tie-keep correct *within the shared queue*, but
  the queue's scope is wrong (spans charges).
- **Layer 2 (Rust):** absent. Same peptide can survive as multiple rows at the same scan.
- **Layer 3 (Rust):** broken. `match_engine.rs:325` picks one `top_charge` GF context
  for the whole spectrum (across all charges).
- **Layer 4 (Rust):** N/A — there is no per-spectrum merge step because the queue is already
  per-spectrum.

### Protein-index aggregation (new finding, not in the 2026-05-16 audit)

A fifth contributor to Rust's row inflation: **Java aggregates multiple protein indices
into ONE `DatabaseMatch` via `addIndex(...)` (`DatabaseMatch.java:75`), then writes them
all in one PIN row** (`DirectPinWriter.java:237` — `for (String acc : proteins.accessions)
row.append(acc)`). Rust resolves exactly one accession from one `Candidate` and writes one
`Proteins` cell (`row_context.rs:47`, `pin.rs:462`). When the same peptide matches multiple
proteins (common for shared tryptic peptides in target+decoy concat), Java emits 1 row
listing all proteins; Rust emits N rows, each with one protein. Some fraction of the 47,270
duplicate (scan, peptide) PSMs in Rust's PIN come from this multiplicity, not just from
R-1's tie keeping.

### The magnitude explained

The 11.6× raw target count (1,042,255 vs Java's 89,479) is not R-1 doing "more than
expected" in terms of correctness. R-1 is the correct fix to the tie-dropping bug. The
over-shoot is the compound effect of:
- Missing layer 2 (pepSeq + score dedup) → same peptide kept at multiple charges
- Missing layer 3 (per-SpecKey GF) → cross-charge candidates all compared against the same
  GF distribution
- Missing layer 4 / cross-charge queue → ties accumulate across all charges, not just
  within one
- Missing protein-index aggregation → same (scan, peptide) emitted once per protein assignment

R-1 fixes ONE divergence faithfully. The remaining four divergences amplify R-1's effect.

## Decision per spec outcome table

The spec's outcome table reads gates on raw target count, T/D ratio, and wall time. All
three pass: 1,042,255 ≥ 86,674; 1.965 ≥ 1.47; 23:52 ≤ 26:00. Per the spec, this is the
"Hypothesis validated" branch.

**But the spec's outcome table is incomplete.** It treats "raw target count moved toward
Java" as sufficient evidence, without checking whether the movement is real coverage or
duplicate inflation. Updated decision logic — empirically driven:

| Empirical | Verdict |
|---|---|
| Raw target count gates pass | ✓ per spec |
| **But distinct (scan, peptide) count is LOWER than Java** | ✗ — R-1 alone doesn't improve coverage |
| **And unique scans + peptides are LOWER than Java** | ✗ — R-1 alone regresses coverage |
| **And Rust Percolator scores are systematically lower** | ✗ — R-1 doesn't fix discrimination |

**R-1 stays committed** (fc16407 is not reverted) because the fix is correct in isolation
and the audit's tie-keep prediction is empirically validated as one component of a larger
problem. But the production-metric framing — "+38K PSMs vs Java" — is **not a real win**
and should not be cited as a success.

## Next iteration — full R-2 sequence (revised)

R-2 is no longer "per-charge GF compute" alone. To match Java, Rust needs all four
retention layers + the protein-index aggregation:

1. **R-2.1 — Per-SpecKey raw-score queue.** Replace the single `TopNQueue` per spectrum
   with one per `(spectrum, charge)`. Insertion happens against the per-charge queue.
   Cite: `DBScanner.java:534`.
2. **R-2.2 — Pre-merge pepSeq + score dedup.** Before the per-spectrum merge, collapse
   degenerate matches (same peptide, same score, possibly different protein indices) into
   one entry. Cite: `DBScanner.java:719`.
3. **R-2.3 — Per-SpecKey GF / SpecEValue.** Compute the GF distribution per
   `(spectrum, charge)`, not once per spectrum. Each per-charge queue gets SpecE calibrated
   against its own GF. Cite: `DBScanner.java:606`, `:779`.
4. **R-2.4 — Spectrum-level merge with SpecE tie keep.** Merge the per-charge queues into
   per-spectrum output, keeping SpecE ties at capacity. Cite: `DBScanner.java:745`.

Plus a separate item closely tied to R-2.2:

5. **R-2.5 — Protein-index aggregation.** When multiple Candidates resolve to the same
   peptide and score (across different proteins), aggregate the protein accessions into one
   PIN row (one `Proteins` cell with tab-separated accessions, like Java's
   `DirectPinWriter.java:237`). Requires `DatabaseMatch.addIndex`-style accumulation in
   Rust's data model.

### Will R-2 alone close the Astral gap?

**No** — the score-comparison evidence above shows that even on PSMs both engines agree on,
Rust's Percolator scores are 99.2% lower than Java's. R-2 will fix the row-inflation
artifact and bring counts toward Java's range, but the **compressed discriminator** is a
separate problem rooted in:

- **enzN/enzC/enzInt zero-stubbed in Rust** (`pin.rs:418` vs `DirectPinWriter.java:195`)
- **longest_y_pct denominator** (Rust uses `n`; Java uses `n − 1`; `match_engine.rs:745` vs
  `PSMFeatureFinder.java:95`)
- **Charge-1 b/y-only feature extraction in Rust** (`match_engine.rs:628` vs Java's full
  scorer ion model `PSMFeatureFinder.java:148`)
- **Rust skips Java's `minDeNovoScore` PIN-row filtering** (`pin.rs:251` vs
  `DirectPinWriter.java:130`)

Per `2026-05-18-piecewise-fixes-dont-work.md` (the iter3-5 lessons-learned doc), these
feature/scoring fixes **regressed Astral by 8K PSMs when applied without the retention
layer**. The hypothesis is that fixing them on top of a correct R-2 architectural baseline
won't regress (because Percolator will be operating on Java-shaped data). That hypothesis
is testable but not validated yet.

### Effort estimate

R-2 is a multi-file refactor across `match_engine.rs`, `psm.rs`, and likely
`scored_spectrum.rs`. Conservative estimate: **1-2 days** for R-2.1 through R-2.4 + R-2.5,
including a strengthened `match_engine_java_parity` test that gates on distinct (scan,
peptide) count agreement with Java (not just top-1 peptide identity).

### Coupling and order

- R-1 (committed, fc16407) and R-2 (planned) are confirmed coupled.
- R-2.5 (protein-index aggregation) is coupled to R-2.2 (pepSeq dedup) — both touch the
  data model for multi-protein matches.
- R-3 (`minDeNovoScore` PIN filter), R-4 (`lnEValue` denominator), F-1 (`matched_ion_ratio`
  denominator) remain deferred until R-2 establishes the architectural baseline.
- The compressed-discriminator feature fixes (audit-tier C-4/C-5/C-5b) require R-2 first;
  attempting them on the current single-queue base recreates iter5's −8K regression.
