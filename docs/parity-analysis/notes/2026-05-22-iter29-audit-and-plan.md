# iter29 audit + next-phase plan

_2026-05-22. Comprehensive code review + perf analysis after iter29's main_ion fix landed (+379 Astral PSMs). Two parallel subagent passes, one full 3-dataset bench. This doc consolidates findings into a prioritized plan for closing the remaining 11.6% Astral gap **and** making Rust faster than Java._

## TL;DR

- **Astral 1% FDR is 31,677 / Java 35,818 (gap 11.6%, was 26% at iter16).**
- **Rust is already faster than Java on PXD001819 (1:07 vs 1:20).** The "Rust 15-45 min" memory entry was stale — Rust is now within 25% of Java on Astral wall time (7:32 vs Java's typical 5:50-6:46).
- Code review surfaced **2 HIGH-severity correctness bugs** (deconvolution gating + prob_peak source) that plausibly account for a chunk of the remaining Astral gap and are <50-line fixes.
- Perf review found a clear **~10-25% wall-time win** in 4 changes totaling <300 LOC, none of which alter bit-identity.

## Section 1 — bench results (3-dataset, iter29 vs Java)

| Dataset | Engine | Wall (8 threads) | targets_total | 1% FDR | gap vs Java |
|---|---|---:|---:|---:|---:|
| PXD001819 | Java | 1:20 | 28,037 | 14,989 | — |
| PXD001819 | Rust iter29 | **1:07** ✓ | 28,038 | 14,751 | -1.6% (-238) |
| Astral | Java | 5:49 | 89,479 | 35,818 | — |
| Astral | Rust iter29 | 7:32 | 92,781 | 31,677 | -11.6% (-4,141) |
| TMT | Java | 3:07 | 28,790 | 10,194 | — |
| TMT | Rust iter29 | 3:26 | 27,683 | **11,091** ✓ | **+8.8% (+897)** |

**Observations:**
- **TMT: Rust BEATS Java by +897 PSMs (+8.8%)** at 1% FDR. Rust enumerates FEWER total targets (27,683 vs 28,790) but gets MORE confident PSMs. This is the first dataset where Rust is materially better than Java post-Percolator.
- **PXD001819: Rust 16% faster wall** (1:07 vs 1:20) but slightly fewer PSMs (-238). The iter29 main_ion fix slightly hurt PXD001819 (low-res CID). Hypothesis: CID_LowRes_Tryp.param's dominant ion may be b-ion (low-res CID has noisier y-ions), so the "fix" went the wrong way for that dataset. Need to verify by checking which ion `main_ion_from_param` picks for PXD001819 partitions.
- **Astral: gap 11.6% (4,141 PSMs)** is the largest single dataset divergence. DeNovoScore is now at agreement-bucket parity (median Δ 0) per the iter29 pin-diff. Remaining gap is in **top-1 peptide selection**: ~40% of scans pick a different peptide across engines (`both_target_diff_peptide=18,893` + `java_target_rust_decoy=16,636` + `rust_target_java_decoy=13,633` = 49,162 / 121,681).
- **Wall summary**: Rust is faster on PXD001819, slower on Astral by 30%, slower on TMT by 10%. The Astral wall gap is the most actionable target.

## Section 2 — code review findings (correctness)

### HIGH — C-1: deconvolution `charge > 2` guard is wrong

`crates/scoring/src/scoring/scored_spectrum.rs:229`:
```rust
if param.apply_deconvolution && charge > 2 {  // ← Java has no charge guard
```

Java's `NewScoredSpectrum.java:76` unconditionally deconvolutes when `applyDeconvolution()` is true. For Astral charge-2 spectra (a large fraction), Rust skips deconvolution and uses the raw peak list for both `prob_peak` and the prefix/suffix node-score cache. The downstream `deconvolute_spectrum` is already safe for charge ≤ 2 (its `2..charge.min(4)` inner loop runs 0 iterations), so removing the guard produces `deconv_peaks` equal to the input — identical to Java's behavior.

**Fix:** drop the `charge > 2` clause. Risk: low (no-op for charge ≤ 2 mathematically, just consistency with Java). LOC: 1.

### HIGH — C-2: `prob_peak` computed BEFORE deconvolution

`crates/scoring/src/scoring/scored_spectrum.rs:196-198`. `prob_peak = kept_count / approx_num_bins` is set from pre-deconv peak count. Java sets `probPeak = spec.size() / approxNumBins` AFTER `spec = spec.getDeconvolutedSpectrum(...)`. For charge ≥ 3 spectra where deconvolution materially changes peak count, `prob_peak` is biased high (more peaks per bin → looks like b-ion present everywhere). This affects every `ion_existence_score` call in the GF DP.

**Fix:** after the deconv block (line 235), if `deconv_peaks.is_some()` recompute `prob_peak = deconv_peaks.len() as f32 / approx_num_bins.max(1.0)`. Risk: low (only fires for charge ≥ 3 with apply_deconvolution=true; matches Java's order). LOC: 5.

### MED — C-4: cleavage credit possible double-count

`crates/search/src/match_engine.rs:261-287`. The per-candidate `compute_cleavage_credit` adds `neighboring_aa_cleavage_credit` and `peptide_cleavage_credit` to every PSM's RawScore. The GF DP separately applies a `neighboring_aa_cleavage_credit/penalty` adjustment to the final score distribution (`generating_function.rs:615-631`). Java's `PrimitiveAminoAcidGraph` applies cleavage on source/sink EDGES only (`addCleavageFromSource = direction == enzyme.isNTerm()`). For C-term enzymes (trypsin), `addCleavageFromSource = false`. The Rust per-candidate cleavage addition may be in the wrong place.

**Fix:** trace per-PSM Java RawScore vs Rust score_psm + cleavage_credit. We already proved score_psm is bit-exact (iter28 trace, scan 47106) and Java's pin RawScore = 73 = 61 + 4 (cleavage) + 8 (edge). Rust pin RawScore = 65 = 61 + 4 (cleavage). So cleavage matches at +4 in both. But for non-tryptic / semi-tryptic cases this may differ — needs a non-tryptic test scan. Risk: medium. LOC: depends on findings.

### LOW — C-5: `dedup_pepseq_score` may collapse protein-terminal flags

`match_engine.rs:338-346`. When peptides at multiple proteins share the same sequence + score, dedup keeps the first one. If a true protein-N-term match is dedup'd against a non-N-term match, the GF DP later uses the wrong source AA list. Latent landmine.

**Fix:** OR `is_protein_n_term` / `is_protein_c_term` across the collapsed candidate set when dedup'ing. LOC: ~20.

### Tests gap

The code-review subagent flagged 5 missing tests — most critically:
- T-1: no test exercises `deconvolute_spectrum` for charge-2 spectrum
- T-2: no test compares `prob_peak` before/after deconvolution
- T-3: `setup_score_threshold` not unit-tested in isolation

## Section 3 — perf review findings

Rust at iter29 Astral runs at **647% CPU on 8 cores** (81% saturation, healthy parallelism). The 25% wall gap vs Java is per-thread cost, not parallelism.

| Opt | LOC | Expected wall reduction | File |
|---|---:|---:|---|
| **#1 hoist `env::var("MSGF_TRACE_PEP")` out of hot path** | ~10 | 3-8% | `psm_score.rs:150` |
| **#2 pipeline mzML read with scoring** | ~80 | 5-12% | `bin/msgf-rust.rs:406-499` |
| **#3 SmallVec for per-PSM matched arrays** | ~20 | 1-3% | `match_engine.rs:733-829` |
| **#4 cache `ions_for_partition_slice` per spectrum** | ~50 | 2-5% | `match_engine.rs:857` |
| **#5 pool `observed_by_mass` Vec into GF arena** | ~30 | 1-2% | `primitive_graph.rs:810` |

Sum: 190 LOC, expected total **12-30% wall reduction** ≈ 5:15-6:35 Astral wall after all five. That brings Rust ahead of Java's 5:50-6:46 typical Astral wall.

## Section 4 — prioritized roadmap

### Phase A: HIGH-impact correctness (iter30, this week)

1. **C-1 + C-2 deconvolution fixes.** Both LOC small, both gate on `apply_deconvolution=true` which is the Astral param. Land together, bench Astral.
2. **C-5 protein-terminal flag OR-in dedup.** Cheap, latent. Land with #1.
3. **PXD001819 regression triage.** Check which ion `main_ion_from_param` picks for CID_LowRes_Tryp.param at each partition. If b-ion is correct (low-res CID y-ions are noisier), restrict the fix to params where the suffix ion truly dominates. Worst case: gate the fix by `param.activation == HCD`.
4. **TMT bench** (when it lands).

Expected outcome: Astral 1% FDR closer to Java; PXD001819 regression closed. Target: gap ≤ 8%.

### Phase B: perf cluster (iter31, day 1-3)

1. **#1 env::var hoist** (30 min, zero risk).
2. **#3 SmallVec matched arrays** (1-2 h, zero risk).
3. **#4 ions_for_partition cache** (~2 h, zero risk — reuses existing `segment_partition_cache`).
4. Bench Astral: target ≤6:30 wall.

### Phase C: deeper perf (iter32, day 4-5)

1. **#2 pipeline parse+score** with `crossbeam::channel`. ~80 LOC, low risk.
2. **#5 arena-pool observed_by_mass.** ~30 LOC, zero risk.
3. Bench all three datasets. Target: Astral ≤6:00 wall, faster than Java on every dataset.

### Phase D: top-1 selection convergence (iter33+, when stable)

After Phase A lands, re-run pin-diff. Focus on the `both_target_diff_peptide=18,893` bucket. The agreement-bucket DeNovoScore is now at parity, so disagreement is driven by:
1. **Tie-breaks** (~25% of cases): multiple peptides at the same RawScore + DeNovoScore. Whoever wins is Percolator-influencing. Java's tie-break order differs from Rust's. Audit `TopNQueue::push` ordering vs Java's `PriorityQueue` insertion order.
2. **score_psm + score_psm-with-edge disagreements**: same node-score but different edge-score path (post-iter29 the direction is right, but per-edge values may still diverge for charge-2 spectra where deconvolution is skipped → C-1 + C-2 may close this).
3. **lnSpecEValue precision differences** affecting Percolator's primary feature.

## Section 5 — work NOT to do

- ❌ Modify RawScore to include edge_score (iter17/18 regressed -8K PSMs; the iter19 EdgeScore separate PIN column is correct).
- ❌ Bit-exact error stats (iter23 regressed -1,404 PSMs).
- ❌ Cleavage efficiency 0.95 → 0.99999 (penalty change only; no DeNovoScore impact per analysis in `2026-05-21-iter27-pin-diff.md`).
- ❌ num_distinct semantic fix (item #2 in known-divergences.md). Item is documented but `lnEValue Δ ≈ -4.56` is already what Percolator was trained against — modifying it likely regresses.

## Section 6 — open questions

1. **Astral wall**: Rust 7:32 vs Java 5:49 = 30% slower. Phase B+C perf cluster targets 6:00 (parity).
2. **TMT direction**: Rust BEATS Java by +897 PSMs (+8.8%) at 1% FDR. Why? Is this stable or an artifact of the param choice? Verify with another TMT dataset (or a replicate of PXD007683).
3. For C-4 (cleavage double-count): need a non-tryptic test scan to confirm/deny.
4. For the PXD001819 regression: which ion type does CID_LowRes_Tryp.param prefer? If `main_ion_from_param` now picks y-ion but Java picks b-ion for CID-LowRes, the fix went the wrong way for that dataset. Quick check: dump `main_ion` for a few CID_LowRes partitions and compare to Java's `getMainIonType`. **Decision tree:**
   - If both engines pick same ion → iter29 fix is right; PXD001819 regression has another cause (look at C-1/C-2 deconvolution path or look at low-res rank-dist tables).
   - If Java picks b-ion and Rust now picks y-ion → guard `main_ion_from_param` by `param.activation` (HCD/CID + High/Low) or by `getDataType().getInstrument()`.

## Section 7 — concrete next actions

1. **Land C-1 + C-2** (deconvolution gate + prob_peak ordering). Single small commit. Bench Astral.
2. **PXD001819 diagnostic**: print `main_ion` for the (charge=2, parent_mass≈1500, seg=0) partition on both engines. 10-min task; informs whether to keep the iter29 fix universally or gate by instrument.
3. **Compose iter30 plan branch**: C-1 + C-2 + C-5 + PXD001819 fix. Bench all three datasets. Goal: Astral gap ≤ 8%, PXD001819 back to ≥14,850, TMT preserved.
4. **iter31 perf cluster**: env::var hoist + SmallVec + ion-cache. Target Astral wall ≤ 6:30.
5. **iter32 deep perf**: pipeline parse+score + arena pool. Target Astral wall ≤ 6:00 (parity).
