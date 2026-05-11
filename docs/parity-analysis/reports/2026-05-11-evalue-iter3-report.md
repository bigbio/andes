# E-Value Iter 3 Results: Mod-Aware Hypothesis Rejected, Ratio Gap is Structural

**Iteration 3 commit:** `47893d7` — replaced `HashSet<Vec<u8>>` with `FxHashSet<u64>` in `SearchIndex.ensure_distinct_peptide_counts`.

## Results

| Metric | Pre-T2 | T2-5 | Iter 3 (47893d7) | Status |
|---|---|---|---|---|
| match_spectra wall | ~64s | 242s | 195s | partial recovery (~20% better than T2-5, still 3× baseline) |
| Total wall | ~5m | 9m17s | 7m53s | partial recovery (~15% improvement) |
| Pin rows | 37,113 | 37,113 | 37,113 | unchanged |
| Percolator @ 1% FDR | 14,798 (Java) | 14,850 | 14,850 | unchanged from T2-5 |
| EValue ratio median (Rust/Java) | — | 0.0368 | 0.0368 | **UNCHANGED** |
| EValue % within ±5% | — | 0.0% | 0.0% | **UNCHANGED** |

## Key Finding: Mod-Aware Hypothesis Rejected

The commit message at `47893d7` claims "Mod-aware counting: seen-set now hashes residues + mod positions/masses". **This is incorrect.**

The actual code hashes only residue bytes. The discovery step found that Java's `CompactSuffixArray.computeNumDistinctPeptides` is also bare-residue-only — therefore, adding mod-awareness was unnecessary.

The code comments correctly describe bare-residue behavior; only the commit message is misleading.

## Structural Divergence: EValue Denominator

The 0.0368 median EValue ratio gap is **not addressable by optimizing the seen-set**. The ratio is unchanged because the root cause is **structural**: Java and Rust use a different peptide-count divisor in the EValue formula.

**Next investigation:** Compare Java's exact EValue computation (`DBScanner.scan` or `MSGFPlusMatch.computeEValue`) to Rust's (`psm.e_value = psm.spec_e_value * num_distinct_peptides_at_length`). If Java uses total-DB peptides (or peptides-at-all-charges, not per-length), a formula fix closes the gap.

## Wall Regression Status

T2-5 → Iter 3 is a ~15% improvement (9m17s → 7m53s), but still 3× the pre-T2 baseline (~5m). The `FxHashSet<u64>` optimization reduces allocations but does not address the load-bearing regression. Further recovery requires identifying why Phase 6's enumerate_candidates pass became 3× slower.
