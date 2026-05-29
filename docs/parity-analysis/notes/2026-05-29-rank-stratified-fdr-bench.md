# Chimeric rank-stratified FDR bench — REFUTED (wide-window inflates rank-1 too)

**Date:** 2026-05-29
**Branch:** `feat/chimeric-dda-plus`
**Design:** `superpowers/specs/2026-05-29-chimeric-rank-stratified-fdr-design.md`
**Method:** split chimeric PIN by rank (SpecId `{spec_id}_{scan}_{rank}_{rowidx}`,
`rank` = 2nd-to-last `_`-field), run Percolator **separately** on rank-1 vs rank≥2,
sum targets@1%. Run on both the Phase-3 rescored PIN (b2) and the non-rescored
Phase-1 PIN (b1, via `MSGF_CHIMERIC_NO_RESCORE=1`).

## Results (@1% FDR, separate Percolator per stratum)

| PIN | dataset | rank-1 (decoy_frac) | rank≥2 (decoy_frac) | total | off | Java |
|---|---|---|---|---:|---:|---:|
| Phase-3 | PXD | 17,473 (0.343) | 778 (0.919) | 18,251 | 14,808 | 14,974 |
| Phase-3 | Astral | **51,579 (0.392)** | 25,535 (0.894) | **77,114** | 36,715 | ~36,271 |
| Phase-1 | PXD | 17,683 (0.360) | 344 (0.921) | 18,027 | 14,808 | 14,974 |

(rank-1 is **rescore-invariant** — it claims peaks first, never gets re-scored —
so Phase-1 and Phase-3 rank-1 agree: PXD 17,683 ≈ 17,473. The b1 Astral rank-1
equals the b2 51,579 by the same argument.)

## Verdict: rank-stratified FDR does NOT make chimeric trustworthy

The total barely changes vs pooled PSM-level (Astral 77,114 ≈ pooled 77,444). The
canary still fires. **The decisive finding is the rank-1 stratum itself:**

- Astral rank-1 = **51,579** (+40% over off 36,715, +42% over Java) — yet its
  **decoy fraction is HEALTHY (0.392, better than off's 0.517).**
- PXD rank-1 = 17,683 (+19% over off 14,808), decoy_frac 0.360 (≈ off 0.396).

**Healthy decoy fraction + implausible absolute count = coincidental-target
inflation that target-decoy FDR cannot see.** A wide isolation window lets a
real-DB-sequence peptide (from elsewhere in the window) win rank-1 by coincidental
fragment matches. It is a *target* with no decoy counterpart, so TDC counts it as
~99% correct — the decoy fraction stays healthy while the count inflates. 51k IDs
on a 15-min Astral run (truth ≈ 36k) is not real.

**The inflation is in the wide-window SELECTION, not just the rank≥2
multi-emission.** Rank-stratification therefore cannot fix it: even restricting to
one best peptide per scan, the wide window admits ~40% spurious-but-real-sequence
rank-1 matches. The rank≥2 strata are separately pure junk (decoy_frac ~0.9).

## Why this closes the post-hoc-FDR line

No post-hoc re-stratification of the chimeric PIN can recover trustworthiness,
because the bad IDs are coincidental *targets* baked into which peptides win the
wide-window search — and TDC is blind to them. The only remaining levers are
**pre-FDR admission control**: hard MS1 precursor-evidence gating (Phase 2 showed
the *soft* feature insufficient; a *hard* gate was the untested mitigation) and/or
an entrapment/competition FDR that models per-scan multiple testing — both
substantial, and per [[merge-gate-beat-java]] chimeric still wouldn't move TMT.

## Chimeric line status (cumulative negative results)

| Phase | Lever | Result |
|---|---|---|
| 1 | wide-window multi-emission | inflates FDR (Astral +94–97%) |
| 2 | MS1 isotope-KL soft feature | insufficient (hard-filter test: still +89%) |
| 3 | fragment competition + residual SpecE re-score | canary fires (+111%); per-PSM rescore moves targets+decoys together |
| — | rank-stratified separate FDR | **this note** — rank-1 itself inflated by coincidental targets (TDC-invisible) |

**Recommendation: shelve chimeric.** All four levers refuted; the remaining
pre-FDR-gating path is large and still gate-irrelevant (no-op Astral, net-negative
TMT). Redirect to the actual gate blocker (TMT). All machinery kept on
`feat/chimeric-dda-plus`; `--chimeric off` byte-identical; nothing shipped.

## References

- `2026-05-29-chimeric-phase3-bench-canary-fails.md`
- `2026-05-28-chimeric-phase2-bench.md`
- `2026-05-29-chimeric-fragment-overlap-diagnostic.md`
