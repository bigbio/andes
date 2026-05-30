# Fragment-index prefilter — PXD (low-res CID) FAILS both gates (coarse-bucket degeneration)

**Date:** 2026-05-30
**Branch:** `feat/chimeric-dda-plus` (impl T1-T3: FragmentIndex+FragmentVoter wired, commit 24228572)

## T4 PXD recall+speed gate result (chimeric NO_RESCORE, normal DB)

| arm | rows | wall | @1% PSMs | T/D |
|---|---:|---:|---:|---:|
| brute (`--chimeric-frag-index off`) | 365,661 | 1:36 | 18,234 | 1.20 |
| index (`on`) | 280,133 | **28:24** | **6,207** | 1.74 |

index entrapment FDP (on): 0.81% raw / 1.61% combined (still ~controlled).

**Both gates FAIL: 18× SLOWER (target: faster), 34% recall (target ≥99.5%).**

## Root cause: coarse-bucket degeneration on low-res CID

PXD001819 = LTQ-Orbitrap Velos, CID ion-trap MS2 → `is_high_resolution()==false` →
fragment bin width 0.5 Da. Coarse bins → each bin holds a huge candidate list →
the voter's `touched` set balloons toward ALL candidates per spectrum → `top_k`'s
collect+filter+sort is O(n·log n) PER spectrum, worse than brute-force scoring of
the mass-window candidates. Same "Tier-1 overhead exceeds Tier-2 savings" failure
the Java fragment-index abandonment hit (ABANDONED-2026-04-20.md), reproduced on
low-res data. The coarse bins also make matched-fragment-COUNT voting
non-selective → the true peptide misses the top-64 → catastrophic recall.

## Key nuance: PXD-low-res is the worst case, NOT the target

The gate-blocking SPEED dataset is **Astral — high-res HCD, 0.02 Da bins (~25×
finer)**. There: small selective buckets, small `touched`, finer voting → the
degeneration may not occur and recall may be far higher. PXD low-res is the
pathological case. Astral (T5) is the real test of whether this approach is viable
at all.

## Options
1. Run Astral (high-res) — different dynamics; the actual speed target.
2. If Astral also degenerates → the simple vote-all-touched prefilter (Approach A)
   is the wrong algorithm; would need Approach B (Sage-style fragment-mass-sorted,
   precursor-filtered index) or shelve.
3. Restrict `--chimeric-frag-index` to high-res only (auto-off on low-res) so the
   PXD path keeps brute force.
