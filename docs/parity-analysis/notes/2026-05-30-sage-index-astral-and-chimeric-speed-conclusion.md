# Sage index Astral result + chimeric-speed conclusion: fundamentally incompatible

**Date:** 2026-05-30
**Branch:** `feat/chimeric-dda-plus`

## Astral gate (high-res HCD)

| arm | rows | wall | MaxRSS | @1% PSMs |
|---|---:|---:|---:|---:|
| brute (off) | 1,242,613 | 18:02 | 12 GB | 77,287 |
| SageIndex (on, TOP_K=128) | 980,140 | 12:03 | 17 GB | 35,836 |
| Java | — | 6:18 | — | 35,818 |

## Verdict: Approach B is refuted by a FUNDAMENTAL recall/speed tension

- **No degeneration** (the Approach-A failure mode is gone — local microbenchmark
  held; 12:03, not >75min). SageIndex is a real 33% speedup over brute.
- **But recall = 46%** (35,836 / 77,287), and the smoking gun: SageIndex @1% =
  35,836 ≈ **Java's 35,818**. The prefilter dropped EXACTLY the chimeric-specific
  gain, collapsing back to Java's primary-ID level.

**Root cause — a top-K fragment prefilter is fundamentally incompatible with
chimeric's goal.** Chimeric's value is the *secondary co-isolated* peptides, which
by construction match FEWER fragments (they share the spectrum). So they rank below
the top-K and are dropped first. The prefilter keeps the dominant peptide (≈ Java)
and discards the co-isolated tail (the +116% gain). Raising TOP_K to recover them
converges to brute's candidate count → brute's speed. You cannot prefilter for
speed AND keep the chimeric gain — they are the same candidates.

Additionally, even brute→12:03 doesn't reach Java's 6:18 because the **GF
SpecEValue scoring** (not candidate enumeration) is the dominant cost, and scoring
the wide candidate set is irreducible if the co-isolated peptides are the goal.

## Both fragment-index approaches now refuted
- **A** (vote-all-touched prefilter): degenerates (touches whole DB) — 28:24 PXD / killed Astral.
- **B** (Sage-style precursor-bounded prefilter): no degeneration, but the top-K
  prefilter drops the co-isolated tail that IS the chimeric gain (46% recall).

## Conclusion for the objective (more PSMs + faster than Java)
**Chimeric delivers real, entrapment-validated PSM gains (PXD +21.6%, Astral +116%)
but is irreducibly ~3× slower than Java, and fragment indexing cannot fix that
without discarding the gains.** So chimeric cannot clear the "more PSMs AND faster"
gate via candidate-generation speedups. The speed lever for chimeric would have to
attack the GF scoring cost itself (a different, harder problem), or chimeric ships
only where its PSM gain outweighs a speed loss (not this gate).

Implementation (B-T1..T3 SageIndex) kept on branch as a reviewed, correct,
non-degenerate record behind `--chimeric-frag-index`; `--chimeric off` bit-identical;
nothing shipped.

## Durable wins from this investigation (unchanged)
1. Chimeric PSM gains are REAL + entrapment-validated (overturned 4 prior count-based "refutations").
2. Broken-ruler lesson: judge FDR with entrapment FDP, not count-vs-narrow-baseline.
3. Chimeric speed is not achievable via fragment indexing (A degenerate, B recall/speed tension).
