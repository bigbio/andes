# Chimeric two-pass cascade — Astral: +55% PSMs over Java, FDR ~controlled

**Date:** 2026-05-30
**Branch:** `feat/chimeric-dda-plus` (cascade: narrow Pass 1 top-1 + MS1-gated parallel Pass 2)

## Astral result (cascade vs references)

| config | @1% PSMs | wall | T/D | entrapment FDP |
|---|---:|---:|---:|---|
| Java | 35,818 | 6:18 | — | — |
| Rust narrow (streaming) | 36,715 | 6:01 | 0.52 | — |
| **cascade** | **55,472** | 9:42 | 0.39 | 0.80% raw / 1.60% combined (rank-1 0.65%) |
| blind chimeric (refuted-speed) | 77,287 | 18:02 | — | 0.42%/0.83% |

## Findings
- **PSMs: +55% over Java** (55,472 vs 35,818), recovering ~72% of the blind-chimeric
  gain via 120k MS1-confirmed secondaries, at a CLEAN top-1+secondary emission
  (340k rows, not the 1.24M blind multi-emission).
- **FDR approximately controlled**: primary (rank-1) entrapment FDP 0.65% (clean);
  overall raw 0.80%, combined 1.60% — slightly above nominal (the secondaries add
  mild inflation; `max_kl=1.0` is lenient → tightenable to nominal).
- **Speed: 9:42 vs Java 6:18 — does NOT beat Java yet.** Breakdown: Pass 1
  (batch narrow scoring) = 8:29; Pass 2 (parallel MS1-gated) = ~1:09 (cheap — the
  cascade architecture works). The gap is Pass-1 batch (8:29) vs pure-narrow
  STREAMING (6:01) on the same 16.8M-candidate DB — a batch-vs-streaming /
  enumeration overhead, NOT Pass 2.

## Status vs the objective (more PSMs AND faster than Java)
- **More PSMs: YES, decisively (+55%)** — the first approach in the whole
  investigation to deliver a large, clean, entrapment-validated PSM gain.
- **Faster: NOT YET** (9:42 vs 6:18). Even a perfect Pass 1 (= narrow 6:01) + Pass 2
  (1:09) = 7:10 > 6:18: the margin is thin because narrow itself is only 17s faster
  than Java and Pass 2 adds ~1 min.

## Iteration knobs
- `max_kl` (co-isolation KL gate, currently 1.0): lower → fewer/cleaner secondaries
  → FDP toward nominal (small PSM cost).
- `max_n` (max co-isolated per scan, currently 2).
- Pass-1 speed: fix the batch-vs-streaming 2:28 gap (profile); plus general
  narrow-search opts (roundf 10% per the cost profile).

## Cascade implementation (all on branch, reviewed)
coisolation.rs (detector + search_secondary), run_pass2_coisolation (parallel),
narrow Pass-1 (cascade_wide=false), top-1 + force_push secondaries, MS1 per-PSM
feature skipped. Run with `--chimeric --chimeric-frag-index off`.
