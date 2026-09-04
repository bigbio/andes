# Decision Trace — RS³ scoring system ABANDONED

**Date:** 2026-06-29
**Status:** ABANDONED (benchmark-refuted on all 3 gold datasets)
**Decision owner:** user; recorded by assistant
**Recoverability:** implementation discarded with branch `feat/rs3-scoring` (head `121d4d1c`, git-reflog-recoverable short-term). Research/design docs retained: `2026-06-29-rs3-spectral-significance-design.md`, `2026-06-29-unique-scoring-campaign-plan.md`, `2026-06-29-literature-review-*.md`.

---

## Decision
Stop pursuing RS³ (per-spectrum decoy-calibrated empirical-null significance, emitted as additive Percolator PIN features `Rs3NegLog10P`/`Rs3StdScore`). It does not clear the merge gate (more PSMs **or** faster than current andes) and is removed from the active line of work.

## Context (how we got here)
- andes leads the field on Astral (high-res) and TMT, trails Java MS-GF+ ~5% only on **low-res UPS1**.
- A multi-agent campaign (literature review + adversarial review + judge) identified **per-spectrum significance calibration** as the textbook low-res lever, and RS³ — a renewal/decoy-calibrated empirical null of andes's *own* emitted score — as a patent-free alternative to MS-GF+'s patented generating function (US 8,639,447).
- The adversarial review + judge **explicitly predicted RS³ might land flat**, for two reasons (see below). The user chose to skip the cheap Gate-0 proxy and test directly on the 3 gold datasets — the authoritative test.

## Evidence (the benchmark that decided it)
Single-variable A/B (baseline vs `--rs3`, identical params, bundled-store model auto-selection, Percolator q≤0.01, on `pride-linux-vm`, results in `/srv/data/msgf-bench/rs3ab/`). Rs3 columns verified populated in `--rs3` arms, exactly 0 in baseline → clean A/B, Percolator genuinely saw the features.

| Dataset | baseline PSMs@1% | `--rs3` PSMs@1% | Δ | wall |
|---|---|---|---|---|
| Astral (high-res) | 41,512 | 41,344 | **−168 (−0.40%)** | ≈ neutral (962→932s) |
| TMT | 11,407 | 11,504 | **+97 (+0.85%)** | ≈ neutral (179→178s) |
| UPS1 (low-res) | 17,476 | 17,469 | **−7 (−0.04%)** | ≈ neutral (154→152s) |

- Flat-to-slightly-negative on all three. The **low-res UPS1 target — the gap RS³ was designed to close — did not move**.
- Speed-neutral (the 256 decoy `score_psm` calls/spectrum are cheap vs the full search), so the "faster" arm is not met either.

## Why it failed (both causes were pre-flagged)
1. **Collinearity.** andes already ships per-spectrum calibration features (`TailorScore`, `RawScoreCal`, `ChanceMatchSurprise`). Percolator already extracts that signal; `Rs3NegLog10P` adds nothing orthogonal → flat. (The deferred "Phase-0 measurement" was meant to catch exactly this before building.)
2. **Closed-search-flat.** Project-standing result: "discriminators/separation features are flat on closed search; they only pay when candidate space is EXPANDED." RS³ is a separation feature on a fixed candidate space.

This re-confirms a pattern seen repeatedly: **RichIonLLR, strong-score-on-low-res, and now RS³ all landed flat/negative as closed-search scoring innovations.** andes's closed-search scoring is at/near ceiling.

## Alternatives considered and NOT taken (and why)
- **RS³ under expansion (`--rs3 --chimeric`)**: would distinguish "RS³ useless" from "RS³ needs expansion." User chose to abandon rather than run it — acceptable, since even if positive it would be a small effect gated behind chimeric, and the cleaner play is to invest expansion effort directly.
- **Collinearity ablation (drop TailorScore)**: would confirm redundancy vs no-signal; not worth the cycle given the flat 3-dataset result.

## Consequence / governing lesson
**Stop investing in new closed-search scoring functions.** The leverage for *more PSMs* is **coverage expansion** (semi-tryptic, open/mass-offset, deeper chimeric) and **model quality** (analyzer-/regime-matched retraining), not new per-PSM scores. The "faster" arm is already won, so a pure speed lever (fragment-ion index) only matters as a *substrate* that makes expansion tractable. Next options surveyed separately.
