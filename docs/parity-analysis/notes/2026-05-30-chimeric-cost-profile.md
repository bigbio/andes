# Chimeric cost profile (perf) — candidate SCORING dominates, not emission/MS1/GF

**Date:** 2026-05-30
**Method:** rebuilt with `RUSTFLAGS=-C debuginfo=1`; `perf record -F 99` over ~4 min
of an Astral chimeric NO_RESCORE (brute, the 18:02 config); `perf report --sort=symbol`.

## Self-time (top)
| function | self% | bucket |
|---|---:|---|
| run_chunk_inner closure (inlined loop) | 34.0 | candidate loop |
| roundf | 10.3 | mass/peak matching (scoring) |
| ScoredSpectrum::edge_score | 10.0 | candidate scoring |
| score_psm | 9.7 | candidate scoring |
| psm_edge_score | 3.3 | candidate scoring |
| observed_node_mass | 3.0 | candidate scoring |
| directional_node_score_inner | 2.5 | candidate scoring |
| GF compute_inner | 3.6 | GF |
| compute_spec_e_values_for_spectrum | 2.3 | GF |

## Findings (overturns the prior hypothesis)
- **Per-candidate SCORING dominates (~60-70%)**: score_psm (node) + edge_score +
  psm_edge_score + observed_node_mass + directional_node_score + the roundf-heavy
  peak matching, plus the closure body (the candidate loop).
- **NOT emission**: `compute_psm_features` does not appear — the "9× emission rows"
  hypothesis is FALSE.
- **NOT MS1**: `precursor_isotope_match` does not appear.
- **GF is only ~6%**.
- Chimeric pays the scoring cost for THOUSANDS of wide-window candidates/scan, and
  the two-stage `could_win` gate (skips edge_score) is far less effective at
  top-N>1 than narrow top-1.

## Implication for the cascade design
The lever is reducing the NUMBER of candidates fully scored, WITHOUT dropping the
real co-isolated tail. The cascade does this via MS1 localization:
- Pass 1 narrow → primary peptide (few candidates, fast).
- Pass 2 → score only candidates within narrow tol of the MS1-detected co-isolated
  precursor mass(es), on the residual spectrum (a handful, not the whole window).
MS1 envelopes ARE the real co-isolated precursors, so this keeps the +116% (unlike
the fragment-vote prefilter that dropped the low-evidence tail).

Secondary micro-opt (flat speedup, both narrow+chimeric): roundf at 10% is high —
integer/cached mass math in the peak-match/node-score hot path could shave it.
