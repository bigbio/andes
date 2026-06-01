# Chimeric fragment-evidence prefilter (fragment index) — design

**Date:** 2026-05-29
**Branch:** `feat/chimeric-dda-plus`
**Type:** Feature design (speed enabler; phased, bench-gated)
**Motivation:** `2026-05-29-gate-chimeric-norescore-vs-java.md` — chimeric NO_RESCORE
delivers real, entrapment-validated PSM gains (PXD +21.6%, Astral +116% vs Java) but
loses on **speed** (wide-window candidate explosion; Astral 17:04 vs Java 6:18). The
fragment index cuts the per-spectrum candidate set so chimeric can also win speed.

## Goal

Make `--chimeric` search **faster than Java on PXD001819 + Astral** while preserving
its entrapment-validated PSM gains. Speed target is modest (Astral ~2.7×: 17:04 →
< 6:18; PXD ~1.15×: 1:33 → < 1:21) — NOT the 10× the abandoned Java plan chased.

## Why the Java fragment-index abandonment does NOT apply here

`~/.claude/plans/msgfplus-fragment-index/ABANDONED-2026-04-20.md` abandoned the
fragment index, but for reasons that are inverted or absent in this context:
- **Root cause #1/#4 (Java can't match Sage's Rust data structures; rewrite-in-Java
  loses 2–5×):** we are now IN Rust. Inverted.
- **Tested on NARROW search** (~5 candidates/spectrum → no Tier-2 savings, Tier-1
  overhead dominates): our case is **wide-window chimeric** (thousands/spectrum) —
  exactly where fragment indexing wins. We scope to chimeric ONLY; narrow keeps the
  current fast bucket scan (already beats Java).
- **Memory (global `Map<SpecKey,Set<String>>` OOM):** v2/this design uses a transient
  per-spectrum vote buffer (no global map). Sidestepped.
- **Recall 95.3% < 99.5% gate:** that was a fingerprint-Hamming Tier-1; here we vote
  by real rank-weighted fragment evidence and gate recall vs brute-force chimeric.

## Scope

Fragment index activates ONLY under `--chimeric` (wide windows). Narrow / `--chimeric
off` is untouched and bit-identical. Approach **A** (in-loop prefilter): the index
changes *which* candidates are scored, never *how* — the validated scoring/emission/
FDR path is preserved by construction.

## Design

### §1 Architecture & data flow
Build a fragment index once per `PreparedSearch`, only when `params.chimeric`, over
the enumerated `candidates` slice (target+decoy, mod-expanded — already built). Per
chimeric spectrum in `run_chunk_inner`: keep the cheap mass-window enumeration, but
accumulate per-candidate rank-weighted fragment votes from observed peaks, keep the
**top-K** in-window candidates, and run the existing `score_psm`+edge+per-charge GF
on only those K (instead of all mass-window candidates).

### §2 Index structure & build
CSR layout: `bucket_offsets: Vec<u32>` (n_buckets+1) + `bucket_candidates: Vec<u32>`
(concatenated per-bucket candidate-id lists). Bucket = quantized charge-1 b/y
theoretical fragment mass at the scorer's matching tolerance (high-res 20 ppm /
low-res 0.5 Da → fixed bin width chosen accordingly; ±tol may map a peak to 1–2
adjacent bins at query time). Two-pass build (count fragments/bucket → prefix-sum →
fill). Memory ≈ total_fragments × 4 B (~1 GB at Astral scale). New module
`crates/search/src/fragment_index.rs`. Built only under `--chimeric` (gate).

### §3 Per-spectrum candidate generation (hot path)
Per-thread reusable scratch: `votes: Vec<f32>` (len = n_candidates) + `touched:
Vec<u32>`; reset is O(touched) (Sage pattern; no global allocation per spectrum).
For each kept (deconvolved, charge-1) observed peak: map m/z → fragment bucket(s) →
for each candidate in those buckets, add a **rank-weighted** vote (weight from the
peak's rank via the scorer, so strong peaks count more) and record `touched`.
Select **top-K** among touched candidates whose precursor mass ∈ the isolation
window (cheap mass check at selection). Feed top-K into the existing scoring path
unchanged. Reset votes via `touched`.

### §4 Recall & correctness gates
- `--chimeric off` / non-chimeric → index unbuilt, path not entered → PIN
  bit-identical (existing gate).
- Chimeric+index must reproduce **≥99.5% of the brute-force chimeric PSMs @1% FDR**,
  with T/D ratio preserved and **entrapment FDP unchanged** (re-run the
  `entrap-pxd` / `entrap-astral` harness with index on). Tune K, bin width, and
  rank-weighting to hit recall.
- CLI `--chimeric-frag-index {auto,on,off}` (default auto = on under `--chimeric`)
  enables direct A/B (brute-force vs index) for the recall gate and a safe revert.

### §5 Speed & memory gates
- Wall: chimeric+index **< Java** (Astral < 6:18, PXD < 1:21). Revert if not met.
- Memory: index ≤ ~1–2 GB; total RSS within the VM budget (~27 GB). Measure; if the
  index blows memory, fall back to bucket sharding by precursor-mass slab.

### §6 Testing (TDD)
1. Unit: hand-built 3-candidate index; a spectrum with peaks at candidate B's
   fragments → B has the top vote → top-K contains B.
2. Unit: in-window restriction — a high-vote candidate outside the mass window is
   excluded from top-K.
3. Unit: rank-weighting — higher-rank peaks contribute more vote.
4. Unit: `touched`-list reset leaves `votes` zeroed between spectra.
5. Integration: chimeric+index PSM count ≈ brute-force chimeric (recall gate) on a
   fixture (BSA `test.mgf`).
6. OFF bit-identity (existing parity test).

## Phase decomposition (each a gated, independently committable milestone)
- **P1** — `fragment_index.rs`: CSR index struct + build over `candidates` (unit tests).
- **P2** — per-spectrum vote/top-K generator (unit tests in isolation).
- **P3** — wire into `run_chunk_inner` under `--chimeric` behind `--chimeric-frag-index`;
  **PXD recall gate** (≥99.5% of brute-force; entrapment FDP unchanged).
- **P4** — Astral speed/memory tune → **< Java wall**, index ≤ budget.
- **P5** — TMT + cross-dataset; confirm entrapment FDP preserved; bench note + PR.

## Risks & mitigations
| Risk | Mitigation |
|---|---|
| Index memory blowup (the Java killer) | u32 CSR layout; measure P1; precursor-slab sharding fallback |
| Recall loss (top-K misses true peptide) | gate ≥99.5% vs brute-force; tune K/bin/weight; A/B flag |
| Perturbing validated FDR | scoring/emission untouched; re-run entrapment harness with index on |
| Narrow-path regression | index gated on `--chimeric`; off path bit-identical |
| Deconvolution/charge mismatch in voting | reuse the same charge-1 deconvolution the scorer/diagnostic use |

## Out of scope
- General (narrow) fragment-index replacement (the abandoned 10× goal) — narrow already beats Java.
- TMT PSM gap (Lever-2a / GF-shape) — orthogonal; chimeric doesn't help TMT.
- MS1/precursor changes; rescore (dead weight — NO_RESCORE is the chimeric path).

## References
- `2026-05-29-gate-chimeric-norescore-vs-java.md` (the speed blocker)
- `2026-05-29-entrapment-fdp-reversal.md` (the validated PSM gains to preserve)
- `~/.claude/plans/msgfplus-fragment-index/{ABANDONED-2026-04-20,speed-rewrite-v2,design}.md` (Java-era algorithm + lessons)
- Sage (lazear/sage) fragment index; MSFragger fragment-ion indexing.
