# Chimeric Sage-style fragment index (Approach B) — design

**Date:** 2026-05-30
**Branch:** `feat/chimeric-dda-plus`
**Type:** Feature design (speed enabler; phased, bench-gated; replaces the abandoned Approach A)
**Supersedes (for the ON path):** `2026-05-29-chimeric-fragment-index-prefilter-design.md` (Approach A)
**Motivation:** chimeric NO_RESCORE delivers real, entrapment-validated PSM gains
(PXD +21.6%, Astral +116% vs Java) but is ~3× slower (wide-window candidate
explosion). Approach A (in-loop vote-all-touched prefilter) was built, reviewed,
and **failed** — it touched the whole DB's fragment space per spectrum, so the
score set ballooned and the per-spectrum collect+sort dominated (PXD 18× slower,
Astral killed >75 min; `2026-05-30-frag-index-pxd-fails-lowres.md`). This design
ports Sage's candidate generator, whose dual-sort bounds the work to the precursor
window — the property Approach A lacked.

## Goal

Make `--chimeric` search faster than Java on Astral (wall < 6:18; current chimeric
~18:30) while preserving its entrapment-validated PSM gains (recall ≥99.5% of
brute, entrapment FDP unchanged). PXD low-res may stay on brute (auto-gated).

## The algorithm (Sage's `IndexedDatabase` / `IndexedQuery::page_search`)

Two structures (verified against lazear/sage `database.rs`):
- **Peptides sorted by precursor (monoisotopic) mass** → a precursor-mass window is
  a *contiguous index range* `[pre_lo, pre_hi]`.
- **Fragments sorted globally by fragment m/z**, divided into fixed buckets
  (Sage `bucket_size = 8192`), and *within each bucket re-sorted by peptide index*;
  a `min_value: Vec<f32>` holds each bucket's min fragment m/z.

Query = dual binary search: (1) precursor filter on the mass-sorted peptides →
`[pre_lo, pre_hi]`; (2) per observed peak, binary-search `min_value` to find buckets
overlapping the peak's m/z±tol, then within each bucket (peptide-index-sorted)
binary-search the `[pre_lo, pre_hi]` sub-range; fragments there within tol increment
that peptide's score. **Only peptides inside the precursor window are ever touched.**

## Design

### §1 Architecture & data flow
New module `crates/search/src/sage_index.rs`, built once per `PreparedSearch` when
`params.frag_index_active()`. Two structures over the enumerated `candidates`
(target+decoy, mod-expanded):
- `mass_order: Vec<u32>` — candidate ids sorted ascending by `candidate.peptide.mass()`.
  A candidate's position is its `pidx`.
- `fragments: Vec<Frag>` where `Frag { mz: f32, pidx: u32 }`, sorted globally by `mz`,
  chunked into buckets of `BUCKET = 8192`, each bucket re-sorted by `pidx`;
  `bucket_min_mz: Vec<f32>` (min mz per bucket).

Per chimeric spectrum, per charge in `charges_to_try`: derive the neutral-mass window
from the isolation window (reuse the existing chimeric window math), binary-search
`mass_order` (by candidate mass) for `[pre_lo, pre_hi]`. Build active peaks
`(mz)` from `scored_spec_for_charge(z).dump_active_peaks()`. For each peak: locate
candidate buckets via `bucket_min_mz` (a peak at `mz±tol` can span buckets — scan the
range), and within each bucket binary-search the `[pre_lo, pre_hi]` pidx sub-range;
for each `Frag` there with `|frag.mz - peak.mz| <= tol`, increment
`scores[pidx - pre_lo]`. Take **top-K** pidx by score → map `mass_order[pidx]` to
candidate ids → feed the **existing** per-charge GF scoring + emission unchanged.

### §2 Why this cannot degenerate (the Approach-A fix)
The per-spectrum score buffer is a **local `Vec<u16>` sized `(pre_hi - pre_lo)`**,
indexed by `pidx - pre_lo` — NOT a global candidate-sized array, and NOT a
touched-set that can grow to the whole DB. Scoring touches only peptides whose
precursor mass is in the window. Work per spectrum = `peaks × log(buckets) +
(matched fragments in window)`, bounded by the window size. No whole-DB sort.

### §3 Build & memory
`mass_order`: `argsort` candidates by mass. `fragments`: for each candidate at its
pidx, `predict_by_ions(peptide, 1..=1)` → push `Frag { mz, pidx }`; sort by mz;
chunk into 8192; re-sort each chunk by pidx; fill `bucket_min_mz`. Compute fragments
ONCE per candidate (avoid the 3× `predict_by_ions` the Approach-A review flagged).
Memory ≈ 8 B/fragment (f32 + u32) × ~280M (Astral) ≈ **2.2 GB** — within the ~27 GB
budget; report + measure at build. Built only under `--chimeric` (gate).

### §4 Correctness & gates
- `--chimeric off` / `--chimeric-frag-index off` → index unbuilt, path not entered →
  PIN bit-identical (existing gate).
- Chimeric+index recall **≥99.5%** of brute chimeric @1% FDR; T/D preserved;
  **entrapment FDP unchanged** (re-run the entrap harness with index on).
- Decisive: Astral wall **< Java 6:18**, index memory within budget.
- Auto-gate to high-res by default if low-res proves unhelpful (PXD keeps brute);
  the `--chimeric-frag-index {auto,on,off}` flag drives A/B + revert.

### §5 Testing (TDD) — with an EARLY degeneration guard
1. Unit (build): a candidate's fragment is retrievable at its pidx; `mass_order` is
   ascending by mass; `bucket_min_mz` monotonic non-decreasing.
2. Unit (query): a spectrum with peaks at candidate B's fragments AND B's precursor
   mass in the window returns B with the top score; a candidate with matching
   fragments but precursor mass OUTSIDE the window is excluded (the core Sage
   property); tolerance edge (a fragment just inside/outside tol).
3. **Local query-cost microbenchmark test:** build an index over a synthetic large
   candidate set, run the query against a synthetic dense spectrum (50+ peaks),
   assert per-query wall stays sub-millisecond. This catches any per-spectrum blowup
   in `cargo test` — BEFORE any VM cycle (the lesson from Approach A).
4. Integration: chimeric+index PSM count ≈ brute on BSA `test.mgf` (recall).
5. Off bit-identity (existing parity test).

### §6 Phase decomposition (each gated)
- **P1** `sage_index.rs`: `mass_order` + sorted/bucketed `fragments` build (unit tests, memory report).
- **P2** query (precursor range + fragment intersection + top-K) **+ the local query-cost microbenchmark** (unit tests).
- **P3** wire into `run_chunk_inner` under `--chimeric` behind the flag; off bit-identical; BSA smoke.
- **P4** PXD recall gate **+ local per-spectrum timing print** before any VM run.
- **P5** Astral speed/mem gate → **< Java wall**.
- **P6** TMT + cross-dataset; entrapment FDP preserved; bench note + PR.

## Risks & mitigations
| Risk | Mitigation |
|---|---|
| Per-spectrum query still too slow | window-bounded score buffer; P2 LOCAL microbenchmark gate catches it before VM |
| Build time (predict_by_ions over all candidates) | compute fragments once per candidate; measure at P1 |
| Memory (~2.2 GB Astral) | f32+u32 Frag; measure; precursor-slab partition fallback if over budget |
| Recall loss (top-K misses true peptide) | gate ≥99.5% vs brute; tune K; the precursor+fragment intersection is far more selective than Approach A's count-vote |
| Perturbing validated FDR | top-K feeds the UNCHANGED GF/emission; re-run entrapment harness |
| f32 m/z precision at 20 ppm high-res | f32 ~7 sig figs ≈ 0.0001 Da at 1000 m/z, finer than 20 ppm (0.02 Da); adequate (Sage uses f32) |

## Out of scope
- Narrow-search replacement (narrow already beats Java; abandonment lesson).
- TMT PSM gap (Lever-2a / GF-shape) — orthogonal.
- Replacing Sage wholesale / calling Sage as an external tool — this is a native port.
- The abandoned Approach-A `fragment_index.rs` stays on the branch as a record; the
  ON path moves to `sage_index.rs` (or `fragment_index.rs` is rewritten — implementer's
  call at plan time; keep the `--chimeric-frag-index` flag).

## References
- lazear/sage `crates/sage/src/database.rs` (`IndexedDatabase`, `IndexedQuery::page_search`, `bucket_size=8192`).
- `2026-05-30-frag-index-pxd-fails-lowres.md` (why Approach A failed — the bound this design adds).
- `2026-05-29-entrapment-fdp-reversal.md` (the validated PSM gains to preserve).
- `2026-05-29-gate-chimeric-norescore-vs-java.md` (the speed blocker + Java wall targets).
