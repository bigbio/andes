# Phase 2 — Clean-room the scoring code (spec)

**2026-06-04** — Make SIMAS's scoring an **independent implementation** of the published
intensity-rank method, removing the derived *expression* of MS-GF+'s Java, so it is no
longer a copyright-derivative. Preserve scoring **semantics** (so the current models +
yield hold; Phase 3 swaps in own models later).

## What "clean-room" means here (and the honest limitation)

Copyright protects the specific **code expression**, *not*: the algorithm, the math, or
physical constants/facts. So we do **not** need to change *what* the scorer computes —
we need the *code, format, and comments* to be independently authored from the published
method + first principles, not a line-by-line translation of their Java.

- **Keep (free — facts/method):** the intensity-rank log-likelihood-ratio node scoring
  (Frank 2005; the academic basis), b/y fragment-ion physics, the nominal-mass integer
  scaling *value* (a measured mass ratio, derivable from amino-acid masses), the
  partition/segment *concept*, the per-rank probability *tables* come from the model
  (Phase 3 replaces the model data).
- **Replace (derived expression):** the MS-GF+ `.param` **binary format** + its loader,
  the line-by-line code structure that mirrors `FastScorer`/`DBScanScorer`/`getNodeScore`,
  the **58 Java-parity comments**, and constants documented as "MS-GF+'s X" (re-derive
  from first principles instead).

Honest limitation: this is *independent reimplementation from public sources*, not a
formal legal clean-room (which walls implementers off from ever seeing the original).
For SIMAS's purpose (not a verbatim derivative) that is the practical standard.

## The method, from public sources (what we reimplement)

Cite/source these, not the Java (`internal-docs/papers/` + acquire-list):
- **Intensity-rank node scoring** (Frank, *J. Proteome Res.* 2005/2009; PMC2738854):
  each fragment ion contributes `log( P(rank | ion) / P(rank | noise) )` where peaks are
  ranked by intensity and the per-rank probabilities are learned (the model). Missing
  expected ions contribute a negative "absent" term.
- **Spectrum preprocessing**: rank peaks by intensity; precursor-peak filtering; optional
  charge-deconvolution; `prob_peak` = peak density used in the existence term.
- **Segments/partitions**: condition the per-rank tables on (charge, parent-mass bin,
  m/z segment). The segment = `floor(peak_mz / parent_mass × num_segments)` — a generic
  binning, document as our own.
- **Node + edge (cleavage) scoring**: `score_psm` = Σ over cleavage sites of
  prefix+suffix node scores; edge/cleavage credit from the amino-acid set.

## Derivation map (current code → action)

| File | lines | parity comments | action |
|---|---|---|---|
| `fragment_ions.rs` | 262 | **0** | already generic (b/y physics) — just confirm/keep; minor doc pass |
| `rank_scorer.rs` | 320 | 5 | reimplement the rank-LLR lookup independently; document `chargeOrSeg` as our own normalization; strip parity comments |
| `psm_score.rs` | 480 | 11 | reimplement node/edge summation; remove the "mirroring Java `DBScanScorer.getScore`" structure + comments |
| `param_model.rs` | 1168 | 4 | **delete the `.param` BigEndian binary loader** (`load_from_bytes`/`load_from_file`) — the Parquet store is canonical; keep only the `Param`/`Partition`/`IonType` types loaded from Parquet; re-derive the nominal-mass scaler from first principles |
| `scored_spectrum.rs` | 2140 | **38** | the big one — reimplement preprocessing/deconvolution/rank-assignment/`prob_peak`/cached prefix-suffix score tables independently; remove `FastScorer` layout mirroring + the 38 parity comments; re-derive constants |

Constants: `INTEGER_MASS_SCALER = 0.999497` and `chargeOrSeg` are **values/methods** (kept),
but must be **documented from first principles** (the scaler = mean nominal/monoisotopic
mass ratio; chargeOrSeg = charge capped at segment count) — not as "MS-GF+'s constant".

## Reimplementation order (each step gated on yield)

Easiest/lowest-risk → hardest. **Build + bench after every step.**
1. **`param_model.rs` — drop the `.param` loader.** Parquet is already canonical; deleting
   the binary-format parser removes the clearest format-derivation with zero scoring
   change. (Update tests that used `.param` fixtures → use the Parquet store / Parquet
   fixtures.) Yield must be **byte-identical**.
2. **`fragment_ions.rs`** — confirm independence (already 0 parity), tidy docs.
3. **`rank_scorer.rs`** — reimplement the LLR lookup from the spec; independent code +
   comments. Yield must hold.
4. **`psm_score.rs`** — reimplement node/edge summation independently. Yield must hold.
5. **`scored_spectrum.rs`** — the 2140-line core. Reimplement section-by-section
   (preprocessing → deconvolution → rank assignment → prob_peak → cached score tables),
   each sub-step bench-gated. Highest regression risk (n=12 lesson).

## Validation gate (every step)

- `cargo build --workspace` + `cargo test --workspace` green.
- **Benchmark yield must hold** vs the current GF-free baselines:
  **Astral 36,868 / 4,746** and **a05058 10,417 / 4,503** (±~0.5% tolerance). A step that
  regresses yield is reverted and reworked. (This is the same hard gate the n=12 parity
  lessons demand for any scoring change.)
- Track derivation removal: `grep -icE "java|msgf|parity|FastScorer|DBScan|getNodeScore"`
  per file should trend to **0**; `.param` loader gone; `git grep "load_from_bytes\|\.param"`
  clean in code.

## Out of scope (later phases)
- Phase 3 replaces the model *data* (the per-rank tables) with own-trained models.
- A *better* scorer (beyond independence) — this phase preserves semantics, doesn't try
  to improve yield (that's Phase R/3 + the chimeric work).
