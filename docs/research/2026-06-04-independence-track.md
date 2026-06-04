# CIMAS independence track — derivation audit + plan

**2026-06-04** — Goal: make CIMAS legally independent of MS-GF+ (no UC-licensed derived
code/models), so it can be relicensed (e.g. Apache-2.0). Keep `NOTICE` as the honest
acknowledgment of origin. No external consultation — replace the derived parts so
there is nothing to argue about. Final relicense gets an internal EBI-legal nod only.

## Derivation audit (what's derived, what's the work)

| Component | Files | Derived? | Replace =  | Effort |
|---|---|---|---|---|
| **Generating function** | `crates/scoring/src/gf/*` (~140 KB) | **Yes + patented** | RawScore ranking (already used) + a target-decoy/EVD calibration to replace `spec_e_value` | **MED** |
| **Core scoring** | `psm_score.rs`, `rank_scorer.rs`, `scored_spectrum.rs` (96 KB), `param_model.rs`, `fragment_ions.rs` | **Yes — line-by-line port** (61 Java-parity comments, copied constants `0.999497`/`chargeOrSeg`, `.param` loader) | clean-room reimplement the intensity-rank-LLR *idea* (Frank 2005 / Kim papers) against the `Param`/`score_psm` interface; drop bit-parity + `.param` loader | **HIGH** |
| **Trained models** | `resources/ionstat/models.parquet` (39 models) | **Yes — MS-GF+'s `.param` files** (slugs = exact lowercased filenames) | retrain all from public/MSnet data via existing `train` engine | **MED-HIGH** |

**Key favorable finding:** the GF does **not** drive candidate ranking — PSMs are
ranked by `rank_score`; the GF only attaches significance afterward (single producer
`compute_spec_e_values_for_spectrum`, `match_engine.rs:1001`). So removing it does
**not** change which peptides win; it removes only the significance score (Percolator +
a replacement calibration cover it).

## Coupling
- GF-removal ⟂ model-retraining (different artifacts).
- GF-removal **simplifies** clean-room scoring (lifts the integer-score-scale parity
  constraint on node scores).
- Clean-room scoring keeps the `Param`/`score_psm` interface → search/queue unchanged,
  retrained models drop in unchanged.

## Plan (dependency order)

### Phase 1 — Remove the generating function *(start here)*
- Add a GF-free path: skip `compute_spec_e_values_for_spectrum`; rank/sort by
  `rank_score` (no ranking change); replace `spec_e_value`/`e_value`/`de_novo_score`
  with a **patent-clean calibration** of `rank_score` (per-spectrum target-decoy
  empirical p-value or fitted EVD/Gumbel tail); handle the two other GF consumers (the
  mass-calibrator PSM gate → use a `rank_score` percentile; the final per-spectrum
  re-sort → by `rank_score`).
- PIN: drop `lnSpecEValue`/`lnDeltaSpecEValue`/`DeNovoScore` (or replace with the new
  calibrated stat); let Percolator lean on RawScore/DeltaRawScore + the rest.
- Ship as a flag first (bit-identical when off), measure, then make it the default for
  high-res.
- **Validation gate:** PSM/protein yield at 1% FDR on Astral (high-res — tonight's
  column-drop proxy showed −1.1%/−0.2%, validate end-to-end) + a low-res set. Caveat
  (from parity lessons): Rust raw TDC is weaker than Java, masked by Percolator — so
  validate the *end-to-end* search, not just a PIN column-drop. Keep only if high-res
  stays ~lossless; low-res CID is allowed to cost ~6% (recovered later / accepted).
- **Outcome:** patent excised; `gf/` subtree removable; biggest single derived chunk gone.

### Phase 2 — Clean-room the scoring code  *(swapped: was Phase 3)*
- Reimplement `rank_scorer` + `psm_score` + the `scored_spectrum` preprocessing from a
  spec written off the *public papers* (not the Java), keeping the `Param`/`score_psm`
  contract; drop the bit-parity goal and the `.param` binary loader. Phase 1 having
  removed the GF means the new node score is free of the integer-scale parity
  constraint.
- Largest single effort (`scored_spectrum.rs` is 96 KB of parity-tuned logic; the
  n=12 lesson warns scoring changes regress Percolator — so gate hard on yield).

### Phase 3 — Retrain models on public/MSnet data  *(swapped: was Phase 2)*
- Harvest high-confidence PSMs from PRIDE/MSnet reanalyses across
  activation×instrument×enzyme×label; train CIMAS-native models via the existing
  engine; write a fresh `models.parquet` with zero MS-GF+-derived rows; retire the
  `.param` fixtures. **Likely also carries the low-res 6%.**
- **Hard part:** bootstrap training *diluted* (−4.3% cross-dataset). Needs a better
  estimator (shrinkage/sharpening toward a non-derived prior) to match curated quality.
  This is the real research risk of the track.

### Phase 4 — Relicense (and drop the NOTICE requirement)
- Once Phases 2+3 are done, CIMAS contains **no MS-GF+-derived code or models** → it is
  no longer a derivative work, so the **UC license + the NOTICE attribution requirement
  no longer apply.** Switch `LICENSE` → Apache-2.0; the legal `NOTICE` can be dropped
  (keep at most a brief *courtesy* acknowledgment of MS-GF+'s intellectual influence —
  optional, good practice, not required).

### Phase R — Literature & strategy review  *(before rebrand finalization)*
- A full review of the papers + strategies (build on `internal-docs/papers`) to:
  (a) **make the analysis even faster** — fragment-ion indexing (MSFragger/Sage),
  tag prefilters, SIMD/bit-parallel matched-peak counting, cache-friendly layout,
  no-GPU tricks; and (b) **improve TMT PSMs** — predicted-intensity rescoring
  (MS2PIP/Prosit-TMT), TMT-aware / complementary-ion scoring, low-res-CID-specific
  methods, and what recovers the low-res 6%.
- Output: a ranked, evidence-backed action list feeding Phases 2/3 + the chimeric
  speed work. (Methodology: the same multi-agent literature fan-out used for the LLR
  charter.)

## Highest-leverage first move
**Phase 1 (remove the GF).** Patented + isolated + doesn't affect ranking + de-risks
Phase 3. Start now.

---

## Phase 1 — DONE + VALIDATED (2026-06-04, commit a8b16788)

`--gf-free` opt-in mode shipped (bit-identical when off, schema-parity tests pass).
End-to-end Percolator validation (same configs as the baselines):

| dataset | default (GF) | `--gf-free` | Δ PSMs / proteins | speed |
|---|---|---|---|---|
| Astral (high-res) | 37,176 / 4,745 | 36,750 / 4,735 | **−1.1% / −0.2%** | **186.9s vs 536.8s = 2.9× faster** |
| a05058 (low-res TMT) | 11,128 / 4,710 | 10,433 / 4,507 | −6.2% / −4.3% | 44.7s vs ~124s = 2.8× faster |

**Verdict: PASS for high-res.** Removing the patented generating function on high-res
is **~lossless AND ~3× faster** (the GF DP was a major runtime cost). Low-res CID keeps
the ~6% gap — recover via Phase 2/3 (own models + clean-room scoring) or keep GF for
low-res only until then.

**Next:** make GF-free the default by resolution (high-res → GF-free; low-res → GF until
recovered), then Phase 2/3 to delete the `gf/` subtree entirely. The speed win also
directly serves the "competitive with Sage/MSFragger" goal.

## Low-res 6% recovery — attempt 1: Tailor calibration (commit dec3bb03)

Added `TailorScore` (per-spectrum RawScore quantile calibration) as an additive PIN feature.

| dataset | GF baseline | GF-free | GF-free + Tailor |
|---|---|---|---|
| Astral (high-res) | 37,176 / 4,745 | 36,750 / 4,735 | **36,868 / 4,746** (≈ GF parity, +118 vs no-Tailor) |
| a05058 (low-res) | 11,128 / 4,710 | 10,433 / 4,507 | 10,417 / 4,503 (flat — NOT recovered) |

**Verdict:** Tailor brings **high-res to GF parity** (proteins ≥ GF, PSMs −0.8%, ~3× faster)
— high-res is effectively DONE, patent-clean. But it does **NOT** recover the low-res 6%:
the gap is not a calibration problem (Tailor normalizes only observed candidate scores;
the GF's low-res value was its full theoretical null over all peptides). Keep Tailor
(helps high-res, harmless low-res). Recover the low-res 6% via richer calibration (EVD)
or, more likely, Phase 2/3 (own models + clean-room scoring) — it looks like scoring
signal, not calibration.
