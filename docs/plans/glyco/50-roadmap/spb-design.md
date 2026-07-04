# SP-B — glyco-specific rank model (design + staged plan, 2026-07-04)

## Why (diagnosis, from the collapse A/B)
Generation finds ~80% of truth backbones; only ~19% rank #1. The top-1-per-scan
collapse (required for honest TDC FDR) picks the winner by `rank_score` — andes'
`RankScorer` b/y score. That scorer is trained on regular high-res DDA, where
peptide backbone b/y are strong. In glyco spectra the precursor energy goes into
glycan cleavage, so backbone b/y are SUPPRESSED — a regime mismatch. The
collapse therefore mis-ranks on noisy b/y and discards correct backbones BEFORE
Percolator. A naive fix (make the core-Y ladder the primary collapse key) was
A/B-refuted (260→197). The real fix: a `rank_score` that models glyco b/y
statistics → better collapse → more correct backbones survive → more IDs.

## Architecture (reverse-engineered)
```
glyco collapse ── rank_score ◄── score_psm(ScoredSpectrum, peptide, RankScorer)
RankScorer ◄── RankScorer::new(&Param)
Param      ◄── model-train::estimate (per-fragment rank-LLR frequency vectors)
estimate   ◄── CountStats ◄── StatsAccumulator.accumulate(spectrum, LabeledMatch)
LabeledMatch { spectrum_index, peptide, charge }   ← THE input we must supply
store      ── Parquet model store; `andes --model-id <id>` loads a specific model
```
`andes train-from-search` already wires estimate→store, but labels via a STANDARD
search (`bootstrap_labels`) — not glyco. The single new piece is a path that
produces `Vec<LabeledMatch>` from GLYCO backbone PSMs and runs the same
accumulate→estimate→store. Everything downstream (RankScorer, --model-id) reuses.

## Staged plan
- **Stage 0 — feasibility gate (leaky, fast falsification).** Label from andes'
  OWN glyco top-1 @train-fdr PSMs on Fc3_r1 (self-training) OR the on-VM
  the reference engine Fc3_r1 IDs (glyco_cmp, 8210 PSMs). Train a glyco rank model,
  re-search Fc3_r1 with it (`--model-id`), measure backbone-correct vs the 101
  baseline. Train==test run ⇒ OPTIMISTIC, not a shippable number — but if a
  glyco-trained rank model can't beat the DDA model even in-domain, the
  direction is dead and we skip the expensive honest data. If it wins, proceed.
- **Stage 1 — honest cross-dataset.** Train on PXD005411 (a glyco search engine2 mouse-brain,
  SAME stepped-HCD regime, 17,188 public PSMs; FTP in 40-data/collection/
  README). Test on PXD025455 Fc3_r1. Clean (independent dataset/organism, same
  fragmentation). This is the shippable number.
- **Stage 2 — integrate.** Default the glyco driver to the glyco model (selection
  key / bundled store entry); keep the DDA model for non-glyco.

## Guardrails
- FDR stays Percolator-only on the 1-PSM/scan PIN (unchanged). SP-B only changes
  `rank_score`, i.e. WHICH backbone the collapse keeps — not the FDR.
- One variable per A/B; label the leaky vs honest numbers explicitly.
- TDD the new labels→accumulate path; reuse estimate/store as-is.

## Success gate
Stage 1 backbone-correct > 101 (and total @1% FDR ≥ 260) on Fc3_r1.

## STAGE 0 RESULT (2026-07-04) — FAILED; direction refuted
Trained a glyco rank model on the on-VM the reference engine Fc3_r1 backbone IDs
(8136/8136 labels parsed, 14 partitions), re-searched Fc3_r1 with `--model`:

| model | @1% FDR | truth scans | backbone-correct |
|---|---|---|---|
| DDA (baseline) | 260 | 237 | 101 |
| glyco-trained (leaky) | 239 | 218 | 94 |

**Worse — even leaky (train==test run).** So retraining the b/y RANK model does
NOT fix the ranking bottleneck. The regime-mismatch hypothesis is refuted: the
issue isn't a miscalibrated DDA model, it's that glyco b/y ions are
fundamentally SPARSE — training on them yields NOISIER per-rank statistics, not
better ones. Gate failed ⇒ stage 1 (PXD005411 download+train) is NOT worth it.

**Combined with the collapse A/B**, two single-signal collapse fixes are now
both dead: Y-ladder-primary (260→197) and glyco b/y rank model (260→239). The
"which backbone" signal is split across b/y + Y-ladder + oxonium; no single one
ranks. Remaining principled option: a LEARNED MULTI-FEATURE collapse score over
ALL glyco features (the 2-pass Percolator re-collapse — use pass-1 Percolator's
weights to re-pick the collapse winner from the top-K accepted candidates). The
`--labels` training path + collapse_cmp refactor are kept (reusable); the SP-B
rank-model arc is closed.

## CORRECTION (2026-07-04) — the stage-0 "refutation" was CONFOUNDED; WITHDRAWN
User pushback ("this has been mainly errors") + a Codex adversarial review +
self-audit found the stage-0 A/B changed TWO variables: the glyco model had
retrained frequencies AND an OWN-DERIVED geometry (train-from-search default),
which COLLAPSED to 14 partitions on the sparse glyco-only corpus vs the seed's
148. Re-ran with `ANDES_SEED_GEOMETRY=1` (frequencies-only):

| model (leaky, the reference engine Fc3_r1 labels) | @1% FDR | backbone-correct | partitions |
|---|---|---|---|
| baseline DDA | 260 | 101 | bundled |
| stage-0 glyco (own geometry) | 239 | 94 | 14 |
| glyco (seed geometry) | 255 | 93 | 148 |

Geometry fixed the COUNT (239→255 ≈ baseline) but backbone-correct stayed ~93
(< 101). So glyco frequencies are NOT dramatically worse (that was a geometry
artifact) but do NOT beat baseline on ranking either — inconclusive-leaning-
negative, NOT a clean refutation. **The stage-0 conclusion is WITHDRAWN.**

Still-uncontrolled confounds (Codex): single-engine (the reference engine) labels,
post-Percolator @1% metric (not decoy-separated ranking separation), same-run
label/eval overlap. A trustworthy call needs multi-dataset eval + multi-tool
consensus truth + a ranking-separation metric. ⇒ Do the dataset harvest
(PXD005411 a glyco search engine2, PXD016175 a glyco search engine2, PXD030670 a commercial glyco engine) FIRST, then re-test.
Do NOT close the direction on current evidence.

## HONEST RE-TEST (2026-07-04) — SP-B rank model refuted, controls CLEAN
Redid the test with every stage-0 confound controlled: seed geometry
(`ANDES_SEED_GEOMETRY=1` → 148 partitions, no degeneration), held-out labels
(the reference engine 8011 labels EXCLUDING the 196 consensus test scans → no test leakage),
measured vs BOTH truths.

| model | @1% FDR | vs 523 correct | vs 196 correct |
|---|---|---|---|
| baseline DDA | 261 | 101 | 51 |
| honest glyco rank model | 235 | 83 | 42 |

Worse on EVERY metric, cleanly. **The SP-B rank-model direction is now honestly
refuted** (the earlier confounded 239 is superseded; this is the trustworthy
verdict). Reason: the DDA rank model carries RICH b/y intensity-rank statistics
from non-glyco spectra; glyco b/y are SPARSE, so training the rank model only on
glyco data yields NOISIER per-rank vectors, not better ones. You cannot learn
better b/y ranking from a signal that is barely present. The discriminative
"which backbone" information is in Y-ladder + oxonium + multi-feature — which a
single b/y rank model structurally cannot capture. ⇒ Only remaining ranking
lever = the 2-pass Percolator re-collapse (learned MULTI-feature collapse score).
The `--labels` path + an open-source glyco engine oracle stay as reusable infra.

Review-found code bugs (fix before trusting any --labels result):
- label ingestion has no duplicate-scan/charge/conflict guard (didn't bite here —
  the reference engine export was 1-per-scan — but poisons on rank-alternative exports)
- collapse_cmp is applied AFTER a rank-only truncation in glyco_search, so the
  true ladder winner can be truncated before the shared comparator sees it
- the fixed-mod decorator is residue-wide only (fine for Cam-C; wrong for TMT)
