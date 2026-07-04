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
