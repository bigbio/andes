# andes glyco — implementation plan: review fixes + SP-B scoring

**Date:** 2026-07-01 · **Branch:** `glyco-phase1` (base of this plan: `ecf0d38b`)
**Inputs consolidated:** Codex adversarial-review (`e8e4c6ca..HEAD`), CodeRabbit review (same), and the SP-B prototype diagnosis on local nglycan labels.

## Where we are (verified this session)
- **Generation is near-ceiling.** Clean-truth backbone find-rate: Y-first cascade @top_k=50 = **80.1%** vs DB-branch ceiling **85.9%** vs de-novo 61.8%.
- **VM end-to-end (new code dc88b819):** peptide-level find-rate **29.4%** (154/523, up from 17.2%); **177 glyco-PSMs @1%FDR** (up from SP-A2's 0) BUT only **66 true** (true-FDP ≈63%) = 12.6% of the reference engine's 523.
- **Bottleneck = ranking + FDR, not generation** (quantified on 523 truth scans): findable 154 → top-1 correct **83** (ranking loses 71) → @1%FDR true 66 (+111 false wrongly pass).
- **Driver speed** ≈ 31 min/run even on the VM (phase-1 b/y-scores every backbone×candidate) — coupled to scoring (a real glyco score is what safely bounds candidates).

## Phase 0 — Correctness fixes from the two reviews (fold in first; same files as SP-B)
- **P0.1 [HIGH, Codex #2] Cross-isotope dedup keeps the wrong residual for novel glycans.** `glyco_search.rs` unions all isotope offsets then dedups DeNovo hits by backbone mass before scoring; offset-0 wins deterministically, so `glycan_mass_residual`/`CalcMass`/`PsmMatch::isotope_offset` can be off by ≥1 isotope on M+1/M+2 picks. Fix: include isotope offset (or the residual) in the dedup key for unannotated (Source::DeNovo) hits, or defer residual emission until the supported isotope is chosen by scoring. **TDD:** novel-glycan M+1 pick → emitted residual == precursor(corrected)−backbone.
- **P0.2 [HIGH, Codex #1] Restore a BOUNDED DB fallback for 0-core-Y spectra.** Data confirms the loss: cascade 80.1% vs DB-ceiling 85.9% = the 5.8% gap is exactly these spectra. Design (graceful degradation): quorum-2 → quorum-1 → **if still empty, `db_branch` candidates fed THROUGH phase-1 b/y ranking + `backbone_top_k`** (Codex's recommendation). Only fires for the rare truly-0-core-Y set (NOT weak-ladder), so it avoids the 40-min brute-force blowup that hit when it fired for every empty-solver spectrum. **TDD:** 0-core-Y known-glycan spectrum → true backbone recovered, candidate count bounded ≤ top_k.
- **P0.3 [MED, Codex #3] Narrow the quorum-1 rescue.** `solve_backbone_min(min_core_y=1)` currently accepts a lone peak as ANY of Y0..Y5 → spurious 1-rung backbones (and the prototype showed quorum-1 added ~nothing: 61.4 vs 61.8). Fix: restrict single-rung rescue to **Y0/Y1 only** (the diagnostic rungs) and/or require complement or known-glycan-residual support. **TDD:** a lone Y3-only noise peak must NOT yield a backbone.
- **P0.4 [MED, Codex #4 + CodeRabbit #3] Probe isotope fidelity + label.** `glyco_probe` calls `hybrid_candidates` once at offset 0; production sweeps `isotope_error_range`. Sweep the same range in the probe (and compute the DB ceiling over corrected neutrals). Update the de-novo baseline label (no longer "oxonium-gated"; it evaluates all matched truth scans).
- **P0.5 [MINOR, CodeRabbit #1/#2] Test robustness.** Novel-glycan test: use an empty glycan list (not `n_glycan_list()`) so the "no annotation" premise is DB-change-proof. Isotope test: assert the wrong-offset case yields NO `Source::Db` hit for the true backbone (suggestion patch provided).

## Phase 1 — SP-B ranking (attack the 71-scan ranking loss: 83 → toward 154)
- **P1.1 Populate the DEAD `YLadderScore` feature.** `y_ladder_intensity_score` is hardcoded `0.0` in `glyco_search.rs` (~L428) → the feature separates 0.00. Populate it with the backbone+glycan Y-ladder INTENSITY match (not just the count `CoreYHits` gives). First indicated, smallest change.
- **P1.2 Add glyco-discriminating features:** full-glycan Y-ladder (Y0→intact, stepping the assigned composition), intact-Y presence, oxonium-composition consistency. Keep them ADDITIVE (per parity-tuning lessons) so Percolator isn't destabilized.
- **P1.3 Evaluate:** top-1 accuracy on the 523 truth scans (target 83 → toward 154), via the local labels — no Codon needed.

## Phase 2 — SP-B FDR (attack the 111 false @1%: true-FDP 63% → ≤ few %)
- **P2.1 Glycan-axis decoy → 2D FDR.** Peptide-reversal decoys don't model wrong-glycan hits. Add a glycan-axis decoy (StrucGP precursor +20–30 Da re-score / a glyco search engine monosaccharide-mass shift / a cross-spectrum glyco engine fragment m/z shift) and control FDR on peptide and glycan axes independently. FDR stays Percolator-only (per the FDR-boundary rule); the glycan-axis decoy is input to it, not a new FDR engine.
- **P2.2 Evaluate:** @1% true-FDP on the 523 truth (overlap script already built).

## Phase 3 — Driver speed (now enabled by a real glyco score)
- Use the glyco score / core-Y evidence as a CHEAP pre-filter to prune phase-1 candidates before full b/y scoring; the P0.1 dedup fix also cuts the isotope×backbone multiplication. Target: iterable (<5 min) VM runs.

## Phase 4 — Scale (only if the Phase-1/2 prototype validates on local labels)
- Codon multi-dataset harvest (PXD005565 dual a glyco search engine+a commercial glyco engine, PXD030622, PXD025859) → train a general glyco model the andes way. Resource-heavy; gated on prototype success.

## Evaluation harness (reusable, already built)
- `glyco_probe` (backbone find-rate, residue-convention) · `findrate.py` (peptide-level) · true/false @1% feature-separation script · glycan-coverage diagnostic · VM `/srv/data/msgf-bench/glyco_bench/` staged (mzML+FASTA+truth) · Percolator via `run_percolator_docker.sh` + top-1/scan collapse.

## Sequencing & gates
1. Phase 0 (correctness) — TDD each; single VM regen to confirm no regression + measure the P0.2 generation gain (80.1→~85.9%).
2. Phase 1+2 (SP-B ranking+FDR) — batch, one VM regen, measure top-1 accuracy + true-FDP on the 523 local labels.
3. Gate before Phase 4: only harvest on Codon if the local prototype lifts top-1 accuracy AND true-FDP materially.
4. No PR until andes beats the reference engine on glyco IDs@1%FDR (standing gate). Milestone commits on `glyco-phase1`, single closing PR.
