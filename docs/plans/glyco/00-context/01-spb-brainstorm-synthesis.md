# SP-B / G2 brainstorm synthesis (2026-07-02) + three-agent review

The user chose to fix the **peptide-axis ranking** (SP-B/G2) via a learned,
regime-matched model trained on a **harvested multi-dataset corpus** with
**multi-engine consensus labels**, reusing the existing `protocol=` model-store
infrastructure. Three independent agents pressure-tested it. Verdicts below.

## Design decisions (user-approved)
1. **Training data:** harvest multiple public glyco datasets on Codon (not just the
   523 truth scans).
2. **Labels:** multi-engine consensus (peptide + glycan agreement).
3. **Model:** retrain the existing own `strong` peptide-conditioned spectral model
   for the glyco regime, registered under a new `protocol=NGlyco` selection key.

## Agent 1 — reuse seam: infra ~90% reusable
`andes train` ingests externally-labeled PSMs from a flat parquet;
`split_store_by_protocol` writes `protocol=NGlyco` generically. New work:
`Protocol::NGlyco` variant (model + CLI enums) + one `build_selection_key` arm;
**glycan-stripped backbone training rows** (the real modeling decision); a
truth-TSV→parquet converter script. (Detail in `00-current-state.md` §5.)

## Agent 2 — premise critique: retraining alone LIKELY-INSUFFICIENT
Stepped-HCD b/y is physically too sparse (a cross-spectrum glyco engine: ~11% direct-b/y ID;
+48% from cross-spectrum, not scoring). Recommendations:
- Kill-gate on **decoy-separated ranking**, not find-rate.
- **Add the Y0/Y1 peptide-mass-anchor feature** (the wasted discriminating signal).
- If the gate fails, the real fix is **cross-spectrum transfer (G4)**.

## Agent 3 — labeling/leakage: three corrections
- **Hold out ALL of PXD025455** (all files — same serum pools/instrument/prep).
- **Multi-engine eval truth** (the reference engine ∩ a glyco search engine/a commercial glyco engine) to break the
  the reference engine-both-sides loop (the dataset was originally a commercial glyco engine; our truth is a
  re-search).
- **Tiered labels, not pure consensus:** Tier A = consensus (gold/calibration);
  Tier B = single-engine, Y-ion/oxonium-confirmed *hard* cases. Consensus alone
  selects easy well-fragmented spectra and *excludes the hard sparse-b/y cases the
  model exists to rank*. Decouple label-trust (via glyco evidence) from difficulty.
- **Datasets (PRIDE-verified), mixed species for anti-leakage:** PXD005411 (a glyco search engine2
  mouse brain), PXD030670 (human saliva, Q-Exactive — closest match), PXD011239
  (serum haptoglobin, QE — same lab, train-only), PXD020254/PXD016175 (Lumos
  stepped-HCD), PXD057219 (venom, QE-HF). Canonical glycan-composition tuple for
  cross-engine matching; join on `(raw_file, scan)`; target ≥100–300k glycoPSMs.

## The reframe (net)
Retraining is **necessary calibration infrastructure but not the ranking fix
alone.** The ranking win comes from **the Y0/Y1 anchor feature + cross-spectrum
transfer**, with the regime-matched model underneath. Sequence to de-risk:
**cheap Y0/Y1 anchor + kill-gate FIRST**, then the harvest+retrain (the chosen
path) as calibration, then cross-spectrum. This is the plan codified in
`50-roadmap/`.
