# Proof experiment: predicted-spectrum rescoring + entrapment-FDP on chimeric

**Date:** 2026-05-29
**Branch:** `feat/chimeric-dda-plus`
**Type:** Experiment (measurement; bundles plan Step 1 entrapment + Step 2 rescoring)
**Follows:** `2026-05-29-chimeric-full-review-and-rethink.md`

## Hypothesis

Coincidental real-DB-sequence targets win chimeric rank-1 by chance and
reversed-decoy TDC is blind to them (our finding). The field's discriminator is
**predicted-spectrum evidence** (MS2PIP predicted fragment-intensity correlation +
DeepLC ΔRT). **Prediction:** adding MS2PIP+DeepLC features before Percolator
collapses the **entrapment FDP** of chimeric IDs toward nominal 1%, while
Percolator-only leaves it inflated.

## Dataset

PXD001819 (LTQ-Orbitrap Velos, UPS1+yeast; wide ~2–3 Da windows). Smaller/faster
than Astral; entrapment-FDP logic holds regardless of real-vs-spurious co-isolation.

## Method

1. **Entrapment DB.** From the target FASTA, generate **shuffled-target entrapment**
   proteins (C-termini fixed, 1:1 target:entrapment; FDRBench-style). Final DB =
   target + entrapment + decoys (reverse of both). Track which accessions are
   entrapment.
2. **Re-search chimeric.** Run msgf-rust `--chimeric` (with `MSGF_CHIMERIC_NO_RESCORE=1`
   for clean Phase-1 emission) on the entrapment DB → chimeric PIN.
3. **Two arms (same PIN):**
   - **(A) baseline:** PIN → Percolator → accepted PSMs @1% q.
   - **(B) rescored:** PIN → **MS²Rescore** (MS2PIP intensity-correlation + DeepLC
     ΔRT features appended) → Percolator → accepted PSMs @1% q.
4. **FDP.** For each arm, FDP = (entrapment accepted) / (target+entrapment accepted),
   adjusted by the target:entrapment ratio (FDRBench combined/paired formula). Also
   report rank-1-only FDP (the coincidental-inflation locus).

## Success criteria

- **Primary:** FDP(B) ≪ FDP(A); ideally FDP(B) ≈ nominal 1% while FDP(A) ≫ 1%.
  → predicted-spectrum rescoring is the discriminator (thesis proven).
- **Secondary:** rank-1 count/FDP deflates under (B); MS2PIP correlation +
  DeepLC ΔRT carry real discriminative weight in the Percolator model.
- **Null result is also informative:** if FDP(B) ≈ FDP(A) and both inflated →
  rescoring is NOT sufficient → coincidental inflation needs a different control.

## Build order (each a gate; stop + report if a gate fails)

1. **Install MS²Rescore stack** (ms2rescore + ms2pip; DeepLC 2.2.38 present). Riskiest
   — do first. Pick MS2PIP model for PXD001819 fragmentation (CID/HCD low-res ion-trap
   → MS2PIP CID model; confirm from the mzML).
2. **Entrapment DB** builder + accession manifest.
3. **Re-search chimeric** on entrapment DB.
4. **Arm A** (Percolator-only) FDP.
5. **Arm B** (MS²Rescore → Percolator) FDP.
6. **FDP compute + verdict** note.

## Tooling notes / risks

- MS²Rescore reads PSMs via `psm_utils`; PIN is the `percolator` format. Must map
  SpecId → mzML scan for MS2PIP (observed spectra) + DeepLC (RT). Glue may be needed.
- MS2PIP model must match fragmentation; wrong model → meaningless correlation.
- Percolator: VM has docker 3.6.5; bench wrapper pulls 3.7.1. MS²Rescore can call its
  own Percolator — keep both arms on the SAME Percolator version for comparability.
- Entrapment ratio 1:1 keeps the search ~2× (target+entrapment) +decoys ≈ 4× DB; PXD
  is small so acceptable.

## Out of scope

- Astral arm (follow-up if PXD proof is positive).
- Productionizing rescoring into msgf-rust (this is an external-handoff measurement).
- Narrow-search gate run (separate, gate-relevant experiment).

## References

- MS²Rescore (compomics); MS2PIP; DeepLC.
- FDRBench (Noble-Lab) / Wen et al. Nat Methods 2025 — entrapment FDP.
- `2026-05-29-chimeric-full-review-and-rethink.md`.
