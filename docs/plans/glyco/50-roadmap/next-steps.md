# andes glyco — next-steps plan + self-review (2026-07-03)

State: honest, competitive, faster. **259 glyco-PSMs @1% FDR** (> published
the reference engine-nglycan ~197), glyco phase **~63→37 min (~1.75×)**. All harness bugs
(top-1 collapse, enumerated-only, mass-preserving-decoy trap) fixed; Codex +
CodeRabbit reviews applied.

## Self-review — issues I see in the current code (beyond Codex/CodeRabbit)

| # | Issue | Severity | Note |
|---|---|---|---|
| SR-1 | **`sialic_consistency` is asymmetric**: penalizes "spectrum HAS NeuAc-oxonium but glycan claims none" (−obs), but does NOT penalize "glycan CLAIMS sialic yet spectrum shows no NeuAc-oxonium" (gives +0). A sialylated glycopeptide MUST show NeuAc oxonium in HCD, so claiming sialic without it should be penalized. | med | fold into GI-2 part 2 |
| SR-2 | **GI-2 features (Y0/Y1 anchor, sialic) are currently INERT for the reported FDR.** On the peptide-reversal-decoy axis a target and its decoy share the glycan → identical feature → Percolator can't use them. They only pay on the GLYCAN axis (2D-FDR). So the 259 came from the evidence prefilter + GI-3, NOT the new features. | high (opportunity) | GI-2 part 2 is what activates them |
| SR-3 | **Evidence prefilter is a per-scan recall tradeoff**: a true weak-core-Y DB backbone is skipped if the scan has ANY other evidence backbone (the non-dropping fallback only triggers when the WHOLE scan is evidence-free). Net-positive on this dataset (259>230), but could drop a true ID on another. | low | monitor; make the gate configurable |
| SR-4 | **Tie-break fix computes `glycan_y_intensity` for ALL winners**, not just those tied at the max b/y rank. Correct but ~O(winners) ladder calls/scan; could restrict to max-rank ties. | low (perf) | micro-opt |
| SR-5 | **`solve_backbone` recomputed per (charge×isotope)** though the Y-ladder bin voting is ~isotope-independent. The identified >2× speed lever. | med (perf) | see speed plan |

## GI-2 part 2 — EMPIRICAL RESULT (2026-07-03): current glycan decoy is NON-FUNCTIONAL

Ran the separate-axis 2D-FDR on the clean top-1-collapsed PIN with the EXISTING
glycan decoy (`ANDES_GLYCO_DECOY=1`, `glycandecoy_` rows, label −1):

```
peptide-axis only (qP<=0.01) : 259 glyco-PSMs, 101 backbone-correct
2D (qP<=0.01 AND qG<=0.01)   : 0 — ALL 259 fail the glycan axis
```

**The glycan axis wipes out every ID.** Root cause: the current glycan decoy
differs from its target in *only* `YLadderScore` (`glyco_pin.rs` deliberately
keeps `SialicConsistency`/`GlycanMass` IDENTICAL on the decoy row — they are
"peptide-axis features, independent of the glycan"). One feature, and it does
NOT separate → Percolator can push no target's `qG` below 0.01 → 0 pass. So the
glycan axis is not a precision filter (the earlier optimistic read); with this
decoy it is an all-or-nothing wipeout, and the glycan-axis machinery is
effectively dead. **This makes SR-1 the blocking prerequisite, not a nicety.**

## Plan A — MORE glycopeptides (identified levers, re-ranked by the above)

1. **GI-2 part 2 — glycan-axis 2D-FDR, REBUILT around an ISOBARIC-COMPOSITION
   decoy (highest leverage, now with a concrete blocker).** The glycan decoy must
   be a DIFFERENT composition of ~equal mass (isobaric), so that ALL
   composition-conditioned features genuinely differ target-vs-decoy:
   `GlycanMass` residual, `SialicConsistency` (stop forcing it identical on the
   decoy row — recompute it for the decoy composition), oxonium-derived signals,
   and `YLadderScore`. Only then does the glycan-axis Percolator have real signal
   to separate glycan-correct from glycan-wrong, and `qG<=0.01` keeps the
   glycan-confident subset instead of zero. Fold in SR-1 (penalize claimed-sialic
   with no NeuAc oxonium). Re-run the 2D-FDR; the target is a NON-zero
   glycan-confident subset of the 259 (correctness), and — if the composition
   features add peptide-axis signal too — more IDs. Without the isobaric decoy
   this whole lever is inert (proven above).
2. **Sequon-preserving glyco decoys (Codex Finding 1).** A decoy that keeps N-X-S/T
   density (shuffle preserving sequon positions) makes the PEPTIDE-axis FDR honest
   AND yields a matched decoy space → more IDs at trustworthy FDR. Replaces the
   external-FASTA requirement for the default path.
3. **Cross-spectrum transfer (G4), redone honestly.** RT-gated transfer with donors
   restricted to FDR-ACCEPTED top-1 GPSMs (Codex's earlier caveat). Now that the
   harness is honest, re-test the sialylated/short-peptide sparse-b/y stratum.
4. **Search-space expansion.** Missed cleavages (glycopeptides often carry them) +
   semi-tryptic; the 22%-covered-only-partly large-glycan tail (GI-3 follow-up).
5. **SP-B learned glyco fragment model.** Train on the harvested consensus corpus;
   the real peptide-axis discriminator for the sparse-b/y stratum.

## Plan B — BETTER speed (toward the full >2×)

1. **Factor `solve_backbone` across isotopes (SR-5).** Split into (a) charge-level
   bin VOTING (isotope-independent, compute once/charge) and (b) cheap per-isotope
   precursor GATE + candidate build. Expected ~another 1.2–1.4× (the bin voting +
   `filter_map` were ~17% after the complement fix, recomputed ~4×). Recall-sensitive
   → re-validate 259.
2. **Trim the new co-bottlenecks.** `gbdt_eval` (17%, the strong-score GBDT in
   `score_psm`) + `compute_psm_features` (16%). Options: a cheaper phase-1 pre-score
   (rank-only) to prune before the GBDT, reserving the strong score for survivors;
   confirm `compute_psm_features` runs once/scan (the phase-2 reduction).
3. **Reduce the isotope range for glyco** (e.g. `0..=2` vs `-1..=2`) IF the −1 offset
   yields ~no IDs — measure before shipping (recall).

## Sequencing
GI-2 part 2 (Plan A #1) is the single highest-leverage next step — it turns the
already-built Y-ladder + sialic features into IDs. Do it on the clean pipeline,
re-validate, then the sequon-preserving decoy (A #2, also fixes Codex Finding 1),
then the `solve_backbone` isotope factoring (B #1) for the full >2×.
