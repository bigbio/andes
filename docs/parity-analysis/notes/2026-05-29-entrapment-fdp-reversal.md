# Entrapment-FDP REVERSES the chimeric verdict — we measured with a broken ruler

**Date:** 2026-05-29
**Branch:** `feat/chimeric-dda-plus`
**Design:** `superpowers/specs/2026-05-29-ms2rescore-entrapment-fdp-proof.md`
**TL;DR:** Built an entrapment database (shuffled targets, K/R fixed, 1:1) and
re-searched chimeric (NO_RESCORE Phase-1 emission) on it. The entrapment
false-discovery proportion (FDP) of the chimeric IDs is **~nominal on BOTH
datasets** — and it AGREES with the reversed-decoy q. The "chimeric inflates FDR"
conclusion (4 refuted levers) was an artifact of comparing chimeric counts to the
*narrow-search* baseline, not measuring true FDP.

## Result (Percolator-only, no rescoring; entrapment = known-false targets)

| Dataset | accepted @1% q | real | entrapment | entrapment FDP (raw / combined) | rank-1 FDP |
|---|---:|---:|---:|---:|---:|
| PXD001819 | 17,036 | 16,978 | 58 | 0.34% / 0.68% | 0.33% |
| Astral | 69,625 | 69,336 | 289 | 0.42% / 0.83% | 0.39% |

Reference counts: PXD narrow-off 14,808, Java 14,974; Astral narrow-off 36,715,
Java 36,271. Search wall (entrapment 2× DB, NO_RESCORE): PXD 22 s, Astral fast.

## The two rulers AGREE → the FDR is honest

At 1% reversed-decoy q we expect ~1% false. With a 1:1 target:entrapment DB,
false-targets ≈ false-entrapment, so combined entrapment-FDP ≈ 2× the raw
entrapment fraction. Astral: 0.83% combined ≈ the 1% decoy-q threshold. PXD: 0.68%
≈ 1%. **Two independent FDR estimators (reversed-decoy and entrapment) agree** —
the standard signal that FDR control is real. The accepted chimeric sets are
~99%+ genuine target peptides.

## Why this overturns the prior verdict

The four "refutations" (phase1 emission, phase2 MS1-KL, phase3 fragment-competition
+ residual rescore, rank-stratified FDR) all judged trustworthiness by **chimeric
count vs the narrow-search count** ("Astral 36,715 → 77,444 = +111% = inflation";
"rank-1 51,579 = +40% inflation"). **The narrow count is not ground truth — it is
just a narrower search.** The entrapment ruler (the field standard, Wen et al. Nat
Methods 2025 / FDRBench) shows the surplus is largely *real* co-isolated peptides,
not coincidental false positives. The "healthy decoy fraction but inflated count"
we flagged as coincidental-target-blindness is, by entrapment, simply **more real
IDs at controlled FDR** — exactly DDA+'s claimed behavior.

This is the broken-ruler failure mode the deep-research survey predicted
(`2026-05-29-chimeric-full-review-and-rethink.md`): we cannot judge wide-window
trustworthiness by count-vs-narrow; only an entrapment ground truth answers it.

## Consequence for the "need a predictor" thesis

The experiment was designed to PROVE predicted-spectrum rescoring (MS2PIP+DeepLC)
is required to make chimeric trustworthy. **It disproved that** for trust:
entrapment FDP is already ~nominal WITHOUT any rescoring. A predictor remains a
plausible *sensitivity* lever (more real IDs, à la MSBooster) but is **not needed
for FDR trust** here. Honest negative result on the stated hypothesis.

## Consequence for the merge gate (potentially large — needs verification)

PXD chimeric **17,036 > Java 14,974 (+14%)**, entrapment-validated; Astral 69,625
≫ Java 36,271. If these hold up under the verification below, chimeric is a *real*
sensitivity gain and the PXD gate blocker may be flipped. NOT yet actionable —
see caveats.

## Caveats / required verification before acting

1. **Entrapment realism.** If shuffled entrapment is easier to reject than real
   coincidental targets, FDP underestimates. Mitigant: it agrees with decoy-q.
   Verify: (a) entrapment vs real are ~50/50 among LOW-score (pre-FDR) target-label
   PSMs (fair competition); (b) try a foreign-species entrapment (e.g. add a
   non-sample proteome) as an independent check.
2. **Peptide-level, not PSM-level.** Counts above are PSMs (chimeric emits multiple
   per scan). Recompute distinct-peptide counts + peptide-level entrapment FDP for
   an honest vs-Java comparison.
3. **vs-Java fairness.** Chimeric searches a larger space; compare like-for-like
   (same FDP ruler on both; ideally entrapment-FDP on the Java run too).
4. **Speed.** Chimeric wall is higher (wide window). The gate needs speed AND PSMs;
   measure wall vs Java at the trustworthy config.
5. **Reproducibility / overlap.** Confirm chimeric IDs ⊇ narrow IDs + plausible
   extras (not a disjoint reshuffle); spot-check extra-ID peptide quality.

## Status

Reopens chimeric from "shelved/refuted" to "promising, pending verification." All
prior refutation notes remain valid as records of the *count-based* analysis but
their trustworthiness conclusions are **superseded by entrapment FDP**. Nothing
shipped; `--chimeric off` byte-identical.

## Artifacts
- Entrapment builder: `benchmark/parity/build_entrapment_db.py`
- FDP compute: `benchmark/parity/compute_entrapment_fdp.py`
- VM: `/srv/data/msgf-bench/entrap-pxd/`, `/srv/data/msgf-bench/entrap-astral/`
