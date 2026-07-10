# SP-A2 Comparative Result — andes --glyco vs the reference glyco engine (PXD025455 HCC_pool_Late_Fc3_r1)

**Date:** 2026-06-30. **andes commit:** 5e922c36. **FASTA (identical both engines):** 9606-reviewed-contam-decoy.fasta (20,716 targets + DECOY_ decoys + contaminants). **FDR:** Percolator 3.7.1 (same container, seed 42), local Docker. **mzML:** HCC_pool_Late_Fc3_r1.

## Headline (raw, as-run)
| Engine | glyco-PSMs@1%FDR | distinct glycopeptides@1%FDR | Percolator mode |
|---|---|---|---|
| the reference glyco engine 4.2 (**labile mode**) | 3,217 | 1,634 | Concatenated |
| andes --glyco (bare-backbone baseline) | **0** | 0 | Concatenated (after top-1/scan collapse) |

## BUT the comparison is INVALID as run — two problems, both must be fixed

### Problem 1 — the reference engine baseline is not a clean N-glyco set (labile ≠ nglycan)
The standalone the reference engine-4.2 jar errored under `nglycan` mode with explicit `mass_offsets`, so the run used `labile_search_mode=labile`. Labile mode does NOT enforce the N-X-S/T sequon. Proof: the reference engine scan 7933's "glyco-PSM" is Haptoglobin (P00738) peptide `AVGDKLPECEADDGCPKPPEIAHGYVEHSVR` — **contains no Asn**, so it cannot carry an N-glycan; it is a peptide + coincidental delta mass, not a glycopeptide. An unknown but large fraction of the 3,217 are such non-glyco false assignments. **Action: re-run the reference engine in proper `nglycan` mode (via FragPipe or corrected params) for a valid baseline.**

### Problem 2 — andes baseline scoring yields 0 IDs at 1% FDR (real andes gap, independent of #1)
After collapsing the andes PIN to top-1 PSM/scan (TDC winner by RankScore → Concatenated mode), the target:decoy ratio is ~1.24:1 (4,452:3,593) and Percolator returns **0 at q≤0.01 AND q≤0.05**. Root cause: the only target/decoy-DISCRIMINATING feature is the backbone b/y match (RankScore), which is weak on glyco spectra (intensity is dominated by oxonium + Y-ladder + glycan peaks; backbone fragmentation is sparse). The glyco features in the PIN (OxoniumScore, YLadderScore, CoreYHits, GlycanMass) are **backbone/spectrum-level — identical for the target and decoy peptide competing at the same backbone mass — so they add NO target/decoy discrimination.** Percolator therefore has only the weak backbone score to separate target from decoy → 0.

### Diagnostic confound
A find-rate check (does andes's candidate set contain the reference engine's peptide?) returned 0.1% — but this is confounded by Problem 1 (the reference engine's peptides are partly non-glyco, which andes correctly does not generate), so it overstates andes's miss rate. The reliable andes-internal signal is the ~1.2:1 top-1 target:decoy ratio.

### Secondary observations
- andes emits ~20 candidates/scan (precursor − each of 2,510 glycans × sequon peptides in each backbone window) → huge coincidental (peptide,glycan) candidate space; the true pair is one needle among many. This is the combinatorial false-match problem; candidate ranking is weak without a glyco-aware score.
- andes run did not set `--max-missed-cleavages`; glycopeptides often carry missed cleavages — verify andes default ≥2 and matches the reference engine's 2 in any re-run.

## Implications for the build
- **SP-B (learned glyco fragment model) is NOT optional** — it is the lever that makes the backbone b/y match discriminate on glyco spectra. The current baseline confirms backbone-only scoring is insufficient.
- A target/decoy-discriminating glyco score is needed (the spectrum-level oxonium/Y-ladder features cannot separate competing peptides at one backbone mass).
- The comparison must be redone with the reference engine in `nglycan` mode AND andes with SP-B scoring before any conclusion about parity.

**Gate verdict: INCONCLUSIVE.** Engine runs end-to-end; baseline scoring insufficient (0@1%FDR); the reference engine baseline mis-configured. No PR until a valid nglycan-mode comparison with andes learned scoring.
