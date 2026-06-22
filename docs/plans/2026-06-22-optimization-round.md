# andes optimization round — 2026-06-22 (/loop)

**Goal:** next-round gains in models / scoring / chimeric / refinement-for-PTMs, beyond v0.2.0.
Metric: PSMs @ 1% TRUE entrapment-FDP (Percolator). One A/B per iteration, banked.
Hard rules: NO MS-GF+ ideas; Percolator for production FDR; experiment-hygiene (binary commit + data SHA + one variable + verify-don't-assume).

## State after v0.2.0
- Closed search beats Comet+Java+MSFragger (Astral +29-31% vs MF). Chimeric ~2x on Astral.
- Native rescorer = single-pass GBDT (validated; Percolator-class high-res). Two "improvements" (semi-sup GBDT, linear) REFUTED — don't re-try.

## Open problems / levers (ranked)
1. **★ Refine (PTM) FDR is UNVALIDATED** — Pass-2 is peptide-anchored, so the pooled q-value is blind to it (ENT_ flat while PSMs rise). MUST measure refine's SUBSET entrapment-FDP (split PIN by IsRefinement) before any refine "gain" is real. Then grouped/subset FDR (Percolator post-process by the IsRefinement group col) to validate.
2. **Chimeric headroom** — N co-isolated (default 4): sweep 4/6/8 at honest FDP; KL gate (0.3). Is there more recoverable co-fragmentation?
3. **Scoring** — strong vs rank per regime; RichIonLLR under expansion; charge-1 fragment blind spot.
4. **Models** — corpus expansion / stronger models (Codon); per-regime model variants.
5. **Search-space expansion** — semi-tryptic (--ntt 1) + group-FDR; open/delta-mass (needs frag-ion index).

## Iteration 1 (running): honest subset-FDP baseline on Astral
Run andes on Astral entrapment, CLOSED vs --chimeric vs --refine; rescore (Percolator); compute entrapment-FDP for the WHOLE set AND the subsets (Pass-1 only / chimeric-secondary only / refinement-Pass-2 only, via the IsRefinement + chimeric markers + SpecId row index). Tells us where the real headroom is + whether refine gains survive subset-FDP. -> picks the next iteration's lever.

## Results (banked)
(none yet)
