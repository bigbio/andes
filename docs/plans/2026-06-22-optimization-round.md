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

### Iter 1 (2026-06-22) — refine Pass-2 subset entrapment-FDP on Astral — ★ REFINE VALIDATED
- closed 38,011 @0.79%; **refine 40,838 @0.62% (+2,827 / +7.4%)**; refine Pass-2 subset 4,020 @ **0.60%** (honest, peptide-remapped).
- ★ The naive Pass-2 FDP=0.00% is a MEASUREMENT ARTIFACT: Pass-2 rewrites accessions to BASEPEP_ namespace (refinement.rs:255) → ENT_ prefix severed → entrapment metric structurally BLIND to Pass-2 (0/74,466 IsRef=1 rows carry ENT_). Honest FDP via bare-peptide re-map vs entrapment FASTA = 0.60% (12 ENT / 4,020). This is the mechanism behind the old "refine unvalidated" note — it was a TOOLING blind spot, NOT an FDR violation.
- VERDICT: refine gains are REAL + honest. Next = OPTIMIZE refine (not grouped-FDR). Carry-forward: (1) entrapment harness MUST peptide-remap for Pass-2 (ENT_ prefix reads 0%); (2) Pass-2 uses BASEPEP self-decoy TDC. Binary 6e896710.

### Iter 2 (2026-06-22) — richer refine tier on Astral — PARTIAL (VPN outage blocked final read)
- default (5 mods, max2): 522,048 Pass-2 candidates -> 48,631 winners; refine 16.9s; total 312.9s.
- expanded (+Methyl/Dimethyl K,R + Trimethyl K, max3): 1,989,209 candidates (3.8x) -> 62,685 winners (+29%, +14,054); refine 27.8s; total 325.6s. Ran clean (exit 0), validated config, under 31GB.
- ★ DECISIVE NUMBER PENDING: the expanded-tier Pass-2 honest FDP (peptide-remap) + percolator total q<=0.01 — NOT yet read (Mac VPN/DNS outage hit at the final percolator step). Job done + persisted on VM. Retrieve: `ssh pride-linux-vm 'tail -32 /srv/data/msgf-bench/refine_iter2.out'` or re-run `python3 /srv/data/msgf-bench/refine_trace_fdp2.py refine-iter2 {default,expanded}`.
- PRELIM read (cautionary): +14,054 winners but Methyl(+14.016 near-isobaric) + Trimethyl(+42.047 ~ Acetyl +42.011) are FDR-inflation risks at max_mods 3 -> likely PARTIAL BLOAT (higher Pass-2 FDP). DON'T flip default until expanded Pass-2 FDP confirmed <=~1%. Provisional: keep 5-mod default; expanded = opt-in --refine-config for methylation studies. Binary a77e2aa3.
