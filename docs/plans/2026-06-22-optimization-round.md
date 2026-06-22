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

### Iter 2 FINAL — richer refine tier REFUTED
- default 40,838 total / 4,020 Pass-2 @0.60%; expanded 40,774 total / 3,586 Pass-2 @0.50%. Expanded = -64 total, -434 Pass-2 confident IDs DESPITE +29% raw winners. Bigger search (3.8x candidates: Methyl/Dimethyl/Trimethyl) -> FDR pressure -> FEWER confident. KEEP the 5-mod default. Expanded = opt-in --refine-config for methylation studies only. Refine track DONE (validated + tier tuned).

### Iter 3 — chimeric N=4 vs N=8 — N=4 CONFIRMED optimal
- N4: 53,068 @1.21% (17,870 sec); N8: 53,207 @1.24% (18,004 sec). +139 PSMs (+0.26%), 8 ENT -> FDP rises. Co-isolation saturates at ~4 on Astral's narrow isolation; 5th-8th slots = FDR-noise. KEEP N=4. Chimeric depth lever tapped. Binary f6805e2a.

### Iter 4 (scoping) — charge-1 fragment blind spot (audit flagged HIGHEST value scoring lever)
The cosine + GBDT intensity prediction only model charge-1 fragments; 2+ precursors carry charge-2 fragments unmodeled. Investigating the code to scope an extension to charge 1..=2.

### Iter 4 (scoping) — charge-1 blind spot ALREADY ADDRESSED
strong_score.rs + ion_features.rs already use predict_by_ions(1..=2); scored_spectrum deconvolves z>1 frags onto charge-1 axis (1760-1810). Only mod_site_features (PTM localization) still 1..=1 (narrow). The audit's "HIGHEST value" charge-1 lever was fixed post-audit. No headroom here. (Pattern this round: the easy parameter/scoring levers are already tuned — refine tier, chimeric N, charge-1 all confirmed optimal — pointing to MODELS/corpus as the next real frontier.)

### Iter 5 (running) — semi-tryptic search-space expansion
--enzyme-specificity semi vs fully on Astral entrapment: does the semi-tryptic expansion (N-term processing, signal peptides) add REAL IDs at honest FDP, or bloat?

### Iter 5 — semi-tryptic BLOCKED by in-RAM enumeration OOM (real engineering finding)
- fully baseline healthy: 38,011 @0.79%. semi OOM'd at ~31.6GB even at 1/32 DB scale -> memory NOT DB-proportional = runaway in-RAM candidate enumeration. Root: build_base_peptide_index/enumerate_candidates materialize the FULL candidate Vec before compacting to 20-byte records; semi emits ~30x candidates. mmap index streams SCORING not the BUILD -> doesn't help.
- VERDICT: semi (+non-specific/open) UNMEASURABLE on 31GB VM. Real optimization = make candidate enumeration STREAMING/out-of-core (emit->compact-record->external-sort, not collect full Vec). Unlocks a whole class of expansion search. Iter 6 = scope it.
