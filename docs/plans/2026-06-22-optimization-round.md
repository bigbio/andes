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

### Iter 6 — model frontier: bundle WINS; gain needs ANALYZER-MATCHED corpus
- shipped bundle 37,521@0.82% (WINS) vs own_winning v2-retrain 36,539@0.90% (-2.6%, MORE QExactive rows REGRESSED) vs own_release 21,656 (store-load issue). Bundled hcd_qexactive_tryp at corpus ceiling for its regime.
- ★ KEY: Astral = Orbitrap-ASTRAL analyzer; the QExactive-trained model only approximates it. The real gain = a DEDICATED `hcd_astral_tryp` slug trained on NATIVE Astral DDA tryptic spectra (analyzer-matched), NOT more QExactive volume. Iter 7 = discover Astral DDA tryptic PRIDE datasets to build that corpus.

## Round summary (so far)
Engine is WELL-TUNED: parameter (chimeric N), scoring (charge-1, strong/rank), and same-regime model levers are at/near ceiling — strong endorsement of v0.2.0. The ONE win: refine VALIDATED (+2,827 real PTM IDs @0.60% FDP; the 'unvalidated' was a BASEPEP measurement artifact). Two concrete next-gain frontiers: (1) dedicated ANALYZER-MATCHED Astral model (new corpus, Codon long-pole); (2) streaming candidate enumeration (unblocks semi/non-specific/open — currently OOMs in-RAM).

### Streaming-index fix (7e1e0e6f) — NECESSARY but INSUFFICIENT for full-DB semi
- Dropped the Vec<Candidate> (real memory win, parity-clean) BUT build_base_peptide_index still sorts the Vec<IndexRecord> IN RAM. Full-DB semi-tryptic Astral = >1.5B records (>31GB) -> still OOMs at ~31.5GB (dmesg-confirmed 3 OOM-kills). The 20-byte-records-fit estimate was wrong (records themselves >31GB).
- Full fix = EXTERNAL/chunked sort (spill records to disk + k-way merge by mass_milli) instead of in-RAM sort_by_key — or DB sharding. The streaming fix stands as a memory improvement (helps giant-mod-space too) but doesn't unblock full-DB semi alone.
- ★ DISCIPLINE: before the external-sort investment, TEST whether semi-tryptic even ADDS IDs (most expansion levers bloated this round). The fix now lets REDUCED-DB semi fit -> A/B semi vs fully on a reduced Astral DB to get the semi value signal cheaply. If semi wins -> external sort worth it; if bloats -> skip it.

### Semi-tryptic — ABANDONED (infeasible + untestable on 31GB)
- Reduced-DB A/B: fully completes (1/4 8,693; 1/8 4,784; 1/16 2,459) but SEMI OOMs at 31.6GB at EVERY fraction down to 1/16 (1994 prot), in the SEARCH LOOP (search_index_build 0.02s then SIGKILL) — a per-spectrum semi-candidate combinatorial blowup that mmap (mod-only deferral) doesn't bound. TWO OOM sites: full-DB index build (>31GB records, in-RAM sort) AND per-spectrum search enumeration.
- The streaming fix (7e1e0e6f) addresses only the build's Candidate Vec; it stands as a small parity-clean memory win for all searches but does NOT unblock semi.
- VERDICT: semi value UNTESTABLE on 31GB; unblocking needs bounded per-spectrum enumeration + external-sort + a bigmem box just to measure — not worth it speculatively (expansion levers bloated this round). ABANDON. Keep the streaming fix.

## ROUND CLOSED. Net: 1 validated win (refine, +2,827 real PTM IDs @0.60% FDP) + 1 small memory fix (streaming index). Everything else confirmed at-ceiling or infeasible -> strong endorsement that v0.2.0 is well-tuned. The ONE remaining real-gain frontier = a dedicated ANALYZER-MATCHED Astral model (hcd_astral_tryp, native Astral DDA corpus) — a fresh multi-session sub-project.

## DEDICATED ASTRAL MODEL CAMPAIGN (2026-06-23, user-driven corpus)
Frontier after the cheap-lever round: the QExactive-trained bundle is corpus-bounded for the Astral regime (more QExactive data REGRESSED -2.6%) -> train analyzer-matched Astral models. Two corpora, two phases:

### Phase 1 — BASE general model (hcd_astral_tryp)
- Corpus = MSV000098998 / PXD067958 (GSK HCP-on-Astral benchmark). Orbitrap Astral, CHO (Cricetulus griseus, no-leakage vs human HeLa benchmark), LABEL-FREE DDA Top80. 21 DDA runs / 29.94 GB (7 load levels x 3 reps), sample-verified genuine Astral raw. DIA (102 GB) excluded. SIL 13C/15N heavy spike-in in every run -> declare as variable mods. Enzyme assumed trypsin (no params in deposit). FTP blocked -> MassIVE HTTPS DownloadResultFile works.
- Harvest LAUNCHED (agent aa21980b) -> $B/astral_corpus/ + CHO FASTA (UP000001075). Next: search -> flat -> andes train -> A/B vs bundled on human Astral benchmark @ honest entrapment-FDP.

### Phase 2 — PTM-AWARE Astral model (the "future", user-endorsed 2026-06-23)
- Corpus = PXD065579 (Kumar et al, MCP 2026:101562). Human HCT116/Jurkat, PTMScan immunoaffinity enrichment (phospho-Y, acetyl-K, ubiquitin K-e-GG, methyl, succinyl), Orbitrap Astral + Lumos, nDIA + DDA mixed. Use ONLY the Astral+DDA subset (inventory agent ae88134e scoping it).
- TIES TO MaxSBM/glyco doc (docs/plans/2026-06-20-glyco-neutral-loss-and-maxsbm.md): train diagnostic/loss-ion (d-ion/p-ion) behavior on REAL Astral PTM spectra — phospho neutral losses, ubiquitin-GG remnant, acetyl. pepXML/mzML in the deposit may give ready PTM labels without re-searching.
- SECONDARY for the base model too (its Astral-DDA runs add analyzer diversity), but PTM-enrichment-biased so not the base.

### Phase 2 inventory result (PXD065579 Astral+DDA subset) — scoped 2026-06-23
~25 Astral-DDA .raw runs (EXCLUDE ~20 Astral-nDIA + ~20 Lumos-DDA). By PTM:
- Ubiquitin K-e-GG: 8 (HCT116, incl MG132/Control) — RICHEST + the +114.0429 K-e-GG remnant IS the MaxSBM d-ion/p-ion case -> BEST PTM to start.
- Phospho-Y: 6 + IMAC global pSTY: 3 = 9 phospho total (classic -98 HPO3 / -80 neutral loss = andes loss_class already handles).
- Acetyl-K: 4. Succinyl-K: 2 (mouse liver). Methyl-R: 2 (mouse liver).
★ READY-MADE LABELS (no re-search): FragPipe interact-*-Astral-DDA-*.mod.pep.xml (PeptideProphet PSMs w/ localized mods) + *_calibrated.mzML + combined_site_{K_114.0429,STY_79.9663,K_42.0106}.tsv -> pair calibrated spectra <-> localized mods directly to build PTM training flats. FTP ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/05/PXD065579/ (.raw not on API page1 -> FTP-list to size).
CAVEAT: thin per-PTM (pilot), single-rep input titrations -> good for diagnostic-ion PRESENCE learning, weak for generality; supplement before claiming generality. Start with ubiquitin KGG.

### PXD065579 COMPLETE FTP inventory (ftp.pride.ebi.ac.uk/.../PXD065579/, 2026-06-23)
679 files / ~733 GB. Astral-DDA .raw = 25 / 139.5 GB (bigger than assumed: KGG 8 [3.5-7.9G], pY 6 [1.7-5.2G], AcK 4 [2.3-4.7G], IMAC-phospho 3 [9.7-18G!], mLiver SuccK 2 + R-Me 2 [4.6-9.8G]). EXCLUDE: Astral-nDIA 120/248GB, Lumos-DDA 164/35GB. Ready labels: Astral-DDA calibrated .mzml (~304GB all) + .pepxml/.mod.pep.xml (66/7.5GB) + site TSVs. FTP works directly (local Mac has egress to PRIDE HTTPS).
DECISION (user asked which to start): START WITH MASSIVE (MSV000098998) for Phase-1 BASE — clean label-free, 30GB (vs 139GB raws), right population for the general Astral benchmark; search step is routine andes. PXD065579 is EASIER (pre-made mod.pep.xml labels + working FTP) but PTM-biased -> Phase-2 only (its ready labels are the head-start there). Lean Phase-2 pull = calibrated mzML + mod.pep.xml (skip the 139GB raws).

### Campaign progress (2026-06-23)
- PRIDE discovery (ad7f4343): ★ BENCHMARK SOURCE = PXD070049 (ProteoBench LFQ Astral DDA 15min/50ng hybrid, Olsen) + DIA sibling PXD071205 -> EXCLUDE from training (leakage). PRIDE has 295 Astral projects but only 51 DDA (Astral is mostly DIA); clean whole-proteome label-free human Astral DDA is SCARCE. BEST base = PXD046453 (Guzman/Olsen NBT2024): HeLa 200ng Astral DDA ~12 runs, label-free trypsin+LysC, SHIPS DDA IDs (top100_DDA.zip/Cycle_DDA.zip). Caveat: SAME LAB as benchmark (not leakage, but mild lab-style optimism) -> use CHO (diff lab+species, already harvested) as INDEPENDENT cross-check. Other: PXD069898 (SpecFormer, diff lab, purpose-built intensity corpus, HeLa+mouse, .sepr2 needs re-search), PXD071864 (insect, non-human diversity).
- DESIGN: train base on PXD046453 (regime-match, easy ready IDs) AND CHO (independent) -> A/B both on benchmark; both-improve = real not lab-fit.
- PXD046453 harvest+label-inspect LAUNCHED (aa34922a). Base CHO harvest job 54144438 still downloading.
- Phase-2 KGG head-start DONE (a53dd208): 8/8 ubiquitin K-e-GG raws (38GB) + site truth combined_site_K_114.0429.tsv + FragPipe labels (6/8 runs) at $B/astral_corpus_ptm/PXD065579_KGG/.

### Phase 3 — Astral TMT model (hcd_astral_tryp_tmt), user-added 2026-06-23
TMT-on-Astral is a distinct regime (reporter ions, TMT/TMTpro fixed on K+Nterm, shifted b/y); bundled TMT models are LOW-RES CID only -> no Astral match. Ties to the historical TMT weak-spot. Slug hcd_astral_tryp_tmt, keyed (beam-CID, Astral, trypsin, TMT). Need a TRAIN + a HELD-OUT VALIDATION Astral-TMT-DDA pair (independent; entrapment-FDP). Candidates from base discovery (flagged as TMT, excluded from label-free base — wanted HERE): PXD058918, PXD060332, PXD062520, PXD063977, PXD055796. Discovery agent aef3296b scoping: confirm Astral+DDA (not DIA, PRIDE DDA facet unreliable), TMT plex, MS2-reporter vs MS3 (prefer MS2), has-IDs, recommend train+holdout. NOTE the low-res TMT benchmark PXD007683 is a DIFFERENT regime (not a validation set for Astral-TMT).
CAMPAIGN = Astral model family: hcd_astral_tryp (LFQ base) + PTM-aware (KGG/phospho diagnostic ions) + hcd_astral_tryp_tmt (TMT).

### Phase 3 scoped (aef3296b) — Astral-TMT viable + MS2-reporter (andes-friendly)
Clean Astral-TMT-DDA exists; Astral has NO SPS-MS3 -> reporters are in the high-res MS2 (the regime andes models). All TMTpro (16/18), LysC+Tryp, TMTpro on K+Nterm + Cam-C static, Ox-M variable.
- TRAIN -> PXD058160 (Tian Zhang, mouse aging, TMTpro-18, ~74 raws, per-raw mzIdentML labels READY, single-instrument, independent). Human alt: PXD063977 (Paulo, TMTpro-18, 52 labeled raws).
- VALIDATE (held-out, independent lab) -> PXD055796 (Emory ALS, human, FragPipe+Percolator results ready). Caveat TMTpro-16 vs train-18: only the reporter cluster (126-135) differs; b/y identical (+304.2071), and the scorer ignores the reporter region -> harmless. Strict-plex alt: filtered-Astral subset of PXD062520 (human TMTpro-18).
- Multi-instrument sets (PXD062520/PXD060332) need Astral-only filtering. PXD058918 mixes DIA+TMT (use TMT raws). Slug hcd_astral_tryp_tmt.
★ HARVEST DEFERRED (~100-150GB) until Phase-1 base PROVES analyzer-matched model beats bundled on benchmark @ honest FDP. Discipline: prove the thesis on the cheap LFQ base before the big TMT/PTM downloads. (KGG 38GB already staged = cheap.)

### Phase-1 A/B VERDICT (2026-06-23) — own-trained hcd_astral_tryp vs bundled
On the ProteoBench LFQ_Astral_DDA_15min_50ng benchmark (ASTRAL_entrapment.fasta, Percolator @ q<=0.01), identical binary/settings:
- bundled hcd_qexactive_tryp (MS-GF+-derived): 37,101 real @ 0.91% FDP, 16:41
- NEW hcd_astral_tryp (FULLY own-trained, 0 seed tables, 218,738 PSMs): 36,685 real @ 0.91% FDP, 17:33
NULL: own-trained is -416 PSMs (-1.1%), effective tie/slight loss vs the bundled MS-GF+-derived model. Thesis (analyzer-matched Astral model beats bundled) NOT supported on N=1. CAVEAT N=1 (one mzML, could flip). The bundled high-res model already serves Astral well.
★ REFRAME (user): the release-relevant gate is NOT vs the bundled MS-GF+-derived model but vs the FIELD (Comet/MSFragger/MS-GF+). If 36,685 beats them -> the own-trained model is FIELD-BEATING + independence-clean = strong release story regardless of the 1.1% vs bundled. Multi-engine benchmark on the SAME file LAUNCHED (ac259ac3). Prior Astral benchmarks (andes +~30% over MSFragger) were on a DIFFERENT file -> measure on this one.
