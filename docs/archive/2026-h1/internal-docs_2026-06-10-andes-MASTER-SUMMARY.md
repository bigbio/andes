# Andes R&D Master Summary — 2026-06-10

> Synthesized from 47 internal research documents (2026-05-28 → 2026-06-10).
> Newer docs supersede older where conflicts exist. All PSM counts at 1% FDR (entrapment-controlled
> unless noted). Numbers are concrete; "refuted" means built + A/B-tested and found harmful or flat.

---

## 1. OVERVIEW

**Andes** is a clean-room Rust port of MS-GF+ (formerly SIMAS, then CIMAS), targeting full
Apache-2.0 relicense. The core scoring algorithm is the **Frank 2005 intensity-rank LLR**:
`RawScore = Σ_i log(P_sig(rank_i) / P_noise(rank_i))` — a difference of two per-ion KL terms.
This is NOT XCorr/Hyperscore/MVH; replacing it without proof is explicitly forbidden.

**Three phases of independence work:**
- Phase 1 (DONE): Remove the MS-GF+ generating function (GF) from IP-critical code paths.
- Phase 2 (DONE): Clean-room scoring pipeline (Rust, additive PIN features, no Java source).
- Phase 3 (~1/39 DONE): Retrain all 39 model slugs from open data. Cannot relicense until complete
  (models.parquet is still MS-GF+-derived for 38/39 — claiming independence would be a false IP claim).

**Canonical benchmarks (2026-06-09, PR #51, all at 1% entrapment-FDR):**

| Dataset | Andes | Next best | Notes |
|---|---|---|---|
| Astral (chimeric) | **79,652** PSMs | — (LEADS all 3) | +31% vs May baseline |
| UPS1 CID LFQ | **17,532** PSMs, **986** proteins | — (LEADS) | |
| TMT PXD016999 | **12,000** PSMs | MSFragger | 2nd overall |
| TMT PXD007683 (a05058 Lumos) | 46,475 PSMs | MSFragger **52,528** | −11.5% gap |

---

## 2. SHIPPED / WORKS

### 2a. Core engine (all merged to dev/main or PR #51)
- **Native Thermo .raw** (PR #44): in-process .NET 8 `thermorawfilereader`; auto-dispatch by extension.
- **timsTOF .d**: also shipped.
- **Chimeric cascade** `--chimeric` (PR #42): two-pass search; Astral +101% PSMs vs single-pass.
  Entrapment-validated (FDP 1.10% → 1.12% — flat).
- **Parquet model store + `msgf-rust train`** (PR #47): 39 `.param` files consolidated into
  `models.parquet` (byte-identical round-trip); `--update --add/--remove-source/--reweight/--decay
  --validate`. Per-source sufficient stats for exact incremental updates.
- **Precursor calibration** `--precursor-cal {auto,on,off}` (PR #33): MassCalibrator.
- **Dense noise model + windowed peak filter**: shipped in production models; key for Astral quality.
- **`keep seed mme`**: critical rule — model fragment tolerance is per-model (e.g.,
  `hcd_qexactive_tryp = 0.5 Da`, NOT ppm). Forcing `--fragment-tol-ppm` on own model collapsed
  Astral from 28,773 → 533 PSMs. Always keep the seed mme when retraining.

### 2b. Additive PIN features (shipped, Percolator reweights these)
- **DeltaRawScore**: +129 PXD001819 / +12 TMT / +104 Astral PSMs. Zero wall cost.
- **PpmGaussianScore** (+0.54 Percolator weight on Astral, top fragment feature).
- **ChanceMatchSurprise** (listwise null, +0.26 weight).
- **LongestComplementaryLadder** (+0.05 weight).
- **Precursor isotope envelope**: chimeric disambiguation.
- **TailorScore**, **StrongScoreCal**: per-spectrum calibration columns.
- **Sanov significance** (PR #55, `feat/sanov-significance`, SHIPPED): two additive features
  — `SanovLogEvalue` (= `m·I(t) − ln N_cand`, Cramér rate × ion count) + `TypeAtypicality`
  (= `D(P̂‖P_sig)`). **+467 PSMs (+0.70%) on Astral at flat FDP.** Non-redundant with RawScore
  (corr = −0.14/−0.35). Neutral on low-res TMT (+15, +0.13%) — TMT is evidence-bound, not
  significance-bound. Mechanism confirmed: recovered PSMs are mean length 9.7, median 8.

### 2c. Speed (already 25-40× faster than Java on Astral)
- No significant new speed work shipped in this period; the roadmap exists (see §7).

### 2d. Training pipeline (operational)
- **MSNet catalog**: 112 datasets, 321.9M PSMs. HCD dominates (86 datasets, 247M PSMs).
- **15-file bootstrap** for `hcd_qexactive_tryp` → 28,773 PSMs, gap −16% vs curated 34,176
  (pipeline proven end-to-end; quality gap = corpus quality, not algorithm).
- **ion_existence bug fixed** (local, uncommitted as of 2026-06-05): `accumulate()` never called
  `bump_existence` → ion_existence stuck at uniform 0.25. Fixed: iter4b 23,949 → 28,859 (+20.5%).

---

## 3. REFUTED / DEAD ENDS

All items below were built and A/B-tested (entrapment FDP controlled). "Refuted" = PSMs flat or
negative, or FDP inflated.

| Idea | Result | Root cause |
|---|---|---|
| **H2 peak-rank modifying** (TMT gap) | REFUTED | Modifying existing distributions regresses (n=12 parity iters) |
| **TMT fragmentation overlay** (tag-loss, low-mass cutoff, existence reweight) | REFUTED (all 3) | 2a tag-loss −8.6%, 2b −1%, 2c ±0 on a05058 |
| **EM/iterative TMT labels** | DEAD END | Training loop; no PSM gain |
| **More QExactive data** (diversity lever) | REFUTED | Data was not the lever |
| **MSnet for low-res CID diversity** | REFUTED | Domain lock-in; no Lumos+CID in MSNet |
| **Fragment-ion inverted index** (top-K prefilter) | REFUTED for chimeric | 33% faster, 46% recall — drops chimeric secondaries |
| **Sage-style top-K prefilter** | REFUTED | Same recall/speed tension |
| **ImmutableCollections** | REFUTED | 2.2× regression |
| **Fragment-vote-all-touched index** | REFUTED | PXD 28min+ (too slow) |
| **Noise sharpening / noise_pseudo knob** | REFUTED | No PSM gain |
| **Isotonic regression calibration** | REFUTED | No PSM gain |
| **Temperature/score sharpening** | REFUTED | No PSM gain |
| **Global peak cap** | REFUTED | No PSM gain |
| **Piecewise scoring fixes** (Java-mirror) | REFUTED | Regresses production; non-additive |
| **Strong-score gating** (H1) | REFUTED | +2.6% PSMs gained by keeping it → ablation overturned |
| **Fragment-ion isotope features** (M+1 presence + ratio) | REFUTED | AUC 0.660; −0 PSMs; redundant with PpmGaussianScore |
| **CalibratedScore** (length-null + per-spectrum LOO-z) | REFUTED | AUC 0.655 vs 0.682 RawScore; −22 PSMs; root cause: top-n 1 → LOO is no-op |
| **Residue-Lattice Rescue** | REFUTED on TMT | a05058: +2 PSMs (flat, 0.02%); density null absorbs all signal on low-res dense spectra |
| **Strong score on CID** | PARKED (much worse) | −23%/−41% on both TMT benchmarks |

**Key pattern (confirmed by 3 independent refutations 2026-06-10):** The per-spectrum additive
feature space is **saturated** on low-res TMT. Sanov, Residue-Lattice, and CalibratedScore all fail
to recover the 1,102 two-engine-agreed a05058 misses. The TMT bottleneck is NOT a per-spectrum
feature problem.

---

## 4. KEY LEARNINGS / RULES

### 4a. The Iron Rule (non-negotiable)
**Additive-only features.** Never modify the RawScore ranking. Percolator reweights additive PIN
columns safely; modifying existing score distributions regresses (empirically confirmed n=12+).

### 4b. Entrapment FDP is the only trusted ruler
`FDP = 2·N_ent/N_total` (paired `ENT_` entrapment DB). Reversed-decoy TDC is blind to coincidental
targets. Never gate on TDC alone. All validation must be entrapment-controlled A/B.

### 4c. Percolator mode caveat
Percolator 3.7.1 auto-detects Concatenated vs Separate from PIN row structure. Cross-mode counts
are not comparable — grep the mode before comparing PSM counts.

### 4d. Training rules
- Keep seed mme: fragment tolerance is per-model, not per-instrument globally.
- Dilution is real: Laplace pseudo + global backoff blends sparse partitions toward bulk datasets
  (MSNet HCD dominates). Phase-3 bootstrap model was −4.3% vs curated `cid_lowres` (15-file test).
- Domain lock-in: zero Lumos+CID, zero TMT+CID in MSNet → `cid_lowres_tryp` suffers.
- Frozen structure: StatsAccumulator must reuse production matchers verbatim.
- ion_existence bug: was unfixed for 38 models (uniform 0.25); fix before any retraining.

### 4e. Score identity
Andes `--score rank` IS the same algorithm as MS-GF+ rank scoring (Frank 2005 LLR). The Lumos
−11.5% gap is a **model-data mismatch**, not a wrong algorithm. Do NOT replace inner-loop scorer
with XCorr/Hyperscore/MVH without proof. Fix the training first.

### 4f. Chimeric cascade rules
- Pass-2: exclude primary peptide, fill secondary Tailor/delta, no spec.clone.
- Entrapment-validated (flat FDP).

### 4g. Additive feature pattern (4-touchpoint)
(1) New `crates/search/src/<feature>.rs` module; (2) call in `compute_psm_features`; (3) new
`PsmFeatures` fields (default 0.0); (4) new PIN columns in `crates/output/src/pin.rs` +
schema-parity test updated.

### 4h. Sanov identity
`RawScore = m·[D(P̂‖P_noise) − D(P̂‖P_sig)]` — point statistic. `SanovLogEvalue = m·I(t) − ln N_cand`
is the missing null tail (large-deviations exponent, NOT the GF; no DP, no graph; IP-distinct).

---

## 5. THE TMT GAP

### 5a. Measured state (2026-06-10)
- PXD016999 (Fusion TMT CID): Andes 20,857 vs MSFragger 20,760 (+97, Andes leads).
- PXD007683 a05058 (Lumos TMT CID): Andes 46,475 vs MSFragger 52,528 (−6,053, **−11.5%**).
- The gap is entirely on Lumos low-res CID, not on Fusion.

### 5b. Root cause (confirmed by diagnostic 2026-06-10)
3-way diagnostic on a05058 (Andes / Sage / MSFragger, same `TMT_entrapment.fasta`, 1% FDR):
- Andes 10,914 · Sage 10,746 · MSFragger 12,037 unique scans.
- Sage ∩ MSFragger agree on 9,897 scans; Andes misses **1,102 (11.1%)**.
- Bucket breakdown: **50% B (scored-but-not-accepted)**, 50% A (candidate-gen/deep-ranking).
- Missed peptides: mean length 10.9, median 9. **Short peptides, few fragment ions.**
- Bucket B: Andes scored the right peptide at RawScore mean 66.6 but ranked below 1%.

### 5c. Evidence ranking for the Lumos gap
1. **Wrong training domain** (no Lumos+CID in MSNet) — primary lever.
2. **Estimator dilution** (Laplace + global backoff; 15-file bootstrap −4.3% vs curated).
3. TMT distorts rank/noise distributions.
4. Flat rank distributions (curated `cid_lowres` vs own model: 0.5 Da now beats 0.05 Da after
   ion_existence fix; gap −27% → −18.7% with fix, residual = structural softness).
5. GF removal ~6% low-res contribution.
6. Algorithm/peak-pick parity — REFUTED.

### 5d. Per-spectrum features are exhausted for TMT (strong evidence, 3 independent failures)
- Sanov (significance recalibration): +0.13% a05058, +0.70% Astral.
- Residue-Lattice (interval-evidence): +0.02% a05058 (flat), density null absorbs all signal.
- CalibratedScore: −22 PSMs. AUC worse than RawScore.

### 5e. Remaining levers for TMT
1. **Phase-3 retraining** (domain-matched Lumos CID labels, anti-dilution estimator) — primary.
2. **Andes-Belief** (cross-spectrum belief propagation; A1 from novel-algorithms synthesis) — TMT's
   dense per-peptide corroboration makes the coupling graph rich; riskier (FDP inflation risk).

---

## 6. MODEL TRAINING & INDEPENDENCE (Phase 3)

### 6a. Status
- 1/39 slugs retrained from open data (`hcd_qexactive_tryp` — pipeline proven, quality gap remains).
- 38/39 models are MS-GF+-derived in `models.parquet`. Cannot relicense Apache-2.0 until all 39 done.
- Branch `chore/cleanup-msgf-heritage` (commits fd1b8086, 067082b1): LICENSE/NOTICE/README credit
  origin + code-independence + models-in-transition + intent to relicense.

### 6b. Replace bar
Own model ≥ curated OR within ≤3% AND entrapment-FDP ≤ curated.

### 6c. Three-tier data strategy
- **T1** (MSNet free sweep): 112 datasets, 321.9M PSMs; HCD dominates; no Lumos+CID, no TMT+CID.
- **T2** (PRIDE harvest via `build_gap_corpus`): disk-bounded, targeted domain-matched datasets.
- **T3** (generic pooled own-model): fallback when T1/T2 insufficient.

### 6d. Model architecture
- Parquet schema: `models` manifest + `tables` bulk + `sources` ledger + `stats` per-source counts.
- Slugs: standard, phospho, tmt, itraq, acetyl, ubiquitin, glyco, immuno.
- Selection backoff: exact → largest subset → labeling only → drop experiment class → instrument
  family → generic.
- Acceptance gate: entrapment-FDP ≤ curated (the provenance-assertion test is red until Plan-4 done).

### 6e. Anti-dilution (key)
- Set `MSGF_DENSE_NOISE=300` per-slug harness.
- Noise-sampling tweak: stop over-smoothing noise rank distributions (dilution is the noise rank
  model 2.5× flatter; signal b/y largely unchanged).
- `retrain_slug.sh` per-slug harness with seed mme preserved.

### 6f. CID/TMT skeleton (for Lumos gap)
- `cid_native_tryp` skeleton builder (Rust): generate model structure without Java seed.
- Structural skeleton = independence from MS-GF+ data; estimator + labels = yield.
- Target: train Lumos TMT CID model; success = Andes ≥ 49,900 PSMs on PXD007683.

---

## 7. OPEN IDEAS / BACKLOG

### 7a. HIGH PRIORITY (next bets, evidence-grounded)

**Andes-Belief (A1)** — cross-spectrum belief propagation.
- Factor graph over the run: PSM → peptide-presence → protein-presence + RT-coelution factor.
- 3 additive PIN features (safe). ~400-700 LOC post-scoring/pre-PIN.
- Kill-switch: RT-shuffle (if shuffled-RT still "helps", fitting noise → reject).
- Targets ALL THREE Andes weaknesses: TMT corroboration, low-abundance recovery, protein consistency.
- Risk: FDP inflation from cross-spectrum coupling — requires decoy-symmetry + RT-shuffle gate.

**Phase-3 retraining** — the primary TMT fix.
- Target: domain-matched Lumos TMT CID corpus via PRIDE harvest.
- PXD014502 (Lamond, CID.ITMS ion-trap-CID-MS2-TMT) identified as a candidate source.
- Discovery method: search PRIDE "CID ion trap TMT" → check description/files.

**Kalman drift filter (D1)** — addresses bucket-A candidate-gen misses.
- Track instrument state (mass bias, drift, log-noise) over RT.
- Static two-pass recal is the degenerate zero-drift case.
- Do in parallel with diagnostic (untested: are the absent peptides' masses inside Andes' window?).

### 7b. MEDIUM PRIORITY

**MinHash-LSH open search (E1)** — weighted-MinHash sketches via CWS (Pr[collision] = weighted
Jaccard). Banded probes → pre-thresholded shortlist. NOT the fragment index; escapes the refuted
top-K recall/speed tension. Cost is window-width-independent → makes open search feasible.

**Andes-Scout** — length-stratified secretary rescue of deep candidates. Partition candidates by
(charge, peptide_length), hash-defined scout sample defines length-local threshold. Only if shadow
diagnostic shows true short peptides have high lattice/feature Z but outside top-N.

**Fragment-Slab BVH** — exact bounding volumes over candidate prefix-mass trajectories. Admissible
upper bounds → only skip candidates that provably cannot beat current queue threshold. Best speed
idea for expanded searches (non-destructive to FDR/top-N).

**Auto-tolerance** — preflight calibration: metadata priors → sampled permissive search → confident
PSMs → mass-error models → validation → main search. Caution: do NOT change core node-scoring
tolerance unless model declares it.

**Multi-enzyme CleavagePolicy** — union cleavage positions; canonical policy key (Chymotrypsin+LysC,
TrypsinP); model store backoff.

**MS1 monoisotopic precursor correction** — documented +85% PSMs in literature; bucket-A lever.

### 7c. LOWER PRIORITY / DEFERRED

**Thermodynamic density-of-states (B2)** — more powerful than Sanov but heavier and GF-adjacent.
Escalate only if Sanov's asymptotic proves too loose.

**Non-negative group-sparse deconvolution (A2)** — convex joint solve (group-LASSO) for chimeric
spectra; generalizes greedy chimeric Pass-2 to a joint solve; peaks partitioned, not double-counted.

**Global submodular PSM↔spectrum assignment (A3)** — subsumes by A1 (Andes-Belief); higher risk.

**X!Tandem staged refinement** — Pass 1 (canonical) → Pass 2 (implicated proteins, budgeted PTMs).
9 lessons documented. Key caution: do not inherit silent invalid params, flat FDR pool for rare PTMs.

**Optimal transport (C1) / HKLS (C2) / TDA (C3)** — elegant mathematics but per-PSM additive
features in the saturated space. Revisit only if diagnostic shows a specific dimension current
features miss.

**Direct feature ablation** (cheaper than CalibratedScore): drop TailorScore / StrongScoreCal /
peplen one at a time and measure. Whichever drops with no loss is redundant → retire it.

**Strong score** `--score strong` — parked for CID. Stage-1 additive bolt-ons (+3.5% Astral via
strong-score track: 35,789 → 37,052 PSMs); ablation weights PpmGaussian +0.54, ChanceMatch +0.26,
Ladder +0.05. Top-K retention under strong score is a known limitation (rank-gated → only re-ranks
what rank_score kept).

### 7d. Backlog tags (from `2026-06-10-andes-scoring-ideas-backlog.md`)
- SHIPPED: dense noise model, windowed peak filter, diverse harvest, chimeric cascade, keep seed mme,
  strong-score features, precursor isotope, native .raw, Parquet model store, Sanov significance.
- REFUTED: noise sharpening, noise_pseudo knob, isotonic regression, temperature/sharpening, more
  QExactive data, MSnet for low-res CID diversity, strong-score gating, global peak cap, TMT
  fragmentation overlay, EM/iterative labels, fragment-index, piecewise scoring fixes,
  fragment-ion isotope features, CalibratedScore, Residue-Lattice.
- HIGH-PRIORITY OPEN: Andes-Belief, Phase-3 harvest, Kalman drift, MinHash-LSH open search.

---

## 8. ROADMAP / PRIORITIES

### Immediate (do first, highest leverage)
1. **Fix ion_existence bug** in production models (commit the local fix; required before any
   Phase-3 retraining or the baseline is wrong).
2. **Phase-3: Tier-1 pilot** — retrain `hcd_qexactive_tryp` with anti-dilution (`MSGF_DENSE_NOISE=300`)
   + ion_existence fix. Gate: own ≥ curated OR within 3% at FDP ≤ curated.
3. **Phase-3: Lumos TMT CID corpus** — identify/harvest PRIDE CID+TMT datasets (PXD014502 first);
   build `cid_lowres_tryp` native skeleton; retrain; target ≥ 49,900 PSMs on PXD007683.

### Near-term (parallel tracks)
4. **Andes-Belief (A1)** — brainstorm → build behind flag → entrapment+RT-shuffle gate.
5. **Bucket-A diagnostic** for a05058 — check if the 548 absent-peptide misses are precursor-window
   failures (→ Kalman/MS1-mono) or deep-ranking (→ Andes-Scout).
6. **Direct feature ablation** — drop TailorScore/StrongScoreCal/peplen one at a time on Astral +
   a05058; retire any that drop with no loss (reduce feature count without CalibratedScore overhead).

### Medium-term
7. **Phase-3: complete 39 slugs** → Apache-2.0 relicense milestone.
8. **Strong score top-K retention fix** — increase per-spectrum retention to top-K under `--score
   strong` (increase K from 1); then revisit the −325 Astral PSM loss.
9. **Auto-tolerance preflight** — precursor auto first, then fragment feature-only, then validation
   grid.
10. **MinHash-LSH open search** — separate capability track (no impact on closed-search quality).

### Architecture decisions (locked)
- Core scorer: Frank 2005 rank-LLR, untouchable without proof.
- Validation ruler: entrapment FDP, always.
- Feature addition: additive-only, never modify RawScore ranking.
- Apache relicense: gated on all 39 slugs passing provenance-assertion test.
- TMT: per-spectrum features exhausted; pivot to model + cross-spectrum (Andes-Belief).

---

## Appendix: Key numbers quick-reference

| Metric | Value | Source |
|---|---|---|
| Astral PSMs @1% (PR #51, chimeric) | 79,652 | 2026-06-09 benchmark |
| Astral PSMs @1% (strong-score track) | 37,052 | 2026-06-06 |
| UPS1 PSMs / proteins @1% | 17,532 / 986 | PR #51 |
| TMT PXD016999 Andes vs MSFragger | 20,857 vs 20,760 | PR #51 |
| TMT PXD007683 (Lumos) Andes vs MSFragger | 46,475 vs 52,528 (−11.5%) | 2026-06-07 |
| ion_existence fix: Astral iter4b | 23,949 → 28,859 (+20.5%) | 2026-06-05 |
| MSNet catalog | 112 datasets, 321.9M PSMs | 2026-06-04 |
| Phase-3 slug retraining | 1/39 complete | 2026-06-09 |
| Sanov Astral gain | +467 PSMs (+0.70%), FDP flat | 2026-06-10 |
| Sanov a05058 gain | +15 PSMs (+0.13%) | 2026-06-10 |
| Residue-Lattice a05058 gain | +2 PSMs (flat, +0.02%) | 2026-06-10 |
| CalibratedScore delta | −22 PSMs, AUC 0.655 vs 0.682 | 2026-06-10 |
| Fragment-isotope delta | −0 PSMs, redundant | 2026-06-10 |
| a05058 headroom (Sage∩MSFragger − Andes) | 1,102 misses (50% B, 50% A) | 2026-06-10 |
| Missed peptide median length | 9 (bucket B) | 2026-06-10 |
| Andes vs Java speed | 25-40× faster (Astral) | 2026-06-04 |
| DeltaRawScore gain (additive) | +129/+12/+104 PSMs (PXD/TMT/Astral) | 2026-05-28 |
