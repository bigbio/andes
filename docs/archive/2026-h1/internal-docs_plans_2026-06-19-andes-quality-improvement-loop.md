# andes Quality-Improvement Loop — living plan

**Started 2026-06-19** (self-paced /loop). Branch `feat/ptm-refinement-cascade` (HEAD a254be65). Goal: systematically hunt bugs/issues, benchmark across configs vs MSFragger (PSMs / complementary IDs / speed / memory), and maintain a prioritized quality plan. Motivation: external review (CodeRabbit) trivially found 7 real bugs our own reviews missed → assume more exist.

## How this loop runs
Each iteration: (a) advance a bug-hunt slice, (b) advance a benchmark slice, (c) fold results into §2/§3 below, (d) reprioritize. Background agents do the heavy work; this doc is the durable accumulator.

---

## 1. Bug / issue findings (status: OPEN / FIXED / WONTFIX)
Seeded from this session (all the CodeRabbit P1–P3 were verified TRUE and FIXED a254be65):
- [FIXED] P1a fixture drop; P1b fixed+variable terminal stacking; P1c TSV merged storage; P2a YAML override; P2b token validation; P3 mod-class location.
- [OPEN] TMT fixed-peptide-N-term mod not applied to protein-N-term peptides (`build_terminal_variants` ProtNTerm-xor-NTerm; pre-existing).
- [OPEN, minor] combined CAM+Acetyl-Cys refine PSM → `mod_class 99` (Δ=99≠42).
- [OPEN, design] RAM `strong_score` order-dependent (listwise features over capacity queue) → mmap not byte-identical.
- [OPEN, perf] heavy-multi-mod mmap >90 min (lazy over-expansion).
- [OPEN, quality] RAM nominal-bucket prefilter is LOSSY (~+1,997 matchable candidates dropped) = Option B, potentially more IDs.
### Iteration-1 bug-hunt (opus; 1 Critical, 3 Important, 9 Minor; ~15 cleared)
- [FIXED 2065b87d] **C1 (Critical)** `candidate_nominal_bounds` (match_engine.rs:256) widened the asymmetric precursor prefilter on the WRONG side (min by `right`, max by `left`) — verified vs `matches_precursor` (left=lighter/lower bound, right=heavier/upper). Inert for symmetric 20/20ppm (tests green) but LIVE for the calibrator's asymmetric tolerance → pruned real candidates → lost PSMs. Affects ram AND mmap (shared fn). Fixed: min↔left, max↔right.
- [FIXED 2065b87d] **I1 (Important)** mmap `index_cache_path` omitted the FIXED-mod set → CAM-C vs CAM-C+TMT collide on one cache file → stale mass-shifted index reused → silent ID loss. Fixed: hash fixed_mod_deltas (sorted) + a regression test.
- [FIXED 2065b87d] **I2 (Important)** same cache key omitted decoy strategy/seed → reseeded shuffle reuses stale decoys → corrupt FDR. Fixed: hash ALL proteins (incl decoy sequences), not just targets + prefix.
- [DEFERRED, downgraded] **I3** out-of-range precursor charge → all-zero PIN charge one-hot (pin.rs:433 iterates `charge_range`). Re-assessed: all-zero IS a valid "charge not in standard set" representation, not a correctness bug; a clamp would also conflate z=6/z=7. Real improvement = widen the PIN charge schema to the max OBSERVED charge (deliberate output-format change), low impact on QE/HCD data (mostly z2-3), matters more on Astral/timsTOF. Not a quick fix — deferred to an output-format pass.
- [OPEN, minor] M1–M9: table-model intensity fragment-charge-blind fallback; mod_site charge-1-only blind spot; unclamped GBDT NaN path; all-or-nothing calibration abort; DeltaRankScore rank-1 attribution under --refine; Peptide::from_str can't parse combined terminal mods (test-only); mmap N-term-Met multiplicity under-count; 0.5Da nominal-window narrowing; non-standard cleavage residue drop. (full detail: `.git/sdd/bughunt-iter1.md`)
- _(new findings appended by the bug-hunt agent each iteration)_

## 2. Benchmark matrix (PXD001468 b1931 + others as added; TRUE entrapment-FDP, best-per-scan)
Metrics per run: total PSMs@true-1%, modified subset, complementary IDs (overlap & unique vs MSFragger), wall-clock, peak RSS, FDP.
Configs to sweep: andes {default, default+refine(gate 0.01), ram vs mmap, ±chimeric, precursor-cal on/off}; MSFragger {closed, open}; (Sage/Comet as added).
### Iteration-1 (b1931, commit a254be65; Percolator Separate mix-max)
| Config | PSMs@true-1% | Modified | True FDP | Wall | Peak RSS |
|---|---|---|---|---|---|
| andes RAM default (CAM+OxM+Acetyl) | 27,291 | 8,836 | 0.50% | 1:41 | 6.83 GB |
| andes mmap default | 27,287 | 8,828 | 0.51% | **29:43** | 4.26 GB |
| andes RAM + --refine 0.01 | **28,926** | **11,411** | **0.29%** | 2:02 | 13.56 GB |
| MSFragger closed (prior) | 28,207 | 9,038 | 1.14% | — | — |
- **No regression** from the bug fixes (27,291 ≥ 26,674 baseline; +2.3% from Acetyl default).
- **--refine 0.01 is the strongest lever:** +5.9% PSMs, +29% modified, HALF the FDP, +21s — but **memory doubles to 13.56 GB** (cost to watch).
- **andes+refine LEADS MSFragger-closed** (+2.5% PSMs, 4× lower FDP, +26% modified) — CAVEAT: FDR methods differ (Percolator-Separate vs TDC) → iter-2 must re-percolate MSFragger + run pep_entrap_curve for an honest cross-engine FDP.
- Complementary (andes-RAM vs MF): both 25,841 · andes-only 1,450 · MF-only 2,366 (86% union overlap).
- ⚠️ **mmap 17.5× slower than RAM at the DEFAULT 2-mod config** (not just heavy-mod) → mmap-with-mods is impractical until the lazy over-expansion + per-candidate multiplicity re-probe (M4) are fixed. RAM default unaffected.

Known baselines (this session): andes default 26,674/8,348/0.49%; andes+Acetyl(C1) 27,625/9,479/0.43%; MSFragger-closed 28,082/8,954; Sage 27,513/3,779; MSFragger-open 37,060/5,397.

## 3. Prioritized quality-improvement plan
Merged with the **auto-backend / one-PIN / PTM-enrichment** roadmap (sibling plan `2026-06-19-andes-auto-backend-one-pin-ptm-enrichment.md`, Cursor multi-agent review — accepted; decisions in §4 below). Ranked:
1. **[DONE 2e9207df] G1** — `--refine-select-psm-fdr` default 0.10→0.01 + hidden (entrapment-leak fix).
2. **Refine memory (real lever).** [DONE aac1c037] Step 1: owned Pass-1 candidates (move, not `to_vec`) → eliminates the duplicated Pass-1 pool (~1.5 GB). [DONE 307a5921] Step 2: PRUNE the Pass-2 pool (up to ~10.8M) to winners-only before merge — compact pool + candidate_idx remap inside `run_refinement`, full pool dropped. Output unchanged (PSMs resolve to same peptide); 30/30 refinement + 140 search tests green. `[REFINE]` log now prints `pass2_candidates=N (pruned to K winners)`. **VERIFY on VM (post usage-limit): RSS drop + PIN byte-identity diff.** Remaining: unidentified-only peak retention (peaks still kept for all spectra under --refine).
3. **Task 4 — official ScanNr best-per-scan collapse.** Scan key already shared (refinement.rs:577); the gap is that `pep_entrap_curve.py` is row-level. Make ScanNr best-per-scan the documented downstream step + fix the bench harness to always collapse before the FDP curve. (No new engine plumbing needed.)
4. **Correctness sweep:** triage iter-2 scoring findings (N1/N2 chimeric feature parity = inherent-to-architecture, low priority; M3 NaN cosine = clamp), then iter-1 M1–M9. Keep hunting until a pass is empty.
5. **VM entrapment gate** (`phaseA_entrapment.sh`) — record Phase-A true-combined-FDP ≤1% on b1931 (the §I ship gate).
6. **`--refine auto` UX collapse** (hide the 5 refine flags; stays OPT-IN — no default-on flip per §4).
7. **Phase A.2** refine-all-scans (precursor-gated, not unidentified-only) — upgrades IDs.
8. **mmap over-expansion fix** (17.5× @ 2 mods; P3) — unblocks low-RAM story. **Option B** prefilter A/B (P5).
9. **Search-quality multipliers** (fragment recal, open-search) + **RT/intensity rescoring** — Phase B, only after the gate passes.

---

### Iteration-2b (b1931, commit 2e9207df) — ⚠️ METRIC CORRECTION
- **Memory:** default 6.74 GB, --refine 13.37 GB (+6.63 GB). **P1 (drop identified-spectra peaks) is DEAD** — saves only ~45 MB (cloned-peak subset is tiny). The driver is the **Pass-2 candidate pool** (10.85M candidates, processed as ONE chunk) + holding Pass-1 AND Pass-2 candidate Vecs simultaneously. Real memory fix = stream/drain Pass-2 candidates per-chunk + don't double-hold. Added `[REFINE]` logging (b4d774f9).
- **🚩 The "+107% / 58,174 PSMs" is a METRIC ARTIFACT, not a win.** 58,174 > 55,352 total MS2 scans → impossible under best-per-scan. The curve script `pep_entrap_curve.py` counts PSM **rows** (Pass-1 + Pass-2 both counted per scan = the double-assignment the doc warned about), NOT scans. The HONEST best-per-scan number is iter-1's `phaseA_collapse_fdp.py`: **refine 28,926 vs default 27,291 ≈ +6%** (plus a large MODIFIED re-explanation of already-identified scans). Default andes already matches MSFragger-closed (~28k).
- **Methodology fix (binding):** all multi-pass (refine/chimeric) counts MUST be best-per-scan (ScanNr-keyed) BEFORE the FDP curve. `pep_entrap_curve.py` is row-level → misleading for refine; always collapse first. The scan key is ALREADY shared (refinement.rs:577 sets global spectrum_idx) → Task-4's remaining work is making ScanNr best-per-scan the official downstream step, not new engine plumbing.
- 83.4% of scans (46,147/55,352) go to Pass-2 because the internal RAW-TDC anchor gate is ~3× pessimistic vs Percolator → wasteful re-search. Lever: gate on a better Pass-1 confidence proxy to shrink the Pass-2 set.

### Iteration-3 results (bf661881, cam_only Pass-1 + M,P,K refine tier; --checksum-verified build)
- **Honest best-per-scan @ ≤1% combined entrapment-FDP:** default 23,259 / 4,650 mod / 0.50% → **refine 25,704 / 8,309 mod / 0.26%** = **+10.5% PSMs, +78.7% modified, FDP IMPROVES** (0.50→0.26%). 2,680 per-scan double-counts removed by the collapse → refine candidates genuinely BEAT Pass-1 winners (not just additive). This is the cleanest honest refine win yet (cam_only isolates Pass-2's contribution; on the full Ox-M default the marginal gain is ~+6% since Pass-1 already finds Ox-M).
- **M,P,K cost (cam_only):** refine +0.57 GB (+13%), +9.8s (+12.5%) vs default — modest. (The iter-2b 13.37 GB was Ox-M-in-Pass-1, a different config — NOT the refine tier.)
- **Collagen test = NEGATIVE on b1931 (wrong sample).** M,P,K found 2,588 ox-P + 860 ox-K real PSMs (0.29% FDP) but ZERO map to the 47 collagen proteins; top sources HNRPU/PTBP1/TXLNA = chemical oxidation on non-ECM proteins. b1931 is HEK293 whole-cell lysate (low collagen). **The collagen protein-recovery hypothesis needs an ECM/skin/bone/tendon dataset to prove** — on common samples M,P,K adds real-but-chemical oxidation PSMs + candidate cost, no collagen payoff. (M,P,K decision stands per user; demonstrate value on an ECM set or reconsider dedicated tier.)
- ⚠️ **PROCESS GOTCHA:** `rsync -az` silently skipped 3 changed files → stale build; fixed with `--checksum`. ALL future VM syncs MUST use `rsync -az --checksum` + verify a known change is in the binary. (Recorded in experiment-hygiene memory.)

### Iteration-4 (in flight): collagen protein-recovery test (user-named datasets)
User pointed at **PXD001765** (mouse lung ECM/matrisome, Q Exactive HCD, trypsin, hydroxyproline var-mod — fits hcd_qexactive_tryp) and **PXD006579** (human/mouse ISDoT decellularized native ECM, LysC+trypsin). Running PXD001765's smallest ECM fraction (1.06 GB .raw, andes reads natively): andes CLOSED vs --refine → collagen proteins/peptides/PSMs recovered ONLY via hydroxy-P/K. This is the RIGHT sample to prove the protein-recovery axis (b1931 was HEK293, no collagen).
- **❌ First attempt INVALID — stale binary.** Agent silently used a June-18 `/tmp/phaseA-target` binary (PRE-dating M,P,K) because my dispatch omitted `--features thermo` (needed for .raw) so a fresh build couldn't open the file. It "concluded" the default tier lacks P/K — true of the OLD binary, false of HEAD. ALL its numbers (closed 6,540 / refine 17,637 / collagen flat / 15.8 GB) are discarded. ROOT CAUSE = same stale-binary class as iter-3 + missing thermo feature. **Lesson: .raw runs MUST build `-p andes --features thermo`; dispatches must HARD-FAIL not fall back to a pre-existing binary.**
- **✅ RE-RUN VALIDATED the insight (bench-collagen2, binary aac1c037 + thermo, `[REFINE] anchors=1672 unident=68975/71144 pass2_candidates=49455`):**
  - **Collagen hydroxy-P/K: CLOSED 0 → REFINE 109 PSMs** @ q≤0.01 (all Pass-2) on CO1A2/CO2A1/CO4A2. Textbook example `TGETGASGPP+15.995GFVGEK` (CO1A2, Gly-X-Y hydroxyprolyl site).
  - Collagen unique peptides **129→275 (+113%)**, 161 refine-only; collagen PSMs 253→433 (+71%). Overall best-per-scan **6,405→7,910 (+23.5%)**.
  - **Memory fix confirmed: 15.8 GB (stale pre-aac1c037) → 4.99 GB** (3.2×; small Pass-2 pool here so mostly the Pass-1-dup elimination).
  - **HONEST NUANCE:** collagen PROTEIN count FLAT (15→14, noise) — alpha-chains already detected closed via non-hydroxy tryptic peptides. The win on THIS sample is COVERAGE DEPTH (the otherwise-invisible hydroxyproline peptidome), NOT net-new proteins. Net-new-protein recovery needs a protein seen ONLY via hydroxy peptides → test on PXD006579 (purer ECM, LysC+trypsin) and/or under protein-level FDR.

### Refinement's PROTEIN-recovery axis (user insight, 2026-06-19)
PTM refinement's value is not only PSM/scan count — it recovers **proteins** a closed search misses entirely. Canonical case: **collagen / ECM**, dominated by hydroxyproline (P +15.995) and hydroxylysine (K +15.995); those peptides never match unmodified → the protein drops out. A scan that "re-explains" from a weak unmodified match to a correct hydroxylated one can also re-assign to a DIFFERENT protein. This rehabilitates the modest "+6% PSM" framing: the real win may be on the protein axis.
- **GAP 1 (eval):** benchmarks measure PSMs/modified only, NEVER protein-level complementary IDs. ADD a metric: proteins identified ONLY via modified/refined PSMs (not present in the closed-search protein set). Run on b1931 + ideally a collagen/ECM-rich set.
- **GAP 2 (tier content) — DECIDED + DONE (bf661881):** user chose to **extend Oxidation residues to M,P,K** in the DEFAULT tier (accepted tradeoffs: ~3-5× more oxidation candidates on every refine run; hydroxy-P/K grouped under the oxidation FDR class). Must MEASURE the Pass-2 memory/time delta + the protein-recovery gain on the next VM run; if the cost is too high on non-ECM samples, revisit (dedicated tier).
- **NEXT VM run (binding):** report (a) protein-complementary IDs (proteins only-via-modified), (b) Pass-2 candidate count + RSS + wall delta from the M→M,P,K change (the `[REFINE]` log now prints pass2_candidates), best-per-scan. Ideally add a collagen/ECM-rich dataset.

## 4. Decisions on the Cursor auto-backend review (2026-06-19)
Reviewed `2026-06-19-andes-auto-backend-one-pin-ptm-enrichment.md`. It is well-grounded (cites this loop's HEAD + numbers) and largely correct. Rulings:
- **§C/G1 — ACCEPT + DONE.** Default 0.01, flag hidden (not removed — kept as a hidden escape hatch rather than a breaking removal). Evidence is decisive.
- **§B/G2 Task 4 — REVISED per user (2026-06-19): NOT in-engine one-winner.** Keep multiple candidates per scan (Pass-1 + Pass-2) as COMPETITORS and let Percolator's rescored score decide the winner — best-per-scan happens DOWNSTREAM (post-Percolator), not pre-Percolator in-engine. Andes stays a pure PIN emitter ("push more than one candidate for percolator to decide"). Task 4 work = (a) ensure Pass-1/Pass-2 candidates for a scan share the scan key so Percolator competes them, (b) make post-Percolator best-per-scan the OFFICIAL documented pipeline (not just `phaseA_collapse_fdp.py`). Verify current PIN SpecId/ScanNr grouping for Pass-1 vs Pass-2.
- **§E P1 (peak drop) — ACCEPT, sequence next.** The refine 2× RSS is real (measured 13.6 GB).
- **§D `--refine auto` collapse — ACCEPT, defer** until after P1+Task4.
- **§F items — ACCEPT into backlog** (BASEPEP accession masking blocks grouped FDR; TMT-skip WARN; mod_class 99; TMT protein-N-term bug).
- **§J1 default-on — REVISED per user (2026-06-19): `--refine` STAYS OPT-IN indefinitely.** No default-on flip. Predictable default behavior; PTM discovery is explicit.
- **Loop focus (user 2026-06-19): EXECUTE THE ROADMAP** (G1✓→P1→Task4→entrapment-gate as code, benchmark per step). Bug-hunting drops to background.
- **§J2 (accept ~14% fewer modified on stress for 15× better FDP) — YES** (honesty > inflated counts; on b1931 production 0.01 gave *more* modified anyway).
- **§J3 (in-engine FDR vs mokapot) — stay downstream for MVP.**
- **Iter-2 scoring bug-hunt (separate):** N1/N2 (chimeric secondary feature parity) judged inherent-to-architecture, low priority; M3 (NaN cosine) → clamp; N3 verify intent. Detail in `.git/sdd/bughunt-iter2.md`.

## Iteration log
- **VM validation (2026-06-19, HEAD a1a4b24b, b1931 default mods):** ✅ refine RSS **7.38 GB** (was 13–15 GB; §I 10 GB gate PASS), `[REFINE] pass2_candidates=393,323 pruned to 131,343`. No regression: refine 29,902 PSMs (+8.3%) / 12,647 mod (+41.8%), no unmodified-anchor bleed. **Entrapment gate PASS** (default 0.43%, refine 0.29%). Safety batch (panic/DoS/cache-validation) confirmed in binary. PXD006579 skipped (VM /tmp disk). The session's memory + safety + FDR-honesty work is validated.

- **Iter 1 (2026-06-19):** bug-hunt → **fixed C1 + I1/I2 (commit 2065b87d, search 140/140)**. Benchmark: andes RAM default 27,291 (no regression), andes+refine 28,926/11,411mod/0.29% (appears to lead MSFragger-closed but FDP cuts unmatched), mmap 17.5× slower at 2 mods. I3 downgraded (output-format, deferred).
- **Iter 2 (2026-06-19):** scoring bug-hunt → core clean; only narrow chimeric/strong feature-parity items (N1/N2 architecture-inherent, M3 NaN clamp). FDR-parity bench died on a session limit. Reviewed Cursor auto-backend doc → decisions §4; **shipped G1** (2e9207df, refine gate 0.10→0.01).
- **Iter 2b (2026-06-19):** VM memory profile + bench re-run. **P1 peak-drop killed** (45 MB; driver = Pass-2 candidate pool) → added `[REFINE]` logging (b4d774f9). **🚩 caught the +107% as a row-counting artifact** (pep_entrap_curve.py counts rows not scans; 58,174 > 55,352 scans) — honest refine gain is ~+6% best-per-scan. Corrected §2 + methodology.
- **Iter 3 (2026-06-19):** user insight → refine has a PROTEIN-recovery axis (collagen hydroxy-P/K). Shipped **Oxidation→M,P,K** (bf661881) + **owned-Pass-1-candidates memory fix** (aac1c037, no `to_vec` dup). Dispatched VM run measuring M,P,K cost (pass2_candidates/RSS) + protein-complementary IDs + honest best-per-scan. NEXT: confirm aac1c037 RSS drop; Pass-2 pool pruning; official ScanNr collapse.
