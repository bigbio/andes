# SP-B — glyco-specific rank model (design + staged plan, 2026-07-04)

## Why (diagnosis, from the collapse A/B)
Generation finds ~80% of truth backbones; only ~19% rank #1. The top-1-per-scan
collapse (required for honest TDC FDR) picks the winner by `rank_score` — andes'
`RankScorer` b/y score. That scorer is trained on regular high-res DDA, where
peptide backbone b/y are strong. In glyco spectra the precursor energy goes into
glycan cleavage, so backbone b/y are SUPPRESSED — a regime mismatch. The
collapse therefore mis-ranks on noisy b/y and discards correct backbones BEFORE
Percolator. A naive fix (make the core-Y ladder the primary collapse key) was
A/B-refuted (260→197). The real fix: a `rank_score` that models glyco b/y
statistics → better collapse → more correct backbones survive → more IDs.

## Architecture (reverse-engineered)
```
glyco collapse ── rank_score ◄── score_psm(ScoredSpectrum, peptide, RankScorer)
RankScorer ◄── RankScorer::new(&Param)
Param      ◄── model-train::estimate (per-fragment rank-LLR frequency vectors)
estimate   ◄── CountStats ◄── StatsAccumulator.accumulate(spectrum, LabeledMatch)
LabeledMatch { spectrum_index, peptide, charge }   ← THE input we must supply
store      ── Parquet model store; `andes --model-id <id>` loads a specific model
```
`andes train-from-search` already wires estimate→store, but labels via a STANDARD
search (`bootstrap_labels`) — not glyco. The single new piece is a path that
produces `Vec<LabeledMatch>` from GLYCO backbone PSMs and runs the same
accumulate→estimate→store. Everything downstream (RankScorer, --model-id) reuses.

## Staged plan
- **Stage 0 — feasibility gate (leaky, fast falsification).** Label from andes'
  OWN glyco top-1 @train-fdr PSMs on Fc3_r1 (self-training) OR the on-VM
  the reference engine Fc3_r1 IDs (glyco_cmp, 8210 PSMs). Train a glyco rank model,
  re-search Fc3_r1 with it (`--model-id`), measure backbone-correct vs the 101
  baseline. Train==test run ⇒ OPTIMISTIC, not a shippable number — but if a
  glyco-trained rank model can't beat the DDA model even in-domain, the
  direction is dead and we skip the expensive honest data. If it wins, proceed.
- **Stage 1 — honest cross-dataset.** Train on PXD005411 (a glyco search engine2 mouse-brain,
  SAME stepped-HCD regime, 17,188 public PSMs; FTP in 40-data/collection/
  README). Test on PXD025455 Fc3_r1. Clean (independent dataset/organism, same
  fragmentation). This is the shippable number.
- **Stage 2 — integrate.** Default the glyco driver to the glyco model (selection
  key / bundled store entry); keep the DDA model for non-glyco.

## Guardrails
- FDR stays Percolator-only on the 1-PSM/scan PIN (unchanged). SP-B only changes
  `rank_score`, i.e. WHICH backbone the collapse keeps — not the FDR.
- One variable per A/B; label the leaky vs honest numbers explicitly.
- TDD the new labels→accumulate path; reuse estimate/store as-is.

## Success gate
Stage 1 backbone-correct > 101 (and total @1% FDR ≥ 260) on Fc3_r1.

## STAGE 0 RESULT (2026-07-04) — FAILED; direction refuted
Trained a glyco rank model on the on-VM the reference engine Fc3_r1 backbone IDs
(8136/8136 labels parsed, 14 partitions), re-searched Fc3_r1 with `--model`:

| model | @1% FDR | truth scans | backbone-correct |
|---|---|---|---|
| DDA (baseline) | 260 | 237 | 101 |
| glyco-trained (leaky) | 239 | 218 | 94 |

**Worse — even leaky (train==test run).** So retraining the b/y RANK model does
NOT fix the ranking bottleneck. The regime-mismatch hypothesis is refuted: the
issue isn't a miscalibrated DDA model, it's that glyco b/y ions are
fundamentally SPARSE — training on them yields NOISIER per-rank statistics, not
better ones. Gate failed ⇒ stage 1 (PXD005411 download+train) is NOT worth it.

**Combined with the collapse A/B**, two single-signal collapse fixes are now
both dead: Y-ladder-primary (260→197) and glyco b/y rank model (260→239). The
"which backbone" signal is split across b/y + Y-ladder + oxonium; no single one
ranks. Remaining principled option: a LEARNED MULTI-FEATURE collapse score over
ALL glyco features (the 2-pass Percolator re-collapse — use pass-1 Percolator's
weights to re-pick the collapse winner from the top-K accepted candidates). The
`--labels` training path + collapse_cmp refactor are kept (reusable); the SP-B
rank-model arc is closed.

## CORRECTION (2026-07-04) — the stage-0 "refutation" was CONFOUNDED; WITHDRAWN
User pushback ("this has been mainly errors") + a Codex adversarial review +
self-audit found the stage-0 A/B changed TWO variables: the glyco model had
retrained frequencies AND an OWN-DERIVED geometry (train-from-search default),
which COLLAPSED to 14 partitions on the sparse glyco-only corpus vs the seed's
148. Re-ran with `ANDES_SEED_GEOMETRY=1` (frequencies-only):

| model (leaky, the reference engine Fc3_r1 labels) | @1% FDR | backbone-correct | partitions |
|---|---|---|---|
| baseline DDA | 260 | 101 | bundled |
| stage-0 glyco (own geometry) | 239 | 94 | 14 |
| glyco (seed geometry) | 255 | 93 | 148 |

Geometry fixed the COUNT (239→255 ≈ baseline) but backbone-correct stayed ~93
(< 101). So glyco frequencies are NOT dramatically worse (that was a geometry
artifact) but do NOT beat baseline on ranking either — inconclusive-leaning-
negative, NOT a clean refutation. **The stage-0 conclusion is WITHDRAWN.**

Still-uncontrolled confounds (Codex): single-engine (the reference engine) labels,
post-Percolator @1% metric (not decoy-separated ranking separation), same-run
label/eval overlap. A trustworthy call needs multi-dataset eval + multi-tool
consensus truth + a ranking-separation metric. ⇒ Do the dataset harvest
(PXD005411 a glyco search engine2, PXD016175 a glyco search engine2, PXD030670 a commercial glyco engine) FIRST, then re-test.
Do NOT close the direction on current evidence.

## HONEST RE-TEST (2026-07-04) — SP-B rank model refuted, controls CLEAN
Redid the test with every stage-0 confound controlled: seed geometry
(`ANDES_SEED_GEOMETRY=1` → 148 partitions, no degeneration), held-out labels
(the reference engine 8011 labels EXCLUDING the 196 consensus test scans → no test leakage),
measured vs BOTH truths.

| model | @1% FDR | vs 523 correct | vs 196 correct |
|---|---|---|---|
| baseline DDA | 261 | 101 | 51 |
| honest glyco rank model | 235 | 83 | 42 |

Worse on EVERY metric, cleanly. **The SP-B rank-model direction is now honestly
refuted** (the earlier confounded 239 is superseded; this is the trustworthy
verdict). Reason: the DDA rank model carries RICH b/y intensity-rank statistics
from non-glyco spectra; glyco b/y are SPARSE, so training the rank model only on
glyco data yields NOISIER per-rank vectors, not better ones. You cannot learn
better b/y ranking from a signal that is barely present. The discriminative
"which backbone" information is in Y-ladder + oxonium + multi-feature — which a
single b/y rank model structurally cannot capture. ⇒ Only remaining ranking
lever = the 2-pass Percolator re-collapse (learned MULTI-feature collapse score).
The `--labels` path + an open-source glyco engine oracle stay as reusable infra.

## 2-PASS PERCOLATOR RE-COLLAPSE (2026-07-04) — worse; re-ranking avenue CLOSED
The "last ranking lever": emit the full multi-row PIN (ANDES_GLYCO_ALL_HITS=1,
~333k candidate rows / 3108 scans, ~107 cand/scan), pass-A Percolator discriminant
over all candidates, re-collapse per scan by that discriminant, pass-B honest FDR.

| collapse selector | @1% FDR | vs 523 correct | vs 196 correct |
|---|---|---|---|
| rank_score (baseline) | 261 | 101 | 51 |
| Percolator multi-feature discriminant | 259 | 85 | 43 |

Worse. And the reason is FUNDAMENTAL, not an implementation detail: within a scan,
ALL backbone candidates for a target peptide are labeled TARGET; the decoys are
reversed PEPTIDES, not wrong BACKBONES. So the target/decoy discriminant learns
"real peptide match vs reversed decoy" — it has NO signal for "which backbone is
correct." Percolator (a target/decoy tool) structurally cannot solve within-scan
backbone selection. Empirically plain rank_score (direct b/y match quality) beats
the multi-feature discriminant as a per-scan selector.

### Ranking-by-re-selection is now EXHAUSTED (all refuted, clean):
- Y-ladder-primary collapse → 260→197
- glyco b/y rank model (honest, controlled) → 51→42 (consensus)
- 2-pass Percolator multi-feature re-collapse → 51→43 (consensus)

rank_score is the best available per-scan selector; you cannot out-rank it by
re-selecting from the SAME candidate set, because the correct-vs-wrong-backbone
signal is neither in the sparse b/y NOR in the peptide target/decoy labels. The
only remaining levers change the INPUTS, not the re-ranking:
1. GENERATION — fewer wrong backbones competing per scan (tighter candidate gen),
   so rank_score has fewer ways to be wrong.
2. GLYCAN-AXIS decoys (GI-2 part 2) — an ISOBARIC-composition glycan decoy makes
   the composition features discriminate glycan-correctness (a DIFFERENT
   target/decoy axis than peptide reversal). The one untested lever.

## GI-2 PART 2 (2026-07-04) — separate glycan axis REFUTED (2 decoy versions)
Added the negated-SialicConsistency glycan decoy (was YLadder-only) so the
composition features discriminate. 2D-FDR still 0: all 261 peptide-axis passers
fail the glycan axis, vs BOTH 523 and 196 truths.

FUNDAMENTAL reason: a glycan-decoy row shares its PEPTIDE BACKBONE with the
target, so all ~40 peptide features are IDENTICAL; only glycan-specific features
can differ, and of those GlycanMass is isobaric (same), core-oxonium is
composition-independent (same) — leaving just YLadderScore + SialicConsistency.
Two (individually weak, sialic-only-for-sialylated) features cannot separate
1127 target/decoy pairs at 1% FDR. A separate glycan axis is structurally
underpowered for andes's feature set.

AND it is redundant: the glycan features (YLadder, sialic, oxonium, GlycanMass)
are already columns in the PEPTIDE-axis PIN, and on that axis a target vs its
REVERSED-PEPTIDE decoy have DIFFERENT backbones → different glycan-by-subtraction
→ those features DO discriminate. So glycan correctness is already partly
controlled by the single peptide-axis FDR (the 261/101). Making the glycan axis
work would need MANY more composition-specific features (per-composition oxonium,
etc.) each with a decoy — a large feature-engineering effort, not a quick lever.

## MORE-IDs INVESTIGATION — every identified lever now tested (summary)
| lever | result |
|---|---|
| Y-ladder-primary collapse | worse (260→197) |
| glyco b/y rank model (honest) | worse (51→42 consensus) |
| 2-pass Percolator re-collapse | worse (51→43); TD can't rank within-scan backbones |
| glycan-axis 2D-FDR (YLadder decoy) | 2D=0 |
| glycan-axis 2D-FDR (+sialic decoy) | 2D=0 (only 2 features differ) |
| generation | near-ceiling (~80% find-rate) |

andes is at an architectural ceiling for this data: 261 @1% FDR (MORE than
an open-source glyco engine's 222), truth validated by 2 engines, backbone-correct 101/523 &
51/196. The remaining paths are LARGE (composition-feature engineering, better
generation, learned models needing more labeled data), not quick levers.

Review-found code bugs (fix before trusting any --labels result):
- label ingestion has no duplicate-scan/charge/conflict guard (didn't bite here —
  the reference engine export was 1-per-scan — but poisons on rank-alternative exports)
- collapse_cmp is applied AFTER a rank-only truncation in glyco_search, so the
  true ladder winner can be truncated before the shared comparator sees it
- the fixed-mod decorator is residue-wide only (fine for Cam-C; wrong for TMT)

## REVIEW ROUND (2026-07-04) — CodeRabbit + Codex + local agent; fixes + robustness
CodeRabbit: 0 findings. Codex: needs-attention (5). Local agent: 1 medium + confirmed
non-bugs (sialic negation algebraically correct; collapse_cmp a proper total order).

FIXED (commit d890fe8f):
- [Codex+agent] `load_labels_from_tsv` marked a scan "seen" BEFORE parsing its
  peptide → an unparseable first row poisoned a later valid row for the same scan.
  Fixed (insert after parse). Did NOT affect SP-B: the the reference engine export was
  verified 1-row-per-scan (0 duplicates), so it never triggered.
- [Codex HIGH] Collapse winner was chosen from a rank-only TRUNCATED subset;
  with exhaustive mode OFF by default (effective_top_k=50) a high-ladder/rank-tie
  winner could be dropped before collapse_cmp. Fixed: select over the FULL accepted
  set. Rank-primary experiments (baseline / SP-B / 2-pass / 2D-FDR) UNAFFECTED (the
  rank-max always survived truncation); the **yladder-primary A/B warrants a re-run**
  (its ladder-max could have been truncated) — pending VM access.

DOCUMENTED (not fixed; assessed impact):
- [Codex HIGH] isotope-factoring top_k truncation at the widest precursor: real but
  the MEASURED real-data impact was ≤1 PSM (261→260 when factoring was introduced),
  and the whole pipeline is capped at effective_top_k=50 anyway, so no NEW loss.
- [Codex HIGH] `--labels` trained model writes SEED data_type (protocol/instrument
  only affect model-id + ledger). Affects model-store AUTO-SELECTION, but every
  experiment loaded via explicit `--model`, so it did NOT affect results. TODO for
  bundled-store integration: set trained_param.data_type + an NGlyco selection key.
- [Codex HIGH] "negated sialic is an artificial mirror decoy" — correct concern, but
  MOOT: the glycan axis yields 0 IDs, so no q-value is ever calibrated from it; the
  finding is that the axis is underpowered, which we already concluded.

ROBUSTNESS VERDICT: the negative conclusions hold for all rank-primary experiments
(baseline, SP-B rank model, 2-pass re-collapse, 2D-FDR). The ONE result to re-run
under the truncation fix was the yladder-primary A/B.

RE-RUN UNDER THE FIX (2026-07-04): confirmed. Rank-primary baseline reproduces
EXACTLY (261 @1% FDR / 101 vs 523 / 51 vs 196) — the fix does not regress the
default. Yladder-primary = 201 / 89 vs 523 / 46 vs 196 — still clearly worse than
baseline (46 < 51 consensus). The fix nudged it (197→201) but did not flip the
verdict. So EVERY negative conclusion now holds under the corrected collapse path;
the review round is closed.
