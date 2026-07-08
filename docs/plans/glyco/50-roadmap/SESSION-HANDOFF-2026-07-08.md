# Glyco Campaign — Session Handoff & Forward Plan (2026-07-08)

> **Purpose:** a clean resume-from-here document so the next session picks up without
> re-deriving the investigation. Test bed throughout: PXD025455 `HCC_pool_Late_Fc3_r1`
> (human serum, stepped-HCD, Q Exactive HF). 523 the reference engine-truth backbones.
> Authoritative measuring stick: `glyco_outrank_audit.py` on the VM (`/srv/data/msgf-bench/glyco_bench`).

---

## 1. Branch state — `glyco-phase1` (UNMERGED, safe to build on)

Tracked tree was deterministic (byte-identical proven) and twice code-reviewed
(Codex + CodeRabbit, all findings resolved). Current worktree also has untracked
planning/audit files from this campaign:
`docs/plans/glyco/50-roadmap/BEAT-the reference engine-PLAN-2026-07-08.md` and
`docs/plans/glyco/scripts/glyco_outrank_audit.py`. Still **beats an open-source glyco engine:
253 @1% / 97 backbone-correct / 1 decoy vs MM ~222**.

Net code delivered this session (on top of the prior 253/97 baseline):
- **RT foundation** — engine-wide + glyco retention-time features
  (`DeltaRT`/`AbsDeltaRT`/`DeltaRTNorm` + glyco `DeltaRTRank`), per-run
  self-calibration. Commits `db85f016` (Commit-1), `08577538` (Commit-2).
  Core in `crates/scoring/src/rt_model.rs` + `crates/andes-glyco/src/rt_offset.rs`.
- **RT review fixes (2 rounds)** — `1fc9c2a1` (determinism: glyco hits total-order
  sort; NaN guards; anchor hygiene) and `db10abd1` (NaN guards completed at the
  caller boundary; glyco anchor decoy-skip parity). Determinism VM-proven.
- **Reverted experiments (net zero, do not resurrect):** `MIN_BY` knob
  (`159a9941`→`4a31ff6d`), isotope-voting fix (`7b765e95`→`357e8532`).

**RT is BUILT, SAFE, but PARKED** — measured a weak discriminator with the untrained
seed model (@1% is byte-identical with RT populated). It does not move the number yet.

**Do NOT open a PR** until we beat the reference engine or have a clean strategy (user directive).

---

## 2. The governing truth after the the reference engine-gap audit

The previous "generation vs scoring" story was too coarse. The authoritative
measurement is still `glyco_outrank_audit.py` on `andes_cap1024_truth_allhits.glyco.pin`,
but interpret its buckets carefully: `truth_absent` means the correct backbone was
not present in the emitted/all-hit candidate surface under that run, NOT necessarily
that the database could never enumerate it.

Decomposition of the 523 truth backbones under the all-hit audit:

> **⚠️ CORRECTED 2026-07-08 (numeric mass matching).** The original counts below came
> from a `glyco_outrank_audit.py` build whose `residue_mass()` STRING parser mis-parses
> **45.9%** of candidate rows (measured over 89,412 rows; e.g. a 3947 Da backbone read
> as 146 Da), inflating `truth_absent` with false absents. Both audit tools are now
> fixed to match numerically via **CalcMass − GlycanMass − H₂O** (repo + VM synced).
> The corrected decomposition changes the whole diagnosis: **generation is ~half the
> believed size; the dominant problem is SCORING.**

| bucket | buggy (string) | **CORRECTED (numeric)** | meaning |
|---|---|---|---|
| truth_absent | 301 (58%) | **148 (28%)** | correct backbone not on the scoreable surface |
| truth_outranked | 106 (20%) | **95 (18%)** | in the pool, but a *wrong* backbone wins top-1 |
| top1_correct (by mass) | 116 (22%) | **280 (54%)** | winner has the correct backbone mass; only **97** survive @1% FDR |

**Headline:** andes already generates AND top-1-selects the correct backbone mass for
**280/523**, but only **97** survive @1%. The **280→97** gap (scoring + FDR-separation,
plus some mass-correct/peptide-wrong) is the real battleground — NOT generation. The
"generation-bound" framing below was largely a parser artifact.
(`top1_correct` is backbone-MASS-only — the truth file has no sequence — so 280 is an
upper bound; the honest peptide floor stays 97 @1%.)

Established beyond doubt this session:
- **It is NOT an FDR threshold problem.** Loosening 1% -> 40% moved
  backbone-correct only 97 -> 103 while decoys exploded. Percolator cannot rescue
  a correct backbone that was not emitted as the one PSM for the scan.
- **It is NOT the mods-file problem.** Source audit shows the no-`--mods` default
  already includes protein N-term acetylation in `default_aa_set_with_tag()`.
  The nearby log string that says "Cam-C fixed, Ox-M variable" is stale/incomplete.
- **It is NOT solved by the current combined selector.** Diagnostic top-1 improved
  modestly, but the honest Percolator run on the existing `honest_comb` PIN gave
  **236 @1% / 86 backbone-correct / 1 decoy**, worse than 253/97/1. Also,
  `ANDES_GLYCO_SELECTOR=combined` implicitly turns `YINDEX` on unless overridden,
  so prior combined-selector A/Bs changed candidate retention and scoring together.
- **The miss pattern is large/high-charge backbones.** Correct IDs skew smaller
  (median backbone mass ~1350 Da); absent/outranked are ~1750-1860 Da. z=2 is
  mostly recoverable, z=5+ is essentially not recovered.
- **Generation and scoring are coupled.** Adding more backbones/glycans without
  better separation adds wrong competitors; full glycans alone crashed 253 -> 119
  @1%.

For the 95 generated-but-outranked cases (corrected numeric audit; the earlier 57/49
split was also a parser artifact):

| reason | count | note |
|---|---:|---|
| truth loses y-ladder | **92** | a wrong (usually implausible) mass-split has a stronger *summed* Y-ladder intensity, even though truth typically has MORE core-Y HITS |
| y-ladder tie loses rank | 3 | glyco evidence tied; sparse b/y rank picks wrong |

**The selector defect is now unambiguous:** collapse is **summed-Y-ladder-INTENSITY-primary**
(b/y rank only a tiebreak), so 92/95 outranked truths lose to a competitor with higher
summed Y intensity despite having more real core-Y hits (e.g. scan 10156: truth coreY=5
loses to winner coreY=1 whose backbone is 534 Da off). Leg-2 fix = count core-Y **hits**
(not raw summed intensity) + peptide-dominant b/y (a glyco search engine 0.65/0.35) + mass-split penalty.

Median candidate count is roughly the same across buckets (~114-118), so the miss is
not simply "too few candidates per scan"; it is **which candidate surface survives**
and **which candidate wins top-1**.

---

## 3. Refuted — DO NOT RETRY (each was a clean A/B or code-verified)

| lever | result |
|---|---|
| loosen FDR threshold | +6 correct, +200 decoys — noise, not IDs |
| full 2510-glycan list | crashes @1% 253→119 (expansion without separation) |
| `charge-expand` (try z+1..) | worse (reported charge is *correct*, not misassigned) |
| `MIN_BY` b/y quorum 6→1 | byte-identical (peptide-first gate is not the bottleneck) |
| isotope-aware backbone voting | z≥5 unchanged 16/18; near-misses were coincidental |
| RT as a standalone lever | weak seed model; @1% flat |
| wider isotope range (−1..4) | no effect |
| current `ANDES_GLYCO_SELECTOR=combined` | honest @1% worse: 253/97/1 -> 236/86/1 |
| glyco-only b/y rank retrain (SP-B) | worse under controlled/honest test |
| two-pass Percolator re-collapse | worse; TD labels cannot learn within-scan backbone correctness |
| unified glycan-decoy Percolator pile | crashed; glycan-axis features too underpowered in one label pile |

Also do not rely on the old VM `gen_audit3.py` result as-is: it double-counts
explicit `C+57.02146` when parsing PIN peptides and reports a false low
`andes_has`. Use `glyco_outrank_audit.py` and numeric mass matching instead.

---

## 4. Field map — what the reference engine/a glyco search engine/a cross-spectrum glyco engine do differently

The useful clean-room lesson is consistent across tools:

| tool | relevant winning idea | andes gap |
|---|---|---|
| the reference glyco engine | peptide-first mass-offset search; score naked b/y + glyco-aware ions with a hyperscore-like selector | andes has strong peptide scoring, but not a robust mass-offset DB enumeration branch feeding a glyco-aware selector |
| a glyco search engine | glycan-first Y-complement indexing, then separate `ScoreG`/`ScoreP`/`ScoreGP` and glycan-level QC | andes has a Y-index and glyco features, but lacks a stored two-axis fused selector |
| a cross-spectrum glyco engine | shared-backbone/cross-spectrum evidence for sparse b/y spectra | andes scaffold exists, but prior transfer payoff is bounded and cannot be the main lever |
| O-Pair/an open-source glyco engine | paired HCD/EThcD localization and graph logic | useful direction for future paired data, not the immediate HCD-only Fc3 gap |

Conclusion: beating the reference engine on this run requires **candidate-surface + selector
coupling**, not a single transfer/FDR tweak.

---

## 5. Forward plan — prioritized, audited

### P0 — Measurement and candidate-level audit (do first, cheap)
Before building, make the candidate-level truth/winner table explicit for the
~106 generated-but-outranked scans plus a sample of `truth_absent` scans:
- wrong winner peptide/glycan/backbone vs truth peptide/glycan/backbone;
- `ScoreP` candidate terms (`rank_score`, edge, b/y coverage/contiguity);
- `ScoreG` candidate terms (`YLadderScore`, `CoreYHits`, Y0/Y1 anchor, sialic);
- mass-split delta (wrong short-backbone/big-glycan patterns);
- whether the truth was excluded by top-k retention, shortlist K=24, glycan list,
  or no precursor-mass DB source.

This decides whether the first patch is a small selector/retention fix or a real
learned/listwise scorer. It also prevents repeating the current `combined` trap:
top-1 can rise while honest @1% falls.

### P1 — Build a real two-axis selector surface
Implement a stored fused score, not a post-collapse PIN-only feature:
- `ScoreP`: peptide-backbone evidence from b/y rank, edge, and preferably an
  additive coverage/contiguity term (hyperscore-like, no `score_psm` rewrite).
- `ScoreG`: glycan-Y evidence from `YLadderScore`, `CoreYHits`, Y0/Y1 anchor, and
  composition-specific terms that actually differ between competitors.
- `ScoreGP`: fixed seed fusion first (a glyco search engine-style starting point: peptide-heavy,
  e.g. ~0.65 peptide / ~0.35 glycan after scale normalization), then learn/tune only
  after the candidate audit shows the axes separate truth from real competitors.

Hard requirements:
- Store the fused scalar on `GlycoPsmKey` or equivalent so driver and PIN writer
  do not recompute different winners.
- Wire the same comparator into both `glyco_search.rs` pre-feature collapse and
  `glyco_pin.rs::select_emitted_hits`.
- Remove or audit the `SELECTOR_SHORTLIST_K=24` shortcut for glyco reranking; a
  weak-b/y truth cannot be rescued if it is outside a bare-rank shortlist.
- Keep selector and retention toggles independent. `ANDES_GLYCO_SELECTOR=...`
  must not silently change `ANDES_GLYCO_YINDEX`.

Gate: top1 audit improves AND honest Percolator @1% improves vs 253/97/1 with decoys
controlled. Top1 alone is not sufficient.

### P2 — Add precursor-mass DB enumeration, but only behind guardrails
The code already has the primitive (`bucket_index`, `db_branch`, mass->peptide lookup),
but the glyco driver still lacks a clean source that enumerates:

```text
precursor_neutral - glycan.mass -> backbone mass -> bucket_index peptides
```

Add this as an opt-in source (for example `ANDES_GLYCO_PRECURSOR_MASS=1`) with:
- sequon filter;
- charge/isotope consistency;
- bounded emission cap;
- glycan-decoy/selector guardrails before full-glycan expansion is trusted.

Do NOT ship precursor-mass expansion by itself. Prior full-list expansion increased
competition and crashed @1%. Generation only pays after P1 can rank the recovered
backbone.

### P3 — FDR/decoys only after stronger candidate scoring exists
Current glycan-axis decoys are implemented but underpowered: glycan-decoy rows share
nearly all peptide features, and only a small set of glycan features differ. A
separate glycan axis is still the field-correct shape, but not a quick fix until
there are richer composition/Y-ladder features.

Near-term FDR rule: one top-1 row per scan, Percolator only, validate with target and
decoy counts plus `glyco_recovery_fdr.py`. Do not use all-hit PINs for FDR except as
diagnostic pass-A inputs.

### P4 — Transfer and RT stay secondary
Cross-spectrum transfer is real but bounded by the measured sibling availability and
did not move the current baseline meaningfully. RT features are built and safe, but
the current seed model is flat. Revisit both after P1/P2 create a candidate surface
where orthogonal evidence can convert IDs.

---

## 6. Concrete first tasks for the next session

1. **Write/run the candidate-level competitor audit** on
   `andes_cap1024_truth_allhits.glyco.pin`: real wrong winner vs truth for the 106
   outranked scans, with the `ScoreP`/`ScoreG` columns above. Include a `truth_absent`
   sample to distinguish no-DB-source vs retention/top-k loss.
2. **Patch the experiment toggles before more A/Bs:** decouple combined selector from
   implicit `YINDEX`, and add a regression that driver collapse and PIN collapse pick
   the same winner under any selector mode.
3. **Prototype P1 stored fused selector** with no learned model first. Gate on the
   all-hit audit and then honest @1% Percolator. If it cannot improve 253/97/1, do
   not proceed to larger generation expansion.
4. **Only then add P2 precursor-mass enumeration** and test joint P1+P2. Never judge
   P2 alone on target count.

---

## 7. Key files, measuring sticks, gotchas

- **Measuring stick:** `glyco_outrank_audit.py --truth truth_nglycan_residue.tsv
  --pin <all-hits PIN> [--out per-scan.tsv]`. Categories: top1_correct /
  truth_outranked / truth_absent. Use `glyco_recovery_fdr.py <truth> <psms> <q> <tol>`
  for @1% backbone-correct.
- **Honest @1% recipe:** `andes --glyco` -> Percolator (`--seed 42 --only-psms`).
  Current baseline 253/97/1. `--mods glyco_mods.txt` should be byte/near-byte
  hygiene only because default code already includes Prot-N-term Acetyl; verify if
  changing binaries.
- **andes mass convention (bit me 3×):** backbone is **residue mass** (no water);
  the PIN writes Cam-C explicitly as `C+57.02146`. Match numerically via
  `CalcMass − GlycanMass = peptide neutral`, NOT by re-parsing the peptide string.
- **Constraints:** FDR = Percolator only (never Mokapot); additive PIN features
  only; deterministic (no HashMap in output paths — a 40% FDR swing came from one);
  model/GBDT changes engine-wide, not glyco-only; validate at @1%+decoys not top1.
- **VM:** `/srv/data/msgf-bench/glyco_bench`; source at `/srv/data/msgf-bench/andes-src`
  (plain synced copy, NOT git — scp individual files, watch zsh word-splitting).
  ~100 artifact PINs/psms accumulated (safe to leave or clean).
- **Memory index:** `[[glyco-hybrid-campaign]]` in the andes-workspace memory has the
  full detail of every thread above.

---

## 8. Housekeeping status (this session end)

- Determinism previously proven; tracked code not edited by this handoff update.
- Untracked roadmap/audit files left as-is (see section 1).
- Memory (`project_glyco_hybrid_campaign.md` + `MEMORY.md` index) fully updated.
- Nothing running on the VM.
