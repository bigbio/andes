# andes glyco — status resume (2026-07-04)

Test bed: PXD025455 `HCC_pool_Late_Fc3_r1` (human serum NASH-HCC, stepped-HCD,
Q Exactive HF). andes numbers use honest defaults (top-1-per-scan collapse,
enumerated-only, Percolator FDR). Truth is 2-engine validated.

## Yield @1% FDR (same run)
| Engine | glyco-PSMs @1% FDR | Role |
|---|---|---|
| the reference engine (raw) | ~1,370 | truth source |
| the reference engine (curated truth) | 523 | eval "truth" set |
| **andes** | **261** | Y-ladder-first generation |
| an open-source glyco engine (O-Pair) | 222 | independent 2nd tool (we ran it) |

andes yields more @1% FDR than an open-source glyco engine.

## Truth validation (measuring stick)
| Quantity | Value |
|---|---|
| the reference engine truth scans | 523 |
| an open-source glyco engine confident N-glyco | 222 |
| Co-identified by both | 210 |
| Backbone agreement | 196/210 = 93.3% |
| Conflicts | 14 |
| an open-source glyco engine-only (truth misses) | 12 |

→ 196-scan the reference engine∩an open-source glyco engine consensus = the trustworthy gold standard.
See `40-data/multitool-truth-validation.md`.

## andes vs 196 consensus
| Metric | andes |
|---|---|
| Scans hit @1% FDR | 112/196 (57.1%) |
| Backbone-correct | 51/196 (26.0%) |
| — wrong backbone (ranking loss) | 61/112 |
| — missed @1% FDR | 84/196 |

## Levers tested (all refuted, clean controls) — see `50-roadmap/spb-design.md`
| Lever | Result | Why |
|---|---|---|
| baseline (rank_score collapse) | 51/196 | — |
| Y-ladder-primary collapse | 260→197 worse | promotes de-novo/offset backbones |
| glyco b/y rank model (seed-geo, held-out) | 51→42 | glyco b/y sparse → noisier stats |
| 2-pass Percolator re-collapse | 51→43 | TD can't rank within-scan backbones |
| glycan-axis 2D-FDR (2 decoy versions) | 2D=0 | only 2 features differ target vs glycan-decoy |
| generation | ~80% find-rate | near-ceiling |

Bottom line: andes leads an open-source glyco engine on count with validated truth; the
backbone-correctness gap is a STRUCTURAL ceiling (correct-vs-wrong-backbone signal
is in neither sparse b/y nor TD labels), not a missing quick lever.

## Ideas not yet explored (ranked by leverage)
- **B. Glyco fragment-intensity model — TESTED & CLOSED (2026-07-05, see
  `50-roadmap/glyco-intensity-model-design.md`).** Feasibility passed (glyco ions
  learnable CV~0.1, model generalizes cross-run cosine 0.989) BUT the decisive
  candidate-level test refuted it for RANKING: core-Y PATTERN-fit vs SUM = identical
  discrimination (0.679 vs 0.682), pattern rescues 1.6% of sum's failures. The
  discriminative core-Y ladder is Y0+Y1 (2 ions, already summed); higher rungs
  ~noise; oxonium is backbone-independent. A richer HCD intensity model does NOT
  beat existing YLadderScore. andes at HCD ceiling.
- **A. Tighter generation** — fewer wrong backbones/scan (stricter core-Y quorum,
  precursor-mass refinement) so rank_score has fewer ways to be wrong. Attacks the
  bottleneck at the source; never A/B'd.
- **C. Orthogonal fragmentation (EThcD/ETD)** — c/z ions preserve the backbone
  that HCD b/y suppression destroys. Root-cause fix on EThcD runs. Unexplored in andes.
- **D. Cross-spectrum / RT transfer (a cross-spectrum glyco engine style)** — propagate a
  confident peptide's backbone to its RT-linked glycoforms; rescue weak-Y spectra.
- **E. Search-space expansion** — missed-cleavage / semi-tryptic backbones andes'
  tryptic search never generates. Untested for glyco.
- **F. Richer glycan-axis features** — per-composition oxonium (Fuc 512, sialyl
  274/292/657…) each with a decoy, to give the 2D glycan axis real power.
- **G. Multi-dataset learned model at scale** — harvest PXD005411 (a glyco search engine2),
  PXD016175 (a glyco search engine2), PXD030670 (a commercial glyco engine); label via the an open-source glyco engine oracle;
  precondition for B. Harvest plan in `40-data/pride-datasets.md`.

## Reusable infra built this session
- an open-source glyco engine (O-Pair) 2nd-tool oracle: conda env `mm` + `mm_nglyco.toml` on VM.
- `train-from-search --labels <tsv>` external-label training path.
- `collapse_cmp` shared collapse comparator; dual-truth recovery (`run_reeval.sh`);
  2-pass re-collapse scaffold (`recollapse.py`).
- `truth_consensus.tsv` / `truth_196.tsv` on VM.

## Foundation audit (2026-07-05) — before building cross-spectrum transfer
Deep 4-agent audit + CodeRabbit + Codex, at the user's request ("verify masses /
ion annotations / candidate generation / glycan DB are correct before a more
complex idea").
- **Mass constants** — ALL CORRECT (max deviation 4.66e-7 Da).
- **Glycan DB** — CORRECT + complete; residue-sum mass formula (no water); all
  common human N-glycans covered (2522 full / 612 common).
- **Ion m/z annotations** — ALL CORRECT (oxonium, Y-ladder, Y0/Y1 anchor, b/y
  complement); no water/proton errors. Coverage note: Y-ladder matched at +1 only
  in 3/4 fns (multiply-charged Y unmatched — a find-rate opportunity, not a bug).
- **Candidate gen + conventions** — emitted glycan-by-subtraction masses CORRECT
  end-to-end. ONE confirmed bug FIXED (f0021d81): off-by-H2O MIN_GLYCAN gate dropped
  de-novo minimal-core (2xHexNAc 406 Da) backbones.
- **Reviews**: Codex reopened the intensity-model conclusion (its "decisive test"
  used synthetic decoys, not real competitors — PROVISIONAL now). Minor: yladder
  A/B pool b/y-censored; isotope not exact superset; decoy final-rung shift; sialic
  +0 case; --labels seed metadata.
Verdict: foundations are SOLID. Fix the H2O gate (done); the "single-spectrum
exhausted" narrative is only PROVISIONAL until the intensity model is retested on
real competitors. Cross-spectrum transfer is a sound next lever regardless.

## ★★ BREAKTHROUGH (2026-07-05) — multi-charge Y-ladder + H2O gate fix
The foundation audit (before building transfer) surfaced two fixes that BROKE the
supposed "HCD ceiling":
- off-by-H2O MIN_GLYCAN gate (f0021d81)
- multiply-charged Y-ladder matching (8463cd7a) — we were matching +1 Y-ions ONLY;
  stepped-HCD deposits 2+/3+ Y-ions (the intact glycopeptide is 2-6+).

| metric | before | after |
|---|---|---|
| @1% FDR | 261 | **319** (+22%) |
| backbone-correct vs 523 | 101 | **117** |
| consensus scans hit | 112 (57.1%) | **150 (76.5%)** |
| consensus backbone-correct | 51 | **60** |

Backbone-correct is vs INDEPENDENT truth (the reference engine 117 / 2-tool consensus 60) → real
IDs, not FDR inflation. andes 319 @1% FDR now far ahead of an open-source glyco engine's 222.
LESSON: the "single-spectrum exhausted / HCD ceiling" conclusion was WRONG — Codex
suspected it; the cause was a coverage bug, not a fundamental limit. Audit foundations
BEFORE concluding a ceiling. Several prior "refuted" ranking conclusions should be
re-examined on this new (larger, correcter) candidate pool.
