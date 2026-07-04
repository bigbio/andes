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
- **B. Glyco fragment-intensity model** — predict expected b/y + Y-ladder +
  oxonium intensities, score by spectral angle (Prosit-style). The principled
  scoring fix the sparse per-fragment *rank* model couldn't be. Needs multi-dataset
  training scale (idea G).
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
