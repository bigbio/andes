# Multi-tool truth validation — Fc3_r1 (2026-07-04)

**Why:** the eval truth (523 scans, `truth_nglycan_residue.tsv`) was derived from
ONE engine (the reference engine). Every andes glyco metric (260 @1% FDR / 101
backbone-correct, and every A/B) is scored against it. A single tool's opinion is
not a gold standard. Also found: the PXD025455 *deposited* result — labelled
"a commercial glyco engine" in the project metadata — is actually **the reference engine-3.2** (`the classic hyperscore engine`
format, 12,994 glyco queries, 0 a commercial glyco engine refs), so PRIDE offers no independent
second tool. We had to RUN one.

**Second tool:** an open-source glyco engine 1.1.5 (O-Pair) N-glyco search on the exact eval run
(`HCC_pool_Late_Fc3_r1.mzML`), installed via bioconda, config
`mm_nglyco.toml` (adapted from the repo's `GlycoSearchTaskconfigNGlycoTest_Run`:
DissociationType→HCD, Protease→trypsin, 2 missed cleavages, Cam-C fixed / Ox-M
var, ±10 ppm MS1 / ±20 ppm MS2, target-only FASTA so an open-source glyco engine adds its own
reversed decoys). Running an external engine for benchmark truth is clean-room
fine (only borrowing its CODE is barred). ~56 min, 222 confident N-glyco PSMs
@1% FDR.

**Result — the truth is CORROBORATED:**

| metric | value |
|---|---|
| the reference engine truth scans | 523 |
| an open-source glyco engine confident N-glyco | 222 |
| scans co-identified by BOTH | 210 (94.6% of an open-source glyco engine's IDs) |
| backbone AGREE (incl. systematic H₂O convention offset) | **196 / 210 = 93.3%** |
| genuine backbone conflicts | 14 |
| an open-source glyco engine-ONLY (truth misses) | 12 |
| the reference engine-only | 313 |

Where two independent engines both fire, they agree on the backbone 93% of the
time. So the single-tool truth is **not** a tool-specific artifact — it is a
trustworthy reference on the overlap. (The an open-source glyco engine "peptide+water" mass vs
the reference engine "residue" mass differ by a systematic 18.0106 Da; those are agreements,
not conflicts — only 14 truly conflict.)

**Artifacts (on VM `glyco_bench/`):**
- `truth_consensus.tsv` — the 523 truth scans tagged `support` =
  `the reference engine+an open-source glyco engine` (196, the high-confidence 2-tool gold standard) /
  `the reference engine-only` / `the reference engine(backbone-conflict)`.
- `mm_out/Task1GlycoSearchTask/…nglyco.psmtsv` — an open-source glyco engine PSMs.
- `mm_nglyco.toml`, `mm_consensus.py`, `mk_consensus_truth.py` — reproducible.

## Baseline re-scored vs both truths (2026-07-04)

andes default (DDA model, honest defaults) on Fc3_r1, scored against each:

| denominator | @1% FDR | scans hit | backbone-correct |
|---|---|---|---|
| full 523 (the reference engine) | 261 | 238 (45.5%) | 101 (19.3%) |
| **196 (2-tool consensus)** | 261 | 112 (57.1%) | **51 (26.0%)** |

andes scores HIGHER against the consensus (26.0% vs 19.3% backbone-correct) —
the full 523 is dragged down by the 313 the reference engine-only scans an open-source glyco engine does
not confirm (harder / less-certain). Baseline reproduces exactly (261/101 vs 523).

**The ranking bottleneck, on trustworthy truth:** of the 112 consensus scans
andes IDs @1% FDR, only 51 carry the correct backbone → **~54% of consensus-set
IDs pick the WRONG backbone**. That (not generation) is the target for any
ranking work. 84/196 consensus scans are missed @1% FDR entirely.

**How to use going forward:** report glyco metrics against BOTH the full
the reference engine truth (523) AND the 2-tool consensus (196). The consensus is the
stricter, more trustworthy denominator for any ranking/SP-B claim. The honest
SP-B target: lift consensus backbone-correct above 51/196 (train on independent
an open-source glyco engine/consensus labels, seed geometry, measure here). Open
questions: are the 313 the reference engine-only scans real (an open-source glyco engine under-sensitivity,
likely — it found only 222 total) or truth inflation? Are the 12
an open-source glyco engine-only real IDs the truth should gain? A third tool (a glyco search engine) would
break ties, but the reference engine∩an open-source glyco engine already gives a solid consensus core.
