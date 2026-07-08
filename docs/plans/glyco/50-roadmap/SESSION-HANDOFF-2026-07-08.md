# Glyco Campaign — Session Handoff & Forward Plan (2026-07-08)

> **Purpose:** a clean resume-from-here document so the next session picks up without
> re-deriving the investigation. Test bed throughout: PXD025455 `HCC_pool_Late_Fc3_r1`
> (human serum, stepped-HCD, Q Exactive HF). 523 the reference engine-truth backbones.
> Authoritative measuring stick: `glyco_outrank_audit.py` on the VM (`/srv/data/msgf-bench/glyco_bench`).

---

## 1. Branch state — `glyco-phase1` (UNMERGED, safe to build on)

Clean working tree, deterministic (byte-identical proven), twice code-reviewed
(Codex + CodeRabbit, all findings resolved). Still **beats an open-source glyco engine: 253 @1% /
97 backbone-correct / 1 decoy vs MM ~222**.

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

## 2. The governing truth — the ~253/97 plateau is a COUPLED generation∧scoring ceiling

Decomposition of the 523 truth backbones (validated tool):

| bucket | count | meaning |
|---|---|---|
| never generated | **301** (58%) | backbone not in andes' candidate pool |
| generated but out-ranked | **106** (20%) | in the pool, but a *wrong* backbone wins the top-1 collapse |
| generated + correct | 116 (22%) | andes wins (97 survive @1% FDR) |

Established beyond doubt this session:
- **It is NOT an FDR problem.** Loosening 1%→40% moved backbone-correct only
  97→103 while decoys went 1→209. The 407 non-correct are never *emitted*.
- **The missing backbones are searchable** — median 4 core-Y rungs, 14 b/y ions,
  7 oxonium; masses/charge/PTMs verified correct (90% close as `precursor =
  backbone + in-list glycan` at the reported charge). The loss is **mechanical**.
- **Recovering generation alone does NOT move @1%.** Even when high-charge (z≥5)
  backbones are forced into the pool, `top1_correct = 0` — they land deeply
  out-ranked. Generation and scoring are coupled: fix one without the other and
  the number does not move.

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

---

## 4. Forward plan — prioritized levers

### P1 (RECOMMENDED FIRST) — Scoring / F1: peptide-b/y ⊕ glycan-Y learned fusion
The binding constraint. The ceiling analysis showed 140/171 out-ranked backbones are
**pareto-blocked** on the current collinear features (rank, anchor, intensity) — they
need a NEW orthogonal axis. F1 scores the **peptide-backbone b/y** and the
**glycan-Y ladder** as *separate* signals, then combines them with a learned fusion
(a glyco search engine-2.0 does w≈0.35 glycan / 0.65 peptide). This directly attacks the 106
out-ranked AND is the prerequisite for ANY generation recovery to convert to IDs.
- **Design doc to start from:** `glyco-scoring-roadmap.md` §6b, `deep-review-synthesis.md`.
- **Validate at @1% + decoys, NOT top1_correct** (the full-glycan trap: top1 rose
  while @1% crashed).
- Needs a brainstorm pass first (separate-score architecture + how the fusion is
  learned within the FDR/Percolator constraint — FDR stays Percolator-only).

### P2 (QUICK, do alongside P1) — Benchmark hygiene: run glyco WITH `--mods`
My glyco runs used the hardcoded default (Cam-C + Met-ox) and dropped the **N-term
acetylation** the reference engine searched. PTMs are already parameterized (same `--mods`
system as normal search). Fix = pass `--mods glyco_mods.txt` (staged on the VM;
= `docs/benchmarks/configs/mods.txt` with Cam-C + Ox-M + Prot-N-term Acetyl).
Re-establish the honest @1% baseline with the fair mod set before/with P1.
No code change.

### P3 (ONLY AFTER scoring can convert) — Generation recovery
Mechanical, but blocked by scoring until P1 lands (recovered backbones land out-ranked):
- **High-charge (z≥4): ~61 backbones.** z≥5 is 89% absent; the Y-ladder voting is
  buried by spurious high-charge bins (voted-but-truncated below top_k). Fix =
  reduce spurious bins (cap/down-weight Y-ion charge interpretations), NOT a giant
  top_k (30 min/18 scans) and NOT charge-expand.
- **Low-charge (z≤3): ~100 backbones.** Have evidence + correct mass but missed;
  cause unknown (likely the voting picking a competing backbone). Un-instrumented.

### PARKED — RT prediction
Revisit only if F1 needs an extra orthogonal axis: swap the seed model for a
trained `RtIndexModel::fit()` (or fold RT-consistency into the pre-collapse
selector, since the post-collapse `DeltaRTRank` is inert on the 1-row-per-scan PIN).

---

## 5. Concrete first tasks for the next session

1. **Re-baseline with correct mods (P2, ~40 min VM):** run the honest @1% with
   `--mods glyco_mods.txt` → the real starting number (should be ≥253/97, likely
   +a few from acetyl coverage).
2. **Brainstorm the F1 fusion design (P1):** separate peptide-b/y and glycan-Y
   scores; where the combined score enters (selector vs additive PIN feature vs
   both); how the fusion weight is set/learned; how it stays FDR-honest
   (Percolator only). → then `writing-plans` → `subagent-driven-development`.
3. Only after P1 shows it can convert out-ranked → pursue P3 high-charge generation.

---

## 6. Key files, measuring sticks, gotchas

- **Measuring stick:** `glyco_outrank_audit.py --truth truth_nglycan_residue.tsv
  --pin <all-hits PIN> [--out per-scan.tsv]`. Categories: top1_correct /
  truth_outranked / truth_absent. Use `glyco_recovery_fdr.py <truth> <psms> <q> <tol>`
  for @1% backbone-correct.
- **Honest @1% recipe:** `andes --glyco` (add `--mods glyco_mods.txt`) → Percolator
  (`--seed 42 --only-psms`). Baseline 253/97/1.
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

## 7. Housekeeping status (this session end)

- Branch clean, all suites green, release build clean, determinism proven.
- Untracked `docs/plans/glyco/scripts/` = analysis scripts (left as-is).
- Memory (`project_glyco_hybrid_campaign.md` + `MEMORY.md` index) fully updated.
- Nothing running on the VM.
