# andes-glyco 2.0 — an architecture to recover the 148 and beat the reference engine (2026-07-08)

> Grounded in this session's measured decomposition of PXD025455 Fc3_r1 (523 the reference engine
> truth). Current honest @1%: gp selector = **287 backbone-correct** (beats an open-source glyco engine
> ~222; the reference engine = 523 truth). The gap is fully attributed below.

---

## 0. The governing principle (why every past attempt failed)

The gap is **FDR-bound, not scoring-bound or generation-bound in isolation.** Three
hard experimental facts from this session define the design constraints:

1. **Expansion dilutes single-Percolator FDR.** Full 2510-glycan list crashed @1%
   253→119. More candidate retention (top-K 50→500) crashed @1% 287→196. Adding
   target-labeled hypotheses without a *matched decoy on the same axis* destroys
   Percolator's target/decoy separation.
2. **Re-ranking the SAME pool is FDR-safe.** The gp selector (`rank + K·ladder`) lifted
   221→287 with decoys 1→2 — because it adds no candidates.
3. **Weak-evidence candidates can't be separated by more search.** The 67
   retention-recoverable backbones are large/high-charge spectra with sparse b/y +
   weak core-Y; surfacing them adds noise Percolator reads as decoy-like.

**Conclusion:** the winning architecture must make comprehensive generation
**self-FDR-controlling** — every added target hypothesis matched 1:1 by a decoy on the
axis it could be wrong on. That is the crux; everything else andes already has.

## 1. The gap, fully attributed (what the architecture must recover)

| bucket | count | current cause | 2.0 lever |
|---|---:|---|---|
| gp-outranked (scored, lost) | 57 | RankScore under-rewards matched-ion COUNT | **hyperscore + gp2 count term** (§3) |
| retention-truncated | 67 | weak b/y AND weak core-Y → dropped by top-K; un-truncating dilutes | **peptide-first + group/2D-FDR + transfer** (§2,§4,§5) |
| glycan-list gap | 81 | glycan not in the 456-mass searched list; full list crashes FDR | **expanded glycan DB + glycan-decoy 2D-FDR** (§2,§4) |

## 2. Stage 1 — Dual-channel comprehensive generation

Replace fragment-gated generation with **precursor-mass-first** (the the reference engine move):

- **Peptide-first / mass-offset (primary):** for each sequon peptide, treat the glycan
  as a labile mass offset — `glycan_mass = precursor − peptide_residue` — and match
  against an **expanded glycan DB** (the full N-glycan space, ~2500+, structured as a
  mass-indexed tree like a glyco search engine's GDB). Reaches ALL of the 67 + 81 regardless of
  fragment evidence. andes primitive: `db_branch` + `bucket_index` (already exist;
  today gated behind fragment evidence + top-K).
- **Glycan-Y-first (secondary):** the existing core-Y index for strong-glycan spectra.
- Union + dedup by (peptide, glycan) within mass tolerance.

Do NOT truncate by fragment rank here — truncation is what dropped the 67. The
explosion is controlled by the FDR structure (§4), not by pre-filtering.

## 3. Stage 2 — Two-channel scoring (separable peptide vs glycan evidence)

Per candidate, compute two INDEPENDENT scores (never fuse early):

- **Peptide score `P` — a hyperscore over backbone b/y:**
  `P = log(Nb! · Ny! · (1+ΣIb) · (1+ΣIy))` on the naked peptide (glycan as offset).
  The factorial terms reward matched-ion COUNT — the exact axis where the 57
  outranked truth beat their competitor (measured: truth +4.2 NumMatchedIons, +1.3
  CoreYHits, but −5 RankScore). This is why the reference engine recovers them and andes doesn't.
  Add andes' existing coverage features (longest_b/y, matchedIonRatio, EdgeScore).
- **Glycan score `G` — core-Y ladder + oxonium:** YLadderScore (intensity), CoreYHits
  (count), Y0/Y1 anchor, sialic consistency, oxonium coverage. All exist.
- **Fused `GP = wP·P + wG·G`** for the per-scan collapse winner (a glyco search engine: wG≈0.35).
  The current gp/gp2 selector (`rank + K·ladder + J·core_y_hits`) is the seed of this;
  replace `rank` with the hyperscore `P` and learn `wP/wG` (§6).

## 4. Stage 3 — Glycan-aware 2D FDR (THE crux — makes expansion pay)

This is the single component that converts §2's comprehensive generation from a
253→119 crash into a net gain. Two decoy axes, competed separately:

- **Peptide-axis decoys:** reversed-sequence sequon peptides (mass-preserving). Exist.
- **Glycan-axis decoys:** decoy glycan compositions / shifted-interior-Y-ladder glycans
  (Y0/Y1 kept, intermediate rungs shifted so a real spectrum scores the decoy below
  the target). andes' G3 machinery (`glycan_y_intensity_decoy`, paired Label-−1 rows)
  exists but is underpowered — **strengthen it with richer glycan features** so the
  target−decoy gap is large (the current gap is too small → the "unified pile crashed").

**FDR estimation (Percolator-only compatible):** emit, for every target glyco-PSM, a
PAIRED glycan-decoy row, so the decoy set = (peptide-decoys ∪ glycan-decoys). One
Percolator run then yields a q-value covering "either axis wrong" — the native
analogue of a glyco search engine's `FDR = FDR^P + FDR^G − FDR^(P∩G)` (slightly conservative, never
optimistic). Mechanism: every added target glycan-hypothesis is matched 1:1 by a
glycan-axis decoy → expansion no longer dilutes, it self-controls. **This is what makes
the 81 glycan-gap recoverable.** Optional refinement: a thin post-Percolator glycan-decoy
competition (a group post-process, not an andes FDR engine — allowed by the FDR boundary).

## 5. Stage 4 — Weak-evidence rescue (the 67, without dilution)

Two orthogonal-evidence levers that add signal WITHOUT adding raw candidates:

- **Group-FDR by glycopeptide:** aggregate evidence across all spectra of the same
  (peptide, glycan) — multiple charges/scans/RT. A weak single spectrum passes as part
  of a confident group. Thin post-process of Percolator (allowed).
- **RT-anchored cross-spectrum transfer:** for a weak spectrum, if a co-eluting sibling
  (same backbone, different glycan) is confidently ID'd, transfer the backbone with the
  sibling's confidence as an additive feature. andes scaffold exists (gated
  experimental; was FDR-unsound — must fix the seed-decoy-label leak first).

## 6. Stage 5 — Learned fusion (replace hand-set weights)

RankSVM / LambdaMART (andes has a GBDT engine) to learn `wP`, `wG`, `K`, `J`, and the
per-feature weights on a held-out split of the truth (avoid leakage on the eval set;
fit the engine-wide ModelStore path, not a glyco-only fork).

## 7. Build order (each gate = honest @1% + decoys, numeric recovery, never top1-alone)

0. **[DONE] gp / gp2 selector** — FDR-safe re-ranking. 221→287→(gp2 validating).
1. **Hyperscore peptide score `P`** — replace RankScore in the fusion. Attacks the 57.
   FDR-safe (re-ranking). Expected: 287 → ~305–320.
2. **Strengthen glycan-axis decoy features** (richer composition/Y-ladder) so the
   glycan-decoy gap is large. Prerequisite for §4; no ID change yet, enables step 3.
3. **Paired glycan-decoy + expanded glycan DB + peptide-first generation, JOINTLY.**
   Never expansion without the decoy channel. This is the make-or-break step for the
   81 + 67. Gate hard: @1%+decoys must exceed 287.
4. **Group-FDR + RT transfer** for the residual weak-evidence 67.
5. **Learn all weights.**

## 8. Honest risk assessment

- **Highest-value, highest-risk = §4 (2D FDR).** If the strengthened glycan-decoy
  channel gives a large enough target−decoy gap, expansion pays and the 81 (+ much of
  the 67) become reachable — that is the path past 287 toward/over 523. If the glycan
  features remain too weak to separate, expansion keeps diluting and the 81 stay lost.
  This is THE experiment that decides whether beating the reference engine's raw count is feasible.
- **§1 (hyperscore) and §5 (learning) are safe, incremental** — re-ranking only,
  bankable gains on the 57 + the already-scored 375.
- **Fair-comparison caveat:** the reference engine's "523" is its own FDR call; a truly fair
  "beat" needs the same numeric backbone matching on both, and ideally an orthogonal
  ground truth (synthetic/entrapment). Re-measure the reference engine's set with
  `glyco_recovery_numeric.py` conventions before declaring victory.

**One-line summary:** andes already has comprehensive generation (`db_branch`), a
FDR-safe fused selector (gp/gp2), and a glycan-decoy scaffold (G3). The missing keystone
is a **glycan-aware 2D FDR with a strong-enough glycan-decoy channel** that lets
comprehensive generation self-control instead of dilute — plus a **count-rewarding
hyperscore** for the peptide axis. Build those two and the 148 + 57 become reachable;
without the 2D FDR, expansion will keep crashing and 287 is near the ceiling.
