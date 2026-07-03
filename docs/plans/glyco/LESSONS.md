# andes glyco — lessons & corrections (READ BEFORE TRUSTING ANY GLYCO NUMBER)

> Written 2026-07-03 after a session that swung from "andes recovers ~30% and
> finds 2,288 glyco-PSMs" to "done correctly, andes yields 0 glyco-PSMs @1% FDR."
> The swing was **not** the algorithm changing — it was the *measurement* being
> wrong, then corrected. Every claim below is something we got wrong and fixed;
> internalize them before running another glyco experiment.

## The arc (what happened, honestly)

1. **"Everything is fine."** Full-file A/Bs showed ~29–30% recovery of a 523-scan
   "truth" set across G1/G2/G4, and andes "reporting 2,288 glyco-PSMs @1% FDR."
2. **"Everything is flat."** Every lever (glycan-Y generation, Y0/Y1 anchor,
   cross-spectrum transfer, 2D-FDR) landed at the same ~30%. Read as a robust
   algorithmic ceiling.
3. **"It was all an artifact."** A cross-check against the *published*
   the reference engine-nglycan result + a proper top-1-per-scan collapse revealed the ~30%
   was a multi-row-PIN artifact. Measured correctly, andes yields **0 glyco-PSMs
   @1% FDR** — reproducing the original SPA2 finding that had been there all along.

## The five measurement errors (each invalidated conclusions)

### L1 — Never trust a glyco PIN that has >1 PSM row per scan
andes emits ~4.5 candidate glyco PSM rows per spectrum (different peptide/glycan
hypotheses). Percolator does **not** collapse; feeding it many correlated rows per
scan fabricates a target/decoy balance and lets ~30% of scans "pass." **Collapse
to the single top-1 PSM per scan (TDC winner) BEFORE Percolator.** Proof: g4off PIN
1,020,645 rows → 5,954 (top-1/scan) → Percolator returns **0** @1% FDR (matches
SPA2 exactly). **Rule: every glyco FDR run starts with a top-1-per-scan collapse.**

### L2 — Count identifications, not PIN rows
"2,288 glyco-PSMs @1% FDR" was 2,288 *rows*, ~4.5 per scan = ~510 unique scans.
Reporting rows as identifications overstated yield ~4.5×. **Rule: report unique
scans / unique glycopeptides, never raw PIN-row counts.**

### L3 — Mass-preserving decoys cannot test mass-conditioned features
A reversed-peptide decoy has the same residue multiset → same peptide mass → same
backbone mass. So **every mass-based glyco feature** (Y0/Y1 anchor, GlycanMass,
backbone core-Y) is *identical* for a target and its reversed decoy (measured:
anchor target-mean 0.0340 vs decoy 0.0342). Percolator learns a ~0 weight and the
feature reads "flat" — not because it lacks signal, but because the decoy defeats
it by construction. **Rule: to test a mass-conditioned feature, decoys must break
that mass (glycan-axis decoy, or entrapment), never mass-preserving reversal.**

### L4 — The "523 truth" was never ground truth
The 523-scan set is a curated backbone-mass subset. It overlaps the *published*
the reference engine-nglycan targets by only 221, and andes's own IDs by ~150–190. Three
"reference" sets (andes ~510, the reference engine-nglycan ~197–716 @1% FDR, truth 523)
mutually overlap only ~190–220 — **none is ground truth.** We also compared against
an *invalid* the reference engine baseline (labile mode, 3,217 PSMs incl. non-glyco hits like
a peptide with no Asn) before finding the deposited **proper nglycan** result
(~197 @1% FDR). **Rule: a valid reference is engine-independent + FDR-controlled
(published proper-mode result, multi-engine consensus, or entrapment) — not a
curated subset, and never a mis-configured run.**

### L5 — G1/G2/G4 verdicts are void (artifact vs artifact)
Because L1–L4 held, the A/B verdicts ("G2 anchor NO-GO," "G4 no gain," "2D-FDR
lowers recall") compared artifacts to artifacts. They tell us nothing reliable.
**Rule: any pre-2026-07-03 glyco recovery number is void; re-measure on
top-1-collapsed PINs with a valid reference before citing it.**

## What is actually TRUE after correction

- **Measured correctly (top-1/scan), andes yields 0 glyco-PSMs @1% FDR** on the
  peptide axis; the glycan axis can't hold a 1% threshold either. This is the real,
  honest state and it matches SPA2.
- The real blocker is **per-scan target/decoy separation**: on stepped-HCD, the
  b/y is sparse and the glyco features are backbone-level, so at the single-best
  per scan the true (peptide, glycan) is not separable from decoys with current
  features/decoys. A single hand-crafted feature (the anchor) provably does not fix
  it.

## What survives (durable, not artifacts)

- **Code**: G0 correctness (6-decimal masses + H2O double-count guard + DET-1),
  glycan-Y index, Y0/Y1 anchor, glycan-decoy scorer, RT-gated cross-spectrum (+
  Codex fixes). All fine as *code*; they were just measured wrong.
- **Knowledge base + roadmap** (`docs/plans/glyco/`), the standardized masses /
  notations, and the harvested collection.
- **The published the reference engine-nglycan reference** (`HCC_pool_Late_Fc3_r1.pepXML`,
  PXD025455) — the first valid external comparator we have.

## The corrected way forward (glyco continues)

1. **Rebuild the eval harness first** (before any more features): top-1-per-scan
   collapse → Percolator → count unique glyco IDs @ a chosen FDR, with entrapment
   for honesty. Ban row-counting and multi-row PINs.
2. **Re-run the SP-B kill-gate the correct way**: decoy-separated AUROC of true vs
   competing-peptide (different-mass) candidates at top-1/scan — the prior NO-GO
   was on the artifact and is void.
3. **Only then** decide whether a learned discriminative per-(peptide, glycan)
   model can achieve honest per-scan separation, or whether the glyco value
   proposition needs to change.

**Meta-lesson: in glyco, the evaluation harness is more dangerous than the search
algorithm. A wrong harness produced months-worth of confident, wrong conclusions in
a single session. Trust no glyco number whose harness you have not audited.**
