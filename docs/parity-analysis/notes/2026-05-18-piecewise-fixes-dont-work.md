# Empirical study: piecewise Rust↔Java fixes regress production Astral

_2026-05-15 → 2026-05-18 investigation. Reverted to b1d45bb baseline; all attempted fixes archived in git reflog._

## TL;DR

Over four iterations between 2026-05-15 and 2026-05-18, we attempted to fix divergences
between msgf-rust and Java MS-GF+ that were measurable on `gf_java_parity` (5 BSA fixture
PSMs, SP-vs-SP). A comprehensive audit identified 17 verified divergences across param
loading, scoring pipeline, Percolator features, and search algorithm.

We applied five audit-driven fixes:

| Commit | Fix | Effect on Astral 1% FDR |
|---|---|---:|
| 331a3ec | score_psm includes per-edge ion_existence + error_score + source/sink cleavage | (combined w/ next two) |
| d4b93df | Per-enzyme cleavage efficiency (0.95 → 0.99999 for Trypsin) | combined Δ = -6,922 vs b1d45bb baseline |
| 43cfca1 | Per-AA FASTA-derived probability priors | combined Δ = -6,922 vs b1d45bb |
| b6844dc | Tier 1: prob_peak − H2O, longest_y_pct denominator | -219 vs prior |
| 893f7bf | Tier 2: multi-charge ion features + enzN/enzC/enzInt | -923 vs prior |

**Net: Astral 1% FDR went from 25,224 (b1d45bb) → 17,160 (after all fixes), a 32% loss.**

Every fix was correct per line-by-line reading of Java source. Net result: large
production regression. **Reverted to b1d45bb.**

## What this empirically demonstrates

1. **The b1d45bb baseline is a happy compensating-bug equilibrium**, not a correct
   implementation. `gf_java_parity` passes at 1.0 OOM (0.07-0.81 OOM range) by
   coincidence: wrong score × wrong distribution × wrong Percolator features happen to
   discriminate target/decoy on Astral about as well as Java's correct pipeline does.

2. **Java's MS-GF+ pipeline is internally consistent.** Score scale, GF distribution
   shape, and PIN features all align. Fixing any one element in Rust toward Java
   semantics breaks the Rust-specific consistency without yet reaching Java's
   actual performance.

3. **Piecewise alignment doesn't work** for closing the Astral gap. The audit's
   tiered fix plan (Tier 1: one-liners; Tier 2: port missing Java features) was
   methodologically sound, but each tier moved production further from working,
   not closer. Even the audit's predicted "5-15% PSM uplift" from Tier 2 (multi-charge
   ion features in Percolator) turned out to be -5% in practice.

## The 17 verified Rust↔Java divergences (still real, still open)

The audit doc was preserved in git reflog (commit `b8f8f77` at HEAD@{...} of
`rust-implement` before the reset). Top items, with status after reset:

### Param loading and usage
- All 14 `Param` fields read/used identically. No divergences. (Verified line-by-line.)

### Scoring pipeline
- **B-1**: `main_ion_from_param` selects prefix ions only; Java's `determineIonTypes`
  aggregates across all segments and considers suffix ions. Dormant for HCD; relevant
  for ETD.
- **B-2**: `prob_peak` numerator uses `parent_mass`; Java uses `parent_mass − H2O`.
  One-line fix. Cascades through every edge in the GF DP. Tried as Tier 1 commit
  `b6844dc`: no production impact (the formula correction is too small).

### Percolator features (.pin)
- **C-1**: `ExpMass`/`CalcMass` neutral in Rust, charge-state in Java. Same `dm`;
  different absolute scale.
- **C-2**: `lnEValue` denominator uses peptide.length() in Rust; Java uses `length + 1`
  for enzymatic searches.
- **C-3**: `isotope_error` from `MassError::isotope_offset` (Rust); from
  `(expMass − theoMass) / ISOTOPE` (Java).
- **C-4**: `enzN`/`enzC`/`enzInt` hardcoded `0` in Rust; Java computes via OpenMS
  enzymatic-boundary rules. Tried as Tier 2 (part of commit `893f7bf`): regressed
  Astral.
- **C-5**: Rust's `compute_psm_features` walks charge-1 b/y only; Java's
  `PSMFeatureFinder` walks all ion types in the partition (multi-charge, multi-segment).
  Tried as Tier 2 (part of commit `893f7bf`): regressed Astral.
- **C-5b**: `longest_y_pct` denominator uses `n` in Rust; Java uses `(n − 1)`. Tried
  as Tier 1 commit `b6844dc`: no impact.
- **C-6**: Rust emits one `Proteins` accession; Java emits one column per match.

### Search algorithm
- All search algorithm aspects (precursor tolerance, isotope iteration, top-N for
  top_n=1, peptide length, charge range) verified parity. Candidate enumeration covers
  all 5 ModLocation contexts (the `known-divergences.md` item #5 entry stating only
  Anywhere is expanded is **out-of-date**).

### Production-only paths
- **E-1**: `-precursorCal auto` not implemented. Confirmed bench uses
  `-precursorCal off`, so this is dormant for the production benchmark.

## What we ALSO confirmed (parity verified — these are not the gap)

- GF DP body, threshold pre-pass, `add_prob_dist`, enzyme adjustment, underflow guard
  — all bit-identical to Java by line-by-line reading.
- Per-node score formula, per-edge score formula.
- Source/sink cleavage accounting in `score_psm` + `compute_cleavage_credit` is
  semantically correct (after the iter3 score_psm refactor, which we reverted).
- Precursor tolerance window, isotope error iteration inclusivity.
- Init order (`set_amino_acid_probabilities` before `register_enzyme` — possible
  in principle, but no production code currently calls `set_amino_acid_probabilities`
  in Rust).

## Iteration arc — full empirical table

| State | gf_java_parity | Astral 1% FDR | PXD001819 1% FDR | TMT 1% FDR |
|---|---|---:|---:|---:|
| b1d45bb (this baseline) | PASS 1.0 OOM (0.07-0.81) | 25,224 | 15,003 | 10,572 |
| Iter3 (cleavage_eff + AA priors + score_psm refactor) | FAIL at 1.0; 2/5 PASS, 3/5 fail 1.7-2.3 OOM | 18,302 | 15,016 | 10,548 (no proper TMT routing) |
| Iter4 (+ Tier 1 prob_peak + longest_y_pct) | unchanged from iter3 | 18,083 | 15,016 | n/m |
| Iter5 (+ Tier 2 multi-charge ions + enzN/enzC/enzInt) | unchanged from iter3 | **17,160** | 15,022 | 0 |

PXD001819 was robust to all changes (CID LowRes is forgiving). TMT requires `-protocol 4`
+ TMT mod auto-application that Rust doesn't implement (bench-config gap, not algorithmic).
Astral is where the regression manifests: HCD HighRes with multi-charge fragments and
complex peptide chemistry is where each piecewise "fix" exposed the next compensating bug.

## What to do next (NOT this iteration)

The empirical lesson is that closing the Astral gap to Java requires:

1. **Either** port Java's entire scoring + feature pipeline as a coherent unit
   (multi-day effort; "rewrite the pipeline to match Java")
2. **Or** accept that Rust has its own scoring model and **tune Rust directly against
   the Astral 1% FDR target** (treat Java as a guide, not a reference; iterate on what
   moves the production metric)
3. **Or** invest in a per-PSM Rust↔Java diff harness with both engines emitting the
   same trace format, and use it to identify the actual divergence sources — not
   formula readings of source files

Path 2 is probably the most pragmatic for production. Path 1 is the most rigorous.
Path 3 is the only way to do path 1 well without more whack-a-mole.

## Status of `gf_java_parity` test

Test tolerance: **1.0 OOM** (the b1d45bb baseline value, restored by the reset).

This passes by coincidence on the 5 BSA fixture PSMs. It is NOT a meaningful gate for
production Astral performance — the audit + empirical results both demonstrated that
moving gf_java_parity SP-vs-SP closer to Java *regresses* Astral 1% FDR. A future
iteration should either:
- Replace gf_java_parity with a test that directly tracks Astral 1% FDR delta on a
  pinned dataset, or
- Document gf_java_parity as a coincidence baseline only, not a parity gate

## Files preserved from the reverted iterations

The audit doc (`docs/parity-analysis/notes/2026-05-16-divergence-audit.md` in
commit `b8f8f77`) is recoverable from git reflog. Local backups at:
- `/tmp/audit-preserved.md`  ← full 4,000-word audit catalog
- `/tmp/iter3-note-preserved.md`  ← iter3 investigation note
- `/tmp/known-divergences-preserved.md`  ← iter3-updated known-divergences

If you reuse them, drop the "iter3 fixes are good" framing and replace with
"these fixes were tried and regressed production; recipe needs whole-pipeline
coordination."

## Commits reverted

Removed from `rust-implement` (recoverable via `git reflog`):

```
893f7bf docs(output): Tier 2 audit fixes — multi-charge ion features + enzN/enzC/enzInt
b6844dc fix(scoring): Tier 1 audit fixes — prob_peak H2O + longest_y_pct denom
b8f8f77 docs(parity): comprehensive Rust↔Java divergence audit
05b8f21 test(parity): loosen gf_java_parity to 2.5 OOM + investigation note
43cfca1 fix(scoring): per-AA probability priors from FASTA frequencies
d4b93df fix(scoring): use per-enzyme cleavage efficiencies
331a3ec fix(score-fix): Task 2 — close score_psm vs GF DP source/sink-edge gap
2693a5c diag(score-fix): Task 1 — failing test exposes score_psm vs GF DP +2 gap
```

Current branch HEAD: `b1d45bb`.
