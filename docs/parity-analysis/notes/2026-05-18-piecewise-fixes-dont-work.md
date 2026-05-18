# Empirical study: scoring + feature fixes regress production Astral; retention layer untested

_2026-05-15 → 2026-05-18 investigation. Reverted to b1d45bb baseline; all attempted fixes archived in git reflog._

> **2026-05-18 follow-up review** (post-revert) flagged that the 2026-05-16 audit and this
> note's original framing missed two HIGH-priority retention-layer divergences. See the
> "Open: retention layer (added 2026-05-18)" section below. The TL;DR has been softened
> accordingly — what we PROVED is that piecewise fixes to scoring + Percolator features
> alone regress Astral. We never tested fixes to retention (TopNQueue tie semantics +
> per-charge GF compute), and the post-review analysis suggests they may be upstream
> prerequisites for any feature-level fix to land cleanly.

## TL;DR

Over four iterations between 2026-05-15 and 2026-05-18, we attempted to fix divergences
between msgf-rust and Java MS-GF+ that were measurable on `gf_java_parity` (5 BSA fixture
PSMs, SP-vs-SP). A comprehensive audit identified 17 divergences across param loading,
scoring pipeline, Percolator features, and search algorithm. **The audit missed two
HIGH-priority retention-layer divergences** (TopNQueue strict-greater eviction vs Java's
tie keeping; single-`top_charge` GF context vs Java's per-(spec, charge) GF). These are
documented in the "Open" section below.

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

1. **Piecewise fixes to scoring + Percolator features REGRESS production Astral.**
   This is empirically established by the iter3 → iter5 arc. Five audit-driven fixes,
   correct per Java source reading, cumulatively lost 8K PSMs on Astral. The compensating-bug
   pattern is real: fixing element X without aligning element Y exposes Y's contribution
   to the previous lucky equilibrium.

2. **The b1d45bb baseline IS a compensating-bug equilibrium**, not a correct
   implementation. `gf_java_parity` passes at 1.0 OOM (0.07-0.81 OOM range) by
   coincidence: wrong score × wrong distribution × wrong Percolator features happen to
   discriminate target/decoy on Astral about as well as Java's correct pipeline does.

## What this does NOT demonstrate

The earlier framing of this note ("piecewise alignment doesn't work, full rewrite needed")
was **overclaiming**. A 2026-05-18 follow-up code review (incorporated below) flagged two
HIGH-priority retention-layer divergences the audit missed. We never tested any retention-layer
fix. The empirical regression we observed is for piecewise fixes to **scoring** and
**Percolator features** only — not retention.

The retention layer (Java keeps tied PSMs; Rust drops them; Java computes per-charge GF;
Rust uses one `top_charge` GF for the whole queue) is **upstream** of features. If retention
is broken, no amount of feature accuracy can rescue Percolator — the right PSMs aren't being
fed to it in the first place. The original audit + the iter3-5 arc never addressed this
layer, so we cannot conclude that fixing scoring + features cleanly (with retention also
fixed) would still regress.

Honest verdict: **piecewise fixes to scoring + features alone regress production. The
retention layer is the most-likely upstream prerequisite and was never tested.**

## Open: retention layer (added 2026-05-18 — audit and iter3-5 work missed this entirely)

Independent code review on 2026-05-18 (two passes) verified the following four
retention-layer + post-merge divergences. All were missed by the 2026-05-16 audit and
were never targeted by any iter3-5 fix. They are LIKELY bug-class, with line-by-line
Java/Rust verification:

### R-1 (HIGH): TopNQueue strict-greater eviction drops ties Java keeps

**Rust** [`psm.rs:163-173`](rust/crates/search/src/psm.rs#L163-L173): `m.cmp(top) == Ordering::Greater`
— strictly greater. With capacity=1 and a tied score, the new PSM is NOT inserted.

**Java** has THREE places where it keeps tied PSMs that Rust drops:

1. [`DBScanner.java:540`](src/main/java/edu/ucsd/msjava/msdbsearch/DBScanner.java#L540) —
   per-SpecKey raw-score retention: `if (queue.size() < n || score == queue.peek().getScore())
   queue.add(...)`. Tied scores added EVEN at capacity.
2. [`DBScanner.java:719-733`](src/main/java/edu/ucsd/msjava/msdbsearch/DBScanner.java#L719-L733) —
   pre-merge dedup keyed by `pepSeq + score`: same peptide-and-score is deduped (protein
   indices merged into one match), but DIFFERENT peptides tying at the same score are NOT
   deduped — both flow to the per-spectrum merge.
3. [`DBScanner.java:745`](src/main/java/edu/ucsd/msjava/msdbsearch/DBScanner.java#L745) —
   per-spectrum SpecE merge: `if (queue.size() < n || curEValue == queue.peek().getSpecEValue())
   queue.add(...)`. Tied SpecE values added EVEN at capacity.

**Effect on `-n 1`:** Java emits multiple tied PSMs per spectrum at every tie point.
Rust emits one. Raw PSM count gap on the no-mods Astral bench was Java 89,479 / Rust 75,457
targets — a **14,022-PSM gap before Percolator even sees the data**. Plausibly explained by
this divergence.

**Audit had this as "Ord-1: dormant for top_n=1" — that conclusion was wrong.** Top_n=1
with ties produces N>1 retained PSMs in Java, not 1.

### R-2 (HIGH): Per-charge GF computation — Rust uses one `top_charge` for the whole spectrum

**Rust** [`match_engine.rs:325-340`](rust/crates/search/src/match_engine.rs#L325-L340):

```rust
let top_charge = queue.iter_psms().max_by(...).map(|p| p.charge_used)...
let scored_spec_for_gf = scored_spec_for_charge(top_charge);
compute_spec_e_values_for_spectrum(... scored_spec_for_gf, top_charge, ...);
```

ONE GF context (`top_charge`) is used for every PSM in the queue. Multi-charge candidates
share one GF distribution.

**Java** [`DBScanner.java:779`](src/main/java/edu/ucsd/msjava/msdbsearch/DBScanner.java#L779):

```java
NewRankScorer scorer = specScanner.getRankScorer(new SpecKey(specIndex, match.getCharge()));
```

Per-match charge lookup — every PSM gets SpecE computed under ITS OWN charge's GF
distribution. SpecKey indexing means each (spectrum, charge) has its own scorer +
ScoredSpectrum + GF distribution; merge happens AFTER per-SpecKey SpecE is computed.

**Interaction with R-1:** For `-n 1` without ties, Rust's queue has 1 PSM whose charge
IS `top_charge`, so R-2 is dormant. But fix R-1 (Rust keeps ties), and PSMs at different
charges land in the queue → R-2's per-charge bug ACTIVATES. R-1 and R-2 are coupled —
fixing one without the other may not help. Need both.

### R-3 (HIGH): Rust PIN writer skips `minDeNovoScore` filtering

**Rust** [`pin.rs:251-268`](rust/crates/output/src/pin.rs#L251-L268) writes every retained
PSM and [`pin.rs:472-483`](rust/crates/output/src/pin.rs#L472-L483) computes rank-2 from
the full queue with no de-novo-score filter.

**Java** [`DirectPinWriter.java:126,130-132`](src/main/java/edu/ucsd/msjava/output/DirectPinWriter.java#L126):

```java
double rank2SpecEValue = findRank2SpecEValue(matchList, params.getMinDeNovoScore());
for (int i = matchList.size() - 1; i >= 0; --i) {
    DatabaseMatch match = matchList.get(i);
    if (match.getDeNovoScore() < params.getMinDeNovoScore()) continue;  // skip
    ...
}
```

`findRank2SpecEValue` also applies the same filter internally
([DirectPinWriter.java:266](src/main/java/edu/ucsd/msjava/output/DirectPinWriter.java#L266)).

**Effect:** Java drops low-de-novo-score PSMs from the .pin; Rust emits them. Affects both
the row count AND `lnDeltaSpecEValue` (which uses rank2 from the filtered set in Java but
the unfiltered set in Rust).

### R-4 (MEDIUM): `lnEValue` denominator length-indexing off by one (sharper than C-2)

**Rust** [`match_engine.rs:589`](rust/crates/search/src/match_engine.rs#L589):
`num_distinct_peptides_at_length(peptide.length())` — passes raw `pepLen`.

**Java** [`DirectPinWriter.java:171`](src/main/java/edu/ucsd/msjava/output/DirectPinWriter.java#L171):

```java
int numPeptides = sa.getNumDistinctPeptides(params.getEnzyme() == null ? length - 2 : length - 1);
```

where `length = match.getLength() = pepLen + 2`. So Java passes:
- `pepLen + 1` for enzymatic searches (Trypsin et al)
- `pepLen` for non-enzymatic

Rust always passes `pepLen`. For enzymatic searches, the lookup is shifted by 1 length,
which is a real number — `num_distinct_peptides_at_length(n+1)` is typically 10-20% larger
than at `n`. Makes Rust's E-values systematically smaller than Java's.

Sharper than the audit's C-2 (which noted the divergence existed but not the precise math).

### F-1 (MEDIUM): `matched_ion_ratio` denominator divergence (NEW — not in original audit)

**Rust** [`match_engine.rs:790`](rust/crates/search/src/match_engine.rs#L790):
`matched_ion_ratio: num_matched as f32 / n as f32` — divides by peptide length.

**Java** [`DirectPinWriter.java:232`](src/main/java/edu/ucsd/msjava/output/DirectPinWriter.java#L232):
`computeMatchedIonRatio(features.get("NumMatchedMainIons"), length)` where
`length = pepLen + 2`. Java divides by `pepLen + 2`.

**Different from C-5b** (`longest_y_pct` denominator), which uses `pepLen − 1` in Java
([PSMFeatureFinder.java:95](src/main/java/edu/ucsd/msjava/msdbsearch/PSMFeatureFinder.java#L95)).
So Java has THREE different denominators in the same .pin row (`n`, `n-1`, `n+2`) for
similar-looking ratios. Rust uses `n` for both `longest_y_pct` (incorrect — should be `n-1`)
and `matched_ion_ratio` (incorrect — should be `n+2`).

### Why these survived

The existing Rust parity test
[`match_engine_java_parity.rs:141`](rust/crates/search/tests/match_engine_java_parity.rs#L141)
checks **scan coverage + top-1 peptide identity** on a tiny BSA fixture. It does NOT check:

- Tied-row retention (would catch R-1)
- Charge-specific GF (would catch R-2)
- .pin row counts after `minDeNovoScore` filtering (would catch R-3)
- E-value denominator parity (would catch R-4)
- Feature-denominator parity (would catch F-1)

So the existing test gives a misleading sense of parity health. R-1 through R-4 + F-1 went
undetected for the entire iter3-5 cycle. A future iteration must add stronger gates before
trusting any Java↔Rust fix as "parity-preserving."

## The 17 verified Rust↔Java divergences from the 2026-05-16 audit

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

The 2026-05-18 follow-up review changes the recommended next steps. The retention layer
(R-1 + R-2 + R-3) was never tested and is the most-likely upstream prerequisite. Suggested
fix order (highest-leverage first), with empirical gates after each:

1. **R-1 (tie semantics).** Make `TopNQueue::push` keep tied PSMs (relax strict-greater
   eviction; mirror Java's three tie-keeping points). Bench Astral. Expected: raw target
   count moves toward Java's 89,479 (we're at 75,457). If true, this alone may account for
   most of the 14K-PSM gap.

2. **R-2 (per-charge GF).** Once R-1 lands and the queue can hold multiple charges per
   spectrum, switch `compute_spec_e_values_for_spectrum` to per-PSM-charge GF lookup. R-1
   and R-2 are coupled — fixing R-1 alone exposes R-2 if queues become multi-charge.

3. **R-3 (`minDeNovoScore` filtering).** Add the filter in `pin.rs` row write AND in
   `find_rank2_spec_e_value`. Bench .pin row counts vs Java for parity.

4. **R-4 + F-1 (denominator math).** Surgical one-liners in `match_engine.rs` for
   `lnEValue` length-indexing and `matched_ion_ratio` divisor. Trivial diff once the upstream
   layers are right.

5. **Then audit-tier fixes** (multi-charge ion features, enzN/enzC/enzInt) — these were
   the ones we tried in iter5 that regressed. Per the new analysis, they may stop
   regressing once R-1 + R-2 are in place, because Percolator would be operating on the
   correct (tied, per-charge-GF-correct) PSM set.

Before any of this:

- **Strengthen `match_engine_java_parity.rs`** to check tied-row retention, .pin row
  counts, and feature-denominator parity. The current test is too weak — it caught nothing
  for the entire iter3-5 cycle.

The original framing of three paths (full port / direct Astral tuning / per-PSM diff
harness) still applies as the strategic backstop, but the retention-first approach is now
the most concrete first move.

---

### Original "three paths" (still applicable as the broader strategic frame)

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
