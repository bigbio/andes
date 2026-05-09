# Java/Rust Parity Investigation — Stratified PIN Analysis

**Status:** drafted 2026-05-09

**Goal:** Identify the source(s) of the remaining ~67% top-1 PSM disagreement
between the Rust port and Java MS-GF+ on PXD001819. Distinguish between (A) a
constant scoring offset and (B) variable per-PSM divergences that flip top-1
ranking, using only the existing PIN files (no Java instrumentation, no Rust
code changes).

## Context — what's already known

After 7 rounds of Rust-side fixes (parent_mass, peak selection,
per-partition ion enumeration, Partition Ord, mme tolerance, cleavage credit,
AA-list cache), Rust↔Java parity sits at 33.6% full match / 30.6% class
flips on full PXD001819. Score distributions are systematically shifted:

- Strong matches (Rust score 50-200): 86% Java/Rust agreement
- Medium matches (0-50): 30% agreement
- Weak matches (<0): 18-22% agreement

Mean ΔRawScore on full-match PSMs is **+42** (Java higher), with stdev ~33,
range −41 to +182. The wide spread argues against "one missing constant" —
Δ is likely a mix of a property-conditioned systematic bias and per-scan
noise. The user's brainstorm refinement made this explicit: investigate
whether Δ is property-correlated, not whether it is constant.

## Investigation strategy: two parallel tracks

### Track A — ΔRawScore decomposition on full-match scans

Goal: characterize Δ as either (i) a missing global additive, (ii)
property-conditioned systematic bias, or (iii) noise that's not systematic at
all. Each of these dictates a different next step.

### Track B — Stratified flip analysis with population denominators

Goal: identify which subpopulation of scans dominates the flips, ideally
narrow enough that we can run `msgf-trace` on a focused bucket rather than
random samples. Use lift (flip_rate / baseline_rate) to find buckets where
disagreement is concentrated, not just frequent.

Both tracks consume the same input artifacts and run as one analysis script.

## Architecture

A single Python analysis script:
`benchmark/parity/analyze_parity.py`

Inputs:
- `benchmark/results/PXD001819-parity/java.pin` (already on disk)
- `benchmark/results/PXD001819-parity/rust.pin` (already on disk)
- Optional: `benchmark/data/PXD001819/*.fasta` for protein-N/C-term feature

Output:
- Markdown report at
  `docs/superpowers/specs/2026-05-09-parity-analysis-report.md`
  (separate from this design doc; the report is generated each run and
  kept for the historical record of the investigation)

The script is read-only. Does not modify any Rust or Java code.

## Components (per single-responsibility units)

Each component is a pure function over its inputs. No global state.

1. `parse_pin(path) -> list[PsmRow]`
   - Reads PIN tab-separated file, returns rows as dicts keyed by column
     name. Already partially implemented in `flip_count.py`.

2. `peptide_features(row) -> dict`
   - Extracts:
     - `length`: from `peplen` PIN column minus 2 (Java convention adds
       flanking residues to peplen)
     - `charge`: derived from one-hot `charge2` / `charge3` / `charge4`
       PIN columns
     - `n_oxidation`: number of `M+15.99491` substrings in `Peptide`
       column (allowing one trailing decimal-digit variant); tolerant
       of `M+15.9949`, `M+15.99491`, `M+15.99492`
     - `n_carbamidomethyl`: number of `C+57.02146` substrings in
       `Peptide`; same tolerance
     - `iso_off`: from `isotope_error` PIN column
     - `last_aa`: last residue character before the trailing `.X`
       flanking marker in `Peptide`
     - `pre_aa`: first character of `Peptide` (the leading flanking
       residue, e.g. `K` in `K.PEPTIDE.D`); special values include `_`
       (protein N-term) and `M` (Met-cleaved or follows-Met)
     - `is_decoy`: derived from PIN `Label` column (-1 = decoy)
     - `score_bucket`: bucket of `RawScore` per the stratification used
       earlier this session: `very_weak` ≤ -10, `weak` (-10, 0],
       `medium` (0, 50], `strong` (50, 200], `very_strong` > 200
   - All purely from the PIN row — no spectrum or FASTA needed.

3. `protein_position(row, fasta_dict) -> str`
   - Optional. Returns "n_term" / "c_term" / "internal" based on
     `start_offset_in_protein` recovered from `SpecId` + protein lookup.
   - Skipped if FASTA missing; report notes the omission.

4. `stratify(rows, feature_fn) -> dict[bucket, stats]`
   - Groups rows by `feature_fn(row)`, computes per-bucket aggregate stats:
     count, mean Δ, median Δ, stdev Δ, agreement rate, flip rate.
   - Pure function. Tested independently.

5. `compute_lift(group_rate, base_rate) -> float`
   - `lift = group_rate / base_rate`. Lift > 1 means the bucket has more
     disagreement than the population baseline; lift > 2 is a strong signal.

6. `classify_ranking_mode(java_row, rust_row) -> str`
   - Returns one of:
     - `"raw_swap"` — Java's top-1 has higher RawScore in Java's view
     - `"spec_e_swap_only"` — RawScore order agrees but SpecE inverts
     - `"both_swap"` — fully inverted on both
     - `"agree"` — same peptide
   - Used in B1 to distinguish per-PSM scoring divergence from GF DP
     (SpecEValue) divergence.

7. `format_report(...) -> str`
   - Renders the markdown report from the computed stats. No business
     logic, just templating.

## Output report sections

### Section 1 — Population overview (sanity check)

- Counts: full-match, class-flip, same-label-different-pep
- Score histograms (text-mode, 10 buckets each) for Rust vs Java
- Distribution shift quantified (median, IQR shift)

### Section 2 — A1: ΔRawScore decomposition

- Full Δ histogram (10 buckets), median, IQR, tail count (|Δ| > 100)
- Stratified Δ table per feature: `length`, `charge`, `n_oxidation`,
  `n_carbamidomethyl`, `iso_off`, `last_aa`, `pre_aa`, `is_decoy`
- Per-feature variance contribution (η²-equivalent): `between-group SS / total SS`
- Cross-tab on top-2 features by variance contribution

**Decision rule** (A2 = "narrowed Java-source line-by-line audit on the
relevant code path", a follow-up investigation step out of scope for this
analysis):

- If any single feature has η² > 0.40 → property-conditioned bias →
  proceed to A2 narrowed to the code path that consumes that feature
- If max η² < 0.20 → Δ is mostly per-scan noise → constant offset
  hypothesis is wrong; the 67% disagreement is genuinely many small
  divergences and Track A may not have a single fix
- If 0.20 ≤ η² < 0.40 for the top feature → mixed: a small systematic
  effect plus noise

### Section 3 — B1: stratified flip analysis with denominators

- For each bucket from the same feature set:
  - `n_total`, `n_flips`, `n_full_match`
  - `flip_rate = n_flips / n_total`
  - `baseline_flip_rate = total_flips / total_scans`
  - `lift = flip_rate / baseline_flip_rate`
- Sort by lift (descending), filter to `n_total >= 50`
- Top 5 buckets reported as "investigation targets for B2"

### Section 4 — Ranking-mode breakdown

- For each flip, classify as `raw_swap` / `spec_e_swap_only` / `both_swap`
- Counts and proportions for each class
- Per-class: median |Δ_RawScore| and median Δ_lnSpecEValue
- **Decision rule:**
  - If `spec_e_swap_only` dominates → GF DP / SpecEValue computation diverges;
    investigate `compute_inner` and `add_prob_dist` paths
  - If `raw_swap` dominates → per-PSM scoring still diverges; investigate
    `score_psm` / `directional_node_score`
  - If `both_swap` dominates → both paths diverge or share a common upstream
    bug

### Section 5 — Generated recommendation

Conditional on Section 2-4 findings, the script emits a one-paragraph
"next-step recommendation" pointing at:
- Specific feature(s) that explain Δ → narrow code audit
- Specific bucket(s) for `msgf-trace` follow-up
- Specific code path(s) for ranking-mode investigation

## Data flow

```
java.pin ──┐
rust.pin ──┼─► parse_pin → PsmRow ──┐
           │                         │
fasta ─────┤                         │
(optional) │                         ▼
           │             match_by_scan ──► (java_row, rust_row) pairs
           │                                          │
           │              ┌──── full-match ───────────┤
           │              ▼                           ▼
           │           A: Δ decomposition           B: flip stratification
           │              │                           │
           │              ▼                           ▼
           │           feature stratification ◄──► same feature_fn library
           │              │                           │
           └──────────────┴────► format_report ───────┘
                                       │
                                       ▼
                          markdown report (saved + printed)
```

## Error handling

- Missing scan on either side → counted separately, excluded from matched
  pool. Reported as "Java only / Rust only" counts.
- Multi-row Java scans (Java keeps tied top-1 PSMs) → analyzed under both
  "first row only" and "any-row tolerance" semantics; both numbers reported
  side by side. Default analysis uses "first row only" for backward compat
  with `flip_count.py` baseline.
- Missing FASTA → skip the protein-position feature; report notes the
  omission in section 1.
- Malformed PIN row (column count mismatch) → log + skip; counted as parse
  failure.
- Empty bucket (n_total < 50) → excluded from lift table to avoid
  statistical artifacts; total excluded count reported.

## Testing strategy

- Smoke test on the 2000-spectrum slice pins (`/tmp/java_slice.pin`,
  `/tmp/rust_slice.pin`) before running on the full 38k-spectrum pin pair.
  Slice runs in <5s and validates the script end-to-end.
- Cross-check: the script's "full match count" must equal
  `flip_count.py`'s baseline (438 / 2000 on the slice; 12,432 / 37,089 on
  the full set). Mismatch indicates a parser bug.
- Median Δ on full-match: must be ~+42. Mismatch indicates a column-index
  bug.
- Unit tests: `peptide_features`, `stratify`, `compute_lift`,
  `classify_ranking_mode` each get a small fixture-based test in the
  script's `if __name__ == "__main__"` block under a `--self-test` flag.

## Out of scope

- No Rust code changes. Track A's eventual fixes (A2/A3) are out of scope
  for this investigation step; they will be follow-up work driven by what
  this analysis reveals.
- No Java instrumentation. The user explicitly chose Rust-side-only
  investigation.
- No spectrum-level diagnostics. `msgf-trace` is the right tool for
  individual scans; this script is the population-level prerequisite that
  decides which scans deserve a trace.

## Success criteria

The investigation is successful if the report unambiguously answers:

1. Is the Δ pattern a missing constant, property-conditioned bias, or noise?
2. Which feature bucket(s) account for the bulk of flips?
3. Is the bulk of disagreement driven by RawScore or by SpecEValue?

Each of these answers maps to a concrete next investigation step (A2, B2,
GF audit, or score_psm audit). The report's "Section 5 recommendation"
formalizes the mapping.

## Why this approach minimizes wasted effort (per refinement feedback)

- A1 first quantifies whether Δ is "constant" before assuming so. The
  user's note about Δ spread (range −41 to +182) suggested a missing
  constant is unlikely; A1 makes that explicit with η² rather than gut
  feel.
- B1 uses lift not raw flip-rate, so we don't chase buckets that are
  large but unconcentrated.
- B2 escalation is targeted on the highest-lift bucket from B1, not
  random samples — avoids the "trace 50 random scans" anti-pattern.
- Ranking-mode classification cleanly separates per-PSM scoring drift
  from GF DP / SpecEValue drift; without this, both look like generic
  "score differs" and we'd waste audits on the wrong path.
