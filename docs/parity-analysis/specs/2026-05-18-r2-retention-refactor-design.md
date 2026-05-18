# R-2 retention-layer refactor — design spec

_2026-05-18. Branch `rust-implement` at HEAD `de77ea9`. Builds on R-1 (commit `fc16407`)._

## Problem

R-1 fixed the `TopNQueue` tie-keeping bug (`psm.rs:163`). Empirically on the Astral no-mods
bench, R-1 alone produced 11.6× more raw target PSMs than Java (1,042,255 vs 89,479) but
**worse coverage** by deduped metrics: 9,439 fewer unique scans, 4,580 fewer unique
peptides, 8,884 fewer distinct (scan, peptide) PSMs. See
[`docs/parity-analysis/notes/2026-05-18-r1-bench-results.md`](../notes/2026-05-18-r1-bench-results.md)
for the full data.

The root cause: Rust has only ONE of Java's FOUR retention layers, plus a missing
protein-index aggregation. R-1 fixed layer 1; the remaining gaps amplify R-1 into an
over-shoot. The four-layer + protein-aggregation sequence is documented in the same note
and labeled R-2.1 through R-2.5:

| Step | What | Java reference |
|---|---|---|
| R-2.1 | per-SpecKey raw-score queue with tie keep | `DBScanner.java:534` |
| R-2.2 | pre-merge `pepSeq + score` dedup | `DBScanner.java:719-733` |
| R-2.3 | per-SpecKey GF / SpecEValue compute | `DBScanner.java:606`, `:779` |
| R-2.4 | spectrum-level merge with SpecE tie keep | `DBScanner.java:745` |
| R-2.5 | protein-index aggregation into one PIN row | `DatabaseMatch.java:75` + `DirectPinWriter.java:237` |

## Goal

Bring Rust's Astral identification metrics within striking distance of Java by closing
all five retention-layer divergences as a coherent unit.

This is **not** an empirical hypothesis test like R-1 was. R-1's role was to verify ties
matter; that's done. R-2 is a full architectural refactor to bring retention semantics
to Java parity.

## Success criteria (4-gate spec)

After R-2 is applied, the Astral no-mods bench (threads=8, same fixture and CLI as Java's
reference bench) must satisfy ALL FOUR of the following gates simultaneously:

| Gate | Threshold | b1d45bb | post-R-1 | Java |
|---|---|---:|---:|---:|
| Raw target count | 0.8× to 1.2× Java (~72K-107K) | 75,457 | 1,042,255 ❌ | 89,479 |
| Distinct (scan, peptide) within 5% | ≥34,027 | unknown | 26,934 ❌ | 35,818 |
| T/D ratio | ≥1.85 | 1.633 | 1.965 ✓ | 1.912 |
| Percolator @ 1% FDR | ≥30,000 | 25,224 | 74,204 (artifact) | 35,818 |

Plus:
- Wall time ≤ 26 min (3× the b1d45bb baseline of ~8:36)
- All existing tests pass: `gf_java_parity`, `score_psm_pxd001819_parity`,
  `score_psm_vs_gf_dp_edge_parity`, `match_engine_java_parity` (3 tests including the new R-1
  regression guard at commit `de77ea9`)

**Rollback criterion:** if ANY of the four gates fails, revert the entire R-2 commit.
Partial retention-layer fixes regress production (proven by iter3-5 and R-1).

## Architecture — Approach A: in-place refactor of `match_engine.rs`

External API unchanged: `match_spectra` still returns `Vec<TopNQueue>` with one merged
queue per spectrum. Per-SpecKey logic is internal to `match_engine.rs`'s per-spectrum loop.

### Per-spectrum loop (revised)

```
For each spectrum:
  per_charge_queues: HashMap<u8, TopNQueue>    // NEW: keyed by charge

  // Scoring loop (R-2.1)
  for each candidate × charge:
    score_psm(...)
    per_charge_queues
      .entry(charge)
      .or_insert_with(|| TopNQueue::new(top_n_psms_per_spectrum))
      .push(psm)

  // Dedup (R-2.2)
  for (charge, queue) in &mut per_charge_queues:
    dedup_pepseq_score_in_place(queue, candidates);

  // Per-charge GF (R-2.3)
  for (charge, queue) in &mut per_charge_queues:
    let scored_spec_charge = scored_spec_for_charge(*charge);
    compute_spec_e_values_for_spectrum(queue, scored_spec_charge, *charge, ...);

  // Spectrum-level merge (R-2.4)
  let merged_queue = merge_per_charge_queues(per_charge_queues, top_n_psms_per_spectrum);

  // Return one merged queue per spectrum, API unchanged
  return merged_queue;
```

### Data model change

`PsmMatch.candidate_idx: u32` → `PsmMatch.candidate_idxs: Vec<u32>`:

```rust
pub struct PsmMatch {
    pub spectrum_idx: usize,

    /// Indices into `PreparedSearch.candidates`. Length ≥ 1; the first is the
    /// "primary" candidate (used by callers that need just one). PIN writer
    /// iterates all to emit a tab-separated `Proteins` column matching Java's
    /// `DirectPinWriter.java:237`.
    ///
    /// Populated by:
    /// 1. Initial `push` from the scoring loop: a single-element Vec.
    /// 2. Dedup (`dedup_pepseq_score_in_place`): union of all candidate_idxs
    ///    sharing the same (peptide_residue_string, rounded_score) key.
    pub candidate_idxs: Vec<u32>,

    // ... other fields unchanged ...
}

impl PsmMatch {
    /// Returns the first (primary) candidate index. Callers that need a single
    /// Candidate (most do) use this; PIN writer iterates `candidate_idxs`.
    pub fn primary_candidate_idx(&self) -> u32 { self.candidate_idxs[0] }
}
```

Callers updated:
- `match_engine.rs::compute_psm_features` and similar: use `primary_candidate_idx()`.
- `pin.rs::write_psm_row`: iterate `psm.candidate_idxs` to emit one Proteins column with
  tab-separated accessions.
- Test files using `.candidate_idx`: migrate to `.primary_candidate_idx()`.

Estimated 6-8 call sites to update.

### Dedup algorithm (R-2.2)

```rust
fn dedup_pepseq_score_in_place(queue: &mut TopNQueue, candidates: &[Candidate]) {
    // Drain into Vec, group by (peptide_residue_string, score.round() as i32),
    // collapse groups, re-push.
    let drained: Vec<PsmMatch> = queue.drain().collect();
    let mut groups: HashMap<(Vec<u8>, i32), PsmMatch> = HashMap::new();
    for mut psm in drained {
        let key = (
            candidates[psm.primary_candidate_idx() as usize].peptide
                .residues.iter().map(|aa| aa.residue).collect::<Vec<_>>(),
            psm.score.round() as i32,
        );
        groups.entry(key)
            .and_modify(|existing| {
                // Aggregate this PSM's indices into the existing entry
                existing.candidate_idxs.extend(psm.candidate_idxs.iter().copied());
            })
            .or_insert(psm);
    }
    for psm in groups.into_values() {
        queue.push(psm);
    }
}
```

Note: `TopNQueue` needs a `drain()` method (currently has `into_sorted_vec()` which
consumes self; a `drain` that takes `&mut self` and clears the heap is needed).

### Per-charge GF (R-2.3)

Existing `compute_spec_e_values_for_spectrum` is called once per non-empty per-charge
queue, with the matching `scored_spec_for_charge(charge)`. No changes to the function
itself; just multiple invocations from the new outer loop.

The `top_charge` heuristic at `match_engine.rs:325` is REMOVED — each per-charge queue
uses its own GF context.

### Spectrum-level merge (R-2.4)

```rust
fn merge_per_charge_queues(
    per_charge_queues: HashMap<u8, TopNQueue>,
    capacity: u32,
) -> TopNQueue {
    let mut merged = TopNQueue::new(capacity);
    for (_charge, queue) in per_charge_queues {
        for psm in queue.into_iter() {  // existing iter_psms or drain
            merged.push(psm);
        }
    }
    merged
}
```

Because R-1 already makes `TopNQueue::push` keep ties at capacity (and `PsmMatch::cmp`
orders by `spec_e_value` first), **the merge step automatically keeps SpecE ties** —
R-2.4 is "free" once the merge loop exists.

### Protein-index aggregation in PIN (R-2.5)

In `pin.rs::write_psm_row` (around line 462), replace the single-accession emit:

```rust
// OLD:
writeln!(writer, "\t{}\t{}", cand.peptide, ctx.accession)
```

with multi-accession iteration:

```rust
// NEW: emit Peptide column, then one tab-separated accession per merged candidate_idx
write!(writer, "\t{}", cand.peptide)?;
for &cidx in &psm.candidate_idxs {
    let acc = resolve_accession(&candidates[cidx as usize], search_index);
    write!(writer, "\t{}", acc)?;
}
writeln!(writer)?;
```

This matches Java's `DirectPinWriter.java:237`: `for (String acc : proteins.accessions)
row.append('\t').append(acc);`.

The `accession` in `RowContext` becomes redundant for the Proteins column (still useful
elsewhere as the primary accession). Consider removing if cleanup is desired; not blocking.

## Testing

### TDD: failing test first

Add to `match_engine.rs` test module (or a new `tests/r2_retention.rs` integration test if
the scenario is too big for a unit test):

```rust
#[test]
fn r2_per_charge_dedup_aggregates_protein_indices() {
    // Synthetic: build a spectrum + two Candidates with the SAME peptide
    // residues but different protein indices (mimicking shared peptide across
    // target proteins). Both should score identically. Expected behavior:
    // - per-charge queue receives both via push
    // - dedup collapses them: queue has 1 PsmMatch
    // - the PsmMatch.candidate_idxs has length 2
    // ...
}
```

The test should FAIL before R-2 is implemented (current code retains 2 separate PsmMatch
entries OR 1 after eviction; either way, candidate_idxs.len() != 2 for the merged result).

### Strengthen integration test (already partially done in commit `de77ea9`)

Update `match_engine_java_parity.rs::r1_tie_retention_active_in_production_pipeline` to
also assert:
- After R-2, per-spectrum queue size is bounded by Java's per-spectrum cardinality
  (within tolerance — exact match unlikely, but the 12× over-shoot from R-1 must close)
- Distinct (scan, peptide) count on BSA fixture matches Java's reference within 5%

### Bench protocol (Astral no-mods)

Same as R-1 (see [`docs/parity-analysis/specs/2026-05-18-r1-tie-retention-test-design.md`](2026-05-18-r1-tie-retention-test-design.md)
section "Bench protocol"). Use the same SSH ControlMaster approach + nohup pattern.

After bench:
1. Count raw targets + decoys from `.pin`
2. Run Percolator (Docker, biocontainers 3.7.1)
3. Compute distinct (scan, peptide) deduped count from Percolator's
   `target.psms.txt` (the methodology established in the R-1 results note's "honest
   comparison" section)
4. Check all four gates simultaneously

### Rollback

If any gate fails: `git revert HEAD~1` to remove the R-2 commit. The R-1 work + the
strengthened test stay; the R-2 branch can be revisited with a different approach.

## Out of scope (NOT this iteration)

- **R-3** (`minDeNovoScore` PIN filter, `pin.rs:251`)
- **R-4** (`lnEValue` denominator length-indexing, `match_engine.rs:589`)
- **F-1** (`matched_ion_ratio` denominator, `match_engine.rs:790`)
- All audit-tier scoring/feature fixes: C-1 (PIN mass scale), C-3 (isotope_error source),
  C-4 (enzN/enzC/enzInt), C-5 (multi-charge ion features), C-5b (longest_y_pct denom),
  C-6 (proteins column count — partially addressed by R-2.5 multi-accession but the
  "one column per match" form remains).
- `gf_java_parity` test tolerance changes
- Performance optimization
- Per-SpecKey architecture exposed in public API (Approach B / C from brainstorming)
- Compressed-discriminator pattern (99.2% lower Rust Percolator scores) — this is a
  separate problem that R-2 alone cannot fix; requires the feature-tier fixes above
  AFTER R-2 establishes a correct retention baseline

## Risks

| Risk | Severity | Mitigation |
|---|---|---|
| `PsmMatch` data model change misses a caller | Medium | `git grep` on `.candidate_idx`; comprehensive test run before bench |
| Per-charge GF compute changes scoring scale unexpectedly | Medium | bench measures empirically; rollback if Percolator < 30K |
| Dedup loses real distinct PSMs (same peptide, different scoring path) | Medium | unit test with synthetic data + spot-check BSA fixture row count post-dedup |
| Protein-aggregation PIN format breaks Percolator parsing | Low | spot-check PIN on BSA; Percolator's input grammar accepts tab-separated proteins per row |
| Memory: dedup HashMap grows per spectrum | Low | bounded by candidates_per_spectrum × charges (~thousands); negligible |
| Hot-path performance hit from HashMap + dedup | Medium | bench wall ≤26 min gate; if exceeded, optimize or rollback |
| TopNQueue `drain` doesn't exist yet | Low | add it; one-line method using `mem::take` |

## Estimated effort

| Phase | Estimate |
|---|---:|
| Data model migration (PsmMatch + 6-8 call sites + tests) | 1 hour |
| TopNQueue `drain` method + dedup function (new code) | 1-2 hours |
| Per-charge queue + per-charge GF in match_engine | 2 hours |
| Spectrum-level merge function | 30 min |
| PIN writer multi-accession iteration | 1 hour |
| Synthetic dedup unit test | 30 min |
| Strengthen match_engine_java_parity (deduped count gate) | 1 hour |
| Sync + bench Astral + Percolator + analysis | 2 hours |
| **Total** | **~9-10 hours, split across 2 sessions** |

## Deliverables

1. One commit (or two if a logical split helps): `fix(search): R-2 retention-layer refactor — per-SpecKey queues + dedup + per-charge GF + protein aggregation`
2. Modified files:
   - `rust/crates/search/src/psm.rs` (PsmMatch data model + TopNQueue drain)
   - `rust/crates/search/src/match_engine.rs` (per-spectrum loop refactor)
   - `rust/crates/output/src/pin.rs` (multi-accession Proteins column)
   - `rust/crates/output/src/row_context.rs` (possibly simplified)
3. New test: synthetic dedup test in match_engine or psm test module
4. Updated test: `match_engine_java_parity.rs` deduped-count gate
5. Astral bench results documented in
   `docs/parity-analysis/notes/2026-05-18-r2-bench-results.md`
6. All four gates pass simultaneously, or full rollback

## Spec → Plan transition

After the user reviews and approves this spec, invoke `superpowers:writing-plans` to
generate the task-by-task implementation plan.
