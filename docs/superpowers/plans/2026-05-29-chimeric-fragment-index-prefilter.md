# Chimeric Fragment-Index Prefilter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `--chimeric` search faster than Java on PXD001819 + Astral by replacing per-spectrum brute-force scoring of thousands of mass-window candidates with a fragment-evidence prefilter that scores only a top-K, preserving the entrapment-validated PSM gains.

**Architecture:** A CSR inverted index `fragment-bin → [candidate_id]` is built once per `PreparedSearch` (only under `--chimeric`) over the enumerated `candidates`. Per chimeric spectrum, observed peaks vote for candidates sharing their fragment bins; the top-K in-window candidates by vote are fed into the **existing, unchanged** `score_psm`/GF/emission path. Off path and narrow search are untouched and bit-identical.

**Tech Stack:** Rust (crates: `search`, `scoring`, `model`, `output`, `msgf-rust`), cargo, Percolator (docker, on VM), entrapment-FDP harness (`benchmark/parity/`).

**Spec:** `docs/superpowers/specs/2026-05-29-chimeric-fragment-index-prefilter-design.md`

---

## File Structure

- **Create** `crates/search/src/fragment_index.rs` — `FragmentIndex` (CSR build + bin query) and `FragmentVoter` (per-spectrum reusable vote/top-K scratch). Single responsibility: candidate generation by fragment evidence.
- **Modify** `crates/search/src/lib.rs` — register the module.
- **Modify** `crates/search/src/match_engine.rs` — build the index in `PreparedSearch` under `--chimeric`; add a `fragment_index: Option<FragmentIndex>` field; in `run_chunk_inner` replace the per-spectrum candidate set with the voter's top-K when the index is present.
- **Modify** `crates/search/src/search_params.rs` — add `chimeric_frag_index: FragIndexMode` (auto/on/off).
- **Modify** `crates/msgf-rust/src/bin/msgf-rust.rs` — `--chimeric-frag-index {auto,on,off}` CLI flag.
- **Bench-only** (local, gitignored `benchmark/parity/`): reuse `gate_chimeric_norescore.sh`, `build_entrapment_db.py`, `compute_entrapment_fdp.py`, `rank_stratified_fdr.sh`.

---

## Task 1: FragmentIndex CSR struct + build

**Files:**
- Create: `crates/search/src/fragment_index.rs`
- Modify: `crates/search/src/lib.rs`

- [ ] **Step 1: Register the module**

In `crates/search/src/lib.rs`, after the `pub(crate) mod shared_fragment;` line add:

```rust
pub(crate) mod fragment_index;
```

- [ ] **Step 2: Write the failing test**

Create `crates/search/src/fragment_index.rs` with only the test module first:

```rust
//! Chimeric Phase-4 fragment-evidence prefilter: an inverted
//! `fragment-bin -> [candidate_id]` index used as a candidate generator under
//! `--chimeric`, so only candidates with real fragment evidence are scored.

#[cfg(test)]
mod tests {
    use super::*;
    use model::amino_acid::AminoAcid;
    use model::peptide::Peptide;
    use crate::candidate_gen::Candidate;

    fn cand(seq: &str) -> Candidate {
        let residues = seq.bytes().map(AminoAcid::from_residue_unmod).collect();
        Candidate {
            peptide: Peptide::new(residues, b'-', b'-'),
            protein_index: 0,
            start_offset_in_protein: 0,
            is_decoy: false,
            is_protein_n_term: false,
            is_protein_c_term: false,
        }
    }

    #[test]
    fn build_indexes_every_candidate_fragment() {
        let cands = vec![cand("PEPTIDEK"), cand("ACDEFGHIK")];
        let idx = FragmentIndex::build(&cands, 0.02);
        // Every charge-1 b/y fragment of every candidate must be retrievable:
        // querying the bin of a known fragment m/z returns that candidate id.
        let frags = scoring_crate::scoring::fragment_ions::predict_by_ions(&cands[0].peptide, 1..=1);
        let probe = frags[0].mz;
        let hits = idx.candidates_in_bin(probe);
        assert!(hits.contains(&0u32), "candidate 0 must be indexed at its own fragment m/z");
    }

    #[test]
    fn unknown_mz_returns_empty() {
        let cands = vec![cand("PEPTIDEK")];
        let idx = FragmentIndex::build(&cands, 0.02);
        assert!(idx.candidates_in_bin(5.0).is_empty());
        assert!(idx.candidates_in_bin(99999.0).is_empty());
    }
}
```

Verify the exact `Candidate` field set and `Peptide::new` / `AminoAcid::from_residue_unmod` signatures against the codebase before running; adjust the `cand` helper if they differ (these are the constructors used elsewhere in `search`/`model` tests).

- [ ] **Step 3: Run the test to verify it fails**

Run: `cargo test -p search --lib fragment_index 2>&1 | tail -15`
Expected: FAIL — `cannot find type FragmentIndex`.

- [ ] **Step 4: Implement `FragmentIndex::build` + `candidates_in_bin`**

Prepend to `crates/search/src/fragment_index.rs` (above the test module):

```rust
use crate::candidate_gen::Candidate;
use scoring_crate::scoring::fragment_ions::predict_by_ions;

/// Inverted index: fragment-m/z bin -> candidate ids that have a charge-1 b/y
/// fragment in that bin. CSR layout (offsets + concatenated ids) keeps it
/// compact (`u32` ids) — memory was the failure mode of the abandoned Java
/// attempt, so this stays as tight as possible.
pub(crate) struct FragmentIndex {
    bin_width: f64,
    min_mz: f64,
    n_bins: usize,
    /// `bucket_offsets[b]..bucket_offsets[b+1]` indexes into `bucket_candidates`.
    bucket_offsets: Vec<u32>,
    bucket_candidates: Vec<u32>,
}

impl FragmentIndex {
    /// Build over the full candidate set (target+decoy, mod-expanded).
    /// `bin_width` is the fragment-m/z bin in Da (caller picks ~tolerance:
    /// 0.02 for high-res, 0.5 for low-res).
    pub(crate) fn build(candidates: &[Candidate], bin_width: f64) -> Self {
        // Pass A: bounds.
        let mut min_mz = f64::INFINITY;
        let mut max_mz = f64::NEG_INFINITY;
        for c in candidates {
            for ion in predict_by_ions(&c.peptide, 1..=1) {
                if ion.mz < min_mz { min_mz = ion.mz; }
                if ion.mz > max_mz { max_mz = ion.mz; }
            }
        }
        if !min_mz.is_finite() {
            return FragmentIndex { bin_width, min_mz: 0.0, n_bins: 0,
                bucket_offsets: vec![0], bucket_candidates: Vec::new() };
        }
        let n_bins = (((max_mz - min_mz) / bin_width).floor() as usize) + 1;

        // Pass B: per-bin counts.
        let mut counts = vec![0u32; n_bins];
        let bin_of = |mz: f64| -> Option<usize> {
            if mz < min_mz { return None; }
            let b = ((mz - min_mz) / bin_width).floor() as usize;
            if b < n_bins { Some(b) } else { None }
        };
        for c in candidates {
            for ion in predict_by_ions(&c.peptide, 1..=1) {
                if let Some(b) = bin_of(ion.mz) { counts[b] += 1; }
            }
        }

        // Prefix sum -> offsets.
        let mut bucket_offsets = vec![0u32; n_bins + 1];
        let mut acc = 0u32;
        for b in 0..n_bins {
            bucket_offsets[b] = acc;
            acc += counts[b];
        }
        bucket_offsets[n_bins] = acc;

        // Pass C: fill via a moving cursor copy of offsets.
        let mut cursor: Vec<u32> = bucket_offsets[..n_bins].to_vec();
        let mut bucket_candidates = vec![0u32; acc as usize];
        for (cid, c) in candidates.iter().enumerate() {
            for ion in predict_by_ions(&c.peptide, 1..=1) {
                if let Some(b) = bin_of(ion.mz) {
                    let pos = cursor[b] as usize;
                    bucket_candidates[pos] = cid as u32;
                    cursor[b] += 1;
                }
            }
        }

        FragmentIndex { bin_width, min_mz, n_bins, bucket_offsets, bucket_candidates }
    }

    #[inline]
    fn bin_index(&self, mz: f64) -> Option<usize> {
        if self.n_bins == 0 || mz < self.min_mz { return None; }
        let b = ((mz - self.min_mz) / self.bin_width).floor() as usize;
        if b < self.n_bins { Some(b) } else { None }
    }

    /// Candidate ids whose charge-1 b/y fragment falls in the bin containing `mz`.
    /// (Callers also probe `mz ± bin_width` to cover tolerance at bin edges.)
    pub(crate) fn candidates_in_bin(&self, mz: f64) -> &[u32] {
        match self.bin_index(mz) {
            Some(b) => {
                let lo = self.bucket_offsets[b] as usize;
                let hi = self.bucket_offsets[b + 1] as usize;
                &self.bucket_candidates[lo..hi]
            }
            None => &[],
        }
    }

    /// Total indexed (fragment, candidate) entries — for memory accounting.
    pub(crate) fn n_entries(&self) -> usize { self.bucket_candidates.len() }
}
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `cargo test -p search --lib fragment_index 2>&1 | tail -15`
Expected: PASS (2 tests).

- [ ] **Step 6: Commit**

```bash
git add crates/search/src/fragment_index.rs crates/search/src/lib.rs
git commit -m "feat(chimeric): FragmentIndex CSR inverted fragment->candidate index (P1)"
```

---

## Task 2: FragmentVoter — per-spectrum vote + top-K

**Files:**
- Modify: `crates/search/src/fragment_index.rs`

- [ ] **Step 1: Write the failing test**

Add to the `tests` module in `crates/search/src/fragment_index.rs`:

```rust
#[test]
fn voter_ranks_candidate_with_most_matched_fragments_first() {
    // B shares 3 fragments with the observed peaks; A shares 1. B must rank above A.
    let cands = vec![cand("PEPTIDEK"), cand("ACDEFGHIK")];
    let idx = FragmentIndex::build(&cands, 0.02);
    let b_frags = predict_by_ions(&cands[1].peptide, 1..=1);
    let a_frags = predict_by_ions(&cands[0].peptide, 1..=1);
    // observed peaks (rank, mz): 3 of B's fragments + 1 of A's.
    let peaks = vec![
        (1u32, b_frags[0].mz), (2, b_frags[1].mz), (3, b_frags[2].mz), (4, a_frags[0].mz),
    ];
    let mut voter = FragmentVoter::new(cands.len());
    // in_window = both candidates eligible.
    let topk = voter.top_k(&idx, &peaks, |_cid| true, 2);
    assert_eq!(topk[0], 1u32, "candidate B (3 matches) ranks first");
    assert!(topk.contains(&0u32));
}

#[test]
fn voter_excludes_out_of_window_candidates() {
    let cands = vec![cand("PEPTIDEK"), cand("ACDEFGHIK")];
    let idx = FragmentIndex::build(&cands, 0.02);
    let b_frags = predict_by_ions(&cands[1].peptide, 1..=1);
    let peaks = vec![(1u32, b_frags[0].mz), (2, b_frags[1].mz)];
    let mut voter = FragmentVoter::new(cands.len());
    // window excludes candidate 1 -> it must not appear even though it has the votes.
    let topk = voter.top_k(&idx, &peaks, |cid| cid == 0, 2);
    assert!(!topk.contains(&1u32));
}

#[test]
fn voter_resets_between_calls() {
    let cands = vec![cand("PEPTIDEK"), cand("ACDEFGHIK")];
    let idx = FragmentIndex::build(&cands, 0.02);
    let b_frags = predict_by_ions(&cands[1].peptide, 1..=1);
    let mut voter = FragmentVoter::new(cands.len());
    let _ = voter.top_k(&idx, &[(1, b_frags[0].mz)], |_| true, 2);
    // second call with NO peaks must yield no votes (scratch cleared).
    let topk = voter.top_k(&idx, &[], |_| true, 2);
    assert!(topk.is_empty());
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p search --lib fragment_index::tests::voter 2>&1 | tail -15`
Expected: FAIL — `cannot find type FragmentVoter`.

- [ ] **Step 3: Implement `FragmentVoter`**

Add to `crates/search/src/fragment_index.rs` (above the tests):

```rust
/// Per-thread reusable scratch for the per-spectrum vote/top-K step. `votes` is
/// sized to the candidate count; `touched` records which entries were written so
/// reset is O(touched), never O(n_candidates) (the Sage pattern — avoids the
/// global per-spectrum allocation that OOM'd the Java attempt).
pub(crate) struct FragmentVoter {
    votes: Vec<f32>,
    touched: Vec<u32>,
}

impl FragmentVoter {
    pub(crate) fn new(n_candidates: usize) -> Self {
        FragmentVoter { votes: vec![0.0; n_candidates], touched: Vec::with_capacity(4096) }
    }

    /// Accumulate one vote per matched fragment bin and return up to `k`
    /// in-window candidate ids ranked by vote (descending; ties broken by
    /// ascending id for determinism). `peaks` are `(rank, mz)` for the active
    /// observed peaks; `in_window(cid)` gates by precursor-mass eligibility.
    /// Probes the peak's bin and both neighbours to cover ±bin_width tolerance.
    pub(crate) fn top_k<F: Fn(u32) -> bool>(
        &mut self,
        idx: &FragmentIndex,
        peaks: &[(u32, f64)],
        in_window: F,
        k: usize,
    ) -> Vec<u32> {
        // Reset prior votes.
        for &c in &self.touched { self.votes[c as usize] = 0.0; }
        self.touched.clear();

        for &(_rank, mz) in peaks {
            // weight = 1.0 (matched-fragment count). Rank-weighting is a P4
            // tuning knob; count is a strong, deterministic baseline.
            for probe in [mz - idx.bin_width, mz, mz + idx.bin_width] {
                for &cid in idx.candidates_in_bin(probe) {
                    let v = &mut self.votes[cid as usize];
                    if *v == 0.0 { self.touched.push(cid); }
                    *v += 1.0;
                }
            }
        }

        // Collect in-window touched candidates with their votes, partial-sort top-k.
        let mut scored: Vec<(f32, u32)> = self.touched.iter()
            .copied()
            .filter(|&c| in_window(c))
            .map(|c| (self.votes[c as usize], c))
            .collect();
        scored.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal)
            .then(a.1.cmp(&b.1)));
        scored.truncate(k);
        scored.into_iter().map(|(_, c)| c).collect()
    }
}
```

Note: `bin_width` is referenced via `idx.bin_width`; it is a private field of the same module, so it is accessible from `FragmentVoter`.

- [ ] **Step 4: Run to verify pass**

Run: `cargo test -p search --lib fragment_index 2>&1 | tail -15`
Expected: PASS (5 tests).

- [ ] **Step 5: Commit**

```bash
git add crates/search/src/fragment_index.rs
git commit -m "feat(chimeric): FragmentVoter per-spectrum vote + top-K with O(touched) reset (P2)"
```

---

## Task 3: Wire into the search — CLI flag, index build, hot-path swap

**Files:**
- Modify: `crates/search/src/search_params.rs`
- Modify: `crates/search/src/match_engine.rs`
- Modify: `crates/msgf-rust/src/bin/msgf-rust.rs`

- [ ] **Step 1: Add the mode enum + param**

In `crates/search/src/search_params.rs`, add near the other enums:

```rust
/// Controls the chimeric fragment-evidence prefilter (Task 4 / fragment index).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FragIndexMode {
    /// On when `--chimeric` is set (the default).
    #[default]
    Auto,
    On,
    Off,
}
```

Add a field to `SearchParams` (next to `chimeric`):

```rust
    pub chimeric_frag_index: FragIndexMode,
```

Set it in every `SearchParams` constructor/default to `FragIndexMode::Auto` (search for existing `chimeric:` initializers and add the field alongside — including test fixtures in `pin.rs`, `match_engine.rs` tests, and `output` tests; the compiler will list each missing-field site).

Helper on `SearchParams`:

```rust
    /// True when the fragment-index prefilter should be active for this search.
    pub fn frag_index_active(&self) -> bool {
        self.chimeric && self.chimeric_frag_index != FragIndexMode::Off
    }
```

- [ ] **Step 2: Run the build to find all missing-field sites**

Run: `cargo build -p search -p output -p msgf-rust 2>&1 | grep -E 'missing field|error' | head`
Expected: a list of struct-literal sites missing `chimeric_frag_index`. Add `chimeric_frag_index: FragIndexMode::Auto,` (or `Off` in pure-narrow test fixtures) at each until it compiles.

- [ ] **Step 3: Add the index field to `PreparedSearch` + build it**

In `crates/search/src/match_engine.rs`, add to `PreparedSearch`:

```rust
    /// Fragment-evidence prefilter index. `Some` only when
    /// `params.frag_index_active()`; `None` keeps the brute-force path (and the
    /// entire `--chimeric off` / narrow path) bit-identical.
    pub fragment_index: Option<fragment_index::FragmentIndex>,
```

Where `PreparedSearch` is constructed (after `candidates` + `bucket_index` are built), add:

```rust
        let fragment_index = if params.frag_index_active() {
            // bin width ~ matching tolerance: high-res 0.02 Da, low-res 0.5 Da.
            let bin_width = if scorer.param().data_type.instrument.is_high_resolution() {
                0.02
            } else {
                0.5
            };
            let fi = fragment_index::FragmentIndex::build(&candidates, bin_width);
            eprintln!("FragmentIndex: {} candidates, {} fragment entries (~{} MB)",
                candidates.len(), fi.n_entries(), fi.n_entries() * 4 / 1_000_000);
            Some(fi)
        } else {
            None
        };
```

and include `fragment_index` in the struct initializer. Add `use crate::fragment_index;` to the imports if not already present.

- [ ] **Step 4: Swap the candidate set in `run_chunk_inner` (the hot path)**

In `run_chunk_inner`, the per-spectrum loop currently iterates `window_cand_indices` (the mass-window candidates) and scores each. Gate it: when `prepared.fragment_index` is `Some` and `params.chimeric`, compute the top-K via the voter and iterate THAT set instead of the full `window_cand_indices`.

Threading: build one `FragmentVoter` per Rayon worker (e.g. via `window_cand_indices`-scope thread-local or construct inside the per-spectrum closure — construct-per-spectrum is simplest and the `vec![0.0; n]` is the cost to watch; if profiling shows it dominates, hoist to a `thread_local!`). Implementation in the closure:

```rust
            // Chimeric fragment-evidence prefilter: replace the brute-force
            // window scan with the top-K candidates by fragment vote. The set
            // fed into scoring shrinks; the scoring/emission path is unchanged.
            let cand_iter: Vec<usize> = if let Some(fi) = prepared.fragment_index.as_ref() {
                // Active observed peaks (rank, mz) for the spectrum's top charge.
                let z0 = charges_to_try[0];
                let active = scored_spec_for_charge(z0).dump_active_peaks();
                let peaks: Vec<(u32, f64)> = active.iter().map(|&(r, mz, _)| (r, mz)).collect();
                // in-window membership: candidate idx present in window_cand_indices.
                // window_cand_indices is sorted+deduped, so binary_search is O(log n).
                let in_window = |cid: u32| window_cand_indices.binary_search(&(cid as usize)).is_ok();
                let mut voter = fragment_index::FragmentVoter::new(prepared.candidates.len());
                const TOP_K: usize = 64; // tuned in P4
                voter.top_k(fi, &peaks, in_window, TOP_K)
                    .into_iter().map(|c| c as usize).collect()
            } else {
                window_cand_indices.clone()
            };
```

Then change the scoring loop to iterate `cand_iter` instead of `window_cand_indices` (`for &cand_idx in &cand_iter { ... }`). Everything downstream (per-charge queues, GF SpecE, fill_post_topn, shared-fragment competition, emission) is unchanged.

`TOP_K` is a compile-time constant for P3; P4 may promote it to a `SearchParams` field for tuning.

- [ ] **Step 5: Add the CLI flag**

In `crates/msgf-rust/src/bin/msgf-rust.rs`, add a clap arg mirroring `--precursor-cal`'s `{auto,on,off}` pattern:

```rust
    /// Chimeric fragment-index prefilter: auto (on under --chimeric), on, off.
    #[arg(long, value_name = "MODE", default_value = "auto")]
    chimeric_frag_index: String,
```

Map it into `SearchParams.chimeric_frag_index`:

```rust
    chimeric_frag_index: match cli.chimeric_frag_index.as_str() {
        "on" => FragIndexMode::On,
        "off" => FragIndexMode::Off,
        _ => FragIndexMode::Auto,
    },
```

(import `FragIndexMode` from `search::search_params`).

- [ ] **Step 6: Build + off-path bit-identity test**

Run: `cargo build -p msgf-rust 2>&1 | tail -3` → Expected: compiles.
Run: `cargo test -p search -p output 2>&1 | grep -E 'test result|FAILED'`
Expected: all pre-existing tests PASS (the `java_fixtures_load` env-fixture test may fail pre-existingly — ignore that one only).

- [ ] **Step 7: Local smoke (BSA fixture) — index path runs + still emits**

Run:
```bash
target/release/msgf-rust --spectrum test-fixtures/test.mgf --database test-fixtures/BSA.fasta \
  --output-pin /tmp/fi_on.pin --chimeric --chimeric-frag-index on 2>/tmp/fi_on.log
grep -c FragmentIndex /tmp/fi_on.log; echo "rows=$(($(wc -l < /tmp/fi_on.pin)-1))"
```
Expected: build first (`cargo build --release -p msgf-rust`); log shows the `FragmentIndex:` line; rows > 0; exit 0.

- [ ] **Step 8: Commit**

```bash
git add crates/search/src/search_params.rs crates/search/src/match_engine.rs crates/msgf-rust/src/bin/msgf-rust.rs
git commit -m "feat(chimeric): wire fragment-index prefilter into search behind --chimeric-frag-index (P3)"
```

---

## Task 4: PXD recall gate (correctness — must reproduce brute-force chimeric)

**Files:** none (bench/measurement on the VM). Reuses `benchmark/parity/` scripts.

This task is empirical. The gate: chimeric+index reproduces **≥99.5%** of brute-force chimeric PSMs @1% FDR, T/D preserved, **entrapment FDP unchanged**. If a gate fails, tune `TOP_K` / `bin_width` and re-run; if it can't be met, that is a finding — report and stop.

- [ ] **Step 1: Ship + rebuild on VM**

```bash
for f in crates/search/src/fragment_index.rs crates/search/src/match_engine.rs crates/search/src/search_params.rs crates/search/src/lib.rs crates/msgf-rust/src/bin/msgf-rust.rs; do
  scp -o ControlPath=/tmp/msgfplus-bench.sock "$f" "pride-linux-vm:/srv/data/msgf-bench/chimeric-build/$f"; done
ssh -S /tmp/msgfplus-bench.sock pride-linux-vm 'cd /srv/data/msgf-bench/chimeric-build && source /root/.cargo/env && cargo build --release -p msgf-rust 2>&1 | tail -3'
```
Expected: `Finished release`.

- [ ] **Step 2: PXD brute-force vs index — PSM count @1% + wall**

Run chimeric NO_RESCORE on PXD with `--chimeric-frag-index off` (brute) and `on` (index); Percolator each; compare @1% counts and wall. Reuse the body of `gate_chimeric_norescore.sh` for the PXD arm with the extra flag. Expected: index @1% ≥ 0.995 × brute @1%; index wall < brute wall.

- [ ] **Step 3: PXD entrapment FDP with index on**

```bash
ssh -S /tmp/msgfplus-bench.sock pride-linux-vm 'cd /srv/data/msgf-bench
RUST=chimeric-build/target/release/msgf-rust
MSGF_CHIMERIC_NO_RESCORE=1 ./$RUST --spectrum data/UPS1_5000amol_R1.mzML --database data/PXD001819_entrapment.fasta --output-pin entrap-pxd/pxd-entrap-fi.pin --mods mods-numeric.txt --enzyme-specificity fully --max-missed-cleavages 2 --min-peaks 10 --min-length 6 --max-length 40 --charge-min 2 --charge-max 4 --threads 8 --precursor-cal auto --chimeric --chimeric-frag-index on --precursor-tol-ppm 5 --isotope-error-min 0 --isotope-error-max 1 >entrap-pxd/fi.log 2>&1
bash run_percolator_docker.sh entrap-pxd/pxd-entrap-fi.pin entrap-pxd pxd-fi
python3 compute_entrapment_fdp.py entrap-pxd/pxd-fi.target.psms.txt 0.01 pxd-fi'
```
Expected: `entrapment_fraction` ≈ the brute-force 0.0034 (FDR still honest). Gate: FDP ≤ ~1.5× brute and < nominal 1%.

- [ ] **Step 4: Tune if needed, then record**

If recall < 99.5%: raise `TOP_K` (e.g. 64→128) or widen `bin_width`; rebuild; re-run Step 2. If wall regresses past brute with higher K, note the recall/speed trade. Record results in `docs/parity-analysis/notes/2026-05-29-frag-index-pxd-recall.md` and commit the note.

```bash
git add docs/parity-analysis/notes/2026-05-29-frag-index-pxd-recall.md && git commit -m "bench(chimeric): P4 PXD recall+FDP gate for fragment-index prefilter"
```

---

## Task 5: Astral speed/memory gate (the objective: beat Java wall)

**Files:** none (bench on VM).

Gate: chimeric+index Astral wall **< Java 6:18** AND index memory within budget, with @1% PSMs ≥ 0.995 × brute-force Astral chimeric (77,287).

- [ ] **Step 1: Astral brute vs index — wall + @1% + RSS**

Run the Astral arm (HCD/QExactive args from `gate_chimeric_norescore.sh`) with `--chimeric-frag-index off` then `on`, capturing `/usr/bin/time -v` wall + MaxRSS, then Percolator @1%.
Expected: index wall < 6:18 (Java) and < brute 17:04; @1% ≥ 0.995 × 77,287; MaxRSS within VM budget. The `FragmentIndex:` log line reports index entries/MB.

- [ ] **Step 2: If memory over budget — slab fallback**

If the index MB is too large, shard `FragmentIndex` by precursor-mass slab (build one sub-index per ~50 Da precursor band; query only the spectrum's band). This is a `fragment_index.rs`-local change behind the same API. Add as a follow-up task only if Step 1 shows the need; otherwise skip (YAGNI).

- [ ] **Step 3: Record**

Write `docs/parity-analysis/notes/2026-05-29-frag-index-astral-speed.md` with the brute-vs-index-vs-Java table (wall, @1%, RSS) and commit.

---

## Task 6: TMT + cross-dataset validation, then PR

**Files:** none (bench) + a closing PR.

- [ ] **Step 1: 3-dataset gate table**

Run `gate_chimeric_norescore.sh` semantics with `--chimeric-frag-index on` for all 3 datasets; build the final table vs Java (PSMs + wall): PXD, Astral, TMT. Note TMT is expected to still trail on PSMs (chimeric doesn't help TMT — Lever-2a territory); the index must not make TMT worse and should speed it up.

- [ ] **Step 2: Confirm entrapment FDP preserved (PXD + Astral)**

Re-run the entrapment harness with index on for Astral (the `entrap-astral` flow) as in Task 4 Step 3. Gate: entrapment FDP ≈ brute (both ~nominal).

- [ ] **Step 3: Final note + open PR**

Write `docs/parity-analysis/notes/2026-05-29-frag-index-gate-result.md` (does chimeric+index now beat Java on PSMs AND speed for PXD+Astral?). Then per [[merge-gate-beat-java]]: only open/merge if the full gate is met on all 3 (it won't be until TMT PSMs are solved separately) — otherwise keep the branch parked and PR open as a milestone. Commit the note; push the branch.

```bash
git add docs/parity-analysis/notes/2026-05-29-frag-index-gate-result.md && git commit -m "bench(chimeric): fragment-index gate result vs Java (P6)"
```

---

## Self-Review notes (author)

- **Spec coverage:** §1 arch → T3; §2 index → T1; §3 voter/top-K → T2 + T3; §4 recall/FDP gate + flag → T1 enum, T3 flag, T4; §5 speed/mem → T5; §6 tests → T1/T2 unit + T3 smoke + T4 recall. Phases P1–P5 → Tasks 1–6.
- **Off-path bit-identity:** enforced by `fragment_index: None` when `!frag_index_active()` (T3 Step 3) + the `cand_iter` else-branch = `window_cand_indices.clone()` (T3 Step 4).
- **Type consistency:** `FragmentIndex::{build, candidates_in_bin, n_entries, bin_width}`; `FragmentVoter::{new, top_k}`; `FragIndexMode::{Auto,On,Off}`; `SearchParams::frag_index_active`; `PreparedSearch.fragment_index` — used consistently across tasks.
- **Known caveat:** the `cand` test helper and the `dump_active_peaks`/`predict_by_ions` paths must be verified against current signatures at execution time (Task 1 Step 2 note); `TOP_K`/`bin_width` are tuning knobs resolved empirically in P4/P5.
