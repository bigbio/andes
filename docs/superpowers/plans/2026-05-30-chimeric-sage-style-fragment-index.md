# Chimeric Sage-style Fragment Index (Approach B) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the chimeric per-spectrum brute-force candidate scoring with a Sage-style fragment index (precursor-mass-sorted peptides + fragment-m/z-sorted buckets, dual binary-search) so wide-window chimeric search beats Java on wall time while preserving the entrapment-validated PSM gains.

**Architecture:** A `SageIndex` built once per `PreparedSearch` (only under `--chimeric`) over the enumerated candidates. Per chimeric spectrum, a precursor-mass binary search yields a contiguous candidate-index range `[pre_lo, pre_hi]`; for each observed peak, bucketed fragment-m/z binary search intersected with that range increments a score buffer sized only to the window; the top-K candidates feed the existing GF scoring/emission unchanged.

**Tech Stack:** Rust (crates: `search`, `scoring`, `model`, `msgf-rust`), cargo, Percolator+entrapment harness on VM.

**Spec:** `docs/superpowers/specs/2026-05-30-chimeric-sage-style-fragment-index-design.md`

**Reused scaffolding (already on the branch from Approach A — do NOT rebuild):** `FragIndexMode {auto,on,off}` enum + `SearchParams.chimeric_frag_index` + `frag_index_active()` (in `search_params.rs`); `--chimeric-frag-index` CLI flag (in `msgf-rust.rs`); the `cand_iter` swap point + `PreparedSearch.fragment_index` field pattern (in `match_engine.rs`). Approach A's `fragment_index.rs` (FragmentIndex/FragmentVoter) is the FAILED algorithm — this plan ADDS `sage_index.rs` and switches the wiring to it; leave `fragment_index.rs` in place (dead) unless P3 cleanly removes it.

---

## File Structure

- **Create** `crates/search/src/sage_index.rs` — `SageIndex` (build + `query`). One responsibility: Sage-style candidate generation.
- **Modify** `crates/search/src/lib.rs` — register `pub(crate) mod sage_index;`.
- **Modify** `crates/search/src/match_engine.rs` — add `PreparedSearch.sage_index: Option<SageIndex>` (built when `frag_index_active()`); in `run_chunk_inner`, when present + chimeric, set `cand_iter` from `sage_index.query(...)` instead of the FragmentVoter / window scan.
- (No change needed to `search_params.rs` / `msgf-rust.rs` — the flag exists.)

---

## Task 1: SageIndex build (mass-sorted candidates + fragment-m/z buckets)

**Files:**
- Create: `crates/search/src/sage_index.rs`
- Modify: `crates/search/src/lib.rs`

- [ ] **Step 1: Register module.** In `crates/search/src/lib.rs`, near the other module decls add: `pub(crate) mod sage_index;`

- [ ] **Step 2: Write the failing build test.** Create `crates/search/src/sage_index.rs`:

```rust
//! Chimeric Sage-style fragment index (Approach B). Peptides sorted by precursor
//! mass (so a mass window is a contiguous index range); fragments sorted by m/z in
//! fixed buckets, each bucket re-sorted by peptide index. The per-spectrum query
//! (Task 2) does a dual binary search bounded to the precursor window — the bound
//! Approach A's vote-all-touched prefilter lacked.

use crate::candidate_gen::Candidate;
use scoring_crate::scoring::fragment_ions::predict_by_ions;

/// Fixed fragment bucket size (Sage uses 8192; power-of-two).
const BUCKET: usize = 8192;

#[derive(Clone, Copy)]
struct Frag {
    mz: f32,
    /// index into `sorted_cand` / `sorted_mass` (Sage's PeptideIx).
    pidx: u32,
}

pub(crate) struct SageIndex {
    /// candidate ids ordered by ascending peptide neutral mass.
    sorted_cand: Vec<u32>,
    /// parallel to `sorted_cand`: ascending neutral masses (for precursor binary search).
    sorted_mass: Vec<f64>,
    /// all candidates' charge-1 b/y fragments; globally m/z-sorted, then each
    /// `BUCKET`-sized chunk re-sorted by `pidx`.
    fragments: Vec<Frag>,
    /// min fragment m/z per bucket (len = ceil(fragments.len()/BUCKET)).
    bucket_min_mz: Vec<f32>,
}

impl SageIndex {
    /// Build over the full candidate set (target+decoy, mod-expanded).
    pub(crate) fn build(candidates: &[Candidate]) -> Self {
        // 1. mass-sorted candidate order.
        let mut order: Vec<u32> = (0..candidates.len() as u32).collect();
        order.sort_by(|&a, &b| {
            candidates[a as usize].peptide.mass()
                .partial_cmp(&candidates[b as usize].peptide.mass())
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let sorted_mass: Vec<f64> = order.iter()
            .map(|&c| candidates[c as usize].peptide.mass())
            .collect();

        // 2. fragments at each candidate's pidx (compute predict_by_ions ONCE/candidate).
        let mut fragments: Vec<Frag> = Vec::new();
        for (pidx, &cid) in order.iter().enumerate() {
            for ion in predict_by_ions(&candidates[cid as usize].peptide, 1..=1) {
                fragments.push(Frag { mz: ion.mz as f32, pidx: pidx as u32 });
            }
        }

        // 3. global m/z sort, then per-bucket re-sort by pidx.
        fragments.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(std::cmp::Ordering::Equal));
        let mut bucket_min_mz = Vec::with_capacity(fragments.len() / BUCKET + 1);
        let mut start = 0;
        while start < fragments.len() {
            let end = (start + BUCKET).min(fragments.len());
            bucket_min_mz.push(fragments[start].mz);
            fragments[start..end].sort_by_key(|f| f.pidx);
            start = end;
        }

        SageIndex { sorted_cand: order, sorted_mass, fragments, bucket_min_mz }
    }

    /// Total indexed fragment entries (memory accounting: 8 B each).
    pub(crate) fn n_fragments(&self) -> usize { self.fragments.len() }
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::amino_acid::AminoAcid;
    use model::peptide::Peptide;

    fn cand(seq: &str) -> Candidate {
        let residues = seq.bytes().map(|b| AminoAcid::standard(b).unwrap()).collect();
        Candidate {
            peptide: Peptide::new(residues, b'-', b'-'),
            protein_index: 0, start_offset_in_protein: 0,
            is_decoy: false, is_protein_n_term: false, is_protein_c_term: false,
        }
    }

    #[test]
    fn build_sorts_by_mass_and_indexes_fragments() {
        // ACDEFGHIK is heavier than AAK -> mass order must be [AAK, ACDEFGHIK].
        let cands = vec![cand("ACDEFGHIK"), cand("AAK")];
        let idx = SageIndex::build(&cands);
        // sorted_mass ascending:
        assert!(idx.sorted_mass[0] <= idx.sorted_mass[1]);
        // lighter peptide (AAK = candidate id 1) is pidx 0:
        assert_eq!(idx.sorted_cand[0], 1u32);
        // total fragments = sum of charge-1 b/y across both candidates:
        let expect: usize = cands.iter()
            .map(|c| predict_by_ions(&c.peptide, 1..=1).len()).sum();
        assert_eq!(idx.n_fragments(), expect);
        // bucket_min_mz is non-decreasing.
        for w in idx.bucket_min_mz.windows(2) { assert!(w[0] <= w[1]); }
    }
}
```

Verify the `AminoAcid::standard` / `Peptide::new` / `Candidate` field usages against the codebase first (these match what Approach A's `fragment_index.rs` tests used — copy that exact pattern). Adjust if signatures differ.

- [ ] **Step 3: Run to verify it fails.** `cargo test -p sage_index 2>&1` won't match; run `cargo test -p search --lib sage_index 2>&1 | tail -15`. Expected: FAIL (compile error / type not found) BEFORE you paste the impl, or PASS after. (Write the test, run, see fail, then the impl above makes it pass — if you pasted both, split: comment out the impl, see the test fail to compile, restore.)

- [ ] **Step 4: Run to verify it passes.** `cargo test -p search --lib sage_index 2>&1 | tail -15` → 1 passed. `cargo build -p search` → clean. If `dead_code` warnings on `sorted_cand`/`fragments`/`bucket_min_mz`/`Frag.mz` (used only in Task 2), add a module-level `#![allow(dead_code)]` with a comment that Task 2 consumes them.

- [ ] **Step 5: Commit.**
```bash
git add crates/search/src/sage_index.rs crates/search/src/lib.rs
git commit -m "feat(chimeric): SageIndex build — mass-sorted candidates + m/z-bucketed fragments (B-P1)"
```

---

## Task 2: SageIndex query (dual binary-search, window-bounded) + local microbenchmark

**Files:**
- Modify: `crates/search/src/sage_index.rs`

- [ ] **Step 1: Write the failing query tests.** Add to the `tests` module:

```rust
fn frag_mzs(c: &Candidate) -> Vec<f64> {
    predict_by_ions(&c.peptide, 1..=1).iter().map(|i| i.mz).collect()
}

#[test]
fn query_returns_in_window_candidate_matching_peaks() {
    let cands = vec![cand("AAK"), cand("ACDEFGHIK"), cand("PEPTIDEK")];
    let idx = SageIndex::build(&cands);
    // Target candidate 2 (PEPTIDEK): observed peaks = its fragments; precursor
    // window = its mass +/- 0.01.
    let m = cands[2].peptide.mass();
    let peaks: Vec<f64> = frag_mzs(&cands[2]);
    let got = idx.query(m - 0.01, m + 0.01, &peaks, 0.02, 5);
    assert!(got.contains(&2u32), "PEPTIDEK must be returned (in window + fragments match)");
}

#[test]
fn query_excludes_out_of_precursor_window_even_if_fragments_match() {
    let cands = vec![cand("AAK"), cand("PEPTIDEK")];
    let idx = SageIndex::build(&cands);
    // Feed PEPTIDEK's fragments but a precursor window around AAK's mass only.
    let m_aak = cands[0].peptide.mass();
    let peaks = frag_mzs(&cands[1]);
    let got = idx.query(m_aak - 0.01, m_aak + 0.01, &peaks, 0.02, 5);
    assert!(!got.contains(&1u32), "PEPTIDEK is out of the precursor window -> excluded");
}

#[test]
fn query_ranks_by_matched_fragment_count() {
    let cands = vec![cand("PEPTIDEK"), cand("ACDEFGHIK")];
    let idx = SageIndex::build(&cands);
    // wide precursor window covering both; peaks = 3 of candidate 0 + 1 of candidate 1.
    let lo = idx.sorted_mass[0] - 1.0;
    let hi = idx.sorted_mass[idx.sorted_mass.len()-1] + 1.0;
    let f0 = frag_mzs(&cands[0]);
    let f1 = frag_mzs(&cands[1]);
    let peaks = vec![f0[0], f0[1], f0[2], f1[0]];
    let got = idx.query(lo, hi, &peaks, 0.02, 2);
    assert_eq!(got[0], 0u32, "candidate 0 (3 matches) ranks first");
}
```

- [ ] **Step 2: Run to verify fail.** `cargo test -p search --lib sage_index::tests::query 2>&1 | tail` → FAIL (`no method query`).

- [ ] **Step 3: Implement `query`.** Add to `impl SageIndex`:

```rust
    /// Candidate ids (top-`k` by matched-fragment count) whose precursor neutral
    /// mass is in `[mass_lo, mass_hi]` and which have charge-1 b/y fragments near
    /// the observed `peaks` (within `tol` Da). The score buffer is sized to the
    /// precursor window only — work is bounded by the window, never the whole DB.
    pub(crate) fn query(
        &self, mass_lo: f64, mass_hi: f64, peaks: &[f64], tol: f64, k: usize,
    ) -> Vec<u32> {
        // 1. precursor window -> contiguous pidx range [pre_lo, pre_hi).
        let pre_lo = self.sorted_mass.partition_point(|&m| m < mass_lo);
        let pre_hi = self.sorted_mass.partition_point(|&m| m <= mass_hi);
        if pre_hi <= pre_lo { return Vec::new(); }
        let mut scores = vec![0u16; pre_hi - pre_lo];

        // 2. per peak: buckets overlapping [mz-tol, mz+tol], intersected with pidx range.
        for &mz in peaks {
            let lo_mz = (mz - tol) as f32;
            let hi_mz = (mz + tol) as f32;
            // first bucket whose min_mz could contain lo_mz: the bucket before the
            // first whose min_mz > lo_mz.
            let b_start = self.bucket_min_mz.partition_point(|&m| m <= lo_mz).saturating_sub(1);
            let b_end = self.bucket_min_mz.partition_point(|&m| m <= hi_mz); // exclusive bucket index
            for b in b_start..b_end.max(b_start + 1) {
                if b >= self.bucket_min_mz.len() { break; }
                let f_lo = b * BUCKET;
                let f_hi = (f_lo + BUCKET).min(self.fragments.len());
                let slice = &self.fragments[f_lo..f_hi]; // sorted by pidx
                // pidx sub-range [pre_lo, pre_hi) via binary search on pidx.
                let s = slice.partition_point(|f| (f.pidx as usize) < pre_lo);
                let e = slice.partition_point(|f| (f.pidx as usize) < pre_hi);
                for f in &slice[s..e] {
                    if (f.mz - mz as f32).abs() <= tol as f32 {
                        scores[f.pidx as usize - pre_lo] += 1;
                    }
                }
            }
        }

        // 3. top-k pidx by score (desc), ties by ascending candidate id; drop zero-score.
        let mut scored: Vec<(u16, u32)> = scores.iter().enumerate()
            .filter(|(_, &s)| s > 0)
            .map(|(i, &s)| (s, self.sorted_cand[pre_lo + i]))
            .collect();
        scored.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
        scored.truncate(k);
        scored.into_iter().map(|(_, c)| c).collect()
    }
```

- [ ] **Step 4: Run to verify pass.** `cargo test -p search --lib sage_index 2>&1 | tail` → 4 passed. `cargo clippy -p search 2>&1 | tail -3` → clean.

- [ ] **Step 5: Write the LOCAL query-cost microbenchmark test (the Approach-A safeguard).** Add:

```rust
    #[test]
    fn query_is_sub_millisecond_on_a_dense_spectrum() {
        // ~2000 synthetic candidates of varied length -> a realistic index;
        // a dense 60-peak spectrum; one query must be well under 1 ms.
        let seqs = ["PEPTIDEK","ACDEFGHIK","SAMPLERK","VWXYTESTR","AAAAAAK",
                    "LLLLLLLR","GGGGTESTK","MNPQRSTK","FFFYYYWK","CCDDEEK"];
        let mut cands = Vec::new();
        for i in 0..2000 { cands.push(cand(seqs[i % seqs.len()])); }
        let idx = SageIndex::build(&cands);
        // dense spectrum: 60 peaks spanning typical fragment m/z.
        let peaks: Vec<f64> = (0..60).map(|i| 200.0 + i as f64 * 18.0).collect();
        let lo = idx.sorted_mass[0] - 5.0;
        let hi = idx.sorted_mass[idx.sorted_mass.len()-1] + 5.0; // worst case: whole window
        let t0 = std::time::Instant::now();
        let iters = 200;
        for _ in 0..iters { let _ = idx.query(lo, hi, &peaks, 0.02, 64); }
        let per = t0.elapsed().as_secs_f64() / iters as f64;
        assert!(per < 1e-3, "per-query {:.3} ms exceeds 1 ms budget", per * 1e3);
    }
```

Note: `std::time::Instant` is allowed in tests. If the assert fails, the algorithm degenerates — STOP and report (do not proceed to VM). The whole-window `lo..hi` here is the worst case; real chimeric windows are far narrower.

- [ ] **Step 6: Run the microbenchmark.** `cargo test -p search --lib sage_index::tests::query_is_sub_millisecond -- --nocapture 2>&1 | tail`. Expected: PASS. If it FAILS, report the per-query time as a finding and STOP.

- [ ] **Step 7: Commit.**
```bash
git add crates/search/src/sage_index.rs
git commit -m "feat(chimeric): SageIndex window-bounded query + local sub-ms microbenchmark (B-P2)"
```

---

## Task 3: Wire SageIndex into the search (replace the chimeric candidate gen)

**Files:**
- Modify: `crates/search/src/match_engine.rs`

- [ ] **Step 1: Add the field + build.** In `PreparedSearch` add `pub(crate) sage_index: Option<sage_index::SageIndex>,` (add `use crate::sage_index;`). Where `PreparedSearch` is built, after `candidates` exist and next to the existing `fragment_index` build, add:

```rust
        let sage_index = if params.frag_index_active() {
            let si = sage_index::SageIndex::build(&candidates);
            eprintln!("SageIndex: {} candidates, {} fragments (~{} MB)",
                candidates.len(), si.n_fragments(), si.n_fragments() * 8 / 1_000_000);
            Some(si)
        } else { None };
```
Include `sage_index` in the struct initializer. (You may set the old `fragment_index` build to `None` unconditionally now, or leave it — but the hot path below must use `sage_index`.)

- [ ] **Step 2: Swap the candidate generation in `run_chunk_inner`.** Locate the `cand_iter` block added by Approach A (the one using `FragmentVoter` / `self.fragment_index`). Replace its `Some(...)` arm so that, when `self.sage_index` is `Some(si)` and `params.chimeric`, it computes the per-charge candidate union from the Sage query. Use the SAME per-charge neutral-mass window the brute path derives (reuse `candidate_nominal_bounds` or the existing chimeric isolation-window math to get `(min_neutral, max_neutral)` per charge — match what `window_cand_indices` covers). For each charge in `charges_to_try`, get active peaks via `scored_spec_for_charge(z).dump_active_peaks().into_iter().map(|(_,mz,_)| mz)`, choose `tol = if scorer.param().data_type.instrument.is_high_resolution() { peak_mz*20e-6 } else { 0.5 }` (compute per-peak for ppm), call `si.query(min_neutral, max_neutral, &peaks, tol, TOP_K)`, and union the results across charges (dedup). Map to `Vec<usize>`:

```rust
            let cand_iter: Vec<usize> = if let (Some(si), true) =
                (self.sage_index.as_ref(), params.chimeric)
            {
                const TOP_K: usize = 64;
                let mut out: Vec<usize> = Vec::new();
                for &z in &charges_to_try {
                    let (min_nom, max_nom) = candidate_nominal_bounds(spec, z, params, shift_ppm);
                    // nominal -> neutral mass bounds for the precursor binary search:
                    let lo = min_nom as f64 / model::mass::INTEGER_MASS_SCALER; // see note
                    let hi = max_nom as f64 / model::mass::INTEGER_MASS_SCALER;
                    let ss = scored_spec_for_charge(z);
                    let high_res = scorer.param().data_type.instrument.is_high_resolution();
                    let peaks: Vec<f64> = ss.dump_active_peaks().into_iter().map(|(_,mz,_)| mz).collect();
                    // ppm tol -> use a representative per-peak tol; query takes a scalar tol,
                    // so pass 0.02 (high-res) / 0.5 (low-res) as a fixed Da window matching
                    // matched_peak_keys semantics. (If per-peak ppm is needed, widen query.)
                    let tol = if high_res { 0.02 } else { 0.5 };
                    out.extend(si.query(lo, hi, &peaks, tol, TOP_K).into_iter().map(|c| c as usize));
                }
                out.sort_unstable(); out.dedup(); out
            } else {
                window_cand_indices.clone()
            };
```

NOTE on `INTEGER_MASS_SCALER`: confirm how nominal mass maps to neutral mass in this codebase (`model::mass::nominal_from` is the forward map; find its inverse/scaler). If a clean nominal→neutral inverse isn't available, instead derive the neutral-mass window directly from `spec.precursor_mz`, the charge, and the isolation offsets the way `compute_spec_e_values_for_spectrum` does (`(precursor_mz - PROTON) * charge - H2O` ± window). Use whichever gives a neutral-mass `[lo, hi]` consistent with `window_cand_indices`. This is the one integration subtlety — verify the window matches the brute path on a sample (the recall gate in Task 4 will catch a mismatch).

Then ensure the scoring loop iterates `&cand_iter` (Approach A already changed this; keep it). Everything downstream unchanged.

- [ ] **Step 3: Build + off-path tests.** `cargo build -p msgf-rust 2>&1 | tail -3` → compiles. `cargo test -p search -p output 2>&1 | grep -E 'test result|FAILED'` → all pass except the pre-existing `java_fixtures_load::tryp_pig_bov_revcat_full_set_loads`. `cargo clippy -p search 2>&1 | tail -3` → clean.

- [ ] **Step 4: BSA smoke.**
```bash
cargo build --release -p msgf-rust 2>&1 | tail -2
target/release/msgf-rust --spectrum test-fixtures/test.mgf --database test-fixtures/BSA.fasta --output-pin /tmp/sage_off.pin --chimeric --chimeric-frag-index off 2>/dev/null; echo off=$?
target/release/msgf-rust --spectrum test-fixtures/test.mgf --database test-fixtures/BSA.fasta --output-pin /tmp/sage_on.pin --chimeric --chimeric-frag-index on 2>/tmp/sage_on.log; echo on=$?
grep SageIndex /tmp/sage_on.log; echo "off_rows=$(($(wc -l</tmp/sage_off.pin)-1)) on_rows=$(($(wc -l</tmp/sage_on.pin)-1))"
```
Expected: both exit 0; `SageIndex:` line in on-run; rows>0 both (on ≈ off if recall is high on this tiny fixture).

- [ ] **Step 5: Commit.**
```bash
git add crates/search/src/match_engine.rs
git commit -m "feat(chimeric): wire SageIndex query as chimeric candidate generator (B-P3)"
```

---

## Task 4: PXD recall gate + local per-spectrum timing (VM)

Empirical. Gate: index@1% ≥ 0.995×brute@1%; index wall < brute; entrapment FDP ~ brute. **The local microbenchmark (Task 2) already proved the query is sub-ms — Task 4 confirms it on real data.**

- [ ] **Step 1: Ship + rebuild on VM.**
```bash
for f in crates/search/src/sage_index.rs crates/search/src/match_engine.rs crates/search/src/lib.rs; do
  scp -o ControlPath=/tmp/msgfplus-bench.sock "$f" "pride-linux-vm:/srv/data/msgf-bench/chimeric-build/$f"; done
ssh -S /tmp/msgfplus-bench.sock pride-linux-vm 'cd /srv/data/msgf-bench/chimeric-build && source /root/.cargo/env && cargo build --release -p msgf-rust 2>&1 | tail -3'
```

- [ ] **Step 2: Reuse `fi_pxd_gate.sh`** (already on the VM from Approach A) — it runs PXD chimeric NO_RESCORE with `--chimeric-frag-index off` and `on`, Percolators each, and computes entrapment FDP. Run it:
```bash
ssh -S /tmp/msgfplus-bench.sock pride-linux-vm 'cd /srv/data/msgf-bench && nohup bash fi_pxd_gate.sh > fi-pxd-sage.out 2>&1 & echo pid=$!'
```
Poll for `T4_DONE`; read `targets_1pct` (off vs on), wall, and `entrapment_fraction`.
Gate: on@1% ≥ 0.995×off@1%; on wall < off wall (1:36); entrapment FDP < ~1%. If recall fails, raise TOP_K and re-run; if wall fails despite the sub-ms microbench, profile (likely the window derivation or build).

- [ ] **Step 3: Record** results in `docs/parity-analysis/notes/2026-05-30-sage-index-pxd-gate.md`, commit.

---

## Task 5: Astral speed/memory gate — beat Java wall (VM)

Empirical. Gate: Astral index wall **< Java 6:18**, index@1% ≥ 0.995×brute (77,287), MaxRSS within budget.

- [ ] **Step 1: Reuse `fi_astral_gate.sh`** (on the VM from Approach A; brute vs index, wall + MaxRSS + @1%, separate search/percolator logs). Run nohup'd; poll for `T5_DONE`.
- [ ] **Step 2:** If wall ≥ 6:18 or recall < 99.5%: profile/tune (TOP_K, the per-charge query union, the build). If index memory over budget: precursor-slab partition the fragments (a `sage_index.rs`-local change). 
- [ ] **Step 3:** Record `docs/parity-analysis/notes/2026-05-30-sage-index-astral-speed.md` (brute vs index vs Java table: wall, @1%, RSS); commit.

---

## Task 6: Cross-dataset + entrapment-preserved + PR (VM)

- [ ] **Step 1:** Run all 3 datasets with `--chimeric-frag-index on` (the `gate_chimeric_norescore.sh` shape) → final PSMs+wall table vs Java.
- [ ] **Step 2:** Re-run the entrapment harness with index on (PXD + Astral) → entrapment FDP ≈ brute (both ~nominal). This confirms the validated FDR is preserved through the index.
- [ ] **Step 3:** Write `docs/parity-analysis/notes/2026-05-30-sage-index-gate-result.md`; commit. Per `merge-gate-beat-java`: only merge/PR if chimeric+index now beats Java on PSMs AND speed on all 3 (TMT PSMs likely still trail → keep parked, PR as milestone). Push the branch.

---

## Self-Review notes (author)

- **Spec coverage:** §1 arch → T1+T3; §2 window-bounded buffer → T2 `query` (scores sized `pre_hi-pre_lo`); §3 build/memory → T1 + T3 eprintln; §4 gates/flag → reused scaffolding + T4/T5/T6; §5 tests incl. microbench → T2 Step 5; §6 phases → Tasks 1-6.
- **Type consistency:** `SageIndex::{build, query, n_fragments}`; `Frag{mz,pidx}`; `sorted_cand`/`sorted_mass`/`fragments`/`bucket_min_mz`; `PreparedSearch.sage_index`. Consistent across tasks.
- **Known integration subtlety (flagged in T3 Step 2):** deriving the neutral-mass `[lo,hi]` window consistent with `window_cand_indices` — the implementer must confirm the nominal↔neutral mapping or derive from precursor_mz directly; the T4 recall gate catches a mismatch.
- **No placeholders:** all code is concrete; the one judgement call (window derivation) is explicitly bounded with two viable methods + a catching gate.
