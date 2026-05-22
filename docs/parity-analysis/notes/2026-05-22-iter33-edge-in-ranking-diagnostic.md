# iter33 diagnostic: top-1 ranking uses node-only score; Java uses node+edge

_2026-05-22. Per-PSM trace on scan 21 (NEEQSR vs TEAPCGK decoy) identified the structural cause of the 40% diff-peptide bucket in the iter32 pin-diff. Closing the remaining 11.4% Astral PSM gap likely needs a top-1-ranking fix that adds edge_score to the queue ordering key — distinct from the iter17 attempt that modified the pin RawScore column directly._

## The trace

**Scan 21**: Java picks `R.NEEQSR.D` (target, RawScore=38); Rust picks `K.TEAPC+57.02146GK.P` (decoy, RawScore=32).

`MSGF_TRACE_PEP=NEEQSR` + `-Dmsgfplus.trace=true` on the single-scan MGF (scan 21 extracted from the Astral mzML):

| | per-split sum | cleavage | edge | total RawScore | source |
|---|---:|---:|---:|---:|---|
| Java NEEQSR | 14 | +4 | **+20** | **38** | DBScanScorer.getScore + DBScanner cleavage |
| Rust NEEQSR | 14 (bit-exact!) | +4 | +20 (computed but NOT in RawScore) | 18 in pin | psm_score node-only + iter19 separate EdgeScore PIN column |

**The score_psm path is BIT-EXACT for NEEQSR.** Both engines compute the per-split node sum = 14. The 24-point divergence is in **what counts toward top-1 ranking**:

- Java's `DBScanScorer.getScore` (override of FastScorer.getScore) returns **node + edge** (= 34). Then DBScanner adds cleavage (+4). So `match.score = 38` — the value that enters `PriorityQueue`'s heap-ordering AND becomes the pin RawScore.
- Rust's `score_psm` returns node only (= 14). `match_engine.rs` adds cleavage (+4). So `PsmMatch.score = 18` — what enters `TopNQueue` AND becomes the pin RawScore. `psm_edge_score = +20` is computed AFTER queue selection (iter19, additive PIN column).

## Why iter17 failed but iter33 should not

iter17 (REVERTED, -8K Astral PSMs) tried to add edge_score to `psm.score`. That changed BOTH the queue ordering AND the pin RawScore column. Percolator's learned weights against Rust's distribution broke.

The proposed iter33 design splits them:

```rust
struct PsmMatch {
    score: f32,        // node + cleavage — UNCHANGED; this is pin RawScore
    rank_score: f32,   // node + cleavage + edge — NEW; queue ordering key
    ...
}

impl Ord for PsmMatch {
    fn cmp(&self, other: &Self) -> Ordering {
        // primary: spec_e_value ASC (unchanged)
        // secondary: rank_score DESC (was score DESC)
    }
}
```

The PIN distribution stays exactly as iter19 emitted: `RawScore = node + cleavage`, `EdgeScore = +edge` (separate column). Percolator's learned weights remain valid. The CHANGE is which PSM ends up at top-1 per spectrum — Java-aligned now.

For scan 21: Rust would still emit `K.NEEQSR.D` with `RawScore = 18`, `EdgeScore = +20` (matching Java's effective +20). Percolator combines them.

## Implementation sketch

Three files touched, ~50 LOC total:

1. `crates/search/src/psm.rs`: add `rank_score: f32` field to `PsmMatch`; change `Ord::cmp` to use `rank_score` as secondary key.
2. `crates/search/src/match_engine.rs`: in the per-candidate scoring loop (~line 295-310), compute
   ```rust
   let node_cleav = score_psm(...) + cleavage_credit;
   let edge = psm_edge_score(scored_spec, &cand.peptide, scorer, z);
   let rank_score = node_cleav + edge as f32;
   ```
   and pass `(node_cleav, rank_score)` to PsmMatch.
3. `compute_psm_features`: reuse the stored `psm.edge_score` instead of recomputing.

## Perf concern

Adding `psm_edge_score` to the per-candidate loop roughly doubles per-PSM scoring work. 16M candidates × ~10 edges × peak-lookup cost. Astral wall could go 5:35 → ~9:00.

Mitigations (in order of complexity):
- **Two-stage gating**: only compute edge_score for candidates where `node + cleavage > queue.worst_score - MAX_EDGE_BONUS`. For Astral top_n=1 this evaluates edge_score on ~5-10 candidates per spectrum (the ones that COULD become top-1) instead of all 16M.
- **Precompute observed_node_mass cache** for the spectrum once; edge_score reads from cache.
- **Inline the edge loop** with node loop sharing prefix/suffix mass accumulators.

## Bench expectation

If the iter17 -8K loss was entirely from breaking Percolator's pin-distribution learning (and iter33 preserves it), and the underlying Java direction is correct (DBScanScorer.getScore includes edge for ranking), then iter33 should land in the same direction as iter17 INTENDED — closing toward Java's 35,818. Target: +500 to +2,000 Astral PSMs.

If iter33 ALSO regresses, the n=10 audit pattern needs a further refinement: "queue-ordering changes that flip top-1 selection are inherently destabilizing even when pin distributions stay intact". That would suggest the 11.4% gap is closeable only via per-feature tuning, not via algorithmic alignment.

## Test plan

1. Implement the 3-file change in a branch off iter32.
2. Build, run scoring + search tests.
3. Sanity-trace scan 21: confirm Rust's top-1 is now `R.NEEQSR.D` (not the decoy).
4. Sanity-trace scan 47106 (the iter28 reference case): confirm `R.HGIPTAQWK.A` still wins (it had +8 edge before, should still).
5. Full 3-dataset bench. Compare 1% FDR + wall.
6. If +PSMs and acceptable wall: ship. If perf regressed too much: add two-stage gating. If PSM regressed: bisect by sample, may need rollback.

## Why not in this session

The change requires PsmMatch struct change + queue ordering + careful re-test. Best done with a 60-90 min focused window, not as a tail action.

## Status

Iter32 perf cluster shipped with Rust now faster than Java on all 3 datasets (Astral 5:35 vs Java 5:49). PSM gap to Java is 11.4% (was 26% at iter16). iter33 is the proposed next algorithmic-fix iteration to close more of the PSM gap.
