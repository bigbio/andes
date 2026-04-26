# Plan: search-path sync cleanup + per-task result buffers

**Status: SHIPPED in PR #25** (https://github.com/bigbio/msgfplus/pull/25)
Branch: `perf/search-sync-cleanup` (worktree at
`/Users/yperez/work/msgfplus-workspace/search-sync-cleanup`).

Successor to PR #24. Pure refactor + instrumentation — no scoring,
parser, or `.pin` feature changes. Output bit-identical to dev's tip
on every measurable axis.

## What shipped (6 commits)

1. **T1 — per-task wall stats + tail-imbalance summary**
   `RunMSGFPlus` captures preprocess / db-search / compute-evalue /
   total wall into a `TaskWallStats` accessor; `MSGFPlus.runMSGFPlus`
   prints a one-line summary at end of search:
   ```
   Task wall summary (n=12): min=101.7s median=224.2s p95=246.4s
     max=246.4s total=2356.7s tail_gap=22.2s (10% of median)
   ```
   On Astral the measured `tail_gap` is **10 % of median**, which means
   T2 and T3 can't deliver substantial wins on this workload.

2. **Drop dead `synchronized` wrappers in DBScanner + ScoredSpectraMap.**
   Each instance is task-local (verified: no internal fork-out in
   `dbSearch`, no shared instance across threads). Plain `HashMap` /
   `TreeMap` replace the `Collections.synchronizedMap` /
   `synchronizedSortedMap` wrappers; `synchronized` modifier dropped
   from `addDBMatches`, `generateSpecIndexDBMatchMap`,
   `addResultsToList`, `addDBSearchResults`. Memory-visibility safety
   preserved via `awaitTermination`'s happens-before.

3. **Per-task local result buffers + final merge.**
   Replaced the global `Collections.synchronizedList<MSGFPlusMatch>`
   with a per-task `ArrayList`. Each `RunMSGFPlus` owns its own buffer;
   main thread drains all buffers after `awaitTermination`.
   `RunMSGFPlus`'s constructor drops the `resultList` parameter; new
   `getResults()` accessor.

4. **T2 — `-Dmsgfplus.numTasksPerThread=N`** (default 3, unchanged).
   Lets operators raise the multiplier on datasets where T1's
   `tail_gap` shows real imbalance.

5. **T3 — `-Dmsgfplus.useForkJoin=true`** (default false, unchanged).
   Opt-in `ForkJoinPool` swap. Default keeps
   `ThreadPoolExecutorWithExceptions` (which retains progress
   reporting + exception-capture-via-afterExecute). FJP path uses
   `Future.get()` for exception propagation.

6. **Polish — tighter result-buffer merge + `drainResultsTo` + reused
   null sink.** Static `NULL_PRINT_STREAM` cached instead of allocated
   per `run()`; `drainResultsTo(dest)` clears per-task buffers
   immediately after merge so heap is collectible; pre-size merged
   `ArrayList` to `sum(t.getResultCount())` to avoid resize-and-copy;
   `submittedTasks.clear()` after summary drops strong refs to all 12
   task instances before the FDR / write phase.

## Validation gate cleared (Astral 3-arm + Percolator)

Astral 3-arm cold, 8 GB heap, 4 threads, default sysprops.
**All 8 parity numbers bit-identical to dev's tip:**

| Metric | dev | this branch |
|---|---:|---:|
| armB raw targets | 89,479 | 89,479 ✓ |
| armB raw decoys | 46,792 | 46,792 ✓ |
| armB 1 % FDR targets | 35,818 | 35,818 ✓ |
| armB 5 % FDR targets | 40,408 | 40,408 ✓ |
| armC raw targets | 89,360 | 89,360 ✓ |
| armC raw decoys | 46,913 | 46,913 ✓ |
| armC 1 % FDR targets | 35,767 | 35,767 ✓ |
| armC 5 % FDR targets | 40,426 | 40,426 ✓ |

Walltime delta vs master in the same run:
- armB: 752.2s vs 848.8s = **−11.4 %**
- armC: 798.2s vs 848.8s = **−5.9 %**

(First run came in with armC at 6298s; root-caused to OS thrashing —
load avg 5-8, ~120 MB free RAM, 165M page reclaims, Rancher VM eating
1 GB. Re-ran after stopping Rancher; wall normalized. Not a code
issue. Documented in PR #25 description.)

## What we learned vs. expected wins

The plan predicted:
- Step 1 (sync removal): 0–2 % wall. Possibly negative if biased
  locking was helping. Code clarity is the more reliable win.
- Step 2 (per-task buffers): 2–8 % wall, scaling with PSM count.
- T2 / T3: only worth doing if profiler shows real tail-imbalance.

What we measured:
- Combined wall improvement: **11.4 % on armB, 5.9 % on armC** —
  better than the upper end of the per-step predictions, suggesting
  the gains compound (less monitor traffic + cheaper drain phase).
- T1's measured tail_gap on Astral: **10 % of median** — small enough
  that T2/T3 default-on would give marginal wins. They ship as opt-in
  knobs precisely so they don't gate the default behavior.

## What this branch is NOT

Not a fragment-index revival. Not a primitive mass-window port. Not
a peak-storage refactor (`Peak` → `float[]`). Not a CLI / format
change. Originated from a third-party review of PR #24.

## Follow-ups (out of scope for this PR)

- **Profile on TMT and a metaproteomic FASTA** with the new T1
  summary. Astral's 10 % tail_gap might not represent uneven
  workloads — homolog-rich DBs are the place T2/T3 should bite.
- **`DatabaseMatch.indices` from `TreeSet<Integer>` to primitive
  `int[]`** (M1 from the broader memory-roadmap discussion). Highest
  expected impact for homolog-heavy databases (5-12× memory reduction
  per match); needs a metaproteomic test fixture to validate.
- **Parser cache stores raw `float[] mz, float[] intensity`** (M3),
  with a fresh `Spectrum` built per `getSpectrumBySpecIndex`. Side
  benefit: cache-layer immutability instead of cloneSpectrum.
- **`Peak`/`Spectrum` storage refactor** (M2). Multi-PR. Big surface
  area. Defer until M1 + M3 land.

## Open questions resolved

- **Did the custom `ThreadPoolExecutorWithExceptions` preserve
  awaitTermination's happens-before on the exception path?** Yes —
  observed bit-identical results in armB / armC across the 3-arm
  benchmark, which would not be the case if visibility were broken.

- **Was HotSpot already eliding the uncontended monitors?** Probably
  partially. Step 2 (sync removal) on its own gives an unmeasured
  delta; combined with steps 3–6 the total is 11.4 %. We can't
  attribute that 11.4 % to any single commit without per-commit
  benchmarks, but the polish commit (#6) likely contributes
  meaningfully via the pre-sized `ArrayList` and immediate
  per-task-buffer release.
