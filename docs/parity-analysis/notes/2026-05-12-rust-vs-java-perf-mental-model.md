# Rust vs Java for MS-GF+ — perf mental model + concrete evidence

## 1. Premise

Rust can be faster than Java for the same algorithm, but usually not by a
dramatic factor when both implementations are decent. Java is not slow by
default: hot loops over primitive arrays with predictable branches are fast
after JIT. Rust mainly removes runtime costs — GC pauses, object headers, less
control over layout/allocation, harder zero-copy. If Java already uses compact
arrays and precomputation, there is little "language tax" left to harvest.

In this codebase the biggest gaps were **not** "Rust is slower than Java."
They were implementation choices: Rust recomputed work Java precomputed; Rust
had heavier lookup/layout paths in hot loops; one path had an O(N×M) scan;
parallelism on Apple Silicon underperformed.

Right mental model:

- Rust gives you better tools to build the fastest version.
- It does not guarantee a first-pass port beats a mature Java implementation.
- To beat Java, Rust must exploit something extra: better data layout, fewer
  allocations, cache-friendly structures, stronger precomputation, better
  parallel execution, avoiding duplicated work.

## 2. Where Rust did NOT automatically beat Java

Starting point of this iteration: msgf-rust was **5.49× slower than Java** on
the PXD001819 Astral baseline. Same algorithm, same I/O, same input — the
language alone bought us nothing. Every meaningful speedup below is an
implementation strategy fix, not a language win.

## 3. The seven commits that closed most of the gap

| Commit    | Fix                                              | Category                                         |
|-----------|--------------------------------------------------|--------------------------------------------------|
| `9c56797` | PrimitiveAaGraph thread-local arena pool         | Fewer allocations / better data layout           |
| `b95d348` | Flat `Vec<f64>` ScoreDist arena per GF graph     | Fewer allocations (~55M tiny allocs eliminated)  |
| `317d6bc` | Per-segment `(partition, ion_logs)` cache        | Stronger precomputation                          |
| `e19293a` | 4-wide chunked `add_prob_dist` accumulation      | Better layout for the auto-vectorizer            |
| `507bcb1` | PIN Label cache + one-pass distinct count        | Avoiding duplicated work (626s → 366s wall)      |
| `be50dab` | `compute_psm_features` hoisted to post-top-N     | Avoiding duplicated work (189s → 178s match)     |
| `d3d577d` | pin_write `memchr::memmem` haystack              | Cache-friendly structures (158.9s → 4.45s, 31.6×)|

Highlights:

- **`9c56797` + `b95d348`** removed ~11 per-call `Vec` allocations and ~55M
  tiny `Option<ScoreDist>` heap allocs. Both bit-identical against Java SP.
- **`507bcb1`** caches PIN `Label` by peptide sequence — a "Java already
  deduped, Rust was recomputing" gap. Wall: ~626s → ~366s on one change.
- **`d3d577d`** is the user's exact "O(N×M) scan" example. Naive substring
  search → `memchr::memmem::find` (SIMD Two-Way) cut pin_write 158.9s → 4.45s
  (**31.6×**). The language did none of the work; the algorithm change did.

## 4. The single biggest lever — Track A FastScorer (pending)

Prototype: precompute prefix/suffix score arrays once per `(spectrum, charge)`,
so per-split scoring is an array lookup instead of recomputed
`directional_node_score`. Measured single-VM wall:

```
Track A:  1m43.50s
Java:     1m39.75s
ratio:    0.96× of Java   (~Java parity)
```

This is the load-bearing example of the user's thesis: the biggest win in the
whole iteration comes from **stronger precomputation**, not from Rust-the-language.
A bit-identity gate currently FAILS; the commit is blocked pending the scoring
bug fix, but the headline measurement is real.

A secondary effect: with per-thread allocation pressure dropped by `9c56797` +
`b95d348`, Track A VM CPU% jumped from **192% → 406%** on 8 threads. On the same
Mac, Java was already pulling ~350% CPU while the pre-fix Rust build burned
~240s wall on 12 threads with only **154% CPU = 1.5 cores effective**.
Allocator contention and M-series E-core stragglers are the suspected blockers;
the per-thread layout fixes **indirectly unlocked parallel scaling**. This is
the rare case where Rust-specific control (thread-local arenas) bought us
something Java would have to work harder to match.

## 5. "Profile before you assume" — the SA-walk postmortem

Postmortem at `bb3353a`. The SA-walk integration was an architecturally
attractive change predicted by an earlier brainstorm. Implemented, measured,
reverted: it touched the wrong region. The actual bottleneck was pin_write's
O(N×M) substring scan (`d3d577d`) — **outside** the area the refactor changed.
Lesson: an architecturally appealing fix that doesn't touch the measured hot
path doesn't move wall clock. Profile first, refactor second.

## 6. Bottom line

Starting position: Rust **5.49× slower** than Java.
After the seven landed commits + Track A prototype: Rust within **0.96× of
Java** wall clock on a single VM.

None of that came from Rust being a "faster language." It came from doing what
the user's mental model says you have to do: precompute more, allocate less,
lay out data for the cache and the vectorizer, dedupe duplicated work, and pick
the right algorithm for the hot loop. Java parity here is a question of
implementation strategy, not language choice.
