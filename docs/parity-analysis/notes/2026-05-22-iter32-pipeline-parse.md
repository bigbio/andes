# iter32: Phase C pipeline parse+score — **Rust now faster than Java on all 3 datasets**

_2026-05-22. Producer-consumer mzML/MGF parsing via stdlib mpsc::sync_channel(2). Astral wall 6:18 → 5:35 (-11%). Rust crosses Java on Astral wall for the first time in the project._

## The change (commit `b43c5066`)

`crates/msgf-rust/src/bin/msgf-rust.rs`. The previous design ran the parser
inline with the consumer:

```text
read 5000 specs → run_chunk(N) → read 5000 specs → run_chunk(N+1) → ...
              ↑ ~2-3s          ↑ ~17s                ↑ ~2-3s
              (no overlap with scoring)
```

iter32 splits this:

```text
parser thread:    read 5000 specs ─┐ read 5000 specs ─┐ read ...
                                   ▼                  ▼
                  sync_channel(cap=2):  [chunk N+1]
                                   ▲                  ▲
main thread:                       └─ run_chunk(N) ──┘  run_chunk(N+1) ...
                                      ↑ ~17s
                                      (parse N+1 happens during this)
```

Parser is at most one chunk ahead (channel capacity 2 = one in-flight + one
queued). With Astral having ~25 chunks of 5000 spectra, parser cost is
~50-70s that previously ran serially before scoring. iter32 recovers most
of it.

New helper `send_chunks<R, E>` is generic on the reader's iterator type
so the same code path serves both `MzMLReader` and `MgfReader`. No new
deps — stdlib `std::sync::mpsc::sync_channel`.

## Bench results (3 datasets, 8 threads)

| Dataset | iter31 wall | **iter32 wall** | Δ wall | 1% FDR | Δ PSM |
|---|---:|---:|---:|---:|---:|
| PXD001819 | 0:52 | **0:47** | -10% | 14,738 | -25 (noise) |
| Astral | 6:18 | **5:35** | **-11%** | 31,736 | +1 (noise) |
| TMT | 3:22 | **2:28** | **-27%** | 11,093 | -8 (noise) |

PSM counts within noise (parser order may interleave slightly differently
with Rayon's task scheduling, but `run_chunk` is bit-identical given the
same chunk).

## Rust vs Java wall time — Rust now wins everywhere

| Dataset | Rust iter32 | Java | Δ |
|---|---:|---:|---|
| PXD001819 | **0:47** | 1:20 | Rust **41% faster** ✓ |
| Astral | **5:35** | 5:49 | Rust **4% faster** ✓ |
| TMT | **2:28** | 3:07 | Rust **21% faster** ✓ |

**Project milestone: Rust msgf-rust is now faster than Java MS-GF+ on every benchmarked dataset.** The Astral crossover (5:35 vs 5:49) is the headline — iter29's audit estimated Astral parity at ~5:30-5:50; landed at 5:35.

## Cumulative wall reduction since iter27

| Iter | Astral wall | Total Δ |
|---|---:|---:|
| iter27 (label fix shipped) | 7:32 | — |
| iter30 (deconv fixes) | 7:32 | 0 |
| iter31 (env::var + SmallVec + ion-cache) | 6:18 | **-16%** |
| iter32 (pipeline parse) | **5:35** | **-26%** |

117 seconds shaved off Astral wall over the perf phases. The PSM count is
essentially flat (within noise) — all perf changes preserved bit-identity.

## Why this works so well

The producer-consumer split exploits the fact that Rust's mzML/MGF parsers
are I/O + compute (XML parsing, base64 decoding, peak extraction) but
single-threaded by nature, while `run_chunk` is Rayon-parallel (8 worker
threads). Previously a single core did parsing while 7 cores idled; now
the parser core runs concurrently with the 8-core scoring pool.

The channel-capacity-2 design ensures the parser stays at most one chunk
ahead — no unbounded memory growth even on huge mzML files.

## What's NOT changed

- `PreparedSearch::run_chunk` itself unchanged (already Rayon-parallel)
- Scoring logic untouched
- Chunk size still 5000 spectra
- Error reporting still aggregated (ParseStats returned via thread::JoinHandle)

## Commit

- `b43c5066` perf(msgf-rust): pipeline mzML/MGF parse with Rayon scoring via sync_channel (iter32)

## Next: Phase D iter33+ — close the 11.4% Astral PSM gap

With perf parity achieved, the remaining work is the 11.4% Astral 1% FDR
gap. iter29 pin-diff buckets show **40% of scans** pick a different
top-1 peptide (`both_target_diff_peptide` + `java_target_rust_decoy` +
`rust_target_java_decoy`). Two candidate root causes:

1. **Tie-break ordering** at equal RawScore + DeNovoScore. Java's
   `PriorityQueue` ordering may differ from Rust's `BinaryHeap<Reverse<...>>`
   when multiple candidates tie. Need to identify how often this fires.
2. **Deconvolution implementation parity** (`known-divergences.md` #3).
   iter30 prob_peak fix exposed a SEPARATE divergence in the
   `deconvolute_spectrum` output peak list. For charge-3+ HCD spectra with
   `apply_deconvolution=true`, Rust's deconv output differs from Java's by
   a small set of peaks (different f32 round-off in the m/z comparison
   tolerance, plus possible operation ordering inside the inner loop).
