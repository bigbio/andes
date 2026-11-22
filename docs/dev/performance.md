# Performance: speed, memory, resources

Speed has priority, but never at the cost of a changed result: a speed-up must keep the PIN
bit-identical (see `bug-prevention`). andes beats Comet on wall time, CPU time and peak
memory on the reference sets (README tables); keep it that way.

## Measure first, every time

```bash
SPECTRUM=path/to/run.mzML DATABASE=path/to/db.fasta
```

```bash
/usr/bin/time -v target/release/andes --spectrum "$SPECTRUM" --database "$DATABASE" --output-pin out.pin --threads 8
```

Record wall, user+sys CPU and `Maximum resident set size`. Run at least twice; wall varies
5-15% between sessions, CPU time is steadier and is the better signal. Profile before
optimising (`perf record -g`, `cargo flamegraph`), and report the hot function's share.

State the before/after numbers, dataset and thread count in the commit message.

## Speed

- No allocation in per-spectrum or per-candidate inner loops. Reuse buffers
  (`clear()` + capacity), pass slices, avoid `clone()` of vectors and strings.
- Parallelise at the spectrum/chunk level with rayon; do not nest `par_iter` inside a
  per-spectrum closure.
- Keep data contiguous (struct-of-arrays, sorted `Vec` + binary search) over `HashMap` in hot
  paths.
- Do not enable FMA or change float evaluation order: `.cargo/config.toml` pins
  `target-cpu=sandybridge` so PINs stay bit-identical.
- Release builds use fat LTO and one codegen unit; benchmark only `--release` builds.

## Memory

- Stream inputs; never hold all spectra of a large run when a chunked pass works (see the
  out-of-core path and its goldens).
- Respect the memory budget in `crates/andes/src/bin/andes/memlimit.rs` (cgroup limit or
  `MemAvailable`). New large buffers must size themselves from that budget, not from the
  machine's total RAM.
- Prefer `f32` and `u32` where precision allows and results stay identical; drop large
  intermediate structures as soon as a stage ends.
- Every new structure that grows with the number of spectra, peptides or fragments needs
  its size estimated in the PR (bytes per item x expected count on the Astral set).

## Resources

- Honour `--threads`; never spawn threads beyond it. Respect cgroup CPU limits on clusters.
- Temporary files go to the output directory or the system temp dir and are removed on
  success and on error.
- No network access at run time.

## Regression gate

A change to a hot path reports wall, CPU and peak RSS on at least one reference set
(UPS1 is quickest) before and after. A regression above noise (about 5% CPU) needs a
reason in the PR.
