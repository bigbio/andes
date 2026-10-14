# Scalable LFQ match-between-runs (P0 + P1) Implementation Plan

> Part A was revised during implementation; where it differs, the code and DOCS.md are
> authoritative.

**Goal:** Give every (precursor, run) pair in a multi-run LFQ experiment one consistent,
FDR-controlled quantity, using per-run stages that run in parallel and global stages that
only touch compact tables.

**Architecture:** `andes search --quant-dir D` writes a persistent run cache and a
per-run PSM table. `andes quant library` builds the precursor library: global confidence,
joint RT alignment and decoys. `andes quant extract` measures every library precursor in
one run by fitting a template to that precursor's expected signal. `andes quant finalize`
trains the transfer FDR and writes the tables. `andes search --lfq --mbr` with several
inputs runs the same four stages in one process.

**Tech Stack:** Rust 1.87 workspace, crates `quant` (pure algorithms), `output`
(writers), `andes` (driver/CLI), `model-train` (GBDT); arrow/parquet 53; zstd.

**Spec:** `docs/plans/2026-10-10-scalable-quant-design.md` (sections 5-8, phases P0 and P1).

## Global Constraints

- Build, test, clippy and benchmark on a dedicated build machine with its own `CARGO_TARGET_DIR`; the local machine does edits and git only.
- CI gate: `RUSTFLAGS="-D warnings" cargo build --release`, `cargo test --release --workspace`, `cargo clippy --workspace --all-targets -- -D warnings`.
- No hostnames, usernames, cluster names, job ids or absolute site paths in code, commits, docs or PR text.
- Code comments describe the code as it is; tried and rejected ideas go in the task report, not in comments.
- andes reports PSM-level identification FDR only through Percolator (or the labelled native rescorer). The library entry gate and the transfer q-value are internal quant-quality scores and must be labelled as such in help text and output column names.
- Existing outputs (`.tmt.tsv`, `.lfq.tsv`, `.lfq_features.tsv`, PIN, idparquet) stay byte-identical when `--quant-dir`/`--mbr` are not given.
- One commit per task on branch `feat/quant-mbr`; no push until the user asks.
- Every benchmark A/B changes one variable; report peptides/precursors along with PSMs.

## Review Focus

- FAIMS runs: a precursor identified at one compensation voltage must be extracted from the MS1 series at that voltage in the target run, never from a series that mixes voltages (Task 8 test `extract_uses_donor_faims_group`).
- A run with very few identifications, or an RT range that barely overlaps the others, must be marked unaligned and skipped as a transfer source and destination, with a warning and no panic (Task 5 test `run_with_few_anchors_is_unaligned`).
- Precursors identified in only one run still need a consensus RT and still get transferred, but they must not steer the alignment (Task 5 test `single_run_precursor_gets_consensus_without_anchoring`).
- A missing, truncated or wrong-version run cache or quant table must fail with an error that names the file, never with a panic or silent zeros (Task 1 test `cache_rejects_wrong_magic_and_version`, Task 2 test `quant_table_rejects_unknown_major_version`).
- Zero transfers, zero decoys, or a fold where all rows share one class must yield `q = None` for those rows and a logged count, not a panic and not q = 0 (Task 10 test `fdr_degenerate_inputs_yield_none`).

---

## Part A: algorithm design (read before Tasks 5-10)

This part defines the method. The tasks implement it, and their tests check it.

### A1. The question per (precursor, run)

For library precursor `p` (charge `z`, monoisotopic m/z `mu`, theoretical isotope fractions
`pi_0..pi_{I-1}` with I = 4, summing to 1) and run `r`, the alignment gives a predicted apex
`t_hat = f_r(t_p)` and an RT uncertainty `s` (A3). The question is whether the MS1 data near
(`mu`, `t_hat`) contains `p`'s co-eluting isotope pattern, and if so how much of it.

We do not pick "the best peak in a window" and then score it with hand-weighted features.
Instead, the hypothesis "p elutes here with this shape and these isotope ratios" is fitted
directly to the data, and one likelihood-ratio statistic summarises the evidence.
Picking-then-scoring has two problems. It mixes interference into the pick. And its decoys
rarely replay the same max-over-window search, so the FDR is miscalibrated.

### A2. The template fit (`quant::evidence`)

Data: for isotope `i`, the trace `y_i[j]` is the most intense centroid within `tol_ppm` of
`mu + i * 1.003355 / z`, at MS1 scans `j` whose RT `t_j` lies in
`[t_hat - W, t_hat + W]`, where `W = 4 s + 2 w_prior`. There is also one extra trace, `y_m1`, at
`mu - 1.003355 / z`. It should be empty when `mu` really is a monoisotopic peak.

Model for isotope `i`: `y_i[j] = b_i + A * pi_i * g(t_j; t0, w) + noise`, where
`g = exp(-(t - t0)^2 / (2 w^2))`.

- Background `b_i` is the 20th percentile of `y_i` over the window.
- Noise `sigma_i` is `1.4826 * median|y_i[j+1] - y_i[j]| / sqrt(2)`, floored at
  `1e-3 * max(y_i)` and at 1.0. Using successive differences makes it robust to the peak itself.
- Fit over `t0` on the scan RTs inside `[t_hat - 3 s, t_hat + 3 s]`, and over
  `w in {0.7, 1.0, 1.4} * w_prior`. For each candidate, `A` is the closed-form weighted
  least-squares amplitude over the used isotopes:
  `num = sum_i pi_i * sum_j (y_i[j] - b_i) g_j / sigma_i^2`,
  `den = sum_i pi_i^2 * sum_j g_j^2 / sigma_i^2`, `A = max(0, num / den)`.
  The log-likelihood ratio against "no peak" is `LLR = A * num - A^2 * den / 2`.
- Objective: `LLR - (t0 - t_hat)^2 / (2 s^2)`. The RT prior is part of the fit, not a
  separate feature applied after the pick.
- Interference exclusion, at most one isotope out of `{1, 2, 3}`. For each used isotope
  compute its own amplitude `A_i = sum_j (y_i - b_i) g_j / (pi_i sum_j g_j^2)` and its
  standard error `se_i = sigma_i / (pi_i sqrt(sum_j g_j^2))`. If
  `(A_i - A) / sqrt(se_i^2 + A_se^2) > 3` for the worst isotope, drop it and refit at the
  same `(t0, w)`. Interference only ever adds signal, so the test is one-sided.

Outputs (`Evidence`):

| field | meaning |
|---|---|
| `amplitude`, `amplitude_se` | `A`, `1/sqrt(den)` |
| `quantity`, `quantity_se` | `A * w * sqrt(2 pi)` and its SE (area of the used-isotope-normalised pattern) |
| `llr` | log-likelihood ratio at the optimum |
| `rt_z` | `(t0 - t_hat) / s` |
| `width_ratio` | `w / w_prior` |
| `envelope_cos` | cosine between `(A_i)` and `pi` over the used isotopes |
| `coherence` | minimum Pearson correlation between `y_i - b_i` and `g` over `[t0 - 2w, t0 + 2w]` for isotopes 0 and 1 |
| `m1_ratio` | M-1 amplitude projected on `g`, divided by `A * pi_0`; near 0 is good |
| `mass_err_ppm` | intensity-weighted ppm error of the isotope-0 centroids over `[t0 - w, t0 + w]` |
| `n_points` | scans inside `[t0 - 2w, t0 + 2w]` |
| `isotopes_used` | bit mask |
| `t0`, `w` | fitted apex and width |

The same estimator measures identified pairs. The only differences are that the prior is
centred on the PSM's MS2 RT and that `s = max(w_prior, 6 s)` (apex vs trigger offset). So
identified and transferred values are the same kind of number.

### A3. Joint RT alignment (`quant::align`)

Model: `rt[p, r] = f_r(t_p) + e`. `f_r` is monotone piecewise-linear and the errors `e`
are heavy-tailed.

- **Anchors:** confident precursors (library q <= 0.01) observed in at least 2 runs. Each
  `(p, r)` observation uses the RT of the best-scoring PSM of `p` in `r`.
- **Initialisation:** `t_p` is the median of `rt[p, .]`.
- **Iteration** (fixed 8 rounds, deterministic):
  1. For each run, form pairs `(t_p, rt[p, r])` and bin them by `t_p` into
     `K = clamp(n / 25, 4, 60)` equal-count bins. Each bin gives the median x and median y.
     Make y non-decreasing with pool-adjacent-violators, then interpolate linearly between
     knots. Outside the knots, extrapolate with the slope of the outermost 20% of anchors.
  2. For each precursor, `t_p` is the median over runs of `f_r^{-1}(rt[p, r])`.
  3. Fix the gauge: the mean map `g(t) = mean_r f_r(t)` is evaluated on the knot grid.
     Set `t_p <- g(t_p)` so the consensus scale is the average run, not an arbitrary one.
- **Residual width `sigma_r(t)`:** `1.4826 * MAD` of `rt - f_r(t_p)` in the same bins,
  smoothed with a running median of 3 and floored at 0.5 times the median MS1 scan spacing.
- **Per-precursor extra dispersion:** `tau_p^2 = max(0, mean_r e^2 - mean_r sigma_r^2) * n_p / (n_p + 3)`.
  The factor shrinks it towards 0 for precursors seen in few runs.
- **Transfer uncertainty into run r:**
  `s = sqrt(sigma_r(t)^2 * (1 + 1 / n_p) + tau_p^2)`.
- **Precursors seen in one run** do not take part in the fit. Their consensus is
  `t_p = f_r^{-1}(rt)` from the final maps.
- **A run is unaligned** when it has fewer than 50 anchors, or when its anchors cover less
  than 30% of the consensus RT range. It is excluded from the fit, is never a transfer
  source or destination, and is logged.

Why no reference run: any reference carries its own noise and gaps into every map. The
joint fit uses all runs symmetrically, and its cost is O(observations x rounds).

### A4. Decoys (`quant::decoys`)

A decoy must be a hypothesis that is certainly false but statistically interchangeable with
a false target:

- the same density of real signal nearby
- the same search breadth (`s` and `W` of the target it shadows)
- the same optimisation (the A2 fit, max over `t0` and `w`, isotope exclusion)

Two families, each one per target per run, so `K = 2`:

1. **Mass-shift decoy.** `mu' = mu + n * 1.000495 / z`, with `n` drawn from
   `{-11, -7, +7, +11}` by a hash of (precursor, run). The factor 1.000495 Da per nominal
   Da keeps the decoy in the peptide mass-defect band. A fractional shift would land where
   peptides rarely occur, so decoys would see less background than targets and the FDR
   would be underestimated. The values 7 and 11 Da sit away from common modification and
   adduct deltas (0.984, 14.016, 15.995, 17.027, 18.011, 21.982, 27.995, 31.990, 37.956) by
   at least 0.9 Da. The isotope fractions `pi`, prior `t_hat`, `s` and `w_prior` stay the
   target's.
2. **RT-shift decoy.** The same `mu` and `pi`, but `t_hat' = t_hat +/- D` with
   `D = max(8 s + 4 w_prior, 0.08 * run RT span)`. The sign comes from the hash and is
   flipped if `t_hat'` would leave the run's RT range; if neither sign fits, no RT decoy is
   made. `s` and `w_prior` are kept from the target.

The two families fail differently. The mass-shift decoy keeps the co-elution neighbourhood
but changes the m/z neighbourhood. The RT-shift decoy does the opposite. Finalize reports
the FDR estimate for each family separately, so a disagreement between them is visible. The
entrapment benchmark (Task 13) decides whether to keep both.

### A5. Transfer FDR (`andes::transfer_fdr`)

- **Rows:** transfer targets (precursor in the library, not identified in the run) and their
  decoys. Identified pairs are not part of the FDR; Percolator already established their
  identity.
- **Features:** `llr`, `rt_z`, `|rt_z|`, `width_ratio`, `envelope_cos`, `coherence`,
  `m1_ratio`, `|mass_err_ppm|`, `n_points`, `popcount(isotopes_used)`,
  `log10(amplitude / amplitude_se + 1)`, `n_donor_runs`, and the library q-value.
- **Never a feature:** any measure of agreement between this run's quantity and other runs'
  quantities. A true 10-fold change would look "inconsistent" and be rejected, which
  compresses the fold changes we are trying to measure.
- **Model:** the existing `model_train::gbdt::train::train_gbdt`, single pass with 3-fold
  CV. The fold is a hash of the precursor id, so a precursor's targets and decoys in every
  run fall in the same fold and nothing leaks across runs.
- **q-value:** rank by held-out score `c`.
  `FDR(c) = (D(c) / K) / max(T(c), 1)`, where `D` and `T` count decoys and targets at or
  above `c`. q is the running minimum from the bottom. This is the K-decoy estimator; the
  competition variant is not used because most transfer targets are true nulls, and
  competition would delete them from the count.
- **Per-family check:** the same estimate is computed with one decoy family at a time
  (`K = 1`), and both are logged next to the pooled one.
- **Pseudo-transfer validation:** every identified pair is also measured as if it were a
  transfer, with the alignment prior instead of its PSM RT, and scored by the same model.
  Its true answer is known: present, with its apex at the identified fit's `t0`. Finalize
  reports:
  - the fraction at q <= 0.01, a lower bound on sensitivity
  - the fraction with `|t0_transfer - t0_identified| > w`, the wrong-peak rate

  The wrong-peak rate is the one error decoys cannot see.

### A6. Library confidence (`quant::library::precursor_qvalues`)

The precursor key is the peptidoform string plus the charge. Take each precursor's best
PSM score over all runs, for targets and decoys alike. Sort descending and set
`q = (D + 1) / max(T, 1)`, then the running minimum. This is best-PSM-per-precursor
target/decoy competition over Percolator's (or the native rescorer's) PSM scores. It is an
internal entry gate (q <= `--mbr-library-fdr`, default 0.01), not a reported identification
FDR.

### A7. Run width prior `w_prior`

Per run, fit all identified pairs with `w` free on the grid `{0.5, 0.7, 1, 1.4, 2} * w0`,
where `w0` is the median of PR #119's per-run half-width (or 6 s when unavailable). Then set
`w_prior(t)` to the running median of the fitted `w` over 20 equal-count RT bins. Transfers
use `w_prior` at `t_hat`.

---

## Part B: tasks

### Task 1: Persistent run cache (P0)

**Files:**
- Modify: `crates/andes/src/bin/andes/ms1_cache.rs`
- Modify: `crates/andes/Cargo.toml` (add `zstd = "0.13"`)

**Interfaces:**
- Consumes: `quant::ms1_index::Ms1RunIndex`, `model::scan::Ms1Scan`.
- Produces: `Ms1Cache::store_at(index: &Ms1RunIndex, path: &Path) -> io::Result<()>`,
  `Ms1Cache::load_from(path: &Path) -> io::Result<Ms1RunIndex>`. The existing
  temp-file `store`/`load` keep their behaviour.

The file layout is the magic `b"ANDESRC1"`, a little-endian `u16` major version = 1, a
`u16` minor = 0, and then a zstd (level 3) stream of the existing body encoding
(`write_body` / `read_body` factored out of `store`/`load`).

- [ ] **Step 1: Write the failing tests** (in `ms1_cache.rs` `mod tests`)

```rust
#[test]
fn persistent_round_trip_is_lossless() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("run.cache");
    let index = Ms1RunIndex::new(vec![
        Ms1Scan { rt: 3.123456789, peaks: vec![(600.123456789, 1234.5678)], id: "scan=2".into(), faims_cv: Some(-45.0) },
        Ms1Scan { rt: 1.0, peaks: vec![], ..Ms1Scan::default() },
    ]);
    Ms1Cache::store_at(&index, &path).unwrap();
    let loaded = Ms1Cache::load_from(&path).unwrap();
    assert_eq!(loaded.len(), 2);
    for i in 0..2 {
        assert_eq!(loaded.rt(i).to_bits(), index.rt(i).to_bits());
        assert_eq!(loaded.peaks(i), index.peaks(i));
        assert_eq!(loaded.scan(i).faims_cv, index.scan(i).faims_cv);
    }
}

#[test]
fn cache_rejects_wrong_magic_and_version() {
    let dir = tempfile::tempdir().unwrap();
    let bad = dir.path().join("bad.cache");
    std::fs::write(&bad, b"NOTACACHE").unwrap();
    let e = Ms1Cache::load_from(&bad).unwrap_err().to_string();
    assert!(e.contains("bad.cache") && e.contains("not an andes run cache"), "{e}");
    let mut v2 = b"ANDESRC1".to_vec();
    v2.extend_from_slice(&2u16.to_le_bytes());
    v2.extend_from_slice(&0u16.to_le_bytes());
    std::fs::write(&bad, &v2).unwrap();
    let e = Ms1Cache::load_from(&bad).unwrap_err().to_string();
    assert!(e.contains("version 2"), "{e}");
}
```

- [ ] **Step 2: Run them on the build machine and see them fail**

Run: `cargo test --release -p andes --bin andes ms1_cache`
Expected: compile error, `store_at`/`load_from` not found.

- [ ] **Step 3: Implement**

Factor the body loop of `store` into `fn write_body(index: &Ms1RunIndex, out: &mut impl Write) -> io::Result<()>`,
and the body of `load` into `fn read_body(input: &mut impl Read, remaining: Option<u64>) -> io::Result<Ms1RunIndex>`.
When `remaining` is `None` (zstd stream), bounded counts use `u32::MAX as u64` as the cap
and trailing-data detection reads one byte past the end. Then:

```rust
const MAGIC: &[u8; 8] = b"ANDESRC1";
const MAJOR: u16 = 1;

pub(super) fn store_at(index: &Ms1RunIndex, path: &Path) -> io::Result<()> {
    let tmp = path.with_extension("cache.partial");
    {
        let mut f = BufWriter::new(std::fs::File::create(&tmp)?);
        f.write_all(MAGIC)?;
        f.write_all(&MAJOR.to_le_bytes())?;
        f.write_all(&0u16.to_le_bytes())?;
        let mut z = zstd::Encoder::new(f, 3)?;
        write_body(index, &mut z)?;
        z.finish()?.flush()?;
    }
    std::fs::rename(&tmp, path)
}

pub(super) fn load_from(path: &Path) -> io::Result<Ms1RunIndex> {
    let named = |e: io::Error| io::Error::new(e.kind(), format!("{}: {e}", path.display()));
    let mut f = BufReader::new(std::fs::File::open(path).map_err(named)?);
    let mut head = [0u8; 12];
    if f.read_exact(&mut head).is_err() || &head[..8] != MAGIC {
        return Err(io::Error::new(io::ErrorKind::InvalidData,
            format!("{}: not an andes run cache", path.display())));
    }
    let major = u16::from_le_bytes([head[8], head[9]]);
    if major != MAJOR {
        return Err(io::Error::new(io::ErrorKind::InvalidData,
            format!("{}: run cache version {major}, this andes reads version {MAJOR}", path.display())));
    }
    let mut z = zstd::Decoder::new(f).map_err(named)?;
    read_body(&mut z, None).map_err(named)
}
```

The temporary `store`/`load` keep writing the uncompressed body, unchanged.

- [ ] **Step 4: Run the tests and see them pass**

Run: `cargo test --release -p andes --bin andes ms1_cache`. Expected: 3 passed.

- [ ] **Step 5: Commit**

```bash
git add crates/andes/Cargo.toml crates/andes/src/bin/andes/ms1_cache.rs Cargo.lock
git commit -m "quant: persistent, versioned, zstd-compressed run cache"
```

### Task 2: Per-run PSM table `<stem>.quant.parquet` (P0)

**Files:**
- Create: `crates/output/src/run_quant_table.rs`
- Modify: `crates/output/src/lib.rs` (add `pub mod run_quant_table;`)
- Test: in-module tests

**Interfaces:**
- Produces:

```rust
pub const RUN_QUANT_SCHEMA_MAJOR: u32 = 1;

#[derive(Debug, Clone, PartialEq)]
pub struct RunPsmRow {
    pub spec_id: String,        // PIN SpecId, joins Percolator PSMId
    pub rt_s: f64,              // MS2 retention time, seconds
    pub precursor_key: String,  // "<peptidoform>/<charge>"
    pub peptidoform: String,
    pub charge: u8,
    pub mz: f64,                // theoretical monoisotopic m/z
    pub is_decoy: bool,
    pub search_score: f32,      // andes RawScore (rank-1)
    pub faims_cv: Option<f32>,
    pub envelope: Vec<f32>,     // theoretical isotope fractions, 4 values, sum 1
    pub proteins: String,       // ';'-joined accessions
}

pub fn write_run_quant_table(path: &Path, run_stem: &str, rows: &[RunPsmRow]) -> Result<(), String>;
pub fn read_run_quant_table(path: &Path) -> Result<(String, Vec<RunPsmRow>), String>; // (run_stem, rows)
```

Parquet key-value metadata: `andes.run_quant.schema = "1.0"`, `andes.run_quant.run = <stem>`.
The reader rejects a different major version with an error naming the file.

- [ ] **Step 1: Write the failing tests**

```rust
#[test]
fn round_trip_preserves_rows_and_run() {
    let dir = tempfile::tempdir().unwrap();
    let p = dir.path().join("a.quant.parquet");
    let rows = vec![RunPsmRow {
        spec_id: "a/scan=7_2_1".into(), rt_s: 1234.5, precursor_key: "PEPTIDEK/2".into(),
        peptidoform: "PEPTIDEK".into(), charge: 2, mz: 465.7380, is_decoy: false,
        search_score: 41.0, faims_cv: Some(-45.0), envelope: vec![0.6, 0.3, 0.08, 0.02],
        proteins: "P1;P2".into(),
    }];
    write_run_quant_table(&p, "a", &rows).unwrap();
    let (stem, back) = read_run_quant_table(&p).unwrap();
    assert_eq!(stem, "a");
    assert_eq!(back, rows);
}

#[test]
fn quant_table_rejects_unknown_major_version() {
    let dir = tempfile::tempdir().unwrap();
    let p = dir.path().join("b.quant.parquet");
    write_with_schema_version(&p, "b", &[], "2.0").unwrap(); // test-only helper
    let e = read_run_quant_table(&p).unwrap_err();
    assert!(e.contains("b.quant.parquet") && e.contains("schema 2.0"), "{e}");
}
```

- [ ] **Step 2: Run and see them fail.** `cargo test --release -p output run_quant_table`

- [ ] **Step 3: Implement** with arrow `RecordBatch` columns: `Utf8, Float64, Utf8, Utf8, UInt8, Float64, Boolean, Float32, Float32 (nullable), List<Float32>, Utf8`. Follow the writer pattern in `crates/output/src/feature_parquet.rs` (ArrowWriter, ZSTD compression, `WriterProperties::builder().set_key_value_metadata(...)`). `write_with_schema_version` is `#[cfg(test)] pub(crate)`, and `write_run_quant_table` delegates to it with `"1.0"`.

- [ ] **Step 4: Run and see them pass.**

- [ ] **Step 5: Commit** `git commit -m "output: versioned per-run PSM table for multi-run quantification"`

### Task 3: `--quant-dir` writes the cache and the table during search (P0)

**Files:**
- Modify: `crates/andes/src/bin/andes/cli.rs` (flag `--quant-dir <DIR>`)
- Modify: `crates/andes/src/config.rs` (`quant.dir`)
- Modify: `crates/andes/src/bin/andes/quant_run.rs` (`QuantSettings.quant_dir: Option<PathBuf>`; `fn write_run_artifacts`)
- Modify: `DOCS.md`, `config.example.yaml`
- Test: `crates/andes/tests/quant_end_to_end.rs`

**Interfaces:**
- Consumes: Task 1 `Ms1Cache::store_at`, Task 2 `write_run_quant_table`/`RunPsmRow`,
  the existing `select_hits`, `quant_id`, `faims_groups` and `formula::Formula::envelope`.
- Produces: for each input file, `<DIR>/<stem>.cache` and `<DIR>/<stem>.quant.parquet`.
  Rows are every rank-1 PSM, targets and decoys, with no FDR filter. `--quant-dir` implies
  MS1 capture (`needs_ms1() == true`) and has the same input restrictions as `--lfq`
  (mzML or `.raw`; not with `--chimeric`, MGF, `.d` or `--glyco`).

- [ ] **Step 1: Write the failing end-to-end test** (uses the existing `lfq_synthetic.mzML.gz` fixture and helpers in `quant_end_to_end.rs`)

```rust
#[test]
fn quant_dir_writes_cache_and_table_with_decoys() {
    let out = tempfile::tempdir().unwrap();
    let qdir = out.path().join("q");
    run_andes(&[ "--lfq", "--quant-dir", qdir.to_str().unwrap() /* + fixture args as in lfq test */ ]);
    let stem = "lfq_synthetic";
    assert!(qdir.join(format!("{stem}.cache")).exists());
    let (run, rows) = output::run_quant_table::read_run_quant_table(&qdir.join(format!("{stem}.quant.parquet"))).unwrap();
    assert_eq!(run, stem);
    assert!(rows.iter().any(|r| r.is_decoy), "decoy PSMs must be written");
    assert!(rows.iter().filter(|r| !r.is_decoy).count() >= 6, "the six designed BSA peptides");
    for r in &rows { assert!((r.envelope.iter().sum::<f32>() - 1.0).abs() < 1e-3); }
}

#[test]
fn outputs_unchanged_without_quant_dir() {
    // Run the existing LFQ fixture twice, with and without --quant-dir; the .lfq.tsv files are byte-identical.
}
```

- [ ] **Step 2: Run and see them fail.**
- [ ] **Step 3: Implement.** In `run_quant`, when `s.quant_dir` is set, then for each run, before the LFQ block:
  1. load the temp cache and `Ms1Cache::store_at(&index, dir/stem.cache)`
  2. build `RunPsmRow`s from every rank-1 PSM of the run (a new `select_rank1_all(collector, spectra, queues, candidates, run_idx)`, the same as `select_hits` without the q and decoy filters)
  3. `write_run_quant_table`

  The envelope is `Formula::peptide(...).envelope(4)` normalised to sum 1. `precursor_key` is `format!("{peptidoform}/{charge}")`.
- [ ] **Step 4: Run the tests and see them pass**, plus the whole `quant_end_to_end` suite.
- [ ] **Step 5: Commit** `git commit -m "quant: --quant-dir writes each run's cache and PSM table"`

### Task 4: Library confidence (P0)

**Files:**
- Create: `crates/quant/src/library.rs` (this task: `precursor_qvalues`)
- Modify: `crates/quant/src/lib.rs`

**Interfaces:**
- Produces: `pub fn precursor_qvalues(rows: &[(u32, bool, f64)]) -> Vec<Option<f64>>`. Input
  is `(precursor_id, is_decoy, psm_score)` with higher score better; output is indexed by
  precursor id, with `None` for ids with no rows. A precursor id is all-target or
  all-decoy (decoy peptidoforms are distinct strings).

- [ ] **Step 1: Write the failing tests**

```rust
#[test]
fn best_psm_per_precursor_then_tdc() {
    // p0 target scores 10 and 3 -> best 10; p1 decoy 9; p2 target 8; p3 target 1
    let rows = [(0, false, 10.0), (0, false, 3.0), (1, true, 9.0), (2, false, 8.0), (3, false, 1.0)];
    let q = precursor_qvalues(&rows);
    // raw (D+1)/T down the list: p0 1/1, p1 2/1 -> 1 (clamp), p2 2/2, p3 2/3;
    // the running minimum from the bottom makes every one 2/3
    for i in 0..4 { assert!((q[i].unwrap() - 2.0 / 3.0).abs() < 1e-12, "{i}: {:?}", q[i]); }
}

#[test]
fn many_targets_reach_low_q() {
    let mut rows: Vec<(u32, bool, f64)> = (0..1000).map(|i| (i, false, 1000.0 - i as f64)).collect();
    rows.push((1000, true, 0.5));
    let q = precursor_qvalues(&rows);
    assert!(q[0].unwrap() <= 0.01 && q[999].unwrap() <= 0.01);
}
```

The first test's expected values were worked out by hand; recompute them if the formula changes.

- [ ] **Step 2: Run and see them fail.**
- [ ] **Step 3: Implement:** keep the max score per id, sort `(score desc, id asc)`, accumulate T/D, `fdr = (D+1)/max(T,1)`, take the running min from the end, and clamp to ≤ 1.
- [ ] **Step 4: Pass.** **Step 5: Commit** `git commit -m "quant: precursor-level library confidence from PSM scores"`

### Task 5: Joint RT alignment (P1)

**Files:**
- Create: `crates/quant/src/align.rs`
- Modify: `crates/quant/src/lib.rs`

**Interfaces:**
- Produces:

```rust
pub struct Observation { pub precursor: u32, pub run: u16, pub rt: f64 }
pub struct AlignParams { pub rounds: usize, pub min_anchors: usize, pub min_coverage: f64, pub min_spacing_s: f64 }
impl Default for AlignParams { /* 8, 50, 0.30, 0.5 */ }
pub struct RunMap { knots_x: Vec<f64>, knots_y: Vec<f64>, slope_lo: f64, slope_hi: f64 }
impl RunMap { pub fn apply(&self, t: f64) -> f64; pub fn invert(&self, rt: f64) -> f64; }
pub struct Alignment {
    pub consensus: Vec<f64>,          // NaN when a precursor has no observation in an aligned run
    pub maps: Vec<Option<RunMap>>,    // None = unaligned run
    pub residual_sd: Vec<Option<Vec<(f64, f64)>>>, // per run (consensus t, sd) knots
    pub tau: Vec<f64>,                // per precursor extra dispersion (A3)
    pub n_runs_seen: Vec<u16>,
}
impl Alignment {
    pub fn predict(&self, p: u32, run: u16) -> Option<(f64, f64)>; // (t_hat, s) per A3
}
pub fn align(obs: &[Observation], n_precursors: usize, n_runs: usize, p: &AlignParams) -> Alignment;
```

- [ ] **Step 1: Write the failing tests**

```rust
fn synth(n_p: usize, warps: &[fn(f64) -> f64], noise: f64, seed: u64) -> Vec<Observation> {
    // true t_p uniform in [600, 3600]; rt = warp_r(t_p) + noise * standard-normal (splitmix64 + Box-Muller)
}

#[test]
fn recovers_monotone_warps() {
    let warps: [fn(f64) -> f64; 3] = [|t| t, |t| 1.1 * t + 30.0, |t| t + 60.0 * (t / 3600.0).powi(2)];
    let obs = synth(2000, &warps, 5.0, 1);
    let a = align(&obs, 2000, 3, &AlignParams::default());
    // prediction error for held-out (p, r): the median |t_hat - truth| is below 1.5x the noise
    let errs: Vec<f64> = obs.iter().map(|o| (a.predict(o.precursor, o.run).unwrap().0 - o.rt).abs()).collect();
    assert!(median(&errs) < 7.5, "{}", median(&errs));
}

#[test]
fn robust_to_ten_percent_outliers() {
    // like the test above, but replace 10% of the rts with uniform random values in [0, 4000]; the median error stays < 7.5
}

#[test]
fn maps_are_monotone_and_invertible() {
    // for every aligned map: apply is non-decreasing on a grid, and invert(apply(t)) == t within 1e-6
}

#[test]
fn run_with_few_anchors_is_unaligned() {
    // run 3 has 20 observations -> maps[3].is_none(), predict(_, 3) == None, the other runs still recover
}

#[test]
fn single_run_precursor_gets_consensus_without_anchoring() {
    // precursor seen only in run 1 at rt 2000: consensus finite; removing it leaves the maps bit-identical
}

#[test]
fn transfer_sd_shrinks_with_more_donors() {
    // the s of a precursor seen in 2 runs is greater than the s of a precursor seen in 6 runs (same run, same t)
}
```

- [ ] **Step 2: Run and see them fail.**
- [ ] **Step 3: Implement A3** exactly: equal-count binning, the median per bin, pool-adjacent-violators on the bin y values, linear interpolation, end slopes from the outer 20% of anchors (least squares), and invert by binary search on the knots. Iterate 8 rounds, fix the gauge, then compute the residual width and `tau`. Use a deterministic order everywhere, so the output is the same for the same input.
- [ ] **Step 4: Pass.** **Step 5: Commit** `git commit -m "quant: joint RT alignment across runs without a reference run"`

### Task 6: Library build (P1)

**Files:**
- Modify: `crates/quant/src/library.rs`

**Interfaces:**
- Consumes: Task 4 `precursor_qvalues`, Task 5 `align`.
- Produces:

```rust
pub struct PsmIn { pub run: u16, pub precursor: u32, pub rt: f64, pub score: f64, pub is_decoy: bool, pub faims_cv: Option<f32> }
pub struct PrecursorMeta { pub key: String, pub mz: f64, pub charge: u8, pub envelope: [f32; 4] }
pub struct LibraryEntry {
    pub precursor: u32, pub q: f64, pub consensus_rt: f64,
    pub faims_cv: Option<f32>,      // the most frequent CV among its identifications
    pub donors: Vec<u16>,           // runs where it is identified at library q (any PSM)
    pub best_rt: Vec<(u16, f64)>,   // per donor run, the RT of its best PSM
}
pub struct Library { pub entries: Vec<LibraryEntry>, pub alignment: Alignment, pub meta: Vec<PrecursorMeta> }
pub fn build_library(psms: &[PsmIn], meta: Vec<PrecursorMeta>, n_runs: usize, max_q: f64, ap: &AlignParams) -> Library;
```

Entries are target precursors with `q <= max_q`, sorted by `mz`. Alignment anchors are the
entries' best-PSM RT per donor run.

- [ ] **Step 1: Failing tests:**
  - `decoys_never_enter_library`
  - `entries_sorted_by_mz`
  - `faims_cv_is_majority_of_donor_ids`
  - `best_rt_uses_best_scoring_psm_per_run`, a precursor with two PSMs in one run where the later one scores higher, so `best_rt` is the later RT
- [ ] **Step 2: See them fail. Step 3: Implement. Step 4: Pass.**
- [ ] **Step 5: Commit** `git commit -m "quant: precursor library with global confidence and joint alignment"`

### Task 7: Decoy hypotheses (P1)

**Files:**
- Create: `crates/quant/src/decoys.rs`

**Interfaces:**
- Produces:

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Hypothesis { Target, MassShiftDecoy(i8), RtShiftDecoy }
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Probe { pub kind: Hypothesis, pub mz: f64, pub t_hat: f64, pub s: f64, pub w_prior: f64 }
pub const MASS_SHIFT_UNITS: [i8; 4] = [-11, -7, 7, 11];
pub const PEPTIDE_MASS_SLOPE: f64 = 1.000495;
pub fn probes(precursor: u32, run: u16, mz: f64, charge: u8, t_hat: f64, s: f64, w_prior: f64, run_rt: (f64, f64)) -> Vec<Probe>;
```

`probes` returns the target and one decoy of each family, per A4 (the RT decoy is omitted
when no sign fits). The hash is splitmix64 of `(precursor as u64) << 16 | run`.

- [ ] **Step 1: Failing tests:**

```rust
#[test]
fn mass_shift_stays_in_mass_defect_band() {
    let ps = probes(5, 2, 700.0, 2, 1500.0, 10.0, 6.0, (0.0, 3600.0));
    let d = ps.iter().find(|p| matches!(p.kind, Hypothesis::MassShiftDecoy(_))).unwrap();
    let Hypothesis::MassShiftDecoy(n) = d.kind else { unreachable!() };
    assert!(((d.mz - 700.0) * 2.0 - n as f64 * PEPTIDE_MASS_SLOPE).abs() < 1e-9);
}
#[test]
fn rt_shift_keeps_search_breadth_and_stays_in_run() { /* s and w_prior equal the target's; t_hat' inside the run range; |shift| >= 8s+4w */ }
#[test]
fn deterministic_per_precursor_and_run() { /* the same inputs give identical probes; different runs can pick different n */ }
#[test]
fn no_rt_decoy_when_run_too_short() { /* run_rt = (1490, 1510) -> only the target and the mass-shift decoy */ }
#[test]
fn shifts_avoid_common_deltas() {
    let deltas = [0.984, 14.016, 15.995, 17.027, 18.011, 21.982, 27.995, 31.990, 37.956];
    for n in MASS_SHIFT_UNITS { for d in deltas { assert!(((n as f64).abs() * PEPTIDE_MASS_SLOPE - d).abs() > 0.9); } }
}
```

- [ ] **Steps 2-4.** **Step 5: Commit** `git commit -m "quant: mass-defect-preserving and RT-shifted decoy hypotheses"`

### Task 8: Template-fit evidence (P1)

**Files:**
- Create: `crates/quant/src/evidence.rs`

**Interfaces:**
- Consumes: `Ms1RunIndex::{scan_window, rt, xic}`, Task 7 `Probe`.
- Produces:

```rust
pub struct EvidenceParams { pub tol_ppm: f64, pub n_isotopes: usize /* 4 */ }
#[derive(Debug, Clone, PartialEq)]
pub struct Evidence {
    pub amplitude: f64, pub amplitude_se: f64, pub quantity: f64, pub quantity_se: f64,
    pub llr: f64, pub rt_z: f64, pub width_ratio: f64, pub envelope_cos: f64, pub coherence: f64,
    pub m1_ratio: f64, pub mass_err_ppm: f64, pub n_points: u32, pub isotopes_used: u8,
    pub t0: f64, pub w: f64,
}
pub fn fit(ms1: &Ms1RunIndex, probe: &Probe, charge: u8, envelope: &[f32; 4], p: &EvidenceParams) -> Option<Evidence>;
// None only when the window holds fewer than 5 MS1 scans.
```

`xic` must return the most intense centroid in tolerance per scan. Check
`ms1_index.rs:187`; if it sums, add `xic_max`. The fit takes the centroid m/z for the mass
error from a new `Ms1RunIndex::apex_centroid(idx, lo, hi) -> Option<(f64, f32)>`.

- [ ] **Step 1: Failing tests.** All use synthetic `Ms1RunIndex` built from Gaussian elution
  profiles, written directly in the test module with a helper
  `synth_run(features: &[(mz, z, envelope, apex_rt, width, amplitude)], noise: f32, scan_dt: f64)`:

```rust
#[test] fn recovers_amplitude_apex_and_width() { /* one feature, noise 1% of the apex: |A/A_true - 1| < 0.05, |t0 - apex| < scan_dt, width_ratio within 0.7..1.4 */ }
#[test] fn absent_precursor_has_small_llr() { /* empty run (noise only): llr < 10, amplitude near 0 */ }
#[test] fn interfered_isotope_is_excluded_and_quantity_unbiased() {
    /* add an unrelated co-eluting peak on isotope 2 at 5x: isotopes_used excludes bit 2; |quantity/true - 1| < 0.1 */
}
#[test] fn rt_prior_prefers_peak_near_prediction() { /* two equal peaks at t_hat+-0 and t_hat+5s: t0 picks the near one */ }
#[test] fn m1_ratio_flags_non_monoisotopic_hypothesis() { /* a probe at the +1 isotope of a real feature: m1_ratio > 0.5; the true mono probe: m1_ratio < 0.1 */ }
#[test] fn decoy_probe_in_empty_region_scores_below_target() { /* target present: llr(target) > 20 * max(1, llr(mass-shift decoy)) */ }
#[test] fn extract_uses_donor_faims_group() {
    /* two FAIMS groups; the feature exists only at -45: fit on the -45 index finds it, fit on the -60 index has llr < 10 */
}
```

- [ ] **Step 2: See them fail. Step 3: Implement A2** exactly (formulas above). Keep each step a small private function:
  - `traces`
  - `background_and_noise`
  - `best_alignment` (the grid over t0 and w)
  - `exclude_interfered`
  - `features`
- [ ] **Step 4: Pass**, then run `cargo bench`-free timing in the test log: fitting 10k probes on a synthetic 3,000-scan run prints the elapsed time (no assert).
- [ ] **Step 5: Commit** `git commit -m "quant: likelihood template fit for (precursor, run) evidence"`

### Task 9: Extract stage (P1)

**Files:**
- Create: `crates/andes/src/bin/andes/quant_mbr.rs` (the stage functions shared by the subcommands and by all-in-one)
- Create: `crates/output/src/extract_table.rs` (writer/reader for `<stem>.extract.parquet`, schema `andes.extract.schema = "1.0"`)
- Modify: `crates/andes/src/bin/andes.rs` (module declaration)

**Interfaces:**
- Consumes: Tasks 1-8.
- Produces:

```rust
// quant_mbr.rs
pub(crate) struct ExtractRow {
    pub precursor: u32, pub run: u16, pub kind: quant::decoys::Hypothesis,
    pub role: Role, // Identified | Transfer | PseudoTransfer
    pub evidence: Option<quant::evidence::Evidence>,
    pub n_donor_runs: u16, pub library_q: f64,
}
pub(crate) enum Role { Identified, Transfer, PseudoTransfer }
pub(crate) fn extract_run(run: u16, cache: &Path, lib: &quant::library::Library, ep: &EvidenceParams) -> Result<Vec<ExtractRow>, String>;
```

Per run:
1. Load the cache and split it by FAIMS (reuse `faims_groups`).
2. Estimate `w_prior(t)` (A7) from the identified pairs.
3. For every library entry:
   - **Identified in this run:** one `Identified` target fit, with the prior at its best
     PSM RT, plus one `PseudoTransfer` target fit, with the prior from `alignment.predict`
     computed with this run's own observation left out. That means
     `t_hat = map_r(consensus)` and `s` from A3 with `n_p - 1`.
   - **Not identified:** a `Transfer` target and its decoy probes.
   - **Run unaligned, or no prediction:** skip and count.
4. Use the entry's FAIMS group; when the run has no such group, use the merged index.
   Fits run in parallel over entries with rayon. The output order is the library order.

- [ ] **Step 1: Failing tests** (unit tests in `quant_mbr.rs` with a synthetic two-run library over `synth_run` caches written to a tempdir):
  - `identified_pairs_get_identified_and_pseudo_rows`
  - `missing_pairs_get_target_and_decoys`
  - `unaligned_run_is_skipped_with_count`
  - `output_order_is_deterministic`, two calls give equal vectors
- [ ] **Steps 2-4.**
- [ ] **Step 5: Commit** `git commit -m "quant: per-run extraction of every library precursor"`

### Task 10: Transfer FDR and finalize (P1)

**Files:**
- Create: `crates/andes/src/transfer_fdr.rs` (in the `andes` lib crate next to `rescore.rs`; the CV pattern follows `rescore.rs::cv_scores`)
- Modify: `crates/andes/src/lib.rs`
- Modify: `crates/andes/src/bin/andes/quant_mbr.rs` (`finalize`)
- Create: `crates/output/src/mbr_out.rs` (`<base>.mbr.tsv`, a long format, one row per (precursor, run))
- Modify: `crates/output/src/feature_parquet.rs` (optional columns `evidence_type`, `transfer_q`, `quantity_se`; absent when not MBR)

**Interfaces:**
- Produces:

```rust
pub struct FdrRow { pub precursor: u32, pub is_decoy: bool, pub family: u8 /* 0 target, 1 mass, 2 rt */, pub x: Vec<f32> }
pub struct FdrResult { pub q: Vec<Option<f64>>, pub q_mass_only: Vec<Option<f64>>, pub q_rt_only: Vec<Option<f64>>, pub score: Vec<f32> }
pub fn transfer_qvalues(rows: &[FdrRow], seed: u64) -> FdrResult;
pub const TRANSFER_FEATURES: &[&str]; // names in A5 order, written to the log
```

`finalize` gathers every run's `ExtractRow`s and computes the A5 features:
- **Transfer and decoy rows** go to `transfer_qvalues`.
- **`PseudoTransfer` rows** are scored by the fitted fold models (the fold of their
  precursor). This gives the sensitivity and wrong-peak rate in A5, both logged and
  written as a two-line `<base>.mbr_validation.tsv`.

Output rows:
- **Identified:** quantity, SE, `evidence_type = identified`, `transfer_q = NA`.
- **Transfer with q ≤ `--mbr-fdr`** (default 0.01): `evidence_type = transfer_ms1`.
- **Other transfers:** written with `evidence_type = below_threshold`, quantity and SE
  kept, so downstream tools can treat them as censored rather than missing.

- [ ] **Step 1: Failing tests** (`transfer_fdr.rs`):

```rust
#[test] fn separable_targets_get_low_q() { /* 500 targets x ~ N(3,1), 1000 decoys x ~ N(0,1) (K=2): >= 400 targets at q <= 0.01 */ }
#[test] fn null_targets_are_not_accepted() { /* targets and decoys from the same distribution: targets at q <= 0.01 are <= 5 */ }
#[test] fn folds_group_by_precursor() { /* every precursor's rows share one fold id */ }
#[test] fn fdr_degenerate_inputs_yield_none() { /* no decoys -> all None; no targets -> empty; one-class fold -> those rows None, no panic */ }
#[test] fn per_family_estimates_reported() { /* family-1-only decoys: q_rt_only is None for all rows */ }
```

- [ ] **Steps 2-4.**
- [ ] **Step 5: Commit** `git commit -m "quant: cross-validated transfer FDR, pseudo-transfer validation, MBR outputs"`

### Task 11: CLI wiring: subcommands and all-in-one (P1)

**Files:**
- Modify: `crates/andes/src/bin/andes/cli.rs`, `crates/andes/src/bin/andes.rs`, `crates/andes/src/bin/andes/search.rs`, `crates/andes/src/config.rs`
- Modify: `DOCS.md`, `README.md`, `config.example.yaml`
- Test: `crates/andes/tests/quant_mbr_end_to_end.rs`

**Interfaces:**
- `andes quant library --quant-dir D --psms T.tsv --decoy-psms D.tsv [--mbr-library-fdr 0.01]` writes `D/library.parquet` (entries + meta) and `D/alignment.parquet`. The PSM tsvs are Percolator outputs, joined to the run tables by `SpecId`. A Percolator `score` column is required; when it is absent, the error says so.
- `andes quant extract --quant-dir D --run <stem>` writes `D/<stem>.extract.parquet`.
- `andes quant finalize --quant-dir D -o <base> [--mbr-fdr 0.01] [--output-parquet DIR]`.
- `andes search ... --lfq --mbr` with two or more inputs: implies `--quant-dir` (a temp dir when not given) and needs rescoring (`--rescore` or `--rescore-native`; error otherwise). After quant, it runs library, extract (runs in parallel) and finalize in process, using the in-memory Percolator map in place of the tsvs.
- Help text labels `--mbr-library-fdr` as "internal library entry gate (precursor-level), not a reported identification FDR" and `--mbr-fdr` as "transfer quality q-value".

- [ ] **Step 1: Failing end-to-end test.**
  - Fixture: a new script, `scripts/make_quant_fixtures.py --mbr`, writes three synthetic runs from the six-peptide BSA design with RT warps `t`, `1.05 t + 20`, `t + 15 sin(t/600)`. In run 3, peptide 4 is never triggered for MS2: its MS1 feature exists but there is no MS2 scan. In all runs, peptide 6 is absent from MS1 and only identified in run 1.
  - Assertions:
    - all-in-one and subcommand paths give byte-identical `.mbr.tsv`
    - peptide 4 in run 3 is `transfer_ms1` with quantity within 10% of the designed value
    - peptide 6 in runs 2 and 3 is not `transfer_ms1`
- [ ] **Steps 2-4.**
- [ ] **Step 5: Commit** `git commit -m "quant: andes quant library/extract/finalize and --mbr all-in-one"`

### Task 12: Speed: RT sweep extraction (P1, only if Task 13 profiling shows extract dominating)

**Files:** `crates/quant/src/evidence.rs` (`fit_many(ms1, probes, ...) -> Vec<Option<Evidence>>`), `quant_mbr.rs`.

Gate: the profile from Task 13 on UPS1 shows `extract` above 30% of the wall time of the
search. If not, skip this task and record the profile.

Method: sort probes by `mz` and keep an active set as the RT cursor advances. For each MS1
scan, a merge-join of its m/z-sorted peaks with the active probes' isotope m/z values
fills the trace buffers. Then fit each probe from its buffer.

- [ ] **Step 1: Failing test:** `fit_many_equals_fit` on 2,000 random probes over a synthetic run (exact equality of `Evidence`).
- [ ] **Steps 2-4**, plus timing before/after on the UPS1 run.
- [ ] **Step 5: Commit** `git commit -m "quant: single RT sweep for all probes of a run"`

### Task 13: Benchmark and gates (P1)

**Files:** benchmark scripts only, on the build machine (scripts are kept; data and outputs are deleted afterwards). Results go into the report, not into the repo.

- [ ] **Step 1: P0 gate.** On the six UPS1 runs (PXD001819), with `--quant-dir` vs without: `.lfq.tsv` byte-identical, and record the cache size per run. Repeat the cache-size measurement on the Astral HYE run. Gate: the cache is at most 50% of the raw file size; otherwise report it and propose f32 m/z.
- [ ] **Step 2: UPS1 MBR.** All-in-one `--lfq --mbr --rescore`.
  - Report:
    - missing fraction (yeast, UPS1)
    - the number of UPS1 precursors with a fold change
    - fold-change median and MAE (true 3.32 / 0)
    - yeast CV, against the PR #119 baseline and Sage 0.14.7 run with its own MBR
  - Report quantity values both ways: identified + `transfer_ms1` only, and also including `below_threshold`.
- [ ] **Step 3: Transfer FDR truth.** Find a DDA dataset with a species missing from some runs by design (search PRIDE for two-species LFQ benchmarks with dropout or blank samples; check the file listing across pages). Measure the true transfer FDP as transfers of the absent species, after removing peptides shared between species, divided by all accepted transfers into those runs. Also report the per-family q estimates and their agreement. Gate: true FDP ≤ 1.5 × nominal at 1%. If no such dataset fits on the build machine, report the per-family agreement and the pseudo-transfer validation only, and state that the true FDP is unmeasured.
- [ ] **Step 4: Validation numbers.** Pseudo-transfer sensitivity at q ≤ 0.01 and the wrong-peak rate.
- [ ] **Step 5: Speed.** Wall time and peak RSS of each stage on UPS1, and the share of extract (decides Task 12).
- [ ] **Step 6: One-variable A/Bs, each kept only if it wins.**
  - the isotope-exclusion rule (on/off)
  - the RT prior inside the objective vs a hard window
  - K = 2 vs mass-shift-only
  - the width grid (fixed `w_prior` vs 3-point)
- [ ] **Step 7: Clean up** the benchmark data and outputs; keep scripts.

---

## Deferred to later plans (from the spec, intentionally not in P0/P1)

- Reduced MS2 in the run cache and MS2-checked transfers: P2. The cache major version moves to 2 then.
- Cross-run signal selection with variance output: P3. A2 already writes `quantity_se`.
- TMT isolation-window capture, reporter noise and unmixing: P4.
- m/z sharding of library/extract, and quantms wiring: P5.
