# Hardening + Unimod-resolver implementation plan

Status: plan. Combines (a) the 7-crate Codex adversarial review (2026-06-24, every
crate `needs-attention`) and (b) the Unimod-accession→mass gap into one prioritized
track. **Independent of the geometry PR (#13)** — that tests a different question and
is gated on the benchmark; this is a parallel correctness/robustness track.

## Guiding principle — order by blast radius

The reviews surfaced three recurring root causes, not a wrong core design:
1. **Fail-silent instead of fail-loud** — corrupted inputs, schema-mismatched models,
   failed Percolator joins, missing artifacts all proceed via NaN / skip / WARN / exit-0.
2. **No machine-enforced contracts** — feature counts, source weights, GBDT quality,
   SpecId uniqueness, CLI value domains are unchecked.
3. **A few genuine correctness holes in opt-in paths** — can change results or bias FDR.

The core algorithm (rank/strong scoring, expand-then-discriminate, Percolator-only FDR,
the model store) was **not** challenged. So this is hardening + one feature, ordered so
the highest-FDR-impact work lands first.

Already done (geometry review): deterministic `main_ion` under tied/zero frequencies
(`7d24bbcc`, scoring) — folded into the same family as Phase 3's finite-LLR guard.

---

## Phase 1 — FDR-boundary correctness (highest blast radius; can corrupt *reported* FDR)

These can silently misreport q-values to users. Do first.

| # | Finding | Location | Fix | Test |
|---|---------|----------|-----|------|
| 1.1 | PIN `SpecId` not globally unique → Percolator `HashMap<PSMId>` overwrite / lost PSMs | `output/src/pin.rs:102` | Add a stable run-local spectrum ordinal to `SpecId`; track a `HashSet` while writing and **error on duplicates**; reject duplicate `PSMId` when parsing Percolator output | dup/missing title+scan, `scan=0_0_1` collapse, dup PSMId |
| 1.2 | QPX silently accepts missing Percolator joins (null PEP == "off" == "join failed") | `output/src/qpx.rs:400` | When `rescore.is_some()`, require every non-decoy emitted row to have a matching Percolator result; else error with missing/extra key counts | partial-join bundle |
| 1.3 | Percolator output stale-read / race (fixed filenames in PIN dir, no unlink/validate) | `output/src/percolator.rs:245` | Write results to a fresh temp/unique-run-id dir; remove expected outputs pre-run; reject user `--percolator-args` that redirect result/weight paths; validate parsed PSMIds ⊆ current PIN | concurrent same-stem, stale file present |
| 1.4 | Malformed q-value/PEP → `NaN` instead of error | `output/src/percolator.rs:331` | Return parse error (line + column) for short rows, non-finite q/PEP, duplicate PSMId | truncated/locale/garbage result row |
| 1.5 | Refinement Pass-2 becomes **target-only** under `DecoyStrategy::None` → no paired decoy competition → FDR bias | `search/src/refinement.rs:244` | For the synthetic anchor mini-DB always generate fresh scoped decoys with a real strategy; OR reject `--refine` + `DecoyStrategy::None` unless anchor decoys are supplied with suffix labeling preserved | refine + external-decoy FASTA |

## Phase 2 — opt-in feature-combination correctness (changes results for specific flags)

| # | Finding | Location | Fix | Test |
|---|---------|----------|-----|------|
| 2.1 | **[critical]** `chimeric + mmap` indexes empty `prepared.candidates`/`bucket_index` → panic or missing secondary candidates | `search/src/match_engine.rs:1292` | Near-term: **reject `chimeric && Mmap` up front** with a clear error. Long-term: lazy `mmap_index` query per co-isolated mass + resolve primary via the mmap accumulator before Pass-2 | chimeric+mmap small DB |
| 2.2 | mmap candidate-backing → order-dependent `strong_score`/primary selection; `--candidate-index auto` makes output machine-dependent | `search/src/match_engine.rs:79` | Make listwise/null features + `could_win` pruning order-independent (compute from the full scored distribution or a deterministic global order); shared stable tie-breakers for RAM and mmap; until then, document that `auto` is not bit-stable and keep FDR-sensitive runs on RAM | RAM vs mmap PSM/PIN equivalence |
| 2.3 | **[high]** Fixed terminal mods don't stack with residue variants → wrong precursor/fragment masses for **TMT/iTRAQ** | `model/src/aa_set.rs:450` | Compose mandatory fixed terminal mods with every applicable base/variable residue variant; exclude unmodified base variants when a fixed terminal mod applies | fixed `* N-term` + variable Ox-M; fixed N-term + fixed K on N-term K |
| 2.4 | Semi-specific (`ntt==1`) digestion ignores `max_missed_cleavages` → larger candidate/null space than intended | `search/src/candidate_gen.rs:170` | Define + enforce missed-cleavage semantics for semi-specific spans (count internal sites in the emitted span); or warn/reject when a finite limit will be ignored | semi-tryptic + finite missed-cleavage |
| 2.5 | `.mzML.gz` mis-routed: global `is_mzml/is_raw/is_d` from `paths[0].extension()` (no `.gz` strip) drives model/chimeric/calibration; multi-file inherits file[0]'s class | `andes/src/bin/andes.rs:1375` | One canonical `.gz`-aware format detector for all paths incl. the global; for multi-file, validate homogeneous routing class or route per-file explicitly | `sample.mzML.gz`; heterogeneous multi-`--spectrum` |

## Phase 3 — machine-enforced contracts (fail-loud on corrupt artifacts)

| # | Finding | Location | Fix | Test |
|---|---------|----------|-----|------|
| 3.1 | GBDT blobs evaluated without feature-contract check → wrong-slot/schema-skew blob scores plausibly via NaN→`default_left` | `scoring/src/gbdt_eval.rs:138` | Validate model kind + exact feature count at load / first predict: peak==`peak_features::N_FEATURES`, frag==`frag_features::N_FRAG_FEATURES`, rich==`ion_features::N_ION_FEATURES`; reject mismatch | mismatched feature count, swapped slot |
| 3.2 | Rank LLR `ln(ion/(noise·norm))` → ±inf from zero/negative/NaN frequencies | `scoring/src/scoring/rank_scorer.rs:94` | Validate rank distributions in `RankScorer::new`: finite, non-negative, positive denominators; documented smoothing floor; fail setup when LLRs can't be finite | zero-noise, zero-ion, NaN freq table |
| 3.3 | Dense spectra → `prob_peak > 1` → negative baseline → NaN edge score → `round() as i32` = 0 (silent loss of edge evidence) | `scoring/src/scoring/rank_scorer.rs:291` | Make `prob_peak` a valid probability before use, or branch on out-of-domain values with a bounded documented score + instrumentation | `prob_peak > 1` edge output |
| 3.4 | Strong/intensity features → NaN/inf from unbounded regressor `exp` and `inf/inf` | `scoring/src/scoring/strong_score.rs:190` | Clamp log-intensity to a finite calibrated range before `exp`; `is_finite` guards on all derived features | extreme leaf values |
| 3.5 | Unchecked `SourceLedger.weight` → NaN→0 (sources erased) or huge→`u64::MAX` saturation | `model-train/src/store/update.rs:140` | Validate `weight.is_finite() && weight >= 0.0` (+ upper bound) at ingestion and before `CountStats::scaled`, from API and Parquet | negative/NaN/inf weights |
| 3.6 | GBDT trainers log failed validation but still **return a deployable model** (empty-tree, negative R², single-class) | `model-train/src/gbdt/train.rs:391` | Return `Result<GbdtPeakModel, TrainError>` with hard gates: min rows, both classes, min AUC/Pearson/R², explicit empty-tree failure unless opted-in fallback | degenerate/no-signal training |
| 3.7 | Partitioned store unions all `*.parquet`; first-match wins on duplicate `model_id` → path-order-dependent load | `model-train/src/store/read.rs:330` | Build a `model_id → part` index at `open`; reject duplicate manifests/rows across parts unless full manifest + table checksums match; ignore temp suffixes / require the Hive `protocol=*/` layout | dup id across parts, stale temp file |
| 3.8 | CLI range/tol parsers validate syntax+ordering but not domain (`--charge 0..0`, `--precursor-tol 0/neg/NaN/inf`, FDR/PEP outside `[0,1]`) — **from the #12 minimization** | `andes/src/bin/andes.rs:4724` | Domain validators at clap-parse time: charge positive + within supported bounds; tolerances finite `> 0`; FDR/PEP/refine-FDR finite in `[0,1]` | each invalid value rejected with message |
| 3.9 | `--output-parquet` write failure downgraded to WARN → **exit 0** with the requested artifact missing | `andes/src/bin/andes.rs:2420` | Treat failure to write any *requested* output as a hard (non-zero) error; expose best-effort side artifacts explicitly if ever wanted | unwritable `--output-parquet` |

## Phase 4 — input fail-loud (corrupt spectra surface as errors, not silent drops)

| # | Finding | Location | Fix | Test |
|---|---------|----------|-----|------|
| 4.1 | Default mzML iterator discards malformed spectra + resyncs, never yields `Err` | `input/src/mzml.rs:1058` | Make the default iterator fail-visible (yield `Err` before resync) or move tolerant resync behind an explicit API returning an error count/sample | bad base64/zlib/length-mismatch scan |
| 4.2 | Truncated mzML binary array silently rounded down (`UnexpectedEof` == clean end) | `input/src/mzml.rs:1439` | Require decoded byte length to be an exact multiple of precision; error on remainder; prefer mzML `arrayLength` | 1–7 byte truncation |
| 4.3 | Thermo `.raw` peaks bypass the finite/positive sanitization other readers apply; `zip` truncates mismatched arrays | `input/src/thermo.rs:238` | Apply the mzML/timsTOF peak validation; make array-length mismatch observable | NaN/negative/ragged Thermo peaks |
| 4.4 | MGF accepts non-finite/negative peak values | `input/src/mgf.rs:178` | Reject/filter unless m/z finite+positive and intensity finite+non-negative | NaN/inf/negative MGF peaks |

## Phase 5 — Unimod-accession → mass resolver (the feature) + peptide mod-identity

Today andes **cannot consume a Unimod accession or mod name**: `--mods` requires a numeric
delta as field 1, and `accession=UNIMOD:NN` is a carried-through *output label* only
(`model/src/modification.rs:138,195`). The bundled `resources/unimod.obo` is **not parsed
anywhere** in `crates/`. So the caller (e.g. the quantms andes module via `andesModLine`)
must pre-resolve every mass. This phase makes andes self-sufficient.

### 5.1 Unimod OBO loader (`crates/model`)
- New `model::unimod` module: parse `unimod.obo` once into a map
  `UNIMOD:NN → { mono_mass, name, sites }` and `name(lowercased) → UNIMOD:NN`.
  Use the bundled path resolution that `models.parquet` already uses (binary-relative
  `resources/`), lazy + cached (`OnceLock`).
- OBO is a simple stanza format (`[Term] id: UNIMOD:35 name: Oxidation
  xref: delta_mono_mass "15.994915" …`); parse `id`, `name`, `delta_mono_mass`, and the
  `delta_composition` for cross-check. No new dependency — hand-rolled line parser.
- Validation test: a curated set (Carbamidomethyl/Ox/Acetyl/TMT/iTRAQ/Phospho/Deamidated/
  GG) resolves to the values pinned in `model/tests/common_mod_masses.rs` (reuse that
  reference). Tolerance-compare at 1e-5; fail on missing/ambiguous.

### 5.2 `--mods` accepts an accession **or** a bare name in field 1 (`model::modification`)
- Extend `Modification::from_mods_txt_line`: if field 1 doesn't parse as `f64`, try
  `UNIMOD:NN` then a bare name against the loader; set `mass_delta` from the resolved mono
  mass and auto-populate `accession`/`name` if absent. A numeric field 1 keeps current
  behavior exactly (back-compat).
- Errors are explicit: `UnknownUnimodAccession`, `AmbiguousModName`, with the offending
  token. No silent fallback.
- The loader is injected (not a global) so tests stay hermetic and a `--no-unimod` / env
  escape hatch can disable it.

### 5.3 Peptide mod-identity round-trip (folds in the `model` review finding)
- `model/src/peptide.rs:99` — `Display` writes 5-dp rounded mass; `from_str` parses only on
  exact bit-match, so `57.021464`→`+57.02146` no longer round-trips. Fix: serialize the
  **accession** when present (`[UNIMOD:4]`, ProForma-style) and parse it back via the
  resolver; for accession-less mods, match deltas with a documented tolerance and reject
  ambiguous matches. Aligns the peptide string with the QPX `[UNIMOD:NN]` peptidoform
  already in `output/src/qpx.rs:524`.
- Add a round-trip property test: build → `Display` → `from_str` → equal, for the curated
  mod set and for accession-less numeric mods.

### 5.4 quantms simplification (follow-up, separate repo)
- Once andes resolves names/accessions, the quantms andes module's `andesModLine`
  name→mass table can be replaced by passing the mod name/accession straight through.
  Track as a quantms-side change; do **not** block this plan on it.

### 5.5 Enzyme bond-context (related `model` medium, optional)
- `model/src/enzyme.rs:31` `is_cleavable_after` sees only the left residue → can't model
  trypsin-not-before-Pro. Change the predicate to `cleaves_between(left, right)` and encode
  K/R-not-before-P where intended. Defaults must preserve current behavior unless a new
  enzyme variant opts in. Lower priority; bundle with Phase 5 since it's also a chemistry-rule
  fidelity item.

---

## Sequencing & packaging

- **PR A — FDR-boundary** (Phase 1): smallest, highest-value, ship first. Pure correctness.
- **PR B — opt-in-combo guards** (Phase 2): the critical `chimeric+mmap` reject is one line and
  should land with A or immediately after; the rest follow.
- **PR C — contracts/fail-loud** (Phase 3 + 4): a consistent "validate-and-error" sweep; larger
  but mechanical. Split 3/4 if review gets heavy.
- **PR D — Unimod resolver** (Phase 5): a self-contained *feature* with its own tests; can land
  any time, independent of A–C.

Each finding gets an adversarial regression test (the reviews' "Next steps" are the test list).
Guard against churn: every fix is additive validation or an explicit error — no change to the
happy-path scoring/FDR numbers, so the bundled-model and parity tests must stay green.

## Explicitly out of scope here
- The geometry derivation (#13) — separate, benchmark-gated.
- Any change to the core scoring math, the expand-then-discriminate design, or Percolator-only
  FDR — the review did not challenge these.
