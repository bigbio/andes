# Native rescoring pipeline for msgf-rust — design (DRAFT, pending review)

> **STATUS: DRAFT for review.** Produced by the brainstorming skill from two
> feasibility investigations (2026-05-31). NOT approved; NO implementation until
> the user reviews and answers the open decisions in §7. Do not transition to
> writing-plans until approved.

## 1. Goal

Let a user go from spectra to FDR-controlled q-values in **one `msgf-rust`
command**, optionally enriched with deep-learning rescoring features — removing
the two external middle layers we use today: the Python **MS²Rescore** stack
(MS²PIP + DeepLC) and the **Docker Percolator** step.

Today: `msgf-rust` → `.pin` file → (Docker `percolator`) → q-values; and any
ML rescoring requires a separate Python MS²Rescore run on the side.

Target: `msgf-rust --rescore --percolator` → q-values, with predicted-spectrum
and predicted-RT features folded in.

## 2. Scope decomposition (two independent subsystems)

These are independent and should each get their own spec + plan:

- **S1 — Native Percolator integration.** Make the PIN flow directly into
  Percolator from `msgf-rust` (one command), no Docker, no temp file required.
  Works on ANY pin; orthogonal to rescoring.
- **S2 — ML rescoring features.** Compute MSBooster/MS²Rescore-style features
  (observed-vs-predicted fragment-intensity similarity; observed-vs-predicted
  RT delta) and emit them as NEW additive `.pin` columns. Feeds ANY Percolator.

This design covers both and a recommended build order; the first implementation
plan should target a single sub-project (see §7 Q1).

## 3. Key findings that shape the design

**Percolator (S1):**
- Apache-2.0. Reads pin-tab from a **file OR stdin** (`percolator -` /
  `--stdinput`). Building with the XML/converter path disabled drops the
  Xerces-C dependency (only Boost + Eigen + internal libs remain).
- `crates/output/src/pin.rs` already writes via a `Write`-generic
  `write_pin_to<W: Write>` — it can stream to a child process's stdin with no
  change, so an in-memory pin → `percolator -` pipe needs no temp file.
- A static `perclibrary` core exists (FFI is *possible*) but exposes C++
  classes, not a C ABI; FFI buys the same parity as a subprocess at much higher
  build cost (Boost/CMake inside `cargo build`). Not worth it.
- A Rust reimplementation (sage ships a built-in LDA "within ~1% of Percolator")
  is the only true zero-external-dep path but will NOT reproduce Percolator
  3.7.1 exact 1%/5% counts — incompatible with our current benchmark gate
  without a deliberate re-baseline.

**ML rescoring (S2):**
- The features (Half B) are independent of how predictions are obtained
  (Half A) and are pure Rust: a spectral-similarity metric (normalized spectral
  angle / dot product over matched b/y intensities) and an RT-delta (with a
  per-run alignment/calibration variant). Emit as additive pin columns — the
  one historically-safe change class (cf. C-4, iter19 EdgeScore; does not
  perturb top-1 selection, so it cannot trigger the T/D-ratio canary).
- Predictions (Half A): **Koina** (wilhelm-lab, Nat Commun 2025) is an
  open, self-hostable Triton inference server exposing HTTP/gRPC (KServe v2),
  serving MS²PIP + DeepLC + Prosit + AlphaPeptDeep. Crucially the **server does
  the peptide feature-engineering** — the client sends sequences/charges/CE and
  gets back mz/intensities or predicted RT. This avoids porting MS²PIP's Cython
  per-fragment encoder (the multi-week, high-parity-risk part).
- Native embedding is feasible but asymmetric in risk: DeepLC is a small CNN
  with a clean atom-composition encoding → ONNX (`tract`/`ort`) is low-risk;
  MS²PIP is XGBoost trees (loadable via pure-Rust `gbdt` or native `xgboost`)
  but its **Cython feature extractor** (AA basicity/hydrophobicity/helicity/pI
  tables + cleavage-position windowing + mod handling) must be ported
  bit-faithfully — the highest-risk piece, deferred to a later phase behind a
  per-fragment parity harness.
- Precedent: Sage deliberately does NOT embed these models; it uses a built-in
  LDA + cheap RT model and delegates heavy ML rescoring to external MS²Rescore.

## 4. Recommended architecture (phased, opt-in)

```
msgf-rust search ─► base PIN (in-memory, Write-generic writer)
                       │
        ┌──────────────┴──────────────┐  S2 (opt: --rescore)
        │ ML rescoring features        │
        │  predictions ◄── Koina HTTP  │  (Phase 1b)
        │  spectral-angle + RT-delta   │  → additive PIN columns (Phase 0)
        └──────────────┬──────────────┘
                       │
        ┌──────────────┴──────────────┐  S1 (opt: --percolator)
        │ stream PIN → `percolator -`  │  (Phase 1a)
        └──────────────┬──────────────┘
                       ▼
                  q-values TSV  ← one command
```

All stages are opt-in flags; with neither flag, behavior is byte-identical to
today (writes the same `.pin`).

## 5. Build order (smallest-viable-first)

- **Phase 0 — Rust feature math (S2 Half B).** Implement spectral-angle and
  RT-delta computation over already-available observed matched ions + a supplied
  prediction, emitted as new additive pin columns. Validate against a Koina
  prediction for a handful of PSMs. Pure Rust, low risk. *The features are the
  point.*
- **Phase 1a — `--percolator` (S1).** Spawn a bundled `percolator -`, stream the
  in-memory pin to its stdin via the existing `write_pin_to`, parse q-value TSV
  back. Perfect parity. Independent — could ship first. Ship the binary as an
  optional separate download initially (don't gate the release matrix on it).
- **Phase 1b — Koina client (S2 Half A).** `reqwest` + `serde_json` against a
  self-hostable Koina endpoint; batch peptides; feed Phase-0 features. Removes
  the Python MS²Rescore layer with minimal parity risk.
- **Phase 2 (later/optional) — native model embedding.** DeepLC via ONNX
  (`tract`) first (low encoder risk); MS²PIP native (`gbdt`/`xgboost`) last,
  gated behind a per-fragment bit-parity harness vs. Python MS²PIP. Enables
  offline single-binary operation. Only if §7 Q4 says offline is required.

## 6. Approaches considered (overall shape)

1. **Phased Koina + bundled Percolator (RECOMMENDED).** Fastest to a working
   one-command pipeline; lowest parity risk; self-hostable for offline/repro.
   Cost: a network/Triton service dependency (self-host for benchmark scale)
   and a per-platform Percolator binary in the release.
2. **Full-native single binary.** Best end-state (offline, one static binary,
   trivial CI). Highest effort/risk — the MS²PIP Cython encoder port and a
   Percolator reimplementation/re-baseline. A destination, not a first step.
3. **Subprocess-everything (interim oracle only).** Shell out to MS²Rescore +
   Percolator Python/Docker. Lowest effort, perfect parity, but reintroduces
   the exact Python stack we want to remove. Use only to validate 1/2.

## 7. Open decisions (need user input before a plan)

1. **First spec scope:** S1/Phase-1a (native Percolator) alone, or the combined
   Phase 0+1a+1b rescoring pipeline?
2. **Percolator distribution:** bundle/vendor the binary (subprocess + stdin), or
   commit to a sage-style Rust reimplementation (accepting a benchmark
   re-baseline)?
3. **Predictions source first:** Koina (self-hosted) or straight to native
   ONNX/XGBoost embedding?
4. **Offline requirement:** is single-binary, no-network operation a hard
   requirement (→ phase 2 native) or nice-to-have (→ Koina is fine)?

## 8. Testing strategy (outline)

- Phase 0: unit tests for spectral-angle/RT-delta on synthetic
  observed-vs-predicted vectors; golden values vs. a manual computation.
- Phase 1a: integration test that `--percolator` on a fixture pin reproduces the
  same target counts as the Docker Percolator 3.7.1 path (parity gate).
- Phase 1b: contract test against a recorded Koina response (no live network in
  CI); a separate opt-in live-endpoint smoke test.
- Bench: re-run the 3-dataset entrapment harness with `--rescore` to measure the
  sensitivity lift and confirm FDP stays nominal.

## 9. Non-goals

- Retraining MS²PIP/DeepLC models (inference/reuse only).
- Ion-mobility features (timsTOF) — out of scope for the first cut.
- Replacing the chimeric cascade (orthogonal; rescoring features would *also*
  feed cascade secondaries).

## References

See the two investigation reports (Percolator integration; MS²PIP+DeepLC Rust
feasibility) — key sources: percolator/percolator (Apache-2.0, `--stdinput`),
mokapot, linfa-svm, sage (built-in LDA), MSBooster (Nat Commun 2023), MS²Rescore
3, Koina (Nat Commun 2025, wilhelm-lab), gbdt/ort/tract crates, compomics ms2pip
v4 (XGBoost models on Zenodo) + DeepLC (Keras HDF5 → ONNX).
