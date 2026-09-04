# andes — stronger own models + public release (design)

Date: 2026-06-15
Status: design (brainstorming) → for review before writing-plans
Location: internal (Claude/planning docs live in `internal-docs/`, never in the public repo)

## Goal

Publish the reanalysis-hardening + multi-enzyme work to the public `main`, **gated on
andes shipping its own, independently-trained models that demonstrably beat the
open-source field**. The model work is the centerpiece: a stronger primary scorer
whose **`RawScore` and `DeltaRawScore` are powerful Percolator inputs**.

**Differentiated thesis (from the research):** andes wins by being the one engine with
strong, CPU-only, **non-NN** intensity/noise modeling that is *tolerance-robust across
low-res CID/TMT and high-res Astral*. Sage sums raw intensity (weak low-res); Comet needs
per-mode bin retuning; ProSE is high-res-only; MSFragger needs a neural net (MSBooster/
Prosit) for its strong intensity feature. A gradient-boosted per-peak signal/noise model
fills that gap. Existence proof: MS2PIP (a GBDT) matched/beat Prosit-rescoring on IDs at
1% FDR.

## Scope boundaries (what this is NOT)

- **Not an E-value / spectral-probability model.** The deliverable is a stronger
  `RawScore`/`DeltaRawScore`, not a per-spectrum significance score.
- **Not rescoring.** Prosit/RT/ion-mobility/intensity-similarity features are Percolator's
  post-processing layer (optionally with external tools) — out of scope here.
- **No neural networks** (explicitly not Prosit/DeepLC-style).
- **No fragment-ion index rewrite.** andes already explored + refuted a fragment-index
  "speed v2"; not revisited.
- **No peptide-conditioned intensity prediction model** (the MS2PIP/Pep2Prob task). The
  GBDT here is **peptide-agnostic per-peak signal/noise** only.

## Decisions (locked with the user, 2026-06-15)

1. The own models are **fully independent** (trained from harvested public PSMs; the
   MS-GF+ seed is an init/prior, not derived IP).
2. The public store ships **own slugs only** ("clean independence"); uncovered combos warn
   + need `--model` (already implemented).
3. Independence claim: **full independence, minimal heritage** (MS-GF+ → brief citation).
4. Model work scope: **GBDT peak signal/noise + cleanly retrained rank core**, focused on a
   strong `RawScore`/`DeltaRawScore`.
5. **Codon = training/data; VM = benchmarking.**

### Decisions — update (2026-06-15, supersedes the LightGBM/Python approach above)

6. **Rust-only trainer (NO Python).** The GBDT is trained in Rust — gradient-boosted
   trees on the logistic loss + isotonic calibration (PAVA), all in `crates/model-train`.
   The shipped binary was already pure-Rust for *inference*; this makes *training* pure-Rust
   too. Consequences: the LightGBM training script, the `model.txt`→AGBD transcoder, AND the
   cross-language feature-parity gate all **disappear** — there is now ONE feature extractor
   (`crates/scoring/src/peak_features.rs`) and ONE ion model (the engine's own), used by both
   training and inference, so the most fragile surface (Rust↔Python parity) is gone. The Rust
   trainer emits a `GbdtPeakModel` (the SoA walker type) directly — no transcode step.
7. **One unified `andes train` command (promoted from `train-from-msnet`).** It is
   **source-agnostic**: training PSMs come from many provenances — PRIDE reanalysis, MSnet,
   quantms, etc. — normalized to one harvested-PSM input contract, with provenance recorded
   per-input via `--source pride|msnet|quantms|...` (replacing the hardcoded `'msnet'` ledger
   tag at andes.rs:2829). A single `andes train` invocation fits BOTH the rank-table core AND
   the GBDT peak model in one pass (sharing the spectrum prep + feature extraction), seeded
   from the archived prior, and writes/updates `models.parquet` (rank tables +
   `gbdt_model_bytes` blob).
   - CLI collision resolution (user, 2026-06-15): andes already had two subcommands —
     `train` (spectra + FASTA, self-labeling) and `train-from-msnet` (externally-labeled gold
     PSMs). **Promote the gold-PSM trainer to `train`**; **rename the existing spectra+FASTA
     `train` → `train-from-search`.** Both kept, clear names.
8. **Rust trainer design (confirmed):** histogram gradient-boosting on the logistic loss +
   isotonic calibration (PAVA), in `crates/model-train`, emitting a `GbdtPeakModel` directly.

## Decomposition

- **Sub-project A — Stronger own models (the gate).** No PR until A4 passes.
- **Sub-project B — Public release.** Consumes A's store.

---

## Phase 0 — hygiene (quick, first)

- **Archive** the MS-GF+-derived `models.parquet` (the 18 MB seed store) to
  `internal-docs/model-archive/msgf-derived-seed.parquet` (outside the repo; kept as the
  training seed/reference, never shipped).
- **Rename** `resources/ionstat/` → `resources/`; the store path becomes
  `resources/models.parquet`. Update the loader path constant + tests. Scrub remaining
  MS-GF+-ish naming ("ionstat", etc.) from code + tool docs. ("As little resemblance to
  MS-GF+ as possible.")

## Phase A — stronger own models (the gate)

### A1 · Rank-based core (keep + retrain cleanly)
andes's per-fragment LLR `log[P(rank|signal)/P(rank|noise)]` is its natural low-res
advantage (ranks robust to scaling/noise where Sage's raw intensity + Comet's coarse bins
lose). Retrain cleanly on more Codon-harvested gold PSMs (MSFragger∪Comet→Percolator — an
internal *data-generation* tool, not model derivation). Keep the existing rank-table
estimator; the only structural change is replacing the flat noise prior (A2).

### A2 · GBDT peak signal/noise model (centerpiece)
**Purpose:** replace the current peptide-agnostic, per-rank *flat noise histogram*
(`rank_dist_table[IonType::Noise]`, consumed once in `RankScorer::new`,
`crates/scoring/src/scoring/rank_scorer.rs:78`) with a learned per-peak
`P(signal | features)` → a sharper `RawScore`.

- **Model:** LightGBM binary classifier (faster than XGBoost on 10⁸ peak rows; text dump
  transcodes cleanly to a flat Rust table).
- **Features (peptide-AGNOSTIC, set "A" only):** local-window intensity rank (top
  discriminator in the denoising literature), top-K-in-window flags, global rank/rank-frac,
  intensity/base-peak, intensity/TIC, m/z + m/z-frac-of-precursor, local peak density,
  spacing to neighbors, mass defect, isotope-partner / n-isotope, complement-in-scan (uses
  precursor mass only), AA-gap-neighbor, water/NH₃-loss partner, charge hint. These are
  computable **once per spectrum** (peptide-independent).
- **Labels:** from strict-FDR (≤0.1–0.5%) confident PSMs — signal = peak matching any
  theoretical b/y(+a/loss/immonium, generously enumerated) ion; noise = the rest.
  Mitigate label noise: exclude high-co-isolation scans; for chimeric IDs label against the
  union of confident peptides; label at the model's exact fragment tolerance.
- **Training (Codon):** negative-undersample to ~4:1; optimize PR-AUC; **isotonic
  calibration** on a run-disjoint split (non-negotiable — the output feeds a log-likelihood);
  group-split by run+peptide (no leakage); shallow trees (depth 4–7, num_leaves ≤64, high
  min_child_samples), early stopping.
- **Export → Rust inference:** transcode LightGBM `model.txt` + the isotonic map (offline,
  Python) into a **flat struct-of-arrays tree table**; evaluate with a **hand-rolled Rust
  SoA tree walker** (`crates/scoring/src/gbdt_eval.rs`, ~300 lines, **zero new native deps**
  — preserves the pure-Rust build). Sub-microsecond/peak.
- **Storage:** a backward-compatible nullable `gbdt_model_bytes` (Arrow `Binary`) column on
  the manifest rows of `resources/models.parquet` (same pattern that added `loss_class`).
- **Scoring integration (parity-safe, additive):** compute features + calibrated
  `P(signal)` **once per peak at spectrum-prep** (`ScoredSpectrum::new`) and cache it
  parallel to the intensity ranks — the GBDT **never runs in the inner candidate loop**, so
  millions of peaks is a non-issue. For a matched fragment, add `log(s/(1−s))` as a **new
  additive LLR term** alongside the existing rank/intensity table term (which is kept). This
  matches andes's parity lessons: additive terms are Percolator-safe; modifying the existing
  scoring regresses. Expected gain concentrated in low-res / wide-tolerance / chimeric — the
  regimes where andes currently trails.
- **Output:** the sharpened **`RawScore`** and the derived **`DeltaRawScore`** (top1 − next)
  are the PIN features Percolator consumes. No E-value, no rescoring.

### A3 · Clean `andes train` pipeline (Codon, documented)
One reproducible flow: harvest strict-FDR PSMs → extract per-peak features + labels →
fit rank tables (existing estimator) + train LightGBM + isotonic → transcode → inject the
blob → write `resources/models.parquet`. Documented as "how to produce a model." Per-slug
training (activation/instrument/enzyme/protocol), specializing the hard regimes (TMT,
low-res CID) where the lever is largest.

### A4 · Competitive-advantage gate (VM)
andes with the new own models must, at **1% true entrapment-FDP** on the VM:
- be **≥ the best open-source engine** (Java MS-GF+/Sage/Comet/ProSE) on PSMs on **all** of
  Astral / TMT / UPS1,
- be **strictly > the field on Astral** (the high-res strength), and
- **beat the field on LysC** (non-tryptic, multi-enzyme showcase).

Baseline (existing own models, 2026-06-15): already > field on Astral (+4.6%) and TMT
(+3.7%); UPS1 top-1 is −1.5% vs Java (chimeric +6%). So A2's targets are concrete: close
the UPS1 top-1 gap and capture the ~9% Astral headroom (the quick own-Astral train trailed
its derived counterpart). **No PR until the gate passes.**

---

## Phase B — public release

Curated branch off `main` (already ~built as `release/multi-enzyme-reanalysis`,
`664e9c09`):
- Features: multi-enzyme `--enzyme`, H5/H6/H1/decoy/M11/Sec, neutral-loss primitive.
- The strong **own-only** `resources/models.parquet` (from Phase A).
- Own-model benchmark (Astral/TMT/UPS1 + LysC) vs Java MS-GF+/Sage/Comet/ProSE, uniform
  Percolator, + the **entrapment-FDP** honesty check re-run with the own models.
- Independence docs: NOTICE/README/HERITAGE → **full independence, minimal heritage**.
- Hygiene: internal docs already in `internal-docs/`; giant fixtures gitignored/uncommitted;
  tool-docs only in the repo.
- → PR into `main`. (`main` and the dev line have unrelated histories; the curated branch
  off `main` is the reviewable vehicle. Tests reconciled — green modulo the pre-existing
  `integer_mass_scaler`.)

## Testing & validation

- A2 is validated as a **pure additive change**: A/B the GBDT term on the FDP-vs-field
  entrapment benchmark across Astral/TMT/UPS (gated on Astral, not TMT-only), milestone-
  committed on a feature branch.
- Rust: unit tests for the SoA tree walker (exact match vs a reference LightGBM prediction
  on fixtures), the store round-trip of `gbdt_model_bytes`, and the once-per-peak caching.
- Parity: the existing parity/golden suites stay green for standard peptides (the GBDT term
  is additive; when absent the path is byte-identical).

## Key risks / open items

- **A2 is genuine research** — feature design, label quality (co-isolation), calibration,
  and beating the field will iterate on Codon+VM across multiple sessions. The gate (A4)
  may take several training rounds.
- **Label noise** (real-but-unmatched ions as "noise") is the main modeling risk — mitigated
  by generous theoretical enumeration + co-isolation exclusion + robust tree settings.
- **Own-slug coverage:** the clean own-only store omits combos without an own model →
  loud WARN + `--model`. Acceptable; expand in later releases.
- **Independence-claim wording** is the user's IP call (the seed-init lineage); the spec
  writes to "full independence, minimal heritage" per the decision.
</content>
