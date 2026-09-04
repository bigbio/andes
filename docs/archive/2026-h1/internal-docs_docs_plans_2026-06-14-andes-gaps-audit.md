# andes gaps & footguns audit (2026-06-14)

Consolidated, deduped, severity-ranked findings from a 4-agent parallel audit of the
andes search engine, motivated by the **parameter-free / reanalysis** goal (zero-config
`andes --database X --spectrum Y --output-pin Z`, robust across heterogeneous public data).
Each finding is verified against real code (file:line). Convergent findings (flagged by
more than one agent) are the strongest signals.

**Already fixed this session (not in the backlog below):**
- `--enzyme` flag for **digestion** (default trypsin) — was hardcoded with no override. Shipped + parity-verified.
- CodeRabbit: chimeric residual tolerance (`param().mme`), loss-ion determinism sort, `emit_loss_ions` mz≤0 guard.

> **Campaign benchmarks are unaffected** by everything below: TMT runs passed explicit
> `--mods` (so C1 didn't bite); a05058/UPS1 are low-res CID (not the HCD/profile/EThcD
> paths); entrapment-FDP is robust to the rare C2 collision. These are reanalysis-hardening
> items, not corrections to the andes-vs-MSFragger results.

---

## CRITICAL — correctness / catastrophic silent loss

| ID | Finding | file:line | Fix |
|---|---|---|---|
| **C1** | **`--protocol TMT`/iTRAQ does NOT inject the tag mod** (+229.163 / +144.10 on K + pep-N-term). Only the reporter peak-filter engages. So `--protocol TMT` without a `--mods` carrying the tag → every labeled peptide is precursor-mass-off → catastrophic loss. Fires on **explicit `--protocol` AND auto-detect.** | `andes.rs:1152-1171` (protocol arms set only `data_type.protocol`); `aa` built at `1007-1053` w/ CAM+Ox only | When protocol resolves to TMT/iTRAQ, inject the plex tag as fixed mod on K + pep-N-term into the `aa` builder **before** `build()`; guard against double-add if `--mods` already supplies it. Needs the aa-build to run **after** protocol detection (or an `inject_fixed_mod` API). |
| **C2** | **Dedup key omits `is_decoy`** → a target and a reversed-decoy peptide with the same sequence+rounded score **merge into one PIN row**; Label (+1/−1) decided by heap/insertion order → a decoy can be emitted as target (and carry mixed accessions). Directly corrupts Percolator TDC FDR + quantms parsimony. **(Confirmed by 2 agents.)** Dead code `peptide_has_target_match` documents the intended "all-decoy→−1" rule but is never called. | `match_engine.rs:2092-2095,2145-2164`; `pin.rs:334,521-525`; dead code `search_index.rs:45-52` | Add `is_decoy: bool` to `DedupMapKey` + its `Ord` (from the primary candidate) so cross-label merges are impossible; `debug_assert!` merged idxs share the label. One-field fix. |

## HIGH

| ID | Finding | file:line | Fix |
|---|---|---|---|
| **H1** | **Profile-mode mzML scored as centroided** — `build_peaks` never checks `MS:1000128` (profile) vs `MS:1000127` (centroid); Thermo path requests centroids, mzML has no equivalent. Profile shoulders flood the intensity-rank scorer → near-unscorable. Common for un-centroided Orbitrap/Astral conversions. | `crates/input/src/mzml.rs:~349-379` | Detect profile representation → centroid in `build_peaks`, or fail **loud** per-file (never silent). |
| **H2** | **MGF multi-value `CHARGE` (`2+ and 3+`) drops the whole spectrum** (peaks+PEPMASS valid). Non-fatal in default search (counts+continues) but **fatal with `--precursor-cal auto/on`**. Extremely common in msconvert MGF. | `crates/input/src/mgf.rs:160-167` | Parse first numeric token; on failure → `None` (let engine sweep `charge_range`), never drop. |
| **H3** | **`--fragment-tol-ppm/-da` silently discarded for mzML/.raw/.d** (dropped with a stderr WARN; high-res feature matching locked at 20 ppm). The flag no-ops on the formats that matter. | discard `andes.rs:1176-1182`; default `rank_scorer.rs:214-223` | Let `--fragment-tol-*` set the override even when instrument is detected, or add `--feature-tol`. |
| **H4** | **Enzyme hardcoded in MODEL SELECTION too** (`build_selection_key` → `enzyme: "Trypsin"`), so even with the new `--enzyme lysc`, digestion is LysC but a **trypsin model** loads — the 13 non-tryptic store models stay unreachable. Completes the enzyme fix. | `andes.rs:3233` | Thread `cli.enzyme` → `SelectionKey.enzyme` via an `Enzyme`→store-string map; default trypsin preserves the equivalence-gate test. |
| **H5** | **Low-res HCD scored with high-res QExactive model** — `("HCD","LowRes") => ("HCD","QExactive")`; ion-trap HCD then matched at 20 ppm on 0.5-Da peaks (~−18% PSMs, silent). No `hcd_lowres` model exists. | `andes.rs:3206` | Route `HCD/LowRes`→`cid_lowres_tryp` (low-res tol, b/y) + WARN; or train `hcd_lowres_tryp`. |
| **H6** | **EThcD/ETciD collapsed to pure ETD model** — c/z-ion model scores b/y-bearing spectra → systematic under-scoring of PTM/phospho data. mzML & Thermo readers both fold supplemental activation into ETD. | `mzml.rs:488-492`, `thermo.rs:259-260` | Detect EThcD (ETD+HCD/CID in one block / Thermo codes 5/6/9/10); route to HCD until an EThcD model exists; log it. |
| **H7** | **Flag-exposure class (capability exists, no CLI flag)** — same pattern as enzyme: **(a)** asymmetric + **(b)** Da precursor tolerance (`PrecursorTolerance::asymmetric`, `Tolerance::Da` exist; CLI forces symmetric ppm); **(c)** `--max-mods` (max variable mods/peptide, hardcoded 3); **(d)** `--decoy-strategy` (reversal-only, **no target-only mode, and a pre-decoyed FASTA gets doubled**). | tol `tolerance.rs:44`+`andes.rs:1188`; max-mods `search_params.rs:99`; decoy `decoy.rs`+`andes.rs:994` | Add the four flags; for decoy, detect pre-existing decoys at load + allow `none`. |

## MEDIUM

| ID | Finding | file:line | Fix |
|---|---|---|---|
| **M1** | Non-standard residues (`U` Sec, `O`, `X`, `B/Z/J`) silently drop the **entire candidate span** → selenoproteins (GPX, SELENOP, TXNRD) + proteogenomics/metaproteomics DBs 100% lost. No counter. | `candidate_gen.rs:221` | Treat `U`=Sec (150.95364); count+log dropped spans. |
| **M2** | **`max_length=40` vs MSFragger/Comet ~50** → long tryptic peptides dropped. **Campaign-relevant** (Astral). | `search_params.rs:99`, `candidate_gen.rs:151-252` | Raise default to 50; re-run Astral A/B per milestone discipline. |
| **M3** | Charge-missing precursors only swept `charge_min..=max` (2..5) → real 6+/7+ charge-missing scans get 0 PSMs. | `match_engine.rs:317-320` | Raise/adapt the missing-charge sweep ceiling. |
| **M4** | Explicit precursor charge **outside** `--charge-min/max` → **all-zero charge one-hot** (degenerate Percolator feature). | `pin.rs:174-176,406-409` | Clamp/skip/widen header, or `chargeN+` overflow column; warn. |
| **M5** | Charge-missing scan can emit **two PIN rows** (z2 & z3) — per-charge dedup, no cross-charge pass after merge. | `match_engine.rs:595-621` | Re-run `dedup_pepseq_score` on the merged spectrum queue. |
| **M6** | `is_decoy = starts_with("XXX")` greedy prefix (a target accession starting `XXX` mislabeled) + decoys generated **unconditionally** (pre-decoyed FASTA → 1:1 broken → FDR mis-estimate). | `candidate_gen.rs:49`, `decoy.rs:11-27` | Token-boundary prefix; detect pre-existing decoys → skip/err. |
| **M7** | Activation/instrument detected from **first 64 scans only** → decision-tree/alternating + lead-in files mis-routed for the whole run. | `andes.rs:3066-3083`, `mzml.rs:994,1157` | Stride-sample across the file (like the isobaric path). |
| **M8** | Isobaric **false-negative** on low labeling + first-1000-contiguous sampling; SPS-MS3 reporters (MS3) unseen at MS2. | `isobaric.rs`, `andes.rs:889` | Stride-sample; add a "weak TMT" log tier; prefer SDRF when available. |
| **M9** | Chimeric gates hardcoded (`max_kl=0.3`, `max_n=2` secondaries/scan). | `match_engine.rs:881-882` | Expose `--chimeric-max-coisolated` / `--chimeric-max-kl`. |
| **M10** | Isobaric peak-filter window (100 Da / K=20) — a documented **~+3.5% PSM lever** — is **env-var-only** (`ANDES_PEAK_WINDOW`). | `scored_spectrum.rs:288,201-208` | Surface `--peak-window-da` / `--peaks-per-window`. |
| **M11** | Precursor-cal skip thresholds (`MIN_SPECKEYS=10_000`, `MIN_CONFIDENT_PSMS=150`) unexposed → small/targeted runs **never calibrate**, silently. | `precursor_cal.rs:76-94` | Expose or scale by dataset size. |
| **M12** | Training charge range hardcoded `2..=5`, ignores `--charge-min/max`; `default_tryptic` doc says `2..=3` — **3 divergent defaults**. | `andes.rs:1805`, `search_params.rs:103`, `andes.rs:142-147` | Honor flags in `build_train_search_params`; reconcile the doc. |
| **M13** | Precursor calibration skipped for **`.raw`/`.d`** (andes's native-format differentiator) — loses the tightened-tolerance gain on exactly those inputs. | `andes.rs` (`is_mzml`-gated) | Add a native-format metadata pre-pass for calibration. |

## LOW

| ID | Finding | file:line |
|---|---|---|
| **L1** | `--fragmentation` lacks PQD; metadata-less MGF can't request HighRes/TOF/TimsTOF/Astral (only QExactive/LowRes). | `andes.rs:55-61,3016-3022` |
| **L2** | Deconvolution fragment-charge cap hardcoded to 3 (z≥5 precursors lose 4+ charged fragments). | `scored_spectrum.rs:1769` |
| **L3** | `--score strong` retention pool `K=25` is a `const` (doc says "raise if it keeps recovering PSMs"). | `search_params.rs:20` |
| **L4** | `write_double` maps NaN/inf **and** exact-zero to `0` (Percolator-safe but masks feature-NaN bugs). | `pin.rs:539-542` |
| **L5** | Multi-file `ScanNr` collisions (SpecId unique, ScanNr restarts per file) → downstream join-by-scan (quantms/SDRF) mis-assigns. | `pin.rs:391` |
| **L6** | Model-store `select()` last-resort fallback (`hcd_qexactive_tryp`) is indistinguishable from a real match in logs. | `andes.rs:3277-3287`, `select.rs:229` |
| **L7** | `--mods` **replaces** the whole mod set (add phospho → silently lose Cam-C/Ox-M); defaults duplicated in 2 places (drift risk). | `andes.rs:1024-1047,2922-2945` |

---

## Implementation batches (proposed)

- **Batch A — Correctness (fix first, with tests):** C1 (TMT/iTRAQ tag-mod injection), C2 (dedup `is_decoy`). Real bugs; C2 is a clean one-field fix.
- **Batch B — Flag-exposure quick wins (additive, low-risk):** H4 (enzyme→model-selection), H7a-d (asymmetric/Da precursor tol, `--max-mods`, `--decoy-strategy` incl. target-only + no-double), H3 (honor `--fragment-tol-*`).
- **Batch C — Silent-loss robustness (reanalysis):** H1 (profile mzML), H2 (MGF multi-charge), M1 (non-standard residues), M2 (max_length→50, re-A/B Astral), M3 (missing-charge sweep).
- **Batch D — Model-resolution + auto-detect hardening (test-guarded; equivalence-gate):** H5 (low-res HCD routing), H6 (EThcD), L6 (match-tier WARN), M7/M8 (stride-sampling).
- **Batch E — Observability + MVP polish:** inferred-params startup banner (reanalysis auditability), L4 (NaN debug-assert), M9/M10/M11 (env-var levers → flags).

**The cheapest universal win** (the auto-detect agent's rec): make `select()` report its match tier and **WARN on family/last-resort fallback** + stride-sample the peek passes — makes every remaining mis-route *visible* in batch logs even before models are added.
