# andes Search-Engine Optimization Campaign (ready-to-`/loop`)

> **How to run this:** start a fresh session and launch the loop, e.g.
> `/loop drive the andes optimization campaign in docs/plans/2026-06-13-andes-search-engine-optimization-campaign.md — one FDP-gated experiment per iteration`.
> This doc is the loop's standing spec: objective, harness, protocol, guardrails, backlog.
> Keep glyco-ID in its own session (see the glyco research plan); this campaign is search-engine only.

## Goal (honest, per-dataset)

Beat **MSFragger** where there's room, hold the lead where andes already has it — all three axes (speed, PSMs, fewer params) optimized in one loop, **gated on entrapment-FDP**.

| Dataset | Current standing (from prior benchmarks) | Campaign target |
|---|---|---|
| **Astral** (high-res) | andes already **leads** PSMs/peptides + ~30× Java speed | **Hold / extend** |
| **UPS1** | andes **leads** PSMs + proteins | **Hold** |
| **TMT** (low-res CID) | andes **mid-pack / 2nd** — the real gap (low-res scoring divergence) | **Win** ← campaign heart |
| **Speed** | fast, but MSFragger fragment-index is the king | **Stay competitive** (don't regress; opportunistic wins) |
| **Params** | `--instrument` removed; mzML/.raw/.d zero-config | **Reduce further** (auto precursor-tol, enzyme, protocol) |

The realistic story: **win TMT, hold Astral/UPS, keep speed competitive, shed params.** Don't over-promise "beat everywhere."

> **⚠ Independence caveat (reframes the "standing" column above).** Those prior
> "andes leads" numbers used the **MS-GF+-derived** models — 38 of 39 bundled
> models in `resources/ionstat/models.parquet` are still MS-GF+ heritage (NOTICE:
> *scoring code independent, models in transition*). The campaign's North Star is
> **beat MSFragger with models andes trained itself, on public/quantms data —
> zero MS-GF+ heritage.** Stripped of that heritage and retrained from scratch,
> andes most likely **starts ~10% behind** MSFragger (the user's expected baseline).
> So the headline below is the gap to close, not the current borrowed lead.
>
> **⇒ The per-dataset "standing" + "target" columns above are PROVISIONAL** — they
> come from inherited-model benchmarks and may not survive the switch to own-data
> models (andes could lose the Astral/UPS lead too, not only trail on TMT).
> **Experiment #0 re-establishes the real per-dataset standing; the actual targets
> are set from that table, not assumed here.**
>
> **Methodology note (NOT circular — clarified).** Training andes on
> MSFragger→Percolator gold PSMs then benchmarking vs MSFragger is **not**
> circular: andes's scoring is fundamentally different — it learns the
> intensity-rank **fragmentation model** from the *(spectrum → true peptide)*
> pairs, which is engine-independent physics (any true-PSM source teaches the
> same; this is how MS-GF+ and all ML rescorers are trained). andes applies that
> model to *all* spectra and can exceed MSFragger on ones it scores better — the
> labels do not cap it. The only residual is minor **training-set selection bias**
> (labels drawn from the regime one engine is confident on); remove it by sourcing
> gold PSMs from **MSFragger ∪ Comet** (consensus, per the user) and report
> **entrapment-FDP** (engine-independent) as the arbiter alongside the head-to-head.

## Experiment #0 — the independence baseline table (run FIRST, it's the scoreboard)

Everything else in this campaign is measured against this table. It is also the
**independence milestone** (it requires andes on its own models), so it is the
gating prerequisite, not an afterthought.

**Prerequisite — own-data models. This GATES the baseline table** (chosen
approach: independence-valid, train-from-public-data first — not a current-models
snapshot). Two stages, Codon → VM:
1. **Codon — gold-standard PSMs:** run **MSFragger ∪ Comet → Percolator
   (consensus)** on the public training corpus (the ~29 experiment-class slugs)
   to produce confident gold-standard PSMs — consensus across the two engines
   removes single-engine selection bias (user-confirmed). Needs the flats +
   PRIDE-curated corpus staged at
   `/hps/.../andes-training`. Access via the **`codon-cluster` skill**. Hand the
   gold PSMs to the VM.
2. **VM — train andes on those labels:** `andes train` with the Codon gold PSMs as
   the confident-label set → an own-data `models.parquet` with **zero MS-GF+
   heritage** (incl. the per-protocol low-res CID-TMT model = backlog experiment
   #3, which feeds sweeps #6/#9/#11). This is what the baseline table benchmarks.
   Update NOTICE / independence status when 39/39 are own-data.

(A current-MS-GF+-heritage-models snapshot is optional reference only — explicitly
**not** the independence baseline; the user chose to gate on the trained table.)

**Training-set cleanliness — signal vs noise have different needs (open Q answered).**
andes scores `LLR = ln(ion_freq[rank] / (noise_freq[rank]·norm))`. The two halves
are learned differently and tolerate dirt differently:
- **Signal (`ion_freq`) is FRAGILE to false PSMs.** A wrong training PSM contributes
  *coincidental* matches as if they were real fragments → flattens/poisons the
  signal distribution → worse discrimination. ⇒ the gold set must be **stringent**
  (≤1% FDR; the **MSFragger ∪ Comet** members each already pass an engine's FDR
  gate, so the union is ~1%-clean *and* large — best clean-and-big balance.
  Intersection is purer but small and re-adds a selection bias).
- **Noise (`noise_freq`) is ROBUST to label dirt.** It is sampled from
  **decoy/background positions** (reversed-peptide decoy ions — `noise_match_facts`
  — or dense random positions — `dense_noise_facts`), which are background
  *regardless of whether the PSM label is right*. So the noise model does **not**
  need an ultra-clean set. Its real limiters are: **(1) corpus SIZE** (too few
  PSMs → noise under-sampled → Laplace smoothing flattens it → the diagnosed
  −4.3% "dilution"), and **(2) over-smoothing** (`noise_pseudo` must stay small so
  the noise stays sharply peaked at the missing-ion slot — guarded by
  `noise_rank_dist_stays_sharp_not_flattened_by_smoothing`).
- **Practical answer:** ~1% FDR consensus gold PSMs are **plenty clean for the
  noise**; what the noise needs is a **large** corpus + minimal smoothing. The
  signal is what demands the stringency. The exact FDR (0.1 / 1 / 5%), union-vs-
  intersection, and noise-smoothing are an **empirical sweep** in the loop (backlog
  noise-shape experiments) — watch the trained noise distribution's peak as the
  canary: if it flattens, the corpus is too small or smoothing too aggressive.

**Protocol (VM).** Matched target+decoy FASTA + foreign-proteome **entrapment**,
matched **1% FDR** (+ glycan-free here), Percolator (grep the mode — Concatenated
vs Separate counts aren't comparable), uniform peptide→protein parsimony
(`protgroups.py`) to avoid the counting-artifact class. MSFragger via FragPipe;
andes native. Same datasets, same DB, same FDR.

**The table the loop fills (per dataset: PSMs / peptides / proteins @1% FDR,
entrapment-FDP, wall-clock):**

| Dataset | MSFragger | andes (own-data models) | Δ% (andes vs MSFragger) | FDP andes / MSFragger | andes speed vs MSFragger |
|---|---|---|---|---|---|
| Astral (high-res) | — | — | — | — | — |
| UPS1 | — | — | — | — | — |
| TMT (low-res CID) | — | — | — | — | — |

WIN = andes ≥ MSFragger on IDs at matched 1% FDR with entrapment-FDP in
tolerance, on **own-data models**. The campaign drives every Δ% from negative
(the ~10% honest start) toward ≥ 0, axis by axis.

### Launch sequence (three hosts)
0a. **Codon:** MSFragger → Percolator on the public corpus → gold-standard PSMs
    (training labels). *(Long pole; gates everything.)* Hand off to the VM.
0b. **VM:** `andes train` on the Codon gold PSMs → own-data `models.parquet`
    (zero MS-GF+ heritage), incl. the low-res CID-TMT model.
1. **VM:** run experiment #0 → fill the baseline table (andes own-models vs
   MSFragger, 3 datasets, matched 1% FDR + entrapment-FDP). Record Δ% per dataset.
2. **Loop (local code ↔ VM benchmark):** round-robin the ranked backlog
   (`2026-06-13-andes-opt-experiment-backlog.md`), one atomic FDP-gated experiment
   per iteration — implement the change **on this machine**, benchmark on the VM,
   re-run the affected baseline row so Δ% stays current. Astral is the mandatory
   high-res canary on every change. (Model-variant experiments loop back through
   0a/0b.)
3. Stop when Δ% ≥ 0 on all three (or the backlog is exhausted); log refutations.

Note: backlog #1 (auto-detect isobaric protocol → the +3.5% TMT win) is a
**peak-filtering** change, independent of the models, so it applies on the
own-data models too — it is not borrowed from MS-GF+ heritage.

## Objective function

Primary, in priority order (from the andes objective: own data / no patent / beat all tools on PSMs / max speed):
1. **PSMs (and peptides/proteins) at 1% FDR — FDP-validated** via entrapment. This is the score. Raw counts without FDP are meaningless.
2. **Wall-clock speed** (and peak RAM) at equal results.
3. **Parameter count / zero-config** (fewer required flags = better), never at a PSM cost.

A change is an **improvement only if** it raises (1) with FDP held ≤ target, OR improves (2)/(3) at byte-identical (1). Everything else is a refutation — log it and revert.

## Model experimentation inside andes (standing design principle)

New scoring models are **tested in-engine** — A/B'd inside andes against the
baseline table — not bolted on externally. "Multiple models based on peak ranks"
generalizes to *any* new model family (retrained rank tables, alternative
rank/noise parameterizations, and later small learned models), under **three
standing requirements**:

1. **Native Rust inference.** A model's *inference path in andes* must be pure
   native Rust — no external runtime (Python/ONNX/JVM), no vendor FFI. This keeps
   the single static binary, the MS-GF+-independence, and "deployable anywhere"
   intact, and keeps the model fully owned. Training may happen anywhere (the VM,
   on Codon gold PSMs); only **inference** must be native Rust. Today that's the
   parquet rank tables + the GF scorer; a future learned model ships its weights +
   a hand-written Rust forward pass, loaded through the same model store.
2. **Percolator-complementary.** A model's output feeds the existing
   PIN → Percolator rescoring step — either as the primary ranking score
   (RawScore) or as **additive feature columns** Percolator weights. New models
   are "next steps into Percolator"; they do **not** replace the rescoring / FDR
   layer. (Additive features must be top-1-preserving and `--`-gated so the
   standard-3 PIN stays byte-identical — parity-tuning lessons.)
3. **Fast enough.** Inference must not regress the speed axis: native Rust + a
   bounded per-spectrum cost (the GF DP is the budget). A model that lifts PSMs
   but tanks throughput **fails** the objective (speed is axis #2) — measure wall
   on every model swap.

This is the rule the model-lever backlog items (#3/#6/#11 + any future learned
model) run under: **native Rust, Percolator-complementary, fast.**

## Hard guardrails (learned the expensive way)

The prior record is littered with plausible ideas that **entrapment-FDP refuted**: chimeric-TMT, the fragmentation overlay, the rank-model ceiling, the speed-v2 fragment index. So:

1. **FDP-gate every experiment.** Target-decoy counts alone lie. Use entrapment FASTA + report FDP. A "gain" with FDP drift is a regression.
2. **Off ⇒ byte-identical.** Any new scoring/feature behind a flag must be bit-identical to baseline when disabled (`if mode==OFF return input_unchanged` at the top — not a deep flag-branch that reorders float ops). Verify on the standard-3 goldens.
3. **One atomic hypothesis per iteration.** "All three axes" = the loop *ranges over* speed/PSMs/params, but each iteration changes ONE thing and benchmarks it. No bundled changes.
4. **No raw-count chasing.** If a change adds PSMs but FDP creeps, it's rejected, full stop.
5. **Milestone commits on measurement.** Gate "improvement" commits on an Astral measurement at minimum (TMT-only is not sufficient — high-res is the canary).
6. **Percolator mode parity.** Grep the Percolator mode (Concatenated vs Separate) before comparing counts — cross-mode counts aren't comparable.
7. **Record refutations** in this doc's log so the loop never re-tries a dead end.

## Harness & infrastructure (THREE-host split — authoritative)

The loop spans three hosts, each with a distinct job:

- **Codon cluster = gold-standard PSM generation (training labels).** Run
  **MSFragger ∪ Comet → Percolator (consensus)** on big *public* datasets, fast,
  to produce high-confidence gold-standard PSMs. These are the **training labels**
  handed to the VM. Consensus across the two engines removes single-engine
  selection bias (user-confirmed); see the methodology note. (Independence-clean: andes's model is learned from public spectra with
  externally-validated labels — no MS-GF+ parameters; MSFragger only supplies the
  label set, not any andes code/model.) Access via the `codon-cluster` skill;
  staging at `/hps/.../andes-training`.
- **VM = training + benchmarking.** (a) **Train** andes models with `andes train`
  using the Codon gold-standard PSMs as the confident-label set → an own-data
  `models.parquet` with zero MS-GF+ heritage. (b) **Benchmark** andes vs MSFragger
  (+ Sage/Comet refs) with Percolator + entrapment-FDP. 3-arm scripts exist
  (gitignored under `benchmark/`): `run_{astral,tmt,pxd001819}_3arm.sh` +
  `compare_*_3arm_percolator.sh`. Gotchas: pre-convert `.raw`→mzML for
  MSFragger/Sage; msgf2pin crashes → `build_pins.py`; target-only FASTA +
  entrapment for FDP; grep the Percolator mode before comparing counts.
- **This machine (local) = code improvements.** After a VM benchmark yields a
  conclusion, the actual andes **source changes** are made and committed here,
  then pushed for the VM to re-benchmark. (Where the `/loop` that edits code runs;
  it dispatches compute to Codon/VM and reads results back.)
- **Datasets (the actual staged VM benchmark set, per `reference_andes_infra_layout`):**
  `astral-data/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML` (high-res),
  `data/UPS1_5000amol_R1.mzML` (low-res CID LFQ), `tmt-data/a05058.mzML` (low-res
  CID-TMT) — each with its FASTA + numeric mods under `/srv/data/msgf-bench/`.
  (Held-out: `qe-holdout/` PXD048603 HeLa.) Use these, not ad-hoc PRIDE pulls.

## Loop protocol (per iteration)

```
1. PICK   the next experiment from the Backlog (highest expected-value, not yet refuted).
2. STATE  the hypothesis + the metric it should move + the FDP gate it must hold.
3. BUILD  the change behind a flag (off ⇒ byte-identical); unit + golden tests green.
4. BENCH  on the VM: 3-arm + Percolator + entrapment-FDP on the relevant dataset(s);
          ALWAYS include Astral as the high-res canary even for a TMT-aimed change.
5. JUDGE  keep iff (PSMs up & FDP held) OR (speed/params better at byte-identical PSMs).
6. RECORD result (kept/refuted + numbers) in the Experiment Log below.
7. COMMIT a milestone if kept; revert if refuted. Then loop.
```

Pacing: each iteration is real compute (search + Percolator = minutes–tens of minutes; gold-PSM generation on Codon + andes training on the VM = longer). Use long/dynamic `/loop` intervals tied to job completion, not tight polling.

## Starter backlog (across all three axes)

**PSMs — TMT/low-res (the lever: the model):**
- Train per-protocol models **on the VM from Codon-generated gold PSMs**: a dedicated **low-res CID-TMT** model (the known divergence) — vary peak-rank distribution, noise-model shape, segment count, training corpus. Benchmark each variant; keep the FDP-clean winner. (Model store + `andes train` machinery already exist; loss-table serialization too.)
- Additive Percolator features (must be *additive* — modifying existing features regresses, per parity-tuning lessons): e.g. a per-ion CID node-score trace, DeltaRawScore, complementary-ion features. Top-1-preserving only.
- Revisit the noise/rank prior shape for low-res (the dilution diagnosis: the noise rank model was over-smoothed/2.5× too flat).

**Speed (stay competitive; opportunistic):**
- Profile candidate generation + the GF DP hot path; SIMD/layout wins that are byte-identical. (Do NOT revisit the fragment-index candidate generator for chimeric speed — built + refuted: irreducible recall/speed tension.)
- Parallelism / IO chunking tuning; peak-RAM reduction.

**Params (zero-config):**
- Auto precursor tolerance (calibrate from high-confidence PSMs), auto enzyme/protocol detection from the data — extend the auto-resolution work. Each must be byte-identical when the user passes the flag explicitly.

## Experiment Log

| # | Axis | Hypothesis | Dataset(s) | PSM Δ (FDP) | Speed Δ | Verdict | Commit |
|---|---|---|---|---|---|---|---|
| 1 | PARAMS+PSMs | auto-detect TMT/iTRAQ from MS2 reporters → engage isobaric filter zero-config | TMT (a05058) | *await VM* | — | **impl done, byte-identical verified locally** (golden parity + real high-res-glyco negative control no-false-positive); FDP-gain pending VM | `b49228b6` |
| 2 | SPEED | single-segment fast path | — | — | — | **REJECTED** — verified all 39 models `num_segments==2`; never triggers | — |
| 10 | PSMs | additive `ComplementaryIonBalance` PIN feature (rank-agreement-weighted complementary b/y) | TMT | *await VM* | — | **impl done, additive byte-identical** (52 existing cols + RawScore + 633 rows unchanged, golden regen, unit test); Percolator gain pending VM A/B | `bf2d1e4e` |

### VM benchmark — a05058 TMT (latest binary #1+#10, entrapment FASTA, --protocol TMT filter, Percolator @1%)
| Model | PSMs@1% | peptides@1% | note |
|---|---|---|---|
| inherited `cid_lowres_tryp` (MS-GF+ heritage) | 11,620 (top1) / 11,678 (chim) | 10,398 | baseline (Jun-8 run) |
| own `cid_lowres_tryp` (7 flats) + #10 | 11,298 | 10,096 | −2.8% vs inherited (tiny corpus) |
| **own `cid_lowres_tryp_tmt` (72 flats) + #10** | **11,539** | **10,333** | own TMT model; q-val@1% (true FDP 2.32% — see below) |
| MSFragger (matched entrapment FASTA, rev_ decoys) | **6,791** | — | q-val@1% (true FDP 1.77%); same 52,820 targets as andes |

### HEAD-TO-HEAD — matched FASTA + 1% TRUE entrapment-FDP (mode-independent, the defensible metric)

> **Methodology fix (2026-06-14):** andes pins percolate as **Separate** mode, MSFragger
> pins as **Concatenated/TDC** — Percolator auto-detects from PIN structure, so the raw
> PSM@1%-*q-value* counts are NOT directly comparable (cross-mode). The fair metric is
> **PSMs at 1% true entrapment-FDP** (entrapment is 1:1 → FDP = 2×ENT-hits/total), which is
> independent of how each engine's q-values were estimated. Numbers below use that.

> **⚠ CORRECTION (2026-06-14): MSFragger TMT numbers were inflated by a config bug.** The TMT
> params (`fragger-tmt2.params`) had `data_type=3` (DDA+/wide-window) instead of `data_type=0`
> (standard narrow DDA). On a05058, a controlled A/B (only `data_type` differs) showed
> **MSFragger 4,965 (dt=3) → 9,369 (dt=0)** — `data_type=3` nearly *halved* MSFragger's IDs.
> So the original "+86% TMT" was an artifact. **Corrected a05058: andes 10,520 vs MSFragger
> 9,369 = +12.3%** (in line with UPS1). UPS1 and Astral-top-1 used `data_type=0` already →
> valid. Astral-chimeric used `data_type=3` → being re-run with `data_type=0 + topN=2`. TMT
> replicates (a05059/a05060) + PXD016999 being re-run with `data_type=0`. Table below corrected.

| Dataset | Regime | andes (own/no-heritage) | MSFragger 4.2 (matched, **dt=0**) | andes margin | calibration (true FDP @ q≤1%) |
|---|---|---|---|---|---|
| **TMT** a05058 | low-res CID | **10,520** | 9,369 | **+12.3%** | andes 2.32% / MF ~1% |
| **UPS1** PXD001819 | low-res CID | **16,132** | 14,282 | **+13%** | andes 1.26% / MF 1.43% |
| **Astral** (LFQ 15min 50ng) | **high-res** | **30,077** | 27,760 | **+8.3%** | andes 1.14% / MF 1.56% |
| **Astral — chimeric** | high-res | **56,243** | 52,073 (dt=3) / 27,722 (dt=0,topN2) | **~+8% vs MF chimeric** | — |

> **Chimeric nuance (2026-06-14, hardened binary):** andes `--chimeric` = 56,243 (its two-pass
> strips the primary's peaks and re-searches the residual → ~2× its top-1 of 30,259). MSFragger's
> chimeric number is **mode-dependent**: its wide-window/chimeric mode `data_type=3` = 52,073
> (andes **+8%**, the fair chimeric-capability comparison), while `data_type=0 + output_report_topN=2`
> (correct narrow-DDA primary mode) = 27,722 — i.e. narrow-mode MSFragger barely recovers
> co-isolated peptides (andes +103% there, but that's not MF's chimeric mode). Honest takeaway:
> andes's co-isolation deconvolution modestly beats MSFragger's chimeric mode and dramatically
> beats its narrow-DDA mode.

**Re-bench with the hardened binary (2026-06-14)** — after Phase 1–3 (enzyme digestion+model-selection,
C1 TMT-mod inject, C2 dedup-is_decoy, CodeRabbit ×3, max_length 40→50, --max-mods, cal=auto default,
MGF multi-charge), at the fair `max_length 50` vs unchanged MSFragger (dt=0): **andes still wins all
three, regression-free** — TMT **10,789** vs 9,369 (+15.2%), UPS1 **16,160** vs 14,282 (+13.1%),
Astral top-1 **30,259** vs 27,760 (+9.0%). Notes: cal-auto **self-gated** (skipped — ~220 confident
PSMs < threshold, so no effect; the small gains are C2 + max_length); `max_length→50` is regime-split
(slightly −1.5% on low-res TMT from FDP pressure, +0.6% on high-res Astral) but is the fair MSFragger
match. The hardening's value is **usability + correctness for zero-config reanalysis**, not raw PSM gain.

### ✅ WIN CONDITION MET (2026-06-14) — andes beats MSFragger on ALL THREE datasets

At matched 1% true entrapment-FDP, with **own-data-trained models only** (no MS-GF+
heritage in the scoring path beyond the train-from-msnet seed), andes wins all three —
**including high-res Astral, MSFragger's strongest regime**. On the two high-res-grade sets
andes is also **better-calibrated** than MSFragger. Models used: TMT `cid_lowres_tryp_tmt`
(72 flats), UPS1 `cid_lowres_tryp`, Astral `hcd_qexactive_tryp_astral` (481,039 own PSMs:
broad high-res HCD + PXD061135). Astral model carried a 0.5 Da seed fragment tol but andes
auto-resolves the high-res matching tol from the data, so it didn't hurt.

**Caveats (for full rigor, not blockers to the directional win):** (1) single representative
file per dataset, not full replicate sets; (2) matched at Cam+Ox mods (andes's optional
Acetyl-N-term dropped to match MSFragger → conservative for andes); (3) fragment tol andes-auto
vs MSFragger 20 ppm (each engine's natural high-res config); (4) andes TMT q-value calibration
is optimistic (2.32%) — a model/score-calibration lever still worth tightening.

**Per the standing rule, this unparks the glyco track** (was gated on beating MSFragger on all
3) — pending user confirmation that the comparison is rigorous enough to call.

### Replicate hardening (2026-06-14) — streaming download→convert→bench

To address the single-file caveat: stream additional replicate runs (download `.raw` from
PRIDE → convert via docker `ThermoRawFileParser` → bench both engines → delete). No mzML on
PRIDE for these (raw-only); no converter was on the VM, so the docker converter is the new
(validated) dependency.

**UPS1 PXD001819 (PSMs @ 1% true entrapment-FDP):**

| Replicate | andes (own `cid_lowres_tryp`) | MSFragger 4.2 | andes margin |
|---|---|---|---|
| R1 (original) | 16,132 | 14,282 | +13.0% |
| R2 | 16,503 | 14,308 | +15.3% |
| R3 | 16,563 | 14,603 | +13.4% |

andes wins **all three replicates** by ~+14% with **<1% variance** on its own counts — the
single-file UPS1 result is robust, not a fluke.

**TMT held-out PXD007683 (PSMs @ 1% true entrapment-FDP):**

**CORRECTED to `data_type=0`** (the dt=3 values, struck, were the config-bug artifacts):

| Replicate | andes (own `cid_lowres_tryp_tmt`) | MSFragger (dt=0, correct) | andes margin | ~~MSFragger (dt=3, WRONG)~~ |
|---|---|---|---|---|
| a05058 | 10,520 | **9,369** | **+12.3%** | ~~5,644~~ |
| a05059 | 11,054 | **9,920** | **+11.4%** | ~~5,702~~ |
| a05060 | 10,656 | **9,244** | **+15.3%** | ~~4,643~~ |

With the corrected `data_type=0`, the TMT margin is a consistent **+11–15%** across all 3
held-out files (the corrected MSFragger counts are ~2× the dt=3 artifacts) — in line with
UPS1's +13–15%. PXD007683 is **not** in the `cid_lowres_tryp_tmt` training set (clean held-out).

**TMT held-out PXD016999 — multi-instrument check (streaming):** GTEx Human Body Map (Jiang
2020), Orbitrap Fusion SPS-MS3 TMT10 (ion-trap CID MS2 for ID). 4 files × each of the two
physical Fusion instruments (`Instrument1_*`, `SecondInstrument_*`), all from **samples NOT
in the 30-file training subset** → genuinely held-out. Verifies the win holds across
instruments + a second TMT dataset. (Results pending.)

**Pipeline note:** no mzML on PRIDE for these (raw-only), no converter on VM → docker
`ThermoRawFileParser` converts each `.raw`→mzML inline. MSFragger writes a ~300 MB internal
`.mzBIN_calibrated` per run (even with `write_calibrated_mzml=0`) that must be cleaned to
avoid filling the 100 GB disk.

**Two gotchas that cost time (in memory + howto doc):** (1) MSFragger only recognises
`rev_` decoys, not `XXX_` → `docs/plans/2026-06-14-msfragger-benchmark-howto.md`. (2) andes
vs MSFragger percolate in different auto-detected modes → always compare at true
entrapment-FDP, never raw q-value counts. Also: andes's q-value calibration is **optimistic
on TMT** (q≤1% → 2.32% true FDP) but **well-calibrated on UPS1** (q≤1% → 1.26%) — a model/
score-calibration lever worth a future experiment.

**Known dead ends (do not re-try):** chimeric-TMT (FDP-flat), fragmentation overlay (all 3 adjustments fail), rank-model "ceiling" via data/recal (no lever), fragment-index speed-v2 (recall/speed tension). See the project memory for details.
