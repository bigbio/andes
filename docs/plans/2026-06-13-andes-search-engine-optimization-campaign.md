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
1. **Codon — gold-standard PSMs:** run MSFragger → Percolator on the public
   training corpus (the ~29 experiment-class slugs) to produce confident
   gold-standard PSMs. Needs the flats + PRIDE-curated corpus staged at
   `/hps/.../andes-training`. Access via the **`codon-cluster` skill**. Hand the
   gold PSMs to the VM.
2. **VM — train andes on those labels:** `andes train` with the Codon gold PSMs as
   the confident-label set → an own-data `models.parquet` with **zero MS-GF+
   heritage** (incl. the per-protocol low-res CID-TMT model = backlog experiment
   #3, which feeds sweeps #6/#9/#11). This is what the baseline table benchmarks.
   Update NOTICE / independence status when 39/39 are own-data.

(A current-MS-GF+-heritage-models snapshot is optional reference only — explicitly
**not** the independence baseline; the user chose to gate on the trained table.)

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
  **MSFragger (∪ Comet for consensus) → Percolator** on big *public* datasets,
  fast, to produce high-confidence gold-standard PSMs. These are the **training
  labels** handed to the VM. (Consensus across two engines removes single-engine
  selection bias; see the methodology note.) (Independence-clean: andes's model is learned from public spectra with
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
| — | — | (loop appends rows here) | | | | | |

**Known dead ends (do not re-try):** chimeric-TMT (FDP-flat), fragmentation overlay (all 3 adjustments fail), rank-model "ceiling" via data/recal (no lever), fragment-index speed-v2 (recall/speed tension). See the project memory for details.
