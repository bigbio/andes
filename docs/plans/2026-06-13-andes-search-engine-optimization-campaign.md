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

## Harness & infrastructure (two-host split)

- **VM = benchmark only.** Runs andes + MSFragger (+ Sage/Comet as references) and Percolator. 3-arm scripts already exist (gitignored under `benchmark/`): `run_{astral,tmt,pxd001819}_3arm.sh` + `compare_*_3arm_percolator.sh`. Gotchas: MSFragger/Sage have no native `.raw` here (pre-convert to mzML); msgf2pin crashes → use `build_pins.py`; target-only FASTA + entrapment for FDP. Access via the VM per `reference_andes_infra_layout`.
- **Codon cluster = generate data + train models.** Phase-3 retraining is staged at `/hps/.../andes-training` (manifest + driver + array ready; needs flats + PRIDE-curated training corpus). Access via the `codon-cluster` skill. This is where "multiple models based on peak ranks" get built.
- **Datasets:** Astral (high-res), UPS1, TMT — low-res CID-TMT (e.g. PXD016999 4-engine TMT set; PXD014502 ion-trap CID-TMT). Pin exact accessions/paths on first run and record them here.

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

Pacing: each iteration is real compute (search + Percolator = minutes–tens of minutes, model training on Codon = longer). Use long/dynamic `/loop` intervals tied to benchmark completion, not tight polling.

## Starter backlog (across all three axes)

**PSMs — TMT/low-res (the lever: the model):**
- Train per-protocol models on Codon: a dedicated **low-res CID-TMT** model (the known divergence) — vary peak-rank distribution, noise-model shape, segment count, training corpus. Benchmark each variant; keep the FDP-clean winner. (Model store + `andes train` machinery already exist; loss-table serialization too.)
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
