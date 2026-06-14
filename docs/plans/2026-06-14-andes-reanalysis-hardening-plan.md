# andes reanalysis-hardening + re-benchmark plan (2026-06-14)

Strategic plan: where we are, what the audits found, the phased fix program, and the
re-benchmark gates. Companion to the finding-level detail in
[`2026-06-14-andes-gaps-audit.md`](2026-06-14-andes-gaps-audit.md) and the campaign in
[`2026-06-13-andes-search-engine-optimization-campaign.md`](2026-06-13-andes-search-engine-optimization-campaign.md).

---

## 1. Where we are

### 1a. Campaign (beat MSFragger) — WON and HARDENED
All matched (same FASTA), own/no-heritage models, **1% true entrapment-FDP**, MSFragger at
the corrected `data_type=0`:

| Dataset | Regime | margin (andes over MSFragger) | replication |
|---|---|---|---|
| **TMT PXD007683** (a05058/9/60) | low-res CID | **+11–15%** | 3 held-out files |
| **TMT PXD016999** (GTEx Fusion SPS-MS3) | low-res CID | **+12–25%** | ≥5 held-out files, **both physical instruments** |
| **UPS1 PXD001819** | low-res CID/LFQ | **+13–15%** | 3 replicates |
| **Astral** (ProteoBench) | high-res | **+8.3%** (top-1) | 1 file; chimeric re-run pending |

Robust, replicated, multi-instrument, multi-dataset, all held-out. **Loose ends:** Astral
chimeric re-run (dt=0 + topN=2) and the `max_length 40→50` lever (M2) that could widen Astral.

> The earlier "+86–130% TMT" was an artifact of MSFragger `data_type=3` (caught by the user);
> corrected to dt=0. All numbers above are post-correction.

### 1b. Fixes already shipped this session (parity-verified, UNCOMMITTED, on `feat/enzyme-support`)
- **`--enzyme` digestion flag** (default trypsin; was hardcoded with no override).
- **C2 — dedup `is_decoy`**: target + reversed-decoy can no longer merge into one mislabeled PIN row (FDR correctness).
- **CodeRabbit ×3**: chimeric residual tolerance (`param().mme`), loss-ion determinism sort, `emit_loss_ions` mz≤0 guard.
- All bit-identical on the default path (java-parity 9/9, decoy-parity, candidate-gen, scoring 91/91 green).

### 1c. Audits (4 parallel agents, ~28 findings) → see the gaps-audit doc
The verdict the user named: andes is **"almost unusable" for heterogeneous zero-config
reanalysis** because real capabilities are hardcoded or silently degrade. It is *correct* for
the campaign's narrow case (tryptic, explicit `--mods`, centroided, correctly-routed
instruments) — which is why the benchmarks are valid — but it silently mis-processes the
diversity of public data.

---

## 2. Decision — tackle it now

**Yes.** The campaign objective (beat MSFragger) is met and hardened, so this is the right
moment to pivot to the parameter-free/reanalysis goal. The blocking gaps are concrete bugs,
not speculative features. Sequence the fixes by **value × (1/risk)**, gate result-changing
phases behind a re-benchmark, and parity-test after every phase.

---

## 3. The phased program

Each phase = a milestone commit. Parity tests (`cargo test`) after every phase; the campaign
numbers (§1a) are the regression baseline. **RE-BENCHMARK gates** marked where results can change.

### Phase 0 — DONE (§1b)
enzyme digestion, C2, CodeRabbit ×3.

### Phase 1 — Correctness + enzyme completion (high value, contained)
- **C1** — `--protocol TMT`/iTRAQ inject the tag mod (+229.163 / +144.10 on K + pep-N-term) into the `aa` set, **with a double-add guard** when `--mods` already supplies it. Requires building `aa` *after* protocol resolution (or an `inject_fixed_mod` API).
- **H4** — thread `--enzyme` into `build_selection_key` (model selection) so non-tryptic data loads the matching store model; default trypsin preserves the equivalence-gate test.
- **→ RE-BENCHMARK gate 1:** (i) TMT a05058 with explicit `--protocol TMT --mods` must be **unchanged** (double-add guard works); (ii) NEW check: `--protocol TMT` *alone* (no `--mods`) now finds labeled peptides; (iii) NEW: a LysC dataset with `--enzyme lysc` now produces real IDs.

### Phase 2 — Flag-exposure quick wins (low-risk, additive; default behavior unchanged)
- Asymmetric + Da precursor tolerance (`--precursor-tol-left/right-ppm`, `--precursor-tol-da`).
- `--max-mods` (max variable mods/peptide; `NumMods=` still overrides).
- `--decoy-strategy {reverse,shuffle,none}` + detect pre-existing decoys at load (no silent doubling; target-only possible).
- Honor `--fragment-tol-*` on mzML/.raw/.d (currently no-op) or add `--feature-tol`.
- **→ No re-benchmark** (defaults unchanged); a parity smoke suffices.

### Phase 3 — Silent-loss robustness (result-affecting → re-benchmark)
- **H1** profile-mode mzML: detect `MS:1000128` → centroid (or fail loud).
- **H2** MGF multi-value CHARGE: parse first token / fall back to `None`, never drop.
- **M1** non-standard residues: treat `U`=Sec, count+log dropped spans (no whole-protein drop).
- **M2** `max_length 40→50` (MSFragger parity; campaign-relevant).
- **M3** missing-charge sweep ceiling (z6/7 charge-missing scans).
- **→ RE-BENCHMARK gate 2:** `max_length→50` can change counts on ALL datasets (hopefully up, esp. Astral) → re-run all 3 datasets matched; confirm no regression + capture gains. Profile-mode only matters if a test dataset is profile (campaign ones are centroided → unchanged).

### Phase 4 — Model-resolution + auto-detect hardening (test-guarded; equivalence gate)
- **H5** route `HCD/LowRes`→`cid_lowres` instead of high-res QExactive (+WARN).
- **H6** detect EThcD/ETciD → route to HCD (not pure ETD) until an EThcD model exists.
- **L6** `select()` reports match tier; WARN on family/last-resort fallback.
- **M7/M8** stride-sample the activation/instrument/isobaric peek passes (not first-N contiguous).
- **→ RE-BENCHMARK gate 3:** campaign datasets are low-res **CID** + high-res HCD (Astral) → their routing is unchanged; verify model_id identical for a05058/UPS1/Astral. Gains accrue on HCD-low-res/EThcD datasets (not in the campaign set).

### Phase 5 — Observability + MVP polish
- **Inferred-params startup banner** (every resolved param tagged inferred/user) — the reanalysis auditability win.
- NaN debug-assert in `write_double` (L4); env-var levers → flags (M10/M11); reconcile training charge-range defaults (M12).
- `--precursor-cal` default `off→auto` — **discuss** (global behavior change; re-benchmark if flipped).

---

## 4. Re-benchmark protocol

Re-benchmarking is a **regression-check + gain-capture**, not a full redo — most fixes are
bit-identical on the campaign datasets (tryptic, explicit-mods, centroided, correctly-routed).

- **After Phase 1 and Phase 3** (the only result-changing phases for the campaign set): re-run
  the matched suite — TMT a05058 + a PXD016999 subset, UPS1, Astral top-1 + chimeric — at 1%
  true entrapment-FDP, MSFragger dt=0. Confirm: (a) no regression vs §1a, (b) the
  `max_length→50` gain.
- **Every phase:** `cargo test` parity suites; assert default-path bit-identical where claimed.
- **New coverage to add:** a non-tryptic (LysC) mini-benchmark and a `--protocol TMT`-without-`--mods` check, so C1/H4 are regression-guarded going forward.

## 5. Sequencing, branch, effort

- **Recommended order:** Phase 1 → re-benchmark gate 1 → Phase 3 → re-benchmark gate 2 (lock
  the new baseline) → Phase 2 (anytime; additive) → Phase 4 → Phase 5. Phase 1+3 first because
  they're the correctness + result-affecting ones; lock the baseline early.
- **Branch:** the done fixes (§1b) sit in `feat/enzyme-support`'s working tree alongside parked
  glyco WIP + pre-existing failing smoke tests. Recommend a **clean `feat/parameter-free`
  branch** off the integration branch for this program (move the §1b fixes onto it), keeping
  glyco separate. One closing PR per the milestone-commit model.
- **Effort:** ~28 findings across 5 phases is a multi-day program. Phases 1–3 (the
  correctness + usability core) are the high-value chunk; 4–5 are hardening/polish.

## 6. Author decisions (2026-06-14) — RESOLVED
- **Branch:** continue on **`feat/enzyme-support`** (this branch was created to support
  multiple enzyme digestions — the audit found trypsin hardcoded on it anyway, so the
  hardening belongs here). No new branch.
- **Re-benchmark scope:** **re-bench everything** — full 3-dataset matched re-run at each
  result-changing gate (not just the regression subset).
- **`--precursor-cal`:** flip default **`off → auto`** (it's better most of the time) — moved
  into Phase 1/2 rather than deferred; re-bench picks up any change.

> **Audit principle (author, 2026-06-14):** every hardcoded value is a suspected missing
> parameter — when the audit (or any code read) finds a literal, ask *why* it's hardcoded and
> whether it should be a CLI flag / config / SDRF-driven option. Apply this lens to all future
> audits and reviews.
