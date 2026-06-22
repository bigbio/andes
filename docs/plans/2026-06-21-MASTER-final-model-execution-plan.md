# MASTER execution plan — andes final fully-own, field-beating model

**Date:** 2026-06-21
**Supersedes/integrates:** `2026-06-21-model-v2-and-training-followups.md` + `2026-06-21-model-improvement-roadmap.md` (rationale lives there; this is the HOW).
**Single source of truth for the campaign.**

---

## ⚠️ 2026-06-21 EXPERIMENT RESULT — Phase-2 levers REFUTED, roadmap redirected

Ran action #1 (F2 + A1) on Astral @1% entrapment-FDP (binary `dbfbe630`, all via `--model-store`):
- shipped/v1 (dormant GBDT): **36,799** PSMs, RichIonLLR 0% live.
- A1-only (fresh train + GBDT packed): **36,878** (+0.2% = **FLAT**), RichIonLLR 99.2% live.
- F2+A1 (dense noise + GBDT): **31,059** (−15.8% = **REGRESSION**).

**Conclusions that change the plan:**
1. **A1 (activate GBDT) is FLAT on closed search** — the discriminative layer is correctly wired and activates, but pays ~nothing on a *fixed* candidate space. Confirms the expand-then-discriminate thesis.
2. **F2 (dense noise) is a regression — DROP it.** Default noise wins; the noise-density diagnosis of the Astral regression was wrong.
3. **Retraining is safe** — a plain fresh train (default noise + GBDT) ≈ the champion (36,878). The earlier "−32% Astral regression" was a v2 corpus/config artifact, not inherent.

**Redirect:** (a) drop F2 everywhere; (b) keep A1 (ship GBDT-packed models — free, correct, own — but not as a standalone lever); (c) **PROMOTE expand-then-discriminate (semi-tryptic + chimeric + group/subset-FDR with the live RichIonLLR/strong-score) to the #1 PSM lever**; (d) final-model recipe per slug = plain fresh train (default noise + GBDT packed). Detail in memory `project_f2a1_gbdt_activation_2026_06_21.md`.

**Corrected next experiment:** on the GBDT-live `hcd_qexactive_tryp`, run **`--ntt semi` + group-FDR** on Astral and measure whether the now-live discriminative score converts the expanded candidate space into PSMs above 36,878 at ≤ FDP. That is the real test of the architecture's value.

---

## 0. Objective & definition of done

A final `resources/models.parquet` that is:
1. **Fully own** — no MS-GF+ table values AND no MS-GF+ partition geometry (independence gate passes; E3 done).
2. **Discriminatively live** — ships the GBDT blob columns (RawScore/RichIonLLR/frag-LLR battery non-zero), not the dormant pre-GBDT schema.
3. **Field-beating** — ≥ the reference engine/a comparison search engine on confident PSMs at **1% true entrapment-FDP** across Astral (high-res HCD), TMT (low-res CID), UPS1 (LFQ), and HLA.
4. **Honestly trained** — dataset-wide q-value labels, density-decoupled noise, capped corpora, per-slug own training (no silent seed fallback; gate-enforced).
5. **Shipped as ONE consolidated PR** — models + pipeline fixes + scoring + gate + NOTICE/license update.

**Hard gates throughout:** every change A/B'd one variable at a time on the VM benchmark at 1% true entrapment-FDP; experiment-hygiene provenance (binary commit + model SHA + flat SHAs); gate keep/drop on Astral, confirm no UPS/TMT regression before banking; `scripts/check_models_independence.py` must pass (zero seed-identical) before any ship.

---

## 1. Infrastructure / harness (where each thing runs)

- **Local repo** (`/Users/yperez/work/msgfplus-workspace/msgf-rust`): all Rust code changes (F2, F3, A1-packing, E3, W1, W3), the gate script, the consolidated PR. Build the release binary here, ship to Codon/VM.
- **Codon** (`/hps/nobackup/juan/pride/reanalysis/andes-training`, codon-cluster skill): harvesting (F1 python edit), training (`cluster_slug.sh`/`submit_array.sh`, fixed), `andes-bin` (rebuild from main after each code change), stores, flats, seed-models.parquet, assembler.
- **VM** (`pride-linux-vm`, `~/.ssh/andes-vm`): the benchmark (`ab_a05058.sh`, Percolator 3.7.1, `compute_entrapment_fdp.py`, `TMPDIR=/srv/data/msgf-bench/abtmp`).
- **Provenance rule:** one binary commit → rsync crates + `cargo build --release -p andes` on Codon AND VM; never let the two hosts drift (experiment-hygiene memory).

---

## Phase 1 — Pipeline foundation (unblocks honest training)

**Goal:** a training pipeline that produces clean, own, discriminatively-live models. Everything downstream depends on this.

### 1a. Bug B durable fix + harness hygiene  *(Codon, small)*
- `cluster_slug.sh`: replace `--seed-model "$SLUG"` with a **slug→valid-seed-base map** (TMT/nocleavage-phospho variants → their base in the 39-seed store, e.g. `cid_lowres_tryp_tmt→cid_lowres_tryp`, `hcd_highres_nocleavage_phosphorylation→hcd_highres_nocleavage`). Add the **exit-code failsafe**: if `train-from-msnet` exits non-zero, `echo FAILED; exit N` — never print "trained OK."
- Add `etd_highres_tryp_phosphorylation` (and any other data-having slug) to `manifest.tsv`.
- Bump training job `--mem` to 90G (fix the `cid_lowres_tryp_tmt` GBDT OOM at 498k PSMs) — superseded by F3 subsampling.
- **Gate:** re-run the 3 gap slugs; all write own stores (tables ≠ seed via `check_models_independence.py`); a deliberately-bad slug errors loudly.

### 1b. F1 — dataset-wide q-value labels in flats  *(Codon python, low effort, highest quality lever)*
- Edit `mzml_pepxml_to_flat.py`: replace the per-PSM `expect ≤ 0.01` cut with a **target/decoy TDC q-value ≤ 0.01** (the reference engine emits decoys; mirror `crates/model-train/src/labeled.rs` monotone q-walk + conservative tie-bucketing). Keep rank-1, charge≥… , peaks≥10.
- **Validation:** regenerate flats for `cid_lowres_tryp` + `hcd_qexactive_tryp`; retrain; A/B old-flats vs new-flats → confident PSMs @1% entrapment-FDP on UPS1/Astral. Expect equal-or-more PSMs at *lower* FDP.
- If neutral-or-positive (expected): regenerate ALL flats with the new cut before the final train.

### 1c. F2 — density-decoupled noise default for high-res  *(Rust, low effort, the Astral −32% lever)*
- In the training noise path (`crates/scoring/src/scoring/scored_spectrum.rs` `noise_match_facts` / `dense_noise_facts`, wired in `crates/model-train/src/accumulate.rs:110-116`): make the density-decoupled sampler the **default for high-res** slugs (or normalize noise counts per-spectrum before accumulation), instead of requiring `ANDES_DENSE_NOISE`. Keep low-res behavior unless the A/B says otherwise.
- **Validation:** retrain `hcd_qexactive_tryp` on the *large* corpus with the new default vs current → Astral PSMs@1% entrapment-FDP. **The −32% regression should shrink or invert.** This is the decisive single-variable test.

### 1d. F3 — per-slug corpus cap / subsampling  *(Rust or harness, low effort)*
- Add a `--max-train-psms` (or per-slug cap) with deterministic subsampling in the train path; default to the learning-curve knee (~100–250k, confirm below). Fixes the GBDT OOM and the over-training tail.
- **Validation:** learning-curve sweep on `hcd_qexactive_tryp_phosphorylation` (50k/100k/250k/500k/1.1M) → PSMs@1% FDP; pick the knee. Also dedup on (seq,charge,mod) and cap per-PXD contribution (harvesting #4).

### 1e. A1-packing prerequisite — ensure the trainer writes GBDT blobs  *(verify, Rust)*
- Confirm `crates/model-train/src/store/write.rs` packs `frag_intensity_model_bytes`/`rich_ion_model_bytes`/`gbdt_model_bytes`, and that `train-from-msnet` runs the GBDT fit and writes them (`gbdt/{frag_dataset,ion_dataset}.rs`). The dormant shipped store predates these columns — the *new* training MUST emit them. (The phospho repro already built a GBDT, so the path works; verify the columns land in `store.parquet`.)

### 1f. Audits  *(cheap, do now)*
- `uvpd_qexactive_tryp`: inspect source mzML activation cvParams — 295k UVPD-bottom-up PSMs is implausible; quarantine if mislabeled HCD.
- Assembler relabel: the current assembler relabels `cid_lowres_tryp`→14 cid slugs + `hcd_qexactive_tryp`→2 phospho. For the final model, train the data-having ones independently; relabel only the truly data-less.

**Phase 1 exit criteria:** clean-label flats, density-decoupled noise default verified to fix Astral, corpus cap set, GBDT blobs confirmed in new stores, Bug B durable, audits done.

---

## Phase 2 — ★ Activate the discriminative layer (A1, the big lever)

**Goal:** measure the real ceiling by turning the dormant GBDT stack live.
- Retrain the **3 benchmark models** (`hcd_qexactive_tryp`, `cid_lowres_tryp`, plus the TMT/UPS-relevant ones) with the full Phase-1 pipeline (F1+F2+F3) AND the GBDT blobs packed.
- Repack into a candidate store; **benchmark vs the shipped store** on Astral/TMT/UPS1 @1% entrapment-FDP. Confirm `RawScore`/`RichIonLLR` columns are **non-zero** in the PIN.
- **Gate:** A1 candidate ≥ shipped on all three (expect a large jump from dormant→live). This is the first true measurement of the architecture's ceiling and tells us how much Tiers 3–4 can add.

---

## Phase 3 — Near-zero-risk wins (parallel, after the A1 baseline)

A/B each on the A1 baseline; keep only what helps at flat-or-better FDP.
- **W1 — group/subset FDR — Percolator-only, RESCORING-LAYER, NOT andes code.** andes does NOT compute production FDR (Percolator owns PSM FDR; we care ONLY about PSM FDR). **FDR tool = Percolator ONLY — no Mokapot/other tool.** andes's only job = emit a group-identifying PIN column (ntt-class / charge / mod-class — mostly already present); group-FDR, if pursued, = a thin post-process of **Percolator's own output**. (`tdc.rs` is training-only — never wire it into production output.) See memory `feedback_andes_fdr_boundary`.
- **W2 — enable 1+ + `--hla` preset** (`andes.rs` defaults): `charge-min 1`, no fixed Cam-C, cysteinylation var-mod, 20/20 ppm, length 8–12/8–25; validate on HLA with NetMHCpan binders.
- **W3 — cheap additive PIN features** (`crates/output/src/pin.rs`): `ln(numCandidates)`, `log(MS2IonCurrent)`, all-ladder ppm-stdev; one-at-a-time, parity discipline; drop `IsolationWindowEfficiency` dead column.

---

## Phase 4 — Train the full data-ready set + own geometry (E3)

**Goal:** every data-ready slug trained own, on own geometry.
- **E3 — re-derive partition geometry from corpus** (`crates/model-train/src/estimate.rs` currently copies `num_segments`/layout from the seed template): derive mass-tier boundaries + charge range + `num_segments` from corpus stats; new `Param` constructor that builds structure from derived geometry (no seed template). Sweep 2 vs 3 vs 4 segments on a high-density regime (Astral) for the resolution/sparsity trade-off.
- **Train all 16+ data-ready slugs** via the fixed `submit_array.sh` with the full pipeline (F1+F2+F3+A1+E3). Re-run `cid_lowres_tryp_tmt` with the cap (no OOM).
- **Gate:** each store own (≠ seed) AND GBDT-packed AND own geometry; per-slug entrapment-FDP gate where a bench is staged.

---

## Phase 5 — Coverage expansion (SAFE only after F2 + A1)

A/B each via the per-source ledger (`train --update --add/--remove-source`, exact/reversible).
- **C1 — cross-species pooling**: add E. coli (PXD018176, 242k) + yeast flats *already on Codon* into `hcd_qexactive_tryp`, weight-swept (0.25/0.5/1.0) → Astral PSMs@1%. Fills high-charge/high-mass tails. (F2 makes this safe — without it, heterogeneous data re-triggers the noise inflation.)
- **C2 — timsTOF curation**: harvest the empty `cid_tof_*`/`hcd_tof_alp` slugs + add a **timsTOF-HLA** slug (the one data-rich gap; PRIDE has plenty). Do NOT chase ETD/UVPD bottom-up (genuinely scarce → rely on backoff from `etd_highres_tryp`).
- **C3 — semi-tryptic (`--ntt semi`) + group-FDR (W1)**: biggest free coverage expansion; must pair with grouped FDR.
- **C4 — length-normalized scoring for HLA** (`rank_scorer.rs`/`pin.rs`): `RawScore/peplen` + length prior so 9mers aren't out-competed.
- Harvest the easy TODO slugs (`hcd_highres_tryp_phosphorylation`, `cid_lowres_aspn` — abundant data).

---

## Phase 6 — Assemble final + benchmark + ONE PR

- **Assemble** best-own-per-slug (all own tables, own geometry, GBDT-packed) + relabel only the truly data-less; the v3 logic generalizes (use real own phospho/TMT now).
- **Independence gate** (`check_models_independence.py` vs seed): MUST pass (zero seed-identical).
- **Full benchmark** vs the reference engine + a comparison search engine on Astral/TMT/UPS1 (+ HLA) at 1% true entrapment-FDP, uniform Percolator, fresh wall + uniform parsimony.
- **NOTICE/LICENSE/README**: now truthfully independent (own tables + own geometry) → update attribution / relicense (the legal step that was deliberately held).
- **ONE consolidated PR**: final `resources/models.parquet` + the pipeline fixes (Bug B, F1/F2/F3, A1-packing, E3) + scoring (W1/W3) + the gate (already PR #5) + NOTICE/license. (Fold or merge the already-open #5/#6 as prerequisites.)

---

## Phase 7 — Longer tail (post-final, future PRs)

- **B1** dedicated HLA model + 1+/13–25mer/class-II corpus (Sarkizova MSV000084172 / Abelin MSV000080527).
- **B3** held-out likelihood + early-stop + seed averaging (regression detectable at train time).
- **B4** self-labeling with andes once it leads a regime (removes the the reference engine-ceiling + independence concern).
- **B5** retention-time features + a learned GBDT re-scorer over the full PIN (biggest missing feature category).
- Harvest the remaining ~20 TODO slugs as data lands.

---

## Dependency graph (critical path)

```
P1 (foundation) ──> P2 (A1 activate, measure ceiling) ──> P3 (cheap wins)
       │                                   │
       └──> P4 (E3 + full train) ──────────┴──> P5 (coverage) ──> P6 (assemble + gate + bench + ONE PR)
```
- P1 gates everything (clean labels + noise fix + GBDT packing).
- P2 must precede P3–P5 (can't A/B levers on a dormant baseline).
- E3 (P4) and the NOTICE update (P6) are the independence long-poles.
- P5 coverage is only safe after F2 (P1).

## Effort / risk summary

| Phase | Effort | Risk | Key gate |
|---|---|---|---|
| P1 foundation | med | low | Astral recovers w/ F2; GBDT cols non-zero; labels cleaner |
| P2 A1 activate | med | low | candidate ≥ shipped on all 3 |
| P3 cheap wins | low | low | each A/B ≥ flat FDP |
| P4 E3 + train | high | med | own geometry ≥ current; all stores own+GBDT |
| P5 coverage | med | med | each source A/B ≥ baseline |
| P6 assemble+PR | med | low | gate passes; beats field; NOTICE truthful |

## Immediate next 3 actions (this week)
1. **P1c (F2) + P2 (A1)** together: rebuild andes, retrain `hcd_qexactive_tryp` with density-decoupled noise + GBDT blobs, benchmark Astral. *One experiment that tests the two biggest levers.*
2. **P1a** Bug B durable `cluster_slug.sh` fix + re-run `cid_lowres_tryp_tmt` with the cap.
3. **P1b (F1)** flat label-cut edit + the `cid_lowres_tryp`/`hcd_qexactive_tryp` A/B.

## Open decisions for the user
- Merge #5 (gate) + #6 (code/heritage removal) now as clean prerequisites, or fold into the one PR? *(rec: merge now.)*
- Final model scope: ship best-own-per-slug for the 16 data-ready + relabel the rest now, then expand coverage in follow-up PRs? Or hold the final PR until timsTOF/HLA coverage lands? *(rec: ship the 16-slug field-beating final first, expand after.)*
