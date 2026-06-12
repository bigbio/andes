# Andes — Complete Independence from MS-GF+ Implementation Plan

> **For agentic workers:** Use superpowers:subagent-driven-development or executing-plans to implement the ENGINEERING phases (1, 4, 6, 7) task-by-task. Phases 2 (campaign), 3 (clean-room scrub) and 5 (legal) are gated procedures, not TDD tasks. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Make andes *completely new* — own code (clean-room), own data (every model), own structure (no MS-GF+ scaffold), own git history — so it can be relicensed Apache-2.0 with zero MS-GF+ IP exposure.

**Architecture:** The only thing andes legitimately keeps is the *published* rank-intensity LLR scoring **concept** (PepNovo/MS-GF+ papers — public-domain science, not protectable). Everything that is *MS-GF+'s specific expression, data, or structure* is removed or independently re-derived: the generating-function patent (already gone), the Java source (gone), the seed model weights (retrain), the model scaffold (build from own data via a new `--from-scratch` path), and the git fork lineage (orphan history).

**Tech Stack:** Rust workspace (`crates/scoring`, `crates/model-train`, `crates/andes`); Codon SLURM for harvest/train; VM for benchmark; Percolator 3.7.1; PRIDE public datasets harvested with pridepy + MSFragger + own combined-FDR + own trainer.

**The five independence buckets (and which phase closes each):**
1. Patent (generating function) — **closed** (Phase 0 verifies).
2. Copyright (Rust expression) — Phase 3 (clean-room scrub + attestation).
3. Data (model weights) — Phase 2 (retrain every model on own public data).
4. Structure (model scaffold) — **Phase 1** (the real code gap: `--from-scratch`).
5. Contract (the UC license terms) + lineage (git history) — Phase 5 + Phase 4.

---

## Phase 0 — Independence baseline & guard rails (do first)

**Files:** Create `docs/independence/2026-06-12-baseline-audit.md`; Create `.github/workflows/independence.yml` (or CI equivalent).

- [ ] **Step 0.1 — Snapshot what is still MS-GF+.** Run the patent grep on the release-candidate branch (`generating|gf|spec_e|evalue|score_dist|denovo|convolv`); dump the `source` column of every model store; list every `--seed-model` call site. Record results in the baseline-audit doc. Expected: patents clean; ~10/39 models own; 100% of scaffold seed-derived.
- [ ] **Step 0.2 — Add CI guard `no-gf`.** Fails build if any live (non-comment, non-test-assertion) symbol from the patent list appears. Run it; confirm green on current tree.
- [ ] **Step 0.3 — Add CI guard `no-java-legacy`.** Fails if `java-legacy*` refs or `.java` files exist in the tree.
- [ ] **Step 0.4 — Add CI guard `no-seed-source`** (stub now, enforced after Phase 2): fails if any shipped model store has a model with `source != "msnet"`/own. Mark `allow-fail` until Phase 2 completes.
- [ ] **Step 0.5 — Commit** the audit + guards.

---

## Phase 1 — `--from-scratch` scaffold builder (THE code gap)

**Goal:** Build a `Param` from (config + own training data) with **no seed input**. Today `train-from-msnet` loads the seed Param and replaces only the 4 learned tables; the scaffold (`data_type, mme, num_segments, max_rank, error_scaling_factor, min/max_charge, partitions, charge_hist, precursor_off_map`) is copied from MS-GF+. After this phase it is derived from your config + data.

**Files:**
- Create: `crates/model-train/src/scaffold.rs` — `build_scratch_param(cfg: &ScaffoldConfig, mass_charges: &[(f64, u8)]) -> Param`
- Modify: `crates/model-train/src/lib.rs` (export `scaffold`)
- Modify: `crates/andes/src/bin/andes.rs` (`TrainFromMsnetArgs` + `run_train_from_msnet` ~2482)
- Test: `crates/model-train/src/scaffold.rs` (`#[cfg(test)] mod tests`)

### Design decision (locked): **data-quantile partitions**
Per-charge parent-mass bins are **equal-population quantiles of the training corpus's precursor-mass distribution** (your data defines your structure). `num_segments`, `mme`, `max_rank`, `error_scaling_factor`, charge range are config scalars (your design choices, not MS-GF+ secrets). `charge_hist` and `precursor_off_map` are counted from the data.

- [ ] **Step 1.1 — Write the failing test for the mass-bin deriver.**
```rust
#[test]
fn quantile_mass_bins_are_equal_population_and_monotonic() {
    // 1000 synthetic precursor masses, charge 2
    let masses: Vec<f64> = (0..1000).map(|i| 500.0 + i as f64).collect();
    let edges = super::quantile_mass_edges(&masses, 4); // 4 bins
    assert_eq!(edges.len(), 3, "4 bins -> 3 interior edges");
    assert!(edges.windows(2).all(|w| w[0] < w[1]), "edges strictly increasing");
    // each bin ~250 masses (equal population)
    let b0 = masses.iter().filter(|&&m| m < edges[0]).count();
    assert!((230..=270).contains(&b0), "first bin ~250, got {b0}");
}
```
- [ ] **Step 1.2 — Run it, confirm FAIL** (`cargo test -p model-train quantile_mass`).
- [ ] **Step 1.3 — Implement `quantile_mass_edges(masses, n_bins) -> Vec<f64>`** (sort, take n_bins-1 quantile cut points). Run, confirm PASS.
- [ ] **Step 1.4 — Write failing test for `ScaffoldConfig` + `build_scratch_param`.** Assert: the returned `Param` has the config's `mme/num_segments/max_rank/esf/charge range`; `partitions` cover `charge × quantile-mass-bins × segments` with no overlap/gap; `charge_hist` reflects the input charge counts; `precursor_off_map` is populated; and the learned tables are *empty/uniform* (to be filled by the estimator). Crucially assert **it takes NO `Param`/seed argument**.
- [ ] **Step 1.5 — Run, confirm FAIL.**
- [ ] **Step 1.6 — Implement `ScaffoldConfig` (struct of the scalar knobs + enzyme/cleavage spec) and `build_scratch_param`.** Build `data_type` from config; `partitions` from `quantile_mass_edges` per charge × `num_segments`; `charge_hist`/`precursor_off_map` from the `mass_charges` input; allocate empty learned tables sized by `max_rank`/`esf`. Run, confirm PASS + `cargo clippy -p model-train -- -D warnings`.
- [ ] **Step 1.7 — Wire the CLI.** Add `--from-scratch` (bool) + scalar args (`--num-segments`, `--max-rank`, `--error-scaling-factor`, `--mass-bins`, `--charge-min/max`, `--enzyme`, etc.) to `TrainFromMsnetArgs`. In `run_train_from_msnet`: if `--from-scratch`, FIRST pass over `psms` to collect `(precursor_mass, charge)`, call `build_scratch_param`, and use it instead of `load_seed_param` (make `--seed-model` optional/mutually-exclusive). The existing accumulate + `Estimator` then fill the tables. Keep the seed path for back-compat until Phase 2 flips the default.
- [ ] **Step 1.8 — Integration test:** train a tiny model `--from-scratch` from a fixture flat; assert it writes a store, has the right model_id, non-empty tables, and `source` = own. `cargo test -p andes from_scratch`.
- [ ] **Step 1.9 — Commit** (`feat(model-train): --from-scratch builds the model scaffold from own config + data-quantile partitions; no seed input`).

---

## Phase 2 — Retrain EVERY model on own data, from scratch (the campaign)

**Gated procedure (uses the harvest/train infra on Codon). Not TDD.**

- [ ] **Step 2.1 — Define your own model taxonomy.** Decide the slug set *you* support based on *your* available public data — do not mechanically reproduce MS-GF+'s 39. Drop exotic combos with no public data (option-B drops from the curation). Record the authoritative list in `docs/independence/model-taxonomy.md`.
- [ ] **Step 2.2 — Harvest each slug's corpus with the combined-FDR pipeline** (public PRIDE only; pridepy + MSFragger + combined Percolator FDR — your pipeline, your decoys). Reuse `harvest_comb.sh`. Log dataset provenance per slug.
- [ ] **Step 2.3 — Train each slug with `--from-scratch`** (Phase 1) + the tight-training matcher (c6084a97). No `--seed-model` anywhere.
- [ ] **Step 2.4 — Per-slug acceptance gate:** benchmark each model; it must **meet-or-beat** the prior best (seed where applicable, else other tools) on a held-out set — premise: no regressions. Resolve the high-res `hcd_qexactive_tryp` gap here (clean corpus + from-scratch is the live hypothesis; the bar is ≥ seed 70,887 on Astral).
- [ ] **Step 2.5 — Assemble one all-own `models.parquet`** from the per-slug stores; every model `source = own/msnet`. **Delete `seed-models.parquet` from the repo/distribution.**
- [ ] **Step 2.6 — Flip the bundled default** (`resources/ionstat/models.parquet`) to the all-own store; regression-test the benchmark board (UPS1/TMT/Astral) still beats all engines.
- [ ] **Step 2.7 — Enforce CI guard `no-seed-source`** (remove the allow-fail from Step 0.4).

---

## Phase 3 — Clean-room code finalization (copyright)

**Gated review. Files: across `crates/`, plus `LICENSE`, `NOTICE`, `README.md`.**

- [ ] **Step 3.1 — Scrub stale MS-GF+ doc-comments** (the `GF`/`SpecEValue`/`DeNovo` descriptive comments the audit flagged in `match_engine.rs`, `coisolation.rs`, `pin.rs`). They describe logic that no longer exists.
- [ ] **Step 3.2 — Comment/identifier sweep:** grep for variable names, comments, or struct layouts that mirror the MS-GF+ Java too closely; rename/rewrite to independent expression. Document anything intentionally parallel (functional necessity) in the dossier.
- [ ] **Step 3.3 — Remove `java-legacy*` / `primitives-gf*` branches from the release remote**; confirm the release branch descends only from clean-room commits.
- [ ] **Step 3.4 — Clean-room attestation:** write `docs/independence/clean-room-attestation.md` — who authored the Rust, from what sources (papers/specs, NOT the Java), and that no Java was translated line-by-line.

---

## Phase 4 — Split the git history (lineage)

**Goal:** no fork lineage to MS-GF+ in the shipped repo.

- [ ] **Step 4.1 — Decide the mechanism** (recommend: fresh **orphan** root). Confirm with the user before destructive history ops.
- [ ] **Step 4.2 — Create an orphan branch** with a single squashed commit of the current clean-room tree (`git checkout --orphan release-clean && git commit`), OR export the tree to a brand-new repo. No MS-GF+ commits, authors, or tags carried over.
- [ ] **Step 4.3 — Verify** `git log` shows only own history; no `.java`, no `java-legacy` refs, no `seed-models.parquet`, no MS-GF+ author lines.
- [ ] **Step 4.4 — Preserve the old repo privately** (internal only, never published) for provenance/audit, separate from the public artifact.

---

## Phase 5 — License & legal (the gate that can't be engineered around)

**Non-engineering. Do in parallel; it gates the actual relicense.**

- [ ] **Step 5.1 — Read the actual MS-GF+/UC license** and record its terms on derivative works, relicensing, and commercial use in the dossier.
- [ ] **Step 5.2 — Hand the independence dossier to tech-transfer / IP counsel** (patent audit + clean-room attestation + model-provenance ledger + scaffold-derivation note from Phase 1) and obtain **written sign-off** that the reauthoring + own models are unencumbered.
- [ ] **Step 5.3 — Update `LICENSE`/`NOTICE`/`README`** to Apache-2.0 with accurate provenance (clean-room code, own-data models, no patented method) **only after** Step 5.2.

---

## Phase 6 — Provenance dossier & final verification

**Files: `docs/independence/dossier.md` + CI.**

- [ ] **Step 6.1 — Model-provenance ledger:** table of every shipped model → public PXD(s) → harvest pipeline → `source=own`. Generated from the store's source ledger, not hand-maintained.
- [ ] **Step 6.2 — Assemble the dossier** (Phase 0 audit + Phase 3 attestation + Step 6.1 ledger + Phase 1 scaffold-derivation note + Phase 5 license read + sign-off).
- [ ] **Step 6.3 — Full independence re-verification** on the release artifact: `no-gf`, `no-java-legacy`, `no-seed-source` all green; `git log` lineage-clean; benchmark board still beats all engines on own data.
- [ ] **Step 6.4 — Tag the Apache release.**

---

## Self-review notes
- The *only* MS-GF+ thing intentionally retained is the published rank-LLR scoring concept (public-domain science) — explicitly justified, not copied expression.
- Phase 1 is the sole net-new feature and the real engineering gap; everything else is retrain (2), scrub (3), history (4), legal (5), docs (6).
- Phases 2 and 5 run longest; start the harvest/train campaign and the license read immediately and in parallel with Phase 1.
- Acceptance premise throughout: **no regression below the prior best** (Phase 2.4 gate).
