# Glyco review brief — algorithmic avenues to close the intact-N-glycopeptide gap

**Audience:** a fresh reviewer (human or AI) tasked with finding *algorithmic* ways to push andes' intact-N-glycopeptide performance beyond the current ~55%.

## State of play
- Benchmark: public AI-ETD dataset **PXD011533** (mouse brain, 6 fractions), scored vs a 5088-glycopeptide reference identification set, **backbone-correct @1% FDR (Percolator)**.
- **andes ≈ 55%** (2820/5088). Reference specialist tools ≈ **71%**. Remaining gap ≈ 815 IDs.
- Trusted diagnosis: the deficit is **broad**, not one hole. Backbone-correct hit-rate ≈ **48% ETD / 59% HCD** (both below reference). Gap split ≈ **45% never generated / 44% generated-but-outranked / 10% FDR** ⇒ anything touching *only* FDR caps near +1–2 pp.
- Verdict: the remaining gap is **not a physics ceiling** (raw c/z presence rises with charge; the reference tool hits ~51% z5 / ~29% z6 on the same data). It is an **architectural ceiling of the current GPSM** plus a small composition-unreachable z6 subset.

## Read first
- Memory: `MEMORY.md`, `reference_glyco_scoring_architecture.md`, `reference_highcharge_rootcause.md`, `project_glyco_hybrid_campaign.md` (in the session memory dir).
- Code: `crates/search/src/glyco_search.rs` (generation → collapse → PIN), `crates/andes-glyco/` (`backbone.rs`, `hybrid.rs`, `glycan_db.rs`, `sequon.rs`, `oxonium.rs`, `glyco_psm.rs`), `crates/scoring/src/scoring/psm_score.rs`, `crates/output/src/glyco_pin.rs`.

## Hard constraints
- FDR is **Percolator PSM-level only**. Express any 2D/glycan-axis idea as **additive PIN features or decoy rows**, never a second FDR engine, never Mokapot.
- **Never modify the existing `rank`/PIN columns** — additive only. Replacing `RankScore` with a c/z model cratered FDR 3× (proven).
- Own training data only; no patented algorithms. Own-data model training *is* allowed (additive only).
- Validate **decoy-safe @1%, pooled fractions, one variable per A/B**. Measure before claiming magnitude.

## Already shipped / decoy-safe (don't redo)
Round-2 reweight + c/z truncation gate (34→51%); intensity into hyperscores (`CzIntensity`, +25); ETD oxonium-gate un-gating (default-on, +19); `chargeHi`; sequon-boundary generation (+6); charge-reconciled + `PAIR_Y_ON_GEN` paired path; paired HCD+ETD (`--glyco-hcd-pair`, +44, biggest single lever).

## Refuted — do NOT retry as-is
Coarse space widening (full 4036-glycan list, isotope-range widening, charge z±1) — all **dilute FDR**. `sialic_consistency` into the collapse (mouse brain is oligomannose-rich). Raising the deconvolution charge range. **Retraining the fragment model as a rank *replacement*** (regresses Percolator every time). The old closed-HCD-only evidence wall.

## Ranked open avenues (build order; magnitudes are hypotheses — measure)
1. **Joint paired HCD+ETD GPSM** (high; partly a wiring bug). Score the HCD glycan channel (Y-ladder, partial-glycan b/y, oxonium fingerprint) on the HCD partner and c/z on the ETD scan, fused into collapse + additive PIN. First step: `ANDES_GLYCO_PAIR_Y_ON_GEN` (glycan-Y read from `gen_peaks` when paired) — directly targets the z5 pairing regression. Precedent: O-Pair (Lu 2020), pGlyco3 (Zeng 2021), MSFragger-Glyco (Polasky 2020).
2. **Promote sequence-specific graded terms into the fused selector** (high). `partial_glycan_by`, `CzIntensity`, `y0y1_anchor`, `strong_score` are computed but *not in the argmax*; the two dominant collapse terms (`K·ladder`, `J·core_y`) are per-backbone (zero isobaric-peptide discrimination — round-2 finding). Offline separability check on `--debug-glyco` PINs first. Precedent: MSFragger multiattribute (Polasky MCP 2022), Byonic intensity + absence penalties (Bern MCP 2021).
3. **Decorated-backbone additive PIN** (medium). Bare-backbone matching hides ~half the ladder (glycosite-spanning ions carry the glycan). Emit `NumMatchedIonsGlyco`/`ExplainedIonCurrentGlyco` via `glyco_aware_peptide`; never replace `RankScore`.
4. **Real glycan-axis competition as Percolator features** (medium-low alone; enabler). Default-on glycan decoys with *recomputed* (not copied) glycan features + `DeltaGlycanScore` + size-normalized Y-completeness + categorized-oxonium log-likelihoods. FDR slice is only ~10% of the gap.
5. **Targeted generation leftovers** (partial). Multi-sequon max-cz, parse `MS:1000633` possible charges (measure first), peptide-first **c/z index** on unpaired ETD (index is b/y-only today), chimeric-in-glyco. Classify the ~45% never-generated with a `--debug-glyco` histogram (sequon / wrong-z / glycan-not-in-DB / peptide-not-in-digest / truncation / oxonium) → one lever per bucket. Part of z6 is a hard composition-alphabet ceiling.
6. **Additive own-data glyco intensity / c/z model** (unknown, high cost). New PIN + fused term only — never replace rank. Do after 1–3.

## Validation rig
`ssh <bench-host>`, dir `<bench-root>/ethcd` — 6-fraction mzML, truth TSVs, `mouse-decoy.fasta`, `eth_bench_eval.py` (@1% evaluator), dockerized Percolator. Build: `cd /srv/data/msgf-bench/cz-src && cargo build --release -p andes`. Run: `andes --spectrum FracN.mzML --database mouse-decoy.fasta --decoy-strategy none --decoy-prefix DECOY_ --glyco --glyco-backbone-top-k 150 --output-pin OUT.glyco.pin` (`--debug-glyco` for all candidate rows → never-generated vs outranked). Access is VPN-dependent and drops.

## Output expected from a review
A ranked list, each avenue with: mechanism, reference-tool precedent + citation, code location, expected magnitude *and why*, additive/FDR-safe check, and a concrete one-variable first experiment. Distinguish likely-ceiling from genuinely-open. Prefer "measure X to decide" over confident claims.
