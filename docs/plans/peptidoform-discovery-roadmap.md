# Andes peptidoform discovery roadmap

**Date:** 2026-10-11.  
**Status:** planning — derived from the Oct 2026 algorithm/strategy review (closed-search
scorecard, glyco structural gaps, proteoform literature, hybrid de novo / DB / open design,
spectral clustering, structured noise).  
**Constraint reminders:** Percolator-only PSM FDR (no Mokapot / no new FDR engine); additive
PIN only; own training data; one-variable A/Bs; pool small searches for stable FDR.

---

## Vision

Keep andes as a fast, intensity-aware, Percolator-first engine. Expand it into an
**evidence-triggered peptidoform discovery system**: the closed fragment-index path stays
primary; tags, open Δm, class templates (glyco), clustering, and structured-noise models
only fire when the spectrum justifies them. Report **peptidoforms** honestly; infer
proteoforms only with multi-peptide or top-down evidence — do not overclaim.

---

## Baseline (2026-10)

| Strength | State |
|---|---|
| Closed search | Faster and richer than Comet (Astral ~+50% PSMs @1%, ~0.7× wall) |
| Retrieval | Fragment-index + out-of-core |
| Expansion knobs | Chimeric, refine Pass-2, glyco class path |
| Gap | No open search; refine does not *discover* mods; glyco selection/pairing still imperfect |
| Empirical scar | Expand the wrong spectra with the wrong scorer → 0 new IDs @1% or fake yield (refine flood; glyco split-election) |

Authoritative numbers: [`docs/benchmarks/README.md`](../benchmarks/README.md).

---

## Design principles

1. **Hypothesis object, not tool silos** — backbone + alterations + support (tags, shifted ions, residual, cluster).
2. **Index first** — spectrum graph / tags aid retrieval and localization; they do **not** replace the fragment index.
3. **Selective expansion** — trigger on twilight + residual + tags; never always-on open/tag/de novo on every scan.
4. **Discovery → closed refine** — open finds Δm; refine locks them in (G-PTM-D pattern).
5. **Grouped FDR** — `HypClass` / `DeltaMassBin` so open/decoy-dense bins do not destroy closed power.
6. **Noise is structured** — use it for nulls, masks, residuals — not as the peptide ID itself.
7. **Clustering recovers repeats** — not singletons; judge **peptides**, not PSM inflation.
8. **Measure before shipping** — Astral do-no-harm; entrapment; peptides @ equal true FDP.

---

## Roadmap

### Phase 0 — Foundations (cheap probes)

| Probe | Question | Go / no-go |
|---|---|---|
| Tag coverage | Do ≥~80% of Astral true PSMs contain a correct length-≥4 tag? | If no, tags are low priority |
| Twilight residual | Among q≈0.01–0.20 closed misses, how often is residual+Δm real vs noise? | Sets Pass-B triggers |
| Pairing Y-on-HCD | Under `--glyco-hcd-pair`, score Y on HCD not ETD | Fix known wiring defect |
| Sibling/cluster ceiling | Fraction of unidentified spectra near a strong same-precursor neighbour | Caps clustering ROI |

**Also ship if green:** glyco Y-ladder / PIN Y\* on the HCD partner when paired (one-variable A/B).

### Phase 1 — Open discovery on the closed backbone

**Goal:** unexpected PTMs without wrecking unmodified power.

```
Pass A  closed fragment-index (unchanged)
Pass B  IF twilight ∨ high unexplained ∨ strong tag:
          open Δm on bounded candidates + shifted-ion score
Pass C  histogram Δm → feed --refine tier on confident proteins / twilight only
```

| Deliverable | Detail |
|---|---|
| Localization-aware open | Unshifted + shifted fragment matching; peptide-level Δm first |
| PIN columns | `HypClass`, `DeltaMass`, `DeltaMassBin`, `ShiftedIonFrac`, `UnexplainedFrac`, `TagCoverage` |
| Refine sink | Discovered frequent Δm → Pass-2 closed (no “refine all unidentified”) |
| FDR | Percolator + group post-process by `HypClass` / bin |
| Gate | Astral unmodified peptides not hurt; mod-rich set gains peptidoforms @ ≤1% true FDP |

**Non-goals for P1:** full de novo paths, multi-unexpected-mod, top-down, TagGraph clone.

### Phase 2 — Tags + selective hybrid retrieval

**Goal:** use de novo *partial* evidence to find DB peptides that the closed mass window misses.

| Deliverable | Detail |
|---|---|
| Cheap tag extractor | From peaks (no full DP rewrite of the hot path) |
| Tag→DB retrieval | Bounded string-index lookup; enzyme-optional |
| Integrate with Pass B | Tags trigger / bound open candidates |
| Gate | Gain on semi/miscleavage or mod-rich sets; Astral do-no-harm |

Full de novo sequence models (e.g. π-PrimeNovo-class) only after tags prove value — see
[De novo arm decision](#de-novo-arm-decision) below.

### Phase 3 — Spectrum neighbourhood (clustering / transfer)

**Goal:** recover weak repeats of peptides already seen well.

| Deliverable | Detail |
|---|---|
| v3a PIN-only | `ClusterSupport` / precursor-neighbour similarity (safest) |
| v3b Transfer | Strong seed → weak cluster members; conflict rules (SIMSI-style) |
| Gate | **Peptides** (not PSMs) ↑ at equal FDP; if only PSMs↑, do not ship hard transfer |

Optional later: spectral-network edges with Δm for mod variants (propagation-aware FDR).

### Phase 4 — Structured noise as null and cleaner

**Goal:** stop treating unexplained peaks as i.i.d. random.

| Deliverable | Detail |
|---|---|
| Run-level background model | Frequent unexplained m/z → down-weight / mask |
| Residual features | After best PSM: contaminant hits, unexplained frac |
| Quality gate | Contaminant-dominated / empty → skip expansion |
| Gate | Better target/decoy separation; fewer false tag/graph edges |

Not a peptide ID from noise alone.

### Phase 5 — Class templates + peptidoform reporting

**Goal:** generalize glyco; honest products.

| Deliverable | Detail |
|---|---|
| Fix glyco pairing/selection debt | Y-on-HCD; graded sequence-specific terms carefully |
| Template API | Class ions + mass residual + backbone channel + class PIN |
| Next classes (pick one) | Ub remnant, ADP-ribose, or middle-down histone combinatorial |
| Reporting | Peptidoform-level export; proteoform inference only with multi-peptide/TD evidence |

Re-benchmark: PXD011533 AI-ETD + current mouse-liver glyco on one binary.

---

## Success metrics (everywhere)

| Metric | Why |
|---|---|
| Peptides / peptidoforms @ ≤1% true FDP | Primary |
| Astral unmodified (do-no-harm) | CONGA criterion |
| Wall / memory vs current main | Do not lose the speed win |
| Entrapment on expanded classes | Especially open / refine / transfer |
| One variable per A/B | Hygiene |

---

## Explicit non-goals

- Global open on every spectrum as default
- Spectrum-graph as replacement for fragment index
- New FDR engine / Mokapot
- Full TopPIC / ProSight competitor
- Claiming “proteoform ID” from single tryptic peptidoforms
- Always-on multi-generator union (refine-flood mode)
- Shipping a deep de novo model as core before P1 open search pays off

---

## Dependency order

```
P0 probes + glyco Y-on-HCD fix
        │
        ▼
P0b ONNX substrate (#92 → #93)          ← shared infra before any DL arm
        │  measure AlphaPeptDeep offline first
        │  then thin andes-onnx (feature-gated)
        │
        ▼
P1 open Δm → refine sink + group PIN     ← main discovery lever
        │
        ├─► P2 tags (accelerate P1 retrieval)
        ├─► P3 clustering (recover weak repeats)
        └─► P4 noise null/mask (sharpen all scores)
                │
                ▼
        P5 class templates + peptidoform product
                │
                ▼
        Optional: deep de novo via ONNX / sidecar — see below
```

---

## Near-term “do next”

1. **[#92](https://github.com/bigbio/andes/issues/92) Step 1** — offline AlphaPeptDeep ONNX → append five PIN columns → Percolator A/B (Astral + UPS1, entrapment, contamination check). No Rust yet.
2. Phase 0 tag + twilight probes on Astral (no engine rewrite).
3. One-variable glyco fix: Y-ladder / PIN Y\* on HCD partner when paired.
4. If #92 pays (or we still want the substrate for other models): **Phase 0b** — thin ONNX crate / feature (see below).
5. Spec + spike for P1: shifted-ion open on twilight only + `DeltaMassBin` PIN + Astral do-no-harm gate.
6. Only after P1 moves peptidoforms on a mod-rich set: tags, then clustering PIN features.

---

## Phase 0b — ONNX substrate (before any de novo crate)

**Scaffold landed (2026-10-11):** workspace crate [`crates/onnx-infer`](../../crates/onnx-infer/)
— tool-agnostic API (`OnnxSession`, `ModelManifest`, `NamedTensor`), `--features tract`
backend (pure Rust), reserved `ort` feature, fixture tests green. Not yet wired into
the andes binary or AlphaPeptDeep PIN columns.

**quantms-rescoring alignment (follow-up):** `rescoring` module shares feature-name
contracts with [quantms-rescoring](https://github.com/bigbio/quantms-rescoring);
`peptdeep` packs OpenMS-compatible ONNX I/O; export scripts port AGBD→TreeEnsemble
and PeptDeep `.pth`→ONNX. Native AGBD remains the in-search hot path; ONNX is the
interchange / quantms path. Top-N-only batching is encoded as `TopNPolicy`.

**ONNX fixtures landed:** `intensity_mlp` (tract), `tree_leq` + `peak_gbdt` AGBD↔ONNX
TreeEnsemble (ort). PeptDeep RT/MS2/CCS export script is ready; full `.pth`→ONNX
run needs a clean env with ~2GB free (blocked once by disk-full on this machine).

**Yes — this should land before a PrimeNovo-style de novo arm.**
Issue [#92](https://github.com/bigbio/andes/issues/92) (AlphaPeptDeep MS2+RT ONNX as additive rescoring) and the linked
porting path (#93) are the right first consumer: additive PIN features, top-N-only
inference, feature-gated / off by default, same shape OpenMS already shipped.

### Why this before de novo

| Reason | Detail |
|---|---|
| Shared infra | Same loader serves AlphaPeptDeep *and* any later exported Transformer (de novo backbone, intensity models) |
| Smaller first model | MS2/RT predictors are simpler than π-PrimeNovo (custom PMC/CTC decode) |
| Fits andes change type | Additive PIN columns — repeatedly measured safe |
| Forces packaging decision | `tract` (pure Rust, Apache-2.0) vs `ort` (ONNX Runtime, broader ops) under the existing apache-core vs vendor-deps split |
| Issue discipline | #92 Step 1 can kill the feature with a measurement; do not build the crate “just in case” if offline A/B is flat **and** no other ONNX consumer is queued |

### What the crate is / next wiring

Crate name: **`onnx-infer`** (generic on purpose — reusable outside andes).

Already:

1. Loads ONNX (+ optional JSON manifest with provenance / task / `extras`).
2. Batched `f32` I/O via `NamedTensor`; tract backend behind `--features tract`.
3. No dependency from `andes` binary yet (opt-in later).

Still to do for #92:

1. Offline AlphaPeptDeep measurement (Step 1 of the issue).
2. Wire top-N-only inference → five additive PIN columns; RankScore untouched.
3. Optional `ort` backend if tract cannot run the AlphaPeptDeep graphs.
4. NCE / RT run-level calibration hooks (manifest `extras` + search-side fit).

**Not in v1:** training, Python, or a PrimeNovo decode loop. Export/custom ops stay outside.

### Sequencing vs #92

```
#92 Step 1 (Python offline) ──► pays? ──yes──► wire the five features into andes on onnx-infer
                         └──no───┬► still want de novo later? ──yes──► keep onnx-infer as a minimal loader
                                 └► else close #92; leave onnx-infer unwired
```

The `onnx-infer` scaffold has already landed, so a flat Step 1 does not remove it: it stops
the AlphaPeptDeep wiring and productization. If we still plan ONNX de novo within a quarter,
the scaffold stays as #93 infrastructure — keep it a **minimal** loader.

---

## De novo arm decision

**Question:** should andes grow a dedicated de novo crate inspired by
[π-PrimeNovo](https://www.nature.com/articles/s41467-024-55021-3)
(Zhang et al., *Nat Commun* 2025 — non-autoregressive Transformer + precise mass control;
~64% peptide recall on nine-species; up to ~89× faster in the paper's best case, PrimeNovo
without precise mass control against Casanovo V2 with beam search)?

### Recommendation: **not as a core crate yet; yes as a staged optional arm after P0b + P1**

| Option | Verdict | Why |
|---|---|---|
| New `andes-denovo` crate that reimplements π-PrimeNovo now | **No** | Wrong dependency order; needs ONNX/runtime substrate first; P1 open→refine is higher ROI for DB-centric workflows |
| ONNX substrate (#92/#93) then optional exported backbone | **Yes, when due** | Same path as AlphaPeptDeep; PMC/CTC decode may stay outside ONNX |
| Python sidecar for a twilight spike | **Yes, anytime** | Validates hybrid value without a crate |
| Own deep de novo crate (own weights, Apache) | **Maybe, after P0b+P1+tag probe** | Only if DB-absent / metaproteome / antibody use cases remain |

### How a de novo arm should plug in (if/when)

```
closed search
    │
    ├─ identified → done
    └─ twilight / unexplained
            │
            ▼
     de novo model (ONNX and/or sidecar)
            │
            ├─ tags / full sequence → tag→DB + open Δm (Pass B)
            └─ DB-absent sequence → report as de novo peptidoform
                 with strict QC (mass error, decoy/entrapment policy)
```

De novo is a **generator for the hybrid hypothesis object**, not a replacement for
Percolator-scored DB search. Metaproteomics / antibody / novel antigen are the
justifying products — not “more Astral PSMs.”

### Practical next step on de novo

1. Land Phase 0b path (measure #92 → ONNX loader).
2. Finish P0 tag-coverage probe (classical tags, no DL).
3. Spike: run π-PrimeNovo (or Casanovo) externally on Astral twilight spectra; measure
   how often the prediction (a) matches a DB peptide outside the closed hit, (b) yields
   a correct open Δm localization, (c) is DB-absent but mass-consistent.
4. If (a)+(b) move peptidoforms under grouped FDR → export what ONNX can hold; wire
   through the same substrate; train/fine-tune on **own** labels before claiming independence.
5. Do **not** vendor third-party weights into the Apache core without license + training-data provenance review.

---

## One-line summary

**Keep the fast closed engine; land a thin ONNX substrate via #92/#93 before any DL arm;
discover unexpected chemistry with selective open→refine; use tags, clusters, and
structured noise as support channels; keep glyco as the first class template; treat deep
de novo as an optional generator on that substrate after P1; always prove peptide-level
gains under grouped Percolator FDR.**
