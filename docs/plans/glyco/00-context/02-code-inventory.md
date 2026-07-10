# andes glyco — code inventory (what exists vs pending)

*Precise source-cited audit of the glyco code. Branch `glyco-phase1`, HEAD `3ce6f6f7`.
Repo: `/Users/yperez/work/msgfplus-workspace/msgf-rust`.
Companion to `00-current-state.md` (measurements) — this file is the CODE map.*

Constraints honored in all recommendations below:
FDR = **Percolator only** (2D-FDR = thin post-process); algorithms from **published
papers only** (a glyco search engine / a cross-spectrum glyco engine = Apache-2.0, O-Pair/an open-source glyco engine = MIT-permissive
— clean-room OK; **no a commercial glyco engine / the reference engine code**); andes differentiates as
**glycan-Y-first + own learned models + in-process cross-spectrum**, not a a comparison search engine/the reference engine clone.

---

## 1. Pipeline stages + env toggles

The glyco path is a **standalone driver** invoked from `--glyco`; it writes its own
`.glyco.pin` and skips the standard PIN/rescore/TSV/Parquet/refine pipeline
(`crates/andes/src/bin/andes.rs:2301` early-return block; guard at `:1762`,
`:2297-2356`). Percolator is run on the `.glyco.pin` separately.

### Stage flow (per spectrum, in `glyco_search_run`)

| # | Stage | Function (file:line) | Notes |
|---|-------|----------------------|-------|
| 0 | CLI entry | `andes.rs:2301` `if cli.glyco` | glycan list = `n_glycan_list_common()` (~600), `tol_ppm=20.0` hardcoded (`andes.rs:2307-2308`) |
| 1 | Oxonium gate | `oxonium.rs:10` `oxonium_gate` | fires iff summed core-oxonium frac ≥ 0.10 AND ≥2 of 5 ions (`glycan_mass.rs:24` `CORE_OXONIUM_MZ`) |
| 2 | Backbone generation (per charge × isotope offset) | `hybrid.rs:159` `hybrid_candidates_with_isotope` | Y-ion-first cascade; see §4 |
| 2a | — de-novo Y-ladder solver | `backbone.rs:109` `solve_backbone_min` | quorum-2 primary → quorum-1 rescue (`hybrid.rs:200-203`) |
| 2b | — DB-union + core-Y truncation | `hybrid.rs:247-324` | union `db_branch` when Y-first has evidence, rank by `core_y_intensity`, truncate to `top_k` |
| 2c | — peptide-first union | `glyco_search.rs:370-416` | `FragmentIndex.query` (≥`MIN_BY_MATCHES=6` b/y) → glycan-by-subtraction; cap `MAX_PEPTIDE_FIRST=64` |
| 2d | — glycan-Y-first index (Phase G1) | `glyco_search.rs:423-450` + `glyco_y_index.rs:91` `GlycanYIndex::query` | opt-in `ANDES_GLYCO_YINDEX` |
| 2e | — cross-spectrum transfer (pass 2) | `glyco_search.rs:801-882` + `crossspectrum.rs:58` `transfer` | opt-in `ANDES_GLYCO_CROSSSPECTRUM` |
| 3 | Dedup backbone hits | `glyco_search.rs:123` `dedup_backbone_hits` | merges only same (backbone, glycan-hypothesis, isotope) |
| 4 | Phase-1 cheap b/y scoring (ALL backbones) | `glyco_search.rs:529-605` | `score_psm` + `psm_edge_score`, sequon filter `has_nxst_sequon`; per-backbone best rank |
| 5 | Two-axis backbone selection + top-k cap | `glyco_search.rs:621-656` | AXIS-1 = b/y rank; AXIS-2 = core-Y (only if `yindex_on`) |
| 6 | Phase-2 expensive features (top-k only) | `glyco_search.rs:669-778` | `compute_psm_features`; builds `GlycoPsmKey` incl. Y-ladder + decoy scores |
| 7 | PIN emission | `glyco_pin.rs:264` `write_glyco_pin` | 6 glyco columns + optional glycan-decoy rows |

### Env toggles (all read in `glyco_search.rs:200-238`, plus PIN decoy in `andes.rs:2337`)

| Env var | Default | Effect | file:line |
|---------|---------|--------|-----------|
| `ANDES_GLYCO_PEPTIDE_FIRST` | **ON** (`!=0`) | peptide-first fragment-index union + build (`FragmentIndex`) | `glyco_search.rs:200`, guard `:253`, use `:370` |
| `ANDES_GLYCO_CROSSSPECTRUM` | OFF (`==1`) | pass-2 sibling backbone transfer | `glyco_search.rs:203`, use `:797` |
| `ANDES_GLYCO_DECOY` | OFF (`==1`) | compute glycan-axis decoy Y-ladder per hit **and** emit paired decoy PIN rows | `glyco_search.rs:209` (compute), `andes.rs:2337` (emit) |
| `ANDES_GLYCO_SCANS=<file>` | unset | dev harness: subset driver to listed scan numbers (truth-scan A/B, ~8 min/arm) | `glyco_search.rs:216`, use `:308` |
| `ANDES_GLYCO_EXHAUSTIVE` | OFF (`==1`) | disable both truncations (`top_k=100_000`) to measure ceiling; SLOW | `glyco_search.rs:227`, use `:230` |
| `ANDES_GLYCO_YINDEX` | OFF (`==1`) | Phase-G1 glycan-Y index generation + two-axis retention | `glyco_search.rs:236`, build `:239`, use `:423`,`:641` |

CLI flags: `--glyco` (`andes.rs:391`), `--glyco-backbone-top-k` (default 50, hidden,
`andes.rs:398`), `--glyco-max-spectra` (dev, hidden, `andes.rs:402`).

---

## 2. The scoring seam (score_psm + SelectionKey + protocol store)

**This is the core of the SP-B gap.** The glyco driver scores every backbone as a
**bare, unmodified peptide** using the SAME model the standard search resolved:

- `glyco_search_run` takes `prepared: &PreparedSearch` and uses `prepared.scorer`
  (`glyco_search.rs:191`), which is the **intact standard model** loaded once at
  `andes.rs:1532` (`load_param_from_store`) → `RankScorer::new` (`andes.rs:1613`),
  then passed to `glyco_search_run` at `andes.rs:2316`.
- Phase-1 ranking = `score_psm(ss, &cand.peptide, scorer, z, frag_tol)`
  (`glyco_search.rs:577`, defined `crates/scoring/src/scoring/psm_score.rs:234`)
  `+ psm_edge_score` (`:578`, def `psm_score.rs:58`). The glycan mass is **not**
  added to Asn — Percolator sees only backbone b/y ions (module doc `glyco_search.rs:8-12`).
- Glycan evidence is **additive PIN columns only** (`GlycoPsmKey`, §5), never fused
  into the score that ranks candidates.

**Model selection** (`build_selection_key`, `andes.rs:4762`) maps a
`model::protocol::Protocol` (`crates/model/src/protocol.rs:4`) to an
`experiment_class` `BTreeSet` via `protocol_to_experiment_class`
(`store/read.rs:255`), then `select` / `select_nearest`
(`select.rs:192`, `:261`) walk the backoff ladder over `SelectionEntry` rows
(`select.rs:31`). The store is partitioned by protocol
(`resources/models/protocol={Automatic,TMT,Phosphorylation,iTRAQ}`).

**There is NO glyco variant anywhere in this seam** — this is the wiring gap:
- `model::protocol::Protocol` has 6 variants, none glyco (`protocol.rs:4-11`).
- The CLI `Protocol` enum (`andes.rs:69`) has no glyco value.
- `protocol_to_experiment_class` (`store/read.rs:255`) has no glyco arm.
- `build_selection_key` (`andes.rs:4809`) has no glyco arm.
- No `resources/models/protocol=NGlyco` partition exists.
- The training catalog DOES know a `"glyco"` slug (`catalog.rs:95`,
  `inference: None` — explicit-tag-only), and the store schema already reserves
  **`loss_class=1` = glyco** (`store/schema.rs:211`, columns `ion_loss_class`
  `:213` + `frag_off_loss_classes` `:222`; read `store/read.rs:491-558`,
  write `store/write.rs:393-470`). So the *store format* supports glyco fragment
  offsets today; nothing *produces or selects* them.

Net: the glyco PSM is scored by a generic intact-peptide model. That intact model
never saw glyco spectra (oxonium-dominated, sparse backbone b/y, Y-ladder), so
top-1 backbone ranking is weak — **the measured 83/154 top-1 failure**
(`00-current-state.md §1`). SP-B/G2 = give glyco its own model + route to it.

---

## 3. What SP-B / G2 needs (concrete wiring)

SP-B = a glyco-specific SCORING/ranking model + its selection. Four wiring points,
each a small, additive change that respects the existing seam:

### (a) `Protocol::NGlyco` variant
- Add `NGlyco` to `model::protocol::Protocol` (`crates/model/src/protocol.rs:4`)
  with `name()`/`from_name()` arms (`:14`, `:26`). (O-glyco later = second variant;
  keep N first — andes is N-glyco today, `sequon.rs` is N-X-S/T only.)
- Add `#[clap(name = "N-glyco")] NGlyco` to the CLI enum (`andes.rs:69`) and map it
  through wherever the CLI `Protocol` → `model::Protocol` conversion happens.
- Note: `--glyco` is a *mode* flag today (`andes.rs:391`); SP-B can either keep
  `--glyco` and force protocol=NGlyco internally, or expose `--protocol N-glyco`.
  Keeping `--glyco` implies the protocol is the cleaner minimal change.

### (b) `build_selection_key` arm
- Add a glyco arm in `build_selection_key` (`andes.rs:4809-4815`):
  `Protocol::NGlyco => "NGlyco"` (or `"glyco"`) as `protocol_for_store`.
- Add the matching arm to `protocol_to_experiment_class` (`store/read.rs:255-268`):
  `"NGlyco" => parse_experiment_class("glyco")` (single opaque slug — mirrors the
  `iTRAQPhospho` single-slug precedent at `:262` so `select` step-1 finds a glyco
  model and falls through cleanly when absent). Catalog slug is already `"glyco"`
  (`catalog.rs:95`).
- `select_nearest` (`select.rs:261`) already routes an absent protocol to the
  standard base with a WARN — so an incremental rollout (variant added, model not
  yet trained) degrades gracefully to today's behaviour.

### (c) Glycan-stripped training rows
- Training accumulates fragment/rank statistics from **labeled (peptide, spectrum,
  charge)** rows via `accumulate` (`accumulate.rs:62`), which calls
  `ScoredSpectrum::new` + `ion_match_facts` (`accumulate.rs:73,80`) — a **bare
  peptide** matcher. For glyco this is exactly right IF the peptide is the
  backbone and the spectrum's glycan Y-ions/oxonium are treated as noise (or as
  `loss_class=1` offset ions). Two sub-tasks:
  1. Produce labeled backbone rows: `(backbone_peptide, glyco_spectrum, charge)`
     where the peptide is the **glycan-stripped backbone** (the Asn carries no
     glycan mass). The spectrum is the raw stepped-HCD/HCD glyco MS2.
  2. Optionally register glyco Y-ion offsets as `loss_class=1` ions so the model
     LEARNS the trimannosyl-core ladder instead of penalizing it as missing
     (schema already supports it, `store/schema.rs:211`; no producer exists).
- No code strips glycans for training today — `labeled.rs:111` `bootstrap_labels`
  and `accumulate.rs` operate on whatever `Peptide` they are handed. The stripper
  + row generator is NEW code (a `train`-side glyco corpus builder), gated by the
  `loss_class=1`/`"glyco"` slug that already exist.

### (d) truth-TSV → parquet converter
- The measurement harness uses truth scans (the reference engine-nglycan re-search of
  PXD025455, `00-current-state.md §1`) as TSV. To TRAIN, these
  `(scan, backbone_peptide, glycan, charge)` truth rows must become the labeled
  input to `accumulate`. No converter exists in-tree (the store writer
  `store/write.rs` writes the *trained model*, not labeled input). NEW code:
  a small truth-TSV → labeled-rows adapter feeding `accumulate` → `ModelStore`
  written under `protocol=NGlyco`. The parquet **model** schema is ready
  (`store/schema.rs`); only the labeled-input path is missing.

**FDR note (hard constraint):** SP-B is ranking only. The 2D-FDR (G3) stays a
Percolator post-process — the glycan-axis decoy rows already emitted by
`write_glyco_pin` (`glyco_pin.rs:320-330`) are the mechanism; do NOT add an
andes-internal FDR. G3 is a separate lever from SP-B.

---

## 4. Cross-spectrum scaffold state

The **andes-unique** lever (a cross-spectrum glyco engine-style backbone transfer;
Apache-2.0 ref, clean-room). Fully built, opt-in, wired into the driver:

- `crossspectrum.rs:25` `GlycoformWhitelist` — sorted/deduped confident backbone
  masses. `new` (`:33`), `transfer` (`:58`, binary-search glycan lookup),
  `nearest_glycan` (`:85`). 3 unit tests pass (`:119-155`).
- Driver integration is a **two-pass** design (`glyco_search.rs:790-887`):
  - PASS 1 (`:791`) = normal generation, no transfer.
  - Build whitelist from pass-1 PSMs with `core_y_hits ≥ CONF_MIN_CORE_Y=3`
    (`:804-814`).
  - PASS 2 (`:831`) = only non-confident oxonium-positive spectra; inject
    transferred backbones via `whitelist.transfer` (`:856`) into the same
    dedup/score path (`process_one(..., &transfer)` at `:876`); pass-2 results
    supersede pass-1 (`:880-882`).
- **State: complete but gated OFF** (`ANDES_GLYCO_CROSSSPECTRUM`, default false,
  `glyco_search.rs:203-205`; early return `:797-799` skips pass 2 entirely).
  Not yet A/B-measured on truth scans (unlike G1, which is verified —
  `00-current-state.md §2`). This is the differentiator to validate once SP-B
  lifts top-1 ranking (transfer only pays when the *donor* glycoforms are
  confidently ID'd, which is exactly what SP-B fixes).

---

## 5. File:line map

### `crates/andes-glyco/src/` (the primitives — mostly complete, well-tested)

| File | Key items (line) | Role |
|------|------------------|------|
| `lib.rs` | modules `:1-9` | crate root; 10 modules |
| `glycan_mass.rs` | `HEXNAC :5`, `HEX :8`, `FUC :11`, `NEUAC :14`, `NEUGC :17`, `PROTON :20`, `CORE_OXONIUM_MZ :24`, `CORE_Y_STEPS :34`, `MONO_STEPS :43` | published monoisotopic constants (Glyco-Fragment/UniMod) |
| `oxonium.rs` | `OxoniumEvidence :3`, `oxonium_gate :10` | glyco-spectrum gate (frac≥0.10, n≥2) |
| `sequon.rs` | `has_nxst_sequon :13` | N-X-S/T (X≠P); **N-glyco only** |
| `glycan_db.rs` | `GlycanComp :27`, `n_glycan_list :40` (~2510), `n_glycan_list_common :111` (~600, default) | clean-room combinatorial enumerator (no vendor list) |
| `backbone.rs` | `BackboneCandidate :13`, `complement_score :42`, `solve_backbone :93`, `solve_backbone_min :109`, `core_y_intensity :357`, `glycan_y_intensity :408`, `glycan_y_intensity_decoy :498` (splitmix64 `:479`), `count_core_y_hits :582` | de-novo Y-ladder solver + Y-ladder scorers + **G3 glycan-axis decoy** |
| `hybrid.rs` | `Source :20`, `BackboneHit :30`, `db_branch :58`, `nearest_glycan :118`, `hybrid_candidates :143`, `hybrid_candidates_with_isotope :159` | DB ∪ de-novo union, cascade, isotope sweep |
| `glyco_y_index.rs` | `GlycanYIndex :51`, `build :66`, `query :91`, `core_partials :34`, `has_core :46` | **Phase-G1** glycan-Y-complementary index (a glyco search engine insight, clean-room) |
| `crossspectrum.rs` | `GlycoformWhitelist :25`, `new :33`, `transfer :58` | **andes-unique** cross-spectrum transfer (see §4) |
| `glyco_psm.rs` | `GlycoPsmKey :36` | carrier of glycan evidence → PIN (all fields listed `:38-59`) |

### `crates/search/src/glyco_search.rs` (the driver — the SP-B locus)

| Item | line | Role |
|------|------|------|
| `FullGlycoPsm` / `GlycoSpectrumResult` | `:54` / `:63` | outputs |
| `MIN_GLYCAN` | `:70` | 406 Da (2×HexNAc core) |
| `nearest_glycan_mass` | `:76` | glycan-by-subtraction lookup |
| `glycan_decoy_seed` | `:115` | stable per-composition G3 seed |
| `dedup_backbone_hits` | `:123` | cross-charge/isotope merge (hypothesis-aware) |
| **`glyco_search_run`** | **`:184`** | main driver |
| env toggles | `:200-238` | see §1 |
| `FragmentIndex::build` (peptide-first) | `:253-282` | `MIN_BY_MATCHES=6 :293`, `MAX_PEPTIDE_FIRST=64 :296` |
| `process_one` closure | `:300` | per-spectrum body (both passes) |
| backbone gen (charge×iso) | `:333-450` | hybrid + peptide-first + Y-index |
| **phase-1 b/y scoring** (`score_psm`/`edge`) | **`:577-578`** | ← SP-B ranking call site |
| two-axis selection + top-k | `:621-656` | |
| **phase-2 features + `GlycoPsmKey` build** | `:691-778` | Y-ladder `:752`, decoy `:760` |
| cross-spectrum pass 2 | `:790-887` | see §4 |

### `crates/output/src/glyco_pin.rs` (PIN writer)

| Item | line | Role |
|------|------|------|
| `write_glyco_header` | `:33` | standard cols + 6 glyco cols `:107-112` (OxoniumScore, NCoreOxoniumIons, YLadderScore, CoreYHits, GlycanMass, IsGlycanDb) |
| `write_glyco_psm_row` | `:124` | `CalcMass = peptide.mass()+glycan_mass :164`; glycan-decoy override `:138` |
| `write_glyco_pin` / `_to` | `:264` / `:287` | G3 paired glycan-decoy emission `:320-330` (Label −1, `glycandecoy_` accession `:247`) |

### `crates/model-train/src/` + `crates/andes/src/bin/andes.rs` (the seam — glyco-blind today)

| Item | file:line | Glyco status |
|------|-----------|--------------|
| `Protocol` (model) | `model/src/protocol.rs:4` | **no glyco variant** — SP-B (a) |
| `Protocol` (CLI) | `andes.rs:69` | **no glyco variant** — SP-B (a) |
| `SelectionKey` / `SelectionEntry` | `model-train/src/select.rs:45` / `:31` | protocol-agnostic; ready |
| `select` / `select_nearest` | `select.rs:192` / `:261` | ready; `select_nearest` WARN-degrades absent protocol |
| `protocol_to_experiment_class` | `store/read.rs:255` | **no glyco arm** — SP-B (b) |
| `build_selection_key` | `andes.rs:4762`, arm `:4809` | **no glyco arm** — SP-B (b) |
| `load_param_from_store` | `andes.rs:4843` | resolves ONE intact model; glyco reuses it |
| `accumulate` (training) | `accumulate.rs:62` | bare-peptide matcher; needs glycan-stripped rows — SP-B (c) |
| `bootstrap_labels` | `labeled.rs:111` | label source; no glyco truth adapter — SP-B (d) |
| store schema `loss_class` | `store/schema.rs:211` (`1=glyco`) | **format ready**, no producer/consumer for glyco |
| catalog `"glyco"` slug | `catalog.rs:95` | slug known (`inference: None`) |
| store partitions | `resources/models/protocol=*` | **no `protocol=NGlyco`** — SP-B (d) output |

---

## Clean-room / license ledger (for the algorithms referenced above)

| Idea in code | Paper / source | License | Where used |
|--------------|----------------|---------|------------|
| Glycan-Y-complementary indexing (`glycan_mass − partial`) | a glyco search engine (Zeng et al., *Nat Methods* 18:1515, 2021; PMC8493561; DOI 10.1038/s41592-021-01306-0). Code github.com/pFindStudio/a glyco search engine (Apache-2.0) | Apache-2.0 (ref only) | `glyco_y_index.rs`, `backbone.rs` |
| Cross-spectrum glycoform transfer | a cross-spectrum glyco engine (Fang et al., *Nat Commun* 13:1900, 2022; PMC8993824; DOI 10.1038/s41467-022-29524-w). Code github.com/pFindStudio/GlycoDecipher (Apache-2.0) | Apache-2.0 (ref only) | `crossspectrum.rs` |
| Decoy-glycan (shifted intermediate Y-rungs) for 2D-FDR | a glyco search engine2/an open-search PTM tool decoy-glycan recipe; O-Pair/an open-source glyco engine (Lu et al., *Nat Commun* 11:4271, 2020; DOI 10.1038/s41467-020-18096-2; MIT) | Apache-2.0 / MIT (ref only) | `backbone.rs:498`, `glyco_pin.rs` |
| Oxonium gating | standard practice (every mature glyco engine) | — | `oxonium.rs` |

No a commercial glyco engine (commercial) or the reference engine (UM-proprietary) code is referenced or copied.
