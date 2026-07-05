# Cross-Spectrum Backbone Transfer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add opt-in cross-spectrum N-glycopeptide backbone transfer that propagates Pass-1 1%-FDR backbones to co-eluting, glycan-delta-linked sibling spectra, raising IDs @1% FDR above the 268 deterministic baseline at honest FDR.

**Architecture:** A two-pass flow inside a single glyco run. Pass 1 is today's glyco search + rescorer. Its 1%-FDR PSMs (target + decoy, symmetric) seed two glycan-delta graphs; each seed backbone is offered to co-eluting acceptor spectra and "sticks" where `precursor − backbone` is a known glycan. Stuck backbones enter Pass-2 scoring as `Source::Transferred` candidates with 5 additive PIN features. Both Percolator and the native GBDT consume the enriched PIN unchanged.

**Tech Stack:** Rust (workspace crates `andes-glyco`, `output`, `search`, `andes`), `cargo test`, Percolator 3.7.1 (biocontainer), Python recovery scripts on the benchmark VM.

## Global Constraints

- **FDR authority = Percolator** (production); native GBDT rescorer is a supported fallback. Feature is rescorer-agnostic.
- **Additive PIN features only** — never modify or reorder existing PIN columns/features.
- **Clean-room** — no borrowed a glyco search engine/a commercial glyco engine/O-Pair/a cross-spectrum glyco engine code.
- **Deterministic** — same input ⇒ byte-identical output. **No `HashMap`/`HashSet` in any output-bearing path.** All sorts carry a total-order tiebreak.
- **Gated** — entire feature behind `--glyco-transfer` (default `false`). With the flag off, output must be byte-identical to the current baseline.
- **Baseline to beat:** 268 PSMs @1% FDR / 96 backbone-correct on PXD025455 Fc3_r1, decoys @1% controlled (~1), deterministic.
- **TDD** — every task: write failing test → watch it fail → minimal code → pass → commit.

Reference design: `docs/plans/glyco/50-roadmap/cross-spectrum-transfer-design.md`.

---

## File Structure

- `crates/andes-glyco/src/hybrid.rs` — add `Source::Transferred` variant.
- `crates/andes-glyco/src/crossspectrum.rs` — extend: `GlycoNode`, `Seed`, `TransferredCandidate`, `build_glyco_nodes`, `propagate_transfers` (builds on existing `GlycoformWhitelist`/`nearest_glycan`).
- `crates/andes-glyco/src/glyco_psm.rs` — add 5 transfer fields to `GlycoPsmKey` (default inert).
- `crates/output/src/glyco_pin.rs` — emit 5 new PIN columns after `SialicConsistency`.
- `crates/andes/src/glyco_seeds.rs` (new) — parse a rescored PIN → 1%-FDR `Seed`s (target + decoy).
- `crates/andes/src/bin/andes.rs` — `--glyco-transfer` flag + two-pass orchestration around the existing glyco block (~line 2349).

---

## Task 1: `Source::Transferred` provenance variant

**Files:**
- Modify: `crates/andes-glyco/src/hybrid.rs:20-28`
- Test: same file `#[cfg(test)]`

**Interfaces:**
- Produces: `Source::Transferred` (enum variant, used by Tasks 4–6).

- [ ] **Step 1: Write the failing test**

```rust
// in crates/andes-glyco/src/hybrid.rs tests module
#[test]
fn source_has_transferred_variant_distinct_from_db_and_denovo() {
    assert_ne!(Source::Transferred, Source::Db);
    assert_ne!(Source::Transferred, Source::DeNovo);
    // clone + eq hold (derives intact)
    assert_eq!(Source::Transferred, Source::Transferred.clone());
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p andes-glyco source_has_transferred_variant -- --nocapture`
Expected: FAIL — `no variant named Transferred found for enum Source`.

- [ ] **Step 3: Add the variant**

```rust
pub enum Source {
    /// Backbone computed as precursor_neutral − known glycan mass.
    Db,
    /// Backbone proposed by the de-novo Y-ladder solver.
    DeNovo,
    /// Backbone borrowed from a confident co-eluting sibling spectrum
    /// (cross-spectrum transfer). Carries no per-spectrum core-Y anchor.
    Transferred,
}
```

- [ ] **Step 4: Run tests**

Run: `cargo test -p andes-glyco` — Expected: PASS (new test + all existing; any `match Source` without a wildcard now needs the arm — fix by adding `Source::Transferred => …` mirroring `Source::Db` where the compiler flags it).

- [ ] **Step 5: Commit**

```bash
git add crates/andes-glyco/src/hybrid.rs
git commit -m "feat(glyco): add Source::Transferred provenance variant"
```

---

## Task 2: Graph node + seed + transferred-candidate types

**Files:**
- Modify: `crates/andes-glyco/src/crossspectrum.rs` (append types near top, after `use`)
- Test: same file `#[cfg(test)]`

**Interfaces:**
- Produces:
  - `pub struct GlycoNode { pub scan: u32, pub precursor_neutral: f64, pub rt_seconds: Option<f64> }`
  - `pub struct Seed { pub scan: u32, pub peptide_idx: u32, pub backbone_mass: f64, pub rt_seconds: Option<f64>, pub seed_score: f64, pub is_decoy: bool }`
  - `pub struct TransferredCandidate { pub acceptor_scan: u32, pub peptide_idx: u32, pub backbone_mass: f64, pub glycan: GlycanComp, pub graph_support: u32, pub seed_score: f64, pub rt_delta: f64, pub ungated: bool, pub is_decoy: bool }`

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn transfer_types_construct_and_clone() {
    let g = crate::glycan_db::n_glycan_list()[0].clone();
    let n = GlycoNode { scan: 5, precursor_neutral: 2000.0, rt_seconds: Some(900.0) };
    let s = Seed { scan: 5, peptide_idx: 3, backbone_mass: 1500.0, rt_seconds: Some(900.0), seed_score: 2.5, is_decoy: false };
    let t = TransferredCandidate { acceptor_scan: 7, peptide_idx: 3, backbone_mass: 1500.0,
        glycan: g, graph_support: 4, seed_score: 2.5, rt_delta: 12.0, ungated: false, is_decoy: false };
    assert_eq!(n.clone().scan, 5);
    assert_eq!(s.clone().peptide_idx, 3);
    assert_eq!(t.clone().graph_support, 4);
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p andes-glyco transfer_types_construct` — Expected: FAIL (`cannot find struct GlycoNode`).

- [ ] **Step 3: Add the types**

```rust
/// A glyco-candidate spectrum as a graph node.
#[derive(Debug, Clone)]
pub struct GlycoNode {
    pub scan: u32,
    pub precursor_neutral: f64,
    pub rt_seconds: Option<f64>,
}

/// A confident Pass-1 seed backbone to propagate. `peptide_idx` indexes the
/// driver's candidate array so Pass-2 can re-score the exact peptide+mods.
#[derive(Debug, Clone)]
pub struct Seed {
    pub scan: u32,
    pub peptide_idx: u32,
    pub backbone_mass: f64,
    pub rt_seconds: Option<f64>,
    pub seed_score: f64,
    pub is_decoy: bool,
}

/// A backbone transferred onto an acceptor spectrum, plus its graph evidence.
#[derive(Debug, Clone)]
pub struct TransferredCandidate {
    pub acceptor_scan: u32,
    pub peptide_idx: u32,
    pub backbone_mass: f64,
    pub glycan: GlycanComp,
    pub graph_support: u32,
    pub seed_score: f64,
    pub rt_delta: f64,
    pub ungated: bool,
    pub is_decoy: bool,
}
```

- [ ] **Step 4: Run tests** — Run: `cargo test -p andes-glyco` — Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add crates/andes-glyco/src/crossspectrum.rs
git commit -m "feat(glyco): add GlycoNode/Seed/TransferredCandidate transfer types"
```

---

## Task 3: `propagate_transfers` — the graph propagation core

**Files:**
- Modify: `crates/andes-glyco/src/crossspectrum.rs`
- Test: same file `#[cfg(test)]`

**Interfaces:**
- Consumes: `Seed`, `GlycoNode`, `TransferredCandidate` (Task 2); `nearest_glycan`, `GlycanComp`.
- Produces:
```rust
pub fn propagate_transfers(
    seeds: &[Seed],
    nodes: &[GlycoNode],
    glycan_sorted: &[(f64, usize)],
    glycans: &[GlycanComp],
    rt_window: f32,
    min_glycan: f64,
    tol: f64,
) -> Vec<TransferredCandidate>
```

**Semantics:** For each seed, offer its backbone to every node co-eluting within `rt_window` (a node with `rt_seconds = None`, or a seed with none, is treated as ungated: it still matches but sets `ungated = true`). A node accepts if `precursor_neutral − backbone` is within `tol` of a known glycan and `≥ min_glycan`. Do not emit a transfer onto the seed's own scan. `graph_support` for an emitted candidate = number of *distinct accepting nodes* for that seed (the glycoform-family size). `rt_delta` = `|acceptor_rt − seed_rt|` (0.0 if either is `None`). `is_decoy` inherited from the seed. Output sorted by `(acceptor_scan, backbone_mass, glycan.mass)` with a final total-order tiebreak — deterministic, no HashMap.

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn propagate_transfers_recovers_sibling_and_counts_support() {
    let glycans = n_glycan_list();
    let sorted = sorted_view(&glycans);
    let backbone = 1500.0_f64;
    let g_a = 2.0 * HEXNAC + 3.0 * HEX;      // core
    let g_b = 2.0 * HEXNAC + 4.0 * HEX;      // +1 Hex sibling
    // seed on scan 1 (a well-fragmented glycoform), two co-eluting acceptors.
    let seeds = vec![Seed { scan: 1, peptide_idx: 9, backbone_mass: backbone,
        rt_seconds: Some(900.0), seed_score: 3.0, is_decoy: false }];
    let nodes = vec![
        GlycoNode { scan: 1, precursor_neutral: backbone + g_a, rt_seconds: Some(900.0) }, // self: skip
        GlycoNode { scan: 2, precursor_neutral: backbone + g_a, rt_seconds: Some(905.0) }, // sibling
        GlycoNode { scan: 3, precursor_neutral: backbone + g_b, rt_seconds: Some(910.0) }, // sibling
    ];
    let out = propagate_transfers(&seeds, &nodes, &sorted, &glycans, 300.0, 406.0, 0.05);
    // Two transfers (scan 2 and 3), NOT onto the seed's own scan 1.
    assert_eq!(out.len(), 2, "got {out:?}");
    assert!(out.iter().all(|t| t.acceptor_scan != 1));
    assert!(out.iter().all(|t| t.peptide_idx == 9 && (t.backbone_mass - backbone).abs() < 1e-6));
    // graph_support = family size (2 accepting non-self nodes).
    assert!(out.iter().all(|t| t.graph_support == 2), "support {out:?}");
}

#[test]
fn propagate_transfers_respects_rt_window_and_marks_ungated() {
    let glycans = n_glycan_list();
    let sorted = sorted_view(&glycans);
    let backbone = 1500.0_f64;
    let g = 2.0 * HEXNAC + 3.0 * HEX;
    let seeds = vec![Seed { scan: 1, peptide_idx: 0, backbone_mass: backbone,
        rt_seconds: Some(1000.0), seed_score: 1.0, is_decoy: false }];
    // one co-eluting (RT 1100), one far (RT 5000), one with no RT (ungated).
    let nodes = vec![
        GlycoNode { scan: 2, precursor_neutral: backbone + g, rt_seconds: Some(1100.0) },
        GlycoNode { scan: 3, precursor_neutral: backbone + g, rt_seconds: Some(5000.0) },
        GlycoNode { scan: 4, precursor_neutral: backbone + g, rt_seconds: None },
    ];
    let out = propagate_transfers(&seeds, &nodes, &sorted, &glycans, 300.0, 406.0, 0.05);
    let scans: Vec<u32> = out.iter().map(|t| t.acceptor_scan).collect();
    assert!(scans.contains(&2), "co-eluting must transfer: {out:?}");
    assert!(!scans.contains(&3), "far RT must NOT transfer: {out:?}");
    // ungated node (no RT) still receives, flagged ungated.
    let u = out.iter().find(|t| t.acceptor_scan == 4).expect("ungated transfer");
    assert!(u.ungated, "no-RT acceptor must be marked ungated");
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p andes-glyco propagate_transfers` — Expected: FAIL (`cannot find function propagate_transfers`).

- [ ] **Step 3: Implement**

```rust
/// Propagate confident seed backbones to co-eluting, glycan-consistent acceptor
/// nodes. See module docs. Deterministic: output sorted by
/// (acceptor_scan, backbone_mass, glycan mass) with an index tiebreak.
pub fn propagate_transfers(
    seeds: &[Seed],
    nodes: &[GlycoNode],
    glycan_sorted: &[(f64, usize)],
    glycans: &[GlycanComp],
    rt_window: f32,
    min_glycan: f64,
    tol: f64,
) -> Vec<TransferredCandidate> {
    let co_elutes = |seed_rt: Option<f64>, node_rt: Option<f64>| -> (bool, bool) {
        // returns (passes_gate, ungated)
        match (seed_rt, node_rt) {
            (Some(s), Some(n)) => (((n - s).abs() as f32) <= rt_window, false),
            _ => (true, true), // missing RT on either side: accept, flag ungated
        }
    };
    let mut out: Vec<TransferredCandidate> = Vec::new();
    for seed in seeds {
        // First pass over nodes: collect accepting (non-self) acceptors so we
        // can attach the family-size support to each emitted candidate.
        let mut accepted: Vec<(u32, GlycanComp, f64, bool)> = Vec::new();
        for node in nodes {
            if node.scan == seed.scan {
                continue;
            }
            let (gate, ungated) = co_elutes(seed.rt_seconds, node.rt_seconds);
            if !gate {
                continue;
            }
            let glycan_mass = node.precursor_neutral - seed.backbone_mass;
            if glycan_mass < min_glycan {
                continue;
            }
            if let Some(g) = nearest_glycan(glycan_sorted, glycans, glycan_mass, tol) {
                let rt_delta = match (seed.rt_seconds, node.rt_seconds) {
                    (Some(s), Some(n)) => (n - s).abs(),
                    _ => 0.0,
                };
                accepted.push((node.scan, g, rt_delta, ungated));
            }
        }
        let support = accepted.len() as u32;
        for (acceptor_scan, glycan, rt_delta, ungated) in accepted {
            out.push(TransferredCandidate {
                acceptor_scan,
                peptide_idx: seed.peptide_idx,
                backbone_mass: seed.backbone_mass,
                glycan,
                graph_support: support,
                seed_score: seed.seed_score,
                rt_delta,
                ungated,
                is_decoy: seed.is_decoy,
            });
        }
    }
    // Deterministic total order (no HashMap anywhere above).
    out.sort_by(|a, b| {
        a.acceptor_scan
            .cmp(&b.acceptor_scan)
            .then(a.backbone_mass.total_cmp(&b.backbone_mass))
            .then(a.glycan.mass.total_cmp(&b.glycan.mass))
            .then(a.peptide_idx.cmp(&b.peptide_idx))
    });
    out
}
```

- [ ] **Step 4: Run tests** — Run: `cargo test -p andes-glyco propagate_transfers` — Expected: PASS (both tests).

- [ ] **Step 5: Commit**

```bash
git add crates/andes-glyco/src/crossspectrum.rs
git commit -m "feat(glyco): propagate_transfers — glycan-delta graph propagation core"
```

---

## Task 4: Determinism guard for propagation

**Files:**
- Modify: `crates/andes-glyco/src/crossspectrum.rs` (test only)

**Interfaces:** Consumes `propagate_transfers` (Task 3).

- [ ] **Step 1: Write the failing test** (fails only if a future edit reintroduces order-dependence)

```rust
#[test]
fn propagate_transfers_is_deterministic_across_input_orders() {
    let glycans = n_glycan_list();
    let sorted = sorted_view(&glycans);
    let bb = 1500.0_f64;
    let g = 2.0 * HEXNAC + 3.0 * HEX;
    let seeds = vec![
        Seed { scan: 1, peptide_idx: 0, backbone_mass: bb, rt_seconds: Some(900.0), seed_score: 1.0, is_decoy: false },
        Seed { scan: 9, peptide_idx: 1, backbone_mass: bb + 100.0, rt_seconds: Some(902.0), seed_score: 1.0, is_decoy: false },
    ];
    let mk = |order: &[usize]| {
        let base = vec![
            GlycoNode { scan: 2, precursor_neutral: bb + g, rt_seconds: Some(901.0) },
            GlycoNode { scan: 3, precursor_neutral: bb + g, rt_seconds: Some(903.0) },
            GlycoNode { scan: 4, precursor_neutral: bb + 100.0 + g, rt_seconds: Some(904.0) },
        ];
        let nodes: Vec<GlycoNode> = order.iter().map(|&i| base[i].clone()).collect();
        propagate_transfers(&seeds, &nodes, &sorted, &glycans, 300.0, 406.0, 0.05)
    };
    let a = mk(&[0, 1, 2]);
    let b = mk(&[2, 0, 1]);
    let c = mk(&[1, 2, 0]);
    let key = |v: &[TransferredCandidate]| -> Vec<(u32, u64, u64)> {
        v.iter().map(|t| (t.acceptor_scan, t.backbone_mass.to_bits(), t.glycan.mass.to_bits())).collect()
    };
    assert_eq!(key(&a), key(&b));
    assert_eq!(key(&a), key(&c));
}
```

- [ ] **Step 2: Run to verify it passes now** (Task 3 already sorts): Run: `cargo test -p andes-glyco propagate_transfers_is_deterministic` — Expected: PASS. (Sanity: temporarily delete the `out.sort_by(...)` block → rerun → observe FAIL, then restore.)

- [ ] **Step 3: (no impl needed — guard test)**

- [ ] **Step 4: Full suite** — Run: `cargo test -p andes-glyco` — Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add crates/andes-glyco/src/crossspectrum.rs
git commit -m "test(glyco): determinism guard for propagate_transfers"
```

---

## Task 5: Transfer fields on `GlycoPsmKey` (inert by default)

**Files:**
- Modify: `crates/andes-glyco/src/glyco_psm.rs` (struct + every constructor/test that builds `GlycoPsmKey`)
- Test: same file

**Interfaces:**
- Produces: `GlycoPsmKey` gains `is_transferred: bool`, `transfer_graph_support: u32`, `transfer_seed_score: f32`, `transfer_rt_delta: f32`, `transfer_ungated: bool`. Defaults (`false/0/0.0`) reproduce current behavior.

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn glyco_psm_key_defaults_to_non_transferred() {
    let key = GlycoPsmKey {
        spectrum_idx: 0, glycan: None, glycan_source: Source::Db,
        oxonium_summed_frac: 0.0, n_core_oxonium_ions: 0,
        y_ladder_intensity_score: 0.0, y_ladder_decoy_score: 0.0,
        y0y1_anchor_score: 0.0, sialic_consistency: 0.0, core_y_hits: 0,
        glycan_mass: 0.0, backbone_mass: 0.0,
        is_transferred: false, transfer_graph_support: 0,
        transfer_seed_score: 0.0, transfer_rt_delta: 0.0, transfer_ungated: false,
    };
    assert!(!key.is_transferred);
    assert_eq!(key.transfer_graph_support, 0);
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p andes-glyco glyco_psm_key_defaults_to_non_transferred` — Expected: FAIL (`missing field is_transferred` at other construction sites → compile error).

- [ ] **Step 3: Add fields + fix all constructors**

Append to the `GlycoPsmKey` struct (after `backbone_mass`):

```rust
    /// Cross-spectrum transfer provenance + evidence (additive PIN features).
    /// All inert (false/0) for natively-generated candidates.
    pub is_transferred: bool,
    /// # co-eluting, glycan-delta-linked sibling spectra corroborating this
    /// backbone (the discriminative transfer signal).
    pub transfer_graph_support: u32,
    /// Pass-1 discriminant of the donor seed.
    pub transfer_seed_score: f32,
    /// |RT(acceptor) − RT(seed)| seconds; 0 = perfect co-elution.
    pub transfer_rt_delta: f32,
    /// RT unavailable ⇒ co-elution gate skipped (distrust signal).
    pub transfer_ungated: bool,
```

Then update EVERY other `GlycoPsmKey { … }` literal (in `glyco_psm.rs` doctests/tests, `glyco_search.rs`, `glyco_pin.rs` tests) to append the five defaults:
```rust
    is_transferred: false, transfer_graph_support: 0,
    transfer_seed_score: 0.0, transfer_rt_delta: 0.0, transfer_ungated: false,
```
Find them with: `grep -rn "GlycoPsmKey {" crates/`.

- [ ] **Step 4: Run tests** — Run: `cargo test -p andes-glyco -p search -p output` — Expected: PASS (all crates compile + green).

- [ ] **Step 5: Commit**

```bash
git add crates/andes-glyco/src/glyco_psm.rs crates/search/src/glyco_search.rs crates/output/src/glyco_pin.rs
git commit -m "feat(glyco): add inert transfer fields to GlycoPsmKey"
```

---

## Task 6: Emit the 5 transfer PIN columns (additive)

**Files:**
- Modify: `crates/output/src/glyco_pin.rs` (header list ~line 111-114 and the row-emission fn)
- Test: same file

**Interfaces:** Consumes `GlycoPsmKey` transfer fields (Task 5). Columns appended AFTER `SialicConsistency`, BEFORE `Peptide` — additive, existing columns unchanged.

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn glyco_pin_header_has_transfer_columns_after_sialic() {
    let header = glyco_pin_header(); // existing header-builder; adjust name to the real one
    for col in ["IsTransferred", "TransferGraphSupport", "TransferSeedScore", "TransferRTDelta", "TransferUngated"] {
        assert!(header.contains(col), "header missing {col}");
    }
    // additive placement: all appear after SialicConsistency, before Peptide
    let pos = |c: &str| header.iter().position(|h| h == c).unwrap();
    assert!(pos("SialicConsistency") < pos("IsTransferred"));
    assert!(pos("TransferUngated") < pos("Peptide"));
}
```
(If the header is built inline rather than via a named fn, assert against the emitted header string from a `write_glyco_pin_to` call on a one-row fixture instead.)

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p output glyco_pin_header_has_transfer_columns` — Expected: FAIL (columns absent).

- [ ] **Step 3: Add columns to header + row**

In the header vec, after `"SialicConsistency".to_string(),`:
```rust
        "IsTransferred".to_string(),          // cross-spectrum transfer provenance (additive)
        "TransferGraphSupport".to_string(),   // # corroborating co-eluting siblings
        "TransferSeedScore".to_string(),      // donor seed Pass-1 discriminant
        "TransferRTDelta".to_string(),        // |RT(acceptor)-RT(seed)| seconds
        "TransferUngated".to_string(),        // 1 = RT gate skipped (no RT)
```
In the row-value emission (mirror where `SialicConsistency` value is written from `key`):
```rust
        fmt_bool(key.is_transferred),               // IsTransferred
        (key.transfer_graph_support as f64),        // TransferGraphSupport
        (key.transfer_seed_score as f64),           // TransferSeedScore
        (key.transfer_rt_delta as f64),             // TransferRTDelta
        fmt_bool(key.transfer_ungated),             // TransferUngated
```
(Use the module's existing numeric/bool feature formatting; match the surrounding code's `Feature`/`format_feature_value` convention. `fmt_bool` = the existing 1.0/0.0 idiom used by `IsGlycanDb`.)

- [ ] **Step 4: Run tests** — Run: `cargo test -p output` — Expected: PASS. Also `cargo test -p output glyco_pin` to confirm existing header tests (OxoniumScore/CoreYHits) still pass — proves additivity.

- [ ] **Step 5: Commit**

```bash
git add crates/output/src/glyco_pin.rs
git commit -m "feat(glyco): emit 5 additive transfer PIN columns"
```

---

## Task 7: Seed extraction from a rescored PIN

**Files:**
- Create: `crates/andes/src/glyco_seeds.rs`
- Modify: `crates/andes/src/lib.rs` (or `main`/module tree) to `mod glyco_seeds;`
- Test: `crates/andes/src/glyco_seeds.rs` `#[cfg(test)]`

**Interfaces:**
- Consumes: rescored PIN/PSM output (reuse the PIN parsing style in `crates/andes/src/rescore.rs:41+` — `SpecId/PSMId`, `Label`, `q-value`/score columns) and a `scan → (peptide_idx, backbone_mass, rt_seconds, score)` map from the Pass-1 glyco results.
- Produces:
```rust
pub struct SeedRow { pub scan: u32, pub is_decoy: bool, pub q_value: f64, pub score: f64 }
pub fn parse_seed_rows(psms_tsv: &str) -> Result<Vec<SeedRow>, String>;
pub fn seeds_at_fdr(rows: &[SeedRow], q_threshold: f64,
    lookup: impl Fn(u32) -> Option<(u32, f64, Option<f64>)>) -> Vec<andes_glyco::crossspectrum::Seed>;
```
`lookup(scan) -> (peptide_idx, backbone_mass, rt_seconds)` bridges the PIN row back to the Pass-1 winner for that scan.

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn seeds_at_fdr_keeps_only_q_below_threshold_and_maps_backbone() {
    let rows = vec![
        SeedRow { scan: 10, is_decoy: false, q_value: 0.004, score: 3.1 },
        SeedRow { scan: 11, is_decoy: false, q_value: 0.05,  score: 0.4 },  // fails FDR
        SeedRow { scan: 12, is_decoy: true,  q_value: 0.008, score: 2.0 },  // decoy seed
    ];
    let lookup = |scan: u32| match scan {
        10 => Some((100u32, 1500.0f64, Some(900.0f64))),
        12 => Some((200u32, 1800.0f64, Some(905.0f64))),
        _ => None,
    };
    let seeds = seeds_at_fdr(&rows, 0.01, lookup);
    assert_eq!(seeds.len(), 2, "only q<=0.01 rows with a backbone: {seeds:?}");
    assert!(seeds.iter().any(|s| s.scan == 10 && !s.is_decoy && (s.backbone_mass - 1500.0).abs() < 1e-6));
    assert!(seeds.iter().any(|s| s.scan == 12 && s.is_decoy));
    assert!(!seeds.iter().any(|s| s.scan == 11));
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p andes seeds_at_fdr` — Expected: FAIL (module/fn missing).

- [ ] **Step 3: Implement**

```rust
//! Extract 1%-FDR Pass-1 glyco seeds (target + decoy) from a rescored PIN.
use andes_glyco::crossspectrum::Seed;

#[derive(Debug, Clone)]
pub struct SeedRow { pub scan: u32, pub is_decoy: bool, pub q_value: f64, pub score: f64 }

/// Parse a Percolator/native `.psms`+`.dpsms` (or combined) TSV. Columns:
/// PSMId/SpecId, Label (1 target / -1 or 0 decoy), score, q-value. Scan is the
/// integer ScanNr embedded in the SpecId (…scan=<N>…), matching glyco PIN ids.
pub fn parse_seed_rows(psms_tsv: &str) -> Result<Vec<SeedRow>, String> {
    let mut lines = psms_tsv.lines();
    let header = lines.next().ok_or("empty psms")?;
    let cols: Vec<&str> = header.split('\t').collect();
    let find = |name: &str| cols.iter().position(|c| c.eq_ignore_ascii_case(name));
    let id_i = find("PSMId").or_else(|| find("SpecId")).ok_or("no PSMId/SpecId")?;
    let q_i = find("q-value").ok_or("no q-value")?;
    let s_i = find("score").unwrap_or(q_i);
    let l_i = find("Label");
    let mut out = Vec::new();
    for line in lines {
        if line.trim().is_empty() { continue; }
        let f: Vec<&str> = line.split('\t').collect();
        let id = f.get(id_i).copied().unwrap_or("");
        let scan = extract_scan(id).ok_or_else(|| format!("no scan in id {id}"))?;
        let q_value: f64 = f.get(q_i).and_then(|v| v.parse().ok()).unwrap_or(1.0);
        let score: f64 = f.get(s_i).and_then(|v| v.parse().ok()).unwrap_or(0.0);
        let is_decoy = l_i
            .and_then(|i| f.get(i))
            .map(|v| v.trim().starts_with('-') || v.trim() == "0")
            .unwrap_or_else(|| id.contains("DECOY_"));
        out.push(SeedRow { scan, is_decoy, q_value, score });
    }
    Ok(out)
}

/// Pull the integer scan from a glyco SpecId like
/// "controllerType=0 controllerNumber=1 scan=3000_glyco_3000_1".
fn extract_scan(id: &str) -> Option<u32> {
    let after = id.split("scan=").nth(1)?;
    let digits: String = after.chars().take_while(|c| c.is_ascii_digit()).collect();
    digits.parse().ok()
}

/// Keep rows at `q <= q_threshold` that map to a Pass-1 backbone, as `Seed`s.
pub fn seeds_at_fdr(
    rows: &[SeedRow],
    q_threshold: f64,
    lookup: impl Fn(u32) -> Option<(u32, f64, Option<f64>)>,
) -> Vec<Seed> {
    let mut seeds: Vec<Seed> = rows
        .iter()
        .filter(|r| r.q_value <= q_threshold)
        .filter_map(|r| {
            let (peptide_idx, backbone_mass, rt_seconds) = lookup(r.scan)?;
            Some(Seed { scan: r.scan, peptide_idx, backbone_mass, rt_seconds,
                seed_score: r.score, is_decoy: r.is_decoy })
        })
        .collect();
    // Deterministic order (target/decoy seeds propagate identically downstream).
    seeds.sort_by(|a, b| a.scan.cmp(&b.scan).then(a.backbone_mass.total_cmp(&b.backbone_mass)));
    seeds
}
```
Add `mod glyco_seeds;` (and `pub use` if needed) to the andes crate module tree.

- [ ] **Step 4: Run tests** — Run: `cargo test -p andes seeds_at_fdr` and `cargo test -p andes` — Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add crates/andes/src/glyco_seeds.rs crates/andes/src/lib.rs
git commit -m "feat(glyco): seed extraction from rescored PIN (1% FDR, target+decoy)"
```

---

## Task 8 — REVISED (2026-07-06, after driver recon)

Recon changed Task 8's shape. Facts: (1) `glyco_search_run` (search crate) already has a Pass-1/Pass-2 transfer using `core_y_hits>=3` donors + flat `GlycoformWhitelist`, gated `ANDES_GLYCO_CROSSSPECTRUM=1` (default off; the 253/97 baseline is transfer-free) — this legacy path is SUPERSEDED. (2) The native GBDT rescorer (`andes::rescore::native_rescore_pin`) lives ABOVE `search` in the crate graph, so seeding must be orchestrated in the DRIVER. (3) `native_rescore_pin` is target-only. (4) The per-candidate scoring is a 400-line closure `process_one` inside `glyco_search_run`; provenance must ride on `BackboneHit`.

Chosen architecture (user: single-invocation, in-process native GBDT): one `andes --glyco --glyco-transfer` run does Pass-1 → build PIN text → native-GBDT rescore (in-proc) → 1%-FDR seeds (target+decoy) → build nodes → `propagate_transfers` → convert to provenance-bearing `BackboneHit`s → Pass-2 injection → final PIN → external Percolator for production FDR. Split into 4 sub-tasks:

### Task 8a: `BackboneHit` provenance fields + wire into the key
**Files:** Modify `crates/andes-glyco/src/hybrid.rs` (BackboneHit struct + all constructors), `crates/search/src/glyco_search.rs` (key build at ~929-933 reads `bb_hit.transfer_*`).
**Produces:** `BackboneHit` gains `is_transferred: bool, transfer_graph_support: u32, transfer_seed_score: f32, transfer_rt_delta: f32, transfer_ungated: bool` (inert defaults). At the `GlycoPsmKey` construction (glyco_search.rs ~929-933), replace the hardcoded inert values with `bb_hit.is_transferred` / `bb_hit.transfer_graph_support` / etc. so a transferred backbone's provenance reaches the PIN.
**Test:** unit test — a `BackboneHit` with `is_transferred=true, transfer_graph_support=5` scored through the key path yields a `GlycoPsmKey` with those values (or a focused test that the key copies them). TDD.

### Task 8b: expose target+decoy q-values from the native rescorer
**Files:** Modify `crates/andes/src/rescore.rs`.
**Produces:** `pub fn native_rescore_qvalues(pin_text: &str, seed: u64) -> Result<Vec<(String /*spec_id*/, bool /*is_decoy*/, f64 /*q*/, f64 /*score*/)>, String>` — same `parse_pin`→`cv_scores`→`qvalues` pipeline as `native_rescore_pin` but returns ALL rows (targets AND decoys), so the symmetric decoy graph can seed from decoy 1%-FDR PSMs. Keep `native_rescore_pin` unchanged.
**Test:** synthetic PIN (a few target + decoy rows with separable scores) → returns all rows, decoys included, q-values monotone non-decreasing by score rank. TDD.

### Task 8c: `glyco_transfer_pass2` entry point (extract `process_one`)
**Files:** Modify `crates/search/src/glyco_search.rs`.
**Produces:** Refactor the `process_one` closure into a standalone `fn` taking its captured context explicitly (frag_index, candidates, glycan_sorted, glycan_list, params, tol_ppm, top_k, flags), so it can be called outside `glyco_search_run`. Add `pub fn glyco_transfer_pass2(spectra, prepared, glycan_list, tol_ppm, backbone_top_k, pass1: Vec<GlycoSpectrumResult>, injected: &std::collections::BTreeMap<usize, Vec<BackboneHit>>) -> Vec<GlycoSpectrumResult>` that, for each spectrum with injected transferred backbones, re-runs scoring with them and supersedes its pass-1 entry (deterministic merge, BTreeMap not HashMap). `glyco_search_run`'s own Pass-1 call uses the same extracted fn (no behavior change — guard with the existing test suite).
**Test:** the existing glyco_search tests must stay green (proves the extraction preserved behavior); add one test that `glyco_transfer_pass2` with an injected transferred `BackboneHit` on a fixture spectrum emits a hit with `is_transferred=true`.

### Task 8d: driver `--glyco-transfer` orchestration
**Files:** Modify `crates/andes/src/bin/andes.rs` (CLI + glyco block ~2351-2392), `crates/andes/src/glyco_seeds.rs` (a `build_seed_lookup` helper if needed).
**Consumes:** 8a/8b/8c + `glyco_seeds::seeds_at_fdr`, `andes_glyco::crossspectrum::{GlycoNode, propagate_transfers}`.
**Produces:** when `cli.glyco && cli.glyco_transfer`: (1) run `glyco_search_run` (Pass-1, legacy xspec OFF); (2) write the Pass-1 glyco PIN to an in-memory `Vec<u8>` via `write_glyco_pin_to`; (3) `native_rescore_qvalues` on it; (4) build `SeedRow`s from the q-values (scan via the typed `ScanNr` column preferred — see carry-forward — else `extract_scan`; **fail loud on ambiguous decoy labels**), and a `lookup(scan) -> (peptide_idx, backbone_mass, rt)` from the Pass-1 winners; (5) `seeds_at_fdr(rows, 0.01, lookup)`; (6) build `Vec<GlycoNode>` from all oxonium-positive spectra, **sorted by scan** (carry-forward: propagate_transfers needs sorted input); (7) `propagate_transfers(...)`; (8) group transferred candidates by acceptor spec_idx into a `BTreeMap<usize, Vec<BackboneHit>>` with provenance set; (9) `glyco_transfer_pass2(...)`; (10) `write_glyco_pin` the final results. Gate: `--glyco-transfer` default false; a `crates/andes/tests/glyco_transfer_gate.rs` byte-identity test (flag off == baseline).
**Test:** byte-identity gate (flag off), plus (VM/fixture-permitting) a transferred-row-present check; the real functional check is Task 9's A/B.

---

## Task 8 (ORIGINAL — superseded by 8a-8d above): `--glyco-transfer` flag + two-pass driver orchestration

**Files:**
- Modify: `crates/andes/src/bin/andes.rs` (CLI struct near line 391; glyco block near 2349-2386)
- Test: `crates/andes/tests/glyco_transfer_gate.rs` (new integration test)

**Interfaces:**
- Consumes: `parse_seed_rows`/`seeds_at_fdr` (Task 7), `propagate_transfers` (Task 3), `glyco_search_run`, `write_glyco_pin`.
- Produces: end-to-end behavior — flag off ⇒ unchanged; flag on ⇒ Pass-2 PIN includes `Source::Transferred` rows with populated transfer features.

**Orchestration (flag on):**
1. Run `glyco_search_run` (Pass 1) → write Pass-1 glyco PIN.
2. Rescore it (existing path: Percolator if available, else native GBDT) → `.psms`/`.dpsms`.
3. `parse_seed_rows` + `seeds_at_fdr(…, 0.01, lookup)` where `lookup` reads the Pass-1 winner per scan (peptide_idx, backbone_mass) and the spectrum RT.
4. Build `Vec<GlycoNode>` from all oxonium-positive spectra (`build_glyco_nodes`).
5. `propagate_transfers(seeds, nodes, …)` → transferred candidates.
6. Pass 2: for each transferred candidate, score its `(peptide_idx → peptide, glycan)` against the acceptor spectrum with the existing glyco scorer, build a `FullGlycoPsm` whose `GlycoPsmKey` has `glycan_source = Source::Transferred` and the transfer fields set; merge into that scan's hit list; re-run the shared `collapse_cmp` top-1 selection.
7. Write the final PIN.

- [ ] **Step 1: Write the failing test** (gate behavior on a tiny fixture)

```rust
// crates/andes/tests/glyco_transfer_gate.rs
// Uses the repo's small glyco fixture (see crates/*/tests for the existing
// mgf/mzML fixture used by glyco tests; reuse that path).
#[test]
fn transfer_flag_off_matches_baseline_pin_bytes() {
    let baseline = run_andes_glyco(&["--glyco"]);            // helper: returns PIN bytes
    let with_flag_off = run_andes_glyco(&["--glyco", "--glyco-transfer=false"]);
    assert_eq!(baseline, with_flag_off, "flag off must be byte-identical to baseline");
}

#[test]
fn transfer_flag_on_emits_transferred_rows() {
    let pin = String::from_utf8(run_andes_glyco(&["--glyco", "--glyco-transfer"])).unwrap();
    assert!(pin.lines().next().unwrap().contains("IsTransferred"));
    // at least one data row flags a transfer (1.0 in the IsTransferred column)
    let hdr: Vec<&str> = pin.lines().next().unwrap().split('\t').collect();
    let it = hdr.iter().position(|c| *c == "IsTransferred").unwrap();
    assert!(pin.lines().skip(1).any(|l| l.split('\t').nth(it) == Some("1")),
        "expected at least one transferred row on the fixture");
}
```
(If the fixture is too small to produce a transfer, mark `transfer_flag_on_emits_transferred_rows` `#[ignore]` with a comment pointing at the VM A/B in Task 9 as the real functional check, and keep the byte-identity gate test as the hard CI gate.)

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p andes --test glyco_transfer_gate` — Expected: FAIL (`--glyco-transfer` unknown arg).

- [ ] **Step 3: Add the flag + orchestration**

CLI (after the `glyco` field ~line 392):
```rust
    /// Enable cross-spectrum backbone transfer (two-pass; glyco mode only).
    /// Off by default — baseline output is unchanged.
    #[arg(long = "glyco-transfer", default_value_t = false)]
    glyco_transfer: bool,
```
Then wrap the existing glyco block: keep Pass-1 exactly as is; when `args.glyco && args.glyco_transfer`, run steps 2-7 above before `write_glyco_pin`. Keep every new sort total-ordered; no HashMap in the merge. Reuse the existing scorer entry point that `glyco_search_run` uses per `(peptide, glycan, spectrum)` — factor it into a small `score_one_glyco(peptide_idx, glycan, spectrum) -> FullGlycoPsm` helper if not already callable, so Pass-2 reuses identical scoring.

- [ ] **Step 4: Run tests**

Run: `cargo test -p andes --test glyco_transfer_gate` then `cargo build --release -p andes` — Expected: PASS + clean release build.

- [ ] **Step 5: Commit**

```bash
git add crates/andes/src/bin/andes.rs crates/andes/tests/glyco_transfer_gate.rs
git commit -m "feat(glyco): --glyco-transfer two-pass orchestration (off by default)"
```

---

## Task 9: End-to-end honest-FDR A/B on Fc3 (VM)

**Files:**
- Create: `docs/plans/glyco/50-roadmap/run_transfer_ab.sh` (VM harness, mirrors `run_cap_sweep.sh`)

**Interfaces:** Consumes the release binary with `--glyco-transfer`. Produces the pass/fail validation numbers.

- [ ] **Step 1: Write the A/B harness** (one build, two arms, both truths, decoy check)

```bash
#!/bin/bash
set -uo pipefail
cd /srv/data/msgf-bench/glyco_bench
BIN=/srv/data/msgf-bench/andes-src/target/release/andes
MZ=HCC_pool_Late_Fc3_r1.mzML; DB=9606-reviewed-contam-decoy.fasta
IMG=quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2; LOG=transfer_ab.log
( cd /srv/data/msgf-bench/andes-src && cargo build --release -p andes 2>/tmp/tr_build.log ) && echo BUILD_OK | tee $LOG
arm(){ local tag=$1; shift
  "$BIN" --spectrum $MZ --database $DB --decoy-strategy none --decoy-prefix DECOY_ --glyco "$@" \
    --output-pin andes_$tag.glyco.pin > andes_$tag.out 2>&1
  echo "$tag rows=$(($(wc -l <andes_$tag.glyco.pin)-1)) md5=$(md5sum andes_$tag.glyco.pin|cut -d' ' -f1)" | tee -a $LOG
  docker run --rm -v "$PWD":/data $IMG percolator --seed 42 --only-psms \
    --results-psms /data/$tag.psms --decoy-results-psms /data/$tag.dpsms /data/andes_$tag.glyco.pin >$tag.perc.log 2>&1
  local t=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="q-value")q=i}NR>1&&$q<=0.01{c++}END{print c+0}' $tag.psms)
  local d=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="q-value")q=i}NR>1&&$q<=0.01{c++}END{print c+0}' $tag.dpsms)
  echo "  $tag targets@1%=$t decoys@1%=$d" | tee -a $LOG
  echo "  vs523:" | tee -a $LOG; python3 glyco_recovery_fdr.py truth_nglycan_residue.tsv $tag.psms 0.01 0.02 | tee -a $LOG
  echo "  vs196:" | tee -a $LOG; python3 glyco_recovery_fdr.py truth_196.tsv $tag.psms 0.01 0.02 | tee -a $LOG
}
arm base                       # transfer OFF = 268 baseline (sanity)
arm xfer --glyco-transfer      # transfer ON
echo "=== TRANSFER_AB_DONE ===" | tee -a $LOG
```

- [ ] **Step 2: Run it on the VM** (background; ~70 min for two arms)

Run: `scp … run_transfer_ab.sh pride-linux-vm:… && ssh pride-linux-vm 'cd … && nohup ./run_transfer_ab.sh &'`

- [ ] **Step 3: Evaluate against the bar**

PASS requires ALL of:
- `base` reproduces 268 @1% (guards against regressions in the shared path).
- `xfer` targets@1% **> 268** AND backbone-correct(523) **≥ 96**.
- `xfer` **decoys@1% still ~1** (honest FDR — a jump in decoys@1% = the feature is inflating ⇒ FAIL, revisit the symmetric decoy graph).
- Determinism: run `xfer` twice, PIN md5 identical.

- [ ] **Step 4: Record the result** in `docs/plans/glyco/STATUS.md` (numbers + verdict) and memory.

- [ ] **Step 5: Commit**

```bash
git add docs/plans/glyco/50-roadmap/run_transfer_ab.sh docs/plans/glyco/STATUS.md
git commit -m "test(glyco): end-to-end transfer A/B harness + Fc3 result"
```

---

## Self-Review

**Spec coverage:** §3 two-pass flow → Tasks 7,8. §4 graph/propagation → Tasks 2,3; features → Tasks 5,6. §4 decoy symmetry → `is_decoy` threaded through Tasks 2,3,7 + Task 9 decoy check. §5 determinism → Tasks 3,4 + total-order sorts throughout. §6 integration points → Tasks 1,5,6,7,8. §7 validation → Task 9. §2 `Source::Transferred` → Task 1. Rescorer-agnostic → Task 8 reuses existing rescore path; features additive (Task 6) so native GBDT consumes them. Gating → Task 8 + byte-identity test. **No gaps.**

**Placeholder scan:** All code steps carry real code. Two soft spots flagged explicitly with fallbacks: the exact glyco header-builder name (Task 6 — assert against emitted string if no named fn) and the fixture's ability to produce a transfer (Task 8 — `#[ignore]` + defer to Task 9). These are honest environment unknowns, not vague instructions.

**Type consistency:** `Seed`/`GlycoNode`/`TransferredCandidate` fields defined in Task 2 are used verbatim in Tasks 3,4,7. `GlycoPsmKey` transfer fields (Task 5) match the PIN emission (Task 6) and column names match Task 9's `IsTransferred` check. `seeds_at_fdr` signature (Task 7) matches its Task 8 call. `propagate_transfers` signature identical across Tasks 3,4,8.
