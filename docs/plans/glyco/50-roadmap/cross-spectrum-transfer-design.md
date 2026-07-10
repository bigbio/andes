# Cross-Spectrum Backbone Transfer — Design

**Status:** implemented as an EXPERIMENTAL prototype (`--glyco-transfer`, off by
default). ⚠ **NOT FDR-VALID — do not use for reported identifications.** See
"Known soundness bugs" below. Implemented Tasks 1–8d (2026-07-05/06); A/B on
Fc3_r1 was net-neutral, but that result is **confounded** by the bugs below.
**Branch:** `glyco-phase1`
**Author:** brainstormed with the user, 2026-07-05

## ⚠ Known soundness bugs (Codex adversarial review, 2026-07-06) — MUST fix before revival

The A/B (net-neutral, decoys@1%=1) looked benign only because transfer barely
moved the top-1 winners; the mechanism that would keep it honest is broken.
Fix in this order before trusting any `--glyco-transfer` number:

1. **[CRITICAL] Transferred seeds lose target/decoy identity** (`andes.rs`, the
   `TransferredCandidate → BackboneHit` injection). `BackboneHit` carries neither
   `peptide_idx` nor `is_decoy`, so Pass-2 scores *any* mass-matching candidate
   and the emitted row's label comes from that candidate — a **decoy** seed can
   emit a **target**-labeled row. Target transfers inflate without matching decoy
   transfers ⇒ the symmetric-decoy graph is invalid ⇒ Percolator q-values are not
   honest. *Fix:* carry the seed `peptide_idx` + `is_decoy` through Pass-2; emit
   label-locked transferred rows; add an end-to-end test that a low-q decoy seed
   produces decoy-labeled transferred PIN rows.
2. **[HIGH] Scan-only q-join can seed the wrong spectrum** (`andes.rs:2407-2507`).
   Lookup keyed by `spec.scan.unwrap_or(0)`; duplicate scans (multi-file) or MGF
   without `SCANS` collapse last-wins. *Fix:* join on exact emitted `SpecId` or a
   unique `spec_idx`; fail loud on duplicate/missing scan ids in transfer mode.
3. **[HIGH] Dedup erases transferred candidates** (`glyco_search.rs:147-183`).
   Dedup prefers `Source::Db` and does not preserve `Transferred`; a transferred
   hit with the same backbone/glycan as a DB hit collapses into the DB hit,
   losing `IsTransferred` + the non-Db scoring exemption — silently defeating
   transfer on the weak-ladder spectra it targets (this is a prime suspect for
   the net-neutral result). *Fix:* preserve/merge `Transferred` provenance in
   dedup; regression test.
4. **[HIGH] Missing RT ⇒ all-run mass transfer** (`crossspectrum.rs:182-187`).
   `co_elutes` accepts when either RT is missing (only flags `ungated`, an output
   feature, not a gate). *Fix:* require RT on both ends by default; explicit
   unsafe opt-in with stricter support otherwise.
5. **[MEDIUM] Fixed 2500 Da tolerance biases edges by mass** (`andes.rs:2523`).
   *Fix:* per-acceptor `(precursor·ppm).max(0.02)` tolerance in
   `propagate_transfers`.

---

## 1. Purpose & context

In stepped-HCD N-glyco data a single peptide backbone appears as many glycoforms
(≈6–7 per peptide in serum). Well-fragmented glycoforms (strong trimannosyl-core
Y-ladder) are confidently identified; their poorly-fragmenting siblings are not,
because per-spectrum candidate generation has no core-Y ladder to anchor the
backbone. Single-spectrum engines (the reference glyco engine, an open-source glyco engine/O-Pair) leave
those siblings on the table.

**Cross-spectrum transfer** borrows a *confident* backbone — learned from the
well-fragmented glycoforms in a first pass — and offers it to co-eluting sibling
spectra whose `precursor − backbone` is a known glycan. This is the
a cross-spectrum glyco engine idea (published +33–178%), and it is a structural advantage
per-spectrum engines cannot replicate.

**Baseline this must beat.** As of 2026-07-05 the deterministic andes glyco
pipeline (yladder-default collapse + no candidate cap) already scores **268 PSMs
@1% FDR / 96 backbone-correct** on PXD025455 Fc3_r1 — *ahead of* an open-source glyco engine
(~222). So transfer is no longer about reaching parity; it must **extend the lead
and beat the reference engine**, at honest FDR.

### Non-negotiable constraints
- **FDR authority = Percolator** (production). andes's native GBDT rescorer
  (mokapot-style 3-fold target/decoy CV over the same PIN) is a supported
  fallback. The design is rescorer-agnostic (see §4).
- **Additive PIN features only** — never modify existing features.
- **Clean-room** — no borrowed a glyco search engine/a commercial glyco engine/O-Pair/a cross-spectrum glyco engine *code*;
  the algorithm is reimplemented from first principles.
- **Deterministic** — same input ⇒ byte-identical output (see §5).

## 2. Design decisions (locked)

| Decision | Choice | Rationale |
|---|---|---|
| v1 scope | Full glycan-delta graph | Structural differentiation for a decisive lead |
| Seed confidence | Pass-1 PSMs @ **1% FDR only** | Seeds are already FDR-controlled; no invented confidence |
| FDR honesty | **Symmetric decoy-transfer graph** | Decoys get identical transfer ⇒ feature is calibrated, not target-inflating |
| Edge width | **Any valid glycan-composition delta** + RT co-elution | Max reach; honesty rests on RT gate + decoy graph |

**Interaction note.** FDR-only seeds + any-glycan-delta edges make propagation
effectively *single-hop*: a seed backbone is offered to every co-eluting node and
"sticks" where `precursor − backbone` is a known glycan. The graph's payoff is
therefore not chained reach but a **discriminative feature** — how many
co-eluting, delta-linked glycoform siblings corroborate a backbone (a real
backbone shows a whole RT-aligned ladder; a coincidence does not).

## 3. Architecture & data flow

Two-pass flow inside a single glyco run, symmetric across target and decoy:

1. **Pass 1** — the current glyco search, unchanged. Run the rescorer once on the
   Pass-1 PIN (Percolator in production; the native GBDT when Percolator is
   unavailable — whichever the run already uses for FDR) → 1% FDR PSMs. Split
   into **target seeds** and **decoy seeds** (decoy-peptide IDs at their own 1%
   FDR). Each seed → `(peptide, backbone_mass, rt_seconds, seed_score)`.
2. **Build two glycan-delta graphs** — one over target seeds, one over decoy
   seeds — by the identical procedure.
3. **Propagate** — offer each seed backbone to co-eluting nodes; it sticks to a
   node when `precursor − backbone` is a known glycan. Record per stuck backbone
   the graph-support count and RT offset to the seed.
4. **Pass 2** — for each acceptor spectrum, add stuck backbones as candidates
   tagged `Source::Transferred` (decoy seeds → decoy candidates), scored by the
   *existing* glyco scorer, collapsed top-1-per-scan as today.
5. **One rescorer run** (Percolator prod / native GBDT fallback) over the
   combined PIN (native + transferred, target + decoy) with the new additive
   features. FDR stays honest because decoys got identical transfer treatment.

The only genuinely new component is the graph + propagation module. Scoring,
collapse, PIN emission, and FDR are reused. Gated behind `--glyco-transfer`
(off by default) so the 268 baseline is untouched until transfer is proven.

## 4. The graph, propagation & features

**Nodes.** Every oxonium-positive (glyco-candidate) spectrum, carrying
`(precursor_neutral, rt_seconds: Option<f64>, scan)`. A node with `rt_seconds =
None` can still *receive* a transfer but cannot be RT-gated — it falls back to a
whole-run window and is flagged `TransferUngated` so the rescorer can distrust it.

**Edges.** Sort nodes by precursor mass. An edge links two co-eluting nodes
(`|Δrt| ≤ rt_window`) whose precursor difference matches any known
glycan-composition delta within tolerance, found by a sorted two-pointer sweep.

**Propagation (single-hop from FDR seeds).** Each target seed offers
`(peptide, backbone)` to co-eluting nodes; it sticks where `precursor_node −
backbone` is a known glycan. Decoy seeds do the identical thing on the decoy
graph. A stuck backbone becomes a Pass-2 candidate.

**Additive PIN features (symmetric on decoys):**
- `IsTransferred` (0/1) — provenance.
- `TransferGraphSupport` — # co-eluting, delta-linked glycoform siblings
  corroborating this backbone. *The key discriminative signal.*
- `TransferSeedScore` — the donor seed's Pass-1 discriminant.
- `TransferRTDelta` — RT offset acceptor↔seed (0 = perfect co-elution).
- `TransferUngated` (0/1) — RT unavailable, co-elution gate skipped.

Features are zero/absent for native candidates, so both Percolator and the
native GBDT see a clean additive schema and *learn* transfer reliability from
this run's target/decoy separation — no retraining, no frozen schema. The native
GBDT folds by `ScanNr`; a transferred candidate carries its **acceptor's**
ScanNr, so it folds with its own spectrum. Any residual cross-spectrum leakage is
symmetric across target/decoy, so q-values stay honest.

## 5. Determinism

- Node list sorted by `(precursor, scan)`; edge sweep is a sorted two-pointer
  with total-order tiebreaks; seeds iterated in sorted order.
- **No `HashMap` in any output-bearing path** — same discipline as the
  `order_peptide_first` determinism fix.
- Test `transfer_is_deterministic`: run the whole graph twice on shuffled input,
  assert identical output.

## 6. Integration points

- `Source::Transferred` variant in `crates/andes-glyco/src/hybrid.rs`.
- New graph code in `crates/andes-glyco/src/crossspectrum.rs` (extends the
  existing, tested `GlycoformWhitelist`): `build_graph(nodes) → edges`,
  `propagate(seeds, graph) → Vec<TransferredCandidate>` with
  `(scan, peptide, backbone, graph_support, seed_score, rt_delta, ungated)`.
- Driver (`crates/search/src/glyco_search.rs` / andes bin): after Pass-1
  Percolator, collect 1% FDR target + decoy seeds, build both graphs, inject
  transferred candidates into Pass-2 scoring. Gated by `--glyco-transfer`.
- PIN writer (`crates/output/src/glyco_pin.rs`, feature list in `pin.rs`): emit
  the 5 new columns.

## 7. Validation plan (the bar this must clear)

1. **Unit** — graph build, edge deltas, propagation stick/no-stick, RT gate,
   decoy symmetry, determinism.
2. **A/B on Fc3** (deterministic, both rescorers) — `--glyco-transfer` on vs off,
   vs both truths. Must show **↑ PSMs @1% AND ↑ backbone-correct with decoys @1%
   still controlled (~1)**. Transfer that only lifts targets is **rejected**.
3. **Truth-anchored** — transferred IDs landing on 523/196-truth scans must be
   backbone-correct at **≥ the native rate**, proving transfer recovers *right*
   backbones, not just more IDs.
4. **Symmetric-decoy sanity** — transferred-decoy count must track
   transferred-target count in the null; large asymmetry ⇒ the graph is leaking
   and the feature would inflate.

**Success** = beats **268** on Fc3 at honest FDR (decoys @1% controlled), then
holds on a second glyco dataset (Fc5_r2 `.raw`, already staged on the VM).

## 8. Out of scope (YAGNI for v1)

- Multi-hop / transitive re-seeding (transferred nodes becoming new seeds) — the
  looser-seed risk we explicitly excluded; revisit only if single-hop proves out.
- Two-tier gold/silver seeds (sub-FDR ladder-strong seeds).
- Entrapment-set validation beyond the symmetric decoy graph — add later if the
  decoy null proves insufficient.
- Orthogonal fragmentation (EThcD c/z) — a different lever, for datasets with ETD
  scans (Fc3 is HCD-only).
