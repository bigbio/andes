# andes N-glyco search — clean implementation roadmap (G0 → G4)

> ⚠️ **CORRECTION (2026-07-03) — read [`../LESSONS.md`](../LESSONS.md) first.** Every
> recovery gate below that cites the **523-scan truth** or the **`154/523 → 83
> top-1 → 66 true`** chain is **VOID** — those numbers came from un-collapsed
> multi-row PINs (fixed in commit 7c269aab). Done correctly (top-1-per-scan), andes
> yields **0 glyco-PSMs @1% FDR**. All gates must be re-expressed as: **unique glyco
> scans/glycopeptides @ FDR, on top-1-collapsed PINs, against a valid FDR-controlled
> reference (the published the reference engine-nglycan result or multi-engine consensus),
> validated by entrapment.** Do not treat the pre-correction numbers as targets.


*Lead-architect roadmap. Authoritative sequencing document; engineers implement from this.
Written 2026-07-02, branch `glyco-phase1` (HEAD ~`35d31bb9`). Synthesizes
`00-context/`, `10-tools/`, `20-theory/`, `30-standards/`, `40-data/`.*

---

## 1. Problem statement

andes finds the glycopeptide **backbone** in ~90% of stepped-HCD N-glyco spectra
(generation is near-ceiling; `PHASE1_RESULT.md`, verified this session) but only
**identifies** 66 of 523 truth PSMs at 1% FDR on the frozen eval file
(PXD025455 `HCC_pool_Late_Fc3_r1`, own the reference engine re-search truth). The loss is
entirely downstream of generation: **154/523 recovered @1% FDR → 83 top-1 correct
(ranking loses 71) → 66 true (111 false pass)**. The two deficits are **ranking**
(peptide-axis discrimination) and **FDR** (a valid 2-dimensional glyco FDR). Both
stem from one physical fact: stepped-HCD deposits collisional energy into labile
glycosidic bonds, so peptide-specific backbone b/y is present in only **~11%** of
glyco spectra (a cross-spectrum glyco engine, Fang et al. *Nat Commun* 2022, PMC8990002; Riley &
Malaker *JPR* 2020, PMC7425838). andes's current glyco PIN features
(`OxoniumScore`, `YLadderScore`, `CoreYHits`, `GlycanMass`) are all
backbone/spectrum-level — identical for two peptides competing at the same
backbone-mass window — so they add **zero** peptide-axis target/decoy separation
(`SPA2_RESULT.md`; top-1 target:decoy ≈ 1.24:1 → Percolator returns 0 IDs @1% on
backbone-only RankScore).

## 2. Guiding thesis

> **Generation is solved. The work is ranking + FDR. Retraining is calibration
> infrastructure, not the ranking fix. The levers are the Y0/Y1 peptide-mass anchor
> feature and in-process cross-spectrum transfer.**

Four consequences that shape every phase:

1. **IDs ≈ Coverage × Separation**, and the deficit is entirely *Separation*. A
   better-calibrated model reshapes a likelihood over the same sparse peaks; it
   **cannot manufacture b/y ions the fragmentation never deposited**
   (`20-theory/why-andes-fails-and-succeed.md`).
2. **Y0 (bare peptide) and Y1 (peptide+HexNAc) are the one recoverable
   peptide-specific signal.** They encode the peptide mass *directly*, are
   high-intensity even when interior b/y is dead, and — unlike oxonium/YLadder —
   they **discriminate competing peptides**. This is the SP-B lever
   (`30-standards/masses.md` Y-ion convention; `40-data/collection/spectra_examples.md`
   takeaway 2).
3. **Single-spectrum scoring is bounded by the ~11% b/y stratum.** The only
   field-validated mechanism that breaks past ~30% is **RT-gated cross-spectrum
   transfer** — a cross-spectrum glyco engine's spectrum-expansion delivered +33.5–178.5%
   glyco-PSMs from *transfer, not scoring* (Fang 2022, PMC8990002). This is G4,
   andes's differentiator, already scaffolded but gated OFF.
4. **FDR is 2-dimensional** (peptide axis and glycan-Y axis fail independently). A
   unified Percolator pile crashed recovery **29.4% → 4.4%** this session
   (glycan-decoy rows differ only in YLadder → flood the −1 pile → Percolator
   over-weights YLadder). 2D-FDR **must** be a thin separate-axis Percolator
   post-process (`20-theory/glyco-fdr.md`).

**Sequencing rule (de-risk cheap-first):** the Y0/Y1 anchor + a decoy-separated
**kill-gate** are cheap and decide whether the expensive harvest is even worth it.
Front-load them **before** the multi-dataset harvest/retrain.

## 3. Phase table G0 → G4

| Phase | Goal | Key work | Gate (go/no-go) | Depends on |
|---|---|---|---|---|
| **G0 — Correctness hygiene** | Remove known defects so later A/Bs measure signal, not bugs. | DET-1: `hybrid.rs:318` sort → `total_cmp` (determinism). P0.3: Y0/Y1-only quorum retention (*measure before shipping*). P0.4: probe isotope fidelity. Confirm H2O convention pinned to `30-standards/masses.md` (residue masses, **no +18.0106** on attached glycan). | All 46 existing glyco tests green; find-rate A/B on 523 truth scans (`ANDES_GLYCO_SCANS`) is unchanged or up; runs bit-reproducible across two invocations. | none (do first) |
| **G1 — Glycan-Y-first generation** | Lock in near-ceiling backbone generation as the spine (the differentiator vs peptide-first the reference engine/a commercial glyco engine). | Already VERIFIED this session: peptide-independent Y-complement index + 2-core-Y-bit retention (`ANDES_GLYCO_YINDEX`), backbone-findability **59.3% → 69.8% @0.05 Da**, +7.2 pts even @0.005 Da. Promote from opt-in toward default after G0; keep peptide-first (`ANDES_GLYCO_PEPTIDE_FIRST`) as a *fallback branch only*. | Find-rate ≥ current on all truth scans with `YINDEX` on; no regression on eval when peptide-first is the fallback, not the spine. | G0 |
| **G2 — Y0/Y1 anchor + kill-gate** ★cheap decision point | Add the one peptide-specific PIN feature and *decide whether retraining can rank at all* before spending the harvest. | Add **additive-only** `Y0Y1AnchorScore` PIN column (peptide-mass-conditioned intensity of Y0=`M_pep+proton`, Y1=`M_pep+203.079373+proton`; parity-safe, never fused into the ranking score). Run the **decoy-separated kill-gate**: true peptide's backbone+anchor score vs a same-backbone reversed-peptide decoy's — measure separation, **not find-rate**. | **GO** if the anchor feature lifts top-1 (83→ up) *and* gives measurable target/decoy separation on the kill-gate. **NO-GO / skip to G4** if separation stays flat (the honest premise: retraining alone is likely insufficient). | G0 (G1 recommended) |
| **G3 — Harvest + regime-matched retrain (calibration)** | Give the peptide axis a learned, stepped-HCD-matched model under `protocol=NGlyco` — the calibration layer the anchor feature rides on. | (a) `Protocol::NGlyco` variant + CLI enum + `build_selection_key` arm → `experiment_class="glyco"` (see §7 code anchors). (b) **Glycan-stripped backbone training rows** — glycan mass must NOT shift backbone b/y; rows carry the bare peptide (the real modeling decision). (c) Truth-TSV→flat-parquet converter script. (d) Harvest on Codon per §5; **tiered labels**. Retrain the own `strong` spectral model; write `resources/models/protocol=NGlyco`. | Retrained `protocol=NGlyco` model + anchor beats the field-default model on eval top-1 **and** clears the G2 kill-gate on held-out data. Entrapment-FDP honest (foreign glycome+proteome padding, §6). | G2 = GO |
| **G4 — In-process cross-spectrum transfer** ★ceiling-breaker | Break past the ~11%/~30% single-spectrum ceiling — the andes-unique differentiator. | Activate the existing scaffold (`ANDES_GLYCO_CROSSSPECTRUM`, `crossspectrum.rs`, two-pass at `glyco_search.rs:790-887`): add **RT gating** + cosine-weighted transfer of a backbone fragment-frequency prior from confident donor PSMs to unassigned same-backbone spectra (a cross-spectrum glyco engine mechanism, Apache-2.0 reference). Pays only once G2/G3 lift donor-glycoform top-1 confidence. | Net glyco-PSM gain on eval (target: recover a meaningful fraction of the sialylated/short-peptide sparse-b/y stratum) with 2D-FDR (§below) held at 1%. | G3 (donor confidence) |
| **G3′ — 2D-FDR post-process** (runs *with* G3/G4, not after) | Valid glyco FDR without a unified pile or a second FDR engine. | Two vanilla Percolator runs + one inclusion-exclusion merge: **(1)** peptide-axis run — reversed-peptide decoys + Y0/Y1 anchor feature → `q_P`; **(2)** glycan-axis run — Y-rung-shifted glycan decoys (shift all peptide+Y ions **except Y0/Y1**, or per-fragment random 1–30 Da; recipe `20-theory/glyco-fdr.md`) with glycan-only features → `q_G`; **(3)** accept where `FDR_P + FDR_G − FDR_{P∩G} ≤ 0.01`. Glycan-decoy rows already emitted at `glyco_pin.rs:320-330`. | Combined 2D q ≤ 1% holds against entrapment-FDP; recovery does **not** collapse (contrast the refuted unified pile 29.4%→4.4%). | G2 (anchor feature); pairs with G3/G4 |

★ = the two levers the thesis rests on: **G2** (cheap kill-gate decision) and **G4**
(ceiling-breaker). Retraining (G3) is the calibration layer between them.

## 4. Standardized-mass + notation decisions (adopt as-is; cite `30-standards/`)

**Residue masses** (free monosaccharide − H2O; sum these, no extra water — the
glycosidic bond to Asn is itself a condensation). Cross-checked to 6 decimals vs
Unimod/ProForma2.0/a glyco search engine/FragPipe (`30-standards/masses.md`,
`40-data/glycan-db.md`):

| Residue | Mass (Da) | | Residue | Mass (Da) |
|---|---|---|---|---|
| Hex | 162.052824 | | NeuGc | 307.090331 |
| HexNAc | 203.079373 | | Pent | 132.042259 |
| Fuc / dHex | 146.057909 | | HexA | 176.032088 |
| NeuAc | 291.095417 | | Kdn | 250.068867 |
| +Phospho | 79.966331 | | +Sulfo | 79.956815 |

Constants: **proton 1.007276** (never neutral-H 1.00783 — else ~0.5 mDa/charge
drift), **water 18.010565**.

**Mass rules (pin to these):**
- Glycopeptide neutral mass = `M_peptide (Σresidues + one H2O + mods) +
  glycan_residue_sum`. Glycan is a **single delta on the sequon residue**.
- **Never add +18.010565 to an attached glycan** — that water is only for
  free/released glycans (GlycoMod frame). This is the #1 divergence risk
  (−18.0106 Da double-count). PHASE1's "H2O convention" fix pins here.
- Store full 6-decimal masses (203.079373, not 203.0793 — 7 mDa can flip
  near-isobaric HexNAc4 vs Hex5 +2 Da; ppmFixer *Glycobiology* 2024).
- **Y-ion convention:** Y0 neutral = `M_peptide`; Y1 = `M_peptide + 203.079373`;
  observed 1+ adds one proton (Y0 m/z = `M_peptide + 1.007276`). These are the
  SP-B anchor targets (G2).
- Ensure NeuGc / Kdn / HexA / sulfo / phospho exist in the DB or non-human /
  sulfated glycopeptides silently drop (coverage hole, not FDR miss). **NeuGc
  307.0903 ≠ NeuAc 291.0954** — mouse data (PXD005411) uses NeuGc; mixing them
  shifts sialylated masses by 16.00 Da.

**Notation** (`30-standards/notations.md`): canonical internal form = **6-tuple
(Hex, HexNAc, Fuc, NeuAc, NeuGc, Other)** with explicit residue masses; **emit**
parenthesized-count composition strings, e.g. `Hex(9)HexNAc(2)Fuc(1)NeuAc(2)` —
which andes already uses and every cited engine consumes. Composition-only (no
topology) is a complete-enough canonical form (a commercial glyco engine states this explicitly).
- Single letters are dialect-overloaded (a glyco search engine H/N/A/F/G ≠ Oxford A/G/S) — parse
  by explicit residue, not letter. Safe aliases: dHex≡Fuc, Neu5Ac≡NeuAc,
  Neu5Gc≡NeuGc.
- **Cross-engine join key = `(raw_file_stem, native scan number)`**; parse
  `scan=N` from the spectrum id (mzML index ≠ scan); charge is a tiebreaker only.
  Fold I→L for the join key only; keep original in the report.
- WURCS/GlycoCT/IUPAC = **read-only, best-effort, flagged low-confidence**; never
  emit as canonical, never reimplement a topology parser.

## 5. Data / harvest + train/eval split policy

**Frozen EVAL holdout (never train on any file from it):** **PXD025455** — human
serum, QE-HF, stepped-HCD. Hold out **all files** (same serum pools / instrument /
prep → leakage risk). Eval truth = **multi-engine** (the reference engine ∩ a glyco search engine/a commercial glyco engine)
to break the the reference engine-both-sides loop, since our 523-scan truth is an own
re-search of a a commercial glyco engine dataset (`00-context/01`, Agent 3).

**TRAIN corpus (harvest on Codon, mixed species for anti-leakage;
`40-data/pride-datasets.md`, all PRIDE-verified):**

| Rank | Dataset | Matrix / species | Instrument / activation | Notes |
|---|---|---|---|---|
| 1 (anchor) | **PXD005411** | mouse brain | LTQ-Orbitrap stepped-energy | a glyco search engine2 reference set; clean per-GPSM FDR tables; **mouse = anti-leakage** vs human eval; sparse-b/y regime. Real-harvested sample already at `40-data/collection/psms_pxd005411.tsv` (45 balanced PSMs, provenance in `README_pxd005411.md`). |
| 2 | PXD016175 | human IgG plasma | Lumos HCD | a glyco search engine2 Results.zip. |
| 3 | PXD030670 | human saliva | QE HCD | a commercial glyco engine xlsx truth; **closest instrument match** to eval; sample spectra at `40-data/collection/spectra_examples.md`. |
| 4 | PXD020254 | human cell/tumor | Lumos stepped-HCD ±10% | rar archives. |

**Demote / drop:** PXD011239 (EThcD, not HCD — different fragmentation) → demote to
train-only if used at all; PXD057219 (venom) has **no deposited glyco results** →
**drop** from labelled train.

**Label policy — tiered, NOT pure consensus** (`00-context/01`, Agent 3):
- **Tier A** = multi-engine consensus (peptide + glycan-composition agreement) →
  gold / calibration.
- **Tier B** = single-engine, **Y-ion/oxonium-confirmed hard cases**. Consensus
  alone selects easy well-fragmented spectra and *excludes the sparse-b/y cases the
  model exists to rank* — decouple label-trust (glyco evidence) from difficulty.
- Join on `(raw_file, scan)`; canonical glycan-composition tuple (§4) for
  cross-engine matching; target **≥100–300k glyco-PSMs**.

**Glycan DB for search** (`40-data/glycan-db.md`): generate from **published
biosynthetic rules** (N 2–8, H 3–12, F 0–3, S 0–5, G 0–2; core floor H≥3/N≥2;
oligomannose N2H5–N2H9 no F/S/G; fuc ≤ hexnac; sialic ≤ hexnac−2). Do **not** import
any engine's list. Validate: ~600 common tier ⊆ (GlyConnect human ∪ a glyco search engine
built-in); all 523 PXD025455 truth glycans present in the 2510 full list (miss =
hard generation gap); |Δmass| < 1 mDa vs GlyConnect canonical strings.

## 6. Hard constraints (non-negotiable)

1. **FDR = Percolator ONLY. Never Mokapot.** 2D-FDR is a **thin Percolator
   post-process** (two vanilla runs + inclusion-exclusion merge, §G3′). No
   andes-internal FDR engine, no finite-mixture-model FDR. Percolator auto-detects
   Concatenated vs Separate from the PIN — grep the mode; cross-mode counts are not
   comparable.
2. **Clean-room, published papers only.** Reference-safe: **a glyco search engine/a glyco search engine**
   (algorithms published; *runtime* is license-gated at i.pfind.net — paper-only,
   do **not** vendor binaries/GUI), **a cross-spectrum glyco engine** (Apache-2.0,
   github.com/DICP-1809/a cross-spectrum glyco engine — the G4 reference), **GlycReSoft**
   (Apache-2.0). **O-Pair/an open-source glyco engine is GPL-3.0** (copyleft — cite the *paper*
   (PMC7606753, PMC8933705); do **NOT** copy code into Apache andes; supersedes the
   earlier "MIT" note). **SugarPy GPL-3.0 — cite paper, do NOT read code.**
   **FORBIDDEN code:** **a commercial glyco engine** (commercial, Protein Metrics/Dotmatics) and
   **the reference glyco engine/FragPipe** (UM-proprietary, academic-only) — labels/notation
   docs only, never algorithms.
3. **Additive-only PIN features.** The Y0/Y1 anchor and all glyco evidence are
   *additive* PIN columns; never fuse glycan evidence into the peptide ranking score
   and never modify existing PIN features (parity-safe rule; modifying-existing
   regresses Percolator).
4. **Differentiate, don't clone.** andes = glycan-Y-first generation + own learned
   regime-matched model + in-process RT-gated cross-spectrum + Percolator-native
   2D-FDR. No competitor combines all four. Not a re-implementation of
   the reference engine/a comparison search engine fragment-index open search.
5. **Kill-gate on separation, not find-rate.** Generation is solved; the honest
   metric for every ranking/scoring change is decoy-separated top-1 ranking.
6. **Entrapment before ship.** Validate 2D-FDR honesty by padding the DB with a
   foreign glycome+proteome (a glyco search engine2's validator, PMC5585273); any GPSM with a
   foreign-only glycan OR foreign-only peptide is a false positive → FDP check.

## 7. What NOT to do (explicit anti-goals)

- **Do NOT feed glycan-decoy rows into one unified Percolator `Label` pile.**
  Measured to crash recovery 29.4% → 4.4% (they differ only in YLadder). Two
  separate axes, always.
- **Do NOT regress to peptide-first / mass-offset generation as the spine**
  (the reference glyco engine, a commercial glyco engine, a cross-spectrum glyco engine, StrucGP, O-Pair all do this). It
  reintroduces the ~20-candidates/scan combinatorial false-match pile behind
  andes's 1.2:1 target:decoy, and gives up the precursor-anchored glycan-Y-first
  edge (+7–10 pts generation). Keep peptide-first as a **fallback branch only**.
- **Do NOT copy static hand-tuned weights** — a glyco search engine's `w≈0.35 / α / β / γ`,
  a commercial glyco engine's `w_peptide ≫ w_Y ≫ w_oxonium`, a cross-spectrum glyco engine's flat `freq^0.3` empirical
  table, O-Pair/an open-search PTM tool Table S1 α constants. **Borrow the score
  *decomposition*** (intensity × quartic-mass-error × ion-ratio × core-ratio) but
  **LEARN the combiner** regime-matched (stepped-HCD, own models).
- **Do NOT treat per-spectrum scoring as the fix for the sparse stepped-HCD b/y
  stratum.** It is physically bounded at ~11% direct-b/y; the escape is
  cross-spectrum transfer (G4).
- **Do NOT build an O-Pair-style localization graph.** Single-sequon N-X-S/T
  N-glyco has the site fixed by the sequon; a graph buys nothing (it solves O-glyco
  multi-site ambiguity).
- **Do NOT build brute-force open/wildcard mass search** or re-implement
  the reference engine/a comparison search engine fragment-index open search — the reference engine/a commercial glyco engine own that niche and it
  inflates the crowded candidate space.
- **Do NOT prefer GBDT for the SP-B fragment-reliability model.** GlycReSoft tested
  gradient-boosted trees, saw substantial overfitting, and chose a regularized
  multinomial/log-linear model (PMC11263600). SP-B should be regularized /
  regime-matched, not GBDT.
- **Do NOT emit topology (WURCS/GlycoCT/IUPAC) as canonical** or reimplement a
  topology parser. Composition-only 6-tuple is the canonical form.
- **Do NOT run the multi-engine benchmark's speed arm as a Codon array** —
  heterogeneous nodes make speed counts non-portable (counts portable, speed not).

## Code anchors (from `00-context/02-code-inventory.md`)

- **G0 DET-1:** `crates/search/src/glyco_search.rs` region + `hybrid.rs:318` → `total_cmp`.
- **G2 anchor / scoring seam:** ranking call `glyco_search.rs:577-578` uses
  `prepared.scorer` (the intact standard model, `andes.rs:1532/:1613`) — glycan
  evidence stays additive PIN (`GlycoPsmKey`), never fused.
- **G3 protocol wiring:** `Protocol::NGlyco` → `crates/model/src/protocol.rs:4`
  (+name/from_name `:14,:26`) + CLI enum `crates/andes/src/bin/andes.rs:69`;
  `build_selection_key` arm `andes.rs:4809`; `protocol_to_experiment_class`
  `crates/model-train/src/store/read.rs:255` (iTRAQPhospho precedent `:262`);
  `select_nearest` `select.rs:261` already WARN-degrades.
- **G3 training rows:** glyco corpus builder feeding **glycan-stripped backbone**
  rows (new; `crates/model-train/src/accumulate.rs:62` is bare-peptide today;
  Y-ions optionally `loss_class=1`, reserved in `store/schema.rs:211`, no producer
  yet). Truth-TSV→parquet adapter is new.
- **G4 scaffold (COMPLETE, gated OFF):** `ANDES_GLYCO_CROSSSPECTRUM`
  (`glyco_search.rs:203`), two-pass `:790-887`,
  `crossspectrum.rs GlycoformWhitelist::transfer :58` — needs RT gating + truth A/B.
- **G3′ 2D-FDR:** glycan-decoy rows already emitted `glyco_pin.rs:320-330`.

## Cited sources

- Fang et al., a cross-spectrum glyco engine, *Nat Commun* 2022 — PMC8990002 (cross-spectrum
  transfer; ~11% direct-b/y ceiling; +33.5–178.5%). Apache-2.0, github.com/DICP-1809/a cross-spectrum glyco engine.
- Riley & Malaker, *JPR* 2020 — PMC7425838 (energy partitioning, glycan-channel
  dominance).
- Liu et al., a glyco search engine2, *Nat Commun* 2017 — PMC5585273 (2D inclusion-exclusion FDR;
  entrapment validator; glycan-decoy recipe). PXD005411.
- Zeng et al., a glyco search engine, *Nat Methods* 2021 — PMC8648562 (glycan-Y-first ion index;
  biosynthetic DB). Repo Apache-2.0, runtime license-gated (i.pfind.net).
- O-Pair, *Nat Methods* 2020 — PMC7606753; Multi-attribute glycan score, *MCP* 2022
  — PMC8933705. an open-source glyco engine **GPL-3.0** (paper-cite only).
- GlycReSoft — PMC11263600 (Apache-2.0; separated axes; GBDT-overfit finding).
- ppmFixer, *Glycobiology* 2024 (6-decimal mass necessity).
- Repo internal: `PHASE1_RESULT.md`, `SPA2_RESULT.md`, `00-context/*`,
  `30-standards/*`, `40-data/*`.
