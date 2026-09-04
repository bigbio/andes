# andes `--glyco` Mode — Design Improvements (v2 Addendum)

**Status:** design addendum. Amends
`docs/specs/2026-06-13-andes-glyco-mode-design.md` (the v1 spec). Synthesized
from 6 per-dimension SOTA research reviews (oxonium trigger, backbone-mass
solver, composition DB / decoy / 2D FDR, scoring integration, output/interop,
Phase-0 validation), reconciled against the codebase and the MVP scope.

**MVP scope is unchanged and binding:** N-linked glyco only, sequon-restricted
(N-X-S/T, X≠P), stepped-HCD / HCD, curated **composition** list (no de-novo, no
structure/linkage), Phase-0 prototype-gated. Anything below that drifts toward
O-glyco, ETD/EThcD, glycan trees, open mass windows, or a hand-tuned mixing
weight is explicitly **deferred** and called out as such. Each subsection says
**what it replaces** in v1.

**Conventions reaffirmed (codebase invariants, do not violate):**
- PIN gets **additive numeric feature columns only** — `Proteins` is
  rest-of-line (`crates/output/src/pin.rs` writes one accession per
  `candidate_idxs` then `writeln!`), so any trailing non-numeric field becomes a
  phantom protein and corrupts Percolator. (n=12: additive features safe;
  modifying existing regresses — see `parity-tuning-lessons`.)
- Human-readable strings (composition, accession, q-value) go to **TSV only**,
  extending the existing additive-last `Modifications` `pos:CURIE` column
  (`crates/output/src/tsv.rs::modifications_field`).
- Optional scoring-path behavior must be **inert / byte-identical when unused**
  (the neutral-loss-disabled invariant; CLAUDE.md).
- Let **Percolator learn** the peptide-vs-glycan weight via additive features —
  do **not** hand-tune pGlyco's `w` / GlycReSoft's `w`.

---

## 0. Reconciliation summary (conflicts resolved)

| Conflict across reviews | Resolution |
|---|---|
| Trigger = peak **count** (v1) vs **summed-intensity** gate (MSFragger/pGlyco) vs **k-of-n peak-presence** (Byonic) | **Two-part AND gate**: summed-oxonium-fraction ≥ floor (primary, calibrated) **AND** ≥2 distinct core ions above a rank/intensity floor (anti-false-trigger). Replaces v1's count-only rule. |
| Backbone solver = **direct Y0/Y1 read** (v1) vs pGlyco3 **Y-complementary voting** vs MSFragger **glycan-stripped b/y** | **Converged solver**: Y-complementary voting is primary; ≥2-core-Y quorum gates it; on marginal quorum emit **top-K** candidates; glycan-stripped backbone **b/y RawScore is the final arbiter**. Direct Y0/Y1 read demoted to "Y0 counts as a core ion," not the sole derivation. |
| 2D FDR = **inclusion-exclusion** `F_GP = F_G + F_P − F_(G∩P)` (v1, pGlyco-style) vs **sequential** TDC (FragPipe/PTM-Shepherd) | **Sequential** TDC: peptide q (Percolator, untouched) then glycan q (new post-Percolator TDC layer), gate on `max(q_pep, q_gly) ≤ α`. Replaces inclusion-exclusion as the primary; the latter kept only as a documented alternative. |
| Hard discard vs soft route at the trigger | **Soft route**: the gate decides only which scans pay for the backbone solve (cost control); the continuous oxonium evidence flows to Percolator as features so a borderline true glyco scan that clears Y-quorum is not lost. |
| "Reuse the neutral-loss primitive scores the whole Y-ladder" (v1) | **Half true.** The primitive scores **only glycosite-spanning backbone b/y ions that retain a partial glycan** (a span containing the loss-bearing Asn). The **Y0/Y1/Y2 whole-backbone ladder** is a precursor-region series `score_psm`'s cleavage loop (sites `1..n`) structurally cannot emit → it needs a **separate cheap peak-lookup matcher** feeding features. |

**Dropped as over-scoped for the MVP** (recorded so they are not silently lost):
NeuGc/phospho/sulfate composition-class branching beyond a single feature flag;
runtime OLS/GlyTouCan accession queries; pGlyco's 3-level finite-mixture FDR;
single-energy-HCD threshold auto-calibration (measure, don't auto-tune in v1);
the multiplicative joint-FDR estimate. Charge-explicit voting is **kept** (it is
a correctness fix, not scope creep) but bounded to the precursor charge range.

---

## 1. Oxonium trigger — replace count-only with a calibrated two-part gate

**Replaces:** v1 Phase-0 "Detector 1" and Phase-1 trigger
(*"flag an MS2 as glyco when ≥2 of the canonical oxonium ions are present above
an intensity floor"*).

### 1.1 The gate (two parts, ANDed)
A scan is **glyco-flagged** iff:
- **(A) Summed-intensity gate (primary):** Σ(intensity of matched core oxonium
  ions) ≥ `--oxonium-min-frac` × base-peak intensity. **Default 0.10** (the
  MSFragger/pGlyco3 `diagnostic_intensity_filter` anchor). Tunable because the
  operating point is activation/instrument dependent (AI-ETD uses 0.025).
- **(B) k-of-n presence guard (anti-false-trigger):** ≥2 distinct **core**
  oxonium ions present, each above an absolute floor — **intensity ≥ 1% of base
  peak** AND (Byonic-style) ranking within the top ~50 peaks. Prevents a single
  coincidental 204.087 noise peak from clearing the gate in a sparse spectrum.

Base-peak normalization **reuses the existing pattern** in
`crates/scoring/src/scoring/strong_score.rs` (`peaks.map(|(_,i)| i).fold(0.0,
f64::max)`) so the trigger and the intensity model share one notion of base peak.

### 1.2 Activation gate (free, reuse existing helper)
Fire the gate **only** when `activation_method` is collisional. Reuse
`ActivationMethod::predicts_neutral_losses()` (already `false` for ETD) — never
fire on pure ETD (no oxonium under electron dissociation). For single-energy
high-NCE HCD note that oxonium can dominate and suppress backbone ions; the
`--oxonium-min-frac` knob covers this per-acquisition rather than a hardcoded
energy branch. (Open question in v1 — partly answered: gate logic is identical
across HCD energies, **only the threshold differs**.)

### 1.3 Pinned oxonium ion set (single source of truth)
Define **one** constant adjacent to `loss_class_id` in
`crates/model/src/modification.rs` (or a new `glyco::oxonium` module) so the
trigger ions and the scored glyco loss ions never drift. Split it:

- **CORE (drives gate parts A+B), HexNAc series + Hex-HexNAc:**
  `138.0550, 168.0655, 186.0761, 204.0867, 366.1395`.
- **EXTENDED (glycan-type feature flags only — NOT required for the gate):**
  `126.0550, 144.0655` (HexNAc fragments), `163.0601` (Hex — **missing in v1,
  add**), `274.0921, 292.1027` (NeuAc → sialylation), `512.1974`
  (HexNAc-Hex-Fuc — fucose marker, **missing in v1, add**), `657.2349`
  (HexNAc-Hex-NeuAc), `290.0870 / 308.0976` (NeuGc).

**Match window:** **±20 ppm** on high-res data (andes targets Q Exactive /
Astral) with a **0.01 Th absolute floor** for the low-m/z ions. **Specify ppm,
not a flat Da** — a flat 0.01 Th mis-matches the heavy markers (366, 512, 657).
This resolves v1's "no tolerance/units stated" gap.

### 1.4 Soft route, not hard discard
The gate's binary decision controls **only** which scans pay for the Y0/Y1
backbone solve (a cost knob). The **continuous** oxonium evidence is emitted as
additive PIN features (§5) so Percolator makes the final identification call. A
hard discard would create an invisible recall ceiling (small high-mannose /
under-fragmented glycans at 8% oxonium); the soft route avoids it.

---

## 2. Backbone-mass solver — converged voting solver with explicit fallbacks

**Replaces:** v1 Phase-0 "Detector 2" and Phase-1 step 1 (a Byonic-style
**direct Y0/Y1 read** — the single most fragile published strategy: Byonic
itself notes absent Y0/Y1 "reduces confidence").

### 2.1 Primary: Y-complementary-mass voting (pGlyco3-style), DB-anchored
Precompute, from the curated ~182-composition list, the table of
`{glycan_total − glycan_Y_rung}` complementary masses (a fixed, offline table).
For each MS2 peak `p` with charge-deconvolved **neutral** mass `m_p`, compute
`c = precursor_neutral − m_p` and look `c` up against the table. Each hit votes
for that glycan composition. Then `backbone = precursor_neutral −
glycan(top-voted)`.

This makes **every rung** of the Y-ladder (Y2, core-trimmed, Hex-trimmed)
contribute evidence — not just Y0 and Y1 — so a single weak/missing Y0/Y1 does
not sink the solve. It is a **pure pre-scoring lookup** over the peak index
(`ScoredSpectrum::nearest_peak_full` / `nearest_peak_rank`,
`parent_mass()`), **zero scoring-core change**, and snaps directly onto v1's
"snap to the composition list" step.

### 2.2 Charge-explicit (correctness fix, not scope creep)
Deconvolve candidate Y-region peaks to **neutral** mass across `z =
1..precursor_z` and vote in **neutral complementary-mass space**. If the
precursor charge is ambiguous/multiple, run voting per candidate precursor
charge and keep the `(charge, glycan, backbone)` triple with the strongest
core-Y quorum — the quorum itself disambiguates charge (wrong charge → votes
land on noise). andes already iterates candidate charges in the search loop and
has the isotope-error window; the solver joins that iteration. Cost is
peak-lookup only.

### 2.3 Core-Y quorum gate
Accept a single backbone solution only when **≥2 core Y-ions** support the
winning composition (pGlyco3's N-glyco threshold; **Y0 always counts when
present**). Define core Y-ions for N-glyco as the trimannosyl-core ladder: Y0
(bare), Y1 (+HexNAc 203.0794), Y2 (+2×HexNAc 406.1587), and inner-core Hex steps
(+162.0528 off Y1). Below quorum → do **not** emit one backbone; emit the
**top-K** voted compositions (K small, **3–5**) as ambiguous candidates.

### 2.4 Bounded Route-(ii) fallback (concrete, not a blind 182-window search)
When quorum fails, restrict per-spectrum windows to the **union** of:
(a) the **top-N compositions by complementary-vote count** (N ≤ **8**), and
(b) **all small glycans (≤3 monosaccharides)** — the HexNAc / HexNAc₂ /
HexNAc₂Hex₁ core series pGlyco3 keeps unconditionally because they can never
accrue enough Y-votes. Open one window per candidate at
`precursor_neutral − mass(glycan)` via `candidate_nominal_bounds`
**unchanged**. **Cap total fallback windows per spectrum** (default ~**10**).

Rationale: the research plan establishes the GF mass-indexed DP node cache
**cannot amortize across compositions** (each backbone mass = a different node
set), so a naïve 182-window fallback is 182× the per-scan cost. The shortlist is
driven by the **same** complementary-vote table — it is the marginal branch of
one solver, not a separate brute-force mode.

### 2.5 Final arbiter: glycan-stripped backbone-b/y cross-check (MSFragger-style)
For each top-K backbone candidate, the deglycosylated backbone **b/y RawScore**
is computed by the **existing `score_psm` unchanged** and selects the winner.
Even when the Y-ladder is too sparse to clear quorum, the peptide's own backbone
fragments (stepped-HCD produces them in the same scan) confirm the right
candidate. This adds **no cost** beyond the bounded top-K windows already opened
(§2.4) and keeps the trained scoring — not a heuristic peak read — as the final
judge.

**Net solver order:** complementary voting → quorum check → (single backbone) OR
(top-K bounded windows + small-glycan safety net) → backbone-b/y RawScore picks
the winner.

---

## 3. Composition DB format, default list, mass-based decoy, 2D FDR

**Replaces:** v1 "Glycan database & decoy strategy" (lines 100–123) — sharpens
the under-specified format, decoy, and FDR.

### 3.1 DB format (drop-in FragPipe/Byonic interop)
- One glycan per non-blank, non-`#`-comment line; **no header**.
- Primary syntax = **Byonic composition string** `Residue(Count)Residue(Count)…`
  with an **optional ` % <mass>` suffix** (literal `%` is the mass separator;
  mass in Da; **computed by andes when absent** — the common case).
- Recognized residue registry (fixed): `Hex, HexNAc, Fuc(=dHex), NeuAc, NeuGc,
  Pentose(=Xyl), Phospho, Sulfo`. Case-insensitive; unknown residue → **hard
  parse error naming the token**.
- Stretch-accept (same loader): MetaMorpheus `Residue-Count_Residue-Count` and
  pGlyco nested `(N(N(H(A))))` as **input aliases only**.
- Extensions: `.txt / .csv / .tsv / .glyc`.
- **Always recompute mass from composition**; if a supplied `% mass` disagrees
  by >1 mDa, **warn (not error)** and use the computed value.

This is a near-verbatim adoption of the FragPipe loader → every published
FragPipe/Byonic/pGlyco list loads unchanged, at zero cost.

### 3.2 Monosaccharide residue masses (water-subtracted, in-chain monoisotopic)
Ship as a `const` table in the new `andes-glyco` crate, mirroring
`crates/model/src/mass.rs`:

| Residue | Mass (Da) | | Residue | Mass (Da) |
|---|---|---|---|---|
| Hex | 162.05282 | | NeuAc | 291.09542 |
| HexNAc | 203.07937 | | NeuGc | 307.09033 |
| Fuc / dHex | 146.05791 | | Pentose / Xyl | 132.04226 |
| Phospho | 79.96633 | | Sulfo | 79.95682 |

`glycan_delta = Σ(count_i × residue_mass_i)` with **NO additional water** (the
reducing-end –OH is part of the Asn-bearing peptide backbone, already in the
andes peptide mass). This note prevents the classic +18 Da glycan-mass bug. The
**Y1 anchor = backbone + HexNAc(203.0794)** used by §2 comes from this same
table.

### 3.3 Default list — provenance (independence-critical)
**Do NOT claim to ship "the published 182" as a single licensed file — none
exists Apache-clean.** Instead:
- Curate `resources/glycans/human_nglycan_182.txt` **clean-room** by intersecting
  the **MIT** MetaMorpheus and **Apache-2.0** Glyco-Decipher N-glycan lists
  (both license-safe to read/adapt) down to the ~180–190 most-frequent human
  serum/tissue N-glycan compositions, cross-checked against the public
  GlyGen/GlyConnect composition browser.
- **Record the derivation + per-composition source note** in a header comment
  (same provenance discipline the independence memo applies to the models — a
  false-provenance claim is the failure mode to avoid).
- Tag it with the existing `glyco` catalog slug
  (`crates/model-train/src/catalog.rs`).
- pGlyco's 1,234 and FragPipe's full 1,670 lists are loadable via `--glycan-db`
  for deep profiling (the 182 is a curated **subset**, not a competing
  authoritative file).

### 3.4 Mass-based glycan decoy (FragPipe 4-type, default Type 1)
Expose `--glycan-decoy-type {0,1,2,3}`, **default 1**:
- **Type 1 (default):** decoy intact mass = random offset within the glycan match
  tolerance **+ a random isotope error** from the search's isotope-error set;
  decoy keeps the target's nominal composition and per-type fragment **count**,
  but **each Y and oxonium fragment is shifted by a unique random value in
  [1,20] Da**.
- **Type 0:** ±3 Da window (looser, discouraged).
- **Type 2:** as Type 1, no isotope error.
- **Type 3:** decoy mass == target mass, mass-error excluded from scoring
  (stringent).
- One decoy per target glycan, **1:1 paired**, labeled `glycan_label = -1`.

The fragment scramble (not composition) is what makes glycan decoys **orthogonal
to andes's reverse-peptide decoys** (`crates/search/src/decoy.rs`) — the two
compose. The per-fragment shift is applied to the **loss-shifted ions
`predict_by_ions_with_losses` emits** so the decoy genuinely degrades the
Y-ladder/glyco score, not just the precursor match.

### 3.5 2D FDR = sequential TDC (NOT inclusion-exclusion)
**Replaces v1's `FDR_GP = FDR_G + FDR_P − FDR_(G∩P)` as the primary.**
- (a) Peptide/PSM FDR **unchanged** — reverse-sequence decoys → Percolator on the
  PIN.
- (b) Glycan FDR = a **new post-Percolator layer** in `andes-glyco`: for every
  PSM passing peptide `q ≤ α`, compute the multi-attribute glycan **absolute
  score** (Y-ladder count + intensity ratio, oxonium count, glycan mass error,
  isotope error), run **target/decoy-glycan competition** over the whole run's
  glyco-PSM population (`q_glycan = #decoy-glycans / #target-glycans` at each
  score cutoff), keep `q_glycan ≤ α`.
- Report **both q-values per row**; gate on **`max(q_pep, q_gly) ≤ α`**
  (conservative upper bound on the joint FDR).
- Inclusion-exclusion and the multiplicative `1−(1−q_p)(1−q_g)` estimate are
  documented as **alternatives only** (no engine validated them; PTM-Shepherd
  filters sequentially).

**Glycan-FDR granularity caveat:** with ~182 targets ⇒ ~182 decoys, glycan
q-values are coarse at small N. Estimate FDR over the **whole run's glyco-PSM
population**, never per-spectrum, or q-values are unstable.

---

## 4. Scoring integration — split the two glyco fragment families

**Replaces:** v1 Phase-1 step 4 (*"the Y-ion ladder rungs … ARE the loss-shifted
b/y ions"* — overstated). The mapping splits cleanly:

### 4.1 Backbone-retaining ions → the neutral-loss primitive (as-is)
The primitive (`predict_by_ions_with_losses` / `loss_node_score`, and
`span_losses` in `scored_spectrum.rs`) emits a loss-shifted partner only for b/y
ions whose residue **span contains the loss-bearing Asn** — i.e. exactly the
**glycosite-crossing backbone b/y ions that retain a partial glycan**. These are
scored by the shipped pooled per-class glyco loss table (`loss_class=1`).

**Declare the per-candidate retained-core loss ladder.** The glyco mode
**synthesizes a `Modification` on the sequon Asn at assignment time** (NOT from a
static `mods.txt` — there is no single loss list covering 182 compositions) with
`neutral_losses = { glycan_total − retained_core }` for `retained_core ∈ {0,
203.0794 (HexNAc), 406.1587 (2×HexNAc), 568.2115, 730.2644, 892.3172 (inner-core
+Hex steps)}`, **intersected with what the composition can actually shed**,
`loss_class = glyco(1)`. The peptide then passes through the **unchanged
`score_psm`**.

**Charge-division correctness note (load-bearing).** In
`visit_directional_loss_ion_matches` the loss is divided by the **ion charge**
(`theo_mz = base_mz − loss/ion_charge`) and the ion is dropped when
`theo_mz ≤ 0`. An **intact** N-glycan (>1500 Da) declared as a single loss on a
charge-1 backbone b/y ion goes negative and is **silently skipped** → zero glyco
evidence on exactly the singly-charged ions stepped-HCD produces. Therefore
declare the **per-rung trim masses** (a Hex 162.0528, a HexNAc 203.0794, core
steps) — modest rungs that keep `theo_mz > 0` — matching the shipped glyco
example (`loss=162.0528;324.1056`), **never the full glycan adduct mass**.

### 4.2 Y0/Y1/Y2 core ladder + oxonium → a separate cheap matcher (NET-NEW)
The whole-backbone Y0/Y1/Y2 ladder is a **precursor-region** series that
`score_psm`'s cleavage loop (sites `1..n` only) **structurally cannot emit**.
Score it with a small Rust port of the Phase-0 reader (pure peak lookups via the
**same `ScoredSpectrum` index** — `nearest_peak_full` / `nearest_peak_rank` at
charge-reduced precursor m/z), producing the features in §5. No new `IonType`,
no scoring-core change.

### 4.3 Combination = additive PIN features, Percolator learns the weight
The backbone **RawScore** (now including the §4.1 glyco loss-ladder contribution
already summed into RawScore) stays the primary ranking score. The §4.2 Y-core /
oxonium evidence is **orthogonal** and is emitted as additive PIN feature columns
(§5). **Load-bearing rationale:** because the loss-ladder is inside RawScore but
the Y-core/oxonium evidence is not, RawScore alone under-ranks glyco PSMs against
non-glyco rows; the additive columns restore separation **without perturbing the
standard-search RawScore distribution**. Do **not** fold glycan score into
RawScore with a hand weight.

**v1 pooling caveat retained:** all composition-specific rungs share one pooled
`loss_class=1` distribution; glyco-PSM gain may be capped until SP3 trains on
real glyco data (deferred pending the user's glyco dataset). Acceptable for v1.

---

## 5. Output + interop columns

**Replaces:** v1 Phase-1 step 6 (*"additive glyco PIN/TSV columns (glycan
composition, glycan score, Y-ion count)"* — conflates the two channels and would
put a non-numeric string into PIN, corrupting Percolator).

### 5.1 PIN — additive **numeric** features only, gated behind `--glyco`
Insert in the existing feature block (between `StrongScoreCal` and
`Peptide`/`Proteins`), populated from a new glyco sub-block on `PsmFeatures`
(`crates/search/src/psm.rs`) via `compute_psm_features`:

| Column | Source |
|---|---|
| `OxoniumIonCount` | # core oxonium ions matched |
| `OxoniumIntensityFraction` | summed oxonium intensity / base peak (the §1 gate signal) |
| `YLadderIonCount` | # matched Y0/Y1/Y2 + Hex-ladder rungs (§4.2) |
| `YLadderIntensityFraction` | summed Y-rung intensity / MS2 TIC |
| `GlycanMassErrorPpm` | ppm(precursor − backbone − assigned-composition mass) |
| `GlycanScore` | the pooled glyco loss-table node-score sum (§4.1) — **not** a bespoke counter |
| `GlycanDecoy` | 0/1 glycan-level TD label (feeds §3.5 2D FDR) |

**Gate the columns behind `--glyco`** so the standard-3 search stays
byte-identical and the existing `pin_header_columns_are_gf_free_schema` golden
test is untouched (a re-baseline would otherwise be forced). **Start with the 3
highest-information features** (`OxoniumIntensityFraction`, `GlycanScore`,
`GlycanMassErrorPpm`) and let the benchmark justify the rest — 7 mostly
correlated features on a sparse glyco PIN risk overfitting.

### 5.2 TSV — human/interop string columns (rest-of-line-free)
Append as additive-last columns, and also add a `pos:CURIE` entry to the existing
`Modifications` column when an accession resolves:

`GlycanComposition` (Byonic `HexNAc(2)Hex(5)NeuAc(0)Fuc(0)`), `GlycanMass`
(4–5 dp), `GlycanScore`, `GlycanQValue` (glycan 2D-FDR q), `GlycanSite` (1-based
sequon-N), `GlycanAccession` (GNO/GlyTouCan CURIE or empty). **All empty for
non-glyco rows → byte-identical to today's TSV.**

### 5.3 Composition string format (canonical)
Byonic `Monosaccharide(count)` using the **ProForma-2.0 token set**
`{Hex, HexNAc, dHex, NeuAc, NeuGc, Pen, HexS, HexP, HexNAcS}` (dHex == Fuc;
zero-count residues may be omitted). Store this **exact string per entry in the
DB file** so output is a lookup, not a re-derivation. pGlyco `H/N/A/F/G`
shorthand is an **input alias only**, never emitted.

### 5.4 Accession scheme — optional, tiered, **GNO not UNIMOD**
Primary identity = composition string + monoisotopic mass (always emitted). The
**`GNO:Gxxxxxxxx`** CURIE is an optional enrichment emitted only when the curated
DB entry carries a pre-resolved accession — **ship the 182-list with GlyTouCan
composition-level accessions pre-populated where known, empty otherwise. Do NOT
runtime-query OLS/GlyTouCan** (zero-config philosophy). GNO is the ProForma-2.0
glycan ontology; **never force a `UNIMOD:` accession** (UNIMOD has no
per-composition glycan terms — it would be semantically wrong).

### 5.5 quantms / PSI interop
Emit a **ProForma-2.0 glyco peptidoform** as the cross-tool key:
`…N[Glycan:HexNAc2Hex2]X(S/T)…` when no accession, `…N[GNO:Gxxxxxxxx]…` when
known. Map each glyco PSM to quantms's structured `modifications` element:
`name` = Byonic composition, `accession` = GNO CURIE-or-null,
`positions=[{position: sequon-N 1-based, amino_acid:'N', scores:[{glycan_score},
{glycan_q_value}]}]`. quantms's `from_proforma` already round-trips bracket tags,
so the output is directly ingestible. Pin the token map (`Fuc == dHex`) and
unit-test the round-trip against pyteomics/spectrum_utils.

---

## 6. Phase-0 go/no-go — tightened, executable

**Replaces:** v1 Phase-0 "Validation metrics" + "Decision".

### 6.1 PRECONDITION (v1 is factually wrong here)
The Byonic `.pepXML` is **NOT staged** — `/tmp/glyco_andes/` currently holds only
`Pool_HCC_early_Fc3_1.mzML`, `.raw`, `human_sp.fasta`, `std.pin`, and
`glyco_mods*.txt`. **Phase 0 cannot run as written.** Make this an explicit,
verified precondition:
- Download `Pool_HCC_early_Fc3_1.pepXML` (or the project's Byonic `.xlsx`/`.mzid`
  results) from PRIDE `…/2021/05/PXD025455/` into `/tmp/glyco_andes/`; **assert
  >0 glyco rows parse** before running.
- **Fallback if no per-scan Byonic file exists:** run pGlyco3 or MSFragger-Glyco
  on the staged `.raw` at 1% peptide × 1% glycan FDR to produce an open,
  per-scan, license-clean truth (also needed for the §6.4 benchmark).

### 6.2 Two truth sources, reported separately
Byonic PXD025455 was filtered at **score ≥ 150** → a high-confidence but
**incomplete, biased** subset (the 2025 comparative eval found all-5-tool overlap
of only 17%; 33.9% tool-unique; Byonic itself made 0.0065 Da peptide-pair
misassignments). Report **both**:
- **(a) Recall vs Byonic-confirmed scans** — a lower bound ("andes misses what
  Byonic found").
- **(b) Orthogonal diagnostic-ion truth** — a scan is "glyco-confirmed" if it has
  **≥4 of 5** HexNAc oxonium ions `{126.055, 138.055, 168.065, 186.076,
  204.087}` **AND** a detectable Y0/Y1 or backbone b/y pair. This is Byonic-
  independent and is the **real** go/no-go for the trigger + Y-reader.

### 6.3 GO/NO-GO decision table
| Component | Metric | GO bar |
|---|---|---|
| Trigger | precision-recall over `--oxonium-min-frac` ∈ {0.02, 0.05, 0.10, 0.20} **and** k-of-n {204.087, 138.055, 366.140} | recall ≥ **0.85** at glyco-precision ≥ **0.90** (the **precision arm is new** — v1 set only recall) |
| Backbone | mass within ±0.02 Da of Byonic, **stratified by Y-ladder richness** (core-Y count 0,1,2,≥3) | ≥ **0.70**, **and the bar must hold within the sparse (≤1 core-Y) stratum** — the regime the direct reader fails |
| Backbone (searchable) | feed derived backbone into a **dry-run of `candidate_nominal_bounds`**; does the true backbone peptide fall in the opened window? | ≥ **0.70** — *this is the metric that actually gates Phase 1; mass accuracy alone is insufficient because the engine searches by nominal window* |
| Glycan-mass closure | `precursor − backbone` residual to nearest 182-list entry | \|residual\| < **50 ppm** (FragPipe's glycan tolerance); **also report the fraction of true glycans OUTSIDE the 182 list** (serum/HCC is heavily sialylated/fucosylated — quantifies the DB-coverage ceiling **before** committing to 182) |
| Charge confusion | fraction where a single-charge direct read picks the wrong precursor charge but neutral-space voting picks the right one | report (proves §2.2 earns its complexity) |
| Fallback | quorum-fail rate + whether the bounded top-K still contains the true backbone | report (proves §2.4 recall ceiling) |

**Decision:** GO only if the trigger and **searchable-backbone** bars clear
(stratified, not aggregate). If the sparse stratum fails → the bounded Route-(ii)
fallback (§2.4) carries the mode, stated explicitly as the fallback. **Record all
numbers back in this doc.**

### 6.4 Phase 0.5 — full-mode benchmark design (reproducible, matched)
Add before Phase 1 ships:
- **Matched FASTA + glycan list:** same target+decoy DB and same 182 composition
  list for andes, MSFragger-Glyco, pGlyco3 (pGlyco3 defaults to 1,234 — restrict
  it to 182 **or** report both and note the asymmetry).
- **Matched FDR:** **1% peptide × 1% glycan**; report glycoPSMs, unique
  glycopeptides, unique (peptide+glycan+site).
- **Primary metric = two-layered entrapment-FDP:** peptide-level via a foreign
  proteome appended to `human_sp.fasta` (Arabidopsis or yeast; target FDP ≤ 1%,
  cf. FragPipe's 0.04%) **and** glycan-level via entrapment glycans absent from
  human serum biology (NeuGc/odd compositions).
- **WIN condition (per the andes objective):** more glycoPSMs/glycopeptides than
  MSFragger-Glyco **and** pGlyco3 at matched 1%×1% FDR, entrapment-FDP within
  tolerance, plus speed. (Avoids the cross-mode / parsimony **counting-artifact**
  class documented in `andes-full-benchmark`.)

---

## 7. Open questions (updated)
- **Single-energy vs stepped-HCD:** gate logic identical, threshold differs;
  measure fallback-trigger rate on any single-energy data before claiming
  generalization. **Do not auto-calibrate the threshold in v1.**
- **Second stepped-HCD dataset:** PXD025455 is one Q Exactive HF file; resolve a
  second dataset before treating the §6.4 benchmark as conclusive.
- **182 coverage:** if §6.3 glycan-mass closure shows high failure, revisit the
  default size **before** Phase 2.
- **Tolerance coupling:** set the ≤0.02 Da backbone target jointly with the
  oxonium floor and the existing isotope-error window — too tight rejects heavy
  isotope-error-prone glyco precursors; too loose lets noise clear quorum.
- **O-glyco later:** the DB format (§3.1) must not preclude it (composition list
  is residue-agnostic) — but O-glyco scoring/quorum (≥1 core Y) is **deferred**.

---

## 8. Residual conflicts / risks (unresolved, flagged)
- **Threshold transfer:** 10% oxonium fraction is validated on Q Exactive
  stepped-HCD; **re-tune for Astral/timsTOF** via the §6.3 sweep — do not assume
  transfer.
- **Co-isolation / chimeric false triggers:** a non-glyco precursor co-isolated
  with a glyco one can inherit oxonium ions. The k≥2 + rank guard + single-
  backbone quorum reduce but don't eliminate it; the §3.5 glycan 2D FDR is the
  backstop and **must not be skipped**. (andes's `--chimeric` two-pass can later
  own the secondary species; Phase 0 just flags/down-weights multi-backbone
  scans.)
- **Near-isobaric glycans:** a confident backbone does **not** imply a confident
  glycan (NeuAc vs Fuc₂+… traps); the 2D FDR treats glycan assignment as
  separately uncertain. Composition (not structure) ambiguity is **reported, not
  silently resolved** — an explicit MVP limitation.
- **No Apache 182 file exists:** the curated list is re-derivation work (§3.3),
  not a download; must be documented per the independence memo.
- **Pooled glyco loss table** caps gain until SP3 trains on real glyco data
  (deferred).
- **Fragment-scramble decoy** can accidentally land a decoy fragment on a real
  Y/oxonium m/z, slightly deflating measured glycan FDR; entrapment-validate
  (§6.4) before claiming FDR control.
