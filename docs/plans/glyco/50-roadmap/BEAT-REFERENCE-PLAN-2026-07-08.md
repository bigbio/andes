# Beating the reference engine on Glyco — Consolidated Plan (2026-07-08)

> Produced by a 6-agent investigation loop (2 code audits, 1 literature review, 1
> hygiene/design review, 2 design agents). Every code claim below is verified against
> the tree on branch `glyco-phase1`. Test bed: PXD025455 `HCC_pool_Late_Fc3_r1`
> (523 the reference engine-truth backbones).
>
> **⚠️ RECONCILED 2026-07-08 with `SESSION-HANDOFF-2026-07-08.md` §2/§3/§5 (rewritten
> after a real the reference engine-gap audit).** The handoff is authoritative where it conflicts
> with this doc. Two corrections were folded in: (1) the glycan-decoy channel is NOT a
> cheap win — it was already run WITH expansion (the "unified glycan-decoy Percolator
> pile") and still crashed because glycan-axis features are too underpowered; decoys are
> demoted to last (P3), blocked on richer glycan features. (2) The first move is a
> **candidate-level competitor audit (P0)**, not a decoy experiment. Also added:
> `SELECTOR_SHORTLIST_K=24` prefilter and the `combined`/`YINDEX` A/B confound.
>
> **⚠️ CORRECTED 2026-07-08 (foundations bug found while running P0).** The
> `glyco_outrank_audit.py` `residue_mass()` string parser mis-parses **45.9%** of
> candidate rows, inflating `truth_absent`. Both audit tools are now fixed to match
> numerically (CalcMass − GlycanMass − H₂O). **Corrected decomposition of 523 truth:
> 148 absent (not 301) / 280 top1-by-mass (not 116) / 95 outranked (not 106); outranked
> reason 92 "truth loses summed Y-ladder" / 3 tie (not 57/49).** Net: generation is ~half
> the believed problem; SCORING dominates (280 top1-by-mass, only 97 @1%). Leg 2 is the
> clear priority; leg 1 drops. Per-leg counts below are updated to match.
>
> Backing detail: `scratchpad/glyco-loop/01..06-*.md` (scoring audit, generation audit,
> literature, hygiene+design, integrated design, FDR-safe expansion).

---

## The one-paragraph answer

Beating the reference engine requires **three interlocked fixes**, not one: (1) **precursor-mass
backbone enumeration** to reach the ~150 "unsearchable" backbones a fragment-driven
engine structurally cannot propose; (2) a **two-channel learned fusion** (peptide-b/y
`ScoreP` ⊕ glycan-Y `ScoreG`, a glyco search engine weighting ≈0.65/0.35) that enters the top-1
**selector** to un-block the 95 generated-but-outranked (corrected numeric audit; and to
close the larger 280→97 top1-survives-FDR gap); and (3) a **glycan-axis
FDR guardrail** (glycan-decoy channel) to make the expansion in (1) FDR-safe rather than
self-defeating — but the *built* version of (3) is underpowered (its decoy rows share
nearly all peptide features) and already crashed with expansion, so it is **blocked on
richer glycan features and demoted to LAST**. They are **coupled** — every past
single-lever A/B came back flat or crashed *because* the other legs were missing. The
code-audit finding: legs (1) and (3) are **already built but dormant/gated** (so this is
wire-up-and-validate, not build-from-scratch), and leg (2) is half-built. **Corrected
order (handoff §5): P0 candidate audit → leg 2 selector → leg 1 precursor-mass → leg 3
decoys last.**

---

## Why every prior attempt failed (the coupling, in one table)

| Past result | Why it failed | The missing leg |
|---|---|---|
| Force high-charge backbones into pool → `top1_correct=0` | recovered backbones land out-ranked | leg 2 (selector fusion) |
| Full 2510-glycan list → @1% crashes 253→119 | target glyco-noise inflates, peptide-decoys are mass-preserving so don't model wrong-glycan error → Percolator separation collapses | leg 3 (glycan-decoy channel) |
| Additive PIN feature for fusion → flat on top-1 | winner is chosen by `collapse_cmp` *upstream* of Percolator; PIN can't resurrect a demoted backbone | leg 2 must enter the **selector**, not just PIN |
| Loosen FDR 1%→40% → +6 correct, +200 decoys | the 407 non-correct are never *emitted* | leg 1 (generation) |

**Refuted / do-not-pursue (re-confirmed this loop):** P2 "unfair mods" is **not real** —
protein N-term acetyl has been in andes' default since commit `0829d621` (2026-06-19),
weeks before any glyco run; `--mods docs/benchmarks/configs/mods.txt` reproduces the
already-active set (≈byte-identical, ~0 ID gain). RT as a standalone lever is weak with
the seed model (parked).

---

## Leg 1 — Precursor-mass backbone enumeration (generation)

**Recovers:** the genuinely-absent backbones — **148 (28%)** by corrected numeric audit,
NOT 301/58% (that count was a parser artifact; see the correction banner). Of these,
a subset are truly unsearchable (no core-Y ladder + no b/y). Leg 1 is real but roughly
HALF the size the handoff claimed, and lower priority than leg 2.

**What's already there (verified):**
- Generation is **fragment-driven**: a backbone is nominated only if it clears a core-Y
  rung quorum (`backbone.rs:287-303`) or a ≥6-b/y peptide-first gate
  (`glyco_search.rs` `MIN_BY_MATCHES`). No source enumerates DB peptides by precursor
  mass alone.
- BUT the primitive exists: `db_branch` computes `bb = precursor_neutral − glycan.mass`
  per composition (`crates/andes-glyco/src/hybrid.rs:65-98`), and the precursor-mass
  peptide index `bucket_index: BTreeMap<i32, Vec<usize>>` is built
  (`crates/search/src/match_engine.rs:294-298`), wired into the glyco ctx
  (`glyco_search.rs:206, 427, 531`), and already used mass→peptide at
  `glyco_search.rs:849` — but only for **already-fragment-nominated** backbones.

**Change:** add ONE ungated `Source::Db` source right after the peptide-first block
(≈`glyco_search.rs:696`), per charge/isotope, sweeping the glycan list → `bb` →
`bucket_index.range(...)` (reusing the `:849` machinery + `has_nxst_sequon`), emitting
`BackboneHit`s field-identical to `:675-687`. Gate behind `ANDES_GLYCO_PRECURSOR_MASS=1`
with `max_precursor_mass` emission cap and `PRECURSOR_MASS_MAX_Z=3` (the paired
high-charge cap the audit asked for — clamps this new source to the charge regime the
core-Y match functions already trust at `backbone.rs:391`, without touching the uncapped
vote loops at `backbone.rs:186` / `glyco_y_index.rs:102`; those are a separate engine fix).

---

## Leg 2 — Two-channel learned fusion into the selector (scoring)

**Recovers:** the **95** generated-but-outranked (corrected numeric audit), of which
**92 are "truth loses the summed Y-ladder"** — the selector is summed-Y-ladder-INTENSITY-
primary, so a wrong (usually implausible, >50 Da off) mass-split with higher summed Y
intensity beats truth even though truth typically has MORE core-Y HITS. Fix = count
core-Y *hits* (not raw summed intensity) + peptide-dominant b/y + mass-split penalty.
**And the bigger prize:** andes top-1-selects the correct backbone mass for 280/523 but
only 97 survive @1% — closing that 280→97 gap (making correct picks survive FDR) is the
same scoring lever.
Also: the existing `ANDES_GLYCO_SELECTOR=combined` was already tried and **regresses honest
@1% 253/97/1 → 236/86/1**, and it silently flips `YINDEX` on (retention+scoring change
together) — the new `gp` selector must be a *clean* variable with `YINDEX` decoupled.
A further trap: `SELECTOR_SHORTLIST_K=24` bare-rank-prefilters candidates before reranking,
so a weak-b/y truth outside the top-24 by bare rank can never be rescued by any fused
score — audit/remove this shortlist as part of leg 2.

**Verified score definitions (both scalars already computed at collapse):**
- **ScoreP** = `rk` at `glyco_search.rs:882-884` = `score_psm(...) + psm_edge_score(...)`
  (naked-backbone b/y rank-LLR).
- **ScoreG** = the `ladder` closure `glyco_search.rs:996-1003` →
  `glycan_y_intensity` (`backbone.rs:451`) or `core_y_intensity` (`backbone.rs:416`)
  (base-peak-normalised summed Y-ladder intensity).
- **ScoreGP** = `W_P·ScoreP + W_G·ScoreG`, seed `W_P=0.65, W_G=0.35` (a glyco search engine2 verified
  `w=0.35`; RankSVM-learned later).

**Change (must enter the selector, not just PIN):** precompute `ScoreGP` into a new
`GlycoPsmKey.gp_fused_score: f32` field (`glyco_psm.rs:147-194`); add a sibling
`gp_collapse_cmp` (next to `glyco_psm.rs:79-113`) that sorts on that one stored scalar
with the legacy lexicographic order as tiebreak; wire it into BOTH collapse sites —
driver `max_by` (`glyco_search.rs:1095-1104`) and PIN `select_emitted_hits`
(`glyco_pin.rs:377-389`) — behind `ANDES_GLYCO_SELECTOR=gp`. Also emit `ScoreP`/`ScoreG`
as additive PIN columns (`glyco_pin.rs` header ~:113-131, row ~:240-281) for Percolator.

**Top-1 correctness risk (the one to watch):** driver/PIN collapse divergence — the
winner is chosen twice and must agree byte-for-byte (documented past bug). Neutralize by
(a) storing the fused scalar so the PIN *reads* it and never recomputes (it has no raw
peaks), and (b) **never per-scan-normalizing** ScoreP/ScoreG — use fixed process-constant
affine rescale + `f32::total_cmp` only (per-scan normalization is the exact class of the
past 40% FDR determinism swing). Guard with a regression test asserting
driver-winner == PIN-winner under `=gp`, plus a Step-0 byte-identical check.

---

## Leg 3 — Glycan-decoy channel makes expansion FDR-safe (guardrail)

**This is the finding that re-frames the campaign.** a glyco search engine controls FDR in 2D
(glycan-decoy + peptide-decoy). andes is Percolator-only — but the **glycan-decoy
channel is already built** (the "G3" machinery) and gated OFF:
- decoy Y-ladder `glycan_y_intensity_decoy` (`backbone.rs:527`, with passing tests at
  `:773, :840`);
- paired Label-−1 glycan-decoy PIN rows (`glyco_pin.rs:446-456`), forced label
  (`:172`), decoy `YLadderScore` (`:245`), **negated** `SialicConsistency` (`:261`),
  `glycandecoy_` accession prefix (`:298`);
- env gates `ANDES_GLYCO_DECOY` (`glyco_search.rs:277`, `andes.rs:2394`) and
  `ANDES_GLYCO_FULL_GLYCANS` (`andes.rs:2357`).

**The mechanism (why crash vs pay):** expansion crashes from **target/decoy imbalance on
the glycan axis** — more isobaric wrong-glycans inflate target-labeled glyco-noise, but
the only decoys are mass-preserving reversed *peptides* that don't model a
wrong-glycan-on-real-peptide error. The glycan-decoy channel restores the balance:
expansion pays **iff every added target glyco-hypothesis is matched 1:1 by a
glycan-axis decoy**. The additive glycan-quality features (`YLadderScore`,
`SialicConsistency`) are necessary but **inert without the decoys**.

**⚠️ CORRECTION (handoff §3):** an earlier version of this doc guessed the "full-glycan
253→119 crash" was confounded by the decoy channel being off, and proposed
`FULL_GLYCANS=1 DECOY=1` as a cheap win. **That is wrong.** The handoff refuted-list now
records BOTH: "full 2510-glycan list crashes 253→119" AND "**unified glycan-decoy
Percolator pile CRASHED — glycan-axis features too underpowered in one label pile.**" So
the decoy channel *was* run with expansion and still crashed. Root cause: glycan-decoy
rows share nearly all peptide features and only a small set of glycan features differ, so
Percolator cannot separate them. **The decoy channel is not a cheap win — it is blocked on
richer composition/Y-ladder features (see leg 2's "49 lose the ladder" cases) and belongs
LAST in the order, after the selector can rank recovered backbones.**

**Single joint Percolator run is defensible and primary:** the decoy set becomes
(peptide-decoys ∪ glycan-decoys), so one q-value covers "either axis wrong" — the native
analogue of a glyco search engine's `FDR^G + FDR^P`. Not subtracting the intersection makes it slightly
conservative (safe, never optimistic). A thin post-Percolator glycan-decoy competition is
optional, build only if measured glycan-FDP stays high.

**Hardening before trusting it (from audit 06):** `CoreYHits` on the decoy row currently
reuses the target's value (`glyco_pin.rs:247`) — recompute against the shifted ladder;
verify the 1–30 Da rung shift reliably misses at 20 ppm; audit 1:1 target:glycan-decoy
balance and 1-target-row-per-scan collapse.

---

## Decisive next step (do this FIRST, no engine code) — the P0 candidate-level audit

Not a decoy experiment. Per handoff §5-P0/§6.1, build the explicit **truth-vs-winner
competitor table** before touching engine code — it decides whether leg 2 is a small
selector fix or a real learned scorer, and prevents the `combined`-selector trap (top-1
rises while honest @1% falls).

The **outranked half already exists**: `glyco_outrank_audit.py --out per-scan.tsv` already
emits, per outranked scan, truth-vs-winner `RankScore`(ScoreP) + `YLadderScore`/`CoreYHits`/
`Y0Y1Anchor`/`SialicConsistency`(ScoreG) with per-term gaps and the 57/49 reason split.
Run it and read the per-scan TSV first.

The **missing half is the `truth_absent` decomposition** — distinguish, for the 148 absent:
*no precursor-mass DB source* vs *top-k/retention loss* vs *`SHORTLIST_K=24` loss* vs
*glycan-list gap*, and characterize by charge and mass-split (handoff: miss skews large /
high-charge, absent ~1750–1860 Da, z≥5 ~unrecovered). This needs (a) the all-hit PIN +
truth TSV on the VM and (b) likely a diff against a max-retention / precursor-mass-on run,
since "excluded by retention" is not observable from a single emitted PIN. See
`glyco_absent_audit.py` (scaffolded alongside) for the PIN-derivable part.

---

## Build sequence (reconciled with handoff §5; each step its own flag; gate @1%+decoys AND top1)

0. **P0 candidate audit (no code):** run `glyco_outrank_audit.py` on the all-hit PIN; add
   the `truth_absent` cause decomposition. Decide: small selector fix vs learned scorer.
1. **Patch toggles before any A/B (small code):** decouple `ANDES_GLYCO_SELECTOR` from the
   implicit `YINDEX` flip; add a regression asserting driver collapse and PIN collapse pick
   the same winner under every selector mode. Audit/remove `SELECTOR_SHORTLIST_K=24`.
2. **Carrier inert (no flag):** add `gp_fused_score` field + `gp_collapse_cmp` +
   precursor-mass ctx fields, all OFF. **Gate: byte-identical to 253/97.**
3. **Selector alone** (`ANDES_GLYCO_SELECTOR=gp`, `YINDEX` fixed): pool unchanged → FDR-safe;
   attacks the 57 ScoreP-fixable outranked. Ships first. **Gate: honest @1% > 253/97/1.**
   Note: the 49 "truth loses the ladder" need a stronger ScoreG (composition terms), not
   just fusion — expect fusion alone to recover the 57, not all 106.
4. **Generation alone** (`ANDES_GLYCO_PRECURSOR_MASS=1`, decoys OFF): adds the ~150
   unsearchable. Judge on top-1 presence + the audit, NOT honest @1% yet (expansion without
   separation crashes — expected). Only proceed if leg 2 proved it can rank.
5. **Joint P1+P2** (selector + precursor-mass): the payoff — new backbones enter AND get
   ranked. **Gate: honest @1% > 253/97/1 with decoys controlled.**
6. **Learn `W_G/W_P`** via RankSVM/LambdaMART on a held-out split (avoid leakage on the
   523-scan Fc3 set; engine-wide ModelStore path, not a glyco-only fork).
7. **Decoys / richer glycan features (LAST, P3):** only after 3–6 land. The glycan-decoy
   pile crashed because glycan-axis features are underpowered — add richer composition/
   Y-ladder features first, then revisit the glycan-decoy channel. Not a near-term lever.

**Constraints (unchanged):** FDR = Percolator only (never Mokapot); additive PIN
features + precomputed-into-key selector scalar only (no `score_psm` rewrite);
deterministic (no HashMap / no per-scan normalization in ordered paths); model/GBDT
changes are engine-wide; validate @1%+decoys, never top1-alone.
