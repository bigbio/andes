# andes PTM-Refinement Cascade — design (X!Tandem-style, FDR-correct)

**Date:** 2026-06-17 · **Branch:** feat/gbdt-stronger-models (fresh-repo line)

**Goal:** A second-pass *refinement* search that expands the candidate space with common PTMs — but **only on the proteins the first pass already showed evidence for** — so the expensive PTM expansion is memory-feasible (global PTM/semi OOMs at 30+ GB), and the newly-found modified peptides are admitted at honest FDR by a discriminative score. This is X!Tandem's refinement model, made FDR-correct with the Crux cascade-search + MaxQuant dependent-peptide + subgroup-FDR construction.

**Hard constraint (the whole point):** do **not** violate the FDR paradigm. X!Tandem refinement is criticized for inflating FDR by re-searching a shrunken DB while reusing the full-DB null. We avoid that by construction.

**Tech stack:** Rust (`search`, `scoring`, `output`, `andes`); reuses the chimeric two-pass machinery (`run_pass2_coisolation` → `force_push`, residual scoring), the internal TDC q-gate (from `mass_calibrator`), RichIonLLR (the discriminative gate), Percolator/mokapot for downstream group-FDR. New `serde_yaml` dep for the refinement config.

---

## Why
- A discriminative score on a FIXED search space is at its ceiling (RichIonLLR flat, +0.13% on Astral); the gains come from **expanding the candidate space** (chimeric: +101%). The score is the *enabler*; the new true candidates are the gain. See `[[search-space-expansion-strategy]]`.
- **Global** PTM/semi expansion is memory-infeasible (andes materializes all candidates: semi OOM'd at 31.6 GB on a 13.5k-protein DB). The fix — proven by X!Tandem refinement, Crux cascade-search, and MaxQuant dependent-peptide — is to expand **only on the confident-protein subset** (~hundreds–low-thousands → fits).

## Prior art (grounding; cited)
- **Crux `cascade-search`** (Kertesz-Farkas, Keich, Noble, *JPR* 2015, 14:3027): ordered DB tiers; each spectrum identified at q≤α is **removed before the next tier**; FDR controlled **per tier against that tier's own paired decoys** → disjoint partition → per-tier α bounds the merged list. **Min-discoveries guard (k≈20)**: don't report a tier with too few decoys (naïve two-stage FDR fails without it).
- **MaxQuant dependent-peptide** (ModifiComb lineage): seeds from the 1%-FDR identified peptides, searches **only spectra the first pass left unidentified**, gives discovered modified peptides their **own target-decoy FDR stratum** + a localization probability.
- **Subgroup FDR** (Bogdanow/Selbach *MCP* 2016; PTMiner mass-shift-grouped FDR; mokapot grouped confidence, Fondrie & Noble 2021): a global FDR is **invalid** for modified peptides (the large unmodified population masks modified false positives). Estimate FDR **within each mod-class group, each with its own targets+decoys**.

## Architecture — the cascade (Design A: disjoint union)
One in-engine invocation (`--refine`), mirroring the existing chimeric two-pass:
1. **Pass-1:** the normal search (full DB, fully-tryptic, light mods). Compute an **internal TDC q-value** (reuse `mass_calibrator`'s q machinery) — used for SCOPING ONLY.
2. **Scope (search-space decision, NOT an FDR claim):**
   - **Confident proteins** = proteins with ≥1 Pass-1 target PSM at internal-TDC-q ≤ `refine-select-psm-fdr` (default **0.10** — permissive, coverage-driven, "proteins with peptide evidence").
   - **Refinement spectra** = the spectra Pass-1 left **unidentified** (no target PSM at the report threshold). Refine only these (cascade/DP style) — avoids re-doing Pass-1 hits and the same-spectrum-two-PSMs problem.
3. **Refinement DB** = {each confident target protein} ∪ {**its paired decoy protein**}, the decoy carried **1:1 by target-membership, regardless of the decoy's own Pass-1 score**. PTM expansion is applied **identically to targets and decoys**.
4. **Pass-2:** search the refinement spectra against the refinement DB with the PTM-rich params; `force_push` winners (like chimeric secondaries), each tagged with its expansion class.
5. **Report = disjoint union:** Pass-1 set (at 1%) ⊎ Pass-2 set (at 1% subgroup-FDR). Because the spectrum sets are disjoint and each is independently decoy-calibrated, the merged list is FDR-bounded.

## FDR design (non-negotiable)
- **Selection ≠ FDR.** The 10% scoping gate only chooses the search space; the reported FDR is re-derived inside the refinement DB. Permissiveness costs coverage/compute, not validity.
- **Subgroup FDR:** Pass-2 PSMs are FDR'd as their **own group(s)**, stratified by mod-class, each against the refinement DB's **paired decoys** — never pooled with Pass-1 unmodified PSMs. Concrete class map: `oxidation`, `deamidation`, `nterm_acetyl`, `nterm_loss` (pyro-Glu Q+E pooled), and `alkyl` (when that tier is on — carbamyl folds here) each get their OWN class (prevalent enough to clear k≈20); `metabolic`, `glycation`, non-canonical-phospho fold into an **"other/rare"** bin unless the sample makes them abundant. A class below k≈20 decoys is folded into "other" or dropped (sparse-class FDR is unreliable).
- **Symmetric expansion:** targets and decoys get the identical PTM expansion (same mods, same `refine-max-mods`) so the modified null is the right size.
- **Min-discoveries guard:** don't report a mod-class group with too few decoys to estimate N_decoy/N_target (k≈20); fold sparse groups into an "other" bin or drop.
- **Score comparability:** emit `is_refinement` (0/1), `num_mods`, and `refine_mod_class` (categorical) as **additive PIN columns** so the downstream rescorer (Percolator/mokapot `--group-column`) scores modified rows on their own scale and computes per-group q-values. (Additive-PIN is the proven-safe integration; never modify existing score columns.)
- **What we explicitly do NOT do (the X!Tandem sin):** re-search a different-size DB while reusing the Pass-1 null. Pass-2's null is always its own paired decoys in its own DB.
- **Validation is mandatory:** entrapment-FDP A/B (below) must confirm the merged 1% is truly ≤1%.

## Refinement PTM set + config (YAML) — tiered (from the 2026-06-17 PTM-selection review)
A `--refine-config <yaml>` defines the refinement spec. **The MVP ships ONE small DEFAULT tier** (the universal head of the open-search Δm distribution — present in every human-proteome tryptic run regardless of biology); everything sample-conditional is an **opt-in tier**, never on by default (each opened variable mod can cost ~16% of base IDs via search-space inflation + raise the modified-class FDR).

**DEFAULT tier (MVP — the only built-in default), 5 mods:**
```yaml
select_psm_fdr: 0.10       # confident-protein SCOPING gate (coverage, not FDR)
max_mods: 2
high_res_only: true        # gates near-isobaric arbitration (below)
mods:
  - {name: Oxidation,       delta:  15.994915, residues: [M],     location: anywhere,       accession: "UNIMOD:35", class: oxidation}
  - {name: Deamidation,     delta:   0.984016, residues: [N, Q],  location: anywhere,       accession: "UNIMOD:7",  class: deamidation}
  - {name: "Gln->pyro-Glu", delta: -17.026549, residues: [Q],     location: n_term,         accession: "UNIMOD:28", class: nterm_loss}
  - {name: "Glu->pyro-Glu", delta: -18.010565, residues: [E],     location: n_term,         accession: "UNIMOD:27", class: nterm_loss}
  - {name: Acetyl,          delta:  42.010565, residues: ["*"],   location: protein_n_term, accession: "UNIMOD:1",  class: nterm_acetyl}
```
Combinatorially cheap: the 3 terminal mods are position-restricted (one site, not every residue), so the real cost under `max_mods:2` is just Ox(M) × Deam(N/Q). (Dropped Pro-oxidation from an earlier draft: collagen-specific + Δm-collides with Met-ox under the cap.)

**Cross-tool grounding (2026-06-17 multi-tool review).** This 5-mod head is a near-exact clone of X!Tandem's *always-on* refinement chemistry — `protein, quick acetyl=yes` (protein-N-term acetyl) + `protein, quick pyrolidone=yes` (pyro-Glu Q −17.0265 / pyro-Glu E −18.0106) — plus oxidation-M, which is the universal default of every engine surveyed (X!Tandem `refine, potential modification mass` convention `15.994915@M`; MaxQuant / MSFragger / Comet default variable mod). MetaMorpheus's **G-PTM-D** curated common-mod list (deamidation, pyro-Glu, oxidation, acetyl, …) is the closest field analog to our tiered design and corroborates each of these five.
- **Deamidation N/Q is the one DEFAULT mod more aggressive than the field** — MaxQuant / MSFragger / Comet treat it as sample-specific (not a default); X!Tandem refinement convention and G-PTM-D's *artifact* bucket do include it. Keeping it is justified (andes is explicitly X!Tandem-refinement-style) but it is near-isobaric with ¹³C/+1 Da (Δ≈19 mDa) → **it is the DEFAULT mod most worth A/B-gating** on the TMT/Astral/UPS entrapment harness before locking in.
- **Carbamyl was REMOVED from DEFAULT** (was in the earlier 6-mod draft): no surveyed engine carries carbamylation in a *default* tier — X!Tandem's always-on N-term chemistry is acetyl + pyro-Glu, not carbamyl. It is a real urea/cyanate sample-prep artifact → demoted to the opt-in `alkylation`/artifact tier (below), not silently dropped.

**Opt-in tiers** (`--refine-config` selects; never stacked): `alkylation` (over-CAM K/H/N-term +57.021, propionamide C +71.037, failed-alkylation, **carbamyl N-term/K +43.0058** — the urea/cyanate artifact demoted out of DEFAULT, no field tool defaults it), `ffpe` (methylol +30.011, methylene/Schiff +12.000, methyl +14.016, formyl +27.995), `metabolic` (succination C +116.011, itaconation C +130.027), `glycation` (MG-H1 R +54.011, carboxyethyl K/R +72.021, hexose +162.053), `phospho` (STY +79.966 with −97.977 loss; only for enrichment runs).

`common-extended` (per **D. Tabb's modification-frequency analysis** — common across storage/handling regardless of sample type; the recommended *first* opt-in for richer reanalysis):
```yaml
mods:
  - {name: "Oxidation(Pro)",  delta: 15.994915, residues: [P],     location: anywhere, accession: "UNIMOD:35", class: oxidation}
  - {name: Dioxidation,       delta: 31.989829, residues: [M, W],  location: anywhere, accession: "UNIMOD:425", class: dioxidation}
  - {name: Formylation,       delta: 27.994915, residues: [K],     location: [anywhere, n_term], accession: "UNIMOD:122", class: formyl}
  - {name: Methylation,       delta: 14.015650, residues: [K, R],  location: anywhere, accession: "UNIMOD:34", class: methyl}  # see caveat
```
Caveats: **Oxidation(Pro)** shares Δm with Met-ox but a different residue (own localization, no collision) — adds oxidation-class sites, so watch the `max_mods` budget. **Methylation +14.0157 is EXACTLY isobaric with the Val↔Leu/Ile SAV (Δ=0 — high-res cannot resolve it)** → only fragment-localization + biological prior separate "methyl-K" from a sequence variant; keep methylation its own subgroup-FDR class and treat any future SAV tier as the alternative hypothesis (see Special Handling). **Dioxidation** and **Formylation** are distinct Δm (no near-isobaric trap). All four of these are independently corroborated by the cross-tool review — they sit in MetaMorpheus G-PTM-D's curated common-mod list and the X!Tandem refinement convention (oxidation M/W/P, dioxidation, formyl, methyl), so `common-extended` is field-precedented as the first opt-in, not andes-invented. (Exact David L. Tabb modification-frequency citation still to confirm with the user — neither WebSearch nor the cross-tool agent could pin the DOI; the mod list itself is sound regardless of which paper formalized it.)

**Combinatorial budget rules:** ≤~6 variable mods per tier; **never stack tiers**; declare **mutually-exclusive groups** (the N-terminal family {pyroGlu-Q, pyroGlu-E, acetyl, +carbamyl when the alkylation tier is on} compete for one position; the Cys family {succination, itaconation, propionamide, over-CAM} — one per Cys). Config guard: warn/cap if active non-terminal variable mods exceed ~8. The cascade's confident-protein scoping protects MEMORY; the `max_mods` cap + tier discipline protect per-peptide combinatorics + FDR.

## Near-isobaric disambiguation (high-res gate)
Several refinement mods are near-isobaric and MUST be gated to high-res + arbitrated at the fragment level:
- **Deamidation +0.98402 vs ¹³C isotope +1.00336** (Δ ≈ 19 mDa) — the headline case. Enumerate **both** (unmodified at isotope-error −1 vs deamidated at offset 0); **fragment-level RichIonLLR arbitrates** (deamidation shifts specific b/y; isotope shifts only the precursor). Emit as distinct PSM rows; never silently merge.
- **Acetyl +42.01057 vs Trimethyl +42.04695** (Δ ≈ 36 mDa) — both on K; needs Orbitrap accuracy + diagnostic ions (acetyl-K immonium m/z 126). Default tier carries only acetyl (protein-N-term, low collision); trimethyl stays out of the MVP.
- **Carboxyethyl +72.021 vs Propionamide +71.037** — not isobaric but aliasable; kept residue-specific (K/R vs C) and `high_res_only` if both tiers ever co-enabled (which the budget rules forbid).
`high_res_only: true` in the DEFAULT tier; any tier adding methylation/trimethyl/glycation/metabolic-Cys mods **hard-requires** `high_res_only` (unresolvable + FDR-inflating at low-res).

## Special handling — NOT plain variable mods
- **SAVs (single-amino-acid variants):** Val↔Leu/Ile (±14.0157) is **exactly isobaric** with methylation; Ala↔Ser (+15.9949) with oxidation. Modeling these as residue mods *creates* the confusion and corrupts the methylation/oxidation FDR classes. → handle as a **separate sequence/variant refinement tier** (expand the confident proteins' sequences with known SAVs, not the mod space), its own `refine_mod_class = "variant"` subgroup-FDR. **Deferred — NOT in the MVP.**
- **Metal adducts (Na +21.9819, K +37.9559):** proton-substitution / desalting artifacts, not residue PTMs (shift the precursor, not the localized backbone). → **exclude** by default (treat as unassigned; flag in QC); if ever modeled, a dedicated precursor-level `adduct` class, never pooled with acidic-residue mod classes.
- **In-source decay / neutral losses (−18.011, −17.027, pyro-Glu):** map to terminal mods (above) so the b/y ladder is consistent; they sit on the ISD-vs-real-mod boundary (PTM-Shepherd: H₂O/NH₃ losses have bimodal RT separating in-source from true) — a later RT-based flag can split the ISD subset. The **−57.021 "failed alkylation"** is decomposition, not a +Δ mod → express as *variable* (not fixed) carbamidomethyl on Cys, never a phantom negative mass.

## Localization (separate axis)
Report a **per-site localization probability** from site-determining ions (Ascore/PTM-score lineage) for modified PSMs; ≥0.75 = "localized". **Localization FLR is a separate axis from identification FDR** (a PSM can be correctly identified but mislocalized) — emit both numbers; do not conflate. Adopt the **localization-aware** approach of Yu et al. (MSFragger, *Nat Commun* 2020): score using **both the modification-shifted AND the regular (unshifted) fragment ions** so the shifted b/y ions that pin the site directly drive both identification and localization — RichIonLLR already scores per-ion, so the refinement just feeds it the shifted theoretical m/z per site hypothesis.

## RichIonLLR as the gate
Refinement PSMs are scored with the existing additive `RichIonLLR` (decoy-aware per-ion LLR) so the rescorer can separate true modified peptides (coherent shifted b/y ladders) from chance matches in the larger refinement space — the score that lets the expansion pay without FDR blowup.

## Code structure
- `crates/search/src/refinement.rs` **(new):** the refinement pass — collect confident proteins (internal-TDC-q), build the refinement `SearchIndex` (confident targets + paired decoys via the existing decoy machinery), run Pass-2 on unidentified spectra, tag PSMs with expansion class. Generalizes the `run_pass2_coisolation` pattern (`match_engine.rs`).
- `crates/search/src/refine_config.rs` **(new):** `serde_yaml` parse of `--refine-config` → refinement mods + params.
- `crates/search/src/search_params.rs`: refinement params (or a `RefineConfig` field).
- `crates/output/src/pin.rs` + `crates/search/src/psm.rs`: additive `is_refinement` / `num_mods` / `refine_mod_class` columns (+ localization-probability column).
- `crates/andes/src/bin/andes.rs`: `--refine`, `--refine-config`, `--refine-select-psm-fdr` (0.10), `--refine-max-mods` (2), `--refine-high-res-only`.
- Reuse: `mass_calibrator` internal-TDC-q; `decoy.rs` paired-decoy generation; `coisolation.rs`/`match_engine.rs` two-pass `force_push`/residual; `RichIonLLR`.

## Validation gates
- **Phase 0 — entrapment-FDP is mandatory:** base (no refine) vs base+refine on Astral entrapment; report PSMs@1% + **true entrapment-FDP** for the unmodified group AND the modified group separately. The merged-list true-FDP must be ≤1% — the only proof the construction didn't inflate.
- Per-mod-class true-FDP must each be ≤ target (subgroup validity).
- Net **new modified IDs** at flat FDP = the win (the thing global semi couldn't even run).
- Sanity: discovered Δm histogram peaks at the configured masses (real chemistry, not noise).

## Phasing
- **MVP (this spec):** PTM refinement (the 5-mod DEFAULT tier: oxidation M, deamidation N/Q, pyro-Glu Q/E, protein-N-term acetyl — a near-exact clone of X!Tandem's always-on refinement chemistry; carbamyl demoted to the opt-in alkylation tier) via the cascade + subgroup-FDR + group-tag PIN features + entrapment validation.
- **Later (same infra):** semi-tryptic and point-mutation refinement tiers (X!Tandem also refines these); open/delta-mass refinement — blueprint = MSFragger's **localization-aware open search** (Yu et al. 2020): a fragment-ion index over **both shifted and regular** ions + fast mass calibration; proper in-engine subgroup-FDR if we stop relying on mokapot.

## Non-goals (deferred)
- Full open-search (arbitrary Δm) — needs the fragment-ion index (separate project).
- In-engine FDR (andes stays a PIN emitter; group-FDR via mokapot `--group-column` downstream for the MVP).
- Refining the *full* spectrum set (we refine only Pass-1-unidentified spectra — cascade/DP style).

## Adjacent directions (researched 2026-06-17 alongside this spec — their OWN specs, NOT in this MVP)
Two pieces of the "expand the search smartly" theme were researched in the same session but are orthogonal to PTM refinement and each warrant a separate spec:
- **Mass calibration / parameter optimization (MaxQuant + MSFragger).** andes today does precursor-only, single global median-ppm calibration (`mass_calibrator.rs` / `precursor_cal.rs`, TDC-q≤1% gated). Peers add (a) **fragment-mass calibration**, (b) **position-dependent error models** (MSFragger: nonparametric RT×m/z grid; MaxQuant: nonlinear m/z/RT/log-intensity, first-search→recalibrate→main-search), (c) **tolerance auto-optimization** (MSFragger 6×5 frag-tol × top-N grid) — together a reported ~5–10% PSM gain. Highest-ROI andes gap = **a fragment-mass shift learned in the existing pre-pass** (cheap, additive, A/B-gated on Astral+TMT). This is a search-quality multiplier that compounds with the cascade; track as its own spec.
- **Glycopeptide refinement (StructGP-informed).** For the parked glyco path (Unimod 393 / per-class `loss_class` IonType), the reusable, low-risk ideas from StructGP's extracted code are: **oxonium triage** (138.055 + 204.087 essential-ion gate to flag glyco spectra), **Y-ion-ladder coverage as additive PIN features** (peptide+stepwise-glycan-loss, maps onto `loss_class`), and a **mass-shifted-glycan decoy** (precursor +20–30 Da) for an honest *glycan-level* FDP axis. Its N-glyco biosynthetic-grammar structure engine is not worth porting. Stays parked behind the campaign; captured here so the analysis isn't lost.

## References
- Kertesz-Farkas, Keich, Noble, "Tandem MS Identification via Cascaded Search," *JPR* 2015, 14:3027 (DOI 10.1021/pr501173s); Crux cascade-search docs (crux.ms).
- Tyanova, Temu, Cox, "MaxQuant computational platform," *Nat. Protoc.* 2016 (dependent-peptide); Savitski et al. ModifiComb, *MCP* 2006.
- Bogdanow, Zauber, Selbach, "Systematic Errors … by Modified Peptides," *MCP* 2016 (subgroup FDR); Fondrie & Noble, "mokapot," *JPR* 2021 (grouped confidence); Freestone, Noble, Keich, "Group-walk," *Bioinformatics* 2022.
- Yu, Teo, Kong, Haynes, Avtonomov, Geiszler, Nesvizhskii, "Identification of modified peptides using localization-aware open search," *Nat Commun* 2020, 11:4065, DOI 10.1038/s41467-020-17921-y (PMC7426425) — MSFragger localization-aware open search: index shifted + regular fragment ions, mass calibration, modified-peptide localization. https://www.nature.com/articles/s41467-020-17921-y
- *Advances in Chemical Proteomics* (X. Yao, ed.), Elsevier 2021, ISBN 978-0-12-821433-6 — Ch. 4 (chemical-proteomics / PTM methods). https://www.sciencedirect.com/science/chapter/edited-volume/abs/pii/B9780128214336000040
- Open-search PTM prevalence + selection (from the 2026-06-17 PTM-selection review): PTM-Shepherd (*MCP* 2021, PMC7950090), Chick et al. (*Nat Biotechnol* 2015), PTMiner mass-shift-grouped FDR (PMC10166509).
- Cross-tool refinement/default mod lists (2026-06-17 multi-tool review): **X!Tandem** API (thegpm.org/TANDEM/api: `pqa.html` quick acetyl, `pqp.html` quick pyrolidone, `refpmm.html` refine potential modification mass), Craig & Beavis, *Bioinformatics* 2004, 20:1466; **MetaMorpheus G-PTM-D** Solntsev, Shortreed, Frey, Smith, *JPR* 2018, 17:1844 (curated common-mod list, closest field analog to our tiers); MaxQuant / MSFragger / Comet default variable-mod sets (Oxidation M + Acetyl protein-N-term). Finding: the 5-mod DEFAULT tier clones X!Tandem's always-on chemistry; deamidation is the one above-field-default mod (A/B-gate); carbamyl has no default-tier precedent (demoted).
- Mass calibration (adjacent direction): Cox & Mann, *Nat Biotechnol* 2008, 26:1367 (MaxQuant individualized ppb recalibration); Cox et al., "Andromeda," *JPR* 2011, 10:1794; Yu et al. 2020 (above; MSFragger fast mass calibration + parameter optimization).
