# andes glyco — tracked issues (from the 2026-07-03 data audit vs published the reference engine-nglycan)

Evidence in `LESSONS.md` + memory. Data files: `scratchpad/{msf_assigned,msf_enum,andes_db_glycans}.txt`, `Fc3_r1.pepXML`, `the reference engine_glyco_targets.tsv`.

| # | Issue | Evidence | Status | Fix |
|---|---|---|---|---|
| **GI-1** | ~50% of andes glyco PSMs are **de-novo mass residuals**, not enumerated glycan compositions | 700,705 target rows: 352,475 `IsGlycanDb=1` vs **348,230 de-novo** (35,882 distinct "GlycanMass" up to 10,698 Da) | **OPEN → fixing** | Emit only enumerated-composition hits into the FDR PIN (default); de-novo = diagnostic-only |
| **GI-2** | Glyco **features do not separate target from (reversed-peptide) decoy** | target vs decoy means IDENTICAL: GlycanMass 1692.8≈1698.8, YLadder 0.043≈0.044, CoreY 0.978≈0.996 | **OPEN — partly structural** | See analysis below: only *composition-specific Y-ladder* features can discriminate a glycan on a fixed spectrum+backbone; Oxonium/CoreY/GlycanMass are spectrum/backbone-level and **cannot**. Real fix = isobaric-composition glycan decoy + stronger composition-specific features, NOT "copy decoy to all columns." |
| **GI-3** | andes's glycan **DB has coverage gaps** vs the reference engine | 27/123 (22%) the reference engine-assigned glycans absent: 8 small <892 Da (HexNAc, HexNAc2, HexNAc2Hex1-2 = paucimannose/truncated) + 19 large (1458–4028) | **OPEN → fixing** | Expand `n_glycan_list` ranges: add paucimannose/truncated (below trimannosyl core) + larger sialyl/fucosyl compositions |

## GI-2 analysis (why "extend decoy to all features" is not literally right)

On a FIXED spectrum + FIXED backbone mass, of the glyco PIN features:
- `OxoniumScore`, `NCoreOxoniumIons` — computed from the **spectrum's oxonium region** (`oxonium_gate(spec.peaks)`), independent of the candidate glycan → identical target/decoy.
- `CoreYHits`, core-Y intensity — the **trimannosyl-core ladder is glycan-INDEPENDENT** (every N-glycan shares it), anchored on the backbone → identical target/decoy.
- `GlycanMass` — the glycan-axis decoy is **isobaric by design** (same total mass so the precursor still matches) → identical by construction.
- `YLadderScore` — the **composition-specific** stepwise Y-ladder → the ONLY feature a glycan decoy can move.

⇒ The current glycan decoy (`glycan_y_intensity_decoy`: same composition, shifted interior rungs) correctly moves only `YLadderScore`. Making *more* features discriminate requires either (a) an **isobaric different-composition** decoy (so a composition-specific *oxonium-ratio* feature and the Y-ladder both differ), and/or (b) **new composition-conditioned features** (e.g. per-composition oxonium-intensity-ratio, sialic-acid diagnostic ratio). This is a research task, tracked here — NOT a one-line "copy to all columns" fix.

## Order of work
1. **GI-1** (enumerated-only PIN filter) — cheapest, removes the ~50% non-ID rows.
2. **GI-3** (DB expansion) — closes the 22% coverage gap.
3. **GI-2** (isobaric decoy + composition-conditioned features) — the hard, real lever for glycan-axis discrimination; do after GI-1/GI-3 so it's measured on a clean, enumerated, top-1-collapsed PIN.
