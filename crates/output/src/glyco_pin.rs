// Glyco PIN writer for andes.
//
// Writes a Percolator-consumable `.glyco.pin` file for glycopeptide PSMs
// produced by `search::glyco_search::glyco_search_run`.
//
// Column layout: standard PIN columns (from `psm_feature_values`) +
// glyco-specific columns appended AFTER `RankScoreFloat`, BEFORE
// `Peptide`/`Proteins`:
//   OxoniumScore, NCoreOxoniumIons, YLadderScore, CoreYHits,
//   GlycanMass, IsGlycanDb
//
// The `Peptide` column appends a glycan composition tag for DB hits:
//   `[HexNAc{n}Hex{m}Fuc{f}NeuAc{a}NeuGc{g}@N{site}]`
// where `@N{site}` is the 1-based position of the sequon Asn in the backbone
// when the backbone has exactly one N-X-S/T, and `@N?` when it has zero or
// several (andes does not localize between sequons by default).

use std::io::{self, BufWriter, Write};
use std::path::Path;

use model::mass::{ISOTOPE, PROTON};
use search::candidate_gen::Candidate;
use search::search_index::SearchIndex;
use search::search_params::SearchParams;
use model::spectrum::Spectrum;

use andes_glyco::glyco_psm::{
    glyco_gp_fused_score, GLYCO_GP_J_DEFAULT, GLYCO_GP_K_DEFAULT,
};
use andes_glyco::hybrid::Source;

use crate::pin::{psm_feature_values, FeatureFmt};
use crate::percolator_enz::{count_internal_enzymatic, is_enzymatic_boundary};
use crate::row_context::RowContext;
use search::glyco_search::{FullGlycoPsm, GlycoSpectrumResult};

// ── header ────────────────────────────────────────────────────────────────────

fn write_glyco_header<W: Write>(
    writer: &mut W,
    min_charge: u8,
    max_charge: u8,
) -> io::Result<()> {
    let mut cols: Vec<String> = vec![
        "SpecId".to_string(),
        "Label".to_string(),
        "ScanNr".to_string(),
        "ExpMass".to_string(),
        "CalcMass".to_string(),
        "mass".to_string(),
        "RankScore".to_string(),
        "isotope_error".to_string(),
        "peplen".to_string(),
        "dm".to_string(),
        "absdm".to_string(),
    ];

    for c in min_charge..=max_charge {
        cols.push(format!("charge{}", c));
    }
    cols.push("chargeHi".to_string()); // overflow flag for z > max_charge (z6/z7)

    cols.extend_from_slice(&[
        "enzN".to_string(),
        "enzC".to_string(),
        "enzInt".to_string(),
        "NumMatchedMainIons".to_string(),
        "longest_b".to_string(),
        "longest_y".to_string(),
        "longest_y_pct".to_string(),
        "ExplainedIonCurrentRatio".to_string(),
        "NTermIonCurrentRatio".to_string(),
        "CTermIonCurrentRatio".to_string(),
        "MS2IonCurrent".to_string(),
        "IsolationWindowEfficiency".to_string(),
        "MeanErrorTop7".to_string(),
        "StdevErrorTop7".to_string(),
        "MeanRelErrorTop7".to_string(),
        "StdevRelErrorTop7".to_string(),
        "matchedIonRatio".to_string(),
        "EdgeScore".to_string(),
        "PrecursorIsotopeKL".to_string(),
        "PrecursorSNR".to_string(),
        "DeltaRankScore".to_string(),
        "TailorScore".to_string(),
        "PpmGaussianScore".to_string(),
        "NeutralLossIonCount".to_string(),
        "LongestComplementaryLadder".to_string(),
        "ComplementaryIonBalance".to_string(),
        "MeanMatchedIntensityRank".to_string(),
        "DoublyChargedMatchedIonCount".to_string(),
        "UniqueMatchFraction".to_string(),
        "ChanceMatchSurprise".to_string(),
        "IntensitySignal".to_string(),
        "FragPredExplained".to_string(),
        "FragPredChanceLLR".to_string(),
        "FragTopKObserved".to_string(),
        "RichIonLLR".to_string(),
        "IsRefinement".to_string(),
        "NumMods".to_string(),
        "RefinementModClass".to_string(),
        "ModSiteShiftedMatched".to_string(),
        "ModSiteShiftedFrac".to_string(),
        "ModSiteIntensFrac".to_string(),
        "ModSiteLocalized".to_string(),
        "ModSiteDetCount".to_string(),
        "MassCompetitionEvidence".to_string(),
        "CandidateRankEntropy".to_string(),
        "ListwiseScoreGap".to_string(),
        "RawScore".to_string(),
        "RawScoreCal".to_string(),
        "RankScoreFloat".to_string(),
        // Engine-wide retention-time columns (ADDITIVE; mirror the shared
        // `psm_feature_values` tail so the glyco PIN carries the SAME RT
        // features as the standard PIN). 0.0 (neutral) without observed RT.
        "DeltaRT".to_string(),
        "AbsDeltaRT".to_string(),
        "DeltaRTNorm".to_string(),
        // Glyco-specific columns (appended after the RT block, before Peptide/Proteins)
        "OxoniumScore".to_string(),
        "NCoreOxoniumIons".to_string(),
        "YLadderScore".to_string(),
        // Completeness of the assigned composition's OWN predicted Y ladder,
        // in [0,1]. Unlike YLadderScore (an unnormalised intensity sum that
        // grows with glycan size), predicting rungs you cannot support LOWERS
        // this. Targets the measured failing stage: choosing the mass split.
        "YHitFrac".to_string(),
        "PartialGlycanBY".to_string(), // idea B: sequence-specific partial-glycan b/y evidence
        "CoreYHits".to_string(),
        "GlycanMass".to_string(),
        "IsGlycanDb".to_string(),
        "Y0Y1Anchor".to_string(), // G2 peptide-mass anchor (additive, peptide-discriminating)
        "SialicConsistency".to_string(), // GI-2 composition-conditioned sialic-oxonium (additive)
        "CzHyperscore".to_string(), // ETD c/z backbone hyperscore (additive; ETD/AI-ETD only, else 0)
        "CzIntensity".to_string(), // ETD c/z matched-intensity fraction (additive; ETD/AI-ETD only, else 0)
        "CzExplained".to_string(), // ETD c/z analytical explained-intensity LLR (additive; graded-model de-risk; ETD only, else 0)
        "CzChanceLlr".to_string(), // ETD c/z prior-weighted chance-match surprise (additive; graded-model de-risk; ETD only, else 0)
        "IsTransferred".to_string(),        // cross-spectrum transfer provenance (additive)
        "TransferGraphSupport".to_string(), // # corroborating co-eluting siblings
        "TransferSeedScore".to_string(),    // donor seed Pass-1 discriminant
        "TransferRTDelta".to_string(),      // |RT(acceptor)-RT(seed)| seconds
        "TransferUngated".to_string(),      // 1 = RT gate skipped (no RT)
        // Glyco RT rank (ADDITIVE, appended LAST): within-scan rank of this
        // candidate's AbsDeltaRT among competing glycoforms, normalized to (0,1].
        // The isobaric-glycoform-disambiguation feature. 0.0 (neutral)
        // when the scan has <2 candidate hits or RT is unavailable.
        "DeltaRTRank".to_string(),
        // Isobaric-glycan RT margin (ADDITIVE): how much better the assigned
        // composition fits the observed RT than the best near-isobaric alternative
        // on the same backbone. >0 ⇒ composition RT-consistent. 0.0
        // when RT unavailable / no isobaric alternative / de-novo hit.
        "IsobaricRTMargin".to_string(),
        // Terminal columns
        "Peptide".to_string(),
        "Proteins".to_string(),
    ]);

    writeln!(writer, "{}", cols.join("\t"))
}

// ── per-row writer ─────────────────────────────────────────────────────────────

#[allow(clippy::too_many_arguments)]
fn write_glyco_psm_row<W: Write>(
    writer: &mut W,
    spec: &Spectrum,
    hit: &FullGlycoPsm,
    cand: &Candidate,
    row_idx: usize,
    min_charge: u8,
    max_charge: u8,
    candidates: &[Candidate],
    search_index: &SearchIndex,
    params: &SearchParams,
    // All of THIS scan's candidate hits (for the within-scan DeltaRTRank).
    scan_hits: &[FullGlycoPsm],
    // When Some, this row is a GLYCAN-AXIS decoy: force Label -1, use the decoy
    // Y-ladder score, and mark the SpecId/accession so it is traceable and
    // Percolator treats it as a decoy (G3 2D-FDR). None = ordinary row.
    glycan_decoy_override: Option<f32>,
) -> io::Result<()> {
    let psm = &hit.psm;
    let key = &hit.glycan_key;

    let is_glycan_decoy = glycan_decoy_override.is_some();
    let ctx = RowContext::new(spec, cand, search_index);
    let scan = ctx.scan;
    let gd_suffix = if is_glycan_decoy { "_gd" } else { "" };
    let spec_id = format!("{}_glyco_{}_{}{}", ctx.spec_id, scan, row_idx + 1, gd_suffix);

    // A glycan-decoy is a decoy on the glycan axis regardless of the peptide's
    // own target/decoy status.
    let label: i32 = if is_glycan_decoy || cand.is_decoy { -1 } else { 1 };

    let charge = psm.charge_used as f64;
    let precursor_mz = spec.precursor_mz;
    let exp_mass = precursor_mz * charge - charge * PROTON;

    // CalcMass: the INTACT glycopeptide neutral mass (peptide backbone +
    // glycan), NOT the bare backbone alone. The b/y scoring against Percolator
    // is bare-backbone (see module doc), but CalcMass/dm/absdm here are
    // compared against the observed PRECURSOR mass, which includes the
    // glycan. Omitting `glycan_mass` (0.0 for de-novo hits without a resolved
    // composition) previously left `dm`/`absdm` off by ~the glycan mass,
    // corrupting the mass-error PIN features Percolator uses for q-values.
    let calc_mass = cand.peptide.mass() + key.glycan_mass;
    let mass = exp_mass;

    let raw_score = psm.score.round() as i32;
    let isotope_error: i32 = psm.isotope_offset as i32;
    let peplen = cand.peptide.length() + 2;

    let theo_mz = calc_mass / charge + PROTON;
    let adjusted_exp_mz = precursor_mz - ISOTOPE * (isotope_error as f64) / charge;
    let dm = adjusted_exp_mz - theo_mz;
    let absdm = dm.abs();

    // SpecId, Label, ScanNr, ExpMass, CalcMass, mass
    write!(writer, "{}\t{}\t{}\t", spec_id, label, scan)?;
    write_double(writer, exp_mass)?;
    writer.write_all(b"\t")?;
    write_double(writer, calc_mass)?;
    writer.write_all(b"\t")?;
    write_double(writer, mass)?;

    // RankScore, isotope_error, peplen, dm, absdm
    write!(writer, "\t{}\t{}\t{}\t", raw_score, isotope_error, peplen)?;
    write_double(writer, dm)?;
    writer.write_all(b"\t")?;
    write_double(writer, absdm)?;

    // Charge one-hot
    for c in min_charge..=max_charge {
        let flag: u8 = if c == psm.charge_used { b'1' } else { b'0' };
        writer.write_all(&[b'\t', flag])?;
    }
    // chargeHi overflow: 1 when the glycopeptide's charge exceeds the one-hot range
    // (glyco tries the spectrum's own charge, which reaches z6/z7 on high-charge
    // glycopeptides — otherwise those get an ALL-ZERO charge vector and Percolator
    // cannot isolate them from lower charges with different target/decoy separation).
    let charge_hi: u8 = if psm.charge_used > max_charge { b'1' } else { b'0' };
    writer.write_all(&[b'\t', charge_hi])?;

    // enzN, enzC, enzInt
    let residues: Vec<u8> = cand.peptide.residues.iter().map(|aa| aa.residue).collect();
    let first = residues.first().copied().unwrap_or(b'-');
    let last = residues.last().copied().unwrap_or(b'-');
    let enz_n: u8 = is_enzymatic_boundary(cand.peptide.pre, first, params.enzyme) as u8;
    let enz_c: u8 = is_enzymatic_boundary(last, cand.peptide.post, params.enzyme) as u8;
    let enz_int = count_internal_enzymatic(&residues, params.enzyme);
    write!(writer, "\t{}\t{}\t{}", enz_n, enz_c, enz_int)?;

    // Standard feature block (reuse psm_feature_values; skip RankScore + isotope_error
    // which were already written above in the leading block).
    let features = psm_feature_values(psm, 1);
    debug_assert_eq!(features[0].0, "RankScore");
    debug_assert_eq!(features[1].0, "isotope_error");
    for &(_name, value, fmt) in &features[2..] {
        writer.write_all(b"\t")?;
        match fmt {
            FeatureFmt::Int => write!(writer, "{}", value.round() as i64)?,
            FeatureFmt::Fixed6 => write!(writer, "{:.6}", value)?,
            FeatureFmt::Double => write_double(writer, value)?,
        }
    }

    // Glyco-specific columns (8, incl. Y0Y1Anchor and YHitFrac).
    let is_glycan_db: u8 = if key.glycan_source == Source::Db { 1 } else { 0 };
    write_double_tab(writer, key.oxonium_summed_frac as f64)?;
    write!(writer, "\t{}", key.n_core_oxonium_ions)?;
    // YLadderScore: the decoy score for a glycan-decoy row, else the target score.
    let y_ladder = glycan_decoy_override.unwrap_or(key.y_ladder_intensity_score);
    write_double_tab(writer, y_ladder as f64)?;
    // YHitFrac: on a glycan-decoy row report the decoy twin's completeness, so the
    // row is a like-for-like control of its target rather than a copy of it.
    let y_frac = if glycan_decoy_override.is_some() {
        key.y_hit_frac_decoy
    } else {
        key.y_hit_frac
    };
    write_double_tab(writer, y_frac as f64)?;
    // PartialGlycanBY (idea B): sequence-specific partial-glycan b/y evidence.
    write_double_tab(writer, key.partial_glycan_by as f64)?;
    write!(writer, "\t{}", key.core_y_hits)?;
    write_double_tab(writer, key.glycan_mass)?;
    write!(writer, "\t{}", is_glycan_db)?;
    // G2 Y0/Y1 anchor (additive, peptide-mass-conditioned). Same value on a
    // glycan-decoy row (it is a peptide-axis feature, independent of the glycan).
    write_double_tab(writer, key.y0y1_anchor_score as f64)?;
    // GI-2 sialic consistency (composition-conditioned). On a glycan-decoy row we
    // emit the NEGATED value: the decoy models an isobaric glycan whose sialic
    // claim is flipped, and flipping each `±obs` term negates the whole score
    // (`a + g` with `a = if neuac>0 {+obs} else {-obs}`). So a sialylated spectrum
    // SUPPORTS the target's sialic claim (+high) and CONTRADICTS the decoy's
    // (−high) — giving the glycan axis a real composition discriminator alongside
    // the shifted Y-ladder (GI-2 part 2). For non-sialylated spectra obs≈0 → ≈0
    // on both (no harm). Additive PIN feature; peptide axis unchanged.
    let sialic = if is_glycan_decoy {
        -key.sialic_consistency
    } else {
        key.sialic_consistency
    };
    write_double_tab(writer, sialic as f64)?;

    // ETD c/z backbone hyperscore (additive; ETD/AI-ETD spectra only, else 0.0).
    // A peptide-axis feature (backbone c/z ladder), so a glycan-decoy row emits the
    // SAME value as its paired target — like Y0Y1Anchor above.
    write_double_tab(writer, key.cz_hyperscore as f64)?;
    write_double_tab(writer, key.cz_intensity as f64)?;
    // Discriminative c/z STRUCTURE features (additive; gated to real values by
    // ANDES_GLYCO_CZ_STRUCT at compute time, else 0.0 = ignored by Percolator).
    write_double_tab(writer, key.cz_explained as f64)?;
    write_double_tab(writer, key.cz_chance_llr as f64)?;

    // Transfer columns (additive; inert 0 for native candidates). Bools mirror
    // the `is_glycan_db` 1/0 idiom above; numerics use write_double_tab.
    write!(writer, "\t{}", if key.is_transferred { 1 } else { 0 })?;
    write_double_tab(writer, key.transfer_graph_support as f64)?;
    write_double_tab(writer, key.transfer_seed_score as f64)?;
    write_double_tab(writer, key.transfer_rt_delta as f64)?;
    write!(writer, "\t{}", if key.transfer_ungated { 1 } else { 0 })?;

    // Glyco RT rank (ADDITIVE, LAST): within-scan AbsDeltaRT rank of THIS hit
    // among the scan's competing glycoforms, normalized to (0,1]. Reuses the
    // `abs_delta_rt` populated by `populate_glyco_rt_features`. A glycan-decoy
    // row shares the row_idx of its paired target, so it takes that hit's rank
    // (its AbsDeltaRT is the same peptide-axis RT delta). 0.0 when <2 candidates.
    write_double_tab(writer, crate::glyco_rt::delta_rt_rank(scan_hits, row_idx))?;

    // Isobaric-glycan RT margin (ADDITIVE): from `psm.features`, populated by
    // `populate_glyco_rt_features`. >0 ⇒ the assigned composition fits the observed
    // RT better than the best near-isobaric alternative. A glycan-decoy
    // row shares its paired target's psm, so it inherits the same value.
    write_double_tab(writer, psm.features.isobaric_rt_margin as f64)?;

    // Peptide column: backbone sequence + optional glycan tag.
    //
    // The tag carries the glycosite as `@N<pos>` (1-based position of the sequon
    // Asn in the backbone) ONLY when the backbone has exactly one N-X-S/T sequon,
    // in which case the site is a fact about the sequence and needs no evidence.
    // With two or more sequons it emits `@N?`: andes does not localize between
    // them by default (the site used internally for c/z scoring is the first
    // sequon, a positional convention, not a localization), and printing that
    // heuristic choice as if it were an assignment would be a false claim. Runs
    // with `--glyco-cz-multisite` pick the site by c/z evidence, but that is
    // opt-in and its choice is not threaded here, so `@N?` stays.
    //
    // Target and glycan-decoy rows are tagged identically, so this is symmetric
    // on the FDR axis; and PSM-level q-values do not depend on the Peptide string.
    let glycan_tag = match &key.glycan {
        Some(g) => glycan_tag_with_site(g, &residues),
        None => String::new(),
    };
    write!(writer, "\t{}{}", cand.peptide, glycan_tag)?;

    // Proteins column(s). A glycan-decoy carries a decoy-prefixed accession so
    // Percolator (and downstream protein grouping) recognise it as a decoy.
    for &cidx in &psm.candidate_idxs {
        let cand_for_acc = &candidates[cidx as usize];
        let accession = crate::row_context::resolve_accession(cand_for_acc, search_index);
        if is_glycan_decoy {
            write!(writer, "\tglycandecoy_{}", accession)?;
        } else {
            write!(writer, "\t{}", accession)?;
        }
    }
    writeln!(writer)
}

/// Render the Peptide-column glycan tag, including the glycosite when the
/// backbone determines it. See the call site for why an ambiguous backbone
/// emits `@N?` rather than andes's internal first-sequon convention.
fn glycan_tag_with_site(g: &andes_glyco::glycan_db::GlycanComp, residues: &[u8]) -> String {
    let site = match andes_glyco::sequon::all_nxst_sites(residues).as_slice() {
        [only] => format!("@N{}", only + 1),
        _ => "@N?".to_string(),
    };
    format!(
        "[HexNAc{}Hex{}Fuc{}NeuAc{}NeuGc{}{}]",
        g.hexnac, g.hex, g.fuc, g.neuac, g.neugc, site
    )
}

// ── public API ────────────────────────────────────────────────────────────────

/// Write a glyco PIN file from the output of `glyco_search_run`.
///
/// Standard PIN columns are written using `psm_feature_values`; six
/// glyco-specific columns are appended after `RankScoreFloat`:
///   OxoniumScore, NCoreOxoniumIons, YLadderScore, CoreYHits,
///   GlycanMass, IsGlycanDb
///
/// `Peptide` column appends `[HexNAc{n}Hex{m}...]` for DB glycan hits.
#[allow(clippy::too_many_arguments)]
pub fn write_glyco_pin(
    path: &Path,
    spectra: &[Spectrum],
    results: &[GlycoSpectrumResult],
    candidates: &[Candidate],
    params: &SearchParams,
    search_index: &SearchIndex,
    emit_glycan_decoy: bool,
    debug: bool,
) -> io::Result<()> {
    let file = std::fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    write_glyco_pin_to(
        &mut writer,
        spectra,
        results,
        candidates,
        params,
        search_index,
        emit_glycan_decoy,
        debug,
    )
}

/// Select which hits of ONE scan to emit into the FDR PIN.
///
/// `collapse=true` (the honest default) keeps only the single best PSM — the
/// top-1-per-scan TDC winner, by `rank_score` then `YLadderScore`, lowest hit
/// index breaking full ties. This is REQUIRED for a valid per-scan FDR:
/// Percolator does NOT do target/decoy competition within a scan, so emitting the
/// ~4.5 correlated (peptide×glycan) alternatives per scan fabricates the
/// target/decoy balance and produces fake recovery (docs/plans/glyco/LESSONS.md
/// L1 — collapsing g4off's 1.02M rows→5,954 changed 155 "recovered" → 0 @1% FDR).
/// The non-glyco path enforces the same one-PSM-per-scan rule
/// (`top_n_psms_per_spectrum = 1`). `collapse=false` dumps all hits — DIAGNOSTIC
/// ONLY, never for FDR.
pub(crate) fn select_emitted_hits(
    hits: &[FullGlycoPsm],
    collapse: bool,
    enumerated_only: bool,
) -> Vec<usize> {
    // GI-1: a glyco IDENTIFICATION requires an ENUMERATED glycan composition. The
    // de-novo branch (glycan = None) reports a bare precursor−backbone mass
    // RESIDUAL that is not in the glycan DB — ~half of raw glyco rows. Those are
    // not glycopeptide IDs and must not enter the FDR PIN (docs/plans/glyco/ISSUES.md
    // GI-1). Filter them out first, THEN collapse.
    if !collapse {
        // Diagnostic dump: all hits (optionally enumerated-filtered).
        return (0..hits.len())
            .filter(|&i| !enumerated_only || hits[i].glycan_key.glycan.is_some())
            .collect();
    }
    if hits.is_empty() {
        return Vec::new();
    }
    // Collapse over ALL hits first — the TDC winner is the best rank_score
    // REGARDLESS of enumeration (Codex: filtering first would promote a
    // 2nd-place enumerated hit when the real winner is de-novo, inflating IDs).
    // Ordering is the `gp` fused selector `rank + K·ladder + J·core_y`, the SHARED
    // single source of truth with the driver's pre-feature reduction: it MUST match
    // glyco_search's gp collapse exactly (same rank/ladder/core-Y/k/j) or the scored
    // winner and the written row diverge (the collapse-parity bug).
    // The PIN-side collapse is over a SINGLE hit under the honest collapse path (the
    // driver already reduced the scan), so these weights are never decisive here; use
    // the default constants (the driver's --glyco-gp-* flags are the real authority).
    let gp_k = GLYCO_GP_K_DEFAULT;
    let gp_j = GLYCO_GP_J_DEFAULT;
    // v2 hyperscore term: the driver (glyco_search) is the SOLE authority for the gp
    // winner — under the honest `collapse=true` path it already reduced each scan to a
    // SINGLE hit before feature extraction, so this max_by is over one element and its
    // score is never decisive. The hyperscore needs raw peaks the PIN writer does not
    // hold, so pass 0.0 here; it cannot change any real output. (rank/ladder/core-Y are
    // stored fields, so the rest of the ordering stays consistent with the driver.)
    let winner = (0..hits.len())
        .max_by(|&a, &b| {
            let sa = glyco_gp_fused_score(
                hits[a].psm.rank_score,
                hits[a].glycan_key.y_ladder_intensity_score,
                hits[a].glycan_key.core_y_hits as f32,
                0.0,
                gp_k,
                gp_j,
                0.0,
            );
            let sb = glyco_gp_fused_score(
                hits[b].psm.rank_score,
                hits[b].glycan_key.y_ladder_intensity_score,
                hits[b].glycan_key.core_y_hits as f32,
                0.0,
                gp_k,
                gp_j,
                0.0,
            );
            sa.total_cmp(&sb).then(b.cmp(&a))
        })
        .expect("non-empty");
    // Emit only if the scan's actual winner has an enumerated glycan; if the
    // winner is de-novo, the scan has NO enumerated ID → drop it by default.
    // Opt-in ANDES_GLYCO_ENUM_FALLBACK mirrors the driver's collapse fallback:
    // emit the best-scoring enumerated hit instead of discarding the scan. Under
    // the honest collapse path the driver already reduced the scan to a single
    // hit, so this is a safety mirror to keep the two collapse sites consistent.
    if enumerated_only && hits[winner].glycan_key.glycan.is_none() {
        // Deliberately OFF here, even though the driver's collapse fallback is ON. The
        // asymmetry is not a split-brain default, it is the FDR boundary: the driver may
        // promote an enumerated candidate while CHOOSING a winner, but once a de-novo
        // candidate has won target-decoy competition the scan has no enumerated
        // identification, and promoting the losing runner-up would add a target that did
        // not win. See `enumerated_only_drops_de_novo_hits`, which asserts exactly this.
        let enum_fallback = false;
        if !enum_fallback {
            return Vec::new();
        }
        return (0..hits.len())
            .filter(|&i| hits[i].glycan_key.glycan.is_some())
            .max_by(|&a, &b| {
                let sa = glyco_gp_fused_score(
                    hits[a].psm.rank_score,
                    hits[a].glycan_key.y_ladder_intensity_score,
                    hits[a].glycan_key.core_y_hits as f32,
                    0.0,
                    gp_k,
                    gp_j,
                    0.0,
                );
                let sb = glyco_gp_fused_score(
                    hits[b].psm.rank_score,
                    hits[b].glycan_key.y_ladder_intensity_score,
                    hits[b].glycan_key.core_y_hits as f32,
                    0.0,
                    gp_k,
                    gp_j,
                    0.0,
                );
                sa.total_cmp(&sb).then(b.cmp(&a))
            })
            .map(|w| vec![w])
            .unwrap_or_default();
    }
    vec![winner]
}

/// Write glyco PIN to an arbitrary writer (useful for testing).
#[allow(clippy::too_many_arguments)]
pub fn write_glyco_pin_to<W: Write>(
    writer: &mut W,
    spectra: &[Spectrum],
    results: &[GlycoSpectrumResult],
    candidates: &[Candidate],
    params: &SearchParams,
    search_index: &SearchIndex,
    emit_glycan_decoy: bool,
    debug: bool,
) -> io::Result<()> {
    let min_charge = *params.charge_range.start();
    let max_charge = *params.charge_range.end();

    write_glyco_header(writer, min_charge, max_charge)?;

    // Top-1-per-scan collapse + enumerated-only is the honest default (see
    // `select_emitted_hits`). `--debug-glyco` restores the full multi-row / de-novo
    // dump for DIAGNOSTICS ONLY (that PIN must never be fed to an FDR tool). MUST
    // mirror the driver's `GlycoConfig.debug` so the kept row == the driver's hit.
    let collapse = !debug;
    let enumerated_only = !debug;

    let mut row_count = 0usize;
    for result in results {
        let spec_idx = result.spectrum_idx;
        if spec_idx >= spectra.len() {
            continue;
        }
        let spec = &spectra[spec_idx];
        for hit_idx in select_emitted_hits(&result.hits, collapse, enumerated_only) {
            let hit = &result.hits[hit_idx];
            let cand_idx = hit.psm.primary_candidate_idx() as usize;
            if cand_idx >= candidates.len() {
                continue;
            }
            let cand = &candidates[cand_idx];
            write_glyco_psm_row(
                writer, spec, hit, cand, hit_idx, min_charge, max_charge, candidates,
                search_index, params, &result.hits, None,
            )?;
            row_count += 1;

            // G3: emit a paired glycan-axis decoy for TARGET-peptide PSMs that
            // have a resolved glycan composition. Skip de-novo (no composition →
            // decoy score 0) and peptide-decoy rows (they are already decoys; a
            // glycan-decoy of a decoy peptide adds no glycan-axis signal).
            if emit_glycan_decoy && !cand.is_decoy && hit.glycan_key.glycan.is_some() {
                write_glyco_psm_row(
                    writer, spec, hit, cand, hit_idx, min_charge, max_charge, candidates,
                    search_index, params, &result.hits, Some(hit.glycan_key.y_ladder_decoy_score),
                )?;
                row_count += 1;
            }
        }
    }

    writer.flush()?;
    // Return row count via a side-channel eprintln so callers can report it.
    eprintln!(
        "[glyco-pin] wrote {} PSM rows ({})",
        row_count,
        if collapse {
            "top-1-per-scan collapse — FDR-valid"
        } else {
            "ALL hits — DIAGNOSTIC ONLY, not FDR-valid"
        }
    );
    Ok(())
}

// ── helpers ───────────────────────────────────────────────────────────────────

/// Write a tab then a double in %.6g style.
fn write_double_tab<W: Write>(writer: &mut W, v: f64) -> io::Result<()> {
    writer.write_all(b"\t")?;
    write_double(writer, v)
}

/// Minimal %.6g double formatter (mirrors pin.rs write_double).
fn write_double<W: Write>(writer: &mut W, v: f64) -> io::Result<()> {
    if v.is_nan() || v.is_infinite() || v == 0.0 {
        return writer.write_all(b"0");
    }
    let mut buf = [0u8; 32];
    let abs = v.abs();
    if !(1e-4..1e6).contains(&abs) {
        let len = {
            let mut cursor = &mut buf[..];
            write!(cursor, "{:.5e}", v)?;
            32 - cursor.len()
        };
        write_trim_scientific(writer, &buf[..len])
    } else {
        let digits_before_decimal = abs.log10().floor() as i32 + 1;
        let decimal_places = (6 - digits_before_decimal).max(0) as usize;
        let len = {
            let mut cursor = &mut buf[..];
            write!(cursor, "{:.prec$}", v, prec = decimal_places)?;
            32 - cursor.len()
        };
        write_trim_fixed(writer, &buf[..len])
    }
}

fn write_trim_fixed<W: Write>(writer: &mut W, s: &[u8]) -> io::Result<()> {
    if !s.contains(&b'.') {
        return writer.write_all(s);
    }
    let mut end = s.len();
    while end > 0 && s[end - 1] == b'0' {
        end -= 1;
    }
    if end > 0 && s[end - 1] == b'.' {
        end -= 1;
    }
    writer.write_all(&s[..end])
}

fn write_trim_scientific<W: Write>(writer: &mut W, s: &[u8]) -> io::Result<()> {
    let pos = match s.iter().position(|&b| b == b'e') {
        Some(p) => p,
        None => return writer.write_all(s),
    };
    let mantissa = &s[..pos];
    let exp_part = &s[pos + 1..];
    let mantissa_end = if mantissa.contains(&b'.') {
        let mut end = mantissa.len();
        while end > 0 && mantissa[end - 1] == b'0' {
            end -= 1;
        }
        if end > 0 && mantissa[end - 1] == b'.' {
            end -= 1;
        }
        end
    } else {
        mantissa.len()
    };
    writer.write_all(&mantissa[..mantissa_end])?;
    // Parse and reformat exponent with explicit sign.
    let exp_str = std::str::from_utf8(exp_part).unwrap_or("0");
    let exp_val: i32 = exp_str.parse().unwrap_or(0);
    write!(writer, "e{:+03}", exp_val)
}

#[cfg(test)]
mod tests {
    use super::*;
    use andes_glyco::glycan_db::GlycanComp;
    use andes_glyco::glyco_psm::GlycoPsmKey;
    use andes_glyco::hybrid::Source;
    use search::glyco_search::{FullGlycoPsm, GlycoSpectrumResult};
    use search::psm::{PsmFeatures, PsmMatch};

    fn make_minimal_psm(spectrum_idx: usize, score: f32) -> PsmMatch {
        PsmMatch {
            spectrum_idx,
            candidate_idxs: vec![0],
            charge_used: 2,
            mass_error_ppm: 0.0,
            score,
            rank_score: score,
            edge_score: 0,
            activation_method: None,
            features: PsmFeatures::default(),
            isotope_offset: 0,
            precursor_mz_override: None,
        }
    }

    fn make_key(is_db: bool, glycan_mass: f64) -> GlycoPsmKey {
        let glycan = if is_db {
            Some(GlycanComp {
                hexnac: 4,
                hex: 5,
                fuc: 1,
                neuac: 2,
                neugc: 0,
                mass: glycan_mass,
            })
        } else {
            None
        };
        GlycoPsmKey {
            spectrum_idx: 0,
            glycan,
            glycan_source: if is_db { Source::Db } else { Source::DeNovo },
            oxonium_summed_frac: 0.25,
            n_core_oxonium_ions: 3,
            y_ladder_intensity_score: 1.2,
            y_ladder_decoy_score: 0.3,
            y_hit_frac: 0.0,
            y_hit_frac_decoy: 0.0,
            partial_glycan_by: 0.0,
            y0y1_anchor_score: 0.6,
            sialic_consistency: 0.15,
            core_y_hits: 5,
            glycan_mass,
            backbone_mass: 1500.0,
            is_transferred: false,
            transfer_graph_support: 0,
            transfer_seed_score: 0.0,
            transfer_rt_delta: 0.0,
            transfer_ungated: false,
            cz_hyperscore: 0.0,
            cz_intensity: 0.0,
            cz_explained: 0.0,
            cz_chance_llr: 0.0,
        }
    }

    #[test]
    fn glycan_tag_reports_an_unambiguous_glycosite_and_refuses_an_ambiguous_one() {
        let g = GlycanComp {
            hexnac: 2,
            hex: 5,
            fuc: 0,
            neuac: 0,
            neugc: 0,
            mass: 0.0,
        };
        // exactly one N-X-S/T: the site is a fact about the sequence.
        assert!(glycan_tag_with_site(&g, b"AANLTK").ends_with("@N3]"));
        // the sequon may sit at the very start; positions are 1-based.
        assert!(glycan_tag_with_site(&g, b"NLTAAK").ends_with("@N1]"));
        // two sequons: andes does not localize by default, so it must not pick.
        assert!(glycan_tag_with_site(&g, b"ANLTANCSK").ends_with("@N?]"));
        // no sequon at all is also "not a site andes assigned".
        assert!(glycan_tag_with_site(&g, b"AAAAK").ends_with("@N?]"));
        // the composition itself is unchanged by the site suffix.
        assert!(glycan_tag_with_site(&g, b"AANLTK").starts_with("[HexNAc2Hex5Fuc0NeuAc0NeuGc0"));
    }

    #[test]
    fn glyco_pin_header_contains_glyco_columns() {
        // Build a tiny result and write to a Vec<u8>.
        // We can't call write_glyco_pin_to without real candidates + index,
        // so we test the header directly.
        let mut buf = Vec::new();
        write_glyco_header(&mut buf, 2, 4).unwrap();
        let header = String::from_utf8(buf).unwrap();
        assert!(header.contains("OxoniumScore"), "header missing OxoniumScore");
        assert!(header.contains("NCoreOxoniumIons"), "header missing NCoreOxoniumIons");
        assert!(header.contains("YLadderScore"), "header missing YLadderScore");
        assert!(header.contains("CoreYHits"), "header missing CoreYHits");
        assert!(header.contains("GlycanMass"), "header missing GlycanMass");
        assert!(header.contains("IsGlycanDb"), "header missing IsGlycanDb");
        assert!(header.contains("Peptide"), "header missing Peptide");
        assert!(header.contains("Proteins"), "header missing Proteins");
    }

    #[test]
    fn glyco_pin_header_has_transfer_columns_after_sialic() {
        let mut buf = Vec::new();
        write_glyco_header(&mut buf, 2, 4).unwrap();
        let header = String::from_utf8(buf).unwrap();
        let cols: Vec<&str> = header.trim().split('\t').collect();
        for c in ["IsTransferred","TransferGraphSupport","TransferSeedScore","TransferRTDelta","TransferUngated"] {
            assert!(cols.contains(&c), "header missing {c}");
        }
        let pos = |c: &str| cols.iter().position(|&h| h == c).unwrap();
        assert!(pos("SialicConsistency") < pos("IsTransferred"));
        assert!(pos("TransferUngated") < pos("Peptide"));
    }

    #[test]
    fn glyco_pin_header_deltartrank_is_last_before_peptide() {
        let mut buf = Vec::new();
        write_glyco_header(&mut buf, 2, 4).unwrap();
        let header = String::from_utf8(buf).unwrap();
        let cols: Vec<&str> = header.trim().split('\t').collect();
        let pos = |c: &str| cols.iter().position(|&h| h == c).unwrap();
        // DeltaRT/AbsDeltaRT/DeltaRTNorm keep their in-place positions; DeltaRTRank
        // is APPENDED last after TransferUngated, immediately before Peptide.
        assert!(cols.contains(&"DeltaRTRank"), "header missing DeltaRTRank");
        assert!(pos("DeltaRTNorm") < pos("DeltaRTRank"));
        assert!(pos("TransferUngated") < pos("DeltaRTRank"));
        // IsobaricRTMargin is appended after DeltaRTRank, immediately before Peptide.
        assert!(cols.contains(&"IsobaricRTMargin"), "header missing IsobaricRTMargin");
        assert!(pos("DeltaRTRank") < pos("IsobaricRTMargin"));
        assert_eq!(pos("IsobaricRTMargin") + 1, pos("Peptide"), "IsobaricRTMargin must be last before Peptide");
        // Header column count == a written data row's column count (consistency).
    }

    #[test]
    fn glyco_pin_header_glyco_columns_before_peptide() {
        let mut buf = Vec::new();
        write_glyco_header(&mut buf, 2, 3).unwrap();
        let header = String::from_utf8(buf).unwrap();
        let cols: Vec<&str> = header.trim().split('\t').collect();
        let glycan_pos = cols.iter().position(|&c| c == "GlycanMass").unwrap();
        let peptide_pos = cols.iter().position(|&c| c == "Peptide").unwrap();
        assert!(
            glycan_pos < peptide_pos,
            "GlycanMass col ({}) must come before Peptide col ({})",
            glycan_pos,
            peptide_pos
        );
    }

    #[test]
    fn is_glycan_db_flag_correct() {
        // DB hit → IsGlycanDb = 1
        let db_key = make_key(true, 1234.0);
        assert_eq!(
            if db_key.glycan_source == Source::Db { 1u8 } else { 0u8 },
            1
        );
        // DeNovo hit → IsGlycanDb = 0
        let dn_key = make_key(false, 0.0);
        assert_eq!(
            if dn_key.glycan_source == Source::Db { 1u8 } else { 0u8 },
            0
        );
    }

    // ── BUG 3 regression: CalcMass must include the glycan mass ──────────────

    use model::aa_set::AminoAcidSetBuilder;
    use model::peptide::Peptide;
    use model::protein::{Protein, ProteinDb};
    use model::tolerance::{PrecursorTolerance, Tolerance};
    use search::candidate_gen::Candidate;
    use search::search_params::SearchParams;

    fn make_glyco_search_index(accession: &str) -> SearchIndex {
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: accession.to_string(),
                description: String::new(),
                sequence: b"MKWVTFISLLNSTK".to_vec(),
            }],
        };
        SearchIndex::from_target_db(&target, "XXX_")
    }

    fn make_glyco_candidate() -> Candidate {
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        // "NSTK" contains an N-X-S/T sequon (N-S-T).
        let peptide = Peptide::from_str("K.NSTK.R", &aa_set).expect("valid peptide");
        Candidate {
            peptide,
            protein_index: 0,
            start_offset_in_protein: 0,
            is_decoy: false,
            is_protein_n_term: false,
            is_protein_c_term: false,
        }
    }

    fn make_glyco_spectrum(precursor_mz: f64) -> Spectrum {
        Spectrum {
            title: "test_spec".to_string(),
            precursor_mz,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: Some(1),
            peaks: vec![],
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
        }
    }

    fn make_glyco_params() -> SearchParams {
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        SearchParams {
            aa_set,
            enzyme: model::enzyme::Enzyme::Trypsin,
            extra_enzymes: Vec::new(),
            min_length: 6,
            max_length: 40,
            max_missed_cleavages: 1,
            max_variable_mods_per_peptide: 3,
            precursor_tolerance: PrecursorTolerance::symmetric(Tolerance::Ppm(20.0)),
            charge_range: 2..=3,
            isotope_error_range: -1..=2,
            top_n_psms_per_spectrum: 10,
            num_tolerable_termini: 2,
            min_peaks: 10,
            precursor_cal_mode: search::PrecursorCalMode::Off,
            cal_min_spec_keys: search::precursor_cal::constants::MIN_SPECKEYS_FOR_PREPASS,
            precursor_mass_shift_ppm: 0.0,
            chimeric: false,
            chimeric_isolation_halfwidth_da: 1.5,
            chimeric_max_coisolated: 2,
            chimeric_max_kl: 0.3,
            chimeric_allow_overlap: false,
            score_mode: search::ScoreMode::Rank,
            refine_select_psm_fdr: 0.01,
            candidate_index: search::CandidateIndexMode::Ram,
        }
    }

    /// BUG 3: `CalcMass` in the glyco PIN row must be the INTACT glycopeptide
    /// neutral mass (`peptide.mass() + glycan_mass`), not the bare backbone
    /// alone. Before the fix, `CalcMass = cand.peptide.mass()` only, which
    /// left `dm`/`absdm` off by ~the glycan mass (here 892.317 Da for
    /// HexNAc2Hex3) since the observed `ExpMass` is the intact precursor.
    #[test]
    fn calc_mass_includes_glycan_mass() {
        let search_index = make_glyco_search_index("sp|P00000|TEST");
        let candidates = vec![make_glyco_candidate()];
        let params = make_glyco_params();

        let glycan_mass = 2.0 * andes_glyco::glycan_mass::HEXNAC
            + 3.0 * andes_glyco::glycan_mass::HEX; // ~892.317 Da
        let peptide_neutral_mass = candidates[0].peptide.mass();

        // Precursor m/z consistent with the INTACT glycopeptide at z=2 so
        // `dm` should be near 0 once CalcMass correctly includes the glycan.
        let intact_neutral = peptide_neutral_mass + glycan_mass;
        let z = 2.0_f64;
        let precursor_mz = (intact_neutral + z * PROTON) / z;
        let spectrum = make_glyco_spectrum(precursor_mz);
        let spectra = vec![spectrum];

        let psm = make_minimal_psm(0, 10.0);
        let glycan_key = GlycoPsmKey {
            spectrum_idx: 0,
            glycan: Some(GlycanComp {
                hexnac: 2,
                hex: 3,
                fuc: 0,
                neuac: 0,
                neugc: 0,
                mass: glycan_mass,
            }),
            glycan_source: Source::Db,
            oxonium_summed_frac: 0.2,
            n_core_oxonium_ions: 2,
            y_ladder_intensity_score: 0.0,
            y_hit_frac: 0.0,
            y_hit_frac_decoy: 0.0,
            y_ladder_decoy_score: 0.0,
            partial_glycan_by: 0.0,
            y0y1_anchor_score: 0.0,
            sialic_consistency: 0.0,
            core_y_hits: 3,
            glycan_mass,
            backbone_mass: peptide_neutral_mass,
            is_transferred: false,
            transfer_graph_support: 0,
            transfer_seed_score: 0.0,
            transfer_rt_delta: 0.0,
            transfer_ungated: false,
            cz_hyperscore: 0.0,
            cz_intensity: 0.0,
            cz_explained: 0.0,
            cz_chance_llr: 0.0,
        };
        let hit = FullGlycoPsm { glycan_key, psm };
        let results = vec![GlycoSpectrumResult { spectrum_idx: 0, hits: vec![hit] }];

        let mut buf = Vec::new();
        write_glyco_pin_to(&mut buf, &spectra, &results, &candidates, &params, &search_index, false, false)
            .expect("write_glyco_pin_to must succeed");
        let text = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert!(lines.len() >= 2, "expected header + at least one PSM row");

        let header_cols: Vec<&str> = lines[0].split('\t').collect();
        let calc_mass_col = header_cols.iter().position(|&c| c == "CalcMass").unwrap();
        let dm_col = header_cols.iter().position(|&c| c == "dm").unwrap();

        let row_cols: Vec<&str> = lines[1].split('\t').collect();
        let calc_mass: f64 = row_cols[calc_mass_col].parse().unwrap();
        let dm: f64 = row_cols[dm_col].parse().unwrap();

        // `write_double` formats with ~6 significant digits (see `write_double`
        // in this file), so a masses-in-the-thousands value round-trips to
        // within ~0.01 Da, not exact float precision. The bug this test
        // guards against (omitting `glycan_mass`, ~892 Da here) is 4+ orders
        // of magnitude larger than this formatting tolerance.
        assert!(
            (calc_mass - intact_neutral).abs() < 0.05,
            "CalcMass must equal peptide.mass() + glycan_mass ({:.4}), got {:.4} \
             (bare peptide mass alone would be {:.4})",
            intact_neutral,
            calc_mass,
            peptide_neutral_mass
        );
        assert!(
            dm.abs() < 0.05,
            "dm should be ~0 for a precursor matched to the intact glycopeptide mass \
             once CalcMass includes the glycan; got dm={:.6} (bug: CalcMass omitting \
             glycan_mass would produce dm ~= -{:.4})",
            dm,
            glycan_mass
        );
    }

    /// L1 (the core harness bug): the FDR PIN must emit ONE PSM per scan — the
    /// top-1 TDC winner by the shared `gp` fused score. Emitting all ~4.5
    /// (peptide×glycan) alternatives fabricates the target/decoy balance and
    /// produced the fake ~30% recovery. `select_emitted_hits(_, true)` must
    /// return only the best. The gp score `rank + K·ladder + J·core_y` (K=50)
    /// is ladder-dominated, so the highest-`y_ladder` hit wins with `rank_score`
    /// as the tiebreak.
    #[test]
    fn top1_collapse_keeps_only_the_best_hit_per_scan() {
        fn make_hit(rank: f32, ladder: f32) -> FullGlycoPsm {
            let mut key = make_key(true, 1000.0);
            key.y_ladder_intensity_score = ladder;
            FullGlycoPsm { glycan_key: key, psm: make_minimal_psm(0, rank) }
        }
        let hits = vec![
            make_hit(5.0, 1.0),
            make_hit(10.0, 0.5),
            make_hit(8.0, 2.0), // highest y_ladder → the winner (ladder-primary default)
        ];
        assert_eq!(select_emitted_hits(&hits, true, false), vec![2], "collapse keeps top-1 by y_ladder (default)");
        assert_eq!(select_emitted_hits(&hits, false, false).len(), 3, "diagnostic mode keeps all");

        // Full tie on y_ladder → break by rank_score, then lowest index.
        let tied = vec![make_hit(1.0, 7.0), make_hit(3.0, 7.0), make_hit(3.0, 7.0)];
        assert_eq!(select_emitted_hits(&tied, true, false), vec![1], "tie broken by rank_score then lowest index");
    }

    /// The `gp` fused `rank + K·ladder` collapse (the shipped default selector)
    /// rescues a rank-strong hit that a bare ladder-primary collapse would reject for
    /// a spurious tiny ladder edge (the P0 failure mode: a wrong mass-split wins on a
    /// coincidental summed Y-ladder intensity despite truth's stronger b/y). This
    /// asserts the exact ordering expression `select_emitted_hits` uses.
    #[test]
    fn gp_fused_collapse_rescues_rank_strong_hit_over_spurious_ladder() {
        use andes_glyco::glyco_psm::{glyco_gp_fused_score, GLYCO_GP_J_DEFAULT, GLYCO_GP_K_DEFAULT};
        fn make_hit(rank: f32, ladder: f32) -> FullGlycoPsm {
            let mut key = make_key(true, 1000.0);
            key.y_ladder_intensity_score = ladder;
            FullGlycoPsm { glycan_key: key, psm: make_minimal_psm(0, rank) }
        }
        // idx 0 = truth (strong b/y rank, modest ladder); idx 1 = wrong split
        // (spurious tiny ladder edge that a bare ladder-primary collapse would pick).
        let hits = vec![make_hit(15.0, 0.05), make_hit(2.0, 0.06)];
        // The shipped gp default rescues the rank-strong truth (idx 0):
        // 15 + 50·0.05 = 17.5 > 2 + 50·0.06 = 5 (a bare ladder collapse would pick idx 1).
        assert_eq!(
            select_emitted_hits(&hits, true, false),
            vec![0],
            "gp default rescues the rank-strong truth over the spurious-ladder split"
        );
        // The gp fused ordering (exactly what select_emitted_hits computes) rescues
        // the rank-strong truth (idx 0): 15 + 50·0.05 = 17.5 > 2 + 50·0.06 = 5.
        let (k, j) = (GLYCO_GP_K_DEFAULT, GLYCO_GP_J_DEFAULT);
        let gp_score = |h: &FullGlycoPsm| {
            glyco_gp_fused_score(
                h.psm.rank_score,
                h.glycan_key.y_ladder_intensity_score,
                h.glycan_key.core_y_hits as f32,
                0.0, // hyperscore: PIN-side single-hit collapse, never decisive (see select_emitted_hits)
                k,
                j,
                0.0,
            )
        };
        let gp_winner = (0..hits.len())
            .max_by(|&a, &b| gp_score(&hits[a]).total_cmp(&gp_score(&hits[b])).then(b.cmp(&a)))
            .unwrap();
        assert_eq!(gp_winner, 0, "gp fusion rescues the rank-strong truth over the spurious-ladder split");
    }

    /// GI-1: a de-novo hit (glycan = None, a bare mass residual) is NOT a glyco
    /// identification and must be dropped by the enumerated-only filter — even
    /// when it has the highest rank_score. The best ENUMERATED hit wins the scan.
    #[test]
    fn enumerated_only_drops_de_novo_hits() {
        fn hit(is_db: bool, rank: f32) -> FullGlycoPsm {
            FullGlycoPsm { glycan_key: make_key(is_db, 1500.0), psm: make_minimal_psm(0, rank) }
        }
        // de-novo hit has the HIGHER rank_score → it is the TDC winner. Under
        // enumerated-only the scan is DROPPED (not a lower enumerated hit promoted).
        let hits = vec![hit(false, 20.0), hit(true, 8.0), hit(true, 5.0)];
        assert!(
            select_emitted_hits(&hits, true, true).is_empty(),
            "de-novo TDC winner → scan dropped, NOT a losing enumerated hit promoted"
        );
        // enumerated_only=false: de-novo (rank 20, idx 0) wins the raw collapse.
        assert_eq!(select_emitted_hits(&hits, true, false), vec![0], "without the filter the de-novo residual wins");
        // When the ENUMERATED hit is the winner, it is emitted.
        let db_wins = vec![hit(true, 20.0), hit(false, 8.0), hit(true, 5.0)];
        assert_eq!(select_emitted_hits(&db_wins, true, true), vec![0], "enumerated winner is emitted");
        // A scan with ONLY de-novo hits yields nothing under enumerated-only.
        let denovo_only = vec![hit(false, 9.0), hit(false, 4.0)];
        assert!(select_emitted_hits(&denovo_only, true, true).is_empty(), "no enumerated glycan → no ID");
    }

    /// G3: with glycan-decoy emission ON, a TARGET glyco-PSM that has a resolved
    /// glycan must emit a PAIRED glycan-decoy row: Label -1, YLadderScore taken
    /// from `y_ladder_decoy_score` (below the target's), so Percolator has a
    /// glycan-axis decoy whose glyco feature actually differs from the target.
    #[test]
    fn emit_glycan_decoy_writes_paired_decoy_row() {
        let search_index = make_glyco_search_index("sp|P00000|TEST");
        let candidates = vec![make_glyco_candidate()]; // is_decoy = false (target)
        let params = make_glyco_params();

        let glycan_mass =
            2.0 * andes_glyco::glycan_mass::HEXNAC + 3.0 * andes_glyco::glycan_mass::HEX;
        let z = 2.0_f64;
        let intact_neutral = candidates[0].peptide.mass() + glycan_mass;
        let precursor_mz = (intact_neutral + z * PROTON) / z;
        let spectra = vec![make_glyco_spectrum(precursor_mz)];

        let psm = make_minimal_psm(0, 10.0);
        let glycan_key = GlycoPsmKey {
            spectrum_idx: 0,
            glycan: Some(GlycanComp {
                hexnac: 2, hex: 3, fuc: 0, neuac: 0, neugc: 0, mass: glycan_mass,
            }),
            glycan_source: Source::Db,
            oxonium_summed_frac: 0.2,
            n_core_oxonium_ions: 2,
            y_ladder_intensity_score: 1.2,
            y_ladder_decoy_score: 0.3,
            y_hit_frac: 0.0,
            y_hit_frac_decoy: 0.0,
            partial_glycan_by: 0.0,
            y0y1_anchor_score: 0.6,
            sialic_consistency: 0.15,
            core_y_hits: 4,
            glycan_mass,
            backbone_mass: candidates[0].peptide.mass(),
            is_transferred: false,
            transfer_graph_support: 0,
            transfer_seed_score: 0.0,
            transfer_rt_delta: 0.0,
            transfer_ungated: false,
            cz_hyperscore: 0.0,
            cz_intensity: 0.0,
            cz_explained: 0.0,
            cz_chance_llr: 0.0,
        };
        let hit = FullGlycoPsm { glycan_key, psm };
        let results = vec![GlycoSpectrumResult { spectrum_idx: 0, hits: vec![hit] }];

        let mut buf = Vec::new();
        // emit_glycan_decoy = true
        write_glyco_pin_to(&mut buf, &spectra, &results, &candidates, &params, &search_index, true, false)
            .expect("write must succeed");
        let text = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        let header: Vec<&str> = lines[0].split('\t').collect();
        let label_col = header.iter().position(|&c| c == "Label").unwrap();
        let yl_col = header.iter().position(|&c| c == "YLadderScore").unwrap();

        let data: Vec<Vec<&str>> = lines[1..].iter().map(|l| l.split('\t').collect()).collect();
        assert_eq!(data.len(), 2, "target + paired glycan-decoy expected, got {}", data.len());

        let targets: Vec<&Vec<&str>> = data.iter().filter(|r| r[label_col] == "1").collect();
        let decoys: Vec<&Vec<&str>> = data.iter().filter(|r| r[label_col] == "-1").collect();
        assert_eq!(targets.len(), 1, "one target row");
        assert_eq!(decoys.len(), 1, "one glycan-decoy row");

        let yt: f64 = targets[0][yl_col].parse().unwrap();
        let yd: f64 = decoys[0][yl_col].parse().unwrap();
        assert!((yt - 1.2).abs() < 1e-6, "target YLadderScore 1.2, got {yt}");
        assert!((yd - 0.3).abs() < 1e-6, "decoy YLadderScore must be the decoy score 0.3, got {yd}");
    }

    /// Every emitted glyco PIN row must have exactly as many tab-separated
    /// fields as the header (header↔row column-count consistency), and the
    /// DeltaRTRank column must be present and parseable. This is the guard that
    /// appending DeltaRTRank did not desynchronise the header and the rows.
    #[test]
    fn header_and_row_column_counts_match_with_delta_rt_rank() {
        let search_index = make_glyco_search_index("sp|P00000|TEST");
        let candidates = vec![make_glyco_candidate()];
        let params = make_glyco_params();

        let glycan_mass =
            2.0 * andes_glyco::glycan_mass::HEXNAC + 3.0 * andes_glyco::glycan_mass::HEX;
        let z = 2.0_f64;
        let intact_neutral = candidates[0].peptide.mass() + glycan_mass;
        let precursor_mz = (intact_neutral + z * PROTON) / z;
        let spectra = vec![make_glyco_spectrum(precursor_mz)];

        let mut key = make_key(true, glycan_mass);
        key.glycan = Some(GlycanComp {
            hexnac: 2, hex: 3, fuc: 0, neuac: 0, neugc: 0, mass: glycan_mass,
        });
        let hit = FullGlycoPsm { glycan_key: key, psm: make_minimal_psm(0, 10.0) };
        let results = vec![GlycoSpectrumResult { spectrum_idx: 0, hits: vec![hit] }];

        let mut buf = Vec::new();
        write_glyco_pin_to(&mut buf, &spectra, &results, &candidates, &params, &search_index, false, false)
            .expect("write must succeed");
        let text = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        let header: Vec<&str> = lines[0].split('\t').collect();
        let ncols = header.len();
        for (i, row) in lines[1..].iter().enumerate() {
            let cols: Vec<&str> = row.split('\t').collect();
            assert_eq!(cols.len(), ncols, "row {i} has {} cols, header has {ncols}", cols.len());
        }
        // DeltaRTRank present; single-hit scan ⇒ neutral 0.0 (parseable).
        let rank_col = header.iter().position(|&c| c == "DeltaRTRank").unwrap();
        let row: Vec<&str> = lines[1].split('\t').collect();
        let rank: f64 = row[rank_col].parse().unwrap();
        assert_eq!(rank, 0.0, "single-hit scan ⇒ DeltaRTRank neutral 0.0");
    }
}
