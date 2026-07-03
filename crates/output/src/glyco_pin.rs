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
//   `[HexNAc{n}Hex{m}Fuc{f}NeuAc{a}NeuGc{g}]`

use std::io::{self, BufWriter, Write};
use std::path::Path;

use model::mass::{ISOTOPE, PROTON};
use search::candidate_gen::Candidate;
use search::search_index::SearchIndex;
use search::search_params::SearchParams;
use model::spectrum::Spectrum;

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
        // Glyco-specific columns (appended after RankScoreFloat, before Peptide/Proteins)
        "OxoniumScore".to_string(),
        "NCoreOxoniumIons".to_string(),
        "YLadderScore".to_string(),
        "CoreYHits".to_string(),
        "GlycanMass".to_string(),
        "IsGlycanDb".to_string(),
        "Y0Y1Anchor".to_string(), // G2 peptide-mass anchor (additive, peptide-discriminating)
        "SialicConsistency".to_string(), // GI-2 composition-conditioned sialic-oxonium (additive)
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

    // Glyco-specific columns (7, incl. Y0Y1Anchor).
    let is_glycan_db: u8 = if key.glycan_source == Source::Db { 1 } else { 0 };
    write_double_tab(writer, key.oxonium_summed_frac as f64)?;
    write!(writer, "\t{}", key.n_core_oxonium_ions)?;
    // YLadderScore: the decoy score for a glycan-decoy row, else the target score.
    let y_ladder = glycan_decoy_override.unwrap_or(key.y_ladder_intensity_score);
    write_double_tab(writer, y_ladder as f64)?;
    write!(writer, "\t{}", key.core_y_hits)?;
    write_double_tab(writer, key.glycan_mass)?;
    write!(writer, "\t{}", is_glycan_db)?;
    // G2 Y0/Y1 anchor (additive, peptide-mass-conditioned). Same value on a
    // glycan-decoy row (it is a peptide-axis feature, independent of the glycan).
    write_double_tab(writer, key.y0y1_anchor_score as f64)?;
    // GI-2 sialic consistency (composition-conditioned; same value on a glycan-decoy
    // row since the shifted-ladder decoy keeps the composition — differs once
    // isobaric-composition decoys land).
    write_double_tab(writer, key.sialic_consistency as f64)?;

    // Peptide column: backbone sequence + optional glycan tag.
    let glycan_tag = match &key.glycan {
        Some(g) => format!(
            "[HexNAc{}Hex{}Fuc{}NeuAc{}NeuGc{}]",
            g.hexnac, g.hex, g.fuc, g.neuac, g.neugc
        ),
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

// ── public API ────────────────────────────────────────────────────────────────

/// Write a glyco PIN file from the output of `glyco_search_run`.
///
/// Standard PIN columns are written using `psm_feature_values`; six
/// glyco-specific columns are appended after `RankScoreFloat`:
///   OxoniumScore, NCoreOxoniumIons, YLadderScore, CoreYHits,
///   GlycanMass, IsGlycanDb
///
/// `Peptide` column appends `[HexNAc{n}Hex{m}...]` for DB glycan hits.
pub fn write_glyco_pin(
    path: &Path,
    spectra: &[Spectrum],
    results: &[GlycoSpectrumResult],
    candidates: &[Candidate],
    params: &SearchParams,
    search_index: &SearchIndex,
    emit_glycan_decoy: bool,
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
    let eligible: Vec<usize> = (0..hits.len())
        .filter(|&i| !enumerated_only || hits[i].glycan_key.glycan.is_some())
        .collect();
    if !collapse || eligible.len() <= 1 {
        return eligible;
    }
    eligible
        .iter()
        .copied()
        .max_by(|&a, &b| {
            hits[a]
                .psm
                .rank_score
                .total_cmp(&hits[b].psm.rank_score)
                .then(
                    hits[a]
                        .glycan_key
                        .y_ladder_intensity_score
                        .total_cmp(&hits[b].glycan_key.y_ladder_intensity_score),
                )
                // lowest index wins a full tie (deterministic)
                .then(b.cmp(&a))
        })
        .into_iter()
        .collect()
}

/// Write glyco PIN to an arbitrary writer (useful for testing).
pub fn write_glyco_pin_to<W: Write>(
    writer: &mut W,
    spectra: &[Spectrum],
    results: &[GlycoSpectrumResult],
    candidates: &[Candidate],
    params: &SearchParams,
    search_index: &SearchIndex,
    emit_glycan_decoy: bool,
) -> io::Result<()> {
    let min_charge = *params.charge_range.start();
    let max_charge = *params.charge_range.end();

    write_glyco_header(writer, min_charge, max_charge)?;

    // Top-1-per-scan collapse is the honest default (see `select_emitted_hits`).
    // ANDES_GLYCO_ALL_HITS=1 restores the full multi-row dump for DIAGNOSTICS ONLY
    // (its PIN must never be fed to Percolator as an FDR input).
    let collapse = std::env::var("ANDES_GLYCO_ALL_HITS")
        .map(|v| v != "1")
        .unwrap_or(true);
    // GI-1: emit only enumerated-composition glyco IDs by default; ANDES_GLYCO_DENOVO=1
    // includes the de-novo mass-residual hits (diagnostic / open-search only).
    let enumerated_only = std::env::var("ANDES_GLYCO_DENOVO")
        .map(|v| v != "1")
        .unwrap_or(true);

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
                search_index, params, None,
            )?;
            row_count += 1;

            // G3: emit a paired glycan-axis decoy for TARGET-peptide PSMs that
            // have a resolved glycan composition. Skip de-novo (no composition →
            // decoy score 0) and peptide-decoy rows (they are already decoys; a
            // glycan-decoy of a decoy peptide adds no glycan-axis signal).
            if emit_glycan_decoy && !cand.is_decoy && hit.glycan_key.glycan.is_some() {
                write_glyco_psm_row(
                    writer, spec, hit, cand, hit_idx, min_charge, max_charge, candidates,
                    search_index, params, Some(hit.glycan_key.y_ladder_decoy_score),
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
            y0y1_anchor_score: 0.6,
            sialic_consistency: 0.15,
            core_y_hits: 5,
            glycan_mass,
            backbone_mass: 1500.0,
        }
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
            y_ladder_decoy_score: 0.0,
            y0y1_anchor_score: 0.0,
            sialic_consistency: 0.0,
            core_y_hits: 3,
            glycan_mass,
            backbone_mass: peptide_neutral_mass,
        };
        let hit = FullGlycoPsm { glycan_key, psm };
        let results = vec![GlycoSpectrumResult { spectrum_idx: 0, hits: vec![hit] }];

        let mut buf = Vec::new();
        write_glyco_pin_to(&mut buf, &spectra, &results, &candidates, &params, &search_index, false)
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
    /// top-1 TDC winner by rank_score. Emitting all ~4.5 (peptide×glycan)
    /// alternatives fabricates the target/decoy balance and produced the fake
    /// ~30% recovery. `select_emitted_hits(_, true)` must return only the best.
    #[test]
    fn top1_collapse_keeps_only_the_best_hit_per_scan() {
        fn make_hit(rank: f32, ladder: f32) -> FullGlycoPsm {
            let mut key = make_key(true, 1000.0);
            key.y_ladder_intensity_score = ladder;
            FullGlycoPsm { glycan_key: key, psm: make_minimal_psm(0, rank) }
        }
        let hits = vec![
            make_hit(5.0, 1.0),
            make_hit(10.0, 0.5), // highest rank_score → the winner
            make_hit(8.0, 2.0),
        ];
        assert_eq!(select_emitted_hits(&hits, true, false), vec![1], "collapse keeps top-1 by rank_score");
        assert_eq!(select_emitted_hits(&hits, false, false).len(), 3, "diagnostic mode keeps all");

        // Full tie on rank_score → break by YLadderScore, then lowest index.
        let tied = vec![make_hit(7.0, 1.0), make_hit(7.0, 3.0), make_hit(7.0, 3.0)];
        assert_eq!(select_emitted_hits(&tied, true, false), vec![1], "tie broken by YLadder then lowest index");
    }

    /// GI-1: a de-novo hit (glycan = None, a bare mass residual) is NOT a glyco
    /// identification and must be dropped by the enumerated-only filter — even
    /// when it has the highest rank_score. The best ENUMERATED hit wins the scan.
    #[test]
    fn enumerated_only_drops_de_novo_hits() {
        fn hit(is_db: bool, rank: f32) -> FullGlycoPsm {
            FullGlycoPsm { glycan_key: make_key(is_db, 1500.0), psm: make_minimal_psm(0, rank) }
        }
        // de-novo hit has the HIGHER rank_score, but no enumerated glycan.
        let hits = vec![hit(false, 20.0), hit(true, 8.0), hit(true, 5.0)];
        // enumerated_only=true: de-novo dropped, best enumerated (rank 8.0, idx 1) wins.
        assert_eq!(select_emitted_hits(&hits, true, true), vec![1], "enumerated-only keeps best DB hit, drops de-novo");
        // enumerated_only=false: de-novo (rank 20, idx 0) wins the raw collapse.
        assert_eq!(select_emitted_hits(&hits, true, false), vec![0], "without the filter the de-novo residual wins");
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
            y0y1_anchor_score: 0.6,
            sialic_consistency: 0.15,
            core_y_hits: 4,
            glycan_mass,
            backbone_mass: candidates[0].peptide.mass(),
        };
        let hit = FullGlycoPsm { glycan_key, psm };
        let results = vec![GlycoSpectrumResult { spectrum_idx: 0, hits: vec![hit] }];

        let mut buf = Vec::new();
        // emit_glycan_decoy = true
        write_glyco_pin_to(&mut buf, &spectra, &results, &candidates, &params, &search_index, true)
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
}
