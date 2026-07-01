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
) -> io::Result<()> {
    let psm = &hit.psm;
    let key = &hit.glycan_key;

    let ctx = RowContext::new(spec, cand, search_index);
    let scan = ctx.scan;
    let spec_id = format!("{}_glyco_{}_{}", ctx.spec_id, scan, row_idx + 1);

    let label: i32 = if cand.is_decoy { -1 } else { 1 };

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

    // Glyco-specific columns (6).
    let is_glycan_db: u8 = if key.glycan_source == Source::Db { 1 } else { 0 };
    write_double_tab(writer, key.oxonium_summed_frac as f64)?;
    write!(writer, "\t{}", key.n_core_oxonium_ions)?;
    write_double_tab(writer, key.y_ladder_intensity_score as f64)?;
    write!(writer, "\t{}", key.core_y_hits)?;
    write_double_tab(writer, key.glycan_mass)?;
    write!(writer, "\t{}", is_glycan_db)?;

    // Peptide column: backbone sequence + optional glycan tag.
    let glycan_tag = match &key.glycan {
        Some(g) => format!(
            "[HexNAc{}Hex{}Fuc{}NeuAc{}NeuGc{}]",
            g.hexnac, g.hex, g.fuc, g.neuac, g.neugc
        ),
        None => String::new(),
    };
    write!(writer, "\t{}{}", cand.peptide, glycan_tag)?;

    // Proteins column(s)
    for &cidx in &psm.candidate_idxs {
        let cand_for_acc = &candidates[cidx as usize];
        let accession = crate::row_context::resolve_accession(cand_for_acc, search_index);
        write!(writer, "\t{}", accession)?;
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
) -> io::Result<()> {
    let file = std::fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    write_glyco_pin_to(&mut writer, spectra, results, candidates, params, search_index)
}

/// Write glyco PIN to an arbitrary writer (useful for testing).
pub fn write_glyco_pin_to<W: Write>(
    writer: &mut W,
    spectra: &[Spectrum],
    results: &[GlycoSpectrumResult],
    candidates: &[Candidate],
    params: &SearchParams,
    search_index: &SearchIndex,
) -> io::Result<()> {
    let min_charge = *params.charge_range.start();
    let max_charge = *params.charge_range.end();

    write_glyco_header(writer, min_charge, max_charge)?;

    let mut row_count = 0usize;
    for result in results {
        let spec_idx = result.spectrum_idx;
        if spec_idx >= spectra.len() {
            continue;
        }
        let spec = &spectra[spec_idx];
        for (hit_idx, hit) in result.hits.iter().enumerate() {
            let cand_idx = hit.psm.primary_candidate_idx() as usize;
            if cand_idx >= candidates.len() {
                continue;
            }
            let cand = &candidates[cand_idx];
            write_glyco_psm_row(
                writer,
                spec,
                hit,
                cand,
                hit_idx,
                min_charge,
                max_charge,
                candidates,
                search_index,
                params,
            )?;
            row_count += 1;
        }
    }

    writer.flush()?;
    // Return row count via a side-channel eprintln so callers can report it.
    eprintln!("[glyco-pin] wrote {} PSM rows", row_count);
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
            core_y_hits: 3,
            glycan_mass,
            backbone_mass: peptide_neutral_mass,
        };
        let hit = FullGlycoPsm { glycan_key, psm };
        let results = vec![GlycoSpectrumResult { spectrum_idx: 0, hits: vec![hit] }];

        let mut buf = Vec::new();
        write_glyco_pin_to(&mut buf, &spectra, &results, &candidates, &params, &search_index)
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
}
