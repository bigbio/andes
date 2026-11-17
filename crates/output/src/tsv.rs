use std::io::{self, BufWriter, Write};

use crate::row_context::{emitted_rows, RowContext};
use model::mass::{ISOTOPE, PROTON};
use model::peptide::Peptide;
use model::protein::{EntryKind, EntryOrigin};
use model::spectrum::Spectrum;
use model::tolerance::Tolerance;
use search::candidate_gen::Candidate;
use search::psm::{PsmMatch, TopNQueue};
use search::search_index::SearchIndex;
use search::search_params::SearchParams;

pub(crate) fn modifications_field(peptide: &Peptide) -> String {
    let mut parts: Vec<String> = peptide
        .residues
        .iter()
        .enumerate()
        .flat_map(|(i, aa)| {
            aa.mod_
                .iter()
                .flat_map(|m| m.components())
                .filter_map(|m| m.accession.as_deref())
                .map(move |acc| format!("{}:{}", i + 1, acc))
                .collect::<Vec<_>>()
        })
        .collect();
    parts.sort_by_key(|s| {
        s.split(':')
            .next()
            .and_then(|n| n.parse::<usize>().ok())
            .unwrap_or(0)
    });
    parts.join(";")
}

#[allow(
    clippy::too_many_arguments,
    reason = "Writer API mirrors PIN writer; grouping into a struct would diverge from the parallel write_pin API"
)]
pub fn write_tsv(
    output_path: &std::path::Path,
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &[Candidate],
    params: &SearchParams,
    search_index: &SearchIndex,
    spec_file_name: &str,
    is_mgf: bool,
) -> io::Result<()> {
    crate::atomic::write_atomic(output_path, |path| {
        let mut writer = BufWriter::new(std::fs::File::create(path)?);
        write_tsv_to(
            &mut writer,
            spectra,
            queues,
            candidates,
            params,
            search_index,
            spec_file_name,
            is_mgf,
        )?;
        writer.flush()
    })
}

#[allow(
    clippy::too_many_arguments,
    reason = "Writer API mirrors PIN writer; grouping into a struct would diverge from the parallel write_pin API"
)]
pub fn write_tsv_to<W: Write>(
    writer: &mut W,
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &[Candidate],
    params: &SearchParams,
    search_index: &SearchIndex,
    spec_file_name: &str,
    is_mgf: bool,
) -> io::Result<()> {
    write_header(writer, params, is_mgf)?;
    for (spec_idx, queue) in queues.iter().enumerate() {
        if queue.is_empty() {
            continue;
        }
        let spec = &spectra[spec_idx];
        write_spectrum_rows(
            writer,
            spec,
            queue,
            candidates,
            params,
            spec_file_name,
            is_mgf,
            search_index,
        )?;
    }
    Ok(())
}

fn write_header<W: Write>(writer: &mut W, params: &SearchParams, is_mgf: bool) -> io::Result<()> {
    let ppm_mode = matches!(params.precursor_tolerance.left, Tolerance::Ppm(_));
    let prec_err_col = if ppm_mode {
        "PrecursorError(ppm)"
    } else {
        "PrecursorError(Da)"
    };

    let mut cols: Vec<&str> = vec!["#SpecFile", "SpecID", "ScanNum"];
    if is_mgf {
        cols.push("Title");
    }
    cols.extend_from_slice(&[
        "FragMethod",
        "Precursor",
        "IsotopeError",
        prec_err_col,
        "Charge",
        "Peptide",
        "Protein",
    ]);
    cols.push("RankScore");
    cols.push("Modifications");
    cols.extend_from_slice(&["ProteinFamily", "EntryKind", "FeatureIds"]);
    if params.open_search.is_some() {
        cols.extend_from_slice(&["DeltaMass", "OpenSite", "OpenLocScore", "OpenLocMargin"]);
    }

    writeln!(writer, "{}", cols.join("\t"))
}

struct RowCtx<'a> {
    spec_file_name: &'a str,
    is_mgf: bool,
    ppm_mode: bool,
    open_search: bool,
}

#[allow(
    clippy::too_many_arguments,
    reason = "Writer API mirrors PIN writer; grouping into a struct would diverge from the parallel write_pin API"
)]
fn write_spectrum_rows<W: Write>(
    writer: &mut W,
    spec: &Spectrum,
    queue: &TopNQueue,
    candidates: &[Candidate],
    params: &SearchParams,
    spec_file_name: &str,
    is_mgf: bool,
    search_index: &SearchIndex,
) -> io::Result<()> {
    let psms = queue.clone().into_rank_sorted_vec();

    let row_ctx = RowCtx {
        spec_file_name,
        is_mgf,
        ppm_mode: matches!(params.precursor_tolerance.left, Tolerance::Ppm(_)),
        open_search: params.open_search.is_some(),
    };

    for (_rank, psm) in emitted_rows(spec, &psms, candidates) {
        let cand = &candidates[psm.primary_candidate_idx() as usize];
        let ctx = RowContext::new(spec, cand, search_index);
        let origin = search_index
            .protein_at(cand.protein_index)
            .and_then(|p| p.origin());
        write_psm_row(writer, spec, &psm, cand, &ctx, origin, &row_ctx)?;
    }
    Ok(())
}

fn write_psm_row<W: Write>(
    writer: &mut W,
    spec: &Spectrum,
    psm: &PsmMatch,
    cand: &Candidate,
    ctx: &RowContext,
    origin: Option<&EntryOrigin>,
    row_ctx: &RowCtx<'_>,
) -> io::Result<()> {
    let is_mgf = row_ctx.is_mgf;
    let ppm_mode = row_ctx.ppm_mode;
    let spec_file_name = row_ctx.spec_file_name;

    let spec_id = &ctx.spec_id;

    let scan_num = ctx.scan;

    let frag_method = psm
        .activation_method
        .map(|m| m.name().to_string())
        .unwrap_or_else(|| "UNKNOWN".to_string());

    let precursor_mz = psm.precursor_mz_override.unwrap_or(spec.precursor_mz);
    let precursor = format!("{:.4}", precursor_mz);

    let isotope_error: i32 = psm.isotope_offset as i32;

    let precursor_error = if ppm_mode {
        format!("{:.4}", psm.mass_error_ppm)
    } else {
        let da = psm.mass_error_ppm * 1e-6 * precursor_mz;
        format!("{:.4}", da)
    };

    let charge = psm.charge_used;

    let peptide = &cand.peptide;
    let protein = &ctx.accession;

    let raw_score = psm.score.round() as i32;

    let lead = if is_mgf {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            spec_file_name,
            spec_id,
            scan_num,
            spec.title,
            frag_method,
            precursor,
            isotope_error,
            precursor_error,
            charge,
            peptide,
            protein,
        )
    } else {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            spec_file_name,
            spec_id,
            scan_num,
            frag_method,
            precursor,
            isotope_error,
            precursor_error,
            charge,
            peptide,
            protein,
        )
    };
    let mods = modifications_field(&cand.peptide);
    let (family, kind, feature) = match origin {
        Some(o) => (&*o.family, o.kind.name(), &*o.feature),
        None => (protein.as_str(), EntryKind::Reference.name(), ""),
    };
    let open_search_tail = if row_ctx.open_search {
        let exp_mass = precursor_mz * f64::from(charge) - f64::from(charge) * PROTON
            - f64::from(psm.isotope_offset) * ISOTOPE;
        let f = &psm.features;
        let delta = if f.open_delta != 0.0 {
            f64::from(f.open_delta)
        } else {
            exp_mass - cand.peptide.mass()
        };
        format!(
            "\t{:.4}\t{}\t{:.4}\t{:.4}",
            delta, f.open_site, f.open_loc_score, f.open_loc_margin
        )
    } else {
        String::new()
    };
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}{}",
        lead, raw_score, mods, family, kind, feature, open_search_tail
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::amino_acid::AminoAcid;
    use model::modification::{ModLocation, Modification, ResidueSpec};
    use model::peptide::Peptide;
    use model::protein::{Protein, ProteinDb};
    use model::tolerance::PrecursorTolerance;
    use search::candidate_gen::Candidate;
    use search::search_index::SearchIndex;

    fn make_search_index(accession: &str) -> SearchIndex {
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: accession.to_string(),
                description: String::new(),
                sequence: b"MKWVTFISLL".to_vec(),
                annotation: None,
            }],
        };
        SearchIndex::from_target_db(&target, "XXX_")
    }

    fn make_empty_search_index() -> SearchIndex {
        let target = ProteinDb { proteins: vec![] };
        SearchIndex::from_target_db(&target, "XXX_")
    }

    fn make_spectrum(title: &str, scan: i32, precursor_mz: f64) -> Spectrum {
        Spectrum {
            title: title.to_string(),
            precursor_mz,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: Some(scan),
            peaks: vec![],
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
        }
    }

    fn make_candidate(protein_index: usize, is_decoy: bool) -> Candidate {
        let aa = AminoAcid::standard(b'A').unwrap();
        let peptide = Peptide::new(vec![aa], b'K', b'S');
        Candidate {
            peptide,
            protein_index,
            start_offset_in_protein: 0,
            is_decoy,
            is_protein_n_term: false,
            is_protein_c_term: false,
        }
    }

    fn make_psm(spectrum_idx: usize, score: f32, rank_score: f32) -> PsmMatch {
        PsmMatch {
            spectrum_idx,
            candidate_idxs: vec![0],
            charge_used: 2,
            mass_error_ppm: 1.5,
            score,
            rank_score,
            edge_score: 0,
            activation_method: Some(model::activation::ActivationMethod::HCD),
            features: search::psm::PsmFeatures::default(),
            isotope_offset: 0,
            precursor_mz_override: None,
            site_evidence: None,
            tie_key: 0,
        }
    }

    fn make_params_ppm() -> SearchParams {
        use model::aa_set::AminoAcidSetBuilder;
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        SearchParams {
            aa_set,
            enzyme: model::enzyme::Enzyme::Trypsin,
            extra_enzymes: Vec::new(),
            min_length: 6,
            max_length: 40,
            max_missed_cleavages: 1,
            max_variable_mods_per_peptide: 3,
            max_label_sites: 0,
            precursor_tolerance: PrecursorTolerance::symmetric(Tolerance::Ppm(20.0)),
            charge_range: 2..=3,
            isotope_error_range: -1..=2,
            top_n_psms_per_spectrum: 10,
            num_tolerable_termini: 2,
            min_peaks: 10,
            precursor_cal_mode: search::PrecursorCalMode::Auto,
            cal_min_spec_keys: search::precursor_cal::constants::MIN_SPECKEYS_FOR_PREPASS,
            fragment_index_top_k: 0,
            precursor_mass_shift_ppm: 0.0,
            chimeric: false,
            dia: false,
            chimeric_isolation_halfwidth_da: 1.5,
            chimeric_max_coisolated: 2,
            chimeric_max_kl: 0.3,
            score_mode: search::ScoreMode::Rank,
            refine_select_psm_fdr: 0.01,
            require_nxst_sequon: false,
            open_search: None,
        }
    }

    fn parse_header(output: &[u8]) -> Vec<String> {
        let text = std::str::from_utf8(output).unwrap();
        let first_line = text.lines().next().unwrap_or("");
        first_line.split('\t').map(|s| s.to_string()).collect()
    }

    fn parse_rows(output: &[u8]) -> Vec<Vec<String>> {
        let text = std::str::from_utf8(output).unwrap();
        text.lines()
            .skip(1)
            .filter(|l| !l.is_empty())
            .map(|l| l.split('\t').map(|s| s.to_string()).collect())
            .collect()
    }

    #[test]
    fn tsv_header_columns_match_expected_when_mgf() {
        let params = make_params_ppm();
        let spectra: Vec<Spectrum> = vec![];
        let queues: Vec<TopNQueue> = vec![];
        let idx = make_empty_search_index();

        let mut buf = Vec::<u8>::new();
        let cands: Vec<search::candidate_gen::Candidate> = vec![];
        write_tsv_to(
            &mut buf, &spectra, &queues, &cands, &params, &idx, "test.mgf", true,
        )
        .unwrap();

        let cols = parse_header(&buf);
        assert_eq!(
            cols,
            vec![
                "#SpecFile",
                "SpecID",
                "ScanNum",
                "Title",
                "FragMethod",
                "Precursor",
                "IsotopeError",
                "PrecursorError(ppm)",
                "Charge",
                "Peptide",
                "Protein",
                "RankScore",
                "Modifications",
                "ProteinFamily",
                "EntryKind",
                "FeatureIds",
            ],
            "Header columns must match expected order when is_mgf=true"
        );
    }

    #[test]
    fn tsv_header_no_title_column_when_not_mgf() {
        let params = make_params_ppm();
        let spectra: Vec<Spectrum> = vec![];
        let queues: Vec<TopNQueue> = vec![];
        let idx = make_empty_search_index();

        let mut buf = Vec::<u8>::new();
        let cands: Vec<search::candidate_gen::Candidate> = vec![];
        write_tsv_to(
            &mut buf,
            &spectra,
            &queues,
            &cands,
            &params,
            &idx,
            "test.mzML",
            false,
        )
        .unwrap();

        let cols = parse_header(&buf);
        assert!(
            !cols.contains(&"Title".to_string()),
            "Title column must be absent when is_mgf=false"
        );
        assert!(cols.contains(&"ScanNum".to_string()));
        assert!(cols.contains(&"SpecID".to_string()));
    }

    #[test]
    fn tsv_handles_empty_queues_gracefully() {
        let params = make_params_ppm();
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];
        let queues = vec![TopNQueue::new(10)];
        let idx = make_empty_search_index();

        let mut buf = Vec::<u8>::new();
        let cands: Vec<search::candidate_gen::Candidate> = vec![];
        write_tsv_to(
            &mut buf, &spectra, &queues, &cands, &params, &idx, "test.mgf", true,
        )
        .unwrap();

        let rows = parse_rows(&buf);
        assert!(rows.is_empty(), "empty queue should produce no data rows");
    }

    #[test]
    fn tsv_writes_one_row_per_psm_in_rank_order() {
        let params = make_params_ppm();
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        let mut queue = TopNQueue::new(10);
        for (cidx, score) in [(0u32, 10.0f32), (1, 8.0), (2, 6.0)] {
            let mut psm = make_psm(0, score, score);
            psm.candidate_idxs = vec![cidx];
            queue.push(psm);
        }
        let queues = vec![queue];
        let idx = make_empty_search_index();

        let mut buf = Vec::<u8>::new();
        let cands: Vec<Candidate> = [b'A', b'G', b'S']
            .iter()
            .map(|&r| {
                let mut c = make_candidate(0, false);
                c.peptide = Peptide::new(vec![AminoAcid::standard(r).unwrap()], b'K', b'S');
                c
            })
            .collect();
        write_tsv_to(
            &mut buf, &spectra, &queues, &cands, &params, &idx, "test.mgf", true,
        )
        .unwrap();

        let rows = parse_rows(&buf);
        assert_eq!(rows.len(), 3, "should have 3 data rows");

        let scores: Vec<&str> = rows.iter().map(|r| r[11].as_str()).collect();

        assert_eq!(scores[0], "10", "first row should be the highest-score PSM");
        assert_eq!(scores[2], "6", "last row should be the lowest-score PSM");
    }

    #[test]
    fn tsv_peptide_column_includes_mods() {
        let params = make_params_ppm();
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        let m_unmod = AminoAcid::standard(b'M').unwrap();
        let ox_mod = Modification {
            name: "Oxidation".to_string(),
            mass_delta: 15.99491,
            residue: model::modification::ResidueSpec::Specific(b'M'),
            location: model::modification::ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        let m_ox = AminoAcid {
            residue: b'M',
            mass: m_unmod.mass,
            mod_: Some(model::modification::leak_mod(ox_mod)),
        };
        let a = AminoAcid::standard(b'A').unwrap();
        let peptide = Peptide::new(vec![a.clone(), m_ox, a], b'K', b'S');

        let psm = PsmMatch {
            spectrum_idx: 0,
            candidate_idxs: vec![0],
            charge_used: 2,
            mass_error_ppm: 0.0,
            score: 10.0,
            rank_score: 10.0,
            edge_score: 0,
            activation_method: None,
            features: search::psm::PsmFeatures::default(),
            isotope_offset: 0,
            precursor_mz_override: None,
            site_evidence: None,
            tie_key: 0,
        };

        let mut queue = TopNQueue::new(10);
        queue.push(psm);
        let queues = vec![queue];
        let idx = make_empty_search_index();

        let mut buf = Vec::<u8>::new();
        let cand = Candidate {
            peptide,
            protein_index: 0,
            start_offset_in_protein: 0,
            is_decoy: false,
            is_protein_n_term: false,
            is_protein_c_term: false,
        };
        let cands = vec![cand];
        write_tsv_to(
            &mut buf, &spectra, &queues, &cands, &params, &idx, "test.mgf", true,
        )
        .unwrap();

        let rows = parse_rows(&buf);
        assert_eq!(rows.len(), 1);
        let peptide_col = &rows[0][9];
        assert!(
            peptide_col.contains("+15.99"),
            "peptide column should contain oxidation mod delta (+15.99...), got: {}",
            peptide_col
        );
    }

    #[test]
    fn tsv_writes_real_accession_when_search_index_provided() {
        let accession = "sp|P02769|ALBU_BOVIN";
        let idx = make_search_index(accession);

        let params = make_params_ppm();
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        let psm = make_psm(0, 10.0, 10.0);
        let mut queue = TopNQueue::new(10);
        queue.push(psm);
        let queues = vec![queue];

        let mut buf = Vec::<u8>::new();
        let cands = vec![make_candidate(0, false)];
        write_tsv_to(
            &mut buf, &spectra, &queues, &cands, &params, &idx, "test.mgf", true,
        )
        .unwrap();

        let cols = parse_header(&buf);
        let rows = parse_rows(&buf);
        assert_eq!(rows.len(), 1);

        let prot_col = cols
            .iter()
            .position(|c| c == "Protein")
            .expect("Protein column missing");
        assert_eq!(
            rows[0][prot_col], accession,
            "Protein column should contain the real accession, not a PROT_N placeholder"
        );
    }

    #[test]
    fn tsv_ends_with_provenance_of_the_named_protein() {
        let accession = "sp|P02769|ALBU_BOVIN";
        let idx = make_search_index(accession);
        let params = make_params_ppm();
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];
        let mut q = TopNQueue::new(10);
        q.push(make_psm(0, 10.0, 10.0));
        let queues = vec![q];
        let mut buf = Vec::<u8>::new();
        let cands = vec![make_candidate(0, false)];
        write_tsv_to(
            &mut buf, &spectra, &queues, &cands, &params, &idx, "test.mgf", true,
        )
        .unwrap();
        let cols = parse_header(&buf);
        assert_eq!(
            &cols[cols.len() - 4..],
            ["Modifications", "ProteinFamily", "EntryKind", "FeatureIds"]
        );
        let row = &parse_rows(&buf)[0];
        assert_eq!(row.len(), cols.len());
        let at = |name: &str| row[cols.iter().position(|c| c == name).unwrap()].as_str();
        assert_eq!(at("ProteinFamily"), accession);
        assert_eq!(at("EntryKind"), "reference");
        assert_eq!(at("FeatureIds"), "");
    }

    #[test]
    fn tsv_writes_decoy_prefix_for_decoy_protein() {
        let accession = "sp|P02769|ALBU_BOVIN";
        let idx = make_search_index(accession);

        let params = make_params_ppm();
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        let psm = make_psm(0, 10.0, 10.0);

        let mut queue = TopNQueue::new(10);
        queue.push(psm);
        let queues = vec![queue];
        let cands = vec![make_candidate(1, true)];

        let mut buf = Vec::<u8>::new();
        write_tsv_to(
            &mut buf, &spectra, &queues, &cands, &params, &idx, "test.mgf", true,
        )
        .unwrap();

        let cols = parse_header(&buf);
        let rows = parse_rows(&buf);
        assert_eq!(rows.len(), 1);

        let prot_col = cols
            .iter()
            .position(|c| c == "Protein")
            .expect("Protein column missing");
        let expected_decoy = format!("XXX_{}", accession);
        assert_eq!(
            rows[0][prot_col], expected_decoy,
            "Protein column should carry decoy prefix for decoy PSM"
        );
    }

    fn make_mod(accession: Option<&str>) -> Modification {
        Modification {
            name: "TestMod".to_string(),
            mass_delta: 1.0,
            residue: ResidueSpec::Specific(b'A'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: accession.map(|s| s.to_string()),
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        }
    }

    #[test]
    fn tsv_modifications_column_precedes_provenance() {
        let params = make_params_ppm();
        let spectra: Vec<Spectrum> = vec![];
        let queues: Vec<TopNQueue> = vec![];
        let idx = make_empty_search_index();

        let mut buf = Vec::<u8>::new();
        let cands: Vec<search::candidate_gen::Candidate> = vec![];
        write_tsv_to(
            &mut buf, &spectra, &queues, &cands, &params, &idx, "test.mgf", true,
        )
        .unwrap();

        let cols = parse_header(&buf);
        assert_eq!(
            &cols[cols.len() - 5..],
            [
                "RankScore",
                "Modifications",
                "ProteinFamily",
                "EntryKind",
                "FeatureIds"
            ],
            "the provenance columns are appended after Modifications"
        );
    }

    #[test]
    fn tsv_modifications_field_single_accession_mod() {
        let params = make_params_ppm();
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        let plain_aa = AminoAcid::standard(b'A').unwrap();
        let mod_aa = AminoAcid {
            residue: b'A',
            mass: plain_aa.mass,
            mod_: Some(model::modification::leak_mod(make_mod(Some("UNIMOD:393")))),
        };
        let residues = vec![
            plain_aa.clone(),
            plain_aa.clone(),
            plain_aa.clone(),
            plain_aa.clone(),
            plain_aa.clone(),
            mod_aa,
        ];
        let peptide = Peptide::new(residues, b'K', b'S');
        let cand = Candidate {
            peptide,
            protein_index: 0,
            start_offset_in_protein: 0,
            is_decoy: false,
            is_protein_n_term: false,
            is_protein_c_term: false,
        };

        let psm = make_psm(0, 10.0, 10.0);
        let mut queue = TopNQueue::new(10);
        queue.push(psm);
        let queues = vec![queue];
        let idx = make_empty_search_index();

        let mut buf = Vec::<u8>::new();
        let cands = vec![cand];
        write_tsv_to(
            &mut buf, &spectra, &queues, &cands, &params, &idx, "test.mgf", true,
        )
        .unwrap();

        let cols = parse_header(&buf);
        let rows = parse_rows(&buf);
        assert_eq!(rows.len(), 1);

        let mod_col = cols
            .iter()
            .position(|c| c == "Modifications")
            .expect("Modifications column missing");
        assert_eq!(
            rows[0][mod_col], "6:UNIMOD:393",
            "single mod at 0-based 5 should emit '6:UNIMOD:393', got: {}",
            rows[0][mod_col]
        );
    }

    #[test]
    fn tsv_modifications_field_empty_when_no_accession_mod() {
        let params = make_params_ppm();
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        let plain_aa = AminoAcid::standard(b'A').unwrap();
        let noaccession_aa = AminoAcid {
            residue: b'A',
            mass: plain_aa.mass,
            mod_: Some(model::modification::leak_mod(make_mod(None))),
        };
        let residues = vec![plain_aa.clone(), noaccession_aa, plain_aa];
        let peptide = Peptide::new(residues, b'K', b'S');
        let cand = Candidate {
            peptide,
            protein_index: 0,
            start_offset_in_protein: 0,
            is_decoy: false,
            is_protein_n_term: false,
            is_protein_c_term: false,
        };

        let psm = make_psm(0, 10.0, 10.0);
        let mut queue = TopNQueue::new(10);
        queue.push(psm);
        let queues = vec![queue];
        let idx = make_empty_search_index();

        let mut buf = Vec::<u8>::new();
        let cands = vec![cand];
        write_tsv_to(
            &mut buf, &spectra, &queues, &cands, &params, &idx, "test.mgf", true,
        )
        .unwrap();

        let cols = parse_header(&buf);
        let rows = parse_rows(&buf);
        assert_eq!(rows.len(), 1);

        let mod_col = cols
            .iter()
            .position(|c| c == "Modifications")
            .expect("Modifications column missing");
        assert_eq!(
            rows[0][mod_col], "",
            "no-accession mod should emit empty Modifications field, got: {:?}",
            rows[0][mod_col]
        );
    }

    #[test]
    fn tsv_modifications_field_multiple_mods_ascending_order() {
        let params = make_params_ppm();
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        let plain_aa = AminoAcid::standard(b'A').unwrap();
        let mod1 = AminoAcid {
            residue: b'A',
            mass: plain_aa.mass,
            mod_: Some(model::modification::leak_mod(make_mod(Some("UNIMOD:4")))),
        };
        let mod2 = AminoAcid {
            residue: b'A',
            mass: plain_aa.mass,
            mod_: Some(model::modification::leak_mod(make_mod(Some("UNIMOD:393")))),
        };
        let residues = vec![
            plain_aa.clone(),
            mod1,
            plain_aa.clone(),
            plain_aa.clone(),
            plain_aa.clone(),
            mod2,
            plain_aa.clone(),
        ];
        let peptide = Peptide::new(residues, b'K', b'S');
        let cand = Candidate {
            peptide,
            protein_index: 0,
            start_offset_in_protein: 0,
            is_decoy: false,
            is_protein_n_term: false,
            is_protein_c_term: false,
        };

        let psm = make_psm(0, 10.0, 10.0);
        let mut queue = TopNQueue::new(10);
        queue.push(psm);
        let queues = vec![queue];
        let idx = make_empty_search_index();

        let mut buf = Vec::<u8>::new();
        let cands = vec![cand];
        write_tsv_to(
            &mut buf, &spectra, &queues, &cands, &params, &idx, "test.mgf", true,
        )
        .unwrap();

        let cols = parse_header(&buf);
        let rows = parse_rows(&buf);
        assert_eq!(rows.len(), 1);

        let mod_col = cols
            .iter()
            .position(|c| c == "Modifications")
            .expect("Modifications column missing");
        assert_eq!(
            rows[0][mod_col], "2:UNIMOD:4;6:UNIMOD:393",
            "multiple mods should be ';'-joined ascending, got: {}",
            rows[0][mod_col]
        );
    }
}
