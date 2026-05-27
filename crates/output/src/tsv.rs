//! TSV output writer.
//!
//! # Column order
//!
//! ```text
//! #SpecFile  SpecID  ScanNum  [Title — only when is_mgf]  FragMethod
//! Precursor  IsotopeError  PrecursorError(ppm|Da)  Charge
//! Peptide  Protein  DeNovoScore  MSGFScore  SpecEValue  EValue
//! ```
//!
//! # Column semantics
//!
//! * **FragMethod**: `ActivationMethod::name()` for the five canonical variants;
//!   `"UNKNOWN"` for unknown / unset activation.
//! * **IsotopeError**: winning isotope offset from the search (`PsmMatch::isotope_offset`).
//! * **Decoy filtering**: decoys are emitted; downstream Percolator labels them.
//! * **SpecID for non-MGF**: `"scan=N"` (mzML convention).

use std::io::{self, BufWriter, Write};

use crate::row_context::{iter_ranked, RowContext};
use search::candidate_gen::Candidate;
use search::psm::{PsmMatch, TopNQueue};
use search::search_index::SearchIndex;
use search::search_params::SearchParams;
use model::spectrum::Spectrum;
use model::tolerance::Tolerance;

// ── public API ──────────────────────────────────────────────────────────────

/// Write all PSMs to a tab-separated file at `output_path`.
///
/// `spectra` and `queues` must be parallel slices (same length): `queues[i]`
/// holds the top-N PSMs for `spectra[i]`.
///
/// `search_index` is used to resolve protein accessions from
/// `psm.candidate.protein_index`.  Decoy accessions already carry the prefix
/// (set by `target_plus_decoy`) — no prefix arithmetic is needed here.
///
/// `spec_file_name` is the bare filename (e.g. `"test.mgf"`) written in the
/// `#SpecFile` column.
///
/// `is_mgf` controls whether a `Title` column is emitted in the header and
/// rows, matching Java's behaviour for MGF vs mzML input.
#[allow(clippy::too_many_arguments, reason = "Writer API mirrors PIN writer; grouping into a struct would diverge from the parallel write_pin API")]
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
    let file = std::fs::File::create(output_path)?;
    let mut writer = BufWriter::new(file);
    write_tsv_to(&mut writer, spectra, queues, candidates, params, search_index, spec_file_name, is_mgf)
}

/// Write all PSMs to an arbitrary writer — useful for testing without temp
/// files.
///
/// See [`write_tsv`] for parameter documentation.
#[allow(clippy::too_many_arguments, reason = "Writer API mirrors PIN writer; grouping into a struct would diverge from the parallel write_pin API")]
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
        write_spectrum_rows(writer, spec, queue, candidates, params, spec_file_name, is_mgf, search_index)?;
    }
    Ok(())
}

// ── header ───────────────────────────────────────────────────────────────────

fn write_header<W: Write>(
    writer: &mut W,
    params: &SearchParams,
    is_mgf: bool,
) -> io::Result<()> {
    let ppm_mode = matches!(params.precursor_tolerance.left, Tolerance::Ppm(_));
    let prec_err_col = if ppm_mode { "PrecursorError(ppm)" } else { "PrecursorError(Da)" };

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
        "DeNovoScore",
        "MSGFScore",
        "SpecEValue",
        "EValue",
    ]);

    writeln!(writer, "{}", cols.join("\t"))
}

// ── per-spectrum rows ─────────────────────────────────────────────────────────

/// Row-writing context: fixed fields derived once per spectrum.
struct RowCtx<'a> {
    spec_file_name: &'a str,
    is_mgf: bool,
    ppm_mode: bool,
}

#[allow(clippy::too_many_arguments, reason = "Writer API mirrors PIN writer; grouping into a struct would diverge from the parallel write_pin API")]
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
    // Sort best-first (lowest spec_e_value first).
    let psms = queue.clone().into_sorted_vec();

    let row_ctx = RowCtx {
        spec_file_name,
        is_mgf,
        ppm_mode: matches!(params.precursor_tolerance.left, Tolerance::Ppm(_)),
    };

    for (_rank, psm) in iter_ranked(&psms) {
        let cand = &candidates[psm.primary_candidate_idx() as usize];
        let ctx = RowContext::new(spec, cand, search_index);
        write_psm_row(writer, spec, psm, cand, &ctx, &row_ctx)?;
    }
    Ok(())
}

fn write_psm_row<W: Write>(
    writer: &mut W,
    spec: &Spectrum,
    psm: &PsmMatch,
    cand: &Candidate,
    ctx: &RowContext,
    row_ctx: &RowCtx<'_>,
) -> io::Result<()> {
    let is_mgf = row_ctx.is_mgf;
    let ppm_mode = row_ctx.ppm_mode;
    let spec_file_name = row_ctx.spec_file_name;

    // SpecID: derived from RowContext (title if non-empty, else "scan=N")
    let spec_id = &ctx.spec_id;

    let scan_num = ctx.scan;

    // FragMethod: use ActivationMethod::name() for known variants, "UNKNOWN" for None
    let frag_method = psm
        .activation_method
        .map(|m| m.name().to_string())
        .unwrap_or_else(|| "UNKNOWN".to_string());

    // Precursor m/z formatted to 4 decimal places
    let precursor = format!("{:.4}", spec.precursor_mz);

    // IsotopeError: winning isotope offset from the search (matches PIN column).
    let isotope_error: i32 = psm.isotope_offset as i32;

    // PrecursorError: mass_error_ppm stored on psm; convert to Da if needed
    let precursor_error = if ppm_mode {
        format!("{:.4}", psm.mass_error_ppm)
    } else {
        // Convert ppm error back to Da using precursor_mz
        let da = psm.mass_error_ppm * 1e-6 * spec.precursor_mz;
        format!("{:.4}", da)
    };

    // Charge
    let charge = psm.charge_used;

    // Peptide: uses the existing Display impl → "pre.SEQ_WITH_MODS.post"
    let peptide = &cand.peptide;
    let protein = &ctx.accession;

    // DeNovoScore
    let de_novo_score = psm.de_novo_score;

    // MSGFScore: integer-rounded raw score
    let msgf_score = psm.score.round() as i32;

    // SpecEValue: format as scientific notation with 6 decimal places
    let spec_e_value = format_e_value(psm.spec_e_value);

    // EValue: same formatting
    let e_value = format_e_value(psm.e_value);

    // Build row
    if is_mgf {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            spec_file_name,
            spec_id,
            scan_num,
            spec.title,   // Title column (MGF only)
            frag_method,
            precursor,
            isotope_error,
            precursor_error,
            charge,
            peptide,
            protein,
            de_novo_score,
            msgf_score,
            spec_e_value,
            e_value,
        )
    } else {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
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
            de_novo_score,
            msgf_score,
            spec_e_value,
            e_value,
        )
    }
}

// ── helpers ───────────────────────────────────────────────────────────────────


/// Format a SpecEValue / EValue in scientific notation.
///
/// Matches Java's `%.6e` formatting: always lowercase `e`, 6 fractional digits.
fn format_e_value(v: f64) -> String {
    format!("{:.6e}", v)
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use model::amino_acid::AminoAcid;
    use search::candidate_gen::Candidate;
    use model::modification::Modification;
    use model::peptide::Peptide;
    use model::protein::{Protein, ProteinDb};
    use search::search_index::SearchIndex;
    use model::tolerance::PrecursorTolerance;

    // ── fixture helpers ─────────────────────────────────────────────────────

    /// Build a minimal `SearchIndex` with one target protein.
    fn make_search_index(accession: &str) -> SearchIndex {
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: accession.to_string(),
                description: String::new(),
                sequence: b"MKWVTFISLL".to_vec(),
            }],
        };
        SearchIndex::from_target_db(&target, "XXX_")
    }

    /// Build an empty `SearchIndex` for tests that don't inspect protein values.
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
        }
    }

    /// Build a single Candidate fixture. Mirrors the make_candidate in pin.rs.
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

    fn make_psm(spectrum_idx: usize, score: f32, spec_e_value: f64) -> PsmMatch {
        PsmMatch {
            spectrum_idx,
            candidate_idxs: vec![0],
            charge_used: 2,
            mass_error_ppm: 1.5,
            score,
            rank_score: score,  // iter33: test fixtures default rank_score = score
            edge_score: 0,
            spec_e_value,
            de_novo_score: 42,
            activation_method: Some(model::activation::ActivationMethod::HCD),
            e_value: spec_e_value * 100.0,
            features: search::psm::PsmFeatures::default(),
            isotope_offset: 0,
        }
    }

    fn make_params_ppm() -> SearchParams {
        use model::aa_set::AminoAcidSetBuilder;
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        SearchParams {
            aa_set,
            enzyme: model::enzyme::Enzyme::Trypsin,
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
            precursor_cal_mode: search::PrecursorCalMode::Auto,
            precursor_mass_shift_ppm: 0.0,
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
            .skip(1) // skip header
            .filter(|l| !l.is_empty())
            .map(|l| l.split('\t').map(|s| s.to_string()).collect())
            .collect()
    }

    // ── Test 1: header columns match expected when MGF ─────────────────────

    #[test]
    fn tsv_header_columns_match_expected_when_mgf() {
        let params = make_params_ppm();
        let spectra: Vec<Spectrum> = vec![];
        let queues: Vec<TopNQueue> = vec![];
        let idx = make_empty_search_index();

        let mut buf = Vec::<u8>::new();
        let cands: Vec<search::candidate_gen::Candidate> = vec![];
        write_tsv_to(&mut buf, &spectra, &queues, &cands, &params, &idx, "test.mgf", true).unwrap();

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
                "DeNovoScore",
                "MSGFScore",
                "SpecEValue",
                "EValue",
            ],
            "Header columns must match expected order when is_mgf=true"
        );
    }

    // ── Test 2: header omits Title when not MGF ────────────────────────────

    #[test]
    fn tsv_header_no_title_column_when_not_mgf() {
        let params = make_params_ppm();
        let spectra: Vec<Spectrum> = vec![];
        let queues: Vec<TopNQueue> = vec![];
        let idx = make_empty_search_index();

        let mut buf = Vec::<u8>::new();
        let cands: Vec<search::candidate_gen::Candidate> = vec![];
        write_tsv_to(&mut buf, &spectra, &queues, &cands, &params, &idx, "test.mzML", false).unwrap();

        let cols = parse_header(&buf);
        assert!(!cols.contains(&"Title".to_string()), "Title column must be absent when is_mgf=false");
        assert!(cols.contains(&"ScanNum".to_string()));
        assert!(cols.contains(&"SpecID".to_string()));
    }

    // ── Test 3: empty queues → only header, no data rows ──────────────────

    #[test]
    fn tsv_handles_empty_queues_gracefully() {
        let params = make_params_ppm();
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];
        let queues = vec![TopNQueue::new(10)]; // empty queue
        let idx = make_empty_search_index();

        let mut buf = Vec::<u8>::new();
        let cands: Vec<search::candidate_gen::Candidate> = vec![];
        write_tsv_to(&mut buf, &spectra, &queues, &cands, &params, &idx, "test.mgf", true).unwrap();

        let rows = parse_rows(&buf);
        assert!(rows.is_empty(), "empty queue should produce no data rows");
    }

    // ── Test 4: PSMs written in rank order (best spec_e_value first) ───────

    #[test]
    fn tsv_writes_one_row_per_psm_in_rank_order() {
        let params = make_params_ppm();
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        let mut queue = TopNQueue::new(10);
        // Push 3 PSMs with descending spec_e_values (best = smallest)
        queue.push(make_psm(0, 10.0, 1e-10)); // best (rank 1)
        queue.push(make_psm(0, 8.0,  1e-8));  // middle (rank 2)
        queue.push(make_psm(0, 6.0,  1e-6));  // worst (rank 3)
        let queues = vec![queue];
        let idx = make_empty_search_index();

        let mut buf = Vec::<u8>::new();
        let cands = vec![make_candidate(0, false)];
        write_tsv_to(&mut buf, &spectra, &queues, &cands, &params, &idx, "test.mgf", true).unwrap();

        let rows = parse_rows(&buf);
        assert_eq!(rows.len(), 3, "should have 3 data rows");

        // Extract SpecEValue column (index 13 when is_mgf=true: 0=#SpecFile 1=SpecID
        // 2=ScanNum 3=Title 4=FragMethod 5=Precursor 6=IsotopeError 7=PrecursorError
        // 8=Charge 9=Peptide 10=Protein 11=DeNovoScore 12=MSGFScore 13=SpecEValue)
        let spec_evalues: Vec<&str> = rows.iter().map(|r| r[13].as_str()).collect();

        // Best PSM (1e-10) should come first
        assert!(
            spec_evalues[0].contains("1.000000e") && spec_evalues[0].contains("-10"),
            "first row should have spec_e_value 1e-10, got: {}",
            spec_evalues[0]
        );
        assert!(
            spec_evalues[2].contains("1.000000e") && spec_evalues[2].contains("-6"),
            "last row should have spec_e_value 1e-6, got: {}",
            spec_evalues[2]
        );
    }

    // ── Test 5: peptide column includes mods ───────────────────────────────

    #[test]
    fn tsv_peptide_column_includes_mods() {
        let params = make_params_ppm();
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        // Build a peptide with an oxidized methionine (+15.99491 Da)
        let m_unmod = AminoAcid::standard(b'M').unwrap();
        let ox_mod = Modification {
            name: "Oxidation".to_string(),
            mass_delta: 15.99491,
            residue: model::modification::ResidueSpec::Specific(b'M'),
            location: model::modification::ModLocation::Anywhere,
            fixed: false,
            accession: None,
        };
        let m_ox = AminoAcid {
            residue: b'M',
            mass: m_unmod.mass,
            mod_: Some(std::sync::Arc::new(ox_mod)),
        };
        let a = AminoAcid::standard(b'A').unwrap();
        // Peptide: K.AM(ox)A.S
        let peptide = Peptide::new(vec![a.clone(), m_ox, a], b'K', b'S');

        let psm = PsmMatch {
            spectrum_idx: 0,
            candidate_idxs: vec![0],
            charge_used: 2,
            mass_error_ppm: 0.0,
            score: 10.0,
            rank_score: 10.0,
            edge_score: 0,
            spec_e_value: 1e-5,
            de_novo_score: 0,
            activation_method: None,
            e_value: 1e-3,
            features: search::psm::PsmFeatures::default(),
            isotope_offset: 0,
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
        write_tsv_to(&mut buf, &spectra, &queues, &cands, &params, &idx, "test.mgf", true).unwrap();

        let rows = parse_rows(&buf);
        assert_eq!(rows.len(), 1);
        // Peptide column is index 9 (0=#SpecFile 1=SpecID 2=ScanNum 3=Title
        // 4=FragMethod 5=Precursor 6=IsotopeError 7=PrecursorError 8=Charge
        // 9=Peptide)
        let peptide_col = &rows[0][9];
        assert!(
            peptide_col.contains("+15.99"),
            "peptide column should contain oxidation mod delta (+15.99...), got: {}",
            peptide_col
        );
    }

    // ── Test 6: real accession emitted for target PSM ─────────────────────────

    #[test]
    fn tsv_writes_real_accession_when_search_index_provided() {
        let accession = "sp|P02769|ALBU_BOVIN";
        let idx = make_search_index(accession);

        let params = make_params_ppm();
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        // protein_index = 0 → first target protein
        let psm = make_psm(0, 10.0, 1e-5);
        let mut queue = TopNQueue::new(10);
        queue.push(psm);
        let queues = vec![queue];

        let mut buf = Vec::<u8>::new();
        let cands = vec![make_candidate(0, false)];
        write_tsv_to(&mut buf, &spectra, &queues, &cands, &params, &idx, "test.mgf", true).unwrap();

        let cols = parse_header(&buf);
        let rows = parse_rows(&buf);
        assert_eq!(rows.len(), 1);

        let prot_col = cols.iter().position(|c| c == "Protein").expect("Protein column missing");
        assert_eq!(
            rows[0][prot_col], accession,
            "Protein column should contain the real accession, not a PROT_N placeholder"
        );
    }

    // ── Test 7: decoy accession carries decoy prefix ──────────────────────────

    #[test]
    fn tsv_writes_decoy_prefix_for_decoy_protein() {
        let accession = "sp|P02769|ALBU_BOVIN";
        let idx = make_search_index(accession);

        let params = make_params_ppm();
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        // SearchIndex: 1 target (idx 0) + 1 decoy (idx 1, accession = "XXX_<base>")
        let psm = make_psm(0, 10.0, 1e-5);

        let mut queue = TopNQueue::new(10);
        queue.push(psm);
        let queues = vec![queue];
        let cands = vec![make_candidate(1, true)]; // decoy candidate at protein_index 1

        let mut buf = Vec::<u8>::new();
        write_tsv_to(&mut buf, &spectra, &queues, &cands, &params, &idx, "test.mgf", true).unwrap();

        let cols = parse_header(&buf);
        let rows = parse_rows(&buf);
        assert_eq!(rows.len(), 1);

        let prot_col = cols.iter().position(|c| c == "Protein").expect("Protein column missing");
        let expected_decoy = format!("XXX_{}", accession);
        assert_eq!(
            rows[0][prot_col], expected_decoy,
            "Protein column should carry decoy prefix for decoy PSM"
        );
    }
}
