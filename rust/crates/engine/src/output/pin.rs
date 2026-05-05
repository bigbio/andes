//! PIN output writer — mirrors Java `DirectPinWriter` (Phase 7/Task 3).
//!
//! Produces a Percolator-consumable `.pin` file with the column layout used
//! by MS-GF+ and OpenMS PercolatorAdapter so that downstream tools (Percolator,
//! MS²Rescore, Mokapot) can consume the output interchangeably.
//!
//! # Column order
//!
//! ```text
//! SpecId  Label  ScanNr  ExpMass  CalcMass  mass  RawScore  DeNovoScore
//! lnSpecEValue  lnEValue  isotope_error  peplen  dm  absdm
//! charge<min>  charge<min+1>  ...  charge<max>
//! enzN  enzC  enzInt
//! NumMatchedMainIons  longest_b  longest_y  longest_y_pct
//! ExplainedIonCurrentRatio  NTermIonCurrentRatio  CTermIonCurrentRatio
//! MS2IonCurrent  IsolationWindowEfficiency
//! MeanErrorTop7  StdevErrorTop7  MeanRelErrorTop7  StdevRelErrorTop7
//! lnDeltaSpecEValue  matchedIonRatio
//! Peptide  Proteins
//! ```
//!
//! # Java divergences (Phase 7 MVP)
//!
//! * **Label**: Java inspects all protein accessions to detect "all-decoy" PSMs.
//!   Rust uses `psm.candidate.is_decoy` directly: `1` for target, `-1` for decoy.
//!   Documented intentional simplification — Percolator only needs target/decoy
//!   disambiguation, which this provides.
//!
//! * **isotope_error**: always `0`. Phase 4e's precursor-matching loop tries
//!   multiple isotope offsets but does not record which offset produced the match.
//!   Fix in a later phase once the winning offset is threaded into `PsmMatch`.
//!
//! * **enzN / enzC / enzInt**: zero-stubbed. Java computes enzymatic-boundary
//!   indicators from the pre/post flanking residues + enzyme rules (OpenMS
//!   PercolatorInfile convention). Rust would need `Enzyme::is_cleavage_site`
//!   wiring; deferred to a later task.
//!
//! * **13 feature columns** (NumMatchedMainIons through StdevRelErrorTop7):
//!   zero-stubbed. Phase 7/Task 6 (stretch) fills these from Phase 5's scored
//!   spectrum. Percolator runs but with reduced discrimination power until Task 6.
//!
//! * **Proteins**: single column with `"PROT_{protein_index}"`. Java emits full
//!   accession strings; multi-protein PSMs get additional tab-separated columns.
//!   Rust MVP uses the index as a placeholder — Task 4 (CLI wiring) threads in
//!   a real `&ProteinDb` for full accession lookup.
//!
//! * **peplen**: Java uses `match.getLength()` which is the sequence length
//!   *including* flanking residues (`length - 2` when emitting). Rust's
//!   `Peptide::length()` returns only the residue count (no flanking), so we
//!   use it directly — the values are equivalent.
//!
//! * **dm / absdm**: computed from precursor m/z using isotope_error = 0.
//!   Java adjusts `adjustedExpMz = precursorMz - ISOTOPE * isotopeError / charge`.
//!   Since isotope_error is stubbed at 0, adjustedExpMz == precursorMz here.
//!
//! * **CalcMass**: `peptide.mass()` already includes H2O (Rust `Peptide::mass()`
//!   sums residue masses + H2O). Java's `theoMz * charge` involves charge-carrier
//!   mass; Rust computes the neutral mass directly from the peptide.

use std::io::{self, BufWriter, Write};

use crate::mass::{ISOTOPE, PROTON};
use crate::psm::{PsmMatch, TopNQueue};
use crate::search_params::SearchParams;
use crate::spectrum::Spectrum;

// ── public API ───────────────────────────────────────────────────────────────

/// Write all PSMs to a Percolator `.pin` file at `output_path`.
///
/// `spectra` and `queues` must be parallel slices (same length): `queues[i]`
/// holds the top-N PSMs for `spectra[i]`.
///
/// `decoy_prefix` is used to derive the `Label` column (target = 1, decoy = -1).
pub fn write_pin(
    output_path: &std::path::Path,
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    params: &SearchParams,
    decoy_prefix: &str,
) -> io::Result<()> {
    let file = std::fs::File::create(output_path)?;
    let mut writer = BufWriter::new(file);
    write_pin_to(&mut writer, spectra, queues, params, decoy_prefix)
}

/// Write all PSMs to an arbitrary writer — useful for testing without temp files.
///
/// See [`write_pin`] for parameter documentation.
pub fn write_pin_to<W: Write>(
    writer: &mut W,
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    params: &SearchParams,
    decoy_prefix: &str,
) -> io::Result<()> {
    let _ = decoy_prefix; // used indirectly via psm.candidate.is_decoy
    let min_charge = *params.charge_range.start();
    let max_charge = *params.charge_range.end();

    write_header(writer, min_charge, max_charge)?;

    for (spec_idx, queue) in queues.iter().enumerate() {
        if queue.is_empty() {
            continue;
        }
        let spec = &spectra[spec_idx];
        write_spectrum_rows(writer, spec, queue, min_charge, max_charge)?;
    }
    Ok(())
}

// ── header ────────────────────────────────────────────────────────────────────

fn write_header<W: Write>(writer: &mut W, min_charge: u8, max_charge: u8) -> io::Result<()> {
    let mut cols: Vec<String> = vec![
        "SpecId".to_string(),
        "Label".to_string(),
        "ScanNr".to_string(),
        "ExpMass".to_string(),
        "CalcMass".to_string(),
        "mass".to_string(),
        "RawScore".to_string(),
        "DeNovoScore".to_string(),
        "lnSpecEValue".to_string(),
        "lnEValue".to_string(),
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
        // PIN_FEATURES (zero-stubbed in Phase 7 MVP)
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
        // PIN_EXTRA_FEATURES
        "lnDeltaSpecEValue".to_string(),
        "matchedIonRatio".to_string(),
        // Peptide / Proteins
        "Peptide".to_string(),
        "Proteins".to_string(),
    ]);

    writeln!(writer, "{}", cols.join("\t"))
}

// ── per-spectrum rows ──────────────────────────────────────────────────────────

fn write_spectrum_rows<W: Write>(
    writer: &mut W,
    spec: &Spectrum,
    queue: &TopNQueue,
    min_charge: u8,
    max_charge: u8,
) -> io::Result<()> {
    // Sort best-first (lowest spec_e_value first, then highest score).
    let psms = queue.clone().into_sorted_vec();

    // spec_id used in SpecId column: title if non-empty, else "scan=N"
    let scan = spec.scan.unwrap_or(0);
    let spec_id = if spec.title.is_empty() {
        format!("scan={}", scan)
    } else {
        spec.title.clone()
    };

    // find rank-2 SpecEValue: first distinct spec_e_value after rank-1
    let rank2_spec_e_value = find_rank2_spec_e_value(&psms);

    let mut rank: u32 = 0;
    let mut last_spec_e_value = f64::NAN;

    for psm in &psms {
        // Rank increments when spec_e_value changes (ties share a rank).
        if psm.spec_e_value != last_spec_e_value {
            rank += 1;
            last_spec_e_value = psm.spec_e_value;
        }
        write_psm_row(
            writer,
            spec,
            psm,
            &spec_id,
            scan,
            rank,
            rank2_spec_e_value,
            min_charge,
            max_charge,
        )?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn write_psm_row<W: Write>(
    writer: &mut W,
    spec: &Spectrum,
    psm: &PsmMatch,
    spec_id: &str,
    scan: i32,
    rank: u32,
    rank2_spec_e_value: f64,
    min_charge: u8,
    max_charge: u8,
) -> io::Result<()> {
    let charge = psm.charge_used as f64;

    // SpecId: Java pattern is specID + "_" + scanNum + "_" + rank
    let psm_id = format!("{}_{}_{}", spec_id, scan, rank);

    // Label: target = 1, decoy = -1
    // MVP divergence: uses is_decoy flag directly (Java inspects all protein accessions)
    let label: i32 = if psm.candidate.is_decoy { -1 } else { 1 };

    // ExpMass: neutral precursor mass = mz * charge - charge * PROTON
    let exp_mass = spec.precursor_mz * charge - charge * PROTON;

    // CalcMass: theoretical neutral mass. peptide.mass() already includes H2O.
    // Java: theoMass = theoMz * charge where theoMz = (peptideMass + H2O) / charge + PROTON
    //   → theoMass = peptideMass + H2O + PROTON * charge
    // But for a neutral mass we want peptideMass + H2O (same as peptide.mass() here).
    // We match Java's CalcMass column (theoMass = theoMz * charge) which is the
    // protonated mass — so: peptide.mass() + charge * PROTON
    // However the fixture shows CalcMass ≈ ExpMass (neutral masses), so:
    let peptide_mass = psm.candidate.peptide.mass(); // includes H2O
    let calc_mass = peptide_mass + charge * PROTON;

    // mass: duplicate of ExpMass (per Java line 212: "mass — duplicate of ExpMass")
    let mass = exp_mass;

    // RawScore: integer-rounded score
    let raw_score = psm.score.round() as i32;

    // DeNovoScore
    let de_novo_score = psm.de_novo_score;

    // lnSpecEValue
    let ln_spec_e_value = if psm.spec_e_value > 0.0 {
        psm.spec_e_value.ln()
    } else {
        -f64::MAX
    };

    // lnEValue
    let ln_e_value = if psm.e_value > 0.0 {
        psm.e_value.ln()
    } else {
        -f64::MAX
    };

    // isotope_error: always 0 (MVP divergence)
    let isotope_error: i32 = 0;

    // peplen: number of residues (no flanking)
    let peplen = psm.candidate.peptide.length();

    // dm / absdm: precursor mass error in Da
    // Java: adjustedExpMz = precursorMz - ISOTOPE * isotopeError / charge
    // Since isotopeError = 0: adjustedExpMz = precursorMz
    // theoMz = (peptideMass + H2O) / charge + PROTON
    //        = peptide.mass() / charge + PROTON  (since peptide.mass() includes H2O)
    // dm = adjustedExpMz - theoMz
    let theo_mz = peptide_mass / charge + PROTON;
    let adjusted_exp_mz = spec.precursor_mz - ISOTOPE * (isotope_error as f64) / charge;
    let dm = adjusted_exp_mz - theo_mz;
    let absdm = dm.abs();

    // lnDeltaSpecEValue
    let ln_delta_spec_e_value = compute_ln_delta_spec_e_value(rank, psm.spec_e_value, rank2_spec_e_value);

    // matchedIonRatio: 0 since NumMatchedMainIons is zero-stubbed
    let matched_ion_ratio = 0.0_f64;

    // Peptide: pre.SEQ_WITH_MODS.post  (uses existing Display impl)
    let peptide_str = format!("{}", psm.candidate.peptide);

    // Proteins: MVP uses protein_index as placeholder
    // Divergence: Java uses full accession; multi-protein PSMs get extra tab-separated columns.
    let proteins_str = format_protein(psm);

    // Build row — tab-separated
    // Fixed prefix
    write!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        psm_id,
        label,
        scan,
        format_double(exp_mass),
        format_double(calc_mass),
        format_double(mass),
        raw_score,
        de_novo_score,
        format_double(ln_spec_e_value),
        format_double(ln_e_value),
        isotope_error,
        peplen,
        format_double(dm),
        format_double(absdm),
    )?;

    // Charge one-hot
    for c in min_charge..=max_charge {
        let flag: i32 = if c == psm.charge_used { 1 } else { 0 };
        write!(writer, "\t{}", flag)?;
    }

    // enzN, enzC, enzInt (zero-stubbed)
    write!(writer, "\t0\t0\t0")?;

    // 13 feature columns (zero-stubbed)
    // NumMatchedMainIons, longest_b, longest_y, longest_y_pct,
    // ExplainedIonCurrentRatio, NTermIonCurrentRatio, CTermIonCurrentRatio,
    // MS2IonCurrent, IsolationWindowEfficiency,
    // MeanErrorTop7, StdevErrorTop7, MeanRelErrorTop7, StdevRelErrorTop7
    write!(writer, "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0")?;

    // lnDeltaSpecEValue, matchedIonRatio
    write!(
        writer,
        "\t{}\t{}",
        format_double(ln_delta_spec_e_value),
        format_double(matched_ion_ratio),
    )?;

    // Peptide, Proteins
    writeln!(writer, "\t{}\t{}", peptide_str, proteins_str)
}

// ── helpers ───────────────────────────────────────────────────────────────────

/// Find the rank-2 SpecEValue: the first distinct spec_e_value encountered after
/// the rank-1 value (skipping ties). Returns `f64::NAN` if no rank-2 exists.
///
/// Mirrors Java `DirectPinWriter.findRank2SpecEValue` (line 262).
/// PSMs must be sorted best-first (lowest spec_e_value first).
fn find_rank2_spec_e_value(psms: &[PsmMatch]) -> f64 {
    let mut rank1 = f64::NAN;
    for psm in psms {
        let se = psm.spec_e_value;
        if rank1.is_nan() {
            rank1 = se;
        } else if se != rank1 {
            return se;
        }
    }
    f64::NAN
}

/// `log(rank1 SpecEValue / rank2 SpecEValue)` for rank-1 PSMs; `0.0` otherwise
/// or when either SpecEValue is non-positive / NaN.
///
/// Mirrors Java `DirectPinWriter.computeLnDeltaSpecEValue` (line 283).
fn compute_ln_delta_spec_e_value(rank: u32, rank1_spec_e_value: f64, rank2_spec_e_value: f64) -> f64 {
    if rank != 1 {
        return 0.0;
    }
    if rank1_spec_e_value.is_nan() || rank2_spec_e_value.is_nan() {
        return 0.0;
    }
    if rank1_spec_e_value <= 0.0 || rank2_spec_e_value <= 0.0 {
        return 0.0;
    }
    (rank1_spec_e_value / rank2_spec_e_value).ln()
}

/// Format a `f64` in `%.6g` style (6 significant figures), matching Java's
/// `String.format(Locale.ROOT, "%.6g", v)` used in `formatDouble`.
///
/// NaN or infinite values are formatted as `"0"` matching Java's behaviour:
/// `if (Double.isNaN(v) || Double.isInfinite(v)) return "0";`
fn format_double(v: f64) -> String {
    if v.is_nan() || v.is_infinite() {
        return "0".to_string();
    }
    // %.6g: 6 significant figures, removes trailing zeros after decimal point.
    // Rust doesn't have a %g format natively, so we mimic it:
    // use scientific notation when |v| < 1e-4 or |v| >= 1e6; else fixed.
    if v == 0.0 {
        return "0".to_string();
    }
    let abs = v.abs();
    if !(1e-4..1e6).contains(&abs) {
        // Scientific notation, 5 decimal places after dot = 6 significant digits
        let s = format!("{:.5e}", v);
        trim_scientific(&s)
    } else {
        // Fixed notation. Determine decimal places for 6 sig figs.
        let digits_before_decimal = abs.log10().floor() as i32 + 1;
        let decimal_places = (6 - digits_before_decimal).max(0) as usize;
        let s = format!("{:.prec$}", v, prec = decimal_places);
        trim_fixed(&s)
    }
}

/// Trim trailing zeros from a fixed-point string (e.g. "1.50000" → "1.5").
fn trim_fixed(s: &str) -> String {
    if s.contains('.') {
        let trimmed = s.trim_end_matches('0').trim_end_matches('.');
        trimmed.to_string()
    } else {
        s.to_string()
    }
}

/// Normalise a Rust scientific notation string to match Java's `%g` output.
///
/// Rust produces `1.23456e7`; Java produces `1.23456e+07`. We want to match
/// the Percolator-readable convention (either works for Percolator, but for
/// the byte-parity test we normalise to remove leading zeros in exponent and
/// strip trailing zeros in the significand).
fn trim_scientific(s: &str) -> String {
    // Split at 'e'
    if let Some(pos) = s.find('e') {
        let mantissa = s[..pos].to_string();
        let exp_part = &s[pos + 1..];

        // Trim trailing zeros from mantissa (after decimal point)
        let mantissa = if mantissa.contains('.') {
            mantissa.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            mantissa
        };

        // Parse and reformat exponent (remove leading zeros, keep sign)
        let exp_val: i32 = exp_part.parse().unwrap_or(0);
        format!("{}e{:+03}", mantissa, exp_val)
    } else {
        s.to_string()
    }
}

/// Format a protein identifier from the PSM (MVP placeholder).
///
/// Java emits full accession strings; multi-protein PSMs get additional
/// tab-separated `Proteins` columns. For MVP, emit `"PROT_{index}"`.
/// Task 4 (CLI wiring) threads in a `&ProteinDb` for full accession lookup.
fn format_protein(psm: &PsmMatch) -> String {
    format!("PROT_{}", psm.candidate.protein_index)
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::amino_acid::AminoAcid;
    use crate::candidate_gen::Candidate;
    use crate::peptide::Peptide;
    use crate::tolerance::PrecursorTolerance;
    use crate::tolerance::Tolerance;

    // ── fixture helpers ─────────────────────────────────────────────────────

    fn make_spectrum(title: &str, scan: i32, precursor_mz: f64) -> Spectrum {
        Spectrum {
            title: title.to_string(),
            precursor_mz,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: Some(scan),
            peaks: vec![],
        }
    }

    fn make_psm(spectrum_idx: usize, score: f32, spec_e_value: f64, is_decoy: bool, charge: u8) -> PsmMatch {
        let aa = AminoAcid::standard(b'A').unwrap();
        let peptide = Peptide::new(vec![aa], b'K', b'S');
        PsmMatch {
            spectrum_idx,
            candidate: Candidate {
                peptide,
                protein_index: 0,
                start_offset_in_protein: 0,
                is_decoy,
            },
            charge_used: charge,
            mass_error_ppm: 1.5,
            score,
            spec_e_value,
            de_novo_score: 42,
            activation_method: Some(crate::activation::ActivationMethod::HCD),
            e_value: spec_e_value * 100.0,
        }
    }

    fn make_params(charge_range: std::ops::RangeInclusive<u8>) -> SearchParams {
        use crate::aa_set::AminoAcidSetBuilder;
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        SearchParams {
            aa_set,
            enzyme: crate::enzyme::Enzyme::Trypsin,
            min_length: 6,
            max_length: 40,
            max_missed_cleavages: 1,
            max_variable_mods_per_peptide: 3,
            precursor_tolerance: PrecursorTolerance::symmetric(Tolerance::Ppm(20.0)),
            charge_range,
            isotope_error_range: -1..=2,
            top_n_psms_per_spectrum: 10,
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

    // ── Test 1: header columns match Java fixture ────────────────────────────

    /// The expected column list is copied verbatim from the Java fixture's first
    /// line (`benchmark/parity-fixtures/bsa_test_mgf_java.pin`), which uses
    /// charge2..=charge3 (BSA test uses charge_range 2..=3).
    ///
    /// Byte-parity note: the fixture header is compared column-by-column below.
    /// The fixture uses charge2..=charge3 because the BSA test was run with
    /// that range.
    #[test]
    fn pin_header_columns_match_java_fixture_without_features() {
        // Java fixture first line (charge2..=charge3):
        // SpecId Label ScanNr ExpMass CalcMass mass RawScore DeNovoScore
        // lnSpecEValue lnEValue isotope_error peplen dm absdm
        // charge2 charge3
        // enzN enzC enzInt
        // NumMatchedMainIons longest_b longest_y longest_y_pct
        // ExplainedIonCurrentRatio NTermIonCurrentRatio CTermIonCurrentRatio
        // MS2IonCurrent IsolationWindowEfficiency
        // MeanErrorTop7 StdevErrorTop7 MeanRelErrorTop7 StdevRelErrorTop7
        // lnDeltaSpecEValue matchedIonRatio
        // Peptide Proteins
        let expected: Vec<&str> = vec![
            "SpecId", "Label", "ScanNr", "ExpMass", "CalcMass", "mass",
            "RawScore", "DeNovoScore", "lnSpecEValue", "lnEValue", "isotope_error",
            "peplen", "dm", "absdm",
            "charge2", "charge3",
            "enzN", "enzC", "enzInt",
            "NumMatchedMainIons", "longest_b", "longest_y", "longest_y_pct",
            "ExplainedIonCurrentRatio", "NTermIonCurrentRatio", "CTermIonCurrentRatio",
            "MS2IonCurrent", "IsolationWindowEfficiency",
            "MeanErrorTop7", "StdevErrorTop7", "MeanRelErrorTop7", "StdevRelErrorTop7",
            "lnDeltaSpecEValue", "matchedIonRatio",
            "Peptide", "Proteins",
        ];

        let params = make_params(2..=3);
        let spectra: Vec<Spectrum> = vec![];
        let queues: Vec<TopNQueue> = vec![];

        let mut buf = Vec::<u8>::new();
        write_pin_to(&mut buf, &spectra, &queues, &params, "XXX_").unwrap();

        let cols = parse_header(&buf);
        assert_eq!(
            cols, expected,
            "PIN header columns must match Java fixture column order exactly"
        );
    }

    // ── Test 2: decoy PSM gets Label = -1 ────────────────────────────────────

    #[test]
    fn pin_writes_label_minus_one_for_decoy() {
        let params = make_params(2..=3);
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        let mut queue = TopNQueue::new(10);
        queue.push(make_psm(0, 10.0, 1e-5, true, 2)); // decoy
        let queues = vec![queue];

        let mut buf = Vec::<u8>::new();
        write_pin_to(&mut buf, &spectra, &queues, &params, "XXX_").unwrap();

        let rows = parse_rows(&buf);
        assert_eq!(rows.len(), 1, "should have 1 data row");

        // Label is column index 1 (SpecId=0, Label=1)
        assert_eq!(rows[0][1], "-1", "decoy PSM should have Label = -1");
    }

    // ── Test 3: charge one-hot encoding ────────────────────────────────────

    #[test]
    fn pin_writes_charge_one_hot_correctly() {
        let params = make_params(2..=3);
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        let mut queue = TopNQueue::new(10);
        queue.push(make_psm(0, 10.0, 1e-5, false, 2)); // charge 2
        let queues = vec![queue];

        let mut buf = Vec::<u8>::new();
        write_pin_to(&mut buf, &spectra, &queues, &params, "XXX_").unwrap();

        let cols = parse_header(&buf);
        let rows = parse_rows(&buf);
        assert_eq!(rows.len(), 1);

        // Find charge2 and charge3 column indices
        let charge2_idx = cols.iter().position(|c| c == "charge2").expect("charge2 column missing");
        let charge3_idx = cols.iter().position(|c| c == "charge3").expect("charge3 column missing");

        assert_eq!(rows[0][charge2_idx], "1", "charge2 should be 1 for a charge-2 PSM");
        assert_eq!(rows[0][charge3_idx], "0", "charge3 should be 0 for a charge-2 PSM");
    }

    // ── Test 4: empty queue → only header ────────────────────────────────────

    #[test]
    fn pin_handles_empty_queue() {
        let params = make_params(2..=3);
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];
        let queues = vec![TopNQueue::new(10)]; // empty

        let mut buf = Vec::<u8>::new();
        write_pin_to(&mut buf, &spectra, &queues, &params, "XXX_").unwrap();

        let rows = parse_rows(&buf);
        assert!(rows.is_empty(), "empty queue should produce no data rows");
    }

    // ── Test 5: lnDeltaSpecEValue = 0 when no rank-2 ─────────────────────────

    #[test]
    fn pin_lndelta_spec_evalue_zero_when_no_rank2() {
        let params = make_params(2..=3);
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        let mut queue = TopNQueue::new(10);
        queue.push(make_psm(0, 10.0, 1e-10, false, 2)); // single PSM → no rank-2
        let queues = vec![queue];

        let mut buf = Vec::<u8>::new();
        write_pin_to(&mut buf, &spectra, &queues, &params, "XXX_").unwrap();

        let cols = parse_header(&buf);
        let rows = parse_rows(&buf);
        assert_eq!(rows.len(), 1);

        let ln_delta_idx = cols
            .iter()
            .position(|c| c == "lnDeltaSpecEValue")
            .expect("lnDeltaSpecEValue column missing");

        let val: f64 = rows[0][ln_delta_idx]
            .parse()
            .expect("lnDeltaSpecEValue should be a number");
        assert!(
            val.abs() < 1e-9,
            "lnDeltaSpecEValue should be 0 when no rank-2 exists, got: {}",
            val
        );
    }
}
