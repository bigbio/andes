//! PIN output writer.
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
//! # Column semantics
//!
//! * **Label**: source-protein TDC rule (iter27, 2026-05-21). `Label = -1`
//!   if the candidate's source protein is a decoy (`cand.is_decoy`), else
//!   `+1`. Matches Java MS-GF+ TDC labeling and avoids inflating Percolator's
//!   target set with peptides whose hit actually came from a decoy protein.
//!
//! * **isotope_error**: threaded from `PsmMatch::isotope_offset`, set by
//!   `match_engine.rs` from `MassError::isotope_offset`.
//!
//! * **enzN / enzC / enzInt**: computed via `crate::percolator_enz`,
//!   mirroring Java's `DirectPinWriter::isEnzymaticBoundary` +
//!   `countInternalEnzymatic` (OpenMS PercolatorInfile rules).
//!
//! * **Proteins**: single column with the real protein accession resolved from
//!   `SearchIndex::protein_at(candidates[psm.primary_candidate_idx() as usize].protein_index)`.
//!   Decoy accessions already carry the decoy prefix. Multi-protein support
//!   (merging Candidates that share pepSeq + score) comes in Task 4 of the R-2 refactor.
//!
//! * **peplen**: residue count + 2 (includes both flanking residues).
//!
//! * **dm / absdm**: mass error in Da using the matched isotope offset.
//!   `adjusted_exp_mz = precursor_mz - ISOTOPE * isotope_error / charge`
//!   (see `write_psm_row`), then `dm = adjusted_exp_mz - theo_mz` and
//!   `absdm = |dm|`. `isotope_error` is the PIN column from
//!   `PsmMatch::isotope_offset`.
//!
//! * **CalcMass**: `peptide.mass()` already includes H2O — neutral mass is
//!   computed directly from the peptide.
//!
//! ## Feature columns
//!
//! All 14 feature columns are filled from `psm.features` (computed by
//! `match_engine::compute_psm_features` at scoring time):
//! - `NumMatchedMainIons` — count of matched charge-1 b/y fragment positions.
//! - `longest_b` — longest contiguous run of matched b-ions.
//! - `longest_y` — longest contiguous run of matched y-ions.
//! - `longest_y_pct` — `longest_y / peptide.length()`.
//! - `ExplainedIonCurrentRatio` — matched b+y intensity / total MS2 intensity.
//! - `NTermIonCurrentRatio` — matched b intensity / total MS2 intensity.
//! - `CTermIonCurrentRatio` — matched y intensity / total MS2 intensity.
//! - `MS2IonCurrent` — raw sum of all MS2 peak intensities (NOT log10).
//! - `IsolationWindowEfficiency` — always 0.0 (not available from the Spectrum object).
//! - `MeanErrorTop7` — mean |Da| error of top-7 most-intense matched ions.
//! - `StdevErrorTop7` — population stdev of |Da| errors for top-7 ions.
//! - `MeanRelErrorTop7` — mean signed ppm error of top-7 ions.
//! - `StdevRelErrorTop7` — population stdev of signed ppm errors for top-7.
//! - `matchedIonRatio` — `NumMatchedMainIons / peptide.length()`.

use std::io::{self, BufWriter, Write};

use model::mass::{ISOTOPE, PROTON};
use crate::percolator_enz::{count_internal_enzymatic, is_enzymatic_boundary};
use crate::row_context::{iter_ranked, RowContext};
use search::candidate_gen::Candidate;
use search::psm::{PsmMatch, TopNQueue};
use search::search_index::SearchIndex;
use search::search_params::SearchParams;
use model::spectrum::Spectrum;

// ── public API ───────────────────────────────────────────────────────────────

/// Write all PSMs to a Percolator `.pin` file at `output_path`.
///
/// `spectra` and `queues` must be parallel slices (same length): `queues[i]`
/// holds the top-N PSMs for `spectra[i]`.
///
/// `candidates` is the per-search candidate pool owned by `PreparedSearch`.
/// PSM-to-candidate resolution goes through `candidates[psm.primary_candidate_idx() as usize]`.
///
/// `search_index` is used to resolve protein accessions from
/// `candidates[psm.primary_candidate_idx() as usize].protein_index`. The combined
/// target+decoy `ProteinDb` inside `search_index` already carries decoy
/// prefixes in the decoy accessions, so no separate prefix string is needed
/// for accession lookup. The `Label` column is derived directly from
/// `cand.is_decoy` (see `write_psm_row`).
pub fn write_pin(
    output_path: &std::path::Path,
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &[Candidate],
    params: &SearchParams,
    search_index: &SearchIndex,
) -> io::Result<()> {
    let file = std::fs::File::create(output_path)?;
    let mut writer = BufWriter::new(file);
    write_pin_to(&mut writer, spectra, queues, candidates, params, search_index)
}

/// Write all PSMs to an arbitrary writer — useful for testing without temp files.
///
/// See [`write_pin`] for parameter documentation.
pub fn write_pin_to<W: Write>(
    writer: &mut W,
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &[Candidate],
    params: &SearchParams,
    search_index: &SearchIndex,
) -> io::Result<()> {
    let min_charge = *params.charge_range.start();
    let max_charge = *params.charge_range.end();

    write_header(writer, min_charge, max_charge)?;

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
            min_charge,
            max_charge,
            search_index,
            params,
        )?;
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
        // Fragment-coverage + ion-current + error-stat features
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
        // ADDITIVE Java-parity feature (2026-05-21 iter19): per-bond
        // DBScanScorer edge sum (IES + error_score), emitted as a NEW
        // column so Percolator can learn weights without disrupting the
        // existing RawScore distribution.
        "EdgeScore".to_string(),
        // Peptide / Proteins
        "Peptide".to_string(),
        "Proteins".to_string(),
    ]);

    writeln!(writer, "{}", cols.join("\t"))
}

// ── per-spectrum rows ──────────────────────────────────────────────────────────

#[allow(clippy::too_many_arguments)]
fn write_spectrum_rows<W: Write>(
    writer: &mut W,
    spec: &Spectrum,
    queue: &TopNQueue,
    candidates: &[Candidate],
    min_charge: u8,
    max_charge: u8,
    search_index: &SearchIndex,
    params: &SearchParams,
) -> io::Result<()> {
    // Sort best-first (lowest spec_e_value first, then highest score).
    let psms = queue.clone().into_sorted_vec();

    // find rank-2 SpecEValue: first distinct spec_e_value after rank-1
    let rank2_spec_e_value = find_rank2_spec_e_value(&psms);

    for (rank, psm) in iter_ranked(&psms) {
        let cand = &candidates[psm.primary_candidate_idx() as usize];
        let ctx = RowContext::new(spec, cand, search_index);
        write_psm_row(
            writer,
            spec,
            psm,
            cand,
            &ctx,
            rank,
            rank2_spec_e_value,
            min_charge,
            max_charge,
            candidates,
            search_index,
            params,
        )?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn write_psm_row<W: Write>(
    writer: &mut W,
    spec: &Spectrum,
    psm: &PsmMatch,
    cand: &Candidate,
    ctx: &RowContext,
    rank: u32,
    rank2_spec_e_value: f64,
    min_charge: u8,
    max_charge: u8,
    candidates: &[Candidate],
    search_index: &SearchIndex,
    params: &SearchParams,
) -> io::Result<()> {
    let charge = psm.charge_used as f64;

    // iter27 (2026-05-21): label by SOURCE PROTEIN accession (standard TDC
    // convention, matches Java MS-GF+). Pre-iter27, Rust used an "any-target-
    // match" rule (Label = 1 if peptide sequence appears in ANY target
    // protein) which inflated target count when a peptide appeared in both
    // target and decoy proteins. Java labels by source: if the source
    // protein is a decoy, label = -1; otherwise +1.
    let label: i32 = if cand.is_decoy { -1 } else { 1 };

    // ExpMass: neutral precursor mass = mz * charge - charge * PROTON
    let exp_mass = spec.precursor_mz * charge - charge * PROTON;

    // CalcMass: theoretical neutral mass. peptide.mass() already includes H2O.
    // ExpMass = mz * charge - charge * PROTON is also a neutral mass.
    // Both columns must be neutral masses so that dm = ExpMass - CalcMass is a
    // true mass error (not a charge-induced offset). Fixture reference:
    // ExpMass=1641.96, CalcMass=1641.95 — both neutral.
    let calc_mass = cand.peptide.mass(); // includes H2O — neutral mass

    // mass: duplicate of ExpMass (column convention).
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

    // isotope_error: from PsmMatch::isotope_offset (threaded from
    // MassError::isotope_offset in match_engine.rs).
    let isotope_error: i32 = psm.isotope_offset as i32;

    // peplen: `residue_count + 2` (counts both flanking residues — the `pre`
    // and `post` characters in the `Peptide` struct). Without the +2, the
    // PIN row count and per-row diff disagree with the reference fixture.
    let peplen = cand.peptide.length() + 2;

    // dm / absdm: precursor mass error in Da.
    //   adjusted_exp_mz = precursor_mz - ISOTOPE * isotope_error / charge
    //   theo_mz         = peptide.mass() / charge + PROTON  (peptide.mass() includes H2O)
    //   dm              = adjusted_exp_mz - theo_mz
    let theo_mz = calc_mass / charge + PROTON;
    let adjusted_exp_mz = spec.precursor_mz - ISOTOPE * (isotope_error as f64) / charge;
    let dm = adjusted_exp_mz - theo_mz;
    let absdm = dm.abs();

    // lnDeltaSpecEValue
    let ln_delta_spec_e_value = compute_ln_delta_spec_e_value(rank, psm.spec_e_value, rank2_spec_e_value);

    // matchedIonRatio: from psm.features.
    let matched_ion_ratio = psm.features.matched_ion_ratio as f64;

    // Build row — tab-separated. We write directly into the BufWriter to
    // avoid heap-allocating each formatted column (the old implementation
    // built ~30 intermediate Strings per row × 37k rows = ~1.1M allocs).
    //
    // SpecId: `specID + "_" + scanNum + "_" + rank` — emitted inline via
    // three `write!` calls so we don't materialise a temporary String.
    write!(writer, "{}_{}_{}", ctx.spec_id, ctx.scan, rank)?;
    write!(writer, "\t{}\t{}\t", label, ctx.scan)?;
    write_double(writer, exp_mass)?;
    writer.write_all(b"\t")?;
    write_double(writer, calc_mass)?;
    writer.write_all(b"\t")?;
    write_double(writer, mass)?;
    write!(writer, "\t{}\t{}\t", raw_score, de_novo_score)?;
    write_double(writer, ln_spec_e_value)?;
    writer.write_all(b"\t")?;
    write_double(writer, ln_e_value)?;
    write!(writer, "\t{}\t{}\t", isotope_error, peplen)?;
    write_double(writer, dm)?;
    writer.write_all(b"\t")?;
    write_double(writer, absdm)?;

    // Charge one-hot
    for c in min_charge..=max_charge {
        let flag: u8 = if c == psm.charge_used { b'1' } else { b'0' };
        writer.write_all(&[b'\t', flag])?;
    }

    // enzN, enzC, enzInt — C-4 (2026-05-19): Java DirectPinWriter.java:199-203
    // emits enzymatic-boundary consistency features. enzN = boundary between
    // protein-pre and peptide[0]; enzC = boundary between peptide[last] and
    // protein-post; enzInt = count of internal positions consistent with the
    // enzyme. Per-rule semantics in crate::percolator_enz, mirroring Java's
    // isEnzymaticBoundary + countInternalEnzymatic (OpenMS PercolatorInfile).
    let residues: Vec<u8> = cand.peptide.residues.iter().map(|aa| aa.residue).collect();
    let first = residues.first().copied().unwrap_or(b'-');
    let last  = residues.last().copied().unwrap_or(b'-');
    let enz_n: u8 = is_enzymatic_boundary(cand.peptide.pre, first, params.enzyme) as u8;
    let enz_c: u8 = is_enzymatic_boundary(last, cand.peptide.post, params.enzyme) as u8;
    let enz_int = count_internal_enzymatic(&residues, params.enzyme);
    write!(writer, "\t{}\t{}\t{}", enz_n, enz_c, enz_int)?;

    // 4 fragment-coverage feature columns:
    // NumMatchedMainIons, longest_b, longest_y, longest_y_pct
    write!(
        writer,
        "\t{}\t{}\t{}\t{:.6}",
        psm.features.num_matched_main_ions,
        psm.features.longest_b,
        psm.features.longest_y,
        psm.features.longest_y_pct,
    )?;
    // 9 feature columns from psm.features:
    // ExplainedIonCurrentRatio, NTermIonCurrentRatio, CTermIonCurrentRatio,
    // MS2IonCurrent, IsolationWindowEfficiency,
    // MeanErrorTop7, StdevErrorTop7, MeanRelErrorTop7, StdevRelErrorTop7
    //
    // IsolationWindowEfficiency is always 0.0 (not available from the Spectrum object).
    writer.write_all(b"\t")?;
    write_double(writer, psm.features.explained_ion_current_ratio as f64)?;
    writer.write_all(b"\t")?;
    write_double(writer, psm.features.n_term_ion_current_ratio as f64)?;
    writer.write_all(b"\t")?;
    write_double(writer, psm.features.c_term_ion_current_ratio as f64)?;
    writer.write_all(b"\t")?;
    write_double(writer, psm.features.ms2_ion_current as f64)?;
    writer.write_all(b"\t")?;
    write_double(writer, psm.features.isolation_window_efficiency as f64)?;
    writer.write_all(b"\t")?;
    write_double(writer, psm.features.mean_error_top7 as f64)?;
    writer.write_all(b"\t")?;
    write_double(writer, psm.features.stdev_error_top7 as f64)?;
    writer.write_all(b"\t")?;
    write_double(writer, psm.features.mean_rel_error_top7 as f64)?;
    writer.write_all(b"\t")?;
    write_double(writer, psm.features.stdev_rel_error_top7 as f64)?;

    // lnDeltaSpecEValue, matchedIonRatio
    writer.write_all(b"\t")?;
    write_double(writer, ln_delta_spec_e_value)?;
    writer.write_all(b"\t")?;
    write_double(writer, matched_ion_ratio)?;

    // EdgeScore: additive Java-parity feature (iter19).
    writer.write_all(b"\t")?;
    write!(writer, "{}", psm.features.edge_score)?;

    // Peptide column (always one).
    // Proteins column(s): one tab-separated accession per candidate_idx.
    // After R-2.2 dedup, a PSM that matches the same peptide across multiple
    // proteins keeps all protein indices in candidate_idxs, and the PIN row
    // emits one accession per index — matching Java DirectPinWriter.java:237.
    // For PSMs with a single candidate_idx (typical), output is identical to
    // the pre-R-2.5 single-accession emit (ctx.accession still used by TSV).
    write!(writer, "\t{}", cand.peptide)?;
    for &cidx in &psm.candidate_idxs {
        let cand_for_acc = &candidates[cidx as usize];
        let accession = crate::row_context::resolve_accession(cand_for_acc, search_index);
        write!(writer, "\t{}", accession)?;
    }
    writeln!(writer)
}

// ── helpers ───────────────────────────────────────────────────────────────────

/// Find the rank-2 SpecEValue: the first distinct spec_e_value encountered after
/// the rank-1 value (skipping ties). Returns `f64::NAN` if no rank-2 exists.
///
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

/// Write a `f64` in `%.6g` style (6 significant figures) directly into
/// `writer`, matching Java's `String.format(Locale.ROOT, "%.6g", v)` used in
/// `formatDouble`.
///
/// NaN, infinite, or zero values are emitted as the single byte `'0'`
/// (matching Java's `if (Double.isNaN(v) || Double.isInfinite(v)) return "0";`).
///
/// This formats into a stack-allocated 32-byte buffer (sufficient for any
/// `%.5e`-style f64) and writes only the trimmed slice — avoiding the
/// per-call `String` allocation that the previous `format_double` returned.
fn write_double<W: Write>(writer: &mut W, v: f64) -> io::Result<()> {
    if v.is_nan() || v.is_infinite() || v == 0.0 {
        return writer.write_all(b"0");
    }

    // Stack buffer — 32 bytes is more than enough for any "%.5e" or
    // "%.prec$" formatting of an f64 (sign + 7 mantissa digits + 'e' +
    // signed 3-digit exponent ≈ 14 bytes worst case).
    let mut buf = [0u8; 32];
    let abs = v.abs();
    if !(1e-4..1e6).contains(&abs) {
        // Scientific notation, 5 decimal places after dot = 6 significant
        // digits. Format into stack buffer, then trim trailing zeros from
        // mantissa and reformat the exponent inline (no heap String).
        let len = {
            let mut cursor = &mut buf[..];
            write!(cursor, "{:.5e}", v)?;
            32 - cursor.len()
        };
        write_trim_scientific(writer, &buf[..len])
    } else {
        // Fixed notation. Determine decimal places for 6 sig figs.
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

/// Write the bytes in `s` to `writer`, trimming any trailing `'0'` (and a
/// dangling `'.'`) from a fixed-point representation. e.g. `"1.50000"` →
/// `"1.5"`. If `s` has no `'.'`, it is written verbatim.
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

/// Write a scientific-notation byte slice to `writer`, normalised to match
/// Java's `%g`-style output.
///
/// Rust formats `1.23456e7`; the reference fixture uses `1.23456e+07`. Trim trailing
/// zeros (and a dangling `.`) from the mantissa, then re-emit the exponent
/// with explicit sign and a minimum width of 2 digits (`e{:+03}` style).
fn write_trim_scientific<W: Write>(writer: &mut W, s: &[u8]) -> io::Result<()> {
    let pos = match s.iter().position(|&b| b == b'e') {
        Some(p) => p,
        None => return writer.write_all(s),
    };
    let mantissa = &s[..pos];
    let exp_part = &s[pos + 1..];

    // Trim trailing zeros (and a dangling '.') from the mantissa if it has
    // a decimal point.
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

    // Parse exponent and re-emit with explicit sign + min width 2. We
    // accept the same `unwrap_or(0)` semantics as the original code.
    let exp_str = std::str::from_utf8(exp_part).unwrap_or("0");
    let exp_val: i32 = exp_str.parse().unwrap_or(0);
    write!(writer, "e{:+03}", exp_val)
}


// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use model::amino_acid::AminoAcid;
    use search::candidate_gen::Candidate;
    use model::peptide::Peptide;
    use model::protein::{Protein, ProteinDb};
    use search::search_index::SearchIndex;
    use model::tolerance::PrecursorTolerance;
    use model::tolerance::Tolerance;

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

    /// Build an empty `SearchIndex` for tests that don't care about protein
    /// accessions (header / label / charge tests).
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

    /// Build a single Candidate for fixture tests. Mirrors the shape that the
    /// real candidate enumerator produces. Tests build a `Vec<Candidate>` from
    /// these and pass it to `write_pin_to`.
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

    fn make_psm(spectrum_idx: usize, score: f32, spec_e_value: f64, candidate_idx: u32, charge: u8) -> PsmMatch {
        PsmMatch {
            spectrum_idx,
            candidate_idxs: vec![candidate_idx],
            charge_used: charge,
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

    fn make_params(charge_range: std::ops::RangeInclusive<u8>) -> SearchParams {
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
            charge_range,
            isotope_error_range: -1..=2,
            top_n_psms_per_spectrum: 10,
            num_tolerable_termini: 2,
            min_peaks: 10,
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

    // ── Test 1: header columns match the reference fixture ──────────────────

    /// The expected column list is copied verbatim from the reference fixture's
    /// first line (`benchmark/parity-fixtures/bsa_test_mgf_java.pin`), which uses
    /// charge2..=charge3 (BSA test uses charge_range 2..=3).
    ///
    /// Byte-parity note: the fixture header is compared column-by-column below.
    #[test]
    fn pin_header_columns_match_java_fixture_without_features() {
        // Reference fixture first line (charge2..=charge3):
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
        // Java-fixture columns followed by Rust-only additive features.
        // `EdgeScore` is an iter19 ADDITIVE Java-parity feature emitted by
        // Rust only (Java doesn't compute it standalone — it's blended into
        // RawScore by DBScanScorer). Lives between matchedIonRatio and
        // Peptide so legacy Percolator readers using column order still
        // parse Peptide/Proteins at the tail.
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
            "EdgeScore",
            "Peptide", "Proteins",
        ];

        let params = make_params(2..=3);
        let spectra: Vec<Spectrum> = vec![];
        let queues: Vec<TopNQueue> = vec![];
        let idx = make_empty_search_index();

        let mut buf = Vec::<u8>::new();
        let cands: Vec<Candidate> = vec![];
        write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

        let cols = parse_header(&buf);
        assert_eq!(
            cols, expected,
            "PIN header columns must match the reference fixture column order exactly"
        );
    }

    // ── Test 2: decoy PSM gets Label = -1 ────────────────────────────────────

    #[test]
    fn pin_writes_label_minus_one_for_decoy() {
        let params = make_params(2..=3);
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        let mut queue = TopNQueue::new(10);
        queue.push(make_psm(0, 10.0, 1e-5, 0, 2)); // decoy
        let queues = vec![queue];
        let idx = make_empty_search_index();

        let mut buf = Vec::<u8>::new();
        let cands = vec![make_candidate(0, true)];
        write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

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
        queue.push(make_psm(0, 10.0, 1e-5, 0, 2)); // charge 2
        let queues = vec![queue];
        let idx = make_empty_search_index();

        let mut buf = Vec::<u8>::new();
        let cands = vec![make_candidate(0, false)];
        write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

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
        let idx = make_empty_search_index();

        let mut buf = Vec::<u8>::new();
        let cands: Vec<Candidate> = vec![];
        write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

        let rows = parse_rows(&buf);
        assert!(rows.is_empty(), "empty queue should produce no data rows");
    }

    // ── Test 5: lnDeltaSpecEValue = 0 when no rank-2 ─────────────────────────

    #[test]
    fn pin_lndelta_spec_evalue_zero_when_no_rank2() {
        let params = make_params(2..=3);
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        let mut queue = TopNQueue::new(10);
        queue.push(make_psm(0, 10.0, 1e-10, 0, 2)); // single PSM → no rank-2
        let queues = vec![queue];
        let idx = make_empty_search_index();

        let mut buf = Vec::<u8>::new();
        let cands = vec![make_candidate(0, false)];
        write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

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

    // ── Test 6: real accession emitted for target PSM ─────────────────────────

    #[test]
    fn pin_writes_real_accession_when_search_index_provided() {
        let accession = "sp|P02769|ALBU_BOVIN";
        let idx = make_search_index(accession);

        let params = make_params(2..=3);
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        // protein_index = 0 → first target protein
        let psm = make_psm(0, 10.0, 1e-5, 0, 2);

        let mut queue = TopNQueue::new(10);
        queue.push(psm);
        let queues = vec![queue];

        let mut buf = Vec::<u8>::new();
        let cands = vec![make_candidate(0, false)];
        write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

        let cols = parse_header(&buf);
        let rows = parse_rows(&buf);
        assert_eq!(rows.len(), 1);

        let prot_idx = cols.iter().position(|c| c == "Proteins").expect("Proteins column missing");
        assert_eq!(
            rows[0][prot_idx], accession,
            "Proteins column should contain the real accession, not a PROT_N placeholder"
        );
    }

    // ── Test 7: decoy accession carries decoy prefix ──────────────────────────

    #[test]
    fn pin_writes_decoy_prefix_for_decoy_protein() {
        let accession = "sp|P02769|ALBU_BOVIN";
        let idx = make_search_index(accession);

        let params = make_params(2..=3);
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        // SearchIndex has 1 target (idx 0) + 1 decoy (idx 1). Decoy accession
        // is set to "XXX_sp|P02769|ALBU_BOVIN" by target_plus_decoy.
        let psm = make_psm(0, 10.0, 1e-5, 0, 2);

        let mut queue = TopNQueue::new(10);
        queue.push(psm);
        let queues = vec![queue];

        let mut buf = Vec::<u8>::new();
        let cands = vec![make_candidate(1, true)]; // protein_index=1 (decoy slot), is_decoy=true
        write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

        let cols = parse_header(&buf);
        let rows = parse_rows(&buf);
        assert_eq!(rows.len(), 1);

        let prot_idx = cols.iter().position(|c| c == "Proteins").expect("Proteins column missing");
        let expected_decoy = format!("XXX_{}", accession);
        assert_eq!(
            rows[0][prot_idx], expected_decoy,
            "Proteins column should carry decoy prefix for decoy PSM"
        );
    }

    // ── Phase 7 followup: PIN emits real feature values ──────────────────────

    /// Verify that `NumMatchedMainIons` is emitted from `psm.features`
    /// rather than always being zero-stubbed.
    #[test]
    fn pin_emits_real_num_matched_main_ions_when_features_populated() {
        let params = make_params(2..=3);
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        let mut psm = make_psm(0, 10.0, 1e-5, 0, 2);
        psm.features.num_matched_main_ions = 5;

        let mut queue = TopNQueue::new(10);
        queue.push(psm);
        let queues = vec![queue];
        let idx = make_empty_search_index();

        let mut buf = Vec::<u8>::new();
        let cands = vec![make_candidate(0, false)];
        write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

        let cols = parse_header(&buf);
        let rows = parse_rows(&buf);
        assert_eq!(rows.len(), 1);

        let col_idx = cols
            .iter()
            .position(|c| c == "NumMatchedMainIons")
            .expect("NumMatchedMainIons column missing");
        assert_eq!(
            rows[0][col_idx], "5",
            "NumMatchedMainIons should be 5, not zero-stubbed"
        );
    }

    /// Verify that `longest_y_pct` is formatted with 6 decimal places.
    #[test]
    fn pin_emits_longest_y_pct_with_six_decimals() {
        let params = make_params(2..=3);
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        let mut psm = make_psm(0, 10.0, 1e-5, 0, 2);
        psm.features.longest_y = 1;
        psm.features.longest_y_pct = 0.5;

        let mut queue = TopNQueue::new(10);
        queue.push(psm);
        let queues = vec![queue];
        let idx = make_empty_search_index();

        let mut buf = Vec::<u8>::new();
        let cands = vec![make_candidate(0, false)];
        write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

        let cols = parse_header(&buf);
        let rows = parse_rows(&buf);
        assert_eq!(rows.len(), 1);

        let col_idx = cols
            .iter()
            .position(|c| c == "longest_y_pct")
            .expect("longest_y_pct column missing");
        assert_eq!(
            rows[0][col_idx], "0.500000",
            "longest_y_pct should be formatted with 6 decimal places"
        );
    }
}
