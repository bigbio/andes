use std::borrow::Cow;
use std::io::{self, BufWriter, Write};

use crate::percolator::enzyme::enzyme_features;
use crate::rows::{emitted_rows, RowContext};
use model::mass::{ISOTOPE, PROTON};
use model::spectrum::Spectrum;
use rayon::prelude::*;
use search::candidates::Candidate;
use search::psm::{PsmMatch, TopNQueue};
use search::index::SearchIndex;
use search::params::SearchParams;

pub fn format_spec_id(
    spec_id: &str,
    scan: i32,
    rank: u32,
    row_idx: usize,
    multi_row: bool,
) -> String {
    if multi_row {
        format!("{spec_id}_{scan}_{rank}_{row_idx}")
    } else {
        format!("{spec_id}_{scan}_{rank}")
    }
}

pub fn spec_id_run(spec_id: &str) -> &str {
    spec_id.split_once('/').map_or("", |(run, _)| run)
}

/// The emitted rows of a spectrum in output order, each with its SpecId.
pub fn emitted_spec_ids<'a>(
    spec: &Spectrum,
    psms_sorted: &'a [PsmMatch],
    candidates: &[Candidate],
) -> Vec<(String, Cow<'a, PsmMatch>)> {
    let ranked = emitted_rows(spec, psms_sorted, candidates);
    let multi_row = ranked.len() > 1;
    ranked
        .into_iter()
        .enumerate()
        .map(|(row_idx, (rank, psm))| (emitted_spec_id(spec, rank, row_idx, multi_row), psm))
        .collect()
}

/// The SpecId of the emitted row `row_idx` (of rank `rank`) of a spectrum; `multi_row` when the
/// spectrum emits more than one row.
pub fn emitted_spec_id(spec: &Spectrum, rank: u32, row_idx: usize, multi_row: bool) -> String {
    let scan = spec.scan.unwrap_or(0);
    if spec.title.is_empty() {
        format_spec_id(&format!("scan={scan}"), scan, rank, row_idx, multi_row)
    } else {
        format_spec_id(&spec.title, scan, rank, row_idx, multi_row)
    }
}

pub fn emitted_row_evidence(
    spec: &Spectrum,
    queue: &TopNQueue,
    candidates: &[Candidate],
) -> Vec<(String, i8, Option<Box<search::sites::SiteEvidence>>)> {
    let psms = queue.clone().into_rank_sorted_vec();
    emitted_spec_ids(spec, &psms, candidates)
        .into_iter()
        .map(|(id, psm)| (id, psm.isotope_offset, psm.site_evidence.clone()))
        .collect()
}

pub fn emitted_row_candidates(
    spec: &Spectrum,
    queue: &TopNQueue,
    candidates: &[Candidate],
) -> Vec<(String, u32, u8)> {
    let psms = queue.clone().into_rank_sorted_vec();
    emitted_spec_ids(spec, &psms, candidates)
        .into_iter()
        .map(|(id, psm)| (id, psm.primary_candidate_idx(), psm.charge_used))
        .collect()
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum FeatureFmt {
    Int,
    Double,
    Fixed6,
}

pub type Feature = (&'static str, f64, FeatureFmt);

pub fn psm_feature_values(psm: &PsmMatch, rank: u32) -> Vec<Feature> {
    use FeatureFmt::{Double, Fixed6, Int};
    let f = &psm.features;
    let rank_score = psm.score.round() as i32 as f64;
    let delta_rank_score = if rank == 1 {
        f.delta_raw_score as f64
    } else {
        0.0
    };
    vec![
        ("RankScore", rank_score, Int),
        ("isotope_error", psm.isotope_offset as f64, Int),
        ("NumMatchedMainIons", f.num_matched_main_ions as f64, Int),
        ("longest_b", f.longest_b as f64, Int),
        ("longest_y", f.longest_y as f64, Int),
        ("longest_y_pct", f.longest_y_pct as f64, Fixed6),
        (
            "ExplainedIonCurrentRatio",
            f.explained_ion_current_ratio as f64,
            Double,
        ),
        (
            "NTermIonCurrentRatio",
            f.n_term_ion_current_ratio as f64,
            Double,
        ),
        (
            "CTermIonCurrentRatio",
            f.c_term_ion_current_ratio as f64,
            Double,
        ),
        ("MS2IonCurrent", f.ms2_ion_current as f64, Double),
        (
            "IsolationWindowEfficiency",
            f.isolation_window_efficiency as f64,
            Double,
        ),
        ("MeanErrorTop7", f.mean_error_top7 as f64, Double),
        ("StdevErrorTop7", f.stdev_error_top7 as f64, Double),
        ("MeanRelErrorTop7", f.mean_rel_error_top7 as f64, Double),
        ("StdevRelErrorTop7", f.stdev_rel_error_top7 as f64, Double),
        ("matchedIonRatio", f.matched_ion_ratio as f64, Double),
        ("EdgeScore", f.edge_score as f64, Int),
        ("PrecursorIsotopeKL", f.precursor_isotope_kl as f64, Double),
        ("PrecursorSNR", f.precursor_snr as f64, Double),
        ("DeltaRankScore", delta_rank_score, Double),
        ("TailorScore", f.tailor_score as f64, Double),
        ("PpmGaussianScore", f.ppm_gaussian_score as f64, Double),
        ("NeutralLossIonCount", f.neutral_loss_ion_count as f64, Int),
        (
            "LongestComplementaryLadder",
            f.longest_complementary_ladder as f64,
            Int,
        ),
        (
            "ComplementaryIonBalance",
            f.complementary_ion_balance as f64,
            Double,
        ),
        (
            "MeanMatchedIntensityRank",
            f.mean_matched_intensity_rank as f64,
            Double,
        ),
        (
            "DoublyChargedMatchedIonCount",
            f.doubly_charged_matched_ion_count as f64,
            Int,
        ),
        (
            "UniqueMatchFraction",
            f.unique_match_fraction as f64,
            Double,
        ),
        (
            "ChanceMatchSurprise",
            f.chance_match_surprise as f64,
            Double,
        ),
        ("IntensitySignal", f.intensity_signal as f64, Double),
        ("FragPredExplained", f.frag_pred_explained as f64, Double),
        ("FragPredChanceLLR", f.frag_pred_chance_llr as f64, Double),
        ("FragTopKObserved", f.frag_topk_observed as f64, Double),
        ("RichIonLLR", f.rich_ion_llr as f64, Double),
        ("IsRefinement", f.is_refinement as f64, Int),
        ("NumMods", f.num_mods as f64, Int),
        ("RefinementModClass", f.refine_mod_class as f64, Int),
        (
            "ModSiteShiftedMatched",
            f.mod_site_shifted_matched as f64,
            Double,
        ),
        ("ModSiteShiftedFrac", f.mod_site_shifted_frac as f64, Double),
        ("ModSiteIntensFrac", f.mod_site_intens_frac as f64, Double),
        ("ModSiteLocalized", f.mod_site_localized as f64, Double),
        ("ModSiteDetCount", f.mod_site_det_count as f64, Double),
        (
            "MassCompetitionEvidence",
            f.mass_competition_evidence as f64,
            Double,
        ),
        (
            "CandidateRankEntropy",
            f.candidate_rank_entropy as f64,
            Double,
        ),
        ("ListwiseScoreGap", f.listwise_score_gap as f64, Double),
        ("RawScore", f.strong_score as f64, Double),
        ("RawScoreCal", f.strong_score_cal as f64, Double),
        ("RankScoreFloat", f.rank_score_float as f64, Double),
        ("DeltaRT", f.delta_rt as f64, Double),
        ("AbsDeltaRT", f.abs_delta_rt as f64, Double),
        ("DeltaRTNorm", f.delta_rt_norm as f64, Double),
        ("IsVariant", f.is_variant as f64, Int),
        ("IsIsoform", f.is_isoform as f64, Int),
        ("IsChain", f.is_chain as f64, Int),
        ("RefCompetitor", f.ref_competitor as f64, Int),
        ("DeltaRefScore", f.delta_ref_score as f64, Double),
        ("IsotopeMimic", f.isotope_mimic as f64, Int),
        (
            "IsotopeErrorAbs",
            f64::from(psm.isotope_offset.unsigned_abs()),
            Int,
        ),
        (
            "IsotopeErrorNeg",
            f64::from(u8::from(psm.isotope_offset < 0)),
            Int,
        ),
    ]
}

pub fn write_pin(
    output_path: &std::path::Path,
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &[Candidate],
    params: &SearchParams,
    search_index: &SearchIndex,
) -> io::Result<()> {
    crate::atomic::write_atomic(output_path, |path| {
        let mut writer = BufWriter::new(std::fs::File::create(path)?);
        write_pin_to(
            &mut writer,
            spectra,
            queues,
            candidates,
            params,
            search_index,
        )?;
        writer.flush()
    })
}

pub fn write_pin_to<W: Write>(
    writer: &mut W,
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &[Candidate],
    params: &SearchParams,
    search_index: &SearchIndex,
) -> io::Result<()> {
    let (min_charge, max_charge) = pin_charge_columns(params);
    check_unique_spec_bases(spectra, queues)?;

    write_header(
        writer,
        min_charge,
        max_charge,
        params.open_search.is_some(),
    )?;

    let emitting: Vec<usize> = (0..queues.len())
        .filter(|&i| !queues[i].is_empty())
        .collect();
    let format_batch = |batch: &[usize]| -> Vec<io::Result<Vec<u8>>> {
        batch
            .par_iter()
            .map(|&spec_idx| {
                let mut buf = Vec::with_capacity(4096);
                write_spectrum_rows(
                    &mut buf,
                    &spectra[spec_idx],
                    &queues[spec_idx],
                    candidates,
                    min_charge,
                    max_charge,
                    search_index,
                    params,
                )?;
                Ok(buf)
            })
            .collect()
    };
    let format_batch = &format_batch;
    let mut batches = emitting.chunks(PIN_FORMAT_BATCH_SPECTRA);
    let mut ready = batches.next().map(format_batch);
    while let Some(current) = ready.take() {
        let next = batches.next();
        let mut formatted = None;
        let slot = &mut formatted;
        rayon::in_place_scope(|s| {
            if let Some(batch) = next {
                s.spawn(move |_| *slot = Some(format_batch(batch)));
            }
            current
                .into_iter()
                .try_for_each(|buf| writer.write_all(&buf?))
        })?;
        ready = formatted;
    }
    Ok(())
}

/// The SpecId base of a spectrum, `<title>_<scan>` (or `scan=<scan>_<scan>` without a title).
fn spec_base(spec: &Spectrum) -> String {
    let scan = spec.scan.unwrap_or(0);
    if spec.title.is_empty() {
        format!("scan={scan}_{scan}")
    } else {
        format!("{}_{scan}", spec.title)
    }
}

/// Fails when two emitting spectra share a SpecId base, naming the base of the first spectrum (in
/// spectrum order) whose base an earlier spectrum already has.
fn check_unique_spec_bases(spectra: &[Spectrum], queues: &[TopNQueue]) -> io::Result<()> {
    let mut bases: Vec<(String, usize)> = queues
        .par_iter()
        .enumerate()
        .filter(|(_, q)| !q.is_empty())
        .map(|(spec_idx, _)| (spec_base(&spectra[spec_idx]), spec_idx))
        .collect();
    bases.par_sort_unstable();
    let first_repeat = bases
        .par_windows(2)
        .filter(|w| w[0].0 == w[1].0)
        .map(|w| w[1].1)
        .min();
    match first_repeat {
        None => Ok(()),
        Some(spec_idx) => {
            let base = spec_base(&spectra[spec_idx]);
            Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "duplicate PIN SpecId base '{base}': two spectra share the same title and \
                     scan, which would collide on the Percolator q-value/PEP join and silently \
                     corrupt FDR. Ensure each spectrum has a unique title or scan number."
                ),
            ))
        }
    }
}

fn pin_charge_columns(params: &SearchParams) -> (u8, u8) {
    (*params.charge_range.start(), *params.charge_range.end())
}

const PIN_FORMAT_BATCH_SPECTRA: usize = 8192;

/// The PIN columns from `SpecId` to `absdm`, shared by the standard and the glyco PIN.
pub(crate) const PIN_LEAD_COLUMNS: [&str; 12] = [
    "SpecId",
    "Label",
    "ScanNr",
    "ExpMass",
    "CalcMass",
    "FileName",
    "mass",
    "RankScore",
    "isotope_error",
    "peplen",
    "dm",
    "absdm",
];

/// The PIN feature columns from `enzN` to `DeltaRTNorm`, shared by the standard and the glyco PIN.
pub(crate) const PIN_FEATURE_COLUMNS: [&str; 52] = [
    "enzN",
    "enzC",
    "enzInt",
    "NumMatchedMainIons",
    "longest_b",
    "longest_y",
    "longest_y_pct",
    "ExplainedIonCurrentRatio",
    "NTermIonCurrentRatio",
    "CTermIonCurrentRatio",
    "MS2IonCurrent",
    "IsolationWindowEfficiency",
    "MeanErrorTop7",
    "StdevErrorTop7",
    "MeanRelErrorTop7",
    "StdevRelErrorTop7",
    "matchedIonRatio",
    "EdgeScore",
    "PrecursorIsotopeKL",
    "PrecursorSNR",
    "DeltaRankScore",
    "TailorScore",
    "PpmGaussianScore",
    "NeutralLossIonCount",
    "LongestComplementaryLadder",
    "ComplementaryIonBalance",
    "MeanMatchedIntensityRank",
    "DoublyChargedMatchedIonCount",
    "UniqueMatchFraction",
    "ChanceMatchSurprise",
    "IntensitySignal",
    "FragPredExplained",
    "FragPredChanceLLR",
    "FragTopKObserved",
    "RichIonLLR",
    "IsRefinement",
    "NumMods",
    "RefinementModClass",
    "ModSiteShiftedMatched",
    "ModSiteShiftedFrac",
    "ModSiteIntensFrac",
    "ModSiteLocalized",
    "ModSiteDetCount",
    "MassCompetitionEvidence",
    "CandidateRankEntropy",
    "ListwiseScoreGap",
    "RawScore",
    "RawScoreCal",
    "RankScoreFloat",
    "DeltaRT",
    "AbsDeltaRT",
    "DeltaRTNorm",
];

fn write_header<W: Write>(
    writer: &mut W,
    min_charge: u8,
    max_charge: u8,
    open_search: bool,
) -> io::Result<()> {
    let mut cols: Vec<String> = PIN_LEAD_COLUMNS.iter().map(|c| c.to_string()).collect();
    for c in min_charge..=max_charge {
        cols.push(format!("charge{c}"));
    }
    cols.extend(PIN_FEATURE_COLUMNS.iter().map(|c| c.to_string()));
    cols.extend(
        [
            "IsVariant",
            "IsIsoform",
            "IsChain",
            "RefCompetitor",
            "DeltaRefScore",
            "IsotopeMimic",
            "IsotopeErrorAbs",
            "IsotopeErrorNeg",
        ]
        .map(String::from),
    );

    if open_search {
        cols.extend_from_slice(&[
            "DeltaMass".to_string(),
            "AbsDeltaMass".to_string(),
            "OpenSite".to_string(),
            "OpenLocScore".to_string(),
            "OpenLocMargin".to_string(),
        ]);
    }

    cols.extend_from_slice(&[
        "Peptide".to_string(),
        "Proteins".to_string(),
    ]);
    writeln!(writer, "{}", cols.join("\t"))
}

pub fn format_feature_value(value: f64, fmt: FeatureFmt) -> String {
    match fmt {
        FeatureFmt::Int => format!("{}", value.round() as i64),
        FeatureFmt::Fixed6 => format!("{value:.6}"),
        FeatureFmt::Double => {
            let mut buf = Vec::with_capacity(16);
            write_double(&mut buf, value).expect("write_double to Vec is infallible");
            String::from_utf8(buf).expect("write_double emits ASCII")
        }
    }
}

pub fn feature_value_type(fmt: FeatureFmt) -> &'static str {
    match fmt {
        FeatureFmt::Int => "int",
        FeatureFmt::Double | FeatureFmt::Fixed6 => "double",
    }
}

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
    let psms = queue.clone().into_rank_sorted_vec();

    let ranked = emitted_rows(spec, &psms, candidates);
    let multi_row = ranked.len() > 1;
    for (row_idx, (rank, psm)) in ranked.into_iter().enumerate() {
        let cand = &candidates[psm.primary_candidate_idx() as usize];
        let ctx = RowContext::new(spec, cand, search_index);
        write_psm_row(
            writer,
            spec,
            &psm,
            cand,
            &ctx,
            rank,
            row_idx,
            multi_row,
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
    row_idx: usize,
    multi_row: bool,
    min_charge: u8,
    max_charge: u8,
    candidates: &[Candidate],
    search_index: &SearchIndex,
    params: &SearchParams,
) -> io::Result<()> {
    let charge = psm.charge_used as f64;

    let label: i32 = if cand.is_decoy { -1 } else { 1 };

    let precursor_mz = psm.precursor_mz_override.unwrap_or(spec.precursor_mz);

    let exp_mass = precursor_mz * charge - charge * PROTON;

    let calc_mass = cand.peptide.mass();

    let mass = exp_mass;

    let raw_score = psm.score.round() as i32;

    let isotope_error: i32 = psm.isotope_offset as i32;

    let peplen = cand.peptide.length() + 2;

    let (dm, _) = precursor_error(precursor_mz, charge, isotope_error, calc_mass);
    let absdm = dm.abs();

    let matched_ion_ratio = psm.features.matched_ion_ratio as f64;

    write!(
        writer,
        "{}",
        format_spec_id(&ctx.spec_id, ctx.scan, rank, row_idx, multi_row)
    )?;
    write!(writer, "\t{}\t{}\t", label, ctx.scan)?;
    write_mass(writer, exp_mass)?;
    writer.write_all(b"\t")?;
    write_mass(writer, calc_mass)?;
    write!(writer, "\t{}\t", spec_id_run(&ctx.spec_id))?;
    write_double(writer, mass)?;
    write!(writer, "\t{raw_score}")?;
    write!(writer, "\t{isotope_error}\t{peplen}\t")?;
    write_double(writer, dm)?;
    writer.write_all(b"\t")?;
    write_double(writer, absdm)?;

    let hot_charge = psm.charge_used.clamp(min_charge, max_charge);
    for c in min_charge..=max_charge {
        let flag: u8 = if c == hot_charge { b'1' } else { b'0' };
        writer.write_all(&[b'\t', flag])?;
    }

    let (enz_n, enz_c, enz_int) = enzyme_features(cand, params.enzyme, &params.extra_enzymes);
    write!(writer, "\t{enz_n}\t{enz_c}\t{enz_int}")?;

    let _ = matched_ion_ratio;
    let features = psm_feature_values(psm, rank);
    debug_assert_eq!(features[0].0, "RankScore");
    debug_assert_eq!(features[1].0, "isotope_error");
    for &(_name, value, fmt) in &features[2..] {
        writer.write_all(b"\t")?;
        match fmt {
            FeatureFmt::Int => write!(writer, "{}", value.round() as i64)?,
            FeatureFmt::Fixed6 => write!(writer, "{value:.6}")?,
            FeatureFmt::Double => write_double(writer, value)?,
        }
    }

    if params.open_search.is_some() {
        let f = &psm.features;
        let delta = if f.open_delta != 0.0 {
            f64::from(f.open_delta)
        } else {
            exp_mass - f64::from(psm.isotope_offset) * ISOTOPE - calc_mass
        };
        writer.write_all(b"\t")?;
        write_mass(writer, delta)?;
        writer.write_all(b"\t")?;
        write_mass(writer, delta.abs())?;
        write!(writer, "\t{}", f.open_site)?;
        writer.write_all(b"\t")?;
        write_double(writer, f64::from(f.open_loc_score))?;
        writer.write_all(b"\t")?;
        write_double(writer, f64::from(f.open_loc_margin))?;
    }

    write!(writer, "\t{}", cand.peptide)?;
    for &cidx in &psm.candidate_idxs {
        let cand_for_acc = &candidates[cidx as usize];
        let accession = crate::rows::resolve_accession(cand_for_acc, search_index);
        write!(writer, "\t{accession}")?;
    }
    writeln!(writer)
}

pub(crate) fn precursor_error(
    precursor_mz: f64,
    charge: f64,
    isotope_error: i32,
    calc_mass: f64,
) -> (f64, f64) {
    let theo_mz = calc_mass / charge + PROTON;
    let adjusted_exp_mz = precursor_mz - ISOTOPE * (isotope_error as f64) / charge;
    (adjusted_exp_mz - theo_mz, theo_mz)
}

pub(crate) fn exp_mass_text(precursor_mz: f64, charge: u8) -> String {
    let z = charge as f64;
    let mut out = Vec::new();
    write_mass(&mut out, precursor_mz * z - z * PROTON).expect("write to a Vec");
    String::from_utf8(out).expect("ASCII mass")
}

fn write_mass<W: Write>(writer: &mut W, v: f64) -> io::Result<()> {
    if !v.is_finite() || v == 0.0 {
        return writer.write_all(b"0");
    }
    let mut buf = [0u8; 32];
    let len = {
        let mut cursor = &mut buf[..];
        write!(cursor, "{v:.6}")?;
        32 - cursor.len()
    };
    write_trim_fixed(writer, &buf[..len])
}

fn write_double<W: Write>(writer: &mut W, v: f64) -> io::Result<()> {
    if v.is_nan() || v.is_infinite() || v == 0.0 {
        return writer.write_all(b"0");
    }

    let mut buf = [0u8; 32];
    let abs = v.abs();
    if !(1e-4..1e6).contains(&abs) {
        let len = {
            let mut cursor = &mut buf[..];
            write!(cursor, "{v:.5e}")?;
            32 - cursor.len()
        };
        write_trim_scientific(writer, &buf[..len])
    } else {
        let digits_before_decimal = abs.log10().floor() as i32 + 1;
        let decimal_places = (6 - digits_before_decimal).max(0) as usize;
        let len = {
            let mut cursor = &mut buf[..];
            write!(cursor, "{v:.decimal_places$}")?;
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

    let exp_str = std::str::from_utf8(exp_part).unwrap_or("0");
    let exp_val: i32 = exp_str.parse().unwrap_or(0);
    write!(writer, "e{exp_val:+03}")
}

#[cfg(test)]
#[path = "pin/tests/pin.rs"]
mod tests;
