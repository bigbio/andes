//! PIN output writer.
//!
//! Produces a Percolator-consumable `.pin` file with the column layout used
//! by Kim et al. (Nat Commun 5:5277, 2014) and OpenMS PercolatorAdapter so
//! that downstream tools (Percolator, MS²Rescore, Mokapot) can consume the
//! output interchangeably.
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
//! * **Label**: source-protein TDC rule. `Label = -1`
//!   if the candidate's source protein is a decoy (`cand.is_decoy`), else
//!   `+1`. Standard target-decoy competition (TDC) labeling avoids inflating
//!   Percolator's target set with peptides whose hit came from a decoy protein.
//!
//! * **isotope_error**: threaded from `PsmMatch::isotope_offset`, set by
//!   `match_engine.rs` from `MassError::isotope_offset`.
//!
//! * **enzN / enzC / enzInt**: computed via `crate::percolator_enz`
//!   (OpenMS PercolatorInfile enzymatic-boundary rules; Kim et al., Nat Commun
//!   5:5277, 2014 Percolator features).
//!
//! * **Proteins**: single column with the real protein accession resolved from
//!   `SearchIndex::protein_at(candidates[psm.primary_candidate_idx() as usize].protein_index)`.
//!   Decoy accessions already carry the decoy prefix. Multi-protein support
//!   merges Candidates that share pepSeq + score.
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
//! - `MeanErrorTop7` — mean |ppm| error of top-7 most-intense matched ions.
//! - `StdevErrorTop7` — population stdev of |ppm| errors for top-7 ions.
//! - `MeanRelErrorTop7` — mean signed ppm error of top-7 ions.
//! - `StdevRelErrorTop7` — population stdev of signed ppm errors for top-7.
//! - `matchedIonRatio` — `NumMatchedMainIons / peptide.length()`.

use std::io::{self, BufWriter, Write};

use model::mass::{ISOTOPE, PROTON};
use crate::percolator_enz::{count_internal_enzymatic, is_enzymatic_boundary};
use crate::row_context::{iter_ranked_by_rank_score, RowContext};
use search::candidate_gen::Candidate;
use search::psm::{PsmMatch, TopNQueue};
use search::search_index::SearchIndex;
use search::search_params::SearchParams;
use model::spectrum::Spectrum;

// ── shared SpecId formatting ───────────────────────────────────────────────────

/// Format the Percolator `SpecId` (== QPX row identity) for one PSM row.
///
/// This is the SINGLE source of truth for the SpecId string, called from BOTH
/// the PIN writer (here) and the QPX writer (`crate::qpx`) so Percolator's
/// `PSMId` join key reconstructs identically on both sides.
///
/// Rule (must stay stable — it is the PIN/QPX/Percolator join contract):
/// - single-row scan → `"{spec_id}_{scan}_{rank}"`
/// - multi-row scan  → `"{spec_id}_{scan}_{rank}_{row_idx}"`
///
/// `spec_id` is `spec.title` (or `"scan={scan}"` when the title is empty — see
/// [`crate::row_context::RowContext`]); `rank` is the 1-based rank from
/// [`crate::row_context::iter_ranked_by_rank_score`]; `row_idx` is the per-scan
/// emission index; `multi_row` is true when the scan emits more than one row
/// (under `--chimeric` or when ties at queue capacity are retained), in which
/// case `_{row_idx}` disambiguates SpecIds that would otherwise collide on equal
/// `rank`. Single-row scans keep the historical `specID_scan_rank` format.
pub fn format_spec_id(spec_id: &str, scan: i32, rank: u32, row_idx: usize, multi_row: bool) -> String {
    if multi_row {
        format!("{spec_id}_{scan}_{rank}_{row_idx}")
    } else {
        format!("{spec_id}_{scan}_{rank}")
    }
}

// ── shared per-PSM feature vector ──────────────────────────────────────────────

/// How a feature value is rendered (PIN bytes + QPX `value_type`).
///
/// This is the SINGLE source of truth for the per-PSM discriminative feature
/// formatting, shared by the PIN writer (here) and the QPX `.idparquet` writer
/// (`crate::qpx`) so the two never diverge in which features they emit or how
/// they format them.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum FeatureFmt {
    /// Integer value. PIN: `{}`. QPX `value_type` = `"int"`.
    Int,
    /// General `%.6g`-style double (see [`write_double`]). PIN: `write_double`.
    /// QPX `value_type` = `"double"`.
    Double,
    /// Fixed 6-decimal double (`{:.6}`). Used only by `longest_y_pct` to match
    /// the historical PIN byte layout. QPX `value_type` = `"double"`.
    Fixed6,
}

/// One per-PSM discriminative feature: `(name, value, fmt)`.
pub type Feature = (&'static str, f64, FeatureFmt);

/// Compute the ordered per-PSM discriminative feature vector for one PSM.
///
/// This is the SINGLE source of truth for andes's Percolator feature set,
/// consumed by BOTH the PIN writer (`write_psm_row`) and the QPX `.idparquet`
/// writer (`crate::qpx::build_psms_batch`, into the `psm_metavalues` column) so
/// the discriminative power of the PIN and the idparquet never diverge.
///
/// The list intentionally EXCLUDES the columns OpenMS' PercolatorAdapter adds
/// generically (SpecId / Label / ScanNr / ExpMass / CalcMass / mass / peplen /
/// dm / absdm / the charge one-hots / enzN / enzC / enzInt) and the
/// non-numeric Peptide / Proteins columns. It DOES include `RankScore` and
/// `isotope_error` (the two numeric per-PSM signals the PIN writes in its
/// leading block) plus the entire fragment/error/strong-score feature block.
///
/// `rank` is the 1-based rank used to gate `DeltaRankScore` to the rank-1 row
/// (mirroring the PIN writer), so the value matches the PIN byte-for-byte.
pub fn psm_feature_values(psm: &PsmMatch, rank: u32) -> Vec<Feature> {
    use FeatureFmt::{Double, Fixed6, Int};
    let f = &psm.features;
    // RankScore: the integer-rounded ranking score (matches the PIN `RankScore`
    // column, which writes `psm.score.round() as i32`).
    let rank_score = psm.score.round() as i32 as f64;
    // DeltaRankScore is per-spectrum and emitted only on the rank-1 row.
    let delta_rank_score = if rank == 1 { f.delta_raw_score as f64 } else { 0.0 };
    vec![
        ("RankScore", rank_score, Int),
        ("isotope_error", psm.isotope_offset as f64, Int),
        ("NumMatchedMainIons", f.num_matched_main_ions as f64, Int),
        ("longest_b", f.longest_b as f64, Int),
        ("longest_y", f.longest_y as f64, Int),
        ("longest_y_pct", f.longest_y_pct as f64, Fixed6),
        ("ExplainedIonCurrentRatio", f.explained_ion_current_ratio as f64, Double),
        ("NTermIonCurrentRatio", f.n_term_ion_current_ratio as f64, Double),
        ("CTermIonCurrentRatio", f.c_term_ion_current_ratio as f64, Double),
        ("MS2IonCurrent", f.ms2_ion_current as f64, Double),
        ("IsolationWindowEfficiency", f.isolation_window_efficiency as f64, Double),
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
        ("LongestComplementaryLadder", f.longest_complementary_ladder as f64, Int),
        ("ComplementaryIonBalance", f.complementary_ion_balance as f64, Double),
        ("MeanMatchedIntensityRank", f.mean_matched_intensity_rank as f64, Double),
        ("DoublyChargedMatchedIonCount", f.doubly_charged_matched_ion_count as f64, Int),
        ("UniqueMatchFraction", f.unique_match_fraction as f64, Double),
        ("ChanceMatchSurprise", f.chance_match_surprise as f64, Double),
        ("IntensitySignal", f.intensity_signal as f64, Double),
        ("FragPredExplained", f.frag_pred_explained as f64, Double),
        ("FragPredChanceLLR", f.frag_pred_chance_llr as f64, Double),
        ("FragTopKObserved", f.frag_topk_observed as f64, Double),
        ("RichIonLLR", f.rich_ion_llr as f64, Double),
        ("IsRefinement", f.is_refinement as f64, Int),
        ("NumMods", f.num_mods as f64, Int),
        ("RefinementModClass", f.refine_mod_class as f64, Int),
        ("ModSiteShiftedMatched", f.mod_site_shifted_matched as f64, Double),
        ("ModSiteShiftedFrac", f.mod_site_shifted_frac as f64, Double),
        ("ModSiteIntensFrac", f.mod_site_intens_frac as f64, Double),
        ("ModSiteLocalized", f.mod_site_localized as f64, Double),
        ("ModSiteDetCount", f.mod_site_det_count as f64, Double),
        ("MassCompetitionEvidence", f.mass_competition_evidence as f64, Double),
        ("CandidateRankEntropy", f.candidate_rank_entropy as f64, Double),
        ("ListwiseScoreGap", f.listwise_score_gap as f64, Double),
        ("RawScore", f.strong_score as f64, Double),
        ("RawScoreCal", f.strong_score_cal as f64, Double),
    ]
}

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

fn write_header<W: Write>(
    writer: &mut W,
    min_charge: u8,
    max_charge: u8,
) -> io::Result<()> {
    // RawScore is the sole score column. The generating function has been
    // removed, so the GF-derived columns (DeNovoScore / lnSpecEValue / lnEValue
    // / lnDeltaSpecEValue) are not emitted: Percolator calibrates FDR from
    // RawScore + the remaining fragment/mass features.
    let mut cols: Vec<String> = vec![
        "SpecId".to_string(),
        "Label".to_string(),
        "ScanNr".to_string(),
        "ExpMass".to_string(),
        "CalcMass".to_string(),
        "mass".to_string(),
        // RankScore = the rank-LLR ranking score (formerly "RawScore"; value is
        // `psm.score`). Renamed so the headline/primary column is "RawScore"
        // (the fused score, formerly "StrongScore").
        "RankScore".to_string(),
    ];
    cols.extend_from_slice(&[
        "isotope_error".to_string(),
        "peplen".to_string(),
        "dm".to_string(),
        "absdm".to_string(),
    ]);

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
    ]);
    cols.extend_from_slice(&[
        "matchedIonRatio".to_string(),
        // ADDITIVE feature: per-bond edge sum (IES + error_score), emitted as a NEW
        // column so Percolator can learn weights without disrupting the
        // existing RawScore distribution (Kim et al., Nat Commun 5:5277, 2014).
        "EdgeScore".to_string(),
        // ADDITIVE chimeric MS1 precursor-envelope features: emitted
        // adjacent to EdgeScore so they sit just before Peptide/Proteins.
        // Both are 0.0 unless `--chimeric` populates them from a linked MS1.
        "PrecursorIsotopeKL".to_string(),
        "PrecursorSNR".to_string(),
        // ADDITIVE top-1 dominance feature: RawScore(best) − RawScore(2nd-best
        // distinct peptide) on the rank-1 row, 0.0 elsewhere. Built on
        // parity-grade RawScore (not the divergent SpecE), so it adds an
        // orthogonal "lead over the runner-up" signal without touching any
        // existing column. Populated only when a distinct runner-up was scored
        // (i.e. effectively needs internal retention ≥ 2 candidates per scan).
        "DeltaRankScore".to_string(),
        // ADDITIVE Tailor per-spectrum calibration (Yang et al., JPR 2020):
        // RawScore / (spectrum's top-1% quantile RawScore). Makes RawScores
        // comparable across spectra — the role the removed generating function
        // used to play — recovering low-res discrimination without the GF.
        "TailorScore".to_string(),
        // ADDITIVE strong-score Stage-1 bolt-ons (deterministic, no model
        // change): PpmGaussianScore = Σ exp(-½(ppm/7)²) over matched ions
        // (turns fragment mass accuracy into evidence the rank model discards);
        "PpmGaussianScore".to_string(),
        // NeutralLossIonCount = matched b/y ions with −H2O/−NH3 partner peaks.
        "NeutralLossIonCount".to_string(),
        // LongestComplementaryLadder = longest consecutive run of complementary
        // cleavage sites (both b and y matched).
        "LongestComplementaryLadder".to_string(),
        // ComplementaryIonBalance = Σ over bonds where both b_i and y_{n-i}
        // matched, weighted by intensity-rank agreement 1/(1+|rank_b−rank_y|)
        // (ADDITIVE; orthogonal to the run-length ladder above).
        "ComplementaryIonBalance".to_string(),
        // MeanMatchedIntensityRank = mean intensity-rank of matched ions (lower
        // = matched dominant peaks).
        "MeanMatchedIntensityRank".to_string(),
        // DoublyChargedMatchedIonCount = matched charge-2 b/y ions.
        "DoublyChargedMatchedIonCount".to_string(),
        // UniqueMatchFraction = within-peptide peak-explanation uniqueness.
        "UniqueMatchFraction".to_string(),
        // ChanceMatchSurprise = strong-score Stage-2 null moat: Σ max(0,
        // -ln(ρ·Δ)) per matched ion — how improbable the matches are by chance.
        "ChanceMatchSurprise".to_string(),
        // IntensitySignal = strong-score S1 numerator: cosine similarity between
        // IntensityModel predictions and observed relative intensities (0 without model).
        "IntensitySignal".to_string(),
        // Tier-2 frag-intensity LLR battery (0 without a frag-intensity model).
        // FragPredExplained = Σ(matched·pred)/Σpred; FragPredChanceLLR =
        // Σ matched·pred·max(0,−ln p_chance); FragTopKObserved = top-K hit rate.
        "FragPredExplained".to_string(),
        "FragPredChanceLLR".to_string(),
        "FragTopKObserved".to_string(),
        // RichIonLLR = decoy-aware per-annotated-ion LLR sum (0 without a rich-ion model).
        "RichIonLLR".to_string(),
        // Refinement-cascade additive columns (0 without --refine):
        // IsRefinement = from Pass-2 search; NumMods = variable-mod count;
        // RefinementModClass = mod-class id for subgroup-FDR grouping.
        "IsRefinement".to_string(),
        "NumMods".to_string(),
        "RefinementModClass".to_string(),
        // Mod-localization site-determining-ion columns (0 for unmodified
        // peptides). ModSiteShiftedMatched = matched mod-bearing b/y ions;
        // ModSiteShiftedFrac = matched/total shifted; ModSiteIntensFrac =
        // shifted/all matched intensity; ModSiteLocalized = 1 if a bracketing
        // ion pair localizes the mod; ModSiteDetCount = # site-determining ions.
        "ModSiteShiftedMatched".to_string(),
        "ModSiteShiftedFrac".to_string(),
        "ModSiteIntensFrac".to_string(),
        "ModSiteLocalized".to_string(),
        "ModSiteDetCount".to_string(),
        // MassCompetitionEvidence = S2 null term 2: Σ 1/(1+ambiguity+ρ).
        "MassCompetitionEvidence".to_string(),
        // CandidateRankEntropy = S2 listwise: softmax entropy over retained top-K.
        "CandidateRankEntropy".to_string(),
        // ListwiseScoreGap = S2 listwise: top-1 − top-2 RawScore in retained queue.
        "ListwiseScoreGap".to_string(),
        // RawScore = S3 fused signal − null (formerly "StrongScore"; the
        // headline/primary score, always emitted; ranks when --score strong).
        "RawScore".to_string(),
        // RawScoreCal = S4 per-spectrum z-scored significance (formerly "StrongScoreCal").
        "RawScoreCal".to_string(),
    ]);

    cols.extend_from_slice(&[
        // Peptide / Proteins
        "Peptide".to_string(),
        "Proteins".to_string(),
    ]);

    writeln!(writer, "{}", cols.join("\t"))
}

/// Format a feature value to a string using the SAME rules the PIN writer uses,
/// so the QPX `psm_metavalues` strings and the PIN columns are byte-consistent.
///
/// - [`FeatureFmt::Int`]: integer (`{}` on the rounded value).
/// - [`FeatureFmt::Double`]: `%.6g`-style via [`write_double`].
/// - [`FeatureFmt::Fixed6`]: `{:.6}` fixed point.
pub fn format_feature_value(value: f64, fmt: FeatureFmt) -> String {
    match fmt {
        FeatureFmt::Int => format!("{}", value.round() as i64),
        FeatureFmt::Fixed6 => format!("{:.6}", value),
        FeatureFmt::Double => {
            let mut buf = Vec::with_capacity(16);
            // write_double is infallible for an in-memory Vec.
            write_double(&mut buf, value).expect("write_double to Vec is infallible");
            String::from_utf8(buf).expect("write_double emits ASCII")
        }
    }
}

/// QPX `value_type` string for a feature format (mirrors a comparison engine's convention:
/// `"int"` for integers, `"double"` for floats).
pub fn feature_value_type(fmt: FeatureFmt) -> &'static str {
    match fmt {
        FeatureFmt::Int => "int",
        FeatureFmt::Double | FeatureFmt::Fixed6 => "double",
    }
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
    // Order by rank_score (RawScore) descending — the sole ranking signal.
    let psms = queue.clone().into_rank_sorted_vec();

    let ranked: Vec<(u32, &PsmMatch)> = iter_ranked_by_rank_score(&psms).collect();
    // When a scan emits more than one row, `rank` alone can collide (ranks tie
    // on equal `rank_score`, and `TopNQueue` retains ties at capacity), so the
    // SpecId must include the per-row index to stay unique. Single-row scans
    // keep the historical `specID_scan_rank` format (schema parity).
    let multi_row = ranked.len() > 1;
    for (row_idx, (rank, psm)) in ranked.into_iter().enumerate() {
        let cand = &candidates[psm.primary_candidate_idx() as usize];
        let ctx = RowContext::new(spec, cand, search_index);
        write_psm_row(
            writer,
            spec,
            psm,
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

    // Label by SOURCE PROTEIN accession (standard TDC convention).
    // An "any-target-match" rule (Label = 1 if peptide sequence appears in
    // ANY target protein) would inflate target count when a peptide appeared
    // in both target and decoy proteins. Label by source: decoy protein → -1,
    // otherwise +1.
    let label: i32 = if cand.is_decoy { -1 } else { 1 };

    // For chimeric Pass-2 secondaries, mass-error columns use the co-isolated
    // precursor m/z (the secondary's true precursor). `None` (every ordinary PSM)
    // falls back to the spectrum's precursor m/z, keeping that path byte-identical.
    let precursor_mz = psm.precursor_mz_override.unwrap_or(spec.precursor_mz);

    // ExpMass: neutral precursor mass = mz * charge - charge * PROTON
    let exp_mass = precursor_mz * charge - charge * PROTON;

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
    let adjusted_exp_mz = precursor_mz - ISOTOPE * (isotope_error as f64) / charge;
    let dm = adjusted_exp_mz - theo_mz;
    let absdm = dm.abs();

    // matchedIonRatio: from psm.features.
    let matched_ion_ratio = psm.features.matched_ion_ratio as f64;

    // Write columns directly into the BufWriter (avoids ~30 String allocs/row).
    //
    // SpecId = `specID_scanNum_rank`. Whenever a scan emits more than one row
    // (under --chimeric, OR because ranks tie on equal `rank_score` and the
    // `TopNQueue` retained the ties), `_{rank}` can collide, producing duplicate
    // SpecIds in the PIN (ambiguous downstream mapping). Append the per-row
    // emission index to disambiguate. Single-row-per-scan keeps the historical
    // `specID_scan_rank` format so the schema/common case is unchanged.
    write!(writer, "{}", format_spec_id(&ctx.spec_id, ctx.scan, rank, row_idx, multi_row))?;
    write!(writer, "\t{}\t{}\t", label, ctx.scan)?;
    write_double(writer, exp_mass)?;
    writer.write_all(b"\t")?;
    write_double(writer, calc_mass)?;
    writer.write_all(b"\t")?;
    write_double(writer, mass)?;
    // RawScore is the sole score column (GF-derived DeNovoScore / lnSpecEValue /
    // lnEValue are no longer emitted).
    write!(writer, "\t{}", raw_score)?;
    write!(writer, "\t{}\t{}\t", isotope_error, peplen)?;
    write_double(writer, dm)?;
    writer.write_all(b"\t")?;
    write_double(writer, absdm)?;

    // Charge one-hot
    for c in min_charge..=max_charge {
        let flag: u8 = if c == psm.charge_used { b'1' } else { b'0' };
        writer.write_all(&[b'\t', flag])?;
    }

    // enzN, enzC, enzInt — Percolator enzymatic-boundary features.
    // enzN = boundary between protein-pre and peptide[0]; enzC = boundary
    // between peptide[last] and protein-post; enzInt = count of internal
    // positions consistent with the enzyme. Per-rule semantics in
    // crate::percolator_enz (OpenMS PercolatorInfile).
    let residues: Vec<u8> = cand.peptide.residues.iter().map(|aa| aa.residue).collect();
    let first = residues.first().copied().unwrap_or(b'-');
    let last  = residues.last().copied().unwrap_or(b'-');
    let enz_n: u8 = is_enzymatic_boundary(cand.peptide.pre, first, params.enzyme) as u8;
    let enz_c: u8 = is_enzymatic_boundary(last, cand.peptide.post, params.enzyme) as u8;
    let enz_int = count_internal_enzymatic(&residues, params.enzyme);
    write!(writer, "\t{}\t{}\t{}", enz_n, enz_c, enz_int)?;

    // Fragment/error/strong-score feature block — written from the SHARED
    // per-PSM feature vector (`psm_feature_values`) so the PIN and the QPX
    // `.idparquet` (`crate::qpx`, `psm_metavalues` column) can never diverge.
    //
    // `psm_feature_values` returns the full ordered feature list; its first two
    // entries (`RankScore`, `isotope_error`) are written inline above in their
    // historical leading-block positions, so this loop skips them and emits the
    // tail (`NumMatchedMainIons` … `RawScoreCal`) — the exact PIN column order.
    //
    // `matched_ion_ratio` is bound above for the doc cross-reference; the shared
    // list re-derives the identical value, so the binding is intentionally inert
    // here (silence the unused-variable lint without changing byte output).
    let _ = matched_ion_ratio;
    let features = psm_feature_values(psm, rank);
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

    // Peptide column (always one).
    // Proteins column(s): one tab-separated accession per candidate_idx.
    // After pepSeq+score dedup, a PSM that matches the same peptide across
    // multiple proteins keeps all protein indices in candidate_idxs, and the
    // PIN row emits one accession per index for multi-protein shared peptides.
    write!(writer, "\t{}", cand.peptide)?;
    for &cidx in &psm.candidate_idxs {
        let cand_for_acc = &candidates[cidx as usize];
        let accession = crate::row_context::resolve_accession(cand_for_acc, search_index);
        write!(writer, "\t{}", accession)?;
    }
    writeln!(writer)
}

// ── helpers ───────────────────────────────────────────────────────────────────

/// Write a `f64` in `%.6g` style (6 significant figures) directly into
/// `writer`.
///
/// NaN, infinite, or zero values are emitted as the single byte `'0'`.
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

/// Write a scientific-notation byte slice to `writer`, normalised to `%g`-style
/// output with explicit signed exponent (`e{:+03}` style).
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
            isolation_lower_offset: None,
            isolation_upper_offset: None,
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

    fn make_psm(spectrum_idx: usize, score: f32, rank_score: f32, candidate_idx: u32, charge: u8) -> PsmMatch {
        PsmMatch {
            spectrum_idx,
            candidate_idxs: vec![candidate_idx],
            charge_used: charge,
            mass_error_ppm: 1.5,
            score,
            rank_score,
            edge_score: 0,
            activation_method: Some(model::activation::ActivationMethod::HCD),
            features: search::psm::PsmFeatures::default(),
            isotope_offset: 0,
            precursor_mz_override: None,
        }
    }

    fn make_params(charge_range: std::ops::RangeInclusive<u8>) -> SearchParams {
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
            precursor_tolerance: PrecursorTolerance::symmetric(Tolerance::Ppm(20.0)),
            charge_range,
            isotope_error_range: -1..=2,
            top_n_psms_per_spectrum: 10,
            num_tolerable_termini: 2,
            min_peaks: 10,
            precursor_cal_mode: search::PrecursorCalMode::Auto,
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

    /// The expected column list is andes's own PIN schema (defined by
    /// `write_header`), here pinned for charge2..=charge3 (the BSA test uses
    /// charge_range 2..=3).
    ///
    /// The header is asserted column-by-column below.
    #[test]
    fn pin_header_columns_are_gf_free_schema() {
        // GF-free schema: RawScore is the sole score column; the GF-derived
        // DeNovoScore / lnSpecEValue / lnEValue / lnDeltaSpecEValue columns are
        // NOT emitted. The additive feature columns (EdgeScore,
        // PrecursorIsotopeKL, PrecursorSNR, DeltaRawScore) sit between
        // matchedIonRatio and Peptide.
        let expected: Vec<&str> = vec![
            "SpecId", "Label", "ScanNr", "ExpMass", "CalcMass", "mass",
            "RankScore", "isotope_error",
            "peplen", "dm", "absdm",
            "charge2", "charge3",
            "enzN", "enzC", "enzInt",
            "NumMatchedMainIons", "longest_b", "longest_y", "longest_y_pct",
            "ExplainedIonCurrentRatio", "NTermIonCurrentRatio", "CTermIonCurrentRatio",
            "MS2IonCurrent", "IsolationWindowEfficiency",
            "MeanErrorTop7", "StdevErrorTop7", "MeanRelErrorTop7", "StdevRelErrorTop7",
            "matchedIonRatio",
            "EdgeScore",
            "PrecursorIsotopeKL", "PrecursorSNR", "DeltaRankScore", "TailorScore",
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

    /// Additive refinement-cascade columns are always present in the header
    /// (they carry 0 without `--refine`). Guards Task 2 of the PTM cascade.
    #[test]
    fn pin_header_has_refinement_columns() {
        let params = make_params(2..=3);
        let spectra: Vec<Spectrum> = vec![];
        let queues: Vec<TopNQueue> = vec![];
        let idx = make_empty_search_index();

        let mut buf = Vec::<u8>::new();
        let cands: Vec<Candidate> = vec![];
        write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

        let cols = parse_header(&buf);
        for col in ["IsRefinement", "NumMods", "RefinementModClass"] {
            assert!(cols.iter().any(|c| c == col), "header missing column {col}");
        }
    }

    // ── Test 2: decoy PSM gets Label = -1 ────────────────────────────────────

    #[test]
    fn pin_writes_label_minus_one_for_decoy() {
        let params = make_params(2..=3);
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        let mut queue = TopNQueue::new(10);
        queue.push(make_psm(0, 10.0, 10.0, 0, 2)); // decoy
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
        queue.push(make_psm(0, 10.0, 10.0, 0, 2)); // charge 2
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

    // ── Test: chimeric SpecId uniqueness for co-fragmented peptides ─────────

    #[test]
    fn chimeric_specids_unique_for_cofragmented_peptides_same_scan() {
        let mut params = make_params(2..=3);
        params.chimeric = true;
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        // Two distinct peptides (candidates 0 and 1) with the SAME rank_score:
        // iter_ranked_by_rank_score assigns them the same rank, so without the
        // chimeric per-row suffix their SpecIds (`spec_scan_rank`) would collide.
        let mut queue = TopNQueue::new(10);
        queue.push(make_psm(0, 10.0, 10.0, 0, 2));
        queue.push(make_psm(0, 9.0, 10.0, 1, 2));
        let queues = vec![queue];
        let idx = make_empty_search_index();
        let cands = vec![make_candidate(0, false), make_candidate(1, false)];

        let mut buf = Vec::<u8>::new();
        write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

        let rows = parse_rows(&buf);
        assert_eq!(rows.len(), 2, "both co-fragmented PSMs should be emitted");
        assert_ne!(rows[0][0], rows[1][0],
            "chimeric SpecIds must be unique per row, got {:?} and {:?}", rows[0][0], rows[1][0]);
    }

    // ── Test: non-chimeric tied-PSM SpecId uniqueness ──────────────────────
    #[test]
    fn non_chimeric_tied_psms_same_scan_get_distinct_specids() {
        // With the generating function removed, ranks tie on equal `rank_score`
        // and `TopNQueue` retains the ties at capacity, so even a *non-chimeric*
        // scan can emit two rows that would share `spec_scan_rank`. Both rows
        // must still get distinct SpecIds.
        let params = make_params(2..=3); // chimeric == false (default)
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        // Two distinct peptides with the SAME rank_score (10.0) → same rank.
        let mut queue = TopNQueue::new(10);
        queue.push(make_psm(0, 10.0, 10.0, 0, 2));
        queue.push(make_psm(0, 10.0, 10.0, 1, 2));
        let queues = vec![queue];
        let idx = make_empty_search_index();
        let cands = vec![make_candidate(0, false), make_candidate(1, false)];

        let mut buf = Vec::<u8>::new();
        write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

        let rows = parse_rows(&buf);
        assert_eq!(rows.len(), 2, "both tied PSMs should be emitted");
        assert_ne!(
            rows[0][0], rows[1][0],
            "non-chimeric tied SpecIds must be unique, got {:?} and {:?}",
            rows[0][0], rows[1][0]
        );
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

    // ── Test 6: real accession emitted for target PSM ─────────────────────────

    #[test]
    fn pin_writes_real_accession_when_search_index_provided() {
        let accession = "sp|P02769|ALBU_BOVIN";
        let idx = make_search_index(accession);

        let params = make_params(2..=3);
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        // protein_index = 0 → first target protein
        let psm = make_psm(0, 10.0, 10.0, 0, 2);

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
        let psm = make_psm(0, 10.0, 10.0, 0, 2);

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

    // ── PIN emits real feature values ────────────────────────────────────────

    /// Verify that `NumMatchedMainIons` is emitted from `psm.features`.
    #[test]
    fn pin_emits_real_num_matched_main_ions_when_features_populated() {
        let params = make_params(2..=3);
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        let mut psm = make_psm(0, 10.0, 10.0, 0, 2);
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

    /// DeltaRawScore is emitted on the rank-1 row from `features.delta_raw_score`
    /// (a per-spectrum scalar the match engine stores on every retained PSM),
    /// and gated to 0.0 on rank > 1 rows — mirroring lnDeltaSpecEValue.
    #[test]
    fn pin_emits_delta_raw_score_on_rank1_only() {
        let params = make_params(2..=3);
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        // Two distinct-rank_score PSMs on one spectrum → rank 1 then rank 2. The
        // engine stores the same spectrum-level delta on both; the writer must
        // emit it for rank 1 and 0.0 for rank 2 (no double-attribution).
        let mut psm1 = make_psm(0, 12.0, 12.0, 0, 2);
        psm1.features.delta_raw_score = 7.0;
        let mut psm2 = make_psm(0, 5.0, 5.0, 1, 2);
        psm2.features.delta_raw_score = 7.0;

        let mut queue = TopNQueue::new(10);
        queue.push(psm1);
        queue.push(psm2);
        let queues = vec![queue];
        let idx = make_empty_search_index();

        let mut buf = Vec::<u8>::new();
        let cands = vec![make_candidate(0, false), make_candidate(1, false)];
        write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

        let cols = parse_header(&buf);
        let rows = parse_rows(&buf);
        assert_eq!(rows.len(), 2, "two PSMs → two rows");

        let col_idx = cols
            .iter()
            .position(|c| c == "DeltaRankScore")
            .expect("DeltaRankScore column missing");

        let r1: f64 = rows[0][col_idx].parse().expect("rank-1 DeltaRankScore numeric");
        let r2: f64 = rows[1][col_idx].parse().expect("rank-2 DeltaRankScore numeric");
        assert!((r1 - 7.0).abs() < 1e-6, "rank-1 DeltaRankScore should be 7.0, got {r1}");
        assert_eq!(r2, 0.0, "rank-2 DeltaRankScore should be gated to 0.0, got {r2}");
    }

    /// Verify that `longest_y_pct` is formatted with 6 decimal places.
    #[test]
    fn pin_emits_longest_y_pct_with_six_decimals() {
        let params = make_params(2..=3);
        let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

        let mut psm = make_psm(0, 10.0, 10.0, 0, 2);
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
