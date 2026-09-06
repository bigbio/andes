//! The `--glyco` standalone driver: glyco scoring and `.glyco.pin` output.

use std::path::{Path, PathBuf};

use crate::cli::{Cli, GlycanListFlag, GlycoTaxonFlag};
use input::ProteinDb;
use model::{activation::ActivationMethod, InstrumentType, Spectrum};
use search::{PreparedSearch, SearchIndex, SearchParams};

/// Glyco mode: run glyco scoring over ALL accumulated spectra (using the
/// `PreparedSearch` from the standard search) and write a separate `.glyco.pin`.
/// The standard PIN is skipped by the caller.
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_glyco(
    cli: &Cli,
    spectra: &[Spectrum],
    prepared: &PreparedSearch,
    params: &SearchParams,
    idx: &SearchIndex,
    output_pin_path: &Path,
    spectrum_paths: &[PathBuf],
    target_db: &ProteinDb,
    detected_activation_instrument: Option<(ActivationMethod, Option<InstrumentType>)>,
    t_total: std::time::Instant,
) -> Result<(), Box<dyn std::error::Error>> {
    let t_glyco = std::time::Instant::now();
    // Use the curated common list (~600 glycans) so that ALL backbone candidates
    // can be b/y-scored in phase-1 (avoids the Y-ladder pre-filter ceiling). The
    // full ~4034-entry list (n_glycan_list()) was A/B-refuted at 1% FDR when
    // the cz collapse scoring was buggy (long/decoy candidates won). The bug
    // hunt (2026-07-16) showed the default ~612 list MISSES the mouse-brain
    // glycome at high charge (z5 69%/z6 38% coverage) — a generation ceiling.
    // ANDES_GLYCO_FULL_GLYCANS re-tests the full list now that cz is fixable.
    // Install scoring settings that reach hot inner functions. Done once here, from
    // validated CLI values, so no scoring code has to read the environment.
    scoring_crate::scoring::init_cz_settings(scoring_crate::scoring::CzSettings {
        zmax_override: cli.glyco_cz_max_charge,
    });
    andes_glyco::backbone::init_y_max_charge(cli.glyco_y_max_charge);

    // Decide whether NeuGc belongs in the search space -- and, if it does, how many per
    // composition the default list may carry. NeuGc is the sole source of isobaric mass
    // degeneracy in this list (Fuc+NeuGc and Hex+NeuAc are the SAME elemental formula),
    // so getting this right is worth more than any downstream scoring fix -- you cannot
    // resolve from fragments what should not be enumerated.
    //
    // The pair is (drop_neugc, neugc_positive). `neugc_positive` is true only when NeuGc
    // is kept ON EVIDENCE -- explicit `mammal`, or an `auto` survey that is conclusive
    // and finds NeuGc -- and false when it is merely kept for lack of signal, or by the
    // FASTA veto alone. It gates the NeuGc bound below, so that widening the search
    // space happens only where the spectra say the compositions are there.
    let (drop_neugc, neugc_positive) =
        if cli.glyco_no_neugc {
            eprintln!("--glyco-no-neugc: NeuGc excluded (explicit).");
            (true, false)
        } else {
            match cli.glyco_taxon {
                GlycoTaxonFlag::Human => {
                    eprintln!("--glyco-taxon human: NeuGc excluded (CMAH-inactivated lineage).");
                    (true, false)
                }
                GlycoTaxonFlag::Mammal => {
                    eprintln!("--glyco-taxon mammal: NeuGc kept (CMAH-competent lineage).");
                    (false, true)
                }
                GlycoTaxonFlag::Auto => {
                    // Signal 1: the spectra themselves. Stronger than the FASTA, because it
                    // measures what is in the tube rather than what was searched.
                    let survey = andes_glyco::oxonium::survey_sialic_oxonium(
                        spectra.iter().map(|s| s.peaks.as_slice()),
                        cli.glyco_tol_ppm,
                        0.10,
                    );
                    // Signal 2: OX= taxon ids in the database.
                    let (taxon, n_hu, n_nh, n_ox) = andes_glyco::glycan_db::taxon_from_headers(
                        target_db.proteins.iter().map(|p| p.description.as_str()),
                    );
                    eprintln!(
                    "--glyco-taxon auto: sialylated spectra {} ({} with NeuGc oxonium >=10% of \
                     NeuAc = {:.2}%, {}), FASTA OX= {:?} ({} CMAH-null / {} competent of {})",
                    survey.neuac_spectra,
                    survey.neugc_spectra,
                    100.0 * survey.neugc_fraction,
                    if survey.conclusive { "conclusive" } else { "INCONCLUSIVE" },
                    taxon,
                    n_hu,
                    n_nh,
                    n_ox
                );
                    // Narrow only when BOTH signals agree, and never on an inconclusive
                    // survey. A CMAH-competent FASTA vetoes: mouse genuinely has NeuGc, and
                    // recombinant human protein from a murine/CHO host does too.
                    let spectra_say_no = survey.conclusive && survey.neugc_fraction < 0.01;
                    let fasta_objects = taxon == andes_glyco::glycan_db::Taxon::CmahCompetent;
                    let decide = spectra_say_no && !fasta_objects;
                    eprintln!(
                        "--glyco-taxon auto: {} (override with --glyco-taxon human|mammal or \
                     --glyco-no-neugc)",
                        if decide {
                            "NeuGc EXCLUDED - no NeuGc oxonium evidence in this run"
                        } else if !survey.conclusive {
                            "NeuGc kept - too little sialic signal to judge"
                        } else if fasta_objects {
                            "NeuGc kept - FASTA is a CMAH-competent organism"
                        } else {
                            "NeuGc kept - NeuGc oxonium evidence present"
                        }
                    );
                    // Positive evidence uses the same 1% boundary that decides exclusion:
                    // conclusive survey, NeuGc present. A keep that only the FASTA veto
                    // produced, or an inconclusive keep, is not evidence.
                    (decide, survey.conclusive && !spectra_say_no)
                }
            }
        };
    // NeuGc bound for the default list. The shipped list is human-tuned (NeuGc <= 1); a
    // CMAH-competent sample routinely carries 2-3 NeuGc per composition, and a
    // composition that cannot be enumerated falls to the de-novo branch, which never
    // reaches the FDR PIN. Measured on pGlyco2 mouse liver T-1 (3,877 reference
    // spectra): 501 carry NeuGc >= 2 (12.9%), and 442 of the run's 838 selection losses
    // (52.7%) were exactly those -- backbone retained, glycan unnameable. Raised to
    // NeuAc's bound (4) only on positive evidence, so human runs are byte-identical: the
    // NeuGc-free subset of the list is the same at every bound.
    let max_neugc: u8 = match cli.glyco_max_neugc {
        Some(n) => n,
        None if !drop_neugc && neugc_positive => 4,
        None => 1,
    };
    let mut glycan_list = match cli.glyco_glycan_list {
        GlycanListFlag::Full => andes_glyco::glycan_db::n_glycan_list(),
        GlycanListFlag::ReferenceHuman => andes_glyco::glycan_db::n_glycan_list_reference_human(),
        GlycanListFlag::Common => {
            andes_glyco::glycan_db::n_glycan_list_common_with_neugc(max_neugc)
        }
    };
    if !drop_neugc && matches!(cli.glyco_glycan_list, GlycanListFlag::Common) {
        let why = match (cli.glyco_max_neugc, max_neugc) {
            (Some(_), _) => "explicit --glyco-max-neugc",
            (None, 1) => "human-validated default; NeuGc kept without positive evidence",
            (None, _) => "raised to NeuAc's bound: NeuGc kept on positive evidence",
        };
        eprintln!(
            "glycan list: NeuGc <= {} per composition ({}); {} compositions",
            max_neugc,
            why,
            glycan_list.len()
        );
    }
    if drop_neugc {
        let before = glycan_list.len();
        glycan_list.retain(|g| g.neugc == 0);
        eprintln!(
            "glycan list: {} -> {} compositions (NeuGc removed; it is the only source of \
             isobaric mass degeneracy here)",
            before,
            glycan_list.len()
        );
    }
    let glycan_list = glycan_list;
    let glyco_tol_ppm = cli.glyco_tol_ppm;
    // Finite and > 0 is enforced by clap (`parse_positive_tol`), so NaN and
    // non-positive values never reach here.
    if glyco_tol_ppm < 20.0 {
        eprintln!(
            "warning: --glyco-tol-ppm {glyco_tol_ppm} is tighter than the 20 ppm the \
             glyco defaults were validated at; oxonium, core-Y and c/z matching may \
             under-fire"
        );
    }
    // Dev cap on glyco scoring uses the global --max-spectra (was the
    // redundant --glyco-max-spectra).
    let spectra_for_glyco: &[_] = if cli.max_spectra > 0 {
        &spectra[..spectra.len().min(cli.max_spectra)]
    } else {
        spectra
    };
    // Peptide-first RETRIEVAL window, resolved from the ACQUISITION, not from the
    // selected scoring model: model routing deliberately sends some high-res CID
    // and ETD acquisitions to a low-res model (`build_selection_key`), and the
    // scoring window that routing chooses is not the retrieval window the index
    // should use. Detected analyzer metadata decides; metadata-less input falls
    // back to the `--fragment-tol-*` unit, the same rule the model resolver uses.
    // An explicit --glyco-retrieval-tol-{ppm,da} overrides either way.
    let acquisition_high_res = match detected_activation_instrument.and_then(|(_, i)| i) {
        Some(i) => i.is_high_resolution(),
        None => cli.fragment_tol_ppm.is_some(),
    };
    let retrieval_ppm: Option<f64> = cli.glyco_retrieval_tol_ppm.or_else(|| {
        (acquisition_high_res && cli.glyco_retrieval_tol_da.is_none()).then_some(glyco_tol_ppm)
    });
    match (cli.glyco_retrieval_tol_da, retrieval_ppm) {
        (Some(da), _) => eprintln!("glyco retrieval window: {da} Da (--glyco-retrieval-tol-da)"),
        (None, Some(ppm)) => eprintln!(
            "glyco retrieval window: {ppm} ppm ({})",
            if cli.glyco_retrieval_tol_ppm.is_some() {
                "--glyco-retrieval-tol-ppm"
            } else {
                "high-resolution acquisition"
            }
        ),
        (None, None) => eprintln!(
            // The Da fallback builds the index with the resolved MODEL fragment
            // tolerance, not a constant, so naming a number here would put a
            // fabricated provenance value into the benchmark record.
            "glyco retrieval window: model fragment tolerance in Da (low-resolution acquisition)"
        ),
    }
    let glyco_cfg = search::glyco_search::GlycoConfig {
        gp_k: cli.glyco_gp_k,
        gp_j: cli.glyco_gp_j,
        gp_h: cli.glyco_gp_h,
        gp_cz: cli.glyco_gp_cz,
        min_core_y: cli.glyco_min_core_y,
        min_raw_score: cli.glyco_min_raw_score,
        diag_splits: cli.glyco_diag_splits.clone(),
        min_matched_by: cli.glyco_min_matched_ions,
        max_gen_peaks: cli.glyco_max_peaks,
        cz_multisite: cli.glyco_cz_multisite,
        sialic_oxonium_min_frac: cli.glyco_sialic_oxonium_min_frac,
        scan_filter_path: cli.glyco_scans.clone(),
        pf_charge: cli.glyco_pf_charge,
        // Peptide-first RETRIEVAL window. High-resolution MS2 defaults to the
        // glyco ppm tolerance; low-resolution keeps the rank model's 0.5 Da.
        // Measured 2026-09-02 (Codon, five seeds): on high-res data the 0.5 Da
        // window admitted b/y matches ~50x wider than every glycan-side matcher
        // and was the dominant glyco cost — 20 ppm was 6.9x faster on mouse
        // PXD011533 and 7x on plasma PXD030622 with identifications neutral
        // (mouse 3198 vs 3183 correct; plasma 399 vs 380). An explicit
        // --glyco-retrieval-tol-ppm overrides the auto default either way.
        retrieval_tol_ppm: retrieval_ppm,
        retrieval_tol_da: cli.glyco_retrieval_tol_da,
        max_pf: cli.glyco_max_pf,
        debug: cli.debug_glyco,
        // Single-file only: cross-file pairing is unsound (see guard above).
        hcd_pair: cli.glyco_hcd_pair && spectrum_paths.len() == 1,
        etd_rank_glycan: cli.glyco_etd_rank_glycan,
        cz_gate: cli.glyco_cz_gate,
        pair_y_on_gen: cli.glyco_pair_y_on_gen,
        enum_fallback: cli.glyco_enum_fallback,
        etd_require_oxonium: cli.glyco_etd_require_oxonium,
        elect_top_k: cli.glyco_elect_top_k,
    };
    let pass1 = search::glyco_search::glyco_search_run(
        spectra_for_glyco,
        prepared,
        &glycan_list,
        glyco_tol_ppm,
        cli.glyco_backbone_top_k,
        glyco_cfg,
    );
    let total_pass1_rows: usize = pass1.iter().map(|r| r.hits.len()).sum();
    eprintln!(
        "[glyco] scored {} spectra → {} glyco-PSM rows [{:.2}s]",
        pass1.len(),
        total_pass1_rows,
        t_glyco.elapsed().as_secs_f64()
    );

    // Derive the glyco PIN path: strip a trailing `.pin` (only), then append
    // `.glyco.pin`.
    //
    // The previous form was `with_extension("").with_extension("glyco.pin")`,
    // which strips at the LAST dot TWICE — so `PXD011533.Frac1.pin` and
    // `PXD011533.Frac2.pin` both collapsed to `PXD011533.glyco.pin` and the
    // second run silently overwrote the first. Pooling fractions is the
    // recommended practice for stable glyco FDR, so `dataset.FracN.pin` is
    // exactly the naming users reach for, and the data loss was silent.
    let glyco_pin_path = {
        let mut s = if output_pin_path
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("pin"))
        {
            output_pin_path.with_extension("").into_os_string()
        } else {
            output_pin_path.to_path_buf().into_os_string()
        };
        s.push(".glyco.pin");
        std::path::PathBuf::from(s)
    };
    eprintln!("Glyco PIN will be written to: {}", glyco_pin_path.display());
    let mut glyco_results = pass1;

    // Populate glyco RT PIN features (DeltaRT/AbsDeltaRT/DeltaRTNorm +
    // predicted_rt_min) in place on each hit, using the engine-wide backbone
    // RT index + per-monosaccharide offset + per-run self-calibration. The
    // glyco PIN writer then also appends the within-scan DeltaRTRank. Neutral
    // 0.0 without observed RT / <MIN_CALIBRATION_ANCHORS anchors (baseline-safe).
    if let Some(q) = cli.glyco_min_raw_score_quantile {
        // Mutual exclusion with --glyco-min-raw-score is enforced by clap.
        let cands = &prepared.candidates;
        let is_decoy = |h: &search::glyco_search::FullGlycoPsm| -> bool {
            h.psm
                .candidate_idxs
                .first()
                .map(|&i| cands[i as usize].is_decoy)
                .unwrap_or(false)
        };
        match search::glyco_search::apply_adaptive_emission_floor(&mut glyco_results, &is_decoy, q)
        {
            Some((floor, before, kept)) => eprintln!(
                "--glyco-min-raw-score-quantile {q}: derived RawScore floor {floor:.3} \
                 from this run's decoy winners; scans {before} -> {kept}"
            ),
            None => eprintln!(
                "WARN: --glyco-min-raw-score-quantile {q} did nothing: this run has no \
                 decoy winners to calibrate on (tiny input?); emitting ungated"
            ),
        }
    }

    // Counted after the emission floor so the reported total matches the rows written.
    let total_glyco_rows: usize = glyco_results.iter().map(|r| r.hits.len()).sum();
    output::populate_glyco_rt_features(
        spectra,
        &mut glyco_results,
        &prepared.candidates,
        &glycan_list,
    );

    output::write_glyco_pin(
        cli.glyco_pin_curated,
        &glyco_pin_path,
        spectra,
        &glyco_results,
        &prepared.candidates,
        params,
        idx,
        cli.debug_glyco,
    )?;
    eprintln!(
        "Wrote glyco PIN: {} ({} PSM rows) [PHASE TOTAL: {:.2}s]",
        glyco_pin_path.display(),
        total_glyco_rows,
        t_total.elapsed().as_secs_f64()
    );
    Ok(())
}
