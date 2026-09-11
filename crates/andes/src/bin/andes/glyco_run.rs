//! The `--glyco` standalone driver: glyco scoring and `.glyco.pin` output.

use std::path::{Path, PathBuf};

use crate::cli::Cli;
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
    detected_activation_instrument: Option<(ActivationMethod, Option<InstrumentType>)>,
    // `--precursor-mono`: per-spectrum MS1 envelope corrections, aligned with
    // `spectra` (`Some` only when the flag is active; the PIN then carries the
    // Mono* columns). `None` keeps the PIN byte-identical.
    mono: Option<&[Option<search::precursor_mono::MonoCorrection>]>,
    // The isotope-error window to search instead of `params.isotope_error_range`
    // (`Some(0..=1)` once `--precursor-mono auto` corrected precursors; see
    // `mono::coupled_isotope_window`). `None` keeps the params' window.
    isotope_error_override: Option<std::ops::RangeInclusive<i8>>,
    t_total: std::time::Instant,
) -> Result<(), Box<dyn std::error::Error>> {
    let t_glyco = std::time::Instant::now();
    // Install scoring settings that reach hot inner functions. Done once here, from
    // validated CLI values, so no scoring code has to read the environment.
    scoring_crate::scoring::init_cz_settings(scoring_crate::scoring::CzSettings {
        zmax_override: cli.glyco_cz_max_charge,
    });
    andes_glyco::backbone::init_y_max_charge(cli.glyco_y_max_charge);

    // The glycan-first search space comes from one of two sources: an explicit
    // pGlyco `.gdb` file (`--glyco-glycan-gdb`, takes precedence) or a bundled
    // species-specific database (`--glyco-species`). Trees are parsed with
    // structure preserved (core- vs antenna-fucose) and NeuGc content is used
    // as-is. The composition enumerator and its NeuGc/taxon tuning are gone.
    let (glycan_list, source_desc) = match (
        cli.glyco_glycan_gdb.as_ref(),
        cli.glyco_species.as_ref(),
    ) {
        (Some(path), _) => {
            let content = std::fs::read_to_string(path)?;
            let list = andes_glyco::glycan_db::load_glycan_gdb(&content)?;
            (list, format!("{} (file)", path.display()))
        }
        (None, Some(species)) => {
            let list = andes_glyco::glycan_db::load_species_glycan_db(species.key())?;
            (list, format!("--glyco-species {}", species.key()))
        }
        (None, None) => {
            return Err(
                "--glyco requires --glyco-glycan-gdb or --glyco-species: the built-in \
                 composition enumerator was removed; supply a pGlyco `.gdb` or a bundled \
                 species"
                    .into(),
            );
        }
    };
    eprintln!(
        "glycan list: {} compositions loaded from {} (structure preserved; NeuGc content used as-is)",
        glycan_list.len(),
        source_desc
    );
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
        isotope_error_override_mask: isotope_error_override
            .as_ref()
            .and(mono)
            .map(|table| table.iter().map(Option::is_some).collect()),
        isotope_error_override,
        pf_charge: cli.glyco_pf_charge,
        // Peptide-first RETRIEVAL window. High-resolution MS2 defaults to the
        // glyco ppm tolerance; low-resolution keeps the rank model's 0.5 Da.
        // Measured 2026-09-02 (five seeds): on high-res data the 0.5 Da
        // window admitted b/y matches ~50x wider than every glycan-side matcher
        // and was the dominant glyco cost — 20 ppm was 6.9x faster on mouse
        // PXD011533 and 7x on plasma PXD030622 with identifications neutral
        // (mouse 3198 vs 3183 correct; plasma 399 vs 380). An explicit
        // --glyco-retrieval-tol-ppm overrides the auto default either way.
        retrieval_tol_ppm: retrieval_ppm,
        retrieval_tol_da: cli.glyco_retrieval_tol_da,
        max_pf: cli.glyco_max_pf,
        full_glycan_db: cli.glyco_full_glycan_db,
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
        mono,
    )?;
    eprintln!(
        "Wrote glyco PIN: {} ({} PSM rows) [PHASE TOTAL: {:.2}s]",
        glyco_pin_path.display(),
        total_glyco_rows,
        t_total.elapsed().as_secs_f64()
    );
    Ok(())
}
