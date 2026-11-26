//! Applying a run config file to the search arguments: one function per config section,
//! each key applied unless the command line set it, then range checks on the result.

use super::*;

pub fn apply(
    mut cfg: RunConfig,
    args: &mut SearchArgs,
    m: &clap::ArgMatches,
    src: &mut ParamSources,
) -> Result<(), String> {
    apply_db(&cfg, &mut args.db, m, src)?;
    let xl = cfg.search.xl.take();
    let xl_linear_only = cfg.search.xl_linear_only.take();
    let mut c = ConfigApply {
        args,
        m,
        src,
        written: &cfg.written,
    };
    apply_io(cfg.io, &mut c)?;
    apply_search(cfg.search, &mut c)?;
    apply_scoring(cfg.scoring, &mut c)?;
    apply_decoys(cfg.decoys, &mut c)?;
    apply_chimeric(cfg.chimeric, &mut c)?;
    apply_refine(cfg.refine, &mut c)?;
    apply_rescoring(cfg.rescoring, &mut c)?;
    apply_quant(cfg.quant, &mut c)?;
    apply_glyco(cfg.glyco, xl, xl_linear_only, &mut c)?;
    check_search_args(c.args)
}

/// What the per-section `apply_*` functions write to: the search arguments, the command line
/// (a key set there wins over the config file), the parameter sources and the config text.
struct ConfigApply<'a> {
    args: &'a mut SearchArgs,
    m: &'a clap::ArgMatches,
    src: &'a mut ParamSources,
    written: &'a serde_yaml::Value,
}

// Each setter applies one config key unless the command line set it, and records the config
// file as its source.
macro_rules! set {
    ($c:ident, $id:literal, $($field:ident).+, $val:expr) => {
        if let Some(v) = $val {
            if !cli_set($c.m, $id) {
                $c.src.record_config($id, $c.written);
                $c.args.$($field).+ = v;
            }
        }
    };
}
macro_rules! set_opt {
    ($c:ident, $id:literal, $($field:ident).+, $val:expr) => {
        if let Some(v) = $val {
            if !cli_set($c.m, $id) {
                $c.src.record_config($id, $c.written);
                $c.args.$($field).+ = Some(v);
            }
        }
    };
}
macro_rules! set_parsed {
    ($c:ident, $id:literal, $($field:ident).+, $val:expr, $parse:expr) => {
        if let Some(s) = $val {
            if !cli_set($c.m, $id) {
                $c.args.$($field).+ = ($parse)(&s)?;
                $c.src.record_config($id, $c.written);
            }
        }
    };
}

/// Applies the config file's input and output paths to the search arguments not set on the command line.
fn apply_io(io: IoCfg, c: &mut ConfigApply<'_>) -> Result<(), String> {
    set_opt!(c, "sdrf", sdrf, io.sdrf);
    set!(c, "spectrum", spectrum, io.spectrum);
    set_opt!(c, "database", database, io.database);
    set_opt!(c, "output_pin", output_pin, io.output_pin);
    set_opt!(c, "output_tsv", output_tsv, io.output_tsv);
    set_opt!(c, "output_parquet", output_parquet, io.output_parquet);
    set_opt!(c, "temp_dir", temp_dir, io.temp_dir);
    Ok(())
}

/// Applies the config file's search settings to the search arguments not set on the command line.
fn apply_search(search: SearchCfg, c: &mut ConfigApply<'_>) -> Result<(), String> {
    set_parsed!(
        c,
        "precursor_tol",
        precursor_tol,
        search.precursor_tol,
        crate::parse_precursor_tol
    );
    set_parsed!(
        c,
        "charge",
        charge,
        search.charge,
        crate::parse_charge_range
    );
    if let Some(s) = search.isotope_error {
        if !cli_set(c.m, "isotope_error") {
            c.src.record_config("isotope_error", c.written);
            c.args.isotope_error = Some(crate::parse_isotope_error_range(&s)?);
        }
    }
    set_parsed!(
        c,
        "enzyme_specificity",
        enzyme_specificity,
        search.enzyme_specificity,
        |s: &str| crate::cli::parse_enum("enzyme_specificity", s, true)
    );
    set!(
        c,
        "max_missed_cleavages",
        max_missed_cleavages,
        search.max_missed_cleavages
    );
    set!(c, "min_length", min_length, search.min_length);
    set!(c, "max_length", max_length, search.max_length);
    set!(c, "max_mods", max_mods, search.max_mods);
    set!(c, "min_peaks", min_peaks, search.min_peaks);
    set!(c, "top_n", top_n, search.top_n);
    set!(c, "dia", dia, search.dia);
    set!(
        c,
        "dia_mobility_slices",
        dia_mobility_slices,
        search.dia_mobility_slices
    );
    set!(
        c,
        "dia_ms1_tol_ppm",
        dia_ms1_tol_ppm,
        search.dia_ms1_tol_ppm
    );
    if let Some(v) = search.localize.filter(|_| !cli_set(c.m, "localize")) {
        if v != "ipf" {
            return Err(format!("search.localize: '{v}' is not one of: ipf"));
        }
        c.args.localize = Some(v);
        c.src.record_config("localize", c.written);
    }
    set_parsed!(
        c,
        "fragment_index",
        fragment_index,
        search.fragment_index,
        parse_value_enum::<crate::cli::FragmentIndexFlag>
    );
    set_opt!(
        c,
        "fragment_index_slice_da",
        fragment_index_slice_da,
        search.fragment_index_slice_da
    );
    set_opt!(
        c,
        "index_cache_dir",
        index_cache_dir,
        search.index_cache_dir
    );
    set!(c, "open_search", open_search, search.open_search);
    set_parsed!(
        c,
        "open_search_delta",
        open_search_delta,
        search.open_search_delta,
        crate::parse_open_search_delta
    );
    set!(
        c,
        "open_search_top_k",
        open_search_top_k,
        search.open_search_top_k
    );
    set!(c, "optimize", optimize, search.optimize);
    set!(
        c,
        "optimize_sample",
        optimize_sample,
        search.optimize_sample
    );
    Ok(())
}

/// Applies the config file's scoring settings to the search arguments not set on the command line.
fn apply_scoring(scoring: ScoringCfg, c: &mut ConfigApply<'_>) -> Result<(), String> {
    set_parsed!(
        c,
        "score",
        score,
        scoring.score,
        parse_value_enum::<crate::ScoreFlag>
    );
    set_parsed!(
        c,
        "precursor_cal",
        precursor_cal,
        scoring.precursor_cal,
        crate::parse_precursor_cal
    );
    set!(c, "threads", threads, scoring.threads);
    set_parsed!(
        c,
        "fragmentation",
        fragmentation,
        scoring.fragmentation,
        |s: &str| crate::cli::parse_enum("fragmentation", s, true)
    );
    set_parsed!(c, "protocol", protocol, scoring.protocol, |s: &str| {
        crate::cli::parse_enum("protocol", s, true)
    });
    if !cli_set(c.m, "fragment_tol_ppm") && !cli_set(c.m, "fragment_tol_da") {
        set_opt!(
            c,
            "fragment_tol_ppm",
            fragment_tol_ppm,
            scoring.fragment_tol_ppm
        );
        set_opt!(
            c,
            "fragment_tol_da",
            fragment_tol_da,
            scoring.fragment_tol_da
        );
    }
    set_opt!(c, "model_store", model_store, scoring.model_store);
    set_opt!(c, "model_id_override", model_id_override, scoring.model);
    set!(c, "ms_level", ms_level, scoring.ms_level);
    set!(c, "max_spectra", max_spectra, scoring.max_spectra);
    set_opt!(c, "gbdt_max_trees", gbdt_max_trees, scoring.gbdt_max_trees);
    set_parsed!(
        c,
        "precursor_mono",
        precursor_mono,
        scoring.precursor_mono,
        parse_value_enum::<crate::cli::PrecursorMonoFlag>
    );
    set_parsed!(
        c,
        "ethcd_activation",
        ethcd_activation,
        scoring.ethcd_activation,
        parse_value_enum::<crate::cli::EthcdActivationFlag>
    );
    set_opt!(c, "peak_filter", peak_filter, scoring.peak_filter);
    Ok(())
}

/// Applies the config file's decoy settings to the search arguments not set on the command line.
fn apply_decoys(decoys: DecoysCfg, c: &mut ConfigApply<'_>) -> Result<(), String> {
    set_parsed!(
        c,
        "decoy_strategy",
        decoy_strategy,
        decoys.strategy,
        crate::cli::parse_decoy_strategy
    );
    set!(c, "decoy_seed", decoy_seed, decoys.seed);
    Ok(())
}

/// Applies the config file's chimeric-search settings to the search arguments not set on the command line.
fn apply_chimeric(chimeric: ChimericCfg, c: &mut ConfigApply<'_>) -> Result<(), String> {
    set!(c, "chimeric", chimeric, chimeric.enabled);
    set!(
        c,
        "chimeric_max_coisolated",
        chimeric_max_coisolated,
        chimeric.max_coisolated
    );
    set!(c, "chimeric_max_kl", chimeric_max_kl, chimeric.max_kl);
    Ok(())
}

/// Applies the config file's refine settings to the search arguments not set on the command line.
fn apply_refine(refine: RefineCfg, c: &mut ConfigApply<'_>) -> Result<(), String> {
    set!(c, "refine", refine, refine.enabled);
    set_opt!(c, "refine_config", refine_config, refine.config);
    Ok(())
}

/// Applies the config file's rescoring and FDR settings to the search arguments not set on the command line.
fn apply_rescoring(rescoring: RescoringCfg, c: &mut ConfigApply<'_>) -> Result<(), String> {
    set!(c, "rescore", rescore, rescoring.enabled);
    set_opt!(c, "psm_fdr", psm_fdr, rescoring.psm_fdr);
    set_opt!(
        c,
        "proteoform_fdr",
        proteoform_fdr,
        rescoring.proteoform_fdr
    );
    set_opt!(c, "protein_fdr", protein_fdr, rescoring.protein_fdr);
    set!(c, "add_entrapment", add_entrapment, rescoring.add_entrapment);
    set!(c, "entrapment_seed", entrapment_seed, rescoring.entrapment_seed);
    set_opt!(c, "pep", pep, rescoring.pep);
    set_opt!(
        c,
        "percolator_bin",
        percolator_bin,
        rescoring.percolator_bin
    );
    set!(
        c,
        "percolator_docker",
        percolator_docker,
        rescoring.percolator_docker
    );
    set!(
        c,
        "percolator_image",
        percolator_image,
        rescoring.percolator_image
    );
    set!(
        c,
        "percolator_args",
        percolator_args,
        rescoring.percolator_args
    );
    set!(c, "keep_pin", keep_pin, rescoring.keep_pin);
    set_opt!(
        c,
        "rescore_memory_gb",
        rescore_memory_gb,
        rescoring.memory_gb
    );
    if let Some(list) = rescoring.predict {
        if !cli_set(c.m, "predict") {
            c.args.predict = list
                .iter()
                .map(|s| crate::cli::parse_enum("predict entry", s, true))
                .collect::<Result<_, _>>()?;
            c.src.record_config("predict", c.written);
        }
    }
    set_opt!(
        c,
        "predictor_models",
        predictor_models,
        rescoring.predictor_models
    );
    set_parsed!(
        c,
        "predictor_backend",
        predictor_backend,
        rescoring.predictor_backend,
        |s: &str| crate::cli::parse_enum("predictor_backend", s, true)
    );
    set_opt!(c, "predictor_nce", predictor_nce, rescoring.predictor_nce);
    set_opt!(
        c,
        "predictor_instrument",
        predictor_instrument,
        rescoring.predictor_instrument
    );
    set!(
        c,
        "predictor_batch",
        predictor_batch,
        rescoring.predictor_batch
    );
    Ok(())
}

/// Applies the config file's quantification settings to the search arguments not set on the command line.
fn apply_quant(quant: QuantCfg, c: &mut ConfigApply<'_>) -> Result<(), String> {
    if let Some(v) = quant.tmt {
        if !cli_set(c.m, "tmt") {
            c.src.record_config("tmt", c.written);
            c.args.tmt = Some(crate::cli::parse_plex(&v)?);
        }
    }
    set_opt!(c, "tmt_level", tmt_level, quant.tmt_level);
    if let Some(v) = quant.tmt_tol {
        if !cli_set(c.m, "tmt_tol") {
            c.src.record_config("tmt_tol", c.written);
            c.args.tmt_tol = Some(parse_tol("tmt_tol", &v)?);
        }
    }
    set_opt!(c, "tmt_correction", tmt_correction, quant.tmt_correction);
    set!(c, "tmt_min_purity", tmt_min_purity, quant.tmt_min_purity);
    set!(c, "tmt_min_sn", tmt_min_sn, quant.tmt_min_sn);
    set_parsed!(
        c,
        "tmt_unmix",
        tmt_unmix,
        quant.tmt_unmix,
        parse_value_enum::<crate::cli::TmtUnmix>
    );
    set!(c, "lfq", lfq, quant.lfq);
    if let Some(v) = quant.lfq_tol {
        if !cli_set(c.m, "lfq_tol") {
            c.src.record_config("lfq_tol", c.written);
            c.args.lfq_tol = Some(parse_tol("lfq_tol", &v)?);
        }
    }
    set_opt!(c, "lfq_rt_window", lfq_rt_window, quant.lfq_rt_window);
    set_opt!(
        c,
        "lfq_gaussian_width",
        lfq_gaussian_width,
        quant.lfq_gaussian_width
    );
    set!(c, "lfq_min_cosine", lfq_min_cosine, quant.lfq_min_cosine);
    set_opt!(c, "label", label, quant.label);
    set!(c, "label_reference", label_reference, quant.label_reference);
    set!(c, "label_arg_pro", label_arg_pro, quant.label_arg_pro);
    set!(
        c,
        "label_2h_enrichment_silac",
        label_2h_enrichment_silac,
        quant.label_2h_enrichment_silac
    );
    set!(
        c,
        "label_2h_enrichment_chemical",
        label_2h_enrichment_chemical,
        quant.label_2h_enrichment_chemical
    );
    set!(
        c,
        "label_heavy_enrichment",
        label_heavy_enrichment,
        quant.label_heavy_enrichment
    );
    set_opt!(c, "quant_fdr", quant_fdr, quant.quant_fdr);
    set_opt!(
        c,
        "quant_transfer_fdr",
        quant_transfer_fdr,
        quant.quant_transfer_fdr
    );
    set_opt!(
        c,
        "quant_channel_fdr",
        quant_channel_fdr,
        quant.quant_channel_fdr
    );
    set_opt!(c, "quant_dir", quant_dir, quant.dir);
    set!(c, "mbr", mbr, quant.mbr);
    set_opt!(c, "mbr_donor_fdr", mbr_donor_fdr, quant.mbr_donor_fdr);
    set!(c, "mbr_wrong_log2", mbr_wrong_log2, quant.mbr_wrong_log2);
    Ok(())
}

/// Applies the config file's glyco settings, with the cross-link keys of `search` that sit among them to the search arguments not set on the command line.
fn apply_glyco(
    glyco: GlycoCfg,
    xl: Option<serde_yaml::Value>,
    xl_linear_only: Option<bool>,
    c: &mut ConfigApply<'_>,
) -> Result<(), String> {
    set!(c, "glyco", glyco, glyco.enabled);
    if let Some(v) = &xl {
        if !cli_set(c.m, "xl") {
            c.src.record_config("xl", c.written);
            c.args.xl = Some(xl_set(v)?);
        }
    }
    set!(c, "xl_linear_only", xl_linear_only, xl_linear_only);
    set_parsed!(c, "glyco_core", glyco_core, glyco.core, parse_glyco_core);
    set!(c, "glyco_tol_ppm", glyco_tol_ppm, glyco.tol_ppm);
    set!(c, "glyco_max_peaks", glyco_max_peaks, glyco.max_peaks);
    set!(c, "debug_glyco", debug_glyco, glyco.debug);
    set!(c, "glyco_hcd_pair", glyco_hcd_pair, glyco.hcd_pair);
    if !cli_set(c.m, "glyco_glycan_gdb")
        && !cli_set(c.m, "glyco_glycan_wurcs")
        && !cli_set(c.m, "glyco_species")
    {
        let n_sources = usize::from(glyco.glycan_gdb.is_some())
            + usize::from(glyco.glycan_wurcs.is_some())
            + usize::from(glyco.species.is_some());
        if n_sources > 1 {
            return Err(
                "config: `glyco.glycan_gdb`, `glyco.glycan_wurcs` and `glyco.species` are \
                 mutually exclusive — set only one"
                    .to_string(),
            );
        }
        if let Some(p) = glyco.glycan_gdb {
            c.src.record_config("glyco_glycan_gdb", c.written);
            c.args.glyco_glycan_gdb = Some(p);
        }
        if let Some(p) = glyco.glycan_wurcs {
            c.src.record_config("glyco_glycan_wurcs", c.written);
            c.args.glyco_glycan_wurcs = Some(p);
        }
        if let Some(s) = glyco.species {
            c.src.record_config("glyco_species", c.written);
            c.args.glyco_species = Some(
                parse_value_enum::<crate::cli::GlycoSpeciesFlag>(&s)
                    .map_err(|e| format!("config: `glyco.species`: {e}"))?,
            );
        }
    }
    set!(c, "glyco_min_core_y", glyco_min_core_y, glyco.min_core_y);
    set_opt!(
        c,
        "glyco_retrieval_tol_ppm",
        glyco_retrieval_tol_ppm,
        glyco.retrieval_tol_ppm
    );
    set!(
        c,
        "glyco_peptide_first",
        glyco_peptide_first,
        glyco.peptide_first
    );
    set!(
        c,
        "glyco_index_sequon_only",
        glyco_index_sequon_only,
        glyco.index_sequon_only
    );
    set_opt!(
        c,
        "glyco_cz_max_charge",
        glyco_cz_max_charge,
        glyco.cz_max_charge
    );
    set!(
        c,
        "glyco_y_max_charge",
        glyco_y_max_charge,
        glyco.y_max_charge
    );
    Ok(())
}

/// Range checks on the merged search arguments, after the config file is applied.
fn check_search_args(args: &SearchArgs) -> Result<(), String> {
    if let Some(v) = args.psm_fdr {
        check_unit_fraction("psm_fdr", v)?;
    }
    if args.rescore_memory_gb == Some(0) {
        return Err("config: `memory_gb` must be at least 1 (GB)".to_string());
    }
    if let Some(v) = args.proteoform_fdr {
        check_unit_fraction("proteoform_fdr", v)?;
    }
    if let Some(v) = args.protein_fdr {
        check_unit_fraction("protein_fdr", v)?;
    }
    if let Some(v) = args.pep {
        check_unit_fraction("pep", v)?;
    }
    for (key, value) in [
        ("mbr_donor_fdr", args.mbr_donor_fdr),
        ("quant_transfer_fdr", args.quant_transfer_fdr),
        ("quant_channel_fdr", args.quant_channel_fdr),
    ] {
        if let Some(v) = value {
            check_unit_fraction(key, v)?;
        }
    }
    check_unit_fraction("tmt_min_purity", args.tmt_min_purity)?;
    if !args.tmt_min_sn.is_finite() || args.tmt_min_sn < 0.0 {
        return Err(format!(
            "config: tmt_min_sn must be >= 0, got {}",
            args.tmt_min_sn
        ));
    }
    check_unit_fraction("lfq_min_cosine", args.lfq_min_cosine)?;
    check_unit_fraction("label_2h_enrichment_silac", args.label_2h_enrichment_silac)?;
    check_unit_fraction(
        "label_2h_enrichment_chemical",
        args.label_2h_enrichment_chemical,
    )?;
    check_unit_fraction("label_heavy_enrichment", args.label_heavy_enrichment)?;
    if let Some(v) = args.quant_fdr {
        check_unit_fraction("quant_fdr", v)?;
    }
    if let Some(v) = args.lfq_rt_window {
        check_positive("lfq_rt_window", v)?;
    }
    check_positive("mbr_wrong_log2", args.mbr_wrong_log2)?;
    if let Some(v) = args.lfq_gaussian_width {
        check_positive("lfq_gaussian_width", v)?;
    }
    check_positive("glyco_tol_ppm", args.glyco_tol_ppm)?;
    if let Some(v) = args.glyco_retrieval_tol_ppm {
        check_positive("glyco_retrieval_tol_ppm", v)?;
    }
    if !args.chimeric_max_kl.is_finite() || args.chimeric_max_kl < 0.0 {
        return Err(format!(
            "config: `chimeric_max_kl` must be finite and >= 0 (got {})",
            args.chimeric_max_kl
        ));
    }
    if let Some(v) = args.fragment_index_slice_da {
        check_positive("fragment_index_slice_da", v)?;
    }
    if let Some(v) = args.fragment_tol_ppm {
        check_positive("fragment_tol_ppm", v)?;
    }
    if let Some(v) = args.fragment_tol_da {
        check_positive("fragment_tol_da", v)?;
    }
    if args.fragment_tol_ppm.is_some() && args.fragment_tol_da.is_some() {
        return Err(
            "config: `scoring.fragment_tol_ppm` and `scoring.fragment_tol_da` are mutually \
             exclusive — set only one"
                .to_string(),
        );
    }

    Ok(())
}
