//! YAML run-configuration file (`--config run.yaml`).
//!
//! A single YAML file can set every search parameter, grouped by experimental
//! design to mirror the DOCS §1a taxonomy (`io`, `search`, `scoring`, `decoys`,
//! `chimeric`, `refine`, `rescoring`, `glyco`). Every key is optional; omitted
//! keys keep their built-in default.
//!
//! PRECEDENCE: an explicit CLI flag ALWAYS wins over the YAML value, which wins
//! over the built-in default. We detect "explicit CLI flag" with clap's
//! `ValueSource::CommandLine`, so the YAML only fills parameters the user did
//! not type on the command line. Unknown keys are a hard error (`deny_unknown_
//! fields`) so a typo never silently no-ops.
//!
//! Every search flag has a key here except four that are CLI-only on purpose:
//! `--config` itself, and the hidden diagnostics `--cal-min-spec-keys` (a test
//! hook that forces the calibration pre-pass on tiny inputs),
//! `--precursor-mono-dump` and `--glyco-scans` (per-run debug dumps). The test
//! `every_search_flag_has_a_config_key` keeps that list honest.
//!
//! Non-scalar values (tolerances, ranges, enums) are written as the SAME strings
//! the CLI accepts (e.g. `precursor_tol: 20ppm`, `charge: "2..5"`, `score: auto`)
//! and are parsed through the identical validators, so the config and the CLI can
//! never diverge in what they accept.

use std::path::PathBuf;

use serde::Deserialize;

use crate::SearchArgs;

/// Top-level run config; sections mirror DOCS §1a.
#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct RunConfig {
    pub io: IoCfg,
    pub search: SearchCfg,
    pub scoring: ScoringCfg,
    pub decoys: DecoysCfg,
    pub chimeric: ChimericCfg,
    pub refine: RefineCfg,
    pub rescoring: RescoringCfg,
    pub glyco: GlycoCfg,
    #[serde(default)]
    pub quant: QuantCfg,
}

#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct IoCfg {
    pub spectrum: Option<Vec<PathBuf>>,
    pub database: Option<PathBuf>,
    pub uniprot_variants: Option<bool>,
    pub uniprot_isoforms: Option<bool>,
    pub uniprot_chains: Option<bool>,
    pub uniprot_evidence: Option<String>,
    pub annotations: Option<Vec<PathBuf>>,
    pub knockoffs: Option<usize>,
    pub output_pin: Option<PathBuf>,
    pub output_tsv: Option<PathBuf>,
    pub output_parquet: Option<PathBuf>,
    pub temp_dir: Option<PathBuf>,
}

#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct SearchCfg {
    pub precursor_tol: Option<String>,
    pub charge: Option<String>,
    pub isotope_error: Option<String>,
    pub enzyme: Option<String>,
    pub enzyme_specificity: Option<String>,
    pub max_missed_cleavages: Option<u32>,
    pub min_length: Option<u32>,
    pub max_length: Option<u32>,
    pub max_mods: Option<u32>,
    pub min_peaks: Option<u32>,
    pub top_n: Option<u32>,
    pub mods: Option<PathBuf>,
    /// `auto`, `on` or `off`, as `--fragment-index`.
    pub fragment_index: Option<String>,
    pub fragment_index_slice_da: Option<f64>,
    pub index_cache_dir: Option<PathBuf>,
}

#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct ScoringCfg {
    pub score: Option<String>,
    pub precursor_cal: Option<String>,
    pub threads: Option<usize>,
    pub fragmentation: Option<String>,
    pub protocol: Option<String>,
    pub fragment_tol_ppm: Option<f64>,
    pub fragment_tol_da: Option<f64>,
    pub model_store: Option<PathBuf>,
    pub model: Option<String>,
    pub ms_level: Option<u8>,
    pub max_spectra: Option<usize>,
    pub gbdt_max_trees: Option<usize>,
    /// `auto` or `off`, as `--precursor-mono`.
    pub precursor_mono: Option<String>,
    /// `hcd` or `etd`, as `--ethcd-activation`.
    pub ethcd_activation: Option<String>,
    /// `WINDOW_DA:PEAKS`, as `--peak-filter`.
    pub peak_filter: Option<String>,
}

#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct DecoysCfg {
    pub strategy: Option<String>,
    pub prefix: Option<String>,
    pub suffix: Option<String>,
    pub seed: Option<u64>,
}

#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct ChimericCfg {
    pub enabled: Option<bool>,
    pub max_coisolated: Option<usize>,
    pub max_kl: Option<f32>,
}

#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct RefineCfg {
    pub enabled: Option<bool>,
    pub config: Option<PathBuf>,
    pub select_psm_fdr: Option<f64>,
}

#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct RescoringCfg {
    pub enabled: Option<bool>,
    pub native: Option<bool>,
    pub fdr: Option<f64>,
    pub proteoform_fdr: Option<f64>,
    pub protein_fdr: Option<f64>,
    pub pep: Option<f64>,
    pub percolator_bin: Option<PathBuf>,
    pub percolator_docker: Option<bool>,
    pub percolator_image: Option<String>,
    pub percolator_args: Option<String>,
    pub keep_pin: Option<bool>,
    pub predict: Option<Vec<String>>,
    pub predictor_models: Option<PathBuf>,
    pub predictor_backend: Option<String>,
    pub predictor_nce: Option<f32>,
    pub predictor_instrument: Option<String>,
    pub predictor_batch: Option<usize>,
}

#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct QuantCfg {
    pub tmt: Option<String>,
    pub tmt_level: Option<u8>,
    pub tmt_tol: Option<String>,
    pub tmt_correction: Option<PathBuf>,
    pub tmt_min_purity: Option<f64>,
    pub tmt_max_pep: Option<f64>,
    pub tmt_min_sn: Option<f64>,
    pub tmt_unmix: Option<bool>,
    pub lfq: Option<bool>,
    pub lfq_tol: Option<String>,
    pub lfq_rt_window: Option<f64>,
    pub lfq_gaussian_width: Option<f64>,
    pub lfq_min_cosine: Option<f64>,
    pub lfq_feature_fdr: Option<f64>,
    pub quant_fdr: Option<f64>,
    pub dir: Option<PathBuf>,
    pub mbr: Option<bool>,
    pub mbr_library_fdr: Option<f64>,
    pub mbr_fdr: Option<f64>,
    /// A number in `(0, 1]` or `off`, as `--quant-feature-fdr`.
    pub quant_feature_fdr: Option<NumberOrString>,
}

/// A YAML scalar written as a number or as a string (e.g. `0.05` or `off`),
/// passed to the CLI validator as its text.
#[derive(Debug, Clone, Deserialize)]
#[serde(untagged)]
pub enum NumberOrString {
    Number(f64),
    String(String),
}

/// `quant_feature_fdr` through the CLI validator, the error naming the key.
fn parse_quant_feature_fdr(s: &str) -> Result<crate::cli::FractionOrOff, String> {
    crate::cli::parse_fraction_or_off(s).map_err(|e| format!("quant_feature_fdr: {e}"))
}

fn parse_glyco_core(s: &str) -> Result<crate::cli::GlycoCoreFlag, String> {
    <crate::cli::GlycoCoreFlag as clap::ValueEnum>::from_str(s, true)
        .map_err(|_| format!("glyco.core: invalid `{s}`: expected n|o"))
}

impl NumberOrString {
    fn text(&self) -> String {
        match self {
            NumberOrString::Number(v) => v.to_string(),
            NumberOrString::String(s) => s.clone(),
        }
    }
}

#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct GlycoCfg {
    pub enabled: Option<bool>,
    pub core: Option<String>,
    pub tol_ppm: Option<f64>,
    pub max_peaks: Option<usize>,
    pub debug: Option<bool>,
    pub hcd_pair: Option<bool>,
    /// External `.gdb` glycan database, as `--glyco-glycan-gdb`.
    pub glycan_gdb: Option<PathBuf>,
    /// WURCS 2.0 glycan database CSV, as `--glyco-glycan-wurcs`.
    pub glycan_wurcs: Option<PathBuf>,
    /// Bundled species database, as `--glyco-species`.
    pub species: Option<String>,
    pub min_core_y: Option<u32>,
    pub retrieval_tol_ppm: Option<f64>,
    pub peptide_first: Option<bool>,
    pub index_sequon_only: Option<bool>,
    pub cz_max_charge: Option<u8>,
    pub y_max_charge: Option<u8>,
}

impl RunConfig {
    /// Parse a YAML config string. Unknown keys are rejected.
    pub fn from_yaml_str(s: &str) -> Result<Self, String> {
        serde_yaml::from_str(s).map_err(|e| format!("--config parse error: {e}"))
    }

    /// Load a YAML config file.
    pub fn load(path: &std::path::Path) -> Result<Self, String> {
        let s = std::fs::read_to_string(path)
            .map_err(|e| format!("--config: cannot read {}: {e}", path.display()))?;
        Self::from_yaml_str(&s)
    }
}

/// True when the flag `id` was set EXPLICITLY on the command line (not a default
/// and not an env fallback), so the YAML must not override it.
fn cli_set(m: &clap::ArgMatches, id: &str) -> bool {
    matches!(
        m.value_source(id),
        Some(clap::parser::ValueSource::CommandLine)
    )
}

/// Apply the YAML config onto the parsed `SearchArgs`, respecting CLI precedence.
///
/// For each parameter: if the user did NOT type the flag (`!cli_set`) and the
/// YAML provides a value, the YAML value replaces the built-in default. Values
/// written as strings (tolerances/ranges/enums) go through the exact CLI parsers.
pub fn apply(cfg: RunConfig, args: &mut SearchArgs, m: &clap::ArgMatches) -> Result<(), String> {
    // scalar field, value type == target type
    macro_rules! set {
        ($id:literal, $field:ident, $val:expr) => {
            if let Some(v) = $val {
                if !cli_set(m, $id) {
                    args.$field = v;
                }
            }
        };
    }
    // scalar into an Option<> target
    macro_rules! set_opt {
        ($id:literal, $field:ident, $val:expr) => {
            if let Some(v) = $val {
                if !cli_set(m, $id) {
                    args.$field = Some(v);
                }
            }
        };
    }
    // string wire value parsed through a CLI validator `fn(&str) -> Result<T,String>`
    macro_rules! set_parsed {
        ($id:literal, $field:ident, $val:expr, $parse:path) => {
            if let Some(s) = $val {
                if !cli_set(m, $id) {
                    args.$field = $parse(&s)?;
                }
            }
        };
    }

    // ── io ──
    set!("spectrum", spectrum, cfg.io.spectrum);
    set_opt!("database", database, cfg.io.database);
    set!(
        "uniprot_variants",
        uniprot_variants,
        cfg.io.uniprot_variants
    );
    set!(
        "uniprot_isoforms",
        uniprot_isoforms,
        cfg.io.uniprot_isoforms
    );
    set!("uniprot_chains", uniprot_chains, cfg.io.uniprot_chains);
    set!("annotations", annotations, cfg.io.annotations);
    set!("knockoffs", knockoffs, cfg.io.knockoffs);
    set_parsed!(
        "uniprot_evidence",
        uniprot_evidence,
        cfg.io.uniprot_evidence,
        crate::parse_uniprot_evidence
    );
    set_opt!("output_pin", output_pin, cfg.io.output_pin);
    set_opt!("output_tsv", output_tsv, cfg.io.output_tsv);
    set_opt!("output_parquet", output_parquet, cfg.io.output_parquet);
    set_opt!("temp_dir", temp_dir, cfg.io.temp_dir);

    // ── search ──
    set_parsed!(
        "precursor_tol",
        precursor_tol,
        cfg.search.precursor_tol,
        crate::parse_precursor_tol
    );
    set_parsed!(
        "charge",
        charge,
        cfg.search.charge,
        crate::parse_charge_range
    );
    // isotope_error is Option<(i8,i8)> on the CLI so an EXPLICIT range is
    // distinguishable from the default (which is mode-dependent: 0..2 under --glyco).
    // A config-file value counts as explicit, so wrap it in Some.
    if let Some(s) = cfg.search.isotope_error {
        if !cli_set(m, "isotope_error") {
            args.isotope_error = Some(crate::parse_isotope_error_range(&s)?);
        }
    }
    set!("enzyme", enzyme, cfg.search.enzyme);
    set_parsed!(
        "enzyme_specificity",
        enzyme_specificity,
        cfg.search.enzyme_specificity,
        crate::parse_enzyme_specificity
    );
    set!(
        "max_missed_cleavages",
        max_missed_cleavages,
        cfg.search.max_missed_cleavages
    );
    set!("min_length", min_length, cfg.search.min_length);
    set!("max_length", max_length, cfg.search.max_length);
    set!("max_mods", max_mods, cfg.search.max_mods);
    set!("min_peaks", min_peaks, cfg.search.min_peaks);
    set!("top_n", top_n, cfg.search.top_n);
    set_opt!("mods", mods, cfg.search.mods);
    set_parsed!(
        "fragment_index",
        fragment_index,
        cfg.search.fragment_index,
        parse_value_enum::<crate::cli::FragmentIndexFlag>
    );
    set_opt!(
        "fragment_index_slice_da",
        fragment_index_slice_da,
        cfg.search.fragment_index_slice_da
    );
    set_opt!(
        "index_cache_dir",
        index_cache_dir,
        cfg.search.index_cache_dir
    );

    // ── scoring ──
    set_parsed!(
        "score",
        score,
        cfg.scoring.score,
        parse_value_enum::<crate::ScoreFlag>
    );
    set_parsed!(
        "precursor_cal",
        precursor_cal,
        cfg.scoring.precursor_cal,
        crate::parse_precursor_cal
    );
    set!("threads", threads, cfg.scoring.threads);
    set_parsed!(
        "fragmentation",
        fragmentation,
        cfg.scoring.fragmentation,
        crate::parse_fragmentation
    );
    set_parsed!(
        "protocol",
        protocol,
        cfg.scoring.protocol,
        crate::parse_protocol
    );
    // The fragment tolerance is one choice: a `--fragment-tol-ppm` or
    // `--fragment-tol-da` on the command line replaces both config keys, so a
    // CLI ppm value never meets a config Da value (or the reverse).
    if !cli_set(m, "fragment_tol_ppm") && !cli_set(m, "fragment_tol_da") {
        set_opt!(
            "fragment_tol_ppm",
            fragment_tol_ppm,
            cfg.scoring.fragment_tol_ppm
        );
        set_opt!(
            "fragment_tol_da",
            fragment_tol_da,
            cfg.scoring.fragment_tol_da
        );
    }
    set_opt!("model_store", model_store, cfg.scoring.model_store);
    set_opt!("model_id_override", model_id_override, cfg.scoring.model);
    set!("ms_level", ms_level, cfg.scoring.ms_level);
    set!("max_spectra", max_spectra, cfg.scoring.max_spectra);
    set_opt!("gbdt_max_trees", gbdt_max_trees, cfg.scoring.gbdt_max_trees);
    set_parsed!(
        "precursor_mono",
        precursor_mono,
        cfg.scoring.precursor_mono,
        parse_value_enum::<crate::cli::PrecursorMonoFlag>
    );
    set_parsed!(
        "ethcd_activation",
        ethcd_activation,
        cfg.scoring.ethcd_activation,
        parse_value_enum::<crate::cli::EthcdActivationFlag>
    );
    set_opt!("peak_filter", peak_filter, cfg.scoring.peak_filter);

    // ── decoys ──
    set!("decoy_strategy", decoy_strategy, cfg.decoys.strategy);
    set!("decoy_prefix", decoy_prefix, cfg.decoys.prefix);
    set_opt!("decoy_suffix", decoy_suffix, cfg.decoys.suffix);
    set!("decoy_seed", decoy_seed, cfg.decoys.seed);

    // ── chimeric ──
    set!("chimeric", chimeric, cfg.chimeric.enabled);
    set!(
        "chimeric_max_coisolated",
        chimeric_max_coisolated,
        cfg.chimeric.max_coisolated
    );
    set!("chimeric_max_kl", chimeric_max_kl, cfg.chimeric.max_kl);

    // ── refine ──
    set!("refine", refine, cfg.refine.enabled);
    set_opt!("refine_config", refine_config, cfg.refine.config);
    set!(
        "refine_select_psm_fdr",
        refine_select_psm_fdr,
        cfg.refine.select_psm_fdr
    );

    // ── rescoring ──
    // The rescoring backend is one choice: `--rescore` or `--rescore-native` on the
    // command line replaces the config's `rescoring.enabled` / `rescoring.native`
    // entirely, so a CLI `--rescore-native` is not overridden by a config that
    // enables Percolator (which wins when both are set).
    if !cli_set(m, "rescore") && !cli_set(m, "rescore_native") {
        set!("rescore", rescore, cfg.rescoring.enabled);
        set!("rescore_native", rescore_native, cfg.rescoring.native);
    }
    set_opt!("fdr", fdr, cfg.rescoring.fdr);
    set_opt!(
        "proteoform_fdr",
        proteoform_fdr,
        cfg.rescoring.proteoform_fdr
    );
    set_opt!("protein_fdr", protein_fdr, cfg.rescoring.protein_fdr);
    set_opt!("pep", pep, cfg.rescoring.pep);
    set_opt!(
        "percolator_bin",
        percolator_bin,
        cfg.rescoring.percolator_bin
    );
    set!(
        "percolator_docker",
        percolator_docker,
        cfg.rescoring.percolator_docker
    );
    set!(
        "percolator_image",
        percolator_image,
        cfg.rescoring.percolator_image
    );
    set!(
        "percolator_args",
        percolator_args,
        cfg.rescoring.percolator_args
    );
    set!("keep_pin", keep_pin, cfg.rescoring.keep_pin);
    if let Some(list) = cfg.rescoring.predict {
        if !cli_set(m, "predict") {
            args.predict = list
                .iter()
                .map(|s| crate::cli::parse_predict_feature(s))
                .collect::<Result<_, _>>()?;
        }
    }
    set_opt!(
        "predictor_models",
        predictor_models,
        cfg.rescoring.predictor_models
    );
    set_parsed!(
        "predictor_backend",
        predictor_backend,
        cfg.rescoring.predictor_backend,
        crate::cli::parse_predictor_backend
    );
    set_opt!("predictor_nce", predictor_nce, cfg.rescoring.predictor_nce);
    set_opt!(
        "predictor_instrument",
        predictor_instrument,
        cfg.rescoring.predictor_instrument
    );
    set!(
        "predictor_batch",
        predictor_batch,
        cfg.rescoring.predictor_batch
    );

    // ── quant ──
    if let Some(v) = cfg.quant.tmt {
        if !cli_set(m, "tmt") {
            args.tmt = Some(crate::cli::parse_plex(&v)?);
        }
    }
    set_opt!("tmt_level", tmt_level, cfg.quant.tmt_level);
    if let Some(v) = cfg.quant.tmt_tol {
        if !cli_set(m, "tmt_tol") {
            args.tmt_tol = Some(parse_tol("tmt_tol", &v)?);
        }
    }
    set_opt!("tmt_correction", tmt_correction, cfg.quant.tmt_correction);
    set!("tmt_min_purity", tmt_min_purity, cfg.quant.tmt_min_purity);
    set!("tmt_max_pep", tmt_max_pep, cfg.quant.tmt_max_pep);
    set!("tmt_min_sn", tmt_min_sn, cfg.quant.tmt_min_sn);
    set!("tmt_unmix", tmt_unmix, cfg.quant.tmt_unmix);
    set!("lfq", lfq, cfg.quant.lfq);
    if let Some(v) = cfg.quant.lfq_tol {
        if !cli_set(m, "lfq_tol") {
            args.lfq_tol = Some(parse_tol("lfq_tol", &v)?);
        }
    }
    set!("lfq_rt_window", lfq_rt_window, cfg.quant.lfq_rt_window);
    set_opt!(
        "lfq_gaussian_width",
        lfq_gaussian_width,
        cfg.quant.lfq_gaussian_width
    );
    set!("lfq_min_cosine", lfq_min_cosine, cfg.quant.lfq_min_cosine);
    set!(
        "lfq_feature_fdr",
        lfq_feature_fdr,
        cfg.quant.lfq_feature_fdr
    );
    set_opt!("quant_fdr", quant_fdr, cfg.quant.quant_fdr);
    set_opt!("quant_dir", quant_dir, cfg.quant.dir);
    set!("mbr", mbr, cfg.quant.mbr);
    set!(
        "mbr_library_fdr",
        mbr_library_fdr,
        cfg.quant.mbr_library_fdr
    );
    set!("mbr_fdr", mbr_fdr, cfg.quant.mbr_fdr);
    set_parsed!(
        "quant_feature_fdr",
        quant_feature_fdr,
        cfg.quant.quant_feature_fdr.map(|v| v.text()),
        parse_quant_feature_fdr
    );

    // ── glyco ──
    set!("glyco", glyco, cfg.glyco.enabled);
    set_parsed!("glyco_core", glyco_core, cfg.glyco.core, parse_glyco_core);
    set!("glyco_tol_ppm", glyco_tol_ppm, cfg.glyco.tol_ppm);
    set!("glyco_max_peaks", glyco_max_peaks, cfg.glyco.max_peaks);
    set!("debug_glyco", debug_glyco, cfg.glyco.debug);
    set!("glyco_hcd_pair", glyco_hcd_pair, cfg.glyco.hcd_pair);
    // The glycan database is one choice made in one place: when any flag is
    // on the command line, the config's choice is ignored entirely, so the CLI
    // never ends up combined with a config value it overrides.
    if !cli_set(m, "glyco_glycan_gdb")
        && !cli_set(m, "glyco_glycan_wurcs")
        && !cli_set(m, "glyco_species")
    {
        let n_sources = usize::from(cfg.glyco.glycan_gdb.is_some())
            + usize::from(cfg.glyco.glycan_wurcs.is_some())
            + usize::from(cfg.glyco.species.is_some());
        if n_sources > 1 {
            return Err(
                "config: `glyco.glycan_gdb`, `glyco.glycan_wurcs` and `glyco.species` are \
                 mutually exclusive — set only one"
                    .to_string(),
            );
        }
        if let Some(p) = cfg.glyco.glycan_gdb {
            args.glyco_glycan_gdb = Some(p);
        }
        if let Some(p) = cfg.glyco.glycan_wurcs {
            args.glyco_glycan_wurcs = Some(p);
        }
        if let Some(s) = cfg.glyco.species {
            args.glyco_species = Some(
                parse_value_enum::<crate::cli::GlycoSpeciesFlag>(&s)
                    .map_err(|e| format!("config: `glyco.species`: {e}"))?,
            );
        }
    }
    set!("glyco_min_core_y", glyco_min_core_y, cfg.glyco.min_core_y);
    set_opt!(
        "glyco_retrieval_tol_ppm",
        glyco_retrieval_tol_ppm,
        cfg.glyco.retrieval_tol_ppm
    );
    set!(
        "glyco_peptide_first",
        glyco_peptide_first,
        cfg.glyco.peptide_first
    );
    set!(
        "glyco_index_sequon_only",
        glyco_index_sequon_only,
        cfg.glyco.index_sequon_only
    );
    set_opt!(
        "glyco_cz_max_charge",
        glyco_cz_max_charge,
        cfg.glyco.cz_max_charge
    );
    set!(
        "glyco_y_max_charge",
        glyco_y_max_charge,
        cfg.glyco.y_max_charge
    );

    // Validate the resolved numeric fields that the CLI enforces via `value_parser`
    // but the YAML path (plain assignment) would otherwise bypass. Re-checking the
    // final value is a no-op for CLI-supplied values (already validated by clap) and
    // rejects out-of-range values injected through `--config`.
    if let Some(v) = args.fdr {
        check_unit_fraction("fdr", v)?;
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
    check_unit_fraction("refine_select_psm_fdr", args.refine_select_psm_fdr)?;
    check_unit_fraction("mbr_library_fdr", args.mbr_library_fdr)?;
    check_unit_fraction("mbr_fdr", args.mbr_fdr)?;
    check_unit_fraction("tmt_min_purity", args.tmt_min_purity)?;
    check_unit_fraction("tmt_max_pep", args.tmt_max_pep)?;
    if !args.tmt_min_sn.is_finite() || args.tmt_min_sn < 0.0 {
        return Err(format!(
            "config: tmt_min_sn must be >= 0, got {}",
            args.tmt_min_sn
        ));
    }
    check_unit_fraction("lfq_min_cosine", args.lfq_min_cosine)?;
    check_unit_fraction("lfq_feature_fdr", args.lfq_feature_fdr)?;
    if let Some(v) = args.quant_fdr {
        check_unit_fraction("quant_fdr", v)?;
    }
    check_positive("lfq_rt_window", args.lfq_rt_window)?;
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
    // `--fragment-tol-ppm` / `--fragment-tol-da` are mutually exclusive on the CLI
    // (clap `conflicts_with`); a config that sets both, with neither on the
    // command line, is ambiguous.
    if args.fragment_tol_ppm.is_some() && args.fragment_tol_da.is_some() {
        return Err(
            "config: `scoring.fragment_tol_ppm` and `scoring.fragment_tol_da` are mutually \
             exclusive — set only one"
                .to_string(),
        );
    }

    Ok(())
}

/// Rules on the merged arguments (CLI, config and defaults together) that clap
/// cannot check because a value may come from the config file. Runs for every
/// search, with or without `--config`.
pub fn check_resolved(args: &SearchArgs) -> Result<(), String> {
    if args.glyco
        && args.glyco_glycan_gdb.is_none()
        && args.glyco_glycan_wurcs.is_none()
        && args.glyco_species.is_none()
    {
        return Err(
            "--glyco requires --glyco-glycan-gdb, --glyco-glycan-wurcs or --glyco-species \
             (config: `glyco.glycan_gdb`, `glyco.glycan_wurcs` or `glyco.species`): supply a \
             `.gdb` glycan database, a WURCS 2.0 CSV, or a bundled species"
                .to_string(),
        );
    }
    Ok(())
}

/// Parse a string into a clap `ValueEnum` (reusing the CLI's accepted names).
fn parse_value_enum<T: clap::ValueEnum>(s: &str) -> Result<T, String> {
    T::from_str(s, false).map_err(|e| e.to_string())
}

/// Same bound the CLI's `parse_unit_fraction` enforces: finite and within [0, 1].
fn check_unit_fraction(field: &str, v: f64) -> Result<(), String> {
    if !v.is_finite() || !(0.0..=1.0).contains(&v) {
        return Err(format!(
            "config: `{field}` must be finite and within [0, 1] (got {v})"
        ));
    }
    Ok(())
}

/// A tolerance string through the CLI's parser (finite and > 0), with the
/// error naming the config key.
fn parse_tol(field: &str, v: &str) -> Result<model::Tolerance, String> {
    crate::parse_precursor_tol(v).map_err(|e| format!("config: `{field}`: {e}"))
}

/// Same bound the CLI's `parse_positive_tol` enforces: finite and > 0.
fn check_positive(field: &str, v: f64) -> Result<(), String> {
    if !v.is_finite() || v <= 0.0 {
        return Err(format!(
            "config: `{field}` must be finite and > 0 (got {v})"
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unknown_key_is_rejected() {
        let y = "search:\n  precursor_tol: 20ppm\n  bogus_key: 1\n";
        assert!(
            RunConfig::from_yaml_str(y).is_err(),
            "unknown key must error"
        );
    }

    #[test]
    fn nested_sections_parse() {
        let y = "\
io:
  database: db.fasta
search:
  precursor_tol: 15ppm
  charge: \"2..4\"
glyco:
  enabled: true
  tol_ppm: 42.0
";
        let c = RunConfig::from_yaml_str(y).expect("valid config");
        assert_eq!(c.search.precursor_tol.as_deref(), Some("15ppm"));
        assert_eq!(c.search.charge.as_deref(), Some("2..4"));
        assert_eq!(c.glyco.enabled, Some(true));
        assert_eq!(c.glyco.tol_ppm, Some(42.0));
        assert_eq!(c.io.database, Some(PathBuf::from("db.fasta")));
    }

    #[test]
    fn empty_config_is_all_none() {
        let c = RunConfig::from_yaml_str("{}").expect("empty ok");
        assert!(c.search.precursor_tol.is_none());
        assert!(c.glyco.enabled.is_none());
    }

    #[test]
    fn unit_fraction_and_positive_bounds() {
        assert!(check_unit_fraction("fdr", 0.01).is_ok());
        assert!(check_unit_fraction("fdr", 0.0).is_ok());
        assert!(check_unit_fraction("fdr", 1.0).is_ok());
        assert!(
            check_unit_fraction("fdr", 5.0).is_err(),
            "fdr>1 must reject"
        );
        assert!(check_unit_fraction("fdr", -0.1).is_err());
        assert!(check_unit_fraction("fdr", f64::NAN).is_err());
        assert!(check_positive("frag", 20.0).is_ok());
        assert!(check_positive("frag", 0.0).is_err(), "tol=0 must reject");
        assert!(check_positive("frag", -1.0).is_err());
    }

    /// YAML glyco and chimeric numbers go through the CLI's bounds.
    #[test]
    fn yaml_glyco_tol_and_chimeric_kl_are_validated() {
        for bad in [
            "glyco:\n  tol_ppm: 0\n",
            "glyco:\n  tol_ppm: -5\n",
            "glyco:\n  tol_ppm: .nan\n",
            "glyco:\n  tol_ppm: .inf\n",
            "chimeric:\n  max_kl: -0.1\n",
            "chimeric:\n  max_kl: .nan\n",
        ] {
            let e = applied(bad, &[]).unwrap_err();
            assert!(
                e.contains("glyco_tol_ppm") || e.contains("chimeric_max_kl"),
                "{bad}: {e}"
            );
        }
        let a = applied("glyco:\n  tol_ppm: 30\nchimeric:\n  max_kl: 0\n", &[]).unwrap();
        assert_eq!(a.glyco_tol_ppm, 30.0);
        assert_eq!(a.chimeric_max_kl, 0.0);
    }

    /// The YAML path must reject values the CLI would reject via `value_parser`.
    #[test]
    fn apply_rejects_out_of_range_and_conflicting_yaml() {
        use clap::{CommandFactory, FromArgMatches};
        let matches = crate::TopCli::command()
            .try_get_matches_from([
                "andes",
                "--spectrum",
                "x",
                "--database",
                "y",
                "--output-pin",
                "o",
            ])
            .expect("parse baseline args");

        // out-of-range fdr from YAML → error
        let mut a = crate::SearchArgs::from_arg_matches(&matches).unwrap();
        let bad = RunConfig::from_yaml_str("rescoring:\n  fdr: 5.0\n").unwrap();
        let e = apply(bad, &mut a, &matches).unwrap_err();
        assert!(e.contains("fdr"), "expected fdr range error, got: {e}");

        // both fragment tolerances from YAML → mutual-exclusion error
        let mut a = crate::SearchArgs::from_arg_matches(&matches).unwrap();
        let both =
            RunConfig::from_yaml_str("scoring:\n  fragment_tol_ppm: 20\n  fragment_tol_da: 0.02\n")
                .unwrap();
        let e = apply(both, &mut a, &matches).unwrap_err();
        assert!(
            e.contains("mutually exclusive"),
            "expected conflict error, got: {e}"
        );

        // a valid config still applies cleanly
        let mut a = crate::SearchArgs::from_arg_matches(&matches).unwrap();
        let ok =
            RunConfig::from_yaml_str("rescoring:\n  fdr: 0.01\nscoring:\n  fragment_tol_ppm: 20\n")
                .unwrap();
        assert!(apply(ok, &mut a, &matches).is_ok());
    }

    /// Out-of-range quant values from YAML fail with an error naming the key,
    /// as the CLI's value parsers would.
    #[test]
    fn apply_rejects_out_of_range_quant_keys() {
        use clap::{CommandFactory, FromArgMatches};
        let matches = crate::TopCli::command()
            .try_get_matches_from([
                "andes",
                "--spectrum",
                "x",
                "--database",
                "y",
                "--output-pin",
                "o",
            ])
            .expect("parse baseline args");
        for (key, value) in [
            ("tmt_min_purity", "1.5"),
            ("tmt_min_purity", "-0.1"),
            ("lfq_min_cosine", "2.0"),
            ("lfq_feature_fdr", "-0.5"),
            ("quant_fdr", "5.0"),
            ("quant_fdr", ".nan"),
            ("quant_feature_fdr", "1.5"),
            ("quant_feature_fdr", "0"),
            ("quant_feature_fdr", "none"),
            ("lfq_rt_window", "0"),
            ("lfq_rt_window", "-30"),
            ("lfq_rt_window", ".inf"),
            ("lfq_tol", "0ppm"),
            ("lfq_tol", "-5ppm"),
            ("tmt_tol", "0da"),
            ("tmt_tol", "-20ppm"),
        ] {
            let mut a = crate::SearchArgs::from_arg_matches(&matches).unwrap();
            let yaml = if key.ends_with("_tol") {
                format!("quant:\n  {key}: \"{value}\"\n")
            } else {
                format!("quant:\n  {key}: {value}\n")
            };
            let cfg = RunConfig::from_yaml_str(&yaml).unwrap();
            let Err(e) = apply(cfg, &mut a, &matches) else {
                panic!("{key}: {value} must be rejected");
            };
            assert!(
                e.contains(key),
                "{key}: {value}: error does not name the key: {e}"
            );
        }

        // In-range values still apply.
        let mut a = crate::SearchArgs::from_arg_matches(&matches).unwrap();
        let ok = RunConfig::from_yaml_str(
            "quant:\n  tmt_min_purity: 0.5\n  lfq_min_cosine: 0.8\n  lfq_feature_fdr: 0.05\n  \
             quant_fdr: 0.01\n  lfq_rt_window: 45\n  lfq_tol: 10ppm\n  tmt_tol: 0.003da\n",
        )
        .unwrap();
        assert!(apply(ok, &mut a, &matches).is_ok());
        assert_eq!(a.quant_feature_fdr.0, None, "the default is off");
        for (value, want) in [("off", None), ("0.01", Some(0.01)), ("0.05", Some(0.05))] {
            let mut a = crate::SearchArgs::from_arg_matches(&matches).unwrap();
            let cfg = RunConfig::from_yaml_str(&format!("quant:\n  quant_feature_fdr: {value}\n"))
                .unwrap();
            apply(cfg, &mut a, &matches).unwrap();
            assert_eq!(a.quant_feature_fdr.0, want, "{value}");
        }
    }

    /// The UniProt XML keys apply through the CLI's validator; an explicit CLI
    /// flag still wins.
    #[test]
    fn uniprot_keys_apply_with_cli_precedence() {
        use clap::{CommandFactory, FromArgMatches};
        let base = [
            "andes",
            "--spectrum",
            "x",
            "--database",
            "y.xml.gz",
            "--output-pin",
            "o",
        ];
        let matches = crate::TopCli::command()
            .try_get_matches_from(base)
            .expect("parse baseline args");
        let mut a = crate::SearchArgs::from_arg_matches(&matches).unwrap();
        assert!(!a.uniprot_variants);
        assert_eq!(a.uniprot_evidence, input::EvidenceLevel::Experimental);
        assert!(!a.uniprot_isoforms && !a.uniprot_chains);
        let cfg = RunConfig::from_yaml_str(
            "io:\n  uniprot_variants: true\n  uniprot_isoforms: true\n  uniprot_chains: true\n  \
             uniprot_evidence: curated\n",
        )
        .unwrap();
        apply(cfg, &mut a, &matches).unwrap();
        assert!(a.uniprot_variants && a.uniprot_isoforms && a.uniprot_chains);
        assert_eq!(a.uniprot_evidence, input::EvidenceLevel::Curated);

        let bad = RunConfig::from_yaml_str("io:\n  uniprot_evidence: strong\n").unwrap();
        let mut a = crate::SearchArgs::from_arg_matches(&matches).unwrap();
        assert!(apply(bad, &mut a, &matches).is_err());

        let cli = crate::TopCli::command()
            .try_get_matches_from(base.into_iter().chain(["--uniprot-evidence", "any"]))
            .unwrap();
        let mut a = crate::SearchArgs::from_arg_matches(&cli).unwrap();
        let cfg = RunConfig::from_yaml_str("io:\n  uniprot_evidence: curated\n").unwrap();
        apply(cfg, &mut a, &cli).unwrap();
        assert_eq!(a.uniprot_evidence, input::EvidenceLevel::Any);
    }

    /// The predictor keys apply through the CLI's validators; an explicit
    /// `--predict` still wins.
    #[test]
    fn predict_keys_apply_with_cli_precedence() {
        use crate::cli::{PredictFeature, PredictorBackend};
        use clap::{CommandFactory, FromArgMatches};
        let base = [
            "andes",
            "--spectrum",
            "x",
            "--database",
            "y.fasta",
            "--output-pin",
            "o",
        ];
        let matches = crate::TopCli::command()
            .try_get_matches_from(base)
            .expect("parse baseline args");
        let mut a = crate::SearchArgs::from_arg_matches(&matches).unwrap();
        let cfg = RunConfig::from_yaml_str(
            "rescoring:\n  predict: [rt, ms2]\n  predictor_models: /m\n  \
             predictor_backend: ort\n  predictor_nce: 27\n",
        )
        .unwrap();
        apply(cfg, &mut a, &matches).unwrap();
        assert_eq!(a.predict, vec![PredictFeature::Rt, PredictFeature::Ms2]);
        assert_eq!(a.predictor_models, Some(PathBuf::from("/m")));
        assert_eq!(a.predictor_backend, PredictorBackend::Ort);
        assert_eq!(a.predictor_nce, Some(27.0));

        let bad = RunConfig::from_yaml_str("rescoring:\n  predict: [ms3]\n").unwrap();
        let mut a = crate::SearchArgs::from_arg_matches(&matches).unwrap();
        assert!(apply(bad, &mut a, &matches).is_err());

        let cli = crate::TopCli::command()
            .try_get_matches_from(base.into_iter().chain(["--predict", "none"]))
            .unwrap();
        let mut a = crate::SearchArgs::from_arg_matches(&cli).unwrap();
        let cfg = RunConfig::from_yaml_str("rescoring:\n  predict: [rt]\n").unwrap();
        apply(cfg, &mut a, &cli).unwrap();
        assert_eq!(a.predict, vec![PredictFeature::None]);
    }

    /// `search.index_cache_dir` applies, and an explicit `--index-cache-dir` wins.
    #[test]
    fn index_cache_dir_key_has_cli_precedence() {
        let yaml = "search:\n  index_cache_dir: from_config\n";
        let a = applied(yaml, &[]).unwrap();
        assert_eq!(a.index_cache_dir, Some(PathBuf::from("from_config")));
        let a = applied(yaml, &["--index-cache-dir", "from_cli"]).unwrap();
        assert_eq!(a.index_cache_dir, Some(PathBuf::from("from_cli")));
        assert_eq!(applied("{}", &[]).unwrap().index_cache_dir, None);

        let yaml = "io:\n  temp_dir: scratch\n";
        assert_eq!(
            applied(yaml, &[]).unwrap().temp_dir,
            Some(PathBuf::from("scratch"))
        );
        assert_eq!(
            applied(yaml, &["--temp-dir", "cli"]).unwrap().temp_dir,
            Some(PathBuf::from("cli"))
        );
    }

    /// The shipped `config.example.yaml` must parse against the current schema, so
    /// the documented template can never drift from the accepted keys.
    #[test]
    fn shipped_example_config_parses() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../config.example.yaml");
        let s = std::fs::read_to_string(path).expect("read config.example.yaml");
        RunConfig::from_yaml_str(&s).expect("config.example.yaml must parse");
    }

    fn baseline_matches(extra: &[&str]) -> clap::ArgMatches {
        use clap::CommandFactory;
        let base = [
            "andes",
            "--spectrum",
            "x",
            "--database",
            "y",
            "--output-pin",
            "o",
        ];
        crate::TopCli::command()
            .try_get_matches_from(base.iter().chain(extra))
            .expect("parse args")
    }

    fn applied(yaml: &str, extra: &[&str]) -> Result<crate::SearchArgs, String> {
        use clap::FromArgMatches;
        let m = baseline_matches(extra);
        let mut a = crate::SearchArgs::from_arg_matches(&m).unwrap();
        apply(RunConfig::from_yaml_str(yaml)?, &mut a, &m)?;
        Ok(a)
    }

    /// Flags deliberately left out of the config file (see the module docs).
    const CLI_ONLY: &[&str] = &[
        "config",
        "cal_min_spec_keys",
        "precursor_mono_dump",
        "glyco_scans",
    ];

    /// Every search flag maps to a config key unless it is listed as CLI-only, so
    /// a new flag cannot silently miss the config file.
    #[test]
    fn every_search_flag_has_a_config_key() {
        use clap::CommandFactory;
        // `apply` reads a flag's config key through one of these macros, or tests
        // `cli_set` itself for keys that need custom handling. A checkout with
        // CRLF line endings (Windows) is normalised before joining the lines.
        let source = include_str!("config.rs")
            .replace("\r\n", "\n")
            .replace("(\n        ", "(");
        let mapped = |id: &str| {
            ["set!", "set_opt!", "set_parsed!"]
                .iter()
                .any(|mac| source.contains(&format!("{mac}(\"{id}\",")))
                || source.contains(&format!("cli_set(m, \"{id}\")"))
        };
        let missing: Vec<String> = crate::TopCli::command()
            .get_arguments()
            .map(|a| a.get_id().to_string())
            .filter(|id| !matches!(id.as_str(), "help" | "version"))
            .filter(|id| !CLI_ONLY.contains(&id.as_str()) && !mapped(id))
            .collect();
        assert!(missing.is_empty(), "flags with no config key: {missing:?}");
    }

    #[test]
    fn new_keys_apply_and_validate() {
        let a = applied(
            "search:\n  fragment_index: \"off\"\n  fragment_index_slice_da: 40\n  index_cache_dir: shared/cache\n\
             scoring:\n  gbdt_max_trees: 0\n  precursor_mono: \"off\"\n  ethcd_activation: etd\n  \
             peak_filter: \"100:20\"\n\
             quant:\n  lfq_gaussian_width: 12.5\n\
             glyco:\n  enabled: true\n  species: mouse-large\n  min_core_y: 2\n  \
             retrieval_tol_ppm: 20\n  peptide_first: true\n  index_sequon_only: true\n  \
             cz_max_charge: 2\n  y_max_charge: 4\n",
            &[],
        )
        .unwrap();
        assert_eq!(a.fragment_index, crate::cli::FragmentIndexFlag::Off);
        assert_eq!(a.fragment_index_slice_da, Some(40.0));
        assert_eq!(a.index_cache_dir, Some(PathBuf::from("shared/cache")));
        assert_eq!(a.gbdt_max_trees, Some(0));
        assert_eq!(a.precursor_mono, crate::cli::PrecursorMonoFlag::Off);
        assert_eq!(a.ethcd_activation, crate::cli::EthcdActivationFlag::Etd);
        assert_eq!(a.peak_filter.as_deref(), Some("100:20"));
        assert_eq!(a.lfq_gaussian_width, Some(12.5));
        assert_eq!(
            a.glyco_species,
            Some(crate::cli::GlycoSpeciesFlag::MouseLarge)
        );
        assert_eq!(a.glyco_min_core_y, 2);
        assert_eq!(a.glyco_retrieval_tol_ppm, Some(20.0));
        assert!(a.glyco_peptide_first && a.glyco_index_sequon_only);
        assert_eq!(a.glyco_cz_max_charge, Some(2));
        assert_eq!(a.glyco_y_max_charge, 4);
        assert!(check_resolved(&a).is_ok());

        // Defaults stay unset, so the search resolves them by mode.
        let d = applied("{}", &[]).unwrap();
        assert_eq!(d.gbdt_max_trees, None);

        for bad in [
            "search:\n  fragment_index: sometimes\n",
            "search:\n  fragment_index_slice_da: 0\n",
            "scoring:\n  precursor_mono: maybe\n",
            "quant:\n  lfq_gaussian_width: -1\n",
            "glyco:\n  retrieval_tol_ppm: 0\n",
            "glyco:\n  species: martian\n",
        ] {
            assert!(applied(bad, &[]).is_err(), "{bad} must be rejected");
        }

        // An explicit CLI flag still wins.
        let a = applied(
            "scoring:\n  gbdt_max_trees: 0\n",
            &["--gbdt-max-trees", "50"],
        )
        .unwrap();
        assert_eq!(a.gbdt_max_trees, Some(50));
    }

    /// A rescoring backend on the command line wins over the config's choice.
    #[test]
    fn cli_rescore_backend_wins_over_yaml() {
        let a = applied("rescoring:\n  enabled: true\n", &["--rescore-native"]).unwrap();
        assert!(
            a.rescore_native && !a.rescore,
            "CLI --rescore-native must win"
        );
        let a = applied("rescoring:\n  native: true\n", &["--rescore"]).unwrap();
        assert!(a.rescore && !a.rescore_native);
        let a = applied("rescoring:\n  enabled: true\n", &[]).unwrap();
        assert!(a.rescore && !a.rescore_native);
        let a = applied("rescoring:\n  native: true\n", &[]).unwrap();
        assert!(a.rescore_native && !a.rescore);
    }

    /// A fragment tolerance on the command line replaces both config keys, so a
    /// CLI ppm value and a config Da value never conflict.
    #[test]
    fn cli_fragment_tol_wins_over_yaml_without_conflict() {
        let a = applied(
            "scoring:\n  fragment_tol_da: 0.5\n",
            &["--fragment-tol-ppm", "20"],
        )
        .unwrap();
        assert_eq!(a.fragment_tol_ppm, Some(20.0));
        assert_eq!(a.fragment_tol_da, None);
        let a = applied(
            "scoring:\n  fragment_tol_ppm: 20\n",
            &["--fragment-tol-da", "0.5"],
        )
        .unwrap();
        assert_eq!(a.fragment_tol_da, Some(0.5));
        assert_eq!(a.fragment_tol_ppm, None);
        let a = applied("scoring:\n  fragment_tol_da: 0.5\n", &[]).unwrap();
        assert_eq!(a.fragment_tol_da, Some(0.5));
    }

    /// `--glyco` needs a glycan database whichever way it was switched on.
    #[test]
    fn glyco_requires_a_glycan_database_after_the_merge() {
        let a = applied("glyco:\n  enabled: true\n", &[]).unwrap();
        let e = check_resolved(&a).unwrap_err();
        assert!(e.contains("glyco.species"), "{e}");
        let a = applied("{}", &["--glyco"]).unwrap();
        assert!(check_resolved(&a).is_err());
        let a = applied("glyco:\n  glycan_gdb: g.gdb\n", &["--glyco"]).unwrap();
        assert!(check_resolved(&a).is_ok());
        let a = applied("glyco:\n  enabled: true\n", &["--glyco-species", "human"]).unwrap();
        assert!(check_resolved(&a).is_ok());
        assert!(check_resolved(&applied("{}", &[]).unwrap()).is_ok());

        // The WURCS CSV source satisfies the requirement like the other two.
        let a = applied(
            "glyco:\n  glycan_wurcs: w.csv\n",
            &["--glyco"],
        )
        .unwrap();
        assert_eq!(a.glyco_glycan_wurcs.as_deref(), Some(std::path::Path::new("w.csv")));
        assert!(check_resolved(&a).is_ok());
        let a = applied("{}", &["--glyco-glycan-wurcs", "w.csv"]).unwrap();
        assert_eq!(a.glyco_glycan_wurcs.as_deref(), Some(std::path::Path::new("w.csv")));

        // Two sources in the config is ambiguous.
        let both = "glyco:\n  enabled: true\n  glycan_gdb: g.gdb\n  species: human\n";
        assert!(applied(both, &[])
            .unwrap_err()
            .contains("mutually exclusive"));
        let both = "glyco:\n  enabled: true\n  glycan_gdb: g.gdb\n  glycan_wurcs: w.csv\n";
        assert!(applied(both, &[])
            .unwrap_err()
            .contains("mutually exclusive"));
        // A source on the command line replaces the config's choice entirely.
        let a = applied(
            "glyco:\n  glycan_gdb: g.gdb\n",
            &["--glyco-species", "mouse"],
        )
        .unwrap();
        assert_eq!(a.glyco_glycan_gdb, None);
        assert_eq!(a.glyco_species, Some(crate::cli::GlycoSpeciesFlag::Mouse));
    }
}
