use std::path::PathBuf;

use serde::Deserialize;

use crate::SearchArgs;

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
    pub sdrf: Option<PathBuf>,
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
    pub fixed_mods: Option<Vec<serde_yaml::Value>>,
    pub variable_mods: Option<Vec<serde_yaml::Value>>,
    pub xl: Option<serde_yaml::Value>,
    pub localize: Option<String>,
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
    pub precursor_mono: Option<String>,
    pub ethcd_activation: Option<String>,
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
}

#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct RescoringCfg {
    pub enabled: Option<bool>,
    pub native: Option<bool>,
    pub psm_fdr: Option<f64>,
    pub proteoform_fdr: Option<f64>,
    pub protein_fdr: Option<f64>,
    pub pep: Option<f64>,
    pub percolator_bin: Option<PathBuf>,
    pub percolator_docker: Option<bool>,
    pub percolator_image: Option<String>,
    pub percolator_args: Option<String>,
    pub keep_pin: Option<bool>,
    pub memory_gb: Option<u64>,
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
    pub tmt_unmix: Option<String>,
    pub lfq: Option<bool>,
    pub lfq_tol: Option<String>,
    pub lfq_rt_window: Option<f64>,
    pub lfq_gaussian_width: Option<f64>,
    pub lfq_min_cosine: Option<f64>,
    pub label: Option<String>,
    pub label_reference: Option<usize>,
    pub label_arg_pro: Option<bool>,
    pub label_2h_enrichment_silac: Option<f64>,
    pub label_2h_enrichment_chemical: Option<f64>,
    pub label_heavy_enrichment: Option<f64>,
    pub quant_fdr: Option<f64>,
    pub quant_transfer_fdr: Option<f64>,
    pub quant_channel_fdr: Option<f64>,
    pub dir: Option<PathBuf>,
    pub mbr: Option<bool>,
    pub mbr_donor_fdr: Option<f64>,
    pub mbr_wrong_log2: Option<f64>,
}

fn parse_glyco_core(s: &str) -> Result<crate::cli::GlycoCoreFlag, String> {
    crate::cli::parse_enum("glyco.core", s, true)
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
    pub glycan_gdb: Option<PathBuf>,
    pub glycan_wurcs: Option<PathBuf>,
    pub species: Option<String>,
    pub min_core_y: Option<u32>,
    pub retrieval_tol_ppm: Option<f64>,
    pub peptide_first: Option<bool>,
    pub index_sequon_only: Option<bool>,
    pub cz_max_charge: Option<u8>,
    pub y_max_charge: Option<u8>,
}

impl RunConfig {
    pub fn from_yaml_str(s: &str) -> Result<Self, String> {
        serde_yaml::from_str(s).map_err(|e| format!("--config parse error: {e}"))
    }

    pub fn load(path: &std::path::Path) -> Result<Self, String> {
        let s = std::fs::read_to_string(path)
            .map_err(|e| format!("--config: cannot read {}: {e}", path.display()))?;
        Self::from_yaml_str(&s)
    }
}

pub(crate) fn cli_set(m: &clap::ArgMatches, id: &str) -> bool {
    matches!(
        m.value_source(id),
        Some(clap::parser::ValueSource::CommandLine)
    )
}

pub fn apply_db(
    cfg: &RunConfig,
    db: &mut crate::cli::DatabaseArgs,
    m: &clap::ArgMatches,
) -> Result<(), String> {
    macro_rules! set {
        ($id:literal, $field:ident, $val:expr) => {
            if let Some(v) = &$val {
                if !cli_set(m, $id) {
                    db.$field = v.clone();
                }
            }
        };
    }
    macro_rules! set_opt {
        ($id:literal, $field:ident, $val:expr) => {
            if let Some(v) = &$val {
                if !cli_set(m, $id) {
                    db.$field = Some(v.clone());
                }
            }
        };
    }
    macro_rules! set_parsed {
        ($id:literal, $field:ident, $val:expr, $parse:expr) => {
            if let Some(s) = &$val {
                if !cli_set(m, $id) {
                    db.$field = ($parse)(s)?;
                }
            }
        };
    }
    set!("uniprot_variants", uniprot_variants, cfg.io.uniprot_variants);
    set!("uniprot_isoforms", uniprot_isoforms, cfg.io.uniprot_isoforms);
    set!("uniprot_chains", uniprot_chains, cfg.io.uniprot_chains);
    set!("annotations", annotations, cfg.io.annotations);
    set!("knockoffs", knockoffs, cfg.io.knockoffs);
    set_parsed!(
        "uniprot_evidence",
        uniprot_evidence,
        cfg.io.uniprot_evidence,
        crate::parse_uniprot_evidence
    );
    set_parsed!(
        "enzyme",
        enzyme,
        cfg.search.enzyme,
        crate::cli::parse_enzyme_list
    );
    if let Some(items) = &cfg.search.fixed_mods {
        if !cli_set(m, "fixed_mods") {
            db.fixed_mods = Some(mod_list("search.fixed_mods", items, true)?);
        }
    }
    if let Some(items) = &cfg.search.variable_mods {
        if !cli_set(m, "variable_mods") {
            db.variable_mods = Some(mod_list("search.variable_mods", items, false)?);
        }
    }
    set!("decoy_prefix", decoy_prefix, cfg.decoys.prefix);
    set_opt!("decoy_suffix", decoy_suffix, cfg.decoys.suffix);
    Ok(())
}

pub fn apply_protein_config(
    p: &mut crate::cli::ProteinArgs,
    m: &clap::ArgMatches,
) -> Result<(), String> {
    if let Some(path) = p.config.clone() {
        apply_db(&RunConfig::load(&path)?, &mut p.db, m)?;
    }
    Ok(())
}

pub fn mod_list(
    key: &str,
    items: &[serde_yaml::Value],
    fixed: bool,
) -> Result<crate::cli::ModList, String> {
    use model::mod_spec::ModDecl;
    let mut mods = Vec::new();
    for (i, item) in items.iter().enumerate() {
        let at = |e: String| format!("{key}[{i}]: {e}");
        let decl = match item {
            serde_yaml::Value::String(s) => ModDecl::parse(s).map_err(|e| at(e.to_string()))?,
            serde_yaml::Value::Mapping(map) => mod_decl(map).map_err(at)?,
            _ => {
                return Err(at(
                    "expected `Name (sites)` or a map with `name` or `mass` and `sites`".into(),
                ))
            }
        };
        mods.extend(
            decl.resolve(fixed, model::unimod::bundled())
                .map_err(|e| at(e.to_string()))?,
        );
    }
    Ok(crate::cli::ModList(mods))
}

fn mod_decl(map: &serde_yaml::Mapping) -> Result<model::mod_spec::ModDecl, String> {
    use serde_yaml::Value;
    let text = |key: &str, v: &Value| -> Result<String, String> {
        match v {
            Value::String(s) => Ok(s.clone()),
            _ => Err(format!("`{key}` must be a string")),
        }
    };
    let numbers = |key: &str, v: &Value| -> Result<Vec<f64>, String> {
        let one = |v: &Value| v.as_f64().ok_or_else(|| format!("`{key}` takes numbers"));
        match v {
            Value::Sequence(seq) => seq.iter().map(one).collect(),
            v => Ok(vec![one(v)?]),
        }
    };
    let mut d = model::mod_spec::ModDecl::default();
    let mut sites = None;
    for (k, v) in map {
        let key = k.as_str().ok_or("modification keys are strings")?;
        match key {
            "name" => d.name = Some(text(key, v)?),
            "mass" => {
                d.mass = Some(
                    v.as_f64()
                        .ok_or("`mass` is a number (the mass shift in Da)")?,
                )
            }
            "sites" => sites = Some(text(key, v)?),
            "neutral_loss" => d.neutral_loss = numbers(key, v)?,
            "loss_class" => d.loss_class = Some(text(key, v)?),
            "restrict_to" => d.restrict_to = Some(text(key, v)?),
            "cleavage" => d.cleavage = Some(text(key, v)?),
            "diagnostic_ions" => d.diagnostic_ions = numbers(key, v)?,
            other => {
                return Err(format!(
                    "unknown key `{other}` (expected name, mass, sites, neutral_loss, \
                     loss_class, restrict_to, cleavage, diagnostic_ions)"
                ))
            }
        }
    }
    match (sites, d.name.clone(), d.mass) {
        (Some(s), _, _) => d.sites = s,
        (None, Some(name), None) => {
            let parsed = model::mod_spec::ModDecl::parse(&name).map_err(|e| e.to_string())?;
            d.name = parsed.name;
            d.sites = parsed.sites;
        }
        (None, _, _) => return Err("`sites` is missing".into()),
    }
    Ok(d)
}

pub fn xl_linker(
    key: &str,
    v: &serde_yaml::Value,
) -> Result<andes_xl::crosslinker::Crosslinker, String> {
    use andes_xl::spec::XlDecl;
    use serde_yaml::Value;
    let at = |e: String| format!("{key}: {e}");
    let db = model::unimod::bundled();
    match v {
        Value::String(s) => XlDecl::parse(s)
            .and_then(|d| d.resolve(db))
            .map_err(|e| at(e.to_string())),
        Value::Mapping(map) if map.contains_key("intact") => {
            let r: model::xl_reagents::Reagent =
                serde_yaml::from_value(v.clone()).map_err(|e| at(e.to_string()))?;
            andes_xl::spec::custom(&r, db).map_err(|e| at(e.to_string()))
        }
        Value::Mapping(map) => {
            let mut name = None;
            let mut sites = None;
            for (k, val) in map {
                let text = || {
                    val.as_str().map(str::to_string).ok_or_else(|| {
                        at(format!("`{}` must be a string", k.as_str().unwrap_or("?")))
                    })
                };
                match k.as_str() {
                    Some("name") => name = Some(text()?),
                    Some("sites") => sites = Some(text()?),
                    other => {
                        return Err(at(format!(
                            "unknown key `{}` (a named reagent takes name and sites; a custom \
                             reagent also takes intact, stubs, doublets, monolinks, aliases, \
                             xlmod, source)",
                            other.unwrap_or("?")
                        )))
                    }
                }
            }
            let name = name.ok_or_else(|| at("`name` is missing".into()))?;
            let mut d = XlDecl::parse(&name).map_err(|e| at(e.to_string()))?;
            if sites.is_some() {
                d.sites = sites;
            }
            d.resolve(db).map_err(|e| at(e.to_string()))
        }
        _ => Err(at("expected `Name (sites)` or a map with `name`".into())),
    }
}

pub fn apply(cfg: RunConfig, args: &mut SearchArgs, m: &clap::ArgMatches) -> Result<(), String> {
    macro_rules! set {
        ($id:literal, $($field:ident).+, $val:expr) => {
            if let Some(v) = $val {
                if !cli_set(m, $id) {
                    args.$($field).+ = v;
                }
            }
        };
    }
    macro_rules! set_opt {
        ($id:literal, $($field:ident).+, $val:expr) => {
            if let Some(v) = $val {
                if !cli_set(m, $id) {
                    args.$($field).+ = Some(v);
                }
            }
        };
    }
    macro_rules! set_parsed {
        ($id:literal, $($field:ident).+, $val:expr, $parse:expr) => {
            if let Some(s) = $val {
                if !cli_set(m, $id) {
                    args.$($field).+ = ($parse)(&s)?;
                }
            }
        };
    }

    apply_db(&cfg, &mut args.db, m)?;

    set_opt!("sdrf", sdrf, cfg.io.sdrf);
    set!("spectrum", spectrum, cfg.io.spectrum);
    set_opt!("database", database, cfg.io.database);
    set_opt!("output_pin", output_pin, cfg.io.output_pin);
    set_opt!("output_tsv", output_tsv, cfg.io.output_tsv);
    set_opt!("output_parquet", output_parquet, cfg.io.output_parquet);
    set_opt!("temp_dir", temp_dir, cfg.io.temp_dir);

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
    if let Some(s) = cfg.search.isotope_error {
        if !cli_set(m, "isotope_error") {
            args.isotope_error = Some(crate::parse_isotope_error_range(&s)?);
        }
    }
    set_parsed!(
        "enzyme_specificity",
        enzyme_specificity,
        cfg.search.enzyme_specificity,
        |s: &str| crate::cli::parse_enum("enzyme_specificity", s, true)
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
    if let Some(v) = cfg.search.localize.filter(|_| !cli_set(m, "localize")) {
        if v != "ipf" {
            return Err(format!("search.localize: '{v}' is not one of: ipf"));
        }
        args.localize = Some(v);
    }
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
        |s: &str| crate::cli::parse_enum("fragmentation", s, true)
    );
    set_parsed!(
        "protocol",
        protocol,
        cfg.scoring.protocol,
        |s: &str| crate::cli::parse_enum("protocol", s, true)
    );
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

    set_parsed!(
        "decoy_strategy",
        decoy_strategy,
        cfg.decoys.strategy,
        crate::cli::parse_decoy_strategy
    );
    set!("decoy_seed", decoy_seed, cfg.decoys.seed);

    set!("chimeric", chimeric, cfg.chimeric.enabled);
    set!(
        "chimeric_max_coisolated",
        chimeric_max_coisolated,
        cfg.chimeric.max_coisolated
    );
    set!("chimeric_max_kl", chimeric_max_kl, cfg.chimeric.max_kl);

    set!("refine", refine, cfg.refine.enabled);
    set_opt!("refine_config", refine_config, cfg.refine.config);

    if !cli_set(m, "rescore") && !cli_set(m, "rescore_native") {
        set!("rescore", rescore, cfg.rescoring.enabled);
        set!("rescore_native", rescore_native, cfg.rescoring.native);
    }
    set_opt!("psm_fdr", psm_fdr, cfg.rescoring.psm_fdr);
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
    set!(
        "rescore_memory_gb",
        rescore_memory_gb,
        cfg.rescoring.memory_gb
    );
    if let Some(list) = cfg.rescoring.predict {
        if !cli_set(m, "predict") {
            args.predict = list
                .iter()
                .map(|s| crate::cli::parse_enum("predict entry", s, true))
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
        |s: &str| crate::cli::parse_enum("predictor_backend", s, true)
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
    set_parsed!(
        "tmt_unmix",
        tmt_unmix,
        cfg.quant.tmt_unmix,
        parse_value_enum::<crate::cli::TmtUnmix>
    );
    set!("lfq", lfq, cfg.quant.lfq);
    if let Some(v) = cfg.quant.lfq_tol {
        if !cli_set(m, "lfq_tol") {
            args.lfq_tol = Some(parse_tol("lfq_tol", &v)?);
        }
    }
    set_opt!("lfq_rt_window", lfq_rt_window, cfg.quant.lfq_rt_window);
    set_opt!(
        "lfq_gaussian_width",
        lfq_gaussian_width,
        cfg.quant.lfq_gaussian_width
    );
    set!("lfq_min_cosine", lfq_min_cosine, cfg.quant.lfq_min_cosine);
    set_opt!("label", label, cfg.quant.label);
    set!(
        "label_reference",
        label_reference,
        cfg.quant.label_reference
    );
    set!("label_arg_pro", label_arg_pro, cfg.quant.label_arg_pro);
    set!(
        "label_2h_enrichment_silac",
        label_2h_enrichment_silac,
        cfg.quant.label_2h_enrichment_silac
    );
    set!(
        "label_2h_enrichment_chemical",
        label_2h_enrichment_chemical,
        cfg.quant.label_2h_enrichment_chemical
    );
    set!(
        "label_heavy_enrichment",
        label_heavy_enrichment,
        cfg.quant.label_heavy_enrichment
    );
    set_opt!("quant_fdr", quant_fdr, cfg.quant.quant_fdr);
    set_opt!(
        "quant_transfer_fdr",
        quant_transfer_fdr,
        cfg.quant.quant_transfer_fdr
    );
    set_opt!(
        "quant_channel_fdr",
        quant_channel_fdr,
        cfg.quant.quant_channel_fdr
    );
    set_opt!("quant_dir", quant_dir, cfg.quant.dir);
    set!("mbr", mbr, cfg.quant.mbr);
    set_opt!("mbr_donor_fdr", mbr_donor_fdr, cfg.quant.mbr_donor_fdr);
    set!("mbr_wrong_log2", mbr_wrong_log2, cfg.quant.mbr_wrong_log2);

    set!("glyco", glyco, cfg.glyco.enabled);
    if let Some(v) = &cfg.search.xl {
        if !cli_set(m, "xl") {
            args.xl = Some(xl_linker("search.xl", v)?);
        }
    }
    set_parsed!("glyco_core", glyco_core, cfg.glyco.core, parse_glyco_core);
    set!("glyco_tol_ppm", glyco_tol_ppm, cfg.glyco.tol_ppm);
    set!("glyco_max_peaks", glyco_max_peaks, cfg.glyco.max_peaks);
    set!("debug_glyco", debug_glyco, cfg.glyco.debug);
    set!("glyco_hcd_pair", glyco_hcd_pair, cfg.glyco.hcd_pair);
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

    if let Some(v) = args.psm_fdr {
        check_unit_fraction("psm_fdr", v)?;
    }
    if args.rescore_memory_gb == 0 {
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
    check_unit_fraction("tmt_max_pep", args.tmt_max_pep)?;
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

fn parse_value_enum<T: clap::ValueEnum>(s: &str) -> Result<T, String> {
    T::from_str(s, false).map_err(|e| e.to_string())
}

fn check_unit_fraction(field: &str, v: f64) -> Result<(), String> {
    if !v.is_finite() || !(0.0..=1.0).contains(&v) {
        return Err(format!(
            "config: `{field}` must be finite and within [0, 1] (got {v})"
        ));
    }
    Ok(())
}

fn parse_tol(field: &str, v: &str) -> Result<model::Tolerance, String> {
    crate::parse_precursor_tol(v).map_err(|e| format!("config: `{field}`: {e}"))
}

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

        let mut a = crate::SearchArgs::from_arg_matches(&matches).unwrap();
        let bad = RunConfig::from_yaml_str("rescoring:\n  psm_fdr: 5.0\n").unwrap();
        let e = apply(bad, &mut a, &matches).unwrap_err();
        assert!(
            e.contains("psm_fdr"),
            "expected psm_fdr range error, got: {e}"
        );

        let mut a = crate::SearchArgs::from_arg_matches(&matches).unwrap();
        let both =
            RunConfig::from_yaml_str("scoring:\n  fragment_tol_ppm: 20\n  fragment_tol_da: 0.02\n")
                .unwrap();
        let e = apply(both, &mut a, &matches).unwrap_err();
        assert!(
            e.contains("mutually exclusive"),
            "expected conflict error, got: {e}"
        );

        let mut a = crate::SearchArgs::from_arg_matches(&matches).unwrap();
        let ok = RunConfig::from_yaml_str(
            "rescoring:\n  psm_fdr: 0.01\nscoring:\n  fragment_tol_ppm: 20\n",
        )
        .unwrap();
        assert!(apply(ok, &mut a, &matches).is_ok());
    }

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
            ("quant_fdr", "5.0"),
            ("quant_fdr", ".nan"),
            ("quant_transfer_fdr", "-0.5"),
            ("quant_channel_fdr", "1.5"),
            ("mbr_donor_fdr", "2"),
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

        let mut a = crate::SearchArgs::from_arg_matches(&matches).unwrap();
        let ok = RunConfig::from_yaml_str(
            "quant:\n  tmt_min_purity: 0.5\n  lfq_min_cosine: 0.8\n  quant_channel_fdr: 0.05\n  \
             quant_fdr: 0.01\n  lfq_rt_window: 45\n  lfq_tol: 10ppm\n  tmt_tol: 0.003da\n",
        )
        .unwrap();
        assert!(apply(ok, &mut a, &matches).is_ok());
    }

    #[test]
    fn fdr_keys_apply_and_inherit() {
        let a = applied(
            "rescoring:\n  psm_fdr: 0.02\n  proteoform_fdr: 0.06\n  protein_fdr: 0.03\n\
             quant:\n  quant_fdr: 0.04\n  mbr_donor_fdr: 0.005\n",
            &[],
        )
        .unwrap();
        assert_eq!(a.psm_fdr(), 0.02);
        assert_eq!(a.proteoform_fdr(), 0.06);
        assert_eq!(a.protein_fdr, Some(0.03));
        assert_eq!(a.mbr_donor_fdr(), 0.005);
        assert_eq!(a.quant_transfer_fdr(), 0.04);
        assert_eq!(a.quant_channel_fdr(), 0.04);
        let a = applied(
            "quant:\n  quant_transfer_fdr: 0.02\n  quant_channel_fdr: 0.03\n",
            &["--quant-transfer-fdr", "0.07"],
        )
        .unwrap();
        assert_eq!(a.quant_transfer_fdr(), 0.07);
        assert_eq!(a.quant_channel_fdr(), 0.03);
        assert_eq!(a.mbr_donor_fdr(), 0.01);
        for old in [
            "rescoring:\n  fdr: 0.01\n",
            "refine:\n  select_psm_fdr: 0.01\n",
            "quant:\n  mbr_fdr: 0.01\n",
            "quant:\n  mbr_library_fdr: 0.01\n",
            "quant:\n  label_fdr: 0.05\n",
            "quant:\n  label_mbr_fdr: 0.05\n",
            "quant:\n  lfq_feature_fdr: 0.05\n",
            "quant:\n  quant_feature_fdr: off\n",
        ] {
            assert!(
                RunConfig::from_yaml_str(old).is_err(),
                "{old} must be rejected"
            );
        }
    }

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
        assert!(!a.db.uniprot_variants);
        assert_eq!(a.db.uniprot_evidence, input::EvidenceLevel::Experimental);
        assert!(!a.db.uniprot_isoforms && !a.db.uniprot_chains);
        let cfg = RunConfig::from_yaml_str(
            "io:\n  uniprot_variants: true\n  uniprot_isoforms: true\n  uniprot_chains: true\n  \
             uniprot_evidence: curated\n",
        )
        .unwrap();
        apply(cfg, &mut a, &matches).unwrap();
        assert!(a.db.uniprot_variants && a.db.uniprot_isoforms && a.db.uniprot_chains);
        assert_eq!(a.db.uniprot_evidence, input::EvidenceLevel::Curated);

        let bad = RunConfig::from_yaml_str("io:\n  uniprot_evidence: strong\n").unwrap();
        let mut a = crate::SearchArgs::from_arg_matches(&matches).unwrap();
        assert!(apply(bad, &mut a, &matches).is_err());

        let cli = crate::TopCli::command()
            .try_get_matches_from(base.into_iter().chain(["--uniprot-evidence", "any"]))
            .unwrap();
        let mut a = crate::SearchArgs::from_arg_matches(&cli).unwrap();
        let cfg = RunConfig::from_yaml_str("io:\n  uniprot_evidence: curated\n").unwrap();
        apply(cfg, &mut a, &cli).unwrap();
        assert_eq!(a.db.uniprot_evidence, input::EvidenceLevel::Any);
    }

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

    const CLI_ONLY: &[&str] = &[
        "config",
        "cal_min_spec_keys",
        "precursor_mono_dump",
        "glyco_scans",
        "glyco_mbr_match_train",
        "ipf_debug",
        "xl_oracle",
        "xl_rescore_intra",
        "xl_rescore_inter",
    ];

    #[test]
    fn every_search_flag_has_a_config_key() {
        use clap::CommandFactory;
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
    fn modification_lists_apply_with_cli_precedence() {
        let yaml = "search:\n  fixed_mods: [\"Carbamidomethyl (C)\", {name: MyTag, mass: 28.0313, sites: N-term}]\n  \
                    variable_mods:\n    - Oxidation (M)\n    \
                    - {name: \"Phospho (ST)\", neutral_loss: [97.976896], loss_class: phospho, restrict_to: annotated}\n    \
                    - {name: GlyGly, sites: K, cleavage: blocked, diagnostic_ions: 215.1502}\n";
        let a = applied(yaml, &[]).unwrap();
        let f = &a.db.fixed_mods.as_ref().unwrap().0;
        assert_eq!(f.len(), 2);
        assert_eq!((f[0].name.as_str(), f[0].mass_delta), ("Carbamidomethyl", 57.021464));
        assert_eq!((f[1].name.as_str(), f[1].mass_delta), ("MyTag", 28.0313));
        assert!(f[1].accession.is_none() && f.iter().all(|m| m.fixed));
        let v = &a.db.variable_mods.as_ref().unwrap().0;
        assert_eq!(v.len(), 4);
        assert_eq!(v[1].neutral_losses, vec![97.976896]);
        assert_eq!(v[1].loss_class, 2);
        assert_eq!(v[1].annotated_site_mass, Some(79.966331));
        assert_eq!(v[3].name, "GG");
        assert!(v[3].blocks_cleavage);
        assert_eq!(v[3].diagnostic_ions, vec![215.1502]);
        let a = applied(yaml, &["--variable-mods", "Phospho (Y)"]).unwrap();
        assert_eq!(a.db.variable_mods.unwrap().0.len(), 1);
        assert_eq!(a.db.fixed_mods.unwrap().0.len(), 2);
        let a = applied("search:\n  variable_mods: []\n", &[]).unwrap();
        assert!(a.db.variable_mods.unwrap().0.is_empty() && a.db.fixed_mods.is_none());
        for (bad, msg) in [
            (
                "search:\n  variable_mods: [{name: Phospho, sites: S, nuetral_loss: [1.0]}]\n",
                "unknown key `nuetral_loss`",
            ),
            ("search:\n  variable_mods: [Phosho (S)]\n", "did you mean Phospho"),
            ("search:\n  variable_mods: [{mass: 1.0}]\n", "`sites` is missing"),
            (
                "search:\n  fixed_mods: [{name: Phospho, sites: S, restrict_to: annotated}]\n",
                "needs a variable",
            ),
            ("search:\n  variable_mods: [42]\n", "expected `Name (sites)`"),
        ] {
            let e = applied(bad, &[]).unwrap_err();
            assert!(e.contains(msg), "{e}");
        }
    }

    #[test]
    fn crosslinker_applies_with_cli_precedence() {
        let a = applied("search:\n  xl: DSSO (K, S, T, Y, Protein N-term)\n", &[]).unwrap();
        let x = a.xl.unwrap();
        assert_eq!(x.name, "DSSO");
        assert_eq!(x.ends[0].residues, b"KSTY");
        let a = applied("search:\n  xl: DSSO (K)\n", &["--xl", "DSBU"]).unwrap();
        assert_eq!(a.xl.unwrap().name, "DSBU");
        let yaml = "search:\n  xl:\n    name: MyLinker\n    sites: K, Protein N-term\n    \
                    intact: {mass: 158.003765}\n    \
                    stubs: [{name: alkene, mass: 54.010565}, {name: thiol, mass: 85.982635}]\n    \
                    doublets: [[alkene, thiol]]\n    monolinks: [{name: water, mass: 176.01433}]\n";
        let x = applied(yaml, &[]).unwrap().xl.unwrap();
        assert_eq!(x.name, "MyLinker");
        let c = x.cleavage.unwrap();
        assert_eq!(c.stubs.len(), 2);
        assert_eq!(c.doublets, vec![(0, 1)]);
        assert_eq!(x.monolinks.len(), 1);
        let x = applied("search:\n  xl: {name: \"BS3 (K)\"}\n", &[])
            .unwrap()
            .xl
            .unwrap();
        assert_eq!(x.name, "DSS");
        assert!(!x.ends[0].protein_n_term);
        let x = applied("search:\n  xl: {name: EDC, sites: K / D}\n", &[])
            .unwrap()
            .xl
            .unwrap();
        assert_eq!(x.ends[1].residues, b"D");
        for (bad, msg) in [
            ("search:\n  xl: DSOS\n", "did you mean DSSO"),
            ("search:\n  xl: {name: X, intact: {mass: 100.0}}\n", "needs `sites`"),
            ("search:\n  xl: {name: X, sites: K, intact: {}}\n", "needs a `unimod`"),
            (
                "search:\n  xl: {name: DSSO, sitess: K}\n",
                "unknown key `sitess`",
            ),
            (
                "search:\n  xl: {name: X, sites: K, intact: {mass: 1.0}, stubz: []}\n",
                "unknown field `stubz`",
            ),
            ("search:\n  xl: 42\n", "expected `Name (sites)`"),
        ] {
            let e = applied(bad, &[]).unwrap_err();
            assert!(e.contains(msg), "{e}");
        }
    }

    #[test]
    fn new_keys_apply_and_validate() {
        let a = applied(
            "search:\n  fragment_index: \"off\"\n  fragment_index_slice_da: 40\n  index_cache_dir: shared/cache\n  \
             localize: ipf\n\
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
        assert_eq!(a.localize.as_deref(), Some("ipf"));
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

        let d = applied("{}", &[]).unwrap();
        assert_eq!(d.gbdt_max_trees, None);

        for bad in [
            "search:\n  fragment_index: sometimes\n",
            "search:\n  fragment_index_slice_da: 0\n",
            "search:\n  localize: ascore\n",
            "scoring:\n  precursor_mono: maybe\n",
            "quant:\n  lfq_gaussian_width: -1\n",
            "glyco:\n  retrieval_tol_ppm: 0\n",
            "glyco:\n  species: martian\n",
        ] {
            assert!(applied(bad, &[]).is_err(), "{bad} must be rejected");
        }
        let o = applied("search:\n  localize: ascore\n", &["--localize", "ipf"]).unwrap();
        assert_eq!(o.localize.as_deref(), Some("ipf"));

        let a = applied(
            "scoring:\n  gbdt_max_trees: 0\n",
            &["--gbdt-max-trees", "50"],
        )
        .unwrap();
        assert_eq!(a.gbdt_max_trees, Some(50));
    }

    #[test]
    fn config_memory_gb_inherits() {
        let a = applied("rescoring:\n  keep_pin: true\n", &[]).unwrap();
        assert_eq!(a.rescore_memory_gb, 16);
        let a = applied("rescoring:\n  memory_gb: 4\n", &[]).unwrap();
        assert_eq!(a.rescore_memory_gb, 4);
        let a = applied(
            "rescoring:\n  memory_gb: 4\n",
            &["--rescore-memory-gb", "8"],
        )
        .unwrap();
        assert_eq!(a.rescore_memory_gb, 8);
        assert!(applied("rescoring:\n  memory_gb: 0\n", &[]).is_err());
    }

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

        let a = applied("glyco:\n  glycan_wurcs: w.csv\n", &["--glyco"]).unwrap();
        assert_eq!(
            a.glyco_glycan_wurcs.as_deref(),
            Some(std::path::Path::new("w.csv"))
        );
        assert!(check_resolved(&a).is_ok());
        let a = applied("{}", &["--glyco-glycan-wurcs", "w.csv"]).unwrap();
        assert_eq!(
            a.glyco_glycan_wurcs.as_deref(),
            Some(std::path::Path::new("w.csv"))
        );

        let both = "glyco:\n  enabled: true\n  glycan_gdb: g.gdb\n  species: human\n";
        assert!(applied(both, &[])
            .unwrap_err()
            .contains("mutually exclusive"));
        let both = "glyco:\n  enabled: true\n  glycan_gdb: g.gdb\n  glycan_wurcs: w.csv\n";
        assert!(applied(both, &[])
            .unwrap_err()
            .contains("mutually exclusive"));
        let a = applied(
            "glyco:\n  glycan_gdb: g.gdb\n",
            &["--glyco-species", "mouse"],
        )
        .unwrap();
        assert_eq!(a.glyco_glycan_gdb, None);
        assert_eq!(a.glyco_species, Some(crate::cli::GlycoSpeciesFlag::Mouse));
    }
}
