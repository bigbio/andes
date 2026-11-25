pub(crate) mod mods;
mod sections;
pub use sections::apply;

use std::path::PathBuf;

use serde::Deserialize;

use crate::params::sources::ParamSources;
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
    /// Read by `andes cohort plan`; a search ignores it.
    pub cohort: Option<crate::cohort::CohortCfg>,
    /// The file as parsed, so a recorded setting keeps the value as written.
    #[serde(skip)]
    pub written: serde_yaml::Value,
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
    pub xl_linear_only: Option<bool>,
    pub localize: Option<String>,
    pub fragment_index: Option<String>,
    pub fragment_index_slice_da: Option<f64>,
    pub index_cache_dir: Option<PathBuf>,
    pub dia: Option<bool>,
    pub dia_mobility_slices: Option<usize>,
    pub dia_ms1_tol_ppm: Option<f64>,
    pub open_search: Option<bool>,
    pub open_search_delta: Option<String>,
    pub open_search_top_k: Option<u32>,
    pub optimize: Option<bool>,
    pub optimize_sample: Option<usize>,
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
        let mut cfg: Self =
            serde_yaml::from_str(s).map_err(|e| format!("--config parse error: {e}"))?;
        cfg.written = serde_yaml::from_str(s).unwrap_or_default();
        Ok(cfg)
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
    src: &mut ParamSources,
) -> Result<(), String> {
    macro_rules! set {
        ($id:literal, $field:ident, $val:expr) => {
            if let Some(v) = &$val {
                if !cli_set(m, $id) {
                    src.record_config($id, &cfg.written);
                    db.$field = v.clone();
                }
            }
        };
    }
    macro_rules! set_opt {
        ($id:literal, $field:ident, $val:expr) => {
            if let Some(v) = &$val {
                if !cli_set(m, $id) {
                    src.record_config($id, &cfg.written);
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
                    src.record_config($id, &cfg.written);
                }
            }
        };
    }
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
    set_parsed!(
        "enzyme",
        enzyme,
        cfg.search.enzyme,
        crate::cli::parse_enzyme_list
    );
    if let Some(items) = &cfg.search.fixed_mods {
        if !cli_set(m, "fixed_mods") {
            db.fixed_mods = Some(mod_list("search.fixed_mods", items, true)?);
            src.record_config("fixed_mods", &cfg.written);
        }
    }
    if let Some(items) = &cfg.search.variable_mods {
        if !cli_set(m, "variable_mods") {
            db.variable_mods = Some(mod_list("search.variable_mods", items, false)?);
            src.record_config("variable_mods", &cfg.written);
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
        let cfg = RunConfig::load(&path)?;
        apply_db(&cfg, &mut p.db, m, &mut ParamSources::default())?;
    }
    Ok(())
}

use crate::config::mods::{mod_list, xl_set};

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

    fn apply(cfg: RunConfig, a: &mut SearchArgs, m: &clap::ArgMatches) -> Result<(), String> {
        super::apply(cfg, a, m, &mut ParamSources::default())
    }

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

    fn cli_only(id: &str) -> bool {
        crate::params::spec(id).is_some_and(|p| p.section == crate::params::Section::CliOnly)
    }

    #[test]
    fn every_search_flag_has_a_config_key() {
        use clap::CommandFactory;
        let source = [
            include_str!("config.rs"),
            include_str!("config/sections.rs"),
        ]
        .concat()
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
        .replace("( ", "(");
        let mapped = |id: &str| {
            ["set!", "set_opt!", "set_parsed!"].iter().any(|mac| {
                source.contains(&format!("{mac}(\"{id}\","))
                    || source.contains(&format!("{mac}(c, \"{id}\","))
            }) || source.contains(&format!("cli_set(m, \"{id}\")"))
                || source.contains(&format!("cli_set(c.m, \"{id}\")"))
        };
        let missing: Vec<String> = crate::TopCli::command()
            .get_arguments()
            .map(|a| a.get_id().to_string())
            .filter(|id| !matches!(id.as_str(), "help" | "version"))
            .filter(|id| !cli_only(id) && !mapped(id))
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
        assert_eq!(
            (f[0].name.as_str(), f[0].mass_delta),
            ("Carbamidomethyl", 57.021464)
        );
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
            (
                "search:\n  variable_mods: [Phosho (S)]\n",
                "did you mean Phospho",
            ),
            (
                "search:\n  variable_mods: [{mass: 1.0}]\n",
                "`sites` is missing",
            ),
            (
                "search:\n  fixed_mods: [{name: Phospho, sites: S, restrict_to: annotated}]\n",
                "needs a variable",
            ),
            (
                "search:\n  variable_mods: [42]\n",
                "expected `Name (sites)`",
            ),
        ] {
            let e = applied(bad, &[]).unwrap_err();
            assert!(e.contains(msg), "{e}");
        }
    }

    #[test]
    fn crosslinker_applies_with_cli_precedence() {
        let a = applied("search:\n  xl: DSSO (K, S, T, Y, Protein N-term)\n", &[]).unwrap();
        let x = a.xl.unwrap().0.remove(0);
        assert_eq!(x.name, "DSSO");
        assert_eq!(x.ends[0].residues, b"KSTY");
        let a = applied("search:\n  xl: DSSO (K)\n", &["--xl", "DSBU"]).unwrap();
        assert_eq!(a.xl.unwrap().0.remove(0).name, "DSBU");
        let yaml = "search:\n  xl:\n    name: MyLinker\n    sites: K, Protein N-term\n    \
                    intact: {mass: 158.003765}\n    \
                    stubs: [{name: alkene, mass: 54.010565}, {name: thiol, mass: 85.982635}]\n    \
                    doublets: [[alkene, thiol]]\n    monolinks: [{name: water, mass: 176.01433}]\n";
        let x = applied(yaml, &[]).unwrap().xl.unwrap().0.remove(0);
        assert_eq!(x.name, "MyLinker");
        let c = x.cleavage.unwrap();
        assert_eq!(c.stubs.len(), 2);
        assert_eq!(c.doublets, vec![(0, 1)]);
        assert_eq!(x.monolinks.len(), 1);
        let x = applied("search:\n  xl: {name: \"BS3 (K)\"}\n", &[])
            .unwrap()
            .xl
            .unwrap()
            .0
            .remove(0);
        assert_eq!(x.name, "DSS");
        assert!(!x.ends[0].protein_n_term);
        let x = applied("search:\n  xl: {name: EDC, sites: K / D}\n", &[])
            .unwrap()
            .xl
            .unwrap()
            .0
            .remove(0);
        assert_eq!(x.ends[1].residues, b"D");
        for (bad, msg) in [
            ("search:\n  xl: DSOS\n", "did you mean DSSO"),
            ("search:\n  xl: []\n", "search.xl: no crosslinker"),
            (
                "search:\n  xl: [DSSO, \"DSSO (K)\"]\n",
                "search.xl: DSSO is given twice",
            ),
            (
                "search:\n  xl: {name: X, intact: {mass: 100.0}}\n",
                "needs `sites`",
            ),
            (
                "search:\n  xl: {name: X, sites: K, intact: {}}\n",
                "needs a `unimod`",
            ),
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
    fn open_search_keys_apply() {
        let a = applied(
            "search:\n  open_search: true\n  open_search_delta: \"-50..200\"\n  open_search_top_k: 50\n",
            &[],
        )
        .unwrap();
        assert!(a.open_search);
        assert_eq!(a.open_search_delta, (-50.0, 200.0));
        assert_eq!(a.open_search_top_k, 50);
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
        assert_eq!(a.rescore_memory_gb, None);
        let a = applied("rescoring:\n  memory_gb: 4\n", &[]).unwrap();
        assert_eq!(a.rescore_memory_gb, Some(4));
        let a = applied(
            "rescoring:\n  memory_gb: 4\n",
            &["--rescore-memory-gb", "8"],
        )
        .unwrap();
        assert_eq!(a.rescore_memory_gb, Some(8));
        assert!(applied("rescoring:\n  memory_gb: 0\n", &[]).is_err());
    }

    #[test]
    fn rescoring_backend_comes_from_the_binary_keys() {
        let a = applied("rescoring:\n  enabled: true\n", &[]).unwrap();
        assert!(a.rescore && !a.external_percolator());
        let a = applied(
            "rescoring:\n  enabled: true\n  percolator_bin: /opt/percolator\n",
            &[],
        )
        .unwrap();
        assert!(a.rescore && a.external_percolator());
        let a = applied(
            "rescoring:\n  enabled: true\n  percolator_docker: true\n",
            &[],
        )
        .unwrap();
        assert!(a.rescore && a.external_percolator());
        let a = applied("rescoring:\n  enabled: false\n", &["--rescore"]).unwrap();
        assert!(a.rescore, "CLI --rescore must win");
        assert!(
            applied("rescoring:\n  native: true\n", &[]).is_err(),
            "rescoring.native is no longer a key"
        );
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
