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
}

#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct IoCfg {
    pub spectrum: Option<Vec<PathBuf>>,
    pub database: Option<PathBuf>,
    pub output_pin: Option<PathBuf>,
    pub output_tsv: Option<PathBuf>,
    pub output_parquet: Option<PathBuf>,
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
    pub intensity_model: Option<PathBuf>,
    pub candidate_index: Option<String>,
    pub ms_level: Option<u8>,
    pub max_spectra: Option<usize>,
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
    pub pep: Option<f64>,
    pub percolator_bin: Option<PathBuf>,
    pub percolator_docker: Option<bool>,
    pub percolator_image: Option<String>,
    pub percolator_args: Option<String>,
    pub keep_pin: Option<bool>,
}

#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct GlycoCfg {
    pub enabled: Option<bool>,
    pub backbone_top_k: Option<usize>,
    pub gp_k: Option<f32>,
    pub gp_j: Option<f32>,
    pub gp_h: Option<f32>,
    pub gp_cz: Option<f32>,
    pub tol_ppm: Option<f64>,
    pub max_peaks: Option<usize>,
    pub pf_charge: Option<u8>,
    pub max_pf: Option<usize>,
    pub debug: Option<bool>,
    pub hcd_pair: Option<bool>,
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
    set_opt!("output_pin", output_pin, cfg.io.output_pin);
    set_opt!("output_tsv", output_tsv, cfg.io.output_tsv);
    set_opt!("output_parquet", output_parquet, cfg.io.output_parquet);

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
    set_opt!("model_store", model_store, cfg.scoring.model_store);
    set_opt!("model_id_override", model_id_override, cfg.scoring.model);
    set_opt!(
        "intensity_model",
        intensity_model,
        cfg.scoring.intensity_model
    );
    set_parsed!(
        "candidate_index",
        candidate_index,
        cfg.scoring.candidate_index,
        parse_value_enum::<crate::CandidateIndexFlag>
    );
    set!("ms_level", ms_level, cfg.scoring.ms_level);
    set!("max_spectra", max_spectra, cfg.scoring.max_spectra);

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
    set!("rescore", rescore, cfg.rescoring.enabled);
    set!("rescore_native", rescore_native, cfg.rescoring.native);
    set_opt!("fdr", fdr, cfg.rescoring.fdr);
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

    // ── glyco ──
    set!("glyco", glyco, cfg.glyco.enabled);
    set!(
        "glyco_backbone_top_k",
        glyco_backbone_top_k,
        cfg.glyco.backbone_top_k
    );
    set!("glyco_gp_k", glyco_gp_k, cfg.glyco.gp_k);
    set!("glyco_gp_j", glyco_gp_j, cfg.glyco.gp_j);
    set!("glyco_gp_h", glyco_gp_h, cfg.glyco.gp_h);
    set!("glyco_gp_cz", glyco_gp_cz, cfg.glyco.gp_cz);
    set!("glyco_tol_ppm", glyco_tol_ppm, cfg.glyco.tol_ppm);
    set!("glyco_max_peaks", glyco_max_peaks, cfg.glyco.max_peaks);
    set!("glyco_pf_charge", glyco_pf_charge, cfg.glyco.pf_charge);
    set!("glyco_max_pf", glyco_max_pf, cfg.glyco.max_pf);
    set!("debug_glyco", debug_glyco, cfg.glyco.debug);
    set!("glyco_hcd_pair", glyco_hcd_pair, cfg.glyco.hcd_pair);

    // Validate the resolved numeric fields that the CLI enforces via `value_parser`
    // but the YAML path (plain assignment) would otherwise bypass. Re-checking the
    // final value is a no-op for CLI-supplied values (already validated by clap) and
    // rejects out-of-range values injected through `--config`.
    if let Some(v) = args.fdr {
        check_unit_fraction("fdr", v)?;
    }
    if let Some(v) = args.pep {
        check_unit_fraction("pep", v)?;
    }
    check_unit_fraction("refine_select_psm_fdr", args.refine_select_psm_fdr)?;
    if let Some(v) = args.fragment_tol_ppm {
        check_positive("fragment_tol_ppm", v)?;
    }
    if let Some(v) = args.fragment_tol_da {
        check_positive("fragment_tol_da", v)?;
    }
    // `--fragment-tol-ppm` / `--fragment-tol-da` are mutually exclusive on the CLI
    // (clap `conflicts_with`); the config could set both, which is ambiguous.
    if args.fragment_tol_ppm.is_some() && args.fragment_tol_da.is_some() {
        return Err(
            "config: `scoring.fragment_tol_ppm` and `scoring.fragment_tol_da` are mutually \
             exclusive — set only one"
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
  gp_k: 42.0
";
        let c = RunConfig::from_yaml_str(y).expect("valid config");
        assert_eq!(c.search.precursor_tol.as_deref(), Some("15ppm"));
        assert_eq!(c.search.charge.as_deref(), Some("2..4"));
        assert_eq!(c.glyco.enabled, Some(true));
        assert_eq!(c.glyco.gp_k, Some(42.0));
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

    /// The shipped `config.example.yaml` must parse against the current schema, so
    /// the documented template can never drift from the accepted keys.
    #[test]
    fn shipped_example_config_parses() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../config.example.yaml");
        let s = std::fs::read_to_string(path).expect("read config.example.yaml");
        RunConfig::from_yaml_str(&s).expect("config.example.yaml must parse");
    }
}
