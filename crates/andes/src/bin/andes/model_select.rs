//! Scoring-model selection: amino-acid set construction and Parquet model-store lookup.

use std::path::{Path, PathBuf};

use crate::cli::{Fragmentation, Protocol};
use model::{
    activation::ActivationMethod, AminoAcidSetBuilder, InstrumentType, ModLocation, Modification,
    ResidueSpec,
};
use model_train::{
    counts::CountStats,
    protocol_to_experiment_class as store_protocol_to_experiment_class,
    select::{select, select_nearest, SelectionKey},
    store::SourceLedger,
    ModelStore,
};
use scoring_crate::Param;

// Type alias to reduce clippy type_complexity warnings in the train path.
pub(crate) type ModelEntryOwned = (String, Param, Vec<(SourceLedger, CountStats)>);

/// Load the seed Param from the optional seed model specifier.
pub(crate) fn load_seed_param(
    seed_model: &Option<String>,
) -> Result<(String, Param), Box<dyn std::error::Error>> {
    match seed_model {
        None => {
            let store_path = bundled_store_path();
            let store =
                ModelStore::open(&store_path).map_err(|e| format!("opening bundled store: {e}"))?;
            let p = store
                .load_param("hcd_qexactive_tryp")
                .map_err(|e| format!("loading seed model: {e}"))?;
            Ok(("hcd_qexactive_tryp".to_string(), p))
        }
        Some(seed) => {
            // Seed by slug from the canonical Parquet store.
            let store_path = bundled_store_path();
            let store =
                ModelStore::open(&store_path).map_err(|e| format!("opening bundled store: {e}"))?;
            let p = store
                .load_param(seed)
                .map_err(|e| format!("loading seed model '{seed}': {e}"))?;
            Ok((seed.clone(), p))
        }
    }
}

/// Build an `AminoAcidSet` from an optional mods file for the training paths,
/// defaulting to Carbamidomethyl-C fixed + Oxidation-M variable. (The search
/// default, `default_aa_set_with_tag`, additionally carries protein-N-term Acetyl.)
pub(crate) fn build_aa_set(
    mods: &Option<PathBuf>,
) -> Result<model::AminoAcidSet, Box<dyn std::error::Error>> {
    match mods {
        Some(path) => {
            let set = AminoAcidSetBuilder::new_standard()
                .add_mods_from_file(path)
                .map_err(|e| format!("loading mods from {}: {e}", path.display()))?
                .build()
                .map_err(|e| format!("building amino-acid set: {e}"))?;
            Ok(set)
        }
        None => {
            let cam = Modification {
                name: "Carbamidomethyl".into(),
                mass_delta: 57.02146,
                residue: ResidueSpec::Specific(b'C'),
                location: ModLocation::Anywhere,
                fixed: true,
                accession: None,
                neutral_losses: Vec::new(),
                loss_class: 0,
            };
            let ox = Modification {
                name: "Oxidation".into(),
                mass_delta: 15.99491,
                residue: ResidueSpec::Specific(b'M'),
                location: ModLocation::Anywhere,
                fixed: false,
                accession: None,
                neutral_losses: Vec::new(),
                loss_class: 0,
            };
            Ok(AminoAcidSetBuilder::new_standard()
                .add_fixed_mod(cam)
                .add_variable_mod(ox)
                .build()?)
        }
    }
}

/// Build the default `AminoAcidSet` (CAM-C fixed + Ox-M variable), optionally
/// with an isobaric tag (TMT/iTRAQ) as a FIXED mod on K + peptide N-term.
///
/// Used by the no-`--mods` (parameter-free) path: when the protocol resolves to
/// TMT/iTRAQ the tag MUST be in the candidate set, or every labeled peptide is
/// `+tag·(nK+1)` Da off at the precursor and silently misses (C1).
pub(crate) fn default_aa_set_with_tag(
    tag: Option<(&str, f64)>,
) -> Result<model::AminoAcidSet, Box<dyn std::error::Error>> {
    let cam = Modification {
        name: "Carbamidomethyl".into(),
        mass_delta: 57.02146,
        residue: ResidueSpec::Specific(b'C'),
        location: ModLocation::Anywhere,
        fixed: true,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
    };
    let ox = Modification {
        name: "Oxidation".into(),
        mass_delta: 15.99491,
        residue: ResidueSpec::Specific(b'M'),
        location: ModLocation::Anywhere,
        fixed: false,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
    };
    // Protein N-terminal acetylation (+42.010565) — a near-universal default in
    // the field (the reference engine / a comparison engine / a quantitation tool). Restricted to the PROTEIN N-term
    // (one site per protein, after Met cleavage), so it is combinatorially cheap
    // (not every peptide N-term) and biologically correct.
    let acetyl = Modification {
        name: "Acetyl".into(),
        mass_delta: 42.010565,
        residue: ResidueSpec::Wildcard,
        location: ModLocation::ProtNTerm,
        fixed: false,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
    };
    let mut b = AminoAcidSetBuilder::new_standard()
        .add_fixed_mod(cam)
        .add_variable_mod(ox)
        .add_variable_mod(acetyl);
    if let Some((name, mass)) = tag {
        let tag_k = Modification {
            name: name.into(),
            mass_delta: mass,
            residue: ResidueSpec::Specific(b'K'),
            location: ModLocation::Anywhere,
            fixed: true,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
        };
        let tag_nterm = Modification {
            name: name.into(),
            mass_delta: mass,
            residue: ResidueSpec::Wildcard,
            location: ModLocation::NTerm,
            fixed: true,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
        };
        b = b.add_fixed_mod(tag_k).add_fixed_mod(tag_nterm);
    }
    Ok(b.build()?)
}

/// Convert the CLI `Fragmentation` enum to `Option<ActivationMethod>`.
///
/// `Fragmentation::Auto` returns `None` (no activation explicitly requested);
/// every concrete variant maps to its `ActivationMethod`. Used by
/// [`resolve_metadataless_selection`] so that an unset `--fragmentation`
/// defers to detection or the class-consistent default.
pub(crate) fn cli_fragmentation_to_activation_opt(f: Fragmentation) -> Option<ActivationMethod> {
    match f {
        Fragmentation::Auto => None,
        Fragmentation::Cid => Some(ActivationMethod::CID),
        Fragmentation::Etd => Some(ActivationMethod::ETD),
        Fragmentation::Hcd => Some(ActivationMethod::HCD),
        Fragmentation::Uvpd => Some(ActivationMethod::UVPD),
    }
}

/// Resolve the CLI fragment-tolerance override (MGF only) into a `Tolerance`.
/// `--fragment-tol-ppm` ⇒ `Ppm`; `--fragment-tol-da` ⇒ `Da`; none ⇒ `None`.
pub(crate) fn cli_fragment_tol_override(
    fragment_tol_ppm: Option<f64>,
    fragment_tol_da: Option<f64>,
) -> Option<model::tolerance::Tolerance> {
    use model::tolerance::Tolerance;
    fragment_tol_ppm
        .map(Tolerance::Ppm)
        .or_else(|| fragment_tol_da.map(Tolerance::Da))
}

/// Parse the `--enzyme` value into `(primary, extras)`.
///
/// Accepts `,` or `+` separators (the latter matching a quantitation tool/the reference engine
/// docs), trims whitespace, drops empty tokens, and de-duplicates by enzyme
/// while preserving order. The first surviving entry is the primary (drives
/// model selection and the cleavage-credit PIN feature); the rest only widen
/// candidate enumeration. Errors on an unknown name or empty input.
pub(crate) fn parse_enzymes(
    spec: &str,
) -> Result<(model::enzyme::Enzyme, Vec<model::enzyme::Enzyme>), Box<dyn std::error::Error>> {
    use model::enzyme::Enzyme;
    let mut all: Vec<Enzyme> = Vec::new();
    for tok in spec
        .split([',', '+'])
        .map(str::trim)
        .filter(|s| !s.is_empty())
    {
        let e = Enzyme::from_name(tok).ok_or_else(|| {
            format!(
                "unknown --enzyme '{tok}' (expected trypsin/chymotrypsin/lysc/aspn/gluc/lysn/argc/\
             alphalp/nocleavage/nonspecific/elastase)"
            )
        })?;
        if !all.contains(&e) {
            all.push(e); // dedup, keep order (first = primary)
        }
    }
    let (primary, extras) = all.split_first().ok_or("empty --enzyme")?;
    Ok((*primary, extras.to_vec()))
}

/// Warn when a universal protease (NonSpecific, or AlphaLP which cleaves every
/// bond) is combined with specific protease(s): the cleavage-site union then
/// covers EVERY position, so the digest is effectively non-specific — usually
/// not what a user adding it to a specific list intends.
pub(crate) fn warn_if_universal_protease_combo(
    primary: model::enzyme::Enzyme,
    extras: &[model::enzyme::Enzyme],
) {
    use model::enzyme::Enzyme;
    let is_universal = |e: &Enzyme| matches!(e, Enzyme::NonSpecific | Enzyme::AlphaLP);
    let is_specific = |e: &Enzyme| {
        !matches!(
            e,
            Enzyme::NonSpecific | Enzyme::AlphaLP | Enzyme::NoCleavage
        )
    };
    let all: Vec<Enzyme> = std::iter::once(primary)
        .chain(extras.iter().copied())
        .collect();
    if all.len() > 1 && all.iter().any(is_universal) && all.iter().any(is_specific) {
        eprintln!(
            "WARN: --enzyme combines a universal protease (NonSpecific/AlphaLP) with specific \
             one(s); their union cleaves at EVERY position, so this is effectively a \
             non-specific digest (much larger search space). Drop the universal protease to \
             keep the specific cleavage rules."
        );
    }
}

/// Resolve (activation, instrument) for model selection on metadata-less input
/// (MGF, or mzML/.raw with no analyzer metadata). Resolution class comes from
/// the `--fragment-tol-*` unit; activation from detected method, else
/// `--fragmentation`, else the class-consistent default. When nothing
/// disambiguates, defaults to CID / LowRes (→ `cid_lowres_tryp`) + a warning.
pub(crate) fn resolve_metadataless_selection(
    detected_activation: Option<ActivationMethod>,
    fragmentation: Fragmentation,
    fragment_tol_ppm: Option<f64>,
    fragment_tol_da: Option<f64>,
) -> (ActivationMethod, Option<InstrumentType>) {
    let instrument: Option<InstrumentType> = if fragment_tol_ppm.is_some() {
        Some(InstrumentType::QExactive)
    } else if fragment_tol_da.is_some() {
        Some(InstrumentType::LowRes)
    } else {
        None
    };
    let explicit = cli_fragmentation_to_activation_opt(fragmentation);
    // Class-consistent default when neither detection nor `--fragmentation`
    // names an activation: high-res classes imply HCD, otherwise CID.
    let class_default = match instrument {
        Some(InstrumentType::QExactive)
        | Some(InstrumentType::HighRes)
        | Some(InstrumentType::TOF) => ActivationMethod::HCD,
        _ => ActivationMethod::CID,
    };
    let activation = detected_activation.or(explicit).unwrap_or(class_default);
    if detected_activation.is_none() && explicit.is_none() && instrument.is_none() {
        eprintln!(
            "WARN: MGF input with no --fragmentation/--fragment-tol; assuming \
             CID / low-res / 0.5 Da. Pass --fragmentation and --fragment-tol-ppm/-da \
             to override."
        );
    }
    (activation, instrument)
}

/// Resolve the path to the bundled model store.
///
/// The store ships as a per-protocol partitioned directory `resources/models/`
/// (Hive-style `protocol=<P>/models.parquet`), which [`ModelStore::open`]
/// reads as one store.
///
/// A packaged release ships `resources/` next to the binary, so prefer
/// `<exe_dir>/resources/...` when it exists — that makes an installed binary
/// self-contained regardless of where it runs. Fall back to the compile-time
/// source tree (`CARGO_MANIFEST_DIR`) for `cargo run` / tests.
pub(crate) fn bundled_store_path() -> PathBuf {
    // Search roots in priority order: next to the binary, then the source tree.
    let mut roots: Vec<PathBuf> = Vec::new();
    if let Ok(exe) = std::env::current_exe() {
        if let Some(dir) = exe.parent() {
            roots.push(dir.join("resources"));
        }
    }
    roots.push(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../resources"));

    for root in &roots {
        let partitioned = root.join("models");
        if partitioned.is_dir() {
            return partitioned;
        }
    }

    // Last-resort default (source tree directory) for error messages.
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../resources/models")
}

/// Build a [`SelectionKey`] from `(activation, instrument, protocol)` applying
/// all old-ladder normalizations. This is the new entry point used by the
/// search binary, replacing the former filename-based resolution ladder.
///
/// `activation`: the detected or explicitly set `ActivationMethod`.
/// `instrument`: the detected or explicitly set `InstrumentType` (None = undetected → LowRes).
/// `protocol`:   the CLI `Protocol` value.
pub(crate) fn build_selection_key(
    activation: ActivationMethod,
    instrument: Option<InstrumentType>,
    protocol: Protocol,
    enzyme: model::enzyme::Enzyme,
) -> SelectionKey {
    use std::collections::BTreeSet;

    // 1. PQD → CID for model routing.
    let act_str: &str = match activation {
        ActivationMethod::PQD => "CID",
        other => other.name(),
    };
    // 2. Apply family fallback (OrbitrapAstral → QExactive, TimsTOF → TOF).
    let inst_after_family: &str = match instrument {
        Some(i) => i.family_fallback().name(),
        None => "LowRes",
    };

    // 3. Apply old-ladder (activation, instrument) normalization.
    //    Because `normalize_for_store` returns `&'static str` only for the
    //    normalizing arms (to avoid lifetime issues), we handle identity
    //    inline here.
    let (final_act, final_inst, drop_protocol): (&str, &str, bool) =
        match (act_str, inst_after_family) {
            // H5: low-res (ion-trap) HCD. Routing this to the high-res
            // QExactive model would match 0.5-Da peaks at 20 ppm and lose
            // ~18% of PSMs silently. No hcd_lowres model exists, so route to the
            // low-res b/y model (cid_lowres_tryp) instead — correct fragment
            // tolerance and ion series. Pinned by store_selection_equivalence.rs.
            ("HCD", "LowRes") => {
                eprintln!(
                    "WARN: low-res (ion-trap) HCD detected — no hcd_lowres model exists; \
                     routing to cid_lowres_tryp (low-res b/y, 0.5-Da tolerance) rather than \
                     the high-res QExactive model. Pass --model to override."
                );
                ("CID", "LowRes", true)
            }
            ("HCD", "TOF") => ("CID", "TOF", true),
            ("CID", "QExactive") => ("CID", "LowRes", true),
            ("ETD", i) if !matches!(i, "LowRes" | "HighRes") => ("ETD", "LowRes", true),
            ("UVPD", i) if i != "QExactive" => ("CID", "LowRes", true),
            _ => (act_str, inst_after_family, false),
        };

    // 4. Build experiment_class from protocol (unless the final fallback dropped it).
    //    Protocol → experiment_class mapping matches the parquet's `protocol` column.
    let protocol_for_store: &str = match protocol {
        Protocol::Auto | Protocol::Standard => "Automatic",
        Protocol::Tmt => "TMT",
        Protocol::Phospho => "Phosphorylation",
        Protocol::Itraq => "iTRAQ",
        Protocol::ItraqPhospho => "iTRAQPhospho",
    };
    let experiment_class: BTreeSet<String> = if drop_protocol {
        BTreeSet::new()
    } else {
        store_protocol_to_experiment_class(protocol_for_store)
    };

    SelectionKey {
        activation: final_act.to_string(),
        instrument: final_inst.to_string(),
        // Parquet stores the enzyme as its `Enzyme::name()` ("Trypsin", "LysC", ...).
        enzyme: enzyme.name().to_string(),
        experiment_class,
    }
}

/// Load the scoring [`Param`] from the bundled Parquet store for the given
/// `(activation, instrument, protocol)` combination.
///
/// The `model_id` selected from the store is guaranteed to match the
/// reference resolution ladder by the equivalence gate test
/// `store_selection_matches_old_ladder_for_all_combos`.
///
/// `custom_store_path`: when `Some`, use that Parquet store instead of the
/// bundled `resources/models/` (honours `--model-store`).
///
/// `model_id_override`: when `Some`, skip automatic selection and load this
/// exact model ID (honours `--model`).
pub(crate) fn load_param_from_store(
    activation: ActivationMethod,
    instrument: Option<InstrumentType>,
    protocol: Protocol,
    enzyme: model::enzyme::Enzyme,
    custom_store_path: Option<&Path>,
    model_id_override: Option<&str>,
) -> Result<(String, Param), Box<dyn std::error::Error>> {
    let store_path = custom_store_path
        .map(|p| p.to_owned())
        .unwrap_or_else(bundled_store_path);
    let store = ModelStore::open(&store_path)
        .map_err(|e| format!("opening model store {}: {e}", store_path.display()))?;

    let model_id: String = if let Some(id) = model_id_override {
        id.to_string()
    } else {
        let entries = store.selection_entries();
        let key = build_selection_key(activation, instrument, protocol, enzyme);

        // Forward-compat: `build_selection_key` collapses instruments with a real
        // family fallback (OrbitrapAstral → QExactive, TimsTOF → TOF) so the
        // bundled models resolve correctly. But that also hides a model trained
        // for the EXACT instrument (e.g. a user-trained OrbitrapAstral model).
        // Try the exact detected instrument FIRST; only when no such model exists
        // (the bundled case) do we fall through to the normalized family ladder —
        // so bundled selection (and the equivalence gate) is unchanged.
        let exact_id: Option<String> = match instrument {
            Some(i) if i.family_fallback().name() != i.name() => {
                let raw_key = SelectionKey {
                    instrument: i.name().to_string(),
                    ..key.clone()
                };
                select(&entries, &raw_key, |s| s.to_string(), None).map(|s| s.to_string())
            }
            _ => None,
        };

        exact_id.unwrap_or_else(|| {
            // `build_selection_key` already applies family fallback + all
            // normalizations, so the family_fn here is the identity. When the
            // exact ladder misses (e.g. a protocol the own-only store doesn't
            // carry), select_nearest routes to the CLOSEST own model — relaxing
            // the enzyme (keeping the activation and instrument) — and only
            // resolves to the standard base as a last resort, WARNing which model
            // it substituted so the user can pin one with --model.
            let (id, substituted) =
                select_nearest(&entries, &key, |i| i.to_string(), "hcd_qexactive_tryp");
            if substituted {
                eprintln!(
                    "WARN: no model matched (activation={}, instrument={}, enzyme={}, class={:?}) \
                     — using the nearest available model '{}'; scores may be mis-calibrated for \
                     this data. Pin a model with --model if this is wrong.",
                    key.activation, key.instrument, key.enzyme, key.experiment_class, id
                );
            }
            id.to_string()
        })
    };

    let param = store
        .load_param(&model_id)
        .map_err(|e| format!("loading model '{model_id}' from store: {e}"))?;

    Ok((model_id, param))
}

#[cfg(test)]
mod param_resolver_tests {
    use super::*;
    use crate::cli::{
        parse_charge_range, parse_enzyme_specificity, parse_fragmentation,
        parse_isotope_error_range, parse_precursor_cal, parse_precursor_tol, parse_protocol,
        EnzymeSpecificity,
    };
    use crate::spectra::title_prefix_for;
    use crate::train_intensity::finalize_intensity_stats;
    use ::search::PrecursorCalMode;
    use model::Tolerance;

    // ── Model resolution is store-based: all bundled models live in
    //    resources/models/ and are selected by `model_id`.
    //    The store_selection_equivalence integration test covers the
    //    activation/instrument/protocol → model_id selection invariant.

    #[test]
    fn parse_fragmentation_accepts_named_rejects_numeric() {
        assert_eq!(parse_fragmentation("HCD").unwrap(), Fragmentation::Hcd);
        assert_eq!(parse_fragmentation("auto").unwrap(), Fragmentation::Auto);
        // Legacy numeric forms are no longer accepted.
        assert!(parse_fragmentation("3").is_err());
        assert!(parse_fragmentation("99").is_err());
    }

    #[test]
    fn parse_protocol_accepts_named_rejects_numeric() {
        assert_eq!(parse_protocol("TMT").unwrap(), Protocol::Tmt);
        assert_eq!(parse_protocol("auto").unwrap(), Protocol::Auto);
        assert!(parse_protocol("4").is_err());
        assert!(parse_protocol("99").is_err());
    }

    #[test]
    fn parse_enzyme_specificity_named_only() {
        assert_eq!(
            parse_enzyme_specificity("fully").unwrap(),
            EnzymeSpecificity::Fully
        );
        assert_eq!(
            parse_enzyme_specificity("semi").unwrap(),
            EnzymeSpecificity::Semi
        );
        assert!(parse_enzyme_specificity("2").is_err());
    }

    #[test]
    fn parse_charge_range_forms() {
        assert_eq!(parse_charge_range("2..5").unwrap(), (2, 5));
        assert_eq!(parse_charge_range("2-5").unwrap(), (2, 5));
        assert_eq!(parse_charge_range("3..3").unwrap(), (3, 3));
        assert!(parse_charge_range("5..2").is_err());
        assert!(parse_charge_range("x..5").is_err());
    }

    #[test]
    fn parse_isotope_error_range_allows_negatives() {
        assert_eq!(parse_isotope_error_range("-1..2").unwrap(), (-1, 2));
        assert_eq!(parse_isotope_error_range("-1-2").unwrap(), (-1, 2));
        assert_eq!(parse_isotope_error_range("0..0").unwrap(), (0, 0));
        assert!(parse_isotope_error_range("2..-1").is_err());
    }

    #[test]
    fn parse_precursor_tol_units() {
        assert_eq!(parse_precursor_tol("20ppm").unwrap(), Tolerance::Ppm(20.0));
        assert_eq!(parse_precursor_tol("0.02da").unwrap(), Tolerance::Da(0.02));
        assert_eq!(parse_precursor_tol("0.02Da").unwrap(), Tolerance::Da(0.02));
        assert!(parse_precursor_tol("20").is_err());
        assert!(parse_precursor_tol("xppm").is_err());
    }

    #[test]
    fn parse_precursor_cal_accepts_named_modes() {
        assert_eq!(parse_precursor_cal("auto").unwrap(), PrecursorCalMode::Auto);
        assert_eq!(parse_precursor_cal("OFF").unwrap(), PrecursorCalMode::Off);
        assert_eq!(parse_precursor_cal("on").unwrap(), PrecursorCalMode::On);
        assert!(parse_precursor_cal("bogus").is_err());
    }

    #[test]
    fn finalize_intensity_stats_drops_zero_count() {
        // A count==0 key (e.g. from a partial aggregation parquet) must not
        // produce NaN mean/var in the finalized model — it carries no signal.
        assert_eq!(finalize_intensity_stats(0.0, 0.0, 0), None);
    }

    #[test]
    fn finalize_intensity_stats_computes_mean_and_clamped_var() {
        // sum=6, sum_sq=14, count=2 -> mean=3, var=max(0, 14/2 - 9)= max(0,-2)=0.
        let (mean, var) = finalize_intensity_stats(6.0, 14.0, 2).unwrap();
        assert_eq!(mean, 3.0);
        assert_eq!(var, 0.0);
    }

    #[test]
    fn single_file_search_does_not_prefix_titles() {
        // Regression (default-path PIN parity): with one input file the
        // SpecId/title must stay unprefixed so PIN output is byte-identical to
        // the single-file behavior that predates multi-`--spectrum` support.
        assert_eq!(title_prefix_for(1, "myfile"), None);
    }

    #[test]
    fn multi_file_search_prefixes_titles_with_file_stem() {
        assert_eq!(title_prefix_for(2, "myfile").as_deref(), Some("myfile/"));
    }
}
