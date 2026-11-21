use std::path::{Path, PathBuf};

use crate::cli::{DatabaseArgs, Fragmentation, ModList, Protocol};
use model::{
    activation::ActivationMethod, AminoAcidSetBuilder, InstrumentType, ModLocation, Modification,
    ResidueSpec,
};
use model_train::{
    counts::CountStats,
    protocol_to_experiment_class as store_protocol_to_experiment_class,
    select::{select_nearest_for, SelectionEntry, SelectionKey},
    store::SourceLedger,
    ModelStore,
};
use scoring_crate::Param;

pub(crate) type ModelEntryOwned = (String, Param, Vec<(SourceLedger, CountStats)>);

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

pub(crate) const DEFAULT_FIXED_MODS: &str = "Carbamidomethyl (C)";

pub(crate) const DEFAULT_VARIABLE_MODS: &str = "Oxidation (M), Acetyl (Protein N-term)";

const TRAINING_VARIABLE_MODS: &str = "Oxidation (M)";

fn builtin_mods(list: &str, fixed: bool) -> Vec<Modification> {
    model::mod_spec::parse_mod_list(list, fixed)
        .unwrap_or_else(|e| panic!("built-in modification list {list:?}: {e}"))
}

pub(crate) fn default_fixed_mods() -> Vec<Modification> {
    builtin_mods(DEFAULT_FIXED_MODS, true)
}

pub(crate) fn default_variable_mods() -> Vec<Modification> {
    builtin_mods(DEFAULT_VARIABLE_MODS, false)
}

pub(crate) fn build_aa_set(
    fixed: Option<&ModList>,
    variable: Option<&ModList>,
) -> Result<model::AminoAcidSet, Box<dyn std::error::Error>> {
    let fixed = fixed.map_or_else(default_fixed_mods, |l| l.0.clone());
    let variable =
        variable.map_or_else(|| builtin_mods(TRAINING_VARIABLE_MODS, false), |l| l.0.clone());
    Ok(AminoAcidSetBuilder::new_standard()
        .add_mods(fixed)
        .add_mods(variable)
        .build()?)
}

pub(crate) fn search_aa_set(
    db: &DatabaseArgs,
    tag: Option<(&str, f64)>,
) -> Result<model::AminoAcidSet, Box<dyn std::error::Error>> {
    let mut fixed = db
        .fixed_mods
        .as_ref()
        .map_or_else(default_fixed_mods, |l| l.0.clone());
    if let (Some((name, mass)), None) = (tag, &db.fixed_mods) {
        for (residue, location) in [
            (ResidueSpec::Specific(b'K'), ModLocation::Anywhere),
            (ResidueSpec::Wildcard, ModLocation::NTerm),
        ] {
            fixed.push(Modification::new(name, mass, residue, location, true));
        }
    }
    let variable = db
        .variable_mods
        .as_ref()
        .map_or_else(default_variable_mods, |l| l.0.clone());
    Ok(AminoAcidSetBuilder::new_standard()
        .add_mods(fixed)
        .add_mods(variable)
        .build()?)
}

pub(crate) fn cli_fragmentation_to_activation_opt(f: Fragmentation) -> Option<ActivationMethod> {
    match f {
        Fragmentation::Auto => None,
        Fragmentation::Cid => Some(ActivationMethod::CID),
        Fragmentation::Etd => Some(ActivationMethod::ETD),
        Fragmentation::Hcd => Some(ActivationMethod::HCD),
        Fragmentation::Uvpd => Some(ActivationMethod::UVPD),
    }
}

pub(crate) fn cli_fragment_tol_override(
    fragment_tol_ppm: Option<f64>,
    fragment_tol_da: Option<f64>,
) -> Option<model::tolerance::Tolerance> {
    use model::tolerance::Tolerance;
    fragment_tol_ppm
        .map(Tolerance::Ppm)
        .or_else(|| fragment_tol_da.map(Tolerance::Da))
}

pub(crate) fn parse_enzymes(
    spec: &str,
) -> Result<(model::enzyme::Enzyme, Vec<model::enzyme::Enzyme>), Box<dyn std::error::Error>> {
    use model::enzyme::Enzyme;
    let mut all: Vec<Enzyme> = Vec::new();
    for tok in spec.split(',').map(str::trim).filter(|s| !s.is_empty()) {
        let e = Enzyme::from_name(tok).ok_or_else(|| {
            let hint = if tok.contains('+') {
                "; separate proteases with a comma, as in trypsin,lysc"
            } else {
                ""
            };
            format!(
                "unknown --enzyme '{tok}' (expected trypsin/chymotrypsin/lysc/aspn/gluc/lysn/argc/\
             alphalp/nocleavage/nonspecific/elastase){hint}"
            )
        })?;
        if !all.contains(&e) {
            all.push(e);
        }
    }
    all.sort_by_key(|&e| primary_rank(e));
    let (primary, extras) = all.split_first().ok_or("empty --enzyme")?;
    Ok((*primary, extras.to_vec()))
}

/// Sort key that puts a protease set in one order whatever order it was listed in; the
/// first one is the primary protease, which selects the scoring model and the cleavage
/// feature.
///
/// Trypsin leads any set that contains it: the scoring models are trained mostly on tryptic
/// digests, and with Lys-C, Arg-C, Glu-C or chymotrypsin beside it the digest is still
/// mostly tryptic. A set without trypsin puts first the protease whose cleavage rule sets
/// the C-terminal residue of most peptides of the union digest: a specific protease before
/// `NoCleavage`/`NonSpecific`, a C-terminal one (cleaving after its residues) before an
/// N-terminal one, then the one with more cleavage residues, a stand-in for the share of
/// peptide C-termini it makes. Ties fall back to the declaration order of
/// [`model::enzyme::Enzyme`].
fn primary_rank(e: model::enzyme::Enzyme) -> (bool, bool, bool, std::cmp::Reverse<usize>, u8) {
    use model::enzyme::Enzyme;
    let specific = !matches!(e, Enzyme::NoCleavage | Enzyme::NonSpecific);
    (
        e != Enzyme::Trypsin,
        !specific,
        !e.is_c_term(),
        std::cmp::Reverse(e.residues().len()),
        e as u8,
    )
}

pub(crate) fn warn_if_universal_protease_combo(
    primary: model::enzyme::Enzyme,
    extras: &[model::enzyme::Enzyme],
) {
    use model::enzyme::Enzyme;
    let is_universal = |e: &Enzyme| matches!(e, Enzyme::NonSpecific);
    let is_specific = |e: &Enzyme| !matches!(e, Enzyme::NonSpecific | Enzyme::NoCleavage);
    let all: Vec<Enzyme> = std::iter::once(primary)
        .chain(extras.iter().copied())
        .collect();
    if all.len() > 1 && all.iter().any(is_universal) && all.iter().any(is_specific) {
        log::warn!(
            target: "model",
            "--enzyme combines a universal protease (NonSpecific) with specific \
             one(s); their union cleaves at EVERY position, so this is effectively a \
             non-specific digest (much larger search space). Drop the universal protease to \
             keep the specific cleavage rules."
        );
    }
}

pub(crate) fn resolve_metadataless_selection(
    detected_activation: Option<ActivationMethod>,
    fragmentation: Fragmentation,
    fragment_tol_ppm: Option<f64>,
    fragment_tol_da: Option<f64>,
    inferred_high_res: Option<bool>,
) -> (ActivationMethod, Option<InstrumentType>) {
    let instrument: Option<InstrumentType> = if fragment_tol_ppm.is_some() {
        Some(InstrumentType::QExactive)
    } else if fragment_tol_da.is_some() {
        Some(InstrumentType::LowRes)
    } else {
        match inferred_high_res {
            Some(true) => Some(InstrumentType::QExactive),
            Some(false) => Some(InstrumentType::LowRes),
            None => None,
        }
    };
    let explicit = cli_fragmentation_to_activation_opt(fragmentation);
    let class_default = match instrument {
        Some(InstrumentType::QExactive)
        | Some(InstrumentType::HighRes)
        | Some(InstrumentType::TOF) => ActivationMethod::HCD,
        _ => ActivationMethod::CID,
    };
    if let (Some(e), Some(d)) = (explicit, detected_activation) {
        if e != d {
            log::info!(
                target: "model",
                "auto-detected activation = {}, OVERRIDDEN by --fragmentation {}",
                d.name(),
                e.name()
            );
        }
    }
    let activation = explicit.or(detected_activation).unwrap_or(class_default);
    if detected_activation.is_none() && explicit.is_none() && instrument.is_none() {
        log::warn!(
            target: "model",
            "no analyzer metadata, no --fragmentation/--fragment-tol, and the peaks \
             do not settle the fragment resolution; assuming CID / low-res / 0.5 Da. Pass \
             --fragmentation and --fragment-tol-ppm/-da to override."
        );
    }
    (activation, instrument)
}

pub(crate) fn bundled_store_path() -> PathBuf {
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

    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../resources/models")
}

pub(crate) fn build_selection_key(
    activation: ActivationMethod,
    instrument: Option<InstrumentType>,
    protocol: Protocol,
    enzyme: model::enzyme::Enzyme,
) -> SelectionKey {
    use std::collections::BTreeSet;

    let inst_after_family: &str = match instrument {
        Some(i) => i.family_fallback().name(),
        None => "LowRes",
    };
    let act_str: &str = match activation {
        ActivationMethod::PQD => "CID",
        ActivationMethod::EThcD if !matches!(inst_after_family, "QExactive" | "HighRes") => "HCD",
        other => other.name(),
    };

    let (final_act, final_inst, drop_protocol): (&str, &str, bool) =
        match (act_str, inst_after_family) {
            ("HCD", "LowRes") => ("CID", "LowRes", true),
            ("HCD" | "CID", "TOF") => ("HCD", "TOF", false),
            ("CID", "QExactive") => ("CID", "LowRes", true),
            ("ETD", i) if !matches!(i, "LowRes" | "HighRes") => ("ETD", "LowRes", true),
            ("UVPD", i) if i != "QExactive" => ("CID", "LowRes", true),
            _ => (act_str, inst_after_family, false),
        };

    let experiment_class: BTreeSet<String> = if drop_protocol {
        BTreeSet::new()
    } else {
        store_protocol_to_experiment_class(protocol_store_name(protocol))
    };

    SelectionKey {
        activation: final_act.to_string(),
        instrument: final_inst.to_string(),
        enzyme: routing_enzyme(enzyme).name().to_string(),
        experiment_class,
    }
}

pub(crate) fn routing_enzyme(enzyme: model::enzyme::Enzyme) -> model::enzyme::Enzyme {
    match enzyme {
        model::enzyme::Enzyme::NonSpecific => model::enzyme::Enzyme::NoCleavage,
        other => other,
    }
}

pub(crate) struct ModelChoice {
    pub(crate) model_id: String,
    pub(crate) key: SelectionKey,
    pub(crate) warning: Option<String>,
}

pub(crate) fn choose_model(
    entries: &[SelectionEntry],
    activation: ActivationMethod,
    instrument: Option<InstrumentType>,
    protocol: Protocol,
    enzyme: model::enzyme::Enzyme,
) -> ModelChoice {
    choose_routed(entries, activation, instrument, protocol, enzyme, true)
}

fn choose_routed(
    entries: &[SelectionEntry],
    activation: ActivationMethod,
    instrument: Option<InstrumentType>,
    protocol: Protocol,
    enzyme: model::enzyme::Enzyme,
    ion_trap_hcd: bool,
) -> ModelChoice {
    if activation == ActivationMethod::EThcD {
        let choice = choose_for_protocol(entries, activation, instrument, protocol, enzyme);
        let fits = entries
            .iter()
            .find(|e| e.model_id == choice.model_id)
            .is_some_and(|e| {
                e.activation == ActivationMethod::EThcD.name()
                    && e.enzyme == choice.key.enzyme
                    && e.experiment_class == choice.key.experiment_class
            });
        if !fits {
            let labelled = protocol == Protocol::Tmt;
            let nonspecific = choice.key.enzyme == model::enzyme::Enzyme::NoCleavage.name();
            if nonspecific || labelled {
                let class = if labelled {
                    store_protocol_to_experiment_class(protocol_store_name(Protocol::Auto))
                } else {
                    choice.key.experiment_class.clone()
                };
                if let Some(e) = entries.iter().find(|e| {
                    e.activation == ActivationMethod::EThcD.name()
                        && e.experiment_class == class
                        && (nonspecific || e.enzyme == choice.key.enzyme)
                        && analyser_of(&e.instrument) == analyser_of(&choice.key.instrument)
                }) {
                    return ModelChoice {
                        warning: untrained_route_warning(
                            entries,
                            &e.model_id,
                            activation,
                            instrument,
                            protocol,
                            enzyme,
                        ),
                        model_id: e.model_id.clone(),
                        key: choice.key,
                    };
                }
            }
            let hcd = choose_routed(
                entries,
                ActivationMethod::HCD,
                instrument,
                protocol,
                enzyme,
                false,
            );
            return ModelChoice {
                warning: untrained_route_warning(
                    entries,
                    &hcd.model_id,
                    activation,
                    instrument,
                    protocol,
                    enzyme,
                ),
                ..hcd
            };
        }
    }
    if ion_trap_hcd
        && activation == ActivationMethod::HCD
        && instrument == Some(InstrumentType::LowRes)
    {
        let enzyme_name = routing_enzyme(enzyme).name();
        let class = store_protocol_to_experiment_class(protocol_store_name(protocol));
        if let Some(e) = entries
            .iter()
            .filter(|e| {
                e.activation == ActivationMethod::HCD.name()
                    && e.instrument == InstrumentType::LowRes.name()
                    && e.enzyme == enzyme_name
                    && e.experiment_class == class
            })
            .min_by(|a, b| a.model_id.cmp(&b.model_id))
        {
            return ModelChoice {
                model_id: e.model_id.clone(),
                key: SelectionKey {
                    activation: e.activation.clone(),
                    instrument: e.instrument.clone(),
                    enzyme: e.enzyme.clone(),
                    experiment_class: class.clone(),
                },
                warning: None,
            };
        }
    }
    let choice = choose_for_protocol(entries, activation, instrument, protocol, enzyme);
    if protocol == Protocol::ItraqPhospho
        && class_of(entries, &choice.model_id) != Some(&choice.key.experiment_class)
    {
        let itraq = choose_for_protocol(entries, activation, instrument, Protocol::Itraq, enzyme);
        if class_of(entries, &itraq.model_id).is_some_and(|c| c.contains("itraq")) {
            return ModelChoice {
                warning: untrained_route_warning(
                    entries,
                    &itraq.model_id,
                    activation,
                    instrument,
                    protocol,
                    enzyme,
                ),
                model_id: itraq.model_id,
                key: choice.key,
            };
        }
    }
    ModelChoice {
        warning: untrained_route_warning(
            entries,
            &choice.model_id,
            activation,
            instrument,
            protocol,
            enzyme,
        ),
        ..choice
    }
}

pub(crate) fn analyser_of(instrument: &str) -> &str {
    match instrument {
        "LowRes" => "ion trap",
        "QExactive" | "HighRes" => "Orbitrap",
        "OrbitrapAstral" => "Astral",
        "TOF" | "TimsTOF" => "TOF",
        other => other,
    }
}

fn describe_cell(
    activation: &str,
    analyser: &str,
    enzyme: &str,
    class: &std::collections::BTreeSet<String>,
) -> String {
    let class = if class.is_empty() {
        "label-free".to_string()
    } else {
        class.iter().cloned().collect::<Vec<_>>().join("+")
    };
    format!("{activation} / {analyser} / {enzyme} / {class}")
}

pub(crate) fn untrained_route_warning(
    entries: &[SelectionEntry],
    model_id: &str,
    activation: ActivationMethod,
    instrument: Option<InstrumentType>,
    protocol: Protocol,
    enzyme: model::enzyme::Enzyme,
) -> Option<String> {
    let act = match activation {
        ActivationMethod::PQD => ActivationMethod::CID,
        other => other,
    }
    .name();
    let analyser = analyser_of(instrument.map_or("LowRes", |i| i.name()));
    let enzyme = routing_enzyme(enzyme).name();
    let class = store_protocol_to_experiment_class(protocol_store_name(protocol));
    let requested = describe_cell(act, analyser, enzyme, &class);
    let Some(e) = entries.iter().find(|e| e.model_id == model_id) else {
        return Some(format!(
            "no model in the store is trained for {requested}; \
             scoring with '{model_id}'. Scores may be mis-calibrated; pass --model to \
             choose explicitly."
        ));
    };
    let trained_analyser = analyser_of(&e.instrument);
    if e.activation == act
        && trained_analyser == analyser
        && e.enzyme == enzyme
        && e.experiment_class == class
    {
        return None;
    }
    let trained = describe_cell(
        &e.activation,
        trained_analyser,
        &e.enzyme,
        &e.experiment_class,
    );
    Some(format!(
        "no model in the store is trained for {requested}; scoring \
         with '{model_id}', trained on {trained}. Scores may be mis-calibrated; pass \
         --model to choose explicitly."
    ))
}

pub(crate) fn label_from_mods(aa: &model::AminoAcidSet) -> Option<input::IsobaricLabel> {
    use input::IsobaricLabel::{Itraq, Itraq8, Tmt, TmtPro};
    let on_label_site = |m: &Modification| {
        matches!(
            (m.residue, m.location),
            (ResidueSpec::Specific(b'K'), _)
                | (
                    ResidueSpec::Wildcard,
                    ModLocation::NTerm | ModLocation::ProtNTerm
                )
        )
    };
    [Tmt, TmtPro, Itraq, Itraq8].into_iter().find(|label| {
        let mass = label.tag().1;
        aa.distinct_mods()
            .iter()
            .any(|m| on_label_site(m) && (m.mass_delta - mass).abs() < 5e-4)
    })
}

pub(crate) fn routing_protocol(
    protocol: Protocol,
    label: Option<input::IsobaricLabel>,
) -> Protocol {
    match (protocol, label) {
        (Protocol::Auto, Some(input::IsobaricLabel::Tmt | input::IsobaricLabel::TmtPro)) => {
            Protocol::Tmt
        }
        (Protocol::Auto, Some(input::IsobaricLabel::Itraq | input::IsobaricLabel::Itraq8)) => {
            Protocol::Itraq
        }
        _ => protocol,
    }
}

fn protocol_store_name(protocol: Protocol) -> &'static str {
    match protocol {
        Protocol::Auto | Protocol::Standard => "Automatic",
        Protocol::Tmt => "TMT",
        Protocol::Phospho => "Phosphorylation",
        Protocol::Itraq => "iTRAQ",
        Protocol::ItraqPhospho => "iTRAQPhospho",
    }
}

fn class_of<'a>(
    entries: &'a [SelectionEntry],
    model_id: &str,
) -> Option<&'a std::collections::BTreeSet<String>> {
    entries
        .iter()
        .find(|e| e.model_id == model_id)
        .map(|e| &e.experiment_class)
}

fn choose_for_protocol(
    entries: &[SelectionEntry],
    activation: ActivationMethod,
    instrument: Option<InstrumentType>,
    protocol: Protocol,
    enzyme: model::enzyme::Enzyme,
) -> ModelChoice {
    let key = build_selection_key(activation, instrument, protocol, enzyme);
    let detected: Option<&str> = instrument
        .filter(|i| i.family_fallback() != *i)
        .map(|i| i.name());
    let (id, _substituted) = select_nearest_for(
        entries,
        &key,
        detected,
        |i| i.to_string(),
        "hcd_qexactive_tryp",
    );
    ModelChoice {
        model_id: id.to_string(),
        key,
        warning: None,
    }
}

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
        let choice = choose_model(&entries, activation, instrument, protocol, enzyme);
        if let Some(w) = &choice.warning {
            log::warn!(target: "model-routing", "{w}");
        }
        choice.model_id
    };

    let param = store
        .load_param(&model_id)
        .map_err(|e| format!("loading model '{model_id}' from store: {e}"))?;

    Ok((model_id, param))
}

pub(crate) fn load_activation_route(
    activation: ActivationMethod,
    instrument: Option<InstrumentType>,
    protocol: Protocol,
    enzyme: model::enzyme::Enzyme,
    store_path: &Path,
    run_model_id: &str,
    run_analyser: &str,
) -> Result<Option<(String, Param)>, Box<dyn std::error::Error>> {
    let store = ModelStore::open(store_path)
        .map_err(|e| format!("opening model store {}: {e}", store_path.display()))?;
    let entries = store.selection_entries();
    let choice = choose_model(&entries, activation, instrument, protocol, enzyme);
    if choice.model_id == run_model_id {
        return Ok(None);
    }
    let param = store
        .load_param(&choice.model_id)
        .map_err(|e| format!("loading model '{}' from store: {e}", choice.model_id))?;
    if analyser_of(param.data_type.instrument.name()) != run_analyser {
        return Ok(None);
    }
    if let Some(w) = &choice.warning {
        log::warn!(target: "model-routing", "{w}");
    }
    Ok(Some((choice.model_id, param)))
}

#[cfg(test)]
mod param_resolver_tests {
    use super::*;

    #[test]
    fn explicit_fragmentation_beats_a_detected_activation_without_instrument() {
        let (act, inst) = resolve_metadataless_selection(
            Some(ActivationMethod::HCD),
            Fragmentation::Etd,
            None,
            None,
            None,
        );
        assert_eq!((act, inst), (ActivationMethod::ETD, None));
        let (act, _) = resolve_metadataless_selection(
            Some(ActivationMethod::HCD),
            Fragmentation::Auto,
            None,
            None,
            None,
        );
        assert_eq!(act, ActivationMethod::HCD);
    }

    #[test]
    fn inferred_label_routes_like_the_explicit_protocol() {
        use input::IsobaricLabel::{Itraq, Itraq8, Tmt, TmtPro};
        let store = ModelStore::open(&bundled_store_path()).expect("bundled store");
        let entries = store.selection_entries();
        let labels = [
            (Tmt, Protocol::Tmt),
            (TmtPro, Protocol::Tmt),
            (Itraq, Protocol::Itraq),
            (Itraq8, Protocol::Itraq),
        ];
        for act in [
            ActivationMethod::HCD,
            ActivationMethod::CID,
            ActivationMethod::ETD,
            ActivationMethod::EThcD,
        ] {
            for inst in [
                None,
                Some(InstrumentType::LowRes),
                Some(InstrumentType::HighRes),
                Some(InstrumentType::QExactive),
                Some(InstrumentType::OrbitrapAstral),
            ] {
                for enz in [
                    model::enzyme::Enzyme::Trypsin,
                    model::enzyme::Enzyme::NoCleavage,
                ] {
                    for (label, explicit) in labels {
                        let routed = routing_protocol(Protocol::Auto, Some(label));
                        assert_eq!(routed, explicit);
                        let got = choose_model(&entries, act, inst, routed, enz);
                        let want = choose_model(&entries, act, inst, explicit, enz);
                        assert_eq!(
                            got.model_id, want.model_id,
                            "{act:?}/{inst:?}/{enz:?}/{label:?}"
                        );
                    }
                }
            }
        }
        let id = choose_model(
            &entries,
            ActivationMethod::EThcD,
            Some(InstrumentType::QExactive),
            routing_protocol(Protocol::Auto, Some(Tmt)),
            model::enzyme::Enzyme::Trypsin,
        )
        .model_id;
        assert_eq!(id, "ethcd_qexactive_tryp");
        let id = choose_model(
            &entries,
            ActivationMethod::EThcD,
            Some(InstrumentType::QExactive),
            Protocol::Tmt,
            model::enzyme::Enzyme::NonSpecific,
        )
        .model_id;
        assert_eq!(id, "ethcd_qexactive_tryp");
        let id = choose_model(
            &entries,
            ActivationMethod::EThcD,
            Some(InstrumentType::LowRes),
            Protocol::Tmt,
            model::enzyme::Enzyme::Trypsin,
        )
        .model_id;
        assert_eq!(
            id,
            choose_model(
                &entries,
                ActivationMethod::HCD,
                Some(InstrumentType::LowRes),
                Protocol::Tmt,
                model::enzyme::Enzyme::Trypsin,
            )
            .model_id
        );
        assert_eq!(
            routing_protocol(Protocol::Phospho, Some(Tmt)),
            Protocol::Phospho
        );
        assert_eq!(routing_protocol(Protocol::Auto, None), Protocol::Auto);
    }

    #[test]
    fn label_is_read_from_the_declared_mods() {
        use input::IsobaricLabel::{Itraq, Tmt, TmtPro};
        let set = |fixed: &str, variable: &str| {
            AminoAcidSetBuilder::new_standard()
                .add_mods(model::mod_spec::parse_mod_list(fixed, true).unwrap())
                .add_mods(model::mod_spec::parse_mod_list(variable, false).unwrap())
                .build()
                .unwrap()
        };
        assert_eq!(label_from_mods(&set("TMT6plex (K)", "")), Some(Tmt));
        assert_eq!(label_from_mods(&set("TMT6plex (N-term)", "")), Some(Tmt));
        assert_eq!(label_from_mods(&set("TMTpro (K)", "")), Some(TmtPro));
        assert_eq!(label_from_mods(&set("", "iTRAQ4plex (K)")), Some(Itraq));
        assert_eq!(
            label_from_mods(&set("Carbamidomethyl (C)", "Oxidation (M)")),
            None
        );
    }

    #[test]
    fn tof_runs_resolve_to_a_high_res_hcd_model() {
        let cases = [
            (Protocol::Auto, "hcd_qexactive_tryp"),
            (Protocol::Standard, "hcd_qexactive_tryp"),
            (Protocol::Tmt, "hcd_qexactive_tryp_tmt"),
            (Protocol::Phospho, "hcd_qexactive_tryp_phosphorylation"),
            (Protocol::Itraq, "hcd_qexactive_tryp_itraq"),
        ];
        for inst in [InstrumentType::TOF, InstrumentType::TimsTOF] {
            for act in [
                ActivationMethod::CID,
                ActivationMethod::HCD,
                ActivationMethod::PQD,
            ] {
                let key = build_selection_key(
                    act,
                    Some(inst),
                    Protocol::Auto,
                    model::enzyme::Enzyme::Trypsin,
                );
                assert_eq!(
                    (key.activation.as_str(), key.instrument.as_str()),
                    ("HCD", "TOF")
                );
                for (protocol, expected) in cases {
                    let (id, param) = load_param_from_store(
                        act,
                        Some(inst),
                        protocol,
                        model::enzyme::Enzyme::Trypsin,
                        None,
                        None,
                    )
                    .expect("bundled store");
                    assert_eq!(id, expected, "{act:?}/{inst:?}/{protocol:?}");
                    assert!(
                        param.data_type.instrument.is_high_resolution(),
                        "{act:?}/{inst:?}/{protocol:?} -> {id} is not high-res"
                    );
                }
            }
        }
    }
    #[test]
    fn a_non_specific_digest_routes_like_nocleavage() {
        use model::enzyme::Enzyme;
        let store = ModelStore::open(&bundled_store_path()).expect("bundled store");
        let entries = store.selection_entries();
        for name in ["nonspecific", "all", "elastase"] {
            let enzyme = Enzyme::from_name(name).expect("known enzyme name");
            for act in [ActivationMethod::HCD, ActivationMethod::CID] {
                for inst in [
                    Some(InstrumentType::QExactive),
                    Some(InstrumentType::HighRes),
                    Some(InstrumentType::OrbitrapAstral),
                    Some(InstrumentType::TOF),
                    Some(InstrumentType::LowRes),
                ] {
                    for prot in [Protocol::Auto, Protocol::Phospho, Protocol::Tmt] {
                        assert_eq!(
                            choose_model(&entries, act, inst, prot, enzyme).model_id,
                            choose_model(&entries, act, inst, prot, Enzyme::NoCleavage).model_id,
                            "{name}/{act:?}/{inst:?}/{prot:?}"
                        );
                    }
                }
            }
            let hla = choose_model(
                &entries,
                ActivationMethod::HCD,
                Some(InstrumentType::QExactive),
                Protocol::Auto,
                enzyme,
            );
            assert_eq!(hla.model_id, "hcd_highres_nocleavage", "{name}");
        }
    }

    #[test]
    fn itraq_phospho_takes_the_itraq_model() {
        use model::enzyme::Enzyme;
        let store = ModelStore::open(&bundled_store_path()).expect("bundled store");
        let entries = store.selection_entries();
        let pick = |act, inst, prot| choose_model(&entries, act, inst, prot, Enzyme::Trypsin);
        for inst in [InstrumentType::QExactive, InstrumentType::HighRes] {
            let c = pick(ActivationMethod::HCD, Some(inst), Protocol::ItraqPhospho);
            assert_eq!(c.model_id, "hcd_qexactive_tryp_itraq", "{inst:?}");
            assert!(c.warning.is_some(), "{inst:?}");
        }
        assert_eq!(
            pick(
                ActivationMethod::CID,
                Some(InstrumentType::LowRes),
                Protocol::ItraqPhospho
            )
            .model_id,
            pick(
                ActivationMethod::CID,
                Some(InstrumentType::LowRes),
                Protocol::Auto
            )
            .model_id
        );
    }

    #[test]
    fn mixed_tmt_run_scores_each_activation_with_its_own_model() {
        use model::enzyme::Enzyme;
        let store_path = bundled_store_path();
        let store = ModelStore::open(&store_path).expect("bundled store");
        let entries = store.selection_entries();
        let inst = Some(InstrumentType::QExactive);
        let run_of = |act| choose_model(&entries, act, inst, Protocol::Tmt, Enzyme::NonSpecific);
        let hcd = run_of(ActivationMethod::HCD).model_id;
        let ethcd = run_of(ActivationMethod::EThcD).model_id;
        assert_eq!(hcd, "hcd_highres_nocleavage");
        assert_eq!(ethcd, "ethcd_qexactive_tryp");
        for (scans, run_model) in [
            (ActivationMethod::EThcD, &hcd),
            (ActivationMethod::HCD, &ethcd),
        ] {
            let route = load_activation_route(
                scans,
                inst,
                Protocol::Tmt,
                Enzyme::NonSpecific,
                &store_path,
                run_model,
                "Orbitrap",
            )
            .expect("route")
            .expect("the scans of the other activation get a model of their own");
            assert_eq!(route.0, run_of(scans).model_id);
            assert_ne!(&route.0, run_model);
        }
    }

    #[test]
    fn bundled_store_selection_grid_is_pinned() {
        use model::enzyme::Enzyme;
        let store = ModelStore::open(&bundled_store_path()).expect("bundled store");
        let entries = store.selection_entries();
        let activations = [
            ActivationMethod::CID,
            ActivationMethod::ETD,
            ActivationMethod::HCD,
            ActivationMethod::UVPD,
            ActivationMethod::EThcD,
        ];
        let instruments = [
            None,
            Some(InstrumentType::LowRes),
            Some(InstrumentType::HighRes),
            Some(InstrumentType::TOF),
            Some(InstrumentType::QExactive),
            Some(InstrumentType::OrbitrapAstral),
            Some(InstrumentType::TimsTOF),
        ];
        let protocols = [
            ("auto", Protocol::Auto),
            ("phospho", Protocol::Phospho),
            ("iTRAQ", Protocol::Itraq),
            ("iTRAQ-phospho", Protocol::ItraqPhospho),
            ("TMT", Protocol::Tmt),
        ];
        let enzymes = [
            Enzyme::Trypsin,
            Enzyme::Chymotrypsin,
            Enzyme::LysC,
            Enzyme::AspN,
            Enzyme::GluC,
            Enzyme::LysN,
            Enzyme::ArgC,
            Enzyme::AlphaLP,
            Enzyme::NoCleavage,
            Enzyme::NonSpecific,
        ];
        let cell = |a, i, p, e| {
            let c = choose_model(&entries, a, i, p, e);
            format!(
                "{}{}",
                c.model_id,
                if c.warning.is_some() { "*" } else { "" }
            )
        };
        let mut table = String::from("activation\tinstrument\tenzyme");
        for (name, _) in protocols {
            table.push('\t');
            table.push_str(name);
        }
        table.push('\n');
        for act in activations {
            for inst in instruments {
                for enz in enzymes {
                    table.push_str(&format!(
                        "{}\t{}\t{}",
                        act.name(),
                        inst.map_or("none", |i| i.name()),
                        enz.name()
                    ));
                    for (_, prot) in protocols {
                        let got = cell(act, inst, prot, enz);
                        if act == ActivationMethod::CID {
                            assert_eq!(
                                cell(ActivationMethod::PQD, inst, prot, enz),
                                got,
                                "PQD/{inst:?}/{prot:?}/{enz:?} must resolve as CID"
                            );
                        }
                        if prot == Protocol::Auto {
                            assert_eq!(
                                cell(act, inst, Protocol::Standard, enz),
                                got,
                                "{act:?}/{inst:?}/standard/{enz:?} must resolve as auto"
                            );
                        }
                        table.push('\t');
                        table.push_str(&got);
                    }
                    table.push('\n');
                }
            }
        }
        let pinned = include_str!("../../../../../test-fixtures/model_selection_grid.tsv")
            .replace("\r\n", "\n");
        if table != pinned {
            let moved: Vec<String> = table
                .lines()
                .zip(pinned.lines())
                .filter(|(got, want)| got != want)
                .map(|(got, want)| format!("  pinned: {want}\n  now:    {got}"))
                .collect();
            panic!(
                "bundled-store selection differs from test-fixtures/model_selection_grid.tsv \
                 ({} row(s) moved):\n{}\n--- full table ---\n{table}",
                moved.len(),
                moved.join("\n")
            );
        }
    }

    #[test]
    fn ion_trap_hcd_takes_its_own_model_only_for_its_cell() {
        use model::enzyme::Enzyme;
        let store = ModelStore::open(&bundled_store_path()).expect("bundled store");
        let mut entries = store.selection_entries();
        let entry = |id: &str, act: &str, inst: &str, enz: &str, class: &[&str]| SelectionEntry {
            model_id: id.to_string(),
            activation: act.to_string(),
            instrument: inst.to_string(),
            enzyme: enz.to_string(),
            experiment_class: class.iter().map(|c| c.to_string()).collect(),
        };
        entries.retain(|e| !(e.activation == "HCD" && e.instrument == "LowRes"));
        let lowres = Some(InstrumentType::LowRes);
        let pick = |entries: &[SelectionEntry], act, prot, enz| {
            choose_model(entries, act, lowres, prot, enz)
        };
        let before = pick(
            &entries,
            ActivationMethod::HCD,
            Protocol::Auto,
            Enzyme::Trypsin,
        );
        assert_eq!(before.model_id, "cid_lowres_tryp");
        assert!(before.warning.is_some());

        entries.push(entry("hcd_lowres_probe", "HCD", "LowRes", "Trypsin", &[]));
        let hit = pick(
            &entries,
            ActivationMethod::HCD,
            Protocol::Auto,
            Enzyme::Trypsin,
        );
        assert_eq!(hit.model_id, "hcd_lowres_probe");
        assert!(hit.warning.is_none());
        let standard = pick(
            &entries,
            ActivationMethod::HCD,
            Protocol::Standard,
            Enzyme::Trypsin,
        );
        assert_eq!(standard.model_id, "hcd_lowres_probe");
        for (act, prot, enz) in [
            (ActivationMethod::HCD, Protocol::Auto, Enzyme::LysC),
            (ActivationMethod::HCD, Protocol::Phospho, Enzyme::Trypsin),
            (ActivationMethod::HCD, Protocol::Tmt, Enzyme::Trypsin),
            (ActivationMethod::CID, Protocol::Auto, Enzyme::Trypsin),
            (ActivationMethod::EThcD, Protocol::Auto, Enzyme::Trypsin),
        ] {
            assert_ne!(
                pick(&entries, act, prot, enz).model_id,
                "hcd_lowres_probe",
                "{act:?}/{prot:?}/{enz:?}"
            );
        }
        for inst in [None, Some(InstrumentType::QExactive)] {
            assert_ne!(
                choose_model(
                    &entries,
                    ActivationMethod::HCD,
                    inst,
                    Protocol::Auto,
                    Enzyme::Trypsin
                )
                .model_id,
                "hcd_lowres_probe"
            );
        }
    }

    #[test]
    fn lysc_and_trypsin_route_alike_in_either_order() {
        let store = ModelStore::open(&bundled_store_path()).expect("bundled store");
        let entries = store.selection_entries();
        let route = |spec: &str, act, inst, prot| {
            let (primary, _) = parse_enzymes(spec).unwrap();
            choose_model(&entries, act, inst, prot, primary)
        };
        for act in [
            ActivationMethod::CID,
            ActivationMethod::HCD,
            ActivationMethod::ETD,
            ActivationMethod::EThcD,
        ] {
            for inst in [
                None,
                Some(InstrumentType::LowRes),
                Some(InstrumentType::QExactive),
                Some(InstrumentType::OrbitrapAstral),
            ] {
                for prot in [Protocol::Auto, Protocol::Tmt] {
                    let lysc_first = route("lysc,trypsin", act, inst, prot);
                    let tryp_first = route("trypsin,lysc", act, inst, prot);
                    let trypsin = route("trypsin", act, inst, prot);
                    assert_eq!(
                        lysc_first.model_id, tryp_first.model_id,
                        "{act:?}/{inst:?}/{prot:?}"
                    );
                    assert_eq!(
                        lysc_first.model_id, trypsin.model_id,
                        "{act:?}/{inst:?}/{prot:?}"
                    );
                }
            }
        }
        let ion_trap = Some(InstrumentType::LowRes);
        let tmt = route(
            "lysc,trypsin",
            ActivationMethod::CID,
            ion_trap,
            Protocol::Tmt,
        );
        assert_eq!(
            (tmt.model_id.as_str(), tmt.warning.is_none()),
            ("cid_lowres_tryp_tmt", true)
        );
        let label_free = route(
            "lysc,trypsin",
            ActivationMethod::CID,
            ion_trap,
            Protocol::Auto,
        );
        assert_eq!(
            (label_free.model_id.as_str(), label_free.warning.is_none()),
            ("cid_lowres_tryp", true)
        );
    }

    #[test]
    fn metadataless_no_flags_default_selects_cid_lowres() {
        let store = ModelStore::open(&bundled_store_path()).expect("bundled store");
        let entries = store.selection_entries();
        let (act, inst) =
            resolve_metadataless_selection(None, Fragmentation::Auto, None, None, None);
        let c = choose_model(
            &entries,
            act,
            inst,
            Protocol::Auto,
            model::enzyme::Enzyme::Trypsin,
        );
        assert_eq!(
            (c.model_id.as_str(), c.warning.is_some()),
            ("cid_lowres_tryp", false)
        );
    }

    #[test]
    fn metadataless_peak_inference_drives_the_resolution_class() {
        let store = ModelStore::open(&bundled_store_path()).expect("bundled store");
        let entries = store.selection_entries();
        let pick = |high: Option<bool>, ppm: Option<f64>, da: Option<f64>| {
            let (act, inst) =
                resolve_metadataless_selection(None, Fragmentation::Auto, ppm, da, high);
            choose_model(
                &entries,
                act,
                inst,
                Protocol::Auto,
                model::enzyme::Enzyme::Trypsin,
            )
            .model_id
        };
        assert_eq!(pick(Some(true), None, None), "hcd_qexactive_tryp");
        assert_eq!(pick(Some(false), None, None), "cid_lowres_tryp");
        assert_eq!(pick(None, None, None), "cid_lowres_tryp");
        assert_eq!(pick(Some(true), None, Some(0.5)), "cid_lowres_tryp");
        assert_eq!(pick(Some(false), Some(20.0), None), "hcd_qexactive_tryp");
    }

    use crate::cli::{
        parse_charge_range, parse_enum, parse_isotope_error_range, parse_precursor_cal,
        parse_precursor_tol, EnzymeSpecificity,
    };
    use crate::spectra::title_prefix_for;
    use ::search::PrecursorCalMode;
    use model::Tolerance;

    #[test]
    fn config_enum_values_are_named() {
        let frag = |s: &str| parse_enum::<Fragmentation>("fragmentation", s, true);
        assert_eq!(frag("HCD").unwrap(), Fragmentation::Hcd);
        assert_eq!(frag("hcd").unwrap(), Fragmentation::Hcd);
        assert_eq!(frag("auto").unwrap(), Fragmentation::Auto);
        assert!(frag("3").is_err());
        let protocol = |s: &str| parse_enum::<Protocol>("protocol", s, true);
        assert_eq!(protocol("TMT").unwrap(), Protocol::Tmt);
        assert_eq!(protocol("itraq-PHOSPHO").unwrap(), Protocol::ItraqPhospho);
        assert!(protocol("4").is_err());
        let ntt = |s: &str| parse_enum::<EnzymeSpecificity>("enzyme_specificity", s, true);
        assert_eq!(ntt("fully").unwrap(), EnzymeSpecificity::Fully);
        assert_eq!(ntt("Semi").unwrap(), EnzymeSpecificity::Semi);
        assert!(ntt("2").is_err());
        assert_eq!(
            frag("x").unwrap_err(),
            "invalid fragmentation `x`: expected auto|CID|ETD|HCD|UVPD"
        );
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
    fn every_search_prefixes_titles_with_file_stem() {
        assert_eq!(title_prefix_for("myfile").as_deref(), Some("myfile/"));
    }
}
