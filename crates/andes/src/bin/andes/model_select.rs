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
    select::{select_nearest_for, SelectionEntry, SelectionKey},
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
                annotated_site_mass: None,
                blocks_cleavage: false,
                parts: Vec::new(),
                diagnostic_ions: Vec::new(),
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
                annotated_site_mass: None,
                blocks_cleavage: false,
                parts: Vec::new(),
                diagnostic_ions: Vec::new(),
            };
            Ok(AminoAcidSetBuilder::new_standard()
                .add_fixed_mod(cam)
                .add_variable_mod(ox)
                .build()?)
        }
    }
}

/// Build the default `AminoAcidSet` (CAM-C fixed, Ox-M variable, protein-N-term
/// Acetyl variable), optionally with an isobaric tag (TMT/iTRAQ) as a FIXED mod on K + peptide N-term.
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
        annotated_site_mass: None,
        blocks_cleavage: false,
        parts: Vec::new(),
        diagnostic_ions: Vec::new(),
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
        annotated_site_mass: None,
        blocks_cleavage: false,
        parts: Vec::new(),
        diagnostic_ions: Vec::new(),
    };
    // Protein N-terminal acetylation (+42.010565) — a near-universal default in
    // the field (MSFragger / Comet / MaxQuant). Restricted to the PROTEIN N-term
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
        annotated_site_mass: None,
        blocks_cleavage: false,
        parts: Vec::new(),
        diagnostic_ions: Vec::new(),
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
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
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
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
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
/// Accepts `,` or `+` separators (the latter matching MaxQuant/MSFragger
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

/// Warn when a universal protease (NonSpecific, which cleaves every bond) is
/// combined with specific protease(s): the cleavage-site union then
/// covers EVERY position, so the digest is effectively non-specific — usually
/// not what a user adding it to a specific list intends.
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
        eprintln!(
            "WARN: --enzyme combines a universal protease (NonSpecific) with specific \
             one(s); their union cleaves at EVERY position, so this is effectively a \
             non-specific digest (much larger search space). Drop the universal protease to \
             keep the specific cleavage rules."
        );
    }
}

/// Resolve (activation, instrument) for model selection on metadata-less input
/// (MGF, or mzML/.raw with no analyzer metadata). Resolution class comes from
/// the `--fragment-tol-*` unit, else from `inferred_high_res` (the class read
/// from the peaks, see `input::infer_fragment_resolution`); activation from an
/// explicit `--fragmentation`, else the detected method, else the
/// class-consistent default (the same precedence as a fully detected input).
/// When nothing disambiguates, defaults to CID / LowRes (→ `cid_lowres_tryp`)
/// + a warning.
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
    // Class-consistent default when neither detection nor `--fragmentation`
    // names an activation: high-res classes imply HCD, otherwise CID.
    let class_default = match instrument {
        Some(InstrumentType::QExactive)
        | Some(InstrumentType::HighRes)
        | Some(InstrumentType::TOF) => ActivationMethod::HCD,
        _ => ActivationMethod::CID,
    };
    if let (Some(e), Some(d)) = (explicit, detected_activation) {
        if e != d {
            eprintln!(
                "Param resolver: auto-detected activation = {}, OVERRIDDEN by --fragmentation {}",
                d.name(),
                e.name()
            );
        }
    }
    let activation = explicit.or(detected_activation).unwrap_or(class_default);
    if detected_activation.is_none() && explicit.is_none() && instrument.is_none() {
        eprintln!(
            "WARN: no analyzer metadata, no --fragmentation/--fragment-tol, and the peaks \
             do not settle the fragment resolution; assuming CID / low-res / 0.5 Da. Pass \
             --fragmentation and --fragment-tol-ppm/-da to override."
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

    // 2. Apply family fallback (OrbitrapAstral → QExactive, TimsTOF → TOF).
    let inst_after_family: &str = match instrument {
        Some(i) => i.family_fallback().name(),
        None => "LowRes",
    };
    // 1. PQD → CID for model routing. EThcD keys its own model on an Orbitrap
    //    read-out; read elsewhere (ion trap, TOF, undetected) it routes like HCD,
    //    the b/y half of its spectrum, as it did before the EThcD model existed.
    let act_str: &str = match activation {
        ActivationMethod::PQD => "CID",
        ActivationMethod::EThcD if !matches!(inst_after_family, "QExactive" | "HighRes") => "HCD",
        other => other.name(),
    };

    // 3. Apply old-ladder (activation, instrument) normalization.
    //    Because `normalize_for_store` returns `&'static str` only for the
    //    normalizing arms (to avoid lifetime issues), we handle identity
    //    inline here.
    let (final_act, final_inst, drop_protocol): (&str, &str, bool) =
        match (act_str, inst_after_family) {
            // H5: HCD read out in the ion trap. A store without a model trained
            // on these spectra routes the key to the low-res CID model of the
            // digest (`cid_lowres_tryp` for trypsin) and the run warns (see
            // `untrained_route_warning`); `choose_model` takes an ion-trap HCD
            // model first when the store has one for the digest and labelling.
            // The fragment tolerance does not depend on this choice: on ion-trap
            // fragments every bundled model scores its rank tables in the stored
            // 0.5 Da window (only high-res runs learn a narrower one, see
            // `RankScorer::node_score_window`), and the PIN features follow the
            // run's acquisition class. The Orbitrap HCD model is no substitute:
            // against it, PSMs at 1 % entrapment FDP were within -9 % to +7 %
            // across four tryptic label-free ion-trap HCD datasets, two in each
            // direction.
            ("HCD", "LowRes") => ("CID", "LowRes", true),
            // TOF analyzers (QTOF, timsTOF) fragment by beam-type collisional
            // activation: the b/y chemistry of HCD, measured at high resolution.
            // Vendors label it "CID" (the timsTOF `.d` reader and Bruker mzML
            // exports both do), so CID/PQD and HCD on TOF share one key, routed
            // to an HCD model. With no TOF-keyed model in the store,
            // `select_nearest` then substitutes the same-resolution high-res HCD
            // model, never the 0.5-Da ion-trap model that the activation alone
            // would pick. The protocol is kept: a TMT/phospho TOF run gets the
            // matching high-res protocol model.
            ("HCD" | "CID", "TOF") => ("HCD", "TOF", false),
            ("CID", "QExactive") => ("CID", "LowRes", true),
            ("ETD", i) if !matches!(i, "LowRes" | "HighRes") => ("ETD", "LowRes", true),
            ("UVPD", i) if i != "QExactive" => ("CID", "LowRes", true),
            _ => (act_str, inst_after_family, false),
        };

    // 4. Build experiment_class from protocol (unless the final fallback dropped it).
    //    Protocol → experiment_class mapping matches the parquet's `protocol` column.
    let experiment_class: BTreeSet<String> = if drop_protocol {
        BTreeSet::new()
    } else {
        store_protocol_to_experiment_class(protocol_store_name(protocol))
    };

    SelectionKey {
        activation: final_act.to_string(),
        instrument: final_inst.to_string(),
        // Parquet stores the enzyme as its `Enzyme::name()` ("Trypsin", "LysC", ...).
        enzyme: routing_enzyme(enzyme).name().to_string(),
        experiment_class,
    }
}

/// The enzyme a digest is routed by. A non-specific digest (`--enzyme
/// nonspecific` and its aliases `all` and `elastase`) is the immunopeptidomics
/// and peptidomics search the no-cleavage models were trained on (their
/// training searches were unspecific HLA ligandomes), so it selects with the
/// `NoCleavage` key, exactly as `--enzyme nocleavage` does.
pub(crate) fn routing_enzyme(enzyme: model::enzyme::Enzyme) -> model::enzyme::Enzyme {
    match enzyme {
        model::enzyme::Enzyme::NonSpecific => model::enzyme::Enzyme::NoCleavage,
        other => other,
    }
}

/// The model [`load_param_from_store`] picks when no `--model` is given.
pub(crate) struct ModelChoice {
    pub(crate) model_id: String,
    /// The key the store was searched with (after the family fallback and the
    /// activation/instrument normalizations of [`build_selection_key`]).
    pub(crate) key: SelectionKey,
    /// Set when the model was not trained on the requested cell (activation,
    /// analyser, enzyme, labelling/enrichment): the one warning the run prints,
    /// naming the cell and the model used. See [`untrained_route_warning`].
    pub(crate) warning: Option<String>,
}

/// Pick the model for `(activation, instrument, protocol, enzyme)` from the
/// store's selection `entries`: [`build_selection_key`], then
/// [`select_nearest_for`] with the detected instrument when the key collapsed
/// it to its family. `build_selection_key` already applies the family fallback
/// and every normalization, so the family function is the identity. When the
/// exact ladder misses (e.g. a protocol the own-only store doesn't carry),
/// `select_nearest_for` routes to the CLOSEST own model and only resolves to
/// the standard base as a last resort. This is the selection the search and
/// the training paths run, and the one the bundled-store grid test pins.
pub(crate) fn choose_model(
    entries: &[SelectionEntry],
    activation: ActivationMethod,
    instrument: Option<InstrumentType>,
    protocol: Protocol,
    enzyme: model::enzyme::Enzyme,
) -> ModelChoice {
    choose_routed(entries, activation, instrument, protocol, enzyme, true)
}

/// [`choose_model`]; `ion_trap_hcd` allows the ion-trap HCD model for HCD read out in the
/// ion trap. An EThcD run read in the ion trap routes as HCD without it: that model is
/// trained and measured on pure HCD scans only.
fn choose_routed(
    entries: &[SelectionEntry],
    activation: ActivationMethod,
    instrument: Option<InstrumentType>,
    protocol: Protocol,
    enzyme: model::enzyme::Enzyme,
    ion_trap_hcd: bool,
) -> ModelChoice {
    // EThcD: the EThcD model when the store has one trained for this cell's enzyme and
    // experiment class (labelling, enrichment). Otherwise the run is routed as HCD, the
    // b/y half of its spectrum, as before the EThcD model existed: the labelled,
    // enriched and enzyme-specific HCD models stay the choice for the cells no EThcD
    // model was trained or measured on, and a store without an EThcD model keeps the HCD
    // route.
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
            // Two cells have no EThcD model of their own and take the tryptic, label-free
            // EThcD model anyway, because it scores the c/z and b/y series of their
            // spectra, which the HCD models cannot:
            // - a non-specific digest (HLA ligands, peptidomics): against the no-cleavage
            //   HCD route, +39 % PSMs and +30 % peptides at 1 % entrapment FDP on HLA
            //   class I EThcD runs;
            // - a TMT run (iTRAQ and phospho are not measured and keep the HCD route): on
            //   TMT10 EThcD scans of a non-tryptic plasma peptidome, +23 % PSMs and +22 %
            //   peptides at 1 % entrapment FDP against the labelled HCD route (12 of 12
            //   draw/seed pairs).
            // The run warns, as every route to a cell the model was not trained on.
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
    // HCD read out in the ion trap: the model trained on those spectra, when the store has
    // one for this digest and experiment class. A cell with no such model keeps the low-res
    // CID route of `build_selection_key`.
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
    // iTRAQ-phospho: without a model trained on that class, the ladder falls to
    // the unlabelled base. The iTRAQ model is the closer one: the tag on every
    // peptide N-terminus and lysine shifts and enhances the b-ion series and adds
    // reporter and tag fragments on every spectrum, while the phosphate touches
    // only the modified peptides (its neutral loss is scored from the
    // modification, not from the table). TMT-phospho already resolves to the
    // TMT model by the same preference (the labelling tag wins the subset
    // tie-break in `select`). The iTRAQ model is taken only when it exists for
    // this activation and instrument; otherwise the first choice stands.
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

/// The fragment analyser an instrument class measures with. The selection
/// keys keep the store's instrument names; a model trained on one Orbitrap
/// (`QExactive`, `HighRes`) is trained for every Orbitrap, and the TOF names
/// share one analyser too.
pub(crate) fn analyser_of(instrument: &str) -> &str {
    match instrument {
        "LowRes" => "ion trap",
        "QExactive" | "HighRes" => "Orbitrap",
        "OrbitrapAstral" => "Astral",
        "TOF" | "TimsTOF" => "TOF",
        other => other,
    }
}

/// `activation / analyser / enzyme / labelling` for a warning.
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

/// The warning for a run scored with a model that was not trained on its cell:
/// the requested (activation, analyser, enzyme, labelling/enrichment) against
/// the cell the store records for `model_id`. `None` when they agree. Every
/// route to another cell warns, whether the ladder substituted the nearest
/// model, a normalization changed the activation or analyser (ion-trap HCD,
/// Orbitrap CID, UVPD off the Orbitrap), the protocol fell to a model without
/// it, or a deliberate route (iTRAQ-phospho to the iTRAQ model) was taken.
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
            "WARN: model routing: no model in the store is trained for {requested}; \
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
        "WARN: model routing: no model in the store is trained for {requested}; scoring \
         with '{model_id}', trained on {trained}. Scores may be mis-calibrated; pass \
         --model to choose explicitly."
    ))
}

/// The isobaric label the declared modifications carry: a mod of a label's
/// mass (to 0.0005 Da) on K or on the peptide N-terminus, fixed or variable.
/// `None` for a set without one.
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

/// The protocol the model is chosen for. An explicit `--protocol` stands; at
/// `auto` a run that carries a label (named by `--tmt`, read from its
/// reporter ions, or declared by its modifications) routes like the explicit
/// protocol of that label, so the two spellings of one run agree on every
/// activation.
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

/// The store's `protocol` name for a CLI protocol.
fn protocol_store_name(protocol: Protocol) -> &'static str {
    match protocol {
        Protocol::Auto | Protocol::Standard => "Automatic",
        Protocol::Tmt => "TMT",
        Protocol::Phospho => "Phosphorylation",
        Protocol::Itraq => "iTRAQ",
        Protocol::ItraqPhospho => "iTRAQPhospho",
    }
}

/// The experiment class the store's `model_id` was trained on.
fn class_of<'a>(
    entries: &'a [SelectionEntry],
    model_id: &str,
) -> Option<&'a std::collections::BTreeSet<String>> {
    entries
        .iter()
        .find(|e| e.model_id == model_id)
        .map(|e| &e.experiment_class)
}

/// [`choose_model`] for one protocol: [`build_selection_key`], then
/// [`select_nearest_for`].
fn choose_for_protocol(
    entries: &[SelectionEntry],
    activation: ActivationMethod,
    instrument: Option<InstrumentType>,
    protocol: Protocol,
    enzyme: model::enzyme::Enzyme,
) -> ModelChoice {
    let key = build_selection_key(activation, instrument, protocol, enzyme);
    // A model trained for the EXACT detected instrument (e.g. a user-trained
    // OrbitrapAstral model) is hidden by the family collapse, so the detected
    // instrument is passed along when the collapse changed it.
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

/// Load the scoring [`Param`] from the bundled Parquet store for the given
/// `(activation, instrument, protocol)` combination.
///
/// The `model_id` is chosen by [`choose_model`]; the bundled-store selection
/// grid test pins its result for every activation × instrument × protocol ×
/// enzyme combination.
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
        let choice = choose_model(&entries, activation, instrument, protocol, enzyme);
        if let Some(w) = &choice.warning {
            eprintln!("{w}");
        }
        choice.model_id
    };

    let param = store
        .load_param(&model_id)
        .map_err(|e| format!("loading model '{model_id}' from store: {e}"))?;

    Ok((model_id, param))
}

/// The model for the scans of a run that were recorded with `activation`, when the
/// run itself was routed for another activation (a decision-tree or alternating
/// acquisition). Same rules as the run-level choice ([`choose_model`] with the
/// run's instrument, protocol and enzyme); the untrained-route warning prints
/// once, here, because the route is built once per activation.
///
/// `None` keeps those scans on the run's model: the choice is the model the run
/// already uses, or it was trained on another analyser than the run's
/// (`run_analyser`: the scans' instrument is the run's, so a model for a
/// different analyser cannot describe them).
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
        eprintln!("{w}");
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

    /// A label named by the run (`--tmt`, reporter ions or modifications)
    /// routes exactly as the explicit `--protocol` of that label: the model
    /// depends on the run, not on how the label was spelled. EThcD, HCD and
    /// CID, every analyser, every enzyme.
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
        // Labelled Orbitrap EThcD takes the EThcD model (measured on TMT10 EThcD scans);
        // the same run read in the ion trap keeps the labelled HCD route.
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
        // An explicit protocol is never overridden by the label.
        assert_eq!(
            routing_protocol(Protocol::Phospho, Some(Tmt)),
            Protocol::Phospho
        );
        assert_eq!(routing_protocol(Protocol::Auto, None), Protocol::Auto);
    }

    /// A label declared by the modifications (fixed or variable, on K or the
    /// peptide N-terminus) is read as that label; other mods are not.
    #[test]
    fn label_is_read_from_the_declared_mods() {
        use input::IsobaricLabel::{Itraq, Tmt, TmtPro};
        let set = |lines: &[&str]| {
            let mut b = AminoAcidSetBuilder::new_standard();
            for l in lines {
                let m = Modification::from_mods_txt_line(l).unwrap();
                b = if m.fixed {
                    b.add_fixed_mod(m)
                } else {
                    b.add_variable_mod(m)
                };
            }
            b.build().unwrap()
        };
        assert_eq!(
            label_from_mods(&set(&["229.162932,K,fix,any,TMT"])),
            Some(Tmt)
        );
        assert_eq!(
            label_from_mods(&set(&["229.162932,*,fix,N-term,TMT"])),
            Some(Tmt)
        );
        assert_eq!(
            label_from_mods(&set(&["304.207146,K,fix,any,TMTpro"])),
            Some(TmtPro)
        );
        assert_eq!(
            label_from_mods(&set(&["144.102063,K,opt,any,iTRAQ"])),
            Some(Itraq)
        );
        assert_eq!(
            label_from_mods(&set(&[
                "57.021464,C,fix,any,Carbamidomethyl",
                "15.994915,M,opt,any,Oxidation"
            ])),
            None
        );
    }

    /// TOF and timsTOF runs fragment by beam-type CID, whether the input labels
    /// it CID (the `.d` reader, Bruker mzML) or HCD. With no TOF model in the
    /// bundled store every combination resolves to a high-resolution HCD model
    /// for the requested protocol, never to the 0.5-Da ion-trap model.
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
    /// A non-specific digest, under every name it accepts, selects exactly the
    /// model `--enzyme nocleavage` selects: on high-res HCD that is the
    /// no-cleavage (HLA) model, never the tryptic base.
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

    /// iTRAQ-phospho with no model trained on it takes the iTRAQ model of the
    /// same activation and instrument (and warns), not the unlabelled base;
    /// where no iTRAQ model exists the first choice stands.
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

    /// A mixed HCD/EThcD TMT run on an Orbitrap, searched non-specifically (PXD004781):
    /// whichever activation the run is routed for, the scans of the other one get
    /// their own model, so the run's choice only decides which of the two is the
    /// run-level model. Each direction of the per-scan route is a model other than the
    /// run's, on the run's analyser.
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

    /// The bundled-store selection, as the search and training paths run it
    /// (`choose_model`: `build_selection_key`, the detected-instrument pre-pass
    /// and `select_nearest_for`), over every activation × detected instrument
    /// (or none) × protocol × enzyme, pinned in
    /// `test-fixtures/model_selection_grid.tsv`. One row per (activation,
    /// instrument, enzyme), one column per protocol; `*` marks a cell the
    /// routed model was not trained on (the run warns, `untrained_route_warning`). PQD must resolve as CID and `standard` as
    /// `auto`, so neither has its own row or column. A change to the store, the
    /// key or the ladder that moves any selection fails here with the new
    /// table; update the fixture only when every moved cell is intended.
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
        // A Windows checkout may give the fixture CRLF line endings.
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

    /// HCD read out in the ion trap takes an ion-trap HCD model of the digest and experiment
    /// class when the store has one, and the low-res CID route otherwise.
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
        // Other digests, labelling, enrichment and activations keep their routes.
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
        // An undetected analyser and the Orbitrap read-out never take it.
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

    /// With no analyzer metadata and no `--fragmentation`/`--fragment-tol-*`,
    /// the resolver yields `(CID, None)`, which selects the low-res b/y model.
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

    /// Metadata-less input whose peaks read as high resolution selects the
    /// high-res HCD model; a low-resolution call keeps the low-res b/y model.
    /// The `--fragment-tol-*` unit still wins over the inferred class.
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

    // ── Model resolution is store-based: all bundled models live in
    //    resources/models/ and are selected by `model_id`.
    //    `bundled_store_selection_grid_is_pinned` covers the
    //    activation/instrument/protocol/enzyme → model_id selection.

    /// The `--config` values of the named enums: their CLI names, in any case;
    /// numbers are not accepted.
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
        // One input file is prefixed too: per-file searches are pooled into one
        // Percolator run, and native titles repeat across files.
        assert_eq!(title_prefix_for("myfile").as_deref(), Some("myfile/"));
    }
}
