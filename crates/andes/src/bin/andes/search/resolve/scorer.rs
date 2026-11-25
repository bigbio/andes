//! The scoring model of the run (bundled, trained or given), the retrieval mode and the
//! model-dependent settings (GBDT truncation, isobaric peak filter).

use std::path::PathBuf;
use std::sync::Arc;

use model::activation::ActivationMethod;
use model::InstrumentType;
use scoring_crate::RankScorer;

use super::{RunInputs, ScoringModel};
use crate::cli::{Cli, FragmentIndexFlag, Protocol};
use crate::models::{
    cli_fragment_tol_override, cli_fragmentation_to_activation_opt, load_param_from_store,
    resolve_metadataless_selection, search_aa_set,
};

/// Chooses the scoring model for the detected activation and instrument, and builds the scorer.
pub(crate) fn resolve_scoring_model(
    cli: &mut Cli,
    inputs: &RunInputs,
    aa: &mut model::AminoAcidSet,
    search_enzyme: model::enzyme::Enzyme,
) -> Result<ScoringModel, Box<dyn std::error::Error>> {
    let spectrum_path: PathBuf = cli.spectrum[0].clone();
    let detected_activation = inputs.detected_activation;
    let detected_instrument = inputs.detected_instrument;
    let instrument_was_detected = inputs.instrument_was_detected;
    let acquisition = &inputs.acquisition;
    let label = inputs.label;
    let route_protocol = inputs.route_protocol;
    let t_phase = std::time::Instant::now();
    let (mut param, routed_activation, routed_instrument, resolved_model_id) = {
        let detected = (detected_activation, detected_instrument);
        let (activation, instrument_opt): (ActivationMethod, Option<InstrumentType>) =
            match detected {
                (Some(method), Some(inst)) => {
                    let chosen =
                        cli_fragmentation_to_activation_opt(cli.fragmentation).unwrap_or(method);
                    if chosen != method {
                        log::info!(
                            target: "model",
                            "auto-detected activation = {} (instrument = {}) from {}, \
                             OVERRIDDEN by --fragmentation {}",
                            method.name(), inst.name(), spectrum_path.display(), chosen.name()
                        );
                    } else {
                        log::info!(
                            target: "model",
                            "auto-detected activation = {} (instrument = {}) from {}",
                            method.name(), inst.name(), spectrum_path.display()
                        );
                    }
                    (chosen, Some(inst))
                }
                (None, Some(inst)) => {
                    let explicit = cli_fragmentation_to_activation_opt(cli.fragmentation);
                    let chosen = explicit.unwrap_or(if inst.is_high_resolution() {
                        ActivationMethod::HCD
                    } else {
                        ActivationMethod::CID
                    });
                    log::info!(
                        target: "model",
                        "instrument = {} from {}, no activation recorded; \
                         activation = {} [{}]",
                        inst.name(),
                        spectrum_path.display(),
                        chosen.name(),
                        if explicit.is_some() {
                            "--fragmentation"
                        } else {
                            "default of the analyzer class"
                        }
                    );
                    (chosen, Some(inst))
                }
                (method, None) => resolve_metadataless_selection(
                    method,
                    cli.fragmentation,
                    cli.fragment_tol_ppm,
                    cli.fragment_tol_da,
                    acquisition.inferred_high_res(),
                ),
            };

        let (model_id, p) = load_param_from_store(
            activation,
            instrument_opt,
            route_protocol,
            search_enzyme,
            cli.model_store.as_deref(),
            cli.model_id_override.as_deref(),
        )?;
        log::info!(target: "model", "{model_id} (from store)");
        (p, activation, instrument_opt, model_id)
    };
    let gbdt_k = cli.gbdt_max_trees.unwrap_or(if cli.glyco {
        0
    } else {
        crate::cli::DEFAULT_GBDT_MAX_TREES
    });
    if gbdt_k > 0 {
        for (slot, name) in [
            (&mut param.frag_intensity_model, "frag-intensity"),
            (&mut param.rich_ion_model, "rich-ion"),
        ] {
            if let Some(before) = truncate_gbdt(slot, gbdt_k) {
                log::info!(target: "search",
                    "--gbdt-max-trees {gbdt_k}: {name} ensemble {before} -> {gbdt_k} trees \
                     (changes predictions; not byte-identical)"
                );
            }
        }
    }
    if cli.no_search_gbdt {
        param.frag_intensity_model = None;
        param.rich_ion_model = None;
        log::info!(target: "search", "--no-search-gbdt: search scoring without the fragment-intensity and rich-ion ensembles");
    }
    match cli.protocol {
        Protocol::Tmt => param.data_type.protocol = model::protocol::Protocol::TMT,
        Protocol::Itraq => param.data_type.protocol = model::protocol::Protocol::ITRAQ,
        Protocol::ItraqPhospho => {
            param.data_type.protocol = model::protocol::Protocol::ITRAQPhospho
        }
        Protocol::Auto => {
            if let Some(l) = label {
                log::info!(
                    target: "protocol",
                    "{} {l:?} → isobaric protocol (windowed peak filter on low-res MS2)",
                    if cli.tmt.is_some() {
                        "--tmt"
                    } else {
                        "detected"
                    }
                );
                param.data_type.protocol = match l {
                    input::IsobaricLabel::Tmt | input::IsobaricLabel::TmtPro => {
                        model::protocol::Protocol::TMT
                    }
                    input::IsobaricLabel::Itraq | input::IsobaricLabel::Itraq8 => {
                        model::protocol::Protocol::ITRAQ
                    }
                };
            }
        }
        _ => {}
    }
    if cli.db.fixed_mods.is_none() {
        let tag = match (param.data_type.protocol, label) {
            (
                model::protocol::Protocol::TMT,
                Some(l @ (input::IsobaricLabel::Tmt | input::IsobaricLabel::TmtPro)),
            ) => Some(l.tag()),
            (model::protocol::Protocol::TMT, _) => Some(input::IsobaricLabel::Tmt.tag()),
            (
                model::protocol::Protocol::ITRAQ | model::protocol::Protocol::ITRAQPhospho,
                Some(l @ (input::IsobaricLabel::Itraq | input::IsobaricLabel::Itraq8)),
            ) => Some(l.tag()),
            (model::protocol::Protocol::ITRAQ | model::protocol::Protocol::ITRAQPhospho, _) => {
                Some(input::IsobaricLabel::Itraq.tag())
            }
            (_, Some(l)) if cli.tmt.is_some() => Some(l.tag()),
            _ => None,
        };
        if let Some((name, mass)) = tag {
            *aa = search_aa_set(&cli.db, Some((name, mass)))?;
            log::info!(
                target: "protocol",
                "injected {name} fixed mod (+{mass:.4} on K + peptide N-term) \
                 into the candidate set (no --fixed-mods given)"
            );
        }
    }
    let mut scorer = RankScorer::new(&param);
    scorer.set_acquisition_high_res(Some(acquisition.high_res));
    let dia_tol_default = cli.dia
        && acquisition.analyzer == Some(InstrumentType::OrbitrapAstral)
        && cli.fragment_tol_ppm.is_none()
        && cli.fragment_tol_da.is_none();
    if dia_tol_default {
        cli.fragment_tol_ppm = Some(crate::dia::FRAGMENT_TOL_PPM);
        log::info!(target: "search",
            "--dia: fragment tolerance {} ppm (Orbitrap Astral MS2)",
            crate::dia::FRAGMENT_TOL_PPM
        );
    }
    let frag_tol_override = cli_fragment_tol_override(cli.fragment_tol_ppm, cli.fragment_tol_da);
    if let Some(t) = frag_tol_override {
        if instrument_was_detected && !dia_tol_default {
            log::info!(target: "search",
                "--fragment-tol-* {t:?} overrides the feature tolerance of the detected {} analyzer",
                acquisition.analyzer.map_or("", |i| i.name())
            );
        }
        scorer.set_fragment_tol_override(frag_tol_override);
    }
    log::info!(
        target: "param_and_scorer",
        "done in {:.2} s",
        t_phase.elapsed().as_secs_f64()
    );
    Ok(ScoringModel {
        param,
        scorer,
        routed_activation,
        routed_instrument,
        resolved_model_id,
        gbdt_k,
        frag_tol_override,
    })
}

pub(crate) struct RetrievalChoice {
    pub use_index: bool,
    pub report: Option<String>,
    pub refused_because: Option<&'static str>,
}

pub(crate) fn retrieval_choice(
    fragment_tol: model::tolerance::Tolerance,
    excluded_because: Option<&'static str>,
    flag: FragmentIndexFlag,
) -> RetrievalChoice {
    use model::tolerance::Tolerance;
    let tol_text = match fragment_tol {
        Tolerance::Ppm(v) => format!("{v} ppm"),
        Tolerance::Da(v) => format!("{v} Da"),
    };
    if let Some(reason) = excluded_because {
        return RetrievalChoice {
            use_index: false,
            report: None,
            refused_because: (flag == FragmentIndexFlag::On).then_some(reason),
        };
    }
    match flag {
        FragmentIndexFlag::Off => RetrievalChoice {
            use_index: false,
            report: Some("out-of-core index; the fragment-ion index was switched off".to_string()),
            refused_because: None,
        },
        FragmentIndexFlag::Auto | FragmentIndexFlag::On => {
            let report = match search::fragments::index::vote_peaks_for(fragment_tol) {
                0 => format!("out-of-core index, high-resolution fragments at {tol_text}"),
                n => format!(
                    "out-of-core index, low-resolution fragments at {tol_text}; \
                     the {n} most intense peaks vote"
                ),
            };
            RetrievalChoice {
                use_index: true,
                report: Some(report),
                refused_because: None,
            }
        }
    }
}

/// Shortens the GBDT ensemble in `slot` to `k` trees and returns its former tree count, or
/// `None` when it is absent or already no longer than `k`.
pub(crate) fn truncate_gbdt(
    slot: &mut Option<Arc<scoring_crate::gbdt::GbdtPeakModel>>,
    k: usize,
) -> Option<usize> {
    let g = slot.as_ref()?;
    let before = g.trees.len();
    if k >= before {
        return None;
    }
    let mut t = (**g).clone();
    t.trees.truncate(k);
    *slot = Some(Arc::new(t));
    Some(before)
}

pub(super) fn isobaric_peak_filter(
    explicit: Option<(f64, usize)>,
    high_res: bool,
) -> Option<(f64, usize)> {
    explicit.or(high_res.then_some((0.0, 0)))
}
