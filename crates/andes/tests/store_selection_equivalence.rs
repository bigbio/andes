//! Equivalence gate: for every (ActivationMethod × InstrumentType × Protocol)
//! combination that the reference resolution ladder handles, assert that the
//! store-based selection returns the same `model_id`.
//!
//! This test is the safety proof that the store-based model selection
//! (`ModelStore::load_param`) returns the correct bundled model for each
//! activation/instrument/protocol combination.

use std::collections::BTreeSet;
use std::path::PathBuf;

use model::{activation::ActivationMethod, InstrumentType};
use model_train::{
    ModelStore,
    select::{select, SelectionEntry, SelectionKey},
};

// ── helpers (mirrors the search binary's ladder) ─────────────────────────────

/// CLI Protocol enum (mirrors the binary's `Protocol`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Protocol { Auto, Phospho, Itraq, ItraqPhospho, Tmt, Standard }

/// Fragmentation enum (mirrors the binary's `Fragmentation`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Fragmentation { Auto, Cid, Etd, Hcd, Uvpd }

/// Instrument CLI enum (mirrors the binary's `Instrument`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Instrument { LowRes, HighRes, Tof, QExactive }

fn protocol_suffix(p: Protocol) -> &'static str {
    match p {
        Protocol::Auto | Protocol::Standard => "",
        Protocol::Phospho      => "_Phosphorylation",
        Protocol::Itraq        => "_iTRAQ",
        Protocol::ItraqPhospho => "_iTRAQPhospho",
        Protocol::Tmt          => "_TMT",
    }
}

/// Reference 5-level resolution ladder, expressed directly as lowercase
/// `model_id` strings (the canonical identifiers in the parquet store).
fn resolve_model_id_old(
    fragmentation: Fragmentation,
    instrument: Instrument,
    protocol: Protocol,
) -> String {
    // Step 0: all-defaults short-circuit.
    if fragmentation == Fragmentation::Auto
        && instrument == Instrument::LowRes
        && protocol == Protocol::Auto
    {
        return "hcd_qexactive_tryp".to_string();
    }

    let frag = match fragmentation {
        Fragmentation::Auto => "CID",
        Fragmentation::Cid  => "CID",
        Fragmentation::Etd  => "ETD",
        Fragmentation::Hcd  => "HCD",
        Fragmentation::Uvpd => "UVPD",
    };
    let inst = match instrument {
        Instrument::LowRes    => "LowRes",
        Instrument::HighRes   => "HighRes",
        Instrument::Tof       => "TOF",
        Instrument::QExactive => "QExactive",
    };
    // H5: low-res (ion-trap) HCD is NOT routed to the high-res model. With no
    // hcd_lowres model bundled, frag="HCD"/inst="LowRes" finds no exact match
    // and falls through to the cid_lowres_tryp final-fallback below — the
    // correct low-res b/y model.

    let prot_suffix = protocol_suffix(protocol);
    let exact = model_id(frag, inst, prot_suffix);
    // The exact protocol-suffixed model is only selectable if it is bundled AND
    // actually keyed for the requested protocol's experiment_class — this
    // mirrors `select()`'s exact-match step, which is protocol-key-aware.
    // (e.g. `etd_highres_tryp_phosphorylation` is keyed protocol=Phosphorylation
    // in the v1 (N=19) store, so ETD/HighRes/Phospho exact-matches it rather
    // than falling back to the base model. This check reads the store's actual
    // experiment_class, so it tracks the manifest protocol key automatically.)
    if bundled_with_protocol(&exact, protocol) { return exact; }

    if !prot_suffix.is_empty() {
        let no_prot = model_id(frag, inst, "");
        if let Some(id) = try_bundled(&no_prot) { return id; }
    }

    // Final fallback ladder. Each fallback target is gated through the bundle:
    // the v1 bundle ships a curated subset (19 own-trained regimes), so a
    // fallback model that is not bundled (e.g. `cid_tof_tryp`, `etd_lowres_tryp`
    // were dropped as un-sourceable) degrades to the global default
    // `hcd_qexactive_tryp` — exactly what the real `select()` does.
    //
    // IMPORTANT: this mirror must match the binary's `build_selection_key`
    // *normalization*, not just bundle membership. The binary only rewrites
    // `(HCD,LowRes)`, `(CID,QExactive)`, `(UVPD,non-QE)`, `(ETD,non-LowRes/HighRes)`
    // and `(HCD,TOF)`; for every other unmatched `(frag,inst)` it keeps the
    // instrument verbatim, so a TOF instrument with no `*_tof_*` model bundled
    // (the v1 case) falls straight through to the global default rather than to
    // `cid_lowres_tryp`. Encode that here so the reference tracks the binary.
    let final_fallback: Option<&str> = match (frag, inst) {
        // TOF/HighRes HCD maps to cid_tof_tryp.
        ("HCD", "TOF") | ("HCD", "HighRes") => Some("cid_tof_tryp"),
        ("ETD", _)                          => Some("etd_lowres_tryp"),
        // The binary rewrites these to cid_lowres_tryp (drop_protocol arms).
        ("CID", "QExactive") | ("UVPD", _)  => Some("cid_lowres_tryp"),
        // CID/HCD on LowRes resolve to the low-res b/y model.
        (_, "LowRes")                       => Some("cid_lowres_tryp"),
        // Anything else (notably CID/PQD on TOF or QExactive-family without a
        // matching model) has no normalization arm → global default.
        _                                   => None,
    };
    match final_fallback {
        Some(id) => try_bundled(id).unwrap_or_else(|| "hcd_qexactive_tryp".to_string()),
        None     => "hcd_qexactive_tryp".to_string(),
    }
}

/// Build a lowercase store `model_id` from the (fragmentation, instrument,
/// protocol-suffix) components, e.g. `("CID","LowRes","_TMT")` → `cid_lowres_tryp_tmt`.
fn model_id(frag: &str, inst: &str, prot_suffix: &str) -> String {
    format!("{frag}_{inst}_Tryp{prot_suffix}").to_ascii_lowercase()
}

/// Reference resolution ladder keyed on (activation, instrument, protocol).
fn resolve_for_activation_old(
    method: ActivationMethod,
    detected_instrument: Option<InstrumentType>,
    protocol: Protocol,
) -> String {
    let frag = match method {
        ActivationMethod::CID  => Fragmentation::Cid,
        ActivationMethod::ETD  => Fragmentation::Etd,
        ActivationMethod::HCD  => Fragmentation::Hcd,
        ActivationMethod::UVPD => Fragmentation::Uvpd,
        ActivationMethod::PQD  => Fragmentation::Cid,
    };
    let inst = match detected_instrument.map(|i| i.family_fallback()) {
        Some(InstrumentType::LowRes)         => Instrument::LowRes,
        Some(InstrumentType::HighRes)        => Instrument::HighRes,
        Some(InstrumentType::TOF)            => Instrument::Tof,
        Some(InstrumentType::QExactive)      => Instrument::QExactive,
        Some(InstrumentType::OrbitrapAstral) => Instrument::QExactive,
        Some(InstrumentType::TimsTOF)        => Instrument::Tof,
        None                                 => Instrument::LowRes,
    };
    resolve_model_id_old(frag, inst, protocol)
}

/// Lazily initialized set of model IDs from the bundled parquet store.
/// Opened once and reused across all `try_bundled` calls in the test.
fn bundled_model_ids() -> &'static std::collections::BTreeSet<String> {
    use std::sync::OnceLock;
    static IDS: OnceLock<std::collections::BTreeSet<String>> = OnceLock::new();
    IDS.get_or_init(|| {
        let store_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../resources/models.parquet");
        let store = model_train::ModelStore::open(&store_path)
            .unwrap_or_else(|e| panic!("failed to open bundled models.parquet: {e}"));
        store.model_ids().into_iter().collect()
    })
}

/// Lazily initialized map `model_id -> experiment_class` from the bundled store.
fn bundled_model_classes() -> &'static std::collections::BTreeMap<String, BTreeSet<String>> {
    use std::sync::OnceLock;
    static MAP: OnceLock<std::collections::BTreeMap<String, BTreeSet<String>>> = OnceLock::new();
    MAP.get_or_init(|| {
        let store_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../resources/models.parquet");
        let store = model_train::ModelStore::open(&store_path)
            .unwrap_or_else(|e| panic!("failed to open bundled models.parquet: {e}"));
        store
            .selection_entries()
            .into_iter()
            .map(|e| (e.model_id, e.experiment_class))
            .collect()
    })
}

/// True iff `model_id` is bundled AND its experiment_class matches the
/// requested protocol (mirrors `select()`'s protocol-aware exact match).
fn bundled_with_protocol(model_id: &str, protocol: Protocol) -> bool {
    let want = protocol_to_experiment_class(protocol);
    match bundled_model_classes().get(model_id) {
        Some(have) => *have == want,
        None => false,
    }
}

/// Return the `model_id` iff it is present in the bundled parquet store.
fn try_bundled(model_id: &str) -> Option<String> {
    if bundled_model_ids().contains(model_id) {
        Some(model_id.to_string())
    } else {
        None
    }
}

// ── new selection helpers ────────────────────────────────────────────────────

/// Convert the CLI Protocol to an experiment_class BTreeSet for the SelectionKey.
///
/// `ItraqPhospho` uses `{"itraqphospho"}` as a single opaque slug (matching
/// how the parquet stores `iTRAQPhospho`), NOT `{"itraq","phospho"}`, so that
/// the exact-match step finds the combo models and the empty-class fallback
/// still fires when no combo model is bundled.
fn protocol_to_experiment_class(p: Protocol) -> BTreeSet<String> {
    fn s(v: &str) -> String { v.to_string() }
    match p {
        Protocol::Auto | Protocol::Standard => BTreeSet::new(),
        Protocol::Tmt          => [s("tmt")].into(),
        Protocol::Phospho      => [s("phospho")].into(),
        Protocol::Itraq        => [s("itraq")].into(),
        // Keep as single "itraqphospho" slug to match the parquet row
        // and avoid spurious phospho-subset matches (step 2) when the
        // combo file is not bundled.
        Protocol::ItraqPhospho => [s("itraqphospho")].into(),
    }
}

/// Apply the same normalization the old ladder performs on (activation, instrument)
/// AFTER the instrument family fallback (`OrbitrapAstral`→`QExactive`,
/// `TimsTOF`→`TOF`) has already been applied.
///
/// Returns `(final_activation, final_instrument, drop_protocol)`.
/// `drop_protocol = true` when this is a "final fallback" normalization that
/// switches to a different activation/instrument class — in that case the old
/// ladder ignores the protocol and returns the base model, so the caller must
/// also clear the experiment_class.
///
/// `drop_protocol = false` means this is the HCD+LowRes upgrade (same model
/// family, different instrument slug only), where the protocol IS preserved.
fn normalize_activation_instrument(act: &str, inst: &str) -> (String, String, bool) {
    match (act, inst) {
        // H5: HCD + LowRes → CID + LowRes. Low-res (ion-trap) HCD is not routed
        // to the high-res model; there is no hcd_lowres model, so route to the
        // low-res b/y model (cid_lowres_tryp). Protocol dropped (final-fallback
        // base model). Mirrored in build_selection_key.
        ("HCD", "LowRes") => ("CID".into(), "LowRes".into(), true),
        // HCD + TOF → CID + TOF (no hcd_tof model bundled; final fallback
        // maps (HCD, TOF|HighRes) → cid_tof_tryp).
        // Protocol dropped: final fallback returns base model only.
        ("HCD", "TOF") => ("CID".into(), "TOF".into(), true),
        // CID + QExactive → CID + LowRes (no cid_qexactive model bundled;
        // final fallback default arm → cid_lowres_tryp).
        // Protocol dropped.
        ("CID", "QExactive") => ("CID".into(), "LowRes".into(), true),
        // ETD + any non-(LowRes|HighRes) → ETD + LowRes (final fallback
        // `("ETD", _)` → etd_lowres_tryp). Protocol dropped.
        ("ETD", i) if !matches!(i, "LowRes" | "HighRes") => {
            ("ETD".into(), "LowRes".into(), true)
        }
        // UVPD + non-QExactive → CID + LowRes (only uvpd_qexactive_tryp
        // is bundled; final fallback default arm → cid_lowres_tryp).
        // Protocol dropped.
        ("UVPD", i) if i != "QExactive" => ("CID".into(), "LowRes".into(), true),
        _ => (act.into(), inst.into(), false),
    }
}

/// Build the SelectionKey from (ActivationMethod, InstrumentType, Protocol),
/// applying all old-ladder normalizations so that `select()` with an identity
/// family_fn performs a direct lookup.
fn build_key(
    method: ActivationMethod,
    instrument: InstrumentType,
    protocol: Protocol,
) -> SelectionKey {
    // 1. PQD → CID (PQD is scored with the CID model).
    let act = match method {
        ActivationMethod::PQD => "CID",
        other                 => other.name(),
    };
    // 2. Apply family fallback (OrbitrapAstral→QExactive, TimsTOF→TOF).
    let inst = instrument.family_fallback().name();
    // 3. Apply old-ladder instrument+activation normalization.
    //    `drop_protocol` is true when this is a final-fallback normalization
    //    that changes the activation/instrument class — in that case the old
    //    ladder returns the base model (no protocol), so we clear the class.
    let (final_act, final_inst, drop_protocol) = normalize_activation_instrument(act, inst);
    let experiment_class = if drop_protocol {
        BTreeSet::new()
    } else {
        protocol_to_experiment_class(protocol)
    };

    SelectionKey {
        activation: final_act,
        instrument: final_inst,
        // Parquet stores enzyme as "Trypsin"; use the same string in the key.
        enzyme: "Trypsin".into(),
        experiment_class,
    }
}

/// The instrument_family closure used by `select()`.
/// Since all normalization is pre-applied in `build_key`, this is identity.
fn instrument_family(inst: &str) -> String {
    inst.to_string()
}

// ── the full matrix test ─────────────────────────────────────────────────────

fn all_activations() -> Vec<ActivationMethod> {
    vec![
        ActivationMethod::CID,
        ActivationMethod::ETD,
        ActivationMethod::HCD,
        ActivationMethod::PQD,
        ActivationMethod::UVPD,
    ]
}

fn all_instruments() -> Vec<InstrumentType> {
    vec![
        InstrumentType::LowRes,
        InstrumentType::HighRes,
        InstrumentType::TOF,
        InstrumentType::QExactive,
        InstrumentType::OrbitrapAstral,
        InstrumentType::TimsTOF,
    ]
}

fn all_protocols() -> Vec<Protocol> {
    vec![
        Protocol::Auto,
        Protocol::Phospho,
        Protocol::Itraq,
        Protocol::ItraqPhospho,
        Protocol::Tmt,
        Protocol::Standard,
    ]
}

/// Open the bundled parquet store and return its selection entries.
/// Note: `ItraqPhospho` entries use experiment_class `{"itraqphospho"}` so
/// that exact-match in select() works for the bundled combo models.
fn bundled_selection_entries() -> Vec<SelectionEntry> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../resources/models.parquet");
    let store = ModelStore::open(&path).expect("failed to open bundled models.parquet");
    store.selection_entries()
}

#[test]
fn store_selection_matches_old_ladder_for_all_combos() {
    let entries = bundled_selection_entries();

    let mut failures: Vec<String> = Vec::new();

    for &act in &all_activations() {
        for &inst in &all_instruments() {
            for &prot in &all_protocols() {
                let old_id = resolve_for_activation_old(act, Some(inst), prot);

                let key = build_key(act, inst, prot);
                let new_id = select(
                    &entries,
                    &key,
                    instrument_family,
                    Some("hcd_qexactive_tryp"),
                )
                .map(|s| s.to_string())
                .unwrap_or_else(|| "hcd_qexactive_tryp".to_string());

                if new_id != old_id {
                    failures.push(format!(
                        "{:?}/{:?}/{:?}: old={old_id} new={new_id}  key=({},{},{},[{:?}])",
                        act, inst, prot,
                        key.activation, key.instrument, key.enzyme,
                        key.experiment_class.iter().collect::<Vec<_>>()
                    ));
                }
            }
        }
    }

    // Also test the `None` instrument (no instrument detected → LowRes).
    for &act in &all_activations() {
        for &prot in &all_protocols() {
            let old_id = resolve_for_activation_old(act, None, prot);

            // None instrument → LowRes (low-res default).
            let key = build_key(act, InstrumentType::LowRes, prot);
            let new_id = select(
                &entries,
                &key,
                instrument_family,
                Some("hcd_qexactive_tryp"),
            )
            .map(|s| s.to_string())
            .unwrap_or_else(|| "hcd_qexactive_tryp".to_string());

            if new_id != old_id {
                failures.push(format!(
                    "{:?}/None/{:?}: old={old_id} new={new_id}  key=({},{},{},[{:?}])",
                    act, prot,
                    key.activation, key.instrument, key.enzyme,
                    key.experiment_class.iter().collect::<Vec<_>>()
                ));
            }
        }
    }

    if !failures.is_empty() {
        panic!(
            "store selection diverges from old ladder for {} combo(s):\n{}",
            failures.len(),
            failures.join("\n")
        );
    }
}

// ── Metadata-less CLI default ─────────────────────────────────────────────────
//
// The equivalence matrix above exercises the *activation-aware* ladder
// (`resolve_for_activation_old`), which always receives a concrete activation
// method, so it never hits the all-defaults short-circuit.
//
// With no analyzer metadata and no `--fragmentation`/`--fragment-tol-*`, the
// binary's `resolve_metadataless_selection` yields `(CID, None)` →
// `cid_lowres_tryp`. This test pins that default against the store directly
// (mirroring the binary's resolver).
#[test]
fn metadataless_no_flags_default_selects_cid_lowres() {
    let entries = bundled_selection_entries();

    // Mirror the binary's `resolve_metadataless_selection` for the no-flags
    // case: no detected activation, Fragmentation::Auto, no fragment-tol.
    // → activation = CID, instrument = None (→ LowRes via the empty-instrument
    //   normalization), protocol = Auto.
    let key = build_key(ActivationMethod::CID, InstrumentType::LowRes, Protocol::Auto);
    let new_id = select(
        &entries,
        &key,
        instrument_family,
        Some("hcd_qexactive_tryp"),
    )
    .expect("cid_lowres_tryp must be present in the bundled store")
    .to_string();

    assert_eq!(
        new_id, "cid_lowres_tryp",
        "decision E: metadata-less no-flags default must resolve to cid_lowres_tryp \
         (not the old hcd_qexactive)"
    );
}
