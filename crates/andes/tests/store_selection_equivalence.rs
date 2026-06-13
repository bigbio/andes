//! Equivalence gate: for every (ActivationMethod × InstrumentType × Protocol)
//! combination that the old `resolve_bundled_param_for_activation` ladder
//! handles, assert that the new store-based selection returns the same
//! `model_id` as the lowercased filename stem of the old path.
//!
//! This test is the safety proof that switching the search binary from
//! `Param::load_from_file` to `ModelStore::load_param` is behavior-preserving
//! for the bundled store.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

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

/// Replicate the OLD 5-level ladder (resolve_bundled_param logic).
fn resolve_bundled_param_old(
    fragmentation: Fragmentation,
    instrument: Instrument,
    protocol: Protocol,
) -> PathBuf {
    // Step 0: all-defaults short-circuit.
    if fragmentation == Fragmentation::Auto
        && instrument == Instrument::LowRes
        && protocol == Protocol::Auto
    {
        return canonicalize_bundled("HCD_QExactive_Tryp.param");
    }

    let frag = match fragmentation {
        Fragmentation::Auto => "CID",
        Fragmentation::Cid  => "CID",
        Fragmentation::Etd  => "ETD",
        Fragmentation::Hcd  => "HCD",
        Fragmentation::Uvpd => "UVPD",
    };
    let mut inst = match instrument {
        Instrument::LowRes    => "LowRes",
        Instrument::HighRes   => "HighRes",
        Instrument::Tof       => "TOF",
        Instrument::QExactive => "QExactive",
    };
    // HCD-upgrade rule.
    if frag == "HCD" && inst == "LowRes" {
        inst = "QExactive";
    }

    let prot_suffix = protocol_suffix(protocol);
    let exact = format!("{frag}_{inst}_Tryp{prot_suffix}.param");
    if let Ok(p) = try_bundled(&exact) { return p; }

    if !prot_suffix.is_empty() {
        let no_prot = format!("{frag}_{inst}_Tryp.param");
        if let Ok(p) = try_bundled(&no_prot) { return p; }
    }

    // Final fallback ladder.
    let final_fallback = match (frag, inst) {
        ("HCD", "TOF") | ("HCD", "HighRes") => "CID_TOF_Tryp.param",
        ("ETD", _)                           => "ETD_LowRes_Tryp.param",
        _                                    => "CID_LowRes_Tryp.param",
    };
    canonicalize_bundled(final_fallback)
}

/// Replicate the OLD resolve_bundled_param_for_activation logic.
fn resolve_for_activation_old(
    method: ActivationMethod,
    detected_instrument: Option<InstrumentType>,
    protocol: Protocol,
) -> PathBuf {
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
    resolve_bundled_param_old(frag, inst, protocol)
}

/// Build a path under resources/ionstat for a given filename WITHOUT requiring
/// the file to exist on disk. Used only to derive the lowercased stem for
/// comparison with the parquet store's model IDs — the .param files themselves
/// are no longer shipped on disk (they live in models.parquet).
fn canonicalize_bundled(filename: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../resources/ionstat")
        .join(filename)
}

/// Lazily initialized set of model IDs from the bundled parquet store.
/// Opened once and reused across all `try_bundled` calls in the test.
fn bundled_model_ids() -> &'static std::collections::BTreeSet<String> {
    use std::sync::OnceLock;
    static IDS: OnceLock<std::collections::BTreeSet<String>> = OnceLock::new();
    IDS.get_or_init(|| {
        let store_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../resources/ionstat/models.parquet");
        let store = model_train::ModelStore::open(&store_path)
            .unwrap_or_else(|e| panic!("failed to open bundled models.parquet: {e}"));
        store.model_ids().into_iter().collect()
    })
}

/// Check whether a bundled .param file WOULD have existed under the old naming
/// scheme. Since the files no longer live on disk, we derive the expected
/// existence from the bundled parquet store: a model exists iff it has an
/// entry in the store.
fn try_bundled(filename: &str) -> Result<PathBuf, ()> {
    // Derive the expected model ID from the filename stem (lowercase).
    let stem = PathBuf::from(filename)
        .file_stem()
        .and_then(|s| s.to_str())
        .map(|s| s.to_ascii_lowercase())
        .unwrap_or_default();
    if bundled_model_ids().contains(&stem) {
        Ok(canonicalize_bundled(filename))
    } else {
        Err(())
    }
}

/// Extract the lowercased filename stem from a `.param` path.
/// e.g. `.../HCD_QExactive_Tryp.param` → `"hcd_qexactive_tryp"`.
fn stem_of(p: &Path) -> String {
    p.file_stem()
        .and_then(|s| s.to_str())
        .map(|s| s.to_ascii_lowercase())
        .unwrap_or_default()
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
        // HCD + LowRes → HCD + QExactive (Java's instrument-upgrade rule).
        // Protocol is preserved: HCD_QExactive_Tryp_TMT etc. exist.
        ("HCD", "LowRes") => ("HCD".into(), "QExactive".into(), false),
        // HCD + TOF → CID + TOF (HCD_TOF_Tryp.param not bundled; old
        // final fallback maps (HCD, TOF|HighRes) → CID_TOF_Tryp).
        // Protocol dropped: final fallback returns base model only.
        ("HCD", "TOF") => ("CID".into(), "TOF".into(), true),
        // CID + QExactive → CID + LowRes (CID_QExactive_Tryp not bundled;
        // old final fallback default arm → CID_LowRes_Tryp).
        // Protocol dropped.
        ("CID", "QExactive") => ("CID".into(), "LowRes".into(), true),
        // ETD + any non-(LowRes|HighRes) → ETD + LowRes (old final fallback
        // `("ETD", _)` → ETD_LowRes_Tryp.param). Protocol dropped.
        ("ETD", i) if !matches!(i, "LowRes" | "HighRes") => {
            ("ETD".into(), "LowRes".into(), true)
        }
        // UVPD + non-QExactive → CID + LowRes (only UVPD_QExactive_Tryp
        // is bundled; old final fallback default arm → CID_LowRes_Tryp).
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
    // 1. PQD → CID (Java's NewScorerFactory rule).
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
        .join("../../resources/ionstat/models.parquet");
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
                let old_path = resolve_for_activation_old(act, Some(inst), prot);
                let old_stem = stem_of(&old_path);

                let key = build_key(act, inst, prot);
                let new_id = select(
                    &entries,
                    &key,
                    instrument_family,
                    Some("hcd_qexactive_tryp"),
                )
                .map(|s| s.to_string())
                .unwrap_or_else(|| "hcd_qexactive_tryp".to_string());

                if new_id != old_stem {
                    failures.push(format!(
                        "{:?}/{:?}/{:?}: old={old_stem} new={new_id}  key=({},{},{},[{:?}])",
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
            let old_path = resolve_for_activation_old(act, None, prot);
            let old_stem = stem_of(&old_path);

            // None instrument → LowRes (Java default LOW_RESOLUTION_LTQ).
            let key = build_key(act, InstrumentType::LowRes, prot);
            let new_id = select(
                &entries,
                &key,
                instrument_family,
                Some("hcd_qexactive_tryp"),
            )
            .map(|s| s.to_string())
            .unwrap_or_else(|| "hcd_qexactive_tryp".to_string());

            if new_id != old_stem {
                failures.push(format!(
                    "{:?}/None/{:?}: old={old_stem} new={new_id}  key=({},{},{},[{:?}])",
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

// ── Decision E: metadata-less CLI default ─────────────────────────────────────
//
// The equivalence matrix above exercises the *activation-aware* ladder
// (`resolve_for_activation_old`), which always receives a concrete activation
// method, so it never hits the historical all-defaults short-circuit
// (`Fragmentation::Auto && Instrument::LowRes && Protocol::Auto → hcd_qexactive`)
// that the removed `--instrument` CLI flag used to reach. That short-circuit
// was a CLI-flag artifact, not a property of the store, so dropping it does not
// affect the matrix above.
//
// Decision E changes that metadata-less, no-flags CLI default from
// `hcd_qexactive` to `cid_lowres`: with no analyzer metadata and no
// `--fragmentation`/`--fragment-tol-*`, the binary's
// `resolve_metadataless_selection` now yields `(CID, None)` → `cid_lowres_tryp`.
// This test pins that new default against the store directly (mirroring the
// binary's resolver), so the behavior change is asserted, not merely implied.
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
