use model::mass::{H2O, NH3, PROTON, Z_DOT_OFFSET};
use scoring_crate::params::IonType;
use scoring_crate::psm::rank::RankScorer;
use std::path::Path;

const KEY_TOLERANCE_DA: f64 = 1e-6;

fn offset_of(ion: &IonType) -> Option<f64> {
    ion.offset().map(f64::from)
}

fn is_ion(ion: &IonType, prefix: bool, offset: f64) -> bool {
    let kind_matches = match ion {
        IonType::Prefix {
            charge: 1,
            loss_class: 0,
            ..
        } => prefix,
        IonType::Suffix {
            charge: 1,
            loss_class: 0,
            ..
        } => !prefix,
        _ => false,
    };
    kind_matches && offset_of(ion).is_some_and(|o| (o - offset).abs() < KEY_TOLERANCE_DA)
}

/// The shipped models carry the f32 offset keys written at training time. Search reads each
/// key from the model, so a change of the physical constants must still find the singly
/// charged b and y (or c and z) tables in every shipped model, and the scorer must hold a
/// table for them.
#[test]
fn shipped_models_keep_their_primary_ion_tables() {
    let bundled = Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../resources/models"
    ));
    let store = model_train::store::ModelStore::open(bundled).expect("open bundled model store");
    let ids = store.model_ids();
    assert!(!ids.is_empty());

    for id in ids {
        let param = store.load_param(&id).expect("load shipped model");
        let scorer = RankScorer::new(&param);
        let pairs: &[(f64, f64)] = &[
            (PROTON, H2O + PROTON),
            (PROTON + NH3, H2O + PROTON + Z_DOT_OFFSET),
        ];
        let mut scored_partitions = 0usize;
        for partition in &param.partitions {
            let logs = scorer.partition_ion_logs(partition);
            for &(prefix_offset, suffix_offset) in pairs {
                let has_prefix = logs.iter().any(|(ion, _)| is_ion(ion, true, prefix_offset));
                let has_suffix = logs
                    .iter()
                    .any(|(ion, _)| is_ion(ion, false, suffix_offset));
                if has_prefix && has_suffix {
                    scored_partitions += 1;
                }
            }
        }
        let mut seen: Vec<String> = Vec::new();
        for partition in &param.partitions {
            for (ion, _) in scorer.partition_ion_logs(partition) {
                let text = format!("{ion:?}={:?}", offset_of(ion));
                if !seen.contains(&text) {
                    seen.push(text);
                }
            }
        }
        assert!(
            scored_partitions > 0,
            "{id}: no partition scores both primary ion types from the stored keys; \
             the model holds {seen:?}"
        );
    }
}
