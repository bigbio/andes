//! Trainer quality gate (replaces the former cross-language parity gate):
//! determinism, blob round-trip at the real 18-feature width, and signal/noise
//! separation. The GBDT is trained in Rust and emitted directly as the SoA
//! `GbdtPeakModel` the scoring walker consumes — this locks the trainer→walker
//! contract.
use model_train::gbdt::train::{train_gbdt, Dataset, TrainParams};
use scoring_crate::gbdt_eval::GbdtPeakModel;

/// Deterministic separable dataset at the production feature width (18).
/// feature 0 drives the label; the rest are seeded noise. No `rand` crate.
fn separable_18(n: usize, seed: u64) -> Dataset {
    const NF: usize = 18;
    let mut s = seed;
    let next = |s: &mut u64| -> f32 {
        *s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*s >> 33) as f32) / ((1u64 << 31) as f32)
    };
    let mut x = vec![0.0f32; n * NF];
    let mut y = vec![0u8; n];
    let mut groups = vec![0u32; n];
    for r in 0..n {
        for f in 0..NF {
            x[r * NF + f] = next(&mut s);
        }
        y[r] = u8::from(x[r * NF] + 0.2 * next(&mut s) > 0.6);
        groups[r] = (r as u32) % 60;
    }
    Dataset {
        x,
        y,
        groups,
        n_features: NF,
    }
}

#[test]
fn trainer_is_deterministic() {
    let ds = separable_18(3000, 11);
    let m1 = train_gbdt(&ds, &TrainParams::default(), 42).expect("separable data passes the gate");
    let m2 = train_gbdt(&ds, &TrainParams::default(), 42).expect("separable data passes the gate");
    assert_eq!(
        m1.to_bytes(),
        m2.to_bytes(),
        "same seed must yield identical model bytes"
    );
}

#[test]
fn trained_model_roundtrips_through_blob_and_separates() {
    let ds = separable_18(3000, 5);
    let model =
        train_gbdt(&ds, &TrainParams::default(), 7).expect("separable data passes the gate");
    assert!(
        !model.trees.is_empty(),
        "should fit at least one tree on separable data"
    );
    assert_eq!(model.n_features, 18);
    assert!(model.apply_sigmoid);

    // Round-trip through the AGBD blob (the exact bytes the store holds and the
    // walker decodes) — must decode and be structurally identical.
    let back = GbdtPeakModel::from_bytes(&model.to_bytes()).expect("blob must decode");
    assert_eq!(
        back, model,
        "model must survive blob round-trip byte-for-byte"
    );

    // Separation holds through the decoded model.
    let mut hi = vec![0.0f32; 18];
    hi[0] = 0.95;
    let mut lo = vec![0.0f32; 18];
    lo[0] = 0.05;
    assert!(
        back.predict_logit(&hi) > back.predict_logit(&lo) + 1.0,
        "decoded model must separate signal from noise: hi={} lo={}",
        back.predict_logit(&hi),
        back.predict_logit(&lo)
    );
}

#[test]
fn different_seeds_can_differ_but_both_separate() {
    // Sanity: two different seeds both yield separating models (not a determinism
    // contradiction — just that training is robust, not seed-fragile).
    let ds = separable_18(3000, 5);
    for seed in [1u64, 2, 999] {
        let m =
            train_gbdt(&ds, &TrainParams::default(), seed).expect("separable data passes the gate");
        let mut hi = vec![0.0f32; 18];
        hi[0] = 0.95;
        let mut lo = vec![0.0f32; 18];
        lo[0] = 0.05;
        assert!(
            m.predict_logit(&hi) > m.predict_logit(&lo),
            "seed {seed} model failed to separate"
        );
    }
}
