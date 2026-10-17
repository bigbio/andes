//! FDR within each evidence class, after rescoring.
//!
//! One q-value over every PSM lets a small class (isoform-, variant- or
//! chain-specific peptides) pass at a much higher error rate than the large
//! canonical class: its decoys are diluted among canonical targets. Each PSM is
//! assigned a class from its proteins' entry kinds and gets a q-value within
//! its class (see [`apply_class_fdr`]); a target's q-value becomes the larger
//! of the global and the class q-value. A PSM then passes a cutoff only if it
//! passes it both overall and within its class.

use std::collections::HashMap;

use model::EntryKind;
use search::SearchIndex;

use crate::percolator::PercolatorPsm;

/// Targets and decoys of one class, as counted for the run log.
#[derive(Debug, Clone, PartialEq)]
pub struct ClassCount {
    pub class: String,
    pub targets: usize,
    pub decoys: usize,
    /// How the class q-values were estimated.
    pub method: ClassMethod,
    /// Storey's estimate of the class's null fraction (1.0 under `ClassTdc`).
    pub pi0: f64,
}

/// How a class's q-values are estimated.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClassMethod {
    /// p-values against every decoy of the run, then Storey q-values within
    /// the class: the class's own decoys score no higher than the pool, so
    /// the pooled null is conservative for it.
    PooledNull,
    /// Target-decoy competition inside the class, q = (decoys + 1) / targets:
    /// the class's decoys score higher than the pool, or are too few to check.
    ClassTdc,
}

impl ClassMethod {
    pub fn name(self) -> &'static str {
        match self {
            ClassMethod::PooledNull => "pooled-null",
            ClassMethod::ClassTdc => "class-TDC",
        }
    }
}

/// The evidence class of a PSM from its `;`-separated protein accessions:
/// `canonical` when any protein is a Reference entry (canonical or FASTA),
/// otherwise the kinds present (`isoform`, `variant`, `chain`), joined by `+`.
/// An accession the index does not know counts as Reference.
pub fn class_of(proteins: &str, kinds: &HashMap<&str, EntryKind>) -> String {
    let (mut isoform, mut variant, mut chain) = (false, false, false);
    for acc in proteins.split(';').filter(|a| !a.is_empty()) {
        match kinds.get(acc).copied().unwrap_or(EntryKind::Reference) {
            EntryKind::Reference => return "canonical".to_string(),
            EntryKind::Isoform => isoform = true,
            EntryKind::Variant => variant = true,
            EntryKind::Chain => chain = true,
        }
    }
    let parts: Vec<&str> = [(chain, "chain"), (isoform, "isoform"), (variant, "variant")]
        .into_iter()
        .filter_map(|(on, name)| on.then_some(name))
        .collect();
    if parts.is_empty() {
        "canonical".to_string()
    } else {
        parts.join("+")
    }
}

/// A class's targets (rank, PSM id) and decoy ranks.
type ClassRows = (Vec<(f64, String)>, Vec<f64>);

/// Fewer decoys than this in a class cannot show its null is no heavier than
/// the pool; the class then uses its own target-decoy competition.
const MIN_GUARD_DECOYS: usize = 20;

/// Raise each target's q-value to its class q-value. Does nothing when every
/// PSM is in one class. Returns the per-class counts, largest first.
///
/// Each target gets an empirical p-value against every decoy of the run (the
/// pooled null), and Storey q-values are computed within its class. This
/// borrows the power of the large canonical decoy set, and it is conservative
/// for a class whose own decoys score no higher than the pool, which is
/// checked: a class whose decoys fall in the pool's top 5% more often than 5%
/// (one-sided z > 2.33), or that has fewer than [`MIN_GUARD_DECOYS`] decoys,
/// uses target-decoy competition inside the class instead.
pub fn apply_class_fdr(
    targets: &mut HashMap<String, PercolatorPsm>,
    decoys: &HashMap<String, PercolatorPsm>,
    idx: &SearchIndex,
) -> Vec<ClassCount> {
    let kinds: HashMap<&str, EntryKind> = idx
        .db
        .proteins
        .iter()
        .map(|p| (p.accession.as_str(), p.entry_kind()))
        .collect();
    // A row's rank: the rescorer's score, or (without one) its q-value inverted.
    let rank = |p: &PercolatorPsm| p.score.unwrap_or(-p.q_value);
    let mut by_class: HashMap<String, ClassRows> = HashMap::new();
    for (id, p) in targets.iter() {
        by_class
            .entry(class_of(&p.proteins, &kinds))
            .or_default()
            .0
            .push((rank(p), id.clone()));
    }
    for p in decoys.values() {
        by_class
            .entry(class_of(&p.proteins, &kinds))
            .or_default()
            .1
            .push(rank(p));
    }
    let single = by_class.len() < 2;
    let mut pooled: Vec<f64> = decoys.values().map(rank).collect();
    pooled.sort_by(f64::total_cmp);
    let n = pooled.len();
    // Share of pooled decoys scoring at least `s`, with the usual +1.
    let p_of = |s: f64| -> f64 {
        let ge = n - pooled.partition_point(|&x| x < s);
        (ge + 1) as f64 / (n + 1) as f64
    };
    let mut counts = Vec::with_capacity(by_class.len());
    for (class, (ts, ds)) in by_class.iter_mut() {
        let method = if ds.len() < MIN_GUARD_DECOYS || heavier_than_pool(ds, &p_of) {
            ClassMethod::ClassTdc
        } else {
            ClassMethod::PooledNull
        };
        let mut pi0 = 1.0;
        if !single {
            let q: Vec<f64> = match method {
                ClassMethod::PooledNull => {
                    let p: Vec<f64> = ts.iter().map(|t| p_of(t.0)).collect();
                    let (q, est) = storey_q(&p);
                    pi0 = est;
                    q
                }
                ClassMethod::ClassTdc => class_tdc_q(ts, ds),
            };
            for (t, qv) in ts.iter().zip(q) {
                if let Some(p) = targets.get_mut(&t.1) {
                    p.q_value = p.q_value.max(qv);
                }
            }
        }
        counts.push(ClassCount {
            class: class.clone(),
            targets: ts.len(),
            decoys: ds.len(),
            method,
            pi0,
        });
    }
    counts.sort_by(|a, b| b.targets.cmp(&a.targets).then(a.class.cmp(&b.class)));
    counts
}

/// Whether a class's decoys land in the pool's top 5% significantly more
/// often than 5% of the time (one-sided z > 2.33).
fn heavier_than_pool(ds: &[f64], p_of: &impl Fn(f64) -> f64) -> bool {
    let m = ds.len() as f64;
    let k = ds.iter().filter(|&&s| p_of(s) <= 0.05).count() as f64;
    let expected = 0.05 * m;
    (k - expected) / (m * 0.05 * 0.95).sqrt() > 2.33
}

/// Storey q-values (lambda = 0.5) of `p`, in the same order, and the null
/// fraction estimate.
fn storey_q(p: &[f64]) -> (Vec<f64>, f64) {
    let m = p.len();
    if m == 0 {
        return (Vec::new(), 1.0);
    }
    let pi0 = (p.iter().filter(|&&x| x > 0.5).count() as f64 / (0.5 * m as f64)).min(1.0);
    let mut order: Vec<usize> = (0..m).collect();
    order.sort_by(|&a, &b| p[a].total_cmp(&p[b]));
    let mut q = vec![1.0; m];
    let mut running = 1.0f64;
    for (r, &i) in order.iter().enumerate().rev() {
        running = running.min(pi0 * p[i] * m as f64 / (r + 1) as f64);
        q[i] = running.min(1.0);
    }
    (q, pi0)
}

/// Target-decoy q-values inside a class, q = (decoys + 1) / targets made
/// monotone, in the order of `ts`. At equal rank a decoy counts first.
fn class_tdc_q(ts: &[(f64, String)], ds: &[f64]) -> Vec<f64> {
    let mut rows: Vec<(f64, bool, usize)> = ts
        .iter()
        .enumerate()
        .map(|(i, t)| (t.0, true, i))
        .chain(ds.iter().map(|&s| (s, false, usize::MAX)))
        .collect();
    rows.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
    let mut q = vec![0.0f64; rows.len()];
    let (mut t, mut d) = (0usize, 0usize);
    for (k, r) in rows.iter().enumerate() {
        if r.1 {
            t += 1;
        } else {
            d += 1;
        }
        q[k] = if t == 0 {
            1.0
        } else {
            ((d + 1) as f64 / t as f64).min(1.0)
        };
    }
    for k in (0..q.len().saturating_sub(1)).rev() {
        q[k] = q[k].min(q[k + 1]);
    }
    let mut out = vec![1.0; ts.len()];
    for (r, qv) in rows.iter().zip(q) {
        if r.1 {
            out[r.2] = qv;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::protein::{EntryOrigin, Protein, ProteinAnnotation, ProteinDb};
    use std::sync::Arc;

    fn protein(acc: &str, kind: EntryKind) -> Protein {
        let annotation = (kind != EntryKind::Reference).then(|| {
            Arc::new(
                ProteinAnnotation::new(Vec::new(), None).with_origin(Some(EntryOrigin {
                    kind,
                    family: Arc::from("P1"),
                    feature: Arc::from("F"),
                    variant_delta_da: 0.0,
                    knockoff_of: None,
                })),
            )
        });
        Protein {
            accession: acc.into(),
            description: String::new(),
            sequence: b"PEPTIDEK".to_vec(),
            annotation,
        }
    }

    /// Canonical P1, its isoform P1-2, a variant P1_VAR_1, and their decoys.
    fn index() -> SearchIndex {
        let target = ProteinDb {
            proteins: vec![
                protein("sp|P1|A", EntryKind::Reference),
                protein("sp|P1-2|A", EntryKind::Isoform),
                protein("sp|P1_VAR_1|A", EntryKind::Variant),
            ],
        };
        SearchIndex::from_target_db(&target, "XXX")
    }

    fn psm(id: &str, score: f64, q: f64, proteins: &str) -> (String, PercolatorPsm) {
        (
            id.to_string(),
            PercolatorPsm {
                psm_id: id.to_string(),
                q_value: q,
                pep: 0.0,
                score: Some(score),
                peptide: "K.PEPTIDEK.A".into(),
                proteins: proteins.into(),
            },
        )
    }

    fn kinds(idx: &SearchIndex) -> HashMap<&str, EntryKind> {
        idx.db
            .proteins
            .iter()
            .map(|p| (p.accession.as_str(), p.entry_kind()))
            .collect()
    }

    #[test]
    fn class_follows_the_most_canonical_protein() {
        let idx = index();
        let k = kinds(&idx);
        assert_eq!(class_of("sp|P1|A;sp|P1-2|A", &k), "canonical");
        assert_eq!(class_of("sp|P1-2|A", &k), "isoform");
        assert_eq!(class_of("sp|P1_VAR_1|A;sp|P1-2|A", &k), "isoform+variant");
        assert_eq!(
            class_of("XXX_sp|P1-2|A", &k),
            "isoform",
            "decoys keep the class"
        );
        assert_eq!(class_of("unknown", &k), "canonical");
    }

    #[test]
    fn one_class_leaves_q_values_alone() {
        let idx = index();
        let mut t: HashMap<_, _> = (0..20)
            .map(|i| psm(&format!("t{i}"), 10.0 - i as f64, 0.001, "sp|P1|A"))
            .collect();
        let d: HashMap<_, _> = [psm("d0", 9.5, 0.5, "XXX_sp|P1|A")].into_iter().collect();
        let before = t.clone();
        let counts = apply_class_fdr(&mut t, &d, &idx);
        assert_eq!(t, before);
        assert_eq!(counts.len(), 1);
    }

    /// `n` scores spread evenly over `[lo, hi]`.
    fn spread(n: usize, lo: f64, hi: f64) -> Vec<f64> {
        (0..n)
            .map(|i| lo + (hi - lo) * i as f64 / n.max(2) as f64)
            .collect()
    }

    #[test]
    fn pooled_null_gives_a_small_class_power() {
        let idx = index();
        // 2000 canonical targets far above 1000 canonical decoys (0..5);
        // isoform: 1000 decoys scoring lower (0..2), 20 targets above every
        // decoy and 80 null-like ones. Globally all pass at q 0.001.
        let mut t: HashMap<_, _> = spread(2000, 10.0, 20.0)
            .into_iter()
            .enumerate()
            .map(|(i, s)| psm(&format!("c{i}"), s, 0.001, "sp|P1|A"))
            .chain(
                spread(20, 6.0, 7.0)
                    .into_iter()
                    .enumerate()
                    .map(|(i, s)| psm(&format!("i{i}"), s, 0.001, "sp|P1-2|A")),
            )
            .chain(
                spread(80, 0.0, 2.0)
                    .into_iter()
                    .enumerate()
                    .map(|(i, s)| psm(&format!("n{i}"), s, 0.5, "sp|P1-2|A")),
            )
            .collect();
        let d: HashMap<_, _> = spread(1000, 0.0, 5.0)
            .into_iter()
            .enumerate()
            .map(|(i, s)| psm(&format!("dc{i}"), s, 0.5, "XXX_sp|P1|A"))
            .chain(
                spread(1000, 0.0, 2.0)
                    .into_iter()
                    .enumerate()
                    .map(|(i, s)| psm(&format!("di{i}"), s, 0.5, "XXX_sp|P1-2|A")),
            )
            .collect();
        let counts = apply_class_fdr(&mut t, &d, &idx);
        let iso = counts.iter().find(|c| c.class == "isoform").unwrap();
        assert_eq!(iso.method, ClassMethod::PooledNull);
        // p = 1/2001 for each of the 20; Storey q = pi0 * p * 100 / rank, at
        // most 100 / (20 * 2001) for the 20th. Class target-decoy would give
        // (0 + 1) / 20 = 5%.
        for k in 0..20 {
            let q = t[&format!("i{k}")].q_value;
            assert!(q <= 0.01, "i{k}: {q}");
        }
        assert!(t["n0"].q_value > 0.01);
    }

    #[test]
    fn a_class_whose_decoys_score_higher_uses_its_own_competition() {
        let idx = index();
        // Isoform decoys (4..6) score above the canonical decoys (0..2): the
        // pooled null would be anti-conservative for the isoform class.
        let mut t: HashMap<_, _> = spread(2000, 10.0, 20.0)
            .into_iter()
            .enumerate()
            .map(|(i, s)| psm(&format!("c{i}"), s, 0.001, "sp|P1|A"))
            .chain(
                spread(20, 7.0, 8.0)
                    .into_iter()
                    .enumerate()
                    .map(|(i, s)| psm(&format!("i{i}"), s, 0.001, "sp|P1-2|A")),
            )
            .collect();
        let d: HashMap<_, _> = spread(2000, 0.0, 2.0)
            .into_iter()
            .enumerate()
            .map(|(i, s)| psm(&format!("dc{i}"), s, 0.5, "XXX_sp|P1|A"))
            .chain(
                spread(200, 4.0, 6.0)
                    .into_iter()
                    .enumerate()
                    .map(|(i, s)| psm(&format!("di{i}"), s, 0.5, "XXX_sp|P1-2|A")),
            )
            .collect();
        let counts = apply_class_fdr(&mut t, &d, &idx);
        let iso = counts.iter().find(|c| c.class == "isoform").unwrap();
        assert_eq!(iso.method, ClassMethod::ClassTdc);
        // No isoform decoy above them: q = (0 + 1) / 20 at best.
        let best = (0..20)
            .map(|k| t[&format!("i{k}")].q_value)
            .fold(f64::INFINITY, f64::min);
        assert!((best - 0.05).abs() < 1e-9, "{best}");
    }
}
