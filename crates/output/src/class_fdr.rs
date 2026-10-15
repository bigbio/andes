//! FDR within each evidence class, after rescoring.
//!
//! One q-value over every PSM lets a small class (isoform-, variant- or
//! chain-specific peptides) pass at a much higher error rate than the large
//! canonical class: its decoys are diluted among canonical targets. Each PSM is
//! assigned a class from its proteins' entry kinds, target-decoy competition is
//! run inside each class on the rescorer's score, and a target's q-value
//! becomes the larger of the global and the class q-value. A PSM then passes a
//! cutoff only if it passes it both overall and within its class.

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

/// Raise each target's q-value to its class q-value. Does nothing when every
/// PSM is in one class. Returns the per-class counts, largest first.
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
    let mut by_class: HashMap<String, Vec<(f64, bool, String)>> = HashMap::new();
    for (id, p) in targets.iter() {
        by_class
            .entry(class_of(&p.proteins, &kinds))
            .or_default()
            .push((rank(p), true, id.clone()));
    }
    for p in decoys.values() {
        by_class
            .entry(class_of(&p.proteins, &kinds))
            .or_default()
            .push((rank(p), false, String::new()));
    }
    let mut counts: Vec<ClassCount> = by_class
        .iter()
        .map(|(c, rows)| ClassCount {
            class: c.clone(),
            targets: rows.iter().filter(|r| r.1).count(),
            decoys: rows.iter().filter(|r| !r.1).count(),
        })
        .collect();
    counts.sort_by(|a, b| b.targets.cmp(&a.targets).then(a.class.cmp(&b.class)));
    if by_class.len() < 2 {
        return counts;
    }
    for rows in by_class.values_mut() {
        // Best first; at equal rank a decoy counts before a target.
        rows.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
        let mut q = vec![0.0f64; rows.len()];
        let (mut t, mut d) = (0usize, 0usize);
        for (i, r) in rows.iter().enumerate() {
            if r.1 {
                t += 1;
            } else {
                d += 1;
            }
            q[i] = if t == 0 {
                1.0
            } else {
                ((d + 1) as f64 / t as f64).min(1.0)
            };
        }
        // Monotone: a row's q-value is the lowest at or below it.
        for i in (0..q.len().saturating_sub(1)).rev() {
            q[i] = q[i].min(q[i + 1]);
        }
        for (r, qv) in rows.iter().zip(q) {
            if r.1 {
                if let Some(p) = targets.get_mut(&r.2) {
                    p.q_value = p.q_value.max(qv);
                }
            }
        }
    }
    counts
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

    #[test]
    fn a_small_class_gets_its_own_error_rate() {
        let idx = index();
        // 200 canonical targets above every decoy; 10 isoform targets with two
        // isoform decoys among them. Globally all pass at q 0.001.
        let mut t: HashMap<_, _> = (0..200)
            .map(|i| psm(&format!("c{i}"), 100.0 - i as f64 * 0.1, 0.001, "sp|P1|A"))
            .chain((0..10).map(|i| psm(&format!("i{i}"), 50.0 - i as f64, 0.001, "sp|P1-2|A")))
            .collect();
        let d: HashMap<_, _> = [
            psm("d1", 47.5, 0.5, "XXX_sp|P1-2|A"),
            psm("d2", 44.5, 0.5, "XXX_sp|P1-2|A"),
            psm("d3", 1.0, 0.5, "XXX_sp|P1|A"),
        ]
        .into_iter()
        .collect();
        let counts = apply_class_fdr(&mut t, &d, &idx);
        // Isoform class: q = (decoys at or above + 1) / targets at or above,
        // made monotone. The lowest value, (2 + 1) / 10 at i9, is at or below
        // every isoform target, so all of them get it.
        assert!(
            (t["i0"].q_value - 0.3).abs() < 1e-9,
            "{}",
            t["i0"].q_value
        );
        assert!(
            (t["i9"].q_value - 3.0 / 10.0).abs() < 1e-9,
            "{}",
            t["i9"].q_value
        );
        // Canonical class: no canonical decoy above them, (0 + 1) / n, but
        // never below the global q-value.
        assert!(
            (t["c199"].q_value - 1.0 / 200.0).abs() < 1e-9,
            "{}",
            t["c199"].q_value
        );
        assert!(t["c0"].q_value >= 0.001);
        assert_eq!(
            counts,
            vec![
                ClassCount {
                    class: "canonical".into(),
                    targets: 200,
                    decoys: 1
                },
                ClassCount {
                    class: "isoform".into(),
                    targets: 10,
                    decoys: 2
                },
            ]
        );
    }
}
