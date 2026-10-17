//! Knockoff q-values for variant PSMs, after rescoring and the class FDR.
//!
//! A knockoff entry carries a real variant's substitution at a position no
//! known variant has, so it is almost certainly absent from the sample. A
//! decoy only models random matches; a knockoff is also matched by a canonical
//! peptide carrying a modification of the substitution's mass or picked at the
//! wrong isotope, which is how most false variant calls arise. The knockoff
//! PSMs above a score estimate the false real variant PSMs above it, within
//! each substitution class (see [`apply_knockoff_fdr`]).

use std::collections::{HashMap, HashSet};

use model::EntryKind;
use search::SearchIndex;

use crate::percolator::PercolatorPsm;

/// Substitution classes, by the substitution's mass: the ones a common
/// modification mimics, then every other.
const STRATA: [(&str, f64); 3] = [
    ("deamidation-like", 0.984016),
    ("methyl-like", 14.01565),
    ("oxidation-like", 15.994915),
];
const OTHER: &str = "other";
const STRATUM_TOL_DA: f64 = 0.02;

fn stratum(delta_da: f64) -> &'static str {
    STRATA
        .iter()
        .find(|(_, m)| (delta_da.abs() - m).abs() < STRATUM_TOL_DA)
        .map_or(OTHER, |(name, _)| name)
}

/// One substitution class, as counted for the run log.
#[derive(Debug, Clone, PartialEq)]
pub struct KnockoffCount {
    pub stratum: &'static str,
    /// Real variant PSMs at q ≤ the cutoff before the knockoff q-values.
    pub real_before: usize,
    /// Of them, the ones whose variant has knockoffs.
    pub tested_before: usize,
    /// Knockoff PSMs at q ≤ the cutoff, each counting 1/k for a variant with
    /// k knockoffs.
    pub knockoff_weighted: f64,
    /// Real variant PSMs still at q ≤ the cutoff after.
    pub real_after: usize,
}

impl KnockoffCount {
    /// Estimated share of the tested real PSMs at the cutoff that are false.
    pub fn false_share(&self) -> f64 {
        if self.tested_before == 0 {
            f64::NAN
        } else {
            self.knockoff_weighted / self.tested_before as f64
        }
    }
}

/// What [`apply_knockoff_fdr`] did.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct KnockoffOutcome {
    /// Per substitution class, in a fixed order.
    pub counts: Vec<KnockoffCount>,
    /// Ids of the knockoff rows: never to be reported or quantified.
    pub knockoff_rows: HashSet<String>,
}

/// How a variant-class target row relates to the knockoffs.
enum Row {
    /// Every protein is a knockoff entry; `weight` is 1/k of the variant it
    /// copies.
    Knockoff { weight: f64 },
    /// A real variant row; `tested` when one of its variants has knockoffs.
    Real { tested: bool },
}

/// Raise each real variant target's q-value to its knockoff q-value, give
/// knockoff rows q = 1 and PEP = 1, and drop knockoff accessions from other
/// rows' protein lists. Does nothing when the database has no knockoff entry.
/// Returns per-class counts at `fdr` and the knockoff rows, which the caller
/// keeps out of every reported output.
///
/// Rows of the variant class only (no Reference, isoform or chain protein).
/// Within a substitution class, at score `t`:
/// FDP(t) = (1/K + Σ 1/k over knockoff rows ≥ t) / #(tested real rows ≥ t),
/// with K the most knockoffs of any variant; a real row's knockoff q-value is
/// the least FDP(t) over t at or below its score. Untested real rows (whose
/// variants have no knockoff) take the q-value at their score too: the error
/// rate is a property of the substitution class, which the tested rows sample.
pub fn apply_knockoff_fdr(
    targets: &mut HashMap<String, PercolatorPsm>,
    idx: &SearchIndex,
    fdr: f64,
) -> KnockoffOutcome {
    // Target knockoffs per real variant, keyed by (family, feature).
    let mut per_variant: HashMap<(&str, &str), usize> = HashMap::new();
    for p in idx
        .db
        .proteins
        .iter()
        .filter(|p| !idx.is_decoy_accession(&p.accession))
    {
        if let Some(o) = p.origin() {
            if let Some(of) = o.knockoff_of.as_deref() {
                *per_variant.entry((&o.family, of)).or_insert(0) += 1;
            }
        }
    }
    if per_variant.is_empty() {
        return KnockoffOutcome::default();
    }
    let k_max = per_variant.values().copied().max().unwrap_or(1) as f64;
    let proteins: HashMap<&str, &model::Protein> = idx
        .db
        .proteins
        .iter()
        .map(|p| (p.accession.as_str(), p))
        .collect();

    let rank = |p: &PercolatorPsm| p.score.unwrap_or(-p.q_value);
    // Per stratum: (rank, row id, row kind).
    let mut by_stratum: HashMap<&'static str, Vec<(f64, String, Row)>> = HashMap::new();
    for (id, psm) in targets.iter() {
        let mut real_delta: Option<f64> = None;
        let mut knockoff: Option<(f64, f64)> = None;
        let mut tested = false;
        let mut variant_only = true;
        let mut all_knockoff = true;
        for acc in psm.proteins.split(';').filter(|a| !a.is_empty()) {
            let Some(p) = proteins.get(acc) else {
                variant_only = false;
                break;
            };
            if p.entry_kind() != EntryKind::Variant {
                variant_only = false;
                break;
            }
            let o = p.origin().expect("a variant entry has an origin");
            match o.knockoff_of.as_deref() {
                Some(of) => {
                    let k = per_variant[&(&*o.family, of)] as f64;
                    let w = 1.0 / k;
                    knockoff = Some(match knockoff {
                        Some((d, prev)) => (d, prev.min(w)),
                        None => (o.variant_delta_da, w),
                    });
                }
                None => {
                    all_knockoff = false;
                    real_delta.get_or_insert(o.variant_delta_da);
                    tested |= per_variant.contains_key(&(&*o.family, &*o.feature));
                }
            }
        }
        if !variant_only {
            continue;
        }
        let (delta, row) = match (all_knockoff, knockoff, real_delta) {
            (true, Some((d, w)), _) => (d, Row::Knockoff { weight: w }),
            (false, _, Some(d)) => (d, Row::Real { tested }),
            _ => continue,
        };
        by_stratum
            .entry(stratum(delta))
            .or_default()
            .push((rank(psm), id.clone(), row));
    }

    let mut counts = Vec::new();
    let mut knockoff_rows = HashSet::new();
    for (name, mut rows) in by_stratum {
        let mut count = KnockoffCount {
            stratum: name,
            real_before: 0,
            tested_before: 0,
            knockoff_weighted: 0.0,
            real_after: 0,
        };
        for (_, id, row) in &rows {
            if targets[id].q_value <= fdr {
                match row {
                    Row::Knockoff { weight } => count.knockoff_weighted += weight,
                    Row::Real { tested } => {
                        count.real_before += 1;
                        count.tested_before += *tested as usize;
                    }
                }
            }
        }
        // Best first; at equal rank a knockoff counts first.
        rows.sort_by(|a, b| {
            b.0.total_cmp(&a.0)
                .then(matches!(b.2, Row::Knockoff { .. }).cmp(&matches!(a.2, Row::Knockoff { .. })))
        });
        let mut fdp = Vec::with_capacity(rows.len());
        let (mut w, mut t) = (0.0f64, 0usize);
        for (_, _, row) in &rows {
            match row {
                Row::Knockoff { weight } => w += weight,
                Row::Real { tested } => t += *tested as usize,
            }
            fdp.push(((1.0 / k_max + w) / t.max(1) as f64).min(1.0));
        }
        // q at a row: the least FDP over thresholds at or below its rank.
        for i in (0..fdp.len().saturating_sub(1)).rev() {
            fdp[i] = fdp[i].min(fdp[i + 1]);
        }
        for ((_, id, row), q) in rows.iter().zip(fdp) {
            let p = targets.get_mut(id).expect("row of the map");
            match row {
                Row::Knockoff { .. } => {
                    p.q_value = 1.0;
                    p.pep = 1.0;
                    knockoff_rows.insert(id.clone());
                }
                Row::Real { .. } => {
                    p.q_value = p.q_value.max(q);
                    count.real_after += (p.q_value <= fdr) as usize;
                }
            }
        }
        counts.push(count);
    }

    // A knockoff accession on any other row (a knockoff peptide another entry
    // also holds) is not reported either.
    for p in targets.values_mut() {
        if p.proteins
            .split(';')
            .any(|a| proteins.get(a).is_some_and(|e| e.is_knockoff()))
        {
            let kept: Vec<&str> = p
                .proteins
                .split(';')
                .filter(|a| !proteins.get(a).is_some_and(|e| e.is_knockoff()))
                .collect();
            if !kept.is_empty() {
                p.proteins = kept.join(";");
            }
        }
    }
    counts.sort_by_key(|c| {
        STRATA
            .iter()
            .position(|(n, _)| *n == c.stratum)
            .unwrap_or(STRATA.len())
    });
    KnockoffOutcome {
        counts,
        knockoff_rows,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::protein::{EntryOrigin, Protein, ProteinAnnotation, ProteinDb};
    use std::sync::Arc;

    fn protein(acc: &str, kind: EntryKind, feature: &str, delta: f64, of: Option<&str>) -> Protein {
        let annotation = (kind != EntryKind::Reference).then(|| {
            Arc::new(
                ProteinAnnotation::new(Vec::new(), None).with_origin(Some(EntryOrigin {
                    kind,
                    family: Arc::from("P1"),
                    feature: Arc::from(feature),
                    variant_delta_da: delta,
                    knockoff_of: of.map(Arc::from),
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

    /// Canonical P1; variant V1 (N→D, deamidation-like) with two knockoffs;
    /// variant V2 (other, +30 Da) with none.
    fn index() -> SearchIndex {
        let target = ProteinDb {
            proteins: vec![
                protein("sp|P1|A", EntryKind::Reference, "", 0.0, None),
                protein("sp|P1_V1|A", EntryKind::Variant, "V1", 0.984016, None),
                protein(
                    "sp|P1_KO1_V1|A",
                    EntryKind::Variant,
                    "KO1_V1",
                    0.984016,
                    Some("V1"),
                ),
                protein(
                    "sp|P1_KO2_V1|A",
                    EntryKind::Variant,
                    "KO2_V1",
                    0.984016,
                    Some("V1"),
                ),
                protein("sp|P1_V2|A", EntryKind::Variant, "V2", 30.0, None),
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

    #[test]
    fn no_knockoff_entries_leave_the_map_alone() {
        let target = ProteinDb {
            proteins: vec![protein("sp|P1_V2|A", EntryKind::Variant, "V2", 30.0, None)],
        };
        let idx = SearchIndex::from_target_db(&target, "XXX");
        let mut t: HashMap<_, _> = [psm("a", 5.0, 0.001, "sp|P1_V2|A")].into_iter().collect();
        let before = t.clone();
        assert_eq!(apply_knockoff_fdr(&mut t, &idx, 0.01), KnockoffOutcome::default());
        assert_eq!(t, before);
    }

    #[test]
    fn knockoffs_above_a_real_row_raise_its_q_and_are_never_accepted() {
        let idx = index();
        // Deamidation-like: 4 real rows (8, 7, 6, 5) and two knockoff rows
        // (7.5, 5.5), each weighing 1/2. K = 2.
        let mut t: HashMap<_, _> = [
            psm("r8", 8.0, 0.001, "sp|P1_V1|A"),
            psm("r7", 7.0, 0.001, "sp|P1_V1|A"),
            psm("r6", 6.0, 0.001, "sp|P1_V1|A"),
            psm("r5", 5.0, 0.001, "sp|P1_V1|A"),
            psm("k75", 7.5, 0.001, "sp|P1_KO1_V1|A"),
            psm("k55", 5.5, 0.001, "sp|P1_KO2_V1|A"),
            psm("c", 9.0, 0.001, "sp|P1|A;sp|P1_KO1_V1|A"),
        ]
        .into_iter()
        .collect();
        let out = apply_knockoff_fdr(&mut t, &idx, 0.01);
        let counts = out.counts;
        assert_eq!(
            out.knockoff_rows,
            HashSet::from(["k75".to_string(), "k55".to_string()])
        );
        // FDP at r8: (0.5 + 0) / 1 = 0.5; at r7: (0.5 + 0.5) / 2 = 0.5;
        // r6: 1.0 / 3; r5: (0.5 + 1.0) / 4 = 0.375. Running minimum from
        // below: r5 0.375, r6 1/3, r7 1/3, r8 1/3.
        assert!((t["r5"].q_value - 0.375).abs() < 1e-9);
        for id in ["r6", "r7", "r8"] {
            assert!((t[id].q_value - 1.0 / 3.0).abs() < 1e-9, "{id}");
        }
        assert_eq!((t["k75"].q_value, t["k75"].pep), (1.0, 1.0));
        assert_eq!(t["c"].q_value, 0.001, "canonical rows keep their q-value");
        assert_eq!(
            t["c"].proteins, "sp|P1|A",
            "knockoff accessions are not reported"
        );
        let c = &counts[0];
        assert_eq!(c.stratum, "deamidation-like");
        assert_eq!((c.real_before, c.tested_before, c.real_after), (4, 4, 0));
        assert!((c.knockoff_weighted - 1.0).abs() < 1e-9);
        assert!((c.false_share() - 0.25).abs() < 1e-9);
    }

    #[test]
    fn real_rows_far_above_the_knockoffs_keep_a_small_q() {
        let idx = index();
        let mut t: HashMap<_, _> = (0..100)
            .map(|i| psm(&format!("r{i}"), 100.0 - i as f64, 0.001, "sp|P1_V1|A"))
            .chain([psm("k", -5.0, 0.5, "sp|P1_KO1_V1|A")])
            .collect();
        apply_knockoff_fdr(&mut t, &idx, 0.01);
        // (1/2) / 100 at the bottom real row, the floor of the estimate.
        assert!((t["r99"].q_value - 0.005).abs() < 1e-9);
        assert!((t["r0"].q_value - 0.005).abs() < 1e-9);
    }

    #[test]
    fn untested_variants_take_their_class_estimate() {
        let idx = index();
        // V2 (other) has no knockoff: no tested row, so FDP = (1/K) / 1.
        let mut t: HashMap<_, _> = [psm("v2", 5.0, 0.001, "sp|P1_V2|A")].into_iter().collect();
        let counts = apply_knockoff_fdr(&mut t, &idx, 0.01).counts;
        assert!((t["v2"].q_value - 0.5).abs() < 1e-9);
        assert_eq!(counts[0].stratum, "other");
        assert_eq!(counts[0].tested_before, 0);
    }
}
