use std::collections::{HashMap, HashSet};

use model::{EntryKind, PriorTier};
use search::SearchIndex;

use crate::percolator::PercolatorPsm;

const STRATA: [(&str, f64); 3] = [
    ("deamidation-like", 0.984016),
    ("methyl-like", 14.01565),
    ("oxidation-like", 15.994915),
];
const OTHER: &str = "other";
const STRATUM_TOL_DA: f64 = 0.02;
const MIN_STRATUM_TESTED: usize = 50;

fn stratum(delta_da: f64) -> &'static str {
    STRATA
        .iter()
        .find(|(_, m)| (delta_da.abs() - m).abs() < STRATUM_TOL_DA)
        .map_or(OTHER, |(name, _)| name)
}

#[derive(Debug, Clone, PartialEq)]
pub struct KnockoffCount {
    pub stratum: String,
    pub real_before: usize,
    pub tested_before: usize,
    pub knockoff_weighted: f64,
    pub real_after: usize,
}

impl KnockoffCount {
    pub fn false_share(&self) -> f64 {
        if self.tested_before == 0 {
            f64::NAN
        } else {
            self.knockoff_weighted / self.tested_before as f64
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct KnockoffOutcome {
    pub counts: Vec<KnockoffCount>,
    pub knockoff_rows: HashSet<String>,
}

type RankedRow = (f64, String, Row);

enum Row {
    Knockoff { weight: f64, k: usize },
    Real { tested: bool, k: usize },
}

impl Row {
    fn k(&self) -> usize {
        match *self {
            Row::Knockoff { k, .. } | Row::Real { k, .. } => k,
        }
    }
}

pub(crate) fn finite_sample_term(ks: impl Iterator<Item = usize>) -> f64 {
    let k_min = ks.filter(|&k| k != usize::MAX).min().unwrap_or(1).max(1);
    1.0 / k_min as f64
}

pub fn apply_knockoff_fdr(
    targets: &mut HashMap<String, PercolatorPsm>,
    idx: &SearchIndex,
    fdr: f64,
    knockoffs_requested: bool,
) -> KnockoffOutcome {
    let mut features: HashMap<(&str, &str), HashSet<&str>> = HashMap::new();
    for p in idx
        .db
        .proteins
        .iter()
        .filter(|p| !idx.is_decoy_accession(&p.accession))
    {
        if let Some(o) = p.origin() {
            if let Some(of) = o.knockoff_of.as_deref() {
                features
                    .entry((&o.family, of))
                    .or_default()
                    .insert(&o.feature);
            }
        }
    }
    let per_variant: HashMap<(&str, &str), usize> =
        features.into_iter().map(|(k, f)| (k, f.len())).collect();
    if per_variant.is_empty() && !knockoffs_requested {
        return KnockoffOutcome::default();
    }
    let proteins: HashMap<&str, &model::Protein> = idx
        .db
        .proteins
        .iter()
        .map(|p| (p.accession.as_str(), p))
        .collect();

    let rank = |p: &PercolatorPsm| p.score.unwrap_or(-p.q_value);
    let mut by_class_tier: HashMap<(&'static str, PriorTier), Vec<RankedRow>> = HashMap::new();
    for (id, psm) in targets.iter() {
        let mut real_delta: Option<f64> = None;
        let mut real_tier: Option<PriorTier> = None;
        let mut knockoff: Option<(f64, PriorTier)> = None;
        let mut copied: HashSet<(&str, &str)> = HashSet::new();
        let mut tested = false;
        let mut k_real = usize::MAX;
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
                    copied.insert((&*o.family, of));
                    knockoff.get_or_insert((o.variant_delta_da, o.tier));
                }
                None => {
                    all_knockoff = false;
                    real_delta.get_or_insert(o.variant_delta_da);
                    real_tier = Some(real_tier.map_or(o.tier, |t| t.min(o.tier)));
                    if let Some(&k) = per_variant.get(&(&*o.family, &*o.feature)) {
                        tested = true;
                        k_real = k_real.min(k);
                    }
                }
            }
        }
        if !variant_only {
            continue;
        }
        let (delta, tier, row) = match (all_knockoff, knockoff, real_delta) {
            (true, Some((d, t)), _) => {
                let w = copied.iter().map(|v| 1.0 / per_variant[v] as f64).sum();
                let k = copied.iter().map(|v| per_variant[v]).min().unwrap_or(1);
                (d, t, Row::Knockoff { weight: w, k })
            }
            (false, _, Some(d)) => (
                d,
                real_tier.unwrap_or_default(),
                Row::Real { tested, k: k_real },
            ),
            _ => continue,
        };
        by_class_tier
            .entry((stratum(delta), tier))
            .or_default()
            .push((rank(psm), id.clone(), row));
    }

    let tested_at_cutoff = |rows: &[RankedRow]| {
        rows.iter()
            .filter(|(_, id, row)| {
                matches!(row, Row::Real { tested: true, .. }) && targets[id].q_value <= fdr
            })
            .count()
    };
    let tiers: std::collections::BTreeSet<PriorTier> = by_class_tier.keys().map(|k| k.1).collect();
    let mut by_stratum: HashMap<String, Vec<RankedRow>> = HashMap::new();
    for ((class, tier), rows) in by_class_tier {
        let name = if tiers.len() < 2 {
            class.to_string()
        } else if tested_at_cutoff(&rows) >= MIN_STRATUM_TESTED {
            format!("{class}/{}", tier.name())
        } else {
            format!("{class}/pooled")
        };
        by_stratum.entry(name).or_default().extend(rows);
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
                    Row::Knockoff { weight, .. } => count.knockoff_weighted += weight,
                    Row::Real { tested, .. } => {
                        count.real_before += 1;
                        count.tested_before += *tested as usize;
                    }
                }
            }
        }
        rows.sort_by(|a, b| {
            b.0.total_cmp(&a.0)
                .then(matches!(b.2, Row::Knockoff { .. }).cmp(&matches!(a.2, Row::Knockoff { .. })))
        });
        let offset = finite_sample_term(rows.iter().map(|(_, _, row)| row.k()));
        let mut fdp = Vec::with_capacity(rows.len());
        let (mut w, mut t) = (0.0f64, 0usize);
        for (_, _, row) in &rows {
            match row {
                Row::Knockoff { weight, .. } => w += weight,
                Row::Real { tested, .. } => t += *tested as usize,
            }
            fdp.push(((offset + w) / t.max(1) as f64).min(1.0));
        }
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
                Row::Real { tested: true, .. } => {
                    p.q_value = p.q_value.max(q);
                    count.real_after += (p.q_value <= fdr) as usize;
                }
                Row::Real { tested: false, .. } => p.q_value = 1.0,
            }
        }
        counts.push(count);
    }

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
        let class = c.stratum.split('/').next().unwrap_or_default().to_string();
        let order = STRATA
            .iter()
            .position(|(n, _)| *n == class)
            .unwrap_or(STRATA.len());
        (order, c.stratum.clone())
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
        tiered(acc, kind, feature, delta, of, PriorTier::Catalogue)
    }

    fn tiered(
        acc: &str,
        kind: EntryKind,
        feature: &str,
        delta: f64,
        of: Option<&str>,
        tier: PriorTier,
    ) -> Protein {
        let annotation = (kind != EntryKind::Reference).then(|| {
            Arc::new(
                ProteinAnnotation::new(Vec::new(), None).with_origin(Some(EntryOrigin {
                    kind,
                    family: Arc::from("P1"),
                    feature: Arc::from(feature),
                    variant_delta_da: delta,
                    knockoff_of: of.map(Arc::from),
                    source: None,
                    tier,
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
                pooled_q_value: None,
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
        assert_eq!(
            apply_knockoff_fdr(&mut t, &idx, 0.01, false),
            KnockoffOutcome::default()
        );
        assert_eq!(t, before);
    }

    #[test]
    fn requested_knockoffs_with_none_built_reject_every_variant() {
        let target = ProteinDb {
            proteins: vec![protein("sp|P1_V2|A", EntryKind::Variant, "V2", 30.0, None)],
        };
        let idx = SearchIndex::from_target_db(&target, "XXX");
        let mut t: HashMap<_, _> = [psm("a", 5.0, 0.001, "sp|P1_V2|A")].into_iter().collect();
        let out = apply_knockoff_fdr(&mut t, &idx, 0.01, true);
        assert_eq!(t["a"].q_value, 1.0);
        assert_eq!(
            (out.counts[0].real_before, out.counts[0].tested_before),
            (1, 0)
        );
    }

    #[test]
    fn knockoff_counts_are_distinct_features_and_shared_rows_sum() {
        let mut proteins = index().db.proteins.clone();
        proteins.retain(|p| !p.accession.starts_with("XXX"));
        proteins.push(protein(
            "sp|P1_KO1_V1|A_refine",
            EntryKind::Variant,
            "KO1_V1",
            0.984016,
            Some("V1"),
        ));
        proteins.push(protein(
            "sp|P1_V3|A",
            EntryKind::Variant,
            "V3",
            0.984016,
            None,
        ));
        proteins.push(protein(
            "sp|P1_KO1_V3|A",
            EntryKind::Variant,
            "KO1_V3",
            0.984016,
            Some("V3"),
        ));
        let idx = SearchIndex::from_target_db(&ProteinDb { proteins }, "XXX");
        let mut t: HashMap<_, _> = [
            psm("r", 8.0, 0.001, "sp|P1_V1|A"),
            psm("k", 7.0, 0.001, "sp|P1_KO1_V1|A;sp|P1_KO1_V3|A"),
            psm("k2", 6.0, 0.001, "sp|P1_KO2_V1|A"),
        ]
        .into_iter()
        .collect();
        let counts = apply_knockoff_fdr(&mut t, &idx, 0.01, true).counts;
        assert_eq!(counts.len(), 1);
        assert!((counts[0].knockoff_weighted - (1.5 + 0.5)).abs() < 1e-9);
    }

    #[test]
    fn knockoffs_above_a_real_row_raise_its_q_and_are_never_accepted() {
        let idx = index();
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
        let out = apply_knockoff_fdr(&mut t, &idx, 0.01, true);
        let counts = out.counts;
        assert_eq!(
            out.knockoff_rows,
            HashSet::from(["k75".to_string(), "k55".to_string()])
        );
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
        apply_knockoff_fdr(&mut t, &idx, 0.01, true);
        assert!((t["r99"].q_value - 0.005).abs() < 1e-9);
        assert!((t["r0"].q_value - 0.005).abs() < 1e-9);
    }

    #[test]
    fn the_finite_sample_term_is_per_stratum() {
        let mut proteins = index().db.proteins.clone();
        proteins.retain(|p| !p.accession.starts_with("XXX"));
        proteins.push(protein(
            "sp|P1_V5|A",
            EntryKind::Variant,
            "V5",
            14.01565,
            None,
        ));
        for k in 1..=4 {
            let f = format!("KO{k}_V5");
            proteins.push(protein(
                &format!("sp|P1_{f}|A"),
                EntryKind::Variant,
                &f,
                14.01565,
                Some("V5"),
            ));
        }
        let idx = SearchIndex::from_target_db(&ProteinDb { proteins }, "XXX");
        let mut t: HashMap<_, _> = (0..100)
            .map(|i| psm(&format!("d{i}"), 100.0 - i as f64, 0.001, "sp|P1_V1|A"))
            .chain((0..100).map(|i| psm(&format!("m{i}"), 100.0 - i as f64, 0.001, "sp|P1_V5|A")))
            .collect();
        apply_knockoff_fdr(&mut t, &idx, 0.01, true);
        assert!(
            (t["d0"].q_value - 0.5 / 100.0).abs() < 1e-12,
            "{}",
            t["d0"].q_value
        );
        assert!(
            (t["m0"].q_value - 0.25 / 100.0).abs() < 1e-12,
            "{}",
            t["m0"].q_value
        );
        assert_eq!(finite_sample_term([4, 2, usize::MAX].into_iter()), 0.5);
        assert_eq!(finite_sample_term([usize::MAX].into_iter()), 1.0);
        assert_eq!(finite_sample_term(std::iter::empty()), 1.0);
    }

    #[test]
    fn untested_variants_are_not_accepted() {
        let idx = index();
        let mut t: HashMap<_, _> = [psm("v2", 5.0, 0.001, "sp|P1_V2|A")].into_iter().collect();
        let counts = apply_knockoff_fdr(&mut t, &idx, 0.01, true).counts;
        assert_eq!(t["v2"].q_value, 1.0);
        assert_eq!(counts[0].stratum, "other");
        assert_eq!(
            (
                counts[0].real_before,
                counts[0].tested_before,
                counts[0].real_after
            ),
            (1, 0, 0)
        );
        assert!(!apply_knockoff_fdr(&mut t, &idx, 0.01, true)
            .knockoff_rows
            .contains("v2"));
    }

    fn tiered_index() -> SearchIndex {
        let mut proteins = vec![protein("sp|P1|A", EntryKind::Reference, "", 0.0, None)];
        for (v, tier) in [
            ("S", PriorTier::Sample),
            ("C", PriorTier::Catalogue),
            ("M", PriorTier::Common),
        ] {
            proteins.push(tiered(
                &format!("sp|P1_{v}|A"),
                EntryKind::Variant,
                v,
                0.984016,
                None,
                tier,
            ));
            for k in 1..=2 {
                let f = format!("KO{k}_{v}");
                proteins.push(tiered(
                    &format!("sp|P1_{f}|A"),
                    EntryKind::Variant,
                    &f,
                    0.984016,
                    Some(v),
                    tier,
                ));
            }
        }
        SearchIndex::from_target_db(&ProteinDb { proteins }, "XXX")
    }

    #[test]
    fn tiers_with_enough_claims_get_their_own_estimate() {
        let idx = tiered_index();
        let mut t: HashMap<_, _> = (0..60)
            .map(|i| psm(&format!("s{i}"), 100.0 - i as f64, 0.001, "sp|P1_S|A"))
            .chain((0..60).map(|i| psm(&format!("c{i}"), 100.0 - i as f64, 0.001, "sp|P1_C|A")))
            .chain((0..30).map(|i| {
                psm(
                    &format!("k{i}"),
                    99.5 - 2.0 * i as f64,
                    0.001,
                    "sp|P1_KO1_C|A",
                )
            }))
            .chain((0..5).map(|i| psm(&format!("m{i}"), 100.0 - i as f64, 0.001, "sp|P1_M|A")))
            .collect();
        let counts = apply_knockoff_fdr(&mut t, &idx, 0.05, true).counts;
        let names: Vec<&str> = counts.iter().map(|c| c.stratum.as_str()).collect();
        assert_eq!(
            names,
            [
                "deamidation-like/catalogue",
                "deamidation-like/pooled",
                "deamidation-like/sample"
            ]
        );
        assert!(t["s59"].q_value < 0.01, "{}", t["s59"].q_value);
        assert!(t["c59"].q_value > 0.2, "{}", t["c59"].q_value);
        let pooled = counts
            .iter()
            .find(|c| c.stratum.ends_with("pooled"))
            .unwrap();
        assert_eq!(pooled.tested_before, 5);
    }

    #[test]
    fn one_tier_keeps_class_strata() {
        let idx = index();
        let mut t: HashMap<_, _> = [psm("r", 5.0, 0.001, "sp|P1_V1|A")].into_iter().collect();
        let counts = apply_knockoff_fdr(&mut t, &idx, 0.01, true).counts;
        assert_eq!(counts[0].stratum, "deamidation-like");
    }
}
