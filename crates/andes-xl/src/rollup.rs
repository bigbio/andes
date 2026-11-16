//! Roll-up of crosslink-spectrum matches to residue pairs and protein pairs,
//! with a class-split `(TD − DD) / TT` FDR at each level.
//!
//! A level groups the items of the level below by a key and keeps the best
//! score of each group. Targets and decoys are grouped by the same rule (the
//! key holds the decoy accessions of decoy peptides, so a group never mixes
//! classes), and q-values are computed separately for intra- and
//! inter-protein groups.

use std::collections::BTreeMap;

use crate::decoy::PairClass;
use crate::fdr::q_values;

/// One item of a level: a CSM, a residue pair or a protein pair.
#[derive(Debug, Clone, PartialEq)]
pub struct Item {
    /// Group key at this level.
    pub key: String,
    /// Score (higher is better).
    pub score: f64,
    pub class: PairClass,
    pub intra: bool,
    /// Number of CSMs behind the item (1 for a CSM).
    pub n_csms: usize,
}

/// A group of items with its best score.
#[derive(Debug, Clone, PartialEq)]
pub struct Group {
    pub key: String,
    /// Best member score.
    pub score: f64,
    pub class: PairClass,
    pub intra: bool,
    /// Members (items of the level below).
    pub n_members: usize,
    /// CSMs behind the members.
    pub n_csms: usize,
    /// Index of the best member in the input.
    pub best: usize,
}

/// Group `items` by `key_of(item)`; each group takes its best member's score,
/// class and link class. Groups come out sorted by key.
pub fn group_by<F: Fn(usize, &Item) -> String>(items: &[Item], key_of: F) -> Vec<Group> {
    let mut groups: BTreeMap<String, Group> = BTreeMap::new();
    for (i, it) in items.iter().enumerate() {
        let key = key_of(i, it);
        let g = groups.entry(key.clone()).or_insert_with(|| Group {
            key,
            score: it.score,
            class: it.class,
            intra: it.intra,
            n_members: 0,
            n_csms: 0,
            best: i,
        });
        g.n_members += 1;
        g.n_csms += it.n_csms;
        if it.score > g.score {
            g.score = it.score;
            g.class = it.class;
            g.intra = it.intra;
            g.best = i;
        }
    }
    groups.into_values().collect()
}

/// Class-split q-values of `groups`: computed among the intra groups and
/// among the inter groups separately, over the groups that `eligible`
/// accepts (`None` for the others).
pub fn split_q_values<F: Fn(&Group) -> bool>(groups: &[Group], eligible: F) -> Vec<Option<f64>> {
    let mut q = vec![None; groups.len()];
    for intra in [true, false] {
        let rows: Vec<usize> = (0..groups.len())
            .filter(|&i| groups[i].intra == intra && eligible(&groups[i]))
            .collect();
        let scores: Vec<f64> = rows.iter().map(|&i| groups[i].score).collect();
        let classes: Vec<PairClass> = rows.iter().map(|&i| groups[i].class).collect();
        for (k, v) in q_values(&scores, &classes).into_iter().enumerate() {
            q[rows[k]] = Some(v);
        }
    }
    q
}

/// The groups as items of the next level up.
pub fn as_items(groups: &[Group], key_of: impl Fn(&Group) -> String) -> Vec<Item> {
    groups
        .iter()
        .map(|g| Item {
            key: key_of(g),
            score: g.score,
            class: g.class,
            intra: g.intra,
            n_csms: g.n_csms,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use PairClass::*;

    fn item(key: &str, score: f64, class: PairClass, intra: bool) -> Item {
        Item {
            key: key.to_string(),
            score,
            class,
            intra,
            n_csms: 1,
        }
    }

    #[test]
    fn groups_keep_the_best_member_and_count_csms() {
        let items = vec![
            item("a", 1.0, TT, true),
            item("b", 5.0, TD, false),
            item("a", 3.0, TT, true),
        ];
        let g = group_by(&items, |_, it| it.key.clone());
        assert_eq!(g.len(), 2);
        assert_eq!(
            (g[0].key.as_str(), g[0].score, g[0].n_members, g[0].best),
            ("a", 3.0, 2, 2)
        );
        assert_eq!((g[1].class, g[1].n_csms), (TD, 1));
        let up = as_items(&g, |_| "p".to_string());
        assert_eq!(up[0].n_csms, 2);
    }

    #[test]
    fn q_values_are_split_by_link_class_and_respect_eligibility() {
        let items = vec![
            item("i1", 10.0, TT, true),
            item("i2", 9.0, TD, true),
            item("x1", 8.0, TT, false),
            item("x2", 7.0, TT, false),
        ];
        let g = group_by(&items, |_, it| it.key.clone());
        let q = split_q_values(&g, |_| true);
        assert_eq!(q, vec![Some(1.0), Some(1.0), Some(1.0), Some(1.0)]);
        let q = split_q_values(&g, |g| g.class == TT);
        assert_eq!(q[1], None);
        assert_eq!(q[0], Some(1.0));
    }
}
