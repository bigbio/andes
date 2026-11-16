//! Roll-up of crosslink-spectrum matches to residue pairs and protein pairs,
//! with a class-split `(TD − DD) / TT` FDR at each level.
//!
//! A level groups the items of the level below by a key and keeps the best
//! score of each group. Targets and decoys are grouped by the same rule (the
//! key holds the decoy accessions of decoy peptides, so a group never mixes
//! classes), and q-values are computed separately for intra-protein, homotypic and
//! inter-protein groups.

use std::collections::BTreeMap;

use crate::decoy::PairClass;
use crate::fdr::{q_values, q_values_homotypic};

/// Link class of a crosslink, the unit the q-values are computed in.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum LinkClass {
    /// Two different peptide sequences of one protein (family).
    Intra,
    /// The same peptide sequence on both sides: a homodimer or a random peptide paired
    /// with itself, whose false matches are TT or DD but never TD.
    Homotypic,
    /// Peptides of different proteins.
    Inter,
}

impl LinkClass {
    pub const ALL: [LinkClass; 3] = [LinkClass::Intra, LinkClass::Homotypic, LinkClass::Inter];

    pub fn of(intra: bool, homotypic: bool) -> Self {
        match (intra, homotypic) {
            (false, _) => LinkClass::Inter,
            (true, true) => LinkClass::Homotypic,
            (true, false) => LinkClass::Intra,
        }
    }

    pub fn is_intra(self) -> bool {
        self != LinkClass::Inter
    }

    /// Output label.
    pub fn name(self) -> &'static str {
        match self {
            LinkClass::Intra => "intra",
            LinkClass::Homotypic => "homotypic",
            LinkClass::Inter => "inter",
        }
    }
}

/// One item of a level: a CSM, a residue pair or a protein pair.
#[derive(Debug, Clone, PartialEq)]
pub struct Item {
    /// Group key at this level.
    pub key: String,
    /// Score (higher is better).
    pub score: f64,
    pub class: PairClass,
    pub link: LinkClass,
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
    pub link: LinkClass,
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
            link: it.link,
            n_members: 0,
            n_csms: 0,
            best: i,
        });
        g.n_members += 1;
        g.n_csms += it.n_csms;
        if it.score > g.score {
            g.score = it.score;
            g.class = it.class;
            g.link = it.link;
            g.best = i;
        }
    }
    groups.into_values().collect()
}

/// Class-split q-values of `groups`: computed separately among the intra, the
/// homotypic and the inter groups, over the groups that `eligible` accepts (`None`
/// for the others).
pub fn split_q_values<F: Fn(&Group) -> bool>(groups: &[Group], eligible: F) -> Vec<Option<f64>> {
    let mut q = vec![None; groups.len()];
    for link in LinkClass::ALL {
        let rows: Vec<usize> = (0..groups.len())
            .filter(|&i| groups[i].link == link && eligible(&groups[i]))
            .collect();
        let scores: Vec<f64> = rows.iter().map(|&i| groups[i].score).collect();
        let classes: Vec<PairClass> = rows.iter().map(|&i| groups[i].class).collect();
        let qv = if link == LinkClass::Homotypic {
            q_values_homotypic(&scores, &classes)
        } else {
            q_values(&scores, &classes)
        };
        for (k, v) in qv.into_iter().enumerate() {
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
            link: g.link,
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
            link: LinkClass::of(intra, false),
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
    fn homotypic_groups_are_a_class_of_their_own() {
        let mut items: Vec<Item> = (0..300)
            .map(|i| item(&format!("h{i}"), 1000.0 - i as f64, TT, true))
            .collect();
        for it in &mut items {
            it.link = LinkClass::Homotypic;
        }
        for k in 0..5 {
            let mut dd = item(&format!("hd{k}"), 900.0 - k as f64, DD, true);
            dd.link = LinkClass::Homotypic;
            items.push(dd);
        }
        items.extend((0..300).map(|i| item(&format!("i{i}"), 1000.0 - i as f64, TT, true)));
        items.extend((0..5).map(|k| item(&format!("id{k}"), 900.0 - k as f64, DD, true)));
        let g = group_by(&items, |_, it| it.key.clone());
        let q = split_q_values(&g, |_| true);
        let q_of = |k: &str| q[g.iter().position(|g| g.key == k).unwrap()].unwrap();
        // A homotypic DD stands for a false homotypic TT ...
        assert!((q_of("h299") - 6.0 / 300.0).abs() < 1e-12);
        // ... while among the other intra links it offsets a TD counted twice.
        assert!((q_of("i299") - 2.0 / 300.0).abs() < 1e-12);
        assert_eq!(LinkClass::of(true, true).name(), "homotypic");
        assert!(LinkClass::Homotypic.is_intra() && !LinkClass::Inter.is_intra());
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
