//! Protein groups and their q-values, after rescoring.
//!
//! Inference works on protein families: a canonical entry with its isoform,
//! variant and chain entries is one family, so alternative entries never split
//! a protein's evidence. On each side (targets, decoys):
//! 1. every peptide gets the best rescored score of its rows;
//! 2. families with the same peptides form one group;
//! 3. a group whose peptides are a strict subset of another group's is
//!    subsumed: listed with that group, not scored;
//! 4. a peptide still shared between groups goes to the group with the most
//!    peptides (razor), ties to the better best score, then the name;
//! 5. a group scores its best unique or razor peptide.
//!
//! Protein q-values come from picked target-decoy competition: each target
//! family competes with the decoy of the same family, and only the better one
//! counts. Surviving families count once per group, and q = decoys / targets
//! at or above each group's score, with a running minimum.
//!
//! A decoy accession `<prefix>_<accession>` is the decoy of `<accession>`'s
//! family, so pairs are named, not guessed.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use model::{EntryKind, Protein};
use search::SearchIndex;

use crate::percolator::PercolatorPsm;
use crate::proteoform_table::bare_sequence;

/// One protein group.
#[derive(Debug, Clone, PartialEq)]
pub struct ProteinGroup {
    /// 1-based, in order of score.
    pub id: usize,
    /// The reference entry of the group's best family.
    pub lead: String,
    /// The families of the group (indistinguishable by their peptides).
    pub families: Vec<String>,
    /// Every target entry of those families that one of its peptides names.
    pub accessions: Vec<String>,
    /// Peptides of the group with a PSM at its cutoff: all, unique to the
    /// group, and taken from a shared peptide.
    pub peptides: usize,
    pub unique: usize,
    pub razor: usize,
    /// Target PSMs at their cutoff whose peptide counts for the group.
    pub psms: usize,
    pub score: f64,
    pub q_value: f64,
    /// Families whose peptides are all in this group.
    pub subsumed: Vec<String>,
    /// Peptides (with a PSM at its cutoff) that name each accession.
    pub peptides_per_accession: Vec<(String, usize)>,
    pub is_decoy: bool,
}

/// What [`infer_protein_groups`] found.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Inference {
    /// Target groups, best score first.
    pub groups: Vec<ProteinGroup>,
    /// Each target peptide (bare sequence) and the index into `groups` it
    /// counts for.
    pub peptide_group: HashMap<String, usize>,
    /// The target peptides that other kept groups also contain: each counts
    /// for its group as a razor peptide, not a unique one.
    pub shared_peptides: HashSet<String>,
    /// Decoy groups that won their family's competition, with the q-value at
    /// their score; their peptide counts include every assigned peptide.
    pub decoy_groups: Vec<ProteinGroup>,
}

impl Inference {
    /// Each scored family's group q-value. A subsumed family has none.
    pub fn family_q(&self) -> HashMap<&str, f64> {
        self.groups
            .iter()
            .flat_map(|g| g.families.iter().map(move |f| (f.as_str(), g.q_value)))
            .collect()
    }

    /// The group a target row counts for.
    pub fn group_of(&self, psm: &PercolatorPsm) -> Option<&ProteinGroup> {
        self.peptide_group
            .get(&bare_sequence(&psm.peptide))
            .map(|&g| &self.groups[g])
    }
}

/// Where a peptide's quantified features go under a protein q-value cut.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GroupAt<'a> {
    /// The group the peptide counts for, at or below the cut; `unique` is
    /// false for a peptide another group also contains (razor).
    Group {
        group: &'a ProteinGroup,
        unique: bool,
    },
    /// The peptide counts for no target group.
    NoGroup,
    /// The peptide's group is above the cut.
    AboveCut,
}

impl Inference {
    /// The group of the target peptide `sequence` (bare residues) at protein
    /// q-value `max_q`.
    pub fn group_at(&self, sequence: &str, max_q: f64) -> GroupAt<'_> {
        let Some(&i) = self.peptide_group.get(sequence) else {
            return GroupAt::NoGroup;
        };
        let group = &self.groups[i];
        if group.q_value > max_q {
            return GroupAt::AboveCut;
        }
        GroupAt::Group {
            group,
            unique: !self.shared_peptides.contains(sequence),
        }
    }
}

/// One side's groups: families, peptide ids, best score.
struct SideGroup {
    families: Vec<String>,
    peptides: BTreeSet<usize>,
    /// The assigned peptides another kept group also has.
    shared: BTreeSet<usize>,
    score: f64,
    subsumed: Vec<String>,
}

/// One side's peptides and families after parsimony.
struct Side {
    peptides: Vec<String>,
    groups: Vec<SideGroup>,
    /// Family to group index, for every family in a scored group.
    family_group: HashMap<String, usize>,
}

/// Parsimony over one side. `rows` is (peptide, score, families).
fn parsimony(rows: Vec<(String, f64, BTreeSet<String>)>) -> Side {
    // Best score per peptide, and the families each peptide names.
    let mut best: BTreeMap<String, (f64, BTreeSet<String>)> = BTreeMap::new();
    for (pep, score, fams) in rows {
        let e = best
            .entry(pep)
            .or_insert_with(|| (f64::NEG_INFINITY, BTreeSet::new()));
        e.0 = e.0.max(score);
        e.1.extend(fams);
    }
    let peptides: Vec<String> = best.keys().cloned().collect();
    let scores: Vec<f64> = best.values().map(|v| v.0).collect();
    let mut fam_peps: BTreeMap<String, BTreeSet<usize>> = BTreeMap::new();
    for (i, (_, fams)) in best.values().enumerate() {
        for f in fams {
            fam_peps.entry(f.clone()).or_default().insert(i);
        }
    }
    // Families with identical peptide sets.
    let mut by_set: BTreeMap<BTreeSet<usize>, Vec<String>> = BTreeMap::new();
    for (f, s) in fam_peps {
        by_set.entry(s).or_default().push(f);
    }
    let sets: Vec<(BTreeSet<usize>, Vec<String>)> = by_set.into_iter().collect();
    // Strict subsets are subsumed by the largest group containing them.
    let mut groups_of_pep: HashMap<usize, Vec<usize>> = HashMap::new();
    for (g, (s, _)) in sets.iter().enumerate() {
        for &p in s {
            groups_of_pep.entry(p).or_default().push(g);
        }
    }
    let mut subsumed_by: Vec<Option<usize>> = vec![None; sets.len()];
    for (g, (s, _)) in sets.iter().enumerate() {
        let first = *s.iter().next().expect("a family names a peptide");
        subsumed_by[g] = groups_of_pep[&first]
            .iter()
            .copied()
            .filter(|&h| h != g && sets[h].0.len() > s.len() && s.is_subset(&sets[h].0))
            .max_by(|&a, &b| {
                sets[a]
                    .0
                    .len()
                    .cmp(&sets[b].0.len())
                    .then_with(|| sets[b].1.cmp(&sets[a].1))
            });
    }
    // Resolve chains: a subsumed group's parent may itself be subsumed.
    let root = |mut g: usize| {
        while let Some(h) = subsumed_by[g] {
            g = h;
        }
        g
    };
    let kept: Vec<usize> = (0..sets.len())
        .filter(|&g| subsumed_by[g].is_none())
        .collect();
    let best_of = |g: usize| {
        sets[g]
            .0
            .iter()
            .map(|&p| scores[p])
            .fold(f64::NEG_INFINITY, f64::max)
    };
    // Razor: each peptide counts for one kept group.
    let mut owner: HashMap<usize, usize> = HashMap::new();
    let mut shared: HashMap<usize, bool> = HashMap::new();
    for (p, gs) in &groups_of_pep {
        let candidates: Vec<usize> = gs
            .iter()
            .copied()
            .filter(|&g| subsumed_by[g].is_none())
            .collect();
        let Some(&win) = candidates.iter().max_by(|&&a, &&b| {
            sets[a]
                .0
                .len()
                .cmp(&sets[b].0.len())
                .then(best_of(a).total_cmp(&best_of(b)))
                .then_with(|| sets[b].1.cmp(&sets[a].1))
        }) else {
            continue;
        };
        owner.insert(*p, win);
        shared.insert(*p, candidates.len() > 1);
    }
    let mut index_of: HashMap<usize, usize> = HashMap::new();
    let mut groups: Vec<SideGroup> = Vec::new();
    for &g in &kept {
        let assigned: BTreeSet<usize> = sets[g]
            .0
            .iter()
            .copied()
            .filter(|p| owner.get(p) == Some(&g))
            .collect();
        if assigned.is_empty() {
            continue;
        }
        index_of.insert(g, groups.len());
        groups.push(SideGroup {
            families: sets[g].1.clone(),
            shared: assigned.iter().copied().filter(|p| shared[p]).collect(),
            score: assigned
                .iter()
                .map(|&p| scores[p])
                .fold(f64::NEG_INFINITY, f64::max),
            peptides: assigned,
            subsumed: Vec::new(),
        });
    }
    // Subsumed sets, and kept sets that lost every peptide to razor, are
    // listed with the group that took their first peptide.
    for (g, (set, fams)) in sets.iter().enumerate() {
        if index_of.contains_key(&g) {
            continue;
        }
        let first = *set.iter().next().expect("non-empty");
        let host = owner.get(&first).copied().unwrap_or_else(|| root(g));
        if let Some(&i) = index_of.get(&host) {
            groups[i].subsumed.extend(fams.iter().cloned());
        }
    }
    let mut family_group = HashMap::new();
    for (i, g) in groups.iter().enumerate() {
        for f in &g.families {
            family_group.insert(f.clone(), i);
        }
    }
    Side {
        peptides,
        groups,
        family_group,
    }
}

/// Protein groups of the rescored rows, with picked-FDR q-values.
/// `accepted(id)` says whether a target row passes its PSM cutoff; it only
/// counts PSMs per group.
pub fn infer_protein_groups(
    targets: &HashMap<String, PercolatorPsm>,
    decoys: &HashMap<String, PercolatorPsm>,
    idx: &SearchIndex,
    accepted: impl Fn(&str) -> bool,
) -> Inference {
    let proteins: HashMap<&str, &Protein> = idx
        .db
        .proteins
        .iter()
        .map(|p| (p.accession.as_str(), p))
        .collect();
    let family_of = |acc: &str| -> String {
        proteins.get(acc).map_or_else(
            || acc.to_string(),
            |p| {
                p.origin()
                    .map_or_else(|| p.accession.clone(), |o| o.family.to_string())
            },
        )
    };
    let prefix = format!("{}_", idx.decoy_prefix);
    let target_of = |acc: &str| -> String {
        if let Some(t) = acc.strip_prefix(&prefix) {
            return t.to_string();
        }
        match idx.decoy_suffix.as_deref() {
            Some(s) => acc.strip_suffix(s).unwrap_or(acc).to_string(),
            None => acc.to_string(),
        }
    };
    // Knockoff entries exist only to estimate false variant claims: their
    // peptides are evidence for no protein, on either side.
    let knockoff = |acc: &str| proteins.get(acc).is_some_and(|p| p.is_knockoff());
    let rows = |map: &HashMap<String, PercolatorPsm>, decoy: bool| {
        let mut out: Vec<(String, f64, BTreeSet<String>)> = map
            .values()
            .filter_map(|p| {
                let score = p.score?;
                let fams: BTreeSet<String> = p
                    .proteins
                    .split(';')
                    .filter(|a| !a.is_empty())
                    .map(|a| if decoy { target_of(a) } else { a.to_string() })
                    .filter(|a| !knockoff(a))
                    .map(|a| family_of(&a))
                    .collect();
                (!fams.is_empty()).then(|| (bare_sequence(&p.peptide), score, fams))
            })
            .collect();
        out.sort_by(|a, b| a.0.cmp(&b.0).then(b.1.total_cmp(&a.1)));
        out
    };
    let target = parsimony(rows(targets, false));
    let decoy = parsimony(rows(decoys, true));

    // Picked competition per family; the winners count once per group.
    let families: BTreeSet<&String> = target
        .family_group
        .keys()
        .chain(decoy.family_group.keys())
        .collect();
    let mut t_win: BTreeSet<usize> = BTreeSet::new();
    let mut d_win: BTreeSet<usize> = BTreeSet::new();
    for f in families {
        let t = target
            .family_group
            .get(f)
            .map(|&g| (g, target.groups[g].score));
        let d = decoy
            .family_group
            .get(f)
            .map(|&g| (g, decoy.groups[g].score));
        match (t, d) {
            (Some((g, ts)), Some((_, ds))) if ts > ds => {
                t_win.insert(g);
            }
            (Some(_), Some((g, _))) | (None, Some((g, _))) => {
                d_win.insert(g);
            }
            (Some((g, _)), None) => {
                t_win.insert(g);
            }
            (None, None) => {}
        }
    }
    // Ties go to the decoy, and a decoy sorts before a target of equal score.
    let mut ranked: Vec<(f64, bool, usize)> = t_win
        .iter()
        .map(|&g| (target.groups[g].score, false, g))
        .chain(d_win.iter().map(|&g| (decoy.groups[g].score, true, g)))
        .collect();
    ranked.sort_by(|a, b| b.0.total_cmp(&a.0).then(b.1.cmp(&a.1)).then(a.2.cmp(&b.2)));
    let mut fdp = Vec::with_capacity(ranked.len());
    let (mut d, mut t) = (0usize, 0usize);
    for &(_, is_decoy, _) in &ranked {
        if is_decoy {
            d += 1;
        } else {
            t += 1;
        }
        fdp.push((d as f64 / t.max(1) as f64).min(1.0));
    }
    for i in (0..fdp.len().saturating_sub(1)).rev() {
        fdp[i] = fdp[i].min(fdp[i + 1]);
    }
    let mut q_of: HashMap<usize, f64> = HashMap::new();
    let mut decoy_q: Vec<(usize, f64)> = Vec::new();
    for (&(_, is_decoy, g), &q) in ranked.iter().zip(&fdp) {
        if is_decoy {
            decoy_q.push((g, q));
        } else {
            q_of.insert(g, q);
        }
    }

    // PSMs and named entries per target group.
    let mut psms = vec![0usize; target.groups.len()];
    let mut accepted_peps: BTreeSet<usize> = BTreeSet::new();
    let mut named: Vec<BTreeSet<String>> = vec![BTreeSet::new(); target.groups.len()];
    let mut per_accession: Vec<BTreeMap<String, BTreeSet<usize>>> =
        vec![BTreeMap::new(); target.groups.len()];
    let pep_index: HashMap<&str, usize> = target
        .peptides
        .iter()
        .enumerate()
        .map(|(i, p)| (p.as_str(), i))
        .collect();
    let mut pep_owner: HashMap<usize, usize> = HashMap::new();
    for (g, grp) in target.groups.iter().enumerate() {
        for &p in &grp.peptides {
            pep_owner.insert(p, g);
        }
    }
    for (id, p) in targets {
        let Some(&g) = pep_index
            .get(bare_sequence(&p.peptide).as_str())
            .and_then(|i| pep_owner.get(i))
        else {
            continue;
        };
        let pep = pep_index[bare_sequence(&p.peptide).as_str()];
        let ok = accepted(id);
        if ok {
            psms[g] += 1;
            accepted_peps.insert(pep);
        }
        for a in p
            .proteins
            .split(';')
            .filter(|a| !a.is_empty() && !knockoff(a))
        {
            if target.groups[g].families.contains(&family_of(a)) {
                named[g].insert(a.to_string());
                if ok {
                    per_accession[g]
                        .entry(a.to_string())
                        .or_default()
                        .insert(pep);
                }
            }
        }
    }

    let lead_of = |g: usize| -> String {
        let fams = &target.groups[g].families;
        let reference = named[g].iter().find(|a| {
            proteins
                .get(a.as_str())
                .is_some_and(|p| p.entry_kind() == EntryKind::Reference)
        });
        reference
            .or_else(|| named[g].iter().next())
            .cloned()
            .unwrap_or_else(|| fams[0].clone())
    };
    let mut order: Vec<usize> = (0..target.groups.len()).collect();
    order.sort_by(|&a, &b| {
        target.groups[b]
            .score
            .total_cmp(&target.groups[a].score)
            .then_with(|| target.groups[a].families.cmp(&target.groups[b].families))
    });
    let mut new_index = vec![0usize; target.groups.len()];
    let groups: Vec<ProteinGroup> = order
        .iter()
        .enumerate()
        .map(|(i, &g)| {
            new_index[g] = i;
            let s = &target.groups[g];
            ProteinGroup {
                id: i + 1,
                lead: lead_of(g),
                families: s.families.clone(),
                accessions: named[g].iter().cloned().collect(),
                peptides: s.peptides.intersection(&accepted_peps).count(),
                unique: s
                    .peptides
                    .iter()
                    .filter(|p| accepted_peps.contains(p) && !s.shared.contains(p))
                    .count(),
                razor: s.shared.intersection(&accepted_peps).count(),
                psms: psms[g],
                score: s.score,
                // A target that lost its family's competition has no q-value
                // of its own; it is no better than its decoy.
                q_value: q_of.get(&g).copied().unwrap_or(1.0),
                subsumed: {
                    let mut v = s.subsumed.clone();
                    v.sort();
                    v.dedup();
                    v
                },
                peptides_per_accession: per_accession[g]
                    .iter()
                    .map(|(a, peps)| (a.clone(), peps.len()))
                    .collect(),
                is_decoy: false,
            }
        })
        .collect();
    let shared_peptides = target
        .groups
        .iter()
        .flat_map(|g| g.shared.iter().map(|&p| target.peptides[p].clone()))
        .collect();
    let peptide_group = pep_owner
        .into_iter()
        .map(|(p, g)| (target.peptides[p].clone(), new_index[g]))
        .collect();
    // A decoy family `F` stands for the accession `<prefix>_F`.
    let decoy_accession = |f: &str| format!("{prefix}{f}");
    let mut decoy_groups: Vec<ProteinGroup> = decoy_q
        .into_iter()
        .map(|(g, q)| {
            let s = &decoy.groups[g];
            let accessions: Vec<String> = s.families.iter().map(|f| decoy_accession(f)).collect();
            ProteinGroup {
                id: 0,
                lead: accessions[0].clone(),
                families: s.families.clone(),
                peptides_per_accession: accessions
                    .iter()
                    .map(|a| (a.clone(), s.peptides.len()))
                    .collect(),
                accessions,
                peptides: s.peptides.len(),
                unique: s.peptides.len() - s.shared.len(),
                razor: s.shared.len(),
                psms: 0,
                score: s.score,
                q_value: q,
                subsumed: Vec::new(),
                is_decoy: true,
            }
        })
        .collect();
    decoy_groups.sort_by(|a, b| {
        b.score
            .total_cmp(&a.score)
            .then_with(|| a.lead.cmp(&b.lead))
    });
    for (i, g) in decoy_groups.iter_mut().enumerate() {
        g.id = groups.len() + i + 1;
    }
    Inference {
        groups,
        peptide_group,
        shared_peptides,
        decoy_groups,
    }
}

/// Write the protein-group table, best score first.
pub fn write_protein_groups_tsv(path: &std::path::Path, inf: &Inference) -> std::io::Result<()> {
    use std::io::Write;
    crate::atomic::write_atomic(path, |path| {
        let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
        writeln!(
            w,
            "protein_group\tlead\taccessions\tfamilies\tpeptides\tunique_peptides\trazor_peptides\tpsms\tbest_score\tprotein_q\tsubsumed"
        )?;
        for g in &inf.groups {
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}",
                g.id,
                g.lead,
                g.accessions.join(";"),
                g.families.join(";"),
                g.peptides,
                g.unique,
                g.razor,
                g.psms,
                g.score,
                g.q_value,
                g.subsumed.join(";")
            )?;
        }
        w.flush()
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::protein::ProteinDb;

    fn psm(id: &str, score: f64, peptide: &str, proteins: &str) -> (String, PercolatorPsm) {
        (
            id.into(),
            PercolatorPsm {
                psm_id: id.into(),
                q_value: 0.001,
                pep: 0.0,
                score: Some(score),
                peptide: peptide.into(),
                proteins: proteins.into(),
            },
        )
    }

    fn index() -> SearchIndex {
        let p = |acc: &str| model::Protein {
            accession: acc.into(),
            description: String::new(),
            sequence: b"PEPTIDEK".to_vec(),
            annotation: None,
        };
        SearchIndex::from_target_db(
            &ProteinDb {
                proteins: vec![p("A"), p("B"), p("C"), p("D"), p("E")],
            },
            "XXX",
        )
    }

    #[test]
    fn parsimony_merges_subsumes_and_assigns_razor_peptides() {
        let targets: HashMap<String, PercolatorPsm> = [
            // A and B share both peptides: one group.
            psm("1", 5.0, "K.AAAK.R", "A;B"),
            psm("2", 4.0, "K.CCCK.R", "A;B"),
            // C has a subset of D's peptides: subsumed by D.
            psm("3", 3.0, "K.DDDK.R", "C;D"),
            psm("4", 2.0, "K.EEEK.R", "D"),
            // E shares one peptide with D and has its own: the shared one is
            // D's (two peptides before against E's one).
            psm("5", 1.0, "K.FFFK.R", "D;E"),
            psm("6", 0.5, "K.GGGK.R", "E"),
        ]
        .into_iter()
        .collect();
        let inf = infer_protein_groups(&targets, &HashMap::new(), &index(), |_| true);
        assert_eq!(inf.groups.len(), 3);
        let ab = &inf.groups[0];
        assert_eq!(ab.families, vec!["A", "B"]);
        assert_eq!((ab.peptides, ab.unique, ab.razor), (2, 2, 0));
        let d = &inf.groups[1];
        assert_eq!(d.families, vec!["D"]);
        assert_eq!(d.subsumed, vec!["C"]);
        assert_eq!((d.peptides, d.unique, d.razor), (3, 2, 1));
        let e = &inf.groups[2];
        assert_eq!((e.peptides, e.unique, e.razor, e.score), (1, 1, 0, 0.5));
        assert_eq!(inf.group_of(&targets["5"]).unwrap().families, vec!["D"]);
        // Only FFFK is shared between kept groups (C is subsumed).
        assert_eq!(inf.shared_peptides, ["FFFK".to_string()].into());
        match inf.group_at("FFFK", 0.01) {
            GroupAt::Group { group, unique } => {
                assert_eq!(
                    (group.families.as_slice(), unique),
                    (&["D".to_string()][..], false)
                );
            }
            other => panic!("{other:?}"),
        }
        assert!(matches!(
            inf.group_at("EEEK", 0.01),
            GroupAt::Group { unique: true, .. }
        ));
        assert_eq!(inf.group_at("WWWK", 0.01), GroupAt::NoGroup);
        // No decoys: every group has q = 0.
        assert!(inf.groups.iter().all(|g| g.q_value == 0.0));
        // A peptide whose PSM fails its cutoff still places the group but is
        // not counted as its evidence.
        let inf = infer_protein_groups(&targets, &HashMap::new(), &index(), |id| id != "6");
        let e = &inf.groups[2];
        assert_eq!((e.peptides, e.unique, e.psms, e.score), (0, 0, 0, 0.5));
    }

    #[test]
    fn knockoff_entries_are_evidence_for_no_protein() {
        use model::protein::{EntryOrigin, ProteinAnnotation};
        use std::sync::Arc;
        let mut db = ProteinDb {
            proteins: vec![model::Protein {
                accession: "A".into(),
                description: String::new(),
                sequence: b"PEPTIDEK".to_vec(),
                annotation: None,
            }],
        };
        let origin = EntryOrigin {
            kind: EntryKind::Variant,
            family: Arc::from("A"),
            feature: Arc::from("VAR_KO"),
            variant_delta_da: 0.0,
            knockoff_of: Some(Arc::from("VAR_1")),
            source: None,
            tier: Default::default(),
        };
        db.proteins.push(model::Protein {
            accession: "A-KO".into(),
            description: String::new(),
            sequence: b"PEPTLDEK".to_vec(),
            annotation: Some(Arc::new(
                ProteinAnnotation::new(Vec::new(), None).with_origin(Some(origin)),
            )),
        });
        let idx = SearchIndex::from_target_db(&db, "XXX");
        let targets: HashMap<String, PercolatorPsm> = [
            psm("1", 2.0, "K.AAAK.R", "A"),
            // A knockoff-only peptide scores higher than A's own.
            psm("2", 9.0, "K.KKKK.R", "A-KO"),
        ]
        .into_iter()
        .collect();
        let inf = infer_protein_groups(&targets, &HashMap::new(), &idx, |_| true);
        assert_eq!(inf.groups.len(), 1);
        let g = &inf.groups[0];
        assert_eq!((g.peptides, g.score), (1, 2.0));
        assert_eq!(g.accessions, vec!["A"]);
        assert!(inf.group_of(&targets["2"]).is_none());
    }

    #[test]
    fn picked_competition_pairs_a_family_with_its_own_decoy() {
        let targets: HashMap<String, PercolatorPsm> = [
            psm("a", 5.0, "K.AAAK.R", "A"),
            psm("b", 4.0, "K.BBBK.R", "B"),
            psm("c", 1.0, "K.CCCK.R", "C"),
            psm("d", 0.5, "K.DDDK.R", "D"),
        ]
        .into_iter()
        .collect();
        let decoys: HashMap<String, PercolatorPsm> = [
            // B's decoy loses to B: it does not count.
            psm("x", 3.0, "K.XBXK.R", "XXX_B"),
            // C's decoy beats C: C drops out, the decoy counts.
            psm("y", 2.0, "K.XCXK.R", "XXX_C"),
        ]
        .into_iter()
        .collect();
        let inf = infer_protein_groups(&targets, &decoys, &index(), |_| true);
        let q: HashMap<&str, f64> = inf
            .groups
            .iter()
            .map(|g| (g.lead.as_str(), g.q_value))
            .collect();
        // Ranked survivors: A 5, B 4, decoy-C 2, D 0.5.
        assert_eq!(q["A"], 0.0);
        assert_eq!(q["B"], 0.0);
        assert_eq!(q["C"], 1.0, "a target that lost to its decoy");
        assert!((q["D"] - 1.0 / 3.0).abs() < 1e-12);
        // The decoy of C won its competition and is reported as a decoy group.
        assert_eq!(inf.decoy_groups.len(), 1);
        let d = &inf.decoy_groups[0];
        assert_eq!((d.lead.as_str(), d.is_decoy, d.id), ("XXX_C", true, 5));
        assert!((d.q_value - 1.0 / 3.0).abs() < 1e-12);
    }
}
