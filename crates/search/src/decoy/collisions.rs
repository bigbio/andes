//! Peptide-level check of generated decoys against the target peptides.
//!
//! The segment-level redraws in the parent module keep each decoy segment from spelling a target
//! segment, but the search enumerates peptides: missed-cleavage joins of short segments,
//! semi-specific sub-spans and, for a non-specific search, every sub-sequence. This module
//! enumerates those peptides exactly as the search does, for the targets and for each decoy, and
//! redraws the part of a decoy that holds a peptide equal to a target peptide (I and L folded).
//! Peptides that still equal a target after the redraws are counted and logged; the search scores
//! such a peptide as a target (`relabel_collision_decoys_with`), so it never counts as a decoy.

use rayon::prelude::*;

use model::hash::{fnv1a64, splitmix64};
use model::protein::{Protein, ProteinDb};

use super::{
    decoy_protein, draw_segment, folded_hash, shuffle_in_place, DecoyStrategy, Proteases,
    MAX_REDRAWS, MUTATE_POOL,
};
use crate::candidates::{
    compute_cleavage_positions, has_met_cleaved_walk, missed_cleavages_between,
};

/// Rounds of redraws for a decoy that still holds a target peptide.
const REPAIR_ROUNDS: u64 = MAX_REDRAWS;

/// Decoy peptides of one decoy protein that equalled a target peptide before the redraws
/// (`found`) and after them (`left`).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(super) struct Repair {
    pub(super) found: usize,
    pub(super) left: usize,
}

/// Hashes of every target peptide that one or more digests enumerate, I and L folded. The first
/// digest is the one the decoys are built for; the others are the further digests of a cohort
/// that searches the same decoys.
pub(super) struct TargetPeptides {
    hashes: Vec<u64>,
    digests: Vec<Proteases>,
}

impl TargetPeptides {
    pub(super) fn union(db: &ProteinDb, digests: &[Proteases]) -> Self {
        let mut hashes: Vec<u64> = db
            .proteins
            .par_iter()
            .flat_map_iter(|p| {
                let mut out = Vec::new();
                for d in digests {
                    for_each_span(&p.sequence, has_met_cleaved_walk(p), d, |s, e| {
                        out.push(peptide_hash(&p.sequence[s..e]));
                    });
                }
                out
            })
            .collect();
        hashes.par_sort_unstable();
        hashes.dedup();
        Self {
            hashes,
            digests: digests.to_vec(),
        }
    }

    pub(super) fn len(&self) -> usize {
        self.hashes.len()
    }

    fn contains(&self, residues: &[u8]) -> bool {
        self.hashes.binary_search(&peptide_hash(residues)).is_ok()
    }

    /// Whether a residue is a cleavage site of any of the digests.
    fn makes_cut(&self, r: u8) -> bool {
        self.digests.iter().any(|d| d.makes_cut(r))
    }

    /// The `[start, end)` spans of the peptides any of the digests enumerates from `decoy` that
    /// equal a target peptide: the build digest's in its enumeration order, then the spans only
    /// a further digest enumerates.
    pub(super) fn collisions(&self, decoy: &Protein) -> Vec<(usize, usize)> {
        let seq = &decoy.sequence;
        let met = has_met_cleaved_walk(decoy);
        let mut out = Vec::new();
        for (k, d) in self.digests.iter().enumerate() {
            let before = out.len();
            for_each_span(seq, met, d, |s, e| {
                if self.contains(&seq[s..e]) && (k == 0 || !out[..before].contains(&(s, e))) {
                    out.push((s, e));
                }
            });
        }
        out
    }
}

/// Distinct peptides of one length that a digest enumerates from the targets and from the
/// decoys, and how many of the decoy ones equal a target one (I and L folded).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LengthOverlap {
    pub length: u32,
    pub target_peptides: u64,
    pub decoy_peptides: u64,
    pub shared_peptides: u64,
}

fn length_keyed_peptides(db: &ProteinDb, digest: &Proteases) -> Vec<(u32, u64)> {
    let mut keys: Vec<(u32, u64)> = db
        .proteins
        .par_iter()
        .flat_map_iter(|p| {
            let mut out = Vec::new();
            for_each_span(&p.sequence, has_met_cleaved_walk(p), digest, |s, e| {
                out.push(((e - s) as u32, peptide_hash(&p.sequence[s..e])));
            });
            out
        })
        .collect();
    keys.par_sort_unstable();
    keys.dedup();
    keys
}

/// Per peptide length, the distinct target and decoy peptides `digest` enumerates and the decoy
/// peptides equal to a target peptide (I and L folded), shortest length first.
pub(super) fn peptide_overlap(
    targets: &ProteinDb,
    decoys: &ProteinDb,
    digest: &Proteases,
) -> Vec<LengthOverlap> {
    let t = length_keyed_peptides(targets, digest);
    let d = length_keyed_peptides(decoys, digest);
    let mut counts: std::collections::BTreeMap<u32, (u64, u64, u64)> = Default::default();
    for &(len, _) in &t {
        counts.entry(len).or_default().0 += 1;
    }
    for &(len, h) in &d {
        let c = counts.entry(len).or_default();
        c.1 += 1;
        c.2 += u64::from(t.binary_search(&(len, h)).is_ok());
    }
    counts
        .into_iter()
        .map(|(length, (t, d, s))| LengthOverlap {
            length,
            target_peptides: t,
            decoy_peptides: d,
            shared_peptides: s,
        })
        .collect()
}

fn peptide_hash(residues: &[u8]) -> u64 {
    splitmix64(folded_hash(residues.iter().copied()))
}

/// Calls `f(start, end)` for every peptide span the search enumerates from a protein with
/// sequence `seq` (`candidates::enumerate_protein`): cleavage sites of the enzymes, the
/// specificity, the length and missed-cleavage limits, the initiator-Met-cleaved walk, and no
/// span holding a residue without a standard mass. Variant and chain windows are not applied, so
/// the spans are a superset of the search's for annotated entries.
pub(super) fn for_each_span(
    seq: &[u8],
    met_walk: bool,
    proteases: &Proteases,
    mut f: impl FnMut(usize, usize),
) {
    let mut unknown = Vec::with_capacity(seq.len() + 1);
    unknown.push(0u32);
    for &r in seq {
        let last = unknown[unknown.len() - 1];
        unknown.push(last + u32::from(model::AminoAcid::standard(r).is_none()));
    }
    let mut emit = |s: usize, e: usize| {
        if unknown[e] == unknown[s] {
            f(s, e);
        }
    };
    spans_from(seq, 0, met_walk, proteases, &mut emit);
    if met_walk && !seq.is_empty() {
        spans_from(seq, 1, false, proteases, &mut emit);
    }
}

/// Spans of `seq[offset..]` (`enumerate_protein_from_offset`); `offset` 1 is the walk after a
/// cleaved initiator Met, which only starts at its first residue. `skip_second` drops the spans
/// starting at position 1 of the offset-0 walk, which the offset-1 walk emits.
fn spans_from(
    seq: &[u8],
    offset: usize,
    skip_second: bool,
    proteases: &Proteases,
    emit: &mut impl FnMut(usize, usize),
) {
    let sub = &seq[offset..];
    let n = sub.len() as u32;
    if n < proteases.min_length {
        return;
    }
    let min = proteases.min_length.max(1);
    let max = proteases.max_length;
    let missed = proteases.max_missed_cleavages;
    let n_term_only = offset > 0;
    let mut out = |s: u32, e: u32| {
        let (s, e) = (s as usize + offset, e as usize + offset);
        if !(skip_second && s == 1) {
            emit(s, e);
        }
    };
    if proteases.is_non_specific() {
        let starts = if n_term_only { 0..n.min(1) } else { 0..n };
        for s in starts {
            for e in (s + min)..=s.saturating_add(max).min(n) {
                out(s, e);
            }
        }
        return;
    }
    let cuts = compute_cleavage_positions(sub, proteases.primary, &proteases.extra);
    for (i, &start) in cuts.iter().enumerate() {
        if n_term_only && start != 0 {
            break;
        }
        for (k, &end) in cuts[i + 1..].iter().enumerate() {
            let len = end - start;
            if len > max {
                break;
            }
            if len < min || k as u32 > missed {
                continue;
            }
            out(start, end);
        }
    }
    if proteases.num_tolerable_termini != 1 {
        return;
    }
    let is_cut = |p: u32| cuts.binary_search(&p).is_ok();
    for &start in &cuts {
        if n_term_only && start != 0 {
            break;
        }
        for end in (start + min)..=start.saturating_add(max).min(n) {
            if is_cut(end) || missed_cleavages_between(&cuts, start, end) > missed {
                continue;
            }
            out(start, end);
        }
    }
    if n_term_only {
        return;
    }
    for &end in &cuts {
        if end < min {
            continue;
        }
        for start in end.saturating_sub(max)..=(end - min) {
            if is_cut(start) || missed_cleavages_between(&cuts, start, end) > missed {
                continue;
            }
            out(start, end);
        }
    }
}

fn repair_seed(seed: u64, residues: impl Iterator<Item = u8>, round: u64) -> u64 {
    splitmix64(fnv1a64(residues) ^ splitmix64(seed ^ (MAX_REDRAWS + 1 + round)))
}

/// Redraws the segments of a peptide-level decoy that hold a peptide equal to a target peptide.
/// Each round redraws one movable segment of every colliding peptide (a different one each
/// round), with draw `MAX_REDRAWS + 1 + round` of that segment, so composition, cleavage sites
/// and terminal residues of every segment stay as they are.
pub(super) fn repair_peptide_level(
    p: &Protein,
    prefix: &str,
    mut source: Vec<u32>,
    strategy: DecoyStrategy,
    seed: u64,
    proteases: &Proteases,
    peptides: &TargetPeptides,
) -> (Protein, Repair) {
    let mut decoy = decoy_protein(p, prefix, source.clone(), proteases);
    let mut hits = peptides.collisions(&decoy);
    let found = hits.len();
    if found == 0 {
        return (decoy, Repair::default());
    }
    let seq = &p.sequence;
    let met = has_met_cleaved_walk(p);
    let cuts = compute_cleavage_positions(seq, proteases.primary, &proteases.extra);
    let segment = |j: usize| (cuts[j] as usize, cuts[j + 1] as usize);
    let draws = |j: usize, k: u64| {
        let (a, b) = segment(j);
        draw_segment(&seq[a..b], met && a == 0, strategy, seed, proteases, k)
    };
    let movable: Vec<bool> = (0..cuts.len() - 1).map(|j| draws(j, 1).is_some()).collect();
    for round in 0..REPAIR_ROUNDS {
        if hits.is_empty() {
            break;
        }
        let mut redrawn = vec![false; movable.len()];
        for &(s, e) in &hits {
            let lo = cuts.partition_point(|&c| c as usize <= s) - 1;
            let hi = cuts.partition_point(|&c| (c as usize) < e);
            let open: Vec<usize> = (lo..hi).filter(|&j| movable[j]).collect();
            if open.is_empty() {
                continue;
            }
            let j = open[(round as usize + s) % open.len()];
            if redrawn[j] {
                continue;
            }
            redrawn[j] = true;
            if let Some(local) = draws(j, MAX_REDRAWS + 1 + round) {
                let a = segment(j).0;
                for (k, from) in local.into_iter().enumerate() {
                    source[a + k] = a as u32 + from;
                }
            }
        }
        decoy = decoy_protein(p, prefix, source.clone(), proteases);
        hits = peptides.collisions(&decoy);
    }
    let left = hits.len();
    (decoy, Repair { found, left })
}

/// Redraws the windows of a whole-sequence decoy (non-specific search) that hold a peptide equal
/// to a target peptide: the residues of the peptide are shuffled among its positions, and each
/// round widens the window by one residue on each side.
pub(super) fn repair_whole_sequence(
    p: &Protein,
    prefix: &str,
    mut source: Vec<u32>,
    seed: u64,
    proteases: &Proteases,
    peptides: &TargetPeptides,
) -> (Protein, Repair) {
    let mut decoy = decoy_protein(p, prefix, source.clone(), proteases);
    let mut hits = peptides.collisions(&decoy);
    let found = hits.len();
    if found == 0 {
        return (decoy, Repair::default());
    }
    let n = source.len();
    for round in 0..REPAIR_ROUNDS {
        if hits.is_empty() {
            break;
        }
        let mut touched = vec![false; n];
        for &(s, e) in &hits {
            let (a, b) = (s.saturating_sub(round as usize), (e + round as usize).min(n));
            if touched[a..b].iter().any(|&t| t) {
                continue;
            }
            let mut window = source[a..b].to_vec();
            let state = repair_seed(seed, window.iter().map(|&i| p.sequence[i as usize]), round);
            shuffle_in_place(&mut window, state);
            source[a..b].copy_from_slice(&window);
            touched[a..b].fill(true);
        }
        decoy = decoy_protein(p, prefix, source.clone(), proteases);
        hits = peptides.collisions(&decoy);
    }
    let left = hits.len();
    (decoy, Repair { found, left })
}

/// Adds a substitution to every peptide of a `mutate` decoy that still equals a target peptide,
/// such as a missed-cleavage join of segments too short to be mutated: an interior residue of the
/// peptide (not a cleavage residue, an annotated site or the residue after an initiator Met)
/// takes a residue of `MUTATE_POOL` whose mass differs from the target residue's by more than
/// 1.1 Da and that makes no cleavage site.
pub(super) fn repair_mutated(
    p: &Protein,
    mut decoy: Protein,
    seed: u64,
    peptides: &TargetPeptides,
) -> (Protein, Repair) {
    let mut hits = peptides.collisions(&decoy);
    let found = hits.len();
    if found == 0 {
        return (decoy, Repair::default());
    }
    let mass = |r: u8| model::AminoAcid::standard(r).map(|a| a.mass);
    let met = has_met_cleaved_walk(p);
    let sites: Vec<usize> = p.mod_sites().iter().map(|s| s.position as usize).collect();
    let n = decoy.sequence.len();
    for round in 0..REPAIR_ROUNDS {
        if hits.is_empty() {
            break;
        }
        let mut changed = vec![false; n];
        for &(s, e) in &hits {
            if changed[s..e].iter().any(|&c| c) {
                continue;
            }
            let seq = &decoy.sequence;
            let open: Vec<usize> = (s + 1..e.saturating_sub(1))
                .filter(|&i| {
                    let fixed = (met && i == 1)
                        || sites.contains(&i)
                        || peptides.makes_cut(seq[i])
                        || mass(p.sequence[i]).is_none();
                    !fixed
                })
                .collect();
            if open.is_empty() {
                continue;
            }
            let h = repair_seed(seed, seq[s..e].iter().copied(), round);
            let at = open[(h % open.len() as u64) as usize];
            let Some(old) = mass(p.sequence[at]) else {
                continue;
            };
            let pool: Vec<u8> = MUTATE_POOL
                .iter()
                .copied()
                .filter(|&r| {
                    r != seq[at]
                        && mass(r).is_some_and(|m| (m - old).abs() > 1.1)
                        && !peptides.makes_cut(r)
                })
                .collect();
            if pool.is_empty() {
                continue;
            }
            decoy.sequence[at] = pool[(splitmix64(h) % pool.len() as u64) as usize];
            changed[at] = true;
        }
        hits = peptides.collisions(&decoy);
    }
    let left = hits.len();
    (decoy, Repair { found, left })
}

fn thousands(n: usize) -> String {
    let s = n.to_string();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// One line with the outcome of the check over all decoy proteins.
pub(super) fn log_repairs(repairs: &[Repair], target_peptides: usize, secs: f64) {
    let found: usize = repairs.iter().map(|r| r.found).sum();
    let left: usize = repairs.iter().map(|r| r.left).sum();
    let proteins = repairs.iter().filter(|r| r.found > 0).count();
    if found == 0 {
        log::info!(target: "database",
            "decoy check: no decoy peptide equals one of the {} target peptides (I = L); {secs:.2} s",
            thousands(target_peptides)
        );
        return;
    }
    log::info!(target: "database",
        "decoy check: {} decoy peptides in {} decoy proteins equalled a target peptide (I = L) \
         and were redrawn; {} still do and are scored as targets; {secs:.2} s",
        thousands(found),
        thousands(proteins),
        thousands(left)
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::decoy::DEFAULT_DECOY_SEED;
    use model::enzyme::Enzyme;
    use std::collections::BTreeSet;

    fn random_db(count: usize) -> ProteinDb {
        let aa = b"ACDEFGHIKLMNPQRSTVWYKRDX";
        let mut state = 0x0DDB_A11C_0FFE_E123u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        ProteinDb {
            proteins: (0..count)
                .map(|i| {
                    let len = (next() % 120) as usize;
                    let mut sequence: Vec<u8> =
                        (0..len).map(|_| aa[(next() % aa.len() as u64) as usize]).collect();
                    if i % 3 == 0 && !sequence.is_empty() {
                        sequence[0] = b'M';
                    }
                    Protein {
                        accession: format!("P{i}"),
                        description: String::new(),
                        sequence,
                        annotation: None,
                    }
                })
                .collect(),
        }
    }

    #[test]
    fn spans_are_the_peptides_the_search_enumerates() {
        let db = random_db(150);
        let base = Proteases {
            min_length: 5,
            max_length: 25,
            ..Proteases::default()
        };
        let digests = [
            base.clone(),
            Proteases {
                num_tolerable_termini: 1,
                max_missed_cleavages: 1,
                ..base.clone()
            },
            Proteases {
                primary: Enzyme::NonSpecific,
                num_tolerable_termini: 0,
                min_length: 8,
                max_length: 12,
                ..base.clone()
            },
            Proteases {
                num_tolerable_termini: 0,
                ..base.clone()
            },
            Proteases {
                primary: Enzyme::AspN,
                ..base.clone()
            },
            Proteases {
                extra: vec![Enzyme::GluC],
                max_missed_cleavages: 0,
                ..base.clone()
            },
            Proteases {
                primary: Enzyme::LysN,
                num_tolerable_termini: 1,
                ..base.clone()
            },
        ];
        let aa = model::aa_set::AminoAcidSetBuilder::new_standard()
            .build()
            .expect("standard amino acids");
        for proteases in &digests {
            let mut params = crate::params::SearchParams::default_tryptic(aa.clone());
            params.enzyme = proteases.primary;
            params.extra_enzymes = proteases.extra.clone();
            params.num_tolerable_termini = proteases.num_tolerable_termini;
            params.min_length = proteases.min_length;
            params.max_length = proteases.max_length;
            params.max_missed_cleavages = proteases.max_missed_cleavages;
            params.max_variable_mods_per_peptide = 0;
            let idx = crate::index::SearchIndex::from_target_db_for(
                &db,
                "NONE",
                None,
                DecoyStrategy::None,
                DEFAULT_DECOY_SEED,
                proteases,
            );
            let want: BTreeSet<(usize, usize, usize)> =
                crate::candidates::enumerate_candidates(&idx, &params, "NONE")
                    .map(|c| {
                        let s = c.start_offset_in_protein;
                        (c.protein_index, s, s + c.peptide.residues.len())
                    })
                    .collect();
            let mut got = BTreeSet::new();
            for (i, p) in db.proteins.iter().enumerate() {
                for_each_span(&p.sequence, has_met_cleaved_walk(p), proteases, |s, e| {
                    assert!(got.insert((i, s, e)), "{proteases:?} span {i} {s}..{e} twice");
                });
            }
            assert!(!want.is_empty(), "{proteases:?}");
            assert_eq!(got, want, "{proteases:?}");
        }
    }

    #[test]
    fn thousands_groups_digits() {
        assert_eq!(thousands(0), "0");
        assert_eq!(thousands(999), "999");
        assert_eq!(thousands(1_000), "1,000");
        assert_eq!(thousands(1_234_567), "1,234,567");
    }
}
