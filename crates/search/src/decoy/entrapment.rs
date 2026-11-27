//! Shuffled entrapment proteins for a search that measures its own false discovery proportion.
//!
//! The construction is the peptide-level shuffled entrapment of Wen et al., Nat. Methods 2025
//! [cite:wen2025_entrapment]: each target protein gets one entrapment protein whose every
//! segment between two cleavage sites is a shuffle of the target segment, with the cleavage
//! residues, the cleavage-side terminal residue and an initiator Met in place. Every enzymatic
//! entrapment peptide therefore has a target peptide of the same length, composition and mass
//! at the same position of its source protein (the pair the paired estimator needs). The
//! shuffle is the decoy shuffle of [`super::DecoyStrategy::Shuffle`] with its own seed, and each
//! entrapment peptide that equals a target peptide (I and L folded) is redrawn as a decoy is.
//! The entrapment proteins are added to the targets before decoys are built, so they get decoys
//! of their own and the search cannot tell them from targets.

use std::collections::{BTreeSet, HashSet};

use aho_corasick::AhoCorasick;
use rayon::prelude::*;

use rustc_hash::FxHashMap;

use model::hash::splitmix64;
use model::protein::{EntryKind, Protein, ProteinDb};

use super::collisions::{repair_peptide_level, TargetPeptides};
use super::{peptide_level_permutation, DecoyStrategy, Proteases, TargetSegments};

/// Accession prefix of an entrapment protein, without its separator: `ENTRAP_<accession>`, the
/// tag the benchmark evaluators under `docs/benchmarks/` read.
pub const ENTRAPMENT_PREFIX: &str = "ENTRAP";

/// Default `--entrapment-seed`.
pub const DEFAULT_ENTRAPMENT_SEED: u64 = 1;

/// Mixed into the entrapment seed so the entrapment shuffle of a segment is never the decoy
/// shuffle of the same segment, whatever the two seeds are.
const ENTRAPMENT_SEED_SALT: u64 = 0xE47A_9C1D_5B3F_2E61;

/// What building the entrapment proteins did.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct EntrapmentStats {
    /// Entrapment proteins built, one per eligible target.
    pub proteins: usize,
    /// Target entries without an entrapment protein: variant, isoform, chain and knockoff
    /// entries, whose peptides the reference entry already pairs.
    pub skipped: usize,
    /// Entrapment peptides that equalled a target peptide and were redrawn.
    pub redrawn: usize,
    /// Entrapment peptides that still equal a target peptide; a match to one maps to a target
    /// protein too and is never counted as entrapment.
    pub left: usize,
    /// Distinct target peptides of the digest, I and L folded.
    pub target_peptides: usize,
    /// Distinct entrapment peptides that are no target peptide, I and L folded.
    pub entrapment_peptides: usize,
}

impl EntrapmentStats {
    /// The entrapment to target ratio `r` of the combined estimator, in peptide space.
    pub fn peptide_ratio(&self) -> f64 {
        self.entrapment_peptides as f64 / self.target_peptides.max(1) as f64
    }
}

/// The entrapment proteins of a target database.
pub struct Entrapment {
    /// `ENTRAP_<accession>` proteins, in the order of their targets.
    pub proteins: ProteinDb,
    pub stats: EntrapmentStats,
}

/// Builds one shuffled entrapment protein for every reference target of `targets` (see the
/// module doc), for the `proteases` digest and seeded by `seed`. An entrapment protein keeps the
/// annotation of its target, its sites moved with its residues, with its own accession as its
/// protein family. `None` for a non-specific digest, which has no segments to shuffle.
pub fn entrapment_db(targets: &ProteinDb, seed: u64, proteases: &Proteases) -> Option<Entrapment> {
    if proteases.is_non_specific() {
        return None;
    }
    let digests = std::slice::from_ref(proteases);
    let peptides = TargetPeptides::union(targets, digests);
    let segments = TargetSegments::new(targets, proteases);
    let shuffle_seed = splitmix64(seed ^ ENTRAPMENT_SEED_SALT);
    let built: Vec<_> = targets
        .proteins
        .par_iter()
        .filter(|p| p.entry_kind() == EntryKind::Reference && !p.is_knockoff())
        .map(|p| {
            let source = peptide_level_permutation(
                p,
                DecoyStrategy::Shuffle,
                shuffle_seed,
                proteases,
                &segments,
            );
            repair_peptide_level(
                p,
                ENTRAPMENT_PREFIX,
                source,
                DecoyStrategy::Shuffle,
                shuffle_seed,
                proteases,
                &peptides,
            )
        })
        .collect();
    let mut stats = EntrapmentStats {
        proteins: built.len(),
        skipped: targets.proteins.len() - built.len(),
        target_peptides: peptides.len(),
        ..Default::default()
    };
    let mut proteins = Vec::with_capacity(built.len());
    for (mut p, repair) in built {
        stats.redrawn += repair.found;
        stats.left += repair.left;
        if let Some(a) = p.annotation.as_mut() {
            if let Some(o) = std::sync::Arc::make_mut(a).origin.as_mut() {
                o.family = std::sync::Arc::from(p.accession.as_str());
            }
        }
        proteins.push(p);
    }
    let proteins = ProteinDb { proteins };
    stats.entrapment_peptides = TargetPeptides::union(&proteins, digests).count_outside(&peptides);
    Some(Entrapment { proteins, stats })
}

/// What reading a searched database's entrapment needs: the original target peptides (I and L
/// folded), each protein by accession, and the entrapment to target ratio of the digest.
pub struct EntrapmentLookup<'a> {
    targets: TargetPeptides,
    proteins: FxHashMap<&'a str, &'a Protein>,
    /// Distinct original target peptides of the digest.
    pub target_peptides: usize,
    /// Distinct entrapment peptides that are no target peptide.
    pub entrapment_peptides: usize,
}

impl<'a> EntrapmentLookup<'a> {
    /// The lookup of `searched`, the target entries of a search (original and entrapment, no
    /// decoys), for the `proteases` digest. `None` when it holds no entrapment protein.
    pub fn new(searched: impl Iterator<Item = &'a Protein>, proteases: &Proteases) -> Option<Self> {
        let searched: Vec<&'a Protein> = searched.collect();
        let (ent, orig): (Vec<&Protein>, Vec<&Protein>) = searched
            .iter()
            .copied()
            .partition(|p| is_entrapment_accession(&p.accession));
        if ent.is_empty() {
            return None;
        }
        let own = |v: Vec<&Protein>| ProteinDb {
            proteins: v.into_iter().cloned().collect(),
        };
        let digests = std::slice::from_ref(proteases);
        let targets = TargetPeptides::union(&own(orig), digests);
        let entrapment_peptides = TargetPeptides::union(&own(ent), digests).count_outside(&targets);
        Some(Self {
            target_peptides: targets.len(),
            entrapment_peptides,
            targets,
            proteins: searched
                .into_iter()
                .map(|p| (p.accession.as_str(), p))
                .collect(),
        })
    }

    /// The entrapment to target ratio `r` in peptide space.
    pub fn ratio(&self) -> f64 {
        self.entrapment_peptides as f64 / self.target_peptides.max(1) as f64
    }

    /// Whether `residues` (unfolded) is an original target peptide, with I and L folded.
    pub fn is_target_peptide(&self, residues: &[u8]) -> bool {
        self.targets.contains(residues)
    }

    /// The target peptide `peptide` of the entrapment protein `accession` was shuffled from: the
    /// residues of its source protein at every place `peptide` occurs. `None` when the
    /// accession is no entrapment protein, its source is not searched, the peptide is absent,
    /// or repeated occurrences map to different target sequences.
    pub fn paired_target(&self, accession: &str, peptide: &[u8]) -> Option<&'a [u8]> {
        let ent = self.proteins.get(accession)?;
        let source = self.proteins.get(entrapment_source(accession)?)?;
        if peptide.is_empty() || ent.sequence.len() != source.sequence.len() {
            return None;
        }
        let mut matches = ent
            .sequence
            .windows(peptide.len())
            .enumerate()
            .filter(|(_, w)| *w == peptide);
        let (at, _) = matches.next()?;
        let paired = source.sequence.get(at..at + peptide.len())?;
        if !same_residue_composition(peptide, paired)
            || matches.any(|(at, _)| {
                let Some(candidate) = source.sequence.get(at..at + peptide.len()) else {
                    return true;
                };
                candidate != paired || !same_residue_composition(peptide, candidate)
            })
        {
            return None;
        }
        Some(paired)
    }
}

fn same_residue_composition(left: &[u8], right: &[u8]) -> bool {
    if left.len() != right.len() {
        return false;
    }
    let mut left = left
        .iter()
        .map(|&residue| if residue == b'I' { b'L' } else { residue })
        .collect::<Vec<_>>();
    let mut right = right
        .iter()
        .map(|&residue| if residue == b'I' { b'L' } else { residue })
        .collect::<Vec<_>>();
    left.sort_unstable();
    right.sort_unstable();
    left == right
}

/// The accession prefix with its separator, `ENTRAP_`.
pub fn entrapment_needle() -> String {
    format!("{ENTRAPMENT_PREFIX}_")
}

/// Whether `accession` is an entrapment protein's.
pub fn is_entrapment_accession(accession: &str) -> bool {
    accession
        .strip_prefix(ENTRAPMENT_PREFIX)
        .is_some_and(|rest| rest.starts_with('_'))
}

/// The accession of the target an entrapment protein was built from.
pub fn entrapment_source(accession: &str) -> Option<&str> {
    accession
        .strip_prefix(ENTRAPMENT_PREFIX)
        .and_then(|rest| rest.strip_prefix('_'))
}

/// Finds normalized query sequences in normalized target proteins, including overlapping matches.
pub fn target_substrings<'a>(
    targets: impl Iterator<Item = &'a [u8]>,
    queries: impl IntoIterator<Item = String>,
) -> Result<HashSet<String>, aho_corasick::BuildError> {
    let queries: Vec<String> = queries.into_iter().collect::<BTreeSet<_>>().into_iter().collect();
    if queries.is_empty() {
        return Ok(HashSet::new());
    }
    let matcher = AhoCorasick::new(&queries)?;
    let mut found = vec![false; queries.len()];
    let mut sequence = Vec::new();
    for target in targets {
        sequence.clear();
        sequence.extend(
            target
                .iter()
                .map(|&residue| if residue == b'I' { b'L' } else { residue }),
        );
        for hit in matcher.find_overlapping_iter(&sequence) {
            found[hit.pattern().as_usize()] = true;
        }
    }
    Ok(queries
        .into_iter()
        .zip(found)
        .filter_map(|(query, found)| found.then_some(query))
        .collect())
}

#[cfg(test)]
#[path = "tests/entrapment.rs"]
mod tests;
