//! Peptide-level dedup of PSMs, where I/L twins are one PSM, and the peptidoform keys it uses;
//! the I/L-folded bare-sequence sets that relabel decoys colliding with a target (I and L are
//! isobaric, so such a decoy is the target's precursor and fragment ladder); the merge of
//! candidate indices. Part of the match engine.

use super::*;

pub(crate) fn dedup_pepseq_score(psms: Vec<PsmMatch>, candidates: &[Candidate]) -> Vec<PsmMatch> {
    use std::collections::btree_map::Entry;
    use std::collections::BTreeMap;

    let mut pep_key_cache: FxHashMap<u32, Arc<PepDedupKey>> = FxHashMap::default();
    let mut groups: BTreeMap<DedupMapKey, PsmMatch> = BTreeMap::new();

    let mut psms = psms;
    psms.sort_by(|a, b| b.cmp(a));
    for psm in psms {
        let primary = psm.primary_candidate_idx();
        let pep_key = pep_key_cache
            .entry(primary)
            .or_insert_with(|| {
                Arc::new(PepDedupKey::from_peptide(
                    &candidates[primary as usize].peptide,
                ))
            })
            .clone();
        let key = DedupMapKey {
            pep: pep_key,
            score: psm.rank_score.round() as i32,
            is_decoy: candidates[primary as usize].is_decoy,
        };

        match groups.entry(key) {
            Entry::Vacant(slot) => {
                slot.insert(psm);
            }
            Entry::Occupied(mut slot) => {
                let existing = slot.get_mut();
                merge_unique_candidate_idxs(&mut existing.candidate_idxs, &psm.candidate_idxs);
                if survives(&psm, existing, candidates) {
                    let merged_idxs = std::mem::take(&mut existing.candidate_idxs);
                    let mut survivor = psm;
                    merge_unique_candidate_idxs(&mut survivor.candidate_idxs, &merged_idxs);
                    *existing = survivor;
                }
            }
        }
    }

    let mut out = Vec::with_capacity(groups.len());
    out.extend(groups.into_values());
    out
}

/// Whether `psm` replaces `existing` as the reported row of their dedup group: a higher rank
/// score wins; at an equal rank score (I/L twins score the same), the lower database sequence
/// wins, so the twins of a precursor report one sequence in every spectrum, the same rule for
/// targets and decoys and independent of the protein order of the database.
fn survives(psm: &PsmMatch, existing: &PsmMatch, candidates: &[Candidate]) -> bool {
    match psm.rank_score.total_cmp(&existing.rank_score) {
        std::cmp::Ordering::Greater => true,
        std::cmp::Ordering::Less => false,
        std::cmp::Ordering::Equal => {
            let residues = |p: &PsmMatch| {
                candidates[p.primary_candidate_idx() as usize]
                    .peptide
                    .residues
                    .iter()
                    .map(|aa| aa.residue)
            };
            residues(psm).lt(residues(existing))
        }
    }
}

#[inline]
pub(crate) fn peptidoform_key(peptide: &model::Peptide) -> u64 {
    peptidoform_hash(peptide, |r| r)
}

/// The key of a spectrum's competing peptidoforms (delta, reference and runner-up scores): I/L
/// twins share it, so a twin is never the runner-up of its own peptide.
#[inline]
pub(crate) fn competitor_key(peptide: &model::Peptide) -> u64 {
    peptidoform_hash(peptide, fold_isoleucine)
}

#[inline]
fn peptidoform_hash(peptide: &model::Peptide, residue: impl Fn(u8) -> u8) -> u64 {
    use std::hash::Hasher;
    let mut h = rustc_hash::FxHasher::default();
    for aa in &peptide.residues {
        h.write_u8(residue(aa.residue));
        let units = aa
            .mod_
            .as_ref()
            .map(|m| (m.mass_delta * 100_000.0).round() as i64)
            .unwrap_or(0);
        h.write_i64(units);
    }
    h.finish()
}

pub(super) use model::isomers::fold_isoleucine;

/// The residues of a candidate with I folded to L, the key of the target peptide set.
pub(super) fn bare_residues(cand: &Candidate) -> Box<[u8]> {
    cand.peptide
        .residues
        .iter()
        .map(|aa| fold_isoleucine(aa.residue))
        .collect()
}

/// A set of I/L-folded peptide sequences, split into shards by hash so that the shards are
/// built in parallel. It only answers membership.
#[derive(Default)]
pub(crate) struct FoldedSet {
    shards: Vec<FxHashSet<Box<[u8]>>>,
}

const FOLDED_SET_SHARD_BITS: u32 = 6;

/// The shard of a sequence. SipHash rather than FxHash: the shard bits must not correlate with
/// the hash a shard's own set probes by.
fn folded_shard(seq: &[u8]) -> usize {
    use std::hash::Hasher;
    let mut h = std::hash::DefaultHasher::new();
    h.write(seq);
    (h.finish() >> (64 - FOLDED_SET_SHARD_BITS)) as usize
}

impl FoldedSet {
    pub(crate) fn is_empty(&self) -> bool {
        self.shards.iter().all(FxHashSet::is_empty)
    }

    pub(crate) fn contains(&self, seq: &[u8]) -> bool {
        self.shards
            .get(folded_shard(seq))
            .is_some_and(|shard| shard.contains(seq))
    }
}

pub(super) fn folded_set(spans: Vec<&[u8]>) -> FoldedSet {
    let n_shards = 1usize << FOLDED_SET_SHARD_BITS;
    let buckets: Vec<Vec<Box<[u8]>>> = spans
        .par_iter()
        .fold(
            || vec![Vec::new(); n_shards],
            |mut buckets, s| {
                let folded: Box<[u8]> = s.iter().map(|&r| fold_isoleucine(r)).collect();
                buckets[folded_shard(&folded)].push(folded);
                buckets
            },
        )
        .reduce(
            || vec![Vec::new(); n_shards],
            |mut a, b| {
                for (into, from) in a.iter_mut().zip(b) {
                    into.extend(from);
                }
                a
            },
        );
    let shards = buckets
        .into_par_iter()
        .map(|bucket| {
            let mut set = FxHashSet::with_capacity_and_hasher(bucket.len(), Default::default());
            set.extend(bucket);
            set
        })
        .collect();
    FoldedSet { shards }
}

pub(super) fn relabel_collision_decoys_with(candidates: &mut [Candidate], target_seqs: &FoldedSet) {
    if target_seqs.is_empty() {
        return;
    }
    for cand in candidates.iter_mut().filter(|c| c.is_decoy) {
        if target_seqs.contains(&bare_residues(cand)) {
            cand.is_decoy = false;
        }
    }
}

pub(super) fn build_target_bare_seqs(mi: &MmapCandidateIndex, db: &SearchIndex) -> FoldedSet {
    use crate::candidates::records::flags as rec_flags;
    let mut spans: Vec<&[u8]> = (0..mi.len())
        .into_par_iter()
        .filter_map(|i| {
            let rec = mi.record_at(i);
            if rec.flags & rec_flags::IS_DECOY != 0 {
                return None;
            }
            let seq = &db.db.proteins.get(rec.protein_index as usize)?.sequence;
            let start = rec.start_offset as usize;
            seq.get(start..start + rec.length as usize)
        })
        .collect();
    spans.par_sort_unstable();
    spans.dedup();
    folded_set(spans)
}

pub(super) fn build_target_bare_seqs_where(
    mi: &MmapCandidateIndex,
    db: &SearchIndex,
    params: &SearchParams,
    keep: &CandidateFilter,
) -> FoldedSet {
    use crate::candidates::records::flags as rec_flags;
    let mut spans: Vec<&[u8]> = (0..mi.len())
        .into_par_iter()
        .filter_map(|i| {
            let rec = mi.record_at(i);
            if rec.flags & rec_flags::IS_DECOY != 0 {
                return None;
            }
            let seq = &db.db.proteins.get(rec.protein_index as usize)?.sequence;
            let start = rec.start_offset as usize;
            let span = seq.get(start..start + rec.length as usize)?;
            expand_base_record(db, params, &rec)
                .iter()
                .any(keep)
                .then_some(span)
        })
        .collect();
    spans.par_sort_unstable();
    spans.dedup();
    folded_set(spans)
}

/// The identity of a peptidoform in the dedup of one spectrum's PSMs: its residues with I read as
/// L and the mass of each residue's modification, so I/L twins are one PSM.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(super) struct PepDedupKey {
    residues: Vec<u8>,
    mod_units: Vec<i32>,
}

impl PepDedupKey {
    pub(super) fn from_peptide(peptide: &model::Peptide) -> Self {
        let len = peptide.residues.len();
        let mut residues = Vec::with_capacity(len);
        let mut mod_units = Vec::with_capacity(len);
        for aa in &peptide.residues {
            residues.push(fold_isoleucine(aa.residue));
            mod_units.push(
                aa.mod_
                    .as_ref()
                    .map(|m| (m.mass_delta * 100_000.0).round() as i32)
                    .unwrap_or(0),
            );
        }
        Self {
            residues,
            mod_units,
        }
    }
}

#[derive(Clone, PartialEq, Eq)]
pub(super) struct DedupMapKey {
    pep: Arc<PepDedupKey>,
    score: i32,
    is_decoy: bool,
}

impl PartialOrd for DedupMapKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for DedupMapKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.pep
            .residues
            .cmp(&other.pep.residues)
            .then_with(|| self.pep.mod_units.cmp(&other.pep.mod_units))
            .then(self.score.cmp(&other.score))
            .then(self.is_decoy.cmp(&other.is_decoy))
    }
}

pub(super) fn merge_unique_candidate_idxs(into: &mut crate::psm::CandidateIdxs, from: &[u32]) {
    if from.is_empty() {
        return;
    }
    let mut seen: FxHashSet<u32> = into.iter().copied().collect();
    for &idx in from {
        if seen.insert(idx) {
            into.push(idx);
        }
    }
}
