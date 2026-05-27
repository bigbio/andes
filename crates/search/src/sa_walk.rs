//! Suffix-array walk that produces `DistinctPeptide`s with LCP-based dedup.
//!
//! Walks `(indices[i], nlcps[i])` in lockstep and, for each peptide length L
//! in `[min, max]`, uses the LCP to decide whether the current suffix shares
//! the same residues (and possibly the same flanks) as the previous suffix:
//!
//! - `lcp >= L + 2`: residues + N-term flank + C-term flank are all shared
//!   with the previous suffix. The previous match's position list gets
//!   another `(protein, offset)` entry; no new distinct peptide is emitted.
//! - `lcp == L + 1`: residues + N-term flank are shared, but the C-term
//!   flank differs. The enzyme decides whether the new C-term flank still
//!   produces a cleavable peptide; if so, append to the previous match;
//!   otherwise start a fresh distinct peptide.
//! - `lcp < L + 1`: residues differ at or before position L. The previous
//!   match (if any) is emitted as a completed `DistinctPeptide`, and a new
//!   match is started at this cursor.
//!
//! This file deliberately implements ONLY the LCP-dedup walk: variable-mod
//! expansion, N-term Met cleavage, and the mass-tolerance filter all live
//! in later layers that consume the stream this iterator produces.
//!
//! ## Residue encoding note
//!
//! `compact.sequence` stores alphabet-indexed bytes (TERMINATOR=0,
//! INVALID=1, 'A'=2, ..., 'Z'=27). The bytes we emit on
//! `DistinctPeptide.residues` are ASCII uppercase residues (decoded via
//! `byte_to_residue`), so downstream consumers can treat them as ordinary
//! AA bytes.
//!
//! ## Simplification
//!
//! The `lcp == L + 1` enzyme-decision branch is currently treated the same
//! as the `lcp < L + 1` "new peptide" branch — i.e., we always start a
//! fresh DistinctPeptide. This costs a small amount of extra emission (the
//! same residue sequence may appear as two adjacent DistinctPeptides
//! differing in C-term flank) but is conservative and never silently
//! merges peptides the enzyme would consider distinct. Porting the full
//! enzyme branch is a follow-up.
//!
//! ## N-terminal Met-cleavage merge
//!
//! For each protein whose first residue is `M`, we run a separate enumeration
//! pass over `sequence[1..]` (the "initial-Met loss" virtual sequence) and
//! emit any peptides that pass the enzyme/length filters with
//! `is_protein_n_term = true` (the post-Met residue is the biological
//! N-terminus). These Met-cleaved variants are always emitted as
//! **separate** `DistinctPeptide`s from the main SA walk: dedup key is
//! `(residues, is_protein_n_term)`. The Met-cleaved pass dedupes among
//! itself by residue bytes (all entries share `is_protein_n_term = true`),
//! so two M-prefixed proteins yielding the same Met-cleaved residue
//! sequence aggregate into one `DistinctPeptide` with two positions —
//! while the same residue sequence appearing elsewhere (non-N-terminal,
//! or non-Met-prefixed N-term) remains a distinct entry from the main
//! pass. See `tests/sa_walk_met_cleavage.rs`.

use std::collections::HashMap;

use model::amino_acid::AminoAcid;
use model::compact_fasta::{byte_to_residue, INVALID_CHAR_CODE, TERMINATOR};
use model::enzyme::Enzyme;
use model::mass::{nominal_from, H2O};

use crate::distinct_peptide::{DistinctPeptide, Position};
use crate::search_index::SearchIndex;
use crate::search_params::SearchParams;

/// Streaming SA-walk iterator over `idx`. Emits one `DistinctPeptide` per
/// unique residue sequence (per peptide length) seen during the walk, with
/// every `(protein, offset)` position accumulated via LCP dedup.
///
/// Stateful: each `next()` call advances the SA cursor until at least one
/// completed `DistinctPeptide` is ready (or the walk ends). Emission order
/// is determined by SA order — same as Java.
pub struct SaPeptideStream<'a> {
    idx: &'a SearchIndex,
    params: &'a SearchParams,
    cursor: usize,
    /// `prev_match[length]` holds the in-progress DistinctPeptide for that
    /// length; `None` if the most recent suffix at that length was invalid
    /// (e.g., contained TERMINATOR) or no match has started yet.
    prev_match: Vec<Option<DistinctPeptide>>,
    /// Completed peptides ready to yield from the next `next()` call.
    pending: Vec<DistinctPeptide>,
    min_length: usize,
    max_length: usize,
    /// Cached per-protein decoy classification (indexed by protein_index).
    /// Avoids a string-prefix check on every emission.
    is_decoy: Vec<bool>,
    /// Set once the main SA-walk is exhausted and the Met-cleavage
    /// finalization pass has been queued into `pending`. Prevents double
    /// emission across repeated `next()` calls after the iterator drains.
    met_cleavage_emitted: bool,
}

impl<'a> SaPeptideStream<'a> {
    pub fn new(idx: &'a SearchIndex, params: &'a SearchParams, decoy_prefix: &'a str) -> Self {
        let min_length = params.min_length as usize;
        let max_length = params.max_length as usize;
        let is_decoy: Vec<bool> = idx
            .db
            .proteins
            .iter()
            .map(|p| p.accession.starts_with(decoy_prefix))
            .collect();
        Self {
            idx,
            params,
            cursor: 0,
            // Indexed 0..=max_length; prev_match[0] unused. +1 slot for ergonomic indexing.
            prev_match: (0..=max_length + 1).map(|_| None).collect(),
            pending: Vec::new(),
            min_length,
            max_length,
            is_decoy,
            met_cleavage_emitted: false,
        }
    }

    /// Resolve the cumulative `(protein_index, offset_in_protein,
    /// is_protein_n_term, is_protein_c_term)` for a suffix starting at
    /// CompactFastaSequence body position `index` and spanning `length`
    /// alphabet-encoded residue bytes. Returns `None` when `index` falls
    /// before the first protein (i.e., on the leading TERMINATOR byte) or
    /// when the span straddles a protein boundary.
    fn make_position(&self, index: usize, length: usize) -> Option<Position> {
        let p_idx = self.idx.compact.protein_index_at(index as u64)?;
        let ann = self.idx.compact.annotations.get(p_idx)?;
        let protein_start = ann.start as usize;
        let offset = index.checked_sub(protein_start)?;
        // The protein's residues are stored from `protein_start` up to (but
        // not including) the next TERMINATOR byte. If the span extends to
        // or past that TERMINATOR, this is not a valid in-protein peptide.
        let protein = self.idx.db.proteins.get(p_idx)?;
        if offset + length > protein.sequence.len() {
            return None;
        }
        let is_protein_n_term = offset == 0;
        let is_protein_c_term = offset + length == protein.sequence.len();
        Some(Position {
            protein_index: p_idx as u32,
            offset: offset as u32,
            is_decoy: self.is_decoy.get(p_idx).copied().unwrap_or(false),
            is_protein_n_term,
            is_protein_c_term,
        })
    }

    /// Build a fresh `DistinctPeptide` at the given SA index for the given
    /// length, applying residue validity + enzyme cleavage checks. Returns
    /// `None` when the peptide is rejected.
    fn build_distinct_peptide(&self, index: usize, length: usize) -> Option<DistinctPeptide> {
        let seq = &self.idx.compact.sequence;
        // Bounds + range guard.
        if index + length > seq.len() {
            return None;
        }
        // Decode the alphabet-encoded residues to ASCII; reject if any byte
        // is TERMINATOR/INVALID or maps outside the 20 standard AAs.
        let mut ascii = Vec::with_capacity(length);
        for &b in &seq[index..index + length] {
            if b == TERMINATOR || b == INVALID_CHAR_CODE {
                return None;
            }
            let aa = byte_to_residue(b);
            AminoAcid::standard(aa)?;
            ascii.push(aa);
        }
        // Position resolution doubles as a protein-boundary check: if the
        // span straddles two proteins, `make_position` returns None.
        let position = self.make_position(index, length)?;

        // Enzyme NTT (num tolerable termini) check. The pre flank is the
        // body byte before `index`; post is the body byte at `index+length`.
        // For protein-terminal positions we treat the flank as cleavable.
        let pre_byte = if index == 0 { TERMINATOR } else { seq[index - 1] };
        let post_byte = seq[index + length]; // safe: index+length <= seq.len()-? — body always ends in TERM, so this is valid for any legal peptide that fits within a protein.
        let pre_ascii = if pre_byte == TERMINATOR {
            None
        } else {
            Some(byte_to_residue(pre_byte))
        };
        let post_ascii = if post_byte == TERMINATOR {
            None
        } else {
            Some(byte_to_residue(post_byte))
        };

        if !self.passes_ntt(&ascii, pre_ascii, post_ascii) {
            return None;
        }

        let nominal_mass = compute_nominal_mass(&ascii);
        let mut dp = DistinctPeptide::new(ascii, nominal_mass);
        dp.add_position(position);
        Some(dp)
    }

    /// Number-of-tolerable-termini check:
    /// - NTT=2 (strict): both ends must be enzyme-cleavable.
    /// - NTT=1 (semi):   at least one end must be cleavable.
    /// - NTT=0 (none):   no constraint.
    ///
    /// For Trypsin-like C-term cutters, "N-term cleavable" means the
    /// preceding residue is K/R (or protein N-term); "C-term cleavable"
    /// means the last residue of the peptide is K/R (or protein C-term).
    fn passes_ntt(&self, residues: &[u8], pre: Option<u8>, post: Option<u8>) -> bool {
        let ntt = self.params.num_tolerable_termini;
        if ntt == 0 {
            return true;
        }
        let enzyme = self.params.enzyme;
        if matches!(enzyme, Enzyme::NonSpecific) {
            return true;
        }
        let n_ok = match pre {
            None => true, // protein N-term: trivially cleavable
            Some(p) => enzyme.is_cleavable_after(p) || enzyme.is_cleavable_before(residues[0]),
        };
        let c_ok = match post {
            None => true, // protein C-term
            Some(post_r) => {
                let last = *residues.last().unwrap();
                enzyme.is_cleavable_after(last) || enzyme.is_cleavable_before(post_r)
            }
        };
        match ntt {
            2 => n_ok && c_ok,
            _ => n_ok || c_ok, // ntt == 1 (or any other non-zero/non-2 value, treated as 1)
        }
    }

    /// Displace the in-progress `prev_match[length]` (push it to pending)
    /// and install a fresh DistinctPeptide for the current cursor at that
    /// length. If the cursor's peptide is invalid, `prev_match[length]` is
    /// left as `None`.
    fn start_new(&mut self, index: usize, length: usize) {
        if let Some(prev) = self.prev_match[length].take() {
            self.pending.push(prev);
        }
        if let Some(fresh) = self.build_distinct_peptide(index, length) {
            self.prev_match[length] = Some(fresh);
        }
    }

    /// Append the cursor's `(protein, offset)` position to
    /// `prev_match[length]` if a match is in progress. If
    /// `prev_match[length]` is `None` (no in-progress match), do nothing —
    /// this can happen when an earlier cursor at the same length was
    /// invalid (e.g., the suffix contained a TERMINATOR). The shared-LCP
    /// guarantee from the SA still holds (suffixes share their first L
    /// characters), but if those characters include a TERMINATOR neither
    /// suffix can produce a valid peptide.
    fn append_position(&mut self, index: usize, length: usize) {
        // Resolve position first to release the immutable self-borrow before
        // taking the mutable borrow on prev_match.
        let pos = self.make_position(index, length);
        if let (Some(prev), Some(p)) = (self.prev_match[length].as_mut(), pos) {
            prev.add_position(p);
        }
    }

    /// Enumerate Met-cleaved peptide variants and append them to
    /// `self.pending`. For each M-prefixed protein, treat `sequence[1..]`
    /// as a virtual protein (post-initial-Met cleavage), enumerate spans
    /// that pass the same residue + enzyme + length filters used by the
    /// main SA pass, and emit them with `is_protein_n_term = true`. The
    /// pre-flank for spans starting at offset 1 of the original protein
    /// is the cleaved `M` itself, so the NTT check uses `pre = Some(b'M')`.
    ///
    /// Multiple M-prefixed proteins producing the same Met-cleaved residue
    /// sequence are aggregated into a single `DistinctPeptide` (positions
    /// vector lists each `(protein, offset=1+..)` site). This matches the
    /// dedup contract for the main SA pass — residue-only identity within
    /// the Met-cleaved sub-pass — while keeping Met-cleaved peptides as
    /// separate `DistinctPeptide`s from non-Met-cleaved peptides with the
    /// same residues (the `is_protein_n_term` axis differs).
    fn enumerate_met_cleaved(&mut self) {
        // Aggregate by residue bytes. All entries here share is_protein_n_term=true.
        let mut by_residues: HashMap<Vec<u8>, DistinctPeptide> = HashMap::new();

        for (p_idx, protein) in self.idx.db.proteins.iter().enumerate() {
            let seq = &protein.sequence;
            if seq.first() != Some(&b'M') || seq.len() <= 1 {
                continue;
            }
            // Met-cleavage's unique contribution: peptides starting at
            // offset 1 of the original protein (the post-Met biological
            // N-terminus). Spans with start > 1 are already enumerated by
            // the main SA walk with is_protein_n_term=false at their
            // native location, so we don't repeat them here.
            let seq_len = seq.len();
            let min_l = self.min_length;
            let max_l = self.max_length;
            if seq_len < 1 + min_l {
                continue;
            }
            let start = 1usize;
            let max_end = seq_len.min(start + max_l);
            for end in (start + min_l)..=max_end {
                let span = &seq[start..end];
                // Residue validity: standard AAs only.
                let mut residues = Vec::with_capacity(span.len());
                let mut ok = true;
                for &b in span {
                    if AminoAcid::standard(b).is_none() {
                        ok = false;
                        break;
                    }
                    residues.push(b);
                }
                if !ok {
                    continue;
                }
                // NTT pre-flank for offset=1 is the cleaved M itself.
                let pre = Some(b'M');
                let post = if end == seq_len { None } else { Some(seq[end]) };
                if !self.passes_ntt(&residues, pre, post) {
                    continue;
                }
                let is_protein_c_term = end == seq_len;
                let position = Position {
                    protein_index: p_idx as u32,
                    offset: start as u32,
                    is_decoy: self.is_decoy.get(p_idx).copied().unwrap_or(false),
                    is_protein_n_term: true, // post-Met biological N-terminus
                    is_protein_c_term,
                };
                let nominal_mass = compute_nominal_mass(&residues);
                let entry = by_residues
                    .entry(residues.clone())
                    .or_insert_with(|| DistinctPeptide::new(residues, nominal_mass));
                entry.add_position(position);
            }
        }

        // Drain into pending. Order is unspecified but deterministic-ish
        // (HashMap iteration); downstream consumers must not rely on order.
        self.pending.extend(by_residues.into_values());
    }
}

impl<'a> Iterator for SaPeptideStream<'a> {
    type Item = DistinctPeptide;

    fn next(&mut self) -> Option<DistinctPeptide> {
        // Drain pending queue first.
        if let Some(dp) = self.pending.pop() {
            return Some(dp);
        }
        let sa_size = self.idx.sa.indices.len();
        while self.cursor < sa_size {
            let index = self.idx.sa.indices[self.cursor] as usize;
            let lcp = if self.cursor == 0 {
                0
            } else {
                self.idx.sa.nlcps[self.cursor] as i64
            };

            for length in self.min_length..=self.max_length {
                let l = length as i64;
                if lcp >= l + 2 {
                    // Shared peptide + flanks: append position to prev_match[length].
                    self.append_position(index, length);
                } else if lcp == l + 1 {
                    // Shared peptide, possibly different C-term flank.
                    // SIMPLIFICATION (see module docs): treat as a new
                    // peptide. Conservative — never silently merges across
                    // a C-term flank change.
                    self.start_new(index, length);
                } else {
                    // Residues differ at or before this length: start a
                    // new distinct peptide. Pre-existing prev_match[length]
                    // is emitted to pending.
                    self.start_new(index, length);
                }
            }

            self.cursor += 1;
            if let Some(dp) = self.pending.pop() {
                return Some(dp);
            }
        }
        // End of walk: flush remaining in-progress matches.
        for length in self.min_length..=self.max_length {
            if let Some(dp) = self.prev_match[length].take() {
                self.pending.push(dp);
            }
        }
        // Met-cleavage finalization: enumerate Met-cleaved peptides for
        // every M-prefixed protein and queue them as separate
        // DistinctPeptides distinguished by (residues, is_protein_n_term=true).
        if !self.met_cleavage_emitted {
            self.met_cleavage_emitted = true;
            self.enumerate_met_cleaved();
        }
        self.pending.pop()
    }
}

/// Compute the unmodified peptide nominal mass from an ASCII residue
/// sequence. Sum residue masses (no mods at this layer) + H2O, then floor
/// via Java's `Constants.INTEGER_MASS_SCALER` conversion.
fn compute_nominal_mass(ascii_residues: &[u8]) -> i32 {
    let residue_sum: f64 = ascii_residues
        .iter()
        .filter_map(|&r| AminoAcid::standard(r).map(|aa| aa.mass))
        .sum();
    nominal_from(residue_sum + H2O)
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::aa_set::AminoAcidSetBuilder;
    use model::protein::ProteinDb;

    fn aa_set() -> model::aa_set::AminoAcidSet {
        AminoAcidSetBuilder::new_standard().build().unwrap()
    }

    #[test]
    fn empty_db_yields_no_peptides() {
        let target = ProteinDb { proteins: vec![] };
        let idx = SearchIndex::from_target_db(&target, "XXX");
        let mut params = SearchParams::default_tryptic(aa_set());
        params.min_length = 6;
        params.max_length = 10;
        let peptides: Vec<_> = SaPeptideStream::new(&idx, &params, "XXX").collect();
        assert!(peptides.is_empty());
    }

    #[test]
    fn nominal_mass_includes_h2o() {
        // GA: G=57, A=71, +H2O ≈ 18 → 146
        let mass = compute_nominal_mass(b"GA");
        assert_eq!(mass, 146);
    }
}
