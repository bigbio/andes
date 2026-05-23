//! Leaf types for SA-walk-based candidate enumeration. No logic; pure data.
//!
//! A `DistinctPeptide` represents a single unique residue sequence (no mods,
//! no flanking context) together with every `(protein, offset)` site where
//! that residue sequence occurs in the target+decoy database. This is the
//! shape produced by walking the suffix array with LCP-based deduplication
//! (`sa_walk::SaPeptideStream`): identical-residue suffixes get collapsed
//! into a single entry whose `positions` accumulate the per-protein
//! occurrences.
//!
//! Each `DistinctPeptide` keeps a single occurrence list keyed by residue
//! identity, with `positions: SmallVec<[Position; 4]>` — most peptides occur
//! in 1-3 proteins so the inline 4-slot smallvec avoids a heap allocation
//! on the common path.

use smallvec::SmallVec;

/// One occurrence of a peptide in the target+decoy database.
///
/// `protein_index` indexes into `SearchIndex.db.proteins` (target half is
/// `[0, target_count)`, decoy half is `[target_count, 2 * target_count)`).
/// `offset` is the start index of this peptide within the protein's residue
/// sequence (ASCII), NOT into the CompactFastaSequence body.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Position {
    pub protein_index: u32,
    pub offset: u32,
    pub is_decoy: bool,
    pub is_protein_n_term: bool,
    pub is_protein_c_term: bool,
}

/// A unique residue sequence and every place it occurs.
///
/// `residues` is the bare residue byte sequence (ASCII uppercase), with no
/// modifications and no flanking context — residue-only identity.
/// `nominal_mass` is the unmodified peptide nominal mass (residue masses +
/// `H2O`); variable-mod expansion happens in a later subtask layered on top
/// of this stream.
#[derive(Debug, Clone)]
pub struct DistinctPeptide {
    pub residues: Vec<u8>,
    pub nominal_mass: i32,
    pub positions: SmallVec<[Position; 4]>,
}

impl DistinctPeptide {
    pub fn new(residues: Vec<u8>, nominal_mass: i32) -> Self {
        Self {
            residues,
            nominal_mass,
            positions: SmallVec::new(),
        }
    }

    pub fn add_position(&mut self, pos: Position) {
        self.positions.push(pos);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_starts_with_no_positions() {
        let dp = DistinctPeptide::new(b"PEPTIDE".to_vec(), 799);
        assert_eq!(dp.residues, b"PEPTIDE");
        assert_eq!(dp.nominal_mass, 799);
        assert!(dp.positions.is_empty());
    }

    #[test]
    fn add_position_accumulates() {
        let mut dp = DistinctPeptide::new(b"PEPTIDE".to_vec(), 799);
        dp.add_position(Position {
            protein_index: 0,
            offset: 5,
            is_decoy: false,
            is_protein_n_term: false,
            is_protein_c_term: false,
        });
        dp.add_position(Position {
            protein_index: 3,
            offset: 12,
            is_decoy: true,
            is_protein_n_term: false,
            is_protein_c_term: false,
        });
        assert_eq!(dp.positions.len(), 2);
        assert_eq!(dp.positions[0].protein_index, 0);
        assert_eq!(dp.positions[1].is_decoy, true);
    }
}
