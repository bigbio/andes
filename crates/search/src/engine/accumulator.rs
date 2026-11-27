//! The candidates a search materialises from the out-of-core index, deduplicated and kept in a
//! compact form: one byte per residue (an index into a table of the distinct residues seen) and
//! a hash-chained index on the candidate's identity. Reading a candidate back rebuilds it equal,
//! bit for bit, to the one interned.

use std::hash::{Hash, Hasher};

use model::amino_acid::AminoAcid;
use model::mass::H2O;
use model::peptide::Peptide;
use rustc_hash::{FxHashMap, FxHasher};

use crate::candidates::Candidate;

const NONE: u32 = u32::MAX;
const DECOY: u8 = 1;
const PROTEIN_N_TERM: u8 = 2;
const PROTEIN_C_TERM: u8 = 4;
/// The candidate did not fit the compact form and is kept whole in `whole`.
const WHOLE: u8 = 8;
/// Scale of the modification mass in the identity key (0.01 mDa units).
const MOD_KEY_SCALE: f64 = 100_000.0;

/// One materialised candidate: its identity fields, its flanks and charge, and where its residue
/// codes start in `codes` (they end where the next candidate's start).
struct Meta {
    hash: u64,
    codes_at: u64,
    protein_index: u32,
    start_offset: u32,
    pre: u8,
    post: u8,
    charge: Option<u8>,
    flags: u8,
}

/// Key of one residue in the table: residue letter, residue mass bits and modification address.
type ResidueKey = (u8, u64, usize);

#[derive(Default)]
pub(super) struct MmapAccumulator {
    metas: Vec<Meta>,
    codes: Vec<u8>,
    residues: Vec<AminoAcid>,
    residue_code: FxHashMap<ResidueKey, u8>,
    whole: FxHashMap<u32, Candidate>,
    heads: FxHashMap<u64, u32>,
    chain: Vec<u32>,
}

fn residue_key(aa: &AminoAcid) -> ResidueKey {
    (
        aa.residue,
        aa.mass.to_bits(),
        aa.mod_.map_or(0, |m| m as *const _ as usize),
    )
}

fn mod_unit(aa: &AminoAcid) -> i32 {
    aa.mod_
        .as_ref()
        .map(|m| (m.mass_delta * MOD_KEY_SCALE).round() as i32)
        .unwrap_or(0)
}

fn flags_of(c: &Candidate) -> u8 {
    (u8::from(c.is_decoy) * DECOY)
        | (u8::from(c.is_protein_n_term) * PROTEIN_N_TERM)
        | (u8::from(c.is_protein_c_term) * PROTEIN_C_TERM)
}

/// Hash of a candidate's identity: protein, start, residue letters with their modification
/// mass at `MOD_KEY_SCALE`, decoy and protein-terminus flags.
fn identity_hash<'r>(
    protein_index: u32,
    start_offset: u32,
    flags: u8,
    residues: impl Iterator<Item = &'r AminoAcid>,
) -> u64 {
    let mut h = FxHasher::default();
    protein_index.hash(&mut h);
    start_offset.hash(&mut h);
    (flags & (DECOY | PROTEIN_N_TERM | PROTEIN_C_TERM)).hash(&mut h);
    for aa in residues {
        aa.residue.hash(&mut h);
        mod_unit(aa).hash(&mut h);
    }
    h.finish()
}

/// The neutral mass [`Peptide::new`] computes from `residues`.
fn neutral_mass_of(residues: &[AminoAcid]) -> f64 {
    let residue_mass: f64 = residues
        .iter()
        .map(|aa| aa.mass + aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta))
        .sum();
    residue_mass + H2O
}

impl MmapAccumulator {
    /// Number of candidates materialised; the index the next new one gets.
    pub(super) fn len(&self) -> usize {
        self.metas.len()
    }

    fn codes_of(&self, i: usize) -> &[u8] {
        let at = self.metas[i].codes_at as usize;
        let end = self
            .metas
            .get(i + 1)
            .map_or(self.codes.len(), |m| m.codes_at as usize);
        &self.codes[at..end]
    }

    fn same_identity(&self, i: usize, c: &Candidate) -> bool {
        if let Some(w) = self.whole.get(&(i as u32)) {
            return w.protein_index as u32 == c.protein_index as u32
                && w.start_offset_in_protein as u32 == c.start_offset_in_protein as u32
                && flags_of(w) == flags_of(c)
                && w.peptide.residues.len() == c.peptide.residues.len()
                && w.peptide
                    .residues
                    .iter()
                    .zip(&c.peptide.residues)
                    .all(|(a, b)| a.residue == b.residue && mod_unit(a) == mod_unit(b));
        }
        let m = &self.metas[i];
        let codes = self.codes_of(i);
        m.protein_index == c.protein_index as u32
            && m.start_offset == c.start_offset_in_protein as u32
            && m.flags & !WHOLE == flags_of(c)
            && codes.len() == c.peptide.residues.len()
            && codes.iter().zip(&c.peptide.residues).all(|(&k, b)| {
                let a = &self.residues[k as usize];
                a.residue == b.residue && mod_unit(a) == mod_unit(b)
            })
    }

    /// Residue codes of `c`, adding new residues to the table; `None` when the table is full.
    fn encode(&mut self, c: &Candidate) -> Option<Vec<u8>> {
        let mut out = Vec::with_capacity(c.peptide.residues.len());
        for aa in &c.peptide.residues {
            let key = residue_key(aa);
            let code = match self.residue_code.get(&key) {
                Some(&k) => k,
                None => {
                    let k = u8::try_from(self.residues.len()).ok()?;
                    self.residues.push(aa.clone());
                    self.residue_code.insert(key, k);
                    k
                }
            };
            out.push(code);
        }
        Some(out)
    }

    /// The global index of `cand`: an earlier candidate with the same identity (protein, start,
    /// residues with their modification mass, decoy and terminus flags) keeps its index, a new
    /// one gets the next index.
    pub(super) fn intern(&mut self, cand: &Candidate) -> u32 {
        let flags = flags_of(cand);
        let hash = identity_hash(
            cand.protein_index as u32,
            cand.start_offset_in_protein as u32,
            flags,
            cand.peptide.residues.iter(),
        );
        let mut at = self.heads.get(&hash).copied().unwrap_or(NONE);
        while at != NONE {
            if self.same_identity(at as usize, cand) {
                return at;
            }
            at = self.chain[at as usize];
        }
        let idx = self.metas.len() as u32;
        let fits = u32::try_from(cand.protein_index).is_ok()
            && u32::try_from(cand.start_offset_in_protein).is_ok()
            && neutral_mass_of(&cand.peptide.residues).to_bits() == cand.peptide.mass().to_bits();
        let codes = if fits { self.encode(cand) } else { None };
        let codes_at = self.codes.len() as u64;
        let flags = match codes {
            Some(codes) => {
                self.codes.extend_from_slice(&codes);
                flags
            }
            None => {
                self.whole.insert(idx, cand.clone());
                flags | WHOLE
            }
        };
        self.metas.push(Meta {
            hash,
            codes_at,
            protein_index: cand.protein_index as u32,
            start_offset: cand.start_offset_in_protein as u32,
            pre: cand.peptide.pre,
            post: cand.peptide.post,
            charge: cand.peptide.charge,
            flags,
        });
        self.chain
            .push(self.heads.insert(hash, idx).unwrap_or(NONE));
        idx
    }

    /// The peptide of candidate `i`, equal to the one interned.
    pub(super) fn peptide(&self, i: usize) -> Option<Peptide> {
        let m = self.metas.get(i)?;
        if m.flags & WHOLE != 0 {
            return self.whole.get(&(i as u32)).map(|c| c.peptide.clone());
        }
        let residues: Vec<AminoAcid> = self
            .codes_of(i)
            .iter()
            .map(|&k| self.residues[k as usize].clone())
            .collect();
        let p = Peptide::new(residues, m.pre, m.post);
        Some(match m.charge {
            Some(z) => p.with_charge(z),
            None => p,
        })
    }

    pub(super) fn is_decoy(&self, i: usize) -> Option<bool> {
        self.metas.get(i).map(|m| m.flags & DECOY != 0)
    }

    /// Candidate `i`, equal to the one interned.
    pub(super) fn get(&self, i: usize) -> Option<Candidate> {
        let m = self.metas.get(i)?;
        if m.flags & WHOLE != 0 {
            return self.whole.get(&(i as u32)).cloned();
        }
        Some(Candidate {
            peptide: self.peptide(i)?,
            protein_index: m.protein_index as usize,
            start_offset_in_protein: m.start_offset as usize,
            is_decoy: m.flags & DECOY != 0,
            is_protein_n_term: m.flags & PROTEIN_N_TERM != 0,
            is_protein_c_term: m.flags & PROTEIN_C_TERM != 0,
        })
    }

    /// Every candidate in index order; the accumulator is left empty.
    pub(super) fn take_all(&mut self) -> Vec<Candidate> {
        use rayon::prelude::*;
        let out = (0..self.len())
            .into_par_iter()
            .filter_map(|i| self.get(i))
            .collect();
        *self = Self::default();
        out
    }

    /// Forgets the candidates from index `len` on and returns their memory.
    pub(super) fn truncate(&mut self, len: usize) {
        if self.metas.len() <= len {
            return;
        }
        self.codes.truncate(self.metas[len].codes_at as usize);
        self.codes.shrink_to_fit();
        self.metas.truncate(len);
        self.metas.shrink_to_fit();
        self.chain.truncate(len);
        self.chain.shrink_to_fit();
        self.whole.retain(|&i, _| (i as usize) < len);
        self.whole.shrink_to_fit();
        self.heads.clear();
        for (i, m) in self.metas.iter().enumerate() {
            self.heads.insert(m.hash, i as u32);
        }
        self.heads.shrink_to_fit();
    }
}
