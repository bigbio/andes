//! Positional-isomer identity of a PSM.
//!
//! Two PSMs on one spectrum are positional isomers when they carry the same
//! residues and the same multiset of modification mass deltas, with the
//! modifications on different sites. Site alternatives of one peptide score
//! almost identically, so the per-spectrum rows keep one of them
//! ([`collapse_isomers`]) and the per-spectrum competition features
//! (DeltaRawScore, ListwiseScoreGap, CandidateRankEntropy) are taken over the
//! same collapsed pool, so a kept row never counts its own dropped isomer as
//! its runner-up.

use crate::candidate_gen::Candidate;
use crate::psm::PsmMatch;
use model::Peptide;

/// A modification mass delta in units of 1e-5 Da, the precision of the PIN
/// `Peptide` column.
#[inline]
fn mod_delta_units(mass_delta: f64) -> i64 {
    (mass_delta * 100_000.0).round() as i64
}

/// Identity of a PSM up to the placement of its modifications: the precursor it
/// was matched at, its charge, its target/decoy side, its residues, and the
/// sorted multiset of its modification mass deltas (terminal mods included,
/// since they are carried on the terminal residue). Two PSMs with equal keys are
/// positional isomers — the same peptide with the same modifications on
/// different sites.
///
/// A target and a decoy never share a key, so a decoy that repeats a target's
/// residues keeps its own row and its `-1` label.
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct IsomerKey {
    precursor_mz_bits: u64,
    charge: u8,
    is_decoy: bool,
    residues: Vec<u8>,
    mod_deltas: Vec<i64>,
}

impl IsomerKey {
    /// The key of `psm` (whose primary candidate is `cand`) on a spectrum with
    /// precursor m/z `spectrum_precursor_mz`. A chimeric secondary is keyed at
    /// its own precursor (`precursor_mz_override`).
    pub fn new(spectrum_precursor_mz: f64, psm: &PsmMatch, cand: &Candidate) -> Self {
        let residues = cand.peptide.residues.iter().map(|aa| aa.residue).collect();
        let mut mod_deltas: Vec<i64> = cand
            .peptide
            .residues
            .iter()
            .filter_map(|aa| aa.mod_.map(|m| mod_delta_units(m.mass_delta)))
            .collect();
        mod_deltas.sort_unstable();
        Self {
            precursor_mz_bits: psm
                .precursor_mz_override
                .unwrap_or(spectrum_precursor_mz)
                .to_bits(),
            charge: psm.charge_used,
            is_decoy: cand.is_decoy,
            residues,
            mod_deltas,
        }
    }
}

/// Allocation-free 64-bit hash of a peptide's residues and its multiset of
/// modification mass deltas (quantised as in [`IsomerKey`]); positional
/// isomers hash equal. The deltas are combined with a commutative sum of
/// mixed values, so the hash does not depend on which residues carry them.
/// Used on every scored candidate, where [`IsomerKey`]'s allocations would
/// cost too much; a collision merges two unrelated peptides, with negligible
/// probability.
#[inline]
pub fn isomer_class_hash(peptide: &Peptide) -> u64 {
    use std::hash::Hasher;
    let mut h = rustc_hash::FxHasher::default();
    let mut mods: u64 = 0;
    let mut n_mods: u64 = 0;
    for aa in &peptide.residues {
        h.write_u8(aa.residue);
        if let Some(m) = aa.mod_.as_ref() {
            mods = mods.wrapping_add(splitmix64(mod_delta_units(m.mass_delta) as u64));
            n_mods += 1;
        }
    }
    h.write_u64(mods);
    h.write_u64(n_mods);
    h.finish()
}

/// SplitMix64 finaliser: a bijective mix that spreads each delta over all 64
/// bits before the commutative sum.
#[inline]
fn splitmix64(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// The PSMs that represent each positional-isomer class of one spectrum, in
/// the order given: of the PSMs sharing an [`IsomerKey`], the first one is
/// kept and the rest are dropped. `psms_best_first` must be in rank order
/// (best first), so the kept PSM is the best-ranked isomer.
///
/// Returns, for each kept PSM, its index in `psms_best_first` and the indices
/// of the PSMs collapsed into it (empty when it has no isomer).
pub fn collapse_isomers<'a, I>(
    spectrum_precursor_mz: f64,
    psms_best_first: I,
    candidates: &[Candidate],
) -> Vec<(usize, Vec<usize>)>
where
    I: IntoIterator<Item = &'a PsmMatch>,
{
    let mut keys: Vec<IsomerKey> = Vec::new();
    let mut kept: Vec<(usize, Vec<usize>)> = Vec::new();
    for (i, psm) in psms_best_first.into_iter().enumerate() {
        let cand = &candidates[psm.primary_candidate_idx() as usize];
        let key = IsomerKey::new(spectrum_precursor_mz, psm, cand);
        match keys.iter().position(|k| *k == key) {
            Some(slot) => kept[slot].1.push(i),
            None => {
                keys.push(key);
                kept.push((i, Vec::new()));
            }
        }
    }
    kept
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::{AminoAcid, ModLocation, Modification, Peptide, ResidueSpec};

    fn aa(c: u8) -> AminoAcid {
        AminoAcid::standard(c).unwrap()
    }

    fn modded(c: u8, delta: f64) -> AminoAcid {
        let mut a = aa(c);
        a.mod_ = Some(model::modification::leak_mod(Modification {
            name: "m".into(),
            mass_delta: delta,
            residue: ResidueSpec::Specific(c),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
        }));
        a
    }

    fn pep(residues: Vec<AminoAcid>) -> Peptide {
        Peptide::new(residues, b'K', b'A')
    }

    #[test]
    fn class_hash_ignores_mod_site_but_not_mod_mass_or_residues() {
        let phospho = 79.96633;
        let a = pep(vec![modded(b'S', phospho), aa(b'T'), aa(b'K')]);
        let b = pep(vec![aa(b'S'), modded(b'T', phospho), aa(b'K')]);
        let unmod = pep(vec![aa(b'S'), aa(b'T'), aa(b'K')]);
        let other = pep(vec![aa(b'T'), modded(b'S', phospho), aa(b'K')]);
        let ox = pep(vec![modded(b'S', 15.99491), aa(b'T'), aa(b'K')]);
        assert_eq!(isomer_class_hash(&a), isomer_class_hash(&b));
        assert_ne!(isomer_class_hash(&a), isomer_class_hash(&unmod));
        assert_ne!(isomer_class_hash(&a), isomer_class_hash(&other));
        assert_ne!(isomer_class_hash(&a), isomer_class_hash(&ox));
    }
}
