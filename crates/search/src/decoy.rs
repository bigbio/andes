//! Decoy database generation via sequence reversal (default), seeded shuffling,
//! or none (target-only / externally-supplied decoys).

use model::protein::{Protein, ProteinDb};

/// Default decoy accession prefix.
pub const DEFAULT_DECOY_PREFIX: &str = "XXX";

/// Default seed for the reproducible shuffle decoy strategy.
pub const DEFAULT_DECOY_SEED: u64 = 42;

/// How decoy proteins are generated for the target+decoy search database.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum DecoyStrategy {
    /// Reverse each target sequence (default; deterministic, MS-GF+ compatible).
    #[default]
    Reverse,
    /// Shuffle each target sequence with a seeded, reproducible RNG. Useful when
    /// reversal produces too-similar decoys (e.g. palindromic/low-complexity).
    Shuffle,
    /// Generate no decoys: the search database is the input FASTA verbatim. For
    /// inputs that already contain decoys (accessions carrying the decoy prefix)
    /// or when FDR is computed externally downstream.
    None,
}

impl DecoyStrategy {
    /// Parse the `--decoy-strategy` CLI value. Case-insensitive.
    pub fn from_name(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "reverse" | "rev"     => Some(Self::Reverse),
            "shuffle" | "shuf"    => Some(Self::Shuffle),
            "none" | "off"        => Some(Self::None),
            _ => None,
        }
    }
}

/// Reverse each protein's sequence and prepend `<prefix>_` to its
/// accession. `prefix` is normalized: trailing `_`s stripped; empty
/// prefix → `DEFAULT_DECOY_PREFIX`.
pub fn reverse_db(db: &ProteinDb, prefix: &str) -> ProteinDb {
    let normalized = normalize_decoy_prefix(prefix);
    let proteins = db.proteins.iter().map(|p| Protein {
        accession: format!("{}_{}", normalized, p.accession),
        description: p.description.clone(),
        sequence: p.sequence.iter().rev().copied().collect(),
    }).collect();
    ProteinDb { proteins }
}

/// Concatenate target + decoy.
pub fn target_plus_decoy(target: &ProteinDb, prefix: &str) -> ProteinDb {
    let decoy = reverse_db(target, prefix);
    let mut proteins = target.proteins.clone();
    proteins.extend(decoy.proteins);
    ProteinDb { proteins }
}

/// Shuffle each protein's sequence with a per-protein, seed-derived RNG, and
/// prepend `<prefix>_` to its accession (same accession scheme as [`reverse_db`]).
///
/// Deterministic: the same `(seed, protein index, sequence)` always yields the
/// same permutation, so repeated runs and re-benchmarks are reproducible. Uses a
/// dependency-free xorshift64 RNG + Fisher–Yates; the per-protein seed mixes the
/// global seed with the protein index so identical sequences don't shuffle alike.
pub fn shuffle_db(db: &ProteinDb, prefix: &str, seed: u64) -> ProteinDb {
    let normalized = normalize_decoy_prefix(prefix);
    let proteins = db
        .proteins
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let mut seq = p.sequence.clone();
            // Per-protein RNG state: mix the global seed with the index via the
            // splitmix64 multiplier so adjacent proteins get well-separated
            // streams. Never let the xorshift state be zero (it would stick).
            let mut state = seed ^ (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
            if state == 0 {
                state = 0xD1B5_4A32_D192_ED03;
            }
            // Fisher–Yates from the back, drawing each index from xorshift64.
            for j in (1..seq.len()).rev() {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                let k = (state % (j as u64 + 1)) as usize;
                seq.swap(j, k);
            }
            Protein {
                accession: format!("{}_{}", normalized, p.accession),
                description: p.description.clone(),
                sequence: seq,
            }
        })
        .collect();
    ProteinDb { proteins }
}

/// Build the searched database from the target proteins according to `strategy`.
/// `Reverse`/`Shuffle` append a 1:1 decoy per target (target proteins occupy the
/// first half); `None` returns the target db unchanged. `seed` is only used by
/// `Shuffle`. Decoy membership downstream is determined by the accession prefix
/// (not index), so all three strategies interoperate with the same labeling.
pub fn build_search_db(
    target: &ProteinDb,
    prefix: &str,
    strategy: DecoyStrategy,
    seed: u64,
) -> ProteinDb {
    match strategy {
        DecoyStrategy::Reverse => target_plus_decoy(target, prefix),
        DecoyStrategy::Shuffle => {
            let decoy = shuffle_db(target, prefix, seed);
            let mut proteins = target.proteins.clone();
            proteins.extend(decoy.proteins);
            ProteinDb { proteins }
        }
        DecoyStrategy::None => target.clone(),
    }
}

/// Normalize a user-supplied decoy accession prefix: trim whitespace,
/// strip trailing `_`, and fall back to [`DEFAULT_DECOY_PREFIX`] when empty.
pub fn normalize_decoy_prefix(prefix: &str) -> String {
    let trimmed = prefix.trim().trim_end_matches('_');
    if trimmed.is_empty() {
        DEFAULT_DECOY_PREFIX.to_string()
    } else {
        trimmed.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_db(proteins: &[(&str, &[u8])]) -> ProteinDb {
        ProteinDb {
            proteins: proteins.iter().map(|(acc, seq)| Protein {
                accession: acc.to_string(),
                description: String::new(),
                sequence: seq.to_vec(),
            }).collect(),
        }
    }

    #[test]
    fn reverse_db_reverses_sequences() {
        let db = make_db(&[("P1", b"MKWV"), ("P2", b"AGCT")]);
        let decoy = reverse_db(&db, "XXX");
        assert_eq!(decoy.len(), 2);
        assert_eq!(decoy.proteins[0].sequence, b"VWKM");
        assert_eq!(decoy.proteins[1].sequence, b"TCGA");
    }

    #[test]
    fn reverse_db_prepends_prefix() {
        let db = make_db(&[("P1", b"AB")]);
        let decoy = reverse_db(&db, "XXX");
        assert_eq!(decoy.proteins[0].accession, "XXX_P1");
    }

    #[test]
    fn reverse_db_strips_trailing_underscores_in_prefix() {
        let db = make_db(&[("P1", b"AB")]);
        let decoy = reverse_db(&db, "XXX_");
        assert_eq!(decoy.proteins[0].accession, "XXX_P1");
    }

    #[test]
    fn reverse_db_empty_prefix_uses_default() {
        let db = make_db(&[("P1", b"AB")]);
        let decoy = reverse_db(&db, "");
        assert_eq!(decoy.proteins[0].accession, "XXX_P1");
    }

    #[test]
    fn reverse_db_preserves_description() {
        let mut db = make_db(&[("P1", b"AB")]);
        db.proteins[0].description = "Some description".into();
        let decoy = reverse_db(&db, "XXX");
        assert_eq!(decoy.proteins[0].description, "Some description");
    }

    #[test]
    fn target_plus_decoy_concats() {
        let target = make_db(&[("P1", b"AB"), ("P2", b"CD")]);
        let combined = target_plus_decoy(&target, "XXX");
        assert_eq!(combined.len(), 4);
        assert_eq!(combined.proteins[0].accession, "P1");
        assert_eq!(combined.proteins[2].accession, "XXX_P1");
        assert_eq!(combined.proteins[2].sequence, b"BA");
    }

    #[test]
    fn shuffle_db_is_deterministic_for_a_fixed_seed() {
        let db = make_db(&[("P1", b"MKWVLPASTNDE"), ("P2", b"AGCTFYHRQGIL")]);
        let a = shuffle_db(&db, "XXX", DEFAULT_DECOY_SEED);
        let b = shuffle_db(&db, "XXX", DEFAULT_DECOY_SEED);
        assert_eq!(a.proteins[0].sequence, b.proteins[0].sequence, "same seed ⇒ same shuffle");
        assert_eq!(a.proteins[1].sequence, b.proteins[1].sequence);
        // A different seed yields a different permutation (overwhelmingly likely
        // for a 12-mer), so shuffles are actually seed-dependent.
        let c = shuffle_db(&db, "XXX", DEFAULT_DECOY_SEED + 1);
        assert_ne!(a.proteins[0].sequence, c.proteins[0].sequence);
        // Shuffle is a permutation: same multiset of residues, prefixed accession.
        let mut orig = db.proteins[0].sequence.clone();
        let mut shuf = a.proteins[0].sequence.clone();
        orig.sort_unstable();
        shuf.sort_unstable();
        assert_eq!(orig, shuf, "shuffle preserves residue composition");
        assert_eq!(a.proteins[0].accession, "XXX_P1");
    }

    #[test]
    fn build_search_db_reverse_matches_target_plus_decoy() {
        let target = make_db(&[("P1", b"ABCD"), ("P2", b"EFGH")]);
        let built = build_search_db(&target, "XXX", DecoyStrategy::Reverse, DEFAULT_DECOY_SEED);
        let legacy = target_plus_decoy(&target, "XXX");
        // Reverse strategy is bit-identical to the legacy path (parity).
        assert_eq!(built.len(), legacy.len());
        for (a, b) in built.proteins.iter().zip(legacy.proteins.iter()) {
            assert_eq!(a.accession, b.accession);
            assert_eq!(a.sequence, b.sequence);
        }
    }

    #[test]
    fn build_search_db_none_is_target_only() {
        let target = make_db(&[("P1", b"ABCD"), ("P2", b"EFGH")]);
        let built = build_search_db(&target, "XXX", DecoyStrategy::None, DEFAULT_DECOY_SEED);
        assert_eq!(built.len(), 2, "None appends no decoys");
        assert_eq!(built.proteins[0].accession, "P1");
        assert_eq!(built.proteins[1].accession, "P2");
    }

    #[test]
    fn decoy_strategy_parses_names() {
        assert_eq!(DecoyStrategy::from_name("reverse"), Some(DecoyStrategy::Reverse));
        assert_eq!(DecoyStrategy::from_name("Shuffle"), Some(DecoyStrategy::Shuffle));
        assert_eq!(DecoyStrategy::from_name("NONE"), Some(DecoyStrategy::None));
        assert_eq!(DecoyStrategy::from_name("bogus"), None);
        assert_eq!(DecoyStrategy::default(), DecoyStrategy::Reverse);
    }
}
