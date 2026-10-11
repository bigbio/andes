//! Minimal PeptDeep-compatible amino-acid indexing (1-based, OpenMS-aligned).
//!
//! Full modification one-hot (`mod_x`) requires `model_const["mod_elements"]`
//! from peptdeep at export time; we only pad zeros here unless the caller
//! supplies a dense mod tensor.

use crate::error::{Error, Result};

/// Default mod-element channel count in PeptDeep pretrained v3 (`mod_elements` len).
/// Overridable via manifest `extras.mod_hidden`.
pub const MOD_ELEMENT_COUNT: usize = 109;

/// 1-based AA vocab used by OpenMS PeptDeep ONNX bindings (index 0 = pad/unknown).
///
/// Order matches common PeptDeep / alphabase residue indexing for the 20 AA +
/// a few specials; unknown residues map to 0.
pub const AA_VOCAB_1BASED: &[(char, i64)] = &[
    ('A', 1),
    ('C', 2),
    ('D', 3),
    ('E', 4),
    ('F', 5),
    ('G', 6),
    ('H', 7),
    ('I', 8),
    ('K', 9),
    ('L', 10),
    ('M', 11),
    ('N', 12),
    ('P', 13),
    ('Q', 14),
    ('R', 15),
    ('S', 16),
    ('T', 17),
    ('V', 18),
    ('W', 19),
    ('Y', 20),
];

fn aa_index(c: char) -> i64 {
    for (ch, idx) in AA_VOCAB_1BASED {
        if *ch == c {
            return *idx;
        }
    }
    0
}

/// Encode a plain amino-acid sequence into 1-based indices (no flanking dots).
///
/// Returns a length-`seq_len` vector, right-padded with 0 when `seq` is shorter,
/// or truncated when longer (PeptDeep dynamic axis usually handles true length;
/// fixed padding is for rectangular batches).
pub fn encode_aa_sequence(seq: &str, seq_len: usize) -> Result<Vec<i64>> {
    if seq_len == 0 {
        return Err(Error::InvalidArg("seq_len must be ≥ 1".into()));
    }
    let mut out = vec![0i64; seq_len];
    for (i, c) in seq.chars().take(seq_len).enumerate() {
        if !c.is_ascii_alphabetic() {
            return Err(Error::InvalidArg(format!(
                "non-alphabetic residue `{c}` in peptide {seq:?}"
            )));
        }
        out[i] = aa_index(c.to_ascii_uppercase());
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encodes_peptide_and_pads() {
        let v = encode_aa_sequence("PEPTIDE", 10).unwrap();
        assert_eq!(v.len(), 10);
        assert_eq!(v[0], 13); // P
        assert_eq!(v[6], 4); // E
        assert_eq!(v[7], 0);
    }
}
