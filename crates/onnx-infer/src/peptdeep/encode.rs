//! PeptDeep amino-acid indexing, matching `peptdeep.model.featurize.get_batch_aa_indices`.
//!
//! PeptDeep indexes residues by ASCII code (`'A'` → 1, `'B'` → 2, …, `'Z'` → 26;
//! 0 is padding) and pads one zero at the N-terminus and one at the C-terminus,
//! so a peptide of `n` residues becomes `n + 2` indices. Those two pad slots are
//! also where N- and C-terminal modification features go in `mod_x`.

use crate::error::{Error, Result};

/// Mod-element channel count in PeptDeep pretrained v3 (`model_const["mod_elements"]`).
/// `export_peptdeep_to_onnx.py` records the exported graph's count in manifest
/// `extras.mod_hidden`; nothing here reads it. For another count use the packers
/// that take `mod_dim` (`pack_rt_batch`, `pack_ms2_batch`, `pack_rt_inputs_with_mod_dim`).
pub const MOD_ELEMENT_COUNT: usize = 109;

/// PeptDeep index of one residue: ASCII position in `A..=Z`, 1-based; a lowercase
/// letter is read as its uppercase. `None` for anything that is not an ASCII letter.
pub fn aa_index(c: char) -> Option<i64> {
    let c = c.to_ascii_uppercase();
    c.is_ascii_uppercase().then(|| c as i64 - 'A' as i64 + 1)
}

/// Encode a plain amino-acid sequence (no flanks, no modification text) into
/// PeptDeep indices of length `seq_len + 2`: a zero N-terminal pad, the
/// residues, a zero C-terminal pad, and zero right-padding up to `seq_len`.
///
/// PeptDeep predicts each batch at one peptide length and does not mask
/// padding, so right-padding a shorter peptide (`seq_len > seq.len()`) changes
/// its prediction. Batch peptides of equal length with `seq_len == seq.len()`
/// to reproduce peptdeep's output.
///
/// A character that is not an ASCII letter, or a peptide longer than `seq_len`,
/// is an error; lowercase letters are read as uppercase.
pub fn encode_aa_sequence(seq: &str, seq_len: usize) -> Result<Vec<i64>> {
    if seq_len == 0 {
        return Err(Error::InvalidArg("seq_len must be ≥ 1".into()));
    }
    let len = seq.chars().count();
    if len > seq_len {
        return Err(Error::InvalidArg(format!(
            "peptide {seq:?} has {len} residues, longer than seq_len {seq_len}"
        )));
    }
    let mut out = vec![0i64; seq_len + 2];
    for (i, c) in seq.chars().enumerate() {
        out[i + 1] = aa_index(c).ok_or_else(|| {
            Error::InvalidArg(format!(
                "residue `{c}` in peptide {seq:?} is not in the PeptDeep vocabulary"
            ))
        })?;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encodes_ascii_indices_with_terminal_pads() {
        let v = encode_aa_sequence("PEPTIDE", 7).unwrap();
        assert_eq!(v, vec![0, 16, 5, 16, 20, 9, 4, 5, 0]);
        let padded = encode_aa_sequence("AC", 4).unwrap();
        assert_eq!(padded, vec![0, 1, 3, 0, 0, 0]);
    }

    #[test]
    fn rejects_non_letters() {
        for seq in ["PEP1IDE", "PEP-IDE", "M+15.995K"] {
            assert!(
                encode_aa_sequence(seq, 16).is_err(),
                "{seq} must be rejected"
            );
        }
    }

    #[test]
    fn rejects_peptides_longer_than_seq_len() {
        assert!(encode_aa_sequence("PEPTIDE", 6).is_err());
        assert!(encode_aa_sequence("PEPTIDE", 7).is_ok());
    }
}
