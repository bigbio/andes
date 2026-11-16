use crate::error::{Error, Result};

pub const MOD_ELEMENT_COUNT: usize = 109;

pub fn aa_index(c: char) -> Option<i64> {
    let c = c.to_ascii_uppercase();
    c.is_ascii_uppercase().then(|| c as i64 - 'A' as i64 + 1)
}

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
