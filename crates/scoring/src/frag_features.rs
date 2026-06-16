//! Peptide-CONDITIONED per-fragment features for the v3 intensity regressor.
//! The SINGLE source of the feature vector — imported by both the trainer's
//! frag-dataset builder and `intensity_signal` — so train/infer features match
//! by construction (mirrors `PeakFeatureCtx::for_spectrum`).

use model::peptide::Peptide;

use crate::scoring::fragment_ions::IonKind;

pub const FEAT_ION_TYPE: usize = 0;
pub const FEAT_CHARGE: usize = 1;
pub const FEAT_NFLANK: usize = 2;
pub const FEAT_CFLANK: usize = 3;
pub const FEAT_PROLINE_FLANK: usize = 4;
pub const FEAT_POS_FRAC: usize = 5;
pub const FEAT_PEP_LEN: usize = 6;
pub const FEAT_NFLANK_MOD: usize = 7;
pub const FEAT_CFLANK_MOD: usize = 8;
pub const FEAT_NCE: usize = 9;
pub const N_FRAG_FEATURES: usize = 10;

/// Resolve the (N-side index, C-side index) of the two residues flanking
/// the cleavage site for a b/y ion at 1-based `position` in a peptide of
/// length `n`. Returns `None` for terminal/degenerate positions.
fn flank_indices(n: usize, kind: IonKind, position: u32) -> Option<(usize, usize)> {
    let i = position as usize;
    if i < 1 || i >= n {
        return None;
    }
    match kind {
        IonKind::B => Some((i - 1, i)),
        IonKind::Y => {
            let left = n - i;
            Some((left - 1, left))
        }
    }
}

/// Map a residue byte to its 0-based integer index (`residue − b'A'`).
#[inline]
fn res_idx(b: u8) -> f32 {
    b.wrapping_sub(b'A') as f32
}

/// Return the modification mass delta (as f32) for residue at `idx`, or 0.0
/// if the residue carries no modification.
#[inline]
fn mod_delta(p: &Peptide, idx: usize) -> f32 {
    p.residues[idx].mod_.as_ref().map_or(0.0, |m| m.mass_delta) as f32
}

/// Feature vector for one annotated b/y ion. `position` is 1-based, `charge`
/// the fragment charge, `nce` the parsed NCE (0.0 when unknown). Returns an
/// all-zero vector for out-of-range positions (terminal/degenerate ions).
pub fn extract_frag_features(
    p: &Peptide,
    kind: IonKind,
    position: u32,
    charge: u8,
    nce: f32,
) -> [f32; N_FRAG_FEATURES] {
    let n = p.residues.len();
    let mut f = [0.0f32; N_FRAG_FEATURES];
    let (ni, ci) = match flank_indices(n, kind, position) {
        Some(v) => v,
        None => return f,
    };
    f[FEAT_ION_TYPE] = match kind {
        IonKind::B => 0.0,
        IonKind::Y => 1.0,
    };
    f[FEAT_CHARGE] = charge as f32;
    f[FEAT_NFLANK] = res_idx(p.residues[ni].residue);
    f[FEAT_CFLANK] = res_idx(p.residues[ci].residue);
    f[FEAT_PROLINE_FLANK] =
        if p.residues[ni].residue == b'P' || p.residues[ci].residue == b'P' {
            1.0
        } else {
            0.0
        };
    f[FEAT_POS_FRAC] = position as f32 / n as f32;
    f[FEAT_PEP_LEN] = n as f32;
    f[FEAT_NFLANK_MOD] = mod_delta(p, ni);
    f[FEAT_CFLANK_MOD] = mod_delta(p, ci);
    f[FEAT_NCE] = nce;
    f
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::amino_acid::AminoAcid;
    use model::modification::{ModLocation, Modification, ResidueSpec};
    use model::peptide::Peptide;

    use crate::scoring::fragment_ions::IonKind;

    fn pep(seq: &str) -> Peptide {
        Peptide::new(
            seq.bytes().map(|b| AminoAcid::standard(b).unwrap()).collect(),
            b'K',
            b'R',
        )
    }

    fn pep_mod(seq: &str, pos1: usize, delta: f64) -> Peptide {
        let mut r: Vec<AminoAcid> =
            seq.bytes().map(|b| AminoAcid::standard(b).unwrap()).collect();
        let res = r[pos1 - 1].residue;
        let m = Modification {
            name: "t".into(),
            mass_delta: delta,
            residue: ResidueSpec::Specific(res),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: vec![],
            loss_class: 0,
        };
        r[pos1 - 1] = r[pos1 - 1].clone().with_mod(m);
        Peptide::new(r, b'K', b'R')
    }

    #[test]
    fn frag_features_stable_shape_and_context() {
        // n=7; b2 cleaves after residue 2 (flanks E|P), P is C-flank of PEPTIDE.
        let p = pep("PEPTIDE");
        let f = extract_frag_features(&p, IonKind::B, 2, 1, 0.0);
        assert_eq!(f.len(), N_FRAG_FEATURES);
        assert_eq!(f[FEAT_ION_TYPE], 0.0);
        assert!((f[FEAT_POS_FRAC] - 2.0 / 7.0).abs() < 1e-6);
        assert_eq!(f[FEAT_PROLINE_FLANK], 1.0);
        assert_eq!(f[FEAT_PEP_LEN], 7.0);
    }

    #[test]
    fn frag_features_reflect_modification() {
        let plain = extract_frag_features(&pep("PEACDEK"), IonKind::B, 3, 1, 0.0);
        let modded =
            extract_frag_features(&pep_mod("PEACDEK", 3, 57.02146), IonKind::B, 3, 1, 0.0);
        assert_ne!(plain[FEAT_NFLANK_MOD], modded[FEAT_NFLANK_MOD]);
    }
}
