//! Peptide-CONDITIONED per-fragment features for the v3 intensity regressor.
//! The SINGLE source of the feature vector — imported by both the trainer's
//! frag-dataset builder and `intensity_signal` — so train/infer features match
//! by construction (mirrors `PeakFeatureCtx::for_spectrum`).

use model::peptide::Peptide;

use crate::scoring::fragment_ions::IonKind;

pub const FEAT_ION_TYPE: usize = 0;
pub const FEAT_PRECURSOR_CHARGE: usize = 1;
pub const FEAT_NFLANK: usize = 2;
pub const FEAT_CFLANK: usize = 3;
pub const FEAT_PROLINE_FLANK: usize = 4;
pub const FEAT_POS_FRAC: usize = 5;
pub const FEAT_PEP_LEN: usize = 6;
pub const FEAT_NFLANK_MOD: usize = 7;
pub const FEAT_CFLANK_MOD: usize = 8;
pub const FEAT_NCE: usize = 9;
// New features (appended; existing indices unchanged)
pub const FEAT_N_BASIC: usize = 10;
pub const FEAT_PROTON_MOBILITY: usize = 11;
pub const FEAT_BASIC_NFLANK: usize = 12;
pub const FEAT_BASIC_CFLANK: usize = 13;
pub const FEAT_HYDRO_NFLANK: usize = 14;
pub const FEAT_HYDRO_CFLANK: usize = 15;
pub const FEAT_TOTAL_MODS: usize = 16;
pub const FEAT_NEAREST_MOD_DIST: usize = 17;
pub const N_FRAG_FEATURES: usize = 18;

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

/// Returns `true` if the residue byte is a basic residue (R, K, or H).
#[inline]
fn is_basic(b: u8) -> bool {
    matches!(b, b'R' | b'K' | b'H')
}

/// Kyte-Doolittle hydropathy for a single-letter residue byte.
/// Unknown residues return 0.0.
fn hydropathy(b: u8) -> f32 {
    match b {
        b'A' => 1.8,
        b'R' => -4.5,
        b'N' => -3.5,
        b'D' => -3.5,
        b'C' => 2.5,
        b'Q' => -3.5,
        b'E' => -3.5,
        b'G' => -0.4,
        b'H' => -3.2,
        b'I' => 4.5,
        b'L' => 3.8,
        b'K' => -3.9,
        b'M' => 1.9,
        b'F' => 2.8,
        b'P' => -1.6,
        b'S' => -0.8,
        b'T' => -0.7,
        b'W' => -0.9,
        b'Y' => -1.3,
        b'V' => 4.2,
        _ => 0.0,
    }
}

/// Feature vector for one annotated b/y ion. `position` is 1-based,
/// `precursor_charge` the PSM precursor charge (not the fragment charge),
/// `nce` the parsed NCE (0.0 when unknown). Returns an all-zero vector for
/// out-of-range positions (terminal/degenerate ions).
pub fn extract_frag_features(
    p: &Peptide,
    kind: IonKind,
    position: u32,
    precursor_charge: u8,
    nce: f32,
) -> [f32; N_FRAG_FEATURES] {
    let n = p.residues.len();
    let mut f = [0.0f32; N_FRAG_FEATURES];
    let (ni, ci) = match flank_indices(n, kind, position) {
        Some(v) => v,
        None => return f,
    };

    // ---- Original features (indices 0-9, semantics unchanged) ----
    f[FEAT_ION_TYPE] = match kind {
        IonKind::B => 0.0,
        IonKind::Y => 1.0,
    };
    f[FEAT_PRECURSOR_CHARGE] = precursor_charge as f32;
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

    // ---- New features (indices 10-17) ----

    // Proton mobility / basic-residue context (peptide-global)
    let n_basic = p.residues.iter().filter(|aa| is_basic(aa.residue)).count();
    f[FEAT_N_BASIC] = n_basic as f32;
    f[FEAT_PROTON_MOBILITY] = n_basic as f32 - precursor_charge as f32;

    // Per-flank chemistry (N-side and C-side residue of the cleavage)
    let nr = p.residues[ni].residue;
    let cr = p.residues[ci].residue;
    f[FEAT_BASIC_NFLANK] = if is_basic(nr) { 1.0 } else { 0.0 };
    f[FEAT_BASIC_CFLANK] = if is_basic(cr) { 1.0 } else { 0.0 };
    f[FEAT_HYDRO_NFLANK] = hydropathy(nr);
    f[FEAT_HYDRO_CFLANK] = hydropathy(cr);

    // PTM context
    let total_mods = p.residues.iter().filter(|aa| aa.mod_.is_some()).count();
    f[FEAT_TOTAL_MODS] = total_mods as f32;

    // Distance from cleavage site (ni, the N-side index) to nearest modified residue.
    // Sentinel = peptide length when there are no modified residues.
    let nearest_mod_dist = if total_mods == 0 {
        n as f32
    } else {
        p.residues
            .iter()
            .enumerate()
            .filter(|(_, aa)| aa.mod_.is_some())
            .map(|(idx, _)| (idx as isize - ni as isize).unsigned_abs() as f32)
            .fold(f32::MAX, f32::min)
    };
    f[FEAT_NEAREST_MOD_DIST] = nearest_mod_dist;

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
        let f = extract_frag_features(&p, IonKind::B, 2, 2, 0.0);
        assert_eq!(f.len(), N_FRAG_FEATURES);
        assert_eq!(f[FEAT_ION_TYPE], 0.0);
        assert!((f[FEAT_POS_FRAC] - 2.0 / 7.0).abs() < 1e-6);
        assert_eq!(f[FEAT_PROLINE_FLANK], 1.0);
        assert_eq!(f[FEAT_PEP_LEN], 7.0);
    }

    #[test]
    fn frag_features_reflect_modification() {
        let plain = extract_frag_features(&pep("PEACDEK"), IonKind::B, 3, 2, 0.0);
        let modded =
            extract_frag_features(&pep_mod("PEACDEK", 3, 57.02146), IonKind::B, 3, 2, 0.0);
        assert_ne!(plain[FEAT_NFLANK_MOD], modded[FEAT_NFLANK_MOD]);
    }

    #[test]
    fn frag_features_chemistry_and_ptm() {
        let p = pep("PEPTIDER"); // R is basic, at C-term
        let f = extract_frag_features(&p, IonKind::B, 2, 2, 0.0);
        assert_eq!(f[FEAT_PRECURSOR_CHARGE], 2.0);
        assert_eq!(f[FEAT_N_BASIC], 1.0);                 // one R
        assert!((f[FEAT_PROTON_MOBILITY] - (1.0 - 2.0)).abs() < 1e-6); // 1 basic - charge 2 = -1
        // b2 flanks E|P (neither basic); hydropathy of E = -3.5, P = -1.6
        assert_eq!(f[FEAT_BASIC_NFLANK], 0.0);
        assert!((f[FEAT_HYDRO_NFLANK] - (-3.5)).abs() < 1e-6);
        assert!((f[FEAT_HYDRO_CFLANK] - (-1.6)).abs() < 1e-6);
        // unmodified peptide → total mods 0, nearest-mod-dist = sentinel (len)
        assert_eq!(f[FEAT_TOTAL_MODS], 0.0);
        assert!((f[FEAT_NEAREST_MOD_DIST] - p.residues.len() as f32).abs() < 1e-6);
    }

    #[test]
    fn frag_features_ptm_distance() {
        // Cam-C on residue 3 of "PEPCDEK": at b2 (cleavage index ni=1, 0-based),
        // nearest mod is residue index 2 → distance 1; total mods 1.
        let p = pep_mod("PEPCDEK", 3, 57.02146);
        let f = extract_frag_features(&p, IonKind::B, 2, 2, 0.0);
        assert_eq!(f[FEAT_TOTAL_MODS], 1.0);
        assert!((f[FEAT_NEAREST_MOD_DIST] - 1.0).abs() < 1e-6);
    }
}
