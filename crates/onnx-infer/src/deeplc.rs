use crate::error::{Error, Result};
use crate::tensor::NamedTensor;

pub const PADDING_LENGTH: usize = 60;
pub const MIN_LENGTH: usize = 4;
pub const N_ATOMS: usize = 6;
pub const N_RESIDUES: usize = 20;
pub const GLOBAL_DIM: usize = 67;

const RESIDUES: &str = "KRPTNAQVSGILCMHFYWED";

fn residue_composition(c: char) -> Option<[i32; N_ATOMS]> {
    Some(match c {
        'A' => [3, 5, 1, 1, 0, 0],
        'R' => [6, 12, 4, 1, 0, 0],
        'N' => [4, 6, 2, 2, 0, 0],
        'D' => [4, 5, 1, 3, 0, 0],
        'C' => [3, 5, 1, 1, 1, 0],
        'E' => [5, 7, 1, 3, 0, 0],
        'Q' => [5, 8, 2, 2, 0, 0],
        'G' => [2, 3, 1, 1, 0, 0],
        'H' => [6, 7, 3, 1, 0, 0],
        'I' | 'L' => [6, 11, 1, 1, 0, 0],
        'K' => [6, 12, 2, 1, 0, 0],
        'M' => [5, 9, 1, 1, 1, 0],
        'F' => [9, 9, 1, 1, 0, 0],
        'P' => [5, 7, 1, 1, 0, 0],
        'S' => [3, 5, 1, 2, 0, 0],
        'T' => [4, 7, 1, 2, 0, 0],
        'W' => [11, 10, 2, 1, 0, 0],
        'Y' => [9, 9, 1, 2, 0, 0],
        'V' => [5, 9, 1, 1, 0, 0],
        _ => return None,
    })
}

#[derive(Debug, Clone, Default)]
pub struct DeepLcPeptide<'a> {
    pub sequence: &'a str,
    pub mods: Vec<(usize, [i32; N_ATOMS])>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct DeepLcFeatures {
    pub x_atom: Vec<f32>,
    pub x_global: Vec<f32>,
    pub x_one_hot: Vec<f32>,
}

fn positional_rows(i: usize, n: usize) -> impl Iterator<Item = usize> {
    let pos = (i <= 3).then_some(i + 4);
    let neg = (n - i <= 4).then(|| 4 - (n - i));
    pos.into_iter().chain(neg)
}

pub fn encode(p: &DeepLcPeptide<'_>) -> Result<DeepLcFeatures> {
    let seq: Vec<char> = p.sequence.chars().collect();
    let n = seq.len();
    if !(MIN_LENGTH..=PADDING_LENGTH).contains(&n) {
        return Err(Error::InvalidArg(format!(
            "peptide {:?} has {n} residues; DeepLC encoding covers {MIN_LENGTH}..={PADDING_LENGTH}",
            p.sequence
        )));
    }
    let mut atom = vec![[0i32; N_ATOMS]; n];
    let mut one_hot = vec![0f32; PADDING_LENGTH * N_RESIDUES];
    for (i, &c) in seq.iter().enumerate() {
        atom[i] = residue_composition(c).ok_or_else(|| {
            Error::InvalidArg(format!(
                "residue `{c}` in {:?} is not in the DeepLC alphabet",
                p.sequence
            ))
        })?;
        let k = RESIDUES.find(c).expect("alphabet");
        one_hot[i * N_RESIDUES + k] = 1.0;
    }
    let mut positional = [[0i32; N_ATOMS]; 8];
    for i in (0..4).chain(n - 4..n) {
        for r in positional_rows(i, n) {
            positional[r] = atom[i];
        }
    }
    let mut terminal = [[0i32; N_ATOMS]; 2];
    for &(site, delta) in &p.mods {
        let i = match site {
            0 => 0,
            s if s <= n => s - 1,
            s if s == n + 1 => n - 1,
            s => {
                return Err(Error::InvalidArg(format!(
                    "mod site {s} outside peptide {:?} (n = {n})",
                    p.sequence
                )));
            }
        };
        let add = |dst: &mut [i32; N_ATOMS]| dst.iter_mut().zip(delta).for_each(|(d, v)| *d += v);
        add(&mut atom[i]);
        for r in positional_rows(i, n) {
            add(&mut positional[r]);
        }
        if site == 0 {
            add(&mut terminal[0]);
        } else if site == n + 1 {
            add(&mut terminal[1]);
        }
    }
    let mut x_atom = vec![0f32; PADDING_LENGTH * N_ATOMS];
    let mut x_global = Vec::with_capacity(GLOBAL_DIM);
    let mut sum = [0i32; N_ATOMS];
    for (i, a) in atom.iter().enumerate() {
        for k in 0..N_ATOMS {
            x_atom[i * N_ATOMS + k] = a[k] as f32;
            sum[k] += a[k];
        }
    }
    x_global.extend(sum.iter().map(|&v| v as f32));
    x_global.push(n as f32);
    x_global.extend(positional.iter().flatten().map(|&v| v as f32));
    x_global.extend(terminal.iter().flatten().map(|&v| v as f32));
    debug_assert_eq!(x_global.len(), GLOBAL_DIM);
    Ok(DeepLcFeatures {
        x_atom,
        x_global,
        x_one_hot: one_hot,
    })
}

pub fn pack_batch(peptides: &[DeepLcPeptide<'_>]) -> Result<[NamedTensor; 3]> {
    if peptides.is_empty() {
        return Err(Error::InvalidArg("empty peptide batch".into()));
    }
    let b = peptides.len();
    let (mut a, mut g, mut h) = (Vec::new(), Vec::new(), Vec::new());
    for p in peptides {
        let f = encode(p)?;
        a.extend(f.x_atom);
        g.extend(f.x_global);
        h.extend(f.x_one_hot);
    }
    Ok([
        NamedTensor::new("x_atom", vec![b, PADDING_LENGTH, N_ATOMS], a)?,
        NamedTensor::new("x_global", vec![b, GLOBAL_DIM], g)?,
        NamedTensor::new("x_one_hot", vec![b, PADDING_LENGTH, N_RESIDUES], h)?,
    ])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pep(s: &str) -> DeepLcPeptide<'_> {
        DeepLcPeptide {
            sequence: s,
            mods: Vec::new(),
        }
    }

    #[test]
    fn global_layout_for_unmodified_peptide() {
        let f = encode(&pep("PEPTIDEK")).unwrap();
        assert_eq!(&f.x_global[..7], &[40.0, 63.0, 9.0, 15.0, 0.0, 0.0, 8.0]);
        assert_eq!(&f.x_global[7..13], &[6.0, 11.0, 1.0, 1.0, 0.0, 0.0]);
        assert_eq!(&f.x_global[7 + 24..7 + 30], &[5.0, 7.0, 1.0, 1.0, 0.0, 0.0]);
        assert_eq!(&f.x_global[7 + 42..7 + 48], &[4.0, 7.0, 1.0, 2.0, 0.0, 0.0]);
        assert!(f.x_global[55..].iter().all(|&v| v == 0.0));
        assert_eq!(f.x_one_hot[2], 1.0, "P is one-hot column 2");
        assert!(f.x_atom[8 * N_ATOMS..].iter().all(|&v| v == 0.0));
    }

    #[test]
    fn terminal_mod_reaches_residue_positional_and_terminal_blocks() {
        let acetyl = [2, 2, 0, 1, 0, 0];
        let p = DeepLcPeptide {
            sequence: "SAEDK",
            mods: vec![(0, acetyl)],
        };
        let f = encode(&p).unwrap();
        let base = encode(&pep("SAEDK")).unwrap();
        assert_eq!(f.x_atom[0] - base.x_atom[0], 2.0);
        assert_eq!(f.x_global[7 + 24] - base.x_global[7 + 24], 2.0);
        assert_eq!(&f.x_global[55..61], &[2.0, 2.0, 0.0, 1.0, 0.0, 0.0]);
    }

    #[test]
    fn short_peptide_residue_fills_both_rows() {
        let ox = [0, 0, 0, 1, 0, 0];
        let p = DeepLcPeptide {
            sequence: "MEEK",
            mods: vec![(1, ox)],
        };
        let f = encode(&p).unwrap();
        let base = encode(&pep("MEEK")).unwrap();
        assert_eq!(
            f.x_global[7 + 3] - base.x_global[7 + 3],
            1.0,
            "row 0 (−4) O"
        );
        assert_eq!(
            f.x_global[7 + 27] - base.x_global[7 + 27],
            1.0,
            "row 4 (0) O"
        );
    }

    #[test]
    fn rejects_unsupported() {
        assert!(encode(&pep("PEK")).is_err());
        assert!(encode(&pep(&"A".repeat(61))).is_err());
        assert!(encode(&pep("PEPXIDE")).is_err());
    }
}
