//! Peptide elemental formula and theoretical precursor isotope envelope.
//!
//! Standard residues use their exact compositions; a modification contributes
//! an averagine-shaped allotment of atoms for its mass delta (the delta is
//! known, its formula is not). The envelope is the exact multinomial fold of
//! the per-element abundance ladders (`model::isotope`).

use model::amino_acid::standard_composition;
use model::isotope::isotope_envelope_from_formula;

/// Atom counts (C, H, N, O, S).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Formula {
    pub c: u32,
    pub h: u32,
    pub n: u32,
    pub o: u32,
    pub s: u32,
}

/// Averagine (Senko 1995): C4.9384 H7.7583 N1.3577 O1.4773 S0.0417 per 111.1254 Da.
const AVERAGINE_MASS: f64 = 111.1254;
const AVERAGINE: [f64; 5] = [4.9384, 7.7583, 1.3577, 1.4773, 0.0417];

impl Formula {
    pub fn water() -> Self {
        Formula {
            c: 0,
            h: 2,
            n: 0,
            o: 1,
            s: 0,
        }
    }

    /// Averagine atoms for a mass (rounded to integers, never negative).
    pub fn averagine(mass: f64) -> Self {
        let units = (mass / AVERAGINE_MASS).max(0.0);
        let r = |k: usize| (AVERAGINE[k] * units).round().max(0.0) as u32;
        Formula {
            c: r(0),
            h: r(1),
            n: r(2),
            o: r(3),
            s: r(4),
        }
    }

    pub fn add(&mut self, other: Formula) {
        self.c += other.c;
        self.h += other.h;
        self.n += other.n;
        self.o += other.o;
        self.s += other.s;
    }

    pub fn sub_saturating(&mut self, other: Formula) {
        self.c = self.c.saturating_sub(other.c);
        self.h = self.h.saturating_sub(other.h);
        self.n = self.n.saturating_sub(other.n);
        self.o = self.o.saturating_sub(other.o);
        self.s = self.s.saturating_sub(other.s);
    }

    /// Formula of a peptide given its residues (one-letter codes) and the
    /// modification mass delta carried by each residue (0.0 when unmodified).
    /// Unknown residues (`U`, `X`, ...) fall back to averagine for a 110 Da
    /// residue.
    pub fn peptide<I>(residues: I) -> Self
    where
        I: IntoIterator<Item = (u8, f64)>,
    {
        let mut f = Formula::water();
        for (res, delta) in residues {
            match standard_composition(res) {
                Some((c, h, n, o, s)) => f.add(Formula { c, h, n, o, s }),
                None => f.add(Formula::averagine(110.0)),
            }
            if delta > 0.0 {
                f.add(Formula::averagine(delta));
            } else if delta < 0.0 {
                f.sub_saturating(Formula::averagine(-delta));
            }
        }
        f
    }

    /// First `n_isotopes` relative abundances of the precursor envelope,
    /// normalised to sum 1.
    pub fn envelope(&self, n_isotopes: usize) -> Vec<f64> {
        isotope_envelope_from_formula(self.c, self.h, self.n, self.o, self.s, n_isotopes)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn peptide_formula_sums_residues_and_water() {
        // G + A: C5 H8 N2 O2 + H2O = C5 H10 N2 O3
        let f = Formula::peptide([(b'G', 0.0), (b'A', 0.0)]);
        assert_eq!(
            f,
            Formula {
                c: 5,
                h: 10,
                n: 2,
                o: 3,
                s: 0
            }
        );
    }

    #[test]
    fn modification_adds_averagine_atoms() {
        let plain = Formula::peptide([(b'K', 0.0)]);
        let tmt = Formula::peptide([(b'K', 229.162932)]);
        assert!(tmt.c > plain.c && tmt.n > plain.n);
        let loss = Formula::peptide([(b'K', -17.0)]);
        assert!(loss.c <= plain.c);
    }

    #[test]
    fn envelope_of_a_small_peptide_is_mono_dominated() {
        let f = Formula::peptide(b"PEPTIDEK".iter().map(|&r| (r, 0.0)));
        let env = f.envelope(4);
        assert_eq!(env.len(), 4);
        assert!((env.iter().sum::<f64>() - 1.0).abs() < 1e-9);
        assert!(env[0] > env[1] && env[1] > env[2]);
    }

    #[test]
    fn envelope_of_a_large_peptide_peaks_later() {
        let f = Formula::averagine(4000.0);
        let env = f.envelope(5);
        assert!(env[1] > env[0] || env[2] > env[0]);
    }
}
