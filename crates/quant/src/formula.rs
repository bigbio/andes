use model::amino_acid::standard_composition;
use model::isotope::{isotope_envelope_and_share_from_formula, isotope_envelope_from_formula};
use model::mass::{C, H, N, O, S};

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Formula {
    pub c: u32,
    pub h: u32,
    pub n: u32,
    pub o: u32,
    pub s: u32,
}

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

    pub fn glycan(hex: u8, hexnac: u8, fuc: u8, neuac: u8, neugc: u8) -> Self {
        let h = hex as u32;
        let n = hexnac as u32;
        let f = fuc as u32;
        let a = neuac as u32;
        let g = neugc as u32;
        Formula {
            c: 6 * h + 8 * n + 6 * f + 11 * a + 11 * g,
            h: 10 * h + 13 * n + 10 * f + 17 * a + 17 * g,
            n: n + a + g,
            o: 5 * h + 5 * n + 4 * f + 8 * a + 9 * g,
            s: 0,
        }
    }

    pub fn monoisotopic_mass(&self) -> f64 {
        self.c as f64 * C
            + self.h as f64 * H
            + self.n as f64 * N
            + self.o as f64 * O
            + self.s as f64 * S
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

    pub fn envelope(&self, n_isotopes: usize) -> Vec<f64> {
        isotope_envelope_from_formula(self.c, self.h, self.n, self.o, self.s, n_isotopes)
    }

    pub fn envelope_and_share(&self, n_isotopes: usize) -> (Vec<f64>, f64) {
        isotope_envelope_and_share_from_formula(self.c, self.h, self.n, self.o, self.s, n_isotopes)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn peptide_formula_sums_residues_and_water() {
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
    fn glycan_formula_mass_matches_the_residue_masses() {
        let f = Formula::glycan(1, 0, 0, 0, 0);
        assert!((f.monoisotopic_mass() - 162.052824).abs() < 5e-6, "{f:?}");
        let f = Formula::glycan(0, 1, 0, 0, 0);
        assert!((f.monoisotopic_mass() - 203.079373).abs() < 5e-6, "{f:?}");
        let f = Formula::glycan(0, 0, 1, 0, 0);
        assert!((f.monoisotopic_mass() - 146.057909).abs() < 5e-6, "{f:?}");
        let f = Formula::glycan(0, 0, 0, 1, 0);
        assert!((f.monoisotopic_mass() - 291.095417).abs() < 5e-6, "{f:?}");
        let f = Formula::glycan(0, 0, 0, 0, 1);
        assert!((f.monoisotopic_mass() - 307.090331).abs() < 5e-6, "{f:?}");
    }

    #[test]
    fn a_glycopeptide_envelope_is_flatter_than_the_peptide_alone() {
        let bare = Formula::peptide(b"IGNTTNA".iter().map(|&r| (r, 0.0)));
        let mut intact = bare;
        intact.add(Formula::glycan(8, 2, 0, 0, 0));
        let a = bare.envelope(6);
        let b = intact.envelope(6);
        assert!(b[1] / b[0] > a[1] / a[0]);
        let mut big = Formula::peptide(b"IGNTTNA".iter().map(|&r| (r, 0.0)));
        big.add(Formula::averagine(2500.0));
        let e = big.envelope(8);
        assert!(e[2] > 0.1 * e[0], "{e:?}");
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
    fn share_of_the_first_four_isotopes_falls_with_mass() {
        let share = |mass: f64| Formula::averagine(mass).envelope_and_share(4).1;
        for (mass, lo, hi) in [
            (1000.0, 0.99, 1.0),
            (2000.0, 0.94, 0.96),
            (3000.0, 0.86, 0.88),
            (4000.0, 0.73, 0.76),
        ] {
            let s = share(mass);
            assert!(s > lo && s < hi, "{mass} Da: {s}");
        }
        let f = Formula::peptide(b"PEPTIDEK".iter().map(|&r| (r, 0.0)));
        let (env, s) = f.envelope_and_share(4);
        assert_eq!(env, f.envelope(4));
        assert!(s > 0.99 && s < 1.0, "{s}");
        assert!(f.envelope_and_share(5).1 > s);
    }

    #[test]
    fn envelope_of_a_large_peptide_peaks_later() {
        let f = Formula::averagine(4000.0);
        let env = f.envelope(5);
        assert!(env[1] > env[0] || env[2] > env[0]);
    }
}
