use model::mass::{H2O, PROTON};
use model::modification::Modification;
use model::peptide::Peptide;
use model::tolerance::Tolerance;

#[derive(Debug, Clone)]
pub struct ModMark {
    pub modification: &'static Modification,
    pub diagnostic: Option<bool>,
    pub loss: Option<bool>,
}

#[derive(Debug, Clone, Default)]
pub struct SiteEvidence {
    pub unique_hit: u16,
    pub unique_miss: u16,
    pub hit_fraction: f32,
    pub marks: Vec<ModMark>,
}

fn fragment_mzs(peptide: &Peptide, charge: u8) -> Vec<f64> {
    let masses: Vec<f64> = peptide
        .residues
        .iter()
        .map(|a| a.mass + a.mod_.map_or(0.0, |m| m.mass_delta))
        .collect();
    let total: f64 = masses.iter().sum();
    let mut out = Vec::with_capacity(4 * masses.len());
    let mut prefix = 0.0;
    for m in &masses[..masses.len().saturating_sub(1)] {
        prefix += m;
        let (b, y) = (prefix + PROTON, total - prefix + H2O + PROTON);
        out.extend([b, y]);
        if charge >= 3 {
            out.extend([(b + PROTON) / 2.0, (y + PROTON) / 2.0]);
        }
    }
    out
}

fn present(peaks: &[(f64, f32)], mz: f64, tol: Tolerance) -> bool {
    let tol = tol.as_da(mz);
    let i = peaks.partition_point(|p| p.0 < mz - tol);
    i < peaks.len() && peaks[i].0 <= mz + tol
}

fn near(ions: &[f64], mz: f64, tol: f64) -> bool {
    let i = ions.partition_point(|&v| v < mz - tol);
    i < ions.len() && ions[i] <= mz + tol
}

pub struct Row<'a> {
    pub peptide: &'a Peptide,
    pub is_decoy: bool,
    pub charge: u8,
    pub precursor_mass: f64,
}

pub fn spectrum_evidence(
    peaks: &[(f64, f32)],
    rows: &[Row<'_>],
    tol: Tolerance,
) -> Vec<Option<Box<SiteEvidence>>> {
    let variable = |p: &Peptide| p.residues.iter().any(|a| a.mod_.is_some_and(|m| !m.fixed));
    if !rows.iter().any(|r| variable(r.peptide)) {
        return vec![None; rows.len()];
    }
    let ions: Vec<Vec<f64>> = rows
        .iter()
        .map(|r| fragment_mzs(r.peptide, r.charge))
        .collect();
    let same_peptide = |a: &Peptide, b: &Peptide| {
        a.residues.len() == b.residues.len()
            && a.residues
                .iter()
                .zip(&b.residues)
                .all(|(x, y)| x.residue == y.residue)
    };
    rows.iter()
        .enumerate()
        .map(|(k, row)| {
            let own = &ions[k];
            let mut others: Vec<f64> = rows
                .iter()
                .enumerate()
                .filter(|&(j, o)| {
                    j != k && o.is_decoy == row.is_decoy && same_peptide(o.peptide, row.peptide)
                })
                .flat_map(|(j, _)| ions[j].iter().copied())
                .collect();
            others.sort_by(f64::total_cmp);
            let (mut hit, mut unique_hit, mut unique_miss) = (0usize, 0u16, 0u16);
            for &mz in own {
                let h = present(peaks, mz, tol);
                hit += h as usize;
                if !near(&others, mz, 2.0 * tol.as_da(mz)) {
                    if h {
                        unique_hit = unique_hit.saturating_add(1);
                    } else {
                        unique_miss = unique_miss.saturating_add(1);
                    }
                }
            }
            let marks = row
                .peptide
                .residues
                .iter()
                .enumerate()
                .filter_map(|(i, a)| {
                    let m = a.mod_.filter(|m| !m.fixed)?;
                    if m.diagnostic_ions.is_empty() && m.neutral_losses.is_empty() {
                        return None;
                    }
                    let diagnostic = (!m.diagnostic_ions.is_empty())
                        .then(|| m.diagnostic_ions.iter().any(|&d| present(peaks, d, tol)));
                    let loss = (!m.neutral_losses.is_empty())
                        .then(|| loss_present(peaks, row, i, &m.neutral_losses, tol));
                    Some(ModMark {
                        modification: m,
                        diagnostic,
                        loss,
                    })
                })
                .collect();
            Some(Box::new(SiteEvidence {
                unique_hit,
                unique_miss,
                hit_fraction: if own.is_empty() {
                    0.0
                } else {
                    hit as f32 / own.len() as f32
                },
                marks,
            }))
        })
        .collect()
}

fn loss_present(
    peaks: &[(f64, f32)],
    row: &Row<'_>,
    at: usize,
    losses: &[f64],
    tol: Tolerance,
) -> bool {
    let masses: Vec<f64> = row
        .peptide
        .residues
        .iter()
        .map(|a| a.mass + a.mod_.map_or(0.0, |m| m.mass_delta))
        .collect();
    let total: f64 = masses.iter().sum();
    losses.iter().any(|&loss| {
        let precursor = (1..=row.charge.max(1)).any(|z| {
            let z = z as f64;
            present(peaks, (row.precursor_mass - loss + z * PROTON) / z, tol)
        });
        precursor || {
            let mut prefix = 0.0;
            (0..masses.len().saturating_sub(1)).any(|i| {
                prefix += masses[i];
                let b = (i >= at).then_some(prefix + PROTON - loss);
                let y = (i < at).then_some(total - prefix + H2O + PROTON - loss);
                b.or(y).is_some_and(|mz| present(peaks, mz, tol))
            })
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::amino_acid::AminoAcid;
    use model::modification::leak_mod;

    fn peptide(seq: &[u8], at: usize, m: &'static Modification) -> Peptide {
        let residues = seq
            .iter()
            .enumerate()
            .map(|(i, &r)| {
                let mut a = AminoAcid::standard(r).unwrap();
                if i == at {
                    a.mod_ = Some(m);
                }
                a
            })
            .collect();
        Peptide::new(residues, b'K', b'G')
    }

    fn lys(name: &str, mass: f64, ions: &[f64], losses: &[f64]) -> &'static Modification {
        leak_mod(Modification {
            diagnostic_ions: ions.to_vec(),
            neutral_losses: losses.to_vec(),
            loss_class: if losses.is_empty() { 0 } else { 255 },
            ..Modification::new(
                name,
                mass,
                model::ResidueSpec::Specific(b'K'),
                model::ModLocation::Anywhere,
                false,
            )
        })
    }

    #[test]
    fn diagnostic_ions_and_losses_mark_their_modification() {
        let ac = lys("Acetyl", 42.010565, &[143.1179], &[]);
        let me3 = lys("Trimethyl", 42.04695, &[143.1543], &[59.0735]);
        let a = peptide(b"PEPKTIDER", 3, ac);
        let t = peptide(b"PEPKTIDER", 3, me3);
        let mass = t.mass();
        let b4: f64 = a.residues[..4]
            .iter()
            .map(|x| x.mass + x.mod_.map_or(0.0, |m| m.mass_delta))
            .sum::<f64>()
            + PROTON;
        let peaks = vec![(143.1179, 1.0), (b4, 1.0)];
        let rows = [
            Row {
                peptide: &a,
                is_decoy: false,
                charge: 2,
                precursor_mass: mass,
            },
            Row {
                peptide: &t,
                is_decoy: false,
                charge: 2,
                precursor_mass: mass,
            },
        ];
        let ev = spectrum_evidence(&peaks, &rows, Tolerance::Da(0.01));
        let (ea, et) = (ev[0].as_ref().unwrap(), ev[1].as_ref().unwrap());
        assert_eq!(ea.marks[0].diagnostic, Some(true));
        assert_eq!(ea.marks[0].loss, None);
        assert_eq!(et.marks[0].diagnostic, Some(false));
        assert_eq!(et.marks[0].loss, Some(false));
        assert_eq!((ea.unique_hit, et.unique_hit), (1, 0));
        assert!(ea.unique_miss >= 4 && et.unique_miss >= 5);

        let peaks = vec![((mass - 59.0735 + 2.0 * PROTON) / 2.0, 1.0)];
        let ev = spectrum_evidence(&peaks, &rows, Tolerance::Da(0.01));
        assert_eq!(ev[1].as_ref().unwrap().marks[0].loss, Some(true));
    }

    #[test]
    fn a_spectrum_without_variable_mods_gets_no_evidence() {
        let p = Peptide::new(
            b"PEPTIDEK"
                .iter()
                .map(|&r| AminoAcid::standard(r).unwrap())
                .collect(),
            b'K',
            b'G',
        );
        let rows = [Row {
            peptide: &p,
            is_decoy: false,
            charge: 2,
            precursor_mass: p.mass(),
        }];
        assert!(spectrum_evidence(&[], &rows, Tolerance::Da(0.01))[0].is_none());
    }

    #[test]
    fn a_ppm_tolerance_separates_near_isobaric_ions() {
        let ac = lys("Acetyl", 42.010565, &[126.0913], &[]);
        let me3 = lys("Trimethyl", 42.04695, &[], &[]);
        let a = peptide(b"PEPKTIDER", 3, ac);
        let t = peptide(b"PEPKTIDER", 3, me3);
        let b4: f64 = a.residues[..4]
            .iter()
            .map(|x| x.mass + x.mod_.map_or(0.0, |m| m.mass_delta))
            .sum::<f64>()
            + PROTON;
        let peaks = vec![(126.1277, 1.0), (b4, 1.0)];
        let rows = [
            Row {
                peptide: &a,
                is_decoy: false,
                charge: 2,
                precursor_mass: a.mass(),
            },
            Row {
                peptide: &t,
                is_decoy: false,
                charge: 2,
                precursor_mass: t.mass(),
            },
        ];
        let wide = spectrum_evidence(&peaks, &rows, Tolerance::Da(0.5));
        let w = wide[0].as_ref().unwrap();
        assert_eq!(w.marks[0].diagnostic, Some(true));
        assert_eq!((w.unique_hit, w.unique_miss), (0, 0));
        let tight = spectrum_evidence(&peaks, &rows, Tolerance::Ppm(20.0));
        let (ta, tt) = (tight[0].as_ref().unwrap(), tight[1].as_ref().unwrap());
        assert_eq!(ta.marks[0].diagnostic, Some(false));
        assert_eq!((ta.unique_hit, tt.unique_hit), (1, 0));
        assert!(ta.unique_miss >= 4);
    }
}
