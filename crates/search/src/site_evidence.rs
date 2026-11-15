//! Evidence on which modification a residue carries, taken from the spectrum
//! while its peaks are at hand.
//!
//! Near-isobaric modifications (acetyl and trimethyl differ by 0.036 Da) or
//! modifications that sum to another (two methyls and one dimethyl) give PSMs
//! whose scores barely differ. Three kinds of evidence separate them:
//! - the fragment ions of a form that no other form of the same peptide among
//!   the spectrum's candidates has (the site-determining ions);
//! - the diagnostic ions a modification declares (Mods.txt `diagnostic=`),
//!   such as the immonium ion of the modified residue;
//! - the neutral losses it declares (Mods.txt `loss=`), from the precursor or
//!   from the fragments that hold the modified residue.
//!
//! Rows are compared after rescoring (see `output::site_identity`); nothing
//! here changes a score.

use model::mass::{H2O, PROTON};
use model::modification::Modification;
use model::peptide::Peptide;
use model::tolerance::Tolerance;

/// One declared signature of a row's modification and whether the spectrum
/// shows it.
#[derive(Debug, Clone)]
pub struct ModMark {
    pub modification: &'static Modification,
    /// One of its diagnostic ions is present; `None` when it declares none.
    pub diagnostic: Option<bool>,
    /// One of its neutral losses is present; `None` when it declares none.
    pub loss: Option<bool>,
}

/// A row's evidence. Every row of a spectrum gets one when any of its rows
/// carries a variable modification, so the unmodified forms compete too.
#[derive(Debug, Clone, Default)]
pub struct SiteEvidence {
    /// Fragment ions of this form that no other form of the same peptide
    /// among the spectrum's rows has, matched and not matched.
    pub unique_hit: u16,
    pub unique_miss: u16,
    /// Share of this form's fragment ions that are matched.
    pub hit_fraction: f32,
    /// One mark per carried modification that declares diagnostic ions or
    /// neutral losses, in residue order.
    pub marks: Vec<ModMark>,
}

/// b and y ions at charge 1, and at charge 2 for a precursor of charge 3 or
/// more.
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

/// Whether a peak lies within `tol` of `mz`. `peaks` is sorted by m/z.
fn present(peaks: &[(f64, f32)], mz: f64, tol: Tolerance) -> bool {
    let tol = tol.as_da(mz);
    let i = peaks.partition_point(|p| p.0 < mz - tol);
    i < peaks.len() && peaks[i].0 <= mz + tol
}

/// Whether `ions` (sorted) holds a value within `tol` Da of `mz`.
fn near(ions: &[f64], mz: f64, tol: f64) -> bool {
    let i = ions.partition_point(|&v| v < mz - tol);
    i < ions.len() && ions[i] <= mz + tol
}

/// One candidate row of a spectrum.
pub struct Row<'a> {
    pub peptide: &'a Peptide,
    pub is_decoy: bool,
    pub charge: u8,
    /// Neutral precursor mass the row was matched at.
    pub precursor_mass: f64,
}

/// The evidence of every row of one spectrum; all `None` when no row carries
/// a variable modification. `tol` is the run's fragment tolerance (ppm on
/// high-resolution data), so near-isobaric ions such as the acetyl-lysine
/// immonium ion (126.091), TMT126 (126.128) and the HexNAc oxonium fragment
/// (126.055) are told apart; an ion is unique when no other form has one
/// within twice that tolerance.
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

/// A neutral loss from the precursor (any charge up to the row's) or from a
/// singly charged b or y ion that holds residue `at`.
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

    fn line(s: &str) -> &'static Modification {
        leak_mod(Modification::from_mods_txt_line(s).unwrap())
    }

    #[test]
    fn diagnostic_ions_and_losses_mark_their_modification() {
        let ac = line("42.010565,K,opt,any,Acetyl,diagnostic=143.1179");
        let me3 = line("42.04695,K,opt,any,Trimethyl,diagnostic=143.1543,loss=59.0735");
        let a = peptide(b"PEPKTIDER", 3, ac);
        let t = peptide(b"PEPKTIDER", 3, me3);
        let mass = t.mass();
        // The acetyl immonium ion and a b4 ion of the acetyl form, nothing of
        // trimethyl (no 143.1543, no loss).
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
        // b4..b8 differ by 0.036 Da (> 2 × 0.01): unique to each form; b4 of
        // the acetyl form is matched.
        assert_eq!((ea.unique_hit, et.unique_hit), (1, 0));
        assert!(ea.unique_miss >= 4 && et.unique_miss >= 5);

        // The precursor minus trimethylamine marks the loss.
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

    /// At 0.5 Da the TMT126 reporter (126.128) counts as the acetyl-lysine
    /// immonium ion (126.091) and the near-isobaric acetyl / trimethyl forms
    /// share every ion; at 20 ppm neither happens.
    #[test]
    fn a_ppm_tolerance_separates_near_isobaric_ions() {
        let ac = line("42.010565,K,opt,any,Acetyl,diagnostic=126.0913");
        let me3 = line("42.04695,K,opt,any,Trimethyl");
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
