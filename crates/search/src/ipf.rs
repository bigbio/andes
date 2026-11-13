//! Site localization by inference of peptidoforms (IPF).
//!
//! Port of the hierarchical Bayesian model of OpenSWATH/PyProphet IPF
//! (Rosenberger et al., Nat. Biotechnol. 2017; `pyprophet/ipf.py`) to DDA
//! PSMs. For one PSM the hypotheses are every placement of its variable
//! modifications (the peptidoforms) plus h0, "the PSM is wrong":
//! - prior: h0 = PSM PEP; each peptidoform = (1 - PSM PEP) / #peptidoforms;
//! - evidence: each site-determining ion (a fragment that not every
//!   peptidoform shares) has a PEP. It is evidence *for* the hypotheses it
//!   maps to (1 - PEP) and *against* every other hypothesis, h0 included
//!   (PEP), exactly as PyProphet's transition bitmask;
//! - posterior: prior x product of evidence, normalised over the PSM's
//!   hypotheses (in log space, as PyProphet's `use_log_space=True`).
//!
//! The peptidoform PEP is 1 - posterior and [`model_fdr`] turns PEPs into
//! q-values as PyProphet's `compute_model_fdr`. A residue's localization
//! probability is the posterior mass of the peptidoforms that modify it.
//!
//! This module is pure: enumeration, ion bookkeeping and the Bayesian step.
//! Matching ions against a spectrum and learning ion PEPs live elsewhere.

use model::mass::{H2O, PROTON};
use model::modification::{ModLocation, Modification};
use model::peptide::Peptide;

/// Hard cap on enumerated placements per PSM; [`placements`] reports when it
/// was hit so callers can flag the PSM instead of silently truncating.
pub const MAX_PLACEMENTS: usize = 1024;

/// Every placement of `peptide`'s variable modifications on the residues each
/// one may occupy. Fixed modifications stay where they are. The input
/// placement is always included. Returns the placements (each a full peptide)
/// and whether [`MAX_PLACEMENTS`] cut the enumeration short.
///
/// Only residue-anywhere variable mods move: terminal mods (protein or peptide
/// N/C-term) have a single possible site and are kept in place.
pub fn placements(peptide: &Peptide) -> (Vec<Peptide>, bool) {
    let n = peptide.residues.len();
    // Movable mods, one entry per occurrence, in residue order.
    let movable: Vec<&'static Modification> = peptide
        .residues
        .iter()
        .filter_map(|a| a.mod_.filter(|m| is_movable(m)))
        .collect();
    if movable.is_empty() {
        return (vec![peptide.clone()], false);
    }
    // Residues free to receive a movable mod: unmodified, or carrying a
    // movable mod (which is lifted off below). Fixed and terminal mods stay.
    let base: Vec<_> = peptide
        .residues
        .iter()
        .map(|a| {
            let mut a = a.clone();
            if a.mod_.is_some_and(is_movable) {
                a.mod_ = None;
            }
            a
        })
        .collect();
    let free: Vec<usize> = (0..n).filter(|&i| base[i].mod_.is_none()).collect();
    // Candidate sites per mod occurrence.
    let sites: Vec<Vec<usize>> = movable
        .iter()
        .map(|m| {
            free.iter()
                .copied()
                .filter(|&i| m.applies_to(base[i].residue, ModLocation::Anywhere))
                .collect()
        })
        .collect();
    let mut e = Enum {
        movable: &movable,
        sites: &sites,
        base: &base,
        peptide,
        chosen: vec![usize::MAX; movable.len()],
        seen: std::collections::HashSet::new(),
        out: Vec::new(),
        truncated: false,
    };
    e.rec(0);
    let (mut out, truncated) = (e.out, e.truncated);
    // Guarantee the observed placement is present (it can be cut by the cap).
    if !out.iter().any(|p| same_form(p, peptide)) {
        out.push(peptide.clone());
    }
    (out, truncated)
}

/// Depth-first enumeration state for [`placements`].
struct Enum<'a> {
    movable: &'a [&'static Modification],
    sites: &'a [Vec<usize>],
    base: &'a [model::amino_acid::AminoAcid],
    peptide: &'a Peptide,
    chosen: Vec<usize>,
    seen: std::collections::HashSet<Vec<(usize, u64)>>,
    out: Vec<Peptide>,
    truncated: bool,
}

impl Enum<'_> {
    fn rec(&mut self, k: usize) {
        if self.out.len() >= MAX_PLACEMENTS {
            self.truncated = true;
            return;
        }
        if k == self.movable.len() {
            // Canonical key: sorted (site, mod mass bits) so permutations of
            // identical mods collapse to one peptidoform.
            let mut key: Vec<(usize, u64)> = self
                .chosen
                .iter()
                .zip(self.movable)
                .map(|(&s, m)| (s, m.mass_delta.to_bits()))
                .collect();
            key.sort_unstable();
            if self.seen.insert(key) {
                let mut res = self.base.to_vec();
                for (&s, m) in self.chosen.iter().zip(self.movable) {
                    res[s].mod_ = Some(m);
                }
                let mut p = Peptide::new(res, self.peptide.pre, self.peptide.post);
                p.charge = self.peptide.charge;
                self.out.push(p);
            }
            return;
        }
        for &s in &self.sites[k] {
            if self.chosen[..k].contains(&s) {
                continue;
            }
            self.chosen[k] = s;
            self.rec(k + 1);
        }
        self.chosen[k] = usize::MAX;
    }
}

fn is_movable(m: &Modification) -> bool {
    !m.fixed && m.location == ModLocation::Anywhere
}

/// Same residues and the same modification (by mass) on each.
pub fn same_form(a: &Peptide, b: &Peptide) -> bool {
    a.residues.len() == b.residues.len()
        && a.residues.iter().zip(&b.residues).all(|(x, y)| {
            x.residue == y.residue
                && x.mod_.map(|m| m.mass_delta.to_bits()) == y.mod_.map(|m| m.mass_delta.to_bits())
        })
}

/// 0-based residue indices carrying a movable (variable, anywhere) mod.
pub fn modified_sites(p: &Peptide) -> Vec<usize> {
    p.residues
        .iter()
        .enumerate()
        .filter(|(_, a)| a.mod_.is_some_and(is_movable))
        .map(|(i, _)| i)
        .collect()
}

/// A fragment ion and the peptidoforms that produce it.
#[derive(Debug, Clone, PartialEq)]
pub struct SiteIon {
    pub mz: f64,
    /// `b` or `y`.
    pub kind: u8,
    /// Ion number (residues covered).
    pub number: u16,
    pub charge: u8,
    /// Indices into the placement list of the forms with an ion at `mz`.
    pub forms: Vec<u32>,
}

/// b/y ions of one form: charge 1, plus charge 2 when the precursor is 3+.
fn by_ions(p: &Peptide, precursor_charge: u8) -> Vec<(f64, u8, u16, u8)> {
    let masses: Vec<f64> = p
        .residues
        .iter()
        .map(|a| a.mass + a.mod_.map_or(0.0, |m| m.mass_delta))
        .collect();
    let n = masses.len();
    let total: f64 = masses.iter().sum();
    let mut out = Vec::with_capacity(4 * n);
    let mut prefix = 0.0;
    for (i, m) in masses[..n.saturating_sub(1)].iter().enumerate() {
        prefix += m;
        let b = prefix + PROTON;
        let y = total - prefix + H2O + PROTON;
        let (bn, yn) = ((i + 1) as u16, (n - i - 1) as u16);
        out.push((b, b'b', bn, 1));
        out.push((y, b'y', yn, 1));
        if precursor_charge >= 3 {
            out.push(((b + PROTON) / 2.0, b'b', bn, 2));
            out.push(((y + PROTON) / 2.0, b'y', yn, 2));
        }
    }
    out
}

/// The site-determining ions of a set of forms: every b/y ion (merged across
/// forms within `tol_da`) that is NOT produced by all forms. Ions every form
/// shares carry no localization information and are dropped, as IPF only
/// scores transitions that discriminate peptidoforms.
pub fn site_ions(
    forms: &[Peptide],
    precursor_charge: u8,
    tol_da: impl Fn(f64) -> f64,
) -> Vec<SiteIon> {
    let mut all: Vec<(f64, u8, u16, u8, u32)> = forms
        .iter()
        .enumerate()
        .flat_map(|(f, p)| {
            by_ions(p, precursor_charge)
                .into_iter()
                .map(move |(mz, k, num, z)| (mz, k, num, z, f as u32))
        })
        .collect();
    all.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut out: Vec<SiteIon> = Vec::new();
    for (mz, kind, number, charge, f) in all {
        match out.last_mut() {
            Some(last) if mz - last.mz <= tol_da(mz) => {
                if !last.forms.contains(&f) {
                    last.forms.push(f);
                }
            }
            _ => out.push(SiteIon {
                mz,
                kind,
                number,
                charge,
                forms: vec![f],
            }),
        }
    }
    let n = forms.len();
    out.retain(|ion| ion.forms.len() < n);
    for ion in &mut out {
        ion.forms.sort_unstable();
    }
    out
}

/// One piece of evidence: an ion's PEP and the hypotheses it supports.
#[derive(Debug, Clone)]
pub struct Evidence {
    pub pep: f64,
    /// Indices of the peptidoforms this ion maps to (bitmask = 1).
    pub forms: Vec<u32>,
}

/// Posterior over one PSM's hypotheses.
#[derive(Debug, Clone, PartialEq)]
pub struct Posterior {
    /// Per peptidoform, in placement order.
    pub forms: Vec<f64>,
    /// h0: the PSM itself is wrong.
    pub h0: f64,
}

/// The IPF Bayesian step for one PSM (PyProphet `prepare_transition_bm` +
/// `apply_bm`, log space). `psm_pep` is the PSM's PEP after rescoring,
/// `n_forms` the number of peptidoforms. `epsilon` clips evidence into
/// `[epsilon, 1 - epsilon]` (PyProphet `bm_evidence_epsilon`; 0 disables).
pub fn posterior(psm_pep: f64, n_forms: usize, evidence: &[Evidence], epsilon: f64) -> Posterior {
    assert!(n_forms > 0, "a PSM has at least its own peptidoform");
    let clip = |e: f64| {
        if epsilon > 0.0 {
            e.clamp(epsilon, 1.0 - epsilon)
        } else {
            e
        }
    };
    let ln = |x: f64| if x > 0.0 { x.ln() } else { f64::NEG_INFINITY };
    let psm_pep = psm_pep.clamp(0.0, 1.0);
    // log prior + sum of log evidence, per hypothesis (index n_forms = h0).
    let mut l: Vec<f64> = vec![ln((1.0 - psm_pep) / n_forms as f64); n_forms];
    l.push(ln(psm_pep));
    for e in evidence {
        let (for_, against) = (ln(clip(1.0 - e.pep)), ln(clip(e.pep)));
        for (h, lh) in l.iter_mut().enumerate() {
            let maps = h < n_forms && e.forms.contains(&(h as u32));
            *lh += if maps { for_ } else { against };
        }
    }
    let mx = l.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !mx.is_finite() {
        // Every hypothesis impossible (all-zero evidence): PyProphet's
        // fillna(0) leaves every posterior at 0.
        return Posterior {
            forms: vec![0.0; n_forms],
            h0: 0.0,
        };
    }
    let z: f64 = l.iter().map(|x| (x - mx).exp()).sum();
    let p: Vec<f64> = l.iter().map(|x| (x - mx).exp() / z).collect();
    Posterior {
        forms: p[..n_forms].to_vec(),
        h0: p[n_forms],
    }
}

/// Model-based q-values from PEPs (PyProphet `compute_model_fdr`): sort
/// ascending, q = cumulative mean of PEP, ties take the value at their
/// highest rank.
pub fn model_fdr(peps: &[f64]) -> Vec<f64> {
    let mut order: Vec<usize> = (0..peps.len()).collect();
    order.sort_by(|&a, &b| peps[a].total_cmp(&peps[b]));
    let mut cum = Vec::with_capacity(peps.len());
    let mut s = 0.0;
    for &i in &order {
        s += peps[i];
        cum.push(s);
    }
    let mut q = vec![0.0; peps.len()];
    let mut start = 0;
    while start < order.len() {
        let mut end = start;
        while end + 1 < order.len() && peps[order[end + 1]] == peps[order[start]] {
            end += 1;
        }
        let v = cum[end] / (end + 1) as f64;
        for &i in &order[start..=end] {
            q[i] = v;
        }
        start = end + 1;
    }
    q
}

/// Localization probability of every residue of the PSM's peptide: the
/// posterior mass of the peptidoforms that put a movable mod there,
/// normalised by the peptidoform mass (so it is conditional on the PSM
/// being right, as a site probability is usually reported).
pub fn site_probabilities(forms: &[Peptide], post: &Posterior) -> Vec<f64> {
    let n = forms.first().map_or(0, |p| p.residues.len());
    let mut out = vec![0.0; n];
    let total: f64 = post.forms.iter().sum();
    if total <= 0.0 {
        return out;
    }
    for (p, &w) in forms.iter().zip(&post.forms) {
        for i in modified_sites(p) {
            out[i] += w / total;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::amino_acid::AminoAcid;
    use model::modification::ResidueSpec;

    fn phospho(r: u8) -> &'static Modification {
        Box::leak(Box::new(Modification {
            name: "Phospho".into(),
            mass_delta: 79.966331,
            residue: ResidueSpec::Specific(r),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: vec![],
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            diagnostic_ions: vec![],
        }))
    }

    /// Peptide from a sequence, phosphorylating the 0-based indices given.
    fn pep(seq: &str, ph: &[usize]) -> Peptide {
        let res = seq
            .bytes()
            .enumerate()
            .map(|(i, r)| {
                let mut a = AminoAcid::standard(r).unwrap();
                if ph.contains(&i) {
                    a.mod_ = Some(phospho(r));
                }
                a
            })
            .collect();
        Peptide::new(res, b'K', b'A')
    }

    #[test]
    fn placements_enumerate_every_site_once() {
        // One phospho, sites S at 1 and 4 (T at 3 needs a T-phospho mod).
        let (forms, cut) = placements(&pep("ASATSK", &[1]));
        assert!(!cut);
        let mut sites: Vec<Vec<usize>> = forms.iter().map(modified_sites).collect();
        sites.sort();
        assert_eq!(sites, vec![vec![1], vec![4]]);
        // Two identical mods on three S: C(3,2) = 3 forms, no duplicates.
        let (forms, _) = placements(&pep("SASAS", &[0, 2]));
        assert_eq!(forms.len(), 3);
        // Unmodified peptide: itself only.
        assert_eq!(placements(&pep("PEPTIDE", &[])).0.len(), 1);
    }

    #[test]
    fn site_ions_drop_shared_ions() {
        let (forms, _) = placements(&pep("ASATSK", &[1]));
        let ions = site_ions(&forms, 2, |_| 0.01);
        assert!(!ions.is_empty());
        // Every kept ion separates the two forms.
        assert!(ions.iter().all(|i| i.forms.len() == 1));
        // b1 (A) is shared by both forms and must be gone.
        assert!(!ions.iter().any(|i| i.kind == b'b' && i.number == 1));
    }

    /// Reference computed with PyProphet's apply_bm (linear space) by hand:
    /// psm_pep 0.1, two forms, one ion for form 0 with pep 0.2.
    #[test]
    fn posterior_matches_pyprophet_arithmetic() {
        let e = [Evidence {
            pep: 0.2,
            forms: vec![0],
        }];
        let p = posterior(0.1, 2, &e, 0.0);
        // likelihood*prior: f0 = 0.45*0.8, f1 = 0.45*0.2, h0 = 0.1*0.2
        let (f0, f1, h0) = (0.45 * 0.8, 0.45 * 0.2, 0.1 * 0.2);
        let z = f0 + f1 + h0;
        assert!((p.forms[0] - f0 / z).abs() < 1e-12);
        assert!((p.forms[1] - f1 / z).abs() < 1e-12);
        assert!((p.h0 - h0 / z).abs() < 1e-12);
    }

    #[test]
    fn no_evidence_gives_the_prior() {
        let p = posterior(0.2, 4, &[], 0.0);
        for f in &p.forms {
            assert!((f - 0.2).abs() < 1e-12);
        }
        assert!((p.h0 - 0.2).abs() < 1e-12);
    }

    #[test]
    fn model_fdr_matches_pyprophet() {
        // compute_model_fdr([0.1, 0.3, 0.1, 0.5]): sorted 0.1,0.1,0.3,0.5;
        // ties at max rank 2 -> 0.2/2 = 0.1; 0.5/3; 1.0/4.
        let q = model_fdr(&[0.1, 0.3, 0.1, 0.5]);
        let want = [0.1, 0.5 / 3.0, 0.1, 0.25];
        for (a, b) in q.iter().zip(want) {
            assert!((a - b).abs() < 1e-12, "{q:?}");
        }
    }

    #[test]
    fn site_probabilities_sum_over_forms() {
        let (forms, _) = placements(&pep("ASATSK", &[1]));
        let post = Posterior {
            forms: vec![0.6, 0.3],
            h0: 0.1,
        };
        let s = site_probabilities(&forms, &post);
        let i1 = modified_sites(&forms[0])[0];
        let i2 = modified_sites(&forms[1])[0];
        assert!((s[i1] - 2.0 / 3.0).abs() < 1e-12);
        assert!((s[i2] - 1.0 / 3.0).abs() < 1e-12);
    }
}
