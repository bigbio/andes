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

/// Shifts of the decoy ions, in Da. Each site ion gets one decoy per shift,
/// at an m/z no b/y ion of the peptide sits at (off the 1.00 Da isotope grid
/// and its halves), so a decoy match is a chance match. Their match rate per
/// intensity tier is the null the ion PEPs are learned against.
pub const DECOY_SHIFTS: [f64; 4] = [-7.37, -3.71, 3.83, 7.59];

/// One site ion, judged against the spectrum.
#[derive(Debug, Clone, PartialEq)]
pub struct IonMatch {
    /// Indices into the placement list of the forms with this ion.
    pub forms: Vec<u32>,
    /// Intensity tier of the matched peak (0 = most intense), `None` if no
    /// peak lies within tolerance.
    pub tier: Option<u8>,
}

/// IPF evidence of one PSM: its peptidoforms (as placement order),
/// the site ions matched against the spectrum and the decoy-ion tiers.
#[derive(Debug, Clone, Default)]
pub struct IpfEvidence {
    /// 0-based modified residue indices of each placement.
    pub forms: Vec<Vec<u16>>,
    /// Index of the PSM's own placement in `forms`.
    pub observed: u32,
    /// Whether [`MAX_PLACEMENTS`] cut the enumeration.
    pub truncated: bool,
    pub ions: Vec<IonMatch>,
    /// Tier of every decoy ion (`None` = unmatched), pooled over shifts.
    pub decoys: Vec<Option<u8>>,
}

/// Number of intensity tiers; peaks are tiered by rank within the spectrum.
pub const TIERS: u8 = 4;

/// Intensity tier of each peak: rank-quartiles of the spectrum (0 = top).
fn peak_tiers(peaks: &[(f64, f32)]) -> Vec<u8> {
    let mut order: Vec<usize> = (0..peaks.len()).collect();
    order.sort_by(|&a, &b| peaks[b].1.total_cmp(&peaks[a].1));
    let mut tier = vec![0u8; peaks.len()];
    let n = peaks.len().max(1);
    for (rank, &i) in order.iter().enumerate() {
        tier[i] = ((rank * TIERS as usize) / n) as u8;
    }
    tier
}

/// Best (lowest) tier of a peak within `tol` Da of `mz`. `peaks` sorted by m/z.
fn match_tier(peaks: &[(f64, f32)], tiers: &[u8], mz: f64, tol: f64) -> Option<u8> {
    let i = peaks.partition_point(|p| p.0 < mz - tol);
    peaks[i..]
        .iter()
        .zip(&tiers[i..])
        .take_while(|(p, _)| p.0 <= mz + tol)
        .map(|(_, &t)| t)
        .min()
}

/// IPF evidence of one PSM's peptide against its spectrum. `None` when the
/// peptide has no movable modification or only one possible placement.
pub fn psm_evidence(
    peaks: &[(f64, f32)],
    peptide: &Peptide,
    precursor_charge: u8,
    tol_da: impl Fn(f64) -> f64 + Copy,
) -> Option<IpfEvidence> {
    let (forms, truncated) = placements(peptide);
    if forms.len() < 2 {
        return None;
    }
    let observed = forms.iter().position(|p| same_form(p, peptide))? as u32;
    let tiers = peak_tiers(peaks);
    let ions = site_ions(&forms, precursor_charge, tol_da);
    let matches = ions
        .iter()
        .map(|ion| IonMatch {
            forms: ion.forms.clone(),
            tier: match_tier(peaks, &tiers, ion.mz, tol_da(ion.mz)),
        })
        .collect();
    let decoys = ions
        .iter()
        .flat_map(|ion| {
            DECOY_SHIFTS
                .iter()
                .map(move |s| (ion.mz + s / ion.charge as f64).max(0.0))
        })
        .map(|mz| match_tier(peaks, &tiers, mz, tol_da(mz)))
        .collect();
    Some(IpfEvidence {
        forms: forms
            .iter()
            .map(|p| modified_sites(p).into_iter().map(|i| i as u16).collect())
            .collect(),
        observed,
        truncated,
        ions: matches,
        decoys,
    })
}

/// PEP of a site-ion observation by intensity tier, learned from one run.
///
/// Decoy ions only match by chance, so their tier distribution (with
/// "unmatched" as a tier of its own) is the null. A target site ion's tier
/// distribution is a mixture of the null and true matches. With pi0 the
/// share of null target ions (estimated from the "unmatched" bin, where true
/// ions are rare, and bounded to [0, 1]) the local FDR of each tier is
/// `pi0 * P(tier | null) / P(tier | target)` (the two-groups model behind
/// PyProphet's PEPs). An unmatched ion is scored like an OpenSWATH transition
/// without signal: a high PEP on the forms that have it, so it counts against
/// them and for h0. Its PEP is 1 - r, with r the share of unmatched target
/// site ions that are not chance (bounded to [0.01, 0.5] so one missing ion
/// never rules a form out).
#[derive(Debug, Clone, PartialEq)]
pub struct IonPepModel {
    /// PEP of a matched ion, per tier.
    pub matched: [f64; TIERS as usize],
    /// PEP of an unmatched site ion (applied to the forms that have it).
    pub unmatched: f64,
    pub pi0: f64,
}

impl IonPepModel {
    /// Fit on every PSM's evidence of the run (targets above a PSM-level
    /// confidence filter are the usual input).
    pub fn fit<'a>(evidence: impl IntoIterator<Item = &'a IpfEvidence>) -> Self {
        let bins = TIERS as usize + 1; // last bin = unmatched
        let (mut t, mut d) = (vec![1.0f64; bins], vec![1.0f64; bins]); // +1 pseudo-counts
        for e in evidence {
            for ion in &e.ions {
                t[ion.tier.map_or(bins - 1, |x| x as usize)] += 1.0;
            }
            for tier in &e.decoys {
                d[tier.map_or(bins - 1, |x| x as usize)] += 1.0;
            }
        }
        let (st, sd): (f64, f64) = (t.iter().sum(), d.iter().sum());
        let pt: Vec<f64> = t.iter().map(|x| x / st).collect();
        let pd: Vec<f64> = d.iter().map(|x| x / sd).collect();
        // Storey-style pi0 from the unmatched bin, where true ions are rare.
        let pi0 = (pt[bins - 1] / pd[bins - 1]).clamp(0.0, 1.0);
        let mut matched = [0.0; TIERS as usize];
        for (k, m) in matched.iter_mut().enumerate() {
            *m = (pi0 * pd[k] / pt[k]).clamp(1e-6, 1.0);
        }
        // Monotone in tier: a weaker peak is never better evidence.
        for k in 1..matched.len() {
            matched[k] = matched[k].max(matched[k - 1]);
        }
        let r = (1.0 - pi0 * pd[bins - 1] / pt[bins - 1]).clamp(0.01, 0.5);
        IonPepModel {
            matched,
            unmatched: 1.0 - r,
            pi0,
        }
    }

    /// IPF evidence rows of one PSM, one per site ion: its tier PEP when
    /// matched, [`Self::unmatched`] when not, always on the forms that have it
    /// (PyProphet's transition bitmask).
    pub fn evidence(&self, e: &IpfEvidence) -> Vec<Evidence> {
        e.ions
            .iter()
            .map(|ion| Evidence {
                pep: ion
                    .tier
                    .map_or(self.unmatched, |t| self.matched[t as usize]),
                forms: ion.forms.clone(),
            })
            .collect()
    }
}

/// Localization result of one PSM.
#[derive(Debug, Clone, PartialEq)]
pub struct Localization {
    /// Index (into `IpfEvidence::forms`) of the most probable placement.
    pub best: u32,
    /// Posterior of every placement (h0 excluded), in placement order.
    pub forms: Vec<f64>,
    pub h0: f64,
    /// PEP of the best placement: 1 - its posterior (PyProphet's peptidoform PEP).
    pub pep: f64,
    /// Localization probability of the best placement given the PSM is right.
    pub best_prob: f64,
}

/// IPF for one PSM: prior from the PSM's PEP, evidence from the ion model.
pub fn localize(psm_pep: f64, e: &IpfEvidence, model: &IonPepModel) -> Localization {
    let post = posterior(psm_pep, e.forms.len(), &model.evidence(e), 0.0);
    let best = (0..post.forms.len())
        .max_by(|&a, &b| post.forms[a].total_cmp(&post.forms[b]).then(b.cmp(&a)))
        .unwrap_or(0);
    let total: f64 = post.forms.iter().sum();
    Localization {
        best: best as u32,
        pep: 1.0 - post.forms[best],
        best_prob: if total > 0.0 {
            post.forms[best] / total
        } else {
            0.0
        },
        h0: post.h0,
        forms: post.forms,
    }
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

    #[test]
    fn localize_recovers_the_true_site_from_a_simulated_spectrum() {
        // Truth: phospho on S4 of ASATSKSAPEK-like peptide; the PSM reported S1.
        let truth = pep("ASATSKSAPEK", &[4]);
        let reported = pep("ASATSKSAPEK", &[1]);
        // Spectrum: b/y of the truth at charge 1 (intense) + uniform noise.
        let mut peaks: Vec<(f64, f32)> = by_ions(&truth, 2)
            .into_iter()
            .map(|(mz, ..)| (mz, 100.0))
            .collect();
        let mut x = 0.123_f64;
        for _ in 0..200 {
            x = (x * 9301.0 + 49297.0) % 233280.0; // LCG
            peaks.push((150.0 + 1200.0 * x / 233280.0, 5.0));
        }
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        let tol = |mz: f64| mz * 10e-6;
        let e = psm_evidence(&peaks, &reported, 2, tol).unwrap();
        assert_eq!(e.forms.len(), 3); // S1, S4, S6
        assert_eq!(e.forms[e.observed as usize], vec![1]);
        // Fit on many copies (one run's worth of PSMs).
        let model = IonPepModel::fit(std::iter::repeat_n(&e, 50));
        assert!(model.matched[0] < 0.2, "{model:?}");
        let loc = localize(0.01, &e, &model);
        assert_eq!(e.forms[loc.best as usize], vec![4], "{loc:?}");
        assert!(loc.best_prob > 0.9, "{loc:?}");
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
