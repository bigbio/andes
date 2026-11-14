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

/// Hard cap on enumerated placements per PSM; [`placements_with`] reports when it
/// was hit so callers can flag the PSM instead of silently truncating.
pub const MAX_PLACEMENTS: usize = 1024;

/// Every placement of `peptide`'s variable modifications on the residues each
/// one may occupy. Fixed modifications stay where they are. The input
/// placement is always included. Returns the placements (each a full peptide)
/// and whether [`MAX_PLACEMENTS`] cut the enumeration short.
///
/// Only residue-anywhere variable mods move: terminal mods (protein or peptide
/// N/C-term) have a single possible site and are kept in place. A mod also
/// moves to residues that a `declared` sibling record (same name and mass,
/// e.g. Mods.txt's separate Phospho S, T and Y lines) accepts; the residue then
/// carries that sibling record.
pub fn placements_with(
    peptide: &Peptide,
    declared: &[&'static Modification],
) -> (Vec<Peptide>, bool) {
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
    let sites: Vec<Vec<(usize, &'static Modification)>> = movable
        .iter()
        .map(|&m| {
            free.iter()
                .filter_map(|&i| {
                    let r = base[i].residue;
                    if m.applies_to(r, ModLocation::Anywhere) {
                        return Some((i, m));
                    }
                    declared
                        .iter()
                        .copied()
                        .find(|d| {
                            is_movable(d)
                                && d.name == m.name
                                && (d.mass_delta - m.mass_delta).abs() < 1e-6
                                && d.applies_to(r, ModLocation::Anywhere)
                        })
                        .map(|d| (i, d))
                })
                .collect()
        })
        .collect();
    // Consecutive mods with the same candidates are interchangeable: visit
    // their sites in increasing order only, so each combination is built once
    // (n identical mods would otherwise be enumerated n! times).
    let same_sites = |a: &[(usize, &'static Modification)],
                      b: &[(usize, &'static Modification)]| {
        a.len() == b.len()
            && a.iter()
                .zip(b)
                .all(|(x, y)| x.0 == y.0 && std::ptr::eq(x.1, y.1))
    };
    let same_as_prev: Vec<bool> = (0..sites.len())
        .map(|k| k > 0 && same_sites(&sites[k], &sites[k - 1]))
        .collect();
    let mut e = Enum {
        sites: &sites,
        same_as_prev: &same_as_prev,
        base: &base,
        peptide,
        chosen: movable.iter().map(|&m| (usize::MAX, m)).collect(),
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

/// Depth-first enumeration state for [`placements_with`].
struct Enum<'a> {
    sites: &'a [Vec<(usize, &'static Modification)>],
    same_as_prev: &'a [bool],
    base: &'a [model::amino_acid::AminoAcid],
    peptide: &'a Peptide,
    chosen: Vec<(usize, &'static Modification)>,
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
        if k == self.sites.len() {
            // Canonical key: sorted (site, mod mass bits) so permutations of
            // identical mods collapse to one peptidoform.
            let mut key: Vec<(usize, u64)> = self
                .chosen
                .iter()
                .map(|&(s, m)| (s, m.mass_delta.to_bits()))
                .collect();
            key.sort_unstable();
            if self.seen.insert(key) {
                let mut res = self.base.to_vec();
                for &(s, m) in &self.chosen {
                    res[s].mod_ = Some(m);
                }
                let mut p = Peptide::new(res, self.peptide.pre, self.peptide.post);
                p.charge = self.peptide.charge;
                self.out.push(p);
            }
            return;
        }
        for &(s, m) in &self.sites[k] {
            if self.chosen[..k].iter().any(|c| c.0 == s)
                || (self.same_as_prev[k] && s < self.chosen[k - 1].0)
            {
                continue;
            }
            self.chosen[k] = (s, m);
            self.rec(k + 1);
        }
        self.chosen[k].0 = usize::MAX;
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
    /// Indices into the placement list of the forms with an ion at `mz`.
    pub forms: Vec<u32>,
}

/// b/y ions of one form: charge 1, plus charge 2 when the precursor is 3+.
/// A fragment that holds a residue whose modification declares neutral losses
/// (Mods.txt `loss=`, e.g. phospho's 97.977 H3PO4) also gives one ion per
/// distinct loss at the same charges: the loss ion pins the mod to the fragment
/// as well as the intact one does. Calls `emit(mz, kind, number)`.
fn by_ions(p: &Peptide, precursor_charge: u8, mut emit: impl FnMut(f64, u8, u16)) {
    let n = p.residues.len();
    if n < 2 {
        return;
    }
    let losses_at = |i: usize| -> &[f64] {
        p.residues[i]
            .mod_
            .filter(|m| !m.fixed)
            .map_or(&[][..], |m| &m.neutral_losses[..])
    };
    // Distinct loss masses carried by residues 0..=i (b) and i+1..n (y).
    fn add(set: &mut Vec<f64>, losses: &[f64]) {
        for &l in losses {
            if !set.iter().any(|x| (x - l).abs() < 1e-6) {
                set.push(l);
            }
        }
    }
    let mut y_losses: Vec<Vec<f64>> = vec![Vec::new(); n];
    for i in (0..n - 1).rev() {
        let mut s = y_losses[i + 1].clone();
        add(&mut s, losses_at(i + 1));
        y_losses[i] = s;
    }
    let mut emit_all = |m: f64, kind: u8, num: u16, losses: &[f64]| {
        for shift in std::iter::once(0.0).chain(losses.iter().map(|l| -l)) {
            emit(m + shift, kind, num);
            if precursor_charge >= 3 {
                emit((m + shift + PROTON) / 2.0, kind, num);
            }
        }
    };
    let total: f64 = p
        .residues
        .iter()
        .map(|a| a.mass + a.mod_.map_or(0.0, |m| m.mass_delta))
        .sum();
    let mut prefix = 0.0;
    let mut b_losses: Vec<f64> = Vec::new();
    for (i, (a, y_loss)) in p.residues.iter().zip(&y_losses).take(n - 1).enumerate() {
        prefix += a.mass + a.mod_.map_or(0.0, |m| m.mass_delta);
        add(&mut b_losses, losses_at(i));
        emit_all(prefix + PROTON, b'b', (i + 1) as u16, &b_losses);
        emit_all(
            total - prefix + H2O + PROTON,
            b'y',
            (n - i - 1) as u16,
            y_loss,
        );
    }
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
    let mut all: Vec<(f64, u8, u16, u32)> = Vec::new();
    for (f, p) in forms.iter().enumerate() {
        by_ions(p, precursor_charge, |mz, k, num| {
            all.push((mz, k, num, f as u32))
        });
    }
    all.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut out: Vec<SiteIon> = Vec::new();
    for (mz, kind, number, f) in all {
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
                forms: vec![f],
            }),
        }
    }
    let n = forms.len();
    out.retain(|ion| ion.forms.len() < n);
    out
}

/// A set of placement indices, one bit per form (an ion's PyProphet
/// transition bitmask). Compact because a multiply phosphorylated peptide can
/// have up to [`MAX_PLACEMENTS`] forms and each site ion maps to many of them.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct FormSet(Box<[u64]>);

impl FormSet {
    pub fn from_indices(n_forms: usize, idx: &[u32]) -> Self {
        let mut w = vec![0u64; n_forms.div_ceil(64)];
        for &i in idx {
            w[i as usize / 64] |= 1 << (i % 64);
        }
        FormSet(w.into_boxed_slice())
    }
    pub fn contains(&self, i: u32) -> bool {
        self.0
            .get(i as usize / 64)
            .is_some_and(|w| w & (1 << (i % 64)) != 0)
    }
}

/// One site ion, judged against the spectrum.
#[derive(Debug, Clone, PartialEq)]
pub struct IonMatch {
    /// The forms with this ion.
    pub forms: FormSet,
    /// Intensity tier of the matched peak (0 = most intense), `None` if no
    /// peak lies within tolerance.
    pub tier: Option<u8>,
}

/// IPF evidence of one PSM: its peptidoforms (as placement order), the site
/// ions matched against the spectrum, and the search score of each form.
#[derive(Debug, Clone, Default)]
pub struct IpfEvidence {
    /// 0-based modified residue indices of each placement.
    pub forms: Vec<Vec<u16>>,
    /// Index of the PSM's own placement in `forms`.
    pub observed: u32,
    /// Whether [`MAX_PLACEMENTS`] cut the enumeration.
    pub truncated: bool,
    pub ions: Vec<IonMatch>,
    /// The search's own score (node + edge, the rank score without cleavage
    /// credit, which every placement shares) of each placement, in `forms`
    /// order. Empty when not computed.
    pub form_scores: Vec<f32>,
}

/// Number of intensity tiers; peaks are tiered by rank within the spectrum.
pub const TIERS: u8 = 4;

/// Intensity tier of each peak: rank-quartiles of the spectrum (0 = top).
/// Computed once per spectrum and shared by its rows.
pub fn peak_tiers(peaks: &[(f64, f32)]) -> Vec<u8> {
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

/// IPF evidence of one PSM's peptide against its spectrum. `tiers` is
/// [`peak_tiers`] of `peaks`; `score`, when given, is the search's score of a
/// placement (fills [`IpfEvidence::form_scores`]). `None` when the peptide has
/// no movable modification or only one possible placement.
pub fn psm_evidence(
    peaks: &[(f64, f32)],
    tiers: &[u8],
    peptide: &Peptide,
    precursor_charge: u8,
    tol_da: impl Fn(f64) -> f64 + Copy,
    declared: &[&'static Modification],
    score: Option<&dyn Fn(&Peptide) -> f32>,
) -> Option<IpfEvidence> {
    let (forms, truncated) = placements_with(peptide, declared);
    if forms.len() < 2 {
        return None;
    }
    let observed = forms.iter().position(|p| same_form(p, peptide))? as u32;
    // `.iter()`, not `.into_iter()`: an in-place collect would keep the larger
    // SiteIon buffer alive for every stored PSM.
    let ions = site_ions(&forms, precursor_charge, tol_da)
        .iter()
        .map(|ion| IonMatch {
            forms: FormSet::from_indices(forms.len(), &ion.forms),
            tier: match_tier(peaks, tiers, ion.mz, tol_da(ion.mz)),
        })
        .collect();
    Some(IpfEvidence {
        form_scores: score.map_or_else(Vec::new, |s| forms.iter().map(s).collect()),
        forms: forms
            .iter()
            .map(|p| modified_sites(p).into_iter().map(|i| i as u16).collect())
            .collect(),
        observed,
        truncated,
        ions,
    })
}

/// PEP of a site-ion observation by intensity tier, learned from one run.
///
/// The IPF evidence of an ion is (1 - PEP) for the hypotheses that produce it
/// and PEP for the others, so the ratio between the two kinds of hypothesis is
/// a likelihood ratio LR = (1 - PEP) / PEP. Here
/// `LR(tier) = P(tier | the form has the ion) / P(tier | it does not)` and
/// `PEP = 1 / (1 + LR)`. Both distributions come from confident PSMs (PSM
/// q-value filter applied by the caller): the site ions of each PSM's own
/// placement give the first, the ions only other placements produce give the
/// second. "Unmatched" is a tier of its own, so a missing ion counts against
/// the forms that should show it only as much as the run says it should: in
/// DDA most theoretical ions are missing even for the right peptidoform.
#[derive(Debug, Clone, PartialEq)]
pub struct IonPepModel {
    /// PEP of a matched ion, per tier.
    pub matched: [f64; TIERS as usize],
    /// PEP of an unmatched site ion (applied to the forms that have it).
    pub unmatched: f64,
}

impl IonPepModel {
    /// Fit on the evidence of confident PSMs of one run, each labelled by its
    /// own placement ([`IpfEvidence::observed`]).
    pub fn fit<'a>(evidence: impl IntoIterator<Item = &'a IpfEvidence>) -> Self {
        Self::fit_labelled(evidence.into_iter().map(|e| (e, e.observed)))
    }

    /// [`Self::fit`] with each PSM labelled by the given placement index (an
    /// EM step relabels by IPF's own best placement).
    pub fn fit_labelled<'a>(evidence: impl IntoIterator<Item = (&'a IpfEvidence, u32)>) -> Self {
        let bins = TIERS as usize + 1; // last bin = unmatched
        let bin = |t: Option<u8>| t.map_or(bins - 1, |x| x as usize);
        // +1 pseudo-counts so an empty bin never gives LR 0 or infinity.
        let (mut own, mut other) = (vec![1.0f64; bins], vec![1.0f64; bins]);
        for (e, label) in evidence {
            for ion in &e.ions {
                if ion.forms.contains(label) {
                    own[bin(ion.tier)] += 1.0;
                } else {
                    other[bin(ion.tier)] += 1.0;
                }
            }
        }
        let (so, sx): (f64, f64) = (own.iter().sum(), other.iter().sum());
        let pep = |k: usize| {
            let lr = (own[k] / so) / (other[k] / sx);
            (1.0 / (1.0 + lr)).clamp(1e-6, 1.0 - 1e-6)
        };
        let mut matched = [0.0; TIERS as usize];
        for (k, m) in matched.iter_mut().enumerate() {
            *m = pep(k);
        }
        // Monotone in tier: a weaker peak is never better evidence.
        for k in 1..matched.len() {
            matched[k] = matched[k].max(matched[k - 1]);
        }
        IonPepModel {
            matched,
            unmatched: pep(bins - 1),
        }
    }

    /// PEP of one site-ion observation: its tier PEP when matched,
    /// [`Self::unmatched`] when not (applied to the forms that have it).
    pub fn ion_pep(&self, ion: &IonMatch) -> f64 {
        ion.tier
            .map_or(self.unmatched, |t| self.matched[t as usize])
    }
}

/// Localization result of one PSM.
#[derive(Debug, Clone, PartialEq)]
pub struct Localization {
    /// Index (into `IpfEvidence::forms`) of the most probable placement.
    pub best: u32,
    /// Posterior of every placement (h0 excluded), in placement order.
    pub forms: Vec<f64>,
    /// PEP of the best placement: 1 - its posterior (PyProphet's peptidoform PEP).
    pub pep: f64,
    /// Localization probability of the best placement given the PSM is right.
    pub best_prob: f64,
}

/// Weight of the search score in the placement prior: prior ∝ exp(w · score).
/// Fitted by maximum likelihood of the reference placement on half the
/// peptides of each benchmark (synthetic ground truth PXD007058, PXD007653 vs
/// MaxQuant, PXD017476 vs Spectronaut), on a grid of factor-1.33 steps; all
/// three picked this grid point, and fitting on any two picks it for the third.
pub const SEARCH_SCORE_WEIGHT: f64 = 0.237;

/// Prior over a PSM's placements from the search's own score of each
/// (softmax of `SEARCH_SCORE_WEIGHT` · score); uniform without scores.
pub fn placement_prior(e: &IpfEvidence) -> Vec<f64> {
    let n = e.forms.len();
    if e.form_scores.len() != n || e.form_scores.iter().any(|s| !s.is_finite()) {
        return vec![1.0 / n as f64; n];
    }
    let m = e
        .form_scores
        .iter()
        .copied()
        .fold(f32::NEG_INFINITY, f32::max) as f64;
    let ex: Vec<f64> = e
        .form_scores
        .iter()
        .map(|&s| (SEARCH_SCORE_WEIGHT * (s as f64 - m)).exp())
        .collect();
    let sum: f64 = ex.iter().sum();
    ex.into_iter().map(|x| x / sum).collect()
}

/// IPF for one PSM: prior from the PSM's PEP (split over placements by
/// [`placement_prior`]), evidence from the ion model.
pub fn localize(psm_pep: f64, e: &IpfEvidence, model: &IonPepModel) -> Localization {
    let post = log_posterior(
        psm_pep,
        &placement_prior(e),
        e.ions.iter().map(|ion| (model.ion_pep(ion), ion)),
        |ion, h| ion.forms.contains(h),
        0.0,
    );
    let best = (0..post.forms.len())
        // Ties (identical evidence) go to the placement with the lowest
        // residues, so the result does not depend on enumeration order.
        .max_by(|&a, &b| {
            post.forms[a]
                .total_cmp(&post.forms[b])
                .then_with(|| e.forms[b].cmp(&e.forms[a]))
        })
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
        forms: post.forms,
    }
}

/// Posterior over the placements from the ion evidence alone (uniform prior,
/// PSM assumed right): its log is the IPF log-likelihood of each placement up
/// to a constant per PSM. For diagnostics (`--ipf-debug`).
pub fn evidence_posterior(e: &IpfEvidence, model: &IonPepModel) -> Posterior {
    let n = e.forms.len();
    log_posterior(
        0.0,
        &vec![1.0 / n as f64; n],
        e.ions.iter().map(|ion| (model.ion_pep(ion), ion)),
        |ion, h| ion.forms.contains(h),
        0.0,
    )
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

/// The IPF Bayesian step for one PSM with PyProphet's uniform prior
/// (`prepare_transition_bm` + `apply_bm`, log space). `psm_pep` is the PSM's
/// PEP after rescoring, `n_forms` the number of peptidoforms. `epsilon` clips
/// evidence into `[epsilon, 1 - epsilon]` (PyProphet `bm_evidence_epsilon`;
/// 0 disables).
pub fn posterior(psm_pep: f64, n_forms: usize, evidence: &[Evidence], epsilon: f64) -> Posterior {
    assert!(n_forms > 0, "a PSM has at least its own peptidoform");
    log_posterior(
        psm_pep,
        &vec![1.0 / n_forms as f64; n_forms],
        evidence.iter().map(|e| (e.pep, e)),
        |e, h| e.forms.contains(&h),
        epsilon,
    )
}

/// The Bayesian step itself: `prior` over the forms (summing to 1; the forms
/// share `1 - psm_pep` in those proportions, h0 gets `psm_pep`), and per
/// evidence item its PEP and whether it maps to form `h` (`maps`).
fn log_posterior<T>(
    psm_pep: f64,
    prior: &[f64],
    evidence: impl Iterator<Item = (f64, T)>,
    maps: impl Fn(&T, u32) -> bool,
    epsilon: f64,
) -> Posterior {
    let n_forms = prior.len();
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
    let mut l: Vec<f64> = prior.iter().map(|p| ln((1.0 - psm_pep) * p)).collect();
    l.push(ln(psm_pep));
    for (pep, item) in evidence {
        let (for_, against) = (ln(clip(1.0 - pep)), ln(clip(pep)));
        for (h, lh) in l[..n_forms].iter_mut().enumerate() {
            *lh += if maps(&item, h as u32) { for_ } else { against };
        }
        l[n_forms] += against;
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
    for x in &mut l {
        *x = (*x - mx).exp();
    }
    let z: f64 = l.iter().sum();
    let h0 = l.pop().unwrap() / z;
    for x in &mut l {
        *x /= z;
    }
    Posterior { forms: l, h0 }
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

#[cfg(test)]
mod tests {
    use super::*;
    use model::amino_acid::AminoAcid;
    use model::modification::ResidueSpec;

    fn placements(peptide: &Peptide) -> (Vec<Peptide>, bool) {
        placements_with(peptide, &[])
    }

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
            parts: Vec::new(),
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
        // Three identical mods on six S: C(6,3) = 20 forms, each once, the
        // input among them (the ordered walk must not lose any combination).
        let input = pep("SSASSASS", &[0, 4, 7]);
        let (forms, cut) = placements(&input);
        assert!(!cut);
        let mut sites: Vec<Vec<usize>> = forms.iter().map(modified_sites).collect();
        sites.sort();
        sites.dedup();
        assert_eq!((forms.len(), sites.len()), (20, 20));
        assert!(forms.iter().any(|p| same_form(p, &input)));
    }

    #[test]
    fn search_scores_set_the_placement_prior() {
        let mut e = IpfEvidence {
            forms: vec![vec![1], vec![4]],
            ..Default::default()
        };
        // No scores: uniform, PyProphet's prior.
        assert_eq!(placement_prior(&e), vec![0.5, 0.5]);
        // A higher search score for form 1 tilts the prior towards it.
        e.form_scores = vec![10.0, 14.0];
        let p = placement_prior(&e);
        assert!((p.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        let want = 1.0 / (1.0 + (-SEARCH_SCORE_WEIGHT * 4.0).exp());
        assert!((p[1] - want).abs() < 1e-12, "{p:?}");
        // No ion evidence: the posterior over forms is the prior.
        let model = IonPepModel::fit([&e]);
        let loc = localize(0.0, &e, &model);
        assert_eq!(loc.best, 1);
        assert!((loc.best_prob - want).abs() < 1e-9);
    }

    #[test]
    fn declared_siblings_extend_sites_to_other_residues() {
        // Mods.txt declares Phospho once per residue; a pS PSM must also be
        // able to move to T and Y, carrying the T/Y record there.
        let declared = [phospho(b'S'), phospho(b'T'), phospho(b'Y')];
        let (forms, _) = placements_with(&pep("ASATYK", &[1]), &declared);
        let mut sites: Vec<Vec<usize>> = forms.iter().map(modified_sites).collect();
        sites.sort();
        assert_eq!(sites, vec![vec![1], vec![3], vec![4]]);
        let y = forms.iter().find(|p| modified_sites(p) == [4]).unwrap();
        assert!(y.residues[4]
            .mod_
            .unwrap()
            .applies_to(b'Y', ModLocation::Anywhere));
    }

    #[test]
    fn neutral_loss_ions_add_site_evidence() {
        let lossy: &'static Modification = Box::leak(Box::new(Modification {
            neutral_losses: vec![97.976896],
            ..phospho(b'S').clone()
        }));
        let with = |i: usize| {
            let res = "ASATSK"
                .bytes()
                .enumerate()
                .map(|(k, r)| {
                    let mut a = AminoAcid::standard(r).unwrap();
                    if k == i {
                        a.mod_ = Some(lossy);
                    }
                    a
                })
                .collect();
            Peptide::new(res, b'K', b'A')
        };
        let (forms, _) = placements(&with(1));
        let plain = site_ions(&placements(&pep("ASATSK", &[1])).0, 2, |_| 0.01).len();
        let lossy_n = site_ions(&forms, 2, |_| 0.01).len();
        assert!(lossy_n > plain, "{lossy_n} <= {plain}");
        // y1 (K) has no phospho in any form and no loss: still dropped.
        assert!(!site_ions(&forms, 2, |_| 0.01)
            .iter()
            .any(|i| i.kind == b'y' && i.number == 1));
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
        let mut peaks: Vec<(f64, f32)> = Vec::new();
        by_ions(&truth, 2, |mz, _, _| peaks.push((mz, 100.0)));
        let mut x = 0.123_f64;
        for _ in 0..200 {
            x = (x * 9301.0 + 49297.0) % 233280.0; // LCG
            peaks.push((150.0 + 1200.0 * x / 233280.0, 5.0));
        }
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        let tol = |mz: f64| mz * 10e-6;
        let tiers = peak_tiers(&peaks);
        let e = psm_evidence(&peaks, &tiers, &reported, 2, tol, &[], None).unwrap();
        assert_eq!(e.forms.len(), 3); // S1, S4, S6
        assert_eq!(e.forms[e.observed as usize], vec![1]);
        // Fit as on a run: most confident PSMs report the right placement.
        let right = psm_evidence(&peaks, &tiers, &truth, 2, tol, &[], None).unwrap();
        let mut train = vec![&right; 45];
        train.extend(std::iter::repeat_n(&e, 5));
        let model = IonPepModel::fit(train);
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
}
