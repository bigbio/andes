use model::mass::{H2O, PROTON};
use model::modification::{ModLocation, Modification};
use model::peptide::Peptide;

pub const MAX_PLACEMENTS: usize = 1024;

pub fn placements_with(
    peptide: &Peptide,
    declared: &[&'static Modification],
) -> (Vec<Peptide>, bool) {
    let n = peptide.residues.len();
    let movable: Vec<&'static Modification> = peptide
        .residues
        .iter()
        .filter_map(|a| a.mod_.filter(|m| is_movable(m)))
        .collect();
    if movable.is_empty() {
        return (vec![peptide.clone()], false);
    }
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
    let input_sites: Vec<usize> = (0..n)
        .filter(|&i| peptide.residues[i].mod_.is_some_and(is_movable))
        .collect();
    let run = |max_displacement: usize| {
        let mut e = Enum {
            sites: &sites,
            same_as_prev: &same_as_prev,
            base: &base,
            peptide,
            input_sites: &input_sites,
            max_displacement,
            chosen: movable.iter().map(|&m| (usize::MAX, m)).collect(),
            seen: std::collections::HashSet::new(),
            out: Vec::new(),
            truncated: false,
        };
        e.rec(0, 0);
        (e.out, e.truncated)
    };
    let (mut out, truncated) = run(usize::MAX);
    if truncated {
        let widest = n * movable.len();
        out = Vec::new();
        for d in 0..=widest {
            let (o, cut) = run(d);
            if cut {
                break;
            }
            out = o;
        }
    }
    if !out.iter().any(|p| same_form(p, peptide)) {
        out.push(peptide.clone());
    }
    (out, truncated)
}

struct Enum<'a> {
    sites: &'a [Vec<(usize, &'static Modification)>],
    same_as_prev: &'a [bool],
    base: &'a [model::amino_acid::AminoAcid],
    peptide: &'a Peptide,
    input_sites: &'a [usize],
    max_displacement: usize,
    chosen: Vec<(usize, &'static Modification)>,
    seen: std::collections::HashSet<Vec<(usize, u64)>>,
    out: Vec<Peptide>,
    truncated: bool,
}

impl Enum<'_> {
    fn rec(&mut self, k: usize, displacement: usize) {
        if self.out.len() >= MAX_PLACEMENTS {
            self.truncated = true;
            return;
        }
        if k == self.sites.len() {
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
            let d = displacement + s.abs_diff(self.input_sites[k]);
            if d > self.max_displacement {
                continue;
            }
            self.chosen[k] = (s, m);
            self.rec(k + 1, d);
        }
        self.chosen[k].0 = usize::MAX;
    }
}

fn is_movable(m: &Modification) -> bool {
    !m.fixed && m.location == ModLocation::Anywhere && m.annotated_site_mass.is_none()
}

pub fn same_form(a: &Peptide, b: &Peptide) -> bool {
    a.residues.len() == b.residues.len()
        && a.residues.iter().zip(&b.residues).all(|(x, y)| {
            x.residue == y.residue
                && x.mod_.map(|m| m.mass_delta.to_bits()) == y.mod_.map(|m| m.mass_delta.to_bits())
        })
}

pub fn modified_sites(p: &Peptide) -> Vec<usize> {
    p.residues
        .iter()
        .enumerate()
        .filter(|(_, a)| a.mod_.is_some_and(is_movable))
        .map(|(i, _)| i)
        .collect()
}

#[derive(Debug, Clone, PartialEq)]
pub struct SiteIon {
    pub mz: f64,
    pub kind: u8,
    pub number: u16,
    pub forms: Vec<u32>,
}

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

#[derive(Debug, Clone, PartialEq)]
pub struct IonMatch {
    pub forms: FormSet,
    pub tier: Option<u8>,
}

#[derive(Debug, Clone, Default)]
pub struct IpfEvidence {
    pub forms: Vec<Vec<u16>>,
    pub mods: Vec<Vec<&'static Modification>>,
    pub observed: u32,
    pub truncated: bool,
    pub ions: Vec<IonMatch>,
    pub form_scores: Vec<f32>,
}

pub const TIERS: u8 = 4;

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

fn match_tier(peaks: &[(f64, f32)], tiers: &[u8], mz: f64, tol: f64) -> Option<u8> {
    let i = peaks.partition_point(|p| p.0 < mz - tol);
    peaks[i..]
        .iter()
        .zip(&tiers[i..])
        .take_while(|(p, _)| p.0 <= mz + tol)
        .map(|(_, &t)| t)
        .min()
}

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
        mods: forms
            .iter()
            .map(|p| {
                modified_sites(p)
                    .into_iter()
                    .filter_map(|i| p.residues[i].mod_)
                    .collect()
            })
            .collect(),
        observed,
        truncated,
        ions,
    })
}

#[derive(Debug, Clone, PartialEq)]
pub struct IonPepModel {
    pub matched: [f64; TIERS as usize],
    pub unmatched: f64,
}

impl IonPepModel {
    pub fn fit<'a>(evidence: impl IntoIterator<Item = &'a IpfEvidence>) -> Self {
        Self::fit_labelled(evidence.into_iter().map(|e| (e, e.observed)))
    }

    pub fn fit_labelled<'a>(evidence: impl IntoIterator<Item = (&'a IpfEvidence, u32)>) -> Self {
        let bins = TIERS as usize + 1;
        let bin = |t: Option<u8>| t.map_or(bins - 1, |x| x as usize);
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
        for k in 1..matched.len() {
            matched[k] = matched[k].max(matched[k - 1]);
        }
        IonPepModel {
            matched,
            unmatched: pep(bins - 1),
        }
    }

    pub fn ion_pep(&self, ion: &IonMatch) -> f64 {
        ion.tier
            .map_or(self.unmatched, |t| self.matched[t as usize])
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Localization {
    pub best: u32,
    pub forms: Vec<f64>,
    pub pep: f64,
    pub best_prob: f64,
}

pub const SEARCH_SCORE_WEIGHT: f64 = 0.237;

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

pub fn localize(psm_pep: f64, e: &IpfEvidence, model: &IonPepModel) -> Localization {
    let post = log_posterior(
        psm_pep,
        &placement_prior(e),
        e.ions.iter().map(|ion| (model.ion_pep(ion), ion)),
        |ion, h| ion.forms.contains(h),
        0.0,
    );
    let best = (0..post.forms.len())
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

#[derive(Debug, Clone)]
pub struct Evidence {
    pub pep: f64,
    pub forms: Vec<u32>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Posterior {
    pub forms: Vec<f64>,
    pub h0: f64,
}

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
    fn over_the_cap_keeps_a_displacement_neighbourhood_of_the_input() {
        let input = pep(&"S".repeat(24), &[8, 10, 12, 14]);
        let (forms, cut) = placements(&input);
        assert!(cut);
        assert!(forms.len() <= MAX_PLACEMENTS && forms.len() > 1);
        assert!(forms.iter().any(|p| same_form(p, &input)));
        let disp = |sites: &[usize]| -> usize {
            sites.iter().zip([8usize, 10, 12, 14]).map(|(s, o)| s.abs_diff(o)).sum()
        };
        let kept: std::collections::HashSet<Vec<usize>> =
            forms.iter().map(modified_sites).collect();
        let bound = kept.iter().map(|s| disp(s)).max().unwrap();
        let mut within = 0;
        for a in 0..24 {
            for b in a + 1..24 {
                for c in b + 1..24 {
                    for d in c + 1..24 {
                        let s = [a, b, c, d];
                        if disp(&s) <= bound {
                            within += 1;
                            assert!(kept.contains(s.as_slice()), "{s:?} missing");
                        }
                    }
                }
            }
        }
        assert_eq!(within, kept.len());
        assert!(kept.iter().any(|s| s[0] < 8));
        assert!(kept.iter().any(|s| s[3] > 14));
    }

    #[test]
    fn placements_enumerate_every_site_once() {
        let (forms, cut) = placements(&pep("ASATSK", &[1]));
        assert!(!cut);
        let mut sites: Vec<Vec<usize>> = forms.iter().map(modified_sites).collect();
        sites.sort();
        assert_eq!(sites, vec![vec![1], vec![4]]);
        let (forms, _) = placements(&pep("SASAS", &[0, 2]));
        assert_eq!(forms.len(), 3);
        assert_eq!(placements(&pep("PEPTIDE", &[])).0.len(), 1);
        let mut annotated = pep("ASATSK", &[1]);
        let mut m = phospho(b'S').clone();
        m.annotated_site_mass = Some(m.mass_delta);
        annotated.residues[1].mod_ = Some(Box::leak(Box::new(m)));
        assert_eq!(placements(&annotated).0.len(), 1);
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
        assert_eq!(placement_prior(&e), vec![0.5, 0.5]);
        e.form_scores = vec![10.0, 14.0];
        let p = placement_prior(&e);
        assert!((p.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        let want = 1.0 / (1.0 + (-SEARCH_SCORE_WEIGHT * 4.0).exp());
        assert!((p[1] - want).abs() < 1e-12, "{p:?}");
        let model = IonPepModel::fit([&e]);
        let loc = localize(0.0, &e, &model);
        assert_eq!(loc.best, 1);
        assert!((loc.best_prob - want).abs() < 1e-9);
    }

    #[test]
    fn declared_siblings_extend_sites_to_other_residues() {
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
        assert!(!site_ions(&forms, 2, |_| 0.01)
            .iter()
            .any(|i| i.kind == b'y' && i.number == 1));
    }

    #[test]
    fn site_ions_drop_shared_ions() {
        let (forms, _) = placements(&pep("ASATSK", &[1]));
        let ions = site_ions(&forms, 2, |_| 0.01);
        assert!(!ions.is_empty());
        assert!(ions.iter().all(|i| i.forms.len() == 1));
        assert!(!ions.iter().any(|i| i.kind == b'b' && i.number == 1));
    }

    #[test]
    fn localize_recovers_the_true_site_from_a_simulated_spectrum() {
        let truth = pep("ASATSKSAPEK", &[4]);
        let reported = pep("ASATSKSAPEK", &[1]);
        let mut peaks: Vec<(f64, f32)> = Vec::new();
        by_ions(&truth, 2, |mz, _, _| peaks.push((mz, 100.0)));
        let mut x = 0.123_f64;
        for _ in 0..200 {
            x = (x * 9301.0 + 49297.0) % 233280.0;
            peaks.push((150.0 + 1200.0 * x / 233280.0, 5.0));
        }
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        let tol = |mz: f64| mz * 10e-6;
        let tiers = peak_tiers(&peaks);
        let e = psm_evidence(&peaks, &tiers, &reported, 2, tol, &[], None).unwrap();
        assert_eq!(e.forms.len(), 3);
        assert_eq!(e.forms[e.observed as usize], vec![1]);
        let right = psm_evidence(&peaks, &tiers, &truth, 2, tol, &[], None).unwrap();
        let mut train = vec![&right; 45];
        train.extend(std::iter::repeat_n(&e, 5));
        let model = IonPepModel::fit(train);
        assert!(model.matched[0] < 0.2, "{model:?}");
        let loc = localize(0.01, &e, &model);
        assert_eq!(e.forms[loc.best as usize], vec![4], "{loc:?}");
        assert!(loc.best_prob > 0.9, "{loc:?}");
    }

    #[test]
    fn posterior_matches_pyprophet_arithmetic() {
        let e = [Evidence {
            pep: 0.2,
            forms: vec![0],
        }];
        let p = posterior(0.1, 2, &e, 0.0);
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
        let q = model_fdr(&[0.1, 0.3, 0.1, 0.5]);
        let want = [0.1, 0.5 / 3.0, 0.1, 0.25];
        for (a, b) in q.iter().zip(want) {
            assert!((a - b).abs() < 1e-12, "{q:?}");
        }
    }
}
