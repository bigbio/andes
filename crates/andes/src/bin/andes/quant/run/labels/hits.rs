//! Labelled hits: identified PSMs assigned to a label channel, grouped per peptide and charge,
//! with the checks for heavy proline and satellites on another channel.

use quant::label::{ChannelEnvelope, LabelSpec, Site};
use quant::Ms1RunIndex;
use search::candidates::Candidate;

pub(super) type RunSeries = Vec<(Option<f32>, Ms1RunIndex)>;

pub(super) struct LabelledHit {
    pub(super) channel: usize,
    pub(super) instances: Vec<Site>,
    pub(super) light_mass: f64,
    pub(super) key: String,
    pub(super) residues: Vec<(u8, f64)>,
}

pub(super) fn labelled(cand: &Candidate, spec: &LabelSpec) -> Option<LabelledHit> {
    let res = &cand.peptide.residues;
    let mut channel: Option<usize> = None;
    let mut instances: Vec<Site> = Vec::new();
    let mut unlabelled: Vec<Site> = Vec::new();
    let mut label_mass = 0.0;
    let mut key = String::new();
    let mut residues = Vec::with_capacity(res.len());
    let spec_sites = spec.sites();
    for (i, aa) in res.iter().enumerate() {
        let mut non_label: Vec<&str> = Vec::new();
        let mut side_taken = false;
        let mut term_taken = false;
        let mut non_label_delta = 0.0;
        let mut used: Vec<Site> = Vec::new();
        if let Some(m) = aa.mod_ {
            let mut label_delta = 0.0;
            let n_parts = m.name.split('+').count();
            for (k, part) in m.name.split('+').enumerate() {
                let terminal =
                    m.location != model::modification::ModLocation::Anywhere && k + 1 == n_parts;
                let site = [Site::Residue(aa.residue), Site::NTerm]
                    .into_iter()
                    .filter(|s| *s != Site::NTerm || i == 0)
                    .filter(|s| !used.contains(s))
                    .find(|s| spec.channel_of_label(*s, part).is_some());
                match site {
                    Some(s) => {
                        let c = spec.channel_of_label(s, part)?;
                        if channel.is_some_and(|prev| prev != c) {
                            return None;
                        }
                        channel = Some(c);
                        used.push(s);
                        instances.push(s);
                        label_delta += spec.channels[c].label_at(s)?.delta;
                    }
                    None => {
                        non_label.push(part);
                        if terminal {
                            term_taken = true;
                        } else {
                            side_taken = true;
                        }
                    }
                }
            }
            label_mass += label_delta;
            non_label_delta = m.mass_delta - label_delta;
        }
        for s in &spec_sites {
            let here = match s {
                Site::Residue(r) => *r == aa.residue && !side_taken,
                Site::NTerm => i == 0 && !term_taken,
            };
            if here && !used.contains(s) {
                unlabelled.push(*s);
            }
        }
        key.push(aa.residue as char);
        if !non_label.is_empty() {
            key.push('[');
            key.push_str(&non_label.join("+"));
            key.push(']');
        }
        residues.push((aa.residue, non_label_delta));
    }
    let channel = match channel {
        Some(c) => {
            if !unlabelled.is_empty() {
                return None;
            }
            c
        }
        None => {
            if !spec.has_unlabelled_channel() {
                return None;
            }
            instances = unlabelled;
            0
        }
    };
    if instances.is_empty() {
        return None;
    }
    Some(LabelledHit {
        channel,
        instances,
        light_mass: cand.peptide.mass() - label_mass,
        key,
        residues,
    })
}

pub(super) struct Group {
    pub(super) key: String,
    pub(super) sequence: String,
    pub(super) charge: u8,
    pub(super) proteins: Vec<String>,
    pub(super) channel_mz: Vec<f64>,
    pub(super) envelopes: Vec<ChannelEnvelope>,
    pub(super) deltas: Vec<f64>,
    pub(super) faims_cv: Option<f32>,
    pub(super) identity_decoy: bool,
    pub(super) satellites: Vec<(usize, f64)>,
}

pub(super) fn heavy_proline_delta(arg: &quant::label::LabelSite) -> Option<f64> {
    const C13_C12: f64 = model::mass::ISOTOPE;
    const N15_N14: f64 = model::mass::N15_N14;
    let n15 = arg.heavy_atoms.saturating_sub(arg.carbon13);
    (arg.carbon13 >= 5).then(|| 5.0 * C13_C12 + f64::from(n15.min(1)) * N15_N14)
}

pub(super) fn satellite_on_another_channel(g: &Group, tol_ppm: f64) -> bool {
    use quant::label::{HIGHEST_OFFSET, LOWEST_OFFSET};
    let z = f64::from(g.charge.max(1));
    let on = |mono: f64, mz: f64| {
        (LOWEST_OFFSET..=HIGHEST_OFFSET).any(|off| {
            let m = mono + f64::from(off) * model::mass::ISOTOPE / z;
            (m - mz).abs() <= m * tol_ppm * 1e-6
        })
    };
    g.satellites.iter().enumerate().any(|(i, &(parent, mz))| {
        g.channel_mz
            .iter()
            .enumerate()
            .any(|(c, &mc)| c != parent && on(mc, mz))
            || g.satellites
                .iter()
                .enumerate()
                .any(|(k, &(_, ms))| k != i && on(ms, mz))
    })
}

#[derive(Clone, Copy)]
pub(super) struct ChannelHit {
    pub(super) hi: usize,
    pub(super) rt: f64,
}
