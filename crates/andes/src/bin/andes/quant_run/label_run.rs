//! `--label`: MS1 labelled-pair quantification (SILAC, dimethyl) of the
//! confident PSMs, after rescoring.
//!
//! Every PSM is mapped onto its *label group* (the peptidoform without its
//! channel labels, at its charge) and the channel its labels name. Per run:
//!
//! 1. the elution width of the run is learned from its identified groups
//!    (rounds of fits with the width free, as [`crate::quant_mbr`] does);
//! 2. with deuterium labels, the RT offset of each channel from the
//!    reference is learned from groups whose channels are all strong when
//!    fitted one at a time ([`quant::label::learn_rt_shifts`]);
//! 3. every identified group is fitted once, all channels jointly
//!    ([`quant::label::fit_group`]), centred on its best PSM. Channels the run
//!    did not identify get a decoy fit at a non-label mass shift; target/decoy
//!    competition of the per-channel likelihood ratios over the run gives
//!    their quality q-values, and [`quant::label::classify_channel`] their
//!    status (match between channels).
//!
//! With `--mbr` (two or more runs), groups identified in some runs are
//! transferred into the others: the runs are aligned jointly on the
//! identified groups' fitted apexes ([`quant::align`]), every group missing
//! from a run is fitted there at its predicted apex, and its whole-group
//! likelihood ratio competes with a mass-shifted decoy of the group
//! ([`quant::decoys::MASS_SHIFT_UNITS`]); transfers at the `--mbr-fdr`
//! q-value are reported.
//!
//! Log-ratios are reported per channel against the reference channel, and
//! normalised by the median log-ratio of the run's identified groups for that
//! channel (a mixing error per run and channel). Transfers and matched
//! channels never enter the alignment or the normalisation reference.

use std::collections::HashMap;

use model::mass::PROTON;
use quant::formula::Formula;
use quant::label::{
    channel_score, classify_channel, classify_identified, classify_transferred, coelution_factor,
    fit_group, learn_rt_shifts, log2_ratio, mass_factor, median, weighted_qvalues, ChannelEnvelope,
    ChannelFit, ChannelStatus, ClassifyParams, FitParams, GroupFit, GroupProbe, LabelSpec, Site,
    FREE_WIDTHS,
};
use quant::Ms1RunIndex;
use rayon::prelude::*;
use search::candidate_gen::Candidate;

use super::{faims_groups, group_of, peptidoform_of, Hit, QuantCollector};
use model::spectrum::Spectrum;

/// `--label` settings.
#[derive(Debug, Clone)]
pub(crate) struct LabelSettings {
    pub spec: LabelSpec,
    /// Quality q-value cut of a matched channel.
    pub fdr: f64,
    /// MS1 tolerance of the traces, ppm.
    pub tol_ppm: f64,
    /// `--mbr`: transfer groups between runs at this q-value.
    pub mbr_fdr: Option<f64>,
    /// `--label-arg-pro`: fit the heavy proline of converted heavy arginine.
    pub arg_pro: bool,
}

impl LabelSettings {
    pub(crate) fn describe(&self) -> String {
        format!(
            "labelled pairs {} (tol {} ppm, match-between-channels q<={}{})",
            self.spec.describe(),
            self.tol_ppm,
            self.fdr,
            match self.mbr_fdr {
                Some(q) => format!(", match between runs at transfer q<={q}"),
                None => String::new(),
            }
        ) + if self.arg_pro {
            ", arginine-to-proline conversion fitted"
        } else {
            ""
        }
    }
}

/// The label mods of the search: one variable modification per (channel,
/// label site), with the channel (1-based, the unlabelled channel has none)
/// that [`model::AminoAcidSetBuilder::with_label_channels`] registers.
pub(crate) fn search_label_mods(spec: &LabelSpec) -> (Vec<model::Modification>, Vec<(String, u8)>) {
    use model::modification::{ModLocation, ResidueSpec};
    let mut mods = Vec::new();
    let mut channels = Vec::new();
    for (c, l) in spec.search_labels() {
        let (residue, location) = match l.site {
            Site::Residue(r) => (ResidueSpec::Specific(r), ModLocation::Anywhere),
            Site::NTerm => (ResidueSpec::Wildcard, ModLocation::NTerm),
        };
        mods.push(model::Modification {
            name: l.name.clone(),
            mass_delta: l.delta,
            residue,
            location,
            fixed: false,
            accession: l.unimod.map(|u| format!("UNIMOD:{u}")),
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            diagnostic_ions: Vec::new(),
            parts: Vec::new(),
        });
        if !channels.iter().any(|(n, _): &(String, u8)| *n == l.name) {
            channels.push((l.name.clone(), c as u8 + 1));
        }
    }
    (mods, channels)
}

/// Most labelled sites a peptide of the search can carry. In a specific or
/// semi-specific digest (`num_tolerable_termini` 2 or 1), a labelled residue
/// that every enzyme cleaves after (SILAC K/R with trypsin) is a cleavage site
/// or the C-terminus, so those number at most the missed cleavages plus one.
/// Labels on other residues (Leu3, or Arg with Lys-C), and every label of a
/// non-specific digest (where missed cleavages do not apply), number at most
/// the most such residues in any `max_length` residues of `db`, which bounds
/// every peptide the search can produce; plus the N-terminus. The labels count
/// apart from the other variable mods (`SearchParams::max_label_sites`), so
/// this bounds the labels only.
pub(crate) fn max_label_sites(
    spec: &LabelSpec,
    enzymes: &[model::enzyme::Enzyme],
    num_tolerable_termini: u8,
    missed_cleavages: u32,
    db: &model::ProteinDb,
    max_length: u32,
) -> u32 {
    use model::enzyme::Enzyme;
    let sites = spec.sites();
    let specific = num_tolerable_termini > 0
        && !enzymes.is_empty()
        && enzymes
            .iter()
            .all(|e| !matches!(e, Enzyme::NonSpecific | Enzyme::NoCleavage));
    let mut cleavage = false;
    let mut free: Vec<u8> = Vec::new();
    for s in &sites {
        if let Site::Residue(r) = s {
            if specific && enzymes.iter().all(|e| e.is_cleavable_after(*r)) {
                cleavage = true;
            } else {
                free.push(*r);
            }
        }
    }
    let nterm = sites.contains(&Site::NTerm);
    let bound = if cleavage { missed_cleavages + 1 } else { 0 };
    let free = if free.is_empty() {
        0
    } else {
        most_in_window(db, &free, max_length as usize)
    };
    bound + free + u32::from(nterm)
}

/// The most residues of `residues` in any `len` consecutive residues of a
/// protein of `db`.
fn most_in_window(db: &model::ProteinDb, residues: &[u8], len: usize) -> u32 {
    let len = len.max(1);
    let hit = |r: &u8| residues.contains(&r.to_ascii_uppercase());
    let mut best = 0u32;
    for p in &db.proteins {
        let s = &p.sequence;
        let mut n = 0u32;
        for (i, r) in s.iter().enumerate() {
            n += u32::from(hit(r));
            if i >= len {
                n -= u32::from(hit(&s[i - len]));
            }
            best = best.max(n);
        }
    }
    best
}

/// A run's MS1 index split by FAIMS compensation voltage ([`faims_groups`]).
type RunSeries = Vec<(Option<f32>, Ms1RunIndex)>;

/// A PSM in label terms.
struct LabelledHit {
    channel: usize,
    /// Site of every label instance of the peptide (labelled or, in the
    /// unlabelled channel, labelable).
    instances: Vec<Site>,
    light_mass: f64,
    key: String,
    /// Residues with their non-label modification deltas (for the formula).
    residues: Vec<(u8, f64)>,
}

/// The label group and channel of a candidate; `None` when its labels mix
/// channels or it has no labelable site.
fn labelled(cand: &Candidate, spec: &LabelSpec) -> Option<LabelledHit> {
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
        // Whether a non-label part holds the side chain / the terminus.
        let mut side_taken = false;
        let mut term_taken = false;
        let mut non_label_delta = 0.0;
        let mut used: Vec<Site> = Vec::new();
        if let Some(m) = aa.mod_ {
            let mut label_delta = 0.0;
            let n_parts = m.name.split('+').count();
            for (k, part) in m.name.split('+').enumerate() {
                // A stack lists its side-chain parts first and its terminal
                // part last, and takes its location from that terminal part.
                let terminal =
                    m.location != model::modification::ModLocation::Anywhere && k + 1 == n_parts;
                // The site of a label part: this residue's, then the N-terminus.
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
        // Labelable sites left without a label: a non-label mod takes a
        // residue's side chain or the N-terminus, not both.
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

/// One label group: a peptide identity at one charge.
struct Group {
    key: String,
    sequence: String,
    charge: u8,
    proteins: Vec<String>,
    channel_mz: Vec<f64>,
    envelopes: Vec<ChannelEnvelope>,
    /// Channel masses relative to channel 1 (Da).
    deltas: Vec<f64>,
    /// FAIMS voltage of its best PSM, if any.
    faims_cv: Option<f32>,
    /// A decoy peptide's group (its decoy PSMs at the targets' cut): measured
    /// like a target, only to count how often the transfer of a false donor
    /// identification is accepted; never reported.
    identity_decoy: bool,
    /// `--label-arg-pro`: per satellite, its parent channel and m/z (the
    /// parent with one heavy proline from converted heavy arginine).
    satellites: Vec<(usize, f64)>,
}

/// Mass of the heavy proline that heavy arginine converts to: 13C5 for Arg6
/// (13C6), 13C5 15N1 for Arg10 (13C6 15N4).
fn heavy_proline_delta(arg: &quant::label::LabelSite) -> Option<f64> {
    const C13_C12: f64 = 1.003_354_838;
    const N15_N14: f64 = 0.997_034_886;
    let n15 = arg.heavy_atoms.saturating_sub(arg.carbon13);
    (arg.carbon13 >= 5).then(|| 5.0 * C13_C12 + f64::from(n15.min(1)) * N15_N14)
}

/// Whether a satellite of `g` lies within `tol_ppm` of an isotope the fit
/// models (M−2..M+5) of another component: a channel other than its parent,
/// or another satellite. The two fit components are then nearly collinear
/// and their split is weakly determined.
fn satellite_on_another_channel(g: &Group, tol_ppm: f64) -> bool {
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

/// The best PSM of a group's channel in a run.
#[derive(Clone, Copy)]
struct ChannelHit {
    hi: usize,
    rt: f64,
}

/// A fitted (group, run).
struct Fitted {
    gid: usize,
    run: usize,
    fit: GroupFit,
    /// Best PSM per channel in this run (identified channels).
    ids: Vec<Option<ChannelHit>>,
    /// Per channel, the decoy fits' scores and monoisotopic mass errors
    /// (missing channels).
    decoy_scores: Vec<Vec<(f64, Option<f64>, f64)>>,
    /// Per channel, another identified peptide of the run sits on it.
    claimed: Vec<bool>,
    transfer_q: Option<f64>,
    rt_shift: Vec<f64>,
    /// Per channel, the converted (heavy proline) fraction (`--label-arg-pro`).
    converted: Vec<Option<f64>>,
}

/// The run-level scales of one MS1 series.
struct Scales {
    dt: f64,
    w: f64,
    shifts: Vec<f64>,
}

const MIN_WIDTH_LLR: f64 = 25.0;
const MIN_WIDTH_FITS: usize = 20;
const MAX_WIDTH_FITS: usize = 500;
const MIN_SHIFT_PAIRS: usize = 20;
const MIN_NORM_RATIOS: usize = 10;
/// Decoys per target of the channel q-values.
const N_DECOYS: usize = 4;
/// Decoys per transfer: a run has fewer transfers than channels, so more
/// decoys resolve the q-value at the cut.
const N_TRANSFER_DECOYS: usize = 8;

fn probe_of(g: &Group, t_hat: f64, s: f64, w: f64, shifts: &[f64]) -> GroupProbe {
    let mut channel_mz = g.channel_mz.clone();
    let mut envelopes = g.envelopes.clone();
    let mut rt_shift = shifts.to_vec();
    for &(parent, mz) in &g.satellites {
        channel_mz.push(mz);
        envelopes.push(g.envelopes[parent]);
        rt_shift.push(shifts[parent]);
    }
    GroupProbe {
        charge: g.charge,
        free: vec![true; channel_mz.len()],
        channel_mz,
        envelopes,
        t_hat,
        s,
        w_prior: w,
        rt_shift,
        mobility: None,
        pinned: false,
    }
}

/// [`fit_group`] with the satellites folded into their parent channels: a
/// parent's amplitude and quantity become the sums, its standard error their
/// root sum of squares. Also returns each channel's converted fraction (the
/// satellites' share; `None` without satellites).
fn fit_folded(
    ms1: &Ms1RunIndex,
    g: &Group,
    probe: &GroupProbe,
    params: &FitParams,
) -> Option<(GroupFit, Vec<Option<f64>>)> {
    let mut fit = fit_group(ms1, probe, params)?;
    let n = g.channel_mz.len();
    let mut converted = vec![None; n];
    if g.satellites.is_empty() {
        return Some((fit, converted));
    }
    let sats: Vec<ChannelFit> = fit.channels.split_off(n);
    for (&(parent, _), sat) in g.satellites.iter().zip(&sats) {
        let c = &mut fit.channels[parent];
        let total = c.amplitude + sat.amplitude;
        converted[parent] = Some(if total > 0.0 {
            sat.amplitude / total
        } else {
            0.0
        });
        c.amplitude = total;
        c.quantity += sat.quantity;
        c.trace_quantity += sat.trace_quantity;
        c.amplitude_se = c.amplitude_se.hypot(sat.amplitude_se);
        c.quantity_se = c.quantity_se.hypot(sat.quantity_se);
        c.llr += sat.llr;
    }
    Some((fit, converted))
}

/// The joint fit at the identified channels' profile (with its satellite
/// conversions) and the pinned probe it was measured with.
type ProfileFit = ((GroupFit, Vec<Option<f64>>), GroupProbe);

/// The fit of a group with unidentified channels at the elution profile of
/// its identified channels: the apex and width are searched with only the
/// identified channels (and their satellites) free, the others held at 0,
/// and every channel is then measured jointly at that profile, pinned. Also
/// returns the pinned probe, from which the decoys of the unidentified
/// channels are made: neither an unidentified channel nor a decoy can move
/// the profile onto a species at its own m/z. `rt_z` is the profile fit's.
fn fit_at_identified_profile(
    ms1: &Ms1RunIndex,
    g: &Group,
    probe: &GroupProbe,
    ids: &[Option<ChannelHit>],
    params: &FitParams,
) -> Option<ProfileFit> {
    let n = g.channel_mz.len();
    let mut shape = probe.clone();
    shape.free = (0..shape.channel_mz.len())
        .map(|k| ids[if k < n { k } else { g.satellites[k - n].0 }].is_some())
        .collect();
    let (profile, _) = fit_folded(ms1, g, &shape, params)?;
    let mut pinned = probe.clone();
    pinned.t_hat = profile.t0;
    pinned.w_prior = profile.w;
    pinned.pinned = true;
    let (mut fit, converted) = fit_folded(ms1, g, &pinned, params)?;
    fit.rt_z = profile.rt_z;
    Some(((fit, converted), pinned))
}

/// Learns the run's elution width from `anchors` (group, RT of its best PSM
/// in reference time): rounds of fits with the width free, the median
/// fitted width of the strong fits, until it moves by less than 10 %.
fn learn_width(ms1: &Ms1RunIndex, groups: &[Group], anchors: &[(usize, f64)], tol: f64) -> f64 {
    let dt = ms1.median_scan_spacing().unwrap_or(1.0).max(1e-3);
    // The sample is chosen by a hash of each group's key, so one group more
    // or fewer changes at most one member (a stride over the list would shift
    // all of them).
    let mut order: Vec<&(usize, f64)> = anchors.iter().collect();
    order.sort_by_key(|(gid, _)| seed_of(&groups[*gid].key, 0));
    let sample: Vec<&(usize, f64)> = order.into_iter().take(MAX_WIDTH_FITS).collect();
    // The coarse grid finds the width's scale; the fine one, once the coarse
    // estimate has settled, resolves it to a few percent.
    let mut params = FitParams {
        tol_ppm: tol,
        widths: FREE_WIDTHS,
    };
    let mut w = 6.0f64.max(dt);
    for _ in 0..16 {
        let s = w.max(dt);
        let mut fitted: Vec<f64> = sample
            .par_iter()
            .filter_map(|(gid, t)| {
                let g = &groups[*gid];
                let p = probe_of(g, *t, s, w, &vec![0.0; g.channel_mz.len()]);
                fit_folded(ms1, g, &p, &params).map(|x| x.0)
            })
            .filter(|f| f.llr >= MIN_WIDTH_LLR)
            .map(|f| f.w.ln())
            .collect();
        if fitted.len() < MIN_WIDTH_FITS {
            return (6.0f64).max(dt);
        }
        // Each fit picks its width from a grid relative to `w`: their median
        // is a grid point, at which the iteration would stop (3.6 or 4.2 s
        // from 6 s, by the sample). The trimmed mean of the log widths moves
        // continuously and settles where the grid choices balance.
        fitted.sort_by(f64::total_cmp);
        let cut = fitted.len() / 10;
        let core = &fitted[cut..fitted.len() - cut];
        let next = (core.iter().sum::<f64>() / core.len() as f64).exp().max(dt);
        let moved = (next - w).abs() / w;
        w = next;
        if params.widths == FREE_WIDTHS {
            if moved < 0.1 {
                params.widths = REFINE_WIDTHS;
            }
        } else if moved < 0.01 {
            break;
        }
    }
    w
}

/// The width grid of [`learn_width`]'s second stage, multiples of the
/// current estimate.
const REFINE_WIDTHS: &[f64] = &[0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.25];

/// Learns the channels' RT offsets from the reference: every channel fitted
/// alone at `anchors`, pairs where both are strong.
fn learn_shifts(
    ms1: &Ms1RunIndex,
    groups: &[Group],
    anchors: &[(usize, f64)],
    scales: (f64, f64),
    spec: &LabelSpec,
    tol: f64,
) -> Vec<f64> {
    let n = spec.len();
    if !spec.has_deuterium() {
        return vec![0.0; n];
    }
    let (dt, w) = scales;
    let params = FitParams {
        tol_ppm: tol,
        widths: quant::label::DEFAULT_WIDTHS,
    };
    let step = (anchors.len() / 2000).max(1);
    let r = spec.reference;
    let pairs: Vec<(usize, f64)> = anchors
        .par_iter()
        .step_by(step)
        .flat_map_iter(|(gid, t)| {
            let g = &groups[*gid];
            let alone = |c: usize| {
                let mut p = probe_of(g, *t, w.max(dt), w, &vec![0.0; n]);
                p.free = (0..p.channel_mz.len())
                    .map(|k| k == c || (k >= n && g.satellites[k - n].0 == c))
                    .collect();
                fit_folded(ms1, g, &p, &params)
                    .map(|x| x.0)
                    .filter(|f| f.llr >= MIN_WIDTH_LLR)
            };
            let Some(fr) = alone(r) else {
                return Vec::new();
            };
            (0..n)
                .filter(|&c| c != r)
                .filter_map(|c| alone(c).map(|f| (c, f.t0 - fr.t0)))
                .collect::<Vec<_>>()
        })
        .collect();
    learn_rt_shifts(&pairs, n, MIN_SHIFT_PAIRS)
}

/// The decoys of `probe` with its channels `moved` (the unidentified
/// channels, or every channel for a transfer) moved by the non-label mass
/// shifts of [`quant::label::decoy_shifts`].
fn decoy_probes(
    g: &Group,
    probe: &GroupProbe,
    moved: &[usize],
    tol_ppm: f64,
    seed: u64,
    n: usize,
) -> Vec<GroupProbe> {
    let z = g.charge.max(1) as f64;
    let mass = (g.channel_mz[0] - PROTON) * z;
    let min_sep = quant::label::decoy_min_sep(mass, tol_ppm);
    quant::label::decoy_shifts(moved, &g.deltas, min_sep, seed, n)
        .into_iter()
        .map(|shift| {
            let mut p = probe.clone();
            for &c in moved {
                p.channel_mz[c] += shift / z;
            }
            // Satellites follow the parents in the probe; they move with them.
            let n_ch = g.channel_mz.len();
            for (k, &(parent, _)) in g.satellites.iter().enumerate() {
                if moved.contains(&parent) {
                    p.channel_mz[n_ch + k] += shift / z;
                }
            }
            p
        })
        .collect()
}

/// The competition score of a group: its strongest channel's
/// [`channel_score`] with the group's total likelihood ratio.
fn group_score(fit: &GroupFit) -> f64 {
    let Some(best) = fit
        .channels
        .iter()
        .max_by(|a, b| a.amplitude.total_cmp(&b.amplitude))
    else {
        return 0.0;
    };
    (2.0 * fit.llr.max(0.0)).sqrt() * best.envelope_cos.max(0.0) * best.coherence.max(0.0)
}

fn seed_of(s: &str, run: usize) -> u64 {
    // FNV-1a: stable across runs of the program.
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for b in s.bytes().chain((run as u64).to_le_bytes()) {
        h ^= b as u64;
        h = h.wrapping_mul(0x0100_0000_01b3);
    }
    h
}

/// Further `label.parquet` metadata keys and values.
type Metadata = Vec<(String, String)>;

/// `label.parquet` rows, QPX feature records and further table metadata.
pub(super) type LabelOutputs = (
    Vec<output::LabelRow>,
    Vec<output::FeatureRecord>,
    Vec<(String, String)>,
);

/// Per-run label quantification. Returns the `label.parquet` rows, the QPX
/// feature records (identified and transferred groups, quantified channels)
/// and the identity-decoy metadata of the transfers. `decoy_hits` are the
/// decoy PSMs at the targets' cut (identity decoys, `--mbr` only).
#[allow(clippy::too_many_arguments)]
pub(super) fn quantify(
    collector: &QuantCollector,
    settings: &LabelSettings,
    spectra: &[Spectrum],
    target_hits: &[Hit<'_>],
    decoy_hits: &[Hit<'_>],
    candidates: &[Candidate],
    index: &search::search_index::SearchIndex,
    with_records: bool,
    next_feature_id: &mut i64,
) -> Result<LabelOutputs, String> {
    let all_hits: Vec<&Hit<'_>> = target_hits.iter().chain(decoy_hits).collect();
    let hits = &all_hits[..];
    let spec = &settings.spec;
    let n_ch = spec.len();
    let sites = spec.sites();
    let tol = settings.tol_ppm;

    // ── Groups from the PSMs ──────────────────────────────────────────────
    let mut groups: Vec<Group> = Vec::new();
    let mut gid_of: HashMap<(String, u8), usize> = HashMap::new();
    // (gid, run) -> best hit per channel.
    let mut ids: HashMap<(usize, usize), Vec<Option<ChannelHit>>> = HashMap::new();
    let mut best_of_group: HashMap<usize, usize> = HashMap::new();
    let (mut mixed, mut no_rt) = (0usize, 0usize);
    let better = |h: &Hit<'_>, o: &Hit<'_>| match (h.q_value, o.q_value) {
        (Some(a), Some(b)) if a != b => a < b,
        _ => h.psm.rank_score > o.psm.rank_score,
    };
    for (hi, h) in hits.iter().enumerate() {
        let Some(mut lh) = labelled(h.cand, spec) else {
            mixed += 1;
            continue;
        };
        let Some(rt) = spectra[h.spec_idx].rt_seconds else {
            no_rt += 1;
            continue;
        };
        let z = h.psm.charge_used.max(1);
        if h.cand.is_decoy {
            lh.key = format!("DECOY:{}", lh.key);
        }
        let gid = *gid_of.entry((lh.key.clone(), z)).or_insert_with(|| {
            let counts: Vec<u32> = sites
                .iter()
                .map(|s| lh.instances.iter().filter(|x| *x == s).count() as u32)
                .collect();
            let light = Formula::peptide(lh.residues.iter().copied());
            let mut channel_mz = Vec::with_capacity(n_ch);
            let mut envelopes = Vec::with_capacity(n_ch);
            let mut deltas = Vec::with_capacity(n_ch);
            let base = spec.channel_delta(0, &sites, &counts);
            for c in 0..n_ch {
                let delta = spec.channel_delta(c, &sites, &counts);
                channel_mz.push((lh.light_mass + delta) / z as f64 + PROTON);
                envelopes.push(spec.channel_envelope(c, &sites, &counts, &light));
                deltas.push(delta - base);
            }
            groups.push(Group {
                key: lh.key.clone(),
                sequence: output::qpx::bare_sequence(h.cand),
                charge: z,
                proteins: super::proteins_of(&h.psm, candidates, index)
                    .into_iter()
                    .map(|(a, ..)| a)
                    .collect(),
                satellites: if settings.arg_pro {
                    let n_pro = lh.residues.iter().filter(|r| r.0 == b'P').count();
                    (0..n_ch)
                        .filter(|_| n_pro > 0)
                        .filter_map(|c| {
                            let arg = spec.channels[c].label_at(Site::Residue(b'R'))?;
                            let d = heavy_proline_delta(arg)?;
                            Some((c, channel_mz[c] + d / z as f64))
                        })
                        .collect()
                } else {
                    Vec::new()
                },
                channel_mz,
                envelopes,
                deltas,
                faims_cv: None,
                identity_decoy: h.cand.is_decoy,
            });
            groups.len() - 1
        });
        let slot = ids
            .entry((gid, h.run_idx))
            .or_insert_with(|| vec![None; n_ch]);
        let replace = slot[lh.channel].is_none_or(|o| better(h, hits[o.hi]));
        if replace {
            slot[lh.channel] = Some(ChannelHit { hi, rt });
        }
        let b = best_of_group.entry(gid).or_insert(hi);
        if better(h, hits[*b]) {
            *b = hi;
        }
    }
    for (gid, &hi) in &best_of_group {
        let h = &hits[hi];
        let run = &collector.runs[h.run_idx];
        groups[*gid].faims_cv = run
            .contexts
            .get(super::native_id(run, &spectra[h.spec_idx].title))
            .and_then(|c| c.faims_cv);
    }
    eprintln!(
        "quant: label: {} label groups from {} PSMs ({} without a consistent label state, \
         {} without a retention time)",
        groups.len(),
        hits.len(),
        mixed,
        no_rt
    );
    let shared = groups
        .iter()
        .filter(|g| !g.identity_decoy && satellite_on_another_channel(g, tol))
        .count();
    if shared > 0 {
        eprintln!(
            "quant: label: {shared} groups have an Arg->Pro satellite within the MS1 tolerance \
             of another channel's isotope; the split between the two (converted_fraction, and \
             the ratio through it) is weakly determined there"
        );
    }

    let fit_params = FitParams {
        tol_ppm: tol,
        ..FitParams::default()
    };
    let classify = ClassifyParams {
        fdr: settings.fdr,
        ..ClassifyParams::default()
    };
    let n_runs = collector.runs.len();
    let mut fitted: Vec<Fitted> = Vec::new();
    // Per (run, series): scales, kept for the transfers.
    let mut run_scales: Vec<Vec<Option<Scales>>> = Vec::with_capacity(n_runs);
    let mut run_series_cv: Vec<Vec<Option<f32>>> = Vec::with_capacity(n_runs);

    // With match between runs, the last run fitted stays loaded so that the
    // transfers into it need not decode it again (one run resident at a time).
    let keep_resident = settings.mbr_fdr.is_some() && n_runs >= 2;
    let mut resident: Option<(usize, RunSeries)> = None;
    for (run_idx, run) in collector.runs.iter().enumerate() {
        let Some(cache) = run.ms1.as_ref() else {
            run_scales.push(Vec::new());
            run_series_cv.push(Vec::new());
            continue;
        };
        resident = None;
        let ms1 = cache
            .load()
            .map_err(|e| format!("load MS1 scans for {}: {e}", run.stem))?;
        let series = faims_groups(ms1);
        let mut scales_here: Vec<Option<Scales>> = Vec::new();
        run_series_cv.push(series.iter().map(|(cv, _)| *cv).collect());
        // Identified groups of this run, per series.
        let mut by_series: Vec<Vec<(usize, Vec<Option<ChannelHit>>)>> =
            vec![Vec::new(); series.len()];
        let mut keys: Vec<&(usize, usize)> = ids.keys().filter(|(_, r)| *r == run_idx).collect();
        keys.sort_unstable();
        for &(gid, _) in keys {
            let chans = ids[&(gid, run_idx)].clone();
            let best = chans
                .iter()
                .flatten()
                .copied()
                .reduce(|a, b| if better(hits[b.hi], hits[a.hi]) { b } else { a })
                .expect("an identified group has a PSM");
            let spec_ = &spectra[hits[best.hi].spec_idx];
            let Some(si) = group_of(&series, run, spec_) else {
                continue;
            };
            by_series[si].push((gid, chans));
        }
        // Identified precursors of the run (claims): (series, charge, mz, rt, gid).
        let mut claims: Vec<(usize, u8, f64, f64, usize)> = Vec::new();
        for (si, list) in by_series.iter().enumerate() {
            for (gid, chans) in list.iter().filter(|(g, _)| !groups[*g].identity_decoy) {
                for (c, ch) in chans.iter().enumerate() {
                    if let Some(ch) = ch {
                        let g = &groups[*gid];
                        claims.push((si, g.charge, g.channel_mz[c], ch.rt, *gid));
                    }
                }
            }
        }
        // Other identified PSMs of the run (not label groups) claim too.
        for h in hits
            .iter()
            .filter(|h| h.run_idx == run_idx && !h.cand.is_decoy)
        {
            if labelled(h.cand, spec).is_none() {
                let spec_ = &spectra[h.spec_idx];
                if let (Some(si), Some(rt)) = (group_of(&series, run, spec_), spec_.rt_seconds) {
                    let z = h.psm.charge_used.max(1);
                    claims.push((
                        si,
                        z,
                        h.cand.peptide.mass() / z as f64 + PROTON,
                        rt,
                        usize::MAX,
                    ));
                }
            }
        }
        // By m/z, for the claim search.
        claims.sort_by(|a, b| a.2.total_cmp(&b.2));
        for (si, (_, ms1)) in series.iter().enumerate() {
            let list = &by_series[si];
            if list.is_empty() {
                scales_here.push(None);
                continue;
            }
            let dt = ms1.median_scan_spacing().unwrap_or(1.0).max(1e-3);
            // Anchors in the reference channel's time (no shifts yet),
            // target groups only.
            let anchors: Vec<(usize, f64)> = list
                .iter()
                .filter(|(gid, _)| !groups[*gid].identity_decoy)
                .map(|(gid, chans)| {
                    let best = chans
                        .iter()
                        .flatten()
                        .copied()
                        .reduce(|a, b| if better(hits[b.hi], hits[a.hi]) { b } else { a })
                        .unwrap();
                    (*gid, best.rt)
                })
                .collect();
            let w = learn_width(ms1, &groups, &anchors, tol);
            let shifts = learn_shifts(ms1, &groups, &anchors, (dt, w), spec, tol);
            eprintln!(
                "quant: label: {}{}: {} identified groups, elution width {:.2} s, MS1 spacing \
                 {:.2} s{}",
                run.stem,
                match series[si].0 {
                    Some(cv) => format!(" (FAIMS {cv} V)"),
                    None => String::new(),
                },
                list.len(),
                w,
                dt,
                if spec.has_deuterium() {
                    format!(
                        ", channel RT offsets {}",
                        shifts
                            .iter()
                            .map(|v| format!("{v:+.2} s"))
                            .collect::<Vec<_>>()
                            .join(" ")
                    )
                } else {
                    String::new()
                }
            );
            let s = w.max(dt);
            let rows: Vec<Fitted> = list
                .par_iter()
                .filter_map(|(gid, chans)| {
                    let g = &groups[*gid];
                    let best = chans
                        .iter()
                        .enumerate()
                        .filter_map(|(c, x)| x.map(|x| (c, x)))
                        .reduce(|a, b| {
                            if better(hits[b.1.hi], hits[a.1.hi]) {
                                b
                            } else {
                                a
                            }
                        })?;
                    let t_hat = best.1.rt - shifts[best.0];
                    let probe = probe_of(g, t_hat, s, w, &shifts);
                    let missing: Vec<usize> = (0..n_ch).filter(|&c| chans[c].is_none()).collect();
                    let ((fit, converted), probe) = if missing.is_empty() {
                        (fit_folded(ms1, g, &probe, &fit_params)?, probe)
                    } else {
                        let ((mut fit, conv), pinned) =
                            fit_at_identified_profile(ms1, g, &probe, chans, &fit_params)?;
                        // The identified channels keep their quantities from the
                        // free joint fit, as a group with every channel
                        // identified does; only the unidentified channels (and
                        // their decoys, made from `pinned`) are read at the
                        // identified channels' profile.
                        if let Some((free, _)) = fit_folded(ms1, g, &probe, &fit_params) {
                            let hits = chans.iter();
                            for ((ch, f), hit) in
                                fit.channels.iter_mut().zip(&free.channels).zip(hits)
                            {
                                if hit.is_some() {
                                    *ch = f.clone();
                                }
                            }
                        }
                        ((fit, conv), pinned)
                    };
                    let mut decoy_scores = vec![Vec::new(); n_ch];
                    if !missing.is_empty() {
                        let seed = seed_of(&g.key, run_idx);
                        // The decoys move the missing channels' masses of the
                        // pinned probe: measured at the same profile as the
                        // channels they stand for, which none of them shaped.
                        for dp in decoy_probes(g, &probe, &missing, tol, seed, N_DECOYS) {
                            let df = fit_folded(ms1, g, &dp, &fit_params).map(|x| x.0);
                            for &c in &missing {
                                decoy_scores[c].push(df.as_ref().map_or((0.0, None, 0.0), |d| {
                                    let ch = &d.channels[c];
                                    (channel_score(ch), ch.mass_err_ppm, ch.apex_offset)
                                }));
                            }
                        }
                    }
                    let claimed = (0..n_ch)
                        .map(|c| {
                            chans[c].is_none() && {
                                let mz = g.channel_mz[c];
                                let tolda = mz * tol * 1e-6;
                                let t = fit.t0 + shifts[c];
                                let lo = claims.partition_point(|c| c.2 < mz - tolda);
                                claims[lo..].iter().take_while(|c| c.2 <= mz + tolda).any(
                                    |&(csi, z, cmz, crt, cg)| {
                                        csi == si
                                            && z == g.charge
                                            && cg != *gid
                                            && (cmz - mz).abs() <= tolda
                                            && (crt - t).abs() <= 3.0 * fit.w + dt
                                    },
                                )
                            }
                        })
                        .collect();
                    Some(Fitted {
                        gid: *gid,
                        run: run_idx,
                        fit,
                        ids: chans.clone(),
                        decoy_scores,
                        claimed,
                        transfer_q: None,
                        rt_shift: shifts.clone(),
                        converted,
                    })
                })
                .collect();
            fitted.extend(rows);
            scales_here.push(Some(Scales { dt, w, shifts }));
        }
        run_scales.push(scales_here);
        if keep_resident {
            resident = Some((run_idx, series));
        }
    }

    // ── Match between runs ────────────────────────────────────────────────
    // ── Mass calibration per run: the median and robust spread of the
    // identified channels' monoisotopic mass errors (strong fits).
    let calib: Vec<(f64, f64)> = (0..n_runs)
        .map(|r| {
            let v: Vec<f64> = fitted
                .iter()
                .filter(|f| f.run == r && f.transfer_q.is_none())
                .flat_map(|f| {
                    (0..n_ch)
                        .filter(|&c| f.ids[c].is_some() && f.fit.channels[c].llr >= MIN_WIDTH_LLR)
                        .filter_map(|c| f.fit.channels[c].mass_err_ppm)
                        .collect::<Vec<_>>()
                })
                .collect();
            if v.len() < MIN_WIDTH_FITS {
                return (0.0, tol / 3.0);
            }
            let m = median(v.clone());
            let mad = median(v.iter().map(|x| (x - m).abs()).collect());
            (m, (1.4826 * mad).max(0.5))
        })
        .collect();
    // Spread of the identified channels' apex offsets per run, floored at a
    // quarter of the MS1 spacing.
    let offset_sd: Vec<f64> = (0..n_runs)
        .map(|r| {
            let v: Vec<f64> = fitted
                .iter()
                .filter(|f| f.run == r && f.transfer_q.is_none())
                .flat_map(|f| {
                    (0..n_ch)
                        .filter(|&c| f.ids[c].is_some() && f.fit.channels[c].llr >= MIN_WIDTH_LLR)
                        .map(|c| f.fit.channels[c].apex_offset)
                        .collect::<Vec<_>>()
                })
                .collect();
            let dt = run_scales
                .get(r)
                .and_then(|s| s.iter().flatten().map(|x| x.dt).reduce(f64::min))
                .unwrap_or(1.0);
            if v.len() < MIN_WIDTH_FITS {
                return dt;
            }
            let m = median(v.clone());
            (1.4826 * median(v.iter().map(|x| (x - m).abs()).collect())).max(0.25 * dt)
        })
        .collect();
    for (r, (m, sd)) in calib.iter().enumerate() {
        eprintln!(
            "quant: label: {}: mass error of the identified channels {m:+.2} ppm, spread {sd:.2} \
             ppm; apex offset spread {:.2} s",
            collector.runs[r].stem, offset_sd[r]
        );
    }

    let mut meta: Vec<(String, String)> = Vec::new();
    if let (Some(mbr_fdr), true) = (settings.mbr_fdr, n_runs >= 2) {
        let (transfers, stats) = transfer(
            collector,
            &groups,
            &fitted,
            &run_scales,
            &run_series_cv,
            &fit_params,
            &calib,
            mbr_fdr,
            resident.take(),
        )?;
        fitted.extend(transfers);
        meta = stats;
    }

    // ── Channel q-values (match between channels), per run ─────────────────
    // (fitted index, channel) of every unidentified channel of an identified
    // group, its score against its decoys'. Each run competes on its own:
    // runs differ in MS1 density and design (a 1:10 run's weak light channels
    // would drag a 1:1 run's q-values down).
    let mut channel_q: HashMap<(usize, usize), f64> = HashMap::new();
    for run_idx in 0..n_runs {
        let mut rows: Vec<(f64, f64)> = Vec::new();
        let mut at: Vec<((usize, usize), usize)> = Vec::new();
        for (fi, f) in fitted.iter().enumerate() {
            if f.run != run_idx || f.transfer_q.is_some() || groups[f.gid].identity_decoy {
                continue;
            }
            for c in 0..n_ch {
                if f.ids[c].is_some() || f.decoy_scores[c].is_empty() {
                    continue;
                }
                at.push(((fi, c), rows.len()));
                // Expected mass error: the group's identified channels'.
                let (center, sigma) = calib[f.run];
                let expected = {
                    let v: Vec<f64> = (0..n_ch)
                        .filter(|&k| f.ids[k].is_some())
                        .filter_map(|k| f.fit.channels[k].mass_err_ppm)
                        .collect();
                    if v.is_empty() {
                        center
                    } else {
                        median(v)
                    }
                };
                // Two channels' errors differ by √2 of one's spread.
                let sd = sigma * std::f64::consts::SQRT_2;
                let ch = &f.fit.channels[c];
                // Two channels' apex offsets differ by √2 of one's spread.
                let st = offset_sd[f.run] * std::f64::consts::SQRT_2;
                rows.push((
                    channel_score(ch)
                        * mass_factor(ch.mass_err_ppm, expected, sd)
                        * coelution_factor(ch.apex_offset, st),
                    0.0,
                ));
                let w = 1.0 / f.decoy_scores[c].len() as f64;
                rows.extend(f.decoy_scores[c].iter().map(|&(d, e, o)| {
                    (
                        d * mass_factor(e, expected, sd) * coelution_factor(o, st),
                        w,
                    )
                }));
            }
        }
        let q = weighted_qvalues(&rows);
        for (key, row) in at {
            channel_q.insert(key, q[row]);
        }
    }

    // ── Statuses, ratios, normalisation ───────────────────────────────────
    let r = spec.reference;
    let status_of = |fi: usize, f: &Fitted| -> Vec<ChannelStatus> {
        let strongest = (0..n_ch)
            .filter(|&c| f.ids[c].is_some() || f.transfer_q.is_some())
            .map(|c| f.fit.channels[c].amplitude)
            .fold(0.0f64, f64::max);
        (0..n_ch)
            .map(|c| {
                let ch: &ChannelFit = &f.fit.channels[c];
                if f.ids[c].is_some() {
                    classify_identified(ch, &classify)
                } else if f.transfer_q.is_some() {
                    classify_transferred(ch, strongest, &classify)
                } else {
                    let q = channel_q.get(&(fi, c)).copied().unwrap_or(1.0);
                    classify_channel(ch, q, strongest, f.claimed[c], &classify)
                }
            })
            .collect()
    };
    let statuses: Vec<Vec<ChannelStatus>> = fitted
        .iter()
        .enumerate()
        .map(|(fi, f)| status_of(fi, f))
        .collect();
    // Normalisation offsets per (run, channel): identified groups with both
    // the channel and the reference measured, not transferred.
    let mut offsets: HashMap<(usize, usize), f64> = HashMap::new();
    for run_idx in 0..n_runs {
        for c in (0..n_ch).filter(|&c| c != r) {
            let v: Vec<f64> = fitted
                .iter()
                .zip(&statuses)
                .filter(|(f, st)| {
                    f.run == run_idx
                        && f.transfer_q.is_none()
                        && !groups[f.gid].identity_decoy
                        && st[c] == ChannelStatus::Identified
                        && st[r] == ChannelStatus::Identified
                })
                .filter_map(|(f, _)| log2_ratio(&f.fit.channels[c], &f.fit.channels[r]))
                .map(|x| x.0)
                .collect();
            if v.len() >= MIN_NORM_RATIOS {
                offsets.insert((run_idx, c), quant::label::median(v));
            }
        }
    }

    // ── Rows ──────────────────────────────────────────────────────────────
    let mut order: Vec<usize> = (0..fitted.len()).collect();
    order.sort_by_key(|&i| (fitted[i].run, fitted[i].gid));
    let mut rows: Vec<output::LabelRow> = Vec::new();
    let mut records: Vec<output::FeatureRecord> = Vec::new();
    let mut counts = vec![[0usize; 5]; n_ch];
    let any_transfer = fitted.iter().any(|f| f.transfer_q.is_some());
    for &fi in &order {
        let f = &fitted[fi];
        let g = &groups[f.gid];
        if g.identity_decoy {
            continue;
        }
        let st = &statuses[fi];
        let run = &collector.runs[f.run];
        for c in 0..n_ch {
            let ch = &f.fit.channels[c];
            let ratio = (c != r && st[c].is_quantified() && st[r].is_quantified())
                .then(|| log2_ratio(ch, &f.fit.channels[r]))
                .flatten();
            let id_hit = f.ids[c].map(|x| &hits[x.hi]);
            counts[c][match st[c] {
                ChannelStatus::Identified => 0,
                ChannelStatus::Matched => 1,
                ChannelStatus::Transferred => 2,
                ChannelStatus::BelowDetection => 3,
                ChannelStatus::NotQuantifiable => 4,
            }] += 1;
            rows.push(output::LabelRow {
                run: run.stem.clone(),
                peptide: g.key.clone(),
                sequence: g.sequence.clone(),
                charge: g.charge,
                proteins: g.proteins.clone(),
                channel: c as u8 + 1,
                channel_name: spec.channels[c].name.clone(),
                label: spec.channels[c].qpx_label.clone(),
                status: st[c].as_str().to_string(),
                spec_id: id_hit.map(|h| h.spec_id.clone()),
                psm_q_value: id_hit.and_then(|h| h.q_value),
                rt_psm: f.ids[c].map(|x| x.rt),
                labelled_peptidoform: id_hit.map(|h| peptidoform_of(h.cand)),
                mz: g.channel_mz[c],
                rt_apex: f.fit.t0 + f.rt_shift[c],
                rt_shift: f.rt_shift[c],
                width: f.fit.w,
                quantity: ch.quantity,
                quantity_se: ch.quantity_se,
                quantity_upper: ch.quantity + 2.0 * ch.quantity_se,
                llr: ch.llr,
                envelope_cos: ch.envelope_cos,
                coherence: ch.coherence,
                m1_ratio: ch.m1_ratio,
                mass_err_ppm: ch.mass_err_ppm,
                channel_q_value: channel_q.get(&(fi, c)).copied(),
                group_llr: f.fit.llr,
                transfer: f.transfer_q.is_some(),
                transfer_q_value: f.transfer_q,
                log2_ratio: ratio.map(|x| x.0),
                log2_ratio_se: ratio.map(|x| x.1),
                log2_ratio_normalized: ratio
                    .and_then(|x| offsets.get(&(f.run, c)).map(|o| x.0 - o)),
                converted_fraction: f.converted[c],
                trace_quantity: ch.trace_quantity,
            });
        }
        // QPX: one feature per (group, run), quantified channels. A transfer
        // carries the identification of the group's best PSM in another run.
        if with_records {
            let best = if f.transfer_q.is_some() {
                best_of_group.get(&f.gid).copied()
            } else {
                f.ids
                    .iter()
                    .flatten()
                    .copied()
                    .reduce(|a, b| if better(hits[b.hi], hits[a.hi]) { b } else { a })
                    .map(|x| x.hi)
            };
            if let Some(hi) = best {
                let hit = hits[hi];
                let spec_ = &spectra[hit.spec_idx];
                let donor_run = &collector.runs[hit.run_idx];
                let id = super::quant_id(hit, donor_run, spec_, candidates, index);
                let proteins = super::proteins_of(&hit.psm, candidates, index);
                let mut rec = super::base_record(&id, proteins, *next_feature_id);
                *next_feature_id += 1;
                if f.transfer_q.is_some() {
                    rec.run_file_name = run.stem.clone();
                    rec.id_run_file_name = None;
                    rec.psm_ids = Vec::new();
                    rec.scan = Vec::new();
                    // The donor's peptidoform and m/z, observed in this run
                    // at that channel's fitted mass error.
                    let c = labelled(hit.cand, spec).map_or(0, |l| l.channel);
                    let calc = g.channel_mz[c];
                    // No centroid there: no observed m/z and no mass error.
                    let err = f.fit.channels[c].mass_err_ppm;
                    rec.calculated_mz = calc as f32;
                    rec.observed_mz = (calc * (1.0 + err.unwrap_or(0.0) * 1e-6)) as f32;
                    rec.mass_error_ppm = err.map(|e| e as f32);
                }
                if any_transfer {
                    rec.mbr = Some(output::feature_parquet::MbrColumns {
                        evidence_type: if f.transfer_q.is_some() {
                            "transfer_ms1"
                        } else {
                            "identified"
                        }
                        .to_string(),
                        transfer_q: f.transfer_q,
                        quant_q: None,
                        quantity_se: None,
                        claimed_by: None,
                        claim_co_apex: None,
                    });
                }
                rec.rt = Some(f.fit.t0 as f32);
                rec.rt_start = Some((f.fit.t0 - 3.0 * f.fit.w) as f32);
                rec.rt_stop = Some((f.fit.t0 + 3.0 * f.fit.w) as f32);
                rec.intensities = (0..n_ch)
                    .filter(|&c| st[c].is_quantified())
                    .filter_map(|c| {
                        spec.channels[c]
                            .qpx_label
                            .clone()
                            .map(|l| (l, f.fit.channels[c].quantity as f32))
                    })
                    .collect();
                rec.additional_scores
                    .push(("label_group_llr".into(), f.fit.llr, true));
                for (c, status) in st.iter().enumerate() {
                    rec.cv_params.push((
                        format!("label_status_{}", c + 1),
                        status.as_str().to_string(),
                    ));
                }
                rec.cv_params.push(("label_peptide".into(), g.key.clone()));
                if !rec.intensities.is_empty() {
                    records.push(rec);
                }
            }
        }
    }
    for (c, k) in counts.iter().enumerate() {
        eprintln!(
            "quant: label: channel {} ({}): {} identified, {} matched between channels, {} \
             transferred, {} below detection, {} not quantifiable",
            c + 1,
            spec.channels[c].name,
            k[0],
            k[1],
            k[2],
            k[3],
            k[4]
        );
    }
    for ((run_idx, c), o) in {
        let mut v: Vec<_> = offsets.iter().collect();
        v.sort_by_key(|x| *x.0);
        v
    } {
        eprintln!(
            "quant: label: {} channel {} median log2 ratio to channel {}: {o:+.3}",
            collector.runs[*run_idx].stem,
            c + 1,
            r + 1
        );
    }
    Ok((rows, records, meta))
}

/// Match between runs: every group of the library (identified in some run)
/// fitted in each run that lacks it, at its aligned apex; the transfers at
/// the q-value cut.
#[allow(clippy::too_many_arguments)]
fn transfer(
    collector: &QuantCollector,
    groups: &[Group],
    fitted: &[Fitted],
    run_scales: &[Vec<Option<Scales>>],
    run_series_cv: &[Vec<Option<f32>>],
    fit_params: &FitParams,
    calib: &[(f64, f64)],
    fdr: f64,
    mut resident: Option<(usize, RunSeries)>,
) -> Result<(Vec<Fitted>, Metadata), String> {
    // The group score with the mass accuracy of its strongest channel
    // against the receiving run's calibration.
    let score = |fit: &GroupFit, run: usize| -> f64 {
        let (center, sigma) = calib[run];
        let err = fit
            .channels
            .iter()
            .max_by(|a, b| a.amplitude.total_cmp(&b.amplitude))
            .and_then(|c| c.mass_err_ppm);
        // The alignment is calibrated (a real transfer's apex is within its
        // uncertainty of the prediction, `rt_z` ~ N(0, 1)); a decoy that hits
        // another species' envelope sits anywhere in the search window.
        group_score(fit)
            * quant::label::mass_factor(err, center, sigma)
            * (-0.5 * fit.rt_z * fit.rt_z).exp()
    };
    let n_runs = collector.runs.len();
    let n_ch = groups.first().map_or(0, |g| g.channel_mz.len());
    // Observations: the fitted apex (reference time) of strong identified fits.
    let observe = |decoy: bool| -> Vec<quant::align::Observation> {
        fitted
            .iter()
            .filter(|f| f.fit.llr >= MIN_WIDTH_LLR && groups[f.gid].identity_decoy == decoy)
            .map(|f| quant::align::Observation {
                precursor: f.gid as u32,
                run: f.run as u16,
                rt: f.fit.t0,
            })
            .collect()
    };
    let obs = observe(false);
    let spacing: Vec<f64> = run_scales
        .iter()
        .map(|s| {
            s.iter()
                .flatten()
                .map(|x| x.dt / 2.0)
                .fold(f64::NAN, f64::min)
        })
        .collect();
    let params = quant::align::AlignParams {
        run_min_spacing_s: spacing,
        ..Default::default()
    };
    let mut alignment = quant::align::align(&obs, groups.len(), n_runs, &params);
    // Identity decoys get their consensus and predictions as a target would,
    // without moving any run map.
    alignment.add_observations(&observe(true))?;
    let mut present = vec![vec![false; n_runs]; groups.len()];
    for f in fitted {
        present[f.gid][f.run] = true;
    }
    let mut out: Vec<Fitted> = Vec::new();
    let mut scored: Vec<(f64, f64)> = Vec::new();
    let mut target_rows: Vec<usize> = Vec::new();
    let mut identity_scores: Vec<f64> = Vec::new();
    // The run still loaded from the fits first; the results are then taken
    // in run order, as if every run had been visited in turn.
    let mut order: Vec<usize> = (0..n_runs).collect();
    if let Some((r, _)) = &resident {
        order.retain(|x| x != r);
        order.insert(0, *r);
    }
    let mut batches: Vec<Vec<(Fitted, Vec<f64>)>> = (0..n_runs).map(|_| Vec::new()).collect();
    for run_idx in order {
        let run = &collector.runs[run_idx];
        let Some(cache) = run.ms1.as_ref() else {
            continue;
        };
        if alignment.component(run_idx as u16).is_none() {
            eprintln!(
                "quant: label: {}: not aligned (too few shared identified groups); no transfers \
                 into it",
                run.stem
            );
            continue;
        }
        let series = match resident.take() {
            Some((r, series)) if r == run_idx => series,
            _ => faims_groups(
                cache
                    .load()
                    .map_err(|e| format!("load MS1 scans for {}: {e}", run.stem))?,
            ),
        };
        // Width scale of a series without identifications: the run's other
        // series, else the median of every run.
        let fallback = || {
            let all: Vec<f64> = run_scales.iter().flatten().flatten().map(|s| s.w).collect();
            quant::label::median(all)
        };
        let jobs: Vec<(usize, usize, f64, f64)> = (0..groups.len())
            .filter(|&gid| !present[gid][run_idx])
            .filter_map(|gid| {
                let (t, s_align) = alignment.predict(gid as u32, run_idx as u16)?;
                let si = run_series_cv
                    .get(run_idx)
                    .and_then(|cvs| {
                        if cvs.len() <= 1 {
                            Some(0)
                        } else {
                            cvs.iter().position(|cv| *cv == groups[gid].faims_cv)
                        }
                    })
                    .filter(|&si| si < series.len())?;
                Some((gid, si, t, s_align))
            })
            .collect();
        let results: Vec<(Fitted, Vec<f64>)> = jobs
            .par_iter()
            .filter_map(|&(gid, si, t, s_align)| {
                let ms1 = &series[si].1;
                let sc = run_scales
                    .get(run_idx)
                    .and_then(|v| v.get(si))
                    .and_then(|x| x.as_ref());
                let dt = sc.map_or_else(|| ms1.median_scan_spacing().unwrap_or(1.0), |s| s.dt);
                let w = sc.map_or_else(fallback, |s| s.w).max(dt);
                let shifts = sc.map_or_else(|| vec![0.0; n_ch], |s| s.shifts.clone());
                let g = &groups[gid];
                let probe = probe_of(g, t, s_align.max(dt), w, &shifts);
                let (fit, converted) = fit_folded(ms1, g, &probe, fit_params)?;
                let all: Vec<usize> = (0..n_ch).collect();
                let dscores: Vec<f64> = decoy_probes(
                    g,
                    &probe,
                    &all,
                    fit_params.tol_ppm,
                    seed_of(&g.key, run_idx),
                    N_TRANSFER_DECOYS,
                )
                .iter()
                .map(|d| fit_folded(ms1, g, d, fit_params).map_or(0.0, |d| score(&d.0, run_idx)))
                .collect();
                Some((
                    Fitted {
                        gid,
                        run: run_idx,
                        fit,
                        ids: vec![None; n_ch],
                        decoy_scores: vec![Vec::new(); n_ch],
                        claimed: vec![false; n_ch],
                        transfer_q: None,
                        rt_shift: shifts,
                        converted,
                    },
                    dscores,
                ))
            })
            .collect();
        batches[run_idx] = results;
    }
    for (f, d) in batches.into_iter().flatten() {
        if groups[f.gid].identity_decoy {
            identity_scores.push(score(&f.fit, f.run));
            continue;
        }
        target_rows.push(scored.len());
        scored.push((score(&f.fit, f.run), 0.0));
        let w = 1.0 / d.len().max(1) as f64;
        scored.extend(d.into_iter().map(|x| (x, w)));
        out.push(f);
    }
    let q = weighted_qvalues(&scored);
    // The q-value of an identity-decoy transfer: that of the weakest target
    // transfer scoring at least as high (the targets' step function); one
    // above every target takes the best target's.
    let mut steps: Vec<(f64, f64)> = target_rows.iter().map(|&r| (scored[r].0, q[r])).collect();
    steps.sort_by(|a, b| b.0.total_cmp(&a.0));
    let q_at = |score: f64| -> f64 {
        steps
            .iter()
            .take_while(|(s, _)| *s >= score)
            .last()
            .or(steps.first())
            .map_or(1.0, |x| x.1)
    };
    let accepted_targets = target_rows.iter().filter(|&&r| q[r] <= fdr).count();
    let accepted_identity = identity_scores.iter().filter(|&&s| q_at(s) <= fdr).count();
    // Accepted transfers of false donor identifications: identity decoys
    // stand for the library's false target entries one for one.
    let identity_rate = accepted_identity as f64 / accepted_targets.max(1) as f64;
    let mut kept = Vec::new();
    for (k, mut f) in out.into_iter().enumerate() {
        let qk = (q[target_rows[k]] + identity_rate).min(1.0);
        if qk <= fdr {
            f.transfer_q = Some(qk);
            kept.push(f);
        }
    }
    eprintln!(
        "quant: label: match between runs: {} of {} transfers at transfer q<={fdr} (lowest \
         q {:.4}; identity decoys: {} of {} transferred, false-donor rate {:.4} added to the \
         q-values)",
        kept.len(),
        target_rows.len(),
        target_rows.iter().map(|&r| q[r]).fold(1.0f64, f64::min),
        accepted_identity,
        identity_scores.len(),
        identity_rate
    );
    let stats = vec![
        (
            "andes.label.identity_decoy_transfers".to_string(),
            format!("{accepted_identity}/{}", identity_scores.len()),
        ),
        (
            "andes.label.identity_decoy_rate".to_string(),
            identity_rate.to_string(),
        ),
    ];
    Ok((kept, stats))
}

/// Fewest peptide ratios a protein ratio needs (MaxQuant's "min. ratio
/// count").
pub(super) const MIN_RATIO_COUNT: usize = 2;

/// One protein ratio of [`protein_ratios`].
#[derive(Debug, Clone, PartialEq)]
pub(super) struct ProteinRatio {
    pub protein: String,
    pub lead: String,
    pub run: String,
    pub channel: u8,
    pub channel_name: String,
    pub log2_ratio: f64,
    pub log2_ratio_normalized: Option<f64>,
    /// Peptide ratios the median is over, and those whose evidence is only
    /// matched between channels or transferred (no identification of the
    /// channel in the run).
    pub ratio_count: usize,
    pub inferred_count: usize,
}

/// Protein ratios as MaxQuant computes them, from `label.parquet` rows: per
/// run, channel and peptide identity (all charges), the median of its
/// evidence log-ratios; per protein, the median of its peptides' ratios, from
/// [`MIN_RATIO_COUNT`] peptides up; the normalised ratio likewise from the
/// normalised evidence. A peptide's protein is its protein group at the
/// protein q-value cut when `group_of` gives one, else its accessions joined.
pub(super) fn protein_ratios(
    rows: &[output::LabelRow],
    group_of: impl Fn(&str) -> Option<(String, String)>,
    min_count: usize,
) -> Vec<ProteinRatio> {
    use std::collections::BTreeMap;
    // (run, protein, lead, channel, name) -> peptide -> (ratios, normalised, identified)
    type Key = (String, String, String, u8, String);
    // Per peptide: its evidence ratios, normalised ratios, and whether the
    // channel is identified in some evidence.
    type Evidence = (Vec<f64>, Vec<f64>, bool);
    let mut by: BTreeMap<Key, BTreeMap<String, Evidence>> = BTreeMap::new();
    for r in rows {
        let Some(v) = r.log2_ratio else { continue };
        let (protein, lead) = group_of(&r.sequence).unwrap_or_else(|| {
            let p = r.proteins.join(";");
            (p.clone(), r.proteins.first().cloned().unwrap_or(p))
        });
        let e = by
            .entry((
                r.run.clone(),
                protein,
                lead,
                r.channel,
                r.channel_name.clone(),
            ))
            .or_default()
            .entry(r.peptide.clone())
            .or_default();
        e.0.push(v);
        if let Some(n) = r.log2_ratio_normalized {
            e.1.push(n);
        }
        e.2 |= r.status == "identified";
    }
    let mut out = Vec::new();
    for ((run, protein, lead, channel, channel_name), peptides) in by {
        let ratios: Vec<f64> = peptides.values().map(|p| median(p.0.clone())).collect();
        if ratios.len() < min_count {
            continue;
        }
        let norm: Vec<f64> = peptides
            .values()
            .filter(|p| !p.1.is_empty())
            .map(|p| median(p.1.clone()))
            .collect();
        out.push(ProteinRatio {
            protein,
            lead,
            run,
            channel,
            channel_name,
            log2_ratio: median(ratios.clone()),
            log2_ratio_normalized: (norm.len() >= min_count).then(|| median(norm)),
            ratio_count: ratios.len(),
            inferred_count: peptides.values().filter(|p| !p.2).count(),
        });
    }
    out
}

/// `<stem>.label_proteins.tsv`: the [`protein_ratios`] as text.
pub(super) fn protein_ratio_tsv(ratios: &[ProteinRatio]) -> String {
    let mut s = String::from(
        "protein_group\tlead\tfilename\tchannel\tchannel_name\tlog2_ratio\t\
         log2_ratio_normalized\tratio_count\tinferred_ratio_count\n",
    );
    for r in ratios {
        s.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\n",
            r.protein,
            r.lead,
            r.run,
            r.channel,
            r.channel_name,
            r.log2_ratio,
            r.log2_ratio_normalized
                .map_or_else(|| "NA".to_string(), |v| format!("{v:.6}")),
            r.ratio_count,
            r.inferred_count
        ));
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A heavy channel 20 % converted to its heavy-proline form: fitted with
    /// the satellite, the channel holds the whole heavy amount.
    #[test]
    fn arg_pro_satellite_is_folded_into_its_channel() {
        use model::scan::Ms1Scan;
        let spec = LabelSpec::parse("silac2").unwrap();
        let arg10 = spec.channels[1].label_at(Site::Residue(b'R')).unwrap();
        let d_pro = heavy_proline_delta(arg10).unwrap();
        assert!((d_pro - 6.013809).abs() < 1e-4, "{d_pro}");
        let light = Formula::peptide(b"PEPTIDER".iter().map(|&r| (r, 0.0)));
        let sites = spec.sites();
        let counts: Vec<u32> = sites
            .iter()
            .map(|s| u32::from(*s == Site::Residue(b'R')))
            .collect();
        let envs: Vec<ChannelEnvelope> = (0..2)
            .map(|c| spec.channel_envelope(c, &sites, &counts, &light))
            .collect();
        let z = 2.0;
        let ml = 478.73;
        let mh = ml + 10.008269 / z;
        let mp = mh + d_pro / z;
        // (m/z, envelope, amplitude) at apex 100 s, sigma 4 s.
        let species = [
            (ml, envs[0], 1e6),
            (mh, envs[1], 0.8e6),
            (mp, envs[1], 0.2e6),
        ];
        let scans: Vec<Ms1Scan> = (0..200)
            .map(|i| {
                let rt = i as f64;
                let g = (-(rt - 100.0f64).powi(2) / 32.0).exp();
                let mut peaks: Vec<(f64, f32)> =
                    (0..40).map(|k| (300.0 + k as f64 * 23.1, 50.0)).collect();
                for (mz, env, amp) in species {
                    for off in -2..=5 {
                        let v = amp * env.at(off) * g;
                        if v > 50.0 {
                            peaks.push((mz + off as f64 * model::mass::ISOTOPE / z, v as f32));
                        }
                    }
                }
                peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
                Ms1Scan {
                    rt,
                    peaks,
                    id: format!("scan={i}"),
                    ..Default::default()
                }
            })
            .collect();
        let ms1 = Ms1RunIndex::new(scans);
        let mut g = Group {
            key: "PEPTIDER".into(),
            sequence: "PEPTIDER".into(),
            charge: 2,
            proteins: Vec::new(),
            channel_mz: vec![ml, mh],
            envelopes: envs.clone(),
            deltas: vec![0.0, 10.008269],
            faims_cv: None,
            satellites: vec![(1, mp)],
            identity_decoy: false,
        };
        // SILAC duplex: the satellite is clear of the light channel.
        assert!(!satellite_on_another_channel(&g, 10.0));
        let p = probe_of(&g, 100.0, 4.0, 4.0, &[0.0, 0.0]);
        let (fit, converted) = fit_folded(&ms1, &g, &p, &FitParams::default()).unwrap();
        let (r, _) = log2_ratio(&fit.channels[1], &fit.channels[0]).unwrap();
        assert!(r.abs() < 0.02, "{r}");
        assert!((converted[1].unwrap() - 0.2).abs() < 0.01, "{converted:?}");
        // Without the satellite, the heavy channel misses the converted part.
        g.satellites.clear();
        let p = probe_of(&g, 100.0, 4.0, 4.0, &[0.0, 0.0]);
        let (fit, _) = fit_folded(&ms1, &g, &p, &FitParams::default()).unwrap();
        let (r, _) = log2_ratio(&fit.channels[1], &fit.channels[0]).unwrap();
        assert!((r - 0.8f64.log2()).abs() < 0.02, "{r}");
    }

    /// A light-only group whose heavy channel is absent, with a strong
    /// species eluting 6 s later either at the heavy channel's m/z (the
    /// target's situation) or at a decoy's (the decoy's): the target's heavy
    /// channel scores no better than that decoy, since both are measured at
    /// the light channel's profile, which neither moves.
    #[test]
    fn missing_channel_scores_like_its_decoys_next_to_an_interferer() {
        use model::scan::Ms1Scan;
        let env = ChannelEnvelope::new(&[0.55, 0.3, 0.11, 0.04], 0, 1.0, 0, 1.0);
        let z = 2.0;
        let ml = 512.27;
        let g = Group {
            key: "PEPTIDEK".into(),
            sequence: "PEPTIDEK".into(),
            charge: 2,
            proteins: Vec::new(),
            channel_mz: vec![ml, ml + 8.014199 / z],
            envelopes: vec![env, env],
            deltas: vec![0.0, 8.014199],
            faims_cv: None,
            satellites: Vec::new(),
            identity_decoy: false,
        };
        let ids = vec![Some(ChannelHit { hi: 0, rt: 100.0 }), None];
        let probe = probe_of(&g, 100.0, 4.0, 4.0, &[0.0, 0.0]);
        let params = FitParams::default();
        // The run: the light channel at apex 100 s, and a species 5x as
        // strong at `mz`, apex 106 s; sigma 4 s.
        let run = |mz: f64| {
            let scans: Vec<Ms1Scan> = (0..600)
                .map(|k| {
                    let rt = k as f64 / 3.0;
                    let mut peaks: Vec<(f64, f32)> =
                        (0..60).map(|j| (350.0 + 11.3 * j as f64, 50.0)).collect();
                    for (m, apex, amp) in [(ml, 100.0, 1e6), (mz, 106.0, 5e6)] {
                        let e = (-(rt - apex) * (rt - apex) / 32.0f64).exp();
                        for off in 0..4 {
                            let v = amp * env.at(off) * e;
                            if v > 50.0 {
                                peaks.push((m + off as f64 * model::mass::ISOTOPE / z, v as f32));
                            }
                        }
                    }
                    peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
                    Ms1Scan {
                        rt,
                        peaks,
                        id: format!("scan={k}"),
                        ..Default::default()
                    }
                })
                .collect();
            Ms1RunIndex::new(scans)
        };
        let decoy_mz = decoy_probes(&g, &probe, &[1], 10.0, 7, 1)[0].channel_mz[1];
        assert!((decoy_mz - g.channel_mz[1]).abs() > 0.01);

        // The target's situation: the species sits on the heavy channel. A
        // free joint fit is pulled onto it; the measurement is not.
        let at_heavy = run(g.channel_mz[1]);
        let (free, _) = fit_folded(&at_heavy, &g, &probe, &params).unwrap();
        assert!(free.t0 > 102.0, "free apex {}", free.t0);
        let ((target, _), _) =
            fit_at_identified_profile(&at_heavy, &g, &probe, &ids, &params).unwrap();
        assert!((target.t0 - 100.0).abs() < 0.5, "apex {}", target.t0);

        // The decoy's situation: the same species on the decoy's m/z.
        let at_decoy = run(decoy_mz);
        let (_, pinned) = fit_at_identified_profile(&at_decoy, &g, &probe, &ids, &params).unwrap();
        let dp = decoy_probes(&g, &pinned, &[1], 10.0, 7, 1).remove(0);
        assert_eq!(dp.channel_mz[1], decoy_mz);
        let (decoy, _) = fit_folded(&at_decoy, &g, &dp, &params).unwrap();

        let (st, sd) = (
            channel_score(&target.channels[1]),
            channel_score(&decoy.channels[1]),
        );
        // Read at the same pinned profile, the target and the decoy see the
        // species alike (6 s off the apex their coherence, and so their score,
        // can be 0): neither scores better than the other.
        assert!(
            (st - sd).abs() <= 0.02 * st.max(sd) + 1e-9,
            "target {st} vs decoy {sd}"
        );
        assert!(
            (target.channels[1].apex_offset - decoy.channels[1].apex_offset).abs() < 0.2,
            "{} vs {}",
            target.channels[1].apex_offset,
            decoy.channels[1].apex_offset
        );
        // The identified channel is measured at its own profile.
        let (r, _) = log2_ratio(&target.channels[0], &decoy.channels[0]).unwrap();
        assert!(r.abs() < 0.02, "{r}");
    }

    /// The learned elution width finds the true width and does not jump when
    /// one group more or fewer is fitted.
    #[test]
    fn learned_width_is_stable() {
        use model::scan::Ms1Scan;
        let env = ChannelEnvelope::new(&[0.55, 0.3, 0.11, 0.04], 0, 1.0, 0, 1.0);
        let n_groups = 60;
        let (sigma, z) = (4.0, 2.0);
        let groups: Vec<Group> = (0..n_groups)
            .map(|i| {
                let ml = 400.0 + 7.3 * i as f64;
                Group {
                    key: format!("PEPTIDE{i}K"),
                    sequence: format!("PEPTIDE{i}K"),
                    charge: 2,
                    proteins: Vec::new(),
                    channel_mz: vec![ml, ml + 8.014199 / z],
                    envelopes: vec![env, env],
                    deltas: vec![0.0, 8.014199],
                    faims_cv: None,
                    satellites: Vec::new(),
                    identity_decoy: false,
                }
            })
            .collect();
        let apex = |i: usize| 40.0 + 9.0 * i as f64;
        let scans: Vec<Ms1Scan> = (0..1800)
            .map(|k| {
                let rt = k as f64 / 3.0;
                let mut peaks: Vec<(f64, f32)> =
                    (0..60).map(|j| (350.0 + 11.3 * j as f64, 50.0)).collect();
                for (i, g) in groups.iter().enumerate() {
                    let e = (-(rt - apex(i)).powi(2) / (2.0 * sigma * sigma)).exp();
                    for (c, amp) in [(0, 1e6), (1, 5e5)] {
                        for off in 0..4 {
                            let v = amp * env.at(off) * e;
                            if v > 50.0 {
                                let m = g.channel_mz[c] + off as f64 * model::mass::ISOTOPE / z;
                                peaks.push((m, v as f32));
                            }
                        }
                    }
                }
                peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
                Ms1Scan {
                    rt,
                    peaks,
                    id: format!("scan={k}"),
                    ..Default::default()
                }
            })
            .collect();
        let ms1 = Ms1RunIndex::new(scans);
        let anchors: Vec<(usize, f64)> = (0..n_groups).map(|i| (i, apex(i))).collect();
        let all = learn_width(&ms1, &groups, &anchors, 10.0);
        let fewer = learn_width(&ms1, &groups, &anchors[1..], 10.0);
        assert!((all - sigma).abs() < 0.04 * sigma, "learned {all}");
        assert!((all - fewer).abs() < 0.02 * all, "{all} vs {fewer}");
    }

    /// The collinearity check covers every isotope the fit models (M−2..M+5)
    /// of the other channels, and the other satellites.
    #[test]
    fn satellite_overlap_check_covers_the_modelled_isotopes() {
        let iso = model::mass::ISOTOPE / 2.0;
        let group = |channel_mz: Vec<f64>, satellites: Vec<(usize, f64)>| Group {
            key: "PEPTIDER".into(),
            sequence: "PEPTIDER".into(),
            charge: 2,
            proteins: Vec::new(),
            envelopes: vec![ChannelEnvelope::new(&[1.0], 0, 1.0, 0, 1.0); channel_mz.len()],
            deltas: vec![0.0; channel_mz.len()],
            channel_mz,
            faims_cv: None,
            satellites,
            identity_decoy: false,
        };
        // On channel 1's M−1 (a lower isotopologue of a 2H label).
        assert!(satellite_on_another_channel(
            &group(vec![500.0, 510.0], vec![(0, 510.0 - iso)]),
            10.0
        ));
        // On channel 1's M+5.
        assert!(satellite_on_another_channel(
            &group(vec![500.0, 510.0], vec![(0, 510.0 + 5.0 * iso)]),
            10.0
        ));
        // On another satellite.
        assert!(satellite_on_another_channel(
            &group(vec![500.0, 520.0], vec![(0, 540.0), (1, 540.0 + iso)]),
            10.0
        ));
        // Clear of everything.
        assert!(!satellite_on_another_channel(
            &group(vec![500.0, 520.0], vec![(0, 540.0), (1, 560.0)]),
            10.0
        ));
    }

    /// Four channels with three heavy-proline satellites take seven fit
    /// components, more than the five channels a specification may have:
    /// the group is still fitted, every channel with its converted part.
    #[test]
    fn four_channels_with_arg_pro_satellites_are_fitted() {
        use model::scan::Ms1Scan;
        let spec = LabelSpec::parse("[][Arg6][R+8.014199:6C2N][Arg10]").unwrap();
        let light = Formula::peptide(b"PEPTIDER".iter().map(|&r| (r, 0.0)));
        let sites = spec.sites();
        let counts: Vec<u32> = sites
            .iter()
            .map(|s| u32::from(*s == Site::Residue(b'R')))
            .collect();
        let envs: Vec<ChannelEnvelope> = (0..4)
            .map(|c| spec.channel_envelope(c, &sites, &counts, &light))
            .collect();
        let z = 2.0;
        let ml = 478.73;
        let deltas = [0.0, 6.020129, 8.014199, 10.008269];
        let mz: Vec<f64> = deltas.iter().map(|d| ml + d / z).collect();
        let sats: Vec<(usize, f64)> = (1..4)
            .map(|c| {
                let arg = spec.channels[c].label_at(Site::Residue(b'R')).unwrap();
                (c, mz[c] + heavy_proline_delta(arg).unwrap() / z)
            })
            .collect();
        // Channel amplitudes 1 : 0.5 : 0.25 : 2, each heavy one 20 % converted.
        let amps = [1e6, 0.5e6, 0.25e6, 2e6];
        let mut species: Vec<(f64, ChannelEnvelope, f64)> = vec![(mz[0], envs[0], amps[0])];
        for &(c, mp) in &sats {
            species.push((mz[c], envs[c], 0.8 * amps[c]));
            species.push((mp, envs[c], 0.2 * amps[c]));
        }
        let scans: Vec<Ms1Scan> = (0..200)
            .map(|i| {
                let rt = i as f64;
                let g = (-(rt - 100.0f64).powi(2) / 32.0).exp();
                let mut peaks: Vec<(f64, f32)> =
                    (0..40).map(|k| (300.0 + k as f64 * 23.1, 50.0)).collect();
                for &(m, env, amp) in &species {
                    for off in -2..=5 {
                        let v = amp * env.at(off) * g;
                        if v > 50.0 {
                            peaks.push((m + off as f64 * model::mass::ISOTOPE / z, v as f32));
                        }
                    }
                }
                peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
                Ms1Scan {
                    rt,
                    peaks,
                    id: format!("scan={i}"),
                    ..Default::default()
                }
            })
            .collect();
        let ms1 = Ms1RunIndex::new(scans);
        let g = Group {
            key: "PEPTIDER".into(),
            sequence: "PEPTIDER".into(),
            charge: 2,
            proteins: Vec::new(),
            channel_mz: mz.clone(),
            envelopes: envs.clone(),
            deltas: deltas.to_vec(),
            faims_cv: None,
            satellites: sats,
            identity_decoy: false,
        };
        // Channel 2's satellite sits on Arg10's M+4 (3 mDa): flagged.
        assert!(satellite_on_another_channel(&g, 10.0));
        let p = probe_of(&g, 100.0, 4.0, 4.0, &[0.0; 4]);
        assert_eq!(p.channel_mz.len(), 7);
        let (fit, converted) = fit_folded(&ms1, &g, &p, &FitParams::default()).unwrap();
        for c in 1..4 {
            let (r, _) = log2_ratio(&fit.channels[c], &fit.channels[0]).unwrap();
            let want = (amps[c] / amps[0]).log2();
            assert!((r - want).abs() < 0.05, "channel {c}: {r} vs {want}");
            // Channel 2's satellite (+7.014 at z 2) sits 3 mDa from Arg10's
            // M+4, inseparable at any MS1 tolerance: a little of it is shared.
            assert!((converted[c].unwrap() - 0.2).abs() < 0.03, "{converted:?}");
        }
    }

    /// A protein N-terminal acetyl holds the N-terminus, not the first
    /// residue's side chain: light and heavy forms of an acetylated
    /// `KPEPTIDER` name the same label sites (K and R).
    #[test]
    fn terminal_mod_leaves_the_residue_label_site() {
        use model::amino_acid::AminoAcid;
        use model::modification::{ModLocation, Modification, ResidueSpec};
        use model::peptide::Peptide;
        let spec = LabelSpec::parse("silac2").unwrap();
        let m = |name: &str, delta: f64, r: u8, loc: ModLocation| -> &'static Modification {
            Box::leak(Box::new(Modification {
                name: name.into(),
                mass_delta: delta,
                residue: ResidueSpec::Specific(r),
                location: loc,
                fixed: false,
                accession: None,
                neutral_losses: Vec::new(),
                loss_class: 0,
                annotated_site_mass: None,
                blocks_cleavage: false,
                diagnostic_ions: Vec::new(),
                parts: Vec::new(),
            }))
        };
        let lys8 = spec.channels[1].label_at(Site::Residue(b'K')).unwrap();
        let arg10 = spec.channels[1].label_at(Site::Residue(b'R')).unwrap();
        let acetyl = m("Acetyl", 42.010565, b'K', ModLocation::ProtNTerm);
        let heavy_k = m(
            &format!("{}+Acetyl", lys8.name),
            lys8.delta + 42.010565,
            b'K',
            ModLocation::ProtNTerm,
        );
        let heavy_r = m(&arg10.name, arg10.delta, b'R', ModLocation::Anywhere);
        let cand = |first: &'static Modification, last: Option<&'static Modification>| {
            let mut res: Vec<AminoAcid> = b"KPEPTIDER"
                .iter()
                .map(|&r| AminoAcid::standard(r).unwrap().clone())
                .collect();
            res[0] = res[0].clone().with_mod(first);
            if let Some(l) = last {
                res[8] = res[8].clone().with_mod(l);
            }
            Candidate {
                peptide: Peptide::new(res, b'_', b'A'),
                protein_index: 0,
                start_offset_in_protein: 0,
                is_decoy: false,
                is_protein_n_term: true,
                is_protein_c_term: false,
            }
        };
        let light = labelled(&cand(acetyl, None), &spec).unwrap();
        let heavy = labelled(&cand(heavy_k, Some(heavy_r)), &spec).unwrap();
        assert_eq!(light.channel, 0);
        assert_eq!(heavy.channel, 1);
        let sorted = |v: &[Site]| {
            let mut v = v.to_vec();
            v.sort_by_key(|s| format!("{s:?}"));
            v
        };
        assert_eq!(sorted(&light.instances), sorted(&heavy.instances));
        assert_eq!(light.instances.len(), 2);
        assert_eq!(light.key, heavy.key);
        assert!((light.light_mass - heavy.light_mass).abs() < 1e-6);
    }

    #[test]
    fn protein_ratio_is_the_median_of_peptide_medians() {
        let row = |pep: &str, z: u8, v: f64, status: &str| output::LabelRow {
            run: "r".into(),
            peptide: pep.into(),
            sequence: pep.into(),
            charge: z,
            proteins: vec!["P1".into()],
            channel: 2,
            channel_name: "heavy".into(),
            label: None,
            status: status.into(),
            spec_id: None,
            psm_q_value: None,
            rt_psm: None,
            labelled_peptidoform: None,
            mz: 0.0,
            rt_apex: 0.0,
            rt_shift: 0.0,
            width: 0.0,
            quantity: 0.0,
            quantity_se: 0.0,
            quantity_upper: 0.0,
            llr: 0.0,
            envelope_cos: 0.0,
            coherence: 0.0,
            m1_ratio: 0.0,
            mass_err_ppm: Some(0.0),
            channel_q_value: None,
            group_llr: 0.0,
            transfer: false,
            transfer_q_value: None,
            log2_ratio: Some(v),
            log2_ratio_se: None,
            log2_ratio_normalized: Some(v - 0.5),
            converted_fraction: None,
            trace_quantity: 0.0,
        };
        // Peptide A: charges 2 and 3 (median 1.0); B: 3.0; C: -1.0.
        let rows = vec![
            row("A", 2, 0.5, "identified"),
            row("A", 3, 1.5, "identified"),
            row("B", 2, 3.0, "match_channel"),
            row("C", 2, -1.0, "identified"),
        ];
        let r = protein_ratios(&rows, |_| None, 2);
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].ratio_count, 3);
        assert_eq!(r[0].inferred_count, 1);
        assert!((r[0].log2_ratio - 1.0).abs() < 1e-12);
        assert!((r[0].log2_ratio_normalized.unwrap() - 0.5).abs() < 1e-12);
        // One peptide is not enough.
        assert!(protein_ratios(&rows[3..], |_| None, 2).is_empty());
        let tsv = protein_ratio_tsv(&r);
        assert!(tsv
            .lines()
            .nth(1)
            .unwrap()
            .starts_with("P1\tP1\tr\t2\theavy\t1.000000"));
    }

    #[test]
    fn max_sites_follow_missed_cleavages() {
        use model::enzyme::Enzyme;
        // Seven L and five R in all; at most four L in 6 residues.
        let db = model::ProteinDb {
            proteins: vec![model::Protein {
                accession: "P1".into(),
                description: String::new(),
                sequence: b"GGLRLRLLRAGLRGLRLGGGGGGGGGGG".to_vec(),
                annotation: None,
            }],
        };
        let s = LabelSpec::parse("silac2").unwrap();
        assert_eq!(max_label_sites(&s, &[Enzyme::Trypsin], 2, 2, &db, 50), 3);
        // Semi-specific spans keep the internal missed-cleavage bound.
        assert_eq!(max_label_sites(&s, &[Enzyme::Trypsin], 1, 2, &db, 50), 3);
        let d = LabelSpec::parse("dimethyl3").unwrap();
        assert_eq!(max_label_sites(&d, &[Enzyme::Trypsin], 2, 1, &db, 50), 3);
        // Leu is no cleavage site: its labels are bound by the database, the
        // most L in any peptide-length window.
        let l = LabelSpec::parse("[][Leu3]").unwrap();
        assert_eq!(max_label_sites(&l, &[Enzyme::Trypsin], 2, 2, &db, 50), 7);
        assert_eq!(max_label_sites(&l, &[Enzyme::Trypsin], 2, 2, &db, 6), 4);
        // Arg labels with Lys-C are not bound by its cleavages either.
        assert_eq!(max_label_sites(&s, &[Enzyme::LysC], 2, 2, &db, 50), 3 + 5);
        // A non-specific digest has no missed cleavages: every K/R counts.
        assert_eq!(max_label_sites(&s, &[Enzyme::Trypsin], 0, 2, &db, 50), 5);
    }

    #[test]
    fn search_mods_name_channels() {
        let d = LabelSpec::parse("dimethyl3").unwrap();
        let (mods, channels) = search_label_mods(&d);
        assert_eq!(mods.len(), 6);
        assert_eq!(channels.len(), 3);
        assert_eq!(channels[0], ("Dimethyl".to_string(), 1));
        let s = LabelSpec::parse("silac2").unwrap();
        let (mods, channels) = search_label_mods(&s);
        assert_eq!(mods.len(), 2);
        assert!(channels.iter().all(|c| c.1 == 2));
    }
}
