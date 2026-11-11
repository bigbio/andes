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
    channel_qvalues, channel_score, classify_channel, classify_transferred, fit_group,
    learn_rt_shifts, log2_ratio, ChannelEnvelope, ChannelFit, ChannelStatus, ClassifyParams,
    FitParams, GroupFit, GroupProbe, LabelSpec, Site, FREE_WIDTHS,
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
        )
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
        });
        if !channels.iter().any(|(n, _): &(String, u8)| *n == l.name) {
            channels.push((l.name.clone(), c as u8 + 1));
        }
    }
    (mods, channels)
}

/// Most labelled sites a peptide of the search can carry: one per residue
/// site up to the missed cleavages (each labelled residue of a tryptic
/// peptide is a cleavage site or its C-terminus), plus the N-terminus.
pub(crate) fn max_label_sites(spec: &LabelSpec, missed_cleavages: u32) -> u32 {
    let sites = spec.sites();
    let residue = sites.iter().any(|s| matches!(s, Site::Residue(_)));
    let nterm = sites.contains(&Site::NTerm);
    u32::from(residue) * (missed_cleavages + 1) + u32::from(nterm)
}

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
        let mut non_label_delta = 0.0;
        let mut used: Vec<Site> = Vec::new();
        if let Some(m) = aa.mod_ {
            let mut label_delta = 0.0;
            for part in m.name.split('+') {
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
                    None => non_label.push(part),
                }
            }
            label_mass += label_delta;
            non_label_delta = m.mass_delta - label_delta;
        }
        // Labelable sites left without a label.
        let nterm_free = aa
            .mod_
            .is_none_or(|m| m.location == model::modification::ModLocation::Anywhere);
        for s in &spec_sites {
            let here = match s {
                Site::Residue(r) => *r == aa.residue && non_label.is_empty(),
                Site::NTerm => i == 0 && nterm_free && non_label.is_empty(),
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
    /// Per channel, the decoy fits' scores (missing channels).
    decoy_scores: Vec<Vec<f64>>,
    /// Per channel, another identified peptide of the run sits on it.
    claimed: Vec<bool>,
    transfer_q: Option<f64>,
    rt_shift: Vec<f64>,
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
/// Decoys per target of the channel and transfer q-values.
const N_DECOYS: usize = 4;

fn probe_of(g: &Group, t_hat: f64, s: f64, w: f64, shifts: &[f64]) -> GroupProbe {
    GroupProbe {
        charge: g.charge,
        channel_mz: g.channel_mz.clone(),
        envelopes: g.envelopes.clone(),
        t_hat,
        s,
        w_prior: w,
        rt_shift: shifts.to_vec(),
        free: vec![true; g.channel_mz.len()],
        mobility: None,
    }
}

/// Learns the run's elution width from `anchors` (group, RT of its best PSM
/// in reference time): rounds of fits with the width free, the median
/// fitted width of the strong fits, until it moves by less than 10 %.
fn learn_width(ms1: &Ms1RunIndex, groups: &[Group], anchors: &[(usize, f64)], tol: f64) -> f64 {
    let dt = ms1.median_scan_spacing().unwrap_or(1.0).max(1e-3);
    let step = (anchors.len() / MAX_WIDTH_FITS).max(1);
    let sample: Vec<&(usize, f64)> = anchors.iter().step_by(step).collect();
    let params = FitParams {
        tol_ppm: tol,
        widths: FREE_WIDTHS,
    };
    let mut w = 6.0f64.max(dt);
    for _ in 0..6 {
        let s = w.max(dt);
        let fitted: Vec<f64> = sample
            .par_iter()
            .filter_map(|(gid, t)| {
                let g = &groups[*gid];
                let p = probe_of(g, *t, s, w, &vec![0.0; g.channel_mz.len()]);
                fit_group(ms1, &p, &params)
            })
            .filter(|f| f.llr >= MIN_WIDTH_LLR)
            .map(|f| f.w)
            .collect();
        if fitted.len() < MIN_WIDTH_FITS {
            return (6.0f64).max(dt);
        }
        let next = quant::label::median(fitted).max(dt);
        let moved = (next - w).abs() / w;
        w = next;
        if moved < 0.1 {
            break;
        }
    }
    w
}

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
                p.free = (0..n).map(|k| k == c).collect();
                fit_group(ms1, &p, &params).filter(|f| f.llr >= MIN_WIDTH_LLR)
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
) -> Vec<GroupProbe> {
    let z = g.charge.max(1) as f64;
    let mass = (g.channel_mz[0] - PROTON) * z;
    let min_sep = quant::label::decoy_min_sep(mass, tol_ppm);
    quant::label::decoy_shifts(moved, &g.deltas, min_sep, seed, N_DECOYS)
        .into_iter()
        .map(|shift| {
            let mut p = probe.clone();
            for &c in moved {
                p.channel_mz[c] += shift / z;
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

/// Per-run label quantification. Returns the `label.parquet` rows and the
/// QPX feature records (identified groups, quantified channels).
#[allow(clippy::too_many_arguments)]
pub(super) fn quantify(
    collector: &QuantCollector,
    settings: &LabelSettings,
    spectra: &[Spectrum],
    hits: &[Hit<'_>],
    candidates: &[Candidate],
    index: &search::search_index::SearchIndex,
    with_records: bool,
    next_feature_id: &mut i64,
) -> Result<(Vec<output::LabelRow>, Vec<output::FeatureRecord>), String> {
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
        if h.cand.is_decoy {
            continue;
        }
        let Some(lh) = labelled(h.cand, spec) else {
            mixed += 1;
            continue;
        };
        let Some(rt) = spectra[h.spec_idx].rt_seconds else {
            no_rt += 1;
            continue;
        };
        let z = h.psm.charge_used.max(1);
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
                channel_mz,
                envelopes,
                deltas,
                faims_cv: None,
            });
            groups.len() - 1
        });
        let slot = ids
            .entry((gid, h.run_idx))
            .or_insert_with(|| vec![None; n_ch]);
        let replace = slot[lh.channel].is_none_or(|o| better(h, &hits[o.hi]));
        if replace {
            slot[lh.channel] = Some(ChannelHit { hi, rt });
        }
        let b = best_of_group.entry(gid).or_insert(hi);
        if better(h, &hits[*b]) {
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

    for (run_idx, run) in collector.runs.iter().enumerate() {
        let Some(cache) = run.ms1.as_ref() else {
            run_scales.push(Vec::new());
            run_series_cv.push(Vec::new());
            continue;
        };
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
                .reduce(|a, b| {
                    if better(&hits[b.hi], &hits[a.hi]) {
                        b
                    } else {
                        a
                    }
                })
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
            for (gid, chans) in list {
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
        for (si, (_, ms1)) in series.iter().enumerate() {
            let list = &by_series[si];
            if list.is_empty() {
                scales_here.push(None);
                continue;
            }
            let dt = ms1.median_scan_spacing().unwrap_or(1.0).max(1e-3);
            // Anchors in the reference channel's time (no shifts yet).
            let anchors: Vec<(usize, f64)> = list
                .iter()
                .map(|(gid, chans)| {
                    let best = chans
                        .iter()
                        .flatten()
                        .copied()
                        .reduce(|a, b| {
                            if better(&hits[b.hi], &hits[a.hi]) {
                                b
                            } else {
                                a
                            }
                        })
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
                            if better(&hits[b.1.hi], &hits[a.1.hi]) {
                                b
                            } else {
                                a
                            }
                        })?;
                    let t_hat = best.1.rt - shifts[best.0];
                    let probe = probe_of(g, t_hat, s, w, &shifts);
                    let fit = fit_group(ms1, &probe, &fit_params)?;
                    let missing: Vec<usize> = (0..n_ch).filter(|&c| chans[c].is_none()).collect();
                    let mut decoy_scores = vec![Vec::new(); n_ch];
                    if !missing.is_empty() {
                        let seed = seed_of(&g.key, run_idx);
                        for dp in decoy_probes(g, &probe, &missing, tol, seed) {
                            let df = fit_group(ms1, &dp, &fit_params);
                            for &c in &missing {
                                decoy_scores[c].push(
                                    df.as_ref().map_or(0.0, |d| channel_score(&d.channels[c])),
                                );
                            }
                        }
                    }
                    let claimed = (0..n_ch)
                        .map(|c| {
                            chans[c].is_none() && {
                                let mz = g.channel_mz[c];
                                let tolda = mz * tol * 1e-6;
                                let t = fit.t0 + shifts[c];
                                claims.iter().any(|&(csi, z, cmz, crt, cg)| {
                                    csi == si
                                        && z == g.charge
                                        && cg != *gid
                                        && (cmz - mz).abs() <= tolda
                                        && (crt - t).abs() <= 3.0 * fit.w + dt
                                })
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
                    })
                })
                .collect();
            fitted.extend(rows);
            scales_here.push(Some(Scales { dt, w, shifts }));
        }
        run_scales.push(scales_here);
    }

    // ── Match between runs ────────────────────────────────────────────────
    if let (Some(mbr_fdr), true) = (settings.mbr_fdr, n_runs >= 2) {
        let transfers = transfer(
            collector,
            &groups,
            &fitted,
            &run_scales,
            &run_series_cv,
            &fit_params,
            mbr_fdr,
        )?;
        fitted.extend(transfers);
    }

    // ── Channel q-values (match between channels), pooled over runs ───────
    // (fitted index, channel) of every unidentified channel of an identified
    // group, its score against its decoys'.
    let mut channel_q: HashMap<(usize, usize), f64> = HashMap::new();
    {
        let mut rows: Vec<(f64, bool)> = Vec::new();
        let mut at: Vec<((usize, usize), usize)> = Vec::new();
        for (fi, f) in fitted.iter().enumerate() {
            if f.transfer_q.is_some() {
                continue;
            }
            for c in 0..n_ch {
                if f.ids[c].is_some() || f.decoy_scores[c].is_empty() {
                    continue;
                }
                at.push(((fi, c), rows.len()));
                rows.push((channel_score(&f.fit.channels[c]), false));
                rows.extend(f.decoy_scores[c].iter().map(|&d| (d, true)));
            }
        }
        let q = channel_qvalues(&rows, N_DECOYS);
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
                    ChannelStatus::Identified
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
    for &fi in &order {
        let f = &fitted[fi];
        let g = &groups[f.gid];
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
            });
        }
        // QPX: one feature per identified (group, run), quantified channels.
        if with_records && f.transfer_q.is_none() {
            let best = f.ids.iter().flatten().copied().reduce(|a, b| {
                if better(&hits[b.hi], &hits[a.hi]) {
                    b
                } else {
                    a
                }
            });
            if let Some(best) = best {
                let hit = &hits[best.hi];
                let spec_ = &spectra[hit.spec_idx];
                let id = super::quant_id(hit, run, spec_, candidates, index);
                let proteins = super::proteins_of(&hit.psm, candidates, index);
                let mut rec = super::base_record(&id, proteins, *next_feature_id);
                *next_feature_id += 1;
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
    Ok((rows, records))
}

/// Match between runs: every group of the library (identified in some run)
/// fitted in each run that lacks it, at its aligned apex; the transfers at
/// the q-value cut.
fn transfer(
    collector: &QuantCollector,
    groups: &[Group],
    fitted: &[Fitted],
    run_scales: &[Vec<Option<Scales>>],
    run_series_cv: &[Vec<Option<f32>>],
    fit_params: &FitParams,
    fdr: f64,
) -> Result<Vec<Fitted>, String> {
    let n_runs = collector.runs.len();
    let n_ch = groups.first().map_or(0, |g| g.channel_mz.len());
    // Observations: the fitted apex (reference time) of strong identified fits.
    let obs: Vec<quant::align::Observation> = fitted
        .iter()
        .filter(|f| f.fit.llr >= MIN_WIDTH_LLR)
        .map(|f| quant::align::Observation {
            precursor: f.gid as u32,
            run: f.run as u16,
            rt: f.fit.t0,
        })
        .collect();
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
    let alignment = quant::align::align(&obs, groups.len(), n_runs, &params);
    let mut present = vec![vec![false; n_runs]; groups.len()];
    for f in fitted {
        present[f.gid][f.run] = true;
    }
    let mut out: Vec<Fitted> = Vec::new();
    let mut scored: Vec<(f64, bool)> = Vec::new();
    let mut target_rows: Vec<usize> = Vec::new();
    for (run_idx, run) in collector.runs.iter().enumerate() {
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
        let ms1 = cache
            .load()
            .map_err(|e| format!("load MS1 scans for {}: {e}", run.stem))?;
        let series = faims_groups(ms1);
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
                let fit = fit_group(ms1, &probe, fit_params)?;
                let all: Vec<usize> = (0..n_ch).collect();
                let dscores: Vec<f64> = decoy_probes(
                    g,
                    &probe,
                    &all,
                    fit_params.tol_ppm,
                    seed_of(&g.key, run_idx),
                )
                .iter()
                .map(|d| fit_group(ms1, d, fit_params).map_or(0.0, |d| group_score(&d)))
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
                    },
                    dscores,
                ))
            })
            .collect();
        for (f, d) in results {
            target_rows.push(scored.len());
            scored.push((group_score(&f.fit), false));
            scored.extend(d.into_iter().map(|x| (x, true)));
            out.push(f);
        }
    }
    let q = channel_qvalues(&scored, N_DECOYS);
    let mut kept = Vec::new();
    for (k, mut f) in out.into_iter().enumerate() {
        let qk = q[target_rows[k]];
        if qk <= fdr {
            f.transfer_q = Some(qk);
            kept.push(f);
        }
    }
    eprintln!(
        "quant: label: match between runs: {} of {} transfers at transfer q<={fdr}",
        kept.len(),
        target_rows.len()
    );
    Ok(kept)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn max_sites_follow_missed_cleavages() {
        let s = LabelSpec::parse("silac2").unwrap();
        assert_eq!(max_label_sites(&s, 2), 3);
        let d = LabelSpec::parse("dimethyl3").unwrap();
        assert_eq!(max_label_sites(&d, 1), 3);
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
