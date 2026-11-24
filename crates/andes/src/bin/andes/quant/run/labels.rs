//! Labelled-pair quantification of a run (SILAC, dimethyl and other mass-shift labels): the
//! channels of each identified group, their joint fit, transfers between runs and protein ratios.

mod fit;
mod hits;
mod proteins;
mod settings;

use std::collections::HashMap;

use model::mass::PROTON;
use model::spectrum::Spectrum;
use model::stats::{finite_vec, median_mad, Even, MAD_TO_SD};
use quant::formula::Formula;
use quant::label::{
    channel_score, classify_channel, classify_identified, classify_transferred, coelution_factor,
    log2_ratio, mass_factor, median, weighted_qvalues, ChannelFit, ChannelStatus, ClassifyParams,
    FitParams, GroupFit, LabelSpec, Site,
};
use rayon::prelude::*;
use search::candidates::Candidate;

use super::{faims_groups, group_of, peptidoform_of, Hit, QuantCollector};
use fit::{
    decoy_probes, fit_at_identified_profile, fit_folded, group_score, learn_shifts, learn_width,
    probe_of, seed_of, Fitted, Scales, MIN_NORM_RATIOS, MIN_WIDTH_FITS, MIN_WIDTH_LLR, N_DECOYS,
    N_TRANSFER_DECOYS,
};
use hits::{
    heavy_proline_delta, labelled, satellite_on_another_channel, ChannelHit, Group, RunSeries,
};
pub(super) use proteins::{protein_ratio_tsv, protein_ratios, MIN_RATIO_COUNT};
pub(crate) use settings::{max_label_sites, search_label_mods, LabelSettings};

type Metadata = Vec<(String, String)>;

pub(super) type LabelOutputs = (
    Vec<output::LabelRow>,
    Vec<output::FeatureRecord>,
    Vec<(String, String)>,
);

#[allow(clippy::too_many_arguments)]
pub(super) fn quantify(
    collector: &QuantCollector,
    settings: &LabelSettings,
    spectra: &[Spectrum],
    target_hits: &[Hit<'_>],
    decoy_hits: &[Hit<'_>],
    candidates: &[Candidate],
    index: &search::index::SearchIndex,
    with_records: bool,
    next_feature_id: &mut i64,
) -> Result<LabelOutputs, String> {
    let all_hits: Vec<&Hit<'_>> = target_hits.iter().chain(decoy_hits).collect();
    let hits = &all_hits[..];
    let spec = &settings.spec;
    let n_ch = spec.len();
    let tol = settings.tol_ppm;

    let (groups, ids, best_of_group) =
        label_groups(hits, settings, spectra, candidates, index, collector);

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
    let mut run_scales: Vec<Vec<Option<Scales>>> = Vec::with_capacity(n_runs);
    let mut run_series_cv: Vec<Vec<Option<f32>>> = Vec::with_capacity(n_runs);

    let keep_resident = settings.mbr_fdr.is_some() && n_runs >= 2;
    let mut resident: Option<(usize, RunSeries)> = None;
    let run_inputs = RunFitInputs {
        hits,
        groups: &groups,
        ids: &ids,
        spectra,
        spec,
        n_ch,
        tol,
        fit_params: &fit_params,
        keep_resident,
    };
    for (run_idx, run) in collector.runs.iter().enumerate() {
        fit_run(
            &run_inputs,
            run_idx,
            run,
            &mut fitted,
            &mut run_scales,
            &mut run_series_cv,
            &mut resident,
        )?;
    }

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
            let (m, mad) = median_mad(&finite_vec(&v), Even::Mean).unwrap_or((f64::NAN, f64::NAN));
            (m, (MAD_TO_SD * mad).max(0.5))
        })
        .collect();
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
            let mad = median_mad(&finite_vec(&v), Even::Mean).map_or(f64::NAN, |(_, mad)| mad);
            (MAD_TO_SD * mad).max(0.25 * dt)
        })
        .collect();
    for (r, (m, sd)) in calib.iter().enumerate() {
        log::info!(
            target: "label",
            "{}: mass error of the identified channels {m:+.2} ppm, spread {sd:.2} \
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
                let sd = sigma * std::f64::consts::SQRT_2;
                let ch = &f.fit.channels[c];
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

    let mut order: Vec<usize> = (0..fitted.len()).collect();
    order.sort_by_key(|&i| (fitted[i].run, fitted[i].gid));
    let any_transfer = fitted.iter().any(|f| f.transfer_q.is_some());
    let emit = EmitInputs {
        collector,
        fitted: &fitted,
        groups: &groups,
        statuses: &statuses,
        hits,
        spectra,
        candidates,
        index,
        spec,
        channel_q: &channel_q,
        offsets: &offsets,
        best_of_group: &best_of_group,
        n_ch,
        r,
        with_records,
        any_transfer,
    };
    let mut out = Emitted {
        rows: Vec::new(),
        records: Vec::new(),
        counts: vec![[0usize; 5]; n_ch],
    };
    for &fi in &order {
        emit_fitted(&emit, fi, &mut out, next_feature_id);
    }
    let Emitted {
        rows,
        records,
        counts,
    } = out;
    for (c, k) in counts.iter().enumerate() {
        log::info!(
            target: "label",
            "channel {} ({}): {} identified, {} matched between channels, {} \
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
        log::debug!(
            target: "label",
            "{} channel {} median log2 ratio to channel {}: {o:+.3}",
            collector.runs[*run_idx].stem,
            c + 1,
            r + 1
        );
    }
    Ok((rows, records, meta))
}

/// The identified hits of a label group in one run, per channel.
type ChannelIds = HashMap<(usize, usize), Vec<Option<ChannelHit>>>;

/// Whether `h` is a better identification than `o`: lower q-value, else higher rank score.
fn better(h: &Hit<'_>, o: &Hit<'_>) -> bool {
    match (h.q_value, o.q_value) {
        (Some(a), Some(b)) if a != b => a < b,
        _ => h.psm.rank_score > o.psm.rank_score,
    }
}

/// The label groups of the hits (one per labelled peptide key and charge), the best identified
/// hit of each group per run and channel, and the best hit of each group overall.
fn label_groups(
    hits: &[&Hit<'_>],
    settings: &LabelSettings,
    spectra: &[Spectrum],
    candidates: &[Candidate],
    index: &search::index::SearchIndex,
    collector: &QuantCollector,
) -> (Vec<Group>, ChannelIds, HashMap<usize, usize>) {
    let spec = &settings.spec;
    let n_ch = spec.len();
    let sites = spec.sites();
    let tol = settings.tol_ppm;
    let mut groups: Vec<Group> = Vec::new();
    let mut gid_of: HashMap<(String, u8), usize> = HashMap::new();
    let mut ids: HashMap<(usize, usize), Vec<Option<ChannelHit>>> = HashMap::new();
    let mut best_of_group: HashMap<usize, usize> = HashMap::new();
    let (mut mixed, mut no_rt) = (0usize, 0usize);
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
    log::info!(
        target: "label",
        "{} label groups from {} PSMs ({} without a consistent label state, \
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
        log::warn!(
            target: "label",
            "{shared} groups have an Arg->Pro satellite within the MS1 tolerance \
             of another channel's isotope; the split between the two (converted_fraction, and \
             the ratio through it) is weakly determined there"
        );
    }
    (groups, ids, best_of_group)
}

/// What the per-run label fit reads: the hits and their label groups, the identified hits per
/// group and run, the spectra, the label spec, the tolerance and the fit parameters.
#[derive(Clone, Copy)]
struct RunFitInputs<'a> {
    hits: &'a [&'a Hit<'a>],
    groups: &'a [Group],
    ids: &'a ChannelIds,
    spectra: &'a [Spectrum],
    spec: &'a LabelSpec,
    n_ch: usize,
    tol: f64,
    fit_params: &'a FitParams,
    keep_resident: bool,
}

/// Fits every label group identified in one run, per FAIMS series: learns the series' elution
/// width and channel RT shifts, fits each group's channels and decoy probes, and appends the fits
/// and the series scales. With `keep_resident`, the run's MS1 series stay loaded for transfers.
fn fit_run(
    inputs: &RunFitInputs<'_>,
    run_idx: usize,
    run: &super::collect::RunCapture,
    fitted: &mut Vec<Fitted>,
    run_scales: &mut Vec<Vec<Option<Scales>>>,
    run_series_cv: &mut Vec<Vec<Option<f32>>>,
    resident: &mut Option<(usize, RunSeries)>,
) -> Result<(), String> {
    let RunFitInputs {
        hits,
        groups,
        ids,
        spectra,
        spec,
        n_ch,
        tol,
        fit_params,
        keep_resident,
    } = *inputs;
    let Some(cache) = run.ms1.as_ref() else {
        run_scales.push(Vec::new());
        run_series_cv.push(Vec::new());
        return Ok(());
    };
    *resident = None;
    let ms1 = cache
        .load()
        .map_err(|e| format!("load MS1 scans for {}: {e}", run.stem))?;
    let series = faims_groups(ms1);
    let mut scales_here: Vec<Option<Scales>> = Vec::new();
    run_series_cv.push(series.iter().map(|(cv, _)| *cv).collect());
    let mut by_series: Vec<Vec<(usize, Vec<Option<ChannelHit>>)>> = vec![Vec::new(); series.len()];
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
    claims.sort_by(|a, b| a.2.total_cmp(&b.2));
    for (si, (_, ms1)) in series.iter().enumerate() {
        let list = &by_series[si];
        if list.is_empty() {
            scales_here.push(None);
            continue;
        }
        let dt = ms1.median_scan_spacing().unwrap_or(1.0).max(1e-3);
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
        let w = learn_width(ms1, groups, &anchors, tol);
        let shifts = learn_shifts(ms1, groups, &anchors, (dt, w), spec, tol);
        log::info!(
            target: "label",
            "{}{}: {} identified groups, elution width {:.2} s, MS1 spacing \
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
                    (fit_folded(ms1, g, &probe, fit_params)?, probe)
                } else {
                    let ((mut fit, conv), pinned) =
                        fit_at_identified_profile(ms1, g, &probe, chans, fit_params)?;
                    if let Some((free, _)) = fit_folded(ms1, g, &probe, fit_params) {
                        let hits = chans.iter();
                        for ((ch, f), hit) in fit.channels.iter_mut().zip(&free.channels).zip(hits)
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
                    for dp in decoy_probes(g, &probe, &missing, tol, seed, N_DECOYS) {
                        let df = fit_folded(ms1, g, &dp, fit_params).map(|x| x.0);
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
        *resident = Some((run_idx, series));
    }
    Ok(())
}

/// What writing the output rows of a fit reads: the fits, their groups and channel statuses, the
/// hits and spectra behind them, the channel q-values and normalisation offsets, the reference
/// channel and whether feature records are written.
#[derive(Clone, Copy)]
struct EmitInputs<'a> {
    collector: &'a QuantCollector,
    fitted: &'a [Fitted],
    groups: &'a [Group],
    statuses: &'a [Vec<ChannelStatus>],
    hits: &'a [&'a Hit<'a>],
    spectra: &'a [Spectrum],
    candidates: &'a [Candidate],
    index: &'a search::index::SearchIndex,
    spec: &'a LabelSpec,
    channel_q: &'a HashMap<(usize, usize), f64>,
    offsets: &'a HashMap<(usize, usize), f64>,
    best_of_group: &'a HashMap<usize, usize>,
    n_ch: usize,
    r: usize,
    with_records: bool,
    any_transfer: bool,
}

/// The rows, feature records and per-channel status counts written so far.
struct Emitted {
    rows: Vec<output::LabelRow>,
    records: Vec<output::FeatureRecord>,
    counts: Vec<[usize; 5]>,
}

/// Appends one fit's label rows (one per channel), its status counts and, with records, its
/// feature record. Identity-decoy groups write nothing.
fn emit_fitted(inputs: &EmitInputs<'_>, fi: usize, out: &mut Emitted, next_feature_id: &mut i64) {
    let EmitInputs {
        collector,
        fitted,
        groups,
        statuses,
        hits,
        spectra,
        candidates,
        index,
        spec,
        channel_q,
        offsets,
        best_of_group,
        n_ch,
        r,
        with_records,
        any_transfer,
    } = *inputs;
    let Emitted {
        rows,
        records,
        counts,
    } = out;
    let f = &fitted[fi];
    let g = &groups[f.gid];
    if g.identity_decoy {
        return;
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
            log2_ratio_normalized: ratio.and_then(|x| offsets.get(&(f.run, c)).map(|o| x.0 - o)),
            converted_fraction: f.converted[c],
            trace_quantity: ch.trace_quantity,
        });
    }
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
                let c = labelled(hit.cand, spec).map_or(0, |l| l.channel);
                let calc = g.channel_mz[c];
                let err = f.fit.channels[c].mass_err_ppm;
                rec.calculated_mz = calc as f32;
                rec.observed_mz = (calc * (1.0 + err.unwrap_or(0.0) * 1e-6)) as f32;
                rec.mass_error_ppm = err.map(|e| e as f32);
            }
            if any_transfer {
                rec.mbr = Some(output::features::MbrColumns {
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
    let score = |fit: &GroupFit, run: usize| -> f64 {
        let (center, sigma) = calib[run];
        let err = fit
            .channels
            .iter()
            .max_by(|a, b| a.amplitude.total_cmp(&b.amplitude))
            .and_then(|c| c.mass_err_ppm);
        group_score(fit)
            * quant::label::mass_factor(err, center, sigma)
            * (-0.5 * fit.rt_z * fit.rt_z).exp()
    };
    let n_runs = collector.runs.len();
    let n_ch = groups.first().map_or(0, |g| g.channel_mz.len());
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
    alignment.add_observations(&observe(true))?;
    let mut present = vec![vec![false; n_runs]; groups.len()];
    for f in fitted {
        present[f.gid][f.run] = true;
    }
    let mut out: Vec<Fitted> = Vec::new();
    let mut scored: Vec<(f64, f64)> = Vec::new();
    let mut target_rows: Vec<usize> = Vec::new();
    let mut identity_scores: Vec<f64> = Vec::new();
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
            log::warn!(
                target: "label",
                "{}: not aligned (too few shared identified groups); no transfers \
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
    let identity_rate = accepted_identity as f64 / accepted_targets.max(1) as f64;
    let mut kept = Vec::new();
    for (k, mut f) in out.into_iter().enumerate() {
        let qk = (q[target_rows[k]] + identity_rate).min(1.0);
        if qk <= fdr {
            f.transfer_q = Some(qk);
            kept.push(f);
        }
    }
    log::info!(
        target: "label",
        "match between runs: {} of {} transfers at transfer q<={fdr} (lowest \
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

#[cfg(test)]
#[path = "labels/tests/labels.rs"]
mod tests;
