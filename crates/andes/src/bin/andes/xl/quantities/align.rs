//! Retention-time alignment of the runs of a cross-link experiment and the transfers it
//! allows. The anchors are the cross-link features and the linear peptides identified in the
//! same spectra: cross-links alone are too few to align a run, and the linear peptides of
//! the same injection elute on the same gradient. Every anchor enters with the apex of its
//! MS1 feature, not the time of its MS2 spectrum, which is taken anywhere along the elution.
//! `quant::align` fits the runs jointly.

use std::collections::{BTreeMap, HashSet};

use model::mass::PROTON;
use quant::align::{align, AlignParams, Alignment, Observation};
use search::tdc::{qvalues, ScoredLabel};

use super::library::{peptide_formula, traced_envelope, Library};
use super::measure::Transfer;
use super::{accepted, Evidence, Row, Search, Settings, ENVELOPE_ISOTOPES};
use crate::quant::run::RunMs1;

/// A linear peptide is an anchor when its best match in a run is a target at or below this
/// q-value of the run's target-decoy competition.
const LINEAR_ANCHOR_Q: f64 = 0.01;
/// Spectra searched for linear anchors at a time.
const LINEAR_CHUNK: usize = 4096;

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(super) struct Counts {
    pub(super) aligned_runs: usize,
    pub(super) anchors: usize,
    pub(super) linear_anchors: usize,
    pub(super) xl_anchors: usize,
    pub(super) attempted: usize,
    pub(super) measured: usize,
}

pub(super) struct Aligned {
    alignment: Alignment,
    pub(super) counts: Counts,
}

/// A linear peptide ion: its peptidoform (with the target or decoy bit) and charge.
pub(super) type LinearKey = (u64, u8);

/// A confident linear peptide ion of one run.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct LinearAnchor {
    pub(super) key: LinearKey,
    /// The retention time of its best match, in seconds.
    pub(super) rt: f64,
    pub(super) mono_mz: f64,
    pub(super) charge: u8,
    /// The traced isotopes, summing to one.
    pub(super) envelope: Vec<f64>,
}

/// The confident linear peptides of one run with the retention time of their best match.
pub(super) fn linear_peptides(search: &Search<'_>, run: &RunMs1<'_>) -> Vec<LinearAnchor> {
    let prepared = search.prepared;
    let spectra = search.spectra.get(run.span.clone()).unwrap_or(&[]);
    let mut tops: Vec<(LinearKey, f64, u32)> = Vec::new();
    let mut labels: Vec<ScoredLabel> = Vec::new();
    for (c, chunk) in spectra.chunks(LINEAR_CHUNK).enumerate() {
        let offset = run.span.start + c * LINEAR_CHUNK;
        let queues = prepared.run_chunk_for_top(chunk, offset, prepared.params);
        for (spec, queue) in chunk.iter().zip(queues) {
            let (Some(top), Some(rt)) = (queue.peek_top(), spec.rt_seconds) else {
                continue;
            };
            if !rt.is_finite() {
                continue;
            }
            tops.push((
                (top.tie_key, top.charge_used),
                rt,
                top.primary_candidate_idx(),
            ));
            labels.push(ScoredLabel {
                score: top.score,
                is_decoy: prepared.candidate_is_decoy(top.primary_candidate_idx()),
            });
        }
    }
    let q = qvalues(&labels);
    let mut best: BTreeMap<LinearKey, (f32, f64, u32)> = BTreeMap::new();
    for (((key, rt, candidate), label), q) in tops.into_iter().zip(&labels).zip(q) {
        if label.is_decoy || q > LINEAR_ANCHOR_Q {
            continue;
        }
        best.entry(key)
            .and_modify(|kept| {
                if label.score > kept.0 || (label.score == kept.0 && rt < kept.1) {
                    *kept = (label.score, rt, candidate);
                }
            })
            .or_insert((label.score, rt, candidate));
    }
    best.into_iter()
        .map(|(key, (_, rt, candidate))| {
            let peptide = prepared.candidate_peptide(candidate);
            let charge = key.1.max(1);
            let (envelope, _) =
                traced_envelope(&peptide_formula(&peptide).envelope(ENVELOPE_ISOTOPES));
            LinearAnchor {
                key,
                rt,
                mono_mz: peptide.mass() / f64::from(charge) + PROTON,
                charge,
                envelope,
            }
        })
        .collect()
}

/// Aligns the runs on the accepted cross-link features of `rows` and on the features of the
/// linear peptides (`linear`, per run) that at least two runs share.
pub(super) fn build(
    lib: &Library,
    rows: &[Row],
    linear: &[BTreeMap<LinearKey, f64>],
    n_runs: usize,
    settings: &Settings,
) -> Aligned {
    let t0 = std::time::Instant::now();
    let n_xl = lib.precursors.len();
    let mut obs: Vec<Observation> = Vec::new();
    let mut xl_runs = vec![0usize; n_xl];
    for row in rows {
        let Some(feature) = &row.feature else {
            continue;
        };
        if row.evidence != Evidence::Identified || !accepted(row, settings) {
            continue;
        }
        xl_runs[row.precursor] += 1;
        obs.push(Observation {
            precursor: row.precursor as u32,
            run: row.run as u16,
            rt: feature.apex_rt,
        });
    }
    let per_run = linear;
    let mut seen: BTreeMap<LinearKey, usize> = BTreeMap::new();
    for peptides in per_run {
        for key in peptides.keys() {
            *seen.entry(*key).or_insert(0) += 1;
        }
    }
    let shared: BTreeMap<LinearKey, u32> = seen
        .into_iter()
        .filter(|(_, n)| *n >= 2)
        .enumerate()
        .map(|(i, (key, _))| (key, (n_xl + i) as u32))
        .collect();
    for (run, peptides) in per_run.iter().enumerate() {
        for (key, &rt) in peptides {
            if let Some(&id) = shared.get(key) {
                obs.push(Observation {
                    precursor: id,
                    run: run as u16,
                    rt,
                });
            }
        }
    }
    let alignment = align(&obs, n_xl + shared.len(), n_runs, &AlignParams::default());
    let xl_anchors = xl_runs.iter().filter(|n| **n >= 2).count();
    let counts = Counts {
        aligned_runs: alignment.maps.iter().flatten().count(),
        anchors: shared.len() + xl_anchors,
        linear_anchors: shared.len(),
        xl_anchors,
        ..Counts::default()
    };
    log::debug!(
        target: "xl-mbr",
        "linear anchors with a feature per run: {:?}; aligned in {:.2}s",
        per_run.iter().map(BTreeMap::len).collect::<Vec<_>>(),
        t0.elapsed().as_secs_f64()
    );
    Aligned { alignment, counts }
}

/// Every precursor in every run that did not identify it and where the alignment predicts a
/// retention time for it.
pub(super) fn transfers(
    lib: &Library,
    rows: &[Row],
    aligned: &Aligned,
    n_runs: usize,
) -> Vec<Transfer> {
    let identified: HashSet<(usize, usize)> = rows
        .iter()
        .filter(|r| r.evidence == Evidence::Identified)
        .map(|r| (r.precursor, r.run))
        .collect();
    let mut out = Vec::new();
    for precursor in 0..lib.precursors.len() {
        for run in 0..n_runs {
            if identified.contains(&(precursor, run)) {
                continue;
            }
            let Some((rt, sd_s)) = aligned.alignment.predict(precursor as u32, run as u16) else {
                continue;
            };
            if rt.is_finite() && sd_s.is_finite() {
                out.push(Transfer {
                    precursor,
                    run,
                    rt,
                    sd_s,
                });
            }
        }
    }
    out
}
