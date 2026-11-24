//! Label-free targets: FAIMS groups, the target of each identification, their quantification,
//! apex competition and isotope envelopes.

use model::mass::PROTON;
use model::spectrum::Spectrum;
use quant::formula::Formula;
use quant::lfq::{competition_score, quantify, quantify_decoy, FeatureQuant, LfqParams, LfqTarget};
use quant::ms1::Ms1RunIndex;
use quant::tdc::Pair;
use rayon::prelude::*;
use search::candidates::Candidate;

use super::collect::RunCapture;
use super::hits::Hit;
use super::scans::native_id;

pub(super) type Quantified = (usize, Option<FeatureQuant>, Option<FeatureQuant>);

pub(super) fn faims_groups(ms1: Ms1RunIndex) -> Vec<(Option<f32>, Ms1RunIndex)> {
    if ms1.faims_voltages().len() > 1 {
        ms1.split_by_faims()
            .into_iter()
            .map(|(cv, idx)| (Some(cv), idx))
            .collect()
    } else {
        vec![(None, ms1)]
    }
}

pub(super) fn group_of(
    groups: &[(Option<f32>, Ms1RunIndex)],
    run: &RunCapture,
    spec: &Spectrum,
) -> Option<usize> {
    match groups.first()?.0 {
        None => Some(0),
        Some(_) => {
            let cv = run
                .contexts
                .get(native_id(run, &spec.title))
                .and_then(|c| c.faims_cv)?;
            groups.iter().position(|(g, _)| *g == Some(cv))
        }
    }
}

pub(super) fn lfq_target(h: &Hit<'_>, anchor_rt: f64, n_iso: usize) -> LfqTarget {
    let z = h.psm.charge_used.max(1);
    let formula = Formula::peptide(
        h.cand
            .peptide
            .residues
            .iter()
            .map(|aa| (aa.residue, aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta))),
    );
    let (envelope, isotope_share) = formula.envelope_and_share(n_iso);
    LfqTarget {
        mono_mz: h.cand.peptide.mass() / z as f64 + PROTON,
        charge: z,
        anchor_rt,
        envelope,
        isotope_share,
    }
}

pub(super) fn quantify_targets(
    groups: &[(Option<f32>, Ms1RunIndex)],
    targets: &[(usize, usize, LfqTarget)],
    params: &LfqParams,
) -> Vec<Quantified> {
    targets
        .par_iter()
        .map(|(hi, g, t)| {
            let ms1 = &groups[*g].1;
            (
                *hi,
                quantify(ms1, t, params),
                quantify_decoy(ms1, t, params),
            )
        })
        .collect()
}

pub(super) fn max_apex_of(quantified: &[Quantified]) -> f32 {
    quantified
        .iter()
        .flat_map(|(_, f, d)| [f.as_ref(), d.as_ref()])
        .flatten()
        .map(|f| f.apex_intensity)
        .fold(0.0f32, f32::max)
}

pub(super) fn competition_pairs(
    quantified: &[Quantified],
    params: &LfqParams,
    max_apex: f32,
) -> Vec<Pair> {
    quantified
        .iter()
        .map(|(_, f, d)| Pair {
            target: f
                .as_ref()
                .and_then(|f| competition_score(f, params, max_apex)),
            decoy: d
                .as_ref()
                .and_then(|f| competition_score(f, params, max_apex)),
        })
        .collect()
}

pub(super) fn envelope_of(cand: &Candidate, n_iso: usize) -> (Vec<f32>, f32) {
    let formula = Formula::peptide(
        cand.peptide
            .residues
            .iter()
            .map(|aa| (aa.residue, aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta))),
    );
    let (env, share) = formula.envelope_and_share(n_iso);
    let total: f64 = env.iter().sum();
    let env = if total > 0.0 {
        env.iter().map(|v| (v / total) as f32).collect()
    } else {
        env.iter().map(|v| *v as f32).collect()
    };
    (env, share as f32)
}
