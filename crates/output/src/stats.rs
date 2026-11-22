use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use model::modification::{ModLocation, ResidueSpec};
use model::spectrum::Spectrum;
use model::stats::{quantile_sorted, MAD_TO_SD};
use model::tolerance::Tolerance;
use rayon::prelude::*;
use search::candidate_gen::Candidate;
use search::psm::TopNQueue;
use search::search_params::SearchParams;

pub struct RunStatistics {
    pub precursor_tolerance: String,
    pub precursor_calibration: String,
    pub fragment_tolerance: String,
    pub node_score_bin: String,
    pub spectra_with_match: usize,
    pub target_psms: usize,
    pub decoy_psms: usize,
    pub unmodified_target_psms: usize,
    pub ptm_counts: Vec<(String, usize)>,
}

fn fmt_tol(t: &Tolerance) -> String {
    match t {
        Tolerance::Ppm(v) => format!("{} ppm", round3(*v)),
        Tolerance::Da(v) => format!("{} Da", round3(*v)),
    }
}

fn round3(v: f64) -> f64 {
    (v * 1000.0).round() / 1000.0
}

fn fmt_fragment_tolerance(matching: &Tolerance, pin_features: &Tolerance) -> String {
    if pin_features == matching {
        format!("{} (retrieval and match features)", fmt_tol(matching))
    } else {
        format!(
            "{} (retrieval and model features); {} (PIN match features, measured)",
            fmt_tol(matching),
            fmt_tol(pin_features)
        )
    }
}

impl RunStatistics {
    pub fn compute(
        queues: &[TopNQueue],
        candidates: &[Candidate],
        params: &SearchParams,
        matching_tolerance: &Tolerance,
        pin_feature_tolerance: &Tolerance,
        node_score_bin: &Tolerance,
    ) -> Self {
        let mut spectra_with_match = 0usize;
        let mut target_psms = 0usize;
        let mut decoy_psms = 0usize;
        let mut unmodified_target_psms = 0usize;
        let mut ptm: BTreeMap<String, usize> = BTreeMap::new();

        for q in queues {
            if q.is_empty() {
                continue;
            }
            spectra_with_match += 1;
            let psms = q.clone().into_rank_sorted_vec();
            let Some(best) = psms.first() else {
                continue;
            };
            let cand = &candidates[best.primary_candidate_idx() as usize];
            if cand.is_decoy {
                decoy_psms += 1;
                continue;
            }
            target_psms += 1;
            let mut names: BTreeSet<&str> = BTreeSet::new();
            for aa in &cand.peptide.residues {
                if let Some(m) = aa.mod_.as_ref() {
                    names.insert(m.name.as_str());
                }
            }
            if names.is_empty() {
                unmodified_target_psms += 1;
            }
            for n in names {
                *ptm.entry(n.to_string()).or_insert(0) += 1;
            }
        }

        let mut ptm_counts: Vec<(String, usize)> = ptm.into_iter().collect();
        ptm_counts.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

        Self {
            precursor_tolerance: match &params.open_search {
                Some(o) => format!(
                    "open-search offset window {:+} … {:+} Da",
                    o.delta_min, o.delta_max
                ),
                None => format!("{:?}", params.precursor_tolerance),
            },
            precursor_calibration: format!("{:?}", params.precursor_cal_mode),
            fragment_tolerance: fmt_fragment_tolerance(matching_tolerance, pin_feature_tolerance),
            node_score_bin: fmt_tol(node_score_bin),
            spectra_with_match,
            target_psms,
            decoy_psms,
            unmodified_target_psms,
            ptm_counts,
        }
    }

    pub fn render(&self) -> String {
        self.render_with(None)
    }

    pub fn render_with(&self, summary: Option<&IdentificationSummary>) -> String {
        let label_w = self
            .ptm_counts
            .iter()
            .map(|(n, _)| n.len())
            .chain(std::iter::once("(unmodified)".len()))
            .max()
            .unwrap_or(12);

        let mut s = String::new();
        s.push_str("──────── andes run summary ────────\n");
        s.push_str(&format!(
            "  Final precursor tolerance : {} (calibration: {})\n",
            self.precursor_tolerance, self.precursor_calibration
        ));
        s.push_str(&format!(
            "  Final fragment tolerance  : {}\n",
            self.fragment_tolerance
        ));
        s.push_str(&format!(
            "  Node-score bin width      : {} (model rank tables, nominal fragment masses)\n",
            self.node_score_bin
        ));
        s.push_str(&format!(
            "  Spectra with a match      : {}\n",
            self.spectra_with_match
        ));
        s.push_str(&format!(
            "  Rank-1 PSMs (pre-FDR)     : {} target, {} decoy\n",
            self.target_psms, self.decoy_psms
        ));
        if summary.is_none() {
            s.push_str("  PTM report (rank-1 target PSMs carrying each modification):\n");
            if self.ptm_counts.is_empty() {
                s.push_str("    (none — no modified target PSMs)\n");
            } else {
                for (name, count) in &self.ptm_counts {
                    s.push_str(&format!("    {name:<label_w$} : {count}\n"));
                }
            }
            s.push_str(&format!(
                "    {:<label_w$} : {}\n",
                "(unmodified)", self.unmodified_target_psms
            ));
        }
        s.push_str("───────────────────────────────────\n");
        if let Some(summary) = summary {
            s.push_str(&summary.render());
        }
        s
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct PsmFacts {
    pub run: usize,
    pub charge: u8,
    pub peptidoform: String,
    pub sequence: String,
    pub ppm: f64,
    pub ppm_calibrated: f64,
    pub spectrum_charge: Option<u8>,
    pub secondary: bool,
    pub missed_cleavages: u32,
    pub variable_mods: Vec<String>,
    pub fixed_carried: Vec<usize>,
    pub fixed_eligible: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct FixedModSpec {
    pub name: String,
    pub residue: Option<char>,
    pub location: ModLocation,
}

impl FixedModSpec {
    pub fn declared(params: &SearchParams) -> Vec<Self> {
        params
            .aa_set
            .distinct_mods()
            .into_iter()
            .filter(|m| m.fixed)
            .map(|m| FixedModSpec {
                name: m.name.clone(),
                residue: match m.residue {
                    ResidueSpec::Specific(r) => Some(r as char),
                    ResidueSpec::Wildcard => None,
                },
                location: m.location,
            })
            .collect()
    }

    fn applies_at(&self, residue: u8, locations: &[ModLocation]) -> bool {
        self.residue.is_none_or(|r| r as u8 == residue)
            && (self.location == ModLocation::Anywhere || locations.contains(&self.location))
    }

    fn place(&self) -> Option<&'static str> {
        match self.location {
            ModLocation::Anywhere => None,
            ModLocation::NTerm => Some("peptide N-terminus"),
            ModLocation::CTerm => Some("peptide C-terminus"),
            ModLocation::ProtNTerm => Some("protein N-terminus"),
            ModLocation::ProtCTerm => Some("protein C-terminus"),
        }
    }

    pub fn target(&self) -> String {
        match (self.residue, self.place()) {
            (Some(r), None) => r.to_string(),
            (Some(r), Some(p)) => format!("{r} at the {p}"),
            (None, Some(p)) => p.to_string(),
            (None, None) => "any residue".to_string(),
        }
    }

    fn site_text(&self) -> String {
        match (self.residue, self.place()) {
            (Some(r), None) => format!("a {r}"),
            (Some(r), Some(p)) => format!("a {r} at the {p}"),
            (None, Some(p)) => format!("the {p}"),
            (None, None) => "a residue".to_string(),
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct PeptideMods {
    pub variable: Vec<String>,
    pub fixed_carried: Vec<usize>,
    pub fixed_eligible: Vec<usize>,
}

fn site_locations(i: usize, n: usize, protein_n: bool, protein_c: bool) -> Vec<ModLocation> {
    let mut locations = Vec::new();
    if i == 0 {
        locations.push(if protein_n {
            ModLocation::ProtNTerm
        } else {
            ModLocation::NTerm
        });
    }
    if i + 1 == n {
        locations.push(if protein_c {
            ModLocation::ProtCTerm
        } else {
            ModLocation::CTerm
        });
    }
    if locations.is_empty() {
        locations.push(ModLocation::Anywhere);
    }
    locations
}

fn components<'a>(record: &'a str, declared: &[&str]) -> Vec<&'a str> {
    fn split<'a>(rest: &'a str, declared: &[&str]) -> Option<Vec<&'a str>> {
        if rest.is_empty() {
            return Some(Vec::new());
        }
        for d in declared {
            let Some(tail) = rest.strip_prefix(d) else {
                continue;
            };
            let tail = if tail.is_empty() {
                tail
            } else if let Some(t) = tail.strip_prefix('+') {
                t
            } else {
                continue;
            };
            if let Some(mut more) = split(tail, declared) {
                more.insert(0, &rest[..d.len()]);
                return Some(more);
            }
        }
        None
    }
    let mut longest_first: Vec<&str> = declared.to_vec();
    longest_first.sort_by_key(|d| std::cmp::Reverse(d.len()));
    split(record, &longest_first).unwrap_or_else(|| vec![record])
}

pub fn classify_mods(
    residues: &[(u8, Option<&str>)],
    protein_n: bool,
    protein_c: bool,
    fixed: &[FixedModSpec],
    declared: &[&str],
) -> PeptideMods {
    let n = residues.len();
    let mut variable = BTreeSet::new();
    let mut carried = BTreeSet::new();
    let mut eligible = BTreeSet::new();
    for (i, &(residue, record)) in residues.iter().enumerate() {
        let locations = site_locations(i, n, protein_n, protein_c);
        for (j, spec) in fixed.iter().enumerate() {
            if spec.applies_at(residue, &locations) {
                eligible.insert(j);
            }
        }
        let Some(record) = record else {
            continue;
        };
        for component in components(record, declared) {
            let matching: Vec<usize> = fixed
                .iter()
                .enumerate()
                .filter(|(_, spec)| spec.name == component && spec.applies_at(residue, &locations))
                .map(|(j, _)| j)
                .collect();
            if matching.is_empty() {
                variable.insert(component.to_string());
            } else {
                carried.extend(matching);
            }
        }
    }
    PeptideMods {
        variable: variable.into_iter().collect(),
        fixed_carried: carried.into_iter().collect(),
        fixed_eligible: eligible.into_iter().collect(),
    }
}

/// The facts of the `wanted` PSMs. `run_shift_ppm` is each run's precursor calibration shift;
/// a run past its end uses `params.precursor_mass_shift_ppm`.
pub fn psm_facts(
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &[Candidate],
    params: &SearchParams,
    run_ends: &[usize],
    run_shift_ppm: &[f64],
    wanted: &HashSet<&str>,
) -> HashMap<String, PsmFacts> {
    use crate::pin::{emitted_spec_ids, precursor_error};
    use crate::qpx::{bare_sequence, peptidoform_string};
    use model::mass::{ISOTOPE, PROTON};
    use search::precursor_cal::adjusted_observed_neutral_mass;

    let fixed = FixedModSpec::declared(params);
    let declared_mods = params.aa_set.distinct_mods();
    let declared: Vec<&str> = declared_mods.iter().map(|m| m.name.as_str()).collect();
    (0..queues.len())
        .into_par_iter()
        .filter(|&i| !queues[i].is_empty())
        .flat_map_iter(|i| {
            let spec = &spectra[i];
            let psms = queues[i].clone().into_rank_sorted_vec();
            let run = run_ends.partition_point(|&end| end <= i);
            let shift_ppm = run_shift_ppm
                .get(run)
                .copied()
                .unwrap_or(params.precursor_mass_shift_ppm);
            let mut out = Vec::new();
            for (id, psm) in emitted_spec_ids(spec, &psms, candidates) {
                if !wanted.contains(id.as_str()) {
                    continue;
                }
                let cand = &candidates[psm.primary_candidate_idx() as usize];
                let precursor_mz = psm.precursor_mz_override.unwrap_or(spec.precursor_mz);
                let (dm, theo_mz) = precursor_error(
                    precursor_mz,
                    f64::from(psm.charge_used),
                    i32::from(psm.isotope_offset),
                    cand.peptide.mass(),
                );
                let (_, _, enz_int) = crate::percolator_enz::enzyme_features(
                    cand,
                    params.enzyme,
                    &params.extra_enzymes,
                );
                let residues: Vec<(u8, Option<&str>)> = cand
                    .peptide
                    .residues
                    .iter()
                    .map(|aa| (aa.residue, aa.mod_.map(|m| m.name.as_str())))
                    .collect();
                let mods = classify_mods(
                    &residues,
                    cand.is_protein_n_term,
                    cand.is_protein_c_term,
                    &fixed,
                    &declared,
                );
                let ppm = dm / theo_mz * 1e6;
                let ppm_calibrated = if shift_ppm == 0.0 {
                    ppm
                } else {
                    let z = f64::from(psm.charge_used);
                    let observed =
                        adjusted_observed_neutral_mass(precursor_mz * z - z * PROTON, shift_ppm);
                    let dm =
                        (observed - ISOTOPE * f64::from(psm.isotope_offset) - cand.peptide.mass())
                            / z;
                    dm / theo_mz * 1e6
                };
                out.push((
                    id,
                    PsmFacts {
                        run,
                        charge: psm.charge_used,
                        peptidoform: peptidoform_string(cand),
                        sequence: bare_sequence(cand),
                        ppm,
                        ppm_calibrated,
                        spectrum_charge: spec
                            .precursor_charge
                            .and_then(|z| u8::try_from(z).ok())
                            .filter(|&z| z > 0),
                        secondary: psm.precursor_mz_override.is_some(),
                        missed_cleavages: enz_int.max(0) as u32,
                        variable_mods: mods.variable,
                        fixed_carried: mods.fixed_carried,
                        fixed_eligible: mods.fixed_eligible,
                    },
                ));
            }
            out
        })
        .collect()
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct SpectraTally {
    pub per_run: Vec<usize>,
    pub by_charge: BTreeMap<u8, usize>,
    pub unknown_charge: usize,
}

impl SpectraTally {
    pub fn of(spectra: &[Spectrum], run_ends: &[usize]) -> Self {
        let mut per_run = Vec::with_capacity(run_ends.len());
        let mut start = 0usize;
        for &end in run_ends {
            per_run.push(end.saturating_sub(start));
            start = start.max(end);
        }
        let mut by_charge: BTreeMap<u8, usize> = BTreeMap::new();
        let mut unknown_charge = 0usize;
        for s in spectra {
            match s.precursor_charge.and_then(|z| u8::try_from(z).ok()) {
                Some(z) if z > 0 => *by_charge.entry(z).or_insert(0) += 1,
                _ => unknown_charge += 1,
            }
        }
        Self {
            per_run,
            by_charge,
            unknown_charge,
        }
    }

    fn total(&self) -> usize {
        self.by_charge.values().sum::<usize>() + self.unknown_charge
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct RateRow {
    pub label: String,
    pub spectra: usize,
    pub psms: usize,
    pub rate: Option<f64>,
    pub z: Option<f64>,
}

const LOW_RATE_Z: f64 = -3.0;
const MIN_LOG_RATE_SCALE: f64 = 0.1;

#[derive(Debug, Clone, Copy)]
pub struct SummaryContext<'a> {
    pub rescorer: &'a str,
    pub cut: &'a str,
    pub decoys_at_cut: usize,
    pub protein_groups: Option<(f64, usize)>,
    pub fixed: &'a [FixedModSpec],
    pub calibration_shift_ppm: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct IdentificationSummary {
    pub rescorer: String,
    pub cut: String,
    pub protein_groups: Option<(f64, usize)>,
    pub spectra: usize,
    pub psms: usize,
    pub unlinked_psms: usize,
    pub precursors: usize,
    pub peptidoforms: usize,
    pub sequences: usize,
    pub decoys_at_cut: usize,
    pub mass_error_ppm: Option<(f64, f64, f64)>,
    pub calibration_shift_ppm: f64,
    pub mass_error_after_shift_ppm: Option<(f64, f64, f64)>,
    pub fixed_mods: Vec<FixedModRow>,
    pub missed_cleavages: Vec<(u32, usize)>,
    pub by_charge: Vec<RateRow>,
    pub charge_rates_comparable: bool,
    pub unknown_charge_spectra: usize,
    pub by_run: Vec<RateRow>,
    pub ptm_counts: Vec<(String, usize)>,
    pub unmodified_psms: usize,
    pub open_search: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct FixedModRow {
    pub name: String,
    pub target: String,
    pub site: String,
    pub psms: usize,
    pub eligible: usize,
}

fn run_rate_z(rates: &[Option<f64>]) -> Vec<Option<f64>> {
    let logs: Vec<f64> = rates.iter().flatten().map(|&r| r.max(1e-6).ln()).collect();
    if logs.len() < 3 {
        return vec![None; rates.len()];
    }
    let mut sorted = logs.clone();
    sorted.sort_by(f64::total_cmp);
    let median = quantile_sorted(&sorted, 0.5, false);
    let mut dev: Vec<f64> = logs.iter().map(|l| (l - median).abs()).collect();
    dev.sort_by(f64::total_cmp);
    let scale = (MAD_TO_SD * quantile_sorted(&dev, 0.5, false)).max(MIN_LOG_RATE_SCALE);
    rates
        .iter()
        .map(|r| r.map(|r| (r.max(1e-6).ln() - median) / scale))
        .collect()
}

fn rate_of(psms: usize, spectra: usize) -> Option<f64> {
    (spectra > 0 && psms <= spectra).then(|| psms as f64 / spectra as f64)
}

impl IdentificationSummary {
    pub fn compute(
        tally: &SpectraTally,
        run_names: &[String],
        facts: &[PsmFacts],
        unlinked_psms: usize,
        ctx: &SummaryContext,
    ) -> Self {
        let precursors: HashSet<(&str, u8)> = facts
            .iter()
            .map(|f| (f.peptidoform.as_str(), f.charge))
            .collect();
        let peptidoforms: HashSet<&str> = facts.iter().map(|f| f.peptidoform.as_str()).collect();
        let sequences: HashSet<&str> = facts.iter().map(|f| f.sequence.as_str()).collect();

        let mut ppm: Vec<f64> = facts
            .iter()
            .map(|f| f.ppm)
            .filter(|v| v.is_finite())
            .collect();
        ppm.sort_by(f64::total_cmp);
        let mass_error_ppm = (!ppm.is_empty()).then(|| {
            (
                quantile_sorted(&ppm, 0.5, false),
                quantile_sorted(&ppm, 0.25, false),
                quantile_sorted(&ppm, 0.75, false),
            )
        });
        let mass_error_after_shift_ppm = (ctx.calibration_shift_ppm != 0.0)
            .then(|| {
                let mut after: Vec<f64> = facts
                    .iter()
                    .map(|f| f.ppm_calibrated)
                    .filter(|v| v.is_finite())
                    .collect();
                after.sort_by(f64::total_cmp);
                (!after.is_empty()).then(|| {
                    (
                        quantile_sorted(&after, 0.5, false),
                        quantile_sorted(&after, 0.25, false),
                        quantile_sorted(&after, 0.75, false),
                    )
                })
            })
            .flatten();

        let mut missed: BTreeMap<u32, usize> = BTreeMap::new();
        let mut ptm: BTreeMap<&str, usize> = BTreeMap::new();
        let mut unmodified_psms = 0usize;
        let mut psm_by_charge: BTreeMap<u8, usize> = BTreeMap::new();
        let mut psm_by_run: Vec<usize> = vec![0; tally.per_run.len().max(run_names.len())];
        for f in facts {
            *missed.entry(f.missed_cleavages).or_insert(0) += 1;
            *psm_by_charge.entry(f.charge).or_insert(0) += 1;
            if let Some(n) = psm_by_run.get_mut(f.run) {
                *n += 1;
            }
            if f.variable_mods.is_empty() {
                unmodified_psms += 1;
            }
            for m in &f.variable_mods {
                *ptm.entry(m.as_str()).or_insert(0) += 1;
            }
        }
        let fixed_mods = ctx
            .fixed
            .iter()
            .enumerate()
            .map(|(j, spec)| FixedModRow {
                name: spec.name.clone(),
                target: spec.target(),
                site: spec.site_text(),
                psms: facts
                    .iter()
                    .filter(|f| f.fixed_carried.contains(&j))
                    .count(),
                eligible: facts
                    .iter()
                    .filter(|f| f.fixed_eligible.contains(&j))
                    .count(),
            })
            .collect();
        let mut ptm_counts: Vec<(String, usize)> =
            ptm.into_iter().map(|(n, c)| (n.to_string(), c)).collect();
        ptm_counts.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

        let charge_rates_comparable = unlinked_psms == 0
            && facts
                .iter()
                .all(|f| !f.secondary && f.spectrum_charge == Some(f.charge));
        let charges: BTreeSet<u8> = psm_by_charge
            .keys()
            .chain(tally.by_charge.keys())
            .copied()
            .collect();
        let by_charge = charges
            .into_iter()
            .map(|z| {
                let spectra = tally.by_charge.get(&z).copied().unwrap_or(0);
                let psms = psm_by_charge.get(&z).copied().unwrap_or(0);
                RateRow {
                    label: format!("{z}+"),
                    spectra,
                    psms,
                    rate: rate_of(psms, spectra).filter(|_| charge_rates_comparable),
                    z: None,
                }
            })
            .collect();

        let by_run = if tally.per_run.len() > 1 {
            let rates: Vec<Option<f64>> = tally
                .per_run
                .iter()
                .enumerate()
                .map(|(i, &s)| rate_of(psm_by_run[i], s))
                .collect();
            let zs = run_rate_z(&rates);
            tally
                .per_run
                .iter()
                .enumerate()
                .map(|(i, &spectra)| RateRow {
                    label: run_names
                        .get(i)
                        .cloned()
                        .unwrap_or_else(|| format!("run {}", i + 1)),
                    spectra,
                    psms: psm_by_run[i],
                    rate: rates[i],
                    z: zs[i],
                })
                .collect()
        } else {
            Vec::new()
        };

        Self {
            rescorer: ctx.rescorer.to_string(),
            cut: ctx.cut.to_string(),
            protein_groups: ctx.protein_groups,
            spectra: tally.total(),
            psms: facts.len() + unlinked_psms,
            unlinked_psms,
            precursors: precursors.len(),
            peptidoforms: peptidoforms.len(),
            sequences: sequences.len(),
            decoys_at_cut: ctx.decoys_at_cut,
            mass_error_ppm,
            calibration_shift_ppm: ctx.calibration_shift_ppm,
            mass_error_after_shift_ppm,
            fixed_mods,
            missed_cleavages: missed.into_iter().collect(),
            by_charge,
            charge_rates_comparable,
            unknown_charge_spectra: tally.unknown_charge,
            by_run,
            ptm_counts,
            unmodified_psms,
            open_search: None,
        }
    }

    pub fn render(&self) -> String {
        const LABEL_W: usize = 27;
        let line = |label: &str, value: String| format!("  {label:<LABEL_W$}: {value}\n");
        let linked = self.psms - self.unlinked_psms;

        let mut s = format!(
            "──────── andes identification summary ({}, {}) ────────\n",
            self.rescorer, self.cut
        );
        s.push_str(&line("Spectra searched", thousands(self.spectra)));
        s.push_str(&line(
            "Identified PSMs",
            format!(
                "{}  ({} of spectra)",
                thousands(self.psms),
                if self.psms <= self.spectra {
                    percent(self.psms, self.spectra)
                } else {
                    "-".to_string()
                }
            ),
        ));
        s.push_str(&line(
            "Precursors / peptidoforms",
            format!(
                "{} / {}   sequences: {}",
                thousands(self.precursors),
                thousands(self.peptidoforms),
                thousands(self.sequences)
            ),
        ));
        if let Some((q, n)) = self.protein_groups {
            s.push_str(&line(&format!("Protein groups (q ≤ {q})"), thousands(n)));
        }
        s.push_str(&line("Decoys at the cut", thousands(self.decoys_at_cut)));
        if let (Some((median, q1, q3)), None) = (self.mass_error_ppm, &self.open_search) {
            s.push_str(&line(
                "Precursor mass error (raw)",
                format!(
                    "median {} ppm (IQR {} … {})",
                    signed(median),
                    signed(q1),
                    signed(q3)
                ),
            ));
        }
        if let (Some((median, q1, q3)), None) = (self.mass_error_after_shift_ppm, &self.open_search) {
            s.push_str(&line(
                "  after calibration shift",
                format!(
                    "median {} ppm (IQR {} … {})  [shift {} ppm]",
                    signed(median),
                    signed(q1),
                    signed(q3),
                    signed(self.calibration_shift_ppm)
                ),
            ));
        }
        if !self.missed_cleavages.is_empty() {
            let parts: Vec<String> = self
                .missed_cleavages
                .iter()
                .map(|&(k, n)| format!("{k}: {}", percent(n, linked)))
                .collect();
            s.push_str(&line("Missed cleavages", parts.join("  ")));
        }
        if self.unlinked_psms > 0 {
            s.push_str(&format!(
                "  ({} identified PSMs could not be matched to a searched spectrum and are left \
                 out of the figures below the PSM count)\n",
                thousands(self.unlinked_psms)
            ));
        }

        s.push_str(&rate_table(
            "By charge",
            &self.by_charge,
            unknown_row(self.unknown_charge_spectra),
        ));
        if !self.charge_rates_comparable {
            s.push_str(
                "    per-charge rates are not shown: some PSMs sit at a charge other than their \
                 spectrum's annotated one (spectra with no annotated charge are searched at every \
                 charge in the range; chimeric secondaries use their own precursor)\n",
            );
        }
        if !self.by_run.is_empty() {
            s.push_str(&rate_table("By run", &self.by_run, None));
            if self
                .by_run
                .iter()
                .any(|r| r.z.is_some_and(|z| z <= LOW_RATE_Z))
            {
                s.push_str(
                    "    '<' marks a run whose rate is far below the other runs (robust z ≤ -3 \
                     on the log rate); it is a note only and changes no result\n",
                );
            }
        }

        for f in &self.fixed_mods {
            s.push_str(&format!(
                "  Fixed: {} ({}) on {} PSMs ({} of PSMs with {})\n",
                f.name,
                f.target,
                thousands(f.psms),
                percent(f.psms, f.eligible),
                f.site
            ));
        }
        if let Some(open) = &self.open_search {
            s.push_str(open);
        } else {
            s.push_str(
                "  Variable modifications on identified PSMs (unmodified: none of them; \
                 share of identified PSMs):\n",
            );
            let label_w = self
                .ptm_counts
                .iter()
                .map(|(n, _)| n.chars().count())
                .chain(std::iter::once("(unmodified)".len()))
                .max()
                .unwrap_or(12);
            let count_w = thousands(
                self.ptm_counts
                    .iter()
                    .map(|&(_, c)| c)
                    .chain(std::iter::once(self.unmodified_psms))
                    .max()
                    .unwrap_or(0),
            )
            .len();
            let ptm_line = |name: &str, count: usize| {
                format!(
                    "    {name:<label_w$} : {:>count_w$}  ({})\n",
                    thousands(count),
                    percent(count, linked)
                )
            };
            for (name, count) in &self.ptm_counts {
                s.push_str(&ptm_line(name, *count));
            }
            s.push_str(&ptm_line("(unmodified)", self.unmodified_psms));
        }
        s.push_str("───────────────────────────────────────────────────────────────\n");
        s
    }
}

fn unknown_row(unknown: usize) -> Option<(String, usize)> {
    (unknown > 0).then(|| ("unknown".to_string(), unknown))
}

fn rate_table(title: &str, rows: &[RateRow], extra: Option<(String, usize)>) -> String {
    let label_w = rows
        .iter()
        .map(|r| r.label.chars().count())
        .chain(extra.iter().map(|(l, _)| l.chars().count()))
        .max()
        .unwrap_or(0);
    let num_w = |f: &dyn Fn(&RateRow) -> usize, head: &str| {
        rows.iter()
            .map(|r| thousands(f(r)).len())
            .chain(std::iter::once(head.len()))
            .max()
            .unwrap_or(0)
    };
    let sp_w = num_w(&|r| r.spectra, "spectra")
        .max(extra.as_ref().map_or(0, |(_, n)| thousands(*n).len()));
    let ps_w = num_w(&|r| r.psms, "PSMs");
    let title_w = label_w.max(title.len());
    let mut s = format!(
        "  {title:<title_w$}  {:>sp_w$}  {:>ps_w$}  {:>6}\n",
        "spectra", "PSMs", "rate"
    );
    for r in rows {
        let rate = r
            .rate
            .map_or("-".to_string(), |v| format!("{:.1} %", v * 100.0));
        let flag = if r.z.is_some_and(|z| z <= LOW_RATE_Z) {
            "  <"
        } else {
            ""
        };
        s.push_str(&format!(
            "  {:<title_w$}  {:>sp_w$}  {:>ps_w$}  {:>6}{flag}\n",
            r.label,
            thousands(r.spectra),
            thousands(r.psms),
            rate
        ));
    }
    if let Some((label, n)) = extra {
        s.push_str(&format!(
            "  {label:<title_w$}  {:>sp_w$}  {:>ps_w$}  {:>6}\n",
            thousands(n),
            "-",
            "-"
        ));
    }
    s
}

/// `n` with a comma between each group of three digits (52,169).
pub fn thousands(n: usize) -> String {
    let digits = n.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(c);
    }
    out
}

fn percent(part: usize, whole: usize) -> String {
    if whole == 0 {
        "-".to_string()
    } else {
        format!("{:.1} %", 100.0 * part as f64 / whole as f64)
    }
}

fn signed(v: f64) -> String {
    let r = (v * 10.0).round() / 10.0;
    let r = if r == 0.0 { 0.0 } else { r };
    format!("{r:+.1}").replace('-', "−")
}

pub fn write_statistics_log(
    stats: &RunStatistics,
    summary: Option<&IdentificationSummary>,
    path: &Path,
) -> io::Result<()> {
    crate::atomic::write_atomic(path, |path| {
        let mut f = File::create(path)?;
        f.write_all(stats.render_with(summary).as_bytes())
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn statistics_log_is_written_through_a_partial_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("run.statistics.log");
        let stats = RunStatistics {
            precursor_tolerance: "10 ppm".to_string(),
            precursor_calibration: "Auto".to_string(),
            fragment_tolerance: "0.02 Da".to_string(),
            node_score_bin: "0.5 Da".to_string(),
            spectra_with_match: 3,
            target_psms: 2,
            decoy_psms: 1,
            unmodified_target_psms: 2,
            ptm_counts: Vec::new(),
        };
        write_statistics_log(&stats, None, &path).unwrap();
        assert_eq!(std::fs::read_to_string(&path).unwrap(), stats.render());
        assert!(!crate::atomic::partial(&path).exists());
    }

    #[test]
    fn render_lists_ptms_desc_and_unmodified() {
        let stats = RunStatistics {
            precursor_tolerance: "Symmetric(10.0 ppm)".to_string(),
            precursor_calibration: "Auto".to_string(),
            fragment_tolerance: "0.5 Da".to_string(),
            node_score_bin: "0.5 Da".to_string(),
            spectra_with_match: 100,
            target_psms: 60,
            decoy_psms: 40,
            unmodified_target_psms: 25,
            ptm_counts: vec![
                ("Carbamidomethyl".to_string(), 30),
                ("Oxidation".to_string(), 10),
            ],
        };
        let out = stats.render();
        assert!(out.contains("Final fragment tolerance  : 0.5 Da"));
        assert!(out.contains("Carbamidomethyl : 30"));
        assert!(out.contains("Oxidation       : 10"));
        assert!(out.contains("(unmodified)    : 25"));
        let cam = out.find("Carbamidomethyl").unwrap();
        let ox = out.find("Oxidation").unwrap();
        assert!(cam < ox, "PTMs should be sorted count-descending");
    }

    #[test]
    fn high_res_run_reports_matching_tolerance_apart_from_node_score_bin() {
        let params = SearchParams::default_tryptic(
            model::aa_set::AminoAcidSetBuilder::new_standard()
                .build()
                .unwrap(),
        );
        let stats = RunStatistics::compute(
            &[],
            &[],
            &params,
            &Tolerance::Ppm(20.0),
            &Tolerance::Ppm(9.87654),
            &Tolerance::Da(0.5),
        );
        let out = stats.render();
        assert!(out.contains(
            "Final fragment tolerance  : 20 ppm (retrieval and model features); \
             9.877 ppm (PIN match features, measured)"
        ));
        assert!(out.contains("Node-score bin width      : 0.5 Da"));
        assert!(!out.contains("Final fragment tolerance  : 0.5 Da"));

        let low = RunStatistics::compute(
            &[],
            &[],
            &params,
            &Tolerance::Da(0.5),
            &Tolerance::Da(0.5),
            &Tolerance::Da(0.5),
        );
        assert!(low
            .render()
            .contains("Final fragment tolerance  : 0.5 Da (retrieval and match features)"));
    }

    fn fact(run: usize, charge: u8, pep: &str, ppm: f64, missed: u32, mods: &[&str]) -> PsmFacts {
        PsmFacts {
            run,
            charge,
            peptidoform: pep.to_string(),
            sequence: pep.chars().filter(|c| c.is_ascii_uppercase()).collect(),
            ppm,
            ppm_calibrated: ppm,
            spectrum_charge: Some(charge),
            secondary: false,
            missed_cleavages: missed,
            variable_mods: mods.iter().map(|m| m.to_string()).collect(),
            fixed_carried: Vec::new(),
            fixed_eligible: Vec::new(),
        }
    }

    fn names(n: usize) -> Vec<String> {
        (1..=n).map(|i| format!("run{i}")).collect()
    }

    fn tally(per_run: Vec<usize>, charge2: usize, unknown: usize) -> SpectraTally {
        SpectraTally {
            per_run,
            by_charge: BTreeMap::from([(2, charge2)]),
            unknown_charge: unknown,
        }
    }

    fn ctx() -> SummaryContext<'static> {
        SummaryContext {
            rescorer: "Percolator",
            cut: "q ≤ 0.01",
            decoys_at_cut: 0,
            protein_groups: None,
            fixed: &[],
            calibration_shift_ppm: 0.0,
        }
    }

    fn summarise(
        tally: &SpectraTally,
        runs: usize,
        facts: &[PsmFacts],
        unlinked: usize,
    ) -> IdentificationSummary {
        IdentificationSummary::compute(tally, &names(runs), facts, unlinked, &ctx())
    }

    fn single_run_summary() -> IdentificationSummary {
        let tally = SpectraTally {
            per_run: vec![20],
            by_charge: BTreeMap::from([(2, 12), (3, 8)]),
            unknown_charge: 0,
        };
        let facts = vec![
            fact(0, 2, "AAK", 1.0, 0, &["Carbamidomethyl"]),
            fact(0, 2, "AAK", 3.0, 0, &["Carbamidomethyl"]),
            fact(0, 3, "AAK", -1.0, 1, &["Carbamidomethyl"]),
            fact(0, 2, "BBR", 2.0, 0, &[]),
        ];
        let ctx = SummaryContext {
            decoys_at_cut: 1,
            protein_groups: Some((0.01, 5)),
            ..ctx()
        };
        IdentificationSummary::compute(&tally, &names(1), &facts, 0, &ctx)
    }

    #[test]
    fn summary_counts_rates_and_distinct_levels() {
        let s = single_run_summary();
        assert_eq!(s.spectra, 20);
        assert_eq!(s.psms, 4);
        assert_eq!((s.precursors, s.peptidoforms, s.sequences), (3, 2, 2));
        assert_eq!(s.decoys_at_cut, 1);
        assert_eq!(s.missed_cleavages, vec![(0, 3), (1, 1)]);
        assert_eq!(s.ptm_counts, vec![("Carbamidomethyl".to_string(), 3)]);
        assert_eq!(s.unmodified_psms, 1);
        let z2 = &s.by_charge[0];
        assert_eq!((z2.label.as_str(), z2.spectra, z2.psms), ("2+", 12, 3));
        assert!((z2.rate.unwrap() - 0.25).abs() < 1e-12);
        let (m, q1, q3) = s.mass_error_ppm.unwrap();
        assert!((m - 1.5).abs() < 1e-12 && (q1 - 0.5).abs() < 1e-12 && (q3 - 2.25).abs() < 1e-12);
    }

    #[test]
    fn single_run_render_has_no_run_table() {
        let out = single_run_summary().render();
        assert!(out.contains("andes identification summary (Percolator, q ≤ 0.01)"));
        assert!(out.contains("  Spectra searched           : 20\n"));
        assert!(!out.contains("MS2"));
        assert!(out.contains("  Identified PSMs            : 4  (20.0 % of spectra)\n"));
        assert!(out.contains("  Precursors / peptidoforms  : 3 / 2   sequences: 2\n"));
        assert!(out.contains("  Protein groups (q ≤ 0.01)  : 5\n"));
        assert!(out.contains("  Decoys at the cut          : 1\n"));
        assert!(out.contains("Precursor mass error (raw)"));
        assert!(out.contains("median +1.5 ppm (IQR +0.5 … +2.3)"));
        assert!(!out.contains("after calibration shift"));
        assert!(!out.contains("Fixed:"));
        assert!(out.contains("  Missed cleavages           : 0: 75.0 %  1: 25.0 %\n"));
        assert!(out.contains("Carbamidomethyl : 3  (75.0 %)"));
        assert!(out.contains("(unmodified)    : 1  (25.0 %)"));
        assert!(!out.contains("By run"));
        assert!(!out.contains("per-charge rates are not shown"));
        assert!(out.contains("25.0 %"));
    }

    #[test]
    fn heading_names_the_rescoring_backend() {
        let native = SummaryContext {
            rescorer: "native rescorer",
            ..ctx()
        };
        let s = IdentificationSummary::compute(&tally(vec![5], 5, 0), &names(1), &[], 0, &native);
        let out = s.render();
        assert!(out.contains("andes identification summary (native rescorer, q ≤ 0.01)"));
        assert!(!out.contains("Percolator"));
    }

    #[test]
    fn unknown_charge_spectra_get_a_row_and_a_note_and_no_rates() {
        let mut psm3 = fact(0, 3, "AAK", 0.0, 0, &[]);
        psm3.spectrum_charge = None;
        let facts = vec![fact(0, 2, "AAK", 0.0, 0, &[]), psm3];
        let s = summarise(&tally(vec![10], 6, 4), 1, &facts, 0);
        assert_eq!(s.spectra, 10);
        assert!(!s.charge_rates_comparable);
        let out = s.render();
        assert!(out.contains("  unknown"));
        assert!(out.contains("per-charge rates are not shown"));
        assert!(out.contains("  3+"));
        assert!(!out.contains("16.7 %"), "per-charge rates are withheld");
        assert!(!out.contains("Protein groups"));
        assert!(out.contains("median +0.0 ppm"));
    }

    #[test]
    fn a_psm_at_another_charge_or_a_secondary_withholds_the_charge_rates() {
        let base = vec![
            fact(0, 2, "AAK", 0.0, 0, &[]),
            fact(0, 2, "BBK", 0.0, 0, &[]),
        ];
        let ok = summarise(&tally(vec![4], 4, 0), 1, &base, 0);
        assert!(ok.charge_rates_comparable);
        assert_eq!(ok.by_charge[0].rate, Some(0.5));

        let mut other = base.clone();
        other[1].charge = 3;
        let s = summarise(&tally(vec![4], 4, 0), 1, &other, 0);
        assert!(!s.charge_rates_comparable);
        assert!(s.by_charge.iter().all(|r| r.rate.is_none()));
        assert!(s.by_charge.iter().any(|r| r.psms > 0));

        let mut secondary = base.clone();
        secondary[1].secondary = true;
        let s = summarise(&tally(vec![4], 4, 0), 1, &secondary, 0);
        assert!(!s.charge_rates_comparable);
        assert!(s.render().contains("per-charge rates are not shown"));

        assert!(!summarise(&tally(vec![4], 4, 0), 1, &base, 1).charge_rates_comparable);
    }

    #[test]
    fn no_rate_exceeds_one_hundred_percent() {
        let facts = vec![
            fact(0, 2, "AAK", 0.0, 0, &[]),
            fact(0, 2, "BBK", 0.0, 0, &[]),
            fact(0, 2, "CCK", 0.0, 0, &[]),
        ];
        let one = summarise(&tally(vec![2], 2, 0), 1, &facts, 0);
        assert_eq!(one.by_charge[0].rate, None);
        let out = one.render();
        assert!(out.contains("Identified PSMs            : 3  (- of spectra)"));
        assert!(!out.contains("150.0 %"));
        let two = summarise(&tally(vec![2, 5], 7, 0), 2, &facts, 0);
        assert_eq!(two.by_run[0].rate, None);
        assert!(two.by_run[1].rate.is_none_or(|v| v <= 1.0));
        assert!(two
            .by_charge
            .iter()
            .all(|r| r.rate.is_none_or(|v| v <= 1.0)));
    }

    #[test]
    fn run_table_sums_and_flags_a_run_far_below_the_others() {
        let mut facts = Vec::new();
        for (run, n) in [(0usize, 50usize), (1, 48), (2, 52), (3, 5)] {
            for i in 0..n {
                facts.push(fact(run, 2, &format!("PEP{i}K"), 0.0, 0, &[]));
            }
        }
        let s = summarise(&tally(vec![100, 100, 100, 100], 400, 0), 4, &facts, 0);
        assert_eq!(s.by_run.iter().map(|r| r.psms).sum::<usize>(), s.psms);
        assert_eq!(s.by_run.iter().map(|r| r.spectra).sum::<usize>(), s.spectra);
        let flagged: Vec<&str> = s
            .by_run
            .iter()
            .filter(|r| r.z.is_some_and(|z| z <= LOW_RATE_Z))
            .map(|r| r.label.as_str())
            .collect();
        assert_eq!(flagged, vec!["run4"]);
        let out = s.render();
        assert!(out.contains("By run"));
        assert!(out.contains("run4") && out.contains("  <\n"));
        assert!(out.contains("changes no result"));
    }

    #[test]
    fn similar_runs_and_short_run_lists_are_not_flagged() {
        let facts: Vec<PsmFacts> = [(0usize, 50usize), (1, 49), (2, 47)]
            .iter()
            .flat_map(|&(r, n)| (0..n).map(move |i| fact(r, 2, &format!("P{i}K"), 0.0, 0, &[])))
            .collect();
        let s = summarise(&tally(vec![100, 100, 100], 300, 0), 3, &facts, 0);
        assert!(s.by_run.iter().all(|r| r.z.is_none_or(|z| z > LOW_RATE_Z)));
        assert!(!s.render().contains("changes no result"));
        let two = summarise(
            &tally(vec![100, 100], 200, 0),
            2,
            &[fact(0, 2, "AAK", 0.0, 0, &[])],
            0,
        );
        assert!(two.by_run.iter().all(|r| r.z.is_none()));
    }

    #[test]
    fn empty_accepted_set_renders_without_panicking() {
        let s = summarise(&tally(vec![5], 5, 0), 1, &[], 0);
        assert_eq!(s.psms, 0);
        assert!(s.mass_error_ppm.is_none());
        let out = s.render();
        assert!(out.contains("Identified PSMs            : 0  (0.0 % of spectra)"));
        assert!(!out.contains("Precursor mass error"));
        assert!(!out.contains("Missed cleavages"));
        assert!(out.contains("(unmodified)"));
    }

    #[test]
    fn ptms_sort_by_count_then_name_and_numbers_get_separators() {
        let facts = vec![
            fact(0, 2, "AK", 0.0, 0, &["Oxidation"]),
            fact(0, 2, "BK", 0.0, 0, &["Acetyl"]),
            fact(0, 2, "CK", 0.0, 0, &["Acetyl", "Oxidation"]),
            fact(0, 2, "DK", 0.0, 0, &["Zeta"]),
        ];
        let s = summarise(&tally(vec![2_000_000], 2_000_000, 0), 1, &facts, 0);
        let order: Vec<&str> = s.ptm_counts.iter().map(|(n, _)| n.as_str()).collect();
        assert_eq!(order, vec!["Acetyl", "Oxidation", "Zeta"]);
        assert!(s.render().contains("2,000,000"));
        assert_eq!(thousands(0), "0");
        assert_eq!(thousands(999), "999");
        assert_eq!(thousands(1000), "1,000");
        assert_eq!(thousands(1234567), "1,234,567");
    }

    #[test]
    fn unlinked_psms_are_counted_and_called_out() {
        let facts = vec![fact(0, 2, "AAK", 0.0, 0, &[])];
        let s = summarise(&tally(vec![10], 10, 0), 1, &facts, 2);
        assert_eq!(s.psms, 3);
        assert!(s
            .render()
            .contains("2 identified PSMs could not be matched"));
    }

    #[test]
    fn spectra_tally_splits_runs_and_charges() {
        let spec = |z: Option<i32>| Spectrum {
            title: String::new(),
            precursor_mz: 500.0,
            precursor_intensity: None,
            precursor_charge: z,
            rt_seconds: None,
            scan: None,
            peaks: Vec::new(),
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
            mobility_window: None,
        };
        let spectra = vec![spec(Some(2)), spec(None), spec(Some(3)), spec(Some(2))];
        let t = SpectraTally::of(&spectra, &[2, 4]);
        assert_eq!(t.per_run, vec![2, 2]);
        assert_eq!(t.by_charge, BTreeMap::from([(2, 2), (3, 1)]));
        assert_eq!(t.unknown_charge, 1);
        assert_eq!(t.total(), 4);
    }

    #[test]
    fn rescored_report_replaces_the_pre_fdr_ptm_list_and_appends_the_summary() {
        let stats = RunStatistics {
            precursor_tolerance: "10 ppm".to_string(),
            precursor_calibration: "Auto".to_string(),
            fragment_tolerance: "0.02 Da".to_string(),
            node_score_bin: "0.5 Da".to_string(),
            spectra_with_match: 3,
            target_psms: 2,
            decoy_psms: 1,
            unmodified_target_psms: 2,
            ptm_counts: Vec::new(),
        };
        let plain = stats.render();
        assert_eq!(plain, stats.render_with(None));
        assert!(plain.contains("PTM report (rank-1 target PSMs"));
        let summary = single_run_summary();
        let out = stats.render_with(Some(&summary));
        assert!(!out.contains("PTM report"));
        assert!(out.contains("Rank-1 PSMs (pre-FDR)     : 2 target, 1 decoy"));
        assert!(out.ends_with(&summary.render()));
    }

    fn spec(name: &str, residue: Option<char>, location: ModLocation) -> FixedModSpec {
        FixedModSpec {
            name: name.to_string(),
            residue,
            location,
        }
    }

    fn peptide<'a>(seq: &str, mods: &[(usize, &'a str)]) -> Vec<(u8, Option<&'a str>)> {
        seq.bytes()
            .enumerate()
            .map(|(i, r)| (r, mods.iter().find(|(p, _)| *p == i).map(|(_, m)| *m)))
            .collect()
    }

    #[test]
    fn components_keep_declared_names_with_a_plus_whole() {
        let declared = ["Oxidation", "Acetyl", "Label+Heavy", "Label"];
        assert_eq!(
            components("Oxidation+Acetyl", &declared),
            vec!["Oxidation", "Acetyl"]
        );
        assert_eq!(components("Label+Heavy", &declared), vec!["Label+Heavy"]);
        assert_eq!(
            components("Label+Heavy+Oxidation", &declared),
            vec!["Label+Heavy", "Oxidation"]
        );
        assert_eq!(
            components("Found+Elsewhere", &declared),
            vec!["Found+Elsewhere"]
        );
    }

    #[test]
    fn classify_splits_fixed_from_variable_by_declaration() {
        let fixed = vec![spec("Carbamidomethyl", Some('C'), ModLocation::Anywhere)];
        let declared = ["Carbamidomethyl", "Oxidation", "Acetyl"];
        let residues = peptide("ACMK", &[(1, "Carbamidomethyl+Acetyl"), (2, "Oxidation")]);
        let m = classify_mods(&residues, false, false, &fixed, &declared);
        assert_eq!(m.variable, vec!["Acetyl", "Oxidation"]);
        assert_eq!(m.fixed_carried, vec![0]);
        assert_eq!(m.fixed_eligible, vec![0]);
        let none = classify_mods(&peptide("AMK", &[]), false, false, &fixed, &declared);
        assert!(none.fixed_eligible.is_empty() && none.fixed_carried.is_empty());
    }

    #[test]
    fn a_name_with_a_plus_is_one_modification_and_one_declaration() {
        let fixed = vec![spec("Label+Heavy", Some('K'), ModLocation::Anywhere)];
        let declared = ["Label+Heavy", "Label", "Heavy"];
        let m = classify_mods(
            &peptide("AKR", &[(1, "Label+Heavy")]),
            false,
            false,
            &fixed,
            &declared,
        );
        assert_eq!(m.fixed_carried, vec![0]);
        assert!(m.variable.is_empty());
    }

    #[test]
    fn terminal_fixed_declarations_are_separate_rows_with_their_own_denominator() {
        let fixed = vec![
            spec("TMT", Some('K'), ModLocation::Anywhere),
            spec("TMT", None, ModLocation::NTerm),
        ];
        let declared = ["TMT"];
        let both = classify_mods(
            &peptide("AKR", &[(0, "TMT"), (1, "TMT")]),
            false,
            false,
            &fixed,
            &declared,
        );
        assert_eq!(both.fixed_carried, vec![0, 1]);
        assert_eq!(both.fixed_eligible, vec![0, 1]);
        let k_first = classify_mods(
            &peptide("KA", &[(0, "TMT+TMT")]),
            false,
            false,
            &fixed,
            &declared,
        );
        assert_eq!(k_first.fixed_carried, vec![0, 1]);
        let prot = classify_mods(
            &peptide("KA", &[(0, "TMT")]),
            true,
            false,
            &fixed,
            &declared,
        );
        assert_eq!(prot.fixed_carried, vec![0]);
        assert_eq!(prot.fixed_eligible, vec![0]);
    }

    #[test]
    fn fixed_rows_count_declarations_against_their_eligible_psms_and_stay_within_100() {
        let fixed = vec![
            spec("Carbamidomethyl", Some('C'), ModLocation::Anywhere),
            spec("TMT", None, ModLocation::NTerm),
        ];
        let mut a = fact(0, 2, "ACK", 0.0, 0, &["Oxidation"]);
        a.fixed_carried = vec![0, 1];
        a.fixed_eligible = vec![0, 1];
        let mut b = fact(0, 2, "AAK", 0.0, 0, &[]);
        b.fixed_carried = vec![1];
        b.fixed_eligible = vec![1];
        let mut c = fact(0, 2, "CCK", 0.0, 0, &[]);
        c.fixed_carried = vec![0, 1];
        c.fixed_eligible = vec![0, 1];
        let ctx = SummaryContext {
            fixed: &fixed,
            ..ctx()
        };
        let s =
            IdentificationSummary::compute(&tally(vec![10], 10, 0), &names(1), &[a, b, c], 0, &ctx);
        assert_eq!(
            s.fixed_mods,
            vec![
                FixedModRow {
                    name: "Carbamidomethyl".to_string(),
                    target: "C".to_string(),
                    site: "a C".to_string(),
                    psms: 2,
                    eligible: 2,
                },
                FixedModRow {
                    name: "TMT".to_string(),
                    target: "peptide N-terminus".to_string(),
                    site: "the peptide N-terminus".to_string(),
                    psms: 3,
                    eligible: 3,
                },
            ]
        );
        assert!(s.fixed_mods.iter().all(|f| f.psms <= f.eligible));
        let out = s.render();
        assert!(out.contains("  Fixed: Carbamidomethyl (C) on 2 PSMs (100.0 % of PSMs with a C)\n"));
        assert!(out.contains(
            "  Fixed: TMT (peptide N-terminus) on 3 PSMs (100.0 % of PSMs with the peptide N-terminus)\n"
        ));
        assert_eq!(s.ptm_counts, vec![("Oxidation".to_string(), 1)]);
        assert_eq!(s.unmodified_psms, 2);
    }

    #[test]
    fn calibration_shift_adds_the_error_after_the_shift() {
        let mut a = fact(0, 2, "AAK", 2.0, 0, &[]);
        a.ppm_calibrated = -1.0;
        let mut b = fact(0, 2, "BBK", 4.0, 0, &[]);
        b.ppm_calibrated = 1.0;
        let with_shift = SummaryContext {
            calibration_shift_ppm: 3.0,
            ..ctx()
        };
        let s = IdentificationSummary::compute(
            &tally(vec![10], 10, 0),
            &names(1),
            &[a.clone(), b.clone()],
            0,
            &with_shift,
        );
        assert_eq!(s.mass_error_ppm.unwrap().0, 3.0);
        assert_eq!(s.mass_error_after_shift_ppm.unwrap().0, 0.0);
        let out = s.render();
        assert!(out.contains("Precursor mass error (raw) : median +3.0 ppm"));
        assert!(out.contains("after calibration shift"));
        assert!(out.contains("median +0.0 ppm (IQR −0.5 … +0.5)  [shift +3.0 ppm]"));
        let none = summarise(&tally(vec![10], 10, 0), 1, &[a, b], 0);
        assert!(none.mass_error_after_shift_ppm.is_none());
        assert!(!none.render().contains("after calibration shift"));
    }

    #[test]
    fn signed_ppm_uses_a_true_minus_and_no_negative_zero() {
        assert_eq!(signed(-1.14), "−1.1");
        assert_eq!(signed(0.4), "+0.4");
        assert_eq!(signed(-0.04), "+0.0");
    }
}
