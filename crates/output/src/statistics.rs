use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use model::modification::{ModLocation, ResidueSpec};
use model::spectrum::Spectrum;
use model::stats::{quantile_sorted, MAD_TO_SD};
use model::tolerance::Tolerance;
use rayon::prelude::*;
use search::candidates::Candidate;
use search::psm::TopNQueue;
use search::params::SearchParams;

pub use crate::spectra::SpectraTally;

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

        let best: Vec<Option<Option<(bool, BTreeSet<&str>)>>> = queues
            .par_iter()
            .map(|q| {
                if q.is_empty() {
                    return None;
                }
                let psms = q.clone().into_rank_sorted_vec();
                let Some(best) = psms.first() else {
                    return Some(None);
                };
                let cand = &candidates[best.primary_candidate_idx() as usize];
                let mut names: BTreeSet<&str> = BTreeSet::new();
                if !cand.is_decoy {
                    for aa in &cand.peptide.residues {
                        if let Some(m) = aa.mod_.as_ref() {
                            names.insert(m.name.as_str());
                        }
                    }
                }
                Some(Some((cand.is_decoy, names)))
            })
            .collect();
        for q in best {
            let Some(q) = q else {
                continue;
            };
            spectra_with_match += 1;
            let Some((is_decoy, names)) = q else {
                continue;
            };
            if is_decoy {
                decoy_psms += 1;
                continue;
            }
            target_psms += 1;
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

/// One emitted PIN row reduced to what [`psm_facts_of_rows`] reads: the parts of its SpecId,
/// and the matched candidate, charge, isotope offset and precursor m/z. 24 bytes a row, where
/// the queue's `PsmMatch` with its features is over 300.
#[derive(Debug, Clone, PartialEq)]
pub struct FactRow {
    rank: u32,
    row_idx: u32,
    candidate: u32,
    charge: u8,
    isotope: i8,
    multi_row: bool,
    has_mz_override: bool,
    mz_override: f64,
}

/// The emitted PIN rows of each spectrum (indexed like the queues), in PIN order.
pub type FactRows = Vec<Box<[FactRow]>>;

/// The emitted PIN rows of every spectrum as [`FactRow`]s.
pub fn fact_rows(
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &[Candidate],
) -> FactRows {
    use crate::rows::emitted_rows;
    queues
        .par_iter()
        .enumerate()
        .map(|(i, q)| {
            if q.is_empty() {
                return Box::default();
            }
            let psms = q.clone().into_rank_sorted_vec();
            let ranked = emitted_rows(&spectra[i], &psms, candidates);
            let multi_row = ranked.len() > 1;
            ranked
                .into_iter()
                .enumerate()
                .map(|(row_idx, (rank, psm))| FactRow {
                    rank,
                    row_idx: row_idx as u32,
                    candidate: psm.primary_candidate_idx(),
                    charge: psm.charge_used,
                    isotope: psm.isotope_offset,
                    multi_row,
                    has_mz_override: psm.precursor_mz_override.is_some(),
                    mz_override: psm.precursor_mz_override.unwrap_or(0.0),
                })
                .collect()
        })
        .collect()
}

pub fn psm_facts(
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &[Candidate],
    params: &SearchParams,
    run_ends: &[usize],
    file_cal: &[search::FileCalibration],
    wanted: &HashSet<&str>,
) -> HashMap<String, PsmFacts> {
    let rows = fact_rows(spectra, queues, candidates);
    psm_facts_of_rows(spectra, &rows, candidates, params, run_ends, file_cal, wanted)
}

/// The [`PsmFacts`] of the `wanted` SpecIds among `rows`. `file_cal` holds each run's precursor
/// calibration; a run past its end (or no per-file calibration) uses `params`' shift.
pub fn psm_facts_of_rows(
    spectra: &[Spectrum],
    rows: &[Box<[FactRow]>],
    candidates: &[Candidate],
    params: &SearchParams,
    run_ends: &[usize],
    file_cal: &[search::FileCalibration],
    wanted: &HashSet<&str>,
) -> HashMap<String, PsmFacts> {
    use crate::pin::{emitted_spec_id, precursor_error};
    use crate::qpx::{bare_sequence, peptidoform_string};
    use model::mass::{ISOTOPE, PROTON};
    use search::precursor::calibration::adjusted_observed_neutral_mass;

    let fixed = FixedModSpec::declared(params);
    let declared_mods = params.aa_set.distinct_mods();
    let declared: Vec<&str> = declared_mods.iter().map(|m| m.name.as_str()).collect();
    (0..rows.len())
        .into_par_iter()
        .flat_map_iter(|i| rows[i].iter().map(move |row| (i, row)))
        .filter_map(|(i, row)| {
            let spec = &spectra[i];
            let id = emitted_spec_id(spec, row.rank, row.row_idx as usize, row.multi_row);
            if !wanted.contains(id.as_str()) {
                return None;
            }
            let run = run_ends.partition_point(|&end| end <= i);
            let shift_ppm = file_cal
                .get(run)
                .map_or(params.precursor_mass_shift_ppm, |f| f.shift_ppm);
            let cand = &candidates[row.candidate as usize];
            let precursor_mz = if row.has_mz_override {
                row.mz_override
            } else {
                spec.precursor_mz
            };
            let (dm, theo_mz) = precursor_error(
                precursor_mz,
                f64::from(row.charge),
                i32::from(row.isotope),
                cand.peptide.mass(),
            );
            let (_, _, enz_int) =
                crate::percolator::enzyme::enzyme_features(cand, params.enzyme, &params.extra_enzymes);
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
                let z = f64::from(row.charge);
                let observed =
                    adjusted_observed_neutral_mass(precursor_mz * z - z * PROTON, shift_ppm);
                let dm = (observed - ISOTOPE * f64::from(row.isotope) - cand.peptide.mass()) / z;
                dm / theo_mz * 1e6
            };
            Some((
                id,
                PsmFacts {
                    run,
                    charge: row.charge,
                    peptidoform: peptidoform_string(cand),
                    sequence: bare_sequence(cand),
                    ppm,
                    ppm_calibrated,
                    spectrum_charge: spec
                        .precursor_charge
                        .and_then(|z| u8::try_from(z).ok())
                        .filter(|&z| z > 0),
                    secondary: row.has_mz_override,
                    missed_cleavages: enz_int.max(0) as u32,
                    variable_mods: mods.variable,
                    fixed_carried: mods.fixed_carried,
                    fixed_eligible: mods.fixed_eligible,
                },
            ))
        })
        .collect()
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
#[path = "stats/tests/stats.rs"]
mod tests;
