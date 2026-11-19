use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use std::path::Path;
use std::sync::Arc;

use arrow::array::{ArrayRef, Float64Array, StringArray, UInt32Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;

use model::stats::{median_mad, median_mut, Even, MAD_TO_SD};
use model::unimod::{Specificity, UnimodDb};

use crate::percolator::PercolatorPsm;
use crate::versioned_parquet::{read_table, write_table, ReadTable, TableVersion, BATCH_ROWS};

const ISOTOPE_DA: f64 = model::mass::ISOTOPE;
const REFINED_DA: f64 = 0.0005;
const MIN_PEAK_PSMS: usize = 2;
const DEFAULT_PRECISION_DA: f64 = 0.002;
const MIN_PRECISION_DA: f64 = 0.0005;
const ZERO_WINDOW_DA: f64 = 0.02;
const MIN_CALIBRATION_PSMS: usize = 20;
const SPLIT_SIGMAS: f64 = 4.0;
const GRID_PER_SIGMA: f64 = 4.0;
const KERNEL_SIGMAS: f64 = 4.0;
const VALLEY_RATIO: f64 = 0.5;
const CANDIDATE_SIGMAS: f64 = 6.0;
const PRIOR_SINGLE: f64 = 1.0;
const PRIOR_RARE: f64 = 0.1;
const PRIOR_ISOTOPE: f64 = 0.2;
const PRIOR_SINGLE_ISOTOPE: f64 = 0.1;
const PRIOR_PAIR: f64 = 0.05;
const PRIOR_UNEXPLAINED: f64 = 0.05;
const RARE_SITE_WEIGHT: f64 = 0.1;
const PAIR_SOURCE_POSTERIOR: f64 = 0.5;
const PAIR_SOURCE_MIN_PSMS: usize = 20;
const PAIR_SOURCE_MAX_DECOY_RATIO: f64 = 0.01;
const TOP_EXPLANATIONS: usize = 3;
const EM_ITERATIONS: usize = 200;
const EM_PSEUDOCOUNT: f64 = 1.0;
const EM_TOLERANCE: f64 = 1e-6;
const UNEXPLAINED: &str = "unexplained";
const SITE_NEUTRAL: f64 = 0.5;
const MIXTURE_WIDTH: f64 = 2.0;

#[derive(Debug, Clone, PartialEq)]
pub struct PinOffset {
    pub psm_id: String,
    pub decoy: bool,
    pub delta: f64,
    pub calc_mass: f64,
    pub refined: bool,
    pub site: usize,
    pub sequence: String,
    pub run: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Explanation {
    pub label: String,
    pub mass: f64,
    pub posterior: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct OffsetPeak {
    pub id: usize,
    pub mass: f64,
    pub width: f64,
    pub psms: usize,
    pub peptides: usize,
    pub decoys: usize,
    pub residues: Vec<(String, usize)>,
    pub explanations: Vec<Explanation>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PsmAssignment {
    pub psm_id: String,
    pub delta: f64,
    pub peak: Option<usize>,
    pub residue: String,
    pub explanation: String,
    pub posterior: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct OffsetReport {
    pub ppm_shift: f64,
    pub run_ppm_shifts: Vec<(String, f64)>,
    pub precision: f64,
    pub refined_precision: f64,
    pub peaks: Vec<OffsetPeak>,
    pub assignments: Vec<PsmAssignment>,
}

pub fn read_pin_offsets(path: &Path) -> std::io::Result<Vec<PinOffset>> {
    let reader = std::io::BufReader::new(std::fs::File::open(path)?);
    parse_pin_offsets(reader)
}

pub fn parse_pin_offsets<R: BufRead>(reader: R) -> std::io::Result<Vec<PinOffset>> {
    let bad = |m: String| std::io::Error::new(std::io::ErrorKind::InvalidData, m);
    let mut lines = reader.lines();
    let header = match lines.next() {
        Some(h) => h?,
        None => return Ok(Vec::new()),
    };
    let cols: Vec<&str> = header.split('\t').collect();
    let col = |name: &str| {
        cols.iter()
            .position(|c| *c == name)
            .ok_or_else(|| bad(format!("PIN has no {name} column")))
    };
    let (id_c, label_c, exp_c, calc_c, iso_c, delta_c, site_c, pep_c) = (
        col("SpecId")?,
        col("Label")?,
        col("ExpMass")?,
        col("CalcMass")?,
        col("isotope_error")?,
        col("DeltaMass")?,
        col("OpenSite")?,
        col("Peptide")?,
    );
    let run_c = cols.iter().position(|c| *c == "FileName");
    let mut out = Vec::new();
    for line in lines {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let get = |i: usize| {
            f.get(i)
                .copied()
                .ok_or_else(|| bad(format!("PIN row has fewer than {} fields", i + 1)))
        };
        let label = get(label_c)?;
        if label != "1" && label != "-1" {
            continue;
        }
        let num = |i: usize, name: &str| -> std::io::Result<f64> {
            get(i)?.parse().map_err(|e| bad(format!("{name}: {e}")))
        };
        let site: usize = get(site_c)?
            .parse()
            .map_err(|e| bad(format!("OpenSite: {e}")))?;
        let delta = num(delta_c, "DeltaMass")?;
        let calc_mass = num(calc_c, "CalcMass")?;
        let precursor_delta =
            num(exp_c, "ExpMass")? - num(iso_c, "isotope_error")? * ISOTOPE_DA - calc_mass;
        out.push(PinOffset {
            psm_id: get(id_c)?.to_string(),
            decoy: label == "-1",
            delta,
            calc_mass,
            refined: (delta - precursor_delta).abs() > REFINED_DA,
            site,
            sequence: stripped_sequence(get(pep_c)?),
            run: match run_c {
                Some(i) => get(i)?.to_string(),
                None => String::new(),
            },
        });
    }
    Ok(out)
}

pub fn stripped_sequence(peptide: &str) -> String {
    let core = match (peptide.find('.'), peptide.rfind('.')) {
        (Some(a), Some(b)) if b > a + 1 && a <= 1 && b + 2 >= peptide.len() => &peptide[a + 1..b],
        _ => peptide,
    };
    let mut out = String::with_capacity(core.len());
    let mut depth = 0usize;
    for c in core.chars() {
        match c {
            '[' | '(' => depth += 1,
            ']' | ')' => depth = depth.saturating_sub(1),
            c if depth == 0 && c.is_ascii_uppercase() => out.push(c),
            _ => {}
        }
    }
    out
}

fn site_label(sequence: &str, site: usize) -> Option<String> {
    let c = sequence.chars().nth(site.checked_sub(1)?)?;
    Some(c.to_string())
}

fn median(values: &mut [f64]) -> f64 {
    median_mut(values, Even::Upper).unwrap_or(f64::NAN)
}

fn mad(values: &[f64]) -> f64 {
    median_mad(values, Even::Upper).map_or(0.0, |(_, mad)| MAD_TO_SD * mad)
}

fn shift_of(rows: &[&PinOffset]) -> Option<f64> {
    let mut ppm: Vec<f64> = rows
        .iter()
        .filter(|r| !r.refined && r.delta.abs() < ZERO_WINDOW_DA && r.calc_mass > 0.0)
        .map(|r| r.delta / r.calc_mass * 1e6)
        .collect();
    (ppm.len() >= MIN_CALIBRATION_PSMS).then(|| median(&mut ppm))
}

pub struct Calibration {
    pub pooled: f64,
    pub runs: Vec<(String, f64)>,
}

impl Calibration {
    pub fn shift(&self, run: &str) -> f64 {
        self.runs
            .iter()
            .find(|(r, _)| r == run)
            .map_or(self.pooled, |(_, s)| *s)
    }
}

pub fn calibration(rows: &[&PinOffset]) -> Calibration {
    let pooled = shift_of(rows).unwrap_or(0.0);
    let mut names: Vec<&str> = rows.iter().map(|r| r.run.as_str()).collect();
    names.sort_unstable();
    names.dedup();
    let runs = names
        .into_iter()
        .map(|name| {
            let own: Vec<&PinOffset> = rows.iter().copied().filter(|r| r.run == name).collect();
            (name.to_string(), shift_of(&own).unwrap_or(pooled))
        })
        .collect();
    Calibration { pooled, runs }
}

pub fn corrected_delta(row: &PinOffset, ppm: f64) -> f64 {
    if row.refined {
        row.delta
    } else {
        row.delta - ppm * 1e-6 * (row.calc_mass + row.delta)
    }
}

pub fn measured_precision(deltas: &[f64]) -> f64 {
    let near: Vec<f64> = deltas
        .iter()
        .copied()
        .filter(|d| d.abs() < ZERO_WINDOW_DA)
        .collect();
    if near.len() < MIN_CALIBRATION_PSMS {
        return DEFAULT_PRECISION_DA;
    }
    mad(&near).max(MIN_PRECISION_DA)
}

fn gap_groups(sorted: &[f64], precision: f64) -> Vec<std::ops::Range<usize>> {
    let mut out = Vec::new();
    let mut start = 0;
    for i in 1..=sorted.len() {
        if i == sorted.len() || sorted[i] - sorted[i - 1] > SPLIT_SIGMAS * precision {
            out.push(start..i);
            start = i;
        }
    }
    out
}

fn valley_splits(values: &[f64], precision: f64) -> Vec<f64> {
    let (lo, hi) = (values[0], values[values.len() - 1]);
    let step = precision / GRID_PER_SIGMA;
    let n = ((hi - lo) / step).ceil() as usize + 1;
    if n < 3 {
        return Vec::new();
    }
    let mut counts = vec![0.0f64; n];
    for v in values {
        counts[(((v - lo) / step).round() as usize).min(n - 1)] += 1.0;
    }
    let radius = (KERNEL_SIGMAS * GRID_PER_SIGMA).ceil() as isize;
    let kernel: Vec<f64> = (-radius..=radius)
        .map(|k| (-0.5 * (k as f64 / GRID_PER_SIGMA).powi(2)).exp())
        .collect();
    let density: Vec<f64> = (0..n as isize)
        .map(|i| {
            (-radius..=radius)
                .filter_map(|k| {
                    let j = i + k;
                    (0..n as isize)
                        .contains(&j)
                        .then(|| counts[j as usize] * kernel[(k + radius) as usize])
                })
                .sum()
        })
        .collect();
    let mut maxima: Vec<usize> = (0..n)
        .filter(|&i| {
            (i == 0 || density[i] > density[i - 1]) && (i + 1 == n || density[i] >= density[i + 1])
        })
        .collect();
    maxima.dedup();
    let mut splits = Vec::new();
    let mut left = match maxima.first() {
        Some(&m) => m,
        None => return splits,
    };
    for &right in &maxima[1..] {
        let (valley, &vmin) = density[left..=right]
            .iter()
            .enumerate()
            .min_by(|a, b| a.1.total_cmp(b.1))
            .map(|(i, v)| (left + i, v))
            .unwrap_or((left, &density[left]));
        if vmin < VALLEY_RATIO * density[left].min(density[right]) {
            splits.push(lo + valley as f64 * step);
            left = right;
        } else if density[right] > density[left] {
            left = right;
        }
    }
    splits
}

pub fn offset_peaks(sorted: &[f64], precision: f64) -> Vec<std::ops::Range<usize>> {
    let mut out = Vec::new();
    for g in gap_groups(sorted, precision) {
        let values = &sorted[g.clone()];
        let mut start = g.start;
        for cut in valley_splits(values, precision) {
            let end = g.start + values.partition_point(|v| *v < cut);
            if end > start {
                out.push(start..end);
                start = end;
            }
        }
        out.push(start..g.end);
    }
    out
}

struct Candidate {
    label: String,
    mass: f64,
    prior: f64,
    specificities: Vec<Specificity>,
}

fn gauss_log(x: f64, mu: f64, s: f64) -> f64 {
    let z = (x - mu) / s;
    -0.5 * z * z - s.ln() - 0.5 * (2.0 * std::f64::consts::PI).ln()
}

fn specificity_matches(spec: &Specificity, (residue, n_term, c_term): &(String, bool, bool)) -> bool {
    let site_ok = spec.site == *residue
        || (spec.site == "N-term" && *n_term)
        || (spec.site == "C-term" && *c_term);
    let position_ok = (!spec.position.needs_n_terminus() || *n_term)
        && (!spec.position.needs_c_terminus() || *c_term);
    site_ok && position_ok
}

fn site_support(c: &Candidate, localised: &[(String, bool, bool)]) -> f64 {
    if localised.is_empty() {
        return 1.0;
    }
    if c.specificities.is_empty() {
        return SITE_NEUTRAL;
    }
    let hits: f64 = localised
        .iter()
        .map(|l| {
            c.specificities
                .iter()
                .filter(|spec| specificity_matches(spec, l))
                .map(|spec| if spec.hidden { RARE_SITE_WEIGHT } else { 1.0 })
                .fold(0.0, f64::max)
        })
        .sum();
    (hits + 0.5) / (localised.len() as f64 + 1.0)
}

fn single_candidates(db: &UnimodDb) -> Vec<Candidate> {
    let mut out: Vec<Candidate> = db
        .mods()
        .filter(|m| m.mono_mass.is_finite() && m.mono_mass != 0.0)
        .map(|m| Candidate {
            label: format!("{} ({})", m.name, m.accession),
            mass: m.mono_mass,
            prior: if m.specificities.iter().all(|spec| spec.hidden) || m.isotopic_label {
                PRIOR_RARE
            } else {
                PRIOR_SINGLE
            },
            specificities: m.specificities.clone(),
        })
        .collect();
    out.sort_by(|a, b| a.mass.total_cmp(&b.mass).then_with(|| a.label.cmp(&b.label)));
    out
}

struct Scored {
    label: String,
    mass: f64,
    prior: f64,
    log_lik: f64,
}

fn scored_candidates(
    delta: f64,
    s: f64,
    span: f64,
    localised: &[(String, bool, bool)],
    singles: &[Candidate],
    pairs: &[Candidate],
) -> Vec<Scored> {
    let reach = CANDIDATE_SIGMAS * s;
    let neutral = if localised.is_empty() { 1.0 } else { SITE_NEUTRAL };
    let mut out: Vec<Scored> = Vec::new();
    let mut push = |label: String, cmass: f64, prior: f64, support: f64| {
        if (delta - cmass).abs() <= reach {
            out.push(Scored {
                label,
                mass: cmass,
                prior,
                log_lik: gauss_log(delta, cmass, s) + support.ln(),
            });
        }
    };
    push("unmodified".to_string(), 0.0, PRIOR_SINGLE, neutral);
    for k in [-2i32, -1, 1, 2] {
        push(
            format!("isotope error {k:+}"),
            f64::from(k) * ISOTOPE_DA,
            PRIOR_ISOTOPE,
            neutral,
        );
    }
    let lo = singles.partition_point(|c| c.mass < delta - reach - 2.0 * ISOTOPE_DA);
    for c in &singles[lo..] {
        if c.mass > delta + reach + 2.0 * ISOTOPE_DA {
            break;
        }
        let support = site_support(c, localised);
        push(c.label.clone(), c.mass, c.prior, support);
        for k in [-1i32, 1] {
            push(
                format!("{} {k:+} isotope", c.label),
                c.mass + f64::from(k) * ISOTOPE_DA,
                c.prior * PRIOR_SINGLE_ISOTOPE,
                support,
            );
        }
    }
    for c in pairs {
        push(c.label.clone(), c.mass, c.prior, site_support(c, localised));
    }
    out.push(Scored {
        label: UNEXPLAINED.to_string(),
        mass: f64::NAN,
        prior: PRIOR_UNEXPLAINED,
        log_lik: neutral.ln() - span.max(1.0).ln(),
    });
    out
}

fn single_species(scored: Vec<Scored>, peak_mass: f64) -> Vec<Explanation> {
    let scores: Vec<f64> = scored.iter().map(|c| c.prior.ln() + c.log_lik).collect();
    let max = scores.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let total: f64 = scores.iter().map(|v| (v - max).exp()).sum();
    let mut out: Vec<Explanation> = scored
        .into_iter()
        .zip(scores)
        .map(|(c, v)| Explanation {
            label: c.label,
            mass: if c.mass.is_nan() { peak_mass } else { c.mass },
            posterior: (v - max).exp() / total,
        })
        .collect();
    out.sort_by(|a, b| {
        b.posterior
            .total_cmp(&a.posterior)
            .then_with(|| a.label.cmp(&b.label))
    });
    out.truncate(TOP_EXPLANATIONS);
    out
}

struct Mixture {
    peak: Vec<Explanation>,
    psms: Vec<Explanation>,
}

fn fit_mixture(rows: Vec<Vec<Scored>>, peak_mass: f64) -> Mixture {
    let mut index: HashMap<String, usize> = HashMap::new();
    let mut comps: Vec<(String, f64, f64)> = Vec::new();
    let lik: Vec<Vec<(usize, f64)>> = rows
        .into_iter()
        .map(|r| {
            r.into_iter()
                .map(|c| {
                    let k = *index.entry(c.label.clone()).or_insert_with(|| {
                        let mass = if c.mass.is_nan() { peak_mass } else { c.mass };
                        comps.push((c.label, mass, c.prior));
                        comps.len() - 1
                    });
                    (k, c.log_lik)
                })
                .collect()
        })
        .collect();
    let total: f64 = comps.iter().map(|c| c.2).sum();
    let prior: Vec<f64> = comps.iter().map(|c| c.2 / total).collect();
    let n = lik.len() as f64;
    let e_step = |pi: &[f64], acc: &mut [f64], tops: &mut Vec<(usize, f64)>| {
        tops.clear();
        for r in &lik {
            let max = r
                .iter()
                .map(|&(k, ll)| pi[k].ln() + ll)
                .fold(f64::NEG_INFINITY, f64::max);
            let sum: f64 = r.iter().map(|&(k, ll)| (pi[k].ln() + ll - max).exp()).sum();
            let mut best = (0usize, -1.0f64);
            for &(k, ll) in r {
                let w = (pi[k].ln() + ll - max).exp() / sum;
                acc[k] += w;
                if w > best.1 || (w == best.1 && comps[k].0 < comps[best.0].0) {
                    best = (k, w);
                }
            }
            tops.push(best);
        }
    };
    let mut pi = prior.clone();
    let mut tops: Vec<(usize, f64)> = Vec::with_capacity(lik.len());
    for _ in 0..EM_ITERATIONS {
        let mut acc = vec![0.0f64; comps.len()];
        e_step(&pi, &mut acc, &mut tops);
        let next: Vec<f64> = acc
            .iter()
            .zip(&prior)
            .map(|(a, p)| (a + EM_PSEUDOCOUNT * p) / (n + EM_PSEUDOCOUNT))
            .collect();
        let change = next
            .iter()
            .zip(&pi)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        pi = next;
        if change < EM_TOLERANCE {
            break;
        }
    }
    let mut acc = vec![0.0f64; comps.len()];
    e_step(&pi, &mut acc, &mut tops);
    let mut peak: Vec<Explanation> = comps
        .iter()
        .zip(&pi)
        .map(|((label, mass, _), p)| Explanation {
            label: label.clone(),
            mass: *mass,
            posterior: *p,
        })
        .collect();
    peak.sort_by(|a, b| {
        b.posterior
            .total_cmp(&a.posterior)
            .then_with(|| a.label.cmp(&b.label))
    });
    peak.truncate(TOP_EXPLANATIONS);
    let psms = tops
        .iter()
        .map(|&(k, w)| Explanation {
            label: comps[k].0.clone(),
            mass: comps[k].1,
            posterior: w,
        })
        .collect();
    Mixture { peak, psms }
}

pub fn decoy_windows(bounds: &[(f64, f64)], precision: f64) -> Vec<(f64, f64)> {
    let mut order: Vec<usize> = (0..bounds.len()).collect();
    order.sort_by(|&a, &b| bounds[a].0.total_cmp(&bounds[b].0));
    let mut windows = vec![(0.0, 0.0); bounds.len()];
    for (j, &k) in order.iter().enumerate() {
        let (lo, hi) = bounds[k];
        let mut wlo = lo - precision;
        let mut whi = hi + precision;
        if j > 0 {
            wlo = wlo.max((bounds[order[j - 1]].1 + lo) / 2.0);
        }
        if let Some(&next) = order.get(j + 1) {
            whi = whi.min((hi + bounds[next].0) / 2.0);
        }
        windows[k] = (wlo, whi);
    }
    windows
}

struct RawPeak {
    members: Vec<usize>,
    mass: f64,
    width: f64,
    localised: Vec<(String, bool, bool)>,
}

pub fn build_report(
    pin: &[PinOffset],
    accepted: &HashSet<String>,
    target_scores: &HashMap<String, f64>,
    decoy_scores: &HashMap<String, f64>,
    db: &UnimodDb,
) -> OffsetReport {
    let targets: Vec<&PinOffset> = pin
        .iter()
        .filter(|r| !r.decoy && accepted.contains(&r.psm_id))
        .collect();
    let cal = calibration(&targets);
    let deltas: Vec<f64> = targets
        .iter()
        .map(|r| corrected_delta(r, cal.shift(&r.run)))
        .collect();
    let precision = measured_precision(&deltas);
    let span = pin
        .iter()
        .map(|r| r.delta)
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), d| (lo.min(d), hi.max(d)));
    let span = if span.0.is_finite() { span.1 - span.0 } else { 1.0 };
    let threshold = targets
        .iter()
        .filter_map(|r| target_scores.get(&r.psm_id))
        .copied()
        .fold(f64::INFINITY, f64::min);
    let mut decoy_deltas: Vec<f64> = pin
        .iter()
        .filter(|r| r.decoy)
        .filter(|r| {
            decoy_scores
                .get(&r.psm_id)
                .is_some_and(|s| threshold.is_finite() && *s >= threshold)
        })
        .map(|r| corrected_delta(r, cal.shift(&r.run)))
        .collect();
    decoy_deltas.sort_by(|a, b| a.total_cmp(b));

    let mut order: Vec<usize> = (0..targets.len()).collect();
    order.sort_by(|&a, &b| deltas[a].total_cmp(&deltas[b]));
    let sorted: Vec<f64> = order.iter().map(|&i| deltas[i]).collect();
    let mut raw: Vec<RawPeak> = offset_peaks(&sorted, precision)
        .into_iter()
        .filter(|r| r.len() >= MIN_PEAK_PSMS)
        .map(|r| {
            let members: Vec<usize> = order[r.clone()].to_vec();
            let values = &sorted[r];
            let width = mad(values).max(precision);
            let mass = median(&mut values.to_vec());
            let localised: Vec<(String, bool, bool)> = members
                .iter()
                .filter_map(|&i| {
                    let t = targets[i];
                    let residue = site_label(&t.sequence, t.site)?;
                    Some((residue, t.site == 1, t.site == t.sequence.len()))
                })
                .collect();
            RawPeak {
                members,
                mass,
                width,
                localised,
            }
        })
        .collect();
    raw.sort_by(|a, b| {
        b.members
            .len()
            .cmp(&a.members.len())
            .then_with(|| a.mass.total_cmp(&b.mass))
    });

    let singles = single_candidates(db);
    let bounds: Vec<(f64, f64)> = raw
        .iter()
        .map(|p| {
            p.members.iter().map(|&i| deltas[i]).fold(
                (f64::INFINITY, f64::NEG_INFINITY),
                |(lo, hi), d| (lo.min(d), hi.max(d)),
            )
        })
        .collect();
    let windows = decoy_windows(&bounds, precision);
    let decoys_in = |k: usize| {
        let (dlo, dhi) = windows[k];
        decoy_deltas.partition_point(|d| *d < dhi) - decoy_deltas.partition_point(|d| *d < dlo)
    };
    let mut refined_widths: Vec<f64> = raw
        .iter()
        .filter(|p| {
            p.members.len() >= PAIR_SOURCE_MIN_PSMS
                && 2 * p.members.iter().filter(|&&i| targets[i].refined).count() >= p.members.len()
        })
        .map(|p| p.width)
        .collect();
    let refined_precision = if refined_widths.is_empty() {
        precision
    } else {
        median(&mut refined_widths).max(precision)
    };
    let explain_peak = |p: &RawPeak, pairs: &[Candidate]| -> (Vec<Explanation>, Vec<Explanation>) {
        let refined = p.members.iter().filter(|&&i| targets[i].refined).count();
        let expected = if 2 * refined >= p.members.len() {
            refined_precision
        } else {
            precision
        };
        if p.width <= MIXTURE_WIDTH * expected {
            let n = p.members.len() as f64;
            let s = (p.width * p.width / n + expected * expected).sqrt();
            let peak = single_species(
                scored_candidates(p.mass, s, span, &p.localised, &singles, pairs),
                p.mass,
            );
            let top = peak.first().cloned().unwrap_or(Explanation {
                label: UNEXPLAINED.to_string(),
                mass: p.mass,
                posterior: 1.0,
            });
            return (peak, vec![top; p.members.len()]);
        }
        let rows: Vec<Vec<Scored>> = p
            .members
            .iter()
            .map(|&i| {
                let t = targets[i];
                let s = if t.refined { refined_precision } else { precision };
                let localised: Vec<(String, bool, bool)> = site_label(&t.sequence, t.site)
                    .map(|r| vec![(r, t.site == 1, t.site == t.sequence.len())])
                    .unwrap_or_default();
                scored_candidates(deltas[i], s, span, &localised, &singles, pairs)
            })
            .collect();
        let m = fit_mixture(rows, p.mass);
        (m.peak, m.psms)
    };
    let mut sources: Vec<&Candidate> = Vec::new();
    for (k, p) in raw.iter().enumerate() {
        if p.members.len() < PAIR_SOURCE_MIN_PSMS
            || decoys_in(k) as f64 > PAIR_SOURCE_MAX_DECOY_RATIO * p.members.len() as f64
        {
            continue;
        }
        let Some(top) = explain_peak(p, &[]).0.into_iter().next() else {
            continue;
        };
        if top.posterior < PAIR_SOURCE_POSTERIOR {
            continue;
        }
        if let Some(c) = singles.iter().find(|c| c.label == top.label) {
            if !sources.iter().any(|s| s.label == c.label) {
                sources.push(c);
            }
        }
    }
    let mut pairs: Vec<Candidate> = Vec::new();
    for (i, a) in sources.iter().enumerate() {
        for b in &sources[i..] {
            pairs.push(Candidate {
                label: format!("{} + {}", a.label, b.label),
                mass: a.mass + b.mass,
                prior: PRIOR_PAIR * a.prior * b.prior,
                specificities: a
                    .specificities
                    .iter()
                    .chain(b.specificities.iter())
                    .cloned()
                    .collect(),
            });
        }
    }

    let mut peaks = Vec::with_capacity(raw.len());
    let mut assignments = Vec::with_capacity(targets.len());
    let mut assigned = vec![false; targets.len()];
    for (k, p) in raw.iter().enumerate() {
        let id = k + 1;
        let (explanations, tops) = explain_peak(p, &pairs);
        let peptides: HashSet<&str> = p
            .members
            .iter()
            .map(|&i| targets[i].sequence.as_str())
            .collect();
        let mut residue_counts: HashMap<&str, usize> = HashMap::new();
        for (r, _, _) in &p.localised {
            *residue_counts.entry(r.as_str()).or_default() += 1;
        }
        let mut residues: Vec<(String, usize)> = residue_counts
            .into_iter()
            .map(|(r, n)| (r.to_string(), n))
            .collect();
        residues.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        for (&i, top) in p.members.iter().zip(tops.iter()) {
            assigned[i] = true;
            let t = targets[i];
            assignments.push(PsmAssignment {
                psm_id: t.psm_id.clone(),
                delta: deltas[i],
                peak: Some(id),
                residue: site_label(&t.sequence, t.site).unwrap_or_default(),
                explanation: top.label.clone(),
                posterior: top.posterior,
            });
        }
        peaks.push(OffsetPeak {
            id,
            mass: p.mass,
            width: p.width,
            psms: p.members.len(),
            peptides: peptides.len(),
            decoys: decoys_in(k),
            residues,
            explanations,
        });
    }
    for (i, t) in targets.iter().enumerate() {
        if !assigned[i] {
            assignments.push(PsmAssignment {
                psm_id: t.psm_id.clone(),
                delta: deltas[i],
                peak: None,
                residue: site_label(&t.sequence, t.site).unwrap_or_default(),
                explanation: String::new(),
                posterior: 0.0,
            });
        }
    }
    assignments.sort_by(|a, b| a.psm_id.cmp(&b.psm_id));
    OffsetReport {
        ppm_shift: cal.pooled,
        run_ppm_shifts: cal.runs,
        precision,
        refined_precision,
        peaks,
        assignments,
    }
}

pub fn report_from_results(
    pin_path: &Path,
    accepted: &[&PercolatorPsm],
    decoys: &HashMap<String, PercolatorPsm>,
    db: &UnimodDb,
) -> std::io::Result<OffsetReport> {
    let pin = read_pin_offsets(pin_path)?;
    let ids: HashSet<String> = accepted.iter().map(|p| p.psm_id.clone()).collect();
    let target_scores: HashMap<String, f64> = accepted
        .iter()
        .filter_map(|p| p.score.map(|s| (p.psm_id.clone(), s)))
        .collect();
    let decoy_scores: HashMap<String, f64> = decoys
        .iter()
        .filter_map(|(id, p)| p.score.map(|s| (id.clone(), s)))
        .collect();
    Ok(build_report(&pin, &ids, &target_scores, &decoy_scores, db))
}

fn residues_field(residues: &[(String, usize)]) -> String {
    residues
        .iter()
        .take(3)
        .map(|(r, n)| format!("{r}:{n}"))
        .collect::<Vec<_>>()
        .join(";")
}

const PEAKS_TABLE: TableVersion = TableVersion {
    key: "andes.open_offsets.schema",
    version: "1.0",
    major: 1,
    what: "open-search mass-offset table",
};
const PSMS_TABLE: TableVersion = TableVersion {
    key: "andes.open_psms.schema",
    version: "1.0",
    major: 1,
    what: "open-search PSM offset table",
};
pub const PPM_SHIFT_KEY: &str = "andes.open_offsets.ppm_shift";
pub const RUN_PPM_SHIFTS_KEY: &str = "andes.open_offsets.run_ppm_shifts";
pub const PRECISION_KEY: &str = "andes.open_offsets.precision_da";
pub const REFINED_PRECISION_KEY: &str = "andes.open_offsets.refined_precision_da";

fn report_meta(report: &OffsetReport) -> Vec<(&'static str, String)> {
    vec![
        (PPM_SHIFT_KEY, report.ppm_shift.to_string()),
        (
            RUN_PPM_SHIFTS_KEY,
            report
                .run_ppm_shifts
                .iter()
                .map(|(run, shift)| format!("{run}\t{shift}"))
                .collect::<Vec<_>>()
                .join("\n"),
        ),
        (PRECISION_KEY, report.precision.to_string()),
        (REFINED_PRECISION_KEY, report.refined_precision.to_string()),
    ]
}

fn peaks_schema() -> Schema {
    let mut fields = vec![
        Field::new("peak", DataType::UInt32, false),
        Field::new("mass_shift", DataType::Float64, false),
        Field::new("width_da", DataType::Float64, false),
        Field::new("psms", DataType::UInt32, false),
        Field::new("peptides", DataType::UInt32, false),
        Field::new("decoys", DataType::UInt32, false),
        Field::new("decoy_ratio", DataType::Float64, false),
        Field::new("residues", DataType::Utf8, false),
    ];
    for k in 1..=TOP_EXPLANATIONS {
        fields.push(Field::new(format!("explanation_{k}"), DataType::Utf8, true));
        fields.push(Field::new(format!("mass_{k}"), DataType::Float64, true));
        fields.push(Field::new(format!("posterior_{k}"), DataType::Float64, true));
    }
    Schema::new(fields)
}

fn psms_schema() -> Schema {
    Schema::new(vec![
        Field::new("psm_id", DataType::Utf8, false),
        Field::new("delta_mass", DataType::Float64, false),
        Field::new("peak", DataType::UInt32, true),
        Field::new("residue", DataType::Utf8, true),
        Field::new("explanation", DataType::Utf8, true),
        Field::new("posterior", DataType::Float64, true),
    ])
}

fn peaks_batch(schema: &Arc<Schema>, peaks: &[OffsetPeak]) -> Result<RecordBatch, String> {
    let u32s = |f: &dyn Fn(&OffsetPeak) -> usize| -> ArrayRef {
        Arc::new(UInt32Array::from_iter_values(peaks.iter().map(|p| f(p) as u32)))
    };
    let f64s = |f: &dyn Fn(&OffsetPeak) -> f64| -> ArrayRef {
        Arc::new(Float64Array::from_iter_values(peaks.iter().map(f)))
    };
    let mut cols: Vec<ArrayRef> = vec![
        u32s(&|p| p.id),
        f64s(&|p| p.mass),
        f64s(&|p| p.width),
        u32s(&|p| p.psms),
        u32s(&|p| p.peptides),
        u32s(&|p| p.decoys),
        f64s(&|p| p.decoys as f64 / p.psms as f64),
        Arc::new(StringArray::from_iter_values(
            peaks.iter().map(|p| residues_field(&p.residues)),
        )),
    ];
    for k in 0..TOP_EXPLANATIONS {
        cols.push(Arc::new(StringArray::from_iter(
            peaks
                .iter()
                .map(|p| p.explanations.get(k).map(|x| x.label.clone())),
        )));
        cols.push(Arc::new(Float64Array::from_iter(
            peaks.iter().map(|p| p.explanations.get(k).map(|x| x.mass)),
        )));
        cols.push(Arc::new(Float64Array::from_iter(
            peaks
                .iter()
                .map(|p| p.explanations.get(k).map(|x| x.posterior)),
        )));
    }
    RecordBatch::try_new(schema.clone(), cols).map_err(|e| e.to_string())
}

fn psms_batch(schema: &Arc<Schema>, rows: &[PsmAssignment]) -> Result<RecordBatch, String> {
    let text = |v: &str| (!v.is_empty()).then(|| v.to_string());
    let cols: Vec<ArrayRef> = vec![
        Arc::new(StringArray::from_iter_values(rows.iter().map(|a| a.psm_id.as_str()))),
        Arc::new(Float64Array::from_iter_values(rows.iter().map(|a| a.delta))),
        Arc::new(UInt32Array::from_iter(rows.iter().map(|a| a.peak.map(|p| p as u32)))),
        Arc::new(StringArray::from_iter(rows.iter().map(|a| text(&a.residue)))),
        Arc::new(StringArray::from_iter(rows.iter().map(|a| text(&a.explanation)))),
        Arc::new(Float64Array::from_iter(
            rows.iter()
                .map(|a| a.peak.map(|_| a.posterior)),
        )),
    ];
    RecordBatch::try_new(schema.clone(), cols).map_err(|e| e.to_string())
}

pub fn write_peaks_table(path: &Path, report: &OffsetReport) -> Result<(), String> {
    let schema = Arc::new(peaks_schema());
    let batches = report
        .peaks
        .chunks(BATCH_ROWS)
        .map(|c| peaks_batch(&schema, c));
    write_table(
        path,
        &PEAKS_TABLE,
        PEAKS_TABLE.version,
        schema.clone(),
        batches,
        &report_meta(report),
    )
}

pub fn write_psms_table(path: &Path, report: &OffsetReport) -> Result<(), String> {
    let schema = Arc::new(psms_schema());
    let batches = report
        .assignments
        .chunks(BATCH_ROWS)
        .map(|c| psms_batch(&schema, c));
    write_table(
        path,
        &PSMS_TABLE,
        PSMS_TABLE.version,
        schema.clone(),
        batches,
        &report_meta(report),
    )
}

pub fn read_peaks_table(path: &Path) -> Result<ReadTable, String> {
    read_table(path, &PEAKS_TABLE)
}

pub fn read_psms_table(path: &Path) -> Result<ReadTable, String> {
    read_table(path, &PSMS_TABLE)
}

const SUMMARY_PEAKS: usize = 10;
const SUMMARY_MODS: usize = 10;
const SUBSTANTIAL_PEAK_PSMS: usize = 20;
const ISOTOPE_PICK: &str = "isotope pick, no modification";

fn thousands(n: usize) -> String {
    let s = n.to_string();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(c);
    }
    out
}

fn share(n: usize, of: usize) -> String {
    if of == 0 {
        "-".to_string()
    } else {
        format!("{:.1} %", 100.0 * n as f64 / of as f64)
    }
}

pub fn modification_of(label: &str) -> String {
    if label.starts_with("isotope error") {
        return ISOTOPE_PICK.to_string();
    }
    for suffix in [" +1 isotope", " -1 isotope"] {
        if let Some(base) = label.strip_suffix(suffix) {
            return base.to_string();
        }
    }
    label.to_string()
}

pub fn summary_block(report: &OffsetReport) -> String {
    let total = report.assignments.len();
    let mut unmodified = 0usize;
    let mut isotope = 0usize;
    let mut unexplained = 0usize;
    let mut outside = 0usize;
    let mut by_mod: HashMap<String, usize> = HashMap::new();
    for a in &report.assignments {
        match a.peak {
            None => outside += 1,
            Some(_) if a.explanation == "unmodified" => unmodified += 1,
            Some(_) if a.explanation == UNEXPLAINED => unexplained += 1,
            Some(_) => {
                let m = modification_of(&a.explanation);
                if m == ISOTOPE_PICK {
                    isotope += 1;
                } else {
                    *by_mod.entry(m).or_default() += 1;
                }
            }
        }
    }
    let explained: usize = by_mod.values().sum();
    let mut mods: Vec<(String, usize)> = by_mod.into_iter().collect();
    mods.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    let substantial = report
        .peaks
        .iter()
        .filter(|p| p.psms >= SUBSTANTIAL_PEAK_PSMS)
        .count();
    let line = |label: &str, value: String| format!("    {label:<29}: {value}\n");
    let mut s = String::from("  Open search (mass offsets of the identified PSMs):\n");
    if report.run_ppm_shifts.len() > 1 {
        let shifts: Vec<f64> = report.run_ppm_shifts.iter().map(|(_, v)| *v).collect();
        let lo = shifts.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = shifts.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        s.push_str(&line(
            "Precursor calibration",
            format!(
                "per run, {} runs, {:+.2} … {:+.2} ppm (measured on each run's zero-offset PSMs)",
                report.run_ppm_shifts.len(),
                lo,
                hi
            ),
        ));
    } else {
        s.push_str(&line(
            "Precursor calibration",
            format!(
                "{:+.2} ppm (measured on zero-offset PSMs, removed from precursor offsets)",
                report.ppm_shift
            ),
        ));
    }
    s.push_str(&line(
        "Offset precision",
        format!(
            "{:.4} Da precursor, {:.4} Da fragment-refined",
            report.precision, report.refined_precision
        ),
    ));
    s.push_str(&line("PSMs", thousands(total)));
    for (label, n) in [
        ("  unmodified (zero offset)", unmodified),
        ("  explained offset", explained),
        ("  isotope pick, no mod.", isotope),
        ("  unexplained offset", unexplained),
        ("  outside every peak", outside),
    ] {
        s.push_str(&line(label, format!("{}  ({})", thousands(n), share(n, total))));
    }
    s.push_str(&line(
        "Offset peaks",
        format!(
            "{}  ({} with at least {SUBSTANTIAL_PEAK_PSMS} PSMs)",
            thousands(report.peaks.len()),
            thousands(substantial)
        ),
    ));
    if !mods.is_empty() {
        s.push_str("    By modification (PSMs; isotope picks counted with their modification):\n");
        let w = mods
            .iter()
            .take(SUMMARY_MODS)
            .map(|(m, _)| m.chars().count())
            .max()
            .unwrap_or(0);
        for (m, n) in mods.iter().take(SUMMARY_MODS) {
            s.push_str(&format!(
                "      {m:<w$}  {:>8}  ({})\n",
                thousands(*n),
                share(*n, total)
            ));
        }
        if mods.len() > SUMMARY_MODS {
            let rest: usize = mods.iter().skip(SUMMARY_MODS).map(|(_, n)| n).sum();
            s.push_str(&format!(
                "      {:<w$}  {:>8}  ({})\n",
                format!("{} others", mods.len() - SUMMARY_MODS),
                thousands(rest),
                share(rest, total)
            ));
        }
    }
    if !report.peaks.is_empty() {
        s.push_str(&format!(
            "    Largest peaks:\n      {:>11}  {:>7}  {:>6}  {:<22}  top explanation (posterior)\n",
            "offset (Da)", "PSMs", "decoys", "residues"
        ));
        for p in report.peaks.iter().take(SUMMARY_PEAKS) {
            let top = p
                .explanations
                .first()
                .map_or_else(String::new, |t| format!("{} ({:.2})", t.label, t.posterior));
            s.push_str(&format!(
                "      {:>+11.4}  {:>7}  {:>6}  {:<22}  {}\n",
                p.mass,
                thousands(p.psms),
                thousands(p.decoys),
                residues_field(&p.residues),
                top
            ));
        }
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(id: &str, decoy: bool, delta: f64, site: usize, seq: &str) -> PinOffset {
        PinOffset {
            psm_id: id.to_string(),
            decoy,
            delta,
            calc_mass: 1500.0,
            refined: false,
            site,
            sequence: seq.to_string(),
            run: "r1".to_string(),
        }
    }

    fn report(rows: &[PinOffset]) -> OffsetReport {
        let accepted: HashSet<String> = rows
            .iter()
            .filter(|r| !r.decoy)
            .map(|r| r.psm_id.clone())
            .collect();
        let target_scores: HashMap<String, f64> = rows
            .iter()
            .filter(|r| !r.decoy)
            .map(|r| (r.psm_id.clone(), 1.0))
            .collect();
        let decoy_scores: HashMap<String, f64> = rows
            .iter()
            .filter(|r| r.decoy)
            .map(|r| (r.psm_id.clone(), 2.0))
            .collect();
        build_report(
            rows,
            &accepted,
            &target_scores,
            &decoy_scores,
            model::unimod::bundled(),
        )
    }

    fn zero_rows(n: usize, ppm: f64) -> Vec<PinOffset> {
        (0..n)
            .map(|i| {
                let jitter = 0.0001 * (i % 5) as f64 - 0.0002;
                row(&format!("z{i}"), false, ppm * 1e-6 * 1500.0 + jitter, 0, "PEPTIDEK")
            })
            .collect()
    }

    #[test]
    fn stripped_sequence_drops_flanks_and_mods() {
        assert_eq!(stripped_sequence("K.PEM[+15.995]TIDEK.R"), "PEMTIDEK");
        assert_eq!(stripped_sequence("-.ACDK.-"), "ACDK");
        assert_eq!(stripped_sequence("PEPTIDE"), "PEPTIDE");
    }

    #[test]
    fn unimod_specificities_keep_hidden_flags_positions_and_label_classes() {
        let db = model::unimod::bundled();
        let ox = db.by_accession("UNIMOD:35").expect("Oxidation");
        assert!(ox.specificities.iter().any(|s| s.site == "M" && !s.hidden), "{:?}", ox.specificities);
        assert!(ox.specificities.iter().any(|s| s.hidden), "{:?}", ox.specificities);
        assert!(!ox.isotopic_label);
        let pyro = db.by_accession("UNIMOD:28").expect("Gln->pyro-Glu");
        assert!(
            pyro.specificities
                .iter()
                .any(|s| s.site == "Q" && s.position.needs_n_terminus()),
            "{:?}",
            pyro.specificities
        );
        assert!(db.by_accession("UNIMOD:994").expect("Label:15N(1)").isotopic_label);
    }

    #[test]
    fn an_n_terminal_specificity_supports_only_an_n_terminal_residue() {
        let pyro = single_candidates(model::unimod::bundled())
            .into_iter()
            .find(|c| c.label.starts_with("Gln->pyro-Glu (UNIMOD:28)"))
            .expect("pyro-Glu candidate");
        let first_q = site_support(&pyro, &[("Q".to_string(), true, false)]);
        let middle_q = site_support(&pyro, &[("Q".to_string(), false, false)]);
        assert!(first_q > middle_q, "{first_q} {middle_q}");
        assert!((middle_q - 0.25).abs() < 1e-12, "{middle_q}");
    }

    #[test]
    fn each_run_is_calibrated_on_its_own_zero_offset_psms() {
        let true_delta = 79.96633;
        let mut rows: Vec<PinOffset> = Vec::new();
        for (run, ppm) in [("a", 5.0), ("b", -5.0)] {
            for i in 0..30 {
                let jitter = 0.0001 * (i % 5) as f64 - 0.0002;
                rows.push(PinOffset {
                    run: run.to_string(),
                    calc_mass: 2000.0,
                    ..row(&format!("{run}z{i}"), false, ppm * 1e-6 * 2000.0 + jitter, 0, "PEPTIDEK")
                });
            }
            for i in 0..10 {
                let observed = true_delta + ppm * 1e-6 * (2000.0 + true_delta) + 0.0001 * (i % 3) as f64;
                rows.push(PinOffset {
                    run: run.to_string(),
                    calc_mass: 2000.0,
                    ..row(&format!("{run}p{i}"), false, observed, 2, "ASPEK")
                });
            }
        }
        rows.push(PinOffset {
            run: "b".to_string(),
            refined: true,
            ..row("bref", false, 15.9949, 3, "PEMTIDEK")
        });
        let r = report(&rows);
        assert_eq!(r.run_ppm_shifts.len(), 2);
        assert!((r.run_ppm_shifts[0].1 - 5.0).abs() < 0.2, "{:?}", r.run_ppm_shifts);
        assert!((r.run_ppm_shifts[1].1 + 5.0).abs() < 0.2, "{:?}", r.run_ppm_shifts);
        let phospho: Vec<&OffsetPeak> = r
            .peaks
            .iter()
            .filter(|p| (p.mass - 79.966).abs() < 0.05)
            .collect();
        assert_eq!(phospho.len(), 1, "{:?}", r.peaks);
        assert_eq!(phospho[0].psms, 20);
        assert!((phospho[0].mass - true_delta).abs() < 0.001, "{}", phospho[0].mass);
        let refined = r.assignments.iter().find(|a| a.psm_id == "bref").expect("refined row");
        assert!((refined.delta - 15.9949).abs() < 1e-12, "{}", refined.delta);
    }

    #[test]
    fn methionine_peak_at_plus_16_is_oxidation() {
        let mut rows = zero_rows(30, 0.0);
        rows.extend((0..12).map(|i| {
            row(&format!("o{i}"), false, 15.9949 + 0.0003 * (i % 3) as f64 - 0.0003, 3, "PEMTIDEK")
        }));
        let r = report(&rows);
        let ox = r
            .peaks
            .iter()
            .find(|p| (p.mass - 15.9949).abs() < 0.01)
            .expect("+16 peak");
        assert_eq!(ox.psms, 12);
        assert_eq!(ox.residues, vec![("M".to_string(), 12)]);
        let top = &ox.explanations[0];
        assert!(top.label.starts_with("Oxidation (UNIMOD:35)"), "{:?}", ox.explanations);
        assert!(top.posterior > 0.5, "{:?}", ox.explanations);
        let zero = r.peaks.iter().find(|p| p.mass.abs() < 0.01).expect("zero peak");
        assert_eq!(zero.explanations[0].label, "unmodified");
    }

    #[test]
    fn a_ppm_calibration_offset_is_removed_before_the_peaks_are_found() {
        let ppm = 4.0;
        let mut rows = zero_rows(40, ppm);
        rows.extend((0..12).map(|i| {
            let true_delta = 79.96633;
            let observed = true_delta + ppm * 1e-6 * (1500.0 + true_delta) + 0.0001 * (i % 3) as f64;
            row(&format!("p{i}"), false, observed, 2, "ASPEK")
        }));
        let r = report(&rows);
        assert!((r.ppm_shift - ppm).abs() < 0.2, "{}", r.ppm_shift);
        let zero = r.peaks.iter().find(|p| p.mass.abs() < 0.01).expect("zero peak");
        assert!(zero.mass.abs() < 0.0005, "{}", zero.mass);
        let phospho = r.peaks.iter().find(|p| (p.mass - 79.966).abs() < 0.01).expect("phospho");
        assert!((phospho.mass - 79.96633).abs() < 0.001, "{}", phospho.mass);
        assert!(phospho.explanations[0].label.starts_with("Phospho (UNIMOD:21)"), "{:?}", phospho.explanations);
    }

    #[test]
    fn fragment_refined_offsets_keep_their_value() {
        let ppm = 4.0;
        let mut rows = zero_rows(40, ppm);
        rows.extend((0..12).map(|i| PinOffset {
            refined: true,
            ..row(&format!("o{i}"), false, 15.9949 + 0.0001 * (i % 3) as f64, 3, "PEMTIDEK")
        }));
        let r = report(&rows);
        let ox = r.peaks.iter().find(|p| (p.mass - 15.99).abs() < 0.05).expect("+16 peak");
        assert!((ox.mass - 15.995).abs() < 0.0005, "{}", ox.mass);
    }

    #[test]
    fn deamidation_and_an_isotope_step_split_at_the_density_valley() {
        let mut sorted: Vec<f64> = Vec::new();
        for i in 0..60 {
            sorted.push(0.9840 + 0.0004 * ((i % 7) as f64 - 3.0));
            sorted.push(1.0034 + 0.0004 * ((i % 7) as f64 - 3.0));
        }
        for i in 0..8 {
            sorted.push(0.9900 + 0.0014 * i as f64);
        }
        sorted.sort_by(|a, b| a.total_cmp(b));
        let groups = offset_peaks(&sorted, 0.002);
        let big: Vec<_> = groups.iter().filter(|g| g.len() >= 40).collect();
        assert_eq!(big.len(), 2, "{groups:?}");
        let centre = |g: &std::ops::Range<usize>| median(&mut sorted[g.clone()].to_vec());
        assert!((centre(big[0]) - 0.9840).abs() < 0.002);
        assert!((centre(big[1]) - 1.0034).abs() < 0.002);
    }

    #[test]
    fn a_mixed_peak_reports_both_explanations_in_proportion() {
        let mut rows = zero_rows(30, 0.0);
        rows.extend((0..30).map(|i| {
            row(&format!("n{i}"), false, 0.9840 + 0.0005 * (i % 5) as f64 - 0.001, 3, "PENTIDEK")
        }));
        rows.extend((0..20).map(|i| {
            row(&format!("j{i}"), false, 1.0034 + 0.0005 * (i % 5) as f64 - 0.001, 0, "PEPTIDEK")
        }));
        rows.extend((0..12).map(|i| row(&format!("b{i}"), false, 0.9870 + 0.0013 * i as f64, 0, "PEPTIDEK")));
        let r = report(&rows);
        let n0_peak = r
            .assignments
            .iter()
            .find(|a| a.psm_id == "n0")
            .and_then(|a| a.peak)
            .expect("n0 is in a peak");
        let peak = r.peaks.iter().find(|p| p.id == n0_peak).expect("its peak");
        let label_of = |e: &Explanation| e.label.clone();
        let labels: Vec<String> = peak.explanations.iter().map(label_of).collect();
        assert!(labels.iter().any(|l| l.starts_with("Deamidated (UNIMOD:7)")), "{:?}", peak.explanations);
        let deamidated: Vec<&PsmAssignment> = r
            .assignments
            .iter()
            .filter(|a| a.psm_id.starts_with('n'))
            .collect();
        assert!(deamidated.iter().all(|a| a.explanation.starts_with("Deamidated (UNIMOD:7)")), "{deamidated:?}");
    }

    #[test]
    fn refined_offsets_are_judged_by_the_spread_of_refined_peaks() {
        let mut rows = zero_rows(30, 0.0);
        rows.extend((0..40).map(|i| PinOffset {
            refined: true,
            ..row(&format!("p{i}"), false, 79.9663 + 0.0012 * ((i % 9) as f64 - 4.0), 2, "ASPEK")
        }));
        let r = report(&rows);
        assert!(r.refined_precision > r.precision, "{} {}", r.refined_precision, r.precision);
        let phospho = r.peaks.iter().find(|p| (p.mass - 79.966).abs() < 0.01).expect("phospho");
        assert!(phospho.explanations[0].label.starts_with("Phospho (UNIMOD:21)"), "{:?}", phospho.explanations);
        assert!(phospho.explanations[0].posterior > 0.8, "{:?}", phospho.explanations);
    }

    #[test]
    fn a_peak_one_neutron_up_is_an_isotope_error() {
        let rows: Vec<PinOffset> = (0..10)
            .map(|i| row(&format!("i{i}"), false, 1.0034 + 0.0002 * (i % 2) as f64, 0, "PEPTIDEK"))
            .collect();
        let r = report(&rows);
        assert_eq!(r.peaks.len(), 1);
        assert!(
            r.peaks[0]
                .explanations
                .iter()
                .any(|e| e.label == "isotope error +1"),
            "{:?}",
            r.peaks[0].explanations
        );
    }

    #[test]
    fn decoys_above_the_threshold_are_counted_in_their_peak() {
        let mut rows: Vec<PinOffset> = (0..5)
            .map(|i| row(&format!("t{i}"), false, 42.0106, 1, "AKDE"))
            .collect();
        rows.push(row("d0", true, 42.0107, 1, "EDKA"));
        rows.push(row("d1", true, 80.0, 1, "EDKA"));
        let r = report(&rows);
        assert_eq!(r.peaks.len(), 1);
        assert_eq!(r.peaks[0].decoys, 1);
    }

    #[test]
    fn neighbouring_decoy_windows_do_not_overlap() {
        let windows = decoy_windows(&[(10.0, 10.0055), (10.0065, 10.012)], 0.002);
        assert!(windows[0].1 <= windows[1].0, "{windows:?}");
        assert!((windows[0].1 - 10.006).abs() < 1e-9, "{windows:?}");
        assert!((windows[0].0 - 9.998).abs() < 1e-9, "{windows:?}");
        assert!((windows[1].1 - 10.014).abs() < 1e-9, "{windows:?}");
    }

    #[test]
    fn singletons_are_not_peaks_but_keep_an_assignment_row() {
        let rows = vec![
            row("a", false, 10.0, 0, "PEPTIDE"),
            row("b", false, 55.5, 0, "PEPTIDE"),
        ];
        let r = report(&rows);
        assert!(r.peaks.is_empty());
        assert_eq!(r.assignments.len(), 2);
        assert!(r.assignments.iter().all(|a| a.peak.is_none()));
    }

    #[test]
    fn report_tables_round_trip_through_parquet() {
        let mut rows = zero_rows(30, 0.0);
        rows.extend((0..12).map(|i| {
            row(&format!("o{i}"), false, 15.9949 + 0.0003 * (i % 3) as f64 - 0.0003, 3, "PEMTIDEK")
        }));
        rows.push(row("lone", false, 55.5, 0, "PEPTIDE"));
        let r = report(&rows);
        let dir = tempfile::tempdir().unwrap();
        let peaks = dir.path().join("x.open_offsets.parquet");
        let psms = dir.path().join("x.open_psms.parquet");
        write_peaks_table(&peaks, &r).unwrap();
        write_psms_table(&psms, &r).unwrap();
        let t = read_peaks_table(&peaks).unwrap();
        let n: usize = t.batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(n, r.peaks.len());
        assert_eq!(t.meta.get(PRECISION_KEY).map(String::as_str), Some(r.precision.to_string().as_str()));
        let t = read_psms_table(&psms).unwrap();
        let n: usize = t.batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(n, r.assignments.len());
        assert!(read_peaks_table(&psms).is_err());
    }

    #[test]
    fn summary_block_counts_each_psm_once_and_folds_isotope_picks() {
        let mut rows = zero_rows(30, 0.0);
        rows.extend((0..12).map(|i| {
            row(&format!("o{i}"), false, 15.9949 + 0.0003 * (i % 3) as f64 - 0.0003, 3, "PEMTIDEK")
        }));
        rows.push(row("lone", false, 55.5, 0, "PEPTIDE"));
        let r = report(&rows);
        let s = summary_block(&r);
        assert!(s.contains("PSMs                         : 43"), "{s}");
        assert!(s.contains("unmodified (zero offset)   : 30"), "{s}");
        assert!(s.contains("outside every peak         : 1"), "{s}");
        assert!(s.contains("Oxidation (UNIMOD:35)"), "{s}");
        assert!(s.contains("Largest peaks:"), "{s}");
        assert_eq!(modification_of("Phospho (UNIMOD:21) +1 isotope"), "Phospho (UNIMOD:21)");
        assert_eq!(modification_of("isotope error -1"), ISOTOPE_PICK);
    }

    #[test]
    fn pin_columns_are_read_by_name() {
        let text = "SpecId\tLabel\tScanNr\tExpMass\tCalcMass\tisotope_error\tDeltaMass\tAbsDeltaMass\tOpenSite\tPeptide\tProteins\n\
                    s1\t1\t5\t1015.99\t1000.0\t0\t15.9949\t15.9949\t3\tK.PEMTIDEK.R\tP1\n\
                    s2\t-1\t6\t1000.001\t1000.0\t0\t0.0010\t0.0010\t0\tK.KEDITMEP.R\tdecoy_P1\n";
        let rows = parse_pin_offsets(std::io::Cursor::new(text)).unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].sequence, "PEMTIDEK");
        assert_eq!(rows[0].site, 3);
        assert_eq!(rows[0].calc_mass, 1000.0);
        assert!(rows[0].refined);
        assert!(rows[1].decoy);
        assert!(!rows[1].refined);
    }
}
