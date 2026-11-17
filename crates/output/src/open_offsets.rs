use std::collections::{HashMap, HashSet};
use std::io::{BufRead, Write};
use std::path::Path;

use model::unimod::UnimodDb;

use crate::percolator::PercolatorPsm;

const ISOTOPE_DA: f64 = 1.003_355;
const MIN_PEAK_PSMS: usize = 2;
const DEFAULT_PRECISION_DA: f64 = 0.002;
const MIN_PRECISION_DA: f64 = 0.0005;
const ZERO_WINDOW_DA: f64 = 0.02;
const SPLIT_SIGMAS: f64 = 4.0;
const CANDIDATE_SIGMAS: f64 = 6.0;
const PRIOR_SINGLE: f64 = 1.0;
const PRIOR_ISOTOPE: f64 = 0.2;
const PRIOR_SINGLE_ISOTOPE: f64 = 0.1;
const PRIOR_PAIR: f64 = 0.05;
const PRIOR_UNEXPLAINED: f64 = 0.05;
const PAIR_SOURCE_POSTERIOR: f64 = 0.5;
const TOP_EXPLANATIONS: usize = 3;

#[derive(Debug, Clone, PartialEq)]
pub struct PinOffset {
    pub psm_id: String,
    pub decoy: bool,
    pub delta: f64,
    pub site: usize,
    pub sequence: String,
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
    pub lo: f64,
    pub hi: f64,
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
    pub precision: f64,
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
    let (id_c, label_c, delta_c, site_c, pep_c) = (
        col("SpecId")?,
        col("Label")?,
        col("DeltaMass")?,
        col("OpenSite")?,
        col("Peptide")?,
    );
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
        let delta: f64 = get(delta_c)?
            .parse()
            .map_err(|e| bad(format!("DeltaMass: {e}")))?;
        let site: usize = get(site_c)?
            .parse()
            .map_err(|e| bad(format!("OpenSite: {e}")))?;
        out.push(PinOffset {
            psm_id: get(id_c)?.to_string(),
            decoy: label == "-1",
            delta,
            site,
            sequence: stripped_sequence(get(pep_c)?),
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

fn mad(values: &mut [f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_by(|a, b| a.total_cmp(b));
    let med = values[values.len() / 2];
    let mut dev: Vec<f64> = values.iter().map(|v| (v - med).abs()).collect();
    dev.sort_by(|a, b| a.total_cmp(b));
    1.4826 * dev[dev.len() / 2]
}

fn median(values: &mut [f64]) -> f64 {
    values.sort_by(|a, b| a.total_cmp(b));
    values[values.len() / 2]
}

pub fn measured_precision(deltas: &[f64]) -> f64 {
    let mut near: Vec<f64> = deltas
        .iter()
        .copied()
        .filter(|d| d.abs() < ZERO_WINDOW_DA)
        .collect();
    if near.len() < 20 {
        return DEFAULT_PRECISION_DA;
    }
    mad(&mut near).max(MIN_PRECISION_DA)
}

struct Candidate {
    label: String,
    mass: f64,
    prior: f64,
    sites: Option<Vec<String>>,
}

fn gauss_log(x: f64, mu: f64, s: f64) -> f64 {
    let z = (x - mu) / s;
    -0.5 * z * z - s.ln() - 0.5 * (2.0 * std::f64::consts::PI).ln()
}

fn site_support(sites: &[String], localised: &[(String, bool, bool)]) -> f64 {
    if localised.is_empty() {
        return 1.0;
    }
    let allowed = localised
        .iter()
        .filter(|(residue, n_term, c_term)| {
            sites.iter().any(|s| {
                s == residue || (s == "N-term" && *n_term) || (s == "C-term" && *c_term)
            })
        })
        .count();
    (allowed as f64 + 0.5) / (localised.len() as f64 + 1.0)
}

fn single_candidates(db: &UnimodDb) -> Vec<Candidate> {
    let mut out: Vec<Candidate> = db
        .mods()
        .filter(|m| m.mono_mass.is_finite() && m.mono_mass != 0.0)
        .map(|m| Candidate {
            label: format!("{} ({})", m.name, m.accession),
            mass: m.mono_mass,
            prior: PRIOR_SINGLE,
            sites: Some(m.sites.clone()),
        })
        .collect();
    out.sort_by(|a, b| a.mass.total_cmp(&b.mass).then_with(|| a.label.cmp(&b.label)));
    out
}

fn explain(
    mass: f64,
    s: f64,
    span: f64,
    localised: &[(String, bool, bool)],
    singles: &[Candidate],
    pairs: &[Candidate],
) -> Vec<Explanation> {
    let reach = CANDIDATE_SIGMAS * s;
    let mut scored: Vec<(String, f64, f64)> = Vec::new();
    let mut push = |label: String, cmass: f64, prior: f64, support: f64| {
        if (mass - cmass).abs() <= reach {
            scored.push((label, cmass, prior.ln() + gauss_log(mass, cmass, s) + support.ln()));
        }
    };
    push("unmodified".to_string(), 0.0, PRIOR_SINGLE, 1.0);
    for k in [-2i32, -1, 1, 2] {
        push(
            format!("isotope error {k:+}"),
            f64::from(k) * ISOTOPE_DA,
            PRIOR_ISOTOPE,
            1.0,
        );
    }
    let lo = singles.partition_point(|c| c.mass < mass - reach - 2.0 * ISOTOPE_DA);
    for c in &singles[lo..] {
        if c.mass > mass + reach + 2.0 * ISOTOPE_DA {
            break;
        }
        let support = c
            .sites
            .as_deref()
            .map_or(1.0, |sites| site_support(sites, localised));
        push(c.label.clone(), c.mass, c.prior, support);
        for k in [-1i32, 1] {
            push(
                format!("{} {k:+} isotope", c.label),
                c.mass + f64::from(k) * ISOTOPE_DA,
                PRIOR_SINGLE_ISOTOPE,
                support,
            );
        }
    }
    for c in pairs {
        let support = c
            .sites
            .as_deref()
            .map_or(1.0, |sites| site_support(sites, localised));
        push(c.label.clone(), c.mass, c.prior, support);
    }
    let unexplained = PRIOR_UNEXPLAINED.ln() - span.max(1.0).ln();
    let max = scored
        .iter()
        .map(|(_, _, l)| *l)
        .fold(unexplained, f64::max);
    let total: f64 = scored.iter().map(|(_, _, l)| (l - max).exp()).sum::<f64>()
        + (unexplained - max).exp();
    let mut out: Vec<Explanation> = scored
        .into_iter()
        .map(|(label, cmass, l)| Explanation {
            label,
            mass: cmass,
            posterior: (l - max).exp() / total,
        })
        .collect();
    out.push(Explanation {
        label: "unexplained".to_string(),
        mass,
        posterior: (unexplained - max).exp() / total,
    });
    out.sort_by(|a, b| {
        b.posterior
            .total_cmp(&a.posterior)
            .then_with(|| a.label.cmp(&b.label))
    });
    out.truncate(TOP_EXPLANATIONS);
    out
}

struct Cluster {
    members: Vec<usize>,
}

struct RawPeak {
    cluster: Cluster,
    mass: f64,
    width: f64,
    localised: Vec<(String, bool, bool)>,
}

fn clusters(rows: &[&PinOffset], precision: f64) -> Vec<Cluster> {
    let mut order: Vec<usize> = (0..rows.len()).collect();
    order.sort_by(|&a, &b| rows[a].delta.total_cmp(&rows[b].delta));
    let mut out: Vec<Cluster> = Vec::new();
    let mut cur: Vec<usize> = Vec::new();
    for i in order {
        if let Some(&last) = cur.last() {
            if rows[i].delta - rows[last].delta > SPLIT_SIGMAS * precision {
                out.push(Cluster {
                    members: std::mem::take(&mut cur),
                });
            }
        }
        cur.push(i);
    }
    if !cur.is_empty() {
        out.push(Cluster { members: cur });
    }
    out
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
    let deltas: Vec<f64> = targets.iter().map(|r| r.delta).collect();
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
        .map(|r| r.delta)
        .collect();
    decoy_deltas.sort_by(|a, b| a.total_cmp(b));

    let singles = single_candidates(db);
    let mut raw: Vec<RawPeak> = clusters(&targets, precision)
        .into_iter()
        .filter(|c| c.members.len() >= MIN_PEAK_PSMS)
        .map(|c| {
            let mut ds: Vec<f64> = c.members.iter().map(|&i| targets[i].delta).collect();
            let width = mad(&mut ds.clone()).max(precision);
            let mass = median(&mut ds);
            let localised: Vec<(String, bool, bool)> = c
                .members
                .iter()
                .filter_map(|&i| {
                    let r = targets[i];
                    let residue = site_label(&r.sequence, r.site)?;
                    Some((residue, r.site == 1, r.site == r.sequence.len()))
                })
                .collect();
            RawPeak {
                cluster: c,
                mass,
                width,
                localised,
            }
        })
        .collect();
    raw.sort_by(|a, b| {
        b.cluster
            .members
            .len()
            .cmp(&a.cluster.members.len())
            .then_with(|| a.mass.total_cmp(&b.mass))
    });

    let explain_peak = |mass: f64, width: f64, n: usize, localised: &[(String, bool, bool)], pairs: &[Candidate]| {
        let s = (width * width / n as f64 + precision * precision).sqrt();
        explain(mass, s, span, localised, &singles, pairs)
    };
    let first: Vec<Vec<Explanation>> = raw
        .iter()
        .map(|p| explain_peak(p.mass, p.width, p.cluster.members.len(), &p.localised, &[]))
        .collect();
    let mut sources: Vec<&Candidate> = Vec::new();
    for ex in &first {
        if let Some(top) = ex.first() {
            if top.posterior >= PAIR_SOURCE_POSTERIOR {
                if let Some(c) = singles.iter().find(|c| c.label == top.label) {
                    if !sources.iter().any(|s| s.label == c.label) {
                        sources.push(c);
                    }
                }
            }
        }
    }
    let mut pairs: Vec<Candidate> = Vec::new();
    for (i, a) in sources.iter().enumerate() {
        for b in &sources[i..] {
            let sites: Vec<String> = a
                .sites
                .iter()
                .flatten()
                .chain(b.sites.iter().flatten())
                .cloned()
                .collect();
            pairs.push(Candidate {
                label: format!("{} + {}", a.label, b.label),
                mass: a.mass + b.mass,
                prior: PRIOR_PAIR,
                sites: Some(sites),
            });
        }
    }

    let mut peaks = Vec::with_capacity(raw.len());
    let mut assignments = Vec::with_capacity(targets.len());
    let mut assigned: HashSet<usize> = HashSet::new();
    for (k, p) in raw.iter().enumerate() {
        let (c, mass, width, localised) = (&p.cluster, &p.mass, &p.width, &p.localised);
        let id = k + 1;
        let explanations = explain_peak(*mass, *width, c.members.len(), localised, &pairs);
        let lo = c.members.iter().map(|&i| targets[i].delta).fold(f64::INFINITY, f64::min);
        let hi = c.members.iter().map(|&i| targets[i].delta).fold(f64::NEG_INFINITY, f64::max);
        let (dlo, dhi) = (lo - SPLIT_SIGMAS * precision / 2.0, hi + SPLIT_SIGMAS * precision / 2.0);
        let decoys = decoy_deltas.partition_point(|d| *d <= dhi) - decoy_deltas.partition_point(|d| *d < dlo);
        let peptides: HashSet<&str> = c.members.iter().map(|&i| targets[i].sequence.as_str()).collect();
        let mut residue_counts: HashMap<&str, usize> = HashMap::new();
        for (r, _, _) in localised {
            *residue_counts.entry(r.as_str()).or_default() += 1;
        }
        let mut residues: Vec<(String, usize)> = residue_counts
            .into_iter()
            .map(|(r, n)| (r.to_string(), n))
            .collect();
        residues.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        let top = explanations.first().cloned();
        for &i in &c.members {
            assigned.insert(i);
            let r = targets[i];
            assignments.push(PsmAssignment {
                psm_id: r.psm_id.clone(),
                delta: r.delta,
                peak: Some(id),
                residue: site_label(&r.sequence, r.site).unwrap_or_default(),
                explanation: top.as_ref().map_or_else(String::new, |t| t.label.clone()),
                posterior: top.as_ref().map_or(0.0, |t| t.posterior),
            });
        }
        peaks.push(OffsetPeak {
            id,
            mass: *mass,
            width: *width,
            lo,
            hi,
            psms: c.members.len(),
            peptides: peptides.len(),
            decoys,
            residues,
            explanations,
        });
    }
    for (i, r) in targets.iter().enumerate() {
        if !assigned.contains(&i) {
            assignments.push(PsmAssignment {
                psm_id: r.psm_id.clone(),
                delta: r.delta,
                peak: None,
                residue: site_label(&r.sequence, r.site).unwrap_or_default(),
                explanation: String::new(),
                posterior: 0.0,
            });
        }
    }
    assignments.sort_by(|a, b| a.psm_id.cmp(&b.psm_id));
    OffsetReport {
        precision,
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

pub fn write_peaks<W: Write>(w: &mut W, report: &OffsetReport) -> std::io::Result<()> {
    write!(
        w,
        "peak\tmass_shift\twidth_da\tpsms\tpeptides\tdecoys\tdecoy_ratio\tresidues"
    )?;
    for k in 1..=TOP_EXPLANATIONS {
        write!(w, "\texplanation_{k}\tmass_{k}\tposterior_{k}")?;
    }
    writeln!(w)?;
    for p in &report.peaks {
        write!(
            w,
            "{}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{:.4}\t{}",
            p.id,
            p.mass,
            p.width,
            p.psms,
            p.peptides,
            p.decoys,
            p.decoys as f64 / p.psms as f64,
            residues_field(&p.residues)
        )?;
        for k in 0..TOP_EXPLANATIONS {
            match p.explanations.get(k) {
                Some(e) => write!(w, "\t{}\t{:.4}\t{:.3}", e.label, e.mass, e.posterior)?,
                None => write!(w, "\t\t\t")?,
            }
        }
        writeln!(w)?;
    }
    Ok(())
}

pub fn write_assignments<W: Write>(w: &mut W, report: &OffsetReport) -> std::io::Result<()> {
    writeln!(w, "PSMId\tDeltaMass\tpeak\tresidue\texplanation\tposterior")?;
    for a in &report.assignments {
        writeln!(
            w,
            "{}\t{:.4}\t{}\t{}\t{}\t{:.3}",
            a.psm_id,
            a.delta,
            a.peak.map_or_else(|| "NA".to_string(), |p| p.to_string()),
            a.residue,
            a.explanation,
            a.posterior
        )?;
    }
    Ok(())
}

pub fn summary_lines(report: &OffsetReport, n: usize) -> Vec<String> {
    report
        .peaks
        .iter()
        .take(n)
        .map(|p| {
            let top = p.explanations.first();
            format!(
                "  {:+10.4} Da  {:>6} PSMs  {:>5} decoys  {}  {}",
                p.mass,
                p.psms,
                p.decoys,
                residues_field(&p.residues),
                top.map_or_else(String::new, |t| format!("{} ({:.2})", t.label, t.posterior))
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(id: &str, decoy: bool, delta: f64, site: usize, seq: &str) -> PinOffset {
        PinOffset {
            psm_id: id.to_string(),
            decoy,
            delta,
            site,
            sequence: seq.to_string(),
        }
    }

    fn db() -> &'static UnimodDb {
        model::unimod::bundled()
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
        build_report(rows, &accepted, &target_scores, &decoy_scores, db())
    }

    #[test]
    fn stripped_sequence_drops_flanks_and_mods() {
        assert_eq!(stripped_sequence("K.PEM[+15.995]TIDEK.R"), "PEMTIDEK");
        assert_eq!(stripped_sequence("-.ACDK.-"), "ACDK");
        assert_eq!(stripped_sequence("PEPTIDE"), "PEPTIDE");
    }

    #[test]
    fn methionine_peak_at_plus_16_is_oxidation() {
        let mut rows: Vec<PinOffset> = (0..30)
            .map(|i| row(&format!("z{i}"), false, 0.0001 * f64::from(i % 5) - 0.0002, 0, "PEPTIDEK"))
            .collect();
        rows.extend((0..12).map(|i| {
            row(&format!("o{i}"), false, 15.9949 + 0.0003 * f64::from(i % 3) - 0.0003, 3, "PEMTIDEK")
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
    fn a_peak_one_neutron_up_is_an_isotope_error() {
        let rows: Vec<PinOffset> = (0..10)
            .map(|i| row(&format!("i{i}"), false, 1.0034 + 0.0002 * f64::from(i % 2), 0, "PEPTIDEK"))
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
    fn pin_columns_are_read_by_name() {
        let text = "SpecId\tLabel\tScanNr\tDeltaMass\tAbsDeltaMass\tOpenSite\tPeptide\tProteins\n\
                    s1\t1\t5\t15.9949\t15.9949\t3\tK.PEMTIDEK.R\tP1\n\
                    s2\t-1\t6\t0.0010\t0.0010\t0\tK.KEDITMEP.R\tdecoy_P1\n";
        let rows = parse_pin_offsets(std::io::Cursor::new(text)).unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0], row("s1", false, 15.9949, 3, "PEMTIDEK"));
        assert!(rows[1].decoy);
    }
}
