//! The precursors that are quantified: one per cross-linked peptide pair, link sites and
//! charge, taken from the target cross-link spectrum matches whose residue pair passes the
//! identification threshold, with the run and retention time of each identification.

use std::collections::{BTreeMap, HashMap};

use model::mass::PROTON;
use model::peptide::Peptide;
use quant::formula::Formula;
use search::crosslinks::{Csm, PeptideHit, XlPeptide};
use xl::crosslinker::Crosslinker;

use super::accept::{self, CohortFdr};
use super::{Search, Settings, ENVELOPE_ISOTOPES, ISOTOPE_COVERAGE, MAX_ISOTOPES, MIN_ISOTOPES};
use crate::quant::run::RunMs1;
use crate::xl::fdr::XlRun;

/// One peptide of a cross-linked pair.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Side {
    /// The peptide with its modifications, without the link.
    pub(crate) peptide: String,
    pub(crate) bare: String,
    /// The linked residue, 1-based within the peptide.
    pub(crate) site: usize,
    /// `nterm` or `cterm` for a link on a terminus, else the linked residue.
    pub(crate) site_kind: String,
}

#[derive(Debug, Clone)]
pub(crate) struct Precursor {
    pub(crate) key: String,
    pub(crate) linker: String,
    pub(crate) sides: [Side; 2],
    pub(crate) charge: u8,
    /// The theoretical monoisotopic m/z.
    pub(crate) mono_mz: f64,
    /// The traced isotopes, summing to one.
    pub(crate) envelope: Vec<f64>,
    /// The share of the whole isotope envelope that the traced isotopes hold.
    pub(crate) isotope_share: f64,
    pub(crate) residue_pair: String,
    pub(crate) protein_pair: String,
    pub(crate) link: &'static str,
    /// The q-value of the residue pair.
    pub(crate) q: f64,
}

/// The best cross-link spectrum match of a precursor in a run.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Identification {
    pub(crate) precursor: usize,
    pub(crate) run: usize,
    pub(crate) rt: f64,
    pub(crate) scan: i32,
    pub(crate) score: f64,
    pub(crate) n_csms: u32,
}

pub(crate) struct Library {
    pub(crate) precursors: Vec<Precursor>,
    /// One per precursor and run, ordered by precursor and run.
    pub(crate) identifications: Vec<Identification>,
}

/// The peptide with its modifications as `to_proforma` writes it, without the flanking residues.
fn peptidoform(p: &Peptide) -> String {
    let s = p.to_proforma();
    s.get(2..s.len().saturating_sub(2))
        .unwrap_or(&s)
        .to_string()
}

fn side_of(hit: &PeptideHit, p: &XlPeptide) -> Side {
    let site_kind = if hit.site.n_term {
        "nterm".to_string()
    } else if hit.site.c_term {
        "cterm".to_string()
    } else {
        p.peptide
            .residues
            .get(hit.site.pos)
            .map_or("NA".to_string(), |aa| (aa.residue as char).to_string())
    };
    Side {
        peptide: peptidoform(&p.peptide),
        bare: String::from_utf8_lossy(&p.bare()).into_owned(),
        site: hit.site.pos + 1,
        site_kind,
    }
}

fn side_key(s: &Side) -> String {
    let terminus = match s.site_kind.as_str() {
        "nterm" => "n",
        "cterm" => "c",
        _ => "",
    };
    format!("{}@{}{terminus}", s.peptide, s.site)
}

pub(super) fn peptide_formula(p: &Peptide) -> Formula {
    Formula::peptide(
        p.residues
            .iter()
            .map(|aa| (aa.residue, aa.mod_.map_or(0.0, |m| m.mass_delta))),
    )
}

/// The elemental formula of a cross-linked pair: both peptides and the intact link. A linker
/// that is defined by a mass alone adds averagine atoms of that mass.
fn pair_formula(a: &Peptide, b: &Peptide, linker: &Crosslinker) -> Formula {
    let mut f = peptide_formula(a);
    f.add(peptide_formula(b));
    match &linker.composition {
        Some(e) => f.add_elements(e),
        None if linker.intact >= 0.0 => f.add(Formula::averagine(linker.intact)),
        None => f.sub_saturating(Formula::averagine(-linker.intact)),
    }
    f
}

/// The leading isotopes of `full` that hold `ISOTOPE_COVERAGE` of it (between `MIN_ISOTOPES`
/// and `MAX_ISOTOPES`), scaled to sum to one, and the share of `full` they hold.
pub(super) fn traced_envelope(full: &[f64]) -> (Vec<f64>, f64) {
    let total: f64 = full.iter().sum();
    if full.is_empty() || total.is_nan() || total <= 0.0 {
        return (vec![1.0], 1.0);
    }
    let mut k = full.len();
    let mut cum = 0.0;
    for (i, &v) in full.iter().enumerate() {
        cum += v;
        if i + 1 >= MIN_ISOTOPES && cum >= ISOTOPE_COVERAGE * total {
            k = i + 1;
            break;
        }
    }
    let k = k.clamp(MIN_ISOTOPES, MAX_ISOTOPES).min(full.len());
    let kept: f64 = full[..k].iter().sum();
    if kept <= 0.0 {
        return (vec![1.0], 1.0);
    }
    (full[..k].iter().map(|v| v / kept).collect(), kept / total)
}

fn precursor_of(
    csm: &Csm,
    run: &XlRun,
    search: &Search<'_>,
    residue_pair: &str,
    link: &'static str,
    q: f64,
) -> Precursor {
    let table = &run.table;
    let pa = &table.peptides()[csm.alpha.pep as usize];
    let pb = &table.peptides()[csm.beta.pep as usize];
    let mut sides = [side_of(&csm.alpha, pa), side_of(&csm.beta, pb)];
    if side_key(&sides[1]) < side_key(&sides[0]) {
        sides.swap(0, 1);
    }
    let charge = csm.charge.max(1);
    let neutral = pa.peptide.mass() + pb.peptide.mass() + run.linker.intact;
    let (envelope, isotope_share) = traced_envelope(
        &pair_formula(&pa.peptide, &pb.peptide, &run.linker).envelope(ENVELOPE_ISOTOPES),
    );
    let (_, protein_pair) = crate::xl::fdr::pair_keys(csm, table, search.prepared);
    Precursor {
        key: format!(
            "{}|{}|{}/{charge}",
            run.linker.name,
            side_key(&sides[0]),
            side_key(&sides[1])
        ),
        linker: run.linker.name.clone(),
        sides,
        charge,
        mono_mz: neutral / f64::from(charge) + PROTON,
        envelope,
        isotope_share,
        residue_pair: residue_pair.to_string(),
        protein_pair,
        link,
        q,
    }
}

pub(super) fn build(
    search: &Search<'_>,
    cohort: Option<&CohortFdr>,
    runs: &[RunMs1<'_>],
    settings: &Settings,
) -> Library {
    let levels = search.levels;
    let qvalues = accept::qvalues(search, cohort);
    let mut by_key: BTreeMap<String, Precursor> = BTreeMap::new();
    let mut found: Vec<(String, Identification)> = Vec::new();
    let (mut no_rt, mut no_run) = (0usize, 0usize);
    for ((item, &row), q) in levels.csm_items.iter().zip(&levels.csm_rows).zip(qvalues) {
        let Some((match_q, q)) = q else {
            continue;
        };
        if match_q > settings.xl_fdr || q > settings.xl_fdr {
            continue;
        }
        let result = &search.results[row];
        let (Some(csm), Some(spec)) = (&result.csm, search.spectra.get(result.spec_idx)) else {
            continue;
        };
        let Some(rt) = spec.rt_seconds.filter(|t| t.is_finite()) else {
            no_rt += 1;
            continue;
        };
        let Some(run) = runs.iter().position(|r| r.span.contains(&result.spec_idx)) else {
            no_run += 1;
            continue;
        };
        let p = precursor_of(
            csm,
            &search.linkers[result.linker],
            search,
            &item.key,
            item.link.name(),
            q,
        );
        let key = p.key.clone();
        by_key.entry(key.clone()).or_insert(p);
        found.push((
            key,
            Identification {
                precursor: 0,
                run,
                rt,
                scan: spec.scan.unwrap_or(0),
                score: search.score[row].unwrap_or(f64::NEG_INFINITY),
                n_csms: 1,
            },
        ));
    }
    if no_rt + no_run > 0 {
        log::warn!(
            target: "xl-quant",
            "{} accepted cross-link spectrum matches are not quantified (no retention time: \
             {no_rt}; no MS1 scans of their run: {no_run})",
            no_rt + no_run
        );
    }
    let index: HashMap<&str, usize> = by_key
        .keys()
        .enumerate()
        .map(|(i, k)| (k.as_str(), i))
        .collect();
    let mut best: BTreeMap<(usize, usize), Identification> = BTreeMap::new();
    for (key, mut id) in found {
        id.precursor = index[key.as_str()];
        best.entry((id.precursor, id.run))
            .and_modify(|kept| {
                let n = kept.n_csms + 1;
                if id.score > kept.score || (id.score == kept.score && id.scan < kept.scan) {
                    *kept = id.clone();
                }
                kept.n_csms = n;
            })
            .or_insert(id);
    }
    drop(index);
    Library {
        precursors: by_key.into_values().collect(),
        identifications: best.into_values().collect(),
    }
}

#[cfg(test)]
#[path = "tests/library.rs"]
mod tests;
