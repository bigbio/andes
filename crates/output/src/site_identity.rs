//! Which explanation of a spectrum is right, after rescoring: a posterior
//! over every row of the spectrum (targets and decoys), so a spectrum's
//! explanations compete with each other.
//!
//! A row's logit combines its rescored score, its precursor isotope offset,
//! its fragment evidence and the signatures of its modifications
//! ([`search::site_evidence`]). The weights of the first three are fixed; they
//! were fitted on synthetic peptides of known modification. How much a
//! diagnostic ion or a neutral loss says is learned from the run itself: among
//! the rows the other evidence already makes confident, how often the
//! signature is present when the row carries the modification and when a
//! competing row of the same spectrum carries it instead.
//!
//! The site tables sum these posteriors per site: the share of the posterior
//! on a site's peptides that places a modification there (occupancy) and,
//! given a modification, the share of each one (identity).

use std::collections::HashMap;

use search::site_evidence::SiteEvidence;

use crate::percolator::PercolatorPsm;

/// Weights of the fixed part of a row's logit, per unit of each feature.
const W_SCORE: f64 = 2.0;
const W_ISOTOPE_1: f64 = -1.0;
const W_ISOTOPE_2: f64 = -5.7;
const W_HIT_FRACTION: f64 = 25.0;
const W_UNIQUE_HIT: f64 = 0.08;
const W_UNIQUE_MISS: f64 = -0.02;
const W_DECOY: f64 = -11.0;
/// Posterior above which a target row teaches the signature likelihoods.
const CONFIDENT: f64 = 0.9;

/// What the search kept of one emitted row.
#[derive(Debug, Clone)]
pub struct RowEvidence {
    /// Spectrum index: rows of one spectrum compete.
    pub spectrum: u32,
    pub isotope_offset: i8,
    pub evidence: Option<Box<SiteEvidence>>,
}

/// A modification signature: (mod name, mass bits, diagnostic or loss).
type Signature = (String, u64, bool);

/// Present and total counts, with one pseudo-count each way.
#[derive(Debug, Clone, Copy)]
struct Rate {
    present: f64,
    n: f64,
}

impl Rate {
    fn new() -> Self {
        Rate {
            present: 1.0,
            n: 2.0,
        }
    }
    fn p(&self) -> f64 {
        self.present / self.n
    }
}

struct Row<'a> {
    id: &'a str,
    decoy: bool,
    base: f64,
    evidence: Option<&'a SiteEvidence>,
}

fn base_logit(score: f64, decoy: bool, row: &RowEvidence) -> f64 {
    let mut l = W_SCORE * score + if decoy { W_DECOY } else { 0.0 };
    l += match row.isotope_offset.unsigned_abs() {
        0 => 0.0,
        1 => W_ISOTOPE_1,
        _ => W_ISOTOPE_2,
    };
    if let Some(e) = row.evidence.as_deref() {
        l += W_HIT_FRACTION * e.hit_fraction as f64
            + W_UNIQUE_HIT * e.unique_hit as f64
            + W_UNIQUE_MISS * e.unique_miss as f64;
    }
    l
}

fn softmax(logits: &[f64]) -> Vec<f64> {
    let mx = logits.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let e: Vec<f64> = logits.iter().map(|l| (l - mx).exp()).collect();
    let z: f64 = e.iter().sum();
    e.into_iter().map(|x| x / z).collect()
}

fn signatures(e: &SiteEvidence) -> impl Iterator<Item = (Signature, bool)> + '_ {
    e.marks.iter().flat_map(|m| {
        let key = |diag: bool| {
            (
                m.modification.name.clone(),
                m.modification.mass_delta.to_bits(),
                diag,
            )
        };
        m.diagnostic
            .map(|p| (key(true), p))
            .into_iter()
            .chain(m.loss.map(|p| (key(false), p)))
    })
}

/// Each target row's posterior among its spectrum's rows. Rows the search
/// kept no evidence for (spectra without a variable modification) are left
/// out.
pub fn row_posteriors(
    targets: &HashMap<String, PercolatorPsm>,
    decoys: &HashMap<String, PercolatorPsm>,
    rows: &HashMap<String, RowEvidence>,
) -> HashMap<String, f64> {
    let mut spectra: HashMap<u32, Vec<Row<'_>>> = HashMap::new();
    for (map, decoy) in [(targets, false), (decoys, true)] {
        for (id, psm) in map {
            let (Some(r), Some(score)) = (rows.get(id), psm.score) else {
                continue;
            };
            spectra.entry(r.spectrum).or_default().push(Row {
                id,
                decoy,
                base: base_logit(score, decoy, r),
                evidence: r.evidence.as_deref(),
            });
        }
    }
    // A fixed row order per spectrum keeps the result independent of hashing.
    let mut spectra: Vec<(u32, Vec<Row<'_>>)> = spectra.into_iter().collect();
    spectra.sort_unstable_by_key(|(s, _)| *s);
    for (_, rows) in &mut spectra {
        rows.sort_unstable_by(|a, b| a.id.cmp(b.id));
    }

    // Signature likelihoods from the confident target rows: P(present | the
    // row carries the mod) and P(present | another row of the spectrum
    // carries it).
    let mut carried: HashMap<Signature, Rate> = HashMap::new();
    let mut other: HashMap<Signature, Rate> = HashMap::new();
    for (_, rows) in &spectra {
        let post = softmax(&rows.iter().map(|r| r.base).collect::<Vec<_>>());
        let shown: HashMap<Signature, bool> = rows
            .iter()
            .filter_map(|r| r.evidence)
            .flat_map(signatures)
            .fold(HashMap::new(), |mut acc, (k, p)| {
                *acc.entry(k).or_insert(false) |= p;
                acc
            });
        for (r, &p) in rows.iter().zip(&post) {
            if r.decoy || p < CONFIDENT {
                continue;
            }
            let own: HashMap<Signature, bool> = r
                .evidence
                .map_or_else(HashMap::new, |e| signatures(e).collect());
            for (k, &present) in &own {
                let c = carried.entry(k.clone()).or_insert_with(Rate::new);
                c.present += present as u8 as f64;
                c.n += 1.0;
            }
            for (k, &present) in &shown {
                if !own.contains_key(k) {
                    let c = other.entry(k.clone()).or_insert_with(Rate::new);
                    c.present += present as u8 as f64;
                    c.n += 1.0;
                }
            }
        }
    }
    let llr = |k: &Signature, present: bool| -> f64 {
        let (p1, p0) = (
            carried.get(k).copied().unwrap_or_else(Rate::new).p(),
            other.get(k).copied().unwrap_or_else(Rate::new).p(),
        );
        if present {
            (p1 / p0).ln()
        } else {
            ((1.0 - p1) / (1.0 - p0)).ln()
        }
    };

    let mut out = HashMap::new();
    for (_, rows) in &spectra {
        let logits: Vec<f64> = rows
            .iter()
            .map(|r| {
                r.base
                    + r.evidence
                        .map_or(0.0, |e| signatures(e).map(|(k, p)| llr(&k, p)).sum::<f64>())
            })
            .collect();
        for (r, p) in rows.iter().zip(softmax(&logits)) {
            if !r.decoy {
                out.insert(r.id.to_string(), p);
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::modification::{leak_mod, Modification};
    use search::site_evidence::ModMark;

    fn psm(id: &str, score: f64) -> (String, PercolatorPsm) {
        (
            id.into(),
            PercolatorPsm {
                psm_id: id.into(),
                q_value: 0.001,
                pep: 0.0,
                score: Some(score),
                peptide: String::new(),
                proteins: String::new(),
            },
        )
    }

    fn ev(spectrum: u32, m: &'static Modification, diagnostic: bool) -> RowEvidence {
        RowEvidence {
            spectrum,
            isotope_offset: 0,
            evidence: Some(Box::new(SiteEvidence {
                unique_hit: 0,
                unique_miss: 0,
                hit_fraction: 0.5,
                marks: vec![ModMark {
                    modification: m,
                    diagnostic: Some(diagnostic),
                    loss: None,
                }],
            })),
        }
    }

    #[test]
    fn rows_of_a_spectrum_share_one_posterior_and_signatures_move_it() {
        let ac = leak_mod(
            Modification::from_mods_txt_line("42.010565,K,opt,any,Acetyl,diagnostic=143.1179")
                .unwrap(),
        );
        let me3 = leak_mod(
            Modification::from_mods_txt_line("42.04695,K,opt,any,Trimethyl,diagnostic=143.1543")
                .unwrap(),
        );
        // 200 spectra teach the likelihoods: a clear acetyl row whose ion is
        // present, against a weak trimethyl row whose ion is absent.
        let mut targets = HashMap::new();
        let mut rows = HashMap::new();
        for s in 0..200u32 {
            let (a, t) = (format!("a{s}"), format!("t{s}"));
            targets.extend([psm(&a, 3.0), psm(&t, 0.0)]);
            rows.insert(a, ev(s, ac, true));
            rows.insert(t, ev(s, me3, false));
        }
        // The spectrum under test: equal scores, the acetyl ion present.
        targets.extend([psm("a", 1.0), psm("t", 1.0)]);
        rows.insert("a".into(), ev(999, ac, true));
        rows.insert("t".into(), ev(999, me3, false));
        let post = row_posteriors(&targets, &HashMap::new(), &rows);
        assert!((post["a"] + post["t"] - 1.0).abs() < 1e-9);
        // Present when carried in all 200, never shown by a confident rival:
        // ln(201/202 / 0.5) each way, so P(acetyl) = 1 / (1 + e^-1.38) ≈ 0.80.
        assert!((post["a"] - 0.799).abs() < 0.01, "{}", post["a"]);
        // Without signatures the two rows tie.
        for r in rows.values_mut() {
            r.evidence.as_mut().unwrap().marks.clear();
        }
        let post = row_posteriors(&targets, &HashMap::new(), &rows);
        assert!((post["a"] - 0.5).abs() < 1e-9);
    }

    #[test]
    fn a_decoy_row_takes_its_share() {
        let mut targets = HashMap::from([psm("t", 1.0)]);
        let decoys = HashMap::from([psm("d", 1.0)]);
        let row = |s| RowEvidence {
            spectrum: s,
            isotope_offset: 0,
            evidence: None,
        };
        let rows = HashMap::from([("t".to_string(), row(0)), ("d".to_string(), row(0))]);
        let p = row_posteriors(&targets, &decoys, &rows)["t"];
        assert!(p > 0.99 && p < 1.0, "{p}");
        targets.insert("x".into(), targets["t"].clone());
        assert!(!row_posteriors(&targets, &decoys, &rows).contains_key("x"));
    }
}
