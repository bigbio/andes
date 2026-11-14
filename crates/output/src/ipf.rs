//! IPF site localization after rescoring (`--localize ipf`).
//!
//! Each identified target PSM (rescored, within the run's PSM cutoff) is
//! localized with the PSM's PEP as the prior (see [`search::ipf`]). Its
//! [`search::ipf::IpfEvidence`] is built after rescoring, from the retained
//! peaks, so every row kind (main, `--refine` Pass 2, `--chimeric`
//! secondaries) gets it the same way. The ion PEP model is fit once per run on
//! the evidence of confident target PSMs. Peptidoform q-values are model-based
//! over every localized PSM, as in PyProphet IPF.

use std::collections::HashMap;
use std::io::Write;

use search::ipf::{localize, model_fdr, IonPepModel, IpfEvidence};

use crate::percolator::PercolatorPsm;

/// PSM q-value below which a target's ion evidence trains the ion PEP model.
pub const TRAIN_Q: f64 = 0.01;

/// One localized PSM.
#[derive(Debug, Clone)]
pub struct IpfRow {
    pub psm_id: String,
    pub peptide: String,
    pub psm_q: f64,
    pub psm_pep: f64,
    /// The reported (search) placement: 1-based residue and its mod name.
    pub reported_sites: Vec<(u16, String)>,
    /// The most probable placement: 1-based residue and its mod name.
    pub best_sites: Vec<(u16, String)>,
    /// Localization probability of the best placement given the PSM is right.
    pub best_prob: f64,
    /// Peptidoform PEP (1 - posterior of the best placement, h0 included).
    pub pep: f64,
    pub q: f64,
    pub n_forms: usize,
    pub site_ions: usize,
    pub matched_ions: usize,
    pub truncated: bool,
    /// Localization probability of each (1-based residue, mod name) that some
    /// placement modifies, in residue order: the posterior mass of the
    /// placements that put that mod on that residue.
    pub site_probs: Vec<(u16, String, f64)>,
    /// Debug (`--ipf-debug`): every placement as
    /// `sites|search_score|ipf_log_evidence`, `;`-separated. Empty otherwise.
    pub forms_debug: String,
}

/// Run IPF on every target PSM of `map` that has evidence in `evidence`
/// (keyed by PIN `SpecId`). Returns rows sorted by PSM id and the model.
/// `debug` fills [`IpfRow::forms_debug`].
pub fn localize_psms(
    map: &HashMap<String, PercolatorPsm>,
    evidence: &HashMap<String, IpfEvidence>,
    debug: bool,
) -> (Vec<IpfRow>, IonPepModel) {
    // Two passes: fit on the search's placements, localize, refit on the
    // placements IPF prefers (so the model does not just learn the search's
    // choice), as an EM step.
    let train: Vec<&IpfEvidence> = map
        .iter()
        .filter(|(_, p)| p.q_value <= TRAIN_Q)
        .filter_map(|(id, _)| evidence.get(id))
        .collect();
    let first = IonPepModel::fit(train.iter().copied());
    let model =
        IonPepModel::fit_labelled(train.iter().map(|&e| (e, localize(0.0, e, &first).best)));
    let mut ids: Vec<&String> = map.keys().filter(|id| evidence.contains_key(*id)).collect();
    ids.sort();
    let mut rows: Vec<IpfRow> = ids
        .into_iter()
        .map(|id| {
            let psm = &map[id];
            let e = &evidence[id];
            let loc = localize(psm.pep, e, &model);
            let total: f64 = loc.forms.iter().sum();
            let mut site_prob: std::collections::BTreeMap<(u16, &str), f64> = Default::default();
            for (k, sites) in e.forms.iter().enumerate() {
                let w = if total > 0.0 { loc.forms[k] / total } else { 0.0 };
                for (&s, m) in sites.iter().zip(&e.mods[k]) {
                    *site_prob.entry((s + 1, m.name.as_str())).or_default() += w;
                }
            }
            let one_based = |k: usize| {
                e.forms[k]
                    .iter()
                    .zip(&e.mods[k])
                    .map(|(s, m)| (s + 1, m.name.clone()))
                    .collect::<Vec<_>>()
            };
            IpfRow {
                psm_id: id.clone(),
                peptide: psm.peptide.clone(),
                psm_q: psm.q_value,
                psm_pep: psm.pep,
                reported_sites: one_based(e.observed as usize),
                best_sites: one_based(loc.best as usize),
                best_prob: loc.best_prob,
                pep: loc.pep,
                q: 0.0,
                n_forms: e.forms.len(),
                site_ions: e.ions.len(),
                matched_ions: e.ions.iter().filter(|i| i.tier.is_some()).count(),
                truncated: e.truncated,
                site_probs: site_prob
                    .into_iter()
                    .map(|((s, m), p)| (s, m.to_string(), p))
                    .collect(),
                forms_debug: if debug {
                    // Evidence-only posterior (no PSM prior): its log is the
                    // IPF log-likelihood up to a constant per PSM.
                    let ll = search::ipf::evidence_posterior(e, &model);
                    e.forms
                        .iter()
                        .enumerate()
                        .map(|(k, s)| {
                            format!(
                                "{}|{:.3}|{:.4}",
                                s.iter()
                                    .map(|x| (x + 1).to_string())
                                    .collect::<Vec<_>>()
                                    .join(","),
                                e.form_scores.get(k).copied().unwrap_or(f32::NAN),
                                ll.forms[k].max(1e-300).ln()
                            )
                        })
                        .collect::<Vec<_>>()
                        .join(";")
                } else {
                    String::new()
                },
            }
        })
        .collect();
    let q = model_fdr(&rows.iter().map(|r| r.pep).collect::<Vec<_>>());
    for (r, q) in rows.iter_mut().zip(q) {
        r.q = q;
    }
    (rows, model)
}

/// What `psms.parquet` carries for one localized row.
#[derive(Debug, Clone, PartialEq)]
pub struct IpfQpx {
    /// Localization probability per (1-based residue, mod name), every
    /// residue some placement modifies (the `positions.scores` of the
    /// reported mods).
    pub site_probs: HashMap<(i32, String), f64>,
    pub best_prob: f64,
    pub pep: f64,
    pub q: f64,
    /// The most probable placement, as [`sites_label`] writes it.
    pub best_sites: String,
    /// The per-residue probabilities, as the TSV `site_probs` column.
    pub site_probs_label: String,
    pub truncated: bool,
}

/// The QPX view of every localized row, keyed by PIN `SpecId`.
pub fn qpx_map(rows: &[IpfRow]) -> HashMap<String, IpfQpx> {
    rows.iter()
        .map(|r| {
            let qpx = IpfQpx {
                site_probs: r
                    .site_probs
                    .iter()
                    .map(|(s, m, p)| ((*s as i32, m.clone()), *p))
                    .collect(),
                best_prob: r.best_prob,
                pep: r.pep,
                q: r.q,
                best_sites: sites_label(&r.best_sites, r),
                site_probs_label: site_probs_label(r),
                truncated: r.truncated,
            };
            (r.psm_id.clone(), qpx)
        })
        .collect()
}

/// Whether the row's placements carry more than one mod name; the labels
/// name the mod only then.
fn mixed_mods(r: &IpfRow) -> bool {
    let first = r.site_probs.first().map(|x| &x.1);
    r.site_probs.iter().any(|x| Some(&x.1) != first)
}

/// `3;7`, or `3[Phospho];7[Oxidation]` when the row mixes mods.
pub fn sites_label(sites: &[(u16, String)], r: &IpfRow) -> String {
    let named = mixed_mods(r);
    sites
        .iter()
        .map(|(s, m)| if named { format!("{s}[{m}]") } else { s.to_string() })
        .collect::<Vec<_>>()
        .join(";")
}

/// Residues at probability 0.01 or more, `3:0.981;5:0.019` (named as in
/// [`sites_label`] when the row mixes mods).
pub fn site_probs_label(r: &IpfRow) -> String {
    let named = mixed_mods(r);
    r.site_probs
        .iter()
        .filter(|x| x.2 >= 0.01)
        .map(|(s, m, p)| {
            if named {
                format!("{s}[{m}]:{p:.3}")
            } else {
                format!("{s}:{p:.3}")
            }
        })
        .collect::<Vec<_>>()
        .join(";")
}

/// Write `<stem>.ipf.tsv`.
pub fn write_ipf_tsv(path: &std::path::Path, rows: &[IpfRow]) -> std::io::Result<()> {
    let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
    let debug = rows.iter().any(|r| !r.forms_debug.is_empty());
    writeln!(
        w,
        "psm_id\tpeptide\tpsm_q\tpsm_pep\treported_sites\tbest_sites\tbest_prob\tipf_pep\tipf_q\tn_forms\tsite_ions\tmatched_ions\ttruncated\tsite_probs{}",
        if debug { "\tforms_debug" } else { "" }
    )?;
    for r in rows {
        write!(
            w,
            "{}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}",
            r.psm_id,
            r.peptide,
            r.psm_q,
            r.psm_pep,
            sites_label(&r.reported_sites, r),
            sites_label(&r.best_sites, r),
            r.best_prob,
            r.pep,
            r.q,
            r.n_forms,
            r.site_ions,
            r.matched_ions,
            r.truncated as u8,
            site_probs_label(r)
        )?;
        if debug {
            write!(w, "\t{}", r.forms_debug)?;
        }
        writeln!(w)?;
    }
    w.flush()
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::modification::{ModLocation, Modification, ResidueSpec};
    use search::ipf::{FormSet, IonMatch};

    fn mod_(name: &str, mass: f64) -> &'static Modification {
        Box::leak(Box::new(Modification {
            name: name.into(),
            mass_delta: mass,
            residue: ResidueSpec::Specific(b'K'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: vec![],
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            diagnostic_ions: vec![],
            parts: Vec::new(),
        }))
    }

    /// Two placements on the same residues that swap which mod sits where are
    /// told apart: probabilities and sites carry the mod, and the best
    /// placement can differ from the reported one on identity alone.
    #[test]
    fn placements_that_swap_mods_on_the_same_residues_stay_apart() {
        let (me, di) = (mod_("Methyl", 14.01565), mod_("Dimethyl", 28.0313));
        // Reported: Methyl on residue 2, Dimethyl on 5. Strong ions favour the swap.
        let e = IpfEvidence {
            forms: vec![vec![1, 4], vec![1, 4]],
            mods: vec![vec![me, di], vec![di, me]],
            observed: 0,
            truncated: false,
            ions: (0..6)
                .map(|_| IonMatch {
                    forms: FormSet::from_indices(2, &[1]),
                    tier: Some(0),
                })
                .chain((0..6).map(|_| IonMatch {
                    forms: FormSet::from_indices(2, &[0]),
                    tier: None,
                }))
                .collect(),
            form_scores: vec![],
        };
        // Several confident PSMs so the ion model learns that tier 0 is evidence.
        let mut map = HashMap::new();
        let mut evidence = HashMap::new();
        for i in 0..20 {
            let id = format!("r_{i}");
            map.insert(
                id.clone(),
                PercolatorPsm {
                    psm_id: id.clone(),
                    q_value: 0.001,
                    pooled_q_value: None,
                    pep: 0.001,
                    score: None,
                    peptide: "K.AKAAKA.A".into(),
                    proteins: "P".into(),
                },
            );
            let mut ei = e.clone();
            if i > 0 {
                // The training PSMs report the placement their ions support.
                ei.observed = 1;
            }
            evidence.insert(id, ei);
        }
        let (rows, _) = localize_psms(&map, &evidence, false);
        let r = rows.iter().find(|r| r.psm_id == "r_0").unwrap();
        assert_eq!(r.reported_sites, vec![(2, "Methyl".into()), (5, "Dimethyl".into())]);
        assert_eq!(r.best_sites, vec![(2, "Dimethyl".into()), (5, "Methyl".into())]);
        let prob = |s: u16, m: &str| {
            r.site_probs.iter().find(|x| x.0 == s && x.1 == m).map(|x| x.2).unwrap()
        };
        assert!(prob(2, "Dimethyl") > 0.9 && prob(2, "Methyl") < 0.1, "{:?}", r.site_probs);
        assert!((prob(2, "Dimethyl") + prob(2, "Methyl") - 1.0).abs() < 1e-9);
        assert_eq!(sites_label(&r.best_sites, r), "2[Dimethyl];5[Methyl]");
        let q = qpx_map(&rows);
        assert!(q["r_0"].site_probs[&(2, "Dimethyl".to_string())] > 0.9);
    }
}
