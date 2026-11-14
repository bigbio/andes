//! IPF site localization after rescoring (`--localize ipf`).
//!
//! Each rescored target PSM whose row carries [`search::ipf::IpfEvidence`] is
//! localized with the PSM's PEP as the prior (see [`search::ipf`]). The ion
//! PEP model is fit once per run on the evidence of confident target PSMs.
//! Peptidoform q-values are model-based over every localized PSM, as in
//! PyProphet IPF.

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
    /// 1-based residue positions of the reported (search) placement.
    pub reported_sites: Vec<u16>,
    /// 1-based residue positions of the most probable placement.
    pub best_sites: Vec<u16>,
    /// Localization probability of the best placement given the PSM is right.
    pub best_prob: f64,
    /// Peptidoform PEP (1 - posterior of the best placement, h0 included).
    pub pep: f64,
    pub q: f64,
    pub n_forms: usize,
    pub site_ions: usize,
    pub matched_ions: usize,
    pub truncated: bool,
    /// Per-residue localization probability (1-based residue, prob), for
    /// residues with prob >= 0.01, in residue order.
    pub site_probs: Vec<(u16, f64)>,
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
            let mut site_prob: std::collections::BTreeMap<u16, f64> = Default::default();
            if total > 0.0 {
                for (sites, w) in e.forms.iter().zip(&loc.forms) {
                    for &s in sites {
                        *site_prob.entry(s + 1).or_default() += w / total;
                    }
                }
            }
            let one_based = |v: &[u16]| v.iter().map(|s| s + 1).collect::<Vec<_>>();
            IpfRow {
                psm_id: id.clone(),
                peptide: psm.peptide.clone(),
                psm_q: psm.q_value,
                psm_pep: psm.pep,
                reported_sites: one_based(&e.forms[e.observed as usize]),
                best_sites: one_based(&e.forms[loc.best as usize]),
                best_prob: loc.best_prob,
                pep: loc.pep,
                q: 0.0,
                n_forms: e.forms.len(),
                site_ions: e.ions.len(),
                matched_ions: e.ions.iter().filter(|i| i.tier.is_some()).count(),
                truncated: e.truncated,
                site_probs: site_prob.into_iter().filter(|x| x.1 >= 0.01).collect(),
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

/// Per-row residue localization probabilities (1-based residue -> prob),
/// keyed by PIN `SpecId`, for the QPX `positions.scores` column.
pub fn site_probability_map(rows: &[IpfRow]) -> HashMap<String, HashMap<i32, f64>> {
    rows.iter()
        .map(|r| {
            let m = r.site_probs.iter().map(|&(s, p)| (s as i32, p)).collect();
            (r.psm_id.clone(), m)
        })
        .collect()
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
    let sites = |v: &[u16]| v.iter().map(u16::to_string).collect::<Vec<_>>().join(";");
    for r in rows {
        write!(
            w,
            "{}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}",
            r.psm_id,
            r.peptide,
            r.psm_q,
            r.psm_pep,
            sites(&r.reported_sites),
            sites(&r.best_sites),
            r.best_prob,
            r.pep,
            r.q,
            r.n_forms,
            r.site_ions,
            r.matched_ions,
            r.truncated as u8,
            r.site_probs
                .iter()
                .map(|(s, p)| format!("{s}:{p:.3}"))
                .collect::<Vec<_>>()
                .join(";")
        )?;
        if debug {
            write!(w, "\t{}", r.forms_debug)?;
        }
        writeln!(w)?;
    }
    w.flush()
}
