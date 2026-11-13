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
    /// Per-residue localization probability, `pos:prob` for prob >= 0.01.
    pub site_probs: String,
}

/// Run IPF on every target PSM of `map` that has evidence in `evidence`
/// (keyed by PIN `SpecId`). Returns rows sorted by PSM id and the model.
pub fn localize_psms(
    map: &HashMap<String, PercolatorPsm>,
    evidence: &HashMap<String, IpfEvidence>,
) -> (Vec<IpfRow>, IonPepModel) {
    let model = IonPepModel::fit(
        map.iter()
            .filter(|(_, p)| p.q_value <= TRAIN_Q)
            .filter_map(|(id, _)| evidence.get(id)),
    );
    let mut ids: Vec<&String> = map.keys().filter(|id| evidence.contains_key(*id)).collect();
    ids.sort();
    let mut rows: Vec<IpfRow> = ids
        .into_iter()
        .map(|id| {
            let psm = &map[id];
            let e = &evidence[id];
            let loc = localize(psm.pep, e, &model);
            let total: f64 = loc.forms.iter().sum();
            let mut site_prob: HashMap<u16, f64> = HashMap::new();
            if total > 0.0 {
                for (sites, w) in e.forms.iter().zip(&loc.forms) {
                    for &s in sites {
                        *site_prob.entry(s).or_default() += w / total;
                    }
                }
            }
            let mut sp: Vec<(u16, f64)> = site_prob.into_iter().filter(|x| x.1 >= 0.01).collect();
            sp.sort_by_key(|x| x.0);
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
                site_probs: sp
                    .iter()
                    .map(|(s, p)| format!("{}:{p:.3}", s + 1))
                    .collect::<Vec<_>>()
                    .join(";"),
            }
        })
        .collect();
    let q = model_fdr(&rows.iter().map(|r| r.pep).collect::<Vec<_>>());
    for (r, q) in rows.iter_mut().zip(q) {
        r.q = q;
    }
    (rows, model)
}

/// Write `<stem>.ipf.tsv`.
pub fn write_ipf_tsv(path: &std::path::Path, rows: &[IpfRow]) -> std::io::Result<()> {
    let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        w,
        "psm_id\tpeptide\tpsm_q\tpsm_pep\treported_sites\tbest_sites\tbest_prob\tipf_pep\tipf_q\tn_forms\tsite_ions\tmatched_ions\ttruncated\tsite_probs"
    )?;
    let sites = |v: &[u16]| v.iter().map(u16::to_string).collect::<Vec<_>>().join(";");
    for r in rows {
        writeln!(
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
        )?;
    }
    w.flush()
}
