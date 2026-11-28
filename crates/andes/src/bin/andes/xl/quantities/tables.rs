//! The cross-link quantification tables: the feature table (one row per precursor and run)
//! and the residue-pair matrix over the runs.

use std::collections::BTreeMap;
use std::io::Write;
use std::path::{Path, PathBuf};

use output::crosslinks::features::{write_xl_features, XlFeatureRow};

use super::library::Library;
use super::{accepted, Row, Settings};

pub(super) struct Written {
    pub(super) features: PathBuf,
    pub(super) residue_pairs: usize,
}

fn sibling(report_base: &Path, suffix: &str) -> PathBuf {
    let mut s = if report_base
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("pin"))
    {
        report_base.with_extension("").into_os_string()
    } else {
        report_base.to_path_buf().into_os_string()
    };
    s.push(suffix);
    PathBuf::from(s)
}

fn feature_row(row: &Row, lib: &Library, run: &str, settings: &Settings) -> XlFeatureRow {
    let p = &lib.precursors[row.precursor];
    let f = row.feature.as_ref();
    let [a, b] = &p.sides;
    XlFeatureRow {
        run: run.to_string(),
        precursor: p.key.clone(),
        residue_pair: p.residue_pair.clone(),
        protein_pair: p.protein_pair.clone(),
        link: p.link.to_string(),
        crosslinker: p.linker.clone(),
        alpha_peptide: a.bare.clone(),
        alpha_peptidoform: a.peptide.clone(),
        alpha_site: u32::try_from(a.site).unwrap_or(u32::MAX),
        alpha_site_kind: a.site_kind.clone(),
        beta_peptide: b.bare.clone(),
        beta_peptidoform: b.peptide.clone(),
        beta_site: u32::try_from(b.site).unwrap_or(u32::MAX),
        beta_site_kind: b.site_kind.clone(),
        charge: u32::from(p.charge),
        mono_mz: p.mono_mz,
        evidence: row.evidence.as_str().to_string(),
        accepted: accepted(row, settings),
        quantity: f.map(|f| f.area),
        quant_q: row.q,
        feature_score: row.score.map(f64::from),
        anchor_rt_s: row.anchor_rt,
        anchor_sd_s: row.anchor_sd_s,
        apex_rt_s: f.map(|f| f.apex_rt),
        rt_start_s: f.map(|f| f.rt_start),
        rt_stop_s: f.map(|f| f.rt_stop),
        apex_intensity: f.map(|f| f64::from(f.apex_intensity)),
        cosine: f.map(|f| f64::from(f.cosine)),
        n_isotopes: u32::try_from(p.envelope.len()).unwrap_or(u32::MAX),
        n_isotopes_found: f.map_or(0, |f| f.n_isotopes_found),
        scan: row.scan,
        csm_score: row.csm_score.is_finite().then_some(row.csm_score),
        n_csms: row.n_csms,
        residue_pair_q: p.q,
    }
}

struct ResiduePair<'a> {
    protein_pair: &'a str,
    link: &'a str,
    q: f64,
    precursors: std::collections::BTreeSet<usize>,
    per_run: Vec<f64>,
}

/// The residue-pair matrix: for each residue pair and run, the sum of the accepted quantities
/// of its precursors; an empty cell where none is accepted.
fn write_residue_pairs(
    path: &Path,
    lib: &Library,
    rows: &[Row],
    runs: &[&str],
    settings: &Settings,
) -> Result<usize, String> {
    let mut pairs: BTreeMap<&str, ResiduePair<'_>> = BTreeMap::new();
    for row in rows {
        let p = &lib.precursors[row.precursor];
        let entry = pairs
            .entry(p.residue_pair.as_str())
            .or_insert_with(|| ResiduePair {
                protein_pair: &p.protein_pair,
                link: p.link,
                q: p.q,
                precursors: std::collections::BTreeSet::new(),
                per_run: vec![0.0; runs.len()],
            });
        if let (true, Some(f)) = (accepted(row, settings), &row.feature) {
            entry.precursors.insert(row.precursor);
            entry.per_run[row.run] += f.area;
        }
    }
    let err = |e: std::io::Error| format!("{}: {e}", path.display());
    let mut out = std::io::BufWriter::new(std::fs::File::create(path).map_err(err)?);
    write!(
        out,
        "residue_pair\tprotein_pair\tlink\tq_value\tn_precursors"
    )
    .map_err(err)?;
    for run in runs {
        write!(out, "\t{run}").map_err(err)?;
    }
    writeln!(out).map_err(err)?;
    for (key, pair) in &pairs {
        write!(
            out,
            "{key}\t{}\t{}\t{:.6}\t{}",
            pair.protein_pair,
            pair.link,
            pair.q,
            pair.precursors.len()
        )
        .map_err(err)?;
        for v in &pair.per_run {
            if *v > 0.0 {
                write!(out, "\t{v:.1}").map_err(err)?;
            } else {
                write!(out, "\t").map_err(err)?;
            }
        }
        writeln!(out).map_err(err)?;
    }
    out.flush().map_err(err)?;
    Ok(pairs.len())
}

pub(super) fn write(
    report_base: &Path,
    lib: &Library,
    rows: &[Row],
    runs: &[&str],
    settings: &Settings,
) -> Result<Written, String> {
    let features = sibling(report_base, ".xl.features.parquet");
    let table: Vec<XlFeatureRow> = rows
        .iter()
        .map(|row| feature_row(row, lib, runs[row.run], settings))
        .collect();
    write_xl_features(&features, &table)?;
    let residue_pairs = write_residue_pairs(
        &sibling(report_base, ".xl.residue_pairs.quant.tsv"),
        lib,
        rows,
        runs,
        settings,
    )?;
    Ok(Written {
        features,
        residue_pairs,
    })
}

#[cfg(test)]
#[path = "tests/tables.rs"]
mod tests;
