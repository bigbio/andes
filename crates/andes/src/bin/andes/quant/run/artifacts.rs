//! Run artifacts on disk: the reduced MS2 table, the run's feature tables and quant weights.

use std::path::Path;

use model::mass::PROTON;
use model::spectrum::Spectrum;
use model::Tolerance;
use output::features::FeatureRecord;
use output::ms2::{write_ms2_table, Ms2Table};
use output::quantities::runs::{write_run_quant_table, RunPsmRow};
use quant::ms2::Ms2Scan;
use search::candidates::Candidate;
use search::index::SearchIndex;
use search::psm::TopNQueue;

use super::collect::{QuantCollector, RunCapture};
use super::hits::{peptidoform_of, proteins_of, select_rank1_all, Hit};
use super::lfq::envelope_of;
use super::scans::native_id;
use super::MS2_TOP_PEAKS;

/// Roll-up weight of a reporter row that passes the gate: its purity squared (1 when unknown)
/// times the probability that its identification is right, 1 - PEP (1 when no PEP is known).
pub(crate) fn quant_weight(purity: Option<f32>, pep: Option<f64>) -> f32 {
    let right = 1.0 - pep.unwrap_or(0.0).clamp(0.0, 1.0);
    purity.unwrap_or(1.0).clamp(0.0, 1.0).powi(2) * right as f32
}

pub(super) fn reduced_ms2(spec: &Spectrum) -> Ms2Scan {
    let scan = spec.scan.unwrap_or(0);
    let base = if spec.title.is_empty() {
        format!("scan={scan}")
    } else {
        spec.title.clone()
    };
    let mut peaks: Vec<(f32, f32)> = spec
        .peaks
        .iter()
        .filter(|p| p.0.is_finite() && p.1.is_finite() && p.1 > 0.0)
        .map(|&(mz, i)| (mz as f32, i))
        .collect();
    if peaks.len() > MS2_TOP_PEAKS {
        peaks.select_nth_unstable_by(MS2_TOP_PEAKS - 1, |a, b| {
            b.1.total_cmp(&a.1).then(a.0.total_cmp(&b.0))
        });
        peaks.truncate(MS2_TOP_PEAKS);
    }
    peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
    let isolation = match (
        spec.isolation_target_mz.or(Some(spec.precursor_mz)),
        spec.isolation_lower_offset,
        spec.isolation_upper_offset,
    ) {
        (Some(t), Some(lo), Some(hi)) => Some((t, lo, hi)),
        _ => None,
    };
    Ms2Scan {
        native_id: spec.title.clone(),
        spec_key: format!("{base}_{scan}"),
        rt: spec.rt_seconds.unwrap_or(f64::NAN),
        precursor_mz: spec.precursor_mz,
        charge: spec
            .precursor_charge
            .and_then(|z| u8::try_from(z).ok())
            .unwrap_or(0),
        isolation,
        faims_cv: None,
        peaks,
        reporters: None,
    }
}

pub(super) fn write_run_ms2(
    dir: &Path,
    stem: &str,
    title_prefix: Option<&str>,
    mut scans: Vec<Ms2Scan>,
    fragment_tol: Option<Tolerance>,
    faims_cv: impl Fn(&str) -> Option<f32>,
) -> Result<(), String> {
    let fragment_tol = fragment_tol
        .ok_or_else(|| format!("--quant-dir: no fragment tolerance for the MS2 table of {stem}"))?;
    std::fs::create_dir_all(dir)
        .map_err(|e| format!("create --quant-dir {}: {e}", dir.display()))?;
    for s in &mut scans {
        if let Some(p) = title_prefix {
            if let Some(id) = s.native_id.strip_prefix(p) {
                s.native_id = id.to_string();
            }
        }
        s.faims_cv = faims_cv(&s.native_id);
    }
    scans.retain(|s| s.rt.is_finite());
    scans.sort_by(|a, b| a.spec_key.cmp(&b.spec_key));
    scans.dedup_by(|a, b| a.spec_key == b.spec_key);
    scans.sort_by(|a, b| a.rt.total_cmp(&b.rt).then(a.spec_key.cmp(&b.spec_key)));
    let path = dir.join(format!("{stem}.ms2.parquet"));
    let n = scans.len();
    write_ms2_table(
        &path,
        &Ms2Table {
            run: stem.to_string(),
            fragment_tol,
            top_peaks: MS2_TOP_PEAKS,
            scans,
        },
    )
    .map_err(|e| format!("write MS2 table: {e}"))?;
    let bytes = std::fs::metadata(&path).map(|m| m.len()).unwrap_or(0);
    log::info!(
        target: "quant",
        "{stem}: wrote {} ({n} MS2 scans, {:.1} MB)",
        path.display(),
        bytes as f64 / 1e6
    );
    Ok(())
}

const TABLE_ENVELOPE_ISOTOPES: usize = 4;

pub(super) fn write_run_artifacts(
    collector: &QuantCollector,
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &[Candidate],
    index: &SearchIndex,
    quant_dir: &Path,
) -> Result<(), String> {
    std::fs::create_dir_all(quant_dir)
        .map_err(|e| format!("create --quant-dir {}: {e}", quant_dir.display()))?;
    let hits = select_rank1_all(collector, spectra, queues, candidates);
    for (run_idx, run) in collector.runs.iter().enumerate() {
        let ms1 = run
            .ms1
            .as_ref()
            .ok_or_else(|| format!("--quant-dir: run {} has no MS1 scans", run.stem))?;
        let ms1_path = quant_dir.join(format!("{}.ms1.parquet", run.stem));
        let table_path = quant_dir.join(format!("{}.quant.parquet", run.stem));
        let run_hits: Vec<&Hit<'_>> = hits.iter().filter(|h| h.run_idx == run_idx).collect();
        let (persisted, written) = rayon::join(
            || ms1.persist(&ms1_path),
            || {
                let (rows, mods, no_rt) = run_psm_rows(run, &run_hits, spectra, candidates, index);
                write_run_quant_table(&table_path, &run.stem, &rows, &mods)
                    .map(|()| (rows.len(), no_rt))
                    .map_err(|e| format!("write {}: {e}", table_path.display()))
            },
        );
        persisted.map_err(|e| format!("write MS1 table: {e}"))?;
        let (n_rows, no_rt) = written?;
        log::info!(
            target: "quant",
            "{}: wrote {} and {} ({} PSMs, one per precursor of a spectrum, targets and decoys{})",
            run.stem,
            ms1_path.display(),
            table_path.display(),
            n_rows,
            if no_rt > 0 {
                format!("; {no_rt} without a retention time left out")
            } else {
                String::new()
            }
        );
    }
    Ok(())
}

/// The PSM rows of one run's quant table, built in parallel in hit order; the modifications
/// they carry (the first mass seen for each accession, in hit order); and how many hits had no
/// retention time and were left out.
fn run_psm_rows(
    run: &RunCapture,
    hits: &[&Hit<'_>],
    spectra: &[Spectrum],
    candidates: &[Candidate],
    index: &SearchIndex,
) -> (Vec<RunPsmRow>, Vec<(String, f64)>, usize) {
    use rayon::prelude::*;
    let mut mods: std::collections::BTreeMap<String, f64> = std::collections::BTreeMap::new();
    let mut no_rt = 0usize;
    for hit in hits {
        if spectra[hit.spec_idx].rt_seconds.is_none() {
            no_rt += 1;
            continue;
        }
        for m in hit.cand.peptide.residues.iter().filter_map(|aa| aa.mod_) {
            if let Some(acc) = m.accession.as_deref().filter(|a| !a.is_empty()) {
                mods.entry(acc.to_string()).or_insert(m.mass_delta);
            }
        }
    }
    let rows: Vec<RunPsmRow> = hits
        .par_iter()
        .filter_map(|hit| {
            let spec = &spectra[hit.spec_idx];
            let rt_s = spec.rt_seconds?;
            let z = hit.psm.charge_used.max(1);
            let peptidoform = peptidoform_of(hit.cand);
            let context = run.contexts.get(native_id(run, &spec.title));
            let (envelope, isotope_share) = envelope_of(hit.cand, TABLE_ENVELOPE_ISOTOPES);
            Some(RunPsmRow {
                spec_id: hit.spec_id.clone(),
                rt_s,
                precursor_key: format!("{peptidoform}/{z}"),
                peptidoform,
                charge: z,
                mz: hit.cand.peptide.mass() / z as f64 + PROTON,
                is_decoy: hit.cand.is_decoy,
                search_score: hit.psm.features.strong_score,
                faims_cv: context.and_then(|c| c.faims_cv),
                envelope,
                proteins: proteins_of(&hit.psm, candidates, index)
                    .into_iter()
                    .map(|(a, ..)| a)
                    .collect::<Vec<_>>()
                    .join(";"),
                isotope_share,
                ion_mobility: context.and_then(|c| c.ion_mobility),
            })
        })
        .collect();
    (rows, mods.into_iter().collect(), no_rt)
}

pub(super) fn write_feature_outputs(
    pdir: &Path,
    groups: Option<&crate::quant::ProteinGroups<'_>>,
    mut records: Vec<FeatureRecord>,
) -> Result<(), String> {
    std::fs::create_dir_all(pdir).map_err(|e| format!("create {}: {e}", pdir.display()))?;
    let mut tally = crate::quant::GroupTally::default();
    if let Some(pg) = groups {
        records.retain_mut(|r| match tally.assign(pg, &r.sequence) {
            Some((proteins, group)) => {
                r.proteins = proteins;
                r.group = Some(group);
                true
            }
            None => false,
        });
    }
    for r in &records {
        tally.quantify(r, &r.run_file_name)?;
    }
    let path = pdir.join("quantms.feature.parquet");
    output::write_feature_parquet(&path, &records)
        .map_err(|e| format!("write {}: {e}", path.display()))?;
    log::info!(
        target: "quant",
        "wrote {} ({} feature rows)",
        path.display(),
        records.len()
    );
    if let Some(pg) = groups {
        let mut runs: Vec<String> = records.iter().map(|r| r.run_file_name.clone()).collect();
        runs.sort();
        runs.dedup();
        tally.write_pg(pdir, pg, &runs)?;
    }
    Ok(())
}
