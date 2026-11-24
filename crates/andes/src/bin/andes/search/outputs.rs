//! Side outputs of a search: the claim database, open-search offsets, filtered TSVs and the
//! crosslink precursor snap.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::path::{Path, PathBuf};

pub(crate) struct ClaimDatabase {
    pub(crate) loaded: input::LoadedDatabase,
    pub(crate) annotations: Option<input::AnnotationLoadStats>,
    pub(crate) site_knockoffs: input::SiteKnockoffStats,
}

pub(crate) fn load_claim_database(
    path: &Path,
    annotations: &[PathBuf],
    opts: &input::UniProtOptions,
) -> Result<ClaimDatabase, Box<dyn std::error::Error>> {
    let mut loaded = input::load_protein_database(path, opts)?;
    let annotations = if annotations.is_empty() {
        None
    } else {
        Some(input::load_annotations(
            annotations,
            &mut loaded.db,
            opts.knockoffs,
        )?)
    };
    let site_knockoffs = opts
        .knockoffs
        .map(|ko| input::add_site_knockoffs(&mut loaded.db, ko))
        .unwrap_or_default();
    Ok(ClaimDatabase {
        loaded,
        annotations,
        site_knockoffs,
    })
}

pub(super) fn release_at_exit(candidates: Vec<search::candidates::Candidate>) {
    std::mem::forget(candidates);
}

pub(super) fn write_open_offsets(
    pin: &Path,
    dir: &Path,
    stem: &str,
    accepted: &[&output::PercolatorPsm],
    decoys: &std::collections::HashMap<String, output::PercolatorPsm>,
) -> Option<String> {
    let report =
        match output::offsets::report_from_results(pin, accepted, decoys, model::unimod::bundled())
        {
            Ok(r) => r,
            Err(e) => {
                log::warn!("open-search offset report not written: {e}");
                return None;
            }
        };
    let peaks_path = dir.join(format!("{stem}.open_offsets.parquet"));
    let psms_path = dir.join(format!("{stem}.open_psms.parquet"));
    let written = output::offsets::write_peaks_table(&peaks_path, &report)
        .and_then(|()| output::offsets::write_psms_table(&psms_path, &report));
    if let Err(e) = written {
        let mut left = Vec::new();
        for stale in [&peaks_path, &psms_path] {
            if let Err(r) = std::fs::remove_file(stale) {
                if r.kind() != std::io::ErrorKind::NotFound {
                    log::warn!(target: "output", "could not remove {}: {r}", stale.display());
                    left.push(stale.display().to_string());
                }
            }
        }
        if left.is_empty() {
            log::warn!(
                "open-search report not written ({e}); neither {} nor {} is left",
                peaks_path.display(),
                psms_path.display()
            );
        } else {
            log::warn!(
                "open-search report not written ({e}); incomplete file left: {}",
                left.join(", ")
            );
        }
        return None;
    }
    log::info!(
        target: "output",
        "Wrote open-search report: {} ({} offset peaks) and {} ({} PSMs)",
        peaks_path.display(),
        report.peaks.len(),
        psms_path.display(),
        report.assignments.len()
    );
    Some(output::offsets::summary_block(&report))
}

pub(crate) fn write_filtered_tsv(
    path: &Path,
    psms: &[&output::PercolatorPsm],
    groups: Option<&output::Inference>,
) -> std::io::Result<()> {
    use std::io::Write;
    output::write_atomic(path, |path| {
        let mut w = std::io::BufWriter::new(File::create(path)?);
        writeln!(
            w,
            "PSMId\tq-value\tposterior_error_probability\tpeptide\tproteinIds\tprotein_group\tprotein_q"
        )?;
        for p in psms {
            let (group, q) = groups
                .and_then(|g| g.group_of(p))
                .map_or(("NA".to_string(), "NA".to_string()), |g| {
                    (g.id.to_string(), g.q_value.to_string())
                });
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                p.psm_id, p.q_value, p.pep, p.peptide, p.proteins, group, q
            )?;
        }
        w.flush()
    })
}

/// Whether a crosslink search moves this file's precursor m/z values to the dominant MS1 peak
/// (`search::crosslinks::PRECURSOR_SNAP_PPM`): an MS2 search of a file that carries its MS1
/// scans (mzML, Thermo `.raw`). Decided per file, so the order of the inputs does not matter.
pub(crate) fn xl_precursor_snap(
    xl: bool,
    ms_level: u8,
    file_is_mzml: bool,
    file_is_raw: bool,
) -> bool {
    xl && ms_level == 2 && (file_is_mzml || file_is_raw)
}

/// The search's PSM tables as the output writers read them.
#[derive(Clone, Copy)]
pub(super) struct PsmTables<'a> {
    pub(super) spectra: &'a [model::Spectrum],
    pub(super) queues: &'a [search::TopNQueue],
    pub(super) candidates: &'a [search::candidates::Candidate],
    pub(super) params: &'a search::SearchParams,
    pub(super) index: &'a search::SearchIndex,
}

/// Rescoring results the QPX bundle reports: target PSMs, decoy ids, rows left out, IPF sites.
pub(super) type QpxRescore<'a> = (
    Option<&'a HashMap<String, output::PercolatorPsm>>,
    Option<&'a HashSet<String>>,
    &'a HashSet<String>,
    Option<&'a HashMap<String, output::ipf::IpfQpx>>,
);

/// Infers the protein groups of the accepted target PSMs (q under their class cutoff, not a
/// knockoff row) and writes them next to the report.
pub(super) fn infer_and_write_groups(
    map: &HashMap<String, output::PercolatorPsm>,
    rescore_decoy_psms: &HashMap<String, output::PercolatorPsm>,
    pin_index: &search::SearchIndex,
    knockoff_rows: &HashSet<String>,
    report_base: &Path,
    (fdr_threshold, proteoform_fdr): (f64, f64),
) -> output::Inference {
    let cutoffs = output::Cutoffs::new(map, pin_index, fdr_threshold, proteoform_fdr);
    let inf = output::infer_protein_groups(map, rescore_decoy_psms, pin_index, |id| {
        map[id].q_value <= cutoffs.of(id) && !knockoff_rows.contains(id)
    });
    let stem = report_base
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "andes".to_string());
    let dir = report_base.parent().unwrap_or_else(|| Path::new("."));
    let path = dir.join(format!("{stem}.protein_groups.tsv"));
    let at = |q: f64| inf.groups.iter().filter(|g| g.q_value <= q).count();
    match output::write_protein_groups_tsv(&path, &inf) {
        Ok(()) => log::info!(
            target: "output",
            "Wrote protein groups: {} ({} groups; {} at protein q<=0.01, {} at q<=0.05)",
            path.display(),
            inf.groups.len(),
            at(0.01),
            at(0.05)
        ),
        Err(e) => log::warn!(target: "output", "could not write {}: {e}", path.display()),
    }
    inf
}

/// Writes the QPX bundle (`--output-parquet`).
pub(super) fn write_qpx_bundle(
    parquet_dir: &Path,
    t: &PsmTables<'_>,
    fragment_tol: &model::tolerance::Tolerance,
    (spectrum_path, original_spectrum_paths, run_ends): (&Path, &[PathBuf], &[usize]),
    (rescore_map, rescore_decoys, knockoff_rows, ipf_qpx): QpxRescore<'_>,
) -> Result<(), Box<dyn std::error::Error>> {
    let run_id = spectrum_path
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "andes_run".to_string());
    let primary_paths: Vec<String> = original_spectrum_paths
        .iter()
        .map(|p| p.display().to_string())
        .collect();
    match output::write_qpx(
        parquet_dir,
        t.spectra,
        t.queues,
        t.candidates,
        t.params,
        t.index,
        fragment_tol,
        &run_id,
        &primary_paths,
        run_ends,
        rescore_map,
        rescore_decoys,
        knockoff_rows,
        ipf_qpx,
    ) {
        Ok(()) => log::info!(target: "output", "Wrote QPX bundle: {}", parquet_dir.display()),
        Err(e) => {
            return Err(format!(
                "failed to write requested --output-parquet {}: {e}",
                parquet_dir.display()
            )
            .into());
        }
    }
    Ok(())
}

/// Writes the proteoform table of the accepted PSMs (families under `--protein-fdr` when set).
pub(super) fn write_proteoform_table(
    map: &HashMap<String, output::PercolatorPsm>,
    pin_index: &search::SearchIndex,
    cutoffs: &output::Cutoffs,
    quant_summary: Option<&output::QuantSummary>,
    protein_fdr_groups: Option<(f64, &output::Inference)>,
    report_base: &Path,
) {
    let mut rows = output::proteoform_table(map, pin_index, cutoffs, quant_summary);
    if let Some((t, inf)) = protein_fdr_groups {
        let fq = inf.family_q();
        rows.retain(|r| fq.get(r.family.as_str()).is_some_and(|&q| q <= t));
    }
    let stem = report_base
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "andes".to_string());
    let dir = report_base.parent().unwrap_or_else(|| Path::new("."));
    let path = dir.join(format!("{stem}.proteoforms.tsv"));
    let columns = quant_summary.map(|q| q.columns.clone()).unwrap_or_default();
    match output::write_proteoform_tsv(&path, &rows, &columns) {
        Ok(()) => log::info!(
            target: "output",
            "Wrote proteoform table: {} ({} family members, {} families)",
            path.display(),
            rows.len(),
            rows.iter()
                .map(|r| r.family.as_str())
                .collect::<std::collections::BTreeSet<_>>()
                .len()
        ),
        Err(e) => log::warn!(target: "output", "could not write {}: {e}", path.display()),
    }
}

/// Writes the search TSV (`--output-tsv`), naming the spectrum file(s) it came from.
pub(super) fn write_search_tsv(
    tsv_path: &Path,
    t: &PsmTables<'_>,
    spectrum_paths: &[PathBuf],
    is_mgf: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let spec_file_name = if spectrum_paths.len() == 1 {
        spectrum_paths[0]
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| spectrum_paths[0].display().to_string())
    } else {
        spectrum_paths
            .iter()
            .filter_map(|p| p.file_name().map(|n| n.to_string_lossy().into_owned()))
            .collect::<Vec<_>>()
            .join("+")
    };
    output::write_tsv(
        tsv_path,
        t.spectra,
        t.queues,
        t.candidates,
        t.params,
        t.index,
        &spec_file_name,
        is_mgf,
    )?;
    log::info!(target: "output", "Wrote TSV: {}", tsv_path.display());
    Ok(())
}
