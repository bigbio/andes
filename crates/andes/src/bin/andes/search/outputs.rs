//! Side outputs of a search: the claim database, open-search offsets, filtered TSVs and the
//! crosslink precursor snap.

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
