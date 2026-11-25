//! Reading and searching one input file of a search (mzML, MGF, Thermo .raw, Bruker .d): the
//! parser threads, the chunked search, monoisotopic correction and MS1/MS3 capture.

use std::path::{Path, PathBuf};
use std::sync::mpsc::sync_channel;
use std::thread;

use input::{MgfReader, Ms1Link, MzMLReader};
use model::Spectrum;
use search::fragments::index::{index_chunk_len, index_order_mass, INDEX_CHUNK_SIZE};
use search::precursor::monoisotopic::{MonoCorrection, MonoParams};
use search::{PreparedSearch, SearchParams, TopNQueue};

use super::calibrate::{calibrate_dia_rt, recalibrate_dia_ms1, FragmentErrorMeasurer};
use super::outputs::xl_precursor_snap;
use super::{Ms1Peaks, CHUNK_SIZE};
use crate::cli::Cli;
use crate::memory::available_memory_budget;
use crate::monoisotopic::{correct_chunk, MonoDump, MonoStats};
use crate::quant::run::{QuantCollector, RunScans};
use crate::report_search_progress;
use crate::spectra::{
    add_mzml_skipped, input_format_flags, merge_parse_stats, prefix_spectrum_titles,
    require_spectra, send_chunks, send_linked_chunks, title_prefix_for, ParseStats,
};

fn index_chunk_span_da(flag: Option<f64>) -> f64 {
    const GIB: f64 = 1024.0 * 1024.0 * 1024.0;
    match flag {
        Some(d) => d.max(1.0),
        None => {
            let budget = available_memory_budget()
                .map(|b| b.bytes as f64)
                .unwrap_or(32.0 * GIB);
            ((budget * 0.6 - 4.0 * GIB) / (0.13 * GIB)).clamp(10.0, 150.0)
        }
    }
}

fn is_gzip_path(path: &Path) -> bool {
    path.extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("gz"))
}

#[allow(clippy::too_many_arguments)]
fn score_mass_ordered(
    prepared: &PreparedSearch,
    mut pending: Vec<Spectrum>,
    span_da: f64,
    all_spectra: &mut Vec<Spectrum>,
    all_queues: &mut Vec<TopNQueue>,
    t_search_start: std::time::Instant,
    keep_peaks: bool,
    mut quant: Option<&mut QuantCollector>,
    measurer: &mut FragmentErrorMeasurer,
) {
    if prepared.params.dia {
        prepared.set_dia_rt_prior(None);
    }
    measurer.measure(prepared, pending.iter(), Some(span_da));
    if prepared.params.dia {
        calibrate_dia_rt(prepared, &pending, span_da);
    }
    pending.sort_by(|a, b| index_order_mass(a).total_cmp(&index_order_mass(b)));
    log::info!(
        target: "fragment-index",
        "{} spectra sorted by precursor mass, chunks of <= {} spectra and <= {:.0} Da",
        pending.len(),
        INDEX_CHUNK_SIZE,
        span_da
    );
    let mut rest = pending;
    while !rest.is_empty() {
        let take = index_chunk_len(&rest, index_order_mass, span_da);
        let chunk: Vec<Spectrum> = rest.drain(..take).collect();
        let offset = all_spectra.len();
        let queues = prepared.run_chunk(&chunk, offset);
        all_queues.extend(queues);
        for mut spec in chunk.into_iter() {
            if let Some(q) = quant.as_deref_mut() {
                q.on_spectrum(&spec);
            }
            if !keep_peaks {
                spec.peaks = Vec::new();
            }
            all_spectra.push(spec);
        }
        report_search_progress(all_spectra.len(), t_search_start);
    }
    if prepared.params.dia {
        prepared.set_dia_rt_prior(None);
    }
}

/// What every input file of a search shares: the flags, the resolved and per-file parameters,
/// the inputs, and the switches `run` derives before the search loop.
#[derive(Clone, Copy)]
pub(super) struct FileInputs<'a> {
    pub(super) cli: &'a Cli,
    pub(super) params: &'a SearchParams,
    pub(super) file_params: &'a [SearchParams],
    pub(super) spectrum_paths: &'a [PathBuf],
    pub(super) bench_cap: usize,
    pub(super) bench_mode: bool,
    pub(super) chimeric_active: bool,
    pub(super) mono_active: bool,
    pub(super) mono_params: &'a MonoParams,
    pub(super) ms1_input: bool,
    pub(super) ms_level_u32: u32,
    pub(super) keep_peaks: bool,
    pub(super) standalone: bool,
    pub(super) t_search_start: std::time::Instant,
    pub(super) xl_ms3_capture: bool,
}

/// What the search loop accumulates over the input files.
pub(super) struct FileOutputs {
    pub(super) all_spectra: Vec<Spectrum>,
    pub(super) all_queues: Vec<TopNQueue>,
    pub(super) run_ends: Vec<usize>,
    pub(super) dia_ms1: Vec<Vec<model::scan::Ms1Scan>>,
    pub(super) xl_ms3: Vec<Vec<model::scan::ProductScan>>,
    pub(super) xl_snapped: usize,
    pub(super) xl_snap_files: usize,
    pub(super) parse_stats: ParseStats,
    pub(super) mono_table: Vec<Option<MonoCorrection>>,
    pub(super) mono_stats: MonoStats,
    pub(super) mono_dump: Option<MonoDump>,
}

impl FileOutputs {
    /// Empty outputs for a search of `n_files` files, writing the monoisotopic dump to `mono_dump`.
    pub(super) fn new(n_files: usize, mono_dump: Option<MonoDump>) -> Self {
        FileOutputs {
            all_spectra: Vec::new(),
            all_queues: Vec::new(),
            run_ends: Vec::with_capacity(n_files),
            dia_ms1: Vec::new(),
            xl_ms3: Vec::new(),
            xl_snapped: 0,
            xl_snap_files: 0,
            parse_stats: ParseStats::default(),
            mono_table: Vec::new(),
            mono_stats: MonoStats::default(),
            mono_dump,
        }
    }
}

/// Reads, searches and records one input file of `run`: its spectra and queues are appended to
/// `outputs`, its quantification scans go to `quant`.
pub(super) fn search_file<'a>(
    file_inputs: &FileInputs<'a>,
    outputs: &mut FileOutputs,
    prepared: &mut PreparedSearch<'a>,
    measurer: &mut FragmentErrorMeasurer,
    quant: &mut Option<QuantCollector>,
    file_idx: usize,
    input_path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let FileInputs {
        cli,
        file_params,
        spectrum_paths,
        bench_cap,
        bench_mode,
        chimeric_active,
        mono_active,
        ms1_input,
        xl_ms3_capture,
        ..
    } = *file_inputs;
    if let Some(p) = file_params.get(file_idx) {
        prepared.params = p;
    }
    measurer.start_file(prepared);
    let remaining_cap = if bench_mode {
        bench_cap.saturating_sub(outputs.all_spectra.len())
    } else {
        usize::MAX
    };

    let (file_is_mzml, file_is_raw, file_is_d, _file_is_mgf) = input_format_flags(input_path);
    let file_stem = input_path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("spectrum");
    let title_prefix = title_prefix_for(file_stem);

    if spectrum_paths.len() > 1 {
        log::info!(
            target: "input",
            "spectrum file {}/{}: {}",
            file_idx + 1,
            spectrum_paths.len(),
            input_path.display()
        );
    }

    let file_ms1 = file_is_mzml || file_is_raw || file_is_d;
    let file_snap = xl_precursor_snap(cli.xl.is_some(), cli.ms_level, file_is_mzml, file_is_raw);
    outputs.xl_snap_files += usize::from(file_snap);
    // This file reads its MS1 scans with its MS2 scans.
    let file_linked = (ms1_input && file_ms1) || file_snap;
    if (chimeric_active || mono_active) && !file_ms1 {
        return Err(format!(
            "{} only supports mzML/.raw/.d inputs, got {}",
            if chimeric_active {
                "--chimeric"
            } else {
                "--precursor-mono"
            },
            input_path.display()
        )
        .into());
    }

    let spectra_before = outputs.all_spectra.len();
    let mut run_scans = RunScans::default();
    let locals = FileLocals {
        input_path,
        remaining_cap,
        file_is_mzml,
        file_is_raw,
        file_is_d,
        title_prefix: &title_prefix,
        file_snap,
    };
    let route = read_route(file_inputs, quant.as_ref(), &locals, file_linked);
    let file_stats = match route {
        ReadRoute::Linked => search_linked(
            file_inputs,
            outputs,
            prepared,
            measurer,
            quant,
            &locals,
            &mut run_scans,
        )?,
        ReadRoute::Indexed => {
            search_indexed(file_inputs, outputs, prepared, measurer, quant, &locals)?
        }
        ReadRoute::IndexedMs1 | ReadRoute::Enumerated => {
            let read_ms1 = route == ReadRoute::IndexedMs1;
            let indexed = if read_ms1 {
                search_indexed_ms1(
                    file_inputs,
                    outputs,
                    prepared,
                    measurer,
                    quant,
                    &locals,
                    &mut run_scans,
                )?
            } else {
                None
            };
            match indexed {
                Some(stats) => stats,
                None => search_enumerated(
                    file_inputs,
                    outputs,
                    prepared,
                    measurer,
                    quant,
                    &locals,
                    &mut run_scans,
                )?,
            }
        }
    };
    let FileOutputs {
        all_spectra,
        run_ends,
        dia_ms1,
        xl_ms3,
        parse_stats,
        ..
    } = outputs;
    require_spectra(
        input_path,
        if file_linked { 2 } else { cli.ms_level },
        all_spectra.len() - spectra_before,
    )
    .map_err(|e| crate::spectra::with_first_read_error(e, &file_stats.first_errors))?;
    merge_parse_stats(parse_stats, file_stats);
    if cli.dia {
        let gated = prepared.take_dia_ms1();
        dia_ms1.push(if !gated.is_empty() {
            gated
        } else if quant.is_some() {
            run_scans.ms1.clone()
        } else {
            std::mem::take(&mut run_scans.ms1)
        });
    }
    if xl_ms3_capture {
        let ms3 = std::mem::take(&mut run_scans.ms3);
        let n = ms3.len();
        let linked = crate::xl::link_ms3(
            &all_spectra[spectra_before..],
            spectra_before,
            title_prefix.as_deref(),
            ms3,
            xl_ms3,
        );
        if n > 0 {
            log::info!(
                target: "xl",
                "{} MS3 scans, {} linked to an MS2 scan of {}",
                n,
                linked,
                input_path.display()
            );
        }
    }
    if let Some(q) = quant.as_mut() {
        q.finish_file(
            input_path,
            title_prefix.as_deref(),
            run_scans,
            all_spectra.len(),
        )?;
    }
    run_ends.push(all_spectra.len());
    Ok(())
}

/// The per-file values the read-and-search helpers of `search_file` share.
#[derive(Clone, Copy)]
struct FileLocals<'l> {
    input_path: &'l Path,
    remaining_cap: usize,
    file_is_mzml: bool,
    file_is_raw: bool,
    file_is_d: bool,
    title_prefix: &'l Option<String>,
    file_snap: bool,
}

/// The parse error count, the first errors and the run scans a parser thread returns.
type Ms1Parsed = (usize, Vec<String>, RunScans);

/// What the MS1-linked parser thread of `search_linked` owns: the file, the spectrum cap, the
/// reader settings and what to capture.
struct LinkedRead {
    spectrum_path: PathBuf,
    cap: usize,
    mslevel: u32,
    capture_ms1: bool,
    capture_context: bool,
    capture_noise: bool,
    parallel_read: bool,
    threads: usize,
    capture_product: bool,
    file_is_mzml: bool,
    file_is_d: bool,
    #[cfg(feature = "timstof")]
    dia_mobility_slices: usize,
}

/// Parses an MS1-linked file (mzML in parallel or streamed, Bruker .d, Thermo .raw) and sends its
/// MS2 chunks, each with its MS1 link, to `tx`.
fn parse_linked(
    read: LinkedRead,
    tx: std::sync::mpsc::SyncSender<(Vec<Spectrum>, Ms1Link)>,
) -> Result<Ms1Parsed, String> {
    let LinkedRead {
        spectrum_path,
        cap,
        mslevel,
        capture_ms1,
        capture_context,
        capture_noise,
        parallel_read,
        threads,
        capture_product,
        file_is_mzml,
        file_is_d,
        #[cfg(feature = "timstof")]
        dia_mobility_slices,
    } = read;
    let linked = if parallel_read {
        input::mzml::read_ms1_linked_parallel(
            &spectrum_path,
            threads,
            capture_context,
            capture_noise,
            false,
        )
        .map_err(|e| format!("read mzML: {e}"))?
    } else {
        None
    };
    if let Some(read) = linked {
        send_linked_chunks(
            read.spectra,
            &read.ms2_to_ms1,
            &read.run_ms1,
            CHUNK_SIZE,
            &tx,
        );
        let scans = RunScans {
            ms1: read.run_ms1,
            ms3: Vec::new(),
            contexts: read.contexts,
        };
        Ok((0, Vec::new(), scans))
    } else if file_is_mzml {
        let mut reader = MzMLReader::new(
            input::open_buf_maybe_gz(&spectrum_path).map_err(|e| format!("open mzML: {e}"))?,
        )
        .with_ms_level_range(mslevel, mslevel)
        .with_ms1_capture(true)
        .with_run_ms1_capture(capture_ms1)
        .with_product_scan_capture(capture_product)
        .with_scan_context_capture(capture_context)
        .with_noise_capture(capture_noise);
        let (errc, errs) = reader.read_with_ms1_chunked(CHUNK_SIZE, cap, |chunk, link| {
            let _ = tx.send((chunk, link));
        });
        let scans = RunScans {
            ms1: reader.take_run_ms1(),
            ms3: reader.take_product_scans(),
            contexts: reader.take_scan_contexts(),
        };
        Ok((errc, errs, scans))
    } else if file_is_d {
        #[cfg(feature = "timstof")]
        {
            let mut reader = input::TimsTofReader::open(&spectrum_path)
                .map_err(|e| format!("open Bruker .d: {e}"))?
                .with_dia_mobility_slices(dia_mobility_slices)
                .with_run_ms1_capture(capture_ms1)
                .with_ms1_mobility(true)
                .with_scan_context_capture(capture_context);
            let (errc, errs) = reader.read_with_ms1_chunked(CHUNK_SIZE, cap, |chunk, link| {
                let _ = tx.send((chunk, link));
            });
            let scans = RunScans {
                ms1: reader.take_run_ms1(),
                ms3: Vec::new(),
                contexts: reader.take_scan_contexts(),
            };
            Ok((errc, errs, scans))
        }
        #[cfg(not(feature = "timstof"))]
        {
            Err("this andes build has no Bruker .d (timsTOF) support; \
                 rebuild with `--features timstof`."
                .to_string())
        }
    } else {
        #[cfg(feature = "thermo")]
        {
            let mut reader = input::ThermoRawReader::open(&spectrum_path)
                .map_err(|e| format!("open Thermo .raw: {e}"))?
                .with_run_ms1_capture(capture_ms1)
                .with_product_scan_capture(capture_product)
                .with_scan_context_capture(capture_context)
                .with_noise_capture(capture_noise);
            let (errc, errs) = reader.read_with_ms1_chunked(CHUNK_SIZE, cap, |chunk, link| {
                let _ = tx.send((chunk, link));
            });
            let scans = RunScans {
                ms1: reader.take_run_ms1(),
                ms3: reader.take_product_scans(),
                contexts: reader.take_scan_contexts(),
            };
            Ok((errc, errs, scans))
        }
        #[cfg(not(feature = "thermo"))]
        {
            Err("this andes build has no Thermo .raw support; \
                 rebuild with `--features thermo`."
                .to_string())
        }
    }
}

/// Reads a file whose MS1 scans are linked to its MS2 scans (MS1 input, or the cross-link
/// precursor snap), on a parser thread, and searches it chunk by chunk.
fn search_linked<'a>(
    file_inputs: &FileInputs<'a>,
    outputs: &mut FileOutputs,
    prepared: &mut PreparedSearch<'a>,
    measurer: &mut FragmentErrorMeasurer,
    quant: &mut Option<QuantCollector>,
    locals: &FileLocals<'_>,
    run_scans: &mut RunScans,
) -> Result<ParseStats, Box<dyn std::error::Error>> {
    let FileInputs {
        cli,
        params,
        chimeric_active,
        mono_active,
        mono_params,
        keep_peaks,
        standalone,
        t_search_start,
        xl_ms3_capture,
        ..
    } = *file_inputs;
    let FileOutputs {
        all_spectra,
        all_queues,
        xl_snapped,
        mono_table,
        mono_stats,
        mono_dump,
        ..
    } = outputs;
    let FileLocals {
        input_path,
        remaining_cap,
        file_is_mzml,
        file_is_d,
        title_prefix,
        file_snap,
        ..
    } = *locals;
    let (tx, rx) = sync_channel::<(Vec<Spectrum>, Ms1Link)>(2);
    let spectrum_path = input_path.to_path_buf();
    let cap = remaining_cap;
    let mslevel = 2;
    let capture_ms1 = quant.as_ref().is_some_and(|q| q.settings.needs_ms1());
    let capture_context = quant.is_some();
    let capture_noise = quant.as_ref().is_some_and(|q| q.settings.plex.is_some());
    let index_ordered = chimeric_active && params.fragment_index_top_k > 0 && !standalone;
    let parallel_read = index_ordered
        && file_is_mzml
        && capture_ms1
        && cap == usize::MAX
        && !is_gzip_path(input_path);
    let threads = cli.threads;
    let capture_product = xl_ms3_capture;
    #[cfg(feature = "timstof")]
    let dia_mobility_slices = cli.dia_mobility_slices;
    let read = LinkedRead {
        spectrum_path,
        cap,
        mslevel,
        capture_ms1,
        capture_context,
        capture_noise,
        parallel_read,
        threads,
        capture_product,
        file_is_mzml,
        file_is_d,
        #[cfg(feature = "timstof")]
        dia_mobility_slices,
    };
    let parser_handle = thread::spawn(move || parse_linked(read, tx));

    let mut file_offset = 0usize;
    let mut ms1_linked = 0usize;
    let mut pending: Vec<(Spectrum, Option<Ms1Peaks>)> = Vec::new();
    for (mut chunk_spectra, chunk_link) in rx {
        if let Some(prefix) = title_prefix {
            prefix_spectrum_titles(&mut chunk_spectra, prefix);
        }
        if file_snap {
            for (i, spec) in chunk_spectra.iter_mut().enumerate() {
                let ms1 = chunk_link
                    .ms2_to_ms1
                    .get(i)
                    .copied()
                    .flatten()
                    .and_then(|m| chunk_link.ms1_peaks.get(m));
                if let Some(mz) = ms1.and_then(|p| {
                    search::precursor::monoisotopic::snap_precursor_mz(
                        p,
                        spec.precursor_mz,
                        search::crosslinks::PRECURSOR_SNAP_PPM,
                        search::crosslinks::PRECURSOR_SNAP_RATIO,
                    )
                }) {
                    spec.precursor_mz = mz;
                    *xl_snapped += 1;
                }
            }
        }
        if mono_active {
            correct_chunk(
                &mut chunk_spectra,
                &chunk_link,
                mono_params,
                mono_table,
                mono_stats,
                mono_dump.as_mut(),
            )
            .map_err(|e| format!("write --precursor-mono-dump: {e}"))?;
        }
        if index_ordered {
            file_offset += chunk_spectra.len();
            ms1_linked += chunk_link.ms1_peaks.len();
            for (i, spec) in chunk_spectra.into_iter().enumerate() {
                let ms1 = chunk_link
                    .ms2_to_ms1
                    .get(i)
                    .copied()
                    .flatten()
                    .and_then(|m| chunk_link.ms1_peaks.get(m))
                    .map(|peaks| {
                        let (lo, hi) =
                            search::engine::coisolation_ms1_bounds(&spec, prepared.params);
                        let a = peaks.partition_point(|p| p.0 < lo);
                        let b = peaks.partition_point(|p| p.0 <= hi);
                        peaks[a..b.max(a)].to_vec()
                    });
                pending.push((spec, ms1));
            }
            continue;
        }
        let offset = all_spectra.len();
        measurer.measure(prepared, chunk_spectra.iter(), None);
        if chimeric_active {
            let mut queues = prepared.run_chunk(&chunk_spectra, offset);
            search::engine::coisolation_pass(
                prepared,
                &chunk_spectra,
                &mut queues,
                prepared.params,
                &chunk_link,
                offset,
            );
            all_queues.extend(queues);
        } else if !standalone {
            let queues = prepared.run_chunk(&chunk_spectra, offset);
            all_queues.extend(queues);
        }
        file_offset += chunk_spectra.len();
        ms1_linked += chunk_link.ms1_peaks.len();
        for mut spec in chunk_spectra.into_iter() {
            if let Some(q) = quant.as_mut() {
                q.on_spectrum(&spec);
            }
            if !keep_peaks && !standalone {
                spec.peaks = Vec::new();
            }
            all_spectra.push(spec);
        }
        report_search_progress(all_spectra.len(), t_search_start);
    }
    if index_ordered {
        let span_da = index_chunk_span_da(cli.fragment_index_slice_da);
        measurer.measure(prepared, pending.iter().map(|e| &e.0), Some(span_da));
        pending.sort_by(|a, b| index_order_mass(&a.0).total_cmp(&index_order_mass(&b.0)));
        log::info!(
            target: "fragment-index",
            "{} spectra sorted by precursor mass, chunks of <= {} spectra and <= {:.0} Da",
            pending.len(),
            INDEX_CHUNK_SIZE,
            span_da
        );
        let base = all_spectra.len();
        let mut scored: Vec<Spectrum> = Vec::with_capacity(pending.len());
        let mut queues: Vec<search::psm::TopNQueue> = Vec::with_capacity(pending.len());
        let mut link = Ms1Link::default();
        let mut rest = pending;
        while !rest.is_empty() {
            let take = index_chunk_len(&rest, |e| index_order_mass(&e.0), span_da);
            let mut chunk: Vec<Spectrum> = Vec::with_capacity(take);
            for (spec, ms1) in rest.drain(..take) {
                chunk.push(spec);
                let ms1_idx = ms1.map(|peaks| {
                    link.ms1_peaks.push(peaks);
                    link.ms1_peaks.len() - 1
                });
                link.ms2_to_ms1.push(ms1_idx);
            }
            queues.extend(prepared.run_chunk(&chunk, base + scored.len()));
            scored.extend(chunk);
            report_search_progress(base + scored.len(), t_search_start);
        }
        search::engine::coisolation_pass(
            prepared,
            &scored,
            &mut queues,
            prepared.params,
            &link,
            base,
        );
        all_queues.extend(queues);
        for mut spec in scored.into_iter() {
            if let Some(q) = quant.as_mut() {
                q.on_spectrum(&spec);
            }
            if !keep_peaks {
                spec.peaks = Vec::new();
            }
            all_spectra.push(spec);
        }
    }
    let (err_count, first_errors, scans) = parser_handle
        .join()
        .map_err(|_| "MS1-linked parser thread panicked".to_string())??;
    *run_scans = scans;
    log::info!(
        target: if chimeric_active {
            "chimeric"
        } else if mono_active {
            "precursor-mono"
        } else {
            "xl"
        },
        "streamed {} MS2 spectra ({} MS1 scans linked) from {}",
        file_offset,
        ms1_linked,
        input_path.display()
    );
    Ok(ParseStats {
        error_count: err_count,
        first_errors,
    })
}

/// Reads an uncompressed mzML file in parallel and searches it with the fragment-ion index.
fn search_indexed<'a>(
    file_inputs: &FileInputs<'a>,
    outputs: &mut FileOutputs,
    prepared: &mut PreparedSearch<'a>,
    measurer: &mut FragmentErrorMeasurer,
    quant: &mut Option<QuantCollector>,
    locals: &FileLocals<'_>,
) -> Result<ParseStats, Box<dyn std::error::Error>> {
    let FileInputs {
        cli,
        ms_level_u32,
        keep_peaks,
        t_search_start,
        ..
    } = *file_inputs;
    let FileOutputs {
        all_spectra,
        all_queues,
        ..
    } = outputs;
    let FileLocals {
        input_path,
        remaining_cap,
        title_prefix,
        ..
    } = *locals;
    let (mut spectra, errors) =
        input::read_spectra_parallel(input_path, ms_level_u32, ms_level_u32, cli.threads)
            .map_err(|e| format!("read mzML {}: {e}", input_path.display()))?;
    spectra.truncate(remaining_cap);
    if let Some(prefix) = title_prefix {
        prefix_spectrum_titles(&mut spectra, prefix);
    }
    score_mass_ordered(
        prepared,
        spectra,
        index_chunk_span_da(cli.fragment_index_slice_da),
        all_spectra,
        all_queues,
        t_search_start,
        keep_peaks,
        quant.as_mut(),
        measurer,
    );
    Ok(ParseStats {
        error_count: errors.len(),
        first_errors: errors.into_iter().take(3).collect(),
    })
}

/// How a file is read and searched.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum ReadRoute {
    /// MS1-linked chunks on a parser thread (`search_linked`).
    Linked,
    /// An uncompressed mzML read in parallel pieces, MS2 only (`search_indexed`).
    Indexed,
    /// An uncompressed mzML read in parallel pieces with its MS1 (and MS3) scans, for
    /// quantification (`search_indexed_ms1`).
    IndexedMs1,
    /// Streamed on a parser thread (`search_enumerated`).
    Enumerated,
}

fn read_route(
    file_inputs: &FileInputs<'_>,
    quant: Option<&QuantCollector>,
    locals: &FileLocals<'_>,
    file_linked: bool,
) -> ReadRoute {
    let FileInputs {
        cli,
        params,
        standalone,
        ms_level_u32,
        xl_ms3_capture,
        ..
    } = *file_inputs;
    let path = locals.input_path;
    let parallel_mzml = params.fragment_index_top_k > 0
        && !standalone
        && !cli.dia
        && locals.file_is_mzml
        && !is_gzip_path(path);
    let needs_scans = quant.is_some_and(|q| q.settings.needs_ms1() || q.settings.needs_ms3());
    if file_linked {
        ReadRoute::Linked
    } else if parallel_mzml && !needs_scans {
        ReadRoute::Indexed
    } else if parallel_mzml
        && ms_level_u32 == 2
        && locals.remaining_cap == usize::MAX
        && !xl_ms3_capture
    {
        ReadRoute::IndexedMs1
    } else {
        ReadRoute::Enumerated
    }
}

/// Reads an uncompressed mzML file in parallel pieces with its MS1 (and MS3) scans and searches
/// it with the fragment-ion index: the spectra, run scans and contexts the streamed reader of
/// `search_enumerated` gives, in the same order. `None` (nothing read or searched) when a piece
/// fails to parse; the streamed reader then reads the file and counts its errors.
fn search_indexed_ms1<'a>(
    file_inputs: &FileInputs<'a>,
    outputs: &mut FileOutputs,
    prepared: &mut PreparedSearch<'a>,
    measurer: &mut FragmentErrorMeasurer,
    quant: &mut Option<QuantCollector>,
    locals: &FileLocals<'_>,
    run_scans: &mut RunScans,
) -> Result<Option<ParseStats>, Box<dyn std::error::Error>> {
    let FileInputs {
        cli,
        keep_peaks,
        t_search_start,
        ..
    } = *file_inputs;
    let FileLocals {
        input_path,
        title_prefix,
        ..
    } = *locals;
    let capture_ms1 = quant.as_ref().is_some_and(|q| q.settings.needs_ms1());
    let capture_ms3 = quant.as_ref().is_some_and(|q| q.settings.needs_ms3());
    let capture_noise = quant.as_ref().is_some_and(|q| q.settings.plex.is_some());
    let Some(read) = input::mzml::read_ms1_linked_parallel(
        input_path,
        cli.threads,
        quant.is_some(),
        capture_noise,
        capture_ms3,
    )
    .map_err(|e| format!("read mzML {}: {e}", input_path.display()))?
    else {
        return Ok(None);
    };
    let mut spectra = read.spectra;
    if let Some(prefix) = title_prefix {
        prefix_spectrum_titles(&mut spectra, prefix);
    }
    *run_scans = RunScans {
        ms1: if capture_ms1 { read.run_ms1 } else { Vec::new() },
        ms3: read.product_scans,
        contexts: read.contexts,
    };
    let FileOutputs {
        all_spectra,
        all_queues,
        ..
    } = outputs;
    score_mass_ordered(
        prepared,
        spectra,
        index_chunk_span_da(cli.fragment_index_slice_da),
        all_spectra,
        all_queues,
        t_search_start,
        keep_peaks,
        quant.as_mut(),
        measurer,
    );
    Ok(Some(ParseStats::default()))
}

/// Reads a file on a parser thread and searches it chunk by chunk.
fn search_enumerated<'a>(
    file_inputs: &FileInputs<'a>,
    outputs: &mut FileOutputs,
    prepared: &mut PreparedSearch<'a>,
    measurer: &mut FragmentErrorMeasurer,
    quant: &mut Option<QuantCollector>,
    locals: &FileLocals<'_>,
    run_scans: &mut RunScans,
) -> Result<ParseStats, Box<dyn std::error::Error>> {
    let FileInputs {
        cli,
        params,
        ms_level_u32,
        keep_peaks,
        standalone,
        t_search_start,
        xl_ms3_capture,
        ..
    } = *file_inputs;
    let FileOutputs {
        all_spectra,
        all_queues,
        ..
    } = outputs;
    let FileLocals {
        input_path,
        remaining_cap,
        file_is_mzml,
        file_is_raw,
        file_is_d,
        title_prefix,
        ..
    } = *locals;
    let (tx, rx) = sync_channel::<Vec<Spectrum>>(2);
    let spectrum_path = input_path.to_path_buf();
    let capture_ms1 = cli.dia || quant.as_ref().is_some_and(|q| q.settings.needs_ms1());
    let capture_ms3 = xl_ms3_capture || quant.as_ref().is_some_and(|q| q.settings.needs_ms3());
    let capture_context = quant.is_some();
    let capture_noise = quant.as_ref().is_some_and(|q| q.settings.plex.is_some());
    #[cfg(feature = "timstof")]
    let dia_mobility_slices = cli.dia_mobility_slices;
    let parser_handle = thread::spawn(
        move || -> Result<(ParseStats, RunScans), Box<dyn std::error::Error + Send + Sync>> {
            if file_is_mzml {
                let mut reader = MzMLReader::new(
                    input::open_buf_maybe_gz(&spectrum_path)
                        .map_err(|e| format!("open mzML: {e}"))?,
                )
                .with_ms_level_range(ms_level_u32, ms_level_u32)
                .with_run_ms1_capture(capture_ms1)
                .with_product_scan_capture(capture_ms3)
                .with_scan_context_capture(capture_context)
                .with_noise_capture(capture_noise);
                let mut stats = send_chunks(reader.by_ref(), CHUNK_SIZE, remaining_cap, tx);
                add_mzml_skipped(&mut stats, &reader);
                let scans = RunScans {
                    ms1: reader.take_run_ms1(),
                    ms3: reader.take_product_scans(),
                    contexts: reader.take_scan_contexts(),
                };
                Ok((stats, scans))
            } else if file_is_raw {
                #[cfg(feature = "thermo")]
                {
                    let mut reader = input::ThermoRawReader::open(&spectrum_path)
                        .map_err(|e| format!("open Thermo .raw: {e}"))?
                        .with_ms_level(Some(2))
                        .with_run_ms1_capture(capture_ms1)
                        .with_product_scan_capture(capture_ms3)
                        .with_scan_context_capture(capture_context)
                        .with_noise_capture(capture_noise);
                    let stats = send_chunks(reader.by_ref(), CHUNK_SIZE, remaining_cap, tx);
                    let scans = RunScans {
                        ms1: reader.take_run_ms1(),
                        ms3: reader.take_product_scans(),
                        contexts: reader.take_scan_contexts(),
                    };
                    Ok((stats, scans))
                }
                #[cfg(not(feature = "thermo"))]
                {
                    Err("this andes build has no Thermo .raw support; \
                         rebuild with `--features thermo` (and run with the \
                         .NET 8 runtime installed). mzML/MGF inputs work without it."
                        .into())
                }
            } else if file_is_d {
                #[cfg(feature = "timstof")]
                {
                    let mut reader = input::TimsTofReader::open(&spectrum_path)
                        .map_err(|e| format!("open Bruker .d: {e}"))?
                        .with_dia_mobility_slices(dia_mobility_slices)
                        .with_run_ms1_capture(capture_ms1)
                        .with_ms1_mobility(true)
                        .with_scan_context_capture(capture_context);
                    let stats = send_chunks(reader.by_ref(), CHUNK_SIZE, remaining_cap, tx);
                    let scans = RunScans {
                        ms1: reader.take_run_ms1(),
                        ms3: Vec::new(),
                        contexts: reader.take_scan_contexts(),
                    };
                    Ok((stats, scans))
                }
                #[cfg(not(feature = "timstof"))]
                {
                    Err("this andes build has no Bruker .d (timsTOF) support; \
                         rebuild with `--features timstof`. mzML/MGF inputs work \
                         without it."
                        .into())
                }
            } else {
                let reader = MgfReader::new(
                    input::open_buf_maybe_gz(&spectrum_path)
                        .map_err(|e| format!("open MGF: {e}"))?,
                );
                Ok((
                    send_chunks(reader, CHUNK_SIZE, remaining_cap, tx),
                    RunScans::default(),
                ))
            }
        },
    );

    let mass_ordered = (params.fragment_index_top_k > 0 || cli.dia) && !standalone;
    let mut pending: Vec<Spectrum> = Vec::new();
    for mut chunk in rx {
        if chunk.is_empty() {
            continue;
        }
        if let Some(prefix) = title_prefix {
            prefix_spectrum_titles(&mut chunk, prefix);
        }
        if mass_ordered {
            pending.extend(chunk);
            continue;
        }
        if !standalone {
            measurer.measure(prepared, chunk.iter(), None);
            let offset = all_spectra.len();
            let queues = prepared.run_chunk(&chunk, offset);
            all_queues.extend(queues);
        }
        for mut spec in chunk.into_iter() {
            if let Some(q) = quant.as_mut() {
                q.on_spectrum(&spec);
            }
            if !keep_peaks && !standalone {
                spec.peaks = Vec::new();
            }
            all_spectra.push(spec);
        }
        report_search_progress(all_spectra.len(), t_search_start);
    }
    let stats = match parser_handle.join() {
        Ok(Ok((stats, scans))) => {
            *run_scans = scans;
            stats
        }
        Ok(Err(e)) => return Err(format!("parser thread error: {e}").into()),
        Err(_) => return Err("parser thread panicked".into()),
    };
    if cli.dia && mass_ordered {
        if run_scans.ms1.is_empty() {
            return Err(format!(
                "--dia needs MS1 survey scans, and {} has none",
                input_path.display()
            )
            .into());
        }
        for s in &mut pending {
            s.precursor_charge = None;
        }
        prepared.set_dia_ms1(if quant.is_some() {
            run_scans.ms1.clone()
        } else {
            std::mem::take(&mut run_scans.ms1)
        });
        if file_is_d {
            recalibrate_dia_ms1(
                prepared,
                &pending,
                index_chunk_span_da(cli.fragment_index_slice_da),
            );
        }
    }
    if mass_ordered {
        score_mass_ordered(
            prepared,
            pending,
            index_chunk_span_da(cli.fragment_index_slice_da),
            all_spectra,
            all_queues,
            t_search_start,
            keep_peaks,
            quant.as_mut(),
            measurer,
        );
    }
    Ok(stats)
}
