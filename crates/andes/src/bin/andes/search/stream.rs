//! The streamed read of an input file on a parser thread (mzML, MGF, Thermo .raw, Bruker .d),
//! and starting the first file's read before the candidate index is built, so that the read
//! and the build overlap.

use std::path::{Path, PathBuf};
use std::sync::mpsc::{sync_channel, Receiver, SyncSender};
use std::thread::{self, JoinHandle};

use input::{MgfReader, MzMLReader};
use model::Spectrum;
use search::SearchParams;

use super::files::file_reads_linked;
#[cfg(feature = "thermo")]
use super::files::THERMO_READERS_MAX;
use super::CHUNK_SIZE;
use crate::cli::Cli;
use crate::quant::run::{QuantCollector, RunScans};
use crate::spectra::{add_mzml_skipped, input_format_flags, send_chunks, ParseStats};

/// What a streamed parser thread returns: its parse statistics and the run scans it captured.
pub(super) type StreamResult =
    Result<(ParseStats, RunScans), Box<dyn std::error::Error + Send + Sync>>;

/// Chunks a streamed read may send ahead of the search.
pub(super) const STREAM_CHANNEL_CHUNKS: usize = 2;

/// Chunks the first file's early read may hold before the search takes them: the search takes
/// them only after the candidate index is built.
const PREFETCH_CHANNEL_CHUNKS: usize = 1024;

/// A parser thread and the spectrum chunks it sends.
pub(super) struct StreamedRead {
    pub(super) handle: JoinHandle<StreamResult>,
    pub(super) rx: Receiver<Vec<Spectrum>>,
}

/// What the streamed parser of one file reads and captures.
pub(super) struct StreamSettings {
    spectrum_path: PathBuf,
    ms_level_u32: u32,
    remaining_cap: usize,
    capture_ms1: bool,
    capture_ms3: bool,
    capture_context: bool,
    capture_noise: bool,
    file_is_mzml: bool,
    file_is_raw: bool,
    file_is_d: bool,
    #[cfg(feature = "thermo")]
    threads: usize,
    #[cfg(feature = "timstof")]
    dia_mobility_slices: usize,
}

impl StreamSettings {
    /// The read of `path` at `ms_level_u32`, at most `remaining_cap` spectra, capturing the scans
    /// the quantification (and the cross-link MS3 search) needs.
    pub(super) fn of(
        cli: &Cli,
        xl_ms3_capture: bool,
        quant: Option<&QuantCollector>,
        (path, ms_level_u32, remaining_cap): (&Path, u32, usize),
    ) -> Self {
        let (file_is_mzml, file_is_raw, file_is_d, _) = input_format_flags(path);
        StreamSettings {
            spectrum_path: path.to_path_buf(),
            ms_level_u32,
            remaining_cap,
            capture_ms1: cli.dia || quant.is_some_and(|q| q.settings.needs_ms1()),
            capture_ms3: xl_ms3_capture || quant.is_some_and(|q| q.settings.needs_ms3()),
            capture_context: quant.is_some(),
            capture_noise: quant.is_some_and(|q| q.settings.plex.is_some()),
            file_is_mzml,
            file_is_raw,
            file_is_d,
            #[cfg(feature = "thermo")]
            threads: cli.threads,
            #[cfg(feature = "timstof")]
            dia_mobility_slices: cli.dia_mobility_slices,
        }
    }
}

/// Starts the streamed read of a file on its parser thread; up to `bound` chunks wait for the
/// search.
pub(super) fn spawn_streamed(settings: StreamSettings, bound: usize) -> StreamedRead {
    let (tx, rx) = sync_channel::<Vec<Spectrum>>(bound);
    let handle = thread::spawn(move || parse_streamed(settings, tx));
    StreamedRead { handle, rx }
}

/// Starts reading the first input when the search will stream it and score it in mass order, as
/// it does a Thermo .raw or Bruker .d file: the read then runs while the candidate index is
/// built, and the search takes the chunks afterwards. `None` for every other input and mode.
pub(super) fn prefetch_first_file(
    cli: &Cli,
    (params, standalone, ms1_input, xl_ms3_capture): (&SearchParams, bool, bool, bool),
    (bench_mode, bench_cap): (bool, usize),
    quant: Option<&QuantCollector>,
) -> Option<StreamedRead> {
    let path = cli.spectrum.first()?;
    let (_, is_raw, is_d, _) = input_format_flags(path);
    let mass_ordered = params.fragment_index_top_k > 0 && !standalone && !cli.dia;
    if !(is_raw || is_d) || !mass_ordered || file_reads_linked(cli, ms1_input, path) {
        return None;
    }
    let remaining_cap = if bench_mode { bench_cap } else { usize::MAX };
    let settings = StreamSettings::of(
        cli,
        xl_ms3_capture,
        quant,
        (path, cli.ms_level as u32, remaining_cap),
    );
    Some(spawn_streamed(settings, PREFETCH_CHANNEL_CHUNKS))
}

fn parse_streamed(settings: StreamSettings, tx: SyncSender<Vec<Spectrum>>) -> StreamResult {
    let StreamSettings {
        spectrum_path,
        ms_level_u32,
        remaining_cap,
        capture_ms1,
        capture_ms3,
        capture_context,
        capture_noise,
        file_is_mzml,
        file_is_raw,
        file_is_d,
        #[cfg(feature = "thermo")]
        threads,
        #[cfg(feature = "timstof")]
        dia_mobility_slices,
    } = settings;
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
            scan_mobility: Vec::new(),
        };
        Ok((stats, scans))
    } else if file_is_raw {
        #[cfg(feature = "thermo")]
        {
            let mut reader = input::ThermoRawReader::open(&spectrum_path)
                .map_err(|e| format!("open Thermo .raw: {e}"))?
                .with_parallel_reads(threads.clamp(1, THERMO_READERS_MAX))
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
                scan_mobility: Vec::new(),
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
            let scan_mobility = reader.dia_scan_mobility();
            let scans = RunScans {
                ms1: reader.take_run_ms1(),
                ms3: Vec::new(),
                contexts: reader.take_scan_contexts(),
                scan_mobility,
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
}
