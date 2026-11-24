//! Parallel mzML reading: the file is cut into byte pieces at `<spectrum` boundaries and
//! each piece is parsed by its own reader; results come back in file order.

use model::scan::{Ms1Scan, ScanContext};
use model::Spectrum;

use super::reader::MzMLReader;

const PARALLEL_PIECE_BYTES: u64 = 16 << 20;

fn next_spectrum_start(file: &std::fs::File, from: u64) -> std::io::Result<Option<u64>> {
    use std::io::{Read, Seek, SeekFrom};
    const PAT: &[u8] = b"<spectrum";
    let mut f = file;
    f.seek(SeekFrom::Start(from))?;
    let mut block = vec![0u8; 1 << 20];
    let keep = PAT.len();
    let mut carried = 0usize;
    let mut block_start = from;
    loop {
        let n = f.read(&mut block[carried..])?;
        let len = carried + n;
        if n == 0 {
            return Ok(None);
        }
        let hay = &block[..len];
        let mut i = 0;
        while i + PAT.len() < len {
            if &hay[i..i + PAT.len()] == PAT && hay[i + PAT.len()].is_ascii_whitespace() {
                return Ok(Some(block_start + i as u64));
            }
            i += 1;
        }
        let tail = keep.min(len);
        block.copy_within(len - tail..len, 0);
        block_start += (len - tail) as u64;
        carried = tail;
    }
}

pub fn read_spectra_parallel(
    path: &std::path::Path,
    ms_level_min: u32,
    ms_level_max: u32,
    threads: usize,
) -> std::io::Result<(Vec<Spectrum>, Vec<String>)> {
    read_spectra_in_pieces(
        path,
        ms_level_min,
        ms_level_max,
        threads,
        PARALLEL_PIECE_BYTES,
    )
}

fn spectrum_pieces(path: &std::path::Path, piece_bytes: u64) -> std::io::Result<Vec<(u64, u64)>> {
    let file = std::fs::File::open(path)?;
    let len = file.metadata()?.len();
    let n_cuts = len / piece_bytes;
    let mut starts: Vec<u64> = vec![0];
    for i in 1..=n_cuts {
        let from = (i * piece_bytes).max(*starts.last().unwrap() + 1);
        match next_spectrum_start(&file, from)? {
            Some(off) if off > *starts.last().unwrap() => starts.push(off),
            Some(_) => {}
            None => break,
        }
    }
    let n_pieces = starts.len();
    Ok((0..n_pieces)
        .map(|i| (starts[i], starts.get(i + 1).copied().unwrap_or(len)))
        .collect())
}

fn piece_bytes(
    path: &std::path::Path,
    pieces: &[(u64, u64)],
    i: usize,
) -> std::io::Result<Vec<u8>> {
    use std::io::{Read, Seek, SeekFrom};
    let (a, b) = pieces[i];
    let mut f = std::fs::File::open(path)?;
    f.seek(SeekFrom::Start(a))?;
    let mut bytes = Vec::with_capacity((b - a) as usize);
    f.take(b - a).read_to_end(&mut bytes)?;
    if i + 1 == pieces.len() {
        const END: &[u8] = b"</spectrumList";
        if let Some(p) = bytes.windows(END.len()).position(|w| w == END) {
            bytes.truncate(p);
        }
    }
    Ok(bytes)
}

fn map_pieces<T: Send>(
    n_pieces: usize,
    threads: usize,
    parse: impl Fn(usize) -> T + Sync,
) -> Vec<T> {
    use std::sync::atomic::{AtomicUsize, Ordering};
    let next = AtomicUsize::new(0);
    let mut results: Vec<(usize, T)> = std::thread::scope(|scope| {
        let workers: Vec<_> = (0..threads.clamp(1, n_pieces.max(1)))
            .map(|_| {
                scope.spawn(|| {
                    let mut out = Vec::new();
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        if i >= n_pieces {
                            break;
                        }
                        out.push((i, parse(i)));
                    }
                    out
                })
            })
            .collect();
        workers
            .into_iter()
            .flat_map(|w| w.join().expect("mzML parse worker panicked"))
            .collect()
    });
    results.sort_by_key(|(i, _)| *i);
    results.into_iter().map(|(_, r)| r).collect()
}

pub(super) fn read_spectra_in_pieces(
    path: &std::path::Path,
    ms_level_min: u32,
    ms_level_max: u32,
    threads: usize,
    piece_bytes_target: u64,
) -> std::io::Result<(Vec<Spectrum>, Vec<String>)> {
    let pieces = spectrum_pieces(path, piece_bytes_target)?;
    type PieceResult = std::io::Result<(Vec<Spectrum>, Vec<String>)>;
    let parse_piece = |i: usize| -> PieceResult {
        let bytes = piece_bytes(path, &pieces, i)?;
        let mut reader =
            MzMLReader::new(bytes.as_slice()).with_ms_level_range(ms_level_min, ms_level_max);
        let mut spectra = Vec::new();
        let mut errors = Vec::new();
        for r in reader.by_ref() {
            match r {
                Ok(s) => spectra.push(s),
                Err(e) => errors.push(e.to_string()),
            }
        }
        let kept = reader.skipped_errors().to_vec();
        let unlisted = reader.skipped_count() - kept.len();
        errors.extend(kept);
        errors.extend(std::iter::repeat_n(
            "malformed spectrum skipped".to_string(),
            unlisted,
        ));
        Ok((spectra, errors))
    };
    let mut spectra = Vec::new();
    let mut errors = Vec::new();
    for r in map_pieces(pieces.len(), threads, parse_piece) {
        let (s, e) = r?;
        spectra.extend(s);
        errors.extend(e);
    }
    Ok((spectra, errors))
}

#[derive(Debug, Default)]
pub struct Ms1LinkedRead {
    pub spectra: Vec<Spectrum>,
    pub ms2_to_ms1: Vec<Option<usize>>,
    pub run_ms1: Vec<Ms1Scan>,
    pub contexts: Vec<ScanContext>,
}

pub fn read_ms1_linked_parallel(
    path: &std::path::Path,
    threads: usize,
    capture_context: bool,
    capture_noise: bool,
) -> std::io::Result<Option<Ms1LinkedRead>> {
    read_ms1_linked_in_pieces(
        path,
        threads,
        capture_context,
        capture_noise,
        PARALLEL_PIECE_BYTES,
    )
}

/// Per MS2 scan of a piece: the MS1 scans before it and its FAIMS compensation voltage.
type Ms2Seen = (usize, Option<f32>);

pub(super) fn read_ms1_linked_in_pieces(
    path: &std::path::Path,
    threads: usize,
    capture_context: bool,
    capture_noise: bool,
    piece_bytes_target: u64,
) -> std::io::Result<Option<Ms1LinkedRead>> {
    let pieces = spectrum_pieces(path, piece_bytes_target)?;
    let parse_piece = |i: usize| -> std::io::Result<Option<(Ms1LinkedRead, Vec<Ms2Seen>)>> {
        let bytes = piece_bytes(path, &pieces, i)?;
        let mut reader = MzMLReader::new(bytes.as_slice())
            .with_ms_level_range(2, 2)
            .with_run_ms1_capture(true)
            .with_scan_context_capture(capture_context)
            .with_noise_capture(capture_noise);
        let mut read = Ms1LinkedRead::default();
        // Per MS2 scan: the MS1 scans of the piece before it, and its compensation voltage.
        let mut seen: Vec<Ms2Seen> = Vec::new();
        loop {
            match reader.pump() {
                Ok(Some(spectrum)) => {
                    seen.push((reader.run_ms1.len(), reader.last_faims_cv));
                    read.spectra.push(spectrum);
                }
                Ok(None) => break,
                Err(_) => return Ok(None),
            }
        }
        read.run_ms1 = reader.take_run_ms1();
        read.contexts = reader.take_scan_contexts();
        Ok(Some((read, seen)))
    };
    let mut out = Ms1LinkedRead::default();
    let mut seen: Vec<Ms2Seen> = Vec::new();
    for r in map_pieces(pieces.len(), threads, parse_piece) {
        let Some((piece, piece_seen)) = r? else {
            return Ok(None);
        };
        let base = out.run_ms1.len();
        seen.extend(piece_seen.into_iter().map(|(n, cv)| (base + n, cv)));
        out.spectra.extend(piece.spectra);
        out.run_ms1.extend(piece.run_ms1);
        out.contexts.extend(piece.contexts);
    }
    // Link in acquisition order across the pieces.
    let mut tracker = crate::ms1::Ms1Tracker::default();
    let mut next = 0usize;
    out.ms2_to_ms1 = seen
        .into_iter()
        .map(|(n, cv)| {
            while next < n {
                tracker.push(next, out.run_ms1[next].faims_cv);
                next += 1;
            }
            tracker.link(cv)
        })
        .collect();
    Ok(Some(out))
}
