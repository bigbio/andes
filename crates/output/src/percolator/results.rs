//! Reading Percolator result files in blocks of whole lines: a reader thread fills the next
//! block while the rows of the current one are parsed in parallel, and the parsed rows reach
//! the caller block by block in file order. Rows, errors and their order are those of
//! [`PsmResultReader`](super::PsmResultReader) over the same input.

use std::collections::HashMap;
use std::io::Read;
use std::sync::mpsc::{sync_channel, SyncSender};

use rayon::prelude::*;

use super::{PercolatorError, PercolatorPsm, PsmResultRow, ResultColumns};

/// Bytes read per block.
const BLOCK_BYTES: usize = 64 << 20;
/// Rows of one parallel work unit.
const ROWS_PER_TASK: usize = 4096;
/// Rows of in-memory text parsed before they are inserted into the result map.
const ROWS_PER_PASS: usize = 1 << 18;

fn utf8_error() -> PercolatorError {
    PercolatorError::Io(std::io::Error::new(
        std::io::ErrorKind::InvalidData,
        "stream did not contain valid UTF-8",
    ))
}

/// The first line of `text` and the text after it; the line ends as
/// [`std::io::BufRead::read_line`] ends it, without its `\n` and a `\r` before that.
fn split_header(text: &str) -> (&str, &str) {
    match text.split_once('\n') {
        Some((line, rest)) => (line.strip_suffix('\r').unwrap_or(line), rest),
        None => (text, ""),
    }
}

/// Pushes the non-empty lines of `text`, each without its `\n` and a `\r` before that (a last
/// line without `\n` keeps a trailing `\r`, as [`std::io::BufRead::read_line`] leaves it).
fn push_lines<'a>(text: &'a str, out: &mut Vec<&'a str>) {
    for raw in text.split_inclusive('\n') {
        let line = match raw.strip_suffix('\n') {
            Some(l) => l.strip_suffix('\r').unwrap_or(l),
            None => raw,
        };
        if !line.is_empty() {
            out.push(line);
        }
    }
}

/// Parses `lines` in parallel into the values `map` makes of them, in order; stops at the first
/// malformed row and returns its error with the values of the rows before it.
fn parse_rows<T: Send>(
    columns: &ResultColumns,
    lines: &[&str],
    map: &(impl Fn(PsmResultRow<'_>) -> T + Sync),
) -> (Vec<T>, Option<PercolatorError>) {
    let parsed: Vec<Result<T, PercolatorError>> = lines
        .par_iter()
        .with_min_len(ROWS_PER_TASK)
        .map(|line| columns.row(line).map(map))
        .collect();
    if parsed.iter().all(Result::is_ok) {
        return (parsed.into_iter().flatten().collect(), None);
    }
    let mut values = Vec::with_capacity(parsed.len());
    for r in parsed {
        match r {
            Ok(v) => values.push(v),
            Err(e) => return (values, Some(e)),
        }
    }
    (values, None)
}

/// Sends `reader`'s bytes in blocks of whole lines (the last block may end without `\n`).
fn fill_blocks(mut reader: impl Read, tx: &SyncSender<std::io::Result<Vec<u8>>>) {
    let mut carry: Vec<u8> = Vec::new();
    loop {
        let mut buf = std::mem::take(&mut carry);
        let start = buf.len();
        buf.resize(start + BLOCK_BYTES, 0);
        let mut filled = start;
        while filled < buf.len() {
            match reader.read(&mut buf[filled..]) {
                Ok(0) => break,
                Ok(n) => filled += n,
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
                Err(e) => {
                    let _ = tx.send(Err(e));
                    return;
                }
            }
        }
        buf.truncate(filled);
        let eof = filled < start + BLOCK_BYTES;
        if !eof {
            match buf.iter().rposition(|&b| b == b'\n') {
                Some(cut) => carry = buf.split_off(cut + 1),
                None => {
                    carry = buf;
                    continue;
                }
            }
        }
        if !buf.is_empty() && tx.send(Ok(buf)).is_err() {
            return;
        }
        if eof {
            return;
        }
    }
}

/// The text of a block: all of it when it is UTF-8, else its lines before the first invalid
/// byte, with `true`.
fn block_text(block: &[u8]) -> (&str, bool) {
    match std::str::from_utf8(block) {
        Ok(text) => (text, false),
        Err(e) => {
            let valid = &block[..e.valid_up_to()];
            let end = valid.iter().rposition(|&b| b == b'\n').map_or(0, |i| i + 1);
            (std::str::from_utf8(&valid[..end]).unwrap_or_default(), true)
        }
    }
}

/// Reads a Percolator result file from `reader`: `map` turns each data row into a `T` (rows in
/// parallel), and `on_block` gets the values of each block of rows in file order. A malformed
/// row, or bytes that are not UTF-8, end the read with the error [`PsmResultReader`] gives on
/// that row; the values of the rows before it reach `on_block` first.
///
/// [`PsmResultReader`]: super::PsmResultReader
pub fn read_psm_results<T: Send>(
    reader: impl Read + Send,
    map: impl Fn(PsmResultRow<'_>) -> T + Sync,
    mut on_block: impl FnMut(Vec<T>) -> Result<(), PercolatorError>,
) -> Result<(), PercolatorError> {
    let (tx, rx) = sync_channel::<std::io::Result<Vec<u8>>>(1);
    std::thread::scope(|scope| {
        scope.spawn(move || fill_blocks(reader, &tx));
        let mut columns: Option<ResultColumns> = None;
        for block in rx {
            let block = block.map_err(PercolatorError::Io)?;
            let (text, invalid) = block_text(&block);
            let (cols, body) = match columns {
                Some(c) => (c, text),
                None if invalid && text.is_empty() => return Err(utf8_error()),
                None => {
                    let (header, body) = split_header(text);
                    let c = ResultColumns::from_header(header)?;
                    columns = Some(c);
                    (c, body)
                }
            };
            let mut lines: Vec<&str> = Vec::new();
            push_lines(body, &mut lines);
            let (values, err) = parse_rows(&cols, &lines, &map);
            on_block(values)?;
            if let Some(e) = err {
                return Err(e);
            }
            if invalid {
                return Err(utf8_error());
            }
        }
        if columns.is_none() {
            ResultColumns::from_header("")?;
        }
        Ok(())
    })
}

/// A result row as a [`PercolatorPsm`].
fn psm_of(row: PsmResultRow<'_>) -> PercolatorPsm {
    PercolatorPsm {
        pooled_q_value: row.pooled_q,
        psm_id: row.psm_id.to_string(),
        q_value: row.q_value,
        pep: row.pep,
        score: row.score,
        peptide: row.peptide.to_string(),
        proteins: row.proteins(),
    }
}

/// [`super::parse_psm_results`]: the rows are parsed in parallel and inserted in file order,
/// so the first duplicate PSMId and the first malformed row are reported as a serial read
/// reports them.
pub(super) fn parse_text(text: &str) -> Result<HashMap<String, PercolatorPsm>, PercolatorError> {
    let (header, body) = split_header(text);
    let columns = ResultColumns::from_header(header)?;
    let mut lines: Vec<&str> = Vec::new();
    push_lines(body, &mut lines);
    let mut map = HashMap::new();
    for pass in lines.chunks(ROWS_PER_PASS) {
        let (psms, err) = parse_rows(&columns, pass, &psm_of);
        for psm in psms {
            if let Some(old) = map.insert(psm.psm_id.clone(), psm) {
                return Err(PercolatorError::DuplicatePsmId { psm_id: old.psm_id });
            }
        }
        if let Some(e) = err {
            return Err(e);
        }
    }
    Ok(map)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::percolator::PsmResultReader;

    /// Every row the serial reader gives, as (PSMId, q, pep, score, peptide, proteins, run).
    type Row = (
        String,
        f64,
        f64,
        Option<f64>,
        String,
        String,
        Option<String>,
    );

    fn owned(r: PsmResultRow<'_>) -> Row {
        (
            r.psm_id.to_string(),
            r.q_value,
            r.pep,
            r.score,
            r.peptide.to_string(),
            r.proteins(),
            r.run.map(str::to_string),
        )
    }

    fn serial(bytes: &[u8]) -> (Vec<Row>, Option<String>) {
        let mut rows = Vec::new();
        let mut reader = match PsmResultReader::new(bytes) {
            Ok(r) => r,
            Err(e) => return (rows, Some(e.to_string())),
        };
        loop {
            match reader.next_row() {
                Ok(Some(r)) => rows.push(owned(r)),
                Ok(None) => return (rows, None),
                Err(e) => return (rows, Some(e.to_string())),
            }
        }
    }

    fn blocks(bytes: &[u8]) -> (Vec<Row>, Option<String>) {
        let mut rows = Vec::new();
        let r = read_psm_results(bytes, owned, |b| {
            rows.extend(b);
            Ok(())
        });
        (rows, r.err().map(|e| e.to_string()))
    }

    const HEADER: &str = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds";

    #[test]
    fn block_reads_give_the_rows_and_errors_of_the_serial_reader() {
        let cases: Vec<Vec<u8>> = vec![
            Vec::new(),
            HEADER.into(),
            format!(
                "{HEADER}\r\na\t1\t0.01\t0.1\tK.PEPK.A\tP1\tP2\r\n\r\nb\t2\t0.02\t0.2\tR.PK.A\tP3"
            )
            .into(),
            format!("{HEADER}\na\t1\t0.01\t0.1\tK.PEPK.A\tP1\nb\tx\t0.5\n").into(),
            format!("{HEADER}\na\t1\tnan\t0.1\tK.PEPK.A\tP1\nb\t1\t0.1\t0.1\tK.K.A\tP1\n").into(),
            format!("{HEADER}\na\t1\t0.01\t0.1\tK.PEPK.A\tP1\r").into(),
            [
                format!("{HEADER}\na\t1\t0.01\t0.1\tK.PEPK.A\tP1\nb\t1\t0.01\t0.1\tK.P").as_bytes(),
                &[0xff, 0xfe],
                b"K.A\tP1\n",
            ]
            .concat(),
            [b"PSMId\t".as_slice(), &[0xff], b"\n"].concat(),
            b"wrong\theader\na\t1\n".to_vec(),
        ];
        for c in &cases {
            assert_eq!(blocks(c), serial(c), "{:?}", String::from_utf8_lossy(c));
        }
    }

    #[test]
    fn rows_across_block_boundaries_arrive_whole_and_in_order() {
        let n = BLOCK_BYTES / 40 + 1000;
        let mut text = format!("{HEADER}\n");
        for i in 0..n {
            text.push_str(&format!("id{i}\t{i}\t0.01\t0.1\tK.PEP{i}K.A\tP{i}\n"));
        }
        let (rows, err) = blocks(text.as_bytes());
        assert!(err.is_none());
        assert_eq!(rows.len(), n);
        assert!(rows
            .iter()
            .enumerate()
            .all(|(i, r)| r.0 == format!("id{i}")));
    }

    #[test]
    fn the_text_parse_reports_the_first_duplicate_before_a_later_malformed_row() {
        let text =
            format!("{HEADER}\na\t1\t0.01\t0.1\tK.PK.A\tP1\na\t1\t0.01\t0.1\tK.PK.A\tP1\nb\n");
        assert!(matches!(
            parse_text(&text),
            Err(PercolatorError::DuplicatePsmId { psm_id }) if psm_id == "a"
        ));
    }
}
