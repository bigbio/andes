//! Reading a large text table (a PIN) in blocks of whole lines: a reader thread fills the next
//! block while the caller works on the current one, so the caller can parse or rewrite each
//! block's lines in parallel and still see them in file order.

use std::fs::File;
use std::io::Read;
use std::path::Path;
use std::sync::mpsc::sync_channel;

/// Bytes read per block.
const BLOCK_BYTES: usize = 64 << 20;
/// Lines per parallel work unit of a block.
pub(crate) const LINES_PER_TASK: usize = 4096;

/// Splits `text` into its lines the way [`std::io::BufRead::lines`] does (a trailing `\n` or
/// `\r\n` taken off each), leaving out empty lines.
fn push_lines<'a>(text: &'a str, out: &mut Vec<&'a str>) {
    for line in text.split_terminator('\n') {
        let line = line.strip_suffix('\r').unwrap_or(line);
        if !line.is_empty() {
            out.push(line);
        }
    }
}

/// Reads `path`: `on_header` gets its first line, then `on_block` gets the following non-empty
/// lines block by block, in file order, with the index of the block's first line among them.
/// Lines end as [`std::io::BufRead::lines`] ends them; a file that is not UTF-8 is an error.
pub(crate) fn read_blocks(
    path: &Path,
    on_header: impl FnOnce(&str) -> Result<(), String>,
    mut on_block: impl FnMut(usize, &[&str]) -> Result<(), String>,
) -> Result<(), String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
    let mut f = File::open(path).map_err(|e| err(&e))?;
    let (tx, rx) = sync_channel::<std::io::Result<Vec<u8>>>(1);
    std::thread::scope(|scope| {
        scope.spawn(move || {
            let mut carry: Vec<u8> = Vec::new();
            loop {
                let mut buf = std::mem::take(&mut carry);
                let start = buf.len();
                buf.resize(start + BLOCK_BYTES, 0);
                let mut filled = start;
                while filled < buf.len() {
                    match f.read(&mut buf[filled..]) {
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
        });
        let mut header: Option<_> = Some(on_header);
        let mut first = 0usize;
        for block in rx {
            let block = block.map_err(|e| err(&e))?;
            let text = std::str::from_utf8(&block)
                .map_err(|_| err(&"stream did not contain valid UTF-8"))?;
            let mut body = text;
            if let Some(h) = header.take() {
                let (line, rest) = text.split_once('\n').unwrap_or((text, ""));
                h(line.strip_suffix('\r').unwrap_or(line))?;
                body = rest;
            }
            let mut lines: Vec<&str> = Vec::new();
            push_lines(body, &mut lines);
            on_block(first, &lines)?;
            first += lines.len();
        }
        match header {
            Some(_) => Err(err(&"empty PIN")),
            None => Ok(()),
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufRead;

    fn lines_of(path: &Path) -> (String, Vec<(usize, String)>) {
        let mut head = String::new();
        let mut rows = Vec::new();
        read_blocks(
            path,
            |h| {
                head = h.to_string();
                Ok(())
            },
            |first, lines| {
                rows.extend(lines.iter().enumerate().map(|(i, l)| (first + i, l.to_string())));
                Ok(())
            },
        )
        .unwrap();
        (head, rows)
    }

    #[test]
    fn blocks_give_the_lines_bufread_gives_without_the_empty_ones() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("t.pin");
        let text = "h1\th2\r\na\t1\n\nb\t2\r\nc\t3";
        std::fs::write(&p, text).unwrap();
        let (head, rows) = lines_of(&p);
        let want: Vec<String> = std::io::BufReader::new(text.as_bytes())
            .lines()
            .map(Result::unwrap)
            .collect();
        assert_eq!(head, want[0]);
        let body: Vec<(usize, String)> = want[1..]
            .iter()
            .filter(|l| !l.is_empty())
            .cloned()
            .enumerate()
            .collect();
        assert_eq!(rows, body);
    }

    #[test]
    fn an_empty_file_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("e.pin");
        std::fs::write(&p, "").unwrap();
        assert!(read_blocks(&p, |_| Ok(()), |_, _| Ok(())).is_err());
    }
}
