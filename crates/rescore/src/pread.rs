//! Positioned reads of spill files on every platform, so several threads can read one file at
//! once: each read names its own offset and never depends on a shared file cursor.

use std::fs::File;
use std::io;

/// Reads exactly `buf.len()` bytes of a file starting at byte `offset`.
pub(crate) trait ReadAt {
    /// Fills `buf` from byte `offset`; fails with `UnexpectedEof` if the file ends first.
    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> io::Result<()>;
}

#[cfg(unix)]
impl ReadAt for File {
    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> io::Result<()> {
        std::os::unix::fs::FileExt::read_exact_at(self, buf, offset)
    }
}

#[cfg(windows)]
impl ReadAt for File {
    fn read_exact_at(&self, mut buf: &mut [u8], mut offset: u64) -> io::Result<()> {
        use std::os::windows::fs::FileExt;
        while !buf.is_empty() {
            match self.seek_read(buf, offset) {
                Ok(0) => {
                    return Err(io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        "failed to fill whole buffer",
                    ))
                }
                Ok(n) => {
                    buf = &mut buf[n..];
                    offset += n as u64;
                }
                Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                Err(e) => return Err(e),
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn reads_at_an_offset_and_fails_past_the_end() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("f");
        File::create(&path).unwrap().write_all(b"0123456789").unwrap();
        let file = File::open(&path).unwrap();
        let mut buf = [0u8; 4];
        file.read_exact_at(&mut buf, 3).unwrap();
        assert_eq!(&buf, b"3456");
        file.read_exact_at(&mut buf, 0).unwrap();
        assert_eq!(&buf, b"0123");
        let err = file.read_exact_at(&mut buf, 8).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    }
}
