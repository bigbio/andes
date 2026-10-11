//! Private, lossless scratch storage for a run's MS1 scans. No run stays
//! resident while Percolator trains; quantification loads one run at a time.
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;

use model::scan::Ms1Scan;
use quant::ms1_index::Ms1RunIndex;
use tempfile::NamedTempFile;

pub(crate) struct Ms1Cache {
    file: NamedTempFile,
}

const MAGIC: &[u8; 8] = b"ANDESRC1";
const MAJOR: u16 = 1;

impl Ms1Cache {
    pub(super) fn store(index: &Ms1RunIndex) -> io::Result<Self> {
        let file = tempfile::Builder::new().prefix("andes-ms1-").tempfile()?;
        let mut out = BufWriter::with_capacity(1024 * 1024, file.as_file());
        write_body(index, &mut out)?;
        out.flush()?;
        drop(out);
        Ok(Self { file })
    }

    pub(super) fn load(&self) -> io::Result<Ms1RunIndex> {
        let file = self.file.reopen()?;
        let remaining = file.metadata()?.len();
        let mut input = BufReader::with_capacity(1024 * 1024, file);
        read_body(&mut input, Some(remaining))
    }

    /// Writes a persistent run cache: magic, major and minor version, then the
    /// zstd-compressed body. The file appears at `path` only once complete.
    pub(crate) fn store_at(index: &Ms1RunIndex, path: &Path) -> io::Result<()> {
        let tmp = path.with_extension("cache.partial");
        {
            let mut f = BufWriter::new(std::fs::File::create(&tmp)?);
            f.write_all(MAGIC)?;
            f.write_all(&MAJOR.to_le_bytes())?;
            f.write_all(&0u16.to_le_bytes())?;
            let mut z = zstd::Encoder::new(f, 3)?;
            write_body(index, &mut z)?;
            z.finish()?.flush()?;
        }
        std::fs::rename(&tmp, path)
    }

    /// Reads a persistent run cache; every error names `path`.
    pub(crate) fn load_from(path: &Path) -> io::Result<Ms1RunIndex> {
        let named = |e: io::Error| io::Error::new(e.kind(), format!("{}: {e}", path.display()));
        let mut f = BufReader::new(std::fs::File::open(path).map_err(named)?);
        let mut head = [0u8; 12];
        if f.read_exact(&mut head).is_err() || &head[..8] != MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("{}: not an andes run cache", path.display()),
            ));
        }
        let major = u16::from_le_bytes([head[8], head[9]]);
        if major != MAJOR {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "{}: run cache version {major}, this andes reads version {MAJOR}",
                    path.display()
                ),
            ));
        }
        let mut z = zstd::Decoder::new(f).map_err(named)?;
        read_body(&mut z, None).map_err(named)
    }
}

fn write_body(index: &Ms1RunIndex, out: &mut impl Write) -> io::Result<()> {
    out.write_all(&(index.len() as u64).to_le_bytes())?;
    for i in 0..index.len() {
        let scan = index.scan(i);
        out.write_all(&index.rt(i).to_le_bytes())?;
        out.write_all(&(scan.id.len() as u64).to_le_bytes())?;
        out.write_all(scan.id.as_bytes())?;
        out.write_all(&scan.faims_cv.map_or(f32::NAN, |v| v).to_le_bytes())?;
        out.write_all(&(index.peaks(i).len() as u64).to_le_bytes())?;
        for &(mz, intensity) in index.peaks(i) {
            out.write_all(&mz.to_le_bytes())?;
            out.write_all(&intensity.to_le_bytes())?;
        }
    }
    Ok(())
}

/// Decodes a body. `remaining` is the exact byte length when known (the
/// temporary file); `None` for a compressed stream, where counts are capped at
/// `u32::MAX` and trailing data is detected by reading past the end.
fn read_body(input: &mut impl Read, remaining: Option<u64>) -> io::Result<Ms1RunIndex> {
    let known = remaining.is_some();
    let mut remaining = remaining.unwrap_or(u64::MAX);
    let cap = |left: u64, unit: u64| if known { left / unit } else { u32::MAX as u64 };
    let count = u64::from_le_bytes(read::<8>(input, &mut remaining)?);
    let count = bounded_count(count, cap(remaining, 16))?;
    let mut scans = Vec::with_capacity(count.min(PREALLOC_LIMIT));
    for _ in 0..count {
        let rt = f64::from_le_bytes(read::<8>(input, &mut remaining)?);
        let id_len = u64::from_le_bytes(read::<8>(input, &mut remaining)?);
        let id_len = bounded_count(id_len, cap(remaining, 1))?;
        let mut id = Vec::new();
        let got = input.take(id_len as u64).read_to_end(&mut id)?;
        if got != id_len {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "truncated MS1 cache",
            ));
        }
        remaining = remaining.saturating_sub(id_len as u64);
        let id = String::from_utf8(id)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "MS1 cache scan id"))?;
        let faims_cv = f32::from_le_bytes(read::<4>(input, &mut remaining)?);
        let faims_cv = (!faims_cv.is_nan()).then_some(faims_cv);
        let n = u64::from_le_bytes(read::<8>(input, &mut remaining)?);
        let n = bounded_count(n, cap(remaining, 12))?;
        let mut peaks = Vec::with_capacity(n.min(PREALLOC_LIMIT));
        for _ in 0..n {
            let mz = f64::from_le_bytes(read::<8>(input, &mut remaining)?);
            let intensity = f32::from_le_bytes(read::<4>(input, &mut remaining)?);
            peaks.push((mz, intensity));
        }
        scans.push(Ms1Scan {
            rt,
            peaks,
            id,
            faims_cv,
        });
    }
    let trailing = if known {
        remaining != 0
    } else {
        input.read(&mut [0u8; 1])? != 0
    };
    if trailing {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "trailing MS1 cache data",
        ));
    }
    Ok(Ms1RunIndex::new(scans))
}

/// Upper bound on an up-front allocation driven by an untrusted length.
const PREALLOC_LIMIT: usize = 1 << 16;

fn bounded_count(n: u64, max: u64) -> io::Result<usize> {
    usize::try_from(n)
        .ok()
        .filter(|_| n <= max)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid MS1 cache length"))
}

fn read<const N: usize>(input: &mut impl Read, remaining: &mut u64) -> io::Result<[u8; N]> {
    let mut bytes = [0; N];
    input.read_exact(&mut bytes)?;
    *remaining = remaining
        .checked_sub(N as u64)
        .ok_or_else(|| io::Error::new(io::ErrorKind::UnexpectedEof, "truncated MS1 cache"))?;
    Ok(bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lossless_repeatable_round_trip_and_cleanup() {
        let index = Ms1RunIndex::new(vec![
            Ms1Scan {
                rt: 3.123456789,
                peaks: vec![(600.123456789, 1234.5678), (700.0, -0.0)],
                id: "scan=2".into(),
                faims_cv: Some(-45.0),
            },
            Ms1Scan {
                rt: 1.0,
                peaks: vec![],
                ..Ms1Scan::default()
            },
        ]);
        let cache = Ms1Cache::store(&index).unwrap();
        let path = cache.file.path().to_owned();
        for _ in 0..2 {
            let loaded = cache.load().unwrap();
            assert_eq!(loaded.len(), index.len());
            for i in 0..index.len() {
                assert_eq!(loaded.rt(i).to_bits(), index.rt(i).to_bits());
                for (a, b) in loaded.peaks(i).iter().zip(index.peaks(i)) {
                    assert_eq!(a.0.to_bits(), b.0.to_bits());
                    assert_eq!(a.1.to_bits(), b.1.to_bits());
                }
                assert_eq!(loaded.peaks(i).len(), index.peaks(i).len());
                assert_eq!(loaded.scan(i).id, index.scan(i).id);
                assert_eq!(loaded.scan(i).faims_cv, index.scan(i).faims_cv);
            }
        }
        drop(cache);
        assert!(!path.exists());
    }

    #[test]
    fn empty_and_truncated_caches() {
        let empty = Ms1Cache::store(&Ms1RunIndex::default()).unwrap();
        assert!(empty.load().unwrap().is_empty());
        let index = Ms1RunIndex::new(vec![Ms1Scan {
            rt: 1.0,
            peaks: vec![(500.0, 1.0)],
            ..Ms1Scan::default()
        }]);
        let cache = Ms1Cache::store(&index).unwrap();
        cache.file.as_file().set_len(25).unwrap();
        assert!(cache.load().is_err());
    }

    #[test]
    fn persistent_round_trip_is_lossless() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("run.cache");
        let index = Ms1RunIndex::new(vec![
            Ms1Scan {
                rt: 3.123456789,
                peaks: vec![(600.123456789, 1234.5678)],
                id: "scan=2".into(),
                faims_cv: Some(-45.0),
            },
            Ms1Scan {
                rt: 1.0,
                peaks: vec![],
                ..Ms1Scan::default()
            },
        ]);
        Ms1Cache::store_at(&index, &path).unwrap();
        let loaded = Ms1Cache::load_from(&path).unwrap();
        assert_eq!(loaded.len(), 2);
        for i in 0..2 {
            assert_eq!(loaded.rt(i).to_bits(), index.rt(i).to_bits());
            assert_eq!(loaded.peaks(i), index.peaks(i));
            assert_eq!(loaded.scan(i).faims_cv, index.scan(i).faims_cv);
        }
    }

    #[test]
    fn cache_rejects_wrong_magic_and_version() {
        let dir = tempfile::tempdir().unwrap();
        let bad = dir.path().join("bad.cache");
        std::fs::write(&bad, b"NOTACACHE").unwrap();
        let e = Ms1Cache::load_from(&bad).unwrap_err().to_string();
        assert!(
            e.contains("bad.cache") && e.contains("not an andes run cache"),
            "{e}"
        );
        let mut v2 = b"ANDESRC1".to_vec();
        v2.extend_from_slice(&2u16.to_le_bytes());
        v2.extend_from_slice(&0u16.to_le_bytes());
        std::fs::write(&bad, &v2).unwrap();
        let e = Ms1Cache::load_from(&bad).unwrap_err().to_string();
        assert!(e.contains("version 2"), "{e}");
    }

    #[test]
    fn truncated_persistent_cache_names_the_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("run.cache");
        let index = Ms1RunIndex::new(vec![Ms1Scan {
            rt: 1.0,
            peaks: vec![(500.0, 1.0); 100],
            ..Ms1Scan::default()
        }]);
        Ms1Cache::store_at(&index, &path).unwrap();
        let bytes = std::fs::read(&path).unwrap();
        std::fs::write(&path, &bytes[..bytes.len() - 6]).unwrap();
        let err = Ms1Cache::load_from(&path).unwrap_err();
        let msg = err.to_string();
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof, "{msg}");
        assert!(msg.contains("run.cache"), "{msg}");
    }
}
