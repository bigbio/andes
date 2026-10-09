//! Private, lossless scratch storage for a run's MS1 scans. No run stays
//! resident while Percolator trains; quantification loads one run at a time.
use std::io::{self, BufReader, BufWriter, Read, Write};

use model::scan::Ms1Scan;
use quant::ms1_index::Ms1RunIndex;
use tempfile::NamedTempFile;

pub(super) struct Ms1Cache {
    file: NamedTempFile,
}

impl Ms1Cache {
    pub(super) fn store(index: &Ms1RunIndex) -> io::Result<Self> {
        let file = tempfile::Builder::new().prefix("andes-ms1-").tempfile()?;
        let mut out = BufWriter::with_capacity(1024 * 1024, file.as_file());
        out.write_all(&(index.len() as u64).to_le_bytes())?;
        for i in 0..index.len() {
            out.write_all(&index.rt(i).to_le_bytes())?;
            out.write_all(&(index.peaks(i).len() as u64).to_le_bytes())?;
            for &(mz, intensity) in index.peaks(i) {
                out.write_all(&mz.to_le_bytes())?;
                out.write_all(&intensity.to_le_bytes())?;
            }
        }
        out.flush()?;
        drop(out);
        Ok(Self { file })
    }

    pub(super) fn load(&self) -> io::Result<Ms1RunIndex> {
        let file = self.file.reopen()?;
        let mut remaining = file.metadata()?.len();
        let mut input = BufReader::with_capacity(1024 * 1024, file);
        let count = u64::from_le_bytes(read::<8>(&mut input, &mut remaining)?);
        let count = bounded_count(count, remaining / 16)?;
        let mut scans = Vec::with_capacity(count);
        for _ in 0..count {
            let rt = f64::from_le_bytes(read::<8>(&mut input, &mut remaining)?);
            let n = u64::from_le_bytes(read::<8>(&mut input, &mut remaining)?);
            let n = bounded_count(n, remaining / 12)?;
            let mut peaks = Vec::with_capacity(n);
            for _ in 0..n {
                let mz = f64::from_le_bytes(read::<8>(&mut input, &mut remaining)?);
                let intensity = f32::from_le_bytes(read::<4>(&mut input, &mut remaining)?);
                peaks.push((mz, intensity));
            }
            scans.push(Ms1Scan { rt, peaks });
        }
        if remaining != 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "trailing MS1 cache data",
            ));
        }
        Ok(Ms1RunIndex::new(scans))
    }
}

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
            },
            Ms1Scan {
                rt: 1.0,
                peaks: vec![],
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
        }]);
        let cache = Ms1Cache::store(&index).unwrap();
        cache.file.as_file().set_len(25).unwrap();
        assert!(cache.load().is_err());
    }
}
