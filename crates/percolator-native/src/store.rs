use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

const BLOCK_BYTES: usize = 8 << 20;

static SPILL_COUNTER: AtomicU64 = AtomicU64::new(0);

enum Data {
    Ram(Vec<f32>),
    Spill(PathBuf),
}

pub struct Store {
    nfeat: usize,
    nrows: usize,
    block_rows: usize,
    data: Data,
}

impl Store {
    pub fn len(&self) -> usize {
        self.nrows
    }

    pub fn is_empty(&self) -> bool {
        self.nrows == 0
    }

    pub fn nfeat(&self) -> usize {
        self.nfeat
    }

    pub fn for_each_block(&self, mut f: impl FnMut(usize, &[f32])) {
        if self.nrows == 0 {
            return;
        }
        let rows = self.block_rows;
        match &self.data {
            Data::Ram(v) => {
                let step = rows * self.nfeat.max(1);
                if self.nfeat == 0 {
                    f(0, &[]);
                    return;
                }
                for (i, chunk) in v.chunks(step).enumerate() {
                    f(i * rows, chunk);
                }
            }
            Data::Spill(path) => {
                let mut file = File::open(path).expect("open spill file");
                let mut start = 0usize;
                let mut bytes = vec![0u8; rows * self.nfeat * 4];
                let mut vals = vec![0f32; rows * self.nfeat];
                while start < self.nrows {
                    let n = rows.min(self.nrows - start);
                    let len = n * self.nfeat;
                    file.read_exact(&mut bytes[..len * 4])
                        .expect("read spill file");
                    for (o, c) in vals[..len].iter_mut().zip(bytes[..len * 4].chunks_exact(4)) {
                        *o = f32::from_le_bytes([c[0], c[1], c[2], c[3]]);
                    }
                    f(start, &vals[..len]);
                    start += n;
                }
            }
        }
    }
}

impl Drop for Store {
    fn drop(&mut self) {
        if let Data::Spill(p) = &self.data {
            let _ = std::fs::remove_file(p);
        }
    }
}

enum Sink {
    Ram(Vec<f32>),
    Spill(PathBuf, BufWriter<File>),
}

pub(crate) struct StoreBuilder {
    nfeat: usize,
    block_rows: usize,
    nrows: usize,
    limit_bytes: u64,
    spill_dir: PathBuf,
    sink: Sink,
}

impl StoreBuilder {
    pub(crate) fn new(nfeat: usize, limit_bytes: u64, spill_dir: &Path) -> Self {
        StoreBuilder {
            nfeat,
            block_rows: (BLOCK_BYTES / (4 * nfeat.max(1))).max(1),
            nrows: 0,
            limit_bytes,
            spill_dir: spill_dir.to_path_buf(),
            sink: Sink::Ram(Vec::new()),
        }
    }

    #[cfg(test)]
    pub(crate) fn with_block_rows(mut self, rows: usize) -> Self {
        self.block_rows = rows.max(1);
        self
    }

    pub(crate) fn push_row(&mut self, row: &[f32]) -> Result<(), String> {
        self.nrows += 1;
        match &mut self.sink {
            Sink::Ram(v) => {
                v.extend_from_slice(row);
                if (v.len() as u64) * 4 > self.limit_bytes {
                    let n = SPILL_COUNTER.fetch_add(1, Ordering::Relaxed);
                    let path = self.spill_dir.join(format!(
                        "percolator-native-{}-{}.f32",
                        std::process::id(),
                        n
                    ));
                    let file = File::create(&path)
                        .map_err(|e| format!("create spill file {}: {e}", path.display()))?;
                    let mut w = BufWriter::with_capacity(1 << 20, file);
                    let res = write_f32(&mut w, v);
                    self.sink = Sink::Spill(path.clone(), w);
                    res.map_err(|e| format!("write spill file {}: {e}", path.display()))?;
                }
            }
            Sink::Spill(path, w) => {
                write_f32(w, row)
                    .map_err(|e| format!("write spill file {}: {e}", path.display()))?;
            }
        }
        Ok(())
    }

    pub(crate) fn finish(mut self) -> Result<Store, String> {
        let sink = std::mem::replace(&mut self.sink, Sink::Ram(Vec::new()));
        let data = match sink {
            Sink::Ram(v) => Data::Ram(v),
            Sink::Spill(path, mut w) => {
                let flushed = w.flush();
                drop(w);
                let store = Store {
                    nfeat: self.nfeat,
                    nrows: self.nrows,
                    block_rows: self.block_rows,
                    data: Data::Spill(path.clone()),
                };
                flushed.map_err(|e| format!("flush spill file {}: {e}", path.display()))?;
                return Ok(store);
            }
        };
        Ok(Store {
            nfeat: self.nfeat,
            nrows: self.nrows,
            block_rows: self.block_rows,
            data,
        })
    }
}

impl Drop for StoreBuilder {
    fn drop(&mut self) {
        if let Sink::Spill(p, _) = &self.sink {
            let _ = std::fs::remove_file(p);
        }
    }
}

fn write_f32(w: &mut impl Write, v: &[f32]) -> std::io::Result<()> {
    for x in v {
        w.write_all(&x.to_le_bytes())?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build(limit: u64, dir: &Path, nrows: usize, nfeat: usize) -> Store {
        let mut b = StoreBuilder::new(nfeat, limit, dir).with_block_rows(3);
        for r in 0..nrows {
            let row: Vec<f32> = (0..nfeat).map(|c| (r * 100 + c) as f32 + 0.5).collect();
            b.push_row(&row).unwrap();
        }
        b.finish().unwrap()
    }

    fn collect(s: &Store) -> Vec<(usize, Vec<f32>)> {
        let mut v = Vec::new();
        s.for_each_block(|start, b| v.push((start, b.to_vec())));
        v
    }

    #[test]
    fn multi_block_spill_equals_ram_with_partial_last_block() {
        let ram_dir = tempfile::tempdir().unwrap();
        let ram = build(u64::MAX, ram_dir.path(), 10, 4);
        let spill_dir = tempfile::tempdir().unwrap();
        let spill = build(1, spill_dir.path(), 10, 4);
        let a = collect(&ram);
        let b = collect(&spill);
        let starts: Vec<usize> = a.iter().map(|x| x.0).collect();
        assert_eq!(starts, vec![0, 3, 6, 9]);
        let lens: Vec<usize> = a.iter().map(|x| x.1.len()).collect();
        assert_eq!(lens, vec![12, 12, 12, 4]);
        assert_eq!(a, b);
        assert_eq!(a[3].1, vec![900.5, 901.5, 902.5, 903.5]);
        assert_eq!(collect(&spill), b);
        drop(spill);
        assert_eq!(std::fs::read_dir(spill_dir.path()).unwrap().count(), 0);
    }
}
