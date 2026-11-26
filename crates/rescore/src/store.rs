use crate::budget::{Budget, Reservation};
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

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
    budget: Arc<Budget>,
    _ram: Reservation,
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

    pub fn is_spilled(&self) -> bool {
        matches!(self.data, Data::Spill(_))
    }

    pub fn block_rows(&self) -> usize {
        self.block_rows
    }

    pub fn par_for_each_block(&self, f: impl Fn(usize, &[f32]) + Sync) {
        if self.nrows == 0 {
            return;
        }
        match &self.data {
            Data::Ram(v) => {
                if self.nfeat == 0 {
                    f(0, &[]);
                    return;
                }
                let rows = self.block_rows;
                v.par_chunks(rows * self.nfeat)
                    .enumerate()
                    .for_each(|(i, c)| f(i * rows, c));
            }
            Data::Spill(path) => self.spilled_groups(path, |group| {
                group.par_iter().for_each(|(s, b)| f(*s, b));
            }),
        }
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
            Data::Spill(path) => self.spilled_groups(path, |group| {
                for (s, b) in group {
                    f(*s, b);
                }
            }),
        }
    }

    /// Reads a spilled store in groups of consecutive blocks, the blocks of a group read at
    /// once on the thread pool with positioned reads, and calls `each` with every group in
    /// file order. A group holds as many blocks as there are threads and the free budget holds
    /// with their read buffers, at least one.
    fn spilled_groups(&self, path: &Path, mut each: impl FnMut(&[(usize, Vec<f32>)])) {
        let rows = self.block_rows;
        let block = (rows * self.nfeat * 4) as u64;
        let threads = rayon::current_num_threads().max(1) as u64;
        let fits = (self.budget.available() / (2 * block).max(1)).min(threads);
        let held = Reservation::take(&self.budget, fits * 2 * block).filter(|_| fits >= 2);
        let group = if held.is_some() { fits as usize } else { 1 };
        let file = File::open(path).expect("open spill file");
        let nblocks = self.nrows.div_ceil(rows);
        let read = |b: usize| {
            let start = b * rows;
            let len = rows.min(self.nrows - start) * self.nfeat;
            let mut bytes = vec![0u8; len * 4];
            file.read_exact_at(&mut bytes, (start * self.nfeat * 4) as u64)
                .expect("read spill file");
            let vals: Vec<f32> = bytes
                .chunks_exact(4)
                .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]))
                .collect();
            (start, vals)
        };
        for first in (0..nblocks).step_by(group) {
            let last = (first + group).min(nblocks);
            let blocks: Vec<(usize, Vec<f32>)> = if group > 1 {
                (first..last).into_par_iter().map(read).collect()
            } else {
                (first..last).map(read).collect()
            };
            each(&blocks);
        }
        drop(held);
    }
}

/// A spilled store's columns, one file per feature, written in one pass over the store so
/// that a column can be read without reading every feature.
pub(crate) struct Columns {
    paths: Vec<PathBuf>,
    nrows: usize,
}

impl Columns {
    /// Calls `f` with every value of column `j`, in row order.
    pub(crate) fn for_each_value(&self, j: usize, mut f: impl FnMut(f32)) -> Result<(), String> {
        let path = &self.paths[j];
        let err = |e: std::io::Error| format!("read column file {}: {e}", path.display());
        let mut file = File::open(path).map_err(err)?;
        let mut bytes = vec![0u8; BLOCK_BYTES];
        let mut left = self.nrows * 4;
        while left > 0 {
            let len = left.min(BLOCK_BYTES);
            file.read_exact(&mut bytes[..len]).map_err(err)?;
            for c in bytes[..len].chunks_exact(4) {
                f(f32::from_le_bytes([c[0], c[1], c[2], c[3]]));
            }
            left -= len;
        }
        Ok(())
    }
}

impl Drop for Columns {
    fn drop(&mut self) {
        for p in &self.paths {
            let _ = std::fs::remove_file(p);
        }
    }
}

impl Store {
    /// Writes every feature column to its own file in the spill directory, in one pass.
    pub(crate) fn columns(&self) -> Result<Columns, String> {
        let nf = self.nfeat;
        let k = SPILL_COUNTER.fetch_add(1, Ordering::Relaxed);
        let mut cols = Columns {
            paths: Vec::with_capacity(nf),
            nrows: self.nrows,
        };
        let mut sinks: Vec<(BufWriter<File>, Vec<u8>)> = Vec::with_capacity(nf);
        for j in 0..nf {
            let path = self
                .budget
                .spill_dir()
                .join(format!("rescore-{}-{k}.col{j}", std::process::id()));
            let file = File::create(&path)
                .map_err(|e| format!("create column file {}: {e}", path.display()))?;
            cols.paths.push(path);
            sinks.push((BufWriter::with_capacity(1 << 20, file), Vec::new()));
        }
        let mut failed: Option<String> = None;
        self.for_each_block(|_, block| {
            if failed.is_some() || nf == 0 {
                return;
            }
            failed = sinks
                .par_iter_mut()
                .enumerate()
                .filter_map(|(j, (w, buf))| {
                    buf.clear();
                    for row in block.chunks_exact(nf) {
                        buf.extend_from_slice(&row[j].to_le_bytes());
                    }
                    w.write_all(buf).err().map(|e| (j, e))
                })
                .min_by_key(|(j, _)| *j)
                .map(|(j, e)| format!("write column file {}: {e}", cols.paths[j].display()));
        });
        for (j, (w, _)) in sinks.iter_mut().enumerate() {
            if failed.is_some() {
                break;
            }
            if let Err(e) = w.flush() {
                failed = Some(format!(
                    "write column file {}: {e}",
                    cols.paths[j].display()
                ));
            }
        }
        drop(sinks);
        match failed {
            Some(e) => Err(e),
            None => Ok(cols),
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
    budget: Arc<Budget>,
    ram: Reservation,
    sink: Sink,
}

impl StoreBuilder {
    pub(crate) fn new(nfeat: usize, budget: &Arc<Budget>) -> Self {
        StoreBuilder {
            nfeat,
            block_rows: (BLOCK_BYTES / (4 * nfeat.max(1))).max(1),
            nrows: 0,
            budget: budget.clone(),
            ram: Reservation::none(),
            sink: Sink::Ram(Vec::new()),
        }
    }

    #[cfg(test)]
    pub(crate) fn with_block_rows(mut self, rows: usize) -> Self {
        self.block_rows = rows.max(1);
        self
    }

    pub(crate) fn push_rows(&mut self, n: usize, data: &[f32]) -> Result<(), String> {
        let nf = self.nfeat;
        for i in 0..n {
            self.push_row(&data[i * nf..(i + 1) * nf])?;
        }
        Ok(())
    }

    pub(crate) fn push_row(&mut self, row: &[f32]) -> Result<(), String> {
        self.nrows += 1;
        match &mut self.sink {
            Sink::Ram(v) => {
                let keep = self.budget.total() / 2;
                if self.ram.grow(&self.budget, row.len() as u64 * 4, keep) {
                    v.extend_from_slice(row);
                } else {
                    let n = SPILL_COUNTER.fetch_add(1, Ordering::Relaxed);
                    let path = self.budget.spill_dir().join(format!(
                        "rescore-{}-{}.f32",
                        std::process::id(),
                        n
                    ));
                    let file = File::create(&path)
                        .map_err(|e| format!("create spill file {}: {e}", path.display()))?;
                    let mut w = BufWriter::with_capacity(1 << 20, file);
                    let res = write_f32(&mut w, v).and_then(|()| write_f32(&mut w, row));
                    self.sink = Sink::Spill(path.clone(), w);
                    self.ram = Reservation::none();
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
                    budget: self.budget.clone(),
                    _ram: Reservation::none(),
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
            budget: self.budget.clone(),
            _ram: std::mem::replace(&mut self.ram, Reservation::none()),
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
    use std::path::Path;

    fn build(limit: u64, dir: &Path, nrows: usize, nfeat: usize) -> Store {
        let mut b = StoreBuilder::new(nfeat, &Budget::new(limit, dir)).with_block_rows(3);
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

    #[test]
    fn spilled_blocks_read_in_groups_come_in_file_order() {
        let ram_dir = tempfile::tempdir().unwrap();
        let ram = build(u64::MAX, ram_dir.path(), 10, 4);
        let dir = tempfile::tempdir().unwrap();
        let spill = build(200, dir.path(), 10, 4);
        assert!(spill.is_spilled());
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap();
        assert_eq!(pool.install(|| collect(&spill)), collect(&ram));
        assert_eq!(spill.budget.available(), 200);
    }

    #[test]
    fn column_files_hold_each_feature_in_row_order() {
        let dir = tempfile::tempdir().unwrap();
        let spill = build(1, dir.path(), 10, 4);
        let cols = spill.columns().unwrap();
        for j in 0..4 {
            let mut got = Vec::new();
            cols.for_each_value(j, |x| got.push(x)).unwrap();
            let want: Vec<f32> = (0..10).map(|r| (r * 100 + j) as f32 + 0.5).collect();
            assert_eq!(got, want);
        }
        drop(cols);
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 1);
    }

    fn block_sums(s: &Store, threads: usize, live: &AtomicU64, peak: &AtomicU64) -> Vec<u64> {
        let br = s.block_rows();
        let nf = s.nfeat();
        let out: Vec<AtomicU64> = (0..s.len().div_ceil(br))
            .map(|_| AtomicU64::new(0))
            .collect();
        let block = (br * nf * 4) as u64;
        let budget = s.budget.clone();
        let spilled = s.is_spilled();
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| {
                s.par_for_each_block(|start, b| {
                    let now = live.fetch_add(1, Ordering::SeqCst) + 1;
                    peak.fetch_max(now, Ordering::SeqCst);
                    if spilled && now > 1 {
                        assert!(budget.total() - budget.available() >= (now + 2) * block);
                    }
                    std::thread::sleep(std::time::Duration::from_millis(2));
                    let mut acc = 0.0f64;
                    for row in b.chunks_exact(nf) {
                        for (j, &x) in row.iter().enumerate() {
                            acc += f64::from(x) * (j as f64 + 0.25);
                        }
                    }
                    out[start / br].store(acc.to_bits(), Ordering::SeqCst);
                    live.fetch_sub(1, Ordering::SeqCst);
                });
            });
        out.into_iter().map(AtomicU64::into_inner).collect()
    }

    #[test]
    fn spilled_parallel_pass_holds_only_the_blocks_the_budget_allows() {
        let (nrows, nfeat) = (100, 4);
        let block = 3 * nfeat as u64 * 4;
        let ram_dir = tempfile::tempdir().unwrap();
        let ram = build(u64::MAX, ram_dir.path(), nrows, nfeat);
        let reference = block_sums(&ram, 16, &AtomicU64::new(0), &AtomicU64::new(0));
        for (limit, most) in [
            (1, 1),
            (block + block / 2, 1),
            (3 * block + block / 2, 1),
            (5 * block + block / 2, 3),
        ] {
            let dir = tempfile::tempdir().unwrap();
            let spill = build(limit, dir.path(), nrows, nfeat);
            assert!(spill.is_spilled());
            let peak = AtomicU64::new(0);
            let got = block_sums(&spill, 16, &AtomicU64::new(0), &peak);
            assert_eq!(got, reference, "limit {limit}");
            let peak = peak.into_inner();
            assert!(peak <= most, "limit {limit}: {peak} blocks in flight");
            assert_eq!(spill.budget.available(), limit);
        }
    }
}
