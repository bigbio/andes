use crate::budget::{Budget, Reservation};
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::os::unix::fs::FileExt;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

static TEXT_COUNTER: AtomicU64 = AtomicU64::new(0);

pub(crate) fn append_record(out: &mut Vec<u8>, psm_id: &str, peptide: &str, proteins: &str) {
    out.extend_from_slice(psm_id.as_bytes());
    out.push(b'\t');
    out.extend_from_slice(peptide.as_bytes());
    out.push(b'\t');
    out.extend_from_slice(proteins.as_bytes());
}

pub struct RowText {
    pub psm_id: String,
    pub peptide: String,
    pub proteins: String,
}

enum Data {
    Ram(Vec<u8>),
    Spill(PathBuf, File),
}

pub struct Texts {
    offsets: Vec<u64>,
    data: Data,
    budget: Arc<Budget>,
    _ram: Reservation,
}

/// Records of a spilled store at most this many bytes apart are read with one positioned
/// read, the bytes between them included.
const TEXT_GAP_BYTES: u64 = 16 << 10;

/// Bytes one positioned read of spilled records reads at most.
const TEXT_SPAN_BYTES: u64 = 4 << 20;

/// Records one task of `get_many` reads.
const TEXT_TASK_ROWS: usize = 1 << 14;

/// Bytes a returned record and the output line made from it hold beyond twice its own length.
const TEXT_ROW_OVERHEAD_BYTES: u64 = 128;

impl Texts {
    /// Bytes of every record together.
    pub fn bytes(&self) -> u64 {
        self.offsets.last().copied().unwrap_or(0)
    }

    pub fn len(&self) -> usize {
        self.offsets.len().saturating_sub(1)
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub(crate) fn spill(&mut self, budget: &Budget) -> Result<(), String> {
        let Data::Ram(v) = &self.data else {
            return Ok(());
        };
        let k = TEXT_COUNTER.fetch_add(1, Ordering::Relaxed);
        let path = budget
            .spill_dir()
            .join(format!("rescore-{}-{}.txt", std::process::id(), k));
        let res = std::fs::write(&path, v).and_then(|()| File::open(&path));
        let file = match res {
            Ok(f) => f,
            Err(e) => {
                let _ = std::fs::remove_file(&path);
                return Err(format!("spill file {}: {e}", path.display()));
            }
        };
        self.data = Data::Spill(path, file);
        self._ram = Reservation::none();
        Ok(())
    }

    pub(crate) fn is_spilled(&self) -> bool {
        matches!(self.data, Data::Spill(..))
    }

    /// How many rows one `get_many` call may ask for: as many as half the free storage budget
    /// holds with their output lines, at least `min`, and the reservation that holds them.
    pub(crate) fn batch_rows(&self, min: usize) -> (usize, Reservation) {
        let per_row = 2 * self.bytes() / self.len().max(1) as u64 + TEXT_ROW_OVERHEAD_BYTES;
        let mut held = Reservation::none();
        held.grow_up_to(&self.budget, self.budget.available() / 2);
        let rows = (held.bytes() / per_row) as usize;
        (rows.max(min), held)
    }

    /// The texts of `rows`, in the order given. A spilled store is read in file order by
    /// several threads, each with positioned reads that take nearby records together.
    pub(crate) fn get_many(&self, rows: &[usize]) -> Result<Vec<RowText>, String> {
        let Data::Spill(path, file) = &self.data else {
            return rows.iter().map(|&i| self.get(i)).collect();
        };
        let mut order: Vec<usize> = (0..rows.len()).collect();
        order.par_sort_unstable_by_key(|&k| rows[k]);
        let parts: Vec<Result<Vec<RowText>, String>> = order
            .par_chunks(TEXT_TASK_ROWS)
            .map(|ks| {
                self.read_spans(file, rows, ks)
                    .map_err(|e| format!("read {}: {e}", path.display()))
            })
            .collect();
        let mut out: Vec<Option<RowText>> = (0..rows.len()).map(|_| None).collect();
        for (ks, part) in order.chunks(TEXT_TASK_ROWS).zip(parts) {
            for (&k, t) in ks.iter().zip(part?) {
                out[k] = Some(t);
            }
        }
        Ok(out
            .into_iter()
            .map(|t| t.expect("every row read"))
            .collect())
    }

    /// The records of `rows[k]` for each `k` of `ks`, whose rows must not decrease.
    fn read_spans(
        &self,
        file: &File,
        rows: &[usize],
        ks: &[usize],
    ) -> std::io::Result<Vec<RowText>> {
        let span = |k: usize| (self.offsets[rows[k]], self.offsets[rows[k] + 1]);
        let mut out = Vec::with_capacity(ks.len());
        let mut buf = Vec::new();
        let mut a = 0usize;
        while a < ks.len() {
            let (lo, mut hi) = span(ks[a]);
            let mut b = a + 1;
            while b < ks.len() {
                let (s, e) = span(ks[b]);
                if s > hi + TEXT_GAP_BYTES || e - lo > TEXT_SPAN_BYTES {
                    break;
                }
                hi = hi.max(e);
                b += 1;
            }
            buf.resize((hi - lo) as usize, 0);
            file.read_exact_at(&mut buf, lo)?;
            for &k in &ks[a..b] {
                let (s, e) = span(k);
                out.push(split_record(&buf[(s - lo) as usize..(e - lo) as usize]));
            }
            a = b;
        }
        Ok(out)
    }

    pub fn get(&self, i: usize) -> Result<RowText, String> {
        let (a, b) = (self.offsets[i], self.offsets[i + 1]);
        let bytes = match &self.data {
            Data::Ram(v) => v[a as usize..b as usize].to_vec(),
            Data::Spill(path, file) => {
                let mut buf = vec![0u8; (b - a) as usize];
                file.read_exact_at(&mut buf, a)
                    .map_err(|e| format!("read {}: {e}", path.display()))?;
                buf
            }
        };
        Ok(split_record(&bytes))
    }
}

fn split_record(bytes: &[u8]) -> RowText {
    let text = String::from_utf8_lossy(bytes);
    let mut parts = text.splitn(3, '\t');
    RowText {
        psm_id: parts.next().unwrap_or("").to_string(),
        peptide: parts.next().unwrap_or("").to_string(),
        proteins: parts.next().unwrap_or("").to_string(),
    }
}

impl Drop for Texts {
    fn drop(&mut self) {
        if let Data::Spill(p, _) = &self.data {
            let _ = std::fs::remove_file(p);
        }
    }
}

enum Sink {
    Ram(Vec<u8>),
    Spill(PathBuf, BufWriter<File>),
}

pub(crate) struct TextsBuilder {
    offsets: Vec<u64>,
    len: u64,
    budget: Arc<Budget>,
    ram: Reservation,
    sink: Sink,
}

impl TextsBuilder {
    pub(crate) fn new(budget: &Arc<Budget>) -> Self {
        TextsBuilder {
            offsets: vec![0],
            len: 0,
            budget: budget.clone(),
            ram: Reservation::none(),
            sink: Sink::Ram(Vec::new()),
        }
    }

    pub(crate) fn push_record(&mut self, record: &[u8]) -> Result<(), String> {
        let n = record.len() as u64;
        if let Sink::Ram(v) = &mut self.sink {
            if self.ram.grow(&self.budget, n, self.budget.total() / 2) {
                v.extend_from_slice(record);
            } else {
                let k = TEXT_COUNTER.fetch_add(1, Ordering::Relaxed);
                let path = self.budget.spill_dir().join(format!(
                    "rescore-{}-{}.txt",
                    std::process::id(),
                    k
                ));
                let file = File::create(&path)
                    .map_err(|e| format!("create spill file {}: {e}", path.display()))?;
                let mut w = BufWriter::with_capacity(1 << 20, file);
                let res = w.write_all(v);
                self.sink = Sink::Spill(path.clone(), w);
                self.ram = Reservation::none();
                res.map_err(|e| format!("write spill file {}: {e}", path.display()))?;
            }
        }
        if let Sink::Spill(path, w) = &mut self.sink {
            w.write_all(record)
                .map_err(|e| format!("write spill file {}: {e}", path.display()))?;
        }
        self.len += n;
        self.offsets.push(self.len);
        Ok(())
    }

    pub(crate) fn finish(mut self) -> Result<Texts, String> {
        let sink = std::mem::replace(&mut self.sink, Sink::Ram(Vec::new()));
        let offsets = std::mem::take(&mut self.offsets);
        let ram = std::mem::replace(&mut self.ram, Reservation::none());
        let data = match sink {
            Sink::Ram(v) => Data::Ram(v),
            Sink::Spill(path, mut w) => {
                let flushed = w.flush();
                drop(w);
                let file = match flushed.and_then(|()| File::open(&path)) {
                    Ok(f) => f,
                    Err(e) => {
                        let _ = std::fs::remove_file(&path);
                        return Err(format!("spill file {}: {e}", path.display()));
                    }
                };
                Data::Spill(path, file)
            }
        };
        Ok(Texts {
            offsets,
            data,
            budget: self.budget.clone(),
            _ram: ram,
        })
    }
}

impl Drop for TextsBuilder {
    fn drop(&mut self) {
        if let Sink::Spill(p, _) = &self.sink {
            let _ = std::fs::remove_file(p);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn texts(budget_bytes: u64, dir: &std::path::Path, n: usize) -> Texts {
        let budget = Budget::new(budget_bytes, dir);
        let mut b = TextsBuilder::new(&budget);
        for i in 0..n {
            let mut rec = Vec::new();
            let pad = "P".repeat(i % 37);
            append_record(&mut rec, &format!("psm{i}"), &format!("K.{pad}R.A"), "p1\tp2");
            b.push_record(&rec).unwrap();
        }
        b.finish().unwrap()
    }

    #[test]
    fn spilled_batch_reads_equal_single_reads_in_any_order() {
        let dir = tempfile::tempdir().unwrap();
        let n = 50_000;
        let ram = texts(1 << 30, dir.path(), n);
        let spilled = texts(0, dir.path(), n);
        assert!(!ram.is_spilled() && spilled.is_spilled());
        let mut rows: Vec<usize> = (0..n).map(|i| (i * 7919) % n).collect();
        rows.extend([3, 3, n - 1, 0, 40_000, 3]);
        let sparse: Vec<usize> = (0..n).rev().step_by(997).collect();
        assert_eq!(spilled.get_many(&sparse).unwrap()[1].psm_id, format!("psm{}", n - 998));
        let got = spilled.get_many(&rows).unwrap();
        assert_eq!(got.len(), rows.len());
        for (&i, t) in rows.iter().zip(&got) {
            let want = ram.get(i).unwrap();
            assert_eq!(
                (&t.psm_id, &t.peptide, &t.proteins),
                (&want.psm_id, &want.peptide, &want.proteins)
            );
            assert_eq!(spilled.get(i).unwrap().psm_id, want.psm_id);
        }
        let (k, held) = spilled.batch_rows(10);
        assert_eq!((k, held.bytes()), (10, 0));
    }
}
