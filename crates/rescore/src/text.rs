use crate::budget::{Budget, Reservation};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

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
    Spill(PathBuf, Mutex<File>),
}

pub struct Texts {
    offsets: Vec<u64>,
    data: Data,
    _ram: Reservation,
}

impl Texts {
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
        self.data = Data::Spill(path, Mutex::new(file));
        self._ram = Reservation::none();
        Ok(())
    }

    pub(crate) fn is_spilled(&self) -> bool {
        matches!(self.data, Data::Spill(..))
    }

    /// The texts of `rows`, in the order given. A spilled store is read in
    /// one locked pass in file order through a buffered reader, instead of a
    /// seek and a read per row.
    pub(crate) fn get_many(&self, rows: &[usize]) -> Result<Vec<RowText>, String> {
        let Data::Spill(path, file) = &self.data else {
            return rows.iter().map(|&i| self.get(i)).collect();
        };
        let err = |e: std::io::Error| format!("read {}: {e}", path.display());
        let mut order: Vec<usize> = (0..rows.len()).collect();
        order.sort_unstable_by_key(|&k| rows[k]);
        let mut out: Vec<Option<RowText>> = (0..rows.len()).map(|_| None).collect();
        let mut f = file
            .lock()
            .map_err(|_| "text spill file lock".to_string())?;
        let mut r = BufReader::with_capacity(1 << 20, &mut *f);
        let mut pos: Option<u64> = None;
        let mut buf = Vec::new();
        for &k in &order {
            let i = rows[k];
            let (a, b) = (self.offsets[i], self.offsets[i + 1]);
            match pos {
                Some(p) if a >= p => r.seek_relative((a - p) as i64).map_err(err)?,
                _ => {
                    r.seek(SeekFrom::Start(a)).map_err(err)?;
                }
            }
            buf.resize((b - a) as usize, 0);
            r.read_exact(&mut buf).map_err(err)?;
            pos = Some(b);
            out[k] = Some(split_record(&buf));
        }
        Ok(out
            .into_iter()
            .map(|t| t.expect("every row read"))
            .collect())
    }

    pub fn get(&self, i: usize) -> Result<RowText, String> {
        let (a, b) = (self.offsets[i], self.offsets[i + 1]);
        let bytes = match &self.data {
            Data::Ram(v) => v[a as usize..b as usize].to_vec(),
            Data::Spill(path, file) => {
                let mut buf = vec![0u8; (b - a) as usize];
                let mut f = file
                    .lock()
                    .map_err(|_| "text spill file lock".to_string())?;
                f.seek(SeekFrom::Start(a))
                    .and_then(|_| f.read_exact(&mut buf))
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
                Data::Spill(path, Mutex::new(file))
            }
        };
        Ok(Texts {
            offsets,
            data,
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
