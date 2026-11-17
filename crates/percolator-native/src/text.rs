use crate::budget::{Budget, Reservation};
use std::fs::File;
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

static TEXT_COUNTER: AtomicU64 = AtomicU64::new(0);

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
        let text = String::from_utf8_lossy(&bytes);
        let mut parts = text.splitn(3, '\t');
        Ok(RowText {
            psm_id: parts.next().unwrap_or("").to_string(),
            peptide: parts.next().unwrap_or("").to_string(),
            proteins: parts.next().unwrap_or("").to_string(),
        })
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

    pub(crate) fn push(
        &mut self,
        psm_id: &str,
        peptide: &str,
        proteins: &str,
    ) -> Result<(), String> {
        let record = [
            psm_id.as_bytes(),
            b"\t",
            peptide.as_bytes(),
            b"\t",
            proteins.as_bytes(),
        ];
        let n: u64 = record.iter().map(|p| p.len() as u64).sum();
        if let Sink::Ram(v) = &mut self.sink {
            if self.ram.grow(&self.budget, n, self.budget.total() / 2) {
                for p in record {
                    v.extend_from_slice(p);
                }
            } else {
                let k = TEXT_COUNTER.fetch_add(1, Ordering::Relaxed);
                let path = self.budget.spill_dir().join(format!(
                    "percolator-native-{}-{}.txt",
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
            for p in record {
                w.write_all(p)
                    .map_err(|e| format!("write spill file {}: {e}", path.display()))?;
            }
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
