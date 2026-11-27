use crate::budget::{Budget, Reservation};
use crate::cv::TASK_BYTES_PER_EXAMPLE;
use crate::store::{Store, StoreBuilder};
use crate::svm::{Rows, GRAM_BYTES};
use crate::text::{append_record, RowText, Texts, TextsBuilder};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

/// How the engine's storage budget is chosen.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Sizing {
    /// The budget is `MemLimit::bytes`, as given.
    Fixed,
    /// `MemLimit::bytes` is the memory the job or machine allows; the budget is what the PIN
    /// needs to stay in memory, never above `cap` nor above the allowance. `why` says where the
    /// allowance comes from.
    FitInput { cap: u64, why: String },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MemLimit {
    pub bytes: u64,
    pub spill_dir: PathBuf,
    pub sizing: Sizing,
}

/// Per-PSM metadata that stays in memory for every row. Packed (25 bytes, not 32 with padding):
/// on 500 M-row studies the padding alone is 3.5 GB of the storage budget. Fields are read by
/// value only.
#[repr(C, packed)]
pub struct RowMeta {
    pub label: i8,
    pub file: u32,
    pub scan: u64,
    pub exp_mass: f64,
    /// Hash of the row's peptide sequence without flanking residues or modifications, with
    /// I and L treated alike; rows that share it share cross-validation folds.
    pub sequence: u32,
}

/// FNV-1a hash of the residue letters of a PIN `Peptide` cell: flanking residues
/// (`K.PEPTIDE.R`), bracketed or parenthesised modification text and everything that is not
/// an upper-case residue letter (terminus markers such as `n`, mass deltas) are dropped, and
/// I is read as L, so every peptidoform of one sequence gets the same key.
pub fn sequence_key(peptide: &str) -> u32 {
    let p = peptide.trim().as_bytes();
    let core = if p.len() > 4 && p[1] == b'.' && p[p.len() - 2] == b'.' {
        &p[2..p.len() - 2]
    } else {
        p
    };
    let mut depth = 0u32;
    let mut h = 0x811c_9dc5u32;
    for &c in core {
        match c {
            b'[' | b'(' => depth += 1,
            b']' | b')' => depth = depth.saturating_sub(1),
            _ if depth == 0 && c.is_ascii_uppercase() => {
                let c = if c == b'I' { b'L' } else { c };
                h = (h ^ u32::from(c)).wrapping_mul(0x0100_0193);
            }
            _ => {}
        }
    }
    h
}

pub struct Pin {
    pub feature_names: Vec<String>,
    pub file_names: Vec<String>,
    pub default_direction: Option<Vec<f64>>,
    pub rows: Store,
    pub meta: Vec<RowMeta>,
    pub texts: Texts,
    pub budget: Arc<Budget>,
    _per_row: Reservation,
}

impl Pin {
    pub fn text(&self, i: usize) -> Result<RowText, String> {
        self.texts.get(i)
    }
}

#[derive(Default)]
struct Cols {
    spec_id: Option<usize>,
    label: Option<usize>,
    scan: Option<usize>,
    exp_mass: Option<usize>,
    file: Option<usize>,
    peptide: Option<usize>,
    proteins: Option<usize>,
    features: Vec<usize>,
}

fn parse_header(cells: &[&str]) -> Result<(Cols, Vec<String>), String> {
    let mut c = Cols::default();
    let mut names = Vec::new();
    for (i, cell) in cells.iter().enumerate() {
        match cell.trim().to_ascii_lowercase().as_str() {
            "specid" => c.spec_id = Some(i),
            "label" => c.label = Some(i),
            "scannr" => c.scan = Some(i),
            "expmass" => c.exp_mass = Some(i),
            "calcmass" => {}
            "filename" => c.file = Some(i),
            "peptide" => c.peptide = Some(i),
            "proteins" | "proteinid" | "proteinids" => c.proteins = Some(i),
            _ => {
                c.features.push(i);
                names.push(cell.trim().to_string());
            }
        }
    }
    if c.label.is_none() {
        return Err("PIN header has no Label column".into());
    }
    if c.peptide.is_none() {
        return Err("PIN header has no Peptide column".into());
    }
    if c.features.is_empty() {
        return Err("PIN header has no feature columns".into());
    }
    Ok((c, names))
}

struct Parsed {
    feats: Vec<f32>,
    meta: Vec<RowMeta>,
    file_names: Vec<String>,
    texts: Vec<u8>,
    text_ends: Vec<usize>,
    default_direction: Option<Vec<f64>>,
}

struct Layout<'a> {
    cols: &'a Cols,
    hcells: &'a [String],
    path: &'a Path,
}

fn cell<'a>(cells: &[&'a str], i: usize) -> &'a str {
    cells.get(i).map(|s| s.trim()).unwrap_or("")
}

fn parse_chunk(
    text: &[u8],
    first_lineno: usize,
    first_chunk: bool,
    l: &Layout,
) -> Result<Parsed, String> {
    let cols = l.cols;
    let nfeat = cols.features.len();
    let peptide_col = cols.peptide.unwrap_or(0);
    let label_col = cols.label.unwrap_or(0);
    let mut out = Parsed {
        feats: Vec::new(),
        meta: Vec::new(),
        file_names: Vec::new(),
        texts: Vec::new(),
        text_ends: Vec::new(),
        default_direction: None,
    };
    let mut local_files: HashMap<&str, u32> = HashMap::new();
    let mut row = vec![0f32; nfeat];
    let mut proteins = String::new();
    for (k, raw) in text.split(|&b| b == b'\n').enumerate() {
        let lineno = first_lineno + k;
        let line = std::str::from_utf8(raw).map_err(|_| {
            format!(
                "read {}: stream did not contain valid UTF-8",
                l.path.display()
            )
        })?;
        let line = line.trim_end_matches(['\r', '\n']);
        if line.is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        if first_chunk
            && k == 0
            && cells
                .first()
                .is_some_and(|c| c.trim().eq_ignore_ascii_case("defaultdirection"))
        {
            let mut w = Vec::with_capacity(nfeat);
            for &i in &cols.features {
                let c = cell(&cells, i);
                w.push(if c.is_empty() || c == "-" {
                    0.0
                } else {
                    c.parse::<f64>()
                        .map_err(|_| format!("line {lineno}: bad DefaultDirection value '{c}'"))?
                });
            }
            out.default_direction = Some(w);
            continue;
        }
        if cells.len() <= peptide_col {
            return Err(format!(
                "line {lineno}: {} cells, Peptide column is {}",
                cells.len(),
                peptide_col + 1
            ));
        }
        let Some(label_cell) = cells.get(label_col) else {
            return Err(format!(
                "line {lineno}: {} cells, Label column is {}",
                cells.len(),
                label_col + 1
            ));
        };
        let label = match label_cell.trim() {
            "1" | "+1" => 1i8,
            "-1" => -1i8,
            other => return Err(format!("line {lineno}: bad Label '{other}'")),
        };
        for (slot, &i) in row.iter_mut().zip(&cols.features) {
            let c = cell(&cells, i);
            let v = c
                .parse::<f64>()
                .map_err(|_| format!("line {lineno}: bad value '{c}' in column {}", i + 1))?;
            if !v.is_finite() || !(v as f32).is_finite() {
                return Err(format!(
                    "line {lineno}: non-finite feature value '{c}' in column {} ({})",
                    i + 1,
                    l.hcells[i].trim()
                ));
            }
            *slot = v as f32;
        }
        let file_idx = match cols.file {
            Some(i) => {
                let name = cells.get(i).copied().unwrap_or("");
                match local_files.get(name) {
                    Some(&p) => p,
                    None => {
                        let id = out.file_names.len() as u32;
                        local_files.insert(name, id);
                        out.file_names.push(name.to_string());
                        id
                    }
                }
            }
            None => 0,
        };
        let scan = match cols.scan {
            Some(i) => {
                let c = cell(&cells, i);
                c.parse::<u64>()
                    .map_err(|_| format!("line {lineno}: bad ScanNr '{c}'"))?
            }
            None => 0,
        };
        let exp_mass = match cols.exp_mass {
            Some(i) => {
                let c = cell(&cells, i);
                c.parse::<f64>()
                    .map_err(|_| format!("line {lineno}: bad ExpMass '{c}'"))?
            }
            None => 0.0,
        };
        proteins.clear();
        if let Some(i) = cols.proteins {
            if i < cells.len() {
                for (j, c) in cells[i..].iter().enumerate() {
                    if j > 0 {
                        proteins.push('\t');
                    }
                    proteins.push_str(c);
                }
            }
        }
        let psm_id = cols
            .spec_id
            .and_then(|i| cells.get(i))
            .copied()
            .unwrap_or("");
        append_record(&mut out.texts, psm_id, cells[peptide_col], &proteins);
        out.text_ends.push(out.texts.len());
        out.meta.push(RowMeta {
            label,
            file: file_idx,
            scan,
            exp_mass,
            sequence: sequence_key(cells[peptide_col]),
        });
        out.feats.extend_from_slice(&row);
    }
    Ok(out)
}

const CHUNK_BYTES: usize = 4 << 20;

/// Bytes each PSM row holds in memory from the read to the end of the run, whatever the
/// budget: its `RowMeta` and its text offset. They are charged to the storage budget as the
/// rows are read, so the rest of the budget is memory the stores, training rows and fits can
/// really use.
const PER_ROW_BYTES: u64 = (std::mem::size_of::<RowMeta>() + std::mem::size_of::<u32>()) as u64;

const IN_FLIGHT_BYTES_PER_CHUNK: u64 = 4 * CHUNK_BYTES as u64;

fn chunks_in_flight(available: u64, threads: usize) -> usize {
    ((available / 4 / IN_FLIGHT_BYTES_PER_CHUNK) as usize).clamp(1, threads.max(1) * 2)
}

fn read_chunk(r: &mut impl BufRead, path: &Path) -> Result<Vec<u8>, String> {
    let mut buf = Vec::with_capacity(CHUNK_BYTES + 4096);
    let err = |e: std::io::Error| format!("read {}: {e}", path.display());
    while buf.len() < CHUNK_BYTES {
        let avail = r.fill_buf().map_err(err)?;
        if avail.is_empty() {
            return Ok(buf);
        }
        let take = avail.len().min(CHUNK_BYTES - buf.len());
        buf.extend_from_slice(&avail[..take]);
        r.consume(take);
    }
    r.read_until(b'\n', &mut buf).map_err(err)?;
    if buf.last() == Some(&b'\n') {
        buf.pop();
    }
    Ok(buf)
}

type Batch = Vec<(Vec<u8>, usize, bool)>;

fn read_batch(
    r: &mut impl BufRead,
    path: &Path,
    batch: usize,
    lineno: &mut usize,
    first: &mut bool,
) -> Result<Batch, String> {
    let mut chunks: Batch = Vec::with_capacity(batch);
    while chunks.len() < batch {
        let c = read_chunk(r, path)?;
        if c.is_empty() {
            break;
        }
        let lines = c.iter().filter(|&&b| b == b'\n').count() + 1;
        chunks.push((c, *lineno, *first));
        *lineno += lines;
        *first = false;
    }
    Ok(chunks)
}

pub fn read_pin(path: &Path, limit: MemLimit) -> Result<Pin, String> {
    let file = File::open(path).map_err(|e| format!("open {}: {e}", path.display()))?;
    let mut reader = BufReader::with_capacity(1 << 20, file);
    let mut header = String::new();
    let got = reader
        .read_line(&mut header)
        .map_err(|e| format!("read {}: {e}", path.display()))?;
    if got == 0 {
        return Err(format!("{} is empty", path.display()));
    }
    let hcells: Vec<String> = header
        .trim_end_matches(['\r', '\n'])
        .split('\t')
        .map(str::to_string)
        .collect();
    let hrefs: Vec<&str> = hcells.iter().map(String::as_str).collect();
    let (cols, feature_names) = parse_header(&hrefs)?;
    let nfeat = cols.features.len();
    let layout = Layout {
        cols: &cols,
        hcells: &hcells,
        path,
    };
    let reading_budget = match &limit.sizing {
        Sizing::Fixed => limit.bytes,
        Sizing::FitInput { cap, .. } => limit.bytes.min(*cap),
    };
    let budget = Budget::new(reading_budget, &limit.spill_dir);
    let mut builder = StoreBuilder::new(nfeat, &budget);
    let mut texts = TextsBuilder::new(&budget);
    let mut meta: Vec<RowMeta> = Vec::new();
    let mut file_ids: HashMap<String, u32> = HashMap::new();
    let mut file_names: Vec<String> = Vec::new();
    let mut default_direction = None;
    let mut per_row = Reservation::none();
    let batch = chunks_in_flight(budget.available(), rayon::current_num_threads());
    let mut in_flight = Reservation::none();
    in_flight.grow_up_to(&budget, batch as u64 * IN_FLIGHT_BYTES_PER_CHUNK);
    let mut lineno = 2usize;
    let mut first = true;
    let mut chunks = read_batch(&mut reader, path, batch, &mut lineno, &mut first)?;
    while !chunks.is_empty() {
        let (next, merged) = rayon::join(
            || read_batch(&mut reader, path, batch, &mut lineno, &mut first),
            || -> Result<(), String> {
                let parsed: Vec<Result<Parsed, String>> = chunks
                    .par_iter()
                    .map(|(c, l, f)| parse_chunk(c, *l, *f, &layout))
                    .collect();
                for p in parsed {
                    let p = p?;
                    if p.default_direction.is_some() {
                        default_direction = p.default_direction;
                    }
                    let map: Vec<u32> = p
                        .file_names
                        .into_iter()
                        .map(|name| match file_ids.get(&name) {
                            Some(&g) => g,
                            None => {
                                let g = file_ids.len() as u32;
                                file_ids.insert(name.clone(), g);
                                file_names.push(name);
                                g
                            }
                        })
                        .collect();
                    let mut start = 0usize;
                    let rows = p.meta.len();
                    per_row.grow_up_to(&budget, rows as u64 * PER_ROW_BYTES);
                    for (i, mut m) in p.meta.into_iter().enumerate() {
                        if cols.file.is_some() {
                            m.file = map[m.file as usize];
                        }
                        texts.push_record(&p.texts[start..p.text_ends[i]])?;
                        start = p.text_ends[i];
                        meta.push(m);
                    }
                    builder.push_rows(rows, &p.feats)?;
                }
                Ok(())
            },
        );
        merged?;
        chunks = next?;
    }
    drop(in_flight);
    log::debug!(
        target: "rescore",
        "after reading the PIN: {:.2} GiB of the {:.2} GiB storage budget free, {:.2} GiB charged to {} rows",
        budget.available() as f64 / (1u64 << 30) as f64,
        budget.total() as f64 / (1u64 << 30) as f64,
        per_row.bytes() as f64 / (1u64 << 30) as f64,
        meta.len()
    );
    let pin = Pin {
        feature_names,
        file_names,
        default_direction,
        rows: builder.finish()?,
        meta,
        texts: texts.finish()?,
        budget,
        _per_row: per_row,
    };
    if let Sizing::FitInput { cap, why } = &limit.sizing {
        fit_budget_to_input(&pin, limit.bytes, *cap, why, &limit.spill_dir);
    }
    Ok(pin)
}

/// SVM fits the automatic budget leaves room for at once. Each fit is itself parallel, so
/// more at once mostly adds memory.
const AUTO_CONCURRENT_FITS: u64 = 4;

/// Bytes per PSM the engine allocates outside the storage budget (fold row lists, label and
/// score vectors); the automatic budget leaves them room in the allowance.
const UNBUDGETED_BYTES_PER_ROW: u64 = 96;

/// The storage budget at which the PIN's stores, the normalised training rows and
/// `AUTO_CONCURRENT_FITS` SVM fits all stay in memory: the stores take at most half of the
/// budget while the PIN is read, and the training rows, the fit reserve and the fits come on
/// top of the stores afterwards.
pub(crate) fn in_memory_need(rows: usize, nfeat: usize, text_bytes: u64, threads: usize) -> u64 {
    let n = rows as u64;
    let stores = n * PER_ROW_BYTES + n * nfeat as u64 * 4 + text_bytes;
    let reading = chunks_in_flight(u64::MAX, threads) as u64 * IN_FLIGHT_BYTES_PER_CHUNK;
    let examples = n * 2 / 3;
    let fit = examples * TASK_BYTES_PER_EXAMPLE + GRAM_BYTES as u64;
    let training = stores
        + Rows::bytes_for(rows, nfeat)
        + examples * TASK_BYTES_PER_EXAMPLE
        + AUTO_CONCURRENT_FITS * fit;
    (2 * (stores + reading)).max(training)
}

fn thousands(n: usize) -> String {
    let digits = n.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(c);
    }
    out
}

static BUDGET_LOGGED: AtomicBool = AtomicBool::new(false);

fn fit_budget_to_input(pin: &Pin, allowed: u64, cap: u64, why: &str, spill_dir: &Path) {
    let n = pin.meta.len();
    let nfeat = pin.rows.nfeat();
    let need = in_memory_need(n, nfeat, pin.texts.bytes(), rayon::current_num_threads());
    let room = allowed.saturating_sub(n as u64 * UNBUDGETED_BYTES_PER_ROW);
    pin.budget.shrink_to(need.min(cap).min(room));
    let budget = pin.budget.total();
    let gb = |b: u64| b as f64 / (1u64 << 30) as f64;
    let spills = need > budget;
    let first = !BUDGET_LOGGED.swap(true, Ordering::Relaxed);
    let level = if spills || first {
        log::Level::Info
    } else {
        log::Level::Debug
    };
    let psms = thousands(n);
    if !spills {
        log::log!(
            target: "rescore",
            level,
            "storage budget {:.1} GB, what {psms} PSMs x {nfeat} features need to stay in memory ({why})",
            gb(budget)
        );
    } else if cap < room {
        log::log!(
            target: "rescore",
            level,
            "storage budget {:.1} GB, the default cap (--rescore-memory-gb sets another): {psms} PSMs x {nfeat} features need {:.1} GB, so part of them spills to {}",
            gb(budget),
            gb(need),
            spill_dir.display()
        );
    } else {
        log::log!(
            target: "rescore",
            level,
            "storage budget {:.1} GB, all the memory allowed ({why}): {psms} PSMs x {nfeat} features need {:.1} GB, so part of them spills to {}",
            gb(budget),
            gb(need),
            spill_dir.display()
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_are_grouped_by_thousands() {
        assert_eq!(thousands(0), "0");
        assert_eq!(thousands(999), "999");
        assert_eq!(thousands(1000), "1,000");
        assert_eq!(thousands(61_234_567), "61,234,567");
    }

    #[test]
    fn the_in_memory_need_covers_twice_the_stores_and_grows_with_the_pin() {
        let small = in_memory_need(1_000, 20, 50_000, 16);
        let large = in_memory_need(1_000_000, 20, 50_000_000, 16);
        assert!(large > small);
        let stores = 1_000_000 * (PER_ROW_BYTES + 80) + 50_000_000;
        assert!(large >= 2 * stores);
    }

    #[test]
    fn reading_charges_each_row_only_what_it_keeps_in_memory() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.pin");
        let mut text = String::from("SpecId\tLabel\tScanNr\tExpMass\tf1\tf2\tPeptide\tProteins\n");
        let rows = 20_000usize;
        for i in 0..rows {
            let label = if i % 2 == 0 { 1 } else { -1 };
            text.push_str(&format!(
                "s{i}\t{label}\t{i}\t{}.5\t{}\t0.5\tK.PEP{i}.R\tP{i}\n",
                1000 + i,
                i % 7
            ));
        }
        std::fs::write(&path, text).unwrap();
        let limit = MemLimit {
            bytes: 64 << 20,
            spill_dir: dir.path().to_path_buf(),
            sizing: Sizing::Fixed,
        };
        let pin = read_pin(&path, limit).unwrap();
        assert!(!pin.rows.is_spilled() && !pin.texts.is_spilled());
        let used = pin.budget.total() - pin.budget.available();
        let stores = (rows * 2 * 4) as u64 + pin.texts.bytes();
        let kept = (std::mem::size_of::<RowMeta>() + std::mem::size_of::<u32>()) as u64;
        assert_eq!(used, stores + rows as u64 * kept);
    }

    #[test]
    fn chunks_in_flight_follow_the_budget() {
        let chunk = IN_FLIGHT_BYTES_PER_CHUNK;
        assert_eq!(chunks_in_flight(0, 128), 1);
        assert_eq!(chunks_in_flight(1 << 30, 128), 16);
        assert_eq!(chunks_in_flight(16 << 30, 16), 32);
        assert_eq!(chunks_in_flight(u64::MAX, 4), 8);
        for threads in [1, 16, 128] {
            for gb in [1u64, 2, 8, 16, 64] {
                let n = chunks_in_flight(gb << 30, threads) as u64;
                assert!(n == 1 || n * chunk <= (gb << 30) / 4);
                assert!(n <= threads as u64 * 2);
            }
        }
    }
}
