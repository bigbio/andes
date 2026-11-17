use crate::budget::Budget;
use crate::store::{Store, StoreBuilder};
use crate::text::{append_record, RowText, Texts, TextsBuilder};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::Arc;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MemLimit {
    pub bytes: u64,
    pub spill_dir: PathBuf,
}

pub struct RowMeta {
    pub label: i8,
    pub file: u32,
    pub scan: u64,
    pub exp_mass: f64,
}

pub struct Pin {
    pub feature_names: Vec<String>,
    pub file_names: Vec<String>,
    pub default_direction: Option<Vec<f64>>,
    pub rows: Store,
    pub meta: Vec<RowMeta>,
    pub texts: Texts,
    pub budget: Arc<Budget>,
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
        let label = match cells[label_col].trim() {
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
        });
        out.feats.extend_from_slice(&row);
    }
    Ok(out)
}

const CHUNK_BYTES: usize = 4 << 20;

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
    let budget = Budget::new(limit.bytes, &limit.spill_dir);
    let mut builder = StoreBuilder::new(nfeat, &budget);
    let mut texts = TextsBuilder::new(&budget);
    let mut meta: Vec<RowMeta> = Vec::new();
    let mut file_ids: HashMap<String, u32> = HashMap::new();
    let mut file_names: Vec<String> = Vec::new();
    let mut default_direction = None;
    let batch = rayon::current_num_threads().max(1) * 2;
    let mut lineno = 2usize;
    let mut first = true;
    loop {
        let mut chunks: Vec<(Vec<u8>, usize, bool)> = Vec::with_capacity(batch);
        while chunks.len() < batch {
            let c = read_chunk(&mut reader, path)?;
            if c.is_empty() {
                break;
            }
            let lines = c.iter().filter(|&&b| b == b'\n').count() + 1;
            chunks.push((c, lineno, first));
            lineno += lines;
            first = false;
        }
        if chunks.is_empty() {
            break;
        }
        let parsed: Vec<Result<Parsed, String>> = chunks
            .par_iter()
            .map(|(c, l, f)| parse_chunk(c, *l, *f, &layout))
            .collect();
        drop(chunks);
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
    }
    Ok(Pin {
        feature_names,
        file_names,
        default_direction,
        rows: builder.finish()?,
        meta,
        texts: texts.finish()?,
        budget,
    })
}
