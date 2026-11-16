use crate::store::{Store, StoreBuilder};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MemLimit {
    pub bytes: u64,
    pub spill_dir: PathBuf,
}

pub struct RowMeta {
    pub psm_id: String,
    pub label: i8,
    pub file: u32,
    pub scan: u64,
    pub exp_mass: f64,
    pub peptide: String,
    pub proteins: String,
}

pub struct Pin {
    pub feature_names: Vec<String>,
    pub file_names: Vec<String>,
    pub default_direction: Option<Vec<f64>>,
    pub rows: Store,
    pub meta: Vec<RowMeta>,
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

pub fn read_pin(path: &Path, limit: MemLimit) -> Result<Pin, String> {
    let file = File::open(path).map_err(|e| format!("open {}: {e}", path.display()))?;
    let mut lines = BufReader::with_capacity(1 << 20, file).lines();
    let header = lines
        .next()
        .ok_or_else(|| format!("{} is empty", path.display()))?
        .map_err(|e| format!("read {}: {e}", path.display()))?;
    let hcells: Vec<&str> = header.trim_end_matches(['\r', '\n']).split('\t').collect();
    let (cols, feature_names) = parse_header(&hcells)?;
    let nfeat = cols.features.len();
    let peptide_col = cols.peptide.unwrap_or(0);
    let label_col = cols.label.unwrap_or(0);
    let mut default_direction = None;
    let mut builder = StoreBuilder::new(nfeat, limit.bytes, &limit.spill_dir);
    let mut meta: Vec<RowMeta> = Vec::new();
    let mut file_ids: HashMap<String, u32> = HashMap::new();
    let mut file_names: Vec<String> = Vec::new();
    let mut row = vec![0f32; nfeat];
    for (n, line) in lines.enumerate() {
        let lineno = n + 2;
        let line = line.map_err(|e| format!("read {}: {e}", path.display()))?;
        let line = line.trim_end_matches(['\r', '\n']);
        if line.is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        if n == 0
            && cells
                .first()
                .is_some_and(|c| c.trim().eq_ignore_ascii_case("defaultdirection"))
        {
            let mut w = Vec::with_capacity(nfeat);
            for &i in &cols.features {
                let cell = cells.get(i).map(|s| s.trim()).unwrap_or("");
                w.push(if cell.is_empty() || cell == "-" {
                    0.0
                } else {
                    cell.parse::<f64>().map_err(|_| {
                        format!("line {lineno}: bad DefaultDirection value '{cell}'")
                    })?
                });
            }
            default_direction = Some(w);
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
            let cell = cells.get(i).map(|s| s.trim()).unwrap_or("");
            let v = cell
                .parse::<f64>()
                .map_err(|_| format!("line {lineno}: bad value '{cell}' in column {}", i + 1))?;
            if !v.is_finite() || !(v as f32).is_finite() {
                return Err(format!(
                    "line {lineno}: non-finite feature value '{cell}' in column {} ({})",
                    i + 1,
                    hcells[i].trim()
                ));
            }
            *slot = v as f32;
        }
        let file_idx = match cols.file {
            Some(i) => {
                let name = cells.get(i).copied().unwrap_or("");
                match file_ids.get(name) {
                    Some(&p) => p,
                    None => {
                        let id = file_ids.len() as u32;
                        file_ids.insert(name.to_string(), id);
                        file_names.push(name.to_string());
                        id
                    }
                }
            }
            None => 0,
        };
        let scan = match cols.scan {
            Some(i) => {
                let cell = cells.get(i).map(|s| s.trim()).unwrap_or("");
                cell.parse::<u64>()
                    .map_err(|_| format!("line {lineno}: bad ScanNr '{cell}'"))?
            }
            None => 0,
        };
        let exp_mass = match cols.exp_mass {
            Some(i) => {
                let cell = cells.get(i).map(|s| s.trim()).unwrap_or("");
                cell.parse::<f64>()
                    .map_err(|_| format!("line {lineno}: bad ExpMass '{cell}'"))?
            }
            None => 0.0,
        };
        let proteins = match cols.proteins {
            Some(i) if i < cells.len() => cells[i..].join("\t"),
            _ => String::new(),
        };
        meta.push(RowMeta {
            psm_id: cols
                .spec_id
                .and_then(|i| cells.get(i))
                .map(|s| s.to_string())
                .unwrap_or_default(),
            label,
            file: file_idx,
            scan,
            exp_mass,
            peptide: cells[peptide_col].to_string(),
            proteins,
        });
        builder.push_row(&row)?;
    }
    Ok(Pin {
        feature_names,
        file_names,
        default_direction,
        rows: builder.finish()?,
        meta,
    })
}
