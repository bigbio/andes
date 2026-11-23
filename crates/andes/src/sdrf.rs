use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use model::activation::ActivationMethod;

use crate::param_sources::{ParamSources, Source};

#[derive(Debug, Clone, PartialEq)]
pub enum Quant {
    Lfq,
    Isobaric(String),
    Label(String),
}

#[derive(Debug, Clone, Default)]
pub struct SdrfSettings {
    pub path: PathBuf,
    pub data_files: Vec<PathBuf>,
    /// The file names the SDRF gives, in the order of `data_files` (which may be a converted
    /// copy found next to the SDRF instead).
    pub data_names: Vec<String>,
    pub enzyme: Option<String>,
    pub precursor_tol: Option<String>,
    pub fragment_tol: Option<String>,
    pub fixed_mods: Option<Vec<String>>,
    pub variable_mods: Option<Vec<String>>,
    pub quant: Option<Quant>,
    pub dissociation: Option<ActivationMethod>,
    pub xl: Vec<String>,
}

const DATA_FILE: &str = "comment[data file]";
const CLEAVAGE: &str = "comment[cleavage agent details]";
const MODS: &str = "comment[modification parameters]";
const PRECURSOR_TOL: &str = "comment[precursor mass tolerance]";
const FRAGMENT_TOL: &str = "comment[fragment mass tolerance]";
const LABEL: &str = "comment[label]";
const DISSOCIATION: &str = "comment[dissociation method]";
const XL: &str = "comment[cross-linker]";

pub fn read(path: &Path) -> Result<SdrfSettings, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("--sdrf: cannot read {}: {e}", path.display()))?;
    let dir = path.parent().unwrap_or_else(|| Path::new("."));
    parse(&text, dir)
        .map(|mut s| {
            s.path = path.to_path_buf();
            s
        })
        .map_err(|e| format!("--sdrf {}: {e}", path.display()))
}

pub fn parse(text: &str, dir: &Path) -> Result<SdrfSettings, String> {
    let (rows, has_mods) = parse_rows(text)?;
    settings_of(&rows, has_mods, dir)
}

/// The SDRF's rows split into groups that agree on every setting andes reads, one
/// [`SdrfSettings`] per group, in the order each group's first row appears. A file whose
/// rows all agree gives one group, the same settings [`parse`] gives.
pub fn split(text: &str, dir: &Path) -> Result<Vec<SdrfSettings>, String> {
    let mut groups: Vec<(Vec<String>, Vec<RowSettings>)> = Vec::new();
    let (rows, has_mods) = parse_rows(text)?;
    for row in rows {
        let key: Vec<String> = SETTINGS.iter().map(|(_, k)| k(&row)).collect();
        match groups.iter_mut().find(|(k, _)| *k == key) {
            Some((_, rows)) => rows.push(row),
            None => groups.push((key, vec![row])),
        }
    }
    groups
        .iter()
        .map(|(_, rows)| settings_of(rows, has_mods, dir))
        .collect()
}

type SettingKey = fn(&RowSettings) -> String;

/// The settings every row of one andes run must share: the SDRF column and the row's value.
const SETTINGS: [(&str, SettingKey); 6] = [
    (CLEAVAGE, |r| r.enzymes.join(",")),
    (MODS, |r| {
        let mut m: Vec<String> = r.mods.iter().map(|m| m.key()).collect();
        m.sort();
        m.join(", ")
    }),
    (PRECURSOR_TOL, |r| r.precursor_tol.clone().unwrap_or_default()),
    (DISSOCIATION, |r| r.dissociation.clone().unwrap_or_default()),
    (LABEL, |r| r.label.as_deref().map(label_kind).unwrap_or("").to_string()),
    (XL, |r| r.xl.clone().unwrap_or_default()),
];

/// The data rows, and whether the file has any modification column.
fn parse_rows(text: &str) -> Result<(Vec<RowSettings>, bool), String> {
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let header: Vec<String> = lines
        .next()
        .ok_or("empty file")?
        .split('\t')
        .map(norm_header)
        .collect();
    let cols = |name: &str| -> Vec<usize> {
        header
            .iter()
            .enumerate()
            .filter(|(_, h)| h.as_str() == name)
            .map(|(i, _)| i)
            .collect()
    };
    let data_col = *cols(DATA_FILE)
        .first()
        .ok_or_else(|| format!("no `{DATA_FILE}` column"))?;
    let (cleave, mods) = (cols(CLEAVAGE), cols(MODS));
    let (ptol, ftol) = (cols(PRECURSOR_TOL), cols(FRAGMENT_TOL));
    let (label_col, diss) = (cols(LABEL), cols(DISSOCIATION));
    let xl_col = cols(XL);

    let mut rows: Vec<RowSettings> = Vec::new();
    for (n, line) in lines.enumerate() {
        let cells: Vec<&str> = line.split('\t').collect();
        if cells.len() != header.len() {
            return Err(format!(
                "row {} has {} fields, the header {}",
                n + 2,
                cells.len(),
                header.len()
            ));
        }
        let cell = |i: usize| value(cells[i]);
        let at = |e: String| format!("row {}: {e}", n + 2);
        let file = cell(data_col).ok_or_else(|| format!("row {}: no data file", n + 2))?;
        let mut enzymes = Vec::new();
        for &i in &cleave {
            if let Some(v) = cell(i) {
                let e = enzyme_of(&v).map_err(at)?;
                if !enzymes.contains(&e) {
                    enzymes.push(e);
                }
            }
        }
        let mut row_mods = Vec::new();
        for &i in &mods {
            if let Some(v) = cell(i) {
                row_mods.push(ModParam::parse(&v).map_err(at)?);
            }
        }
        rows.push(RowSettings {
            file,
            enzymes,
            mods: row_mods,
            precursor_tol: ptol
                .first()
                .and_then(|&i| cell(i))
                .map(|v| tolerance(&v))
                .transpose()?,
            fragment_tol: ftol.first().and_then(|&i| cell(i)),
            label: label_col
                .first()
                .and_then(|&i| cell(i))
                .map(|v| nt(&v).map_err(|e| at(format!("`{LABEL}` {e}"))))
                .transpose()?,
            dissociation: diss
                .first()
                .and_then(|&i| cell(i))
                .map(|v| nt(&v).map_err(|e| at(format!("`{DISSOCIATION}` {e}"))))
                .transpose()?,
            xl: xl_col.first().and_then(|&i| cell(i)),
        });
    }
    if rows.is_empty() {
        return Err("no data rows".into());
    }
    Ok((rows, !mods.is_empty()))
}

fn settings_of(rows: &[RowSettings], has_mods: bool, dir: &Path) -> Result<SdrfSettings, String> {
    for (column, key) in &SETTINGS[..5] {
        agree(rows, column, key)?;
    }
    let mut linkers: Vec<String> = Vec::new();
    for v in rows.iter().filter_map(|r| r.xl.clone()) {
        if !linkers.contains(&v) {
            xl::spec::from_sdrf(&v, model::unimod::bundled())
                .map_err(|e| format!("`{XL}`: {e}"))?;
            linkers.push(v);
        }
    }
    if !linkers.is_empty() && rows.iter().any(|r| r.xl.is_none()) {
        return Err(format!(
            "`{XL}` is set for some files and not others; split the SDRF into crosslinked and \
             linear runs"
        ));
    }

    let first = &rows[0];
    let mut data_files: Vec<PathBuf> = Vec::new();
    let mut data_names: Vec<String> = Vec::new();
    for r in rows {
        let p = resolve_data_file(dir, &r.file);
        if !data_files.contains(&p) {
            data_files.push(p);
            data_names.push(r.file.clone());
        }
    }
    let labels: BTreeSet<String> = rows.iter().filter_map(|r| r.label.clone()).collect();
    let quant = quant_of(&labels, &first.mods)?;
    let searched = |m: &&ModParam| match &quant {
        Some(Quant::Label(_)) => !m.is_label_channel(),
        _ => true,
    };
    let (fixed_mods, variable_mods) = if !has_mods {
        (None, None)
    } else {
        let list = |fixed: bool| -> Result<Vec<String>, String> {
            first
                .mods
                .iter()
                .filter(searched)
                .filter(|m| m.fixed == fixed)
                .map(ModParam::item)
                .collect()
        };
        (Some(list(true)?), Some(list(false)?))
    };
    Ok(SdrfSettings {
        path: PathBuf::new(),
        data_files,
        data_names,
        enzyme: (!first.enzymes.is_empty()).then(|| first.enzymes.join(",")),
        precursor_tol: first.precursor_tol.clone(),
        fragment_tol: first.fragment_tol.clone(),
        fixed_mods,
        variable_mods,
        quant,
        dissociation: first.dissociation.as_deref().and_then(dissociation_of),
        xl: linkers,
    })
}

const SPECTRUM_EXTENSIONS: &[&str] = &["mzML", "mzml", "raw", "d", "mgf"];

fn resolve_data_file(dir: &Path, file: &str) -> PathBuf {
    let p = dir.join(file);
    if p.exists() {
        return p;
    }
    SPECTRUM_EXTENSIONS
        .iter()
        .map(|ext| p.with_extension(ext))
        .find(|c| c.exists())
        .unwrap_or(p)
}

struct RowSettings {
    file: String,
    enzymes: Vec<String>,
    mods: Vec<ModParam>,
    precursor_tol: Option<String>,
    fragment_tol: Option<String>,
    label: Option<String>,
    dissociation: Option<String>,
    xl: Option<String>,
}

fn agree(
    rows: &[RowSettings],
    column: &str,
    key: impl Fn(&RowSettings) -> String,
) -> Result<(), String> {
    let first = key(&rows[0]);
    if let Some(r) = rows.iter().find(|r| key(r) != first) {
        let show = |s: &str| {
            if s.is_empty() {
                "(none)".to_string()
            } else {
                s.to_string()
            }
        };
        return Err(format!(
            "`{column}` differs between files: {} has {}, {} has {}; one andes run searches every \
             file with one setting, so split the SDRF by setting",
            rows[0].file,
            show(&first),
            r.file,
            show(&key(r))
        ));
    }
    Ok(())
}

fn norm_header(h: &str) -> String {
    let h = h.trim().to_ascii_lowercase();
    let h = h.split_whitespace().collect::<Vec<_>>().join(" ");
    h.replace(" [", "[").replace("[ ", "[").replace(" ]", "]")
}

fn value(cell: &str) -> Option<String> {
    let v = cell.trim();
    let lower = v.to_ascii_lowercase();
    if v.is_empty() || lower == "not available" || lower == "not applicable" {
        None
    } else {
        Some(v.to_string())
    }
}

/// The `KEY=value` pairs of an SDRF cell. A key given twice is an error: the cell holds two
/// terms (`NT=Trypsin/P;AC=MS:1001313;NT=Lys-C;AC=MS:1001309`), and each needs its own
/// repeated column; reading only one of them would search with the wrong settings.
fn key_values(v: &str) -> Result<BTreeMap<String, String>, String> {
    let mut out = BTreeMap::new();
    for part in v.split(';') {
        match part.split_once('=') {
            Some((k, val)) => {
                let key = k.trim().to_ascii_uppercase();
                if out.insert(key.clone(), val.trim().to_string()).is_some() {
                    return Err(format!(
                        "`{v}`: {key}= appears twice in one cell; give each term its own \
                         column (repeat the column header)"
                    ));
                }
            }
            None if !part.trim().is_empty() => {
                if out.insert("NT".to_string(), part.trim().to_string()).is_some() {
                    return Err(format!(
                        "`{v}`: two terms in one cell; give each term its own column (repeat \
                         the column header)"
                    ));
                }
            }
            None => {}
        }
    }
    Ok(out)
}

fn nt(v: &str) -> Result<String, String> {
    Ok(key_values(v)?
        .remove("NT")
        .unwrap_or_else(|| v.trim().to_string()))
}

fn enzyme_of(v: &str) -> Result<String, String> {
    let given = nt(v).map_err(|e| format!("`{CLEAVAGE}` {e}"))?;
    let name = given.to_ascii_lowercase();
    let andes = match name.as_str() {
        "trypsin" | "trypsin/p" => "trypsin",
        "lys-c" | "lys-c/p" | "lysc" | "endoproteinase lys-c" => "lysc",
        "lys-n" | "lysn" => "lysn",
        "asp-n" | "asp-n/b" | "aspn" | "asp-n_ambic" => "aspn",
        "glu-c" | "gluc" | "glutamyl endopeptidase" | "v8-de" | "v8-e" => "gluc",
        "chymotrypsin" | "chymotrypsin/p" => "chymotrypsin",
        "arg-c" | "arg-c/p" | "argc" => "argc",
        "alpha-lytic protease" | "alpha-lp" | "alp" => "alphalp",
        "unspecific cleavage" | "nonspecific" | "non-specific" => "nonspecific",
        "no cleavage" | "nocleavage" => "nocleavage",
        _ => {
            return Err(format!(
                "`{CLEAVAGE}`: protease `{given}` has no andes equivalent; give --enzyme"
            ))
        }
    };
    Ok(andes.to_string())
}

fn tolerance(v: &str) -> Result<String, String> {
    let lower = v.trim().to_ascii_lowercase();
    for (unit, scale, out) in [("ppm", 1.0, "ppm"), ("mmu", 0.001, "da"), ("da", 1.0, "da")] {
        if let Some(num) = lower.strip_suffix(unit) {
            let x: f64 = num
                .trim()
                .parse()
                .map_err(|_| format!("tolerance `{v}`: `{}` is not a number", num.trim()))?;
            return Ok(format!("{}{out}", x * scale));
        }
    }
    Err(format!(
        "tolerance `{v}`: expected a value in ppm, Da or mmu"
    ))
}

#[derive(Debug, Clone, PartialEq)]
struct ModParam {
    names: Vec<String>,
    fixed: bool,
    position: String,
    residues: String,
}

impl ModParam {
    fn parse(v: &str) -> Result<Self, String> {
        let kv = key_values(v).map_err(|e| format!("`{MODS}` {e}"))?;
        let mut names: Vec<String> = kv.get("NT").cloned().into_iter().collect();
        if let Some(ac) = kv
            .get("AC")
            .filter(|a| a.to_ascii_uppercase().starts_with("UNIMOD:"))
        {
            names.push(ac.to_ascii_uppercase());
        }
        if names.is_empty() {
            return Err(format!("`{MODS}` `{v}`: no NT= or UNIMOD AC="));
        }
        let fixed = match kv.get("MT").map(|s| s.to_ascii_lowercase()) {
            Some(t) if t == "fixed" => true,
            Some(t) if t == "variable" => false,
            _ => return Err(format!("`{MODS}` `{v}`: MT= must be Fixed or Variable")),
        };
        let position = match kv.get("PP").map(|s| s.to_ascii_lowercase()).as_deref() {
            None | Some("anywhere") => "Anywhere",
            Some("protein n-term") => "Protein N-term",
            Some("protein c-term") => "Protein C-term",
            Some("any n-term") | Some("n-term") => "N-term",
            Some("any c-term") | Some("c-term") => "C-term",
            Some(other) => return Err(format!("`{MODS}` `{v}`: unknown PP= `{other}`")),
        }
        .to_string();
        let residues: String = kv
            .get("TA")
            .map(|ta| {
                ta.split(',')
                    .map(str::trim)
                    .filter(|r| {
                        r.len() == 1 && r.chars().all(|c| c.is_ascii_alphabetic()) && *r != "X"
                    })
                    .collect::<String>()
                    .to_ascii_uppercase()
            })
            .unwrap_or_default();
        if position == "Anywhere" && residues.is_empty() {
            return Err(format!("`{MODS}` `{v}`: PP=Anywhere needs TA= residues"));
        }
        Ok(ModParam {
            names,
            fixed,
            position,
            residues,
        })
    }

    fn item(&self) -> Result<String, String> {
        let sites = match (self.position.as_str(), self.residues.is_empty()) {
            ("Anywhere", _) => self.residues.clone(),
            (pos, true) => pos.to_string(),
            (pos, false) => format!("{pos} {}", self.residues),
        };
        let mut last = String::new();
        for name in &self.names {
            let item = format!("{name} ({sites})");
            match model::mod_spec::parse_mod_list(&item, self.fixed) {
                Ok(_) => return Ok(item),
                Err(e) => last = e.to_string(),
            }
        }
        Err(format!(
            "`{MODS}`: {} ({sites}) does not resolve: {last}",
            self.names.join(" / ")
        ))
    }

    fn key(&self) -> String {
        format!(
            "{} ({} {}) {}",
            self.names.join("/"),
            self.position,
            self.residues,
            if self.fixed { "fixed" } else { "variable" }
        )
    }

    /// A SILAC label on every residue it names (`TA=K,R` holds Lys6 and Arg6), or dimethyl.
    fn is_label_channel(&self) -> bool {
        let silac = !self.residues.is_empty()
            && self.residues.chars().all(|r| self.silac_token(&r.to_string()).is_some());
        silac || self.names.iter().any(|n| dimethyl_name(n))
    }

    fn silac_token(&self, residues: &str) -> Option<(&'static str, f64)> {
        self.names.iter().find_map(|n| silac_token(n, residues))
    }
}

fn silac_token(name: &str, residues: &str) -> Option<(&'static str, f64)> {
    let n = name.to_ascii_lowercase();
    let r = residues;
    let tok = match (n.as_str(), r) {
        ("label:2h(4)" | "unimod:481", "K") => ("Lys4", 4.0),
        ("label:13c(6)" | "unimod:188", "K") => ("Lys6", 6.0),
        ("label:13c(6)15n(2)" | "unimod:259", "K") => ("Lys8", 8.0),
        ("label:13c(6)" | "unimod:188", "R") => ("Arg6", 6.0),
        ("label:13c(6)15n(4)" | "unimod:267", "R") => ("Arg10", 10.0),
        ("label:2h(3)" | "unimod:262", "L") => ("Leu3", 3.0),
        _ => return None,
    };
    Some(tok)
}

fn dimethyl_name(name: &str) -> bool {
    let n = name.to_ascii_lowercase();
    n.starts_with("dimethyl") || matches!(n.as_str(), "unimod:36" | "unimod:199" | "unimod:510")
}

fn label_kind(label: &str) -> &str {
    let l = label.to_ascii_lowercase();
    if l == "label free sample" {
        "lfq"
    } else if l.starts_with("tmt") {
        "tmt"
    } else if l.starts_with("itraq") {
        "itraq"
    } else if l.starts_with("silac") {
        "silac"
    } else if l.starts_with("dimethyl") {
        "dimethyl"
    } else {
        label
    }
}

const ISOBARIC: &[(&str, &[&str])] = &[
    ("tmt6", &["126", "127", "128", "129", "130", "131"]),
    (
        "tmt10",
        &[
            "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131",
        ],
    ),
    (
        "tmt11",
        &[
            "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C",
        ],
    ),
    (
        "tmt16",
        &[
            "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C",
            "132N", "132C", "133N", "133C", "134N",
        ],
    ),
    (
        "tmt18",
        &[
            "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C",
            "132N", "132C", "133N", "133C", "134N", "134C", "135N",
        ],
    ),
    ("itraq4", &["114", "115", "116", "117"]),
    (
        "itraq8",
        &["113", "114", "115", "116", "117", "118", "119", "121"],
    ),
];

fn isobaric_channel(label: &str) -> String {
    let l = label.to_ascii_uppercase();
    let tail = l.rsplit(['-', ' ']).next().unwrap_or(&l);
    let tail = tail
        .trim_start_matches("TMTPRO")
        .trim_start_matches("TMT")
        .trim_start_matches("ITRAQ");
    tail.to_string()
}

fn quant_of(labels: &BTreeSet<String>, mods: &[ModParam]) -> Result<Option<Quant>, String> {
    let Some(first) = labels.iter().next() else {
        return Ok(None);
    };
    let kind = label_kind(first).to_string();
    match kind.as_str() {
        "lfq" => Ok(Some(Quant::Lfq)),
        "tmt" | "itraq" => {
            let channels: BTreeSet<String> = labels.iter().map(|l| isobaric_channel(l)).collect();
            for (plex, set) in ISOBARIC
                .iter()
                .filter(|(p, _)| p.starts_with(kind.as_str()))
            {
                let fits = channels
                    .iter()
                    .all(|c| set.contains(&c.as_str()) || set.contains(&format!("{c}N").as_str()));
                if fits {
                    return Ok(Some(Quant::Isobaric(plex.to_string())));
                }
            }
            Err(format!(
                "`{LABEL}`: no {} plex holds the channels {}; give --tmt",
                kind.to_ascii_uppercase(),
                channels.iter().cloned().collect::<Vec<_>>().join(", ")
            ))
        }
        "silac" => {
            let n = labels.len();
            if !(2..=3).contains(&n) {
                return Err(format!("`{LABEL}`: {n} SILAC channels; andes takes 2 or 3"));
            }
            let mut per_residue: BTreeMap<String, Vec<(&str, f64)>> = BTreeMap::new();
            for m in mods {
                for r in m.residues.chars() {
                    if let Some(t) = m.silac_token(&r.to_string()) {
                        let v = per_residue.entry(r.to_string()).or_default();
                        if !v.contains(&t) {
                            v.push(t);
                        }
                    }
                }
            }
            if per_residue.is_empty() {
                return Err(format!(
                    "`{LABEL}`: SILAC labels but no Label: modification in `{MODS}`; give --label"
                ));
            }
            let mut medium = Vec::new();
            let mut heavy = Vec::new();
            for v in per_residue.values_mut() {
                v.sort_by(|a, b| a.1.total_cmp(&b.1));
                heavy.push(v[v.len() - 1].0);
                if n == 3 {
                    if v.len() < 2 {
                        return Err(format!(
                            "`{LABEL}`: three SILAC channels but one label on a residue; give --label"
                        ));
                    }
                    medium.push(v[v.len() - 2].0);
                }
            }
            let spec = if n == 3 {
                format!("[][{}][{}]", medium.join(","), heavy.join(","))
            } else {
                format!("[][{}]", heavy.join(","))
            };
            Ok(Some(Quant::Label(spec)))
        }
        "dimethyl" => {
            let set: BTreeSet<String> = labels.iter().map(|l| l.to_ascii_uppercase()).collect();
            let two: BTreeSet<String> = ["DIMETHYL0", "DIMETHYL4"].map(String::from).into();
            let three: BTreeSet<String> = ["DIMETHYL0", "DIMETHYL4", "DIMETHYL8"]
                .map(String::from)
                .into();
            if set == two {
                Ok(Some(Quant::Label("dimethyl2".into())))
            } else if set == three {
                Ok(Some(Quant::Label("dimethyl3".into())))
            } else {
                Err(format!(
                    "`{LABEL}`: dimethyl channels {set:?}; give --label"
                ))
            }
        }
        other => Err(format!(
            "`{LABEL}`: label `{other}` has no andes quantification; give --tmt, --lfq or --label"
        )),
    }
}

fn dissociation_of(v: &str) -> Option<ActivationMethod> {
    match v.to_ascii_lowercase().as_str() {
        "hcd"
        | "beam-type collision-induced dissociation"
        | "higher energy beam-type collision-induced dissociation" => Some(ActivationMethod::HCD),
        "cid" | "collision-induced dissociation" => Some(ActivationMethod::CID),
        "etd" | "electron transfer dissociation" => Some(ActivationMethod::ETD),
        "ethcd" | "electron transfer/higher-energy collision dissociation" => {
            Some(ActivationMethod::EThcD)
        }
        _ => None,
    }
}

/// Two crosslinker definitions that link the same way: equal once each end's residues are
/// taken as a set (`TA=K,S` and `TA=S,K`).
fn same_linker(a: &xl::crosslinker::Crosslinker, b: &xl::crosslinker::Crosslinker) -> bool {
    let canonical = |x: &xl::crosslinker::Crosslinker| {
        let mut x = x.clone();
        for e in x.ends.iter_mut() {
            e.residues.sort_unstable();
            e.residues.dedup();
        }
        x
    };
    canonical(a) == canonical(b)
}

pub fn apply(
    s: &SdrfSettings,
    args: &mut crate::SearchArgs,
    m: &clap::ArgMatches,
    src: &mut ParamSources,
) -> Result<Vec<String>, String> {
    let held = |src: &ParamSources, id: &str| {
        crate::config::cli_set(m, id) || src.source(id) == Source::Config
    };
    let mut applied = Vec::new();
    let mut kept = Vec::new();
    let mut warnings = Vec::new();
    let names = |v: &[String]| {
        if v.is_empty() {
            "none".to_string()
        } else {
            v.join(", ")
        }
    };

    if held(src, "spectrum") {
        let sdrf_stems: BTreeSet<_> = s.data_files.iter().filter_map(|p| p.file_stem()).collect();
        let extra: Vec<String> = args
            .spectrum
            .iter()
            .filter(|p| p.file_stem().is_some_and(|n| !sdrf_stems.contains(n)))
            .map(|p| p.display().to_string())
            .collect();
        if !extra.is_empty() {
            warnings.push(format!(
                "WARN: --spectrum files not in the SDRF (its settings still apply): {}",
                extra.join(", ")
            ));
        }
        kept.push("spectrum");
    } else {
        let missing: Vec<String> = s
            .data_files
            .iter()
            .filter(|p| !p.exists())
            .map(|p| p.display().to_string())
            .collect();
        if !missing.is_empty() {
            return Err(format!(
                "--sdrf: data files not found next to the SDRF: {} (give --spectrum to point at them)",
                missing.join(", ")
            ));
        }
        args.spectrum = s.data_files.clone();
        src.record_as("spectrum", Source::Sdrf, &s.data_files);
        applied.push(format!("{} data files", s.data_files.len()));
    }
    if let Some(e) = &s.enzyme {
        if held(src, "enzyme") {
            kept.push("enzyme");
        } else {
            args.db.enzyme = crate::cli::parse_enzyme_list(e)?;
            src.record_as("enzyme", Source::Sdrf, e);
            applied.push(format!("enzyme {e}"));
        }
    }
    if let Some(t) = &s.precursor_tol {
        if held(src, "precursor_tol") {
            kept.push("precursor_tol");
        } else {
            args.precursor_tol = crate::parse_precursor_tol(t)?;
            src.record_as("precursor_tol", Source::Sdrf, t);
            applied.push(format!("precursor tolerance {t}"));
        }
    }
    if let Some(f) = &s.fixed_mods {
        if held(src, "fixed_mods") {
            kept.push("fixed_mods");
        } else {
            args.db.fixed_mods = Some(crate::cli::parse_fixed_mods(&f.join(", "))?);
            src.record_as("fixed_mods", Source::Sdrf, f);
            applied.push(format!("fixed mods {}", names(f)));
        }
    }
    if let Some(v) = &s.variable_mods {
        if held(src, "variable_mods") {
            kept.push("variable_mods");
        } else {
            args.db.variable_mods = Some(crate::cli::parse_variable_mods(&v.join(", "))?);
            src.record_as("variable_mods", Source::Sdrf, v);
            applied.push(format!("variable mods {}", names(v)));
        }
    }
    if let Some(q) = &s.quant {
        if ["tmt", "lfq", "label"].iter().any(|&id| held(src, id)) {
            kept.push("quantification");
        } else {
            match q {
                Quant::Lfq => {
                    args.lfq = true;
                    src.record_as("lfq", Source::Sdrf, &true);
                    applied.push("label-free quantification (--lfq)".into());
                }
                Quant::Isobaric(p) => {
                    args.tmt = Some(quant::Plex::parse(p)?);
                    src.record_as("tmt", Source::Sdrf, p);
                    applied.push(format!("--tmt {p}"));
                }
                Quant::Label(spec) => {
                    args.label = Some(spec.clone());
                    src.record_as("label", Source::Sdrf, spec);
                    applied.push(format!("--label {spec}"));
                }
            }
        }
    }
    if !s.xl.is_empty() {
        if held(src, "xl") {
            kept.push("xl");
        } else {
            let mut set: Vec<xl::crosslinker::Crosslinker> = Vec::new();
            for v in &s.xl {
                let x = xl::spec::from_sdrf(v, model::unimod::bundled())
                    .map_err(|e| format!("--sdrf `{XL}`: {e}"))?;
                match set.iter().find(|y| y.name == x.name) {
                    // The same definition written differently (fields or sites in another
                    // order).
                    Some(y) if same_linker(y, &x) => continue,
                    Some(_) => {
                        return Err(format!(
                            "--sdrf `{XL}`: {} appears with different settings in different rows",
                            x.name
                        ))
                    }
                    None => set.push(x),
                }
            }
            let names: Vec<String> = set
                .iter()
                .map(|x| format!("{} ({})", x.name, x.describe_sites()))
                .collect();
            applied.push(format!("crosslink search --xl {}", names.join("; ")));
            args.xl = Some(crate::cli::XlSet(set));
            src.record_without_value("xl", Source::Sdrf);
        }
    }
    args.sdrf_dissociation = s.dissociation;
    let mut out = vec![format!("SDRF {}: {}", s.path.display(), applied.join("; "))];
    out.extend(warnings);
    if !kept.is_empty() {
        out.push(format!(
            "SDRF: the command line or --config set {}; the SDRF values are not used",
            kept.join(", ")
        ));
    }
    if let Some(f) = &s.fragment_tol {
        out.push(format!(
            "SDRF: fragment mass tolerance {f} not applied; andes measures each run's fragment \
             error (--fragment-tol-ppm / --fragment-tol-da fix one)"
        ));
    }
    Ok(out)
}

const SOURCE_NAME: &str = "source name";
const FRACTION: &str = "comment[fraction identifier]";
const TECHNICAL_REPLICATE: &str = "comment[technical replicate]";
const INSTRUMENT: &str = "comment[instrument]";
const POOLED_SAMPLE: &str = "characteristics[pooled sample]";
const BIOLOGICAL_REPLICATE: &str = "characteristics[biological replicate]";
const CHARACTERISTICS: &str = "characteristics[";
/// Cell values (lower case) that mark a channel without sample.
const BLANK_VALUES: [&str; 3] = ["blank", "empty", "empty channel"];

/// Role of an isobaric channel in its plex.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChannelRole {
    /// A sample of the study.
    Sample,
    /// A pooled reference: `characteristics[pooled sample]` is `pooled` or `SN=...`, or
    /// `characteristics[biological replicate]` is `pooled`.
    Reference,
    /// A channel without sample: a `characteristics[...]` cell, or the last `-`/`_`/space
    /// separated word of the source name, is `blank`, `empty` or `empty channel`.
    Blank,
}

impl ChannelRole {
    pub fn as_str(self) -> &'static str {
        match self {
            ChannelRole::Sample => "sample",
            ChannelRole::Reference => "reference",
            ChannelRole::Blank => "blank",
        }
    }
}

/// The sample of one channel of a run.
#[derive(Debug, Clone, PartialEq)]
pub struct ChannelSample {
    /// The SDRF source name.
    pub sample: String,
    pub role: ChannelRole,
}

/// One run (data file) of an isobaric SDRF.
#[derive(Debug, Clone, PartialEq)]
pub struct RunDesign {
    /// The plex (one labelled sample set, 1-based in SDRF order): the files whose channels
    /// have the same samples and technical replicates (`comment[technical replicate]` may
    /// differ between the channels of a file, as when it counts a sample's measurements).
    pub plex: u32,
    pub fraction: Option<String>,
    /// The NT= name of `comment[instrument]`.
    pub instrument: Option<String>,
    /// SDRF channel (the label without its kit prefix, e.g. `127N`) to its sample.
    pub channels: BTreeMap<String, ChannelSample>,
}

/// The channel design of an isobaric SDRF, by run name (the data file name without its
/// extension, and without `.gz`).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct IsobaricDesign {
    pub path: PathBuf,
    pub runs: BTreeMap<String, RunDesign>,
    pub plexes: u32,
}

/// Run name of an SDRF data file, as andes names a run after its spectrum file.
fn run_name(file: &str) -> String {
    let name = file.rsplit(['/', '\\']).next().unwrap_or(file);
    let name = name.strip_suffix(".gz").unwrap_or(name);
    match name.rfind('.') {
        Some(i) if i > 0 => name[..i].to_string(),
        _ => name.to_string(),
    }
}

fn is_blank_word(v: &str) -> bool {
    BLANK_VALUES.contains(&v.trim().to_ascii_lowercase().as_str())
}

/// A plex: per channel its source name and technical replicate.
type PlexKey = Vec<(String, String, Option<String>)>;

/// Reads the channel design of an isobaric SDRF (`andes quant finalize --sdrf`).
pub fn isobaric_design(path: &Path) -> Result<IsobaricDesign, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("--sdrf: cannot read {}: {e}", path.display()))?;
    parse_isobaric_design(&text)
        .map(|mut d| {
            d.path = path.to_path_buf();
            d
        })
        .map_err(|e| format!("--sdrf {}: {e}", path.display()))
}

pub fn parse_isobaric_design(text: &str) -> Result<IsobaricDesign, String> {
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let header: Vec<String> = lines
        .next()
        .ok_or("empty file")?
        .split('\t')
        .map(norm_header)
        .collect();
    let first = |name: &str| header.iter().position(|h| h == name);
    let data_col = first(DATA_FILE).ok_or_else(|| format!("no `{DATA_FILE}` column"))?;
    let label_col = first(LABEL).ok_or_else(|| format!("no `{LABEL}` column"))?;
    let source_col = first(SOURCE_NAME).ok_or_else(|| format!("no `{SOURCE_NAME}` column"))?;
    let (fraction_col, replicate_col) = (first(FRACTION), first(TECHNICAL_REPLICATE));
    let (instrument_col, pooled_col) = (first(INSTRUMENT), first(POOLED_SAMPLE));
    let bio_col = first(BIOLOGICAL_REPLICATE);
    let characteristics: Vec<usize> = header
        .iter()
        .enumerate()
        .filter(|(_, h)| h.starts_with(CHARACTERISTICS))
        .map(|(i, _)| i)
        .collect();

    let mut order: Vec<String> = Vec::new();
    let mut runs: BTreeMap<String, RunDesign> = BTreeMap::new();
    let mut replicates: BTreeMap<(String, String), Option<String>> = BTreeMap::new();
    for (n, line) in lines.enumerate() {
        let row = n + 2;
        let cells: Vec<&str> = line.split('\t').collect();
        if cells.len() != header.len() {
            return Err(format!(
                "row {row} has {} fields, the header {}",
                cells.len(),
                header.len()
            ));
        }
        let cell = |i: usize| value(cells[i]);
        let file = cell(data_col).ok_or_else(|| format!("row {row}: no data file"))?;
        let label = cell(label_col).ok_or_else(|| format!("row {row}: no `{LABEL}`"))?;
        let kind = label_kind(&label);
        if kind != "tmt" && kind != "itraq" {
            return Err(format!(
                "row {row}: `{LABEL}` is `{label}`, not an isobaric (TMT or iTRAQ) channel"
            ));
        }
        let sample = cells[source_col].trim().to_string();
        if sample.is_empty() {
            return Err(format!("row {row}: no `{SOURCE_NAME}`"));
        }
        let last_word = sample.rsplit(['-', '_', ' ']).next().unwrap_or(&sample);
        let blank =
            is_blank_word(last_word) || characteristics.iter().any(|&i| is_blank_word(cells[i]));
        let pooled = pooled_col.and_then(cell).is_some_and(|v| {
            let l = v.to_ascii_lowercase();
            l == "pooled" || l.starts_with("sn=")
        }) || bio_col
            .and_then(cell)
            .is_some_and(|v| v.eq_ignore_ascii_case("pooled"));
        let role = if blank {
            ChannelRole::Blank
        } else if pooled {
            ChannelRole::Reference
        } else {
            ChannelRole::Sample
        };
        let name = run_name(&file);
        let fraction = fraction_col.and_then(cell);
        let technical_replicate = replicate_col.and_then(cell);
        let instrument = instrument_col
            .and_then(cell)
            .map(|v| nt(&v).map_err(|e| format!("row {row}: `{INSTRUMENT}` {e}")))
            .transpose()?;
        let run = runs.entry(name.clone()).or_insert_with(|| {
            order.push(name.clone());
            RunDesign {
                plex: 0,
                fraction: fraction.clone(),
                instrument: instrument.clone(),
                channels: BTreeMap::new(),
            }
        });
        if run.fraction != fraction || run.instrument != instrument {
            return Err(format!(
                "row {row}: {file} has another fraction or instrument than its first row"
            ));
        }
        let channel = isobaric_channel(&label);
        replicates.insert((name.clone(), channel.clone()), technical_replicate);
        let entry = ChannelSample { sample, role };
        if let Some(prev) = run.channels.insert(channel.clone(), entry.clone()) {
            if prev != entry {
                return Err(format!(
                    "row {row}: {file} channel {channel} is given twice, as {} and {}",
                    prev.sample, entry.sample
                ));
            }
        }
    }
    if runs.is_empty() {
        return Err("no data rows".into());
    }

    let mut plex_of: BTreeMap<PlexKey, u32> = BTreeMap::new();
    let mut seen: BTreeMap<(u32, Option<String>), String> = BTreeMap::new();
    for name in &order {
        let run = runs.get_mut(name).expect("every ordered run was inserted");
        let key: PlexKey = run
            .channels
            .iter()
            .map(|(c, s)| {
                let rep = replicates
                    .get(&(name.clone(), c.clone()))
                    .cloned()
                    .flatten();
                (c.clone(), s.sample.clone(), rep)
            })
            .collect();
        let next = plex_of.len() as u32 + 1;
        let plex = *plex_of.entry(key).or_insert(next);
        run.plex = plex;
        if let Some(other) = seen.insert((plex, run.fraction.clone()), name.clone()) {
            return Err(format!(
                "{other} and {name} have the same channel samples, technical replicate and \
                 fraction ({}), so their plexes cannot be told apart; give each file its \
                 `{FRACTION}` or `{TECHNICAL_REPLICATE}`",
                run.fraction.as_deref().unwrap_or("none")
            ));
        }
    }
    Ok(IsobaricDesign {
        path: PathBuf::new(),
        plexes: plex_of.len() as u32,
        runs,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sdrf(header: &[&str], rows: &[&[&str]]) -> String {
        let mut s = header.join("\t");
        for r in rows {
            s.push('\n');
            s.push_str(&r.join("\t"));
        }
        s
    }

    const H: &[&str] = &[
        "source name",
        "Comment[data file]",
        "comment[label]",
        "comment[cleavage agent details]",
        "comment[modification parameters]",
        "comment[modification parameters]",
        "comment[precursor mass tolerance]",
        "comment[fragment mass tolerance]",
        "comment[dissociation method]",
    ];

    #[test]
    fn label_free_trypsin_phospho() {
        let t = sdrf(
            H,
            &[
                &[
                    "s1",
                    "a.raw",
                    "label free sample",
                    "NT=Trypsin;AC=MS:1001251",
                    "NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed;PP=Anywhere",
                    "NT=Phospho;TA=S,T,Y;MT=Variable;PP=Anywhere",
                    "10 ppm",
                    "0.02 Da",
                    "NT=HCD;AC=PRIDE:0000590",
                ],
                &[
                    "s2",
                    "b.raw",
                    "label free sample",
                    "NT=Trypsin;AC=MS:1001251",
                    "NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed;PP=Anywhere",
                    "NT=Phospho;TA=S,T,Y;MT=Variable;PP=Anywhere",
                    "10 ppm",
                    "0.02 Da",
                    "NT=HCD;AC=PRIDE:0000590",
                ],
            ],
        );
        let s = parse(&t, Path::new("/d")).unwrap();
        assert_eq!(
            s.data_files,
            vec![PathBuf::from("/d/a.raw"), PathBuf::from("/d/b.raw")]
        );
        assert_eq!(s.enzyme.as_deref(), Some("trypsin"));
        assert_eq!(s.precursor_tol.as_deref(), Some("10ppm"));
        assert_eq!(s.fixed_mods, Some(vec!["Carbamidomethyl (C)".to_string()]));
        assert_eq!(s.variable_mods, Some(vec!["Phospho (STY)".to_string()]));
        assert_eq!(s.quant, Some(Quant::Lfq));
        assert_eq!(s.dissociation, Some(ActivationMethod::HCD));
        assert_eq!(s.fragment_tol.as_deref(), Some("0.02 Da"));
        model::mod_spec::parse_mod_list(&s.fixed_mods.unwrap().join(", "), true).unwrap();
        model::mod_spec::parse_mod_list(&s.variable_mods.unwrap().join(", "), false).unwrap();
    }

    #[test]
    fn tmt_plex_is_the_smallest_holding_the_labels() {
        let row = |f: &str, l: &str| -> Vec<String> {
            vec![
                "s".into(),
                f.into(),
                l.into(),
                "NT=Trypsin".into(),
                "NT=TMT6plex;TA=K;MT=Fixed;PP=Anywhere".into(),
                "NT=TMT6plex;PP=Any N-term;MT=Fixed".into(),
                "20 ppm".into(),
                "not available".into(),
                "HCD".into(),
            ]
        };
        let labels10 = [
            "TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C", "TMT129N", "TMT129C", "TMT130N",
            "TMT130C", "TMT131",
        ];
        let rows: Vec<Vec<String>> = labels10.iter().map(|l| row("a.raw", l)).collect();
        let refs: Vec<Vec<&str>> = rows
            .iter()
            .map(|r| r.iter().map(String::as_str).collect())
            .collect();
        let slices: Vec<&[&str]> = refs.iter().map(|r| r.as_slice()).collect();
        let s = parse(&sdrf(H, &slices), Path::new(".")).unwrap();
        assert_eq!(s.quant, Some(Quant::Isobaric("tmt10".into())));
        assert_eq!(s.data_files.len(), 1);
        assert_eq!(
            s.fixed_mods,
            Some(vec![
                "TMT6plex (K)".to_string(),
                "TMT6plex (N-term)".to_string()
            ])
        );
        let mut rows11 = rows.clone();
        rows11.push(row("a.raw", "TMT131C"));
        let refs: Vec<Vec<&str>> = rows11
            .iter()
            .map(|r| r.iter().map(String::as_str).collect())
            .collect();
        let slices: Vec<&[&str]> = refs.iter().map(|r| r.as_slice()).collect();
        let s = parse(&sdrf(H, &slices), Path::new(".")).unwrap();
        assert_eq!(s.quant, Some(Quant::Isobaric("tmt11".into())));
    }

    #[test]
    fn silac_channels_come_from_the_label_mods_and_leave_the_search_lists() {
        let silac = |mods: [&str; 2]| {
            let rows = ["SILAC light", "SILAC heavy"]
                .map(|l| ["s", "a.raw", l, "NT=Trypsin", mods[0], mods[1], "10 ppm", "", ""]);
            let slices: Vec<&[&str]> = rows.iter().map(|r| r.as_slice()).collect();
            parse(&sdrf(H, &slices), Path::new(".")).unwrap()
        };
        let s = silac([
            "NT=Label:13C(6)15N(2);TA=K;MT=Variable;PP=Anywhere",
            "NT=Label:13C(6)15N(4);TA=R;MT=Variable;PP=Anywhere",
        ]);
        assert_eq!(s.quant, Some(Quant::Label("[][Lys8,Arg10]".into())));
        assert_eq!(s.variable_mods, Some(vec![]));
        let s = silac(["NT=Label:13C(6);TA=K,R;MT=Variable;PP=Anywhere", ""]);
        assert_eq!(s.quant, Some(Quant::Label("[][Lys6,Arg6]".into())));
        assert_eq!(s.variable_mods, Some(vec![]));
    }

    #[test]
    fn rows_that_disagree_are_an_error_naming_the_column_and_files() {
        let a: &[&str] = &[
            "s",
            "a.raw",
            "label free sample",
            "Trypsin",
            "",
            "",
            "10 ppm",
            "",
            "",
        ];
        let b: &[&str] = &[
            "s",
            "b.raw",
            "label free sample",
            "Lys-C",
            "",
            "",
            "10 ppm",
            "",
            "",
        ];
        let e = parse(&sdrf(H, &[a, b]), Path::new(".")).unwrap_err();
        assert!(
            e.contains("cleavage agent details") && e.contains("a.raw") && e.contains("b.raw"),
            "{e}"
        );
    }

    #[test]
    fn the_name_or_the_accession_whichever_resolves() {
        let item = |v: &str| ModParam::parse(v).unwrap().item();
        assert_eq!(
            item("NT=Acetyl;AC=UNIMOD:67;PP=Protein N-term;MT=Variable").unwrap(),
            "Acetyl (Protein N-term)"
        );
        assert_eq!(
            item("NT=Phosphorylation;AC=UNIMOD:21;TA=S,T;MT=Variable").unwrap(),
            "UNIMOD:21 (ST)"
        );
        assert!(item("NT=Nothing;TA=S;MT=Variable").is_err());
    }

    #[test]
    fn a_converted_data_file_is_found_by_its_stem() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("run1.mzML"), "").unwrap();
        assert_eq!(
            resolve_data_file(dir.path(), "run1.raw"),
            dir.path().join("run1.mzML")
        );
        assert_eq!(
            resolve_data_file(dir.path(), "run2.raw"),
            dir.path().join("run2.raw")
        );
    }

    #[test]
    fn terminal_mods_and_reserved_words() {
        assert_eq!(
            ModParam::parse("NT=Acetyl;AC=UNIMOD:1;PP=Protein N-term;MT=Variable")
                .unwrap()
                .item()
                .unwrap(),
            "Acetyl (Protein N-term)"
        );
        assert_eq!(
            ModParam::parse("NT=Gln->pyro-Glu;TA=Q;PP=Any N-term;MT=Variable")
                .unwrap()
                .item()
                .unwrap(),
            "Gln->pyro-Glu (N-term Q)"
        );
        assert!(
            ModParam::parse("NT=Oxidation;TA=M").is_err(),
            "MT is required"
        );
        assert_eq!(value(" not applicable "), None);
        assert_eq!(tolerance("20 mmu").unwrap(), "0.02da");
        assert_eq!(norm_header("Comment [ Data File ]"), "comment[data file]");
        assert!(enzyme_of("NT=Pepsin").is_err());
        let two = enzyme_of("NT=Trypsin/P;AC=MS:1001313;NT=Lys-C;AC=MS:1001309").unwrap_err();
        assert!(
            two.contains("NT= appears twice") && two.contains(CLEAVAGE),
            "{two}"
        );
        assert_eq!(enzyme_of("NT=Trypsin/P;AC=MS:1001313").unwrap(), "trypsin");
        for both in ["NT=Trypsin;Lys-C", "Lys-C;NT=Trypsin"] {
            assert!(enzyme_of(both).is_err(), "{both}");
        }
        assert_eq!(enzyme_of("Trypsin").unwrap(), "trypsin");
    }

    const HX: &[&str] = &[
        "source name",
        "comment[data file]",
        "comment[label]",
        "comment[cross-linker]",
    ];

    #[test]
    fn the_crosslinkers_are_read_checked_and_listed() {
        let row = |f: &'static str, x: &'static str| -> Vec<&'static str> {
            vec!["s", f, "label free sample", x]
        };
        let dsso = "NT=DSSO;AC=XLMOD:02126;CL=yes;TA=K,S,T,Y,nterm;MH=54.01;ML=85.98";
        let t = sdrf(HX, &[&row("a.raw", dsso), &row("b.raw", dsso)]);
        let s = parse(&t, Path::new(".")).unwrap();
        assert_eq!(s.xl, vec![dsso.to_string()]);
        let t = sdrf(
            HX,
            &[&row("a.raw", dsso), &row("b.raw", "NT=DSS;AC=XLMOD:02001")],
        );
        let s = parse(&t, Path::new(".")).unwrap();
        assert_eq!(
            s.xl,
            vec![dsso.to_string(), "NT=DSS;AC=XLMOD:02001".to_string()]
        );
        let t = sdrf(HX, &[&row("a.raw", dsso), &row("b.raw", "not applicable")]);
        let e = parse(&t, Path::new(".")).unwrap_err();
        assert!(
            e.contains("cross-linker") && e.contains("some files"),
            "{e}"
        );
        let t = sdrf(HX, &[&row("a.raw", "NT=DSS;CL=yes")]);
        let e = parse(&t, Path::new(".")).unwrap_err();
        assert!(e.contains("CL=yes but DSS is not MS-cleavable"), "{e}");
        let t = sdrf(HX, &[&row("a.raw", "not applicable")]);
        assert!(parse(&t, Path::new(".")).unwrap().xl.is_empty());
    }

    #[test]
    fn the_sdrf_crosslinker_applies_below_the_command_line_and_the_config() {
        use clap::{CommandFactory, FromArgMatches};
        let s = SdrfSettings {
            xl: vec!["NT=DSBU;AC=XLMOD:02120;TA=K,nterm".into()],
            ..Default::default()
        };
        let run = |extra: &[&str], yaml: Option<&str>| -> crate::SearchArgs {
            let base = [
                "andes",
                "--spectrum",
                "x",
                "--database",
                "y",
                "--output-pin",
                "o",
            ];
            let m = crate::TopCli::command()
                .try_get_matches_from(base.iter().chain(extra))
                .unwrap();
            let mut a = crate::SearchArgs::from_arg_matches(&m).unwrap();
            let mut src = ParamSources::default();
            if let Some(y) = yaml {
                let c = crate::config::RunConfig::from_yaml_str(y).unwrap();
                crate::config::apply(c, &mut a, &m, &mut src).unwrap();
            }
            apply(&s, &mut a, &m, &mut src).unwrap();
            a
        };
        let x = run(&[], None).xl.unwrap().0.remove(0);
        assert_eq!(x.name, "DSBU");
        assert_eq!(x.describe_sites(), "K, Protein N-term");
        assert_eq!(
            run(&["--xl", "DSS"], None).xl.unwrap().0.remove(0).name,
            "DSS"
        );
        let cfg = Some("search:\n  xl: DSSO\n");
        assert_eq!(run(&[], cfg).xl.unwrap().0.remove(0).name, "DSSO");
        assert_eq!(
            run(&["--xl", "DSS"], cfg).xl.unwrap().0.remove(0).name,
            "DSS"
        );
    }

    #[test]
    fn equal_sdrf_crosslinkers_written_differently_are_one_linker() {
        use clap::{CommandFactory, FromArgMatches};
        let apply_xl = |cells: &[&str]| -> Result<Vec<String>, String> {
            let s = SdrfSettings {
                xl: cells.iter().map(|c| c.to_string()).collect(),
                ..Default::default()
            };
            let base = [
                "andes",
                "--spectrum",
                "x",
                "--database",
                "y",
                "--output-pin",
                "o",
            ];
            let m = crate::TopCli::command().try_get_matches_from(base).unwrap();
            let mut a = crate::SearchArgs::from_arg_matches(&m).unwrap();
            apply(&s, &mut a, &m, &mut ParamSources::default())?;
            Ok(a.xl.unwrap().0.into_iter().map(|x| x.name).collect())
        };
        assert_eq!(
            apply_xl(&["NT=DSSO;AC=XLMOD:02126", "AC=XLMOD:02126;NT=DSSO"]).unwrap(),
            vec!["DSSO"]
        );
        assert_eq!(
            apply_xl(&["NT=DSSO;TA=K,S", "NT=DSSO;TA=S,K"]).unwrap(),
            vec!["DSSO"]
        );
        let e = apply_xl(&["NT=DSSO;TA=K", "NT=DSSO;TA=K,S"]).unwrap_err();
        assert!(e.contains("different settings"), "{e}");
    }

    #[test]
    fn a_spectrum_missing_from_the_sdrf_is_a_line_of_its_own_marked_as_a_warning() {
        use clap::{CommandFactory, FromArgMatches};
        let s = SdrfSettings {
            data_files: vec![PathBuf::from("run1.mzML")],
            ..Default::default()
        };
        let base = [
            "andes",
            "--spectrum",
            "other.mzML",
            "--database",
            "y",
            "--output-pin",
            "o",
        ];
        let m = crate::TopCli::command().try_get_matches_from(base).unwrap();
        let mut a = crate::SearchArgs::from_arg_matches(&m).unwrap();
        let lines = apply(&s, &mut a, &m, &mut ParamSources::default()).unwrap();
        assert!(!lines[0].contains("WARN"), "{lines:?}");
        assert!(
            lines
                .iter()
                .any(|l| l.starts_with("WARN: --spectrum files not in the SDRF")
                    && l.ends_with("other.mzML")),
            "{lines:?}"
        );
    }

    const DESIGN_H: &[&str] = &[
        "source name",
        "characteristics[organism part]",
        "characteristics[biological replicate]",
        "comment[instrument]",
        "comment[label]",
        "comment[fraction identifier]",
        "comment[technical replicate]",
        "comment[data file]",
    ];

    fn design_row<'a>(
        source: &'a str,
        part: &'a str,
        bio: &'a str,
        label: &'a str,
        fraction: &'a str,
        file: &'a str,
    ) -> [&'a str; 8] {
        [
            source,
            part,
            bio,
            "NT=Orbitrap Fusion;AC=MS:1002416",
            label,
            fraction,
            "1",
            file,
        ]
    }

    #[test]
    fn reference_and_blank_channels_get_their_roles_and_fractions_one_plex() {
        let rows = [
            design_row("S1", "Lung", "1", "TMT126", "1", "p1_Fr01.raw"),
            design_row("S2", "Liver", "1", "TMT127N", "1", "p1_Fr01.raw"),
            design_row(
                "PXD-GIS",
                "Global Internal Standard",
                "pooled",
                "TMT131",
                "1",
                "p1_Fr01.raw",
            ),
            design_row("PXD-blank", "blank", "1", "TMT127C", "1", "p1_Fr01.raw"),
            design_row("S1", "Lung", "1", "TMT126", "2", "p1_Fr02.raw"),
            design_row("S2", "Liver", "1", "TMT127N", "2", "p1_Fr02.raw"),
            design_row(
                "PXD-GIS",
                "Global Internal Standard",
                "pooled",
                "TMT131",
                "2",
                "p1_Fr02.raw",
            ),
            design_row("PXD-blank", "blank", "1", "TMT127C", "2", "p1_Fr02.raw"),
            design_row("S3", "Lung", "1", "TMT126", "1", "p2_Fr01.mzML.gz"),
        ];
        let rows: Vec<&[&str]> = rows.iter().map(|r| &r[..]).collect();
        let d = parse_isobaric_design(&sdrf(DESIGN_H, &rows)).unwrap();
        assert_eq!(d.plexes, 2);
        let r1 = &d.runs["p1_Fr01"];
        let r2 = &d.runs["p1_Fr02"];
        assert_eq!((r1.plex, r2.plex, d.runs["p2_Fr01"].plex), (1, 1, 2));
        assert_eq!(r2.fraction.as_deref(), Some("2"));
        assert_eq!(r1.instrument.as_deref(), Some("Orbitrap Fusion"));
        assert_eq!(r1.channels["126"].role, ChannelRole::Sample);
        assert_eq!(r1.channels["131"].role, ChannelRole::Reference);
        assert_eq!(r1.channels["131"].sample, "PXD-GIS");
        assert_eq!(r1.channels["127C"].role, ChannelRole::Blank);
    }

    #[test]
    fn the_pooled_sample_column_and_an_empty_source_name_word_set_the_roles() {
        let h = [
            "source name",
            "characteristics[pooled sample]",
            "comment[label]",
            "comment[data file]",
        ];
        let rows: [&[&str]; 3] = [
            &["ref", "SN=a,b", "TMT126", "x.raw"],
            &["s", "not pooled", "TMT127N", "x.raw"],
            &["chan_empty", "not pooled", "TMT127C", "x.raw"],
        ];
        let d = parse_isobaric_design(&sdrf(&h, &rows)).unwrap();
        let c = &d.runs["x"].channels;
        assert_eq!(c["126"].role, ChannelRole::Reference);
        assert_eq!(c["127N"].role, ChannelRole::Sample);
        assert_eq!(c["127C"].role, ChannelRole::Blank);
        assert_eq!(d.runs["x"].fraction, None);
    }

    #[test]
    fn two_files_of_one_channel_map_and_fraction_are_refused() {
        let rows = [
            design_row("S1", "Lung", "1", "TMT126", "1", "a.raw"),
            design_row("S1", "Lung", "1", "TMT126", "1", "b.raw"),
        ];
        let rows: Vec<&[&str]> = rows.iter().map(|r| &r[..]).collect();
        let e = parse_isobaric_design(&sdrf(DESIGN_H, &rows)).unwrap_err();
        assert!(e.contains("cannot be told apart"), "{e}");
        let rows = [design_row(
            "S1",
            "Lung",
            "1",
            "label free sample",
            "1",
            "a.raw",
        )];
        let rows: Vec<&[&str]> = rows.iter().map(|r| &r[..]).collect();
        let e = parse_isobaric_design(&sdrf(DESIGN_H, &rows)).unwrap_err();
        assert!(e.contains("not an isobaric"), "{e}");
    }

    #[test]
    fn a_technical_replicate_per_channel_keeps_the_fractions_of_a_plex_together() {
        let h = [
            "source name",
            "comment[label]",
            "comment[technical replicate]",
            "comment[fraction identifier]",
            "comment[data file]",
        ];
        let rows: [&[&str]; 8] = [
            &["S1", "TMT126", "2", "1", "p_Fr01.raw"],
            &["S2", "TMT127N", "3", "1", "p_Fr01.raw"],
            &["S1", "TMT126", "2", "2", "p_Fr02.raw"],
            &["S2", "TMT127N", "3", "2", "p_Fr02.raw"],
            &["S1", "TMT126", "4", "1", "q_Fr01.raw"],
            &["S2", "TMT127N", "4", "1", "q_Fr01.raw"],
            &["S1", "TMT126", "2", "3", "p_Fr03.raw"],
            &["S2", "TMT127N", "3", "3", "p_Fr03.raw"],
        ];
        let d = parse_isobaric_design(&sdrf(&h, &rows)).unwrap();
        let plex = |r: &str| d.runs[r].plex;
        assert_eq!((plex("p_Fr01"), plex("p_Fr02"), plex("p_Fr03")), (1, 1, 1));
        assert_eq!(
            plex("q_Fr01"),
            2,
            "a repeat injection of the plex is its own plex"
        );
    }
}
