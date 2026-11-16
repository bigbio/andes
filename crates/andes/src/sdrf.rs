use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use model::activation::ActivationMethod;

/// The quantification an SDRF's labels name.
#[derive(Debug, Clone, PartialEq)]
pub enum Quant {
    Lfq,
    /// `--tmt` plex name (`tmt10`, `itraq4`).
    Isobaric(String),
    /// `--label` spec (`[][Lys8,Arg10]`, `dimethyl2`).
    Label(String),
}

/// The settings read from one SDRF.
#[derive(Debug, Clone, Default)]
pub struct SdrfSettings {
    pub path: PathBuf,
    /// Spectrum files, in order of first appearance, resolved next to the SDRF.
    pub data_files: Vec<PathBuf>,
    /// Protease list as `--enzyme` takes it (`trypsin`, `trypsin+lysc`).
    pub enzyme: Option<String>,
    /// `--precursor-tol` value (`10ppm`).
    pub precursor_tol: Option<String>,
    /// The fragment tolerance as written, reported but not applied.
    pub fragment_tol: Option<String>,
    /// `Name (sites)` items; `None` when the SDRF has no modification column.
    pub fixed_mods: Option<Vec<String>>,
    pub variable_mods: Option<Vec<String>>,
    pub quant: Option<Quant>,
    pub dissociation: Option<ActivationMethod>,
    pub xl: Option<String>,
}

const DATA_FILE: &str = "comment[data file]";
const CLEAVAGE: &str = "comment[cleavage agent details]";
const MODS: &str = "comment[modification parameters]";
const PRECURSOR_TOL: &str = "comment[precursor mass tolerance]";
const FRAGMENT_TOL: &str = "comment[fragment mass tolerance]";
const LABEL: &str = "comment[label]";
const DISSOCIATION: &str = "comment[dissociation method]";
const XL: &str = "comment[cross-linker]";

/// Read and check an SDRF file.
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

/// Parse SDRF text; data files resolve against `dir`.
pub fn parse(text: &str, dir: &Path) -> Result<SdrfSettings, String> {
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
        let file = cell(data_col).ok_or_else(|| format!("row {}: no data file", n + 2))?;
        let mut enzymes = Vec::new();
        for &i in &cleave {
            if let Some(v) = cell(i) {
                let e = enzyme_of(&v)?;
                if !enzymes.contains(&e) {
                    enzymes.push(e);
                }
            }
        }
        let mut row_mods = Vec::new();
        for &i in &mods {
            if let Some(v) = cell(i) {
                row_mods.push(ModParam::parse(&v)?);
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
            label: label_col.first().and_then(|&i| cell(i)).map(|v| nt(&v)),
            dissociation: diss.first().and_then(|&i| cell(i)).map(|v| nt(&v)),
            xl: xl_col.first().and_then(|&i| cell(i)),
        });
    }
    if rows.is_empty() {
        return Err("no data rows".into());
    }

    agree(&rows, CLEAVAGE, |r| r.enzymes.join("+"))?;
    agree(&rows, MODS, |r| {
        let mut m: Vec<String> = r.mods.iter().map(|m| m.key()).collect();
        m.sort();
        m.join(", ")
    })?;
    agree(&rows, PRECURSOR_TOL, |r| {
        r.precursor_tol.clone().unwrap_or_default()
    })?;
    agree(&rows, DISSOCIATION, |r| {
        r.dissociation.clone().unwrap_or_default()
    })?;
    agree(&rows, LABEL, |r| {
        r.label.as_deref().map(label_kind).unwrap_or("").to_string()
    })?;
    agree(&rows, XL, |r| r.xl.clone().unwrap_or_default())?;
    if let Some(v) = &rows[0].xl {
        andes_xl::spec::from_sdrf(v, model::unimod::bundled())
            .map_err(|e| format!("`{XL}`: {e}"))?;
    }

    let first = &rows[0];
    let mut data_files: Vec<PathBuf> = Vec::new();
    for r in &rows {
        let p = resolve_data_file(dir, &r.file);
        if !data_files.contains(&p) {
            data_files.push(p);
        }
    }
    let labels: BTreeSet<String> = rows.iter().filter_map(|r| r.label.clone()).collect();
    let quant = quant_of(&labels, &first.mods)?;
    // Labels the quantification searches itself leave the modification lists.
    let searched = |m: &&ModParam| match &quant {
        Some(Quant::Label(_)) => !m.is_label_channel(),
        _ => true,
    };
    let (fixed_mods, variable_mods) = if mods.is_empty() {
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
        enzyme: (!first.enzymes.is_empty()).then(|| first.enzymes.join("+")),
        precursor_tol: first.precursor_tol.clone(),
        fragment_tol: first.fragment_tol.clone(),
        fixed_mods,
        variable_mods,
        quant,
        dissociation: first.dissociation.as_deref().and_then(dissociation_of),
        xl: first.xl.clone(),
    })
}

/// Spectrum formats andes reads, tried in this order when the SDRF's file is
/// not there under its own name (an SDRF names the vendor files; pipelines often
/// convert them).
const SPECTRUM_EXTENSIONS: &[&str] = &["mzML", "mzml", "raw", "d", "mgf"];

/// `file` in `dir`, or the same stem with another spectrum extension when only
/// that exists; the SDRF's own name when neither does (reported as missing).
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

/// Every row must give the same `key` for `column`.
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

/// `NT=a;AC=b` -> {NT: a, AC: b}; a cell without `=` is its own NT.
fn key_values(v: &str) -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    for part in v.split(';') {
        match part.split_once('=') {
            Some((k, val)) => {
                out.insert(k.trim().to_ascii_uppercase(), val.trim().to_string());
            }
            None if !part.trim().is_empty() => {
                out.entry("NT".to_string())
                    .or_insert_with(|| part.trim().to_string());
            }
            None => {}
        }
    }
    out
}

fn nt(v: &str) -> String {
    key_values(v)
        .remove("NT")
        .unwrap_or_else(|| v.trim().to_string())
}

/// The andes protease of an SDRF cleavage agent.
fn enzyme_of(v: &str) -> Result<String, String> {
    let name = nt(v).to_ascii_lowercase();
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
                "`{CLEAVAGE}`: protease `{}` has no andes equivalent; give --enzyme",
                nt(v)
            ))
        }
    };
    Ok(andes.to_string())
}

/// `10 ppm`, `0.02 Da`, `20 mmu` -> the `--precursor-tol` form.
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

/// One `comment[modification parameters]` cell.
#[derive(Debug, Clone, PartialEq)]
struct ModParam {
    /// The name (`NT=`), then the UNIMOD accession (`AC=`) when given: the
    /// first that the UNIMOD table resolves is used, so a wrong accession next
    /// to a right name (or the reverse) still reads.
    names: Vec<String>,
    fixed: bool,
    /// `Anywhere`, `Protein N-term`, `Protein C-term`, `N-term`, `C-term`.
    position: String,
    residues: String,
}

impl ModParam {
    fn parse(v: &str) -> Result<Self, String> {
        let kv = key_values(v);
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

    /// A SILAC or dimethyl channel label, which `--label` searches itself.
    fn is_label_channel(&self) -> bool {
        self.silac_token(&self.residues).is_some() || self.names.iter().any(|n| dimethyl_name(n))
    }

    /// The SILAC `--label` token of this mod on `residues`, by any of its names.
    fn silac_token(&self, residues: &str) -> Option<(&'static str, f64)> {
        self.names.iter().find_map(|n| silac_token(n, residues))
    }
}

/// The andes `--label` token of a SILAC label modification on one residue.
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

/// `lfq`, `tmt`, `itraq`, `silac`, `dimethyl` or the label itself.
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

/// Channels of each isobaric plex, smallest first.
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

/// The channel of an isobaric label: `TMT127N` -> `127N`, `iTRAQ4plex-114` -> `114`.
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
                // A channel without N/C matches its N form when the plex has no bare one.
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
            // Per residue, the label tokens by mass: heavy is the heaviest.
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

#[derive(Debug, Default, Clone, Copy)]
pub struct ConfigKeys {
    spectrum: bool,
    enzyme: bool,
    fixed_mods: bool,
    variable_mods: bool,
    precursor_tol: bool,
    quant: bool,
    xl: bool,
}

impl ConfigKeys {
    pub fn of(cfg: &crate::config::RunConfig) -> Self {
        ConfigKeys {
            spectrum: cfg.io.spectrum.is_some(),
            enzyme: cfg.search.enzyme.is_some(),
            fixed_mods: cfg.search.fixed_mods.is_some(),
            variable_mods: cfg.search.variable_mods.is_some(),
            precursor_tol: cfg.search.precursor_tol.is_some(),
            quant: cfg.quant.tmt.is_some() || cfg.quant.lfq.is_some() || cfg.quant.label.is_some(),
            xl: cfg.search.xl.is_some(),
        }
    }
}

/// Apply `s` to the search arguments where neither the command line nor the
/// config set the value. Returns the summary lines to print.
pub fn apply(
    s: &SdrfSettings,
    args: &mut crate::SearchArgs,
    m: &clap::ArgMatches,
    cfg: ConfigKeys,
) -> Result<Vec<String>, String> {
    let cli = |id: &str| crate::config::cli_set(m, id);
    let mut applied = Vec::new();
    let mut kept = Vec::new();
    let names = |v: &[String]| {
        if v.is_empty() {
            "none".to_string()
        } else {
            v.join(", ")
        }
    };

    if cli("spectrum") || cfg.spectrum {
        // By stem: a converted file (x.mzML for the SDRF's x.raw) is the same run.
        let sdrf_stems: BTreeSet<_> = s.data_files.iter().filter_map(|p| p.file_stem()).collect();
        let extra: Vec<String> = args
            .spectrum
            .iter()
            .filter(|p| p.file_stem().is_some_and(|n| !sdrf_stems.contains(n)))
            .map(|p| p.display().to_string())
            .collect();
        if !extra.is_empty() {
            applied.push(format!(
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
        applied.push(format!("{} data files", s.data_files.len()));
    }
    if let Some(e) = &s.enzyme {
        if cli("enzyme") || cfg.enzyme {
            kept.push("enzyme");
        } else {
            args.db.enzyme = crate::cli::parse_enzyme_list(e)?;
            applied.push(format!("enzyme {e}"));
        }
    }
    if let Some(t) = &s.precursor_tol {
        if cli("precursor_tol") || cfg.precursor_tol {
            kept.push("precursor_tol");
        } else {
            args.precursor_tol = crate::parse_precursor_tol(t)?;
            applied.push(format!("precursor tolerance {t}"));
        }
    }
    if let Some(f) = &s.fixed_mods {
        if cli("fixed_mods") || cfg.fixed_mods {
            kept.push("fixed_mods");
        } else {
            args.db.fixed_mods = Some(crate::cli::parse_fixed_mods(&f.join(", "))?);
            applied.push(format!("fixed mods {}", names(f)));
        }
    }
    if let Some(v) = &s.variable_mods {
        if cli("variable_mods") || cfg.variable_mods {
            kept.push("variable_mods");
        } else {
            args.db.variable_mods = Some(crate::cli::parse_variable_mods(&v.join(", "))?);
            applied.push(format!("variable mods {}", names(v)));
        }
    }
    if let Some(q) = &s.quant {
        if cli("tmt") || cli("lfq") || cli("label") || cfg.quant {
            kept.push("quantification");
        } else {
            match q {
                Quant::Lfq => {
                    args.lfq = true;
                    applied.push("label-free quantification (--lfq)".into());
                }
                Quant::Isobaric(p) => {
                    args.tmt = Some(quant::Plex::parse(p)?);
                    applied.push(format!("--tmt {p}"));
                }
                Quant::Label(spec) => {
                    args.label = Some(spec.clone());
                    applied.push(format!("--label {spec}"));
                }
            }
        }
    }
    if let Some(v) = &s.xl {
        if cli("xl") || cfg.xl {
            kept.push("xl");
        } else {
            let x = andes_xl::spec::from_sdrf(v, model::unimod::bundled())
                .map_err(|e| format!("--sdrf `{XL}`: {e}"))?;
            applied.push(format!(
                "crosslink search --xl {} ({})",
                x.name,
                x.describe_sites()
            ));
            args.xl = Some(x);
        }
    }
    args.sdrf_dissociation = s.dissociation;
    let mut out = vec![format!("SDRF {}: {}", s.path.display(), applied.join("; "))];
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
        let row = |l: &str| -> Vec<String> {
            vec![
                "s".into(),
                "a.raw".into(),
                l.into(),
                "NT=Trypsin".into(),
                "NT=Label:13C(6)15N(2);TA=K;MT=Variable;PP=Anywhere".into(),
                "NT=Label:13C(6)15N(4);TA=R;MT=Variable;PP=Anywhere".into(),
                "10 ppm".into(),
                "".into(),
                "".into(),
            ]
        };
        let rows = [row("SILAC light"), row("SILAC heavy")];
        let refs: Vec<Vec<&str>> = rows
            .iter()
            .map(|r| r.iter().map(String::as_str).collect())
            .collect();
        let slices: Vec<&[&str]> = refs.iter().map(|r| r.as_slice()).collect();
        let s = parse(&sdrf(H, &slices), Path::new(".")).unwrap();
        assert_eq!(s.quant, Some(Quant::Label("[][Lys8,Arg10]".into())));
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

    /// Real SDRFs carry wrong accessions next to right names (PXD014415:
    /// Acetyl with UNIMOD:67) and the reverse; the first that resolves wins.
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
        // Missing under every extension: the SDRF's own name, reported later.
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
    }

    const HX: &[&str] = &[
        "source name",
        "comment[data file]",
        "comment[label]",
        "comment[cross-linker]",
    ];

    #[test]
    fn the_crosslinker_is_read_checked_and_must_agree() {
        let row = |f: &'static str, x: &'static str| -> Vec<&'static str> {
            vec!["s", f, "label free sample", x]
        };
        let dsso = "NT=DSSO;AC=XLMOD:02126;CL=yes;TA=K,S,T,Y,nterm;MH=54.01;ML=85.98";
        let t = sdrf(HX, &[&row("a.raw", dsso), &row("b.raw", dsso)]);
        let s = parse(&t, Path::new(".")).unwrap();
        assert_eq!(s.xl.as_deref(), Some(dsso));
        let t = sdrf(
            HX,
            &[&row("a.raw", dsso), &row("b.raw", "NT=DSS;AC=XLMOD:02001")],
        );
        let e = parse(&t, Path::new(".")).unwrap_err();
        assert!(e.contains("cross-linker") && e.contains("differs"), "{e}");
        let t = sdrf(HX, &[&row("a.raw", "NT=DSS;CL=yes")]);
        let e = parse(&t, Path::new(".")).unwrap_err();
        assert!(e.contains("CL=yes but DSS is not MS-cleavable"), "{e}");
        let t = sdrf(HX, &[&row("a.raw", "not applicable")]);
        assert_eq!(parse(&t, Path::new(".")).unwrap().xl, None);
    }

    #[test]
    fn the_sdrf_crosslinker_applies_below_the_command_line_and_the_config() {
        use clap::{CommandFactory, FromArgMatches};
        let s = SdrfSettings {
            xl: Some("NT=DSBU;AC=XLMOD:02120;TA=K,nterm".into()),
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
            let mut keys = ConfigKeys::default();
            if let Some(y) = yaml {
                let c = crate::config::RunConfig::from_yaml_str(y).unwrap();
                keys = ConfigKeys::of(&c);
                crate::config::apply(c, &mut a, &m).unwrap();
            }
            apply(&s, &mut a, &m, keys).unwrap();
            a
        };
        let x = run(&[], None).xl.unwrap();
        assert_eq!(x.name, "DSBU");
        assert_eq!(x.describe_sites(), "K, Protein N-term");
        assert_eq!(run(&["--xl", "DSS"], None).xl.unwrap().name, "DSS");
        let cfg = Some("search:\n  xl: DSSO\n");
        assert_eq!(run(&[], cfg).xl.unwrap().name, "DSSO");
        assert_eq!(run(&["--xl", "DSS"], cfg).xl.unwrap().name, "DSS");
    }
}
