use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::OnceLock;

#[derive(Debug, Clone, PartialEq)]
pub struct UnimodMod {
    pub accession: String,
    pub name: String,
    pub mono_mass: f64,
    pub sites: Vec<String>,
    pub common_sites: Vec<String>,
    pub isotopic_label: bool,
}

#[derive(thiserror::Error, Debug, PartialEq)]
pub enum UnimodResolveError {
    #[error("unknown Unimod accession {token:?}")]
    UnknownUnimodAccession { token: String },
    #[error("ambiguous modification name {token:?} (matches {accessions:?})")]
    AmbiguousModName {
        token: String,
        accessions: Vec<String>,
    },
}

#[derive(Debug, Default)]
pub struct UnimodDb {
    by_accession: HashMap<String, UnimodMod>,
    by_name: HashMap<String, Vec<String>>,
    aliases: Vec<(String, String)>,
}

impl UnimodDb {
    pub fn parse_obo(text: &str) -> Self {
        let mut db = UnimodDb::default();

        let mut in_term = false;
        let mut cur_id: Option<String> = None;
        let mut cur_name: Option<String> = None;
        let mut cur_mass: Option<f64> = None;
        let mut cur_sites: Vec<String> = Vec::new();
        let mut cur_specs: HashMap<String, SpecFields> = HashMap::new();

        let commit = |db: &mut UnimodDb,
                      id: &mut Option<String>,
                      name: &mut Option<String>,
                      mass: &mut Option<f64>,
                      sites: &mut Vec<String>,
                      specs: &mut HashMap<String, SpecFields>| {
            let mut common_sites: Vec<String> = Vec::new();
            let isotopic_label = !specs.is_empty()
                && specs
                    .values()
                    .all(|f| f.classification.as_deref() == Some("Isotopic label"));
            for f in specs.drain().map(|(_, v)| v) {
                if let (Some(site), Some(false)) = (f.site, f.hidden) {
                    if !common_sites.contains(&site) {
                        common_sites.push(site);
                    }
                }
            }
            common_sites.sort();
            if let (Some(id_v), Some(name_v), Some(mass_v)) = (id.take(), name.take(), mass.take())
            {
                let entry = UnimodMod {
                    accession: id_v.clone(),
                    name: name_v.clone(),
                    mono_mass: mass_v,
                    sites: std::mem::take(sites),
                    common_sites,
                    isotopic_label,
                };
                db.by_name
                    .entry(name_v.to_ascii_lowercase())
                    .or_default()
                    .push(id_v.clone());
                db.by_accession.insert(id_v, entry);
            } else {
                sites.clear();
            }
        };

        for raw in text.lines() {
            let line = raw.trim();
            if line == "[Term]" {
                commit(
                    &mut db,
                    &mut cur_id,
                    &mut cur_name,
                    &mut cur_mass,
                    &mut cur_sites,
                    &mut cur_specs,
                );
                in_term = true;
                continue;
            }
            if line.starts_with('[') && line.ends_with(']') {
                commit(
                    &mut db,
                    &mut cur_id,
                    &mut cur_name,
                    &mut cur_mass,
                    &mut cur_sites,
                    &mut cur_specs,
                );
                in_term = false;
                continue;
            }
            if !in_term {
                continue;
            }

            if let Some(rest) = line.strip_prefix("id:") {
                let id = rest.trim();
                if id.starts_with("UNIMOD:") {
                    cur_id = Some(id.to_string());
                }
            } else if let Some(rest) = line.strip_prefix("name:") {
                cur_name = Some(rest.trim().to_string());
            } else if let Some(rest) = line.strip_prefix("xref:") {
                let rest = rest.trim();
                if let Some(val) = rest.strip_prefix("delta_mono_mass") {
                    if let Some(v) = parse_quoted_f64(val) {
                        cur_mass = Some(v);
                    }
                } else if let Some(site) = parse_spec_site(rest) {
                    if let Some(idx) = spec_index(rest) {
                        cur_specs.entry(idx).or_default().site = Some(site.clone());
                    }
                    if !cur_sites.contains(&site) {
                        cur_sites.push(site);
                    }
                } else if let Some((idx, hidden)) = parse_spec_hidden(rest) {
                    cur_specs.entry(idx).or_default().hidden = Some(hidden);
                } else if let Some((idx, class)) = parse_spec_value(rest, "_classification") {
                    cur_specs.entry(idx).or_default().classification = Some(class);
                }
            }
        }
        commit(
            &mut db,
            &mut cur_id,
            &mut cur_name,
            &mut cur_mass,
            &mut cur_sites,
            &mut cur_specs,
        );

        db
    }

    pub fn resolve(&self, token: &str) -> Result<Option<&UnimodMod>, UnimodResolveError> {
        let t = token.trim();
        if let Some(num) = t
            .get(..7)
            .filter(|p| p.eq_ignore_ascii_case("UNIMOD:"))
            .map(|_| &t[7..])
        {
            let key = format!("UNIMOD:{}", num.trim());
            return match self.by_accession.get(&key) {
                Some(m) => Ok(Some(m)),
                None => Err(UnimodResolveError::UnknownUnimodAccession {
                    token: token.to_string(),
                }),
            };
        }

        match self.by_name.get(&t.to_ascii_lowercase()) {
            None => Ok(None),
            Some(accs) if accs.len() == 1 => Ok(self.by_accession.get(&accs[0])),
            Some(accs) => Err(UnimodResolveError::AmbiguousModName {
                token: token.to_string(),
                accessions: accs.clone(),
            }),
        }
    }

    pub fn with_aliases(mut self, text: &str) -> Self {
        for line in text.lines() {
            let line = line.split('#').next().unwrap_or("").trim();
            let Some((alias, accession)) = line.split_once('\t') else {
                continue;
            };
            let (alias, accession) = (alias.trim(), accession.trim());
            if !alias.is_empty() && self.by_accession.contains_key(accession) {
                self.aliases.push((alias.to_string(), accession.to_string()));
            }
        }
        self
    }

    pub fn resolve_alias(&self, token: &str) -> Option<&UnimodMod> {
        let t = token.trim();
        self.aliases
            .iter()
            .find(|(alias, _)| alias.eq_ignore_ascii_case(t))
            .and_then(|(_, accession)| self.by_accession.get(accession))
    }

    pub fn mods(&self) -> impl Iterator<Item = &UnimodMod> {
        self.by_accession.values()
    }

    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.by_accession
            .values()
            .map(|m| m.name.as_str())
            .chain(self.aliases.iter().map(|(a, _)| a.as_str()))
    }

    pub fn by_accession(&self, accession: &str) -> Option<&UnimodMod> {
        self.by_accession.get(accession.trim())
    }

    pub fn len(&self) -> usize {
        self.by_accession.len()
    }

    pub fn is_empty(&self) -> bool {
        self.by_accession.is_empty()
    }
}

fn parse_quoted_f64(s: &str) -> Option<f64> {
    let start = s.find('"')? + 1;
    let end = s[start..].find('"')? + start;
    s[start..end].trim().parse().ok()
}

#[derive(Default)]
struct SpecFields {
    site: Option<String>,
    hidden: Option<bool>,
    classification: Option<String>,
}

fn parse_spec_value(s: &str, field: &str) -> Option<(String, String)> {
    let idx = spec_index(s)?;
    let rest = s.strip_prefix("spec_")?;
    let val = rest[idx.len()..].strip_prefix(field)?;
    let start = val.find('"')? + 1;
    let end = val[start..].find('"')? + start;
    Some((idx, val[start..end].trim().to_string()))
}

fn spec_index(s: &str) -> Option<String> {
    let rest = s.strip_prefix("spec_")?;
    let underscore = rest.find('_')?;
    Some(rest[..underscore].to_string())
}

fn parse_spec_hidden(s: &str) -> Option<(String, bool)> {
    parse_spec_value(s, "_hidden").map(|(idx, v)| (idx, v == "1"))
}

fn parse_spec_site(s: &str) -> Option<String> {
    let rest = s.strip_prefix("spec_")?;
    let underscore = rest.find('_')?;
    let tail = &rest[underscore..];
    let val = tail.strip_prefix("_site")?;
    let start = val.find('"')? + 1;
    let end = val[start..].find('"')? + start;
    let site = val[start..end].trim();
    if site.is_empty() {
        None
    } else {
        Some(site.to_string())
    }
}

pub(crate) fn bundled_path(file: &str) -> PathBuf {
    let mut roots: Vec<PathBuf> = Vec::new();
    if let Ok(exe) = std::env::current_exe() {
        if let Some(dir) = exe.parent() {
            roots.push(dir.join("resources"));
        }
    }
    roots.push(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../resources"));

    for root in &roots {
        let p = root.join(file);
        if p.exists() {
            return p;
        }
    }
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../resources")
        .join(file)
}

pub fn bundled() -> &'static UnimodDb {
    static DB: OnceLock<UnimodDb> = OnceLock::new();
    DB.get_or_init(|| {
        let db = match std::fs::read_to_string(bundled_path("unimod.obo")) {
            Ok(text) => UnimodDb::parse_obo(&text),
            Err(_) => UnimodDb::default(),
        };
        match std::fs::read_to_string(bundled_path("unimod_aliases.tsv")) {
            Ok(text) => db.with_aliases(&text),
            Err(_) => db,
        }
    })
}

pub fn resolve(token: &str) -> Result<Option<&'static UnimodMod>, UnimodResolveError> {
    bundled().resolve(token)
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = r#"
format-version: 1.2

[Term]
id: UNIMOD:35
name: Oxidation
xref: delta_mono_mass "15.994915"
xref: delta_composition "O"
xref: spec_1_site "M"
xref: spec_1_position "Anywhere"
xref: spec_2_site "C"

[Term]
id: UNIMOD:4
name: Carbamidomethyl
xref: delta_mono_mass "57.021464"
xref: spec_1_site "C"

[Term]
id: UNIMOD:9999
name: NoMass
xref: delta_composition "X"
"#;

    fn db() -> UnimodDb {
        UnimodDb::parse_obo(SAMPLE)
    }

    #[test]
    fn parses_accession_and_name() {
        let db = db();
        let ox = db.resolve("UNIMOD:35").unwrap().unwrap();
        assert_eq!(ox.name, "Oxidation");
        assert!((ox.mono_mass - 15.994915).abs() < 1e-9);
        assert_eq!(ox.sites, vec!["M".to_string(), "C".to_string()]);
    }

    #[test]
    fn resolves_bare_name_case_insensitive() {
        let db = db();
        let m = db.resolve("oxidation").unwrap().unwrap();
        assert_eq!(m.accession, "UNIMOD:35");
        let m2 = db.resolve("CARBAMIDOMETHYL").unwrap().unwrap();
        assert_eq!(m2.accession, "UNIMOD:4");
    }

    #[test]
    fn unknown_accession_is_error() {
        let db = db();
        assert_eq!(
            db.resolve("UNIMOD:123456"),
            Err(UnimodResolveError::UnknownUnimodAccession {
                token: "UNIMOD:123456".to_string()
            })
        );
    }

    #[test]
    fn unknown_name_is_none() {
        let db = db();
        assert_eq!(db.resolve("Frobnicate"), Ok(None));
    }

    #[test]
    fn term_without_mass_is_skipped() {
        let db = db();
        assert!(db.resolve("UNIMOD:9999").is_err());
        assert_eq!(db.resolve("NoMass"), Ok(None));
    }

    #[test]
    fn ambiguous_name_errors() {
        let dup = r#"
[Term]
id: UNIMOD:1
name: Dup
xref: delta_mono_mass "1.0"

[Term]
id: UNIMOD:2
name: dup
xref: delta_mono_mass "2.0"
"#;
        let db = UnimodDb::parse_obo(dup);
        match db.resolve("Dup") {
            Err(UnimodResolveError::AmbiguousModName { accessions, .. }) => {
                assert_eq!(accessions.len(), 2);
            }
            other => panic!("expected ambiguous, got {other:?}"),
        }
    }

    #[test]
    fn bundled_obo_resolves_curated_set() {
        let curated: &[(&str, &str, f64)] = &[
            ("UNIMOD:4", "Carbamidomethyl", 57.021464),
            ("UNIMOD:35", "Oxidation", 15.994915),
            ("UNIMOD:21", "Phospho", 79.966331),
            ("UNIMOD:1", "Acetyl", 42.010565),
            ("UNIMOD:737", "TMT6plex", 229.162932),
            ("UNIMOD:214", "iTRAQ4plex", 144.102063),
            ("UNIMOD:730", "iTRAQ8plex", 304.205360),
            ("UNIMOD:7", "Deamidated", 0.984016),
            ("UNIMOD:121", "GG", 114.042927),
        ];
        let db = bundled();
        assert!(!db.is_empty(), "bundled OBO did not parse");
        for (acc, name, mass) in curated {
            let by_acc = db.resolve(acc).unwrap().unwrap_or_else(|| {
                panic!("accession {acc} did not resolve");
            });
            assert_eq!(&by_acc.name, name, "name mismatch for {acc}");
            assert!(
                (by_acc.mono_mass - mass).abs() < 1e-5,
                "mass mismatch for {acc}: obo={}, ref={}",
                by_acc.mono_mass,
                mass
            );
            let by_name = db.resolve(name).unwrap().unwrap();
            assert_eq!(by_name.accession, *acc, "name {name} resolved elsewhere");
        }
    }
}
