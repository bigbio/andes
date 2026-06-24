//! Unimod accession / name → monoisotopic mass resolver.
//!
//! Parses the bundled `resources/unimod.obo` (a simple OBO stanza format) into
//! two lookup maps so that `--mods` and the peptide round-trip can accept a
//! `UNIMOD:NN` accession or a bare modification name instead of a numeric mass
//! delta.
//!
//! # Mass provenance
//!
//! **No modification masses are hardcoded.** Every resolved mass comes from the
//! `delta_mono_mass` Unimod cross-reference inside the OBO term, e.g.
//!
//! ```text
//! [Term]
//! id: UNIMOD:35
//! name: Oxidation
//! xref: delta_mono_mass "15.994915"
//! ```
//!
//! The loader is lazy and cached behind a [`OnceLock`]; the bundled path is
//! resolved the same way the model store is (binary-relative `resources/`,
//! falling back to the compile-time source tree).

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::OnceLock;

/// A resolved Unimod modification record.
#[derive(Debug, Clone, PartialEq)]
pub struct UnimodMod {
    /// Canonical accession, e.g. `UNIMOD:35`.
    pub accession: String,
    /// Canonical Unimod name, e.g. `Oxidation`.
    pub name: String,
    /// Monoisotopic mass delta in Da, taken verbatim from the OBO
    /// `delta_mono_mass` cross-reference.
    pub mono_mass: f64,
    /// Residue sites this modification targets (e.g. `["M"]`, `["K", "N-term"]`),
    /// collected from the `spec_<n>_site` cross-references.
    pub sites: Vec<String>,
}

/// Errors returned when resolving a token fails.
#[derive(thiserror::Error, Debug, PartialEq)]
pub enum UnimodResolveError {
    /// A `UNIMOD:NN` accession that is not present in the OBO.
    #[error("unknown Unimod accession {token:?}")]
    UnknownUnimodAccession { token: String },
    /// A bare name that maps to more than one accession (case-insensitive).
    #[error("ambiguous modification name {token:?} (matches {accessions:?})")]
    AmbiguousModName {
        token: String,
        accessions: Vec<String>,
    },
}

/// Parsed Unimod database: accession and name indexes.
#[derive(Debug, Default)]
pub struct UnimodDb {
    by_accession: HashMap<String, UnimodMod>,
    /// Lowercased name → list of accessions (a name *could* be non-unique).
    by_name: HashMap<String, Vec<String>>,
}

impl UnimodDb {
    /// Parse an OBO document from a string. Hand-rolled line parser — no OBO
    /// dependency. Recognises `[Term]` stanzas with `id:`, `name:`, and the
    /// `xref: delta_mono_mass "<f64>"` / `xref: spec_<n>_site "<site>"`
    /// cross-references. Terms without both an `id`/`name` and a parseable
    /// `delta_mono_mass` are skipped.
    pub fn parse_obo(text: &str) -> Self {
        let mut db = UnimodDb::default();

        let mut in_term = false;
        let mut cur_id: Option<String> = None;
        let mut cur_name: Option<String> = None;
        let mut cur_mass: Option<f64> = None;
        let mut cur_sites: Vec<String> = Vec::new();

        // Commit the in-progress term (if complete) into the db.
        let commit = |db: &mut UnimodDb,
                      id: &mut Option<String>,
                      name: &mut Option<String>,
                      mass: &mut Option<f64>,
                      sites: &mut Vec<String>| {
            if let (Some(id_v), Some(name_v), Some(mass_v)) =
                (id.take(), name.take(), mass.take())
            {
                let entry = UnimodMod {
                    accession: id_v.clone(),
                    name: name_v.clone(),
                    mono_mass: mass_v,
                    sites: std::mem::take(sites),
                };
                db.by_name
                    .entry(name_v.to_ascii_lowercase())
                    .or_default()
                    .push(id_v.clone());
                db.by_accession.insert(id_v, entry);
            } else {
                // Incomplete term — drop any partial site list.
                sites.clear();
            }
        };

        for raw in text.lines() {
            let line = raw.trim();
            if line == "[Term]" {
                // Starting a new term: commit the previous one.
                commit(
                    &mut db,
                    &mut cur_id,
                    &mut cur_name,
                    &mut cur_mass,
                    &mut cur_sites,
                );
                in_term = true;
                continue;
            }
            if line.starts_with('[') && line.ends_with(']') {
                // A different stanza type (e.g. `[Typedef]`) — commit & leave.
                commit(
                    &mut db,
                    &mut cur_id,
                    &mut cur_name,
                    &mut cur_mass,
                    &mut cur_sites,
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
                    if !cur_sites.contains(&site) {
                        cur_sites.push(site);
                    }
                }
            }
        }
        // Commit the final term at EOF.
        commit(
            &mut db,
            &mut cur_id,
            &mut cur_name,
            &mut cur_mass,
            &mut cur_sites,
        );

        db
    }

    /// Resolve a token: either `UNIMOD:NN` (exact, case-insensitive prefix) or
    /// a bare name (case-insensitive). Returns `None` if the token is neither
    /// form found; returns `Err` for a known-but-unresolvable token (unknown
    /// accession, ambiguous name).
    pub fn resolve(&self, token: &str) -> Result<Option<&UnimodMod>, UnimodResolveError> {
        let t = token.trim();
        if let Some(num) = t
            .strip_prefix("UNIMOD:")
            .or_else(|| t.strip_prefix("unimod:"))
            .or_else(|| t.strip_prefix("Unimod:"))
        {
            // Canonicalise to `UNIMOD:<num>`.
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

    /// Look up a record by exact `UNIMOD:NN` accession.
    pub fn by_accession(&self, accession: &str) -> Option<&UnimodMod> {
        self.by_accession.get(accession.trim())
    }

    /// Number of parsed terms (test/diagnostic helper).
    pub fn len(&self) -> usize {
        self.by_accession.len()
    }

    /// Whether the db parsed no usable terms.
    pub fn is_empty(&self) -> bool {
        self.by_accession.is_empty()
    }
}

/// Parse a `"<float>"` (possibly with surrounding text after) into an f64.
/// Accepts the OBO form `delta_mono_mass "15.994915"`.
fn parse_quoted_f64(s: &str) -> Option<f64> {
    let start = s.find('"')? + 1;
    let end = s[start..].find('"')? + start;
    s[start..end].trim().parse().ok()
}

/// Parse `spec_<n>_site "<site>"` → the site string. Returns `None` for any
/// other `spec_*` or non-spec xref.
fn parse_spec_site(s: &str) -> Option<String> {
    let rest = s.strip_prefix("spec_")?;
    // rest looks like `1_site "K"`. Require the `_site` segment.
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

/// Resolve the bundled `unimod.obo` path the same way the model store does:
/// prefer `<exe_dir>/resources/unimod.obo` (packaged release next to the
/// binary), fall back to the compile-time source tree.
fn bundled_obo_path() -> PathBuf {
    let mut roots: Vec<PathBuf> = Vec::new();
    if let Ok(exe) = std::env::current_exe() {
        if let Some(dir) = exe.parent() {
            roots.push(dir.join("resources"));
        }
    }
    roots.push(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../resources"));

    for root in &roots {
        let p = root.join("unimod.obo");
        if p.exists() {
            return p;
        }
    }
    // Last-resort default (source tree) for error messages.
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../resources/unimod.obo")
}

/// Process-wide cached parse of the bundled OBO. Returns an empty db if the
/// file cannot be read (callers treat "unresolved" as an explicit error, so a
/// missing OBO surfaces as `UnknownUnimodAccession`/unresolved rather than a
/// panic).
pub fn bundled() -> &'static UnimodDb {
    static DB: OnceLock<UnimodDb> = OnceLock::new();
    DB.get_or_init(|| {
        let path = bundled_obo_path();
        match std::fs::read_to_string(&path) {
            Ok(text) => UnimodDb::parse_obo(&text),
            Err(_) => UnimodDb::default(),
        }
    })
}

/// Convenience: resolve a token against the bundled OBO.
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
        // Reference values pinned to the same Unimod numbers as
        // tests/common_mod_masses.rs. Production masses come from the OBO.
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
            // Bare name resolves to the same record.
            let by_name = db.resolve(name).unwrap().unwrap();
            assert_eq!(by_name.accession, *acc, "name {name} resolved elsewhere");
        }
    }
}
