//! The crosslinker of a search, from three sources like the modification
//! lists: `--xl` on the command line, `search.xl` in the config file, and the
//! SDRF column `comment[cross-linker]`.
//!
//! The command-line spec is written like a modification list item,
//! `Name (sites)`:
//!
//! - **Name**: a reagent of the reagent table (`resources/crosslinkers.yaml`,
//!   [`model::xl_reagents`]) by name or alias (`DSSO`, `BS3`, `DSBU`, `PhoX`,
//!   any case), its XLMOD accession (`XLMOD:02126`), or the UNIMOD name or
//!   accession of its intact cross-link record (`Xlink:DSS[138]`,
//!   `UNIMOD:1898`).
//! - **Sites** (optional; the reagent's usual reactivity when omitted):
//!   comma-separated residues and termini, each in the modification site
//!   grammar (`K, S, T, Y, Protein N-term`; `KSTY` also works). A
//!   heterobifunctional reagent gives its two ends apart with `/`:
//!   `EDC (K, Protein N-term / D, E, Protein C-term)`. `N-term` is the
//!   N-terminus of every peptide, for synthetic peptide libraries or samples
//!   linked after digestion (a protein linked before digestion has only its
//!   protein N-terminus free).
//!
//! The config `search.xl` is the same string, or a reagent entry of the
//! table's schema ([`model::xl_reagents::Reagent`]) for a reagent the table
//! lacks. The SDRF cell is read by [`from_sdrf`].

use model::modification::{ModLocation, ResidueSpec};
use model::unimod::UnimodDb;
use model::xl_reagents::{Reagent, ReagentTable};

use crate::crosslinker::{Cleavage, Crosslinker, LinkEnd, Monolink, Stub, MAX_STUBS};

/// A command-line crosslinker spec before it is resolved.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct XlDecl {
    /// Reagent name, alias, XLMOD accession or UNIMOD record.
    pub name: String,
    /// Site spec; the reagent's default when `None`.
    pub sites: Option<String>,
}

/// Why a crosslinker spec was rejected. Each message names the spec.
#[derive(thiserror::Error, Debug, PartialEq)]
pub enum XlSpecError {
    #[error("`{spec}`: unknown crosslinker `{name}`{hint}")]
    Unknown {
        spec: String,
        name: String,
        hint: String,
    },
    #[error("`{spec}`: `{name}` is a {role} of {reagent}; give `{reagent}` as the crosslinker{}",
            monolink_hint(.role))]
    NotALinker {
        spec: String,
        name: String,
        role: &'static str,
        reagent: String,
    },
    #[error("`{spec}`: {why}")]
    BadMass { spec: String, why: String },
    #[error("`{spec}`: no crosslinker name")]
    MissingName { spec: String },
    #[error("`{spec}`: unbalanced parentheses")]
    Unbalanced { spec: String },
    #[error("`{spec}`: bad sites `{sites}`: {why}")]
    BadSites {
        spec: String,
        sites: String,
        why: String,
    },
    #[error("`{spec}`: a custom reagent needs `sites`")]
    CustomWithoutSites { spec: String },
    #[error("`{spec}`: {why}")]
    BadCleavage { spec: String, why: String },
    #[error("`{spec}`: {why}")]
    Sdrf { spec: String, why: String },
}

fn monolink_hint(role: &str) -> &'static str {
    if role == "monolink" {
        " (and the monolink to --variable-mods)"
    } else {
        ""
    }
}

impl XlDecl {
    /// Parse `Name` or `Name (sites)`. The sites are the last parenthesised
    /// group.
    pub fn parse(spec: &str) -> Result<Self, XlSpecError> {
        let s = spec.trim();
        let err = || s.to_string();
        let (name, sites) = if s.ends_with(')') {
            let mut depth = 0usize;
            let mut open = None;
            for (i, c) in s.char_indices().rev() {
                match c {
                    ')' => depth += 1,
                    '(' => {
                        depth = depth
                            .checked_sub(1)
                            .ok_or(XlSpecError::Unbalanced { spec: err() })?;
                        if depth == 0 {
                            open = Some(i);
                            break;
                        }
                    }
                    _ => {}
                }
            }
            let open = open.ok_or(XlSpecError::Unbalanced { spec: err() })?;
            (s[..open].trim(), Some(s[open + 1..s.len() - 1].trim()))
        } else {
            (s, None)
        };
        if name.is_empty() {
            return Err(XlSpecError::MissingName { spec: err() });
        }
        if let Some(sites) = sites {
            parse_ends(s, sites)?;
        }
        Ok(XlDecl {
            name: name.to_string(),
            sites: sites.map(str::to_string),
        })
    }

    /// How the spec reads in messages.
    pub fn label(&self) -> String {
        match &self.sites {
            Some(s) => format!("{} ({s})", self.name),
            None => self.name.clone(),
        }
    }

    /// The crosslinker of the bundled reagent table this spec names.
    pub fn resolve(&self, db: &UnimodDb) -> Result<Crosslinker, XlSpecError> {
        self.resolve_in(model::xl_reagents::bundled(), db)
    }

    /// [`Self::resolve`] against a given reagent table.
    pub fn resolve_in(
        &self,
        table: &ReagentTable,
        db: &UnimodDb,
    ) -> Result<Crosslinker, XlSpecError> {
        let spec = self.label();
        let reagent = lookup(table, db, &self.name, &spec)?;
        from_reagent(reagent, db, &spec, self.sites.as_deref())
    }
}

/// The crosslinker of a custom reagent entry (the config file's `search.xl`
/// map): the table's schema, with `sites` required.
pub fn custom(r: &Reagent, db: &UnimodDb) -> Result<Crosslinker, XlSpecError> {
    if r.sites.trim().is_empty() {
        return Err(XlSpecError::CustomWithoutSites {
            spec: r.name.clone(),
        });
    }
    from_reagent(r, db, &r.name, None)
}

/// The reagent `name` stands for: the table by name, alias or intact UNIMOD
/// record; then by XLMOD accession; then a UNIMOD accession or name whose
/// record is a reagent's intact cross-link.
fn lookup<'t>(
    table: &'t ReagentTable,
    db: &UnimodDb,
    name: &str,
    spec: &str,
) -> Result<&'t Reagent, XlSpecError> {
    if let Some(r) = table.find(name).or_else(|| table.find_xlmod(name)) {
        return Ok(r);
    }
    if let Ok(Some(rec)) = db.resolve(name) {
        if let Some(r) = table.find(&rec.name) {
            return Ok(r);
        }
        if let Some((r, role)) = table.owner_of_record(&rec.name) {
            return Err(XlSpecError::NotALinker {
                spec: spec.to_string(),
                name: rec.name.clone(),
                role,
                reagent: r.name.clone(),
            });
        }
    }
    Err(XlSpecError::Unknown {
        spec: spec.to_string(),
        name: name.trim().to_string(),
        hint: unknown_hint(table, name),
    })
}

fn unknown_hint(table: &ReagentTable, name: &str) -> String {
    if let Some(e) = model::xl_reagents::bundled_error() {
        return format!(" (the reagent table could not be loaded: {e})");
    }
    let s = suggest(table, name);
    if s.is_empty() {
        format!(" (known: {})", table.names().collect::<Vec<_>>().join(", "))
    } else {
        format!("; did you mean {}?", s.join(", "))
    }
}

/// Up to three reagent names within a small edit distance of `token`,
/// closest first; on a tie, the name of the token's length first.
fn suggest(table: &ReagentTable, token: &str) -> Vec<String> {
    let t = token.trim().to_ascii_lowercase();
    let mut scored: Vec<(usize, usize, &str)> = table
        .names()
        .map(|n| {
            let d = edit_distance(&t, &n.to_ascii_lowercase());
            (d, n.len().abs_diff(t.len()), n)
        })
        .filter(|(d, _, _)| *d <= 2)
        .collect();
    scored.sort();
    scored
        .into_iter()
        .take(3)
        .map(|(_, _, n)| n.to_string())
        .collect()
}

fn edit_distance(a: &str, b: &str) -> usize {
    let a: Vec<char> = a.chars().collect();
    let b: Vec<char> = b.chars().collect();
    let mut d = vec![vec![0usize; b.len() + 1]; a.len() + 1];
    for (i, row) in d.iter_mut().enumerate() {
        row[0] = i;
    }
    for (j, v) in d[0].iter_mut().enumerate() {
        *v = j;
    }
    for i in 1..=a.len() {
        for j in 1..=b.len() {
            let cost = usize::from(a[i - 1] != b[j - 1]);
            let mut v = (d[i - 1][j] + 1)
                .min(d[i][j - 1] + 1)
                .min(d[i - 1][j - 1] + cost);
            if i > 1 && j > 1 && a[i - 1] == b[j - 2] && a[i - 2] == b[j - 1] {
                v = v.min(d[i - 2][j - 2] + 1);
            }
            d[i][j] = v;
        }
    }
    d[a.len()][b.len()]
}

/// The crosslinker of reagent entry `r`, with `sites` in place of its default
/// sites when given.
pub fn from_reagent(
    r: &Reagent,
    db: &UnimodDb,
    spec: &str,
    sites: Option<&str>,
) -> Result<Crosslinker, XlSpecError> {
    let mass = |p: &model::xl_reagents::Part| {
        p.resolve(db).map_err(|why| XlSpecError::BadMass {
            spec: spec.to_string(),
            why: format!("{}: {why}", r.name),
        })
    };
    let bad = |why: String| XlSpecError::BadCleavage {
        spec: spec.to_string(),
        why,
    };
    let (intact, accession) = mass(&r.intact)?;
    let cleavage = if r.stubs.is_empty() {
        if !r.doublets.is_empty() {
            return Err(bad("`doublets` without `stubs`".into()));
        }
        None
    } else {
        let mut stubs: Vec<Stub> = Vec::new();
        for s in &r.stubs {
            if stubs.iter().any(|t| t.name == s.name) {
                return Err(bad(format!("stub `{}` given twice", s.name)));
            }
            stubs.push(Stub {
                name: s.name.clone(),
                mass: mass(s)?.0,
            });
        }
        let index = |n: &str| {
            stubs
                .iter()
                .position(|s| s.name == n)
                .ok_or_else(|| bad(format!("doublet stub `{n}` is not one of the stubs")))
        };
        let mut doublets = Vec::new();
        for [a, b] in &r.doublets {
            doublets.push((index(a)?, index(b)?));
        }
        let c = Cleavage { stubs, doublets };
        check_cleavage(&c).map_err(bad)?;
        Some(order_doublets(c))
    };
    let mut monolinks = Vec::new();
    for m in &r.monolinks {
        let (mass, accession) = mass(m)?;
        monolinks.push(Monolink {
            name: m
                .unimod
                .clone()
                .unwrap_or_else(|| format!("{}:{}", r.name, m.name)),
            accession,
            mass,
        });
    }
    Ok(Crosslinker {
        name: r.name.clone(),
        accession,
        intact,
        cleavage,
        ends: parse_ends(spec, sites.unwrap_or(&r.sites))?,
        monolinks,
    })
}

/// Doublets with the lighter stub first.
fn order_doublets(mut c: Cleavage) -> Cleavage {
    for d in &mut c.doublets {
        if c.stubs[d.0].mass > c.stubs[d.1].mass {
            *d = (d.1, d.0);
        }
    }
    c
}

fn check_cleavage(c: &Cleavage) -> Result<(), String> {
    if c.stubs.len() > MAX_STUBS {
        return Err(format!("at most {MAX_STUBS} stubs"));
    }
    if c.doublets.is_empty() {
        return Err(
            "an MS-cleavable linker needs at least one reporter doublet (`doublets`)".into(),
        );
    }
    for &(a, b) in &c.doublets {
        if (c.stubs[a].mass - c.stubs[b].mass).abs() < 1e-3 {
            return Err(format!(
                "doublet `{}`/`{}` has no mass spacing",
                c.stubs[a].name, c.stubs[b].name
            ));
        }
    }
    Ok(())
}

/// Tolerance of the SDRF stub masses `MH` / `ML` against the reagent's stubs,
/// in Da: the template writes them to two decimals.
pub const SDRF_STUB_TOL_DA: f64 = 0.02;

/// The crosslinker of an SDRF `comment[cross-linker]` cell, e.g.
/// `NT=DSSO;AC=XLMOD:02126;CL=yes;TA=K,S,T,Y,nterm;MH=54.01;ML=85.98`.
///
/// - The reagent is found by `AC` (XLMOD, or a UNIMOD intact record) first,
///   then by `NT`; when both resolve to different reagents it is an error.
/// - `TA` replaces the sites. `nterm` is the protein N-terminus: the SDRF
///   describes the reagent on the sample, and a protein linked before
///   digestion has only its protein N-terminal amine free (the peptide
///   N-termini appear at digestion, after the linking). `cterm` is the
///   protein C-terminus likewise. On a heterobifunctional reagent each end
///   keeps the `TA` targets it reacts with.
/// - `CL`, `MH` and `ML` are checked against the reagent: `CL` against
///   whether it has stubs, `MH` / `ML` each against a stub mass within
///   [`SDRF_STUB_TOL_DA`] (the template's names do not say which is the light
///   stub, so either order is accepted).
/// - A reagent the table lacks is a custom MS-cleavable reagent when `TA`,
///   `MH` and `ML` are all given: the two stubs are one reporter doublet and
///   their sum is the intact cross-link (the split of a urea such as DSBU or
///   CDI); any other unknown reagent is an error that names the fix.
pub fn from_sdrf(value: &str, db: &UnimodDb) -> Result<Crosslinker, XlSpecError> {
    from_sdrf_in(model::xl_reagents::bundled(), value, db)
}

/// [`from_sdrf`] against a given reagent table.
pub fn from_sdrf_in(
    table: &ReagentTable,
    value: &str,
    db: &UnimodDb,
) -> Result<Crosslinker, XlSpecError> {
    let spec = value.trim().to_string();
    let err = |why: String| XlSpecError::Sdrf {
        spec: spec.clone(),
        why,
    };
    let mut kv: std::collections::BTreeMap<String, String> = Default::default();
    for part in value.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        match part.split_once('=') {
            Some((k, v)) => {
                kv.insert(k.trim().to_ascii_uppercase(), v.trim().to_string());
            }
            None => {
                kv.entry("NT".into()).or_insert_with(|| part.to_string());
            }
        }
    }
    let num = |key: &str| -> Result<Option<f64>, XlSpecError> {
        kv.get(key)
            .map(|v| {
                v.parse::<f64>()
                    .map_err(|_| err(format!("{key}={v} is not a number")))
            })
            .transpose()
    };
    let (mh, ml) = (num("MH")?, num("ML")?);
    let by_ac = kv
        .get("AC")
        .and_then(|ac| lookup(table, db, ac, &spec).ok());
    let by_nt = kv
        .get("NT")
        .and_then(|nt| lookup(table, db, nt, &spec).ok());
    let sites = kv
        .get("TA")
        .map(|ta| ta_sites(ta))
        .transpose()
        .map_err(err)?;
    let mut x = match (by_ac, by_nt) {
        (Some(a), Some(n)) if a.name != n.name => {
            return Err(err(format!(
                "AC={} is {} but NT={} is {}; fix one of them",
                kv["AC"], a.name, kv["NT"], n.name
            )))
        }
        (Some(r), _) | (None, Some(r)) => from_reagent(r, db, &spec, None)?,
        (None, None) => {
            let (Some(ta), Some(mh), Some(ml)) = (kv.get("TA"), mh, ml) else {
                return Err(err(format!(
                    "unknown crosslinker (NT={}, AC={}); give a known name or XLMOD accession, \
                     or define it with TA, MH and ML (an MS-cleavable reagent whose two stubs \
                     add up to the cross-link), or define it in the config file (search.xl)",
                    kv.get("NT").map_or("-", String::as_str),
                    kv.get("AC").map_or("-", String::as_str)
                )));
            };
            let (light, heavy) = if mh <= ml { (mh, ml) } else { (ml, mh) };
            let r = Reagent {
                name: kv
                    .get("NT")
                    .cloned()
                    .unwrap_or_else(|| "custom".to_string()),
                xlmod: kv.get("AC").cloned().into_iter().collect(),
                sites: ta_sites(ta).map_err(err)?,
                intact: model::xl_reagents::Part {
                    mass: Some(light + heavy),
                    ..Default::default()
                },
                stubs: vec![
                    model::xl_reagents::Part {
                        name: "light".into(),
                        mass: Some(light),
                        ..Default::default()
                    },
                    model::xl_reagents::Part {
                        name: "heavy".into(),
                        mass: Some(heavy),
                        ..Default::default()
                    },
                ],
                doublets: vec![["light".into(), "heavy".into()]],
                ..Default::default()
            };
            return from_reagent(&r, db, &spec, None);
        }
    };
    if let Some(cl) = kv.get("CL") {
        let says = match cl.to_ascii_lowercase().as_str() {
            "yes" | "true" | "1" => true,
            "no" | "false" | "0" => false,
            other => return Err(err(format!("CL={other}: expected yes or no"))),
        };
        if says != x.cleavage.is_some() {
            return Err(err(format!(
                "CL={cl} but {} is {}MS-cleavable",
                x.name,
                if x.cleavage.is_some() { "" } else { "not " }
            )));
        }
    }
    for (key, v) in [("MH", mh), ("ML", ml)] {
        let Some(v) = v else { continue };
        let stubs: Vec<f64> = x
            .cleavage
            .as_ref()
            .map(|c| c.stub_masses().collect())
            .unwrap_or_default();
        if !stubs.iter().any(|s| (s - v).abs() <= SDRF_STUB_TOL_DA) {
            return Err(err(format!(
                "{key}={v} matches no cleavage stub of {} ({}) within {SDRF_STUB_TOL_DA} Da",
                x.name,
                if stubs.is_empty() {
                    "it has none".to_string()
                } else {
                    stubs
                        .iter()
                        .map(|s| format!("{s:.4}"))
                        .collect::<Vec<_>>()
                        .join(", ")
                }
            )));
        }
    }
    if let Some(sites) = sites {
        let ta = parse_ends(&spec, &sites)?[0].clone();
        if x.is_homobifunctional() {
            x.ends = [ta.clone(), ta];
        } else {
            let keep = |e: &LinkEnd| LinkEnd {
                residues: e
                    .residues
                    .iter()
                    .copied()
                    .filter(|r| ta.residues.contains(r))
                    .collect(),
                protein_n_term: e.protein_n_term && ta.protein_n_term,
                protein_c_term: e.protein_c_term && ta.protein_c_term,
                peptide_n_term: false,
            };
            let ends = [keep(&x.ends[0]), keep(&x.ends[1])];
            let covered = |r: u8| ends.iter().any(|e| e.residues.contains(&r));
            if let Some(&r) = ta.residues.iter().find(|&&r| !covered(r)) {
                return Err(err(format!(
                    "TA target {} is not a target of either end of {} ({})",
                    r as char,
                    x.name,
                    x.describe_sites()
                )));
            }
            if ends
                .iter()
                .any(|e| e.residues.is_empty() && !e.protein_n_term && !e.protein_c_term)
            {
                return Err(err(format!(
                    "TA leaves an end of {} ({}) without targets",
                    x.name,
                    x.describe_sites()
                )));
            }
            x.ends = ends;
        }
    }
    Ok(x)
}

/// An SDRF `TA` list (`K,S,T,Y,nterm`) as a site spec (`K, S, T, Y, Protein
/// N-term`): `nterm` / `cterm` are the protein termini (see [`from_sdrf`]).
fn ta_sites(ta: &str) -> Result<String, String> {
    let mut out = Vec::new();
    for t in ta.split(',').map(str::trim).filter(|t| !t.is_empty()) {
        match t.to_ascii_lowercase().as_str() {
            "nterm" | "n-term" | "protein n-term" => out.push("Protein N-term".to_string()),
            "cterm" | "c-term" | "protein c-term" => out.push("Protein C-term".to_string()),
            _ if t.len() == 1 && t.as_bytes()[0].is_ascii_alphabetic() => {
                out.push(t.to_ascii_uppercase())
            }
            _ => {
                return Err(format!(
                    "TA target `{t}`: expected one-letter residues, nterm or cterm"
                ))
            }
        }
    }
    if out.is_empty() {
        return Err("TA lists no targets".into());
    }
    Ok(out.join(", "))
}

/// The two ends of a site spec: one end for both, or `end / end`.
pub fn parse_ends(spec: &str, sites: &str) -> Result<[LinkEnd; 2], XlSpecError> {
    let parts: Vec<&str> = sites.split('/').collect();
    match parts.as_slice() {
        [one] => {
            let e = parse_end(spec, sites, one)?;
            Ok([e.clone(), e])
        }
        [a, b] => Ok([parse_end(spec, sites, a)?, parse_end(spec, sites, b)?]),
        _ => Err(XlSpecError::BadSites {
            spec: spec.to_string(),
            sites: sites.to_string(),
            why: "at most two ends, separated by one `/`".into(),
        }),
    }
}

fn parse_end(spec: &str, sites: &str, end: &str) -> Result<LinkEnd, XlSpecError> {
    let bad = |why: String| XlSpecError::BadSites {
        spec: spec.to_string(),
        sites: sites.to_string(),
        why,
    };
    let mut out = LinkEnd::default();
    for item in end.split(',') {
        let item = item.trim();
        if item.is_empty() {
            return Err(bad("an empty site".into()));
        }
        let parsed = model::mod_spec::parse_sites(spec, item).map_err(|e| bad(e.to_string()))?;
        for (residue, location) in parsed {
            match (residue, location) {
                (ResidueSpec::Specific(r), ModLocation::Anywhere) => {
                    if !out.residues.contains(&r) {
                        out.residues.push(r);
                    }
                }
                (ResidueSpec::Wildcard, ModLocation::ProtNTerm) => out.protein_n_term = true,
                (ResidueSpec::Wildcard, ModLocation::ProtCTerm) => out.protein_c_term = true,
                (ResidueSpec::Wildcard, ModLocation::NTerm) => out.peptide_n_term = true,
                (ResidueSpec::Wildcard, ModLocation::CTerm) => {
                    return Err(bad(format!(
                        "`{item}`: the peptide C-terminus is not a link site; write Protein \
                         C-term"
                    )))
                }
                _ => {
                    return Err(bad(format!(
                        "`{item}`: a terminus takes no residue; write the residue and the \
                         terminus apart"
                    )))
                }
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn db() -> &'static UnimodDb {
        model::unimod::bundled()
    }

    fn resolve(s: &str) -> Result<Crosslinker, XlSpecError> {
        XlDecl::parse(s).and_then(|d| d.resolve(db()))
    }

    #[test]
    fn every_table_reagent_resolves_with_its_default_sites() {
        for r in &model::xl_reagents::bundled().reagents {
            let x = resolve(&r.name).unwrap_or_else(|e| panic!("{}: {e}", r.name));
            assert_eq!(x.name, r.name);
            assert_eq!(x.cleavage.is_some(), !r.stubs.is_empty());
            assert_eq!(x.monolinks.len(), r.monolinks.len());
            for a in &r.aliases {
                assert_eq!(resolve(a).unwrap(), x, "{a}");
            }
            for a in &r.xlmod {
                assert_eq!(resolve(a).unwrap(), x, "{a}");
            }
            if let Some(u) = &r.intact.unimod {
                assert_eq!(resolve(u).unwrap(), x, "{u}");
                let acc = x.accession.clone().unwrap();
                assert_eq!(resolve(&acc).unwrap(), x, "{acc}");
            }
        }
    }

    #[test]
    fn sites_override_the_default() {
        let x = resolve("DSSO (K, S, T, Y, Protein N-term)").unwrap();
        assert_eq!(x.ends[0].residues, b"KSTY");
        assert!(x.ends[0].protein_n_term && x.is_homobifunctional());
        let y = resolve("dsso (KSTY, Protein N-term)").unwrap();
        assert_eq!(x, y);
        assert_eq!(x.describe_sites(), "K, S, T, Y, Protein N-term");
        let z = resolve("DSS (K)").unwrap();
        assert!(!z.ends[0].protein_n_term);
        let n = resolve("DSBU (K, N-term)").unwrap();
        assert!(n.ends[0].peptide_n_term && !n.ends[0].protein_n_term);
        assert_eq!(n.describe_sites(), "K, N-term");
    }

    #[test]
    fn heterobifunctional_ends() {
        let x = resolve("EDC").unwrap();
        assert!(!x.is_homobifunctional());
        assert_eq!(x.ends[0].residues, b"K");
        assert!(x.ends[0].protein_n_term);
        assert_eq!(x.ends[1].residues, b"DE");
        assert!(x.ends[1].protein_c_term);
        assert!((x.intact + 18.010_565).abs() < 1e-5);
        assert_eq!(
            x.describe_sites(),
            "K, Protein N-term / D, E, Protein C-term"
        );
    }

    #[test]
    fn bad_specs_name_the_problem() {
        let e = resolve("DSOS").unwrap_err().to_string();
        assert!(e.contains("did you mean DSSO"), "{e}");
        let e = resolve("Xlink:DSSO[176]").unwrap_err().to_string();
        assert!(e.contains("monolink of DSSO"), "{e}");
        assert!(e.contains("--variable-mods"), "{e}");
        let e = resolve("Xlink:DSSO[54]").unwrap_err().to_string();
        assert!(e.contains("cleavage stub of DSSO"), "{e}");
        let e = resolve("EDC (K / D, C-term)").unwrap_err().to_string();
        assert!(e.contains("C-terminus is not a link site"), "{e}");
        let e = resolve("DSSO (N-term K)").unwrap_err().to_string();
        assert!(e.contains("takes no residue"), "{e}");
        let e = resolve("DSSO (Protein N-term K)").unwrap_err().to_string();
        assert!(e.contains("takes no residue"), "{e}");
        let e = resolve("DSSO (K / D / E)").unwrap_err().to_string();
        assert!(e.contains("at most two ends"), "{e}");
        assert!(matches!(
            resolve(" (K)"),
            Err(XlSpecError::MissingName { .. })
        ));
        let e = resolve("DSSO (k)").unwrap_err().to_string();
        assert!(e.contains("bad sites"), "{e}");
    }

    fn entry(yaml: &str) -> Reagent {
        serde_yaml::from_str(yaml).unwrap()
    }

    #[test]
    fn a_custom_entry_has_the_table_schema() {
        let r = entry(
            "name: MyDSSO\nsites: K, Protein N-term\nintact: {mass: 158.003765}\n\
             stubs: [{name: thiol, mass: 85.982635}, {name: alkene, composition: H(2) C(3) O}, \
             {name: sulfenic, mass: 103.9932}]\ndoublets: [[thiol, alkene]]\n\
             monolinks: [{name: water, mass: 176.01433}]\n",
        );
        let x = custom(&r, db()).unwrap();
        let c = x.cleavage.as_ref().unwrap();
        // The doublet is ordered light (alkene) then heavy (thiol).
        assert_eq!(c.doublets, vec![(1, 0)]);
        assert!((c.doublet_delta(0) - 31.972_07).abs() < 1e-4);
        assert_eq!(x.monolinks[0].name, "MyDSSO:water");
        let no_doublets = Reagent {
            doublets: vec![],
            ..r.clone()
        };
        assert!(custom(&no_doublets, db())
            .unwrap_err()
            .to_string()
            .contains("reporter doublet"));
        let no_sites = Reagent {
            sites: String::new(),
            ..r.clone()
        };
        assert!(matches!(
            custom(&no_sites, db()),
            Err(XlSpecError::CustomWithoutSites { .. })
        ));
        let bad_stub = Reagent {
            doublets: vec![["alkene".into(), "thiolx".into()]],
            ..r.clone()
        };
        assert!(custom(&bad_stub, db()).is_err());
        let zl = entry("name: ZL\nsites: K / D, E\nintact: {composition: H(-2) O(-1)}\n");
        let x = custom(&zl, db()).unwrap();
        assert!(x.cleavage.is_none() && !x.is_homobifunctional());
    }

    #[test]
    fn transpositions_are_one_edit() {
        assert_eq!(edit_distance("dsos", "dsso"), 1);
        assert_eq!(edit_distance("dss", "dsso"), 1);
        assert_eq!(edit_distance("", "abc"), 3);
    }

    fn sdrf(v: &str) -> Result<Crosslinker, XlSpecError> {
        from_sdrf(v, db())
    }

    #[test]
    fn sdrf_resolves_by_accession_then_name() {
        let x = sdrf("NT=DSS;AC=XLMOD:02001").unwrap();
        assert_eq!(x.name, "DSS");
        assert_eq!(sdrf("NT=BS3;AC=XLMOD:02000").unwrap().name, "DSS");
        // An accession the table lacks falls back to the name.
        assert_eq!(sdrf("NT=DSSO;AC=XLMOD:99999").unwrap().name, "DSSO");
        assert_eq!(sdrf("AC=XLMOD:02126").unwrap().name, "DSSO");
        assert_eq!(sdrf("DSBU").unwrap().name, "DSBU");
        // Name and accession of two different reagents.
        let e = sdrf("NT=DSSO;AC=XLMOD:02010").unwrap_err().to_string();
        assert!(
            e.contains("AC=XLMOD:02010 is EDC but NT=DSSO is DSSO"),
            "{e}"
        );
    }

    #[test]
    fn sdrf_targets_replace_the_sites() {
        let x = sdrf("NT=DSSO;AC=XLMOD:02126;CL=yes;TA=K,S,T,Y,nterm;MH=54.01;ML=85.98").unwrap();
        assert_eq!(x.ends[0].residues, b"KSTY");
        // `nterm` is the protein N-terminus, not every peptide's.
        assert!(x.ends[0].protein_n_term && !x.ends[0].peptide_n_term);
        // Heterobifunctional: each end keeps the targets it reacts with.
        let e = sdrf("NT=EDC;AC=XLMOD:02010;CL=no;TA=K,D,E").unwrap();
        assert_eq!(e.ends[0].residues, b"K");
        assert!(!e.ends[0].protein_n_term);
        assert_eq!(e.ends[1].residues, b"DE");
        assert!(!e.ends[1].protein_c_term);
        let err = sdrf("NT=EDC;TA=K,D,E,Y").unwrap_err().to_string();
        assert!(err.contains("TA target Y"), "{err}");
        let err = sdrf("NT=EDC;TA=K").unwrap_err().to_string();
        assert!(err.contains("without targets"), "{err}");
        let err = sdrf("NT=DSS;TA=K,xyz").unwrap_err().to_string();
        assert!(err.contains("TA target `xyz`"), "{err}");
    }

    #[test]
    fn sdrf_cleavability_and_stub_masses_are_checked() {
        let e = sdrf("NT=DSS;CL=yes").unwrap_err().to_string();
        assert!(e.contains("CL=yes but DSS is not MS-cleavable"), "{e}");
        let e = sdrf("NT=DSSO;CL=no").unwrap_err().to_string();
        assert!(e.contains("CL=no but DSSO is MS-cleavable"), "{e}");
        // Either order of the two stub masses.
        assert!(sdrf("NT=DSSO;MH=85.98;ML=54.01").is_ok());
        let e = sdrf("NT=DSSO;MH=54.01;ML=86.10").unwrap_err().to_string();
        assert!(
            e.contains("ML=86.1 matches no cleavage stub of DSSO"),
            "{e}"
        );
        let e = sdrf("NT=DSS;MH=54.01").unwrap_err().to_string();
        assert!(e.contains("it has none"), "{e}");
    }

    #[test]
    fn sdrf_defines_an_unknown_cleavable_reagent_only_with_targets_and_stubs() {
        let x = sdrf("NT=NewUrea;CL=yes;TA=K,nterm;MH=111.03;ML=85.05").unwrap();
        assert_eq!(x.name, "NewUrea");
        assert!((x.intact - 196.08).abs() < 1e-9);
        let c = x.cleavage.unwrap();
        assert!((c.doublet_delta(0) - 25.98).abs() < 1e-9);
        assert_eq!(x.ends[0].residues, b"K");
        let e = sdrf("NT=NewUrea;MH=111.03;ML=85.05")
            .unwrap_err()
            .to_string();
        assert!(e.contains("define it with TA, MH and ML"), "{e}");
        assert!(sdrf("NT=NewUrea;TA=K").is_err());
    }
}
