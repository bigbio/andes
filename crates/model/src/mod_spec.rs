//! The modification lists of a search: `--fixed-mods` / `--variable-mods` on
//! the command line and `search.fixed_mods` / `search.variable_mods` in the
//! config file.
//!
//! A list is comma-separated items `Name (sites)`, e.g.
//! `Oxidation (M), Phospho (STY), Acetyl (Protein N-term), Gln->pyro-Glu (N-term Q)`.
//!
//! - **Name**: a UNIMOD name (`Phospho`), a UNIMOD accession (`UNIMOD:21`) or
//!   an alias from `resources/unimod_aliases.tsv` (`GlyGly`), any case. The
//!   mass and the accession come from the bundled UNIMOD table
//!   ([`crate::unimod`]); the record's UNIMOD name is the one written out.
//! - **Sites**: residues (`STY`); a terminus (`N-term`, `C-term`,
//!   `Protein N-term`, `Protein C-term`, or UNIMOD's `Any N-term`,
//!   `Any C-term`); or a terminus followed by residues (`N-term Q`: Q at the
//!   peptide N-terminus).
//!
//! A config item can also be a [`ModDecl`] map, which adds the attributes the
//! command line does not take (neutral losses, loss class, annotated-only
//! sites, cleavage blocking, diagnostic ions) and custom masses.

use crate::modification::{loss_class_id, ModLocation, Modification, ResidueSpec};
use crate::unimod::{UnimodDb, UnimodMod, UnimodResolveError};

/// One declared modification before it is resolved against UNIMOD and
/// expanded into one [`Modification`] per site.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ModDecl {
    /// UNIMOD name, accession or alias; with `mass`, a free label.
    pub name: Option<String>,
    /// Custom mass shift in Da; the name is then not looked up.
    pub mass: Option<f64>,
    /// Site spec, as inside the parentheses of `Name (sites)`.
    pub sites: String,
    /// Neutral-loss masses (Da) of the mod's fragment ions.
    pub neutral_loss: Vec<f64>,
    /// `glyco`, `phospho`, `sulfo` or `generic` (the default with losses).
    pub loss_class: Option<String>,
    /// `annotated` or `any` (the default).
    pub restrict_to: Option<String>,
    /// `blocked` or `kept` (the default).
    pub cleavage: Option<String>,
    /// m/z of singly charged ions the mod gives off when fragmented.
    pub diagnostic_ions: Vec<f64>,
}

/// Why a modification list or item was rejected. Each message names the item.
#[derive(thiserror::Error, Debug, PartialEq)]
pub enum ModSpecError {
    #[error("`{item}`: unknown modification `{name}`{}", suggestion(.suggestions, *.table_empty))]
    UnknownName {
        item: String,
        name: String,
        suggestions: Vec<String>,
        table_empty: bool,
    },
    #[error("`{item}`: unknown UNIMOD accession `{token}`")]
    UnknownAccession { item: String, token: String },
    #[error("`{item}`: `{token}` names several UNIMOD records ({}); give the accession", .accessions.join(", "))]
    AmbiguousName {
        item: String,
        token: String,
        accessions: Vec<String>,
    },
    #[error("`{item}`: no sites; write `Name (sites)`, e.g. `Oxidation (M)` or `Acetyl (Protein N-term)`")]
    MissingSites { item: String },
    #[error("`{item}`: no modification name before the sites")]
    MissingName { item: String },
    #[error("`{item}`: unbalanced parentheses")]
    Unbalanced { item: String },
    #[error(
        "`{item}`: bad sites `{sites}`; expected residues (STY), N-term, C-term, \
         Protein N-term, Protein C-term, or a terminus then residues (N-term Q)"
    )]
    BadSites { item: String, sites: String },
    #[error("an empty item in the modification list `{list}`")]
    EmptyItem { list: String },
    #[error("`{item}`: give a `name` or a `mass`")]
    NoNameOrMass { item: String },
    #[error("`{item}`: mass {mass} is not a finite number")]
    BadMass { item: String, mass: f64 },
    #[error("`{item}`: neutral loss {value} must be in (0, 5000) Da")]
    BadNeutralLoss { item: String, value: f64 },
    #[error("`{item}`: diagnostic ion {value} must be in (0, 5000) m/z")]
    BadDiagnosticIon { item: String, value: f64 },
    #[error("`{item}`: loss_class `{value}`; expected glyco, phospho, sulfo or generic")]
    BadLossClass { item: String, value: String },
    #[error("`{item}`: restrict_to `{value}`; expected annotated or any")]
    BadRestrictTo { item: String, value: String },
    #[error("`{item}`: cleavage `{value}`; expected blocked or kept")]
    BadCleavage { item: String, value: String },
    #[error("`{item}`: restrict_to: annotated needs a variable modification; a fixed one applies everywhere")]
    AnnotatedFixed { item: String },
}

fn suggestion(suggestions: &[String], table_empty: bool) -> String {
    if table_empty {
        " (the bundled UNIMOD table resources/unimod.obo could not be read)".to_string()
    } else if suggestions.is_empty() {
        " (not a UNIMOD name, accession or alias)".to_string()
    } else {
        format!("; did you mean {}?", suggestions.join(", "))
    }
}

impl ModDecl {
    /// Parse one `Name (sites)` item. The sites are the last parenthesised
    /// group, so names with parentheses work: `Label:13C(6) (K)`.
    pub fn parse(item: &str) -> Result<Self, ModSpecError> {
        let s = item.trim();
        let err_item = || s.to_string();
        if !s.ends_with(')') {
            return Err(ModSpecError::MissingSites { item: err_item() });
        }
        let mut depth = 0usize;
        let mut open = None;
        for (i, c) in s.char_indices().rev() {
            match c {
                ')' => depth += 1,
                '(' => {
                    depth = depth
                        .checked_sub(1)
                        .ok_or(ModSpecError::Unbalanced { item: err_item() })?;
                    if depth == 0 {
                        open = Some(i);
                        break;
                    }
                }
                _ => {}
            }
        }
        let open = open.ok_or(ModSpecError::Unbalanced { item: err_item() })?;
        let name = s[..open].trim();
        let sites = s[open + 1..s.len() - 1].trim();
        if name.is_empty() {
            return Err(ModSpecError::MissingName { item: err_item() });
        }
        parse_sites(s, sites)?;
        Ok(ModDecl {
            name: Some(name.to_string()),
            sites: sites.to_string(),
            ..ModDecl::default()
        })
    }

    /// How the item reads in messages: `Name (sites)`.
    pub fn label(&self) -> String {
        let name = match (&self.name, self.mass) {
            (Some(n), _) => n.clone(),
            (None, Some(m)) => format!("{m:+}"),
            (None, None) => String::new(),
        };
        format!("{name} ({})", self.sites)
    }

    /// The modifications this declaration stands for, one per site, in the
    /// order the sites are written. A UNIMOD name gives the mass, name and
    /// accession; a custom `mass` keeps the declared name (or `+<mass>`) and
    /// has no accession.
    pub fn resolve(&self, fixed: bool, db: &UnimodDb) -> Result<Vec<Modification>, ModSpecError> {
        let item = self.label();
        let (name, mass_delta, accession) = match (self.mass, &self.name) {
            (Some(mass), name) => {
                if !mass.is_finite() {
                    return Err(ModSpecError::BadMass { item, mass });
                }
                let name = name.clone().unwrap_or_else(|| format!("{mass:+}"));
                (name, mass, None)
            }
            (None, Some(name)) => {
                let u = lookup(db, name, &item)?;
                (u.name.clone(), u.mono_mass, Some(u.accession.clone()))
            }
            (None, None) => return Err(ModSpecError::NoNameOrMass { item }),
        };
        for &value in &self.neutral_loss {
            if !(value > 0.0 && value < 5000.0) {
                return Err(ModSpecError::BadNeutralLoss { item, value });
            }
        }
        for &value in &self.diagnostic_ions {
            if !(value > 0.0 && value < 5000.0) {
                return Err(ModSpecError::BadDiagnosticIon { item, value });
            }
        }
        let loss_class = match &self.loss_class {
            Some(v) => loss_class_id(v).ok_or_else(|| ModSpecError::BadLossClass {
                item: item.clone(),
                value: v.clone(),
            })?,
            // Losses without a class pool as generic.
            None if !self.neutral_loss.is_empty() => 255,
            None => 0,
        };
        let annotated = match self.restrict_to.as_deref().map(str::trim) {
            None => false,
            Some(v) if v.eq_ignore_ascii_case("any") => false,
            Some(v) if v.eq_ignore_ascii_case("annotated") => true,
            Some(v) => {
                return Err(ModSpecError::BadRestrictTo {
                    item,
                    value: v.to_string(),
                })
            }
        };
        if annotated && fixed {
            return Err(ModSpecError::AnnotatedFixed { item });
        }
        let blocks_cleavage = match self.cleavage.as_deref().map(str::trim) {
            None => false,
            Some(v) if v.eq_ignore_ascii_case("kept") => false,
            Some(v) if v.eq_ignore_ascii_case("blocked") => true,
            Some(v) => {
                return Err(ModSpecError::BadCleavage {
                    item,
                    value: v.to_string(),
                })
            }
        };
        let sites = parse_sites(&item, &self.sites)?;
        Ok(sites
            .into_iter()
            .map(|(residue, location)| Modification {
                accession: accession.clone(),
                neutral_losses: self.neutral_loss.clone(),
                loss_class,
                annotated_site_mass: annotated.then_some(mass_delta),
                blocks_cleavage,
                diagnostic_ions: self.diagnostic_ions.clone(),
                ..Modification::new(name.clone(), mass_delta, residue, location, fixed)
            })
            .collect())
    }
}

/// Split a comma-separated list at the commas outside parentheses. An empty
/// or blank list has no items.
pub fn split_list(list: &str) -> Result<Vec<&str>, ModSpecError> {
    if list.trim().is_empty() {
        return Ok(Vec::new());
    }
    let mut items = Vec::new();
    let (mut depth, mut start) = (0i32, 0usize);
    for (i, c) in list.char_indices() {
        match c {
            '(' => depth += 1,
            ')' => depth -= 1,
            ',' if depth == 0 => {
                items.push(&list[start..i]);
                start = i + 1;
            }
            _ => {}
        }
    }
    items.push(&list[start..]);
    if items.iter().any(|i| i.trim().is_empty()) {
        return Err(ModSpecError::EmptyItem {
            list: list.to_string(),
        });
    }
    Ok(items.into_iter().map(str::trim).collect())
}

/// Parse a `Name (sites), ...` list against the bundled UNIMOD table.
pub fn parse_mod_list(list: &str, fixed: bool) -> Result<Vec<Modification>, ModSpecError> {
    parse_mod_list_with(list, fixed, crate::unimod::bundled())
}

/// [`parse_mod_list`] against a given UNIMOD table.
pub fn parse_mod_list_with(
    list: &str,
    fixed: bool,
    db: &UnimodDb,
) -> Result<Vec<Modification>, ModSpecError> {
    let mut mods = Vec::new();
    for item in split_list(list)? {
        mods.extend(ModDecl::parse(item)?.resolve(fixed, db)?);
    }
    Ok(mods)
}

/// `Name (site)` of one expanded modification, as a list item would name it.
pub fn describe(m: &Modification) -> String {
    let terminus = match m.location {
        ModLocation::Anywhere => "",
        ModLocation::NTerm => "N-term",
        ModLocation::CTerm => "C-term",
        ModLocation::ProtNTerm => "Protein N-term",
        ModLocation::ProtCTerm => "Protein C-term",
    };
    let residue = match m.residue {
        ResidueSpec::Specific(r) => (r as char).to_string(),
        ResidueSpec::Wildcard => String::new(),
    };
    let sep = if terminus.is_empty() || residue.is_empty() {
        ""
    } else {
        " "
    };
    format!("{} ({terminus}{sep}{residue})", m.name)
}

/// The sites of `spec`: `(residue, location)` per site, in written order.
fn parse_sites(item: &str, spec: &str) -> Result<Vec<(ResidueSpec, ModLocation)>, ModSpecError> {
    const TERMINI: &[(&str, ModLocation)] = &[
        ("protein n-term", ModLocation::ProtNTerm),
        ("protein c-term", ModLocation::ProtCTerm),
        ("any n-term", ModLocation::NTerm),
        ("any c-term", ModLocation::CTerm),
        ("n-term", ModLocation::NTerm),
        ("c-term", ModLocation::CTerm),
    ];
    let bad = || ModSpecError::BadSites {
        item: item.to_string(),
        sites: spec.to_string(),
    };
    let spec_t = spec.trim();
    let lower = spec_t.to_ascii_lowercase();
    let (location, residues) = match TERMINI.iter().find(|(word, _)| {
        lower.starts_with(word)
            && lower[word.len()..]
                .chars()
                .next()
                .is_none_or(char::is_whitespace)
    }) {
        Some((word, loc)) => (*loc, spec_t[word.len()..].trim()),
        None => (ModLocation::Anywhere, spec_t),
    };
    if !residues.bytes().all(|b| b.is_ascii_uppercase()) {
        return Err(bad());
    }
    if residues.is_empty() {
        return match location {
            ModLocation::Anywhere => Err(bad()),
            loc => Ok(vec![(ResidueSpec::Wildcard, loc)]),
        };
    }
    let mut sites: Vec<(ResidueSpec, ModLocation)> = Vec::new();
    for r in residues.bytes() {
        let site = (ResidueSpec::Specific(r), location);
        if !sites.contains(&site) {
            sites.push(site);
        }
    }
    Ok(sites)
}

/// The UNIMOD record `token` names: an accession, a name, then an alias.
fn lookup<'a>(db: &'a UnimodDb, token: &str, item: &str) -> Result<&'a UnimodMod, ModSpecError> {
    match db.resolve(token) {
        Ok(Some(m)) => Ok(m),
        Ok(None) => db
            .resolve_alias(token)
            .ok_or_else(|| ModSpecError::UnknownName {
                item: item.to_string(),
                name: token.trim().to_string(),
                suggestions: suggest(db, token),
                table_empty: db.is_empty(),
            }),
        Err(UnimodResolveError::UnknownUnimodAccession { token }) => {
            Err(ModSpecError::UnknownAccession {
                item: item.to_string(),
                token,
            })
        }
        Err(UnimodResolveError::AmbiguousModName { token, accessions }) => {
            Err(ModSpecError::AmbiguousName {
                item: item.to_string(),
                token,
                accessions,
            })
        }
    }
}

/// Up to three names close to `token`: by edit distance (case-insensitive),
/// at most a third of the token's length (at least 2), or names that start
/// with the token.
pub fn suggest(db: &UnimodDb, token: &str) -> Vec<String> {
    let t = token.trim().to_ascii_lowercase();
    let limit = (t.chars().count() / 3).max(2);
    let mut scored: Vec<(usize, &str)> = db
        .names()
        .filter_map(|name| {
            let n = name.to_ascii_lowercase();
            let d = edit_distance(&t, &n);
            if d <= limit {
                Some((d, name))
            } else if t.len() >= 3 && n.starts_with(&t) {
                Some((limit + 1 + n.len() - t.len(), name))
            } else {
                None
            }
        })
        .collect();
    scored.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.len().cmp(&b.1.len())).then(a.1.cmp(b.1)));
    let mut out: Vec<String> = Vec::new();
    for (_, name) in scored {
        if !out.iter().any(|o| o == name) {
            out.push(name.to_string());
        }
        if out.len() == 3 {
            break;
        }
    }
    out
}

/// Levenshtein distance over chars.
fn edit_distance(a: &str, b: &str) -> usize {
    let b: Vec<char> = b.chars().collect();
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    let mut cur = vec![0; b.len() + 1];
    for (i, ca) in a.chars().enumerate() {
        cur[0] = i + 1;
        for (j, cb) in b.iter().enumerate() {
            let sub = prev[j] + usize::from(ca != *cb);
            cur[j + 1] = sub.min(prev[j + 1] + 1).min(cur[j] + 1);
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    prev[b.len()]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn db() -> UnimodDb {
        UnimodDb::parse_obo(
            r#"
[Term]
id: UNIMOD:4
name: Carbamidomethyl
xref: delta_mono_mass "57.021464"

[Term]
id: UNIMOD:35
name: Oxidation
xref: delta_mono_mass "15.994915"

[Term]
id: UNIMOD:21
name: Phospho
xref: delta_mono_mass "79.966331"

[Term]
id: UNIMOD:1
name: Acetyl
xref: delta_mono_mass "42.010565"

[Term]
id: UNIMOD:28
name: Gln->pyro-Glu
xref: delta_mono_mass "-17.026549"

[Term]
id: UNIMOD:121
name: GG
xref: delta_mono_mass "114.042927"

[Term]
id: UNIMOD:188
name: Label:13C(6)
xref: delta_mono_mass "6.020129"
"#,
        )
        .with_aliases("# comment\nGlyGly\tUNIMOD:121\nMissing\tUNIMOD:999\n")
    }

    fn sites(mods: &[Modification]) -> Vec<(ResidueSpec, ModLocation)> {
        mods.iter().map(|m| (m.residue, m.location)).collect()
    }

    #[test]
    fn names_resolve_mass_name_and_accession() {
        let mods = parse_mod_list_with("phospho (STY)", false, &db()).unwrap();
        assert_eq!(mods.len(), 3);
        for m in &mods {
            assert_eq!(m.name, "Phospho");
            assert_eq!(m.mass_delta, 79.966331);
            assert_eq!(m.accession.as_deref(), Some("UNIMOD:21"));
            assert!(!m.fixed);
            assert!(m.neutral_losses.is_empty() && m.loss_class == 0);
            assert!(m.annotated_site_mass.is_none() && !m.blocks_cleavage);
        }
        use ModLocation::Anywhere;
        assert_eq!(
            sites(&mods),
            vec![
                (ResidueSpec::Specific(b'S'), Anywhere),
                (ResidueSpec::Specific(b'T'), Anywhere),
                (ResidueSpec::Specific(b'Y'), Anywhere),
            ]
        );
    }

    #[test]
    fn accessions_and_aliases_resolve() {
        let m = parse_mod_list_with("UNIMOD:21 (S)", false, &db()).unwrap();
        assert_eq!(m[0].name, "Phospho");
        let m = parse_mod_list_with("unimod:4 (C)", true, &db()).unwrap();
        assert_eq!((m[0].name.as_str(), m[0].fixed), ("Carbamidomethyl", true));
        // An alias resolves to its record, under the record's UNIMOD name.
        let m = parse_mod_list_with("GlyGly (K)", false, &db()).unwrap();
        assert_eq!(m[0].name, "GG");
        assert_eq!(m[0].accession.as_deref(), Some("UNIMOD:121"));
        assert!(matches!(
            parse_mod_list_with("UNIMOD:999 (K)", false, &db()),
            Err(ModSpecError::UnknownAccession { .. })
        ));
    }

    #[test]
    fn site_specs() {
        let db = db();
        let one = |s: &str| sites(&parse_mod_list_with(s, false, &db).unwrap());
        use ModLocation::*;
        use ResidueSpec::*;
        assert_eq!(one("Acetyl (Protein N-term)"), vec![(Wildcard, ProtNTerm)]);
        assert_eq!(one("Acetyl (protein c-term)"), vec![(Wildcard, ProtCTerm)]);
        assert_eq!(one("Acetyl (N-term)"), vec![(Wildcard, NTerm)]);
        assert_eq!(one("Acetyl (Any N-term)"), vec![(Wildcard, NTerm)]);
        assert_eq!(one("Acetyl (C-term)"), vec![(Wildcard, CTerm)]);
        assert_eq!(one("Gln->pyro-Glu (N-term Q)"), vec![(Specific(b'Q'), NTerm)]);
        assert_eq!(
            one("Acetyl (Protein N-term MS)"),
            vec![(Specific(b'M'), ProtNTerm), (Specific(b'S'), ProtNTerm)]
        );
        // Repeated residues count once; names with parentheses keep them.
        assert_eq!(one("Oxidation (MM)"), vec![(Specific(b'M'), Anywhere)]);
        let m = parse_mod_list_with("Label:13C(6) (K)", false, &db).unwrap();
        assert_eq!(m[0].name, "Label:13C(6)");
    }

    #[test]
    fn lists_split_outside_parentheses() {
        let db = db();
        let mods = parse_mod_list_with(
            "Oxidation (M), Phospho (STY),Acetyl (Protein N-term), Gln->pyro-Glu (N-term Q)",
            false,
            &db,
        )
        .unwrap();
        assert_eq!(mods.len(), 6);
        assert!(parse_mod_list_with("", false, &db).unwrap().is_empty());
        assert!(parse_mod_list_with("  ", false, &db).unwrap().is_empty());
        assert!(matches!(
            parse_mod_list_with("Oxidation (M),, Phospho (S)", false, &db),
            Err(ModSpecError::EmptyItem { .. })
        ));
        assert!(matches!(
            parse_mod_list_with("Oxidation (M),", false, &db),
            Err(ModSpecError::EmptyItem { .. })
        ));
    }

    #[test]
    fn malformed_items_are_errors() {
        let db = db();
        let err = |s: &str| parse_mod_list_with(s, false, &db).unwrap_err();
        assert!(matches!(err("Oxidation"), ModSpecError::MissingSites { .. }));
        assert!(matches!(err("(M)"), ModSpecError::MissingName { .. }));
        assert!(matches!(err("Oxidation M)"), ModSpecError::Unbalanced { .. }));
        for bad in [
            "Oxidation (m)",
            "Oxidation ()",
            "Oxidation (Q N-term)",
            "Oxidation (N-termQ)",
            "Oxidation (anywhere)",
        ] {
            assert!(
                matches!(err(bad), ModSpecError::BadSites { .. }),
                "{bad}: {:?}",
                err(bad)
            );
        }
    }

    #[test]
    fn unknown_names_suggest_close_ones() {
        let db = db();
        match parse_mod_list_with("Phosho (STY)", false, &db).unwrap_err() {
            ModSpecError::UnknownName {
                name, suggestions, ..
            } => {
                assert_eq!(name, "Phosho");
                assert_eq!(suggestions, vec!["Phospho".to_string()]);
            }
            other => panic!("{other:?}"),
        }
        let msg = parse_mod_list_with("Oxidaton (M)", false, &db)
            .unwrap_err()
            .to_string();
        assert!(msg.contains("did you mean Oxidation?"), "{msg}");
        // A prefix suggests the names it starts.
        assert_eq!(suggest(&db, "Carbamido"), vec!["Carbamidomethyl".to_string()]);
        assert_eq!(suggest(&db, "GlyGl"), vec!["GlyGly".to_string()]);
        let msg = parse_mod_list_with("Frobnicate (K)", false, &db)
            .unwrap_err()
            .to_string();
        assert!(msg.contains("not a UNIMOD name"), "{msg}");
        let msg = parse_mod_list_with("Oxidation (M)", false, &UnimodDb::default())
            .unwrap_err()
            .to_string();
        assert!(msg.contains("could not be read"), "{msg}");
    }

    #[test]
    fn config_attributes_apply_to_every_site() {
        let decl = ModDecl {
            name: Some("Phospho".into()),
            sites: "ST".into(),
            neutral_loss: vec![97.976896],
            restrict_to: Some("annotated".into()),
            diagnostic_ions: vec![216.0426],
            ..ModDecl::default()
        };
        let mods = decl.resolve(false, &db()).unwrap();
        assert_eq!(mods.len(), 2);
        for m in &mods {
            assert_eq!(m.neutral_losses, vec![97.976896]);
            assert_eq!(m.loss_class, 255, "losses without a class are generic");
            assert_eq!(m.annotated_site_mass, Some(79.966331));
            assert_eq!(m.diagnostic_ions, vec![216.0426]);
        }
        let gg = ModDecl {
            name: Some("GlyGly".into()),
            sites: "K".into(),
            cleavage: Some("Blocked".into()),
            loss_class: Some("generic".into()),
            ..ModDecl::default()
        };
        let m = &gg.resolve(false, &db()).unwrap()[0];
        assert!(m.blocks_cleavage);
        assert_eq!(m.loss_class, 255);
        let glyco = ModDecl {
            loss_class: Some("glyco".into()),
            neutral_loss: vec![162.0528],
            ..gg.clone()
        };
        assert_eq!(glyco.resolve(false, &db()).unwrap()[0].loss_class, 1);
    }

    #[test]
    fn config_attribute_errors() {
        let base = ModDecl {
            name: Some("Phospho".into()),
            sites: "S".into(),
            ..ModDecl::default()
        };
        let err = |d: ModDecl, fixed: bool| d.resolve(fixed, &db()).unwrap_err();
        assert!(matches!(
            err(
                ModDecl {
                    neutral_loss: vec![0.0],
                    ..base.clone()
                },
                false
            ),
            ModSpecError::BadNeutralLoss { .. }
        ));
        assert!(matches!(
            err(
                ModDecl {
                    neutral_loss: vec![5000.0],
                    ..base.clone()
                },
                false
            ),
            ModSpecError::BadNeutralLoss { .. }
        ));
        assert!(matches!(
            err(
                ModDecl {
                    diagnostic_ions: vec![-1.0],
                    ..base.clone()
                },
                false
            ),
            ModSpecError::BadDiagnosticIon { .. }
        ));
        assert!(matches!(
            err(
                ModDecl {
                    loss_class: Some("bogus".into()),
                    ..base.clone()
                },
                false
            ),
            ModSpecError::BadLossClass { .. }
        ));
        assert!(matches!(
            err(
                ModDecl {
                    restrict_to: Some("some".into()),
                    ..base.clone()
                },
                false
            ),
            ModSpecError::BadRestrictTo { .. }
        ));
        assert!(matches!(
            err(
                ModDecl {
                    restrict_to: Some("annotated".into()),
                    ..base.clone()
                },
                true
            ),
            ModSpecError::AnnotatedFixed { .. }
        ));
        assert!(matches!(
            err(
                ModDecl {
                    cleavage: Some("maybe".into()),
                    ..base.clone()
                },
                false
            ),
            ModSpecError::BadCleavage { .. }
        ));
        assert!(matches!(
            err(
                ModDecl {
                    name: None,
                    ..base.clone()
                },
                false
            ),
            ModSpecError::NoNameOrMass { .. }
        ));
    }

    #[test]
    fn custom_masses_keep_their_name_and_have_no_accession() {
        let decl = ModDecl {
            name: Some("MyTag".into()),
            mass: Some(12.3456),
            sites: "N-term".into(),
            ..ModDecl::default()
        };
        let m = &decl.resolve(true, &UnimodDb::default()).unwrap()[0];
        assert_eq!((m.name.as_str(), m.mass_delta), ("MyTag", 12.3456));
        assert_eq!((m.residue, m.location), (ResidueSpec::Wildcard, ModLocation::NTerm));
        assert!(m.fixed && m.accession.is_none());
        // A UNIMOD name with a mass is a label, not a lookup.
        let decl = ModDecl {
            name: Some("Phospho".into()),
            mass: Some(80.0),
            sites: "S".into(),
            ..ModDecl::default()
        };
        let m = &decl.resolve(false, &db()).unwrap()[0];
        assert_eq!((m.mass_delta, m.accession.clone()), (80.0, None));
        let unnamed = ModDecl {
            mass: Some(-1.5),
            sites: "K".into(),
            ..ModDecl::default()
        };
        assert_eq!(unnamed.resolve(false, &db()).unwrap()[0].name, "-1.5");
        assert!(matches!(
            ModDecl {
                mass: Some(f64::NAN),
                ..unnamed
            }
            .resolve(false, &db()),
            Err(ModSpecError::BadMass { .. })
        ));
    }

    #[test]
    fn describe_names_one_site() {
        let mods = parse_mod_list_with(
            "Phospho (S), Acetyl (Protein N-term), Gln->pyro-Glu (N-term Q)",
            false,
            &db(),
        )
        .unwrap();
        let d: Vec<String> = mods.iter().map(describe).collect();
        assert_eq!(
            d,
            vec![
                "Phospho (S)",
                "Acetyl (Protein N-term)",
                "Gln->pyro-Glu (N-term Q)"
            ]
        );
    }

    #[test]
    fn bundled_table_knows_the_common_names() {
        for (list, name, mass) in [
            ("Carbamidomethyl (C)", "Carbamidomethyl", 57.021464),
            ("Oxidation (M)", "Oxidation", 15.994915),
            ("Phospho (S)", "Phospho", 79.966331),
            ("Acetyl (Protein N-term)", "Acetyl", 42.010565),
            ("Gln->pyro-Glu (N-term Q)", "Gln->pyro-Glu", -17.026549),
            ("Deamidated (N)", "Deamidated", 0.984016),
            ("TMT6plex (K)", "TMT6plex", 229.162932),
            ("TMT10plex (K)", "TMT6plex", 229.162932),
            ("TMTpro (K)", "TMTpro", 304.207146),
            ("iTRAQ4plex (K)", "iTRAQ4plex", 144.102063),
            ("iTRAQ8plex (K)", "iTRAQ8plex", 304.20536),
            ("GlyGly (K)", "GG", 114.042927),
            ("GG (K)", "GG", 114.042927),
            ("Label:13C(6)15N(2) (K)", "Label:13C(6)15N(2)", 8.014199),
        ] {
            let m = &parse_mod_list(list, false).unwrap()[0];
            assert_eq!(m.name, name, "{list}");
            assert!((m.mass_delta - mass).abs() < 1e-6, "{list}: {}", m.mass_delta);
        }
    }
}
