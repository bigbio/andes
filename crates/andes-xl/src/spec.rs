//! The crosslinker of a search: `--xl` on the command line and `search.xl`
//! in the config file.
//!
//! The spec is written like a modification list item, `Name (sites)`:
//!
//! - **Name**: a reagent of the built-in table ([`crate::reagents`]) by name
//!   or alias (`DSSO`, `BS3`, `DSBU`, `PhoX`, any case), or the UNIMOD name or
//!   accession of its intact cross-link record (`Xlink:DSS[138]`,
//!   `UNIMOD:1898`). Masses come from UNIMOD when it has the record and from
//!   the cited composition otherwise.
//! - **Sites** (optional; the reagent's usual reactivity when omitted):
//!   comma-separated residues and protein termini, each in the modification
//!   site grammar (`K, S, T, Y, Protein N-term`; `KSTY` also works). A
//!   heterobifunctional reagent gives its two ends apart with `/`:
//!   `EDC (K, Protein N-term / D, E, Protein C-term)`. `N-term` is the
//!   N-terminus of every peptide, for synthetic peptide libraries or samples
//!   linked after digestion (a protein linked before digestion has only its
//!   protein N-terminus free).
//!
//! A config `search.xl` can also be a map ([`XlDecl`]) that defines a linker
//! the table lacks: intact `mass`, `sites`, and for an MS-cleavable linker its
//! `stubs` and reporter `doublets`.

use model::modification::{ModLocation, ResidueSpec};
use model::unimod::UnimodDb;

use crate::crosslinker::{Cleavage, Crosslinker, LinkEnd, Monolink, Stub, MAX_STUBS};
use crate::reagents::{self, Part, Reagent};

/// One declared crosslinker before it is resolved.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct XlDecl {
    /// Reagent name, alias or UNIMOD record; with `mass`, a free label.
    pub name: Option<String>,
    /// Custom intact mass in Da; the name is then not looked up.
    pub mass: Option<f64>,
    /// Site spec; the reagent's default when `None`.
    pub sites: Option<String>,
    /// Custom MS-cleavage stubs, `(name, mass)`.
    pub stubs: Vec<(String, f64)>,
    /// Custom reporter doublets as `(light, heavy)` stub names.
    pub doublets: Vec<(String, String)>,
    /// Custom monolinks, `(name, mass)`.
    pub monolinks: Vec<(String, f64)>,
}

/// Why a crosslinker spec was rejected. Each message names the spec.
#[derive(thiserror::Error, Debug, PartialEq)]
pub enum XlSpecError {
    #[error("`{spec}`: unknown crosslinker `{name}`{}", hint(.suggestions))]
    Unknown {
        spec: String,
        name: String,
        suggestions: Vec<String>,
    },
    #[error("`{spec}`: `{name}` is a {role} of {reagent}; give `{reagent}` as the crosslinker{}",
            monolink_hint(.role))]
    NotALinker {
        spec: String,
        name: String,
        role: &'static str,
        reagent: &'static str,
    },
    #[error("`{spec}`: the UNIMOD record {name} is not available (the bundled UNIMOD table resources/unimod.obo could not be read)")]
    MissingUnimod { spec: String, name: String },
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
    #[error("`{spec}`: a custom linker (with `mass`) needs `sites`")]
    CustomWithoutSites { spec: String },
    #[error("`{spec}`: mass {mass} is not a finite number")]
    BadMass { spec: String, mass: f64 },
    #[error("`{spec}`: {why}")]
    BadCleavage { spec: String, why: String },
}

fn monolink_hint(role: &str) -> &'static str {
    if role == "monolink" {
        " (and the monolink to --variable-mods)"
    } else {
        ""
    }
}

fn hint(suggestions: &[String]) -> String {
    if suggestions.is_empty() {
        format!(
            " (known: {})",
            reagents::names().collect::<Vec<_>>().join(", ")
        )
    } else {
        format!("; did you mean {}?", suggestions.join(", "))
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
            name: Some(name.to_string()),
            sites: sites.map(str::to_string),
            ..XlDecl::default()
        })
    }

    /// How the spec reads in messages.
    pub fn label(&self) -> String {
        let name = match (&self.name, self.mass) {
            (Some(n), _) => n.clone(),
            (None, Some(m)) => format!("{m:+}"),
            (None, None) => String::new(),
        };
        match &self.sites {
            Some(s) => format!("{name} ({s})"),
            None => name,
        }
    }

    /// The crosslinker this declaration stands for.
    pub fn resolve(&self, db: &UnimodDb) -> Result<Crosslinker, XlSpecError> {
        let spec = self.label();
        match self.mass {
            Some(mass) => self.resolve_custom(&spec, mass),
            None => {
                if !(self.stubs.is_empty() && self.doublets.is_empty() && self.monolinks.is_empty())
                {
                    return Err(XlSpecError::BadCleavage {
                        spec,
                        why: "`stubs`, `doublets` and `monolinks` define a custom linker and need \
                              its intact `mass`"
                            .into(),
                    });
                }
                let name = self
                    .name
                    .as_deref()
                    .ok_or(XlSpecError::MissingName { spec: spec.clone() })?;
                let reagent = lookup(db, name, &spec)?;
                let mut x = from_reagent(reagent, db, &spec)?;
                if let Some(sites) = &self.sites {
                    x.ends = parse_ends(&spec, sites)?;
                }
                Ok(x)
            }
        }
    }

    fn resolve_custom(&self, spec: &str, mass: f64) -> Result<Crosslinker, XlSpecError> {
        let finite = |m: f64| {
            if m.is_finite() {
                Ok(m)
            } else {
                Err(XlSpecError::BadMass {
                    spec: spec.to_string(),
                    mass: m,
                })
            }
        };
        finite(mass)?;
        let sites = self
            .sites
            .as_deref()
            .ok_or(XlSpecError::CustomWithoutSites {
                spec: spec.to_string(),
            })?;
        let bad = |why: String| XlSpecError::BadCleavage {
            spec: spec.to_string(),
            why,
        };
        let cleavage = if self.stubs.is_empty() {
            if !self.doublets.is_empty() {
                return Err(bad("`doublets` without `stubs`".into()));
            }
            None
        } else {
            let mut stubs = Vec::new();
            for (name, m) in &self.stubs {
                finite(*m)?;
                if stubs.iter().any(|s: &Stub| s.name == *name) {
                    return Err(bad(format!("stub `{name}` given twice")));
                }
                stubs.push(Stub {
                    name: name.clone(),
                    mass: *m,
                });
            }
            let index = |n: &str| {
                stubs
                    .iter()
                    .position(|s| s.name == n)
                    .ok_or_else(|| bad(format!("doublet stub `{n}` is not one of the stubs")))
            };
            let mut doublets = Vec::new();
            if self.doublets.is_empty() {
                if stubs.len() != 2 {
                    return Err(bad(
                        "give `doublets` (pairs of stub names) when there are not exactly two stubs"
                            .into(),
                    ));
                }
                doublets.push((0, 1));
            } else {
                for (a, b) in &self.doublets {
                    doublets.push((index(a)?, index(b)?));
                }
            }
            let c = Cleavage { stubs, doublets };
            check_cleavage(&c).map_err(bad)?;
            Some(order_doublets(c))
        };
        let mut monolinks = Vec::new();
        for (name, m) in &self.monolinks {
            monolinks.push(Monolink {
                name: name.clone(),
                accession: None,
                mass: finite(*m)?,
            });
        }
        Ok(Crosslinker {
            name: self.name.clone().unwrap_or_else(|| format!("{mass:+}")),
            accession: None,
            intact: mass,
            cleavage,
            ends: parse_ends(spec, sites)?,
            monolinks,
        })
    }
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
        return Err("an MS-cleavable linker needs at least one reporter doublet".into());
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

/// The reagent `name` stands for: the table by name or alias, then UNIMOD.
fn lookup(db: &UnimodDb, name: &str, spec: &str) -> Result<&'static Reagent, XlSpecError> {
    if let Some(r) = reagents::find(name) {
        return Ok(r);
    }
    // A UNIMOD accession or name: its record must be a reagent's intact link.
    if let Ok(Some(rec)) = db.resolve(name) {
        if let Some(r) = reagents::find(&rec.name) {
            return Ok(r);
        }
        if let Some((r, role)) = reagents::owner_of_record(&rec.name) {
            return Err(XlSpecError::NotALinker {
                spec: spec.to_string(),
                name: rec.name.clone(),
                role,
                reagent: r.name,
            });
        }
    }
    Err(XlSpecError::Unknown {
        spec: spec.to_string(),
        name: name.trim().to_string(),
        suggestions: suggest(name),
    })
}

/// Up to three reagent names within a small edit distance of `token`.
fn suggest(token: &str) -> Vec<String> {
    let t = token.trim().to_ascii_lowercase();
    let mut scored: Vec<(usize, &str)> = reagents::names()
        .map(|n| (edit_distance(&t, &n.to_ascii_lowercase()), n))
        .filter(|(d, _)| *d <= 2)
        .collect();
    scored.sort();
    scored
        .into_iter()
        .take(3)
        .map(|(_, n)| n.to_string())
        .collect()
}

fn edit_distance(a: &str, b: &str) -> usize {
    let b: Vec<char> = b.chars().collect();
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    let mut cur = vec![0; b.len() + 1];
    for (i, ca) in a.chars().enumerate() {
        cur[0] = i + 1;
        for (j, cb) in b.iter().enumerate() {
            cur[j + 1] = (prev[j] + usize::from(ca != *cb))
                .min(prev[j + 1] + 1)
                .min(cur[j] + 1);
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    prev[b.len()]
}

/// The mass of one reagent part: its UNIMOD record, or its composition.
fn part_mass(part: &Part, db: &UnimodDb, spec: &str) -> Result<(f64, Option<String>), XlSpecError> {
    match part.unimod {
        Some(name) => match db.resolve(name) {
            Ok(Some(rec)) => Ok((rec.mono_mass, Some(rec.accession.clone()))),
            _ => Err(XlSpecError::MissingUnimod {
                spec: spec.to_string(),
                name: name.to_string(),
            }),
        },
        None => Ok((
            reagents::composition_mass(part.composition)
                .expect("reagent table compositions parse (tested)"),
            None,
        )),
    }
}

fn from_reagent(r: &Reagent, db: &UnimodDb, spec: &str) -> Result<Crosslinker, XlSpecError> {
    let (intact, accession) = part_mass(&r.intact, db, spec)?;
    let cleavage = if r.stubs.is_empty() {
        None
    } else {
        let mut stubs = Vec::new();
        for s in r.stubs {
            stubs.push(Stub {
                name: s.name.to_string(),
                mass: part_mass(s, db, spec)?.0,
            });
        }
        Some(Cleavage {
            stubs,
            doublets: r.doublets.to_vec(),
        })
    };
    let mut monolinks = Vec::new();
    for m in r.monolinks {
        let (mass, accession) = part_mass(m, db, spec)?;
        monolinks.push(Monolink {
            name: m
                .unimod
                .map_or_else(|| format!("{}-{}", r.name, m.name), str::to_string),
            accession,
            mass,
        });
    }
    Ok(Crosslinker {
        name: r.name.to_string(),
        accession,
        intact,
        cleavage,
        ends: parse_ends(spec, r.sites).expect("reagent table sites parse (tested)"),
        monolinks,
    })
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
        for r in reagents::REAGENTS {
            let x = resolve(r.name).unwrap_or_else(|e| panic!("{}: {e}", r.name));
            assert_eq!(x.name, r.name);
            assert_eq!(x.cleavage.is_some(), !r.stubs.is_empty());
            assert_eq!(x.monolinks.len(), r.monolinks.len());
            assert!(!x.ends[0].residues.is_empty() || x.ends[0].protein_n_term);
            // Every alias resolves to the same linker.
            for a in r.aliases {
                assert_eq!(resolve(a).unwrap(), x, "{a}");
            }
            if let Some(u) = r.intact.unimod {
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
            resolve("DSSO (K"),
            Err(XlSpecError::Unbalanced { .. }) | Err(XlSpecError::Unknown { .. })
        ));
        assert!(matches!(
            resolve(" (K)"),
            Err(XlSpecError::MissingName { .. })
        ));
        let e = resolve("DSSO (k)").unwrap_err().to_string();
        assert!(e.contains("bad sites"), "{e}");
    }

    #[test]
    fn custom_cleavable_linker() {
        let d = XlDecl {
            name: Some("MyDSSO".into()),
            mass: Some(158.003_765),
            sites: Some("K, Protein N-term".into()),
            stubs: vec![
                ("thiol".into(), 85.982_635),
                ("alkene".into(), 54.010_565),
                ("sulfenic".into(), 103.993_2),
            ],
            doublets: vec![("thiol".into(), "alkene".into())],
            monolinks: vec![("water".into(), 176.014_33)],
        };
        let x = d.resolve(db()).unwrap();
        let c = x.cleavage.as_ref().unwrap();
        // The doublet is ordered light (alkene) then heavy (thiol).
        assert_eq!(c.doublets, vec![(1, 0)]);
        assert!((c.doublet_delta(0) - 31.972_07).abs() < 1e-4);
        assert_eq!(x.monolinks[0].mass, 176.014_33);
        // Two stubs need no `doublets`.
        let two = XlDecl {
            stubs: vec![("a".into(), 10.0), ("b".into(), 20.0)],
            doublets: vec![],
            ..d.clone()
        };
        assert_eq!(
            two.resolve(db()).unwrap().cleavage.unwrap().doublets,
            vec![(0, 1)]
        );
        // Three stubs do.
        let three = XlDecl {
            doublets: vec![],
            ..d.clone()
        };
        assert!(three
            .resolve(db())
            .unwrap_err()
            .to_string()
            .contains("give `doublets`"));
        let no_sites = XlDecl {
            sites: None,
            ..d.clone()
        };
        assert!(matches!(
            no_sites.resolve(db()),
            Err(XlSpecError::CustomWithoutSites { .. })
        ));
        let unknown_stub = XlDecl {
            doublets: vec![("alkene".into(), "thiolx".into())],
            ..d
        };
        assert!(unknown_stub.resolve(db()).is_err());
    }

    #[test]
    fn custom_non_cleavable_linker() {
        let d = XlDecl {
            name: Some("ZL".into()),
            mass: Some(-18.010_565),
            sites: Some("K / D, E".into()),
            ..XlDecl::default()
        };
        let x = d.resolve(db()).unwrap();
        assert!(x.cleavage.is_none() && !x.is_homobifunctional());
    }
}
