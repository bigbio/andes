//! The crosslinking reagents andes knows by name, loaded from the bundled
//! `resources/crosslinkers.yaml`: the crosslinker table of `--xl` (resolved
//! in the `andes-xl` crate) and the monolinks that `--variable-mods` accepts
//! by name ([`monolink`]).
//!
//! The file holds no code, only chemistry: per reagent its name, aliases,
//! XLMOD accessions, default sites, and its intact, stub and monolink masses.
//! Each mass is a UNIMOD record (then the UNIMOD mass is the one used) with its
//! elemental composition, or a composition alone for a reagent UNIMOD lacks;
//! a custom reagent of the config file may give a plain `mass`. The schema is
//! [`Reagent`], shared by the bundled file and the config file.

use std::sync::OnceLock;

use serde::Deserialize;

use crate::unimod::{UnimodDb, UnimodMod};

/// Most MS-cleavage stubs a reagent can have.
pub const MAX_STUBS: usize = 4;

/// One mass of a reagent: the intact cross-link, a stub or a monolink.
#[derive(Debug, Clone, Default, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Part {
    /// Column or display name (a stub's short name, a monolink's quencher).
    #[serde(default)]
    pub name: String,
    /// The UNIMOD record of this mass, if UNIMOD has it.
    #[serde(default)]
    pub unimod: Option<String>,
    /// Elemental composition in UNIMOD notation (`H(6) C(6) O(3) S`;
    /// negative counts allowed).
    #[serde(default)]
    pub composition: Option<String>,
    /// A mass in Da, for a custom reagent without a composition.
    #[serde(default)]
    pub mass: Option<f64>,
}

impl Part {
    /// The mass and UNIMOD accession of this part: the UNIMOD record's mass
    /// when it names one, else the composition's, else `mass`.
    pub fn resolve(&self, db: &UnimodDb) -> Result<(f64, Option<String>), String> {
        if let Some(name) = &self.unimod {
            return match db.resolve(name) {
                Ok(Some(rec)) => Ok((rec.mono_mass, Some(rec.accession.clone()))),
                _ => Err(format!(
                    "UNIMOD record {name} not found (the bundled UNIMOD table \
                     resources/unimod.obo could not be read, or lacks it)"
                )),
            };
        }
        if let Some(c) = &self.composition {
            return composition_mass(c)
                .map(|m| (m, None))
                .ok_or_else(|| format!("bad composition `{c}`"));
        }
        match self.mass {
            Some(m) if m.is_finite() => Ok((m, None)),
            Some(m) => Err(format!("mass {m} is not a finite number")),
            None => Err(format!(
                "`{}` needs a `unimod` record, a `composition` or a `mass`",
                if self.name.is_empty() {
                    "intact"
                } else {
                    &self.name
                }
            )),
        }
    }
}

/// A crosslinking reagent: one entry of `resources/crosslinkers.yaml`, or a
/// custom reagent of the config file.
#[derive(Debug, Clone, Default, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Reagent {
    /// Canonical name.
    pub name: String,
    /// Other names accepted for it.
    #[serde(default)]
    pub aliases: Vec<String>,
    /// XLMOD accessions of the reagent and its aliases.
    #[serde(default)]
    pub xlmod: Vec<String>,
    /// Default sites, as in a spec (`K, Protein N-term`; `end / end` for a
    /// heterobifunctional reagent).
    #[serde(default)]
    pub sites: String,
    /// The intact cross-link.
    pub intact: Part,
    /// MS-cleavage stubs; empty for a non-cleavable reagent.
    #[serde(default)]
    pub stubs: Vec<Part>,
    /// Reporter doublets as `[light, heavy]` stub names.
    #[serde(default)]
    pub doublets: Vec<[String; 2]>,
    /// Quenched monolinks, named by quencher.
    #[serde(default)]
    pub monolinks: Vec<Part>,
    /// Where the chemistry and masses come from.
    #[serde(default)]
    pub source: Option<String>,
}

impl Reagent {
    /// Whether `name` is this reagent's name or an alias (any case).
    pub fn is_named(&self, name: &str) -> bool {
        let n = name.trim();
        self.name.eq_ignore_ascii_case(n) || self.aliases.iter().any(|a| a.eq_ignore_ascii_case(n))
    }
}

/// The reagent table.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ReagentTable {
    pub reagents: Vec<Reagent>,
}

impl ReagentTable {
    /// Parse a YAML list of [`Reagent`] entries.
    pub fn parse_yaml(text: &str) -> Result<Self, String> {
        let reagents: Vec<Reagent> = serde_yaml::from_str(text).map_err(|e| e.to_string())?;
        Ok(ReagentTable { reagents })
    }

    /// The reagent named `name`: a name or alias (any case), or the UNIMOD
    /// name of its intact record.
    pub fn find(&self, name: &str) -> Option<&Reagent> {
        let n = name.trim();
        self.reagents.iter().find(|r| {
            r.is_named(n)
                || r.intact
                    .unimod
                    .as_deref()
                    .is_some_and(|u| u.eq_ignore_ascii_case(n))
        })
    }

    /// The reagent with XLMOD accession `accession` (`XLMOD:02001`, any case).
    pub fn find_xlmod(&self, accession: &str) -> Option<&Reagent> {
        let a = accession.trim();
        self.reagents
            .iter()
            .find(|r| r.xlmod.iter().any(|x| x.eq_ignore_ascii_case(a)))
    }

    /// The reagent with a monolink or stub record named `unimod_name`, and
    /// the role of that record (`monolink` or `cleavage stub`).
    pub fn owner_of_record(&self, unimod_name: &str) -> Option<(&Reagent, &'static str)> {
        let is = |p: &Part| {
            p.unimod
                .as_deref()
                .is_some_and(|u| u.eq_ignore_ascii_case(unimod_name))
        };
        self.reagents.iter().find_map(|r| {
            if r.monolinks.iter().any(is) {
                Some((r, "monolink"))
            } else if r.stubs.iter().any(is) {
                Some((r, "cleavage stub"))
            } else {
                None
            }
        })
    }

    /// Every name [`Self::find`] accepts other than UNIMOD records.
    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.reagents.iter().flat_map(|r| {
            std::iter::once(r.name.as_str()).chain(r.aliases.iter().map(String::as_str))
        })
    }
}

/// The bundled reagent table (`resources/crosslinkers.yaml`), parsed once.
/// An unreadable file gives an empty table, so every lookup fails with the
/// caller's "unknown crosslinker" error; [`bundled_error`] says why.
pub fn bundled() -> &'static ReagentTable {
    load().as_ref().unwrap_or_else(|_| {
        static EMPTY: OnceLock<ReagentTable> = OnceLock::new();
        EMPTY.get_or_init(ReagentTable::default)
    })
}

/// Why the bundled reagent table could not be loaded, if it could not.
pub fn bundled_error() -> Option<&'static str> {
    load().as_ref().err().map(String::as_str)
}

fn load() -> &'static Result<ReagentTable, String> {
    static TABLE: OnceLock<Result<ReagentTable, String>> = OnceLock::new();
    TABLE.get_or_init(|| {
        let path = crate::unimod::bundled_path("crosslinkers.yaml");
        let text = std::fs::read_to_string(&path)
            .map_err(|e| format!("cannot read {}: {e}", path.display()))?;
        ReagentTable::parse_yaml(&text).map_err(|e| format!("{}: {e}", path.display()))
    })
}

/// The monolink `token` names, written `<reagent>:<quencher>` (`DSBSO:water`,
/// `PhoX:tris`, `DSSO:ammonia`; any case, reagent aliases accepted), as a
/// modification record: the UNIMOD record when UNIMOD has it, otherwise the
/// name `<Reagent>:<quencher>` with the mass of its composition and an empty
/// accession. `None` when `token` is not of that form or names no monolink of
/// the bundled table.
pub fn monolink(db: &UnimodDb, token: &str) -> Option<UnimodMod> {
    monolink_in(bundled(), db, token)
}

fn monolink_in(table: &ReagentTable, db: &UnimodDb, token: &str) -> Option<UnimodMod> {
    let (reagent, quencher) = token.trim().split_once(':')?;
    let r = table.find(reagent)?;
    let p = r
        .monolinks
        .iter()
        .find(|p| p.name.eq_ignore_ascii_case(quencher.trim()))?;
    match &p.unimod {
        Some(name) => db.resolve(name).ok().flatten().cloned(),
        None => Some(UnimodMod {
            accession: String::new(),
            name: format!("{}:{}", r.name, p.name),
            mono_mass: p.resolve(db).ok()?.0,
            sites: Vec::new(),
        }),
    }
}

/// Monoisotopic masses of the elements reagent compositions use (AME2020).
const ELEMENTS: &[(&str, f64)] = &[
    ("H", 1.007_825_032),
    ("C", 12.0),
    ("N", 14.003_074_004),
    ("O", 15.994_914_620),
    ("S", 31.972_071_174),
    ("P", 30.973_761_998),
];

/// Monoisotopic mass of a UNIMOD-style composition (`H(6) C(6) O(3) S`),
/// or `None` for an unknown element or a malformed count. The empty
/// composition weighs 0.
pub fn composition_mass(composition: &str) -> Option<f64> {
    let mut total = 0.0;
    for tok in composition.split_whitespace() {
        let (el, n) = match tok.find('(') {
            Some(i) => {
                let n: i32 = tok[i + 1..].strip_suffix(')')?.parse().ok()?;
                (&tok[..i], n)
            }
            None => (tok, 1),
        };
        let m = ELEMENTS.iter().find(|(e, _)| *e == el)?.1;
        total += m * n as f64;
    }
    Some(total)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn db() -> &'static UnimodDb {
        let db = crate::unimod::bundled();
        assert!(!db.is_empty(), "bundled UNIMOD table not found");
        db
    }

    fn table() -> &'static ReagentTable {
        assert_eq!(bundled_error(), None);
        let t = bundled();
        assert!(!t.reagents.is_empty());
        t
    }

    fn parts(r: &Reagent) -> impl Iterator<Item = &Part> {
        std::iter::once(&r.intact)
            .chain(r.stubs.iter())
            .chain(r.monolinks.iter())
    }

    fn mass_of(p: &Part) -> f64 {
        composition_mass(p.composition.as_deref().expect("composition")).unwrap()
    }

    #[test]
    fn the_bundled_file_parses_and_every_mass_has_a_composition() {
        for r in &table().reagents {
            assert!(!r.sites.is_empty(), "{}: no sites", r.name);
            assert!(r.source.is_some(), "{}: no source", r.name);
            for p in parts(r) {
                let c = p.composition.as_deref();
                assert!(
                    c.is_some_and(|c| composition_mass(c).is_some()),
                    "{} {}: composition {c:?}",
                    r.name,
                    p.name
                );
                assert!(
                    p.mass.is_none(),
                    "{} {}: bundled masses are compositions",
                    r.name,
                    p.name
                );
            }
        }
    }

    #[test]
    fn every_reagent_has_unique_xlmod_accessions_and_names() {
        let mut acc = std::collections::HashSet::new();
        let mut names = std::collections::HashSet::new();
        for r in &table().reagents {
            assert!(!r.xlmod.is_empty(), "{}: no XLMOD accession", r.name);
            for a in &r.xlmod {
                assert!(
                    a.strip_prefix("XLMOD:")
                        .is_some_and(|n| n.len() == 5 && n.bytes().all(|b| b.is_ascii_digit())),
                    "{}: bad accession {a}",
                    r.name
                );
                assert!(acc.insert(a.clone()), "duplicate XLMOD accession {a}");
            }
        }
        for n in table().names() {
            assert!(names.insert(n.to_ascii_lowercase()), "duplicate name {n}");
        }
    }

    #[test]
    fn every_unimod_backed_mass_matches_its_composition() {
        let db = db();
        let mut checked = 0;
        for r in &table().reagents {
            for p in parts(r) {
                let Some(name) = &p.unimod else { continue };
                let rec =
                    db.resolve(name).ok().flatten().unwrap_or_else(|| {
                        panic!("{} {}: no UNIMOD record {name}", r.name, p.name)
                    });
                let m = mass_of(p);
                assert!(
                    (rec.mono_mass - m).abs() < 2e-5,
                    "{} {} {name}: UNIMOD {} vs composition {m}",
                    r.name,
                    p.name,
                    rec.mono_mass
                );
                checked += 1;
            }
        }
        assert_eq!(checked, 36, "UNIMOD-backed masses checked");
    }

    #[test]
    fn every_intact_unimod_xlink_record_is_in_the_table_or_named_here() {
        // UNIMOD Xlink records that are not a reagent's intact cross-link,
        // stub or monolink, and why.
        const LEFT_OUT: &[(&str, &str)] = &[
            (
                "Xlink:B10621",
                "iodoacetyl reagent whose composition keeps the iodine",
            ),
            ("Xlink:DTSSP[88]", "chemical cleavage product (reduction)"),
            ("Xlink:DTBP[87]", "chemical cleavage product (reduction)"),
            (
                "Xlink:EGS[115]",
                "chemical cleavage product (hydroxylamine)",
            ),
            ("Xlink:DST[56]", "chemical cleavage product (periodate)"),
            (
                "Xlink:DMP[154]",
                "free monolink other than the quenched forms",
            ),
            ("Xlink:SMCC[321]", "monolink quenched by an amine reagent"),
        ];
        let db = db();
        for acc in 1..=3000 {
            let Some(rec) = db.by_accession(&format!("UNIMOD:{acc}")) else {
                continue;
            };
            if !rec.name.starts_with("Xlink:") {
                continue;
            }
            let listed = table()
                .reagents
                .iter()
                .flat_map(parts)
                .any(|p| p.unimod.as_deref() == Some(rec.name.as_str()));
            assert!(
                listed || LEFT_OUT.iter().any(|(n, _)| *n == rec.name),
                "UNIMOD {} ({}) is neither in the reagent table nor left out with a reason",
                rec.name,
                rec.accession
            );
        }
    }

    #[test]
    fn cleavable_stubs_split_the_intact_linker() {
        let h2o = composition_mass("H(2) O").unwrap();
        for r in table().reagents.iter().filter(|r| !r.stubs.is_empty()) {
            let intact = mass_of(&r.intact);
            let stubs: Vec<f64> = r.stubs.iter().map(mass_of).collect();
            // Some pair of stubs (the light doublet stub and a heavy form, or
            // the heavy form plus water) adds up to the intact linker.
            let splits = stubs.iter().enumerate().any(|(i, a)| {
                stubs[i..]
                    .iter()
                    .any(|b| (a + b - intact).abs() < 1e-6 || (a + b + h2o - intact).abs() < 1e-6)
            });
            assert!(
                splits,
                "{}: no two stubs add up to the intact linker",
                r.name
            );
            assert!(
                !r.doublets.is_empty(),
                "{}: cleavable without doublets",
                r.name
            );
            let idx = |n: &str| r.stubs.iter().position(|s| s.name == n).unwrap();
            for [l, h] in &r.doublets {
                assert!(
                    stubs[idx(l)] < stubs[idx(h)],
                    "{}: doublet not light-heavy",
                    r.name
                );
            }
            assert!(r.stubs.len() <= MAX_STUBS);
        }
    }

    #[test]
    fn monolinks_are_the_intact_linker_plus_the_quencher() {
        let q = |name: &str| match name {
            "water" => composition_mass("H(2) O"),
            "ammonia" => composition_mass("H(3) N"),
            "tris" => composition_mass("H(11) C(4) N O(3)"),
            _ => None,
        };
        for r in &table().reagents {
            let intact = mass_of(&r.intact);
            for p in &r.monolinks {
                let quencher =
                    q(&p.name).unwrap_or_else(|| panic!("{}: quencher {}", r.name, p.name));
                let m = mass_of(p);
                assert!(
                    (m - intact - quencher).abs() < 1e-6,
                    "{} {}: {m} != {intact} + {quencher}",
                    r.name,
                    p.name
                );
            }
        }
    }

    #[test]
    fn published_masses_of_reagents_unimod_lacks() {
        // Masses as the cited papers give them (rounded there).
        let t = table();
        let m = |name: &str| mass_of(&t.find(name).unwrap().intact);
        assert!((m("DSBSO") - 308.0388).abs() < 1e-3);
        assert!((m("PhoX") - 209.9718).abs() < 1e-3);
        assert!((m("CDI") - 25.9793).abs() < 1e-3);
        assert!((m("SDA") - 82.0419).abs() < 1e-3);
        assert!((m("EDC") + 18.0106).abs() < 1e-3);
        let dsbso = t.find("dsbso").unwrap();
        let stub = |k: usize| mass_of(&dsbso.stubs[k]);
        // Kaake et al.: alkene +54, unsaturated thiol +236, sulfenic acid +254.
        assert!((stub(0) - 54.0106).abs() < 1e-3);
        assert!((stub(1) - 236.0177).abs() < 1e-3);
        assert!((stub(2) - 254.0283).abs() < 1e-3);
    }

    #[test]
    fn lookup_by_name_alias_unimod_record_and_xlmod() {
        let t = table();
        assert_eq!(t.find("dsso").unwrap().name, "DSSO");
        assert_eq!(t.find("BS3").unwrap().name, "DSS");
        assert_eq!(t.find("buurbu").unwrap().name, "DSBU");
        assert_eq!(t.find("Xlink:DSS[138]").unwrap().name, "DSS");
        assert!(t.find("Xlink:DSS[156]").is_none());
        assert_eq!(t.owner_of_record("Xlink:DSS[156]").unwrap().0.name, "DSS");
        assert_eq!(t.find_xlmod("XLMOD:02126").unwrap().name, "DSSO");
        assert_eq!(t.find_xlmod("xlmod:02000").unwrap().name, "DSS");
        assert_eq!(t.find_xlmod("XLMOD:02010").unwrap().name, "EDC");
        assert!(t.find_xlmod("XLMOD:99999").is_none());
    }

    #[test]
    fn monolinks_resolve_by_reagent_and_quencher() {
        let db = db();
        let m = monolink(db, "dsbso:WATER").unwrap();
        assert_eq!(m.name, "DSBSO:water");
        assert!(m.accession.is_empty());
        assert!((m.mono_mass - (308.038_831 + 18.010_565)).abs() < 1e-5);
        let m = monolink(db, "PhoX:tris").unwrap();
        assert_eq!(m.name, "PhoX:tris");
        // A UNIMOD-backed monolink is its UNIMOD record.
        let m = monolink(db, "BS3:water").unwrap();
        assert_eq!(
            (m.name.as_str(), m.accession.as_str()),
            ("Xlink:DSS[156]", "UNIMOD:1020")
        );
        assert!(monolink(db, "DSSO:acid").is_none());
        assert!(monolink(db, "Oxidation").is_none());
        assert!(monolink(db, "NoSuch:water").is_none());
    }

    #[test]
    fn a_custom_entry_parses_with_plain_masses() {
        let t = ReagentTable::parse_yaml(
            "- name: MyLinker\n  sites: K\n  intact: {mass: 100.5}\n  \
             stubs: [{name: a, mass: 10.0}, {name: b, mass: 20.0}]\n  doublets: [[a, b]]\n",
        )
        .unwrap();
        let r = &t.reagents[0];
        assert_eq!(r.intact.resolve(db()).unwrap(), (100.5, None));
        assert_eq!(r.doublets, vec![["a".to_string(), "b".to_string()]]);
        assert!(ReagentTable::parse_yaml("- name: X\n  intact: {}\n  stubz: []\n").is_err());
        let p = Part::default();
        assert!(p.resolve(db()).unwrap_err().contains("needs a `unimod`"));
    }

    #[test]
    fn composition_parser() {
        assert_eq!(composition_mass(""), Some(0.0));
        assert!((composition_mass("H(2) O").unwrap() - 18.010_565).abs() < 1e-6);
        assert!((composition_mass("H(-2) O(-1)").unwrap() + 18.010_565).abs() < 1e-6);
        assert_eq!(composition_mass("Xx(2)"), None);
        assert_eq!(composition_mass("H(2"), None);
    }
}
