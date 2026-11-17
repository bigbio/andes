use std::sync::OnceLock;

use serde::Deserialize;

use crate::unimod::{UnimodDb, UnimodMod};

pub const MAX_STUBS: usize = 4;

#[derive(Debug, Clone, Default, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Part {
    #[serde(default)]
    pub name: String,
    #[serde(default)]
    pub unimod: Option<String>,
    #[serde(default)]
    pub composition: Option<String>,
    #[serde(default)]
    pub mass: Option<f64>,
}

impl Part {
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

#[derive(Debug, Clone, Default, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Reagent {
    pub name: String,
    #[serde(default)]
    pub aliases: Vec<String>,
    #[serde(default)]
    pub xlmod: Vec<String>,
    #[serde(default)]
    pub sites: String,
    pub intact: Part,
    #[serde(default)]
    pub stubs: Vec<Part>,
    #[serde(default)]
    pub doublets: Vec<[String; 2]>,
    #[serde(default)]
    pub monolinks: Vec<Part>,
    #[serde(default)]
    pub source: Option<String>,
}

impl Reagent {
    pub fn is_named(&self, name: &str) -> bool {
        let n = name.trim();
        self.name.eq_ignore_ascii_case(n) || self.aliases.iter().any(|a| a.eq_ignore_ascii_case(n))
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct ReagentTable {
    pub reagents: Vec<Reagent>,
}

impl ReagentTable {
    pub fn parse_yaml(text: &str) -> Result<Self, String> {
        let reagents: Vec<Reagent> = serde_yaml::from_str(text).map_err(|e| e.to_string())?;
        Ok(ReagentTable { reagents })
    }

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

    pub fn find_xlmod(&self, accession: &str) -> Option<&Reagent> {
        let a = accession.trim();
        self.reagents
            .iter()
            .find(|r| r.xlmod.iter().any(|x| x.eq_ignore_ascii_case(a)))
    }

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

    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.reagents.iter().flat_map(|r| {
            std::iter::once(r.name.as_str()).chain(r.aliases.iter().map(String::as_str))
        })
    }
}

pub fn bundled() -> &'static ReagentTable {
    load().as_ref().unwrap_or_else(|_| {
        static EMPTY: OnceLock<ReagentTable> = OnceLock::new();
        EMPTY.get_or_init(ReagentTable::default)
    })
}

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
            common_sites: Vec::new(),
            isotopic_label: false,
        }),
    }
}

const ELEMENTS: &[(&str, f64)] = &[
    ("H", 1.007_825_032),
    ("C", 12.0),
    ("N", 14.003_074_004),
    ("O", 15.994_914_620),
    ("S", 31.972_071_174),
    ("P", 30.973_761_998),
];

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
        let t = table();
        let m = |name: &str| mass_of(&t.find(name).unwrap().intact);
        assert!((m("DSBSO") - 308.0388).abs() < 1e-3);
        assert!((m("PhoX") - 209.9718).abs() < 1e-3);
        assert!((m("CDI") - 25.9793).abs() < 1e-3);
        assert!((m("SDA") - 82.0419).abs() < 1e-3);
        assert!((m("EDC") + 18.0106).abs() < 1e-3);
        let dsbso = t.find("dsbso").unwrap();
        let stub = |k: usize| mass_of(&dsbso.stubs[k]);
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
