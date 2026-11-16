//! The crosslinking reagents andes knows by name.
//!
//! Every mass of a reagent (intact cross-link, MS-cleavage stubs, quenched
//! monolinks) is an elemental composition. When UNIMOD has the record
//! (`Xlink:<NAME>[<nominal>]`), the record names it and its UNIMOD mass is
//! the one used; a test checks that each composition gives the UNIMOD mass.
//! Reagents UNIMOD lacks carry the composition from the cited paper, and their
//! mass is computed from it.
//!
//! Sites are the reagent's usual reactivity, written as in a spec (`K,
//! Protein N-term`; `end / end` for a heterobifunctional reagent). NHS esters
//! default to Lys and the protein N-terminus; Ser, Thr and Tyr, which they
//! also acylate at a lower rate, are opt-in (`DSSO (K, S, T, Y, Protein
//! N-term)`).

/// One mass of a reagent.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Part {
    /// Column or display name (a stub's short name, a monolink's quencher).
    pub name: &'static str,
    /// Elemental composition in UNIMOD notation (`H(6) C(6) O(3) S`;
    /// negative counts allowed).
    pub composition: &'static str,
    /// The UNIMOD record of this mass, if UNIMOD has it.
    pub unimod: Option<&'static str>,
}

/// A crosslinking reagent of the built-in table.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Reagent {
    /// Canonical name.
    pub name: &'static str,
    /// Other names accepted for it.
    pub aliases: &'static [&'static str],
    /// The intact cross-link.
    pub intact: Part,
    /// MS-cleavage stubs; empty for a non-cleavable reagent.
    pub stubs: &'static [Part],
    /// Reporter doublets as `(light, heavy)` stub indices.
    pub doublets: &'static [(usize, usize)],
    /// Quenched monolinks.
    pub monolinks: &'static [Part],
    /// Default sites, as in a spec.
    pub sites: &'static str,
    /// Where the chemistry and masses come from.
    pub source: &'static str,
}

const fn p(name: &'static str, composition: &'static str, unimod: Option<&'static str>) -> Part {
    Part {
        name,
        composition,
        unimod,
    }
}

const NHS: &str = "K, Protein N-term";

/// The built-in reagents.
pub const REAGENTS: &[Reagent] = &[
    // MS-cleavable.
    Reagent {
        name: "DSSO",
        aliases: &[],
        intact: p("intact", "H(6) C(6) O(3) S", Some("Xlink:DSSO[158]")),
        stubs: &[
            p("alkene", "H(2) C(3) O", Some("Xlink:DSSO[54]")),
            p("thiol", "H(2) C(3) O S", Some("Xlink:DSSO[86]")),
            p("sulfenic", "H(4) C(3) O(2) S", Some("Xlink:DSSO[104]")),
        ],
        doublets: &[(0, 1)],
        monolinks: &[
            p("water", "H(8) C(6) O(4) S", Some("Xlink:DSSO[176]")),
            p("ammonia", "H(9) C(6) N O(3) S", Some("Xlink:DSSO[175]")),
            p("tris", "H(17) C(10) N O(6) S", Some("Xlink:DSSO[279]")),
        ],
        sites: NHS,
        source: "Kao et al. 2011, Mol Cell Proteomics 10:M110.002212; UNIMOD",
    },
    Reagent {
        name: "DSBU",
        aliases: &["BuUrBu"],
        intact: p("intact", "H(12) C(9) N(2) O(3)", Some("Xlink:BuUrBu[196]")),
        stubs: &[
            p("Bu", "H(7) C(4) N O", Some("Xlink:BuUrBu[85]")),
            p("BuUr", "H(5) C(5) N O(2)", Some("Xlink:BuUrBu[111]")),
        ],
        doublets: &[(0, 1)],
        monolinks: &[
            p("water", "H(14) C(9) N(2) O(4)", Some("Xlink:BuUrBu[214]")),
            p("ammonia", "H(15) C(9) N(3) O(3)", Some("Xlink:BuUrBu[213]")),
            p("tris", "H(23) C(13) N(3) O(6)", Some("Xlink:BuUrBu[317]")),
        ],
        sites: NHS,
        source: "Mueller et al. 2010, Anal Chem 82:6958; UNIMOD",
    },
    Reagent {
        name: "DSBSO",
        aliases: &["Azide-A-DSBSO"],
        intact: p("intact", "H(16) C(11) O(6) S(2)", None),
        stubs: &[
            p("alkene", "H(2) C(3) O", Some("Xlink:DSSO[54]")),
            p("thiol", "H(12) C(8) O(4) S(2)", None),
            p("sulfenic", "H(14) C(8) O(5) S(2)", None),
        ],
        doublets: &[(0, 1)],
        monolinks: &[
            p("water", "H(18) C(11) O(7) S(2)", None),
            p("ammonia", "H(19) C(11) N O(6) S(2)", None),
        ],
        sites: NHS,
        source: "Kaake et al. 2014, Mol Cell Proteomics 13:3533; compositions as in \
                 Matzinger et al. 2020, J Proteome Res 19:2071",
    },
    Reagent {
        name: "CDI",
        aliases: &["carbonyldiimidazole"],
        intact: p("intact", "H(-2) C O", None),
        stubs: &[p("amine", "", None), p("isocyanate", "H(-2) C O", None)],
        doublets: &[(0, 1)],
        monolinks: &[],
        sites: "K, S, T, Y, Protein N-term",
        source: "Hage et al. 2017, Angew Chem Int Ed 56:14551",
    },
    // Non-cleavable, amine-reactive.
    Reagent {
        name: "DSS",
        aliases: &["BS3"],
        intact: p("intact", "H(10) C(8) O(2)", Some("Xlink:DSS[138]")),
        stubs: &[],
        doublets: &[],
        monolinks: &[
            p("water", "H(12) C(8) O(3)", Some("Xlink:DSS[156]")),
            p("ammonia", "H(13) C(8) N O(2)", Some("Xlink:DSS[155]")),
            p("tris", "H(21) C(12) N O(5)", Some("Xlink:DSS[259]")),
        ],
        sites: NHS,
        source: "UNIMOD",
    },
    Reagent {
        name: "BS2G",
        aliases: &["DSG"],
        intact: p("intact", "H(4) C(5) O(2)", Some("Xlink:BS2G[96]")),
        stubs: &[],
        doublets: &[],
        monolinks: &[
            p("water", "H(6) C(5) O(3)", Some("Xlink:BS2G[114]")),
            p("ammonia", "H(7) C(5) N O(2)", Some("Xlink:BS2G[113]")),
            p("tris", "H(15) C(9) N O(5)", Some("Xlink:BS2G[217]")),
        ],
        sites: NHS,
        source: "UNIMOD",
    },
    Reagent {
        name: "EGS",
        aliases: &["sulfo-EGS"],
        intact: p("intact", "H(10) C(10) O(6)", Some("Xlink:EGS[226]")),
        stubs: &[],
        doublets: &[],
        monolinks: &[p("water", "H(12) C(10) O(7)", Some("Xlink:EGS[244]"))],
        sites: NHS,
        source: "UNIMOD",
    },
    Reagent {
        name: "DST",
        aliases: &["sulfo-DST"],
        intact: p("intact", "H(2) C(4) O(4)", Some("Xlink:DST[114]")),
        stubs: &[],
        doublets: &[],
        monolinks: &[p("water", "H(4) C(4) O(5)", Some("Xlink:DST[132]"))],
        sites: NHS,
        source: "UNIMOD",
    },
    Reagent {
        name: "DSP",
        aliases: &["DTSSP"],
        intact: p("intact", "H(6) C(6) O(2) S(2)", Some("Xlink:DTSSP[174]")),
        stubs: &[],
        doublets: &[],
        monolinks: &[p("water", "H(8) C(6) O(3) S(2)", Some("Xlink:DTSSP[192]"))],
        sites: NHS,
        source: "UNIMOD",
    },
    Reagent {
        name: "PhoX",
        aliases: &["DSPP"],
        intact: p("intact", "H(3) C(8) O(5) P", None),
        stubs: &[],
        doublets: &[],
        monolinks: &[
            p("water", "H(5) C(8) O(6) P", None),
            p("ammonia", "H(6) C(8) N O(5) P", None),
            p("tris", "H(14) C(12) N O(8) P", None),
        ],
        sites: NHS,
        source: "Steigenberger et al. 2019, ACS Cent Sci 5:1514",
    },
    Reagent {
        name: "DMP",
        aliases: &[],
        intact: p("intact", "H(10) C(7) N(2)", Some("Xlink:DMP[122]")),
        stubs: &[],
        doublets: &[],
        monolinks: &[
            p("water", "H(12) C(7) N(2) O", Some("Xlink:DMP[140]")),
            p("ammonia", "H(13) C(7) N(3)", Some("Xlink:DMP[139]")),
        ],
        sites: NHS,
        source: "UNIMOD",
    },
    Reagent {
        name: "DTBP",
        aliases: &[],
        intact: p("intact", "H(8) C(6) N(2) S(2)", Some("Xlink:DTBP[172]")),
        stubs: &[],
        doublets: &[],
        monolinks: &[],
        sites: NHS,
        source: "UNIMOD",
    },
    Reagent {
        name: "DFDNB",
        aliases: &[],
        intact: p("intact", "C(6) N(2) O(4)", Some("Xlink:DFDNB")),
        stubs: &[],
        doublets: &[],
        monolinks: &[],
        sites: NHS,
        source: "UNIMOD",
    },
    // Non-cleavable, other chemistries.
    Reagent {
        name: "BMOE",
        aliases: &[],
        intact: p("intact", "H(8) C(10) N(2) O(4)", Some("Xlink:BMOE")),
        stubs: &[],
        doublets: &[],
        monolinks: &[],
        sites: "C",
        source: "UNIMOD",
    },
    Reagent {
        name: "SMCC",
        aliases: &["sulfo-SMCC"],
        intact: p("intact", "H(13) C(12) N O(3)", Some("Xlink:SMCC[219]")),
        stubs: &[],
        doublets: &[],
        monolinks: &[p("water", "H(15) C(12) N O(4)", Some("Xlink:SMCC[237]"))],
        sites: "K, Protein N-term / C",
        source: "UNIMOD",
    },
    Reagent {
        name: "SDA",
        aliases: &["sulfo-SDA"],
        intact: p("intact", "H(6) C(5) O", None),
        stubs: &[],
        doublets: &[],
        monolinks: &[p("water", "H(8) C(5) O(2)", None)],
        sites: "K, S, T, Y, Protein N-term / ACDEFGHIKLMNPQRSTVWY",
        source: "Belsom et al. 2016, Mol Cell Proteomics 15:1105",
    },
    Reagent {
        name: "EDC",
        aliases: &["DMTMM"],
        intact: p("intact", "H(-2) O(-1)", None),
        stubs: &[],
        doublets: &[],
        monolinks: &[],
        sites: "K, Protein N-term / D, E, Protein C-term",
        source: "Leitner et al. 2014, Proc Natl Acad Sci USA 111:9455",
    },
];

/// The reagent named `name`: a canonical name or alias (any case), or the
/// UNIMOD name or accession of its intact record.
pub fn find(name: &str) -> Option<&'static Reagent> {
    let n = name.trim();
    REAGENTS.iter().find(|r| {
        r.name.eq_ignore_ascii_case(n)
            || r.aliases.iter().any(|a| a.eq_ignore_ascii_case(n))
            || r.intact.unimod.is_some_and(|u| u.eq_ignore_ascii_case(n))
    })
}

/// The reagent with a monolink or stub record named `unimod_name`, and the
/// role of that record (`monolink` or `stub`).
pub fn owner_of_record(unimod_name: &str) -> Option<(&'static Reagent, &'static str)> {
    let is = |p: &Part| {
        p.unimod
            .is_some_and(|u| u.eq_ignore_ascii_case(unimod_name))
    };
    REAGENTS.iter().find_map(|r| {
        if r.monolinks.iter().any(is) {
            Some((r, "monolink"))
        } else if r.stubs.iter().any(is) {
            Some((r, "cleavage stub"))
        } else {
            None
        }
    })
}

/// Every name [`find`] accepts other than UNIMOD records, for suggestions.
pub fn names() -> impl Iterator<Item = &'static str> {
    REAGENTS
        .iter()
        .flat_map(|r| std::iter::once(r.name).chain(r.aliases.iter().copied()))
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

    fn db() -> &'static model::unimod::UnimodDb {
        let db = model::unimod::bundled();
        assert!(!db.is_empty(), "bundled UNIMOD table not found");
        db
    }

    fn parts(r: &Reagent) -> impl Iterator<Item = &Part> {
        std::iter::once(&r.intact)
            .chain(r.stubs.iter())
            .chain(r.monolinks.iter())
    }

    #[test]
    fn every_composition_parses() {
        for r in REAGENTS {
            for p in parts(r) {
                assert!(
                    composition_mass(p.composition).is_some(),
                    "{} {}: {}",
                    r.name,
                    p.name,
                    p.composition
                );
            }
        }
    }

    #[test]
    fn every_unimod_backed_mass_matches_its_composition() {
        let db = db();
        let mut checked = 0;
        for r in REAGENTS {
            for p in parts(r) {
                let Some(name) = p.unimod else { continue };
                let rec =
                    db.resolve(name).ok().flatten().unwrap_or_else(|| {
                        panic!("{} {}: no UNIMOD record {name}", r.name, p.name)
                    });
                let m = composition_mass(p.composition).unwrap();
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
        // UNIMOD intact cross-link records andes does not list, and why.
        const LEFT_OUT: &[&str] = &[
            // Iodoacetyl reagent whose UNIMOD composition keeps the iodine.
            "Xlink:B10621",
        ];
        let db = db();
        for acc in 1..=3000 {
            let Some(rec) = db.by_accession(&format!("UNIMOD:{acc}")) else {
                continue;
            };
            if !rec.name.starts_with("Xlink:") {
                continue;
            }
            let listed = REAGENTS
                .iter()
                .flat_map(parts)
                .any(|p| p.unimod == Some(rec.name.as_str()));
            // Chemical cleavage products and other non-cross-link forms.
            let other_product = [
                "Xlink:DTSSP[88]",
                "Xlink:DTBP[87]",
                "Xlink:EGS[115]",
                "Xlink:DST[56]",
                "Xlink:DMP[154]",
                "Xlink:SMCC[321]",
            ]
            .contains(&rec.name.as_str());
            assert!(
                listed || other_product || LEFT_OUT.contains(&rec.name.as_str()),
                "UNIMOD {} ({}) is neither in the reagent table nor left out on purpose",
                rec.name,
                rec.accession
            );
        }
    }

    #[test]
    fn cleavable_stubs_split_the_intact_linker() {
        let h2o = composition_mass("H(2) O").unwrap();
        for r in REAGENTS.iter().filter(|r| !r.stubs.is_empty()) {
            let m = |p: &Part| composition_mass(p.composition).unwrap();
            let intact = m(&r.intact);
            let stubs: Vec<f64> = r.stubs.iter().map(m).collect();
            // Some pair of stubs (the light doublet stub and a heavy form,
            // or the heavy form plus water) adds up to the intact linker.
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
            for &(l, h) in r.doublets {
                assert!(
                    stubs[l] < stubs[h],
                    "{}: doublet ({l}, {h}) not light-heavy",
                    r.name
                );
            }
            assert!(r.stubs.len() <= crate::crosslinker::MAX_STUBS);
        }
    }

    #[test]
    fn nhs_monolinks_are_the_intact_linker_plus_the_quencher() {
        let q = |name: &str| match name {
            "water" => composition_mass("H(2) O"),
            "ammonia" => composition_mass("H(3) N"),
            "tris" => composition_mass("H(11) C(4) N O(3)"),
            _ => None,
        };
        for r in REAGENTS {
            let intact = composition_mass(r.intact.composition).unwrap();
            for p in r.monolinks {
                let quencher =
                    q(p.name).unwrap_or_else(|| panic!("{}: quencher {}", r.name, p.name));
                let m = composition_mass(p.composition).unwrap();
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
        let m = |name: &str| composition_mass(find(name).unwrap().intact.composition).unwrap();
        assert!((m("DSBSO") - 308.0385).abs() < 1e-3);
        assert!((m("PhoX") - 209.9718).abs() < 1e-3);
        assert!((m("CDI") - 25.9793).abs() < 1e-3);
        assert!((m("SDA") - 82.0419).abs() < 1e-3);
        assert!((m("EDC") + 18.0106).abs() < 1e-3);
        let dsbso = find("dsbso").unwrap();
        let stub = |k: usize| composition_mass(dsbso.stubs[k].composition).unwrap();
        // Kaake et al.: alkene +54, unsaturated thiol +236, sulfenic acid +254.
        assert!((stub(0) - 54.0106).abs() < 1e-3);
        assert!((stub(1) - 236.0177).abs() < 1e-3);
        assert!((stub(2) - 254.0283).abs() < 1e-3);
    }

    #[test]
    fn lookup_by_name_alias_and_unimod_record() {
        assert_eq!(find("dsso").unwrap().name, "DSSO");
        assert_eq!(find("BS3").unwrap().name, "DSS");
        assert_eq!(find("buurbu").unwrap().name, "DSBU");
        assert_eq!(find("Xlink:DSS[138]").unwrap().name, "DSS");
        assert!(find("Xlink:DSS[156]").is_none());
        assert_eq!(owner_of_record("Xlink:DSS[156]").unwrap().0.name, "DSS");
        let mut seen = std::collections::HashSet::new();
        for n in names() {
            assert!(seen.insert(n.to_ascii_lowercase()), "duplicate name {n}");
        }
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
