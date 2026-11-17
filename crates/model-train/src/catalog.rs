struct CatalogEntry {
    slug: &'static str,
    aliases: &'static [&'static str],
    inference: Option<(f64, f64, &'static str)>,
}

pub struct ProtocolCatalog {
    entries: Vec<CatalogEntry>,
}

impl ProtocolCatalog {
    pub fn builtin() -> Self {
        Self {
            entries: vec![
                CatalogEntry {
                    slug: "standard",
                    aliases: &["unlabeled", "none"],
                    inference: None,
                },
                CatalogEntry {
                    slug: "phospho",
                    aliases: &["phosphorylation", "phospho-enrichment", "phosphorylated"],
                    inference: Some((79.96633, 0.01, "STY")),
                },
                CatalogEntry {
                    slug: "tmt",
                    aliases: &[
                        "tmt6",
                        "tmt10",
                        "tmt11",
                        "tmt16",
                        "tmtpro",
                        "tandem mass tag",
                    ],
                    inference: Some((229.1629, 0.01, "K")),
                },
                CatalogEntry {
                    slug: "itraq",
                    aliases: &["itraq4", "itraq8"],
                    inference: Some((144.1021, 0.01, "K")),
                },
                CatalogEntry {
                    slug: "itraq",
                    aliases: &[],
                    inference: Some((304.2054, 0.01, "K")),
                },
                CatalogEntry {
                    slug: "acetyl",
                    aliases: &["acetylation", "acetylome"],
                    inference: Some((42.0106, 0.01, "K")),
                },
                CatalogEntry {
                    slug: "ubiquitin",
                    aliases: &[
                        "ubiquitination",
                        "diglycine",
                        "di-gly",
                        "gg",
                        "ubiquitinome",
                    ],
                    inference: Some((114.0429, 0.01, "K")),
                },
                CatalogEntry {
                    slug: "glyco",
                    aliases: &["glycosylation", "glycoproteomics", "glycan"],
                    inference: None,
                },
                CatalogEntry {
                    slug: "immuno",
                    aliases: &["immunopeptidomics", "hla", "mhc"],
                    inference: None,
                },
            ],
        }
    }

    pub fn canonical(&self, name: &str) -> Option<String> {
        let key = name.trim().to_lowercase();
        for entry in &self.entries {
            if entry.slug == key {
                return Some(entry.slug.to_string());
            }
            if entry.aliases.iter().any(|a| *a == key) {
                return Some(entry.slug.to_string());
            }
        }
        None
    }

    pub fn canonical_set<I, S>(&self, slugs: I) -> String
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        let mut set: Vec<String> = slugs
            .into_iter()
            .filter_map(|s| self.canonical(s.as_ref()))
            .collect();

        set.sort_unstable();
        set.dedup();

        if set.len() > 1 {
            set.retain(|s| s != "standard");
        }

        if set.is_empty() {
            "standard".to_string()
        } else {
            set.join("+")
        }
    }

    pub fn infer_from_mods(&self, mods: &[(f64, char)]) -> String {
        let matched: Vec<String> = self
            .entries
            .iter()
            .filter_map(|entry| {
                let (sig_mass, tol, residues) = entry.inference?;
                let matches = mods.iter().any(|(mass, residue)| {
                    (mass - sig_mass).abs() <= tol
                        && (residues == "*" || residues.contains(*residue))
                });
                if matches {
                    Some(entry.slug.to_string())
                } else {
                    None
                }
            })
            .collect();

        self.canonical_set(matched)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn folds_aliases_to_canonical_slugs() {
        let c = ProtocolCatalog::builtin();
        assert_eq!(c.canonical("Phosphorylation"), Some("phospho".to_string()));
        assert_eq!(
            c.canonical("phospho-enrichment"),
            Some("phospho".to_string())
        );
        assert_eq!(c.canonical("TMT"), Some("tmt".to_string()));
        assert_eq!(c.canonical("phospho"), Some("phospho".to_string()));
        assert_eq!(c.canonical("nonsense"), None);
    }

    #[test]
    fn canonical_set_is_sorted_and_joined() {
        let c = ProtocolCatalog::builtin();
        assert_eq!(c.canonical_set(["tmt", "phospho"]), "phospho+tmt");
        assert_eq!(c.canonical_set(["phospho", "tmt"]), "phospho+tmt");
        assert_eq!(c.canonical_set::<[&str; 0], &str>([]), "standard");
    }

    #[test]
    fn infers_classes_from_mod_masses() {
        let c = ProtocolCatalog::builtin();
        let mods = [(229.1629_f64, 'K'), (79.96633_f64, 'S')];
        assert_eq!(c.infer_from_mods(&mods), "phospho+tmt");
        let none = [(15.9949_f64, 'M')];
        assert_eq!(c.infer_from_mods(&none), "standard");
    }

    #[test]
    fn canonical_set_drops_standard_when_other_present() {
        let c = ProtocolCatalog::builtin();
        assert_eq!(c.canonical_set(["phospho", "unlabeled"]), "phospho");
    }

    #[test]
    fn canonical_set_deduplicates() {
        let c = ProtocolCatalog::builtin();
        assert_eq!(
            c.canonical_set(["phospho", "phosphorylation", "tmt"]),
            "phospho+tmt"
        );
    }

    #[test]
    fn infers_itraq_from_light_mass() {
        let c = ProtocolCatalog::builtin();
        let mods = [(144.1021_f64, 'K')];
        assert_eq!(c.infer_from_mods(&mods), "itraq");
    }

    #[test]
    fn infers_ubiquitin_from_diglycine() {
        let c = ProtocolCatalog::builtin();
        let mods = [(114.0429_f64, 'K')];
        assert_eq!(c.infer_from_mods(&mods), "ubiquitin");
    }

    #[test]
    fn canonical_aliases_case_insensitive() {
        let c = ProtocolCatalog::builtin();
        assert_eq!(c.canonical("GLYCO"), Some("glyco".to_string()));
        assert_eq!(c.canonical("Glycosylation"), Some("glyco".to_string()));
        assert_eq!(c.canonical("HLA"), Some("immuno".to_string()));
        assert_eq!(c.canonical("  tmt  "), Some("tmt".to_string()));
    }
}
