//! Curated registry of experiment-class slugs with alias folding and mod-inference.
//!
//! A model is tagged with a **set** of slugs (e.g. `phospho+tmt`) produced by
//! [`ProtocolCatalog::canonical_set`].  The registry also supports inferring the
//! experiment class from configured modification masses via
//! [`ProtocolCatalog::infer_from_mods`].

/// One entry in the built-in protocol catalog.
struct CatalogEntry {
    /// Canonical lowercase slug (e.g. `"phospho"`).
    slug: &'static str,
    /// Extra lowercase aliases (the slug itself is always accepted).
    aliases: &'static [&'static str],
    /// Optional mod-inference signature: `(mass, tolerance, residue_chars)`.
    /// `residue_chars` is a `&str` whose chars are accepted residue letters;
    /// `'*'` means any residue (wildcard).
    inference: Option<(f64, f64, &'static str)>,
}

/// Curated registry of experiment-class slugs.
///
/// # Slug sets
/// Slugs can be combined; the canonical representation is a sorted, `+`-joined
/// string like `"phospho+tmt"`.  The empty set maps to `"standard"`.
pub struct ProtocolCatalog {
    entries: Vec<CatalogEntry>,
}

impl ProtocolCatalog {
    /// Returns the built-in catalog covering the most common proteomics
    /// experiment classes.
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
                    aliases: &[
                        "phosphorylation",
                        "phospho-enrichment",
                        "phosphorylated",
                    ],
                    // +79.96633 Da on Ser / Thr / Tyr
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
                    // +229.1629 Da on Lys (or peptide N-term — matched by 'K')
                    inference: Some((229.1629, 0.01, "K")),
                },
                CatalogEntry {
                    slug: "itraq",
                    aliases: &["itraq4", "itraq8"],
                    // +144.1021 or +304.2054 on K/N-term.
                    // We register the lighter variant here; the heavier one is
                    // handled by a second entry below.
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
                    // +42.0106 Da on Lys
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
                    // +114.0429 Da on Lys (di-Gly remnant)
                    inference: Some((114.0429, 0.01, "K")),
                },
                CatalogEntry {
                    slug: "glyco",
                    aliases: &["glycosylation", "glycoproteomics", "glycan"],
                    // No single-mass signature — explicit tag only.
                    inference: None,
                },
                CatalogEntry {
                    slug: "immuno",
                    aliases: &["immunopeptidomics", "hla", "mhc"],
                    // No mass signature — explicit tag only.
                    inference: None,
                },
            ],
        }
    }

    /// Returns the canonical slug for `name`, or `None` if unrecognised.
    ///
    /// Matching is case-insensitive and trims leading/trailing whitespace.
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

    /// Folds each item in `slugs` through [`canonical`](Self::canonical),
    /// drops unknowns, deduplicates, sorts, and joins with `+`.
    ///
    /// * If the resulting set is empty, returns `"standard"`.
    /// * If the set contains `"standard"` alongside other slugs, `"standard"`
    ///   is dropped (it is the no-op base class).
    pub fn canonical_set<I, S>(&self, slugs: I) -> String
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        let mut set: Vec<String> = slugs
            .into_iter()
            .filter_map(|s| self.canonical(s.as_ref()))
            .collect();

        // Dedup + sort (stable).
        set.sort_unstable();
        set.dedup();

        // Drop the "standard" no-op when other classes are present.
        if set.len() > 1 {
            set.retain(|s| s != "standard");
        }

        if set.is_empty() {
            "standard".to_string()
        } else {
            set.join("+")
        }
    }

    /// Infers the experiment class from a slice of `(mass, residue)` pairs
    /// (typically parsed from a `mods.txt`).
    ///
    /// Each pair is matched against catalog entries that carry a mass-inference
    /// signature.  Matching slugs are passed through [`canonical_set`](Self::canonical_set).
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
        assert_eq!(c.canonical("phospho-enrichment"), Some("phospho".to_string()));
        assert_eq!(c.canonical("TMT"), Some("tmt".to_string()));
        assert_eq!(c.canonical("phospho"), Some("phospho".to_string())); // already canonical
        assert_eq!(c.canonical("nonsense"), None);
    }

    #[test]
    fn canonical_set_is_sorted_and_joined() {
        let c = ProtocolCatalog::builtin();
        assert_eq!(c.canonical_set(["tmt", "phospho"]), "phospho+tmt");
        assert_eq!(c.canonical_set(["phospho", "tmt"]), "phospho+tmt"); // order-independent
        assert_eq!(c.canonical_set::<[&str; 0], &str>([]), "standard"); // empty -> "standard"
    }

    #[test]
    fn infers_classes_from_mod_masses() {
        let c = ProtocolCatalog::builtin();
        // (mass, residue) pairs from a configured mods.txt
        let mods = [(229.1629_f64, 'K'), (79.96633_f64, 'S')];
        assert_eq!(c.infer_from_mods(&mods), "phospho+tmt");
        let none = [(15.9949_f64, 'M')]; // oxidation -> not an experiment class
        assert_eq!(c.infer_from_mods(&none), "standard");
    }

    #[test]
    fn canonical_set_drops_standard_when_other_present() {
        let c = ProtocolCatalog::builtin();
        // "unlabeled" folds to "standard"; combined with "phospho" it should be dropped.
        assert_eq!(c.canonical_set(["phospho", "unlabeled"]), "phospho");
    }

    #[test]
    fn canonical_set_deduplicates() {
        let c = ProtocolCatalog::builtin();
        assert_eq!(c.canonical_set(["phospho", "phosphorylation", "tmt"]), "phospho+tmt");
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
        assert_eq!(c.canonical("  tmt  "), Some("tmt".to_string())); // trims whitespace
    }
}
