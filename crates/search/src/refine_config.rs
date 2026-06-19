//! `--refine-config` YAML -> refinement modification spec. A built-in DEFAULT
//! tier (5 mods, X!Tandem always-on chemistry) is used when no config is given.

use serde::Deserialize;

/// One refinement variable modification.
#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct RefineMod {
    pub name: String,
    pub delta: f64,
    pub residues: Vec<String>, // e.g. ["M"], ["N","Q"], ["*"]
    #[serde(default = "loc_anywhere")]
    pub location: String, // anywhere | n_term | c_term | protein_n_term
    pub class: String,    // oxidation | deamidation | nterm_acetyl | nterm_loss | alkyl | other
}
fn loc_anywhere() -> String {
    "anywhere".into()
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct RefineConfig {
    #[serde(default = "default_max_mods")]
    pub max_mods: u32,
    #[serde(default = "default_high_res_only")]
    pub high_res_only: bool,
    pub mods: Vec<RefineMod>,
}
fn default_max_mods() -> u32 {
    2
}
fn default_high_res_only() -> bool {
    true
}

impl RefineConfig {
    /// The built-in DEFAULT tier: the X!Tandem always-on chemistry. Used when
    /// `--refine` is given without `--refine-config`.
    pub fn default_tier() -> Self {
        let m = |name: &str, delta: f64, res: &[&str], loc: &str, class: &str| RefineMod {
            name: name.into(),
            delta,
            residues: res.iter().map(|s| s.to_string()).collect(),
            location: loc.into(),
            class: class.into(),
        };
        RefineConfig {
            max_mods: 2,
            high_res_only: true,
            mods: vec![
                // Oxidation/hydroxylation (+15.995, Unimod 35). M is the common
                // artifact; P and K cover collagen/ECM hydroxyproline &
                // hydroxylysine — heavily-modified proteins whose peptides never
                // match unmodified, so a closed search loses the whole protein.
                // (P+K widen the oxidation candidate space ~3-5× vs M-only.)
                m("Oxidation", 15.994915, &["M", "P", "K"], "anywhere", "oxidation"),
                m("Deamidation", 0.984016, &["N", "Q"], "anywhere", "deamidation"),
                m("Gln->pyro-Glu", -17.026549, &["Q"], "n_term", "nterm_loss"),
                m("Glu->pyro-Glu", -18.010565, &["E"], "n_term", "nterm_loss"),
                m("Acetyl", 42.010565, &["*"], "protein_n_term", "nterm_acetyl"),
            ],
        }
    }

    pub fn from_yaml_str(s: &str) -> Result<Self, String> {
        let cfg: RefineConfig =
            serde_yaml::from_str(s).map_err(|e| format!("refine-config parse: {e}"))?;
        cfg.validate()?;
        Ok(cfg)
    }

    /// Reject invalid residue/location tokens that would otherwise SILENTLY become
    /// broad defaults (a residue typo → wildcard = every residue; an unknown
    /// location → Anywhere), massively (and wrongly) expanding the candidate space.
    /// A residue must be `"*"` or a single uppercase amino-acid letter; a location
    /// must be one of the accepted spellings (those `refinement::parse_location`
    /// recognizes). Fails fast with a clear message.
    pub fn validate(&self) -> Result<(), String> {
        // Accepted location spellings (lowercased), matching parse_location +
        // the `anywhere` default.
        const VALID_LOCATIONS: &[&str] = &[
            "anywhere",
            "n_term", "n-term", "nterm",
            "c_term", "c-term", "cterm",
            "protein_n_term", "prot-n-term", "prot_n_term",
            "protein_c_term", "prot-c-term", "prot_c_term",
        ];
        for m in &self.mods {
            for r in &m.residues {
                let ok = r == "*"
                    || (r.len() == 1 && r.as_bytes()[0].is_ascii_uppercase());
                if !ok {
                    return Err(format!(
                        "refine-config mod '{}': invalid residue token {:?} \
                         (expected \"*\" or a single uppercase amino-acid letter)",
                        m.name, r
                    ));
                }
            }
            let loc = m.location.trim().to_ascii_lowercase();
            if !VALID_LOCATIONS.contains(&loc.as_str()) {
                return Err(format!(
                    "refine-config mod '{}': unknown location {:?} (expected one of \
                     anywhere | n_term | c_term | protein_n_term | protein_c_term)",
                    m.name, m.location
                ));
            }
        }
        Ok(())
    }

    /// Map a class string to the PIN `refine_mod_class` id (Task 2 encoding).
    pub fn class_id(class: &str) -> u32 {
        match class {
            "oxidation" => 1,
            "deamidation" => 2,
            "nterm_acetyl" => 3,
            "nterm_loss" => 4,
            "alkyl" => 5,
            _ => 99,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_tier_has_five_mods() {
        let t = RefineConfig::default_tier();
        assert_eq!(t.mods.len(), 5);
        assert_eq!(t.max_mods, 2);
        assert!(t.high_res_only);
        assert!(t.mods.iter().any(|m| m.name == "Deamidation" && m.class == "deamidation"));
    }

    #[test]
    fn yaml_round_trip() {
        let y = "max_mods: 2\nhigh_res_only: true\nmods:\n  - {name: Oxidation, delta: 15.994915, residues: [M], location: anywhere, class: oxidation}\n";
        let c = RefineConfig::from_yaml_str(y).unwrap();
        assert_eq!(c.mods.len(), 1);
        assert_eq!(c.mods[0].delta, 15.994915);
    }

    #[test]
    fn class_ids() {
        assert_eq!(RefineConfig::class_id("deamidation"), 2);
        assert_eq!(RefineConfig::class_id("unknown"), 99);
    }

    #[test]
    fn default_tier_validates() {
        assert!(RefineConfig::default_tier().validate().is_ok());
    }

    #[test]
    fn validate_rejects_bad_residue_token() {
        let y = "mods:\n  - {name: Oxidation, delta: 15.994915, residues: [Met], location: anywhere, class: oxidation}\n";
        let err = RefineConfig::from_yaml_str(y).unwrap_err();
        assert!(err.contains("invalid residue token"), "got: {err}");
    }

    #[test]
    fn validate_rejects_unknown_location() {
        let y = "mods:\n  - {name: Acetyl, delta: 42.010565, residues: [\"*\"], location: ProteinNTerm, class: nterm_acetyl}\n";
        let err = RefineConfig::from_yaml_str(y).unwrap_err();
        assert!(err.contains("unknown location"), "got: {err}");
    }

    #[test]
    fn validate_accepts_star_and_single_letter() {
        let y = "mods:\n  - {name: Acetyl, delta: 42.010565, residues: [\"*\"], location: protein_n_term, class: nterm_acetyl}\n  - {name: Oxidation, delta: 15.994915, residues: [M], location: anywhere, class: oxidation}\n";
        assert!(RefineConfig::from_yaml_str(y).is_ok());
    }
}
