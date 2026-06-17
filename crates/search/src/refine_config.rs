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
                m("Oxidation", 15.994915, &["M"], "anywhere", "oxidation"),
                m("Deamidation", 0.984016, &["N", "Q"], "anywhere", "deamidation"),
                m("Gln->pyro-Glu", -17.026549, &["Q"], "n_term", "nterm_loss"),
                m("Glu->pyro-Glu", -18.010565, &["E"], "n_term", "nterm_loss"),
                m("Acetyl", 42.010565, &["*"], "protein_n_term", "nterm_acetyl"),
            ],
        }
    }

    pub fn from_yaml_str(s: &str) -> Result<Self, String> {
        serde_yaml::from_str(s).map_err(|e| format!("refine-config parse: {e}"))
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
}
