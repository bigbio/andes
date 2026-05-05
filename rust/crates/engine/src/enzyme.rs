//! Enzymatic cleavage rules. Mirrors Java
//! `edu.ucsd.msjava.msutil.Enzyme`. The 8 canonical variants are pinned by
//! `tests/enzyme_rules_match_java.rs`. Custom enzymes are deferred (see
//! Phase 1 spec § Open decisions deferred to implementation).

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Enzyme {
    Trypsin,
    Chymotrypsin,
    LysC,
    AspN,
    GluC,
    LysN,
    ArgC,
    AlphaLP,
    NoCleavage,
    NonSpecific,
}

/// Cleavage rule table — one per `Enzyme` variant.
///
/// `after`: residues whose C-terminal peptide bond is cleaved.
/// `before`: residues whose N-terminal peptide bond is cleaved.
struct EnzymeRules {
    after:  &'static [u8],
    before: &'static [u8],
    /// Special flag: NonSpecific cleaves between any pair, NoCleavage never.
    universal: Option<bool>, // Some(true) = always, Some(false) = never
}

impl Enzyme {
    fn rules(self) -> EnzymeRules {
        match self {
            Enzyme::Trypsin       => EnzymeRules { after: b"KR",      before: b"",  universal: None },
            Enzyme::Chymotrypsin  => EnzymeRules { after: b"FYWL",    before: b"",  universal: None },
            Enzyme::LysC          => EnzymeRules { after: b"K",       before: b"",  universal: None },
            Enzyme::AspN          => EnzymeRules { after: b"",        before: b"D", universal: None },
            Enzyme::GluC          => EnzymeRules { after: b"E",       before: b"",  universal: None },
            Enzyme::LysN          => EnzymeRules { after: b"",        before: b"K", universal: None },
            Enzyme::ArgC          => EnzymeRules { after: b"R",       before: b"",  universal: None },
            Enzyme::AlphaLP       => EnzymeRules { after: b"",        before: b"",  universal: Some(true) },
            Enzyme::NoCleavage    => EnzymeRules { after: b"",        before: b"",  universal: Some(false) },
            Enzyme::NonSpecific   => EnzymeRules { after: b"",        before: b"",  universal: Some(true) },
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Enzyme::Trypsin      => "Trypsin",
            Enzyme::Chymotrypsin => "Chymotrypsin",
            Enzyme::LysC         => "LysC",
            Enzyme::AspN         => "AspN",
            Enzyme::GluC         => "GluC",
            Enzyme::LysN         => "LysN",
            Enzyme::ArgC         => "ArgC",
            Enzyme::AlphaLP      => "aLP",
            Enzyme::NoCleavage   => "NoCleavage",
            Enzyme::NonSpecific  => "NonSpecific",
        }
    }

    /// Case-insensitive name lookup. Aliases mirror Java's
    /// `Enzyme.getEnzymeByName()` (e.g. "Tryp"→Trypsin, "Asp-N"→AspN).
    pub fn from_name(s: &str) -> Option<Self> {
        let n = s.trim().to_ascii_lowercase();
        match n.as_str() {
            "trypsin" | "tryp"        => Some(Enzyme::Trypsin),
            "chymotrypsin" | "chymo"  => Some(Enzyme::Chymotrypsin),
            "lysc" | "lys-c"          => Some(Enzyme::LysC),
            "aspn" | "asp-n"          => Some(Enzyme::AspN),
            "gluc" | "glu-c"          => Some(Enzyme::GluC),
            "lysn" | "lys-n"          => Some(Enzyme::LysN),
            "argc" | "arg-c"          => Some(Enzyme::ArgC),
            "alp" | "alpha-lp" | "alphalp" => Some(Enzyme::AlphaLP),
            "nocleavage" | "none"     => Some(Enzyme::NoCleavage),
            "nonspecific" | "all"     => Some(Enzyme::NonSpecific),
            _                         => None,
        }
    }

    pub fn is_cleavable_after(self, residue: u8) -> bool {
        match self.rules().universal {
            Some(b) => b,
            None    => self.rules().after.contains(&residue),
        }
    }

    pub fn is_cleavable_before(self, residue: u8) -> bool {
        match self.rules().universal {
            Some(b) => b,
            None    => self.rules().before.contains(&residue),
        }
    }

    /// Required by Phase 4's candidate-generation walk. For builtin
    /// enzymes this is always `true`: any residue is allowed *inside* a
    /// peptide. The hook exists for future custom-enzyme support that
    /// might forbid certain residues internally.
    pub fn allows_internal(self, _residue: u8) -> bool {
        true
    }

    // -----------------------------------------------------------------------
    // Phase 6 / Task 5 GF helpers — mirroring Java Enzyme.isNTerm(),
    // isCTerm(), isCleavable(char), and getResidues().
    // -----------------------------------------------------------------------

    /// Returns `true` for N-terminal enzymes (cleavage before the target
    /// residue: LysN, AspN). `false` for C-terminal enzymes (Trypsin, LysC,
    /// ArgC, Chymotrypsin, GluC) and for AlphaLP / NoCleavage /
    /// NonSpecific.
    ///
    /// Java: `Enzyme.isNTerm()` — the flag is set at construction time and
    /// hard-coded per variant. LysN and AspN are the only two builtins
    /// with `isNTerm = true`.
    pub fn is_n_term(self) -> bool {
        matches!(self, Enzyme::LysN | Enzyme::AspN)
    }

    /// `true` for C-terminal enzymes. Mirrors Java `Enzyme.isCTerm() = !isNTerm`.
    pub fn is_c_term(self) -> bool {
        !self.is_n_term()
    }

    /// Direction-agnostic cleavability: returns `true` if `residue` is a
    /// cleavage-target for this enzyme.
    ///
    /// For C-terminal enzymes (`after` list) this is equivalent to
    /// `is_cleavable_after`. For N-terminal enzymes (`before` list) this is
    /// equivalent to `is_cleavable_before`. For NoCleavage always `false`; for
    /// AlphaLP / NonSpecific always `true`. Mirrors Java `Enzyme.isCleavable(char)`.
    pub fn is_cleavable(self, residue: u8) -> bool {
        match self.rules().universal {
            Some(b) => b,
            None => {
                if self.is_n_term() {
                    self.rules().before.contains(&residue)
                } else {
                    self.rules().after.contains(&residue)
                }
            }
        }
    }

    /// The residues targeted by this enzyme's primary cleavage rule.
    ///
    /// For C-terminal enzymes: the `after` residues (e.g. `[b'K', b'R']` for
    /// Trypsin). For N-terminal enzymes: the `before` residues (e.g. `[b'K']`
    /// for LysN). For NoCleavage / NonSpecific / AlphaLP: `&[]` (the
    /// `universal` flag handles cleavability; there are no specific residues).
    ///
    /// Java: `Enzyme.getResidues()` returns a `char[]` that is `null` for
    /// universal enzymes and the target residues otherwise. We return `&[]`
    /// in place of `null`.
    pub fn residues(self) -> &'static [u8] {
        if self.rules().universal.is_some() {
            return &[];
        }
        if self.is_n_term() {
            self.rules().before
        } else {
            self.rules().after
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trypsin_cleaves_after_k_and_r() {
        assert!(Enzyme::Trypsin.is_cleavable_after(b'K'));
        assert!(Enzyme::Trypsin.is_cleavable_after(b'R'));
        assert!(!Enzyme::Trypsin.is_cleavable_after(b'A'));
        assert!(!Enzyme::Trypsin.is_cleavable_before(b'K'));
    }

    #[test]
    fn aspn_cleaves_before_d() {
        assert!(Enzyme::AspN.is_cleavable_before(b'D'));
        assert!(!Enzyme::AspN.is_cleavable_after(b'D'));
        assert!(!Enzyme::AspN.is_cleavable_before(b'A'));
    }

    #[test]
    fn lysc_cleaves_after_k_only() {
        assert!(Enzyme::LysC.is_cleavable_after(b'K'));
        assert!(!Enzyme::LysC.is_cleavable_after(b'R'));
    }

    #[test]
    fn lysn_cleaves_before_k() {
        assert!(Enzyme::LysN.is_cleavable_before(b'K'));
        assert!(!Enzyme::LysN.is_cleavable_after(b'K'));
    }

    #[test]
    fn gluc_cleaves_after_e() {
        assert!(Enzyme::GluC.is_cleavable_after(b'E'));
        assert!(!Enzyme::GluC.is_cleavable_after(b'D'));
    }

    #[test]
    fn no_cleavage_never_cleaves() {
        for r in b'A'..=b'Z' {
            assert!(!Enzyme::NoCleavage.is_cleavable_after(r));
            assert!(!Enzyme::NoCleavage.is_cleavable_before(r));
        }
    }

    #[test]
    fn nonspecific_always_cleaves() {
        for r in b'A'..=b'Z' {
            assert!(Enzyme::NonSpecific.is_cleavable_after(r));
            assert!(Enzyme::NonSpecific.is_cleavable_before(r));
        }
    }

    #[test]
    fn from_name_aliases() {
        assert_eq!(Enzyme::from_name("Trypsin"), Some(Enzyme::Trypsin));
        assert_eq!(Enzyme::from_name("trypsin"), Some(Enzyme::Trypsin));
        assert_eq!(Enzyme::from_name("Tryp"),    Some(Enzyme::Trypsin));
        assert_eq!(Enzyme::from_name("Asp-N"),   Some(Enzyme::AspN));
        assert_eq!(Enzyme::from_name("AspN"),    Some(Enzyme::AspN));
        assert_eq!(Enzyme::from_name("garbage"), None);
    }

    #[test]
    fn argc_cleaves_after_r() {
        assert!(Enzyme::ArgC.is_cleavable_after(b'R'));
        assert!(!Enzyme::ArgC.is_cleavable_after(b'K'));
        assert!(!Enzyme::ArgC.is_cleavable_before(b'R'));
    }

    #[test]
    fn alphalp_is_universal() {
        for r in b'A'..=b'Z' {
            assert!(Enzyme::AlphaLP.is_cleavable_after(r));
            assert!(Enzyme::AlphaLP.is_cleavable_before(r));
        }
    }

    #[test]
    fn from_name_argc_and_alphalp() {
        assert_eq!(Enzyme::from_name("ArgC"), Some(Enzyme::ArgC));
        assert_eq!(Enzyme::from_name("Arg-C"), Some(Enzyme::ArgC));
        assert_eq!(Enzyme::from_name("aLP"), Some(Enzyme::AlphaLP));
        assert_eq!(Enzyme::from_name("AlphaLP"), Some(Enzyme::AlphaLP));
    }

    // Phase 6 / Task 5a: GF helper tests
    #[test]
    fn trypsin_is_c_term_and_cleaves_after_kr() {
        assert!(!Enzyme::Trypsin.is_n_term());
        assert!(Enzyme::Trypsin.is_c_term());
        assert!(Enzyme::Trypsin.is_cleavable(b'K'));
        assert!(Enzyme::Trypsin.is_cleavable(b'R'));
        assert!(!Enzyme::Trypsin.is_cleavable(b'A'));
        let res = Enzyme::Trypsin.residues();
        assert!(res.contains(&b'K'));
        assert!(res.contains(&b'R'));
    }

    #[test]
    fn lysc_is_c_term_and_cleaves_after_k_only() {
        assert!(!Enzyme::LysC.is_n_term());
        assert!(Enzyme::LysC.is_c_term());
        assert!(Enzyme::LysC.is_cleavable(b'K'));
        assert!(!Enzyme::LysC.is_cleavable(b'R'));
        assert_eq!(Enzyme::LysC.residues(), b"K");
    }

    #[test]
    fn nocleavage_residues_is_empty() {
        assert_eq!(Enzyme::NoCleavage.residues(), &[] as &[u8]);
        // NoCleavage.isCleavable should return false for all residues.
        assert!(!Enzyme::NoCleavage.is_cleavable(b'K'));
        assert!(!Enzyme::NoCleavage.is_cleavable(b'R'));
        assert!(!Enzyme::NoCleavage.is_cleavable(b'A'));
    }

    #[test]
    fn lysn_is_n_term_cleaves_before_k() {
        assert!(Enzyme::LysN.is_n_term());
        assert!(!Enzyme::LysN.is_c_term());
        assert!(Enzyme::LysN.is_cleavable(b'K'));
        assert!(!Enzyme::LysN.is_cleavable(b'R'));
        assert_eq!(Enzyme::LysN.residues(), b"K");
    }

    #[test]
    fn aspn_is_n_term_cleaves_before_d() {
        assert!(Enzyme::AspN.is_n_term());
        assert!(!Enzyme::AspN.is_c_term());
        assert!(Enzyme::AspN.is_cleavable(b'D'));
        assert!(!Enzyme::AspN.is_cleavable(b'K'));
        assert_eq!(Enzyme::AspN.residues(), b"D");
    }

    #[test]
    fn nonspecific_residues_is_empty_but_always_cleavable() {
        assert_eq!(Enzyme::NonSpecific.residues(), &[] as &[u8]);
        assert!(Enzyme::NonSpecific.is_cleavable(b'K'));
        assert!(Enzyme::NonSpecific.is_cleavable(b'A'));
    }

    #[test]
    fn name_round_trips() {
        for e in [
            Enzyme::Trypsin, Enzyme::Chymotrypsin, Enzyme::LysC,
            Enzyme::AspN, Enzyme::GluC, Enzyme::LysN,
            Enzyme::ArgC, Enzyme::AlphaLP,
            Enzyme::NoCleavage, Enzyme::NonSpecific,
        ] {
            let n = e.name();
            assert_eq!(Enzyme::from_name(n), Some(e), "round-trip failed for {n}");
        }
    }
}
