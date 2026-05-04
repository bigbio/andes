//! Search protocol categories. Mirrors Java
//! `edu.ucsd.msjava.msutil.Protocol`.

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Protocol {
    Automatic,
    Phosphorylation,
    ITRAQ,
    ITRAQPhospho,
    TMT,
    Standard,
}

impl Protocol {
    pub fn name(self) -> &'static str {
        match self {
            Protocol::Automatic       => "Automatic",
            Protocol::Phosphorylation => "Phosphorylation",
            Protocol::ITRAQ           => "iTRAQ",
            Protocol::ITRAQPhospho    => "iTRAQPhospho",
            Protocol::TMT             => "TMT",
            Protocol::Standard        => "Standard",
        }
    }

    /// Case-sensitive lookup — matches Java's `Protocol.get(name)`.
    pub fn from_name(s: &str) -> Option<Self> {
        match s {
            "Automatic"       => Some(Protocol::Automatic),
            "Phosphorylation" => Some(Protocol::Phosphorylation),
            "iTRAQ"           => Some(Protocol::ITRAQ),
            "iTRAQPhospho"    => Some(Protocol::ITRAQPhospho),
            "TMT"             => Some(Protocol::TMT),
            "Standard"        => Some(Protocol::Standard),
            _                 => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn name_round_trips() {
        for p in [
            Protocol::Automatic, Protocol::Phosphorylation,
            Protocol::ITRAQ,     Protocol::ITRAQPhospho,
            Protocol::TMT,       Protocol::Standard,
        ] {
            assert_eq!(Protocol::from_name(p.name()), Some(p));
        }
    }

    #[test]
    fn from_name_known_variants() {
        assert_eq!(Protocol::from_name("Automatic"),       Some(Protocol::Automatic));
        assert_eq!(Protocol::from_name("Phosphorylation"), Some(Protocol::Phosphorylation));
        assert_eq!(Protocol::from_name("iTRAQ"),           Some(Protocol::ITRAQ));
        assert_eq!(Protocol::from_name("iTRAQPhospho"),    Some(Protocol::ITRAQPhospho));
        assert_eq!(Protocol::from_name("TMT"),             Some(Protocol::TMT));
        assert_eq!(Protocol::from_name("Standard"),        Some(Protocol::Standard));
    }

    #[test]
    fn from_name_case_sensitive() {
        assert_eq!(Protocol::from_name("itraq"), None);
        assert_eq!(Protocol::from_name("automatic"), None);
    }

    #[test]
    fn from_name_unknown() {
        assert_eq!(Protocol::from_name("garbage"), None);
        assert_eq!(Protocol::from_name(""), None);
    }
}
