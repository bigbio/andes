#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum PairClass {
    TT,
    TD,
    DD,
}

impl PairClass {
    pub fn of(a_decoy: bool, b_decoy: bool) -> Self {
        match (a_decoy, b_decoy) {
            (false, false) => PairClass::TT,
            (true, true) => PairClass::DD,
            _ => PairClass::TD,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            PairClass::TT => "TT",
            PairClass::TD => "TD",
            PairClass::DD => "DD",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn class_is_symmetric() {
        assert_eq!(PairClass::of(false, false), PairClass::TT);
        assert_eq!(PairClass::of(true, false), PairClass::TD);
        assert_eq!(PairClass::of(false, true), PairClass::TD);
        assert_eq!(PairClass::of(true, true), PairClass::DD);
    }
}
