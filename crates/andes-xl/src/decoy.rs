//! Target/decoy classes of a crosslink-spectrum match.
//!
//! Both peptides come from the same target+decoy peptide table and are looked
//! up and scored the same way, so a CSM is target-target (TT), target-decoy
//! (TD, either order) or decoy-decoy (DD).

/// Target/decoy class of a peptide pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum PairClass {
    /// Both peptides are targets.
    TT,
    /// One target and one decoy peptide.
    TD,
    /// Both peptides are decoys.
    DD,
}

impl PairClass {
    /// The class of a pair whose peptides have decoy flags `a` and `b`.
    pub fn of(a_decoy: bool, b_decoy: bool) -> Self {
        match (a_decoy, b_decoy) {
            (false, false) => PairClass::TT,
            (true, true) => PairClass::DD,
            _ => PairClass::TD,
        }
    }

    /// Output label.
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
