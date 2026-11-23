use std::io::Cursor;

#[test]
fn glyco_pin_default_header_carries_no_dead_columns() {
    let mut buf = Cursor::new(Vec::new());
    output::glycopeptides::pin::write_glyco_header_for_test(&mut buf, 2, 4).unwrap();
    let hdr = String::from_utf8(buf.into_inner()).unwrap();
    let cols: Vec<&str> = hdr.trim_end().split('\t').collect();
    for dead in [
        "mass",
        "IsolationWindowEfficiency",
        "PrecursorIsotopeKL",
        "PrecursorSNR",
        "IsRefinement",
        "NumMods",
        "RefinementModClass",
        "ModSiteShiftedMatched",
        "ModSiteShiftedFrac",
        "ModSiteIntensFrac",
        "ModSiteLocalized",
        "ModSiteDetCount",
        "RawScoreCal",
        "DeltaRTRank",
    ] {
        assert!(
            !cols.contains(&dead),
            "structurally dead column {dead} is back in the default glyco PIN"
        );
    }
    for live in ["RawScore", "RankScore", "YHitFrac", "CzHyperscore"] {
        assert!(
            cols.contains(&live),
            "live column {live} missing from default header"
        );
    }
}

#[test]
fn retired_research_columns_are_gone() {
    let retired = [
        "YTreeLLR",
        "YTreeHitFrac",
        "YTreeHighPriorMissing",
        "YTreeDecoyGap",
        "OxoniumCompLLR",
        "RankScoreMasked",
        "MaskedPeakCount",
        "ChanceLlrMasked",
        "ExplainedMasked",
    ];
    let mut buf = Cursor::new(Vec::new());
    output::glycopeptides::pin::write_glyco_header_for_test(&mut buf, 2, 4).unwrap();
    let hdr = String::from_utf8(buf.into_inner()).unwrap();
    let cols: Vec<&str> = hdr.trim_end().split('\t').collect();
    for c in retired {
        assert!(!cols.contains(&c), "retired column {c} is back");
    }
}
