//! Guard: the glyco PIN must not regrow its structurally dead columns.
//!
//! Measured 2026-08-29 on pooled human plasma: 19 of 91 columns were
//! structurally constant or byte-identical duplicates (RawScoreCal == RawScore
//! under the 1-hit collapse), and the small-sample Percolator SVM fed on them
//! produced arbitrary-signed weights. Removing them plus the per-scan
//! spectrum-level columns and retraining at --trainFDR 0.05 took the pooled
//! yield from 256.8 +/- 16.5 to 384.6 +/- 23 glycoPSMs @1% with entrapment FDP
//! 0.00% on all five seeds. This test pins the column policy so the dead set
//! cannot silently return -- the repo's recurring defect class.
use std::io::Cursor;

#[test]
fn glyco_pin_default_header_carries_no_dead_columns() {
    let mut buf = Cursor::new(Vec::new());
    output::glyco_pin::write_glyco_header_for_test(&mut buf, 2, 4, false).unwrap();
    let hdr = String::from_utf8(buf.into_inner()).unwrap();
    let cols: Vec<&str> = hdr.trim_end().split('\t').collect();
    for dead in [
        "mass", "IsolationWindowEfficiency", "PrecursorIsotopeKL", "PrecursorSNR",
        "IsRefinement", "NumMods", "RefinementModClass",
        "ModSiteShiftedMatched", "ModSiteShiftedFrac", "ModSiteIntensFrac",
        "ModSiteLocalized", "ModSiteDetCount", "RawScoreCal", "DeltaRTRank",
    ] {
        assert!(
            !cols.contains(&dead),
            "structurally dead column {dead} is back in the default glyco PIN"
        );
    }
    // and the live ones the policy must NOT have taken with it
    for live in ["RawScore", "RankScore", "YHitFrac", "CzHyperscore", "IsTransferred"] {
        assert!(cols.contains(&live), "live column {live} missing from default header");
    }
}

#[test]
fn glyco_pin_curated_header_is_exactly_the_validated_set() {
    let mut buf = Cursor::new(Vec::new());
    output::glyco_pin::write_glyco_header_for_test(&mut buf, 2, 4, true).unwrap();
    let hdr = String::from_utf8(buf.into_inner()).unwrap();
    let cols: Vec<&str> = hdr.trim_end().split('\t').collect();
    // The validated set is defined by MEMBERSHIP, not a count: charge one-hot
    // columns vary with the run's charge range (the plasma validation run
    // carried 52 columns at its range). Every emitted column must be either a
    // charge one-hot or in the validated keep-list, and no keep-listed
    // non-charge column may be missing.
    let keep: &[&str] = &[
        "SpecId", "Label", "ScanNr", "ExpMass", "CalcMass",
        "RankScore", "RankScoreFloat", "RawScore", "TailorScore", "EdgeScore",
        "NumMatchedMainIons", "matchedIonRatio", "longest_y_pct",
        "ExplainedIonCurrentRatio", "NTermIonCurrentRatio", "CTermIonCurrentRatio",
        "ComplementaryIonBalance", "MeanMatchedIntensityRank", "PpmGaussianScore",
        "ChanceMatchSurprise", "MassCompetitionEvidence", "RichIonLLR",
        "FragPredExplained", "FragPredChanceLLR", "IntensitySignal",
        "dm", "absdm", "peplen", "isotope_error",
        "enzN", "enzC", "enzInt",
        "DeltaRT", "AbsDeltaRT", "DeltaRTNorm", "IsobaricRTMargin",
        "OxoniumScore", "NCoreOxoniumIons", "YLadderScore", "YHitFrac", "CoreYHits",
        "PartialGlycanBY", "Y0Y1Anchor", "SialicConsistency", "GlycanMass",
        "Peptide", "Proteins",
    ];
    for c in &cols {
        assert!(
            c.starts_with("charge") || keep.contains(c),
            "unexpected column {c} in curated glyco PIN"
        );
    }
    for k in keep {
        assert!(cols.contains(k), "validated column {k} missing from curated PIN");
    }
    for gone in ["CzHyperscore", "IsTransferred", "MS2IonCurrent",
                 "CandidateRankEntropy", "ListwiseScoreGap", "DeltaRankScore"] {
        assert!(!cols.contains(&gone), "{gone} must not appear in curated mode");
    }
}

/// The four ID-neutral research emitters (`--glyco-y-tree`, `--glyco-oxonium-llr`,
/// `--glyco-rank-masked`, `--glyco-chance-llr-masked`) were measured
/// identification-neutral and removed. Their columns must not reappear in
/// either header, curated or default.
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
    for curated in [false, true] {
        let mut buf = Cursor::new(Vec::new());
        output::glyco_pin::write_glyco_header_for_test(&mut buf, 2, 4, curated).unwrap();
        let hdr = String::from_utf8(buf.into_inner()).unwrap();
        let cols: Vec<&str> = hdr.trim_end().split('\t').collect();
        for c in retired {
            assert!(!cols.contains(&c), "retired column {c} is back (curated={curated})");
        }
    }
}
