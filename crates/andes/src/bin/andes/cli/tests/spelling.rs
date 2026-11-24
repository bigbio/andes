//! Tests of flag spelling: names, aliases and value lists.

use super::flags::{parse_enum, Fragmentation, PredictFeature, PredictorBackend, ScoreFlag};
use super::options::SearchArgs;
use search::decoy::DecoyStrategy;
use search::PrecursorCalMode;

use clap::Parser;

#[derive(Parser, Debug)]
struct Wrap {
    #[command(flatten)]
    args: SearchArgs,
}

fn parse(flag: &str, value: &str) -> Result<SearchArgs, clap::Error> {
    Wrap::try_parse_from(["andes", flag, value]).map(|w| w.args)
}

fn accepts(flag: &str, values: &[&str]) {
    for v in values {
        if let Err(e) = parse(flag, v) {
            panic!("{flag} {v:?} rejected: {e}");
        }
    }
}

fn rejects(flag: &str, values: &[&str]) {
    for v in values {
        assert!(parse(flag, v).is_err(), "{flag} {v:?} accepted");
    }
}

#[test]
fn listed_flags_accept_the_spellings_of_their_parsers() {
    accepts(
        "--decoy-strategy",
        &[
            "shuffle",
            "Shuffle",
            "shuf",
            "SHUF",
            "reverse",
            "rev",
            "REV",
            " reverse ",
            "sequon",
            "SEQUON",
            "none",
            "off",
            "NONE",
        ],
    );
    rejects("--decoy-strategy", &["bogus", "", "random"]);
    assert_eq!(
        parse("--decoy-strategy", "rev").unwrap().decoy_strategy,
        DecoyStrategy::Reverse
    );
    assert_eq!(
        parse("--decoy-strategy", "off").unwrap().decoy_strategy,
        DecoyStrategy::None
    );

    accepts(
        "--enzyme",
        &[
            "trypsin",
            "Trypsin",
            "tryp",
            "chymotrypsin",
            "chymo",
            "lysc",
            "lys-c",
            "aspn",
            "Asp-N",
            "gluc",
            "glu-c",
            "lysn",
            "lys-n",
            "argc",
            "arg-c",
            "alp",
            "alpha-lp",
            "alphalp",
            "nocleavage",
            "none",
            "nonspecific",
            "all",
            "elastase",
            "gluc,trypsin",
            "  lysc , , trypsin ",
        ],
    );
    rejects(
        "--enzyme",
        &["bogus", "", " , ", "trypsin,bogus", "gluc+trypsin"],
    );
    assert_eq!(
        parse("--enzyme", "gluc,trypsin").unwrap().db.enzyme,
        "gluc,trypsin",
        "the value is kept as given"
    );

    accepts(
        "--tmt",
        &[
            "tmt6",
            "TMT10",
            "tmt10plex",
            "TMT11plex",
            "tmt16",
            "tmtpro",
            "tmtpro16",
            "tmt18",
            "tmtpro18",
            "tmt32",
            "tmtpro32",
            "tmt35",
            "tmtpro35",
            "itraq4",
            "iTRAQ8plex",
        ],
    );
    rejects("--tmt", &["tmt9", "", "tmt"]);

    accepts("--uniprot-evidence", &["experimental", "Curated", "ANY"]);
    rejects("--uniprot-evidence", &["strong", ""]);

    accepts("--precursor-cal", &["off", "OFF", "Auto", "on"]);
    rejects("--precursor-cal", &["bogus", "1", ""]);
    assert_eq!(
        parse("--precursor-cal", "On").unwrap().precursor_cal,
        PrecursorCalMode::On
    );
}

#[test]
fn case_insensitive_enum_flags() {
    accepts(
        "--fragmentation",
        &[
            "auto", "AUTO", "CID", "cid", "ETD", "etd", "HCD", "Hcd", "UVPD", "uvpd",
        ],
    );
    rejects("--fragmentation", &["3", "99", "", "EThcD"]);
    assert_eq!(
        parse("--fragmentation", "hcd").unwrap().fragmentation,
        Fragmentation::Hcd
    );
    accepts(
        "--protocol",
        &[
            "auto",
            "phospho",
            "PHOSPHO",
            "iTRAQ",
            "itraq",
            "iTRAQ-phospho",
            "ITRAQ-PHOSPHO",
            "TMT",
            "tmt",
            "standard",
        ],
    );
    rejects("--protocol", &["4", "42", "", "silac"]);
    accepts(
        "--enzyme-specificity",
        &[
            "non-specific",
            "NON-SPECIFIC",
            "semi",
            "Semi",
            "fully",
            "Fully",
        ],
    );
    rejects("--enzyme-specificity", &["2", "nonspecific", ""]);
}

#[test]
fn exact_enum_flags_are_unchanged() {
    accepts("--score", &["auto", "rank", "strong"]);
    rejects("--score", &["Rank", "STRONG"]);
    accepts("--fragment-index", &["auto", "on", "off"]);
    rejects("--fragment-index", &["OFF"]);
    accepts("--tmt-unmix", &["auto", "off", "on"]);
    rejects("--tmt-unmix", &["Auto"]);
    accepts(
        "--predict",
        &[
            "rt",
            "ms2",
            "peptdeep-ms2",
            "peptdeep-rt",
            "deeplc-rt",
            "none",
            "rt,ms2",
        ],
    );
    rejects("--predict", &["RT", "ms3"]);
    accepts("--predictor-backend", &["tract", "ort"]);
    rejects("--predictor-backend", &["ORT"]);
    accepts("--glyco-core", &["n", "o"]);
    rejects("--glyco-core", &["N"]);
    accepts(
        "--glyco-species",
        &["human", "mouse", "mouse-large", "high-mannose"],
    );
    rejects("--glyco-species", &["Human", "rat"]);
    accepts("--precursor-mono", &["auto", "off"]);
    rejects("--precursor-mono", &["on"]);
    accepts("--ethcd-activation", &["ethcd", "etd"]);
    rejects("--ethcd-activation", &["hcd"]);
}

#[test]
fn config_enum_parser_names_the_values() {
    let e = parse_enum::<PredictFeature>("predict entry", "ms3", true).unwrap_err();
    assert_eq!(
        e,
        "invalid predict entry `ms3`: expected rt|ms2|peptdeep-ms2|peptdeep-rt|deeplc-rt|none"
    );
    assert_eq!(
        parse_enum::<PredictorBackend>("predictor_backend", "ORT", true).unwrap(),
        PredictorBackend::Ort
    );
    assert!(parse_enum::<ScoreFlag>("score", "Rank", false).is_err());
}
