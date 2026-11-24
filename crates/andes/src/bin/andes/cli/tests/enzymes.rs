//! Tests of enzyme flags: names, lists and specificity.

use crate::models::parse_enzymes;
use model::enzyme::Enzyme;

#[test]
fn single_enzyme_has_no_extras() {
    let (primary, extras) = parse_enzymes("trypsin").unwrap();
    assert_eq!(primary, Enzyme::Trypsin);
    assert!(
        extras.is_empty(),
        "single enzyme ⇒ empty extras (bit-identical path)"
    );
}

#[test]
fn proteases_are_a_comma_list() {
    let (p, e) = parse_enzymes("gluc,trypsin").unwrap();
    assert_eq!(
        (p, e.as_slice()),
        (Enzyme::Trypsin, [Enzyme::GluC].as_slice())
    );
    let err = parse_enzymes("gluc+trypsin").unwrap_err().to_string();
    assert!(err.contains("separate proteases with a comma"), "{err}");
}

#[test]
fn the_listed_order_of_proteases_does_not_change_the_parsed_set() {
    for (a, b) in [
        ("trypsin,lysc", "lysc,trypsin"),
        ("trypsin,gluc", "gluc,trypsin"),
        ("trypsin,argc", "argc,trypsin"),
        ("lysn,trypsin", "trypsin,lysn"),
        ("nonspecific,trypsin", "trypsin,nonspecific"),
        ("lysc,argc,trypsin", "argc,trypsin,lysc"),
        ("trypsin,chymotrypsin", "chymotrypsin,trypsin"),
    ] {
        assert_eq!(
            parse_enzymes(a).unwrap(),
            parse_enzymes(b).unwrap(),
            "{a} vs {b}"
        );
    }
}

#[test]
fn the_primary_protease_sets_most_peptide_c_termini() {
    assert_eq!(parse_enzymes("lysc,trypsin").unwrap().0, Enzyme::Trypsin);
    assert_eq!(parse_enzymes("argc,trypsin").unwrap().0, Enzyme::Trypsin);
    assert_eq!(parse_enzymes("lysn,trypsin").unwrap().0, Enzyme::Trypsin);
    assert_eq!(parse_enzymes("aspn,lysc").unwrap().0, Enzyme::LysC);
    assert_eq!(parse_enzymes("nonspecific,lysc").unwrap().0, Enzyme::LysC);
    assert_eq!(parse_enzymes("gluc,lysc").unwrap().0, Enzyme::LysC);
    assert_eq!(
        parse_enzymes("argc,chymotrypsin").unwrap().0,
        Enzyme::Chymotrypsin
    );
}

#[test]
fn trypsin_is_the_primary_of_any_set_that_contains_it() {
    for spec in [
        "trypsin,chymotrypsin",
        "chymotrypsin,trypsin",
        "gluc,trypsin",
        "trypsin,gluc",
        "lysc,trypsin",
        "alphalp,trypsin",
        "nonspecific,trypsin",
        "aspn,chymotrypsin,trypsin",
    ] {
        assert_eq!(parse_enzymes(spec).unwrap().0, Enzyme::Trypsin, "{spec}");
    }
}

#[test]
fn a_single_protease_is_its_own_primary() {
    use model::enzyme::Enzyme::*;
    for e in [
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
    ] {
        assert_eq!(parse_enzymes(e.name()).unwrap(), (e, Vec::new()), "{e:?}");
    }
}

#[test]
fn duplicate_tokens_are_deduped_preserving_order() {
    let (primary, extras) = parse_enzymes("trypsin, tryp , gluc, gluc").unwrap();
    assert_eq!(primary, Enzyme::Trypsin);
    assert_eq!(extras, vec![Enzyme::GluC]);
}

#[test]
fn whitespace_and_empties_tolerated() {
    let (primary, extras) = parse_enzymes("  lysc , , trypsin ").unwrap();
    assert_eq!(primary, Enzyme::Trypsin);
    assert_eq!(extras, vec![Enzyme::LysC]);
}

#[test]
fn unknown_enzyme_errors_and_empty_errors() {
    assert!(parse_enzymes("bogus").is_err());
    assert!(parse_enzymes("").is_err());
    assert!(parse_enzymes("  , ").is_err());
}

#[test]
fn elastase_aliases_to_nonspecific() {
    assert_eq!(parse_enzymes("elastase").unwrap().0, Enzyme::NonSpecific);
}
