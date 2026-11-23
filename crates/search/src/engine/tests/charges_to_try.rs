use super::charges_to_try;
use crate::params::SearchParams;
use model::aa_set::AminoAcidSetBuilder;
use model::Spectrum;

fn params() -> SearchParams {
    let aa = AminoAcidSetBuilder::new_standard().build().unwrap();
    SearchParams::default_tryptic(aa)
}

fn spec(charge: Option<i32>) -> Spectrum {
    Spectrum {
        precursor_charge: charge,
        ..Default::default()
    }
}

#[test]
fn reported_charge_is_searched_alone() {
    let p = params();
    assert_eq!(charges_to_try(&spec(Some(5)), &p).as_slice(), &[5]);
}

#[test]
fn out_of_range_charge_falls_back_to_the_range() {
    let p = params();
    assert_eq!(charges_to_try(&spec(Some(300)), &p).as_slice(), &[2, 3]);
    assert_eq!(charges_to_try(&spec(Some(-2)), &p).as_slice(), &[2, 3]);
    assert_eq!(charges_to_try(&spec(Some(255)), &p).as_slice(), &[255]);
}

#[test]
fn missing_charge_still_falls_back_to_the_range() {
    let p = params();
    assert_eq!(charges_to_try(&spec(None), &p).as_slice(), &[2, 3]);
}
