//! Tests of the value parsers' domain checks (ranges, tolerances, fractions).

use crate::cli::parse::{parse_non_negative_f32, parse_positive_tol};

use super::*;

#[test]
fn charge_range_rejects_zero_and_out_of_bounds() {
    assert!(parse_charge_range("2..5").is_ok());
    assert!(parse_charge_range("1..1").is_ok());
    assert!(
        parse_charge_range("0..5").is_err(),
        "charge 0 must be rejected"
    );
    assert!(
        parse_charge_range("0..0").is_err(),
        "charge 0..0 must be rejected"
    );
    assert!(
        parse_charge_range("2..60").is_err(),
        "charge above bound must be rejected"
    );
}

#[test]
fn open_search_delta_parses_a_range() {
    assert_eq!(parse_open_search_delta("-100..300"), Ok((-100.0, 300.0)));
    assert_eq!(
        parse_open_search_delta(" -20.5..79.966 "),
        Ok((-20.5, 79.966))
    );
    assert!(parse_open_search_delta("-100-300").is_err());
}

#[test]
fn open_search_delta_rejects_bad_windows() {
    assert!(parse_open_search_delta("300..-100").is_err(), "min >= max");
    assert!(
        parse_open_search_delta("10..300").is_err(),
        "0 not in window"
    );
    assert!(
        parse_open_search_delta("-300..-1").is_err(),
        "0 not in window"
    );
    assert!(parse_open_search_delta("0..0").is_err(), "empty window");
    assert!(parse_open_search_delta("nan..300").is_err(), "non-finite");
    assert!(parse_open_search_delta("-1200..1200").is_err(), "width cap");
    assert!(parse_open_search_delta("abc").is_err(), "no separator");
}

#[test]
fn precursor_tol_rejects_nonpositive_and_nonfinite() {
    assert!(parse_precursor_tol("20ppm").is_ok());
    assert!(parse_precursor_tol("0.02da").is_ok());
    assert!(parse_precursor_tol("0ppm").is_err(), "zero tol rejected");
    assert!(
        parse_precursor_tol("-5ppm").is_err(),
        "negative tol rejected"
    );
    assert!(parse_precursor_tol("infda").is_err(), "inf tol rejected");
    assert!(parse_precursor_tol("nanppm").is_err(), "nan tol rejected");
}

#[test]
fn unit_fraction_rejects_out_of_range() {
    assert_eq!(parse_unit_fraction("0.01").unwrap(), 0.01);
    assert_eq!(parse_unit_fraction("0").unwrap(), 0.0);
    assert_eq!(parse_unit_fraction("1").unwrap(), 1.0);
    assert!(parse_unit_fraction("-0.1").is_err());
    assert!(parse_unit_fraction("1.5").is_err());
    assert!(parse_unit_fraction("nan").is_err());
    assert!(parse_unit_fraction("inf").is_err());
}

#[test]
fn non_negative_f32_rejects_negative_and_nonfinite() {
    assert_eq!(parse_non_negative_f32("0.3").unwrap(), 0.3);
    assert_eq!(parse_non_negative_f32("0").unwrap(), 0.0);
    assert!(parse_non_negative_f32("-0.1").is_err());
    assert!(parse_non_negative_f32("nan").is_err());
    assert!(parse_non_negative_f32("inf").is_err());
    assert!(parse_non_negative_f32("1e300").is_err());
}

#[test]
fn positive_tol_rejects_nonpositive_and_nonfinite() {
    assert_eq!(parse_positive_tol("20").unwrap(), 20.0);
    assert!(parse_positive_tol("0").is_err());
    assert!(parse_positive_tol("-1").is_err());
    assert!(parse_positive_tol("nan").is_err());
    assert!(parse_positive_tol("inf").is_err());
}
