//! Tests of the generated Methods paragraph and its references.

use super::methods::{methods_paragraph, references, MethodsFacts};

fn facts<'a>(fixed: &'a [String], variable: &'a [String], rescored: bool) -> MethodsFacts<'a> {
    MethodsFacts {
        version: "0.4.0",
        database: "ecoli.fasta",
        enzyme: "Trypsin",
        missed_cleavages: 1,
        precursor_tolerance: "20 ppm",
        fragment_tolerance: "20 ppm",
        fixed_mods: fixed,
        variable_mods: variable,
        decoys: "shuffling each target peptide (seed 42)",
        rescored,
        psm_fdr: rescored.then_some(0.01),
        protein_fdr: rescored.then_some(0.01),
    }
}

#[test]
fn methods_without_rescoring_has_no_percolator_sentence() {
    let (text, keys) = methods_paragraph(&facts(&[], &[], false));
    assert!(!text.contains("Percolator"), "{text}");
    assert_eq!(keys, ["kim2014_msgfplus", "elias2007_target_decoy"]);
}

#[test]
fn methods_numbers_citations_in_order() {
    let fixed = vec!["Carbamidomethyl (C)".to_string()];
    let variable = vec!["Oxidation (M)".to_string()];
    let (text, keys) = methods_paragraph(&facts(&fixed, &variable, true));
    let at = |m: &str| text.find(m).unwrap_or_else(|| panic!("no {m} in {text}"));
    assert!(at("[1]") < at("[2]") && at("[2]") < at("[3]") && at("[3]") < at("[4]"));
    assert!(!text.contains("[cite"), "{text}");
    assert!(text.contains("Percolator [3] and filtered at 0.01 PSM-level FDR."), "{text}");
    assert!(text.contains("Fixed modifications: Carbamidomethyl (C)."), "{text}");
    assert_eq!(keys.len(), 4);
}

#[test]
fn methods_skips_empty_mod_lists() {
    let variable = vec!["Oxidation (M)".to_string()];
    let (text, _) = methods_paragraph(&facts(&[], &variable, false));
    assert!(!text.contains("Fixed modifications:"), "{text}");
    assert!(text.contains("Variable modifications: Oxidation (M)."), "{text}");
}

#[test]
fn every_methods_key_is_in_the_registry() {
    let (_, keys) = methods_paragraph(&facts(&[], &[], true));
    let refs = references(&keys);
    assert_eq!(refs.len(), 4);
    assert!(refs.iter().all(|r| r.text != "unknown citation key"), "{refs:?}");
    assert!(refs[0].text.starts_with("Kim S, Pevzner PA (2014). "), "{:?}", refs[0]);
    assert!(refs[0].text.ends_with("doi:10.1038/ncomms6277"), "{:?}", refs[0]);
    assert_eq!(references(&["no_such_key"])[0].text, "unknown citation key");
}
