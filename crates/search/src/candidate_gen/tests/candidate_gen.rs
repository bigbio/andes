#[test]
fn decoy_membership_requires_underscore_delimiter() {
    use crate::decoy::is_decoy_accession;
    assert!(is_decoy_accession("XXX_protein1", "XXX"));
    assert!(
        !is_decoy_accession("XXXprotein1", "XXX"),
        "no '_' delimiter ⇒ a real 'XXX...' target must NOT be read as decoy"
    );
    assert!(!is_decoy_accession("DECOY_protein1", "XXX"));
    assert!(is_decoy_accession("rev_sp|P1", "rev"));
    assert!(!is_decoy_accession("reverse_kinase", "rev"));
}

#[test]
fn multi_enzyme_cleavage_is_union_of_rules() {
    use super::compute_cleavage_positions;
    use model::enzyme::Enzyme;
    let seq = b"MKAESRPE";

    let trypsin = compute_cleavage_positions(seq, Enzyme::Trypsin, &[]);
    assert_eq!(trypsin, vec![0, 2, 6, 8], "Trypsin alone: after K(2), R(6)");

    let gluc = compute_cleavage_positions(seq, Enzyme::GluC, &[]);
    assert_eq!(gluc, vec![0, 4, 8], "GluC alone: after first E(4)");

    let both = compute_cleavage_positions(seq, Enzyme::GluC, &[Enzyme::Trypsin]);
    assert_eq!(
        both,
        vec![0, 2, 4, 6, 8],
        "GluC+Trypsin = union of cut sites"
    );
    let both_rev = compute_cleavage_positions(seq, Enzyme::Trypsin, &[Enzyme::GluC]);
    assert_eq!(both, both_rev, "union is order-independent");
}

#[test]
fn single_enzyme_path_unchanged_with_empty_extras() {
    use super::compute_cleavage_positions;
    use model::enzyme::Enzyme;
    let seq = b"MKAESRPEDK";
    for e in [
        Enzyme::Trypsin,
        Enzyme::LysC,
        Enzyme::GluC,
        Enzyme::AspN,
        Enzyme::NonSpecific,
        Enzyme::NoCleavage,
    ] {
        let single = compute_cleavage_positions(seq, e, &[]);
        let with_nocleave = compute_cleavage_positions(seq, e, &[Enzyme::NoCleavage]);
        assert_eq!(
            single, with_nocleave,
            "{e:?}: NoCleavage extra adds no sites"
        );
    }
    let ns = compute_cleavage_positions(seq, Enzyme::Trypsin, &[Enzyme::NonSpecific]);
    assert_eq!(ns, (0..=seq.len() as u32).collect::<Vec<_>>());
}

use super::{
    enumerate_candidates, lazy_candidates_for_nominal_window, lazy_candidates_for_precursor,
    Candidate,
};
use crate::candidate_index::{build_base_peptide_index, MmapCandidateIndex};
use crate::search_index::SearchIndex;
use crate::search_params::SearchParams;
use model::aa_set::AminoAcidSetBuilder;
use model::modification::{ModLocation, Modification, ResidueSpec};
use model::protein::{Protein, ProteinDb};

fn canonical_key(c: &Candidate) -> (bool, usize, usize, String, u64) {
    (
        c.is_decoy,
        c.protein_index,
        c.start_offset_in_protein,
        c.peptide.to_string(),
        (c.peptide.mass() * 1000.0).round() as u64,
    )
}

fn base_mass_of(db: &SearchIndex, params: &SearchParams, seq: &[u8]) -> f64 {
    let base_params = SearchParams {
        max_variable_mods_per_peptide: 0,
        ..params.clone()
    };
    let cands: Vec<_> = enumerate_candidates(db, &base_params, "XXX").collect();
    cands
        .into_iter()
        .filter(|c| !c.is_decoy)
        .find(|c| {
            c.peptide
                .residues
                .iter()
                .map(|aa| aa.residue)
                .eq(seq.iter().copied())
        })
        .unwrap_or_else(|| {
            panic!(
                "no unmodified candidate spelling {:?}",
                std::str::from_utf8(seq)
            )
        })
        .peptide
        .mass()
}

fn build_and_open(
    db: &SearchIndex,
    params: &SearchParams,
) -> (tempfile::NamedTempFile, MmapCandidateIndex) {
    use std::io::BufWriter;
    let tmp = tempfile::NamedTempFile::new().expect("tempfile");
    {
        let mut bw = BufWriter::new(tmp.as_file());
        build_base_peptide_index(db, params, "XXX", &mut bw)
            .expect("build_base_peptide_index failed");
    }
    let mi = MmapCandidateIndex::open(tmp.path()).expect("open index");
    (tmp, mi)
}

fn oxidation_m() -> Modification {
    Modification {
        name: "Oxidation".to_string(),
        mass_delta: 15.994915,
        residue: ResidueSpec::Specific(b'M'),
        location: ModLocation::Anywhere,
        fixed: false,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
        annotated_site_mass: None,
        blocks_cleavage: false,
        parts: Vec::new(),
        diagnostic_ions: Vec::new(),
    }
}

fn label_mod(name: &str, delta: f64, residue: u8, location: ModLocation) -> Modification {
    use model::modification::ResidueSpec;
    Modification {
        name: name.into(),
        mass_delta: delta,
        residue: if residue == b'*' {
            ResidueSpec::Wildcard
        } else {
            ResidueSpec::Specific(residue)
        },
        location,
        fixed: false,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
        annotated_site_mass: None,
        blocks_cleavage: false,
        diagnostic_ions: Vec::new(),
        parts: Vec::new(),
    }
}

fn label_names(c: &Candidate) -> Vec<String> {
    c.peptide
        .residues
        .iter()
        .filter_map(|aa| aa.mod_.map(|m| m.name.clone()))
        .collect()
}

#[test]
fn non_specific_digest_labels_every_site() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "silac".into(),
            sequence: b"AKKKKA".to_vec(),
            annotation: None,
        }],
    };
    let db = SearchIndex::from_target_db(&target, "XXX");
    let k8 = label_mod("Label:13C(6)15N(2)", 8.014199, b'K', ModLocation::Anywhere);
    let aa_set = AminoAcidSetBuilder::new_standard()
        .add_variable_mod(k8)
        .with_label_channels(vec![("Label:13C(6)15N(2)".into(), 1)])
        .build()
        .unwrap();
    let mut params = SearchParams::default_tryptic(aa_set);
    params.min_length = 6;
    params.num_tolerable_termini = 0;
    params.max_missed_cleavages = 0;
    params.max_variable_mods_per_peptide = 0;
    params.max_label_sites = 8;
    let full: Vec<_> = enumerate_candidates(&db, &params, "XXX")
        .filter(|c| c.peptide.residues.len() == 6)
        .map(|c| label_names(&c).len())
        .collect();
    assert!(full.contains(&0), "{full:?}");
    assert!(full.contains(&4), "{full:?}");
    assert!(full.iter().all(|&n| n == 0 || n == 4), "{full:?}");
}

#[test]
fn fixed_mod_on_a_labelable_site_is_neutral() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "fixed".into(),
            sequence: b"MAGLKPEPTIDEK".to_vec(),
            annotation: None,
        }],
    };
    let db = SearchIndex::from_target_db(&target, "XXX");
    let mut b = AminoAcidSetBuilder::new_standard();
    let mut channels = Vec::new();
    for (c, (name, d)) in [("Dimethyl", 28.0313), ("Dimethyl:2H(4)", 32.056407)]
        .into_iter()
        .enumerate()
    {
        for (res, loc, suffix) in [
            (b'K', ModLocation::Anywhere, ""),
            (b'*', ModLocation::NTerm, "-N"),
        ] {
            let nm = format!("{name}{suffix}");
            b = b.add_variable_mod(label_mod(&nm, d, res, loc));
            channels.push((nm, c as u8 + 1));
        }
    }
    let mut acetyl = label_mod("Acetyl", 42.010565, b'*', ModLocation::ProtNTerm);
    acetyl.fixed = true;
    let aa_set = b
        .add_fixed_mod(acetyl)
        .with_label_channels(channels)
        .build()
        .unwrap();
    let mut params = SearchParams::default_tryptic(aa_set);
    params.min_length = 3;
    params.max_label_sites = 3;
    let kept: Vec<_> = enumerate_candidates(&db, &params, "XXX").collect();
    let maglk: Vec<Vec<String>> = kept
        .iter()
        .filter(|c| c.peptide.residues.iter().map(|a| a.residue).eq(*b"MAGLK"))
        .map(label_names)
        .collect();
    for ch in ["Dimethyl", "Dimethyl:2H(4)"] {
        assert!(
            maglk.iter().any(|n| n == &["Acetyl", ch]),
            "{ch}: {maglk:?}"
        );
    }
}

#[test]
fn mass_offsets_of_a_large_label_allowance() {
    let mut b = AminoAcidSetBuilder::new_standard();
    let mut channels = Vec::new();
    for (r, d) in [
        (b'K', 8.0),
        (b'R', 9.0),
        (b'L', 10.0),
        (b'I', 11.0),
        (b'V', 12.0),
    ] {
        let name = format!("Label:{}", r as char);
        b = b.add_variable_mod(label_mod(&name, d, r, ModLocation::Anywhere));
        channels.push((name, 1));
    }
    let mut params =
        SearchParams::default_tryptic(b.with_label_channels(channels).build().unwrap());
    params.max_variable_mods_per_peptide = 0;
    params.max_label_sites = 100;
    let mut got: Vec<i64> = crate::candidate_gen::distinct_mod_mass_offsets(&params)
        .iter()
        .map(|&(_, d)| d.round() as i64)
        .collect();
    got.sort_unstable();
    let mut want: Vec<i64> = (0..=100i64).flat_map(|k| 8 * k..=12 * k).collect();
    want.sort_unstable();
    want.dedup();
    assert_eq!(got, want);
    assert!(std::sync::Arc::ptr_eq(
        &crate::candidate_gen::distinct_mod_mass_offsets(&params),
        &crate::candidate_gen::distinct_mod_mass_offsets(&params),
    ));
}

#[test]
fn memoized_offsets_equal_the_plain_walk() {
    fn plain(params: &SearchParams) -> Vec<u64> {
        let deltas = crate::candidate_gen::variable_mod_deltas(params);
        fn rec(
            items: &[crate::candidate_gen::FormDelta],
            start: usize,
            sum: f64,
            budget: (u32, u32),
            out: &mut std::collections::BTreeSet<u64>,
        ) {
            for (i, f) in items.iter().enumerate().skip(start) {
                if f.mods > budget.0 || f.labels > budget.1 {
                    continue;
                }
                let s = sum + f.delta;
                out.insert(crate::candidate_gen::milli_of(s));
                rec(items, i, s, (budget.0 - f.mods, budget.1 - f.labels), out);
            }
        }
        let mut out = std::collections::BTreeSet::from([crate::candidate_gen::milli_of(0.0)]);
        let mut channels: Vec<u8> = deltas.iter().map(|f| f.channel).collect();
        channels.push(0);
        channels.sort_unstable();
        channels.dedup();
        for c in channels {
            let items: Vec<_> = deltas
                .iter()
                .filter(|f| f.channel == 0 || f.channel == c)
                .copied()
                .collect();
            rec(
                &items,
                0,
                0.0,
                (params.max_variable_mods_per_peptide, params.max_label_sites),
                &mut out,
            );
        }
        out.into_iter().collect()
    }
    let memo = |params: &SearchParams| -> Vec<u64> {
        let mut v: Vec<u64> = crate::candidate_gen::distinct_mod_mass_offsets(params)
            .iter()
            .map(|&(k, _)| k)
            .collect();
        v.sort_unstable();
        v.dedup();
        v
    };
    let build = |labels: &[(&str, f64, u8, u8)], ox: bool, mods: u32, sites: u32| {
        let mut b = AminoAcidSetBuilder::new_standard();
        for &(name, d, r, _) in labels {
            b = b.add_variable_mod(label_mod(name, d, r, ModLocation::Anywhere));
        }
        if ox {
            b = b.add_variable_mod(oxidation_m());
        }
        let channels = labels
            .iter()
            .map(|&(n, _, _, c)| (n.to_string(), c))
            .collect();
        let mut p =
            SearchParams::default_tryptic(b.with_label_channels(channels).build().unwrap());
        p.max_variable_mods_per_peptide = mods;
        p.max_label_sites = sites;
        p
    };
    let coderabbit = build(
        &[
            ("A", 1.0, b'K', 1),
            ("B", 1.000_000_24, b'L', 1),
            ("C", 10.000_249_9, b'R', 1),
        ],
        false,
        0,
        4,
    );
    let got = memo(&coderabbit);
    assert_eq!(got, plain(&coderabbit));
    assert!(got.contains(&22001), "{got:?}");
    let silac = build(
        &[
            ("Label:2H(4)", 4.025107, b'K', 1),
            ("Label:13C(6)", 6.020129, b'R', 1),
            ("Label:13C(6)15N(2)", 8.014199, b'K', 2),
            ("Label:13C(6)15N(4)", 10.008269, b'R', 2),
        ],
        true,
        1,
        6,
    );
    assert_eq!(memo(&silac), plain(&silac));
}

#[test]
fn mass_offsets_keep_channels_apart() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "offsets".into(),
            sequence: b"PEMKTKIDERAAKGGLK".to_vec(),
            annotation: None,
        }],
    };
    let db = SearchIndex::from_target_db(&target, "XXX");
    let aa = AminoAcidSetBuilder::new_standard()
        .add_variable_mod(label_mod(
            "Label:2H(4)",
            4.025107,
            b'K',
            ModLocation::Anywhere,
        ))
        .add_variable_mod(label_mod(
            "Label:13C(6)15N(2)",
            8.014199,
            b'K',
            ModLocation::Anywhere,
        ))
        .add_variable_mod(oxidation_m())
        .with_label_channels(vec![
            ("Label:2H(4)".into(), 1),
            ("Label:13C(6)15N(2)".into(), 2),
        ])
        .build()
        .unwrap();
    let mut params = SearchParams::default_tryptic(aa);
    params.min_length = 3;
    params.max_variable_mods_per_peptide = 1;
    params.max_label_sites = 3;
    let offsets: Vec<f64> = crate::candidate_gen::distinct_mod_mass_offsets(&params)
        .iter()
        .map(|&(_, d)| d)
        .collect();
    let has = |d: f64| offsets.iter().any(|o| (o - d).abs() < 1e-6);
    assert!(!has(4.025107 + 8.014199), "{offsets:?}");
    assert!(has(3.0 * 8.014199 + 15.994915));
    assert!(!has(2.0 * 15.994915));
    for c in enumerate_candidates(&db, &params, "XXX") {
        let d: f64 = c
            .peptide
            .residues
            .iter()
            .filter_map(|a| a.mod_.map(|m| m.mass_delta))
            .sum();
        assert!(has(d), "{d} missing from {offsets:?}");
    }
}

#[test]
fn silac_labels_are_channel_consistent() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "silac".into(),
            sequence: b"PEPKTIDERAAGGLK".to_vec(),
            annotation: None,
        }],
    };
    let db = SearchIndex::from_target_db(&target, "XXX");
    let k8 = label_mod("Label:13C(6)15N(2)", 8.014199, b'K', ModLocation::Anywhere);
    let r10 = label_mod("Label:13C(6)15N(4)", 10.008269, b'R', ModLocation::Anywhere);
    let build = |labels: bool| {
        let mut b = AminoAcidSetBuilder::new_standard()
            .add_variable_mod(k8.clone())
            .add_variable_mod(r10.clone())
            .add_variable_mod(oxidation_m());
        if labels {
            b = b.with_label_channels(vec![
                ("Label:13C(6)15N(2)".into(), 1),
                ("Label:13C(6)15N(4)".into(), 1),
            ]);
        }
        let mut params = SearchParams::default_tryptic(b.build().unwrap());
        params.min_length = 3;
        params.max_variable_mods_per_peptide = 3;
        params.max_label_sites = 3;
        params
    };
    let free = build(false);
    let labelled = build(true);
    let all: Vec<_> = enumerate_candidates(&db, &free, "XXX").collect();
    let kept: Vec<_> = enumerate_candidates(&db, &labelled, "XXX").collect();
    let span = |c: &Candidate| {
        c.peptide
            .residues
            .iter()
            .map(|a| a.residue as char)
            .collect::<String>()
    };
    let n = |v: &[Candidate], s: &str| v.iter().filter(|c| span(c) == s).count();
    assert_eq!(n(&all, "PEPKTIDER"), 4);
    assert_eq!(n(&kept, "PEPKTIDER"), 2);
    for c in kept.iter().filter(|c| span(c) == "PEPKTIDER") {
        let names = label_names(c);
        assert!(names.is_empty() || names.len() == 2, "{names:?}");
    }
    assert_eq!(n(&kept, "TIDER"), 2);
    let mut no_mods = build(true);
    no_mods.max_variable_mods_per_peptide = 0;
    let strict: Vec<_> = enumerate_candidates(&db, &no_mods, "XXX").collect();
    assert_eq!(n(&strict, "PEPKTIDER"), 2);
    assert!(strict.iter().all(|c| c
        .peptide
        .residues
        .iter()
        .all(|a| a.mod_.is_none_or(|m| !m.name.contains("Oxidation")))));
    let tables = crate::candidate_gen::FormWalkTables::new(&labelled);
    let mut walked = 0;
    crate::candidate_gen::walk_forms_bounded(
        b"PEPKTIDER",
        &labelled,
        &tables,
        true,
        false,
        crate::candidate_gen::SpanSites::new(&target.proteins[0], 0),
        0.0,
        f64::MAX,
        |_, _, _| walked += 1,
    );
    assert_eq!(walked, 2);
}

#[test]
fn max_mods_bounds_only_non_label_mods() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "silac".into(),
            sequence: b"PEPMKTIDERAAGGLK".to_vec(),
            annotation: None,
        }],
    };
    let db = SearchIndex::from_target_db(&target, "XXX");
    let k8 = label_mod("Label:13C(6)15N(2)", 8.014199, b'K', ModLocation::Anywhere);
    let r10 = label_mod("Label:13C(6)15N(4)", 10.008269, b'R', ModLocation::Anywhere);
    let aa_set = AminoAcidSetBuilder::new_standard()
        .add_variable_mod(k8)
        .add_variable_mod(r10)
        .add_variable_mod(oxidation_m())
        .with_label_channels(vec![
            ("Label:13C(6)15N(2)".into(), 1),
            ("Label:13C(6)15N(4)".into(), 1),
        ])
        .build()
        .unwrap();
    let oxidised = |c: &[&model::amino_acid::AminoAcid]| {
        c.iter()
            .any(|a| a.mod_.is_some_and(|m| m.name.contains("Oxidation")))
    };
    for (max_mods, expect, with_ox) in [(1u32, 4usize, 2usize), (0, 2, 0)] {
        let mut params = SearchParams::default_tryptic(aa_set.clone());
        params.min_length = 3;
        params.max_variable_mods_per_peptide = max_mods;
        params.max_label_sites = 3;
        let forms: Vec<Vec<model::amino_acid::AminoAcid>> =
            enumerate_candidates(&db, &params, "XXX")
                .filter(|c| {
                    c.peptide
                        .residues
                        .iter()
                        .map(|a| a.residue)
                        .collect::<Vec<_>>()
                        == b"PEPMKTIDER"
                })
                .map(|c| c.peptide.residues)
                .collect();
        assert_eq!(forms.len(), expect, "max_mods {max_mods}");
        let refs: Vec<Vec<&model::amino_acid::AminoAcid>> =
            forms.iter().map(|f| f.iter().collect()).collect();
        assert_eq!(refs.iter().filter(|f| oxidised(f)).count(), with_ox);
        let tables = crate::candidate_gen::FormWalkTables::new(&params);
        let (mut walked, mut walked_ox) = (0, 0);
        crate::candidate_gen::walk_forms_bounded(
            b"PEPMKTIDER",
            &params,
            &tables,
            true,
            false,
            crate::candidate_gen::SpanSites::new(&target.proteins[0], 0),
            0.0,
            f64::MAX,
            |_, chosen, _| {
                walked += 1;
                walked_ox += usize::from(oxidised(chosen));
            },
        );
        assert_eq!(
            (walked, walked_ox),
            (expect, with_ox),
            "max_mods {max_mods}"
        );
    }
}

#[test]
fn dimethyl_labels_cover_n_terminus_and_lysine() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "dimethyl".into(),
            sequence: b"MAGLKPEPTIDEK".to_vec(),
            annotation: None,
        }],
    };
    let db = SearchIndex::from_target_db(&target, "XXX");
    let mut b = AminoAcidSetBuilder::new_standard();
    let mut channels = Vec::new();
    for (c, (name, d)) in [
        ("Dimethyl", 28.0313),
        ("Dimethyl:2H(4)", 32.056407),
        ("Dimethyl:2H(6)13C(2)", 36.07567),
    ]
    .into_iter()
    .enumerate()
    {
        for (res, loc, suffix) in [
            (b'K', ModLocation::Anywhere, ""),
            (b'*', ModLocation::NTerm, "-N"),
        ] {
            let nm = format!("{name}{suffix}");
            b = b.add_variable_mod(label_mod(&nm, d, res, loc));
            channels.push((nm, c as u8 + 1));
        }
    }
    let acetyl = label_mod("Acetyl", 42.010565, b'*', ModLocation::ProtNTerm);
    let aa_set = b
        .add_variable_mod(acetyl)
        .with_label_channels(channels)
        .build()
        .unwrap();
    let mut params = SearchParams::default_tryptic(aa_set);
    params.min_length = 3;
    params.max_variable_mods_per_peptide = 1;
    params.max_label_sites = 3;
    let kept: Vec<_> = enumerate_candidates(&db, &params, "XXX").collect();
    let of = |s: &str| -> Vec<Vec<String>> {
        kept.iter()
            .filter(|c| {
                c.peptide
                    .residues
                    .iter()
                    .map(|a| a.residue as char)
                    .collect::<String>()
                    == s
            })
            .map(label_names)
            .collect()
    };
    let pk = of("PEPTIDEK");
    assert_eq!(pk.len(), 4, "{pk:?}");
    for names in &pk {
        assert!(names.is_empty() || names.len() == 2, "{names:?}");
    }
    let mk = of("MAGLK");
    assert_eq!(mk.len(), 8, "{mk:?}");
}

fn deamidation_n() -> Modification {
    Modification {
        name: "Deamidation".to_string(),
        mass_delta: 0.984016,
        residue: ResidueSpec::Specific(b'N'),
        location: ModLocation::Anywhere,
        fixed: false,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
        annotated_site_mass: None,
        blocks_cleavage: false,
        parts: Vec::new(),
        diagnostic_ions: Vec::new(),
    }
}

fn oxm_fixture() -> (SearchIndex, SearchParams) {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "ox-m fixture".into(),
            sequence: b"PEPTMKDALMRQPK".to_vec(),
            annotation: None,
        }],
    };
    let db = SearchIndex::from_target_db(&target, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard()
        .add_variable_mod(oxidation_m())
        .build()
        .unwrap();
    let mut params = SearchParams::default_tryptic(aa_set);
    params.min_length = 3;
    params.max_variable_mods_per_peptide = 1;
    (db, params)
}

#[test]
fn fixed_n_term_tag_applies_to_protein_n_term_peptides() {
    use model::modification::ResidueSpec;
    let tmt = Modification {
        name: "TMT6plex".into(),
        mass_delta: 229.162932,
        residue: ResidueSpec::Wildcard,
        location: ModLocation::NTerm,
        fixed: true,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
        annotated_site_mass: None,
        blocks_cleavage: false,
        parts: Vec::new(),
        diagnostic_ions: Vec::new(),
    };
    let acetyl = Modification {
        name: "Acetyl".into(),
        mass_delta: 42.010565,
        residue: ResidueSpec::Wildcard,
        location: ModLocation::ProtNTerm,
        fixed: false,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
        annotated_site_mass: None,
        blocks_cleavage: false,
        parts: Vec::new(),
        diagnostic_ions: Vec::new(),
    };
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "tmt protein-nterm".into(),
            sequence: b"AAAPEPTIDEKQQQR".to_vec(),
            annotation: None,
        }],
    };
    let db = SearchIndex::from_target_db(&target, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard()
        .add_fixed_mod(tmt)
        .add_variable_mod(acetyl)
        .build()
        .unwrap();
    let mut params = SearchParams::default_tryptic(aa_set);
    params.min_length = 3;
    params.max_variable_mods_per_peptide = 1;

    let pnt: Vec<_> = enumerate_candidates(&db, &params, "XXX")
        .filter(|c| !c.is_decoy && c.start_offset_in_protein == 0)
        .collect();
    assert!(
        !pnt.is_empty(),
        "expected a protein-N-term peptide candidate"
    );
    let has_tmt = pnt.iter().any(|c| {
        c.peptide
            .residues
            .first()
            .and_then(|aa| aa.mod_.as_ref())
            .map(|m| (m.mass_delta - 229.162932).abs() < 1e-3)
            .unwrap_or(false)
    });
    assert!(
        has_tmt,
        "protein-N-term peptide must carry the fixed N-term TMT tag (finding 2.3)"
    );
}

#[test]
fn a_mod_that_blocks_cleavage_never_sits_on_the_cleaved_residue() {
    let seq = b"AAAAKGGGGKCCCCK".to_vec();
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: String::new(),
            sequence: seq.clone(),
            annotation: None,
        }],
    };
    let db = SearchIndex::from_target_db(&target, "XXX");
    let run = |blocks: bool| {
        let acetyl = Modification {
            name: "Acetyl".into(),
            mass_delta: 42.010565,
            residue: ResidueSpec::Specific(b'K'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: blocks,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        let aa = AminoAcidSetBuilder::new_standard()
            .add_variable_mod(acetyl)
            .build()
            .unwrap();
        let mut params = SearchParams::default_tryptic(aa);
        params.min_length = 3;
        params.max_variable_mods_per_peptide = 2;
        enumerate_candidates(&db, &params, "XXX")
            .filter(|c| !c.is_decoy)
            .map(|c| {
                let mods: Vec<usize> = c
                    .peptide
                    .residues
                    .iter()
                    .enumerate()
                    .filter(|(_, a)| a.mod_.is_some())
                    .map(|(i, _)| i)
                    .collect();
                (c.start_offset_in_protein, c.peptide.length(), mods)
            })
            .collect::<Vec<_>>()
    };
    let blocked = run(true);
    for (start, len, mods) in &blocked {
        let protein_end = start + len == seq.len();
        assert!(
            protein_end || !mods.contains(&(len - 1)),
            "acetyl on the cleaved K of {start}+{len}"
        );
    }
    assert!(
        blocked
            .iter()
            .any(|(s, l, m)| *s == 0 && *l == 10 && m == &vec![4]),
        "acetyl on the internal K of a missed cleavage"
    );
    assert!(
        blocked
            .iter()
            .any(|(s, l, m)| *s == 10 && *l == 5 && m == &vec![4]),
        "acetyl on the protein's own C-terminal K"
    );
    assert!(run(false)
        .iter()
        .any(|(s, l, m)| *s == 0 && *l == 5 && m == &vec![4]));
}

#[test]
fn substitution_at_a_peptide_terminus_uses_the_terminal_forms() {
    use model::protein::{ProteinAnnotation, Substitution, VariantRegion};
    let annotation = ProteinAnnotation {
        variant: vec![VariantRegion { start: 4, end: 5 }],
        ..Default::default()
    }
    .with_substitutions(vec![Substitution {
        position: 4,
        original: b'N',
    }]);
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "variant fixture".into(),
            sequence: b"PEPKDAGTRK".to_vec(),
            annotation: Some(std::sync::Arc::new(annotation)),
        }],
    };
    let db = SearchIndex::from_target_db(&target, "XXX");
    let has_dagtr = |deamidation: bool| {
        let mut b = AminoAcidSetBuilder::new_standard().add_fixed_mod(Modification {
            name: "TMT6plex".into(),
            mass_delta: 229.162932,
            residue: ResidueSpec::Wildcard,
            location: ModLocation::NTerm,
            fixed: true,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        });
        if deamidation {
            b = b.add_variable_mod(deamidation_n());
        }
        let mut params = SearchParams::default_tryptic(b.build().unwrap());
        params.min_length = 3;
        params.max_variable_mods_per_peptide = 1;
        let found = enumerate_candidates(&db, &params, "XXX")
            .any(|c| !c.is_decoy && c.start_offset_in_protein == 4 && c.peptide.length() == 5);
        found
    };
    assert!(
        has_dagtr(false),
        "without deamidation D is evidence for the variant"
    );
    assert!(
        !has_dagtr(true),
        "D+TMT is the canonical N deamidated + TMT"
    );
}

#[test]
fn mature_chain_starting_with_met_keeps_its_second_residue_spans() {
    use model::protein::{ChainEnds, ProteinAnnotation};
    let annotation = ProteinAnnotation {
        chain: Some(ChainEnds {
            n_term: false,
            c_term: true,
        }),
        ..Default::default()
    };
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "chain fixture".into(),
            sequence: b"MAKAGERLLK".to_vec(),
            annotation: Some(std::sync::Arc::new(annotation)),
        }],
    };
    let db = SearchIndex::from_target_db(&target, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let mut params = SearchParams::default_tryptic(aa_set);
    params.min_length = 2;
    params.num_tolerable_termini = 0;
    let starting_at_1: Vec<_> = enumerate_candidates(&db, &params, "XXX")
        .filter(|c| !c.is_decoy && c.start_offset_in_protein == 1)
        .collect();
    assert_eq!(starting_at_1.len(), 1, "AKAGERLLK once");
    assert_eq!(starting_at_1[0].peptide.residues.len(), 9);
    assert!(!starting_at_1[0].is_protein_n_term);
}

#[test]
fn met_cleaved_walk_adds_only_post_met_n_term_peptides() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "met fixture".into(),
            sequence: b"MAKAGERLLKPEPTIDER".to_vec(),
            annotation: None,
        }],
    };
    let db = SearchIndex::from_target_db(&target, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let mut params = SearchParams::default_tryptic(aa_set);
    params.min_length = 2;
    params.max_missed_cleavages = 1;
    for ntt in [2u8, 1, 0] {
        params.num_tolerable_termini = ntt;
        let cands: Vec<_> = enumerate_candidates(&db, &params, "XXX")
            .filter(|c| !c.is_decoy)
            .collect();
        let mut keys: Vec<(usize, usize)> = cands
            .iter()
            .map(|c| (c.start_offset_in_protein, c.peptide.residues.len()))
            .collect();
        let n = keys.len();
        keys.sort_unstable();
        keys.dedup();
        assert_eq!(keys.len(), n, "ntt={ntt}: a span is enumerated twice");
        assert!(
            cands
                .iter()
                .any(|c| c.start_offset_in_protein == 1 && c.is_protein_n_term),
            "ntt={ntt}: missing the Met-cleaved N-terminal peptide"
        );
        assert!(
            cands
                .iter()
                .all(|c| c.is_protein_n_term == (c.start_offset_in_protein <= 1)),
            "ntt={ntt}: only spans starting at 0 or 1 are protein N-terminal"
        );
    }
}

#[test]
fn stacked_terminal_forms_enumerate_and_match_lazy() {
    use model::mod_spec::parse_mod_list;
    let aa_set = AminoAcidSetBuilder::new_standard()
        .add_mods(
            parse_mod_list(
                "TMT6plex (K), TMT6plex (N-term), Carbamidomethyl (C)",
                true,
            )
            .unwrap(),
        )
        .add_mods(parse_mod_list("Oxidation (M), Acetyl (Protein N-term)", false).unwrap())
        .build()
        .unwrap();
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "tmt stacking fixture".into(),
            sequence: b"MCKAGERKLLRMAEPK".to_vec(),
            annotation: None,
        }],
    };
    let db = SearchIndex::from_target_db(&target, "XXX");
    let mut params = SearchParams::default_tryptic(aa_set);
    params.min_length = 2;
    params.max_variable_mods_per_peptide = 2;
    let all: Vec<_> = enumerate_candidates(&db, &params, "XXX").collect();
    let first_delta = |c: &Candidate| {
        c.peptide.residues[0]
            .mod_
            .as_ref()
            .map_or(0.0, |m| m.mass_delta)
    };
    let has = |start: usize, len: usize, delta: f64| {
        all.iter().any(|c| {
            !c.is_decoy
                && c.start_offset_in_protein == start
                && c.peptide.residues.len() == len
                && (first_delta(c) - delta).abs() < 1e-3
        })
    };
    assert!(has(0, 3, 15.994915 + 42.010565));
    assert!(has(1, 2, 57.021464 + 229.162932));
    assert!(has(1, 2, 57.021464 + 42.010565));
    assert!(!has(1, 2, 42.010565));
    assert!(has(7, 4, 2.0 * 229.162932));
    assert!(has(11, 5, 229.162932 + 15.994915));

    let (_tmp, mi) = build_and_open(&db, &params);
    let mut masses: Vec<f64> = all.iter().map(|c| c.peptide.mass()).collect();
    masses.sort_by(|a, b| a.partial_cmp(b).unwrap());
    masses.dedup_by(|a, b| (*a - *b).abs() < 1e-6);
    for &m in &masses {
        assert_lazy_matches_inram(&db, &params, &mi, m, 0.01);
    }
}

fn assert_lazy_matches_inram(
    db: &SearchIndex,
    params: &SearchParams,
    mi: &MmapCandidateIndex,
    precursor: f64,
    tol: f64,
) {
    use std::collections::HashSet;
    let lazy: HashSet<_> = lazy_candidates_for_precursor(mi, db, params, precursor, tol, tol)
        .iter()
        .map(canonical_key)
        .collect();
    let inram: HashSet<_> = enumerate_candidates(db, params, "XXX")
        .filter(|c| (c.peptide.mass() - precursor).abs() <= tol)
        .map(|c| canonical_key(&c))
        .collect();
    assert_eq!(
        lazy,
        inram,
        "lazy enumeration must match the in-RAM candidate set\n\
         only-in-lazy: {:?}\n only-in-inram: {:?}",
        lazy.difference(&inram).collect::<Vec<_>>(),
        inram.difference(&lazy).collect::<Vec<_>>(),
    );
}

fn ram_window_keys(
    db: &SearchIndex,
    params: &SearchParams,
    spec: &model::Spectrum,
    z: u8,
) -> std::collections::HashSet<(bool, usize, usize, String, u64)> {
    let (min_nominal, max_nominal) =
        crate::match_engine::candidate_nominal_bounds(spec, z, params, 0.0);
    enumerate_candidates(db, params, "XXX")
        .filter(|c| {
            let k = c.peptide.nominal_residue_mass();
            k >= min_nominal && k <= max_nominal
        })
        .map(|c| canonical_key(&c))
        .collect()
}

fn spectrum_at_mass(neutral_mass: f64, z: u8) -> model::Spectrum {
    let mz = (neutral_mass + z as f64 * model::mass::PROTON) / z as f64;
    model::Spectrum {
        title: "boundary".into(),
        precursor_mz: mz,
        precursor_intensity: None,
        precursor_charge: Some(z as i32),
        rt_seconds: None,
        scan: None,
        peaks: Vec::new(),
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
        isolation_target_mz: None,
        mobility_window: None,
    }
}

#[test]
fn mmap_nominal_window_equals_ram_bucket_window_on_boundary() {
    use model::tolerance::{PrecursorTolerance, Tolerance};
    use std::collections::HashSet;

    let (db, mut params) = oxm_fixture();
    params.isotope_error_range = 0..=0;
    params.precursor_tolerance = PrecursorTolerance::symmetric(Tolerance::Da(0.5));
    let (_tmp, mi) = build_and_open(&db, &params);

    let z = 1u8;

    let masses: Vec<f64> = enumerate_candidates(&db, &params, "XXX")
        .map(|c| c.peptide.mass())
        .collect();
    assert!(!masses.is_empty(), "fixture produced no candidates");

    let mut boundary_seen = false;
    for &m in &masses {
        for step in -8i32..=8 {
            let precursor = m + (step as f64) * 0.125;
            let spec = spectrum_at_mass(precursor, z);

            let ram = ram_window_keys(&db, &params, &spec, z);

            let (min_nominal, max_nominal) =
                crate::match_engine::candidate_nominal_bounds(&spec, z, &params, 0.0);
            let mmap_nominal: HashSet<_> =
                lazy_candidates_for_nominal_window(&mi, &db, &params, min_nominal, max_nominal)
                    .iter()
                    .map(canonical_key)
                    .collect();
            assert_eq!(
                mmap_nominal,
                ram,
                "NEW mmap nominal window must equal RAM bucket window \
                 (precursor={precursor})\n only-in-mmap: {:?}\n only-in-ram: {:?}",
                mmap_nominal.difference(&ram).collect::<Vec<_>>(),
                ram.difference(&mmap_nominal).collect::<Vec<_>>(),
            );

            let neutral_obs = precursor;
            let tol_left = params.precursor_tolerance.left.as_da(neutral_obs);
            let tol_right = params.precursor_tolerance.right.as_da(neutral_obs);
            let old_f64: HashSet<_> = lazy_candidates_for_precursor(
                &mi,
                &db,
                &params,
                neutral_obs,
                tol_left,
                tol_right,
            )
            .iter()
            .map(canonical_key)
            .collect();
            if old_f64 != ram {
                boundary_seen = true;
            }
        }
    }
    assert!(
        boundary_seen,
        "fixture never exercised a bucket boundary where the OLD f64 window \
         diverged from RAM — the regression guard would not have caught the bug"
    );
}

#[test]
fn lazy_enum_matches_inram_candidate_set_for_window() {
    let (db, params) = oxm_fixture();
    let (_tmp, mi) = build_and_open(&db, &params);

    let base = base_mass_of(&db, &params, b"PEPTMK");
    let p_oxm = base + 15.994915;
    assert_lazy_matches_inram(&db, &params, &mi, p_oxm, 0.01);

    let lazy = lazy_candidates_for_precursor(&mi, &db, &params, p_oxm, 0.01, 0.01);
    assert!(
        lazy.iter().any(|c| {
            !c.is_decoy
                && c.peptide
                    .residues
                    .iter()
                    .map(|aa| aa.residue)
                    .eq(b"PEPTMK".iter().copied())
                && c.peptide.residues.iter().any(|aa| aa.is_modified())
        }),
        "expected the Ox-M PEPTMK peptidoform"
    );

    assert_lazy_matches_inram(&db, &params, &mi, base, 0.01);
    let lazy_base = lazy_candidates_for_precursor(&mi, &db, &params, base, 0.01, 0.01);
    assert!(
        lazy_base.iter().any(|c| {
            !c.is_decoy
                && c.peptide
                    .residues
                    .iter()
                    .map(|aa| aa.residue)
                    .eq(b"PEPTMK".iter().copied())
                && c.peptide.residues.iter().all(|aa| !aa.is_modified())
        }),
        "expected the unmodified PEPTMK peptidoform"
    );

    assert_lazy_matches_inram(&db, &params, &mi, base, 200.0);
}

#[test]
fn lazy_enum_matches_inram_two_variable_mods() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P2".into(),
            description: "two-mod fixture".into(),
            sequence: b"MNPEMTKNMDALERQNMPK".to_vec(),
            annotation: None,
        }],
    };
    let db = SearchIndex::from_target_db(&target, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard()
        .add_variable_mod(oxidation_m())
        .add_variable_mod(deamidation_n())
        .build()
        .unwrap();
    let mut params = SearchParams::default_tryptic(aa_set);
    params.min_length = 3;
    params.max_variable_mods_per_peptide = 2;

    let (_tmp, mi) = build_and_open(&db, &params);

    let all: Vec<_> = enumerate_candidates(&db, &params, "XXX").collect();
    assert!(!all.is_empty(), "fixture must enumerate candidates");
    let mut masses: Vec<f64> = all.iter().map(|c| c.peptide.mass()).collect();
    masses.sort_by(|a, b| a.partial_cmp(b).unwrap());
    masses.dedup_by(|a, b| (*a - *b).abs() < 1e-6);
    for &m in &masses {
        assert_lazy_matches_inram(&db, &params, &mi, m, 0.01);
    }

    let max_mass = masses.last().copied().unwrap();
    let min_mass = masses.first().copied().unwrap();
    let mid = (min_mass + max_mass) / 2.0;
    let half = (max_mass - min_mass) / 2.0 + 1.0;
    assert_lazy_matches_inram(&db, &params, &mi, mid, half);
}
#[test]
fn lazy_enum_matches_inram_at_tolerance_edge() {
    let (db, params) = oxm_fixture();
    let (_tmp, mi) = build_and_open(&db, &params);

    let base = base_mass_of(&db, &params, b"PEPTMK");

    let frac_tol = 0.40_f64;
    let tol = (10.0 + frac_tol) / 1000.0;

    let ox_delta = 15.994915_f64;
    let center = base - tol * 0.999;
    let precursor = center + ox_delta;

    {
        let tol_milli = (tol * 1000.0).round() as i64;
        let center_milli = (center * 1000.0).round() as i64;
        let base_key = (base * 1000.0).round() as i64;
        let hi = center_milli + tol_milli;
        assert!(
            base_key > hi,
            "expected hi-boundary mismatch (base_key={base_key} > hi={hi})                  — rounding conditions may not hold for this mass; verify frac_tol"
        );
    }

    let inram_count: usize = enumerate_candidates(&db, &params, "XXX")
        .filter(|c| (c.peptide.mass() - precursor).abs() <= tol)
        .count();
    assert!(
        inram_count > 0,
        "fixture must yield at least one in-RAM candidate at the engineered precursor              (base={base:.6}, tol={tol:.6}, center={center:.6}, precursor={precursor:.6})"
    );

    assert_lazy_matches_inram(&db, &params, &mi, precursor, tol);
}
