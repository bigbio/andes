use super::*;
use model::{AminoAcidSetBuilder, Protein, ProteinDb};

#[test]
fn sync_resets_the_identity_map() {
    let db = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"MKWVTFISLLRKDAFLGSFLYEYSRK".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&db, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let params = SearchParams::default_tryptic(aa_set);
    let scorer = crate::engine::coisolation::tests::tiny_scorer();
    let mut prepared =
        PreparedSearch::prepare_in_memory(&idx, &params, &scorer, 0.02, "XXX").unwrap();
    let cands = prepared.candidates_where(|_, _| true);
    assert!(cands.len() >= 2);
    {
        let mut accum = prepared.mmap_accum.lock().unwrap();
        assert_eq!(accum.intern(&cands[0]), 0);
        assert_eq!(accum.intern(&cands[1]), 1);
    }
    prepared.sync_materialized_candidates();
    assert_eq!(prepared.candidates.len(), 2);
    let mut accum = prepared.mmap_accum.lock().unwrap();
    assert_eq!(accum.len(), 0);
    assert_eq!(accum.intern(&cands[1]), 0);
    assert_eq!(accum.len(), 1);
}

#[test]
fn truncate_forgets_only_the_later_candidates() {
    let db = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"MKWVTFISLLRKDAFLGSFLYEYSRK".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&db, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let params = SearchParams::default_tryptic(aa_set);
    let scorer = crate::engine::coisolation::tests::tiny_scorer();
    let prepared = PreparedSearch::prepare_in_memory(&idx, &params, &scorer, 0.02, "XXX").unwrap();
    let cands = prepared.candidates_where(|_, _| true);
    assert!(cands.len() >= 3);
    {
        let mut accum = prepared.mmap_accum.lock().unwrap();
        accum.intern(&cands[0]);
        accum.intern(&cands[1]);
        accum.intern(&cands[2]);
    }
    prepared.truncate_materialized(1);
    assert_eq!(prepared.materialized_len(), 1);
    let mut accum = prepared.mmap_accum.lock().unwrap();
    assert_eq!(accum.len(), 1);
    assert_eq!(accum.intern(&cands[0]), 0);
    assert_eq!(accum.intern(&cands[2]), 1);
    assert_eq!(accum.intern(&cands[1]), 2);
}

#[test]
fn dia_ms1_gate_reads_neighbouring_survey_scans() {
    let db = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"MKWVTFISLLRK".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&db, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let params = SearchParams::default_tryptic(aa_set);
    let scorer = crate::engine::coisolation::tests::tiny_scorer();
    let mut prepared =
        PreparedSearch::prepare_in_memory(&idx, &params, &scorer, 0.02, "XXX").unwrap();
    let scan = |rt: f64, mz: f64| model::scan::Ms1Scan {
        rt,
        peaks: vec![(mz, 1.0)],
        ..Default::default()
    };
    prepared.set_dia_ms1(vec![
        scan(30.0, 700.0),
        scan(10.0, 500.0),
        scan(20.0, 600.0),
    ]);
    assert_eq!(prepared.dia_ms1_scans(None), None);
    let r = prepared.dia_ms1_scans(Some(9.0)).unwrap();
    assert_eq!(r, 0..2);
    assert!(prepared.dia_ms1_has(r.clone(), 600.0 * (1.0 + 4e-6), 5.0));
    assert!(!prepared.dia_ms1_has(r.clone(), 600.0 * (1.0 + 6e-6), 5.0));
    assert!(prepared.dia_ms1_has(r.clone(), 600.0 * (1.0 + 6e-6), 25.0));
    assert!(!prepared.dia_ms1_has(r, 700.0, 25.0));
    assert_eq!(prepared.dia_ms1_scans(Some(26.0)).unwrap(), 1..3);

    let e = prepared
        .dia_ms1_error_ppm(Some(20.0), 600.0 * (1.0 - 7e-6), 25.0)
        .unwrap();
    assert!((e - 7.0).abs() < 1e-3, "{e}");
    assert_eq!(prepared.dia_ms1_error_ppm(Some(20.0), 650.0, 25.0), None);
    prepared.shift_dia_ms1_ppm(7.0);
    let e = prepared
        .dia_ms1_error_ppm(Some(20.0), 600.0 * (1.0 - 7e-6), 25.0)
        .unwrap();
    assert!(e.abs() < 1e-3, "{e}");
}

#[test]
fn windowed_top_k_keeps_the_k_most_intense_per_window_in_mz_order() {
    let peaks = [
        (101.0, 1.0),
        (150.0, 5.0),
        (180.0, 3.0),
        (199.0, 4.0),
        (250.0, 2.0),
    ];
    assert_eq!(
        windowed_top_k(&peaks, 100.0, 2),
        vec![(150.0, 5.0), (199.0, 4.0), (250.0, 2.0)]
    );
    assert_eq!(windowed_top_k(&peaks, 100.0, 9), peaks.to_vec());
}

#[test]
fn compact_store_returns_each_candidate_bit_for_bit() {
    use model::modification::{ModLocation, Modification, ResidueSpec};
    let target = ProteinDb {
        proteins: vec![
            Protein {
                accession: "P1".into(),
                description: "".into(),
                sequence: b"MKWVTFISLLRKDAFLGSFLYEYSRKMPEMLLKHPMAGWSTK".to_vec(),
                annotation: None,
            },
            Protein {
                accession: "P2".into(),
                description: "".into(),
                sequence: b"MDEMAKGLLMRNPQSTVWYKMMCHERK".to_vec(),
                annotation: None,
            },
        ],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let ox = Modification {
        name: "Oxidation".to_string(),
        mass_delta: 15.99491,
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
    };
    let aa_set = AminoAcidSetBuilder::new_standard()
        .add_variable_mod(ox)
        .build()
        .unwrap();
    let mut params = SearchParams::default_tryptic(aa_set);
    params.min_length = 3;
    params.max_variable_mods_per_peptide = 2;
    let scorer = crate::engine::coisolation::tests::tiny_scorer();
    let mut prepared =
        PreparedSearch::prepare_in_memory(&idx, &params, &scorer, 0.02, "XXX").unwrap();
    let mut cands = prepared.candidates_where(|_, _| true);
    assert!(cands
        .iter()
        .any(|c| c.peptide.residues.iter().any(|a| a.mod_.is_some())));
    cands[0].peptide = cands[0].peptide.clone().with_charge(3);
    let same = |a: &Candidate, b: &Candidate| {
        a.peptide == b.peptide
            && a.peptide.mass().to_bits() == b.peptide.mass().to_bits()
            && a.peptide.nominal_mass() == b.peptide.nominal_mass()
            && a.peptide.nominal_residue_mass() == b.peptide.nominal_residue_mass()
            && a.peptide
                .residues
                .iter()
                .zip(&b.peptide.residues)
                .all(|(x, y)| x.mod_.map(|m| m as *const _) == y.mod_.map(|m| m as *const _))
            && a.protein_index == b.protein_index
            && a.start_offset_in_protein == b.start_offset_in_protein
            && a.is_decoy == b.is_decoy
            && a.is_protein_n_term == b.is_protein_n_term
            && a.is_protein_c_term == b.is_protein_c_term
    };
    {
        let mut accum = prepared.mmap_accum.lock().unwrap();
        for (i, c) in cands.iter().enumerate() {
            assert_eq!(accum.intern(c), i as u32);
        }
        for (i, c) in cands.iter().enumerate().rev() {
            assert_eq!(accum.intern(c), i as u32);
            assert_eq!(accum.is_decoy(i), Some(c.is_decoy));
            assert!(same(&accum.get(i).unwrap(), c), "candidate {i}");
            assert_eq!(accum.peptide(i).unwrap().charge, c.peptide.charge);
        }
        assert_eq!(accum.get(cands.len()).map(|_| ()), None);
    }
    prepared.sync_materialized_candidates();
    assert_eq!(prepared.candidates.len(), cands.len());
    for (a, b) in prepared.candidates.iter().zip(&cands) {
        assert!(same(a, b));
    }
}

#[test]
fn compact_store_keeps_candidates_whole_when_the_residue_table_is_full() {
    let db = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"MKWVTFISLLRK".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&db, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let params = SearchParams::default_tryptic(aa_set);
    let scorer = crate::engine::coisolation::tests::tiny_scorer();
    let mut prepared =
        PreparedSearch::prepare_in_memory(&idx, &params, &scorer, 0.02, "XXX").unwrap();
    let cand = |k: u32| {
        let residues: Vec<model::AminoAcid> = (0..3)
            .map(|j| model::AminoAcid {
                residue: b'G',
                mass: 57.0 + f64::from(3 * k + j) * 1e-6,
                mod_: None,
            })
            .collect();
        Candidate {
            peptide: Peptide::new(residues, b'K', b'K'),
            protein_index: k as usize,
            start_offset_in_protein: 0,
            is_decoy: k % 2 == 1,
            is_protein_n_term: false,
            is_protein_c_term: false,
        }
    };
    let cands: Vec<Candidate> = (0..120).map(cand).collect();
    {
        let mut accum = prepared.mmap_accum.lock().unwrap();
        for (i, c) in cands.iter().enumerate() {
            assert_eq!(accum.intern(c), i as u32);
        }
        for (i, c) in cands.iter().enumerate() {
            assert_eq!(accum.intern(c), i as u32);
            let got = accum.get(i).unwrap();
            assert_eq!(got.peptide, c.peptide, "candidate {i}");
            assert_eq!(got.peptide.mass().to_bits(), c.peptide.mass().to_bits());
            assert_eq!(got.protein_index, c.protein_index);
            assert_eq!(got.is_decoy, c.is_decoy);
        }
    }
    prepared.sync_materialized_candidates();
    assert_eq!(prepared.candidates.len(), cands.len());
    for (a, b) in prepared.candidates.iter().zip(&cands) {
        assert_eq!(a.peptide, b.peptide);
        assert_eq!(a.peptide.mass().to_bits(), b.peptide.mass().to_bits());
        assert_eq!(a.protein_index, b.protein_index);
    }
}

#[test]
fn the_per_spectrum_ms1_gate_answers_as_the_per_scan_search() {
    let db = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"MKWVTFISLLRK".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&db, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let params = SearchParams::default_tryptic(aa_set);
    let scorer = crate::engine::coisolation::tests::tiny_scorer();
    let mut prepared =
        PreparedSearch::prepare_in_memory(&idx, &params, &scorer, 0.02, "XXX").unwrap();
    let mut rng = 0x9e37_79b9_7f4a_7c15_u64;
    let mut next = move || {
        rng ^= rng << 13;
        rng ^= rng >> 7;
        rng ^= rng << 17;
        rng
    };
    let mut unif = move || (next() >> 11) as f64 / (1u64 << 53) as f64;
    let mut scans: Vec<model::scan::Ms1Scan> = (0..12)
        .map(|i| {
            let mut rows: Vec<(f64, f32, f32)> = (0..3000)
                .map(|_| (400.0 + unif() * 200.0, 1.0, 0.6 + unif() as f32))
                .collect();
            rows.sort_by(|a, b| a.0.total_cmp(&b.0));
            rows[11].0 = rows[10].0;
            model::scan::Ms1Scan {
                rt: f64::from(i),
                peaks: rows.iter().map(|r| (r.0, r.1)).collect(),
                mobility: if i % 3 == 1 {
                    Vec::new()
                } else {
                    rows.iter().map(|r| r.2).collect()
                },
                ..Default::default()
            }
        })
        .collect();
    scans[9].peaks.swap(5, 2000);
    prepared.set_dia_ms1(scans);
    let tol_ppm = 8.0;
    let mut checked = 0usize;
    for rt in 0..12 {
        let rt = Some(f64::from(rt) + 0.3);
        let scans = prepared.dia_ms1_scans(rt).unwrap();
        for k0 in [None, Some((0.9_f32, 1.1_f32))] {
            for _ in 0..4 {
                let lo_mz = 410.0 + unif() * 160.0;
                let windows: Vec<(u8, f64, f64)> = [2u8, 3]
                    .iter()
                    .map(|&z| {
                        let zf = f64::from(z);
                        (z, (lo_mz - PROTON) * zf, (lo_mz + 12.0 - PROTON) * zf)
                    })
                    .collect();
                let gate = prepared.dia_ms1_gate(rt, &windows, tol_ppm, k0).unwrap();
                for &(z, lo, hi) in &windows {
                    let zf = f64::from(z);
                    let mut masses: Vec<f64> = (0..400).map(|_| lo + unif() * (hi - lo)).collect();
                    masses.extend([lo, hi, lo.next_up(), hi.next_down()]);
                    for s in &prepared.dia_ms1[scans.clone()] {
                        for p in &s.peaks {
                            for mz in [
                                p.0,
                                p.0 / (1.0 + tol_ppm * 1e-6),
                                p.0 / (1.0 - tol_ppm * 1e-6),
                            ] {
                                let m = (mz - PROTON) * zf;
                                masses.extend([m, m.next_up(), m.next_down()]);
                            }
                        }
                    }
                    for mass in masses.into_iter().filter(|m| (lo..=hi).contains(m)) {
                        let mz = mass / zf + PROTON;
                        assert_eq!(
                            prepared.dia_ms1_gate_has(&gate, mz, tol_ppm),
                            prepared.dia_ms1_has_at(scans.clone(), mz, tol_ppm, k0),
                            "rt={rt:?} k0={k0:?} z={z} mz={mz}"
                        );
                        checked += 1;
                    }
                }
            }
        }
    }
    assert!(checked > 100_000, "{checked}");
}
