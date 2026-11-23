use super::*;
use model::isotope::max_isolated_isotope;

pub(crate) fn linker(spec: &str) -> Crosslinker {
    xl::spec::XlDecl::parse(spec)
        .and_then(|d| d.resolve(model::unimod::bundled()))
        .unwrap()
}

fn pep(seq: &[u8]) -> Peptide {
    let residues: Vec<AminoAcid> = seq
        .iter()
        .map(|&r| AminoAcid::standard(r).unwrap())
        .collect();
    Peptide::new(residues, b'K', b'A')
}

#[test]
fn a_peptide_of_two_proteins_lists_both_whatever_its_flanks() {
    use model::aa_set::AminoAcidSetBuilder;
    use model::protein::{Protein, ProteinDb};
    let prot = |acc: &str, seq: &[u8]| Protein {
        accession: acc.into(),
        description: String::new(),
        sequence: seq.to_vec(),
        annotation: None,
    };
    let db = ProteinDb {
        proteins: vec![
            prot("CETN2", b"MSEKAFKLFDDDETGKISFKNLK"),
            prot("CETN1", b"MSEKAFRLFDDDETGKISFKNLK"),
        ],
    };
    let idx = crate::index::SearchIndex::from_target_db(&db, "XXX_");
    let aa = AminoAcidSetBuilder::new_standard().build().unwrap();
    let params = crate::params::SearchParams::default_tryptic(aa);
    let cands: Vec<Candidate> =
        crate::candidates::enumerate_candidates(&idx, &params, "XXX_").collect();
    let table = XlCandidates::from_candidates(cands, &linker("DSSO"));
    let shared = |t: &XlCandidates| {
        t.peptides()
            .iter()
            .find(|p| p.bare() == b"LFDDDETGKISFK" && !p.is_decoy)
            .map(|p| p.proteins.clone())
            .unwrap()
    };
    assert_eq!(shared(&table).len(), 2);
    let scorer = crate::engine::coisolation::tests::tiny_scorer();
    let prepared =
        PreparedSearch::prepare_in_memory(&idx, &params, &scorer, 0.02, "XXX_").unwrap();
    let table = XlCandidates::from_index(&prepared, &linker("DSSO"));
    assert_eq!(shared(&table).len(), 2, "{:?}", shared(&table));
}

#[test]
fn only_the_doublets_in_excess_of_chance_are_paired() {
    let mass = |m: f64, i: f64| DoubletMass {
        mass: m,
        intensity: i,
        n_charges: 1,
        isotope_confirmed: false,
        min_charge: 2,
    };
    let summary = |n: usize, chance: f32| {
        let mut d = DoubletSummary {
            masses: (0..n)
                .map(|k| mass(500.0 + k as f64, (n - k) as f64))
                .collect(),
            ..Default::default()
        };
        d.probe[0][0] = [n as f32, chance];
        d
    };
    let clean = summary(30, 0.0);
    assert_eq!(paired_masses(&clean).len(), 30);
    let noisy = summary(540, 500.5);
    let kept = paired_masses(&noisy);
    assert_eq!(kept.len(), 40);
    assert!(kept
        .iter()
        .all(|m| m.intensity >= noisy.masses[39].intensity));
    assert_eq!(
        paired_masses(&summary(300, 299.0)).len(),
        MIN_PAIRED_DOUBLET_MASSES
    );
    assert_eq!(paired_masses(&summary(5, 9.0)).len(), 5);
}

#[test]
fn an_ms3_precursor_matches_its_stub_peptide_at_its_charge() {
    let dsso = linker("DSSO");
    let alkene = dsso.cleavage.as_ref().unwrap().stubs[0].mass;
    let stub_peptide = 1000.5 + alkene;
    let mz2 = stub_peptide / 2.0 + PROTON;
    assert!(ms3_precursor_matches(
        mz2,
        2,
        stub_peptide,
        max_isolated_isotope(stub_peptide)
    ));
    assert!(ms3_precursor_matches(
        mz2 + 0.009,
        2,
        stub_peptide,
        max_isolated_isotope(stub_peptide)
    ));
    assert!(!ms3_precursor_matches(
        mz2 + 0.02,
        2,
        stub_peptide,
        max_isolated_isotope(stub_peptide)
    ));
    assert!(!ms3_precursor_matches(
        mz2,
        3,
        stub_peptide,
        max_isolated_isotope(stub_peptide)
    ));
    assert!(!ms3_precursor_matches(
        mz2,
        0,
        stub_peptide,
        max_isolated_isotope(stub_peptide)
    ));
    let heavy = 2496.2496 + alkene;
    let m1 = (heavy + ISOTOPE_SPACING) / 2.0 + PROTON;
    assert!(
        ms3_precursor_matches(m1, 2, heavy, max_isolated_isotope(heavy)),
        "M+1 of a 2.5 kDa stub peptide"
    );
    assert!(!ms3_precursor_matches(
        m1 + 0.02,
        2,
        heavy,
        max_isolated_isotope(heavy)
    ));
    let light = 600.0 + alkene;
    let far = (light + 3.0 * ISOTOPE_SPACING) / 2.0 + PROTON;
    assert!(
        !ms3_precursor_matches(far, 2, light, max_isolated_isotope(light)),
        "M+3 of a small stub peptide"
    );
    let s = ms3_spectrum(&ProductScan {
        id: "scan=9".into(),
        precursor_mz: mz2,
        peaks: vec![(300.0, 1.0), (200.0, 2.0)],
        ..Default::default()
    });
    assert_eq!(s.peaks, vec![(200.0, 2.0), (300.0, 1.0)]);
    assert_eq!(s.precursor_charge, None);
}

#[test]
fn by_ions_and_peak_removal() {
    let p = pep(b"GAK");
    let ions = by_ions(&p, 1);
    let g = AminoAcid::standard(b'G').unwrap().mass;
    let k = AminoAcid::standard(b'K').unwrap().mass;
    assert!(ions.iter().any(|&m| (m - (g + PROTON)).abs() < 1e-9));
    assert!(ions
        .iter()
        .any(|&m| (m - (k + model::mass::H2O + PROTON)).abs() < 1e-9));
    assert_eq!(ions.len(), 4);
    let spec = Spectrum {
        title: "t".into(),
        precursor_mz: 0.0,
        precursor_intensity: None,
        precursor_charge: Some(3),
        rt_seconds: None,
        scan: None,
        peaks: vec![(g + PROTON + 0.001, 1.0), (500.0, 2.0)],
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
        isolation_target_mz: None,
        mobility_window: None,
    };
    let kept = without_peaks(&spec, &ions, 0.01);
    assert_eq!(kept.peaks, vec![(500.0, 2.0)]);
}

#[test]
fn the_csm_score_adds_a_share_of_the_gap_to_the_next_pair() {
    assert_eq!(combined_pair_score(40.0, None), 40.0);
    assert!(
        (combined_pair_score(40.0, Some(30.0)) - (40.0 + PAIR_DELTA_WEIGHT * 10.0)).abs()
            < 1e-6
    );
    assert_eq!(combined_pair_score(40.0, Some(45.0)), 40.0);
    assert!(combined_pair_score(40.0, Some(10.0)) > combined_pair_score(40.0, Some(39.0)));
}

fn trial(z: u8, pair: &str, score: f32, second: Option<(&str, f32)>) -> XlSpectrumResult {
    let key = |p: &str| [p.as_bytes().to_vec(), b"PEPTIDEK".to_vec()];
    let hit = PeptideHit {
        pep: 0,
        site: xl::sites::LinkSite {
            pos: 0,
            n_term: false,
            c_term: false,
        },
        form: xl::crosslinker::LinkForm::Intact,
        score,
        form_scores: Default::default(),
        n_sites: 1,
        gap_next: None,
        mass_err_ppm: 0.0,
        ms3: None,
    };
    let mut pairs = vec![(key(pair), score)];
    pairs.extend(second.map(|(p, s)| (key(p), s)));
    XlSpectrumResult {
        spec_idx: 0,
        charge: z,
        precursor_mass: 0.0,
        doublets: Default::default(),
        csm: Some(Csm {
            alpha: hit,
            beta: hit,
            pair_key: key(pair),
            pair_score: score,
            second_pair_score: second.map(|(_, s)| s),
            kind: PairKind::Full,
            iso: 0,
            class: PairClass::TT,
            intra: true,
            homotypic: false,
            precursor_err_ppm: 0.0,
            b_stub_ions: 0,
            n_mass_pairs: 1,
            alpha_doublet: None,
            beta_doublet: None,
            rescores: Default::default(),
            charge: z,
        }),
        oracle: None,
        linker: 0,
        pairs,
    }
}

#[test]
fn an_uncharged_spectrum_keeps_the_runner_up_of_its_chosen_charge() {
    // Charge 3 (100 / 50, combined 115) beats charge 4 (110 / 109, combined 110.3).
    let trials = || {
        [
            trial(3, "AK", 100.0, Some(("BK", 50.0))),
            trial(4, "CK", 110.0, Some(("DK", 109.0))),
        ]
        .into_iter()
    };
    // Several linkers: the charge of the best pair score (110 at charge 4).
    let r = choose_charge(trials(), true).unwrap();
    assert_eq!((r.charge, r.csm.as_ref().unwrap().pair_score), (4, 110.0));
    let r = choose_charge(trials(), false).unwrap();
    let c = r.csm.as_ref().unwrap();
    assert_eq!(
        (r.charge, c.pair_score, c.second_pair_score),
        (3, 100.0, Some(50.0))
    );
    // Every charge's pairs stay for a merge over several linkers.
    assert_eq!(r.pairs.len(), 4);
    assert_eq!(runner_up(&r.pairs, &c.pair_key), Some(110.0));
}

#[test]
fn a_half_pair_of_a_triply_charged_precursor_pays_a_penalty() {
    assert_eq!(
        half_pair_penalty(PairKind::Half, 3),
        LOW_CHARGE_HALF_PAIR_PENALTY
    );
    assert_eq!(half_pair_penalty(PairKind::Half, 4), 0.0);
    assert_eq!(half_pair_penalty(PairKind::Full, 3), 0.0);
    assert_eq!(half_pair_penalty(PairKind::Open, 3), 0.0);
}

#[test]
fn link_mass_is_on_the_site_and_in_the_peptide_mass() {
    let x = Crosslinker::dsso();
    let mods = StubMods::new(&x);
    let p = pep(b"PEPKTIDE");
    let stub = LinkForm::Stub(0).delta(&x, 0.0).unwrap();
    let q = with_link_mass(&p, 3, LinkForm::Stub(0), stub, &mods);
    assert_eq!(q.residues[3].mod_.unwrap().name, "XL:alkene");
    assert!((q.mass() - p.mass() - stub).abs() < 1e-9);
    assert!((q.residues[3].mod_.unwrap().mass_delta - stub).abs() < 1e-12);
    let r = with_link_mass(&p, 3, LinkForm::Stub(0), stub, &mods);
    assert!(std::ptr::eq(
        q.residues[3].mod_.unwrap(),
        r.residues[3].mod_.unwrap()
    ));
    let intact = LinkForm::Intact.delta(&x, 1000.0).unwrap();
    let s = with_link_mass(&p, 3, LinkForm::Intact, intact, &mods);
    assert!((s.mass() - p.mass() - intact).abs() < 1e-9);
    assert!(s.residues[3].mod_.is_none());
    assert_eq!(
        s.nominal_residue_mass(),
        model::mass::nominal_from(s.mass() - model::mass::H2O)
    );
}

#[test]
fn a_terminal_site_is_reported_only_where_the_peptide_is_terminal() {
    let cand = |protein: usize, start: usize, n_term: bool| Candidate {
        peptide: pep(b"AKR"),
        protein_index: protein,
        start_offset_in_protein: start,
        is_decoy: false,
        is_protein_n_term: n_term,
        is_protein_c_term: false,
    };
    let t = XlCandidates::from_candidates(
        vec![cand(0, 0, true), cand(1, 7, false)],
        &Crosslinker::dsso(),
    );
    assert_eq!(t.peptides().len(), 1);
    let p = &t.peptides()[0];
    let n = p.sites.iter().find(|s| s.n_term).unwrap();
    let k = p.sites.iter().find(|s| !s.n_term).unwrap();
    assert_eq!(p.occurrences_at(n), &[(0, 0)]);
    assert_eq!(p.occurrences_at(k), &[(0, 0), (1, 7)]);
}

#[test]
fn defaults_search_two_missing_isotopes_and_credit_stub_ions() {
    let cfg = XlConfig::new(Crosslinker::dsso());
    assert_eq!(cfg.isotope_offsets, 0..=2);
    assert_eq!(cfg.stub_weight, STUB_ION_WEIGHT);
    assert_eq!(cfg.secondary_site_penalty, SECONDARY_SITE_PENALTY);
}

#[test]
fn s_t_y_sites_of_an_nhs_ester_pay_the_secondary_site_penalty() {
    let x = crate::crosslinks::tests::linker("DSSO (K, S, T, Y, Protein N-term)");
    let p = pep(b"SAKTY");
    let at = |pos: usize, n_term: bool| LinkSite {
        pos,
        n_term,
        c_term: false,
    };
    assert_eq!(site_penalty(12.0, &p, &at(2, false), &x, 0), 0.0);
    assert_eq!(site_penalty(12.0, &p, &at(0, false), &x, 0), 12.0);
    assert_eq!(site_penalty(12.0, &p, &at(3, false), &x, 0), 12.0);
    assert_eq!(site_penalty(12.0, &p, &at(0, true), &x, 0), 0.0);
    assert_eq!(site_penalty(0.0, &p, &at(0, false), &x, 0), 0.0);
    let sda = crate::crosslinks::tests::linker("SDA");
    let any_end = (0..2).find(|&e| sda.ends[e].reacts_with(b'G')).unwrap();
    assert_eq!(site_penalty(12.0, &p, &at(0, false), &sda, any_end), 0.0);
}

#[test]
fn mass_pairs_find_full_and_half_pairs() {
    let cfg = XlConfig::new(Crosslinker::dsso());
    let thiol = cfg.linker.cleavage.as_ref().unwrap().stubs[1].mass;
    let (ma, mb) = (1500.0, 900.0);
    let prec = ma + mb + cfg.linker.intact;
    let dm = |m: f64| DoubletMass {
        mass: m,
        intensity: 1.0,
        n_charges: 1,
        isotope_confirmed: false,
        min_charge: 1,
    };
    let full = mass_pairs(&[], &[dm(ma), dm(mb)], prec, 0.01, 4, &cfg);
    assert_eq!(full.len(), 1);
    assert_eq!(full[0].kind, PairKind::Full);
    assert_eq!(full[0].b_stub_ions, 3);
    assert!((full[0].a_mass - ma).abs() < 1e-9 && (full[0].b_mass - mb).abs() < 1e-9);
    assert!(mass_pairs(&[], &[dm(ma)], prec, 0.01, 4, &cfg).is_empty());
    let stub_peak = [((mb + thiol) / 2.0 + PROTON, 10.0_f32)];
    let half = mass_pairs(&stub_peak, &[dm(ma)], prec, 0.01, 4, &cfg);
    assert_eq!(half.len(), 1);
    assert_eq!(half[0].kind, PairKind::Half);
    assert_eq!(half[0].iso, 0);
    assert_eq!(half[0].b_stub_ions, 1);
}
