//! Tests of `andes quant`: score loading, the library join, donors and protein groups.

use super::extract::{evidence_params, scan_of, split_key};
use super::groups::{feature_records, GroupTally, ProteinGroups};
use super::library::build_and_write_library;
use super::scores::{percolator_scores, rescore_scores, score_source, ScoreLoader, ScoreSource};
use crate::quant::finalize::MbrFeature;
use model::isomers::bundled_unimod_mass;
use model::{IsomerClass, Tolerance};
use output::features::FeatureGroup;
use output::mbr::EvidenceType;
use output::{FeatureCounts, Inference, PercolatorPsm};
use quant::align::AlignParams;
use quant::library::PrecursorMeta;
use quant::ms2::{fragment_ladder, residue_masses};
use std::collections::HashMap;

use super::*;
use quant::library::build_library;

#[test]
fn a_precursor_key_ends_in_its_charge() {
    assert_eq!(split_key("LC[UNIMOD:4]AK/3"), "LC[UNIMOD:4]AK");
}

#[test]
fn library_runs_checks_artifacts_and_duplicates() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    for f in [
        "b.quant.parquet",
        "b.ms1.parquet",
        "b.ms2.parquet",
        "a.quant.parquet",
        "a.ms2.parquet",
    ] {
        std::fs::write(d.join(f), b"x").unwrap();
    }
    let e = library_runs(d, &[]).unwrap_err();
    assert!(e.contains("a.ms1.parquet") && !e.contains("b.ms"), "{e}");
    std::fs::remove_file(d.join("a.ms2.parquet")).unwrap();
    std::fs::write(d.join("a.ms1.parquet"), b"x").unwrap();
    let e = library_runs(d, &[]).unwrap_err();
    assert!(e.contains("a.ms2.parquet"), "{e}");
    std::fs::write(d.join("a.ms2.parquet"), b"x").unwrap();
    assert_eq!(library_runs(d, &[]).unwrap(), vec!["a", "b"]);
    assert_eq!(
        library_runs(d, &["b".into(), "a".into()]).unwrap(),
        vec!["a", "b"]
    );
    let e = library_runs(d, &["a".into(), "a".into()]).unwrap_err();
    assert!(e.contains("twice"), "{e}");
    let e = library_runs(d, &["c".into()]).unwrap_err();
    assert!(
        e.contains("c.quant.parquet") && e.contains("c.ms1.parquet") && e.contains("c.ms2.parquet"),
        "{e}"
    );
}

#[test]
fn score_files_are_merged_or_one_pair_per_run() {
    let stems: Vec<String> = vec!["a".into(), "b".into()];
    let v = |x: &[&str]| x.iter().map(|s| s.to_string()).collect::<Vec<_>>();
    assert_eq!(
        score_source(&v(&["t.tsv"]), &v(&["d.tsv"]), &stems).unwrap(),
        ScoreSource::Merged("t.tsv".into(), "d.tsv".into())
    );
    let s = score_source(&v(&["b=tb", "a=ta"]), &v(&["a=da", "b=db"]), &stems).unwrap();
    let ScoreSource::PerRun(m) = s else {
        panic!("{s:?}")
    };
    assert_eq!(m["a"], ("ta".into(), "da".into()));
    assert_eq!(m["b"], ("tb".into(), "db".into()));
    let e = score_source(&v(&["a=ta"]), &v(&["a=da"]), &stems).unwrap_err();
    assert!(e.contains("no file for runs b"), "{e}");
    let e = score_source(&v(&["a=ta", "b.tsv"]), &v(&["a=da", "b=db"]), &stems).unwrap_err();
    assert!(e.contains("cannot be mixed"), "{e}");
    let e = score_source(&v(&["a=ta", "a=tb"]), &v(&["a=da", "b=db"]), &stems).unwrap_err();
    assert!(e.contains("twice"), "{e}");
    assert!(score_source(&v(&["t1", "t2"]), &v(&["d"]), &stems).is_err());
}

#[test]
fn extraction_tolerance_is_ppm() {
    assert_eq!(evidence_params(None).unwrap().tol_ppm, 10.0);
    assert_eq!(
        evidence_params(Some(Tolerance::Ppm(5.0))).unwrap().tol_ppm,
        5.0
    );
    assert!(evidence_params(Some(Tolerance::Da(0.01))).is_err());
}

#[test]
fn scores_need_the_score_column() {
    let mut t = HashMap::new();
    t.insert(
        "s1".to_string(),
        PercolatorPsm {
            pooled_q_value: None,
            psm_id: "s1".into(),
            q_value: 0.0,
            pep: 0.0,
            score: None,
            peptide: "K.A.R".into(),
            proteins: "P".into(),
        },
    );
    let none = HashMap::new();
    let e = rescore_scores(&t, &none, "t.tsv").err().unwrap();
    assert!(e.contains("t.tsv") && e.contains("score"), "{e}");
    t.get_mut("s1").unwrap().score = Some(2.0);
    let e = rescore_scores(&t, &t.clone(), "t.tsv").err().unwrap();
    assert!(e.contains("both"), "{e}");
    let mut d = t.clone();
    let p = d.remove("s1").unwrap();
    d.insert(
        "s2".into(),
        PercolatorPsm {
            pooled_q_value: None,
            psm_id: "s2".into(),
            ..p
        },
    );
    let m = rescore_scores(&t, &d, "t.tsv").unwrap();
    assert_eq!(m.len(), 2);
    assert!(!m.get("s1").unwrap().is_decoy() && m.get("s2").unwrap().is_decoy());
    assert!(m.get("s3").is_none());
}

#[test]
fn percolator_files_are_read_line_by_line() {
    let dir = tempfile::tempdir().unwrap();
    let write = |name: &str, rows: &[&str]| {
        let path = dir.path().join(name);
        let mut text =
            "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\r\n".to_string();
        for r in rows {
            text.push_str(r);
            text.push_str("\r\n");
        }
        std::fs::write(&path, text).unwrap();
        path
    };
    let t = write(
        "t.tsv",
        &[
            "a_1_1\t2.5\t0.001\t0.01\tK.PEPC+57.02146TIDEK.A\tP1\tP2",
            "",
            "b_2_1\t1.5\t0.002\t0.02\t-.PEPM[+15.99491]K.-\tP3",
            "e_4_1\t1.0\t0.003\t0.03\tR.PEPS+79.96633TK.A\tP4",
        ],
    );
    let d = write("d.tsv", &["c_3_1\t-1.0\t0.5\t0.9\tK.KEDITPEP.A\tDECOY_P1"]);
    let set = percolator_scores(&t, &d, 0).unwrap();
    assert_eq!(set.len(), 4);
    let a = *set.get("a_1_1").unwrap();
    assert_eq!(
        (a.score, a.q, a.pep, a.is_decoy()),
        (2.5, 0.001, 0.01, false)
    );
    let c = *set.get("c_3_1").unwrap();
    assert!(c.is_decoy() && c.score == -1.0);
    assert!(set.get("a_1").is_none());
    let cls = |t: &str| IsomerClass::parse(t, bundled_unimod_mass).unwrap();
    assert!(set.agrees(&a, false, &cls("PEPC[UNIMOD:4]TIDEK")));
    assert!(set.agrees(&a, false, &cls("PEPC[+57.021464]TIDEK")));
    assert!(
        !set.agrees(&a, false, &cls("PEPCTIDEK")),
        "another modification composition"
    );
    assert!(!set.agrees(&a, true, &cls("PEPC[UNIMOD:4]TIDEK")));
    assert!(!set.agrees(&a, false, &cls("PEPC[UNIMOD:4]TIDEKK")));
    let b = *set.get("b_2_1").unwrap();
    assert!(set.agrees(&b, false, &cls("PEPM[UNIMOD:35]K")));
    assert!(!set.agrees(&b, false, &cls("PEPMK")));
    let e = *set.get("e_4_1").unwrap();
    assert!(
        set.agrees(&e, false, &cls("PEPST[UNIMOD:21]K")),
        "a positional isomer"
    );
    assert!(!set.agrees(&e, false, &cls("PEPS[UNIMOD:21]T[UNIMOD:21]K")));
    assert!(!set.agrees(&e, false, &cls("PEPST[+79.95682]K")));

    let dup = write(
        "dup.tsv",
        &[
            "a_1_1\t2.5\t0.001\t0.01\tK.PEPTIDEK.A\tP1",
            "a_1_1\t2.0\t0.001\t0.01\tK.PEPTIDEK.A\tP1",
        ],
    );
    let e = percolator_scores(&dup, &d, 0).err().unwrap();
    assert!(e.contains("dup.tsv") && e.contains("a_1_1"), "{e}");
    let e = percolator_scores(&t, &t, 0).err().unwrap();
    assert!(e.contains("both a target and a decoy"), "{e}");
    let no_score = write("ns.tsv", &["a_1_1\tnan\t0.001\t0.01\tK.PEPTIDEK.A\tP1"]);
    let e = percolator_scores(&no_score, &d, 0).err().unwrap();
    assert!(e.contains("ns.tsv") && e.contains("score"), "{e}");
    let bad = write("bad.tsv", &["a_1_1\t2.5\tx\t0.01\tK.PEPTIDEK.A\tP1"]);
    let e = percolator_scores(&bad, &d, 0).err().unwrap();
    assert!(e.contains("bad.tsv") && e.contains("q-value"), "{e}");
    let e = percolator_scores(&dir.path().join("missing.tsv"), &d, 0)
        .err()
        .unwrap();
    assert!(e.contains("missing.tsv"), "{e}");
}

#[test]
fn feature_records_leave_out_below_threshold_rows() {
    let meta: Vec<PrecursorMeta> = (0..2)
        .map(|i| PrecursorMeta {
            key: format!("PEPT{}IDEK/2", ["", "M[UNIMOD:35]"][i]),
            mz: 500.0 + i as f64,
            charge: 2,
            envelope: [0.55, 0.3, 0.1, 0.05],
            isotope_share: 1.0,
        })
        .collect();
    let library = build_library(&[], meta, 2, 0.01, 0.01, None, &AlignParams::default());
    let stored = StoredLibrary {
        runs: vec!["a".into(), "b".into()],
        max_q: 0.01,
        donor_q: 0.01,
        donor_pep: None,
        library,
        proteins: vec!["P1".into(), "P2".into()],
        spectra: vec![None; 2],
        donor_scans: vec![Vec::new(); 2],
    };
    let feature = |precursor: u32, run: u16, evidence_type, transfer_q| MbrFeature {
        precursor,
        run,
        evidence_type,
        quantity: Some(1e6),
        quantity_se: Some(1e4),
        quantity_normalized: Some(1e6),
        transfer_q,
        quant_q: transfer_q.or(Some(0.003)),
        apex_rt_s: Some(100.0),
        library_q: 0.001,
        mass_err_ppm: Some(0.5),
        claim: quant::claims::Claim::default(),
    };
    let mut features = [
        feature(0, 0, EvidenceType::Identified, None),
        feature(0, 1, EvidenceType::BelowThreshold, Some(0.2)),
        feature(1, 0, EvidenceType::TransferMs1, Some(0.004)),
        feature(1, 1, EvidenceType::IdentifiedBelowThreshold, None),
    ];
    features[2].claim = quant::claims::Claim {
        claimant: Some(0),
        dt_s: Some(0.5),
        identified: true,
        co_apex: true,
        n: 1,
        n_co_apex: 1,
        n_identified_co_apex: 1,
        ..quant::claims::Claim::default()
    };
    let records = feature_records(&stored, &features, &mut 0, None, &mut GroupTally::default());
    let got: Vec<(i64, &str, &str)> = records
        .iter()
        .map(|r| {
            (
                r.feature_id,
                r.run_file_name.as_str(),
                r.mbr.as_ref().unwrap().evidence_type.as_str(),
            )
        })
        .collect();
    assert_eq!(
        got,
        [(1, "a", "identified"), (2, "a", "transfer_ms1")],
        "below-threshold transfers and identified rows are not feature rows"
    );
    assert_eq!(records[1].peptidoform, "PEPTM[UNIMOD:35]IDEK");
    let quant_q: Vec<Option<f64>> = records
        .iter()
        .map(|r| r.mbr.as_ref().unwrap().quant_q)
        .collect();
    assert_eq!(quant_q, [Some(0.003), Some(0.004)]);
    let claims: Vec<(Option<&str>, Option<bool>)> = records
        .iter()
        .map(|r| {
            let m = r.mbr.as_ref().unwrap();
            (m.claimed_by.as_deref(), m.claim_co_apex)
        })
        .collect();
    assert_eq!(
        claims,
        [(None, None), (Some("PEPTIDEK/2"), Some(true))],
        "a transfer reports its claimant's key and co-apex flag, an identified row neither"
    );
    assert_eq!(records[0].proteins[0].0, "P1");
    assert!(records[0].group.is_none());

    let group = |id: usize, lead: &str, q: f64| output::ProteinGroup {
        id,
        lead: lead.into(),
        families: vec![lead.into()],
        accessions: vec![lead.into(), format!("{lead}-2")],
        peptides: 1,
        unique: 0,
        razor: 1,
        psms: 1,
        score: 3.0,
        q_value: q,
        subsumed: Vec::new(),
        peptides_per_accession: Vec::new(),
        is_decoy: false,
    };
    let inf = Inference {
        groups: vec![group(1, "Q1", 0.0), group(2, "Q2", 0.05)],
        peptide_group: [("PEPTIDEK".into(), 0), ("PEPTMIDEK".into(), 1)].into(),
        shared_peptides: ["PEPTIDEK".to_string()].into(),
        decoy_groups: Vec::new(),
    };
    let groups = ProteinGroups {
        inf: &inf,
        max_q: 0.01,
    };
    let mut tally = GroupTally::default();
    let records = feature_records(&stored, &features, &mut 0, Some(&groups), &mut tally);
    assert_eq!(records.len(), 1);
    let r = &records[0];
    assert_eq!(r.feature_id, 1);
    let accs: Vec<&str> = r.proteins.iter().map(|p| p.0.as_str()).collect();
    assert_eq!(accs, ["Q1", "Q1-2"]);
    assert_eq!(
        r.group,
        Some(FeatureGroup {
            anchor: "Q1".into(),
            unique: false,
            global_q: 0.0
        })
    );
    assert_eq!(
        tally.counts[&1],
        FeatureCounts {
            unique: 0,
            total: 1
        }
    );
    assert_eq!((tally.no_group, tally.above_cut), (0, 1));
}

fn one_run_library(
    dir: &Path,
    pep: f64,
    wrong: Option<(&str, &str)>,
) -> (Vec<String>, HashMap<String, PercolatorPsm>) {
    use output::quantities::runs::{write_run_quant_table, RunPsmRow};
    let peptide = |i: usize| -> String {
        let aa = b"ACDEFGHIKLMNPQRSTVWY";
        let tail: String = [aa[i % 20], aa[(i / 20) % 20]]
            .iter()
            .map(|&c| c as char)
            .collect();
        format!("PEPTIDE{tail}K")
    };
    let rows: Vec<RunPsmRow> = (0..200)
        .map(|i| RunPsmRow {
            spec_id: format!("s{i}_{i}_1"),
            rt_s: 100.0 + i as f64,
            precursor_key: format!("{}/2", peptide(i)),
            peptidoform: peptide(i),
            charge: 2,
            mz: 500.0 + i as f64,
            is_decoy: false,
            search_score: 10.0,
            faims_cv: None,
            envelope: vec![0.55, 0.3, 0.1, 0.05],
            isotope_share: 1.0,
            proteins: "P1".into(),
            ion_mobility: None,
        })
        .collect();
    write_run_quant_table(&artifact(dir, "a", TABLE_SUFFIX), "a", &rows, &[]).unwrap();
    let scans: Vec<quant::ms2::Ms2Scan> = rows
        .iter()
        .enumerate()
        .map(|(i, r)| {
            let ladder =
                fragment_ladder(&residue_masses(&r.peptidoform, &HashMap::new()).unwrap(), 2);
            let peaks: Vec<(f32, f32)> = if i % 2 == 0 {
                ladder
                    .iter()
                    .enumerate()
                    .map(|(k, f)| (f.mz as f32, 100.0 + k as f32))
                    .collect()
            } else {
                vec![(150.05, 10.0), (250.5, 20.0)]
            };
            quant::ms2::Ms2Scan {
                native_id: format!("scan={i}"),
                spec_key: format!("s{i}_{i}"),
                rt: r.rt_s,
                precursor_mz: r.mz,
                charge: 2,
                isolation: Some((r.mz, 0.7, 0.7)),
                faims_cv: None,
                peaks,
                reporters: None,
            }
        })
        .collect();
    output::ms2::write_ms2_table(
        &artifact(dir, "a", MS2_SUFFIX),
        &output::ms2::Ms2Table {
            run: "a".into(),
            fragment_tol: Tolerance::Ppm(20.0),
            top_peaks: 100,
            scans,
        },
    )
    .unwrap();
    let scores = rows
        .iter()
        .map(|r| {
            let peptide = if let Some((_, p)) = wrong.filter(|w| w.0 == r.spec_id) {
                p.to_string()
            } else {
                format!("K.{}.A", r.peptidoform)
            };
            let psm = PercolatorPsm {
                pooled_q_value: None,
                psm_id: r.spec_id.clone(),
                q_value: 0.001,
                pep,
                score: Some(100.0 + r.rt_s),
                peptide,
                proteins: "P1".into(),
            };
            (r.spec_id.clone(), psm)
        })
        .collect();
    (vec!["a".to_string()], scores)
}

fn open_peptide(i: usize) -> String {
    let aa = b"ACDEFGHIKLMNPQRSTVWY";
    let tail: String = [aa[i % 20], aa[(i / 20) % 20]]
        .iter()
        .map(|&c| c as char)
        .collect();
    format!("PEPTIDE{tail}K")
}

#[test]
fn scan_of_strips_rank_and_row() {
    let keys: HashMap<&str, usize> = [("a/scan=5_5", 0usize), ("x_1", 1)].into_iter().collect();
    assert_eq!(scan_of(&keys, "a/scan=5_5_1"), Some(0));
    assert_eq!(scan_of(&keys, "a/scan=5_5_1_2"), Some(0));
    assert_eq!(scan_of(&keys, "x_1_1"), Some(1));
    assert_eq!(scan_of(&keys, "a/scan=6_6_1"), None);
    assert_eq!(scan_of(&keys, "nounderscore"), None);
}

#[test]
fn donors_honour_the_pep_cap() {
    let dir = tempfile::tempdir().unwrap();
    let (stems, t) = one_run_library(dir.path(), 0.5, None);
    let build = |pep: Option<f64>| {
        let set = rescore_scores(&t, &HashMap::new(), "t.tsv").unwrap();
        build_and_write_library(
            dir.path(),
            &stems,
            ScoreLoader::Merged(Box::new(move |_| Ok(set))),
            0.01,
            0.01,
            pep,
        )
        .unwrap()
    };
    let open = build(None);
    assert_eq!(open.library.entries.len(), 200);
    let with: Vec<usize> = (0..open.spectra.len())
        .filter(|&id| open.spectra[id].is_some())
        .collect();
    assert_eq!(with.len(), 100);
    for &id in &with {
        let key = &open.library.meta[id].key;
        let i = (0..200)
            .find(|&i| key == &format!("{}/2", open_peptide(i)))
            .unwrap();
        assert_eq!(i % 2, 0, "{key}");
        let s = open.spectra[id].as_ref().unwrap();
        assert_eq!(s.fragments.len(), quant::ms2::LIBRARY_FRAGMENTS);
    }
    assert!(open.donor_scans.iter().all(|d| d.len() == 1));
    assert_eq!(read_library(dir.path()).unwrap().spectra, open.spectra);
    let capped = build(Some(0.01));
    assert!(
        capped.library.entries.is_empty(),
        "donors above the PEP cap"
    );
    assert_eq!(capped.donor_pep, Some(0.01));
    assert_eq!(read_library(dir.path()).unwrap().donor_pep, Some(0.01));
    assert_ne!(library_fingerprint(&open), library_fingerprint(&capped));
}

#[test]
fn mismatched_joins_are_not_donors() {
    for (wrong, joins) in [
        ("K.WRONGPEPTIDEK.A", false),
        ("K.PEPTIDEIAK+8.01420.A", false),
        ("K.PEPTIDEIAK.A", true),
    ] {
        let dir = tempfile::tempdir().unwrap();
        let (stems, t) = one_run_library(dir.path(), 0.001, Some(("s7_7_1", wrong)));
        let set = rescore_scores(&t, &HashMap::new(), "t.tsv").unwrap();
        let stored = build_and_write_library(
            dir.path(),
            &stems,
            ScoreLoader::Merged(Box::new(move |_| Ok(set))),
            0.01,
            0.01,
            None,
        )
        .unwrap();
        let key_of = |e: &quant::library::LibraryEntry| {
            stored.library.meta[e.precursor as usize].key.clone()
        };
        let keys: Vec<String> = stored.library.entries.iter().map(key_of).collect();
        assert_eq!(keys.len(), 199 + usize::from(joins), "{wrong}");
        assert_eq!(
            keys.iter().any(|k| k == "PEPTIDEIAK/2"),
            joins,
            "{wrong}: the PSM is a donor"
        );
    }
}

#[test]
fn a_run_with_mismatched_modifications_fails_the_join_check() {
    let dir = tempfile::tempdir().unwrap();
    let (stems, mut t) = one_run_library(dir.path(), 0.001, None);
    for p in t.values_mut() {
        p.peptide = p.peptide.replace("K.A", "K+8.01420.A");
    }
    let set = rescore_scores(&t, &HashMap::new(), "t.tsv").unwrap();
    let e = build_and_write_library(
        dir.path(),
        &stems,
        ScoreLoader::Merged(Box::new(move |_| Ok(set))),
        0.01,
        0.01,
        None,
    )
    .err()
    .unwrap();
    assert!(
        e.contains("200 of 200") && e.contains("modifications"),
        "{e}"
    );
}

#[test]
fn site_forms_of_one_peptide_are_one_precursor() {
    use output::quantities::runs::{write_run_quant_table, RunPsmRow};
    let dir = tempfile::tempdir().unwrap();
    let s_form = "PEPS[UNIMOD:21]TIDEK";
    let t_form = "PEPST[UNIMOD:21]IDEK";
    let row = |spec_id: String, peptidoform: &str, rt_s: f64| RunPsmRow {
        spec_id,
        rt_s,
        precursor_key: format!("{peptidoform}/2"),
        peptidoform: peptidoform.to_string(),
        charge: 2,
        mz: 520.0,
        is_decoy: false,
        search_score: 10.0,
        faims_cv: None,
        envelope: vec![0.55, 0.3, 0.1, 0.05],
        isotope_share: 1.0,
        proteins: "P1".into(),
        ion_mobility: None,
    };
    let mut scores = HashMap::new();
    let stems: Vec<String> = ["a", "b", "c"].map(String::from).to_vec();
    for (r, stem) in stems.iter().enumerate() {
        let shift = 5.0 * r as f64;
        let mut rows: Vec<RunPsmRow> = (0..150)
            .map(|i| {
                row(
                    format!("{stem}{i}_{i}_1"),
                    &open_peptide(i),
                    100.0 + 10.0 * i as f64 + shift,
                )
            })
            .collect();
        let (form, score) = match r {
            1 => (t_form, 600.0),
            _ => (s_form, 500.0 - 100.0 * r as f64),
        };
        rows.push(row(format!("{stem}p_900_1"), form, 900.0 + shift));
        for (i, x) in rows.iter().enumerate() {
            let pin = x.peptidoform.replace("[UNIMOD:21]", "+79.96633");
            let score = if i == 150 { score } else { 100.0 + x.rt_s };
            scores.insert(
                x.spec_id.clone(),
                PercolatorPsm {
                    pooled_q_value: None,
                    psm_id: x.spec_id.clone(),
                    q_value: 0.001,
                    pep: 0.001,
                    score: Some(score),
                    peptide: format!("K.{pin}.A"),
                    proteins: "P1".into(),
                },
            );
        }
        write_run_quant_table(
            &artifact(dir.path(), stem, TABLE_SUFFIX),
            stem,
            &rows,
            &[("UNIMOD:21".to_string(), 79.966331)],
        )
        .unwrap();
        output::ms2::write_ms2_table(
            &artifact(dir.path(), stem, MS2_SUFFIX),
            &output::ms2::Ms2Table {
                run: stem.clone(),
                fragment_tol: Tolerance::Ppm(20.0),
                top_peaks: 100,
                scans: Vec::new(),
            },
        )
        .unwrap();
    }
    let set = rescore_scores(&scores, &HashMap::new(), "t.tsv").unwrap();
    let stored = build_and_write_library(
        dir.path(),
        &stems,
        ScoreLoader::Merged(Box::new(move |_| Ok(set))),
        0.01,
        0.01,
        None,
    )
    .unwrap();
    let lib = &stored.library;
    assert_eq!(lib.meta.len(), 151, "150 peptides and one phosphopeptide");
    let phospho: Vec<usize> = (0..lib.meta.len())
        .filter(|&id| lib.meta[id].key.contains("UNIMOD:21"))
        .collect();
    assert_eq!(phospho.len(), 1, "one precursor for both site forms");
    let id = phospho[0];
    assert_eq!(
        lib.meta[id].key,
        format!("{t_form}/2"),
        "keyed by its best PSM"
    );
    let entries: Vec<&quant::library::LibraryEntry> = lib
        .entries
        .iter()
        .filter(|e| e.precursor as usize == id)
        .collect();
    assert_eq!(entries.len(), 1);
    assert_eq!(entries[0].donors, vec![0, 1, 2], "a donor in every run");
    let key_of = |id: usize| lib.meta[id].key.clone();
    let read = read_library(dir.path()).unwrap();
    assert_eq!(read.library.meta.len(), 151);
    assert_eq!(read.library.meta[id].key, key_of(id));
}

#[test]
fn a_rescoring_peptide_without_a_mass_fails() {
    let mut t = HashMap::new();
    t.insert(
        "s1".to_string(),
        PercolatorPsm {
            pooled_q_value: None,
            psm_id: "s1".into(),
            q_value: 0.0,
            pep: 0.0,
            score: Some(1.0),
            peptide: "K.PEPC[NotAModification]K.A".into(),
            proteins: "P".into(),
        },
    );
    let e = rescore_scores(&t, &HashMap::new(), "t.tsv").err().unwrap();
    assert!(e.contains("NotAModification") && e.contains("s1"), "{e}");
}
