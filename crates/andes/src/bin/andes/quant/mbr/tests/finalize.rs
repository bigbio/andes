use output::mbr::{
    write_mbr_validation, Accepted, EvidenceType, MbrClaim, MbrParquetWriter, MbrRow,
    MbrValidation, VALIDATION_Q,
};
use quant::align::AlignParams;
use quant::library::{build_library, PrecursorMeta, PsmIn};
use quant::normalize::{fit as fit_normalization, NormParams};

use super::*;

type OutEntry<'a> = (u32, u16, MbrRow<'a>, (Option<f64>, claims::Claim));
use crate::quant::finalize::{
    finalize, finalize_stream, for_each_scored, mismatch, pseudo_example, running_mean_q,
    table_features, transfer_quant_q, wrong_features, wrong_label, BlockContext, EntryIndex,
    ExtractRuns, FinalizeLimits, FinalizeParams, FinalizeSummary, MbrFeature, Mismatch, Table,
    WrongCalibration, SCORE_CHUNK,
};
use crate::transfers::{
    row_unit, train_bagged_folds, train_folds, transfer_train_params, Estimate, FdrRow,
    FoldScores, QValues, TrainingSample, TARGET,
};

const N: u32 = 900;
const WRONG_EVERY: u32 = 10;

struct Rng(u64);
impl Rng {
    fn unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((self.0 >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }
    fn normal(&mut self) -> f64 {
        let (u1, u2) = (self.unit(), self.unit());
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

fn evidence(g: &mut Rng, present: bool, t0: f64) -> Evidence {
    let (llr, cos, coh) = if present {
        (60.0 + 10.0 * g.normal(), 0.97 + 0.01 * g.normal(), 0.9)
    } else {
        (
            3.0 * g.normal().abs(),
            0.6 + 0.2 * g.normal(),
            0.2 * g.normal(),
        )
    };
    let amplitude = if present { 1e6 } else { 2e3 * g.unit() };
    Evidence {
        amplitude,
        amplitude_se: 1e4,
        quantity: amplitude * 10.0,
        quantity_se: 1e5,
        llr,
        rt_z: if present {
            0.3 * g.normal()
        } else {
            2.0 * g.normal()
        },
        width_ratio: 1.0 + 0.1 * g.normal(),
        envelope_cos: cos,
        coherence: coh,
        m1_ratio: 0.05 * g.unit(),
        mass_err_ppm: if present {
            g.normal()
        } else {
            5.0 * g.normal()
        },
        n_points: 12,
        isotopes_used: 0b1111,
        t0,
        w: 4.0,
        iso_heights: [0.55, 0.3, 0.11, 0.04].map(|v| v * amplitude),
        tailing: 0.0,
        trace_quantity: amplitude * 12.5,
        trace_quantity_se: 1.25e5,
    }
}

fn row(p: u32, run: u16, kind: Hypothesis, role: Role, e: Option<Evidence>) -> ExtractRow {
    ExtractRow {
        precursor: p,
        run,
        kind,
        role,
        evidence: e,
        n_donor_runs: 1,
        library_q: 0.001,
        ms2: None,
        claim: claims::Claim::default(),
        shared: None,
    }
}

fn extract() -> Vec<Vec<ExtractRow>> {
    let mut g = Rng(5);
    let mut run0 = Vec::new();
    let mut run1 = Vec::new();
    for p in 0..N {
        let apex = 100.0 + p as f64;
        run0.push(row(
            p,
            0,
            Hypothesis::Target,
            Role::Identified,
            Some(evidence(&mut g, true, apex)),
        ));
        let pseudo_t0 = if p % WRONG_EVERY == 0 {
            apex + 60.0
        } else {
            apex + 0.5
        };
        run0.push(row(
            p,
            0,
            Hypothesis::Target,
            Role::PseudoTransfer,
            Some(evidence(&mut g, true, pseudo_t0)),
        ));
        let target = (p != 1).then(|| evidence(&mut g, p % 2 == 0, apex));
        run1.push(row(p, 1, Hypothesis::Target, Role::Transfer, target));
        run1.push(row(
            p,
            1,
            Hypothesis::MassShiftDecoy(7),
            Role::Transfer,
            Some(evidence(&mut g, false, apex)),
        ));
        run1.push(row(
            p,
            1,
            Hypothesis::RtShiftDecoy,
            Role::Transfer,
            Some(evidence(&mut g, false, apex + 300.0)),
        ));
    }
    vec![run0, run1]
}

fn library() -> Library {
    let meta = (0..N)
        .map(|p| PrecursorMeta {
            key: format!("PEP{p}/2"),
            mz: 500.0 + p as f64,
            charge: 2,
            envelope: [0.55, 0.30, 0.11, 0.04],
            isotope_share: 1.0,
        })
        .collect();
    build_library(&[], meta, 2, 0.01, 0.01, None, &AlignParams::default())
}

fn stems() -> Vec<String> {
    vec!["a".into(), "b".into()]
}

fn table(path: &Path) -> Vec<Vec<String>> {
    use output::mbr::{read_mbr_parquet, read_mbr_validation, ValidationCell, MBR_COLUMNS};
    let num = |v: Option<f64>| v.map_or_else(|| "NA".to_string(), |v| v.to_string());
    if path.to_string_lossy().ends_with(".mbr_validation.parquet") {
        let cells = read_mbr_validation(path).unwrap();
        let header = cells.iter().map(|(n, _)| n.clone()).collect();
        let values = cells
            .iter()
            .map(|(_, c)| match *c {
                ValidationCell::Count(n) => n.to_string(),
                ValidationCell::Value(v) => num(v),
            })
            .collect();
        return vec![header, values];
    }
    let mut t = vec![MBR_COLUMNS.iter().map(|c| c.to_string()).collect()];
    for r in read_mbr_parquet(path).unwrap() {
        t.push(vec![
            r.precursor_key,
            r.run,
            r.evidence_type,
            num(r.quantity),
            num(r.quantity_se),
            num(r.quantity_normalized),
            num(r.transfer_q),
            num(r.quant_q),
            num(r.apex_rt_s),
            num(Some(r.library_q)),
            num(r.measurement_q),
        ]);
    }
    t
}

fn library_with(identity: &dyn Fn(u32) -> bool, donor_q: &dyn Fn(u32) -> f64) -> Library {
    let mut lib = library();
    lib.entries = (0..N)
        .map(|p| LibraryEntry {
            precursor: p,
            identity_decoy: identity(p),
            q: 0.001,
            donor_psm: quant::library::DonorPsm {
                score: 10.0,
                q: donor_q(p),
                pep: 1e-4,
            },
            consensus_rt: f64::NAN,
            faims_cv: None,
            ion_mobility: None,
            donors: vec![0],
            best_rt: vec![(0, 100.0 + p as f64)],
        })
        .collect();
    lib
}

fn on_real_features(p: u32) -> bool {
    p % 50 == 2
}

fn confident(_: u32) -> f64 {
    1e-4
}

#[test]
fn identity_decoys_on_real_features_count_in_the_transfer_q() {
    let dir = tempfile::tempdir().unwrap();
    let (mbr, val) = mbr_paths(&dir.path().join("x"));
    let lib = library_with(&on_real_features, &confident);
    let s = finalize(&stems(), &lib, &extract(), 0.2, None, 7, &mbr, &val).unwrap();
    let v = &s.validation;
    assert_eq!((v.target_entries, v.identity_entries), (882, 18));
    assert_eq!(v.decoy_transfers[2], 18);
    assert!(v.at_threshold.decoys[2] >= 12, "{v:?}");
    assert!(v.at_q01.wrong_peak_only >= 380, "{v:?}");
    assert!(2 * v.at_q01.pooled < v.at_q01.wrong_peak_only, "{v:?}");
    assert!(v.at_threshold.pooled >= 380, "{v:?}");
    let plain = library_with(&|_| false, &confident);
    let s = finalize(&stems(), &plain, &extract(), 0.2, None, 7, &mbr, &val).unwrap();
    let v = &s.validation;
    assert_eq!(v.decoy_transfers[2], 0);
    assert_eq!(v.at_q01.pooled, v.at_q01.wrong_peak_only);
}

#[test]
fn identity_decoy_rows_never_reach_the_outputs() {
    let dir = tempfile::tempdir().unwrap();
    let (mbr, val) = mbr_paths(&dir.path().join("x"));
    let lib = library_with(&on_real_features, &confident);
    let s = finalize(&stems(), &lib, &extract(), 0.5, None, 7, &mbr, &val).unwrap();
    let t = table(&mbr);
    let keys: Vec<&str> = t[1..].iter().map(|r| r[0].as_str()).collect();
    assert_eq!(keys.len(), 2 * (N as usize - 18));
    for p in (0..N).filter(|&p| on_real_features(p)) {
        let key = format!("PEP{p}/2");
        assert!(!keys.contains(&key.as_str()), "{key}");
    }
    assert_eq!(s.rows_written, keys.len());
    assert!(s.features.iter().all(|f| !on_real_features(f.precursor)));
    assert_eq!(s.validation.identified, N as usize - 18);
}

#[test]
fn entries_with_a_marginal_best_donor_psm_are_still_scored() {
    let dir = tempfile::tempdir().unwrap();
    let (mbr, val) = mbr_paths(&dir.path().join("x"));
    let lib = library_with(&|_| false, &|p| if p > 400 { 0.005 } else { 0.001 });
    let s = finalize(&stems(), &lib, &extract(), 0.01, None, 7, &mbr, &val).unwrap();
    let scored_high = s
        .features
        .iter()
        .filter(|f| f.run == 1 && f.precursor > 400 && f.transfer_q.is_some())
        .count();
    assert!(scored_high > 0);
    assert_eq!(s.validation.identified, N as usize);
}

#[test]
fn finalize_writes_one_row_per_target_pair_with_evidence_types() {
    let dir = tempfile::tempdir().unwrap();
    let (mbr, val) = mbr_paths(&dir.path().join("x"));
    let s = finalize(&stems(), &library(), &extract(), 0.01, None, 7, &mbr, &val).unwrap();
    let t = table(&mbr);
    assert_eq!(t[0][2], "evidence_type");
    let body = &t[1..];
    assert_eq!(body.len(), 2 * N as usize);
    assert_eq!(s.rows_written, 2 * N as usize);
    let (mut true_acc, mut false_acc) = (0, 0);
    for (k, r) in body.iter().enumerate() {
        let p = (k / 2) as u32;
        assert_eq!(r[0], format!("PEP{p}/2"));
        if r[1] == "a" {
            assert_eq!((r[2].as_str(), r[5].as_str()), ("identified", "NA"));
            continue;
        }
        assert_eq!(r[1], "b");
        assert!(r[2] == "transfer_ms1" || r[2] == "below_threshold", "{r:?}");
        if p == 1 {
            assert_eq!(&r[2..7], ["below_threshold", "NA", "NA", "NA", "NA"]);
        } else {
            assert_ne!(r[3], "NA", "a below-threshold row keeps its quantity");
        }
        if r[2] == "transfer_ms1" {
            if p % 2 == 0 {
                true_acc += 1;
            } else {
                false_acc += 1;
            }
        }
    }
    assert!(
        true_acc >= 400,
        "{true_acc} of 450 present transfers accepted"
    );
    assert!(false_acc <= 5, "{false_acc} absent transfers accepted");
    assert_eq!(s.no_features, 1);
}

#[test]
fn finalize_reports_pseudo_transfer_validation() {
    let dir = tempfile::tempdir().unwrap();
    let (mbr, val) = mbr_paths(&dir.path().join("x"));
    let s = finalize(&stems(), &library(), &extract(), 0.01, None, 7, &mbr, &val).unwrap();
    let v = &s.validation;
    assert_eq!(v.pseudo_transfers, N as usize);
    assert_eq!(v.wrong_peak_pairs, N as usize);
    assert_eq!(v.wrong_peak, (N / WRONG_EVERY) as usize);
    assert_eq!(v.wrong_measurement, v.wrong_peak);
    assert_eq!(v.at_q01.wrong_measurement, v.at_q01.wrong_peak);
    assert!(v.pseudo_acceptance(&v.at_q01).unwrap() > 0.9, "{v:?}");
    assert_eq!(v.at_q01, v.at_threshold);
    assert_eq!(v.transfers, N as usize);
    assert_eq!(v.transfers_scored, N as usize - 1);
    assert!(v.at_q01.mass_only > 0 && v.at_q01.rt_only > 0);
    let t = table(&val);
    assert_eq!(t.len(), 2);
    let get = |name: &str| &t[1][t[0].iter().position(|h| h == name).unwrap()];
    assert_eq!(get("wrong_peak_rate"), "0.1");
    assert_eq!(get("wrong_measurement_rate"), "0.1");
    assert_eq!(get("transfers"), &N.to_string());
    assert_eq!(get("accepted_pooled_q01"), &v.at_q01.pooled.to_string());

    let loose = finalize(&stems(), &library(), &extract(), 0.5, None, 7, &mbr, &val).unwrap();
    let w = &loose.validation;
    assert_eq!(w.q_threshold, 0.5);
    assert_eq!(w.at_q01, v.at_q01);
    assert!(w.at_threshold.pooled > w.at_q01.pooled, "{w:?}");
    let t = table(&val);
    let get = |name: &str| &t[1][t[0].iter().position(|h| h == name).unwrap()];
    assert_eq!(get("accepted_pooled_q01"), &v.at_q01.pooled.to_string());
    assert_eq!(
        get("accepted_pooled_at_threshold"),
        &w.at_threshold.pooled.to_string()
    );
}

#[test]
fn finalize_is_deterministic_for_a_seed() {
    let dir = tempfile::tempdir().unwrap();
    let runs = extract();
    let (m1, v1) = mbr_paths(&dir.path().join("one"));
    let (m2, v2) = mbr_paths(&dir.path().join("two"));
    finalize(&stems(), &library(), &runs, 0.01, None, 3, &m1, &v1).unwrap();
    finalize(&stems(), &library(), &runs, 0.01, None, 3, &m2, &v2).unwrap();
    assert_eq!(std::fs::read(&m1).unwrap(), std::fs::read(&m2).unwrap());
    assert_eq!(std::fs::read(&v1).unwrap(), std::fs::read(&v2).unwrap());
}

fn library_with_entries() -> Library {
    let psms: Vec<PsmIn> = (0..N)
        .map(|p| PsmIn {
            q_gate: f64::NAN,
            run: 0,
            precursor: p,
            rt: 100.0 + p as f64,
            score: 10.0 + p as f64,
            q: 0.0001,
            pep: 0.0001,
            is_decoy: false,
            faims_cv: None,
            ion_mobility: None,
        })
        .collect();
    build_library(
        &psms,
        library().meta,
        2,
        1.0,
        0.01,
        None,
        &AlignParams::default(),
    )
}

struct One<'a>(&'a [ExtractRow]);

impl One<'_> {
    fn for_each(
        &self,
        mut f: impl FnMut(&[ExtractRow]) -> Result<(), String>,
    ) -> Result<(), String> {
        f(self.0)
    }
}

#[allow(clippy::too_many_arguments)]
fn reference(
    stems: &[String],
    lib: &Library,
    runs: &[Vec<ExtractRow>],
    mbr_fdr: f64,
    quant_feature_fdr: Option<f64>,
    seed: u64,
    mbr_path: &Path,
    validation_path: &Path,
) -> Result<FinalizeSummary, String> {
    let mut all_features = Vec::new();
    let on_features = &mut |f: Vec<MbrFeature>| -> Result<(), String> {
        all_features.extend(f);
        Ok(())
    };
    if !(mbr_fdr > 0.0 && mbr_fdr <= 1.0) {
        return Err(format!(
            "transfer q-value threshold {mbr_fdr} is outside (0, 1]"
        ));
    }
    if let Some(t) = quant_feature_fdr.filter(|t| !(*t > 0.0 && *t <= 1.0)) {
        return Err(format!(
            "identified quality q-value threshold {t} is outside (0, 1]"
        ));
    }
    let maps = &lib.alignment.maps;
    let idx = EntryIndex::new(lib, stems.len());
    let all: Vec<ExtractRow> = runs.concat();
    let shards = One(&all);
    let mut summary = FinalizeSummary::default();

    type Key = (u16, u32, u64);
    let mut t_sample: TrainingSample<Key> = TrainingSample::new(usize::MAX, seed);
    let mut i_sample: TrainingSample<Key> = TrainingSample::new(usize::MAX, seed);
    let mut seq = 0u64;
    shards.for_each(|rows| {
        if let Some(r) = rows
            .iter()
            .find(|r| r.run as usize >= stems.len() || r.precursor as usize >= lib.meta.len())
        {
            return Err(format!(
                "extract row of precursor {} in run {} is outside the {} precursors and {} runs",
                r.precursor,
                r.run,
                lib.meta.len(),
                stems.len()
            ));
        }
        let shapes = DonorShapes::from_rows(rows.iter(), maps);
        for (c, chunk) in rows.chunks(SCORE_CHUNK).enumerate() {
            let feats = table_features(chunk, &idx, &shapes, maps, |t: Table| t != Table::Pseudo);
            for (i, t, fam, x) in feats {
                seq += 1;
                let r = &rows[c * SCORE_CHUNK + i];
                let Some(x) = x else {
                    match t {
                        Table::Identified => summary.identified_no_features += 1,
                        _ => summary.no_features += 1,
                    }
                    continue;
                };
                let key = (0u16, 0u32, seq);
                let row = FdrRow {
                    precursor: r.precursor,
                    is_decoy: fam != 0,
                    family: fam,
                    x,
                };
                match t {
                    Table::Identified => i_sample.push(key, row_unit(r.precursor, r.run), row),
                    _ => t_sample.push(key, row_unit(r.precursor, r.run), row),
                }
            }
        }
        Ok(())
    })?;
    let t_models =
        train_bagged_folds(&[&t_sample.into_items()], seed, &transfer_train_params());
    let i_models = train_folds(&i_sample.into_items(), seed);
    let mut w_sample: TrainingSample<u64> = TrainingSample::new(usize::MAX, seed);
    let mut w_seq = 0u64;
    shards.for_each(|rows| {
        let ctx = BlockContext::new(rows, maps)?;
        for r in rows {
            w_seq += 1;
            if let Some((y, x)) = pseudo_example(r, &ctx, &idx) {
                w_sample.push(
                    w_seq,
                    row_unit(r.precursor, r.run),
                    FdrRow {
                        precursor: r.precursor,
                        is_decoy: y,
                        family: 0,
                        x,
                    },
                );
            }
        }
        Ok(())
    })?;
    let w_models = train_folds(&w_sample.into_items(), seed);
    let mut w_held: Vec<(f32, bool)> = Vec::new();
    shards.for_each(|rows| {
        let ctx = BlockContext::new(rows, maps)?;
        let ex: Vec<(bool, Vec<f32>, u32)> = rows
            .iter()
            .filter_map(|r| pseudo_example(r, &ctx, &idx).map(|(y, x)| (y, x, r.precursor)))
            .collect();
        let precursors: Vec<u32> = ex.iter().map(|e| e.2).collect();
        let xs: Vec<&[f32]> = ex.iter().map(|e| e.1.as_slice()).collect();
        let mut margins = Vec::new();
        w_models.score(&precursors, &xs, &mut margins);
        w_held.extend(
            margins
                .iter()
                .zip(&ex)
                .filter_map(|(m, e)| m.map(|m| (m, e.0))),
        );
        Ok(())
    })?;
    let wrong_cal = WrongCalibration::fit(w_held);

    let mut t_scores = FoldScores::with_strata(seed, idx.strata());
    let mut i_scores = FoldScores::new(seed);
    shards.for_each(|rows| {
        let shapes = DonorShapes::from_rows(rows.iter(), maps);
        let want = |t: Table| t != Table::Pseudo;
        for_each_scored(
            rows,
            &idx,
            &shapes,
            maps,
            &t_models,
            &i_models,
            want,
            |i, t, fam, s| {
                let p = rows[i].precursor;
                match t {
                    Table::Identified => i_scores.add(p, fam, s),
                    _ => t_scores.add(p, fam, s),
                }
            },
        );
        Ok(())
    })?;
    let fdr = t_scores.finish(&t_models.failed_folds);
    let id_fdr = i_scores.finish(&i_models.failed_folds);

    let mut cells: TrainingSample<u64, (NormCell, Option<f64>)> =
        TrainingSample::new(usize::MAX, seed);
    let mut cell_seq = 0u64;
    let empty_shapes = DonorShapes::default();
    shards.for_each(|rows| {
        let k = relative_widths(rows.iter());
        let want = |t: Table| t == Table::Identified;
        for_each_scored(
            rows,
            &idx,
            &empty_shapes,
            maps,
            &t_models,
            &i_models,
            want,
            |i, _, fam, s| {
                let r = &rows[i];
                if r.kind != Hypothesis::Target
                    || !id_fdr
                        .q(Estimate::Pooled, r.precursor, fam, s)
                        .is_some_and(|q| q <= VALIDATION_Q)
                {
                    return;
                }
                let Some(e) = r.evidence.as_ref() else { return };
                let k = k.get(&r.precursor).copied();
                let q = reported_quantity(e, k, false, 0.0).0;
                if q > 0.0 && q.is_finite() {
                    cell_seq += 1;
                    cells.push(
                        cell_seq,
                        u64::from(r.precursor),
                        (
                            NormCell {
                                precursor: r.precursor,
                                run: r.run,
                                log2_quantity: q.log2(),
                                rt: e.t0,
                            },
                            width_deviation(e, k, false),
                        ),
                    );
                }
            },
        );
        Ok(())
    })?;
    let (cells, weights) = weigh_widths(cells.into_items(), stems.len());
    let norm = fit_normalization(&cells, stems.len(), &NormParams::default());

    let accepted = |q: Option<f64>| q.is_some_and(|q| q <= mbr_fdr);
    let thresholds = [VALIDATION_Q, mbr_fdr];
    let mut at = [Accepted::default(), Accepted::default()];
    let identity_entries = lib.entries.iter().filter(|e| e.identity_decoy).count();
    let mut v = MbrValidation {
        q_threshold: mbr_fdr,
        quant_feature_threshold: quant_feature_fdr,
        target_entries: lib.entries.len() - identity_entries,
        identity_entries,
        ..MbrValidation::default()
    };
    let count_decoy =
        |v: &mut MbrValidation, at: &mut [Accepted; 2], q: Option<f64>, f: usize| {
            let Some(q) = q else { return };
            v.decoy_transfers[f] += 1;
            for (a, &t) in at.iter_mut().zip(&thresholds) {
                a.decoys[f] += usize::from(q <= t);
            }
        };
    let mut writer = MbrParquetWriter::create(mbr_path)?;
    shards.for_each(|rows| {
        let shapes = DonorShapes::from_rows(rows.iter(), maps);
        let k = relative_widths(rows.iter());
        let mut identified: HashMap<(u32, u16), &Evidence> = HashMap::new();
        for r in rows {
            if let (Role::Identified, Hypothesis::Target, Some(e)) =
                (r.role, r.kind, &r.evidence)
            {
                let previous = identified.insert((r.precursor, r.run), e);
                debug_assert!(
                    previous.is_none(),
                    "two identified rows of precursor {} in run {}",
                    r.precursor,
                    r.run
                );
            }
        }
        let mut scored: Vec<Option<(u32, Option<f32>)>> = vec![None; rows.len()];
        for_each_scored(
            rows,
            &idx,
            &shapes,
            maps,
            &t_models,
            &i_models,
            |_: Table| true,
            |i, _, fam, s| {
                scored[i] = Some((fam, s));
            },
        );
        let q_of = |q: &QValues, e: Estimate, i: usize| {
            scored[i].and_then(|(fam, s)| q.q(e, rows[i].precursor, fam, s))
        };
        let reported = |r: &ExtractRow| {
            r.evidence.as_ref().filter(|e| e.amplitude > 0.0).map(|e| {
                whole_pattern(
                    reported_quantity(
                        e,
                        k.get(&r.precursor).copied(),
                        false,
                        row_width_weight(r, &weights),
                    ),
                    lib.meta[r.precursor as usize].isotope_share,
                )
            })
        };
        let mut out: Vec<OutEntry<'_>> = Vec::new();
        let mut p_wrong: Vec<f64> = vec![1.0; rows.len()];
        {
            let at: Vec<usize> = (0..rows.len())
                .filter(|&i| {
                    rows[i].role == Role::Transfer
                        && matches!(idx.table(&rows[i]), Some((Table::Transfer, TARGET)))
                })
                .collect();
            let xs: Vec<Option<Vec<f32>>> = at
                .iter()
                .map(|&i| rows[i].evidence.as_ref().and_then(wrong_features))
                .collect();
            let ok: Vec<usize> = (0..at.len()).filter(|&k| xs[k].is_some()).collect();
            let precursors: Vec<u32> = ok.iter().map(|&k| rows[at[k]].precursor).collect();
            let x_refs: Vec<&[f32]> = ok.iter().map(|&k| xs[k].as_deref().unwrap()).collect();
            let mut margins = Vec::new();
            w_models.score(&precursors, &x_refs, &mut margins);
            for (&k, m) in ok.iter().zip(margins) {
                if let (Some(m), Some(c)) = (m, &wrong_cal) {
                    p_wrong[at[k]] = c.at(m);
                }
            }
        }
        let mut wrongs: Vec<(usize, f64)> = Vec::new();
        for (i, r) in rows.iter().enumerate() {
            if idx.is_identity(r.precursor) {
                if (r.role, r.kind) == (Role::Transfer, Hypothesis::Target) {
                    count_decoy(&mut v, &mut at, q_of(&fdr, Estimate::Pooled, i), 2);
                }
                continue;
            }
            if r.kind != Hypothesis::Target {
                if r.role == Role::Transfer {
                    let f = family(r.kind) as usize;
                    count_decoy(&mut v, &mut at, q_of(&fdr, Estimate::Pooled, i), f - 1);
                    if let Some(m) = &r.ms2 {
                        v.ms2_rows[f] += 1;
                        v.ms2_supported[f] += usize::from(m.supports());
                    }
                }
                continue;
            }
            let e = r.evidence.as_ref();
            let quantity = reported(r);
            let quantity_normalized = quantity
                .zip(e)
                .and_then(|((q, _), e)| norm.normalize(r.run, r.precursor, e.t0, q));
            let (evidence_type, transfer_q, quant_q) = match r.role {
                Role::Identified => {
                    let q = q_of(&id_fdr, Estimate::Pooled, i);
                    v.identified += 1;
                    v.identified_scored += usize::from(q.is_some());
                    v.identified_q01 += usize::from(q.is_some_and(|q| q <= VALIDATION_Q));
                    let t = match quant_feature_fdr {
                        Some(t) if !q.is_some_and(|q| q <= t) => {
                            v.identified_below_threshold += 1;
                            EvidenceType::IdentifiedBelowThreshold
                        }
                        _ => EvidenceType::Identified,
                    };
                    (t, None, q)
                }
                Role::Transfer => {
                    let q = q_of(&fdr, Estimate::Pooled, i);
                    v.transfers += 1;
                    v.transfers_scored += usize::from(q.is_some());
                    let (mass, rt, wrong) = (
                        q_of(&fdr, Estimate::MassOnly, i),
                        q_of(&fdr, Estimate::RtOnly, i),
                        q_of(&fdr, Estimate::WrongPeak, i),
                    );
                    for (a, &t) in at.iter_mut().zip(&thresholds) {
                        let ok = |q: Option<f64>| usize::from(q.is_some_and(|q| q <= t));
                        a.pooled += ok(q);
                        a.mass_only += ok(mass);
                        a.rt_only += ok(rt);
                        a.wrong_peak_only += ok(wrong);
                    }
                    let supported = r.ms2.is_some_and(|m| m.supports());
                    if let Some(m) = &r.ms2 {
                        v.ms2_rows[0] += 1;
                        v.ms2_supported[0] += usize::from(m.supports());
                    }
                    let t = match (accepted(q), supported) {
                        (true, true) => {
                            v.transfer_ms2 += 1;
                            EvidenceType::TransferMs2
                        }
                        (true, false) => EvidenceType::TransferMs1,
                        (false, _) => EvidenceType::BelowThreshold,
                    };
                    (t, q, q)
                }
                Role::PseudoTransfer => {
                    let q = scored[i].and_then(|(fam, s)| fdr.pseudo_q(r.precursor, fam, s));
                    v.pseudo_transfers += 1;
                    if let Some(m) = &r.ms2 {
                        v.ms2_pseudo_rows += 1;
                        v.ms2_pseudo_supported += usize::from(m.supports());
                    }
                    v.pseudo_scored += usize::from(q.is_some());
                    let wrong = e
                        .zip(identified.get(&(r.precursor, r.run)))
                        .and_then(|(p, id)| mismatch(p, id));
                    if let Some(m) = wrong {
                        v.wrong_peak_pairs += 1;
                        v.wrong_peak += usize::from(m.apex);
                        v.wrong_measurement += usize::from(m.wrong());
                    }
                    for (a, &t) in at.iter_mut().zip(&thresholds) {
                        if !q.is_some_and(|q| q <= t) {
                            continue;
                        }
                        a.pseudo += 1;
                        if let Some(m) = wrong {
                            a.pairs += 1;
                            a.wrong_peak += usize::from(m.apex);
                            a.wrong_measurement += usize::from(m.wrong());
                        }
                    }
                    continue;
                }
            };
            if r.role == Role::Transfer {
                wrongs.push((out.len(), p_wrong[i]));
            }
            out.push((
                r.precursor,
                r.run,
                MbrRow {
                    precursor_key: &lib.meta[r.precursor as usize].key,
                    run: &stems[r.run as usize],
                    evidence_type,
                    quantity: quantity.map(|q| q.0),
                    quantity_se: quantity.map(|q| q.1),
                    quantity_normalized,
                    transfer_q,
                    quant_q,
                    apex_rt_s: e.map(|e| e.t0),
                    library_q: r.library_q,
                    claim: (r.role == Role::Transfer).then(|| MbrClaim {
                        claimed_by: r.claim.claimant.map(|p| lib.meta[p as usize].key.as_str()),
                        dt_s: r.claim.dt_s.map(f64::from),
                        identified: r.claim.identified,
                        co_apex: r.claim.co_apex,
                        n_co_apex: u32::from(r.claim.n_co_apex),
                        n_identified_co_apex: u32::from(r.claim.n_identified_co_apex),
                    }),
                    measurement_q: None,
                    transfer_pep: if r.role == Role::Transfer {
                        scored[i].and_then(|(fam, s)| fdr.pep(r.precursor, fam, s))
                    } else {
                        None
                    },
                },
                (
                    e.map(|e| e.mass_err_ppm),
                    if r.role == Role::Transfer {
                        r.claim
                    } else {
                        claims::Claim::default()
                    },
                ),
            ));
        }
        let mut by_run: HashMap<u16, Vec<(f64, usize)>> = HashMap::new();
        for &(k, p) in &wrongs {
            if matches!(
                out[k].2.evidence_type,
                EvidenceType::TransferMs2 | EvidenceType::TransferMs1
            ) {
                by_run.entry(out[k].1).or_default().push((p, k));
            }
        }
        for pool in by_run.values_mut() {
            running_mean_q(pool, |k, q| {
                let r = &mut out[k].2;
                r.measurement_q = Some(q);
                r.quant_q = Some(transfer_quant_q(r.transfer_q, q));
            });
        }
        out.sort_by_key(|&(p, run, ..)| (p, run));
        let mut features = Vec::with_capacity(out.len());
        for (p, run, r, (mass_err, claim)) in &out {
            writer.push(r)?;
            features.push(MbrFeature {
                precursor: *p,
                run: *run,
                evidence_type: r.evidence_type,
                quantity: r.quantity,
                quantity_se: r.quantity_se,
                quantity_normalized: r.quantity_normalized,
                transfer_q: r.transfer_q,
                quant_q: r.quant_q,
                apex_rt_s: r.apex_rt_s,
                library_q: r.library_q,
                mass_err_ppm: *mass_err,
                claim: *claim,
            });
        }
        summary.rows_written += out.len();
        on_features(features)
    })?;
    let [at_q01, at_threshold] = at;
    v.at_q01 = at_q01;
    v.at_threshold = at_threshold;
    writer.finish()?;
    write_mbr_validation(validation_path, &v)?;
    summary.validation = v;
    summary.fdr = fdr.counts;
    summary.identified_fdr = id_fdr.counts;
    summary.features = all_features;
    Ok(summary)
}

fn unlimited() -> FinalizeLimits {
    FinalizeLimits {
        max_train_rows: usize::MAX,
        transfer_models: 1,
        max_norm_cells: usize::MAX,
        max_donor_fits: usize::MAX,
        max_held_scores: usize::MAX,
        temp_dir: None,
    }
}

#[allow(clippy::too_many_arguments)]
fn stream(
    lib: &Library,
    runs: &[Vec<ExtractRow>],
    limits: FinalizeLimits,
    batch_rows: usize,
    base: &Path,
    mbr_fdr: f64,
    quant_feature_fdr: Option<f64>,
) -> (FinalizeSummary, Vec<MbrFeature>) {
    let (mbr, val) = mbr_paths(base);
    let params = FinalizeParams {
        mbr_fdr,
        quant_feature_fdr,
        seed: 3,
        limits,
        mbr_path: &mbr,
        validation_path: &val,
        trace_fraction: 0.0,
    };
    let source = ExtractRuns { runs, batch_rows };
    let mut features = Vec::new();
    let s = finalize_stream(&stems(), lib, &source, &params, &mut |f| {
        features.extend(f);
        Ok(())
    })
    .unwrap();
    (s, features)
}

#[test]
fn map_reduce_matches_the_in_memory_reference() {
    let dir = tempfile::tempdir().unwrap();
    let spill_dir = tempfile::tempdir().unwrap();
    let lib = library_with_entries();
    assert_eq!(lib.entries.len(), N as usize);
    for (name, runs) in [
        ("plain", extract()),
        ("decoys", extract_with_identified_decoys()),
    ] {
        for (mbr_fdr, quant_feature_fdr) in [(0.01, None), (0.05, Some(0.01))] {
            let (m1, v1) = mbr_paths(&dir.path().join(format!("{name}-{mbr_fdr}-ref")));
            let want = reference(
                &stems(),
                &lib,
                &runs,
                mbr_fdr,
                quant_feature_fdr,
                3,
                &m1,
                &v1,
            )
            .unwrap();
            assert!(want.rows_written > 0 && want.validation.transfers_scored > 0);
            let arms = [
                ("one", usize::MAX, unlimited()),
                ("batches", 7, unlimited()),
                (
                    "groups",
                    97,
                    FinalizeLimits {
                        max_donor_fits: 100,
                        ..unlimited()
                    },
                ),
                (
                    "spill",
                    1000,
                    FinalizeLimits {
                        max_held_scores: 50,
                        temp_dir: Some(spill_dir.path().to_path_buf()),
                        ..unlimited()
                    },
                ),
            ];
            for (arm, batch_rows, limits) in arms {
                let base = dir.path().join(format!("{name}-{mbr_fdr}-{arm}"));
                let (s, features) = stream(
                    &lib,
                    &runs,
                    limits,
                    batch_rows,
                    &base,
                    mbr_fdr,
                    quant_feature_fdr,
                );
                let (m2, v2) = mbr_paths(&base);
                let what = format!("{name} {mbr_fdr} {arm}");
                assert_eq!(
                    std::fs::read(&m1).unwrap(),
                    std::fs::read(&m2).unwrap(),
                    "{what}"
                );
                assert_eq!(
                    std::fs::read(&v1).unwrap(),
                    std::fs::read(&v2).unwrap(),
                    "{what}"
                );
                assert_eq!(features, want.features, "{what}");
                assert_eq!(s.fdr, want.fdr, "{what}");
                assert_eq!(s.identified_fdr, want.identified_fdr, "{what}");
                assert_eq!(s.validation, want.validation, "{what}");
                assert_eq!(
                    (s.no_features, s.identified_no_features, s.rows_written),
                    (
                        want.no_features,
                        want.identified_no_features,
                        want.rows_written
                    ),
                    "{what}"
                );
            }
        }
    }
    for d in [dir.path(), spill_dir.path()] {
        let left: Vec<_> = std::fs::read_dir(d)
            .unwrap()
            .filter_map(|e| e.ok())
            .filter(|e| {
                e.file_name()
                    .to_string_lossy()
                    .starts_with(".andes-finalize-")
            })
            .collect();
        assert!(left.is_empty(), "temporary tables are removed");
    }
}

#[test]
fn sampled_training_reports_every_row() {
    let dir = tempfile::tempdir().unwrap();
    let runs = extract();
    let lib = library_with_entries();
    let (m1, v1) = mbr_paths(&dir.path().join("one"));
    let one = finalize(&stems(), &lib, &runs, 0.01, None, 3, &m1, &v1).unwrap();
    let limits = FinalizeLimits {
        max_train_rows: 600,
        max_norm_cells: 100,
        ..unlimited()
    };
    let (s, features) = stream(
        &lib,
        &runs,
        limits,
        1000,
        &dir.path().join("sampled"),
        0.01,
        None,
    );
    assert_eq!(s.rows_written, one.rows_written);
    assert_eq!(s.validation.transfers, one.validation.transfers);
    assert_eq!(s.validation.identified, one.validation.identified);
    let keys = |f: &[MbrFeature]| f.iter().map(|f| (f.precursor, f.run)).collect::<Vec<_>>();
    assert_eq!(keys(&features), keys(&one.features));
    assert!(s.validation.transfers_scored > 0);
    let accepted = features
        .iter()
        .filter(|f| f.evidence_type == EvidenceType::TransferMs1)
        .count();
    assert!(accepted > 0, "{accepted}");
}

#[test]
fn bagged_transfer_classifier_reports_every_row_and_is_reproducible() {
    let dir = tempfile::tempdir().unwrap();
    let runs = extract_with_identified_decoys();
    let lib = library_with_entries();
    let run = |name: &str, members: usize, batch_rows: usize| {
        let limits = FinalizeLimits {
            transfer_models: members,
            max_train_rows: 900,
            ..unlimited()
        };
        let base = dir.path().join(name);
        let (s, features) = stream(&lib, &runs, limits, batch_rows, &base, 0.05, None);
        let (m, v) = mbr_paths(&base);
        (
            s,
            features,
            std::fs::read(m).unwrap(),
            std::fs::read(v).unwrap(),
        )
    };
    let (one, one_f, _, _) = run("one", 1, 1000);
    let (a, a_f, a_m, a_v) = run("a", 3, 1000);
    let (_, _, b_m, b_v) = run("b", 3, 7);
    assert_eq!((a_m, a_v), (b_m, b_v));
    assert_eq!(a.rows_written, one.rows_written);
    assert_eq!(a.validation.transfers, one.validation.transfers);
    assert_eq!(
        a.validation.transfers_scored,
        one.validation.transfers_scored
    );
    let keys = |f: &[MbrFeature]| f.iter().map(|f| (f.precursor, f.run)).collect::<Vec<_>>();
    assert_eq!(keys(&a_f), keys(&one_f));
    assert!(a.fdr.failed_folds.is_empty(), "{:?}", a.fdr.failed_folds);
    let accepted = a_f
        .iter()
        .filter(|f| f.evidence_type == EvidenceType::TransferMs1)
        .count();
    assert!(accepted > 0, "{accepted}");
    assert!(a_f
        .iter()
        .zip(&one_f)
        .any(|(x, y)| x.transfer_q != y.transfer_q));
}

#[test]
fn finalize_rejects_rows_of_a_precursor_apart() {
    let dir = tempfile::tempdir().unwrap();
    let (mbr, val) = mbr_paths(&dir.path().join("x"));
    let mut runs = extract();
    let first = runs[1].remove(0);
    runs[1].push(first);
    let e = finalize(&stems(), &library(), &runs, 0.01, None, 1, &mbr, &val).unwrap_err();
    assert!(e.contains("not contiguous"), "{e}");
}

#[test]
fn finalize_rejects_out_of_range_rows_and_thresholds() {
    let dir = tempfile::tempdir().unwrap();
    let (mbr, val) = mbr_paths(&dir.path().join("x"));
    let mut runs = extract();
    runs[1].push(row(0, 2, Hypothesis::Target, Role::Transfer, None));
    let e = finalize(&stems(), &library(), &runs, 0.01, None, 1, &mbr, &val).unwrap_err();
    assert!(e.contains("run 2"), "{e}");
    runs[1].pop();
    runs[1].push(row(N, 1, Hypothesis::Target, Role::Transfer, None));
    let e = finalize(&stems(), &library(), &runs, 0.01, None, 1, &mbr, &val).unwrap_err();
    assert!(e.contains(&format!("precursor {N}")), "{e}");
    for bad in [0.0, -0.1, 1.5, f64::NAN] {
        assert!(finalize(&stems(), &library(), &extract(), bad, None, 1, &mbr, &val).is_err());
    }
    assert!(!mbr.exists());
}

fn jittered_rows(p: u32, jitter: [f64; 3]) -> Vec<ExtractRow> {
    let (a0, w0) = (1e6, 5.0);
    let scales = [1.0, 1.25, 2.0];
    let mut g = Rng(p as u64 + 1);
    (0..3u16)
        .flat_map(|run| {
            let scale = scales[run as usize];
            let amplitude = a0 / scale;
            let mut e = evidence(&mut g, true, 500.0 * scale);
            e.amplitude = amplitude;
            e.amplitude_se = 0.01 * amplitude;
            e.width_ratio = K_TRUE * jitter[run as usize];
            e.w = w0 * scale * e.width_ratio;
            e.quantity = amplitude * e.w * (2.0 * std::f64::consts::PI).sqrt();
            e.quantity_se = 0.01 * e.quantity;
            let identified = row(
                p,
                run,
                Hypothesis::Target,
                Role::Identified,
                Some(e.clone()),
            );
            let pseudo = row(p, run, Hypothesis::Target, Role::PseudoTransfer, Some(e));
            [identified, pseudo]
        })
        .collect()
}

const K_TRUE: f64 = 1.2;

#[test]
fn relative_width_gives_one_quantity_across_runs() {
    let truth = 1e6 * 5.0 * K_TRUE * (2.0 * std::f64::consts::PI).sqrt();
    for jitter in [[0.7, 1.4, 1.0], [1.4, 0.7, 1.4], [1.0, 1.0, 0.7]] {
        let rows = jittered_rows(3, jitter);
        let k = relative_widths(&rows);
        let per_run: Vec<(f64, f64)> = rows
            .iter()
            .filter(|r| r.role == Role::Identified)
            .map(|r| {
                let e = r.evidence.as_ref().unwrap();
                reported_quantity(e, k.get(&3).copied(), false, 0.0)
            })
            .collect();
        for &(q, se) in &per_run {
            assert!((q / per_run[0].0 - 1.0).abs() < 1e-12, "{per_run:?}");
            assert!((se / q - 0.01).abs() < 1e-12);
        }
        let mut sorted = jitter;
        sorted.sort_by(f64::total_cmp);
        let want = truth * sorted[1];
        assert!((per_run[0].0 / want - 1.0).abs() < 1e-12, "{per_run:?}");
    }
    let rows = jittered_rows(3, [0.7, 1.4, 1.0]);
    let raw = |run: usize| rows[2 * run].evidence.as_ref().unwrap().quantity;
    assert!((raw(1) / raw(0) - 2.0).abs() < 1e-9);
}

#[test]
fn relative_width_skips_unfitted_rows() {
    let mut rows = jittered_rows(4, [1.0, 0.7, 1.4]);
    rows[2].evidence.as_mut().unwrap().amplitude = 0.0;
    rows[4].evidence.as_mut().unwrap().width_ratio = f64::NAN;
    let k = relative_widths(&rows);
    assert!((k[&4] - K_TRUE).abs() < 1e-12, "{k:?}");
    let e = rows[4].evidence.as_ref().unwrap();
    assert_eq!(
        reported_quantity(e, Some(k[&4]), false, 0.0),
        (e.quantity, e.quantity_se)
    );
    let mut rows = jittered_rows(5, [1.0, 1.0, 1.0]);
    for r in rows.iter_mut().filter(|r| r.role == Role::Identified) {
        r.evidence = None;
    }
    assert!(relative_widths(&rows).is_empty());
    let e = rows[1].evidence.as_ref().unwrap();
    assert_eq!(
        reported_quantity(e, None, false, 0.0),
        (e.quantity, e.quantity_se)
    );
}

#[test]
fn width_weight_moves_the_quantity_to_the_own_fitted_width() {
    let rows = jittered_rows(3, [0.7, 1.4, 1.0]);
    let k = relative_widths(&rows)[&3];
    let e = rows[0].evidence.as_ref().unwrap();
    let at = |w: f64| reported_quantity(e, Some(k), false, w);
    assert!((at(1.0).0 / e.quantity - 1.0).abs() < 1e-12);
    assert!((at(1.0).1 / e.quantity_se - 1.0).abs() < 1e-12);
    let mid = (at(0.0).0 * at(1.0).0).sqrt();
    assert!((at(0.5).0 / mid - 1.0).abs() < 1e-12);
    let d = width_deviation(e, Some(k), false).unwrap();
    assert!((d - (0.7f64 / 1.0).log2()).abs() < 1e-12, "{d}");
    assert_eq!(width_deviation(e, None, false), None);
    let mut t = e.clone();
    t.trace_quantity = 3.0;
    assert_eq!(width_deviation(&t, Some(k), true), None);
    assert_eq!(reported_quantity(&t, Some(k), true, 1.0).0, 3.0);
}

fn width_cells(precursors: u32) -> Vec<WidthCell> {
    let mut g = Rng(11);
    let mut cells = Vec::new();
    for p in 0..precursors {
        let truth = 20.0 + 3.0 * g.normal();
        for run in 0..4u16 {
            let d = if run < 2 { 0.4 * g.normal() } else { 0.0 };
            let missed = if run == 0 { d } else { 0.0 };
            cells.push(WidthCell {
                precursor: p,
                run,
                log2_quantity: truth - missed + 0.05 * g.normal() + f64::from(run),
                deviation: d,
            });
        }
    }
    cells
}

#[test]
fn run_width_weights_keep_real_width_changes_and_drop_fit_noise() {
    let cells = width_cells(300);
    assert_eq!(run_width_weights(&cells, 4), vec![1.0, 0.0, 0.0, 0.0]);
    assert_eq!(run_width_weights(&cells, 5), vec![1.0, 0.0, 0.0, 0.0, 0.0]);
    let mut reversed = cells.clone();
    reversed.reverse();
    assert_eq!(run_width_weights(&reversed, 4), vec![1.0, 0.0, 0.0, 0.0]);
}

#[test]
fn run_width_weights_need_enough_cells() {
    let few = width_cells(WIDTH_WEIGHT_MIN_RUN_CELLS as u32 - 1);
    assert_eq!(run_width_weights(&few, 4), vec![0.0; 4]);
    let three: Vec<WidthCell> = width_cells(300).into_iter().filter(|c| c.run < 3).collect();
    assert_eq!(run_width_weights(&three, 4), vec![0.0; 4]);
    assert_eq!(run_width_weights(&[], 2), vec![0.0; 2]);
}

#[test]
fn only_identified_rows_take_the_run_width_weight() {
    let rows = jittered_rows(3, [0.7, 1.4, 1.0]);
    let weights = [0.6, 0.3, 1.0];
    for r in &rows {
        let want = if r.role == Role::Identified {
            weights[r.run as usize]
        } else {
            0.0
        };
        assert_eq!(row_width_weight(r, &weights), want);
    }
    assert_eq!(row_width_weight(&rows[0], &[]), 0.0);
    let mut transfer = rows[0].clone();
    transfer.role = Role::Transfer;
    assert_eq!(row_width_weight(&transfer, &weights), 0.0);
}

#[test]
fn weigh_widths_moves_only_cells_with_a_deviation() {
    let cells: Vec<(NormCell, Option<f64>)> = width_cells(300)
        .into_iter()
        .enumerate()
        .map(|(i, c)| {
            let cell = NormCell {
                precursor: c.precursor,
                run: c.run,
                log2_quantity: c.log2_quantity,
                rt: 100.0,
            };
            (cell, (i % 10 != 0).then_some(c.deviation))
        })
        .collect();
    let (out, weights) = weigh_widths(cells.clone(), 4);
    assert_eq!(weights, vec![1.0, 0.0, 0.0, 0.0]);
    for ((before, d), after) in cells.iter().zip(&out) {
        let want = before.log2_quantity + d.map_or(0.0, |d| weights[before.run as usize] * d);
        assert_eq!(after.log2_quantity, want);
        assert_eq!(
            (after.precursor, after.run, after.rt),
            (before.precursor, before.run, before.rt)
        );
    }
    assert_eq!(describe_width_weights(&weights), "0.0: 3, 1.0: 1");
}

#[test]
fn trace_fraction_reports_the_trace_area_of_the_brightest_precursors() {
    use crate::quant::finalize::{ExtractFiles, ExtractSource};
    let lib = library();
    let mut runs = extract();
    for r in runs[0].iter_mut().filter(|r| r.role == Role::Identified) {
        let e = r.evidence.as_mut().unwrap();
        e.amplitude_se = e.amplitude / (10.0 + r.precursor as f64);
        if r.precursor == 0 {
            e.trace_quantity = 0.0;
        }
    }
    let dir = tempfile::tempdir().unwrap();
    let files = ExtractFiles(
        runs.iter()
            .enumerate()
            .map(|(r, rows)| {
                let path = dir.path().join(format!("{r}.extract.parquet"));
                output::extract::write_extract_table(&path, rows).unwrap();
                path
            })
            .collect(),
    );
    let memory = ExtractRuns {
        runs: &runs,
        batch_rows: usize::MAX,
    };
    let sources: [(&str, &dyn ExtractSource); 2] = [("memory", &memory), ("files", &files)];
    let target: HashMap<(u32, u16), &Evidence> = runs
        .iter()
        .flatten()
        .filter(|r| r.kind == Hypothesis::Target && r.role != Role::PseudoTransfer)
        .filter_map(|r| r.evidence.as_ref().map(|e| ((r.precursor, r.run), e)))
        .collect();
    for (name, source) in sources {
        let finalize_at = |fraction: f64| {
            let (mbr, val) = mbr_paths(&dir.path().join(format!("{name}{fraction}")));
            let params = FinalizeParams {
                mbr_fdr: 0.01,
                quant_feature_fdr: None,
                seed: 7,
                limits: FinalizeLimits::default(),
                mbr_path: &mbr,
                validation_path: &val,
                trace_fraction: fraction,
            };
            let mut features = Vec::new();
            finalize_stream(&stems(), &lib, source, &params, &mut |f| {
                features.extend(f);
                Ok(())
            })
            .unwrap();
            features
        };
        let off = finalize_at(0.0);
        let all = finalize_at(1.0);
        let half = finalize_at(0.5);
        assert_eq!(off.len(), all.len(), "{name}");
        assert_eq!(off.len(), half.len(), "{name}");
        let mut traced = 0;
        for ((o, a), h) in off.iter().zip(&all).zip(&half) {
            let key = (o.precursor, o.run);
            assert_eq!(key, (a.precursor, a.run));
            assert_eq!(key, (h.precursor, h.run));
            for f in [a, h] {
                assert_eq!(
                    (f.evidence_type, f.transfer_q, f.quant_q),
                    (o.evidence_type, o.transfer_q, o.quant_q),
                    "{name}"
                );
            }
            let Some(q) = o.quantity else {
                assert!(a.quantity.is_none() && h.quantity.is_none());
                continue;
            };
            let e = target[&key];
            assert_ne!(q, e.trace_quantity);
            if e.trace_quantity <= 0.0 {
                for f in [a, h] {
                    assert_eq!(
                        (f.quantity, f.quantity_se),
                        (o.quantity, o.quantity_se),
                        "{name} {key:?}"
                    );
                }
                continue;
            }
            assert_eq!(
                (a.quantity, a.quantity_se),
                (Some(e.trace_quantity), Some(e.trace_quantity_se)),
                "{name} {key:?}"
            );
            if o.precursor >= N / 2 {
                assert_eq!(h.quantity, Some(e.trace_quantity), "{name} {key:?}");
                traced += 1;
            } else {
                assert_eq!((h.quantity, h.quantity_se), (o.quantity, o.quantity_se));
            }
        }
        assert!(traced > 0, "{name}");
    }
}

#[test]
fn finalize_reports_the_quantity_at_the_relative_width() {
    let lib = library();
    let mut runs = extract();
    let mut g = Rng(9);
    for (run, amplitude, w, ratio) in [(0u16, 1e6, 4.0, 0.8), (1, 5e5, 12.0, 1.2)] {
        let mut e = evidence(&mut g, true, 100.0);
        e.amplitude = amplitude;
        e.amplitude_se = 1e3;
        e.w = w;
        e.width_ratio = ratio;
        e.quantity = e.amplitude * w * (2.0 * std::f64::consts::PI).sqrt();
        let rows = &mut runs[run as usize];
        rows.retain(|r| !(r.precursor == 0 && r.role != Role::PseudoTransfer));
        let at = rows
            .iter()
            .position(|r| r.precursor != 0)
            .unwrap_or(rows.len());
        rows.insert(
            at,
            row(0, run, Hypothesis::Target, Role::Identified, Some(e)),
        );
    }
    let dir = tempfile::tempdir().unwrap();
    let (mbr, val) = mbr_paths(&dir.path().join("x"));
    let s = finalize(&stems(), &lib, &runs, 0.01, None, 7, &mbr, &val).unwrap();
    let mut p0 = s.features.iter().filter(|f| f.precursor == 0);
    let sqrt_2pi = (2.0 * std::f64::consts::PI).sqrt();
    let a = p0.next().unwrap();
    let b = p0.next().unwrap();
    assert_eq!((a.run, b.run), (0, 1));
    assert!((a.quantity.unwrap() / (1e6 * 5.0 * sqrt_2pi) - 1.0).abs() < 1e-12);
    assert!((b.quantity.unwrap() / (5e5 * 10.0 * sqrt_2pi) - 1.0).abs() < 1e-12);
    assert!((a.quantity_se.unwrap() / (1e3 * 5.0 * sqrt_2pi) - 1.0).abs() < 1e-12);
    let id0 = runs[0]
        .iter()
        .find(|r| r.precursor == 0 && r.role == Role::Identified);
    assert!(
        (id0.unwrap().evidence.as_ref().unwrap().quantity / (4e6 * sqrt_2pi) - 1.0).abs()
            < 1e-12
    );
}

fn foreign(g: &mut Rng, t0: f64) -> Evidence {
    let mut e = evidence(g, true, t0);
    e.llr = 50.0 + 15.0 * g.normal();
    e.envelope_cos = 0.7 + 0.1 * g.normal();
    e.coherence = 0.4 + 0.15 * g.normal();
    e.rt_z = 2.0 * g.normal();
    e.m1_ratio = 0.5 + 0.2 * g.unit();
    e.isotopes_used = 0b0111;
    e.mass_err_ppm = 5.0 * g.normal();
    e
}

fn identified_kind(p: u32) -> &'static str {
    match p % 10 {
        0 => "wrong",
        5 => "interfered",
        _ => "clean",
    }
}

fn extract_with_identified_decoys() -> Vec<Vec<ExtractRow>> {
    let mut runs = extract();
    let mut g = Rng(11);
    let mut run0 = Vec::new();
    for r in runs[0].drain(..) {
        if r.role != Role::Identified {
            run0.push(r);
            continue;
        }
        let (p, apex) = (r.precursor, r.evidence.as_ref().unwrap().t0);
        let target = match identified_kind(p) {
            "wrong" => evidence(&mut g, false, apex),
            "interfered" => foreign(&mut g, apex),
            _ => evidence(&mut g, true, apex),
        };
        run0.push(row(
            p,
            0,
            Hypothesis::Target,
            Role::Identified,
            Some(target),
        ));
        for (k, kind) in [Hypothesis::MassShiftDecoy(-7), Hypothesis::RtShiftDecoy]
            .into_iter()
            .enumerate()
        {
            let t0 = apex + 300.0 * k as f64;
            let e = if (p + k as u32) % 3 == 0 {
                foreign(&mut g, t0)
            } else {
                evidence(&mut g, false, t0)
            };
            run0.push(row(p, 0, kind, Role::Identified, Some(e)));
        }
    }
    runs[0] = run0;
    runs
}

#[test]
fn identified_quality_q_separates_clean_from_wrong_and_interfered_peaks() {
    let dir = tempfile::tempdir().unwrap();
    let (mbr, val) = mbr_paths(&dir.path().join("x"));
    let runs = extract_with_identified_decoys();
    let s = finalize(&stems(), &library(), &runs, 0.01, None, 7, &mbr, &val).unwrap();
    let mut q: HashMap<&str, Vec<f64>> = HashMap::new();
    for f in s.features.iter().filter(|f| f.run == 0) {
        assert_eq!(f.evidence_type, EvidenceType::Identified);
        assert_eq!(f.transfer_q, None);
        q.entry(identified_kind(f.precursor))
            .or_default()
            .push(f.quant_q.unwrap());
    }
    let share = |kind: &str, ok: &dyn Fn(f64) -> bool| {
        let v = &q[kind];
        v.iter().filter(|&&x| ok(x)).count() as f64 / v.len() as f64
    };
    let clean = share("clean", &|x| x <= 0.01);
    let wrong = share("wrong", &|x| x > 0.05);
    let interfered = share("interfered", &|x| x > 0.05);
    assert!(clean >= 0.95, "{clean}");
    assert!(wrong >= 0.95, "{wrong}");
    assert!(interfered >= 0.8, "{interfered}");
    let v = &s.validation;
    assert_eq!(
        (v.identified, v.identified_scored),
        (N as usize, N as usize)
    );
    assert_eq!(v.identified_below_threshold, 0);
    assert_eq!(s.identified_no_features, 0);
    for f in s.features.iter().filter(|f| f.run == 1) {
        if f.evidence_type == EvidenceType::BelowThreshold {
            assert_eq!(f.quant_q, f.transfer_q);
        } else {
            assert!(f.quant_q.is_some_and(|q| (0.0..=1.0).contains(&q)));
            assert!(f.quant_q >= f.transfer_q, "{f:?}");
        }
    }
    let rows = output::mbr::read_mbr_parquet(&mbr).unwrap();
    let accepted = |t: &str| matches!(t, "transfer_ms2" | "transfer_ms1" | "shared_feature");
    assert!(rows.iter().any(|r| accepted(r.evidence_type.as_str())));
    for r in &rows {
        if accepted(r.evidence_type.as_str()) {
            let (t, m) = (r.transfer_q.unwrap(), r.measurement_q.unwrap());
            assert!((0.0..=1.0).contains(&m), "{r:?}");
            assert_eq!(r.quant_q, Some(t.max(m)), "{r:?}");
        } else {
            assert_eq!(r.measurement_q, None, "{r:?}");
        }
    }
}

fn extract_two_identified_runs(gain: f64) -> Vec<Vec<ExtractRow>> {
    let mut runs = extract_with_identified_decoys();
    let copies: Vec<ExtractRow> = runs[0]
        .iter()
        .filter(|r| r.precursor % 2 == 0 && r.role == Role::Identified)
        .map(|r| {
            let mut r = r.clone();
            r.run = 1;
            if let Some(e) = r.evidence.as_mut() {
                e.amplitude *= gain;
                e.amplitude_se *= gain;
                e.quantity *= gain;
                e.quantity_se *= gain;
            }
            r
        })
        .collect();
    let run1: Vec<ExtractRow> = runs[1]
        .drain(..)
        .filter(|r| r.precursor % 2 == 1)
        .chain(copies)
        .collect();
    runs[1] = run1;
    runs
}

#[test]
fn finalize_reports_run_normalized_quantities() {
    let dir = tempfile::tempdir().unwrap();
    let (mbr, val) = mbr_paths(&dir.path().join("x"));
    let runs = extract_two_identified_runs(2.0);
    let s = finalize(&stems(), &library(), &runs, 0.01, None, 7, &mbr, &val).unwrap();
    let get = |p: u32, run: u16| {
        s.features
            .iter()
            .find(|f| f.precursor == p && f.run == run)
            .unwrap()
    };
    let mut clean = 0;
    for p in (0..N).filter(|p| p % 2 == 0 && identified_kind(*p) == "clean") {
        let (a, b) = (get(p, 0), get(p, 1));
        assert!((b.quantity.unwrap() / a.quantity.unwrap() - 2.0).abs() < 1e-9);
        let (na, nb) = (
            a.quantity_normalized.unwrap(),
            b.quantity_normalized.unwrap(),
        );
        assert!((nb / na - 1.0).abs() < 1e-9, "{p}: {na} {nb}");
        assert!((na / (a.quantity.unwrap() * 2f64.sqrt()) - 1.0).abs() < 1e-9);
        clean += 1;
    }
    assert!(clean > 300, "{clean}");
    for f in s.features.iter().filter(|f| f.quantity.is_some()) {
        let want = f.quantity.unwrap()
            * if f.run == 0 {
                2f64.sqrt()
            } else {
                0.5f64.sqrt()
            };
        let got = f.quantity_normalized.unwrap();
        assert!((got / want - 1.0).abs() < 1e-9, "{f:?}");
    }
    let t = table(&mbr);
    let col = |name: &str| t[0].iter().position(|h| h == name).unwrap();
    assert_eq!(col("quantity_normalized"), col("quantity_se") + 1);
}

#[test]
fn whole_pattern_divides_quantity_and_error_by_the_share() {
    assert_eq!(whole_pattern((3.0, 0.5), 1.0), (3.0, 0.5));
    assert_eq!(whole_pattern((3.0, 0.5), 0.5), (6.0, 1.0));
    let share = 0.871_108_8_f32;
    let (q, se) = whole_pattern((1e7, 2e5), share);
    assert_eq!((q, se), (1e7 / f64::from(share), 2e5 / f64::from(share)));
    assert!((q / 1e7 - 1.0 / 0.871_108_8).abs() < 1e-7);
}

#[test]
fn finalize_reports_whole_pattern_quantities() {
    let share = |p: u32| [1.0_f32, 0.948_998_3, 0.871_108_8, 0.5][p as usize % 4];
    let mut lib = library();
    for (p, m) in lib.meta.iter_mut().enumerate() {
        m.isotope_share = share(p as u32);
    }
    let dir = tempfile::tempdir().unwrap();
    let runs = extract_two_identified_runs(2.0);
    let run = |name: &str, lib: &Library| {
        let (mbr, val) = mbr_paths(&dir.path().join(name));
        let s = finalize(&stems(), lib, &runs, 0.01, None, 7, &mbr, &val).unwrap();
        (s.features, table(&mbr), table(&val))
    };
    let (four, four_table, four_validation) = run("four", &library());
    let (whole, whole_table, whole_validation) = run("whole", &lib);
    assert_eq!(four_validation, whole_validation);
    assert_eq!(four.len(), whole.len());
    let (mut scaled, mut normalized) = (0, 0);
    for (f, w) in four.iter().zip(&whole) {
        assert_eq!((f.precursor, f.run), (w.precursor, w.run));
        assert_eq!(
            (f.evidence_type, f.transfer_q, f.quant_q, f.apex_rt_s),
            (w.evidence_type, w.transfer_q, w.quant_q, w.apex_rt_s)
        );
        let s = f64::from(share(f.precursor));
        assert_eq!(w.quantity, f.quantity.map(|q| q / s), "{f:?}");
        assert_eq!(w.quantity_se, f.quantity_se.map(|q| q / s), "{f:?}");
        assert_eq!(
            w.quantity_normalized.is_some(),
            f.quantity_normalized.is_some()
        );
        if let (Some(a), Some(b)) = (f.quantity_normalized, w.quantity_normalized) {
            assert!((b * s / a - 1.0).abs() < 1e-12, "{f:?}");
            normalized += 1;
        }
        if f.precursor % 4 == 0 {
            assert_eq!(w.quantity_normalized, f.quantity_normalized);
        } else if f.quantity.is_some() {
            scaled += 1;
        }
    }
    assert!(scaled > 300 && normalized > 300, "{scaled} {normalized}");
    assert_eq!(four_table.len(), whole_table.len());
    let quantity_columns = ["quantity", "quantity_se", "quantity_normalized"];
    for (a, b) in four_table.iter().zip(&whole_table) {
        for ((name, a), b) in four_table[0].iter().zip(a).zip(b) {
            if !quantity_columns.contains(&name.as_str()) {
                assert_eq!(a, b, "{name}");
            }
        }
    }
    let quantity = four_table[0].iter().position(|h| h == "quantity").unwrap();
    assert!(four_table
        .iter()
        .zip(&whole_table)
        .skip(1)
        .any(|(a, b)| a[quantity] != b[quantity]));
}

#[test]
fn finalize_leaves_quantities_unnormalized_without_two_reference_runs() {
    let dir = tempfile::tempdir().unwrap();
    let (mbr, val) = mbr_paths(&dir.path().join("x"));
    let runs = extract_with_identified_decoys();
    let s = finalize(&stems(), &library(), &runs, 0.01, None, 7, &mbr, &val).unwrap();
    assert!(s.features.iter().all(|f| f.quantity_normalized.is_none()));
    let t = table(&mbr);
    let col = t[0]
        .iter()
        .position(|h| h == "quantity_normalized")
        .unwrap();
    assert!(t[1..].iter().all(|r| r[col] == "NA"));
}

#[test]
fn identified_quality_leaves_the_transfer_q_values_unchanged() {
    let dir = tempfile::tempdir().unwrap();
    let (m1, v1) = mbr_paths(&dir.path().join("plain"));
    let (m2, v2) = mbr_paths(&dir.path().join("decoys"));
    let plain = finalize(&stems(), &library(), &extract(), 0.01, None, 7, &m1, &v1).unwrap();
    let runs = extract_with_identified_decoys();
    let with = finalize(&stems(), &library(), &runs, 0.01, None, 7, &m2, &v2).unwrap();
    let transfers = |s: &FinalizeSummary| -> Vec<(u32, Option<f64>, EvidenceType)> {
        s.features
            .iter()
            .filter(|f| f.run == 1)
            .map(|f| (f.precursor, f.transfer_q, f.evidence_type))
            .collect()
    };
    assert_eq!(transfers(&plain), transfers(&with));
    assert_eq!(plain.validation.at_q01, with.validation.at_q01);
    assert!(plain
        .features
        .iter()
        .filter(|f| f.run == 0)
        .all(|f| f.quant_q.is_none() && f.evidence_type == EvidenceType::Identified));
}

#[test]
fn quant_feature_fdr_flags_identified_rows_above_it() {
    let dir = tempfile::tempdir().unwrap();
    let (mbr, val) = mbr_paths(&dir.path().join("x"));
    let runs = extract_with_identified_decoys();
    let s = finalize(&stems(), &library(), &runs, 0.01, Some(0.01), 7, &mbr, &val).unwrap();
    let (mut flagged, mut kept) = (0, 0);
    for f in s.features.iter().filter(|f| f.run == 0) {
        let below = f.quant_q.is_none_or(|q| q > 0.01);
        let want = if below {
            EvidenceType::IdentifiedBelowThreshold
        } else {
            EvidenceType::Identified
        };
        assert_eq!(f.evidence_type, want, "{f:?}");
        assert!(f.quantity.is_some());
        if below {
            flagged += 1;
        } else {
            kept += 1;
        }
    }
    assert_eq!(s.validation.identified_below_threshold, flagged);
    assert!(
        flagged >= 160 && kept >= 680,
        "{flagged} flagged, {kept} kept"
    );
    let t = table(&mbr);
    assert_eq!(t[0][7], "quant_q");
    let flagged_rows = t
        .iter()
        .filter(|r| r[2] == "identified_below_threshold")
        .count();
    assert_eq!(flagged_rows, flagged);
    let t = table(&val);
    let get = |name: &str| &t[1][t[0].iter().position(|h| h == name).unwrap()];
    assert_eq!(get("quant_feature_threshold"), "0.01");
    assert_eq!(get("identified_below_threshold"), &flagged.to_string());
    for bad in [0.0, 1.5, f64::NAN] {
        assert!(finalize(&stems(), &library(), &runs, 0.01, Some(bad), 7, &mbr, &val).is_err());
    }
}

#[test]
fn transfer_features_follow_the_feature_names() {
    let mut g = Rng(1);
    let mut r = row(
        0,
        1,
        Hypothesis::Target,
        Role::Transfer,
        Some(evidence(&mut g, true, 50.0)),
    );
    r.n_donor_runs = 3;
    let e = r.evidence.clone().unwrap();
    let x = transfer_features(&r, &DonorShapes::default(), &[], 0.0, 1e-4).unwrap();
    assert_eq!(x.len(), TRANSFER_FEATURES.len());
    let at = |name: &str| x[TRANSFER_FEATURES.iter().position(|n| *n == name).unwrap()];
    assert_eq!(at("abs_rt_z"), e.rt_z.abs() as f32);
    assert_eq!(at("n_isotopes_used"), 4.0);
    assert_eq!(at("n_donor_runs"), 3.0);
    assert_eq!(
        at("log10_amplitude_snr"),
        (1e6f64 / 1e4 + 1.0).log10() as f32
    );
    for se in [0.0, -1.0, f64::NAN] {
        r.evidence.as_mut().unwrap().amplitude_se = se;
        assert!(transfer_features(&r, &DonorShapes::default(), &[], 0.0, 1e-4).is_none());
        assert!(identified_features(&r).is_none());
        assert!(amplitude_snr(r.evidence.as_ref().unwrap()).is_nan());
    }
    r.evidence.as_mut().unwrap().amplitude_se = f64::INFINITY;
    let x = transfer_features(&r, &DonorShapes::default(), &[], 0.0, 1e-4).unwrap();
    let snr = TRANSFER_FEATURES
        .iter()
        .position(|n| *n == "log10_amplitude_snr")
        .unwrap();
    assert_eq!(x[snr], 0.0);
    assert!(identified_features(&r).is_some());
    r.evidence.as_mut().unwrap().amplitude_se = 1e4;
    r.evidence.as_mut().unwrap().llr = f64::NAN;
    assert!(transfer_features(&r, &DonorShapes::default(), &[], 0.0, 1e-4).is_none());
    r.evidence = None;
    assert!(transfer_features(&r, &DonorShapes::default(), &[], 0.0, 1e-4).is_none());
}

#[test]
fn a_wrong_measurement_is_a_wrong_peak_or_a_wrong_fitted_quantity() {
    let mut g = Rng(3);
    let id = evidence(&mut g, true, 50.0);
    let mut p = id.clone();
    p.w = 2.5 * id.w;
    p.width_ratio = 2.5 * id.width_ratio;
    p.quantity = 2.5 * id.quantity;
    assert_eq!(
        mismatch(&p, &id),
        Some(Mismatch {
            apex: false,
            quantity: Some(true)
        })
    );
    assert_eq!(wrong_label(&p, &id), Some(true));
    p.quantity = 1.9 * id.quantity;
    p.trace_quantity = 3.0 * id.trace_quantity;
    assert!(!mismatch(&p, &id).unwrap().wrong());
    assert_eq!(wrong_label(&p, &id), Some(false));
    p.t0 = id.t0 + id.w;
    assert!(!mismatch(&p, &id).unwrap().apex);
    p.t0 = id.t0 + id.w + 0.1;
    let m = mismatch(&p, &id).unwrap();
    assert!(m.apex && m.quantity == Some(false) && m.wrong());
    let no_quantity = |e: &Evidence| {
        let mut e = e.clone();
        e.amplitude = 0.0;
        e.quantity = 0.0;
        e
    };
    for (a, b) in [(no_quantity(&p), id.clone()), (p.clone(), no_quantity(&id))] {
        let m = mismatch(&a, &b).unwrap();
        assert!(m.apex && m.quantity.is_none() && m.wrong());
        assert_eq!(wrong_label(&a, &b), None);
        let mut a = a;
        a.t0 = b.t0;
        let m = mismatch(&a, &b).unwrap();
        assert!(!m.apex && m.quantity.is_none() && !m.wrong());
    }
    let (mut a, mut b) = (p.clone(), id.clone());
    a.t0 = f64::NAN;
    b.t0 = f64::NAN;
    assert_eq!(mismatch(&a, &id), None);
    assert_eq!(mismatch(&p, &b), None);
    let mut b = id.clone();
    b.w = f64::NAN;
    assert_eq!(mismatch(&p, &b), None);
}

#[test]
fn a_pseudo_transfer_copied_into_its_identified_row_has_no_label() {
    let mut g = Rng(21);
    let e = evidence(&mut g, true, 500.0);
    assert_eq!(mismatch(&e, &e), None);
    assert_eq!(wrong_label(&e, &e), None);
    let mut n = e.clone();
    n.m1_ratio = f64::NAN;
    assert_eq!(mismatch(&n, &n.clone()), None);
    let mut id = e.clone();
    id.t0 += 3.0 * id.w;
    assert!(mismatch(&e, &id).unwrap().apex);
}

#[test]
fn wrong_features_take_the_fitted_quantity_without_floors() {
    let mut g = Rng(5);
    let mut e = evidence(&mut g, true, 50.0);
    e.quantity = 2e7;
    e.quantity_se = 4e5;
    let x = wrong_features(&e).unwrap();
    assert_eq!(x.len(), 11);
    assert_eq!(x[2], e.width_ratio.ln().abs() as f32);
    assert_eq!(x[8], (1e6f64 / 1e4 + 1.0).log10() as f32);
    assert_eq!(x[9], (4e5f64 / 2e7).ln() as f32);
    assert_eq!(x[10], 2e7f64.log2() as f32);
    let unusable: [fn(&mut Evidence); 9] = [
        |e| e.quantity = 0.0,
        |e| e.quantity = -1.0,
        |e| e.quantity = f64::NAN,
        |e| e.quantity_se = 0.0,
        |e| e.quantity_se = -1.0,
        |e| e.quantity_se = f64::NAN,
        |e| e.amplitude_se = 0.0,
        |e| e.amplitude_se = f64::NAN,
        |e| e.width_ratio = 0.0,
    ];
    for (i, f) in unusable.iter().enumerate() {
        let mut e = e.clone();
        f(&mut e);
        assert!(wrong_features(&e).is_none(), "{i}");
    }
}

#[test]
fn an_accepted_transfer_has_the_larger_of_its_two_q_values() {
    assert_eq!(transfer_quant_q(Some(0.004), 0.002), 0.004);
    assert_eq!(transfer_quant_q(Some(0.004), 0.02), 0.02);
    assert_eq!(transfer_quant_q(None, 0.3), 0.3);
}

use quant::ms2::Ms2Evidence;

fn supporting() -> Ms2Evidence {
    Ms2Evidence {
        n_scans: 2,
        best: quant::ms2::ScanMatch {
            n_matched: 9,
            angle: 0.9,
            explained: 0.95,
        },
    }
}

fn empty_ms2() -> Ms2Evidence {
    Ms2Evidence {
        n_scans: 1,
        best: quant::ms2::ScanMatch::default(),
    }
}

fn identity_maps(n: usize) -> Vec<Option<RunMap>> {
    (0..n)
        .map(|_| {
            Some(RunMap::from_parts(vec![0.0, 1000.0], vec![0.0, 1000.0], 1.0, 1.0).unwrap())
        })
        .collect()
}

fn shaped(w: f64, tail: f64, iso: [f64; 4]) -> Evidence {
    let mut e = evidence(&mut Rng(9), true, 100.0);
    e.w = w;
    e.tailing = tail;
    e.iso_heights = iso;
    e
}

const DONOR_ISO: [f64; 4] = [5.5e5, 3.0e5, 1.1e5, 0.4e5];

#[test]
fn shape_features_compare_with_the_donor_fits_of_other_runs() {
    let maps = identity_maps(3);
    let id = |run: u16, w: f64, tail: f64, iso| {
        row(
            7,
            run,
            Hypothesis::Target,
            Role::Identified,
            Some(shaped(w, tail, iso)),
        )
    };
    let donors = [id(0, 4.0, 0.10, DONOR_ISO), id(2, 4.0, 0.10, DONOR_ISO)];
    let shapes = DonorShapes::from_rows(&donors, &maps);
    let like = row(
        7,
        1,
        Hypothesis::Target,
        Role::Transfer,
        Some(shaped(4.0, 0.10, DONOR_ISO.map(|v| 2.0 * v))),
    );
    let f = shapes.features(&like, &maps);
    assert!(
        f[0].abs() < 1e-12 && f[1].abs() < 1e-12 && (f[2] - 1.0).abs() < 1e-12,
        "{f:?}"
    );
    let other = row(
        7,
        1,
        Hypothesis::MassShiftDecoy(7),
        Role::Transfer,
        Some(shaped(8.0, -0.3, [2.0e5, 4.0e5, 3.0e5, 1.0e5])),
    );
    let g = shapes.features(&other, &maps);
    assert!((g[0] - 2f64.ln()).abs() < 1e-12, "{g:?}");
    assert!((g[1] + 0.4).abs() < 1e-12, "{g:?}");
    assert!(g[2] < f[2] - 0.1, "{g:?}");
    let x = transfer_features(&other, &shapes, &maps, 0.0, 1e-4).unwrap();
    let n = TRANSFER_FEATURES.len();
    assert_eq!(x.len(), n);
    assert_eq!(
        &TRANSFER_FEATURES[n - 4..],
        ["shape_log_width", "shape_tailing_diff", "shape_iso_cos", "log10_donor_pep"]
    );
    assert_eq!(&x[n - 4..n - 1], &g.map(|v| v as f32));
    assert_eq!(x[n - 1], -4.0);
}

#[test]
fn shape_features_leave_out_the_rows_own_run() {
    let maps = identity_maps(2);
    let donors = [
        row(
            3,
            0,
            Hypothesis::Target,
            Role::Identified,
            Some(shaped(4.0, 0.0, DONOR_ISO)),
        ),
        row(
            3,
            1,
            Hypothesis::Target,
            Role::Identified,
            Some(shaped(9.0, 0.5, [1.0, 1.0, 1.0, 1.0])),
        ),
    ];
    let shapes = DonorShapes::from_rows(&donors, &maps);
    let pseudo = row(
        3,
        1,
        Hypothesis::Target,
        Role::PseudoTransfer,
        Some(shaped(9.0, 0.5, [1.0, 1.0, 1.0, 1.0])),
    );
    let f = shapes.features(&pseudo, &maps);
    assert!((f[0] - (9.0f64 / 4.0).ln()).abs() < 1e-12, "{f:?}");
    assert!((f[1] - 0.5).abs() < 1e-12, "{f:?}");
    assert!(f[2] < 0.99, "{f:?}");
    let only = DonorShapes::from_rows(&donors[1..], &maps);
    assert_eq!(only.features(&pseudo, &maps), NO_SHAPE);
    assert_eq!(shapes.features(&pseudo, &[None, None]), NO_SHAPE);
    let mut empty = pseudo.clone();
    empty.evidence.as_mut().unwrap().amplitude = 0.0;
    assert_eq!(shapes.features(&empty, &maps), NO_SHAPE);
    let decoy = row(
        3,
        0,
        Hypothesis::RtShiftDecoy,
        Role::Identified,
        Some(shaped(1.0, 0.0, DONOR_ISO)),
    );
    let with_decoy = DonorShapes::from_rows(donors.iter().chain([&decoy]), &maps);
    assert_eq!(with_decoy.features(&pseudo, &maps), f);
}

#[test]
fn transfer_features_leave_out_the_ms2_evidence() {
    let mut g = Rng(3);
    let mut r = row(
        0,
        1,
        Hypothesis::Target,
        Role::Transfer,
        Some(evidence(&mut g, true, 100.0)),
    );
    let x = transfer_features(&r, &DonorShapes::default(), &[], 0.0, 1e-4).unwrap();
    r.ms2 = Some(supporting());
    assert_eq!(
        transfer_features(&r, &DonorShapes::default(), &[], 0.0, 1e-4).unwrap(),
        x
    );
    assert!(TRANSFER_FEATURES.iter().all(|f| !f.starts_with("ms2")));
}

#[test]
fn supported_accepted_transfers_are_transfer_ms2() {
    let mut runs = extract();
    for r in runs[1].iter_mut() {
        r.ms2 = Some(if r.kind == Hypothesis::Target && r.precursor % 4 == 0 {
            supporting()
        } else {
            empty_ms2()
        });
    }
    let dir = tempfile::tempdir().unwrap();
    let (mbr, val) = mbr_paths(&dir.path().join("x"));
    let s = finalize(&stems(), &library(), &runs, 0.01, None, 7, &mbr, &val).unwrap();
    let t = table(&mbr);
    let ms2_rows: Vec<&Vec<String>> =
        t[1..].iter().filter(|r| r[2] == "transfer_ms2").collect();
    assert!(!ms2_rows.is_empty());
    for r in &ms2_rows {
        let p: u32 = r[0]
            .trim_start_matches("PEP")
            .trim_end_matches("/2")
            .parse()
            .unwrap();
        assert_eq!(p % 4, 0, "{r:?}");
    }
    let v = &s.validation;
    assert_eq!(v.transfer_ms2, ms2_rows.len());
    assert_eq!(v.ms2_rows[0], N as usize);
    assert_eq!(v.ms2_supported[0], (N as usize).div_ceil(4));
    assert_eq!((v.ms2_supported[1], v.ms2_supported[2]), (0, 0));
    assert!(v.ms2_rows[1] > 0 && v.ms2_rows[2] > 0);
    let vt = table(&val);
    let get = |name: &str| &vt[1][vt[0].iter().position(|h| h == name).unwrap()];
    assert_eq!(get("transfer_ms2"), &ms2_rows.len().to_string());
}

fn ambiguous_ms1(with_ms2: bool) -> Vec<Vec<ExtractRow>> {
    let mut runs = extract();
    let mut g = Rng(11);
    let mut run1 = Vec::new();
    for p in 0..N {
        let apex = 100.0 + p as f64;
        for kind in [
            Hypothesis::Target,
            Hypothesis::MassShiftDecoy(7),
            Hypothesis::RtShiftDecoy,
        ] {
            let mut r = row(
                p,
                1,
                kind,
                Role::Transfer,
                Some(evidence(&mut g, true, apex)),
            );
            if with_ms2 {
                r.ms2 = Some(if kind == Hypothesis::Target && p % 2 == 0 {
                    supporting()
                } else {
                    empty_ms2()
                });
            }
            run1.push(r);
        }
    }
    runs[1] = run1;
    runs
}

#[test]
fn ms2_evidence_does_not_change_the_transfer_q_values() {
    let accepted = |with_ms2: bool| {
        let dir = tempfile::tempdir().unwrap();
        let (mbr, val) = mbr_paths(&dir.path().join("x"));
        let runs = ambiguous_ms1(with_ms2);
        finalize(&stems(), &library(), &runs, 0.01, None, 7, &mbr, &val).unwrap();
        let (mut present, mut absent) = (0usize, 0usize);
        for r in &table(&mbr)[1..] {
            if r[1] != "b" || !r[2].starts_with("transfer_ms") {
                continue;
            }
            let p: u32 = r[0]
                .trim_start_matches("PEP")
                .trim_end_matches("/2")
                .parse()
                .unwrap();
            if p % 2 == 0 {
                present += 1;
            } else {
                absent += 1;
            }
        }
        (present, absent)
    };
    let (p0, a0) = accepted(false);
    let (p1, a1) = accepted(true);
    assert_eq!((p0, a0), (p1, a1));
    assert!(p0 + a0 < 50, "{p0} + {a0}");
}

#[test]
fn attach_ms2_scores_each_probe_on_its_own_window() {
    use quant::ms2::{fragment_ladder, library_spectrum, residue_masses, Ms2Scan};
    let lib = library();
    let tol = model::Tolerance::Ppm(20.0);
    let ladder = fragment_ladder(
        &residue_masses("LGEYGFQNALIVR", &std::collections::HashMap::new()).unwrap(),
        2,
    );
    let donor: Vec<(f32, f32)> = ladder
        .iter()
        .enumerate()
        .map(|(k, f)| (f.mz as f32, 100.0 + k as f32))
        .collect();
    let spectrum = library_spectrum(&ladder, &donor, tol).unwrap();
    let mut spectra = vec![None; lib.meta.len()];
    spectra[0] = Some(spectrum);
    let mut donor_scans = vec![Vec::new(); lib.meta.len()];
    donor_scans[0] = vec![(0u16, 50.0), (1u16, 100.0)];
    let noise = vec![(321.1f32, 10.0f32), (654.3, 20.0)];
    let target_mz = lib.meta[0].mz;
    let decoy_mz = probe_mz(&lib.meta[0], Hypothesis::MassShiftDecoy(7));
    assert!((decoy_mz - target_mz - 7.0 * PEPTIDE_MASS_SLOPE / 2.0).abs() < 1e-12);
    let scan = |rt: f64, mz: f64, peaks: &Vec<(f32, f32)>| Ms2Scan {
        native_id: format!("scan={rt}"),
        spec_key: format!("k{rt}"),
        rt,
        precursor_mz: mz,
        charge: 2,
        isolation: Some((mz, 0.7, 0.7)),
        faims_cv: None,
        peaks: peaks.clone(),
        reporters: None,
    };
    let index = Ms2RunIndex::new(vec![
        scan(100.0, target_mz, &donor),
        scan(101.0, target_mz, &noise),
        scan(200.0, decoy_mz, &noise),
        scan(400.0, target_mz, &donor),
    ]);
    let ev = |t0: f64| {
        let mut g = Rng(1);
        Some(evidence(&mut g, true, t0))
    };
    let mut rows = vec![
        row(0, 1, Hypothesis::Target, Role::Transfer, ev(100.0)),
        row(
            0,
            1,
            Hypothesis::MassShiftDecoy(7),
            Role::Transfer,
            ev(200.0),
        ),
        row(0, 1, Hypothesis::RtShiftDecoy, Role::Transfer, ev(400.0)),
        row(0, 1, Hypothesis::Target, Role::PseudoTransfer, ev(100.0)),
        row(0, 1, Hypothesis::Target, Role::Identified, ev(100.0)),
        row(0, 1, Hypothesis::Target, Role::Transfer, None),
        row(1, 1, Hypothesis::Target, Role::Transfer, ev(100.0)),
    ];
    let c = attach_ms2(&mut rows, &lib, &spectra, &donor_scans, &index, tol);
    let m = |i: usize| rows[i].ms2.unwrap();
    assert_eq!(m(0).n_scans, 2);
    assert!(m(0).supports());
    assert_eq!(m(1).n_scans, 1);
    assert!(!m(1).supports());
    assert_eq!(m(2).n_scans, 1);
    assert!(m(2).supports());
    assert_eq!(m(3).n_scans, 1);
    assert!(!m(3).supports());
    assert!(rows[4].ms2.is_none() && rows[5].ms2.is_none() && rows[6].ms2.is_none());
    assert_eq!(c.rows[0], [1, 1, 1]);
    assert_eq!(c.supported[0], [1, 0, 1]);
    assert_eq!((c.rows[1][0], c.supported[1][0]), (1, 0));
    assert_eq!(c.no_spectrum, 1);
}
