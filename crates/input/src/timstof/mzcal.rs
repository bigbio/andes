//! TOF index to m/z from the MzCalibration table of a Bruker `.d`: ModelType 1 and 2 rows,
//! the frame temperature correction, the ModelType 2 correction polynomial, and the row
//! each frame uses.

use std::collections::HashMap;

/// The instrument's own TOF to m/z model, a ModelType 1 or 2 row of the MzCalibration table.
/// The flight time of a TOF index is `delay + timebase * tof` ns, and with the temperature
/// factor `cf = 1 - 1e-6 * (dC1 (T1 - T1cal) + dC2 (T2 - T2cal))` of the frame,
/// `t = C0 + b s + C2 / cf * s^2` with `b = 1e6 / sqrt(C1 cf)` and `s = sqrt(mz + C4)`.
/// A ModelType 2 row has C4 = 0 and may add a correction polynomial. This reproduces the
/// Bruker SDK to about 1e-5 ppm.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct TdfMzModel {
    pub timebase: f64,
    pub delay: f64,
    pub t_cal: [f64; 2],
    pub dc: [f64; 2],
    pub c0: f64,
    pub c1: f64,
    pub c2: f64,
    pub c4: f64,
    pub correction: Option<MzCorrection>,
}

/// The ModelType 2 correction in C5..C14: inside `[low, high]` the base curve m/z `m` becomes
/// `m - sum_k coef[k] m^k`, constant first.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct MzCorrection {
    pub low: f64,
    pub high: f64,
    pub n: usize,
    pub coef: [f64; 7],
}

impl MzCorrection {
    fn at(&self, m: f64) -> f64 {
        if m < self.low || m > self.high {
            return 0.0;
        }
        self.coef[..self.n]
            .iter()
            .rev()
            .fold(0.0, |acc, c| acc * m + c)
    }

    /// C5..C14 of a row, or `None` when the range, the count (1 to 7, the C8..C14 columns)
    /// or a coefficient is unusable.
    fn from_columns(v: [Option<f64>; 10]) -> Option<Self> {
        let (low, high, n) = (v[0]?, v[1]?, v[2]?);
        let usable = low.is_finite() && high.is_finite() && low < high;
        if !usable || n.fract() != 0.0 || !(1.0..=7.0).contains(&n) {
            return None;
        }
        let n = n as usize;
        let mut coef = [0.0; 7];
        for k in 0..n {
            coef[k] = v[3 + k].filter(|x| x.is_finite())?;
        }
        Some(Self { low, high, n, coef })
    }
}

impl TdfMzModel {
    pub fn mz(&self, tof: f64, t: [f64; 2]) -> f64 {
        let drift = self.dc[0] * (t[0] - self.t_cal[0]) + self.dc[1] * (t[1] - self.t_cal[1]);
        let cf = 1.0 - 1e-6 * drift;
        let a = self.c2 / cf;
        let b = 1e6 / (self.c1 * cf).sqrt();
        // a u^2 + b u = y with u = sqrt(mz + C4); this root is stable for either sign of a.
        let y = self.delay + self.timebase * tof - self.c0;
        let u = 2.0 * y / (b + (b * b + 4.0 * a * y).sqrt());
        let m = u * u - self.c4;
        m - self.correction.map_or(0.0, |c| c.at(m))
    }
}

/// A frame's (T1, T2) for the m/z model. A sensor with no reading takes the calibration row's
/// own temperature, so it adds no correction (a reading of 0 would shift m/z by hundreds of
/// ppm). Without a TDF model it is 0, as before, and the linear model ignores it.
pub(crate) fn frame_temperatures(t: [Option<f64>; 2], model: Option<&TdfMzModel>) -> [f64; 2] {
    let fallback = model.map_or([0.0; 2], |m| m.t_cal);
    [0, 1].map(|k| t[k].filter(|x| x.is_finite()).unwrap_or(fallback[k]))
}

/// The usable rows of the MzCalibration table, the reason each other row is not usable, and
/// the row each frame uses.
#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct TdfMzModels {
    rows: HashMap<i64, TdfMzModel>,
    unusable: HashMap<i64, String>,
    /// The frame's row, `None` for a NULL Frames.MzCalibration.
    frame_row: HashMap<usize, Option<i64>>,
    main: Option<i64>,
    /// Frames whose row is NULL, missing or unusable; they use the main row.
    pub unusable_frames: usize,
    /// ModelType 2 rows in use without a usable C5..C14 correction.
    pub without_correction: usize,
}

impl TdfMzModels {
    /// The row most frames use, or `None` when that row is unusable or there is no table.
    pub fn main(&self) -> Option<TdfMzModel> {
        self.main.and_then(|id| self.rows.get(&id).copied())
    }

    /// The Id of the main row.
    pub fn main_id(&self) -> Option<i64> {
        self.main
    }

    /// The frame's own row, or the main row when the frame's row is unusable or unknown.
    pub fn for_frame(&self, frame: usize) -> Option<TdfMzModel> {
        let main = self.main()?;
        Some(
            self.frame_row
                .get(&frame)
                .and_then(|id| self.rows.get(id.as_ref()?))
                .copied()
                .unwrap_or(main),
        )
    }

    /// The usable rows at least one frame uses.
    pub fn rows_in_use(&self) -> usize {
        let mut ids: Vec<i64> = self
            .frame_row
            .values()
            .flatten()
            .filter(|id| self.rows.contains_key(id))
            .copied()
            .collect();
        ids.sort_unstable();
        ids.dedup();
        ids.len()
    }

    /// The frames on an unusable row, per row and reason, and the frames with no row, in one
    /// line, e.g. "12 frames on row 2 (ModelType 3), 1 frame with no row". `None` when every
    /// frame has a usable row.
    pub fn unusable_summary(&self) -> Option<String> {
        let mut per_row: HashMap<i64, usize> = HashMap::new();
        let mut no_row = 0;
        for id in self.frame_row.values() {
            match id {
                Some(id) if self.rows.contains_key(id) => {}
                Some(id) => *per_row.entry(*id).or_default() += 1,
                None => no_row += 1,
            }
        }
        let frames = |n: usize| if n == 1 { "frame" } else { "frames" };
        let mut rows: Vec<(i64, usize)> = per_row.into_iter().collect();
        rows.sort_by_key(|&(id, n)| (std::cmp::Reverse(n), id));
        let mut parts: Vec<String> = rows
            .iter()
            .take(3)
            .map(|&(id, n)| {
                let why = self.unusable.get(&id).map_or("no such row", String::as_str);
                format!("{n} {} on row {id} ({why})", frames(n))
            })
            .collect();
        if rows.len() > 3 {
            let n: usize = rows[3..].iter().map(|r| r.1).sum();
            parts.push(format!(
                "{n} {} on {} other rows",
                frames(n),
                rows.len() - 3
            ));
        }
        if no_row > 0 {
            parts.push(format!("{no_row} {} with no row", frames(no_row)));
        }
        (!parts.is_empty()).then(|| parts.join(", "))
    }
}

/// A row from ModelType and DigitizerTimebase..C4, or the reason it is not usable. A NULL C3
/// or C4 reads as 0, as OpenMS does. In ModelType 2, C3 and C4 repeat C0 and C2, so they are
/// ignored there.
fn tdf_row(kind: Option<i64>, v: [Option<f64>; 11]) -> Result<TdfMzModel, String> {
    match kind {
        Some(1) | Some(2) => {}
        Some(k) => return Err(format!("ModelType {k}")),
        None => return Err("NULL ModelType".into()),
    }
    let mut base = [0.0; 9];
    for (x, y) in base.iter_mut().zip(v) {
        *x = y
            .filter(|y| y.is_finite())
            .ok_or("NULL or non-finite coefficients")?;
    }
    let [timebase, delay, t1, t2, dc1, dc2, c0, c1, c2] = base;
    let (c3, c4) = match kind {
        Some(1) => (v[9].unwrap_or(0.0), v[10].unwrap_or(0.0)),
        _ => (0.0, 0.0),
    };
    if !c4.is_finite() {
        return Err("NULL or non-finite coefficients".into());
    }
    if c3 != 0.0 {
        return Err("a set C3".into());
    }
    if !(timebase > 0.0 && c1 > 0.0) {
        return Err("C1 or DigitizerTimebase not positive".into());
    }
    Ok(TdfMzModel {
        timebase,
        delay,
        t_cal: [t1, t2],
        dc: [dc1, dc2],
        c0,
        c1,
        c2,
        c4,
        correction: None,
    })
}

/// The C5..C14 correction of each row that has one. These columns exist only in files with
/// ModelType 2 rows, so a failed query means no correction.
fn read_corrections(conn: &rusqlite::Connection) -> HashMap<i64, MzCorrection> {
    let sql = "SELECT Id, C5, C6, C7, C8, C9, C10, C11, C12, C13, C14 FROM MzCalibration";
    let Ok(mut stmt) = conn.prepare(sql) else {
        return HashMap::new();
    };
    let rows = stmt.query_map([], |r| {
        let mut v = [None; 10];
        for (k, x) in v.iter_mut().enumerate() {
            *x = r.get::<_, Option<f64>>(k + 1)?;
        }
        Ok((r.get::<_, i64>(0)?, v))
    });
    rows.into_iter()
        .flatten()
        .flatten()
        .filter_map(|(id, v)| Some((id, MzCorrection::from_columns(v)?)))
        .collect()
}

/// Every usable ModelType 1 or 2 row, the reason each other row is not usable, and each
/// frame's row. A missing table or column gives no rows, so m/z keeps the linear model.
pub(crate) fn read_tdf_mz_models(conn: &rusqlite::Connection) -> TdfMzModels {
    let mut out = TdfMzModels::default();
    let sql = "SELECT Id, ModelType, DigitizerTimebase, DigitizerDelay, T1, T2, dC1, dC2, C0, \
               C1, C2, C3, C4 FROM MzCalibration";
    let Ok(mut stmt) = conn.prepare(sql) else {
        return out;
    };
    let rows = stmt.query_map([], |r| {
        let mut v = [None; 11];
        for (k, x) in v.iter_mut().enumerate() {
            *x = r.get::<_, Option<f64>>(k + 2)?;
        }
        Ok((r.get::<_, i64>(0)?, r.get::<_, Option<i64>>(1)?, v))
    });
    let corrections = read_corrections(conn);
    let mut type2 = Vec::new();
    for (id, kind, v) in rows.into_iter().flatten().flatten() {
        match tdf_row(kind, v) {
            Ok(mut model) => {
                if kind == Some(2) {
                    model.correction = corrections.get(&id).copied();
                    if model.correction.is_none() {
                        type2.push(id);
                    }
                }
                out.rows.insert(id, model);
            }
            Err(why) => {
                out.unusable.insert(id, why);
            }
        }
    }
    let Ok(mut stmt) = conn.prepare("SELECT Id, MzCalibration FROM Frames") else {
        return out;
    };
    let frames = stmt.query_map([], |r| {
        Ok((r.get::<_, i64>(0)?, r.get::<_, Option<i64>>(1)?))
    });
    let mut count: HashMap<i64, usize> = HashMap::new();
    for (frame, row) in frames.into_iter().flatten().flatten() {
        out.frame_row.insert(frame as usize, row);
        if let Some(row) = row.filter(|r| out.rows.contains_key(r)) {
            *count.entry(row).or_default() += 1;
        }
    }
    // With no frame on a usable row, a usable row is still closer than the acquisition-range
    // line, so take the lowest Id.
    out.main = count
        .iter()
        .max_by_key(|&(id, n)| (*n, std::cmp::Reverse(*id)))
        .map(|(id, _)| *id)
        .or_else(|| out.rows.keys().min().copied());
    out.unusable_frames = out
        .frame_row
        .values()
        .filter(|id| !id.is_some_and(|id| out.rows.contains_key(&id)))
        .count();
    out.without_correction = type2.iter().filter(|id| count.contains_key(id)).count();
    out
}

#[cfg(test)]
// The SDK and OpenMS values are kept as printed.
#[allow(clippy::excessive_precision)]
mod tests {
    use super::*;
    use crate::timstof::ms1::MzCalibration;

    /// A model, the frame T1 and T2, and TOF index -> m/z pairs of the Bruker SDK.
    type SdkCase = (TdfMzModel, [f64; 2], Vec<(f64, f64)>);

    fn type1(v: [f64; 9]) -> TdfMzModel {
        let [timebase, delay, t1, t2, dc1, dc2, c0, c1, c2] = v;
        TdfMzModel {
            timebase,
            delay,
            t_cal: [t1, t2],
            dc: [dc1, dc2],
            c0,
            c1,
            c2,
            c4: 0.0,
            correction: None,
        }
    }

    /// SDK cases for frame 101 of PXD028735 Condition A Alpha 01 (C2 > 0), of a 21 min HeLa
    /// run (C2 < 0), and the ModelType 1 row of the OpenMS tests with dC2 0 (A) and 25 (B).
    fn sdk_cases() -> Vec<SdkCase> {
        let pxd = TdfMzModel {
            c4: -0.021456006606015886,
            ..type1([
                0.2,
                25279.6,
                25.803492210895893,
                26.16594468081322,
                13.4,
                0.0,
                317.79692276420536,
                152422.6786479605,
                0.00044948185059656226,
            ])
        };
        let hela = TdfMzModel {
            c4: 0.045331,
            ..type1([
                0.19999999999999998,
                25581.8,
                25.304233,
                22.467882,
                21.0,
                0.0,
                308.058693,
                156618.35376,
                -0.003637,
            ])
        };
        let openms = |dc2: f64| TdfMzModel {
            c4: -0.068565,
            ..type1([
                0.2,
                24833.0,
                25.727448,
                28.353124,
                27.0,
                dc2,
                320.480481,
                157992.031315,
                0.002115,
            ])
        };
        let tofs = [1.0, 50000.0, 150000.0, 250000.0, 350000.0];
        let a = [
            95.00007900187677,
            188.25046260982114,
            469.54221185610277,
            877.2135834525166,
            1411.2595103197655,
        ];
        let b = [
            95.0039059785749,
            188.25804879176218,
            469.56113775995885,
            877.2489438234296,
            1411.31639969861,
        ];
        vec![
            (
                pxd,
                [25.759546401430235, 25.415852921226346],
                vec![
                    (20000.5, 127.87192146547858),
                    (100000.0, 308.1518988296756),
                    (250000.0, 856.5176411319167),
                    (380000.0, 1553.688785182179),
                ],
            ),
            (
                hela,
                [25.327912672795613, 22.637715307372783],
                vec![
                    (20000.5, 134.1743960310776),
                    (100000.0, 320.99353264624693),
                    (250000.0, 887.4512495268721),
                    (380000.0, 1606.4751195124252),
                ],
            ),
            (
                openms(0.0),
                [25.783137561298357, 28.876605339990224],
                tofs.into_iter().zip(a).collect(),
            ),
            (
                openms(25.0),
                [26.583137561298358, 25.876605339990224],
                tofs.into_iter().zip(b).collect(),
            ),
        ]
    }

    fn ppm(mz: f64, sdk: f64) -> f64 {
        (mz - sdk) / sdk * 1e6
    }

    #[test]
    fn the_tdf_model_reproduces_the_bruker_sdk() {
        for (model, t, pairs) in sdk_cases() {
            let cal = MzCalibration::tdf(model, 400_000.0);
            for (tof, sdk) in pairs {
                let ppm = ppm(cal.mz(tof, [t[0], t[1], 0.0, 0.0, 0.0]), sdk);
                assert!(ppm.abs() < 1e-4, "tof {tof}: {ppm} ppm");
            }
        }
    }

    #[test]
    fn the_tdf_model_drops_the_temperature_term_only_at_a_cost() {
        let (model, t, pairs) = sdk_cases().swap_remove(0);
        let (tof, sdk) = pairs[2];
        let cold = ppm(model.mz(tof, model.t_cal), sdk);
        assert!((cold - -0.5889).abs() < 0.001, "{cold}");
        assert!(ppm(model.mz(tof, t), sdk).abs() < 1e-4);
    }

    #[test]
    fn the_t2_term_counts() {
        let (model, t, pairs) = sdk_cases().swap_remove(3);
        let without = TdfMzModel {
            dc: [model.dc[0], 0.0],
            ..model
        };
        let (tof, sdk) = pairs[2];
        let lost = ppm(without.mz(tof, t), sdk);
        assert!((lost - -62.0).abs() < 1.0, "{lost}");
    }

    const C8_C14: [f64; 7] = [
        0.058519441907584555,
        -0.0005411344569873044,
        1.8634987450927556e-06,
        -3.1297916646308927e-09,
        2.7471965457754794e-12,
        -1.207997933655641e-15,
        2.095848318930801e-19,
    ];

    /// The ModelType 2 row of the OpenMS tests, its frame T1 and T2 and the SDK m/z.
    fn type2_case() -> SdkCase {
        let model = TdfMzModel {
            correction: Some(MzCorrection {
                low: 225.951491,
                high: 1519.712539,
                n: 7,
                coef: C8_C14,
            }),
            ..type1([
                0.19999999999999998,
                25726.199999999997,
                25.66893460672872,
                25.10638483329421,
                20.0,
                0.0,
                313.51167546836524,
                154833.4542880268,
                -7.3593189552069465e-06,
            ])
        };
        let tofs = [
            1.0, 20000.0, 50000.0, 150000.0, 250000.0, 385000.0, 396000.0,
        ];
        let sdk = [
            99.9937348434843,
            133.94736657275706,
            194.17018217157815,
            475.42674239965703,
            880.5493429108868,
            1623.9488391924265,
            1694.4686480927487,
        ];
        (
            model,
            [25.67976562246388, 25.350275695575547],
            tofs.into_iter().zip(sdk).collect(),
        )
    }

    #[test]
    fn a_model_type_2_row_with_its_correction_reproduces_the_sdk() {
        let (model, t, pairs) = type2_case();
        for &(tof, sdk) in &pairs {
            let ppm = ppm(model.mz(tof, t), sdk);
            assert!(ppm.abs() < 1e-4, "tof {tof}: {ppm} ppm");
        }
        let base = TdfMzModel {
            correction: None,
            ..model
        };
        let (tof, sdk) = pairs[3];
        let off = ppm(base.mz(tof, t), sdk);
        assert!((off - -0.927).abs() < 0.001, "{off}");
    }

    const TYPE2_ROW: &str = "1, 2, 0.19999999999999998, 25726.199999999997, \
        25.66893460672872, 25.10638483329421, 20.0, 0.0, 313.51167546836524, \
        154833.4542880268, -7.3593189552069465e-06, 313.51167546836524, \
        -7.3593189552069465e-06";

    #[test]
    fn a_model_type_2_row_reads_its_correction_when_the_columns_exist() {
        let (model, t, pairs) = type2_case();
        let conn = calibration_db(&[TYPE2_ROW], &[1]);
        let m = read_tdf_mz_models(&conn);
        let base = m
            .main()
            .expect("a ModelType 2 row is usable without C5..C14");
        assert_eq!(base.c4, 0.0, "C3 and C4 of ModelType 2 are ignored");
        assert_eq!(base.correction, None);
        assert_eq!(m.without_correction, 1);
        let mut cols = String::new();
        for c in 5..=14 {
            conn.execute(
                &format!("ALTER TABLE MzCalibration ADD COLUMN C{c} REAL"),
                [],
            )
            .unwrap();
            cols.push_str(&format!("C{c} = ?{}, ", c - 4));
        }
        let mut v = vec![225.951491, 1519.712539, 7.0];
        v.extend(C8_C14);
        let sql = format!(
            "UPDATE MzCalibration SET {} WHERE Id = 1",
            cols.trim_end_matches(", ")
        );
        conn.execute(&sql, rusqlite::params_from_iter(v)).unwrap();
        let m = read_tdf_mz_models(&conn);
        assert_eq!(m.main(), Some(model));
        assert_eq!(m.without_correction, 0);
        for (tof, sdk) in pairs {
            assert!(ppm(m.main().unwrap().mz(tof, t), sdk).abs() < 1e-4);
        }
    }

    #[test]
    fn a_correction_needs_a_range_a_count_and_its_coefficients() {
        let mut v = [None; 10];
        v[..3].copy_from_slice(&[Some(200.0), Some(1500.0), Some(2.0)]);
        v[3] = Some(0.01);
        v[4] = Some(1e-5);
        let c = MzCorrection::from_columns(v).unwrap();
        assert_eq!(c.at(1000.0), 0.01 + 1e-5 * 1000.0);
        assert_eq!(c.at(100.0), 0.0, "outside the calibrant range");
        let mut no_coef = v;
        no_coef[4] = None;
        assert_eq!(MzCorrection::from_columns(no_coef), None);
        let mut bad_n = v;
        bad_n[2] = Some(11.0);
        assert_eq!(MzCorrection::from_columns(bad_n), None);
        let mut bad_range = v;
        bad_range[1] = Some(100.0);
        assert_eq!(MzCorrection::from_columns(bad_range), None);
    }

    fn calibration_db(rows: &[&str], frames: &[i64]) -> rusqlite::Connection {
        let conn = rusqlite::Connection::open_in_memory().unwrap();
        conn.execute_batch(
            "CREATE TABLE MzCalibration (Id INTEGER, ModelType INTEGER, DigitizerTimebase REAL, \
               DigitizerDelay REAL, T1 REAL, T2 REAL, dC1 REAL, dC2 REAL, C0 REAL, C1 REAL, \
               C2 REAL, C3 REAL, C4 REAL);
             CREATE TABLE Frames (Id INTEGER, MzCalibration INTEGER);",
        )
        .unwrap();
        for r in rows {
            conn.execute(&format!("INSERT INTO MzCalibration VALUES ({r})"), [])
                .unwrap();
        }
        for (id, c) in frames.iter().enumerate() {
            let c = (*c > 0).then_some(*c);
            conn.execute(
                "INSERT INTO Frames VALUES (?1, ?2)",
                rusqlite::params![id as i64 + 1, c],
            )
            .unwrap();
        }
        conn
    }

    #[test]
    fn each_frame_uses_its_own_calibration_row() {
        let rows = [
            "1, 1, 0.2, 25279.6, 25.8, 26.1, 13.4, 0, 317.8, 152422.7, 0.00045, 0, -0.02",
            "2, 1, 0.2, 25279.6, 25.8, 26.1, 13.4, 0, 320.5, 151916.0, 0.0021, 0, -0.08",
        ];
        let m = read_tdf_mz_models(&calibration_db(&rows, &[1, 2, 2]));
        assert_eq!(m.main().unwrap().c0, 320.5);
        assert_eq!(m.main_id(), Some(2));
        assert_eq!(m.for_frame(1).unwrap().c0, 317.8);
        assert_eq!(m.for_frame(3).unwrap().c0, 320.5);
        assert_eq!(
            m.for_frame(99).unwrap().c0,
            320.5,
            "an unknown frame takes the main row"
        );
        assert_eq!(m.rows_in_use(), 2);
        assert_eq!(m.unusable_frames, 0);
        assert_eq!(m.unusable_summary(), None);
        let m = read_tdf_mz_models(&calibration_db(&rows, &[1, 1, 2]));
        assert_eq!(m.main().unwrap().c0, 317.8);
        let base = MzCalibration::tdf(m.main().unwrap(), 400_000.0);
        let own = base.for_frame(m.for_frame(3));
        let t = [25.8, 26.1, 0.0, 0.0, 0.0];
        assert_ne!(base.mz(250_000.0, t), own.mz(250_000.0, t));
        let linear = MzCalibration::linear(10.0, 1e-4, 1.0);
        assert_eq!(linear.for_frame(m.for_frame(3)), linear);
    }

    #[test]
    fn a_frame_on_an_unusable_row_uses_the_main_row() {
        let rows = [
            "1, 1, 0.2, 25279.6, 25.8, 26.1, 13.4, 0, 317.8, 152422.7, 0.00045, 0, -0.02",
            "2, 3, 0.2, 25279.6, 25.8, 26.1, 13.4, 0, 320.5, 151916.0, 0.0021, 0, -0.08",
        ];
        let m = read_tdf_mz_models(&calibration_db(&rows, &[1, 1, 2]));
        assert_eq!(m.for_frame(3).unwrap().c0, 317.8);
        assert_eq!(m.unusable_frames, 1);
        assert_eq!(m.rows_in_use(), 1, "only usable rows count");
        let m = read_tdf_mz_models(&calibration_db(&rows, &[1, 2, 2]));
        assert_eq!(
            m.main().unwrap().c0,
            317.8,
            "an unusable row is never the main row, however many frames use it"
        );
        assert_eq!(m.for_frame(1).unwrap().c0, 317.8);
        assert_eq!(m.unusable_frames, 2);
    }

    #[test]
    fn a_frame_with_no_row_uses_the_main_row() {
        let rows = ["1, 1, 0.2, 25279.6, 25.8, 26.1, 13.4, 0, 317.8, 152422.7, 0.00045, 0, -0.02"];
        let m = read_tdf_mz_models(&calibration_db(&rows, &[1, 0]));
        assert_eq!(m.for_frame(2).unwrap().c0, 317.8);
        assert_eq!(m.unusable_frames, 1);
        assert_eq!(m.unusable_summary().unwrap(), "1 frame with no row");
    }

    #[test]
    fn with_no_frame_on_a_usable_row_the_lowest_usable_row_is_the_main_row() {
        let rows = [
            "2, 1, 0.2, 25279.6, 25.8, 26.1, 13.4, 0, 320.5, 151916.0, 0.0021, 0, -0.08",
            "1, 1, 0.2, 25279.6, 25.8, 26.1, 13.4, 0, 317.8, 152422.7, 0.00045, 0, -0.02",
            "3, 3, 0.2, 25279.6, 25.8, 26.1, 13.4, 0, 320.5, 151916.0, 0.0021, 0, -0.08",
        ];
        for frames in [[0, 0], [3, 3]] {
            let m = read_tdf_mz_models(&calibration_db(&rows, &frames));
            assert_eq!(m.main_id(), Some(1), "frames {frames:?}");
            assert_eq!(m.for_frame(1).unwrap().c0, 317.8);
            assert_eq!(m.unusable_frames, 2);
        }
    }

    #[test]
    fn the_warning_lists_frames_per_row_and_reason() {
        let ok = "1, 1, 0.2, 25279.6, 25.8, 26.1, 13.4, 0, 317.8, 152422.7, 0.00045, 0, -0.02";
        let rows = [
            ok,
            "2, 3, 0.2, 25279.6, 25.8, 26.1, 13.4, 0, 317.8, 152422.7, 0.00045, 0, 0",
            "3, 1, 0.2, 25279.6, 25.8, 26.1, 13.4, 0, 317.8, 152422.7, 0.00045, 1, 0",
            "4, 1, 0.2, 25279.6, NULL, 26.1, 13.4, 0, 317.8, 152422.7, 0.00045, 0, 0",
            "5, 1, 0, 25279.6, 25.8, 26.1, 13.4, 0, 317.8, 152422.7, 0.00045, 0, 0",
        ];
        let frames = [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 0];
        let m = read_tdf_mz_models(&calibration_db(&rows, &frames));
        assert_eq!(
            m.unusable_summary().unwrap(),
            "3 frames on row 2 (ModelType 3), 2 frames on row 3 (a set C3), \
             1 frame on row 4 (NULL or non-finite coefficients), 1 frame with no row"
        );
        let mut frames = frames.to_vec();
        frames.extend([5, 9, 9]);
        let m = read_tdf_mz_models(&calibration_db(&rows, &frames));
        assert_eq!(
            m.unusable_summary().unwrap(),
            "3 frames on row 2 (ModelType 3), 2 frames on row 3 (a set C3), \
             2 frames on row 9 (no such row), 2 frames on 2 other rows, 1 frame with no row"
        );
        assert_eq!(
            tdf_row(Some(1), [Some(0.0); 11]).unwrap_err(),
            "C1 or DigitizerTimebase not positive"
        );
        assert_eq!(m.unusable_frames, 10);
    }

    #[test]
    fn a_null_c3_or_c4_reads_as_zero() {
        let rows =
            ["1, 1, 0.2, 25279.6, 25.8, 26.1, 13.4, 0, 317.8, 152422.7, 0.00045, NULL, NULL"];
        let m = read_tdf_mz_models(&calibration_db(&rows, &[1]))
            .main()
            .unwrap();
        assert_eq!(m.c4, 0.0);
    }

    #[test]
    fn an_unknown_calibration_model_falls_back() {
        let other = ["1, 3, 0.2, 25279.6, 25.8, 26.1, 13.4, 0, 317.8, 152422.7, 0.00045, 0, 0"];
        assert_eq!(
            read_tdf_mz_models(&calibration_db(&other, &[1])).main(),
            None
        );
        let c3 = ["1, 1, 0.2, 25279.6, 25.8, 26.1, 13.4, 0, 317.8, 152422.7, 0.00045, 1, 0"];
        assert_eq!(read_tdf_mz_models(&calibration_db(&c3, &[1])).main(), None);
        let conn = rusqlite::Connection::open_in_memory().unwrap();
        assert_eq!(read_tdf_mz_models(&conn).main(), None);
    }

    #[test]
    fn a_missing_temperature_takes_the_calibration_temperature() {
        let (model, t, pairs) = sdk_cases().swap_remove(0);
        let (tof, sdk) = pairs[2];
        let at_cal = model.mz(tof, model.t_cal);
        for reading in [[None, None], [None, Some(t[1])]] {
            let ft = frame_temperatures(reading, Some(&model));
            assert_eq!(ft[0], model.t_cal[0]);
            assert!(ppm(model.mz(tof, ft), at_cal).abs() < 1e-4);
        }
        let ft = frame_temperatures([Some(t[0]), Some(f64::NAN)], Some(&model));
        assert_eq!(ft, [t[0], model.t_cal[1]]);
        assert!(ppm(model.mz(tof, ft), sdk).abs() < 1e-4);
        let zero = frame_temperatures([Some(0.0), Some(0.0)], Some(&model));
        assert_eq!(
            zero,
            [0.0, 0.0],
            "a measured 0 is a reading, not a missing one"
        );
        assert_eq!(frame_temperatures([None, None], None), [0.0, 0.0]);
    }
}
