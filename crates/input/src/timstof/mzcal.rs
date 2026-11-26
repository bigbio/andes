use std::collections::HashMap;

/// The instrument's own TOF to m/z model, a ModelType 1 row of the MzCalibration table.
/// The flight time of a TOF index is `delay + timebase * tof` ns, and with the temperature
/// factor `cf = 1 - 1e-6 * (dC1 (T1 - T1cal) + dC2 (T2 - T2cal))` of the frame,
/// `t = C0 + 1e6 / sqrt(C1 cf) * sqrt(mz + C4) + C2 / cf * mz`. This agrees with the
/// Bruker SDK to 0.02 ppm, the rounding of the stored coefficients.
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
}

impl TdfMzModel {
    pub fn mz(&self, tof: f64, t: [f64; 2]) -> f64 {
        let drift = self.dc[0] * (t[0] - self.t_cal[0]) + self.dc[1] * (t[1] - self.t_cal[1]);
        let cf = 1.0 - 1e-6 * drift;
        let a = self.c2 / cf;
        let b = 1e6 / (self.c1 * cf).sqrt();
        // a u^2 + b u = y with u = sqrt(mz + C4); this root is stable for either sign of a.
        let y = self.delay + self.timebase * tof - self.c0 + a * self.c4;
        let u = 2.0 * y / (b + (b * b + 4.0 * a * y).sqrt());
        u * u - self.c4
    }
}

/// A frame's (T1, T2) for the m/z model. A sensor with no reading takes the calibration row's
/// own temperature, so it adds no correction (a reading of 0 would shift m/z by hundreds of
/// ppm). Without a TDF model it is 0, as before, and the linear model ignores it.
pub(crate) fn frame_temperatures(t: [Option<f64>; 2], model: Option<&TdfMzModel>) -> [f64; 2] {
    let fallback = model.map_or([0.0; 2], |m| m.t_cal);
    [0, 1].map(|k| t[k].filter(|x| x.is_finite()).unwrap_or(fallback[k]))
}

/// The usable ModelType 1 rows of the MzCalibration table and the row each frame uses.
#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct TdfMzModels {
    rows: HashMap<i64, TdfMzModel>,
    frame_row: HashMap<usize, i64>,
    main: Option<i64>,
    /// Frames whose row is missing or unusable; they use the main row.
    pub unusable_frames: usize,
}

impl TdfMzModels {
    /// The row most frames use, or `None` when that row is unusable or there is no table.
    pub fn main(&self) -> Option<TdfMzModel> {
        self.main.and_then(|id| self.rows.get(&id).copied())
    }

    /// The frame's own row, or the main row when the frame's row is unusable or unknown.
    pub fn for_frame(&self, frame: usize) -> Option<TdfMzModel> {
        let main = self.main()?;
        Some(
            self.frame_row
                .get(&frame)
                .and_then(|id| self.rows.get(id))
                .copied()
                .unwrap_or(main),
        )
    }

    pub fn rows_in_use(&self) -> usize {
        let mut ids: Vec<i64> = self.frame_row.values().copied().collect();
        ids.sort_unstable();
        ids.dedup();
        ids.len()
    }
}

fn tdf_row(kind: Option<i64>, v: [f64; 11]) -> Option<TdfMzModel> {
    let [timebase, delay, t1, t2, dc1, dc2, c0, c1, c2, c3, c4] = v;
    let usable = kind == Some(1)
        && c3 == 0.0
        && timebase > 0.0
        && c1 > 0.0
        && [delay, t1, t2, dc1, dc2, c0, c2, c4]
            .iter()
            .all(|x| x.is_finite());
    usable.then_some(TdfMzModel {
        timebase,
        delay,
        t_cal: [t1, t2],
        dc: [dc1, dc2],
        c0,
        c1,
        c2,
        c4,
    })
}

/// Every usable ModelType 1 row (another model type or a set C3 is not usable) and each
/// frame's row. A missing table or column gives no rows, so m/z keeps the linear model.
pub(crate) fn read_tdf_mz_models(conn: &rusqlite::Connection) -> TdfMzModels {
    let mut out = TdfMzModels::default();
    let sql = "SELECT Id, ModelType, DigitizerTimebase, DigitizerDelay, T1, T2, dC1, dC2, C0, \
               C1, C2, C3, C4 FROM MzCalibration";
    let Ok(mut stmt) = conn.prepare(sql) else {
        return out;
    };
    let rows = stmt.query_map([], |r| {
        let mut v = [0.0f64; 11];
        for (k, x) in v.iter_mut().enumerate() {
            *x = r.get::<_, Option<f64>>(k + 2)?.unwrap_or(f64::NAN);
        }
        Ok((r.get::<_, i64>(0)?, r.get::<_, Option<i64>>(1)?, v))
    });
    for (id, kind, v) in rows.into_iter().flatten().flatten() {
        if let Some(model) = tdf_row(kind, v) {
            out.rows.insert(id, model);
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
        if let Some(row) = row {
            out.frame_row.insert(frame as usize, row);
            if out.rows.contains_key(&row) {
                *count.entry(row).or_default() += 1;
            }
        }
    }
    out.main = count
        .iter()
        .max_by_key(|&(id, n)| (*n, std::cmp::Reverse(*id)))
        .map(|(id, _)| *id);
    out.unusable_frames = out
        .frame_row
        .values()
        .filter(|id| !out.rows.contains_key(id))
        .count();
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::timstof::ms1::MzCalibration;

    /// A model, the frame T1 and T2, and TOF index -> m/z pairs of the Bruker SDK.
    type SdkCase = (TdfMzModel, [f64; 2], [(f64, f64); 4]);

    /// SDK cases for frame 101 of PXD028735 Condition A Alpha 01 (C2 > 0) and of a 21 min
    /// HeLa run (C2 < 0).
    fn sdk_cases() -> [SdkCase; 2] {
        let pxd = TdfMzModel {
            timebase: 0.2,
            delay: 25279.6,
            t_cal: [25.803492210895893, 26.16594468081322],
            dc: [13.4, 0.0],
            c0: 317.79692276420536,
            c1: 152422.6786479605,
            c2: 0.00044948185059656226,
            c4: -0.021456006606015886,
        };
        let hela = TdfMzModel {
            timebase: 0.19999999999999998,
            delay: 25581.8,
            t_cal: [25.304233, 22.467882],
            dc: [21.0, 0.0],
            c0: 308.058693,
            c1: 156618.35376,
            c2: -0.003637,
            c4: 0.045331,
        };
        [
            (
                pxd,
                [25.759546401430235, 25.415852921226346],
                [
                    (20000.5, 127.87192146547858),
                    (100000.0, 308.1518988296756),
                    (250000.0, 856.5176411319167),
                    (380000.0, 1553.688785182179),
                ],
            ),
            (
                hela,
                [25.327912672795613, 22.637715307372783],
                [
                    (20000.5, 134.1743960310776),
                    (100000.0, 320.99353264624693),
                    (250000.0, 887.4512495268721),
                    (380000.0, 1606.4751195124252),
                ],
            ),
        ]
    }

    #[test]
    fn the_tdf_model_reproduces_the_bruker_sdk() {
        for (model, t, pairs) in sdk_cases() {
            let cal = MzCalibration::tdf(model, 400_000.0);
            for (tof, sdk) in pairs {
                let ppm = (cal.mz(tof, [t[0], t[1], 0.0, 0.0, 0.0]) - sdk) / sdk * 1e6;
                assert!(ppm.abs() < 0.05, "tof {tof}: {ppm} ppm");
            }
        }
    }

    #[test]
    fn the_tdf_model_drops_the_temperature_term_only_at_a_cost() {
        let (model, t, pairs) = sdk_cases()[0];
        let (tof, sdk) = pairs[2];
        let cold = model.mz(tof, model.t_cal);
        let ppm = (cold - sdk) / sdk * 1e6;
        assert!((ppm - -0.59).abs() < 0.05, "{ppm}");
        assert!(((model.mz(tof, t) - sdk) / sdk * 1e6).abs() < 0.05);
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
            conn.execute("INSERT INTO Frames VALUES (?1, ?2)", [id as i64 + 1, *c])
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
        assert_eq!(m.for_frame(1).unwrap().c0, 317.8);
        assert_eq!(m.for_frame(3).unwrap().c0, 320.5);
        assert_eq!(
            m.for_frame(99).unwrap().c0,
            320.5,
            "an unknown frame takes the main row"
        );
        assert_eq!(m.rows_in_use(), 2);
        assert_eq!(m.unusable_frames, 0);
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
            "2, 2, 0.2, 25279.6, 25.8, 26.1, 13.4, 0, 320.5, 151916.0, 0.0021, 0, -0.08",
        ];
        let m = read_tdf_mz_models(&calibration_db(&rows, &[1, 1, 2]));
        assert_eq!(m.for_frame(3).unwrap().c0, 317.8);
        assert_eq!(m.unusable_frames, 1);
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
    fn an_unknown_calibration_model_falls_back() {
        let other = ["1, 2, 0.2, 25279.6, 25.8, 26.1, 13.4, 0, 317.8, 152422.7, 0.00045, 0, 0"];
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
        let (model, t, pairs) = sdk_cases()[0];
        let (tof, sdk) = pairs[2];
        let at_cal = model.mz(tof, model.t_cal);
        for reading in [[None, None], [None, Some(t[1])]] {
            let ft = frame_temperatures(reading, Some(&model));
            assert_eq!(ft[0], model.t_cal[0]);
            let ppm = (model.mz(tof, ft) - at_cal) / at_cal * 1e6;
            assert!(ppm.abs() < 0.05, "{ppm}");
        }
        let ft = frame_temperatures([Some(t[0]), Some(f64::NAN)], Some(&model));
        assert_eq!(ft, [t[0], model.t_cal[1]]);
        assert!(((model.mz(tof, ft) - sdk) / sdk * 1e6).abs() < 0.05);
        let zero = frame_temperatures([Some(0.0), Some(0.0)], Some(&model));
        assert_eq!(
            zero,
            [0.0, 0.0],
            "a measured 0 is a reading, not a missing one"
        );
        assert_eq!(frame_temperatures([None, None], None), [0.0, 0.0]);
    }
}
