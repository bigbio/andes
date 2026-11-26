//! Scan number to inverse reduced ion mobility from a Bruker `TimsCalibration` table.
//!
//! ModelType 2 uses the vendor rational calibration, while runs without a usable row keep the
//! acquisition-range line supplied by timsrust. The frame foreign key is resolved once per frame
//! so the peak loops do not perform database or hash-map lookups.

use std::collections::HashMap;

use timsrust::converters::{ConvertableDomain, Scan2ImConverter};

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum TimsImModel {
    Linear {
        intercept: f64,
        slope: f64,
    },
    Rational {
        c0: f64,
        c1: f64,
        c2: f64,
        c3: f64,
        c4: f64,
        c6: f64,
        c7: f64,
    },
}

impl TimsImModel {
    pub(crate) fn linear(converter: Scan2ImConverter) -> Self {
        let intercept = converter.convert(0.0);
        Self::Linear {
            intercept,
            slope: converter.convert(1.0) - intercept,
        }
    }

    pub(crate) fn convert(self, scan: f64) -> f64 {
        match self {
            Self::Linear { intercept, slope } => intercept + slope * scan,
            Self::Rational {
                c0,
                c1,
                c2,
                c3,
                c4,
                c6,
                c7,
            } => {
                let v = c2 + ((c3 - c2) / c1) * (scan - c4 - c0);
                let v = if v == 0.0 { 1e-10 } else { v };
                let denominator = c6 + c7 / v.abs();
                let denominator = if denominator == 0.0 {
                    1e-10
                } else {
                    denominator
                };
                1.0 / denominator
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct TimsImModels {
    linear: TimsImModel,
    rows: HashMap<i64, TimsImModel>,
    unusable: HashMap<i64, String>,
    frame_row: HashMap<usize, Option<i64>>,
    main: Option<i64>,
}

impl TimsImModels {
    pub(crate) fn new(linear: TimsImModel) -> Self {
        Self {
            linear,
            rows: HashMap::new(),
            unusable: HashMap::new(),
            frame_row: HashMap::new(),
            main: None,
        }
    }

    pub(crate) fn main(&self) -> Option<TimsImModel> {
        self.main.and_then(|id| self.rows.get(&id).copied())
    }

    pub(crate) fn main_id(&self) -> Option<i64> {
        self.main
    }

    pub(crate) fn for_frame(&self, frame: usize) -> TimsImModel {
        self.frame_row
            .get(&frame)
            .and_then(|id| self.rows.get(id.as_ref()?))
            .copied()
            .or_else(|| self.main())
            .unwrap_or(self.linear)
    }

    pub(crate) fn rows_in_use(&self) -> usize {
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

    pub(crate) fn unusable_summary(&self) -> Option<String> {
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

fn tims_row(kind: Option<i64>, v: [Option<f64>; 10]) -> Result<TimsImModel, String> {
    if kind != Some(2) {
        return Err(kind.map_or_else(|| "NULL ModelType".into(), |k| format!("ModelType {k}")));
    }
    let value = |k: usize| {
        v[k].filter(|x| x.is_finite())
            .ok_or("NULL or non-finite coefficients")
    };
    let (c0, c1, c2, c3, c4, c6, c7) = (
        value(0)?,
        value(1)?,
        value(2)?,
        value(3)?,
        value(4)?,
        value(6)?,
        value(7)?,
    );
    if c1 == 0.0 {
        return Err("C1 is zero".into());
    }
    if c7 == 0.0 {
        return Err("C7 is zero".into());
    }
    Ok(TimsImModel::Rational {
        c0,
        c1,
        c2,
        c3,
        c4,
        c6,
        c7,
    })
}

pub(crate) fn read_tims_im_models(
    conn: &rusqlite::Connection,
    linear: TimsImModel,
) -> TimsImModels {
    let mut out = TimsImModels::new(linear);
    let sql = "SELECT Id, ModelType, C0, C1, C2, C3, C4, C5, C6, C7, C8, C9 FROM TimsCalibration";
    let Ok(mut stmt) = conn.prepare(sql) else {
        return out;
    };
    let rows = stmt.query_map([], |r| {
        let mut v = [None; 10];
        for (k, x) in v.iter_mut().enumerate() {
            *x = r.get::<_, Option<f64>>(k + 2)?;
        }
        Ok((r.get::<_, i64>(0)?, r.get::<_, Option<i64>>(1)?, v))
    });
    for (id, kind, v) in rows.into_iter().flatten().flatten() {
        match tims_row(kind, v) {
            Ok(model) => {
                out.rows.insert(id, model);
            }
            Err(why) => {
                out.unusable.insert(id, why);
            }
        }
    }
    if let Ok(mut stmt) = conn.prepare("SELECT Id, TimsCalibration FROM Frames") {
        if let Ok(frames) = stmt.query_map([], |r| {
            Ok((r.get::<_, i64>(0)?, r.get::<_, Option<i64>>(1)?))
        }) {
            for (frame, row) in frames.into_iter().flatten() {
                out.frame_row.insert(frame as usize, row);
            }
        }
    }
    let mut count: HashMap<i64, usize> = HashMap::new();
    for row in out.frame_row.values().flatten() {
        if out.rows.contains_key(row) {
            *count.entry(*row).or_default() += 1;
        }
    }
    out.main = count
        .iter()
        .max_by_key(|&(id, n)| (*n, std::cmp::Reverse(*id)))
        .map(|(id, _)| *id)
        .or_else(|| out.rows.keys().min().copied());
    out
}

pub(crate) fn log_tims_im_calibration(models: &TimsImModels) {
    match models.main() {
        Some(_) => {
            if models.rows_in_use() > 1 {
                log::info!(
                    "frames use {} TimsCalibration rows; each frame's mobility uses its own row",
                    models.rows_in_use()
                );
            }
            if let Some(summary) = models.unusable_summary() {
                log::warn!(
                    "TimsCalibration rows this reader cannot use: {summary}. They use row {}, the row most frames use",
                    models.main_id().unwrap_or_default()
                );
            }
            log::info!("1/K0 from the instrument TimsCalibration rational model (ModelType 2)");
        }
        None => {
            if let Some(summary) = models.unusable_summary() {
                log::warn!(
                    "TimsCalibration rows this reader cannot use: {summary}; 1/K0 uses the acquisition-range model"
                );
            } else {
                log::warn!(
                    "no usable TimsCalibration row in analysis.tdf; 1/K0 uses the acquisition-range model"
                );
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn linear() -> TimsImModel {
        TimsImModel::Linear {
            intercept: 1.0,
            slope: -0.001,
        }
    }

    fn db() -> rusqlite::Connection {
        let conn = rusqlite::Connection::open_in_memory().unwrap();
        conn.execute_batch(
            "CREATE TABLE TimsCalibration (
                Id INTEGER, ModelType INTEGER, C0 REAL, C1 REAL, C2 REAL, C3 REAL,
                C4 REAL, C5 REAL, C6 REAL, C7 REAL, C8 REAL, C9 REAL
            );
            CREATE TABLE Frames (Id INTEGER, TimsCalibration INTEGER);
            INSERT INTO TimsCalibration VALUES
                (1, 2, 1, 909, 195.764817832641, 88.27535992154185, 33, 1,
                 0.03392383094913572, 128.36912660898196, 12.770219268312093,
                 3391.9821773547433),
                (2, 3, 1, 909, 1, 2, 3, 4, 5, 6, 7, 8);
            INSERT INTO Frames VALUES (10, 1), (11, 2), (12, NULL);",
        )
        .unwrap();
        conn
    }

    #[test]
    fn rational_model_matches_hand_calculation() {
        let model = TimsImModel::Rational {
            c0: 1.0,
            c1: 909.0,
            c2: 195.764817832641,
            c3: 88.27535992154185,
            c4: 33.0,
            c6: 0.03392383094913572,
            c7: 128.36912660898196,
        };
        let scan = 450.0;
        let v: f64 = 195.764817832641
            + ((88.27535992154185 - 195.764817832641) / 909.0) * (scan - 33.0 - 1.0);
        let expected = 1.0 / (0.03392383094913572 + 128.36912660898196 / v.abs());
        assert!((model.convert(scan) - expected).abs() < 1e-15);
    }

    #[test]
    #[allow(clippy::excessive_precision)]
    fn rational_model_matches_sdk_samples() {
        let cases = [
            (
                TimsImModel::Rational {
                    c0: 1.0,
                    c1: 909.0,
                    c2: 195.764817832641,
                    c3: 88.27535992154185,
                    c4: 33.0,
                    c6: 0.03392383094913572,
                    c7: 128.36912660898196,
                },
                [
                    (0.0, 1.4782858799),
                    (455.0, 1.0949592212),
                    (910.0, 0.7009910276),
                ],
            ),
            (
                TimsImModel::Rational {
                    c0: 1.0,
                    c1: 935.0,
                    c2: 222.18858702860615,
                    c3: 111.17632025272353,
                    c4: 33.64485981308411,
                    c6: 0.059633991196902866,
                    c7: 157.66430547139623,
                },
                [
                    (0.0, 1.3221692962),
                    (468.0, 1.0172209066),
                    (936.0, 0.6999895306),
                ],
            ),
            (
                TimsImModel::Rational {
                    c0: 1.0,
                    c1: 935.0,
                    c2: 217.04878684486127,
                    c3: 110.01572505598676,
                    c4: 33.64485981308411,
                    c6: 0.04158491601354941,
                    c7: 157.93464969349154,
                },
                [
                    (0.0, 1.3224473855),
                    (468.0, 1.0154233984),
                    (936.0, 0.6999896557),
                ],
            ),
            (
                TimsImModel::Rational {
                    c0: 1.0,
                    c1: 909.0,
                    c2: 196.54809695879558,
                    c3: 89.36111344088289,
                    c4: 33.0,
                    c6: 0.02398539468508476,
                    c7: 130.83612801556208,
                },
                [
                    (0.0, 1.4785282477),
                    (455.0, 1.0933698281),
                    (910.0, 0.7007620917),
                ],
            ),
        ];
        for (model, samples) in cases {
            for (scan, expected) in samples {
                assert!((model.convert(scan) - expected).abs() < 1e-9);
            }
        }
    }

    #[test]
    fn rational_model_uses_absolute_velocity() {
        let model = TimsImModel::Rational {
            c0: 0.0,
            c1: 1.0,
            c2: -2.0,
            c3: -3.0,
            c4: 0.0,
            c6: 0.5,
            c7: 4.0,
        };
        assert!((model.convert(1.0) - 1.0 / (0.5 + 4.0 / 3.0)).abs() < 1e-15);
    }

    #[test]
    fn frames_use_usable_rows_and_linear_fallback() {
        let models = read_tims_im_models(&db(), linear());
        assert_eq!(models.main_id(), Some(1));
        assert!(matches!(models.for_frame(10), TimsImModel::Rational { .. }));
        assert_eq!(models.for_frame(11), models.for_frame(12));
        assert_eq!(models.for_frame(99), models.main().unwrap());
        assert_eq!(models.rows_in_use(), 1);
        assert_eq!(
            models.unusable_summary().as_deref(),
            Some("1 frame on row 2 (ModelType 3), 1 frame with no row")
        );
    }

    #[test]
    fn missing_table_keeps_linear_model() {
        let conn = rusqlite::Connection::open_in_memory().unwrap();
        let models = read_tims_im_models(&conn, linear());
        assert_eq!(models.for_frame(1), linear());
        assert_eq!(models.main(), None);
    }
}
