use std::collections::HashMap;

use crate::timstof_ms1::{
    centroid_profile, conditions, sum_over_mobility, FrameConditions, MzCalibration,
};

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Selection {
    pub frame_id: usize,
    pub scan_begin: usize,
    pub scan_end: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PasefPrecursor {
    pub id: usize,
    pub mz: f64,
    pub charge: Option<i64>,
    pub intensity: Option<f64>,
    pub parent: Option<usize>,
    pub scan_number: Option<f64>,
    pub isolation_width: f64,
    pub collision_energy: f64,
    pub selections: Vec<Selection>,
    pub dia: bool,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct FrameRow {
    pub time: f64,
    pub num_peaks: u64,
    pub t: [f64; 2],
}

pub(crate) fn read_frame_rows(
    conn: &rusqlite::Connection,
) -> Result<HashMap<usize, FrameRow>, rusqlite::Error> {
    let with_temperatures = conn.prepare("SELECT T1, T2 FROM Frames LIMIT 1").is_ok();
    let sql = if with_temperatures {
        "SELECT Id, Time, NumPeaks, T1, T2 FROM Frames"
    } else {
        "SELECT Id, Time, NumPeaks, NULL, NULL FROM Frames"
    };
    let mut stmt = conn.prepare(sql)?;
    let rows = stmt.query_map([], |r| {
        Ok((
            r.get::<_, i64>(0)?,
            FrameRow {
                time: r.get::<_, Option<f64>>(1)?.unwrap_or(f64::NAN),
                num_peaks: r.get::<_, Option<i64>>(2)?.unwrap_or(0).max(0) as u64,
                t: [
                    r.get::<_, Option<f64>>(3)?.unwrap_or(0.0),
                    r.get::<_, Option<f64>>(4)?.unwrap_or(0.0),
                ],
            },
        ))
    })?;
    let mut out = HashMap::new();
    for row in rows {
        let (id, frame) = row?;
        out.insert(id as usize, frame);
    }
    Ok(out)
}

pub(crate) fn has_table(conn: &rusqlite::Connection, name: &str) -> Result<bool, rusqlite::Error> {
    conn.query_row(
        "SELECT count(*) FROM sqlite_master WHERE type = 'table' AND name = ?1",
        [name],
        |r| r.get::<_, i64>(0),
    )
    .map(|n| n > 0)
}

type DiaWindow = (usize, usize, f64, f64, f64);

pub(crate) fn mobility_slices(begin: usize, end: usize, n: usize) -> Vec<(usize, usize)> {
    if n <= 1 || end - begin < 2 * n {
        return vec![(begin, end)];
    }
    let step = (end - begin).div_ceil(n + 1);
    (0..n)
        .map(|k| {
            let lo = begin + k * step;
            let hi = if k + 1 == n {
                end
            } else {
                (lo + 2 * step).min(end)
            };
            (lo, hi)
        })
        .collect()
}

pub(crate) fn read_dia_windows(
    conn: &rusqlite::Connection,
) -> Result<Vec<PasefPrecursor>, rusqlite::Error> {
    if !has_table(conn, "DiaFrameMsMsInfo")? || !has_table(conn, "DiaFrameMsMsWindows")? {
        return Ok(Vec::new());
    }
    let mut groups: HashMap<i64, Vec<DiaWindow>> = HashMap::new();
    let mut stmt = conn.prepare(
        "SELECT WindowGroup, ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth, \
         CollisionEnergy FROM DiaFrameMsMsWindows ORDER BY WindowGroup, ScanNumBegin",
    )?;
    let rows = stmt.query_map([], |r| {
        Ok((
            r.get::<_, i64>(0)?,
            r.get::<_, i64>(1)?,
            r.get::<_, i64>(2)?,
            r.get::<_, Option<f64>>(3)?,
            r.get::<_, Option<f64>>(4)?,
            r.get::<_, Option<f64>>(5)?,
        ))
    })?;
    for row in rows {
        let (group, begin, end, mz, width, energy) = row?;
        let mz = mz.filter(|m| m.is_finite() && *m > 0.0);
        if let (Some(mz), true) = (mz, begin >= 0 && end > begin) {
            groups.entry(group).or_default().push((
                begin as usize,
                end as usize,
                mz,
                width.unwrap_or(0.0),
                energy.unwrap_or(0.0),
            ));
        }
    }
    let mut parent_of: HashMap<i64, usize> = HashMap::new();
    let mut stmt = conn.prepare("SELECT Id, MsMsType FROM Frames ORDER BY Id")?;
    let rows = stmt.query_map([], |r| {
        Ok((r.get::<_, i64>(0)?, r.get::<_, Option<i64>>(1)?))
    })?;
    let mut last_ms1: Option<usize> = None;
    for row in rows {
        let (id, kind) = row?;
        if kind == Some(0) {
            last_ms1 = Some(id as usize);
        } else if let Some(p) = last_ms1 {
            parent_of.insert(id, p);
        }
    }
    let mut out = Vec::new();
    let mut stmt =
        conn.prepare("SELECT Frame, WindowGroup FROM DiaFrameMsMsInfo ORDER BY Frame")?;
    let rows = stmt.query_map([], |r| Ok((r.get::<_, i64>(0)?, r.get::<_, i64>(1)?)))?;
    for row in rows {
        let (frame, group) = row?;
        for &(begin, end, mz, width, energy) in groups.get(&group).into_iter().flatten() {
            out.push(PasefPrecursor {
                id: out.len() + 1,
                mz,
                charge: None,
                intensity: None,
                parent: parent_of.get(&frame).copied(),
                scan_number: None,
                isolation_width: width,
                collision_energy: energy,
                selections: vec![Selection {
                    frame_id: frame as usize,
                    scan_begin: begin,
                    scan_end: end,
                }],
                dia: true,
            });
        }
    }
    Ok(out)
}

pub(crate) fn slice_dia_windows(precursors: Vec<PasefPrecursor>, n: usize) -> Vec<PasefPrecursor> {
    if n <= 1 || !precursors.iter().all(|p| p.dia) {
        return precursors;
    }
    let mut out = Vec::with_capacity(precursors.len() * n);
    for p in precursors {
        let parts = match p.selections.as_slice() {
            [s] => mobility_slices(s.scan_begin, s.scan_end, n),
            _ => Vec::new(),
        };
        if parts.len() < 2 {
            out.push(PasefPrecursor {
                id: out.len() + 1,
                ..p
            });
            continue;
        }
        for (lo, hi) in parts {
            out.push(PasefPrecursor {
                id: out.len() + 1,
                selections: vec![Selection {
                    scan_begin: lo,
                    scan_end: hi,
                    ..p.selections[0]
                }],
                ..p.clone()
            });
        }
    }
    out
}

pub(crate) fn read_pasef_precursors(
    conn: &rusqlite::Connection,
) -> Result<Vec<PasefPrecursor>, rusqlite::Error> {
    if !has_table(conn, "Precursors")? || !has_table(conn, "PasefFrameMsMsInfo")? {
        return Ok(Vec::new());
    }
    let mut by_id: HashMap<usize, PasefPrecursor> = HashMap::new();
    let mut stmt = conn.prepare(
        "SELECT Id, MonoisotopicMz, LargestPeakMz, Charge, Intensity, Parent, ScanNumber \
         FROM Precursors",
    )?;
    let rows = stmt.query_map([], |r| {
        Ok((
            r.get::<_, i64>(0)?,
            r.get::<_, Option<f64>>(1)?,
            r.get::<_, Option<f64>>(2)?,
            r.get::<_, Option<i64>>(3)?,
            r.get::<_, Option<f64>>(4)?,
            r.get::<_, Option<i64>>(5)?,
            r.get::<_, Option<f64>>(6)?,
        ))
    })?;
    for row in rows {
        let (id, mono, largest, charge, intensity, parent, scan) = row?;
        let mono = mono.filter(|m| m.is_finite() && *m > 0.0);
        let (mz, charge) = match mono {
            Some(m) => (m, charge.filter(|z| *z > 0)),
            None => (largest.unwrap_or(0.0), None),
        };
        let id = id as usize;
        by_id.insert(
            id,
            PasefPrecursor {
                id,
                mz,
                charge,
                intensity,
                parent: parent.map(|p| p as usize),
                scan_number: scan,
                isolation_width: 0.0,
                collision_energy: 0.0,
                selections: Vec::new(),
                dia: false,
            },
        );
    }
    let mut stmt = conn.prepare(
        "SELECT Frame, ScanNumBegin, ScanNumEnd, IsolationWidth, CollisionEnergy, Precursor \
         FROM PasefFrameMsMsInfo ORDER BY Frame, ScanNumBegin",
    )?;
    let rows = stmt.query_map([], |r| {
        Ok((
            r.get::<_, i64>(0)?,
            r.get::<_, i64>(1)?,
            r.get::<_, i64>(2)?,
            r.get::<_, Option<f64>>(3)?,
            r.get::<_, Option<f64>>(4)?,
            r.get::<_, Option<i64>>(5)?,
        ))
    })?;
    for row in rows {
        let (frame, begin, end, width, energy, precursor) = row?;
        let Some(p) = precursor.and_then(|p| by_id.get_mut(&(p as usize))) else {
            continue;
        };
        if begin < 0 || end <= begin {
            continue;
        }
        p.isolation_width = width.unwrap_or(0.0);
        p.collision_energy = energy.unwrap_or(0.0);
        p.selections.push(Selection {
            frame_id: frame as usize,
            scan_begin: begin as usize,
            scan_end: end as usize,
        });
    }
    let mut out: Vec<PasefPrecursor> = by_id
        .into_values()
        .filter(|p| !p.selections.is_empty())
        .collect();
    out.sort_by_key(|p| p.id);
    Ok(out)
}

pub(crate) fn push_scans(
    scan_offsets: &[usize],
    frame_tof: &[u32],
    frame_intensity: &[u32],
    begin: usize,
    end: usize,
    tof: &mut Vec<u32>,
    intensity: &mut Vec<u32>,
) {
    let n_points = frame_tof.len().min(frame_intensity.len());
    let offset = |k: usize| {
        scan_offsets
            .get(k)
            .or(scan_offsets.last())
            .copied()
            .unwrap_or(0)
            .min(n_points)
    };
    let (lo, hi) = (offset(begin), offset(end));
    if lo < hi {
        tof.extend_from_slice(&frame_tof[lo..hi]);
        intensity.extend_from_slice(&frame_intensity[lo..hi]);
    }
}

pub(crate) fn ms2_peaks(
    tof: &[u32],
    intensity: &[u32],
    calibration: &MzCalibration,
    frame: FrameConditions,
    n_frames: usize,
) -> Vec<(f64, f32)> {
    let per_frame = 1.0 / n_frames.max(1) as f64;
    crate::timstof::finite_sorted_peaks(
        centroid_profile(&sum_over_mobility(tof, intensity))
            .iter()
            .map(|c| {
                let cond = conditions(frame, c.intensity * per_frame);
                (calibration.mz(c.tof, cond), c.intensity as f32)
            }),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tdf_with(precursors: &[&str], pasef: &[&str], frames: &[&str]) -> rusqlite::Connection {
        let conn = rusqlite::Connection::open_in_memory().unwrap();
        conn.execute_batch(
            "CREATE TABLE Precursors (Id INTEGER, LargestPeakMz REAL, AverageMz REAL, \
               MonoisotopicMz REAL, Charge INTEGER, ScanNumber REAL, Intensity REAL, Parent INTEGER);
             CREATE TABLE PasefFrameMsMsInfo (Frame INTEGER, ScanNumBegin INTEGER, \
               ScanNumEnd INTEGER, IsolationMz REAL, IsolationWidth REAL, \
               CollisionEnergy REAL, Precursor INTEGER);
             CREATE TABLE Frames (Id INTEGER, Time REAL, NumPeaks INTEGER, T1 REAL, T2 REAL);",
        )
        .unwrap();
        for p in precursors {
            conn.execute(&format!("INSERT INTO Precursors VALUES ({p})"), [])
                .unwrap();
        }
        for p in pasef {
            conn.execute(&format!("INSERT INTO PasefFrameMsMsInfo VALUES ({p})"), [])
                .unwrap();
        }
        for f in frames {
            conn.execute(&format!("INSERT INTO Frames VALUES ({f})"), [])
                .unwrap();
        }
        conn
    }

    #[test]
    fn precursors_collect_their_selections_and_fall_back_to_the_largest_peak() {
        let conn = tdf_with(
            &[
                "1, 500.30, 500.4, 500.25, 2, 700.5, 1000.0, 1",
                "2, 600.31, 600.4, NULL, NULL, 650.0, 800.0, 1",
                "3, 700.0, 700.0, 699.9, 3, 600.0, 10.0, 1",
            ],
            &[
                "4, 710, 735, 600.3, 2.0, 30.0, 2",
                "2, 690, 715, 500.3, 2.0, 28.0, 1",
                "5, 690, 715, 500.3, 2.0, 28.5, 1",
            ],
            &["1, 0.1, 100, 25.0, 26.0", "2, 0.2, 0, 25.0, 26.0"],
        );
        let p = read_pasef_precursors(&conn).unwrap();
        assert_eq!(p.iter().map(|p| p.id).collect::<Vec<_>>(), vec![1, 2]);
        assert_eq!(p[0].mz, 500.25);
        assert_eq!(p[0].charge, Some(2));
        assert_eq!(p[0].parent, Some(1));
        assert_eq!(
            p[0].selections
                .iter()
                .map(|s| s.frame_id)
                .collect::<Vec<_>>(),
            vec![2, 5]
        );
        assert_eq!(p[0].selections[0].scan_begin, 690);
        assert_eq!(p[0].selections[0].scan_end, 715);
        assert_eq!(p[1].mz, 600.31);
        assert_eq!(p[1].charge, None);
        assert_eq!(p[1].collision_energy, 30.0);

        let frames = read_frame_rows(&conn).unwrap();
        assert_eq!(frames[&1].num_peaks, 100);
        assert_eq!(frames[&2].num_peaks, 0);
        assert_eq!(frames[&2].t, [25.0, 26.0]);
    }

    #[test]
    fn a_tdf_without_pasef_tables_has_no_precursors() {
        let conn = rusqlite::Connection::open_in_memory().unwrap();
        conn.execute_batch("CREATE TABLE Frames (Id INTEGER, Time REAL, NumPeaks INTEGER);")
            .unwrap();
        assert!(read_pasef_precursors(&conn).unwrap().is_empty());
        let frames = read_frame_rows(&conn).unwrap();
        assert!(frames.is_empty());
    }

    #[test]
    fn dia_windows_expand_to_one_precursor_per_frame_and_window() {
        let conn = rusqlite::Connection::open_in_memory().unwrap();
        conn.execute_batch(
            "CREATE TABLE Frames (Id INTEGER, Time REAL, NumPeaks INTEGER, MsMsType INTEGER);
             INSERT INTO Frames VALUES (1, 0.1, 10, 0), (2, 0.2, 10, 9), (3, 0.3, 10, 9),
               (4, 0.4, 10, 0), (5, 0.5, 10, 9);
             CREATE TABLE DiaFrameMsMsInfo (Frame INTEGER, WindowGroup INTEGER);
             INSERT INTO DiaFrameMsMsInfo VALUES (2, 1), (3, 2), (5, 1);
             CREATE TABLE DiaFrameMsMsWindows (WindowGroup INTEGER, ScanNumBegin INTEGER, \
               ScanNumEnd INTEGER, IsolationMz REAL, IsolationWidth REAL, CollisionEnergy REAL);
             INSERT INTO DiaFrameMsMsWindows VALUES (1, 400, 600, 745.0, 45.0, 40.0),
               (1, 100, 300, 417.0, 35.0, 25.0), (2, 100, 300, 452.0, 37.0, 26.0);",
        )
        .unwrap();
        assert!(read_pasef_precursors(&conn).unwrap().is_empty());
        let w = read_dia_windows(&conn).unwrap();
        let got: Vec<(usize, f64, Option<usize>, usize)> = w
            .iter()
            .map(|p| {
                (
                    p.selections[0].frame_id,
                    p.mz,
                    p.parent,
                    p.selections[0].scan_begin,
                )
            })
            .collect();
        assert_eq!(
            got,
            vec![
                (2, 417.0, Some(1), 100),
                (2, 745.0, Some(1), 400),
                (3, 452.0, Some(1), 100),
                (5, 417.0, Some(4), 100),
                (5, 745.0, Some(4), 400),
            ]
        );
        assert!(w
            .iter()
            .all(|p| p.dia && p.charge.is_none() && p.scan_number.is_none()));
        assert_eq!(w[1].isolation_width, 45.0);
        assert_eq!(w[1].collision_energy, 40.0);
        let ids: Vec<usize> = w.iter().map(|p| p.id).collect();
        assert_eq!(ids, vec![1, 2, 3, 4, 5]);
    }

    #[test]
    fn mobility_slices_overlap_by_half_and_cover_the_window() {
        assert_eq!(mobility_slices(100, 400, 1), vec![(100, 400)]);
        assert_eq!(mobility_slices(100, 105, 4), vec![(100, 105)]);
        assert_eq!(mobility_slices(100, 400, 2), vec![(100, 300), (200, 400)]);
        let s = mobility_slices(34, 538, 8);
        assert_eq!(s.len(), 8);
        assert_eq!((s[0].0, s[7].1), (34, 538));
        for w in s.windows(2) {
            assert_eq!(w[1].0 - w[0].0, 56);
            assert!(w[1].0 < w[0].1);
        }
    }

    #[test]
    fn slice_dia_windows_splits_only_dia_precursors() {
        let conn = rusqlite::Connection::open_in_memory().unwrap();
        conn.execute_batch(
            "CREATE TABLE Frames (Id INTEGER, Time REAL, NumPeaks INTEGER, MsMsType INTEGER);
             INSERT INTO Frames VALUES (1, 0.1, 10, 0), (2, 0.2, 10, 9);
             CREATE TABLE DiaFrameMsMsInfo (Frame INTEGER, WindowGroup INTEGER);
             INSERT INTO DiaFrameMsMsInfo VALUES (2, 1);
             CREATE TABLE DiaFrameMsMsWindows (WindowGroup INTEGER, ScanNumBegin INTEGER, \
               ScanNumEnd INTEGER, IsolationMz REAL, IsolationWidth REAL, CollisionEnergy REAL);
             INSERT INTO DiaFrameMsMsWindows VALUES (1, 100, 400, 745.0, 45.0, 40.0),
               (1, 400, 700, 417.0, 35.0, 25.0);",
        )
        .unwrap();
        let w = slice_dia_windows(read_dia_windows(&conn).unwrap(), 2);
        let got: Vec<(usize, f64, usize, usize)> = w
            .iter()
            .map(|p| {
                (
                    p.id,
                    p.mz,
                    p.selections[0].scan_begin,
                    p.selections[0].scan_end,
                )
            })
            .collect();
        assert_eq!(
            got,
            vec![
                (1, 745.0, 100, 300),
                (2, 745.0, 200, 400),
                (3, 417.0, 400, 600),
                (4, 417.0, 500, 700),
            ]
        );
        assert!(w.iter().all(|p| p.dia && p.parent == Some(1)));
        let dda = vec![PasefPrecursor {
            id: 7,
            dia: false,
            ..w[0].clone()
        }];
        assert_eq!(slice_dia_windows(dda.clone(), 4), dda);
    }

    #[test]
    fn push_scans_takes_the_selected_scan_range_only() {
        let offsets = [0, 2, 3, 6];
        let tof = [10, 11, 20, 30, 31, 32];
        let intensity = [1, 2, 3, 4, 5, 6];
        let (mut t, mut i) = (Vec::new(), Vec::new());
        push_scans(&offsets, &tof, &intensity, 1, 3, &mut t, &mut i);
        assert_eq!(t, vec![20, 30, 31, 32]);
        assert_eq!(i, vec![3, 4, 5, 6]);
        let (mut t, mut i) = (Vec::new(), Vec::new());
        push_scans(&offsets, &tof, &intensity, 2, 9, &mut t, &mut i);
        assert_eq!(t, vec![30, 31, 32]);
        let (mut t, mut i) = (Vec::new(), Vec::new());
        push_scans(&offsets, &tof, &intensity, 7, 9, &mut t, &mut i);
        assert!(t.is_empty() && i.is_empty());
    }

    #[test]
    fn selections_of_one_ion_merge_into_one_centroid() {
        let tof = [1000, 1001, 1002, 1000, 1001, 1002, 5000];
        let intensity = [10, 40, 10, 10, 40, 10, 7];
        let cal = MzCalibration::linear(10.0, 1e-4, 1.0);
        let cond = FrameConditions {
            t: [0.0; 2],
            rt_s: 0.0,
        };
        let peaks = ms2_peaks(&tof, &intensity, &cal, cond, 2);
        assert_eq!(peaks.len(), 2);
        let expected = cal.mz(1001.0, conditions(cond, 1.0));
        assert!((peaks[0].0 - expected).abs() < 1e-9, "{peaks:?}");
        assert_eq!(peaks[0].1, 120.0);
        assert_eq!(peaks[1].1, 7.0);
        assert!(peaks[0].0 < peaks[1].0);
    }
}
