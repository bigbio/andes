//! PASEF MS2 spectra of a Bruker timsTOF `.d`, built from the raw frames.
//!
//! In DDA-PASEF the quadrupole isolates a precursor during a slice of the
//! mobility scans of an MS2 frame (`PasefFrameMsMsInfo`: frame, scan range,
//! isolation window, collision energy, precursor id), and a precursor below
//! its target intensity is selected again in later frames. One fragment
//! spectrum per precursor (`Precursors` row) is the sum of all its selections:
//! the TOF points of every selected scan of every frame are summed per TOF
//! index and centroided like the MS1 survey frames (see
//! [`crate::timstof_ms1::centroid_profile`]), then converted to m/z with the
//! run's calibration (the TOF axis and its calibration are shared by MS1 and
//! MS2 frames). Frames the TDF records with no peaks hold no data blob and are
//! skipped.
//!
//! The precursor m/z is the vendor `MonoisotopicMz`; where the acquisition
//! software did not determine the monoisotopic peak (and so no charge), the m/z
//! of the largest isotope peak (`LargestPeakMz`) is used with an unknown
//! charge, as mzML converters write it.

use std::collections::HashMap;

use crate::timstof_ms1::{
    centroid_profile, conditions, sum_over_mobility, FrameConditions, MzCalibration,
};

/// One quadrupole selection of a precursor: the scans `[scan_begin,
/// scan_end)` of a frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Selection {
    pub frame_id: usize,
    pub scan_begin: usize,
    pub scan_end: usize,
}

/// A DDA-PASEF precursor with every selection that fragmented it.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PasefPrecursor {
    pub id: usize,
    /// Monoisotopic m/z, else the largest isotope peak's m/z.
    pub mz: f64,
    /// `None` when the acquisition software assigned no charge.
    pub charge: Option<i64>,
    pub intensity: Option<f64>,
    /// Parent MS1 frame id.
    pub parent: Option<usize>,
    /// Mobility scan at the precursor's apex (fractional).
    pub scan_number: Option<f64>,
    pub isolation_width: f64,
    pub collision_energy: f64,
    pub selections: Vec<Selection>,
}

/// The `Frames` columns the MS2 builder needs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct FrameRow {
    /// Retention time, seconds.
    pub time: f64,
    pub num_peaks: u64,
    /// `T1`, `T2` temperatures (0 when the schema lacks them).
    pub t: [f64; 2],
}

/// Read the `Frames` rows keyed by frame id.
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

/// Read every DDA-PASEF precursor that has at least one selection, ordered by
/// precursor id, its selections in frame order. A `.d` without the PASEF
/// tables (not DDA-PASEF) has none.
pub(crate) fn read_pasef_precursors(
    conn: &rusqlite::Connection,
) -> Result<Vec<PasefPrecursor>, rusqlite::Error> {
    let has_table = |name: &str| -> Result<bool, rusqlite::Error> {
        conn.query_row(
            "SELECT count(*) FROM sqlite_master WHERE type = 'table' AND name = ?1",
            [name],
            |r| r.get::<_, i64>(0),
        )
        .map(|n| n > 0)
    };
    if !has_table("Precursors")? || !has_table("PasefFrameMsMsInfo")? {
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

/// The TOF points of the scans `[begin, end)` of a decoded frame, appended to
/// `tof` / `intensity`. `scan_offsets[k]` is the first point of scan `k`; scan
/// indices past the frame's last scan are clipped.
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

/// Centroided `(m/z, intensity)` peaks of a precursor's summed TOF points,
/// m/z ascending. `frame` holds the mean conditions of the selected frames and
/// `n_frames` their number: a centroid's m/z depends on the intensity its ions
/// reach the detector with in one frame, so the calibration's intensity term
/// is evaluated at the centroid's mean intensity per frame.
pub(crate) fn ms2_peaks(
    tof: &[u32],
    intensity: &[u32],
    calibration: &MzCalibration,
    frame: FrameConditions,
    n_frames: usize,
) -> Vec<(f64, f32)> {
    let per_frame = 1.0 / n_frames.max(1) as f64;
    let mut peaks: Vec<(f64, f32)> = centroid_profile(&sum_over_mobility(tof, intensity))
        .iter()
        .map(|c| {
            let cond = conditions(frame, c.intensity * per_frame);
            (calibration.mz(c.tof, cond), c.intensity as f32)
        })
        .filter(|&(mz, i)| mz.is_finite() && mz > 0.0 && i.is_finite() && i > 0.0)
        .collect();
    if peaks.windows(2).any(|w| w[0].0 > w[1].0) {
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
    }
    peaks
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
        // Precursor 3 has no selection and gets no spectrum.
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
        // No monoisotopic peak: largest peak, charge unknown.
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
    fn push_scans_takes_the_selected_scan_range_only() {
        // Three scans: [10, 11], [20], [30, 31, 32].
        let offsets = [0, 2, 3, 6];
        let tof = [10, 11, 20, 30, 31, 32];
        let intensity = [1, 2, 3, 4, 5, 6];
        let (mut t, mut i) = (Vec::new(), Vec::new());
        push_scans(&offsets, &tof, &intensity, 1, 3, &mut t, &mut i);
        assert_eq!(t, vec![20, 30, 31, 32]);
        assert_eq!(i, vec![3, 4, 5, 6]);
        // A range past the last scan is clipped to the frame.
        let (mut t, mut i) = (Vec::new(), Vec::new());
        push_scans(&offsets, &tof, &intensity, 2, 9, &mut t, &mut i);
        assert_eq!(t, vec![30, 31, 32]);
        let (mut t, mut i) = (Vec::new(), Vec::new());
        push_scans(&offsets, &tof, &intensity, 7, 9, &mut t, &mut i);
        assert!(t.is_empty() && i.is_empty());
    }

    #[test]
    fn selections_of_one_ion_merge_into_one_centroid() {
        // The same fragment seen in two frames at slightly different TOF
        // bins: one centroid at the intensity-weighted position.
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
