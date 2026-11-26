use std::collections::HashMap;

use model::stats::{median_mad, Even, MAD_TO_SD};
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct TofCentroid {
    pub tof: f64,
    pub intensity: f64,
}

pub(crate) fn sum_over_mobility(tof: &[u32], intensity: &[u32]) -> Vec<(u32, u64)> {
    let n = tof.len().min(intensity.len());
    let mut pairs: Vec<(u32, u32)> = (0..n).map(|i| (tof[i], intensity[i])).collect();
    pairs.sort_unstable_by_key(|p| p.0);
    let mut out: Vec<(u32, u64)> = Vec::with_capacity(pairs.len() / 2 + 1);
    for (t, i) in pairs {
        match out.last_mut() {
            Some(last) if last.0 == t => last.1 += u64::from(i),
            _ => out.push((t, u64::from(i))),
        }
    }
    out
}

const MAX_TOF_GAP: u32 = 2;

const VALLEY_SPLIT_FRACTION: f64 = 0.5;

const SMOOTHING: [f64; 5] = [1.0 / 16.0, 4.0 / 16.0, 6.0 / 16.0, 4.0 / 16.0, 1.0 / 16.0];

const CORE_FRACTION: f64 = 0.5;

pub(crate) fn centroid_profile(profile: &[(u32, u64)]) -> Vec<TofCentroid> {
    centroid_profile_with_bounds(profile, None)
}

pub(crate) fn centroid_profile_with_bounds(
    profile: &[(u32, u64)],
    mut bounds: Option<&mut Vec<(u32, u32)>>,
) -> Vec<TofCentroid> {
    let mut out = Vec::new();
    let mut start = 0;
    while start < profile.len() {
        let mut end = start + 1;
        while end < profile.len() && profile[end].0 - profile[end - 1].0 <= MAX_TOF_GAP {
            end += 1;
        }
        centroid_run(&profile[start..end], &mut out, bounds.as_deref_mut());
        start = end;
    }
    out
}

fn centroid_run(
    run: &[(u32, u64)],
    out: &mut Vec<TofCentroid>,
    mut bounds: Option<&mut Vec<(u32, u32)>>,
) {
    let first = run[0].0;
    let n = (run[run.len() - 1].0 - first) as usize + 1;
    let mut raw = vec![0.0f64; n];
    for &(t, i) in run {
        raw[(t - first) as usize] = i as f64;
    }
    let half = SMOOTHING.len() / 2;
    let smooth: Vec<f64> = (0..n)
        .map(|i| {
            SMOOTHING
                .iter()
                .enumerate()
                .filter_map(|(k, w)| {
                    let j = (i + k).checked_sub(half)?;
                    raw.get(j).map(|v| v * w)
                })
                .sum()
        })
        .collect();
    let mut segments: Vec<(usize, usize, f64, f64)> = Vec::new();
    let mut lo = 0;
    let mut valley = 0.0;
    for i in 1..n.saturating_sub(1) {
        if smooth[i] < smooth[i - 1] && smooth[i] <= smooth[i + 1] {
            let max = smooth[lo..i].iter().copied().fold(0.0, f64::max);
            segments.push((lo, i, max, valley));
            lo = i;
            valley = smooth[i];
        }
    }
    let max = smooth[lo..].iter().copied().fold(0.0, f64::max);
    segments.push((lo, n, max, valley));
    let mut merged: Vec<(usize, usize, f64)> = Vec::with_capacity(segments.len());
    for (lo, hi, max, valley) in segments {
        match merged.last_mut() {
            Some(last) if valley >= VALLEY_SPLIT_FRACTION * last.2.min(max) => {
                last.1 = hi;
                last.2 = last.2.max(max);
            }
            _ => merged.push((lo, hi, max)),
        }
    }
    for (lo, hi, max) in merged {
        let total: f64 = raw[lo..hi].iter().sum();
        if total <= 0.0 {
            continue;
        }
        let (mut w, mut wt) = (0.0, 0.0);
        for i in lo..hi {
            if smooth[i] >= CORE_FRACTION * max {
                w += raw[i];
                wt += raw[i] * (first as f64 + i as f64);
            }
        }
        if w <= 0.0 {
            w = total;
            wt = (lo..hi).map(|i| raw[i] * (first as f64 + i as f64)).sum();
        }
        out.push(TofCentroid {
            tof: wt / w,
            intensity: total,
        });
        if let Some(b) = bounds.as_deref_mut() {
            b.push((first + lo as u32, first + hi as u32 - 1));
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct MobilityCentroid {
    pub tof: f64,
    pub intensity: f64,
    pub scan: f64,
}

const MOBILITY_SMOOTHING_SCANS: f64 = 3.0;

const MOBILITY_GAP_SCANS: u32 = 9;

const MIN_COMPONENT_POINTS: usize = 4;

pub(crate) fn centroid_frame(
    scan_offsets: &[usize],
    tof: &[u32],
    intensity: &[u32],
    resolve: bool,
) -> (Vec<TofCentroid>, Vec<MobilityCentroid>) {
    let summed_profile = sum_over_mobility(tof, intensity);
    if !resolve {
        return (centroid_profile(&summed_profile), Vec::new());
    }
    let mut bounds = Vec::new();
    let summed = centroid_profile_with_bounds(&summed_profile, Some(&mut bounds));
    let n = tof.len().min(intensity.len());
    let mut points: Vec<(u32, u32, u32)> = Vec::with_capacity(n);
    for s in 0..scan_offsets.len().saturating_sub(1) {
        let (a, b) = (scan_offsets[s].min(n), scan_offsets[s + 1].min(n));
        for i in a..b {
            points.push((tof[i], s as u32, intensity[i]));
        }
    }
    points.sort_unstable_by_key(|p| (p.0, p.1));
    let mut out = Vec::with_capacity(summed.len());
    let mut at = 0;
    let mut buf = Vec::new();
    for &(lo, hi) in &bounds {
        while at < points.len() && points[at].0 < lo {
            at += 1;
        }
        let start = at;
        while at < points.len() && points[at].0 <= hi {
            at += 1;
        }
        buf.clear();
        buf.extend_from_slice(&points[start..at]);
        split_by_mobility(&mut buf, &mut out);
    }
    (summed, out)
}

#[derive(Debug, Clone)]
struct Component {
    lo: u32,
    hi: u32,
    points: usize,
    total: f64,
    max: f64,
    valley: f64,
}

fn split_by_mobility(pts: &mut [(u32, u32, u32)], out: &mut Vec<MobilityCentroid>) {
    if pts.is_empty() {
        return;
    }
    pts.sort_unstable_by_key(|p| (p.1, p.0));
    let mut by_scan: Vec<(u32, f64, usize)> = Vec::new();
    for &(_, s, i) in pts.iter() {
        match by_scan.last_mut() {
            Some(last) if last.0 == s => {
                last.1 += i as f64;
                last.2 += 1;
            }
            _ => by_scan.push((s, i as f64, 1)),
        }
    }
    let half = (3.0 * MOBILITY_SMOOTHING_SCANS).ceil() as i64;
    let kernel: Vec<f64> = (-half..=half)
        .map(|d| (-0.5 * (d as f64 / MOBILITY_SMOOTHING_SCANS).powi(2)).exp())
        .collect();
    let mut parts: Vec<Component> = Vec::new();
    let mut start = 0;
    while start < by_scan.len() {
        let mut end = start + 1;
        while end < by_scan.len() && by_scan[end].0 - by_scan[end - 1].0 <= MOBILITY_GAP_SCANS {
            end += 1;
        }
        mobility_run(&by_scan[start..end], &kernel, half, &mut parts);
        start = end;
    }
    let mut merged: Vec<Component> = Vec::with_capacity(parts.len());
    for c in parts {
        match merged.last_mut() {
            Some(last)
                if c.valley > 0.0 && c.valley >= VALLEY_SPLIT_FRACTION * last.max.min(c.max) =>
            {
                last.hi = c.hi;
                last.points += c.points;
                last.total += c.total;
                last.max = last.max.max(c.max);
            }
            _ => merged.push(c),
        }
    }
    for c in merged.iter().filter(|c| c.points >= MIN_COMPONENT_POINTS) {
        emit(pts, c.lo, c.hi, out);
    }
}

fn mobility_run(run: &[(u32, f64, usize)], kernel: &[f64], half: i64, parts: &mut Vec<Component>) {
    let first = run[0].0;
    let n = (run[run.len() - 1].0 - first) as usize + 1;
    let mut raw = vec![0.0f64; n];
    let mut count = vec![0usize; n];
    for &(s, i, c) in run {
        raw[(s - first) as usize] = i;
        count[(s - first) as usize] = c;
    }
    let smooth: Vec<f64> = (0..n as i64)
        .map(|i| {
            kernel
                .iter()
                .enumerate()
                .filter_map(|(k, w)| {
                    let j = i + k as i64 - half;
                    (0..n as i64).contains(&j).then(|| raw[j as usize] * w)
                })
                .sum()
        })
        .collect();
    let mut lo = 0;
    let mut valley = 0.0;
    let push = |lo: usize, hi: usize, valley: f64, parts: &mut Vec<Component>| {
        let own: f64 = raw[lo..hi].iter().sum();
        let points: usize = count[lo..hi].iter().sum();
        let max = smooth[lo..hi].iter().copied().fold(0.0, f64::max);
        let a = (lo..hi).find(|&k| count[k] > 0).unwrap_or(lo);
        let b = (lo..hi).rev().find(|&k| count[k] > 0).unwrap_or(hi - 1);
        parts.push(Component {
            lo: first + a as u32,
            hi: first + b as u32,
            points,
            total: own,
            max,
            valley,
        });
    };
    for i in 1..n.saturating_sub(1) {
        if smooth[i] < smooth[i - 1] && smooth[i] <= smooth[i + 1] {
            push(lo, i, valley, parts);
            lo = i;
            valley = smooth[i];
        }
    }
    push(lo, n, valley, parts);
}

fn emit(pts: &[(u32, u32, u32)], lo: u32, hi: u32, out: &mut Vec<MobilityCentroid>) {
    let a = pts.partition_point(|p| p.1 < lo);
    let b = pts.partition_point(|p| p.1 <= hi);
    let own = &pts[a..b];
    if own.is_empty() {
        return;
    }
    let (mut sum, mut sum_scan) = (0.0f64, 0.0f64);
    let (mut tmin, mut tmax) = (u32::MAX, 0u32);
    for &(t, s, i) in own {
        sum += i as f64;
        sum_scan += i as f64 * s as f64;
        tmin = tmin.min(t);
        tmax = tmax.max(t);
    }
    if sum <= 0.0 {
        return;
    }
    let width = (tmax - tmin) as usize + 1;
    let mut raw = vec![0.0f64; width];
    for &(t, _, i) in own {
        raw[(t - tmin) as usize] += i as f64;
    }
    let half = SMOOTHING.len() / 2;
    let smooth: Vec<f64> = (0..width)
        .map(|i| {
            SMOOTHING
                .iter()
                .enumerate()
                .filter_map(|(k, w)| {
                    let j = (i + k).checked_sub(half)?;
                    raw.get(j).map(|v| v * w)
                })
                .sum()
        })
        .collect();
    let max = smooth.iter().copied().fold(0.0, f64::max);
    let (mut w, mut wt) = (0.0, 0.0);
    for (i, &r) in raw.iter().enumerate() {
        if smooth[i] >= CORE_FRACTION * max {
            w += r;
            wt += r * (tmin as f64 + i as f64);
        }
    }
    if w <= 0.0 {
        w = sum;
        wt = raw
            .iter()
            .enumerate()
            .map(|(i, r)| r * (tmin as f64 + i as f64))
            .sum();
    }
    out.push(MobilityCentroid {
        tof: wt / w,
        intensity: sum,
        scan: sum_scan / sum,
    });
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct FrameConditions {
    pub t: [f64; 2],
    pub rt_s: f64,
}

pub(crate) type Conditions = [f64; 5];

pub(crate) fn conditions(f: FrameConditions, intensity: f64) -> Conditions {
    let h = f.rt_s / 3600.0;
    [f.t[0], f.t[1], intensity.max(1.0).ln(), h, h * h]
}

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

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct MzCalibration {
    scale: f64,
    poly: [f64; 4],
    cond: [f64; 5],
    c_ref: Conditions,
    tdf: Option<TdfMzModel>,
}

impl MzCalibration {
    pub fn linear(intercept: f64, slope: f64, scale: f64) -> Self {
        Self {
            scale,
            poly: [intercept, slope * scale, 0.0, 0.0],
            cond: [0.0; 5],
            c_ref: [0.0; 5],
            tdf: None,
        }
    }

    pub fn tdf(model: TdfMzModel, scale: f64) -> Self {
        Self {
            tdf: Some(model),
            ..Self::linear(0.0, 0.0, scale)
        }
    }

    pub fn is_tdf(&self) -> bool {
        self.tdf.is_some()
    }

    /// This calibration with a frame's own TDF row in place of the run's main row. A linear
    /// or fitted calibration is returned unchanged.
    pub fn for_frame(&self, model: Option<TdfMzModel>) -> Self {
        match (self.tdf, model) {
            (Some(_), Some(m)) => Self {
                tdf: Some(m),
                ..self.clone()
            },
            _ => self.clone(),
        }
    }

    pub fn mz(&self, tof: f64, c: Conditions) -> f64 {
        if let Some(model) = &self.tdf {
            return model.mz(tof, [c[0], c[1]]);
        }
        let x = tof / self.scale;
        let p = &self.poly;
        let terms: f64 = (0..5).map(|k| self.cond[k] * (c[k] - self.c_ref[k])).sum();
        let root = p[0] + x * (p[1] + x * (p[2] + x * p[3])) + terms;
        root * root
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CalibrationPoint {
    pub tof: f64,
    pub mz: f64,
    pub t: Conditions,
    /// The base calibration's m/z for this point, from its own frame's MzCalibration row.
    /// Base residuals use it, so a point on another row is not charged with the row
    /// difference. `NaN` means "evaluate the run base" (points built without a frame).
    pub base_mz: f64,
}

impl CalibrationPoint {
    fn base_ppm(&self, base: &MzCalibration) -> f64 {
        let m = if self.base_mz.is_finite() {
            self.base_mz
        } else {
            base.mz(self.tof, self.t)
        };
        (m - self.mz) / self.mz * 1e6
    }
}

pub(crate) const WIDE_WINDOW_PPM: f64 = 200.0;
pub(crate) const NARROW_WINDOW_PPM: f64 = 15.0;
pub(crate) const MIN_CALIBRATION_POINTS: usize = 200;
const OUTLIER_SIGMAS: f64 = 4.0;

const MATCH_MARGIN_PPM: f64 = 50.0;

/// Match each vendor m/z to the most intense centroid inside the window. `calibration` places
/// the centroids, and `base` (the run base with this frame's MzCalibration row) gives each
/// point's `base_mz`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn match_precursors(
    centroids: &[TofCentroid],
    vendor_mz: &[f64],
    frame: FrameConditions,
    calibration: &MzCalibration,
    base: &MzCalibration,
    offset: &OffsetCurve,
    window_ppm: f64,
    out: &mut Vec<CalibrationPoint>,
) {
    let nominal = conditions(frame, 1.0);
    for &mz in vendor_mz {
        if !(mz.is_finite() && mz > 0.0) {
            continue;
        }
        let offset_ppm = offset.at(mz);
        let reach = (offset_ppm.abs() + window_ppm + MATCH_MARGIN_PPM) * 1e-6;
        let i = centroids.partition_point(|c| calibration.mz(c.tof, nominal) < mz * (1.0 - reach));
        let j = centroids.partition_point(|c| calibration.mz(c.tof, nominal) <= mz * (1.0 + reach));
        let best = centroids[i..j.max(i)]
            .iter()
            .map(|c| (c, conditions(frame, c.intensity)))
            .filter(|(c, cond)| {
                let ppm = (calibration.mz(c.tof, *cond) - mz) / mz * 1e6;
                (ppm - offset_ppm).abs() <= window_ppm
            })
            .max_by(|a, b| a.0.intensity.total_cmp(&b.0.intensity));
        if let Some((c, cond)) = best {
            out.push(CalibrationPoint {
                tof: c.tof,
                mz,
                t: cond,
                base_mz: base.mz(c.tof, cond),
            });
        }
    }
}

pub(crate) const REFINE_ROUNDS: usize = 2;
pub(crate) const REFINE_SIGMAS: f64 = 3.0;
pub(crate) const MIN_REFINE_WINDOW_PPM: f64 = 8.0;

const OFFSET_BIN_PPM: f64 = 1.0;

pub(crate) fn mode_offset_ppm(points: &[CalibrationPoint], base: &MzCalibration) -> Option<f64> {
    if points.is_empty() {
        return None;
    }
    let half = (WIDE_WINDOW_PPM / OFFSET_BIN_PPM).ceil() as i64;
    let mut hist = vec![0.0f64; (2 * half + 1) as usize];
    for p in points {
        let ppm = p.base_ppm(base);
        let k = (ppm / OFFSET_BIN_PPM).round() as i64 + half;
        if let Some(h) = usize::try_from(k).ok().and_then(|k| hist.get_mut(k)) {
            *h += 1.0;
        }
    }
    let smooth = |k: usize| {
        hist[k.saturating_sub(1)..(k + 2).min(hist.len())]
            .iter()
            .sum::<f64>()
    };
    let k = (0..hist.len()).max_by(|&a, &b| smooth(a).total_cmp(&smooth(b)))?;
    let centre = (k as i64 - half) as f64 * OFFSET_BIN_PPM;
    let near: Vec<f64> = points
        .iter()
        .map(|p| p.base_ppm(base))
        .filter(|ppm| (ppm - centre).abs() <= 1.5 * OFFSET_BIN_PPM)
        .collect();
    Some(near.iter().sum::<f64>() / near.len().max(1) as f64)
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct OffsetCurve {
    knots: Vec<(f64, f64)>,
}

impl OffsetCurve {
    pub fn constant(ppm: f64) -> Self {
        Self {
            knots: vec![(0.0, ppm)],
        }
    }

    pub fn at(&self, mz: f64) -> f64 {
        let k = &self.knots;
        let i = k.partition_point(|&(x, _)| x <= mz);
        if i == 0 {
            return k[0].1;
        }
        if i == k.len() {
            return k[k.len() - 1].1;
        }
        let ((x0, y0), (x1, y1)) = (k[i - 1], k[i]);
        y0 + (y1 - y0) * (mz - x0) / (x1 - x0)
    }
}

const OFFSET_RANGE_POINTS: usize = 400;
const OFFSET_MAX_RANGES: usize = 10;

pub(crate) fn offset_curve(
    points: &[CalibrationPoint],
    base: &MzCalibration,
) -> Option<OffsetCurve> {
    let mut sorted = points.to_vec();
    sorted.sort_by(|a, b| a.mz.total_cmp(&b.mz));
    let ranges = (sorted.len() / OFFSET_RANGE_POINTS).clamp(1, OFFSET_MAX_RANGES);
    let mut knots = Vec::with_capacity(ranges);
    for r in 0..ranges {
        let part = &sorted[r * sorted.len() / ranges..(r + 1) * sorted.len() / ranges];
        if let Some(ppm) = mode_offset_ppm(part, base) {
            knots.push((part[part.len() / 2].mz, ppm));
        }
    }
    (!knots.is_empty()).then_some(OffsetCurve { knots })
}

fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        let pivot = (col..n).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
        if a[pivot][col].abs() < 1e-300 {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        for row in col + 1..n {
            let f = a[row][col] / a[col][col];
            for k in col..n {
                a[row][k] -= f * a[col][k];
            }
            b[row] -= f * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let s: f64 = (row + 1..n).map(|k| a[row][k] * x[k]).sum();
        x[row] = (b[row] - s) / a[row][row];
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

const MIN_CONDITION_SPREAD: f64 = 1e-3;

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CalibrationFit {
    pub calibration: MzCalibration,
    pub n_used: usize,
    pub median_ppm: f64,
    pub sigma_ppm: f64,
    pub base_median_ppm: f64,
    pub base_sigma_ppm: f64,
}

fn residual_stats(ppm: &[f64]) -> (f64, f64) {
    let (med, mad) = median_mad(ppm, Even::Upper).unwrap_or((0.0, 0.0));
    (med, MAD_TO_SD * mad)
}

/// Fraction of the base model's residual (offset and spread together) that a precursor-peak
/// fit must stay under to replace an instrument (TDF) calibration.
const FIT_OVER_TDF: f64 = 0.8;

/// Whether a precursor-peak fit is clearly closer to the vendor m/z than its base model. The
/// vendor m/z come from the same instrument calibration, so on a TDF base the fit should only
/// win when the calibration does not describe the frames (for example an unusable row).
pub(crate) fn fit_beats_base(fit: &CalibrationFit) -> bool {
    let fit_err = fit.median_ppm.hypot(fit.sigma_ppm);
    let base_err = fit.base_median_ppm.hypot(fit.base_sigma_ppm);
    fit_err < FIT_OVER_TDF * base_err
}

pub(crate) fn fit_calibration(
    points: &[CalibrationPoint],
    base: &MzCalibration,
) -> Option<CalibrationFit> {
    if points.len() < MIN_CALIBRATION_POINTS {
        return None;
    }
    let scale = base.scale;
    let mean_t = |k: usize| points.iter().map(|p| p.t[k]).sum::<f64>() / points.len() as f64;
    let t_ref: Conditions = [0, 1, 2, 3, 4].map(mean_t);
    let spread = |k: usize| {
        (points
            .iter()
            .map(|p| (p.t[k] - t_ref[k]).powi(2))
            .sum::<f64>()
            / points.len() as f64)
            .sqrt()
    };
    let use_temp = [0, 1, 2, 3, 4].map(|k| spread(k) > MIN_CONDITION_SPREAD);
    let features = |p: &CalibrationPoint| -> Vec<f64> {
        let x = p.tof / scale;
        let mut f = vec![1.0, x, x * x, x * x * x];
        for k in 0..5 {
            if use_temp[k] {
                f.push(p.t[k] - t_ref[k]);
            }
        }
        f
    };
    let base_ppm: Vec<f64> = points.iter().map(|p| p.base_ppm(base)).collect();
    let (base_median_ppm, base_sigma_ppm) = residual_stats(&base_ppm);
    let mut keep: Vec<bool> = vec![true; points.len()];
    let mut result: Option<CalibrationFit> = None;
    for _ in 0..3 {
        let used: Vec<&CalibrationPoint> = points
            .iter()
            .zip(&keep)
            .filter(|(_, &k)| k)
            .map(|(p, _)| p)
            .collect();
        if used.len() < MIN_CALIBRATION_POINTS {
            break;
        }
        let n = features(used[0]).len();
        let mut ata = vec![vec![0.0; n]; n];
        let mut atb = vec![0.0; n];
        for p in &used {
            let f = features(p);
            let y = p.mz.sqrt();
            for i in 0..n {
                atb[i] += f[i] * y;
                for j in 0..n {
                    ata[i][j] += f[i] * f[j];
                }
            }
        }
        let coef = solve(ata, atb)?;
        let mut temp = [0.0; 5];
        let mut k = 4;
        for (t, &on) in temp.iter_mut().zip(&use_temp) {
            if on {
                *t = coef[k];
                k += 1;
            }
        }
        let calibration = MzCalibration {
            scale,
            poly: [coef[0], coef[1], coef[2], coef[3]],
            cond: temp,
            c_ref: t_ref,
            tdf: None,
        };
        let ppm: Vec<f64> = points
            .iter()
            .map(|p| (calibration.mz(p.tof, p.t) - p.mz) / p.mz * 1e6)
            .collect();
        let used_ppm: Vec<f64> = ppm
            .iter()
            .zip(&keep)
            .filter(|(_, &k)| k)
            .map(|(r, _)| *r)
            .collect();
        let (median_ppm, sigma_ppm) = residual_stats(&used_ppm);
        result = Some(CalibrationFit {
            calibration,
            n_used: used_ppm.len(),
            median_ppm,
            sigma_ppm,
            base_median_ppm,
            base_sigma_ppm,
        });
        let limit = OUTLIER_SIGMAS * sigma_ppm.max(1e-3);
        keep = ppm
            .iter()
            .map(|r| (r - median_ppm).abs() <= limit)
            .collect();
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn summing_merges_equal_tof_indices_across_scans() {
        let tof = [10, 12, 10, 11, 12];
        let inten = [1, 2, 3, 4, 5];
        assert_eq!(
            sum_over_mobility(&tof, &inten),
            vec![(10, 4), (11, 4), (12, 7)]
        );
    }

    #[test]
    fn one_peak_is_one_weighted_centroid() {
        let profile = vec![(100, 1), (101, 3), (102, 5), (103, 3), (104, 1)];
        let c = centroid_profile(&profile);
        assert_eq!(c.len(), 1);
        assert!((c[0].tof - 102.0).abs() < 1e-12);
        assert_eq!(c[0].intensity, 13.0);
    }

    fn gaussian(
        lo: u32,
        hi: u32,
        centre: f64,
        sigma: f64,
        height: f64,
        jitter: f64,
    ) -> Vec<(u32, u64)> {
        (lo..hi)
            .map(|t| {
                let x = (t as f64 - centre) / sigma;
                let sign = if t % 2 == 0 { 1.0 } else { -1.0 };
                let v = height * (-0.5 * x * x).exp() * (1.0 + sign * jitter);
                (t, v.round() as u64)
            })
            .filter(|p| p.1 > 0)
            .collect()
    }

    fn add(a: &[(u32, u64)], b: &[(u32, u64)]) -> Vec<(u32, u64)> {
        let mut m = std::collections::BTreeMap::new();
        for &(t, i) in a.iter().chain(b) {
            *m.entry(t).or_insert(0) += i;
        }
        m.into_iter().collect()
    }

    #[test]
    fn bin_to_bin_noise_does_not_split_one_peak() {
        let p = gaussian(990, 1011, 1000.3, 2.0, 40.0, 0.3);
        let c = centroid_profile(&p);
        assert_eq!(c.len(), 1, "{c:?}");
        assert!((c[0].tof - 1000.3).abs() < 0.25, "{c:?}");
    }

    #[test]
    fn resolved_peaks_split_and_an_unresolved_shoulder_does_not() {
        let p = add(
            &gaussian(990, 1011, 1000.0, 2.0, 100.0, 0.0),
            &gaussian(1000, 1021, 1010.0, 2.0, 60.0, 0.0),
        );
        let c = centroid_profile(&p);
        assert_eq!(c.len(), 2, "{c:?}");
        assert!(
            (c[0].tof - 1000.0).abs() < 0.2 && (c[1].tof - 1010.0).abs() < 0.2,
            "{c:?}"
        );
        let q = add(
            &gaussian(990, 1011, 1000.0, 2.0, 100.0, 0.0),
            &gaussian(993, 1014, 1003.0, 2.0, 60.0, 0.0),
        );
        assert_eq!(centroid_profile(&q).len(), 1);
    }

    fn frame(n_scans: usize, mut pts: Vec<(u32, u32, u32)>) -> (Vec<usize>, Vec<u32>, Vec<u32>) {
        pts.sort_by_key(|p| (p.1, p.0));
        let mut offsets = vec![0usize];
        let (mut tof, mut inten) = (Vec::new(), Vec::new());
        let mut k = 0;
        for s in 0..n_scans as u32 {
            while k < pts.len() && pts[k].1 == s {
                tof.push(pts[k].0);
                inten.push(pts[k].2);
                k += 1;
            }
            offsets.push(tof.len());
        }
        (offsets, tof, inten)
    }

    fn ion(tof: f64, scan: f64, height: f64) -> Vec<(u32, u32, u32)> {
        let mut out = Vec::new();
        for s in (scan as i64 - 30)..(scan as i64 + 30) {
            for t in (tof as i64 - 6)..(tof as i64 + 7) {
                let ds = (s as f64 - scan) / 8.0;
                let dt = (t as f64 - tof) / 1.5;
                let v = (height * (-0.5 * (ds * ds + dt * dt)).exp()).round() as u32;
                if v > 0 {
                    out.push((t as u32, s as u32, v));
                }
            }
        }
        out
    }

    fn merge(parts: &[Vec<(u32, u32, u32)>]) -> Vec<(u32, u32, u32)> {
        let mut m = std::collections::BTreeMap::new();
        for &(t, s, i) in parts.iter().flatten() {
            *m.entry((t, s)).or_insert(0) += i;
        }
        m.into_iter().map(|((t, s), i)| (t, s, i)).collect()
    }

    #[test]
    fn without_resolution_the_frame_is_summed_over_mobility() {
        let (off, tof, inten) = frame(600, merge(&[ion(1000.0, 300.0, 50.0)]));
        let (summed, resolved) = centroid_frame(&off, &tof, &inten, false);
        assert_eq!(summed, centroid_profile(&sum_over_mobility(&tof, &inten)));
        assert!(resolved.is_empty());
    }

    #[test]
    fn ions_at_one_tof_are_split_by_mobility() {
        let pts = merge(&[ion(1000.0, 300.0, 20.0), ion(1000.4, 340.0, 80.0)]);
        let (off, tof, inten) = frame(600, pts);
        let (summed, resolved) = centroid_frame(&off, &tof, &inten, true);
        assert_eq!(summed, centroid_profile(&sum_over_mobility(&tof, &inten)));
        assert_eq!(summed.len(), 1, "{summed:?}");
        assert_eq!(resolved.len(), 2, "{resolved:?}");
        assert!((resolved[0].scan - 300.0).abs() < 1.5, "{resolved:?}");
        assert!((resolved[1].scan - 340.0).abs() < 1.5, "{resolved:?}");
        assert!((resolved[0].tof - 1000.0).abs() < 0.3, "{resolved:?}");
        assert!((resolved[1].tof - 1000.4).abs() < 0.3, "{resolved:?}");
        let total: f64 = resolved.iter().map(|c| c.intensity).sum();
        assert!(total <= summed[0].intensity && total > 0.99 * summed[0].intensity);
        assert!(resolved[1].intensity > 3.0 * resolved[0].intensity);
    }

    #[test]
    fn scattered_counts_do_not_make_components_of_their_own() {
        let mut noise = Vec::new();
        for s in [20u32, 80, 150, 450, 520] {
            noise.push((1000u32, s, 1u32));
        }
        let pts = merge(&[ion(1000.0, 300.0, 30.0), noise]);
        let (off, tof, inten) = frame(600, pts);
        let (summed, resolved) = centroid_frame(&off, &tof, &inten, true);
        assert_eq!(summed.len(), 1);
        assert_eq!(resolved.len(), 1, "{resolved:?}");
        assert_eq!(resolved[0].intensity, summed[0].intensity - 5.0);
        assert!((resolved[0].scan - 300.0).abs() < 1.0, "{resolved:?}");
        let sparse: Vec<(u32, u32, u32)> = [10u32, 90, 200, 333]
            .iter()
            .map(|&s| (2000u32, s, 1u32))
            .collect();
        let (off, tof, inten) = frame(600, sparse);
        let (summed, resolved) = centroid_frame(&off, &tof, &inten, true);
        assert_eq!(summed.len(), 1);
        assert!(resolved.is_empty(), "{resolved:?}");
    }

    #[test]
    fn a_gap_wider_than_the_limit_separates_peaks() {
        let profile = vec![(1, 5), (2, 5), (3 + MAX_TOF_GAP + 1, 5)];
        assert_eq!(centroid_profile(&profile).len(), 2);
        let close = vec![(1, 5), (1 + MAX_TOF_GAP, 5)];
        assert_eq!(centroid_profile(&close).len(), 1);
    }

    fn truth(tof: f64, t: Conditions) -> f64 {
        let x = tof / 400_000.0;
        let root = 10.0 + 31.0 * x + 0.02 * x * x - 0.01 * x * x * x
            + 2e-4 * (t[1] - 27.0)
            + 4e-5 * (t[2] - 5.0)
            + 3e-4 * t[3]
            - 1e-4 * t[4];
        root * root
    }

    #[test]
    fn the_fit_recovers_a_cubic_temperature_dependent_calibration() {
        let base = MzCalibration::linear(10.0, 31.0 / 400_000.0, 400_000.0);
        let mut points = Vec::new();
        for i in 0..2000 {
            let tof = 20_000.0 + (i as f64) * 180.0;
            let t = conditions(
                FrameConditions {
                    t: [25.7, 27.0 + 0.3 * ((i % 7) as f64 / 6.0)],
                    rt_s: ((i * 7919) % 2000) as f64 * 3.6,
                },
                (3.0 + 0.5 * (i % 11) as f64).exp(),
            );
            points.push(CalibrationPoint {
                tof,
                mz: truth(tof, t),
                t,
                base_mz: f64::NAN,
            });
        }
        let fit = fit_calibration(&points, &base).expect("enough pairs");
        assert!(fit.sigma_ppm < 0.01, "{fit:?}");
        assert!(fit.median_ppm.abs() < 0.01, "{fit:?}");
        assert!(fit.base_sigma_ppm > fit.sigma_ppm);
        for &(tof, t) in &[
            (50_000.0, [25.7, 27.1, 4.0, 0.3, 0.09]),
            (300_000.0, [25.7, 27.25, 7.5, 1.7, 2.89]),
        ] {
            let ppm = (fit.calibration.mz(tof, t) - truth(tof, t)) / truth(tof, t) * 1e6;
            assert!(ppm.abs() < 0.01, "{ppm}");
        }
    }

    #[test]
    fn outliers_are_left_out_of_the_refit() {
        let base = MzCalibration::linear(10.0, 31.0 / 400_000.0, 400_000.0);
        let mut points: Vec<CalibrationPoint> = (0..1000)
            .map(|i| {
                let tof = 20_000.0 + (i as f64) * 350.0;
                let t = [25.7, 27.0, 5.0, 1.0, 1.0];
                CalibrationPoint {
                    tof,
                    mz: truth(tof, t),
                    t,
                    base_mz: f64::NAN,
                }
            })
            .collect();
        for p in points.iter_mut().step_by(20) {
            p.mz *= 1.0 + 30e-6;
        }
        let fit = fit_calibration(&points, &base).unwrap();
        assert!(fit.median_ppm.abs() < 0.05, "{fit:?}");
        assert!(fit.n_used <= 950 + 1, "{fit:?}");
    }

    #[test]
    fn too_few_pairs_give_no_fit() {
        let base = MzCalibration::linear(10.0, 1e-4, 1.0);
        let points = vec![
            CalibrationPoint {
                tof: 1.0,
                mz: 100.0,
                t: [0.0; 5],
                base_mz: f64::NAN,
            };
            MIN_CALIBRATION_POINTS - 1
        ];
        assert!(fit_calibration(&points, &base).is_none());
    }

    #[test]
    fn matching_picks_the_most_intense_centroid_inside_the_window() {
        let base = MzCalibration::linear(10.0, 31.0 / 400_000.0, 400_000.0);
        let t = FrameConditions {
            t: [0.0; 2],
            rt_s: 0.0,
        };
        let target = base.mz(100_000.0, [0.0; 5]);
        let c = |tof: f64, intensity: f64| TofCentroid { tof, intensity };
        let centroids = vec![
            c(99_990.0, 50.0),
            c(100_000.4, 1.0),
            c(100_001.0, 20.0),
            c(100_200.0, 99.0),
        ];
        let mut out = Vec::new();
        match_precursors(
            &centroids,
            &[target],
            t,
            &base,
            &base,
            &OffsetCurve::constant(0.0),
            20.0,
            &mut out,
        );
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].tof, 100_001.0);
        out.clear();
        match_precursors(
            &centroids,
            &[target],
            t,
            &base,
            &base,
            &OffsetCurve::constant(500.0),
            0.001,
            &mut out,
        );
        assert!(out.is_empty());
    }

    #[test]
    fn the_offset_is_the_mode_not_the_median() {
        let base = MzCalibration::linear(10.0, 31.0 / 400_000.0, 400_000.0);
        let t = [0.0; 5];
        let mut points = Vec::new();
        for i in 0..1000 {
            let tof = 50_000.0 + i as f64 * 200.0;
            let m = base.mz(tof, t);
            let ppm = if i < 300 {
                8.0
            } else {
                -150.0 + 300.0 * ((i * 37 % 700) as f64 / 700.0)
            };
            points.push(CalibrationPoint {
                tof,
                mz: m / (1.0 + ppm * 1e-6),
                t,
                base_mz: f64::NAN,
            });
        }
        let off = mode_offset_ppm(&points, &base).unwrap();
        assert!((off - 8.0).abs() < 0.2, "{off}");
    }

    #[test]
    fn the_offset_curve_follows_an_offset_that_changes_with_mz() {
        let base = MzCalibration::linear(10.0, 31.0 / 400_000.0, 400_000.0);
        let t = [0.0; 5];
        let (lo, hi) = (base.mz(10_000.0, t), base.mz(390_000.0, t));
        let truth = |m: f64| 80.0 * (hi - m) / (hi - lo);
        let mut points = Vec::new();
        for i in 0..4000 {
            let tof = 10_000.0 + i as f64 * 95.0;
            let m = base.mz(tof, t);
            let ppm = if i % 4 == 3 {
                -150.0 + 300.0 * ((i * 37 % 1000) as f64 / 1000.0)
            } else {
                truth(m)
            };
            points.push(CalibrationPoint {
                tof,
                mz: m / (1.0 + ppm * 1e-6),
                t,
                base_mz: f64::NAN,
            });
        }
        let curve = offset_curve(&points, &base).unwrap();
        for m in [lo + 0.2 * (hi - lo), 0.5 * (lo + hi), lo + 0.8 * (hi - lo)] {
            assert!((curve.at(m) - truth(m)).abs() < 3.0, "{m}: {}", curve.at(m));
        }
        assert_eq!(OffsetCurve::constant(4.0).at(123.0), 4.0);
    }

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
    fn base_residuals_use_each_points_own_row() {
        let base = MzCalibration::linear(10.0, 31.0 / 400_000.0, 400_000.0);
        let t = [25.7, 27.0, 5.0, 1.0, 1.0];
        let own = |tof: f64| base.mz(tof, t) * (1.0 + 40e-6);
        let points: Vec<CalibrationPoint> = (0..400)
            .map(|i| {
                let tof = 20_000.0 + (i as f64) * 900.0;
                CalibrationPoint {
                    tof,
                    mz: own(tof),
                    t,
                    base_mz: own(tof),
                }
            })
            .collect();
        assert!(mode_offset_ppm(&points, &base).unwrap().abs() < 1e-6);
        let fit = fit_calibration(&points, &base).unwrap();
        assert!(fit.base_median_ppm.abs() < 1e-6, "{fit:?}");
        assert!(!fit_beats_base(&fit));
        let run_base: Vec<CalibrationPoint> = points
            .iter()
            .map(|p| CalibrationPoint {
                base_mz: f64::NAN,
                ..*p
            })
            .collect();
        assert!((mode_offset_ppm(&run_base, &base).unwrap() + 40.0).abs() < 0.5);
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

    fn fit_with(median: f64, sigma: f64, base_median: f64, base_sigma: f64) -> CalibrationFit {
        CalibrationFit {
            calibration: MzCalibration::linear(10.0, 1e-4, 1.0),
            n_used: 1000,
            median_ppm: median,
            sigma_ppm: sigma,
            base_median_ppm: base_median,
            base_sigma_ppm: base_sigma,
        }
    }

    #[test]
    fn a_fit_replaces_the_tdf_model_only_when_clearly_closer() {
        assert!(!fit_beats_base(&fit_with(-0.1, 2.0, 0.1, 2.1)));
        assert!(!fit_beats_base(&fit_with(0.0, 1.8, 0.3, 2.0)));
        assert!(fit_beats_base(&fit_with(-0.18, 2.07, 16.78, 7.63)));
    }
}
