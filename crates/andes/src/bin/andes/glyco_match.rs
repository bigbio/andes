//! Cross-run match scoring for `--glyco --mbr` (a 10 → 64 → 32 → 8 → 2 MLP).
//!
//! MBR does not transfer a glycopeptide by trusting an RT prediction alone: it
//! extracts the identifying run's isotope XIC (the "reference"), slides that
//! window across the target run's predicted RT, and scores every candidate
//! with a small MLP over ten similarity features. The best
//! candidate's own bounds become the quantification window, and
//! `0.85 · p + 0.15 · cosine` is the match score that separates true transfers
//! from decoys. The method is the published one of Zeng et al., *Nat. Commun.*
//! 2022; this module keeps its feature definitions, candidate enumeration and
//! two-score combination, and re-trains the weights on andes benchmark data
//! instead of the authors' private training set (the shipped weights carry
//! their provenance in `resources/quant/glyco_mbr_model.txt`). The differences
//! are adaptations to andes' data model, and each is deliberate:
//!
//! * **The reference is a fixed 1 s RT grid**, not a copy of the donor's scan
//!   indices. The original compares equal-length slices of two runs whose scan
//!   rates happen to agree; resampling both to seconds makes the features
//!   independent of the instrument's scan rate.
//! * **The candidate shift is measured from the RT-map prediction**, and the
//!   search is bounded by the alignment uncertainty (60–120 s) instead of ±4
//!   min from the raw seed RT. The RT map is a better prior than "same RT", and
//!   the shift feature becomes the deviation from that prediction.
//! * **Feature 3 repeats feature 2.** Both M+1 trace cosines are computed from
//!   the same reference trace; the duplicate is kept so the retrained model
//!   keeps the ten-feature input the method defines.

use std::collections::HashMap;

use model::mass::ISOTOPE;
use model::Tolerance;
use quant::ms1_index::Ms1RunIndex;
use quant::xic::cosine;

/// Features the model consumes, in the published order.
pub(crate) const N_FEATURES: usize = 10;

/// Candidates the model scores per (donor, target) pair (the top 10).
const TOP_CANDIDATES: usize = 10;

/// RT grid every matrix is resampled to, seconds.
const GRID_S: f64 = 1.0;

/// The identifying run's match reference: the isotope XIC of one feature, on a
/// uniform RT grid, plus the peak bounds the quantification found in it.
#[derive(Debug, Clone)]
pub(crate) struct Reference {
    /// Isotope m/z of row 0 (the observed monoisotope).
    mono_mz: f64,
    charge: f64,
    n_iso: usize,
    /// Isotope traces over the extraction window, row 0 = monoisotope.
    window: Vec<Vec<f32>>,
    /// Peak bounds as `[start, stop)` grid indices.
    peak: (usize, usize),
    /// Per-isotope summed intensity over the peak window.
    peak_sums: Vec<f64>,
    /// Width of the extraction window, seconds.
    width_s: f64,
    /// Donor anchor RT minus the window start, seconds.
    anchor_offset_s: f64,
}

/// The window chosen in the target run, with its match score.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Match {
    /// RT of the chosen peak window's start, seconds.
    pub(crate) start_rt: f64,
    /// RT of the chosen peak window's end, seconds.
    pub(crate) stop_rt: f64,
    /// Half-window for the re-extraction, seconds (`quantify` searches
    /// around the window centre).
    pub(crate) window_s: f64,
    /// The match score `0.85 · p + 0.15 · cosine`.
    pub(crate) score: f32,
    /// The ten model features of the chosen candidate.
    pub(crate) features: [f64; N_FEATURES],
}

impl Reference {
    /// Build the reference of a quantified feature: its isotope XIC over the
    /// extraction window (`anchor_rt ± rt_window_s`), resampled to [`GRID_S`],
    /// with `peak_start_rt..=peak_stop_rt` as the peak window.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn build(
        ms1: &Ms1RunIndex,
        mono_mz: f64,
        charge: u8,
        n_iso: usize,
        anchor_rt: f64,
        rt_window_s: f64,
        peak_start_rt: f64,
        peak_stop_rt: f64,
        tol: Tolerance,
    ) -> Option<Self> {
        let z = charge.max(1) as f64;
        let window = ms1.scan_window(anchor_rt - rt_window_s, anchor_rt + rt_window_s);
        if window.len() < 3 {
            return None;
        }
        let rt: Vec<f64> = window.clone().map(|i| ms1.rt(i)).collect();
        let width_s = rt[rt.len() - 1] - rt[0];
        if width_s < 2.0 * GRID_S {
            return None;
        }
        let n_grid = (width_s / GRID_S).round() as usize + 1;
        let grid: Vec<f64> = (0..n_grid).map(|i| rt[0] + i as f64 * GRID_S).collect();
        let traces: Vec<Vec<f32>> = (0..n_iso.max(1))
            .map(|k| {
                let mz = mono_mz + k as f64 * ISOTOPE / z;
                let raw = ms1.xic(mz, tol.as_da(mz), window.clone());
                resample(&rt, &raw, &grid)
            })
            .collect();
        let p0 = (((peak_start_rt - rt[0]) / GRID_S).round() as isize).clamp(0, n_grid as isize - 1)
            as usize;
        let mut p1 = (((peak_stop_rt - rt[0]) / GRID_S).round() as isize + 1)
            .clamp(0, n_grid as isize) as usize;
        if p1 <= p0 + 1 {
            p1 = (p0 + 2).min(n_grid);
        }
        let peak_sums: Vec<f64> = traces
            .iter()
            .map(|t| t[p0..p1].iter().map(|&v| v as f64).sum())
            .collect();
        Some(Self {
            mono_mz,
            charge: z,
            n_iso: traces.len(),
            window: traces,
            peak: (p0, p1),
            peak_sums,
            width_s,
            anchor_offset_s: anchor_rt - rt[0],
        })
    }

    /// The window start the donor's anchor maps to, given the target run's
    /// predicted anchor RT.
    pub(crate) fn window_start(&self, predicted_anchor_rt: f64) -> f64 {
        predicted_anchor_rt - self.anchor_offset_s
    }
}

/// Score every candidate shift of `reference` in `ms1` and return the best by
/// match score. `predicted_start_rt` is the target-run RT of the
/// reference window's start (the donor's `window_start` mapped through the RT
/// map); `search_s` bounds how far the window may slide.
pub(crate) fn best_match(
    ms1: &Ms1RunIndex,
    reference: &Reference,
    predicted_start_rt: f64,
    search_s: f64,
    tol: Tolerance,
) -> Option<Match> {
    let candidates = score_candidates(ms1, reference, predicted_start_rt, search_s, tol)?;
    scored_top(&candidates).into_iter().next()
}

/// Every candidate's features, in RT order. Only the top [`TOP_CANDIDATES`]
/// reach the model, but the training dump needs them all for labelling.
pub(crate) struct CandidateSet {
    /// Candidate window starts, RT order.
    pub(crate) start_rt: Vec<f64>,
    /// Candidate peak-window stop, RT order.
    pub(crate) stop_rt: Vec<f64>,
    /// The ten features of each candidate, RT order.
    pub(crate) features: Vec<[f64; N_FEATURES]>,
    /// The peak-window half-width used for re-extraction, seconds.
    pub(crate) window_s: f64,
    /// Peak-window length in grid points (the non-maximum suppression radius).
    peak_pts: usize,
}

/// Enumerate and score every candidate window shift.
pub(crate) fn score_candidates(
    ms1: &Ms1RunIndex,
    reference: &Reference,
    predicted_start_rt: f64,
    search_s: f64,
    tol: Tolerance,
) -> Option<CandidateSet> {
    let region_start = predicted_start_rt - search_s;
    let region_stop = predicted_start_rt + reference.width_s + search_s;
    let window = ms1.scan_window(region_start, region_stop);
    if window.len() < 3 {
        return None;
    }
    let rt: Vec<f64> = window.clone().map(|i| ms1.rt(i)).collect();
    let span_s = rt[rt.len() - 1] - rt[0];
    if span_s < reference.width_s + GRID_S {
        return None;
    }
    let n_grid = (span_s / GRID_S).round() as usize + 1;
    let grid: Vec<f64> = (0..n_grid).map(|i| rt[0] + i as f64 * GRID_S).collect();
    let traces: Vec<Vec<f32>> = (0..reference.n_iso)
        .map(|k| {
            let mz = reference.mono_mz + k as f64 * ISOTOPE / reference.charge;
            let raw = ms1.xic(mz, tol.as_da(mz), window.clone());
            resample(&rt, &raw, &grid)
        })
        .collect();

    let n_win = (reference.width_s / GRID_S).round() as usize + 1;
    if n_win + 2 > n_grid {
        return None;
    }
    let (p0, p1) = reference.peak;
    let peak_pts = p1 - p0;
    // Candidate window starts: one per grid point within ±search of the
    // prediction, while the whole reference window still fits the region.
    let first = 0usize;
    let last = (n_grid - n_win).min(((2.0 * search_s) / GRID_S).round() as usize);
    if last <= first {
        return None;
    }
    let ref_window_flat = flatten(&reference.window);
    let ref_peak_flat = flatten_peak(&reference.window, p0, p1);

    let mut set = CandidateSet {
        start_rt: Vec::with_capacity(last - first + 1),
        stop_rt: Vec::with_capacity(last - first + 1),
        features: Vec::with_capacity(last - first + 1),
        window_s: ((p1 - p0) as f64 * GRID_S / 2.0).max(5.0) + GRID_S,
        peak_pts,
    };
    for (k, &start) in grid.iter().enumerate().skip(first).take(last + 1 - first) {
        let sam_window: Vec<f64> = (0..reference.n_iso)
            .flat_map(|i| traces[i][k..k + n_win].iter().map(|&v| v as f64))
            .collect();
        let sam_peak: Vec<Vec<f64>> = (0..reference.n_iso)
            .map(|i| {
                traces[i][k + p0..k + p1]
                    .iter()
                    .map(|&v| v as f64)
                    .collect()
            })
            .collect();
        let sam_sums: Vec<f64> = sam_peak.iter().map(|t| t.iter().sum()).collect();
        let sam_peak_flat: Vec<f64> = sam_peak.iter().flatten().copied().collect();
        let mono = 0usize;
        let mut f = [0.0f64; N_FEATURES];
        // [0] peak-window envelope cosine over all isotopes.
        f[0] = cosine(&sam_peak_flat, &ref_peak_flat);
        // [1] monoisotopic trace, [2] and [3] the M+1 trace (COS_2 is
        // computed from the M+1 trace as well).
        f[1] = trace_cosine(&sam_peak, mono, &reference.window, p0, p1);
        f[2] = trace_cosine(&sam_peak, mono + 1, &reference.window, p0, p1);
        f[3] = f[2];
        // [4] and [5]: the candidate's own mono-vs-M+1/M+2 coherence.
        f[4] = if mono + 1 < sam_peak.len() {
            cosine(&sam_peak[mono], &sam_peak[mono + 1])
        } else {
            0.5
        };
        f[5] = if mono + 2 < sam_peak.len() {
            cosine(&sam_peak[mono], &sam_peak[mono + 2])
        } else {
            0.5
        };
        // [6] the isotope-intensity distribution.
        f[6] = cosine(&sam_sums, &reference.peak_sums);
        // [7] the shift from the RT-map prediction, in window widths.
        let shift = start - predicted_start_rt;
        f[7] = (shift.abs() / reference.width_s).min(1.0);
        // [8] the monoisotopic intensity ratio (log10, quartered).
        let ratio = reference.peak_sums.first().copied().unwrap_or(0.0)
            / sam_sums.first().copied().unwrap_or(0.0).max(1.0);
        f[8] = ratio.clamp(0.001, 1000.0).log10().abs() / 4.0;
        // [9] the extraction-window envelope cosine.
        f[9] = cosine(&sam_window, &ref_window_flat);
        set.start_rt.push(start + p0 as f64 * GRID_S);
        set.stop_rt.push(start + (p1 - 1) as f64 * GRID_S);
        set.features.push(f);
    }
    Some(set)
}

/// The candidates the model scores (top [`TOP_CANDIDATES`] by the screening
/// score, non-maximum suppressed around the peak width) with their final score,
/// best first.
pub(crate) fn scored_top(candidates: &CandidateSet) -> Vec<Match> {
    let radius = (candidates.peak_pts / 2).max(1) as isize;
    let mut order: Vec<usize> = (0..candidates.features.len()).collect();
    order.sort_by(|&a, &b| {
        let sa = screen(&candidates.features[a]);
        let sb = screen(&candidates.features[b]);
        sb.total_cmp(&sa).then(a.cmp(&b))
    });
    let mut picked: Vec<usize> = Vec::with_capacity(TOP_CANDIDATES);
    for i in order {
        if picked
            .iter()
            .all(|&j| (i as isize - j as isize).abs() > radius)
        {
            picked.push(i);
            if picked.len() == TOP_CANDIDATES {
                break;
            }
        }
    }
    if picked.is_empty() {
        return Vec::new();
    }
    let model = Model::shipped();
    let dims = model.dims;
    let mut logits = Vec::with_capacity(picked.len() * dims);
    for &i in &picked {
        let out = model.forward(&candidates.features[i]);
        logits.extend_from_slice(&out);
    }
    // Softmax over candidates of the class-1 logit (the method's
    // `torch.softmax(prediction, 0)`).
    let n = picked.len();
    let max_logit = (0..n)
        .map(|c| logits[c * dims + 1])
        .fold(f64::NEG_INFINITY, f64::max);
    let exp_sum: f64 = (0..n)
        .map(|c| (logits[c * dims + 1] - max_logit).exp())
        .sum();
    let mut out: Vec<Match> = picked
        .iter()
        .enumerate()
        .map(|(c, &i)| {
            let p = if exp_sum > 0.0 {
                (logits[c * dims + 1] - max_logit).exp() / exp_sum
            } else {
                0.0
            };
            let sim = candidates.features[i][0].clamp(0.0, 1.0);
            Match {
                start_rt: candidates.start_rt[i],
                stop_rt: candidates.stop_rt[i],
                window_s: candidates.window_s,
                score: (0.85 * p + 0.15 * sim).min(1.0) as f32,
                features: candidates.features[i],
            }
        })
        .collect();
    out.sort_by(|a, b| b.score.total_cmp(&a.score));
    out
}

/// The screening score: the small-window cosine plus the
/// extraction-window cosine.
fn screen(f: &[f64; N_FEATURES]) -> f64 {
    f[0] + f[9]
}

/// Cosine of one isotope trace against the reference's same trace.
fn trace_cosine(sam: &[Vec<f64>], k: usize, reference: &[Vec<f32>], p0: usize, p1: usize) -> f64 {
    match (sam.get(k), reference.get(k)) {
        (Some(s), Some(r)) => {
            let r: Vec<f64> = r[p0..p1].iter().map(|&v| v as f64).collect();
            cosine(s, &r)
        }
        _ => 0.5,
    }
}

fn flatten(traces: &[Vec<f32>]) -> Vec<f64> {
    traces
        .iter()
        .flat_map(|t| t.iter().map(|&v| v as f64))
        .collect()
}

fn flatten_peak(traces: &[Vec<f32>], p0: usize, p1: usize) -> Vec<f64> {
    traces
        .iter()
        .flat_map(|t| t[p0..p1].iter().map(|&v| v as f64))
        .collect()
}

/// Linear resampling of one trace onto a uniform grid; 0 outside the raw range.
fn resample(rt: &[f64], y: &[f32], grid: &[f64]) -> Vec<f32> {
    let mut out = Vec::with_capacity(grid.len());
    let mut j = 0usize;
    for &t in grid {
        while j + 1 < rt.len() && rt[j + 1] < t {
            j += 1;
        }
        if t < rt[0] || t > rt[rt.len() - 1] {
            out.push(0.0);
            continue;
        }
        if j + 1 >= rt.len() {
            out.push(y[y.len() - 1]);
            continue;
        }
        let dt = rt[j + 1] - rt[j];
        let v = if dt <= 0.0 {
            y[j]
        } else {
            let w = ((t - rt[j]) / dt) as f32;
            y[j] * (1.0 - w) + y[j + 1] * w
        };
        out.push(v);
    }
    out
}

/// The shipped match weights (10 → 64 → 32 → 8 → 2, BatchNorm + ReLU).
pub(crate) struct Model {
    dims: usize,
    fc1: Vec<f64>,
    b1: Vec<f64>,
    bn_w: Vec<f64>,
    bn_b: Vec<f64>,
    bn_mean: Vec<f64>,
    bn_var: Vec<f64>,
    fc2: Vec<f64>,
    b2: Vec<f64>,
    fc3: Vec<f64>,
    b3: Vec<f64>,
    fc4: Vec<f64>,
    b4: Vec<f64>,
}

impl Model {
    /// The model compiled into the binary. See
    /// `resources/quant/glyco_mbr_model.txt` for the training provenance.
    pub(crate) fn shipped() -> &'static Model {
        static MODEL: std::sync::OnceLock<Model> = std::sync::OnceLock::new();
        MODEL.get_or_init(|| {
            Model::parse(include_str!(
                "../../../../../resources/quant/glyco_mbr_model.txt"
            ))
            .expect("the shipped glyco MBR match model parses")
        })
    }

    fn parse(text: &str) -> Result<Model, String> {
        let mut fields: HashMap<&str, Vec<f64>> = HashMap::new();
        let mut arch: Vec<usize> = Vec::new();
        for line in text.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut it = line.split_whitespace();
            let key = it.next().ok_or("empty model line")?;
            let values: Result<Vec<f64>, _> = it.map(str::parse).collect();
            let values = values.map_err(|e| format!("model {key}: {e}"))?;
            if key == "arch" {
                arch = values.iter().map(|&v| v as usize).collect();
            } else {
                fields.insert(key, values);
            }
        }
        if arch.len() != 5 {
            return Err("model arch must be `input h1 h2 h3 out`".into());
        }
        let (n, h1, h2, h3, out) = (arch[0], arch[1], arch[2], arch[3], arch[4]);
        let mut get = |name: &str, len: usize| -> Result<Vec<f64>, String> {
            let v = fields
                .remove(name)
                .ok_or_else(|| format!("model field {name} is missing"))?;
            if v.len() != len {
                return Err(format!(
                    "model field {name}: {} values, want {len}",
                    v.len()
                ));
            }
            Ok(v)
        };
        let model = Model {
            dims: out,
            fc1: get("fc1.weight", h1 * n)?,
            b1: get("fc1.bias", h1)?,
            bn_w: get("bn1.weight", h1)?,
            bn_b: get("bn1.bias", h1)?,
            bn_mean: get("bn1.running_mean", h1)?,
            bn_var: get("bn1.running_var", h1)?,
            fc2: get("fc2.weight", h2 * h1)?,
            b2: get("fc2.bias", h2)?,
            fc3: get("fc3.weight", h3 * h2)?,
            b3: get("fc3.bias", h3)?,
            fc4: get("fc4.weight", out * h3)?,
            b4: get("fc4.bias", out)?,
        };
        if n != N_FEATURES {
            return Err(format!("model input is {n} features, want {N_FEATURES}"));
        }
        Ok(model)
    }

    /// My_NetNN.forward in eval mode: ReLU → BatchNorm(after the first layer),
    /// ReLU on every later layer including the output.
    fn forward(&self, x: &[f64; N_FEATURES]) -> Vec<f64> {
        let h1 = self.b1.len();
        let mut a = vec![0.0f64; h1];
        for (o, a_o) in a.iter_mut().enumerate() {
            let row = &self.fc1[o * N_FEATURES..(o + 1) * N_FEATURES];
            let mut s = self.b1[o];
            for (w, v) in row.iter().zip(x.iter()) {
                s += w * v;
            }
            *a_o = s.max(0.0);
        }
        for (o, a_o) in a.iter_mut().enumerate() {
            *a_o = (*a_o - self.bn_mean[o]) / (self.bn_var[o] + 1e-5).sqrt() * self.bn_w[o]
                + self.bn_b[o];
        }
        let h2 = self.b2.len();
        let mut b = vec![0.0f64; h2];
        for (o, b_o) in b.iter_mut().enumerate() {
            let row = &self.fc2[o * h1..(o + 1) * h1];
            let mut s = self.b2[o];
            for (w, v) in row.iter().zip(a.iter()) {
                s += w * v;
            }
            *b_o = s.max(0.0);
        }
        let h3 = self.b3.len();
        let mut c = vec![0.0f64; h3];
        for (o, c_o) in c.iter_mut().enumerate() {
            let row = &self.fc3[o * h2..(o + 1) * h2];
            let mut s = self.b3[o];
            for (w, v) in row.iter().zip(b.iter()) {
                s += w * v;
            }
            *c_o = s.max(0.0);
        }
        let out = self.b4.len();
        let mut d = vec![0.0f64; out];
        for (o, d_o) in d.iter_mut().enumerate() {
            let row = &self.fc4[o * h3..(o + 1) * h3];
            let mut s = self.b4[o];
            for (w, v) in row.iter().zip(c.iter()) {
                s += w * v;
            }
            *d_o = s.max(0.0);
        }
        d
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::scan::Ms1Scan;

    const MONO: f64 = 1000.5;
    const CHARGE: u8 = 2;
    const RATIOS: [f32; 4] = [1.0, 0.6, 0.3, 0.15];

    fn gauss(rt: f64, centre: f64, sigma: f64) -> f32 {
        (-0.5 * ((rt - centre) / sigma).powi(2)).exp() as f32
    }

    /// A synthetic run: one 4-isotope envelope eluting at `centre`, sampled at
    /// 0.5 s (finer than the 1 s grid, like a real QExactive survey).
    fn run(centre: f64) -> Ms1RunIndex {
        let scans = (0..601)
            .map(|i| {
                let rt = i as f64 * 0.5;
                let amplitude = gauss(rt, centre, 5.0);
                let peaks = RATIOS
                    .iter()
                    .enumerate()
                    .map(|(k, ratio)| {
                        (
                            MONO + k as f64 * ISOTOPE / CHARGE as f64,
                            1000.0 * ratio * amplitude,
                        )
                    })
                    .collect();
                Ms1Scan {
                    rt,
                    peaks,
                    ..Ms1Scan::default()
                }
            })
            .collect();
        Ms1RunIndex::new(scans)
    }

    fn reference(donor: &Ms1RunIndex, mono_mz: f64) -> Reference {
        Reference::build(
            donor,
            mono_mz,
            CHARGE,
            RATIOS.len(),
            100.0,
            60.0,
            90.0,
            110.0,
            Tolerance::Ppm(10.0),
        )
        .expect("the synthetic donor window has enough scans")
    }

    #[test]
    fn the_reference_window_start_round_trips_through_the_anchor() {
        let donor = run(100.0);
        let r = reference(&donor, MONO);
        assert!((r.window_start(100.0) - 40.0).abs() < 0.6);
    }

    #[test]
    fn the_match_chooses_the_window_where_the_peak_actually_is() {
        let donor = run(100.0);
        let target = run(150.0);
        let r = reference(&donor, MONO);
        let tol = Tolerance::Ppm(10.0);
        // The target apex is 50 s later; the RT map predicts exactly that.
        let predicted = r.window_start(150.0);
        let m = best_match(&target, &r, predicted, 60.0, tol).expect("a scored window");
        let centre = 0.5 * (m.start_rt + m.stop_rt);
        assert!(
            (centre - 150.0).abs() < 4.0,
            "chosen window centre {centre} is not at the elution apex"
        );
        assert_eq!(m.features.len(), N_FEATURES);
        assert!(m.score > 0.2, "match score {} is implausibly low", m.score);
    }

    #[test]
    fn the_match_recovers_a_peak_far_from_the_prediction() {
        let donor = run(100.0);
        let target = run(170.0);
        let r = reference(&donor, MONO);
        let tol = Tolerance::Ppm(10.0);
        let predicted = r.window_start(150.0); // 20 s off
        let m = best_match(&target, &r, predicted, 60.0, tol).expect("a scored window");
        let centre = 0.5 * (m.start_rt + m.stop_rt);
        assert!(
            (centre - 170.0).abs() < 4.0,
            "chosen window centre {centre}"
        );
    }

    #[test]
    fn the_decoy_reference_scores_below_the_target_reference() {
        let donor = run(100.0);
        let target = run(150.0);
        let tol = Tolerance::Ppm(10.0);
        let target_ref = reference(&donor, MONO);
        let decoy_ref = reference(&donor, MONO + 11.0);
        let m_target = best_match(
            &target,
            &target_ref,
            target_ref.window_start(150.0),
            60.0,
            tol,
        )
        .expect("a target window");
        let m_decoy = best_match(
            &target,
            &decoy_ref,
            decoy_ref.window_start(150.0),
            60.0,
            tol,
        )
        .expect("a decoy window");
        assert!(
            m_decoy.score < m_target.score,
            "decoy {} should lose to target {}",
            m_decoy.score,
            m_target.score
        );
    }

    #[test]
    fn the_model_forward_matches_its_shipped_dimensions() {
        let out = Model::shipped().forward(&[0.0; N_FEATURES]);
        assert_eq!(out.len(), 2);
        assert!(out.iter().all(|v| v.is_finite()));
    }
}
