//! Joint retention-time alignment across runs, without a reference run.
//!
//! Model: `rt[p, r] = f_r(t_p) + e`, where `t_p` is a consensus RT per
//! precursor, `f_r` is a monotone piecewise-linear map per run and the errors
//! `e` are heavy-tailed. All runs enter the fit symmetrically:
//!
//! 1. Anchors are precursors observed in at least two aligned runs of one
//!    alignment component (below), once per such component; only those
//!    observations enter the fit. `t_p` starts as the median of their
//!    observed RTs.
//! 2. Each round fits every run's map to its `(t_p, rt)` pairs (equal-count
//!    bins, the median x and y per bin, pool-adjacent-violators on y, linear
//!    interpolation, least-squares end slopes from the outer 20% of anchors),
//!    sets `t_p` to the median over runs of `f_r^{-1}(rt)`, and fixes the
//!    gauge by mapping `t_p` through the mean map `g = mean_r f_r` of the
//!    anchor's component, so each component's consensus scale is its average
//!    run.
//! 3. After the last round the maps are fitted once more to the gauge-fixed
//!    `t_p`, and every precursor's consensus is the median over its aligned
//!    runs of `f_r^{-1}(rt)` under those final maps (for a precursor seen in
//!    one run, just `f_r^{-1}(rt)`).
//! 4. Residuals are leave-one-out: for an observation of `p` in run `r`,
//!    `e = (rt - f_r(t_p^(-r))) / sqrt(1 + 1 / (n_p - 1))`, where `t_p^(-r)`
//!    is the consensus without run `r`; the scale makes `e` estimate the
//!    per-observation noise (an in-sample residual against a median over few
//!    runs collapses towards 0). Each run's residual width `sigma_r(t)` is
//!    `1.4826 * MAD` of the anchors' `e` in the same bins as the map, smoothed
//!    by a running median of 3 and floored at [`AlignParams::min_spacing_s`].
//!    Each precursor's extra dispersion is
//!    `tau_p^2 = max(0, mean_r e^2 - mean_r sigma_r^2) * n_p / (n_p + 3)`
//!    over the same residuals (0 when `n_p = 1`).
//!
//! The transfer uncertainty into run `q` combines the recipient's residual
//! width, the donors' residual widths and the extra dispersion:
//! `s^2 = sigma_q(t_p)^2 + f_q'(t_p)^2 * V_p + tau_p^2`, with
//! `V_p = sum_{r in D_p} (sigma_r(t_p) / f_r'(t_p))^2 / n_p^2` the variance
//! of the consensus RT on the consensus scale: each donor run's width is
//! carried into consensus time by its map's local slope `f_r'`, averaged over
//! the `n_p` donor runs `D_p`, and carried into run `q` by `f_q'`. A
//! precursor seen in one noisy run thus gets that run's width in every
//! recipient. With equal widths and unit slopes this is
//! `sigma^2 * (1 + 1 / n_p) + tau_p^2`. The local slope of a map at `t` is
//! the chord over the knot segment holding `t` and one segment on each side
//! (the end slope beyond the knots).
//!
//! Two aligned runs that share at least [`AlignParams::min_shared`]
//! precursors are joined; the connected components of that graph are the
//! alignment components, numbered in the order of their lowest run. Only
//! shared precursors fix the relative RT scale of two runs, so each
//! component is fitted and gauged on its own anchors (with one component,
//! this is the plain joint fit), and a transfer stays inside a component: the
//! consensus RT, `n_p`, `tau_p` and `V_p` of a prediction into run `q` use
//! only `p`'s observations in `q`'s component, and there is no prediction
//! when `p` has none there. [`Alignment::consensus`] is taken in the
//! component holding most of `p`'s observations (the lowest component id on
//! ties). A run alone in its component is unaligned.
//!
//! A run is unaligned when it has fewer than [`AlignParams::min_anchors`]
//! anchors, or when its anchors span less than [`AlignParams::min_coverage`]
//! of the anchors' consensus RT range (1st to 99th percentile of the initial
//! `t_p`). Removing a run can drop others below the thresholds, so the check
//! repeats until no run changes. A run that passes these checks but shares
//! fewer than [`AlignParams::min_shared`] precursors with every other aligned
//! run is then unaligned too. An unaligned run is excluded from the fit and
//! has no map; its observations are ignored everywhere, and [`Alignment`]
//! records its anchor count and coverage so the caller can log why.
//!
//! The output is a deterministic function of the set of observations: they are
//! processed in `(precursor, run)` order and every loop runs over sorted ids.

/// One precursor observed in one run.
///
/// `rt` is the retention time, in seconds, of the precursor's best-scoring PSM
/// in that run.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Observation {
    pub precursor: u32,
    pub run: u16,
    pub rt: f64,
}

/// Alignment settings.
#[derive(Debug, Clone, PartialEq)]
pub struct AlignParams {
    /// Number of alternating rounds (maps, then consensus, then gauge).
    pub rounds: usize,
    /// A run with fewer anchors than this is unaligned.
    pub min_anchors: usize,
    /// A run whose anchors span less than this fraction of the consensus RT
    /// range is unaligned.
    pub min_coverage: f64,
    /// Floor of the residual width `sigma_r(t)`, in seconds; the caller sets it
    /// to half the median MS1 scan spacing.
    pub min_spacing_s: f64,
    /// Two aligned runs that share at least this many precursors are in one
    /// alignment component; transfers never cross components.
    pub min_shared: usize,
}

impl Default for AlignParams {
    fn default() -> Self {
        Self {
            rounds: 8,
            min_anchors: 50,
            min_coverage: 0.30,
            min_spacing_s: 0.5,
            min_shared: 20,
        }
    }
}

/// Monotone piecewise-linear map from the consensus RT scale to one run.
///
/// Knot x and y are strictly increasing and both end slopes are positive, so
/// the map is strictly increasing and [`RunMap::invert`] is its exact inverse.
#[derive(Debug, Clone, PartialEq)]
pub struct RunMap {
    knots_x: Vec<f64>,
    knots_y: Vec<f64>,
    slope_lo: f64,
    slope_hi: f64,
}

impl RunMap {
    /// Run RT of consensus RT `t`.
    pub fn apply(&self, t: f64) -> f64 {
        interp(
            &self.knots_x,
            &self.knots_y,
            self.slope_lo,
            self.slope_hi,
            t,
        )
    }

    /// Consensus RT of run RT `rt` (binary search on the knots).
    pub fn invert(&self, rt: f64) -> f64 {
        interp(
            &self.knots_y,
            &self.knots_x,
            1.0 / self.slope_lo,
            1.0 / self.slope_hi,
            rt,
        )
    }

    /// Local slope `d run RT / d consensus RT` at consensus RT `t`: the chord
    /// over the knot segment holding `t` and one segment on each side, clamped
    /// to the knots; the end slope below the first or above the last knot.
    pub fn slope_at(&self, t: f64) -> f64 {
        let n = self.knots_x.len();
        if t.is_nan() {
            return f64::NAN;
        }
        if n < 2 || t <= self.knots_x[0] {
            return self.slope_lo;
        }
        if t >= self.knots_x[n - 1] {
            return self.slope_hi;
        }
        // knots_x[i - 1] <= t < knots_x[i], with 1 <= i <= n - 1.
        let i = self.knots_x.partition_point(|&k| k <= t);
        let (a, b) = (i.saturating_sub(2), (i + 1).min(n - 1));
        (self.knots_y[b] - self.knots_y[a]) / (self.knots_x[b] - self.knots_x[a])
    }

    /// Knot x (consensus RT) and y (run RT).
    pub fn knots(&self) -> (&[f64], &[f64]) {
        (&self.knots_x, &self.knots_y)
    }

    /// Slopes below the first and above the last knot.
    pub fn end_slopes(&self) -> (f64, f64) {
        (self.slope_lo, self.slope_hi)
    }

    /// A map from stored knots and end slopes, as returned by
    /// [`RunMap::knots`] and [`RunMap::end_slopes`].
    ///
    /// Errors: no knots, knot vectors of different lengths, a non-finite knot,
    /// or an end slope that is not finite and positive.
    pub fn from_parts(
        knots_x: Vec<f64>,
        knots_y: Vec<f64>,
        slope_lo: f64,
        slope_hi: f64,
    ) -> Result<RunMap, String> {
        if knots_x.is_empty() || knots_x.len() != knots_y.len() {
            return Err(format!(
                "run map needs equal, non-zero numbers of knots (x {}, y {})",
                knots_x.len(),
                knots_y.len()
            ));
        }
        if knots_x.iter().chain(&knots_y).any(|v| !v.is_finite()) {
            return Err("run map has a non-finite knot".into());
        }
        if !(slope_lo.is_finite() && slope_lo > 0.0 && slope_hi.is_finite() && slope_hi > 0.0) {
            return Err(format!(
                "run map end slopes {slope_lo} and {slope_hi} must be finite and positive"
            ));
        }
        Ok(RunMap {
            knots_x,
            knots_y,
            slope_lo,
            slope_hi,
        })
    }
}

/// The stored state of one run in an [`Alignment`]: everything
/// [`Alignment::from_runs`] needs to rebuild it exactly.
#[derive(Debug, Clone, PartialEq)]
pub struct RunAlignment {
    /// `None` for an unaligned run.
    pub map: Option<RunMap>,
    /// Residual width knots `(consensus t, sd)`; `None` for an unaligned run.
    pub residual_sd: Option<Vec<(f64, f64)>>,
    pub anchors: usize,
    pub coverage: f64,
    /// Alignment component; `None` for an unaligned run.
    pub component: Option<u32>,
    /// The run's observations `(precursor, rt)`, ascending by precursor;
    /// empty for an unaligned run.
    pub observations: Vec<(u32, f64)>,
}

/// Linear interpolation through `(xs, ys)` (xs strictly increasing, not
/// empty), extrapolated with `slope_lo` / `slope_hi` beyond the ends.
fn interp(xs: &[f64], ys: &[f64], slope_lo: f64, slope_hi: f64, t: f64) -> f64 {
    if t.is_nan() {
        return f64::NAN;
    }
    let n = xs.len();
    if t <= xs[0] {
        return ys[0] + slope_lo * (t - xs[0]);
    }
    if t >= xs[n - 1] {
        return ys[n - 1] + slope_hi * (t - xs[n - 1]);
    }
    // xs[i - 1] <= t < xs[i], with 1 <= i <= n - 1.
    let i = xs.partition_point(|&k| k <= t);
    let (x0, x1, y0, y1) = (xs[i - 1], xs[i], ys[i - 1], ys[i]);
    y0 + (y1 - y0) * (t - x0) / (x1 - x0)
}

/// Result of [`align`]. Vectors are indexed by precursor id or run id.
#[derive(Debug, Clone)]
pub struct Alignment {
    /// Consensus RT per precursor, in the component holding most of its
    /// observations; NaN when it has no observation in an aligned run.
    pub consensus: Vec<f64>,
    /// Map per run; `None` for an unaligned run.
    pub maps: Vec<Option<RunMap>>,
    /// Residual width per run as `(consensus t, sd)` knots, interpolated
    /// linearly and held constant beyond the ends; `None` for an unaligned run.
    pub residual_sd: Vec<Option<Vec<(f64, f64)>>>,
    /// Extra RT dispersion per precursor, in the same component as
    /// `consensus` (0 without observations).
    pub tau: Vec<f64>,
    /// Number of aligned runs in which each precursor is observed (`n_p`).
    pub n_runs_seen: Vec<u16>,
    /// Anchors per run (observations of precursors seen in at least two
    /// aligned runs), as counted when the run was last checked.
    pub anchors: Vec<usize>,
    /// Fraction of the consensus RT range covered by each run's anchors, as
    /// measured when the run was last checked.
    pub coverage: Vec<f64>,
    /// Alignment component per run; `None` for an unaligned run.
    components: Vec<Option<u32>>,
    /// Per precursor, the component of `consensus` and `tau` (`u32::MAX`
    /// without observations) and the consensus variance `V_p` there.
    home: Vec<u32>,
    consensus_var: Vec<f64>,
    /// Observations in aligned runs grouped by precursor: those of `p` are
    /// `obs_runs[obs_offsets[p]..obs_offsets[p + 1]]` (runs ascending).
    obs_offsets: Vec<usize>,
    obs_runs: Vec<u16>,
    obs_rts: Vec<f64>,
}

impl Alignment {
    /// Predicted apex RT of precursor `p` in `run` and its uncertainty `s`,
    /// from `p`'s observations in `run`'s alignment component.
    ///
    /// `None` when the run is unaligned or `p` has no observation in an
    /// aligned run of `run`'s component.
    pub fn predict(&self, p: u32, run: u16) -> Option<(f64, f64)> {
        let map = self.maps.get(run as usize)?.as_ref()?;
        let c = self.component(run)?;
        let pi = p as usize;
        if pi >= self.consensus.len() {
            return None;
        }
        if self.home[pi] == c {
            let t = self.consensus[pi];
            return Some(self.transfer(map, run, t, self.consensus_var[pi], self.tau[pi]));
        }
        let e = self.estimate(pi, None, c)?;
        Some(self.transfer(map, run, e.t, e.var, e.tau))
    }

    /// Like [`Alignment::predict`], but the consensus RT, `n_p`, `tau_p` and
    /// `V_p` are recomputed without `run`'s own observation of `p`. `None`
    /// when the run is unaligned or `p` has no observation in another aligned
    /// run of `run`'s component.
    pub fn predict_excluding(&self, p: u32, run: u16) -> Option<(f64, f64)> {
        let map = self.maps.get(run as usize)?.as_ref()?;
        let c = self.component(run)?;
        if p as usize >= self.consensus.len() {
            return None;
        }
        let e = self.estimate(p as usize, Some(run), c)?;
        Some(self.transfer(map, run, e.t, e.var, e.tau))
    }

    /// Alignment component of `run`; `None` for an unaligned run or a run
    /// outside the alignment.
    pub fn component(&self, run: u16) -> Option<u32> {
        self.components.get(run as usize).copied().flatten()
    }

    /// Number of alignment components.
    pub fn n_components(&self) -> usize {
        self.components
            .iter()
            .flatten()
            .max()
            .map_or(0, |&c| c as usize + 1)
    }

    /// True when `run` is aligned and `p` is observed in an aligned run
    /// other than `exclude`, but in none of `run`'s component: the
    /// prediction is refused only because it would cross components.
    pub fn blocked_by_component(&self, p: u32, run: u16, exclude: Option<u16>) -> bool {
        let Some(c) = self.component(run) else {
            return false;
        };
        let pi = p as usize;
        if pi >= self.consensus.len() {
            return false;
        }
        let mut other = false;
        for i in self.obs_offsets[pi]..self.obs_offsets[pi + 1] {
            let r = self.obs_runs[i];
            if Some(r) == exclude {
                continue;
            }
            if self.components[r as usize] == Some(c) {
                return false;
            }
            other = true;
        }
        other
    }

    /// The per-run state, in run order, from which [`Alignment::from_runs`]
    /// rebuilds this alignment.
    pub fn runs(&self) -> Vec<RunAlignment> {
        let n_runs = self.maps.len();
        let mut observations: Vec<Vec<(u32, f64)>> = vec![Vec::new(); n_runs];
        for p in 0..self.consensus.len() {
            for i in self.obs_offsets[p]..self.obs_offsets[p + 1] {
                observations[self.obs_runs[i] as usize].push((p as u32, self.obs_rts[i]));
            }
        }
        (0..n_runs)
            .zip(observations)
            .map(|(r, observations)| RunAlignment {
                map: self.maps[r].clone(),
                residual_sd: self.residual_sd[r].clone(),
                anchors: self.anchors[r],
                coverage: self.coverage[r],
                component: self.components[r],
                observations,
            })
            .collect()
    }

    /// Adds observations that take no part in the fit: the run maps, residual
    /// widths, anchors and components stay as they are, and each observed
    /// precursor gets its consensus RT, `n_p`, `tau_p` and `V_p` from its
    /// observations under the fitted maps, by the code every precursor's
    /// estimate uses. Observations in unaligned runs and with a non-finite RT
    /// are dropped.
    ///
    /// Errors: a precursor outside the alignment, or a second observation of
    /// one precursor in one run (counting the existing observations).
    pub fn add_observations(&mut self, obs: &[Observation]) -> Result<(), String> {
        if obs.is_empty() {
            return Ok(());
        }
        let n = self.consensus.len();
        let mut runs = self.runs();
        for o in obs {
            let Some(run) = runs.get_mut(o.run as usize) else {
                continue;
            };
            if run.map.is_some() && o.rt.is_finite() {
                run.observations.push((o.precursor, o.rt));
            }
        }
        *self = Alignment::from_runs(runs, n)?;
        Ok(())
    }

    /// Rebuilds an alignment of `n_precursors` precursors from the per-run
    /// state of [`Alignment::runs`]. The consensus RT, `tau_p` and `n_p` are
    /// recomputed from the observations by the same code [`align`] uses, so
    /// every prediction equals the original's bit for bit.
    ///
    /// Errors: a map without residual widths or the reverse, an aligned run
    /// without a component or the reverse, an observation in an unaligned
    /// run, a precursor outside `0..n_precursors`, a non-finite RT, or two
    /// observations of one precursor in one run.
    pub fn from_runs(runs: Vec<RunAlignment>, n_precursors: usize) -> Result<Alignment, String> {
        let n_runs = runs.len();
        if n_runs > u16::MAX as usize + 1 {
            return Err(format!("{n_runs} runs exceed the run id range"));
        }
        let mut all: Vec<(u32, u16, f64)> = Vec::new();
        for (r, run) in runs.iter().enumerate() {
            if run.map.is_some() != run.residual_sd.is_some() {
                return Err(format!(
                    "run {r}: a run map needs residual widths and the reverse"
                ));
            }
            if run.map.is_some() != run.component.is_some() {
                return Err(format!(
                    "run {r}: an aligned run needs a component and the reverse"
                ));
            }
            if run.map.is_none() && !run.observations.is_empty() {
                return Err(format!("run {r} is unaligned but has observations"));
            }
            for &(p, rt) in &run.observations {
                if p as usize >= n_precursors {
                    return Err(format!(
                        "run {r}: observation of precursor {p} outside the {n_precursors} precursors"
                    ));
                }
                if !rt.is_finite() {
                    return Err(format!("run {r}: precursor {p} has a non-finite RT"));
                }
                all.push((p, r as u16, rt));
            }
        }
        all.sort_by_key(|o| (o.0, o.1));
        if let Some(w) = all
            .windows(2)
            .find(|w| w[0].0 == w[1].0 && w[0].1 == w[1].1)
        {
            return Err(format!(
                "run {}: precursor {} is observed twice",
                w[0].1, w[0].0
            ));
        }
        let mut obs_offsets = vec![0usize; n_precursors + 1];
        for &(p, ..) in &all {
            obs_offsets[p as usize + 1] += 1;
        }
        for i in 0..n_precursors {
            obs_offsets[i + 1] += obs_offsets[i];
        }
        let n_runs_seen: Vec<u16> = (0..n_precursors)
            .map(|p| (obs_offsets[p + 1] - obs_offsets[p]).min(u16::MAX as usize) as u16)
            .collect();
        let mut maps = Vec::with_capacity(n_runs);
        let mut residual_sd = Vec::with_capacity(n_runs);
        let mut anchors = Vec::with_capacity(n_runs);
        let mut coverage = Vec::with_capacity(n_runs);
        let mut components = Vec::with_capacity(n_runs);
        for run in runs {
            maps.push(run.map);
            residual_sd.push(run.residual_sd);
            anchors.push(run.anchors);
            coverage.push(run.coverage);
            components.push(run.component);
        }
        let mut out = Alignment {
            consensus: vec![f64::NAN; n_precursors],
            maps,
            residual_sd,
            tau: vec![0.0; n_precursors],
            n_runs_seen,
            anchors,
            coverage,
            components,
            home: vec![u32::MAX; n_precursors],
            consensus_var: vec![0.0; n_precursors],
            obs_offsets,
            obs_runs: all.iter().map(|o| o.1).collect(),
            obs_rts: all.iter().map(|o| o.2).collect(),
        };
        out.fill_consensus();
        Ok(out)
    }

    /// Sets `consensus`, `tau`, `home` and `consensus_var` of every precursor
    /// from its observations in its home component.
    fn fill_consensus(&mut self) {
        let mut counts: Vec<(u32, usize)> = Vec::new();
        for p in 0..self.consensus.len() {
            counts.clear();
            for i in self.obs_offsets[p]..self.obs_offsets[p + 1] {
                let Some(c) = self.components[self.obs_runs[i] as usize] else {
                    continue;
                };
                match counts.iter_mut().find(|e| e.0 == c) {
                    Some(e) => e.1 += 1,
                    None => counts.push((c, 1)),
                }
            }
            // Most observations, then the lowest component id.
            let Some(&(c, _)) = counts
                .iter()
                .max_by(|a, b| a.1.cmp(&b.1).then(b.0.cmp(&a.0)))
            else {
                continue;
            };
            if let Some(e) = self.estimate(p, None, c) {
                self.consensus[p] = e.t;
                self.tau[p] = e.tau;
                self.home[p] = c;
                self.consensus_var[p] = e.var;
            }
        }
    }

    /// `(f_q(t), s)` for recipient `run` with map `map`, consensus `t`,
    /// consensus variance `var` and extra dispersion `tau`.
    fn transfer(&self, map: &RunMap, run: u16, t: f64, var: f64, tau: f64) -> (f64, f64) {
        let sigma = self.sigma_at(run, t);
        let slope = map.slope_at(t);
        let s = (sigma * sigma + slope * slope * var + tau * tau).sqrt();
        (map.apply(t), s)
    }

    /// Residual width of `run` at consensus RT `t` (0 for an unaligned run).
    fn sigma_at(&self, run: u16, t: f64) -> f64 {
        let Some(Some(knots)) = self.residual_sd.get(run as usize) else {
            return 0.0;
        };
        let n = knots.len();
        if t.is_nan() || n == 0 {
            return f64::NAN;
        }
        if t <= knots[0].0 {
            return knots[0].1;
        }
        if t >= knots[n - 1].0 {
            return knots[n - 1].1;
        }
        let i = knots.partition_point(|k| k.0 <= t);
        let ((x0, s0), (x1, s1)) = (knots[i - 1], knots[i]);
        s0 + (s1 - s0) * (t - x0) / (x1 - x0)
    }

    /// Consensus RT, `n_p`, `tau_p` and `V_p` of precursor `p` from its
    /// observations in the aligned runs of `component` other than `exclude`,
    /// under the final maps. `None` when none is left.
    fn estimate(&self, p: usize, exclude: Option<u16>, component: u32) -> Option<Estimate> {
        let range = self.obs_offsets[p]..self.obs_offsets[p + 1];
        let used = || {
            range
                .clone()
                .filter(move |&i| {
                    let r = self.obs_runs[i];
                    Some(r) != exclude && self.components[r as usize] == Some(component)
                })
                .map(|i| (self.obs_runs[i], self.obs_rts[i]))
        };
        let inv: Vec<f64> = used()
            .map(|(r, rt)| {
                self.maps[r as usize]
                    .as_ref()
                    .map_or(f64::NAN, |m| m.invert(rt))
            })
            .collect();
        if inv.is_empty() {
            return None;
        }
        let n = inv.len();
        let mut sorted = inv.clone();
        let t = median_in_place(&mut sorted);
        let nf = n as f64;
        let var = used()
            .map(|(r, _)| {
                let sd = self.sigma_at(r, t)
                    / self.maps[r as usize]
                        .as_ref()
                        .map_or(f64::NAN, |m| m.slope_at(t));
                sd * sd
            })
            .sum::<f64>()
            / (nf * nf);
        if n < 2 {
            return Some(Estimate { t, tau: 0.0, var });
        }
        let pairs: Vec<(u16, f64)> = used().collect();
        let res = loo_residuals(&self.maps, &pairs, &inv);
        let (mut e2, mut s2) = (0.0, 0.0);
        for (&(r, _), &e) in pairs.iter().zip(&res) {
            let s = self.sigma_at(r, t);
            e2 += e * e;
            s2 += s * s;
        }
        let tau2 = (e2 / nf - s2 / nf).max(0.0) * nf / (nf + 3.0);
        Some(Estimate {
            t,
            tau: tau2.sqrt(),
            var,
        })
    }
}

/// Consensus RT `t`, extra dispersion `tau` and consensus variance `var`
/// (`V_p`) of one precursor in one component.
struct Estimate {
    t: f64,
    tau: f64,
    var: f64,
}

/// Jointly aligns all runs; see the module documentation for the method.
///
/// Preconditions: `obs` holds observations of confident (library) precursors,
/// at most one per `(precursor, run)` — the RT of the best-scoring PSM — with
/// `precursor < n_precursors` and `run < n_runs`. In debug builds a violation
/// panics; in release builds a repeated `(precursor, run)` keeps the first
/// observation given, and out-of-range ids and non-finite RTs are ignored.
pub fn align(
    obs: &[Observation],
    n_precursors: usize,
    n_runs: usize,
    p: &AlignParams,
) -> Alignment {
    // Valid observations in (precursor, run) order; the stable sort keeps the
    // first one given for a repeated pair.
    let mut order: Vec<usize> = (0..obs.len())
        .filter(|&i| {
            let o = &obs[i];
            debug_assert!(
                (o.precursor as usize) < n_precursors && (o.run as usize) < n_runs,
                "observation out of range: precursor {} run {}",
                o.precursor,
                o.run
            );
            (o.precursor as usize) < n_precursors && (o.run as usize) < n_runs && o.rt.is_finite()
        })
        .collect();
    order.sort_by_key(|&i| (obs[i].precursor, obs[i].run));
    let mut all: Vec<(u32, u16, f64)> = Vec::with_capacity(order.len());
    for i in order {
        let o = &obs[i];
        let repeated = all
            .last()
            .is_some_and(|l| l.0 == o.precursor && l.1 == o.run);
        debug_assert!(
            !repeated,
            "precursor {} has several observations in run {}",
            o.precursor, o.run
        );
        if !repeated {
            all.push((o.precursor, o.run, o.rt));
        }
    }

    let (mut aligned, anchors, coverage) = select_runs(&all, n_precursors, n_runs, p);
    let mut components = find_components(&all, &aligned, p.min_shared);
    // A run alone in its component shares no anchor with a run of its own
    // component, so nothing fixes its map: it is unaligned.
    let mut size = vec![0usize; n_runs];
    for &c in components.iter().flatten() {
        size[c as usize] += 1;
    }
    let isolated: Vec<usize> = (0..n_runs)
        .filter(|&r| components[r].is_some_and(|c| size[c as usize] == 1))
        .collect();
    if !isolated.is_empty() {
        for &r in &isolated {
            aligned[r] = false;
        }
        components = find_components(&all, &aligned, p.min_shared);
    }
    let n_components = components
        .iter()
        .flatten()
        .max()
        .map_or(0, |&c| c as usize + 1);

    // Observations in aligned runs, grouped by precursor.
    let mut obs_offsets = vec![0usize; n_precursors + 1];
    let mut obs_runs = Vec::new();
    let mut obs_rts = Vec::new();
    for &(pr, r, rt) in &all {
        if aligned[r as usize] {
            obs_offsets[pr as usize + 1] += 1;
            obs_runs.push(r);
            obs_rts.push(rt);
        }
    }
    for i in 0..n_precursors {
        obs_offsets[i + 1] += obs_offsets[i];
    }
    let group = |pr: usize| obs_offsets[pr]..obs_offsets[pr + 1];

    // Anchors: a precursor seen in at least two runs of one component, once
    // per such component (ascending precursor, then component). Only those
    // observations enter the fit; anchor `a` has id `a` in `t`.
    let mut anchor_ids: Vec<(usize, u32)> = Vec::new();
    let mut counts: Vec<(u32, usize)> = Vec::new();
    for pr in 0..n_precursors {
        counts.clear();
        for i in group(pr) {
            let c = components[obs_runs[i] as usize].expect("aligned run has a component");
            match counts.iter_mut().find(|e| e.0 == c) {
                Some(e) => e.1 += 1,
                None => counts.push((c, 1)),
            }
        }
        counts.sort_unstable();
        anchor_ids.extend(counts.iter().filter(|e| e.1 >= 2).map(|e| (pr, e.0)));
    }
    let (offsets, runs_of, comp_of) = (&obs_offsets, &obs_runs, &components);
    let anchor_obs = move |&(pr, c): &(usize, u32)| {
        (offsets[pr]..offsets[pr + 1]).filter(move |&i| comp_of[runs_of[i] as usize] == Some(c))
    };
    let mut t = vec![f64::NAN; anchor_ids.len()];
    let mut buf = Vec::new();
    for (a, anchor) in anchor_ids.iter().enumerate() {
        buf.clear();
        buf.extend(anchor_obs(anchor).map(|i| obs_rts[i]));
        t[a] = median_in_place(&mut buf);
    }
    // Per run, the anchors' (anchor id, rt) pairs.
    let mut run_pairs: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n_runs];
    for (a, anchor) in anchor_ids.iter().enumerate() {
        for i in anchor_obs(anchor) {
            run_pairs[obs_runs[i] as usize].push((a, obs_rts[i]));
        }
    }

    let mut maps: Vec<Option<RunMap>> = vec![None; n_runs];
    let fit_all = |maps: &mut Vec<Option<RunMap>>, t: &[f64]| {
        for r in 0..n_runs {
            maps[r] = if aligned[r] {
                Some(fit_map(&sorted_pairs(&run_pairs[r], t)))
            } else {
                None
            };
        }
    };
    for _ in 0..p.rounds {
        fit_all(&mut maps, &t);
        for (a, anchor) in anchor_ids.iter().enumerate() {
            buf.clear();
            buf.extend(anchor_obs(anchor).map(|i| {
                maps[obs_runs[i] as usize]
                    .as_ref()
                    .map_or(f64::NAN, |m| m.invert(obs_rts[i]))
            }));
            t[a] = median_in_place(&mut buf);
        }
        // The gauge of each component is the mean map of its runs.
        let gauges: Vec<Option<RunMap>> = (0..n_components as u32)
            .map(|c| {
                mean_map(
                    (0..n_runs)
                        .filter(|&r| components[r] == Some(c))
                        .filter_map(|r| maps[r].as_ref()),
                )
            })
            .collect();
        for (a, &(_, c)) in anchor_ids.iter().enumerate() {
            if let Some(g) = &gauges[c as usize] {
                t[a] = g.apply(t[a]);
            }
        }
    }
    fit_all(&mut maps, &t);

    // Leave-one-out residuals of the anchors, as (consensus t_p, e) per run.
    let mut run_res: Vec<Vec<(f64, f64)>> = vec![Vec::new(); n_runs];
    let mut pairs: Vec<(u16, f64)> = Vec::new();
    for anchor in &anchor_ids {
        pairs.clear();
        pairs.extend(anchor_obs(anchor).map(|i| (obs_runs[i], obs_rts[i])));
        let inv: Vec<f64> = pairs
            .iter()
            .map(|&(r, rt)| maps[r as usize].as_ref().map_or(f64::NAN, |m| m.invert(rt)))
            .collect();
        buf.clear();
        buf.extend_from_slice(&inv);
        let tc = median_in_place(&mut buf);
        for (&(r, _), e) in pairs.iter().zip(loo_residuals(&maps, &pairs, &inv)) {
            run_res[r as usize].push((tc, e));
        }
    }
    let residual_sd: Vec<Option<Vec<(f64, f64)>>> = (0..n_runs)
        .map(|r| {
            maps[r].as_ref().map(|_| {
                let v = &mut run_res[r];
                v.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
                residual_width(v, p.min_spacing_s)
            })
        })
        .collect();

    let n_runs_seen: Vec<u16> = (0..n_precursors)
        .map(|pr| group(pr).len().min(u16::MAX as usize) as u16)
        .collect();
    let mut out = Alignment {
        consensus: vec![f64::NAN; n_precursors],
        maps,
        residual_sd,
        tau: vec![0.0; n_precursors],
        n_runs_seen,
        anchors,
        coverage,
        components,
        home: vec![u32::MAX; n_precursors],
        consensus_var: vec![0.0; n_precursors],
        obs_offsets,
        obs_runs,
        obs_rts,
    };
    out.fill_consensus();
    out
}

/// Alignment component per run: the connected components of the graph whose
/// edges join two aligned runs sharing at least `min_shared` precursors,
/// numbered from 0 in the order of their lowest run; `None` for an unaligned
/// run. `all` is sorted by precursor.
fn find_components(
    all: &[(u32, u16, f64)],
    aligned: &[bool],
    min_shared: usize,
) -> Vec<Option<u32>> {
    let n_runs = aligned.len();
    let mut per_run: Vec<Vec<u32>> = vec![Vec::new(); n_runs];
    for &(pr, r, _) in all {
        if aligned[r as usize] {
            per_run[r as usize].push(pr);
        }
    }
    let need = min_shared.max(1);
    let shares = |a: &[u32], b: &[u32]| {
        let (mut i, mut j, mut k) = (0, 0, 0);
        while i < a.len() && j < b.len() {
            match a[i].cmp(&b[j]) {
                std::cmp::Ordering::Less => i += 1,
                std::cmp::Ordering::Greater => j += 1,
                std::cmp::Ordering::Equal => {
                    k += 1;
                    if k >= need {
                        return true;
                    }
                    i += 1;
                    j += 1;
                }
            }
        }
        false
    };
    let mut parent: Vec<usize> = (0..n_runs).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    for r in 0..n_runs {
        if !aligned[r] {
            continue;
        }
        for q in 0..r {
            if !aligned[q] {
                continue;
            }
            let (a, b) = (find(&mut parent, r), find(&mut parent, q));
            if a != b && shares(&per_run[r], &per_run[q]) {
                parent[a.max(b)] = a.min(b);
            }
        }
    }
    let mut label: Vec<Option<u32>> = vec![None; n_runs];
    let mut out = vec![None; n_runs];
    let mut next = 0u32;
    for r in 0..n_runs {
        if !aligned[r] {
            continue;
        }
        let root = find(&mut parent, r);
        let c = *label[root].get_or_insert_with(|| {
            next += 1;
            next - 1
        });
        out[r] = Some(c);
    }
    out
}

/// Decides which runs are aligned: returns the flags, and each run's anchor
/// count and coverage when it was last checked.
fn select_runs(
    all: &[(u32, u16, f64)],
    n_precursors: usize,
    n_runs: usize,
    p: &AlignParams,
) -> (Vec<bool>, Vec<usize>, Vec<f64>) {
    let mut aligned = vec![true; n_runs];
    let mut anchors = vec![0usize; n_runs];
    let mut coverage = vec![0.0f64; n_runs];
    let mut t0 = vec![f64::NAN; n_precursors];
    let mut buf = Vec::new();
    loop {
        // `all` is sorted by precursor: walk each precursor's group.
        let mut start = 0;
        while start < all.len() {
            let pr = all[start].0;
            let mut end = start;
            while end < all.len() && all[end].0 == pr {
                end += 1;
            }
            buf.clear();
            buf.extend(
                all[start..end]
                    .iter()
                    .filter(|o| aligned[o.1 as usize])
                    .map(|o| o.2),
            );
            t0[pr as usize] = if buf.len() >= 2 {
                median_in_place(&mut buf)
            } else {
                f64::NAN
            };
            start = end;
        }
        let mut run_t: Vec<Vec<f64>> = vec![Vec::new(); n_runs];
        for &(pr, r, _) in all {
            let tp = t0[pr as usize];
            if aligned[r as usize] && tp.is_finite() {
                run_t[r as usize].push(tp);
            }
        }
        let mut every: Vec<f64> = t0.iter().copied().filter(|x| x.is_finite()).collect();
        let total = robust_span(&mut every);
        let mut changed = false;
        for r in 0..n_runs {
            if !aligned[r] {
                continue;
            }
            anchors[r] = run_t[r].len();
            let span = robust_span(&mut run_t[r]);
            coverage[r] = if total > 0.0 { span / total } else { 0.0 };
            if anchors[r] < p.min_anchors.max(1) || coverage[r] < p.min_coverage {
                aligned[r] = false;
                changed = true;
            }
        }
        if !changed {
            return (aligned, anchors, coverage);
        }
    }
}

/// Spread between the 1st and 99th percentiles (nearest rank); 0 when empty.
fn robust_span(v: &mut [f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    v.sort_by(f64::total_cmp);
    let q = |f: f64| v[((v.len() - 1) as f64 * f).round() as usize];
    q(0.99) - q(0.01)
}

/// A run's anchors as `(t_p, rt)`, sorted by `t_p` then `rt`.
fn sorted_pairs(pairs: &[(usize, f64)], t: &[f64]) -> Vec<(f64, f64)> {
    let mut v: Vec<(f64, f64)> = pairs.iter().map(|&(pr, rt)| (t[pr], rt)).collect();
    v.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    v
}

/// Equal-count bins over `n` sorted pairs: `K = clamp(n / 25, 4, 60)`, at most
/// `n`, each as an index range.
fn bins(n: usize) -> Vec<std::ops::Range<usize>> {
    let k = (n / 25).clamp(4, 60).min(n);
    (0..k).map(|b| b * n / k..(b + 1) * n / k).collect()
}

/// Fits one run's map to its sorted, non-empty `(t_p, rt)` pairs.
fn fit_map(pairs: &[(f64, f64)]) -> RunMap {
    let n = pairs.len();
    debug_assert!(n > 0);
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    // (x, y, weight) per bin; bins with the same median x are pooled.
    let mut knots: Vec<(f64, f64, f64)> = Vec::new();
    for b in bins(n) {
        let w = b.len() as f64;
        xs.clear();
        ys.clear();
        xs.extend(pairs[b.clone()].iter().map(|q| q.0));
        ys.extend(pairs[b].iter().map(|q| q.1));
        let x = median_in_place(&mut xs);
        let y = median_in_place(&mut ys);
        match knots.last_mut() {
            Some(last) if last.0 == x => {
                last.1 = (last.1 * last.2 + y * w) / (last.2 + w);
                last.2 += w;
            }
            _ => knots.push((x, y, w)),
        }
    }
    // Pool-adjacent-violators: blocks of (sum w*x, sum w*y, sum w), pooled
    // until the block means of y are strictly increasing.
    let mut blocks: Vec<(f64, f64, f64)> = Vec::with_capacity(knots.len());
    for &(x, y, w) in &knots {
        blocks.push((w * x, w * y, w));
        while blocks.len() >= 2 {
            let b = blocks[blocks.len() - 1];
            let a = blocks[blocks.len() - 2];
            if a.1 / a.2 < b.1 / b.2 {
                break;
            }
            blocks.pop();
            let last = blocks.last_mut().expect("two blocks");
            *last = (a.0 + b.0, a.1 + b.1, a.2 + b.2);
        }
    }
    let knots_x: Vec<f64> = blocks.iter().map(|b| b.0 / b.2).collect();
    let knots_y: Vec<f64> = blocks.iter().map(|b| b.1 / b.2).collect();

    let k = knots_x.len();
    let chord = if k >= 2 {
        (knots_y[k - 1] - knots_y[0]) / (knots_x[k - 1] - knots_x[0])
    } else {
        1.0
    };
    let m = n.div_ceil(5).max(2).min(n);
    let slope = |s: &[(f64, f64)]| {
        let v = ls_slope(s);
        if v.is_finite() && v > 0.0 {
            v
        } else {
            chord
        }
    };
    RunMap {
        slope_lo: slope(&pairs[..m]),
        slope_hi: slope(&pairs[n - m..]),
        knots_x,
        knots_y,
    }
}

/// Least-squares slope of y on x; NaN when x has no spread.
fn ls_slope(pairs: &[(f64, f64)]) -> f64 {
    let n = pairs.len() as f64;
    let mx = pairs.iter().map(|q| q.0).sum::<f64>() / n;
    let my = pairs.iter().map(|q| q.1).sum::<f64>() / n;
    let (mut sxx, mut sxy) = (0.0, 0.0);
    for &(x, y) in pairs {
        sxx += (x - mx) * (x - mx);
        sxy += (x - mx) * (y - my);
    }
    if sxx > 0.0 {
        sxy / sxx
    } else {
        f64::NAN
    }
}

/// Mean of the maps, evaluated on the union of their knots (exact for
/// piecewise-linear maps); `None` without maps.
fn mean_map<'a>(maps: impl Iterator<Item = &'a RunMap>) -> Option<RunMap> {
    let ms: Vec<&RunMap> = maps.collect();
    if ms.is_empty() {
        return None;
    }
    let k = ms.len() as f64;
    let mut grid: Vec<f64> = ms.iter().flat_map(|m| m.knots_x.iter().copied()).collect();
    grid.sort_by(f64::total_cmp);
    grid.dedup();
    let knots_y = grid
        .iter()
        .map(|&x| ms.iter().map(|m| m.apply(x)).sum::<f64>() / k)
        .collect();
    Some(RunMap {
        knots_x: grid,
        knots_y,
        slope_lo: ms.iter().map(|m| m.slope_lo).sum::<f64>() / k,
        slope_hi: ms.iter().map(|m| m.slope_hi).sum::<f64>() / k,
    })
}

/// Leave-one-out residuals of one precursor's observations `(run, rt)` (at
/// least two), given their consensus-scale values `inv = f_r^{-1}(rt)`:
/// `(rt - f_r(median of the other inv)) / sqrt(1 + 1 / (n - 1))`.
fn loo_residuals(maps: &[Option<RunMap>], obs: &[(u16, f64)], inv: &[f64]) -> Vec<f64> {
    let n = inv.len();
    debug_assert!(n >= 2 && obs.len() == n);
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| inv[a].total_cmp(&inv[b]).then(a.cmp(&b)));
    let mut rank = vec![0usize; n];
    for (k, &i) in order.iter().enumerate() {
        rank[i] = k;
    }
    let m = n - 1;
    let scale = (1.0 + 1.0 / m as f64).sqrt();
    (0..n)
        .map(|i| {
            let k = rank[i];
            // j-th smallest of the others.
            let get = |j: usize| inv[order[if j < k { j } else { j + 1 }]];
            let t_loo = if m % 2 == 1 {
                get(m / 2)
            } else {
                0.5 * (get(m / 2 - 1) + get(m / 2))
            };
            let (r, rt) = obs[i];
            let fitted = maps[r as usize]
                .as_ref()
                .map_or(f64::NAN, |mp| mp.apply(t_loo));
            (rt - fitted) / scale
        })
        .collect()
}

/// Residual width knots `(median t_p, 1.4826 * MAD of e)` over equal-count
/// bins of `pairs = (t_p, e)` sorted by `t_p`, smoothed by a running median of
/// 3 (ends kept) and floored at `floor`.
fn residual_width(pairs: &[(f64, f64)], floor: f64) -> Vec<(f64, f64)> {
    let mut xs = Vec::new();
    let mut res = Vec::new();
    let mut knots: Vec<(f64, f64)> = bins(pairs.len())
        .into_iter()
        .map(|b| {
            xs.clear();
            res.clear();
            xs.extend(pairs[b.clone()].iter().map(|q| q.0));
            res.extend(pairs[b].iter().map(|q| q.1));
            let x = median_in_place(&mut xs);
            let med = median_in_place(&mut res);
            for e in res.iter_mut() {
                *e = (*e - med).abs();
            }
            (x, 1.4826 * median_in_place(&mut res))
        })
        .collect();
    let raw: Vec<f64> = knots.iter().map(|k| k.1).collect();
    for i in 1..raw.len().saturating_sub(1) {
        let mut w = [raw[i - 1], raw[i], raw[i + 1]];
        knots[i].1 = median_in_place(&mut w);
    }
    for k in knots.iter_mut() {
        k.1 = k.1.max(floor);
    }
    knots
}

/// Median (mean of the middle two for an even count); sorts `v`. NaN when
/// empty.
fn median_in_place(v: &mut [f64]) -> f64 {
    let n = v.len();
    if n == 0 {
        return f64::NAN;
    }
    v.sort_by(f64::total_cmp);
    if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Rng(u64);

    impl Rng {
        /// splitmix64.
        fn next_u64(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        }

        /// Uniform in (0, 1].
        fn uniform(&mut self) -> f64 {
            ((self.next_u64() >> 11) as f64 + 1.0) / (1u64 << 53) as f64
        }

        /// Standard normal by Box-Muller.
        fn normal(&mut self) -> f64 {
            let u1 = self.uniform();
            let u2 = self.uniform();
            (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
        }
    }

    /// True `t_p` per precursor and every `(p, r)` observation.
    fn synth_truth(
        n_p: usize,
        warps: &[fn(f64) -> f64],
        noise: f64,
        seed: u64,
    ) -> (Vec<f64>, Vec<Observation>) {
        let mut rng = Rng(seed);
        let truth: Vec<f64> = (0..n_p).map(|_| 600.0 + 3000.0 * rng.uniform()).collect();
        let mut obs = Vec::with_capacity(n_p * warps.len());
        for (p, &t) in truth.iter().enumerate() {
            for (r, w) in warps.iter().enumerate() {
                obs.push(Observation {
                    precursor: p as u32,
                    run: r as u16,
                    rt: w(t) + noise * rng.normal(),
                });
            }
        }
        (truth, obs)
    }

    fn synth(n_p: usize, warps: &[fn(f64) -> f64], noise: f64, seed: u64) -> Vec<Observation> {
        synth_truth(n_p, warps, noise, seed).1
    }

    fn median(v: &[f64]) -> f64 {
        let mut v = v.to_vec();
        median_in_place(&mut v)
    }

    fn three_warps() -> [fn(f64) -> f64; 3] {
        [
            |t| t,
            |t| 1.1 * t + 30.0,
            |t| t + 60.0 * (t / 3600.0).powi(2),
        ]
    }

    #[test]
    fn recovers_monotone_warps() {
        let warps = three_warps();
        let obs = synth(2000, &warps, 5.0, 1);
        let a = align(&obs, 2000, 3, &AlignParams::default());
        let errs: Vec<f64> = obs
            .iter()
            .map(|o| (a.predict(o.precursor, o.run).unwrap().0 - o.rt).abs())
            .collect();
        let m = median(&errs);
        // Held out: the run's own observation does not enter its prediction.
        let held: Vec<f64> = obs
            .iter()
            .map(|o| (a.predict_excluding(o.precursor, o.run).unwrap().0 - o.rt).abs())
            .collect();
        let mh = median(&held);
        eprintln!("recovers_monotone_warps: median |t_hat - rt| = {m:.3}, held out = {mh:.3}");
        let sds: Vec<f64> = a.residual_sd[0]
            .as_ref()
            .unwrap()
            .iter()
            .map(|k| k.1)
            .collect();
        eprintln!(
            "recovers_monotone_warps: run 0 median sigma_r = {:.3} (noise 5)",
            median(&sds)
        );
        assert!(m < 7.5, "{m}");
        assert!(mh < 7.5, "{mh}");
    }

    #[test]
    fn robust_to_ten_percent_outliers() {
        let warps = three_warps();
        let (truth, mut obs) = synth_truth(2000, &warps, 5.0, 2);
        let mut rng = Rng(99);
        let mut is_outlier = vec![false; obs.len()];
        for (i, o) in obs.iter_mut().enumerate() {
            if rng.uniform() <= 0.10 {
                o.rt = 4000.0 * rng.uniform();
                is_outlier[i] = true;
            }
        }
        let a = align(&obs, 2000, 3, &AlignParams::default());
        let errs: Vec<f64> = obs
            .iter()
            .map(|o| (a.predict(o.precursor, o.run).unwrap().0 - o.rt).abs())
            .collect();
        let m = median(&errs);
        // Against the noise-free truth, over the observations that are not outliers.
        let errs_truth: Vec<f64> = obs
            .iter()
            .zip(&is_outlier)
            .filter(|(_, &out)| !out)
            .map(|(o, _)| {
                let t = truth[o.precursor as usize];
                (a.predict(o.precursor, o.run).unwrap().0 - warps[o.run as usize](t)).abs()
            })
            .collect();
        let mt = median(&errs_truth);
        eprintln!(
            "robust_to_ten_percent_outliers: median |t_hat - rt| = {m:.3}, \
             median |t_hat - truth| (inliers) = {mt:.3}"
        );
        assert!(m < 7.5, "{m}");
        assert!(mt < 5.0, "{mt}");
    }

    #[test]
    fn maps_are_monotone_and_invertible() {
        let warps = three_warps();
        let mut obs = synth(2000, &warps, 5.0, 3);
        let mut rng = Rng(7);
        for o in obs.iter_mut() {
            if rng.uniform() <= 0.10 {
                o.rt = 4000.0 * rng.uniform();
            }
        }
        let a = align(&obs, 2000, 3, &AlignParams::default());
        let mut n_aligned = 0;
        for m in a.maps.iter().flatten() {
            n_aligned += 1;
            let mut prev = f64::NEG_INFINITY;
            for i in 0..=5000 {
                let t = -500.0 + i as f64;
                let y = m.apply(t);
                assert!(y.is_finite());
                assert!(y >= prev, "apply decreases at {t}: {prev} -> {y}");
                prev = y;
                let back = m.invert(y);
                assert!((back - t).abs() < 1e-6, "invert(apply({t})) = {back}");
            }
        }
        assert_eq!(n_aligned, 3);
    }

    fn bits(v: Option<(f64, f64)>) -> Option<(u64, u64)> {
        v.map(|(a, b)| (a.to_bits(), b.to_bits()))
    }

    #[test]
    fn rebuilt_alignment_predicts_bit_for_bit() {
        let warps = three_warps();
        let mut obs = synth(600, &warps, 5.0, 8);
        // An unaligned fourth run and precursors seen in one run only.
        obs.extend((0..10u32).map(|p| Observation {
            precursor: p,
            run: 3,
            rt: 50.0 + p as f64,
        }));
        obs.retain(|o| !(o.precursor >= 590 && o.run > 0));
        let a = align(&obs, 610, 4, &AlignParams::default());
        assert!(a.maps[3].is_none());
        let runs = a.runs();
        assert!(runs[3].observations.is_empty());
        let b = Alignment::from_runs(runs.clone(), 610).unwrap();
        assert_eq!(b.runs(), runs);
        assert_eq!(a.n_runs_seen, b.n_runs_seen);
        assert_eq!(a.anchors, b.anchors);
        for p in 0..610u32 {
            let i = p as usize;
            assert_eq!(a.consensus[i].to_bits(), b.consensus[i].to_bits(), "{p}");
            assert_eq!(a.tau[i].to_bits(), b.tau[i].to_bits(), "{p}");
            for r in 0..4u16 {
                assert_eq!(bits(a.predict(p, r)), bits(b.predict(p, r)), "{p} {r}");
                assert_eq!(
                    bits(a.predict_excluding(p, r)),
                    bits(b.predict_excluding(p, r)),
                    "{p} {r}"
                );
            }
        }
    }

    #[test]
    fn added_observations_are_predicted_without_moving_the_fit() {
        let warps = three_warps();
        let (truth, mut obs) = synth_truth(620, &warps, 5.0, 9);
        // Precursors 600..620 are observed after the fit, in two runs each,
        // at RTs far off their truth: anchors would bend the maps.
        let extra: Vec<Observation> = obs
            .iter()
            .filter(|o| o.precursor >= 600 && o.run < 2)
            .map(|o| Observation {
                rt: o.rt + 400.0,
                ..*o
            })
            .collect();
        obs.retain(|o| o.precursor < 600);
        let fitted = align(&obs, 620, 3, &AlignParams::default());
        let mut a = fitted.clone();
        a.add_observations(&extra).unwrap();
        for r in 0..3 {
            let (x, y) = (fitted.maps[r].as_ref(), a.maps[r].as_ref());
            assert_eq!(x.map(|m| m.knots()), y.map(|m| m.knots()));
        }
        assert_eq!(fitted.anchors, a.anchors);
        for p in 0..600u32 {
            for r in 0..3u16 {
                assert_eq!(bits(fitted.predict(p, r)), bits(a.predict(p, r)));
            }
        }
        // An added precursor is predicted from its own observations as a
        // fitted one would be: the median of its RTs carried onto the
        // consensus scale, mapped into run 2.
        let map = |r: usize| a.maps[r].as_ref().unwrap();
        for p in 600..620u32 {
            assert!(fitted.predict(p, 2).is_none());
            assert_eq!(a.n_runs_seen[p as usize], 2);
            let rt = |r: u16| {
                extra
                    .iter()
                    .find(|o| o.precursor == p && o.run == r)
                    .unwrap()
                    .rt
            };
            let tc = 0.5 * (map(0).invert(rt(0)) + map(1).invert(rt(1)));
            let (t, s) = a.predict(p, 2).unwrap();
            assert!((t - map(2).apply(tc)).abs() < 1e-6, "{p}");
            assert!(t > warps[2](truth[p as usize]) + 200.0, "{p}");
            assert!(s.is_finite() && s > 0.0);
            assert!(a.predict_excluding(p, 0).is_some());
        }
        // Round trip, and a second observation of a pair is refused.
        let b = Alignment::from_runs(a.runs(), 620).unwrap();
        assert_eq!(bits(a.predict(605, 2)), bits(b.predict(605, 2)));
        assert!(a.add_observations(&extra[..1]).is_err());
    }

    #[test]
    fn from_runs_rejects_inconsistent_state() {
        let a = align(
            &synth(300, &three_warps(), 5.0, 2),
            300,
            3,
            &AlignParams::default(),
        );
        let mut runs = a.runs();
        let first = runs[0].observations[0];
        runs[0].observations.push(first);
        let e = Alignment::from_runs(runs, 300).unwrap_err();
        assert!(e.contains("observed twice"), "{e}");
        let e = Alignment::from_runs(a.runs(), 10).unwrap_err();
        assert!(e.contains("outside"), "{e}");
        let mut runs = a.runs();
        runs[1].map = None;
        assert!(Alignment::from_runs(runs, 300).is_err());
        assert!(RunMap::from_parts(vec![], vec![], 1.0, 1.0).is_err());
        assert!(RunMap::from_parts(vec![1.0], vec![2.0], 0.0, 1.0).is_err());
        assert!(RunMap::from_parts(vec![1.0], vec![f64::NAN], 1.0, 1.0).is_err());
    }

    #[test]
    fn run_with_few_anchors_is_unaligned() {
        let warps = three_warps();
        let mut obs = synth(2000, &warps, 5.0, 4);
        // Run 3 sees only 20 of the precursors.
        for p in 0..20u32 {
            obs.push(Observation {
                precursor: p,
                run: 3,
                rt: 1000.0 + 10.0 * p as f64,
            });
        }
        let a = align(&obs, 2000, 4, &AlignParams::default());
        assert!(a.maps[3].is_none());
        assert!(a.residual_sd[3].is_none());
        assert_eq!(a.anchors[3], 20);
        for p in 0..2000u32 {
            assert!(a.predict(p, 3).is_none());
            assert!(a.predict_excluding(p, 3).is_none());
        }
        // Run 3's observations do not count towards n_p.
        assert_eq!(a.n_runs_seen[0], 3);
        let errs: Vec<f64> = obs
            .iter()
            .filter(|o| o.run < 3)
            .map(|o| (a.predict(o.precursor, o.run).unwrap().0 - o.rt).abs())
            .collect();
        let m = median(&errs);
        assert!(m < 7.5, "{m}");
    }

    #[test]
    fn run_with_narrow_rt_range_is_unaligned() {
        let warps = three_warps();
        let (truth, mut obs) = synth_truth(2000, &warps, 5.0, 5);
        // Run 3 sees 200 precursors, all eluting in the first ~7% of the gradient.
        let mut added = 0;
        for (p, &t) in truth.iter().enumerate() {
            if t < 800.0 && added < 200 {
                obs.push(Observation {
                    precursor: p as u32,
                    run: 3,
                    rt: t,
                });
                added += 1;
            }
        }
        assert!(added >= 100);
        let a = align(&obs, 2000, 4, &AlignParams::default());
        assert!(a.maps[3].is_none());
        assert!(a.coverage[3] < 0.30, "{}", a.coverage[3]);
        assert!(a.maps[..3].iter().all(Option::is_some));
    }

    #[test]
    fn single_run_precursor_gets_consensus_without_anchoring() {
        let warps = three_warps();
        let mut obs = synth(2000, &warps, 5.0, 6);
        let base = align(&obs, 2001, 3, &AlignParams::default());
        obs.push(Observation {
            precursor: 2000,
            run: 1,
            rt: 2000.0,
        });
        let a = align(&obs, 2001, 3, &AlignParams::default());
        assert!(base.consensus[2000].is_nan());
        assert!(base.predict(2000, 0).is_none());
        assert!(a.consensus[2000].is_finite());
        assert_eq!(a.n_runs_seen[2000], 1);
        assert_eq!(a.maps, base.maps);
        assert_eq!(a.residual_sd, base.residual_sd);
        for p in 0..2000 {
            assert_eq!(a.consensus[p].to_bits(), base.consensus[p].to_bits());
        }
        // Its own run maps it back onto the observed RT; other runs get a transfer.
        let (t1, _) = a.predict(2000, 1).unwrap();
        assert!((t1 - 2000.0).abs() < 1e-6, "{t1}");
        assert!(a.predict(2000, 0).is_some());
        assert!(a.predict_excluding(2000, 1).is_none());
    }

    #[test]
    fn transfer_sd_shrinks_with_more_donors() {
        let warps: [fn(f64) -> f64; 7] = [
            |t| t,
            |t| 1.05 * t + 20.0,
            |t| 0.95 * t - 10.0,
            |t| t + 40.0,
            |t| t + 60.0 * (t / 3600.0).powi(2),
            |t| 1.02 * t,
            |t| t - 30.0,
        ];
        let mut obs = synth(1500, &warps, 5.0, 8);
        // Precursor 1500 is seen in runs 1..=2, precursor 1501 in runs 1..=6, both at t = 2000.
        for (r, w) in warps.iter().enumerate().take(3).skip(1) {
            obs.push(Observation {
                precursor: 1500,
                run: r as u16,
                rt: w(2000.0),
            });
        }
        for (r, w) in warps.iter().enumerate().skip(1) {
            obs.push(Observation {
                precursor: 1501,
                run: r as u16,
                rt: w(2000.0),
            });
        }
        let a = align(&obs, 1502, 7, &AlignParams::default());
        assert_eq!(a.n_runs_seen[1500], 2);
        assert_eq!(a.n_runs_seen[1501], 6);
        let (t2, s2) = a.predict(1500, 0).unwrap();
        let (t6, s6) = a.predict(1501, 0).unwrap();
        eprintln!(
            "transfer_sd_shrinks_with_more_donors: s(2 donors) = {s2:.3}, s(6 donors) = {s6:.3}"
        );
        assert!(
            (t2 - 2000.0).abs() < 5.0 && (t6 - 2000.0).abs() < 5.0,
            "{t2} {t6}"
        );
        assert!(s2 > s6, "{s2} <= {s6}");
    }

    /// Four runs over 300 shared anchors: run 0 has PSM RT jitter of +-20 or
    /// +-40 s, runs 1 to 3 are exact. Precursor 300 is seen only in run 0,
    /// 47 s after its true apex at 2000 s; precursor 301 is seen in run 0 at
    /// 2047 s and in run 1 at 2000 s. Returns the alignment and the noisy
    /// run's median residual width.
    fn noisy_singleton_donor() -> (Alignment, f64) {
        let mut rng = Rng(21);
        let mut obs = Vec::new();
        for p in 0..300u32 {
            let t = 600.0 + 3000.0 * rng.uniform();
            for r in 0..4u16 {
                let jitter = if r == 0 {
                    [-40.0, -20.0, 20.0, 40.0][(rng.next_u64() % 4) as usize]
                } else {
                    0.0
                };
                obs.push(Observation {
                    precursor: p,
                    run: r,
                    rt: t + jitter,
                });
            }
        }
        for (precursor, run, rt) in [(300, 0, 2047.0), (301, 0, 2047.0), (301, 1, 2000.0)] {
            obs.push(Observation { precursor, run, rt });
        }
        let a = align(&obs, 302, 4, &AlignParams::default());
        let sds: Vec<f64> = a.residual_sd[0]
            .as_ref()
            .unwrap()
            .iter()
            .map(|k| k.1)
            .collect();
        let sd0 = median(&sds);
        (a, sd0)
    }

    #[test]
    fn transfer_window_includes_the_donor_uncertainty() {
        use crate::decoys::{Hypothesis, Probe};
        use crate::evidence::{fit, EvidenceParams, Prior};
        use crate::ms1_index::{Ms1RunIndex, Ms1Scan};
        use model::mass::ISOTOPE;

        let (a, sd0) = noisy_singleton_donor();
        let (t_hat, s) = a.predict(300, 1).unwrap();
        eprintln!(
            "transfer_window_includes_the_donor_uncertainty: donor sigma = {sd0:.2}, \
             t_hat = {t_hat:.2}, s = {s:.3}"
        );
        let probe = Probe {
            kind: Hypothesis::Target,
            mz: 600.0,
            t_hat,
            s,
            w_prior: 5.0,
        };
        let prior = Prior::of(&probe);

        // A strong peak at the true apex in the recipient is recovered.
        let env = [0.55f32, 0.30, 0.11, 0.04];
        let scans: Vec<Ms1Scan> = (1700..2300)
            .map(|k| {
                let t = k as f64;
                let g = (-(t - 2000.0).powi(2) / (2.0 * 25.0)).exp();
                let peaks = (0..4)
                    .map(|i| {
                        (
                            600.0 + i as f64 * ISOTOPE / 2.0,
                            (100.0 + 1e6 * env[i] as f64 * g) as f32,
                        )
                    })
                    .collect();
                Ms1Scan {
                    rt: t,
                    peaks,
                    ..Ms1Scan::default()
                }
            })
            .collect();
        let run = Ms1RunIndex::new(scans);
        let e = fit(&run, &probe, 2, &env, &EvidenceParams::default()).unwrap();
        eprintln!(
            "transfer_window_includes_the_donor_uncertainty: quantity = {:.0}, t0 = {:.2}",
            e.quantity, e.t0
        );
        assert!(e.quantity > 0.0, "{e:?}");
        assert!((e.t0 - 2000.0).abs() <= 1.0, "{e:?}");
        // The true apex lies inside the search and extraction windows, and
        // the donor's own width enters the transfer width.
        assert!((2000.0 - t_hat).abs() <= 3.0 * s, "{t_hat} {s}");
        assert!(t_hat - prior.half_window() < 2000.0);
        assert!(sd0 > 20.0, "{sd0}");
        assert!(s >= 0.8 * sd0, "s = {s}, donor sigma = {sd0}");
    }

    #[test]
    fn excluded_prediction_carries_the_donor_uncertainty() {
        let (a, sd0) = noisy_singleton_donor();
        // Anchor 0 is seen in every run: without exact run 1, two exact
        // donors and the noisy one remain. The noisy donor's width enters,
        // shrunk by the three donors.
        let (_, s) = a.predict_excluding(0, 1).unwrap();
        eprintln!(
            "excluded_prediction_carries_the_donor_uncertainty: s = {s:.3}, \
             donor sigma = {sd0:.2}"
        );
        assert!(s > 0.2 * sd0 && s < 0.6 * sd0, "{s}");
        // Precursor 301 without run 1 has only the noisy donor left.
        let (_, s) = a.predict_excluding(301, 1).unwrap();
        assert!(s >= 0.8 * sd0, "{s}");
        // Predicting the singleton into a run without it excludes nothing.
        assert_eq!(a.predict_excluding(300, 1), a.predict(300, 1));
    }

    /// Runs 0 and 1 share precursors 0..100, runs 2 and 3 share 100..200,
    /// and no precursor is seen in both pairs; runs 2 and 3 elute 120 s
    /// later. Precursor 200 is seen in run 0 only, 201 in run 2 only.
    fn two_pairs() -> Alignment {
        let mut rng = Rng(31);
        let mut obs = Vec::new();
        for p in 0..200u32 {
            let t = 600.0 + 3000.0 * rng.uniform();
            let (runs, shift) = if p < 100 {
                (0..2u16, 0.0)
            } else {
                (2..4u16, 120.0)
            };
            for r in runs {
                obs.push(Observation {
                    precursor: p,
                    run: r,
                    rt: t + shift + 0.5 * rng.normal(),
                });
            }
        }
        obs.push(Observation {
            precursor: 200,
            run: 0,
            rt: 1500.0,
        });
        obs.push(Observation {
            precursor: 201,
            run: 2,
            rt: 1620.0,
        });
        align(&obs, 202, 4, &AlignParams::default())
    }

    #[test]
    fn no_transfer_between_disconnected_runs() {
        let a = two_pairs();
        assert!(a.maps.iter().all(Option::is_some));
        assert_eq!(a.component(0), a.component(1));
        assert_eq!(a.component(2), a.component(3));
        assert_ne!(a.component(0), a.component(2));
        if let Some((t, s)) = a.predict(200, 2) {
            panic!("precursor 200 transferred across components: {t:.2} +- {s:.3}");
        }
        for (p, r) in [(200u32, 3u16), (201, 0), (201, 1), (5, 2), (150, 0)] {
            assert!(a.predict(p, r).is_none(), "{p} into {r}");
            assert!(a.predict_excluding(p, r).is_none(), "{p} into {r}");
            assert!(a.blocked_by_component(p, r, None), "{p} into {r}");
        }
        // Within a pair transfers still work.
        let (t, _) = a.predict(200, 1).unwrap();
        assert!((t - 1500.0).abs() < 3.0, "{t}");
        let (t, _) = a.predict(201, 3).unwrap();
        assert!((t - 1620.0).abs() < 3.0, "{t}");
        assert!(a.predict_excluding(5, 1).is_some());
        assert!(a.predict_excluding(150, 3).is_some());
        assert!(!a.blocked_by_component(200, 1, None));
        assert!(!a.blocked_by_component(5, 1, Some(1)));
    }

    #[test]
    fn precursor_seen_in_two_components_uses_the_destination_component() {
        // As `two_pairs`, with 600 precursors per pair.
        let mut rng = Rng(31);
        let mut obs = Vec::new();
        for p in 0..1200u32 {
            let t = 600.0 + 3000.0 * rng.uniform();
            let (runs, shift) = if p < 600 {
                (0..2u16, 0.0)
            } else {
                (2..4u16, 120.0)
            };
            for r in runs {
                obs.push(Observation {
                    precursor: p,
                    run: r,
                    rt: t + shift + 0.5 * rng.normal(),
                });
            }
        }
        // Precursor 1200 elutes at 1500 s in the first pair, 1620 s in the
        // second.
        obs.push(Observation {
            precursor: 1200,
            run: 0,
            rt: 1500.0,
        });
        obs.push(Observation {
            precursor: 1200,
            run: 2,
            rt: 1620.0,
        });
        let a = align(&obs, 1201, 4, &AlignParams::default());
        assert_eq!(a.n_components(), 2);
        let (t1, _) = a.predict(1200, 1).unwrap();
        let (t3, _) = a.predict(1200, 3).unwrap();
        eprintln!(
            "precursor_seen_in_two_components_uses_the_destination_component: \
             run 1 {t1:.2} (1500), run 3 {t3:.2} (1620)"
        );
        assert!((t1 - 1500.0).abs() < 3.0, "{t1}");
        assert!((t3 - 1620.0).abs() < 3.0, "{t3}");
        // Without run 2, the second pair has no donor.
        assert!(a.predict_excluding(1200, 2).is_none());
        assert!(a.blocked_by_component(1200, 2, Some(2)));
    }

    #[test]
    fn run_sharing_too_few_precursors_with_any_run_is_unaligned() {
        // Runs 0/1 and 2/3 as in `two_pairs`; run 4 sees 12 precursors of
        // each pair, so it shares 12 with every run.
        let mut rng = Rng(41);
        let mut obs = Vec::new();
        for p in 0..200u32 {
            let t = 600.0 + 3000.0 * rng.uniform();
            let mut runs: Vec<u16> = if p < 100 { vec![0, 1] } else { vec![2, 3] };
            if p % 100 < 12 {
                runs.push(4);
            }
            for r in runs {
                obs.push(Observation {
                    precursor: p,
                    run: r,
                    rt: t + 0.5 * rng.normal(),
                });
            }
        }
        let ap = AlignParams {
            min_anchors: 20,
            ..AlignParams::default()
        };
        let a = align(&obs, 200, 5, &ap);
        assert_eq!(a.anchors[4], 24);
        assert!(a.maps[4].is_none() && a.component(4).is_none());
        assert!(a.maps[..4].iter().all(Option::is_some));
        assert_eq!(a.n_components(), 2);
        assert!(a.predict(0, 4).is_none());
        assert_eq!(a.n_runs_seen[0], 2);
    }

    #[test]
    fn predict_excluding_ignores_own_observation() {
        let warps = three_warps();
        let (truth, mut obs) = synth_truth(2000, &warps, 1.0, 9);
        // Precursor 2000 is seen in runs 0 and 1; run 0's RT is 300 s off.
        let t = 2000.0;
        obs.push(Observation {
            precursor: 2000,
            run: 0,
            rt: warps[0](t) + 300.0,
        });
        obs.push(Observation {
            precursor: 2000,
            run: 1,
            rt: warps[1](t),
        });
        let a = align(&obs, 2001, 3, &AlignParams::default());
        let (own, _) = a.predict(2000, 0).unwrap();
        let (excl, s_excl) = a.predict_excluding(2000, 0).unwrap();
        eprintln!(
            "predict_excluding_ignores_own_observation: predict = {own:.3}, \
             predict_excluding = {excl:.3}, truth = {:.3}",
            warps[0](t)
        );
        assert!((own - warps[0](t)).abs() > 100.0, "{own}");
        assert!((excl - warps[0](t)).abs() < 3.0, "{excl}");
        // Excluding run 0 leaves a single donor, so no extra dispersion from the outlier.
        let (_, s_other) = a.predict(1, 0).unwrap();
        assert!(s_excl < 3.0 * s_other.max(1.0), "{s_excl}");
        // An anchor seen in all three runs: excluding its own run stays on the truth.
        let p = 5u32;
        let (e, _) = a.predict_excluding(p, 2).unwrap();
        assert!((e - warps[2](truth[p as usize])).abs() < 5.0, "{e}");
        // Excluding a run the precursor is not seen in equals predict.
        assert!(a.predict_excluding(2000, 2).is_some());
        assert_eq!(a.predict_excluding(2000, 2), a.predict(2000, 2));
    }

    #[test]
    fn mean_map_is_the_identity_after_the_gauge_fix() {
        let warps = three_warps();
        let obs = synth(2000, &warps, 5.0, 10);
        let a = align(&obs, 2000, 3, &AlignParams::default());
        let maps: Vec<&RunMap> = a.maps.iter().flatten().collect();
        assert_eq!(maps.len(), 3);
        // g(t) - t on the consensus range: no offset or scale drift, only the
        // bin-median noise of the maps.
        let dev: Vec<f64> = (0..=30)
            .map(|i| {
                let t = 700.0 + 100.0 * i as f64;
                maps.iter().map(|m| m.apply(t)).sum::<f64>() / maps.len() as f64 - t
            })
            .collect();
        let mean = dev.iter().sum::<f64>() / dev.len() as f64;
        let worst = dev.iter().fold(0.0f64, |w, d| w.max(d.abs()));
        eprintln!(
            "mean_map_is_the_identity_after_the_gauge_fix: mean (g(t) - t) = {mean:.3}, \
             max |g(t) - t| = {worst:.3}"
        );
        assert!(mean.abs() < 0.5, "{mean}");
        assert!(worst < 4.0, "{worst}");
    }

    #[test]
    fn sigma_tracks_true_noise() {
        let three = three_warps();
        let six: [fn(f64) -> f64; 6] = [
            |t| t,
            |t| 1.1 * t + 30.0,
            |t| t + 60.0 * (t / 3600.0).powi(2),
            |t| 0.95 * t - 10.0,
            |t| t + 40.0,
            |t| 1.02 * t,
        ];
        for (name, warps) in [("3 runs", &three[..]), ("6 runs", &six[..])] {
            let obs = synth(2000, warps, 5.0, 12);
            let a = align(&obs, 2000, warps.len(), &AlignParams::default());
            for (r, knots) in a.residual_sd.iter().enumerate() {
                let sds: Vec<f64> = knots.as_ref().unwrap().iter().map(|k| k.1).collect();
                let m = median(&sds);
                eprintln!("sigma_tracks_true_noise: {name}, run {r}: median sigma_r = {m:.3}");
                assert!((m - 5.0).abs() < 1.25, "{name} run {r}: {m}");
            }
        }
    }

    #[test]
    fn output_is_independent_of_input_order() {
        let warps = three_warps();
        let obs = synth(2000, &warps, 5.0, 11);
        let a = align(&obs, 2000, 3, &AlignParams::default());
        let mut rev = obs.clone();
        rev.reverse();
        let b = align(&rev, 2000, 3, &AlignParams::default());
        assert_eq!(a.maps, b.maps);
        assert_eq!(a.residual_sd, b.residual_sd);
        assert_eq!(a.n_runs_seen, b.n_runs_seen);
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(&a.consensus), bits(&b.consensus));
        assert_eq!(bits(&a.tau), bits(&b.tau));
    }

    #[test]
    fn degenerate_inputs_do_not_panic() {
        let p = AlignParams::default();
        let a = align(&[], 0, 0, &p);
        assert!(a.maps.is_empty() && a.consensus.is_empty());
        let a = align(&[], 3, 2, &p);
        assert!(a.maps.iter().all(Option::is_none));
        assert!(a.predict(0, 0).is_none() && a.predict(9, 9).is_none());
        // One run only: nothing is seen twice, so the run is unaligned.
        let obs: Vec<Observation> = (0..100)
            .map(|i| Observation {
                precursor: i,
                run: 0,
                rt: 600.0 + i as f64,
            })
            .collect();
        let a = align(&obs, 100, 1, &p);
        assert!(a.maps[0].is_none());
        assert!(a.consensus.iter().all(|c| c.is_nan()));
        // Identical RTs everywhere and a non-finite RT.
        let mut obs: Vec<Observation> = (0..200u32)
            .flat_map(|i| {
                (0..2u16).map(move |r| Observation {
                    precursor: i,
                    run: r,
                    rt: 1000.0,
                })
            })
            .collect();
        obs.push(Observation {
            precursor: 200,
            run: 0,
            rt: f64::NAN,
        });
        let a = align(&obs, 201, 2, &p);
        assert!(a.consensus[200].is_nan());
    }
}
