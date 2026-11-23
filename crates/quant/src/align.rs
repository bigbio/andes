use model::stats::{median_mad_mut, median_mut, Even, MAD_TO_SD};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Observation {
    pub precursor: u32,
    pub run: u16,
    pub rt: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AlignParams {
    pub rounds: usize,
    pub min_anchors: usize,
    pub min_coverage: f64,
    pub min_spacing_s: f64,
    pub run_min_spacing_s: Vec<f64>,
    pub min_shared: usize,
}

impl AlignParams {
    pub fn floor_s(&self, run: usize) -> f64 {
        self.run_min_spacing_s
            .get(run)
            .copied()
            .filter(|f| f.is_finite() && *f > 0.0)
            .unwrap_or(self.min_spacing_s)
    }
}

impl Default for AlignParams {
    fn default() -> Self {
        Self {
            rounds: 8,
            min_anchors: 50,
            min_coverage: 0.30,
            min_spacing_s: 0.5,
            run_min_spacing_s: Vec::new(),
            min_shared: 20,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct RunMap {
    knots_x: Vec<f64>,
    knots_y: Vec<f64>,
    slope_lo: f64,
    slope_hi: f64,
}

impl RunMap {
    pub fn apply(&self, t: f64) -> f64 {
        interp(
            &self.knots_x,
            &self.knots_y,
            self.slope_lo,
            self.slope_hi,
            t,
        )
    }

    pub fn invert(&self, rt: f64) -> f64 {
        interp(
            &self.knots_y,
            &self.knots_x,
            1.0 / self.slope_lo,
            1.0 / self.slope_hi,
            rt,
        )
    }

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
        let i = self.knots_x.partition_point(|&k| k <= t);
        let (a, b) = (i.saturating_sub(2), (i + 1).min(n - 1));
        (self.knots_y[b] - self.knots_y[a]) / (self.knots_x[b] - self.knots_x[a])
    }

    pub fn knots(&self) -> (&[f64], &[f64]) {
        (&self.knots_x, &self.knots_y)
    }

    pub fn end_slopes(&self) -> (f64, f64) {
        (self.slope_lo, self.slope_hi)
    }

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

#[derive(Debug, Clone, PartialEq)]
pub struct RunAlignment {
    pub map: Option<RunMap>,
    pub residual_sd: Option<Vec<(f64, f64)>>,
    pub anchors: usize,
    pub coverage: f64,
    pub component: Option<u32>,
    pub observations: Vec<(u32, f64)>,
}

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
    let i = xs.partition_point(|&k| k <= t);
    let (x0, x1, y0, y1) = (xs[i - 1], xs[i], ys[i - 1], ys[i]);
    y0 + (y1 - y0) * (t - x0) / (x1 - x0)
}

#[derive(Debug, Clone)]
pub struct Alignment {
    pub consensus: Vec<f64>,
    pub maps: Vec<Option<RunMap>>,
    pub residual_sd: Vec<Option<Vec<(f64, f64)>>>,
    pub tau: Vec<f64>,
    pub n_runs_seen: Vec<u16>,
    pub anchors: Vec<usize>,
    pub coverage: Vec<f64>,
    components: Vec<Option<u32>>,
    home: Vec<u32>,
    consensus_var: Vec<f64>,
    obs_offsets: Vec<usize>,
    obs_runs: Vec<u16>,
    obs_rts: Vec<f64>,
}

impl Alignment {
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

    pub fn predict_excluding(&self, p: u32, run: u16) -> Option<(f64, f64)> {
        let map = self.maps.get(run as usize)?.as_ref()?;
        let c = self.component(run)?;
        if p as usize >= self.consensus.len() {
            return None;
        }
        let e = self.estimate(p as usize, Some(run), c)?;
        Some(self.transfer(map, run, e.t, e.var, e.tau))
    }

    pub fn component(&self, run: u16) -> Option<u32> {
        self.components.get(run as usize).copied().flatten()
    }

    pub fn n_components(&self) -> usize {
        self.components
            .iter()
            .flatten()
            .max()
            .map_or(0, |&c| c as usize + 1)
    }

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

    fn transfer(&self, map: &RunMap, run: u16, t: f64, var: f64, tau: f64) -> (f64, f64) {
        let sigma = self.sigma_at(run, t);
        let slope = map.slope_at(t);
        let s = (sigma * sigma + slope * slope * var + tau * tau).sqrt();
        (map.apply(t), s)
    }

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

struct Estimate {
    t: f64,
    tau: f64,
    var: f64,
}

pub fn align(
    obs: &[Observation],
    n_precursors: usize,
    n_runs: usize,
    p: &AlignParams,
) -> Alignment {
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
                residual_width(v, p.floor_s(r))
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

fn robust_span(v: &mut [f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    v.sort_by(f64::total_cmp);
    let q = |f: f64| v[((v.len() - 1) as f64 * f).round() as usize];
    q(0.99) - q(0.01)
}

fn sorted_pairs(pairs: &[(usize, f64)], t: &[f64]) -> Vec<(f64, f64)> {
    let mut v: Vec<(f64, f64)> = pairs.iter().map(|&(pr, rt)| (t[pr], rt)).collect();
    v.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    v
}

fn bins(n: usize) -> Vec<std::ops::Range<usize>> {
    let k = (n / 25).clamp(4, 60).min(n);
    (0..k).map(|b| b * n / k..(b + 1) * n / k).collect()
}

fn fit_map(pairs: &[(f64, f64)]) -> RunMap {
    let n = pairs.len();
    debug_assert!(n > 0);
    let mut xs = Vec::new();
    let mut ys = Vec::new();
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
            let mad = median_mad_mut(&mut res, Even::Mean).map_or(f64::NAN, |(_, mad)| mad);
            (x, MAD_TO_SD * mad)
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

fn median_in_place(v: &mut [f64]) -> f64 {
    median_mut(v, Even::Mean).unwrap_or(f64::NAN)
}

#[cfg(test)]
#[path = "align/tests/align.rs"]
mod tests;
