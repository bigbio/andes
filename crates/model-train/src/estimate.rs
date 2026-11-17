use rustc_hash::FxHashMap;

use scoring_crate::param_model::{IonType, Param, Partition};

use crate::counts::CountStats;

#[derive(Debug, Clone)]
pub struct EstimatorConfig {
    pub pseudo: f32,
    pub noise_pseudo: f32,
    pub min_count: u64,
    pub backoff_weight: f32,
    pub error_scaling_factor_override: Option<i32>,
    pub rank_smoothing: bool,
}

impl Default for EstimatorConfig {
    fn default() -> Self {
        Self {
            pseudo: 1.0,
            noise_pseudo: 0.05,
            min_count: 50,
            backoff_weight: 20.0,
            error_scaling_factor_override: None,
            rank_smoothing: false,
        }
    }
}

pub struct Estimator {
    cfg: EstimatorConfig,
}

impl Estimator {
    pub fn new(cfg: EstimatorConfig) -> Self {
        Self { cfg }
    }

    pub fn estimate(&self, counts: &CountStats, template: &Param) -> Param {
        self.estimate_with_prior(counts, template, None)
    }

    pub fn estimate_with_prior(
        &self,
        counts: &CountStats,
        template: &Param,
        prior: Option<&Param>,
    ) -> Param {
        assert!(
            self.cfg.pseudo > 0.0,
            "EstimatorConfig.pseudo must be > 0 (got {})",
            self.cfg.pseudo
        );
        assert!(
            self.cfg.noise_pseudo > 0.0,
            "EstimatorConfig.noise_pseudo must be > 0 (got {})",
            self.cfg.noise_pseudo
        );
        let max_rank = template.max_rank;
        let esf = self
            .cfg
            .error_scaling_factor_override
            .unwrap_or(template.error_scaling_factor);

        let rank_dist_table = self.build_rank_dist_table(counts, template, max_rank, prior);
        let (ion_err_dist_table, noise_err_dist_table) =
            self.build_error_tables(counts, template, esf, prior);
        let ion_existence_table = self.build_existence_table(counts, template, prior);
        let charge_hist = build_charge_hist(&counts.charge);
        let (min_charge, max_charge) = charge_range(&charge_hist, template);

        let mut partitions: Vec<Partition> = rank_dist_table.keys().copied().collect();
        partitions.sort();

        let mut param = Param {
            version: template.version,
            data_type: template.data_type.clone(),
            mme: template.mme,
            apply_deconvolution: template.apply_deconvolution,
            deconvolution_error_tolerance: template.deconvolution_error_tolerance,
            charge_hist,
            min_charge,
            max_charge,
            num_segments: template.num_segments,
            partitions,
            num_precursor_off: template.num_precursor_off,
            precursor_off_map: template.precursor_off_map.clone(),
            frag_off_table: template.frag_off_table.clone(),
            max_rank,
            rank_dist_table,
            error_scaling_factor: esf,
            ion_err_dist_table,
            noise_err_dist_table,
            ion_existence_table,
            partition_ion_types_cache: FxHashMap::default(),
            gbdt_peak_model: None,
            frag_intensity_model: None,
            rich_ion_model: None,
        };
        param.rebuild_cache();
        param
    }

    fn build_rank_dist_table(
        &self,
        counts: &CountStats,
        template: &Param,
        max_rank: i32,
        prior: Option<&Param>,
    ) -> FxHashMap<Partition, FxHashMap<IonType, Vec<f32>>> {
        let n_slots = (max_rank + 1) as usize;
        let pseudo = self.cfg.pseudo;
        let noise_pseudo = self.cfg.noise_pseudo;
        let min_count = self.cfg.min_count;
        let w = self.cfg.backoff_weight;
        let rank_smoothing = self.cfg.rank_smoothing;

        let mut all_partitions: std::collections::HashSet<Partition> =
            template.frag_off_table.keys().copied().collect();
        for &(part, _) in counts.rank.keys() {
            all_partitions.insert(part);
        }

        let ion_lists: FxHashMap<Partition, Vec<IonType>> = template
            .frag_off_table
            .iter()
            .map(|(&part, frags)| {
                let ions: Vec<IonType> = frags.iter().map(|f| f.ion_type).collect();
                (part, ions)
            })
            .collect();

        let global_pool = build_global_pool(&counts.rank, n_slots);

        let seg_collapsed = build_seg_collapsed(&counts.rank, n_slots);

        let charge_neighbor = build_charge_neighbor_pool(&counts.rank, n_slots);

        let mut out: FxHashMap<Partition, FxHashMap<IonType, Vec<f32>>> = FxHashMap::default();

        for &part in &all_partitions {
            let ions = match ion_lists.get(&part) {
                Some(v) if !v.is_empty() => v.clone(),
                _ => continue,
            };

            let mut ion_table: FxHashMap<IonType, Vec<f32>> = FxHashMap::default();

            let seg_key = (part.charge, part.parent_mass.to_bits());
            let seg_parent: Option<&FxHashMap<IonType, Vec<u64>>> = seg_collapsed.get(&seg_key);

            let parent_vec = |ion: IonType, ps: f32| -> Vec<f32> {
                if let Some(pr) = prior {
                    if let Some(dist) = pr.rank_dist_table.get(&part).and_then(|m| m.get(&ion)) {
                        if dist.len() == n_slots {
                            return dist.clone();
                        }
                    }
                }
                if let Some(seg_map) = seg_parent {
                    if let Some(raw) = seg_map.get(&ion) {
                        let n: u64 = raw.iter().sum();
                        if n >= min_count {
                            return normalize_with_pseudo(raw, n_slots, ps);
                        }
                    }
                }
                if charge_bucket(part.charge) == 4 {
                    if let Some(raw) = charge_neighbor.get(&(4, ion)) {
                        let n: u64 = raw.iter().sum();
                        if n >= min_count {
                            return normalize_with_pseudo(raw, n_slots, ps);
                        }
                    }
                }
                let graw = global_pool.get(&ion).map(|v| v.as_slice()).unwrap_or(&[]);
                normalize_with_pseudo(graw, n_slots, ps)
            };

            for &ion in &ions {
                let raw = counts
                    .rank
                    .get(&(part, ion))
                    .map(|v| v.as_slice())
                    .unwrap_or(&[]);
                let n: u64 = raw.iter().sum();
                let emp = normalize_with_pseudo(raw, n_slots, pseudo);
                let blended = if n < min_count {
                    blend(&emp, &parent_vec(ion, pseudo), n as f32, w)
                } else {
                    emp
                };
                let final_dist = if rank_smoothing {
                    smooth_rank_window(&blended, n_slots - 1)
                } else {
                    blended
                };
                ion_table.insert(ion, final_dist);
            }

            for (&(p2, ion), raw_counts) in &counts.rank {
                if p2 != part || !ion.is_loss() || ion_table.contains_key(&ion) {
                    continue;
                }
                let raw = raw_counts.as_slice();
                let n: u64 = raw.iter().sum();
                let emp = normalize_with_pseudo(raw, n_slots, pseudo);
                let blended = if n < min_count {
                    blend(&emp, &parent_vec(ion, pseudo), n as f32, w)
                } else {
                    emp
                };
                let final_dist = if rank_smoothing {
                    smooth_rank_window(&blended, n_slots - 1)
                } else {
                    blended
                };
                ion_table.insert(ion, final_dist);
            }

            let noise_raw = counts
                .rank
                .get(&(part, IonType::Noise))
                .map(|v| v.as_slice())
                .unwrap_or(&[]);
            let noise_n: u64 = noise_raw.iter().sum();
            let noise_emp = normalize_with_pseudo(noise_raw, n_slots, noise_pseudo);
            let noise_dist = if noise_n < min_count {
                blend(
                    &noise_emp,
                    &parent_vec(IonType::Noise, noise_pseudo),
                    noise_n as f32,
                    w,
                )
            } else {
                noise_emp
            };
            ion_table.insert(IonType::Noise, noise_dist);

            out.insert(part, ion_table);
        }

        out
    }

    fn build_error_tables(
        &self,
        counts: &CountStats,
        template: &Param,
        esf: i32,
        prior: Option<&Param>,
    ) -> (
        FxHashMap<Partition, Vec<f32>>,
        FxHashMap<Partition, Vec<f32>>,
    ) {
        if esf <= 0 {
            return (FxHashMap::default(), FxHashMap::default());
        }
        let dist_len = (esf as usize) * 2 + 1;
        let pseudo = self.cfg.pseudo;
        let min_count = self.cfg.min_count;
        let w = self.cfg.backoff_weight;

        let mut global_ion_raw = vec![0u64; dist_len];
        let mut global_noise_raw = vec![0u64; dist_len];
        for v in counts.error.values() {
            for (i, &c) in v.iter().enumerate() {
                if i < dist_len {
                    global_ion_raw[i] = global_ion_raw[i].saturating_add(c);
                }
            }
        }
        for v in counts.noise_error.values() {
            for (i, &c) in v.iter().enumerate() {
                if i < dist_len {
                    global_noise_raw[i] = global_noise_raw[i].saturating_add(c);
                }
            }
        }
        let global_ion_norm = normalize_with_pseudo(&global_ion_raw, dist_len, pseudo);
        let global_noise_norm = normalize_with_pseudo(&global_noise_raw, dist_len, pseudo);

        let mut ion_out: FxHashMap<Partition, Vec<f32>> = FxHashMap::default();
        let mut noise_out: FxHashMap<Partition, Vec<f32>> = FxHashMap::default();

        for &part in &template.partitions {
            let ion_raw = counts.error.get(&part).map(|v| v.as_slice()).unwrap_or(&[]);
            let ion_n: u64 = ion_raw.iter().sum();
            let ion_emp = normalize_with_pseudo(ion_raw, dist_len, pseudo);
            let ion_parent: Vec<f32> = prior
                .and_then(|p| p.ion_err_dist_table.get(&part))
                .filter(|d| d.len() == dist_len)
                .cloned()
                .unwrap_or_else(|| global_ion_norm.clone());
            let ion_dist = if ion_n < min_count {
                blend(&ion_emp, &ion_parent, ion_n as f32, w)
            } else {
                ion_emp
            };
            ion_out.insert(part, ion_dist);

            let noise_raw = counts
                .noise_error
                .get(&part)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);
            let noise_n: u64 = noise_raw.iter().sum();
            let noise_emp = normalize_with_pseudo(noise_raw, dist_len, pseudo);
            let noise_parent: Vec<f32> = prior
                .and_then(|p| p.noise_err_dist_table.get(&part))
                .filter(|d| d.len() == dist_len)
                .cloned()
                .unwrap_or_else(|| global_noise_norm.clone());
            let noise_dist = if noise_n < min_count {
                blend(&noise_emp, &noise_parent, noise_n as f32, w)
            } else {
                noise_emp
            };
            noise_out.insert(part, noise_dist);
        }

        (ion_out, noise_out)
    }

    fn build_existence_table(
        &self,
        counts: &CountStats,
        template: &Param,
        prior: Option<&Param>,
    ) -> FxHashMap<Partition, Vec<f32>> {
        const N_EX: usize = 4;
        let pseudo = self.cfg.pseudo;
        let min_count = self.cfg.min_count;
        let w = self.cfg.backoff_weight;

        let mut global_raw = [0u64; N_EX];
        for (&(_part, idx), &c) in &counts.existence {
            if (idx as usize) < N_EX {
                global_raw[idx as usize] = global_raw[idx as usize].saturating_add(c);
            }
        }
        let global_norm = normalize_with_pseudo(&global_raw, N_EX, pseudo);

        let mut out: FxHashMap<Partition, Vec<f32>> = FxHashMap::default();
        for &part in &template.partitions {
            let raw: Vec<u64> = (0..N_EX as u32)
                .map(|idx| counts.existence.get(&(part, idx)).copied().unwrap_or(0))
                .collect();
            let n: u64 = raw.iter().sum();
            let emp = normalize_with_pseudo(&raw, N_EX, pseudo);
            let parent: Vec<f32> = prior
                .and_then(|p| p.ion_existence_table.get(&part))
                .filter(|d| d.len() == N_EX)
                .cloned()
                .unwrap_or_else(|| global_norm.clone());
            let dist = if n < min_count {
                blend(&emp, &parent, n as f32, w)
            } else {
                emp
            };
            out.insert(part, dist);
        }
        out
    }
}

pub fn smooth_rank_window(dist: &[f32], max_rank: usize) -> Vec<f32> {
    let n = dist.len();
    let last = max_rank.min(n.saturating_sub(1));
    let mut out = dist.to_vec();
    let halfwidth = |r: usize| -> usize {
        if r < 3 {
            0
        } else if r < 5 {
            1
        } else if r < 10 {
            2
        } else if r < 20 {
            3
        } else if r < 50 {
            4
        } else {
            5
        }
    };
    for (i, slot) in out.iter_mut().enumerate().take(last) {
        let hw = halfwidth(i);
        if hw == 0 {
            continue;
        }
        let lo = i.saturating_sub(hw);
        let hi = (i + hw + 1).min(last);
        let mut s = 0.0f32;
        let mut c = 0usize;
        for v in dist.iter().take(hi).skip(lo) {
            s += *v;
            c += 1;
        }
        if c > 0 {
            *slot = s / c as f32;
        }
    }
    let tot: f32 = out.iter().sum();
    if tot > 0.0 {
        for x in &mut out {
            *x /= tot;
        }
    }
    out
}

fn normalize_with_pseudo(raw: &[u64], n_slots: usize, pseudo: f32) -> Vec<f32> {
    let mut v: Vec<f32> = Vec::with_capacity(n_slots);
    for i in 0..n_slots {
        let c = raw.get(i).copied().unwrap_or(0) as f32;
        v.push(c + pseudo);
    }
    let total: f32 = v.iter().sum();
    if total > 0.0 {
        for x in &mut v {
            *x /= total;
        }
    }
    v
}

fn blend(emp: &[f32], prior: &[f32], n: f32, w: f32) -> Vec<f32> {
    let denom = n + w;
    if denom <= 0.0 {
        return prior.to_vec();
    }
    emp.iter()
        .zip(prior.iter())
        .map(|(&e, &p)| (n * e + w * p) / denom)
        .collect()
}

fn build_global_pool(
    rank: &FxHashMap<(Partition, IonType), Vec<u64>>,
    n_slots: usize,
) -> FxHashMap<IonType, Vec<u64>> {
    let mut pool: FxHashMap<IonType, Vec<u64>> = FxHashMap::default();
    for (&(_part, ion), v) in rank {
        let entry = pool.entry(ion).or_insert_with(|| vec![0u64; n_slots]);
        if entry.len() < n_slots {
            entry.resize(n_slots, 0);
        }
        for (i, &c) in v.iter().enumerate() {
            if i < n_slots {
                entry[i] = entry[i].saturating_add(c);
            }
        }
    }
    pool
}

fn build_seg_collapsed(
    rank: &FxHashMap<(Partition, IonType), Vec<u64>>,
    n_slots: usize,
) -> FxHashMap<(i32, u32), FxHashMap<IonType, Vec<u64>>> {
    let mut out: FxHashMap<(i32, u32), FxHashMap<IonType, Vec<u64>>> = FxHashMap::default();
    for (&(part, ion), v) in rank {
        let key = (part.charge, part.parent_mass.to_bits());
        let ion_map = out.entry(key).or_default();
        let entry = ion_map.entry(ion).or_insert_with(|| vec![0u64; n_slots]);
        if entry.len() < n_slots {
            entry.resize(n_slots, 0);
        }
        for (i, &c) in v.iter().enumerate() {
            if i < n_slots {
                entry[i] = entry[i].saturating_add(c);
            }
        }
    }
    out
}

fn charge_bucket(charge: i32) -> i32 {
    match charge {
        c if c <= 2 => 2,
        3 => 3,
        _ => 4,
    }
}

fn build_charge_neighbor_pool(
    rank: &FxHashMap<(Partition, IonType), Vec<u64>>,
    n_slots: usize,
) -> FxHashMap<(i32, IonType), Vec<u64>> {
    let mut out: FxHashMap<(i32, IonType), Vec<u64>> = FxHashMap::default();
    for (&(part, ion), v) in rank {
        let entry = out
            .entry((charge_bucket(part.charge), ion))
            .or_insert_with(|| vec![0u64; n_slots]);
        if entry.len() < n_slots {
            entry.resize(n_slots, 0);
        }
        for (i, &c) in v.iter().enumerate() {
            if i < n_slots {
                entry[i] = entry[i].saturating_add(c);
            }
        }
    }
    out
}

fn build_charge_hist(charge: &FxHashMap<i32, u64>) -> Vec<(i32, i32)> {
    let mut hist: Vec<(i32, i32)> = charge
        .iter()
        .map(|(&ch, &n)| (ch, n.min(i32::MAX as u64) as i32))
        .collect();
    hist.sort_by_key(|(ch, _)| *ch);
    hist
}

fn charge_range(hist: &[(i32, i32)], template: &Param) -> (i32, i32) {
    if hist.is_empty() {
        return (template.min_charge, template.max_charge);
    }
    let min = hist
        .iter()
        .map(|(c, _)| *c)
        .min()
        .unwrap_or(template.min_charge);
    let max = hist
        .iter()
        .map(|(c, _)| *c)
        .max()
        .unwrap_or(template.max_charge);
    (min, max)
}
