//! Ranked target-decoy q-values.
//!
//! Callers rank their rows with their own ordering and tie-break, then hand the ranked
//! rows to [`ranked_q`], which walks them once keeping a running state, evaluates a value
//! at each row or each tied group, and applies the suffix (cumulative) minimum. Every
//! behaviour that can change a bit of the result is an explicit argument: when the value
//! is evaluated ([`Eval`]), how ties are detected, the value rule, and the starting cap
//! of the suffix minimum.

/// Running target and decoy counts, the state used by the standard value rules.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Counts {
    pub targets: u64,
    pub decoys: u64,
}

impl Counts {
    pub fn add(&mut self, is_decoy: bool) {
        if is_decoy {
            self.decoys += 1;
        } else {
            self.targets += 1;
        }
    }

    /// `D / max(T, 1)`.
    pub fn ratio(&self) -> f64 {
        self.decoys as f64 / self.targets.max(1) as f64
    }

    /// `(D + 1) / max(T, 1)`.
    pub fn ratio_plus_one(&self) -> f64 {
        (self.decoys + 1) as f64 / self.targets.max(1) as f64
    }

    /// `((D + 1) / k) / max(T, 1)` for `k` decoys per target.
    pub fn ratio_plus_one_over(&self, k: f64) -> f64 {
        ((self.decoys + 1) as f64 / k) / self.targets.max(1) as f64
    }
}

/// When the value is evaluated.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Eval {
    /// After every row; tied rows are not grouped.
    PerRow,
    /// After every row, then every tied group takes its largest value after the suffix
    /// minimum. Within-group order of targets and decoys changes the result.
    PerRowGroupMax,
    /// Once per tied group, after all its rows are counted.
    GroupEnd,
}

/// The suffix minimum: every entry becomes the smallest of itself and all later entries,
/// starting from `cap`.
pub fn suffix_min(values: &mut [f64], cap: f64) {
    let mut running = cap;
    for v in values.iter_mut().rev() {
        if *v < running {
            running = *v;
        }
        *v = running;
    }
}

/// Q-values in rank order.
///
/// `add` folds row `rank` into the running state; `value` turns the state into a value;
/// `same_group(start, rank)` tells whether the row at `rank` ties with the first row of
/// its group at `start` (it is never called for `Eval::PerRow`); `cap` starts the suffix
/// minimum.
pub fn ranked_q<S: Default>(
    n: usize,
    mut add: impl FnMut(&mut S, usize),
    same_group: impl Fn(usize, usize) -> bool,
    eval: Eval,
    value: impl Fn(&S) -> f64,
    cap: f64,
) -> Vec<f64> {
    let mut q = vec![0.0_f64; n];
    let mut state = S::default();
    match eval {
        Eval::PerRow | Eval::PerRowGroupMax => {
            for (rank, qi) in q.iter_mut().enumerate() {
                add(&mut state, rank);
                *qi = value(&state);
            }
            suffix_min(&mut q, cap);
            if eval == Eval::PerRowGroupMax {
                let mut start = 0usize;
                while start < n {
                    let mut end = start + 1;
                    while end < n && same_group(start, end) {
                        end += 1;
                    }
                    let mut worst = q[start];
                    for &qi in &q[start + 1..end] {
                        if qi > worst {
                            worst = qi;
                        }
                    }
                    q[start..end].fill(worst);
                    start = end;
                }
            }
        }
        Eval::GroupEnd => {
            let mut start = 0usize;
            while start < n {
                add(&mut state, start);
                let mut end = start + 1;
                while end < n && same_group(start, end) {
                    add(&mut state, end);
                    end += 1;
                }
                q[start..end].fill(value(&state));
                start = end;
            }
            suffix_min(&mut q, cap);
        }
    }
    q
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(rows: &[(f32, bool)], eval: Eval, value: fn(&Counts) -> f64, cap: f64) -> Vec<f64> {
        ranked_q::<Counts>(
            rows.len(),
            |c, i| c.add(rows[i].1),
            |a, b| rows[a].0.to_bits() == rows[b].0.to_bits(),
            eval,
            value,
            cap,
        )
    }

    #[test]
    fn decoy_ratio_is_the_suffix_minimum() {
        let r = [
            (5.0_f32, false),
            (4.0, false),
            (3.0, true),
            (2.0, false),
            (1.0, true),
        ];
        let q = run(&r, Eval::PerRow, Counts::ratio, 1.0);
        assert_eq!(q, vec![0.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0]);
    }

    #[test]
    fn plus_one_value_is_capped_at_one() {
        let r = [(3.0_f32, true), (2.0, false), (1.0, false)];
        let q = run(&r, Eval::PerRow, Counts::ratio_plus_one, 1.0);
        assert_eq!(q, vec![1.0, 1.0, 1.0]);
        let r = [(3.0_f32, false), (2.0, false), (1.0, false)];
        let q = run(&r, Eval::PerRow, Counts::ratio_plus_one, 1.0);
        assert_eq!(q, vec![1.0 / 3.0; 3]);
    }

    #[test]
    fn an_infinite_cap_leaves_values_above_one() {
        let r = [(3.0_f32, true), (2.0, false)];
        let q = run(&r, Eval::PerRow, Counts::ratio_plus_one, f64::INFINITY);
        assert_eq!(q, vec![2.0, 2.0]);
    }

    #[test]
    fn per_row_evaluation_depends_on_order_within_a_tie() {
        let a = [(2.0_f32, false), (1.0, false), (1.0, true)];
        let b = [(2.0_f32, false), (1.0, true), (1.0, false)];
        assert_ne!(
            run(&a, Eval::PerRow, Counts::ratio, 1.0),
            run(&b, Eval::PerRow, Counts::ratio, 1.0)
        );
    }

    #[test]
    fn group_maximum_makes_a_tie_independent_of_order() {
        let a = [(2.0_f32, false), (1.0, false), (1.0, true)];
        let b = [(2.0_f32, false), (1.0, true), (1.0, false)];
        let qa = run(&a, Eval::PerRowGroupMax, Counts::ratio, 1.0);
        assert_eq!(qa, run(&b, Eval::PerRowGroupMax, Counts::ratio, 1.0));
        assert_eq!(qa, vec![0.0, 0.5, 0.5]);
    }

    #[test]
    fn group_end_counts_the_whole_tie_before_evaluating() {
        let r = [(2.0_f32, false), (1.0, true), (1.0, false)];
        let q = run(&r, Eval::GroupEnd, Counts::ratio, 1.0);
        assert_eq!(q, vec![0.0, 0.5, 0.5]);
    }

    #[test]
    fn group_end_differs_from_per_row_group_maximum_when_a_later_group_dips() {
        let r = [
            (3.0_f32, false),
            (2.0, true),
            (1.0, false),
            (1.0, false),
            (1.0, true),
        ];
        let per_row = run(&r, Eval::PerRowGroupMax, Counts::ratio, 1.0);
        let end = run(&r, Eval::GroupEnd, Counts::ratio, 1.0);
        assert_ne!(per_row, end);
        assert_eq!(end, vec![0.0, 2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0]);
    }

    #[test]
    fn scaled_plus_one_divides_the_decoys_by_k() {
        let c = Counts {
            targets: 4,
            decoys: 3,
        };
        assert_eq!(c.ratio_plus_one_over(1.0), c.ratio_plus_one());
        assert_eq!(c.ratio_plus_one_over(4.0), 0.25);
    }

    #[test]
    fn a_custom_state_can_weight_the_decoys() {
        let rows = [(false, 0.0_f64), (true, 0.5), (false, 0.0)];
        let q = ranked_q::<(usize, f64)>(
            rows.len(),
            |(t, d), i| {
                if rows[i].0 {
                    *d += rows[i].1;
                } else {
                    *t += 1;
                }
            },
            |_, _| false,
            Eval::PerRow,
            |&(t, d)| d / t.max(1) as f64,
            1.0,
        );
        assert_eq!(q, vec![0.0, 0.25, 0.25]);
    }

    #[test]
    fn empty_input_gives_no_values() {
        let q = ranked_q::<Counts>(0, |_, _| {}, |_, _| true, Eval::GroupEnd, Counts::ratio, 1.0);
        assert!(q.is_empty());
    }
}
