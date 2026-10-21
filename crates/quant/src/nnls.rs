//! Non-negative least squares (Lawson–Hanson active-set algorithm).
//!
//! Solves `min ‖A·x − b‖₂` subject to `x ≥ 0` for the small dense systems the
//! isobaric impurity correction needs (one reporter channel per unknown, at
//! most 18). The same solver OpenMS's `IsobaricIsotopeCorrector` uses
//! (`NonNegativeLeastSquaresSolver`), so corrected channel intensities agree
//! with IsobaricAnalyzer / IsobaricWorkflow for the same matrix.

/// Solve `min ‖A·x − b‖` with `x ≥ 0`. `a` is row-major `m × n`.
///
/// Returns `x` (length `n`). A degenerate system (empty, mismatched sizes)
/// returns all zeros. `max_iter` bounds the outer iterations (3·n is plenty for
/// well-conditioned correction matrices; the caller passes a small multiple).
pub fn nnls(a: &[f64], m: usize, n: usize, b: &[f64], max_iter: usize) -> Vec<f64> {
    let mut x = vec![0.0f64; n];
    if m == 0 || n == 0 || a.len() != m * n || b.len() != m {
        return x;
    }
    // Tolerance on the dual variable: relative to the problem scale so a
    // channel vector in the 1e6 range and one in the 1e2 range converge alike.
    let scale = b.iter().fold(0.0f64, |acc, v| acc.max(v.abs())).max(1.0);
    let tol = 1e-10 * scale;

    let mut passive = vec![false; n]; // P: free (positive) set
    let mut w = vec![0.0f64; n];
    let mut s = vec![0.0f64; n];

    let residual_grad = |x: &[f64], w: &mut [f64]| {
        // w = Aᵀ (b − A x)
        let mut r = vec![0.0f64; m];
        for i in 0..m {
            let mut acc = b[i];
            for j in 0..n {
                acc -= a[i * n + j] * x[j];
            }
            r[i] = acc;
        }
        for j in 0..n {
            let mut acc = 0.0;
            for i in 0..m {
                acc += a[i * n + j] * r[i];
            }
            w[j] = acc;
        }
    };

    residual_grad(&x, &mut w);
    let mut outer = 0usize;
    loop {
        outer += 1;
        if outer > max_iter {
            break;
        }
        // Pick the most positive gradient among the active (zero) variables.
        let mut best: Option<(usize, f64)> = None;
        for j in 0..n {
            if !passive[j] && w[j] > tol && best.is_none_or(|(_, bw)| w[j] > bw) {
                best = Some((j, w[j]));
            }
        }
        let Some((j_in, _)) = best else { break };
        passive[j_in] = true;

        // Inner loop: unconstrained solve on the passive set, then step back
        // toward the feasible region if any passive variable went negative.
        let mut inner = 0usize;
        loop {
            inner += 1;
            if inner > 3 * n + 3 {
                break;
            }
            if !solve_passive(a, m, n, b, &passive, &mut s) {
                // Singular subsystem: drop the variable just added and stop.
                passive[j_in] = false;
                break;
            }
            let all_positive = (0..n).all(|j| !passive[j] || s[j] > 0.0);
            if all_positive {
                for j in 0..n {
                    x[j] = if passive[j] { s[j] } else { 0.0 };
                }
                break;
            }
            // alpha = min over passive j with s_j <= 0 of x_j / (x_j − s_j)
            let mut alpha = f64::INFINITY;
            for j in 0..n {
                if passive[j] && s[j] <= 0.0 {
                    let denom = x[j] - s[j];
                    let a_j = if denom.abs() > 0.0 { x[j] / denom } else { 0.0 };
                    if a_j < alpha {
                        alpha = a_j;
                    }
                }
            }
            if !alpha.is_finite() {
                alpha = 0.0;
            }
            for j in 0..n {
                if passive[j] {
                    x[j] += alpha * (s[j] - x[j]);
                    if x[j] <= tol {
                        x[j] = 0.0;
                        passive[j] = false;
                    }
                }
            }
        }
        residual_grad(&x, &mut w);
    }
    x
}

/// Unconstrained least squares restricted to the passive columns, by the
/// normal equations with partial pivoting. Returns false if singular.
fn solve_passive(
    a: &[f64],
    m: usize,
    n: usize,
    b: &[f64],
    passive: &[bool],
    s: &mut [f64],
) -> bool {
    let cols: Vec<usize> = (0..n).filter(|&j| passive[j]).collect();
    let k = cols.len();
    for v in s.iter_mut() {
        *v = 0.0;
    }
    if k == 0 {
        return true;
    }
    // Normal equations: (AᵀA) y = Aᵀb over the passive columns.
    let mut g = vec![0.0f64; k * (k + 1)]; // augmented [AᵀA | Aᵀb]
    for (ci, &cj) in cols.iter().enumerate() {
        for (ri, &rj) in cols.iter().enumerate() {
            let mut acc = 0.0;
            for i in 0..m {
                acc += a[i * n + rj] * a[i * n + cj];
            }
            g[ri * (k + 1) + ci] = acc;
        }
        let mut acc = 0.0;
        for i in 0..m {
            acc += a[i * n + cj] * b[i];
        }
        g[ci * (k + 1) + k] = acc;
    }
    // Gaussian elimination with partial pivoting.
    for p in 0..k {
        let mut piv = p;
        for r in (p + 1)..k {
            if g[r * (k + 1) + p].abs() > g[piv * (k + 1) + p].abs() {
                piv = r;
            }
        }
        if g[piv * (k + 1) + p].abs() < 1e-14 {
            return false;
        }
        if piv != p {
            for c in 0..=k {
                g.swap(p * (k + 1) + c, piv * (k + 1) + c);
            }
        }
        for r in (p + 1)..k {
            let f = g[r * (k + 1) + p] / g[p * (k + 1) + p];
            if f != 0.0 {
                for c in p..=k {
                    g[r * (k + 1) + c] -= f * g[p * (k + 1) + c];
                }
            }
        }
    }
    for p in (0..k).rev() {
        let mut acc = g[p * (k + 1) + k];
        for c in (p + 1)..k {
            acc -= g[p * (k + 1) + c] * s[cols[c]];
        }
        s[cols[p]] = acc / g[p * (k + 1) + p];
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mat_vec(a: &[f64], m: usize, n: usize, x: &[f64]) -> Vec<f64> {
        (0..m)
            .map(|i| (0..n).map(|j| a[i * n + j] * x[j]).sum())
            .collect()
    }

    #[test]
    fn identity_returns_input() {
        let a = [1.0, 0.0, 0.0, 1.0];
        let b = [3.0, 4.0];
        let x = nnls(&a, 2, 2, &b, 20);
        assert!((x[0] - 3.0).abs() < 1e-9 && (x[1] - 4.0).abs() < 1e-9);
    }

    #[test]
    fn recovers_true_intensities_through_impurity_matrix() {
        // 3 channels, 5% spill-over to the next channel.
        let a = [0.95, 0.0, 0.0, 0.05, 0.95, 0.0, 0.0, 0.05, 1.0];
        let truth = [1000.0, 500.0, 2000.0];
        let b = mat_vec(&a, 3, 3, &truth);
        let x = nnls(&a, 3, 3, &b, 20);
        for (xi, ti) in x.iter().zip(truth.iter()) {
            assert!((xi - ti).abs() < 1e-6, "{x:?} vs {truth:?}");
        }
    }

    #[test]
    fn clamps_negative_naive_solution_to_zero() {
        // A naive solve would need a negative x[1] to reproduce b exactly.
        let a = [1.0, 0.5, 0.0, 1.0];
        let b = [1.0, -0.4];
        let x = nnls(&a, 2, 2, &b, 20);
        assert!(x.iter().all(|&v| v >= 0.0), "{x:?}");
        assert!((x[0] - 1.0).abs() < 1e-9 && x[1] == 0.0, "{x:?}");
    }

    #[test]
    fn zero_vector_gives_zero_solution() {
        let a = [0.9, 0.1, 0.1, 0.9];
        let x = nnls(&a, 2, 2, &[0.0, 0.0], 20);
        assert_eq!(x, vec![0.0, 0.0]);
    }
}
