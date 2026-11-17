pub fn nnls(a: &[f64], m: usize, n: usize, b: &[f64], max_iter: usize) -> Vec<f64> {
    let mut x = vec![0.0f64; n];
    if m == 0 || n == 0 || a.len() != m * n || b.len() != m {
        return x;
    }
    let scale = b.iter().fold(0.0f64, |acc, v| acc.max(v.abs())).max(1.0);
    let tol = 1e-10 * scale;

    let mut passive = vec![false; n];
    let mut w = vec![0.0f64; n];
    let mut s = vec![0.0f64; n];

    let residual_grad = |x: &[f64], w: &mut [f64]| {
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
        let mut best: Option<(usize, f64)> = None;
        for j in 0..n {
            if !passive[j] && w[j] > tol && best.is_none_or(|(_, bw)| w[j] > bw) {
                best = Some((j, w[j]));
            }
        }
        let Some((j_in, _)) = best else { break };
        passive[j_in] = true;

        let mut inner = 0usize;
        loop {
            inner += 1;
            if inner > 3 * n + 3 {
                break;
            }
            if !solve_passive(a, m, n, b, &passive, &mut s) {
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
    let mut g = vec![0.0f64; k * (k + 1)];
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
