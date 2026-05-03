//! Row-level diff with per-field relative tolerances.

use crate::PinFile;
use std::collections::HashMap;

/// Map of column name to relative tolerance threshold.
pub type Tolerance = HashMap<String, f64>;

/// Result of a tolerance compare. `Ok(())` means all rows pass; `Err` carries
/// a human-readable message with the first batch of failures.
pub fn compare_with_tolerance(
    a: &PinFile,
    b: &PinFile,
    tolerance: &Tolerance,
) -> Result<(), String> {
    let join_keys = ["SpecID", "ScanNum", "Charge", "Peptide"];
    let key_idx: Vec<usize> = join_keys
        .iter()
        .filter_map(|k| a.columns.iter().position(|c| c == *k))
        .collect();
    if key_idx.is_empty() {
        return Err("no join-key columns (SpecID/ScanNum/Charge/Peptide) in header".into());
    }

    let mut a_by_key: HashMap<String, Vec<&str>> = HashMap::new();
    for line in a.body.lines().filter(|l| !l.is_empty()) {
        let cells: Vec<&str> = line.split('\t').collect();
        let key = key_idx
            .iter()
            .map(|i| cells.get(*i).copied().unwrap_or(""))
            .collect::<Vec<_>>()
            .join("|");
        a_by_key.insert(key, cells);
    }

    let mut failures: Vec<String> = Vec::new();
    let mut seen_in_b: usize = 0;
    for line in b.body.lines().filter(|l| !l.is_empty()) {
        let cells_b: Vec<&str> = line.split('\t').collect();
        let key = key_idx
            .iter()
            .map(|i| cells_b.get(*i).copied().unwrap_or(""))
            .collect::<Vec<_>>()
            .join("|");
        seen_in_b += 1;
        let cells_a = match a_by_key.remove(&key) {
            Some(v) => v,
            None => {
                failures.push(format!("only in B: {key}"));
                continue;
            }
        };
        for (col_name, tol) in tolerance {
            let col_idx = a.columns.iter().position(|c| c == col_name);
            let Some(ci) = col_idx else {
                continue;
            };
            let va: f64 = cells_a
                .get(ci)
                .copied()
                .unwrap_or("0")
                .parse()
                .unwrap_or(0.0);
            let vb: f64 = cells_b
                .get(ci)
                .copied()
                .unwrap_or("0")
                .parse()
                .unwrap_or(0.0);
            if va == 0.0 && vb == 0.0 {
                continue;
            }
            let denom = va.abs().max(vb.abs());
            let drift = (va - vb).abs() / denom;
            if drift > *tol {
                failures.push(format!(
                    "{col_name} drift {drift:.6e} exceeds tolerance {tol:.6e} at key={key}"
                ));
            }
        }
    }
    for (k, _) in a_by_key {
        failures.push(format!("only in A: {k}"));
    }
    if failures.is_empty() {
        Ok(())
    } else {
        let preview = failures
            .iter()
            .take(10)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n  ");
        Err(format!(
            "{} row failures (showing first 10) over {} compared rows:\n  {}",
            failures.len(),
            seen_in_b,
            preview
        ))
    }
}
