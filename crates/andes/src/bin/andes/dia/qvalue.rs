//! Precursor-level q-values for DIA: target-decoy competition over each precursor's best row, overall and per class.

use super::*;

/// Reports precursor-level q-values for DIA: target-decoy competition over the best
/// Percolator score of each precursor (run, peptide form, charge), with q = (D + 1) / T
/// and monotone in the score. A row receives the q of the threshold at its own score, so
/// keeping rows with q <= x keeps the precursors whose best row passes, at precursor FDR x.
/// Rows Percolator did not report keep no entry. Returns the precursor counts (targets,
/// decoys, target precursors at q <= 0.01) and the best `(score, is decoy)` of every precursor,
/// for comparing the run's precursor-level target-decoy curve at other q-values. A precursor
/// PIN (one with the `score::PRESENT_COLUMN` column) is keyed by its Peptide string and
/// charge columns and needs no spectra. With `interference` (--dia-interference) and a PSM
/// PIN, [`interference_removed`] drops precursors from the competition first; their target rows
/// get q = 1.
pub(crate) fn precursor_q(
    pin: &Path,
    spectra: &[Spectrum],
    targets: &mut HashMap<String, output::PercolatorPsm>,
    decoys: &HashMap<String, output::PercolatorPsm>,
    interference: Option<Interference>,
) -> Result<PrecursorQ, String> {
    let header = {
        use std::io::BufRead;
        let f = std::fs::File::open(pin).map_err(|e| format!("{}: {e}", pin.display()))?;
        std::io::BufReader::new(f)
            .lines()
            .find(|l| !matches!(l, Ok(l) if l.is_empty()))
            .transpose()
            .map_err(|e| format!("{}: {e}", pin.display()))?
            .unwrap_or_default()
    };
    let h: Vec<&str> = header.split('\t').collect();
    let at = |n: &str| {
        h.iter()
            .position(|c| *c == n)
            .ok_or_else(|| format!("{}: PIN has no {n} column", pin.display()))
    };
    let (c_id, c_pep, c_label) = (at("SpecId")?, at("Peptide")?, at("Label")?);
    let c_file = h.iter().position(|c| *c == "FileName");
    let charge_cols: Vec<(usize, u8)> = h
        .iter()
        .enumerate()
        .filter_map(|(i, c)| Some((i, c.strip_prefix("charge")?.parse().ok()?)))
        .collect();
    let table = if h.contains(&crate::dia::score::PRESENT_COLUMN) {
        None
    } else {
        Some(read_pin(pin, spectra, &[])?)
    };
    let mut keys: HashMap<(String, String, u8), u32> = HashMap::new();
    let mut rows: Vec<(u32, f64, bool, String)> = Vec::new();
    let mut raw: Vec<(String, String, u8)> = Vec::new();
    let mut line_of: Vec<usize> = Vec::new();
    type Parsed<'a> = Option<(bool, String, u8, &'a str, f64, &'a str, &'a str, usize)>;
    let targets_ro: &HashMap<String, output::PercolatorPsm> = targets;
    crate::dia::blocks::read_blocks(
        pin,
        |_| Ok(()),
        |first, lines| {
            let parsed: Vec<Result<Parsed<'_>, String>> = lines
                .par_iter()
                .enumerate()
                .map(|(j, line)| {
                    let i = first + j;
                    let x: Vec<&str> = line.split('\t').collect();
                    if x.len() <= c_pep.max(c_label).max(c_id) {
                        return Err(format!("{}: short PIN row {}", pin.display(), i + 2));
                    }
                    let (decoy, peptide, charge) = match &table {
                        Some(t) => {
                            let r = t.rows.get(i).ok_or_else(|| {
                                format!("{}: row {} not read", pin.display(), i + 2)
                            })?;
                            let peptide = if r.form == NONE {
                                x[c_pep].to_string()
                            } else {
                                format!("#{}", r.form)
                            };
                            (r.decoy, peptide, r.charge)
                        }
                        None => (
                            x[c_label].trim() == "-1",
                            x[c_pep].to_string(),
                            charge_cols
                                .iter()
                                .find(|(c, _)| x.get(*c) == Some(&"1"))
                                .map_or(0, |c| c.1),
                        ),
                    };
                    let id = x[c_id];
                    let psm = if decoy {
                        decoys.get(id)
                    } else {
                        targets_ro.get(id)
                    };
                    let Some(score) = psm.and_then(|p| p.score) else {
                        return Ok(None);
                    };
                    Ok(Some((
                        decoy,
                        peptide,
                        charge,
                        id,
                        score,
                        c_file.map_or("", |i| x[i]),
                        x[c_pep],
                        i,
                    )))
                })
                .collect();
            for p in parsed {
                let Some((decoy, peptide, charge, id, score, run, raw_pep, i)) = p? else {
                    continue;
                };
                let n = keys.len() as u32;
                let k = *keys.entry((run.to_string(), peptide, charge)).or_insert(n);
                raw.push((run.to_string(), raw_pep.to_string(), charge));
                rows.push((k, score, decoy, id.to_string()));
                line_of.push(i);
            }
            Ok(())
        },
    )?;
    let removed = match (&table, interference) {
        (Some(t), Some(choice)) => {
            let q0 = precursor_q_values(keys.len(), &rows);
            removed_keys(&keys, &rows, &line_of, &q0, t, spectra, choice)
        }
        _ => vec![false; keys.len()],
    };
    for (k, _, decoy, id) in &rows {
        if removed[*k as usize] && !*decoy {
            if let Some(p) = targets.get_mut(id) {
                p.q_value = 1.0;
            }
        }
    }
    let removed_precursors = rows
        .iter()
        .zip(raw)
        .filter(|(r, _)| removed[r.0 as usize])
        .map(|(_, n)| n)
        .collect();
    rows.retain(|r| !removed[r.0 as usize]);
    let q = precursor_q_values(keys.len(), &rows);
    let n_targets = q.iter().filter(|p| p.is_some_and(|p| !p.1)).count();
    let n_decoys = q.iter().filter(|p| p.is_some_and(|p| p.1)).count();
    let passing = q.iter().flatten().filter(|p| !p.1 && p.2 <= 0.01).count();
    let curve = q_curve(&q);
    for (_, score, decoy, id) in &rows {
        if *decoy {
            continue;
        }
        if let Some(p) = targets.get_mut(id) {
            p.q_value = curve_at(&curve, *score);
        }
    }
    let best = q.iter().flatten().map(|p| (p.0, p.1)).collect();
    Ok(((n_targets, n_decoys, passing), best, removed_precursors))
}

/// Best score of each precursor and its precursor-level q: `(best score, is decoy, q)`,
/// `None` for a key no row reached. Ties in score count decoys first.
pub(super) fn precursor_q_values(
    n_keys: usize,
    rows: &[(u32, f64, bool, String)],
) -> Vec<Option<(f64, bool, f64)>> {
    let mut best: Vec<Option<(f64, bool)>> = vec![None; n_keys];
    for (k, score, decoy, _) in rows {
        let b = &mut best[*k as usize];
        if b.is_none_or(|(s, _)| *score > s) {
            *b = Some((*score, *decoy));
        }
    }
    let mut order: Vec<usize> = (0..n_keys).filter(|&i| best[i].is_some()).collect();
    order.sort_by(|&a, &b| {
        let (sa, da) = best[a].unwrap_or((f64::NEG_INFINITY, true));
        let (sb, db) = best[b].unwrap_or((f64::NEG_INFINITY, true));
        sb.total_cmp(&sa).then(db.cmp(&da))
    });
    let mut q = vec![1.0f64; n_keys];
    let (mut t, mut d) = (0usize, 0usize);
    let mut i = 0;
    while i < order.len() {
        let s = best[order[i]].map_or(f64::NEG_INFINITY, |b| b.0);
        let mut j = i;
        while j < order.len() && best[order[j]].map_or(f64::NEG_INFINITY, |b| b.0) == s {
            if best[order[j]].is_some_and(|b| b.1) {
                d += 1;
            } else {
                t += 1;
            }
            j += 1;
        }
        let fdr = ((d + 1) as f64 / t.max(1) as f64).min(1.0);
        for &k in &order[i..j] {
            q[k] = fdr;
        }
        i = j;
    }
    let mut running = 1.0f64;
    for &k in order.iter().rev() {
        running = running.min(q[k]);
        q[k] = running;
    }
    best.iter()
        .zip(q)
        .map(|(b, qv)| b.map(|(s, d)| (s, d, qv)))
        .collect()
}

/// Score thresholds in descending order with the precursor q at each.
pub(super) fn q_curve(q: &[Option<(f64, bool, f64)>]) -> Vec<(f64, f64)> {
    let mut c: Vec<(f64, f64)> = q.iter().flatten().map(|p| (p.0, p.2)).collect();
    c.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.total_cmp(&b.1)));
    c
}

/// The precursor q a row with this score gets: the q of the lowest precursor best score
/// that is still at or above it; 1 below every precursor.
pub(super) fn curve_at(curve: &[(f64, f64)], score: f64) -> f64 {
    let n = curve.partition_point(|p| p.0 >= score);
    if n == 0 {
        1.0
    } else {
        curve[n - 1].1
    }
}

/// Precursor-level q of every row computed within its class only: one target-decoy competition
/// per class (rows of class `false` against each other, rows of class `true` against each
/// other), with the rules of [`precursor_q`]. Returns the q of each row, 1 for a row whose
/// class holds no precursor at or above its score.
pub(crate) fn classwise_q_values(
    n_keys: usize,
    rows: &[(u32, f64, bool, String)],
    class: &[bool],
) -> Vec<f64> {
    let mut q = vec![1.0; rows.len()];
    for c in [false, true] {
        let idx: Vec<usize> = (0..rows.len()).filter(|&i| class[i] == c).collect();
        let sub: Vec<(u32, f64, bool, String)> = idx.iter().map(|&i| rows[i].clone()).collect();
        let curve = q_curve(&precursor_q_values(n_keys, &sub));
        for &i in &idx {
            q[i] = curve_at(&curve, rows[i].1);
        }
    }
    q
}

/// Replaces the targets' q-values of a precursor PIN by the class-wise q
/// ([`classwise_q_values`], classes by `class_col` 0 or 1) and counts the target precursors at
/// q <= 0.01 by class under the pooled q they had and under the class-wise q. Returns
/// `[[pooled class 0, pooled class 1], [class-wise class 0, class-wise class 1]]`.
pub(crate) fn precursor_q_by_class(
    pin: &Path,
    class_col: &str,
    targets: &mut HashMap<String, output::PercolatorPsm>,
    decoys: &HashMap<String, output::PercolatorPsm>,
) -> Result<[[usize; 2]; 2], String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", pin.display());
    let f = std::fs::File::open(pin).map_err(|e| err(&e))?;
    let mut lines = std::io::BufRead::lines(std::io::BufReader::new(f));
    let header = lines
        .next()
        .ok_or_else(|| err(&"empty PIN"))?
        .map_err(|e| err(&e))?;
    let h: Vec<&str> = header.split('\t').collect();
    let at = |n: &str| {
        h.iter()
            .position(|c| *c == n)
            .ok_or_else(|| format!("{}: PIN has no {n} column", pin.display()))
    };
    let (c_id, c_pep, c_label, c_class) = (at("SpecId")?, at("Peptide")?, at("Label")?, at(class_col)?);
    let charge_cols: Vec<(usize, u8)> = h
        .iter()
        .enumerate()
        .filter_map(|(i, c)| Some((i, c.strip_prefix("charge")?.parse().ok()?)))
        .collect();
    let mut keys: HashMap<(String, u8), u32> = HashMap::new();
    let mut rows: Vec<(u32, f64, bool, String)> = Vec::new();
    let mut class: Vec<bool> = Vec::new();
    for line in lines {
        let line = line.map_err(|e| err(&e))?;
        if line.is_empty() {
            continue;
        }
        let x: Vec<&str> = line.split('\t').collect();
        if x.len() <= c_pep.max(c_label).max(c_id).max(c_class) {
            return Err(err(&"short PIN row"));
        }
        let decoy = x[c_label].trim() == "-1";
        let psm = if decoy { decoys.get(x[c_id]) } else { targets.get(x[c_id]) };
        let Some(score) = psm.and_then(|p| p.score) else {
            continue;
        };
        let charge = charge_cols
            .iter()
            .find(|(c, _)| x.get(*c) == Some(&"1"))
            .map_or(0, |c| c.1);
        let n = keys.len() as u32;
        let k = *keys.entry((x[c_pep].to_string(), charge)).or_insert(n);
        rows.push((k, score, decoy, x[c_id].to_string()));
        class.push(x[c_class].trim().parse::<f32>().unwrap_or(0.0) > 0.5);
    }
    let q = classwise_q_values(keys.len(), &rows, &class);
    let mut counts = [[0usize; 2]; 2];
    for (i, (_, _, decoy, id)) in rows.iter().enumerate() {
        if *decoy {
            continue;
        }
        let c = usize::from(class[i]);
        if let Some(p) = targets.get_mut(id) {
            counts[0][c] += usize::from(p.q_value <= 0.01);
            counts[1][c] += usize::from(q[i] <= 0.01);
            p.q_value = q[i];
        }
    }
    Ok(counts)
}
