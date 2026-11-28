//! Entrapment estimates of the false discovery proportion (FDP) of a search's identifications,
//! after Wen et al., Nat. Methods 2025 [cite:wen2025_entrapment].
//!
//! A search run with `--add-entrapment` holds one shuffled entrapment protein per target
//! (`search::decoy::entrapment`). An identification is an entrapment discovery when every
//! protein it maps to is an entrapment protein and its peptide (I and L folded) is not a
//! substring of any original target protein, independently of the search digest. With `N_T`
//! target and `N_E` entrapment discoveries above a cut, the estimators are
//! the lower bound `N_E / (N_T + N_E)`, the combined estimate `N_E (1 + 1/r) / (N_T + N_E)` with
//! `r` the entrapment to target ratio, and the paired estimate
//! `(N_E + N_{E>=s>T} + 2 N_{E>T>=s}) / (N_T + N_E)`, which counts each entrapment discovery
//! against the target it was shuffled from: once more when that target is below the cut and
//! twice more when it is above the cut with a lower score. A known but unobserved target is
//! below the cut. A missing pair identity is different and makes the paired estimate unavailable.
//! The lower bound can show a failure to control the FDR. Upper-bound interpretations of the
//! combined and paired estimates require the method's assumptions, including exchangeability
//! of incorrect target and entrapment matches, an appropriate database ratio and, for pairing,
//! genuine one-to-one source pairs. Finite estimates alone do not establish those assumptions
//! or prove FDR control for an arbitrary acceptance rule. The levels are PSMs, precursors
//! (best PSM per peptidoform and charge), peptides (best PSM per bare sequence) and protein groups.

use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;

use arrow::array::{ArrayRef, BooleanArray, Float64Array, StringArray, UInt64Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;

use crate::versioned::{write_table, TableVersion};

const CURVE_VERSION: &str = "1.0";
const CURVE: TableVersion = TableVersion {
    key: "andes.entrapment.schema",
    version: CURVE_VERSION,
    major: 1,
    what: "entrapment FDP table",
};

/// One identification of the level being estimated.
#[derive(Debug, Clone, PartialEq)]
pub struct Discovery {
    /// The key of the discovery at its level: the PSM id, the bare sequence, the precursor or the
    /// group lead.
    pub key: String,
    /// The pairing key: for a target, its own; for an entrapment discovery, that of the target it
    /// was shuffled from, even when that target has no observed discovery. `None` means the
    /// pairing identity is unknown or ambiguous, not that its target is merely unobserved.
    pub pair: Option<String>,
    pub entrapment: bool,
    /// The rescorer's score, larger is better.
    pub score: f64,
    pub q_value: f64,
}

/// The estimates over the discoveries down to one rank.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Estimate {
    pub targets: usize,
    pub entrapment: usize,
    /// Entrapment discoveries strictly outscoring their paired target, which is not accepted.
    pub pair_below_cut: usize,
    /// Entrapment discoveries strictly outscoring their accepted paired target.
    pub pair_lower: usize,
    pub lower_bound: f64,
    pub combined: f64,
    /// The paired estimate, or NaN when pairing is missing, ambiguous or not one-to-one.
    pub paired: f64,
    ratio: f64,
}

impl Estimate {
    fn of(targets: usize, entrapment: usize, below_cut: usize, lower: usize, r: f64) -> Self {
        let n = (targets + entrapment).max(1) as f64;
        let e = entrapment as f64;
        Estimate {
            targets,
            entrapment,
            pair_below_cut: below_cut,
            pair_lower: lower,
            lower_bound: e / n,
            combined: if r > 0.0 {
                e * (1.0 + 1.0 / r) / n
            } else {
                f64::NAN
            },
            paired: (e + below_cut as f64 + 2.0 * lower as f64) / n,
            ratio: r,
        }
    }
}

/// Sorts `rows` best first (score descending, then q-value, then key) and returns the estimates
/// down to each rank, with `r` the entrapment to target ratio. Paired estimates are NaN unless
/// all rows have finite scores and unambiguous one-to-one pair keys. A known pair key whose
/// target has no observed discovery is treated as below the cut, not as missing pairing.
/// Callers must supply the true source pairing and ensure that omitted targets are unobserved,
/// rather than removed by an unrelated filter.
pub fn curve(rows: &mut [Discovery], r: f64) -> Vec<Estimate> {
    rows.sort_by(|a, b| {
        b.score
            .total_cmp(&a.score)
            .then(a.q_value.total_cmp(&b.q_value))
            .then_with(|| a.key.cmp(&b.key))
            .then(a.entrapment.cmp(&b.entrapment))
            .then_with(|| a.pair.cmp(&b.pair))
    });
    let (target_scores, applicable) = pairing(rows);
    let mut seen_target: HashMap<&str, ()> = HashMap::new();
    let mut waiting: HashMap<&str, usize> = HashMap::new();
    let (mut t, mut e, mut below_cut, mut lower) = (0, 0, 0, 0);
    let mut out = Vec::with_capacity(rows.len());
    for d in rows.iter() {
        let pair = d.pair.as_deref();
        if d.entrapment {
            e += 1;
            match pair {
                Some(p) if seen_target.contains_key(p) => {}
                Some(p) if target_scores.get(p).is_some_and(|score| d.score <= *score) => {}
                Some(p) => {
                    below_cut += 1;
                    *waiting.entry(p).or_default() += 1;
                }
                None => below_cut += 1,
            }
        } else {
            t += 1;
            if let Some(p) = pair {
                if seen_target.insert(p, ()).is_none() {
                    if let Some(k) = waiting.remove(p) {
                        below_cut -= k;
                        lower += k;
                    }
                }
            }
        }
        let mut estimate = Estimate::of(t, e, below_cut, lower, r);
        if !applicable {
            estimate.paired = f64::NAN;
        }
        out.push(estimate);
    }
    out
}

/// The best-scoring discovery of each key (ties to the lower q-value, then target before
/// entrapment, then lexicographic pair key), in key order: a peptide, a precursor or a group
/// is ranked by its best PSM.
pub fn best_per_key(rows: impl IntoIterator<Item = Discovery>) -> Vec<Discovery> {
    let mut best: std::collections::BTreeMap<String, Discovery> = Default::default();
    for d in rows {
        let better = best.get(&d.key).is_none_or(|b| {
            d.score
                .total_cmp(&b.score)
                .then(b.q_value.total_cmp(&d.q_value))
                .then(b.entrapment.cmp(&d.entrapment))
                .then_with(|| b.pair.cmp(&d.pair))
                .is_gt()
        });
        if better {
            best.insert(d.key.clone(), d);
        }
    }
    best.into_values().collect()
}

/// The estimate over exactly the rows whose q-value is at most `q`.
/// `curve` must be the estimates returned for these rows, so the original ratio is retained.
pub fn at_q(rows: &[Discovery], curve: &[Estimate], q: f64) -> Estimate {
    let (target_scores, applicable) = pairing(rows);
    let accepted_targets: std::collections::HashSet<&str> = rows
        .iter()
        .filter(|d| !d.entrapment && d.q_value <= q)
        .filter_map(|d| d.pair.as_deref())
        .collect();
    let (mut t, mut e, mut below, mut lower) = (0, 0, 0, 0);
    for d in rows.iter().filter(|d| d.q_value <= q) {
        if !d.entrapment {
            t += 1;
            continue;
        }
        e += 1;
        match d.pair.as_deref() {
            Some(p) if target_scores.get(p).is_some_and(|score| d.score <= *score) => {}
            Some(p) if accepted_targets.contains(p) => lower += 1,
            _ => below += 1,
        }
    }
    let r = curve.first().map_or(1.0, |estimate| estimate.ratio);
    let mut estimate = Estimate::of(t, e, below, lower, r);
    if !applicable {
        estimate.paired = f64::NAN;
    }
    estimate
}

fn pairing(rows: &[Discovery]) -> (HashMap<&str, f64>, bool) {
    let mut targets = HashMap::new();
    let mut entrapments = std::collections::HashSet::new();
    let mut applicable = true;
    for d in rows {
        let Some(pair) = d.pair.as_deref() else {
            applicable = false;
            continue;
        };
        if pair.is_empty() || !d.score.is_finite() {
            applicable = false;
        }
        if d.entrapment {
            applicable &= entrapments.insert(pair);
        } else if targets.insert(pair, d.score).is_some() {
            applicable = false;
        }
    }
    (targets, applicable)
}

/// One level of the table: its name, the ratio used for the combined estimate, its rows best
/// first and their estimates (from [`curve`]).
pub struct Level<'a> {
    pub name: &'a str,
    pub r: f64,
    pub rows: &'a [Discovery],
    pub curve: &'a [Estimate],
}

/// Writes the entrapment table `path` (`andes.entrapment.schema`): one row per discovery of
/// each level, best first, with `level`, `key`, `entrapment`, `score`, `q_value` and, down to
/// that row, `targets`, `entrapment_discoveries`, `lower_bound_fdp`, `combined_fdp` and
/// `paired_fdp`. Nullable `pair`, `pair_below_cut` and `pair_lower` preserve the pairing
/// inputs and counts for independent reproduction. The metadata holds the seed and the
/// ratio of each level (`r_<level>`).
pub fn write_curve(path: &Path, levels: &[Level<'_>], seed: u64) -> Result<(), String> {
    let schema = Arc::new(Schema::new(vec![
        Field::new("level", DataType::Utf8, false),
        Field::new("key", DataType::Utf8, false),
        Field::new("pair", DataType::Utf8, true),
        Field::new("entrapment", DataType::Boolean, false),
        Field::new("score", DataType::Float64, false),
        Field::new("q_value", DataType::Float64, false),
        Field::new("targets", DataType::UInt64, false),
        Field::new("entrapment_discoveries", DataType::UInt64, false),
        Field::new("pair_below_cut", DataType::UInt64, false),
        Field::new("pair_lower", DataType::UInt64, false),
        Field::new("lower_bound_fdp", DataType::Float64, false),
        Field::new("combined_fdp", DataType::Float64, false),
        Field::new("paired_fdp", DataType::Float64, false),
    ]));
    let batches = levels.iter().map(|l| {
        let n = l.rows.len();
        let cols: Vec<ArrayRef> = vec![
            Arc::new(StringArray::from(vec![l.name; n])),
            Arc::new(StringArray::from_iter_values(
                l.rows.iter().map(|d| d.key.as_str()),
            )),
            Arc::new(StringArray::from_iter(
                l.rows.iter().map(|d| d.pair.as_deref()),
            )),
            Arc::new(BooleanArray::from_iter(
                l.rows.iter().map(|d| Some(d.entrapment)),
            )),
            Arc::new(Float64Array::from_iter_values(
                l.rows.iter().map(|d| d.score),
            )),
            Arc::new(Float64Array::from_iter_values(
                l.rows.iter().map(|d| d.q_value),
            )),
            Arc::new(UInt64Array::from_iter_values(
                l.curve.iter().map(|c| c.targets as u64),
            )),
            Arc::new(UInt64Array::from_iter_values(
                l.curve.iter().map(|c| c.entrapment as u64),
            )),
            Arc::new(UInt64Array::from_iter_values(
                l.curve.iter().map(|c| c.pair_below_cut as u64),
            )),
            Arc::new(UInt64Array::from_iter_values(
                l.curve.iter().map(|c| c.pair_lower as u64),
            )),
            Arc::new(Float64Array::from_iter_values(
                l.curve.iter().map(|c| c.lower_bound),
            )),
            Arc::new(Float64Array::from_iter_values(
                l.curve.iter().map(|c| c.combined),
            )),
            Arc::new(Float64Array::from_iter_values(
                l.curve.iter().map(|c| c.paired),
            )),
        ];
        RecordBatch::try_new(schema.clone(), cols).map_err(|e| e.to_string())
    });
    let mut meta = vec![("entrapment_seed", seed.to_string())];
    let names: Vec<String> = levels.iter().map(|l| format!("r_{}", l.name)).collect();
    for (l, name) in levels.iter().zip(&names) {
        meta.push((name.as_str(), l.r.to_string()));
    }
    write_table(path, &CURVE, CURVE_VERSION, schema.clone(), batches, &meta)
}

#[cfg(test)]
#[path = "tests/entrapment.rs"]
mod tests;
