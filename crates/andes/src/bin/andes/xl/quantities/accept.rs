//! Which cross-link spectrum matches are quantified: the target matches whose own q-value and
//! whose residue pair's q-value are both at or below `--xl-fdr`. The q-values are those of
//! this search, or those of the pooled tables of `andes cohort xl-fdr` (`--xl-cohort-fdr`),
//! which rescores the matches of all runs together.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use xl::decoy::PairClass;

use super::Search;

/// The q-values `andes cohort xl-fdr` wrote: of each residue pair, and of each match by the
/// title of its spectrum.
pub(crate) struct CohortFdr {
    residue_pairs: HashMap<String, f64>,
    matches: HashMap<String, f64>,
}

fn with_suffix(base: &Path, suffix: &str) -> PathBuf {
    let mut s = base.to_path_buf().into_os_string();
    s.push(suffix);
    PathBuf::from(s)
}

/// The `key` and `q_value` columns of the target rows of a table of `andes cohort xl-fdr`.
fn target_qvalues(path: &Path, key: &str) -> Result<HashMap<String, f64>, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("--xl-cohort-fdr: cannot read {}: {e}", path.display()))?;
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
    let column = |name: &str| {
        header.iter().position(|h| *h == name).ok_or_else(|| {
            format!(
                "--xl-cohort-fdr: {} has no `{name}` column; it is not a table of \
                 `andes cohort xl-fdr`",
                path.display()
            )
        })
    };
    let (key, class, q) = (column(key)?, column("class")?, column("q_value")?);
    let mut out = HashMap::new();
    for line in lines {
        let cells: Vec<&str> = line.split('\t').collect();
        let (Some(k), Some(c), Some(v)) = (cells.get(key), cells.get(class), cells.get(q)) else {
            continue;
        };
        if *c != PairClass::TT.name() {
            continue;
        }
        if let Ok(v) = v.parse::<f64>() {
            out.insert((*k).to_string(), v);
        }
    }
    Ok(out)
}

impl CohortFdr {
    /// Reads `<base>.xl.residue_pairs.tsv` and `<base>.xl.tsv`.
    pub(crate) fn read(base: &Path) -> Result<Self, String> {
        Ok(Self {
            residue_pairs: target_qvalues(
                &with_suffix(base, ".xl.residue_pairs.tsv"),
                "residue_pair",
            )?,
            matches: target_qvalues(&with_suffix(base, ".xl.tsv"), "spectrum")?,
        })
    }
}

/// For each match of `search.levels` (in the order of its `csm_items`), the q-value of the
/// match and of its residue pair when it is a target with both; `None` otherwise.
pub(super) fn qvalues(search: &Search<'_>, cohort: Option<&CohortFdr>) -> Vec<Option<(f64, f64)>> {
    let levels = search.levels;
    let own_pairs: HashMap<&str, f64> = levels
        .rp
        .iter()
        .zip(&levels.rp_q)
        .filter_map(|(g, q)| q.map(|q| (g.key.as_str(), q)))
        .collect();
    levels
        .csm_items
        .iter()
        .zip(&levels.csm_rows)
        .zip(&levels.csm_q)
        .map(|((item, &row), own_q)| {
            if item.class != PairClass::TT {
                return None;
            }
            match cohort {
                None => Some(((*own_q)?, *own_pairs.get(item.key.as_str())?)),
                Some(c) => {
                    let spec = search.spectra.get(search.results[row].spec_idx)?;
                    Some((
                        *c.matches.get(spec.title.as_str())?,
                        *c.residue_pairs.get(item.key.as_str())?,
                    ))
                }
            }
        })
        .collect()
}

#[cfg(test)]
#[path = "tests/accept.rs"]
mod tests;
