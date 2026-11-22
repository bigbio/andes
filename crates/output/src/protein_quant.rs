//! Default protein-group quantities of `andes quant finalize`: total intensity and Top3.
//!
//! Both are computed from the rows of the feature view, the same way for label-free and
//! isobaric data:
//! - **Unit and label.** A quantity belongs to one protein group, one quantification unit and one
//!   label. Label-free: the unit is a run and the label is `LFQ`. Isobaric: the unit is a plex
//!   (every fraction of the plex together; a run when the design gives no plex) and the label is a
//!   reporter channel, reference channels included. Values are the feature intensities as written
//!   (label-free: the feature quantity; isobaric: the reporter intensities of the run tables, raw
//!   or impurity-corrected as searched). Nothing is normalised and no channel is divided by the
//!   reference: both quantities are raw sums for downstream normalisation.
//! - **Rows used.** Only features of peptides unique to the group; shared and razor peptides are
//!   left out. Isobaric rows that failed the quantification gate (weight 0) are left out.
//! - **Peptide.** A peptide is a peptidoform (modified sequence); its intensity in a unit and
//!   label is the sum over its charge states, fractions and rows.
//! - **Total intensity** is the sum over the group's peptides of their intensity in the label.
//! - **Top3** is the mean intensity in the label of the group's three most intense peptides in the
//!   unit, ranked by their intensity summed over the unit's labels (so on a plex every channel
//!   averages the same peptides and the channel ratios stay those of the same ions). A group with
//!   fewer than three peptides averages the ones it has. Ties in the ranking are broken by the
//!   peptidoform, ascending; tied peptides have equal intensity, so the choice does not change
//!   the mean.
//! - **Missing values.** A value that is not finite or not above zero is missing. Total
//!   intensity is missing when every peptide is missing in that label; Top3 averages the ranked
//!   peptides with a value and is missing when none of the three has one.

use std::collections::{BTreeMap, HashMap};

/// `intensity_name` of the summed quantity in `additional_intensities`; also the row's
/// primary `intensity`.
pub const TOTAL_INTENSITY: &str = "total_intensity";
/// `intensity_name` of the Top3 quantity in `additional_intensities`.
pub const TOP3: &str = "top3";
/// Version of the default protein quantities, written in the pg view's footer.
pub const PROTEIN_QUANT_VERSION: &str = "1.0";
/// Number of peptides averaged by Top3.
pub const TOP_N: usize = 3;

/// One quantification unit: its name, its runs and its labels in channel order.
#[derive(Debug, Clone, PartialEq)]
pub struct QuantUnit {
    pub name: String,
    pub runs: Vec<String>,
    pub labels: Vec<String>,
}

/// The quantities of one group in one unit and label.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LabelQuantity {
    pub total_intensity: Option<f64>,
    pub top3: Option<f64>,
}

/// Accumulates unique-peptide intensities per (group, unit) and rolls them up.
#[derive(Debug, Default)]
pub struct ProteinQuant {
    units: Vec<QuantUnit>,
    unit_index: HashMap<String, usize>,
    cells: HashMap<(String, usize), BTreeMap<String, Vec<f64>>>,
}

fn present(v: f64) -> bool {
    v.is_finite() && v > 0.0
}

impl ProteinQuant {
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds one feature row of `group` (its anchor accession): its peptidoform and its
    /// intensity per label, in the unit `unit` from the run `run`. The first row of a unit sets
    /// its labels; a later row with other labels is an error.
    pub fn add(
        &mut self,
        group: &str,
        unit: &str,
        run: &str,
        peptidoform: &str,
        labels: &[String],
        values: &[f32],
    ) -> Result<(), String> {
        if labels.len() != values.len() {
            return Err(format!(
                "protein quantity: {} labels for {} values",
                labels.len(),
                values.len()
            ));
        }
        let u = match self.unit_index.get(unit) {
            Some(&u) => {
                if self.units[u].labels != labels {
                    return Err(format!(
                        "protein quantity: unit {unit} has labels {:?} and {:?}",
                        self.units[u].labels, labels
                    ));
                }
                u
            }
            None => {
                self.units.push(QuantUnit {
                    name: unit.to_string(),
                    runs: Vec::new(),
                    labels: labels.to_vec(),
                });
                self.unit_index.insert(unit.to_string(), self.units.len() - 1);
                self.units.len() - 1
            }
        };
        let runs = &mut self.units[u].runs;
        if !runs.iter().any(|r| r == run) {
            runs.push(run.to_string());
        }
        let peptide = self
            .cells
            .entry((group.to_string(), u))
            .or_default()
            .entry(peptidoform.to_string())
            .or_insert_with(|| vec![0.0; labels.len()]);
        for (acc, &v) in peptide.iter_mut().zip(values) {
            let v = f64::from(v);
            if present(v) {
                *acc += v;
            }
        }
        Ok(())
    }

    pub fn is_empty(&self) -> bool {
        self.cells.is_empty()
    }

    /// The units, each with its runs sorted.
    pub fn units(&self) -> Vec<QuantUnit> {
        let mut out = self.units.clone();
        for u in &mut out {
            u.runs.sort();
        }
        out
    }

    /// The quantities of `group` per unit (index into [`Self::units`]) and label; units where
    /// the group has no row are absent.
    pub fn group(&self, group: &str) -> Vec<(usize, Vec<LabelQuantity>)> {
        let mut out: Vec<(usize, Vec<LabelQuantity>)> = (0..self.units.len())
            .filter_map(|u| {
                self.cells
                    .get(&(group.to_string(), u))
                    .map(|peptides| (u, roll_up(peptides, self.units[u].labels.len())))
            })
            .collect();
        out.sort_by_key(|(u, _)| *u);
        out
    }
}

/// Total intensity and Top3 per label from peptide intensities (module doc).
pub fn roll_up(peptides: &BTreeMap<String, Vec<f64>>, n_labels: usize) -> Vec<LabelQuantity> {
    let mut ranked: Vec<(&String, f64)> = peptides
        .iter()
        .map(|(p, v)| (p, v.iter().copied().filter(|&x| present(x)).sum::<f64>()))
        .filter(|(_, s)| *s > 0.0)
        .collect();
    ranked.sort_by(|a, b| b.1.total_cmp(&a.1).then_with(|| a.0.cmp(b.0)));
    let top: Vec<&Vec<f64>> = ranked
        .iter()
        .take(TOP_N)
        .map(|(p, _)| &peptides[*p])
        .collect();
    (0..n_labels)
        .map(|l| {
            let all: Vec<f64> = peptides
                .values()
                .map(|v| v[l])
                .filter(|&x| present(x))
                .collect();
            let best: Vec<f64> = top.iter().map(|v| v[l]).filter(|&x| present(x)).collect();
            LabelQuantity {
                total_intensity: (!all.is_empty()).then(|| all.iter().sum()),
                top3: (!best.is_empty()).then(|| best.iter().sum::<f64>() / best.len() as f64),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn labels(n: &[&str]) -> Vec<String> {
        n.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn total_is_the_sum_and_top3_the_mean_of_the_three_most_intense_peptides() {
        let mut q = ProteinQuant::new();
        let l = labels(&["LFQ"]);
        for (p, v) in [("A", 10.0), ("B", 40.0), ("C", 30.0), ("D", 20.0)] {
            q.add("P1", "run1", "run1", p, &l, &[v]).unwrap();
        }
        let g = q.group("P1");
        assert_eq!(g.len(), 1);
        assert_eq!(g[0].1[0].total_intensity, Some(100.0));
        assert_eq!(g[0].1[0].top3, Some(30.0));
    }

    #[test]
    fn a_peptide_sums_its_charge_states_and_fractions_before_ranking() {
        let mut q = ProteinQuant::new();
        let l = labels(&["LFQ"]);
        q.add("P1", "u", "fr1", "A", &l, &[5.0]).unwrap();
        q.add("P1", "u", "fr2", "A", &l, &[50.0]).unwrap();
        q.add("P1", "u", "fr1", "B", &l, &[30.0]).unwrap();
        let g = q.group("P1");
        assert_eq!(g[0].1[0].total_intensity, Some(85.0));
        assert_eq!(g[0].1[0].top3, Some(42.5));
        assert_eq!(q.units()[0].runs, vec!["fr1".to_string(), "fr2".to_string()]);
    }

    #[test]
    fn fewer_than_three_peptides_average_the_ones_there_are() {
        let mut q = ProteinQuant::new();
        let l = labels(&["LFQ"]);
        q.add("P1", "r", "r", "A", &l, &[8.0]).unwrap();
        assert_eq!(q.group("P1")[0].1[0].top3, Some(8.0));
    }

    #[test]
    fn plex_channels_average_the_same_peptides_ranked_on_the_plex_total() {
        let mut q = ProteinQuant::new();
        let l = labels(&["TMT126", "TMT127N"]);
        q.add("P1", "plex1", "f1", "A", &l, &[100.0, 1.0]).unwrap();
        q.add("P1", "plex1", "f1", "B", &l, &[10.0, 60.0]).unwrap();
        q.add("P1", "plex1", "f1", "C", &l, &[10.0, 50.0]).unwrap();
        q.add("P1", "plex1", "f1", "D", &l, &[5.0, 55.0]).unwrap();
        let g = &q.group("P1")[0].1;
        assert_eq!(g[0].top3, Some(40.0));
        assert_eq!(g[1].top3, Some(37.0));
        assert_eq!(g[0].total_intensity, Some(125.0));
        assert_eq!(g[1].total_intensity, Some(166.0));
    }

    #[test]
    fn zero_and_non_finite_values_are_missing() {
        let mut q = ProteinQuant::new();
        let l = labels(&["TMT126", "TMT127N"]);
        q.add("P1", "p", "f", "A", &l, &[0.0, f32::NAN]).unwrap();
        q.add("P1", "p", "f", "B", &l, &[4.0, 0.0]).unwrap();
        let g = &q.group("P1")[0].1;
        assert_eq!(g[0].total_intensity, Some(4.0));
        assert_eq!(g[0].top3, Some(4.0));
        assert_eq!(g[1].total_intensity, None);
        assert_eq!(g[1].top3, None);
    }

    #[test]
    fn ties_are_broken_by_peptidoform_and_do_not_change_the_mean() {
        let mut q = ProteinQuant::new();
        let l = labels(&["LFQ"]);
        for p in ["D", "C", "B", "A"] {
            q.add("P1", "r", "r", p, &l, &[7.0]).unwrap();
        }
        assert_eq!(q.group("P1")[0].1[0].top3, Some(7.0));
    }

    #[test]
    fn units_are_separate_and_labels_must_agree_within_a_unit() {
        let mut q = ProteinQuant::new();
        let l = labels(&["LFQ"]);
        q.add("P1", "r1", "r1", "A", &l, &[1.0]).unwrap();
        q.add("P1", "r2", "r2", "A", &l, &[3.0]).unwrap();
        let g = q.group("P1");
        assert_eq!(g.len(), 2);
        assert_eq!(g[1].1[0].total_intensity, Some(3.0));
        assert!(q.group("P2").is_empty());
        assert!(q
            .add("P1", "r1", "r1", "A", &labels(&["TMT126"]), &[1.0])
            .is_err());
        assert!(q.add("P1", "r1", "r1", "A", &l, &[1.0, 2.0]).is_err());
    }
}
