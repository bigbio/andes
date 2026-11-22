//! The spectra a search read, per run and per precursor charge (spectra without a charge apart),
//! for the run summary's identification rates.

use std::collections::BTreeMap;

use model::spectrum::Spectrum;

#[derive(Debug, Clone, Default, PartialEq)]
pub struct SpectraTally {
    pub per_run: Vec<usize>,
    pub by_charge: BTreeMap<u8, usize>,
    pub unknown_charge: usize,
}

impl SpectraTally {
    pub fn of(spectra: &[Spectrum], run_ends: &[usize]) -> Self {
        let mut per_run = Vec::with_capacity(run_ends.len());
        let mut start = 0usize;
        for &end in run_ends {
            per_run.push(end.saturating_sub(start));
            start = start.max(end);
        }
        let mut by_charge: BTreeMap<u8, usize> = BTreeMap::new();
        let mut unknown_charge = 0usize;
        for s in spectra {
            match s.precursor_charge.and_then(|z| u8::try_from(z).ok()) {
                Some(z) if z > 0 => *by_charge.entry(z).or_insert(0) += 1,
                _ => unknown_charge += 1,
            }
        }
        Self {
            per_run,
            by_charge,
            unknown_charge,
        }
    }

    pub(crate) fn total(&self) -> usize {
        self.by_charge.values().sum::<usize>() + self.unknown_charge
    }
}
