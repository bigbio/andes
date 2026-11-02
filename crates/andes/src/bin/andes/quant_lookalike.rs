//! Look-alike precursors of match-between-runs: transfers that land on the
//! identified peak of another precursor of the same m/z.
//!
//! Two target library entries are look-alikes when their m/z differ by at
//! most the extract's MS1 tolerance and their consensus RTs (the median RT
//! of the donor PSMs for an entry without one) by at most
//! [`LOOKALIKE_RT_S`]: in a run, the template fits of both read the same
//! isotope-0 trace over the same stretch of elution time. I/L swaps, N vs GG,
//! Q vs AG, reordered sequences and conserved homologues of another species
//! are such pairs. The mass- and RT-shift decoys of the transfer FDR model a
//! transfer that lands on a peak by chance; a look-alike's peak is there
//! because of the sequence, so they cannot model it, and a transfer of an
//! absent precursor onto it looks like a good measurement to every
//! classifier.
//!
//! A transfer (or pseudo-transfer) of precursor `P` in run `r` *shares a
//! peak* with `Q` when `Q` is a look-alike of `P` identified in `r` with a
//! fitted peak (positive amplitude) and the transfer's fitted apex lies
//! within `Q`'s identified width of `Q`'s identified apex,
//! `|t0_P - t0_Q| <= w_Q` (the wrong-peak rule of the pseudo-transfers). With
//! several such `Q`, the one whose apex is nearest is reported (lowest id on a
//! tie). The identification of `Q` is direct MS2 evidence that the peak is
//! `Q`'s; an MS2 scan that supports `P` is evidence for `P` too, unless the
//! two sequences differ only by I/L, which MS2 cannot tell apart
//! ([`Lookalikes::ms2_tells_apart`]).

use std::collections::{HashMap, HashSet};

use quant::library::Library;

/// Largest difference of consensus RTs, seconds, between two look-alikes: a
/// pre-filter only (the test is on the fitted apexes in each run), wide
/// enough for the transfer's RT search and the alignment error.
pub(crate) const LOOKALIKE_RT_S: f64 = 120.0;

/// The look-alike pairs of a library's target entries.
#[derive(Debug, Default)]
pub(crate) struct Lookalikes {
    /// Per target precursor with look-alikes, their ids, ascending.
    partners: HashMap<u32, Vec<u32>>,
    /// Pairs `(min, max)` that MS2 cannot tell apart (I/L only).
    ms2_blind: HashSet<(u32, u32)>,
}

/// The RT of an entry for the look-alike pre-filter: its consensus RT, or
/// without one (no observation in an aligned run) the median RT of its
/// donor runs' best PSMs; `None` when neither is finite.
fn entry_rt(e: &quant::library::LibraryEntry) -> Option<f64> {
    if e.consensus_rt.is_finite() {
        return Some(e.consensus_rt);
    }
    let mut v: Vec<f64> = e
        .best_rt
        .iter()
        .map(|r| r.1)
        .filter(|t| t.is_finite())
        .collect();
    if v.is_empty() {
        return None;
    }
    v.sort_unstable_by(f64::total_cmp);
    let n = v.len();
    Some(if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    })
}

/// The residues of a `<peptidoform>/<charge>` key with I read as L; the
/// modification brackets and the charge are kept as they are.
fn il_blind(key: &str) -> String {
    let mut depth = 0i32;
    key.chars()
        .map(|c| {
            match c {
                '[' | '(' | '{' => depth += 1,
                ']' | ')' | '}' => depth -= 1,
                'I' if depth == 0 => return 'L',
                _ => {}
            }
            c
        })
        .collect()
}

impl Lookalikes {
    /// The look-alike pairs of `lib`'s target entries at m/z tolerance
    /// `tol_ppm` (relative to the higher m/z of the pair); identity decoys and
    /// entries without a finite consensus RT have none.
    pub(crate) fn new(lib: &Library, tol_ppm: f64) -> Self {
        let mut out = Self::default();
        if !(tol_ppm.is_finite() && tol_ppm > 0.0) {
            return out;
        }
        // Target entries with an RT, ascending by m/z.
        let mut ents: Vec<(f64, f64, u32)> = lib
            .entries
            .iter()
            .filter(|e| !e.identity_decoy)
            .filter_map(|e| {
                let mz = lib.meta.get(e.precursor as usize)?.mz;
                let rt = entry_rt(e)?;
                mz.is_finite().then_some((mz, rt, e.precursor))
            })
            .collect();
        ents.sort_unstable_by(|a, b| a.0.total_cmp(&b.0).then(a.2.cmp(&b.2)));
        for i in 0..ents.len() {
            let (mz_i, rt_i, p) = ents[i];
            for &(mz_j, rt_j, q) in &ents[i + 1..] {
                if mz_j - mz_i > tol_ppm * 1e-6 * mz_j {
                    break;
                }
                if (rt_j - rt_i).abs() > LOOKALIKE_RT_S || p == q {
                    continue;
                }
                out.partners.entry(p).or_default().push(q);
                out.partners.entry(q).or_default().push(p);
                let (kp, kq) = (&lib.meta[p as usize].key, &lib.meta[q as usize].key);
                if il_blind(kp) == il_blind(kq) {
                    out.ms2_blind.insert((p.min(q), p.max(q)));
                }
            }
        }
        for v in out.partners.values_mut() {
            v.sort_unstable();
            v.dedup();
        }
        out
    }

    /// Precursors with at least one look-alike.
    pub(crate) fn len(&self) -> usize {
        self.partners.len()
    }

    /// Whether `p` is some precursor's look-alike (its identified fits are
    /// needed). Look-alikes are symmetric, so these are the precursors with
    /// look-alikes.
    pub(crate) fn wants(&self, p: u32) -> bool {
        self.partners.contains_key(&p)
    }

    /// The look-alike `Q` whose identified peak in `run` the apex `t0` of a
    /// (pseudo-)transfer of `p` lies on (see the module docs); `fit(q, run)`
    /// is `Q`'s identified `(t0, w)` in `run`, `None` without a fitted peak.
    pub(crate) fn shared(
        &self,
        p: u32,
        run: u16,
        t0: f64,
        fit: impl Fn(u32, u16) -> Option<(f64, f64)>,
    ) -> Option<u32> {
        if !t0.is_finite() {
            return None;
        }
        let mut best: Option<(f64, u32)> = None;
        for &q in self.partners.get(&p)? {
            let Some((t_q, w_q)) = fit(q, run) else {
                continue;
            };
            let d = (t0 - t_q).abs();
            if d <= w_q && best.is_none_or(|(b, _)| d < b) {
                best = Some((d, q));
            }
        }
        best.map(|(_, q)| q)
    }

    /// Whether an MS2 scan that supports `p` is evidence for `p` against
    /// its look-alike `q`: false when the two differ only by I/L.
    pub(crate) fn ms2_tells_apart(&self, p: u32, q: u32) -> bool {
        !self.ms2_blind.contains(&(p.min(q), p.max(q)))
    }
}

/// A look-alike's identified fit usable by [`Lookalikes::shared`]: a
/// positive amplitude and a finite apex and width.
pub(crate) fn usable_fit(e: &quant::evidence::Evidence) -> Option<(f64, f64)> {
    (e.amplitude > 0.0 && e.t0.is_finite() && e.w.is_finite() && e.w > 0.0).then_some((e.t0, e.w))
}

#[cfg(test)]
mod tests {
    use super::*;
    use quant::align::AlignParams;
    use quant::library::{build_library, DonorPsm, LibraryEntry, PrecursorMeta};

    fn lib(items: &[(&str, f64, f64, bool)]) -> Library {
        let meta: Vec<PrecursorMeta> = items
            .iter()
            .map(|&(key, mz, _, _)| PrecursorMeta {
                key: key.to_string(),
                mz,
                charge: 2,
                envelope: [1.0, 0.5, 0.2, 0.05],
            })
            .collect();
        let entries = items
            .iter()
            .enumerate()
            .map(|(i, &(_, _, rt, identity))| LibraryEntry {
                precursor: i as u32,
                identity_decoy: identity,
                q: 0.0,
                donor_psm: DonorPsm {
                    score: 1.0,
                    q: 0.0,
                    pep: 0.0,
                },
                consensus_rt: rt,
                faims_cv: None,
                ion_mobility: None,
                donors: vec![0],
                best_rt: vec![(0, rt)],
            })
            .collect();
        let mut lib = build_library(&[], meta, 1, 0.01, 0.01, None, &AlignParams::default());
        lib.entries = entries;
        lib
    }

    #[test]
    fn look_alikes_are_close_in_mz_and_consensus_rt() {
        let l = lib(&[
            ("PEPTLDE/2", 500.0, 1000.0, false),
            ("PEPTIDE/2", 500.003, 1030.0, false), // 6 ppm, 30 s: look-alike of 0
            ("AAAK/2", 500.004, 1500.0, false),    // RT too far
            ("CCCK/2", 500.02, 1000.0, false),     // 40 ppm
            ("DECOY/2", 500.001, 1000.0, true),    // identity decoy
            ("NORT/2", 500.002, f64::NAN, false),  // no RT at all
        ]);
        let la = Lookalikes::new(&l, 10.0);
        assert_eq!(la.partners.get(&0), Some(&vec![1]));
        assert_eq!(la.partners.get(&1), Some(&vec![0]));
        assert!(la.wants(0) && la.wants(1) && !la.wants(2) && !la.wants(3));
        assert!(!la.wants(4) && !la.wants(5));
        assert_eq!(la.len(), 2);
        // I/L only: MS2 cannot tell them apart.
        assert!(!la.ms2_tells_apart(0, 1) && !la.ms2_tells_apart(1, 0));
        assert!(Lookalikes::new(&l, 0.0).partners.is_empty());
    }

    #[test]
    fn an_entry_without_a_consensus_rt_uses_its_donor_psms() {
        let mut l = lib(&[
            ("AK/2", 500.0, 1000.0, false),
            ("BK/2", 500.001, 1000.0, false),
        ]);
        l.entries[1].consensus_rt = f64::NAN;
        l.entries[1].best_rt = vec![(0, 990.0), (1, 1010.0), (2, 2000.0)];
        assert_eq!(Lookalikes::new(&l, 10.0).partners.get(&0), Some(&vec![1]));
        l.entries[1].best_rt = vec![(0, 1500.0)];
        assert!(Lookalikes::new(&l, 10.0).partners.is_empty());
    }

    #[test]
    fn i_l_blindness_leaves_modifications_and_other_residues() {
        assert_eq!(il_blind("PEPTIDE/2"), "PEPTLDE/2");
        assert_eq!(il_blind("M[UNIMOD:35]IK/2"), "M[UNIMOD:35]LK/2");
        assert_ne!(il_blind("PEPTIDE/2"), il_blind("PEPTIDE/3"));
        let l = lib(&[
            ("NGGK/2", 500.0, 1000.0, false),
            ("NNK/2", 500.0, 1000.0, false),
        ]);
        let la = Lookalikes::new(&l, 10.0);
        assert!(la.ms2_tells_apart(0, 1));
    }

    #[test]
    fn a_transfer_shares_the_nearest_identified_look_alike_peak() {
        let l = lib(&[
            ("AK/2", 500.0, 1000.0, false),
            ("BK/2", 500.001, 1000.0, false),
            ("CK/2", 500.002, 1000.0, false),
        ]);
        let la = Lookalikes::new(&l, 10.0);
        let fits = |q: u32, run: u16| match (q, run) {
            (1, 0) => Some((1004.0, 5.0)),
            (2, 0) => Some((1001.0, 2.0)),
            _ => None,
        };
        // On both peaks: the nearer apex wins.
        assert_eq!(la.shared(0, 0, 1000.5, fits), Some(2));
        // Outside 2's width, inside 1's.
        assert_eq!(la.shared(0, 0, 1006.0, fits), Some(1));
        // Outside both.
        assert_eq!(la.shared(0, 0, 1010.0, fits), None);
        // Another run: no identified look-alike.
        assert_eq!(la.shared(0, 1, 1000.5, fits), None);
        assert_eq!(la.shared(0, 0, f64::NAN, fits), None);
    }
}
