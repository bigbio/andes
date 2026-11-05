//! Shared-feature evidence and downgrade: does the MS1 feature of a transfer
//! belong to another compound?
//!
//! A transfer whose fitted feature is also claimed by another probed target
//! of the run ([`crate::claims`]) is either a second real compound on that
//! feature (isobaric residue combinations, permutations) or a feature that
//! belongs to the claimant alone. The MS2 scans of the receiving run that
//! cover the feature tell them apart in part: fragment ions that only the
//! transfer's sequence has (and not the claimant's) are present or absent,
//! and the covering scans were identified as the transfer's peptide, as a
//! claimant's, as something else, or not at all ([`SharedEvidence`]).
//!
//! A small gradient-boosted classifier ([`Model`]) turns these facts and the
//! claim geometry into the probability that an accepted transfer is a shared
//! feature; finalize reports such a transfer as `shared_feature` instead of
//! `transfer_ms1` / `transfer_ms2`. The model is a shipped resource with its
//! provenance (`resources/quant/shared_feature_model.txt`).

use std::collections::HashMap;

use model::Tolerance;

use crate::claims::Claim;
use crate::ms2::{Fragment, Ms2RunIndex};

/// Covering scans of a feature: within this many fitted widths of the apex
/// (and not closer than [`COVER_MIN_S`]).
pub const COVER_WIDTHS: f64 = 2.0;
/// Least half-width of the covering window, seconds.
pub const COVER_MIN_S: f64 = 5.0;

/// The MS2 evidence of one transfer row against its claimants.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct SharedEvidence {
    /// Scans covering the feature.
    pub n_scans: u16,
    /// Of them: no identification, identified as the transfer's peptide, as
    /// a claimant's, as another peptide.
    pub n_unid: u16,
    pub n_id_target: u16,
    pub n_id_claimant: u16,
    pub n_id_other: u16,
    /// Fragment ions of the transfer's sequence absent from the claimant
    /// with the best-scoring pair: ions found in the covering scans, and
    /// ions looked for (ions x scans); the same for the claimant's own.
    pub hit_target: u32,
    pub n_target: u32,
    pub hit_claimant: u32,
    pub n_claimant: u32,
    /// Summed share of the scan intensity carried by the matched ions.
    pub int_target: f32,
    pub int_claimant: f32,
}

impl SharedEvidence {
    /// A claimant with distinguishing ions on both sides was scored.
    pub fn scored(&self) -> bool {
        self.n_target > 0 && self.n_claimant > 0
    }

    /// Share of the transfer's own ions found minus the claimant's; 0 when
    /// not scored.
    pub fn dfrac(&self) -> f64 {
        if !self.scored() {
            return 0.0;
        }
        f64::from(self.hit_target) / f64::from(self.n_target)
            - f64::from(self.hit_claimant) / f64::from(self.n_claimant)
    }

    /// Intensity share of the transfer's own ions minus the claimant's.
    pub fn dint(&self) -> f64 {
        if !self.scored() {
            return 0.0;
        }
        f64::from(self.int_target) - f64::from(self.int_claimant)
    }
}

/// The peptides identified at the MS2 scans of one run, by scan retention
/// time (the exact value of the donor PSM).
#[derive(Debug, Clone, Default)]
pub struct ScanIdentity {
    by_rt: HashMap<u64, Vec<u32>>,
}

impl ScanIdentity {
    /// Records that peptide `peptide` was identified at the scan at `rt`.
    pub fn add(&mut self, rt: f64, peptide: u32) {
        let v = self.by_rt.entry(rt.to_bits()).or_default();
        if !v.contains(&peptide) {
            v.push(peptide);
        }
    }

    pub fn at(&self, rt: f64) -> &[u32] {
        self.by_rt.get(&rt.to_bits()).map_or(&[], Vec::as_slice)
    }
}

/// The fragment ladder of a precursor: its peptide id and the sorted m/z of
/// its b/y ions.
#[derive(Debug, Clone, PartialEq)]
pub struct Ions {
    pub peptide: u32,
    pub mz: Vec<f64>,
}

impl Ions {
    pub fn new(peptide: u32, ladder: &[Fragment]) -> Self {
        let mut mz: Vec<f64> = ladder.iter().map(|f| f.mz).collect();
        mz.sort_by(f64::total_cmp);
        Self { peptide, mz }
    }
}

/// Ions of `a` with no ion of `b` within `tol`.
fn distinguishing(a: &[f64], b: &[f64], tol: Tolerance) -> Vec<f64> {
    a.iter()
        .copied()
        .filter(|&x| {
            let d = tol.as_da(x);
            let i = b.partition_point(|&y| y < x - d);
            !b.get(i).is_some_and(|&y| y <= x + d)
        })
        .collect()
}

/// Ions found in the sorted peak list, and the summed intensity of the most
/// intense peak within `tol` of each.
fn match_ions(peaks: &[(f32, f32)], ions: &[f64], tol: Tolerance) -> (u32, f64) {
    let (mut hits, mut sum) = (0u32, 0.0f64);
    for &x in ions {
        let d = tol.as_da(x);
        let start = peaks.partition_point(|p| f64::from(p.0) < x - d);
        let best = peaks[start..]
            .iter()
            .take_while(|p| f64::from(p.0) <= x + d)
            .map(|p| p.1)
            .filter(|i| i.is_finite() && *i > 0.0)
            .max_by(f32::total_cmp);
        if let Some(i) = best {
            hits += 1;
            sum += f64::from(i);
        }
    }
    (hits, sum)
}

/// A transfer probe and its claimants.
#[derive(Debug, Clone, Copy)]
pub struct Query<'a> {
    pub mz: f64,
    /// Fitted apex and width, seconds.
    pub t0: f64,
    pub w: f64,
    pub faims: Option<f32>,
    pub target: &'a Ions,
    pub claimants: &'a [&'a Ions],
}

/// The evidence of probe `q` in the MS2 scans of run `index`, whose scans
/// `identity` says were identified.
pub fn evidence(
    index: &Ms2RunIndex,
    identity: &ScanIdentity,
    q: &Query<'_>,
    tol: Tolerance,
) -> SharedEvidence {
    let Query {
        mz,
        t0,
        w,
        faims,
        target,
        claimants,
    } = *q;
    let half = (COVER_WIDTHS * w).max(COVER_MIN_S);
    let scans = index.covering(mz, (t0 - half, t0 + half), faims);
    let mut out = SharedEvidence {
        n_scans: u16::try_from(scans.len()).unwrap_or(u16::MAX),
        ..SharedEvidence::default()
    };
    if scans.is_empty() {
        return out;
    }
    for s in &scans {
        let ids = identity.at(s.rt);
        if ids.is_empty() {
            out.n_unid += 1;
        } else if ids.contains(&target.peptide) {
            out.n_id_target += 1;
        } else if ids
            .iter()
            .any(|i| claimants.iter().any(|c| c.peptide == *i))
        {
            out.n_id_claimant += 1;
        } else {
            out.n_id_other += 1;
        }
    }
    let totals: Vec<f64> = scans
        .iter()
        .map(|s| s.peaks.iter().map(|p| f64::from(p.1)).sum())
        .collect();
    let n = scans.len() as u32;
    let mut best: Option<(f64, SharedEvidence)> = None;
    for c in claimants {
        let dt = distinguishing(&target.mz, &c.mz, tol);
        let dc = distinguishing(&c.mz, &target.mz, tol);
        if dt.is_empty() || dc.is_empty() {
            continue;
        }
        let (mut ht, mut hc, mut it, mut ic) = (0u32, 0u32, 0.0f64, 0.0f64);
        for (s, total) in scans.iter().zip(&totals) {
            let (a, ai) = match_ions(&s.peaks, &dt, tol);
            let (b, bi) = match_ions(&s.peaks, &dc, tol);
            ht += a;
            hc += b;
            if *total > 0.0 {
                it += ai / total;
                ic += bi / total;
            }
        }
        let cand = SharedEvidence {
            hit_target: ht,
            n_target: dt.len() as u32 * n,
            hit_claimant: hc,
            n_claimant: dc.len() as u32 * n,
            int_target: it as f32,
            int_claimant: ic as f32,
            ..out
        };
        let d = cand.dfrac();
        if best.as_ref().is_none_or(|(b, _)| d < *b) {
            best = Some((d, cand));
        }
    }
    best.map_or(out, |(_, e)| e)
}

/// The features of the shared-feature model, in the order of the shipped
/// model's `features` line (a model names the ones it uses).
pub const FEATURES: [&str; 21] = [
    "isid", "co", "min_dt", "n_claim", "n_id_co", "n_any_co", "n_tr_co", "n_scans", "scored",
    "dfrac", "dint", "fT", "fP", "n_unid", "n_idT", "n_idP", "n_idO", "A", "logq", "loglr", "w",
];

/// What the features are computed from: one accepted transfer row.
#[derive(Debug, Clone, Copy)]
pub struct Row<'a> {
    pub claim: &'a Claim,
    pub shared: Option<&'a SharedEvidence>,
    /// An MS2 scan of the run supports the transfer (`Ms2Evidence::supports`).
    pub ms2_supported: bool,
    pub transfer_q: f64,
    pub llr: f64,
    /// Fitted width, seconds.
    pub width: f64,
}

/// The value of feature `name` of `row`; `None` for an unknown name.
pub fn feature(name: &str, row: &Row<'_>) -> Option<f64> {
    let c = row.claim;
    let e = row.shared.copied().unwrap_or_default();
    let ratio = |h: u32, n: u32| {
        if n == 0 {
            0.0
        } else {
            f64::from(h) / f64::from(n)
        }
    };
    Some(match name {
        "isid" => f64::from(u8::from(c.identified)),
        "co" => f64::from(u8::from(c.co_apex)),
        "min_dt" => c.min_abs_dt_s.map_or(0.0, f64::from),
        "n_claim" => f64::from(c.n),
        "n_id_co" => f64::from(c.n_identified_co_apex),
        "n_any_co" => f64::from(c.n_co_apex),
        "n_tr_co" => f64::from(c.n_co_apex.saturating_sub(c.n_identified_co_apex)),
        "n_scans" => f64::from(e.n_scans),
        "scored" => f64::from(u8::from(e.scored())),
        "dfrac" => e.dfrac(),
        "dint" => e.dint(),
        "fT" => ratio(e.hit_target, e.n_target),
        "fP" => ratio(e.hit_claimant, e.n_claimant),
        "n_unid" => f64::from(e.n_unid),
        "n_idT" => f64::from(e.n_id_target),
        "n_idP" => f64::from(e.n_id_claimant),
        "n_idO" => f64::from(e.n_id_other),
        "A" => f64::from(u8::from(row.ms2_supported)),
        "logq" => row.transfer_q.clamp(1e-6, 1.0).log10(),
        "loglr" => (1.0 + row.llr.max(0.0)).log10(),
        "w" => row.width,
        _ => return None,
    })
}

/// One regression tree: split nodes `(feature, threshold, left, right)`
/// (a value goes left when `x <= threshold`) and leaves.
#[derive(Debug, Clone, PartialEq)]
enum Node {
    Split(usize, f64, usize, usize),
    Leaf(f64),
}

/// The shipped shared-feature classifier.
#[derive(Debug, Clone, PartialEq)]
pub struct Model {
    /// The model's own description (the `#` lines of the resource).
    pub provenance: Vec<String>,
    features: Vec<String>,
    trees: Vec<Vec<Node>>,
    /// A transfer is a shared feature when its probability exceeds this.
    pub threshold: f64,
}

impl Model {
    /// Parses the text resource. Errors name the offending line.
    pub fn parse(text: &str) -> Result<Self, String> {
        let mut provenance = Vec::new();
        let (mut features, mut threshold) = (Vec::new(), None);
        let mut trees: Vec<Vec<Node>> = Vec::new();
        for (n, line) in text.lines().enumerate() {
            let at = |m: &str| format!("shared-feature model line {}: {m}", n + 1);
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            if let Some(p) = line.strip_prefix('#') {
                provenance.push(p.trim().to_string());
            } else if let Some(f) = line.strip_prefix("features:") {
                features = f.split_whitespace().map(str::to_string).collect();
            } else if let Some(t) = line.strip_prefix("threshold:") {
                threshold = Some(t.trim().parse::<f64>().map_err(|e| at(&e.to_string()))?);
            } else if line.starts_with("tree") {
                trees.push(Vec::new());
            } else {
                let tree = trees.last_mut().ok_or_else(|| at("node before any tree"))?;
                let f: Vec<&str> = line.split_whitespace().collect();
                let num = |i: usize| -> Result<f64, String> {
                    f.get(i)
                        .ok_or_else(|| at("too few fields"))?
                        .parse::<f64>()
                        .map_err(|e| at(&e.to_string()))
                };
                tree.push(match f[0] {
                    "S" => Node::Split(
                        num(1)? as usize,
                        num(2)?,
                        num(3)? as usize,
                        num(4)? as usize,
                    ),
                    "L" => Node::Leaf(num(1)?),
                    other => return Err(at(&format!("unknown node kind {other:?}"))),
                });
            }
        }
        let threshold = threshold.ok_or("shared-feature model: no threshold line")?;
        for f in &features {
            let probe = Row {
                claim: &Claim::default(),
                shared: None,
                ms2_supported: false,
                transfer_q: 0.01,
                llr: 0.0,
                width: 1.0,
            };
            if feature(f, &probe).is_none() {
                return Err(format!("shared-feature model: unknown feature {f:?}"));
            }
        }
        for (k, t) in trees.iter().enumerate() {
            let ok = !t.is_empty()
                && t.iter().all(|n| match n {
                    Node::Split(f, _, l, r) => *f < features.len() && *l < t.len() && *r < t.len(),
                    Node::Leaf(_) => true,
                });
            if !ok {
                return Err(format!("shared-feature model: tree {k} is malformed"));
            }
        }
        Ok(Self {
            provenance,
            features,
            trees,
            threshold,
        })
    }

    /// The model compiled into the binary.
    pub fn shipped() -> &'static Model {
        static MODEL: std::sync::OnceLock<Model> = std::sync::OnceLock::new();
        MODEL.get_or_init(|| {
            Model::parse(include_str!(
                "../../../resources/quant/shared_feature_model.txt"
            ))
            .expect("the shipped shared-feature model parses")
        })
    }

    /// Probability that the transfer `row` is a shared feature.
    pub fn probability(&self, row: &Row<'_>) -> f64 {
        let x: Vec<f64> = self
            .features
            .iter()
            .map(|f| feature(f, row).expect("features are checked at parse"))
            .collect();
        let raw: f64 = self
            .trees
            .iter()
            .map(|t| {
                let mut i = 0;
                loop {
                    match t[i] {
                        Node::Leaf(v) => break v,
                        Node::Split(f, th, l, r) => i = if x[f] <= th { l } else { r },
                    }
                }
            })
            .sum();
        1.0 / (1.0 + (-raw).exp())
    }

    /// The transfer is a shared feature: it has a claimant and the model's
    /// probability exceeds the threshold.
    pub fn is_shared(&self, row: &Row<'_>) -> bool {
        row.claim.n > 0 && !self.trees.is_empty() && self.probability(row) > self.threshold
    }

    /// Whether the model reads the MS2 distinguishing-ion evidence (else the
    /// extract can skip computing it).
    pub fn needs_ms2_evidence(&self) -> bool {
        self.features.iter().any(|f| {
            matches!(
                f.as_str(),
                "n_scans"
                    | "scored"
                    | "dfrac"
                    | "dint"
                    | "fT"
                    | "fP"
                    | "n_unid"
                    | "n_idT"
                    | "n_idP"
                    | "n_idO"
            )
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ms2::{IonKind, Ms2Scan};

    fn scan(rt: f64, peaks: &[(f32, f32)]) -> Ms2Scan {
        Ms2Scan {
            native_id: format!("scan={rt}"),
            spec_key: format!("t_{rt}"),
            rt,
            precursor_mz: 500.0,
            charge: 2,
            isolation: Some((500.0, 1.0, 1.0)),
            faims_cv: None,
            peaks: peaks.to_vec(),
            reporters: None,
        }
    }

    fn frag(mz: f64) -> Fragment {
        Fragment {
            mz,
            ion: crate::ms2::Ion {
                kind: IonKind::B,
                length: 3,
                charge: 1,
            },
        }
    }

    #[test]
    fn target_ions_present_and_claimant_ions_absent() {
        let target = Ions::new(1, &[frag(300.0), frag(400.0), frag(500.0)]);
        let claimant = Ions::new(2, &[frag(300.0), frag(450.0), frag(550.0)]);
        let tol = Tolerance::Ppm(20.0);
        // Distinguishing: target 400 and 500; claimant 450 and 550.
        assert_eq!(
            distinguishing(&target.mz, &claimant.mz, tol),
            [400.0, 500.0]
        );
        assert_eq!(
            distinguishing(&claimant.mz, &target.mz, tol),
            [450.0, 550.0]
        );
        let index = Ms2RunIndex::new(vec![
            scan(100.0, &[(300.0, 5.0), (400.0, 50.0), (500.0, 50.0)]),
            scan(101.0, &[(400.0, 100.0), (450.0, 100.0)]),
            scan(900.0, &[(450.0, 100.0)]),
        ]);
        let mut id = ScanIdentity::default();
        id.add(100.0, 1);
        id.add(101.0, 7);
        let q = Query {
            mz: 500.0,
            t0: 100.5,
            w: 2.0,
            faims: None,
            target: &target,
            claimants: &[&claimant],
        };
        let e = evidence(&index, &id, &q, tol);
        assert_eq!(e.n_scans, 2);
        assert_eq!(
            (e.n_id_target, e.n_id_other, e.n_unid, e.n_id_claimant),
            (1, 1, 0, 0)
        );
        assert_eq!((e.hit_target, e.n_target), (3, 4));
        assert_eq!((e.hit_claimant, e.n_claimant), (1, 4));
        assert!(e.scored());
        assert!((e.dfrac() - (0.75 - 0.25)).abs() < 1e-12);
        // Matched intensity shares: scan 100: 100/105; scan 101: 100/200.
        assert!((f64::from(e.int_target) - (100.0 / 105.0 + 0.5)).abs() < 1e-5);
        assert!((f64::from(e.int_claimant) - 0.5).abs() < 1e-5);
    }

    #[test]
    fn no_distinguishing_ions_is_not_scored() {
        let a = Ions::new(1, &[frag(300.0), frag(400.0)]);
        let b = Ions::new(2, &[frag(300.0), frag(400.0), frag(500.0)]);
        let index = Ms2RunIndex::new(vec![scan(100.0, &[(300.0, 5.0)])]);
        let q = Query {
            mz: 500.0,
            t0: 100.0,
            w: 2.0,
            faims: None,
            target: &a,
            claimants: &[&b],
        };
        let e = evidence(&index, &ScanIdentity::default(), &q, Tolerance::Ppm(20.0));
        assert_eq!(e.n_scans, 1);
        assert_eq!(e.n_unid, 1);
        assert!(!e.scored());
        assert_eq!(e.dfrac(), 0.0);
    }

    const MODEL: &str = "# test model\nfeatures: co logq\nthreshold: 0.5\ntree 0\nS 0 0.5 1 2\nL -2.0\nL 2.0\ntree 1\nL 0.5\n";

    #[test]
    fn the_model_follows_its_trees() {
        let m = Model::parse(MODEL).unwrap();
        let claim = Claim {
            claimant: Some(3),
            co_apex: true,
            n: 1,
            ..Claim::default()
        };
        let row = Row {
            claim: &claim,
            shared: None,
            ms2_supported: false,
            transfer_q: 0.001,
            llr: 3.0,
            width: 5.0,
        };
        let p = m.probability(&row);
        assert!((p - 1.0 / (1.0 + (-2.5f64).exp())).abs() < 1e-12);
        assert!(m.is_shared(&row));
        let alone = Row {
            claim: &Claim::default(),
            ..row
        };
        assert!(
            !m.is_shared(&alone),
            "a transfer without a claimant is never shared"
        );
        assert!((m.probability(&alone) - 1.0 / (1.0 + 1.5f64.exp())).abs() < 1e-12);
        assert_eq!(m.provenance, ["test model"]);
        assert!(!m.needs_ms2_evidence());
    }

    #[test]
    fn the_shipped_model_is_well_formed() {
        let m = Model::shipped();
        assert!(m.threshold > 0.0 && m.threshold < 1.0, "{}", m.threshold);
        assert!(m.trees.len() > 1);
        assert!(m.needs_ms2_evidence());
        assert!(m.provenance.iter().any(|l| l.contains("trained")));
        assert!(m.provenance.iter().any(|l| l.contains("held-out")));
        // A transfer with no claimant at all is never downgraded, whatever
        // else the row says.
        let row = Row {
            claim: &Claim::default(),
            shared: None,
            ms2_supported: false,
            transfer_q: 0.0001,
            llr: 100.0,
            width: 30.0,
        };
        assert!(!m.is_shared(&row));
        let p = m.probability(&row);
        assert!((0.0..=1.0).contains(&p));
    }

    #[test]
    fn malformed_models_are_rejected() {
        assert!(Model::parse("features: co\ntree\nL 1.0\n")
            .unwrap_err()
            .contains("threshold"));
        assert!(Model::parse("features: nope\nthreshold: 0.5\n")
            .unwrap_err()
            .contains("unknown feature"));
        assert!(
            Model::parse("features: co\nthreshold: 0.5\ntree\nS 3 0.5 1 2\nL 1\nL 2\n")
                .unwrap_err()
                .contains("malformed")
        );
    }
}
