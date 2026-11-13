//! MS1 labelled-pair quantification (SILAC, reductive dimethyl).
//!
//! A *label group* is one peptide identity (the peptidoform without its
//! channel labels) at one charge, measured in every channel of the label
//! specification at once. The channels of a group co-elute (13C/15N labels) or
//! elute with a small, learned offset (deuterium labels), so they are fitted
//! jointly: one Gaussian elution profile (apex `t0`, width `w`) shared by all
//! channels, one amplitude per channel.
//!
//! Model. For every isotope trace `t` read for the group (isotopes −1..3 of
//! every channel, merged when two read the same centroid):
//!
//! `y_t[j] = b_t + Σ_c M[t][c] · A_c · g(t_j; t0 + δ_c, w)`
//!
//! where `M[t][c]` is the part of channel `c`'s isotope pattern that falls on
//! trace `t` (see [`ChannelEnvelope`]): its own isotopes, its lower
//! isotopologues from incomplete label enrichment, and the isotopes of
//! *another* channel that land on the trace (the light M+4 lies about 9 ppm
//! below a Lys4 or Dimethyl4 channel's monoisotope, so the two are one
//! centroid). Interference from the group's own channels is therefore part of
//! the model, never a reason to drop a channel. With the backgrounds `b_t`
//! profiled out and each trace weighted by its noise `1 / σ_t²`, the
//! amplitudes for a fixed `(t0, w)` solve a small non-negative least-squares
//! problem, solved exactly over the subsets of free channels. The profile
//! log-likelihood ratio against "background only" is `LLR = A·b / 2` at the
//! optimum. `(t0, w)` maximise `LLR − (t0 − t_hat)² / (2 s²)` over the scan
//! RTs within `t_hat ± 3 s` and the widths `{0.7, 1, 1.4} · w_prior`, as in
//! [`crate::evidence`], then refined between scans (the vertex of the
//! parabola through the best apex and its neighbours). The noise `σ_t` of a trace is `1.4826 · median|Δy| /
//! √2` over the window's scans away from the expected peak (beyond
//! `3 · max(s, w)` of `t_hat`), floored at half the median detection limit.
//!
//! Per channel at the optimum: its amplitude and standard error, its own
//! likelihood ratio (the LLR lost when the channel is forced to 0), the
//! cosine between its isotope heights (other channels' modelled parts
//! removed) and its pattern, the coherence of its monoisotopic trace with the
//! profile, and the residual at its M−1 position (after its own lower
//! isotopologue is accounted for) relative to its monoisotope.
//!
//! Match between channels. A channel without an identification is measured
//! at its exact m/z with the identified channels anchoring the profile; its
//! decoys are the same fit with the channel moved by non-label mass shifts
//! ([`decoy_shifts`]), and target/decoy competition of the channel
//! scores ([`channel_score`]) gives the channel's quality q-value
//! ([`channel_qvalues`]). [`classify_channel`] turns the fit into a status.

use model::mass::ISOTOPE;

use crate::ms1_index::{MobilityGate, Ms1RunIndex};

/// Mass difference of deuterium and protium (2H − 1H), Da.
pub const D_H: f64 = 2.014_101_778 - 1.007_825_032;
/// Isotope offsets of [`ChannelEnvelope`]: `LOWEST_OFFSET ..= HIGHEST_OFFSET`
/// relative to the channel's monoisotope.
pub const LOWEST_OFFSET: i32 = -2;
/// Isotopes up to M+5 are modelled so that a lighter channel's M+4..M+5
/// lands on a heavier channel 4 Da up (Lys4, Dimethyl4).
pub const HIGHEST_OFFSET: i32 = 5;
const N_OFFSETS: usize = (HIGHEST_OFFSET - LOWEST_OFFSET + 1) as usize;
/// Isotopes read per channel: M−1 (the check for another envelope) to M+3.
const READ_OFFSETS: [i32; 5] = [-1, 0, 1, 2, 3];
/// Most channels a specification may have (the fit enumerates their subsets).
pub const MAX_CHANNELS: usize = 5;

/// Where a label sits.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Site {
    /// The side chain of a residue (one-letter code).
    Residue(u8),
    /// The peptide N-terminal amine (any residue).
    NTerm,
}

impl Site {
    /// `K`, `R`, ... or `N-term`.
    pub fn describe(self) -> String {
        match self {
            Site::Residue(r) => (r as char).to_string(),
            Site::NTerm => "N-term".into(),
        }
    }
}

/// One label of a channel: where it sits, its mass and its isotope
/// composition (for the enrichment model of [`ChannelEnvelope`]).
#[derive(Debug, Clone, PartialEq)]
pub struct LabelSite {
    pub site: Site,
    /// Mass added by the label, Da.
    pub delta: f64,
    /// Number of 2H atoms of the label.
    pub deuterium: u32,
    /// Number of 13C and 15N atoms of the label.
    pub heavy_atoms: u32,
    /// Of those, the 13C atoms (the rest are 15N).
    pub carbon13: u32,
    /// C and H atoms the label adds to the peptide, its own 13C and 2H
    /// included (a chemical label adds atoms, 2 C and 4 H per dimethyl site;
    /// a metabolic one substitutes them and adds none). The envelope's
    /// natural part is the peptide's atoms plus these, minus the label's 13C,
    /// 15N and 2H.
    pub natural_c: i32,
    pub natural_h: i32,
    /// Unimod name (the modification name the search uses).
    pub name: String,
    pub unimod: Option<u32>,
    /// A chemical label (dimethyl, ICPL) rather than a metabolic one (SILAC):
    /// selects the 2H enrichment.
    pub chemical: bool,
}

/// One channel: its labels (none for the unlabelled SILAC light channel).
#[derive(Debug, Clone, PartialEq)]
pub struct Channel {
    /// Short name (`light`, `Lys8+Arg10`, `Dimethyl4`).
    pub name: String,
    /// The SDRF/QPX label (`SILAC heavy`, `DIMETHYL4`); `None` when the
    /// channel is outside that vocabulary.
    pub qpx_label: Option<String>,
    pub sites: Vec<LabelSite>,
}

impl Channel {
    /// The label of `site` in this channel, if any.
    pub fn label_at(&self, site: Site) -> Option<&LabelSite> {
        self.sites.iter().find(|l| l.site == site)
    }
}

/// A label specification: the channels in increasing mass order, the
/// reference channel of the ratios and the label enrichments.
#[derive(Debug, Clone, PartialEq)]
pub struct LabelSpec {
    pub channels: Vec<Channel>,
    /// Index of the reference channel of the ratios.
    pub reference: usize,
    /// 2H enrichment of metabolic labels (Lys4, Leu3). OpenMS measured 0.96
    /// for Lys4 on PXD003327.
    pub enrichment_2h_metabolic: f64,
    /// 2H enrichment of chemical labels (dimethyl). 0.993–0.999 measured for
    /// the dimethyl channels of PXD000714; 0.99 by default.
    pub enrichment_2h_chemical: f64,
    /// 13C / 15N enrichment of heavy labels, per atom (99 % for the common
    /// SILAC amino acids).
    pub enrichment_heavy: f64,
}

/// The label atoms of one channel of one peptide ([`LabelSpec::channel_atoms`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LabelAtoms {
    pub deuterium: u32,
    /// 2H enrichment that applies to them.
    pub p_2h: f64,
    pub carbon13: u32,
    pub nitrogen15: u32,
    /// 13C / 15N enrichment.
    pub p_heavy: f64,
    /// Natural C and H atoms added (see [`LabelSite::natural_c`]).
    pub natural_c: i32,
    pub natural_h: i32,
}

/// A label the specification syntax knows by name.
struct NamedLabel {
    token: &'static str,
    sites: &'static [Site],
    delta: f64,
    deuterium: u32,
    heavy_atoms: u32,
    carbon13: u32,
    natural_c: i32,
    natural_h: i32,
    name: &'static str,
    unimod: u32,
    chemical: bool,
}

const LYS: &[Site] = &[Site::Residue(b'K')];
const ARG: &[Site] = &[Site::Residue(b'R')];
const LEU: &[Site] = &[Site::Residue(b'L')];
const AMINES: &[Site] = &[Site::Residue(b'K'), Site::NTerm];

const NAMED_LABELS: &[NamedLabel] = &[
    NamedLabel {
        token: "Lys4",
        sites: LYS,
        delta: 4.025107,
        deuterium: 4,
        heavy_atoms: 0,
        carbon13: 0,
        natural_c: 0,
        natural_h: 0,
        name: "Label:2H(4)",
        unimod: 481,
        chemical: false,
    },
    NamedLabel {
        token: "Lys6",
        sites: LYS,
        delta: 6.020129,
        deuterium: 0,
        heavy_atoms: 6,
        carbon13: 6,
        natural_c: 0,
        natural_h: 0,
        name: "Label:13C(6)",
        unimod: 188,
        chemical: false,
    },
    NamedLabel {
        token: "Lys8",
        sites: LYS,
        delta: 8.014199,
        deuterium: 0,
        heavy_atoms: 8,
        carbon13: 6,
        natural_c: 0,
        natural_h: 0,
        name: "Label:13C(6)15N(2)",
        unimod: 259,
        chemical: false,
    },
    NamedLabel {
        token: "Arg6",
        sites: ARG,
        delta: 6.020129,
        deuterium: 0,
        heavy_atoms: 6,
        carbon13: 6,
        natural_c: 0,
        natural_h: 0,
        name: "Label:13C(6)",
        unimod: 188,
        chemical: false,
    },
    NamedLabel {
        token: "Arg10",
        sites: ARG,
        delta: 10.008269,
        deuterium: 0,
        heavy_atoms: 10,
        carbon13: 6,
        natural_c: 0,
        natural_h: 0,
        name: "Label:13C(6)15N(4)",
        unimod: 267,
        chemical: false,
    },
    NamedLabel {
        token: "Leu3",
        sites: LEU,
        delta: 3.01883,
        deuterium: 3,
        heavy_atoms: 0,
        carbon13: 0,
        natural_c: 0,
        natural_h: 0,
        name: "Label:2H(3)",
        unimod: 262,
        chemical: false,
    },
    NamedLabel {
        token: "Dimethyl0",
        sites: AMINES,
        delta: 28.0313,
        deuterium: 0,
        heavy_atoms: 0,
        carbon13: 0,
        natural_c: 2,
        natural_h: 4,
        name: "Dimethyl",
        unimod: 36,
        chemical: true,
    },
    NamedLabel {
        token: "Dimethyl4",
        sites: AMINES,
        delta: 32.056407,
        deuterium: 4,
        heavy_atoms: 0,
        carbon13: 0,
        natural_c: 2,
        natural_h: 4,
        name: "Dimethyl:2H(4)",
        unimod: 199,
        chemical: true,
    },
    NamedLabel {
        token: "Dimethyl6",
        sites: AMINES,
        delta: 34.063117,
        deuterium: 4,
        heavy_atoms: 2,
        carbon13: 2,
        natural_c: 2,
        natural_h: 4,
        name: "Dimethyl:2H(4)13C(2)",
        unimod: 510,
        chemical: true,
    },
    NamedLabel {
        token: "Dimethyl8",
        sites: AMINES,
        delta: 36.07567,
        deuterium: 6,
        heavy_atoms: 2,
        carbon13: 2,
        natural_c: 2,
        natural_h: 4,
        name: "Dimethyl:2H(6)13C(2)",
        unimod: 330,
        chemical: true,
    },
];

/// The presets of `--label`.
pub const PRESETS: &[(&str, &str)] = &[
    ("silac2", "[][Lys8,Arg10]"),
    ("silac3", "[][Lys4,Arg6][Lys8,Arg10]"),
    ("dimethyl2", "[Dimethyl0][Dimethyl4]"),
    ("dimethyl3", "[Dimethyl0][Dimethyl4][Dimethyl8]"),
];

impl LabelSpec {
    /// Parses a preset (`silac2`, `silac3`, `dimethyl2`, `dimethyl3`) or a
    /// channel list in the syntax of OpenMS FeatureFinderMultiplex: one
    /// bracket per channel, in increasing mass order, each holding
    /// comma-separated labels, e.g. `[][Lys8,Arg10]`. A label is a known name
    /// (`Lys4`, `Lys6`, `Lys8`, `Arg6`, `Arg10`, `Leu3`, `Dimethyl0`,
    /// `Dimethyl4`, `Dimethyl6`, `Dimethyl8`) or a custom `<site>+<mass>`
    /// with an optional `:<n>D` for its 2H atoms (`K+8.014199`,
    /// `Nterm+32.056407:4D`).
    pub fn parse(text: &str) -> Result<Self, String> {
        let t = text.trim();
        let body = PRESETS
            .iter()
            .find(|(name, _)| name.eq_ignore_ascii_case(t))
            .map_or(t, |(_, spec)| spec);
        let mut channels = Vec::new();
        let mut rest = body.trim();
        while !rest.is_empty() {
            let Some(open) = rest.strip_prefix('[') else {
                return Err(format!(
                    "label specification {text:?}: expected '[' at {rest:?} (one bracket per \
                     channel, e.g. [][Lys8,Arg10], or a preset: {})",
                    PRESETS.iter().map(|p| p.0).collect::<Vec<_>>().join(", ")
                ));
            };
            let Some(close) = open.find(']') else {
                return Err(format!("label specification {text:?}: unclosed '['"));
            };
            let inner = &open[..close];
            rest = open[close + 1..].trim_start();
            let mut sites = Vec::new();
            for token in inner.split(',').map(str::trim).filter(|s| !s.is_empty()) {
                sites.extend(parse_label_token(token)?);
            }
            channels.push(sites);
        }
        Self::from_channels(channels)
    }

    fn from_channels(channels: Vec<Vec<LabelSite>>) -> Result<Self, String> {
        if channels.len() < 2 {
            return Err("a label specification needs at least two channels".into());
        }
        if channels.len() > MAX_CHANNELS {
            return Err(format!(
                "a label specification has at most {MAX_CHANNELS} channels, got {}",
                channels.len()
            ));
        }
        // Every channel labels the same sites (or none: the unlabelled one).
        let site_set = |c: &[LabelSite]| {
            let mut s: Vec<Site> = c.iter().map(|l| l.site).collect();
            s.sort_by_key(|s| match s {
                Site::Residue(r) => *r as i32,
                Site::NTerm => -1,
            });
            s
        };
        for c in &channels {
            let s = site_set(c);
            if s.windows(2).any(|w| w[0] == w[1]) {
                return Err("a channel labels one site twice".into());
            }
        }
        let labelled: Vec<Vec<Site>> = channels
            .iter()
            .filter(|c| !c.is_empty())
            .map(|c| site_set(c))
            .collect();
        if labelled.windows(2).any(|w| w[0] != w[1]) {
            return Err(
                "every labelled channel must label the same sites (e.g. K and R in each)".into(),
            );
        }
        if channels.iter().skip(1).any(|c| c.is_empty()) {
            return Err("only the first channel may be unlabelled".into());
        }
        // Increasing mass for one of each site.
        let mass = |c: &[LabelSite]| c.iter().map(|l| l.delta).sum::<f64>();
        if channels.windows(2).any(|w| mass(&w[1]) <= mass(&w[0])) {
            return Err("channels must be given in increasing mass order".into());
        }
        // One modification name names one (site, channel).
        for (i, a) in channels.iter().enumerate() {
            for (j, b) in channels.iter().enumerate().skip(i + 1) {
                for la in a {
                    if b.iter().any(|lb| lb.name == la.name) {
                        return Err(format!(
                            "label {} appears in channels {} and {}",
                            la.name,
                            i + 1,
                            j + 1
                        ));
                    }
                }
            }
        }
        let n = channels.len();
        let all_silac = channels
            .iter()
            .flatten()
            .all(|l| matches!(l.site, Site::Residue(b'K' | b'R')) && !l.chemical)
            && channels[0].is_empty();
        let all_dimethyl = channels.iter().all(|c| {
            !c.is_empty()
                && c.iter()
                    .all(|l| l.name.starts_with("Dimethyl") && l.chemical)
        });
        let channels = channels
            .into_iter()
            .enumerate()
            .map(|(i, sites)| {
                let name = if sites.is_empty() {
                    "light".to_string()
                } else {
                    let mut names: Vec<String> = Vec::new();
                    for l in &sites {
                        let token = NAMED_LABELS
                            .iter()
                            .find(|nl| {
                                nl.name == l.name
                                    && nl.sites.contains(&l.site)
                                    && (nl.delta - l.delta).abs() < 1e-4
                            })
                            .map_or_else(
                                || format!("{}+{:.4}", l.site.describe(), l.delta),
                                |nl| nl.token.to_string(),
                            );
                        if !names.contains(&token) {
                            names.push(token);
                        }
                    }
                    names.join("+")
                };
                let qpx_label = if all_silac {
                    Some(
                        match (n, i) {
                            (_, 0) => "SILAC light",
                            (3, 1) => "SILAC medium",
                            _ => "SILAC heavy",
                        }
                        .to_string(),
                    )
                } else if all_dimethyl {
                    let d = sites[0].delta;
                    NAMED_LABELS
                        .iter()
                        .find(|nl| nl.chemical && (nl.delta - d).abs() < 1e-4)
                        .map(|nl| nl.token.to_ascii_uppercase())
                } else {
                    None
                };
                Channel {
                    name,
                    qpx_label,
                    sites,
                }
            })
            .collect();
        Ok(LabelSpec {
            channels,
            reference: 0,
            enrichment_2h_metabolic: 0.96,
            enrichment_2h_chemical: 0.99,
            enrichment_heavy: 0.99,
        })
    }

    /// Number of channels.
    pub fn len(&self) -> usize {
        self.channels.len()
    }

    pub fn is_empty(&self) -> bool {
        self.channels.is_empty()
    }

    /// Every site some channel labels.
    pub fn sites(&self) -> Vec<Site> {
        let mut out: Vec<Site> = Vec::new();
        for l in self.channels.iter().flatten_sites() {
            if !out.contains(&l.site) {
                out.push(l.site);
            }
        }
        out
    }

    /// The labels as variable modifications of the search: `(channel index,
    /// label)` for every label of every channel. The unlabelled channel has
    /// none: it is the search's light form.
    pub fn search_labels(&self) -> Vec<(usize, &LabelSite)> {
        self.channels
            .iter()
            .enumerate()
            .flat_map(|(i, c)| c.sites.iter().map(move |l| (i, l)))
            .collect()
    }

    /// The channel whose label at `site` has modification name `name`.
    pub fn channel_of_label(&self, site: Site, name: &str) -> Option<usize> {
        self.channels
            .iter()
            .position(|c| c.label_at(site).is_some_and(|l| l.name == name))
    }

    /// Whether the first channel is the unlabelled one (SILAC): a labelable
    /// site without a label is then that channel's.
    pub fn has_unlabelled_channel(&self) -> bool {
        self.channels[0].sites.is_empty()
    }

    /// Whether some channel carries deuterium (its RT may shift).
    pub fn has_deuterium(&self) -> bool {
        self.channels
            .iter()
            .flat_map(|c| &c.sites)
            .any(|l| l.deuterium > 0)
    }

    /// Mass added to the unlabelled peptide by channel `c`, for `counts[k]`
    /// labelled positions of site `sites[k]` (see [`Self::sites`]).
    pub fn channel_delta(&self, c: usize, sites: &[Site], counts: &[u32]) -> f64 {
        sites
            .iter()
            .zip(counts)
            .map(|(s, &n)| self.channels[c].label_at(*s).map_or(0.0, |l| l.delta) * n as f64)
            .sum()
    }

    /// The label atoms channel `c` puts on a peptide with `counts[k]`
    /// labelled positions of site `sites[k]`.
    pub fn channel_atoms(&self, c: usize, sites: &[Site], counts: &[u32]) -> LabelAtoms {
        let mut a = LabelAtoms {
            deuterium: 0,
            p_2h: self.enrichment_2h_metabolic,
            carbon13: 0,
            nitrogen15: 0,
            p_heavy: self.enrichment_heavy,
            natural_c: 0,
            natural_h: 0,
        };
        let mut chemical = false;
        for (s, &n) in sites.iter().zip(counts) {
            if let Some(l) = self.channels[c].label_at(*s) {
                a.deuterium += l.deuterium * n;
                a.carbon13 += l.carbon13 * n;
                a.nitrogen15 += (l.heavy_atoms - l.carbon13.min(l.heavy_atoms)) * n;
                a.natural_c += l.natural_c * n as i32;
                a.natural_h += l.natural_h * n as i32;
                chemical |= l.chemical && l.deuterium > 0;
            }
        }
        if chemical {
            a.p_2h = self.enrichment_2h_chemical;
        }
        a
    }

    /// The isotope pattern of channel `c` of a peptide whose unlabelled
    /// formula is `light`: the natural envelope of the atoms the labels do
    /// not fix (a Lys8 fixes 6 C and 2 N as 13C/15N, which then take no part
    /// in the natural 13C/15N distribution) plus the natural atoms a chemical
    /// label adds (two C per dimethyl site), shifted by the label, with the
    /// lower isotopologues of incomplete enrichment.
    pub fn channel_envelope(
        &self,
        c: usize,
        sites: &[Site],
        counts: &[u32],
        light: &crate::formula::Formula,
    ) -> ChannelEnvelope {
        let a = self.channel_atoms(c, sites, counts);
        let mut f = *light;
        let add = |x: u32, d: i32| (x as i64 + d as i64).max(0) as u32;
        f.c = add(f.c, a.natural_c).saturating_sub(a.carbon13);
        f.n = f.n.saturating_sub(a.nitrogen15);
        f.h = add(f.h, a.natural_h).saturating_sub(a.deuterium);
        let env = f.envelope((HIGHEST_OFFSET + 1) as usize);
        ChannelEnvelope::new(
            &env,
            a.deuterium,
            a.p_2h,
            a.carbon13 + a.nitrogen15,
            a.p_heavy,
        )
    }

    /// A one-line description.
    pub fn describe(&self) -> String {
        self.channels
            .iter()
            .enumerate()
            .map(|(i, c)| {
                let r = if i == self.reference {
                    " (reference)"
                } else {
                    ""
                };
                format!("{}: {}{r}", i + 1, c.name)
            })
            .collect::<Vec<_>>()
            .join(", ")
    }
}

trait FlattenSites<'a> {
    fn flatten_sites(self) -> Box<dyn Iterator<Item = &'a LabelSite> + 'a>;
}

impl<'a, I: Iterator<Item = &'a Channel> + 'a> FlattenSites<'a> for I {
    fn flatten_sites(self) -> Box<dyn Iterator<Item = &'a LabelSite> + 'a> {
        Box::new(self.flat_map(|c| c.sites.iter()))
    }
}

fn parse_label_token(token: &str) -> Result<Vec<LabelSite>, String> {
    if let Some(nl) = NAMED_LABELS
        .iter()
        .find(|nl| nl.token.eq_ignore_ascii_case(token))
    {
        return Ok(nl
            .sites
            .iter()
            .map(|&site| LabelSite {
                site,
                delta: nl.delta,
                deuterium: nl.deuterium,
                heavy_atoms: nl.heavy_atoms,
                carbon13: nl.carbon13,
                natural_c: nl.natural_c,
                natural_h: nl.natural_h,
                name: nl.name.to_string(),
                unimod: Some(nl.unimod),
                chemical: nl.chemical,
            })
            .collect());
    }
    let (body, d) = match token.split_once(':') {
        Some((b, d)) => {
            let n = d
                .trim()
                .trim_end_matches(['D', 'd'])
                .parse::<u32>()
                .map_err(|_| format!("label {token:?}: expected ':<n>D' for its 2H atoms"))?;
            (b, n)
        }
        None => (token, 0),
    };
    let Some((site, mass)) = body.split_once('+') else {
        return Err(format!(
            "unknown label {token:?}: use a name (Lys4, Lys6, Lys8, Arg6, Arg10, Leu3, \
             Dimethyl0, Dimethyl4, Dimethyl6, Dimethyl8) or <site>+<mass>, e.g. K+8.014199"
        ));
    };
    let delta: f64 = mass
        .trim()
        .parse()
        .map_err(|_| format!("label {token:?}: bad mass {mass:?}"))?;
    if !(delta > 0.0 && delta < 100.0) {
        return Err(format!("label {token:?}: mass must be in (0, 100) Da"));
    }
    let site = match site.trim() {
        s if s.eq_ignore_ascii_case("nterm") || s.eq_ignore_ascii_case("n-term") => Site::NTerm,
        s if s.len() == 1 && s.as_bytes()[0].is_ascii_uppercase() => Site::Residue(s.as_bytes()[0]),
        s => {
            return Err(format!(
                "label {token:?}: unknown site {s:?} (a residue or Nterm)"
            ))
        }
    };
    // Nominal-mass shift not carried by 2H is carried by 13C/15N.
    let heavy = ((delta - d as f64 * D_H) / 1.0).round().max(0.0) as u32;
    Ok(vec![LabelSite {
        site,
        delta,
        deuterium: d,
        heavy_atoms: heavy,
        carbon13: heavy,
        natural_c: 0,
        natural_h: 0,
        name: format!("Label:{}:{delta:.6}", site.describe()),
        unimod: None,
        chemical: false,
    }])
}

/// The isotope pattern of one channel of a peptide, relative to the
/// channel's monoisotopic m/z, offsets [`LOWEST_OFFSET`]..=[`HIGHEST_OFFSET`]
/// (index 0 is M−2), summing to 1.
///
/// It is the unlabelled peptide's natural envelope shifted by the label, not
/// an averagine of the labelled mass: a 13C/15N/2H label moves the pattern
/// without widening it. Incomplete enrichment adds lower isotopologues: each
/// of the `n` label atoms at enrichment `p` is light with probability `1 − p`,
/// which moves the molecule one isotope down (by `D_H` for 2H, about
/// `ISOTOPE` for 13C, within any MS1 tolerance of it). With two 2H-labelled
/// Lys4 at p = 0.96, the M−1 peak is about 0.3 of the monoisotope.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChannelEnvelope(pub [f64; N_OFFSETS]);

impl ChannelEnvelope {
    /// `light`: the unlabelled peptide's envelope from its monoisotope (any
    /// length; entries beyond M+5 are ignored). `deuterium` 2H atoms at
    /// enrichment `p_2h`, `heavy` 13C/15N atoms at `p_heavy`.
    pub fn new(light: &[f64], deuterium: u32, p_2h: f64, heavy: u32, p_heavy: f64) -> Self {
        // Lower-isotopologue distribution: shift 0, -1, -2.
        let lower = |n: u32, p: f64| -> [f64; 3] {
            if n == 0 || p >= 1.0 {
                return [1.0, 0.0, 0.0];
            }
            let q = (1.0 - p).clamp(0.0, 1.0);
            let nf = n as f64;
            let p0 = p.powf(nf);
            let p1 = nf * q * p.powf(nf - 1.0);
            let p2 = if n >= 2 {
                nf * (nf - 1.0) / 2.0 * q * q * p.powf(nf - 2.0)
            } else {
                0.0
            };
            [p0, p1, p2]
        };
        let a = lower(deuterium, p_2h);
        let b = lower(heavy, p_heavy);
        let mut shift = [0.0; 3];
        for (i, x) in a.iter().enumerate() {
            for (j, y) in b.iter().enumerate() {
                if i + j < 3 {
                    shift[i + j] += x * y;
                }
            }
        }
        let mut out = [0.0; N_OFFSETS];
        for (k, &v) in light.iter().take((HIGHEST_OFFSET + 1) as usize).enumerate() {
            for (j, &s) in shift.iter().enumerate() {
                // offset k - j, index (k - j) - LOWEST_OFFSET
                let idx = k as i32 - j as i32 - LOWEST_OFFSET;
                out[idx as usize] += v * s;
            }
        }
        let total: f64 = out.iter().sum();
        if total > 0.0 {
            for v in &mut out {
                *v /= total;
            }
        }
        ChannelEnvelope(out)
    }

    /// Relative abundance at isotope `offset` (0 = monoisotope).
    pub fn at(&self, offset: i32) -> f64 {
        if (LOWEST_OFFSET..=HIGHEST_OFFSET).contains(&offset) {
            self.0[(offset - LOWEST_OFFSET) as usize]
        } else {
            0.0
        }
    }
}

/// One group hypothesis to fit.
#[derive(Debug, Clone)]
pub struct GroupProbe {
    pub charge: u8,
    /// Monoisotopic m/z of every channel.
    pub channel_mz: Vec<f64>,
    pub envelopes: Vec<ChannelEnvelope>,
    /// Expected apex (seconds) of the reference profile and its prior width.
    pub t_hat: f64,
    pub s: f64,
    /// Expected elution width (Gaussian sigma), seconds.
    pub w_prior: f64,
    /// RT offset of each channel's elution from the shared profile, seconds
    /// (0 for 13C/15N labels; see [`learn_rt_shifts`]).
    pub rt_shift: Vec<f64>,
    /// Channels whose amplitude is fitted; the others are held at 0 (used to
    /// measure one channel alone).
    pub free: Vec<bool>,
    pub mobility: Option<MobilityGate>,
}

/// Fit settings.
#[derive(Debug, Clone, Copy)]
pub struct FitParams {
    /// MS1 tolerance of the traces, ppm.
    pub tol_ppm: f64,
    /// Width grid, multiples of `w_prior`.
    pub widths: &'static [f64],
}

pub const DEFAULT_WIDTHS: &[f64] = &[0.7, 1.0, 1.4];
pub const FREE_WIDTHS: &[f64] = &[0.5, 0.6, 0.7, 0.85, 1.0, 1.2, 1.4, 1.7, 2.0];

impl Default for FitParams {
    fn default() -> Self {
        Self {
            tol_ppm: 10.0,
            widths: DEFAULT_WIDTHS,
        }
    }
}

/// One channel at the optimum of [`fit_group`].
#[derive(Debug, Clone, PartialEq)]
pub struct ChannelFit {
    /// Amplitude of the whole isotope pattern at the apex.
    pub amplitude: f64,
    pub amplitude_se: f64,
    /// `amplitude · w · √(2π)`: the area of the fitted profile.
    pub quantity: f64,
    pub quantity_se: f64,
    /// Log-likelihood ratio lost when the channel is held at 0.
    pub llr: f64,
    /// Cosine between the channel's own isotope heights (other channels'
    /// modelled parts removed) and its pattern; 0 when either is zero.
    pub envelope_cos: f64,
    /// Pearson correlation of the channel's monoisotopic trace (others'
    /// modelled parts removed) with its profile over `t0 ± 2w`.
    pub coherence: f64,
    /// Residual at the channel's M−1 (after every modelled part, its own
    /// lower isotopologues included) over `π_0 · max(A, se)`.
    pub m1_ratio: f64,
    /// Intensity-weighted ppm error of the monoisotopic centroids within
    /// `w` of the channel's apex; `None` when no centroid fell within the
    /// tolerance there (nothing was measured at the channel's m/z).
    pub mass_err_ppm: Option<f64>,
    /// Elution centroid of the channel's own signal (monoisotope and M+1,
    /// the other channels' modelled parts removed, positive values over
    /// `t0 ± 2w` of the channel's profile) minus the profile's apex, seconds;
    /// 0 when it has no signal there. A real channel elutes with its partner
    /// channels; another species' envelope at its position, at a random
    /// offset.
    pub apex_offset: f64,
    /// Shape-free quantity: the channel's own isotope traces (M..M+3, the
    /// other channels' modelled parts and the fitted background removed)
    /// integrated over `t0 − 3w … t0 + 4w` (trapezoids), divided by the
    /// pattern share they hold. Not clamped (an absent channel averages 0).
    pub trace_quantity: f64,
}

/// The joint fit of a group.
#[derive(Debug, Clone, PartialEq)]
pub struct GroupFit {
    pub t0: f64,
    pub w: f64,
    /// Total log-likelihood ratio of the channels against background only.
    pub llr: f64,
    /// `(t0 − t_hat) / s`.
    pub rt_z: f64,
    pub n_points: u32,
    pub channels: Vec<ChannelFit>,
}

/// The traces of one group: unique m/z values and the mixing matrix.
struct Traces {
    mz: Vec<f64>,
    /// `mix[t * n_ch + c]`.
    mix: Vec<f64>,
    /// Per channel, the trace of each read offset (`READ_OFFSETS`).
    own: Vec<[usize; READ_OFFSETS.len()]>,
}

fn build_traces(probe: &GroupProbe, tol_ppm: f64) -> Traces {
    let z = probe.charge.max(1) as f64;
    let n_ch = probe.channel_mz.len();
    let mut mz: Vec<f64> = Vec::new();
    let tol = |m: f64| m * tol_ppm * 1e-6;
    let mut own = vec![[0usize; READ_OFFSETS.len()]; n_ch];
    for (c, &mc) in probe.channel_mz.iter().enumerate() {
        for (k, &off) in READ_OFFSETS.iter().enumerate() {
            let m = mc + off as f64 * ISOTOPE / z;
            let t = match mz.iter().position(|&x| (x - m).abs() <= tol(m)) {
                Some(t) => t,
                None => {
                    mz.push(m);
                    mz.len() - 1
                }
            };
            own[c][k] = t;
        }
    }
    let n_t = mz.len();
    let mut mix = vec![0.0; n_t * n_ch];
    for (t, &mt) in mz.iter().enumerate() {
        for (c, &mc) in probe.channel_mz.iter().enumerate() {
            for off in LOWEST_OFFSET..=HIGHEST_OFFSET {
                let m = mc + off as f64 * ISOTOPE / z;
                if (m - mt).abs() <= tol(mt) {
                    mix[t * n_ch + c] += probe.envelopes[c].at(off);
                }
            }
        }
    }
    Traces { mz, mix, own }
}

/// Solves `max_A A·b − AᵀHA/2` over `A ≥ 0` with `A_c = 0` where `!free[c]`,
/// exactly, by enumerating the free subsets. Returns `(A, LLR, passive set)`.
fn nn_quadratic(h: &[f64], b: &[f64], n: usize, free: &[bool]) -> (Vec<f64>, f64, Vec<bool>) {
    let mut best = (vec![0.0; n], 0.0, vec![false; n]);
    for mask in 1u32..(1 << n) {
        if (0..n).any(|c| mask & (1 << c) != 0 && !free[c]) {
            continue;
        }
        let idx: Vec<usize> = (0..n).filter(|c| mask & (1 << c) != 0).collect();
        let k = idx.len();
        let mut m = vec![0.0; k * k];
        let mut v = vec![0.0; k];
        for (i, &ci) in idx.iter().enumerate() {
            v[i] = b[ci];
            for (j, &cj) in idx.iter().enumerate() {
                m[i * k + j] = h[ci * n + cj];
            }
        }
        let Some(x) = solve(&mut m, &mut v, k) else {
            continue;
        };
        if x.iter().any(|&a| a < 0.0) {
            continue;
        }
        let llr: f64 = idx.iter().zip(&x).map(|(&c, &a)| a * b[c]).sum::<f64>() / 2.0;
        if llr > best.1 {
            let mut a = vec![0.0; n];
            let mut p = vec![false; n];
            for (&c, &xa) in idx.iter().zip(&x) {
                a[c] = xa;
                p[c] = true;
            }
            best = (a, llr, p);
        }
    }
    best
}

/// Gaussian elimination with partial pivoting; `None` when singular.
fn solve(m: &mut [f64], v: &mut [f64], k: usize) -> Option<Vec<f64>> {
    for col in 0..k {
        let piv =
            (col..k).max_by(|&a, &b| m[a * k + col].abs().total_cmp(&m[b * k + col].abs()))?;
        if m[piv * k + col].abs() < 1e-300 {
            return None;
        }
        if piv != col {
            for j in 0..k {
                m.swap(piv * k + j, col * k + j);
            }
            v.swap(piv, col);
        }
        for r in col + 1..k {
            let f = m[r * k + col] / m[col * k + col];
            for j in col..k {
                m[r * k + j] -= f * m[col * k + j];
            }
            v[r] -= f * v[col];
        }
    }
    let mut x = vec![0.0; k];
    for r in (0..k).rev() {
        let mut s = v[r];
        for j in r + 1..k {
            s -= m[r * k + j] * x[j];
        }
        x[r] = s / m[r * k + r];
        if !x[r].is_finite() {
            return None;
        }
    }
    Some(x)
}

/// Inverse of a small symmetric positive matrix (the free block of `H`).
fn invert(h: &[f64], k: usize) -> Option<Vec<f64>> {
    let mut inv = vec![0.0; k * k];
    for col in 0..k {
        let mut m = h.to_vec();
        let mut e = vec![0.0; k];
        e[col] = 1.0;
        let x = solve(&mut m, &mut e, k)?;
        for r in 0..k {
            inv[r * k + col] = x[r];
        }
    }
    Some(inv)
}

fn robust_sigma(y: &[f64], floor: f64) -> f64 {
    if y.len() < 2 {
        return floor.max(1.0);
    }
    let mut d: Vec<f64> = y.windows(2).map(|w| (w[1] - w[0]).abs()).collect();
    let mid = d.len() / 2;
    let med = *d.select_nth_unstable_by(mid, f64::total_cmp).1;
    (1.4826 * med / std::f64::consts::SQRT_2)
        .max(floor)
        .max(1e-9)
}

/// Fits `probe` on `ms1`. `None` when the window holds fewer than 3 scans.
pub fn fit_group(ms1: &Ms1RunIndex, probe: &GroupProbe, params: &FitParams) -> Option<GroupFit> {
    let n_ch = probe.channel_mz.len();
    if n_ch == 0 || n_ch > MAX_CHANNELS {
        return None;
    }
    let s = probe.s.max(1e-3);
    let w_prior = probe.w_prior.max(1e-3);
    let max_shift = probe.rt_shift.iter().fold(0.0f64, |a, &d| a.max(d.abs()));
    let w_max = w_prior * params.widths.iter().fold(1.0f64, |a, &b| a.max(b));
    // Enough baseline on both sides for the noise estimate whatever `s`.
    let half = (4.0 * s + 4.0 * w_max).max(8.0 * w_max) + max_shift;
    let window = ms1.scan_window(probe.t_hat - half, probe.t_hat + half);
    if window.len() < 3 {
        return None;
    }
    let rt: Vec<f64> = window.clone().map(|i| ms1.rt(i)).collect();
    let traces = build_traces(probe, params.tol_ppm);
    let n_t = traces.mz.len();
    let read = |m: f64| -> Vec<f64> {
        let tol = m * params.tol_ppm * 1e-6;
        window
            .clone()
            .map(|i| ms1.slice_in(i, m - tol, m + tol, probe.mobility).1 as f64)
            .collect()
    };
    let y: Vec<Vec<f64>> = traces.mz.iter().map(|&m| read(m)).collect();
    let mut limits: Vec<f64> = window
        .clone()
        .map(|i| ms1.detection_limit(i) as f64)
        .collect();
    let mid = limits.len() / 2;
    let limit = *limits.select_nth_unstable_by(mid, f64::total_cmp).1;
    // Noise from the scans away from the expected peak (all of them when
    // too few are), so that the peak's own slopes do not inflate it.
    let away: Vec<bool> = rt
        .iter()
        .map(|&t| (t - probe.t_hat).abs() > 3.0 * s.max(w_max) + max_shift)
        .collect();
    let n_away = away.iter().filter(|&&a| a).count();
    let sigma: Vec<f64> = y
        .iter()
        .map(|t| {
            if n_away >= 10 {
                let base: Vec<f64> = t
                    .iter()
                    .zip(&away)
                    .filter(|(_, &a)| a)
                    .map(|(v, _)| *v)
                    .collect();
                robust_sigma(&base, 0.5 * limit)
            } else {
                robust_sigma(t, 0.5 * limit)
            }
        })
        .collect();
    let w_inv: Vec<f64> = sigma.iter().map(|s| 1.0 / (s * s)).collect();
    // The amplitudes are fitted over the peak's own neighbourhood only
    // (`t_hat ± (3 s + 3 w)`): the wide window serves the noise estimate, and
    // other species eluting in it would otherwise pull on the background.
    let reach = 3.0 * s + 3.0 * w_max + max_shift;
    let j0 = rt.partition_point(|&t| t < probe.t_hat - reach);
    let j1 = rt.partition_point(|&t| t <= probe.t_hat + reach);
    let (rt, y, window) = if j1 >= j0 + 5 {
        (
            rt[j0..j1].to_vec(),
            y.iter().map(|t| t[j0..j1].to_vec()).collect::<Vec<_>>(),
            window.start + j0..window.start + j1,
        )
    } else {
        (rt, y, window)
    };
    let n_j = rt.len();
    let ybar: Vec<f64> = y
        .iter()
        .map(|t| t.iter().sum::<f64>() / n_j as f64)
        .collect();
    let yc: Vec<Vec<f64>> = y
        .iter()
        .zip(&ybar)
        .map(|(t, m)| t.iter().map(|v| v - m).collect())
        .collect();

    // Per (t0, w): centred profiles per channel and the normal equations.
    let profiles = |t0: f64, w: f64| -> Vec<Vec<f64>> {
        (0..n_ch)
            .map(|c| {
                let tc = t0 + probe.rt_shift.get(c).copied().unwrap_or(0.0);
                let g: Vec<f64> = rt
                    .iter()
                    .map(|&t| (-(t - tc) * (t - tc) / (2.0 * w * w)).exp())
                    .collect();
                let gbar = g.iter().sum::<f64>() / n_j as f64;
                g.into_iter().map(|v| v - gbar).collect()
            })
            .collect()
    };
    let normal =
        |gc: &[Vec<f64>], w_inv: &[f64]| -> (Vec<f64>, Vec<f64>, Vec<Vec<f64>>, Vec<f64>) {
            // proj[t][c] = Σ_j yc_t g_c ; sgg[c][c'] = Σ_j g_c g_c'
            let mut sgg = vec![0.0; n_ch * n_ch];
            for a in 0..n_ch {
                for b in a..n_ch {
                    let v: f64 = gc[a].iter().zip(&gc[b]).map(|(x, y)| x * y).sum();
                    sgg[a * n_ch + b] = v;
                    sgg[b * n_ch + a] = v;
                }
            }
            let proj: Vec<Vec<f64>> = yc
                .iter()
                .map(|t| {
                    (0..n_ch)
                        .map(|c| t.iter().zip(&gc[c]).map(|(x, y)| x * y).sum())
                        .collect()
                })
                .collect();
            let mut h = vec![0.0; n_ch * n_ch];
            let mut b = vec![0.0; n_ch];
            for t in 0..n_t {
                let mt = &traces.mix[t * n_ch..(t + 1) * n_ch];
                for a in 0..n_ch {
                    if mt[a] == 0.0 {
                        continue;
                    }
                    b[a] += w_inv[t] * mt[a] * proj[t][a];
                    for c in 0..n_ch {
                        h[a * n_ch + c] += w_inv[t] * mt[a] * mt[c] * sgg[a * n_ch + c];
                    }
                }
            }
            (h, b, proj, sgg)
        };

    let lo = probe.t_hat - 3.0 * s;
    let hi = probe.t_hat + 3.0 * s;
    let mut t0s: Vec<f64> = rt.iter().copied().filter(|&t| t >= lo && t <= hi).collect();
    if t0s.is_empty() {
        let k = rt
            .iter()
            .enumerate()
            .min_by(|a, b| {
                (a.1 - probe.t_hat)
                    .abs()
                    .total_cmp(&(b.1 - probe.t_hat).abs())
            })
            .map(|x| x.0)?;
        t0s.push(rt[k]);
    }
    // The apex and width that maximise the likelihood ratio under the trace
    // weights `wt` (with the RT prior), refined to a sub-scan apex by the
    // vertex of the parabola through the best grid apex and its neighbours.
    let search_apex = |wt: &[f64]| -> Option<(f64, f64)> {
        let mut best: Option<(f64, f64, f64)> = None; // (objective, t0, w)
        for &t0 in &t0s {
            for &wm in params.widths {
                let w = wm * w_prior;
                let gc = profiles(t0, w);
                if gc
                    .iter()
                    .all(|g| g.iter().map(|v| v * v).sum::<f64>() < 1e-9)
                {
                    continue;
                }
                let (h, b, _, _) = normal(&gc, wt);
                let (_, llr, _) = nn_quadratic(&h, &b, n_ch, &probe.free);
                let obj = llr - (t0 - probe.t_hat).powi(2) / (2.0 * s * s);
                if best.is_none_or(|(o, _, _)| obj > o) {
                    best = Some((obj, t0, w));
                }
            }
        }
        let (best_obj, t0_grid, w) = best?;
        let objective = |t0: f64| -> f64 {
            let gc = profiles(t0, w);
            let (h, b, _, _) = normal(&gc, wt);
            let (_, llr, _) = nn_quadratic(&h, &b, n_ch, &probe.free);
            llr - (t0 - probe.t_hat).powi(2) / (2.0 * s * s)
        };
        let mut t0 = t0_grid;
        if let Some(k) = t0s.iter().position(|&x| x == t0_grid) {
            if k > 0 && k + 1 < t0s.len() {
                let (x0, x1, x2) = (t0s[k - 1], t0_grid, t0s[k + 1]);
                let (y0, y1, y2) = (objective(x0), best_obj, objective(x2));
                let den = (x0 - x1) * (x0 - x2) * (x1 - x2);
                if den.abs() > 0.0 {
                    let a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / den;
                    let b = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / den;
                    if a < 0.0 {
                        let v = (-b / (2.0 * a)).clamp(x0, x2);
                        if v.is_finite() && objective(v) > best_obj {
                            t0 = v;
                        }
                    }
                }
            }
        }
        Some((t0, w))
    };
    let (mut t0, mut w) = search_apex(&w_inv)?;
    // Interference: a trace whose residual at the optimum is far above its
    // noise carries another species (an interfered isotope of one channel);
    // its weight is divided by its excess (robust reweighting, twice), so
    // the amplitudes rest on the traces that follow the model. An interferer
    // strong enough to have pulled the apex to itself loses its weight here,
    // so the apex is searched again under the robust weights (and the weights
    // re-derived there), until the apex stays put.
    let mut w_fit = w_inv.clone();
    let mut gc = profiles(t0, w);
    let (mut h, mut b, mut proj, mut sgg) = normal(&gc, &w_fit);
    let (mut a, mut llr, mut passive) = nn_quadratic(&h, &b, n_ch, &probe.free);
    for round in 0..3 {
        let mut reweighted = false;
        for _ in 0..2 {
            let mut changed = false;
            for t in 0..n_t {
                let mt = &traces.mix[t * n_ch..(t + 1) * n_ch];
                if mt.iter().all(|&m| m == 0.0) {
                    continue;
                }
                let mut ss = 0.0;
                let mut n = 0usize;
                for (j, &tj) in rt.iter().enumerate() {
                    if (tj - t0).abs() > 2.0 * w {
                        continue;
                    }
                    let model: f64 = (0..n_ch).map(|c| mt[c] * a[c] * gc[c][j]).sum();
                    let r = yc[t][j] - model;
                    ss += r * r;
                    n += 1;
                }
                if n < 3 {
                    continue;
                }
                let chi = ss / n as f64 * w_inv[t];
                let target = if chi > 9.0 {
                    w_inv[t] * 9.0 / chi
                } else {
                    w_inv[t]
                };
                if (target - w_fit[t]).abs() > 1e-12 * w_inv[t] {
                    w_fit[t] = target;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
            reweighted = true;
            (h, b, proj, sgg) = normal(&gc, &w_fit);
            (a, llr, passive) = nn_quadratic(&h, &b, n_ch, &probe.free);
        }
        if !reweighted || round == 2 {
            break;
        }
        let Some((t1, w1)) = search_apex(&w_fit) else {
            break;
        };
        if t1 == t0 && w1 == w {
            break;
        }
        (t0, w) = (t1, w1);
        // Weights are re-derived at the new apex, from the noise alone.
        w_fit.clone_from(&w_inv);
        gc = profiles(t0, w);
        (h, b, proj, sgg) = normal(&gc, &w_fit);
        (a, llr, passive) = nn_quadratic(&h, &b, n_ch, &probe.free);
    }
    let area = w * (2.0 * std::f64::consts::PI).sqrt();

    // Standard errors: inverse of the passive block; 1/sqrt(H_cc) otherwise.
    let pidx: Vec<usize> = (0..n_ch).filter(|&c| passive[c]).collect();
    let mut se = vec![0.0; n_ch];
    let k = pidx.len();
    let mut block = vec![0.0; k * k];
    for (i, &ci) in pidx.iter().enumerate() {
        for (j, &cj) in pidx.iter().enumerate() {
            block[i * k + j] = h[ci * n_ch + cj];
        }
    }
    let inv = if k > 0 { invert(&block, k) } else { None };
    for c in 0..n_ch {
        se[c] = match (inv.as_ref(), pidx.iter().position(|&p| p == c)) {
            (Some(inv), Some(i)) => inv[i * k + i].max(0.0).sqrt(),
            _ => 1.0 / h[c * n_ch + c].max(1e-300).sqrt(),
        };
    }

    let n_points = rt.iter().filter(|&&t| (t - t0).abs() <= 2.0 * w).count() as u32;
    let channels = (0..n_ch)
        .map(|c| {
            // LLR lost when c is held at 0.
            let mut free = probe.free.clone();
            free[c] = false;
            let (_, llr0, _) = nn_quadratic(&h, &b, n_ch, &free);
            let llr_c = (llr - llr0).max(0.0);
            // Own isotope heights, the others' modelled parts removed.
            let s_cc = sgg[c * n_ch + c].max(1e-300);
            let height = |t: usize| -> f64 {
                let mut p = proj[t][c];
                for c2 in 0..n_ch {
                    if c2 != c {
                        p -= traces.mix[t * n_ch + c2] * a[c2] * sgg[c * n_ch + c2];
                    }
                }
                p / s_cc
            };
            let mut obs = Vec::new();
            let mut exp = Vec::new();
            for (k, &off) in READ_OFFSETS.iter().enumerate() {
                if off < 0 {
                    continue;
                }
                let t = traces.own[c][k];
                obs.push(height(t));
                exp.push(traces.mix[t * n_ch + c]);
            }
            let envelope_cos = cosine(&obs, &exp);
            // Coherence of the monoisotopic trace with the profile, around
            // the channel's own (shifted) apex.
            let t_ch = t0 + probe.rt_shift.get(c).copied().unwrap_or(0.0);
            let mono = traces.own[c][1];
            let mut xs = Vec::new();
            let mut gs = Vec::new();
            for (j, &t) in rt.iter().enumerate() {
                if (t - t_ch).abs() > 2.0 * w {
                    continue;
                }
                let mut v = yc[mono][j];
                for c2 in 0..n_ch {
                    if c2 != c {
                        v -= traces.mix[mono * n_ch + c2] * a[c2] * gc[c2][j];
                    }
                }
                xs.push(v);
                gs.push(gc[c][j]);
            }
            let coherence = if xs.len() >= 3 {
                pearson(&xs, &gs)
            } else {
                0.0
            };
            // Residual at M−1 after every modelled part.
            let m1 = traces.own[c][0];
            let mut r = proj[m1][c];
            for c2 in 0..n_ch {
                r -= traces.mix[m1 * n_ch + c2] * a[c2] * sgg[c * n_ch + c2];
            }
            let pi0 = probe.envelopes[c].at(0).max(1e-9);
            let m1_ratio = (r / s_cc).max(0.0) / (pi0 * a[c].max(se[c]).max(1e-300));
            // Elution centroid of the channel's own monoisotope and M+1.
            let tc = t0 + probe.rt_shift.get(c).copied().unwrap_or(0.0);
            let (mut sy, mut syt) = (0.0, 0.0);
            for k in [1usize, 2] {
                let tr = traces.own[c][k];
                for (j, &t) in rt.iter().enumerate() {
                    if (t - tc).abs() > 2.0 * w {
                        continue;
                    }
                    let mut v = y[tr][j] - ybar[tr];
                    for c2 in 0..n_ch {
                        if c2 != c {
                            v -= traces.mix[tr * n_ch + c2] * a[c2] * gc[c2][j];
                        }
                    }
                    if v > 0.0 {
                        sy += v;
                        syt += v * (t - tc);
                    }
                }
            }
            let apex_offset = if sy > 0.0 { syt / sy } else { 0.0 };
            // Shape-free quantity over t0 - 3w .. t0 + 4w.
            let trace_quantity = {
                let (lo, hi) = (tc - 3.0 * w, tc + 4.0 * w);
                let mut area = 0.0;
                let mut share = 0.0;
                for k in 1..READ_OFFSETS.len() {
                    let tr = traces.own[c][k];
                    let pi = traces.mix[tr * n_ch + c];
                    if pi <= 0.0 {
                        continue;
                    }
                    share += pi;
                    let v: Vec<f64> = (0..rt.len())
                        .map(|j| {
                            let mut x = yc[tr][j];
                            for c2 in 0..n_ch {
                                if c2 != c {
                                    x -= traces.mix[tr * n_ch + c2] * a[c2] * gc[c2][j];
                                }
                            }
                            x
                        })
                        .collect();
                    // Background: the trace's mean outside the window.
                    let out: Vec<f64> = (0..rt.len())
                        .filter(|&j| rt[j] < lo || rt[j] > hi)
                        .map(|j| v[j])
                        .collect();
                    let b = if out.is_empty() {
                        0.0
                    } else {
                        out.iter().sum::<f64>() / out.len() as f64
                    };
                    for j in 1..rt.len() {
                        if rt[j - 1] >= lo && rt[j] <= hi {
                            area += 0.5 * (v[j - 1] + v[j] - 2.0 * b) * (rt[j] - rt[j - 1]);
                        }
                    }
                }
                // The traces hold `share` of the channel's pattern; scale to
                // the whole pattern.
                if share > 0.0 {
                    area / share
                } else {
                    0.0
                }
            };
            // Mass error of the monoisotopic centroids near the apex.
            let mz = probe.channel_mz[c];
            let tol = mz * params.tol_ppm * 1e-6;
            let (mut sw, mut swe) = (0.0, 0.0);
            let t_ch = t0 + probe.rt_shift.get(c).copied().unwrap_or(0.0);
            for (j, i) in window.clone().enumerate() {
                if (rt[j] - t_ch).abs() > w {
                    continue;
                }
                if let Some((m, inten)) = ms1.apex_centroid(i, mz - tol, mz + tol) {
                    sw += inten as f64;
                    swe += inten as f64 * (m - mz) / mz * 1e6;
                }
            }
            ChannelFit {
                amplitude: a[c],
                amplitude_se: se[c],
                quantity: a[c] * area,
                quantity_se: se[c] * area,
                llr: llr_c,
                envelope_cos,
                coherence,
                m1_ratio,
                mass_err_ppm: (sw > 0.0).then(|| swe / sw),
                apex_offset,
                trace_quantity,
            }
        })
        .collect();
    Some(GroupFit {
        t0,
        w,
        llr,
        rt_z: (t0 - probe.t_hat) / s,
        n_points,
        channels,
    })
}

fn cosine(a: &[f64], b: &[f64]) -> f64 {
    let dot: f64 = a.iter().zip(b).map(|(x, y)| x * y).sum();
    let na: f64 = a.iter().map(|x| x * x).sum::<f64>().sqrt();
    let nb: f64 = b.iter().map(|x| x * x).sum::<f64>().sqrt();
    if na <= 0.0 || nb <= 0.0 {
        0.0
    } else {
        (dot / (na * nb)).clamp(-1.0, 1.0)
    }
}

fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    let ma = a.iter().sum::<f64>() / n;
    let mb = b.iter().sum::<f64>() / n;
    let (mut sab, mut saa, mut sbb) = (0.0, 0.0, 0.0);
    for (x, y) in a.iter().zip(b) {
        sab += (x - ma) * (y - mb);
        saa += (x - ma) * (x - ma);
        sbb += (y - mb) * (y - mb);
    }
    if saa <= 0.0 || sbb <= 0.0 {
        0.0
    } else {
        sab / (saa * sbb).sqrt()
    }
}

/// Nominal decoy shifts (Da) tried by [`decoy_shifts`], scaled by
/// [`crate::decoys::PEPTIDE_MASS_SLOPE`] to stay in the peptide mass-defect
/// band.
const DECOY_UNITS: [f64; 12] = [
    7.0, -7.0, 11.0, -11.0, 13.0, -13.0, 17.0, -17.0, 19.0, -19.0, 23.0, -23.0,
];

/// Neutral mass shifts (Da, before division by the charge) of the decoys of a
/// group whose channels `moved` are moved, for channel masses `deltas`
/// (relative to the first channel). A shift is admissible when no moved
/// channel lands within `min_sep` Da of an isotope (M−3..M+6) of any channel
/// at its true position: a decoy must not read the group's own signal (a
/// heavy channel moved by −7 Da sits on the light M+1 of a Lys8 pair within
/// a few ppm). `seed` rotates the admissible shifts; the first `n` are
/// returned (fewer when fewer are admissible).
pub fn decoy_shifts(
    moved: &[usize],
    deltas: &[f64],
    min_sep: f64,
    seed: u64,
    n: usize,
) -> Vec<f64> {
    let ok: Vec<f64> = DECOY_UNITS
        .iter()
        .map(|u| u * crate::decoys::PEPTIDE_MASS_SLOPE)
        .filter(|&s| {
            moved.iter().all(|&m| {
                deltas.iter().all(|&d| {
                    (-3..=6).all(|k| (deltas[m] + s - d - k as f64 * ISOTOPE).abs() >= min_sep)
                })
            })
        })
        .collect();
    let k = ok.len().max(1);
    (0..n.min(ok.len()))
        .map(|i| ok[(seed as usize).wrapping_add(i) % k])
        .collect()
}

/// The decoy separation of [`decoy_shifts`] for a peptide of neutral mass
/// `mass` at an MS1 tolerance of `tol_ppm`: three tolerances.
pub fn decoy_min_sep(mass: f64, tol_ppm: f64) -> f64 {
    3.0 * tol_ppm * 1e-6 * mass
}

/// The mass-accuracy factor of a measured channel: `exp(−½ (Δ / σ)²)` for a
/// mass error `err_ppm` `Δ` away from where the channel is expected
/// (`expected_ppm`, the run's calibration or the group's identified
/// channels), at the run's mass-error spread `sigma_ppm`. A decoy that lands
/// on another species' envelope sits anywhere in the tolerance window; a real
/// channel sits where its identified channels do.
pub fn mass_factor(err_ppm: Option<f64>, expected_ppm: f64, sigma_ppm: f64) -> f64 {
    // No centroid at the channel's m/z near its apex: no evidence for it.
    let Some(err_ppm) = err_ppm else {
        return 0.0;
    };
    let d = (err_ppm - expected_ppm) / sigma_ppm.max(1e-3);
    (-0.5 * d * d).exp()
}

/// The co-elution factor of a measured channel: `exp(−½ (offset / σ)²)` for
/// its [`ChannelFit::apex_offset`] against the run's spread `sigma_s` of the
/// identified channels' offsets.
pub fn coelution_factor(offset: f64, sigma_s: f64) -> f64 {
    let d = offset / sigma_s.max(1e-3);
    (-0.5 * d * d).exp()
}

/// The competition score of a channel (or a group's strongest channel):
/// `√(2·LLR)` scaled by its envelope cosine and its coherence (both clamped
/// at 0), so that a decoy that lands on another species' isotope, with the
/// wrong pattern or elution, scores low whatever its intensity.
pub fn channel_score(fit: &ChannelFit) -> f64 {
    (2.0 * fit.llr.max(0.0)).sqrt() * fit.envelope_cos.max(0.0) * fit.coherence.max(0.0)
}

/// Target/decoy q-values of scores: `rows` are `(score, is_decoy)`, higher is
/// better, with `decoys_per_target` decoys drawn for every target. The
/// q-value of each row, in input order, is the running minimum from the worst
/// rank of `((D + 1) / k) / max(T, 1)`, clamped to 1; decoys get theirs too.
/// Ties rank decoys first.
pub fn channel_qvalues(rows: &[(f64, bool)], decoys_per_target: usize) -> Vec<f64> {
    let k = decoys_per_target.max(1) as f64;
    let mut order: Vec<usize> = (0..rows.len()).collect();
    let key = |s: f64| if s.is_nan() { f64::NEG_INFINITY } else { s };
    order.sort_by(|&a, &b| {
        key(rows[b].0)
            .total_cmp(&key(rows[a].0))
            .then(rows[b].1.cmp(&rows[a].1))
    });
    let mut q = vec![1.0; rows.len()];
    let (mut t, mut d) = (0usize, 0usize);
    let mut raw = Vec::with_capacity(rows.len());
    for &i in &order {
        if rows[i].1 {
            d += 1;
        } else {
            t += 1;
        }
        raw.push((((d + 1) as f64 / k) / t.max(1) as f64).min(1.0));
    }
    let mut run = 1.0f64;
    for (r, &i) in order.iter().enumerate().rev() {
        run = run.min(raw[r]);
        q[i] = run;
    }
    q
}

/// Target/decoy q-values when targets have different numbers of decoys:
/// `rows` are `(score, decoy_weight)`, a target with weight 0 and each of a
/// target's `k` decoys with weight `1 / k`, so the weighted decoy count above
/// a score estimates the false targets there. The q-value of each row, in
/// input order, is the running minimum from the worst rank of `(D_w + w̄) /
/// max(T, 1)` (`w̄` the mean decoy weight: the "+1" of one decoy), clamped to
/// 1. Ties rank decoys first.
pub fn weighted_qvalues(rows: &[(f64, f64)]) -> Vec<f64> {
    let decoys: Vec<f64> = rows.iter().filter(|r| r.1 > 0.0).map(|r| r.1).collect();
    let w_bar = if decoys.is_empty() {
        1.0
    } else {
        decoys.iter().sum::<f64>() / decoys.len() as f64
    };
    let key = |s: f64| if s.is_nan() { f64::NEG_INFINITY } else { s };
    let mut order: Vec<usize> = (0..rows.len()).collect();
    order.sort_by(|&a, &b| {
        key(rows[b].0)
            .total_cmp(&key(rows[a].0))
            .then((rows[b].1 > 0.0).cmp(&(rows[a].1 > 0.0)))
    });
    let (mut t, mut d) = (0usize, 0.0f64);
    let mut raw = Vec::with_capacity(rows.len());
    for &i in &order {
        if rows[i].1 > 0.0 {
            d += rows[i].1;
        } else {
            t += 1;
        }
        raw.push(((d + w_bar) / t.max(1) as f64).min(1.0));
    }
    let mut q = vec![1.0; rows.len()];
    let mut run = 1.0f64;
    for (r, &i) in order.iter().enumerate().rev() {
        run = run.min(raw[r]);
        q[i] = run;
    }
    q
}

/// What a channel of a group reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChannelStatus {
    /// A PSM of the run identifies the group in this channel.
    Identified,
    /// Not identified, measured at its position (match between channels or
    /// a transfer between runs) at the quality q-value cut.
    Matched,
    /// The group was not identified in this run; it was transferred from
    /// another run (match between runs) at the transfer q-value cut and this
    /// channel shows its pattern there.
    Transferred,
    /// No signal: the quantity is an upper bound, not a measurement.
    BelowDetection,
    /// Something is there that is not this channel (another identified
    /// peptide, a signal that does not co-elute or does not show the
    /// channel's pattern, or a quality q-value above the cut).
    NotQuantifiable,
}

impl ChannelStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            ChannelStatus::Identified => "identified",
            ChannelStatus::Matched => "match_channel",
            ChannelStatus::Transferred => "transfer",
            ChannelStatus::BelowDetection => "below_detection",
            ChannelStatus::NotQuantifiable => "not_quantifiable",
        }
    }

    /// Whether the quantity is a measurement of this channel.
    pub fn is_quantified(self) -> bool {
        matches!(
            self,
            ChannelStatus::Identified | ChannelStatus::Matched | ChannelStatus::Transferred
        )
    }
}

/// Status thresholds of [`classify_channel`].
#[derive(Debug, Clone, Copy)]
pub struct ClassifyParams {
    /// Quality q-value cut of a matched channel.
    pub fdr: f64,
    /// Per-channel LLR below which a channel shows no signal.
    pub min_signal_llr: f64,
    /// Minimum envelope cosine of a matched channel.
    pub min_cosine: f64,
    /// A channel weaker than this fraction of the strongest identified one…
    pub weak_fraction: f64,
    /// …must correlate with its profile by at least this much.
    pub min_weak_coherence: f64,
    /// Residual M−1 above this fraction of the monoisotope marks another
    /// envelope (one isotope lower) under the channel.
    pub max_m1_ratio: f64,
}

impl Default for ClassifyParams {
    fn default() -> Self {
        Self {
            fdr: 0.01,
            min_signal_llr: 4.0,
            min_cosine: 0.7,
            weak_fraction: 0.2,
            min_weak_coherence: 0.5,
            max_m1_ratio: 0.3,
        }
    }
}

/// The status of a channel without an identification. `q` is its quality
/// q-value ([`channel_qvalues`]), `anchor_amplitude` the largest amplitude of
/// the group's identified channels, `claimed` whether another identified
/// peptide of the run sits on the channel's position.
pub fn classify_channel(
    fit: &ChannelFit,
    q: f64,
    anchor_amplitude: f64,
    claimed: bool,
    p: &ClassifyParams,
) -> ChannelStatus {
    if claimed {
        return ChannelStatus::NotQuantifiable;
    }
    if fit.amplitude <= 0.0 || fit.llr < p.min_signal_llr {
        return ChannelStatus::BelowDetection;
    }
    let weak = fit.amplitude < p.weak_fraction * anchor_amplitude;
    if q > p.fdr
        || fit.envelope_cos < p.min_cosine
        || fit.m1_ratio > p.max_m1_ratio
        || (weak && fit.coherence < p.min_weak_coherence)
    {
        return ChannelStatus::NotQuantifiable;
    }
    ChannelStatus::Matched
}

/// The status of an identified channel: identified, unless its own pattern
/// fails the quality checks of [`classify_channel`] (envelope cosine, a
/// strong residual M−1 of another envelope), when its quantity is not this
/// channel's alone and it is not quantifiable. Its identification stays.
pub fn classify_identified(fit: &ChannelFit, p: &ClassifyParams) -> ChannelStatus {
    if fit.amplitude > 0.0 && (fit.envelope_cos < p.min_cosine || fit.m1_ratio > p.max_m1_ratio) {
        ChannelStatus::NotQuantifiable
    } else {
        ChannelStatus::Identified
    }
}

/// The status of a channel of a group transferred from another run (its
/// transfer accepted): [`ChannelStatus::Transferred`] when the channel shows
/// its signal and pattern, else below detection or not quantifiable as in
/// [`classify_channel`].
pub fn classify_transferred(fit: &ChannelFit, strongest: f64, p: &ClassifyParams) -> ChannelStatus {
    match classify_channel(fit, 0.0, strongest, false, p) {
        ChannelStatus::Matched => ChannelStatus::Transferred,
        other => other,
    }
}

/// `log2(A_c / A_ref)` with its standard error by the delta method; `None`
/// when either amplitude is not positive.
pub fn log2_ratio(c: &ChannelFit, reference: &ChannelFit) -> Option<(f64, f64)> {
    if c.amplitude <= 0.0 || reference.amplitude <= 0.0 {
        return None;
    }
    let r = (c.amplitude / reference.amplitude).log2();
    let rel = (c.amplitude_se / c.amplitude).powi(2)
        + (reference.amplitude_se / reference.amplitude).powi(2);
    Some((r, rel.sqrt() / std::f64::consts::LN_2))
}

/// Per-channel RT offsets (seconds) of a run from `pairs`, each `(channel,
/// apex of the channel fitted alone − apex of the reference fitted alone)`
/// for groups where both are strong: the median per channel when it has at
/// least `min_pairs`, else 0. The reference's offset is 0.
pub fn learn_rt_shifts(pairs: &[(usize, f64)], n_channels: usize, min_pairs: usize) -> Vec<f64> {
    (0..n_channels)
        .map(|c| {
            let mut v: Vec<f64> = pairs
                .iter()
                .filter(|p| p.0 == c && p.1.is_finite())
                .map(|p| p.1)
                .collect();
            if v.len() < min_pairs.max(1) {
                return 0.0;
            }
            v.sort_by(f64::total_cmp);
            let n = v.len();
            if n % 2 == 1 {
                v[n / 2]
            } else {
                (v[n / 2 - 1] + v[n / 2]) / 2.0
            }
        })
        .collect()
}

/// The median of `v` (NaN when empty).
pub fn median(mut v: Vec<f64>) -> f64 {
    v.retain(|x| x.is_finite());
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(f64::total_cmp);
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ms1_index::Ms1Scan;

    #[test]
    fn presets_parse() {
        let s = LabelSpec::parse("silac2").unwrap();
        assert_eq!(s.len(), 2);
        assert!(s.has_unlabelled_channel());
        assert_eq!(s.channels[1].name, "Lys8+Arg10");
        assert_eq!(s.channels[1].qpx_label.as_deref(), Some("SILAC heavy"));
        let s = LabelSpec::parse("silac3").unwrap();
        assert_eq!(s.channels[1].qpx_label.as_deref(), Some("SILAC medium"));
        assert!(s.has_deuterium());
        let s = LabelSpec::parse("dimethyl3").unwrap();
        assert_eq!(s.len(), 3);
        assert!(!s.has_unlabelled_channel());
        assert_eq!(s.channels[2].qpx_label.as_deref(), Some("DIMETHYL8"));
        assert_eq!(s.sites().len(), 2);
        let s = LabelSpec::parse("[][K+8.014199,R+10.008269]").unwrap();
        assert_eq!(s.channels[1].sites.len(), 2);
        assert_eq!(s.channels[1].sites[0].heavy_atoms, 8);
        assert!(s.channels[1].qpx_label.is_some());
        // A stacked mod's name joins its parts with '+': a custom label's
        // own name must not contain it, or the label is lost when split.
        for l in &s.channels[1].sites {
            assert!(!l.name.contains('+'), "{}", l.name);
            assert_eq!(s.channel_of_label(l.site, &l.name), Some(1));
        }
    }

    #[test]
    fn bad_specs_are_refused() {
        assert!(LabelSpec::parse("[Lys8]").is_err());
        assert!(LabelSpec::parse("[Lys8][]").is_err());
        assert!(LabelSpec::parse("[][Lys8][Arg10]").is_err());
        assert!(LabelSpec::parse("[][Lys8,Arg10][Lys4,Arg6]").is_err());
        assert!(LabelSpec::parse("[][Foo]").is_err());
    }

    #[test]
    fn channel_delta_counts_sites() {
        let s = LabelSpec::parse("silac2").unwrap();
        let sites = s.sites();
        let k = sites
            .iter()
            .position(|x| *x == Site::Residue(b'K'))
            .unwrap();
        let mut counts = vec![0; sites.len()];
        counts[k] = 2;
        assert!((s.channel_delta(1, &sites, &counts) - 16.028398).abs() < 1e-6);
        assert_eq!(s.channel_delta(0, &sites, &counts), 0.0);
    }

    #[test]
    fn deuterium_enrichment_gives_an_m_minus_1() {
        let light = [0.6, 0.3, 0.08, 0.02];
        let e = ChannelEnvelope::new(&light, 8, 0.96, 0, 1.0);
        // 8 x 0.04 / 0.96 ≈ 0.33 of the monoisotope one isotope down.
        let r = e.at(-1) / e.at(0);
        assert!(r > 0.2 && r < 0.6, "{r}");
        let total: f64 = e.0.iter().sum();
        assert!((total - 1.0).abs() < 1e-12);
        let plain = ChannelEnvelope::new(&light, 0, 1.0, 0, 1.0);
        assert_eq!(plain.at(-1), 0.0);
        assert!((plain.at(0) - 0.6).abs() < 1e-12);
    }

    /// A synthetic run: Gaussian peaks of the given (m/z, envelope, amplitude,
    /// apex) channels at charge 2, 1 s scans, plus a flat noise floor.
    pub(crate) fn synthetic_run(
        species: &[(f64, ChannelEnvelope, f64, f64)],
        w: f64,
        noise: f32,
    ) -> Ms1RunIndex {
        let scans = (0..240)
            .map(|i| {
                let rt = i as f64;
                let mut peaks: Vec<(f64, f32)> = Vec::new();
                for &(mz, env, amp, apex) in species {
                    let g = (-(rt - apex).powi(2) / (2.0 * w * w)).exp();
                    for off in LOWEST_OFFSET..=HIGHEST_OFFSET {
                        let v = amp * env.at(off) * g;
                        if v > noise as f64 {
                            peaks.push((mz + off as f64 * ISOTOPE / 2.0, v as f32));
                        }
                    }
                }
                // Background centroids everywhere else.
                for k in 0..50 {
                    peaks.push((300.0 + k as f64 * 17.3, noise));
                }
                peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
                // Peaks closer than the resolution are one centroid (summed).
                let mut merged: Vec<(f64, f32)> = Vec::new();
                for p in peaks {
                    match merged.last_mut() {
                        Some(l) if (p.0 - l.0).abs() < l.0 * 8e-6 => l.1 += p.1,
                        _ => merged.push(p),
                    }
                }
                Ms1Scan {
                    rt,
                    peaks: merged,
                    id: format!("scan={i}"),
                    ..Default::default()
                }
            })
            .collect();
        Ms1RunIndex::new(scans)
    }

    fn probe(mzs: &[f64], envs: &[ChannelEnvelope], t_hat: f64) -> GroupProbe {
        GroupProbe {
            charge: 2,
            channel_mz: mzs.to_vec(),
            envelopes: envs.to_vec(),
            t_hat,
            s: 3.0,
            w_prior: 4.0,
            rt_shift: vec![0.0; mzs.len()],
            free: vec![true; mzs.len()],
            mobility: None,
        }
    }

    #[test]
    fn duplex_ratio_is_recovered() {
        let light = [0.55, 0.3, 0.11, 0.04];
        let el = ChannelEnvelope::new(&light, 0, 1.0, 0, 1.0);
        let eh = ChannelEnvelope::new(&light, 0, 1.0, 8, 0.99);
        let ml = 600.3;
        let mh = ml + 8.014199 / 2.0;
        for ratio in [1.0, 0.25, 0.05] {
            let run = synthetic_run(
                &[(ml, el, 1e6, 120.0), (mh, eh, 1e6 * ratio, 120.0)],
                4.0,
                50.0,
            );
            let fit = fit_group(
                &run,
                &probe(&[ml, mh], &[el, eh], 121.0),
                &FitParams::default(),
            )
            .unwrap();
            let (r, se) = log2_ratio(&fit.channels[1], &fit.channels[0]).unwrap();
            assert!((r - ratio.log2()).abs() < 0.05, "ratio {ratio}: {r} ± {se}");
            assert!((fit.t0 - 120.0).abs() <= 1.0);
            assert!(fit.channels[1].envelope_cos > 0.95);
            assert!(
                fit.channels[1].m1_ratio < 0.1,
                "{}",
                fit.channels[1].m1_ratio
            );
        }
    }

    #[test]
    fn absent_channel_is_below_detection() {
        let light = [0.55, 0.3, 0.11, 0.04];
        let el = ChannelEnvelope::new(&light, 0, 1.0, 0, 1.0);
        let ml = 600.3;
        let mh = ml + 8.014199 / 2.0;
        let run = synthetic_run(&[(ml, el, 1e6, 120.0)], 4.0, 50.0);
        let fit = fit_group(
            &run,
            &probe(&[ml, mh], &[el, el], 121.0),
            &FitParams::default(),
        )
        .unwrap();
        let h = &fit.channels[1];
        assert!(h.amplitude < 1e3, "{}", h.amplitude);
        let status = classify_channel(
            h,
            0.0,
            fit.channels[0].amplitude,
            false,
            &ClassifyParams::default(),
        );
        assert_eq!(status, ChannelStatus::BelowDetection);
    }

    // Lys4 medium: the light M+4 is ~9 ppm below the medium monoisotope, so
    // the two share a centroid. The joint fit attributes it to both.
    #[test]
    fn overlapping_lys4_channel_is_unmixed() {
        // A large peptide, so that the light M+4 is substantial.
        let light = [0.25, 0.30, 0.22, 0.12, 0.07, 0.04];
        let el = ChannelEnvelope::new(&light, 0, 1.0, 0, 1.0);
        let em = ChannelEnvelope::new(&light, 4, 0.96, 0, 1.0);
        let z = 2.0;
        let ml = 1200.4;
        let mm = ml + 4.025107 / z;
        for ratio in [1.0, 0.2, 0.05] {
            let run = synthetic_run(
                &[(ml, el, 1e6, 100.0), (mm, em, 1e6 * ratio, 100.0)],
                4.0,
                50.0,
            );
            let fit = fit_group(
                &run,
                &probe(&[ml, mm], &[el, em], 100.0),
                &FitParams::default(),
            )
            .unwrap();
            let (r, _) = log2_ratio(&fit.channels[1], &fit.channels[0]).unwrap();
            assert!((r - ratio.log2()).abs() < 0.05, "ratio {ratio}: {r}");
            // The medium's own 2H lower isotopologue (and the light M+3) at
            // its M−1 are modelled, not taken for another envelope.
            assert!(
                fit.channels[1].m1_ratio < 0.3,
                "{}",
                fit.channels[1].m1_ratio
            );
            assert!(fit.channels[1].envelope_cos > 0.95);
        }
    }

    #[test]
    fn deuterium_shift_is_learned_and_fitted() {
        let light = [0.55, 0.3, 0.11, 0.04];
        let el = ChannelEnvelope::new(&light, 0, 1.0, 0, 1.0);
        let eh = ChannelEnvelope::new(&light, 6, 0.99, 2, 0.99);
        let ml = 700.0;
        let mh = ml + (36.07567 - 28.0313) / 2.0;
        let run = synthetic_run(&[(ml, el, 1e6, 120.0), (mh, eh, 5e5, 117.0)], 4.0, 50.0);
        // Each channel alone.
        let mut alone = probe(&[ml, mh], &[el, eh], 120.0);
        alone.free = vec![true, false];
        let a0 = fit_group(&run, &alone, &FitParams::default()).unwrap().t0;
        alone.free = vec![false, true];
        let a1 = fit_group(&run, &alone, &FitParams::default()).unwrap().t0;
        let shifts = learn_rt_shifts(&[(1, a1 - a0)], 2, 1);
        assert!((shifts[1] + 3.0).abs() <= 1.0, "{shifts:?}");
        let mut p = probe(&[ml, mh], &[el, eh], 120.0);
        p.rt_shift = shifts;
        let fit = fit_group(&run, &p, &FitParams::default()).unwrap();
        let (r, _) = log2_ratio(&fit.channels[1], &fit.channels[0]).unwrap();
        assert!((r + 1.0).abs() < 0.1, "{r}");
    }

    /// A single-trace interferer strong enough to win the unweighted apex
    /// search inside the search range: the robust weights take it out and
    /// the apex is searched again, back at the peptide.
    #[test]
    fn single_trace_interferer_does_not_move_the_apex() {
        let light = [0.55, 0.3, 0.11, 0.04];
        let el = ChannelEnvelope::new(&light, 0, 1.0, 0, 1.0);
        let eh = ChannelEnvelope::new(&light, 0, 1.0, 8, 0.99);
        let lone = ChannelEnvelope::new(&[1.0], 0, 1.0, 0, 1.0);
        let ml = 600.3;
        let mh = ml + 8.014199 / 2.0;
        // On the heavy channel's M+2 only, 7 s from the peptide's apex.
        let other = mh + 2.0 * ISOTOPE / 2.0;
        let run = synthetic_run(
            &[
                (ml, el, 1e6, 100.0),
                (mh, eh, 1e6, 100.0),
                (other, lone, 8e6, 107.0),
            ],
            4.0,
            50.0,
        );
        let fit = fit_group(
            &run,
            &probe(&[ml, mh], &[el, eh], 100.0),
            &FitParams::default(),
        )
        .unwrap();
        assert!((fit.t0 - 100.0).abs() < 1.0, "apex {}", fit.t0);
        let (r, _) = log2_ratio(&fit.channels[1], &fit.channels[0]).unwrap();
        assert!(r.abs() < 0.15, "{r}");
    }

    /// A channel with no centroid at its m/z near its apex has no mass
    /// evidence: it earns no mass credit, unlike a measured one at 0 ppm.
    #[test]
    fn unmeasured_channel_earns_no_mass_credit() {
        assert_eq!(mass_factor(None, 0.0, 2.0), 0.0);
        assert!((mass_factor(Some(0.0), 0.0, 2.0) - 1.0).abs() < 1e-12);
        assert!(mass_factor(Some(4.0), 0.0, 2.0) < 0.2);
    }

    #[test]
    fn interferer_at_channel_is_flagged() {
        let light = [0.55, 0.3, 0.11, 0.04];
        let el = ChannelEnvelope::new(&light, 0, 1.0, 0, 1.0);
        let ml = 600.3;
        let mh = ml + 8.014199 / 2.0;
        // Another species one isotope below the heavy, eluting elsewhere.
        let other = mh - ISOTOPE / 2.0;
        let run = synthetic_run(&[(ml, el, 1e6, 120.0), (other, el, 2e6, 121.0)], 4.0, 50.0);
        let fit = fit_group(
            &run,
            &probe(&[ml, mh], &[el, el], 121.0),
            &FitParams::default(),
        )
        .unwrap();
        let h = &fit.channels[1];
        // The other species' M+1 sits on the heavy monoisotope; its own
        // monoisotope sits on the heavy M−1.
        let status = classify_channel(
            h,
            0.0,
            fit.channels[0].amplitude,
            false,
            &ClassifyParams::default(),
        );
        assert_ne!(status, ChannelStatus::Matched, "{h:?}");
    }

    #[test]
    fn qvalues_rank_targets_over_decoys() {
        let rows = vec![
            (10.0, false),
            (9.0, false),
            (8.0, true),
            (7.0, false),
            (1.0, true),
        ];
        let q = channel_qvalues(&rows, 1);
        assert!((q[0] - 0.5).abs() < 1e-12);
        let q4 = channel_qvalues(&rows, 4);
        assert!((q4[0] - 0.125).abs() < 1e-12);
        // The weighted form with equal weights is the same estimate.
        let w: Vec<(f64, f64)> = rows
            .iter()
            .map(|&(s, d)| (s, if d { 0.25 } else { 0.0 }))
            .collect();
        let qw = weighted_qvalues(&w);
        for (a, b) in qw.iter().zip(&q4) {
            assert!((a - b).abs() < 1e-12);
        }
        assert!(q[3] >= q[0]);
        assert_eq!(q.len(), 5);
    }

    #[test]
    fn decoy_shifts_avoid_every_isotope() {
        // SILAC Lys8: the heavy channel moved by -7 Da would sit on the
        // light M+1 (8.0142 - 7.0035 = 1.0107 vs 1.0034).
        let deltas = [0.0, 8.014199];
        let min_sep = decoy_min_sep(1200.0, 10.0);
        let all = decoy_shifts(&[0, 1], &deltas, min_sep, 0, 12);
        assert!(!all.is_empty());
        for s in &all {
            for m in [0, 1] {
                for d in deltas {
                    for k in -3..=6 {
                        assert!(
                            (deltas[m] + s - d - k as f64 * ISOTOPE).abs() >= min_sep,
                            "{s}"
                        );
                    }
                }
            }
        }
        assert!(!all
            .iter()
            .any(|s| (s + 7.0 * crate::decoys::PEPTIDE_MASS_SLOPE).abs() < 1e-9));
        let heavy_only = decoy_shifts(&[1], &deltas, min_sep, 3, 4);
        assert_eq!(heavy_only.len(), 4);
    }
}
