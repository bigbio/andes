//! Label specifications: sites, channels and their atom deltas, the named labels and presets,
//! and the parsers for label tokens and compositions.

use super::fit::ChannelEnvelope;
use super::{D_H, HIGHEST_OFFSET, MAX_CHANNELS};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Site {
    Residue(u8),
    NTerm,
}

impl Site {
    pub fn describe(self) -> String {
        match self {
            Site::Residue(r) => (r as char).to_string(),
            Site::NTerm => "N-term".into(),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct LabelSite {
    pub site: Site,
    pub delta: f64,
    pub deuterium: u32,
    pub heavy_atoms: u32,
    pub carbon13: u32,
    pub natural_c: i32,
    pub natural_h: i32,
    pub name: String,
    pub unimod: Option<u32>,
    pub chemical: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Channel {
    pub name: String,
    pub qpx_label: Option<String>,
    pub sites: Vec<LabelSite>,
}

impl Channel {
    pub fn label_at(&self, site: Site) -> Option<&LabelSite> {
        self.sites.iter().find(|l| l.site == site)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct LabelSpec {
    pub channels: Vec<Channel>,
    pub reference: usize,
    pub enrichment_2h_metabolic: f64,
    pub enrichment_2h_chemical: f64,
    pub enrichment_heavy: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LabelAtoms {
    pub deuterium: u32,
    pub p_2h: f64,
    pub carbon13: u32,
    pub nitrogen15: u32,
    pub p_heavy: f64,
    pub natural_c: i32,
    pub natural_h: i32,
}

pub(super) struct NamedLabel {
    token: &'static str,
    pub(super) sites: &'static [Site],
    delta: f64,
    pub(super) deuterium: u32,
    pub(super) heavy_atoms: u32,
    pub(super) carbon13: u32,
    pub(super) natural_c: i32,
    pub(super) natural_h: i32,
    pub(super) name: &'static str,
    unimod: u32,
    pub(super) chemical: bool,
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

pub const PRESETS: &[(&str, &str)] = &[
    ("silac2", "[][Lys8,Arg10]"),
    ("silac3", "[][Lys4,Arg6][Lys8,Arg10]"),
    ("dimethyl2", "[Dimethyl0][Dimethyl4]"),
    ("dimethyl3", "[Dimethyl0][Dimethyl4][Dimethyl8]"),
];

impl LabelSpec {
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
        let mass = |c: &[LabelSite]| c.iter().map(|l| l.delta).sum::<f64>();
        if channels.windows(2).any(|w| mass(&w[1]) <= mass(&w[0])) {
            return Err("channels must be given in increasing mass order".into());
        }
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
                let qpx_label = if all_silac && n <= 3 {
                    Some(
                        match (n, i) {
                            (_, 0) => "SILAC light",
                            (3, 1) => "SILAC medium",
                            _ => "SILAC heavy",
                        }
                        .to_string(),
                    )
                } else if all_silac {
                    None
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

    pub fn len(&self) -> usize {
        self.channels.len()
    }

    pub fn is_empty(&self) -> bool {
        self.channels.is_empty()
    }

    pub fn sites(&self) -> Vec<Site> {
        let mut out: Vec<Site> = Vec::new();
        for l in self.channels.iter().flatten_sites() {
            if !out.contains(&l.site) {
                out.push(l.site);
            }
        }
        out
    }

    pub fn search_labels(&self) -> Vec<(usize, &LabelSite)> {
        self.channels
            .iter()
            .enumerate()
            .flat_map(|(i, c)| c.sites.iter().map(move |l| (i, l)))
            .collect()
    }

    pub fn channel_of_label(&self, site: Site, name: &str) -> Option<usize> {
        self.channels
            .iter()
            .position(|c| c.label_at(site).is_some_and(|l| l.name == name))
    }

    pub fn has_unlabelled_channel(&self) -> bool {
        self.channels[0].sites.is_empty()
    }

    pub fn has_deuterium(&self) -> bool {
        self.channels
            .iter()
            .flat_map(|c| &c.sites)
            .any(|l| l.deuterium > 0)
    }

    pub fn channel_delta(&self, c: usize, sites: &[Site], counts: &[u32]) -> f64 {
        sites
            .iter()
            .zip(counts)
            .map(|(s, &n)| self.channels[c].label_at(*s).map_or(0.0, |l| l.delta) * n as f64)
            .sum()
    }

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

pub(super) trait FlattenSites<'a> {
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
    let (body, comp) = match token.split_once(':') {
        Some((b, a)) => (b, Some(parse_composition(token, a)?)),
        None => (token, None),
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
    let comp = comp.unwrap_or(Composition {
        chemical: true,
        c13: 0,
        n15: 0,
        d: 0,
    });
    let shift = comp.c13 as f64 * C13_C12 + comp.n15 as f64 * N15_N14 + comp.d as f64 * D_H;
    let (natural_c, natural_h) = if comp.chemical {
        let light = delta - shift;
        if light < 0.5 * CH2 {
            return Err(format!(
                "label {token:?}: a chemical label adds atoms, but {delta} Da less its heavy \
                 isotopes leaves {light:.4} Da"
            ));
        }
        let c = (light / CH2).round() as i32;
        if comp.c13 as i32 > c {
            return Err(format!(
                "label {token:?}: {} 13C, but the label adds about {c} C",
                comp.c13
            ));
        }
        (c, 2 * c)
    } else {
        if (shift - delta).abs() > 0.01 {
            return Err(format!(
                "label {token:?}: {}C {}N {}D add {shift:.4} Da, not {delta} (for a chemical \
                 label that adds atoms write ':chem…', e.g. Nterm+32.056407:chem4D)",
                comp.c13, comp.n15, comp.d
            ));
        }
        (0, 0)
    };
    Ok(vec![LabelSite {
        site,
        delta,
        deuterium: comp.d,
        heavy_atoms: comp.c13 + comp.n15,
        carbon13: comp.c13,
        natural_c,
        natural_h,
        name: format!("Label:{}:{delta:.6}", site.describe()),
        unimod: None,
        chemical: comp.chemical,
    }])
}

const C13_C12: f64 = model::mass::ISOTOPE;
const N15_N14: f64 = model::mass::N15_N14;
const CH2: f64 = model::mass::C + 2.0 * model::mass::H;

#[derive(Debug, Clone, Copy)]
pub(super) struct Composition {
    pub(super) chemical: bool,
    c13: u32,
    n15: u32,
    d: u32,
}

fn parse_composition(token: &str, text: &str) -> Result<Composition, String> {
    let bad = || {
        format!(
            "label {token:?}: expected ':' then [chem] and the heavy atoms as <n>C <n>N <n>D, \
             e.g. K+8.014199:6C2N or Nterm+32.056407:chem4D"
        )
    };
    let t = text.trim();
    let (chemical, mut rest) = match t.get(..4) {
        Some(p) if p.eq_ignore_ascii_case("chem") => (true, &t[4..]),
        _ => (false, t),
    };
    let mut comp = Composition {
        chemical,
        c13: 0,
        n15: 0,
        d: 0,
    };
    let mut seen = Vec::new();
    while !rest.is_empty() {
        let digits = rest.find(|ch: char| !ch.is_ascii_digit()).ok_or_else(bad)?;
        if digits == 0 {
            return Err(bad());
        }
        let n: u32 = rest[..digits].parse().map_err(|_| bad())?;
        let element = rest[digits..]
            .chars()
            .next()
            .ok_or_else(bad)?
            .to_ascii_uppercase();
        if seen.contains(&element) {
            return Err(bad());
        }
        seen.push(element);
        match element {
            'C' => comp.c13 = n,
            'N' => comp.n15 = n,
            'D' => comp.d = n,
            _ => return Err(bad()),
        }
        rest = &rest[digits + 1..];
    }
    if !chemical && seen.is_empty() {
        return Err(bad());
    }
    Ok(comp)
}
