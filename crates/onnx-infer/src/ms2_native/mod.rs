use std::path::Path;

use crate::rt_native::mods::N_COMP;
use crate::rt_native::{read_tensors, silu, Conv, RtError, RtPeptide, Tensor};

const PHOSPHO: u32 = 21;
const OXIDATION: u32 = 35;
const N_TOK: usize = 23;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NeutralLoss {
    None,
    H3po4,
    Ch4so,
}

/// Fragment ion series. `A`, `B` and `C` are N-terminal (prefix) fragments; `Y` and `Z`
/// (z-dot, `y - NH3 + H`) are C-terminal. `C` and `Z` come from electron-based activation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Series {
    A,
    B,
    Y,
    C,
    Z,
}

impl Series {
    /// True for the N-terminal series (a, b, c).
    pub const fn is_prefix(self) -> bool {
        matches!(self, Series::A | Series::B | Series::C)
    }

    /// True for the series that only electron-based activation (ETD, EThcD) produces.
    pub const fn is_electron(self) -> bool {
        matches!(self, Series::C | Series::Z)
    }
}

/// One output channel: ion series, fragment charge and neutral loss.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IonType {
    pub series: Series,
    pub charge: u8,
    pub loss: NeutralLoss,
}

impl IonType {
    pub const fn new(series: Series, charge: u8, loss: NeutralLoss) -> Self {
        Self {
            series,
            charge,
            loss,
        }
    }

    /// True for an N-terminal fragment; for the v1 channels this is the old `is_b` flag.
    pub const fn is_prefix(&self) -> bool {
        self.series.is_prefix()
    }
}

/// The 12 channels of an `AMS20001` (v1) model, per cleavage site.
pub const ION_TYPES: [IonType; 12] = [
    IonType::new(Series::B, 1, NeutralLoss::None),
    IonType::new(Series::B, 2, NeutralLoss::None),
    IonType::new(Series::Y, 1, NeutralLoss::None),
    IonType::new(Series::Y, 2, NeutralLoss::None),
    IonType::new(Series::B, 1, NeutralLoss::H3po4),
    IonType::new(Series::B, 2, NeutralLoss::H3po4),
    IonType::new(Series::Y, 1, NeutralLoss::H3po4),
    IonType::new(Series::Y, 2, NeutralLoss::H3po4),
    IonType::new(Series::B, 1, NeutralLoss::Ch4so),
    IonType::new(Series::B, 2, NeutralLoss::Ch4so),
    IonType::new(Series::Y, 1, NeutralLoss::Ch4so),
    IonType::new(Series::Y, 2, NeutralLoss::Ch4so),
];
pub const N_ION: usize = ION_TYPES.len();
pub const N_BASE_ION: usize = 4;

/// The channels of an `AMS20002` (v2) model: the 12 v1 channels in the same order, then
/// b/y 3+, a 1+, c 1+/2+ and z-dot 1+/2+ (`ms2_common.ION_TYPES_V2`).
pub const ION_TYPES_V2: [IonType; 19] = [
    ION_TYPES[0],
    ION_TYPES[1],
    ION_TYPES[2],
    ION_TYPES[3],
    ION_TYPES[4],
    ION_TYPES[5],
    ION_TYPES[6],
    ION_TYPES[7],
    ION_TYPES[8],
    ION_TYPES[9],
    ION_TYPES[10],
    ION_TYPES[11],
    IonType::new(Series::B, 3, NeutralLoss::None),
    IonType::new(Series::Y, 3, NeutralLoss::None),
    IonType::new(Series::A, 1, NeutralLoss::None),
    IonType::new(Series::C, 1, NeutralLoss::None),
    IonType::new(Series::C, 2, NeutralLoss::None),
    IonType::new(Series::Z, 1, NeutralLoss::None),
    IonType::new(Series::Z, 2, NeutralLoss::None),
];
pub const N_ION_V2: usize = ION_TYPES_V2.len();

/// Activation (dissociation + MS2 analyser); index 0 is unknown. The v2 entries come only
/// after the 6 v1 entries, so v1 indices never move.
pub const ACTIVATIONS: [&str; 9] = [
    "unknown",
    "HCD-FTMS",
    "HCD-ITMS",
    "CID-ITMS",
    "CID-FTMS",
    "HCD-ASTMS",
    "ETD",
    "EThcD",
    "HCD-TOF",
];
/// Instrument family; index 0 is unknown. The v2 entries come only after the 7 v1 entries.
pub const INSTRUMENTS: [&str; 9] = [
    "unknown",
    "QExactive",
    "Exploris",
    "Tribrid",
    "VelosElite",
    "Astral",
    "LTQ",
    "timsTOF",
    "TOF",
];
/// Isobaric label condition of a v2 model; index 0 is none.
pub const LABELS: [&str; 3] = ["none", "TMT", "iTRAQ"];
/// Activations whose spectra carry c and z-dot ions.
pub const ELECTRON_ACTIVATIONS: [&str; 2] = ["ETD", "EThcD"];
/// Rows of the activation and instrument embeddings in a v1 (`AMS20001`) model.
const V1_N_ACTIVATIONS: usize = 6;
const V1_N_INSTRUMENTS: usize = 7;
const MAGIC_V1: &[u8; 8] = b"AMS20001";
const MAGIC_V2: &[u8; 8] = b"AMS20002";
pub const MAX_CHARGE: usize = 6;

pub fn activation_index(name: &str) -> usize {
    ACTIVATIONS.iter().position(|a| *a == name).unwrap_or(0)
}

/// Index in [`LABELS`] for a label name (`none`, `TMT`, `TMTpro`, `iTRAQ8plex`, ...);
/// anything else is 0 (none).
pub fn label_index(name: &str) -> usize {
    let s = name.to_lowercase();
    if s.contains("itraq") {
        2
    } else if s.contains("tmt") {
        1
    } else {
        0
    }
}

pub fn instrument_class(name: &str) -> usize {
    let s: String = name
        .to_lowercase()
        .chars()
        .filter(|c| *c != ' ' && *c != '-')
        .collect();
    let has = |k: &str| s.contains(k);
    if s.is_empty() {
        0
    } else if has("astral") {
        5
    } else if has("exploris") {
        2
    } else if has("exactive") {
        1
    } else if ["fusion", "lumos", "eclipse", "ascend", "tribrid"]
        .iter()
        .any(|k| has(k))
    {
        3
    } else if has("velos") || has("elite") {
        4
    } else if has("ltq") {
        6
    } else if has("timstof") {
        7
    } else if has("tof") {
        8
    } else {
        0
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ms2Conditions {
    pub charge: u8,
    pub nce: Option<f32>,
    pub activation: usize,
    pub instrument: usize,
    /// Index in [`LABELS`] (0 none, 1 TMT, 2 iTRAQ); a v1 model ignores it.
    pub label: usize,
}

impl Ms2Conditions {
    pub fn new(
        charge: u8,
        nce: Option<f32>,
        activation: &str,
        instrument: &str,
        label: &str,
    ) -> Self {
        Self {
            charge,
            nce: nce.filter(|v| v.is_finite()),
            activation: activation_index(activation),
            instrument: instrument_class(instrument),
            label: label_index(label),
        }
    }

    fn numeric(&self) -> [f32; MAX_CHARGE + 2] {
        let mut v = [0.0f32; MAX_CHARGE + 2];
        v[(self.charge as usize).clamp(1, MAX_CHARGE) - 1] = 1.0;
        if let Some(nce) = self.nce {
            v[MAX_CHARGE] = nce / 100.0;
            v[MAX_CHARGE + 1] = 1.0;
        }
        v
    }
}

/// Valid v1 channels ([`ION_TYPES`]) per site, site-major.
pub fn valid_mask(p: &RtPeptide, charge: u8) -> Vec<bool> {
    channel_mask(p, charge, 0, &ION_TYPES)
}

/// Valid channels of `ions` per site, site-major. A channel is invalid when its fragment
/// charge exceeds the precursor charge, when it is a loss channel and the fragment holds no
/// residue that can lose that group, or when it is a c / z-dot channel and either
/// `activation` (an [`ACTIVATIONS`] index) is not electron-based or the cleavage is
/// N-terminal to proline.
pub fn channel_mask(p: &RtPeptide, charge: u8, activation: usize, ions: &[IonType]) -> Vec<bool> {
    let electron = ACTIVATIONS
        .get(activation)
        .is_some_and(|a| ELECTRON_ACTIVATIONS.contains(a));
    let seq = p.sequence.as_bytes();
    let n = seq.len();
    let (mut ph, mut ox) = (vec![0u32; n + 1], vec![0u32; n + 1]);
    for &(pos, id) in &p.mods {
        if (1..=n).contains(&pos) {
            let aa = seq[pos - 1];
            ph[pos] += u32::from(id == PHOSPHO && (aa == b'S' || aa == b'T'));
            ox[pos] += u32::from(id == OXIDATION && aa == b'M');
        }
    }
    let prefix = |v: &mut [u32]| {
        let mut acc = 0;
        for x in v.iter_mut() {
            acc += *x;
            *x = acc;
        }
    };
    prefix(&mut ph);
    prefix(&mut ox);
    let n_ion = ions.len();
    let mut v = vec![false; n.saturating_sub(1) * n_ion];
    for i in 0..n.saturating_sub(1) {
        let before_proline = seq[i + 1] == b'P';
        for (k, ion) in ions.iter().enumerate() {
            if ion.charge > charge {
                continue;
            }
            if ion.series.is_electron() && (!electron || before_proline) {
                continue;
            }
            let prefix = ion.is_prefix();
            let count = |c: &[u32]| if prefix { c[i + 1] } else { c[n] - c[i + 1] };
            v[i * n_ion + k] = match ion.loss {
                NeutralLoss::None => true,
                NeutralLoss::H3po4 => count(&ph) > 0,
                NeutralLoss::Ch4so => count(&ox) > 0,
            };
        }
    }
    v
}

struct Linear {
    din: usize,
    dout: usize,
    w: Vec<f32>,
    b: Vec<f32>,
}

impl Linear {
    fn from_torch(w: &Tensor, b: &Tensor) -> Result<Self, RtError> {
        let [dout, din] = w.dims[..] else {
            return Err(RtError::Format("linear weight must be 2-D".into()));
        };
        if b.dims != [dout] {
            return Err(RtError::Format("linear bias shape".into()));
        }
        Ok(Self {
            din,
            dout,
            w: w.data.clone(),
            b: b.data.clone(),
        })
    }

    fn apply(&self, x: &[f32], y: &mut Vec<f32>) {
        debug_assert_eq!(x.len(), self.din);
        y.clear();
        for o in 0..self.dout {
            let row = &self.w[o * self.din..(o + 1) * self.din];
            y.push(self.b[o] + row.iter().zip(x).map(|(a, b)| a * b).sum::<f32>());
        }
    }

    fn n_params(&self) -> usize {
        self.w.len() + self.b.len()
    }
}

struct Attention {
    heads: usize,
    qkv: Conv,
    proj: Conv,
}

pub struct Ms2Model {
    version: u8,
    emb: Vec<f32>,
    emb_dim: usize,
    ch: usize,
    act_emb: Vec<f32>,
    instr_emb: Vec<f32>,
    label_emb: Option<Vec<f32>>,
    n_act: usize,
    n_instr: usize,
    cond_emb: usize,
    film1: Linear,
    film2: Linear,
    inp: Conv,
    blocks: Vec<(Conv, Conv)>,
    attn: Option<Attention>,
    head1: Conv,
    head_cond: Linear,
    head2: Conv,
    head_ext: Option<Conv>,
}

#[derive(Default)]
pub struct Ms2Scratch {
    input: Vec<f32>,
    x: Vec<f32>,
    h: Vec<f32>,
    z: Vec<f32>,
    qkv: Vec<f32>,
    att: Vec<f32>,
    cond: Vec<f32>,
    film_h: Vec<f32>,
    film: Vec<f32>,
    hc: Vec<f32>,
    site: Vec<f32>,
    ext: Vec<f32>,
    cond_cache: Vec<CondEntry>,
    cond_model: usize,
}

struct CondEntry {
    key: Vec<u32>,
    film: Vec<f32>,
    hc: Vec<f32>,
}

const COND_CACHE: usize = 64;

impl Ms2Model {
    pub fn load(path: impl AsRef<Path>) -> Result<Self, RtError> {
        Self::from_bytes(&std::fs::read(path)?)
    }

    /// Parses an `AMS20001` (v1) or `AMS20002` (v2) weight file.
    pub fn from_bytes(buf: &[u8]) -> Result<Self, RtError> {
        let bad = |m: &str| RtError::Format(m.to_string());
        let (version, magic) = if buf.starts_with(MAGIC_V1) {
            (1u8, MAGIC_V1)
        } else if buf.starts_with(MAGIC_V2) {
            (2u8, MAGIC_V2)
        } else {
            return Err(bad("missing AMS20001 or AMS20002 magic"));
        };
        let t = read_tensors(buf, magic)?;
        let get = |k: &str| {
            t.get(k)
                .ok_or_else(|| RtError::Format(format!("missing tensor {k}")))
        };
        let conv = |name: &str, dil: usize| {
            Conv::from_torch(
                get(&format!("{name}.weight"))?,
                get(&format!("{name}.bias"))?,
                dil,
            )
        };
        let lin = |name: &str| {
            Linear::from_torch(
                get(&format!("{name}.weight"))?,
                get(&format!("{name}.bias"))?,
            )
        };
        let dil: Vec<usize> = get("dilations")?.data.iter().map(|&d| d as usize).collect();
        let heads = get("n_heads")?.data.first().copied().unwrap_or(0.0) as usize;
        let emb = get("emb.weight")?;
        let [ntok, emb_dim] = emb.dims[..] else {
            return Err(bad("emb shape"));
        };
        if ntok != N_TOK {
            return Err(bad("vocabulary size"));
        }
        let act = get("act_emb.weight")?;
        let ins = get("instr_emb.weight")?;
        let cond_emb = act.dims.get(1).copied().unwrap_or(0);
        let (n_act, n_instr) = if version >= 2 {
            (ACTIVATIONS.len(), INSTRUMENTS.len())
        } else {
            (V1_N_ACTIVATIONS, V1_N_INSTRUMENTS)
        };
        if act.dims != [n_act, cond_emb] || ins.dims != [n_instr, cond_emb] {
            return Err(bad("activation / instrument embedding shapes"));
        }
        let label_emb = if version >= 2 {
            let lab = get("label_emb.weight")?;
            if lab.dims != [LABELS.len(), cond_emb] {
                return Err(bad("label embedding shape"));
            }
            Some(lab.data.clone())
        } else {
            None
        };
        let n_cond = if label_emb.is_some() { 3 } else { 2 };
        let inp = conv("inp", 1)?;
        if inp.cin != emb_dim + N_COMP || inp.k != 1 {
            return Err(bad("inp shape"));
        }
        let ch = inp.cout;
        let mut blocks = Vec::new();
        for (i, &d) in dil.iter().enumerate() {
            let c = conv(&format!("convs.{i}"), d)?;
            let m = conv(&format!("mix.{i}"), 1)?;
            if c.cin != ch
                || c.cout != ch
                || m.cin != ch
                || m.cout != ch
                || c.k % 2 != 1
                || m.k != 1
            {
                return Err(bad("block shape"));
            }
            blocks.push((c, m));
        }
        let cond_dim = MAX_CHARGE + 2 + n_cond * cond_emb;
        let film1 = lin("film1")?;
        let film2 = lin("film2")?;
        if film1.din != cond_dim || film2.din != film1.dout || film2.dout != 2 * ch * dil.len() {
            return Err(bad("film shapes"));
        }
        let attn = if heads > 0 {
            let qkv = conv("qkv", 1)?;
            let proj = conv("proj", 1)?;
            if ch % heads != 0
                || qkv.cin != ch
                || qkv.cout != 3 * ch
                || proj.cin != ch
                || proj.cout != ch
            {
                return Err(bad("attention shapes"));
            }
            Some(Attention { heads, qkv, proj })
        } else {
            None
        };
        let head1 = conv("head1", 1)?;
        let head_cond = lin("head_cond")?;
        let head2 = conv("head2", 1)?;
        if head1.k != 2 || head1.cin != ch || head1.cout != ch {
            return Err(bad("head1 shape"));
        }
        if head_cond.din != cond_dim || head_cond.dout != ch {
            return Err(bad("head_cond shape"));
        }
        if head2.k != 1 || head2.cin != ch || head2.cout != N_ION {
            return Err(bad("head2 shape"));
        }
        let head_ext = if version >= 2 {
            let h = conv("head_ext", 1)?;
            if h.k != 1 || h.cin != ch || h.cout != N_ION_V2 - N_ION {
                return Err(bad("head_ext shape"));
            }
            Some(h)
        } else {
            None
        };
        Ok(Self {
            version,
            emb: emb.data.clone(),
            emb_dim,
            ch,
            act_emb: act.data.clone(),
            instr_emb: ins.data.clone(),
            label_emb,
            n_act,
            n_instr,
            cond_emb,
            film1,
            film2,
            inp,
            blocks,
            attn,
            head1,
            head_cond,
            head2,
            head_ext,
        })
    }

    /// Format version of the loaded file: 1 for `AMS20001`, 2 for `AMS20002`.
    pub fn version(&self) -> u8 {
        self.version
    }

    /// Output channels per site, in order: [`ION_TYPES`] for v1, [`ION_TYPES_V2`] for v2.
    pub fn ion_types(&self) -> &'static [IonType] {
        if self.head_ext.is_some() {
            &ION_TYPES_V2
        } else {
            &ION_TYPES
        }
    }

    /// Output channels per site: 12 for v1, `ION_TYPES_V2.len()` for v2.
    pub fn n_ion(&self) -> usize {
        self.ion_types().len()
    }

    pub fn n_params(&self) -> usize {
        let conv = |c: &Conv| c.w.len() + c.b.len();
        self.emb.len()
            + self.act_emb.len()
            + self.instr_emb.len()
            + self.label_emb.as_ref().map_or(0, Vec::len)
            + self.film1.n_params()
            + self.film2.n_params()
            + conv(&self.inp)
            + self
                .blocks
                .iter()
                .map(|(a, b)| conv(a) + conv(b))
                .sum::<usize>()
            + self
                .attn
                .as_ref()
                .map_or(0, |a| conv(&a.qkv) + conv(&a.proj))
            + conv(&self.head1)
            + self.head_cond.n_params()
            + conv(&self.head2)
            + self.head_ext.as_ref().map_or(0, conv)
    }

    pub fn predict(&self, p: &RtPeptide, c: &Ms2Conditions) -> Option<Vec<f32>> {
        self.predict_with(p, c, &mut Ms2Scratch::default())
    }

    pub fn predict_with(
        &self,
        p: &RtPeptide,
        c: &Ms2Conditions,
        s: &mut Ms2Scratch,
    ) -> Option<Vec<f32>> {
        let mut out = self.predict_raw_with(p, c, s)?;
        let mask = channel_mask(p, c.charge, c.activation, self.ion_types());
        for (v, ok) in out.iter_mut().zip(mask) {
            *v = if ok { v.max(0.0) } else { 0.0 };
        }
        Some(out)
    }

    pub fn predict_raw_with(
        &self,
        p: &RtPeptide,
        c: &Ms2Conditions,
        s: &mut Ms2Scratch,
    ) -> Option<Vec<f32>> {
        if c.charge == 0 || p.sequence.len() < 2 {
            return None;
        }
        let (tok, comp) = p.encode()?;
        let np = tok.len();
        let n = np - 2;
        let (ed, ch) = (self.emb_dim, self.ch);
        let known = |i: usize, rows: usize| if i < rows { i } else { 0 };
        let act = known(c.activation, self.n_act);
        let ins = known(c.instrument, self.n_instr);

        let ce = self.cond_emb;
        s.cond.clear();
        s.cond.extend_from_slice(&c.numeric());
        s.cond.extend_from_slice(&self.act_emb[act * ce..(act + 1) * ce]);
        s.cond.extend_from_slice(&self.instr_emb[ins * ce..(ins + 1) * ce]);
        if let Some(lab_emb) = &self.label_emb {
            let lab = known(c.label, LABELS.len());
            s.cond.extend_from_slice(&lab_emb[lab * ce..(lab + 1) * ce]);
        }
        let model = self as *const Self as usize;
        if s.cond_model != model {
            s.cond_cache.clear();
            s.cond_model = model;
        }
        let mut cached = s.cond_cache.iter().position(|e| {
            e.key.len() == s.cond.len() && e.key.iter().zip(&s.cond).all(|(k, c)| *k == c.to_bits())
        });
        if cached.is_none() {
            self.film1.apply(&s.cond, &mut s.film_h);
            for v in s.film_h.iter_mut() {
                *v = silu(*v);
            }
            self.film2.apply(&s.film_h, &mut s.film);
            self.head_cond.apply(&s.cond, &mut s.hc);
            if s.cond_cache.len() < COND_CACHE {
                s.cond_cache.push(CondEntry {
                    key: s.cond.iter().map(|c| c.to_bits()).collect(),
                    film: s.film.clone(),
                    hc: s.hc.clone(),
                });
                cached = Some(s.cond_cache.len() - 1);
            }
        }
        let (film, hc): (&[f32], &[f32]) = match cached {
            Some(i) => (&s.cond_cache[i].film, &s.cond_cache[i].hc),
            None => (&s.film, &s.hc),
        };

        s.input.clear();
        for (t, cv) in tok.iter().zip(&comp) {
            s.input.extend_from_slice(&self.emb[t * ed..(t + 1) * ed]);
            s.input.extend_from_slice(cv);
        }
        s.x.resize(np * ch, 0.0);
        s.h.resize(np * ch, 0.0);
        s.z.resize(np * ch, 0.0);
        self.inp.apply(&s.input, np, &mut s.x);

        for (b, (conv, mix)) in self.blocks.iter().enumerate() {
            conv.apply(&s.x, np, &mut s.h);
            let gamma = &film[(2 * b) * ch..(2 * b + 1) * ch];
            let beta = &film[(2 * b + 1) * ch..(2 * b + 2) * ch];
            for row in s.h.chunks_exact_mut(ch) {
                for ((v, g), be) in row.iter_mut().zip(gamma).zip(beta) {
                    *v = silu(*v * (1.0 + g) + be);
                }
            }
            mix.apply(&s.h, np, &mut s.z);
            for (x, z) in s.x.iter_mut().zip(&s.z) {
                *x += z;
            }
        }

        if let Some(a) = &self.attn {
            let dh = ch / a.heads;
            let scale = 1.0 / (dh as f32).sqrt();
            s.qkv.resize(np * 3 * ch, 0.0);
            a.qkv.apply(&s.x, np, &mut s.qkv);
            s.att.resize(np, 0.0);
            for v in s.h.iter_mut() {
                *v = 0.0;
            }
            for hd in 0..a.heads {
                let off = hd * dh;
                for i in 0..np {
                    let qi = &s.qkv[i * 3 * ch + off..i * 3 * ch + off + dh];
                    let mut mx = f32::NEG_INFINITY;
                    for (j, a) in s.att.iter_mut().enumerate() {
                        let kj = &s.qkv[j * 3 * ch + ch + off..j * 3 * ch + ch + off + dh];
                        *a = qi.iter().zip(kj).map(|(x, y)| x * y).sum::<f32>() * scale;
                        mx = mx.max(*a);
                    }
                    let mut sum = 0.0f32;
                    for w in s.att.iter_mut() {
                        *w = (*w - mx).exp();
                        sum += *w;
                    }
                    let oi = &mut s.h[i * ch + off..i * ch + off + dh];
                    for (j, &a) in s.att.iter().enumerate() {
                        let w = a / sum;
                        let vj = &s.qkv[j * 3 * ch + 2 * ch + off..j * 3 * ch + 2 * ch + off + dh];
                        for (o, v) in oi.iter_mut().zip(vj) {
                            *o += w * v;
                        }
                    }
                }
            }
            a.proj.apply(&s.h, np, &mut s.z);
            for (x, z) in s.x.iter_mut().zip(&s.z) {
                *x += z;
            }
        }

        self.head1.apply(&s.x, np, &mut s.h);
        s.site.clear();
        for j in 0..n - 1 {
            let row = &s.h[(j + 2) * ch..(j + 3) * ch];
            s.site.extend(row.iter().zip(hc).map(|(a, b)| silu(a + b)));
        }
        let mut out = vec![0.0f32; (n - 1) * N_ION];
        self.head2.apply(&s.site, n - 1, &mut out);
        let Some(ext) = &self.head_ext else {
            return Some(out);
        };
        let n_ext = ext.cout;
        s.ext.resize((n - 1) * n_ext, 0.0);
        ext.apply(&s.site, n - 1, &mut s.ext);
        let mut full = Vec::with_capacity((n - 1) * (N_ION + n_ext));
        for (base, more) in out.chunks_exact(N_ION).zip(s.ext.chunks_exact(n_ext)) {
            full.extend_from_slice(base);
            full.extend_from_slice(more);
        }
        Some(full)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> std::path::PathBuf {
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("fixtures/ms2_native")
            .join(name)
    }

    /// FNV-1a 64 of a fixture file, carriage returns dropped so a CRLF checkout hashes the same.
    fn fixture_hash(name: &str) -> u64 {
        std::fs::read(fixture(name))
            .unwrap()
            .into_iter()
            .filter(|&b| b != b'\r')
            .fold(0xcbf2_9ce4_8422_2325u64, |h, b| {
                (h ^ u64::from(b)).wrapping_mul(0x0000_0100_0000_01b3)
            })
    }

    /// Checks every row of a torch reference TSV against `model`; returns (encodable rows,
    /// rows neither side encodes, largest relative difference). A missing `label` column
    /// means label "none".
    fn check_reference(model: &Ms2Model, reference: &str) -> (usize, usize, f32) {
        let text = std::fs::read_to_string(fixture(reference)).unwrap();
        let mut lines = text.lines();
        let h: Vec<&str> = lines.next().unwrap().split('\t').collect();
        let col = |n: &str| h.iter().position(|x| *x == n);
        let need = |n: &str| col(n).unwrap();
        let (cf, cp, cz, cn, ca, ci, cv, cr) = (
            need("format"),
            need("peptide"),
            need("charge"),
            need("nce"),
            need("activation"),
            need("instrument"),
            need("values"),
            need("raw"),
        );
        let cl = col("label");
        let parse = |s: &str| {
            s.split(',')
                .map(|v| v.parse::<f32>().unwrap())
                .collect::<Vec<_>>()
        };
        let (mut n_ok, mut n_none, mut max_err) = (0, 0, 0f32);
        let mut s = Ms2Scratch::default();
        for line in lines {
            let x: Vec<&str> = line.split('\t').collect();
            let p = match x[cf] {
                "proforma" => RtPeptide::from_proforma(x[cp]),
                _ => RtPeptide::from_andes(x[cp]),
            };
            let label = cl.map_or("none", |i| x[i]);
            let c = Ms2Conditions::new(
                x[cz].parse().unwrap(),
                x[cn].parse().ok(),
                x[ca],
                x[ci],
                label,
            );
            let got = p.as_ref().and_then(|p| model.predict_with(p, &c, &mut s));
            let raw = p
                .as_ref()
                .and_then(|p| model.predict_raw_with(p, &c, &mut s));
            if x[cv] == "NA" {
                assert!(got.is_none(), "{} should not encode", x[cp]);
                n_none += 1;
                continue;
            }
            let (got, raw) = (got.expect(x[cp]), raw.unwrap());
            let (want, want_raw) = (parse(x[cv]), parse(x[cr]));
            let sites = p.as_ref().unwrap().sequence.len() - 1;
            assert_eq!(got.len(), want.len(), "{}", x[cp]);
            assert_eq!(got.len(), sites * model.n_ion(), "{}", x[cp]);
            assert_eq!(raw.len(), want_raw.len(), "{}", x[cp]);
            for (g, w) in got.iter().zip(&want).chain(raw.iter().zip(&want_raw)) {
                max_err = max_err.max((g - w).abs() / (1.0 + w.abs()));
            }
            n_ok += 1;
        }
        (n_ok, n_none, max_err)
    }

    fn bits(v: &[f32]) -> Vec<u32> {
        v.iter().map(|x| x.to_bits()).collect()
    }

    fn index_zero(charge: u8) -> Ms2Conditions {
        Ms2Conditions {
            charge,
            nce: None,
            activation: 0,
            instrument: 0,
            label: 0,
        }
    }

    #[test]
    fn instrument_and_activation_vocab() {
        assert_eq!(instrument_class("Q Exactive HF-X"), 1);
        assert_eq!(instrument_class("Orbitrap Exploris 480"), 2);
        assert_eq!(instrument_class("Orbitrap Fusion Lumos"), 3);
        assert_eq!(instrument_class("LTQ Orbitrap Velos"), 4);
        assert_eq!(instrument_class("Orbitrap Astral"), 5);
        assert_eq!(instrument_class("LTQ"), 6);
        assert_eq!(instrument_class("timsTOF Pro 2"), 7);
        assert_eq!(instrument_class("TripleTOF 6600"), 8);
        assert_eq!(instrument_class(""), 0);
        assert_eq!(activation_index("CID-ITMS"), 3);
        assert_eq!(activation_index("ETD"), 6);
        assert_eq!(activation_index("EThcD"), 7);
        assert_eq!(activation_index("HCD-TOF"), 8);
        assert_eq!(activation_index("XYZ"), 0);
        assert_eq!(label_index("none"), 0);
        assert_eq!(label_index("TMTpro"), 1);
        assert_eq!(label_index("iTRAQ8plex"), 2);
        assert_eq!(
            ACTIVATIONS[..V1_N_ACTIVATIONS],
            ["unknown", "HCD-FTMS", "HCD-ITMS", "CID-ITMS", "CID-FTMS", "HCD-ASTMS"]
        );
        assert_eq!(
            INSTRUMENTS[..V1_N_INSTRUMENTS],
            ["unknown", "QExactive", "Exploris", "Tribrid", "VelosElite", "Astral", "LTQ"]
        );
    }

    #[test]
    fn loss_channels_follow_the_modified_residues() {
        let p = RtPeptide::from_proforma("M[UNIMOD:35]PS[UNIMOD:21]K").unwrap();
        let v = valid_mask(&p, 1);
        let col = |k: usize| (0..3).map(|i| v[i * N_ION + k]).collect::<Vec<_>>();
        assert_eq!(col(1), vec![false; 3]);
        assert_eq!(col(4), vec![false, false, true]);
        assert_eq!(col(6), vec![true, true, false]);
        assert_eq!(col(8), vec![true, true, true]);
        assert_eq!(col(10), vec![false, false, false]);
        assert!(valid_mask(&p, 2)[1]);
    }

    #[test]
    fn v2_channels_keep_the_v1_order_and_mask_impossible_ions() {
        assert_eq!(ION_TYPES_V2[..N_ION], ION_TYPES);
        let p = RtPeptide::from_proforma("PEPTIDE").unwrap();
        let etd = channel_mask(&p, 3, activation_index("ETD"), &ION_TYPES_V2);
        let col = |v: &[bool], k: usize| (0..6).map(|i| v[i * N_ION_V2 + k]).collect::<Vec<_>>();
        assert_eq!(col(&etd, 12), vec![true; 6]);
        assert_eq!(col(&etd, 14), vec![true; 6]);
        assert_eq!(col(&etd, 15), vec![true, false, true, true, true, true]);
        assert_eq!(col(&etd, 17), vec![true, false, true, true, true, true]);
        let hcd = channel_mask(&p, 3, activation_index("HCD-FTMS"), &ION_TYPES_V2);
        assert!((15..N_ION_V2).all(|k| col(&hcd, k) == vec![false; 6]));
        let z2 = channel_mask(&p, 2, activation_index("ETD"), &ION_TYPES_V2);
        assert_eq!(col(&z2, 12), vec![false; 6]);
        assert_eq!(col(&z2, 13), vec![false; 6]);
        assert_eq!(channel_mask(&p, 2, 0, &ION_TYPES), valid_mask(&p, 2));
    }

    #[test]
    fn v1_fixture_bit_identical() {
        assert_eq!(fixture_hash("tiny.ams2"), 0xfc2c_f35f_cbc2_0fae, "tiny.ams2 changed");
        assert_eq!(
            fixture_hash("tiny_reference.tsv"),
            0x96f2_a0f2_83ad_2435,
            "tiny_reference.tsv changed"
        );
        let model = Ms2Model::load(fixture("tiny.ams2")).unwrap();
        assert_eq!(model.version(), 1);
        assert_eq!(model.n_params(), 1614);
        let (n_ok, n_none, max_err) = check_reference(&model, "tiny_reference.tsv");
        assert_eq!((n_ok, n_none), (21, 5));
        assert!(max_err < 1e-4, "max relative error {max_err}");
        let tmt = "[UNIMOD:737]-M[UNIMOD:35]PEPS[UNIMOD:21]K[UNIMOD:737]";
        let p = RtPeptide::from_proforma(tmt).unwrap();
        let mut s = Ms2Scratch::default();
        for z in 1..=4u8 {
            let plain = Ms2Conditions::new(z, Some(28.0), "HCD-FTMS", "Q Exactive", "none");
            let want = model.predict_with(&p, &plain, &mut s).unwrap();
            for label in ["TMT", "iTRAQ", "bogus"] {
                let c = Ms2Conditions::new(z, Some(28.0), "HCD-FTMS", "Q Exactive", label);
                let got = model.predict_with(&p, &c, &mut s).unwrap();
                assert_eq!(bits(&got), bits(&want), "v1 must ignore label {label}");
            }
        }
    }

    #[test]
    fn v2_fixture_parity() {
        let model = Ms2Model::load(fixture("tiny_v2.ams2")).unwrap();
        assert_eq!(model.version(), 2);
        assert_eq!(model.n_params(), 1717);
        let (n_ok, n_none, max_err) = check_reference(&model, "tiny_v2_reference.tsv");
        assert_eq!((n_ok, n_none), (23, 5));
        assert!(max_err < 1e-4, "max relative error {max_err}");
    }

    #[test]
    fn unknown_condition_maps_to_zero() {
        let p = RtPeptide::from_proforma("PEPTIDEK").unwrap();
        for name in ["tiny.ams2", "tiny_v2.ams2"] {
            let model = Ms2Model::load(fixture(name)).unwrap();
            let mut s = Ms2Scratch::default();
            for z in 1..=4u8 {
                let zero = index_zero(z);
                let want_raw = model.predict_raw_with(&p, &zero, &mut s).unwrap();
                let want = model.predict_with(&p, &zero, &mut s).unwrap();
                let unknown = Ms2Conditions::new(z, None, "XYZ", "", "none");
                assert_eq!(unknown, zero);
                let out_of_range = Ms2Conditions {
                    activation: usize::MAX,
                    instrument: INSTRUMENTS.len(),
                    label: usize::MAX,
                    ..zero
                };
                for c in [unknown, out_of_range] {
                    let raw = model.predict_raw_with(&p, &c, &mut s).unwrap();
                    let got = model.predict_with(&p, &c, &mut s).unwrap();
                    assert_eq!(bits(&raw), bits(&want_raw), "{name} z{z} {c:?}");
                    assert_eq!(bits(&got), bits(&want), "{name} z{z} {c:?}");
                }
            }
        }
        let v1 = Ms2Model::load(fixture("tiny.ams2")).unwrap();
        let mut s = Ms2Scratch::default();
        let want = v1.predict_raw_with(&p, &index_zero(2), &mut s).unwrap();
        let unseen = Ms2Conditions::new(2, None, "EThcD", "timsTOF Pro 2", "TMT");
        assert_eq!((unseen.activation, unseen.instrument), (7, 7));
        let got = v1.predict_raw_with(&p, &unseen, &mut s).unwrap();
        assert_eq!(bits(&got), bits(&want), "v1 maps v2-only indices to 0");
    }

    #[test]
    fn v1_model_reports_12_channels() {
        let v1 = Ms2Model::load(fixture("tiny.ams2")).unwrap();
        assert_eq!(v1.n_ion(), 12);
        assert_eq!(v1.ion_types(), &ION_TYPES[..]);
        let v2 = Ms2Model::load(fixture("tiny_v2.ams2")).unwrap();
        assert_eq!(v2.n_ion(), ION_TYPES_V2.len());
        assert_eq!(v2.ion_types(), &ION_TYPES_V2[..]);
        let p = RtPeptide::from_proforma("PEPTIDEK").unwrap();
        let c = Ms2Conditions::new(3, Some(30.0), "HCD-FTMS", "Q Exactive", "none");
        assert_eq!(v1.predict(&p, &c).unwrap().len(), 7 * 12);
        assert_eq!(v2.predict(&p, &c).unwrap().len(), 7 * ION_TYPES_V2.len());
    }

    #[test]
    fn reused_scratch_matches_fresh_scratch() {
        for name in ["tiny.ams2", "tiny_v2.ams2"] {
            let model = Ms2Model::load(fixture(name)).unwrap();
            let peptides = ["PEPTIDEK", "AC[UNIMOD:4]DEFMK", "LLGGSR"];
            let labels = ["none", "TMT", "iTRAQ"];
            let mut reused = Ms2Scratch::default();
            for round in 0..2 {
                for (i, z) in [2u8, 3, 1, 2, 4, 3].into_iter().enumerate() {
                    let nce = if i % 2 == 0 { Some(27.0) } else { None };
                    let c = Ms2Conditions::new(z, nce, "HCD", "Q Exactive", labels[i % 3]);
                    let p =
                        RtPeptide::from_proforma(peptides[(i + round) % peptides.len()]).unwrap();
                    let fresh = model
                        .predict_raw_with(&p, &c, &mut Ms2Scratch::default())
                        .unwrap();
                    let again = model.predict_raw_with(&p, &c, &mut reused).unwrap();
                    assert_eq!(bits(&fresh), bits(&again), "{name} charge {z} nce {nce:?}");
                }
            }
        }
    }

    #[test]
    fn rejects_other_formats() {
        assert!(Ms2Model::from_bytes(b"ARTN0001\0\0\0\0").is_err());
        assert!(Ms2Model::from_bytes(b"AMS20003\0\0\0\0").is_err());
        assert!(Ms2Model::from_bytes(b"AMS2").is_err());
    }
}
