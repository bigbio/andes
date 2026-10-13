//! andes-ms2: a small fragment-intensity CNN evaluated in plain Rust (no ONNX runtime).
//!
//! The network is defined in `scripts/student/ms2/ms2_model.py`; this module
//! reproduces its forward pass from an `AMS20001` weight file. Peptides are
//! encoded as in [`rt_native`](crate::rt_native) (tokens with I as L, plus the
//! elemental composition of each position's modifications), so any
//! modification of the built-in Unimod table is accepted.
//!
//! The model is conditioned on precursor charge, NCE (or "unknown"), activation
//! ([`ACTIVATIONS`]) and instrument family ([`INSTRUMENTS`]). It predicts
//! [`N_ION`] channels per cleavage site ([`ION_TYPES`]); site `i` is the cleavage
//! after residue `i + 1`, its b ion the N-terminal fragment and its y ion the
//! complementary C-terminal fragment. The first four channels are b 1+, b 2+,
//! y 1+, y 2+ (the PeptDeep `b_z1, b_z2, y_z1, y_z2` order).

use std::path::Path;

use crate::rt_native::mods::N_COMP;
use crate::rt_native::{read_tensors, silu, Conv, RtError, RtPeptide, Tensor};

/// Unimod accession of phosphorylation.
const PHOSPHO: u32 = 21;
/// Unimod accession of oxidation.
const OXIDATION: u32 = 35;
const N_TOK: usize = 23;

/// Neutral loss of an ion channel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NeutralLoss {
    /// Intact fragment.
    None,
    /// −H3PO4 (97.9769 Da), from phosphorylated S/T.
    H3po4,
    /// −CH4SO (63.9983 Da), from oxidised M.
    Ch4so,
}

/// One output channel: (b ion?, fragment charge, neutral loss).
pub type IonType = (bool, u8, NeutralLoss);

/// Output channels per cleavage site, in model order.
pub const ION_TYPES: [IonType; 12] = [
    (true, 1, NeutralLoss::None),
    (true, 2, NeutralLoss::None),
    (false, 1, NeutralLoss::None),
    (false, 2, NeutralLoss::None),
    (true, 1, NeutralLoss::H3po4),
    (true, 2, NeutralLoss::H3po4),
    (false, 1, NeutralLoss::H3po4),
    (false, 2, NeutralLoss::H3po4),
    (true, 1, NeutralLoss::Ch4so),
    (true, 2, NeutralLoss::Ch4so),
    (false, 1, NeutralLoss::Ch4so),
    (false, 2, NeutralLoss::Ch4so),
];
/// Channels per cleavage site.
pub const N_ION: usize = ION_TYPES.len();
/// Channels without a neutral loss (b 1+, b 2+, y 1+, y 2+).
pub const N_BASE_ION: usize = 4;
/// Activation vocabulary (dissociation + MS2 analyser); index 0 = unknown.
pub const ACTIVATIONS: [&str; 6] = [
    "unknown",
    "HCD-FTMS",
    "HCD-ITMS",
    "CID-ITMS",
    "CID-FTMS",
    "HCD-ASTMS",
];
/// Instrument families; index 0 = unknown.
pub const INSTRUMENTS: [&str; 7] = [
    "unknown",
    "QExactive",
    "Exploris",
    "Tribrid",
    "VelosElite",
    "Astral",
    "LTQ",
];
/// Precursor charges are one-hot encoded up to this value; higher charges share the last slot.
pub const MAX_CHARGE: usize = 6;

/// Index in [`ACTIVATIONS`] (0 when not listed).
pub fn activation_index(name: &str) -> usize {
    ACTIVATIONS.iter().position(|a| *a == name).unwrap_or(0)
}

/// Index in [`INSTRUMENTS`] for a free-text instrument model name (0 = unknown).
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
    } else {
        0
    }
}

/// Acquisition conditions the prediction is conditioned on.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ms2Conditions {
    /// Precursor charge (≥ 1).
    pub charge: u8,
    /// Normalised collision energy, `None` when unknown.
    pub nce: Option<f32>,
    /// Index in [`ACTIVATIONS`].
    pub activation: usize,
    /// Index in [`INSTRUMENTS`].
    pub instrument: usize,
}

impl Ms2Conditions {
    /// Conditions from an activation name ([`ACTIVATIONS`]) and an instrument model name.
    pub fn new(charge: u8, nce: Option<f32>, activation: &str, instrument: &str) -> Self {
        Self {
            charge,
            nce: nce.filter(|v| v.is_finite()),
            activation: activation_index(activation),
            instrument: instrument_class(instrument),
        }
    }

    /// Charge one-hot followed by `[NCE / 100, NCE known]`.
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

/// Which channels exist for `p` at precursor `charge`: `(n − 1) × N_ION`,
/// site-major. A channel is absent when its fragment charge exceeds the
/// precursor charge or, for a loss channel, when the fragment holds no pS/pT
/// (H3PO4) or oxidised M (CH4SO).
pub fn valid_mask(p: &RtPeptide, charge: u8) -> Vec<bool> {
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
    for i in 1..=n {
        ph[i] += ph[i - 1];
        ox[i] += ox[i - 1];
    }
    let mut v = vec![false; n.saturating_sub(1) * N_ION];
    for i in 0..n.saturating_sub(1) {
        for (k, &(b, z, loss)) in ION_TYPES.iter().enumerate() {
            if z > charge {
                continue;
            }
            let count = |c: &[u32]| if b { c[i + 1] } else { c[n] - c[i + 1] };
            v[i * N_ION + k] = match loss {
                NeutralLoss::None => true,
                NeutralLoss::H3po4 => count(&ph) > 0,
                NeutralLoss::Ch4so => count(&ox) > 0,
            };
        }
    }
    v
}

/// Dense layer, torch layout `w[out][in]`.
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

/// The andes-ms2 network.
pub struct Ms2Model {
    emb: Vec<f32>,
    emb_dim: usize,
    ch: usize,
    act_emb: Vec<f32>,
    instr_emb: Vec<f32>,
    cond_emb: usize,
    film1: Linear,
    film2: Linear,
    inp: Conv,
    blocks: Vec<(Conv, Conv)>,
    attn: Option<Attention>,
    head1: Conv,
    head_cond: Linear,
    head2: Conv,
}

/// Reusable per-thread buffers for [`Ms2Model::predict_with`].
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
}

impl Ms2Model {
    /// Load an `AMS20001` weight file written by `ms2_model.export_weights`.
    pub fn load(path: impl AsRef<Path>) -> Result<Self, RtError> {
        Self::from_bytes(&std::fs::read(path)?)
    }

    /// Parse an `AMS20001` blob.
    pub fn from_bytes(buf: &[u8]) -> Result<Self, RtError> {
        let bad = |m: &str| RtError::Format(m.to_string());
        let t = read_tensors(buf, b"AMS20001")?;
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
            Linear::from_torch(get(&format!("{name}.weight"))?, get(&format!("{name}.bias"))?)
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
        if act.dims != [ACTIVATIONS.len(), cond_emb] || ins.dims != [INSTRUMENTS.len(), cond_emb] {
            return Err(bad("activation / instrument embedding shapes"));
        }
        let inp = conv("inp", 1)?;
        if inp.cin != emb_dim + N_COMP || inp.k != 1 {
            return Err(bad("inp shape"));
        }
        let ch = inp.cout;
        let mut blocks = Vec::new();
        for (i, &d) in dil.iter().enumerate() {
            let c = conv(&format!("convs.{i}"), d)?;
            let m = conv(&format!("mix.{i}"), 1)?;
            if c.cin != ch || c.cout != ch || m.cin != ch || m.cout != ch || c.k % 2 != 1 || m.k != 1 {
                return Err(bad("block shape"));
            }
            blocks.push((c, m));
        }
        let cond_dim = MAX_CHARGE + 2 + 2 * cond_emb;
        let film1 = lin("film1")?;
        let film2 = lin("film2")?;
        if film1.din != cond_dim || film2.din != film1.dout || film2.dout != 2 * ch * dil.len() {
            return Err(bad("film shapes"));
        }
        let attn = if heads > 0 {
            let qkv = conv("qkv", 1)?;
            let proj = conv("proj", 1)?;
            if ch % heads != 0 || qkv.cin != ch || qkv.cout != 3 * ch || proj.cin != ch || proj.cout != ch {
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
        Ok(Self {
            emb: emb.data.clone(),
            emb_dim,
            ch,
            act_emb: act.data.clone(),
            instr_emb: ins.data.clone(),
            cond_emb,
            film1,
            film2,
            inp,
            blocks,
            attn,
            head1,
            head_cond,
            head2,
        })
    }

    /// Number of parameters.
    pub fn n_params(&self) -> usize {
        let conv = |c: &Conv| c.w.len() + c.b.len();
        self.emb.len()
            + self.act_emb.len()
            + self.instr_emb.len()
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
    }

    /// Predicted intensities, `(n − 1) × N_ION` site-major, clipped at 0 and with
    /// absent channels ([`valid_mask`]) set to 0. `None` when the peptide cannot
    /// be encoded, has fewer than two residues, or `charge` is 0.
    pub fn predict(&self, p: &RtPeptide, c: &Ms2Conditions) -> Option<Vec<f32>> {
        self.predict_with(p, c, &mut Ms2Scratch::default())
    }

    /// As [`predict`](Self::predict), reusing `scratch` across calls.
    pub fn predict_with(&self, p: &RtPeptide, c: &Ms2Conditions, s: &mut Ms2Scratch) -> Option<Vec<f32>> {
        let mut out = self.predict_raw_with(p, c, s)?;
        for (v, ok) in out.iter_mut().zip(valid_mask(p, c.charge)) {
            *v = if ok { v.max(0.0) } else { 0.0 };
        }
        Some(out)
    }

    /// The unprocessed network output (same layout as [`predict`](Self::predict)).
    pub fn predict_raw_with(&self, p: &RtPeptide, c: &Ms2Conditions, s: &mut Ms2Scratch) -> Option<Vec<f32>> {
        if c.charge == 0 || p.sequence.len() < 2 {
            return None;
        }
        let (tok, comp) = p.encode()?;
        let np = tok.len();
        let n = np - 2;
        let (ed, ch) = (self.emb_dim, self.ch);
        if c.activation >= ACTIVATIONS.len() || c.instrument >= INSTRUMENTS.len() {
            return None;
        }

        // conditioning
        let ce = self.cond_emb;
        s.cond.clear();
        s.cond.extend_from_slice(&c.numeric());
        s.cond.extend_from_slice(&self.act_emb[c.activation * ce..(c.activation + 1) * ce]);
        s.cond.extend_from_slice(&self.instr_emb[c.instrument * ce..(c.instrument + 1) * ce]);
        self.film1.apply(&s.cond, &mut s.film_h);
        for v in s.film_h.iter_mut() {
            *v = silu(*v);
        }
        self.film2.apply(&s.film_h, &mut s.film);

        // input projection
        s.input.clear();
        for (t, cv) in tok.iter().zip(&comp) {
            s.input.extend_from_slice(&self.emb[t * ed..(t + 1) * ed]);
            s.input.extend_from_slice(cv);
        }
        s.x.resize(np * ch, 0.0);
        s.h.resize(np * ch, 0.0);
        s.z.resize(np * ch, 0.0);
        self.inp.apply(&s.input, np, &mut s.x);

        // FiLM-modulated residual dilated convolutions
        for (b, (conv, mix)) in self.blocks.iter().enumerate() {
            conv.apply(&s.x, np, &mut s.h);
            let gamma = &s.film[(2 * b) * ch..(2 * b + 1) * ch];
            let beta = &s.film[(2 * b + 1) * ch..(2 * b + 2) * ch];
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

        // multi-head self-attention
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
                    for j in 0..np {
                        let kj = &s.qkv[j * 3 * ch + ch + off..j * 3 * ch + ch + off + dh];
                        let sc = qi.iter().zip(kj).map(|(a, b)| a * b).sum::<f32>() * scale;
                        s.att[j] = sc;
                        mx = mx.max(sc);
                    }
                    let mut sum = 0.0f32;
                    for w in s.att.iter_mut() {
                        *w = (*w - mx).exp();
                        sum += *w;
                    }
                    let oi = &mut s.h[i * ch + off..i * ch + off + dh];
                    for j in 0..np {
                        let w = s.att[j] / sum;
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

        // site head: head1 is a k = 2 conv, so its output at token p reads tokens
        // p − 1 and p; site j (between residues j + 1 and j + 2) is token j + 2.
        self.head1.apply(&s.x, np, &mut s.h);
        self.head_cond.apply(&s.cond, &mut s.hc);
        s.site.clear();
        for j in 0..n - 1 {
            let row = &s.h[(j + 2) * ch..(j + 3) * ch];
            s.site.extend(row.iter().zip(&s.hc).map(|(a, b)| silu(a + b)));
        }
        let mut out = vec![0.0f32; (n - 1) * N_ION];
        self.head2.apply(&s.site, n - 1, &mut out);
        Some(out)
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

    #[test]
    fn instrument_and_activation_vocab() {
        assert_eq!(instrument_class("Q Exactive HF-X"), 1);
        assert_eq!(instrument_class("Orbitrap Exploris 480"), 2);
        assert_eq!(instrument_class("Orbitrap Fusion Lumos"), 3);
        assert_eq!(instrument_class("LTQ Orbitrap Velos"), 4);
        assert_eq!(instrument_class("Orbitrap Astral"), 5);
        assert_eq!(instrument_class("LTQ"), 6);
        assert_eq!(instrument_class(""), 0);
        assert_eq!(activation_index("CID-ITMS"), 3);
        assert_eq!(activation_index("ETD"), 0);
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

    /// Parity with the PyTorch model on the committed random-weight fixture
    /// (`ms2_export.py --random`), in ProForma and andes PIN notation.
    #[test]
    fn matches_pytorch_reference() {
        let model = Ms2Model::load(fixture("tiny.ams2")).unwrap();
        assert_eq!(model.n_params(), 1614);
        let text = std::fs::read_to_string(fixture("tiny_reference.tsv")).unwrap();
        let mut lines = text.lines();
        let h: Vec<&str> = lines.next().unwrap().split('\t').collect();
        let col = |n: &str| h.iter().position(|x| *x == n).unwrap();
        let (cf, cp, cz, cn, ca, ci, cv, cr) = (
            col("format"),
            col("peptide"),
            col("charge"),
            col("nce"),
            col("activation"),
            col("instrument"),
            col("values"),
            col("raw"),
        );
        let (mut n_ok, mut n_none, mut max_err) = (0, 0, 0f32);
        let mut s = Ms2Scratch::default();
        for line in lines {
            let x: Vec<&str> = line.split('\t').collect();
            let p = match x[cf] {
                "proforma" => RtPeptide::from_proforma(x[cp]),
                _ => RtPeptide::from_andes(x[cp]),
            };
            let c = Ms2Conditions::new(x[cz].parse().unwrap(), x[cn].parse().ok(), x[ca], x[ci]);
            let got = p.as_ref().and_then(|p| model.predict_with(p, &c, &mut s));
            let raw = p.as_ref().and_then(|p| model.predict_raw_with(p, &c, &mut s));
            if x[cv] == "NA" {
                assert!(got.is_none(), "{} should not encode", x[cp]);
                n_none += 1;
                continue;
            }
            let parse = |s: &str| s.split(',').map(|v| v.parse::<f32>().unwrap()).collect::<Vec<_>>();
            let (got, raw) = (got.expect(x[cp]), raw.unwrap());
            for (g, w) in got.iter().zip(parse(x[cv])).chain(raw.iter().zip(parse(x[cr]))) {
                max_err = max_err.max((g - w).abs() / (1.0 + w.abs()));
            }
            assert_eq!(got.len(), parse(x[cv]).len(), "{}", x[cp]);
            n_ok += 1;
        }
        assert_eq!((n_ok, n_none), (21, 5));
        assert!(max_err < 1e-4, "max relative error {max_err}");
    }

    #[test]
    fn rejects_other_formats() {
        assert!(Ms2Model::from_bytes(b"ARTN0001\0\0\0\0").is_err());
    }
}
