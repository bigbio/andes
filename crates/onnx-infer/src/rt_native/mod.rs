//! andes-rt: a small retention-time CNN evaluated in plain Rust (no ONNX runtime).
//!
//! The network and the encoding are defined in `scripts/student/rt_model.py`;
//! this module reproduces its forward pass from an `ARTN0001` weight file.
//!
//! Encoding: tokens `[N-term, residues…, C-term]` (N-term = 1, C-term = 2,
//! residue = 3 + index in [`AA`]; I is encoded as L) and, per position, the summed elemental
//! composition of its modifications (see [`mods`]) scaled by [`COMP_SCALE`].
//! N-terminal mods sit on the N-term token, C-terminal mods on the C-term token,
//! so an unseen PTM is described by its composition rather than by a learned id.
//!
//! The output is a dimensionless retention index; map it to minutes per run
//! with a monotone calibration on confident PSMs.

pub mod mods;

use std::collections::HashMap;
use std::io::Read;
use std::path::Path;

use mods::N_COMP;

/// Residue alphabet, in token order.
pub const AA: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
/// Longest supported peptide.
pub const MAX_LEN: usize = 50;
/// Scale applied to the composition features before the network.
pub const COMP_SCALE: [f32; N_COMP] = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 1.0];
const N_TOK: usize = 23;

/// Errors from loading weights or encoding peptides.
#[derive(Debug, thiserror::Error)]
pub enum RtError {
    /// I/O failure reading the weight file.
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
    /// Malformed or inconsistent weight file.
    #[error("weights: {0}")]
    Format(String),
}

/// A peptidoform: residues plus Unimod modifications by position
/// (0 = N-terminus, 1..=n residues, n + 1 = C-terminus).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct RtPeptide {
    /// One-letter residues.
    pub sequence: String,
    /// (position, Unimod id), sorted.
    pub mods: Vec<(usize, u32)>,
}

impl RtPeptide {
    /// Parse ProForma with Unimod accessions, e.g. `[UNIMOD:737]-PEPS[UNIMOD:21]K`.
    pub fn from_proforma(s: &str) -> Option<Self> {
        let b = s.as_bytes();
        let mut mods = Vec::new();
        let mut seq = String::new();
        let mut i = 0;
        let read_mod = |i: usize| -> Option<(u32, usize)> {
            let rest = s.get(i..)?;
            let rest = rest.strip_prefix("[UNIMOD:")?;
            let end = rest.find(']')?;
            Some((rest[..end].parse().ok()?, i + 8 + end + 1))
        };
        // N-terminal mods: one or more [..] followed by '-'
        if b.first() == Some(&b'[') {
            let mut ids = Vec::new();
            while b.get(i) == Some(&b'[') {
                let (id, j) = read_mod(i)?;
                ids.push(id);
                i = j;
            }
            if b.get(i) != Some(&b'-') {
                return None;
            }
            i += 1;
            mods.extend(ids.into_iter().map(|id| (0, id)));
        }
        while i < b.len() {
            match b[i] {
                b'[' => {
                    let (id, j) = read_mod(i)?;
                    if seq.is_empty() {
                        return None;
                    }
                    mods.push((seq.len(), id));
                    i = j;
                }
                b'-' => {
                    i += 1;
                    while i < b.len() {
                        let (id, j) = read_mod(i)?;
                        mods.push((seq.len() + 1, id));
                        i = j;
                    }
                }
                c if c.is_ascii_uppercase() => {
                    seq.push(c as char);
                    i += 1;
                }
                _ => return None,
            }
        }
        mods.sort_unstable();
        (!seq.is_empty()).then_some(Self {
            sequence: seq,
            mods,
        })
    }

    /// Parse an andes PIN peptide (`K.M+42.01056+15.99491PEPS+79.96633K.D`).
    ///
    /// Mass deltas follow their residue; N-terminal deltas follow the first
    /// residue, and the first one that is N-term-capable goes to the N-terminus
    /// (`K+42.01` reads as N-terminal acetyl, the common case). `None` when a
    /// delta matches no Unimod entry within 2 mDa.
    pub fn from_andes(s: &str) -> Option<Self> {
        let b = s.as_bytes();
        let core = if b.len() > 4 && b[1] == b'.' && b[b.len() - 2] == b'.' {
            &s[2..s.len() - 2]
        } else {
            s
        };
        let c = core.as_bytes();
        let mut seq = String::new();
        let mut mods: Vec<(usize, u32)> = Vec::new();
        let mut nterm_used = false;
        let mut i = 0;
        while i < c.len() {
            if c[i].is_ascii_uppercase() {
                seq.push(c[i] as char);
                i += 1;
                continue;
            }
            if c[i] != b'+' && c[i] != b'-' {
                return None;
            }
            let st = i;
            i += 1;
            while i < c.len() && (c[i].is_ascii_digit() || c[i] == b'.') {
                i += 1;
            }
            let d: f64 = core[st..i].parse().ok()?;
            let pos = seq.len();
            if pos == 0 {
                return None;
            }
            let aa = seq.as_bytes()[pos - 1];
            let res_id = mods::match_mass(d, Some(aa), false, 2e-3);
            let nt_id = if pos == 1 && !nterm_used {
                mods::match_mass(d, None, true, 2e-3)
            } else {
                None
            };
            match (nt_id, res_id) {
                (Some(nt), _) => {
                    mods.push((0, nt));
                    nterm_used = true;
                }
                (None, Some(r)) => mods.push((pos, r)),
                (None, None) => return None,
            }
        }
        mods.sort_unstable();
        (!seq.is_empty()).then_some(Self {
            sequence: seq,
            mods,
        })
    }

    /// ProForma with Unimod accessions (inverse of [`from_proforma`](Self::from_proforma)).
    pub fn to_proforma(&self) -> String {
        let n = self.sequence.len();
        let tag = |pos: usize| -> String {
            self.mods
                .iter()
                .filter(|(p, _)| *p == pos)
                .map(|(_, id)| format!("[UNIMOD:{id}]"))
                .collect()
        };
        let mut out = tag(0);
        if !out.is_empty() {
            out.push('-');
        }
        for (i, c) in self.sequence.chars().enumerate() {
            out.push(c);
            out.push_str(&tag(i + 1));
        }
        let ct = tag(n + 1);
        if !ct.is_empty() {
            out.push('-');
            out.push_str(&ct);
        }
        out
    }

    /// Tokens and scaled composition rows, or `None` if unsupported.
    pub fn encode(&self) -> Option<(Vec<usize>, Vec<[f32; N_COMP]>)> {
        let n = self.sequence.len();
        if n == 0 || n > MAX_LEN {
            return None;
        }
        let mut tok = Vec::with_capacity(n + 2);
        tok.push(1);
        for c in self.sequence.bytes() {
            // I is encoded as L: the training sequences carry no I.
            let c = if c == b'I' { b'L' } else { c };
            tok.push(3 + AA.iter().position(|&a| a == c)?);
        }
        tok.push(2);
        let mut comp = vec![[0.0f32; N_COMP]; n + 2];
        for &(pos, id) in &self.mods {
            if pos > n + 1 {
                return None;
            }
            let f = mods::features(mods::by_id(id)?);
            for (c, v) in comp[pos].iter_mut().zip(f) {
                *c += v;
            }
        }
        for row in &mut comp {
            for (c, k) in row.iter_mut().zip(COMP_SCALE) {
                *c *= k;
            }
        }
        Some((tok, comp))
    }
}

/// Conv layer with weights stored `[k][in][out]` so the inner loop runs over
/// contiguous output channels.
struct Conv {
    k: usize,
    cin: usize,
    cout: usize,
    dil: usize,
    w: Vec<f32>,
    b: Vec<f32>,
}

impl Conv {
    fn from_torch(w: &Tensor, b: &Tensor, dil: usize) -> Result<Self, RtError> {
        let [cout, cin, k] = w.dims[..] else {
            return Err(RtError::Format("conv weight must be 3-D".into()));
        };
        if b.dims != [cout] {
            return Err(RtError::Format("conv bias shape".into()));
        }
        let mut t = vec![0.0; k * cin * cout];
        for o in 0..cout {
            for i in 0..cin {
                for kk in 0..k {
                    t[(kk * cin + i) * cout + o] = w.data[(o * cin + i) * k + kk];
                }
            }
        }
        Ok(Self {
            k,
            cin,
            cout,
            dil,
            w: t,
            b: b.data.clone(),
        })
    }

    /// `y[p] = b + Σ_k Σ_i w[k][i] · x[p + (k - k/2)·dil]` with zero padding.
    fn apply(&self, x: &[f32], n: usize, y: &mut [f32]) {
        let (cin, cout) = (self.cin, self.cout);
        for p in 0..n {
            let yp = &mut y[p * cout..(p + 1) * cout];
            yp.copy_from_slice(&self.b);
            for kk in 0..self.k {
                let q = p as isize + (kk as isize - (self.k / 2) as isize) * self.dil as isize;
                if q < 0 || q >= n as isize {
                    continue;
                }
                let xq = &x[q as usize * cin..(q as usize + 1) * cin];
                let wk = &self.w[kk * cin * cout..(kk + 1) * cin * cout];
                for (i, &xv) in xq.iter().enumerate() {
                    let wi = &wk[i * cout..(i + 1) * cout];
                    for (yo, &wo) in yp.iter_mut().zip(wi) {
                        *yo += wo * xv;
                    }
                }
            }
        }
    }
}

struct Tensor {
    dims: Vec<usize>,
    data: Vec<f32>,
}

#[inline]
fn silu(v: f32) -> f32 {
    v / (1.0 + (-v).exp())
}

/// The andes-rt network.
pub struct RtModel {
    emb: Vec<f32>,
    emb_dim: usize,
    ch: usize,
    inp: Conv,
    blocks: Vec<(Conv, Conv)>,
    att_w: Vec<f32>,
    att_b: f32,
    fc1_w: Vec<f32>,
    fc1_b: Vec<f32>,
    fc2_w: Vec<f32>,
    fc2_b: f32,
}

/// Reusable per-thread buffers for [`RtModel::predict_with`].
#[derive(Default)]
pub struct Scratch {
    input: Vec<f32>,
    x: Vec<f32>,
    h: Vec<f32>,
    z: Vec<f32>,
}

impl RtModel {
    /// Load an `ARTN0001` weight file written by `rt_model.export_weights`.
    pub fn load(path: impl AsRef<Path>) -> Result<Self, RtError> {
        let mut buf = Vec::new();
        std::fs::File::open(path)?.read_to_end(&mut buf)?;
        Self::from_bytes(&buf)
    }

    /// Parse an `ARTN0001` blob.
    pub fn from_bytes(buf: &[u8]) -> Result<Self, RtError> {
        let bad = |m: &str| RtError::Format(m.to_string());
        if buf.len() < 12 || &buf[..8] != b"ARTN0001" {
            return Err(bad("missing ARTN0001 magic"));
        }
        let mut off = 8;
        let u32_at = |off: &mut usize| -> Result<usize, RtError> {
            let v = buf.get(*off..*off + 4).ok_or_else(|| bad("truncated"))?;
            *off += 4;
            Ok(u32::from_le_bytes(v.try_into().unwrap()) as usize)
        };
        let n = u32_at(&mut off)?;
        let mut t: HashMap<String, Tensor> = HashMap::new();
        for _ in 0..n {
            let len = u32_at(&mut off)?;
            let name = std::str::from_utf8(
                buf.get(off..off + len)
                    .ok_or_else(|| bad("truncated name"))?,
            )
            .map_err(|_| bad("name not UTF-8"))?
            .to_string();
            off += len;
            let nd = u32_at(&mut off)?;
            let mut dims = Vec::with_capacity(nd);
            for _ in 0..nd {
                dims.push(u32_at(&mut off)?);
            }
            let count: usize = dims.iter().product();
            let raw = buf
                .get(off..off + 4 * count)
                .ok_or_else(|| bad("truncated data"))?;
            off += 4 * count;
            let data = raw
                .chunks_exact(4)
                .map(|c| f32::from_le_bytes(c.try_into().unwrap()))
                .collect();
            t.insert(name, Tensor { dims, data });
        }
        if off != buf.len() {
            return Err(bad("trailing bytes"));
        }
        let get = |k: &str| {
            t.get(k)
                .ok_or_else(|| RtError::Format(format!("missing tensor {k}")))
        };
        let dil: Vec<usize> = get("dilations")?.data.iter().map(|&d| d as usize).collect();
        let emb = get("emb.weight")?;
        let [ntok, emb_dim] = emb.dims[..] else {
            return Err(bad("emb shape"));
        };
        if ntok != N_TOK {
            return Err(bad("vocabulary size"));
        }
        let inp = Conv::from_torch(get("inp.weight")?, get("inp.bias")?, 1)?;
        if inp.cin != emb_dim + N_COMP || inp.k != 1 {
            return Err(bad("inp shape"));
        }
        let ch = inp.cout;
        let mut blocks = Vec::new();
        for (i, &d) in dil.iter().enumerate() {
            let c = Conv::from_torch(
                get(&format!("convs.{i}.weight"))?,
                get(&format!("convs.{i}.bias"))?,
                d,
            )?;
            let m = Conv::from_torch(
                get(&format!("mix.{i}.weight"))?,
                get(&format!("mix.{i}.bias"))?,
                1,
            )?;
            if c.cin != ch || c.cout != ch || m.cin != ch || m.cout != ch || c.k % 2 != 1 {
                return Err(bad("block shape"));
            }
            blocks.push((c, m));
        }
        let att = get("att.weight")?;
        let fc1 = get("fc1.weight")?;
        let fc2 = get("fc2.weight")?;
        if att.dims != [1, ch, 1] || fc1.dims != [ch, 2 * ch] || fc2.dims != [1, ch] {
            return Err(bad("head shapes"));
        }
        Ok(Self {
            emb: emb.data.clone(),
            emb_dim,
            ch,
            inp,
            blocks,
            att_w: att.data.clone(),
            att_b: get("att.bias")?.data[0],
            fc1_w: fc1.data.clone(),
            fc1_b: get("fc1.bias")?.data.clone(),
            fc2_w: fc2.data.clone(),
            fc2_b: get("fc2.bias")?.data[0],
        })
    }

    /// Number of parameters.
    pub fn n_params(&self) -> usize {
        let conv = |c: &Conv| c.w.len() + c.b.len();
        self.emb.len()
            + conv(&self.inp)
            + self
                .blocks
                .iter()
                .map(|(a, b)| conv(a) + conv(b))
                .sum::<usize>()
            + self.att_w.len()
            + 1
            + self.fc1_w.len()
            + self.fc1_b.len()
            + self.fc2_w.len()
            + 1
    }

    /// Retention index for one peptidoform; `None` if it cannot be encoded.
    pub fn predict(&self, p: &RtPeptide) -> Option<f32> {
        self.predict_with(p, &mut Scratch::default())
    }

    /// As [`predict`](Self::predict), reusing `scratch` across calls.
    pub fn predict_with(&self, p: &RtPeptide, s: &mut Scratch) -> Option<f32> {
        let (tok, comp) = p.encode()?;
        let n = tok.len();
        let (ed, ch) = (self.emb_dim, self.ch);
        let din = ed + N_COMP;
        s.input.clear();
        for (t, c) in tok.iter().zip(&comp) {
            s.input.extend_from_slice(&self.emb[t * ed..(t + 1) * ed]);
            s.input.extend_from_slice(c);
        }
        debug_assert_eq!(s.input.len(), n * din);
        s.x.resize(n * ch, 0.0);
        s.h.resize(n * ch, 0.0);
        s.z.resize(n * ch, 0.0);
        self.inp.apply(&s.input, n, &mut s.x);
        for (conv, mix) in &self.blocks {
            conv.apply(&s.x, n, &mut s.h);
            for v in s.h.iter_mut() {
                *v = silu(*v);
            }
            mix.apply(&s.h, n, &mut s.z);
            for (x, z) in s.x.iter_mut().zip(&s.z) {
                *x += z;
            }
        }
        // attention + mean pooling
        let mut logits = Vec::with_capacity(n);
        for p in 0..n {
            let xp = &s.x[p * ch..(p + 1) * ch];
            logits.push(self.att_b + xp.iter().zip(&self.att_w).map(|(a, b)| a * b).sum::<f32>());
        }
        let mx = logits.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
        let mut wsum = 0.0f32;
        for l in &mut logits {
            *l = (*l - mx).exp();
            wsum += *l;
        }
        let mut pooled = vec![0.0f32; 2 * ch];
        for (p, l) in logits.iter().enumerate() {
            let a = l / wsum;
            let xp = &s.x[p * ch..(p + 1) * ch];
            for c in 0..ch {
                pooled[c] += a * xp[c];
                pooled[ch + c] += xp[c] / n as f32;
            }
        }
        let mut out = self.fc2_b;
        for o in 0..ch {
            let row = &self.fc1_w[o * 2 * ch..(o + 1) * 2 * ch];
            let h = self.fc1_b[o] + row.iter().zip(&pooled).map(|(a, b)| a * b).sum::<f32>();
            out += self.fc2_w[o] * silu(h);
        }
        Some(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_proforma() {
        let p = RtPeptide::from_proforma("[UNIMOD:737]-PEPS[UNIMOD:21]K[UNIMOD:737]").unwrap();
        assert_eq!(p.sequence, "PEPSK");
        assert_eq!(p.mods, vec![(0, 737), (4, 21), (5, 737)]);
        let p = RtPeptide::from_proforma("AC[UNIMOD:4]K-[UNIMOD:2]").unwrap();
        assert_eq!(p.mods, vec![(2, 4), (4, 2)]);
        assert!(RtPeptide::from_proforma("PEP[Oxidation]K").is_none());
        for s in [
            "[UNIMOD:737]-PEPS[UNIMOD:21]K[UNIMOD:737]",
            "AC[UNIMOD:4]K-[UNIMOD:2]",
            "PEPTIDE",
        ] {
            assert_eq!(RtPeptide::from_proforma(s).unwrap().to_proforma(), s);
        }
    }

    #[test]
    fn parses_andes() {
        let p = RtPeptide::from_andes("K.M+42.01056+15.99491PEPS+79.96633K.D").unwrap();
        assert_eq!(
            (p.sequence.as_str(), p.mods.clone()),
            ("MPEPSK", vec![(0, 1), (1, 35), (5, 21)])
        );
        let p = RtPeptide::from_andes("K.K+229.16293+229.16293PEPK+229.16293.D").unwrap();
        assert_eq!(p.mods, vec![(0, 737), (1, 737), (5, 737)]);
        let p = RtPeptide::from_andes("-.K+42.01056PEPK.-").unwrap();
        assert_eq!(p.mods, vec![(0, 1)]);
        let p = RtPeptide::from_andes("K.Q-17.02655PEPK.D").unwrap();
        assert_eq!(p.mods, vec![(0, 28)]);
        assert!(RtPeptide::from_andes("K.PEPXK.D")
            .unwrap()
            .encode()
            .is_none());
        assert!(RtPeptide::from_andes("K.PEP+1.23456K.D").is_none());
    }

    #[test]
    fn tmtpro_is_not_itraq8plex() {
        let p = RtPeptide::from_andes("K.P+304.20715EPTIDEK+304.20715.R").unwrap();
        assert_eq!(p.to_proforma(), "[UNIMOD:2016]-PEPTIDEK[UNIMOD:2016]");
    }

    #[test]
    fn encodes_composition_on_terminal_tokens() {
        let p = RtPeptide::from_proforma("[UNIMOD:1]-AK").unwrap();
        let (tok, comp) = p.encode().unwrap();
        assert_eq!(tok, vec![1, 3, 3 + 8, 2]);
        assert!((comp[0][0] - 0.2).abs() < 1e-6 && (comp[0][7] - 0.42010565).abs() < 1e-6);
        assert!(comp[1].iter().all(|v| *v == 0.0));
        let i = RtPeptide::from_proforma("PIK").unwrap().encode().unwrap().0;
        let l = RtPeptide::from_proforma("PLK").unwrap().encode().unwrap().0;
        assert_eq!(i, l);
    }
}
