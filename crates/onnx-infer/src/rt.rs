pub mod mods;

use std::collections::HashMap;
use std::io::Read;
use std::ops::Range;
use std::path::Path;

use mods::N_COMP;

pub const AA: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
pub const MAX_LEN: usize = 50;
pub const COMP_SCALE: [f32; N_COMP] = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 1.0];
const N_TOK: usize = 23;

#[derive(Debug, thiserror::Error)]
pub enum RtError {
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
    #[error("weights: {0}")]
    Format(String),
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct RtPeptide {
    pub sequence: String,
    pub mods: Vec<(usize, u32)>,
}

impl RtPeptide {
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
                (None, None) if pos == 1 && !nterm_used && mods.is_empty() => {
                    let (nt, r) = mods::split_nterm_residue(d, aa, 2e-3)?;
                    mods.push((0, nt));
                    mods.push((1, r));
                    nterm_used = true;
                }
                (None, None) => return None,
            }
        }
        mods.sort_unstable();
        (!seq.is_empty()).then_some(Self {
            sequence: seq,
            mods,
        })
    }

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

    pub fn encode(&self) -> Option<(Vec<usize>, Vec<[f32; N_COMP]>)> {
        let n = self.sequence.len();
        if n == 0 || n > MAX_LEN {
            return None;
        }
        let mut tok = Vec::with_capacity(n + 2);
        tok.push(1);
        for c in self.sequence.bytes() {
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

pub(crate) struct Conv {
    pub(crate) k: usize,
    pub(crate) cin: usize,
    pub(crate) cout: usize,
    dil: usize,
    pub(crate) w: Weights,
    pub(crate) b: Vec<f32>,
}

/// Conv weights in 32-byte aligned storage: with a width that is a multiple of 8, every
/// 8-channel group of a weight row starts on a 32-byte boundary, so the AVX kernel loads it
/// with one aligned load.
#[derive(Default)]
pub(crate) struct Weights {
    lanes: Vec<Lane>,
    len: usize,
}

#[derive(Clone, Copy, Default)]
#[repr(C, align(32))]
struct Lane([f32; 8]);

impl FromIterator<f32> for Weights {
    fn from_iter<I: IntoIterator<Item = f32>>(values: I) -> Self {
        let mut w = Self::default();
        for v in values {
            if w.len % 8 == 0 {
                w.lanes.push(Lane::default());
            }
            w.lanes[w.len / 8].0[w.len % 8] = v;
            w.len += 1;
        }
        w
    }
}

impl std::ops::Deref for Weights {
    type Target = [f32];

    fn deref(&self) -> &[f32] {
        // SAFETY: `Lane` is `repr(C)` around `[f32; 8]` with no padding, so the lanes are
        // `8 * lanes.len()` contiguous initialised `f32`, and `len` never exceeds that.
        unsafe { std::slice::from_raw_parts(self.lanes.as_ptr().cast::<f32>(), self.len) }
    }
}

impl Conv {
    pub(crate) fn from_torch(w: &Tensor, b: &Tensor, dil: usize) -> Result<Self, RtError> {
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
            w: t.into_iter().collect(),
            b: b.data.clone(),
        })
    }

    pub(crate) fn apply(&self, x: &[f32], n: usize, y: &mut [f32]) {
        self.apply_rows(x, n, 0..n, y);
    }

    /// The output rows `rows` of [`Self::apply`] over the `n` rows of `x`, written to the first
    /// `rows.len()` rows of `y`. Each output sums its bias, then the taps in order, then the
    /// input channels in order, whatever rows are asked for.
    pub(crate) fn apply_rows(&self, x: &[f32], n: usize, rows: Range<usize>, y: &mut [f32]) {
        let (cin, cout) = (self.cin, self.cout);
        let half = (self.k / 2) as isize;
        let dil = self.dil as isize;
        let tiled = cout - cout % TILE_O;
        let aligned = cout % 8 == 0 && (self.w.as_ptr() as usize) % 32 == 0;
        let mut p0 = rows.start;
        while p0 < rows.end {
            let np = TILE_P.min(rows.end - p0);
            let y0 = p0 - rows.start;
            for o0 in (0..tiled).step_by(TILE_O) {
                let mut acc = [[0.0f32; TILE_O]; TILE_P];
                for a in acc.iter_mut().take(np) {
                    a.copy_from_slice(&self.b[o0..o0 + TILE_O]);
                }
                for kk in 0..self.k {
                    let off = (kk as isize - half) * dil;
                    let lo = (-(p0 as isize) - off).clamp(0, np as isize) as usize;
                    let hi = (n as isize - p0 as isize - off).clamp(0, np as isize) as usize;
                    if lo >= hi {
                        continue;
                    }
                    let q0 = ((p0 + lo) as isize + off) as usize;
                    let xs = &x[q0 * cin..(q0 + hi - lo) * cin];
                    let wk = &self.w[kk * cin * cout..(kk + 1) * cin * cout];
                    tap_rows(&mut acc[lo..hi], xs, wk, (cin, cout, o0), aligned);
                }
                for (r, a) in acc.iter().take(np).enumerate() {
                    y[(y0 + r) * cout + o0..][..TILE_O].copy_from_slice(a);
                }
            }
            if tiled < cout {
                for p in p0..p0 + np {
                    let yr = p - rows.start;
                    let yp = &mut y[yr * cout + tiled..(yr + 1) * cout];
                    yp.copy_from_slice(&self.b[tiled..]);
                    for kk in 0..self.k {
                        let q = p as isize + (kk as isize - half) * dil;
                        if q < 0 || q >= n as isize {
                            continue;
                        }
                        let xq = &x[q as usize * cin..(q as usize + 1) * cin];
                        let wk = &self.w[kk * cin * cout..(kk + 1) * cin * cout];
                        for (i, &xv) in xq.iter().enumerate() {
                            let wi = &wk[i * cout + tiled..(i + 1) * cout];
                            for (yo, &wo) in yp.iter_mut().zip(wi) {
                                *yo += wo * xv;
                            }
                        }
                    }
                }
            }
            p0 += np;
        }
    }
}

const TILE_P: usize = 6;
const TILE_O: usize = 16;

/// `(cin, cout, o0)` of a tap: the input width, the output width and the first output channel
/// of the tile.
type TapShape = (usize, usize, usize);

fn tap_rows(acc: &mut [[f32; TILE_O]], xs: &[f32], wk: &[f32], shape: TapShape, aligned: bool) {
    #[cfg(all(target_arch = "x86_64", target_feature = "avx"))]
    use avx::tap;
    match acc.len() {
        1 => tap::<1>(acc, xs, wk, shape, aligned),
        2 => tap::<2>(acc, xs, wk, shape, aligned),
        3 => tap::<3>(acc, xs, wk, shape, aligned),
        4 => tap::<4>(acc, xs, wk, shape, aligned),
        5 => tap::<5>(acc, xs, wk, shape, aligned),
        6 => tap::<6>(acc, xs, wk, shape, aligned),
        _ => unreachable!("tile holds at most TILE_P rows"),
    }
}

/// Adds one tap to the `R` rows of `acc`: `xs` holds the `R` input rows and `wk` the tap's
/// weights, `cin` rows of `cout`. `aligned` only selects the load instruction of the AVX build.
#[cfg_attr(
    all(target_arch = "x86_64", target_feature = "avx", not(test)),
    allow(dead_code)
)]
#[inline(always)]
fn tap<const R: usize>(
    acc: &mut [[f32; TILE_O]],
    xs: &[f32],
    wk: &[f32],
    (cin, cout, o0): TapShape,
    _aligned: bool,
) {
    let Some(acc) = acc.first_chunk_mut::<R>() else {
        unreachable!("the caller passes R rows")
    };
    let mut a = *acc;
    let xs = &xs[..R * cin];
    for i in 0..cin {
        let w = &wk[i * cout + o0..][..TILE_O];
        for (r, ar) in a.iter_mut().enumerate() {
            let xv = xs[r * cin + i];
            for (ao, &wo) in ar.iter_mut().zip(w) {
                *ao += wo * xv;
            }
        }
    }
    *acc = a;
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx"))]
mod avx {
    use super::{TapShape, TILE_O};
    use std::arch::x86_64::{
        __m256, _mm256_add_ps, _mm256_load_ps, _mm256_loadu_ps, _mm256_mul_ps, _mm256_set1_ps,
        _mm256_storeu_ps,
    };

    const V: usize = TILE_O / 8;
    const _: () = assert!(TILE_O == 16, "the AVX tap holds 2 × 8 channels per row");

    #[inline(always)]
    unsafe fn load(row: &[f32; TILE_O]) -> [__m256; V] {
        let p = row.as_ptr();
        [_mm256_loadu_ps(p), _mm256_loadu_ps(p.add(8))]
    }

    #[inline(always)]
    unsafe fn store(row: &mut [f32; TILE_O], a: [__m256; V]) {
        let p = row.as_mut_ptr();
        _mm256_storeu_ps(p, a[0]);
        _mm256_storeu_ps(p.add(8), a[1]);
    }

    #[inline(always)]
    unsafe fn step(a: [__m256; V], w: [__m256; V], x: f32) -> [__m256; V] {
        let b = _mm256_set1_ps(x);
        [
            _mm256_add_ps(a[0], _mm256_mul_ps(w[0], b)),
            _mm256_add_ps(a[1], _mm256_mul_ps(w[1], b)),
        ]
    }

    /// The AVX form of [`super::tap`], with the same sums in the same order.
    #[inline(always)]
    pub(super) fn tap<const R: usize>(
        acc: &mut [[f32; TILE_O]],
        xs: &[f32],
        wk: &[f32],
        shape: TapShape,
        aligned: bool,
    ) {
        if aligned {
            rows::<R, true>(acc, xs, wk, shape)
        } else {
            rows::<R, false>(acc, xs, wk, shape)
        }
    }

    #[inline(always)]
    fn rows<const R: usize, const ALIGNED: bool>(
        acc: &mut [[f32; TILE_O]],
        xs: &[f32],
        wk: &[f32],
        (cin, cout, o0): TapShape,
    ) {
        let Some(acc) = acc.first_chunk_mut::<R>() else {
            unreachable!("the caller passes R rows")
        };
        if cin == 0 {
            return;
        }
        assert!(xs.len() >= R * cin && wk.len() >= (cin - 1) * cout + o0 + TILE_O);
        assert!(!ALIGNED || (cout % 8 == 0 && (wk.as_ptr() as usize + 4 * o0) % 32 == 0));
        let (xp, wp) = (xs.as_ptr(), wk.as_ptr());
        // SAFETY: AVX is enabled for this build (module `cfg`). The asserts above keep every
        // weight row `i * cout + o0 .. + TILE_O` (i < cin) and every input `r * cin + i`
        // (r < R) in bounds, and make each aligned load start on a 32-byte boundary.
        unsafe {
            let mut a: [[__m256; V]; R] = std::array::from_fn(|r| load(&acc[r]));
            for i in 0..cin {
                let row = wp.add(i * cout + o0);
                let w = if ALIGNED {
                    [_mm256_load_ps(row), _mm256_load_ps(row.add(8))]
                } else {
                    [_mm256_loadu_ps(row), _mm256_loadu_ps(row.add(8))]
                };
                for (r, ar) in a.iter_mut().enumerate() {
                    *ar = step(*ar, w, *xp.add(r * cin + i));
                }
            }
            for (row, ar) in acc.iter_mut().zip(a) {
                store(row, ar);
            }
        }
    }
}

pub(crate) struct Tensor {
    pub(crate) dims: Vec<usize>,
    pub(crate) data: Vec<f32>,
}

pub(crate) fn read_tensors(
    buf: &[u8],
    magic: &[u8; 8],
) -> Result<HashMap<String, Tensor>, RtError> {
    let bad = |m: &str| RtError::Format(m.to_string());
    if buf.len() < 12 || &buf[..8] != magic {
        return Err(RtError::Format(format!(
            "missing {} magic",
            String::from_utf8_lossy(magic)
        )));
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
        let bytes = crate::tensor::element_count(&dims)
            .and_then(|count| count.checked_mul(4))
            .ok_or_else(|| bad("tensor size overflows usize"))?;
        let raw = buf
            .get(off..off.saturating_add(bytes))
            .ok_or_else(|| bad("truncated data"))?;
        off += bytes;
        let data = raw
            .chunks_exact(4)
            .map(|c| f32::from_le_bytes(c.try_into().unwrap()))
            .collect();
        t.insert(name, Tensor { dims, data });
    }
    if off != buf.len() {
        return Err(bad("trailing bytes"));
    }
    Ok(t)
}

#[inline]
pub(crate) fn silu(v: f32) -> f32 {
    v / (1.0 + (-v).exp())
}

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

#[derive(Default)]
pub struct Scratch {
    input: Vec<f32>,
    x: Vec<f32>,
    h: Vec<f32>,
    z: Vec<f32>,
}

impl RtModel {
    pub fn load(path: impl AsRef<Path>) -> Result<Self, RtError> {
        let mut buf = Vec::new();
        std::fs::File::open(path)?.read_to_end(&mut buf)?;
        Self::from_bytes(&buf)
    }

    pub fn from_bytes(buf: &[u8]) -> Result<Self, RtError> {
        let bad = |m: &str| RtError::Format(m.to_string());
        let t = read_tensors(buf, b"ARTN0001")?;
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

    pub fn predict(&self, p: &RtPeptide) -> Option<f32> {
        self.predict_with(p, &mut Scratch::default())
    }

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
    fn conv_apply_matches_reference_bitwise() {
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut rnd = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            ((state >> 40) as f32 / (1u64 << 24) as f32) * 2.0 - 1.0
        };
        for &(k, cin, cout, dil) in &[
            (1, 40, 96, 1),
            (3, 96, 96, 8),
            (2, 96, 96, 1),
            (3, 64, 64, 4),
            (1, 96, 12, 1),
            (3, 8, 40, 2),
            (5, 7, 33, 3),
        ] {
            let conv = Conv {
                k,
                cin,
                cout,
                dil,
                w: (0..k * cin * cout).map(|_| rnd()).collect(),
                b: (0..cout).map(|_| rnd()).collect(),
            };
            for n in 1..=20 {
                let x: Vec<f32> = (0..n * cin).map(|_| rnd()).collect();
                let mut y = vec![f32::NAN; n * cout];
                conv.apply(&x, n, &mut y);
                for p in 0..n {
                    for (o, &bias) in conv.b.iter().enumerate() {
                        let mut want = bias;
                        for kk in 0..k {
                            let q = p as isize + (kk as isize - (k / 2) as isize) * dil as isize;
                            if q < 0 || q >= n as isize {
                                continue;
                            }
                            for i in 0..cin {
                                let wv = conv.w[(kk * cin + i) * cout + o];
                                want += wv * x[q as usize * cin + i];
                            }
                        }
                        assert_eq!(
                            y[p * cout + o].to_bits(),
                            want.to_bits(),
                            "k{k} cin{cin} cout{cout} dil{dil} n{n} p{p} o{o}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn conv_rows_match_the_full_apply_bitwise() {
        let mut state = 0x6a09_e667_f3bc_c908u64;
        let mut rnd = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            ((state >> 40) as f32 / (1u64 << 24) as f32) * 2.0 - 1.0
        };
        for &(k, cin, cout, dil) in &[
            (2, 96, 96, 1),
            (3, 64, 64, 8),
            (1, 96, 12, 1),
            (3, 7, 33, 2),
        ] {
            let conv = Conv {
                k,
                cin,
                cout,
                dil,
                w: (0..k * cin * cout).map(|_| rnd()).collect(),
                b: (0..cout).map(|_| rnd()).collect(),
            };
            for n in 1..=24 {
                let x: Vec<f32> = (0..n * cin).map(|_| rnd()).collect();
                let mut full = vec![f32::NAN; n * cout];
                conv.apply(&x, n, &mut full);
                for (lo, hi) in [(0, n), (n.min(2), n), (n / 3, n - n / 4), (n - 1, n)] {
                    let mut part = vec![f32::NAN; (hi - lo) * cout];
                    conv.apply_rows(&x, n, lo..hi, &mut part);
                    let bits = |v: &[f32]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
                    assert_eq!(
                        bits(&part),
                        bits(&full[lo * cout..hi * cout]),
                        "k{k} n{n} {lo}..{hi}"
                    );
                }
            }
        }
    }

    #[test]
    fn conv_weights_are_32_byte_aligned() {
        for len in [1usize, 7, 8, 9, 96 * 96 * 3] {
            let w: Weights = (0..len).map(|i| i as f32).collect();
            assert_eq!(w.as_ptr() as usize % 32, 0);
            assert_eq!(w.len(), len);
            assert!(w.iter().enumerate().all(|(i, &v)| v == i as f32));
        }
    }

    #[cfg(all(target_arch = "x86_64", target_feature = "avx"))]
    #[test]
    fn avx_tap_matches_portable_tap() {
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut rnd = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            ((state >> 40) as f32 / (1u64 << 24) as f32) * 2.0 - 1.0
        };
        let bits =
            |a: &[[f32; TILE_O]]| a.iter().flatten().map(|v| v.to_bits()).collect::<Vec<_>>();
        for (cin, cout, o0) in [(37, 96, 32), (96, 96, 80), (5, 40, 16), (8, 33, 16)] {
            let wk: Weights = (0..cin * cout).map(|_| rnd()).collect();
            let xs: Vec<f32> = (0..TILE_P * cin).map(|_| rnd()).collect();
            let mut init = [[0.0f32; TILE_O]; TILE_P];
            for v in init.iter_mut().flatten() {
                *v = rnd();
            }
            let aligned = cout % 8 == 0;
            for r in 1..=TILE_P {
                let (mut p, mut v) = (init, init);
                tap_rows_portable(&mut p[..r], &xs, &wk, (cin, cout, o0));
                tap_rows(&mut v[..r], &xs, &wk, (cin, cout, o0), aligned);
                assert_eq!(bits(&p), bits(&v), "cin{cin} cout{cout} o0{o0} rows{r}");
            }
        }
    }

    #[cfg(all(target_arch = "x86_64", target_feature = "avx"))]
    fn tap_rows_portable(acc: &mut [[f32; TILE_O]], xs: &[f32], wk: &[f32], shape: TapShape) {
        match acc.len() {
            1 => tap::<1>(acc, xs, wk, shape, false),
            2 => tap::<2>(acc, xs, wk, shape, false),
            3 => tap::<3>(acc, xs, wk, shape, false),
            4 => tap::<4>(acc, xs, wk, shape, false),
            5 => tap::<5>(acc, xs, wk, shape, false),
            _ => tap::<6>(acc, xs, wk, shape, false),
        }
    }

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
