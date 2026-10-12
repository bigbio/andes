"""andes-rt v1: peptidoform encoding, the CNN, and the weight export format.

Encoding (mirrored by ``onnx-infer::rt_native``):
- tokens: [N-term, residues..., C-term], PAD = 0, N-term = 1, C-term = 2,
  residue = 3 + index in ``AA``; at most ``MAX_LEN`` residues. I is encoded
  as L (isobaric; the training sequences carry no I).
- per position a composition vector ``[C, H, N, O, S, P, other, mass/100]``
  (``mod_composition``), summed over the mods at that position and multiplied
  by ``COMP_SCALE``. N-terminal mods sit on the N-term token, C-terminal mods
  on the C-term token.

Network: embedding(32) ++ composition(8) -> 1x1 conv to 64 channels -> four
residual blocks ``x + mix(silu(conv_k3_dilated(x)))`` with dilations 1, 2, 4, 8
-> attention pooling ++ mean pooling -> linear(128, 64) -> silu -> linear(64, 1).
Padded positions are zeroed after every layer, so a peptide's output does not
depend on what it was batched with.

Export (``ARTN0001``): magic, u32 tensor count, then per tensor u32 name length,
UTF-8 name, u32 ndim, u32 dims, little-endian f32 data.
"""
from __future__ import annotations

import os
import re
import struct
import sys

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mod_composition import N_COMP, UNIMOD, VECTORS, match_mass  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_TOK = {a: i + 3 for i, a in enumerate(AA)}
N_TOK = 3 + len(AA)
MAX_LEN = 50
COMP_SCALE = np.array([0.1] * 7 + [1.0], dtype=np.float32)
DILATIONS = (1, 2, 4, 8)

_PF = re.compile(r"\[UNIMOD:(\d+)\]")


def parse_proforma(s: str):
    """``[UNIMOD:737]-PEPS[UNIMOD:21]K`` -> (sequence, [(pos, unimod id)]); pos 0 = N-term."""
    mods = []
    i = 0
    if s.startswith("["):
        j = s.index("]-")
        mods += [(0, int(m)) for m in _PF.findall(s[: j + 1])]
        i = j + 2
    seq = []
    while i < len(s):
        c = s[i]
        if c == "-" and s[i + 1] == "[":
            mods += [(len(seq) + 1, int(m)) for m in _PF.findall(s[i + 1:])]
            break
        if c == "[":
            j = s.index("]", i)
            mods.append((len(seq), int(_PF.match(s[i: j + 1]).group(1))))
            i = j + 1
            continue
        seq.append(c)
        i += 1
    return "".join(seq), mods


def parse_andes(s: str):
    """andes PIN peptide ``K.M+15.99491PEPS+79.96633K.D`` -> (sequence, [(pos, unimod id)]).

    N-terminal deltas follow the first residue; the first N-term-capable delta
    there goes to the N-term token (``K+42.01`` is read as N-terminal acetyl,
    the common case). Returns None for an unmatched delta."""
    b = s
    if len(b) > 4 and b[1] == "." and b[-2] == ".":
        b = b[2:-2]
    seq, mods = [], []
    nterm_used = False
    for tok in re.finditer(r"([A-Z])|([+-]\d+(?:\.\d+)?)", b):
        if tok.group(1):
            seq.append(tok.group(1))
            continue
        d = float(tok.group(2))
        pos = len(seq)
        if pos == 0:
            return None
        aa = seq[-1]
        res_id = match_mass(d, aa, False)
        nt_id = match_mass(d, None, True) if pos == 1 and not nterm_used else None
        if nt_id is not None:
            mods.append((0, nt_id))
            nterm_used = True
        elif res_id is not None:
            mods.append((pos, res_id))
        else:
            return None
    return "".join(seq), sorted(mods)


def encode(seq: str, mods):
    """-> (tokens int64 [L+2], comp float32 [L+2, 8]) or None if unsupported.

    I is encoded as L: the training sequences (quantms reanalysis) carry no I."""
    seq = seq.replace("I", "L")
    if not (1 <= len(seq) <= MAX_LEN) or any(c not in AA_TOK for c in seq):
        return None
    n = len(seq)
    tok = np.zeros(n + 2, dtype=np.int64)
    tok[0], tok[-1] = 1, 2
    tok[1:-1] = [AA_TOK[c] for c in seq]
    comp = np.zeros((n + 2, N_COMP), dtype=np.float32)
    for pos, uid in mods:
        if uid not in VECTORS or not 0 <= pos <= n + 1:
            return None
        comp[pos] += np.asarray(VECTORS[uid], dtype=np.float32)
    comp *= COMP_SCALE
    return tok, comp


def batch(encoded):
    """List of (tok, comp) -> padded tensors (B, Lmax), (B, Lmax, 8)."""
    L = max(len(t) for t, _ in encoded)
    T = np.zeros((len(encoded), L), dtype=np.int64)
    C = np.zeros((len(encoded), L, N_COMP), dtype=np.float32)
    for i, (t, c) in enumerate(encoded):
        T[i, : len(t)] = t
        C[i, : len(t)] = c
    return torch.from_numpy(T), torch.from_numpy(C)


class RtCnn(nn.Module):
    def __init__(self, emb: int = 32, ch: int = 64, dilations=DILATIONS):
        super().__init__()
        self.dilations = tuple(dilations)
        self.emb = nn.Embedding(N_TOK, emb, padding_idx=0)
        self.inp = nn.Conv1d(emb + N_COMP, ch, 1)
        self.convs = nn.ModuleList(nn.Conv1d(ch, ch, 3, padding=d, dilation=d) for d in self.dilations)
        self.mix = nn.ModuleList(nn.Conv1d(ch, ch, 1) for _ in self.dilations)
        self.att = nn.Conv1d(ch, 1, 1)
        self.fc1 = nn.Linear(2 * ch, ch)
        self.fc2 = nn.Linear(ch, 1)

    def forward(self, tok, comp):
        m = (tok > 0).unsqueeze(1).to(comp.dtype)
        x = torch.cat([self.emb(tok), comp], -1).transpose(1, 2)
        x = self.inp(x) * m
        for conv, mix in zip(self.convs, self.mix):
            x = x + mix(F.silu(conv(x))) * m
        a = self.att(x).masked_fill(m == 0, float("-inf"))
        a = torch.softmax(a, -1)
        pooled = torch.cat([(a * x).sum(-1), x.sum(-1) / m.sum(-1)], -1)
        return self.fc2(F.silu(self.fc1(pooled))).squeeze(-1)


def export_weights(model: RtCnn, path: str):
    sd = {k: v.detach().cpu().float().numpy() for k, v in model.state_dict().items()}
    sd["dilations"] = np.asarray(model.dilations, dtype=np.float32)
    with open(path, "wb") as f:
        f.write(b"ARTN0001")
        f.write(struct.pack("<I", len(sd)))
        for name, arr in sd.items():
            nb = name.encode()
            f.write(struct.pack("<I", len(nb)))
            f.write(nb)
            f.write(struct.pack("<I", arr.ndim))
            f.write(struct.pack(f"<{arr.ndim}I", *arr.shape))
            f.write(np.ascontiguousarray(arr, dtype="<f4").tobytes())


def predict(model: RtCnn, forms, device="cpu", bs: int = 4096, fmt: str = "proforma"):
    """Index per peptidoform string (NaN when unsupported)."""
    parse = parse_proforma if fmt == "proforma" else parse_andes
    out = np.full(len(forms), np.nan, dtype=np.float64)
    enc, idx = [], []
    for i, s in enumerate(forms):
        p = parse(s)
        e = encode(*p) if p else None
        if e is not None:
            enc.append(e)
            idx.append(i)
    model.eval()
    order = np.argsort([len(e[0]) for e in enc], kind="stable")
    with torch.no_grad():
        for k in range(0, len(order), bs):
            sel = order[k: k + bs]
            T, C = batch([enc[j] for j in sel])
            v = model(T.to(device), C.to(device)).double().cpu().numpy()
            out[np.asarray(idx)[sel]] = v
    return out


if __name__ == "__main__":
    assert parse_proforma("[UNIMOD:737]-PEPS[UNIMOD:21]K[UNIMOD:737]") == ("PEPSK", [(0, 737), (4, 21), (5, 737)])
    assert parse_proforma("AC[UNIMOD:4]K") == ("ACK", [(2, 4)])
    assert parse_andes("K.M+42.01056+15.99491PEPS+79.96633K.D") == ("MPEPSK", [(0, 1), (1, 35), (5, 21)])
    assert parse_andes("K.K+229.16293+229.16293PEPK+229.16293.D") == ("KPEPK", [(0, 737), (1, 737), (5, 737)])
    assert parse_andes("K.Q-17.02655PEPK.D") == ("QPEPK", [(0, 28)])
    assert parse_andes("-.K+42.01056PEPK.-") == ("KPEPK", [(0, 1)])
    assert parse_andes("K.PEPXK.D") is None
    m = RtCnn()
    print("params", sum(p.numel() for p in m.parameters()))
    e = [encode("PEPTIDEK", []), encode("AMCK", [(2, 35), (3, 4)])]
    T, C = batch(e)
    full = m(T, C)
    solo = torch.stack([m(*batch([x]))[0] for x in e])
    assert torch.allclose(full, solo, atol=1e-5), (full, solo)
    print("ok", full.detach().numpy())
