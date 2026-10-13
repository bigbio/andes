"""andes-ms2 v1: fragment-intensity CNN, its loss, and the weight export format.

Peptide encoding is andes-rt's (``rt_model.encode``): tokens ``[N-term, residues,
C-term]`` (I as L) plus a per-position modification composition vector.

Conditioning vector (``cond``, length ``MAX_CHARGE + 2 + 2 * cond_emb``):
precursor charge one-hot (charges above ``MAX_CHARGE`` use the last slot),
NCE / 100 and an NCE-known flag (0 / 0 when unknown), an activation embedding
(``ACTIVATIONS``) and an instrument-family embedding (``INSTRUMENTS``).

Network::

    x = conv1x1(emb(tok) ++ comp)                         (B, ch, L+2), padding zeroed
    film = lin2(silu(lin1(cond)))  -> per block (gamma, beta), each ch wide
    block b: x += mix_b(silu(conv_b(x) * (1 + gamma_b) + beta_b))      k3, dilation d_b
    optional: x += proj(multi-head self-attention(x))   (padding keys masked)
    site j (cleavage after residue j+1): h = silu(conv_k2(x)[j+1] + cond_lin(cond))
    out = conv1x1(h) -> N_ION values per site (linear; clipped at 0 when used)

Every layer output is multiplied by the padding mask, so a peptide's prediction
does not depend on what it was batched with.

Export (``AMS20001``): the ``ARTN0001`` tensor layout (magic, u32 count, then per
tensor u32 name length, name, u32 ndim, u32 dims, little-endian f32 data) with
extra tensors ``dilations`` and ``n_heads``.
"""
from __future__ import annotations

import math
import os
import struct
import sys

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
sys.path.insert(0, _HERE)
from ms2_common import (  # noqa: E402
    ACTIVATIONS, INSTRUMENTS, MAX_CHARGE, N_ION, activation_index, instrument_class,
    parse_proforma, valid_mask,
)
from rt_model import N_TOK, encode  # noqa: E402
from mod_composition import N_COMP  # noqa: E402

MAGIC = b"AMS20001"


def cond_numeric(charge: int, nce: float) -> np.ndarray:
    """Charge one-hot ++ [NCE/100, NCE known] (the non-embedded part of ``cond``)."""
    v = np.zeros(MAX_CHARGE + 2, np.float32)
    v[min(max(int(charge), 1), MAX_CHARGE) - 1] = 1.0
    if nce is not None and np.isfinite(nce):
        v[MAX_CHARGE] = float(nce) / 100.0
        v[MAX_CHARGE + 1] = 1.0
    return v


class Ms2Cnn(nn.Module):
    def __init__(self, emb: int = 32, ch: int = 96, dilations=(1, 2, 4, 8, 1, 2), n_heads: int = 4,
                 cond_emb: int = 8, film_hidden: int = 64):
        super().__init__()
        self.dilations = tuple(dilations)
        self.n_heads = n_heads
        self.ch = ch
        self.emb = nn.Embedding(N_TOK, emb, padding_idx=0)
        self.act_emb = nn.Embedding(len(ACTIVATIONS), cond_emb)
        self.instr_emb = nn.Embedding(len(INSTRUMENTS), cond_emb)
        cond_dim = MAX_CHARGE + 2 + 2 * cond_emb
        self.film1 = nn.Linear(cond_dim, film_hidden)
        self.film2 = nn.Linear(film_hidden, 2 * ch * len(self.dilations))
        nn.init.zeros_(self.film2.weight)
        nn.init.zeros_(self.film2.bias)
        self.inp = nn.Conv1d(emb + N_COMP, ch, 1)
        self.convs = nn.ModuleList(nn.Conv1d(ch, ch, 3, padding=d, dilation=d) for d in self.dilations)
        self.mix = nn.ModuleList(nn.Conv1d(ch, ch, 1) for _ in self.dilations)
        if n_heads > 0:
            assert ch % n_heads == 0, "ch must be divisible by n_heads"
            self.qkv = nn.Conv1d(ch, 3 * ch, 1)
            self.proj = nn.Conv1d(ch, ch, 1)
        self.head1 = nn.Conv1d(ch, ch, 2)
        self.head_cond = nn.Linear(cond_dim, ch)
        self.head2 = nn.Conv1d(ch, N_ION, 1)

    def cond(self, num, act, instr):
        return torch.cat([num, self.act_emb(act), self.instr_emb(instr)], -1)

    def forward(self, tok, comp, num, act, instr):
        """-> (B, L, N_ION) for L = padded residue count; site j valid for j < n - 1."""
        B, Lp = tok.shape
        m = (tok > 0).unsqueeze(1).to(comp.dtype)
        c = self.cond(num, act, instr)
        film = self.film2(F.silu(self.film1(c))).view(B, len(self.dilations), 2, self.ch, 1)
        x = torch.cat([self.emb(tok), comp], -1).transpose(1, 2)
        x = self.inp(x) * m
        for b, (conv, mix) in enumerate(zip(self.convs, self.mix)):
            h = conv(x) * (1 + film[:, b, 0]) + film[:, b, 1]
            x = x + mix(F.silu(h)) * m
        if self.n_heads > 0:
            H, dh = self.n_heads, self.ch // self.n_heads
            q, k, v = self.qkv(x).chunk(3, dim=1)
            q, k, v = (t.reshape(B, H, dh, Lp) for t in (q, k, v))
            s = torch.einsum("bhdi,bhdj->bhij", q, k) / math.sqrt(dh)
            s = s.masked_fill(m.unsqueeze(1) == 0, float("-inf"))
            o = torch.einsum("bhij,bhdj->bhdi", torch.softmax(s, -1), v).reshape(B, self.ch, Lp)
            x = x + self.proj(o) * m
        h = F.silu(self.head1(x[:, :, 1:]) + self.head_cond(c).unsqueeze(-1))
        return self.head2(h).transpose(1, 2)


def masked_sa_loss(pred, target, mask):
    """Mean (1 - spectral angle) over the batch; ``mask`` selects valid channels."""
    p = (pred * mask).flatten(1)
    t = (target.clamp(min=0) * mask).flatten(1)
    cos = (p * t).sum(1) / (p.norm(dim=1) * t.norm(dim=1) + 1e-8)
    sa = 1 - 2 * torch.acos(cos.clamp(-1 + 1e-7, 1 - 1e-7)) / math.pi
    return (1 - sa).mean()


def teacher_kl(pred, teacher, mask):
    """KL(teacher || student) between sum-normalised intensity distributions over valid channels."""
    m = mask * (teacher >= 0)
    t = teacher.clamp(min=0) * m
    t = t.flatten(1) / (t.flatten(1).sum(1, keepdim=True) + 1e-8)
    s = (F.relu(pred) * m).flatten(1) + 1e-6 * m.flatten(1)
    s = s / (s.sum(1, keepdim=True) + 1e-8)
    return (t * (torch.log(t + 1e-8) - torch.log(s + 1e-8))).sum(1).mean()


def spectral_angle(pred: np.ndarray, obs: np.ndarray, valid: np.ndarray) -> float:
    p = np.where(valid, np.maximum(pred, 0), 0).ravel()
    o = np.where(valid, np.maximum(obs, 0), 0).ravel()
    np_, no = np.linalg.norm(p), np.linalg.norm(o)
    if np_ == 0 or no == 0:
        return 0.0
    cos = float(np.clip(p @ o / (np_ * no), -1, 1))
    return 1 - 2 * math.acos(cos) / math.pi


def export_weights(model: Ms2Cnn, path: str):
    sd = {k: v.detach().cpu().float().numpy() for k, v in model.state_dict().items()}
    sd["dilations"] = np.asarray(model.dilations, dtype=np.float32)
    sd["n_heads"] = np.asarray([model.n_heads], dtype=np.float32)
    with open(path, "wb") as f:
        f.write(MAGIC)
        f.write(struct.pack("<I", len(sd)))
        for name, arr in sd.items():
            nb = name.encode()
            f.write(struct.pack("<I", len(nb)))
            f.write(nb)
            f.write(struct.pack("<I", arr.ndim))
            f.write(struct.pack(f"<{arr.ndim}I", *arr.shape))
            f.write(np.ascontiguousarray(arr, dtype="<f4").tobytes())


def config_of(model: Ms2Cnn) -> dict:
    return {"emb": model.emb.embedding_dim, "ch": model.ch, "dilations": list(model.dilations),
            "n_heads": model.n_heads, "cond_emb": model.act_emb.embedding_dim,
            "film_hidden": model.film1.out_features}


def make_batch(items):
    """items: list of (tok, comp, charge, nce, activation idx, instrument idx) -> tensors."""
    L = max(len(t) for t, *_ in items)
    T = np.zeros((len(items), L), np.int64)
    C = np.zeros((len(items), L, N_COMP), np.float32)
    num = np.zeros((len(items), MAX_CHARGE + 2), np.float32)
    act = np.zeros(len(items), np.int64)
    ins = np.zeros(len(items), np.int64)
    for i, (t, c, z, nce, a, k) in enumerate(items):
        T[i, : len(t)] = t
        C[i, : len(t)] = c
        num[i] = cond_numeric(z, nce)
        act[i], ins[i] = a, k
    return (torch.from_numpy(T), torch.from_numpy(C), torch.from_numpy(num),
            torch.from_numpy(act), torch.from_numpy(ins))


def predict(model: Ms2Cnn, requests, bs: int = 1024, device="cpu", raw: bool = False):
    """requests: (proforma, charge, nce, activation name, instrument name) ->
    list of (n-1, N_ION) arrays clipped at 0 with invalid channels set to 0 (the
    unprocessed network output with ``raw``), or None when the peptidoform cannot
    be encoded."""
    out = [None] * len(requests)
    enc = []
    for i, (pf, z, nce, act, instr) in enumerate(requests):
        p = parse_proforma(pf)
        e = encode(*p) if p else None
        if e is None or len(p[0]) < 2 or z < 1:
            continue
        enc.append((i, p, (e[0], e[1], z, nce, activation_index(act), instrument_class(instr))))
    model.eval()
    with torch.no_grad():
        for k in range(0, len(enc), bs):
            chunk = enc[k: k + bs]
            y = model(*(t.to(device) for t in make_batch([c[2] for c in chunk]))).cpu().numpy()
            for b, (i, (seq, mods), it) in enumerate(chunk):
                n = len(seq)
                v = valid_mask(seq, mods, it[2])
                y_b = y[b, : n - 1]
                out[i] = (y_b if raw else np.where(v, np.maximum(y_b, 0), 0)).astype(np.float32)
    return out


if __name__ == "__main__":
    torch.manual_seed(0)
    m = Ms2Cnn()
    print("params", sum(p.numel() for p in m.parameters()))
    reqs = [("PEPTIDEK", 2, 30.0, "HCD-FTMS", "Q Exactive"),
            ("[UNIMOD:737]-AMS[UNIMOD:21]C[UNIMOD:4]K[UNIMOD:737]", 3, float("nan"), "CID-ITMS", "LTQ Orbitrap Velos"),
            ("LLLLLLLLLLLLLLLLLLLLLLLLLLLLLK", 4, 27.0, "HCD-ASTMS", "Orbitrap Astral")]
    with torch.no_grad():
        nn.init.normal_(m.film2.weight, std=0.1)
    full = predict(m, reqs)
    solo = [predict(m, [r])[0] for r in reqs]
    for a_, b_ in zip(full, solo):
        assert np.allclose(a_, b_, atol=1e-5)
    assert full[1].shape == (4, N_ION) and (full[0][:, 4:] == 0).all()
    print("ok", [f.shape for f in full])
