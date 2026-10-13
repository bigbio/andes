#!/usr/bin/env python3
"""Train andes-ms2 on annotated spectra (``ms2_extract_raw.py`` / ``ms2_extract_flats.py`` shards).

Splits: train / val / test by an md5 of the stripped sequence (``ms2_common.split_of``,
the andes-rt rule, so a sequence is in the same split for both models); rows of
--holdout-datasets never train or validate and are reported separately.

Loss per batch: mean (1 - spectral angle) over all valid channels, plus
--teacher-weight x KL(teacher || student) when the shards carry a ``teacher``
column (``ms2_teacher.py``). Validation metric: median spectral angle of the
plain b/y 1+/2+ channels (the ``rescore_pin`` feature) on a fixed val subset.
The best checkpoint is exported as ``andes_ms2.ams2`` whenever it improves.

At the end the best model is scored on the test split and the held-out datasets,
stratified by activation, precursor charge, length, dataset and PTM; everything
goes to ``train.json``.

Usage::

  python ms2_train.py --shards SHARDS/ --out RUN/ [--steps 200000 --batch 512]
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ms2_common import N_BASE_ION, N_ION, activation_index, instrument_class, parse_proforma, split_of  # noqa: E402
from ms2_model import Ms2Cnn, config_of, export_weights, make_batch, masked_sa_loss, teacher_kl  # noqa: E402
from rt_model import encode  # noqa: E402

COLS = ["dataset", "peptidoform", "sequence", "charge", "nce", "activation", "instrument",
        "n_matched", "explained", "intensities"]


def shard_files(paths):
    out = []
    for p in paths:
        out += sorted(glob.glob(os.path.join(p, "**", "*.parquet"), recursive=True)) if os.path.isdir(p) else [p]
    return out


class Data:
    """All rows in memory: per-row metadata arrays plus flat float16 intensities."""

    def __init__(self, files, min_matched, min_explained, max_rows, seed, use_teacher):
        tables = []
        has_teacher = use_teacher
        for f in files:
            names = pq.ParquetFile(f).schema_arrow.names
            has_teacher = has_teacher and "teacher" in names
        for f in files:
            t = pq.read_table(f, columns=COLS + (["teacher"] if has_teacher else []))
            t = t.filter(pc.and_(pc.greater_equal(t["n_matched"], min_matched),
                                 pc.greater_equal(t["explained"], min_explained)))
            tables.append(t)
        t = pa.concat_tables(tables)
        if max_rows and t.num_rows > max_rows:
            idx = np.sort(np.random.default_rng(seed).choice(t.num_rows, max_rows, replace=False))
            t = t.take(pa.array(idx))
        self.n = t.num_rows
        self.has_teacher = has_teacher
        self.dataset = np.asarray(t["dataset"].to_pylist(), dtype=object)
        forms = pc.dictionary_encode(t["peptidoform"]).combine_chunks()
        self.form_id = forms.indices.to_numpy().astype(np.int64)
        self.forms = forms.dictionary.to_pylist()
        self.enc = []
        self.seqmods = []
        for f in self.forms:
            p = parse_proforma(f)
            self.seqmods.append(p)
            self.enc.append(encode(*p) if p else None)
        self.charge = t["charge"].to_numpy().astype(np.int64)
        self.nce = t["nce"].to_numpy(zero_copy_only=False).astype(np.float32)
        acts = pc.dictionary_encode(t["activation"]).combine_chunks()
        self.act = np.array([activation_index(x) for x in acts.dictionary.to_pylist()])[acts.indices.to_numpy()]
        ins = pc.dictionary_encode(t["instrument"]).combine_chunks()
        self.instr = np.array([instrument_class(x) for x in ins.dictionary.to_pylist()])[ins.indices.to_numpy()]
        self.act_name = np.asarray(acts.dictionary.to_pylist(), dtype=object)[acts.indices.to_numpy()]

        def flat(col):
            arr = t[col].combine_chunks()
            off = arr.offsets.to_numpy().astype(np.int64)
            return off, arr.values.to_numpy(zero_copy_only=False).astype(np.float16)

        self.off, self.vals = flat("intensities")
        self.teacher = flat("teacher")[1] if has_teacher else None
        seqs = pc.dictionary_encode(t["sequence"]).combine_chunks()
        sp = np.array([split_of(s) for s in seqs.dictionary.to_pylist()], dtype=object)
        self.split = sp[seqs.indices.to_numpy()]
        self.length = np.array([len(p[0]) if p else 0 for p in self.seqmods])[self.form_id]
        ok = np.array([e is not None for e in self.enc])[self.form_id]
        ok &= (self.off[1:] - self.off[:-1]) == (self.length - 1) * N_ION
        self.ok = ok

    def batch(self, idx, device):
        items = []
        for i in idx:
            e = self.enc[self.form_id[i]]
            items.append((e[0], e[1], int(self.charge[i]), float(self.nce[i]), int(self.act[i]), int(self.instr[i])))
        T, C, num, act, ins = (x.to(device) for x in make_batch(items))
        S = T.shape[1] - 2  # sites in the padded batch (model output rows)
        tgt = np.full((len(idx), S, N_ION), -1.0, np.float32)
        tch = np.full((len(idx), S, N_ION), -1.0, np.float32) if self.teacher is not None else None
        for b, i in enumerate(idx):
            n = self.length[i] - 1
            tgt[b, :n] = self.vals[self.off[i]: self.off[i + 1]].reshape(n, N_ION)
            if tch is not None:
                tch[b, :n] = self.teacher[self.off[i]: self.off[i + 1]].reshape(n, N_ION)
        tgt = torch.from_numpy(tgt).to(device)
        tch = torch.from_numpy(tch).to(device) if tch is not None else None
        return (T, C, num, act, ins), tgt, tch


def sa_rows(model, data, idx, device, bs=1024):
    """Per-row spectral angle on (all valid channels, plain b/y 1+/2+)."""
    model.eval()
    out_all, out_base = np.zeros(len(idx)), np.zeros(len(idx))
    with torch.no_grad():
        for k in range(0, len(idx), bs):
            sel = idx[k: k + bs]
            x, tgt, _ = data.batch(sel, device)
            pred = model(*x).clamp(min=0)
            mask = (tgt >= 0).float()
            t = tgt.clamp(min=0) * mask
            p = pred * mask
            for name, out, sl in (("all", out_all, slice(None)), ("base", out_base, slice(0, N_BASE_ION))):
                pp, tt = p[:, :, sl].flatten(1), t[:, :, sl].flatten(1)
                cos = (pp * tt).sum(1) / (pp.norm(dim=1) * tt.norm(dim=1) + 1e-12)
                sa = 1 - 2 * torch.acos(cos.clamp(-1, 1)) / np.pi
                sa[(pp.norm(dim=1) == 0) | (tt.norm(dim=1) == 0)] = 0
                out[k: k + len(sel)] = sa.cpu().numpy()
    model.train()
    return out_all, out_base


def strata(data, idx, sa):
    def med(mask):
        return {"n": int(mask.sum()), "median_sa": round(float(np.median(sa[mask])), 4) if mask.any() else None}

    res = {"all": med(np.ones(len(idx), bool))}
    for name, key in (("activation", data.act_name[idx]), ("charge", data.charge[idx]),
                      ("dataset", data.dataset[idx])):
        res[name] = {str(k): med(key == k) for k in np.unique(key)}
    L = data.length[idx]
    res["length"] = {f"{lo}-{hi}": med((L >= lo) & (L <= hi)) for lo, hi in ((6, 10), (11, 15), (16, 20), (21, 50))}
    forms = np.asarray(data.forms, dtype=object)[data.form_id[idx]]
    res["ptm"] = {
        "phospho": med(np.array(["UNIMOD:21]" in f for f in forms])),
        "oxidation": med(np.array(["UNIMOD:35]" in f for f in forms])),
        "tmt": med(np.array([("UNIMOD:737]" in f or "UNIMOD:2016]" in f) for f in forms])),
    }
    return res


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--shards", nargs="+", required=True, help="parquet files or directories")
    p.add_argument("--out", required=True)
    p.add_argument("--holdout-datasets", nargs="*", default=["PXD019263", "PXD005235", "PXD012255"])
    p.add_argument("--min-matched", type=int, default=4)
    p.add_argument("--min-explained", type=float, default=0.0)
    p.add_argument("--max-rows", type=int, default=0, help="subsample the shards (0 = all)")
    p.add_argument("--emb", type=int, default=32)
    p.add_argument("--ch", type=int, default=96)
    p.add_argument("--dilations", type=int, nargs="+", default=[1, 2, 4, 8, 1, 2])
    p.add_argument("--n-heads", type=int, default=4, help="0 disables the attention layer")
    p.add_argument("--cond-emb", type=int, default=8)
    p.add_argument("--film-hidden", type=int, default=64)
    p.add_argument("--steps", type=int, default=200000)
    p.add_argument("--batch", type=int, default=512)
    p.add_argument("--lr", type=float, default=2e-3)
    p.add_argument("--teacher-weight", type=float, default=0.0)
    p.add_argument("--balance-activation", action="store_true",
                   help="sample rows with weight 1/sqrt(rows of their activation)")
    p.add_argument("--val-n", type=int, default=20000)
    p.add_argument("--eval-every", type=int, default=2000)
    p.add_argument("--init", help="checkpoint directory to start from")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    torch.manual_seed(a.seed)
    rng = np.random.default_rng(a.seed)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    if dev == "cpu":
        torch.set_num_threads(a.threads)
    t0 = time.time()
    data = Data(shard_files(a.shards), a.min_matched, a.min_explained, a.max_rows, a.seed, a.teacher_weight > 0)
    hold = np.isin(data.dataset, a.holdout_datasets)
    train = np.flatnonzero(data.ok & ~hold & (data.split == "train"))
    val = np.flatnonzero(data.ok & ~hold & (data.split == "val"))
    test = np.flatnonzero(data.ok & ~hold & (data.split == "test"))
    held = np.flatnonzero(data.ok & hold)
    if len(train) == 0 or len(val) == 0:
        sys.exit("no train or val rows")
    val_idx = rng.choice(val, min(a.val_n, len(val)), replace=False)
    print(f"rows {data.n} ({int((~data.ok).sum())} unusable), forms {len(data.forms)}, train {len(train)}, "
          f"val {len(val)}, test {len(test)}, held-out {len(held)}, teacher {data.has_teacher}, device {dev}, "
          f"load {time.time() - t0:.0f}s", flush=True)
    for k in np.unique(data.act_name[train]):
        print(f"  train activation {k}: {int((data.act_name[train] == k).sum())}", flush=True)

    weights = None
    if a.balance_activation:
        counts = {k: (data.act[train] == k).sum() for k in np.unique(data.act[train])}
        w = np.array([1 / np.sqrt(counts[k]) for k in data.act[train]])
        weights = w / w.sum()

    model = Ms2Cnn(a.emb, a.ch, a.dilations, a.n_heads, a.cond_emb, a.film_hidden).to(dev)
    if a.init:
        model.load_state_dict(torch.load(os.path.join(a.init, "model.pt"), map_location=dev))
    n_params = sum(p.numel() for p in model.parameters())
    print(f"model params {n_params} config {config_of(model)}", flush=True)
    opt = torch.optim.AdamW(model.parameters(), lr=a.lr, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=a.lr, total_steps=a.steps, pct_start=0.05)
    best, hist, run_loss = -1.0, [], []
    meta = {"args": vars(a), "params": n_params, "config": config_of(model)}

    def save(step, val_sa):
        torch.save(model.state_dict(), os.path.join(a.out, "model.pt"))
        export_weights(model, os.path.join(a.out, "andes_ms2.ams2"))
        json.dump({**meta, "best_val_median_sa": val_sa, "best_step": step, "history": hist},
                  open(os.path.join(a.out, "train.json"), "w"), indent=1)

    _, sa0 = sa_rows(model, data, val_idx, dev)
    hist.append({"step": 0, "train_loss": None, "val_median_sa_base": round(float(np.median(sa0)), 5),
                 "seconds": round(time.time() - t0, 1)})
    print(f"step 0 val median SA (b/y 1-2+) {np.median(sa0):.4f}", flush=True)
    for step in range(1, a.steps + 1):
        idx = rng.choice(train, a.batch, replace=True, p=weights)
        x, tgt, tch = data.batch(idx, dev)
        pred = model(*x)
        mask = (tgt >= 0).float()
        loss = masked_sa_loss(pred, tgt, mask)
        if tch is not None:
            loss = loss + a.teacher_weight * teacher_kl(pred, tch, mask)
        opt.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step()
        sched.step()
        run_loss.append(float(loss))
        if step % a.eval_every == 0 or step == a.steps:
            _, sa_b = sa_rows(model, data, val_idx, dev)
            v = float(np.median(sa_b))
            tl = float(np.mean(run_loss))
            run_loss = []
            hist.append({"step": step, "train_loss": round(tl, 5), "val_median_sa_base": round(v, 5),
                         "seconds": round(time.time() - t0, 1)})
            print(f"step {step} train loss {tl:.4f} val median SA (b/y 1-2+) {v:.4f} ({time.time() - t0:.0f}s)",
                  flush=True)
            if v > best:
                best = v
                save(step, v)
    model.load_state_dict(torch.load(os.path.join(a.out, "model.pt"), map_location=dev))
    report = {}
    for name, idx in (("test", test), ("held_out_datasets", held)):
        if len(idx):
            sa_all, sa_b = sa_rows(model, data, idx, dev)
            report[name] = {"base": strata(data, idx, sa_b), "all_channels_median_sa": round(float(np.median(sa_all)), 4)}
            print(f"{name}: n {len(idx)} median SA b/y 1-2+ {np.median(sa_b):.4f}, all channels {np.median(sa_all):.4f}",
                  flush=True)
    model.cpu()
    export_weights(model, os.path.join(a.out, "andes_ms2.ams2"))
    json.dump({**meta, "best_val_median_sa": best, "history": hist, "report": report},
              open(os.path.join(a.out, "train.json"), "w"), indent=1)
    print(f"done; best val median SA {best:.4f}; exported {a.out}/andes_ms2.ams2", flush=True)


if __name__ == "__main__":
    main()
