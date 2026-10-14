#!/usr/bin/env python3
"""Train andes-rt v1 on per-(run, peptidoform) median RTs (``rt_extract.py``).

Loss, per mini-batch of R runs x K peptidoforms from the same run:
- pairwise logistic ranking inside each run (pairs closer than 1% of the run's
  RT span are ignored), so runs with different gradients need no alignment;
- plus MSE of a per-run affine head ``a_r * index + b_r`` against the run's
  RT scaled to [0, 1] by its 1st/99th percentiles (weight --affine-weight).

Rows of --holdout-datasets and the "test" split never train. Ablations:
--exclude-phospho drops every phospho peptidoform; --exclude-tmt drops every
TMT-labelled one. Few-shot: --init a checkpoint, --fewshot-ptm phospho|tmt and
--fewshot-n N add N train-split peptidoforms of that PTM (all their runs) and
fine-tune on half PTM / half general batches; --fewshot-ptm implies the
matching --exclude-* flag, so no other peptidoform of that PTM is trained on.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pyarrow.parquet as pq
import torch
import torch.nn.functional as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rt_model import RtCnn, encode, export_weights, parse_proforma  # noqa: E402


def load_table(path, holdout):
    """Rows as categoricals; ``run_id`` numbers (dataset, run) pairs, ``form_id`` peptidoforms."""
    cats = ["dataset", "run", "peptidoform", "split"]
    t = pq.read_table(path, columns=cats + ["rt_min", "has_phospho", "has_tmt"], read_dictionary=cats).to_pandas()
    for c in cats:
        t[c] = t[c].astype("category")
    t["run_id"] = t.groupby(["dataset", "run"], observed=True, sort=True).ngroup()
    t["form_id"] = t.peptidoform.cat.codes.astype(np.int64)
    t["holdout"] = t.dataset.isin(holdout)
    return t


def encode_forms(forms):
    toks, comps, ok = [], [], np.zeros(len(forms), bool)
    for i, f in enumerate(forms):
        e = encode(*parse_proforma(f))
        if e is None:
            toks.append(None)
            comps.append(None)
            continue
        toks.append(e[0])
        comps.append(e[1])
        ok[i] = True
    return toks, comps, ok


class Batcher:
    def __init__(self, rows, toks, comps, runs_per_batch, per_run, rng, min_rows=16):
        self.rng = rng
        self.toks, self.comps = toks, comps
        self.runs_per_batch, self.per_run = runs_per_batch, per_run
        g = rows.groupby("run_id").indices
        self.groups = [(r, idx) for r, idx in g.items() if len(idx) >= min_rows]
        if not self.groups:
            raise ValueError("no run has enough rows to sample")
        self.runs_per_batch = min(runs_per_batch, len(self.groups))
        w = np.array([np.sqrt(len(i)) for _, i in self.groups])
        self.p = w / w.sum()
        self.form = rows.form_id.values
        self.y = rows.rt_norm.values.astype(np.float32)

    def sample(self):
        sel = self.rng.choice(len(self.groups), self.runs_per_batch, replace=False, p=self.p)
        rid, items = [], []
        for s in sel:
            r, idx = self.groups[s]
            k = min(self.per_run, len(idx))
            items.append(self.rng.choice(idx, k, replace=False))
            rid.append(np.full(k, r))
        items = np.concatenate(items)
        rid = np.concatenate(rid)
        f = self.form[items]
        L = max(len(self.toks[j]) for j in f)
        T = np.zeros((len(f), L), np.int64)
        C = np.zeros((len(f), L, 8), np.float32)
        for i, j in enumerate(f):
            n = len(self.toks[j])
            T[i, :n] = self.toks[j]
            C[i, :n] = self.comps[j]
        return torch.from_numpy(T), torch.from_numpy(C), torch.from_numpy(rid), torch.from_numpy(self.y[items])


def rank_loss(s, y, rid):
    same = rid[:, None] == rid[None, :]
    dy = y[:, None] - y[None, :]
    valid = same & (dy.abs() > 0.01)
    ds = s[:, None] - s[None, :]
    return F.softplus(-torch.sign(dy) * ds)[valid].mean()


def select_rows(t, a, rng):
    """(train, val, few) row sets of the RT table for the parsed arguments ``a``.

    ``few`` is None unless ``a.fewshot_ptm`` is set. A few-shot PTM implies its
    exclusion flag: every peptidoform carrying it is dropped from ``train`` and
    ``val``, so the ``a.fewshot_n`` forms in ``few`` are the only rows of that
    PTM the run sees (the checkpoint arguments are not restored from --init)."""
    base = (~t.holdout) & (t.split != "test")
    if a.exclude_phospho or a.fewshot_ptm == "phospho":
        base &= ~t.has_phospho
    if a.exclude_tmt or a.fewshot_ptm == "tmt":
        base &= ~t.has_tmt
    train = t[base & (t.split == "train")]
    val = t[base & (t.split == "val")]
    few = None
    if a.fewshot_ptm:
        col = "has_phospho" if a.fewshot_ptm == "phospho" else "has_tmt"
        pool = t[(~t.holdout) & (t.split == "train") & t[col]]
        pick = rng.choice(np.asarray(pool.peptidoform.unique()), a.fewshot_n, replace=False)
        picked = pool.peptidoform.isin(set(pick))
        runs = set(pool[picked].run_id)
        # the picked PTM peptidoforms plus every non-PTM train row of the same runs
        cand = t[t.run_id.isin(runs) & (~t.holdout) & (t.split == "train")]
        few = cand[(~cand[col]) | cand.peptidoform.isin(set(pick))]
        few = few[few.run_id.isin(set(few[few[col]].run_id))]
        print(f"few-shot rows: {int(few[col].sum())} PTM rows of {few[col].sum() and few[few[col]].peptidoform.nunique()} "
              f"peptidoforms in {few.run_id.nunique()} runs", flush=True)
        print(f"few-shot {a.fewshot_ptm}: {a.fewshot_n} peptidoforms, {len(few)} rows incl. partners", flush=True)
    return train, val, few


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--table", required=True)
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--holdout-datasets", nargs="*", default=["PXD019263", "PXD005235", "PXD012255"])
    p.add_argument("--exclude-phospho", action="store_true")
    p.add_argument("--exclude-tmt", action="store_true")
    p.add_argument("--init")
    p.add_argument("--fewshot-ptm", choices=["phospho", "tmt"])
    p.add_argument("--fewshot-n", type=int, default=0)
    p.add_argument("--steps", type=int, default=40000)
    p.add_argument("--lr", type=float, default=2e-3)
    p.add_argument("--runs-per-batch", type=int, default=16)
    p.add_argument("--per-run", type=int, default=64)
    p.add_argument("--affine-weight", type=float, default=1.0)
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    torch.manual_seed(a.seed)
    rng = np.random.default_rng(a.seed)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    if dev == "cpu":
        torch.set_num_threads(int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
    t0 = time.time()
    t = load_table(a.table, set(a.holdout_datasets))
    forms = t.peptidoform.cat.categories
    toks, comps, ok = encode_forms(list(forms))
    t = t[ok[t.form_id.values]].copy()
    lo = t.groupby("run_id").rt_min.transform(lambda v: v.quantile(0.01))
    hi = t.groupby("run_id").rt_min.transform(lambda v: v.quantile(0.99))
    t["rt_norm"] = ((t.rt_min - lo) / (hi - lo).clip(lower=1e-3)).clip(-0.2, 1.2)
    n_runs = int(t.run_id.max()) + 1
    print(f"table {len(t)} rows, {len(forms)} forms ({(~ok).sum()} unsupported), {n_runs} runs, "
          f"device {dev}, load {time.time() - t0:.0f}s", flush=True)

    train, val, few = select_rows(t, a, rng)
    print(f"train rows {len(train)} val rows {len(val)}", flush=True)

    model = RtCnn().to(dev)
    affine = torch.nn.Embedding(n_runs, 2).to(dev)
    with torch.no_grad():
        affine.weight[:, 0] = 0.1
        affine.weight[:, 1] = 0.5
    if a.init:
        model.load_state_dict(torch.load(os.path.join(a.init, "model.pt"), map_location=dev))
        affine.load_state_dict(torch.load(os.path.join(a.init, "affine.pt"), map_location=dev))
    opt = torch.optim.AdamW(list(model.parameters()) + list(affine.parameters()), lr=a.lr, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=a.lr, total_steps=a.steps, pct_start=0.05)
    bt = Batcher(train, toks, comps, a.runs_per_batch, a.per_run, rng)
    bf = Batcher(few, toks, comps, a.runs_per_batch, a.per_run, rng, min_rows=4) if few is not None else None
    bv = Batcher(val, toks, comps, a.runs_per_batch, a.per_run, np.random.default_rng(123))
    val_batches = [bv.sample() for _ in range(50)]

    def evaluate():
        model.eval()
        rs = []
        with torch.no_grad():
            for T, C, rid, y in val_batches:
                s = model(T.to(dev), C.to(dev)).cpu()
                for r in rid.unique():
                    m = rid == r
                    if m.sum() > 8:
                        rs.append(np.corrcoef(s[m].numpy(), y[m].numpy())[0, 1])
        model.train()
        return float(np.nanmean(rs))

    best, hist = -1.0, []
    for step in range(1, a.steps + 1):
        if bf is not None and step % 2 == 0:
            T, C, rid, y = bf.sample()
        else:
            T, C, rid, y = bt.sample()
        T, C, rid, y = T.to(dev), C.to(dev), rid.to(dev), y.to(dev)
        s = model(T, C)
        ab = affine(rid)
        loss = rank_loss(s, y, rid) + a.affine_weight * F.mse_loss(ab[:, 0] * s + ab[:, 1], y)
        opt.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step()
        sched.step()
        if step % 1000 == 0 or step == a.steps:
            r = evaluate()
            hist.append((step, float(loss), r))
            print(f"step {step} loss {loss.item():.4f} val per-run pearson {r:.4f} ({time.time() - t0:.0f}s)", flush=True)
            if r > best:
                best = r
                torch.save(model.state_dict(), os.path.join(a.out, "model.pt"))
                torch.save(affine.state_dict(), os.path.join(a.out, "affine.pt"))
                export_weights(model, os.path.join(a.out, "andes_rt.artn"))
                json.dump({"args": vars(a), "best_val_pearson": best, "best_step": step, "history": hist},
                          open(os.path.join(a.out, "train.json"), "w"), indent=1)
    model.load_state_dict(torch.load(os.path.join(a.out, "model.pt"), map_location=dev))
    export_weights(model.cpu(), os.path.join(a.out, "andes_rt.artn"))
    json.dump({"args": vars(a), "best_val_pearson": best, "history": hist,
               "params": sum(p.numel() for p in model.parameters())}, open(os.path.join(a.out, "train.json"), "w"), indent=1)
    print(f"done best val per-run pearson {best:.4f}; exported {a.out}/andes_rt.artn", flush=True)


if __name__ == "__main__":
    main()
