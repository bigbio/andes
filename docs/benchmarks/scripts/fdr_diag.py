#!/usr/bin/env python3
# FDR sanity diagnostic on andes Astral Percolator output.
# Columns: PSMId  score  q-value  posterior_error_prob  peptide  proteinIds
import sys

def load(path):
    rows = []
    with open(path) as f:
        next(f)  # header
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            try:
                rows.append((float(p[1]), float(p[2]), float(p[3])))  # score, q, pep
            except ValueError:
                continue
    return rows

def pct(xs, q):
    if not xs:
        return float("nan")
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(q * len(xs)))]

for tag in ["astral_top1", "astral_chim"]:
    base = f"/srv/data/msgf-bench/bench-astral-mzml/{tag}"
    T = load(base + ".t.psms")
    D = load(base + ".d.psms")
    tq = [s for (s, q, p) in T if q <= 0.01]
    dq = [s for (s, q, p) in D if q <= 0.01]
    # score threshold at the 1% q boundary = lowest target score still q<=0.01
    thr = min((s for (s, q, p) in T if q <= 0.01), default=float("inf"))
    # PEP at the boundary = max PEP among accepted targets
    pep_at_bound = max((p for (s, q, p) in T if q <= 0.01), default=float("nan"))
    # decoys scoring above the 1% threshold (independent of q labelling)
    dec_above = sum(1 for (s, q, p) in D if s >= thr)
    tgt_above = sum(1 for (s, q, p) in T if s >= thr)
    # separation
    maxdec = max((s for (s, q, p) in D), default=float("nan"))
    tgt_gt_maxdec = sum(1 for (s, q, p) in T if s > maxdec)
    print(f"### {tag} ###")
    print(f"  total: targets={len(T)} decoys={len(D)}")
    print(f"  accepted @q<=0.01: targets={len(tq)}  decoys={len(dq)}  "
          f"(decoy/target = {100*len(dq)/max(len(tq),1):.2f}%)")
    print(f"  score@1%-boundary = {thr:.4f}   PEP@boundary = {pep_at_bound:.4g}")
    print(f"  by score>=boundary: targets={tgt_above} decoys={dec_above}  "
          f"=> empirical decoy FDR = {100*dec_above/max(tgt_above,1):.2f}%")
    print(f"  score sep: target p50={pct([s for s,q,p in T],0.5):.3f} "
          f"p90={pct([s for s,q,p in T],0.9):.3f} | "
          f"decoy p50={pct([s for s,q,p in D],0.5):.3f} "
          f"p90={pct([s for s,q,p in D],0.9):.3f} max={maxdec:.3f}")
    print(f"  targets scoring above ALL decoys (unambiguous): {tgt_gt_maxdec}")
    print()
