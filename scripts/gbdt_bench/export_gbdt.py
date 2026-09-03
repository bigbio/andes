#!/usr/bin/env python3
"""Export one shipped andes GBDT blob for the compiled-ensemble benchmark.

Writes, into --out:
  frag.agbd          the raw blob (the Rust micro-benchmark loads it with GbdtPeakModel::from_bytes)
  frag_lgbm.txt      the same trees in LightGBM model-text form (lleaves / treelite load this)
  rows.f32           N x n_features little-endian f32 feature rows, sampled per feature
                     UNIFORMLY BETWEEN THAT FEATURE'S MIN AND MAX SPLIT THRESHOLD, so every
                     internal node branches both ways at a realistic rate (rows drawn from
                     [0,1) would mostly fall one side of every split and flatter any evaluator)
  meta.json          n_trees, n_nodes, n_features, N, leaf statistics

Exactness is checked downstream by comparing each evaluator's predictions on rows.f32
against the Rust reference to 1e-6 (f32 accumulation order differs between evaluators,
so bit-identity is not expected across implementations - only across OUR implementations).
"""
import argparse, glob, io, json, os, struct
import numpy as np
import pyarrow.parquet as pq

def parse(blob):
    c = io.BytesIO(bytes(blob))
    magic = c.read(4); assert magic == b"AGBD", magic
    ver, nfeat, flags, ntree = struct.unpack("<IIII", c.read(16))
    trees = []
    for _ in range(ntree):
        n = struct.unpack("<I", c.read(4))[0]
        feat = struct.unpack("<%di" % n, c.read(4 * n))
        thr = struct.unpack("<%df" % n, c.read(4 * n))
        left = struct.unpack("<%di" % n, c.read(4 * n))
        right = struct.unpack("<%di" % n, c.read(4 * n))
        val = struct.unpack("<%df" % n, c.read(4 * n))
        dl = c.read(n)
        trees.append(dict(feature=feat, threshold=thr, left=left, right=right, value=val, default_left=dl))
    return nfeat, flags, trees

def to_lgbm(trees, nfeat):
    """LightGBM model text. Internal nodes get indices 0..num_leaves-2 in preorder of OUR
    node array; leaves get ~index (LightGBM encodes a leaf child as -(leaf_idx+1))."""
    out = ["tree", "version=v4", "num_class=1", "num_tree_per_iteration=1", "label_index=0",
           "max_feature_idx=%d" % (nfeat - 1), "objective=regression",
           "feature_names=" + " ".join("f%d" % i for i in range(nfeat)),
           "feature_infos=" + " ".join("[-inf:inf]" for _ in range(nfeat)),
           ""]  # no tree_sizes line: LightGBM then splits on "Tree=" sequentially; a
                #  zero-filled tree_sizes made it mis-locate every tree after the first
    for ti, t in enumerate(trees):
        n = len(t["feature"])
        internal = [i for i in range(n) if t["feature"][i] >= 0]
        leaves = [i for i in range(n) if t["feature"][i] < 0]
        imap = {node: k for k, node in enumerate(internal)}
        lmap = {node: k for k, node in enumerate(leaves)}
        def child(node):
            return imap[node] if t["feature"][node] >= 0 else -(lmap[node] + 1)
        if not internal:  # single-leaf tree
            out += ["Tree=%d" % ti, "num_leaves=1", "num_cat=0", "leaf_value=%.9g" % t["value"][leaves[0]], "shrinkage=1", ""]
            continue
        out += ["Tree=%d" % ti, "num_leaves=%d" % len(leaves), "num_cat=0",
                "split_feature=" + " ".join(str(t["feature"][i]) for i in internal),
                "split_gain=" + " ".join("0" for _ in internal),
                "threshold=" + " ".join("%.9g" % t["threshold"][i] for i in internal),
                # bit1 (=2): default_left; bits 2-3 (=8): missing type NaN
                "decision_type=" + " ".join(str(8 | (2 if t["default_left"][i] else 0)) for i in internal),
                "left_child=" + " ".join(str(child(t["left"][i])) for i in internal),
                "right_child=" + " ".join(str(child(t["right"][i])) for i in internal),
                "leaf_value=" + " ".join("%.9g" % t["value"][i] for i in leaves),
                "leaf_weight=" + " ".join("1" for _ in leaves),
                "leaf_count=" + " ".join("1" for _ in leaves),
                "internal_value=" + " ".join("0" for _ in internal),
                "internal_weight=" + " ".join("1" for _ in internal),
                "internal_count=" + " ".join("1" for _ in internal),
                "is_linear=0", "shrinkage=1", ""]
    out += ["end of trees", "", "feature_importances:", "", "parameters:", "[objective: regression]", "end of parameters", "", "pandas_categorical:null"]
    return "\n".join(out) + "\n"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", default="resources/models/protocol=Automatic/models.parquet")
    ap.add_argument("--model-id", default="hcd_astral_tryp")
    ap.add_argument("--column", default="frag_intensity_model_bytes")
    ap.add_argument("--rows", type=int, default=36864)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    t = pq.ParquetFile(a.store).read()
    ids = t.column("model_id").to_pylist(); blobs = t.column(a.column).to_pylist()
    blob = next(b for i, b in zip(ids, blobs) if i == a.model_id and b)
    nfeat, flags, trees = parse(blob)
    os.makedirs(a.out, exist_ok=True)
    open(os.path.join(a.out, "frag.agbd"), "wb").write(bytes(blob))
    open(os.path.join(a.out, "frag_lgbm.txt"), "w").write(to_lgbm(trees, nfeat))
    # per-feature threshold ranges -> realistic rows
    lo = np.full(nfeat, np.inf); hi = np.full(nfeat, -np.inf)
    for tr in trees:
        for f, th in zip(tr["feature"], tr["threshold"]):
            if f >= 0:
                lo[f] = min(lo[f], th); hi[f] = max(hi[f], th)
    for f in range(nfeat):
        if not np.isfinite(lo[f]): lo[f], hi[f] = 0.0, 1.0
        span = hi[f] - lo[f]; lo[f] -= 0.1 * span + 1e-3; hi[f] += 0.1 * span + 1e-3
    rng = np.random.default_rng(a.seed)
    rows = (lo + rng.random((a.rows, nfeat)) * (hi - lo)).astype("<f4")
    rows.tofile(os.path.join(a.out, "rows.f32"))
    leaves = [sum(1 for f in tr["feature"] if f < 0) for tr in trees]
    meta = dict(model_id=a.model_id, column=a.column, n_trees=len(trees), n_nodes=sum(len(t["feature"]) for t in trees),
                n_features=nfeat, sigmoid=bool(flags & 1), rows=a.rows, leaves_max=max(leaves), leaves_median=float(np.median(leaves)))
    json.dump(meta, open(os.path.join(a.out, "meta.json"), "w"), indent=1)
    print(json.dumps(meta))

if __name__ == "__main__":
    main()
