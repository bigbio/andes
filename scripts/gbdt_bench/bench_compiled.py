#!/usr/bin/env python3
"""Time compiled tree-ensemble evaluators on the exported model and compare to the Rust
reference. Each evaluator is optional (skipped with a note if its package is missing).

  lleaves   - LLVM-compiled LightGBM trees (pip install lleaves)
  tl2cgen   - treelite -> C -> shared library (pip install treelite tl2cgen; needs a C compiler)
  lightgbm  - the reference LightGBM predictor (pip install lightgbm), as a sanity baseline

Reports ns/row for per-PSM-sized batches (48 rows) and for the whole set, plus the max
absolute deviation from preds_rust.f32. A deviation above 1e-4 means the LightGBM export
is wrong, not that the evaluator is inexact.
"""
import argparse, json, os, sys, time
import numpy as np

def timeit(fn, rows, batch, reps=5):
    n = rows.shape[0]; best = float("inf")
    for _ in range(reps):
        t = time.perf_counter()
        if batch >= n:
            fn(rows)
        else:
            for i in range(0, n, batch):
                fn(rows[i:i + batch])
        best = min(best, time.perf_counter() - t)
    return best * 1e9 / n

def report(name, fn, rows, ref):
    pred = np.asarray(fn(rows), dtype=np.float64).ravel()
    dev = float(np.max(np.abs(pred - ref)))
    r48 = timeit(fn, rows, 48); rall = timeit(fn, rows, rows.shape[0])
    print(f"{name:10s} batch=48: {r48:8.1f} ns/row   batch=all: {rall:8.1f} ns/row   max|dev vs rust|={dev:.2e}")

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--dir", required=True); a = ap.parse_args()
    meta = json.load(open(os.path.join(a.dir, "meta.json")))
    nf = meta["n_features"]
    rows = np.fromfile(os.path.join(a.dir, "rows.f32"), dtype="<f4").reshape(-1, nf)
    refp = os.path.join(a.dir, "preds_rust.f32")
    if not os.path.exists(refp):
        sys.exit("run the Rust bench first (writes preds_rust.f32)")
    ref = np.fromfile(refp, dtype="<f4").astype(np.float64)
    print(f"model {meta['model_id']}: {meta['n_trees']} trees, {meta['n_nodes']} nodes, {nf} features; rows {rows.shape[0]}")
    txt = os.path.join(a.dir, "frag_lgbm.txt")
    try:
        import lightgbm as lgb
        b = lgb.Booster(model_file=txt)
        report("lightgbm", lambda x: b.predict(x, num_threads=1), rows, ref)
    except Exception as e:
        print(f"lightgbm: skipped ({type(e).__name__}: {str(e)[:120]})")
    try:
        import lleaves
        m = lleaves.Model(model_file=txt); m.compile()
        report("lleaves", lambda x: m.predict(x, n_jobs=1), rows, ref)
    except Exception as e:
        print(f"lleaves: skipped ({type(e).__name__}: {str(e)[:120]})")
    try:
        import treelite, tl2cgen
        model = treelite.frontend.load_lightgbm_model(txt) if hasattr(treelite.frontend, "load_lightgbm_model") else treelite.Model.load(txt, model_format="lightgbm")
        lib = os.path.join(a.dir, "tl2cgen_model.so")
        tl2cgen.export_lib(model, toolchain="gcc", libpath=lib, params={"parallel_comp": 8}, verbose=False)
        pred = tl2cgen.Predictor(lib, nthread=1)
        report("tl2cgen", lambda x: pred.predict(tl2cgen.DMatrix(x)), rows, ref)
    except Exception as e:
        print(f"tl2cgen: skipped ({type(e).__name__}: {str(e)[:120]})")

if __name__ == "__main__":
    main()
