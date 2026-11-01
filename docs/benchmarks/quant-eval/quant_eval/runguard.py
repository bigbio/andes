"""Lineage guard: refuse to join quantities that come from different pipeline runs.

Every source table carries a run tag (hash of the producing file plus a lineage label). `join` asserts equal tags and
one-to-one keys; `assert_cols_agree` adds a value-level check on a column that two sources both carry, so a table labelled as
one run that is not one still fails."""
import hashlib
import os

import numpy as np
import pandas as pd


def file_hash(path, nbytes=8 << 20):
    h = hashlib.sha256()
    h.update(str(os.path.getsize(path)).encode())
    with open(path, "rb") as f:
        h.update(f.read(nbytes))
    return h.hexdigest()[:12]


def run_tag(label, path):
    return f"{label}:{file_hash(path)}"


def tagged(df, tag):
    df = df.copy()
    df.attrs["run_tag"] = tag
    df["run_tag"] = tag
    return df


def tag_of(df):
    t = df.attrs.get("run_tag")
    if t is None and "run_tag" in df.columns:
        u = df["run_tag"].dropna().unique()
        assert len(u) == 1, f"mixed run tags in one table: {u[:5]}"
        t = u[0]
    assert t is not None, "table has no run tag"
    return t


def assert_same_run(*dfs):
    tags = [tag_of(d) for d in dfs]
    assert len(set(tags)) == 1, f"cross-run join refused: {tags}"
    return tags[0]


def assert_cols_agree(a, b, tol=1e-6, log2=False, min_share=0.999, what=""):
    x = np.asarray(a, float)
    y = np.asarray(b, float)
    ok = np.isfinite(x) & np.isfinite(y)
    d = np.abs(np.log2(x[ok] / y[ok])) if log2 else np.abs(x[ok] - y[ok])
    share = float((d <= tol).mean()) if ok.any() else 0.0
    assert share >= min_share, f"value check failed ({what}): only {share:.4f} of {ok.sum()} cells agree within {tol}"
    return share


def join(left, right, on, how="inner", **kw):
    t = assert_same_run(left, right)
    out = left.drop(columns=["run_tag"], errors="ignore").merge(
        right.drop(columns=["run_tag"], errors="ignore"), on=on, how=how, validate="one_to_one", **kw)
    return tagged(out, t)
