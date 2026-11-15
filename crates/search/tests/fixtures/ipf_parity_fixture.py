#!/usr/bin/env python3
"""Generate the IPF parity fixture for crates/search/tests/ipf_pyprophet_parity.rs.

Runs PyProphet's own `apply_bm` (log space and linear space) and `compute_model_fdr`,
copied verbatim from pyprophet/ipf.py (loguru-free), on random PSMs: 1-6 peptidoforms,
0-12 ions each mapping to a random non-empty proper subset of forms, plus h0.

Output (JSON): [{"psm_pep", "n_forms", "evidence": [[pep, [forms]]], "posterior": [...forms, h0]}],
and {"peps": [...], "qvalues": [...]} for compute_model_fdr.
Usage: python3 ipf_parity_fixture.py PYPROPHET_SRC_DIR OUT.json
"""
import json, random, sys, ast
import numpy as np, pandas as pd
from scipy.special import logsumexp
from scipy.stats import rankdata

src = open(f"{sys.argv[1]}/pyprophet/ipf.py").read()
tree = ast.parse(src)
ns = {"np": np, "pd": pd, "logsumexp": logsumexp, "rankdata": rankdata}
for node in tree.body:
    if isinstance(node, ast.FunctionDef) and node.name in ("apply_bm", "compute_model_fdr"):
        exec(compile(ast.Module([node], []), "ipf.py", "exec"), ns)
apply_bm, compute_model_fdr = ns["apply_bm"], ns["compute_model_fdr"]

rng = random.Random(20261113)
cases = []
for fid in range(300):
    n = rng.randint(1, 6)
    psm_pep = rng.choice([0.0, 1e-6, rng.random() * 0.2, rng.random()])
    ev = []
    for t in range(rng.randint(0, 12)):
        k = rng.randint(1, max(1, n - 1))
        forms = sorted(rng.sample(range(n), k)) if n > 1 else [0]
        ev.append([rng.choice([1e-4, rng.random() * 0.1, rng.random()]), forms])
    rows = []
    hyps = list(range(n)) + [-1]
    # prepare_transition_bm: prior per hypothesis, evidence per (transition, hypothesis)
    prior = {h: (psm_pep if h == -1 else (1 - psm_pep) / n) for h in hyps}
    if not ev:  # a feature with no transitions still has its priors
        for h in hyps:
            rows.append(dict(feature_id=fid, hypothesis=h, prior=prior[h], evidence=1.0))
    for pep, forms in ev:
        for h in hyps:
            bm = h != -1 and h in forms
            rows.append(dict(feature_id=fid, hypothesis=h, prior=prior[h],
                             evidence=(1 - pep) if bm else pep))
    df = pd.DataFrame(rows)
    out = {}
    for space, flag in (("log", True), ("linear", False)):
        pp = apply_bm(df, use_log_space=flag).set_index("hypothesis")["posterior"]
        out[space] = [float(pp.loc[h]) for h in hyps]
    cases.append({"psm_pep": psm_pep, "n_forms": n, "evidence": ev,
                  "posterior_log": out["log"], "posterior_linear": out["linear"]})

peps = [round(rng.random(), 3) for _ in range(200)] + [0.5] * 5 + [0.0] * 3
q = compute_model_fdr(np.array(peps)).tolist()
json.dump({"cases": cases, "fdr": {"peps": peps, "qvalues": q}}, open(sys.argv[2], "w"))
print("cases", len(cases), "fdr", len(peps))
