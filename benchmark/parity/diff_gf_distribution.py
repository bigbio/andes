#!/usr/bin/env python3
"""Aligns Rust and Java GF score-distribution dumps and reports the first node
where the distributions diverge by more than a tolerance. Emits a localization
hint per the decision tree in
docs/superpowers/specs/2026-05-10-gf-tail-reconciliation-design.md."""
import re, sys, math
from collections import defaultdict

NODE_RE = re.compile(r"GF_NODE: scan=(\d+) pep=(\S+) node_idx=(\d+) mass=(-?\d+) min_score=(-?\d+) max_score=(-?\d+)")
PROB_RE = re.compile(r"GF_PROB: scan=(\d+) pep=(\S+) node_idx=(\d+) score=(-?\d+) prob=(\S+)")
TAIL_RE = re.compile(r"GF_TAIL: scan=(\d+) pep=(\S+) matched_score=(-?\d+) spec_prob=(\S+) tail_sum=(\S+)")

def parse(path):
    nodes = {}
    probs = defaultdict(dict)
    tail = None
    with open(path) as f:
        for line in f:
            m = NODE_RE.search(line)
            if m: nodes[int(m.group(3))] = (int(m.group(4)), int(m.group(5)), int(m.group(6))); continue
            m = PROB_RE.search(line)
            if m: probs[int(m.group(3))][int(m.group(4))] = float(m.group(5)); continue
            m = TAIL_RE.search(line)
            if m: tail = (int(m.group(3)), float(m.group(4)), float(m.group(5)))
    return nodes, probs, tail

def kl(p, q):
    eps = 1e-30
    keys = set(p) | set(q)
    s = 0.0
    for k in keys:
        pk = max(p.get(k, 0.0), eps); qk = max(q.get(k, 0.0), eps)
        s += pk * math.log(pk / qk)
    return s

def main(rust_log, java_log, prob_tol=1e-6, kl_tol=1e-3):
    rn, rp, rt = parse(rust_log)
    jn, jp, jt = parse(java_log)
    print(f"Rust: {len(rn)} nodes, {sum(len(v) for v in rp.values())} prob entries; tail={rt}")
    print(f"Java: {len(jn)} nodes, {sum(len(v) for v in jp.values())} prob entries; tail={jt}")
    if set(rn) != set(jn):
        only_r = sorted(set(rn) - set(jn))[:5]
        only_j = sorted(set(jn) - set(rn))[:5]
        print(f"NODE_INDEX_MISMATCH: rust_only={only_r}  java_only={only_j}")
        print("HINT: graph-construction divergence (PrimitiveAaGraph::new vs Java AminoAcidGraph). Investigate before per-score diff.")
        return 2
    print("\nPer-node first-divergence scan:")
    print(f"{'node':>5}  {'mass':>6}  {'r_range':>11}  {'j_range':>11}  {'KL':>9}  {'max|dP|':>9}  flag")
    first_div = None
    for ni in sorted(rn):
        rmass, rmin, rmax = rn[ni]
        jmass, jmin, jmax = jn[ni]
        rdist, jdist = rp.get(ni, {}), jp.get(ni, {})
        d = kl(rdist, jdist)
        max_dp = max((abs(rdist.get(k, 0) - jdist.get(k, 0)) for k in (set(rdist) | set(jdist))), default=0.0)
        flag = ""
        if (rmin, rmax) != (jmin, jmax): flag = "RANGE_DIFF"
        elif d > kl_tol: flag = "DIST_DIFF"
        elif max_dp > prob_tol: flag = "PROB_DIFF"
        if flag and first_div is None: first_div = (ni, flag)
        print(f"{ni:>5}  {rmass:>6}  ({rmin:>3},{rmax:>3})  ({jmin:>3},{jmax:>3})  {d:>9.3e}  {max_dp:>9.3e}  {flag}")
    if first_div:
        ni, flag = first_div
        print(f"\nFIRST_DIVERGENCE: node_idx={ni}  flag={flag}")
        print("LIKELY ROOT CAUSE (per spec decision tree):")
        if flag == "RANGE_DIFF":
            print("  -> mass-bin window rounding (PrimitiveAaGraph::new mass extents) OR graph topology")
        elif flag == "DIST_DIFF":
            print("  -> edge prob accumulation (add_prob_dist) OR AA probability lookup OR upstream node carrying error")
        elif flag == "PROB_DIFF":
            print("  -> small per-score precision drift (f32 vs f64) -- check the cumulative pattern in following nodes")
        return 1
    if rt and jt and abs(math.log10(max(rt[1], 1e-300)) - math.log10(max(jt[1], 1e-300))) > 0.3:
        print("\nALL_NODES_MATCH but TAIL_MISMATCH:")
        print("LIKELY ROOT CAUSE: underflow guard or score-range cutoff at max_score (compute_inner / spectral_probability)")
        return 1
    print("\nALL_MATCH within tolerances.")
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
