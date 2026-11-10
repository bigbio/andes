#!/usr/bin/env python3
"""Short comparison of one andes build with the other tools, from the suite's result JSONs: one row per metric.

  python -m quant_eval.compare --results results/ [--build cand] [--arm V] [--tools MQ,FP] [--seed 42] [--only ups27,qe]

One table per dataset. Per metric: the value of the andes arm, and per tool its value and the difference of andes to it (relative to the
tool's value; a factor from 2x on). With several seeds in the results the printed value is the one of `--seed` and the range over the seeds
follows in brackets. Marker: `✅` andes better, `❌` andes worse, `≈` the two values are equal or their ranges over the seeds overlap.
Every value is read from the result JSON as the board reads it; nothing is recomputed here. The two rows on the brightest cells depend on
the reference build of the suite run: its quantities, with those of the tools, rank the cells."""
import argparse
import sys

from .board import g, label, load

# (row label, getter(result, arm), higher is better, format)
METRICS = (
    ("precursors", lambda r, n: g(r, "arms", n, "precursors"), True, "{:,.0f}"),
    ("proteins (>= 2 precursors)", lambda r, n: g(r, "arms", n, "proteins_ge2"), True, "{:,.0f}"),
    ("cells off by > 2x, % (shared cells)", lambda r, n: g(r, "cell_accuracy", "arms", n, "bad"), False, "{:.2f}"),
    ("brightest 10 % of cells off by > 2x, %", lambda r, n: g(r, "cell_accuracy", "arms", n, "bad_bright"), False, "{:.2f}"),
    ("replicate SD, brightest 10 % (log2)", lambda r, n: g(r, "cell_accuracy", "arms", n, "rep_sd_top"), False, "{:.3f}"),
    ("protein CV, %", lambda r, n: g(r, "cv", "protein", n, "all", "med_cv_pct"), False, "{:.1f}"),
    ("DE true calls, common family", lambda r, n: g(r, "arms", n, "de_family", "common", "TP"), True, "{:,.0f}"),
    ("DE false calls, common family", lambda r, n: g(r, "arms", n, "de_family", "common", "FP"), False, "{:,.0f}"),
    ("DE true calls at 10 % false proportion, own rows", lambda r, n: g(r, "arms", n, "de_family", "inter", "tp_at_fdp", "10"), True, "{:,.0f}"),
    ("DE true calls at absolute log2FC >= 0.5, own rows", lambda r, n: g(r, "arms", n, "de_family", "inter", "TP5"), True, "{:,.0f}"),
    ("DE false calls at absolute log2FC >= 0.5, own rows", lambda r, n: g(r, "arms", n, "de_family", "inter", "FP5"), False, "{:,.0f}"),
    ("spiked precursors within 2x, three lowest levels", lambda r, n: low_levels(r, n), True, "{:,.0f}"),
)


def low_levels(r, n):
    ud, lv = g(r, "arms", n, "ups_design"), r.get("levels") or []
    return sum(g(ud, str(L), "t1.0") or 0 for L in lv[:3]) if ud else None


def default_tools(r):
    """the MaxQuant and FragPipe arms that use match between runs (the arms of the common DE family)."""
    return [n for n, k in r["arm_roles"].items() if k["kind"] in ("maxquant", "fragpipe") and not n.lower().replace("_", "").endswith("nombr")]


def span(S, get, n):
    v = [get(S[s], n) for s in sorted(S)]
    v = [x for x in v if x is not None and x == x]
    return (min(v), max(v)) if v else None


def cell(x, rng, f):
    if x is None or x != x:
        return "-"
    return f.format(x) + (f" ({f.format(rng[0])} to {f.format(rng[1])})" if rng and rng[0] != rng[1] else "")


def marker(a, ra, t, rt, higher):
    """`≈` when the two ranges over the seeds overlap (equal values included), else better or worse by the metric's direction."""
    if a is None or t is None or a != a or t != t:
        return "-"
    ra, rt = ra or (a, a), rt or (t, t)
    if ra[0] <= rt[1] and rt[0] <= ra[1]:
        return "≈"
    return "✅" if (a > t) == higher else "❌"


def difference(a, t):
    if a is None or t is None or a != a or t != t:
        return "-"
    if t == 0:
        return f"{a - t:+g}"
    return f"{a / t:.1f}x" if a / t >= 2 else f"{100 * (a - t) / t:+.1f} %"


def table(S, seed, build=None, arm="", tools=None):
    sd = seed if seed in S else sorted(S)[0]
    r = S[sd]
    builds = r.get("builds") or []
    b = build or (builds[-1] if builds else None)
    n = b + ("_" + arm if arm else "") if b else None
    if n not in r["arms"]:
        return f"\n### {r['dataset']}: no arm {n} in the results"
    T = [t for t in (tools or default_tools(r)) if t in r["arms"]]
    seeds = sorted(S)
    out = [f"\n### {r['dataset']}: {label(n, builds)} against {', '.join(label(t, builds) for t in T) or 'no tool'}",
           f"\nSeed {sd}" + (f"; in brackets the range over seeds {', '.join(seeds)}" if len(seeds) > 1 else "")
           + ". ✅ andes better, ❌ andes worse, ≈ equal or the ranges over the seeds overlap. Difference = andes relative to the tool.\n",
           "| metric | " + label(n, builds) + " | " + " | ".join(f"{label(t, builds)} | vs {label(t, builds)}" for t in T) + " |",
           "|---|---|" + "---|---|" * len(T)]
    for name, get, higher, f in METRICS:
        a, ra = get(r, n), span(S, get, n)
        if a is None and all(get(r, t) is None for t in T):
            continue
        row = [name, cell(a, ra, f)]
        for t in T:
            v, rv = get(r, t), span(S, get, t)
            row += [cell(v, rv, f), f"{marker(a, ra, v, rv, higher)} {difference(a, v)}".strip()]
        out.append("| " + " | ".join(row) + " |")
    return "\n".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True)
    ap.add_argument("--build", help="andes build to compare (default: the last build of the results)")
    ap.add_argument("--arm", default="", help="arm of the build: empty (transfer_q <= 0.01), V, H, I or M")
    ap.add_argument("--tools", help="comma-separated tool arms (default: the MaxQuant and FragPipe arms with match between runs)")
    ap.add_argument("--seed", default="42")
    ap.add_argument("--only", help="comma-separated dataset keys (default: every dataset of the results)")
    a = ap.parse_args(argv)
    R = load(a.results)
    if not R:
        raise SystemExit(f"no suite results in {a.results}")
    for dn in (a.only.split(",") if a.only else sorted(R)):
        if dn not in R:
            print(f"skip {dn}: not in the results", file=sys.stderr)
            continue
        print(table(R[dn], a.seed, a.build, a.arm, a.tools.split(",") if a.tools else None))


if __name__ == "__main__":
    main()
