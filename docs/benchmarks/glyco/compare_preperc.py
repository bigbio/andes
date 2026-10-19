#!/usr/bin/env python3
"""
Pre-Percolator A/B comparison of .glyco.pin files.

No Percolator. Reports, per run:
  * counts (total / target / decoy / unique scans)
  * RawScore target-vs-decoy separation (AUC = P(target > decoy), plus stats)
  * decoy-estimated FDR yield (2x rule, 1:1 entrapment) on RawScore
  * entrapment (ENTRAP_) rate at the 1% FDR threshold
  * per-feature AUC (target > decoy), top discriminators

Usage: compare_preperc.py <label:pin> [<label:pin> ...]
"""
import sys
import csv


def load(pin):
    with open(pin) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def auc_target_gt_decoy(ts, ds):
    """AUC = P(random target > random decoy), tie-split. 0.5 = random, 1 = perfect."""
    ts = sorted(v for v in ts if v is not None)
    ds = sorted(v for v in ds if v is not None)
    n, m = len(ts), len(ds)
    if n == 0 or m == 0:
        return float("nan"), 0, 0
    # rank each decoy against targets: count targets strictly greater, tie-split
    total = 0.0
    for d in ds:
        # number of targets > d
        import bisect
        gt = n - bisect.bisect_right(ts, d)   # targets strictly > d
        eq = bisect.bisect_right(ts, d) - bisect.bisect_left(ts, d)  # targets == d
        total += gt + 0.5 * eq
    return total / (n * m), n, m


def decoy_fdr_yield(rows):
    scored = [r for r in rows if num(r["RawScore"]) is not None]
    scored.sort(key=lambda r: num(r["RawScore"]), reverse=True)
    t = d = 0
    best = 0
    kept_rows = []
    for r in scored:
        if r["Label"] == "1":
            t += 1
        else:
            d += 1
        fdr = 2 * d / (d + t) if (d + t) else 0.0
        if fdr <= 0.01:
            best = t
            kept_rows = scored[: t + d]
    ent_at_thr = sum(1 for r in kept_rows if r["Label"] == "1" and "ENTRAP_" in r["Proteins"])
    return best, kept_rows, ent_at_thr


def pct(a, q):
    if not a:
        return float("nan")
    a = sorted(a)
    return a[min(len(a) - 1, int(len(a) * q))]


def summarize(name, pin):
    rows = load(pin)
    tgt = [r for r in rows if r["Label"] == "1"]
    dec = [r for r in rows if r["Label"] == "-1"]
    ts = [num(r["RawScore"]) for r in tgt]
    ds = [num(r["RawScore"]) for r in dec]
    ts_v = [v for v in ts if v is not None]
    ds_v = [v for v in ds if v is not None]

    a, nt, nd = auc_target_gt_decoy(ts_v, ds_v)
    best, kept, ent_thr = decoy_fdr_yield(rows)

    print(f"\n===== {name} =====")
    print(f"  file           : {pin}")
    print(f"  total PSMs     : {len(rows)}  (target {len(tgt)}, decoy {len(dec)})")
    print(f"  unique scans   : {len(set(r['ScanNr'] for r in rows))}")
    print(f"  --- RawScore target vs decoy ---")
    print(f"    target n={nt:6d}  min={min(ts_v):.3f}  med={pct(ts_v,.5):.3f}  mean={sum(ts_v)/len(ts_v):.3f}  p90={pct(ts_v,.9):.3f}  max={max(ts_v):.3f}")
    print(f"    decoy  n={nd:6d}  min={min(ds_v):.3f}  med={pct(ds_v,.5):.3f}  mean={sum(ds_v)/len(ds_v):.3f}  p90={pct(ds_v,.9):.3f}  max={max(ds_v):.3f}")
    print(f"    RawScore AUC (target>decoy) = {a:.4f}")
    print(f"    targets > max(decoy) = {sum(1 for v in ts_v if v > max(ds_v))}")
    print(f"  --- decoy-FDR (2x rule, RawScore desc) ---")
    print(f"    targets @ <=1% FDR = {best}")
    print(f"    ENTRAP_ among those targets = {ent_thr}  ({100*ent_thr/max(1,best):.2f}%)")

    # per-feature AUC: which single feature best separates target/decoy?
    feats = ["RawScore", "RankScore", "RankScoreFloat", "OxoniumScore", "YLadderScore",
             "YHitFrac", "CoreYHits", "GlycanMass", "SialicConsistency", "matchedIonRatio",
             "MeanErrorTop7", "EdgeScore", "DeltaRankScore", "TailorScore", "IntensitySignal",
             "NumMatchedMainIons", "longest_y_pct", "ExplainedIonCurrentRatio"]
    aucs = []
    for f_ in feats:
        if f_ not in rows[0]:
            continue
        tv = [num(r[f_]) for r in tgt]
        dv = [num(r[f_]) for r in dec]
        av, _, _ = auc_target_gt_decoy(tv, dv)
        # AUC can be inverted; report the max distance from 0.5
        aucs.append((abs(av - 0.5), av, f_))
    aucs.sort(reverse=True)
    print(f"  --- top single-feature separators (|AUC-0.5| desc) ---")
    for _, av, f_ in aucs[:10]:
        print(f"    {f_:24s} AUC(target>decoy)={av:.3f}")


def main():
    args = sys.argv[1:]
    if len(args) < 1:
        print(__doc__)
        sys.exit(1)
    for a in args:
        label, _, pin = a.partition(":")
        summarize(label, pin)


if __name__ == "__main__":
    main()
