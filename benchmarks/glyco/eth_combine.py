#!/usr/bin/env python3
# Pool 3 fractions into one PIN (frac-tagged SpecId) for a STABLE @1% Percolator
# (per-fraction had only 0-2 decoys => noisy). Concatenate existing per-frac PINs.
import sys
mode = sys.argv[1]  # "concat" or "eval"

# per-weight fraction PIN files
PINS = {
    "0": {1: "eth_czoff.glyco.pin", 2: "sc_f2_off.glyco.pin", 3: "sc_f3_off.glyco.pin"},
    "2": {1: "wt_f1_w2.glyco.pin", 2: "wt_f2_w2.glyco.pin", 3: "wt_f3_w2.glyco.pin"},
    "3": {1: "wt_f1_w3.glyco.pin", 2: "wt_f2_w3.glyco.pin", 3: "wt_f3_w3.glyco.pin"},
    "5": {1: "eth_czon.glyco.pin", 2: "sc_f2_on.glyco.pin", 3: "sc_f3_on.glyco.pin"},
}

if mode == "concat":
    W = sys.argv[2]
    out = open(f"comb_w{W}.pin", "w")
    wrote_hdr = False
    for frac, pin in PINS[W].items():
        with open(pin) as fh:
            hdr = fh.readline()
            if not wrote_hdr:
                out.write(hdr); wrote_hdr = True
            for ln in fh:
                # tag SpecId (col0) with fraction so scan numbers don't collide
                i = ln.find("\t")
                out.write(f"f{frac}_" + ln)
    out.close()
    print(f"wrote comb_w{W}.pin")

elif mode == "eval":
    W = sys.argv[2]; H2O = 18.010565
    # combined truth keyed by (frac, scan)
    truth = {}
    for frac in (1, 2, 3):
        for ln in open(f"ethcd_truth_frac{frac}.tsv"):
            if ln.startswith("scan"): continue
            c = ln.rstrip("\n").split("\t")
            truth[(frac, int(c[0]))] = (float(c[4]), int(c[2]))
    pin = f"comb_w{W}.pin"; psms = f"comb_w{W}.psms"
    hdr = None; info = {}
    for ln in open(pin):
        f = ln.rstrip("\n").split("\t")
        if hdr is None:
            hdr = f; iS = hdr.index("ScanNr"); iC = hdr.index("CalcMass"); iG = hdr.index("GlycanMass"); continue
        sid = f[0]  # fN_controllerType=...scan=...
        try:
            frac = int(sid[1]); scan = int(float(f[iS])); calc = float(f[iC]); gly = float(f[iG])
        except Exception:
            continue
        info[sid] = (frac, scan, calc, gly)
    from collections import defaultdict
    corr = defaultdict(set)
    ntgt = 0
    for ln in open(psms):
        f = ln.rstrip("\n").split("\t")
        if f[0] == "PSMId" or len(f) < 3: continue
        try:
            if float(f[2]) > 0.01: continue
        except Exception:
            continue
        ntgt += 1
        if f[0] not in info: continue
        frac, scan, calc, gly = info[f[0]]
        k = (frac, scan)
        if k in truth:
            tb, ch = truth[k]
            if abs(calc - gly - H2O - tb) < 0.1 or abs(calc - gly - tb) < 0.1:
                corr[ch].add(k)
    tot = sum(len(v) for v in corr.values())
    ttot = defaultdict(int)
    for (fr, s), (tb, ch) in truth.items(): ttot[ch] += 1
    print(f"W={W}: target@1%={ntgt}  BACKBONE-CORRECT={tot} / {len(truth)}  " +
          "  ".join(f"z{ch}:{len(corr[ch])}/{ttot[ch]}" for ch in sorted(ttot) if ttot[ch] > 3))
