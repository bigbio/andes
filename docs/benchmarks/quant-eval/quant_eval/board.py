#!/usr/bin/env python3
"""ONE scoreboard table from the suite result JSONs (suite.py): per set and tool (arm).

  python -m quant_eval.board --results results/ [--seed 42] [--build main] > scoreboard.md

Columns (seed `--seed`; DE also as the range over all seeds present):
  prec / prot         precursors and proteins with >= 2 precursors; rows identified / transferred
  acc                 identified cells, every arm on the same cells: share |err| > 1 log2, within-2x share, bright-decile gross share (accuracy.py)
  prec.SD / CV        replicate SD of the top predicted-intensity decile (log2); protein replicate CV, median % (matched features)
  L1                  UPS1 precursors within 2x of their design reference, summed over the three lowest levels
  absence             PXD014415: yeast share of the transferred rows in the 20 human-only runs; single-cell set: rows in the Blank run
  wrong-peak          masked transfers off the identified apex at the arm's quant_q cut, pooled (level B) (andes hard-cut arms)
  DE                  intersection family, BH 5 % per pair: TP / FP without a fold-change cut, with |log2FC| >= 0.5, TP at pooled FDP 10 %
New columns (right of the legacy ones; `--legacy` prints the legacy columns only; `--detail` adds the per-level, absence-class and blank tables):
  DE columns          inter, signed: TP with the expected sign / FP plus wrong-sign calls (wrong-sign calls in brackets), intersection family;
                      common: the same on the family of proteins testable in every match-between-runs arm (identical tested rows), unsigned and signed;
                      tested rows: (protein, pair) rows each arm tests in its own / the intersection / the common family. The legacy DE column is the
                      intersection family, unsigned, on each arm's own testable rows.
  within 2x %         dropped from the default board (it is 100 - err>1 %); `--legacy` keeps it
  L1 levels           the amol levels summed in L1 (the three lowest); a level that is also a reference level is flagged `(self-ref)`
  L1 rows (+ unscored)   UPS1 rows the arm reports at those levels for precursors with a design reference (+ the rows of precursors without one,
                      which leave the count and the denominator of the scored share); within-2x % of the scored rows, and with the unscored rows counted as not within 2x
  L1 med abs err      median |log2 error| of the scored rows
  L1 @ matched rows   L1 count when every arm is limited, per level, to the smallest row count among the plain arms (andes build, MaxQuant, FragPipe).
                      Headline: the expected count for the same number of rows drawn at random (order-free). Second column: best rows first, which orders andes
                      by quality column and the arms without one by intensity, so the order differs per arm and is unfavourable to the latter
  absence             PXD014415: noise-removed (lower bound) rate of the species' rows in the runs without it, transfers / identified: yeast-enriched rows, and in
                      brackets enriched + only-absent rows (no yeast-run reference). Classes come from one reference table (the andes every-row table) for every arm
  blank id / tr / carry   single-cell set: identified and transferred rows in the blank run and the rows below the learned carry-over fraction
Arm names: `X` transfer_q <= 0.01, `X_V` all accepted transfers with learned precision weights, `X_H` hard cut quant_q <= 0.01, `X_I` identified rows only (no match between runs), `X_M` hard cut at the
matched true-error cut of the reference build (matched.py). Tools are context only."""
import argparse
import glob
import json
import os

from . import lowlevel as LL

NAME = {"MQ": "MaxQuant (MBR)", "MQ_noMBR": "MaxQuant no MBR", "FP": "FragPipe", "Sage": "Sage"}
# labels of the default board: the context arms say what they are (FlashLFQ runs per peptide, Sage has no transfer flag)
NAME_NEW = {"Sage": "Sage (no transfer flag: MBR fills counted as identified)"}
_STATE = {"legacy": False}


def load(d):
    R = {}
    for f in sorted(glob.glob(os.path.join(d, "*_s*.json"))):
        j = json.load(open(f))
        if "suite_dataset" in j:
            R.setdefault(j["suite_dataset"], {})[j["seed"]] = j
    return R


def g(d, *k):
    for x in k:
        d = d.get(x) if isinstance(d, dict) else None
        if d is None:
            return None
    return d


def fmt(x, f="{:.2f}", none="-"):
    return none if x is None else f.format(x)


def de(A, n, key, seed_range=None):
    f = g(A, n, "de_family", "inter")
    if not f:
        return "-"
    return f"{f[key[0]]:,} / {f[key[1]]:,}"


def label(n, builds):
    for b in builds:
        if n == b:
            return f"andes {b}"
        for suf, nm in (("_V", "V"), ("_H", "H"), ("_I", "no MBR"), ("_M", "M (matched)")):
            if n == b + suf:
                return f"andes {b} {nm}"
    if not _STATE["legacy"]:
        if n.startswith("FL"):
            return f"{n} (FlashLFQ, peptide level)"
        return NAME_NEW.get(n, NAME.get(n, n))
    return NAME.get(n, n)


NEW_COLS = (" L1 levels (amol) | L1 rows (+ unscored) | L1 within-2x % | L1 within-2x %, unscored counted as not within | L1 med abs err |"
            " L1 @ matched rows, random rows (headline) | L1 @ matched rows, best first (order differs per arm) |"
            " absence, noise removed (lower bound): enriched (enriched + only-absent), tr / id | blank id / tr / carry |"
            " DE inter, signed: TP / false (wrong sign) | DE common: TP / FP | DE common, signed: TP / false (wrong sign) | DE tested rows own / inter / common |")
N_NEW = 13


def plain_arms(r, builds):
    """Arms of the matched-row reading: the andes builds' default arms, MaxQuant (MBR) and FragPipe."""
    return [n for n in r["arms"] if n in builds or n in ("MQ", "FP")]


def rate_pct(num, den, f):
    """num / den in percent, '-' when the denominator is zero or missing."""
    return "-" if not den else f.format(100 * num / den)


def absence_cell(o):
    """Noise-removed absence rate of one absent species: enriched (enriched + only-absent) for transferred / identified rows, '%'."""
    dn = o["denominators"]
    out = []
    for kind, den, f in (("transferred", dn["transferred_rows_absent_runs"], "{:.3f}"), ("identified", dn["identified_rows_absent_runs"], "{:.4f}")):
        c = o["learned"][kind]
        out.append(f"{rate_pct(c['yeast_enriched'], den, f)} ({rate_pct(c['yeast_enriched'] + c['only_absent'], den, f)})")
    return f"{out[0]} / {out[1]} %"


def new_cells(r, n, builds, matched):
    d = g(r, "arms", n, "l1_detail") or {}
    c = ["-"] * N_NEW
    if d:
        l1 = d["l1_levels"]
        c[0] = "/".join(l1) + (" (self-ref)" if d["self_reference"] else "")
        p = d["l1_pooled"]
        c[1] = f"{p['rows']:,} (+{p.get('unscored_rows', 0):,})"
        c[2] = "-" if p["share_within2x"] is None else f"{100 * p['share_within2x']:.1f}"
        c[3] = "-" if p.get("within2x_incl_unscored") is None else f"{100 * p['within2x_incl_unscored']:.1f}"
        c[4] = fmt(p["median_abs_err"])
        if matched is not None and n in matched:
            c[5] = f"{matched[n][1]:.1f}"
            c[6] = f"{matched[n][0]:,}"
    ac = (r.get("absence_classes") or {}).get(n) or {}
    for s_, o in ac.items():
        if "learned" in o:
            c[7] = absence_cell(o)
    bs = (r.get("blank_split") or {}).get(n)
    if bs:
        c[8] = f"{bs['identified']['rows']} / {bs['transferred']['rows']} / {bs['carry_over_learned']}"
    fi, fc = g(r, "arms", n, "de_family", "inter"), g(r, "arms", n, "de_family", "common")
    if fi and "WS" in fi:
        c[9] = f"{fi['TPs']:,} / {fi['FP'] + fi['WS']:,} ({fi['WS']})"
    if fc:
        c[10] = f"{fc['TP']:,} / {fc['FP']:,}"
        c[11] = f"{fc['TPs']:,} / {fc['FP'] + fc['WS']:,} ({fc['WS']})"
    fo = g(r, "arms", n, "de_family", "own")
    if fo and fi and fc:
        c[12] = f"{fo['tested']:,} / {fi['tested']:,} / {fc['tested']:,}"
    return c


def matched_l1(r, builds):
    """{arm: L1 count at the matched rows} of one result (sum over the L1 levels of the keys within 2x among the first k rows)."""
    det = {n: g(r, "arms", n, "l1_detail") for n in r["arms"]}
    det = {n: d for n, d in det.items() if d}
    if not det:
        return None
    ref = [n for n in plain_arms(r, builds) if n in det]
    tot = {n: [0, 0.0] for n in det}
    for L in next(iter(det.values()))["l1_levels"]:
        k, cs = LL.matched_reading(det, ref, L)
        for n, c in cs.items():
            tot[n][0] += c["keys"]
            tot[n][1] += c["keys_random"]
    return tot


def build_rows(R, seed, builds, order, legacy=False):
    hdr = ("| set | arm | prec | prot | rows id / tr | err>1 % | within 2x % | bright gross % | prec.SD top | prot CV % | L1 | absence | wrong-peak % (B) | DE TP / FP | DE FC>=0.5 | TP@FDP10 | DE TP/FP range over seeds |\n"
           "|" + "---|" * 17)
    _STATE["legacy"] = legacy
    if not legacy:
        h = hdr.split("\n")[0].replace(" within 2x % |", "", 1).replace("| DE TP / FP |", "| DE TP / FP (inter, unsigned, arm's own testable rows) |", 1)
        hdr = h + NEW_COLS + "\n|" + "---|" * (16 + N_NEW)
    out = [hdr]
    for dn in order:
        if dn not in R:
            continue
        S = R[dn]
        sd = seed if seed in S else sorted(S)[0]
        r = S[sd]
        A = r["arms"]
        mt = matched_l1(r, builds)
        ca = r.get("cell_accuracy") or {}
        for n in A:
            a = A[n]
            plain = n in builds or r["arm_roles"][n]["kind"] != "andes"
            acc = (ca.get("arms") or {}).get(n) if plain else None
            ud = a.get("ups_design") or {}
            lv = r.get("levels") or []
            l1 = sum(g(ud, str(L), "t1.0") or 0 for L in lv[:3]) if ud else None
            ab = None
            if dn_is_p14(r):
                ab = g(r, "rates", n, "YEAST", "rate")
                ab = None if ab is None else f"{100 * ab:.2f} %"
            elif r.get("carry_over"):
                co = r["carry_over"].get(n)
                ab = None if not co else f"{co['rows']} rows"
            wp = "-"
            for b in builds:
                if n in (b + "_H", b + "_M"):
                    cu = g(r, "wrong_peak", b)
                    tau = 0.01 if n.endswith("_H") else g(r, "matched_cuts", b, "tau")
                    if cu and tau:
                        p = min(cu, key=lambda x: abs(x["tau"] - tau))
                        wp = f"{100 * p['wrong'] / max(p['n'], 1):.2f} ({100 * p['B']['wrong'] / max(p['B']['n'], 1):.2f})" + (f" @{tau:g}" if n.endswith("_M") else "")
            f = g(A, n, "de_family", "inter")
            rng = "-"
            if f:
                tps = [g(S[s], "arms", n, "de_family", "inter", "TP") for s in S]
                fps = [g(S[s], "arms", n, "de_family", "inter", "FP") for s in S]
                tps, fps = [x for x in tps if x is not None], [x for x in fps if x is not None]
                if len(tps) > 1:
                    rng = f"{min(tps):,}-{max(tps):,} / {min(fps):,}-{max(fps):,}"
            cvp = g(r, "cv", "protein", n, "all", "med_cv_pct")
            trs = f"{a['transfer_rows']:,}"
            if not legacy and r["arm_roles"][n]["kind"] == "sage":
                trs = "n/a"
            if not legacy and ab and ab.startswith("nan"):
                ab = None
            w2 = [] if not legacy else [fmt(acc and acc["within2x"])]   # within 2x % = 100 - err>1 %: dropped from the default board
            out.append("| " + " | ".join([
                r["dataset"] + (f" s{sd}" if False else ""), label(n, builds), f"{a['precursors']:,}", f"{a['proteins_ge2']:,}", f"{a['ident_rows']:,} / {trs}",
                fmt(acc and acc["bad"])] + w2 + [fmt(acc and acc["bad_bright"]), fmt(acc and acc.get("rep_sd_top"), "{:.3f}"), fmt(cvp, "{:.1f}"),
                fmt(l1, "{:,}"), ab or "-", wp, de(A, n, ("TP", "FP")), de(A, n, ("TP5", "FP5")),
                (f"{f['tp_at_fdp']['10']:,}" if f else "-"), rng] + ([] if legacy else new_cells(r, n, builds, mt))) + " |")
    return "\n".join(out)


def _pct(x, f="{:.1f}"):
    return "-" if x is None else f.format(100 * x)


def detail_tables(R, seed, builds, order, arms_only=None):
    """Per-level L1 table with precision and matched rows, A2 at the lowest level, absence classes and blank split, per set."""
    out = []
    for dn in order:
        if dn not in R:
            continue
        S = R[dn]
        r = S[seed if seed in S else sorted(S)[0]]
        A = r["arms"]
        det = {n: g(r, "arms", n, "l1_detail") for n in A}
        det = {n: d for n, d in det.items() if d}
        names = [n for n in A if (arms_only is None or n in arms_only or n in builds)]
        if det:
            d0 = next(iter(det.values()))
            ref = [n for n in plain_arms(r, builds) if n in det]
            out.append(f"\n### {r['dataset']}: UPS1 levels (amol), L1 = {'/'.join(d0['l1_levels'])}; references built from {'/'.join(d0['ref_levels'])}; "
                       f"self-reference: {'/'.join(d0['self_reference']) or 'none'}; row order: andes and quality-column arms by quality, others by intensity\n")
            out.append("| level | role | arm | rows | keys within 2x | rows within 2x % | med abs err | med err | matched rows | keys @ matched, random rows (headline) | keys @ matched, best first (order differs per arm) | rows within 2x % @ matched, best first |\n|" + "---|" * 12)
            for L in d0["levels"]:
                k, mc = LL.matched_reading(det, ref, L) if L in d0["l1_levels"] else (0, {})
                for n in names:
                    if n not in det:
                        continue
                    p = det[n]["per_level"][L]
                    c = mc.get(n)
                    out.append("| " + " | ".join([L, p["role"] + (" (self-ref)" if L in d0["self_reference"] else ""), label(n, builds), f"{p['rows']:,}", f"{p['keys_within2x']:,}",
                                                  _pct(p["share_within2x"]), fmt(p["median_abs_err"]), fmt(p["median_err"], "{:+.2f}"),
                                                  f"{k:,}" if c else "-", f"{c['keys_random']:.1f}" if c else "-", f"{c['keys']:,}" if c else "-", _pct(c["share_within2x"]) if c else "-"]) + " |")
            L0 = d0["levels"][0]
            out.append(f"\nA2 at the lowest level ({L0} amol): error of the transferred rows (log2 against the design reference), identified rows next to them\n")
            out.append("| arm | transfer rows | median err | median abs err | within 2x % | identified rows | median err | within 2x % |\n|" + "---|" * 8)
            for n in names:
                if n in det:
                    p = det[n]["per_level"][L0]
                    t, i = p["transfer"], p["identified"]
                    out.append("| " + " | ".join([label(n, builds), f"{t['rows']:,}", fmt(t["median_err"], "{:+.2f}"), fmt(t["median_abs_err"]), _pct(t["share_within2x"]),
                                                   f"{i['rows']:,}", fmt(i["median_err"], "{:+.2f}"), _pct(i["share_within2x"])]) + " |")
        ac = r.get("absence_classes")
        if ac:
            for sp_ in sorted({s_ for a in ac.values() for s_ in a}):
                first = next(a[sp_] for a in ac.values() if sp_ in a and "cutoffs" in a[sp_])
                ref = first.get("reference") or first
                src = f"one reference table ({(r.get('absence_reference') or {}).get('table', '?')}, {(r.get('absence_reference') or {}).get('rows', 0):,} rows) for every arm" if first.get("reference") else "each arm's own rows"
                out.append(f"\n### {r['dataset']}: {sp_} rows in the runs without it, by class of the precursor's cross-run profile from {src} (control {ref['control']}, "
                           f"boundary = {100 * (1 - ref['alpha']):g} % quantile of the control's d per bin of absent-run rows; fixed boundary {ref['fixed_log2']:g} log2). "
                           "Yeast-enriched is a lower bound of the species-specific errors.\n")
                out.append("Learned boundary per bin of absent-run rows (n control precursors; * = fewer than 100, pooled cutoff used): "
                           + "; ".join(f"{x['n_abs_from']}+: {x['cutoff']:.2f} (n={x['n_control']:,}{'*' if x.get('pooled_fallback') else ''})" for x in ref["cutoffs"]) + "\n")
                out.append("| arm | class source | transfers: enriched / not specific / only absent / not in reference | rate all % | rate enriched % | rate enriched + only-absent % | identified: enriched / not specific / only absent / not in reference | rate all % | rate enriched % | rate enriched + only-absent % |\n|" + "---|" * 10)
                for bl in ("learned", "fixed", "own_learned", "own_fixed"):
                    for n in names:
                        o = (ac.get(n) or {}).get(sp_)
                        if not o or bl not in o:
                            continue
                        dn = o["denominators"]
                        cells = [label(n, builds), bl]
                        for kind, den in (("transferred", dn["transferred_rows_absent_runs"]), ("identified", dn["identified_rows_absent_runs"])):
                            c = o[bl][kind]
                            f = "{:.4f}" if kind == "identified" else "{:.3f}"
                            cells += [f"{c['yeast_enriched']:,} / {c['not_specific']:,} / {c['only_absent']:,} / {c.get('not_in_reference', 0):,}",
                                      rate_pct(sum(c.values()), den, f), rate_pct(c["yeast_enriched"], den, f), rate_pct(c["yeast_enriched"] + c["only_absent"], den, f)]
                        out.append("| " + " | ".join(cells) + " |")
        bs = r.get("blank_split")
        if bs:
            b0 = next(b for b in bs.values() if b)
            out.append(f"\n### {r['dataset']}: blank-run rows. Carry-over-like = raw intensity below the learned fraction {b0['learned_fraction']:.4f} of the precursor's median in the other runs "
                       f"(the {100 * b0['quantile']:g} % quantile of ordinary rows); fixed fraction {b0['fixed_fraction']:g} next to it\n")
            out.append("| arm | blank rows | identified | of which carry-over learned / fixed / no reference | transferred | of which carry-over learned / fixed / no reference |\n|" + "---|" * 6)
            for n in names:
                b = bs.get(n)
                if b:
                    i, t = b["identified"], b["transferred"]
                    out.append("| " + " | ".join([label(n, builds), f"{b['rows']:,}", f"{i['rows']:,}", f"{i['carry_over_learned']} / {i['carry_over_fixed']} / {i['no_reference']}",
                                                   f"{t['rows']:,}", f"{t['carry_over_learned']} / {t['carry_over_fixed']} / {t['no_reference']}"]) + " |")
    return "\n".join(out)


def dn_is_p14(r):
    return bool(r.get("rates"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True)
    ap.add_argument("--seed", default="42")
    ap.add_argument("--order", default="u9740,u9815,ups27,qe,astral,pxd030780,p14,sc")
    ap.add_argument("--legacy", action="store_true", help="print the legacy columns only (the board of earlier runs)")
    ap.add_argument("--detail", action="store_true", help="also print the per-level L1, A2, absence-class and blank tables")
    a = ap.parse_args(argv)
    R = load(a.results)
    builds = sorted({b for S in R.values() for r in S.values() for b in r["builds"]})
    order = a.order.split(",")
    order = order + [d for d in sorted(R) if d not in order]
    print(build_rows(R, a.seed, builds, order, a.legacy))
    if not a.legacy:
        for dn in order:
            if dn in R:
                S = R[dn]
                r = S[a.seed if a.seed in S else sorted(S)[0]]
                fam = r.get("families") or {}
                if fam.get("common_arms"):
                    print(f"\nCommon DE family of {r['dataset']}: arms {', '.join(label(x, builds) for x in fam['common_arms'])}; every other arm is scored on the part of it that it can test.")
    if a.detail:
        print(detail_tables(R, a.seed, builds, order))


if __name__ == "__main__":
    main()
