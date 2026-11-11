"""Low-level UPS1 detection read three ways: with its levels named, with a precision companion, and at a matched number of rows.

The L1 count of the board (`design_detection`) is the number of UPS1 precursors with one row within 2x of the design at the three
lowest levels. It grows with the number of rows an arm reports (a wrong row has a chance of landing within 2x), so every level
here also carries what the arm reports and how accurate it is:

  rows               scored UPS1 rows of the arm at the level (precursors with a reference). `unscored_rows` are the UPS1 rows of precursors without a
                     reference row (seen only at the lowest levels): they leave the count and its denominator, so the board also prints the share
                     within 2x with them counted as not within 2x (`within2x_incl_unscored` of `l1_pooled`)
  keys_within2x      distinct precursors with a row within 2x (the L1 count of the level)
  share_within2x     rows within 2x / rows
  median_abs_err     median |log2 error| of the rows;  median_err the signed median (positive = over-quantified)
  transfer, identified   the same statistics on the transferred rows and on the identified rows alone (A2: error of transfers)
  l1_pooled          the same statistics over the rows of all L1 levels together
  key_correct_rows   (L1 levels) {c: number of precursors with exactly c rows within 2x}, for the order-free reading
  key_ranks, row_ranks   (L1 levels) 1-based ranks, in the arm's own quality order, at which a new precursor / a row within 2x appears

Role of a level: `L1` (one of the `n_l1` lowest), `ref` (one of the levels the per-precursor references are built from, those above
the `ref_skip` lowest) or `mid`. A level that is both is reported in `self_reference`; with the harness defaults (3 and 3) none is.

Order of rows for the matched reading: an arm that reports a quality column (andes: identified rows first, then transfers by quant_q;
FlashLFQ: PIP q-value) is read in that order, best first; an arm without one (MaxQuant, FragPipe) by offset-corrected intensity,
brightest first. `matched_reading` limits every arm to the same number of rows at a level (the smallest row count among the reference arms)
and counts what lies within 2x among those rows, in that order (`keys`) and for rows drawn at random (`keys_random`, expected value; an
order-free reading, because an intensity order puts the brightest, most over-quantified rows first at the lowest levels)."""
import bisect

import numpy as np

from .metrics import ups_rows

WITHIN = 1.0   # log2: within 2x


def has_quality_order(raw):
    """True when the arm table carries a usable quality column (more than one distinct finite value among its transfers)."""
    if "thr" not in raw.columns:
        return False
    t = raw.loc[raw.is_transfer.astype(bool), "thr"]
    t = t[np.isfinite(t)]
    return t.nunique() > 1


def _stats(e):
    """Row statistics of an error vector (log2)."""
    n = int(len(e))
    if n == 0:
        return {"rows": 0, "within2x": 0, "share_within2x": None, "median_err": None, "median_abs_err": None}
    w = int((np.abs(e) <= WITHIN).sum())
    return {"rows": n, "within2x": w, "share_within2x": w / n, "median_err": float(np.median(e)), "median_abs_err": float(np.median(np.abs(e)))}


def order_rows(u, quality):
    """Rows of one level, best first: identified before transferred, then quant_q ascending (NaN last); intensity (descending) when the arm has no quality column."""
    if quality:
        thr = u.thr.where(np.isfinite(u.thr), np.inf)
        return u.assign(_t=u.is_transfer.astype(int), _q=thr).sort_values(["_t", "_q", "n"], ascending=[True, True, False], kind="stable")
    return u.sort_values("n", ascending=False, kind="stable")


def ranks(u_ordered):
    """(ranks at which a new precursor within 2x appears, ranks of the rows within 2x), 1-based, over rows in the given order."""
    ok = (u_ordered.err.abs() <= WITHIN).values
    keys = u_ordered.key.values
    seen, kr, rr = set(), [], []
    for i in range(len(ok)):
        if ok[i]:
            rr.append(i + 1)
            if keys[i] not in seen:
                seen.add(keys[i])
                kr.append(i + 1)
    return kr, rr


def l1_detail(raw, design, sp, ref_skip=3, n_l1=3):
    """Per-level detail of one arm table (see the module docstring). {} when the design has no UPS1 levels."""
    r = ups_rows(raw, design, sp, ref_skip)
    if r is None:
        return {}
    U, unscored = r
    levels = [str(x) for x in design.levels]
    l1, ref = levels[:n_l1], levels[ref_skip:]
    quality = has_quality_order(raw)
    out = {"levels": levels, "l1_levels": l1, "ref_levels": ref, "ref_skip": ref_skip, "self_reference": [x for x in l1 if x in ref],
           "order": "quality" if quality else "intensity", "per_level": {}}
    out["l1_pooled"] = _stats(U[U.lev.isin([float(x) for x in l1])].err.values)
    out["l1_pooled"]["unscored_rows"] = int(sum(unscored.get(L, 0) for L in l1))
    out["l1_pooled"]["within2x_incl_unscored"] = out["l1_pooled"]["within2x"] / max(out["l1_pooled"]["rows"] + out["l1_pooled"]["unscored_rows"], 1)
    for L in levels:
        u = U[U.lev == float(L)]
        ok = u.err.abs() <= WITHIN
        o = {"role": "L1" if L in l1 else ("ref" if L in ref else "mid"), "unscored_rows": unscored.get(L, 0),
             "keys": int(u.key.nunique()), "keys_within2x": int(u[ok].key.nunique())}
        o.update(_stats(u.err.values))
        o["transfer"] = _stats(u[u.is_transfer.astype(bool)].err.values)
        o["identified"] = _stats(u[~u.is_transfer.astype(bool)].err.values)
        if L in l1:
            o["key_ranks"], o["row_ranks"] = ranks(order_rows(u, quality))
            cc = u[ok].groupby("key").size().value_counts()
            o["key_correct_rows"] = {str(int(c)): int(v) for c, v in cc.items()}
        out["per_level"][L] = o
    return out


def counts_at(detail_level, k):
    """Precursors and rows within 2x among the first `k` rows of an arm at one level (the arm's own order), and the rows considered."""
    n = min(k, detail_level["rows"])
    return {"keys": bisect.bisect_right(detail_level["key_ranks"], k), "rows_within2x": bisect.bisect_right(detail_level["row_ranks"], k),
            "rows": n}


def expected_keys_random(detail_level, k):
    """Expected number of precursors with a row within 2x among `k` rows drawn at random (without replacement) from the arm's `rows` at one level:
    a precursor with c correct rows is counted with probability 1 - C(n - c, k) / C(n, k) = 1 - prod_{i<c} (n - k - i) / (n - i). An order-free reading of the count at a row budget."""
    n = detail_level["rows"]
    k = min(k, n)
    tot = 0.0
    for c, cnt in detail_level["key_correct_rows"].items():
        c = int(c)
        p0 = 0.0 if n - k < c else float(np.prod([(n - k - i) / (n - i) for i in range(c)]))
        tot += cnt * (1.0 - p0)
    return tot


def matched_reading(details, ref_arms, level):
    """Every arm of `details` ({arm: l1_detail}) limited to the same number of rows at `level`: k = the smallest row count among
    `ref_arms` (arms that have the level). Returns (k, {arm: counts_at}) or (0, {})."""
    have = [a for a in ref_arms if details.get(a, {}).get("per_level", {}).get(level, {}).get("rows")]
    if not have:
        return 0, {}
    k = min(details[a]["per_level"][level]["rows"] for a in have)
    out = {}
    for a, d in details.items():
        pl = d.get("per_level", {}).get(level)
        if pl and "key_ranks" in pl:
            c = counts_at(pl, k)
            c["share_within2x"] = c["rows_within2x"] / c["rows"] if c["rows"] else None
            c["keys_random"] = expected_keys_random(pl, k)
            out[a] = c
    return k, out
