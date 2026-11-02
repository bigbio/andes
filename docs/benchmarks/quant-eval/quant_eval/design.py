"""Dataset design from two small TSV files.

design.tsv   columns: run, condition (optional), blank (optional, true for blank runs), absent_species (optional: a species
             that is absent from that run, e.g. yeast in human-only runs), order (optional: integer acquisition order; used to
             find the run acquired immediately before each blank for the carry-over count). Run names are the andes run stems.
amounts.tsv  columns: species, condition, amount. `condition` may be `*` (every condition); `amount` may be `=condition`
             (the amount equals the numeric condition value, e.g. UPS1 amol per level). Species without an amount are not
             scored in DE or accuracy; an amount of 0 means absent at that condition.
The expected log2 fold change of a species between two conditions is log2(amount_high) - log2(amount_low)."""
import math

import pandas as pd


def _truthy(x):
    return str(x).strip().lower() in ("true", "1", "yes", "t")


def _order(d, path):
    """run -> integer acquisition order (empty cells skipped); a non-integer or repeated order is refused."""
    o = {}
    for r, x in zip(d.run, d["order"]):
        if x.strip() == "":
            continue
        try:
            o[r] = int(x)
        except ValueError:
            raise SystemExit(f"{path}: order {x!r} of run {r!r} is not an integer")
    if len(set(o.values())) != len(o):
        raise SystemExit(f"{path}: repeated values in column 'order'")
    return o


class Design:
    def __init__(self, design_tsv, amounts_tsv=None, null=None, classes=(), low_max=500.0):
        d = pd.read_csv(design_tsv, sep="\t", dtype=str).fillna("")
        if "run" not in d.columns:
            raise SystemExit(f"{design_tsv}: column 'run' required")
        self.runs = sorted(d.run)
        self.cond = {r: c for r, c in zip(d.run, d["condition"])} if "condition" in d.columns else {}
        self.cond = {r: c for r, c in self.cond.items() if c != ""}
        self.blank = {r for r, b in zip(d.run, d["blank"]) if _truthy(b)} if "blank" in d.columns else set()
        self.absent = {r: s for r, s in zip(d.run, d["absent_species"]) if s != ""} if "absent_species" in d.columns else {}
        self.order = _order(d, design_tsv) if "order" in d.columns else {}
        self.null = null
        self.classes = tuple(classes)
        self.low_max = low_max
        self.levels = []
        if self.cond:
            lv = set(self.cond.values())
            try:
                self.levels = sorted(lv, key=float)
            except ValueError:
                self.levels = sorted(lv)
        self._amount = {}
        if amounts_tsv:
            a = pd.read_csv(amounts_tsv, sep="\t", dtype=str)
            for s, c, v in zip(a.species, a.condition, a.amount):
                for cc in (self.levels if c == "*" else [c]):
                    self._amount[(s, cc)] = float(cc) if v == "=condition" else float(v)
        self.has_de = bool(self._amount) and len(self.levels) >= 2 and null is not None

    def run_before(self, r):
        """The run acquired immediately before run r (largest order below r's order); None without an order or a predecessor."""
        if r not in self.order:
            return None
        prev = [x for x, o in self.order.items() if o < self.order[r]]
        return max(prev, key=self.order.get) if prev else None

    def amount(self, s, c):
        v = self._amount.get((s, c))
        return v if (v is not None and v > 0) else None

    def expected(self, s, c):
        """Expected log2 amount of species s at condition c; None when not defined."""
        v = self.amount(s, c)
        return math.log2(v) if v is not None else None

    def expected_fc(self, s, lo, hi):
        a, b = self.amount(s, lo), self.amount(s, hi)
        return math.log2(b) - math.log2(a) if (a is not None and b is not None) else float("nan")

    def is_low(self, level):
        try:
            return float(level) <= self.low_max
        except ValueError:
            return False
