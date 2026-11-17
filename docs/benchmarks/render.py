#!/usr/bin/env python3
"""Render the headline benchmark tables from docs/benchmarks/results/*.json.

Each table is stored once, in results/. This script writes it into every page that shows it,
between the markers `<!-- BEGIN generated: <id> -->` and `<!-- END generated: <id> -->`, and
writes website/benchmarks.json, the input of the site figures. With --check it writes nothing
and exits non-zero when any copy differs from what the JSON gives.
"""
import argparse
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
RESULTS = HERE / "results"
SITE_JSON = ROOT / "website" / "benchmarks.json"
PAGES = [ROOT / "README.md", ROOT / "website" / "pages" / "quantification.md", *sorted(HERE.glob("*.md"))]
BLOCK = re.compile(
    r"(<!-- BEGIN generated: (?P<id>[\w-]+) \(docs/benchmarks/render\.py\) -->\n)"
    r"(?P<body>.*?)"
    r"(<!-- END generated: (?P=id) -->)",
    re.S,
)


def load():
    return {p.stem: json.loads(p.read_text(encoding="utf-8")) for p in sorted(RESULTS.glob("*.json"))}


def fmt(v, spec=""):
    if v is None:
        return "—"
    if isinstance(v, list):
        return "–".join(fmt(x, spec) for x in v)
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, int):
        return f"{v:,}"
    if isinstance(v, float):
        return format(v, spec) if spec else f"{v:g}"
    return str(v)


def lookup(obj, path):
    for key in path.split("."):
        if isinstance(obj, list):
            obj = next(x for x in obj if x.get("id") == key)
        elif "tools" in obj and key in {t["id"] for t in obj["tools"]}:
            obj = next(t for t in obj["tools"] if t["id"] == key)["values"]
        elif "modules" in obj and key in obj["modules"]:
            obj = obj["modules"][key]["rows"]
        else:
            obj = obj[key]
    return obj


def resolve(text, data):
    sources = {"engines": data["identification"]["engines"], "pb": data["proteobench"]}

    def repl(m):
        src, path, spec = m.group(1), m.group(2), m.group(3) or ""
        return fmt(lookup(sources[src], path), spec.lstrip("|"))

    return re.sub(r"\{\{(\w+):([^|}]+)(\|[^}]+)?\}\}", repl, text)


def table(header, align, rows):
    sep = ["---:" if a == "right" else "---" for a in align]
    def line(cells):
        return "|" + "".join(f" {c} |" if c else " |" for c in cells)

    lines = [line(header), "|" + "|".join(sep) + "|"] + [line(r) for r in rows]
    return "\n".join(lines) + "\n"


def summary(data, name):
    s = data[name]["summary"]
    return table(s["header"], s["align"], [[resolve(c, data) for c in r] for r in s["rows"]])


def tool_cell(tool, ds, key, unit="", spec=""):
    v = tool["values"].get(ds["id"], {}).get(key)
    if v is None:
        return "—"
    mark = " ‡" if tool["values"][ds["id"]].get("build", tool.get("build")) != tool.get("build") else ""
    return f"{fmt(v, spec)}{unit}{mark}"


def engines_psms(data):
    e = data["identification"]["engines"]
    header = ["Engine"] + [d["long"] for d in e["datasets"] if d["id"] != "phospho"]
    sets = [d for d in e["datasets"] if d["id"] != "phospho"]
    rows = []
    for t in e["tools"]:
        label = t["name"]
        if t.get("build"):
            label = f"{t['name']} `{t['build']}` ({t['date']})"
        elif t.get("date"):
            label = f"{t['name']} ({t['date']})"
        else:
            label = f"{t['name']} †"
        rows.append([label] + [tool_cell(t, d, "psms") for d in sets])
    for t in e["tools"]:
        if any("wall" in t["values"].get(d["id"], {}) for d in sets):
            rows.append([f"*{t['name'].split()[0]} wall time*"] + [f"*{tool_cell(t, d, 'wall', ' s')}*" for d in sets])
    return table(header, ["left"] + ["right"] * len(sets), rows)


def engines_cost(data):
    e = data["identification"]["engines"]
    rows = []
    for d in e["datasets"]:
        first = True
        for t in e["tools"]:
            if "wall" not in t["values"].get(d["id"], {}):
                continue
            rows.append([
                d["name"] if first else "",
                t["name"],
                tool_cell(t, d, "psms"),
                tool_cell(t, d, "wall", " s"),
                tool_cell(t, d, "cpu", " s"),
                tool_cell(t, d, "mem", " GB", ".1f"),
            ])
            first = False
    header = ["dataset", "engine", "PSMs @ q≤0.01", "wall", "CPU time", "peak memory"]
    return table(header, ["left", "left", "right", "right", "right", "right"], rows)


def opt_in(data):
    rows = []
    for r in data["opt_in"]["rows"]:
        cells = [r["mode"]]
        for side in ("andes", "without"):
            v = r[side]
            cells.append(v if isinstance(v, str) else f"{fmt(v['psms'])} / {fmt(v['peptides'])}")
        rows.append(cells)
    return table(["mode", "andes", "without the mode"], ["left", "left", "left"], rows)


def leaderboard(data, module):
    rows = [[r["submission"], r["mbr"], fmt(r["precursors"]), fmt(r["eps"], ".3f"), fmt(r["ecoli_eps"], ".3f"),
             fmt(r["cv"], ".3f"), r["transfer_control"]] for r in data["proteobench"]["modules"][module]["rows"]]
    header = ["tool (submission)", "MBR", "precursors", "median \\|ε\\|", "E. coli \\|ε\\|", "CV", "transfer control"]
    return table(header, ["left", "left", "right", "right", "right", "right", "left"], rows)


def proteobench_plotted(data):
    rows = []
    for module, m in data["proteobench"]["modules"].items():
        first = True
        for r in m["rows"]:
            if "plot" not in r:
                continue
            name = f"**{r['plot']}**" if r.get("andes") else r["plot"]
            rows.append([f"{module} ({m['dataset']})" if first else "", name, fmt(r["precursors"]),
                         fmt(r["eps"], ".3f"), fmt(r["cv"], ".3f")])
            first = False
    header = ["module", "submission", "precursors in ≥ 3 runs", "median \\|ε\\|", "CV"]
    return table(header, ["left", "left", "right", "right", "right"], rows)


TABLES = {
    "identification-summary": lambda d: summary(d, "identification"),
    "quantification-summary": lambda d: summary(d, "quantification"),
    "opt-in": opt_in,
    "engines-psms": engines_psms,
    "engines-cost": engines_cost,
    "proteobench-qexactive": lambda d: leaderboard(d, "QExactive"),
    "proteobench-astral": lambda d: leaderboard(d, "Astral"),
    "proteobench-plotted": proteobench_plotted,
}


def mid(v):
    if not isinstance(v, list):
        return v
    m = (v[0] + v[1]) / 2
    return int(m) if m == int(m) else m


def site_json(data):
    e = data["identification"]["engines"]
    sets = e["datasets"]
    andes = next(t for t in e["tools"] if t["id"] == "andes")
    out = {
        "_provenance": "Generated by docs/benchmarks/render.py from docs/benchmarks/results/*.json; do not edit. "
                       + e["_note"],
        "andes_commit": andes["build"],
        "andes_date": andes["date"],
        "datasets": [d["chart"] for d in sets],
        "psms": {t["name"]: [mid(t["values"].get(d["id"], {}).get("psms")) for d in sets] for t in e["tools"]},
    }
    for key, name in (("wall", "wall_s"), ("cpu", "cpu_s"), ("mem", "peak_gb")):
        cost = {}
        for t in e["tools"]:
            vals = []
            for d in sets:
                v = t["values"].get(d["id"], {})
                same = v.get("host", t["host"]) == "vm"
                vals.append(mid(v[key]) if key in v and same else None)
            if any(x is not None for x in vals):
                cost[t["name"]] = vals
        out[name] = cost
    out["wall_note"] = ("Bars: midpoint of a reported range; cost only where both engines ran on the benchmark VM "
                        "(andes UPS1 was measured on another host and is left out).")
    out.update({k: v for k, v in data["figures"].items() if not k.startswith("_")})
    rows = data["opt_in"]["rows"]
    chim = [r for r in rows if "chart" in r]
    out["opt_in"] = {
        "_note": data["opt_in"]["_provenance"],
        "labels": ["default"] + [r["chart"] for r in chim],
        "psms_lo": [chim[0]["without"]["psms"][0]] + [r["andes"]["psms"][0] for r in chim],
        "psms_hi": [chim[0]["without"]["psms"][1]] + [r["andes"]["psms"][1] for r in chim],
    }
    out["quant_accuracy"] = {"_note": data["proteobench"]["_provenance"]}
    for module, m in data["proteobench"]["modules"].items():
        out["quant_accuracy"][module] = [[r["plot"], r["precursors"], r["eps"], 1 if r.get("andes") else 0]
                                         for r in m["rows"] if "plot" in r]
    return json.dumps(out, ensure_ascii=False, indent=1) + "\n"


def render_page(text, data, seen):
    def repl(m):
        tid = m.group("id")
        if tid not in TABLES:
            raise SystemExit(f"unknown generated table id '{tid}'")
        seen.add(tid)
        return m.group(1) + TABLES[tid](data) + m.group(4)

    return BLOCK.sub(repl, text)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="write nothing; exit 1 if any copy differs from the JSON")
    args = ap.parse_args()
    data = load()
    stale, seen = [], set()
    outputs = {p: render_page(p.read_text(encoding="utf-8"), data, seen) for p in PAGES}
    outputs[SITE_JSON] = site_json(data)
    missing = set(TABLES) - seen
    if missing:
        raise SystemExit(f"tables with no markers in any page: {', '.join(sorted(missing))}")
    for path, new in outputs.items():
        old = path.read_text(encoding="utf-8") if path.exists() else ""
        if old != new:
            stale.append(path.relative_to(ROOT).as_posix())
            if not args.check:
                path.write_text(new, encoding="utf-8")
    if args.check:
        if stale:
            print("out of date with docs/benchmarks/results (run docs/benchmarks/render.py):", *stale, sep="\n  ")
            sys.exit(1)
        print(f"{len(outputs)} files match docs/benchmarks/results")
    else:
        print("rewrote:", *stale, sep="\n  ") if stale else print("nothing to rewrite")


if __name__ == "__main__":
    main()
