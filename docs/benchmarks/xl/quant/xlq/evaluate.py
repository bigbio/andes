"""Score cross-link quantification tables on one dataset of datasets.json.

    python3 -m xlq.evaluate --dataset hsa_dilution --data $XLBENCH_DATA \\
        --arm skyline_published:skyline:ref/HSA_BS3_DDA_dilution.sky \\
        --arm andes:table:/path/to/xl_quant.tsv --out results/hsa_dilution.json

An arm is NAME:KIND:PATH. KIND `table` is the long table of xlq/table.py (TSV or Parquet); KIND
`skyline` is a Skyline document with pseudo-peptide cross-links, read with the modification names
of the dataset's `published.skyline_document` entry. A relative PATH is taken under --data. Every
arm is scored on its own residue pairs and, with more than one arm, also on the residue pairs every
arm quantifies (`common`). The `own` scores of an arm depend on that arm alone: adding or removing
another arm does not change them. Only the choice of the `common` residue pairs looks across arms.
There, keys that differ only in the protein positions listed for a peptide are matched by
`keys.comparison_ids`; the keys themselves are never changed, so the truth of an arm reads the
positions and proteins that arm reported. The result says how many sets of keys were matched and
how many were ambiguous.

An arm that keeps no row (its run names do not match the design, or a filter removes everything)
is an error. Runs of the table that the design does not have, and runs of the design without a row
in an arm, are named on stderr and in the result.

KIND takes row filters after `?`, joined by `&`: `evidence=identified` keeps the identified rows
(the arm without transfers); `min_COL=X` and `max_COL=X` keep the rows whose column COL is at least
or at most X. A filter on COL applies only to rows with a value in COL, so `max_quant_q=0.01` keeps
identified rows that carry no quant_q. Example: `skyline?min_isotope_dotp=0.9`.
"""
import argparse
import hashlib
import json
import sys
import warnings
from pathlib import Path

import pandas as pd

from . import design, keys, metrics, skyline, table
from .keys import read_fasta


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, result):
    """Write a result file. NaN is refused: it is not JSON, and a metric without a value is null."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(result, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    print(f"wrote {path}", file=sys.stderr)


def parse_kind(kind):
    """`skyline?min_isotope_dotp=0.9&evidence=identified` -> ('skyline', {filters})."""
    base, _, opts = kind.partition("?")
    filters = {}
    for opt in filter(None, opts.split("&")):
        key, sep, value = opt.partition("=")
        if not sep:
            raise ValueError(f"filter {opt!r} needs KEY=VALUE")
        filters[key] = value
    return base, filters


def apply_filters(df, filters):
    for key, value in filters.items():
        if key == "evidence":
            if value not in table.EVIDENCE:
                raise ValueError(f"filter evidence={value}: not one of {', '.join(table.EVIDENCE)}")
            df = df[df["evidence"] == value]
        elif key[:4] in ("min_", "max_"):
            col = key[4:]
            if col not in df.columns:
                raise ValueError(f"filter {key}: the table has no column {col!r}")
            x = pd.to_numeric(df[col], errors="coerce")
            ok = x >= float(value) if key.startswith("min_") else x <= float(value)
            df = df[x.isna() | ok]
        else:
            raise ValueError(f"unknown filter {key!r} (evidence, min_COL, max_COL)")
    return df


def read_arm(kind, path, dataset, data_dir):
    base, filters = parse_kind(kind)
    if base == "table":
        df, counts = table.read(path), {}
    elif base == "skyline":
        doc = dataset.get("published", {}).get("skyline_document")
        if doc is None:
            raise ValueError(f"dataset {dataset['name']} has no published.skyline_document entry, "
                             "which names the modifications of a Skyline arm")
        proteins = read_fasta(Path(data_dir) / dataset["fasta"]["file"])
        df, counts = skyline.read_document(path, proteins, doc["spacer_mod"], doc["site_mod"])
    else:
        raise ValueError(f"unknown arm kind {base!r} (table, skyline)")
    kept = apply_filters(df, filters)
    if filters:
        counts = dict(counts, rows_before_filters=int(len(df)), rows_after_filters=int(len(kept)))
    return kept, counts


def score(df, dataset, data_dir, rollup, keep=None):
    """All metrics the dataset's design supports, for one arm's table."""
    runs = [r["run"] for r in dataset["runs"]]
    if keep is not None:
        df = df[df["residue_pair"].isin(keep)]
    M = metrics.log_matrix(df, runs, rollup)
    P = metrics.precursor_matrix(df, runs)
    if dataset.get("normalise", False):
        shift = metrics.run_shifts(M, runs)
        M, P = M - shift, P - shift
    out = {"depth": metrics.depth(df, M)}
    if dataset["design"] == "dilution":
        level_of = {r["run"]: r["level"] for r in dataset["runs"]}
        out["dilution"] = metrics.dilution(M, level_of, dataset["reference_level"])
        groups = {f"{lvl:g}": [r for r in runs if level_of[r] == lvl] for lvl in sorted(set(level_of.values()))}
        out["precision"] = metrics.precision(M, groups)
        out["precision_precursor"] = metrics.precision(P, groups)
        out["run_offsets"] = metrics.run_offsets(P, groups)
    elif dataset["design"] == "replicates":
        groups = {}
        for r in dataset["runs"]:
            groups.setdefault(r["group"], []).append(r["run"])
        out["precision"] = metrics.precision(M, groups)
        out["precision_precursor"] = metrics.precision(P, groups)
        out["run_offsets"] = metrics.run_offsets(P, groups)
        out["null_calls"] = metrics.null_calls(M, groups)
    elif dataset["design"] == "absence":
        out["transfers"] = metrics.transfers(df, design.absence(dataset, data_dir))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--data", required=True, help="the benchmark data directory ($XLBENCH_DATA)")
    ap.add_argument("--arm", action="append", required=True, help="NAME:KIND:PATH")
    ap.add_argument("--rollup", default="median-polish", choices=["median-polish", "sum"])
    ap.add_argument("--manifest", default=str(design.HERE / "datasets.json"))
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    dataset = design.load(a.dataset, a.manifest)
    design_runs = [r["run"] for r in dataset["runs"]]
    tables, result = {}, {"dataset": a.dataset, "design": dataset["design"], "rollup": a.rollup, "arms": {}}
    for spec in a.arm:
        if spec.count(":") < 2:
            sys.exit(f"--arm {spec}: expected NAME:KIND:PATH")
        name, kind, path = spec.split(":", 2)
        path = Path(path) if Path(path).is_absolute() else Path(a.data) / path
        try:
            df, counts = read_arm(kind, path, dataset, a.data)
            read_rows = len(df)
            df, dropped = design.match_runs(df, dataset)
        except ValueError as e:
            sys.exit(f"arm {name}: {e}")
        if df.empty:
            why = (f"none of its runs is a run of {a.dataset} (e.g. {dropped[0]!r} against {design_runs[0]!r})"
                   if read_rows else "its filters keep no row")
            sys.exit(f"arm {name}: no row to score, {why}")
        missing = [r for r in design_runs if r not in set(df["run"])]
        if dropped:
            print(f"arm {name}: {len(dropped)} run(s) of the table are not in the design and were "
                  f"dropped, e.g. {dropped[0]!r}", file=sys.stderr)
        if missing:
            print(f"arm {name}: no row for {len(missing)} of the {len(design_runs)} runs of the "
                  f"design, e.g. {missing[0]!r}", file=sys.stderr)
        tables[name] = df
        result["arms"][name] = {"kind": kind, "file": path.name, "sha256": sha256(path),
                                "read": counts, "runs_dropped": dropped, "runs_without_rows": missing}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        try:
            for name, df in tables.items():
                result["arms"][name]["own"] = score(df, dataset, a.data, a.rollup)
            if len(tables) > 1:
                arm_keys = {name: set(df["residue_pair"]) for name, df in tables.items()}
                ids, report = keys.comparison_ids(arm_keys)
                result["key_matching"] = report
                if report["ambiguous_sets"]:
                    print(f"{report['ambiguous_sets']} set(s) of residue-pair keys overlap between arms "
                          "but do not match one to one; they are not treated as shared", file=sys.stderr)
                common = set.intersection(*({ids[k] for k in arm} for arm in arm_keys.values()))
                result["common_residue_pairs"] = len(common)
                if not common:
                    print("the arms share no residue pair: the `common` view is empty", file=sys.stderr)
                for name, df in tables.items():
                    keep = {k for k in arm_keys[name] if ids[k] in common}
                    result["arms"][name]["common"] = score(df, dataset, a.data, a.rollup, keep=keep)
        except ValueError as e:
            sys.exit(f"{a.dataset}: {e}")
    write_json(a.out, result)


if __name__ == "__main__":
    main()
