"""Summarise the published Skyline CV reports of a replicate dataset (they hold CVs, not areas).

    python3 -m xlq.published --dataset hsa_replicates --data $XLBENCH_DATA --out results/hsa_replicates.published.json

For each `published` entry named `skyline_cv_<group>`: the number of precursors and the median CV
per Skyline transition (`precursor` is the monoisotopic peak, then `[M+1]`, `[M+2]`). The areas
behind these CVs are not published, so this is the only reference for that dataset until the
reference workflow is run again.
"""
import argparse
import json
import sys
from pathlib import Path

from . import design, skyline
from .evaluate import sha256

PREFIX = "skyline_cv_"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--data", required=True)
    ap.add_argument("--manifest", default=str(design.HERE / "datasets.json"))
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    dataset = design.load(a.dataset, a.manifest)
    out = {"dataset": a.dataset, "groups": {}}
    for key, entry in dataset.get("published", {}).items():
        if not key.startswith(PREFIX):
            continue
        path = Path(a.data) / entry["file"]
        cv = skyline.read_cv_report(path)
        out["groups"][key[len(PREFIX):]] = {
            "file": path.name,
            "sha256": sha256(path),
            "precursors": int(cv["precursor"].nunique()),
            "median_cv": {ion: round(float(v), 4) for ion, v in cv.groupby("ion")["cv"].median().items()},
        }
    if not out["groups"]:
        sys.exit(f"{a.dataset}: no {PREFIX}* entries in its `published` block")
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
