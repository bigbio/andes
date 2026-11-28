"""Summarise the published reference outputs of a dataset that the evaluator does not read.

    python3 -m xlq.published --dataset hsa_dilution --data $XLBENCH_DATA --out hsa_dilution_published.json
    python3 -m xlq.published --dataset hsa_replicates --data $XLBENCH_DATA --out hsa_replicates_published.json

For each `published` entry of the dataset in datasets.json:

- `skyline_cv_<group>`, a Skyline report of CVs: the number of precursors and the median CV per
  Skyline transition (`precursor` is the monoisotopic peak, then `[M+1]`, `[M+2]`). Where only
  these reports are published, they are the one reference until the workflow is run again.
- `skyline_document` with `xifdr_psms`: whether the residue pairs read from the Skyline document
  are those of the target-target matches of the xiFDR PSM file (a check of the reader's mapping
  of pseudo-peptides to protein positions), and the document's site isomers
  (`skyline.site_isomers`).
"""
import argparse
import sys
from pathlib import Path

from . import design, skyline
from .evaluate import sha256, write_json
from .keys import read_fasta

PREFIX = "skyline_cv_"
EXAMPLES = 5


def cv_reports(dataset, data):
    out = {}
    for key, entry in dataset.get("published", {}).items():
        if not key.startswith(PREFIX):
            continue
        path = Path(data) / entry["file"]
        cv = skyline.read_cv_report(path)
        out[key[len(PREFIX):]] = {
            "file": path.name,
            "sha256": sha256(path),
            "precursors": int(cv["precursor"].nunique()),
            "median_cv": {ion: round(float(v), 4) for ion, v in cv.groupby("ion")["cv"].median().items()},
        }
    return out


def document(dataset, data):
    published = dataset.get("published", {})
    doc = published.get("skyline_document")
    if doc is None:
        return None
    path = Path(data) / doc["file"]
    proteins = read_fasta(Path(data) / dataset["fasta"]["file"])
    df, counts = skyline.read_document(path, proteins, doc["spacer_mod"], doc["site_mod"])
    out = {"file": path.name, "sha256": sha256(path), "read": counts,
           "site_isomers": skyline.site_isomers(path, doc["spacer_mod"], doc["site_mod"])}
    psms = published.get("xifdr_psms")
    if psms is not None:
        psm_path = Path(data) / psms["file"]
        ours, theirs = set(df["residue_pair"]), skyline.xifdr_residue_pairs(psm_path)
        out["xifdr_psms"] = {
            "file": psm_path.name,
            "sha256": sha256(psm_path),
            "residue_pairs_document": len(ours),
            "residue_pairs_xifdr": len(theirs),
            "shared": len(ours & theirs),
            "only_document": sorted(ours - theirs)[:EXAMPLES],
            "only_xifdr": sorted(theirs - ours)[:EXAMPLES],
        }
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--data", required=True)
    ap.add_argument("--manifest", default=str(design.HERE / "datasets.json"))
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    dataset = design.load(a.dataset, a.manifest)
    out = {"dataset": a.dataset}
    reports = cv_reports(dataset, a.data)
    if reports:
        out["cv_reports"] = reports
    doc = document(dataset, a.data)
    if doc is not None:
        out["skyline_document"] = doc
    if len(out) == 1:
        sys.exit(f"{a.dataset}: its `published` block has no {PREFIX}* and no skyline_document entry")
    write_json(a.out, out)


if __name__ == "__main__":
    main()
