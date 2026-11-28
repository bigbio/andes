"""Download the inputs of the cross-link quantification benchmark into a data directory.

    python3 -m xlq.fetch $XLBENCH_DATA [--only hsa_dilution,hsa_replicates] [--references-only]

Raw files go to raw/, FASTAs to db/ and the published reference outputs to ref/, next to the files
of docs/benchmarks/xl/fetch.py (the false-transfer designs reuse its raws, FASTAs and truth tables;
run that script first for them). Large files are fetched in parallel ranges by the same code.
"""
import argparse
import importlib.util
import json
import subprocess
from pathlib import Path

from .design import HERE, XL, load


def _xl_fetch():
    spec = importlib.util.spec_from_file_location("xl_fetch", XL / "fetch.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def small(url, dest):
    if dest.exists():
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".tmp")
    subprocess.run(["curl", "-sS", "-fL", "--retry", "4", "-o", str(tmp), url], check=True)
    tmp.rename(dest)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("data")
    ap.add_argument("--only", help="comma-separated dataset names")
    ap.add_argument("--references-only", action="store_true", help="FASTAs and published outputs, no raw files")
    a = ap.parse_args(argv)
    data = Path(a.data)
    names = [d["name"] for d in json.loads((HERE / "datasets.json").read_text(encoding="utf-8"))]
    if a.only:
        names = [n for n in names if n in a.only.split(",")]
    xl_fetch = None if a.references_only else _xl_fetch()
    for name in names:
        d = load(name)
        if d["fasta"].get("url"):
            small(d["fasta"]["url"], data / d["fasta"]["file"])
        for entry in d.get("published", {}).values():
            small(entry["url"], data / entry["file"])
        if a.references_only:
            continue
        (data / "raw").mkdir(parents=True, exist_ok=True)
        for r in d["runs"]:
            raw = r.get("raw") or Path(r["url"]).name
            xl_fetch.fetch_large(r["url"], str(data / "raw" / raw))


if __name__ == "__main__":
    main()
