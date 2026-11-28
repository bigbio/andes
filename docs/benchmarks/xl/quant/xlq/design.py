"""Datasets of the cross-link quantification benchmark (datasets.json) and their truth."""
import importlib.util
import json
import os
import re
from pathlib import Path
from types import SimpleNamespace

from .keys import SEP, accessions

HERE = Path(__file__).resolve().parent.parent
XL = HERE.parent


def load(name, manifest=HERE / "datasets.json"):
    """The dataset `name` of the manifest, with the runs of the false-transfer designs resolved to
    the raw-file stems of docs/benchmarks/xl/datasets.json."""
    sets = {d["name"]: d for d in json.loads(Path(manifest).read_text(encoding="utf-8"))}
    if name not in sets:
        raise KeyError(f"unknown dataset {name!r}; known: {', '.join(sets)}")
    d = dict(sets[name])
    if d["design"] == "absence":
        xl = {x["name"]: x for x in json.loads((XL / "datasets.json").read_text(encoding="utf-8"))}
        d["runs"] = [dict(r, run=Path(xl[r["dataset"]]["raw"]).stem, raw=xl[r["dataset"]]["raw"],
                          url=xl[r["dataset"]]["url"]) for r in d["runs"]]
    return d


def run_key(name):
    """Run names compared without directory (`/` or `\\`), extension, case, `-` and `_`."""
    stem = re.split(r"[\\/]", str(name))[-1]
    for ext in (".raw", ".mzml", ".mgf", ".d"):
        if stem.lower().endswith(ext):
            stem = stem[: -len(ext)]
    return stem.lower().replace("-", "").replace("_", "")


def match_runs(df, dataset):
    """Rename the table's runs to the design's run names; rows of other runs are dropped.

    Returns (table, dropped run names). Two run names of the table that stand for one run of the
    design are an error: their rows would be duplicates of one (precursor, run)."""
    by_key = {run_key(r["run"]): r["run"] for r in dataset["runs"]}
    names = {}
    for name in sorted(set(df["run"])):
        target = by_key.get(run_key(name))
        if target is not None:
            names.setdefault(target, []).append(name)
    clash = {t: n for t, n in names.items() if len(n) > 1}
    if clash:
        target, both = sorted(clash.items())[0]
        raise ValueError(f"{len(clash)} run(s) of the design are named more than once in the table, "
                         f"e.g. {both} are all {target!r}")
    mapped = df["run"].map(lambda r: by_key.get(run_key(r)))
    dropped = sorted(set(df.loc[mapped.isna(), "run"]))
    return df[mapped.notna()].assign(run=mapped[mapped.notna()]), dropped


def _xl_eval(data_dir):
    os.environ["XLBENCH_DATA"] = str(data_dir)
    spec = importlib.util.spec_from_file_location("xl_eval", XL / "eval.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def absence(dataset, data_dir):
    """`absent(row)` for the transfer designs, from the truth of docs/benchmarks/xl/eval.py.

    A residue pair is absent from a run when that run's truth does not hold for it: for a
    peptide-library truth (`lib:`), its peptide pair is not a library pair linked at the
    designated residues; for the Clasen truth (`clasen:`), its proteins are not in one interaction
    group of the run's batch. A run whose truth is null is not scored.

    What a transferred row of an absent residue pair means depends on the design. Where the runs
    hold different real cross-links (Clasen batches), it is a false transfer: a real precursor
    placed in a sample that lacks it. Where every run has the same truth (the Matzinger library
    alone and in a HEK background that carries no cross-links), it is a wrong identity carried into
    another run.
    """
    ev = _xl_eval(data_dir)
    truths = {r["run"]: ev.truth_for(r["truth"]) for r in dataset["runs"] if r.get("truth")}

    def absent(row):
        tr = truths.get(row["run"])
        if tr is None:
            return None
        a, b = str(row["residue_pair"]).split(SEP)
        if isinstance(tr, ev.LibTruth):
            ap, bp = row.get("alpha_peptide"), row.get("beta_peptide")
            if not isinstance(ap, str) or not isinstance(bp, str):
                raise ValueError("a peptide-library truth needs alpha_peptide/beta_peptide columns")
            a_site, b_site = int(row["alpha_site"]), int(row["beta_site"])
            r = SimpleNamespace(alpha_peptide=ap, alpha_site=a_site, alpha_site_kind=ap[a_site - 1],
                                beta_peptide=bp, beta_site=b_site, beta_site_kind=bp[b_site - 1])
        else:
            r = SimpleNamespace(alpha_proteins=";".join(sorted(accessions(a))),
                                beta_proteins=";".join(sorted(accessions(b))))
        return not tr.is_true(r)

    return absent
