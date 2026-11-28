"""Run andes on one dataset of datasets.json and write its cross-link quantities.

    python3 -m xlq.andes --dataset hsa_dilution --data $XLBENCH_DATA --andes /path/to/andes --tag mytag

Three steps, each with its wall time, CPU time and peak resident memory in `run.json`:

1. `andes --xl` on all runs of the dataset in one process (the linker and its sites are the
   dataset's `search.xl`);
2. `andes cohort xl-fdr` on that search's table: the q-values of all runs together, after
   rescoring;
3. `andes --xl --lfq --mbr --xl-cohort-fdr`: the same search with quantification of the
   cross-links that step 2 accepts at the dataset's `search.xl_fdr`.

The outputs go to `$XLBENCH_DATA/results/<tag>/<dataset>/`; the feature table that
`xlq.evaluate` reads is `quant.xl.features.parquet` there. Thermo .raw input needs an andes
built with `--features thermo` and the .NET 8 runtime (DOTNET_ROOT).
"""
import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

from . import design


def run_measured(cmd, log):
    """Run cmd with its output in the file `log`; return its exit code, wall time, CPU time
    (user + system) and peak resident memory."""
    t = time.time()
    with open(log, "w") as out:
        if not hasattr(os, "wait4"):
            code = subprocess.run(cmd, stdout=out, stderr=subprocess.STDOUT).returncode
            return {"exit": code, "wall_s": round(time.time() - t, 1)}
        p = subprocess.Popen(cmd, stdout=out, stderr=subprocess.STDOUT)
        _, status, ru = os.wait4(p.pid, 0)
        p.returncode = os.waitstatus_to_exitcode(status)
    rss_bytes = ru.ru_maxrss if sys.platform == "darwin" else ru.ru_maxrss * 1024
    return {"exit": p.returncode, "wall_s": round(time.time() - t, 1),
            "cpu_s": round(ru.ru_utime + ru.ru_stime, 1), "max_rss_mb": round(rss_bytes / 2**20, 1)}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--data", required=True, help="the benchmark data directory ($XLBENCH_DATA)")
    ap.add_argument("--andes", required=True, help="the andes binary")
    ap.add_argument("--tag", required=True, help="name of this result under <data>/results/")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--manifest", default=str(design.HERE / "datasets.json"))
    ap.add_argument("extra", nargs="*", help="further arguments of the quantification step, after --")
    a = ap.parse_args(argv)

    dataset = design.load(a.dataset, a.manifest)
    search = dataset.get("search")
    if search is None:
        sys.exit(f"{a.dataset}: datasets.json gives it no `search` settings")
    data = Path(a.data)
    out = data / "results" / a.tag / a.dataset
    out.mkdir(parents=True, exist_ok=True)
    raws = [data / "raw" / (r.get("raw") or f"{r['run']}.raw") for r in dataset["runs"]]
    missing = [str(p) for p in raws if not p.exists()]
    if missing:
        sys.exit(f"{a.dataset}: {len(missing)} raw files are missing, e.g. {missing[0]}; run xlq.fetch")
    spectra = [arg for p in raws for arg in ("--spectrum", str(p))]
    common = [a.andes, *spectra, "--database", str(data / dataset["fasta"]["file"]),
              "--xl", search["xl"], "--threads", str(a.threads)]
    steps = [
        ("search", [*common, "--output-pin", str(out / "search.pin")]),
        ("fdr", [a.andes, "cohort", "xl-fdr", str(out / "search.xl.tsv"), "--output", str(out / "cohort")]),
        ("quant", [*common, "--lfq", "--mbr", "--xl-fdr", str(search["xl_fdr"]),
                   "--xl-cohort-fdr", str(out / "cohort"), "--output-pin", str(out / "quant.pin"), *a.extra]),
    ]
    version = subprocess.run([a.andes, "--version"], capture_output=True, text=True).stdout.strip()
    record = {"dataset": a.dataset, "andes": version, "threads": a.threads, "steps": {}}
    for name, cmd in steps:
        usage = run_measured(cmd, out / f"{name}.log")
        relative = [str(c).replace(str(data), "$XLBENCH_DATA") for c in cmd[1:]]
        record["steps"][name] = dict(usage, command="andes " + " ".join(relative))
        print(name, json.dumps(usage), file=sys.stderr, flush=True)
        (out / "run.json").write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
        if usage["exit"] != 0:
            sys.exit(f"{name} failed with exit code {usage['exit']}; see {out / (name + '.log')}")
    print(out / "quant.xl.features.parquet")


if __name__ == "__main__":
    main()
