#!/usr/bin/env python3
"""Run the CI benchmark datasets and append one record each to a JSON-lines results file.

    run_ci.py --data <data-dir> --andes <binary> --out results.jsonl [dataset ...]

<data-dir> holds what ../reproduce/fetch_spectra.sh and build_databases.sh write. Datasets are
the keys of datasets.json (default: all). Each dataset gives one line: commit, IDs (from the run's
statistics log), cost (wall, CPU, peak RSS of the andes process tree), the true FDP where the
database carries an entrapment, and the precursor ratio accuracy of an --lfq --mbr design. A
failed run is recorded too (status "failed"), so a regression shows up on the page. The CI page
(website/pages/ci.md) reads these lines.

Needs python3, pyarrow for quant datasets, and Percolator (on PATH, else its Docker image).
"""
import argparse
import glob
import json
import math
import os
import re
import statistics as st
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import entrap_common as E  # noqa: E402

SPECIES = ("HUMAN", "YEAST", "ECOLI")


def git(*a):
    return subprocess.run(["git", "-C", str(HERE), *a], capture_output=True, text=True).stdout.strip()


def spectra(d):
    """Native .raw first (see ../reproduce/README.md), else mzML."""
    for ext in ("raw", "mzML", "mzml", "mzML.gz", "d"):
        f = sorted(glob.glob(f"{d}/*.{ext}"))
        if f:
            return f
    return []


def stat(log, label):
    """First integer after `label :` in the statistics log (thousands separators allowed)."""
    m = re.search(rf"{re.escape(label)}\s*:\s*([\d,]+)(?:\s*/\s*([\d,]+))?", log)
    return [int(g.replace(",", "")) for g in m.groups() if g] if m else None


def run(cmd, log):
    """(exit code, wall s, CPU s, peak RSS GB) of cmd and every process it waited for."""
    t0 = time.time()
    with open(log, "w") as fh:
        p = subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT)
        _, status, ru = os.wait4(p.pid, 0)
    rss_kb = ru.ru_maxrss * (1 / 1024 if sys.platform == "darwin" else 1)  # macOS reports bytes
    return os.waitstatus_to_exitcode(status), time.time() - t0, ru.ru_utime + ru.ru_stime, rss_kb / 1024**2


def quantiles(xs):
    q = st.quantiles(xs, n=4, method="inclusive") if len(xs) > 1 else [xs[0]] * 3
    iqr = q[2] - q[0]
    lo = min(x for x in xs if x >= q[0] - 1.5 * iqr)
    hi = max(x for x in xs if x <= q[2] + 1.5 * iqr)
    return q, lo, hi


def quant(mbr, quant_dir, design, amounts):
    """Precursor log2 ratios, condition / reference, per species; replicate CVs.

    Rows: identified, plus transfers at transfer_q <= 0.01 (the reported view), library_q <= 0.01,
    positive quantity_normalized. A precursor needs >= 2 observed runs per condition; one whose
    proteins span several species is dropped. Ratio = difference of mean log2 per condition.
    """
    import pyarrow.parquet as pq

    runs = dict(r.split("\t") for r in Path(design).read_text().split("\n")[1:] if r.strip())
    conds = sorted(set(runs.values()))
    ref, others = conds[-1], conds[:-1]
    amt = {}
    for r in Path(amounts).read_text().split("\n")[1:]:
        if r.strip():
            sp, c, a = r.split("\t")
            amt[(sp, c)] = float(a)

    def expected(sp, c):
        a = lambda cc: amt.get((sp, cc), amt.get((sp, "*")))  # noqa: E731
        return math.log2(a(c) / a(ref)) if a(c) and a(ref) else None

    species = {}
    for f in glob.glob(f"{quant_dir}/*.quant.parquet"):
        t = pq.read_table(f, columns=["precursor_key", "proteins", "is_decoy"]).to_pydict()
        for k, prot, dec in zip(t["precursor_key"], t["proteins"], t["is_decoy"]):
            if not dec and k not in species:
                s = {sp for sp in SPECIES if f"_{sp}" in (prot or "")}
                species[k] = s.pop() if len(s) == 1 else None

    t = pq.read_table(mbr, columns=["precursor_key", "run", "evidence_type", "quantity_normalized",
                                    "transfer_q", "library_q"]).to_pydict()
    vals = {}
    for k, run_, ev, v, tq, lq in zip(*t.values()):
        cond = runs.get(os.path.basename(run_).split(".")[0])
        transfer = not ev.startswith("identified")
        if (cond is None or ev.endswith("below_threshold") or ev == "shared_feature" or v is None
                or not v > 0 or (lq or 1) > 0.01 or (transfer and (tq is None or tq > 0.01))):
            continue
        vals.setdefault(k, {}).setdefault(cond, []).append(v)

    groups, cvs = [], []
    for sp in SPECIES:
        rows = [v for k, v in vals.items() if species.get(k) == sp]
        for c in conds:
            cv = [st.stdev(x) / st.mean(x) for r in rows if len(x := r.get(c, [])) >= 2]
            if cv:
                cvs.append({"group": sp, "condition": c, "n": len(cv), "median_cv": st.median(cv)})
        for c in others:
            lr = [st.mean(map(math.log2, r[c])) - st.mean(map(math.log2, r[ref]))
                  for r in rows if len(r.get(c, [])) >= 2 and len(r.get(ref, [])) >= 2]
            if not lr:
                continue
            q, lo, hi = quantiles(lr)
            e = expected(sp, c)
            groups.append({
                "group": sp, "comparison": f"{c}/{ref}", "n": len(lr), "median": q[1], "q1": q[0],
                "q3": q[2], "lo": lo, "hi": hi, "expected": e,
                "bias": None if e is None else q[1] - e,
                "mae": None if e is None else st.median(abs(x - e) for x in lr),
            })
    return {"quantified": len(vals), "ratios": groups, "cv": cvs,
            "level": "precursor", "rows": "identified + transfers at transfer_q <= 0.01"}


def bench(name, ds, data, andes, threads, workdir):
    out = Path(workdir) / name
    out.mkdir(parents=True, exist_ok=True)
    rec = {"dataset": name, "description": ds.get("description", "")}
    specs, fasta = spectra(Path(data) / ds["fetch"]), Path(data) / "databases" / ds["database"]
    if not specs or not fasta.is_file():
        return {**rec, "status": "missing", "error": f"no spectra in {data}/{ds['fetch']} or no {fasta}"}
    stem = out / "run"
    cmd = [andes, *sum((["--spectrum", s] for s in specs), []), "--database", str(fasta),
           "--threads", str(threads), "--rescore", "--output-pin", f"{stem}.pin", *ds.get("args", [])]
    if "quant" in ds:
        cmd += ["--quant-dir", str(out / "quant")]
    code, wall, cpu, rss = run(cmd, out / "andes.log")
    rec.update(wall_s=round(wall, 1), cpu_s=round(cpu, 1), peak_rss_gb=round(rss, 2), threads=threads)
    log_path = Path(f"{stem}.statistics.log")
    if code != 0 or not log_path.is_file():
        tail = (out / "andes.log").read_text(errors="replace").splitlines()[-5:]
        return {**rec, "status": "failed", "error": f"exit {code}: " + " | ".join(tail)}
    log = log_path.read_text()
    psms, prec, groups = stat(log, "Identified PSMs"), stat(log, "Precursors / peptidoforms"), \
        stat(log, "Protein groups (q ≤ 0.01)")
    seqs = re.search(r"sequences:\s*([\d,]+)", log)
    n_spec = stat(log, "Spectra searched")
    rec.update(status="ok", spectra=n_spec and n_spec[0], psms=psms and psms[0],
               precursors=prec and prec[0], peptides=seqs and int(seqs.group(1).replace(",", "")),
               protein_groups=groups and groups[0])
    rec["throughput"] = round(rec["spectra"] / wall, 1) if rec["spectra"] else None
    if ds.get("entrapment"):
        n, hits, _ = E.count_percolator_psms(f"{stem}.percolator.target.psms.txt", 0.01,
                                             E.SampleSequences(str(fasta)))
        t, e, _ = E.te_ratio(str(fasta))
        rec["true_fdp"] = round(hits / max(n, 1) * (1 + t / max(e, 1)), 5)
    if "quant" in ds:
        q = ds["quant"]
        rec["quant"] = quant(f"{stem}.mbr.parquet", out / "quant", HERE / q["design"], HERE / q["amounts"])
    return rec


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("datasets", nargs="*")
    ap.add_argument("--data", required=True)
    ap.add_argument("--andes", default=os.environ.get("ANDES", "andes"))
    ap.add_argument("--out", required=True, help="JSON-lines file to append to")
    ap.add_argument("--workdir", default="ci-work")
    ap.add_argument("--threads", type=int, default=int(os.environ.get("THREADS", os.cpu_count() or 4)))
    ap.add_argument("--host", default=os.environ.get("BENCH_HOST", ""), help="label of the machine")
    a = ap.parse_args(argv)
    all_ds = {k: v for k, v in json.loads((HERE / "datasets.json").read_text()).items() if not k.startswith("_")}
    names = a.datasets or list(all_ds)
    unknown = set(names) - set(all_ds)
    if unknown:
        ap.error(f"unknown dataset(s) {sorted(unknown)}; known: {sorted(all_ds)}")
    commit = {"commit": git("rev-parse", "HEAD"), "date": git("log", "-1", "--format=%cI"),
              "message": git("log", "-1", "--format=%s"), "run_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "host": a.host}
    bad = 0
    for n in names:
        rec = {**commit, **bench(n, all_ds[n], a.data, a.andes, a.threads, a.workdir)}
        bad += rec["status"] != "ok"
        with open(a.out, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), file=sys.stderr)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
