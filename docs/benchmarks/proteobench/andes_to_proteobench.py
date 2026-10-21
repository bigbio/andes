#!/usr/bin/env python3
"""Convert andes match-between-runs output to the ProteoBench "Custom" table.

Reads:
  --mbr        <stem>.mbr.parquet from `andes quant finalize -o <stem>` (or `andes --lfq --mbr`):
               one row per (precursor, run) with evidence_type, quantity and quantity_normalized.
  --quant-dir  the quant directory of that run (`--quant-dir`): every <run>.quant.parquet holds
               the peptidoform, charge and proteins of each precursor key. The in-process
               `andes --lfq --mbr` run keeps it only when `--quant-dir` is given.

Writes one tab-separated row per precursor (peptidoform and charge):
  Sequence, Modified sequence (ProForma, Unimod names), Proteins (';'-separated), Charge,
  then one intensity column per run, named as the ProteoBench module expects (empty = missing).

Rows kept: the evidence types in --evidence (default: identified, transfer_ms2, transfer_ms1,
the rows that passed andes's transfer q-value); `below_threshold` and
`identified_below_threshold` rows are never reported. A non-positive or missing quantity is
missing. ProteoBench itself drops contaminants (`Cont_`) and precursors whose proteins span
more than one species, so every target precursor is written.

Requirements: python3 and pyarrow. Example:
  python3 andes_to_proteobench.py --mbr pb.mbr.parquet --quant-dir q/ --module qe \
      --column quantity_normalized -o custom_input.tsv
"""
import argparse
import csv
import glob
import math
import os
import re
import sys

# Run names (file stems) of the two ProteoBench DDA ion modules, in column order.
MODULE_RUNS = {
    "qe": [f"LFQ_Orbitrap_DDA_Condition_{c}_Sample_Alpha_0{i}" for c in "AB" for i in (1, 2, 3)],
    "astral": [f"LFQ_Astral_DDA_15min_50ng_Condition_{c}_REP{i}" for c in "AB" for i in (1, 2, 3)],
}
REPORTED = ("identified", "transfer_ms2", "transfer_ms1")

# andes writes a modification as [UNIMOD:<n>] when the mods file gives an accession, else as a
# signed mass delta. A protein N-terminal acetyl on a residue that also carries a residue mod is
# written as one summed delta on the first residue (e.g. C[+99.03203] = Acetyl + Carbamidomethyl).
UNIMOD = {"UNIMOD:1": "Acetyl", "UNIMOD:4": "Carbamidomethyl", "UNIMOD:35": "Oxidation"}
DELTAS = {  # mass delta -> the Unimod names it stands for
    42.010565: ["Acetyl"],
    57.021464: ["Carbamidomethyl"],
    15.994915: ["Oxidation"],
    99.032029: ["Acetyl", "Carbamidomethyl"],
    58.005480: ["Acetyl", "Oxidation"],
}
TOKEN = re.compile(r"([A-Z])((?:\[[^\]]+\])*)")


def mod_names(tag):
    """Unimod names for one bracketed tag of an andes peptidoform."""
    if tag.upper() in UNIMOD:
        return [UNIMOD[tag.upper()]]
    try:
        delta = float(tag)
    except ValueError:
        raise ValueError(f"unknown modification tag [{tag}]") from None
    for mass, names in DELTAS.items():
        if abs(mass - delta) < 0.002:
            return names
    raise ValueError(f"unknown modification delta [{tag}]")


def to_proforma(peptidoform):
    """andes peptidoform -> (bare sequence, ProForma with Unimod names).

    Acetyl on the first residue is the protein N-terminal acetylation, written as the
    ProForma N-terminal prefix `[Acetyl]-`.
    """
    toks = TOKEN.findall(peptidoform)
    if "".join(a + m for a, m in toks) != peptidoform or not toks:
        raise ValueError(f"cannot parse peptidoform {peptidoform!r}")
    nterm, out = "", []
    for i, (aa, mods) in enumerate(toks):
        res = aa
        for tag in re.findall(r"\[([^\]]+)\]", mods):
            for name in mod_names(tag):
                if name == "Acetyl" and i == 0:
                    nterm = "[Acetyl]-"
                else:
                    res += f"[{name}]"
        out.append(res)
    return "".join(a for a, _ in toks), nterm + "".join(out)


def run_name(run):
    """mbr.parquet `run` -> file stem (drops a directory and every extension)."""
    return os.path.basename(run).split(".")[0]


def read_metadata(quant_dir):
    """precursor_key -> (peptidoform, charge, proteins) from the target PSMs of every run."""
    import pyarrow.parquet as pq

    files = sorted(glob.glob(os.path.join(quant_dir, "*.quant.parquet")))
    if not files:
        raise SystemExit(f"no *.quant.parquet in {quant_dir}")
    meta = {}
    cols = ["precursor_key", "peptidoform", "charge", "proteins", "is_decoy"]
    for f in files:
        for batch in pq.ParquetFile(f).iter_batches(columns=cols, batch_size=200_000):
            d = batch.to_pydict()
            for k, pf, z, prot, dec in zip(*(d[c] for c in cols)):
                if not dec and k not in meta:
                    meta[k] = (pf, int(z), prot)
    return meta


def read_quantities(mbr_path, column, evidence):
    """precursor_key -> {run stem: value} for the kept rows with a positive finite value."""
    import pyarrow.parquet as pq

    t = pq.read_table(mbr_path, columns=["precursor_key", "run", "evidence_type", column]).to_pydict()
    vals, counts = {}, {}
    for k, run, ev, v in zip(t["precursor_key"], t["run"], t["evidence_type"], t[column]):
        counts[ev] = counts.get(ev, 0) + 1
        if ev not in evidence or v is None or not (math.isfinite(v) and v > 0):
            continue
        vals.setdefault(k, {})[run_name(run)] = v
    return vals, counts


def write_custom(out, vals, meta, runs):
    """Write the Custom table; returns (precursors written, keys without PSM metadata)."""
    missing = written = 0
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["Sequence", "Modified sequence", "Proteins", "Charge"] + runs)
        for k in sorted(vals):
            row_vals = [vals[k].get(r, "") for r in runs]
            if all(v == "" for v in row_vals):
                continue
            if k not in meta:
                missing += 1
                continue
            pf, z, prot = meta[k]
            seq, proforma = to_proforma(pf)
            proteins = ";".join(p for p in re.split(r"[;,]", prot) if p)
            w.writerow([seq, proforma, proteins, z] + row_vals)
            written += 1
    return written, missing


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--mbr", required=True, help="<stem>.mbr.parquet")
    ap.add_argument("--quant-dir", required=True, help="the --quant-dir with <run>.quant.parquet")
    ap.add_argument("--module", choices=sorted(MODULE_RUNS), help="ProteoBench module (sets the run columns)")
    ap.add_argument("--runs", help="comma-separated run stems, in column order (instead of --module)")
    ap.add_argument("--column", default="quantity_normalized", choices=["quantity", "quantity_normalized"])
    ap.add_argument("--evidence", default=",".join(REPORTED),
                    help="comma-separated evidence types to keep (default: %(default)s)")
    ap.add_argument("-o", "--output", required=True)
    a = ap.parse_args(argv)
    if bool(a.module) == bool(a.runs):
        ap.error("give exactly one of --module and --runs")
    runs = MODULE_RUNS[a.module] if a.module else a.runs.split(",")
    evidence = set(a.evidence.split(","))
    vals, counts = read_quantities(a.mbr, a.column, evidence)
    seen = {r for v in vals.values() for r in v}
    if not seen & set(runs):
        raise SystemExit(f"no run of {a.mbr} matches the columns {runs}; found {sorted(seen)[:6]}")
    written, missing = write_custom(a.output, vals, read_metadata(a.quant_dir), runs)
    print(f"{a.output}: {written} precursors, column {a.column}, evidence {sorted(evidence)}; "
          f"mbr rows by evidence {dict(sorted(counts.items()))}; {missing} keys without PSM metadata",
          file=sys.stderr)
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
