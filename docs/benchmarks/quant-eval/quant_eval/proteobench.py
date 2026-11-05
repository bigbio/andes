#!/usr/bin/env python3
"""ProteoBench submission files (DDA quantification, precursor-ion modules) from an andes finalize table.

  python -m quant_eval.proteobench --finalize fin.parquet --design design.tsv --fasta hye.fasta --rule hard-cut --out-prefix sub/QE \
      --version 0.4.x --meta key=value ...

Writes `<prefix>_custom_input.tsv` (the module's Custom input: one row per precursor ion, wide format: Sequence, Proforma, Charge, Proteins
(semicolon-separated accessions of every database entry containing the peptide, I = L), then one RAW intensity column per run) and
`<prefix>_params.json` (the free-text metadata fields of the submission form). The module's documented column names and long / wide rule
(a `Raw file` column means long format) were read from its public documentation for the 0.2 - 0.5 series; check them against the module
version at submission time. Rows: identified rows and the transfers allowed by --rule (default: transfer_q <= 0.01; hard-cut: quant_q <= --hard-cut),
library_q <= 0.01. Intensities are the raw `quantity`; ProteoBench's own normalisation applies downstream, as for every tool."""
import argparse
import json
import re
import sys

import pandas as pd

from . import readers as R
from .design import Design
from .fasta import read_fasta
from .keys import bare

MOD = {57.02146: "Carbamidomethyl", 15.99491: "Oxidation", 42.010565: "Acetyl"}


def proforma(key):
    pep, z = key.rsplit("/", 1)
    out, nterm, first = [], "", True
    for tok in re.findall(r"[A-Z]|\[[^\]]*\]", pep):
        if tok.startswith("["):
            d = float(tok[1:-1])
            names = []
            for m, n in MOD.items():
                if abs(d - m) < 0.002:
                    names = [n]
            if not names:      # sums of two or three of the three modifications
                items = list(MOD.items())
                for i, (m1, n1) in enumerate(items):
                    for m2, n2 in items[i:]:
                        if abs(d - m1 - m2) < 0.002:
                            names = [n1, n2]
            if not names:
                raise ValueError(f"unmapped delta {d} in {key}")
            # an N-terminal acetyl is folded onto residue 1 in the key; write it as the ProForma N-terminal modification
            if len(out) == 1 and "Acetyl" in names:
                names.remove("Acetyl")
                nterm = "[Acetyl]-"
            out.append("".join(f"[{n}]" for n in names))
        else:
            out.append(tok)
    return nterm + "".join(out), int(z)


def proteins_of(keys, fasta, k=6):
    ent = [(n, s.replace("I", "L")) for n, s in read_fasta(fasta)]
    idx = {}
    for j, (_, s) in enumerate(ent):
        for i in range(len(s) - k + 1):
            idx.setdefault(s[i:i + k], set()).add(j)
    out, cache = {}, {}
    for key in keys:
        b = bare(key).replace("I", "L")
        if b not in cache:
            cache[b] = ";".join(sorted(ent[j][0] for j in idx.get(b[:k], ()) if b in ent[j][1]))
        out[key] = cache[b]
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--finalize", required=True)
    ap.add_argument("--design", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--rule", default="default", choices=("default", "hard-cut", "all"))
    ap.add_argument("--hard-cut", type=float, default=0.01)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--version", default="")
    ap.add_argument("--meta", action="append", default=[], metavar="KEY=VALUE")
    a = ap.parse_args(argv)
    des = Design(a.design)
    raw = R.read_andes(a.finalize, des.runs, a.rule, a.hard_cut)
    keys = sorted(set(raw.key))
    prot = proteins_of(keys, a.fasta)
    W = raw.pivot_table(index="key", columns="run", values="quantity", aggfunc="first").reindex(columns=des.runs)
    rows, bad = [], 0
    for key in W.index:
        try:
            pf, z = proforma(key)
        except ValueError:
            bad += 1
            continue
        if not prot[key]:
            continue
        rows.append([bare(key), pf, z, prot[key]] + list(W.loc[key].values))
    df = pd.DataFrame(rows, columns=["Sequence", "Proforma", "Charge", "Proteins"] + list(des.runs))
    df.to_csv(a.out_prefix + "_custom_input.tsv", sep="\t", index=False, na_rep="")
    meta = {"software_name": "andes", "software_version": a.version, "search_engine": "andes", "ident_fdr_psm": 0.01, "mbr": True,
            "quantification": "andes quant (MS1 fit, match between runs with decoy-controlled transfer)", "intensity": f"raw precursor quantity, rule {a.rule}"}
    meta.update(dict(m.split("=", 1) for m in a.meta))
    json.dump(meta, open(a.out_prefix + "_params.json", "w"), indent=1)
    print(f"{a.out_prefix}_custom_input.tsv: {len(df)} ions x {len(des.runs)} runs ({bad} unmapped, {len(W) - len(df) - bad} without a database protein)", file=sys.stderr)


if __name__ == "__main__":
    main()
