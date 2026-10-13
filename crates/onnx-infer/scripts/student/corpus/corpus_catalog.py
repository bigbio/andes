#!/usr/bin/env python3
"""Training-corpus catalog for andes-rt / andes-ms2: schema + fillers.

Two tables (parquet, with a TSV twin for reading by eye):

* ``datasets``  one row per (source, accession): provenance, licence, biology,
  acquisition, identification provenance, totals, PTM inventory, distributions, QC.
* ``runs``      one row per (source, accession, run): acquisition as read from the
  run itself (analyser / fragmentation / NCE inferred or parsed), counts and QC.

Sources (``source`` column):
  ``msnet``         MsNet-format flats (peak lists attached; the cluster MsNet collection
                    and the PXD009875 flats are this kind)
  ``qpx``           quantms reanalysis ``psm.parquet`` (IDs + RT, no spectra; spectra via raw)
  ``pride``         a PRIDE submission's own identifications (mzIdentML / mzTab / engine output)
  ``proteometools`` ProteomeTools synthetic-peptide sets (PRIDE submissions, special-cased)
  ``other``         anything else (MassIVE, iProX, jPOST, ...), filled from a manual TSV/JSON

The schema lives in ``FIELDS`` (name, type, level, description); ``spec`` prints it so
downstream scripts can read it. Every filler returns plain dicts keyed by the schema;
missing values are None, never guessed. Values that came from text heuristics carry
the suffix ``_src`` = ``cv`` | ``text`` | ``qc`` | ``manual``.

Usage::

  python corpus_catalog.py spec > SCHEMA.tsv
  python corpus_catalog.py msnet --accession PXD009875 --qc-dir QC/PXD009875 --out-dir CAT [FLAT_DIR_OR_FILES...]
  python corpus_catalog.py qpx --accession PXD012255 --qc-dir QC/PXD012255 --out-dir CAT --sdrf PXD012255.sdrf.tsv
  python corpus_catalog.py pride --classified projects_classified.tsv --files files_summary.tsv \\
         --accessions PXD021013 PXD004732 --out-dir CAT
  python corpus_catalog.py merge --out-dir CAT CAT_PART...       # concatenate + validate
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re
import sys

# Benchmarks / test sets: never trained on.
EXCLUDED = {"PXD070049", "PXD001819", "PXD007683", "PXD007653", "PXD005553", "PXD002047"}
# Accessions in the bundled andes rank/GBDT model store's source ledger (2026-09-05). Training on
# them is fine; *testing* a bundled model on them is circular.
MODEL_STORE_LEDGER = {
    "PXD000900", "PXD002214", "PXD003109", "PXD003702", "PXD004732", "PXD004894", "PXD006939",
    "PXD007058", "PXD007740", "PXD009875", "PXD010429", "PXD010595", "PXD014502", "PXD014525",
    "PXD016766", "PXD016999", "PXD017092", "PXD018176", "PXD019646", "PXD028735", "PXD046453",
    "PXD061135", "PXD065289",
}
# Held out by andes-rt v1 / andes-ms2 v1 (dataset-level test sets).
PREDICTOR_HELDOUT = {"PXD019263", "PXD005235", "PXD012255"}

# (name, type, level, description); level = D (datasets), R (runs), DR (both)
FIELDS = [
    ("source", "str", "DR", "msnet | qpx | pride | proteometools | other"),
    ("accession", "str", "DR", "PXD / MSV / IPX / PMID-derived id of the dataset"),
    ("run", "str", "R", "run (raw file stem)"),
    ("title", "str", "D", "project title"),
    ("year", "int", "D", "publication year"),
    ("licence", "str", "D", "repository licence string (PRIDE: 'EBI terms of use' or CC*)"),
    ("doi", "str", "D", "paper DOI(s), '|'-separated"),
    ("pubmed", "str", "D", "PubMed id(s)"),
    ("organism", "str", "D", "organism name(s) ('|')"),
    ("taxid", "str", "D", "NCBI taxon id(s) ('|')"),
    ("sample_type", "str", "D", "cell line | tissue | plasma/serum | HLA immunoprecipitate | synthetic | microbe | ..."),
    ("enrichment", "str", "DR", "none | phospho | glyco | ubiquitin(GG) | acetyl | HLA-I | HLA-II | ... ('|')"),
    ("labelling", "str", "DR", "LFQ | TMT6/10/11 | TMTpro | iTRAQ4/8 | SILAC | dimethyl"),
    ("enzyme", "str", "DR", "trypsin | LysC | GluC | AspN | chymotrypsin | ArgC | LysN | none (HLA) ('|')"),
    ("enzyme_specificity", "str", "DR", "specific | semi | nonspecific"),
    ("instrument_model", "str", "DR", "instrument model(s) (PSI-MS name)"),
    ("analyser", "str", "DR", "MS2 analyser: FTMS (Orbitrap) | ASTMS (Astral) | ITMS | TOF | timsTOF"),
    ("fragmentation", "str", "DR", "HCD | CID | ETD | EThcD | ETciD | UVPD ('|')"),
    ("activation", "str", "R", "andes-ms2 activation token, e.g. HCD-FTMS, CID-ITMS"),
    ("nce", "float", "DR", "normalised collision energy (median over MS2 scans if read from runs)"),
    ("ms2_resolution", "float", "DR", "MS2 resolving power (Orbitrap, at m/z 200) if known"),
    ("gradient_min", "float", "DR", "LC gradient length in minutes (run: max RT)"),
    ("acquisition", "str", "DR", "DDA | DIA | PRM"),
    ("id_engine", "str", "D", "search engine(s) of the provided identifications"),
    ("id_fdr_method", "str", "D", "e.g. Percolator PEP<=0.01 | target-decoy q<=0.01 | MSFragger expect<=0.01"),
    ("id_format", "str", "D", "mzIdentML | mzTab | PRIDE XML | MaxQuant txt | pepXML | qpx parquet | MsNet flat"),
    ("spectra_available", "str", "DR", "attached | mzML/MGF in submission | raw only | none"),
    ("n_runs", "int", "D", "runs with PSMs"),
    ("n_psms", "int", "DR", "target PSMs after the source's own filter"),
    ("n_unique_peptidoforms", "int", "DR", "distinct ProForma peptidoforms (I kept as given)"),
    ("n_unique_peptidoforms_il", "int", "D", "distinct peptidoforms with I->L merged"),
    ("ptm_inventory", "json", "D", "{'UNIMOD:21': n_psms, ..., 'mass:+x': n} incl. fixed mods"),
    ("n_ptm_types", "int", "D", "distinct non-fixed modifications (excl. UNIMOD:4 and labels)"),
    ("length_hist", "json", "D", "{length: n_psms}, lengths >= 40 pooled at 40"),
    ("charge_hist", "json", "D", "{charge: n_psms}"),
    ("tryptic_cterm_frac", "float", "DR", "fraction of PSMs whose sequence ends in K/R"),
    ("isoleucine_written_as_L", "bool", "D", "no I in any sequence while L present (qpx convention)"),
    # QC (filled from corpus_qc.py outputs)
    ("qc_ann_frac_median", "float", "DR", "median annotated-intensity fraction, b/y 1+/2+ at analyser tolerance"),
    ("qc_ann_frac_p10", "float", "D", "10th percentile of the same"),
    ("qc_low_ann_frac", "float", "DR", "fraction of PSMs with annotated fraction < 0.2"),
    ("qc_n_peaks_median", "float", "DR", "median peaks per MS2"),
    ("qc_prec_err_ppm_median", "float", "DR", "median precursor error (ppm, best isotope)"),
    ("qc_prec_fit_frac", "float", "DR", "fraction of PSMs whose precursor fits the peptidoform within tolerance"),
    ("qc_frag_err_da_median", "float", "DR", "median signed fragment error (Da) of b/y 1+ in 0.5 Da"),
    ("qc_charge_alt_fit_frac", "float", "DR", "fraction whose precursor fits better at z +/- 1"),
    ("qc_decoy_auc", "float", "D", "target-vs-decoy score AUC where decoys exist"),
    ("qc_pep_p99", "float", "D", "99th percentile PEP of retained targets"),
    ("qc_outlier_runs", "int", "D", "runs flagged by robust z / analyser mismatch"),
    ("qc_outlier_flags", "str", "R", "flags of the run"),
    ("analyser_src", "str", "DR", "cv | text | qc | manual"),
    ("fragmentation_src", "str", "DR", "cv | text | qc | manual"),
    # governance
    ("excluded_benchmark", "bool", "D", "accession is an andes benchmark/test set: never train"),
    ("in_model_store_ledger", "bool", "D", "trained a bundled andes model: never use as a test set for it"),
    ("predictor_heldout", "bool", "D", "held out by andes-rt/andes-ms2 v1"),
    ("contamination_note", "str", "D", "free text"),
    ("priority", "str", "D", "user-priority stratum: HLA | PTM | activation | instrument | none"),
    ("notes", "str", "D", "free text"),
]
DATASET_COLS = [f[0] for f in FIELDS if "D" in f[2]]
RUN_COLS = [f[0] for f in FIELDS if "R" in f[2]]
FIXED_OR_LABEL = {"UNIMOD:4", "UNIMOD:737", "UNIMOD:2016", "UNIMOD:214", "UNIMOD:730", "UNIMOD:259", "UNIMOD:267"}


def governance(acc: str) -> dict:
    base = acc.split("-", 1)[0]
    notes = []
    if base in EXCLUDED:
        notes.append("andes benchmark/test set")
    if base in MODEL_STORE_LEDGER:
        notes.append("in bundled model-store ledger")
    if base in PREDICTOR_HELDOUT:
        notes.append("andes-rt/ms2 v1 held-out")
    return {"excluded_benchmark": base in EXCLUDED, "in_model_store_ledger": base in MODEL_STORE_LEDGER,
            "predictor_heldout": base in PREDICTOR_HELDOUT, "contamination_note": "; ".join(notes) or None}


def empty(cols):
    return {c: None for c in cols}


def from_qc(acc: str, source: str, qc_dir: str) -> tuple[dict, list[dict]]:
    """Fill the composition + QC part of a dataset row and its run rows from corpus_qc.py output."""
    d = empty(DATASET_COLS)
    d.update(source=source, accession=acc, **governance(acc))
    S = json.load(open(os.path.join(qc_dir, "dataset_qc.json")))
    inv = S.get("mod_psm_counts", {})
    d.update(
        n_psms=S.get("n_psms"), n_runs=S.get("n_runs"), n_unique_peptidoforms=S.get("unique_peptidoforms"),
        n_unique_peptidoforms_il=S.get("unique_peptidoforms_IL_merged"), ptm_inventory=json.dumps(inv),
        n_ptm_types=len([k for k in inv if k not in FIXED_OR_LABEL]),
        length_hist=json.dumps(S.get("length_hist")), charge_hist=json.dumps(S.get("charge_hist")),
        tryptic_cterm_frac=S.get("tryptic_cterm_frac"), isoleucine_written_as_L=S.get("isoleucine_written_as_L"),
        qc_prec_fit_frac=S.get("prec_fit_frac"), qc_charge_alt_fit_frac=S.get("charge_alt_fit_frac"),
        qc_outlier_runs=len(S.get("outlier_runs") or []),
    )
    if isinstance(S.get("prec_err_ppm_quantiles"), dict):
        d["qc_prec_err_ppm_median"] = S["prec_err_ppm_quantiles"].get("0.5")
    auc = S.get("target_decoy_auc")
    if isinstance(auc, dict) and auc:
        d["qc_decoy_auc"] = max(auc.values())
    if isinstance(S.get("pep_quantiles_targets"), dict):
        d["qc_pep_p99"] = S["pep_quantiles_targets"].get("0.99")
    by = S.get("by_inferred_analyser") or {}
    if by:
        main = max(by, key=lambda k: by[k]["psms"])
        d["analyser"] = "|".join(sorted(by))
        d["analyser_src"] = "qc"
        d["qc_ann_frac_median"] = by[main]["ann_frac_quantiles"]["0.5"]
        d["qc_ann_frac_p10"] = by[main]["ann_frac_quantiles"]["0.1"]
        d["qc_low_ann_frac"] = by[main]["low_ann_frac_lt_0p2"]
        d["qc_n_peaks_median"] = by[main]["n_peaks_quantiles"]["0.5"]
        d["qc_frag_err_da_median"] = by[main]["frag_err_da_quantiles"]["0.5"]
    runs = []
    rp = os.path.join(qc_dir, "run_qc.tsv")
    if os.path.exists(rp):
        for r in csv.DictReader(open(rp), delimiter="\t"):
            row = empty(RUN_COLS)
            run = r.get("run") or r.get("run_file_name")
            an = r.get("analyser_inferred")
            row.update(source=source, accession=acc, run=run, n_psms=int(float(r["n_psms"])),
                       n_unique_peptidoforms=int(float(r.get("unique_peptidoforms") or 0)) or None,
                       analyser=an, analyser_src="qc" if an else None,
                       qc_ann_frac_median=_f(r.get("ann_frac_median")), qc_low_ann_frac=_f(r.get("low_ann_frac")),
                       qc_n_peaks_median=_f(r.get("n_peaks_median")), qc_prec_err_ppm_median=_f(r.get("prec_err_ppm_median")),
                       qc_prec_fit_frac=_f(r.get("prec_fit_frac")), qc_frag_err_da_median=_f(r.get("frag_err_da_median")),
                       qc_charge_alt_fit_frac=_f(r.get("charge_alt_fit_frac")),
                       tryptic_cterm_frac=_f(r.get("tryptic_frac")), qc_outlier_flags=r.get("outlier_flags") or "",
                       gradient_min=(_f(r.get("rt_max")) / 60.0) if r.get("rt_max") else None)
            runs.append(row)
    return d, runs


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


# ----------------------------------------------------------------------------- PRIDE metadata
def pride_rows(classified_tsv: str, files_tsv: str | None, accessions) -> dict:
    want = set(accessions or [])
    files = {}
    if files_tsv and os.path.exists(files_tsv):
        files = {r["accession"]: r for r in csv.DictReader(open(files_tsv), delimiter="\t")}
    out = {}
    for r in csv.DictReader(open(classified_tsv), delimiter="\t"):
        if want and r["accession"] not in want:
            continue
        d = empty(DATASET_COLS)
        d.update(governance(r["accession"]))
        sample = []
        if r["hla"]:
            sample.append("HLA immunoprecipitate")
        if r["plasma"] == "1":
            sample.append("plasma/serum")
        if r["synthetic"] == "1":
            sample.append("synthetic")
        f = files.get(r["accession"], {})
        usable = f.get("usable_without_research")
        d.update(
            source="pride", accession=r["accession"], title=r["title"], year=int(r["year"]) if r["year"] else None,
            licence=r["licence"], doi=r["doi"] or None, pubmed=r["pubmed"] or None, organism=r["organisms"],
            sample_type="|".join(sample) or None,
            enrichment="|".join(x for x in (r["enrichment"], r["hla"].replace("HLA|", "") if r["hla"] else "") if x) or None,
            labelling=r["labelling"], enzyme=r["enzymes_text"] or None,
            enzyme_specificity="nonspecific" if (r["hla"] or "nonspecific" in r["enzymes_text"]) else None,
            instrument_model=r["instruments"], analyser=r["analyser"] or None, analyser_src="cv",
            fragmentation=r["fragmentation_text"] or None, fragmentation_src="text",
            acquisition="DIA" if r["dia"] == "1" else "DDA", id_engine=r["software"] or None,
            id_format=usable, spectra_available=("mzML/MGF in submission" if usable and "+peaks" in usable
                                                 else "raw only" if usable else None),
            ptm_inventory=json.dumps({p.split(":", 2)[0] + ":" + p.split(":", 2)[1]: None for p in r["ptms_cv"].split("|") if p.count(":") >= 2}),
            priority="HLA" if r["hla"] else "PTM" if r["enrichment"] else None,
        )
        out[r["accession"]] = d
    return out


# ----------------------------------------------------------------------------- writers
def write(out_dir, datasets, runs, tag):
    import pandas as pd

    os.makedirs(out_dir, exist_ok=True)
    D = pd.DataFrame(datasets, columns=DATASET_COLS)
    R = pd.DataFrame(runs, columns=RUN_COLS)
    for df, name in ((D, "datasets"), (R, "runs")):
        if len(df):
            df.to_parquet(os.path.join(out_dir, f"{name}.{tag}.parquet"), index=False)
            df.to_csv(os.path.join(out_dir, f"{name}.{tag}.tsv"), sep="\t", index=False)
    print(f"{len(D)} dataset rows, {len(R)} run rows -> {out_dir} ({tag})")


def cmd_spec(a):
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(["field", "type", "level", "description"])
    for f in FIELDS:
        w.writerow(f)


def cmd_from_qc(a, source):
    d, runs = from_qc(a.accession, source, a.qc_dir)
    meta = {}
    if a.meta:
        meta = json.load(open(a.meta))
    if source == "msnet":
        d.update(id_format="MsNet flat", spectra_available="attached")
    if source == "qpx":
        d.update(id_format="qpx parquet", id_engine="quantms (SAGE/Comet/MSGF+ + Percolator)",
                 id_fdr_method="Percolator PEP<=0.01 (as filtered here)", spectra_available="raw only")
    d.update({k: v for k, v in meta.items() if k in d})
    for k in meta:
        if k in d and f"{k}_src" in d:
            d[f"{k}_src"] = "manual"
    for r in runs:
        r.update({k: v for k, v in meta.items() if k in r and r[k] is None})
    write(a.out_dir, [d], runs, f"{source}.{a.accession}")


def cmd_pride(a):
    rows = pride_rows(a.classified, a.files, a.accessions)
    if a.proteometools:
        for acc in a.proteometools:
            if acc in rows:
                rows[acc].update(source="proteometools", sample_type="synthetic",
                                 contamination_note=(rows[acc]["contamination_note"] or "") + " synthetic: no biological overlap with benchmarks")
    write(a.out_dir, list(rows.values()), [], "pride")


def cmd_merge(a):
    import pandas as pd

    for name, cols in (("datasets", DATASET_COLS), ("runs", RUN_COLS)):
        parts = [pd.read_parquet(p) for d in a.parts for p in glob.glob(os.path.join(d, f"{name}.*.parquet"))]
        if not parts:
            continue
        df = pd.concat(parts, ignore_index=True)[cols]
        key = ["source", "accession"] + (["run"] if name == "runs" else [])
        dup = df.duplicated(key, keep="last")
        if dup.any():
            print(f"{name}: {int(dup.sum())} duplicate keys, keeping the last", file=sys.stderr)
        df = df[~dup]
        bad = df[df.accession.str.split("-").str[0].isin(EXCLUDED)] if name == "datasets" else df.iloc[:0]
        if len(bad):
            print(f"WARNING benchmark accessions present (excluded_benchmark=True): {sorted(bad.accession)}", file=sys.stderr)
        os.makedirs(a.out_dir, exist_ok=True)
        df.to_parquet(os.path.join(a.out_dir, f"{name}.parquet"), index=False)
        df.to_csv(os.path.join(a.out_dir, f"{name}.tsv"), sep="\t", index=False)
        print(f"{name}: {len(df)} rows")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("spec")
    for src in ("msnet", "qpx", "other"):
        s = sub.add_parser(src)
        s.add_argument("--accession", required=True)
        s.add_argument("--qc-dir", required=True, help="corpus_qc.py output directory")
        s.add_argument("--meta", help="JSON of schema fields known from the paper/SDRF (instrument_model, nce, ...)")
        s.add_argument("--out-dir", required=True)
    s = sub.add_parser("pride")
    s.add_argument("--classified", required=True)
    s.add_argument("--files")
    s.add_argument("--accessions", nargs="*")
    s.add_argument("--proteometools", nargs="*", default=[])
    s.add_argument("--out-dir", required=True)
    s = sub.add_parser("merge")
    s.add_argument("--out-dir", required=True)
    s.add_argument("parts", nargs="+")
    a = p.parse_args()
    if a.cmd == "spec":
        cmd_spec(a)
    elif a.cmd in ("msnet", "qpx", "other"):
        cmd_from_qc(a, a.cmd)
    elif a.cmd == "pride":
        cmd_pride(a)
    else:
        cmd_merge(a)


if __name__ == "__main__":
    main()
