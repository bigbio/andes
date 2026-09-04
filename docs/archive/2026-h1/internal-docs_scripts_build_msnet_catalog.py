#!/usr/bin/env python3
"""Build a first-pass MSNet dataset catalog for model selection.

The PRIDE MSNet collection stores one large PSM parquet plus small metadata parquets
per dataset folder. This script avoids downloading the large parquet files; DuckDB
uses HTTP range reads and only scans the small columns needed for counts.
"""

from __future__ import annotations

import argparse
import csv
import html.parser
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

import duckdb


DEFAULT_BASE_URL = (
    "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/"
    "quantms-collections/msnet/"
)
DEFAULT_OUT_DIR = "internal-docs/msnet-model-catalog"


class IndexParser(html.parser.HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag != "a":
            return
        attrs_dict = dict(attrs)
        href = attrs_dict.get("href")
        if href:
            self.links.append(href)


def fetch_text(url: str, timeout: int = 60) -> str:
    with urllib.request.urlopen(url, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="replace")


def list_links(url: str) -> list[str]:
    parser = IndexParser()
    parser.feed(fetch_text(url))
    return parser.links


def list_dataset_ids(base_url: str) -> list[str]:
    ids = []
    for href in list_links(base_url):
        if (
            not href.endswith("/")
            or href.startswith("?")
            or href.startswith("/")
            or href == "../"
        ):
            continue
        name = urllib.parse.unquote(href.rstrip("/"))
        if name and name not in {"Parent Directory"}:
            ids.append(name)
    return sorted(set(ids))


def list_dataset_files(dataset_url: str) -> list[str]:
    files = []
    for href in list_links(dataset_url):
        if href.startswith("?") or href.endswith("/"):
            continue
        name = urllib.parse.unquote(href)
        if name.endswith(".parquet"):
            files.append(name)
    return sorted(set(files))


def parquet_url(dataset_url: str, dataset_id: str, files: list[str], suffix: str) -> str | None:
    def file_url(name: str) -> str:
        return urllib.parse.urljoin(dataset_url, urllib.parse.quote(name))

    if suffix == "MSNet":
        exact = f"{dataset_id}-MSNet.parquet"
        if exact in files:
            return file_url(exact)
        psm_exact = f"{dataset_id}.psm.parquet"
        if psm_exact in files:
            return file_url(psm_exact)
        candidates = [name for name in files if name.endswith("-MSNet.parquet")]
        if not candidates:
            candidates = [name for name in files if name.endswith(".psm.parquet")]
        return file_url(candidates[0]) if candidates else None

    exact = f"{dataset_id}.{suffix}.parquet"
    if exact in files:
        return file_url(exact)
    for name in files:
        if name.endswith(f".{suffix}.parquet"):
            return file_url(name)
    return None


def rows_as_dicts(con: duckdb.DuckDBPyConnection, url: str, columns: str = "*") -> list[dict[str, Any]]:
    cur = con.execute(f"SELECT {columns} FROM read_parquet(?)", [url])
    names = [desc[0] for desc in cur.description]
    return [dict(zip(names, row)) for row in cur.fetchall()]


def one_row(con: duckdb.DuckDBPyConnection, url: str) -> dict[str, Any]:
    cur = con.execute("SELECT * FROM read_parquet(?) LIMIT 1", [url])
    names = [desc[0] for desc in cur.description]
    row = cur.fetchone()
    return dict(zip(names, row)) if row else {}


def unique_join(values: list[Any]) -> str:
    flat: list[str] = []
    for value in values:
        if value is None:
            continue
        if isinstance(value, (list, tuple)):
            flat.extend(str(v) for v in value if v is not None)
        else:
            flat.append(str(value))
    seen = []
    for item in flat:
        item = item.strip()
        if item and item not in seen:
            seen.append(item)
    return "; ".join(seen)


def normalize_mod(mod: Any) -> dict[str, Any]:
    if mod is None:
        return {}
    if isinstance(mod, dict):
        return mod
    # DuckDB structs can arrive as tuple-like values in some versions. Keep a
    # fallback text representation so the catalog remains useful.
    return {"name": str(mod)}


def format_mod(mod: dict[str, Any]) -> str:
    name = str(mod.get("name") or "").strip()
    accession = str(mod.get("accession") or "").strip()
    position = str(mod.get("position") or "").strip()
    aa = str(mod.get("target_amino_acid") or "").strip()
    label = name or accession or "unknown"
    details = []
    if accession and accession != label:
        details.append(accession)
    if position and position != "None":
        details.append(position)
    if aa and aa != "None":
        details.append(aa)
    return f"{label} ({', '.join(details)})" if details else label


def collect_mods(run_rows: list[dict[str, Any]]) -> tuple[str, str, str]:
    fixed: list[str] = []
    variable: list[str] = []
    all_mods: list[str] = []
    for row in run_rows:
        for raw_mod in row.get("modification_parameters") or []:
            mod = normalize_mod(raw_mod)
            label = format_mod(mod)
            all_mods.append(label)
            if bool(mod.get("fixed")):
                fixed.append(label)
            else:
                variable.append(label)
    return unique_join(fixed), unique_join(variable), classify_ptms(all_mods)


def classify_ptms(mod_labels: list[str]) -> str:
    classes = []

    def add(label: str) -> None:
        if label not in classes:
            classes.append(label)

    for raw_label in mod_labels:
        label = raw_label.lower()
        if "tmt" in label or "isobaric" in label:
            add("tmt")
        if "itraq" in label:
            add("itraq")
        if "phospho" in label:
            add("phospho")
        if "acetyl" in label:
            add("acetyl")
        if "oxidation" in label:
            add("oxidation")
        if "deamid" in label:
            add("deamidation")
        if "carbamidomethyl" in label:
            add("carbamidomethyl")
        elif "methyl" in label:
            add("methyl")
        if any(term in label for term in ["gly-gly", "glygly", "digly"]):
            add("ubiquitin/glygly")
        if any(term in label for term in ["glycan", "glyco", "hexnac"]):
            add("glyco")
        if "pyro" in label:
            add("pyro-glu")
        if "deamid" not in label and any(term in label for term in ["amidated", "amidation"]):
            add("amidation")
    return "; ".join(classes)


def base_accession(dataset_id: str) -> str:
    match = re.search(r"(PXD\d+|IPX\d+)", dataset_id)
    if match:
        return match.group(1)
    return dataset_id


def infer_species_hint(dataset_id: str, species: str) -> str:
    if species:
        return species
    hints = {
        "Homo": "Homo sapiens",
        "Mouse": "Mus musculus",
        "Rat": "Rattus norvegicus",
        "Pig": "Sus scrofa",
        "Horse": "Equus caballus",
        "Danio": "Danio rerio",
        "Ecoli": "Escherichia coli",
    }
    for needle, value in hints.items():
        if needle.lower() in dataset_id.lower():
            return value
    return ""


def normalize_species_value(species: str) -> str:
    if not species:
        return ""
    known = {
        "homo sapiens": "Homo sapiens",
        "mus musculus": "Mus musculus",
        "rattus norvegicus": "Rattus norvegicus",
        "sus scrofa": "Sus scrofa",
        "sus scrofa domesticus": "Sus scrofa domesticus",
        "danio rerio": "Danio rerio",
        "gallus gallus": "Gallus gallus",
        "arabidopsis thaliana": "Arabidopsis thaliana",
        "saccharomyces cerevisiae": "Saccharomyces cerevisiae",
        "escherichia coli": "Escherichia coli",
        "unknown": "unknown",
    }
    parts = []
    for raw_part in species.split(";"):
        part = raw_part.strip()
        lower = part.lower()
        if not part:
            continue
        if lower in known:
            parts.append(known[lower])
            continue
        words = part.split()
        if len(words) >= 2:
            parts.append(" ".join([words[0].capitalize(), *[word.lower() for word in words[1:]]]))
        else:
            parts.append(part)
    return unique_join(parts)


def schema_columns(con: duckdb.DuckDBPyConnection, url: str) -> set[str]:
    rows = con.execute("DESCRIBE SELECT * FROM read_parquet(?)", [url]).fetchall()
    return {str(row[0]) for row in rows}


def charge_counts(con: duckdb.DuckDBPyConnection, msnet_url: str, with_pep_counts: bool) -> dict[str, Any]:
    columns = schema_columns(con, msnet_url)
    if "charge" not in columns:
        return psm_counts_without_charge(con, msnet_url, columns)

    decoy_expr = "coalesce(is_decoy, false)" if "is_decoy" in columns else "false"
    pep_expr = (
        "sum(CASE WHEN posterior_error_probability <= 0.01 THEN 1 ELSE 0 END)::BIGINT AS pep_le_0_01"
        if with_pep_counts and "posterior_error_probability" in columns
        else "0::BIGINT AS pep_le_0_01"
    )
    query = f"""
        SELECT
            charge,
            {decoy_expr} AS is_decoy,
            count(*)::BIGINT AS n,
            {pep_expr}
        FROM read_parquet(?)
        GROUP BY charge, {decoy_expr}
        ORDER BY charge, is_decoy
    """
    rows = con.execute(query, [msnet_url]).fetchall()
    hist: dict[str, int] = {}
    target_hist: dict[str, int] = {}
    decoy_hist: dict[str, int] = {}
    pep_le_0_01 = 0
    for charge, is_decoy, n, pep_count in rows:
        key = "null" if charge is None else str(int(charge))
        hist[key] = hist.get(key, 0) + int(n)
        if is_decoy:
            decoy_hist[key] = decoy_hist.get(key, 0) + int(n)
        else:
            target_hist[key] = target_hist.get(key, 0) + int(n)
        pep_le_0_01 += int(pep_count or 0)

    def c(histogram: dict[str, int], key: str) -> int:
        return int(histogram.get(key, 0))

    def x(histogram: dict[str, int]) -> int:
        return sum(v for k, v in histogram.items() if k not in {"2", "3"})

    return {
        "psm_rows": sum(hist.values()),
        "target_psm_rows": sum(target_hist.values()),
        "decoy_psm_rows": sum(decoy_hist.values()),
        "charge_2_psms": c(hist, "2"),
        "charge_3_psms": c(hist, "3"),
        "charge_x_psms": x(hist),
        "target_charge_2_psms": c(target_hist, "2"),
        "target_charge_3_psms": c(target_hist, "3"),
        "target_charge_x_psms": x(target_hist),
        "pep_le_0_01_psms": pep_le_0_01,
        "charge_histogram": json.dumps(hist, sort_keys=True),
        "charge_status": "ok",
    }


def psm_counts_without_charge(
    con: duckdb.DuckDBPyConnection,
    msnet_url: str,
    columns: set[str],
) -> dict[str, Any]:
    if "is_decoy" in columns:
        rows = con.execute(
            """
            SELECT coalesce(is_decoy, false) AS is_decoy, count(*)::BIGINT AS n
            FROM read_parquet(?)
            GROUP BY coalesce(is_decoy, false)
            """,
            [msnet_url],
        ).fetchall()
        target = sum(int(n) for is_decoy, n in rows if not is_decoy)
        decoy = sum(int(n) for is_decoy, n in rows if is_decoy)
    else:
        total = con.execute("SELECT count(*)::BIGINT FROM read_parquet(?)", [msnet_url]).fetchone()[0]
        target = int(total)
        decoy = 0
    return {
        "psm_rows": target + decoy,
        "target_psm_rows": target,
        "decoy_psm_rows": decoy,
        "charge_2_psms": "",
        "charge_3_psms": "",
        "charge_x_psms": "",
        "target_charge_2_psms": "",
        "target_charge_3_psms": "",
        "target_charge_x_psms": "",
        "pep_le_0_01_psms": "",
        "charge_histogram": "{}",
        "charge_status": "missing_scalar_charge",
    }


def summarize_dataset(
    con: duckdb.DuckDBPyConnection,
    base_url: str,
    dataset_id: str,
    with_psm_counts: bool,
    with_pep_counts: bool,
) -> dict[str, Any]:
    dataset_url = urllib.parse.urljoin(base_url, urllib.parse.quote(dataset_id) + "/")
    files = list_dataset_files(dataset_url)
    urls = {
        suffix: parquet_url(dataset_url, dataset_id, files, suffix)
        for suffix in ["dataset", "run", "sample", "provenance", "MSNet"]
    }

    out: dict[str, Any] = {
        "dataset_id": dataset_id,
        "base_accession": base_accession(dataset_id),
        "dataset_url": dataset_url,
        "status": "ok",
        "error": "",
    }

    try:
        dataset = one_row(con, urls["dataset"]) if urls.get("dataset") else {}
        run_rows = rows_as_dicts(
            con,
            urls["run"],
            "instrument, enzymes, dissociation_method, modification_parameters",
        ) if urls.get("run") else []
        sample_rows = rows_as_dicts(con, urls["sample"], "organism, organism_part") if urls.get("sample") else []

        fixed_mods, variable_mods, ptm_classes = collect_mods(run_rows)
        species = unique_join([row.get("organism") for row in sample_rows])
        species = infer_species_hint(dataset_id, species)

        out.update(
            {
                "project_accession": dataset.get("project_accession") or "",
                "project_title": dataset.get("project_title") or "",
                "software_name": dataset.get("software_name") or "",
                "software_version": dataset.get("software_version") or "",
                "run_count": len(run_rows),
                "sample_count": len(sample_rows),
                "species": species,
                "species_normalized": normalize_species_value(species),
                "organism_parts": unique_join([row.get("organism_part") for row in sample_rows]),
                "fragmentation_methods": unique_join([row.get("dissociation_method") for row in run_rows]),
                "instruments": unique_join([row.get("instrument") for row in run_rows]),
                "enzymes": unique_join([row.get("enzymes") for row in run_rows]),
                "fixed_mods": fixed_mods,
                "variable_mods": variable_mods,
                "ptm_classes": ptm_classes,
                "msnet_parquet_url": urls.get("MSNet") or "",
            }
        )

        if with_psm_counts and urls.get("MSNet"):
            out.update(charge_counts(con, urls["MSNet"], with_pep_counts))
        else:
            out.update(
                {
                    "psm_rows": "",
                    "target_psm_rows": "",
                    "decoy_psm_rows": "",
                    "charge_2_psms": "",
                    "charge_3_psms": "",
                    "charge_x_psms": "",
                    "target_charge_2_psms": "",
                    "target_charge_3_psms": "",
                    "target_charge_x_psms": "",
                    "pep_le_0_01_psms": "",
                    "charge_histogram": "",
                    "charge_status": "not_counted",
                }
            )
    except Exception as exc:  # Keep the catalog moving if one remote object is bad.
        out["status"] = "error"
        out["error"] = str(exc).replace("\n", " ")
    return out


FIELDNAMES = [
    "dataset_id",
    "base_accession",
    "project_accession",
    "project_title",
    "species",
    "species_normalized",
    "organism_parts",
    "fragmentation_methods",
    "instruments",
    "enzymes",
    "fixed_mods",
    "variable_mods",
    "ptm_classes",
    "run_count",
    "sample_count",
    "psm_rows",
    "target_psm_rows",
    "decoy_psm_rows",
    "charge_2_psms",
    "charge_3_psms",
    "charge_x_psms",
    "target_charge_2_psms",
    "target_charge_3_psms",
    "target_charge_x_psms",
    "pep_le_0_01_psms",
    "charge_histogram",
    "charge_status",
    "software_name",
    "software_version",
    "dataset_url",
    "msnet_parquet_url",
    "status",
    "error",
]


def write_outputs(rows: list[dict[str, Any]], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "msnet_model_catalog.csv"
    tsv_path = out_dir / "msnet_model_catalog.tsv"
    jsonl_path = out_dir / "msnet_model_catalog.jsonl"

    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    with tsv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    with jsonl_path.open("w") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")

    ok_rows = [row for row in rows if row.get("status") == "ok"]
    total_psms = sum(int(row.get("psm_rows") or 0) for row in ok_rows)
    total_targets = sum(int(row.get("target_psm_rows") or 0) for row in ok_rows)
    total_decoys = sum(int(row.get("decoy_psm_rows") or 0) for row in ok_rows)
    fragmentations = sorted({v for row in ok_rows for v in row.get("fragmentation_methods", "").split("; ") if v})
    enzymes = sorted({v for row in ok_rows for v in row.get("enzymes", "").split("; ") if v})
    species = sorted({v for row in ok_rows for v in row.get("species", "").split("; ") if v})

    report_path = out_dir / "README.md"
    report_path.write_text(
        "\n".join(
            [
                "# MSNet Model Catalog",
                "",
                f"Generated rows: {len(rows)}",
                f"Successful rows: {len(ok_rows)}",
                f"Total PSM rows: {total_psms:,}",
                f"Target PSM rows: {total_targets:,}",
                f"Decoy PSM rows: {total_decoys:,}",
                "",
                "## Dimensions Observed",
                "",
                f"- Fragmentation methods: {', '.join(fragmentations) if fragmentations else 'pending'}",
                f"- Enzymes: {', '.join(enzymes) if enzymes else 'pending'}",
                f"- Species count: {len(species)}",
                "",
                "## Files",
                "",
                f"- `{csv_path.name}`",
                f"- `{tsv_path.name}`",
                f"- `{jsonl_path.name}`",
                "",
                "Columns include fragmentation, species, enzyme, fixed/variable PTMs, total PSM rows,",
                "target/decoy rows, charge-2 rows, charge-3 rows, and charge-x rows.",
                "",
            ]
        )
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-url", default=DEFAULT_BASE_URL)
    parser.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    parser.add_argument("--limit", type=int, default=0, help="Limit datasets for probing.")
    parser.add_argument(
        "--dataset-id",
        action="append",
        default=[],
        help="Dataset folder to process. Can be passed multiple times.",
    )
    parser.add_argument("--skip-psm-counts", action="store_true", help="Only read small metadata parquets.")
    parser.add_argument(
        "--with-pep-counts",
        action="store_true",
        help="Also count PSMs with posterior_error_probability <= 0.01. Slower.",
    )
    args = parser.parse_args()

    base_url = args.base_url.rstrip("/") + "/"
    out_dir = Path(args.out_dir)

    dataset_ids = args.dataset_id or list_dataset_ids(base_url)
    if args.limit:
        dataset_ids = dataset_ids[: args.limit]

    con = duckdb.connect()
    rows: list[dict[str, Any]] = []
    started = time.time()
    for idx, dataset_id in enumerate(dataset_ids, start=1):
        print(f"[{idx}/{len(dataset_ids)}] {dataset_id}", file=sys.stderr, flush=True)
        row = summarize_dataset(
            con,
            base_url,
            dataset_id,
            with_psm_counts=not args.skip_psm_counts,
            with_pep_counts=args.with_pep_counts,
        )
        rows.append({field: row.get(field, "") for field in FIELDNAMES})

    write_outputs(rows, out_dir)
    print(f"Wrote {len(rows)} rows to {out_dir} in {time.time() - started:.1f}s", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
