#!/usr/bin/env python3
"""MSNet RAW parquet -> intensity context model (Phase T1).

Reads annotated b/y peaks from MSNet parquets (needs ion_type_array; NOT flat3/).
Accumulates per-context-key log relative intensity stats for the strong-score
numerator. Reuses the quality filters and reader from msnet_to_flat.py.

Usage:
  python3 msnet_intensity_agg.py <in_parquet> <out_partial.parquet> [pep_max=0.001] [sample_rows]
  python3 msnet_intensity_agg.py --merge <partial1> <partial2> ... -o intensity_model.parquet

  sample_rows: optional DuckDB USING SAMPLE N ROWS cap (e.g. 300000). Omit for full scan.
  THREADS env var overrides DuckDB thread count (default: min(8, cpu_count)).

Partial output columns: ion_type, flank_n, flank_c, pos_bin, charge, nce_bin,
count, sum_log_rel, sum_log_rel_sq.

Final output adds mean_log_rel and var_log_rel (population variance).
"""
from __future__ import annotations

import os
import re
import sys
from collections import defaultdict
from typing import DefaultDict, Iterable, Optional, Tuple

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

# Reuse reader + mod resolution from the flat converter.
from msnet_to_flat import build_query, resolve_mods

_KEY = Tuple[str, str, str, int, int, str]
Agg = DefaultDict[_KEY, list]  # [count, sum_log, sum_sq]

PARTIAL_SCHEMA = pa.schema([
    ("ion_type", pa.string()),
    ("flank_n", pa.string()),
    ("flank_c", pa.string()),
    ("pos_bin", pa.int32()),
    ("charge", pa.int32()),
    ("nce_bin", pa.string()),
    ("count", pa.int64()),
    ("sum_log_rel", pa.float64()),
    ("sum_log_rel_sq", pa.float64()),
])

FINAL_SCHEMA = pa.schema([
    ("ion_type", pa.string()),
    ("flank_n", pa.string()),
    ("flank_c", pa.string()),
    ("pos_bin", pa.int32()),
    ("charge", pa.int32()),
    ("nce_bin", pa.string()),
    ("count", pa.int64()),
    ("mean_log_rel", pa.float64()),
    ("var_log_rel", pa.float64()),
])

_ION_RE = re.compile(r"^([by])(\d+)")


class ProgressBar:
    """Lightweight stderr progress bar (TTY redraw; periodic lines when logged)."""

    def __init__(self, total: int, label: str = "", width: int = 36) -> None:
        self.total = max(int(total), 1)
        self.label = label
        self.width = width
        self.current = 0
        self._tty = sys.stderr.isatty()
        self._last_logged = -1

    def update(self, n: int = 1) -> None:
        self.current = min(self.current + n, self.total)
        pct = 100.0 * self.current / self.total
        if self._tty:
            filled = int(self.width * self.current / self.total)
            bar = "=" * filled + "-" * (self.width - filled)
            print(
                f"\r{self.label} [{bar}] {pct:5.1f}% {self.current:,}/{self.total:,}",
                end="",
                file=sys.stderr,
                flush=True,
            )
            return
        step = int(pct // 5) * 5
        if step > self._last_logged:
            self._last_logged = step
            print(
                f"{self.label} {step:3d}% ({self.current:,}/{self.total:,})",
                file=sys.stderr,
                flush=True,
            )

    def finish(self) -> None:
        self.current = self.total
        if self._tty:
            self.update(0)
            print(file=sys.stderr, flush=True)
        else:
            print(
                f"{self.label} 100% ({self.total:,}/{self.total:,})",
                file=sys.stderr,
                flush=True,
            )


def nce_bin_from_cvp(cvp) -> str:
    """Bucket normalized collision energy from cv_params; unknown if absent."""
    if cvp is None:
        return "unknown"
    try:
        if isinstance(cvp, dict):
            for k, v in cvp.items():
                kn = (k or "").lower()
                if "nce" in kn or "collision energy" in kn:
                    val = float(v)
                    return f"{int(round(val / 5) * 5)}"
        if isinstance(cvp, list):
            for kv in cvp:
                name = (kv.get("cv_name") or "").lower()
                if "nce" in name or "collision energy" in name:
                    val = float(kv.get("cv_value") or 0)
                    return f"{int(round(val / 5) * 5)}"
    except (TypeError, ValueError):
        return "unknown"
    return "unknown"


def flank_pair(seq: str, ion_kind: str, idx: int) -> Optional[Tuple[str, str]]:
    """Residues N- and C-side of the cleavage for b_idx or y_idx (1-based)."""
    n = len(seq)
    if idx < 1 or idx >= n:
        return None
    if ion_kind == "b":
        # Cleavage after residue idx.
        return seq[idx - 1], seq[idx]
    # y_idx: cleavage between residue (n-idx) and (n-idx+1).
    left = n - idx
    return seq[left - 1], seq[left]


def pos_bin(idx: int, pep_len: int) -> int:
    if pep_len <= 0:
        return 0
    return int(round(10.0 * idx / pep_len))


def accumulate_psm(
    agg: Agg,
    seq: str,
    charge: int,
    nce_bin: str,
    intensities: Iterable[float],
    iontypes: Iterable[Optional[str]],
) -> int:
    """Add annotated b/y peaks from one PSM; return number of peaks kept."""
    n = len(seq)
    ints = list(intensities or [])
    its = list(iontypes or [])
    if not ints:
        return 0
    base_peak = max(float(x) for x in ints)
    if base_peak <= 0:
        return 0
    kept = 0
    for t, raw_i in zip(its, ints):
        if not t:
            continue
        m = _ION_RE.match(t)
        if not m:
            continue
        kind, idx_s = m.group(1), int(m.group(2))
        # Start with plain b/y only (skip -NH3/-H2O variants for aggregation).
        if "-" in t:
            continue
        flank = flank_pair(seq, kind, idx_s)
        if flank is None:
            continue
        log_r = __import__("math").log(max(float(raw_i), 1e-12) / base_peak)
        key = (
            kind,
            flank[0],
            flank[1],
            pos_bin(idx_s, n),
            int(charge),
            nce_bin,
        )
        slot = agg[key]
        slot[0] += 1
        slot[1] += log_r
        slot[2] += log_r * log_r
        kept += 1
    return kept


def agg_to_table(agg: Agg, schema: pa.Schema) -> pa.Table:
    cols = {name: [] for name in schema.names}
    for key, (cnt, s, ss) in sorted(agg.items()):
        ion, fn, fc, pb, ch, nce = key
        cols["ion_type"].append(ion)
        cols["flank_n"].append(fn)
        cols["flank_c"].append(fc)
        cols["pos_bin"].append(pb)
        cols["charge"].append(ch)
        cols["nce_bin"].append(nce)
        cols["count"].append(cnt)
        if "sum_log_rel" in cols:
            cols["sum_log_rel"].append(s)
            cols["sum_log_rel_sq"].append(ss)
        else:
            mean = s / cnt
            cols["mean_log_rel"].append(mean)
            cols["var_log_rel"].append(max(ss / cnt - mean * mean, 0.0))
    return pa.table(cols, schema=schema)


def _row_total(con: duckdb.DuckDBPyConnection, inner: str, sample_rows: int | None) -> int:
    count_q = f"SELECT count(*)::BIGINT FROM ({inner})"
    if sample_rows:
        count_q += f" USING SAMPLE {int(sample_rows)} ROWS"
    return int(con.execute(count_q).fetchone()[0])


def run_aggregate(
    url: str,
    out: str,
    pep_max: float,
    sample_rows: int | None = None,
    frag_filter: str = "hcd",
) -> None:
    inner = build_query(url, pep_max, 7, 40, min_consensus=0.0, frag_filter=frag_filter)
    q = f"SELECT * FROM ({inner})"
    if sample_rows:
        q += f" USING SAMPLE {int(sample_rows)} ROWS"
    threads = int(os.environ.get("THREADS", min(8, os.cpu_count() or 4)))
    batch_size = int(os.environ.get("BATCH_SIZE", 20000))
    con = duckdb.connect()
    con.execute("INSTALL httpfs; LOAD httpfs;")
    con.execute(
        f"SET memory_limit='6GB'; SET threads={threads}; SET preserve_insertion_order=false;"
    )
    sample_note = f" sample={sample_rows:,}" if sample_rows else " sample=full"
    print(
        f"aggregating intensity keys from {url} "
        f"(PEP<={pep_max} frag={frag_filter}{sample_note} threads={threads})"
    )
    total_rows = _row_total(con, inner, sample_rows)
    print(f"rows_to_process={total_rows:,}", file=sys.stderr, flush=True)
    rel = con.execute(q)
    agg: Agg = defaultdict(lambda: [0, 0.0, 0.0])
    psms = peaks = skipped = 0
    progress = ProgressBar(total_rows, label="aggregate")
    while True:
        batch = rel.fetchmany(batch_size)
        if not batch:
            break
        for seq, charge, _prec_mz, mods, _mz, inten, iontypes, _charge_arr, cvp in batch:
            if resolve_mods(seq, mods) is None:
                skipped += 1
            else:
                nce = nce_bin_from_cvp(cvp)
                n = accumulate_psm(agg, seq, int(charge), nce, inten, iontypes)
                if n:
                    psms += 1
                    peaks += n
                else:
                    skipped += 1
            progress.update(1)
    progress.finish()
    table = agg_to_table(agg, PARTIAL_SCHEMA)
    pq.write_table(table, out, compression="zstd")
    print(f"psms={psms:,} peaks={peaks:,} keys={len(agg):,} skipped={skipped:,} -> {out}")


def merge_partials(paths: list[str], out: str) -> None:
    agg: Agg = defaultdict(lambda: [0, 0.0, 0.0])
    merge_progress = ProgressBar(len(paths), label="merge")
    for path in paths:
        t = pq.read_table(path, columns=PARTIAL_SCHEMA.names)
        for row in t.to_pylist():
            key = (
                row["ion_type"],
                row["flank_n"],
                row["flank_c"],
                int(row["pos_bin"]),
                int(row["charge"]),
                row["nce_bin"],
            )
            slot = agg[key]
            slot[0] += int(row["count"])
            slot[1] += float(row["sum_log_rel"])
            slot[2] += float(row["sum_log_rel_sq"])
        print(f"merged {path} ({t.num_rows:,} keys)")
        merge_progress.update(1)
    merge_progress.finish()
    table = agg_to_table(agg, FINAL_SCHEMA)
    pq.write_table(table, out, compression="zstd")
    print(f"final keys={len(agg):,} -> {out}")


def main() -> None:
    if len(sys.argv) >= 2 and sys.argv[1] == "--merge":
        if "-o" not in sys.argv:
            sys.exit("usage: msnet_intensity_agg.py --merge <partial...> -o <out.parquet>")
        o_idx = sys.argv.index("-o")
        out = sys.argv[o_idx + 1]
        paths = [p for p in sys.argv[2:o_idx] if p != "--merge"]
        if not paths:
            sys.exit("no partial parquets given")
        merge_partials(paths, out)
        return
    if len(sys.argv) < 3:
        sys.exit(
            "usage: msnet_intensity_agg.py <in> <out_partial> [pep_max] [sample_rows]\n"
            "   or: msnet_intensity_agg.py --merge <partial...> -o <final>"
        )
    url, out = sys.argv[1], sys.argv[2]
    pep_max = float(sys.argv[3]) if len(sys.argv) > 3 else 0.001
    sample_rows = int(sys.argv[4]) if len(sys.argv) > 4 else None
    frag_filter = (
        sys.argv[5] if len(sys.argv) > 5 else os.environ.get("FRAG_FILTER", "hcd")
    )
    run_aggregate(url, out, pep_max, sample_rows, frag_filter=frag_filter)


if __name__ == "__main__":
    main()
