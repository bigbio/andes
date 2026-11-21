"""arms.py <dia_quant.parquet> <out_dir>: per-arm variants of the DIA quantity table for scoring.
Each arm is normalised by the rule of `andes quant finalize` (per-run factor 2^median over the precursors
quantified in every run of log2 q - the precursor's mean); the Python rule is first checked against the
Rust quantity_normalized of the default (fragment) arm. Writes <out>/<arm>.dia_quant.parquet with
quantity / quantity_normalized set to that arm (raw and normalised)."""
import os, sys
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
src, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
t = pq.read_table(src); d = t.to_pandas()
runs = sorted(d.filename.unique())
def factors(col):
    w = d.pivot_table(index=["peptide", "charge"], columns="filename", values=col, aggfunc="first")
    full = w[(w > 0).all(axis=1)].dropna()
    l = np.log2(full)
    dev = l.sub(l.mean(axis=1), axis=0)
    return {r: float(2 ** dev[r].median()) for r in runs}, len(full)
ff, nf = factors("quantity_frag")
chk = d.quantity_frag / d.filename.map(ff)
m = d.quantity_normalized.notna()
rel = float((chk[m] / d.quantity_normalized[m] - 1).abs().max())
print(f"fragment factors {[round(ff[r], 4) for r in runs]} over {nf} shared; max rel diff vs Rust {rel:.2e}")
assert rel < 1e-6, "python normalisation does not reproduce the Rust one"
for arm in ("frag", "ms1"):
    f, n = factors(f"quantity_{arm}")
    e = d.copy()
    e["quantity"] = e[f"quantity_{arm}"]
    e["quantity_normalized"] = e["quantity"] / e.filename.map(f)
    print(f"{arm}: factors {[round(f[r], 4) for r in runs]} over {n} shared; rows with a quantity {int(e.quantity.notna().sum())}")
    pq.write_table(pa.Table.from_pandas(e, preserve_index=False), os.path.join(out, f"{arm}.dia_quant.parquet"))
