"""The long table every tool's cross-link quantities are brought to before scoring.

One row per (precursor, run), TSV or Parquet:

| column | meaning |
|---|---|
| `run` | run name: the stem of the raw file (andes writes the input file stem) |
| `precursor` | the tool's own key of the cross-link precursor (peptide pair, sites, modifications, charge) |
| `residue_pair` | residue-pair key in the andes form (`ACC:pos -- ACC:pos`, see keys.py) |
| `quantity` | the raw (not normalised) MS1 quantity; rows without a value are left out |
| `evidence` | `identified` (an MS2 identification in this run) or `transferred` (measured without one) |
| `quant_q` | optional: the tool's quantification q-value of a transferred row |
| `alpha_peptide`, `alpha_site`, `beta_peptide`, `beta_site` | optional: bare peptides and 1-based link sites, needed by the peptide-library truth of the false-transfer designs |
"""
import pandas as pd

from .keys import normalise_key

REQUIRED = ("run", "precursor", "residue_pair", "quantity", "evidence")
EVIDENCE = ("identified", "transferred")


def validate(df, source="table"):
    """Check the contract and return the table with canonical residue-pair keys."""
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f"{source}: missing columns {missing}")
    df = df.copy()
    df["quantity"] = pd.to_numeric(df["quantity"], errors="coerce")
    df = df[df["quantity"] > 0]
    bad = sorted(set(df["evidence"]) - set(EVIDENCE))
    if bad:
        raise ValueError(f"{source}: unknown evidence values {bad}")
    dup = df.duplicated(["precursor", "run"])
    if dup.any():
        first = df.loc[dup, ["precursor", "run"]].iloc[0].tolist()
        raise ValueError(f"{source}: {int(dup.sum())} duplicate (precursor, run) rows, e.g. {first}")
    df["residue_pair"] = [normalise_key(k) for k in df["residue_pair"]]
    return df.reset_index(drop=True)


def read(path):
    """Read and validate a long table (Parquet when the name ends in .parquet, else TSV)."""
    path = str(path)
    df = pd.read_parquet(path) if path.endswith(".parquet") else pd.read_csv(path, sep="\t")
    return validate(df, path)
