#!/usr/bin/env python3
"""Independence gate: fail if any shipped model's tables are the MS-GF+ seed.

WHY THIS EXISTS
---------------
The model-training pipeline pre-seeds each slug store with the 39 MS-GF+ seed
models (version 7061) and overlays one slug per training run. If a training run
errors (e.g. 0 confident labels) the slug keeps its *seed* tables, and the
assembler then re-stamps every row to an "own" version (1/2). The version field
therefore CANNOT be trusted as proof of independence — a slug can carry the
original MS-GF+ table values under an "own" version label.

This gate is version-stamp-independent: it hashes each model's actual table
VALUES and fails if any model in the candidate bundle is byte-identical to any
model in the MS-GF+ seed. Run it at assembly/release time against the seed
(which lives on the training host — it is deliberately NOT shipped in this repo).

Note: the geometry (num_segments / partition layout) is a separate axis and is
NOT checked here (see the partition-structure re-derivation track / Proposal B).

USAGE
-----
    check_models_independence.py --models own_models.parquet --seed seed-models.parquet
    check_models_independence.py --models own_models.parquet --seed seed.parquet \
        --allow hcd_qexactive_tryp_phosphorylation,hcd_highres_tryp_phosphorylation

Exit code 0 = clean (no model matches the seed). Exit code 1 = leak detected
(one or more models are byte-identical to an MS-GF+ seed model).
"""
import argparse
import collections
import hashlib
import sys

# Columns that key a table record; the `values` list is hashed per (model, key).
_KEY_COLS = [
    "table_kind", "ion_kind", "ion_charge", "ion_offset_bits",
    "part_charge", "part_mass_bits", "part_seg",
]


def table_fingerprints(path):
    """Return {model_id: hex16} — a canonical hash of every table's `values`.

    The hash is independent of row order and of the `version` stamp, so two
    stores with identical table content hash identically regardless of how
    their version field was set.
    """
    import pyarrow.parquet as pq
    import pyarrow.compute as pc

    t = pq.read_table(path)
    tt = t.filter(pc.equal(t["record_kind"], "table")).to_pydict()
    n = len(tt["model_id"])
    rows = collections.defaultdict(list)
    for i in range(n):
        key = tuple(str(tt[c][i]) for c in _KEY_COLS if c in tt)
        vals = tuple(round(float(x), 9) for x in (tt["values"][i] or []))
        rows[tt["model_id"][i]].append((key, vals))
    out = {}
    for mid, recs in rows.items():
        recs.sort(key=lambda r: r[0])
        out[mid] = hashlib.sha256(repr(recs).encode()).hexdigest()[:16]
    return out


def find_leaks(models_path, seed_path):
    """Return list of (model_id, matched_seed_id, hash) for seed-identical models."""
    model_fp = table_fingerprints(models_path)
    seed_fp = table_fingerprints(seed_path)
    # Reverse index: hash -> seed model_id(s) that produce it.
    seed_by_hash = collections.defaultdict(list)
    for sid, h in seed_fp.items():
        seed_by_hash[h].append(sid)
    leaks = []
    for mid, h in sorted(model_fp.items()):
        if h in seed_by_hash:
            leaks.append((mid, ",".join(sorted(seed_by_hash[h])), h))
    return leaks


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", required=True, help="candidate bundle .parquet to check")
    ap.add_argument("--seed", required=True, help="MS-GF+ seed .parquet to check against")
    ap.add_argument("--allow", default="",
                    help="comma-separated model_ids permitted to match the seed "
                         "(temporary waiver; each one is a known MS-GF+ leak still to fix)")
    args = ap.parse_args(argv)

    allow = {s.strip() for s in args.allow.split(",") if s.strip()}
    leaks = find_leaks(args.models, args.seed)

    blocking = [l for l in leaks if l[0] not in allow]
    waived = [l for l in leaks if l[0] in allow]

    for mid, sid, h in waived:
        print(f"WAIVED  {mid:40s} == MS-GF+ seed {sid} (hash {h})")
    for mid, sid, h in blocking:
        print(f"LEAK    {mid:40s} == MS-GF+ seed {sid} (hash {h})")

    if blocking:
        print(f"\nFAIL: {len(blocking)} model(s) are byte-identical to the MS-GF+ seed "
              f"(version stamp is not independence). Retrain or drop them before shipping.",
              file=sys.stderr)
        return 1
    print(f"\nOK: no model matches the MS-GF+ seed "
          f"({len(waived)} waived).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
