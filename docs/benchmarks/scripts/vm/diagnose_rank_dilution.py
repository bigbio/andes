#!/usr/bin/env python3
"""Compare rank_dist_table sharpness: trained cid_msnet_tryp vs seed cid_lowres_tryp."""
import sys
from pathlib import Path

import pyarrow.parquet as pq


def load_param_tables(store: Path, model_id: str):
    t = pq.read_table(store, filters=[("model_id", "=", model_id)])
    if t.num_rows == 0:
        raise SystemExit(f"model {model_id!r} not in {store}")
    row = t.to_pydict()
    # rank_dist_table is nested; use first partition slice mean max prob as sharpness proxy
    rdt = row["rank_dist_table"][0]
    if not rdt:
        return 0.0, 0
    maxes = []
    for ion_vec in rdt:
        if not ion_vec:
            continue
        for dist in ion_vec:
            if dist:
                maxes.append(max(dist))
    if not maxes:
        return 0.0, 0
    return sum(maxes) / len(maxes), len(maxes)


def main() -> None:
    store = Path(sys.argv[1] if len(sys.argv) > 1 else "/srv/data/msnet/models_cid_msnet.parquet")
    seed = Path(sys.argv[2] if len(sys.argv) > 2 else str(
        Path(__file__).resolve().parents[2] / "resources/ionstat/models.parquet"
    ))
    trained_id = sys.argv[3] if len(sys.argv) > 3 else "cid_msnet_tryp"
    seed_id = sys.argv[4] if len(sys.argv) > 4 else "cid_lowres_tryp"

    t_mean, t_n = load_param_tables(store, trained_id)
    s_mean, s_n = load_param_tables(seed, seed_id)
    print(f"store={store}")
    print(f"{trained_id}: mean max rank-bin prob = {t_mean:.4f} ({t_n} dists)")
    print(f"{seed_id}: mean max rank-bin prob = {s_mean:.4f} ({s_n} dists)")
    if s_mean > 0:
        print(f"ratio trained/seed = {t_mean / s_mean:.3f}  (<1 => softer/diluted)")


if __name__ == "__main__":
    main()
