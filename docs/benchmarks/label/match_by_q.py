import sys, pyarrow.parquet as pq, pandas as pd, numpy as np
t = pq.read_table(sys.argv[1]).to_pandas()
for run, g in t.groupby("filename"):
    L = g[g.channel == 1].set_index(["peptide", "charge"]); H = g[g.channel == 2].set_index(["peptide", "charge"])
    j = L[["status", "channel_q_value", "quantity", "transfer"]].join(H[["status", "channel_q_value", "quantity"]], lsuffix="_l", rsuffix="_h")
    j = j[~j.transfer]
    j["r"] = np.log2(j.quantity_h / j.quantity_l)
    j = j.replace([np.inf, -np.inf], np.nan)
    both = j[(j.status_l == "identified") & (j.status_h == "identified")]
    print(run, "both identified", len(both), round(both.r.median(), 3), round(both.r.quantile(.75) - both.r.quantile(.25), 3))
    one = j[(j.status_l == "identified") ^ (j.status_h == "identified")].copy()
    one["q"] = np.where(one.status_l == "identified", one.channel_q_value_h, one.channel_q_value_l)
    one["status"] = np.where(one.status_l == "identified", one.status_h, one.status_l)
    one["qbin"] = pd.cut(one.q, [0, 0.01, 0.02, 0.05, 0.1, 0.2, 1.0], include_lowest=True)
    x = one.dropna(subset=["r"]).groupby("qbin", observed=False).r
    print(pd.DataFrame({"n": x.size(), "median": x.median().round(3), "iqr": (x.quantile(.75) - x.quantile(.25)).round(3)}))
    print(one.status.value_counts().to_dict())
