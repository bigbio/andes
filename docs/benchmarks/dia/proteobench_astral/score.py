"""score.py <custom.tsv> <label>: score a Custom table on the ProteoBench DIA Astral module with the installed
proteobench package (on PYTHONPATH); print metrics at min runs 1..6 and save them as <tsv>.metrics.json."""
import json, sys
import proteobench
from proteobench.modules.quant.quant_lfq_ion_DIA_Astral import DIAQuantIonModuleAstral
path, label = sys.argv[1:3]
m = DIAQuantIonModuleAstral(token="")
ui = dict(software_name="andes", software_version="dev", search_engine="andes", search_engine_version="dev",
          ident_fdr_psm=0.01, ident_fdr_peptide=None, ident_fdr_protein=None, enable_match_between_runs=True,
          precursor_mass_tolerance="learned", fragment_mass_tolerance="7 ppm", enzyme="Trypsin/P",
          allowed_miscleavages=1, min_peptide_length=6, max_peptide_length=30)
inter, alldp, inp = m.benchmarking(path, "Custom", ui, None)
cur = alldp[alldp["old_new"] == "new"].iloc[0] if "old_new" in alldp else alldp.iloc[-1]
res = cur["results"]
print(f"== {label}  proteobench {getattr(proteobench, '__version__', '?')}  rows in {len(inp)}  intermediate {len(inter)}")
for k in sorted(res, key=int):
    r = res[k]
    print(f"  min_runs {k}: nr_prec {r['nr_feature']}  median|eps| {r['median_abs_epsilon_global']:.4f}  "
          f"HUMAN {r['median_abs_epsilon_HUMAN']:.4f}  YEAST {r['median_abs_epsilon_YEAST']:.4f}  ECOLI {r['median_abs_epsilon_ECOLI']:.4f}  "
          f"CV_median {r['CV_median']:.4f}  log2 H/Y/E {r['median_log2_empirical_HUMAN']:.3f}/{r['median_log2_empirical_YEAST']:.3f}/{r['median_log2_empirical_ECOLI']:.3f}")
json.dump({str(k): {kk: (float(v) if isinstance(v, (int, float)) else str(v)) for kk, v in res[k].items()} for k in res},
          open(path + ".metrics.json", "w"), indent=1)
