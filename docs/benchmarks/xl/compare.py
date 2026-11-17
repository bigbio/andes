"""compare.py tagA tagB ...: per dataset, true CSMs at andes q<=1% (all classes), actual false share,
worst per-class false share (classes with >=50 accepted), and true CSMs at 1% actual FDP."""
import json, os, sys
tags = sys.argv[1:]
S = {t: json.load(open(os.path.join(os.environ.get('XLBENCH_DATA', '.'), 'results', t, 'summary.json'))) for t in tags}
names = [n for n in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'datasets.json'))) if all(n['name'] in S[t] and 'all' in S[t][n['name']] for t in tags)]
print(f"{'dataset':18}" + ''.join(f"| {t[:22]:>22} " for t in tags))
print(f"{'':18}" + ''.join(f"| {'true@1%  fdp  wcls  @act':>22} " for t in tags))
tot = {t: [0, 0] for t in tags}
for n in names:
    row = f"{n['name']:18}"
    for t in tags:
        r = S[t][n['name']]; a = r['all']
        worst = max([r[c]['fdp'] for c in ('intra', 'homotypic', 'inter') if c in r and r[c]['csm'] >= 50] or [0])
        true = a['csm'] - a['false']
        tot[t][0] += true; tot[t][1] += a['true_at_1pct_actual']
        row += f"| {true:6d} {a['fdp']*100:4.1f} {worst*100:5.1f} {a['true_at_1pct_actual']:6d} "
    print(row)
print(f"{'TOTAL':18}" + ''.join(f"| {tot[t][0]:6d}            {tot[t][1]:6d} " for t in tags))
