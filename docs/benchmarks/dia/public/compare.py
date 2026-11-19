"""compare.py <base tag> <variant tag> [results.tsv]: paired-by-seed comparison of two runs in
RESULTS.tsv (the last line of each tag and raw): per raw, the mean and per-seed difference of
N@FDP_c 1 % and 2 % and N at q <= 0.01, the number of seeds where the variant is higher, and
FDP_c at q <= 0.01 / 0.05 of both."""
import os, sys, csv

base, var = sys.argv[1:3]
path = sys.argv[3] if len(sys.argv) > 3 else os.environ.get('BENCH_ROOT', '/home/user/dia-bench') + '/runs/RESULTS.tsv'
rows = {}
with open(path) as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        rows[(r['tag'], r['raw'])] = r


def per(r, k):
    return dict(zip(r['seeds'].split(','), (float(x) for x in r[k + '_seeds'].split(','))))


for raw in sorted({k[1] for k in rows}):
    a, b = rows.get((base, raw)), rows.get((var, raw))
    if not a or not b:
        continue
    out = [f'raw {raw}']
    for k in ('N_at_FDPc1', 'N_at_FDPc2', 'N_q01'):
        pa, pb = per(a, k), per(b, k)
        s = [x for x in pa if x in pb]
        d = [pb[x] - pa[x] for x in s]
        up = sum(x > 0 for x in d)
        out.append(f'{k} {sum(pa[x] for x in s) / len(s):.0f} -> {sum(pb[x] for x in s) / len(s):.0f} '
                   f'({sum(d) / len(d):+.0f}; per seed {",".join(f"{x:+.0f}" for x in d)}; higher {up}/{len(d)})')
    out.append(f'FDP_c q01 {a["FDPc_q01_mean"]} -> {b["FDPc_q01_mean"]}; q05 {a["FDPc_q05_mean"]} -> {b["FDPc_q05_mean"]}')
    print(' | '.join(out))
