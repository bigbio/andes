"""summarize.py <raw run dir> <tag> <commit> <raw> <tier> <extra args> <perc cmd> [results.tsv]

Reads <dir>/seed<S>/eval.txt (eval_dia.py output) and the /usr/bin/time -v files of the
search (<dir>/search.time, minus the in-andes Percolator calls in <dir>/perc_calls/) and of the
per-seed rescoring (<dir>/seed<S>/perc.time). Prints one summary line and, with a results file,
appends it (writing the header first if the file is new).
"""
import os, re, subprocess, sys, statistics, datetime

d, tag, commit, raw, tier, extra, perc = sys.argv[1:8]
out = sys.argv[8] if len(sys.argv) > 8 else None
here = os.path.dirname(os.path.abspath(__file__))


def timev(paths):
    paths = [p for p in paths if os.path.exists(p)]
    if not paths:
        return (float('nan'),) * 3
    r = subprocess.run([sys.executable, '-I', os.path.join(here, 'timev.py'), *paths],
                       capture_output=True, text=True, check=True).stdout.split()
    return tuple(float(x) for x in r)


pat = re.compile(r'precursors (\d+) \(not in library (\d+)\) \| q<=0.01: N (\d+) E (\d+) FDP_c ([\d.]+)'
                 r' \| q<=0.05: N (\d+) E (\d+) FDP_c ([\d.]+) \| N@FDP_c 1%/2%: (\d+)/(\d+)')
seeds = sorted((int(s[4:]) for s in os.listdir(d) if s.startswith('seed') and
                os.path.exists(os.path.join(d, s, 'eval.txt'))),
               key=lambda s: [42, 1, 2, 3, 7].index(s) if s in (42, 1, 2, 3, 7) else 100 + s)
rows = {}
for s in seeds:
    m = pat.search(open(os.path.join(d, f'seed{s}', 'eval.txt')).read())
    if m:
        g = m.groups()
        rows[s] = dict(unmapped=int(g[1]), n01=int(g[2]), f01=float(g[4]), n05=int(g[5]),
                       f05=float(g[7]), a1=int(g[8]), a2=int(g[9]))
seeds = [s for s in seeds if s in rows]


def agg(key, fmt):
    v = [rows[s][key] for s in seeds]
    if not v:
        return 'NA', 'NA'
    per = ','.join(fmt.format(x) for x in v)
    return per, f'{fmt.format(statistics.mean(v))} ({fmt.format(min(v))}-{fmt.format(max(v))})'


calls = sorted(os.path.join(d, 'perc_calls', f) for f in os.listdir(os.path.join(d, 'perc_calls'))
               if f.endswith('.time')) if os.path.isdir(os.path.join(d, 'perc_calls')) else []
aw, ac, ar = timev([os.path.join(d, 'search.time')])
pw, pc, pr = timev(calls)
if not calls:
    pw = pc = 0.0
rw, rc, rr = timev([os.path.join(d, f'seed{s}', 'perc.time') for s in seeds])
cols = [('date', datetime.datetime.now().strftime('%Y-%m-%dT%H:%M')), ('tag', tag),
        ('commit', commit), ('raw', raw), ('tier', tier), ('seeds', ','.join(map(str, seeds)))]
for key, name, fmt in (('n01', 'N_q01', '{:.0f}'), ('f01', 'FDPc_q01', '{:.4f}'),
                       ('f05', 'FDPc_q05', '{:.4f}'), ('a1', 'N_at_FDPc1', '{:.0f}'),
                       ('a2', 'N_at_FDPc2', '{:.0f}')):
    per, mean = agg(key, fmt)
    cols += [(name + '_mean', mean), (name + '_seeds', per)]
cols += [('unmapped', ','.join(str(rows[s]['unmapped']) for s in seeds) or 'NA'),
         ('search_wall_s', f'{aw - pw:.0f}'), ('search_cpu_s', f'{ac - pc:.0f}'),
         ('search_rss_gb', f'{ar:.2f}'),
         ('andes_perc_calls', str(len(calls))), ('andes_perc_wall_s', f'{pw:.0f}'),
         ('andes_perc_cpu_s', f'{pc:.0f}'),
         ('rescore_wall_s', f'{rw:.0f}'), ('rescore_cpu_s', f'{rc:.0f}'),
         ('rescore_rss_gb', f'{rr:.2f}'), ('total_wall_s', f'{aw + (rw if rw == rw else 0):.0f}'),
         ('total_cpu_s', f'{ac + (rc if rc == rc else 0):.0f}'), ('extra_args', extra or '-'), ('perc', perc)]
line = '\t'.join(v for _, v in cols)
print(line)
if out:
    new = not os.path.exists(out)
    with open(out, 'a') as fh:
        if new:
            fh.write('\t'.join(k for k, _ in cols) + '\n')
        fh.write(line + '\n')
