"""timev.py <time -v file>...: wall seconds, CPU seconds (user + system) and the largest peak RSS
in GB, summed (wall, CPU) or maxed (RSS) over the files, tab-separated."""
import sys

wall = cpu = 0.0
rss = 0
for path in sys.argv[1:]:
    for line in open(path):
        k, _, v = line.strip().rpartition(': ')
        if k.startswith('User time') or k.startswith('System time'):
            cpu += float(v)
        elif k.startswith('Elapsed (wall clock) time'):
            parts = [float(x) for x in v.split(':')]
            s = 0.0
            for x in parts:
                s = s * 60 + x
            wall += s
        elif k.startswith('Maximum resident set size'):
            rss = max(rss, int(v))
print(f'{wall:.1f}\t{cpu:.1f}\t{rss / 1048576:.2f}')
