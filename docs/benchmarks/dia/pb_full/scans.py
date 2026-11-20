"""Scan table of an mzML: MS2 scan number, retention time (min) and isolation window.

scans.py <mzML> <out.parquet>
"""
import re
import sys

import numpy as np
import pandas as pd

scan_re = re.compile(rb'scan=(\d+)')
val = re.compile(rb'value="([^"]*)"')
rows, cur = [], None
with open(sys.argv[1], 'rb') as fh:
    for line in fh:
        s = line.lstrip()
        if s.startswith(b'<spectrum '):
            m = scan_re.search(s)
            cur = {'scan': int(m.group(1)) if m else -1, 'ms': 0, 'rt': np.nan, 't': np.nan, 'lo': 0.0, 'hi': 0.0}
        elif cur is None:
            continue
        elif b'MS:1000511' in s:
            cur['ms'] = int(val.search(s).group(1))
        elif b'MS:1000016' in s:
            v = float(val.search(s).group(1))
            cur['rt'] = v if b'minute' in s else v / 60.0
        elif b'MS:1000827' in s:
            cur['t'] = float(val.search(s).group(1))
        elif b'MS:1000828' in s:
            cur['lo'] = float(val.search(s).group(1))
        elif b'MS:1000829' in s:
            cur['hi'] = float(val.search(s).group(1))
        elif s.startswith(b'</spectrum>'):
            if cur['ms'] == 2:
                rows.append(cur)
            cur = None
d = pd.DataFrame(rows)
d['wlo'], d['whi'] = d.t - d.lo, d.t + d.hi
d[['scan', 'rt', 't', 'wlo', 'whi']].to_parquet(sys.argv[2])
print(len(d), 'MS2 scans')
