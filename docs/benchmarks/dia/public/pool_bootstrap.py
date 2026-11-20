"""pool_bootstrap.py <name> <npz>...: sum of the paired differences over slices; the slices' bootstrap
replicates are independent, so replicate i of the sum is the sum of every slice's replicate i."""
import sys

import numpy as np

name, paths = sys.argv[1], sys.argv[2:]
z = [np.load(p) for p in paths]
point = sum(x['point'] for x in z)
reps = sum(x['reps'] for x in z)
for j, lab in enumerate(('N@FDP_c 1%', 'N@FDP_c 2%')):
    lo, hi = np.quantile(reps[:, j], [0.025, 0.975])
    print(f'{name} pooled over {len(paths)} slices {lab}: sum {point[j]:+.0f}, boot sd '
          f'{reps[:, j].std():.0f}, 95% [{lo:+.0f}, {hi:+.0f}]')
