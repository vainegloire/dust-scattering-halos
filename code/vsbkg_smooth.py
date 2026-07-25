"""Accumulate squared localization errors vs background to smooth Fig 4.
Usage: python3 vsbkg_smooth.py <seed> [K]
Accumulates into vsbkg_acc.json (sum of err^2 and count per background).
"""
import sys, json, os
import numpy as np
from halo import Cloud, simulate_halo, localize

seed = int(sys.argv[1])
K = int(sys.argv[2]) if len(sys.argv) > 2 else 80
FILE = "../results/vsbkg_acc.json"
cl = [Cloud(100.0, 1.0, 0.1)]
bkg_list = [0.0, 0.03, 0.07, 0.15, 0.3, 0.6]

acc = json.load(open(FILE)) if os.path.exists(FILE) else \
    {"bkg": bkg_list, "sumsq": [0.0]*len(bkg_list), "count": [0]*len(bkg_list)}
r = np.random.default_rng(seed)
for i, b in enumerate(bkg_list):
    for _ in range(K):
        d = simulate_halo(cl, 6*3600, 1300, 1.0, n_signal=300,
                          bkg_per_arcmin2=b, fov_arcmin=20.0,
                          source_xy=(0.0, 0.0), rng=r)
        bfrac = min(0.6, (b*1600)/(b*1600+300))
        xc, yc, _, _ = localize(d, half_width=5, n_grid=41, bkg_frac_guess=bfrac)
        acc["sumsq"][i] += float((np.hypot(xc, yc)*60.0)**2)
        acc["count"][i] += 1
json.dump(acc, open(FILE, "w"), indent=2)
rms = [np.sqrt(s/c) for s, c in zip(acc["sumsq"], acc["count"])]
print("counts:", acc["count"])
print("rms:", [round(x, 2) for x in rms])
