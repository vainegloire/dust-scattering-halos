"""Localization-success fraction vs N for one background. Saves to detect.json.
Usage: python3 detect_mc.py <bkg> [K]
"""
import sys, json, os
import numpy as np
from halo import Cloud, simulate_halo, localize_selfcal

bkg = float(sys.argv[1])
K = int(sys.argv[2]) if len(sys.argv) > 2 else 20
FILE = "../results/detect.json"
cl = [Cloud(100.0, 1.0, 0.1)]
N_list = [20, 30, 45, 60, 120, 300]
thresh = 5.0
r = np.random.default_rng(int(1000*bkg) + 7)
frac = []
for n in N_list:
    ok = 0
    for _ in range(K):
        d = simulate_halo(cl, 6*3600, 1300, 1.0, n_signal=n,
                          bkg_per_arcmin2=bkg, fov_arcmin=20.0,
                          source_xy=(0.0, 0.0), rng=r)
        xc, yc, _, _ = localize_selfcal(d, half_width=5, n_grid=31)
        if np.hypot(xc, yc) * 60.0 < thresh:
            ok += 1
    frac.append(ok / K)
res = json.load(open(FILE)) if os.path.exists(FILE) else {"N_list": N_list, "thresh": thresh}
res[f"bkg_{bkg}"] = frac
json.dump(res, open(FILE, "w"), indent=2)
print(f"bkg={bkg} K={K}:", frac)
