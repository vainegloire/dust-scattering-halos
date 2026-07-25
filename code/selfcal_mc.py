"""Monte-Carlo the self-calibrating (distance-agnostic) localizer.
Usage: python3 selfcal_mc.py <single|multi> [K]
"""
import sys, json, os, time
import numpy as np
from halo import Cloud, simulate_halo, localize_selfcal

CFG = sys.argv[1]
K = int(sys.argv[2]) if len(sys.argv) > 2 else 25
FILE = "../results/selfcal_results.json"

cl1 = [Cloud(100.0, 1.0, 0.1)]
clm = [Cloud(60, 1.2, 0.1), Cloud(150, 1.3, 0.1), Cloud(400, 1.5, 0.1)]
n_list = [30, 60, 120, 250, 500, 1000]
clouds = cl1 if CFG == "single" else clm


def scan(seed):
    r = np.random.default_rng(seed)
    out = {}
    for n in n_list:
        errs = []
        for _ in range(K):
            d = simulate_halo(clouds, 6*3600, 1300, 1.0, n_signal=n,
                              bkg_per_arcmin2=0.07, fov_arcmin=20.0,
                              source_xy=(0.0, 0.0), rng=r)
            xc, yc, _, _ = localize_selfcal(d, half_width=5, n_grid=41)
            errs.append(np.hypot(xc, yc) * 60.0)
        errs = np.array(errs)
        out[str(n)] = dict(rms=float(np.sqrt(np.mean(errs**2))),
                           median=float(np.median(errs)))
    return out


t0 = time.time()
res = json.load(open(FILE)) if os.path.exists(FILE) else {}
res[CFG] = scan(11 if CFG == "single" else 12)
json.dump(res, open(FILE, "w"), indent=2)
print(f"[selfcal:{CFG}] {time.time()-t0:.1f}s  K={K}")
print(json.dumps(res[CFG], indent=2))
