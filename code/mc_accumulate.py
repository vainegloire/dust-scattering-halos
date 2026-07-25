"""Accumulate Monte-Carlo localization trials in chunks (so long runs can be
split across short executions) and merge them into mc_results.json.

Usage:
  python3 mc_accumulate.py run <ideal|bkg|multi> <K_chunk> <seed>
  python3 mc_accumulate.py finalize
"""
import sys, json, os
import numpy as np
from halo import Cloud, simulate_halo, localize

cl1 = [Cloud(100.0, 1.0, 0.1)]
clm = [Cloud(60, 1.0, 0.1), Cloud(150, 1.3, 0.1), Cloud(400, 1.6, 0.1)]
n_list = [30, 60, 120, 250, 500, 1000]
BKG = {"ideal": 0.0, "bkg": 0.07, "multi": 0.07}
CLOUDS = {"ideal": cl1, "bkg": cl1, "multi": clm}


def one_error(clouds, n, bkg, r):
    d = simulate_halo(clouds, 6*3600, 1300, 1.0, n_signal=n,
                      bkg_per_arcmin2=bkg, fov_arcmin=20.0,
                      source_xy=(0.0, 0.0), rng=r)
    bfrac = min(0.6, (bkg*1600) / (bkg*1600 + n))
    xc, yc, _, _ = localize(d, half_width=5, n_grid=41, bkg_frac_guess=bfrac)
    return float(np.hypot(xc, yc) * 60.0)


if sys.argv[1] == "run":
    cfg, K, seed = sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    acc_file = f"../results/mc_acc_{cfg}.json"
    acc = (json.load(open(acc_file)) if os.path.exists(acc_file)
           else {str(n): [] for n in n_list})
    r = np.random.default_rng(seed)
    for n in n_list:
        acc[str(n)] += [one_error(CLOUDS[cfg], n, BKG[cfg], r)
                        for _ in range(K)]
    json.dump(acc, open(acc_file, "w"))
    print(cfg, {n: len(v) for n, v in acc.items()})

elif sys.argv[1] == "finalize":
    res = json.load(open("../results/mc_results.json"))
    for cfg in ("ideal", "bkg", "multi"):
        acc = json.load(open(f"../results/mc_acc_{cfg}.json"))
        res[cfg] = {n: dict(rms=float(np.sqrt(np.mean(np.array(e) ** 2))),
                            median=float(np.median(e)), K=len(e))
                    for n, e in acc.items()}
    json.dump(res, open("../results/mc_results.json", "w"), indent=2)
    for cfg in ("ideal", "bkg", "multi"):
        print(cfg, {n: round(v["rms"], 2) for n, v in res[cfg].items()})
