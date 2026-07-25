"""Ensemble of random realistic sightlines (disk-and-log-normal model of
Section 2.3, as in fig7): localize each halo with the self-calibrating
estimator and accumulate the errors, to turn the single-sightline example of
Section 3.7 into a statistic.

Usage:
  python3 ensemble_sightlines.py <n_chunk> <seed>   (accumulates)
  python3 ensemble_sightlines.py report
"""
import sys, json, os
import numpy as np
from halo import simulate_halo, sample_disk_sightline, localize_selfcal

FILE = "../results/sightline_ensemble.json"

if sys.argv[1] == "report":
    acc = json.load(open(FILE))
    e = np.array(acc["err_arcsec"])
    nc = np.array(acc["n_clouds"])
    print(f"n_sightlines={e.size}  clouds/sightline: {nc.min()}-{nc.max()} "
          f"(median {int(np.median(nc))})")
    print(f"median={np.median(e):.2f}\"  p90={np.percentile(e, 90):.2f}\"  "
          f"max={e.max():.2f}\"  fraction<5\"={np.mean(e < 5.0):.2f}")
else:
    K, seed = int(sys.argv[1]), int(sys.argv[2])
    acc = (json.load(open(FILE)) if os.path.exists(FILE)
           else {"err_arcsec": [], "n_clouds": []})
    rng = np.random.default_rng(seed)
    for _ in range(K):
        clouds = sample_disk_sightline(rng, n_clouds_mean=12, d_scale_pc=280,
                                       d_max_pc=1200, a_um=0.1)
        clouds = [c for c in clouds if c.d_pc > 40]
        data = simulate_halo(clouds, 6 * 3600, 1300, E_keV=1.0, n_signal=900,
                             bkg_per_arcmin2=0.07, fov_arcmin=20.0,
                             grain_spread=(0.05, 0.25), source_xy=(0.0, 0.0),
                             rng=rng)
        xc, yc, _, _ = localize_selfcal(data, half_width=5, n_grid=41)
        acc["err_arcsec"].append(float(np.hypot(xc, yc) * 60.0))
        acc["n_clouds"].append(len(clouds))
    json.dump(acc, open(FILE, "w"))
    print(f"total {len(acc['err_arcsec'])} sightlines")
