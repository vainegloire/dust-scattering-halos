"""Monte-Carlo the self-calibrating (distance-agnostic) localizer.

Stores the RAW per-trial errors, as mc_accumulate.py does for Table 1, so
that derived quantities -- in particular the catastrophic-failure rate quoted
in Section 3.5 -- are auditable from the repository rather than only from the
summary statistics.

Usage:  python3 selfcal_mc.py run <single|multi> <K_chunk> <seed>
        python3 selfcal_mc.py finalize
"""
import sys
import json
import os
import time

import numpy as np
from halo import Cloud, fiducial_clouds, simulate_halo, localize_selfcal

N_LIST = [30, 60, 120, 250, 500, 1000]
CLOUDS = {"single": [Cloud(100.0, 1.0, 0.1)], "multi": fiducial_clouds()}
FAIL_ARCSEC = 30.0          # "catastrophic": locked onto a spurious centre


def _acc_path(cfg):
    return f"../results/selfcal_acc_{cfg}.json"


def one_error(clouds, n, rng):
    d = simulate_halo(clouds, 6 * 3600, 1300, 1.0, n_signal=n,
                      bkg_per_arcmin2=0.07, fov_arcmin=20.0,
                      source_xy=(0.0, 0.0), rng=rng)
    xc, yc, _, _ = localize_selfcal(d, half_width=5, n_grid=41)
    return float(np.hypot(xc, yc) * 60.0)


if __name__ == "__main__":
    mode = sys.argv[1]

    if mode == "run":
        cfg, K, seed = sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
        path = _acc_path(cfg)
        acc = (json.load(open(path)) if os.path.exists(path)
               else {str(n): [] for n in N_LIST})
        rng = np.random.default_rng(seed)
        t0 = time.time()
        for n in N_LIST:
            acc[str(n)] += [one_error(CLOUDS[cfg], n, rng) for _ in range(K)]
        json.dump(acc, open(path, "w"))
        print(f"[{cfg}] {time.time()-t0:.0f}s  "
              f"{ {n: len(v) for n, v in acc.items()} }")

    elif mode == "finalize":
        out = {}
        for cfg in ("single", "multi"):
            acc = json.load(open(_acc_path(cfg)))
            out[cfg] = {}
            for n, e in acc.items():
                e = np.array(e)
                out[cfg][n] = dict(
                    rms=float(np.sqrt(np.mean(e ** 2))),
                    median=float(np.median(e)),
                    fail_frac=float(np.mean(e > FAIL_ARCSEC)),
                    K=int(e.size))
        json.dump(out, open("../results/selfcal_results.json", "w"), indent=2)
        for cfg in out:
            print(cfg, {n: (round(v["median"], 2), round(v["fail_frac"], 3))
                        for n, v in out[cfg].items()})
