import numpy as np, json
from halo import Cloud, simulate_halo, localize

rng = np.random.default_rng(101)
cl = [Cloud(100.0, 1.0, 0.1)]
d0 = simulate_halo(cl, 6*3600, 1300, 1.0, n_signal=20000, bkg_per_arcmin2=0,
                   fov_arcmin=20, source_xy=(0, 0), rng=rng)
s = d0['is_signal']
sig_r = np.hypot(d0['x'][s], d0['y'][s]).std()      # arcmin
out = {"sigma_r_arcsec": sig_r*60}
for N in [200]:
    K = 220
    xs, ys = [], []
    for _ in range(K):
        d = simulate_halo(cl, 6*3600, 1300, 1.0, n_signal=N,
                          bkg_per_arcmin2=0.07, fov_arcmin=20,
                          source_xy=(0, 0), rng=rng)
        xc, yc, _, _ = localize(d, half_width=5, n_grid=41)
        xs.append(xc*60); ys.append(yc*60)
    xs, ys = np.array(xs), np.array(ys)
    bias = np.hypot(xs.mean(), ys.mean())
    sd = 0.5*(xs.std()+ys.std())
    sig_pred = sig_r*60*np.sqrt(2.0/N)
    r68 = 1.5136*sig_pred
    cov = float(np.mean(np.hypot(xs, ys) < r68))
    out[str(N)] = dict(bias=float(bias), sd=float(sd), sig_pred=float(sig_pred),
                       ratio=float(sd/sig_pred), coverage68=cov, K=K)
    print(f"N={N} K={K}: bias={bias:.3f}\"  sd={sd:.3f}\"  pred_sigma={sig_pred:.3f}\"  "
          f"sd/pred={sd/sig_pred:.2f}  coverage(1sigma,68%)={cov:.2f}")
json.dump(out, open("../results/calib.json", "w"), indent=2)
