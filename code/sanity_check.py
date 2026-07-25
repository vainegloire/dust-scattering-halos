import numpy as np
from halo import (theta_ring_arcmin, theta0_arcmin, Cloud,
                  simulate_halo, localize)

print("=== Ring geometry vs Nederlander & Paerels (2020) ===")
for t_hr, d in [(6.0, 100.0), (1.9, 100.0), (6.0, 50.0), (6.0, 200.0)]:
    th = theta_ring_arcmin(t_hr * 3600.0, d)
    print(f"  t={t_hr:>4} hr, d={d:>5.0f} pc  ->  ring radius = {th:6.2f} arcmin")
print("  (paper: 7 arcmin at 6 hr for d=100 pc; ~4 arcmin at 1.9 hr)")

print("\n=== theta0 (Gaussian cross-section width) ===")
for E, a in [(1.0, 0.1), (1.0, 0.2), (2.0, 0.1)]:
    print(f"  E={E} keV, a={a} um -> theta0 = {theta0_arcmin(E, a):.2f} arcmin")

print("\n=== Localization sanity: single-cloud, no background ===")
rng = np.random.default_rng(0)
cl = [Cloud(d_pc=100.0, tau=1.0, a_um=0.1)]
true_xy = (1.2, -0.8)
data = simulate_halo(cl, t_start_s=6*3600, t_exp_s=1300, E_keV=1.0,
                     n_signal=300, bkg_per_arcmin2=0.0, fov_arcmin=20.0,
                     source_xy=true_xy, rng=rng)
xc, yc, L, ext = localize(data, half_width=6, n_grid=61)
err = np.hypot(xc - true_xy[0], yc - true_xy[1]) * 60.0
print(f"  true center = {true_xy},  recovered = ({xc:.3f}, {yc:.3f}) arcmin")
print(f"  localization error = {err:.1f} arcsec  (with {int(data['is_signal'].sum())} halo photons)")

print("\n=== Multi-cloud nested rings ===")
clm = [Cloud(60, tau=1.0, a_um=0.1),
       Cloud(150, tau=1.2, a_um=0.1),
       Cloud(400, tau=1.5, a_um=0.1)]
tmid = 6*3600 + 650
for c in clm:
    print(f"  cloud d={c.d_pc:>4.0f} pc -> ring radius {theta_ring_arcmin(tmid, c.d_pc):5.2f} arcmin")
