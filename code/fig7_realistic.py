import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from halo import (simulate_halo, sample_disk_sightline, localize_selfcal,
                  recover_ring_radii, theta_ring_arcmin)

plt.rcParams.update({"font.size": 11, "figure.dpi": 130})
rng = np.random.default_rng(15)

# a realistic sightline: ~12 clouds, exponential distances, lognormal columns
clouds = sample_disk_sightline(rng, n_clouds_mean=12, d_scale_pc=280,
                               d_max_pc=1200, a_um=0.1)
clouds = [c for c in clouds if c.d_pc > 40]          # keep rings within FoV
t0, texp = 6*3600, 1300
t_obs = t0 + texp/2
data = simulate_halo(clouds, t0, texp, E_keV=1.0, n_signal=900,
                     bkg_per_arcmin2=0.07, fov_arcmin=20.0,
                     grain_spread=(0.05, 0.25), source_xy=(0.0, 0.0), rng=rng)

xc, yc, _, _ = localize_selfcal(data, half_width=6, n_grid=51)
err = np.hypot(xc, yc) * 60.0
peaks, centers, dens = recover_ring_radii(data, xc, yc, t_obs,
                                          min_prominence=0.20)

true = sorted([(theta_ring_arcmin(t_obs, c.d_pc), c.d_pc, c.tau)
               for c in clouds if theta_ring_arcmin(t_obs, c.d_pc) < 20],
              key=lambda z: z[0])

fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
s = data['is_signal']
ax[0].scatter(data['x'][s], data['y'][s], s=6, c='k')
ax[0].scatter(data['x'][~s], data['y'][~s], s=6, c='0.75', marker='x')
ax[0].plot(0, 0, 'r+', ms=12, mew=2)
ax[0].plot(xc, yc, 'bx', ms=9, mew=2)
ax[0].set_aspect('equal'); ax[0].set_xlim(-18, 18); ax[0].set_ylim(-18, 18)
ax[0].set_xlabel('offset (arcmin)'); ax[0].set_ylabel('offset (arcmin)')
ax[0].set_title(f'Realistic sightline: {len(clouds)} clouds\ncenter recovered to {err:.1f}"')

ax[1].plot(centers, dens/dens.max(), color='C0')
for th, dpc, tau in true:
    ax[1].axvline(th, color='0.6', ls=':', lw=0.9)
for th, dpc in peaks:
    ax[1].axvline(th, color='C3', ls='--', lw=1.1)
ax[1].set_xlim(0, 16)
ax[1].set_xlabel('ring radius (arcmin)')
ax[1].set_ylabel('normalized radial density')
ax[1].set_title('Radial profile (dotted = true clouds,\ndashed = recovered rings)')
fig.tight_layout(); fig.savefig("../figures/fig7_realistic_sightline.png"); plt.close(fig)

print(f"n_clouds within FoV: {len(true)}, total clouds: {len(clouds)}")
print(f"self-cal center error: {err:.2f} arcsec")
print("true clouds (theta_arcmin, d_pc, tau):")
for th, dpc, tau in true:
    print(f"   {th:5.2f}'  {dpc:6.0f} pc  tau={tau:.2f}")
print("recovered rings (theta_arcmin, d_pc):")
for th, dpc in peaks:
    print(f"   {th:5.2f}'  {dpc:6.0f} pc")
