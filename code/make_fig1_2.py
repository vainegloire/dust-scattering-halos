"""Generate Figure 1 (single-screen halo) and Figure 2 (multi-cloud nested
rings), each with the localization log-likelihood map."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from halo import Cloud, simulate_halo, localize, theta_ring_arcmin

plt.rcParams.update({"font.size": 11, "figure.dpi": 130})
rng = np.random.default_rng(42)

# ---- Figure 1: single screen at 100 pc -------------------------------
cl1 = [Cloud(100.0, 1.0, 0.1)]
d1 = simulate_halo(cl1, 6*3600, 1300, 1.0, n_signal=60, bkg_per_arcmin2=0.07,
                   fov_arcmin=20.0, source_xy=(0.0, 0.0), rng=rng)
xc, yc, L, ext = localize(d1, half_width=6, n_grid=61)
fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.6))
s = d1['is_signal']
ax[0].scatter(d1['x'][s], d1['y'][s], s=9, c='k', label='scattered')
ax[0].scatter(d1['x'][~s], d1['y'][~s], s=9, c='0.7', marker='x', label='background')
ax[0].plot(0, 0, 'r+', ms=13, mew=2, label='true source')
ax[0].plot(xc, yc, 'bx', ms=10, mew=2, label='recovered')
ax[0].set_aspect('equal'); ax[0].set_xlim(-15, 15); ax[0].set_ylim(-15, 15)
ax[0].set_xlabel('offset (arcmin)'); ax[0].set_ylabel('offset (arcmin)')
ax[0].set_title('Single screen at 100 pc, 6 hr delay')
ax[0].legend(loc='upper right', fontsize=8, framealpha=0.9)
im = ax[1].imshow(L, origin='lower', extent=ext, cmap='viridis', aspect='equal')
ax[1].plot(0, 0, 'r+', ms=13, mew=2); ax[1].plot(xc, yc, 'wx', ms=9, mew=2)
ax[1].set_xlabel('center offset (arcmin)'); ax[1].set_ylabel('center offset (arcmin)')
ax[1].set_title('Localization log-likelihood')
fig.colorbar(im, ax=ax[1], fraction=0.046, pad=0.04, label='log L')
fig.tight_layout(); fig.savefig("../figures/fig1_single_cloud.png"); plt.close(fig)

# ---- Figure 2: three screens -----------------------------------------
clm = [Cloud(60, 1.0, 0.1), Cloud(150, 1.3, 0.1), Cloud(400, 1.6, 0.1)]
d2 = simulate_halo(clm, 6*3600, 1300, 1.0, n_signal=600, bkg_per_arcmin2=0.07,
                   fov_arcmin=20.0, source_xy=(0.0, 0.0), rng=rng)
xc2, yc2, L2, ext2 = localize(d2, half_width=6, n_grid=61)
fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.6))
s = d2['is_signal']
ax[0].scatter(d2['x'][s], d2['y'][s], s=7, c='k')
ax[0].scatter(d2['x'][~s], d2['y'][~s], s=7, c='0.75', marker='x')
tmid = 6*3600 + 650
for c in clm:
    rr = theta_ring_arcmin(tmid, c.d_pc); th = np.linspace(0, 2*np.pi, 200)
    ax[0].plot(rr*np.cos(th), rr*np.sin(th), 'r--', lw=0.7, alpha=0.6)
ax[0].plot(0, 0, 'r+', ms=12, mew=2)
ax[0].set_aspect('equal'); ax[0].set_xlim(-15, 15); ax[0].set_ylim(-15, 15)
ax[0].set_xlabel('offset (arcmin)'); ax[0].set_ylabel('offset (arcmin)')
ax[0].set_title('Three screens: 60, 150, 400 pc')
im = ax[1].imshow(L2, origin='lower', extent=ext2, cmap='viridis', aspect='equal')
ax[1].plot(0, 0, 'r+', ms=12, mew=2); ax[1].plot(xc2, yc2, 'wx', ms=9, mew=2)
ax[1].set_xlabel('center offset (arcmin)'); ax[1].set_ylabel('center offset (arcmin)')
ax[1].set_title('Localization log-likelihood')
fig.colorbar(im, ax=ax[1], fraction=0.046, pad=0.04, label='log L')
fig.tight_layout(); fig.savefig("../figures/fig2_multi_cloud.png"); plt.close(fig)
print("Saved fig1, fig2 to ../figures/")
