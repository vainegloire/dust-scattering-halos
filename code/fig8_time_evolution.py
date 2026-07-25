import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from halo import (Cloud, simulate_halo, localize_selfcal, recover_ring_radii,
                  theta_ring_arcmin)

plt.rcParams.update({"font.size": 11, "figure.dpi": 130})
rng = np.random.default_rng(5)

fig, ax = plt.subplots(1, 2, figsize=(11, 4.7))

# --- Panel A: one screen imaged at two epochs -> the ring expands ------
cl = [Cloud(100.0, 1.0, 0.1)]
texp = 1300
for (t_hr, col, lab) in [(3.0, 'C0', '3 hr'), (9.0, 'C3', '9 hr')]:
    d = simulate_halo(cl, t_hr*3600, texp, 1.0, n_signal=120,
                      bkg_per_arcmin2=0.0, fov_arcmin=20.0,
                      source_xy=(0.0, 0.0), rng=rng)
    s = d['is_signal']
    ax[0].scatter(d['x'][s], d['y'][s], s=10, c=col, alpha=0.8,
                  label=f'{lab}  ({theta_ring_arcmin(t_hr*3600,100):.1f}$\'$)')
# static background points (same sky positions at both epochs)
bx = rng.uniform(-18, 18, 40); by = rng.uniform(-18, 18, 40)
ax[0].scatter(bx, by, s=12, c='0.6', marker='x', label='background (static)')
ax[0].plot(0, 0, 'k+', ms=12, mew=2)
ax[0].set_aspect('equal'); ax[0].set_xlim(-16, 16); ax[0].set_ylim(-16, 16)
ax[0].set_xlabel('offset (arcmin)'); ax[0].set_ylabel('offset (arcmin)')
ax[0].set_title('Single screen (100 pc) at two epochs:\nthe ring expands as $\\sqrt{t}$')
ax[0].legend(loc='upper right', fontsize=8)

# --- Panel B: ring radius vs time for three screens -------------------
clm = [Cloud(60, 1.4, 0.1), Cloud(150, 1.3, 0.1), Cloud(400, 1.5, 0.1)]
epochs_hr = [2, 3, 4, 6, 9, 12]
meas = {c.d_pc: [] for c in clm}
meas_t = {c.d_pc: [] for c in clm}
for th_hr in epochs_hr:
    t0 = th_hr*3600
    d = simulate_halo(clm, t0, texp, 1.0, n_signal=800,
                      bkg_per_arcmin2=0.05, fov_arcmin=30.0,
                      source_xy=(0.0, 0.0), rng=rng)
    xc, yc, _, _ = localize_selfcal(d, half_width=5, n_grid=41)
    peaks, _, _ = recover_ring_radii(d, xc, yc, t0 + texp/2,
                                     min_prominence=0.15)
    for th, dpc in peaks:
        # match to nearest true cloud
        k = int(np.argmin([abs(dpc - c.d_pc) for c in clm]))
        d_true = clm[k].d_pc
        meas[d_true].append(th)
        meas_t[d_true].append(th_hr)

tt = np.linspace(1.5, 13, 100)
colors = {60: 'C0', 150: 'C1', 400: 'C2'}
for c in clm:
    ax[1].plot(tt, theta_ring_arcmin(tt*3600, c.d_pc), '-', color=colors[c.d_pc],
               label=f'{int(c.d_pc)} pc (theory)')
    ax[1].plot(meas_t[c.d_pc], meas[c.d_pc], 'o', color=colors[c.d_pc], ms=7,
               mec='k', mew=0.5)
ax[1].set_xlabel('delay time (hr)')
ax[1].set_ylabel('ring radius (arcmin)')
ax[1].set_title('Recovered ring radius vs time\n(points = measured, lines = $\\sqrt{2ct/d}$)')
ax[1].grid(True, alpha=0.3); ax[1].legend(fontsize=8.5)
fig.tight_layout(); fig.savefig("../figures/fig8_time_evolution.png"); plt.close(fig)

# quantify expansion for the single screen
r3 = theta_ring_arcmin(3*3600, 100); r9 = theta_ring_arcmin(9*3600, 100)
print(f"single screen 100pc: ring {r3:.2f}' at 3hr -> {r9:.2f}' at 9hr  (ratio {r9/r3:.3f}, sqrt3={np.sqrt(3):.3f})")
print("Saved fig8.")
