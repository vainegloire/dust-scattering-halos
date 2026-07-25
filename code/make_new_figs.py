import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from halo import (Cloud, simulate_halo, localize_selfcal,
                  recover_ring_radii, theta_ring_arcmin)

plt.rcParams.update({"font.size": 11, "figure.dpi": 130})
n_list = [30, 60, 120, 250, 500, 1000]
tmpl = json.load(open("../results/mc_results.json"))
self_ = json.load(open("../results/selfcal_results.json"))


def med(d, cfg):
    return np.array([d[cfg][str(n)]["median"] for n in n_list])


# ---- Figure 5: template (known distances) vs self-calibrating --------
fig, ax = plt.subplots(figsize=(6.6, 5.0))
ns = np.array(n_list, float)
ax.plot(ns, med(tmpl, "multi"),  'o-', color='C0',
        label='known dust distances (template)')
ax.plot(ns, med(self_, "multi"), '^--', color='C3',
        label='unknown distances (self-calibrating)')
ax.plot(ns, med(tmpl, "multi")[2] * np.sqrt(ns[2] / ns), 'k:',
        label=r'$\propto N^{-1/2}$ reference')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('number of detected halo photons  $N$')
ax.set_ylabel('median localization error (arcsec)')
ax.set_title('Localizing without knowing the dust distribution\n(three-screen halo)')
ax.grid(True, which='both', alpha=0.3); ax.legend(fontsize=9)
fig.tight_layout(); fig.savefig("../figures/fig5_selfcal_vs_template.png"); plt.close(fig)

# ---- Figure 6: halo tomography ---------------------------------------
rng = np.random.default_rng(7)
clm = [Cloud(60, 1.2, 0.1), Cloud(150, 1.3, 0.1), Cloud(400, 1.5, 0.1)]
t_obs = 6*3600 + 650
d = simulate_halo(clm, 6*3600, 1300, 1.0, n_signal=800,
                  bkg_per_arcmin2=0.07, fov_arcmin=20.0,
                  source_xy=(0.0, 0.0), rng=rng)
xc, yc, _, _ = localize_selfcal(d, half_width=6, n_grid=51)
peaks, centers, dens = recover_ring_radii(d, xc, yc, t_obs, min_prominence=0.15)
true_theta = [theta_ring_arcmin(t_obs, c.d_pc) for c in clm]
true_d = [c.d_pc for c in clm]

fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
ax[0].plot(centers, dens / dens.max(), color='C0')
for th in true_theta:
    ax[0].axvline(th, color='0.6', ls=':', lw=1)
for th, dd in peaks:
    ax[0].axvline(th, color='C3', ls='--', lw=1.2)
ax[0].set_xlim(0, 14)
ax[0].set_xlabel('ring radius (arcmin)')
ax[0].set_ylabel('normalized radial density')
ax[0].set_title('Radial profile about recovered center\n(dotted = truth, dashed = detected peaks)')

rec_d = [dd for _, dd in peaks]
matched_true = []
for th, dd in peaks:                       # match each detected ring to nearest truth
    k = int(np.argmin([abs(th - tt) for tt in true_theta]))
    matched_true.append(true_d[k])
ax[1].plot([0, 450], [0, 450], 'k:', lw=1)
ax[1].scatter(matched_true, rec_d, s=80, color='C3', zorder=3)
for mt, rd in zip(matched_true, rec_d):
    ax[1].annotate(f'{rd:.0f} pc', (mt, rd), textcoords='offset points',
                   xytext=(8, -4), fontsize=9)
ax[1].set_xlabel('true screen distance (pc)')
ax[1].set_ylabel('recovered distance (pc)')
ax[1].set_title('Dust-screen distances recovered from ring radii')
ax[1].set_xlim(0, 450); ax[1].set_ylim(0, 450)
ax[1].grid(True, alpha=0.3)
fig.tight_layout(); fig.savefig("../figures/fig6_tomography.png"); plt.close(fig)

# ---- catastrophic-failure fraction at low N (single & multi) ---------
def failure_fraction(clouds, n, thresh=20.0, K=60, seed=99):
    r = np.random.default_rng(seed)
    fails = 0
    for _ in range(K):
        dd = simulate_halo(clouds, 6*3600, 1300, 1.0, n_signal=n,
                           bkg_per_arcmin2=0.07, fov_arcmin=20.0,
                           source_xy=(0.0, 0.0), rng=r)
        xx, yy, _, _ = localize_selfcal(dd, half_width=5, n_grid=41)
        if np.hypot(xx, yy) * 60.0 > thresh:
            fails += 1
    return fails / K


f30 = failure_fraction(clm, 30)
f60 = failure_fraction(clm, 60)
print("tomography detected peaks (theta_arcmin, d_pc):", [(round(t,2), round(dd)) for t,dd in peaks])
print("true (theta, d):", [(round(tt,2), dd) for tt, dd in zip(true_theta, true_d)])
print(f"self-cal catastrophic-failure fraction (>20 arcsec): N=30 -> {f30:.2f}, N=60 -> {f60:.2f}")
print("Saved fig5, fig6.")
