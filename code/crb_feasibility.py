import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from halo import Cloud, simulate_halo

plt.rcParams.update({"font.size": 11, "figure.dpi": 130})
n_list = [30, 60, 120, 250, 500, 1000]
mc = json.load(open("../results/mc_results.json"))

# --- measure the ring radial width sigma_r (single 100 pc screen) -----
rng = np.random.default_rng(3)
d = simulate_halo([Cloud(100.0, 1.0, 0.1)], 6*3600, 1300, 1.0, n_signal=20000,
                  bkg_per_arcmin2=0.0, fov_arcmin=20.0, source_xy=(0.0, 0.0),
                  rng=rng)
s = d['is_signal']
r = np.hypot(d['x'][s], d['y'][s])
sigma_r = r.std()                       # arcmin
sigma_r_as = sigma_r * 60.0             # arcsec

# --- Fig 3 (regenerated): MC vs Cramer-Rao bound 2 sigma_r / sqrt(N) --
fig, ax = plt.subplots(figsize=(6.8, 5.1))
ns = np.array(n_list, float)
def rms(cfg): return np.array([mc[cfg][str(n)]["rms"] for n in n_list])
ax.plot(ns, rms("ideal"), 'o', color='C0', label='simulation, no background')
ax.plot(ns, rms("bkg"),   's', color='C1', label='simulation, bkg 0.07 arcmin$^{-2}$')
ax.plot(ns, rms("multi"), '^', color='C2', label='simulation, three screens')
Ngrid = np.logspace(np.log10(25), np.log10(1100), 100)
ax.plot(Ngrid, 2*sigma_r_as/np.sqrt(Ngrid), 'k--',
        label='Cramér–Rao bound  $2\\sigma_r/\\sqrt{N}$')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('number of detected halo photons  $N$')
ax.set_ylabel('localization error (arcsec, RMS)')
ax.set_title(f'Localization: simulation vs analytic bound\n($\\sigma_r$ = {sigma_r_as:.1f}$^{{\\prime\\prime}}$ measured ring width)')
ax.grid(True, which='both', alpha=0.3); ax.legend(fontsize=8.5)
fig.tight_layout(); fig.savefig("../figures/fig3_localization_scaling.png"); plt.close(fig)

# --- Fig 10: localization vs X-ray efficiency for a fiducial event ----
# NP20 fiducial: source at 400 Mpc, 100 pc dust path, 100 cm^2 effective
# area -> N ~ 4500 * eps_-5 halo photons over the halo lifetime (eps_-5 = eps/1e-5).
N_of_eps = lambda e5: 4500.0 * e5
eps5 = np.logspace(-2, 1, 100)
err = 2*sigma_r_as/np.sqrt(N_of_eps(eps5))

fig, ax = plt.subplots(figsize=(7.0, 5.1))
ax.plot(eps5, err, 'C3-', lw=2, label=r'$2\sigma_r/\sqrt{N},\ N=4500\,\epsilon_{-5}$')
# validation: MC points mapped to their equivalent efficiency
mc_eps = np.array(n_list)/4500.0
ax.plot(mc_eps, rms("ideal"), 'ko', ms=6, label='simulation')
ax.axhline(35, color='0.5', ls=':', lw=1.2)
ax.text(1.1e-2, 38, 'host-galaxy separation ($\\sim$35$^{\\prime\\prime}$)', fontsize=8.5, color='0.4')
ax.axhline(5, color='C0', ls=':', lw=1.2)
ax.text(1.1e-2, 5.4, '5$^{\\prime\\prime}$', fontsize=9, color='C0')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel(r'X-ray efficiency  $\epsilon_{-5} = \epsilon/10^{-5}$')
ax.set_ylabel('localization error (arcsec, RMS)')
ax.set_title('Expected localization for a fiducial event\n(400 Mpc, 100 pc dust, 100 cm$^2$ effective area)')
secax = ax.secondary_xaxis('top', functions=(N_of_eps, lambda N: N/4500.0))
secax.set_xlabel('detected halo photons  $N$')
ax.grid(True, which='both', alpha=0.3); ax.legend(fontsize=9, loc='lower left')
fig.tight_layout(); fig.savefig("../figures/fig11_feasibility.png"); plt.close(fig)

print(f"measured sigma_r = {sigma_r:.4f} arcmin = {sigma_r_as:.2f} arcsec")
print("N   MC_ideal(\")   CRB(\")")
for n in n_list:
    print(f"{n:5d}  {mc['ideal'][str(n)]['rms']:8.2f}   {2*sigma_r_as/np.sqrt(n):6.2f}")
for e5 in [0.01, 0.1, 1.0, 10.0]:
    print(f"eps_-5={e5:>5}: N={N_of_eps(e5):.0f}, error={2*sigma_r_as/np.sqrt(N_of_eps(e5)):.2f} arcsec")
print("Saved fig3 (regenerated) and fig11.")
