import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from halo import theta_ring_arcmin

plt.rcParams.update({"font.size": 11, "figure.dpi": 130})
det = json.load(open("../results/detect.json"))
N_list = det["N_list"]; thresh = det["thresh"]

fig, ax = plt.subplots(1, 2, figsize=(11, 4.7))

# Panel A: analytic detection significance in the ring annulus
theta = theta_ring_arcmin(6*3600, 100.0)
width = 1.0
annulus_area = 2*np.pi*theta*width
N = np.logspace(1, 3, 60)
for bkg, col in [(0.03, 'C0'), (0.07, 'C1'), (0.15, 'C2'), (0.30, 'C3')]:
    B = bkg * annulus_area
    ax[0].plot(N, N/np.sqrt(N + B), color=col, label=f'bkg {bkg} arcmin$^{{-2}}$')
ax[0].axhline(5, color='k', ls=':', lw=1); ax[0].text(11, 5.4, '5$\\sigma$', fontsize=9)
ax[0].set_xscale('log')
ax[0].set_xlabel('halo photons in annulus  $N$')
ax[0].set_ylabel('detection significance  $S/\\sqrt{S+B}$')
ax[0].set_title(f'Halo detectability\n(ring {theta:.0f}$\'$, annulus {annulus_area:.0f} arcmin$^2$)')
ax[0].grid(True, which='both', alpha=0.3); ax[0].legend(fontsize=8.5)

# Panel B: localization-success fraction (self-cal MC)
for bkg, col in [(0.07, 'C1'), (0.3, 'C3')]:
    key = f"bkg_{bkg}"
    if key in det:
        ax[1].plot(N_list, det[key], 'o-', color=col, label=f'bkg {bkg} arcmin$^{{-2}}$')
ax[1].set_xlabel('detected halo photons  $N$')
ax[1].set_ylabel(f'fraction localized to < {thresh:.0f}$^{{\\prime\\prime}}$')
ax[1].set_ylim(0, 1.05)
ax[1].set_title('Self-calibrating localization success\n(single screen)')
ax[1].grid(True, alpha=0.3); ax[1].legend(fontsize=8.5)
fig.tight_layout(); fig.savefig("../figures/fig9_detectability.png"); plt.close(fig)
print(f"annulus area {annulus_area:.1f} arcmin^2, ring {theta:.2f}'")
for sky in [100, 1000]:
    print(f"GW {sky} deg^2 -> 1\": area x{sky/(np.pi*(1/3600)**2):.1e} smaller")
print("Saved fig9.")
