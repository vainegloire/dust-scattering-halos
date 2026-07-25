"""Generate Figure 4 (localization error vs background) from the accumulated
Monte-Carlo results in ../results/vsbkg_acc.json (produced by vsbkg_smooth.py)."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 11, "figure.dpi": 130})
a = json.load(open("../results/vsbkg_acc.json"))
b = a["bkg"]
rms = [np.sqrt(s/c) for s, c in zip(a["sumsq"], a["count"])]

fig, ax = plt.subplots(figsize=(6.4, 4.8))
ax.plot(b, rms, 'o-')
ax.set_ylim(0, 1.6)
ax.set_xlabel('background surface density (photons arcmin$^{-2}$)')
ax.set_ylabel('localization error (arcsec, RMS)')
ax.set_title(f'Localization vs background ($N=300$, {a["count"][0]} trials/point)')
ax.grid(True, alpha=0.3)
fig.tight_layout(); fig.savefig("../figures/fig4_error_vs_background.png"); plt.close(fig)
print("rms:", [round(x, 2) for x in rms])
print("Saved fig4 to ../figures/")
