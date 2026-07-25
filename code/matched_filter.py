"""Matched-filter detection statistic for the scattering halo, with a
false-alarm / trials-factor treatment for tiling many fields.

The statistic at a trial center is  D(c) = sum_j T(|x_j - c|),  where T(r) is
the (zero-mean over the field) radial ring template.  Under background only,
E[D]=0; a halo produces a positive peak at its center.  The detection
statistic for a field is the maximum of D over trial centers.
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from halo import Cloud, simulate_halo, _ring_radial_pdf

plt.rcParams.update({"font.size": 11, "figure.dpi": 130})
rng = np.random.default_rng(2024)

FOV = 20.0
CL = [Cloud(100.0, 1.0, 0.1)]
T0, TEXP, E = 6*3600, 1300, 1.0
BKG = 0.07

# --- build the zero-mean radial matched-filter template ---------------
rg = np.linspace(0.0, FOV*np.sqrt(2.0), 1400)
s = _ring_radial_pdf(rg, CL, T0, TEXP, E, 0.15, n_r_samples=32)
within = rg <= FOV
w = 2*np.pi*rg
mean_s = np.sum(s[within]*w[within]) / np.sum(w[within])   # area-weighted mean
Tr = s - mean_s
Tr = Tr / np.abs(Tr).max()


def mf_max(data, half_width=6.0, n_grid=25):
    x, y = data['x'], data['y']
    gx = np.linspace(-half_width, half_width, n_grid)
    best = -np.inf
    for yc in gx:
        yy = (y - yc)**2
        for xc in gx:
            r = np.sqrt((x - xc)**2 + yy)
            v = np.interp(r, rg, Tr).sum()
            if v > best:
                best = v
    return best


def make(n_signal, seed):
    r = np.random.default_rng(seed)
    return simulate_halo(CL, T0, TEXP, E, n_signal=n_signal,
                         bkg_per_arcmin2=BKG, fov_arcmin=FOV,
                         source_xy=(0.0, 0.0), rng=r)


# --- null distribution (background only) ------------------------------
Knull = 400
Dnull = np.array([mf_max(make(0, 5000+i)) for i in range(Knull)])
mu0, sd0 = Dnull.mean(), Dnull.std()
# Gumbel fit (method of moments) to the max statistic
beta = sd0*np.sqrt(6)/np.pi
mu_g = mu0 - 0.5772*beta

# --- signal: mean statistic vs N --------------------------------------
Nsig = [30, 60, 120, 240]
Ksig = 40
Dsig_mean = []
Dsig_all = {}
for N in Nsig:
    vals = np.array([mf_max(make(N, 9000+N*10+i)) for i in range(Ksig)])
    Dsig_mean.append(vals.mean())
    Dsig_all[N] = vals
Dsig_mean = np.array(Dsig_mean)
# linear fit E[D] = c0 + a N
A = np.vstack([np.ones_like(Nsig, float), Nsig]).T
c0, a = np.linalg.lstsq(A, Dsig_mean, rcond=None)[0]

# --- required photons for a GLOBAL 5-sigma detection over M fields ----
alpha_global = 2.87e-7                 # one-sided 5 sigma
Mfields = np.logspace(0, 4, 60)
Tthr = mu_g + beta*np.log(Mfields/alpha_global)   # Gumbel tail threshold
Nreq = np.clip((Tthr - c0)/a, 0, None)

# single-field 5 sigma threshold and the significance of a faint halo
Tthr_1 = mu_g + beta*np.log(1.0/alpha_global)
sig45 = (np.mean([mf_max(make(45, 20000+i)) for i in range(30)]) - mu0)/sd0

# ---------------- figure ----------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
ax[0].hist(Dnull, bins=22, color='0.6', alpha=0.8, label='background only')
ax[0].hist(Dsig_all[60], bins=22, color='C3', alpha=0.7, label='halo, $N=60$')
ax[0].axvline(Tthr_1, color='k', ls='--', lw=1.3, label='single-field 5$\\sigma$')
ax[0].set_xlabel('matched-filter statistic  $D_{\\max}$')
ax[0].set_ylabel('number of trials')
ax[0].set_title('Detection statistic:\nbackground vs a faint halo')
ax[0].legend(fontsize=8.5)

ax[1].plot(Mfields, Nreq, 'C0-', lw=2)
ax[1].set_xscale('log')
ax[1].set_xlabel('number of tiled fields  $M$ (trials)')
ax[1].set_ylabel('halo photons for global $5\\sigma$')
ax[1].set_title('Detection threshold vs survey size')
ax[1].grid(True, which='both', alpha=0.3)
fig.tight_layout(); fig.savefig("../figures/fig10_matched_filter.png"); plt.close(fig)

res = dict(mu0=float(mu0), sd0=float(sd0), c0=float(c0), a=float(a),
           sig45=float(sig45),
           Nreq_M1=float(np.clip((Tthr_1-c0)/a, 0, None)),
           Nreq_M100=float(np.clip((mu_g+beta*np.log(100/alpha_global)-c0)/a, 0, None)),
           Nreq_M1000=float(np.clip((mu_g+beta*np.log(1000/alpha_global)-c0)/a, 0, None)))
json.dump(res, open("../results/matched_filter.json", "w"), indent=2)
print(f"null: mu0={mu0:.2f} sd0={sd0:.2f}  Gumbel(mu={mu_g:.2f},beta={beta:.2f})")
print(f"signal fit E[D]=c0+a*N: c0={c0:.2f}, a={a:.4f}")
print(f"single faint halo N=45 significance = {sig45:.1f} sigma")
print(f"required halo photons for global 5sigma: M=1 -> {res['Nreq_M1']:.0f}, "
      f"M=100 -> {res['Nreq_M100']:.0f}, M=1000 -> {res['Nreq_M1000']:.0f}")
print("Saved fig10.")
