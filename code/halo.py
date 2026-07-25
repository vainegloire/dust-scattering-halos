"""
halo.py
--------
Forward model and localization tools for Galactic dust-scattering halos
produced by X-ray flashes from gravitational-wave sources.

Physics follows Nederlander & Paerels (2020, ApJ 890, 135):

  * A prompt X-ray flash from a source at (effectively) infinity illuminates
    Galactic dust clouds along the line of sight.
  * A photon scattered by angle theta at a thin dust screen a distance d from
    the observer is seen at angular offset theta from the source direction and
    arrives with a time delay

        t = (d / 2c) * theta**2          (small-angle, source at infinity)

    so that at observation time t the screen produces a ring of angular radius

        theta_ring(t, d) = sqrt(2 c t / d).

  * The differential scattering cross section is approximately Gaussian in the
    scattering angle (Mauche & Gorenstein 1986; Rayleigh-Gans),

        dsigma/dOmega ~ exp(-theta**2 / (2 theta0**2)),
        theta0 = 10.4 * (E/keV)**-1 * (a/0.1um)**-1 arcmin.

This module keeps ONLY a dependence on numpy (matplotlib is used in the figure
scripts, not here).  No scipy required.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Physical constants (cgs) and unit helpers
# ---------------------------------------------------------------------------
C_CGS = 2.99792458e10          # speed of light, cm / s
PC_CGS = 3.0856775814913673e18  # parsec, cm
RAD_TO_ARCMIN = (180.0 / np.pi) * 60.0
ARCMIN_TO_RAD = 1.0 / RAD_TO_ARCMIN


def theta_ring_arcmin(t_s, d_pc):
    """Angular radius (arcmin) of the scattering ring from a screen at
    distance d_pc (parsec) observed at delay t_s (seconds)."""
    t_s = np.asarray(t_s, dtype=float)
    d_cm = np.asarray(d_pc, dtype=float) * PC_CGS
    theta_rad = np.sqrt(2.0 * C_CGS * t_s / d_cm)
    return theta_rad * RAD_TO_ARCMIN


def theta0_arcmin(E_keV=1.0, a_um=0.1):
    """Characteristic (Gaussian) scattering angle theta0 in arcmin."""
    return 10.4 * (E_keV ** -1) * ((a_um / 0.1) ** -1)


def scattering_weight(theta_arcmin, E_keV=1.0, a_um=0.1):
    """Relative Gaussian scattering cross section evaluated at scattering
    angle theta (arcmin)."""
    th0 = theta0_arcmin(E_keV, a_um)
    return np.exp(-(np.asarray(theta_arcmin, float) ** 2) / (2.0 * th0 ** 2))


# ---------------------------------------------------------------------------
# Cloud / observation description
# ---------------------------------------------------------------------------
class Cloud:
    """A thin dust screen along the line of sight.

    Parameters
    ----------
    d_pc   : distance from observer (parsec)
    tau    : scattering optical depth (sets relative photon yield)
    a_um   : characteristic grain radius (micron)
    """
    def __init__(self, d_pc, tau=1.0, a_um=0.1):
        self.d_pc = float(d_pc)
        self.tau = float(tau)
        self.a_um = float(a_um)


# ---------------------------------------------------------------------------
# Forward model: simulate a detected photon list
# ---------------------------------------------------------------------------
def sample_grain_sizes(n, rng, amin=0.05, amax=0.25, power=0.5):
    """Draw grain radii (micron) from a scattering-weighted size distribution.
    An MRN number distribution n(a) ~ a^-3.5 times an X-ray scattering cross
    section sigma_scat ~ a^4 gives an effective p(a) ~ a^{+0.5}; larger grains
    dominate the scattered flux.  Sampled by inverse-CDF of a power law.
    """
    p = power + 1.0
    u = rng.random(n)
    return (amin**p + u * (amax**p - amin**p)) ** (1.0 / p)


def sample_disk_sightline(rng, n_clouds_mean=10, d_scale_pc=150.0,
                          d_max_pc=1000.0, a_um=0.1):
    """Generate a realistic set of dust screens along a line of sight.
    Cloud distances follow an exponential (disk-like) profile and column
    densities (-> relative optical depths) are log-normal.
    """
    n = max(1, int(rng.poisson(n_clouds_mean)))
    d = rng.exponential(d_scale_pc, n)
    d = d[(d > 20.0) & (d < d_max_pc)]
    if d.size == 0:
        d = np.array([d_scale_pc])
    tau = rng.lognormal(mean=0.0, sigma=0.8, size=d.size)
    return [Cloud(float(di), tau=float(wi), a_um=a_um) for di, wi in zip(d, tau)]


def simulate_halo(clouds, t_start_s, t_exp_s, E_keV=1.0,
                  n_signal=300, bkg_per_arcmin2=0.0, fov_arcmin=20.0,
                  psf_sigma_arcmin=0.15, source_xy=(0.0, 0.0),
                  time_weight_index=0.0, grain_spread=None, rng=None):
    """Monte-Carlo a detected photon list (offsets in arcmin) for an exposure
    that runs from t_start_s to t_start_s + t_exp_s after the prompt flash.

    Returns
    -------
    dict with keys:
       'x', 'y'          : arcmin offsets of ALL detected photons
       'is_signal'       : boolean mask, True for scattered (halo) photons
       'clouds'          : the input clouds
       't_start','t_exp' : exposure window (s)
       'E_keV','fov'     : passthrough
    """
    if rng is None:
        rng = np.random.default_rng()
    x0, y0 = source_xy

    # ---- cloud-choice probabilities -------------------------------------
    # Candidates are drawn per cloud proportional to optical depth ONLY;
    # the angular cross section enters once, through the acceptance step
    # below, so that accepted photons per cloud scale as tau_i * w(theta_i).
    # (Weighting the choice by the cross section as well would double-count
    # it and suppress outer rings as w^2.)
    yields = np.asarray([c.tau for c in clouds], float)
    yields = yields / yields.sum()

    xs, ys = [], []
    # oversample then rejection-accept on the (radius-dependent) cross section
    n_left = n_signal
    guard = 0
    while n_left > 0 and guard < 200:
        guard += 1
        m = int(np.ceil(n_left * 1.6)) + 8
        # choose cloud per candidate
        ci = rng.choice(len(clouds), size=m, p=yields)
        # choose scatter time in window; optional t^-index brightness weighting
        u = rng.random(m)
        if time_weight_index == 0.0:
            t = t_start_s + u * t_exp_s
        else:
            # sample t within window with pdf ~ t^{-index}
            a, b = t_start_s, t_start_s + t_exp_s
            p = 1.0 - time_weight_index
            if abs(p) < 1e-9:
                t = a * (b / a) ** u
            else:
                t = (a ** p + u * (b ** p - a ** p)) ** (1.0 / p)
        d = np.array([clouds[k].d_pc for k in ci])
        if grain_spread is not None:
            a_um = sample_grain_sizes(len(ci), rng, *grain_spread)
            # dsigma/dOmega at fixed theta carries a 1/theta0^2 ~ a^2
            # normalization on top of the Gaussian shape; normalize by the
            # largest grain so the acceptance probability stays <= 1.
            size_fac = (a_um / grain_spread[1]) ** 2
        else:
            a_um = np.array([clouds[k].a_um for k in ci])
            size_fac = 1.0
        th = theta_ring_arcmin(t, d)
        # rejection against the differential cross section at the actual radius
        acc_prob = size_fac * scattering_weight(th, E_keV, a_um)
        keep = rng.random(m) < acc_prob
        th = th[keep]
        if th.size == 0:
            continue
        phi = rng.random(th.size) * 2.0 * np.pi
        cx = x0 + th * np.cos(phi)
        cy = y0 + th * np.sin(phi)
        # PSF blur
        cx = cx + rng.normal(0.0, psf_sigma_arcmin, cx.size)
        cy = cy + rng.normal(0.0, psf_sigma_arcmin, cy.size)
        # keep only those inside FoV
        inside = (np.abs(cx) <= fov_arcmin) & (np.abs(cy) <= fov_arcmin)
        cx, cy = cx[inside], cy[inside]
        take = min(cx.size, n_left)
        xs.append(cx[:take]); ys.append(cy[:take])
        n_left -= take

    xs = np.concatenate(xs) if xs else np.array([])
    ys = np.concatenate(ys) if ys else np.array([])
    is_sig = np.ones(xs.size, bool)

    # ---- background: uniform over the square field ---------------------
    area = (2.0 * fov_arcmin) ** 2
    n_bkg = rng.poisson(bkg_per_arcmin2 * area)
    if n_bkg > 0:
        bx = rng.uniform(-fov_arcmin, fov_arcmin, n_bkg)
        by = rng.uniform(-fov_arcmin, fov_arcmin, n_bkg)
        xs = np.concatenate([xs, bx])
        ys = np.concatenate([ys, by])
        is_sig = np.concatenate([is_sig, np.zeros(n_bkg, bool)])

    return dict(x=xs, y=ys, is_signal=is_sig, clouds=clouds,
                t_start=t_start_s, t_exp=t_exp_s, E_keV=E_keV,
                fov=fov_arcmin, psf_sigma=psf_sigma_arcmin)


# ---------------------------------------------------------------------------
# Localization: maximum-likelihood centroiding
# ---------------------------------------------------------------------------
def _ring_radial_pdf(r, clouds, t_start_s, t_exp_s, E_keV, psf_sigma,
                     n_r_samples=48):
    """Un-normalized radial intensity (per unit area) at radius r produced by
    all clouds, given the exposure window.  Built by integrating the thin ring
    over the exposure time (top-hat in delay), weighting each cloud by its
    scattering yield, and smearing radially by the PSF.
    """
    r = np.asarray(r, float)
    ts = np.linspace(t_start_s, t_start_s + t_exp_s, n_r_samples)
    t_mid = t_start_s + 0.5 * t_exp_s

    # Build a flat list of thin-ring radii and their weights across all clouds
    radii = []
    weights = []
    for c in clouds:
        th_mid = theta_ring_arcmin(t_mid, c.d_pc)
        yield_w = c.tau * scattering_weight(th_mid, E_keV, c.a_um)
        th_of_t = theta_ring_arcmin(ts, c.d_pc)          # ring radius vs time
        radii.append(th_of_t)
        weights.append(np.full(th_of_t.shape, yield_w))
    radii = np.concatenate(radii)                        # (M,)
    weights = np.concatenate(weights)                    # (M,)

    # Vectorized PSF-smeared sum of thin rings: (N, M) -> (N,)
    diff = r[:, None] - radii[None, :]
    g = np.exp(-(diff ** 2) / (2.0 * psf_sigma ** 2)) * weights[None, :]
    intensity = g.sum(axis=1) / (2.0 * np.pi * np.maximum(r, 1e-3))
    intensity /= len(ts)
    return intensity


def build_log_template(data, bkg_frac_guess=0.3, clouds=None,
                       n_r_samples=32, r_max=None, n_r_grid=1200):
    """Precompute log(mixture intensity) on a 1-D radial grid so that the
    log-likelihood at any candidate center is just an interpolation in r.

    Returns (r_grid, log_lam_grid) for use with np.interp.
    """
    clouds = clouds if clouds is not None else data['clouds']
    fov = data['fov']
    if r_max is None:
        r_max = 2.0 * fov * np.sqrt(2.0)
    r_grid = np.linspace(0.0, r_max, n_r_grid)
    sig = _ring_radial_pdf(r_grid, clouds, data['t_start'], data['t_exp'],
                           data['E_keV'], data['psf_sigma'], n_r_samples)
    s = sig / (sig.mean() + 1e-12)
    bkg = 1.0                                    # flat, arbitrary units
    lam = (1.0 - bkg_frac_guess) * s + bkg_frac_guess * bkg
    lam = np.maximum(lam, 1e-12)
    return r_grid, np.log(lam)


def loglike_center(xc, yc, data, template):
    """Log-likelihood of a candidate center using a precomputed template."""
    r_grid, log_lam = template
    r = np.sqrt((data['x'] - xc) ** 2 + (data['y'] - yc) ** 2)
    return np.sum(np.interp(r, r_grid, log_lam))


def localize(data, half_width=6.0, n_grid=61, refine=True,
             bkg_frac_guess=0.3, n_r_samples=32):
    """Grid-search the center that maximizes the log-likelihood.

    Returns (xc, yc, loglike_grid, extent).
    """
    template = build_log_template(data, bkg_frac_guess=bkg_frac_guess,
                                  n_r_samples=n_r_samples)
    gx = np.linspace(-half_width, half_width, n_grid)
    gy = np.linspace(-half_width, half_width, n_grid)
    L = np.empty((n_grid, n_grid))
    for i, yc in enumerate(gy):
        for j, xc in enumerate(gx):
            L[i, j] = loglike_center(xc, yc, data, template)
    i0, j0 = np.unravel_index(np.argmax(L), L.shape)
    xc, yc = gx[j0], gy[i0]
    if refine:
        step = gx[1] - gx[0]
        for _stage in range(2):                  # two-stage local refinement
            fx = np.linspace(xc - step, xc + step, 21)
            fy = np.linspace(yc - step, yc + step, 21)
            bestL, bxy = -np.inf, (xc, yc)
            for yc2 in fy:
                for xc2 in fx:
                    v = loglike_center(xc2, yc2, data, template)
                    if v > bestL:
                        bestL, bxy = v, (xc2, yc2)
            xc, yc = bxy
            step = (fx[1] - fx[0])               # zoom in for next stage
    return xc, yc, L, (-half_width, half_width, -half_width, half_width)


# ---------------------------------------------------------------------------
# Self-calibrating localization: no knowledge of dust distances required
# ---------------------------------------------------------------------------
def _radial_concentration_score(r, r_max, bw=0.30, nbin=400):
    """Score how tightly a set of radii clusters into rings.  Builds a
    Gaussian-smoothed radial density from the data itself and returns
    sum_j log p(r_j).  Divides out the geometric r d r area element so that a
    uniform (background-only) field is flat and does not masquerade as rings.
    """
    edges = np.linspace(0.0, r_max, nbin + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    hist, _ = np.histogram(r, bins=edges)
    # Gaussian smoothing kernel in radius
    dr = centers[1] - centers[0]
    half = int(np.ceil(3 * bw / dr))
    k = np.exp(-(np.arange(-half, half + 1) * dr) ** 2 / (2 * bw ** 2))
    k /= k.sum()
    dens = np.convolve(hist, k, mode='same')
    dens = dens / np.maximum(centers, 0.5 * dr)      # remove r d r area factor
    dens = dens / (dens.sum() * dr + 1e-12) + 1e-9   # normalize + floor
    idx = np.clip((r / r_max * nbin).astype(int), 0, nbin - 1)
    return np.sum(np.log(dens[idx]))


def localize_selfcal(data, half_width=6.0, n_grid=61, refine=True,
                     bw=0.30):
    """Locate the source as the common center of concentric rings WITHOUT
    knowing the dust screen distances or how many there are.  Maximizes the
    radial concentration score over trial centers.

    Returns (xc, yc, score_grid, extent).
    """
    x, y = data['x'], data['y']
    r_max = data['fov'] * np.sqrt(2.0)
    gx = np.linspace(-half_width, half_width, n_grid)
    gy = np.linspace(-half_width, half_width, n_grid)
    S = np.empty((n_grid, n_grid))
    for i, yc in enumerate(gy):
        rr_y = (y - yc) ** 2
        for j, xc in enumerate(gx):
            r = np.sqrt((x - xc) ** 2 + rr_y)
            S[i, j] = _radial_concentration_score(r, r_max, bw)
    i0, j0 = np.unravel_index(np.argmax(S), S.shape)
    xc, yc = gx[j0], gy[i0]
    if refine:
        step = gx[1] - gx[0]
        for _stage in range(2):
            fx = np.linspace(xc - step, xc + step, 21)
            fy = np.linspace(yc - step, yc + step, 21)
            best, bxy = -np.inf, (xc, yc)
            for yc2 in fy:
                rr_y = (y - yc2) ** 2
                for xc2 in fx:
                    r = np.sqrt((x - xc2) ** 2 + rr_y)
                    v = _radial_concentration_score(r, r_max, bw)
                    if v > best:
                        best, bxy = v, (xc2, yc2)
            xc, yc = bxy
            step = fx[1] - fx[0]
    return xc, yc, S, (-half_width, half_width, -half_width, half_width)


def recover_ring_radii(data, xc, yc, t_obs_s, bw=0.25, nbin=600,
                       min_prominence=0.25):
    """Given a recovered center, find the ring radii (radial-density peaks)
    and convert each to a dust-screen distance via d = 2 c t / theta^2.

    Returns list of (theta_arcmin, d_pc) for detected peaks.
    """
    r = np.sqrt((data['x'] - xc) ** 2 + (data['y'] - yc) ** 2)
    r_max = data['fov'] * np.sqrt(2.0)
    edges = np.linspace(0.0, r_max, nbin + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    hist, _ = np.histogram(r, bins=edges)
    dr = centers[1] - centers[0]
    half = int(np.ceil(3 * bw / dr))
    k = np.exp(-(np.arange(-half, half + 1) * dr) ** 2 / (2 * bw ** 2))
    k /= k.sum()
    dens = np.convolve(hist, k, mode='same') / np.maximum(centers, 0.5 * dr)
    if dens.max() <= 0:
        return [], centers, dens
    d = dens / dens.max()
    peaks = []
    for i in range(2, len(d) - 2):
        if d[i] >= d[i-1] and d[i] > d[i+1] and d[i] >= min_prominence:
            peaks.append(i)
    # merge peaks closer than the bandwidth
    merged = []
    for p in peaks:
        if merged and (centers[p] - centers[merged[-1]]) < bw:
            if d[p] > d[merged[-1]]:
                merged[-1] = p
        else:
            merged.append(p)
    out = []
    for p in merged:
        th = centers[p]
        th_rad = th * ARCMIN_TO_RAD
        d_pc = 2.0 * C_CGS * t_obs_s / (th_rad ** 2) / PC_CGS
        out.append((float(th), float(d_pc)))
    return out, centers, dens
