"""Unit test for the relative photon yields of multi-screen halos.

A screen at distance d delivers photons at a rate proportional to
tau * w(theta) * c/d, because theta dtheta/dt = c/d.  Until 28 Sept 2026 the
simulator drew screens in proportion to tau alone; these checks would have
failed then (test 1 would have given 1:1 instead of 2:1).

    python3 test_yields.py        # prints PASS or raises AssertionError
"""
import numpy as np
from halo import (Cloud, simulate_halo, theta_ring_arcmin, scattering_weight,
                  _ring_radial_pdf)

T0, TEXP = 6 * 3600.0, 1300.0
trapezoid = getattr(np, "trapezoid", None) or np.trapz     # NumPy 1.x and 2.x


def split_counts(d, r_cut):
    r = np.hypot(d["x"], d["y"])[d["is_signal"]]
    return np.sum(r > r_cut), np.sum(r <= r_cut)       # near screen, far screen


def expected_ratio(c1, c2, E_keV):
    t = np.linspace(T0, T0 + TEXP, 2001)
    rate = lambda c: c.tau / c.d_pc * trapezoid(
        scattering_weight(theta_ring_arcmin(t, c.d_pc), E_keV, c.a_um), t)
    return rate(c1) / rate(c2)


def check(c1, c2, E_keV, n, r_cut, seed, label):
    d = simulate_halo([c1, c2], T0, TEXP, E_keV=E_keV, n_signal=n,
                      rng=np.random.default_rng(seed))
    n1, n2 = split_counts(d, r_cut)
    p = n1 / (n1 + n2)
    ratio, want = n1 / n2, expected_ratio(c1, c2, E_keV)
    se = np.sqrt(p * (1 - p) / (n1 + n2)) / (1 - p) ** 2
    ok = abs(ratio - want) < 4 * se
    print(f"{label}: sampler {ratio:.3f} +/- {se:.3f}, expected {want:.3f}"
          f"  {'ok' if ok else 'FAIL'}")
    assert ok, label
    return want


# 1. Flat cross section (theta0 >> theta): equal-tau screens at d and 2d
#    must give exactly 2 : 1 per unit time.
want = check(Cloud(100, 1.0), Cloud(200, 1.0), 0.01, 40000, 6.0, 1,
             "flat w, 100 vs 200 pc")
assert abs(want - 2.0) < 1e-3

# 2. The paper's cross section at 1 keV: equal-tau screens at 60 and 400 pc
#    (from first principles about 4.7 : 1; the old sampler gave about 0.7).
check(Cloud(60, 1.0), Cloud(400, 1.0), 1.0, 40000, 6.0, 2,
      "1 keV, 60 vs 400 pc")

# 3. The likelihood template must use the same yields as the sampler:
#    integrate its radial intensity over each ring and compare.
c1, c2 = Cloud(60, 1.0), Cloud(400, 1.0)
r = np.linspace(0.01, 15.0, 30000)
lam = _ring_radial_pdf(r, [c1, c2], T0, TEXP, 1.0, 0.15)
dn = lam * 2 * np.pi * r                     # photons per unit radius
near, far = dn[r > 6.0].sum(), dn[r <= 6.0].sum()
want = expected_ratio(c1, c2, 1.0)
print(f"template: {near / far:.3f}, expected {want:.3f}")
assert abs(near / far / want - 1) < 0.01

print("PASS")
