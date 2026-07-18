#!/usr/bin/env python3
"""Yeh-Hummer finite-size analysis of the Madrid-2019 NaCl(aq) sweep.

Reads the LAMMPS `fix ave/time` MSD files (one per cubic box), fits D from the
Einstein slope, checks the MSD is genuinely diffusive over the fit window, and
extrapolates D(L) -> D_0 via  D_pbc = D_0 - kT*xi/(6 pi eta L).

The eta that comes out of the D-vs-1/L slope is a prediction, not an input: no
viscosity calculation is needed to place the intercept.
"""
import argparse
import glob
import os
import re

import numpy as np

KB = 1.380649e-23      # J/K
XI = 2.837297          # Yeh-Hummer lattice constant (CUBIC cells only)
M2S_TO_A2PS = 1e8      # 1 m^2/s = 1e20 A^2 / 1e12 ps
SPECIES = ["O", "Na", "Cl"]


def load_msd(path, dt_fs=2.0):
    """-> t[ps], msd[nsp, nt] (A^2). LAMMPS ave/time: step c_msdO c_msdNa c_msdCl."""
    d = np.loadtxt(path, comments="#")
    return d[:, 0] * dt_fs / 1000.0, d[:, 1:4].T


def fit_D(t, msd, tmin, tmax):
    """Einstein fit MSD = 6 D t + c over [tmin,tmax]. Returns D (A^2/ps)."""
    m = (t >= tmin) & (t <= tmax)
    slope, _ = np.polyfit(t[m], msd[m], 1)
    return slope / 6.0


def loglog_slope(t, msd, tmin, tmax):
    """d(log MSD)/d(log t): 1.0 = diffusive, 2.0 = ballistic, <1 = subdiffusive."""
    m = (t >= tmin) & (t <= tmax) & (msd > 0) & (t > 0)
    return np.polyfit(np.log(t[m]), np.log(msd[m]), 1)[0]


def window_scatter(t, msd, tmin, tmax, nwin=4):
    """Spread of D over sub-windows: a cheap stability check, NOT a true error bar
    (compute msd uses a single time origin, so sub-windows are correlated).
    A real uncertainty needs independent seeds."""
    edges = np.linspace(tmin, tmax, nwin + 1)
    Ds = [fit_D(t, msd, edges[i], edges[i + 1]) for i in range(nwin)]
    return float(np.std(Ds))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/data/cerberus1/crm98/nacl_diffusion/"
                                      "lammps_madrid/finite_size")
    ap.add_argument("--tmin", type=float, default=200.0, help="fit start (ps)")
    ap.add_argument("--tmax", type=float, default=2000.0, help="fit end (ps)")
    ap.add_argument("--T", type=float, default=298.15)
    args = ap.parse_args()

    rows = []
    for d in sorted(glob.glob(os.path.join(args.root, "L*"))):
        L = float(re.search(r"L([\d.]+)$", d).group(1))
        f = os.path.join(d, f"L{L:.2f}.msd")
        if not os.path.exists(f):
            continue
        t, msd = load_msd(f)
        rec = {"L": L, "tmax_run": t[-1]}
        for i, s in enumerate(SPECIES):
            rec[f"D_{s}"] = fit_D(t, msd[i], args.tmin, args.tmax)
            rec[f"a_{s}"] = loglog_slope(t, msd[i], args.tmin, args.tmax)
            rec[f"s_{s}"] = window_scatter(t, msd[i], args.tmin, args.tmax)
        rows.append(rec)
    rows.sort(key=lambda r: r["L"])

    print(f"Einstein fit over {args.tmin:.0f}-{args.tmax:.0f} ps  "
          f"(D in 1e-9 m^2/s = 0.1 A^2/ps units shown as A^2/ps)\n")
    hdr = f"{'L(A)':>7} {'run(ps)':>8}"
    for s in SPECIES:
        hdr += f" | {'D_'+s:>8} {'+-':>7} {'loglog':>6}"
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        line = f"{r['L']:7.2f} {r['tmax_run']:8.0f}"
        for s in SPECIES:
            line += (f" | {r['D_'+s]:8.4f} {r['s_'+s]:7.4f} "
                     f"{r['a_'+s]:6.2f}")
        print(line)

    # ---- Yeh-Hummer extrapolation: D_pbc vs 1/L ----
    print(f"\nYeh-Hummer  D_pbc = D_0 - kT*xi/(6 pi eta L)   (xi={XI}, cubic)\n")
    print(f"{'sp':>3} {'D_0(A^2/ps)':>12} {'D_0(1e-9 m2/s)':>15} "
          f"{'R^2':>7} {'eta_fit(mPa s)':>15}")
    invL = np.array([1.0 / r["L"] for r in rows])
    for s in SPECIES:
        D = np.array([r[f"D_{s}"] for r in rows])
        slope, D0 = np.polyfit(invL, D, 1)
        pred = slope * invL + D0
        ss_res = np.sum((D - pred) ** 2)
        ss_tot = np.sum((D - D.mean()) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
        # slope = -kT xi/(6 pi eta)  [A^2/ps * A]  -> eta in Pa s
        # kT*xi/(6 pi) in J*m ; divide by |slope| converted to m^3/s
        slope_si = abs(slope) * 1e-20 / 1e-12 * 1e-10   # A^3/ps -> m^3/s
        eta = KB * args.T * XI / (6 * np.pi * slope_si) if slope < 0 else float("nan")
        print(f"{s:>3} {D0:12.4f} {D0*10:15.3f} {r2:7.3f} {eta*1e3:15.3f}")

    print("\nexperiment @298 K, 1 m NaCl:  D_water ~ 0.21 A^2/ps, "
          "D_Na ~ 0.117, D_Cl ~ 0.177 (1e-9 m2/s: 2.1 / 1.17 / 1.77)")
    print("eta(1 m NaCl, expt) ~ 0.95 mPa s")


if __name__ == "__main__":
    main()
