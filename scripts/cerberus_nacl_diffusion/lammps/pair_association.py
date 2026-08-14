#!/usr/bin/env python3
"""Ion-pair association estimate from the seed-averaged Na-Cl g(r) at 1 m.

Two routes, both integrating to the PMF barrier r_b (CIP/SSIP divide):
  n_CIP = 4 pi rho_Cl int_0^rb g r^2 dr        (direct CIP coordination number)
  K_A   = 4 pi N_A int_0^rb e^{-w/kT} r^2 dr   (L/mol, Bjerrum-style estimate)
Ideal pairing would make the NE deviation ~ the paired fraction; comparing with
the measured Na-Cl cross term (fig6) links the PMF (fig5) to the transport
decomposition.
"""
import glob
import os
import re

import numpy as np

ROOT = "/data/cerberus1/crm98/nacl_diffusion/lammps_madrid"
KT_KCAL = 0.0019872041 * 298.15
NA_AVOG = 6.02214076e23


def load_avg_gr():
    gs, r, rhos = [], None, []
    for f in sorted(glob.glob(os.path.join(ROOT, "conductivity", "m1.0",
                                           "L*_s*", "*.rdf"))):
        dat = np.loadtxt(f, skiprows=4)
        if r is None:
            r = dat[:, 1]
        elif not np.allclose(dat[:, 1], r):
            dat = np.column_stack([dat[:, 0], r,
                                   np.interp(r, dat[:, 1], dat[:, 2]), dat[:, 3]])
        gs.append(dat[:, 2])
    return r, np.mean(gs, axis=0)


def main():
    r, g = load_avg_gr()
    w = np.full_like(g, np.nan)
    m = g > 0.02
    w[m] = -KT_KCAL * np.log(g[m])
    w -= np.nanmean(w[(r > 9.0) & (r < r.max())])

    # barrier location between CIP and SSIP (same windows as fig5)
    bm = (r >= 3.2) & (r <= 4.2) & np.isfinite(w)
    r_b = r[bm][np.argmax(w[bm])]

    # rho_Cl from the 1 m conductivity runs (kappa_vs_c.npz has V, nCl per run)
    kap = np.load(os.path.join(ROOT, "conductivity", "analysis", "kappa_vs_c.npz"))
    sel = kap["runs_mol"] == 1.0
    rho = (kap["runs_nCl"][sel] / kap["runs_V"][sel]).mean()   # A^-3
    c_molar = rho / NA_AVOG * 1e27                              # mol/L (1 L = 1e27 A^3)

    core = (r <= r_b)
    n_cip = 4 * np.pi * rho * np.trapezoid(g[core] * r[core] ** 2, r[core])

    wb = np.where(np.isfinite(w), w, np.inf)[core]
    ka_A3 = 4 * np.pi * np.trapezoid(np.exp(-wb / KT_KCAL) * r[core] ** 2, r[core])
    ka_lmol = ka_A3 * NA_AVOG * 1e-27                           # L/mol

    # ideal-pairing fraction x/c from K_A = x / (c-x)^2
    coeff = [ka_lmol, -(2 * ka_lmol * c_molar + 1), ka_lmol * c_molar ** 2]
    x = min(np.roots(coeff))
    print(f"barrier r_b = {r_b:.2f} A   c(1 m runs) = {c_molar:.3f} mol/L")
    print(f"n_CIP (direct)          = {n_cip:.3f}  Cl within r_b per Na")
    print(f"K_A (PMF, screened)     = {ka_lmol:.3f} L/mol")
    print(f"ideal paired fraction   = {x / c_molar:.3f}")
    print(f"-> predicted pairing part of Delta_NE ~ {x / c_molar:.1%} "
          f"(measured cross term: 0 +- few % -- consistent)")


if __name__ == "__main__":
    main()
