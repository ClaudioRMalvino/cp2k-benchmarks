#!/usr/bin/env python3
"""Seed-replica analysis of the Madrid-2019 NaCl(aq) finite-size study.

Reads the multi-time-origin MSDs (.msdmo) from the 4-box x 5-seed replica
array, fits D per species per replica, and does two Yeh-Hummer treatments:

  free fit    : slope and intercept fitted per species (water should be clean;
                ions historically too noisy for a meaningful slope)
  shared slope: the YH correction k_B T xi/(6 pi eta L) contains no property
                of the diffusing species, so all species share the water slope;
                only the ion intercepts D_0 are fitted.

Error bars are standard errors over the 5 independent seeds -- real,
independent-sample errors, unlike the single-origin window scatter.
"""
import glob
import os
import re

import numpy as np

KB = 1.380649e-23
XI = 2.837297
T_K = 298.15
ROOT = "/data/cerberus1/crm98/nacl_diffusion/lammps_madrid/replicas"
OUT = "/data/cerberus1/crm98/nacl_diffusion/lammps_madrid/replicas/analysis"
SPECIES = ["O", "Na", "Cl"]
# Multi-origin MSD: statistics degrade at lags approaching the 2 ns run length
# (few independent origins), so fit early-mid lags, well past the sub-ps
# ballistic/caging regime.
TMIN, TMAX = 20.0, 400.0


def fit_D(t, msd):
    m = (t >= TMIN) & (t <= TMAX)
    sl, _ = np.polyfit(t[m], msd[m], 1)
    return sl / 6.0


def loglog_slope(t, msd):
    m = (t >= TMIN) & (t <= TMAX) & (t > 0) & (msd > 0)
    return np.polyfit(np.log(t[m]), np.log(msd[m]), 1)[0]


def main():
    os.makedirs(OUT, exist_ok=True)
    data = {}  # L -> species -> [D per seed]
    diag = {}  # L -> species -> [loglog slope per seed]
    for d in sorted(glob.glob(os.path.join(ROOT, "L*_s*"))):
        mo = re.search(r"L([\d.]+)_s(\d+)$", d)
        L, seed = float(mo.group(1)), mo.group(2)
        f = os.path.join(d, f"L{L:.2f}_s{seed}.msdmo")
        raw = np.loadtxt(f, comments="#")
        t = raw[:, 0]
        for i, s in enumerate(SPECIES):
            data.setdefault(L, {}).setdefault(s, []).append(fit_D(t, raw[:, i + 1]))
            diag.setdefault(L, {}).setdefault(s, []).append(loglog_slope(t, raw[:, i + 1]))

    Ls = np.array(sorted(data))
    invL = 1.0 / Ls
    print(f"fit window {TMIN}-{TMAX} ps, {len(Ls)} boxes, "
          f"{len(data[Ls[0]]['O'])} seeds/box\n")

    mean, sem = {}, {}
    print(f"{'L':>7} {'spec':>4} {'D mean':>9} {'SEM':>8} {'loglog':>7}")
    for s in SPECIES:
        mean[s] = np.array([np.mean(data[L][s]) for L in Ls])
        sem[s] = np.array([np.std(data[L][s], ddof=1) / np.sqrt(len(data[L][s]))
                           for L in Ls])
        for j, L in enumerate(Ls):
            print(f"{L:7.2f} {s:>4} {mean[s][j]:9.4f} {sem[s][j]:8.4f} "
                  f"{np.mean(diag[L][s]):7.3f}")

    # --- free YH fit per species (weighted by seed SEM) ---
    print("\nfree fits (slope + intercept per species):")
    results = {}
    for s in SPECIES:
        w = 1.0 / sem[s] ** 2
        (sl, D0), cov = np.polyfit(invL, mean[s], 1, w=np.sqrt(w), cov=True)
        pred = sl * invL + D0
        r2 = 1 - np.sum((mean[s] - pred) ** 2) / np.sum((mean[s] - mean[s].mean()) ** 2)
        results[s] = (sl, D0, np.sqrt(cov[1, 1]), r2)
        print(f"  {s:>2}: slope={sl:8.3f}  D0={D0:.4f} +- {np.sqrt(cov[1,1]):.4f}"
              f"  R2={r2:.3f}")

    sl_w = results["O"][0]
    eta = KB * T_K * XI / (6 * np.pi * abs(sl_w) * 1e-18) * 1e3
    print(f"  water-slope viscosity eta_fit = {eta:.3f} mPa*s")

    # --- shared-slope fit: ions inherit the water slope, fit intercept only ---
    print("\nshared-slope fits (slope fixed to water's):")
    for s in ["Na", "Cl"]:
        resid = mean[s] - sl_w * invL          # D0 estimate per box
        w = 1.0 / sem[s] ** 2
        D0 = np.sum(w * resid) / np.sum(w)
        D0_err = 1.0 / np.sqrt(np.sum(w))
        scatter = np.std(resid, ddof=1)
        print(f"  {s:>2}: D0={D0:.4f} +- {D0_err:.4f} (weighted)"
              f"   per-box scatter {scatter:.4f}")
        results[s + "_shared"] = (sl_w, D0, D0_err, np.nan)

    np.savez(os.path.join(OUT, "replica_D.npz"),
             Ls=Ls, species=SPECIES,
             D_mean=np.array([mean[s] for s in SPECIES]),
             D_sem=np.array([sem[s] for s in SPECIES]),
             free_fits=np.array([results[s][:3] for s in SPECIES]),
             shared_fits=np.array([results[s + "_shared"][:3] for s in ["Na", "Cl"]]))
    print(f"\nsaved {OUT}/replica_D.npz")
    print("units: D in A^2/ps (1 A^2/ps = 1e-4 cm^2/s = 1e-8 m^2/s)")


if __name__ == "__main__":
    main()
