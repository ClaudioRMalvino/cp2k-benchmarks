#!/usr/bin/env python3
"""MSD vs VACF estimator cross-check for the Madrid-2019 replica state points.

Both estimate the same D_PBC at each box size:
  MSD  : multi-origin Einstein slope (replicas/analysis/replica_D.npz, NVT)
  VACF : Green-Kubo running-integral plateau (vacf_replica/L*_s*/, NVE)
Green-Kubo and Einstein are mathematically equivalent, so agreement within
seed error bars validates both pipelines (supervisor's cross-check).

Plateau = mean of D(t) over --plat-min..--plat-max ps (VACF decays in <2 ps;
the tail of the running integral is the diffusive limit).
"""
import argparse
import glob
import os
import re

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vacf-root",
                    default="/data/cerberus1/crm98/nacl_diffusion/lammps_madrid/vacf_replica")
    ap.add_argument("--replica-npz",
                    default="/data/cerberus1/crm98/nacl_diffusion/lammps_madrid/replicas/analysis/replica_D.npz")
    ap.add_argument("--plat-min", type=float, default=3.0)
    ap.add_argument("--plat-max", type=float, default=5.0)
    args = ap.parse_args()

    z = np.load(args.replica_npz, allow_pickle=True)
    Ls, species = list(z["Ls"]), list(z["species"])
    msd_mean, msd_sem = z["D_mean"], z["D_sem"]

    # collect VACF plateau D per (L, seed, species)
    vacf = {}   # (L, sp) -> [D per seed]
    for f in sorted(glob.glob(os.path.join(args.vacf_root, "L*_s*", "*.vacf"))):
        L = float(re.search(r"L([\d.]+)_s", f).group(1))
        dat = np.loadtxt(f, comments="#")
        with open(f) as fh:
            names = fh.readline().lstrip("# ").split()
        t = dat[:, names.index("t_ps")]
        m = (t >= args.plat_min) & (t <= args.plat_max)
        for sp in species:
            col = f"D_{sp}"
            if col in names:
                vacf.setdefault((L, sp), []).append(dat[m, names.index(col)].mean())

    nseed = max((len(v) for v in vacf.values()), default=0)
    print(f"VACF runs found: {len(vacf)//max(len(species),1)} boxes x {nseed} seeds; "
          f"plateau window {args.plat_min}-{args.plat_max} ps")
    print(f"\n{'L':>7} {'sp':>3} | {'D_MSD':>7} {'+-':>6} | {'D_VACF':>7} {'+-':>6} | "
          f"{'ratio':>6} {'sigma':>6}")
    for j, L in enumerate(Ls):
        for i, sp in enumerate(species):
            v = np.array(vacf.get((L, sp), []))
            if v.size == 0:
                continue
            vm, vs = v.mean(), v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else np.nan
            mm, ms = msd_mean[i, j], msd_sem[i, j]
            pull = (vm - mm) / np.hypot(vs, ms)
            print(f"{L:7.2f} {sp:>3} | {mm:7.4f} {ms:6.4f} | {vm:7.4f} {vs:6.4f} | "
                  f"{vm/mm:6.3f} {pull:+6.1f}")
    print("\nratio = VACF/MSD (expect ~1); sigma = discrepancy / combined error "
          "(|sigma| < ~2 means the estimators agree within noise)")


if __name__ == "__main__":
    main()
