#!/usr/bin/env python3
"""Conductivity vs concentration for the Madrid-2019 NaCl(aq) series.

Reads <ROOT>/m<MOLALITY>/L<L>_s<SEED>/<tag>.onsager (from onsager_multiorigin.py)
plus <tag>.final.data (box volume + ion counts), and computes per run:

  D_i        self-diffusion per species (Einstein fit)
  kappa_NE   Nernst-Einstein: e^2 z^2 (N_Na D_Na + N_Cl D_Cl) / (kB T V)
  kappa_Ons  Einstein-Helfand: e^2 z^2 slope(C_NaNa + C_ClCl - 2 C_NaCl)/(6 kB T V)
  t_Na       (slope C_NaNa - slope C_NaCl) / slope(C_sum)   [charge-free]
  Delta_NE   1 - kappa_Ons/kappa_NE  (ion-correlation deviation)

Both charge conventions reported (displacement correlations are charge-free):
z = 0.85 (model) vs z = 1 (formal) differ by exactly 0.7225. Mean +- SEM over seeds.
"""
import argparse
import glob
import os
import re

import numpy as np

E2 = (1.602176634e-19) ** 2          # C^2
KBT = 1.380649e-23 * 298.15          # J
A2PS_TO_M2S = 1.0e-8                 # A^2/ps -> m^2/s
A3_TO_M3 = 1.0e-30


def box_and_counts(final_data):
    """Volume (A^3) and (N_Na, N_Cl) from a LAMMPS data file (atom_style full)."""
    lo_hi, n_atoms = [], None
    with open(final_data) as fh:
        lines = fh.readlines()
    for i, ln in enumerate(lines):
        if re.search(r"^\s*\d+\s+atoms\s*$", ln):
            n_atoms = int(ln.split()[0])
        if re.search(r"[xyz]lo .*[xyz]hi", ln):
            a, b = map(float, ln.split()[:2])
            lo_hi.append(b - a)
        if ln.strip().startswith("Atoms"):
            start = i + 2
            break
    types = np.array([int(lines[start + j].split()[2]) for j in range(n_atoms)])
    V = lo_hi[0] * lo_hi[1] * lo_hi[2]
    return V, int((types == 3).sum()), int((types == 4).sum())


def slope(t, y, tmin, tmax):
    m = (t >= tmin) & (t <= tmax)
    return np.polyfit(t[m], y[m], 1)[0]


def analyze_run(d, tmin, tmax):
    tag = os.path.basename(d).replace("_", "_")  # L<L>_s<seed>; tag has m too
    mo = re.search(r"L([\d.]+)_s(\d+)$", d)
    L, seed = mo.group(1), mo.group(2)
    mol = re.search(r"/m([\d.]+)/", d + "/").group(1)
    tag = f"L{L}_m{mol}_s{seed}"
    ons = np.loadtxt(os.path.join(d, f"{tag}.onsager"), comments="#")
    V, nNa, nCl = box_and_counts(os.path.join(d, f"{tag}.final.data"))
    t = ons[:, 0]
    # columns: tau msd_O msd_Na msd_Cl coll_NaNa coll_ClCl coll_NaCl
    D = {s: slope(t, ons[:, i + 1], tmin, tmax) / 6.0
         for i, s in enumerate(["O", "Na", "Cl"])}
    sNaNa = slope(t, ons[:, 4], tmin, tmax)
    sClCl = slope(t, ons[:, 5], tmin, tmax)
    sNaCl = slope(t, ons[:, 6], tmin, tmax)
    s_sum = sNaNa + sClCl - 2.0 * sNaCl

    pref = E2 / (KBT * V * A3_TO_M3) * A2PS_TO_M2S      # z=1 prefactor, S/m
    kNE = pref * (nNa * D["Na"] + nCl * D["Cl"])
    kOns = pref * s_sum / 6.0
    tNa = (sNaNa - sNaCl) / s_sum if s_sum != 0 else np.nan
    return dict(L=float(L), mol=float(mol), seed=int(seed), V=V,
                nNa=nNa, nCl=nCl, D=D, kNE=kNE, kOns=kOns, tNa=tNa,
                sNaNa=sNaNa, sClCl=sClCl, sNaCl=sNaCl)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root",
                    default="/data/cerberus1/crm98/nacl_diffusion/lammps_madrid/conductivity")
    ap.add_argument("--tmin", type=float, default=20.0)
    ap.add_argument("--tmax", type=float, default=400.0)
    args = ap.parse_args()

    runs = []
    for d in sorted(glob.glob(os.path.join(args.root, "m*", "L*_s*"))):
        mo = re.search(r"L([\d.]+)_s(\d+)$", d)
        mol = re.search(r"/m([\d.]+)/", d + "/").group(1)
        tag = f"L{mo.group(1)}_m{mol}_s{mo.group(2)}"
        if not os.path.exists(os.path.join(d, f"{tag}.onsager")):
            continue
        runs.append(analyze_run(d, args.tmin, args.tmax))
    if not runs:
        print("no completed runs found"); return

    print(f"{len(runs)} runs, fit window {args.tmin}-{args.tmax} ps")
    print("kappa in S/m; z=1 (formal) column, and z=0.85 (model) = x0.7225")
    print(f"\n{'m':>5} {'L':>7} {'#s':>3} | {'kNE(z1)':>8} {'kOns(z1)':>9} "
          f"{'+-':>6} | {'kOns(.85)':>9} | {'t_Na':>6} {'+-':>6} | {'1-Ons/NE':>8}")
    combos = sorted({(r["mol"], r["L"]) for r in runs})
    agg = {}
    for mol, L in combos:
        sel = [r for r in runs if r["mol"] == mol and r["L"] == L]
        kNE = np.array([r["kNE"] for r in sel])
        kOns = np.array([r["kOns"] for r in sel])
        tNa = np.array([r["tNa"] for r in sel])
        n = len(sel)
        sem = lambda a: np.std(a, ddof=1) / np.sqrt(n) if n > 1 else np.nan
        agg[(mol, L)] = (kNE.mean(), kOns.mean(), sem(kOns), tNa.mean(), sem(tNa))
        print(f"{mol:5.2f} {L:7.2f} {n:3d} | {kNE.mean():8.3f} {kOns.mean():9.3f} "
              f"{sem(kOns):6.3f} | {kOns.mean()*0.7225:9.3f} | "
              f"{tNa.mean():6.3f} {sem(tNa):6.3f} | {1-kOns.mean()/kNE.mean():8.3f}")

    # Onsager decomposition: kappa = k_NaNa + k_ClCl - 2 k_NaCl, k_ss' = e^2 z^2
    # slope(C_ss')/(6 kB T V); diagonals split into self (6 N_i D_i) + distinct.
    # Positive distinct NaCl (co-motion, e.g. pairing) REDUCES kappa. S/m at z=1.
    print(f"\nOnsager decomposition (z=1, S/m): kappa = kNaNa + kClCl - 2 kNaCl")
    print(f"{'m':>5} {'L':>7} | {'kNaNa':>7} {'kClCl':>7} {'-2kNaCl':>8} | "
          f"{'dNaNa':>7} {'dClCl':>7} {'dNaCl':>7} | {'2kNaCl/kNE':>10}")
    for mol, L in combos:
        sel = [r for r in runs if r["mol"] == mol and r["L"] == L]
        pref = lambda r: E2 / (KBT * r["V"] * A3_TO_M3) * A2PS_TO_M2S / 6.0
        kNaNa = np.mean([pref(r) * r["sNaNa"] for r in sel])
        kClCl = np.mean([pref(r) * r["sClCl"] for r in sel])
        kNaCl = np.mean([pref(r) * r["sNaCl"] for r in sel])
        dNaNa = np.mean([pref(r) * (r["sNaNa"] - 6 * r["nNa"] * r["D"]["Na"])
                         for r in sel])
        dClCl = np.mean([pref(r) * (r["sClCl"] - 6 * r["nCl"] * r["D"]["Cl"])
                         for r in sel])
        kNE_m = np.mean([r["kNE"] for r in sel])
        print(f"{mol:5.2f} {L:7.2f} | {kNaNa:7.3f} {kClCl:7.3f} {-2*kNaCl:8.3f} | "
              f"{dNaNa:7.3f} {dClCl:7.3f} {kNaCl:7.3f} | {2*kNaCl/kNE_m:10.3f}")

    out = os.path.join(args.root, "analysis")
    os.makedirs(out, exist_ok=True)
    np.savez(os.path.join(out, "kappa_vs_c.npz"),
             combos=np.array(combos),
             table=np.array([agg[c] for c in combos]),
             runs_mol=np.array([r["mol"] for r in runs]),
             runs_L=np.array([r["L"] for r in runs]),
             runs_kNE=np.array([r["kNE"] for r in runs]),
             runs_kOns=np.array([r["kOns"] for r in runs]),
             runs_tNa=np.array([r["tNa"] for r in runs]),
             runs_DNa=np.array([r["D"]["Na"] for r in runs]),
             runs_DCl=np.array([r["D"]["Cl"] for r in runs]),
             runs_DO=np.array([r["D"]["O"] for r in runs]),
             runs_seed=np.array([r["seed"] for r in runs]),
             runs_V=np.array([r["V"] for r in runs]),
             runs_nNa=np.array([r["nNa"] for r in runs]),
             runs_nCl=np.array([r["nCl"] for r in runs]),
             runs_sNaNa=np.array([r["sNaNa"] for r in runs]),
             runs_sClCl=np.array([r["sClCl"] for r in runs]),
             runs_sNaCl=np.array([r["sNaCl"] for r in runs]))
    print(f"\nsaved {out}/kappa_vs_c.npz")


if __name__ == "__main__":
    main()
