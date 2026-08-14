#!/usr/bin/env python3
"""Aggregate NaCl(aq) NVE production segments into diffusion/viscosity CSVs.

Adapted from figS4/aggregate_figS4.py: cells labelled by tiling ("111", ...) of
the 12.42 A base cell; diffusion per species (O/Na/Cl). Reuses the figS4
compute_viscosity.py / compute_diffusion.py unchanged (Yeh-Hummer:
D_0 = D_PBC + kB*T*xi/(6*pi*eta*L), xi = 2.837297, T = 300 K, L = V^(1/3)).
Usage: aggregate_nacl.py --prod-root <.../production> --out-dir <.../analysis>
"""
import argparse
import glob
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS4 = os.path.normpath(os.path.join(
    HERE, "..", "CSD3_benchmark_scripts", "figS4"))
VISC = os.path.join(FIGS4, "compute_viscosity.py")
DIFF = os.path.join(FIGS4, "compute_diffusion.py")
BASE_L = 12.42
SPECIES = ["O", "Na", "Cl"]


def read_header(path):
    meta = {}
    with open(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            s = line.lstrip("#").strip()
            if ":" in s:
                k, _, v = s.partition(":")
                meta[k.strip()] = v.strip()
    return meta


def read_curve(path, ncols):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.strip().split(",")
            if len(p) != ncols:
                continue
            try:
                rows.append([float(x) for x in p])
            except ValueError:
                continue
    return np.asarray(rows)


def stack_mean_sem(curves, ycol):
    n = min(len(c) for c in curves)
    x = curves[0][:n, 0]
    ys = np.array([c[:n, ycol] for c in curves])
    mean = ys.mean(axis=0)
    sem = ys.std(axis=0, ddof=1) / np.sqrt(len(ys)) if len(ys) > 1 else np.zeros(n)
    return x, mean, sem


def msem(a):
    a = np.asarray(a, float)
    if len(a) == 0:
        return float("nan"), float("nan")
    m = a.mean()
    s = a.std(ddof=1) / np.sqrt(len(a)) if len(a) > 1 else 0.0
    return m, s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prod-root", required=True,
                    help=".../nacl_diffusion/<model>/production")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--cells", nargs="+", default=["111", "211", "221", "222"])
    ap.add_argument("--discard-ps", type=float, default=10.0)
    ap.add_argument("--int-limit-ps", type=float, default=3.0)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    summary = ["cell,species,L_ang,n_segments,eta_mPa_s,eta_sem,"
               "D_PBC_ang2_ps,D_PBC_sem,D_0_ang2_ps,D_0_sem"]

    for cell_label in args.cells:
        mult = tuple(int(c) for c in cell_label)
        cell = [BASE_L * m for m in mult]
        volume = float(np.prod(cell))
        L_eff = volume ** (1.0 / 3.0)
        seg_dirs = sorted(glob.glob(
            os.path.join(args.prod_root, f"cell{cell_label}", "seg*")))

        etas, acf_curves, reta_curves = [], [], []
        dpbc = {s: [] for s in SPECIES}
        d0 = {s: [] for s in SPECIES}
        msd_curves = {s: [] for s in SPECIES}

        for sd in seg_dirs:
            stress = glob.glob(os.path.join(sd, "*-1.stress"))
            traj = glob.glob(os.path.join(sd, "*-pos-1.xyz"))
            if not stress or not traj:
                print(f"  WARN cell{cell_label} {os.path.basename(sd)}: "
                      f"missing stress/traj, skipping", file=sys.stderr)
                continue

            vcsv = os.path.join(sd, "viscosity.csv")
            subprocess.run([sys.executable, VISC, stress[0],
                            "--n-molecules", "1",
                            "--volume-ang3", f"{volume:.6f}",
                            "--discard-ps", str(args.discard_ps),
                            "--int-limit-ps", str(args.int_limit_ps),
                            "--out", vcsv], check=True)
            eta = float(read_header(vcsv)["eta_mPa_s"])
            etas.append(eta)
            acf_curves.append(read_curve(vcsv, 3)[:, [0, 1]])
            reta_curves.append(read_curve(vcsv, 3)[:, [0, 2]])

            for sp in SPECIES:
                dcsv = os.path.join(sd, f"diffusion_{sp}.csv")
                subprocess.run(
                    [sys.executable, DIFF, traj[0],
                     "--cell-ang", f"{cell[0]}", f"{cell[1]}", f"{cell[2]}",
                     "--viscosity-mpa-s", f"{eta}",
                     "--species", sp,
                     "--out", dcsv], check=True)
                dmeta = read_header(dcsv)
                dpbc[sp].append(float(dmeta["D_PBC_ang2_ps"]))
                d0[sp].append(float(dmeta["D_0_ang2_ps"]))
                msd_curves[sp].append(read_curve(dcsv, 2))

        n = len(etas)
        if n == 0:
            print(f"  cell{cell_label}: no usable segments", file=sys.stderr)
            continue

        e_m, e_s = msem(etas)
        for sp in SPECIES:
            dp_m, dp_s = msem(dpbc[sp])
            d0_m, d0_s = msem(d0[sp])
            summary.append(f"{cell_label},{sp},{L_eff:.4f},{n},"
                           f"{e_m:.6f},{e_s:.6f},{dp_m:.8f},{dp_s:.8f},"
                           f"{d0_m:.8f},{d0_s:.8f}")
            print(f"  cell{cell_label} {sp:2s}  n={n}  eta={e_m:.3f} mPa.s  "
                  f"D_PBC={dp_m:.4f}  D_0={d0_m:.4f} A^2/ps")
            x, m, s = stack_mean_sem(msd_curves[sp], 1)
            np.savetxt(os.path.join(args.out_dir,
                                    f"msd_{sp}_cell{cell_label}.csv"),
                       np.c_[x, m, s], delimiter=",",
                       header="lag_ps,msd_mean,msd_sem", comments="")

        x, m, s = stack_mean_sem(acf_curves, 1)
        np.savetxt(os.path.join(args.out_dir, f"acf_cell{cell_label}.csv"),
                   np.c_[x, m, s], delimiter=",",
                   header="lag_fs,acf_norm_mean,acf_norm_sem", comments="")
        x, m, s = stack_mean_sem(reta_curves, 1)
        np.savetxt(os.path.join(args.out_dir,
                                f"running_eta_cell{cell_label}.csv"),
                   np.c_[x, m, s], delimiter=",",
                   header="lag_fs,running_eta_mean,running_eta_sem",
                   comments="")

    out = os.path.join(args.out_dir, "nacl_summary.csv")
    with open(out, "w") as f:
        f.write("\n".join(summary) + "\n")
    print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
