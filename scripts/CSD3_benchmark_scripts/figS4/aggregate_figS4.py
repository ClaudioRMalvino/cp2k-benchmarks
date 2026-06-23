#!/usr/bin/env python3
"""Aggregate Fig. S4 NVE production segments into panel-ready CSVs."""
import argparse
import glob
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
VISC = os.path.join(HERE, "compute_viscosity.py")
DIFF = os.path.join(HERE, "compute_diffusion.py")
BASE_L = 12.42

MULT = {64: (1, 1, 1), 128: (2, 1, 1), 256: (2, 2, 1),
        512: (2, 2, 2), 1024: (4, 2, 2)}


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prod-root", required=True,
                    help="<...>/results/figS4/production")
    ap.add_argument("--out-dir", required=True)
    # Report 2 validates chebyshev only: master's N>=128 production is stale
    # (May, pre-fix snapshots) and master==chebyshev physics is established by
    # the single-point/force equivalence table instead, not by re-run MD.
    ap.add_argument("--branches", nargs="+",
                    default=["feature-nnp-chebyshev"])
    ap.add_argument("--sizes", nargs="+", type=int,
                    default=[64, 128, 256, 512, 1024])
    ap.add_argument("--discard-ps", type=float, default=10.0)
    ap.add_argument("--int-limit-ps", type=float, default=3.0)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    summary = [("branch,n_molecules,n_segments,"
                "eta_mPa_s,eta_sem,D_PBC_ang2_ps,D_PBC_sem,"
                "D_0_ang2_ps,D_0_sem,time_per_step_s,time_per_step_sem,"
                "seg_walltime_s,seg_walltime_sem")]

    for branch in args.branches:
        for size in args.sizes:
            mult = MULT[size]
            cell = [BASE_L * m for m in mult]
            volume = float(np.prod(cell))
            seg_dirs = sorted(glob.glob(
                os.path.join(args.prod_root, branch, f"N{size}", "seg*")))
            etas, dpbc, d0, tps, wts = [], [], [], [], []
            acf_curves, reta_curves, msd_curves = [], [], []

            for sd in seg_dirs:
                stress = glob.glob(os.path.join(sd, "*-1.stress"))
                traj = glob.glob(os.path.join(sd, "*-pos-1.xyz"))
                timing = os.path.join(sd, "timing.csv")
                if not stress or not traj:
                    print(f"  WARN {branch} N{size} {os.path.basename(sd)}: "
                          f"missing stress/traj, skipping", file=sys.stderr)
                    continue

                vcsv = os.path.join(sd, "viscosity.csv")
                subprocess.run([sys.executable, VISC, stress[0],
                                "--n-molecules", str(size),
                                "--volume-ang3", f"{volume:.6f}",
                                "--discard-ps", str(args.discard_ps),
                                "--int-limit-ps", str(args.int_limit_ps),
                                "--out", vcsv], check=True)
                vmeta = read_header(vcsv)
                eta = float(vmeta["eta_mPa_s"])
                etas.append(eta)
                acf_curves.append(read_curve(vcsv, 3)[:, [0, 1]])
                reta_curves.append(read_curve(vcsv, 3)[:, [0, 2]])

                dcsv = os.path.join(sd, "diffusion.csv")
                subprocess.run([sys.executable, DIFF, traj[0],
                                "--cell-ang", f"{cell[0]}", f"{cell[1]}", f"{cell[2]}",
                                "--viscosity-mpa-s", f"{eta}",
                                "--out", dcsv], check=True)
                dmeta = read_header(dcsv)
                dpbc.append(float(dmeta["D_PBC_ang2_ps"]))
                d0.append(float(dmeta["D_0_ang2_ps"]))
                msd_curves.append(read_curve(dcsv, 2))

                if os.path.exists(timing):
                    with open(timing) as f:
                        for line in f:
                            p = line.strip().split(",")
                            if len(p) >= 9 and p[0] == branch:
                                try:
                                    tps.append(float(p[6]))
                                    wts.append(float(p[8]))
                                except ValueError:
                                    pass

            n = len(etas)
            if n == 0:
                print(f"  {branch} N{size}: no usable segments", file=sys.stderr)
                continue

            def msem(a):
                a = np.asarray(a, float)
                m = a.mean()
                s = a.std(ddof=1) / np.sqrt(len(a)) if len(a) > 1 else 0.0
                return m, s

            e_m, e_s = msem(etas)
            dp_m, dp_s = msem(dpbc)
            d0_m, d0_s = msem(d0)
            t_m, t_s = msem(tps) if tps else (float("nan"), float("nan"))
            w_m, w_s = msem(wts) if wts else (float("nan"), float("nan"))
            summary.append(f"{branch},{size},{n},"
                           f"{e_m:.6f},{e_s:.6f},{dp_m:.8f},{dp_s:.8f},"
                           f"{d0_m:.8f},{d0_s:.8f},{t_m:.6f},{t_s:.6f},"
                           f"{w_m:.2f},{w_s:.2f}")
            print(f"  {branch:28s} N{size:<5d} n={n}  "
                  f"eta={e_m:.3f}+/-{e_s:.3f} mPa.s  "
                  f"D_PBC={dp_m:.4f}  D_0={d0_m:.4f} A^2/ps  "
                  f"t/step={t_m:.4f}s")

            x, m, s = stack_mean_sem(acf_curves, 1)
            np.savetxt(os.path.join(args.out_dir, f"acf_{branch}_N{size}.csv"),
                       np.c_[x, m, s], delimiter=",",
                       header="lag_fs,acf_norm_mean,acf_norm_sem", comments="")
            x, m, s = stack_mean_sem(reta_curves, 1)
            np.savetxt(os.path.join(args.out_dir, f"running_eta_{branch}_N{size}.csv"),
                       np.c_[x, m, s], delimiter=",",
                       header="lag_fs,running_eta_mean,running_eta_sem", comments="")
            x, m, s = stack_mean_sem(msd_curves, 1)
            np.savetxt(os.path.join(args.out_dir, f"msd_{branch}_N{size}.csv"),
                       np.c_[x, m, s], delimiter=",",
                       header="lag_ps,msd_mean,msd_sem", comments="")

    with open(os.path.join(args.out_dir, "figS4_summary.csv"), "w") as f:
        f.write("\n".join(summary) + "\n")
    print(f"\nWrote {args.out_dir}/figS4_summary.csv and per-(branch,size) curves")
    return 0


if __name__ == "__main__":
    sys.exit(main())
