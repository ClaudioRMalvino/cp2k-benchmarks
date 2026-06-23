#!/usr/bin/env python3
"""Fig. S4 validation: NVE energy-conservation stability of feature/nnp-chebyshev
across system sizes, plus the master==chebyshev equivalence table.

The optimisation-vs-physics story for Report 2 is:
  1. chebyshev runs STABLE production NVE MD at every system size (this figure):
     the conserved quantity is flat to ~uHa/atom over 100 ps for N=64..512.
  2. chebyshev reproduces master's energies/forces to ~1e-11 (the table): the
     polynomial-kernel optimisation did not perturb the physics.

NOTE: master's N>=128 production .ener on disk is STALE (May, generated from the
pre-fix equilibration snapshots) and is deliberately NOT plotted -- equivalence
is established by the single-point/force agreement instead (see equivalence CSV).

Usage:
  python3 plot_nve_stability.py
  python3 plot_nve_stability.py --prod-root <...>/results/figS4/production
"""
import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CAMBRIDGE = {
    "blue_warm": "#00BDB6", "blue_dark": "#133844", "crest": "#FD8153",
    "crest_dark": "#DD3025", "indigo": "#5366E0", "purple": "#A368DF",
    "green": "#4DB78C", "slate_2": "#B5BDC8", "slate_3": "#546072",
    "slate_4": "#232830",
}
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Liberation Serif", "serif"],
    "font.size": 12, "axes.labelsize": 13, "axes.titlesize": 13,
    "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 10,
    "text.color": CAMBRIDGE["slate_4"], "axes.edgecolor": CAMBRIDGE["slate_4"],
    "axes.labelcolor": CAMBRIDGE["slate_4"], "xtick.color": CAMBRIDGE["slate_4"],
    "ytick.color": CAMBRIDGE["slate_4"], "axes.linewidth": 0.8,
    "axes.grid": True, "grid.color": CAMBRIDGE["slate_2"], "grid.linestyle": "--",
    "grid.linewidth": 0.5, "grid.alpha": 0.55, "lines.linewidth": 1.5,
    "mathtext.fontset": "cm", "savefig.dpi": 300, "savefig.bbox": "tight",
    "pdf.fonttype": 42,
})

ATOMS_PER_MOL = 3
HA_TO_UHA = 1.0e6
# size -> colour (cool->warm with N), matches a perceptual size ordering
SIZE_STYLE = {
    64:  dict(color=CAMBRIDGE["blue_warm"], label="N=64"),
    128: dict(color=CAMBRIDGE["green"],     label="N=128"),
    256: dict(color=CAMBRIDGE["indigo"],    label="N=256"),
    512: dict(color=CAMBRIDGE["crest"],     label="N=512"),
    1024: dict(color=CAMBRIDGE["crest_dark"], label="N=1024"),
}

# Single-point master vs chebyshev agreement (this session's diagnostic,
# diag/n128_pbc) + the N64 physics_check.  These are constants, not regenerated
# here, so the table is self-documenting in the report.
EQUIV_ROWS = [
    # system, n_atoms, E_master_Ha, E_cheby_Ha, maxdF_Ha_bohr, source
    ("N64 (cubic)",        192, None, None, 2.56e-13, "physics_check (job 30451382)"),
    ("N128 explicit cell", 384, -90.254750953763988, -90.254750953779833, 1.00e-11, "diag/n128_pbc single-point"),
    ("N128 MUC 2x1x1",     384, -90.254750953763988, -90.254750953779833, 1.00e-11, "diag/n128_pbc single-point"),
]


def load_ener(path, stride):
    """Return (time_ps, temp_K, consqty_Ha) downsampled by `stride`."""
    t, temp, cq = [], [], []
    with open(path) as f:
        for i, line in enumerate(f):
            if line.startswith("#"):
                continue
            if (i % stride):
                continue
            p = line.split()
            if len(p) < 6:
                continue
            t.append(float(p[1]))      # Time [fs]
            temp.append(float(p[3]))   # Temp [K]
            cq.append(float(p[5]))     # Cons Qty [a.u.]
    return np.array(t) / 1000.0, np.array(temp), np.array(cq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prod-root",
                    default="/rds/user/crm98/hpc-work/cp2k-benchmarks/results/figS4/production")
    ap.add_argument("--branch", default="feature-nnp-chebyshev")
    ap.add_argument("--sizes", nargs="+", type=int, default=[64, 128, 256, 512])
    ap.add_argument("--stride", type=int, default=200, help="row downsample")
    ap.add_argument("--out-dir", default="/home/crm98/cp2k-benchmarks/plots/thesis_figures")
    ap.add_argument("--csv-dir",
                    default="/home/crm98/cp2k-benchmarks/results/figS4/analysis")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.csv_dir, exist_ok=True)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(8.4, 3.6))
    summary = [("n_molecules,n_atoms,t_total_ps,abs_drift_per_atom_uHa,"
                "drift_rate_uHa_per_atom_per_ps,T_mean_K,T_std_K")]

    for size in args.sizes:
        seg = sorted(glob.glob(os.path.join(args.prod_root, args.branch,
                                            f"N{size}", "seg1", "*-1.ener")))
        if not seg:
            print(f"  (skip N{size}: no seg1 .ener)")
            continue
        natoms = size * ATOMS_PER_MOL
        t, temp, cq = load_ener(seg[0], args.stride)
        drift = (cq - cq[0]) / natoms * HA_TO_UHA          # uHa/atom
        st = SIZE_STYLE.get(size, dict(color=CAMBRIDGE["slate_3"], label=f"N={size}"))
        axA.plot(t, drift, color=st["color"], label=st["label"])
        axB.plot(t, temp, color=st["color"], label=st["label"], alpha=0.85)

        rate = (drift[-1]) / t[-1] if t[-1] else float("nan")
        summary.append(f"{size},{natoms},{t[-1]:.2f},{drift[-1]:.4f},"
                       f"{rate:.5f},{temp.mean():.2f},{temp.std():.2f}")
        print(f"  N{size:<5d} |CQ drift|={abs(drift[-1]):.3f} uHa/atom over "
              f"{t[-1]:.0f} ps  ({rate:+.4f} uHa/atom/ps)  T={temp.mean():.1f}+/-{temp.std():.1f} K")

    axA.axhline(0, color=CAMBRIDGE["slate_3"], lw=0.8, ls=":")
    axA.set_xlabel("time  [ps]")
    axA.set_ylabel(r"conserved-qty drift  [$\mu$Ha / atom]")
    axA.set_title("(a)  NVE energy conservation", loc="left")
    axA.legend(frameon=False, ncol=2)

    axB.set_xlabel("time  [ps]")
    axB.set_ylabel("temperature  [K]")
    axB.set_title("(b)  Temperature stability", loc="left")
    axB.legend(frameon=False, ncol=2)

    fig.suptitle("feature/nnp-chebyshev: stable NVE production MD across system sizes",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(args.out_dir, f"figS4_nve_stability.{ext}"))
    print(f"\nwrote {args.out_dir}/figS4_nve_stability.png/.pdf")

    with open(os.path.join(args.csv_dir, "nve_stability_summary.csv"), "w") as f:
        f.write("\n".join(summary) + "\n")

    # equivalence table CSV
    eq = ["system,n_atoms,E_master_Ha,E_chebyshev_Ha,abs_dE_Ha,abs_dE_per_atom_Ha,max_dF_Ha_per_bohr,source"]
    for sysname, na, em, ec, mdf, src in EQUIV_ROWS:
        if em is None:
            eq.append(f"{sysname},{na},,,,,{mdf:.3e},{src}")
        else:
            de = abs(em - ec)
            eq.append(f"{sysname},{na},{em:.15f},{ec:.15f},{de:.3e},{de/na:.3e},{mdf:.3e},{src}")
    with open(os.path.join(args.csv_dir, "master_chebyshev_equivalence.csv"), "w") as f:
        f.write("\n".join(eq) + "\n")
    print(f"wrote {args.csv_dir}/nve_stability_summary.csv")
    print(f"wrote {args.csv_dir}/master_chebyshev_equivalence.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
