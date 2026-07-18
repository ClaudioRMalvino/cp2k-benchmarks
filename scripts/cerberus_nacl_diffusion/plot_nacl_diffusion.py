#!/usr/bin/env python3
"""Yeh-Hummer finite-size plot for NaCl(aq) diffusion (Morawietz-figS4 style).

Reads nacl_summary.csv from aggregate_nacl.py and produces
nacl_diffusion_replication.png:
  left panel  - D_PBC vs 1/L per species with linear fits extrapolating to
                1/L = 0, plus the Yeh-Hummer-corrected D_0 points (which
                should scatter around a horizontal line if the correction
                works for every species, since the YH slope
                -kB*T*xi/(6*pi*eta) is species-independent).
  right panel - Green-Kubo viscosity vs box size (sanity check).
"""
import argparse
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

STYLE = {
    "O":  dict(color="#133844", marker="o", label=r"H$_2$O (O)"),
    "Na": dict(color="#5366E0", marker="^", label=r"Na$^+$"),
    "Cl": dict(color="#4DB78C", marker="s", label=r"Cl$^-$"),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", required=True, help="nacl_summary.csv")
    ap.add_argument("--out", default="nacl_diffusion_replication.png")
    args = ap.parse_args()

    rows = []
    with open(args.summary) as f:
        for r in csv.DictReader(f):
            rows.append(r)

    fig, (ax, axv) = plt.subplots(1, 2, figsize=(11, 4.5))

    for sp, st in STYLE.items():
        sub = [r for r in rows if r["species"] == sp]
        if not sub:
            continue
        invL = np.array([1.0 / float(r["L_ang"]) for r in sub])
        dpbc = np.array([float(r["D_PBC_ang2_ps"]) for r in sub])
        dpbc_e = np.array([float(r["D_PBC_sem"]) for r in sub])
        d0 = np.array([float(r["D_0_ang2_ps"]) for r in sub])
        d0_e = np.array([float(r["D_0_sem"]) for r in sub])

        ax.errorbar(invL, dpbc, yerr=dpbc_e, ls="none", capsize=3,
                    color=st["color"], marker=st["marker"],
                    label=st["label"] + r" $D_\mathrm{PBC}$")
        ax.errorbar(invL, d0, yerr=d0_e, ls="none", capsize=3, mfc="white",
                    color=st["color"], marker=st["marker"], alpha=0.75,
                    label=st["label"] + r" $D_0$ (YH)")
        if len(invL) >= 2:
            b, a = np.polyfit(invL, dpbc, 1)
            xs = np.linspace(0.0, invL.max() * 1.05, 50)
            ax.plot(xs, a + b * xs, color=st["color"], lw=1.0, ls="--")
            ax.scatter([0.0], [a], color=st["color"], marker="*", s=90,
                       zorder=5)

    ax.set_xlabel(r"$1/L$ [$\mathrm{\AA}^{-1}$]")
    ax.set_ylabel(r"$D$ [$\mathrm{\AA}^2$/ps]")
    ax.set_title("NaCl(aq) diffusion: finite-size extrapolation")
    ax.set_xlim(left=0.0)
    ax.legend(fontsize=8, ncol=1)

    seen = {}
    for r in rows:
        seen[r["cell"]] = (float(r["L_ang"]),
                           float(r["eta_mPa_s"]), float(r["eta_sem"]))
    Ls = sorted(v[0] for v in seen.values())
    etas = [v[1] for v in sorted(seen.values())]
    errs = [v[2] for v in sorted(seen.values())]
    axv.errorbar(Ls, etas, yerr=errs, marker="o", capsize=3, color="#DD3025")
    axv.set_xlabel(r"$L$ [$\mathrm{\AA}$]")
    axv.set_ylabel(r"$\eta$ [mPa$\cdot$s]")
    axv.set_title("Green-Kubo shear viscosity")

    fig.tight_layout()
    fig.savefig(args.out, dpi=300)
    print(f"wrote {os.path.abspath(args.out)}")


if __name__ == "__main__":
    main()
