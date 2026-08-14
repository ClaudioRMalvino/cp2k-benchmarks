#!/usr/bin/env python3
"""Diffusion vs NaCl concentration (O'Neill SI Fig S10 style).

Takes molality=summary.csv pairs (all the SAME cell size, so finite-size effects
cancel in the ratio); plots per-species YH-corrected D_0 vs molality and water
D(c)/D(c=0) (the quantity the O'Neill SI compares to experiment).
Example: plot_conc.py --cell 211 0.0=np0/analysis/nacl_summary.csv 0.90=np1/analysis/nacl_summary.csv
"""
import argparse
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

STYLE = {
    "O":  dict(color="#133844", marker="o", label=r"H$_2$O (O)"),
    "Na": dict(color="#5366E0", marker="^", label=r"Na$^+$"),
    "Cl": dict(color="#4DB78C", marker="s", label=r"Cl$^-$"),
}


def load(summary, cell, species):
    with open(summary) as f:
        for r in csv.DictReader(f):
            if r["cell"] == cell and r["species"] == species:
                return float(r["D_0_ang2_ps"]), float(r["D_0_sem"])
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("series", nargs="+",
                    help="molality=path/to/nacl_summary.csv")
    ap.add_argument("--cell", default="211",
                    help="cell label to read from each summary")
    ap.add_argument("--out", default="nacl_D_vs_concentration.png")
    args = ap.parse_args()

    pts = {sp: [] for sp in STYLE}
    for item in args.series:
        mol, _, path = item.partition("=")
        for sp in STYLE:
            v = load(path, args.cell, sp)
            if v is not None and np.isfinite(v[0]):
                pts[sp].append((float(mol), v[0], v[1]))

    fig, (ax, axr) = plt.subplots(1, 2, figsize=(11, 4.5))
    for sp, st in STYLE.items():
        if not pts[sp]:
            continue
        arr = np.array(sorted(pts[sp]))
        ax.errorbar(arr[:, 0], arr[:, 1], yerr=arr[:, 2],
                    capsize=3, **st)
    ax.set_xlabel("NaCl molality [mol/kg]")
    ax.set_ylabel(r"$D_0$ [$\mathrm{\AA}^2$/ps]")
    ax.set_title(f"Diffusion vs concentration (cell{args.cell})")
    ax.legend()

    if pts["O"]:
        arr = np.array(sorted(pts["O"]))
        ref_rows = arr[arr[:, 0] == arr[:, 0].min()]
        d_ref = ref_rows[0, 1]
        rel = arr[:, 1] / d_ref
        rel_err = rel * np.sqrt((arr[:, 2] / arr[:, 1]) ** 2
                                + (ref_rows[0, 2] / d_ref) ** 2)
        axr.errorbar(arr[:, 0], rel, yerr=rel_err, capsize=3,
                     **STYLE["O"])
        axr.axhline(1.0, color="gray", lw=0.8, ls=":")
        axr.set_xlabel("NaCl molality [mol/kg]")
        axr.set_ylabel(r"$D_{\mathrm{H_2O}}(c)\,/\,D_{\mathrm{H_2O}}(0)$")
        axr.set_title("Water diffusion relative to pure water")

    fig.tight_layout()
    fig.savefig(args.out, dpi=300)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
