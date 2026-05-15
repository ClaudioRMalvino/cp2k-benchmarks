#!/usr/bin/env python3
"""
Plot the Fig. S4 replication (Morawietz et al., PNAS 2016) -- system-size
dependence of the shear viscosity and self-diffusion coefficient of liquid
water, comparing two CP2K builds (upstream master vs feature/nnp-native-spline)
over a 64 -> 1024 H2O size sweep with the reconstructed RPBE-vdW NNP.

Reads the aggregated CSVs from aggregate_figS4.py and writes:
  figS4_replication.png  4 panels: A stress ACF, B running Green-Kubo
                         viscosity, C eta vs 1/L, D D_PBC & Yeh-Hummer D_0
                         vs 1/L -- both branches overlaid.
  figS4_performance.png  master vs native-spline: wall-time per MD step,
                         wall-time per production segment, and the speedup.
  figS4_accuracy.png     eta and D_0 master vs native-spline side by side --
                         the "do the two branches give comparable physics?"
                         check (they share identical NVE starting configs).

Usage:
  plot_figS4.py [--analysis-dir DIR] [--plot-dir DIR]
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE_L = 12.42
MULT = {64: (1, 1, 1), 128: (2, 1, 1), 256: (2, 2, 1),
        512: (2, 2, 2), 1024: (4, 2, 2)}
SIZES = [64, 128, 256, 512, 1024]

BRANCHES = {
    "master":                    dict(label="upstream master",  color="tab:blue",  marker="o"),
    "feature-nnp-native-spline": dict(label="nnp-native-spline", color="tab:green", marker="^"),
}


def inv_L(size):
    """1/L with L = V^(1/3); V from the tiled cell."""
    cell = [BASE_L * m for m in MULT[size]]
    return 1.0 / (np.prod(cell) ** (1.0 / 3.0))


def load_summary(path):
    """figS4_summary.csv -> {branch: {size: row-dict}}."""
    out = {}
    with open(path) as f:
        header = f.readline().strip().split(",")
        for line in f:
            p = line.strip().split(",")
            if len(p) != len(header):
                continue
            row = dict(zip(header, p))
            b, s = row["branch"], int(row["n_molecules"])
            out.setdefault(b, {})[s] = {k: (v if k == "branch" else float(v))
                                        for k, v in row.items()}
    return out


def load_curve(path):
    if not os.path.exists(path):
        return None
    return np.genfromtxt(path, delimiter=",", names=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--analysis-dir",
                    default=os.path.expanduser("~/cp2k-benchmarks/results/figS4/analysis"))
    ap.add_argument("--plot-dir",
                    default=os.path.expanduser("~/cp2k-benchmarks/plots"))
    args = ap.parse_args()
    os.makedirs(args.plot_dir, exist_ok=True)
    A = args.analysis_dir

    summ = load_summary(os.path.join(A, "figS4_summary.csv"))

    # =====================================================================
    # Figure 1 -- the Fig. S4 replication (4 panels), both branches overlaid
    # =====================================================================
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle("Fig. S4 replication -- size dependence of viscosity & diffusion\n"
                 "RPBE-vdW NNP (Morawietz 2016), CP2K master vs nnp-native-spline",
                 fontsize=12)

    # -- Panel A: normalised stress ACF, N = 64 ---------------------------
    axA = ax[0, 0]
    for b, st in BRANCHES.items():
        c = load_curve(os.path.join(A, f"acf_{b}_N64.csv"))
        if c is None:
            continue
        axA.plot(c["lag_fs"] / 1000.0, c["acf_norm_mean"],
                 color=st["color"], label=st["label"], lw=1.5)
    axA.axhline(0, color="k", lw=0.6)
    axA.set_xlabel("lag time (ps)")
    axA.set_ylabel("normalised stress ACF")
    axA.set_title("A  stress autocorrelation (N = 64)")
    axA.legend(fontsize=8)

    # -- Panel B: running viscosity, N = 64/512/1024 (branch=colour, size=line)
    axB = ax[0, 1]
    for size, ls in zip((64, 512, 1024), ("-", "--", ":")):
        for b, st in BRANCHES.items():
            c = load_curve(os.path.join(A, f"running_eta_{b}_N{size}.csv"))
            if c is None:
                continue
            axB.plot(c["lag_fs"] / 1000.0, c["running_eta_mean"],
                     color=st["color"], ls=ls, lw=1.4,
                     label=f"{st['label']}, N={size}")
    axB.set_xlabel("Green-Kubo integration limit (ps)")
    axB.set_ylabel(r"running viscosity $\eta$ (mPa$\cdot$s)")
    axB.set_title("B  running Green-Kubo viscosity")
    axB.legend(fontsize=7)

    # -- Panel C: eta vs 1/L ----------------------------------------------
    axC = ax[1, 0]
    for b, st in BRANCHES.items():
        if b not in summ:
            continue
        ss = [s for s in SIZES if s in summ[b]]
        axC.errorbar([inv_L(s) for s in ss],
                     [summ[b][s]["eta_mPa_s"] for s in ss],
                     yerr=[summ[b][s]["eta_sem"] for s in ss],
                     color=st["color"], marker=st["marker"],
                     label=st["label"], capsize=3, lw=1.5)
    axC.set_xlabel(r"$1/L$ ($\mathrm{\AA}^{-1}$)")
    axC.set_ylabel(r"viscosity $\eta$ (mPa$\cdot$s)")
    axC.set_title("C  viscosity vs system size")
    axC.legend(fontsize=8)

    # -- Panel D: D_PBC and D_0 vs 1/L ------------------------------------
    axD = ax[1, 1]
    for b, st in BRANCHES.items():
        if b not in summ:
            continue
        ss = [s for s in SIZES if s in summ[b]]
        xs = [inv_L(s) for s in ss]
        axD.errorbar(xs, [summ[b][s]["D_PBC_ang2_ps"] for s in ss],
                     yerr=[summ[b][s]["D_PBC_sem"] for s in ss],
                     color=st["color"], marker=st["marker"], mfc="none",
                     ls="--", capsize=3, lw=1.2, label=f"{st['label']}  $D_{{PBC}}$")
        axD.errorbar(xs, [summ[b][s]["D_0_ang2_ps"] for s in ss],
                     yerr=[summ[b][s]["D_0_sem"] for s in ss],
                     color=st["color"], marker=st["marker"],
                     ls="-", capsize=3, lw=1.5, label=f"{st['label']}  $D_0$")
    axD.set_xlabel(r"$1/L$ ($\mathrm{\AA}^{-1}$)")
    axD.set_ylabel(r"diffusion coefficient ($\mathrm{\AA}^2$/ps)")
    axD.set_title("D  diffusion: uncorrected $D_{PBC}$ and Yeh-Hummer $D_0$")
    axD.legend(fontsize=7)

    for a in ax.flat:
        a.grid(True, ls="--", alpha=0.4)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(args.plot_dir, "figS4_replication.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out}")

    # =====================================================================
    # Figure 2 -- performance: master vs native-spline
    #   wall-time per MD step | wall-time per production segment | speedup
    # =====================================================================
    fig, (axp, axw, axr) = plt.subplots(1, 3, figsize=(15, 4.5))
    for b, st in BRANCHES.items():
        if b not in summ:
            continue
        ss = [s for s in SIZES if s in summ[b]]
        atoms = [s * 3 for s in ss]
        axp.errorbar(atoms, [summ[b][s]["time_per_step_s"] for s in ss],
                     yerr=[summ[b][s]["time_per_step_sem"] for s in ss],
                     color=st["color"], marker=st["marker"], label=st["label"],
                     capsize=3, lw=1.5)
        if all("seg_walltime_s" in summ[b][s] for s in ss):
            axw.errorbar(atoms, [summ[b][s]["seg_walltime_s"] / 3600.0 for s in ss],
                         yerr=[summ[b][s]["seg_walltime_sem"] / 3600.0 for s in ss],
                         color=st["color"], marker=st["marker"], label=st["label"],
                         capsize=3, lw=1.5)
    axp.set_ylabel("wall-time per MD step (s)")
    axp.set_title("per-step cost vs system size")
    axw.set_ylabel("wall-time per production segment (h)")
    axw.set_title("segment wall-time vs system size")
    for a in (axp, axw):
        a.set_xlabel("number of atoms")
        a.set_xscale("log", base=2)
        a.set_yscale("log")
        a.set_xticks([s * 3 for s in SIZES])
        a.set_xticklabels([str(s * 3) for s in SIZES], rotation=45, fontsize=7)
        a.grid(True, which="both", ls="--", alpha=0.4)
        a.legend(fontsize=8)

    # speedup = master time-per-step / native-spline time-per-step, per size
    m = summ.get("master", {})
    o = summ.get("feature-nnp-native-spline", {})
    ss = [s for s in SIZES if s in m and s in o]
    if ss:
        axr.plot([s * 3 for s in ss],
                 [m[s]["time_per_step_s"] / o[s]["time_per_step_s"] for s in ss],
                 "k-o", lw=1.5)
        axr.axhline(1.0, color="grey", ls="--", lw=1)
    axr.set_xlabel("number of atoms")
    axr.set_ylabel("speedup (master / native-spline, per step)")
    axr.set_title("native-spline speedup over master")
    axr.set_xscale("log", base=2)
    axr.set_xticks([s * 3 for s in SIZES])
    axr.set_xticklabels([str(s * 3) for s in SIZES], rotation=45, fontsize=7)
    axr.grid(True, which="both", ls="--", alpha=0.4)

    fig.suptitle("NNP MD performance -- master vs feature/nnp-native-spline "
                 "(one icelake node)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = os.path.join(args.plot_dir, "figS4_performance.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out}")

    # =====================================================================
    # Figure 3 -- accuracy: do the two branches give comparable physics?
    # eta and D_0 side by side; bars should agree within error bars since the
    # two branches run from identical NVE starting configurations.
    # =====================================================================
    fig, (axe, axd) = plt.subplots(1, 2, figsize=(12, 4.5))
    width = 0.38
    ss_common = [s for s in SIZES if all(s in summ.get(b, {}) for b in BRANCHES)]
    xpos = np.arange(len(ss_common))
    for i, (b, st) in enumerate(BRANCHES.items()):
        off = (i - 0.5) * width
        axe.bar(xpos + off, [summ[b][s]["eta_mPa_s"] for s in ss_common], width,
                yerr=[summ[b][s]["eta_sem"] for s in ss_common],
                color=st["color"], label=st["label"], capsize=3)
        axd.bar(xpos + off, [summ[b][s]["D_0_ang2_ps"] for s in ss_common], width,
                yerr=[summ[b][s]["D_0_sem"] for s in ss_common],
                color=st["color"], label=st["label"], capsize=3)
    for a, ylab, title in ((axe, r"$\eta$ (mPa$\cdot$s)", "viscosity"),
                           (axd, r"$D_0$ ($\mathrm{\AA}^2$/ps)", "diffusion D_0")):
        a.set_xticks(xpos)
        a.set_xticklabels([f"N={s}" for s in ss_common])
        a.set_ylabel(ylab)
        a.set_title(f"accuracy cross-check: {title}")
        a.legend(fontsize=8)
        a.grid(True, axis="y", ls="--", alpha=0.4)
    fig.suptitle("Branches should agree within error bars "
                 "(identical NVE starting configs)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.join(args.plot_dir, "figS4_accuracy.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out}")


if __name__ == "__main__":
    main()
