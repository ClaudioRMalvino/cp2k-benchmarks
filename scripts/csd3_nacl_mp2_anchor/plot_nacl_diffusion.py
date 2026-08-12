#!/usr/bin/env python3
"""NaCl(aq) MP2 C-NNP anchor: water-oxygen MSD and Yeh-Hummer finite-size
extrapolation, D(L) = D_0 - kB T xi / (6 pi eta L), xi = 2.837297 (cubic).

Reads results/nacl_mp2_anchor/{msd_cube2.csv, msd_cube3.csv,
diffusion_summary.csv} written by scripts/csd3_nacl_mp2_anchor/
analyze_diffusion.py. Two cells (cube2 24.84 A / 1500 atoms, cube3 37.26 A /
5064 atoms) x 5 NVE segments each -> D vs 1/L with SEM error bars, linear
extrapolation to 1/L = 0 gives the infinite-box D_0 and, from the slope, the
implied shear viscosity eta.

Styled to match scripts/benchmark_figures/thesis_figures.py.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/nacl_mp2_anchor")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../plots/nacl_mp2_anchor")
os.makedirs(OUT, exist_ok=True)

CAM = {
    "blue_warm": "#00BDB6", "blue_dark": "#133844",
    "crest": "#FD8153", "crest_dark": "#DD3025",
    "indigo": "#5366E0", "green": "#4DB78C", "purple": "#A368DF",
    "slate_2": "#B5BDC8", "slate_3": "#546072", "slate_4": "#232830",
}
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Liberation Serif", "serif"],
    "font.size": 15, "axes.labelsize": 18, "axes.titlesize": 17,
    "xtick.labelsize": 15, "ytick.labelsize": 15, "legend.fontsize": 13,
    "text.color": CAM["slate_4"], "axes.edgecolor": CAM["slate_4"],
    "axes.labelcolor": CAM["slate_4"], "xtick.color": CAM["slate_4"],
    "ytick.color": CAM["slate_4"], "axes.linewidth": 0.8,
    "axes.grid": True, "grid.color": CAM["slate_2"], "grid.linestyle": "--",
    "grid.linewidth": 0.5, "grid.alpha": 0.55,
    "lines.linewidth": 1.7, "lines.markersize": 6.5, "lines.markeredgewidth": 0.8,
    "mathtext.fontset": "cm", "figure.dpi": 100, "savefig.dpi": 300,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.05,
})

KB = 1.380649e-23
XI = 2.837297
CELLS = {"cube2": (24.84, CAM["crest"], "s"),
         "cube3": (37.26, CAM["blue_warm"], "o")}

summ = np.genfromtxt(os.path.join(R, "diffusion_summary.csv"),
                     delimiter=",", names=True, dtype=None, encoding=None)

fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))

# ---- (a) MSD curves --------------------------------------------------------
# The CSVs span half of each cell's trajectory (cube2 50 ps, cube3 80 ps);
# show both over the common lag range.
tabs = {cell: np.loadtxt(os.path.join(R, f"msd_{cell}.csv"), delimiter=",",
                         skiprows=1) for cell in CELLS}
lag_max = min(tab[-1, 0] for tab in tabs.values())
for cell, (L, c, m) in CELLS.items():
    tab = tabs[cell][tabs[cell][:, 0] <= lag_max]
    lag, curves = tab[:, 0], tab[:, 1:]
    mean, lo, hi = curves.mean(1), curves.min(1), curves.max(1)
    ax[0].fill_between(lag, lo, hi, color=c, alpha=0.18, lw=0)
    ax[0].plot(lag, mean, color=c, label=f"L = {L:.1f} $\\AA$")
sel = summ["cell"] == "cube2"
flo, fhi = summ["fit_lo_ps"][sel][0], summ["fit_hi_ps"][sel][0]
ax[0].axvspan(flo, fhi, color=CAM["slate_2"], alpha=0.25, lw=0)
ax[0].set_xlabel("lag time [ps]")
ax[0].set_ylabel("water-O MSD [$\\AA^2$]")
ax[0].set_title("(a)")
ax[0].legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.09),
             ncol=2, columnspacing=1.2, handlelength=1.6, handletextpad=0.4)

# ---- (b) Yeh-Hummer D vs 1/L -----------------------------------------------
pts = {}
for cell, (L, c, m) in CELLS.items():
    sel = summ["cell"] == cell
    ds = summ["D_1e9_m2_s"][sel]
    ts = summ["T_K"][sel]
    d_mean, d_sem = ds.mean(), ds.std(ddof=1) / np.sqrt(len(ds))
    pts[cell] = (L, d_mean, d_sem, ts.mean())
    ax[1].errorbar(1.0 / L, d_mean, yerr=d_sem, marker=m, color=c,
                   capsize=3.5, ms=8, ls="none", label=f"L = {L:.1f} $\\AA$")

(l2, d2, s2, t2), (l3, d3, s3, t3) = pts["cube2"], pts["cube3"]
x2, x3 = 1.0 / l2, 1.0 / l3
slope = (d2 - d3) / (x2 - x3)
d0 = d2 - slope * x2
d0_err = np.sqrt((s3 * x2 / (x2 - x3)) ** 2 + (s2 * x3 / (x2 - x3)) ** 2)
tbar = 0.5 * (t2 + t3)
# slope [1e-9 m^2/s per 1/A] = -kB T xi/(6 pi eta) * 1e19  ->  eta [Pa s]
eta = -KB * tbar * XI / (6 * np.pi * slope * 1e-19) if slope < 0 else np.nan

xs = np.linspace(0, x2 * 1.12, 50)
ax[1].plot(xs, d0 + slope * xs, color=CAM["slate_3"], lw=1.2, ls="--")
ax[1].errorbar([0], [d0], yerr=[d0_err], marker="*", color=CAM["green"],
               ms=15, capsize=3.5, ls="none",
               label="$L \\to \\infty$ extrapolation")
ax[1].set_xlim(-0.003, x2 * 1.12)
ax[1].set_xlabel("1/L [$\\AA^{-1}$]")
ax[1].set_ylabel("$D_{\\mathrm{H_2O}}$ [$10^{-9}$ m$^2$ s$^{-1}$]")
ax[1].set_title("(b)")
ax[1].legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.09),
             ncol=3, columnspacing=1.0, handlelength=1.4, handletextpad=0.35)

for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/nacl_diffusion_yh.{ext}")
plt.close(fig)

print(f"cube2: D = {d2:.4f} +/- {s2:.4f} e-9 m^2/s  (T = {t2:.1f} K)")
print(f"cube3: D = {d3:.4f} +/- {s3:.4f} e-9 m^2/s  (T = {t3:.1f} K)")
print(f"D_0 (1/L -> 0) = {d0:.4f} +/- {d0_err:.4f} e-9 m^2/s")
print(f"implied viscosity eta = {eta*1e3:.3f} mPa s  (expt. water ~0.85 at 300 K)")
print(f"wrote {OUT}/nacl_diffusion_yh.[png|pdf]")
