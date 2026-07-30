"""One-off supervisor figure: self-diffusion vs concentration, Madrid-2019.

Same data (kappa_vs_c.npz), palette and axis conventions as make_figures.py
(fig3): x = molarity via MOL_TO_M, errors = SEM over 5 seeds per (m, L).
D stored in A^2/ps -> x10 = 1e-9 m^2/s.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NPZ = sys.argv[1]
OUT = sys.argv[2]

C_NA, C_CL, C_O = "#5366E0", "#4DB78C", "#CD3572"  # Cambridge indigo/green/cherry
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
MOL_TO_M = {0.25: 0.231, 0.5: 0.491, 1.0: 0.960, 2.0: 1.915, 4.0: 3.666}
# infinite-dilution experiment: Robinson & Stokes (ions), Holz 2000 (water)
EXPT_D0 = {"Na": 1.334, "Cl": 2.032, "O": 2.299}
LABEL = {"Na": r"Na$^+$", "Cl": r"Cl$^-$", "O": r"H$_2$O"}

d = np.load(NPZ, allow_pickle=True)
mol, L = d["runs_mol"], d["runs_L"]
mols = sorted(set(mol))
c = np.array([MOL_TO_M[m] for m in mols])

fig, ax = plt.subplots(figsize=(5.2, 4.0))
for sp, col in (("O", C_O), ("Cl", C_CL), ("Na", C_NA)):
    y = d[f"runs_D{sp}"] * 10.0          # -> 1e-9 m^2/s
    for Lval, fmt, mfc, lw in ((31.05, "o--", "white", 1.2),
                               (43.47, "o-", None, 1.6)):
        sel_L = L == Lval
        mean = np.array([y[(mol == m) & sel_L].mean() for m in mols])
        sem = np.array([y[(mol == m) & sel_L].std(ddof=1) / np.sqrt(5)
                        for m in mols])
        ax.errorbar(c, mean, yerr=sem, fmt=fmt, ms=5.5, color=col,
                    mfc=mfc if mfc else col, mew=1.4, lw=lw,
                    elinewidth=1.0, capsize=2)
    ax.plot(0, EXPT_D0[sp], marker="_", ms=12, mew=2.2, ls="", color=col,
            zorder=5, clip_on=False)

hs = ([plt.Line2D([], [], color=col, lw=2.2, label=LABEL[sp])
       for sp, col in (("O", C_O), ("Na", C_NA), ("Cl", C_CL))]
      + [plt.Line2D([], [], marker="o", ls="--", color=MUTED, mfc="white",
                    mew=1.4, label=r"$L=31.1$ Å"),
         plt.Line2D([], [], marker="o", ls="-", color=MUTED,
                    label=r"$L=43.5$ Å"),
         plt.Line2D([], [], marker="_", ls="", ms=12, mew=2.2, color=MUTED,
                    label=r"Experiment, $c\to0$")])
fig.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, 1.06),
           ncol=2, frameon=False, fontsize=9, columnspacing=2.0,
           handlelength=2.4, handletextpad=0.6, labelcolor=INK2)
ax.set_xlim(-0.08, 3.9)
ax.set_ylim(0, 2.6)
ax.grid(True, color=GRID, linewidth=0.7, alpha=0.9)
ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_color(MUTED)
ax.tick_params(colors=INK2, labelsize=9)
ax.set_xlabel("c (mol/L)", color=INK, fontsize=10)
ax.set_ylabel(r"$D_{\mathrm{PBC}}$ ($10^{-9}$ m$^2$ s$^{-1}$)", color=INK,
              fontsize=10)

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, f"fig9_D_vs_c.{ext}"),
                dpi=400 if ext == "png" else None,
                facecolor="white", bbox_inches="tight")
print("wrote fig9_D_vs_c.png/.pdf")
