#!/usr/bin/env python3
"""Head-to-head: Dhruv's pristine PR head (8be1dfa50f) vs the SAME head with
the bit-exact atom-level OpenMP graft (omp-atom-bitexact-v2). Both binaries
built from the SAME configured CMake tree (-O3 -xCORE-AVX2 -fp-model=precise
-qopenmp, Intel 2022.1.0), run on the SAME icelake node in the SAME job
(31837144) -> the ONLY variable is the graft.

Two axes, H2O-1024 NNP MD (3072 atoms), 100 MD steps, 5 reps + warm-up,
t/step from qs_mol_dyn_low:
  (1) OMP thread scaling   1 MPI rank, threads {1..76}
  (2) pure-MPI strong scale OMP=1, ranks {1..76}
Accuracy leg of the same job: forces bit-identical (max abs diff 0.0) between
the two binaries at OMP=1 and OMP=76; energy diff 0.0 Ha.
Styled to match scripts/benchmark_figures/thesis_figures.py.
"""
import glob
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/cp2k_omp_headtohead/NNP")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../plots/omp_headtohead")
os.makedirs(OUT, exist_ok=True)

CAM = {
    "blue_warm": "#00BDB6", "blue_dark": "#133844",
    "crest": "#FD8153", "crest_dark": "#DD3025",
    "indigo": "#5366E0", "green_dark": "#13553A", "purple": "#A368DF",
    "slate_2": "#B5BDC8", "slate_3": "#546072", "slate_4": "#232830",
}
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Liberation Serif", "Bitstream Vera Serif", "serif"],
    "font.size": 12, "axes.labelsize": 13, "axes.titlesize": 13,
    "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 11,
    "figure.titlesize": 14,
    "text.color": CAM["slate_4"], "axes.edgecolor": CAM["slate_4"],
    "axes.labelcolor": CAM["slate_4"], "xtick.color": CAM["slate_4"],
    "ytick.color": CAM["slate_4"], "axes.linewidth": 0.8,
    "axes.grid": True, "grid.color": CAM["slate_2"], "grid.linestyle": "--",
    "grid.linewidth": 0.5, "grid.alpha": 0.55,
    "lines.linewidth": 1.7, "lines.markersize": 6.5, "lines.markeredgewidth": 0.8,
    "mathtext.fontset": "cm", "figure.dpi": 100, "savefig.dpi": 300,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.05,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
W_TEXT = 6.3


def panel_letter(ax, letter):
    # top-left, above the axes, so it never collides with the centred title
    ax.text(-0.02, 1.07, f"({letter})", transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="bottom", ha="left", color=CAM["slate_4"])

BASE, OMP = "dhruv-head-8be1dfa", "dhruv-omp-bitexact-v2"


def load(kind, label):
    """kind in {omp_thread, core}; returns array cols of the summary CSV."""
    pat = os.path.join(R, f"NNP_{kind}_scaling_{label}_*",
                       f"results_{kind}_scaling_{label}_*.csv")
    fs = [f for f in sorted(glob.glob(pat)) if not f.endswith("_raw.csv")]
    if not fs:
        raise SystemExit(f"no CSV for {kind}/{label}: {pat}")
    return np.loadtxt(fs[-1], delimiter=",", comments="#")


# col layout (both drivers): 0 sweep-var, 2 total_cores, 4 t/step mean,
# 5 t/step std, 10 speedup, 11 efficiency
omp_b, omp_g = load("omp_thread", BASE), load("omp_thread", OMP)
mpi_b, mpi_g = load("core", BASE), load("core", OMP)

# Thesis colour semantics (thesis_figures.py BRANCH_STYLE): the unmodified
# reference is neutral slate with marker "o" (like master); the optimised
# variant is Cambridge warm blue with marker "^" (like native-spline).
# Linestyle separates the parallel mode: OMP solid, pure MPI dotted.
SER = {
    "base_omp": (omp_b, CAM["slate_3"], "o", "-", "PR head, OMP (1 rank)"),
    "graft_omp": (omp_g, CAM["blue_warm"], "^", "-", "PR head + graft, OMP (1 rank)"),
    "base_mpi": (mpi_b, CAM["slate_3"], "o", ":", "PR head, pure MPI"),
    "graft_mpi": (mpi_g, CAM["blue_warm"], "^", ":", "PR head + graft, pure MPI"),
}

fig, ax = plt.subplots(1, 2, figsize=(W_TEXT * 1.85, 4.7))
for a, c, m, ls, lab in SER.values():
    ax[0].errorbar(a[:, 2], a[:, 4], yerr=a[:, 5], marker=m, color=c,
                   ls=ls, capsize=2.5)
ax[0].set_xscale("log", base=2); ax[0].set_yscale("log")
ax[0].set_xticks(omp_b[:, 2])
ax[0].set_xticklabels([f"{int(x)}" for x in omp_b[:, 2]], rotation=45)
ax[0].set_xlabel("Cores")
ax[0].set_ylabel("Time per MD step (s)")
ax[0].set_title(r"$N = 1024$ H$_2$O, time per step")
panel_letter(ax[0], "a")

for a, c, m, ls, lab in SER.values():
    ax[1].plot(a[:, 2], a[:, 10], marker=m, color=c, ls=ls)
ideal = omp_b[:, 2]
h_ideal, = ax[1].plot(ideal, ideal, color=CAM["slate_3"], lw=1.0, ls="--")
ax[1].set_xscale("log", base=2); ax[1].set_yscale("log", base=2)
ax[1].set_xticks(ideal)
ax[1].set_xticklabels([f"{int(x)}" for x in ideal], rotation=45)
ax[1].set_yticks([1, 2, 4, 8, 16, 32, 64]); ax[1].set_yticklabels(
    ["1", "2", "4", "8", "16", "32", "64"])
ax[1].set_xlabel("Cores")
ax[1].set_ylabel("Speedup vs 1 core")
ax[1].set_title("Strong-scaling speedup")
ax[1].legend([h_ideal], ["Ideal (linear)"],
             loc="upper left", fontsize=10, frameon=True, framealpha=0.93)
panel_letter(ax[1], "b")

# One frameless figure-level legend, thesis-style: ncol=2 column-major gives
# row 1 = PR head (OMP | MPI), row 2 = PR head + graft (OMP | MPI).
handles = [plt.Line2D([0], [0], color=c, marker=m, ls=ls, lw=1.7, ms=6.5,
                      label=lab) for _, c, m, ls, lab in SER.values()]
fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.00),
           ncol=2, frameon=False, fontsize=12, columnspacing=2.0,
           handlelength=2.6, handletextpad=0.6)


fig.tight_layout(rect=[0, 0, 1, 0.86])
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/omp_headtohead_scaling.{ext}")
plt.close(fig)

# ---- console summary table (markdown, ready for the PR / report) -----------
def row(a, i):
    return a[i, 4], a[i, 5], a[i, 10], a[i, 11]

print("### Head-to-head, N = 1024 H2O (3072 atoms), 100 MD steps, "
      "icelake, job 31837144\n")
print("| cores | head OMP t/step [s] | graft OMP t/step [s] | OMP ratio "
      "| head MPI t/step [s] | graft MPI t/step [s] |")
print("|---:|---:|---:|---:|---:|---:|")
for i in range(omp_b.shape[0]):
    n = int(omp_b[i, 2])
    tb, sb, _, _ = row(omp_b, i)
    tg, sg, _, _ = row(omp_g, i)
    mb = row(mpi_b, i)[0]
    mg = row(mpi_g, i)[0]
    print(f"| {n} | {tb:.4f} ± {sb:.4f} | {tg:.4f} ± {sg:.4f} "
          f"| {tb/tg:.2f}× | {mb:.4f} | {mg:.4f} |")

print(f"graft OMP-76 vs best full-node baseline (head, 76 MPI): "
      f"{mpi_b[-1, 4] / omp_g[-1, 4]:.2f}x")
print(f"serial overhead of graft (OMP=1): "
      f"{(omp_g[0, 4] / omp_b[0, 4] - 1) * 100:+.1f}%")
print(f"wrote {OUT}/omp_headtohead_scaling.[png|pdf]")
