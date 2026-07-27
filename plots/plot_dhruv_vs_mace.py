#!/usr/bin/env python3
"""Head-to-head: Dhruv's linked-cell NNP (pr/nnp-cpu-cell-list) vs my
feature/nnp-mace (= chebyshev + idle MACE backend). Both binaries built at the
SAME AVX-2 profile (-O3 -xCORE-AVX2 -fp-model=precise) + sequential MKL on the
same toolchain, so the only difference is the NNP source.

Three axes, H2O-N NNP MD, CSD3 icelake (76 cores/node), 100 MD steps,
5 timed reps + 1 warm-up, t/step from qs_mol_dyn_low:
  (1) SIZE scaling          76 MPI x 1 OMP, N = 64..4096
  (2) PURE-MPI strong scale N = 1024, cores {1,2,4,8,16,19,32,38,76}
  (3) OMP thread scaling    1 MPI rank, threads {1..76}, N = 1024
Jobs 31038227 (axes 1-2) + 31038352 (axis 3).
Styled to match plots/thesis_figures.py.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../results")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../plots/dhruv_vs_mace")
os.makedirs(OUT, exist_ok=True)

CAM = {
    "blue_warm": "#00BDB6", "blue_dark": "#133844",
    "crest": "#FD8153", "crest_dark": "#DD3025",
    "indigo": "#5366E0", "green_dark": "#13553A", "purple": "#A368DF",
    "slate_2": "#B5BDC8", "slate_3": "#546072", "slate_4": "#232830",
}
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Liberation Serif", "serif"],
    "font.size": 12, "axes.labelsize": 13, "axes.titlesize": 13,
    "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 10,
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

C = {"dhruv": CAM["slate_3"], "mace": CAM["crest_dark"]}
M = {"dhruv": "o", "mace": "D"}
LAB = {"dhruv": "Dhruv  (linked-cell, spline radials)",
       "mace": "nnp-mace  (table-free minimax)"}


def load(path):
    return np.genfromtxt(path, delimiter=",", comments="#")


SZ_D = R + "/cp2k_dhruv_cell_list/NNP/NNP_size_scaling_dhruv-cell-list_25-06_14-50/results_size_scaling_dhruv-cell-list_25-06_14-50.csv"
SZ_M = R + "/cp2k_feature_mace/NNP/NNP_size_scaling_feature-nnp-mace_25-06_15-07/results_size_scaling_feature-nnp-mace_25-06_15-07.csv"
CO_D = R + "/cp2k_dhruv_cell_list/NNP/NNP_core_scaling_dhruv-cell-list_25-06_15-21-34_31038227/results_core_scaling_dhruv-cell-list_25-06_15-21-34_31038227.csv"
CO_M = R + "/cp2k_feature_mace/NNP/NNP_core_scaling_feature-nnp-mace_25-06_16-16-09_31038227/results_core_scaling_feature-nnp-mace_25-06_16-16-09_31038227.csv"
OM_D = R + "/cp2k_dhruv_cell_list/NNP/NNP_omp_thread_scaling_dhruv-cell-list_25-06_14-50/results_omp_thread_scaling_dhruv-cell-list_25-06_14-50.csv"
OM_M = R + "/cp2k_feature_mace/NNP/NNP_omp_thread_scaling_feature-nnp-mace_25-06_18-39/results_omp_thread_scaling_feature-nnp-mace_25-06_18-39.csv"


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/{name}.{ext}")
    plt.close(fig)
    print("wrote", name)


# ============================ FIG 1: SIZE SCALING ============================
d, m = load(SZ_D), load(SZ_M)
# cols: n_molecules,n_reps,tps_mean,tps_std,tps_min,wall_mean,wall_std,wall_min
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
for k, a in [("dhruv", d), ("mace", m)]:
    ax[0].errorbar(a[:, 0], a[:, 2], yerr=a[:, 3], marker=M[k], color=C[k],
                   label=LAB[k], capsize=3)
ax[0].set_xscale("log", base=2); ax[0].set_yscale("log")
ax[0].set_xlabel("system size  (H$_2$O molecules)")
ax[0].set_ylabel("time per MD step  (s)")
ax[0].set_title("Size scaling  (76 MPI $\\times$ 1 OMP)")
ax[0].set_xticks(d[:, 0]); ax[0].set_xticklabels([int(x) for x in d[:, 0]])
ax[0].legend(frameon=False)
# relative: mace/dhruv ratio of t/step (>1 = mace slower)
ratio = m[:, 2] / d[:, 2]
ax[1].axhline(1.0, color=CAM["slate_3"], lw=1.0, ls="--")
ax[1].plot(d[:, 0], ratio, marker="D", color=CAM["crest_dark"])
for x, r in zip(d[:, 0], ratio):
    ax[1].annotate(f"{r:.2f}", (x, r), textcoords="offset points",
                   xytext=(0, 7), ha="center", fontsize=9)
ax[1].set_xscale("log", base=2)
ax[1].set_xlabel("system size  (H$_2$O molecules)")
ax[1].set_ylabel("nnp-mace t/step  $\\div$  Dhruv t/step")
ax[1].set_title("Relative cost  ($<1$: nnp-mace faster)")
ax[1].set_xticks(d[:, 0]); ax[1].set_xticklabels([int(x) for x in d[:, 0]])
ax[1].set_ylim(0.7, 1.15)
save(fig, "fig1_size_scaling")

# ============================ FIG 2: PURE-MPI ===============================
d, m = load(CO_D), load(CO_M)
# cols: mpi,omp,cores,reps,tps_mean,tps_std,tps_min,wall_mean,wall_std,wall_min,speedup,eff
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
for k, a in [("dhruv", d), ("mace", m)]:
    ax[0].errorbar(a[:, 2], a[:, 4], yerr=a[:, 5], marker=M[k], color=C[k],
                   label=LAB[k], capsize=3)
ax[0].set_xscale("log", base=2); ax[0].set_yscale("log")
ax[0].set_xlabel("MPI ranks  (= cores, 1 OMP)")
ax[0].set_ylabel("time per MD step  (s)")
ax[0].set_title("Pure-MPI strong scaling  (N = 1024)")
ax[0].set_xticks(d[:, 2]); ax[0].set_xticklabels([int(x) for x in d[:, 2]])
ax[0].legend(frameon=False)
for k, a in [("dhruv", d), ("mace", m)]:
    ax[1].plot(a[:, 2], a[:, 11], marker=M[k], color=C[k], label=LAB[k])
ax[1].axhline(100, color=CAM["slate_3"], lw=0.8, ls=":")
ax[1].set_xscale("log", base=2)
ax[1].set_xlabel("MPI ranks")
ax[1].set_ylabel("parallel efficiency  (%)")
ax[1].set_title("Strong-scaling efficiency")
ax[1].set_xticks(d[:, 2]); ax[1].set_xticklabels([int(x) for x in d[:, 2]])
ax[1].legend(frameon=False)
save(fig, "fig2_pure_mpi")

# ============================ FIG 3: OMP THREADS ============================
d, m = load(OM_D), load(OM_M)
# cols: omp,mpi,cores,reps,tps_mean,tps_std,tps_min,wall_mean,wall_std,wall_min,speedup,eff
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
for k, a in [("dhruv", d), ("mace", m)]:
    ax[0].plot(a[:, 0], a[:, 10], marker=M[k], color=C[k], label=LAB[k])
ax[0].plot(d[:, 0], d[:, 0], color=CAM["slate_2"], lw=1.0, ls="--", label="ideal")
ax[0].set_xscale("log", base=2); ax[0].set_yscale("log", base=2)
ax[0].set_xlabel("OMP threads  (1 MPI rank)")
ax[0].set_ylabel("speedup vs 1 thread")
ax[0].set_title("OpenMP thread scaling  (N = 1024)")
ax[0].set_xticks(d[:, 0]); ax[0].set_xticklabels([int(x) for x in d[:, 0]])
ax[0].legend(frameon=False)
for k, a in [("dhruv", d), ("mace", m)]:
    ax[1].plot(a[:, 0], a[:, 4], marker=M[k], color=C[k], label=LAB[k])
ax[1].set_xscale("log", base=2); ax[1].set_yscale("log")
ax[1].set_xlabel("OMP threads  (1 MPI rank)")
ax[1].set_ylabel("time per MD step  (s)")
ax[1].set_title("Absolute t/step under threading")
ax[1].set_xticks(d[:, 0]); ax[1].set_xticklabels([int(x) for x in d[:, 0]])
ax[1].legend(frameon=False)
save(fig, "fig3_omp_threads")

print("\nall figures ->", OUT)
