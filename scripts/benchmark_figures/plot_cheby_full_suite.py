#!/usr/bin/env python3
"""Complete chebyshev benchmark figure suite (size / multi-node / core scaling).
CSD3 icelake, H2O-N NNP MD, 76 cores/node.  Size: cheby/master 76x1 (jobs 20-06 /
10-06), native-spline (15-05); core: cheby+master same job 30565007, native-spline
30901621; multi-node: cheby+master same job 30565008.  Styled per thesis_figures.py."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../plots/cheby_benchmark_figs")
os.makedirs(OUT, exist_ok=True)
DATA_OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/cheby_benchmark")
os.makedirs(DATA_OUT, exist_ok=True)

# ---------------- thesis aesthetic (mirrors thesis_figures.py) ----------------
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
C = {"master": CAM["slate_3"], "native-spline": CAM["blue_warm"], "chebyshev": CAM["purple"]}
M = {"master": "o", "native-spline": "^", "chebyshev": "D"}

LBL = {"master": "master", "native-spline": "feature/nnp-native-spline",
       "chebyshev": "feature/nnp-chebyshev"}


def branch_handles(branches, ls="-"):
    return [plt.Line2D([0], [0], color=C[b], marker=M[b], ls=ls, lw=1.7,
                       ms=6.5, label=LBL[b]) for b in branches]


def top_legend(fig, handles, ncol, top=0.90):
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.00),
               ncol=ncol, frameon=False, fontsize=12, columnspacing=2.0,
               handlelength=2.6, handletextpad=0.6)
    fig.tight_layout(rect=[0, 0, 1, top])


def panel_letter(ax, letter):
    # above the axes so it clears the centred title
    ax.text(-0.02, 1.07, f"({letter})", transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="bottom", ha="left", color=CAM["slate_4"])


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/{name}.{ext}")
    plt.close(fig)


# ---------- SIZE SCALING (76x1) ----------
Nsz = np.array([64, 256, 512, 1024, 2048, 4096])
sz = {
    "master":        np.array([0.036160, 0.045658, 0.078658, 0.196540, 0.807022, 3.479442]),
    "native-spline": np.array([0.039248, 0.039694, 0.046292, 0.072146, 0.111614, 0.204228]),
    "chebyshev":     np.array([0.035932, 0.050208, 0.062922, 0.092486, 0.155230, 0.293732]),
}
# ---------- CORE SCALING (pure MPI, N=1024) ----------
cores = np.array([1, 2, 4, 8, 16, 32, 64, 76])
cs = {
    "master":        np.array([8.000887, 4.263673, 2.125633, 1.084407, 0.575807, 0.296140, 0.182683, 0.195253]),
    "native-spline": np.array([1.297553, 0.739707, 0.408620, 0.209600, 0.108280, 0.062020, 0.040693, 0.068900]),
    "chebyshev":     np.array([3.007183, 1.686183, 0.945683, 0.477153, 0.241927, 0.126557, 0.074240, 0.092067]),
}
# ---------- MULTI-NODE (nodes 1/2/4 = 76/152/304 ranks) ----------
nodes = np.array([1, 2, 4])
mn = {
    1024: {"master":    np.array([0.180377, 0.106880, 0.071587]),
           "chebyshev": np.array([0.092053, 0.061330, 0.046870])},
    4096: {"master":    np.array([3.258930, 1.715460, 0.885183]),
           "chebyshev": np.array([0.288730, 0.164613, 0.106697])},
}

# ===== FIG 5: size scaling =====
fig, (a1, a2) = plt.subplots(1, 2, figsize=(W_TEXT * 1.85, 4.7))
for b in ("master", "native-spline", "chebyshev"):
    a1.plot(Nsz, sz[b], M[b] + "-", color=C[b])
a1.set_xscale("log", base=2); a1.set_yscale("log")
a1.set_xticks(Nsz); a1.set_xticklabels(Nsz)
a1.set_xlabel(r"System size ($N$ water molecules)"); a1.set_ylabel("Time per MD step (s)")
a1.set_title(r"Size scaling, 76 MPI $\times$ 1 OMP, 1 node")
a1.grid(True, which="both"); panel_letter(a1, "a")
for b in ("native-spline", "chebyshev"):
    a2.plot(Nsz, sz["master"] / sz[b], M[b] + "-", color=C[b])
baseline = a2.axhline(1, color=CAM["slate_3"], ls="--", lw=1.0)
a2.legend([baseline], [r"Master baseline ($1\times$)"],
          loc="upper left", fontsize=10, frameon=True, framealpha=0.93)
a2.set_xscale("log", base=2); a2.set_xticks(Nsz); a2.set_xticklabels(Nsz)
a2.set_xlabel(r"System size ($N$ water molecules)"); a2.set_ylabel("Speedup vs master")
a2.set_title("Speedup over master grows with system size")
a2.grid(True); panel_letter(a2, "b")
top_legend(fig, branch_handles(("master", "native-spline", "chebyshev")), ncol=3)
save(fig, "fig5_size_scaling")

# ===== FIG 6: multi-node =====
fig, (a1, a2) = plt.subplots(1, 2, figsize=(W_TEXT * 1.85, 4.7))
for N in (1024, 4096):
    for b in ("master", "chebyshev"):
        ls = "-" if N == 1024 else "--"
        a1.plot(nodes, mn[N][b], M[b] + ls, color=C[b])
a1.set_yscale("log"); a1.set_xticks(nodes); a1.set_xticklabels([f"{n}\n({n*76})" for n in nodes])
a1.set_xlabel("Nodes (MPI ranks)"); a1.set_ylabel("Time per MD step (s)")
a1.set_title("Multi-node strong scaling")
a1.grid(True, which="both"); panel_letter(a1, "a")
for N in (1024, 4096):
    for b in ("master", "chebyshev"):
        eff = 100 * mn[N][b][0] / (mn[N][b] * nodes)
        ls = "-" if N == 1024 else "--"
        a2.plot(nodes, eff, M[b] + ls, color=C[b])
ideal_eff = a2.axhline(100, color=CAM["slate_3"], ls=":", lw=1.0)
a2.legend([ideal_eff], ["Ideal (100%)"],
          loc="lower left", fontsize=10, frameon=True, framealpha=0.93)
a2.set_xticks(nodes); a2.set_xticklabels(nodes)
a2.set_xlabel("Nodes"); a2.set_ylabel("Parallel efficiency (%)")
a2.set_title("Multi-node parallel efficiency")
a2.grid(True); panel_letter(a2, "b")
# ncol=2 column-major -> row per branch, column per system size.
mn_handles = [plt.Line2D([0], [0], color=C[b], marker=M[b],
                         ls=("-" if N == 1024 else "--"), lw=1.7, ms=6.5,
                         label=f"{LBL[b]}, $N$={N}")
              for N in (1024, 4096) for b in ("master", "chebyshev")]
top_legend(fig, mn_handles, ncol=2, top=0.87)
save(fig, "fig6_multinode")

# ===== FIG 7: pure-MPI strong scaling, N=1024 =====
fig, (a1, a2) = plt.subplots(1, 2, figsize=(W_TEXT * 1.85, 4.7))
for b in ("master", "native-spline", "chebyshev"):
    a1.plot(cores, cs[b][0] / cs[b], M[b] + "-", color=C[b])
h_ideal, = a1.plot(cores, cores, "--", color=CAM["slate_3"], lw=1.0)
a1.legend([h_ideal], ["Ideal (linear)"],
          loc="upper left", fontsize=10, frameon=True, framealpha=0.93)
a1.set_xlabel("MPI ranks (1 OMP each)"); a1.set_ylabel("Speedup vs 1 core")
a1.set_title(r"Pure-MPI strong scaling, $N=1024$")
a1.grid(True); panel_letter(a1, "a")
for b in ("master", "native-spline", "chebyshev"):
    eff = 100 * cs[b][0] / (cs[b] * cores)
    a2.plot(cores, eff, M[b] + "-", color=C[b])
ideal_100 = a2.axhline(100, color=CAM["slate_3"], ls=":", lw=1.0)
a2.legend([ideal_100], ["Ideal (100%)"],
          loc="lower left", fontsize=10, frameon=True, framealpha=0.93)
a2.set_xlabel("MPI ranks"); a2.set_ylabel("Parallel efficiency (%)")
a2.set_title(r"Pure-MPI parallel efficiency, $N=1024$")
a2.grid(True); panel_letter(a2, "b")
top_legend(fig, branch_handles(("master", "native-spline", "chebyshev")), ncol=3)
save(fig, "fig7_core_scaling")

# ===== combined CSV =====
with open(f"{DATA_OUT}/combined_full_suite.csv", "w") as f:
    f.write("dataset,branch,N_molecules,nodes,mpi_ranks,omp_threads,tps_s\n")
    for b in sz:
        for n, t in zip(Nsz, sz[b]): f.write(f"size_scaling,{b},{n},1,76,1,{t}\n")
    for b in cs:
        for c, t in zip(cores, cs[b]): f.write(f"core_scaling,{b},1024,1,{c},1,{t}\n")
    for N in mn:
        for b in mn[N]:
            for nd, t in zip(nodes, mn[N][b]): f.write(f"multinode,{b},{N},{nd},{nd*76},1,{t}\n")

print("Wrote suite (png+pdf) to", OUT)
print("Size scaling cheby-vs-master:", {int(n): round(sz['master'][i]/sz['chebyshev'][i], 2) for i, n in enumerate(Nsz)})
