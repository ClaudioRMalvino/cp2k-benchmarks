#!/usr/bin/env python3
"""Chebyshev hybrid-MPI/OMP + memory benchmark figures.
Data: SLURM job 30901621 (NNP_hybrid_mem), CSD3 icelake, 76 cores/node, H2O-N NNP MD.
Styled to match plots/thesis_figures.py (Cambridge palette, serif, grid choices).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.expanduser("~/cp2k-benchmarks/plots/cheby_benchmark_figs")
os.makedirs(OUT, exist_ok=True)

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
# master = neutral slate, native-spline = Cambridge warm blue, chebyshev = Cambridge red
C = {"master": CAM["slate_3"], "native-spline": CAM["blue_warm"], "chebyshev": CAM["crest_dark"]}
M = {"master": "o", "native-spline": "^", "chebyshev": "D"}
# Full branch names in legends, as in thesis_figures.py BRANCH_STYLE.
LBL = {"master": "master", "native-spline": "feature/nnp-native-spline",
       "chebyshev": "feature/nnp-chebyshev"}


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/{name}.{ext}")
    plt.close(fig)


# ----- decomposition order: MPI-heavy -> OMP-heavy -----
decomp_labels = ["76x1", "38x2", "19x4", "4x19", "2x38", "1x76"]
x = np.arange(len(decomp_labels))

# ----- t/step (s), N=1024 -----
tps_1024 = {
    "master":        [0.191500, 0.284450, 0.473875, 2.203925, 4.470600, 8.153175],
    "native-spline": [0.065900, 0.055525, 0.096725, 0.417450, 0.747650, 1.316400],
    "chebyshev":     [0.098175, 0.066950, 0.065050, 0.065475, 0.061250, 0.057025],
}
# ----- aggregate memory (MiB), N=1024 -----
mem_1024 = {
    "master":        [23978, 10241,  4588,  970,  505,  254],
    "native-spline": [24358, 10450,  4636,  972,  493,  251],
    "chebyshev":     [24092, 10450,  4740, 1172,  718,  476],
}
# ----- t/step (s), N=4096 (master/native-spline only high-MPI) -----
tps_4096 = {
    "master":        [3.304400, 4.421033, 7.312967, None, None, None],
    "native-spline": [0.190200, 0.235467, 0.428233, None, None, None],
    "chebyshev":     [0.295375, 0.251750, 0.246825, 0.255075, 0.231575, 0.212475],
}
# ----- chebyshev OMP ladder (1 MPI x T threads), N=1024 -----
omp_threads = np.array([1, 2, 4, 8, 16, 76])
omp_tps     = np.array([3.097925, 1.543250, 0.777600, 0.394250, 0.204025, 0.057025])
# ----- native-spline pure-MPI core scaling, N=1024 -----
ns_cores    = np.array([1, 2, 4, 8, 16, 32, 64, 76])
ns_speedup  = np.array([1.000, 1.754, 3.175, 6.191, 11.983, 20.922, 31.886, 18.832])
ns_eff      = np.array([100.0, 87.7, 79.4, 77.4, 74.9, 65.4, 49.8, 24.8])

# ===== FIG 1: hybrid decomposition sweep, t/step, N=1024 (headline) =====
fig, ax = plt.subplots(figsize=(W_TEXT, 4.4))
for b in ("master", "native-spline", "chebyshev"):
    ax.plot(x, tps_1024[b], M[b] + "-", color=C[b], label=LBL[b])
ax.set_yscale("log")
ax.set_xticks(x); ax.set_xticklabels(decomp_labels)
ax.set_xlabel(r"Decomposition (MPI ranks $\times$ OMP threads, 76 cores fixed)")
ax.set_ylabel("Time per MD step (s)")
ax.set_title(r"Hybrid MPI/OMP scaling, 1024 H$_2$O (3072 atoms), 1 node")
ax.axvspan(-0.3, 0.3, color=CAM["blue_dark"], alpha=0.05)
ax.annotate("pure MPI", (0, ax.get_ylim()[0]), color=CAM["slate_4"], ha="center", va="bottom", fontsize=9)
ax.annotate("pure OMP", (5, ax.get_ylim()[0]), color=CAM["slate_4"], ha="center", va="bottom", fontsize=9)
ax.grid(True, which="both")
ax.legend(title="Branch", title_fontsize=9, frameon=True, framealpha=0.93)
fig.tight_layout(); save(fig, "fig1_hybrid_sweep_N1024")

# ===== FIG 2: speed-memory tradeoff, N=1024 =====
fig, ax = plt.subplots(figsize=(W_TEXT, 4.6))
for b in ("master", "native-spline", "chebyshev"):
    ax.plot(tps_1024[b], mem_1024[b], M[b] + "-", color=C[b], label=LBL[b], alpha=0.95)
ax.annotate("chebyshev 1$\\times$76\n0.057 s, 476 MiB", (0.057025, 476), textcoords="offset points",
            xytext=(42, 12), fontsize=9, color=C["chebyshev"],
            arrowprops=dict(arrowstyle="->", color=C["chebyshev"], lw=1.0))
ax.annotate("native-spline 38$\\times$2\n0.056 s, 10450 MiB", (0.055525, 10450), textcoords="offset points",
            xytext=(42, -6), fontsize=9, color=C["native-spline"],
            arrowprops=dict(arrowstyle="->", color=C["native-spline"], lw=1.0))
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"Time per MD step (s) $\longleftarrow$ faster")
ax.set_ylabel(r"Aggregate node memory (MiB) $\longleftarrow$ lower")
ax.set_title("Speed vs node memory, 1024 H$_2$O, 1 node\n(each point = one MPI$\\times$OMP decomposition)")
ax.grid(True, which="both")
ax.legend(title="Branch", title_fontsize=9, frameon=True, framealpha=0.93)
fig.tight_layout(); save(fig, "fig2_speed_vs_memory_N1024")

# ===== FIG 3: chebyshev OMP scaling ladder =====
fig, ax = plt.subplots(figsize=(W_TEXT, 4.4))
speedup = omp_tps[0] / omp_tps
ax.plot(omp_threads, speedup, M["chebyshev"] + "-", color=C["chebyshev"],
        label="feature/nnp-chebyshev (measured)")
ax.plot(omp_threads, omp_threads, "--", color=CAM["slate_3"], lw=1.0, label="Ideal (linear)")
for t, s in zip(omp_threads, speedup):
    ax.annotate(f"{s:.1f}$\\times$", (t, s), textcoords="offset points", xytext=(7, -3),
                fontsize=9, color=CAM["slate_4"])
ax.set_xlabel("OpenMP threads (1 MPI rank)")
ax.set_ylabel(r"Speedup vs 1 thread")
ax.set_title("Chebyshev centre-level OpenMP scaling, 1024 H$_2$O")
ax.grid(True); ax.legend(frameon=True, framealpha=0.93)
fig.tight_layout(); save(fig, "fig3_cheby_omp_scaling")

# ===== FIG 4: speedup vs master bar chart (pure MPI 76x1, both sizes) =====
fig, ax = plt.subplots(figsize=(W_TEXT, 4.4))
sizes = [r"$N=1024$" "\n" r"(76$\times$1)", r"$N=4096$" "\n" r"(76$\times$1)"]
xb = np.arange(len(sizes)); w = 0.36
cheby_vs_master = [tps_1024["master"][0]/tps_1024["chebyshev"][0],
                   tps_4096["master"][0]/tps_4096["chebyshev"][0]]
nspline_vs_master = [tps_1024["master"][0]/tps_1024["native-spline"][0],
                     tps_4096["master"][0]/tps_4096["native-spline"][0]]
ax.bar(xb - w/2, nspline_vs_master, w, color=C["native-spline"],
       label=LBL["native-spline"], edgecolor=CAM["slate_4"], linewidth=0.8)
ax.bar(xb + w/2, cheby_vs_master, w, color=C["chebyshev"],
       label=LBL["chebyshev"], edgecolor=CAM["slate_4"], linewidth=0.8)
for i, v in enumerate(nspline_vs_master): ax.text(i - w/2, v+0.12, f"{v:.1f}$\\times$", ha="center", fontsize=10)
for i, v in enumerate(cheby_vs_master):  ax.text(i + w/2, v+0.12, f"{v:.1f}$\\times$", ha="center", fontsize=10)
ax.axhline(1.0, color=CAM["slate_3"], ls="--", lw=1.0)
ax.set_xticks(xb); ax.set_xticklabels(sizes)
ax.set_ylabel("Speedup vs upstream master")
ax.set_title("Throughput vs master at matched pure-MPI config")
ax.legend(frameon=True, framealpha=0.93); ax.grid(True, axis="y")
fig.tight_layout(); save(fig, "fig4_speedup_vs_master")

# ===== combined tidy CSV =====
with open(f"{OUT}/combined_hybrid_mem_30901621.csv", "w") as f:
    f.write("dataset,branch,N_molecules,mpi_ranks,omp_threads,tps_s,aggregate_mem_mib\n")
    decomp_mpi = [76, 38, 19, 4, 2, 1]; decomp_omp = [1, 2, 4, 19, 38, 76]
    for b in tps_1024:
        for i in range(6):
            f.write(f"decomp_sweep,{b},1024,{decomp_mpi[i]},{decomp_omp[i]},{tps_1024[b][i]},{mem_1024[b][i]}\n")
    for b in tps_4096:
        for i in range(6):
            if tps_4096[b][i] is not None:
                f.write(f"decomp_sweep,{b},4096,{decomp_mpi[i]},{decomp_omp[i]},{tps_4096[b][i]},\n")
    for t, s in zip(omp_threads, omp_tps):
        f.write(f"cheby_omp_ladder,chebyshev,1024,1,{t},{s},\n")
    for c, sp, ef in zip(ns_cores, ns_speedup, ns_eff):
        f.write(f"nspline_core_scaling,native-spline,1024,{c},1,,\n")

print("Wrote figures (png+pdf) + CSV to", OUT)
for fn in sorted(os.listdir(OUT)):
    print("  ", fn)
print("\n=== HEADLINE NUMBERS ===")
print(f"N=1024 chebyshev vs master @76x1:   {tps_1024['master'][0]/tps_1024['chebyshev'][0]:.2f}x")
print(f"N=4096 chebyshev vs master @76x1:   {tps_4096['master'][0]/tps_4096['chebyshev'][0]:.2f}x")
print(f"chebyshev best(1x76)={min(tps_1024['chebyshev']):.4f}s vs master best(76x1)={min(tps_1024['master']):.4f}s -> {min(tps_1024['master'])/min(tps_1024['chebyshev']):.2f}x")
print(f"chebyshev OMP speedup @76 threads:  {omp_tps[0]/omp_tps[-1]:.1f}x (eff {100*(omp_tps[0]/omp_tps[-1])/76:.0f}%)")
print(f"memory @ equal speed: cheby 1x76={mem_1024['chebyshev'][-1]} MiB vs nspline 38x2={mem_1024['native-spline'][1]} MiB -> {mem_1024['native-spline'][1]/mem_1024['chebyshev'][-1]:.1f}x less")
