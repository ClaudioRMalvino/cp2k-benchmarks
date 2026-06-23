#!/usr/bin/env python3
"""Chebyshev hybrid-MPI/OMP + memory benchmark figures for supervisor meeting.
Data: SLURM job 30901621 (NNP_hybrid_mem), CSD3 icelake, 76 cores/node, H2O-N NNP MD.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.expanduser("~/cp2k-benchmarks/plots/cheby_benchmark_figs")
os.makedirs(OUT, exist_ok=True)

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

colors = {"master": "#777777", "native-spline": "#1f77b4", "chebyshev": "#d62728"}
mark   = {"master": "s", "native-spline": "o", "chebyshev": "^"}

# ===== FIG 1: hybrid decomposition sweep, t/step, N=1024 (headline) =====
fig, ax = plt.subplots(figsize=(7.2, 5.0))
for b, y in tps_1024.items():
    ax.plot(x, y, mark[b]+"-", color=colors[b], lw=2, ms=8, label=b)
ax.set_yscale("log")
ax.set_xticks(x); ax.set_xticklabels(decomp_labels)
ax.set_xlabel("Decomposition  (MPI ranks x OMP threads,  76 cores fixed)")
ax.set_ylabel("Time per MD step (s, log scale)")
ax.set_title("Hybrid MPI/OMP scaling, 1024 water (N=3072 atoms), 1 node")
ax.axvspan(-0.3, 0.3, color="navy", alpha=0.05)
ax.annotate("pure MPI", (0, ax.get_ylim()[0]), color="navy", ha="center", va="bottom", fontsize=8)
ax.annotate("pure OMP", (5, ax.get_ylim()[0]), color="navy", ha="center", va="bottom", fontsize=8)
ax.grid(True, which="both", ls=":", alpha=0.4)
ax.legend(title="branch", frameon=True)
fig.tight_layout(); fig.savefig(f"{OUT}/fig1_hybrid_sweep_N1024.png", dpi=160); plt.close(fig)

# ===== FIG 2: speed-memory tradeoff, N=1024 =====
fig, ax = plt.subplots(figsize=(7.2, 5.0))
for b in tps_1024:
    ax.plot(tps_1024[b], mem_1024[b], mark[b]+"-", color=colors[b], lw=1.5, ms=8, label=b, alpha=0.9)
# annotate chebyshev sweet spot (1x76) and native-spline best (38x2)
ax.annotate("cheby 1x76\n0.057 s, 476 MiB", (0.057025, 476), textcoords="offset points",
            xytext=(40, 10), fontsize=8, color=colors["chebyshev"],
            arrowprops=dict(arrowstyle="->", color=colors["chebyshev"]))
ax.annotate("n-spline 38x2\n0.056 s, 10450 MiB", (0.055525, 10450), textcoords="offset points",
            xytext=(40, -5), fontsize=8, color=colors["native-spline"],
            arrowprops=dict(arrowstyle="->", color=colors["native-spline"]))
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("Time per MD step (s, log scale)  --  faster left")
ax.set_ylabel("Aggregate node memory (MiB, log scale)  --  lower better")
ax.set_title("Speed vs node memory, 1024 water, 1 node\n(each point = one MPIxOMP decomposition)")
ax.grid(True, which="both", ls=":", alpha=0.4)
ax.legend(title="branch"); fig.tight_layout()
fig.savefig(f"{OUT}/fig2_speed_vs_memory_N1024.png", dpi=160); plt.close(fig)

# ===== FIG 3: chebyshev OMP scaling ladder =====
fig, ax = plt.subplots(figsize=(7.2, 5.0))
speedup = omp_tps[0] / omp_tps
ax.plot(omp_threads, speedup, "^-", color=colors["chebyshev"], lw=2, ms=9, label="chebyshev (measured)")
ax.plot(omp_threads, omp_threads, "k--", lw=1, alpha=0.6, label="ideal linear")
for t, s in zip(omp_threads, speedup):
    ax.annotate(f"{s:.1f}x", (t, s), textcoords="offset points", xytext=(6, -2), fontsize=8)
ax.set_xlabel("OpenMP threads (1 MPI rank)")
ax.set_ylabel("Speedup vs 1 thread")
ax.set_title("Chebyshev centre-level OpenMP scaling, 1024 water\n(this is the parallelism master/native-spline lack)")
ax.grid(True, ls=":", alpha=0.4); ax.legend()
fig.tight_layout(); fig.savefig(f"{OUT}/fig3_cheby_omp_scaling.png", dpi=160); plt.close(fig)

# ===== FIG 4: speedup vs master bar chart (pure MPI 76x1, both sizes) =====
fig, ax = plt.subplots(figsize=(7.2, 5.0))
sizes = ["N=1024\n(76x1)", "N=4096\n(76x1)"]
xb = np.arange(len(sizes)); w = 0.35
cheby_vs_master = [tps_1024["master"][0]/tps_1024["chebyshev"][0],
                   tps_4096["master"][0]/tps_4096["chebyshev"][0]]
nspline_vs_master = [tps_1024["master"][0]/tps_1024["native-spline"][0],
                     tps_4096["master"][0]/tps_4096["native-spline"][0]]
ax.bar(xb - w/2, nspline_vs_master, w, color=colors["native-spline"], label="native-spline")
ax.bar(xb + w/2, cheby_vs_master, w, color=colors["chebyshev"], label="chebyshev")
for i, v in enumerate(nspline_vs_master): ax.text(i - w/2, v+0.1, f"{v:.1f}x", ha="center", fontsize=9)
for i, v in enumerate(cheby_vs_master):  ax.text(i + w/2, v+0.1, f"{v:.1f}x", ha="center", fontsize=9)
ax.axhline(1.0, color="k", ls="--", lw=1, alpha=0.6)
ax.set_xticks(xb); ax.set_xticklabels(sizes)
ax.set_ylabel("Speedup vs upstream master (pure MPI)")
ax.set_title("Throughput vs master at matched pure-MPI config")
ax.legend(); ax.grid(True, axis="y", ls=":", alpha=0.4)
fig.tight_layout(); fig.savefig(f"{OUT}/fig4_speedup_vs_master.png", dpi=160); plt.close(fig)

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

print("Wrote figures + CSV to", OUT)
for fn in sorted(os.listdir(OUT)):
    print("  ", fn)

# headline numbers for the talk track
print("\n=== HEADLINE NUMBERS ===")
print(f"N=1024 chebyshev vs master @76x1:   {tps_1024['master'][0]/tps_1024['chebyshev'][0]:.2f}x")
print(f"N=4096 chebyshev vs master @76x1:   {tps_4096['master'][0]/tps_4096['chebyshev'][0]:.2f}x")
print(f"chebyshev best(1x76)={min(tps_1024['chebyshev']):.4f}s vs master best(76x1)={min(tps_1024['master']):.4f}s -> {min(tps_1024['master'])/min(tps_1024['chebyshev']):.2f}x")
print(f"chebyshev OMP speedup @76 threads:  {omp_tps[0]/omp_tps[-1]:.1f}x (eff {100*(omp_tps[0]/omp_tps[-1])/76:.0f}%)")
print(f"memory @ equal speed: cheby 1x76={mem_1024['chebyshev'][-1]} MiB vs nspline 38x2={mem_1024['native-spline'][1]} MiB -> {mem_1024['native-spline'][1]/mem_1024['chebyshev'][-1]:.1f}x less")
