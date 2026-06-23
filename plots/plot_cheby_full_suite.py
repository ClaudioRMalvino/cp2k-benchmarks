#!/usr/bin/env python3
"""Complete chebyshev benchmark figure suite for supervisor meeting.
Sources (all CSD3 icelake, H2O-N NNP MD, 76 cores/node):
  size scaling   : chebyshev/master 76x1 (jobs 20-06 / 10-06); native-spline 76x1 (15-05)
  core scaling   : chebyshev + master SAME JOB 30565007 (20-06); native-spline job 30901621
  multi-node     : chebyshev + master SAME JOB 30565008 (14-06)
"""
import os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.expanduser("~/cp2k-benchmarks/results/cheby_benchmark_figs")
os.makedirs(OUT, exist_ok=True)
C = {"master": "#777777", "native-spline": "#1f77b4", "chebyshev": "#d62728"}
M = {"master": "s", "native-spline": "o", "chebyshev": "^"}

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

# ===== FIG 5: size scaling (2 panels: t/step + speedup vs master) =====
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5))
for b in ["master", "native-spline", "chebyshev"]:
    a1.plot(Nsz, sz[b], M[b]+"-", color=C[b], lw=2, ms=8, label=b)
a1.set_xscale("log", base=2); a1.set_yscale("log")
a1.set_xticks(Nsz); a1.set_xticklabels(Nsz)
a1.set_xlabel("System size (water molecules)"); a1.set_ylabel("Time per MD step (s)")
a1.set_title("Size scaling, 76 MPI x 1 OMP, 1 node"); a1.grid(True, which="both", ls=":", alpha=0.4); a1.legend()
for b in ["native-spline", "chebyshev"]:
    a2.plot(Nsz, sz["master"]/sz[b], M[b]+"-", color=C[b], lw=2, ms=8, label=b)
a2.axhline(1, color="k", ls="--", lw=1, alpha=0.6)
a2.set_xscale("log", base=2); a2.set_xticks(Nsz); a2.set_xticklabels(Nsz)
a2.set_xlabel("System size (water molecules)"); a2.set_ylabel("Speedup vs master")
a2.set_title("Speedup over master grows with system size"); a2.grid(True, ls=":", alpha=0.4); a2.legend()
for n, b in zip(Nsz, sz["master"]/sz["chebyshev"]):
    a2.annotate(f"{b:.1f}x", (n, b), textcoords="offset points", xytext=(0, 6), fontsize=8, ha="center", color=C["chebyshev"])
fig.tight_layout(); fig.savefig(f"{OUT}/fig5_size_scaling.png", dpi=160); plt.close(fig)

# ===== FIG 6: multi-node (t/step + parallel efficiency) =====
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5))
for N in (1024, 4096):
    for b in ("master", "chebyshev"):
        ls = "-" if N == 1024 else "--"
        a1.plot(nodes, mn[N][b], M[b]+ls, color=C[b], lw=2, ms=8, label=f"{b} N={N}")
a1.set_yscale("log"); a1.set_xticks(nodes); a1.set_xticklabels([f"{n}\n({n*76})" for n in nodes])
a1.set_xlabel("Nodes (MPI ranks)"); a1.set_ylabel("Time per MD step (s)")
a1.set_title("Multi-node strong scaling (chebyshev vs master)"); a1.grid(True, which="both", ls=":", alpha=0.4); a1.legend(fontsize=8)
for N in (1024, 4096):
    for b in ("master", "chebyshev"):
        eff = 100 * mn[N][b][0] / (mn[N][b] * nodes)
        ls = "-" if N == 1024 else "--"
        a2.plot(nodes, eff, M[b]+ls, color=C[b], lw=2, ms=8, label=f"{b} N={N}")
a2.axhline(100, color="k", ls=":", lw=1); a2.set_xticks(nodes); a2.set_xticklabels(nodes)
a2.set_xlabel("Nodes"); a2.set_ylabel("Parallel efficiency (%)")
a2.set_title("Multi-node parallel efficiency"); a2.grid(True, ls=":", alpha=0.4); a2.legend(fontsize=8)
fig.tight_layout(); fig.savefig(f"{OUT}/fig6_multinode.png", dpi=160); plt.close(fig)

# ===== FIG 7: core scaling (strong scaling speedup + efficiency), pure MPI N=1024 =====
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5))
for b in ("master", "native-spline", "chebyshev"):
    sp = cs[b][0] / cs[b]
    a1.plot(cores, sp, M[b]+"-", color=C[b], lw=2, ms=7, label=b)
a1.plot(cores, cores, "k--", lw=1, alpha=0.5, label="ideal")
a1.set_xlabel("MPI ranks (1 OMP each)"); a1.set_ylabel("Speedup vs 1 core")
a1.set_title("Pure-MPI strong scaling, N=1024"); a1.grid(True, ls=":", alpha=0.4); a1.legend()
for b in ("master", "native-spline", "chebyshev"):
    eff = 100 * cs[b][0] / (cs[b] * cores)
    a2.plot(cores, eff, M[b]+"-", color=C[b], lw=2, ms=7, label=b)
a2.axhline(100, color="k", ls=":", lw=1)
a2.set_xlabel("MPI ranks"); a2.set_ylabel("Parallel efficiency (%)")
a2.set_title("Pure-MPI parallel efficiency, N=1024"); a2.grid(True, ls=":", alpha=0.4); a2.legend()
fig.tight_layout(); fig.savefig(f"{OUT}/fig7_core_scaling.png", dpi=160); plt.close(fig)

# ===== combined CSV =====
with open(f"{OUT}/combined_full_suite.csv", "w") as f:
    f.write("dataset,branch,N_molecules,nodes,mpi_ranks,omp_threads,tps_s\n")
    for b in sz:
        for n, t in zip(Nsz, sz[b]): f.write(f"size_scaling,{b},{n},1,76,1,{t}\n")
    for b in cs:
        for c, t in zip(cores, cs[b]): f.write(f"core_scaling,{b},1024,1,{c},1,{t}\n")
    for N in mn:
        for b in mn[N]:
            for nd, t in zip(nodes, mn[N][b]): f.write(f"multinode,{b},{N},{nd},{nd*76},1,{t}\n")

print("Wrote suite to", OUT)
for fn in sorted(os.listdir(OUT)):
    if fn.endswith(".png") or fn.endswith(".csv"): print("  ", fn)
print("\n=== KEY NUMBERS ===")
print("Size scaling cheby-vs-master:", {int(n): round(sz['master'][i]/sz['chebyshev'][i],2) for i,n in enumerate(Nsz)})
print("Multinode N=4096 cheby vs master @4 nodes:", round(mn[4096]['master'][2]/mn[4096]['chebyshev'][2],2), "x")
print("Core-scaling serial (1 core) cheby vs nspline:", round(cs['chebyshev'][0]/cs['native-spline'][0],2),
      "x slower (minimax overhead at OMP=1)")
