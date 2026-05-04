#!/usr/bin/env python3
import glob
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd

RESULTS_ROOT = "/home/raid/crm98/cp2k-benchmarks/results"
# Plots are written to local scratch (no per-user quota).
# After the job syncs CSVs to home, run this script from cerberus1/athena;
# retrieve PNGs with: scp athena:/local/data/public/crm98/cp2k-benchmarks/plots/*.png .
PLOT_DIR     = "/local/data/public/crm98/cp2k-benchmarks/plots"
os.makedirs(PLOT_DIR, exist_ok=True)

BRANCHES = {
    "Upstream Master": {
        "dir":    "cp2k_master",
        "label":  "upstream-master",
        "color":  "tab:blue",
        "marker": "o",
    },
    "NNP Verlet Cells": {
        "dir":    "cp2k_feature_verlet_cells",
        "label":  "feature-nnp-verlet-cells",
        "color":  "tab:orange",
        "marker": "s",
    },
    "NNP Native Spline": {
        "dir":    "cp2k_feature_native_spline",
        "label":  "feature-nnp-native-spline",
        "color":  "tab:green",
        "marker": "^",
    },
}

SIZE_COLS = ["n_molecules", "walltime_s"]
CORE_COLS = ["mpi_ranks", "omp_threads", "total_cores", "walltime_s", "speedup"]


def latest_csv(pattern):
    matches = glob.glob(pattern, recursive=True)
    if not matches:
        print(f"  WARNING: no files matched: {pattern}")
        return None
    return max(matches, key=os.path.getmtime)


def load_csv(path, columns):
    if path is None:
        return None
    df = pd.read_csv(path, comment="#", header=None, names=columns)
    print(f"  Loaded: {path}")
    return df


# ── Discover and load data ────────────────────────────────────────────────────
data = {}
for name, info in BRANCHES.items():
    d, lbl = info["dir"], info["label"]
    size_pat = f"{RESULTS_ROOT}/{d}/NNP/NNP_size_scaling_{lbl}_*/results_size_scaling_{lbl}_*.csv"
    core_pat = f"{RESULTS_ROOT}/{d}/NNP/NNP_core_scaling_{lbl}_*/results_core_scaling_{lbl}_*.csv"
    data[name] = {
        "size":   load_csv(latest_csv(size_pat), SIZE_COLS),
        "core":   load_csv(latest_csv(core_pat), CORE_COLS),
        "color":  info["color"],
        "marker": info["marker"],
    }

# H₂O: 3 atoms per molecule
N_ATOMS   = [64*3, 256*3, 512*3, 1024*3, 2048*3, 4096*3]
N_CORES   = [1, 2, 4, 8, 16, 32]

exact_int  = ticker.FuncFormatter(lambda x, _: f"{int(x):,}")
exact_core = ticker.FixedFormatter([str(c) for c in N_CORES])


# ── Plot 1: Size Scaling — Walltime vs Number of Atoms ────────────────────────
print("\nPlotting size scaling (walltime vs atoms)...")
fig, ax = plt.subplots(figsize=(9, 6))

for name, d in data.items():
    if d["size"] is not None:
        atoms = d["size"]["n_molecules"] * 3
        ax.plot(atoms, d["size"]["walltime_s"],
                marker=d["marker"], color=d["color"],
                label=name, linewidth=1.8, markersize=6)

ax.set_title("NNP Size Scaling: Walltime vs Number of Atoms\n(36 MPI ranks, 1 OMP thread)")
ax.set_xlabel("Number of Atoms")
ax.set_ylabel("Walltime (s)")
ax.set_xscale("log", base=2)
ax.set_yscale("log")
ax.set_xticks(N_ATOMS)
ax.xaxis.set_major_formatter(exact_int)
ax.grid(True, which="both", ls="--", alpha=0.5)
ax.legend(framealpha=0.9)
fig.tight_layout()
out = os.path.join(PLOT_DIR, "size_scaling_walltime.png")
fig.savefig(out, dpi=300, bbox_inches="tight")
print(f"  Saved: {out}")
plt.close(fig)


# ── Plot 2: Core Scaling — Walltime vs Cores ─────────────────────────────────
print("\nPlotting core scaling (walltime vs cores)...")
fig, ax = plt.subplots(figsize=(9, 6))

for name, d in data.items():
    if d["core"] is not None:
        ax.plot(d["core"]["total_cores"], d["core"]["walltime_s"],
                marker=d["marker"], color=d["color"],
                label=name, linewidth=1.8, markersize=6)

ax.set_title("NNP Core Scaling: Walltime vs Cores\n(64 H₂O molecules, 192 atoms)")
ax.set_xlabel("Total Cores (MPI ranks)")
ax.set_ylabel("Walltime (s)")
ax.set_xscale("log", base=2)
ax.set_yscale("log")
ax.set_xticks(N_CORES)
ax.xaxis.set_major_formatter(ticker.FixedFormatter([str(c) for c in N_CORES]))
ax.grid(True, which="both", ls="--", alpha=0.5)
ax.legend(framealpha=0.9)
fig.tight_layout()
out = os.path.join(PLOT_DIR, "core_scaling_walltime.png")
fig.savefig(out, dpi=300, bbox_inches="tight")
print(f"  Saved: {out}")
plt.close(fig)


# ── Plot 3: Core Scaling — Speedup vs Cores ──────────────────────────────────
print("\nPlotting core scaling (speedup vs cores)...")
fig, ax = plt.subplots(figsize=(9, 6))

for name, d in data.items():
    if d["core"] is not None:
        ax.plot(d["core"]["total_cores"], d["core"]["speedup"],
                marker=d["marker"], color=d["color"],
                label=name, linewidth=1.8, markersize=6)

ax.plot(N_CORES, N_CORES, "k--", linewidth=1.2, label="Ideal (linear)")

ax.set_title("NNP Core Scaling: Speedup vs Cores\n(64 H₂O molecules, 192 atoms)")
ax.set_xlabel("Total Cores (MPI ranks)")
ax.set_ylabel("Speedup (relative to 1 core)")
ax.set_xticks(N_CORES)
ax.xaxis.set_major_formatter(ticker.FixedFormatter([str(c) for c in N_CORES]))
ax.grid(True, ls="--", alpha=0.5)
ax.legend(framealpha=0.9)
fig.tight_layout()
out = os.path.join(PLOT_DIR, "core_scaling_speedup.png")
fig.savefig(out, dpi=300, bbox_inches="tight")
print(f"  Saved: {out}")
plt.close(fig)


print("\nAll plots saved to:", PLOT_DIR)
