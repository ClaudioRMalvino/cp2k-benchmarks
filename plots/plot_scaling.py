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


def parse_header(path):
    # The run scripts emit `# key: value` lines before the data rows.  Pull
    # them into a dict so titles can show the actual run parameters instead
    # of hardcoded values that drift if STEPS/MPI_RANKS are overridden.
    if path is None:
        return {}
    meta = {}
    with open(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            stripped = line.lstrip("#").strip()
            if ":" in stripped:
                key, _, val = stripped.partition(":")
                meta[key.strip()] = val.strip()
    return meta


# ── Discover and load data ────────────────────────────────────────────────────
data = {}
for name, info in BRANCHES.items():
    d, lbl = info["dir"], info["label"]
    size_pat = f"{RESULTS_ROOT}/{d}/NNP/NNP_size_scaling_{lbl}_*/results_size_scaling_{lbl}_*.csv"
    core_pat = f"{RESULTS_ROOT}/{d}/NNP/NNP_core_scaling_{lbl}_*/results_core_scaling_{lbl}_*.csv"
    size_csv = latest_csv(size_pat)
    core_csv = latest_csv(core_pat)
    data[name] = {
        "size":      load_csv(size_csv, SIZE_COLS),
        "size_meta": parse_header(size_csv),
        "core":      load_csv(core_csv, CORE_COLS),
        "core_meta": parse_header(core_csv),
        "color":     info["color"],
        "marker":    info["marker"],
    }


def consensus(kind, key, fallback="?"):
    # Return the value of `key` from `<kind>_meta` across branches.  Warn if
    # branches disagree (meaning runs were launched with different params and
    # shouldn't really be plotted on the same axes) but proceed with the
    # first non-empty value so plotting still works.
    seen = {}
    for branch, d in data.items():
        v = d[f"{kind}_meta"].get(key)
        if v:
            seen.setdefault(v, []).append(branch)
    if not seen:
        return fallback
    if len(seen) > 1:
        print(f"  WARNING: {kind} '{key}' differs across branches: {seen}")
    return next(iter(seen))

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

size_mpi   = consensus("size", "MPI ranks")
size_omp   = consensus("size", "OMP threads")
size_steps = consensus("size", "steps")
ax.set_title(
    "NNP Size Scaling: Walltime vs Number of Atoms\n"
    f"({size_mpi} MPI ranks, {size_omp} OMP threads, {size_steps} MD steps)"
)
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

core_nmol  = consensus("core", "N_molecules")
core_steps = consensus("core", "steps")
core_atoms = int(core_nmol) * 3 if core_nmol.isdigit() else "?"
ax.set_title(
    "NNP Core Scaling: Walltime vs Cores\n"
    f"({core_nmol} H₂O molecules, {core_atoms} atoms, {core_steps} MD steps)"
)
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

ax.set_title(
    "NNP Core Scaling: Speedup vs Cores\n"
    f"({core_nmol} H₂O molecules, {core_atoms} atoms, {core_steps} MD steps)"
)
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
