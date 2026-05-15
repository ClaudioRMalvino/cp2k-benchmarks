#!/usr/bin/env python3
"""
CSD3 / Peta4-IceLake scaling plots --
upstream master  vs  feature/nnp-native-spline  vs  feature/nnp-native-spline-omp.

Reads the rep-aware CSVs produced by the scripts under
scripts/CSD3_benchmark_scripts/scaling/ (the *_raw.csv files are ignored
here -- the summary CSVs already carry mean / std / min over reps).

Figures written to  plots/scaling_csd3/  :
  size_scaling_time_per_step.png    t / MD step vs N_atoms, three branches
  size_scaling_speedup.png          native-spline & native-spline-omp speedup over master
  core_scaling_time_per_step.png    strong scaling: t / step vs cores  (master, native-spline)
  core_scaling_speedup.png          speedup vs cores + ideal-linear reference
  core_scaling_efficiency.png       parallel efficiency vs cores
  omp_thread_scaling.png            OMP=1..76 at MPI=1, N=64: t/step + speedup
  omp_size_scaling_2d.png           t/step vs N_atoms, line per OMP thread count
"""
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

RESULTS_ROOT = "/home/crm98/rds/hpc-work/cp2k-benchmarks/results"
PLOT_DIR     = "/home/crm98/cp2k-benchmarks/plots/scaling_csd3"
os.makedirs(PLOT_DIR, exist_ok=True)

BRANCHES = {
    "upstream master":   dict(dir="cp2k_master",
                              label="upstream-master",
                              decomp="76 MPI x 1 OMP",
                              color="tab:blue",   marker="o"),
    "native-spline":     dict(dir="cp2k_feature_native_spline",
                              label="feature-nnp-native-spline",
                              decomp="76 MPI x 1 OMP",
                              color="tab:green",  marker="^"),
    "native-spline-omp": dict(dir="cp2k_feature_native_spline_omp",
                              label="feature-nnp-native-spline-omp",
                              decomp="38 MPI x 2 OMP",
                              color="tab:orange", marker="s"),
}

SIZE_COLS = ["n_molecules", "n_reps",
             "tps_mean", "tps_std", "tps_min",
             "wt_mean",  "wt_std",  "wt_min"]
CORE_COLS = ["mpi_ranks", "omp_threads", "total_cores", "n_reps",
             "tps_mean", "tps_std", "tps_min",
             "wt_mean",  "wt_std",  "wt_min",
             "speedup", "efficiency"]
OMP_THREAD_COLS = ["omp_threads", "mpi_ranks", "total_cores", "n_reps",
                   "tps_mean", "tps_std", "tps_min",
                   "wt_mean",  "wt_std",  "wt_min",
                   "speedup", "efficiency"]
OMP_2D_COLS = ["omp_threads", "n_molecules", "mpi_ranks", "total_cores", "n_reps",
               "tps_mean", "tps_std", "tps_min",
               "wt_mean",  "wt_std",  "wt_min"]

ATOMS_PER_MOL = 3   # H2O


# --- helpers ---------------------------------------------------------------
def latest(pattern):
    # Filter out the per-rep *_raw.csv files -- they share the same directory
    # and the same prefix as the summary CSV, so a naive glob picks them up.
    matches = [m for m in glob.glob(pattern) if not m.endswith("_raw.csv")]
    return max(matches, key=os.path.getmtime) if matches else None


def load(path, cols):
    if path is None or not os.path.exists(path):
        return None
    df = pd.read_csv(path, comment="#", header=None, names=cols, na_values=["NA"])
    print(f"  loaded {os.path.basename(path)} ({len(df)} rows)")
    return df


def parse_header(path):
    meta = {}
    if path is None:
        return meta
    with open(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            s = line.lstrip("#").strip()
            if ":" in s:
                k, _, v = s.partition(":")
                meta[k.strip()] = v.strip()
    return meta


def atoms(n_mol):
    return n_mol * ATOMS_PER_MOL


def save(fig, name):
    path = os.path.join(PLOT_DIR, name + ".pdf")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {path}")


# --- discover --------------------------------------------------------------
print(f"Reading CSVs from {RESULTS_ROOT}")
data = {}
for name, b in BRANCHES.items():
    d, lbl = b["dir"], b["label"]
    size = latest(f"{RESULTS_ROOT}/{d}/NNP/NNP_size_scaling_{lbl}_*/results_size_scaling_{lbl}_*.csv")
    core = latest(f"{RESULTS_ROOT}/{d}/NNP/NNP_core_scaling_{lbl}_*/results_core_scaling_{lbl}_*.csv")
    data[name] = {
        "size":      load(size, SIZE_COLS),
        "size_meta": parse_header(size),
        "core":      load(core, CORE_COLS),
        "core_meta": parse_header(core),
        **b,
    }

omp_thread_csv = latest(f"{RESULTS_ROOT}/cp2k_feature_native_spline_omp/NNP/"
                        f"NNP_omp_thread_scaling_*/results_omp_thread_scaling_*.csv")
omp_2d_csv     = latest(f"{RESULTS_ROOT}/cp2k_feature_native_spline_omp/NNP/"
                        f"NNP_omp_size_scaling_*/results_omp_size_scaling_*.csv")
omp_thread     = load(omp_thread_csv, OMP_THREAD_COLS)
omp_thread_meta = parse_header(omp_thread_csv)
omp_2d         = load(omp_2d_csv, OMP_2D_COLS)
omp_2d_meta    = parse_header(omp_2d_csv)


# --- Plot 1: size scaling -- time per MD step vs N_atoms --------------------
print("\n[1/7] size scaling -- time per MD step vs N_atoms")
fig, ax = plt.subplots(figsize=(8, 5.5))
for name, d in data.items():
    s = d["size"]
    if s is None or s.empty:
        continue
    ax.errorbar(atoms(s["n_molecules"]), s["tps_mean"], yerr=s["tps_std"],
                marker=d["marker"], color=d["color"], capsize=3,
                lw=1.6, ms=6, label=f"{name} ({d['decomp']})")
ax.set_xscale("log", base=2)
ax.set_yscale("log")
ax.set_xlabel("number of atoms")
ax.set_ylabel("time per MD step (s)")
steps_label = data["upstream master"]["size_meta"].get("steps", "?")
ax.set_title("NNP size scaling -- time per MD step vs system size\n"
             f"(one Peta4-IceLake node = 76 cores, "
             f"{steps_label} steps x 5 reps, error bars = std)")
ax.grid(True, which="both", ls="--", alpha=0.4)
ax.legend(loc="upper left")
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
fig.tight_layout()
save(fig, "size_scaling_time_per_step")


# --- Plot 2: native-spline / omp speedup over master per size ---------------
print("[2/7] size scaling -- optimised branches' speedup over master")
fig, ax = plt.subplots(figsize=(8, 5.5))
m = data["upstream master"]["size"]
if m is not None:
    for name in ("native-spline", "native-spline-omp"):
        d = data[name]
        s = d["size"]
        if s is None:
            continue
        merged = pd.merge(m[["n_molecules", "tps_mean"]],
                          s[["n_molecules", "tps_mean"]],
                          on="n_molecules", suffixes=("_m", "_o"))
        speedup = merged["tps_mean_m"] / merged["tps_mean_o"]
        ax.plot(atoms(merged["n_molecules"]), speedup,
                marker=d["marker"], color=d["color"], lw=1.8, ms=7,
                label=f"{name} / master")
    ax.axhline(1.0, color="grey", ls="--", lw=1, label="master (= 1x)")
ax.set_xscale("log", base=2)
ax.set_xlabel("number of atoms")
ax.set_ylabel("speedup over upstream master  (master t/step  /  branch t/step)")
ax.set_title("Speedup of the optimised branches over upstream master\n"
             "(per-system-size, one Peta4-IceLake node)")
ax.grid(True, which="both", ls="--", alpha=0.4)
ax.legend()
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
fig.tight_layout()
save(fig, "size_scaling_speedup")


# --- Plot 3: strong scaling -- time per MD step vs cores --------------------
# omp branch's core_scaling completed only 1 point, so this curve is master +
# native-spline (the apples-to-apples pure-MPI comparison).
print("[3/7] strong scaling -- time per MD step vs cores")
fig, ax = plt.subplots(figsize=(8, 5.5))
N_CORES = [1, 2, 4, 8, 16, 32, 64, 76]
for name in ("upstream master", "native-spline"):
    d = data[name]
    c = d["core"]
    if c is None or c.empty:
        continue
    ax.errorbar(c["total_cores"], c["tps_mean"], yerr=c["tps_std"],
                marker=d["marker"], color=d["color"], capsize=3,
                lw=1.6, ms=6, label=name)
ax.set_xscale("log", base=2)
ax.set_yscale("log")
ax.set_xticks(N_CORES)
ax.xaxis.set_major_formatter(ticker.FixedFormatter([str(c) for c in N_CORES]))
ax.set_xlabel("total cores (MPI ranks)")
ax.set_ylabel("time per MD step (s)")
nmol  = data["upstream master"]["core_meta"].get("N_molecules", "?")
steps = data["upstream master"]["core_meta"].get("steps", "?")
ax.set_title("NNP strong scaling -- time per MD step vs cores\n"
             f"(N = {nmol} H2O = {int(nmol)*3 if str(nmol).isdigit() else '?'} atoms, "
             f"{steps} steps x 5 reps, pure MPI)")
ax.grid(True, which="both", ls="--", alpha=0.4)
ax.legend()
fig.tight_layout()
save(fig, "core_scaling_time_per_step")


# --- Plot 4: strong scaling -- speedup vs cores -----------------------------
print("[4/7] strong scaling -- speedup vs cores")
fig, ax = plt.subplots(figsize=(8, 5.5))
ax.plot(N_CORES, N_CORES, "k--", lw=1.2, label="ideal (linear)")
for name in ("upstream master", "native-spline"):
    d = data[name]
    c = d["core"]
    if c is None or c.empty:
        continue
    ax.plot(c["total_cores"], c["speedup"],
            marker=d["marker"], color=d["color"], lw=1.6, ms=7, label=name)
ax.set_xscale("log", base=2)
ax.set_xticks(N_CORES)
ax.xaxis.set_major_formatter(ticker.FixedFormatter([str(c) for c in N_CORES]))
ax.set_xlabel("total cores (MPI ranks)")
ax.set_ylabel("speedup (relative to 1 core)")
ax.set_title("NNP strong scaling -- speedup vs cores\n"
             f"(N = {nmol} H2O, pure MPI, baseline = 1 rank)")
ax.grid(True, which="both", ls="--", alpha=0.4)
ax.legend()
fig.tight_layout()
save(fig, "core_scaling_speedup")


# --- Plot 5: parallel efficiency vs cores -----------------------------------
print("[5/7] strong scaling -- parallel efficiency vs cores")
fig, ax = plt.subplots(figsize=(8, 5.5))
ax.axhline(100, color="k", ls="--", lw=1.2, label="ideal (100%)")
for name in ("upstream master", "native-spline"):
    d = data[name]
    c = d["core"]
    if c is None or c.empty:
        continue
    ax.plot(c["total_cores"], c["efficiency"],
            marker=d["marker"], color=d["color"], lw=1.6, ms=7, label=name)
ax.set_xscale("log", base=2)
ax.set_xticks(N_CORES)
ax.xaxis.set_major_formatter(ticker.FixedFormatter([str(c) for c in N_CORES]))
ax.set_xlabel("total cores (MPI ranks)")
ax.set_ylabel("parallel efficiency (%)")
ax.set_ylim(0, 110)
ax.set_title("NNP strong scaling -- parallel efficiency\n"
             f"(N = {nmol} H2O, pure MPI)")
ax.grid(True, which="both", ls="--", alpha=0.4)
ax.legend()
fig.tight_layout()
save(fig, "core_scaling_efficiency")


# --- Plot 6: OMP thread scaling (native-spline-omp only) --------------------
print("[6/7] OMP thread scaling -- native-spline-omp, MPI=1")
if omp_thread is not None and not omp_thread.empty:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    # left: t / step vs OMP
    axes[0].errorbar(omp_thread["omp_threads"], omp_thread["tps_mean"],
                     yerr=omp_thread["tps_std"], marker="s", color="tab:orange",
                     capsize=3, lw=1.6, ms=6)
    axes[0].set_xscale("log", base=2)
    axes[0].set_yscale("log")
    axes[0].set_xlabel("OMP threads")
    axes[0].set_ylabel("time per MD step (s)")
    axes[0].set_title("time per step vs OMP threads")
    axes[0].grid(True, which="both", ls="--", alpha=0.4)
    # right: speedup vs OMP
    ax2 = axes[1]
    ax2.plot(omp_thread["omp_threads"], omp_thread["omp_threads"],
             "k--", lw=1.2, label="ideal (linear)")
    ax2.plot(omp_thread["omp_threads"], omp_thread["speedup"],
             marker="s", color="tab:orange", lw=1.6, ms=7,
             label="native-spline-omp")
    ax2.set_xscale("log", base=2)
    ax2.set_xlabel("OMP threads")
    ax2.set_ylabel("speedup (vs OMP=1)")
    ax2.set_title("OMP speedup")
    ax2.grid(True, which="both", ls="--", alpha=0.4)
    ax2.legend()
    omp_ticks = [int(v) for v in omp_thread["omp_threads"]]
    for a in axes:
        a.set_xticks(omp_ticks)
        a.xaxis.set_major_formatter(ticker.FixedFormatter([str(c) for c in omp_ticks]))
    nmol_omp  = omp_thread_meta.get("N_molecules", "?")
    steps_omp = omp_thread_meta.get("steps", "?")
    fig.suptitle("OMP thread scaling -- feature/nnp-native-spline-omp  "
                 f"(MPI=1, N={nmol_omp} H2O, {steps_omp} steps x 5 reps)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    save(fig, "omp_thread_scaling")
else:
    print("  no OMP thread-scaling data")


# --- Plot 7: OMP 2D -- t / step vs N_atoms, line per OMP count --------------
print("[7/7] OMP size scaling 2D -- t/step vs N_atoms, line per OMP count")
if omp_2d is not None and not omp_2d.empty:
    fig, ax = plt.subplots(figsize=(8, 5.5))
    omps = sorted(omp_2d["omp_threads"].unique())
    cmap = plt.cm.viridis(np.linspace(0.15, 0.85, len(omps)))
    for c, omp in zip(cmap, omps):
        sub = omp_2d[omp_2d["omp_threads"] == omp].sort_values("n_molecules")
        ax.errorbar(atoms(sub["n_molecules"]), sub["tps_mean"],
                    yerr=sub["tps_std"],
                    marker="o", color=c, capsize=3, lw=1.5, ms=5,
                    label=f"OMP = {int(omp)}")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("number of atoms")
    ax.set_ylabel("time per MD step (s)")
    steps_2d = omp_2d_meta.get("steps", "?")
    ax.set_title("OMP size scaling -- feature/nnp-native-spline-omp\n"
                 f"(MPI=1, sweep OMP x system size, {steps_2d} steps x 5 reps)")
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend(title="OMP threads")
    fig.tight_layout()
    save(fig, "omp_size_scaling_2d")
else:
    print("  no OMP 2D scaling data")


print(f"\nDone. PNGs in {PLOT_DIR}")
