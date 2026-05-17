#!/usr/bin/env python3
"""Consolidated CSD3 / Peta4-IceLake scaling benchmark plotting & analysis.

Error bars are 95% CIs from the per-rep std via Student-t (df = n_reps - 1),
following Hoefler & Belli SC'15 Rule 5.

Usage:
  python3 plot_csd3_scaling.py                     # produce everything (default)
  python3 plot_csd3_scaling.py --only main         # 7 base scaling figures
  python3 plot_csd3_scaling.py --only enhanced     # 4 enhancement outputs
  python3 plot_csd3_scaling.py --only stats        # 2 statistical-rigour outputs
  python3 plot_csd3_scaling.py --only <section>    # a single section by id
"""
import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from scipy import stats as sstats


RESULTS_ROOT  = "/home/crm98/rds/hpc-work/cp2k-benchmarks/results"
PLOT_DIR      = "/home/crm98/cp2k-benchmarks/plots/scaling_csd3"
ATOMS_PER_MOL = 3   # H2O

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
OMP_2D_COLS     = ["omp_threads", "n_molecules", "mpi_ranks", "total_cores", "n_reps",
                   "tps_mean", "tps_std", "tps_min",
                   "wt_mean",  "wt_std",  "wt_min"]

N_CORES_TICKS = [1, 2, 4, 8, 16, 32, 64, 76]


def _ensure_dirs():
    os.makedirs(PLOT_DIR, exist_ok=True)


def latest(pattern):
    # *_raw.csv shares the prefix; exclude it so glob only picks the summary CSV.
    matches = [m for m in glob.glob(pattern) if not m.endswith("_raw.csv")]
    return max(matches, key=os.path.getmtime) if matches else None


def load_csv(path, cols):
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


# 95% CI half-width via Student-t: t(n-1, 0.025) * s / sqrt(n).
# Hoefler & Belli SC'15 Rule 5.
def add_ci95(df, std_col="tps_std", n_col="n_reps", out_col="tps_ci95"):
    if df is None or std_col not in df.columns or n_col not in df.columns:
        return df
    n = df[n_col].astype(float)
    s = df[std_col].astype(float)
    # t.ppf is undefined for df <= 0; mask those rows.
    tcrit = np.where(n > 1, sstats.t.ppf(0.975, np.maximum(n - 1, 1)), np.nan)
    df = df.copy()
    df[out_col] = tcrit * s / np.sqrt(np.maximum(n, 1))
    return df


def save_fig(fig, name):
    os.makedirs(PLOT_DIR, exist_ok=True)
    path = os.path.join(PLOT_DIR, name + ".pdf")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {path}")
    return path


def save_table(df_or_text, basename, header=None):
    """Save DataFrame as .txt + .csv (string: .txt only)."""
    os.makedirs(PLOT_DIR, exist_ok=True)
    txt = os.path.join(PLOT_DIR, basename + ".txt")
    csv = os.path.join(PLOT_DIR, basename + ".csv")
    if isinstance(df_or_text, pd.DataFrame):
        with open(txt, "w") as f:
            if header:
                f.write(header)
            f.write(df_or_text.to_string(index=False))
            f.write("\n")
        df_or_text.to_csv(csv, index=False)
        print(f"  saved {txt}\n         {csv}")
        return txt, csv
    with open(txt, "w") as f:
        if header:
            f.write(header)
        f.write(df_or_text)
    print(f"  saved {txt}")
    return txt, None


def load_all_data():
    print(f"Reading CSVs from {RESULTS_ROOT}")
    data = {}
    for name, b in BRANCHES.items():
        d, lbl = b["dir"], b["label"]
        size = latest(f"{RESULTS_ROOT}/{d}/NNP/NNP_size_scaling_{lbl}_*/results_size_scaling_{lbl}_*.csv")
        core = latest(f"{RESULTS_ROOT}/{d}/NNP/NNP_core_scaling_{lbl}_*/results_core_scaling_{lbl}_*.csv")
        size_df = load_csv(size, SIZE_COLS)
        if size_df is not None:
            size_df = size_df.sort_values("n_molecules").reset_index(drop=True)
            size_df = add_ci95(size_df)
        core_df = add_ci95(load_csv(core, CORE_COLS))
        data[name] = {
            "size":      size_df,
            "size_meta": parse_header(size),
            "size_path": size,
            "core":      core_df,
            "core_meta": parse_header(core),
            "core_path": core,
            **b,
        }

    omp_thread_csv = latest(f"{RESULTS_ROOT}/cp2k_feature_native_spline_omp/NNP/"
                            f"NNP_omp_thread_scaling_*/results_omp_thread_scaling_*.csv")
    omp_2d_csv     = latest(f"{RESULTS_ROOT}/cp2k_feature_native_spline_omp/NNP/"
                            f"NNP_omp_size_scaling_*/results_omp_size_scaling_*.csv")
    omp = {
        "thread":      add_ci95(load_csv(omp_thread_csv, OMP_THREAD_COLS)),
        "thread_meta": parse_header(omp_thread_csv),
        "two_d":       add_ci95(load_csv(omp_2d_csv, OMP_2D_COLS)),
        "two_d_meta":  parse_header(omp_2d_csv),
    }
    return data, omp


def fig_size_scaling(data, omp):
    print("\n[main 1/7] size scaling -- time per MD step vs N_atoms")
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for name, d in data.items():
        s = d["size"]
        if s is None or s.empty:
            continue
        ax.errorbar(atoms(s["n_molecules"]), s["tps_mean"], yerr=s["tps_ci95"],
                    marker=d["marker"], color=d["color"], capsize=3,
                    lw=1.6, ms=6, label=f"{name} ({d['decomp']})")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("number of atoms")
    ax.set_ylabel("time per MD step (s)")
    steps = data["upstream master"]["size_meta"].get("steps", "?")
    ax.set_title("NNP size scaling -- time per MD step vs system size\n"
                 f"(one Peta4-IceLake node = 76 cores, "
                 f"{steps} steps x 5 reps, 95% CI (Student-t))")
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend(loc="upper left")
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    fig.tight_layout()
    save_fig(fig, "size_scaling_time_per_step")


def fig_size_scaling_speedup(data, omp):
    print("\n[main 2/7] size scaling -- optimised branches' speedup over master")
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
    save_fig(fig, "size_scaling_speedup")


def fig_core_scaling_time(data, omp):
    print("\n[main 3/7] strong scaling -- time per MD step vs cores")
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for name in ("upstream master", "native-spline"):
        d = data[name]
        c = d["core"]
        if c is None or c.empty:
            continue
        ax.errorbar(c["total_cores"], c["tps_mean"], yerr=c["tps_ci95"],
                    marker=d["marker"], color=d["color"], capsize=3,
                    lw=1.6, ms=6, label=name)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks(N_CORES_TICKS)
    ax.xaxis.set_major_formatter(ticker.FixedFormatter([str(c) for c in N_CORES_TICKS]))
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
    save_fig(fig, "core_scaling_time_per_step")


def fig_core_scaling_speedup(data, omp):
    print("\n[main 4/7] strong scaling -- speedup vs cores")
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(N_CORES_TICKS, N_CORES_TICKS, "k--", lw=1.2, label="ideal (linear)")
    for name in ("upstream master", "native-spline"):
        d = data[name]
        c = d["core"]
        if c is None or c.empty:
            continue
        ax.plot(c["total_cores"], c["speedup"],
                marker=d["marker"], color=d["color"], lw=1.6, ms=7, label=name)
    ax.set_xscale("log", base=2)
    ax.set_xticks(N_CORES_TICKS)
    ax.xaxis.set_major_formatter(ticker.FixedFormatter([str(c) for c in N_CORES_TICKS]))
    ax.set_xlabel("total cores (MPI ranks)")
    ax.set_ylabel("speedup (relative to 1 core)")
    nmol = data["upstream master"]["core_meta"].get("N_molecules", "?")
    ax.set_title("NNP strong scaling -- speedup vs cores\n"
                 f"(N = {nmol} H2O, pure MPI, baseline = 1 rank)")
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend()
    fig.tight_layout()
    save_fig(fig, "core_scaling_speedup")


def fig_core_scaling_efficiency(data, omp):
    print("\n[main 5/7] strong scaling -- parallel efficiency vs cores")
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
    ax.set_xticks(N_CORES_TICKS)
    ax.xaxis.set_major_formatter(ticker.FixedFormatter([str(c) for c in N_CORES_TICKS]))
    ax.set_xlabel("total cores (MPI ranks)")
    ax.set_ylabel("parallel efficiency (%)")
    ax.set_ylim(0, 110)
    nmol = data["upstream master"]["core_meta"].get("N_molecules", "?")
    ax.set_title("NNP strong scaling -- parallel efficiency\n"
                 f"(N = {nmol} H2O, pure MPI)")
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend()
    fig.tight_layout()
    save_fig(fig, "core_scaling_efficiency")


def fig_omp_thread_scaling(data, omp):
    print("\n[main 6/7] OMP thread scaling -- native-spline-omp, MPI=1")
    t = omp["thread"]
    if t is None or t.empty:
        print("  no OMP thread-scaling data; skipped")
        return
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    axes[0].errorbar(t["omp_threads"], t["tps_mean"], yerr=t["tps_ci95"],
                     marker="s", color="tab:orange",
                     capsize=3, lw=1.6, ms=6)
    axes[0].set_xscale("log", base=2)
    axes[0].set_yscale("log")
    axes[0].set_xlabel("OMP threads")
    axes[0].set_ylabel("time per MD step (s)")
    axes[0].set_title("time per step vs OMP threads")
    axes[0].grid(True, which="both", ls="--", alpha=0.4)

    ax2 = axes[1]
    ax2.plot(t["omp_threads"], t["omp_threads"], "k--", lw=1.2, label="ideal (linear)")
    ax2.plot(t["omp_threads"], t["speedup"], marker="s", color="tab:orange",
             lw=1.6, ms=7, label="native-spline-omp")
    ax2.set_xscale("log", base=2)
    ax2.set_xlabel("OMP threads")
    ax2.set_ylabel("speedup (vs OMP=1)")
    ax2.set_title("OMP speedup")
    ax2.grid(True, which="both", ls="--", alpha=0.4)
    ax2.legend()

    ticks = [int(v) for v in t["omp_threads"]]
    for a in axes:
        a.set_xticks(ticks)
        a.xaxis.set_major_formatter(ticker.FixedFormatter([str(c) for c in ticks]))
    nmol  = omp["thread_meta"].get("N_molecules", "?")
    steps = omp["thread_meta"].get("steps", "?")
    fig.suptitle("OMP thread scaling -- feature/nnp-native-spline-omp  "
                 f"(MPI=1, N={nmol} H2O, {steps} steps x 5 reps)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    save_fig(fig, "omp_thread_scaling")


def fig_omp_size_scaling_2d(data, omp):
    print("\n[main 7/7] OMP size scaling 2D -- t/step vs N_atoms, line per OMP count")
    td = omp["two_d"]
    if td is None or td.empty:
        print("  no OMP 2D scaling data; skipped")
        return
    fig, ax = plt.subplots(figsize=(8, 5.5))
    omps = sorted(td["omp_threads"].unique())
    cmap = plt.cm.viridis(np.linspace(0.15, 0.85, len(omps)))
    for c, omp_val in zip(cmap, omps):
        sub = td[td["omp_threads"] == omp_val].sort_values("n_molecules")
        ax.errorbar(atoms(sub["n_molecules"]), sub["tps_mean"], yerr=sub["tps_ci95"],
                    marker="o", color=c, capsize=3, lw=1.5, ms=5,
                    label=f"OMP = {int(omp_val)}")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("number of atoms")
    ax.set_ylabel("time per MD step (s)")
    steps = omp["two_d_meta"].get("steps", "?")
    ax.set_title("OMP size scaling -- feature/nnp-native-spline-omp\n"
                 f"(MPI=1, sweep OMP x system size, {steps} steps x 5 reps)")
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend(title="OMP threads")
    fig.tight_layout()
    save_fig(fig, "omp_size_scaling_2d")


def fig_size_scaling_with_complexity_refs(data, omp):
    """Size scaling with O(N) and O(N^2) reference lines anchored to master."""
    print("\n[enhanced 1/4] size scaling with complexity reference overlays")
    fig, ax = plt.subplots(figsize=(9, 6))
    for name, d in data.items():
        s = d["size"]
        if s is None:
            continue
        x = atoms(s["n_molecules"])
        ax.errorbar(x, s["tps_mean"], yerr=s["tps_ci95"],
                    marker=d["marker"], color=d["color"], capsize=3, lw=1.6, ms=6,
                    label=f"{name} ({d['decomp']})")
    m = data["upstream master"]["size"]
    x_ref = m["n_molecules"].iloc[0] * ATOMS_PER_MOL
    t_ref = m["tps_mean"].iloc[0]
    x_arr = np.array([x_ref, m["n_molecules"].iloc[-1] * ATOMS_PER_MOL])
    ax.plot(x_arr, t_ref * (x_arr / x_ref), "--", color="gray", lw=1.4,
            label=r"$\propto N$ reference (O$(N)$ scaling)")
    ax.plot(x_arr, t_ref * (x_arr / x_ref) ** 2, ":", color="black", lw=1.4,
            label=r"$\propto N^2$ reference (O$(N^2)$ scaling)")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("number of atoms")
    ax.set_ylabel("time per MD step (s)")
    ax.set_title("NNP size scaling vs algorithmic-complexity references\n"
                 "(one Peta4-IceLake node, 100 steps x 5 reps, 95% CI (Student-t))")
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend(loc="upper left", fontsize=9)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    fig.tight_layout()
    save_fig(fig, "size_scaling_with_complexity_refs")


def write_power_law_fits(data, omp):
    """Weighted log-log linear fit t = a * N^b per branch."""
    print("\n[enhanced 2/4] power-law fit per branch")
    rows = []
    for name, d in data.items():
        s = d["size"]
        if s is None or len(s) < 3:
            continue
        x = (s["n_molecules"] * ATOMS_PER_MOL).to_numpy(dtype=float)
        y = s["tps_mean"].to_numpy(dtype=float)
        yerr = s["tps_std"].to_numpy(dtype=float)
        logx, logy = np.log(x), np.log(y)
        logy_err = np.where(y > 0, yerr / y, 1.0)
        (b, log_a), cov = np.polyfit(logx, logy, 1, w=1.0 / np.maximum(logy_err, 1e-9), cov=True)
        b_err = np.sqrt(cov[0, 0])
        a_err = np.sqrt(cov[1, 1])
        rows.append(dict(branch=name, a=np.exp(log_a), a_err=np.exp(log_a) * a_err,
                         b=b, b_err=b_err, n_points=len(x)))
    fit_df = pd.DataFrame(rows)
    display = fit_df.copy()
    display["exponent (b)"]      = [f"{r.b:.3f} +/- {r.b_err:.3f}" for r in fit_df.itertuples()]
    display["prefactor (a, sec)"] = [f"{r.a:.2e}" for r in fit_df.itertuples()]
    display = display[["branch", "exponent (b)", "prefactor (a, sec)", "n_points"]]
    print()
    print("Fit:  t(N) = a * N^b   (N = number of atoms)")
    print()
    print(display.to_string(index=False))

    header = ("Power-law fit of time-per-step vs system size\n"
              "  Model:   t(N) = a * N^b   (N = number of atoms)\n"
              "  Method:  weighted log-log linear regression, weights = (yerr / y)\n\n")
    footer = ("\n\nInterpretation:\n"
              "  b = 1.0  -> O(N) linear scaling  (cell-list neighbour search)\n"
              "  b = 2.0  -> O(N^2) pair enumeration\n"
              "  b between 1 and 2 -> mixed regime (e.g. communication overhead masking O(N))\n")
    txt, _ = save_table(display, "power_law_fits", header=header)
    with open(txt, "a") as f:
        f.write(footer)
    fit_df.to_csv(os.path.join(PLOT_DIR, "power_law_fits.csv"), index=False)


def fig_size_scaling_per_atom(data, omp):
    """Time per atom per MD step vs N (flat curve = O(N))."""
    print("\n[enhanced 3/4] time per atom per MD step vs N")
    fig, ax = plt.subplots(figsize=(9, 6))
    for name, d in data.items():
        s = d["size"]
        if s is None:
            continue
        n_atoms = s["n_molecules"] * ATOMS_PER_MOL
        t_per_atom = s["tps_mean"] / n_atoms
        t_per_atom_err = s["tps_ci95"] / n_atoms
        ax.errorbar(n_atoms, t_per_atom * 1e6, yerr=t_per_atom_err * 1e6,
                    marker=d["marker"], color=d["color"], capsize=3, lw=1.6, ms=6,
                    label=f"{name} ({d['decomp']})")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("number of atoms")
    ax.set_ylabel(r"time per MD step per atom ($\mu$s)")
    ax.set_title("Per-atom cost: how flat the curve is, is how O(N) the algorithm is\n"
                 "(one Peta4-IceLake node, 76 cores total)")
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend(loc="upper left", fontsize=10)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    fig.tight_layout()
    save_fig(fig, "size_scaling_per_atom_per_step")


def write_summary_table(data, omp):
    """Headline workload summary at N = 1024 H2O."""
    print("\n[enhanced 4/4] headline-workload summary table (N = 1024 H2O, full node)")
    m_t1024 = data["upstream master"]["size"][
        data["upstream master"]["size"]["n_molecules"] == 1024]["tps_mean"].iloc[0]
    rows = []
    for name, d in data.items():
        s = d["size"]
        if s is None:
            continue
        sel = s[s["n_molecules"] == 1024]
        if len(sel) == 0:
            continue
        sel = sel.iloc[0]
        rows.append(dict(
            branch=name,
            decomposition=d["decomp"],
            t_step_mean_s=f"{sel['tps_mean']:.4f}",
            t_step_std_s=f"{sel['tps_std']:.4f}",
            t_step_ci95_s=f"+/-{sel['tps_ci95']:.4f}",
            walltime_s_for_100_steps=f"{sel['wt_mean']:.2f}",
            speedup_over_master=f"{m_t1024 / sel['tps_mean']:.2f}x",
        ))
    table = pd.DataFrame(rows)
    print()
    print(table.to_string(index=False))
    header = "Headline-workload comparison: N = 1024 H2O (3072 atoms), one Peta4-IceLake node\n\n"
    footer = ("\n\nColumn key:\n"
              "  t_step_mean_s            : sample mean of 5 timed reps (warm-up discarded)\n"
              "  t_step_std_s             : sample std deviation of those 5 reps\n"
              "  t_step_ci95_s            : 95% confidence interval half-width via Student-t,\n"
              "                             df = n - 1 = 4, t(4, 0.025) = 2.776 (Hoefler & Belli SC'15 Rule 5)\n"
              "  walltime_s_for_100_steps : end-to-end run wall time for 100 MD steps (mean)\n"
              "  speedup_over_master      : ratio of master mean / branch mean (cf. Rule 1, also see\n"
              "                             statistical_comparison_N1024 for significance testing)\n")
    txt, _ = save_table(table, "size_scaling_summary_table", header=header)
    with open(txt, "a") as f:
        f.write(footer)


# Statistical-rigour sections operate on per-rep raw CSVs (..._raw.csv) rather
# than summary tables, to surface distributions and run significance tests.
HEADLINE_N = 1024


def _raw_path_from_summary(summary_path):
    if summary_path is None:
        return None
    base, ext = os.path.splitext(summary_path)
    cand = base + "_raw" + ext
    return cand if os.path.exists(cand) else None


def load_raw_size_reps(summary_path, n_molecules):
    raw_path = _raw_path_from_summary(summary_path)
    if raw_path is None:
        return None
    df = pd.read_csv(raw_path, comment="#", header=None,
                     names=["n_molecules", "rep", "time_per_step_s", "walltime_s"],
                     na_values=["NA", "FAILED"])
    sel = df[df["n_molecules"] == n_molecules].dropna(subset=["time_per_step_s"])
    return sel["time_per_step_s"].astype(float).to_numpy() if len(sel) else None


def fig_headline_distribution(data, omp):
    """Box + violin of the raw reps at the headline N (Hoefler Rule 8/12)."""
    print(f"\n[stats 1/2] headline-workload raw-rep distribution (box+violin, N = {HEADLINE_N})")
    reps_per_branch = {}
    for name, d in data.items():
        reps = load_raw_size_reps(d.get("size_path"), HEADLINE_N)
        if reps is None or len(reps) == 0:
            print(f"  !! no raw reps for {name} at N={HEADLINE_N}; skipped from distribution plot")
            continue
        reps_per_branch[name] = reps
        print(f"    {name:20s} n_reps={len(reps):d}  mean={reps.mean():.4f}  std={reps.std(ddof=1):.4f}")

    if not reps_per_branch:
        print("  no usable data; skipping section")
        return

    fig, ax = plt.subplots(figsize=(9, 6))
    labels  = list(reps_per_branch.keys())
    samples = [reps_per_branch[k] for k in labels]
    colors  = [data[k]["color"] for k in labels]
    positions = np.arange(len(labels)) + 1

    parts = ax.violinplot(samples, positions=positions, widths=0.6,
                          showmeans=False, showmedians=False, showextrema=False)
    for pc, col in zip(parts['bodies'], colors):
        pc.set_facecolor(col); pc.set_alpha(0.30); pc.set_edgecolor(col); pc.set_linewidth(1.2)

    bp = ax.boxplot(samples, positions=positions, widths=0.25, patch_artist=True,
                    showfliers=True,
                    boxprops=dict(linewidth=1.0),
                    medianprops=dict(color='black', linewidth=1.5),
                    whiskerprops=dict(linewidth=1.0),
                    capprops=dict(linewidth=1.0))
    for patch, col in zip(bp['boxes'], colors):
        patch.set_facecolor(col); patch.set_alpha(0.55); patch.set_edgecolor('black')

    for pos, samp, col in zip(positions, samples, colors):
        jitter = np.random.RandomState(42).uniform(-0.08, 0.08, size=len(samp))
        ax.scatter(pos + jitter, samp, s=28, color=col, edgecolor='black',
                   linewidth=0.5, alpha=0.9, zorder=3)

    ax.set_xticks(positions)
    ax.set_xticklabels([f"{l}\n({data[l]['decomp']})" for l in labels], fontsize=9)
    ax.set_ylabel("time per MD step (s)")
    ax.set_title(f"Headline-workload distribution at N = {HEADLINE_N} H$_2$O ({HEADLINE_N*ATOMS_PER_MOL} atoms)\n"
                 f"Box: 25/50/75 %ile. Violin: kernel density. Dots: individual reps (n = 5)")
    ax.grid(axis='y', ls='--', alpha=0.4)
    fig.tight_layout()
    save_fig(fig, f"headline_distribution_N{HEADLINE_N}")


def write_statistical_comparison(data, omp):
    """Pairwise inter-branch hypothesis tests at the headline N (Hoefler Rule 7)."""
    print(f"\n[stats 2/2] inter-branch significance test (N = {HEADLINE_N})")
    samples = {}
    for name, d in data.items():
        reps = load_raw_size_reps(d.get("size_path"), HEADLINE_N)
        if reps is not None and len(reps) > 1:
            samples[name] = reps
    if len(samples) < 2:
        print("  not enough branches with raw data to compare; skipped")
        return

    rows = []
    names = list(samples.keys())
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            xa, xb = samples[a], samples[b]
            mean_a, mean_b = xa.mean(), xb.mean()
            std_a,  std_b  = xa.std(ddof=1), xb.std(ddof=1)
            n_a,   n_b    = len(xa), len(xb)

            t_stat, p_t = sstats.ttest_ind(xa, xb, equal_var=False)
            u_stat, p_u = sstats.mannwhitneyu(xa, xb, alternative='two-sided')
            _, p_norm_a = sstats.shapiro(xa) if n_a >= 3 else (np.nan, np.nan)
            _, p_norm_b = sstats.shapiro(xb) if n_b >= 3 else (np.nan, np.nan)

            pooled = np.sqrt(((n_a - 1)*std_a**2 + (n_b - 1)*std_b**2) / (n_a + n_b - 2))
            cohens_d = (mean_a - mean_b) / pooled if pooled > 0 else np.nan
            mag = abs(cohens_d)
            interp = ("negligible" if mag < 0.2 else
                      "small"      if mag < 0.5 else
                      "medium"     if mag < 0.8 else
                      "large"      if mag < 1.2 else "very large")

            rows.append(dict(
                pair=f"{a}  vs  {b}",
                mean_a=f"{mean_a:.5f}",
                mean_b=f"{mean_b:.5f}",
                mean_diff_pct=f"{100*(mean_a - mean_b)/mean_b:+.1f}",
                shapiro_p_a=f"{p_norm_a:.3f}",
                shapiro_p_b=f"{p_norm_b:.3f}",
                welch_t=f"{t_stat:+.3f}",
                welch_p=f"{p_t:.2e}",
                mannwhitney_p=f"{p_u:.2e}",
                cohens_d=f"{cohens_d:+.2f}",
                effect=interp,
                sig_at_5pct=("YES" if p_t < 0.05 else "no"),
            ))
    table = pd.DataFrame(rows)
    print()
    print(table.to_string(index=False))

    header = (f"Inter-branch statistical comparison at N = {HEADLINE_N} H2O on one full Ice Lake node\n"
              f"(5 timed reps per branch; warm-up discarded)\n\n")
    footer = ("\n\nColumn key:\n"
              "  mean_a, mean_b   : sample mean of branch A / B (s per MD step)\n"
              "  mean_diff_pct    : (mean_a - mean_b) / mean_b * 100\n"
              "  shapiro_p_a/b    : Shapiro-Wilk p-value for normality of each sample.\n"
              "                     p < 0.05 -> reject normality, prefer Mann-Whitney over Welch's t.\n"
              "                     With only n=5, the Shapiro-Wilk test has low power; results are\n"
              "                     reported for transparency, not as definitive evidence of normality.\n"
              "  welch_t          : Welch's two-sample t-statistic (unequal variances allowed)\n"
              "  welch_p          : two-sided p-value of Welch's test (parametric, assumes ~normal)\n"
              "  mannwhitney_p    : two-sided Mann-Whitney U p-value (non-parametric backup)\n"
              "  cohens_d         : standardised effect size, (mean_a - mean_b) / pooled_std\n"
              "  effect           : Cohen's conventional bin: small (0.2), medium (0.5), large (0.8),\n"
              "                     very large (1.2). All inter-branch effects here are 'very large'\n"
              "                     because the per-rep std is on the order of 10^-3 s while branch\n"
              "                     differences are 10^-2 s.\n"
              "  sig_at_5pct      : Welch test rejects equal-means null at alpha = 0.05\n"
              "\n"
              "Method (Hoefler & Belli SC'15 Rule 7): for each pair of branches we test the\n"
              "null hypothesis 'mean time per step is equal' using Welch's t-test (parametric,\n"
              "robust to unequal variances), with a non-parametric Mann-Whitney U as a fallback\n"
              "in case the per-rep data is non-normal.  Effect size reported via Cohen's d.\n")
    txt, _ = save_table(table, f"statistical_comparison_N{HEADLINE_N}", header=header)
    with open(txt, "a") as f:
        f.write(footer)


MAIN_SECTIONS = [
    ("size",         "size scaling -- t/step vs N_atoms",            fig_size_scaling),
    ("size_speedup", "size scaling -- speedup over master",          fig_size_scaling_speedup),
    ("core_t",       "strong scaling -- t/step vs cores",            fig_core_scaling_time),
    ("core_sp",      "strong scaling -- speedup vs cores",           fig_core_scaling_speedup),
    ("core_eff",     "strong scaling -- parallel efficiency",        fig_core_scaling_efficiency),
    ("omp_thread",   "OMP thread scaling (MPI=1)",                   fig_omp_thread_scaling),
    ("omp_2d",       "OMP 2D size sweep",                            fig_omp_size_scaling_2d),
]
ENHANCED_SECTIONS = [
    ("complexity_refs", "size scaling + O(N)/O(N^2) ref overlays",   fig_size_scaling_with_complexity_refs),
    ("power_law",       "weighted log-log power-law fit per branch", write_power_law_fits),
    ("per_atom",        "time per atom per MD step vs N",            fig_size_scaling_per_atom),
    ("summary_table",   "headline-workload summary at N=1024 H2O",   write_summary_table),
]
STATS_SECTIONS = [
    ("distribution",    f"box+violin of raw 5 reps at N={HEADLINE_N}",   fig_headline_distribution),
    ("comparison",      f"Welch's t-test + Mann-Whitney U at N={HEADLINE_N}", write_statistical_comparison),
]
ALL_SECTIONS = MAIN_SECTIONS + ENHANCED_SECTIONS + STATS_SECTIONS
SECTION_IDS = [s[0] for s in ALL_SECTIONS]


def main():
    ap = argparse.ArgumentParser(description="Consolidated CSD3 scaling plots + analysis.")
    ap.add_argument("--only",
                    help="Run a subset: 'main' (7 base figs), 'enhanced' (4 enhancements), "
                         "'stats' (2 statistical-rigour outputs), or a single section id ("
                         + ", ".join(SECTION_IDS) + ")")
    args = ap.parse_args()

    _ensure_dirs()
    data, omp = load_all_data()

    if args.only == "main":
        sections = MAIN_SECTIONS
    elif args.only == "enhanced":
        sections = ENHANCED_SECTIONS
    elif args.only == "stats":
        sections = STATS_SECTIONS
    elif args.only:
        sections = [s for s in ALL_SECTIONS if s[0] == args.only]
        if not sections:
            raise SystemExit(f"unknown section id {args.only!r}; "
                             f"valid: main, enhanced, stats, or one of {SECTION_IDS}")
    else:
        sections = ALL_SECTIONS

    for _, _, fn in sections:
        fn(data, omp)

    print(f"\nDone. Outputs in {PLOT_DIR}")


if __name__ == "__main__":
    main()
