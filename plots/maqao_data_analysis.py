#!/usr/bin/env python3
"""Consolidated MAQAO-derived analysis for the thesis report.

Numbers are from the MAQAO ONE View run of 2026-05-14; subsequent re-runs were
blocked by the CSD3 ptrace mitigation applied on 2026-05-15 (the
feature/nnp-native-spline May-14 measurement predates cache-port commit
7cbc8b3008 and is labelled "(pre-port)" throughout).

Usage:
  python3 maqao_data_analysis.py            # produce everything (default)
  python3 maqao_data_analysis.py --only N   # produce only one section (1..5)
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


OUTDIR = "/home/crm98/cp2k-benchmarks/plots/maqao_plots"
THREE_BRANCHES = ["master", "native-spline (pre-port)", "native-spline-omp"]


GLOBAL_METRICS = pd.DataFrame([
    dict(branch="master",                   application_time_s=31.74,
         array_access_efficiency_pct=52.92, vec_headroom_speedup=1.50,
         runtime_overhead_speedup=1.29,     user_time_pct=50.83),
    dict(branch="native-spline (pre-port)", application_time_s=28.21,
         array_access_efficiency_pct=52.52, vec_headroom_speedup=1.45,
         runtime_overhead_speedup=1.32,     user_time_pct=55.01),
    dict(branch="native-spline-omp",        application_time_s=43.54,
         array_access_efficiency_pct=78.17, vec_headroom_speedup=1.30,
         runtime_overhead_speedup=1.95,     user_time_pct=37.16),
])

HOTSPOTS = {
    "master": {
        "MPI broadcast":                33.25,
        "libm tanh / exp (NN activ.)":   5.86 + 3.05,
        "nnp_calc_ang":                  3.54,
        "pbc2 (O(N^2) neighbour cost)":  2.81,
        "__powr8i4":                     2.55,
        "OMP runtime (atomics / barriers)": 0.0,
        "hermite_spline_value":          0.0,
        "other / user code":            100 - (33.25 + 5.86 + 3.05 + 3.54 + 2.81 + 2.55),
    },
    "native-spline (pre-port)": {
        "MPI broadcast":               34.77,
        "libm tanh / exp (NN activ.)":  0.14,
        "nnp_calc_ang":                 2.64,
        "pbc2 (O(N^2) neighbour cost)": 0.00,
        "__powr8i4":                    5.24,
        "OMP runtime (atomics / barriers)": 0.0,
        "hermite_spline_value":        10.88,
        "other / user code":           100 - (34.77 + 0.14 + 2.64 + 5.24 + 10.88),
    },
    "native-spline-omp": {
        "MPI broadcast":               27.21,
        "libm tanh / exp (NN activ.)":  0.21,
        "nnp_calc_ang":                 1.65,
        "pbc2 (O(N^2) neighbour cost)": 0.00,
        "__powr8i4":                    3.22,
        "OMP runtime (atomics / barriers)":
            25.67 + 1.26 + 1.09 + 0.17 + 0.13 + 0.12,
        "hermite_spline_value":         5.98,
        "other / user code":
            100 - (27.21 + 0.21 + 1.65 + 3.22 + 25.67 + 1.26 + 1.09 + 0.17 + 0.13 + 0.12 + 5.98),
    },
}

CQA_LOOPS = pd.DataFrame([
    dict(loop="nnp_calc_ang inner loop", branch="master",
         coverage_pct=13.8, cqa_cycles_per_iter=52, vec_ratio_pct=8,
         speedup_if_fully_vectorised="8.40x"),
    dict(loop="nnp_calc_ang inner loop", branch="native-spline (pre-port)",
         coverage_pct=21.1, cqa_cycles_per_iter=75, vec_ratio_pct=26,
         speedup_if_fully_vectorised="6.68x"),
    dict(loop="nnp_calc_ang inner loop", branch="native-spline-omp",
         coverage_pct=13.1, cqa_cycles_per_iter=75, vec_ratio_pct=26,
         speedup_if_fully_vectorised="6.68x"),
    dict(loop="nnp_calc_acsf (rad fctn body)", branch="master",
         coverage_pct=5.8,  cqa_cycles_per_iter=28, vec_ratio_pct=1,
         speedup_if_fully_vectorised="9.58x"),
    dict(loop="nnp_calc_acsf (rad fctn body)", branch="native-spline (pre-port)",
         coverage_pct=4.1,  cqa_cycles_per_iter=31, vec_ratio_pct=0,
         speedup_if_fully_vectorised="8.00x"),
    dict(loop="nnp_calc_acsf (rad fctn body)", branch="native-spline-omp",
         coverage_pct=3.8,  cqa_cycles_per_iter=17, vec_ratio_pct=0,
         speedup_if_fully_vectorised="8.00x"),
    dict(loop="nnp_scale_acsf", branch="master",
         coverage_pct=3.2,  cqa_cycles_per_iter=8,  vec_ratio_pct=100,
         speedup_if_fully_vectorised="1.00x"),
    dict(loop="nnp_scale_acsf", branch="native-spline (pre-port)",
         coverage_pct=2.0,  cqa_cycles_per_iter=8,  vec_ratio_pct=100,
         speedup_if_fully_vectorised="1.00x"),
    dict(loop="nnp_scale_acsf", branch="native-spline-omp",
         coverage_pct=0.5,  cqa_cycles_per_iter=1,  vec_ratio_pct=0,
         speedup_if_fully_vectorised="8.00x"),
])


def _ensure_outdir():
    os.makedirs(OUTDIR, exist_ok=True)


def _save_table(df_or_text, basename, header=None):
    """Save a DataFrame or text block as .txt (+ .csv if DataFrame)."""
    txt_path = os.path.join(OUTDIR, basename + ".txt")
    csv_path = os.path.join(OUTDIR, basename + ".csv")
    if isinstance(df_or_text, pd.DataFrame):
        with open(txt_path, "w") as f:
            if header:
                f.write(header)
            f.write(df_or_text.to_string(index=False))
            f.write("\n")
        df_or_text.to_csv(csv_path, index=False)
        return txt_path, csv_path
    with open(txt_path, "w") as f:
        if header:
            f.write(header)
        f.write(df_or_text)
    return txt_path, None


def plot_hotspot_stackbar():
    print("\n[4] hotspot time-fraction stacked bar chart")
    categories = list(next(iter(HOTSPOTS.values())).keys())
    colors = ["#7fb3d5", "#e74c3c", "#2ecc71", "#f39c12",
              "#8e44ad", "#c0392b", "#3498db", "#bdc3c7"]

    fig, ax = plt.subplots(figsize=(9, 6))
    bottoms = np.zeros(len(THREE_BRANCHES))
    x = np.arange(len(THREE_BRANCHES))
    for cat, col in zip(categories, colors):
        vals = np.array([HOTSPOTS[b][cat] for b in THREE_BRANCHES])
        ax.bar(x, vals, bottom=bottoms, label=cat, color=col,
               edgecolor="white", linewidth=0.7)
        for xi, vi, bi in zip(x, vals, bottoms):
            if vi >= 4.0:
                ax.text(xi, bi + vi/2, f"{vi:.1f}%", ha="center", va="center",
                        fontsize=8, color="white" if vi >= 10 else "black")
        bottoms += vals

    ax.set_xticks(x)
    ax.set_xticklabels(THREE_BRANCHES, fontsize=10)
    ax.set_ylabel("share of application time (%)")
    ax.set_title("MAQAO hotspot decomposition by branch\n"
                 "(N = 64 H2O, 16 cores, top-15 leaf functions grouped)")
    ax.set_ylim(0, 105)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9)
    ax.grid(axis="y", ls="--", alpha=0.3)
    fig.tight_layout()
    path = os.path.join(OUTDIR, "maqao_hotspot_stackbar.pdf")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {path}")


def write_vectorisation_table():
    print("\n[5] per-loop CQA vectorisation table")
    print()
    print(CQA_LOOPS.to_string(index=False))
    notes = (
        "\n\nNotes\n"
        "  * coverage_pct: percent of total application time spent inside the loop\n"
        "  * cqa_cycles_per_iter: minimum cycles per loop iteration the assembly can achieve\n"
        "  * vec_ratio_pct: percentage of arithmetic instructions issued as SIMD\n"
        "  * speedup_if_fully_vectorised: optimistic upper-bound speedup if ALL\n"
        "    arithmetic were vectorised; large value -> compiler left vec opportunity on the table\n"
    )
    txt, csv = _save_table(
        CQA_LOOPS, "maqao_vectorisation_table",
        header="MAQAO CQA -- per-hot-loop vectorisation analysis (N = 64 H2O, 16 cores)\n\n",
    )
    with open(txt, "a") as f:
        f.write(notes)
    print(f"  saved {txt}\n         {csv}")


def write_global_metrics_table():
    print("\n[global] MAQAO global-metrics summary")
    print()
    print(GLOBAL_METRICS.to_string(index=False))
    notes = (
        "\n\nColumns\n"
        "  array_access_efficiency_pct : MAQAO static stride-friendliness metric\n"
        "                                (proxy for cache locality, no PMU needed)\n"
        "  vec_headroom_speedup        : 'speedup_if_fully_vectorised' -- lower = closer to peak\n"
        "  runtime_overhead_speedup    : 'speedup_if_perfect_MPI_OMP' -- higher = more runtime cost\n"
        "  user_time_pct               : percent of time in user code (vs MPI / OMP / libm)\n\n"
        "Headline finding\n"
        "  array-access efficiency jumps from 52.9% (master) to 78.2% (omp) due to the\n"
        "  cache-conscious optimisations (precomputed scale_factor, sorted active_atoms,\n"
        "  heap sort).  Pure-MPI native-spline at the May-14 measurement did NOT yet\n"
        "  contain the cache opts (subsequently ported in commit 7cbc8b3008), so its\n"
        "  53% number reflects only the cell-list + native-spline subset.\n"
    )
    txt, csv = _save_table(
        GLOBAL_METRICS, "maqao_global_metrics_table",
        header=("MAQAO global metrics  (N = 64 H2O, 16 cores: 16 MPI for master/native-spline;"
                " 8 MPI x 2 OMP for omp)\n\n"),
    )
    with open(txt, "a") as f:
        f.write(notes)
    print(f"  saved {txt}\n         {csv}")


def plot_array_access_focused():
    print("\n[focused] array-access-efficiency, master vs native-spline-omp")
    labels = ["upstream master\n(pure MPI, no cache opts)",
              "feature/nnp-native-spline-omp\n(cell list + Horner spline + cache opts)"]
    values = [
        float(GLOBAL_METRICS.loc[GLOBAL_METRICS.branch == "master",
                                 "array_access_efficiency_pct"].iloc[0]),
        float(GLOBAL_METRICS.loc[GLOBAL_METRICS.branch == "native-spline-omp",
                                 "array_access_efficiency_pct"].iloc[0]),
    ]
    colors = ["tab:blue", "tab:orange"]
    delta = values[1] - values[0]

    fig, ax = plt.subplots(figsize=(8, 5.5))
    bars = ax.bar(range(len(values)), values, color=colors,
                  edgecolor="black", linewidth=0.8, width=0.55)
    for b, v in zip(bars, values):
        ax.text(b.get_x() + b.get_width()/2, v + 1.2, f"{v:.1f}%",
                ha="center", va="bottom", fontsize=12, fontweight="bold")
    ax.text(0.5, max(values) + 6, f"+{delta:.1f} pp",
            ha="center", va="center", fontsize=14, color="darkgreen", fontweight="bold")
    ax.annotate("", xy=(0.95, max(values) + 3), xytext=(0.05, max(values) + 3),
                arrowprops=dict(arrowstyle="-|>", color="darkgreen", lw=1.5))
    ax.set_xticks(range(len(values)))
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("array-access efficiency (%)", fontsize=11)
    ax.set_ylim(0, 100)
    ax.axhline(100, color="grey", ls=":", lw=1, alpha=0.6)
    ax.text(1.5, 99, "theoretical max", color="grey", fontsize=9,
            ha="left", va="top", style="italic")
    ax.set_title("Cache locality of the optimised NNP code path\n"
                 "(MAQAO static stride-friendliness metric, N = 64 H$_2$O on 16 cores)")
    ax.grid(axis="y", ls="--", alpha=0.35)

    fig.text(0.5, -0.03,
             "MAQAO measurement from 14 May 2026. The post-port feature/nnp-native-spline branch\n"
             "could not be re-measured: CSD3 applied a kernel-level ptrace mitigation on 15 May 2026\n"
             "in response to a Red Hat / Rocky Linux zero-day, blocking MAQAO's sampler. The OMP\n"
             "branch shown is the same cache-relevant algorithm as the post-port pure-MPI branch.",
             ha="center", va="top", fontsize=8, style="italic", wrap=True,
             bbox=dict(facecolor="#f5f5f5", edgecolor="#cccccc", boxstyle="round,pad=0.4"))
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    path = os.path.join(OUTDIR, "array_access_efficiency_focused.pdf")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {path}")

    relative = delta / values[0] * 100
    table = (
        "MAQAO array-access-efficiency comparison\n"
        "========================================\n\n"
        "Branch                                  array_access_efficiency\n"
        "--------------------------------------  -----------------------\n"
        f"upstream master                                          {values[0]:.2f} %\n"
        f"feature/nnp-native-spline-omp                            {values[1]:.2f} %\n"
        "--------------------------------------  -----------------------\n"
        f"Improvement                                              +{delta:.2f} pp ({relative:.1f}% relative)\n\n"
        "Workload   : 64 H2O molecules (192 atoms), NNP MD, 2000 MD steps\n"
        "Parallelism: 16 cores total per branch (16x1 MPI/OMP for master, 8x2 for omp)\n"
        "Measured   : 2026-05-14 with MAQAO ONE View, --create-report=one (R1)\n"
        "Metric     : array_access_efficiency from shared/run_0/global_metrics.csv\n"
        "             (MAQAO's static stride-friendliness score, derived from CQA\n"
        "              analysis of the hot loops' assembly -- a proxy for cache locality\n"
        "              that does NOT require hardware performance counters.)\n"
    )
    txt, _ = _save_table(table, "array_access_efficiency_focused")
    print(f"  saved {txt}")


def write_disclaimer():
    print("\n[disclaimer] thesis-ready ptrace-mitigation disclaimer paragraph")
    disclaimer = (
        "**Disclaimer (CSD3 ptrace mitigation).** The MAQAO measurements presented in\n"
        "this report were collected on 14 May 2026, prior to a kernel-level ptrace\n"
        "mitigation applied across CSD3 on 15 May 2026 in response to a Red Hat / Rocky\n"
        "Linux zero-day security advisory (HPC support communication, S. Rankin, 15 May\n"
        "2026). The mitigation prevents `ptrace(PTRACE_ATTACH, ...)` -- the syscall\n"
        "MAQAO's LProf sampler relies on -- from being used by unprivileged processes,\n"
        "and a permanent fix awaits a Red Hat / Rocky Linux update at the time of\n"
        "writing. As a consequence, MAQAO could not be re-run against the post-port\n"
        "build of `feature/nnp-native-spline`. The cache-locality and hotspot\n"
        "measurements reported below for \"cache-conscious\" code are therefore taken\n"
        "from `feature/nnp-native-spline-omp`, which already contained the cache-conscious\n"
        "optimisations (precomputed `scale_factor`/`inv_sigma_factor`, sorted\n"
        "`active_atoms` via in-place heap sort) at the time of the May 14 measurement,\n"
        "and is algorithmically identical to the post-port pure-MPI branch on the\n"
        "cache-relevant code paths (only the OpenMP threading layer differs).\n\n"
        "The static `array_access_efficiency` metric used here is derived by MAQAO's\n"
        "CQA module from inspection of load/store stride patterns at the assembly level\n"
        "and is independent of the blocked hardware performance counters; the\n"
        "absolute number is therefore unaffected by the mitigation and is directly\n"
        "comparable across the two branches shown.\n"
    )
    path = os.path.join(OUTDIR, "array_access_efficiency_disclaimer.md")
    with open(path, "w") as f:
        f.write(disclaimer)
    print(f"  saved {path}")


SECTIONS = [
    ("4",          "Hotspot stacked bar (3 branches)",          plot_hotspot_stackbar),
    ("5",          "Per-loop vectorisation table",              write_vectorisation_table),
    ("global",     "Global-metrics summary table",              write_global_metrics_table),
    ("focused",    "Array-access-efficiency focused comparison", plot_array_access_focused),
    ("disclaimer", "Thesis-ready disclaimer (Markdown)",        write_disclaimer),
]


def main():
    p = argparse.ArgumentParser(description="Consolidated MAQAO analysis script.")
    p.add_argument("--only", help="Run only one section: " + ", ".join(s[0] for s in SECTIONS))
    args = p.parse_args()

    _ensure_outdir()
    for key, _, fn in SECTIONS:
        if args.only and args.only != key:
            continue
        fn()
    print(f"\nDone. Outputs in {OUTDIR}")


if __name__ == "__main__":
    main()
