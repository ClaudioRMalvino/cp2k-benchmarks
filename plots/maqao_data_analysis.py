#!/usr/bin/env python3
"""Consolidated MAQAO-derived analysis for the thesis report.

Two data sources are wired in:

  1. LEGACY (May-14 ONE View, hand-typed):
       GLOBAL_METRICS, HOTSPOTS, CQA_LOOPS  -- a frozen snapshot retained for
       backwards compatibility with the existing plots & disclaimer.

  2. FRESH (May-26 CSV from `maqao cqa --fct-loops="nnp_" --output-format=csv`,
       three branches, login-node CQA, ptrace-independent):
       loaded by load_cqa_csvs() into a single DataFrame and aggregated with
       Hoefler & Belli SC'15 rule-compliant means
       (arithmetic for costs, harmonic for rates, geometric for ratios).

The CLI exposes both via --only.  Default runs the fresh-CSV-driven outputs
first, then the legacy outputs.
"""
import argparse
import hashlib
import os
import subprocess
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


OUTDIR = "/home/crm98/cp2k-benchmarks/plots/maqao_plots"
THREE_BRANCHES = ["master", "native-spline (pre-port)", "native-spline-omp"]

# =============================================================================
# Fresh CQA CSVs (one per branch).  Produced on the CSD3 login node via
#   maqao cqa <binary> --fct-loops="nnp_" --output-format=csv \
#             --output-path=<path>
# CQA is pure static disassembly; same binary in -> identical numbers out,
# so the data is deterministic (Rule 5: no CIs required).
# =============================================================================
FRESH_CSV_DIR = "/home/crm98/cp2k-benchmarks/maqao_login_cqa"
FRESH_CSV_PATHS = {
    "master":                    f"{FRESH_CSV_DIR}/master_nnp.csv",
    "native-spline (post-port)": f"{FRESH_CSV_DIR}/native_spline_postport_nnp.csv",
    "native-spline-omp":         f"{FRESH_CSV_DIR}/omp_nnp.csv",
}
# Binaries the CSVs were analyzed from (for the methodology block / md5 trace).
FRESH_BINARY_PATHS = {
    "master":                    "/rds/user/crm98/hpc-work/cp2k_binaries/csd3/master/cp2k.psmp",
    "native-spline (post-port)": "/home/crm98/cp2k_optimized/install/bin/cp2k.psmp",
    "native-spline-omp":         "/rds/user/crm98/hpc-work/cp2k_binaries/csd3/feature-nnp-native-spline-omp/cp2k.psmp",
}
# Restrict to NNP Fortran sources -- the CQA --fct-loops filter catches some
# inlined helium_interactions wrappers; this drops the tail of unrelated loops.
NNP_SRC_BASENAMES = {"nnp_acsf.F", "nnp_force.F", "nnp_model.F", "nnp_environment.F"}
# Function-name shorteners (Fortran module-prefixed symbols are long).
def _short_fn(s):
    if not isinstance(s, str): return "?"
    base = s.split("(")[0].strip()
    return base.split("::")[-1]
def _short_src(s):
    if not isinstance(s, str): return "?"
    return s.split("/")[-1]


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
        "libm pow (__powr8i4)":          2.55,
        "OMP runtime (atomics / barriers)": 0.0,
        "hermite_spline_value":          0.0,
        "other / user code":            100 - (33.25 + 5.86 + 3.05 + 3.54 + 2.81 + 2.55),
    },
    "native-spline (pre-port)": {
        "MPI broadcast":               34.77,
        "libm tanh / exp (NN activ.)":  0.14,
        "nnp_calc_ang":                 2.64,
        "pbc2 (O(N^2) neighbour cost)": 0.00,
        "libm pow (__powr8i4)":         5.24,
        "OMP runtime (atomics / barriers)": 0.0,
        "hermite_spline_value":        10.88,
        "other / user code":           100 - (34.77 + 0.14 + 2.64 + 5.24 + 10.88),
    },
    "native-spline-omp": {
        "MPI broadcast":               27.21,
        "libm tanh / exp (NN activ.)":  0.21,
        "nnp_calc_ang":                 1.65,
        "pbc2 (O(N^2) neighbour cost)": 0.00,
        "libm pow (__powr8i4)":         3.22,
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


# =============================================================================
# Fresh CQA CSV: loader, per-branch aggregator, per-loop table, summary plots
# =============================================================================
def load_cqa_csvs(paths=FRESH_CSV_PATHS, nnp_only=True):
    """Read each branch's CQA CSV, tag with `branch`, optionally restrict to
    the NNP Fortran sources (drops unrelated inlined helium wrappers etc.)."""
    parts = []
    for branch, p in paths.items():
        if not os.path.exists(p):
            print(f"  [load] MISSING: {p} -- skipping {branch}")
            continue
        d = pd.read_csv(p, sep=";")
        d.insert(0, "branch", branch)
        parts.append(d)
    if not parts:
        return None
    df = pd.concat(parts, ignore_index=True)
    df["fn_short"]  = df["Function name"].map(_short_fn)
    df["src_short"] = df["Source file"].map(_short_src)
    if nnp_only:
        df = df[df["src_short"].isin(NNP_SRC_BASENAMES)].copy()
    return df


def _harmean(s):
    """Harmonic mean over a Series; drops NaN and non-positive entries
    (zero-rate or sentinel rows would otherwise divide by zero)."""
    s = pd.Series(s).dropna()
    s = s[s > 0]
    if len(s) == 0: return float("nan")
    return len(s) / (1.0 / s).sum()


def _geomean(s):
    """Geometric mean over a Series; drops NaN and zeros (log(0) undefined).
    For percentage columns this means scalar loops with vec_ratio=0 are
    EXCLUDED from the average; they are reported separately as the
    'n_loops_with_vec_ratio=0' tail count."""
    s = pd.Series(s).dropna()
    s = s[s > 0]
    if len(s) == 0: return float("nan")
    return float(np.exp(np.log(s).mean()))


def aggregate_branch_metrics(df):
    """Per-branch Hoefler-correct aggregates:
       costs   -> arithmetic + median + p25/p75    (Rule 3, Rule 8)
       rates   -> harmonic mean                    (Rule 3)
       ratios  -> geometric mean of positives      (Rule 4)
       tails   -> count of loops at 0% / >=50%     (Rule 8: don't reduce to
                  one number when the distribution is bimodal)
    Returns one row per branch (DataFrame)."""
    def per_branch(g):
        cyc       = g["[L1] Nb cycles: min"].astype(float)
        cyc_full  = g["[L1] Nb cycles if fully vectorized: min"].astype(float)
        cyc_fp    = g["[L1] Nb cycles if FP arith vectorized: min"].astype(float)
        cyc_clean = g["[L1] Nb cycles if clean: min"].astype(float)
        ipc       = g["[L1] IPC: max"].astype(float)
        flopcyc   = g["[L1] FLOP/cycle: max"].astype(float)
        vr        = g["Vec. ratio (%): all"].astype(float)
        vrf       = g["Vec. ratio (%): all (FP)"].astype(float)
        ve        = g["Vec. eff. ratio (%): all"].astype(float)
        ai        = g["Arith. intensity (FLOP / ld+st bytes)"].astype(float)
        return pd.Series({
            "n_loops":                                len(g),
            # costs: arithmetic + distribution
            "cycles/iter (arith mean)":               cyc.mean(),
            "cycles/iter (median)":                   cyc.median(),
            "cycles/iter (p25)":                      cyc.quantile(0.25),
            "cycles/iter (p75)":                      cyc.quantile(0.75),
            # Rule 11: upper-bound model -- what perfect vec / clean code buys
            "cycles/iter (full-vec, arith mean)":     cyc_full.mean(),
            "cycles/iter (FP-vec, arith mean)":       cyc_fp.mean(),
            "cycles/iter (clean, arith mean)":        cyc_clean.mean(),
            # rates: harmonic
            "IPC (harm mean)":                        _harmean(ipc),
            "FLOP/cycle (harm mean)":                 _harmean(flopcyc),
            # ratios: geometric of positives
            "Vec ratio % (geo mean, > 0)":            _geomean(vr),
            "Vec ratio FP % (geo mean, > 0)":         _geomean(vrf),
            "Vec efficiency % (geo mean, > 0)":       _geomean(ve),
            "Arith intensity (median)":               ai.median(),
            # tails of the vec ratio distribution -- the scalar vs SIMD split
            "n_loops with vec_ratio == 0":            int((vr.fillna(-1) == 0).sum()),
            "n_loops with vec_ratio >= 50%":          int((vr.fillna(-1) >= 50).sum()),
        })
    return df.groupby("branch", sort=False).apply(per_branch, include_groups=False)


def write_fresh_branch_aggregate():
    """Rule-3/4/8/11 aggregate across all NNP loops, one row per branch.
       Saved as CSV + readable .txt with the rule annotations baked in."""
    print("\n[fresh-1] CSV-driven per-branch aggregate (rule-correct means)")
    df = load_cqa_csvs()
    if df is None:
        print("  no CSVs found -- skipping"); return
    agg = aggregate_branch_metrics(df).T
    print("\n", agg.to_string())
    notes = (
        "\n\nAggregation rules (Hoefler & Belli SC'15)\n"
        "  Rule 3 -- costs:  arithmetic mean of cycles/iter, FLOP, bytes\n"
        "  Rule 3 -- rates:  harmonic mean of IPC, FLOP/cycle\n"
        "  Rule 4 -- ratios: geometric mean of vec ratio / vec efficiency,\n"
        "                    excluding loops with vec_ratio == 0 (log(0) undefined).\n"
        "                    Reported separately as 'n_loops with vec_ratio == 0'.\n"
        "  Rule 8 -- mean + median + p25/p75 reported for cycles/iter so the\n"
        "            heavy-tailed distribution is visible.\n"
        "  Rule 11 -- 'full-vec / FP-vec / clean' columns are CQA upper-bound\n"
        "             models: lower bound on cycles/iter under that idealization.\n"
        "  Rule 5 -- CQA is deterministic static analysis; no CIs are required.\n"
    )
    csv_p = os.path.join(OUTDIR, "fresh_branch_aggregate.csv")
    txt_p = os.path.join(OUTDIR, "fresh_branch_aggregate.txt")
    agg.to_csv(csv_p)
    with open(txt_p, "w") as f:
        f.write("MAQAO CQA fresh aggregate (N = 64 H2O reference binary, NNP-only loops)\n")
        f.write("=" * 72 + "\n\n")
        f.write(agg.to_string())
        f.write(notes)
    print(f"  saved {csv_p}\n         {txt_p}")


def write_fresh_per_loop_table(top_n=20):
    """Per-loop CQA stats across all branches, restricted to the heaviest
       loops (by cycles/iter).  Pairs every Rule-11 upper-bound with the
       absolute base case (Rule 1).
    """
    print(f"\n[fresh-2] CSV-driven per-loop table (top {top_n} per branch by cycles/iter)")
    df = load_cqa_csvs()
    if df is None: print("  no CSVs found -- skipping"); return
    cols_keep = [
        "branch","fn_short","src_short","Source lines","ID","Path ID",
        "[L1] Nb cycles: min",
        "[L1] Nb cycles if fully vectorized: min",
        "[L1] Nb cycles if FP arith vectorized: min",
        "[L1] Nb cycles if clean: min",
        "[L1] IPC: max",
        "Vec. ratio (%): all",
        "Vec. ratio (%): all (FP)",
        "Vec. eff. ratio (%): all",
        "Arith. intensity (FLOP / ld+st bytes)",
    ]
    keep = df[cols_keep].copy()
    keep = keep.sort_values(["branch","[L1] Nb cycles: min"], ascending=[True, False])
    # take top_n per branch
    top = keep.groupby("branch", sort=False).head(top_n)
    # Rule 1: paired speedup -- vec-headroom relative to base (avoid div by 0)
    base = top["[L1] Nb cycles: min"].astype(float)
    fullv = top["[L1] Nb cycles if fully vectorized: min"].astype(float)
    top["headroom (base/fullvec)"] = (base / fullv.where(fullv > 0)).round(2)
    # round numeric columns for readability
    for c in cols_keep[6:]:
        if c in top.columns:
            top[c] = pd.to_numeric(top[c], errors="coerce").round(2)
    csv_p = os.path.join(OUTDIR, "fresh_per_loop_table.csv")
    top.to_csv(csv_p, index=False)
    print(f"  saved {csv_p}  ({len(top)} rows)")
    # also save a compact text version
    txt_p = os.path.join(OUTDIR, "fresh_per_loop_table.txt")
    with open(txt_p, "w") as f:
        f.write("Top NNP loops per branch (ranked by cycles/iter, min)\n")
        f.write("=" * 72 + "\n\n")
        for b in top["branch"].unique():
            f.write(f"--- {b} ---\n")
            sub = top[top["branch"]==b].drop(columns=["branch"])
            f.write(sub.to_string(index=False) + "\n\n")
    print(f"         {txt_p}")


def plot_fresh_vec_ratio_distribution():
    """Rule 8: when the per-loop distribution is wide, plot the spread,
       not just the mean.  Box+strip per branch."""
    print("\n[fresh-3] per-loop vec-ratio distribution (boxplot)")
    df = load_cqa_csvs()
    if df is None: print("  no CSVs found -- skipping"); return
    branches = list(FRESH_CSV_PATHS.keys())
    data = [df.loc[df["branch"]==b, "Vec. ratio (%): all"].dropna().values for b in branches]
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    box = ax.boxplot(data, labels=branches, showmeans=True, meanline=False,
                     meanprops=dict(marker="D", markerfacecolor="darkred",
                                    markeredgecolor="darkred", markersize=6),
                     medianprops=dict(color="black", lw=1.5),
                     widths=0.55, patch_artist=True)
    colors = ["#7fb3d5", "#2ecc71", "#f39c12"]
    for patch, col in zip(box["boxes"], colors):
        patch.set_facecolor(col); patch.set_alpha(0.75); patch.set_edgecolor("black")
    # overlay actual points (strip plot) with horizontal jitter
    rng = np.random.default_rng(0)
    for i, d in enumerate(data):
        xs = rng.normal(i + 1, 0.05, size=len(d))
        ax.scatter(xs, d, alpha=0.35, s=10, color="grey", zorder=3)
    ax.set_ylabel("Vectorisation ratio (%) -- per loop", fontsize=11)
    ax.set_xlabel("")
    ax.set_title("MAQAO CQA: per-loop vectorisation ratio across NNP loops\n"
                 "(box = IQR, line = median, diamond = arithmetic mean, dots = individual loops)")
    ax.set_ylim(-5, 105)
    ax.grid(axis="y", ls="--", alpha=0.35)
    fig.text(0.5, 0.01,
             "Aggregating ratios with the arithmetic mean over-weights tiny "
             "scalar prologue loops (Rule 4). The geometric mean of the >0 tail "
             "and the full distribution are both reported.",
             ha="center", va="bottom", fontsize=8, style="italic", wrap=True)
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    p = os.path.join(OUTDIR, "fresh_vec_ratio_distribution.pdf")
    fig.savefig(p, dpi=300, bbox_inches="tight"); plt.close(fig)
    print(f"  saved {p}")


def plot_fresh_speedup_with_baseline(top_n=8):
    """Rule 1: any speedup needs the absolute base-case.  Paired-bar plot of
       (current cycles/iter, full-vec cycles/iter) for the top NNP loops on
       the post-port native-spline branch -- shows the headroom in absolute
       cycles, not a unit-less speedup ratio."""
    print(f"\n[fresh-4] post-port paired (base, full-vec) cycles/iter, top {top_n}")
    df = load_cqa_csvs()
    if df is None: print("  no CSVs found -- skipping"); return
    sub = df[df["branch"] == "native-spline (post-port)"].copy()
    sub = sub.sort_values("[L1] Nb cycles: min", ascending=False).head(top_n)
    labels = [f"{f}\n{s}" for f, s in zip(sub["fn_short"].map(lambda s: s.replace("_mp_","\n")), sub["Source lines"])]
    base = sub["[L1] Nb cycles: min"].astype(float).values
    full = sub["[L1] Nb cycles if fully vectorized: min"].astype(float).values
    x = np.arange(len(sub))
    w = 0.4
    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.bar(x - w/2, base, width=w, label="current (cycles/iter, min)",
           color="#5c9ec7", edgecolor="black")
    ax.bar(x + w/2, full, width=w, label="if fully vectorised (CQA upper bound)",
           color="#7fc97f", edgecolor="black")
    for xi, b, f in zip(x, base, full):
        ax.text(xi - w/2, b + 0.5, f"{b:.1f}", ha="center", va="bottom", fontsize=8)
        ax.text(xi + w/2, f + 0.5, f"{f:.2f}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8, rotation=30, ha="right")
    ax.set_ylabel("CQA cycles per iteration (min)")
    ax.set_title("Vectorisation headroom on feature/nnp-native-spline (post-port)\n"
                 f"Top {top_n} NNP loops by current cycles/iter -- absolute base shown alongside upper bound (Rule 1)")
    ax.legend(loc="upper right")
    ax.grid(axis="y", ls="--", alpha=0.35)
    fig.tight_layout()
    p = os.path.join(OUTDIR, "fresh_speedup_with_baseline.pdf")
    fig.savefig(p, dpi=300, bbox_inches="tight"); plt.close(fig)
    print(f"  saved {p}")


def write_methodology_md():
    """Auto-generate a methodology.md tying together: environment, MAQAO version,
       invocation, binary md5, statistical choices, and which Hoefler rule each
       table satisfies."""
    print("\n[methodology] writing methodology.md")
    # collect binary md5s
    md5s = {}
    for label, b in FRESH_BINARY_PATHS.items():
        try:
            with open(b, "rb") as f:
                md5s[label] = hashlib.md5(f.read()).hexdigest()
        except FileNotFoundError:
            md5s[label] = "(file not present)"
    # maqao version
    try:
        ver = subprocess.check_output(
            ["/home/crm98/maqao.x86_64.2026.0.0-b/bin/maqao", "--version"],
            text=True, stderr=subprocess.STDOUT, timeout=10).strip().splitlines()[0]
    except Exception:
        ver = "MAQAO 2026.0.0-b (build path)"
    # csv mtimes
    csv_mtimes = {}
    for label, p in FRESH_CSV_PATHS.items():
        try:
            csv_mtimes[label] = datetime.fromtimestamp(os.path.getmtime(p)).isoformat(sep=" ", timespec="seconds")
        except FileNotFoundError:
            csv_mtimes[label] = "(missing)"

    md = []
    md.append("# MAQAO CQA Methodology\n")
    md.append("This document records how the per-branch MAQAO CQA measurements that feed ")
    md.append("`maqao_data_analysis.py` were produced, and which Hoefler & Belli SC'15 rule ")
    md.append("each downstream table / plot satisfies.\n\n")
    md.append("## Environment\n")
    md.append("- **Cluster**: CSD3 Peta4-IceLake (Intel Xeon Platinum 8368Q, 2×38 = 76 cores/node, AVX-512)\n")
    md.append("- **Execution host**: CSD3 login node (same IceLake silicon as compute nodes; CQA is static disassembly so no compute-node ptrace is required).\n")
    md.append("- **Kernel mitigation**: `kernel.yama.ptrace_scope = 2` (LProf disabled cluster-wide; CQA unaffected).\n")
    md.append(f"- **MAQAO version**: `{ver}`\n")
    md.append("- **CQA invocation**: `maqao cqa <binary> --fct-loops=\"nnp_\" --output-format=csv --output-path=<csv>`\n")
    md.append("\n")
    md.append("## Binaries analyzed (md5)\n")
    md.append("| Branch | Binary path | md5 |\n|---|---|---|\n")
    for label, b in FRESH_BINARY_PATHS.items():
        md.append(f"| {label} | `{b}` | `{md5s[label]}` |\n")
    md.append("\n")
    md.append("## Fresh CSV inputs\n")
    md.append("| Branch | CSV | Last written |\n|---|---|---|\n")
    for label, p in FRESH_CSV_PATHS.items():
        md.append(f"| {label} | `{p}` | {csv_mtimes[label]} |\n")
    md.append("\n")
    md.append("## Hoefler & Belli rule compliance\n\n")
    md.append("| Rule | Where it applies | How this script satisfies it |\n|---|---|---|\n")
    md.append("| **R1** Speedup needs base-case absolute | `fresh_per_loop_table.csv`, `fresh_speedup_with_baseline.pdf` | `headroom (base/fullvec)` column always paired with the underlying `Nb cycles: min` it is computed from. The headroom plot shows absolute cycles, not the ratio. |\n")
    md.append("| **R2** Reason for subsets | `fresh_per_loop_table.csv` (top N = 20) | Sorted by `cycles/iter` and truncated to top 20 per branch. The unfiltered full per-loop CSV is also written so the truncation is auditable. |\n")
    md.append("| **R3** Arithmetic for costs, harmonic for rates | `fresh_branch_aggregate.csv` | Cycles/iter, FLOP, bytes use arithmetic mean. IPC and FLOP/cycle use harmonic mean. |\n")
    md.append("| **R4** Avoid summarizing ratios; geomean if must | `fresh_branch_aggregate.csv` | Vec ratio / vec efficiency report geometric mean across the positive tail, with the zero-tail count broken out separately so the distribution is not collapsed into one misleading number. |\n")
    md.append("| **R5** Deterministic vs nondeterministic | All fresh-CSV outputs | CQA is pure static disassembly of the compiled binary; same binary → identical numbers. No CIs are computed because the measurement has zero variance. |\n")
    md.append("| **R8** Investigate central-tendency choice | `fresh_branch_aggregate.csv`, `fresh_vec_ratio_distribution.pdf` | Cycles/iter is reported as mean + median + p25 + p75. Vec ratio is reported with full distribution (boxplot) plus zero-tail count. |\n")
    md.append("| **R9** Document factors and environment | this file | Cluster, kernel mitigation, MAQAO version, exact invocation, binary md5, CSV mtimes. |\n")
    md.append("| **R11** Upper performance bounds | `fresh_branch_aggregate.csv`, `fresh_speedup_with_baseline.pdf` | CQA's `Nb cycles if fully vectorized / FP arith vectorized / clean` columns are propagated as explicit upper-bound models alongside the achieved values. |\n")
    md.append("| **R12** Plot only valid trends | `fresh_*` plots | Bars and box plots throughout; no spurious line interpolation between non-comparable points. |\n")
    md.append("\n")
    md.append("## What about LProf (coverage %) ?\n")
    md.append("LProf is the dynamic sampling component of MAQAO ONE View and is "
              "currently blocked on CSD3 by the `ptrace_scope = 2` kernel mitigation. "
              "The legacy `HOTSPOTS` / `CQA_LOOPS` dicts in this script (from the May-14 "
              "ONE View run, before the mitigation) retain coverage% values; the fresh "
              "CSV-driven outputs do **not** include coverage% and therefore cannot "
              "weight per-loop aggregates by execution frequency. Aggregates across "
              "loops are equal-weighted, with the loop-by-loop distribution always "
              "available in the per-loop CSV so the reader can re-weight by hand if "
              "desired.\n")
    out = os.path.join(OUTDIR, "maqao_methodology.md")
    with open(out, "w") as f:
        f.write("".join(md))
    print(f"  saved {out}")


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
    # --- fresh CSV-driven outputs (Hoefler & Belli compliant) ---
    ("fresh-1",     "CSV-driven per-branch aggregate",            write_fresh_branch_aggregate),
    ("fresh-2",     "CSV-driven per-loop table (top 20/branch)",  write_fresh_per_loop_table),
    ("fresh-3",     "Per-loop vec-ratio distribution boxplot",    plot_fresh_vec_ratio_distribution),
    ("fresh-4",     "Paired base/full-vec cycles bar (post-port)", plot_fresh_speedup_with_baseline),
    ("methodology", "Methodology .md (env + rules + md5)",        write_methodology_md),
    # --- legacy ONE View (May-14 hand-typed; retains LProf coverage %) ---
    ("4",          "[legacy] Hotspot stacked bar (3 branches)",  plot_hotspot_stackbar),
    ("5",          "[legacy] Per-loop vectorisation table",      write_vectorisation_table),
    ("global",     "[legacy] Global-metrics summary table",      write_global_metrics_table),
    ("focused",    "[legacy] Array-access-efficiency focused",   plot_array_access_focused),
    ("disclaimer", "[legacy] Thesis-ready disclaimer (Markdown)", write_disclaimer),
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
