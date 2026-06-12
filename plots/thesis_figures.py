#!/usr/bin/env python3
"""Multi-panel composite figures for the MPhil report.

Uses the University of Cambridge brand palette + LaTeX-friendly serif font.

Usage:
  python3 thesis_figures.py                # produce all 5
  python3 thesis_figures.py --only N       # produce only fig N (1..5)
"""
import argparse
import glob
import os
import re
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from scipy import stats as sstats


# Cambridge brand palette (https://www.cam.ac.uk/brand-resources/colours)
CAMBRIDGE = {
    "blue_light":   "#D1F9F1",
    "blue":         "#8EE8D8",
    "blue_warm":    "#00BDB6",
    "blue_dark":    "#133844",
    "crest_light":  "#FFE2C8",
    "crest_warm":   "#FFC392",
    "crest":        "#FD8153",
    "crest_dark":   "#DD3025",
    "cherry":       "#CD3572",
    "indigo":       "#5366E0",
    "indigo_warm":  "#B0B9F1",
    "indigo_dark":  "#29347A",
    "green":        "#4DB78C",
    "green_dark":   "#13553A",
    "purple":       "#A368DF",
    "slate_1":      "#ECEEF1",
    "slate_2":      "#B5BDC8",
    "slate_3":      "#546072",
    "slate_4":      "#232830",
    "white":        "#FFFFFF",
}

# Master = neutral slate, native-spline = Cambridge Warm Blue, omp = Crest
# orange (blue vs orange is the safest pair for common forms of colour-blindness).
BRANCH_STYLE = {
    "upstream master":   dict(color=CAMBRIDGE["slate_3"],  marker="o", label="master"),
    "native-spline":     dict(color=CAMBRIDGE["blue_warm"], marker="^", label="feature/nnp-native-spline"),
    "native-spline-omp": dict(color=CAMBRIDGE["crest"],     marker="s", label="feature/nnp-native-spline-omp"),
}


plt.rcParams.update({
    "font.family":        "serif",
    "font.serif":         ["DejaVu Serif", "Liberation Serif", "Bitstream Vera Serif", "serif"],
    "font.size":          12,
    "axes.labelsize":     13,
    "axes.titlesize":     13,
    "xtick.labelsize":    11,
    "ytick.labelsize":    11,
    "legend.fontsize":    11,
    "figure.titlesize":   14,
    "text.color":         CAMBRIDGE["slate_4"],
    "axes.edgecolor":     CAMBRIDGE["slate_4"],
    "axes.labelcolor":    CAMBRIDGE["slate_4"],
    "xtick.color":        CAMBRIDGE["slate_4"],
    "ytick.color":        CAMBRIDGE["slate_4"],
    "axes.linewidth":     0.8,
    "axes.grid":          True,
    "grid.color":         CAMBRIDGE["slate_2"],
    "grid.linestyle":     "--",
    "grid.linewidth":     0.5,
    "grid.alpha":         0.55,
    "lines.linewidth":    1.7,
    "lines.markersize":   6.5,
    "lines.markeredgewidth": 0.8,
    "mathtext.fontset":   "cm",
    "figure.dpi":         100,
    "savefig.dpi":        300,
    "savefig.bbox":       "tight",
    "savefig.pad_inches": 0.05,
    "pdf.fonttype":       42,    # TrueType (avoids font subset issues)
    "ps.fonttype":        42,
})


RESULTS_ROOT  = "/home/crm98/rds/hpc-work/cp2k-benchmarks/results"
OUT_DIR       = "/home/crm98/cp2k-benchmarks/plots/thesis_figures"
ATOMS_PER_MOL = 3   # H2O

# PhDThesisPSnPDF textwidth (a4paper, 12pt body, default margins) is ~6.0".
W_TEXT = 6.3

BRANCH_DIRS = {
    "upstream master":   ("cp2k_master",                    "upstream-master"),
    "native-spline":     ("cp2k_feature_native_spline",     "feature-nnp-native-spline"),
    "native-spline-omp": ("cp2k_feature_native_spline_omp", "feature-nnp-native-spline-omp"),
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


def _latest(pattern):
    matches = [m for m in glob.glob(pattern) if not m.endswith("_raw.csv")]
    return max(matches, key=os.path.getmtime) if matches else None


def _load_csv(path, cols):
    if path is None or not os.path.exists(path):
        return None, None
    df = pd.read_csv(path, comment="#", header=None, names=cols, na_values=["NA"])
    return df, path


def _add_ci95(df):
    if df is None or "tps_std" not in df.columns or "n_reps" not in df.columns:
        return df
    n, s = df["n_reps"].astype(float), df["tps_std"].astype(float)
    tcrit = np.where(n > 1, sstats.t.ppf(0.975, np.maximum(n - 1, 1)), np.nan)
    df = df.copy()
    df["tps_ci95"] = tcrit * s / np.sqrt(np.maximum(n, 1))
    return df


def load_all_data():
    data = {}
    for branch, (d, lbl) in BRANCH_DIRS.items():
        size_csv = _latest(f"{RESULTS_ROOT}/{d}/NNP/NNP_size_scaling_{lbl}_*/results_size_scaling_{lbl}_*.csv")
        core_csv = _latest(f"{RESULTS_ROOT}/{d}/NNP/NNP_core_scaling_{lbl}_*/results_core_scaling_{lbl}_*.csv")
        size_df, _ = _load_csv(size_csv, SIZE_COLS)
        if size_df is not None:
            size_df = size_df.sort_values("n_molecules").reset_index(drop=True)
            size_df = _add_ci95(size_df)
        core_df, _ = _load_csv(core_csv, CORE_COLS)
        core_df = _add_ci95(core_df)
        data[branch] = dict(size=size_df, size_path=size_csv, core=core_df, **BRANCH_STYLE[branch])
    omp_thread_csv = _latest(f"{RESULTS_ROOT}/cp2k_feature_native_spline_omp/NNP/"
                              f"NNP_omp_thread_scaling_*/results_omp_thread_scaling_*.csv")
    # Merge ALL omp_size_scaling CSVs into one dataframe.  The 2D OMP x N
    # grid was filled in across multiple submissions (OMP=1,2,4,8 originally;
    # OMP=16 added later as a focused re-run), so a single _latest() picks
    # only one batch and drops the others.  Concatenating + de-duplicating
    # on (omp_threads, n_molecules) keeps the union.
    omp_2d_paths = sorted(
        p for p in glob.glob(f"{RESULTS_ROOT}/cp2k_feature_native_spline_omp/NNP/"
                              f"NNP_omp_size_scaling_*/results_omp_size_scaling_*.csv")
        if not p.endswith("_raw.csv"))
    omp_2d_frames = []
    for p in omp_2d_paths:
        df_p, _ = _load_csv(p, OMP_2D_COLS)
        if df_p is not None: omp_2d_frames.append(df_p)
    if omp_2d_frames:
        omp_2d_df = pd.concat(omp_2d_frames, ignore_index=True)
        omp_2d_df = omp_2d_df.drop_duplicates(
            subset=["omp_threads", "n_molecules"], keep="last")
        omp_2d_df = omp_2d_df.sort_values(["omp_threads", "n_molecules"]).reset_index(drop=True)
    else:
        omp_2d_df = None
    omp = dict(
        thread=_add_ci95(_load_csv(omp_thread_csv, OMP_THREAD_COLS)[0]),
        two_d= _add_ci95(omp_2d_df),
    )
    return data, omp


def _atoms(n_mol):
    return n_mol * ATOMS_PER_MOL


def _panel_letter(ax, letter, xoff=0.5, yoff=1.03, va="bottom"):
    """Place a bold panel letter (a, b, c, ...).  Default position is above
    the axes; pass yoff < 0 and va="top" to drop the letter below the x-axis
    label instead."""
    ax.text(xoff, yoff, f"({letter})", transform=ax.transAxes,
            fontsize=12, fontweight="bold", va=va, ha="center",
            color=CAMBRIDGE["slate_4"])


def _row_major(items, ncols):
    """Reorder items so matplotlib's column-major legend renders them
    in row-major (left-to-right, top-down) order at the given ncols."""
    n = len(items)
    nrows = (n + ncols - 1) // ncols
    out = []
    for col in range(ncols):
        for row in range(nrows):
            idx = row * ncols + col
            if idx < n:
                out.append(items[idx])
    return out


def _save(fig, name):
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, name + ".pdf")
    fig.savefig(path)
    plt.close(fig)
    print(f"  saved {path}")


def _log_decade_ticks(ax, axis, atom_fmt=False, expand_to_decades=False):
    """Force a major tick at every decade (works around LogLocator collapsing
    to one tick when the visible range spans only ~1.5 decades).
    expand_to_decades pushes limits out to enclosing decade boundaries."""
    target = ax.xaxis if axis == "x" else ax.yaxis
    lo_lim, hi_lim = ax.get_xlim() if axis == "x" else ax.get_ylim()
    lo = int(np.floor(np.log10(lo_lim)))
    hi = int(np.ceil(np.log10(hi_lim)))
    if expand_to_decades:
        new_lo, new_hi = 10.0 ** lo, 10.0 ** hi
        (ax.set_xlim if axis == "x" else ax.set_ylim)(new_lo, new_hi)
    ticks = [10.0 ** k for k in range(lo, hi + 1)]
    target.set_major_locator(ticker.FixedLocator(ticks))
    if atom_fmt:
        target.set_major_formatter(ticker.FuncFormatter(
            lambda v, _: f"{int(v):,}"))
    else:
        target.set_major_formatter(ticker.LogFormatterSciNotation(base=10))
    target.set_minor_locator(ticker.LogLocator(
        base=10, subs=np.arange(2, 10), numticks=100))
    target.set_minor_formatter(ticker.NullFormatter())


def fig1_algorithmic_complexity(data, omp):
    print("\n[fig 1] algorithmic complexity (linear + log size scaling + per-atom + speedup)")
    fig, axes = plt.subplots(2, 2, figsize=(W_TEXT * 1.45, 6.8))
    ax_a, ax_b = axes[0]
    ax_c, ax_d = axes[1]

    m = data["upstream master"]["size"]

    for name, d in data.items():
        s = d["size"]
        if s is None:
            continue
        ax_a.errorbar(_atoms(s["n_molecules"]), s["tps_mean"], yerr=s["tps_ci95"],
                      marker=d["marker"], color=d["color"], capsize=2, lw=1.6, ms=5.5)
    if m is not None:
        x_ref, t_ref = m["n_molecules"].iloc[0] * ATOMS_PER_MOL, m["tps_mean"].iloc[0]
        x_max = m["n_molecules"].iloc[-1] * ATOMS_PER_MOL
        x_smooth = np.linspace(x_ref, x_max, 50)
        ax_a.plot(x_smooth, t_ref * (x_smooth / x_ref),       "--",
                  color=CAMBRIDGE["slate_3"],  lw=1.0)
        ax_a.plot(x_smooth, t_ref * (x_smooth / x_ref) ** 2,  ":",
                  color=CAMBRIDGE["blue_dark"], lw=1.0)
        # Clip y to the measured master max so the N^2 parabola extends off-panel.
        master_max = m["tps_mean"].max()
        ax_a.set_ylim(0, master_max * 1.10)
    ax_a.set_xlabel(r"Number of atoms, $N$")
    ax_a.set_ylabel(r"Time per MD step (s)")
    _panel_letter(ax_a, "a")

    for name, d in data.items():
        s = d["size"]
        if s is None:
            continue
        ax_b.errorbar(_atoms(s["n_molecules"]), s["tps_mean"], yerr=s["tps_ci95"],
                      marker=d["marker"], color=d["color"], capsize=2, lw=1.6, ms=5.5)
    if m is not None:
        x_ref, t_ref = m["n_molecules"].iloc[0] * ATOMS_PER_MOL, m["tps_mean"].iloc[0]
        x_arr = np.array([x_ref, m["n_molecules"].iloc[-1] * ATOMS_PER_MOL])
        ax_b.plot(x_arr, t_ref * (x_arr / x_ref),       "--",
                  color=CAMBRIDGE["slate_3"],  lw=1.0, label=r"$\propto N$")
        ax_b.plot(x_arr, t_ref * (x_arr / x_ref) ** 2,  ":",
                  color=CAMBRIDGE["blue_dark"], lw=1.0, label=r"$\propto N^2$")
    ax_b.set(xscale="log", yscale="log")
    ax_b.set_xlabel(r"Number of atoms, $N$")
    ax_b.set_ylabel(r"Time per MD step (s)")
    _log_decade_ticks(ax_b, "x", atom_fmt=True)
    _panel_letter(ax_b, "b")

    for name, d in data.items():
        s = d["size"]
        if s is None:
            continue
        n_atoms = s["n_molecules"] * ATOMS_PER_MOL
        ax_c.errorbar(n_atoms, s["tps_mean"] / n_atoms,
                      yerr=s["tps_ci95"] / n_atoms,
                      marker=d["marker"], color=d["color"], capsize=2, lw=1.6, ms=5.5)
    ax_c.set(xscale="log", yscale="log")
    ax_c.set_xlabel(r"Number of atoms, $N$")
    ax_c.set_ylabel(r"Time per MD step per atom (s)")
    _log_decade_ticks(ax_c, "x", atom_fmt=True)
    _log_decade_ticks(ax_c, "y", expand_to_decades=True)
    _panel_letter(ax_c, "c")

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
            ax_d.plot(_atoms(merged["n_molecules"]), speedup,
                      marker=d["marker"], color=d["color"], lw=1.7, ms=6.5)
        baseline = ax_d.axhline(1.0, color=CAMBRIDGE["slate_3"], ls="--", lw=1.0)
    ax_d.set(xscale="log")
    ax_d.set_xlabel(r"Number of atoms, $N$")
    ax_d.set_ylabel("Speedup")
    if m is not None:
        ax_d.legend([baseline], [r"Master baseline ($1\times$)"],
                    loc="upper left", fontsize=8, frameon=True, framealpha=0.93)
    _log_decade_ticks(ax_d, "x", atom_fmt=True)
    _panel_letter(ax_d, "d")

    handles = []
    for name, d in data.items():
        if d["size"] is None:
            continue
        handles.append(plt.Line2D([0], [0], color=d["color"], marker=d["marker"],
                                  lw=1.6, ms=5.5, label=d["label"]))
    handles += [
        plt.Line2D([0], [0], color=CAMBRIDGE["slate_3"],  ls="--", lw=1.0,
                   label=r"$\propto N$"),
        plt.Line2D([0], [0], color=CAMBRIDGE["blue_dark"], ls=":",  lw=1.0,
                   label=r"$\propto N^{2}$"),
    ]
    # Two rows of three for legibility: 3 branches on top, 2 power-law
    # references on the bottom.  _row_major restores left-to-right reading
    # order against matplotlib's column-major legend layout.
    fig.legend(handles=_row_major(handles, 3), loc="upper center",
               bbox_to_anchor=(0.5, 1.00), ncol=3,
               frameon=False, fontsize=12, columnspacing=2.0,
               handlelength=2.6, handletextpad=0.6)

    fig.tight_layout(rect=[0, 0, 1, 0.91], h_pad=2.0, w_pad=2.5)
    _save(fig, "fig1_algorithmic_complexity")


def fig2_strong_scaling(data, omp):
    print("\n[fig 2] strong scaling (t/step + speedup + efficiency + overhead)")
    fig, axes = plt.subplots(2, 2, figsize=(W_TEXT * 1.45, 6.8))
    ax_a, ax_b = axes[0]
    ax_c, ax_d = axes[1]
    # 64 omitted: log_2(64)=6 collides with log_2(76)=6.25, and 76 is the
    # physical cap (full Peta4-IceLake node).
    N_TICKS = [1, 2, 4, 8, 16, 32, 76]

    # Common-baseline speedup: master's 1-core time divided by each branch's
    # N-core time. Avoids the per-branch-baseline artefact where the better-
    # optimised branch shows a lower speedup ratio just because its 1-core
    # baseline is already faster.
    master_t1 = None
    master_t1_ci95 = 0.0
    if data["upstream master"]["core"] is not None and not data["upstream master"]["core"].empty:
        master_t1 = data["upstream master"]["core"]["tps_mean"].iloc[0]
        master_t1_ci95 = data["upstream master"]["core"]["tps_ci95"].iloc[0]

    drawn_branches = []
    for name in ("upstream master", "native-spline", "native-spline-omp"):
        d = data[name]; c = d["core"]
        if c is None or c.empty:
            continue
        drawn_branches.append((name, d))
        ax_a.errorbar(c["total_cores"], c["tps_mean"], yerr=c["tps_ci95"],
                      marker=d["marker"], color=d["color"], capsize=2, lw=1.6, ms=5.5)
        # Index of the baseline (lowest core count) row; robust to row order.
        base_idx = int(np.asarray(c["total_cores"]).argmin())
        # Speedup uses a common (master 1-core) baseline.  For N>1 the
        # numerator and denominator are independent runs, so propagate CI95 of
        # both in quadrature (delta method).  master's own 1-core point is a
        # self-ratio (speedup == 1 exactly), so its interval is zero there.
        common_speedup = (master_t1 / c["tps_mean"]) if master_t1 is not None else c["speedup"]
        if master_t1 is not None and master_t1 > 0:
            rel_speedup = np.sqrt((master_t1_ci95/master_t1)**2 +
                                  (c["tps_ci95"]/c["tps_mean"])**2)
            speedup_ci = np.array(common_speedup*rel_speedup, dtype=float)
            if name == "upstream master":
                speedup_ci[base_idx] = 0.0   # master_t1 / master_t1 is exact
        else:
            speedup_ci = None
        ax_b.errorbar(c["total_cores"], common_speedup, yerr=speedup_ci,
                      marker=d["marker"], color=d["color"], capsize=2, lw=1.7, ms=6.5)
        # Efficiency uses the per-branch 1-core baseline.  At the baseline core
        # count efficiency is exactly 100% (self-ratio, zero variance); for
        # N>1 the 1-core and N-core runs are independent.
        branch_t1 = c["tps_mean"].iloc[base_idx]
        branch_t1_ci95 = c["tps_ci95"].iloc[base_idx]
        rel_eff = np.sqrt((branch_t1_ci95/branch_t1)**2 +
                          (c["tps_ci95"]/c["tps_mean"])**2)
        eff_ci = np.array(c["efficiency"]*rel_eff, dtype=float)
        eff_ci[base_idx] = 0.0           # baseline efficiency == 100% exactly
        ax_c.errorbar(c["total_cores"], c["efficiency"], yerr=eff_ci,
                      marker=d["marker"], color=d["color"], capsize=2, lw=1.7, ms=6.5)
        # Estimated overhead per step = t(N) - t(1)/N: absolute residual above
        # ideal strong-scaling (communication + serial + NUMA noise). Not
        # strictly t_comm but the dominant term at high N is collective MPI.
        # Skip the baseline point (trivially overhead=0 by construction); on
        # a log-y axis it would render as a spurious vertical drop.  This is a
        # difference of two independent runs, so the interval propagates as
        # sqrt(CI_tN^2 + (CI_t1/N)^2) (absolute, not relative).
        t1 = c["tps_mean"].iloc[base_idx]
        t1_ci95 = c["tps_ci95"].iloc[base_idx]
        c_tail = c.drop(index=c.index[base_idx])
        overhead = c_tail["tps_mean"] - t1 / c_tail["total_cores"]
        overhead_ci = np.sqrt(c_tail["tps_ci95"]**2 +
                              (t1_ci95 / c_tail["total_cores"])**2)
        ax_d.errorbar(c_tail["total_cores"], overhead.clip(lower=1e-4),
                      yerr=np.array(overhead_ci, dtype=float),
                      marker=d["marker"], color=d["color"], capsize=2, lw=1.7, ms=6.5)
    # Dense sampling: on log-x linear-y axes y=x is a curve; 7 points would
    # render as visible chord segments.
    x_ideal = np.linspace(1, max(N_TICKS), 200)
    ideal_lin, = ax_b.plot(x_ideal, x_ideal, "--", color=CAMBRIDGE["slate_3"], lw=1.0)
    ideal_eff  = ax_c.axhline(100, ls="--", color=CAMBRIDGE["slate_3"], lw=1.0)

    for ax in axes.flat:
        ax.set_xscale("log", base=2)
        ax.set_xticks(N_TICKS)
        ax.xaxis.set_major_formatter(ticker.FixedFormatter([str(c) for c in N_TICKS]))
        ax.set_xlim(0.85, 90)
    ax_a.set_yscale("log")
    ax_a.set_ylabel("Time per MD step (s)")
    _log_decade_ticks(ax_a, "y", expand_to_decades=True)
    ax_b.set_ylabel("Speedup")
    ax_c.set_ylabel("Parallel efficiency (%)")
    ax_c.set_ylim(0, 110)
    ax_d.set_yscale("log")
    ax_d.set_ylabel(r"Overhead per step (s)" "\n" r"($t_N - t_1/N$)")
    _log_decade_ticks(ax_d, "y", expand_to_decades=True)
    ax_b.legend([ideal_lin], ["Ideal (linear)"],
                loc="upper left", fontsize=10, frameon=True, framealpha=0.93)
    ax_c.legend([ideal_eff], ["Ideal (100%)"],
                loc="lower left", fontsize=10, frameon=True, framealpha=0.93)
    for ax, ltr in zip(axes.flat, "abcd"):
        _panel_letter(ax, ltr)

    fig.supxlabel("Total cores", fontsize=12, y=0.02)

    handles = [plt.Line2D([0], [0], color=d["color"], marker=d["marker"],
                          lw=1.9, ms=8.0, label=d["label"])
               for _, d in drawn_branches]
    fig.legend(handles=handles, loc="upper center",
               bbox_to_anchor=(0.5, 1.00), ncol=len(handles),
               frameon=False, fontsize=13, columnspacing=2.2, handlelength=3.0,
               handletextpad=0.7)

    fig.tight_layout(rect=[0, 0.03, 1, 0.96], h_pad=2.5, w_pad=2.5)
    _save(fig, "fig2_strong_scaling")


def fig3_openmp_threading(data, omp):
    print("\n[fig 3] OpenMP threading characterisation")
    fig, axes = plt.subplots(1, 2, figsize=(W_TEXT * 1.55, 4.2),
                             gridspec_kw=dict(wspace=0.55))
    ax_a, ax_b = axes
    omp_col = BRANCH_STYLE["native-spline-omp"]["color"]

    t = omp["thread"]
    panel_a_handles = []
    if t is not None and not t.empty:
        ticks = [int(v) for v in t["omp_threads"]]
        h_time = ax_a.errorbar(t["omp_threads"], t["tps_mean"], yerr=t["tps_ci95"],
                               marker="s", color=omp_col, capsize=2, lw=1.6, ms=5.5,
                               label="Time / step")
        ax_a.set_xscale("log", base=2); ax_a.set_yscale("log")
        ax_a.set_xticks(ticks); ax_a.xaxis.set_major_formatter(ticker.FixedFormatter([str(c) for c in ticks]))
        ax_a.set_xlabel("OMP threads")
        ax_a.set_ylabel("Time per MD step (s)")
        ax_a2 = ax_a.twinx()
        h_ideal,    = ax_a2.plot(t["omp_threads"], t["omp_threads"], "--",
                                 color=CAMBRIDGE["slate_3"], lw=1.0, label="Ideal (linear)")
        h_observed, = ax_a2.plot(t["omp_threads"], t["speedup"], marker="o",
                                 color=CAMBRIDGE["indigo"], lw=1.5, ms=5.5,
                                 label="Observed speedup")
        ax_a2.set_ylabel("Speedup vs OMP = 1", color=CAMBRIDGE["indigo"])
        ax_a2.tick_params(axis="y", labelcolor=CAMBRIDGE["indigo"])
        ax_a2.grid(False)
        panel_a_handles = [h_time, h_ideal, h_observed]
    _panel_letter(ax_a, "a", yoff=-0.22, va="top")

    panel_b_handles = []; panel_b_labels = []
    td = omp["two_d"]
    if td is not None and not td.empty:
        omps = sorted(td["omp_threads"].unique())
        # Cool-color gradient for OMP=1..8 (progressively darker, fits the
        # "more threads → darker" visual story).  OMP=16 breaks to crest_dark
        # since it sits in the diminishing-returns regime (≤10% gain over OMP=8)
        # -- the hue change deliberately flags this as a different operating point.
        ramp_colors = [CAMBRIDGE["blue"],       CAMBRIDGE["blue_warm"],
                       CAMBRIDGE["indigo"],     CAMBRIDGE["blue_dark"],
                       CAMBRIDGE["crest_dark"]]
        ramp = (ramp_colors * (len(omps) // len(ramp_colors) + 1))[:len(omps)]
        for col, ompv in zip(ramp, omps):
            sub = td[td["omp_threads"] == ompv].sort_values("n_molecules")
            x = _atoms(sub["n_molecules"])
            y = sub["tps_mean"]; ci = sub["tps_ci95"]
            # Shaded mean ± CI95 band (Student's-t, n=5) — variance now visible
            # at a glance instead of via thin error-bar caps.
            ax_b.fill_between(x, y - ci, y + ci, color=col, alpha=0.18, lw=0)
            ax_b.plot(x, y, marker="o", color=col, lw=1.5, ms=5)
            panel_b_handles.append(plt.Line2D([0], [0], color=col, marker="o",
                                              lw=1.5, ms=5, label=f"OMP = {int(ompv)}"))
            panel_b_labels.append(f"OMP = {int(ompv)}")
        ax_b.set_xscale("log", base=2); ax_b.set_yscale("log")
        ax_b.set_xlabel(r"Number of atoms, $N$")
        ax_b.set_ylabel("Time per MD step (s)")
        ax_b.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    _panel_letter(ax_b, "b", yoff=-0.22, va="top")

    if panel_a_handles:
        ax_a.legend(panel_a_handles,
                    [h.get_label() for h in panel_a_handles],
                    loc="lower center", bbox_to_anchor=(0.5, 1.04),
                    ncol=3, frameon=False, fontsize=11,
                    columnspacing=1.2, handlelength=2.0, handletextpad=0.5)
    if panel_b_handles:
        # Reorder so matplotlib's column-major legend renders row-major:
        # top row OMP=1,2,4; bottom row OMP=8,16.
        ncol_b = 3
        ax_b.legend(_row_major(panel_b_handles, ncol_b),
                    _row_major(panel_b_labels,  ncol_b),
                    loc="lower center", bbox_to_anchor=(0.5, 1.04),
                    ncol=ncol_b, frameon=False, fontsize=11,
                    columnspacing=1.4, handlelength=2.0, handletextpad=0.5)

    fig.tight_layout(rect=[0, 0.06, 1, 0.84])
    _save(fig, "fig3_openmp_threading")


HEADLINE_N = 1024
def fig4_statistical_significance(data, omp):
    print(f"\n[fig 4] statistical significance distribution at N = {HEADLINE_N}")
    samples = {}
    for name, d in data.items():
        if d.get("size_path") is None: continue
        raw_p = d["size_path"].replace(".csv", "_raw.csv")
        if not os.path.exists(raw_p): continue
        rdf = pd.read_csv(raw_p, comment="#", header=None,
                          names=["n_molecules", "rep", "time_per_step_s", "walltime_s"],
                          na_values=["NA", "FAILED"])
        sel = rdf[rdf["n_molecules"] == HEADLINE_N].dropna(subset=["time_per_step_s"])
        if len(sel):
            samples[name] = sel["time_per_step_s"].astype(float).to_numpy()
    if not samples:
        print("  no raw rep data found, skipping")
        return

    fig, ax = plt.subplots(figsize=(W_TEXT, 4.2))
    labels = list(samples.keys())
    pos    = np.arange(len(labels)) + 1
    colors = [data[k]["color"] for k in labels]

    # Violin/KDE dropped: at n=5 the KDE bandwidth dominates the curve shape
    # and isn't statistically meaningful; raw dots + box carry the information.
    bp = ax.boxplot([samples[k] for k in labels], positions=pos, widths=0.35,
                    patch_artist=True, showfliers=True,
                    boxprops=dict(linewidth=1.0),
                    medianprops=dict(color=CAMBRIDGE["slate_4"], linewidth=1.6),
                    whiskerprops=dict(linewidth=1.0),
                    capprops=dict(linewidth=1.0))
    for patch, col in zip(bp["boxes"], colors):
        patch.set_facecolor(col); patch.set_alpha(0.55); patch.set_edgecolor(CAMBRIDGE["slate_4"])

    rng = np.random.RandomState(42)
    for p, k, col in zip(pos, labels, colors):
        x = p + rng.uniform(-0.08, 0.08, size=len(samples[k]))
        ax.scatter(x, samples[k], s=34, color=col, edgecolor=CAMBRIDGE["slate_4"],
                   linewidth=0.7, alpha=0.95, zorder=3)

    ax.set_xticks(pos)
    ax.set_xticklabels([f"{l}\n({data[l].get('decomp', '')})".strip()
                        if False else f"{l}" for l in labels])
    ax.set_ylabel("Time per MD step (s)")
    ax.set_title(f"Headline-workload distribution at $N = {HEADLINE_N}$ H$_2$O "
                 f"({HEADLINE_N * ATOMS_PER_MOL} atoms)")
    ymax = max(np.max(samples[k]) for k in labels)
    ax.set_ylim(top=ymax * 1.08)

    fig.tight_layout()
    _save(fig, "fig4_statistical_significance")


# MAQAO data from the May-14 measurement (re-runs blocked by the post-2026-05-15
# CSD3 ptrace mitigation; see methods.md).
MAQAO_BRANCHES = ["master", "native-spline (pre-port)", "native-spline-omp"]
MAQAO_HOTSPOTS = {
    "master": {
        "MPI broadcast":                   33.25,
        "libm tanh / exp (NN activ.)":      8.91,
        "nnp_calc_ang":                     3.54,
        "pbc2 (O(N$^2$) neighbour cost)":   2.81,
        "libm pow (__powr8i4)":             2.55,
        "OMP runtime (atomics / barriers)": 0.00,
        "hermite_spline_value":             0.00,
        "other / user code":               48.94,
    },
    "native-spline (pre-port)": {
        "MPI broadcast":                   34.77,
        "libm tanh / exp (NN activ.)":      0.14,
        "nnp_calc_ang":                     2.64,
        "pbc2 (O(N$^2$) neighbour cost)":   0.00,
        "libm pow (__powr8i4)":             5.24,
        "OMP runtime (atomics / barriers)": 0.00,
        "hermite_spline_value":            10.88,
        "other / user code":               46.33,
    },
    "native-spline-omp": {
        "MPI broadcast":                   27.21,
        "libm tanh / exp (NN activ.)":      0.21,
        "nnp_calc_ang":                     1.65,
        "pbc2 (O(N$^2$) neighbour cost)":   0.00,
        "libm pow (__powr8i4)":             3.22,
        "OMP runtime (atomics / barriers)": 28.44,
        "hermite_spline_value":             5.98,
        "other / user code":               33.29,
    },
}
MAQAO_AAEFF = {"master": 52.92, "native-spline-omp": 78.17}


# ============================================================================
# Fresh CQA per-loop data (May-26, post-port aware).  Static disassembly only;
# numbers are deterministic per binary.  Used by figs 7-9.
# Produced via:
#   maqao cqa <binary> --fct-loops="nnp_" --output-format=csv \
#                      --output-path=<csv>
# ============================================================================
FRESH_CQA_DIR = "/home/crm98/cp2k-benchmarks/maqao_login_cqa"
FRESH_CQA_PATHS = {
    "master":            f"{FRESH_CQA_DIR}/master_nnp.csv",
    "native-spline":     f"{FRESH_CQA_DIR}/native_spline_postport_nnp.csv",
    "native-spline-omp": f"{FRESH_CQA_DIR}/omp_nnp.csv",
}
NNP_SRC_BASENAMES = {"nnp_acsf.F", "nnp_force.F", "nnp_model.F", "nnp_environment.F"}

def _load_fresh_cqa():
    """Read the three CQA CSVs, tag each row with branch, restrict to NNP
       Fortran sources.  Returns combined DataFrame or None if any CSV missing."""
    parts = []
    for branch, p in FRESH_CQA_PATHS.items():
        if not os.path.exists(p):
            print(f"  [_load_fresh_cqa] missing CSV: {p}")
            return None
        d = pd.read_csv(p, sep=";")
        d.insert(0, "branch", branch)
        parts.append(d)
    df = pd.concat(parts, ignore_index=True)
    df["src_short"] = df["Source file"].apply(
        lambda s: s.split("/")[-1] if isinstance(s, str) else "?")
    return df[df["src_short"].isin(NNP_SRC_BASENAMES)].copy()


def _branch_style_for(short):
    """Map a fresh-CSV branch key to the BRANCH_STYLE entry."""
    return BRANCH_STYLE[{"master": "upstream master",
                         "native-spline": "native-spline",
                         "native-spline-omp": "native-spline-omp"}[short]]


# ============================================================================
# ONE View HTML report parsers (May-14 measurement, pre-port).  Source pages
# come from `maqao oneview --create-report=one` runs.  Used by figs 7-10.
# ============================================================================
import re as _re
import html as _html
from pathlib import Path as _Path

ONE_VIEW_HTML_ROOT = "/home/crm98/cp2k-benchmarks/results/maqao"
ONE_VIEW_HTML_DIRS = {
    "master":            f"{ONE_VIEW_HTML_ROOT}/xp_master/cp2k.psmp_one_html",
    "native-spline":     f"{ONE_VIEW_HTML_ROOT}/xp_feature-nnp-native-spline/cp2k.psmp_one_html",
    "native-spline-omp": f"{ONE_VIEW_HTML_ROOT}/xp_feature-nnp-native-spline-omp/cp2k.psmp_one_html",
}

# Column order in the "Detailed Application Categorization" table after the
# row label.  Empirically observed; matches MAQAO ONE View 2026.0.0-b.
APP_CATEGORIES = ["Time(s)", "Binary", "MPI", "OMP", "TBB", "Math",
                  "System", "Pthread", "IO", "String", "Memory", "Others"]

# Column indices in the expert_summary table (excluding row label).
EXPERT_HDR_INDEX = {
    "Source Function":      3,
    "Coverage (% app. time)": 7,
    "Speedup if no scalar integer":     8,
    "Speedup if FP arith vectorized":   9,
    "Speedup if fully vectorized":     10,
    "Speedup if FP only":              11,
    "Vectorization Ratio (%)":         13,
    "CQA cycles":                      15,
}


def _ov_strip(s):
    return _html.unescape(_re.sub(r"<[^>]+>", " ", s).strip())


def _ov_cells(row_html):
    return [_ov_strip(c) for c in _re.findall(r"<t[hd][^>]*>([\s\S]*?)</t[hd]>",
                                              row_html)]


def parse_application_categorization(branch):
    """Read `application.html` for `branch` and return a dict with two keys:
        aggregate : {category -> value} for the whole run_0 row
        processes : list[(pid:str, dict{category -> value})] for each MPI rank
       Categories follow APP_CATEGORIES ordering.  Time(s) is absolute, all
       others are percentages of the per-row total."""
    path = _Path(ONE_VIEW_HTML_DIRS[branch]) / "application.html"
    if not path.exists(): return None
    rows = _re.findall(r"<tr[^>]*>([\s\S]*?)</tr>", path.read_text())
    out = {"aggregate": None, "processes": []}
    for r in rows:
        cs = _ov_cells(r)
        if not cs: continue
        nums_raw = cs[1:1 + len(APP_CATEGORIES)]
        if len(nums_raw) < len(APP_CATEGORIES): continue
        try:
            nums = [float(v.replace(",", "")) for v in nums_raw]
        except ValueError:
            continue
        d = dict(zip(APP_CATEGORIES, nums))
        label = cs[0]
        if "run_0" in label and out["aggregate"] is None:
            out["aggregate"] = d
        elif label.startswith("○ Thread") or "Thread" in label:
            out["processes"].append((label.replace("○", "").strip(), d))
    return out


def parse_expert_summary(branch):
    """Read `expert_summary.html` for `branch`.  Returns list of dicts, one
       per hot loop, with the EXPERT_HDR_INDEX-named keys parsed to float
       where possible (string passed through for the function-name column)."""
    path = _Path(ONE_VIEW_HTML_DIRS[branch]) / "expert_summary.html"
    if not path.exists(): return []
    rows = _re.findall(r"<tr[^>]*>([\s\S]*?)</tr>", path.read_text())
    loops = []
    for r in rows:
        cs = _ov_cells(r)
        if not cs or "Loop" not in cs[0]: continue
        if len(cs) < max(EXPERT_HDR_INDEX.values()) + 1: continue
        d = {}
        for key, idx in EXPERT_HDR_INDEX.items():
            v = cs[idx]
            if key == "Source Function":
                d[key] = v
                continue
            # speedup cols can be "1.06" or "8.68 - 7.23" (range across paths)
            v = v.split(" - ")[0].split(",")[0].strip()
            try:
                d[key] = float(v)
            except ValueError:
                d[key] = float("nan")
        loops.append(d)
    return loops


def _amdahl_speedup(coverage_pct, speedup):
    """Application-wide speedup factor if a set of hot loops with the given
       per-loop coverages (as % of app time) each speed up by speedup_i:
           S_app = 1 / [ (1 - sum c_i) + sum c_i / s_i ]
       coverage_pct / speedup may be array-like; NaNs are dropped."""
    c = np.asarray(coverage_pct, dtype=float) / 100.0
    s = np.asarray(speedup,      dtype=float)
    mask = np.isfinite(c) & np.isfinite(s) & (s > 0)
    c = c[mask]; s = s[mask]
    if len(c) == 0: return 1.0
    # ensure sum(c) <= 1 (defensive against double-counted outer loops)
    if c.sum() > 1.0: c = c / c.sum()
    denom = (1.0 - c.sum()) + (c / s).sum()
    return 1.0 / denom if denom > 0 else float("inf")


def fig7_per_process_load_balance(data, omp):
    """Three-panel figure: one panel per branch.  Each panel shows a stacked
       bar per MPI rank with the Binary/MPI/OMP/Math/Other share of the
       rank's wall-clock time.  Exposes load imbalance that the aggregate
       hotspot bar (fig 5) averages out."""
    print("\n[fig 7] Per-process load balance (ONE View 'application' page)")
    branches = ["master", "native-spline", "native-spline-omp"]
    parsed = {b: parse_application_categorization(b) for b in branches}
    if any(v is None for v in parsed.values()):
        print("  one or more application.html pages missing -- skipping"); return

    # Five visible categories; rest collapsed into "Other".
    SHOW = ["Binary", "MPI", "OMP", "Math"]
    fig, axes = plt.subplots(1, 3, figsize=(W_TEXT * 1.5, 4.4), sharey=True)
    colors = [CAMBRIDGE["slate_3"], CAMBRIDGE["crest"],
              CAMBRIDGE["indigo"],  CAMBRIDGE["cherry"], CAMBRIDGE["slate_2"]]

    for ax, branch in zip(axes, branches):
        procs = parsed[branch]["processes"]
        n = len(procs)
        if n == 0:
            ax.set_title(f"{branch} (no per-rank data)"); continue
        x = np.arange(n)
        bottoms = np.zeros(n)
        for cat, col in zip(SHOW + ["Other"], colors):
            if cat == "Other":
                vals = np.array([100.0 - sum(d[c] for c in SHOW)
                                 for _, d in procs])
            else:
                vals = np.array([d[cat] for _, d in procs])
            ax.bar(x, vals, bottom=bottoms, color=col, edgecolor=CAMBRIDGE["white"],
                   linewidth=0.4, label=cat)
            bottoms += vals
        ax.set_xticks(x)
        ax.set_xticklabels([str(i) for i in range(n)], fontsize=7)
        ax.set_xlabel(f"MPI rank (1...{n})", fontsize=9)
        ax.set_ylim(0, 102)
        ax.set_title(branch, fontsize=10)
        ax.grid(axis="y", ls="--", alpha=0.35)
        # MPI-share spread annotation
        mpi_vals = [d["MPI"] for _, d in procs]
        ax.text(0.02, 0.97,
                f"MPI share: {min(mpi_vals):.0f}--{max(mpi_vals):.0f}% "
                f"(spread {max(mpi_vals)-min(mpi_vals):.0f} pp)",
                transform=ax.transAxes, fontsize=7.5, va="top",
                color=CAMBRIDGE["slate_4"],
                bbox=dict(facecolor="white",
                          edgecolor=CAMBRIDGE["slate_2"],
                          alpha=0.85, pad=2))
    axes[0].set_ylabel("Share of rank's wall-clock time (%)")
    # Single legend on the right
    handles, labs = axes[0].get_legend_handles_labels()
    fig.legend(handles, labs, loc="center right", fontsize=9,
               bbox_to_anchor=(0.995, 0.5), frameon=True,
               title="Time category", title_fontsize=9)
    fig.suptitle("Per-MPI-rank load balance (LProf, N = 64 H$_2$O, 16 cores)",
                 fontsize=11.5, y=0.99)
    fig.tight_layout(rect=[0, 0, 0.90, 0.96])
    _save(fig, "fig7_per_process_load_balance")


def fig8_headroom_speedup(data, omp):
    """Coverage-weighted application-wide speedup if every hot loop achieved
       its CQA upper bound for each optimisation class.  Four classes:
         - no scalar integer
         - FP arithmetic vectorized
         - fully vectorized
         - FP-only (CQA-projected scenario)
       Shown as grouped bars: one cluster per branch, four bars per cluster."""
    print("\n[fig 8] Coverage-weighted Amdahl headroom (ONE View 'expert_summary')")
    branches = ["master", "native-spline", "native-spline-omp"]
    loops_by_branch = {b: parse_expert_summary(b) for b in branches}
    if not any(loops_by_branch.values()):
        print("  no expert_summary data -- skipping"); return

    OPTS = [
        ("Speedup if no scalar integer",   "No scalar-int penalty"),
        ("Speedup if FP arith vectorized", "FP arith vectorised"),
        ("Speedup if fully vectorized",    "Fully vectorised"),
        ("Speedup if FP only",             "FP-only kernel"),
    ]

    rows = []
    for b in branches:
        loops = loops_by_branch[b]
        cov = [L["Coverage (% app. time)"] for L in loops]
        d = {"branch": b}
        for k, _ in OPTS:
            sp = [L[k] for L in loops]
            d[k] = _amdahl_speedup(cov, sp)
        rows.append(d)
    df = pd.DataFrame(rows).set_index("branch")

    fig, ax = plt.subplots(figsize=(W_TEXT, 4.4))
    x = np.arange(len(OPTS))
    width = 0.26
    colors = [CAMBRIDGE["slate_3"], CAMBRIDGE["blue_warm"], CAMBRIDGE["crest"]]
    for i, b in enumerate(branches):
        offsets = (i - 1) * width
        vals = [df.at[b, k] for k, _ in OPTS]
        bars = ax.bar(x + offsets, vals, width, color=colors[i],
                      edgecolor=CAMBRIDGE["slate_4"], linewidth=0.6,
                      label=_branch_style_for(b)["label"])
        for xi, v, bar in zip(x, vals, bars):
            ax.text(bar.get_x() + bar.get_width()/2, v + 0.02,
                    f"{v:.2f}x", ha="center", va="bottom",
                    fontsize=8, color=CAMBRIDGE["slate_4"])
    ax.set_xticks(x)
    ax.set_xticklabels([lbl for _, lbl in OPTS], fontsize=9, rotation=10)
    ax.axhline(1.0, color=CAMBRIDGE["slate_3"], ls=":", lw=0.8)
    ax.text(0.01, 1.02, "current performance (1.00x)", transform=ax.get_yaxis_transform(),
            fontsize=8, color=CAMBRIDGE["slate_3"], style="italic")
    ax.set_ylabel("Application-wide speedup factor (CQA upper bound, Amdahl-weighted)")
    ax.set_title("Headroom: theoretical app speedup if every hot loop hit each CQA bound")
    ax.legend(loc="upper left", fontsize=8, frameon=True)
    ax.grid(axis="y", ls="--", alpha=0.4)
    ax.set_ylim(0.9, max(df.values.max() * 1.15, 1.4))
    fig.tight_layout()
    _save(fig, "fig8_headroom_speedup")


def fig9_library_category_share(data, omp):
    """Aggregate library-category time share per branch -- MAQAO's official
       categorisation rather than the hand-curated HOTSPOTS dict.  Stacked
       horizontal bar so the long category names sit cleanly outside."""
    print("\n[fig 9] Aggregate library-category time share (ONE View 'application')")
    branches = ["master", "native-spline", "native-spline-omp"]
    parsed = {b: parse_application_categorization(b) for b in branches}
    if any(v is None for v in parsed.values()):
        print("  one or more application.html pages missing -- skipping"); return

    # Collapse the noise categories into "Other".
    KEEP = ["Binary", "MPI", "OMP", "Math", "Memory"]
    OTHER_KEYS = [k for k in APP_CATEGORIES if k not in KEEP and k != "Time(s)"]

    fig, ax = plt.subplots(figsize=(W_TEXT, 3.6))
    y = np.arange(len(branches))
    bottoms = np.zeros(len(branches))
    colors = [CAMBRIDGE["slate_3"], CAMBRIDGE["crest"],
              CAMBRIDGE["indigo"],  CAMBRIDGE["cherry"],
              CAMBRIDGE["green"],   CAMBRIDGE["slate_2"]]
    for cat, col in zip(KEEP + ["Other"], colors):
        if cat == "Other":
            vals = np.array([sum(parsed[b]["aggregate"][k] for k in OTHER_KEYS)
                             for b in branches])
        else:
            vals = np.array([parsed[b]["aggregate"][cat] for b in branches])
        ax.barh(y, vals, left=bottoms, color=col, edgecolor=CAMBRIDGE["white"],
                linewidth=0.6, label=cat, height=0.62)
        for yi, vi, bi in zip(y, vals, bottoms):
            if vi >= 3:
                ax.text(bi + vi/2, yi, f"{vi:.0f}%",
                        ha="center", va="center", fontsize=8, color="white",
                        fontweight="bold")
        bottoms += vals
    ax.set_yticks(y)
    ax.set_yticklabels([_branch_style_for(b)["label"] for b in branches], fontsize=9)
    ax.set_xlim(0, 102)
    ax.set_xlabel("Share of application wall-clock time (%)")
    ax.set_title("Library-category time share (MAQAO ONE View 'Detailed Application Categorization')")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9,
              frameon=True, title="Category", title_fontsize=9)
    ax.grid(axis="x", ls="--", alpha=0.35)
    fig.tight_layout()
    _save(fig, "fig9_library_category_share")


def fig10_top_loops_by_gain(data, omp):
    """Per-loop *opportunity* ranking: for each hot loop, the application-wide
       gain you'd unlock by fully vectorising IT ALONE.  Computed via Amdahl
       on each single loop, so it weights both 'how much time does this loop
       cost?' and 'how much faster could it be?'.  Top 10 per branch in
       three panels."""
    print("\n[fig 10] Top NNP loops ranked by potential vectorisation gain")
    branches = ["master", "native-spline", "native-spline-omp"]
    loops_by_branch = {b: parse_expert_summary(b) for b in branches}
    if not any(loops_by_branch.values()):
        print("  no expert_summary data -- skipping"); return

    fig, axes = plt.subplots(1, 3, figsize=(W_TEXT * 1.55, 4.6))
    for ax, b in zip(axes, branches):
        loops = loops_by_branch[b]
        gains = []
        for L in loops:
            c = L["Coverage (% app. time)"]; s = L["Speedup if fully vectorized"]
            if not (np.isfinite(c) and np.isfinite(s) and s > 1.0 and c > 0):
                continue
            # single-loop Amdahl: 1 / [(1-c/100) + (c/100)/s] -- the speedup
            # you'd get if ONLY this loop were optimised to its CQA bound
            s_app = 1.0 / ((1.0 - c/100) + (c/100)/s)
            gain_pct = (s_app - 1.0) * 100.0
            fn = L["Source Function"][:24]
            gains.append((gain_pct, fn, c, s, L["Vectorization Ratio (%)"]))
        gains.sort(reverse=True)
        gains = gains[:10]
        if not gains:
            ax.set_title(f"{b} (no loops with measurable gain)"); continue

        ys = np.arange(len(gains))[::-1]
        vals = [g[0] for g in gains]
        labs = [f"{g[1]}  ({g[2]:.1f}% cov, {g[3]:.1f}x bound)" for g in gains]
        ax.barh(ys, vals, color=_branch_style_for(b)["color"],
                edgecolor=CAMBRIDGE["slate_4"], linewidth=0.5)
        for yi, v in zip(ys, vals):
            ax.text(v + 0.05, yi, f"+{v:.1f}%", va="center",
                    fontsize=7.5, color=CAMBRIDGE["slate_4"])
        ax.set_yticks(ys)
        ax.set_yticklabels(labs, fontsize=7)
        ax.set_xlabel("App-wide speedup gain (%)")
        ax.set_title(_branch_style_for(b)["label"], fontsize=9.5)
        ax.grid(axis="x", ls="--", alpha=0.35)
        ax.set_xlim(0, max(vals) * 1.4)
    fig.suptitle("Top-10 NNP loops by single-loop full-vectorisation opportunity",
                 fontsize=11.5, y=1.01)
    fig.tight_layout()
    _save(fig, "fig10_top_loops_by_gain")


def _unused_fig_old_stub(data, omp):
    """(deprecated; kept as placeholder so editor diffs are minimal)"""
    fig, ax = plt.subplots(figsize=(W_TEXT, 4.4))

    THRESH = [50.0, 80.0, 95.0]
    branches = ["master", "native-spline", "native-spline-omp"]
    legend_lines = []
    max_rank = 0

    for branch in branches:
        style = _branch_style_for(branch)
        sub = df[df["branch"] == branch]
        cyc = (sub["[L1] Nb cycles: min"].astype(float)
                  .dropna().sort_values(ascending=False).values)
        if len(cyc) == 0: continue
        cum = np.cumsum(cyc) / cyc.sum() * 100.0
        rank = np.arange(1, len(cyc) + 1)
        max_rank = max(max_rank, len(cyc))

        # ranks at which the curve first crosses each threshold
        k = [int(np.searchsorted(cum, t) + 1) for t in THRESH]
        ax.plot(rank, cum, color=style["color"], lw=1.8)
        ax.scatter(k, THRESH, s=42, marker=style["marker"],
                   facecolor=style["color"],
                   edgecolor=CAMBRIDGE["slate_4"], linewidth=0.7, zorder=5)

        # legend entry includes the per-threshold rank counts
        n_total = len(cyc)
        legend_lines.append(
            (style["label"], style["color"], style["marker"],
             k[0], k[1], k[2], n_total))

    # Reference dotted lines at the threshold percentages
    for t in THRESH:
        ax.axhline(t, color=CAMBRIDGE["slate_2"], ls=":", lw=0.7, zorder=1)
        ax.text(max_rank * 0.995, t + 0.6, f"{int(t)}%",
                fontsize=7.5, ha="right", va="bottom",
                color=CAMBRIDGE["slate_3"], style="italic")

    # Custom legend with the threshold-rank info baked in
    from matplotlib.lines import Line2D
    handles = []
    labels = []
    for lbl, col, mk, k50, k80, k95, n_total in legend_lines:
        handles.append(Line2D([0], [0], color=col, marker=mk, lw=1.8,
                               markersize=7, markeredgecolor=CAMBRIDGE["slate_4"],
                               markeredgewidth=0.7))
        labels.append(f"{lbl}\n  50% at rank {k50}, 80% at {k80}, "
                      f"95% at {k95}  (of {n_total} loops)")
    leg = ax.legend(handles, labels, loc="lower right", fontsize=8,
                    frameon=True, handlelength=2.4, labelspacing=0.9,
                    title="Branch  —  rank at which top-N loops account for ...",
                    title_fontsize=8.5)
    leg.get_title().set_color(CAMBRIDGE["slate_4"])

    ax.set_xlabel("Loop rank  (k = top-$k$ slowest NNP loops, sorted by cycles/iter)")
    ax.set_ylabel("Cumulative share of static cycle budget (%)")
    ax.set_title("Where the static cycle budget lives across NNP loops")
    ax.set_xlim(0, max_rank); ax.set_ylim(0, 102)
    ax.grid(ls="--", alpha=0.3)
    fig.tight_layout()
    _save(fig, "fig7_work_concentration")

def fig5_maqao_microarchitectural(data, omp):
    print("\n[fig 5] MAQAO microarchitectural analysis")
    fig, axes = plt.subplots(1, 2, figsize=(W_TEXT * 1.18, 4.5),
                             gridspec_kw={"width_ratios": [1.7, 1.0]})
    ax_a, ax_b = axes

    categories = list(next(iter(MAQAO_HOTSPOTS.values())).keys())
    cat_colors = [
        CAMBRIDGE["blue_warm"],
        CAMBRIDGE["cherry"],
        CAMBRIDGE["green"],
        CAMBRIDGE["crest"],
        CAMBRIDGE["purple"],
        CAMBRIDGE["crest_dark"],
        CAMBRIDGE["indigo"],
        CAMBRIDGE["slate_2"],
    ]
    x = np.arange(len(MAQAO_BRANCHES))
    bottoms = np.zeros(len(MAQAO_BRANCHES))
    for cat, col in zip(categories, cat_colors):
        vals = np.array([MAQAO_HOTSPOTS[b][cat] for b in MAQAO_BRANCHES])
        ax_a.bar(x, vals, bottom=bottoms, label=cat, color=col,
                 edgecolor=CAMBRIDGE["white"], linewidth=0.6)
        for xi, vi, bi in zip(x, vals, bottoms):
            if vi >= 4:
                ax_a.text(xi, bi + vi/2, f"{vi:.0f}%", ha="center", va="center",
                          fontsize=8, color="white", fontweight="bold")
        bottoms += vals
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(MAQAO_BRANCHES, fontsize=9, rotation=25, ha="right")
    ax_a.set_ylim(0, 105)
    ax_a.set_ylabel("Share of application time (%)")
    ax_a.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8,
                title="Leaf-function category", title_fontsize=9, frameon=True)
    ax_a.grid(axis="y", ls="--", alpha=0.4)
    _panel_letter(ax_a, "a")

    labels = ["master", "native-spline-omp"]
    vals   = [MAQAO_AAEFF["master"], MAQAO_AAEFF["native-spline-omp"]]
    bcols  = [CAMBRIDGE["slate_3"], CAMBRIDGE["blue_warm"]]
    bars = ax_b.bar(np.arange(2), vals, color=bcols, width=0.55,
                    edgecolor=CAMBRIDGE["slate_4"], linewidth=0.8)
    for bar, v in zip(bars, vals):
        ax_b.text(bar.get_x() + bar.get_width()/2, v + 1.5, f"{v:.1f}%",
                  ha="center", va="bottom", fontsize=9, fontweight="bold",
                  color=CAMBRIDGE["slate_4"])
    ax_b.axhline(100, color=CAMBRIDGE["slate_3"], ls=":", lw=0.8)
    ax_b.text(0.5, 99, "Theoretical max", color=CAMBRIDGE["slate_3"],
              ha="center", va="top", fontsize=8, style="italic")
    ax_b.set_xticks(np.arange(2))
    ax_b.set_xticklabels(labels, fontsize=9, rotation=25, ha="right")
    ax_b.set_ylim(0, 110)
    ax_b.set_ylabel("Array-access efficiency (%)")
    _panel_letter(ax_b, "b")

    fig.tight_layout()
    _save(fig, "fig5_maqao_microarchitectural")


# A true MPI x OMP hybrid sweep at fixed N was not run; this renders the
# omp_size_scaling grid (native-spline-omp, MPI=1, OMP in {1,2,4,8,16,32,76},
# N in {64,256,512,1024,2048}) as a heatmap for picking the best OMP setting
# per workload.
def fig6_omp_size_heatmap(data, omp):
    print("\n[fig 6] OMP x N hybrid sweep heatmap (preview)")
    td = omp["two_d"]
    if td is None or td.empty:
        print("  no two_d omp data available, skipping")
        return

    n_atoms_col = td["n_molecules"] * ATOMS_PER_MOL
    pivot = (td.assign(n_atoms=n_atoms_col)
               .pivot(index="omp_threads", columns="n_atoms", values="tps_mean")
               .sort_index(axis=0)
               .sort_index(axis=1))

    fig, ax = plt.subplots(figsize=(W_TEXT * 1.05, W_TEXT * 0.85))

    # Cambridge-palette sequential cmap: blue_light (fast) -> blue_warm ->
    # blue_dark (slow), so good values appear pale and bad values saturated.
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list(
        "cambridge_seq",
        [CAMBRIDGE["blue_light"], CAMBRIDGE["blue_warm"], CAMBRIDGE["blue_dark"]],
        N=256,
    )

    im = ax.imshow(pivot.values, aspect="auto", cmap=cmap,
                   origin="lower", interpolation="nearest")

    best_val = pivot.values[~np.isnan(pivot.values)].min()
    for i, omp_v in enumerate(pivot.index):
        for j, n_v in enumerate(pivot.columns):
            v = pivot.values[i, j]
            if np.isnan(v):
                continue
            rel = (v - pivot.values[~np.isnan(pivot.values)].min()) / \
                  max(1e-9, np.nanmax(pivot.values) - np.nanmin(pivot.values))
            color = "white" if rel > 0.55 else CAMBRIDGE["slate_4"]
            weight = "bold" if v == best_val else "normal"
            ax.text(j, i, f"{v:.3f}", ha="center", va="center",
                    color=color, fontsize=9, fontweight=weight)
            if v == best_val:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                           fill=False, edgecolor=CAMBRIDGE["crest_dark"],
                                           linewidth=2.5))

    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels([f"{int(c):,}" for c in pivot.columns])
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels([int(t) for t in pivot.index])
    ax.set_xlabel("Number of atoms, $N$")
    ax.set_ylabel("OMP threads (MPI = 1)")
    cbar = fig.colorbar(im, ax=ax, label="Time per MD step (s)", pad=0.02)
    cbar.outline.set_edgecolor(CAMBRIDGE["slate_4"])

    best_omp = pivot.values.argmin() // pivot.shape[1]
    best_n_idx = pivot.values.argmin() % pivot.shape[1]
    ax.set_title(f"Native-spline-omp hybrid sweep — best ({pivot.index[best_omp]} OMP, "
                 f"N={int(pivot.columns[best_n_idx]):,}) = {best_val:.3f} s/step",
                 loc="center", pad=8, fontsize=11)

    fig.tight_layout()
    _save(fig, "fig6_omp_size_heatmap")


# ============================================================================
# Perf-stat throughput + DRAM summary (login-node perf run, no SLURM burn)
# Reads results/perf_cache/<TIMESTAMP>/{master,feature-nnp-native-spline}/
# perf_loads.txt — see scripts/CSD3_benchmark_scripts/perf_cache/.
# ============================================================================
PERF_RESULTS_ROOT = "/home/crm98/cp2k-benchmarks/results/perf_cache"
PERF_BRANCHES = [
    ("master",                    "master",                    CAMBRIDGE["slate_3"]),
    ("feature-nnp-native-spline", "feature/nnp-native-spline", CAMBRIDGE["blue_warm"]),
]
_PERF_EVENT_LINE   = re.compile(r"^\s*([\d,]+)\s+([A-Za-z0-9_\.\-]+):u\b")
_PERF_ELAPSED_LINE = re.compile(r"^\s*([0-9.]+)\s+seconds time elapsed")


def _parse_perf_stat(path):
    """Parse one `perf stat -o` file into {event: int_count, 'wall_s': float}."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            m = _PERF_EVENT_LINE.match(line)
            if m:
                out[m.group(2)] = int(m.group(1).replace(",", ""))
                continue
            m = _PERF_ELAPSED_LINE.match(line)
            if m:
                out["wall_s"] = float(m.group(1))
    return out


def _load_perf_cache():
    """Locate the most recent results/perf_cache/<TIMESTAMP>/ and merge each
    branch's loads + backing perf-stat passes.  Returns (ts, {branch: dict})
    or None if no results are present."""
    if not os.path.isdir(PERF_RESULTS_ROOT):
        return None
    ts_dirs = sorted(d for d in os.listdir(PERF_RESULTS_ROOT)
                     if os.path.isdir(os.path.join(PERF_RESULTS_ROOT, d)))
    if not ts_dirs:
        return None
    ts  = ts_dirs[-1]
    out = {}
    for key, _, _ in PERF_BRANCHES:
        bdir    = os.path.join(PERF_RESULTS_ROOT, ts, key)
        loads   = _parse_perf_stat(os.path.join(bdir, "perf_loads.txt"))
        backing = _parse_perf_stat(os.path.join(bdir, "perf_backing.txt"))
        merged  = dict(backing)
        merged.update(loads)
        merged["wall_s"] = loads.get("wall_s", backing.get("wall_s"))
        out[key] = merged
    return ts, out


def fig7_perf_summary(data, omp):
    """Throughput + memory-pressure summary from perf-stat hardware counters.

    Four panels in a single row:
      (a) wall time          — lower better
      (b) retired instructions — lower better
      (c) IPC                — higher better
      (d) DRAM-resident load traffic (load-retired L3 misses) — lower better

    All four favour the native-spline branch, framed as the compute-for-memory
    trade described in §4: the L1/L2 spike from spline-table lookups is fully
    absorbed inside the 1.25 MiB L2, and DRAM traffic itself drops.
    """
    print("\n[fig 7] perf-stat throughput + DRAM summary")
    loaded = _load_perf_cache()
    if loaded is None:
        print("  no perf_cache results found, skipping")
        return
    ts, perf = loaded
    if not all("instructions" in perf[k] for k, _, _ in PERF_BRANCHES):
        print(f"  perf_cache/{ts} is missing events for one or more branches, skipping")
        return
    print(f"  reading perf_cache/{ts}")

    fig, axes = plt.subplots(1, 4, figsize=(W_TEXT * 1.30, 3.2),
                             gridspec_kw={"wspace": 0.60})

    # (panel letter, y-axis label, value extractor, format-string)
    metrics = [
        ("a", "Wall time (s)",
            lambda d: d["wall_s"],                            "{:.2f}"),
        ("b", "Retired instructions ($\\times 10^{9}$)",
            lambda d: d["instructions"] / 1.0e9,              "{:.1f}"),
        ("c", "Instructions per cycle",
            lambda d: d["instructions"] / d["cycles"],        "{:.2f}"),
        ("d", "Load-retired L3 misses ($\\times 10^{3}$)",
            lambda d: d["mem_load_retired.l3_miss"] / 1.0e3,  "{:.0f}"),
    ]

    labels  = [lab for _, lab, _ in PERF_BRANCHES]
    colours = [c   for _, _,   c in PERF_BRANCHES]

    for ax, (letter, ylabel, getter, fmt) in zip(axes, metrics):
        values = [getter(perf[k]) for k, _, _ in PERF_BRANCHES]

        bars = ax.bar(np.arange(2), values, color=colours, width=0.55,
                      edgecolor=CAMBRIDGE["slate_4"], linewidth=0.8)
        for bar, v in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width()/2, v * 1.02,
                    fmt.format(v),
                    ha="center", va="bottom",
                    fontsize=9, fontweight="bold",
                    color=CAMBRIDGE["slate_4"])

        ax.set_xticks(np.arange(2))
        ax.set_xticklabels(labels, fontsize=8.5, rotation=25, ha="right")
        ax.set_ylim(0, max(values) * 1.18)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(f"({letter})", fontsize=12, fontweight="bold", pad=6)
        ax.grid(axis="y", ls="--", alpha=0.4)

    fig.tight_layout()
    _save(fig, "fig7_perf_summary")


# ============================================================================
# Compiler-flags appendix table
# ============================================================================
# Compiler flags passed by the build scripts (cp2k_CSD3_master_build.sh and
# cp2k_CSD3_opt_build.sh).  Both scripts pass the same set; CMake's Release
# preset and CP2K's CMakeLists.txt add further flags downstream, but those
# are build-system plumbing rather than choices made for this work.
COMPILER_FLAGS = [
    ("-O2",
        "C, C++, Fortran",
        "Base optimisation level: aggressive inlining, scalar replacement "
        "of aggregates, and inner-loop transforms enabled."),
    ("-g",
        "C, C++, Fortran",
        "Emit DWARF debug symbols so MAQAO and perf can resolve samples "
        "back to source lines."),
    ("-xCORE-AVX512",
        "C, C++, Fortran",
        "Target the Intel Ice Lake-server AVX-512 ISA; implies all earlier "
        "vector extensions (SSE 4.2, AVX, AVX2)."),
    ("-qopenmp",
        "C, C++, Fortran",
        "Enable OpenMP directives.  Required for the OMP-threaded NN "
        "evaluation in the feature/nnp-native-spline-omp branch and for "
        "the threaded MKL paths used by CP2K's solver."),
    ("-funroll-loops",
        "Fortran",
        "Unroll loops with statically known trip counts; complements the "
        "vectoriser on the symmetry-function kernels."),
    ("-ftree-vectorize",
        "Fortran",
        "Run the inner-loop vectoriser independently of the optimisation "
        "level."),
]


def _wrap_flag_token(text, width):
    """Wrap a flag-style string at width, breaking at `_`, `=`, `,`, `/`
    in addition to spaces/hyphens that textwrap already understands."""
    if len(text) <= width:
        return text
    breakables = "_=,/ "
    tokens, current = [], ""
    for ch in text:
        current += ch
        if ch in breakables:
            tokens.append(current)
            current = ""
    if current:
        tokens.append(current)
    lines, line = [], ""
    for t in tokens:
        if line and len(line) + len(t) > width:
            lines.append(line.rstrip())
            line = t
        else:
            line += t
    if line:
        lines.append(line.rstrip())
    return "\n".join(lines)


def fig8_compiler_flags(data, omp):
    """Compiler-flags appendix table rendered with the same palette as the
    rest of the thesis figures.  Reads its data from COMPILER_FLAGS above
    (no external inputs)."""
    print("\n[fig 8] compiler-flags appendix table")

    HEADERS    = ("Flag", "Languages", "Purpose / notes")
    COL_WIDTHS = [0.22, 0.20, 0.58]
    WRAP_CHARS = (22, 16, 45)

    n_rows = len(COMPILER_FLAGS)
    n_cols = len(HEADERS)

    cell_text = [
        [
            _wrap_flag_token(flag, WRAP_CHARS[0]),
            "\n".join(textwrap.wrap(langs,   WRAP_CHARS[1])) or langs,
            "\n".join(textwrap.wrap(purpose, WRAP_CHARS[2])) or purpose,
        ]
        for flag, langs, purpose in COMPILER_FLAGS
    ]
    row_lines = [max(c.count("\n") + 1 for c in row) for row in cell_text]

    line_inches   = 0.18
    header_inches = 0.40
    total_h = header_inches + sum(r * line_inches for r in row_lines) + 0.40

    fig, ax = plt.subplots(figsize=(8.5, total_h))
    ax.axis("off")

    table = ax.table(cellText=cell_text,
                     colLabels=list(HEADERS),
                     cellLoc="left",
                     loc="upper center",
                     colWidths=COL_WIDTHS)
    table.auto_set_font_size(False)
    table.set_fontsize(9)

    base_h = 1.0 / (n_rows + 1)
    for col in range(n_cols):
        table[(0, col)].set_height(base_h)
    for i, lines in enumerate(row_lines, start=1):
        h = base_h * (lines * 0.80 + 0.4)
        for col in range(n_cols):
            table[(i, col)].set_height(h)

    # Header row: navy fill, white bold text.
    for col in range(n_cols):
        cell = table[(0, col)]
        cell.set_facecolor(CAMBRIDGE["blue_dark"])
        cell.set_text_props(weight="bold", color=CAMBRIDGE["white"],
                            size=10, ha="left")
        cell.set_edgecolor("white")
        cell.set_linewidth(0.5)
        cell.PAD = 0.04

    # Body rows: alternating slate/white stripes.
    body_stripe = [CAMBRIDGE["slate_1"], CAMBRIDGE["white"]]
    for i in range(1, n_rows + 1):
        for col in range(n_cols):
            cell = table[(i, col)]
            cell.set_facecolor(body_stripe[(i - 1) % 2])
            cell.set_text_props(color=CAMBRIDGE["slate_4"], ha="left",
                                va="center")
            cell.set_edgecolor(CAMBRIDGE["slate_2"])
            cell.set_linewidth(0.4)
            cell.PAD = 0.04

    _save(fig, "fig8_compiler_flags")


FIGURES = [
    (1, fig1_algorithmic_complexity),
    (2, fig2_strong_scaling),
    (3, fig3_openmp_threading),
    (4, fig4_statistical_significance),
    (5, fig5_maqao_microarchitectural),
    (6, fig6_omp_size_heatmap),
    (7, fig7_perf_summary),
    (8, fig8_compiler_flags),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", type=int, choices=[1, 2, 3, 4, 5, 6, 7, 8], default=None,
                    help="only produce figure N")
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"Loading data from {RESULTS_ROOT}")
    data, omp = load_all_data()
    for n in BRANCH_DIRS:
        s = data[n].get("size")
        rows = len(s) if s is not None else 0
        print(f"  {n:22s}  size_rows={rows}")

    for num, fn in FIGURES:
        if args.only and args.only != num:
            continue
        fn(data, omp)

    print(f"\nDone.  PDFs in {OUT_DIR}")


if __name__ == "__main__":
    main()
