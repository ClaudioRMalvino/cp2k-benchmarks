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
    "font.size":          11,
    "axes.labelsize":     11,
    "axes.titlesize":     11.5,
    "xtick.labelsize":    9.5,
    "ytick.labelsize":    9.5,
    "legend.fontsize":    9,
    "figure.titlesize":   12.5,
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


def _panel_letter(ax, letter, xoff=0.5, yoff=1.03):
    """Place a bold panel letter (a, b, c, ...) above the axes."""
    ax.text(xoff, yoff, f"({letter})", transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="bottom", ha="center",
            color=CAMBRIDGE["slate_4"])


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
    ax_d.set_ylabel(r"Speedup ($t_\mathrm{master} / t_\mathrm{branch}$)")
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
    fig.legend(handles=handles, loc="upper center",
               bbox_to_anchor=(0.5, 1.00), ncol=len(handles),
               frameon=False, fontsize=10, columnspacing=1.8, handlelength=2.6)

    fig.tight_layout(rect=[0, 0, 1, 0.96], h_pad=2.0, w_pad=2.5)
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
    if data["upstream master"]["core"] is not None and not data["upstream master"]["core"].empty:
        master_t1 = data["upstream master"]["core"]["tps_mean"].iloc[0]

    drawn_branches = []
    for name in ("upstream master", "native-spline", "native-spline-omp"):
        d = data[name]; c = d["core"]
        if c is None or c.empty:
            continue
        drawn_branches.append((name, d))
        ax_a.errorbar(c["total_cores"], c["tps_mean"], yerr=c["tps_ci95"],
                      marker=d["marker"], color=d["color"], capsize=2, lw=1.6, ms=5.5)
        common_speedup = (master_t1 / c["tps_mean"]) if master_t1 is not None else c["speedup"]
        ax_b.plot(c["total_cores"], common_speedup,
                  marker=d["marker"], color=d["color"], lw=1.7, ms=6.5)
        ax_c.plot(c["total_cores"], c["efficiency"],
                  marker=d["marker"], color=d["color"], lw=1.7, ms=6.5)
        # Estimated overhead per step = t(N) - t(1)/N: absolute residual above
        # ideal strong-scaling (communication + serial + NUMA noise). Not
        # strictly t_comm but the dominant term at high N is collective MPI.
        # Skip the baseline point (trivially overhead=0 by construction); on
        # a log-y axis it would render as a spurious vertical drop.
        t1 = c["tps_mean"].iloc[0]
        c_tail = c.iloc[1:]
        overhead = c_tail["tps_mean"] - t1 / c_tail["total_cores"]
        ax_d.plot(c_tail["total_cores"], overhead.clip(lower=1e-4),
                  marker=d["marker"], color=d["color"], lw=1.7, ms=6.5)
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
    ax_b.set_ylabel(r"Speedup ($t_\mathrm{master,1} / t_\mathrm{branch,N}$)")
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
               frameon=False, fontsize=12, columnspacing=2.2, handlelength=3.0,
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
        ax_a.set_xlabel("OMP threads (MPI = 1, N = 64 H$_2$O)")
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
    _panel_letter(ax_a, "a", yoff=1.18)

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
    _panel_letter(ax_b, "b", yoff=1.18)

    if panel_a_handles:
        ax_a.legend(panel_a_handles,
                    [h.get_label() for h in panel_a_handles],
                    loc="lower center", bbox_to_anchor=(0.5, 1.04),
                    ncol=3, frameon=False, fontsize=9,
                    columnspacing=1.0, handlelength=2.0, handletextpad=0.5)
    if panel_b_handles:
        ax_b.legend(panel_b_handles, panel_b_labels,
                    loc="lower center", bbox_to_anchor=(0.5, 1.04),
                    ncol=len(panel_b_handles), frameon=False, fontsize=9,
                    columnspacing=1.0, handlelength=2.0, handletextpad=0.5)

    fig.tight_layout(rect=[0, 0, 1, 0.88])
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
        "__powr8i4":                        2.55,
        "OMP runtime (atomics / barriers)": 0.00,
        "hermite_spline_value":             0.00,
        "other / user code":               48.94,
    },
    "native-spline (pre-port)": {
        "MPI broadcast":                   34.77,
        "libm tanh / exp (NN activ.)":      0.14,
        "nnp_calc_ang":                     2.64,
        "pbc2 (O(N$^2$) neighbour cost)":   0.00,
        "__powr8i4":                        5.24,
        "OMP runtime (atomics / barriers)": 0.00,
        "hermite_spline_value":            10.88,
        "other / user code":               46.33,
    },
    "native-spline-omp": {
        "MPI broadcast":                   27.21,
        "libm tanh / exp (NN activ.)":      0.21,
        "nnp_calc_ang":                     1.65,
        "pbc2 (O(N$^2$) neighbour cost)":   0.00,
        "__powr8i4":                        3.22,
        "OMP runtime (atomics / barriers)": 28.44,
        "hermite_spline_value":             5.98,
        "other / user code":               33.29,
    },
}
MAQAO_AAEFF = {"master": 52.92, "native-spline-omp": 78.17}

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


FIGURES = [
    (1, fig1_algorithmic_complexity),
    (2, fig2_strong_scaling),
    (3, fig3_openmp_threading),
    (4, fig4_statistical_significance),
    (5, fig5_maqao_microarchitectural),
    (6, fig6_omp_size_heatmap),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", type=int, choices=[1, 2, 3, 4, 5], default=None,
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
