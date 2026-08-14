#!/usr/bin/env python3
"""Single-panel speed-up vs system size for the A0 poster.
Same data and style as thesis fig 4.1(d); writes plots/poster/poster_speedup.pdf."""
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "benchmark_figures"))
import thesis_figures as tf

REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
tf.RESULTS_ROOT = os.path.join(REPO, "results")
OUT_DIR = os.path.join(REPO, "plots", "poster")

import matplotlib.pyplot as plt


def main():
    data, _ = tf.load_all_data()
    m = data["upstream master"]["size"]

    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    for name in ("native-spline", "native-spline-omp"):
        d = data[name]
        merged = pd.merge(m[["n_molecules", "tps_mean"]],
                          d["size"][["n_molecules", "tps_mean"]],
                          on="n_molecules", suffixes=("_m", "_o"))
        speedup = merged["tps_mean_m"] / merged["tps_mean_o"]
        ax.plot(tf._atoms(merged["n_molecules"]), speedup,
                marker=d["marker"], color=d["color"], lw=2.0, ms=7.5,
                label=d["label"])
    baseline = ax.axhline(1.0, color=tf.CAMBRIDGE["slate_3"], ls="--", lw=1.2,
                          label=r"master baseline ($1\times$)")
    ax.set(xscale="log")
    ax.set_xlabel(r"Number of atoms, $N$")
    ax.set_ylabel("Speed-up vs master")
    tf._log_decade_ticks(ax, "x", atom_fmt=True)

    fig.legend(loc="upper center", bbox_to_anchor=(0.5, 1.00), ncol=3,
               frameon=False, fontsize=10.5, columnspacing=1.2,
               handlelength=2.0, handletextpad=0.5)
    fig.tight_layout(rect=[0, 0, 1, 0.90])

    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, "poster_speedup.pdf")
    fig.savefig(path)
    print(f"saved {path}")


if __name__ == "__main__":
    main()
