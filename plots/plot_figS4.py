#!/usr/bin/env python3
"""Fig. S4 (Morawietz et al., PNAS 2016) replication plots for the CP2K NNP
chebyshev branch.

Consumes the analysis CSVs written by aggregate_figS4.py:
  figS4_summary.csv                      (eta, D_PBC, D_0, t/step per size)
  msd_<branch>_N<size>.csv               (lag_ps, msd_mean, msd_sem)
  running_eta_<branch>_N<size>.csv       (lag_fs, running_eta_mean, sem)
  acf_<branch>_N<size>.csv               (lag_fs, acf_norm_mean, sem)

Produces three figures the STAGE-3 wrapper expects:
  figS4_replication.png   D_PBC vs 1/L finite-size scaling -> D_0; eta vs size
  figS4_accuracy.png      MSD(t) diffusive regime + running eta(t) convergence
  figS4_performance.png   chebyshev production t/step vs system size

Yeh-Hummer convention matches compute_diffusion.py: L = V^(1/3),
xi = 2.837297, D_0 = D_PBC + k_B T xi / (6 pi eta L).
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CAMBRIDGE = {
    "blue_warm": "#00BDB6", "blue_dark": "#133844", "crest": "#FD8153",
    "crest_dark": "#DD3025", "indigo": "#5366E0", "purple": "#A368DF",
    "green": "#4DB78C", "slate_2": "#B5BDC8", "slate_3": "#546072",
    "slate_4": "#232830", "cherry": "#CD3572",
}
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Liberation Serif", "serif"],
    "font.size": 12, "axes.labelsize": 13, "axes.titlesize": 13,
    "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 10,
    "text.color": CAMBRIDGE["slate_4"], "axes.edgecolor": CAMBRIDGE["slate_4"],
    "axes.labelcolor": CAMBRIDGE["slate_4"], "xtick.color": CAMBRIDGE["slate_4"],
    "ytick.color": CAMBRIDGE["slate_4"], "axes.linewidth": 0.8,
    "axes.grid": True, "grid.color": CAMBRIDGE["slate_2"], "grid.linestyle": "--",
    "grid.linewidth": 0.5, "grid.alpha": 0.55, "lines.linewidth": 1.6,
    "lines.markersize": 6.5, "mathtext.fontset": "cm",
    "savefig.dpi": 300, "savefig.bbox": "tight", "pdf.fonttype": 42,
})

BASE_L = 12.42
MULT = {64: (1, 1, 1), 128: (2, 1, 1), 256: (2, 2, 1),
        512: (2, 2, 2), 1024: (4, 2, 2)}
SIZE_STYLE = {
    64:  dict(color=CAMBRIDGE["blue_warm"], marker="o"),
    128: dict(color=CAMBRIDGE["green"],     marker="^"),
    256: dict(color=CAMBRIDGE["indigo"],    marker="s"),
    512: dict(color=CAMBRIDGE["crest"],     marker="D"),
    1024: dict(color=CAMBRIDGE["crest_dark"], marker="v"),
}


def L_of(size):
    return BASE_L * float(np.prod(MULT[size])) ** (1.0 / 3.0)


def read_summary(path):
    rows = []
    with open(path) as f:
        hdr = f.readline().strip().split(",")
        for line in f:
            p = line.strip().split(",")
            if len(p) != len(hdr):
                continue
            d = dict(zip(hdr, p))
            try:
                d["n_molecules"] = int(d["n_molecules"])
                for k in ("eta_mPa_s", "eta_sem", "D_PBC_ang2_ps", "D_PBC_sem",
                          "D_0_ang2_ps", "D_0_sem", "time_per_step_s",
                          "time_per_step_sem"):
                    d[k] = float(d[k]) if d.get(k) not in (None, "", "nan") else float("nan")
            except (ValueError, KeyError):
                continue
            rows.append(d)
    return rows


def read_curve(path):
    if not os.path.exists(path):
        return None
    try:
        return np.genfromtxt(path, delimiter=",", names=True)
    except Exception:
        return None


def wls_fit(x, y, yerr):
    """Weighted linear fit y = a + b x; returns a, b, sigma_a."""
    w = 1.0 / np.where(yerr > 0, yerr, np.nanmean(yerr) or 1.0) ** 2
    X = np.vstack([np.ones_like(x), x]).T
    WX = X * w[:, None]
    cov = np.linalg.inv(X.T @ WX)
    beta = cov @ (WX.T @ y)
    return beta[0], beta[1], np.sqrt(cov[0, 0])


def fig_replication(rows, branch, out_dir, ref_eta, ref_d0):
    rows = [r for r in rows if r["branch"] == branch]
    rows.sort(key=lambda r: r["n_molecules"])
    sizes = [r["n_molecules"] for r in rows]
    invL = np.array([1.0 / L_of(s) for s in sizes])
    dpbc = np.array([r["D_PBC_ang2_ps"] for r in rows])
    dpbc_e = np.array([r["D_PBC_sem"] for r in rows])
    d0 = np.array([r["D_0_ang2_ps"] for r in rows])
    d0_e = np.array([r["D_0_sem"] for r in rows])
    eta = np.array([r["eta_mPa_s"] for r in rows])
    eta_e = np.array([r["eta_sem"] for r in rows])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.6, 3.7))

    # (a) finite-size scaling D_PBC vs 1/L, with WLS fit -> D_0 intercept
    for s, x, y, ye in zip(sizes, invL, dpbc, dpbc_e):
        st = SIZE_STYLE.get(s, dict(color=CAMBRIDGE["slate_3"], marker="o"))
        ax1.errorbar(x, y, yerr=ye, color=st["color"], marker=st["marker"],
                     ls="none", capsize=3, label=f"N={s}")
    ax1.errorbar(invL, d0, yerr=d0_e, color=CAMBRIDGE["slate_3"], marker="x",
                 ls="none", capsize=3, label=r"$D_0$ (YH-corrected)")
    if len(invL) >= 2 and np.all(np.isfinite(dpbc)):
        a, b, sa = wls_fit(invL, dpbc, dpbc_e)
        xs = np.linspace(0, invL.max() * 1.05, 50)
        ax1.plot(xs, a + b * xs, color=CAMBRIDGE["blue_dark"], lw=1.3, ls="-",
                 label=rf"fit: $D_0$={a:.3f}$\pm${sa:.3f}")
        ax1.axhline(a, color=CAMBRIDGE["blue_dark"], lw=0.7, ls=":")
    if ref_d0:
        ax1.axhline(ref_d0, color=CAMBRIDGE["cherry"], lw=1.2, ls="--",
                    label=f"ref $D_0$={ref_d0:.3f}")
    ax1.set_xlabel(r"$1/L$  [$\mathrm{\AA}^{-1}$]")
    ax1.set_ylabel(r"$D$  [$\mathrm{\AA}^2$/ps]")
    ax1.set_xlim(left=0)
    ax1.set_title("(a)  finite-size diffusion", loc="left")
    ax1.legend(frameon=False, fontsize=8.5)

    # (b) viscosity vs system size (should be ~size-independent)
    for s, y, ye in zip(sizes, eta, eta_e):
        st = SIZE_STYLE.get(s, dict(color=CAMBRIDGE["slate_3"], marker="o"))
        ax2.errorbar(s, y, yerr=ye, color=st["color"], marker=st["marker"],
                     ls="none", capsize=3)
    if np.any(np.isfinite(eta)):
        m = np.nanmean(eta)
        ax2.axhline(m, color=CAMBRIDGE["slate_3"], lw=0.9, ls=":",
                    label=rf"mean={m:.3f} mPa$\cdot$s")
    if ref_eta:
        ax2.axhline(ref_eta, color=CAMBRIDGE["cherry"], lw=1.2, ls="--",
                    label=f"ref={ref_eta:.3f} mPa$\\cdot$s")
    ax2.set_xscale("log", base=2)
    ax2.set_xlabel("number of molecules")
    ax2.set_ylabel(r"shear viscosity $\eta$  [mPa$\cdot$s]")
    ax2.set_title("(b)  shear viscosity", loc="left")
    ax2.legend(frameon=False, fontsize=8.5)

    fig.suptitle("Fig. S4 replication (feature/nnp-chebyshev, RPBE-vdW water, 300 K)",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out_dir, f"figS4_replication.{ext}"))
    plt.close(fig)


def fig_accuracy(rows, branch, analysis_dir, out_dir):
    rows = [r for r in rows if r["branch"] == branch]
    sizes = sorted(r["n_molecules"] for r in rows)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.6, 3.7))
    for s in sizes:
        st = SIZE_STYLE.get(s, dict(color=CAMBRIDGE["slate_3"], marker="o"))
        msd = read_curve(os.path.join(analysis_dir, f"msd_{branch}_N{s}.csv"))
        if msd is not None:
            ax1.plot(msd["lag_ps"], msd["msd_mean"], color=st["color"], label=f"N={s}")
        reta = read_curve(os.path.join(analysis_dir, f"running_eta_{branch}_N{s}.csv"))
        if reta is not None:
            ax2.plot(reta["lag_fs"] / 1000.0, reta["running_eta_mean"],
                     color=st["color"], label=f"N={s}")
    ax1.set_xlabel("time  [ps]")
    ax1.set_ylabel(r"MSD  [$\mathrm{\AA}^2$]")
    ax1.set_title("(a)  mean-squared displacement", loc="left")
    ax1.legend(frameon=False, fontsize=9)
    ax2.set_xlabel("integration limit  [ps]")
    ax2.set_ylabel(r"running $\eta$  [mPa$\cdot$s]")
    ax2.set_title("(b)  Green-Kubo viscosity convergence", loc="left")
    ax2.legend(frameon=False, fontsize=9)
    fig.suptitle("Fig. S4 observables (feature/nnp-chebyshev)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out_dir, f"figS4_accuracy.{ext}"))
    plt.close(fig)


def fig_performance(rows, branch, out_dir):
    rows = [r for r in rows if r["branch"] == branch]
    rows.sort(key=lambda r: r["n_molecules"])
    sizes = np.array([r["n_molecules"] for r in rows])
    tps = np.array([r["time_per_step_s"] for r in rows])
    tps_e = np.array([r["time_per_step_sem"] for r in rows])
    m = np.isfinite(tps)
    if not np.any(m):
        return
    fig, ax = plt.subplots(figsize=(4.6, 3.7))
    ax.errorbar(sizes[m], tps[m], yerr=tps_e[m], color=CAMBRIDGE["blue_warm"],
                marker="o", capsize=3, label="chebyshev production")
    if m.sum() >= 2:
        b = np.polyfit(np.log(sizes[m]), np.log(tps[m]), 1)[0]
        xs = np.array([sizes[m].min(), sizes[m].max()], float)
        ax.plot(xs, tps[m][0] * (xs / sizes[m][0]) ** b, color=CAMBRIDGE["slate_3"],
                lw=1.0, ls="--", label=rf"$\propto N^{{{b:.2f}}}$")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("number of molecules")
    ax.set_ylabel("time / step  [s]")
    ax.set_title("Production cost vs system size", loc="left")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out_dir, f"figS4_performance.{ext}"))
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--analysis-dir", required=True)
    ap.add_argument("--plot-dir", required=True)
    ap.add_argument("--branch", default="feature-nnp-chebyshev")
    # Defaults are liquid-water EXPERIMENT (298 K): D ~ 0.23 Ang^2/ps,
    # eta ~ 0.896 mPa.s.  Override with the Morawietz RPBE-vdW model values.
    ap.add_argument("--ref-d0", type=float, default=0.23)
    ap.add_argument("--ref-eta", type=float, default=0.896)
    args = ap.parse_args()
    os.makedirs(args.plot_dir, exist_ok=True)

    summary_path = os.path.join(args.analysis_dir, "figS4_summary.csv")
    rows = read_summary(summary_path)
    if not rows:
        print(f"!! no usable rows in {summary_path}; nothing to plot")
        return 1
    print(f"  {len(rows)} (branch,size) rows from {summary_path}")
    fig_replication(rows, args.branch, args.plot_dir, args.ref_eta, args.ref_d0)
    fig_accuracy(rows, args.branch, args.analysis_dir, args.plot_dir)
    fig_performance(rows, args.branch, args.plot_dir)
    print(f"  wrote figS4_replication / figS4_accuracy / figS4_performance "
          f"(.png/.pdf) to {args.plot_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
