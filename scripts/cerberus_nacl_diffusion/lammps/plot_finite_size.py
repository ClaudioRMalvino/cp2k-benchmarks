#!/usr/bin/env python3
"""Supervisor-ready figures for the Madrid-2019 NaCl(aq) finite-size study.

Fig 1  D_PBC vs 1/L with the Yeh-Hummer fit per species (small multiples,
       one shared y per panel; water tight, ions honest about their noise).
Fig 2  MSD diagnostics: log-log local slope vs t (the "flat derivative" check)
       and the MSD curves themselves.

Reads the same .msd files as analyze_finite_size.py.
"""
import glob
import os
import re

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

KB = 1.380649e-23
XI = 2.837297
T_K = 298.15
ROOT = "/data/cerberus1/crm98/nacl_diffusion/lammps_madrid/finite_size"
OUT = "/data/cerberus1/crm98/nacl_diffusion/lammps_madrid/figures"
TMIN, TMAX = 200.0, 2000.0

SPECIES = ["O", "Na", "Cl"]
NICE = {"O": "water (O)", "Na": "Na$^+$", "Cl": "Cl$^-$"}
COL = {"O": "#2a78d6", "Na": "#1baf7a", "Cl": "#eda100"}
INK, INK2, GRID = "#1a1a19", "#5f5e56", "#e5e4dd"

plt.rcParams.update({
    "font.size": 11, "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "axes.linewidth": 0.8, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "axes.facecolor": "white",
})


def load_all():
    runs = []
    for d in sorted(glob.glob(os.path.join(ROOT, "L*"))):
        L = float(re.search(r"L([\d.]+)$", d).group(1))
        f = os.path.join(d, f"L{L:.2f}.msd")
        if not os.path.exists(f):
            continue
        raw = np.loadtxt(f, comments="#")
        t = raw[:, 0] * 2.0 / 1000.0          # 2 fs steps -> ps
        runs.append((L, t, raw[:, 1:4].T))
    runs.sort(key=lambda r: r[0])
    return runs


def fit_D(t, msd):
    m = (t >= TMIN) & (t <= TMAX)
    sl, _ = np.polyfit(t[m], msd[m], 1)
    return sl / 6.0


def local_slope(t, msd, nbin=16):
    """d log MSD / d log t in log-spaced windows."""
    m = (t > 1.0) & (msd > 0)
    lt, lm = np.log(t[m]), np.log(msd[m])
    edges = np.linspace(lt[0], lt[-1], nbin + 1)
    xc, sl = [], []
    for i in range(nbin):
        w = (lt >= edges[i]) & (lt < edges[i + 1])
        if w.sum() > 10:
            xc.append(np.exp(edges[i:i + 2].mean()))
            sl.append(np.polyfit(lt[w], lm[w], 1)[0])
    return np.array(xc), np.array(sl)


def main():
    os.makedirs(OUT, exist_ok=True)
    runs = load_all()
    Ls = np.array([r[0] for r in runs])
    invL = 1.0 / Ls
    D = {s: np.array([fit_D(t, msd[i]) for _, t, msd in runs])
         for i, s in enumerate(SPECIES)}

    # ---------- Fig 1: Yeh-Hummer extrapolation, small multiples ----------
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6), constrained_layout=True)
    x = np.linspace(0, invL.max() * 1.12, 50)
    for ax, s in zip(axes, SPECIES):
        y = D[s]
        sl, D0 = np.polyfit(invL, y, 1)
        pred = sl * invL + D0
        r2 = 1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)
        ax.grid(True, color=GRID, lw=0.7)
        ax.set_axisbelow(True)
        if r2 > 0.5 and sl < 0:                       # fit is meaningful
            # raw D_PBC fit line sloping down to the infinite-size intercept
            ax.plot(x, sl * x + D0, color=COL[s], lw=1.6, zorder=2)
            ax.scatter([0], [D0], marker="D", s=48, facecolor="white",
                       edgecolor=COL[s], lw=1.6, zorder=5)
            # YH-corrected: add the finite-size drag back using the fit slope
            # itself (no external viscosity). D_YH = D_PBC - slope*(1/L), which
            # by construction collapses onto the flat D_0 line.
            y_corr = y - sl * invL
            ax.axhline(D0, color=INK2, ls="--", lw=1.0, zorder=1)
            ax.scatter(invL, y_corr, s=46, facecolor="white",
                       edgecolor=COL[s], lw=1.6, zorder=4)
            slope_si = abs(sl) * 1e-18                # A^3/ps -> m^3/s
            eta = KB * T_K * XI / (6 * np.pi * slope_si) * 1e3
            note = (f"$D_0$ = {D0:.3f} Å$^2$/ps\n"
                    f"$R^2$ = {r2:.2f},  $\\eta_{{fit}}$ = {eta:.2f} mPa·s")
        else:
            note = f"no reliable fit\n($R^2$ = {max(r2,0):.2f})"
        # raw D_PBC measured points (filled)
        ax.scatter(invL, y, s=46, color=COL[s], zorder=3)
        ax.text(0.04, 0.05, note, transform=ax.transAxes, fontsize=9.5,
                color=INK, va="bottom")
        ax.set_title(NICE[s], color=INK, fontsize=12)
        ax.set_xlabel("1/L  (Å$^{-1}$)")
        ax.set_xlim(-0.002, x.max())
    axes[0].set_ylabel("$D$  (Å$^2$/ps)")
    # legend: filled = raw, open = YH-corrected, diamond = infinite-size intercept
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", ls="none", mfc=COL["O"], mec=COL["O"],
               ms=7, label=r"$D_{\mathrm{PBC}}$ (raw, per box)"),
        Line2D([0], [0], marker="o", ls="none", mfc="white", mec=COL["O"],
               mew=1.6, ms=7, label=r"$D_{\mathrm{YH}}$ (corrected)"),
        Line2D([0], [0], marker="D", ls="none", mfc="white", mec=COL["O"],
               mew=1.6, ms=7, label=r"$D_0$ (intercept)"),
    ]
    axes[0].legend(handles=handles, loc="upper right", fontsize=8.5,
                   frameon=False)
    fig.suptitle("Yeh–Hummer finite-size extrapolation — Madrid-2019, 1.0 mol/kg NaCl(aq), 298.15 K",
                 fontsize=12, color=INK)
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/fig1_yehhummer.{ext}", dpi=200)
    plt.close(fig)

    # ---------- Fig 2: MSD + flat-derivative check (largest box) ----------
    L, t, msd = runs[-1]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.2, 3.6), constrained_layout=True)
    for ax in (a1, a2):
        ax.grid(True, color=GRID, lw=0.7)
        ax.set_axisbelow(True)
    label_shift = {"O": 1.6, "Na": 0.55, "Cl": 1.0}   # stagger colliding end labels
    for i, s in enumerate(SPECIES):
        a1.loglog(t[t > 0.1], msd[i][t > 0.1], color=COL[s], lw=1.6)
        a1.text(t[-1] * 1.12, msd[i][-1] * label_shift[s], NICE[s],
                color=COL[s], fontsize=10, va="center")
        xc, sl = local_slope(t, msd[i])
        a2.semilogx(xc, sl, color=COL[s], lw=1.6, marker="o", ms=3.5)
    a1.set_xlabel("t  (ps)"); a1.set_ylabel("MSD  (Å$^2$)")
    a1.set_xlim(0.1, t[-1] * 3.2)
    a1.set_title(f"MSD, L = {L:.2f} Å box", color=INK, fontsize=12)
    a2.axhline(1.0, color=INK2, lw=0.9, ls="--")
    a2.axvspan(TMIN, TMAX, color=GRID, alpha=0.5, zorder=0)
    a2.text(np.sqrt(TMIN * TMAX), 1.62, "fit window", ha="center",
            fontsize=9, color=INK2)
    a2.set_ylim(0.4, 1.9)
    a2.set_xlabel("t  (ps)")
    a2.set_ylabel(r"d log MSD / d log t")
    a2.set_title("diffusive check (slope → 1)", color=INK, fontsize=12)
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/fig2_msd_check.{ext}", dpi=200)
    plt.close(fig)

    print(f"wrote {OUT}/fig1_yehhummer.png/.pdf and fig2_msd_check.png/.pdf")


if __name__ == "__main__":
    main()
