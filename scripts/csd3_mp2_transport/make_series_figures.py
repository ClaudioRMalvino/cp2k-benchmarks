#!/usr/bin/env python3
"""RPA concentration-series pooling + report figures (fig9-fig12).

Pools the per-segment CP2K analysis products in results/mp2_transport/
(written on CSD3 by analyze_onsager_cp2k.py / analyze_viscosity_cp2k.py) and
overlays the RPA series on the Madrid-2019 curves of make_figures.py, same
style and experimental anchors. Headline datasets:

  1 m: rpa_gpu_cube3_full_w1040_kappa.npz    + rpa_gpu_cube3_final_eta.npz
  2 m: rpa_gpu_cube3_2m_5x80_w1040_kappa.npz + rpa_gpu_cube3_2m_final_eta.npz
  4 m: rpa_gpu_cube3_4m_5x80_w1040_kappa.npz + rpa_gpu_cube3_4m_final_eta.npz
  MP2 anchor at 1 m, same pipeline: mp2_anchor_cube3_{kappa,eta}.npz

Model encoding across figures: Madrid-2019 = hollow markers + solid lines,
RPA = filled markers + dashed lines (x-offset so error bars stay legible),
MP2 anchor = filled indigo marker at 1 m only. Channel/species hues follow
make_figures.py so a channel keeps its colour across classical and RPA.

kappa and t_Na experimental anchors are the pinned Chambers-Stokes 1956 /
Smits-Duyvis 1966 values of make_figures.py. The eta anchors are UNPINNED
placeholders (Kestin-Khalifa-Correia 1981 ballpark, 25 C) and say so in the
legend; pin them from the primary source before any report use.

Also prints: fit-window sensitivity (10-30 vs 10-40 ps), the 1 m three-way
Madrid/MP2/RPA block, the D_w timescale check, and the timings-ledger
aggregation. Writes results/mp2_transport/rpa_series_summary.csv.
"""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_REPO = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                      "..", ".."))
MP2T = os.path.join(_REPO, "results", "mp2_transport")
MADRID = os.path.join(_REPO, "results", "lammps_madrid",
                      "conductivity", "analysis", "kappa_vs_c.npz")
OUT = os.path.join(MP2T, "figures")

# style block of make_figures.py (kept verbatim so the figures match)
CAM_BLUE, SLATE3 = "#00BDB6", "#546072"
C_NA, C_CL, C_O = CAM_BLUE, "#4DB78C", "#CD3572"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
INDIGO = "#5366E0"                       # MP2-anchor entity colour (renders' Na)

MOL_TO_M = {0.25: 0.231, 0.5: 0.491, 1.0: 0.960, 2.0: 1.915, 4.0: 3.666}
EXPT_C = {0.25: 0.248, 0.5: 0.494, 1.0: 0.979, 2.0: 1.920, 4.0: 3.686}
EXPT_KAPPA = {0.25: 2.48, 0.5: 4.63, 1.0: 8.42, 2.0: 14.50, 4.0: 22.04}
SD66_M = np.array([0.024, 0.048, 0.096, 0.192, 0.384, 0.768,
                   1.536, 3.072, 6.144])
SD66_LOGA = np.array([-1.683, -1.403, -1.126, -0.848, -0.573, -0.292,
                      0.0036, 0.3436, 0.7897])
M_NA, M_CL = 0.0229898, 0.0354530
E2, KB = (1.602176634e-19) ** 2, 1.380649e-23
CONV = 1.0e-8 * 1e9                      # A^2/ps -> 1e-9 m^2/s
DW_EXP = 2.299                           # pure water, 1e-9 m^2/s, 298.15 K (Holz 2000)

# UNPINNED eta placeholders (mPa s, ~25 C, Kestin-Khalifa-Correia 1981
# ballpark) -- pin from the primary source before report use.
ETA_EXP_UNPINNED = {1.0: 0.973, 2.0: 1.085, 4.0: 1.38}

# GK plateau window. The pipeline default (2-10 ps, stored in the npz) is
# biased low at RPA 1 m: the running integral is still rising there, the
# per-segment late-early difference is +0.31 mPa s at 7.9 sigma. Over
# 10-18 ps the slope of the segment-mean curve is zero within segment noise
# for 2 m / 4 m / the anchor and the anchor value moves by less than its
# SEM, so 10-18 ps is the headline window for all four datasets; the 2-10 ps
# value is printed alongside as the sensitivity.
ETA_WINDOW = (10.0, 18.0)

SERIES = {1.0: ("rpa_gpu_cube3_full_w1040_kappa.npz", "rpa_gpu_cube3_final_eta.npz"),
          2.0: ("rpa_gpu_cube3_2m_5x80_w1040_kappa.npz", "rpa_gpu_cube3_2m_final_eta.npz"),
          4.0: ("rpa_gpu_cube3_4m_5x80_w1040_kappa.npz", "rpa_gpu_cube3_4m_final_eta.npz")}
WINDOWS = {1.0: ("rpa_gpu_cube3_5x80_w1030_kappa.npz", "rpa_gpu_cube3_5x80_w1040_kappa.npz"),
           2.0: ("rpa_gpu_cube3_2m_5x80_w1030_kappa.npz", "rpa_gpu_cube3_2m_5x80_w1040_kappa.npz"),
           4.0: ("rpa_gpu_cube3_4m_5x80_w1030_kappa.npz", "rpa_gpu_cube3_4m_5x80_w1040_kappa.npz")}
ANCHOR = ("mp2_anchor_cube3_kappa.npz", "mp2_anchor_cube3_eta.npz")
MASTER_S_PER_STEP = None                 # read from timings.csv (master_ref row)


def tna_expt_hittorf(mol):
    loga = np.interp(np.log(np.asarray(mol, float)), np.log(SD66_M), SD66_LOGA)
    return 0.3720 - 0.0118 * loga


def style(ax, xlabel=None, ylabel=None, title=None):
    ax.grid(True, color=GRID, linewidth=0.7, alpha=0.9)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=INK2, labelsize=9)
    if xlabel:
        ax.set_xlabel(xlabel, color=INK, fontsize=10)
    if ylabel:
        ax.set_ylabel(ylabel, color=INK, fontsize=10)
    if title:
        ax.set_title(title, color=INK, fontsize=11, loc="center", pad=8)


def savefig(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"),
                    dpi=400 if ext == "png" else None,
                    facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {name}.png/.pdf")


def sem(a):
    a = np.asarray(a)
    return a.std(ddof=1) / np.sqrt(len(a)) if len(a) > 1 else np.nan


def pooled(runs_mol, y):
    mols = sorted(set(runs_mol))
    mean = np.array([y[runs_mol == m].mean() for m in mols])
    err = np.array([y[runs_mol == m].std(ddof=1) / np.sqrt((runs_mol == m).sum())
                    for m in mols])
    return np.array(mols), mean, err


def load_kappa(fname, mol):
    """Per-segment arrays + derived channels for one (level, molality)."""
    z = np.load(os.path.join(MP2T, fname))
    d = {k: z[k] for k in z.files}
    pref = E2 / (KB * float(d["temp"]) * float(d["V"]) * 1e-30) * 1e-8 / 6.0
    d["pref"], d["mol"] = pref, mol
    d["dNE"] = 1.0 - d["kOns"] / d["kNE"]
    d["kNaNa"], d["kClCl"], d["kNaCl"] = (pref * d[s] for s in ("sNaNa", "sClCl", "sNaCl"))
    d["dNaNa"] = d["kNaNa"] - pref * 6.0 * int(d["nNa"]) * d["DNa"]
    d["dClCl"] = d["kClCl"] - pref * 6.0 * int(d["nCl"]) * d["DCl"]
    d["chNaNa"] = -d["dNaNa"] / d["kNE"]
    d["chClCl"] = -d["dClCl"] / d["kNE"]
    d["chNaCl"] = 2.0 * d["kNaCl"] / d["kNE"]
    s_sum = d["sNaNa"] + d["sClCl"] - 2.0 * d["sNaCl"]
    d["tH"] = d["tNa"] + mol * (M_NA * (d["sNaNa"] - d["sNaCl"])
                                - M_CL * (d["sClCl"] - d["sNaCl"])) / s_sum
    d["rNa"], d["rCl"] = d["DNa"] / d["DO"], d["DCl"] / d["DO"]
    return d


def load_eta(fname):
    z = np.load(os.path.join(MP2T, fname))
    d = {k: z[k] for k in z.files}
    d["etas_default"], d["plateau_default"] = d["etas"], d["plateau"]
    m = (d["t_ps"] >= ETA_WINDOW[0]) & (d["t_ps"] <= ETA_WINDOW[1])
    d["etas"] = d["curves"][:, m].mean(axis=1)
    d["plateau"] = np.array(ETA_WINDOW)
    return d


def report_block(name, d, eta):
    n = len(d["kOns"])
    print(f"\n--- {name} ({n} segments, fit {d['tmin']:g}-{d['tmax']:g} ps, "
          f"T={float(d['temp']):g} K, {int(d['nNa'])} ion pairs) ---")
    for lab, arr, c in (("D_w ", d["DO"], CONV), ("D_Na", d["DNa"], CONV),
                        ("D_Cl", d["DCl"], CONV)):
        print(f"  {lab} = {arr.mean()*c:6.3f} +- {sem(arr)*c:5.3f}  x1e-9 m^2/s")
    for lab, arr in (("D_Na/D_w", d["rNa"]), ("D_Cl/D_w", d["rCl"])):
        print(f"  {lab} = {arr.mean():6.3f} +- {sem(arr):5.3f}")
    for lab, arr in (("kappa_NE ", d["kNE"]), ("kappa_Ons", d["kOns"]),
                     ("Delta_NE ", d["dNE"]), ("t_Na bary", d["tNa"]),
                     ("t_Na Hitt", d["tH"]), ("ch Na-Na ", d["chNaNa"]),
                     ("ch Cl-Cl ", d["chClCl"]), ("ch Na-Cl ", d["chNaCl"])):
        print(f"  {lab} = {arr.mean():6.3f} +- {sem(arr):5.3f}")
    if eta is not None:
        e = eta["etas"] * 1e3
        e0 = eta["etas_default"] * 1e3
        lo, hi = eta["plateau"]
        l0, h0 = eta["plateau_default"]
        print(f"  eta      = {e.mean():6.3f} +- {sem(e):5.3f}  mPa s  "
              f"(plateau {lo:g}-{hi:g} ps; {e0.mean():.3f} +- {sem(e0):.3f} "
              f"over the {l0:g}-{h0:g} ps default)")


def window_sensitivity():
    print("\n=== fit-window sensitivity (5x80 segment sets, 10-30 vs 10-40 ps) ===")
    for mol, (f30, f40) in sorted(WINDOWS.items()):
        d30, d40 = load_kappa(f30, mol), load_kappa(f40, mol)
        for lab, k in (("kOns", "kOns"), ("dNE", "dNE"), ("t_Na", "tNa")):
            print(f"  m={mol:g} {lab:>5}: {d30[k].mean():+6.3f}+-{sem(d30[k]):.3f} (w1030)"
                  f"  vs {d40[k].mean():+6.3f}+-{sem(d40[k]):.3f} (w1040)"
                  f"   shift {(d40[k].mean()-d30[k].mean()):+.3f}")


def ledger():
    """Aggregate timings.csv; rows with wall_s <= 10 are requeue bookkeeping."""
    rows = []
    with open(os.path.join(MP2T, "timings.csv")) as fh:
        for r in csv.DictReader(fh):
            r["wall_s"], r["steps"] = float(r["wall_s"]), float(r["steps"])
            r["s_per_step"] = float(r["s_per_step"])
            rows.append(r)
    master = [r for r in rows if r["stage"] == "master_ref"]
    sps_master = master[0]["s_per_step"] if master else np.nan
    keep = [r for r in rows if r["wall_s"] > 10
            and r["stage"] not in ("master_ref", "smoke", "gpu_smoke")]
    dropped = len(rows) - len(keep) - len(master) - sum(
        r["stage"] in ("smoke", "gpu_smoke") for r in rows)
    print(f"\n=== cost ledger (timings.csv; {dropped} requeue-bookkeeping rows "
          f"with wall<=10 s dropped; master s/step = {sps_master:g}) ===")
    print(f"  {'conc':9s} {'partition':9s} {'jobs':>4s} {'steps':>9s} "
          f"{'wall h':>8s} {'proj. master h':>14s} {'ratio':>6s}")
    tot_w = tot_p = 0.0
    for conc in sorted({r["conc"] for r in keep}):
        for part in sorted({r["partition"] for r in keep if r["conc"] == conc}):
            g = [r for r in keep if r["conc"] == conc and r["partition"] == part]
            wall_h = sum(r["wall_s"] for r in g) / 3600.0
            steps = sum(r["steps"] for r in g)
            proj_h = steps * sps_master / 3600.0
            tot_w, tot_p = tot_w + wall_h, tot_p + proj_h
            print(f"  {conc:9s} {part:9s} {len(g):4d} {steps:9.0f} "
                  f"{wall_h:8.1f} {proj_h:14.1f} {proj_h/wall_h:6.1f}x")
    print(f"  {'TOTAL':9s} {'':9s} {'':4s} {'':9s} {tot_w:8.1f} {tot_p:14.1f} "
          f"{tot_p/tot_w:6.1f}x")
    print("  (wall-clock basis, one node per job: icelake rows are CPU node-"
          "hours; ampere rows used 1xA100 + 3 SPME ranks of a 4-GPU node. "
          "The projected column prices the same steps at the measured master "
          "s/step on one icelake node.)")


def fig_kappa_series(mad, rpa, mp2):
    mols, kNE, kNE_s = pooled(mad["runs_mol"], mad["runs_kNE"])
    _, kOns, kOns_s = pooled(mad["runs_mol"], mad["runs_kOns"])
    c = np.array([MOL_TO_M[m] for m in mols])
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    ax.errorbar(c, kNE, yerr=kNE_s, fmt="^-", ms=6, color=C_O, mfc="white",
                mew=1.5, lw=1.4, elinewidth=1.0, capsize=2,
                label="Madrid NE")
    ax.errorbar(c, kOns, yerr=kOns_s, fmt="o-", ms=6.5, color=CAM_BLUE,
                mfc="white", mew=1.6, lw=1.8, elinewidth=1.0, capsize=2,
                label="Madrid Onsager")
    cr = np.array([MOL_TO_M[m] for m in sorted(rpa)])
    for key, color, mk, lab in (("kNE", C_O, "^", "RPA NE"),
                                ("kOns", CAM_BLUE, "o", "RPA Onsager")):
        y = [rpa[m][key].mean() for m in sorted(rpa)]
        ys = [sem(rpa[m][key]) for m in sorted(rpa)]
        ax.errorbar(cr + 0.07, y, yerr=ys, fmt=mk + "--", ms=6.5, color=color,
                    mew=0, lw=1.3, elinewidth=1.0, capsize=2, label=lab)
    ax.errorbar([MOL_TO_M[1.0] - 0.07], [mp2["kOns"].mean()],
                yerr=[sem(mp2["kOns"])], fmt="s", ms=6.5, color=INDIGO, mew=0,
                elinewidth=1.0, capsize=2, label="MP2 anchor Onsager")
    ce = [EXPT_C[m] for m in EXPT_KAPPA]
    ax.plot(ce, list(EXPT_KAPPA.values()), marker="o", ms=5, ls="",
            color=SLATE3, label="Experiment", zorder=5)
    h, l = ax.get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 1.09), ncol=3,
               frameon=False, fontsize=8.5, columnspacing=1.4,
               handlelength=2.2, handletextpad=0.5, labelcolor=INK2)
    ax.set_xlim(0, 3.9)
    ax.set_ylim(0, 24)
    style(ax, "c (mol/L)", r"$\kappa$ (S/m)")
    savefig(fig, "fig9_kappa_series")


def fig_dne_series(mad, rpa, mp2):
    E2KBT = E2 / (KB * 298.15)
    pref = E2KBT / (mad["runs_V"] * 1.0e-30) * 1.0e-8 / 6.0
    kNaNa, kClCl, kNaCl = (pref * mad[f"runs_s{s}"] for s in ("NaNa", "ClCl", "NaCl"))
    dNaNa = kNaNa - pref * 6.0 * mad["runs_nNa"] * mad["runs_DNa"]
    dClCl = kClCl - pref * 6.0 * mad["runs_nCl"] * mad["runs_DCl"]
    kNE, mol = mad["runs_kNE"], mad["runs_mol"]
    c = np.array([MOL_TO_M[m] for m in sorted(set(mol))])
    cr = np.array([MOL_TO_M[m] for m in sorted(rpa)])

    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.2), sharex=True, sharey=True)
    parts = ((-dNaNa / kNE, "chNaNa", CAM_BLUE, "s", "Na–Na distinct"),
             (-dClCl / kNE, "chClCl", C_CL, "^", "Cl–Cl distinct"),
             (2 * kNaCl / kNE, "chNaCl", C_O, "D", "Na–Cl distinct (pairing)"),
             (1.0 - mad["runs_kOns"] / kNE, "dNE", SLATE3, "o",
              r"total $\Delta_{\mathrm{NE}}$"))
    for ax, (y_mad, key, color, mk, lab) in zip(axes.flat, parts):
        _, ym, ys = pooled(mol, y_mad)
        ax.errorbar(c, ym, yerr=ys, fmt=mk + "-", ms=5.5, color=color,
                    mfc="white", mew=1.5, lw=1.4, elinewidth=1.0, capsize=2,
                    label="Madrid-2019")
        y = [rpa[m][key].mean() for m in sorted(rpa)]
        ys = [sem(rpa[m][key]) for m in sorted(rpa)]
        ax.errorbar(cr + 0.07, y, yerr=ys, fmt=mk, ms=6, color=color,
                    mew=0, elinewidth=1.0, capsize=2, label="RPA")
        ax.errorbar([MOL_TO_M[1.0] - 0.07], [mp2[key].mean()],
                    yerr=[sem(mp2[key])], fmt="s", ms=6, color=INDIGO, mew=0,
                    elinewidth=1.0, capsize=2, label="MP2 anchor")
        ax.axhline(0, color=MUTED, lw=0.8)
        ax.set_xlim(0, 3.9)
        style(ax, None, None, lab)
    for ax in axes[1]:
        ax.set_xlabel("c (mol/L)", color=INK, fontsize=10)
    for ax in axes[:, 0]:
        ax.set_ylabel(r"contribution to $\Delta_{\mathrm{NE}}$",
                      color=INK, fontsize=10)
    h, l = axes.flat[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 1.02), ncol=3,
               frameon=False, fontsize=9.5, columnspacing=2.0,
               handlelength=2.2, handletextpad=0.5, labelcolor=INK2)
    fig.text(0.5, 0.955, "hollow: Madrid-2019 (lines)   filled: RPA / MP2 (points)",
             ha="center", fontsize=8, color=INK2)
    savefig(fig, "fig10_dne_series")


def fig_tna_series(mad, rpa, mp2):
    mols = np.array(sorted(set(mad["runs_mol"])))
    c = np.array([MOL_TO_M[m] for m in mols])
    sNN, sCC, sNC = (mad[f"runs_s{s}"] for s in ("NaNa", "ClCl", "NaCl"))
    s_sum = sNN + sCC - 2.0 * sNC
    tH_runs = (mad["runs_tNa"]
               + mad["runs_mol"] * (M_NA * (sNN - sNC) - M_CL * (sCC - sNC)) / s_sum)
    _, tH, tH_s = pooled(mad["runs_mol"], tH_runs)
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    ax.errorbar(c, tH, yerr=tH_s, fmt="D-", ms=5.5, color=CAM_BLUE,
                mfc="white", mew=1.3, lw=1.3, elinewidth=1.0, capsize=2,
                label="Madrid, Hittorf")
    cr = np.array([MOL_TO_M[m] for m in sorted(rpa)])
    y = [rpa[m]["tH"].mean() for m in sorted(rpa)]
    ys = [sem(rpa[m]["tH"]) for m in sorted(rpa)]
    ax.errorbar(cr + 0.07, y, yerr=ys, fmt="D--", ms=5.5, color=CAM_BLUE,
                mew=0, lw=1.1, elinewidth=1.0, capsize=2, label="RPA, Hittorf")
    ax.errorbar([MOL_TO_M[1.0] - 0.07], [mp2["tH"].mean()],
                yerr=[sem(mp2["tH"])], fmt="s", ms=6, color=INDIGO, mew=0,
                elinewidth=1.0, capsize=2, label="MP2 anchor, Hittorf")
    ax.plot([EXPT_C[m] for m in mols], tna_expt_hittorf(mols), marker="o",
            ms=5, ls="", color=SLATE3, label="Experiment", zorder=5)
    ax.set_xlim(0, 3.9)
    h, l = ax.get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 1.05), ncol=2,
               frameon=False, fontsize=9, columnspacing=2.0,
               handlelength=2.4, handletextpad=0.6, labelcolor=INK2)
    style(ax, "c (mol/L)", r"$t_{\mathrm{Na}^+}$ (Hittorf frame)")
    savefig(fig, "fig11_tna_series")


def fig_eta_series(rpa_eta, mp2_eta):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.8, 4.0))
    # (a) eta(c)
    cr = np.array([MOL_TO_M[m] for m in sorted(rpa_eta)])
    y = [rpa_eta[m]["etas"].mean() * 1e3 for m in sorted(rpa_eta)]
    ys = [sem(rpa_eta[m]["etas"] * 1e3) for m in sorted(rpa_eta)]
    ax1.errorbar(cr, y, yerr=ys, fmt="o--", ms=6.5, color=CAM_BLUE, mew=0,
                 lw=1.3, elinewidth=1.0, capsize=2, label="RPA")
    ax1.errorbar([MOL_TO_M[1.0] - 0.07], [mp2_eta["etas"].mean() * 1e3],
                 yerr=[sem(mp2_eta["etas"] * 1e3)], fmt="s", ms=6.5,
                 color=INDIGO, mew=0, elinewidth=1.0, capsize=2,
                 label="MP2 anchor")
    ax1.plot([EXPT_C[m] for m in ETA_EXP_UNPINNED],
             list(ETA_EXP_UNPINNED.values()), marker="o", ms=5, ls="",
             color=SLATE3, label="Experiment (UNPINNED)", zorder=5)
    ax1.set_xlim(0, 3.9)
    ax1.legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="upper left")
    style(ax1, "c (mol/L)", r"$\eta$ (mPa s)", "(a)")
    # (b) GK running integrals, per-molality segment mean
    for m, color in ((1.0, CAM_BLUE), (2.0, C_CL), (4.0, C_O)):
        e = rpa_eta[m]
        ax2.plot(e["t_ps"], e["curves"].mean(axis=0) * 1e3, color=color,
                 lw=1.6, label=f"RPA {m:g} mol/kg")
        for cv in e["curves"]:
            ax2.plot(e["t_ps"], cv * 1e3, color=color, lw=0.5, alpha=0.25)
    ax2.plot(mp2_eta["t_ps"], mp2_eta["curves"].mean(axis=0) * 1e3,
             color=INDIGO, lw=1.4, ls="--", label="MP2 anchor 1 mol/kg")
    lo, hi = rpa_eta[1.0]["plateau"]
    ax2.axvspan(lo, hi, color=GRID, alpha=0.5, zorder=0)
    ax2.set_xlim(0, float(rpa_eta[1.0]["t_ps"].max()))
    ax2.legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="lower right")
    style(ax2, "t (ps)", r"$\eta_{\mathrm{GK}}(t)$ (mPa s)", "(b)")
    savefig(fig, "fig12_eta_series")


def write_summary(mad, rpa, rpa_eta, mp2, mp2_eta):
    path = os.path.join(MP2T, "rpa_series_summary.csv")
    cols = ("level,m_molkg,n_seg,fit_ps,Dw_1e9,Dw_sem,DNa_1e9,DNa_sem,"
            "DCl_1e9,DCl_sem,rNa,rNa_sem,rCl,rCl_sem,kNE_Sm,kNE_sem,"
            "kOns_Sm,kOns_sem,dNE,dNE_sem,chNaNa,chNaNa_sem,chClCl,"
            "chClCl_sem,chNaCl,chNaCl_sem,tNa_bary,tNa_bary_sem,"
            "tNa_hittorf,tNa_hittorf_sem,eta_mPas,eta_sem\n")
    def row(level, mol, d, eta):
        f = [f"{level}", f"{mol:g}", f"{len(d['kOns'])}",
             f"{d['tmin']:g}-{d['tmax']:g}"]
        for k in ("DO", "DNa", "DCl"):
            f += [f"{d[k].mean()*CONV:.4f}", f"{sem(d[k])*CONV:.4f}"]
        for k in ("rNa", "rCl", "kNE", "kOns", "dNE", "chNaNa", "chClCl",
                  "chNaCl", "tNa", "tH"):
            f += [f"{d[k].mean():.4f}", f"{sem(d[k]):.4f}"]
        if eta is not None:
            f += [f"{eta['etas'].mean()*1e3:.4f}", f"{sem(eta['etas']*1e3):.4f}"]
        else:
            f += ["", ""]
        return ",".join(f) + "\n"
    with open(path, "w") as fh:
        fh.write(cols)
        for m in sorted(rpa):
            fh.write(row("RPA", m, rpa[m], rpa_eta[m]))
        fh.write(row("MP2_anchor", 1.0, mp2, mp2_eta))
        # Madrid rows, pooled from kappa_vs_c.npz in the same shape
        for m in sorted(set(mad["runs_mol"])):
            g = mad["runs_mol"] == m
            d = {"tmin": 20, "tmax": 400, "kOns": mad["runs_kOns"][g],
                 "kNE": mad["runs_kNE"][g], "DO": mad["runs_DO"][g],
                 "DNa": mad["runs_DNa"][g], "DCl": mad["runs_DCl"][g],
                 "tNa": mad["runs_tNa"][g]}
            pref = (E2 / (KB * 298.15) / (mad["runs_V"][g] * 1e-30) * 1e-8 / 6.0)
            kNaNa, kClCl, kNaCl = (pref * mad[f"runs_s{s}"][g]
                                   for s in ("NaNa", "ClCl", "NaCl"))
            d["chNaNa"] = -(kNaNa - pref * 6 * mad["runs_nNa"][g]
                            * mad["runs_DNa"][g]) / d["kNE"]
            d["chClCl"] = -(kClCl - pref * 6 * mad["runs_nCl"][g]
                            * mad["runs_DCl"][g]) / d["kNE"]
            d["chNaCl"] = 2 * kNaCl / d["kNE"]
            d["dNE"] = 1.0 - d["kOns"] / d["kNE"]
            d["rNa"] = d["DNa"] / d["DO"]
            d["rCl"] = d["DCl"] / d["DO"]
            sNN, sCC, sNC = (mad[f"runs_s{s}"][g] for s in ("NaNa", "ClCl", "NaCl"))
            s_sum = sNN + sCC - 2 * sNC
            d["tH"] = d["tNa"] + m * (M_NA * (sNN - sNC) - M_CL * (sCC - sNC)) / s_sum
            fh.write(row("Madrid2019", m, d, None))
    print(f"wrote {path}")


def main():
    os.makedirs(OUT, exist_ok=True)
    mad = {k: v for k, v in np.load(MADRID, allow_pickle=True).items()}
    rpa = {m: load_kappa(f, m) for m, (f, _) in SERIES.items()}
    rpa_eta = {m: load_eta(e) for m, (_, e) in SERIES.items()}
    mp2 = load_kappa(ANCHOR[0], 1.0)
    mp2_eta = load_eta(ANCHOR[1])

    print("=== RPA series (headline datasets) ===")
    for m in sorted(rpa):
        report_block(f"RPA {m:g} mol/kg", rpa[m], rpa_eta[m])
    report_block("MP2 anchor 1 mol/kg", mp2, mp2_eta)

    print("\n=== three-way at 1 mol/kg (Madrid pooled n=10, MP2/RPA n=5) ===")
    g = mad["runs_mol"] == 1.0
    print("  (Madrid D ratios here are finite-L run values; the report's "
          "headline Madrid ratios use the Yeh-Hummer D_0)")
    mad_dne = 1.0 - mad["runs_kOns"][g] / mad["runs_kNE"][g]
    print(f"  kOns : Madrid {mad['runs_kOns'][g].mean():.3f}+-{sem(mad['runs_kOns'][g]):.3f}"
          f"  MP2 {mp2['kOns'].mean():.3f}+-{sem(mp2['kOns']):.3f}"
          f"  RPA {rpa[1.0]['kOns'].mean():.3f}+-{sem(rpa[1.0]['kOns']):.3f}"
          f"  expt {EXPT_KAPPA[1.0]}")
    print(f"  dNE  : Madrid {mad_dne.mean():.3f}+-{sem(mad_dne):.3f}"
          f"  MP2 {mp2['dNE'].mean():.3f}+-{sem(mp2['dNE']):.3f}"
          f"  RPA {rpa[1.0]['dNE'].mean():.3f}+-{sem(rpa[1.0]['dNE']):.3f}")
    for lab, k in (("rNa", "rNa"), ("rCl", "rCl")):
        madr = mad[f"runs_D{'Na' if k=='rNa' else 'Cl'}"][g] / mad["runs_DO"][g]
        print(f"  {lab}  : Madrid {madr.mean():.3f}+-{sem(madr):.3f}"
              f"  MP2 {mp2[k].mean():.3f}+-{sem(mp2[k]):.3f}"
              f"  RPA {rpa[1.0][k].mean():.3f}+-{sem(rpa[1.0][k]):.3f}")
    print(f"  D_w  : MP2 {mp2['DO'].mean()*CONV:.3f}  RPA {rpa[1.0]['DO'].mean()*CONV:.3f}"
          f"  x1e-9 m^2/s;  pure-water expt {DW_EXP} (298.15 K, Holz 2000)")
    print(f"  timescale factor vs pure-water D_w: "
          f"MP2 {DW_EXP/(mp2['DO'].mean()*CONV):.2f}x  "
          f"RPA {DW_EXP/(rpa[1.0]['DO'].mean()*CONV):.2f}x  "
          f"(1 m solution value slightly lower -> both factors upper bounds)")

    window_sensitivity()
    ledger()
    print("\nNOT in this data drop: RPA g_NaCl(r) / PMF / n_CIP products "
          "(figure v and the tab:pmf RPA row still need the cluster-side RDFs).")
    print("UNPINNED: the eta experimental anchors in fig12a.")

    fig_kappa_series(mad, rpa, mp2)
    fig_dne_series(mad, rpa, mp2)
    fig_tna_series(mad, rpa, mp2)
    fig_eta_series(rpa_eta, mp2_eta)
    write_summary(mad, rpa, rpa_eta, mp2, mp2_eta)


if __name__ == "__main__":
    main()
