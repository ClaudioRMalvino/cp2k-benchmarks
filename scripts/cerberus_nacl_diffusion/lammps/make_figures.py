#!/usr/bin/env python3
"""Paper/supervisor figures for the Madrid-2019 NaCl(aq) transport study.

Reads only finished analysis products + raw .rdf/.vacf files:
  replicas/analysis/replica_D.npz      (MSD self-diffusion, Yeh-Hummer)
  vacf_replica/L*_s*/*.vacf            (Green-Kubo D, corrected campaign)
  conductivity/analysis/kappa_vs_c.npz (Onsager/NE conductivity, 50 runs)
  conductivity/m*/L*_s*/*.rdf          (Na-Cl g(r) -> PMF)

Outputs PNG (400 dpi, slides) + PDF (vector, paper) per figure into --out.

Experimental anchors (298.15 K), pinned 2026-07-27 from the primary sources:
kappa at m = 0.25-4 from Chambers, Stokes & Stokes 1956 Table II (the actual
source behind the Blazquez 2023 compilation = their ref 143); Hittorf t_Na
from Smits & Duyvis 1966 (their eq 3 + Table I activities). Fig. 4b also
shows our t_Na transformed to the solvent (Hittorf) frame so the comparison
is like-for-like (cf. Shao et al. 2022, "Mind the Reference-Frame Gap").
"""
import argparse
import glob
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Data root: cerberus live tree if present, else the analysis products bundled
# in the repo (results/lammps_madrid), so figures regenerate on any machine.
# Override with LAMMPS_MADRID_ROOT.
_REPO = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                      "..", "..", ".."))
ROOT = next(p for p in (
    os.environ.get("LAMMPS_MADRID_ROOT"),
    "/data/cerberus1/crm98/nacl_diffusion/lammps_madrid",
    os.path.join(_REPO, "results", "lammps_madrid"),
) if p and os.path.isdir(p))
KT_KCAL = 0.0019872041 * 298.15          # k_B T in kcal/mol at 298.15 K

# Cambridge palette for non-species curves (cf. plots/ CAM dict); two-series
# figures pair cambridge blue with slate_3 (slate_4 #232830 reads as black
# at line weight, so the visible mid slate is used for dark series)
CAM_BLUE, SLATE3 = "#00BDB6", "#546072"
C_NA, C_CL, C_O = CAM_BLUE, "#4DB78C", "#CD3572"   # slots 1-3 (non-species curves)
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
# element colours: Cambridge indigo / green / cherry (cf. plots/ CAM palette)
SPECIES_COLOR = {"Na": "#5366E0", "Cl": "#4DB78C", "O": "#CD3572"}

MOL_TO_M = {0.25: 0.231, 0.5: 0.491, 1.0: 0.960, 2.0: 1.915, 4.0: 3.666}  # sim c

# --- pinned experimental anchors, 298.15 K ---------------------------------
# kappa: Chambers, Stokes & Stokes 1956 (J. Phys. Chem. 60, 985), Table II
# equivalent conductances Lambda(c) (Int-ohm units, stated error <=0.03%);
# kappa[S/m] = Lambda*c/10. c(m) from experimental densities 997.04/1036.21/
# 1072.27/1136.91 kg/m^3 at m=0/1/2/4 (expt column, Blazquez 2023 Table 1;
# m=0.25/0.5 interpolated, <0.1% error in c). Blazquez's own m=2/4 quotes
# (14.49/22.04) reproduce exactly; their m=1 quote (8.48) sits 0.7% above
# this like-for-like conversion (8.42).
EXPT_C     = {0.25: 0.248, 0.5: 0.494, 1.0: 0.979, 2.0: 1.920, 4.0: 3.686}
EXPT_KAPPA = {0.25: 2.48, 0.5: 4.63, 1.0: 8.42, 2.0: 14.50, 4.0: 22.04}

# t_Na, Hittorf (solvent) frame: Smits & Duyvis 1966 (J. Phys. Chem. 70,
# 2747) eq 3, t_Na = 0.3720 - 0.0118*log10(a+-), valid 0.024-6.144 mol/kg,
# with log10(a+-) from their Table I (Robinson-Stokes activities).
SD66_M    = np.array([0.024, 0.048, 0.096, 0.192, 0.384, 0.768,
                      1.536, 3.072, 6.144])
SD66_LOGA = np.array([-1.683, -1.403, -1.126, -0.848, -0.573, -0.292,
                      0.0036, 0.3436, 0.7897])
M_NA, M_CL = 0.0229898, 0.0354530        # kg/mol (frame transform)


def tna_expt_hittorf(mol):
    loga = np.interp(np.log(np.asarray(mol, float)), np.log(SD66_M), SD66_LOGA)
    return 0.3720 - 0.0118 * loga


def style(ax, xlabel=None, ylabel=None, title=None, fs=10):
    ax.grid(True, color=GRID, linewidth=0.7, alpha=0.9)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=INK2, labelsize=fs - 1)
    if xlabel:
        ax.set_xlabel(xlabel, color=INK, fontsize=fs)
    if ylabel:
        ax.set_ylabel(ylabel, color=INK, fontsize=fs)
    if title:
        ax.set_title(title, color=INK, fontsize=fs + 1, loc="center", pad=8)


def savefig(fig, out, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"),
                    dpi=400 if ext == "png" else None,
                    facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {name}.png/.pdf")


def load_vacf_D(plat_min=3.0, plat_max=5.0):
    """(L, sp) -> list of plateau D per seed; plus one averaged VACF/D(t) set."""
    D = {}
    curves = {}   # sp -> list of (t, vacf_norm, Dt) from the largest box
    for f in sorted(glob.glob(os.path.join(ROOT, "vacf_replica", "L*_s*", "*.vacf"))):
        L = float(re.search(r"L([\d.]+)_s", f).group(1))
        dat = np.loadtxt(f, comments="#")
        with open(f) as fh:
            names = fh.readline().lstrip("# ").split()
        t = dat[:, names.index("t_ps")]
        m = (t >= plat_min) & (t <= plat_max)
        for sp in ("O", "Na", "Cl"):
            if f"D_{sp}" not in names:
                continue
            D.setdefault((L, sp), []).append(dat[m, names.index(f"D_{sp}")].mean())
            if L == 43.47:
                curves.setdefault(sp, []).append(
                    (t, dat[:, names.index(f"vacf_{sp}")], dat[:, names.index(f"D_{sp}")]))
    return D, curves


def fig_msd_vs_vacf(out, rep, vacf_D, curves):
    Ls = list(rep["Ls"])
    species = list(rep["species"])
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.0, 4.0))

    # (a) seed-averaged VACF (normalized) + running integral, largest box
    for sp in ("O", "Na", "Cl"):
        t = curves[sp][0][0]
        v = np.mean([c[1] for c in curves[sp]], axis=0)
        a.plot(t, v, color=SPECIES_COLOR[sp], lw=1.8)
    a.axhline(0, color=MUTED, lw=0.8)
    a.set_xlim(0, 2.0)
    style(a, "t (ps)", r"$\langle v(0)\cdot v(t)\rangle$ / $\langle v^2\rangle$",
          "(a)", fs=13)

    # (b) parity: D_VACF vs D_MSD, all boxes x species
    lims = [0.06, 0.23]
    b.plot(lims, lims, color=MUTED, lw=1.0, ls="--", zorder=1)
    b.fill_between(lims, [l * 0.9 for l in lims], [l * 1.1 for l in lims],
                   color=GRID, alpha=0.5, zorder=0, label=r"$\pm$10%")
    markers = {24.84: "o", 31.05: "s", 37.26: "^", 43.47: "D"}
    for j, L in enumerate(Ls):
        for i, sp in enumerate(species):
            v = np.array(vacf_D[(L, sp)])
            b.errorbar(rep["D_mean"][i, j], v.mean(),
                       xerr=rep["D_sem"][i, j],
                       yerr=v.std(ddof=1) / np.sqrt(len(v)),
                       fmt=markers[L], ms=6, color=SPECIES_COLOR[sp],
                       mfc="white", mew=1.6, elinewidth=1.0, capsize=2, zorder=3)
    b.set_xlim(*lims)
    b.set_ylim(*lims)
    hs = [plt.Line2D([], [], marker=markers[L], ls="", mfc="white", mew=1.6,
                     color=INK2, ms=6, label=f"L = {L:.1f} Å") for L in Ls]
    b.legend(handles=hs, fontsize=10.5, frameon=False, loc="upper left",
             labelcolor=INK2)
    style(b, r"$D_{E}$ from MSD slope (Å$^2$/ps)",
          r"$D_{GK}$ from VACF integral (Å$^2$/ps)", "(b)", fs=13)
    # one species legend above the figure instead of per-panel annotations
    sp_hs = [plt.Line2D([], [], color=SPECIES_COLOR[sp], lw=2.5, label=sp)
             for sp in ("O", "Na", "Cl")]
    fig.legend(handles=sp_hs, loc="upper center", bbox_to_anchor=(0.5, 1.08),
               ncol=3, frameon=False, fontsize=13, columnspacing=2.0,
               handlelength=1.8, labelcolor=INK2)
    fig.tight_layout(w_pad=3)
    savefig(fig, out, "fig1_msd_vs_vacf")


def fig_yeh_hummer(out, rep, vacf_D):
    Ls = np.array(rep["Ls"])
    species = list(rep["species"])
    invL = 1.0 / Ls
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    x = np.linspace(0, invL.max() * 1.06, 50)
    # O: free fit (sets the YH slope / viscosity); ions: shared-slope YH fits —
    # a free per-ion slope is noise (Cl even comes out positive).
    fits = {"O": rep["free_fits"][0][:2],
            "Na": rep["shared_fits"][0][:2], "Cl": rep["shared_fits"][1][:2]}
    label_y = {"O": 0.0045, "Na": -0.013, "Cl": 0.0045}
    for i, sp in enumerate(species):
        c = SPECIES_COLOR[sp]
        m, D0 = fits[sp]
        ax.plot(x, D0 + m * x, color=c, lw=1.6, ls="--", zorder=2)
        ax.errorbar(invL, rep["D_mean"][i], yerr=rep["D_sem"][i], fmt="o", ms=6,
                    color=c, mfc="white", mew=1.6, elinewidth=1.0, capsize=2, zorder=3)
        vm = [np.mean(vacf_D[(L, sp)]) for L in Ls]
        vs = [np.std(vacf_D[(L, sp)], ddof=1) / np.sqrt(len(vacf_D[(L, sp)])) for L in Ls]
        ax.errorbar(invL + 0.0006, vm, yerr=vs, fmt="x", ms=6, color=c,
                    elinewidth=0.9, capsize=2, alpha=0.75, zorder=3)
        # per-box YH-corrected values D_PBC + |slope|/L: flat by construction
        # where the raw points are linear, so they replot the fit residuals
        ax.errorbar(invL - 0.0006, rep["D_mean"][i] - m * invL,
                    yerr=rep["D_sem"][i], fmt="D", ms=4.5, color=c,
                    elinewidth=0.9, capsize=2, alpha=0.85, zorder=3)
        ax.axhline(D0, color=c, lw=0.9, ls=":", alpha=0.55, zorder=1)
        ax.errorbar([0], [D0], fmt="*", ms=11, color=c, zorder=4)
        ax.annotate(f"{sp}:  $D_0$ = {D0:.3f}", (0.0015, D0 + label_y[sp]),
                    color=c, fontsize=9, fontweight="bold")
    hs = [plt.Line2D([], [], marker="o", ls="--", mfc="white", mew=1.6, color=INK2,
                     ms=6, label="MSD"),
          plt.Line2D([], [], marker="x", ls="", color=INK2, ms=6, label="VACF"),
          plt.Line2D([], [], marker="D", ls=":", color=INK2, ms=4.5,
                     label="YH-corrected"),
          plt.Line2D([], [], marker="*", ls="", color=INK2, ms=10,
                     label=r"$L\to\infty$ extrapolation")]
    # frameless figure-level legend above the axes (thesis style, cf.
    # plot_omp_headtohead.py); series names only, methodology in the caption
    fig.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, 1.03),
               ncol=4, frameon=False, fontsize=8.5, columnspacing=1.4,
               handlelength=2.4, handletextpad=0.6, labelcolor=INK2)
    ax.set_xlim(0, x.max())
    style(ax, r"$1/L$ (Å$^{-1}$)", r"$D_{\mathrm{PBC}}$ (Å$^2$/ps)")
    savefig(fig, out, "fig2_yeh_hummer")


def pooled(runs_mol, y):
    mols = sorted(set(runs_mol))
    mean = np.array([y[runs_mol == m].mean() for m in mols])
    sem = np.array([y[runs_mol == m].std(ddof=1) / np.sqrt((runs_mol == m).sum())
                    for m in mols])
    return np.array(mols), mean, sem


def fig_kappa(out, kap):
    mol = kap["runs_mol"]
    c = np.array([MOL_TO_M[m] for m in sorted(set(mol))])
    _, kOns, kOns_s = pooled(mol, kap["runs_kOns"])
    _, kNE, kNE_s = pooled(mol, kap["runs_kNE"])
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    ax.errorbar(c, kNE, yerr=kNE_s, fmt="^--", ms=6, color=C_O, mfc="white",
                mew=1.5, lw=1.4, elinewidth=1.0, capsize=2,
                label="Nernst–Einstein")
    ax.errorbar(c, kOns, yerr=kOns_s, fmt="o-", ms=6.5, color=C_NA, mfc="white",
                mew=1.6, lw=1.8, elinewidth=1.0, capsize=2,
                label=r"Onsager, $z=\pm1$")
    ax.errorbar(c, kOns * 0.7225, yerr=kOns_s * 0.7225, fmt="s-", ms=6, color=C_CL,
                mfc="white", mew=1.6, lw=1.4, elinewidth=1.0, capsize=2,
                label=r"Onsager, $z=\pm0.85$")
    ce = [EXPT_C[m] for m in EXPT_KAPPA]
    ax.plot(ce, list(EXPT_KAPPA.values()), marker="o", ms=5, ls="", color=SLATE3,
            label="Experiment", zorder=5)
    print("  fig3 expt anchors (Chambers-Stokes 1956): " +
          "  ".join(f"m={m:g}: {k:.2f} S/m @ c={EXPT_C[m]:.3f}"
                    for m, k in EXPT_KAPPA.items()))
    # handles come back lines-first: Experiment, NE, Onsager z=1, z=0.85 —
    # column-major that puts Experiment/NE left, the two Onsager right
    h, l = ax.get_legend_handles_labels()
    fig.legend(h, l,
               loc="upper center", bbox_to_anchor=(0.5, 1.06), ncol=2,
               frameon=False, fontsize=8.5, columnspacing=1.6,
               handlelength=2.4, handletextpad=0.6, labelcolor=INK2)
    ax.set_xlim(0, 3.9)
    ax.set_ylim(0, 24)
    style(ax, "c (mol/L)", r"$\kappa$ (S/m)")
    savefig(fig, out, "fig3_kappa_vs_c")


def fig_kappa_tna(out, kap):
    """fig35: conductivity and transport number against experiment in one
    two-panel figure (report request 2026-08-12, replaces the separate
    fig3/fig5 in the validation section)."""
    mol = kap["runs_mol"]
    mols = np.array(sorted(set(mol)))
    c = np.array([MOL_TO_M[m] for m in mols])
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.8, 4.4))

    _, kOns, kOns_s = pooled(mol, kap["runs_kOns"])
    _, kNE, kNE_s = pooled(mol, kap["runs_kNE"])
    a.errorbar(c, kNE, yerr=kNE_s, fmt="^--", ms=6, color=C_O, mfc="white",
               mew=1.5, lw=1.4, elinewidth=1.0, capsize=2,
               label="Nernst-Einstein")
    a.errorbar(c, kOns, yerr=kOns_s, fmt="o-", ms=6.5, color=C_NA,
               mfc="white", mew=1.6, lw=1.8, elinewidth=1.0, capsize=2,
               label=r"Onsager, $z=\pm1$")
    a.errorbar(c, kOns * 0.7225, yerr=kOns_s * 0.7225, fmt="s-", ms=6,
               color=C_CL, mfc="white", mew=1.6, lw=1.4, elinewidth=1.0,
               capsize=2, label=r"Onsager, $z=\pm0.85$")
    a.plot([EXPT_C[m] for m in EXPT_KAPPA], list(EXPT_KAPPA.values()),
           marker="o", ms=5, ls="", color=SLATE3, label="Experiment",
           zorder=5)
    a.set_xlim(0, 3.9)
    a.set_ylim(0, 24)
    a.legend(frameon=False, fontsize=10, labelcolor=INK2, ncol=2,
             loc="lower center", bbox_to_anchor=(0.5, 1.0),
             columnspacing=1.2, handlelength=1.9, handletextpad=0.5)
    style(a, "c (mol/L)", r"$\kappa$ (S/m)", fs=13)
    a.set_title("(a)", color=INK, fontsize=14, pad=52)

    _, tNa, tNa_s = pooled(mol, kap["runs_tNa"])
    sNN, sCC, sNC = kap["runs_sNaNa"], kap["runs_sClCl"], kap["runs_sNaCl"]
    s_sum = sNN + sCC - 2.0 * sNC
    tH_runs = (kap["runs_tNa"]
               + mol * (M_NA * (sNN - sNC) - M_CL * (sCC - sNC)) / s_sum)
    _, tH, tH_s = pooled(mol, tH_runs)
    b.errorbar(c, tNa, yerr=tNa_s, fmt="o-", ms=6.5, color=CAM_BLUE,
               mfc="white", mew=1.6, lw=1.8, elinewidth=1.0, capsize=2,
               label="Barycentric")
    b.errorbar(c, tH, yerr=tH_s, fmt="D--", ms=5.5, color=SLATE3,
               mfc="white", mew=1.3, lw=1.3, elinewidth=1.0, capsize=2,
               label="Hittorf frame")
    b.plot([EXPT_C[m] for m in mols], tna_expt_hittorf(mols), marker="o",
           ms=5, ls="", color=SLATE3, label="Experiment", zorder=5)
    b.set_xlim(0, 3.9)
    b.set_ylim(0.28, 0.56)
    b.legend(frameon=False, fontsize=10, labelcolor=INK2, ncol=2,
             loc="lower center", bbox_to_anchor=(0.5, 1.0),
             columnspacing=1.2, handlelength=1.9, handletextpad=0.5)
    style(b, "c (mol/L)", r"$t_{\mathrm{Na}^+}$", fs=13)
    b.set_title("(b)", color=INK, fontsize=14, pad=52)
    fig.tight_layout(w_pad=3)
    savefig(fig, out, "fig35_kappa_tna")


def fig_ne_deviation(out, kap):
    mol = kap["runs_mol"]
    c = np.array([MOL_TO_M[m] for m in sorted(set(mol))])
    ratio = 1.0 - kap["runs_kOns"] / kap["runs_kNE"]
    _, dNE, dNE_s = pooled(mol, ratio)
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.errorbar(c, dNE, yerr=dNE_s, fmt="o-", ms=6.5, color=CAM_BLUE,
                mfc="white", mew=1.6, lw=1.8, elinewidth=1.0, capsize=2)
    ax.set_ylim(-0.05, 0.32)
    ax.set_xlim(0, 3.9)
    style(ax, "c (mol/L)", r"$\Delta_{\mathrm{NE}}$")
    savefig(fig, out, "fig4_ne_deviation")


def fig_tna(out, kap):
    mol = kap["runs_mol"]
    mols = np.array(sorted(set(mol)))
    c = np.array([MOL_TO_M[m] for m in mols])
    _, tNa, tNa_s = pooled(mol, kap["runs_tNa"])
    # barycentric -> Hittorf (solvent frame), per run. The solvent Onsager
    # row follows from the mass constraint sum_i M_i L^ij = 0 (Fong 2020
    # Eq. 130), so with collective slopes s_ij (conductivity prefactors
    # cancel) and molality in mol/kg, masses in kg/mol:
    #   t^H = t^b + m*[M_Na(s++ - s+-) - M_Cl(s-- - s+-)] / s_sum
    # cf. Shao et al. 2022 (JACS 144, 7583) on the reference-frame gap.
    sNN, sCC, sNC = kap["runs_sNaNa"], kap["runs_sClCl"], kap["runs_sNaCl"]
    s_sum = sNN + sCC - 2.0 * sNC
    tH_runs = (kap["runs_tNa"]
               + mol * (M_NA * (sNN - sNC) - M_CL * (sCC - sNC)) / s_sum)
    _, tH, tH_s = pooled(mol, tH_runs)

    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    ax.errorbar(c, tNa, yerr=tNa_s, fmt="o-", ms=6.5, color=CAM_BLUE,
                mfc="white", mew=1.6, lw=1.8, elinewidth=1.0, capsize=2,
                label="Barycentric")
    ax.errorbar(c, tH, yerr=tH_s, fmt="D--", ms=5.5, color=SLATE3,
                mfc="white", mew=1.3, lw=1.3, elinewidth=1.0, capsize=2,
                label="Hittorf frame")
    ce = [EXPT_C[m] for m in mols]
    ax.plot(ce, tna_expt_hittorf(mols), marker="o", ms=5, ls="", color=SLATE3,
            label="Experiment", zorder=5)
    ax.set_ylim(0.28, 0.56)
    ax.set_xlim(0, 3.9)
    style(ax, "c (mol/L)", r"$t_{\mathrm{Na}^+}$")
    # handles come back lines-first: Experiment, then the two errorbars
    h, l = ax.get_legend_handles_labels()
    order = [1, 2, 0]
    fig.legend([h[i] for i in order], [l[i] for i in order],
               loc="upper center", bbox_to_anchor=(0.5, 1.03), ncol=3,
               frameon=False, fontsize=9, columnspacing=2.0,
               handlelength=2.4, handletextpad=0.6, labelcolor=INK2)
    savefig(fig, out, "fig5_tNa")
    print("  fig5 t_Na (pooled):")
    for m_, tb, tbs, th, ths in zip(mols, tNa, tNa_s, tH, tH_s):
        print(f"    m={m_:4.2f}  bary {tb:.3f}±{tbs:.3f}  "
              f"Hittorf {th:.3f}±{ths:.3f}  expt(SD66) "
              f"{tna_expt_hittorf(m_):.3f}")


def fig_onsager_decomp(out, kap):
    """kappa = kNaNa + kClCl - 2 kNaCl (z=1), and the distinct-correlation
    origin of the Nernst-Einstein deviation. Pooled over both boxes (n=10/m)."""
    E2 = (1.602176634e-19) ** 2
    KBT = 1.380649e-23 * 298.15
    # per-run prefactor: S/m per unit collective-MSD slope (A^2/ps)
    pref = E2 / (KBT * kap["runs_V"] * 1.0e-30) * 1.0e-8 / 6.0
    kNaNa = pref * kap["runs_sNaNa"]
    kClCl = pref * kap["runs_sClCl"]
    kNaCl = pref * kap["runs_sNaCl"]
    dNaNa = kNaNa - pref * 6.0 * kap["runs_nNa"] * kap["runs_DNa"]
    dClCl = kClCl - pref * 6.0 * kap["runs_nCl"] * kap["runs_DCl"]
    kNE = kap["runs_kNE"]

    mol = kap["runs_mol"]
    c = np.array([MOL_TO_M[m] for m in sorted(set(mol))])

    # fig6: conductivity contributions in absolute units
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    for y, color, fmt, lw, label in (
            (kap["runs_kOns"], SLATE3, "o-", 1.8, r"total $\kappa_{\mathrm{Ons}}$"),
            (kNaNa, C_NA, "s-", 1.4, r"$\kappa_{\mathrm{NaNa}}$"),
            (kClCl, C_CL, "^-", 1.4, r"$\kappa_{\mathrm{ClCl}}$"),
            (-2 * kNaCl, C_O, "D-", 1.4, r"$-2\,\kappa_{\mathrm{NaCl}}$")):
        _, ym, ys = pooled(mol, y)
        ax.errorbar(c, ym, yerr=ys, fmt=fmt, ms=5.5, color=color, mfc="white",
                    mew=1.5, lw=lw, elinewidth=1.0, capsize=2, label=label)
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_xlim(0, 3.9)
    h, l = ax.get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 1.06), ncol=2,
               frameon=False, fontsize=9, columnspacing=2.0,
               handlelength=2.4, handletextpad=0.6, labelcolor=INK2)
    style(ax, "c (mol/L)", r"$\kappa$ contribution (S/m)")
    savefig(fig, out, "fig6_kappa_decomposition")

    # fig7: contributions to the NE deviation: Delta = 1 - kOns/kNE
    #       = [-dNaNa - dClCl + 2 dNaCl] / kNE   (per run, then pooled)
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    parts = ((-dNaNa / kNE, C_NA, "s-", "Na–Na distinct"),
             (-dClCl / kNE, C_CL, "^-", "Cl–Cl distinct"),
             (2 * kNaCl / kNE, C_O, "D-", "Na–Cl distinct"))
    for y, color, fmt, label in parts:
        _, ym, ys = pooled(mol, y)
        ax.errorbar(c, ym, yerr=ys, fmt=fmt, ms=5.5, color=color, mfc="white",
                    mew=1.5, lw=1.4, elinewidth=1.0, capsize=2, label=label)
    _, tot, tot_s = pooled(mol, 1.0 - kap["runs_kOns"] / kNE)
    ax.errorbar(c, tot, yerr=tot_s, fmt="o-", ms=5.5, color=SLATE3, mfc="white",
                mew=1.5, lw=1.8, elinewidth=1.0, capsize=2,
                label=r"total $\Delta_{\mathrm{NE}}$")
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_xlim(0, 3.9)
    ax.set_ylim(-0.19, 0.31)
    h, l = ax.get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 1.06), ncol=2,
               frameon=False, fontsize=9, columnspacing=2.0,
               handlelength=2.4, handletextpad=0.6, labelcolor=INK2)
    style(ax, "c (mol/L)", r"contribution to $\Delta_{\mathrm{NE}}$")
    savefig(fig, out, "fig7_ne_origin")

    # fig67: both decompositions as one two-panel figure (report request
    # 2026-08-12, replaces the separate fig6/fig7 in section 9.1)
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.8, 4.6))
    for y, color, fmt, lw, label in (
            (kap["runs_kOns"], SLATE3, "o-", 1.8, r"total $\kappa_{\mathrm{Ons}}$"),
            (kNaNa, C_NA, "s-", 1.4, r"$\kappa_{\mathrm{NaNa}}$"),
            (kClCl, C_CL, "^-", 1.4, r"$\kappa_{\mathrm{ClCl}}$"),
            (-2 * kNaCl, C_O, "D-", 1.4, r"$-2\,\kappa_{\mathrm{NaCl}}$")):
        _, ym, ys = pooled(mol, y)
        a.errorbar(c, ym, yerr=ys, fmt=fmt, ms=5.5, color=color, mfc="white",
                   mew=1.5, lw=lw, elinewidth=1.0, capsize=2, label=label)
    a.axhline(0, color=MUTED, lw=0.8)
    a.set_xlim(0, 3.9)
    a.legend(frameon=False, fontsize=10, labelcolor=INK2, ncol=2,
             loc="lower center", bbox_to_anchor=(0.5, 1.0),
             columnspacing=1.2, handlelength=1.9, handletextpad=0.5)
    style(a, "c (mol/L)", r"$\kappa$ contribution (S/m)", fs=13)
    a.set_title("(a)", color=INK, fontsize=14, pad=52)
    for y, color, fmt, label in ((-dNaNa / kNE, C_NA, "s-", "Na-Na distinct"),
                                 (-dClCl / kNE, C_CL, "^-", "Cl-Cl distinct"),
                                 (2 * kNaCl / kNE, C_O, "D-", "Na-Cl distinct")):
        _, ym, ys = pooled(mol, y)
        b.errorbar(c, ym, yerr=ys, fmt=fmt, ms=5.5, color=color, mfc="white",
                   mew=1.5, lw=1.4, elinewidth=1.0, capsize=2, label=label)
    _, tot, tot_s = pooled(mol, 1.0 - kap["runs_kOns"] / kNE)
    b.errorbar(c, tot, yerr=tot_s, fmt="o-", ms=5.5, color=SLATE3,
               mfc="white", mew=1.5, lw=1.8, elinewidth=1.0, capsize=2,
               label=r"total $\Delta_{\mathrm{NE}}$")
    b.axhline(0, color=MUTED, lw=0.8)
    b.set_xlim(0, 3.9)
    b.set_ylim(-0.19, 0.31)
    b.legend(frameon=False, fontsize=10, labelcolor=INK2, ncol=2,
             loc="lower center", bbox_to_anchor=(0.5, 1.0),
             columnspacing=1.2, handlelength=1.9, handletextpad=0.5)
    style(b, "c (mol/L)", r"contribution to $\Delta_{\mathrm{NE}}$", fs=13)
    b.set_title("(b)", color=INK, fontsize=14, pad=52)
    fig.tight_layout(w_pad=3)
    savefig(fig, out, "fig67_decomposition")

    print("  pooled decomposition (z=1, S/m):")
    for name, y in (("kNaNa", kNaNa), ("kClCl", kClCl), ("-2kNaCl", -2 * kNaCl),
                    ("dNaNa", dNaNa), ("dClCl", dClCl), ("dNaCl", kNaCl)):
        _, ym, ys = pooled(mol, y)
        print(f"    {name:>8}: " + "  ".join(
            f"{m:4.2f}m {v:+6.3f}±{s:5.3f}" for m, v, s in zip(sorted(set(mol)), ym, ys)))


# fig5 (standalone PMF-vs-benchmark) retired 2026-07-30: redundant with
# fig7 panel (a), which carries the same comparison plus the small MP2 cell.


ANCHOR = os.path.join(_REPO, "results", "nacl_mp2_anchor")


def _madrid_pmf():
    """Seed-averaged Madrid w(r) at 1 m (same loader as fig_pmf)."""
    gs, r = [], None
    for f in sorted(glob.glob(os.path.join(ROOT, "conductivity", "m1.0", "L*_s*", "*.rdf"))):
        dat = np.loadtxt(f, skiprows=4)
        if r is None:
            r = dat[:, 1]
        elif not np.allclose(dat[:, 1], r):
            dat = np.column_stack([dat[:, 0], r,
                                   np.interp(r, dat[:, 1], dat[:, 2]), dat[:, 3]])
        gs.append(dat[:, 2])
    g = np.mean(gs, axis=0)
    w = np.full_like(g, np.nan)
    m = g > 0.02
    w[m] = -KT_KCAL * np.log(g[m])
    w -= np.nanmean(w[(r > 9.0) & (r < r.max())])
    return r, w


def _extrema(r, w):
    def ext(lo, hi, kind):
        m = (r >= lo) & (r <= hi) & np.isfinite(w)
        i = (np.argmin if kind == "min" else np.argmax)(w[m])
        return r[m][i], w[m][i]
    return ext(2.4, 3.2, "min"), ext(3.2, 4.2, "max"), ext(4.2, 5.6, "min")


def fig_nnp_anchor(out, rep):
    """fig7: MP2 C-NNP anchor vs Madrid-2019 — pairing PMF + D ratios."""
    r_mad, w_mad = _madrid_pmf()
    pmf = {c: np.loadtxt(os.path.join(ANCHOR, f"rdf_nacl_{c}.csv"),
                         delimiter=",", skiprows=1).T for c in ("cube2", "cube3")}

    # diffusion ratios D_ion / D_water: MP2 per box (uncorrected), Madrid D0
    wat = np.genfromtxt(os.path.join(ANCHOR, "diffusion_summary.csv"),
                        delimiter=",", names=True, dtype=None, encoding=None)
    ion = np.genfromtxt(os.path.join(ANCHOR, "ion_diffusion_summary.csv"),
                        delimiter=",", names=True, dtype=None, encoding=None)

    def pool(rows, key):
        d = np.array([x["D_1e9_m2_s"] for x in rows])
        return d.mean(), d.std(ddof=1) / np.sqrt(len(d))

    mp2 = {}
    for cell in ("cube2", "cube3"):
        dO, sO = pool([x for x in wat if x["cell"] == cell], None)
        mp2[cell] = {"O": (dO, sO)}
        for sp in ("Na", "Cl"):
            mp2[cell][sp] = pool([x for x in ion if x["cell"] == cell
                                  and x["species"] == sp], None)

    def ratio(num, den):
        v = num[0] / den[0]
        return v, v * np.hypot(num[1] / num[0], den[1] / den[0])

    # Madrid D0 (YH shared-slope fit): free_fits row 0 = O, shared_fits = Na, Cl
    mad_O = (rep["free_fits"][0][1], rep["free_fits"][0][2])
    mad = {"Na": ratio((rep["shared_fits"][0][1], rep["shared_fits"][0][2]), mad_O),
           "Cl": ratio((rep["shared_fits"][1][1], rep["shared_fits"][1][2]), mad_O)}
    # expt infinite-dilution: D0_ion from limiting conductivities (Robinson &
    # Stokes: lambda0 50.10 / 76.35 S cm2/mol -> 1.334 / 2.032 e-9 m2/s via
    # Nernst), water self-D 2.299e-9 (Holz 2000)
    expt = {"Na": 1.334 / 2.299, "Cl": 2.032 / 2.299}

    fig, (a, b) = plt.subplots(1, 2, figsize=(9.8, 4.1))

    a.axhline(0, color=MUTED, lw=0.8)
    a.plot(r_mad, w_mad, color=C_NA, lw=2.0, label="Madrid-2019")
    r3, g3, w3 = pmf["cube3"]
    r2, g2, w2 = pmf["cube2"]
    a.plot(r3, w3, color=SLATE3, lw=2.0, label="MP2 C-NNP, $L$=37.3 Å")
    a.plot(r2, w2, color=MUTED, lw=1.1, ls="--",
           label="MP2 C-NNP, $L$=24.8 Å")
    a.set_xlim(2.2, 8.0)
    a.set_ylim(-1.3, 1.5)
    a.legend(fontsize=8.5, frameon=False, loc="lower right", labelcolor=INK2)
    style(a, r"$r_{\mathrm{Na-Cl}}$ (Å)", r"$w(r)$ (kcal/mol)", "(a)")

    xs = {"Na": 0.0, "Cl": 1.0}
    for sp in ("Na", "Cl"):
        x = xs[sp]
        v, e = mad[sp]
        b.errorbar(x - 0.18, v, yerr=e, fmt="s", ms=6, color=C_NA,
                   capsize=3, label="Madrid-2019" if sp == "Na" else None)
        for off, cell, mfc in ((-0.02, "cube2", "white"), (0.14, "cube3", SLATE3)):
            v, e = ratio(mp2[cell][sp], mp2[cell]["O"])
            b.errorbar(x + off, v, yerr=e, fmt="o", ms=6, color=SLATE3, mfc=mfc,
                       capsize=3,
                       label=(f"MP2 C-NNP, $L$={'24.8' if cell == 'cube2' else '37.3'} Å"
                              if sp == "Na" else None))
        b.plot(x + 0.30, expt[sp], marker="_", ms=14, mew=2.2, color=INK2, ls="",
               label="Experiment" if sp == "Na" else None)
    b.set_xticks([0.06, 1.06])
    b.set_xticklabels([r"$D_{\mathrm{Na}}/D_{\mathrm{w}}$",
                       r"$D_{\mathrm{Cl}}/D_{\mathrm{w}}$"])
    b.set_xlim(-0.5, 1.6)
    b.set_ylim(0.35, 1.0)
    b.legend(fontsize=8, frameon=False, loc="upper left", labelcolor=INK2)
    style(b, "", r"$D_{\mathrm{ion}}\,/\,D_{\mathrm{water}}$", "(b)")
    savefig(fig, out, "fig8_nnp_anchor")

    for tag, rr, ww in (("Madrid", r_mad, w_mad), ("MP2 cube3", r3, w3),
                        ("MP2 cube2", r2, w2)):
        (rc, wc), (rt, wt), (rs, ws) = _extrema(rr, ww)
        print(f"  PMF {tag}: CIP {wc:+.3f} @ {rc:.2f} | TS {wt:+.3f} @ {rt:.2f} | "
              f"SSIP {ws:+.3f} @ {rs:.2f} | dW(CIP-SSIP) {wc - ws:+.3f} kcal/mol")

    print("  fig7 table: D (1e-9 m2/s) and ratios")
    for cell in ("cube2", "cube3"):
        row = mp2[cell]
        print(f"    MP2 {cell}: O {row['O'][0]:.3f}±{row['O'][1]:.3f}  "
              f"Na {row['Na'][0]:.3f}±{row['Na'][1]:.3f}  "
              f"Cl {row['Cl'][0]:.3f}±{row['Cl'][1]:.3f}  "
              f"Na/O {ratio(row['Na'], row['O'])[0]:.3f}  "
              f"Cl/O {ratio(row['Cl'], row['O'])[0]:.3f}")
    print(f"    Madrid D0 ratios: Na/O {mad['Na'][0]:.3f}±{mad['Na'][1]:.3f}  "
          f"Cl/O {mad['Cl'][0]:.3f}±{mad['Cl'][1]:.3f}   "
          f"expt(inf dil): Na 0.580  Cl 0.884")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(ROOT, "figures"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    plt.rcParams.update({"font.size": 9.5, "axes.edgecolor": MUTED,
                         "text.color": INK, "axes.labelcolor": INK})

    rep = np.load(os.path.join(ROOT, "replicas", "analysis", "replica_D.npz"),
                  allow_pickle=True)
    kap = np.load(os.path.join(ROOT, "conductivity", "analysis", "kappa_vs_c.npz"),
                  allow_pickle=True)
    vacf_D, curves = load_vacf_D()

    fig_msd_vs_vacf(args.out, rep, vacf_D, curves)
    fig_yeh_hummer(args.out, rep, vacf_D)
    fig_kappa(args.out, kap)
    fig_ne_deviation(args.out, kap)
    fig_tna(args.out, kap)
    fig_kappa_tna(args.out, kap)
    fig_onsager_decomp(args.out, kap)
    fig_nnp_anchor(args.out, rep)

    # console check: pooled kappa table vs headline numbers
    mol = kap["runs_mol"]
    mols, kOns, kOns_s = pooled(mol, kap["runs_kOns"])
    for m, k, s in zip(mols, kOns, kOns_s):
        print(f"  m={m:4.2f}  kappa_Ons(z=1) = {k:5.2f} ± {s:4.2f} S/m")


if __name__ == "__main__":
    main()
