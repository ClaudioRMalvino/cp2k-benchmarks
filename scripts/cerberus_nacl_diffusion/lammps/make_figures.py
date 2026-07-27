#!/usr/bin/env python3
"""Paper/supervisor figures for the Madrid-2019 NaCl(aq) transport study.

Reads only finished analysis products + raw .rdf/.vacf files:
  replicas/analysis/replica_D.npz      (MSD self-diffusion, Yeh-Hummer)
  vacf_replica/L*_s*/*.vacf            (Green-Kubo D, corrected campaign)
  conductivity/analysis/kappa_vs_c.npz (Onsager/NE conductivity, 50 runs)
  conductivity/m*/L*_s*/*.rdf          (Na-Cl g(r) -> PMF)

Outputs PNG (400 dpi, slides) + PDF (vector, paper) per figure into --out.

Experimental kappa anchors (298.15 K): m = 1, 2, 4 -> 8.48, 14.49, 22.04 S/m
as compiled in Blazquez et al., J. Phys. Chem. B 2023 ("Two Surfaces, One
Property"); dilute molalities left unpinned until we fix a citable source.
"""
import argparse
import glob
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/data/cerberus1/crm98/nacl_diffusion/lammps_madrid"
KT_KCAL = 0.0019872041 * 298.15          # k_B T in kcal/mol at 298.15 K

# validated categorical palette (dataviz reference, light mode, fixed order)
C_NA, C_CL, C_O = "#2a78d6", "#008300", "#e87ba4"   # slots 1-3
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
SPECIES_COLOR = {"Na": C_NA, "Cl": C_CL, "O": C_O}

MOL_TO_M = {0.25: 0.231, 0.5: 0.491, 1.0: 0.960, 2.0: 1.915, 4.0: 3.666}
EXPT_KAPPA = {1.0: 8.48, 2.0: 14.49, 4.0: 22.04}    # S/m, Blazquez 2023 compilation


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
        ax.set_title(title, color=INK, fontsize=11, loc="left", pad=8)


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
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.0, 3.8))

    # (a) seed-averaged VACF (normalized) + running integral, largest box
    for sp in ("O", "Na", "Cl"):
        t = curves[sp][0][0]
        v = np.mean([c[1] for c in curves[sp]], axis=0)
        a.plot(t, v, color=SPECIES_COLOR[sp], lw=1.8)
        a.annotate(sp, (t[np.argmin(np.abs(t - {"O": 0.12, "Na": 0.32, "Cl": 0.55}[sp]))],
                        v[np.argmin(np.abs(t - {"O": 0.12, "Na": 0.32, "Cl": 0.55}[sp]))]),
                   textcoords="offset points", xytext=(6, 4),
                   color=SPECIES_COLOR[sp], fontsize=9, fontweight="bold")
    a.axhline(0, color=MUTED, lw=0.8)
    a.set_xlim(0, 2.0)
    style(a, "t (ps)", r"$\langle v(0)\cdot v(t)\rangle$ / $\langle v^2\rangle$",
          "a  Velocity autocorrelation (L = 43.5 Å, seed avg.)")

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
    for sp, xy in (("O", (0.192, 0.207)), ("Na", (0.089, 0.077)), ("Cl", (0.135, 0.122))):
        b.annotate(sp, xy, color=SPECIES_COLOR[sp], fontsize=9, fontweight="bold")
    b.set_xlim(*lims)
    b.set_ylim(*lims)
    hs = [plt.Line2D([], [], marker=markers[L], ls="", mfc="white", mew=1.6,
                     color=INK2, ms=6, label=f"L = {L:.1f} Å") for L in Ls]
    b.legend(handles=hs, fontsize=8, frameon=False, loc="upper left",
             labelcolor=INK2)
    style(b, r"$D$ from MSD slope (Å$^2$/ps)", r"$D$ from VACF integral (Å$^2$/ps)",
          "b  Einstein vs Green–Kubo (all 12 within 2σ)")
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
        ax.errorbar([0], [D0], fmt="*", ms=11, color=c, zorder=4)
        ax.annotate(f"{sp}:  $D_0$ = {D0:.3f}", (0.0015, D0 + label_y[sp]),
                    color=c, fontsize=9, fontweight="bold")
    hs = [plt.Line2D([], [], marker="o", ls="--", mfc="white", mew=1.6, color=INK2,
                     ms=6, label="MSD (NVT) + YH fit (shared slope)"),
          plt.Line2D([], [], marker="x", ls="", color=INK2, ms=6, label="VACF (NVE)"),
          plt.Line2D([], [], marker="*", ls="", color=INK2, ms=10,
                     label=r"$L\to\infty$ extrapolation")]
    ax.legend(handles=hs, fontsize=8, frameon=False, loc="lower left", labelcolor=INK2)
    ax.set_xlim(0, x.max())
    style(ax, r"$1/L$ (Å$^{-1}$)", r"$D_{\mathrm{PBC}}$ (Å$^2$/ps)",
          "Yeh–Hummer finite-size scaling, 1 m NaCl, 298.15 K")
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
                label=r"Nernst–Einstein, $z=\pm1$")
    ax.errorbar(c, kOns, yerr=kOns_s, fmt="o-", ms=6.5, color=C_NA, mfc="white",
                mew=1.6, lw=1.8, elinewidth=1.0, capsize=2,
                label=r"Onsager (Einstein–Helfand), $z=\pm1$")
    ax.errorbar(c, kOns * 0.7225, yerr=kOns_s * 0.7225, fmt="s-", ms=6, color=C_CL,
                mfc="white", mew=1.6, lw=1.4, elinewidth=1.0, capsize=2,
                label=r"Onsager, $z=\pm0.85$ (scaled)")
    ce = [MOL_TO_M[m] for m in EXPT_KAPPA]
    ax.plot(ce, list(EXPT_KAPPA.values()), marker="*", ms=13, ls="", color=INK,
            label="Experiment (298 K)", zorder=5)
    ax.legend(fontsize=8.5, frameon=False, loc="upper left", labelcolor=INK2)
    ax.set_xlim(0, 3.9)
    ax.set_ylim(0, 24)
    style(ax, "c (mol/L)", r"$\kappa$ (S/m)",
          "Electrical conductivity of NaCl(aq), Madrid-2019")
    savefig(fig, out, "fig3_kappa_vs_c")


def fig_transport(out, kap):
    mol = kap["runs_mol"]
    c = np.array([MOL_TO_M[m] for m in sorted(set(mol))])
    ratio = 1.0 - kap["runs_kOns"] / kap["runs_kNE"]
    _, dNE, dNE_s = pooled(mol, ratio)
    _, tNa, tNa_s = pooled(mol, kap["runs_tNa"])
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.0, 3.6))

    a.axhline(0, color=MUTED, lw=0.8)
    a.errorbar(c, dNE, yerr=dNE_s, fmt="o-", ms=6.5, color=C_NA, mfc="white",
               mew=1.6, lw=1.8, elinewidth=1.0, capsize=2)
    a.set_ylim(-0.05, 0.32)
    a.set_xlim(0, 3.9)
    style(a, "c (mol/L)", r"$1-\kappa_{\mathrm{Ons}}/\kappa_{\mathrm{NE}}$",
          "a  Ion-correlation reduction of κ")

    b.errorbar(c, tNa, yerr=tNa_s, fmt="o-", ms=6.5, color=C_NA, mfc="white",
               mew=1.6, lw=1.8, elinewidth=1.0, capsize=2, label="Madrid-2019 (this work)")
    b.axhline(0.39, color=INK2, lw=1.2, ls=":")
    b.annotate("expt ≈ 0.39", (3.85, 0.393), color=INK2, fontsize=8.5, ha="right")
    b.axhline(0.37, color=MUTED, lw=1.0, ls="--")
    b.annotate("Gullbrekken 2023 (sim) 0.37", (3.85, 0.352), color=MUTED,
               fontsize=8.5, ha="right")
    b.set_ylim(0.25, 0.60)
    b.set_xlim(0, 3.9)
    b.legend(fontsize=8.5, frameon=False, loc="upper left", labelcolor=INK2)
    style(b, "c (mol/L)", r"$t_{\mathrm{Na}^+}$", "b  Cation transport number")
    fig.tight_layout(w_pad=3)
    savefig(fig, out, "fig4_ne_deviation_tNa")


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
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.0, 3.6))

    # (a) conductivity contributions
    for y, color, fmt, lw, label in (
            (kap["runs_kOns"], INK, "o-", 1.8, r"total $\kappa_{\mathrm{Ons}}$"),
            (kNaNa, C_NA, "s-", 1.4, r"$\kappa_{\mathrm{NaNa}}$"),
            (kClCl, C_CL, "^-", 1.4, r"$\kappa_{\mathrm{ClCl}}$"),
            (-2 * kNaCl, C_O, "D-", 1.4, r"$-2\,\kappa_{\mathrm{NaCl}}$ (cross)")):
        _, ym, ys = pooled(mol, y)
        a.errorbar(c, ym, yerr=ys, fmt=fmt, ms=5.5, color=color, mfc="white",
                   mew=1.5, lw=lw, elinewidth=1.0, capsize=2, label=label)
    a.axhline(0, color=MUTED, lw=0.8)
    a.legend(fontsize=8.5, frameon=False, loc="upper left", labelcolor=INK2)
    a.set_xlim(0, 3.9)
    style(a, "c (mol/L)", r"$\kappa$ contribution (S/m)",
          "a  Onsager decomposition of κ  (z = ±1)")

    # (b) contributions to the NE deviation: Delta = 1 - kOns/kNE
    #     = [-dNaNa - dClCl + 2 dNaCl] / kNE   (per run, then pooled)
    parts = ((-dNaNa / kNE, C_NA, "s-", "Na–Na distinct"),
             (-dClCl / kNE, C_CL, "^-", "Cl–Cl distinct"),
             (2 * kNaCl / kNE, C_O, "D-", "Na–Cl distinct (pairing)"))
    for y, color, fmt, label in parts:
        _, ym, ys = pooled(mol, y)
        b.errorbar(c, ym, yerr=ys, fmt=fmt, ms=5.5, color=color, mfc="white",
                   mew=1.5, lw=1.4, elinewidth=1.0, capsize=2, label=label)
    _, tot, tot_s = pooled(mol, 1.0 - kap["runs_kOns"] / kNE)
    b.errorbar(c, tot, yerr=tot_s, fmt="o-", ms=5.5, color=INK, mfc="white",
               mew=1.5, lw=1.8, elinewidth=1.0, capsize=2,
               label=r"total $\Delta_{\mathrm{NE}}$")
    b.axhline(0, color=MUTED, lw=0.8)
    b.annotate("pairing term ≈ 0 within noise\n(under-paired PMF, cf. Fig. 5)",
               (3.75, -0.155), color=C_O, fontsize=8.5, ha="right")
    b.legend(fontsize=8.5, frameon=False, loc="upper left", labelcolor=INK2, ncol=2)
    b.set_xlim(0, 3.9)
    b.set_ylim(-0.19, 0.31)
    style(b, "c (mol/L)", r"contribution to $1-\kappa_{\mathrm{Ons}}/\kappa_{\mathrm{NE}}$",
          "b  Origin of the Nernst–Einstein deviation")
    fig.tight_layout(w_pad=3)
    savefig(fig, out, "fig6_onsager_decomposition")

    print("  pooled decomposition (z=1, S/m):")
    for name, y in (("kNaNa", kNaNa), ("kClCl", kClCl), ("-2kNaCl", -2 * kNaCl),
                    ("dNaNa", dNaNa), ("dClCl", dClCl), ("dNaCl", kNaCl)):
        _, ym, ys = pooled(mol, y)
        print(f"    {name:>8}: " + "  ".join(
            f"{m:4.2f}m {v:+6.3f}±{s:5.3f}" for m, v, s in zip(sorted(set(mol)), ym, ys)))


def fig_pmf(out):
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
    mask = g > 0.02
    w = np.full_like(g, np.nan)
    w[mask] = -KT_KCAL * np.log(g[mask])
    w -= np.nanmean(w[(r > 9.0) & (r < r.max())])          # zero at large r

    def extremum(lo, hi, kind):
        m = (r >= lo) & (r <= hi) & np.isfinite(w)
        i = (np.argmin if kind == "min" else np.argmax)(w[m])
        return r[m][i], w[m][i]

    r_cip, w_cip = extremum(2.5, 3.2, "min")
    r_ts, w_ts = extremum(3.2, 4.2, "max")
    r_ssip, w_ssip = extremum(4.2, 5.6, "min")

    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.plot(r, w, color=C_NA, lw=2.0, label="Madrid-2019, 1 m (this work)")
    offsets = {"CIP": (0, -26), "barrier": (0, 10), "SSIP": (12, -14)}
    for rx, wx, name in ((r_cip, w_cip, "CIP"), (r_ts, w_ts, "barrier"),
                         (r_ssip, w_ssip, "SSIP")):
        ax.plot(rx, wx, "o", ms=6, color=C_NA, mfc="white", mew=1.6)
        ax.annotate(f"{name}\n{wx:+.2f}", (rx, wx),
                    textcoords="offset points", xytext=offsets[name],
                    ha="left" if name == "SSIP" else "center",
                    color=INK2, fontsize=8.5)
    # O'Neill 2024 MP2/RPA benchmark: CIP and SSIP equistable within ~0.2 kcal/mol
    # -> gray band = where the CIP minimum would sit if CIP ~ SSIP
    ax.fill_between([r_cip - 0.3, r_cip + 0.3], w_ssip - 0.2, w_ssip + 0.2,
                    color=GRID, alpha=0.9, zorder=0)
    ax.annotate("gray band: CIP depth if CIP ≈ SSIP\n(MP2/RPA benchmark, O'Neill 2024)",
                (r_cip - 0.35, w_ssip - 0.42), color=INK2, fontsize=8.5)
    dW = w_cip - w_ssip
    ax.annotate(rf"Madrid: $\Delta W_{{\mathrm{{CIP-SSIP}}}}$ = {dW:+.2f} kcal/mol"
                "\n→ under-paired vs benchmark",
                (0.97, 0.70), xycoords="axes fraction", ha="right",
                color=INK, fontsize=9)
    ax.set_xlim(2.2, 8.0)
    ax.set_ylim(-1.0, 1.6)
    ax.legend(fontsize=8.5, frameon=False, loc="upper right", labelcolor=INK2)
    style(ax, r"$r_{\mathrm{Na-Cl}}$ (Å)", r"$w(r) = -k_BT\,\ln g(r)$ (kcal/mol)",
          "Na–Cl pairing free energy vs first-principles benchmark")
    savefig(fig, out, "fig5_pmf_vs_benchmark")
    print(f"  PMF: CIP {w_cip:+.3f} @ {r_cip:.2f} | TS {w_ts:+.3f} @ {r_ts:.2f} | "
          f"SSIP {w_ssip:+.3f} @ {r_ssip:.2f} | dW(CIP-SSIP) {dW:+.3f} kcal/mol")


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
    fig_transport(args.out, kap)
    fig_pmf(args.out)
    fig_onsager_decomp(args.out, kap)

    # console check: pooled kappa table vs headline numbers
    mol = kap["runs_mol"]
    mols, kOns, kOns_s = pooled(mol, kap["runs_kOns"])
    for m, k, s in zip(mols, kOns, kOns_s):
        print(f"  m={m:4.2f}  kappa_Ons(z=1) = {k:5.2f} ± {s:4.2f} S/m")


if __name__ == "__main__":
    main()
