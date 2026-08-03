#!/usr/bin/env python3
"""RPA concentration-series pooling + report figures (fig9-fig14).

Pools the per-segment CP2K analysis products (written on CSD3 by
analyze_onsager_cp2k.py / analyze_viscosity_cp2k.py; RPA-series files in
results/rpa_transport/, MP2-anchor files in results/mp2_transport/) and
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

All experimental anchors are pinned from primary sources (2026-08-02):
kappa Chambers-Stokes 1956, t_Na Smits-Duyvis 1966, eta Kestin-Khalifa-
Correia 1981, D_w Muller-Hertz 1996 via Blazquez 2023 Table V.

Also prints: fit-window sensitivity (10-30 vs 10-40 ps), the 1 m three-way
Madrid/MP2/RPA block, the D_w timescale check, and the timings-ledger
aggregation. Cross-level outputs (figures, transport_series_summary.csv)
go to results/transport_comparison/.
"""
import csv
import glob
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_REPO = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                      "..", ".."))
MP2T = os.path.join(_REPO, "results", "mp2_transport")      # MP2-anchor products
RPAT = os.path.join(_REPO, "results", "rpa_transport")      # RPA-series products
CMP = os.path.join(_REPO, "results", "transport_comparison")  # cross-level outputs
MADRID = os.path.join(_REPO, "results", "lammps_madrid",
                      "conductivity", "analysis", "kappa_vs_c.npz")
OUT = os.path.join(CMP, "figures")


def res_path(fname):
    """Level routing, same rule as the CSD3 analyzers: mp2* files live in
    results/mp2_transport/, everything else in results/rpa_transport/."""
    return os.path.join(MP2T if fname.startswith("mp2") else RPAT, fname)

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
FARADAY, NA = 96485.33212, 6.02214076e23
CONV = 1.0e-8 * 1e9                      # A^2/ps -> 1e-9 m^2/s

# water self-diffusion (1e-9 m^2/s, 298.15 K), PINNED 2026-08-02:
# m = 0 from Holz 2000 (PFG-NMR); m = 1/2/4 from the Muller & Hertz 1996
# NMR dataset (J. Phys. Chem. 100, 1256) as tabulated in Blazquez 2023
# Table V (J. Chem. Phys. 158, 054505) -- the same chain O'Neill's SI
# Fig. S10 cites for its experimental curve.
DW_EXP = {0.0: 2.299, 1.0: 2.17, 2.0: 2.02, 4.0: 1.71}

# eta (mPa s, 25 C, p* = 0.1 MPa), PINNED 2026-08-02 from Kestin, Khalifa
# & Correia 1981 (J. Phys. Chem. Ref. Data 10, 71; DOI 10.1063/1.555641),
# Tables 1/3/5/9, stated accuracy +-0.5%; m = 3 for reference: 1.1998.
ETA_EXP = {0.0: 0.8901, 1.0: 0.9723, 2.0: 1.0745, 4.0: 1.3504}

# GK plateau window. The pipeline default (2-10 ps, stored in the npz) is
# biased low at RPA 1 m: the running integral is still rising there, the
# per-segment late-early difference is +0.31 mPa s at 7.9 sigma. Over
# 10-18 ps the slope of the segment-mean curve is zero within segment noise
# for 2 m / 4 m / the anchor and the anchor value moves by less than its
# SEM, so 10-18 ps is the headline window for all four datasets; the 2-10 ps
# value is printed alongside as the sensitivity.
ETA_WINDOW = (10.0, 18.0)

PAIRING = {1.0: "rpa_gpu_cube3_1M", 2.0: "rpa_gpu_cube3_2m",
           4.0: "rpa_gpu_cube3_4m", "mp2": "mp2_anchor_cube3"}
MADRID_RDF = os.path.join(_REPO, "results", "lammps_madrid", "conductivity")
KT_KCAL = 0.0019872041 * 298.15
MADRID_NCIP_1M = 0.026   # report value: coordination integral of the Madrid
                         # g_NaCl at 1 m to the 3.38 A barrier (verified)

# 1 m uses the uniform 5x80 truncation, not the *_full_* set (segs 4-5 there
# run to 125/155 ps): user decision 2026-08-03, keep all three molalities
# like-for-like until the 160 ps extension suite lands.
SERIES = {1.0: ("rpa_gpu_cube3_5x80_w1040_kappa.npz", "rpa_gpu_cube3_final_eta.npz"),
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
    z = np.load(res_path(fname))
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
    # Fong-style observables. Lambda = kappa/c; electrophoretic mobilities
    # from the Hittorf (solvent-fixed) partition sigma_Na = t^H kappa, the
    # frame experimental mobilities are defined in. u_Cl is the magnitude
    # (-u for the anion, as Fong plots it).
    d["c_molar"] = int(d["nNa"]) / (NA * float(d["V"]) * 1e-27)
    d["Lambda"] = 10.0 * d["kOns"] / d["c_molar"]
    fc = FARADAY * d["c_molar"] * 1e3
    d["uNa"] = d["tH"] * d["kOns"] / fc * 1e8
    d["uCl"] = (1.0 - d["tH"]) * d["kOns"] / fc * 1e8
    return d


def load_eta(fname):
    z = np.load(res_path(fname))
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


def load_pairing(tag):
    """CSD3 pairing products: summary row + w(r) with report-window extrema."""
    row = next(csv.DictReader(open(res_path(f"{tag}_pairing.csv"))))
    rdf = np.loadtxt(res_path(f"{tag}_rdf_ions.csv"),
                     delimiter=",", skiprows=1)
    r, w = rdf[:, 0], rdf[:, 2]
    d = {k: float(v) for k, v in row.items() if k not in ("label",)}
    d["r"], d["w"], d["g"] = r, w, rdf[:, 1]
    # tab:pmf conventions: CIP min 2.4-3.2, barrier max 3.2-4.2, SSIP min 4.2-5.6
    for key, lo, hi, kind in (("CIP", 2.4, 3.2, "min"), ("TS", 3.2, 4.2, "max"),
                              ("SSIP", 4.2, 5.6, "min")):
        m = (r >= lo) & (r <= hi) & np.isfinite(w)
        i = (np.argmin if kind == "min" else np.argmax)(w[m])
        d[f"r_{key}"], d[f"w_{key}"] = r[m][i], w[m][i]
    d["dW_cip_ssip"] = d["w_CIP"] - d["w_SSIP"]
    return d


def madrid_pairing(mol):
    """Seed-averaged Madrid pairing observables at one molality, from the
    conductivity campaign's per-run LAMMPS rdf fixes (compute rdf 200 for
    Na-Cl; columns: bin, r, g, coord). w(r) extrema use the tab:pmf windows
    of load_pairing; n_CIP per seed is the LAMMPS running coordination
    number interpolated at the pooled barrier radius (identical integral to
    the CSD3 pipeline's 4 pi rho_Cl int g r^2 dr). Returns None when the
    molality's rdf files are not synced into the repo (they exist on
    cerberus for every conductivity run; only m1.0 was synced originally)."""
    # dir names are m0.25 / m0.5 / m1.0 / m2.0 / m4.0 (mixed precision)
    for pat in (f"m{mol:g}", f"m{mol:.1f}"):
        files = sorted(glob.glob(os.path.join(MADRID_RDF, pat,
                                              "L*_s*", "*.rdf")))
        if files:
            break
    if not files:
        return None
    r, gs, coords = None, [], []
    for f in files:
        dat = np.loadtxt(f, skiprows=4)
        if r is None:
            r = dat[:, 1]
        elif not np.allclose(dat[:, 1], r):
            dat = np.column_stack([dat[:, 0], r,
                                   np.interp(r, dat[:, 1], dat[:, 2]),
                                   np.interp(r, dat[:, 1], dat[:, 3])])
        gs.append(dat[:, 2])
        coords.append(dat[:, 3])
    g = np.mean(gs, axis=0)
    w = np.full_like(g, np.nan)
    m = g > 0.02
    w[m] = -KT_KCAL * np.log(g[m])
    w -= np.nanmean(w[(r > 9.0) & (r < r.max())])
    d = {"mol": mol, "n_seeds": len(files), "r": r, "g": g, "w": w}
    for key, lo, hi, kind in (("CIP", 2.4, 3.2, "min"), ("TS", 3.2, 4.2, "max"),
                              ("SSIP", 4.2, 5.6, "min")):
        mask = (r >= lo) & (r <= hi) & np.isfinite(w)
        i = (np.argmin if kind == "min" else np.argmax)(w[mask])
        d[f"r_{key}"], d[f"w_{key}"] = r[mask][i], w[mask][i]
    d["dW_cip_ssip"] = d["w_CIP"] - d["w_SSIP"]
    n_seed = np.array([np.interp(d["r_TS"], r, c) for c in coords])
    d["n_CIP"], d["n_CIP_sem"] = n_seed.mean(), sem(n_seed)
    return d


def report_pairing(pair, mad_pair):
    print("\n=== Na-Cl association (CSD3 round-2 P0 products; tab:pmf "
          "windows) ===")
    print("  level      CIP (kcal/mol)   SSIP            dW_CIP-SSIP  "
          "n_CIP          n_hyd(Na)  n_hyd(Cl)")
    for key, name in ((1.0, "RPA 1 m"), (2.0, "RPA 2 m"), (4.0, "RPA 4 m"),
                      ("mp2", "MP2 1 m")):
        d = pair[key]
        print(f"  {name:9s}  {d['w_CIP']:+.2f} @ {d['r_CIP']:.2f}    "
              f"{d['w_SSIP']:+.2f} @ {d['r_SSIP']:.2f}   {d['dW_cip_ssip']:+.2f}"
              f"        {d['n_CIP']:.3f}+-{d['n_CIP_sem']:.3f}"
              f"    {d['n_hyd_Na']:.2f}       {d['n_hyd_Cl']:.2f}")
    for m in sorted(mad_pair):
        d = mad_pair.get(m)
        if d is None:
            print(f"  Madrid {m:g} m  -- rdf files not synced (fetch "
                  f"results/lammps_madrid/conductivity/m{m:g}* from cerberus)")
            continue
        print(f"  Madrid {m:g} m  {d['w_CIP']:+.2f} @ {d['r_CIP']:.2f}    "
              f"{d['w_SSIP']:+.2f} @ {d['r_SSIP']:.2f}   {d['dW_cip_ssip']:+.2f}"
              f"        {d['n_CIP']:.3f}+-{d['n_CIP_sem']:.3f}"
              f"    (n={d['n_seeds']} runs, barrier @ {d['r_TS']:.2f} A)")
    print(f"  (cross-check, report Table tab:pmf Madrid 1 m: CIP +0.11, "
          f"SSIP -0.57, dW +0.68, n_CIP {MADRID_NCIP_1M:.3f})")
    print("  regression vs O'Neill 2024 (pinned from the paper): published "
          "RPA/MP2 have CIP and SSIP degenerate within ~0.2 kcal/mol and a "
          "CIP-SSIP barrier of ~1.6 kcal/mol (their Figs 3a/4b); CIP depth "
          "~ -1.1 to -1.2 (Fig 4a).")
    print(f"    ours at RPA 1 m: CIP-SSIP {pair[1.0]['dW_cip_ssip']:+.2f} "
          f"(degenerate, PASS), barrier {pair[1.0]['dW_kcalmol']:.2f} "
          f"(vs ~1.6), CIP depth {pair[1.0]['w_CIP']:+.2f} (shallower; their "
          "PMF cells run at the NpT density of each level -- RPA +5% dense "
          "-- in 24.8 A cells with ~20 ns sampling, vs our fixed "
          "experimental-density 37.26 A transport cells and ~0.4 ns).")


def fig_pairing_series(pair, mad_pair):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.8, 4.0))
    ax1.plot(mad_pair[1.0]["r"], mad_pair[1.0]["w"], color=SLATE3, lw=1.6,
             label="Madrid 1 mol/kg")
    if mad_pair.get(4.0) is not None:
        ax1.plot(mad_pair[4.0]["r"], mad_pair[4.0]["w"], color=SLATE3,
                 lw=1.3, ls="--", label="Madrid 4 mol/kg")
    ax1.plot(pair["mp2"]["r"], pair["mp2"]["w"], color=INDIGO, lw=1.4,
             ls="--", label="MP2 anchor 1 mol/kg")
    for m, color in ((1.0, CAM_BLUE), (2.0, C_CL), (4.0, C_O)):
        ax1.plot(pair[m]["r"], pair[m]["w"], color=color, lw=1.6,
                 label=f"RPA {m:g} mol/kg")
    ax1.axhline(0, color=MUTED, lw=0.8)
    ax1.set_xlim(2.2, 7.0)
    ax1.set_ylim(-1.3, 2.2)
    ax1.legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="upper right")
    style(ax1, r"$r_{\mathrm{NaCl}}$ (Å)", r"$w(r)$ (kcal/mol)", "(a)")
    cr = [MOL_TO_M[m] for m in (1.0, 2.0, 4.0)]
    y = [pair[m]["n_CIP"] for m in (1.0, 2.0, 4.0)]
    ys = [pair[m]["n_CIP_sem"] for m in (1.0, 2.0, 4.0)]
    ax2.errorbar(cr, y, yerr=ys, fmt="o--", ms=6.5, color=CAM_BLUE, mew=0,
                 lw=1.3, elinewidth=1.0, capsize=2, label="RPA")
    ax2.errorbar([MOL_TO_M[1.0] - 0.07], [pair["mp2"]["n_CIP"]],
                 yerr=[pair["mp2"]["n_CIP_sem"]], fmt="s", ms=6.5,
                 color=INDIGO, mew=0, elinewidth=1.0, capsize=2,
                 label="MP2 anchor")
    mols_mad = [m for m in sorted(mad_pair) if mad_pair.get(m) is not None]
    ax2.errorbar([MOL_TO_M[m] for m in mols_mad],
                 [mad_pair[m]["n_CIP"] for m in mols_mad],
                 yerr=[mad_pair[m]["n_CIP_sem"] for m in mols_mad],
                 fmt="o-", ms=6.5, color=SLATE3, mfc="white", mew=1.6,
                 lw=1.4, elinewidth=1.0, capsize=2, label="Madrid-2019")
    ax2.set_xlim(0, 3.9)
    ax2.set_ylim(0, 0.45)
    ax2.legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="upper left")
    style(ax2, "c (mol/L)", r"$n_{\mathrm{CIP}}$ (contact Cl per Na)", "(b)")
    savefig(fig, "fig13_pairing_series")


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
    with open(os.path.join(RPAT, "timings.csv")) as fh:
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


def madrid_fong(mad):
    """Per-seed Lambda / uNa / uCl for the Madrid series (Hittorf partition)."""
    sNN, sCC, sNC = (mad[f"runs_s{s}"] for s in ("NaNa", "ClCl", "NaCl"))
    s_sum = sNN + sCC - 2.0 * sNC
    tH = (mad["runs_tNa"]
          + mad["runs_mol"] * (M_NA * (sNN - sNC) - M_CL * (sCC - sNC)) / s_sum)
    c = mad["runs_nNa"] / (NA * mad["runs_V"] * 1e-27)
    fc = FARADAY * c * 1e3
    return {"mol": mad["runs_mol"], "c": c,
            "Lambda": 10.0 * mad["runs_kOns"] / c,
            "uNa": tH * mad["runs_kOns"] / fc * 1e8,
            "uCl": (1.0 - tH) * mad["runs_kOns"] / fc * 1e8}


def expt_fong(mols):
    """Lambda / uNa / uCl derived from the pinned kappa (Chambers-Stokes
    1956) and Hittorf t_Na (Smits-Duyvis 1966) anchors -- derived, not
    independently measured."""
    mols = np.asarray(mols, float)
    tH = tna_expt_hittorf(mols)
    k = np.array([EXPT_KAPPA[m] for m in mols])
    c = np.array([EXPT_C[m] for m in mols])
    fc = FARADAY * c * 1e3
    return c, 10.0 * k / c, tH * k / fc * 1e8, (1.0 - tH) * k / fc * 1e8


def report_fong(mad, rpa, mp2):
    mf = madrid_fong(mad)
    print("\n=== Fong-style observables (Lambda = kOns/c; u from the Hittorf "
          "partition t^H kOns / F c; expt derived from the pinned "
          "Chambers-Stokes kappa + Smits-Duyvis t^H) ===")
    for key, unit in (("Lambda", "S cm^2/mol"), ("uNa", "1e-8 m^2/Vs"),
                      ("uCl", "1e-8 m^2/Vs")):
        print(f"  {key} ({unit}):")
        for m in (1.0, 2.0, 4.0):
            g = mf["mol"] == m
            _, lam_e, una_e, ucl_e = expt_fong([m])
            exp_v = {"Lambda": lam_e, "uNa": una_e, "uCl": ucl_e}[key][0]
            line = (f"    m={m:g}: Madrid {mf[key][g].mean():6.2f}+-{sem(mf[key][g]):.2f}"
                    f"  RPA {rpa[m][key].mean():6.2f}+-{sem(rpa[m][key]):.2f}")
            if m == 1.0:
                line += f"  MP2 {mp2[key].mean():6.2f}+-{sem(mp2[key]):.2f}"
            print(line + f"  expt {exp_v:6.2f}")


def fig_fong_series(mad, rpa, mp2):
    mf = madrid_fong(mad)
    mols = np.array(sorted(set(mf["mol"])))
    ce, lam_e, una_e, ucl_e = expt_fong(mols)
    cr = np.array([MOL_TO_M[m] for m in sorted(rpa)])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.8, 4.0))
    # (a) molar conductivity vs sqrt(c), Kohlrausch presentation
    _, lam, lam_s = pooled(mf["mol"], mf["Lambda"])
    cm = np.array([MOL_TO_M[m] for m in mols])
    ax1.errorbar(np.sqrt(cm), lam, yerr=lam_s, fmt="o-", ms=6, color=CAM_BLUE,
                 mfc="white", mew=1.5, lw=1.4, elinewidth=1.0, capsize=2,
                 label="Madrid-2019")
    y = [rpa[m]["Lambda"].mean() for m in sorted(rpa)]
    ys = [sem(rpa[m]["Lambda"]) for m in sorted(rpa)]
    ax1.errorbar(np.sqrt(cr) + 0.03, y, yerr=ys, fmt="o--", ms=6.5,
                 color=CAM_BLUE, mew=0, lw=1.3, elinewidth=1.0, capsize=2,
                 label="RPA")
    ax1.errorbar([np.sqrt(MOL_TO_M[1.0]) - 0.03], [mp2["Lambda"].mean()],
                 yerr=[sem(mp2["Lambda"])], fmt="s", ms=6.5, color=INDIGO,
                 mew=0, elinewidth=1.0, capsize=2, label="MP2 anchor")
    ax1.plot(np.sqrt(ce), lam_e, marker="o", ms=5, ls="", color=SLATE3,
             label="Experiment (derived)", zorder=5)
    ax1.set_xlim(0.3, 2.05)
    ax1.legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="lower left")
    style(ax1, r"$\sqrt{c}$ ((mol/L)$^{1/2}$)",
          r"$\Lambda$ (S cm$^2$ mol$^{-1}$)", "(a)")
    # (b) electrophoretic mobilities (Hittorf frame; -u shown for Cl-)
    for key, color, mk, lab in (("uNa", CAM_BLUE, "o", r"Na$^+$"),
                                ("uCl", C_CL, "^", r"$-$Cl$^-$")):
        _, um, us = pooled(mf["mol"], mf[key])
        ax2.errorbar(cm, um, yerr=us, fmt=mk + "-", ms=5.5, color=color,
                     mfc="white", mew=1.5, lw=1.4, elinewidth=1.0, capsize=2,
                     label=f"Madrid {lab}")
        y = [rpa[m][key].mean() for m in sorted(rpa)]
        ys = [sem(rpa[m][key]) for m in sorted(rpa)]
        ax2.errorbar(cr + 0.07, y, yerr=ys, fmt=mk + "--", ms=6, color=color,
                     mew=0, lw=1.1, elinewidth=1.0, capsize=2,
                     label=f"RPA {lab}")
        ax2.errorbar([MOL_TO_M[1.0] - 0.07], [mp2[key].mean()],
                     yerr=[sem(mp2[key])], fmt="s", ms=6, color=INDIGO, mew=0,
                     elinewidth=1.0, capsize=2,
                     label="MP2 anchor" if key == "uNa" else None)
    ax2.plot(ce, una_e, marker="o", ms=5, ls="", color=SLATE3,
             label="Experiment (derived)", zorder=5)
    ax2.plot(ce, ucl_e, marker="^", ms=5, ls="", color=SLATE3, zorder=5)
    ax2.set_xlim(0, 3.9)
    ax2.set_ylim(0, None)
    h, l = ax2.get_legend_handles_labels()
    ax2.legend(h, l, frameon=False, fontsize=8, labelcolor=INK2,
               loc="upper right", ncol=2, columnspacing=1.2)
    style(ax2, "c (mol/L)", r"$u$ ($10^{-8}$ m$^2$ V$^{-1}$ s$^{-1}$)", "(b)")
    savefig(fig, "fig14_fong_series")


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
    ce = [0.0] + [EXPT_C[m] for m in (1.0, 2.0, 4.0)]
    ax1.plot(ce, [ETA_EXP[m] for m in (0.0, 1.0, 2.0, 4.0)], marker="o",
             ms=5, ls="", color=SLATE3, label="Experiment", zorder=5)
    ax1.set_xlim(-0.12, 3.9)
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
    path = os.path.join(CMP, "transport_series_summary.csv")
    cols = ("level,m_molkg,n_seg,fit_ps,Dw_1e9,Dw_sem,DNa_1e9,DNa_sem,"
            "DCl_1e9,DCl_sem,rNa,rNa_sem,rCl,rCl_sem,kNE_Sm,kNE_sem,"
            "kOns_Sm,kOns_sem,dNE,dNE_sem,chNaNa,chNaNa_sem,chClCl,"
            "chClCl_sem,chNaCl,chNaCl_sem,tNa_bary,tNa_bary_sem,"
            "tNa_hittorf,tNa_hittorf_sem,Lambda_Scm2mol,Lambda_sem,"
            "uNa_1e8m2Vs,uNa_sem,uCl_1e8m2Vs,uCl_sem,eta_mPas,eta_sem\n")
    def row(level, mol, d, eta):
        f = [f"{level}", f"{mol:g}", f"{len(d['kOns'])}",
             f"{d['tmin']:g}-{d['tmax']:g}"]
        for k in ("DO", "DNa", "DCl"):
            f += [f"{d[k].mean()*CONV:.4f}", f"{sem(d[k])*CONV:.4f}"]
        for k in ("rNa", "rCl", "kNE", "kOns", "dNE", "chNaNa", "chClCl",
                  "chNaCl", "tNa", "tH", "Lambda", "uNa", "uCl"):
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
            mf = madrid_fong(mad)
            for k in ("Lambda", "uNa", "uCl"):
                d[k] = mf[k][g]
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
          f"  x1e-9 m^2/s;  expt at 1 m {DW_EXP[1.0]} "
          "(Muller-Hertz 1996 via Blazquez 2023 Table V)")
    print(f"  timescale factor vs expt D_w(1 m): "
          f"MP2 {DW_EXP[1.0]/(mp2['DO'].mean()*CONV):.2f}x  "
          f"RPA {DW_EXP[1.0]/(rpa[1.0]['DO'].mean()*CONV):.2f}x")
    r_exp = DW_EXP[4.0] / DW_EXP[1.0]
    r_rpa = rpa[4.0]["DO"].mean() / rpa[1.0]["DO"].mean()
    print(f"  water-slowdown response D_w(4m)/D_w(1m): expt {r_exp:.3f}  "
          f"RPA {r_rpa:.3f}  (Madrid-2019 published, Blazquez Table V: "
          f"{1.29/2.02:.3f})")

    report_fong(mad, rpa, mp2)

    pair = {k: load_pairing(tag) for k, tag in PAIRING.items()}
    mad_pair = {m: madrid_pairing(m) for m in (0.25, 0.5, 1.0, 2.0, 4.0)}
    report_pairing(pair, mad_pair)

    window_sensitivity()
    ledger()
    print("\nAll experimental anchors pinned 2026-08-02: kappa Chambers-"
          "Stokes 1956, t_Na Smits-Duyvis 1966, eta Kestin 1981 (Zotero "
          "8R3TJWV9), D_w Muller-Hertz 1996 via Blazquez 2023 Table V "
          "(Zotero VMZZ3IN9), O'Neill 2024 PMF benchmarks from the paper.")

    fig_kappa_series(mad, rpa, mp2)
    fig_dne_series(mad, rpa, mp2)
    fig_tna_series(mad, rpa, mp2)
    fig_eta_series(rpa_eta, mp2_eta)
    fig_pairing_series(pair, mad_pair)
    fig_fong_series(mad, rpa, mp2)
    write_summary(mad, rpa, rpa_eta, mp2, mp2_eta)


if __name__ == "__main__":
    main()
