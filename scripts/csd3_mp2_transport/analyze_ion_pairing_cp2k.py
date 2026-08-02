#!/usr/bin/env python3
"""Na-Cl association analysis (g(r) / PMF / n_CIP) for CP2K NaCl(aq) segments.

Port of lammps/pair_association.py + the RDF half of analyze_ions_rdf.py to
the transport-campaign layout. Reads the reduced_traj_e{K}.npz caches written
by analyze_onsager_cp2k.py (O/Na/Cl every K dumped frames), so it never
re-streams the pos files; run the Onsager analysis first if a cache is
missing.

Per concentration (segments pooled for curves, per-segment for SEMs):

  g_NaCl(r), w(r) = -kT ln g referenced to the 9 A..L/2 tail   (PMF)
  r_peak / w_min (contact basin), r_b / w_b (CIP-SSIP barrier), dW = w_b - w_min
  n_CIP = 4 pi rho_Cl int_0^rb g r^2 dr    per-segment -> mean +- SEM
  K_A   = 4 pi N_A int_0^rb e^{-w/kT} r^2 dr  (L/mol, Bjerrum-style) and the
          ideal paired fraction x/c from K_A = x/(c-x)^2
  g_NaO(r), g_ClO(r) + first-shell hydration numbers n_hyd (same integral to
          the first minimum of the pooled g)

Outputs (results/mp2_transport/):
  {label}_rdf_ions.csv   r_A, g_NaCl, w_NaCl_kcalmol, g_NaO, g_ClO
  {label}_pairing.csv    one summary row (assemble rows into the PMF table)
  {label}_pairing.npz    per-segment histograms + all scalars

Usage:
  analyze_ion_pairing_cp2k.py --segdirs runs_gpu/RPA/cubic_4m/production/cube3/seg* \
      --box-a 37.26 --label rpa_gpu_cube3_4m [--temp 298.15] [--every 20]
  analyze_ion_pairing_cp2k.py --selftest
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_onsager_cp2k import reduce_segment, EL_CODE  # noqa: E402

RESULTS = os.path.expanduser("~/cp2k-benchmarks/results/mp2_transport")
NA_AVOG = 6.02214076e23
RDF_DR = 0.05  # A


def pair_hist(A, B, box_l, edges):
    """Min-image A-B distance histogram over frames. A: (T,nA,3), B: (T,nB,3)."""
    h = np.zeros(len(edges) - 1)
    if A.shape[1] == 0 or B.shape[1] == 0:
        return h
    chunk = max(1, int(2e7 / (A.shape[1] * B.shape[1])))
    for i in range(0, A.shape[0], chunk):
        d = A[i:i + chunk, :, None, :] - B[i:i + chunk, None, :, :]
        d -= box_l * np.round(d / box_l)
        h += np.histogram(np.sqrt(np.square(d).sum(-1)).ravel(), bins=edges)[0]
    return h


def coord_number(r, g, rho, r_cut):
    """n(r_cut) = 4 pi rho int_0^rcut g r^2 dr (trapezoid on the binned g)."""
    m = r <= r_cut
    return 4.0 * np.pi * rho * np.trapezoid(g[m] * r[m] ** 2, r[m])


def first_extrema(r, g, peak_win, min_span=1.6):
    """(r_peak, r_min): first maximum of g inside peak_win, then the first
    minimum within min_span A after it."""
    pm = (r >= peak_win[0]) & (r <= peak_win[1])
    r_peak = r[pm][np.argmax(g[pm])]
    mm = (r >= r_peak) & (r <= r_peak + min_span)
    r_min = r[mm][np.argmin(g[mm])]
    return r_peak, r_min


def _selftest():
    """Ideal gas: g = 1 everywhere, n(rc) = 4/3 pi rho rc^3."""
    rng = np.random.default_rng(11)
    L, n = 20.0, 400
    A = rng.uniform(0, L, size=(60, n, 3))
    B = rng.uniform(0, L, size=(60, n, 3))
    edges = np.arange(0.0, L / 2 + RDF_DR, RDF_DR)
    h = pair_hist(A, B, L, edges)
    r = 0.5 * (edges[:-1] + edges[1:])
    shell = 4.0 / 3.0 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    g = h / (60 * n * n / L ** 3 * shell)
    m = r > 2.0
    assert abs(np.mean(g[m]) - 1.0) < 0.01, "ideal-gas g(r) != 1"
    n_c = coord_number(r, g, n / L ** 3, 6.0)
    assert abs(n_c / (4 / 3 * np.pi * n / L ** 3 * 6.0 ** 3) - 1) < 0.02, \
        "coordination integral off on ideal gas"
    print("selftest OK")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--segdirs", nargs="+")
    ap.add_argument("--box-a", type=float, help="cubic box edge in A")
    ap.add_argument("--temp", type=float, default=298.15,
                    help="kT for the PMF (298.15 matches the MP2 anchor rows)")
    ap.add_argument("--every", type=int, default=20,
                    help="cache stride: reads reduced_traj_e{K}.npz")
    ap.add_argument("--tmax-ps", type=float, default=80.0,
                    help="truncate segments to a common length (overruns)")
    ap.add_argument("--rb-win", type=float, nargs=2, default=(3.2, 4.2),
                    help="search window for the CIP/SSIP barrier [A]")
    ap.add_argument("--label", default="run")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        _selftest()
        return
    if not (args.segdirs and args.box_a):
        ap.error("--segdirs and --box-a required")

    L = args.box_a
    V = L ** 3
    kt = 0.0019872041 * args.temp
    edges = np.arange(0.0, L / 2 + RDF_DR, RDF_DR)
    r = 0.5 * (edges[:-1] + edges[1:])
    shell = 4.0 / 3.0 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)

    hists = {k: [] for k in ("NaCl", "NaO", "ClO")}
    frames = []
    for sd in args.segdirs:
        steps, codes, traj = reduce_segment(sd.rstrip("/"), args.every)
        keep = (steps - steps[0]) * 0.0005 <= args.tmax_ps  # 0.5 fs steps
        traj = traj[keep]
        sel = {n: traj[:, codes == EL_CODE[n], :] for n in ("Na", "Cl", "O")}
        hists["NaCl"].append(pair_hist(sel["Na"], sel["Cl"], L, edges))
        hists["NaO"].append(pair_hist(sel["Na"], sel["O"], L, edges))
        hists["ClO"].append(pair_hist(sel["Cl"], sel["O"], L, edges))
        frames.append(traj.shape[0])
        print(f"{sd}: {traj.shape[0]} frames", flush=True)
    n_na, n_cl, n_o = (int((codes == EL_CODE[n]).sum()) for n in ("Na", "Cl", "O"))
    frames = np.array(frames)

    def g_of(key, n_a, n_b, fr, hs):
        return np.sum(hs, axis=0) / (fr.sum() * n_a * n_b / V * shell)

    g_pool = {k: g_of(k, *ab, frames, hists[k])
              for k, ab in (("NaCl", (n_na, n_cl)), ("NaO", (n_na, n_o)),
                            ("ClO", (n_cl, n_o)))}

    # PMF from the pooled Na-Cl g(r), tail-referenced (as the MP2 anchor)
    g = g_pool["NaCl"]
    w = np.full_like(g, np.nan)
    m = g > 0.02
    w[m] = -kt * np.log(g[m])
    tail = (r > 9.0) & (r < r.max()) & np.isfinite(w)
    w -= np.nanmean(w[tail])

    pk = (r >= 2.0) & (r <= 3.5)
    r_peak = r[pk][np.argmax(g[pk])]
    w_min = np.nanmin(w[pk])
    bm = (r >= args.rb_win[0]) & (r <= args.rb_win[1]) & np.isfinite(w)
    r_b = r[bm][np.argmax(w[bm])]
    w_b = np.nanmax(w[bm])

    # n_CIP per segment (pooled r_b) -> mean +- SEM
    rho_cl = n_cl / V
    n_cip_seg = np.array([
        coord_number(r, h / (f * n_na * n_cl / V * shell), rho_cl, r_b)
        for h, f in zip(hists["NaCl"], frames)])
    n_cip, n_cip_sem = n_cip_seg.mean(), n_cip_seg.std(ddof=1) / np.sqrt(len(n_cip_seg))

    # Bjerrum-style K_A from the pooled PMF + ideal paired fraction
    core = r <= r_b
    wb = np.where(np.isfinite(w), w, np.inf)[core]
    ka_lmol = (4 * np.pi * np.trapezoid(np.exp(-wb / kt) * r[core] ** 2, r[core])
               * NA_AVOG * 1e-27)
    c_molar = rho_cl / NA_AVOG * 1e27
    coeff = [ka_lmol, -(2 * ka_lmol * c_molar + 1), ka_lmol * c_molar ** 2]
    paired_frac = min(np.roots(coeff)) / c_molar

    # hydration shells from the pooled ion-O g(r), same integral per segment
    hyd = {}
    for key, n_ion, rho_o_pairs in (("NaO", n_na, n_o), ("ClO", n_cl, n_o)):
        r_pk, r_min = first_extrema(r, g_pool[key], (2.0, 3.6))
        n_seg = np.array([
            coord_number(r, h / (f * n_ion * n_o / V * shell), n_o / V, r_min)
            for h, f in zip(hists[key], frames)])
        hyd[key] = (r_pk, r_min, n_seg.mean(),
                    n_seg.std(ddof=1) / np.sqrt(len(n_seg)))

    os.makedirs(RESULTS, exist_ok=True)
    np.savetxt(os.path.join(RESULTS, f"{args.label}_rdf_ions.csv"),
               np.column_stack([r, g, w, g_pool["NaO"], g_pool["ClO"]]),
               delimiter=",", comments="",
               header="r_A,g_NaCl,w_NaCl_kcalmol,g_NaO,g_ClO")
    hdr = ("label,n_NaCl_pairs,c_mol_L,frames,r_peak_A,w_min_kcalmol,r_b_A,"
           "w_b_kcalmol,dW_kcalmol,n_CIP,n_CIP_sem,K_A_L_mol,paired_frac,"
           "r_NaO_peak_A,r_NaO_min_A,n_hyd_Na,n_hyd_Na_sem,"
           "r_ClO_peak_A,r_ClO_min_A,n_hyd_Cl,n_hyd_Cl_sem")
    row = (f"{args.label},{n_na},{c_molar:.4f},{frames.sum()},{r_peak:.3f},"
           f"{w_min:.4f},{r_b:.3f},{w_b:.4f},{w_b - w_min:.4f},{n_cip:.4f},"
           f"{n_cip_sem:.4f},{ka_lmol:.4f},{paired_frac:.4f},"
           f"{hyd['NaO'][0]:.3f},{hyd['NaO'][1]:.3f},{hyd['NaO'][2]:.4f},"
           f"{hyd['NaO'][3]:.4f},{hyd['ClO'][0]:.3f},{hyd['ClO'][1]:.3f},"
           f"{hyd['ClO'][2]:.4f},{hyd['ClO'][3]:.4f}")
    with open(os.path.join(RESULTS, f"{args.label}_pairing.csv"), "w") as f:
        f.write(hdr + "\n" + row + "\n")
    np.savez(os.path.join(RESULTS, f"{args.label}_pairing.npz"),
             r=r, g_NaCl=g, w_NaCl=w, g_NaO=g_pool["NaO"], g_ClO=g_pool["ClO"],
             hist_NaCl_seg=np.array(hists["NaCl"]), frames=frames,
             n_na=n_na, n_cl=n_cl, n_o=n_o, box_a=L, temp=args.temp,
             r_peak=r_peak, w_min=w_min, r_b=r_b, w_b=w_b,
             n_cip_seg=n_cip_seg, ka_lmol=ka_lmol, c_molar=c_molar,
             paired_frac=paired_frac)

    print(f"\n=== {args.label}: {n_na} Na / {n_cl} Cl, c = {c_molar:.3f} mol/L, "
          f"{frames.sum()} frames pooled ===")
    print(f"contact peak r = {r_peak:.2f} A, w_min = {w_min:.3f} kcal/mol")
    print(f"barrier r_b = {r_b:.2f} A, w_b = {w_b:.3f}, dW = {w_b - w_min:.3f}")
    print(f"n_CIP = {n_cip:.3f} +- {n_cip_sem:.3f} Cl within r_b per Na")
    print(f"K_A = {ka_lmol:.3f} L/mol, ideal paired fraction = {paired_frac:.3f}")
    print(f"n_hyd Na-O = {hyd['NaO'][2]:.2f} +- {hyd['NaO'][3]:.2f} "
          f"(min {hyd['NaO'][1]:.2f} A), Cl-O = {hyd['ClO'][2]:.2f} +- "
          f"{hyd['ClO'][3]:.2f} (min {hyd['ClO'][1]:.2f} A)")
    print(f"saved {args.label}_rdf_ions.csv / _pairing.csv / _pairing.npz")


if __name__ == "__main__":
    main()
