#!/usr/bin/env python3
"""Independent audit of every scientific claim in the NaCl(aq) paper notes.

Recomputes each headline number from the data products bundled in the repo
(results/lammps_madrid/*, results/nacl_mp2_anchor/*) with fresh code —
NOT by rerunning make_figures.py — and prints claim vs. recomputed value.
"""
import glob
import os
import re

import numpy as np

REPO = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     "..", "..", ".."))
ROOT = os.path.join(REPO, "results", "lammps_madrid")
ANCHOR = os.path.join(REPO, "results", "nacl_mp2_anchor")

KB = 1.380649e-23
T = 298.15
KT_KCAL = 0.0019872041 * T
E2 = (1.602176634e-19) ** 2
XI = 2.837297
NA_AVOG = 6.02214076e23

MOL_TO_M = {0.25: 0.231, 0.5: 0.491, 1.0: 0.960, 2.0: 1.915, 4.0: 3.666}
EXPT_C = {0.25: 0.248, 0.5: 0.494, 1.0: 0.979, 2.0: 1.920, 4.0: 3.686}
EXPT_KAPPA = {0.25: 2.48, 0.5: 4.63, 1.0: 8.42, 2.0: 14.50, 4.0: 22.04}
SD66_M = np.array([0.024, 0.048, 0.096, 0.192, 0.384, 0.768, 1.536, 3.072, 6.144])
SD66_LOGA = np.array([-1.683, -1.403, -1.126, -0.848, -0.573, -0.292,
                      0.0036, 0.3436, 0.7897])
M_NA, M_CL = 0.0229898, 0.0354530

kap = np.load(os.path.join(ROOT, "conductivity", "analysis", "kappa_vs_c.npz"))
rep = np.load(os.path.join(ROOT, "replicas", "analysis", "replica_D.npz"),
              allow_pickle=True)

mol = kap["runs_mol"]
mols = np.array(sorted(set(mol)))


def pooled(y):
    mean = np.array([y[mol == m].mean() for m in mols])
    sem = np.array([y[mol == m].std(ddof=1) / np.sqrt((mol == m).sum())
                    for m in mols])
    return mean, sem


print("=" * 78)
print("CLAIM 1 (fig6): Na-Cl Onsager cross term = 0 within error at every m;")
print("NE deviation dominated by like-ion (mostly Na-Na) anticorrelation")
print("=" * 78)
pref = E2 / (KB * T * kap["runs_V"] * 1e-30) * 1e-8 / 6.0   # S/m per A^2/ps slope
kNaCl = pref * kap["runs_sNaCl"]
kNE = kap["runs_kNE"]
dNaNa = pref * (kap["runs_sNaNa"] - 6.0 * kap["runs_nNa"] * kap["runs_DNa"])
dClCl = pref * (kap["runs_sClCl"] - 6.0 * kap["runs_nCl"] * kap["runs_DCl"])
cross_m, cross_s = pooled(-2 * kNaCl)             # cross contribution to kappa
tot_m, tot_s = pooled(1.0 - kap["runs_kOns"] / kNE)
pNa_m, pNa_s = pooled(-dNaNa / kNE)
pCl_m, pCl_s = pooled(-dClCl / kNE)
pX_m, pX_s = pooled(2 * kNaCl / kNE)
print(f"{'m':>5} {'-2k_NaCl (S/m)':>18} {'|z|':>5} | NE-dev total    "
      f"Na-Na part      Cl-Cl part      Na-Cl(pairing) part")
for i, m in enumerate(mols):
    z = abs(cross_m[i]) / cross_s[i]
    print(f"{m:5.2f} {cross_m[i]:+8.3f} ± {cross_s[i]:5.3f} {z:5.1f} | "
          f"{tot_m[i]:+.3f}±{tot_s[i]:.3f}  {pNa_m[i]:+.3f}±{pNa_s[i]:.3f}  "
          f"{pCl_m[i]:+.3f}±{pCl_s[i]:.3f}  {pX_m[i]:+.3f}±{pX_s[i]:.3f}")
print("-> cross term zero within ~2 SEM at all m?  "
      + str(all(abs(cross_m) < 2 * cross_s)))
print("-> Na-Na part > Cl-Cl part at all m?  "
      + str(all(pNa_m > pCl_m)))

print()
print("=" * 78)
print("CLAIM 2 (fig3): kappa vs experiment (Chambers-Stokes 1956 anchors)")
print("=" * 78)
kOns_m, kOns_s = pooled(kap["runs_kOns"])
kNE_m, kNE_s = pooled(kNE)
print(f"{'m':>5} {'expt':>6} {'Ons z=.85':>10} {'dev%':>6} "
      f"{'NE z=1':>8} {'dev%':>6} {'Ons z=1':>8} {'dev%':>6}")
for i, m in enumerate(mols):
    ke = EXPT_KAPPA[m]
    k85 = kOns_m[i] * 0.85 ** 2
    print(f"{m:5.2f} {ke:6.2f} {k85:7.2f}±{kOns_s[i]*0.7225:4.2f} "
          f"{100*(k85/ke-1):+6.1f} {kNE_m[i]:8.2f} {100*(kNE_m[i]/ke-1):+6.1f} "
          f"{kOns_m[i]:8.2f} {100*(kOns_m[i]/ke-1):+6.1f}")

print()
print("=" * 78)
print("CLAIM 3 (fig4b): t_Na barycentric + Hittorf-transformed vs Smits-Duyvis")
print("=" * 78)
sNN, sCC, sNC = kap["runs_sNaNa"], kap["runs_sClCl"], kap["runs_sNaCl"]
s_sum = sNN + sCC - 2.0 * sNC
tH_runs = kap["runs_tNa"] + mol * (M_NA * (sNN - sNC) - M_CL * (sCC - sNC)) / s_sum
tb_m, tb_s = pooled(kap["runs_tNa"])
tH_m, tH_s = pooled(tH_runs)
loga = np.interp(np.log(mols), np.log(SD66_M), SD66_LOGA)
t_expt = 0.3720 - 0.0118 * loga
print(f"{'m':>5} {'bary':>13} {'Hittorf':>13} {'expt(SD66)':>10} {'|z| (H vs expt)':>15}")
for i, m in enumerate(mols):
    z = abs(tH_m[i] - t_expt[i]) / tH_s[i]
    print(f"{m:5.2f} {tb_m[i]:.3f}±{tb_s[i]:.3f} {tH_m[i]:.3f}±{tH_s[i]:.3f} "
          f"{t_expt[i]:10.3f} {z:15.1f}")

print()
print("=" * 78)
print("CLAIM 4 (fig2): Yeh-Hummer D0 values + implied viscosity (Madrid, 1 m)")
print("=" * 78)
Ls = np.array(rep["Ls"])
species = list(rep["species"])
# independent refit: free fit for O; Na/Cl with slope fixed to O's
invL = 1.0 / Ls
D = rep["D_mean"]
w = 1.0 / rep["D_sem"] ** 2
iO = species.index("O")
mO, bO = np.polyfit(invL, D[iO], 1, w=np.sqrt(w[iO]))
print(f"O : free fit slope {mO:+.4f} A^2/ps*A, D0 = {bO:.4f} A^2/ps = "
      f"{bO*1e-8*1e9:.3f} e-9 m2/s   (npz free_fits: {rep['free_fits'][0][:2]})")
for sp, row in (("Na", 0), ("Cl", 1)):
    i = species.index(sp)
    b = np.average(D[i] - mO * invL, weights=w[i])
    print(f"{sp:2}: shared-slope D0 = {b:.4f} A^2/ps = {b*1e-8*1e9:.3f} e-9 m2/s"
          f"   (npz shared_fits: {rep['shared_fits'][row][1]:.4f})")
slope_SI = mO * 1e-18            # A^2/ps per 1/A  ->  m^3/s
eta = -KB * T * XI / (6 * np.pi * slope_SI)
print(f"implied viscosity from O slope: eta = {eta*1e3:.3f} mPa s "
      f"(TIP4P/2005 pure water ~0.83; expt 1 m NaCl ~0.97)")
print(f"Madrid water D0 = {bO*10:.3f} e-9 vs expt pure-water 2.299 e-9 (Holz 2000)"
      f" -> {100*(bO*10/2.299-1):+.1f}%")

print()
print("=" * 78)
print("CLAIM 5 (fig1): Einstein (MSD) vs Green-Kubo (VACF) D agree within 2 sigma")
print("=" * 78)
vacf = {}
for f in sorted(glob.glob(os.path.join(ROOT, "vacf_replica", "L*_s*", "*.vacf"))):
    L = float(re.search(r"L([\d.]+)_s", f).group(1))
    dat = np.loadtxt(f, comments="#")
    with open(f) as fh:
        names = fh.readline().lstrip("# ").split()
    t = dat[:, names.index("t_ps")]
    msel = (t >= 3.0) & (t <= 5.0)
    for sp in ("O", "Na", "Cl"):
        vacf.setdefault((L, sp), []).append(dat[msel, names.index(f"D_{sp}")].mean())
worst = 0.0
for j, L in enumerate(Ls):
    for i, sp in enumerate(species):
        v = np.array(vacf[(L, sp)])
        vm, vs = v.mean(), v.std(ddof=1) / np.sqrt(len(v))
        z = abs(D[i, j] - vm) / np.hypot(rep["D_sem"][i, j], vs)
        worst = max(worst, z)
        print(f"L={L:5.2f} {sp:2}: MSD {D[i,j]:.4f}±{rep['D_sem'][i,j]:.4f}  "
              f"VACF {vm:.4f}±{vs:.4f}  |z|={z:.2f}")
print(f"-> worst |z| over all 12 = {worst:.2f}  (claim: all within 2)")

print()
print("=" * 78)
print("CLAIM 6 (fig5/fig7a): PMF  Madrid dW(CIP-SSIP) ~ +0.7, MP2 ~ +0.1 kcal/mol")
print("=" * 78)


def extrema(r, wv):
    def ext(lo, hi, kind):
        m = (r >= lo) & (r <= hi) & np.isfinite(wv)
        i = (np.argmin if kind == "min" else np.argmax)(wv[m])
        return r[m][i], wv[m][i]
    return ext(2.4, 3.2, "min"), ext(3.2, 4.2, "max"), ext(4.2, 5.6, "min")


gs, r = [], None
for f in sorted(glob.glob(os.path.join(ROOT, "conductivity", "m1.0", "L*_s*", "*.rdf"))):
    dat = np.loadtxt(f, skiprows=4)
    if r is None:
        r = dat[:, 1]
    gs.append(np.interp(r, dat[:, 1], dat[:, 2]) if not np.allclose(dat[:, 1], r)
              else dat[:, 2])
g = np.mean(gs, axis=0)
wmad = np.full_like(g, np.nan)
mask = g > 0.02
wmad[mask] = -KT_KCAL * np.log(g[mask])
wmad -= np.nanmean(wmad[(r > 9.0) & (r < r.max())])
(rc, wc), (rt, wt), (rs, ws) = extrema(r, wmad)
print(f"Madrid (10 seed-avg RDFs, 1 m): CIP {wc:+.3f} @ {rc:.2f} A, "
      f"barrier {wt:+.3f} @ {rt:.2f}, SSIP {ws:+.3f} @ {rs:.2f}, "
      f"dW = {wc-ws:+.3f} kcal/mol")
for cell in ("cube2", "cube3"):
    rr, gg, ww = np.loadtxt(os.path.join(ANCHOR, f"rdf_nacl_{cell}.csv"),
                            delimiter=",", skiprows=1).T
    (rc2, wc2), _, (rs2, ws2) = extrema(rr, ww)
    print(f"MP2 {cell}: CIP {wc2:+.3f} @ {rc2:.2f}, SSIP {ws2:+.3f} @ {rs2:.2f}, "
          f"dW = {wc2-ws2:+.3f} kcal/mol")

print()
print("=" * 78)
print("CLAIM 7 (pair_association): ~2.5% contact pairs at 1 m")
print("=" * 78)
bm = (r >= 3.2) & (r <= 4.2) & np.isfinite(wmad)
r_b = r[bm][np.argmax(wmad[bm])]
sel = mol == 1.0
rho = (kap["runs_nCl"][sel] / kap["runs_V"][sel]).mean()      # A^-3
core = r <= r_b
n_cip = 4 * np.pi * rho * np.trapezoid(g[core] * r[core] ** 2, r[core])
wint = np.where(np.isfinite(wmad), wmad, np.inf)
ka = 4 * np.pi * NA_AVOG * 1e-27 * np.trapezoid(
    np.exp(-wint[core] / KT_KCAL) * r[core] ** 2, r[core])
c_molar = rho / NA_AVOG * 1e27
print(f"barrier r_b = {r_b:.2f} A, rho_Cl = {c_molar:.3f} mol/L")
print(f"n_CIP (direct coordination) = {n_cip:.4f}  -> {100*n_cip:.2f}% of Na+ paired")
print(f"K_A (Bjerrum-style) = {ka:.3f} L/mol -> paired fraction ~ "
      f"{100*ka*c_molar:.2f}% (self-consistency ignored)")
print(f"compare: total NE deviation at 1 m = {tot_m[list(mols).index(1.0)]*100:.1f}%")

print()
print("=" * 78)
print("CLAIM 8 (fig7b): D_Cl/D_w  MP2 0.76 | Madrid 0.68 | expt 0.88; MP2 ~3x slow")
print("=" * 78)
wat = np.genfromtxt(os.path.join(ANCHOR, "diffusion_summary.csv"),
                    delimiter=",", names=True, dtype=None, encoding=None)
ion = np.genfromtxt(os.path.join(ANCHOR, "ion_diffusion_summary.csv"),
                    delimiter=",", names=True, dtype=None, encoding=None)
for cell in ("cube2", "cube3"):
    dO = wat["D_1e9_m2_s"][wat["cell"] == cell]
    out = [f"MP2 {cell}: O {dO.mean():.3f}"]
    for sp in ("Na", "Cl"):
        di = ion["D_1e9_m2_s"][(ion["cell"] == cell) & (ion["species"] == sp)]
        ratio = di.mean() / dO.mean()
        sem = ratio * np.hypot(di.std(ddof=1) / np.sqrt(len(di)) / di.mean(),
                               dO.std(ddof=1) / np.sqrt(len(dO)) / dO.mean())
        out.append(f"{sp} {di.mean():.3f}  {sp}/O = {ratio:.3f}±{sem:.3f}")
    print("  ".join(out))
madNa = rep["shared_fits"][0][1] / rep["free_fits"][0][1]
madCl = rep["shared_fits"][1][1] / rep["free_fits"][0][1]
print(f"Madrid D0 ratios: Na/O = {madNa:.3f}  Cl/O = {madCl:.3f}")
print(f"expt (inf. dilution): Na/O = {1.334/2.299:.3f}  Cl/O = {2.032/2.299:.3f}")
# absolute slowdown of MP2 water
d2 = wat["D_1e9_m2_s"][wat["cell"] == "cube2"].mean()
d3 = wat["D_1e9_m2_s"][wat["cell"] == "cube3"].mean()
x2, x3 = 1 / 24.84, 1 / 37.26
slope = (d2 - d3) / (x2 - x3)
d0 = d2 - slope * x2
print(f"MP2 water: cube2 {d2:.3f}, cube3 {d3:.3f}, YH-extrapolated D0 {d0:.3f} e-9")
print(f"slowdown vs expt 2.299: cube2 {2.299/d2:.2f}x, cube3 {2.299/d3:.2f}x, "
      f"D0 {2.299/d0:.2f}x")
