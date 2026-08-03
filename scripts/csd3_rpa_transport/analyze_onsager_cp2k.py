#!/usr/bin/env python3
"""Self + collective (Onsager) transport analysis for CP2K NaCl(aq) segments.

Port of cerberus onsager_multiorigin.py + kappa_analysis.py to the CP2K
campaign layout: the FFT multi-origin machinery is verbatim (validated by
--selftest against brute force); only the trajectory reader changed.

Per segment (an awk fast-path reduces the ~10 GB pos file to ions + O at a
coarser time grid, with a monotonic step-number slice for walltime-resumed
segments; the reduced array is cached as reduced_traj.npz in the segdir):

  self MSDs  -> D_Na, D_Cl, D_O   (Einstein fit over the lag window)
  collective C_ss'(tau) = <[R_s(t+tau)-R_s(t)].[R_s'(t+tau)-R_s'(t)]>
             for R_s = sum of species positions, ss' in NaNa/ClCl/NaCl

Aggregate over segments (mean +- SEM, each segment one independent sample):

  kappa_NE   = e^2 (N_Na D_Na + N_Cl D_Cl) / (kB T V)        [z = +-1]
  kappa_Ons  = e^2 slope(C_NaNa + C_ClCl - 2 C_NaCl) / (6 kB T V)
  Delta_NE   = 1 - kappa_Ons / kappa_NE
  t_Na       = (slope C_NaNa - slope C_NaCl) / slope(C_sum)
  channel decomposition kNaNa/kClCl/kNaCl incl. self/distinct split

Usage:
  analyze_onsager_cp2k.py --segdirs runs/MP2/cubic_1M/production/cube3/seg* \
      --box-a 37.26 --label mp2_anchor_cube3 [--every 20] [--tmin 10 --tmax 40]
  analyze_onsager_cp2k.py --selftest
"""
import argparse
import glob
import os
import subprocess
import sys

import numpy as np

E2 = (1.602176634e-19) ** 2      # C^2
KB = 1.380649e-23                # J/K
A2PS_TO_M2S = 1.0e-8             # A^2/ps -> m^2/s
A3_TO_M3 = 1.0e-30
RESULTS_ROOT = os.path.expanduser("~/cp2k-benchmarks/results")


def results_dir(label):
    """Level-separated outputs: mp2* labels go to results/mp2_transport/,
    everything else (rpa_*) to results/rpa_transport/."""
    sub = "mp2_transport" if label.startswith("mp2") else "rpa_transport"
    return os.path.join(RESULTS_ROOT, sub)

EL_CODE = {"Na": 1, "Cl": 2, "O": 3}

# awk reducer: keep every Kth frame (monotonic in step number), emit
# "step elcode x y z" for Na/Cl/O only. Comment line: " i = N, time = ..."
AWK = r'''
{
  r = (NR - 1) % (NAT + 2);
  if (r == 1) {
    for (i = 1; i <= NF; i++) if ($i == "i") { s = $(i+2); break }
    gsub(",", "", s); step = s + 0;
    if (nf == 0 || step > seen) { seen = step; nf++; keep = (nf % K == 1 || K == 1) }
    else keep = 0;
  } else if (r >= 2 && keep) {
    c = ($1 == "Na") ? 1 : ($1 == "Cl") ? 2 : ($1 == "O") ? 3 : 0;
    if (c) print step, c, $2, $3, $4;
  }
}
'''


def _prefix(p):
    P = np.zeros(len(p) + 1)
    P[1:] = np.cumsum(p)
    return P


def _crosscorr_fft(a, b):
    n = len(a)
    Fa = np.fft.fft(a, 2 * n)
    Fb = np.fft.fft(b, 2 * n)
    res = np.fft.ifft(Fa.conjugate() * Fb)[:n].real
    return res / (n - np.arange(n))


def cross_msd_fft(A, B):
    """< [A(t+tau)-A(t)] . [B(t+tau)-B(t)] >_t, A/B of shape (T, 3)."""
    T = len(A)
    tau = np.arange(T)
    out = np.zeros(T)
    for d in range(3):
        a, b = A[:, d], B[:, d]
        P = _prefix(a * b)
        S1 = (P[T - tau] + P[T] - P[tau]) / (T - tau)
        S2 = _crosscorr_fft(a, b) + _crosscorr_fft(b, a)
        out += S1 - S2
    return out


def cross_msd_brute(A, B):
    T = len(A)
    out = np.zeros(T)
    for tau in range(1, T):
        dA = A[tau:] - A[: T - tau]
        dB = B[tau:] - B[: T - tau]
        out[tau] = (dA * dB).sum(axis=1).mean()
    return out


def selftest():
    rng = np.random.default_rng(7)
    A = np.cumsum(rng.normal(size=(300, 3)), axis=0)
    B = np.cumsum(rng.normal(size=(300, 3)), axis=0) + 0.4 * A
    for X, Y, name in [(A, A, "auto"), (A, B, "cross")]:
        err = np.max(np.abs(cross_msd_fft(X, Y) - cross_msd_brute(X, Y)))
        print(f"selftest {name}: max|fft - brute| = {err:.3e}")
        assert err < 1e-8, name
    print("selftest OK")


def reduce_segment(segdir, every, refresh=False):
    """awk-reduce the pos file -> (steps[T], codes[N], traj[T, N, 3]); cached."""
    cache = os.path.join(segdir, f"reduced_traj_e{every}.npz")
    if os.path.exists(cache) and not refresh:
        z = np.load(cache)
        return z["steps"], z["codes"], z["traj"]
    pos = sorted(glob.glob(os.path.join(segdir, "*-pos-1.xyz")))
    if not pos:
        raise SystemExit(f"no *-pos-1.xyz in {segdir}")
    pos = pos[0]
    with open(pos) as fh:
        nat = int(fh.readline().strip())
    cmd = ["awk", "-v", f"NAT={nat}", "-v", f"K={every}", AWK, pos]
    raw = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    data = np.loadtxt(raw.splitlines())
    n_per = int(np.sum(data[:, 0] == data[0, 0]))   # rows in the first frame
    n_full = len(data) // n_per                     # drop a truncated tail frame
    data = data[: n_full * n_per]
    steps = data[::n_per, 0].astype(int)
    if not np.all(np.diff(steps) > 0):
        raise SystemExit(f"{segdir}: non-monotonic frames after reduction")
    codes = data[:n_per, 1].astype(int)
    traj = data[:, 2:5].reshape(n_full, n_per, 3)
    np.savez_compressed(cache, steps=steps, codes=codes, traj=traj)
    return steps, codes, traj


def segment_onsager(segdir, every, dt_dump_ps, tmax_traj_ps=None, refresh=False):
    """Per-segment self/collective MSD table -> written to onsager_cp2k.dat."""
    steps, codes, traj = reduce_segment(segdir, every, refresh)
    dt = dt_dump_ps * every
    if tmax_traj_ps is not None:
        Tcap = int(tmax_traj_ps / dt) + 1
        traj = traj[:Tcap]
    T = len(traj)
    t = np.arange(T) * dt
    cols, header = [t], ["tau_ps"]
    for name in ("O", "Na", "Cl"):
        sel = traj[:, codes == EL_CODE[name], :]
        msd = np.mean([cross_msd_fft(sel[:, i], sel[:, i])
                       for i in range(sel.shape[1])], axis=0)
        cols.append(msd)
        header.append(f"msd_{name}")
    R = {name: traj[:, codes == EL_CODE[name], :].sum(axis=1)
         for name in ("Na", "Cl")}
    for a, b in [("Na", "Na"), ("Cl", "Cl"), ("Na", "Cl")]:
        cols.append(cross_msd_fft(R[a], R[b]))
        header.append(f"coll_{a}{b}")
    out = os.path.join(segdir, "onsager_cp2k.dat")
    np.savetxt(out, np.column_stack(cols), header=" ".join(header))
    nNa = int((codes == EL_CODE["Na"]).sum())
    nCl = int((codes == EL_CODE["Cl"]).sum())
    return out, nNa, nCl, T, dt


def slope(t, y, tmin, tmax):
    m = (t >= tmin) & (t <= tmax)
    return np.polyfit(t[m], y[m], 1)[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--segdirs", nargs="+")
    ap.add_argument("--box-a", type=float, help="cubic box edge in A")
    ap.add_argument("--temp", type=float, default=300.0)
    ap.add_argument("--every", type=int, default=20,
                    help="keep every Kth dumped frame (dump=5 fs -> K=20 = 0.1 ps)")
    ap.add_argument("--dt-dump-ps", type=float, default=0.005,
                    help="time between dumped frames (10 steps x 0.5 fs)")
    ap.add_argument("--tmin", type=float, default=10.0, help="fit window lo [ps]")
    ap.add_argument("--tmax", type=float, default=40.0, help="fit window hi [ps]")
    ap.add_argument("--traj-max-ps", type=float, default=160.0,
                    help="truncate segments to a common length (overrun segments)")
    ap.add_argument("--label", default="run")
    ap.add_argument("--refresh", action="store_true", help="ignore npz caches")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return
    if not (args.segdirs and args.box_a):
        ap.error("--segdirs and --box-a required")

    V = args.box_a ** 3
    pref = E2 / (KB * args.temp * V * A3_TO_M3) * A2PS_TO_M2S  # z=1, S/m per A^2/ps

    rows = []
    for sd in args.segdirs:
        sd = sd.rstrip("/")
        out, nNa, nCl, T, dt = segment_onsager(
            sd, args.every, args.dt_dump_ps, args.traj_max_ps, args.refresh)
        ons = np.loadtxt(out)
        t = ons[:, 0]
        D = {s: slope(t, ons[:, i + 1], args.tmin, args.tmax) / 6.0
             for i, s in enumerate(["O", "Na", "Cl"])}
        sNaNa = slope(t, ons[:, 4], args.tmin, args.tmax)
        sClCl = slope(t, ons[:, 5], args.tmin, args.tmax)
        sNaCl = slope(t, ons[:, 6], args.tmin, args.tmax)
        s_sum = sNaNa + sClCl - 2.0 * sNaCl
        kNE = pref * (nNa * D["Na"] + nCl * D["Cl"])
        kOns = pref * s_sum / 6.0
        tNa = (sNaNa - sNaCl) / s_sum if s_sum else np.nan
        rows.append(dict(seg=os.path.basename(sd), nNa=nNa, nCl=nCl, T=T,
                         DO=D["O"], DNa=D["Na"], DCl=D["Cl"], kNE=kNE,
                         kOns=kOns, tNa=tNa, sNaNa=sNaNa, sClCl=sClCl,
                         sNaCl=sNaCl))
        print(f"{rows[-1]['seg']}: T={T} frames  D_O={D['O']/6*6*A2PS_TO_M2S*1e9:.3f}e-9  "
              f"kNE={kNE:.3f}  kOns={kOns:.3f} S/m  t_Na={tNa:.3f}")

    n = len(rows)
    get = lambda k: np.array([r[k] for r in rows])
    sem = lambda a: (np.std(a, ddof=1) / np.sqrt(n)) if n > 1 else np.nan
    DNa, DCl, DO = get("DNa"), get("DCl"), get("DO")
    kNE, kOns, tNa = get("kNE"), get("kOns"), get("tNa")
    dlt = 1.0 - kOns / kNE
    nNa, nCl = rows[0]["nNa"], rows[0]["nCl"]

    print(f"\n=== {args.label}: {n} segments, box {args.box_a} A, "
          f"{nNa} Na / {nCl} Cl, T={args.temp} K, fit {args.tmin}-{args.tmax} ps ===")
    conv = A2PS_TO_M2S * 1e9  # A^2/ps -> 1e-9 m^2/s
    for nm, arr in [("D_Na", DNa), ("D_Cl", DCl), ("D_w(O)", DO)]:
        print(f"{nm:7s} = {arr.mean()*conv:6.3f} +- {sem(arr)*conv:5.3f}  x1e-9 m^2/s")
    print(f"kappa_NE  = {kNE.mean():6.3f} +- {sem(kNE):5.3f}  S/m   (z=+-1)")
    print(f"kappa_Ons = {kOns.mean():6.3f} +- {sem(kOns):5.3f}  S/m")
    print(f"Delta_NE  = {dlt.mean():6.3f} +- {sem(dlt):5.3f}")
    print(f"t_Na      = {tNa.mean():6.3f} +- {sem(tNa):5.3f}")
    kNaNa = pref * get("sNaNa") / 6.0
    kClCl = pref * get("sClCl") / 6.0
    kNaCl = pref * get("sNaCl") / 6.0
    dNaNa = kNaNa - pref * nNa * DNa
    dClCl = kClCl - pref * nCl * DCl
    print("channel decomposition (S/m, z=1): "
          f"kNaNa={kNaNa.mean():.3f} kClCl={kClCl.mean():.3f} "
          f"kNaCl={kNaCl.mean():.3f} (distinct: dNaNa={dNaNa.mean():.3f} "
          f"dClCl={dClCl.mean():.3f}; -2kNaCl contributes {-2*kNaCl.mean():.3f})")

    results = results_dir(args.label)
    os.makedirs(results, exist_ok=True)
    csv = os.path.join(results, f"{args.label}_kappa.csv")
    with open(csv, "w") as fh:
        fh.write("seg,nNa,nCl,frames,D_O_A2ps,D_Na_A2ps,D_Cl_A2ps,"
                 "kNE_Sm,kOns_Sm,t_Na,sNaNa,sClCl,sNaCl\n")
        for r in rows:
            fh.write(f"{r['seg']},{r['nNa']},{r['nCl']},{r['T']},{r['DO']:.6e},"
                     f"{r['DNa']:.6e},{r['DCl']:.6e},{r['kNE']:.4f},"
                     f"{r['kOns']:.4f},{r['tNa']:.4f},{r['sNaNa']:.6e},"
                     f"{r['sClCl']:.6e},{r['sNaCl']:.6e}\n")
    np.savez(os.path.join(results, f"{args.label}_kappa.npz"),
             **{k: get(k) for k in ("DO", "DNa", "DCl", "kNE", "kOns", "tNa",
                                    "sNaNa", "sClCl", "sNaCl")},
             nNa=nNa, nCl=nCl, V=V, temp=args.temp,
             tmin=args.tmin, tmax=args.tmax)
    print(f"saved {csv} (+ .npz)")


if __name__ == "__main__":
    sys.exit(main())
