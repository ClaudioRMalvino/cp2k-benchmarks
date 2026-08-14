#!/usr/bin/env python3
"""Green-Kubo shear viscosity from CP2K per-step stress files.

Daivis-Evans symmetric-traceless estimator: eta = V/(10 kB T) int sum_ab
<P^os_ab(0) P^os_ab(t)> dt; ACFs via FFT, plateau = mean of the running integral
over [plateau-lo, plateau-hi]; segments -> mean +- SEM. Stress rows monotonic-
sliced on step (walltime resumes duplicate up to one checkpoint interval).
Writes {label}_eta.npz to results/{mp2,rpa}_transport/.
Usage: --segdirs ... --box-a 37.26 --label X [--plateau-lo 2 --plateau-hi 10]
"""
import argparse
import glob
import os
import sys

import numpy as np

KB = 1.380649e-23
BAR = 1.0e5                      # bar -> Pa
RESULTS_ROOT = os.path.expanduser("~/cp2k-benchmarks/results")


def results_dir(label):
    """Level-separated outputs: mp2* labels go to results/mp2_transport/,
    everything else (rpa_*) to results/rpa_transport/."""
    sub = "mp2_transport" if label.startswith("mp2") else "rpa_transport"
    return os.path.join(RESULTS_ROOT, sub)


def acf_fft(a):
    """<a(0) a(t)> with (n - t) normalisation, via FFT."""
    n = len(a)
    F = np.fft.fft(a, 2 * n)
    res = np.fft.ifft(F.conjugate() * F)[:n].real
    return res / (n - np.arange(n))


def segment_eta(segdir, V_m3, temp, dt_fs, acf_max_ps):
    st = sorted(glob.glob(os.path.join(segdir, "*-1.stress")))
    if not st:
        raise SystemExit(f"no *-1.stress in {segdir}")
    d = np.loadtxt(st[0], comments="#")
    # monotonic slice on the step column
    keep = np.zeros(len(d), dtype=bool)
    seen = -1
    for i, s in enumerate(d[:, 0].astype(int)):
        if s > seen:
            keep[i] = True
            seen = s
    d = d[keep]
    P = d[:, 2:11].reshape(-1, 3, 3) * BAR          # Pa
    Ps = 0.5 * (P + np.transpose(P, (0, 2, 1)))     # symmetrize
    tr = np.trace(Ps, axis1=1, axis2=2) / 3.0
    for k in range(3):
        Ps[:, k, k] -= tr                            # traceless
    nacf = int(acf_max_ps * 1000.0 / dt_fs)
    C = np.zeros(nacf)
    for a in range(3):
        for b in range(3):
            C += acf_fft(Ps[:, a, b])[:nacf]
    dt_s = dt_fs * 1e-15
    t_ps = np.arange(nacf) * dt_fs / 1000.0
    run = V_m3 / (10.0 * KB * temp) * np.concatenate(
        ([0.0], np.cumsum(0.5 * (C[1:] + C[:-1]) * dt_s)))
    return t_ps, run, len(d)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--segdirs", nargs="+", required=True)
    ap.add_argument("--box-a", type=float, required=True)
    ap.add_argument("--temp", type=float, default=300.0)
    ap.add_argument("--dt-fs", type=float, default=0.5)
    ap.add_argument("--acf-max-ps", type=float, default=20.0)
    ap.add_argument("--plateau-lo", type=float, default=2.0)
    ap.add_argument("--plateau-hi", type=float, default=10.0)
    ap.add_argument("--label", default="run")
    args = ap.parse_args()

    V = (args.box_a * 1e-10) ** 3
    etas, curves = [], []
    for sd in args.segdirs:
        sd = sd.rstrip("/")
        t, run, nrow = segment_eta(sd, V, args.temp, args.dt_fs, args.acf_max_ps)
        m = (t >= args.plateau_lo) & (t <= args.plateau_hi)
        eta = run[m].mean()
        etas.append(eta)
        curves.append(run)
        print(f"{os.path.basename(sd)}: {nrow} stress rows, "
              f"eta_plateau = {eta*1e3:.3f} mPa s")

    etas = np.array(etas)
    n = len(etas)
    sem = np.std(etas, ddof=1) / np.sqrt(n) if n > 1 else np.nan
    print(f"\n=== {args.label}: {n} segments, box {args.box_a} A, "
          f"T={args.temp} K, plateau {args.plateau_lo}-{args.plateau_hi} ps ===")
    print(f"eta = {etas.mean()*1e3:.3f} +- {sem*1e3:.3f} mPa s")

    results = results_dir(args.label)
    os.makedirs(results, exist_ok=True)
    np.savez(os.path.join(results, f"{args.label}_eta.npz"),
             t_ps=t, curves=np.array(curves), etas=etas,
             plateau=(args.plateau_lo, args.plateau_hi), temp=args.temp,
             box_a=args.box_a)
    print(f"saved {results}/{args.label}_eta.npz")


if __name__ == "__main__":
    sys.exit(main())
