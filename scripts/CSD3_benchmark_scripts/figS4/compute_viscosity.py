#!/usr/bin/env python3
"""Green-Kubo shear viscosity from a CP2K MD .stress file (Morawietz et al.,
PNAS 113 (2016) 8368, SI section B / Fig. S4 A,B,C).

ACF is averaged over the five independent off-diagonal / traceless components
(Pxy, Pxz, Pyz, (Pxx-Pyy)/2, (Pyy-Pzz)/2); GK integral truncated at
--int-limit ps.
"""
import argparse
import sys

import numpy as np

KB = 1.380649e-23
BAR_TO_PA = 1.0e5
ANG3_TO_M3 = 1.0e-30
VOL_PER_MOLECULE_ANG3 = 12.42 ** 3 / 64.0


def read_stress(path):
    steps, times, comps = [], [], []
    with open(path) as f:
        for line in f:
            if line.lstrip().startswith("#") or not line.strip():
                continue
            p = line.split()
            if len(p) < 11:
                continue
            steps.append(int(p[0]))
            times.append(float(p[1]))
            comps.append([float(x) for x in p[2:11]])
    if not comps:
        raise RuntimeError(f"no stress data found in {path}")
    t = np.asarray(times)
    P = np.asarray(comps).reshape(-1, 3, 3)
    return t, P


def fft_acf(x, max_lag):
    """Unbiased autocorrelation of a 1-D series via FFT. No mean subtraction:
    the off-diagonal stress already fluctuates about 0."""
    n = len(x)
    fsize = 1
    while fsize < 2 * n:
        fsize *= 2
    fx = np.fft.rfft(x, n=fsize)
    acf = np.fft.irfft(fx * np.conj(fx))[:n]
    acf /= (n - np.arange(n))
    return acf[:max_lag]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stress_file")
    ap.add_argument("--n-molecules", type=int, required=True)
    ap.add_argument("--volume-ang3", type=float, default=None,
                    help="cell volume; default = n_molecules * (12.42^3/64)")
    ap.add_argument("--temperature", type=float, default=300.0)
    ap.add_argument("--timestep-fs", type=float, default=0.5,
                    help="spacing between stress samples (stress printed every step)")
    ap.add_argument("--discard-ps", type=float, default=10.0,
                    help="initial transient discarded before the ACF")
    ap.add_argument("--int-limit-ps", type=float, default=3.0,
                    help="Green-Kubo integral upper limit (paper: 3 ps)")
    ap.add_argument("--out", default="viscosity.csv")
    args = ap.parse_args()

    V_ang3 = args.volume_ang3 or args.n_molecules * VOL_PER_MOLECULE_ANG3
    V = V_ang3 * ANG3_TO_M3
    kT = KB * args.temperature

    t_fs, P = read_stress(args.stress_file)
    dt_fs = args.timestep_fs
    n_discard = int(round(args.discard_ps * 1000.0 / dt_fs))
    P = P[n_discard:]
    if len(P) < 100:
        raise RuntimeError(f"too few stress frames after discarding "
                            f"{args.discard_ps} ps ({len(P)} left)")

    comps = np.array([
        P[:, 0, 1],
        P[:, 0, 2],
        P[:, 1, 2],
        0.5 * (P[:, 0, 0] - P[:, 1, 1]),
        0.5 * (P[:, 1, 1] - P[:, 2, 2]),
    ]) * BAR_TO_PA

    max_lag = int(round(args.int_limit_ps * 1000.0 / dt_fs)) + 1
    max_lag = min(max_lag, len(P) // 2)

    acf = np.mean([fft_acf(c, max_lag) for c in comps], axis=0)

    lag_fs = np.arange(max_lag) * dt_fs
    dt_s = dt_fs * 1e-15
    running_integral = np.concatenate([[0.0],
                                       np.cumsum(0.5 * (acf[1:] + acf[:-1]) * dt_s)])
    running_eta = (V / kT) * running_integral * 1.0e3

    eta_final = running_eta[-1]
    acf_norm = acf / acf[0]

    with open(args.out, "w") as f:
        f.write(f"# Green-Kubo viscosity from {args.stress_file}\n")
        f.write(f"# n_molecules: {args.n_molecules}\n")
        f.write(f"# volume_ang3: {V_ang3:.4f}\n")
        f.write(f"# temperature_K: {args.temperature}\n")
        f.write(f"# timestep_fs: {dt_fs}\n")
        f.write(f"# discard_ps: {args.discard_ps}\n")
        f.write(f"# int_limit_ps: {args.int_limit_ps}\n")
        f.write(f"# n_frames_used: {len(P)}\n")
        f.write(f"# eta_mPa_s: {eta_final:.6f}\n")
        f.write("lag_fs,acf_normalized,running_eta_mPa_s\n")
        for lf, an, re_ in zip(lag_fs, acf_norm, running_eta):
            f.write(f"{lf:.3f},{an:.8f},{re_:.8f}\n")

    print(f"{args.stress_file}: n={args.n_molecules}  frames={len(P)}  "
          f"eta({args.int_limit_ps}ps) = {eta_final:.4f} mPa.s  -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
