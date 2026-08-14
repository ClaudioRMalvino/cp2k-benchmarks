#!/usr/bin/env python3
"""Multi-time-origin self AND collective (Onsager) MSDs from a LAMMPS dump.

Self MSDs -> D_i (Nernst-Einstein). Collective MSDs of R_s(t) = sum_i r_i(t)
give the Onsager coefficients via Einstein-Helfand:
    C_ss'(tau) = <[R_s(t+tau)-R_s(t)] . [R_s'(t+tau)-R_s'(t)]>_t
kappa = e^2 z^2 / (6 kB T V) * d/dt (C_NaNa + C_ClCl - 2 C_NaCl), z applied at
analysis time (the correlations are charge-free; z = 0.85 model, 1 formal).
FFT time-origin trick (Calandrini et al. 2011) extended to cross-correlations;
verify with --selftest.
"""
import argparse
import sys

import numpy as np

TYPE_NAME = {1: "O", 3: "Na", 4: "Cl"}


def read_dump(path):
    """Parse ITEM: blocks -> (types, traj[T, natoms, 3]) sorted by atom id."""
    frames, types = [], None
    with open(path) as fh:
        while True:
            line = fh.readline()
            if not line:
                break
            if line.startswith("ITEM: NUMBER OF ATOMS"):
                n = int(fh.readline())
            elif line.startswith("ITEM: ATOMS"):
                block = np.loadtxt((fh.readline() for _ in range(n)))
                block = block[np.argsort(block[:, 0])]
                if types is None:
                    types = block[:, 1].astype(int)
                frames.append(block[:, 2:5])
    return types, np.array(frames)


def _prefix(p):
    P = np.zeros(len(p) + 1)
    P[1:] = np.cumsum(p)
    return P


def _crosscorr_fft(a, b):
    """C(tau) = (1/(n-tau)) sum_t a(t) b(t+tau), via FFT."""
    n = len(a)
    Fa = np.fft.fft(a, 2 * n)
    Fb = np.fft.fft(b, 2 * n)
    res = np.fft.ifft(Fa.conjugate() * Fb)[:n].real
    return res / (n - np.arange(n))


def cross_msd_fft(A, B):
    """< [A(t+tau)-A(t)] . [B(t+tau)-B(t)] >_t, A/B (T, 3), via the S1 - S2
    decomposition (prefix sums + two FFT cross-correlations); A=B = standard FFT MSD."""
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
    for tau in range(T):
        dA = A[tau:] - A[: T - tau] if tau else np.zeros_like(A)
        dB = B[tau:] - B[: T - tau] if tau else np.zeros_like(B)
        out[tau] = (dA * dB).sum(axis=1).mean() if tau else 0.0
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump")
    ap.add_argument("--dt-ps", type=float, default=1.0)
    ap.add_argument("--out")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return
    if not (args.dump and args.out):
        ap.error("--dump and --out required")

    types, traj = read_dump(args.dump)
    T = len(traj)
    t = np.arange(T) * args.dt_ps
    cols, header = [t], ["tau_ps"]

    # self MSDs (per-atom average), same as msd_multiorigin
    for tp, name in TYPE_NAME.items():
        sel = traj[:, types == tp, :]
        msd = np.mean([cross_msd_fft(sel[:, i], sel[:, i])
                       for i in range(sel.shape[1])], axis=0)
        cols.append(msd)
        header.append(f"msd_{name}")

    # collective displacement correlations for the ion Onsager block
    R = {name: traj[:, types == tp, :].sum(axis=1)
         for tp, name in TYPE_NAME.items() if name in ("Na", "Cl")}
    for a, b in [("Na", "Na"), ("Cl", "Cl"), ("Na", "Cl")]:
        cols.append(cross_msd_fft(R[a], R[b]))
        header.append(f"coll_{a}{b}")

    np.savetxt(args.out, np.column_stack(cols), header=" ".join(header))
    nNa = int((types == 3).sum())
    nCl = int((types == 4).sum())
    print(f"wrote {args.out}  (T={T} frames, {nNa} Na, {nCl} Cl)")


if __name__ == "__main__":
    sys.exit(main())
