#!/usr/bin/env python3
"""Green-Kubo diffusion from a LAMMPS velocity dump.

D = (1/3) int_0^t <v(0).v(t')> dt', VACF averaged over atoms and all time
origins (FFT, exact). Outputs the normalized VACF and the running integral
D(t) per species; D is read off the plateau of D(t). Velocities are in
`real` units (A/fs); D is reported in A^2/ps.

Usage:
  vacf_analyze.py --dump L24.84_vacf.vel.lammpstrj --dt-fs 4.0 --out L24.84.vacf
"""
import argparse

import numpy as np

TYPE_NAME = {1: "O", 3: "Na", 4: "Cl"}


def read_dump(path):
    """-> types[N], vel[T, N, 3], atoms sorted by id."""
    frames, types = [], None
    with open(path) as f:
        while True:
            line = f.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            f.readline()
            assert f.readline().startswith("ITEM: NUMBER")
            n = int(f.readline())
            f.readline()
            for _ in range(3):
                f.readline()
            assert f.readline().startswith("ITEM: ATOMS")
            block = np.array([f.readline().split() for _ in range(n)], dtype=float)
            block = block[np.argsort(block[:, 0])]
            if types is None:
                types = block[:, 1].astype(int)
            frames.append(block[:, 2:5])
    return types, np.asarray(frames)


def acf_fft(x):
    """Autocorrelation over all origins for series x[T]; normalized by overlap."""
    n = len(x)
    F = np.fft.fft(x, 2 * n)
    r = np.fft.ifft(F * F.conjugate())[:n].real
    return r / (n - np.arange(n))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True)
    ap.add_argument("--dt-fs", type=float, default=4.0, help="fs between frames")
    ap.add_argument("--tmax-ps", type=float, default=5.0, help="truncate VACF here")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    types, vel = read_dump(args.dump)
    T = len(vel)
    nlag = min(T, int(args.tmax_ps * 1000 / args.dt_fs))
    t_ps = np.arange(nlag) * args.dt_fs / 1000.0

    cols, names = [t_ps], ["t_ps"]
    print(f"{T} frames, dt={args.dt_fs} fs, VACF to {t_ps[-1]:.2f} ps")
    for ty, name in sorted(TYPE_NAME.items()):
        idx = np.where(types == ty)[0]
        if len(idx) == 0:
            continue
        # <v(0).v(t)> averaged over atoms, xyz, and origins  [A^2/fs^2]
        vacf = np.zeros(nlag)
        for i in idx:
            for k in range(3):
                vacf += acf_fft(vel[:, i, k])[:nlag]
        vacf /= len(idx)
        # D(t) = (1/3) int_0^t vacf dt'   [A^2/fs] -> *1000 = A^2/ps
        Dt = np.concatenate([[0.0], np.cumsum((vacf[1:] + vacf[:-1]) / 2)]) \
            * args.dt_fs / 3.0 * 1000.0
        cols += [vacf / vacf[0], Dt]
        names += [f"vacf_{name}", f"D_{name}"]
        # plateau estimate: average D(t) over the last 40% of the window
        plat = Dt[int(0.6 * nlag):].mean()
        print(f"  {name:>2}: {len(idx):4d} atoms  D_plateau = {plat:.4f} A^2/ps")

    np.savetxt(args.out, np.column_stack(cols), header=" ".join(names))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
