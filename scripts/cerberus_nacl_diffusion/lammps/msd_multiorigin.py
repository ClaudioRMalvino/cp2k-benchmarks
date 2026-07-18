#!/usr/bin/env python3
"""Multi-time-origin MSD from a LAMMPS unwrapped-coordinate dump.

`compute msd` in LAMMPS measures displacement from one fixed origin, so a 2 ns
run gives a single sample per atom. Diffusion is time-translation invariant, so
we can instead average |r(t0+tau) - r(t0)|^2 over every origin t0. That is ~1000x
more samples for free -- decisive for the ions, where we only have 9-48 of them.

Done naively this is O(T^2) per atom; we use the FFT method of Calandrini et al.
(Collection SFN 12, 201 (2011)), which is O(T log T) and exact, not approximate.

Usage:
  msd_multiorigin.py --dump L37.26.pos.lammpstrj --dt-ps 1.0 --out L37.26.msdmo
"""
import argparse

import numpy as np

TYPE_NAME = {1: "O", 3: "Na", 4: "Cl"}


def read_dump(path):
    """-> types[N], pos[T, N, 3] of unwrapped coords, atoms sorted by id."""
    frames, types = [], None
    with open(path) as f:
        while True:
            line = f.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            f.readline()                                  # timestep value
            assert f.readline().startswith("ITEM: NUMBER")
            n = int(f.readline())
            f.readline()                                  # ITEM: BOX BOUNDS
            for _ in range(3):
                f.readline()
            assert f.readline().startswith("ITEM: ATOMS")
            block = np.array([f.readline().split() for _ in range(n)], dtype=float)
            block = block[np.argsort(block[:, 0])]        # by id
            if types is None:
                types = block[:, 1].astype(int)
            frames.append(block[:, 2:5])
    return types, np.asarray(frames)


def _autocorr_fft(x):
    """Autocorrelation of a 1-D series via FFT, normalised by overlap count."""
    n = len(x)
    F = np.fft.fft(x, 2 * n)
    res = np.fft.ifft(F * F.conjugate())[:n].real
    return res / (n - np.arange(n))


def msd_fft(r):
    """MSD(tau) averaged over all time origins. r: (T,3) -> (T,)"""
    T = len(r)
    sq = np.square(r).sum(axis=1)
    s2 = sum(_autocorr_fft(r[:, i]) for i in range(3))
    padded = np.append(sq, 0.0)
    q = 2.0 * sq.sum()
    s1 = np.zeros(T)
    for m in range(T):
        q -= padded[m - 1] + padded[T - m]
        s1[m] = q / (T - m)
    return s1 - 2.0 * s2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True)
    ap.add_argument("--dt-ps", type=float, default=1.0, help="ps between frames")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    types, pos = read_dump(args.dump)
    T = len(pos)
    tau = np.arange(T) * args.dt_ps

    cols, names = [tau], ["tau_ps"]
    for t, name in sorted(TYPE_NAME.items()):
        sel = types == t
        if not sel.any():
            continue
        m = np.mean([msd_fft(pos[:, i, :]) for i in np.where(sel)[0]], axis=0)
        cols.append(m)
        names.append(f"msd_{name}")
        print(f"{name}: {sel.sum()} atoms, {T} frames, "
              f"{T} origins at tau=0 -> {1} at tau_max")

    np.savetxt(args.out, np.column_stack(cols), header=" ".join(names))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
