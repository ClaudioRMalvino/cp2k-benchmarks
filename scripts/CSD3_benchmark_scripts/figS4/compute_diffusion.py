#!/usr/bin/env python3
"""Self-diffusion coefficient from a CP2K MD position trajectory with
Yeh-Hummer finite-size correction (Fig. S4 D of Morawietz et al., PNAS 113
(2016) 8368).

For non-cubic tiled cells (128/256/1024) the Yeh-Hummer L is taken as V^(1/3);
xi=2.837297 is strictly the cubic value, so D_0 at those sizes carries a small
box-shape approximation.
"""
import argparse
import re
import sys

import numpy as np

KB = 1.380649e-23
XI = 2.837297

_TIME_RE = re.compile(r"time\s*=\s*([-\d.Ee+]+)")


def parse_xyz(path, species):
    times, frames = [], []
    with open(path) as f:
        while True:
            header = f.readline()
            if not header:
                break
            try:
                n_atoms = int(header.split()[0])
            except (ValueError, IndexError):
                break
            comment = f.readline()
            m = _TIME_RE.search(comment)
            times.append(float(m.group(1)) if m else None)
            pos = []
            for _ in range(n_atoms):
                p = f.readline().split()
                if p[0] == species:
                    pos.append((float(p[1]), float(p[2]), float(p[3])))
            frames.append(pos)
    if not frames:
        raise RuntimeError(f"no frames parsed from {path}")
    times = None if any(t is None for t in times) else np.asarray(times)
    return times, np.asarray(frames, dtype=np.float64)


def unwrap(pos, cell):
    d = np.diff(pos, axis=0)
    d -= cell * np.round(d / cell)
    out = np.empty_like(pos)
    out[0] = pos[0]
    out[1:] = pos[0] + np.cumsum(d, axis=0)
    return out


def msd_fft(r):
    """FFT MSD estimator for one atom, r shape (nframes,3)."""
    n = len(r)
    d = np.sum(r ** 2, axis=1)
    d = np.append(d, 0.0)
    s2 = np.zeros(n)
    for ax in range(3):
        x = r[:, ax]
        fsize = 1
        while fsize < 2 * n:
            fsize *= 2
        fx = np.fft.rfft(x, n=fsize)
        s2 += np.fft.irfft(fx * np.conj(fx))[:n]
    q = 2.0 * d.sum()
    s1 = np.zeros(n)
    for m in range(n):
        q -= (d[m - 1] + d[n - m]) if m > 0 else 0.0
        s1[m] = q / (n - m)
    return s1 - 2.0 * s2 / (n - np.arange(n))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("trajectory")
    ap.add_argument("--cell-ang", type=float, nargs=3, required=True,
                    metavar=("AX", "AY", "AZ"))
    ap.add_argument("--viscosity-mpa-s", type=float, required=True,
                    help="Green-Kubo viscosity of this segment, for Yeh-Hummer")
    ap.add_argument("--species", default="O")
    ap.add_argument("--timestep-fs", type=float, default=0.5)
    ap.add_argument("--traj-every", type=int, default=10,
                    help="MD steps between trajectory frames")
    ap.add_argument("--fit-frac", type=float, nargs=2, default=(0.2, 0.8),
                    metavar=("LO", "HI"),
                    help="MSD lag window used for the linear fit, as fractions")
    ap.add_argument("--out", default="diffusion.csv")
    args = ap.parse_args()

    cell = np.asarray(args.cell_ang)
    times_fs, pos = parse_xyz(args.trajectory, args.species)
    nframes, natoms, _ = pos.shape

    if times_fs is not None and nframes > 1:
        dt_ps = (times_fs[1] - times_fs[0]) / 1000.0
    else:
        dt_ps = args.timestep_fs * args.traj_every / 1000.0

    pos = unwrap(pos, cell)
    msd = np.mean([msd_fft(pos[:, i, :]) for i in range(natoms)], axis=0)
    lag_ps = np.arange(nframes) * dt_ps

    lo = int(args.fit_frac[0] * nframes)
    hi = int(args.fit_frac[1] * nframes)
    slope, intercept = np.polyfit(lag_ps[lo:hi], msd[lo:hi], 1)
    d_pbc = slope / 6.0

    L_m = (np.prod(cell)) ** (1.0 / 3.0) * 1.0e-10
    eta = args.viscosity_mpa_s * 1.0e-3
    kT = KB * 300.0
    correction_si = kT * XI / (6.0 * np.pi * eta * L_m)
    correction = correction_si * 1.0e8
    d_0 = d_pbc + correction

    with open(args.out, "w") as f:
        f.write(f"# diffusion from {args.trajectory}\n")
        f.write(f"# species: {args.species}   n_atoms: {natoms}   n_frames: {nframes}\n")
        f.write(f"# cell_ang: {cell[0]:.4f} {cell[1]:.4f} {cell[2]:.4f}\n")
        f.write(f"# dt_ps: {dt_ps:.5f}\n")
        f.write(f"# fit_window_ps: {lag_ps[lo]:.3f} {lag_ps[hi-1]:.3f}\n")
        f.write(f"# viscosity_mPa_s: {args.viscosity_mpa_s:.6f}\n")
        f.write(f"# D_PBC_ang2_ps: {d_pbc:.8f}\n")
        f.write(f"# YH_correction_ang2_ps: {correction:.8f}\n")
        f.write(f"# D_0_ang2_ps: {d_0:.8f}\n")
        f.write("lag_ps,msd_ang2\n")
        for lp, mv in zip(lag_ps, msd):
            f.write(f"{lp:.5f},{mv:.8f}\n")

    print(f"{args.trajectory}: D_PBC={d_pbc:.5f}  +YH={correction:.5f}  "
          f"D_0={d_0:.5f} Ang^2/ps  -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
