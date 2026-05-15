#!/usr/bin/env python3
"""
Compute O-O / O-H / H-H radial distribution functions from a CP2K NVT
trajectory of liquid water in a cubic periodic box.  Produces a CSV with
columns: r_angstrom, gOO, gOH, gHH — directly plottable to replicate Fig 2b
of Iannuzzi et al., J. Phys. Chem. B 2026, 130, 1237.

Usage:
    python3 compute_rdf.py <trajectory.xyz> --box <length_angstrom> --out <rdf.csv>

Single-frame XYZ format (CP2K default):
    <n_atoms>
    i = <step>, ... <comment>
    <element> <x> <y> <z>
    ...

Distances use minimum-image convention.  Bins are 0..6 Å in 200 steps to match
the resolution of the published figure.
"""
import argparse
import sys

import numpy as np


def parse_xyz_frames(path):
    """Generator yielding (atom_symbols, positions) for each frame.  Reads in a
    streaming fashion so trajectories of any length fit in memory."""
    with open(path, "r") as f:
        while True:
            header = f.readline()
            if not header:
                return
            try:
                n_atoms = int(header.strip())
            except ValueError:
                return
            f.readline()                                       # comment line
            symbols = []
            positions = np.empty((n_atoms, 3), dtype=np.float64)
            for i in range(n_atoms):
                parts = f.readline().split()
                symbols.append(parts[0])
                positions[i] = (float(parts[1]), float(parts[2]), float(parts[3]))
            yield symbols, positions


def pair_distances_pbc(pos_a, pos_b, box, same_set):
    """All pairwise distances between two atom subsets under cubic PBC.
    same_set=True drops self-pairs and double-counting (i<j only)."""
    dr = pos_a[:, None, :] - pos_b[None, :, :]
    dr -= box * np.round(dr / box)
    d = np.sqrt(np.sum(dr * dr, axis=-1))
    if same_set:
        n = d.shape[0]
        iu = np.triu_indices(n, k=1)
        return d[iu]
    return d.ravel()


def accumulate_rdf(trajectory_path, box, n_bins=200, r_max=6.0):
    edges = np.linspace(0.0, r_max, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    dr = r_max / n_bins

    hist_oo = np.zeros(n_bins, dtype=np.float64)
    hist_oh = np.zeros(n_bins, dtype=np.float64)
    hist_hh = np.zeros(n_bins, dtype=np.float64)

    n_frames = 0
    n_O = n_H = 0
    for symbols, positions in parse_xyz_frames(trajectory_path):
        symbols = np.array(symbols)
        idx_O = np.where(symbols == "O")[0]
        idx_H = np.where(symbols == "H")[0]
        if n_frames == 0:
            n_O, n_H = len(idx_O), len(idx_H)

        # O-O (same set)
        d = pair_distances_pbc(positions[idx_O], positions[idx_O], box, same_set=True)
        h, _ = np.histogram(d[d < r_max], bins=edges)
        hist_oo += h

        # O-H (different sets)
        d = pair_distances_pbc(positions[idx_O], positions[idx_H], box, same_set=False)
        h, _ = np.histogram(d[d < r_max], bins=edges)
        hist_oh += h

        # H-H (same set)
        d = pair_distances_pbc(positions[idx_H], positions[idx_H], box, same_set=True)
        h, _ = np.histogram(d[d < r_max], bins=edges)
        hist_hh += h

        n_frames += 1

    if n_frames == 0:
        raise RuntimeError(f"No frames found in {trajectory_path}")

    volume = box ** 3
    shell_vol = 4.0 * np.pi * centers ** 2 * dr

    # g_AB(r) = N_pairs(r) / [n_frames * shell_vol * rho_B * (effective N_A)]
    # For same-species AA (i<j only): rho_A * (N_A - 1) / 2 effective pair count per atom
    rho_O = n_O / volume
    rho_H = n_H / volume
    g_oo = hist_oo / (n_frames * shell_vol * rho_O * (n_O - 1) / 2.0)
    g_oh = hist_oh / (n_frames * shell_vol * rho_H * n_O)
    g_hh = hist_hh / (n_frames * shell_vol * rho_H * (n_H - 1) / 2.0)

    return centers, g_oo, g_oh, g_hh, n_frames, n_O, n_H


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("trajectory", help="CP2K *-pos-1.xyz file")
    ap.add_argument("--box", type=float, required=True,
                    help="cubic cell length in Angstrom")
    ap.add_argument("--out", default="rdf.csv", help="output CSV path")
    ap.add_argument("--bins", type=int, default=200)
    ap.add_argument("--rmax", type=float, default=6.0)
    args = ap.parse_args()

    centers, g_oo, g_oh, g_hh, n_frames, n_O, n_H = accumulate_rdf(
        args.trajectory, args.box, n_bins=args.bins, r_max=args.rmax
    )

    with open(args.out, "w") as f:
        f.write(f"# RDF from {args.trajectory}\n")
        f.write(f"# box = {args.box} A,  n_frames = {n_frames},  N_O = {n_O},  N_H = {n_H}\n")
        f.write("r_angstrom,gOO,gOH,gHH\n")
        for r, oo, oh, hh in zip(centers, g_oo, g_oh, g_hh):
            f.write(f"{r:.4f},{oo:.6f},{oh:.6f},{hh:.6f}\n")
    print(f"Wrote {args.out} ({n_frames} frames, {n_O} O, {n_H} H)")


if __name__ == "__main__":
    sys.exit(main())
