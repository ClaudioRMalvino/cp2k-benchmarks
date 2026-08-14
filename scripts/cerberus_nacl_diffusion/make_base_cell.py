#!/usr/bin/env python3
"""Build the NaCl(aq) base simulation cell for the diffusion study.

"substitute" (default): first 64-water frame of the O'Neill training data
(data-sets/<model>/input_water.xyz, 12.42 A), remove 2*n_pairs waters at
maximally separated O sites, place Na+/Cl- there (n_pairs=1: 188 atoms, ~0.90
mol/kg). "extract": first config_type=NaCl frame of input.xyz verbatim.
Writes a plain "El x y z" file for @INCLUDE in &COORD. Equilibrate before
production: substituted ions start at water O sites with no solvation shell.
"""
import argparse
import sys

import numpy as np

BASE_L = 12.42


def read_first_frame(path, require_na=False):
    """Return (symbols, coords) of the first frame (or first frame with Na)."""
    with open(path) as f:
        while True:
            line = f.readline()
            if not line:
                raise SystemExit(f"no suitable frame found in {path}")
            try:
                nat = int(line.strip())
            except ValueError:
                continue
            f.readline()  # comment / Lattice line
            syms, xyz = [], []
            for _ in range(nat):
                p = f.readline().split()
                syms.append(p[0])
                xyz.append([float(p[1]), float(p[2]), float(p[3])])
            if not require_na or "Na" in syms:
                return np.array(syms), np.array(xyz)


def min_image(d, L):
    return d - L * np.round(d / L)


def farthest_point_sites(opos, k, L):
    """Pick k oxygen indices that are mutually well separated (greedy FPS)."""
    chosen = [0]
    while len(chosen) < k:
        dmin = np.full(len(opos), np.inf)
        for c in chosen:
            d = np.linalg.norm(min_image(opos - opos[c], L), axis=1)
            dmin = np.minimum(dmin, d)
        dmin[chosen] = -1.0
        chosen.append(int(np.argmax(dmin)))
    return chosen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", required=True,
                    help=".../nacl-water/data-sets directory")
    ap.add_argument("--model", default="revPBE-D3")
    ap.add_argument("--mode", choices=["substitute", "extract"],
                    default="substitute")
    ap.add_argument("--n-pairs", type=int, default=1,
                    help="NaCl pairs per cell (substitute mode); ignored when "
                         "--target-molality is set")
    ap.add_argument("--cubic-mult", type=int, default=1,
                    help="tile the 12.42 A water frame N x N x N into a CUBIC "
                         "supercell of side 12.42*N before substituting ions "
                         "(N=2 -> 24.84 A, N=3 -> 37.26 A). Keeps the box cubic "
                         "so the Yeh-Hummer xi=2.837297 constant stays valid.")
    ap.add_argument("--target-molality", type=float, default=None,
                    help="if set (e.g. 1.0), pick the integer pair count that "
                         "best matches this molality in mol/kg for the chosen "
                         "box, overriding --n-pairs")
    ap.add_argument("--water-frame", default=None,
                    help="explicit pure-water .xyz seed frame; overrides "
                         "<datasets>/<model>/input_water.xyz. Use when a model "
                         "(e.g. MP2) ships no water frame - the seed geometry is "
                         "erased by equilibration, so another method's frame is "
                         "fine.")
    ap.add_argument("--out", required=True,
                    help="coordinate file to write (CP2K @INCLUDE format)")
    args = ap.parse_args()

    water_frame = args.water_frame or \
        f"{args.datasets}/{args.model}/input_water.xyz"

    # cubic supercell side used for tiling / wrapping / ion placement
    cell_L = BASE_L * args.cubic_mult

    if args.mode == "extract":
        syms, xyz = read_first_frame(
            f"{args.datasets}/{args.model}/input.xyz", require_na=True)
    elif args.n_pairs == 0 and args.target_molality is None:
        # pure-water reference cell (needed for D(c)/D(c=0) ratios)
        syms, xyz = read_first_frame(water_frame)
    else:
        syms, xyz = read_first_frame(water_frame)
        # wrap the base frame into [0, BASE_L) then tile N x N x N into a cubic
        # supercell so bigger cubic boxes keep the Yeh-Hummer xi valid
        xyz = xyz - BASE_L * np.floor(xyz / BASE_L)
        if args.cubic_mult > 1:
            N = args.cubic_mult
            shifts = np.array([[i, j, k] for i in range(N)
                               for j in range(N) for k in range(N)]) * BASE_L
            xyz = np.vstack([xyz + s for s in shifts])
            syms = np.tile(syms, N ** 3)
        # decide pair count: match a target molality if requested
        n_w0 = int(np.sum(syms == "O"))
        if args.target_molality is not None:
            a = args.target_molality * 18.015e-3
            n_pairs = int(round(a * n_w0 / (1.0 + 2.0 * a)))
        else:
            n_pairs = args.n_pairs
        omask = syms == "O"
        oidx = np.where(omask)[0]
        opos = xyz[omask]
        # group each O with its 2 nearest H (minimum image)
        hidx = np.where(syms == "H")[0]
        hpos = xyz[hidx]
        mol_h = {}
        for o_i, op in zip(oidx, opos):
            d = np.linalg.norm(min_image(hpos - op, cell_L), axis=1)
            mol_h[o_i] = hidx[np.argsort(d)[:2]]
        sites = farthest_point_sites(opos, 2 * n_pairs, cell_L)
        remove = set()
        ions = []
        for j, s in enumerate(sites):
            o_i = oidx[s]
            remove.add(o_i)
            remove.update(mol_h[o_i].tolist())
            ions.append(("Na" if j % 2 == 0 else "Cl", xyz[o_i]))
        keep = [i for i in range(len(syms)) if i not in remove]
        syms = np.concatenate([[el for el, _ in ions],
                               syms[keep]])
        xyz = np.vstack([[pos for _, pos in ions], xyz[keep]])

    n_na = int(np.sum(syms == "Na"))
    n_cl = int(np.sum(syms == "Cl"))
    n_w = int(np.sum(syms == "O"))
    if n_na != n_cl:
        sys.exit(f"unbalanced ions: {n_na} Na vs {n_cl} Cl")
    molality = n_na / (n_w * 18.015e-3) if n_w else 0.0

    # wrap into [0, L) so MULTIPLE_UNIT_CELL tiling replicates a clean cell
    xyz = xyz - cell_L * np.floor(xyz / cell_L)

    with open(args.out, "w") as f:
        for s, (x, y, z) in zip(syms, xyz):
            f.write(f"{s:2s} {x:18.10f} {y:18.10f} {z:18.10f}\n")

    print(f"wrote {args.out}: {len(syms)} atoms "
          f"({n_na} Na, {n_cl} Cl, {n_w} H2O) in {cell_L:.2f} A "
          f"{'cubic ' if args.cubic_mult > 1 else ''}cell "
          f"-> {molality:.3f} mol/kg")
    return 0


if __name__ == "__main__":
    sys.exit(main())
