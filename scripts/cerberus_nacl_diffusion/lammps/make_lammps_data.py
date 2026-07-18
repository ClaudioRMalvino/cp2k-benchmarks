#!/usr/bin/env python3
"""Build a LAMMPS data file for NaCl(aq) with the Madrid-2019 force field
(TIP4P/2005 water + scaled-charge ions).

Generates a cubic box of ideal-geometry TIP4P/2005 water molecules on a jittered
grid at ~0.997 g/cm^3, then substitutes n_pairs waters with Na+/Cl-. The M-site
is implicit (LAMMPS lj/cut/tip4p/long computes it), so each water is 3 sites
(O, H, H) and the O atom carries the M charge -1.1128 e.

Atom types: 1=O, 2=H, 3=Na, 4=Cl.  Bond type 1 = O-H.  Angle type 1 = H-O-H.

Usage:
  make_lammps_data.py --side 24.84 --target-molality 1.0 --out data.nacl
  make_lammps_data.py --side 34.0  --target-molality 1.0 --out data.big
"""
import argparse
import numpy as np

# TIP4P/2005 geometry
DOH = 0.9572
ANG = np.deg2rad(104.52)
# charges (M charge sits on O for LAMMPS TIP4P)
Q = {"O": -1.1128, "H": 0.5564, "Na": 0.85, "Cl": -0.85}
MASS = {"O": 15.9994, "H": 1.008, "Na": 22.98977, "Cl": 35.453}
TYPE = {"O": 1, "H": 2, "Na": 3, "Cl": 4}
RHO = 0.03334          # water molecules per A^3 (~0.997 g/cm^3)


def water(o, rng):
    """Return (O,H,H) coords for an ideal TIP4P/2005 water at O with random rot."""
    h1 = np.array([DOH * np.sin(ANG / 2), DOH * np.cos(ANG / 2), 0.0])
    h2 = np.array([-DOH * np.sin(ANG / 2), DOH * np.cos(ANG / 2), 0.0])
    # random rotation (uniform via QR of a Gaussian matrix)
    R, _ = np.linalg.qr(rng.standard_normal((3, 3)))
    if np.linalg.det(R) < 0:
        R[:, 0] = -R[:, 0]
    return o, o + h1 @ R.T, o + h2 @ R.T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--side", type=float, required=True, help="cubic box side (A)")
    ap.add_argument("--target-molality", type=float, default=1.0)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    L = args.side
    n_w0 = int(round(RHO * L ** 3))
    a = args.target_molality * 18.015e-3
    n_pairs = int(round(a * n_w0 / (1.0 + 2.0 * a)))
    n_w = n_w0 - 2 * n_pairs

    # grid of molecule sites
    ns = int(np.ceil(n_w0 ** (1 / 3)))
    spacing = L / ns
    sites = [np.array([i, j, k]) * spacing + spacing / 2
             for i in range(ns) for j in range(ns) for k in range(ns)]
    rng.shuffle(sites)
    sites = sites[:n_w0]

    # first 2*n_pairs sites -> ions (alternating Na/Cl), rest -> water
    ion_sites = sites[:2 * n_pairs]
    wat_sites = sites[2 * n_pairs:]

    atoms = []   # (mol, elem, x,y,z)
    bonds = []   # (a,b)
    angles = []  # (a,o,b)
    mol = 0
    for j, s in enumerate(ion_sites):
        mol += 1
        atoms.append((mol, "Na" if j % 2 == 0 else "Cl", *(s % L)))
    for s in wat_sites:
        mol += 1
        o, h1, h2 = water(s, rng)
        io = len(atoms) + 1
        atoms.append((mol, "O", *(o % L)))
        atoms.append((mol, "H", *(h1 % L)))
        atoms.append((mol, "H", *(h2 % L)))
        bonds.append((io, io + 1)); bonds.append((io, io + 2))
        angles.append((io + 1, io, io + 2))

    with open(args.out, "w") as f:
        f.write(f"Madrid-2019 NaCl(aq): {n_w} H2O + {n_pairs} NaCl, "
                f"{args.target_molality:.2f} mol/kg, {L:.4f} A cubic\n\n")
        f.write(f"{len(atoms)} atoms\n{len(bonds)} bonds\n{len(angles)} angles\n\n")
        f.write("4 atom types\n1 bond types\n1 angle types\n\n")
        f.write(f"0.0 {L:.6f} xlo xhi\n0.0 {L:.6f} ylo yhi\n0.0 {L:.6f} zlo zhi\n\n")
        f.write("Masses\n\n")
        for el, t in sorted(TYPE.items(), key=lambda x: x[1]):
            f.write(f"{t} {MASS[el]:.5f}  # {el}\n")
        f.write("\nAtoms  # full\n\n")
        for i, (m, el, x, y, z) in enumerate(atoms, 1):
            f.write(f"{i} {m} {TYPE[el]} {Q[el]:+.4f} "
                    f"{x:.6f} {y:.6f} {z:.6f}\n")
        f.write("\nBonds\n\n")
        for i, (a1, a2) in enumerate(bonds, 1):
            f.write(f"{i} 1 {a1} {a2}\n")
        f.write("\nAngles\n\n")
        for i, (a1, ao, a2) in enumerate(angles, 1):
            f.write(f"{i} 1 {a1} {ao} {a2}\n")

    print(f"wrote {args.out}: {len(atoms)} atoms ({n_pairs} Na, {n_pairs} Cl, "
          f"{n_w} H2O) in {L:.3f} A cubic -> "
          f"{n_pairs / (n_w * 18.015e-3):.3f} mol/kg")


if __name__ == "__main__":
    main()
