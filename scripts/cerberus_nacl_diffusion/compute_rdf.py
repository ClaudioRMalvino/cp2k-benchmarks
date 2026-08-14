#!/usr/bin/env python3
"""Radial distribution functions from CP2K XYZ trajectories.

Computes g(r) for O-O, Na-O, Cl-O (the pairs in Fig. 1 of O'Neill et al. JPCL
2024) to check production trajectories against the published curves. Accepts
several trajectory files (concatenated for statistics); orthorhombic min-image.
Example: compute_rdf.py seg*/*-pos-1.xyz --cell-ang 24.84 24.84 12.42 --out-prefix rdf_cell221 --plot
"""
import argparse

import numpy as np

PAIR_DEFAULT = ["O-O", "Na-O", "Cl-O"]


def frames(path):
    """Yield (symbols, coords) per frame; symbols only on first frame."""
    with open(path) as f:
        syms = None
        while True:
            line = f.readline()
            if not line:
                return
            try:
                nat = int(line.strip())
            except ValueError:
                continue
            f.readline()
            s, c = [], np.empty((nat, 3))
            for i in range(nat):
                p = f.readline().split()
                s.append(p[0])
                c[i] = [float(p[1]), float(p[2]), float(p[3])]
            if syms is None:
                syms = np.array(s)
            yield syms, c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("trajectories", nargs="+")
    ap.add_argument("--cell-ang", type=float, nargs=3, required=True)
    ap.add_argument("--pairs", nargs="+", default=PAIR_DEFAULT,
                    help="element pairs like Na-O")
    ap.add_argument("--rmax", type=float, default=6.0)
    ap.add_argument("--dr", type=float, default=0.02)
    ap.add_argument("--stride", type=int, default=5,
                    help="use every Nth frame (RDF frames are correlated "
                         "anyway; stride 5 = every 25 fs)")
    ap.add_argument("--out-prefix", default="rdf")
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()

    cell = np.array(args.cell_ang)
    rmax = min(args.rmax, cell.min() / 2.0)
    nbins = int(rmax / args.dr)
    edges = np.linspace(0.0, rmax, nbins + 1)
    vol = float(np.prod(cell))

    pairs = [tuple(p.split("-")) for p in args.pairs]
    hists = {p: np.zeros(nbins) for p in pairs}
    counts = {p: None for p in pairs}          # (n_A, n_B) per pair
    nframes = 0

    for traj in args.trajectories:
        for k, (syms, xyz) in enumerate(frames(traj)):
            if k % args.stride:
                continue
            nframes += 1
            idx = {el: np.where(syms == el)[0] for el in
                   set(e for p in pairs for e in p)}
            for (a, b) in pairs:
                ia, ib = idx[a], idx[b]
                if len(ia) == 0 or len(ib) == 0:
                    continue
                counts[(a, b)] = (len(ia), len(ib))
                d = xyz[ia][:, None, :] - xyz[ib][None, :, :]
                d -= cell * np.round(d / cell)
                r = np.linalg.norm(d, axis=2).ravel()
                if a == b:
                    r = r[r > 1e-6]            # drop self-pairs
                hists[(a, b)] += np.histogram(r, bins=edges)[0]

    r_mid = 0.5 * (edges[1:] + edges[:-1])
    shell = (4.0 / 3.0) * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    results = {}
    for p, h in hists.items():
        if counts[p] is None or nframes == 0:
            continue
        na, nb = counts[p]
        npairs = na * (nb - 1) if p[0] == p[1] else na * nb
        g = h * vol / (npairs * shell * nframes)
        # running coordination number of B around A
        rho_b = (nb - (1 if p[0] == p[1] else 0)) / vol
        n_r = np.cumsum(g * shell * rho_b)
        results[p] = (g, n_r)
        out = f"{args.out_prefix}_{p[0]}-{p[1]}.csv"
        np.savetxt(out, np.c_[r_mid, g, n_r], delimiter=",",
                   header=f"r_ang,g_r,coord_n  ({nframes} frames, "
                          f"{na}x{nb} atoms)", comments="# ")
        imax = np.argmax(g)
        print(f"{p[0]}-{p[1]}: first peak g={g[imax]:.2f} at "
              f"r={r_mid[imax]:.2f} A  -> {out}")

    if args.plot and results:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, len(results), figsize=(4 * len(results), 3.6))
        axes = np.atleast_1d(axes)
        for ax, (p, (g, _)) in zip(axes, results.items()):
            ax.plot(r_mid, g, color="#133844")
            ax.axhline(1.0, color="gray", lw=0.6, ls=":")
            ax.set_xlabel(r"r [$\mathrm{\AA}$]")
            ax.set_ylabel(f"g(r) {p[0]}-{p[1]}")
            ax.set_xlim(0, rmax)
        fig.tight_layout()
        fig.savefig(f"{args.out_prefix}.png", dpi=300)
        print(f"wrote {args.out_prefix}.png")


if __name__ == "__main__":
    main()
