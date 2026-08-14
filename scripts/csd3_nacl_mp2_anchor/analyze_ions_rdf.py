#!/usr/bin/env python3
"""Ion self-diffusion + Na-Cl RDF/PMF from the NaCl(aq) MP2 C-NNP anchor runs.

Companion to analyze_diffusion.py (water O). One streaming pass per segment: awk
prefilters Na/Cl rows of -pos-1.xyz; per-species FFT multi-origin MSD (D over the
same 10-25 ps window); pooled Na-Cl histogram -> g(r) -> w(r) = -kT ln g (tail-ref).
Writes results/nacl_mp2_anchor/{ion_diffusion_summary, msd_ions_<cell>,
rdf_nacl_<cell>}.csv. Usage: analyze_ions_rdf.py [cube2] [cube3] (default: both)
"""
import io
import os
import subprocess
import sys

import numpy as np

# rsync'd copy of /rds/user/crm98/hpc-work/nacl_mp2_anchor (CSD3 original)
RUNS = "/data/phy-cerberus6/crm98/nacl_diffusion/csd3_mp2/runs/MP2/cubic_1M/production"
OUT = os.path.expanduser("~/cp2k-benchmarks/results/nacl_mp2_anchor")
os.makedirs(OUT, exist_ok=True)

DT_FS = 0.5          # MD timestep
TRAJ_EVERY = 10      # &TRAJECTORY EACH MD 10  -> 5 fs between frames
STRIDE = 2           # analyse every 2nd frame -> 10 fs resolution (MSD)
FIT_LO_PS, FIT_HI_PS = 10.0, 25.0   # match analyze_diffusion.py
KT_KCAL = 0.0019872041 * 298.15     # kcal/mol at 298.15 K (PMF reference)
RDF_DR = 0.05        # A
RDF_STRIDE = 1       # every stored frame (5 fs) for the histogram

CELLS = {
    #        L [A]   frames used (incl. frame 0; same truncation as water)
    "cube2": (24.84, 20001),
    "cube3": (37.26, 32001),
}
SEGMENTS = (1, 2, 3, 4, 5)


CACHE = "/data/phy-cerberus6/crm98/nacl_diffusion/csd3_mp2/analysis_cache"
os.makedirs(CACHE, exist_ok=True)


def read_ions(xyz, n_frames_max):
    """Stream Na/Cl rows of the first n_frames_max frames -> (labels, pos (F,nIon,3)).
    Assumes fixed atom order (CP2K prints topology order); cached as npz so
    re-analysis never re-streams the multi-GB trajectories."""
    cache = os.path.join(
        CACHE, os.path.basename(xyz).replace(".xyz", f"_ions_{n_frames_max}.npz"))
    if os.path.exists(cache):
        z = np.load(cache)
        return z["labels"], z["pos"]
    with open(xyz) as f:
        n_atoms = int(f.readline().split()[0])
    max_lines = n_frames_max * (n_atoms + 2)
    awk = subprocess.Popen(
        ["awk", f'NR>{max_lines}{{exit}} $1=="Na"||$1=="Cl"{{print}}', xyz],
        stdout=subprocess.PIPE)
    rows = np.loadtxt(io.TextIOWrapper(awk.stdout),
                      dtype={"names": ("el", "x", "y", "z"),
                             "formats": ("U2", "f8", "f8", "f8")})
    if awk.wait() != 0:
        raise RuntimeError(f"awk failed on {xyz}")
    n_total = rows.shape[0]
    n_ion = n_total // n_frames_max
    if n_ion * n_frames_max != n_total:
        raise RuntimeError(f"{xyz}: {n_total} ion rows not divisible by "
                           f"{n_frames_max} frames")
    labels = rows["el"][:n_ion].copy()
    if not (rows["el"].reshape(n_frames_max, n_ion) == labels).all():
        raise RuntimeError(f"{xyz}: ion order changes between frames")
    pos = np.stack([rows["x"], rows["y"], rows["z"]], axis=-1)
    pos = pos.reshape(n_frames_max, n_ion, 3)
    np.savez_compressed(cache, labels=labels, pos=pos)
    return labels, pos


def unwrap(pos, box_l):
    d = np.diff(pos, axis=0)
    d -= box_l * np.round(d / box_l)
    return np.concatenate([pos[:1], pos[:1] + np.cumsum(d, axis=0)])


def autocorr_fft(x):
    t = x.shape[0]
    f = np.fft.fft(x, n=2 * t, axis=0)
    res = np.fft.ifft(f * f.conjugate(), axis=0)[:t].real
    return res / (t - np.arange(t))[:, None]


def msd_fft(pos):
    """MSD(lag) averaged over atoms; pos (T, N, 3) unwrapped. Vectorized form of
    the verified per-atom algorithm in lammps/msd_multiorigin.py."""
    t, n, _ = pos.shape
    d = np.square(pos).sum(axis=2)
    csum = np.cumsum(d, axis=0)
    total = csum[-1]
    tail = total[None, :] - np.vstack([np.zeros((1, n)), csum[:-1]])
    head = csum[::-1]
    s1 = (tail + head) / (t - np.arange(t))[:, None]
    s2 = sum(autocorr_fft(pos[:, :, k]) for k in range(3))
    return (s1 - 2 * s2).mean(axis=1)


def _selftest():
    """Verify msd_fft against the brute-force multi-origin MSD."""
    rng = np.random.default_rng(7)
    pos = np.cumsum(rng.normal(size=(400, 3, 3)), axis=0)
    fast = msd_fft(pos)
    brute = np.array([
        np.mean(np.square(pos[m:] - pos[:len(pos) - m]).sum(axis=2))
        if m else 0.0 for m in range(len(pos))])
    if not np.allclose(fast[:200], brute[:200], atol=1e-8):
        raise AssertionError("msd_fft disagrees with brute force")


def rdf_accumulate(pos_na, pos_cl, box_l, hist, edges):
    """Add min-image Na-Cl pair distances to hist (frames strided)."""
    for f in range(0, pos_na.shape[0], RDF_STRIDE):
        d = pos_na[f][:, None, :] - pos_cl[f][None, :, :]
        d -= box_l * np.round(d / box_l)
        r = np.sqrt(np.square(d).sum(axis=2)).ravel()
        hist += np.histogram(r, bins=edges)[0]
    return pos_na.shape[0] // RDF_STRIDE + (pos_na.shape[0] % RDF_STRIDE > 0)


def main():
    _selftest()
    cells = sys.argv[1:] or list(CELLS)
    dt_ps = DT_FS * TRAJ_EVERY * STRIDE / 1000.0
    summary = []
    for cell in cells:
        box_l, n_frames = CELLS[cell]
        edges = np.arange(0.0, box_l / 2 + RDF_DR, RDF_DR)
        hist = np.zeros(len(edges) - 1)
        n_rdf_frames = 0
        msd_cols, msd_hdr = [], []
        for seg in SEGMENTS:
            xyz = os.path.join(RUNS, cell, f"seg{seg}",
                               f"NaCl_MP2_{cell}_seg{seg}-pos-1.xyz")
            labels, pos = read_ions(xyz, n_frames)
            na, cl = pos[:, labels == "Na"], pos[:, labels == "Cl"]
            n_rdf_frames += rdf_accumulate(na, cl, box_l, hist, edges)
            for species, p in (("Na", na), ("Cl", cl)):
                uw = unwrap(p[::STRIDE], box_l)
                msd = msd_fft(uw)
                lags = np.arange(len(msd)) * dt_ps
                fit = (lags >= FIT_LO_PS) & (lags <= FIT_HI_PS)
                slope = np.polyfit(lags[fit], msd[fit], 1)[0]
                d_si = slope / 6.0 * 1.0e-8 * 1.0e9        # 1e-9 m2/s
                summary.append((cell, box_l, seg, species, p.shape[1], d_si))
                msd_cols.append(msd)
                msd_hdr.append(f"msd_{species}_seg{seg}_A2")
            print(f"{cell} seg{seg}: " + "  ".join(
                f"D_{s[3]}={s[5]:.4f}e-9" for s in summary[-2:]), flush=True)

        lags = np.arange(len(msd_cols[0])) * dt_ps
        np.savetxt(os.path.join(OUT, f"msd_ions_{cell}.csv"),
                   np.column_stack([lags] + msd_cols),
                   delimiter=",", comments="",
                   header="lag_ps," + ",".join(msd_hdr))

        # g(r): ideal Na-Cl pair count per frame per shell = nNa*nCl/V * 4 pi r^2 dr
        n_na = int((labels == "Na").sum())
        n_cl = int((labels == "Cl").sum())
        r_mid = 0.5 * (edges[:-1] + edges[1:])
        shell = 4.0 / 3.0 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
        ideal = n_na * n_cl / box_l ** 3 * shell
        g = hist / (n_rdf_frames * ideal)
        w = np.full_like(g, np.nan)
        m = g > 0.02
        w[m] = -KT_KCAL * np.log(g[m])
        tail = (r_mid > 9.0) & (r_mid < r_mid.max()) & np.isfinite(w)
        w -= np.nanmean(w[tail])
        np.savetxt(os.path.join(OUT, f"rdf_nacl_{cell}.csv"),
                   np.column_stack([r_mid, g, w]),
                   delimiter=",", comments="",
                   header="r_A,g_NaCl,w_kcalmol")
        print(f"{cell}: RDF pooled over {n_rdf_frames} frames "
              f"({n_na}x{n_cl} pairs)", flush=True)

    with open(os.path.join(OUT, "ion_diffusion_summary.csv"), "w") as f:
        f.write("cell,L_A,segment,species,n_ions,D_1e-9_m2_s\n")
        for row in summary:
            f.write(f"{row[0]},{row[1]},{row[2]},{row[3]},{row[4]},"
                    f"{row[5]:.6f}\n")
    print("done")


if __name__ == "__main__":
    main()
