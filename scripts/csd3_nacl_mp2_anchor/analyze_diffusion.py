#!/usr/bin/env python3
"""Water-oxygen self-diffusion from the NaCl(aq) MP2 C-NNP anchor runs.

For each cell (cube2 24.84 A, cube3 37.26 A) and each of the 5 NVE production
segments: stream the -pos-1.xyz (O atoms only, via awk prefilter), unwrap
minimum-image, MSD by the FFT multiple-origin algorithm, D from a linear fit
of MSD over a common lag window. Segments are truncated to a common length
per cell so all five are statistically equivalent.

Outputs (data -> results/, per benchmarks layout):
  results/nacl_mp2_anchor/diffusion_summary.csv   one row per (cell, segment)
  results/nacl_mp2_anchor/msd_<cell>.csv          lag_ps + one MSD column/segment

Usage: analyze_diffusion.py [cube2] [cube3]   (default: both)
"""
import io
import os
import subprocess
import sys

import numpy as np

RUNS = "/rds/user/crm98/hpc-work/nacl_mp2_anchor/runs/MP2/cubic_1M/production"
OUT = os.path.expanduser("~/cp2k-benchmarks/results/nacl_mp2_anchor")
os.makedirs(OUT, exist_ok=True)

DT_FS = 0.5          # MD timestep
TRAJ_EVERY = 10      # &TRAJECTORY EACH MD 10  -> 5 fs between frames
STRIDE = 2           # analyse every 2nd frame -> 10 fs resolution
FIT_LO_PS, FIT_HI_PS = 10.0, 25.0   # common diffusive-regime fit window

CELLS = {
    #        L [A]   frames used (incl. frame 0; common truncation per cell)
    "cube2": (24.84, 20001),   # 100 ps, segments completed exactly
    "cube3": (37.26, 32001),   # 160 ps, min over segments (overran 100 ps target)
}
SEGMENTS = (1, 2, 3, 4, 5)


def read_oxygens(xyz, n_frames_max):
    """Stream O-atom coordinates of the first n_frames_max frames -> (F, nO, 3)."""
    with open(xyz) as f:
        n_atoms = int(f.readline().split()[0])
    max_lines = n_frames_max * (n_atoms + 2)
    awk = subprocess.Popen(
        ["awk", f'NR>{max_lines}{{exit}} $1=="O"{{print $2, $3, $4}}', xyz],
        stdout=subprocess.PIPE)
    raw = np.loadtxt(io.TextIOWrapper(awk.stdout), dtype=np.float64)
    if awk.wait() != 0:
        raise RuntimeError(f"awk failed on {xyz}")
    n_o_total = raw.shape[0]
    # nO per frame from frame count (all frames complete by construction)
    n_o = n_o_total // n_frames_max
    if n_o * n_frames_max != n_o_total:
        raise RuntimeError(f"{xyz}: {n_o_total} O rows not divisible by "
                           f"{n_frames_max} frames")
    return raw.reshape(n_frames_max, n_o, 3)


def unwrap(pos, box_l):
    d = np.diff(pos, axis=0)
    d -= box_l * np.round(d / box_l)
    return np.concatenate([pos[:1], pos[:1] + np.cumsum(d, axis=0)])


def autocorr_fft(x):
    """Column-wise autocorrelation sum S2 for the FFT MSD algorithm. x: (T, M)."""
    t = x.shape[0]
    f = np.fft.fft(x, n=2 * t, axis=0)
    res = np.fft.ifft(f * f.conjugate(), axis=0)[:t].real
    return res / (t - np.arange(t))[:, None]


def msd_fft(pos):
    """MSD(lag) averaged over atoms. pos: (T, N, 3) unwrapped."""
    t, n, _ = pos.shape
    d = np.square(pos).sum(axis=2)                    # (T, N)
    csum = np.cumsum(d, axis=0)
    total = csum[-1]
    tail = total[None, :] - np.vstack([np.zeros((1, n)), csum[:-1]])
    s1 = (tail + csum[::-1]) / (t - np.arange(t))[:, None]
    s2 = np.zeros((t, n))
    for lo in range(0, n, 128):                       # chunk atoms for FFT memory
        hi = min(lo + 128, n)
        blk = pos[:, lo:hi, :].reshape(t, -1)
        s2[:, lo:hi] = autocorr_fft(blk).reshape(t, hi - lo, 3).sum(axis=2)
    return (s1 - 2 * s2).mean(axis=1)                 # (T,) A^2


def mean_temperature(rundir, proj):
    ener = os.path.join(rundir, f"{proj}-1.ener")
    data = np.loadtxt(ener, comments="#", ndmin=2)
    return data[:, 3].mean()


def analyse_cell(cell):
    box_l, n_frames = CELLS[cell]
    dt_frame_ps = DT_FS * TRAJ_EVERY * STRIDE / 1000.0
    rows, msd_cols = [], {}
    for seg in SEGMENTS:
        rundir = os.path.join(RUNS, cell, f"seg{seg}")
        proj = f"NaCl_MP2_{cell}_seg{seg}"
        pos = read_oxygens(os.path.join(rundir, f"{proj}-pos-1.xyz"), n_frames)
        pos = unwrap(pos, box_l)[::STRIDE]
        msd = msd_fft(pos)
        lag = np.arange(len(msd)) * dt_frame_ps
        sel = (lag >= FIT_LO_PS) & (lag <= FIT_HI_PS)
        slope, icpt = np.polyfit(lag[sel], msd[sel], 1)
        d_si = slope / 6.0 * 1e-8                     # A^2/ps -> m^2/s
        temp = mean_temperature(rundir, proj)
        rows.append((cell, box_l, seg, d_si * 1e9, temp, pos.shape[0],
                     pos.shape[1], FIT_LO_PS, FIT_HI_PS))
        msd_cols[f"seg{seg}"] = msd
        print(f"{cell} seg{seg}: D = {d_si*1e9:.4f} e-9 m^2/s   "
              f"T = {temp:.1f} K   ({pos.shape[0]} frames x {pos.shape[1]} O)",
              flush=True)
    # MSD curves, downsampled x10, up to lag = half trajectory
    half = len(msd_cols["seg1"]) // 2
    lag = (np.arange(len(msd_cols["seg1"])) * dt_frame_ps)[:half:10]
    table = np.column_stack([lag] + [msd_cols[f"seg{s}"][:half:10]
                                     for s in SEGMENTS])
    np.savetxt(os.path.join(OUT, f"msd_{cell}.csv"), table, delimiter=",",
               header="lag_ps," + ",".join(f"msd_seg{s}_A2" for s in SEGMENTS),
               comments="")
    return rows


def main():
    cells = sys.argv[1:] or list(CELLS)
    all_rows = []
    for cell in cells:
        all_rows += analyse_cell(cell)
    path = os.path.join(OUT, "diffusion_summary.csv")
    with open(path, "w") as f:
        f.write("cell,L_A,segment,D_1e-9_m2_s,T_K,frames_used,n_O,"
                "fit_lo_ps,fit_hi_ps\n")
        for r in all_rows:
            f.write(f"{r[0]},{r[1]},{r[2]},{r[3]:.6f},{r[4]:.2f},"
                    f"{r[5]},{r[6]},{r[7]},{r[8]}\n")
    print(f"wrote {path}")
    for cell in cells:
        ds = [r[3] for r in all_rows if r[0] == cell]
        print(f"{cell}: D = {np.mean(ds):.4f} +/- "
              f"{np.std(ds, ddof=1)/np.sqrt(len(ds)):.4f} e-9 m^2/s "
              f"(mean +/- SEM, {len(ds)} segments)")


if __name__ == "__main__":
    main()
