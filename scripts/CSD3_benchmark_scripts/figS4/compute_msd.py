#!/usr/bin/env python3
"""Read an NVE positions trajectory (XYZ) and compute the self-diffusion
coefficient of water (oxygen atoms only) by a linear fit of the mean
squared displacement.

Usage:
    compute_msd.py <traj.xyz> <dt_fs> <box_Lx_A> <box_Ly_A> <box_Lz_A>

dt_fs is the time between frames in femtoseconds (frame interval × MD timestep).
box dimensions used only for unwrapping (NVE positions printed by CP2K are
typically wrapped into the cell).
"""
import sys
import numpy as np

if len(sys.argv) != 6:
    print("usage: compute_msd.py traj.xyz dt_fs Lx_A Ly_A Lz_A", file=sys.stderr)
    sys.exit(2)

path  = sys.argv[1]
dt_fs = float(sys.argv[2])
box   = np.array([float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])])

# --- read trajectory ---------------------------------------------------------
frames = []
syms   = None
with open(path) as f:
    while True:
        hdr = f.readline()
        if not hdr: break
        n = int(hdr.strip())
        _ = f.readline()        # comment line
        pos  = np.empty((n, 3))
        elem = []
        for i in range(n):
            p = f.readline().split()
            elem.append(p[0])
            pos[i] = [float(p[1]), float(p[2]), float(p[3])]
        if syms is None: syms = elem
        frames.append(pos)
coords = np.array(frames)                # (n_frames, n_atoms, 3)  Angstrom
n_frames, n_atoms, _ = coords.shape

# --- oxygen atoms only -------------------------------------------------------
o_idx = [i for i, s in enumerate(syms) if s.upper() == 'O']
pos_o = coords[:, o_idx, :]              # (n_frames, n_O, 3)
n_O   = pos_o.shape[1]

# --- unwrap (minimum-image, walking forward) ---------------------------------
unwrap = pos_o.copy()
for t in range(1, n_frames):
    d = unwrap[t] - unwrap[t-1]
    unwrap[t] = unwrap[t-1] + d - box * np.round(d / box)

# --- MSD by FFT trick (Allen & Tildesley appendix; O(N log N) per atom) ----
# msd(tau) = <|r(t+tau) - r(t)|^2>_t,atoms
def msd_fft(x):
    """x : (n_frames, n_atoms, 3) → MSD : (n_frames,)"""
    N = x.shape[0]
    # square term: 2*<r_t·r_t> averaged
    D = np.square(x).sum(axis=2)         # (N, n_atoms)
    D2 = np.append(D, np.zeros(D.shape[1])[None, :], axis=0)
    Q = 2 * D.sum(axis=0)                 # (n_atoms,)
    S1 = np.zeros((N, x.shape[1]))
    for t in range(N):
        Q  = Q - D2[t-1] - D2[N-t]
        S1[t] = Q / (N - t)
    # autocorrelation S2[t] = <r_{i+t}·r_i>_i via FFT
    S2 = np.zeros((N, x.shape[1]))
    for dim in range(3):
        f = np.fft.fft(x[:, :, dim], n=2*N, axis=0)
        psd = f * f.conjugate()
        s = np.fft.ifft(psd, axis=0).real
        s = s[:N]
        norm = (N - np.arange(N))[:, None]
        S2 += s / norm
    return (S1 - 2 * S2).mean(axis=1)     # average over atoms

msd = msd_fft(unwrap)                    # (n_frames,)  Å²
t_ps = np.arange(n_frames) * dt_fs / 1000.0

# --- linear fit over the middle of the trajectory (10%..50%) ----------------
i0, i1 = int(n_frames * 0.10), int(n_frames * 0.50)
slope, intercept = np.polyfit(t_ps[i0:i1], msd[i0:i1], 1)

D_A2_ps    = slope / 6.0                 # Å²/ps
D_cm2_s_e5 = D_A2_ps * 1.0e-1            # 1 Å²/ps = 1e-5 cm²/s   (1e-16 cm²/1e-12 s = 1e-4 cm²/s; the ×0.1 converts to 10⁻⁵ cm²/s)

# Yeh-Hummer (need viscosity, which we don't have — print the structure only):
#   D_inf = D_PBC + k_B T * xi / (6 pi eta L) , xi = 2.837297
L_min  = box.min()
yh_pre = 8.314 * 300.0 * 2.837297 / (6 * np.pi * L_min * 1e-10) / 1000  # using L in m
# practically: D_inf [10⁻⁵ cm²/s] ≈ D_PBC + 0.0247/(eta_mPa_s · L_A/10)
yh_const = 2.837297 * 1.38064852e-23 * 300.0 / (6 * np.pi)
yh_corr_per_invvisc_per_invL = yh_const / (L_min * 1.0e-10) / 1.0e-3      # m²/s per (mPa·s)⁻¹
# Convert to 10⁻⁵ cm²/s:  1 m²/s = 1e4 cm²/s = 1e9 ·10⁻⁵ cm²/s
yh_corr_per_invvisc_units = yh_corr_per_invvisc_per_invL * 1.0e9

print(f"trajectory: {path}")
print(f"  frames                : {n_frames}    (dt = {dt_fs} fs → total {t_ps[-1]:.2f} ps)")
print(f"  oxygen atoms used     : {n_O}")
print(f"  cell                  : {box[0]:.3f} x {box[1]:.3f} x {box[2]:.3f} Å  (Lmin = {L_min:.3f})")
print(f"  fit range             : {t_ps[i0]:.2f} .. {t_ps[i1]:.2f} ps  (indices {i0}..{i1})")
print(f"  MSD slope             : {slope:.5f} Å²/ps   (intercept {intercept:+.3f} Å²)")
print(f"")
print(f"  ==> D_PBC  =  {D_A2_ps:.5f} Å²/ps  =  {D_cm2_s_e5:.4f} × 10⁻⁵ cm²/s")
print(f"")
print(f"  Yeh-Hummer correction would add:  D_inf - D_PBC  =  {yh_corr_per_invvisc_units:.4f} / η[mPa·s]   (×10⁻⁵ cm²/s)")
print(f"  Morawietz cites experimental η ≈ 0.85 mPa·s at 300 K")
print(f"  reported NNP D at Morawietz equilibrium (Fig 1D): D_RPBE-vdW ≈ 0.3 × 10⁻⁵ cm²/s")
