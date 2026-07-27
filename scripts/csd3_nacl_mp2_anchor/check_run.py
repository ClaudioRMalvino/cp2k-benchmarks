#!/usr/bin/env python3
"""Physics sanity check for NaCl(aq) MD run directories.

Reads the CP2K .ener file (+ committee energies + trajectory) of each run
directory and reports the three canaries that catch broken forces or an
exploding simulation:

  1. temperature   - mean/max over the last half of the run (NVT should sit
                     near 300 K; NVE should stay near where equilibration
                     left it; > TEMP_HARD K anywhere = exploded)
  2. energy drift  - linear drift of the conserved quantity per atom per ps
                     (THE standard test that integration + forces are
                     consistent; NVE production must pass this)
  3. committee     - mean disagreement (sigma) of the 8 C-NNP models per
                     atom; large values mean the trajectory left the
                     models' training distribution (results untrustworthy
                     even if nothing "crashed")

Usage:
  check_run.py <rundir> [<rundir> ...]
  check_run.py /data/.../np1/production/cell*/seg*     # glob from shell

Exit code 0 = all PASS/WARN, 1 = any FAIL.
"""
import glob
import os
import sys

import numpy as np

# thresholds (per atom, atomic units)
DRIFT_PASS = 2e-6    # a.u./atom/ps
DRIFT_WARN = 2e-5
SIGMA_PASS = 1.5e-4  # a.u./atom  (~4 meV/atom)
SIGMA_WARN = 5e-4
TEMP_HARD = 500.0    # K, anywhere in the run -> exploded


def first(pattern):
    m = sorted(glob.glob(pattern))
    return m[0] if m else None


def natoms_of(rundir):
    xyz = first(os.path.join(rundir, "*-pos-1.xyz"))
    if xyz:
        with open(xyz) as f:
            try:
                return int(f.readline().strip())
            except ValueError:
                pass
    return None


def check(rundir):
    ener = first(os.path.join(rundir, "*-1.ener"))
    if ener is None:
        return "FAIL", f"{rundir}: no .ener file found"
    data = np.loadtxt(ener, comments="#", ndmin=2)
    if len(data) < 10:
        return "WARN", f"{rundir}: only {len(data)} energy records"
    t_ps = data[:, 1] / 1000.0
    temp = data[:, 3]
    cons = data[:, 5]
    if not np.all(np.isfinite(cons)):
        return "FAIL", f"{rundir}: NaN/Inf in conserved quantity - exploded"

    n = natoms_of(rundir)
    half = len(data) // 2
    msgs, level = [], 0

    tmax = temp.max()
    msgs.append(f"T={temp[half:].mean():6.1f}+/-{temp[half:].std():4.1f}K "
                f"(max {tmax:.0f})")
    if tmax > TEMP_HARD:
        level = max(level, 2)
        msgs[-1] += " EXPLODED?"

    # conserved-quantity drift, last 80% (skip initial transient)
    lo = len(data) // 5
    slope = np.polyfit(t_ps[lo:], cons[lo:], 1)[0]  # a.u./ps
    if n:
        d = abs(slope) / n
        msgs.append(f"drift={d:.2e} au/atom/ps")
        level = max(level, 0 if d < DRIFT_PASS else 1 if d < DRIFT_WARN else 2)
    else:
        msgs.append(f"drift={slope:.2e} au/ps (atom count unknown)")

    comm = first(os.path.join(rundir, "*nnp-energies-1.data"))
    if comm and n:
        cdat = np.loadtxt(comm, comments="#", ndmin=2)
        sig = cdat[len(cdat) // 2:, 1].mean() / n
        msgs.append(f"committee sigma={sig:.2e} au/atom")
        level = max(level, 0 if sig < SIGMA_PASS else 1 if sig < SIGMA_WARN else 2)

    verdict = ["PASS", "WARN", "FAIL"][level]
    return verdict, f"{rundir}: {'  '.join(msgs)}"


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    worst = 0
    for rd in sys.argv[1:]:
        verdict, msg = check(rd.rstrip("/"))
        print(f"[{verdict}] {msg}")
        worst = max(worst, ["PASS", "WARN", "FAIL"].index(verdict))
    return 1 if worst == 2 else 0


if __name__ == "__main__":
    sys.exit(main())
