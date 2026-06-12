#!/usr/bin/env python3
"""Fit the empirical size-scaling power law  t(N) = a * N^b  per branch.

This is the *single committed source of truth* for the size-scaling
exponents quoted in the thesis text, the presentation, and the
`power_law_fits.pdf` table.  It loads the same size-scaling CSVs that
`thesis_figures.py` plots, fits a power law in log-log space over one
declared N-range (the SAME range for every branch), and overwrites
`scaling_csd3/power_law_fits.csv`.

Why a fixed range for every branch:
    A single power law does not describe the whole sweep, because at
    small N the 76-core node is underutilised and the step time is
    dominated by fixed overhead (flat region).  The exponent is only
    meaningful in the larger-N regime where the cost actually scales.
    Fitting *every* branch over the same N >= MIN_MOLECULES window keeps
    the comparison honest -- the master and optimised exponents are then
    directly comparable instead of being read off different ranges.

    The resulting exponents are EMPIRICAL over the tested range, not
    asymptotic complexities; the optimised branch is sub-O(N) because
    communication / fixed overhead grows more slowly than N across the
    boxes measured here, not forever.

Run:
    source ~/.fortran_env/bin/activate
    python3 plots/compute_power_law_fits.py
"""
import glob
import os

import numpy as np
import pandas as pd

# --- match thesis_figures.py data loading -----------------------------------
RESULTS_ROOT = "/home/crm98/rds/hpc-work/cp2k-benchmarks/results"
if not os.path.isdir(RESULTS_ROOT):
    RESULTS_ROOT = "/home/crm98/cp2k-benchmarks/results"

OUT_CSV = "/home/crm98/cp2k-benchmarks/plots/scaling_csd3/power_law_fits.csv"

# Atoms per water molecule (N in the fit is the atom count).
ATOMS_PER_MOLECULE = 3

# Fit window: only boxes at or above this many H2O molecules are used, so
# the flat overhead-dominated small-N points do not bias the slope.  The
# SAME window is applied to every branch.  512 mol = 1536 atoms = >=20
# atoms/core on the 76-core node, i.e. the compute-dominated regime.
MIN_MOLECULES = 512

SIZE_COLS = ["n_molecules", "n_reps",
             "tps_mean", "tps_std", "tps_min",
             "wt_mean", "wt_std", "wt_min"]

# branch label -> glob (sub-directory + filename stem), in the CSV order.
BRANCHES = {
    "upstream master":
        "cp2k_master/NNP/NNP_size_scaling_upstream-master_*",
    "native-spline":
        "cp2k_feature_native_spline/NNP/NNP_size_scaling_feature-nnp-native-spline_*",
    "native-spline-omp":
        "cp2k_feature_native_spline_omp/NNP/NNP_size_scaling_feature-nnp-native-spline-omp_*",
}


def _latest(pattern):
    """Most recently modified non-raw CSV matching pattern (figure logic)."""
    matches = [m for m in glob.glob(pattern) if not m.endswith("_raw.csv")]
    return max(matches, key=os.path.getmtime) if matches else None


def _load_size(branch_glob):
    path = _latest(f"{RESULTS_ROOT}/{branch_glob}/results_size_scaling_*.csv")
    if path is None:
        raise FileNotFoundError(f"no size-scaling CSV under {branch_glob}")
    df = pd.read_csv(path, comment="#", header=None,
                     names=SIZE_COLS, na_values=["NA"])
    df = df.sort_values("n_molecules").reset_index(drop=True)
    return df, path


def _fit_power_law(n_atoms, t):
    """Weighted-free OLS of log t on log N.  Returns a, a_err, b, b_err.

    b is the power-law exponent (slope in log-log space); a is the
    prefactor t = a * N^b.  Errors are the 1-sigma standard errors from
    the regression covariance (a_err via the delta method on a=exp(c)).
    """
    x = np.log(n_atoms.astype(float))
    y = np.log(t.astype(float))
    # degree-1 fit with covariance; needs n >= 3 for a meaningful cov.
    (b, c), cov = np.polyfit(x, y, 1, cov=True)
    b_err = float(np.sqrt(cov[0, 0]))
    c_err = float(np.sqrt(cov[1, 1]))
    a = float(np.exp(c))
    a_err = a * c_err          # delta method: d(exp c) = exp(c) dc
    return a, a_err, float(b), b_err


def main():
    rows = []
    print(f"RESULTS_ROOT : {RESULTS_ROOT}")
    print(f"Fit window   : N >= {MIN_MOLECULES} H2O "
          f"(>= {MIN_MOLECULES * ATOMS_PER_MOLECULE} atoms), same for all branches")
    print(f"Model        : t(N) = a * N^b,  N = atom count, OLS in log-log\n")
    print(f"{'branch':<20} {'n':>2}  {'exponent b':>16}  {'prefactor a (s)':>20}")
    print("-" * 64)

    for branch, branch_glob in BRANCHES.items():
        df, path = _load_size(branch_glob)
        sel = df[df["n_molecules"] >= MIN_MOLECULES]
        n_atoms = sel["n_molecules"].to_numpy() * ATOMS_PER_MOLECULE
        t = sel["tps_mean"].to_numpy()
        a, a_err, b, b_err = _fit_power_law(n_atoms, t)
        rows.append({
            "branch": branch, "a": a, "a_err": a_err,
            "b": b, "b_err": b_err, "n_points": len(sel),
        })
        print(f"{branch:<20} {len(sel):>2}  "
              f"{b:>7.3f} +/- {b_err:<5.3f}  "
              f"{a:>9.2e} +/- {a_err:<8.2e}  [{os.path.basename(path)}]")

    out = pd.DataFrame(rows, columns=["branch", "a", "a_err",
                                      "b", "b_err", "n_points"])
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    out.to_csv(OUT_CSV, index=False)
    print(f"\nwrote {OUT_CSV}")
    print("Re-render the table with:  python3 plots/csd3_tables.py")


if __name__ == "__main__":
    main()
