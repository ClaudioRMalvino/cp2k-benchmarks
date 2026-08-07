#!/usr/bin/env python3
"""
Fit the measured on-the-fly DFT cost ladder and build the full cost table.

Reads results/dft_cost/dft_ladder_timings.csv (four measured rungs: 192, 384,
768, 1536 atoms of liquid water at revPBE-D3, the level the NNP was trained
to), fits cost per MD step against system size in log-log space, and
extrapolates to the production cells (5064 / 4952 / 4748 atoms).

The output table places that extrapolation alongside the three MEASURED NNP
rungs (master CP2K, optimised CPU, GPU) and the analytic RI-RPA estimate, and
converts each into the cost of the round-1 production campaign.

    source ~/.fortran_env/bin/activate
    python analyze_dft_ladder.py
"""
import csv
import os

import numpy as np

REPO = "/home/crm98/cp2k-benchmarks"
LADDER = f"{REPO}/results/dft_cost/dft_ladder_timings.csv"
OUT = f"{REPO}/results/dft_cost/dft_cost_summary.csv"

# --- campaign definition (round 1 production) ------------------------------
N_SEG, STEPS_PER_SEG, N_CONC = 5, 160_000, 3
CAMPAIGN_STEPS = N_SEG * STEPS_PER_SEG * N_CONC        # 2.4e6 MD steps
FS_PER_STEP = 0.5
CORES_PER_NODE = 76
SL3_ALLOCATION_CORE_H = 84_000
# CSD3 icelake partition, for the "does it even fit" check: 552 nodes x 256 GB
ICELAKE_TOTAL_RAM_BYTES = 552 * 256_120 * 1024**2

# production cells (atoms); 1M is the reference for the headline number
CELLS = {"cubic_1M": 5064, "cubic_2m": 4952, "cubic_4m": 4748}
N_PROD = CELLS["cubic_1M"]

# --- measured NNP rungs (this campaign, icelake 1 node / 76 ranks; GPU =
#     1xA100 + 3 SPME ranks). These are measurements, not estimates. -------
NNP = {
    "NNP, master CP2K (757bb76a80)": (4.0473, "node"),
    "NNP, optimised CPU (PR #5295)": (0.2554, "node"),
    "NNP, GPU (PR #5295 + nnp_gpu)": (0.1467, "gpu"),
}

# --- analytic RI-RPA rung ---------------------------------------------------
# O'Neill's RPA deck: RI-RPA on PBE, cc-TZ + RI_AUX, 20 minimax quadrature
# points, truncated-Coulomb exact exchange, run on 168-atom cells.
#   time   scales O(N^4)  (standard RI-RPA)
#   memory scales O(N^3)  (the (ia|P) three-centre integrals)
# Both are anchored at the 168-atom training configuration.
RPA_REF_ATOMS = 168
RPA_TIME_EXPONENT = 4.0
RPA_MEM_EXPONENT = 3.0
# RI integral count at the 168-atom reference: N_occ * N_virt * N_aux
RPA_REF_INTEGRALS = 224 * 3024 * 11_000        # ~7.4e9 doubles ~ 60 GB
# RI-RPA single point vs a GGA SCF on the same 168-atom cell (order of
# magnitude; 20 quadrature points each requiring the full RI contraction)
RPA_OVER_DFT_AT_REF = 300.0


def fit_ladder(rows):
    """Fit log(s/step) = log(a) + alpha*log(N). Returns (alpha, a, r2)."""
    n = np.array([float(r["natoms"]) for r in rows])
    t = np.array([float(r["s_per_step"]) for r in rows])
    order = np.argsort(n)
    n, t = n[order], t[order]
    alpha, loga = np.polyfit(np.log(n), np.log(t), 1)
    pred = np.exp(loga) * n**alpha
    ss_res = np.sum((np.log(t) - np.log(pred)) ** 2)
    ss_tot = np.sum((np.log(t) - np.log(t).mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return alpha, np.exp(loga), r2, n, t


def node_hours(s_per_step, steps):
    return s_per_step * steps / 3600.0


def main():
    if not os.path.exists(LADDER):
        raise SystemExit(f"no ladder data yet: {LADDER}\n"
                         "run sbatch_dft_ladder.sh rungs 1-4 first")
    rows = [r for r in csv.DictReader(open(LADDER)) if r.get("s_per_step")]
    # keep the newest row per system size
    latest = {}
    for r in rows:
        latest[int(r["natoms"])] = r
    rows = [latest[k] for k in sorted(latest)]
    if len(rows) < 3:
        raise SystemExit(f"need >=3 rungs to fit, have {len(rows)}")

    alpha, a, r2, n, t = fit_ladder(rows)
    print(f"=== measured on-the-fly DFT ladder (revPBE-D3, TZV2P-GTH, 1200 Ry) ===")
    for ni, ti in zip(n, t):
        print(f"  {int(ni):5d} atoms  {ti:10.3f} s/step")
    print(f"\n  fit: s/step = {a:.3e} * N^{alpha:.3f}   (log-log R^2 = {r2:.4f})")

    dft_prod = a * N_PROD**alpha
    print(f"  extrapolated to {N_PROD} atoms: {dft_prod:,.0f} s/step "
          f"= {dft_prod/3600:.2f} h per MD step")

    # RPA, analytic
    scale = N_PROD / RPA_REF_ATOMS
    rpa_prod = dft_prod * RPA_OVER_DFT_AT_REF * scale**RPA_TIME_EXPONENT / scale**alpha
    rpa_bytes = RPA_REF_INTEGRALS * scale**RPA_MEM_EXPONENT * 8
    print(f"\n=== analytic RI-RPA at {N_PROD} atoms ===")
    print(f"  size factor vs {RPA_REF_ATOMS}-atom training cell: {scale:.1f}x")
    print(f"  time   O(N^{RPA_TIME_EXPONENT:.0f}) -> {scale**RPA_TIME_EXPONENT:,.0f}x  "
          f"=> ~{rpa_prod:.3e} s/step")
    print(f"  memory O(N^{RPA_MEM_EXPONENT:.0f}) -> {scale**RPA_MEM_EXPONENT:,.0f}x  "
          f"=> {rpa_bytes/1e15:.2f} PB for the (ia|P) integrals alone")
    print(f"  that is {rpa_bytes/ICELAKE_TOTAL_RAM_BYTES:.1f}x the TOTAL RAM of the "
          f"entire 552-node icelake partition ({ICELAKE_TOTAL_RAM_BYTES/1e12:.0f} TB)")
    print("  -> on-the-fly RPA at production size is memory-infeasible, not merely slow")

    # assemble the table
    table = []
    table.append(("On-the-fly RPA (analytic)", rpa_prod, "node", "infeasible: memory"))
    table.append(("On-the-fly revPBE-D3 (extrapolated)", dft_prod, "node", "measured ladder + fit"))
    for label, (sps, unit) in NNP.items():
        table.append((label, sps, unit, "measured"))

    ref_gpu = NNP["NNP, GPU (PR #5295 + nnp_gpu)"][0]
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rung", "s_per_step_5064atoms", "unit", "campaign_node_or_gpu_h",
                    "campaign_core_h", "ratio_vs_gpu", "basis"])
        for label, sps, unit, basis in table:
            h = node_hours(sps, CAMPAIGN_STEPS)
            core_h = h * CORES_PER_NODE if unit == "node" else ""
            w.writerow([label, f"{sps:.6g}", unit, f"{h:.6g}",
                        f"{core_h:.6g}" if core_h != "" else "",
                        f"{sps/ref_gpu:.6g}", basis])

    print(f"\n=== campaign cost ({CAMPAIGN_STEPS:,} MD steps = "
          f"{CAMPAIGN_STEPS*FS_PER_STEP/1e6:.1f} ns aggregate) ===")
    for label, sps, unit, _ in table:
        h = node_hours(sps, CAMPAIGN_STEPS)
        tag = "GPU-h" if unit == "gpu" else "node-h"
        extra = f" = {h*CORES_PER_NODE:,.0f} core-h" if unit == "node" else ""
        print(f"  {label:38s} {h:14,.1f} {tag}{extra}")

    # inverse framing: what the whole allocation buys on the fly
    steps_afforded = SL3_ALLOCATION_CORE_H / (dft_prod * CORES_PER_NODE / 3600)
    ps_afforded = steps_afforded * FS_PER_STEP / 1000
    print(f"\n  the entire {SL3_ALLOCATION_CORE_H:,} core-h SL3 allocation buys "
          f"{steps_afforded:,.0f} on-the-fly DFT steps = {ps_afforded:.3f} ps")
    print(f"  the campaign actually produced "
          f"{CAMPAIGN_STEPS*FS_PER_STEP/1000:,.0f} ps  "
          f"({CAMPAIGN_STEPS*FS_PER_STEP/1000/ps_afforded:,.0f}x more trajectory)")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
