#!/usr/bin/env python3
"""Campaign cost (GBP) and energy (kWh) ledger for the RPA transport campaign.

Reads sacct, not timings.csv: log_timing.sh appends at job END, so wall-killed
jobs never log a row and timings.csv undercounts the CPU campaign by ~10%.
kWh is a TDP-based ESTIMATE: CSD3 has AcctGatherEnergyType=(null) and
ConsumedEnergyRaw is 0 for every cluster job. GBP rates are CLI inputs (rate
card is behind Raven; record --accessed date for the citation).
Usage: ./energy_cost_ledger.py [--gbp-core-hour X --gbp-gpu-hour Y --accessed DATE] [--pue P]
"""
import argparse
import re
import subprocess
import sys

# hardware pinned from scontrol/lscpu, 2026-08-10
W_ICELAKE_NODE = 700.0    # 2 x Xeon 8368Q (270 W TDP each) + DRAM/board; mid of 650-750 W band
W_AMPERE_PER_GPU = 550.0  # A100-SXM4-80GB 400 W + quarter-node host share

CAMPAIGN_STEPS = 2_400_000        # 3 molalities x 5 segments x 320,000
S_PER_STEP = {"optimised CPU (PR #5295)": 0.2554,
              "master 2026.1 (757bb76a80)": 4.0473}

GROUPS = [("CPU campaign",    r"^rpa_(equil|prod|ext|vacf)$"),
          ("GPU campaign",    r"^rpa_gpu_"),
          ("DFT cost ladder", r"^(dft_ladder|dft_probe|equil_cube1)$"),
          ("analysis",        r"^rpa_ana160$")]


def elapsed_s(t):
    d, _, rest = t.partition("-")
    if rest:
        hh, mm, ss = rest.split(":")
        return int(d) * 86400 + int(hh) * 3600 + int(mm) * 60 + int(ss)
    hh, mm, ss = t.split(":")
    return int(hh) * 3600 + int(mm) * 60 + int(ss)


def collect(since):
    out = subprocess.run(
        ["sacct", "-u", subprocess.run(["whoami"], capture_output=True, text=True).stdout.strip(),
         "-S", since, "-X", "--noheader", "--parsable2",
         "--format=JobName,Partition,State,Elapsed,AllocNodes,CPUTimeRAW,AllocTRES"],
        capture_output=True, text=True).stdout

    acc = {g: dict(jobs=0, node_s=0.0, core_s=0.0, gpu_s=0.0) for g, _ in GROUPS}
    for line in out.splitlines():
        f = line.split("|")
        if len(f) < 7:
            continue
        name, _part, _state, el, nodes, cput, tres = f[:7]
        try:
            el = elapsed_s(el)
        except ValueError:
            continue
        # sub-10 s rows = insurance twins that found nothing to do
        if el <= 10:
            continue
        grp = next((g for g, pat in GROUPS if re.search(pat, name)), None)
        if grp is None:
            continue
        m = re.search(r"gres/gpu=(\d+)", tres)
        a = acc[grp]
        a["jobs"] += 1
        a["node_s"] += el * int(nodes or 0)
        a["core_s"] += int(cput or 0)
        a["gpu_s"] += el * (int(m.group(1)) if m else 0)
    return acc


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--since", default="2026-07-29")
    p.add_argument("--gbp-core-hour", type=float, help="SL3/SL2 GBP per CPU core-hour (icelake)")
    p.add_argument("--gbp-gpu-hour", type=float, help="GBP per GPU-hour (ampere)")
    p.add_argument("--pue", type=float, default=1.0,
                   help="data-centre PUE multiplier; 1.0 = IT load only")
    p.add_argument("--accessed", help="rate-card access date, for the citation")
    a = p.parse_args()

    acc = collect(a.since)
    money = a.gbp_core_hour is not None and a.gbp_gpu_hour is not None

    print(f"=== measured consumption (sacct, since {a.since}) ===\n")
    hdr = f"{'group':<18}{'jobs':>6}{'node-h':>10}{'core-h':>12}{'GPU-h':>9}{'kWh*':>10}"
    if money:
        hdr += f"{'GBP':>10}"
    print(hdr)
    tot = dict(node=0.0, core=0.0, gpu=0.0, kwh=0.0, gbp=0.0)
    for g, _ in GROUPS:
        d = acc[g]
        nh, ch, gh = d["node_s"] / 3600, d["core_s"] / 3600, d["gpu_s"] / 3600
        # ampere jobs hold one GPU (charged per GPU); icelake jobs take the whole node
        kwh = ((gh * W_AMPERE_PER_GPU if gh else nh * W_ICELAKE_NODE) / 1000.0) * a.pue
        row = f"{g:<18}{d['jobs']:>6}{nh:>10.1f}{ch:>12.1f}{gh:>9.1f}{kwh:>10.1f}"
        if money:
            gbp = gh * a.gbp_gpu_hour if gh else ch * a.gbp_core_hour
            row += f"{gbp:>10.2f}"
            tot["gbp"] += gbp
        print(row)
        tot["node"] += nh; tot["core"] += ch; tot["gpu"] += gh; tot["kwh"] += kwh
    row = f"{'TOTAL':<18}{'':>6}{tot['node']:>10.1f}{tot['core']:>12.1f}{tot['gpu']:>9.1f}{tot['kwh']:>10.1f}"
    if money:
        row += f"{tot['gbp']:>10.2f}"
    print(row)

    print(f"\n* kWh is a TDP-BASED ESTIMATE, not a measurement: CSD3 has "
          f"AcctGatherEnergyType=(null),\n  so ConsumedEnergyRaw is 0 for every job on the "
          f"cluster. Assumes {W_ICELAKE_NODE:.0f} W per\n  icelake node "
          f"(2 x Xeon Platinum 8368Q, 270 W TDP each, plus DRAM/board) and "
          f"{W_AMPERE_PER_GPU:.0f} W\n  per A100 including host share"
          + (f", scaled by PUE {a.pue}." if a.pue != 1.0 else "; PUE NOT applied (IT load only)."))

    print("\n=== idealised production-only cost (like-for-like with the DFT ladder) ===")
    print(f"  basis: {CAMPAIGN_STEPS:,} MD steps, 1 node / 76 ranks, no equilibration or re-runs")
    for lab, sps in S_PER_STEP.items():
        nh = CAMPAIGN_STEPS * sps / 3600
        line = f"  {lab:<30}{nh:>10.1f} node-h{nh * 76:>12.1f} core-h"
        if money:
            line += f"{nh * 76 * a.gbp_core_hour:>10.2f} GBP"
        line += f"{nh * W_ICELAKE_NODE / 1000 * a.pue:>10.1f} kWh*"
        print(line)
    r = S_PER_STEP["master 2026.1 (757bb76a80)"] / S_PER_STEP["optimised CPU (PR #5295)"]
    print(f"  measured engine speedup: {r:.2f}x")

    print("\n  NOTE: 'idealised production-only' and 'measured consumption' are "
          "different quantities\n  and must not be mixed in one sentence. The first is what the "
          "science needed and is the\n  fair basis for the DFT comparison (that figure is also "
          "idealised). The second is what\n  the machine actually spent, including equilibration, "
          "insurance arrays and re-runs,\n  and is the honest basis for a GBP or kWh claim.")

    if a.accessed:
        print(f"\n  Cite rates as: University of Cambridge Research Computing Services, "
              f"internal\n  service charges page, accessed {a.accessed}.")
    elif not money:
        print("\n  Pass --gbp-core-hour / --gbp-gpu-hour / --accessed once you have the "
              "rate card\n  from hpc.cam.ac.uk (Raven login) to add the GBP column.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
