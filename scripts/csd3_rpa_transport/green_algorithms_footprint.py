#!/usr/bin/env python3
"""Carbon footprint of the RPA NaCl(aq) transport campaign, Green Algorithms
methodology (Lannelongue, Grealey & Inouye, Adv. Sci. 2021, 8:2100707).

Not the GreenAlgorithms4HPC tool: it parses sacct, but the master and 200 ns
rows are PROJECTIONS with no sacct records; report Table 9.4 holds the hours.
(Tool installed at /rds/user/crm98/hpc-work/GreenAlgorithms4HPC for cross-checks.)

THE FORMULA (identical to the tool's):

    E [kWh] = t[h] x (n_cores x TDP_core x usage + n_GPU x TDP_GPU
                      + mem_GB x 0.3725 W/GB) / 1000 x PUE x PSF
    carbon  = E x CI

PUE and CI are pure multipliers; see the sensitivity block in the output.

RESERVED, NOT USED: GPU jobs reserved a quarter ampere node (AllocTRES
billing=32,cpu=32,gres/gpu=1,mem=250G) though only ~3 cores carry SPME; carbon
charges what is reserved. Table 9.4's "3 SPME cores" fits the GBP column, not here.
"""
import argparse
import json
import sys
import urllib.request
from datetime import date, timedelta

# --- Green Algorithms fixed parameters (data/fixed_parameters.yaml) ---------
W_PER_GB          = 0.3725     # W/GB of memory reserved
TREE_MONTH_G      = 917.0      # gCO2e absorbed by a mature tree in a month
CAR_EU_G_PER_KM   = 175.0      # gCO2e/km, average EU passenger car
FLIGHT_PAR_LON_G  = 50_000.0   # gCO2e, one passenger, Paris-London

# --- hardware, from the nodes (scontrol/lscpu/nvidia-smi, 2026-08-11) ------
# icelake: 2 x Xeon 8368Q, 270 W / 38-core socket; AllocTRES billing=76,cpu=76,mem=256120M
ICELAKE = dict(cores=76, w_per_core=270.0 / 38, gpus=0, w_per_gpu=0.0, mem_gb=256120 / 1024)
# ampere: EPYC 7763, 280 W / 64-core socket; A100-SXM4-80GB 400 W (250-300 W is
# the PCIe card, not SXM4); AllocTRES billing=32,cpu=32,gres/gpu=1,mem=250G
AMPERE_QUARTER = dict(cores=32, w_per_core=280.0 / 64, gpus=1, w_per_gpu=400.0, mem_gb=250.0)

# campaign hours from report Table 9.4 (canon after the 2026-08-11 factor-2
# correction: 340 / 5,396 / 227)
SCENARIOS = {
    "as run: 2.4 ns (3 conc x 5 x 160 ps, 4.8M steps)": [
        # label,                      hours, profile,          unit
        ("master (non-optimised)",    5396,  ICELAKE,          "node-h"),
        ("optimised CPU",              340,  ICELAKE,          "node-h"),
        ("optimised GPU",              227,  AMPERE_QUARTER,   "GPU-h"),
    ],
    "projected: 200 ns (400M steps, O'Neill-scale sampling)": [
        ("master (non-optimised)",  450000,  ICELAKE,          "node-h"),
        ("optimised CPU",            28000,  ICELAKE,          "node-h"),
        ("optimised GPU",            19000,  AMPERE_QUARTER,   "GPU-h"),
    ],
}

CAMPAIGN_START, CAMPAIGN_END = date(2026, 5, 10), date(2026, 8, 11)
CI_FALLBACK = 231.12      # Green Algorithms UK value, used only if the API fails


def fetch_uk_ci(start, end, postcode=None, verbose=True):
    """Mean UK grid CI over [start, end] from the National Grid API; returns
    (gCO2e/kWh, description) or (None, reason). Measured daily means for the run
    dates: UK CI has ~halved since the 2021 GA default (would overstate ~2x)."""
    base = "https://api.carbonintensity.org.uk/intensity/stats"
    vals, cur = [], start
    while cur < end:
        chunk_end = min(cur + timedelta(days=29), end)
        url = f"{base}/{cur:%Y-%m-%d}T00:00Z/{chunk_end:%Y-%m-%d}T00:00Z/24"
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                data = json.load(r)["data"]
            vals += [d["intensity"]["average"] for d in data]
        except Exception as exc:                       # offline, API down, etc.
            return None, f"API unavailable ({type(exc).__name__}: {exc})"
        cur = chunk_end
    if not vals:
        return None, "API returned no data"
    mean = sum(vals) / len(vals)
    if verbose:
        print(f"  UK grid CI {start} to {end}: mean {mean:.1f} gCO2e/kWh "
              f"over {len(vals)} daily blocks (min {min(vals)}, max {max(vals)})")
    return mean, f"UK National Grid API, {start} to {end}, {len(vals)} daily means"


def watts(profile, usage=1.0):
    """Instantaneous draw of one allocation unit, before PUE."""
    return (profile["cores"] * profile["w_per_core"] * usage
            + profile["gpus"] * profile["w_per_gpu"]
            + profile["mem_gb"] * W_PER_GB)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pue", type=float, default=1.67,
                   help="Power Usage Effectiveness (default 1.67, the Green "
                        "Algorithms global default; CSD3 publishes none)")
    p.add_argument("--pue-alt", type=float, default=1.20,
                   help="second PUE for the sensitivity line (default 1.20, "
                        "typical modern-facility target)")
    p.add_argument("--ci", type=float, default=None,
                   help="carbon intensity gCO2e/kWh; default is fetched live")
    p.add_argument("--usage", type=float, default=1.0,
                   help="CPU usage factor 0-1 (default 1.0: these are "
                        "compute-bound MD runs holding every reserved core)")
    p.add_argument("--psf", type=float, default=1.0, help="pragmatic scaling factor")
    args = p.parse_args()

    print(__doc__.split("THE FORMULA")[0].rstrip())
    print("=" * 78)

    if args.ci is not None:
        ci, ci_src = args.ci, "supplied on the command line"
    else:
        print("Carbon intensity:")
        ci, ci_src = fetch_uk_ci(CAMPAIGN_START, CAMPAIGN_END)
        if ci is None:
            print(f"  ! {ci_src}; falling back to Green Algorithms UK default")
            ci, ci_src = CI_FALLBACK, "Green Algorithms UK default (2021 vintage)"
    print(f"  using CI = {ci:.1f} gCO2e/kWh  [{ci_src}]")
    print(f"  using PUE = {args.pue}, PSF = {args.psf}, usage factor = {args.usage}")

    print("\nPer-unit draw before PUE (Green Algorithms model):")
    for name, prof in (("icelake node   ", ICELAKE), ("ampere 1/4 node", AMPERE_QUARTER)):
        w = watts(prof, args.usage)
        detail = (f"{prof['cores']}c x {prof['w_per_core']:.3f} W"
                  + (f" + {prof['gpus']} GPU x {prof['w_per_gpu']:.0f} W" if prof["gpus"] else "")
                  + f" + {prof['mem_gb']:.0f} GB x {W_PER_GB} W/GB")
        print(f"  {name}: {w:7.1f} W   ({detail})")

    for scenario, rows in SCENARIOS.items():
        print(f"\n{'=' * 78}\n{scenario}\n{'-' * 78}")
        print(f"{'code':<26}{'hours':>10}  {'unit':<7}{'kWh':>12}{'kgCO2e':>11}"
              f"{'tree-mo':>10}{'car-km':>11}")
        total_kwh = 0.0
        for label, hours, prof, unit in rows:
            kwh = hours * watts(prof, args.usage) / 1000.0 * args.pue * args.psf
            g = kwh * ci
            total_kwh += kwh
            print(f"{label:<26}{hours:>10,}  {unit:<7}{kwh:>12,.0f}{g/1000:>11,.1f}"
                  f"{g/TREE_MONTH_G:>10,.0f}{g/CAR_EU_G_PER_KM:>11,.0f}")
        # avoided vs master (the report's number)
        m_kwh = rows[0][1] * watts(rows[0][2], args.usage) / 1000.0 * args.pue * args.psf
        for label, hours, prof, unit in rows[1:]:
            kwh = hours * watts(prof, args.usage) / 1000.0 * args.pue * args.psf
            saved = m_kwh - kwh
            print(f"  -> {label} avoids {saved:,.0f} kWh / "
                  f"{saved * ci / 1e6:,.2f} tCO2e vs master "
                  f"({saved * ci / FLIGHT_PAR_LON_G:,.0f} Paris-London flights)")

    print(f"\n{'=' * 78}\nSENSITIVITY (both parameters are pure multipliers)")
    print(f"  PUE {args.pue} -> {args.pue_alt}: multiply every number by "
          f"{args.pue_alt / args.pue:.2f}")
    print(f"  CI {ci:.1f} -> {CI_FALLBACK} (GA UK default): multiply carbon by "
          f"{CI_FALLBACK / ci:.2f}")
    print("  Equilibration (30 ps x 3) adds ~4% to the as-run rows if wanted all-in.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
