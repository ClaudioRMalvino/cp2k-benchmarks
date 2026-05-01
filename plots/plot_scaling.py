#!/usr/bin/env python3
import os
import pandas as pd
import matplotlib

matplotlib.use('Agg') 
import matplotlib.pyplot as plt

# --- Define File Paths ---
# Output directory
PLOT_DIR = "/home/raid/crm98/cp2k_benchmarks/plots"

# Size Scaling CSVs
SIZE_MASTER = "/home/raid/crm98/cp2k_benchmarks/cp2k_master/NNP/NNP_size_scaling_upstream-master_30-04_04-01/results_size_scaling_upstream-master_30-04_04-01.csv"
SIZE_OPT = "/home/raid/crm98/cp2k_benchmarks/cp2k_optimized/NNP/NNP_size_scaling_feature-nnp-verlet-cells_30-04_04-02/results_size_scaling_feature-nnp-verlet-cells_30-04_04-02.csv"

# Core Scaling CSVs
CORE_MASTER = "/home/raid/crm98/cp2k_benchmarks/cp2k_master/NNP/NNP_core_scaling_upstream-master_30-04_04-04/results_core_scaling_upstream-master_30-04_04-04.csv"
CORE_OPT = "/home/raid/crm98/cp2k_benchmarks/cp2k_optimized/NNP/NNP_core_scaling_feature-nnp-verlet-cells_30-04_04-04/results_core_scaling_feature-nnp-verlet-cells_30-04_04-04.csv"

# Ensure output directory exists
os.makedirs(PLOT_DIR, exist_ok=True)

# --- Helper Function to Read Data ---
def load_csv(filepath, columns):
    if not os.path.exists(filepath):
        print(f"WARNING: File not found: {filepath}")
        return None
    return pd.read_csv(filepath, comment='#', header=None, names=columns)

# --- 1. Plot Size Scaling ---
print("Plotting Size Scaling...")
df_size_master = load_csv(SIZE_MASTER, ['n_molecules', 'walltime_s'])
df_size_opt = load_csv(SIZE_OPT, ['n_molecules', 'walltime_s'])

if df_size_master is not None and df_size_opt is not None:
    plt.figure(figsize=(8, 6))
    plt.plot(df_size_master['n_molecules'], df_size_master['walltime_s'], marker='o', label='Upstream Master', linestyle='-')
    plt.plot(df_size_opt['n_molecules'], df_size_opt['walltime_s'], marker='s', label='Feature: Verlet Cells', linestyle='-')
    
    plt.title("NNP Size Scaling: Walltime vs. Number of Molecules")
    plt.xlabel("Number of Molecules")
    plt.ylabel("Walltime (seconds)")
    plt.xscale('log', base=2)
    plt.yscale('log')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
    
    out_file = os.path.join(PLOT_DIR, "size_scaling_walltime.png")
    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {out_file}")


# --- 2. Plot Core Scaling ---
print("\nPlotting Core Scaling...")
df_core_master = load_csv(CORE_MASTER, ['mpi_ranks', 'omp_threads', 'total_cores', 'walltime_s', 'speedup'])
df_core_opt = load_csv(CORE_OPT, ['mpi_ranks', 'omp_threads', 'total_cores', 'walltime_s', 'speedup'])

if df_core_master is not None and df_core_opt is not None:
    # Plot A: Walltime vs Cores
    plt.figure(figsize=(8, 6))
    plt.plot(df_core_master['total_cores'], df_core_master['walltime_s'], marker='o', label='Upstream Master')
    plt.plot(df_core_opt['total_cores'], df_core_opt['walltime_s'], marker='s', label='Feature: Verlet Cells')
    
    plt.title("NNP Core Scaling: Walltime vs. Total Cores (64 Molecules)")
    plt.xlabel("Total Cores (MPI * OMP)")
    plt.ylabel("Walltime (seconds)")
    plt.xscale('log', base=2)
    plt.yscale('log')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
    
    out_file_time = os.path.join(PLOT_DIR, "core_scaling_walltime.png")
    plt.savefig(out_file_time, dpi=300, bbox_inches='tight')
    print(f"Saved: {out_file_time}")

    # Plot B: Speedup vs Cores
    plt.figure(figsize=(8, 6))
    plt.plot(df_core_master['total_cores'], df_core_master['speedup'], marker='o', label='Upstream Master')
    plt.plot(df_core_opt['total_cores'], df_core_opt['speedup'], marker='s', label='Feature: Verlet Cells')
    
    # Add Ideal Scaling Line (y = x)
    ideal_cores = df_core_master['total_cores']
    plt.plot(ideal_cores, ideal_cores, 'k--', label='Ideal Linear Scaling')
    
    plt.title("NNP Core Scaling: Speedup vs. Total Cores")
    plt.xlabel("Total Cores")
    plt.ylabel("Speedup (Relative to 1 Core)")
    plt.grid(True, ls="--", alpha=0.5)
    plt.legend()
    
    out_file_speedup = os.path.join(PLOT_DIR, "core_scaling_speedup.png")
    plt.savefig(out_file_speedup, dpi=300, bbox_inches='tight')
    print(f"Saved: {out_file_speedup}")

print("\nAll plotting complete!")