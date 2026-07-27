#!/usr/bin/env bash
# Shared CSD3 runtime environment for the NaCl(aq) MP2 anchor runs.
# Sourced by run_equil.sh / run_production.sh / the sbatch wrappers.
# Mirrors the proven figS4 job environment (modules + toolchain setup +
# TMPDIR-on-scratch fix from job 30395309's ENOSPC).

. /etc/profile.d/modules.sh
module purge
# module env + master-toolchain setup reference unbound vars -> relax -u
set +u
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
BIN_LABEL="${BIN_LABEL:-dhruv-cell-list}"   # the PR #5295 binary (validated)
BIN="${BIN:-$BIN_ROOT/$BIN_LABEL/cp2k.psmp}"
source "$BIN_ROOT/setup"
set -u

export LD_LIBRARY_PATH="$BIN_ROOT/$BIN_LABEL/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

# bash expands here-docs via TMPDIR; node-local /tmp can be full
export TMPDIR=/rds/user/$USER/hpc-work/tmp
mkdir -p "$TMPDIR"
