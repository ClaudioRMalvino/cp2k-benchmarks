#!/bin/bash
#! Rebuild all three CP2K branches with -g -xCORE-AVX512.
#! Caches each binary into $BIN_ROOT immediately after the build so the master
#! and optimized builds don't trample each other's $INSTALL_DIR (both scripts
#! write to /home/crm98/cp2k_optimized/install).

set -e
exec 2>&1

BIN_ROOT=/rds/user/$USER/hpc-work/cp2k_binaries/csd3
mkdir -p "$BIN_ROOT/master/lib" \
         "$BIN_ROOT/feature-nnp-native-spline/lib" \
         "$BIN_ROOT/feature-nnp-native-spline-omp/lib"

CP2K_OPT=/home/crm98/cp2k_optimized
INSTALL_BIN=$CP2K_OPT/install/bin/cp2k.psmp

log() { echo "[$(date +%H:%M:%S)] $*"; }

cache_binary() {
   local label=$1
   if [[ ! -x "$INSTALL_BIN" ]]; then
      log "!! $label: $INSTALL_BIN missing - build failed?"
      exit 1
   fi
   cp "$INSTALL_BIN" "$BIN_ROOT/$label/cp2k.psmp"
   log "    cached $label : $(md5sum $BIN_ROOT/$label/cp2k.psmp | cut -c1-32)"
}

log "==== building master ===="
SKIP_REGTEST=1 bash /home/crm98/cp2k-benchmarks/cp2k_master/CSD3_build_scripts/cp2k_CSD3_master_build.sh
cache_binary master

log "==== building feature/nnp-native-spline (with regtest) ===="
cd $CP2K_OPT
git checkout feature/nnp-native-spline
bash /home/crm98/cp2k-benchmarks/cp2k_optimized/CSD3_build_scripts/cp2k_CSD3_opt_build.sh
cache_binary feature-nnp-native-spline

log "==== building feature/nnp-native-spline-omp (incremental, SKIP_REGTEST) ===="
cd $CP2K_OPT
git checkout feature/nnp-native-spline-omp
source /home/crm98/cp2k-benchmarks/scripts/CSD3_benchmark_scripts/cp2k_CSD3_env.sh
cmake --build build -j 16
cmake --install build
cache_binary feature-nnp-native-spline-omp

cp /home/crm98/cp2k_master/tools/toolchain/install/setup "$BIN_ROOT/setup"

log "==== ALL THREE BUILDS DONE ===="
echo
md5sum "$BIN_ROOT/master/cp2k.psmp" \
       "$BIN_ROOT/feature-nnp-native-spline/cp2k.psmp" \
       "$BIN_ROOT/feature-nnp-native-spline-omp/cp2k.psmp"
echo
echo "=== branch-specific symbols in each binary ==="
for b in master feature-nnp-native-spline feature-nnp-native-spline-omp; do
   exe="$BIN_ROOT/$b/cp2k.psmp"
   printf "%-32s  hermite_spline_value=%d  scale_factor=%d  sort_active_atoms=%d\n" \
      "$b" \
      "$(nm -a $exe 2>/dev/null | grep -c 'hermite_spline_value')" \
      "$(nm -a $exe 2>/dev/null | grep -c 'scale_factor')" \
      "$(nm -a $exe 2>/dev/null | grep -c 'sort_active_atoms')"
done
