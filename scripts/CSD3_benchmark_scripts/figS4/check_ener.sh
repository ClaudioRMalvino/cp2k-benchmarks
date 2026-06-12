#!/usr/bin/env bash
# Trajectory health check for CP2K .ener files (login-node, free).
#
#   ./check_ener.sh file1.ener [file2.ener ...]
#   ./check_ener.sh /path/to/equil/N*/equil_*.ener
#
# Reports, per file: steps, T range, conserved-quantity drift, and the
# count of single-step jumps: Pot > 0.10 Ha (the May neighbour-list bug
# signature; legitimate NVT relaxation can reach ~0.05) and ConsQty
# > 0.05 Ha (CQ should be near-flat regardless of relaxation).  Exit
# code 1 if any file shows jumps or T > 1000 K, to gate pipelines.

set -u
fail=0
for f in "$@"; do
    [[ -f "$f" ]] || { echo "MISSING  $f"; fail=1; continue; }
    awk -v fname="$f" '
        NR==1 { next }
        {
            t=$2; T=$4; pot=$5; cq=$6; n++
            if (n==1) { cq0=cq; Tmin=T; Tmax=T }
            if (T<Tmin) Tmin=T
            if (T>Tmax) Tmax=T
            if (n>1) {
                dp = pot-prevpot; if (dp<0) dp=-dp
                dc = cq-prevcq;   if (dc<0) dc=-dc
                if (dp>0.10) pj++
                if (dc>0.05) cj++
            }
            prevpot=pot; prevcq=cq
        }
        END {
            if (n==0) { printf "EMPTY    %s\n", fname; exit 1 }
            drift = cq-cq0; if (drift<0) drift=-drift
            bad = (pj>0 || cj>0 || Tmax>1000.0)
            printf "%s  %s\n", (bad ? "*** BAD" : "OK     "), fname
            printf "         steps=%d  T=[%.1f, %.1f] K  |CQ drift|=%.4f Ha  Pot jumps>0.10Ha=%d  CQ jumps>0.05Ha=%d\n", \
                   n, Tmin, Tmax, drift, pj+0, cj+0
            exit bad ? 1 : 0
        }' "$f" || fail=1
done
exit $fail
