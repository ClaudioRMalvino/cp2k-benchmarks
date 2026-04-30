#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SCRATCH_USER="crm98"
SCRATCH_ROOT="${SCRATCH_ROOT:-/local/data/public}"
SRC_HOME="/home/raid/crm98/cp2k"
SCRATCH_REPO="/local/data/public/crm98/original_cp2k"
TMPDIR_BASE="${SCRATCH_ROOT}/${SCRATCH_USER}/tmp"
PIP_CACHE_BASE="${SCRATCH_ROOT}/${SCRATCH_USER}/pip-cache"

echo "SCRATCH_USER:  $SCRATCH_USER"
echo "SCRATCH_ROOT:  $SCRATCH_ROOT"
echo "SCRATCH_REPO:  $SCRATCH_REPO"
echo "SRC_HOME:      $SRC_HOME"
echo

mkdir -p "$TMPDIR_BASE" "$PIP_CACHE_BASE"
mkdir -p "${SCRATCH_ROOT}/${SCRATCH_USER}"

if [ ! -d "$SCRATCH_REPO/.git" ]; then
  git clone "$SRC_HOME" "$SCRATCH_REPO"
else
  rsync -a --delete --exclude '.git' "$SRC_HOME/" "$SCRATCH_REPO/"
fi

export TMPDIR="$TMPDIR_BASE"
export PIP_CACHE_DIR="$PIP_CACHE_BASE"

cd "$SCRATCH_REPO/tools/toolchain"

if [ -f "$SCRATCH_REPO/tools/toolchain/install/setup" ]; then
  echo "Toolchain already installed in scratch."
  exit 0
fi

./install_cp2k_toolchain.sh \
  --install-all \
  --with-gcc=system \
  --with-openmpi=system \
  --with-openblas=system \
  --math-mode=openblas \
  --with-elpa=no \
  --with-sirius=no
