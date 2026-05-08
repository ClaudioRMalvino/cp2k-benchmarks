#!/usr/bin/env bash

cp2k_known_scratch_roots() {
  local scratch_user="crm98"
  local scratch_root="${SCRATCH_ROOT:-/local/data/public}"

  printf '%s\n' \
    "$scratch_root/$scratch_user/new_cp2k" \
    "$scratch_root/$scratch_user/original_cp2k" \
    "$scratch_root/$scratch_user/cp2k-buildtree" \
    "/local/data/public/$scratch_user/new_cp2k" \
    "/local/data/public/$scratch_user/original_cp2k" \
    "/local/data/public/$scratch_user/cp2k-buildtree" \
    "/data/lovelace/public/$scratch_user/new_cp2k" \
    "/data/lovelace/public/$scratch_user/original_cp2k" \
    "/data/lovelace/public/$scratch_user/cp2k-buildtree" \
    "/data/lovelace/$scratch_user/new_cp2k" \
    "/data/lovelace/$scratch_user/original_cp2k" \
    "/data/lovelace/$scratch_user/cp2k-buildtree"
}

cp2k_source_toolchain_setup() {
  local cp2k_root="$1"
  local setup_path="$cp2k_root/tools/toolchain/install/setup"
  local temp_setup=""
  local needs_rewrite=0
  local known_root
  local -a known_roots=()

  if [ ! -f "$setup_path" ]; then
    echo "Missing scratch toolchain setup:" >&2
    echo "  $setup_path" >&2
    return 1
  fi

  while IFS= read -r known_root; do
    known_roots+=("$known_root")
  done < <(cp2k_known_scratch_roots)

  for known_root in "${known_roots[@]}"; do
    if [ "$cp2k_root" != "$known_root" ] && grep -Fq "$known_root" "$setup_path"; then
      needs_rewrite=1
      break
    fi
  done

  if [ "$needs_rewrite" -eq 1 ]; then
    temp_setup="$(mktemp "${TMPDIR:-/tmp}/cp2k_toolchain_setup.XXXXXX")"
    {
      printf 'sed'
      for known_root in "${known_roots[@]}"; do
        printf ' -e %q' "s|$known_root|$cp2k_root|g"
      done
      printf ' %q\n' "$setup_path"
    } | sh > "$temp_setup"

    set +u
    source "$temp_setup"
    set -u
    rm -f "$temp_setup"
    return 0
  fi

  set +u
  source "$setup_path"
  set -u
}
