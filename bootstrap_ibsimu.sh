#!/bin/bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IBSIMU_DIR="${ROOT_DIR}/libibsimu_patched"
PATCH_DIR="${ROOT_DIR}/patches/ibsimu"
STAMP_FILE="${IBSIMU_DIR}/.cpi-bootstrap-stamp"

UPSTREAM_URL="https://git.code.sf.net/p/ibsimu/code"
UPSTREAM_COMMIT="e8500a2"
PATCHES=(
    "0001-local-ibsimu-customizations.patch"
    "0002-upstream-missing-source-files.patch"
)

require_tool() {
    local tool_name="$1"
    if ! command -v "$tool_name" >/dev/null 2>&1; then
        echo "Missing required tool: $tool_name" >&2
        exit 1
    fi
}

bootstrap_manifest() {
    printf 'upstream_url=%s\n' "$UPSTREAM_URL"
    printf 'upstream_commit=%s\n' "$UPSTREAM_COMMIT"
    printf 'patches=%s\n' "$(IFS=,; echo "${PATCHES[*]}")"
}

current_manifest="$(bootstrap_manifest)"

if [[ -f "$STAMP_FILE" ]] && [[ "$(cat "$STAMP_FILE")" == "$current_manifest" ]]; then
    echo "libibsimu_patched is already bootstrapped"
    exit 0
fi

require_tool git
require_tool autoreconf
require_tool autoconf
require_tool automake
require_tool aclocal
require_tool libtoolize

for patch_file in "${PATCHES[@]}"; do
    if [[ ! -f "${PATCH_DIR}/${patch_file}" ]]; then
        echo "Missing IBSimu patch file: ${PATCH_DIR}/${patch_file}" >&2
        exit 1
    fi
done

echo "Bootstrapping libibsimu_patched from upstream ${UPSTREAM_COMMIT}"
rm -rf "$IBSIMU_DIR"

git clone "$UPSTREAM_URL" "$IBSIMU_DIR"

pushd "$IBSIMU_DIR" >/dev/null
git checkout "$UPSTREAM_COMMIT"

# Upstream automake expects this file to exist before autoreconf.
touch ChangeLog

for patch_file in "${PATCHES[@]}"; do
    echo "Applying ${patch_file}"
    git apply "${PATCH_DIR}/${patch_file}"
done

autoreconf -fi
./configure
rm -rf .git

printf '%s\n' "$current_manifest" > "$STAMP_FILE"
popd >/dev/null

echo "Finished bootstrapping libibsimu_patched"