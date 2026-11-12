#!/usr/bin/env bash
# Configure is done; run ci-tidy build in parallel with prebuild lint, then clang-tidy.
set -uo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_root"

FAILURES_FILE="${repo_root}/.ci-static-failures"
export FAILURES_FILE
export CI_AGGREGATE=1
export LLVM_VERSION="${LLVM_VERSION:-23}"
export BUILD_DIR="${BUILD_DIR:-build-ci-tidy}"

ncpu() {
    if command -v nproc >/dev/null 2>&1; then
        nproc
    else
        sysctl -n hw.ncpu 2>/dev/null || echo 4
    fi
}

: >"$FAILURES_FILE"

cmake --build "$BUILD_DIR" --parallel "$(ncpu)" &
build_pid=$!

PHASE=prebuild scripts/ci_static_analysis.sh

if ! wait "$build_pid"; then
    printf 'cmake build (%s)\n' "$BUILD_DIR" >>"$FAILURES_FILE"
fi

PHASE=postbuild scripts/ci_static_analysis.sh
rc=$?
rm -f "$FAILURES_FILE"
exit $rc
