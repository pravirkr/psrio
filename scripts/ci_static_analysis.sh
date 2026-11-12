#!/usr/bin/env bash
# Architecture, pre-commit lint/format hooks, and clang-tidy.
#
# PHASE:
#   prebuild  — architecture + pre-commit (no compile DB)
#   postbuild — clang-tidy (needs BUILD_DIR)
#   all       — default; run everything sequentially
#
# CI_AGGREGATE=1: prebuild/postbuild record failures and only exit non-zero on
# postbuild or all after printing the combined summary (for parallel lint job).
#
# Usage:
#   LLVM_VERSION=23 BUILD_DIR=build-ci-tidy scripts/ci_static_analysis.sh

set -uo pipefail

LLVM_VERSION=${LLVM_VERSION:-23}
BUILD_DIR=${BUILD_DIR:-build-ci-tidy}
PHASE=${PHASE:-all}
CI_AGGREGATE=${CI_AGGREGATE:-0}
FAILURES_FILE=${FAILURES_FILE:-}

# Resolve LLVM binaries across Linux (clang-tidy-23) and macOS Homebrew
resolve_llvm_tool() {
    local base=$1
    if command -v "${base}-${LLVM_VERSION}" >/dev/null 2>&1; then
        echo "${base}-${LLVM_VERSION}"
    elif [[ -x "/opt/homebrew/opt/llvm/bin/${base}" ]]; then
        echo "/opt/homebrew/opt/llvm/bin/${base}"
    elif command -v "${base}" >/dev/null 2>&1; then
        echo "${base}"
    else
        echo "${base}-${LLVM_VERSION}"
    fi
}

CLANG_TIDY=${CLANG_TIDY:-$(resolve_llvm_tool clang-tidy)}
RUN_CLANG_TIDY=${RUN_CLANG_TIDY:-$(resolve_llvm_tool run-clang-tidy)}
ARCH_CXX=${ARCH_CXX:-$(resolve_llvm_tool clang++)}

ncpu() {
    if command -v nproc >/dev/null 2>&1; then
        nproc
    else
        sysctl -n hw.ncpu 2>/dev/null || echo 4
    fi
}

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_root"

failures=()

record_failure() {
    local name=$1
    failures+=("$name")
    if [[ -n $FAILURES_FILE ]]; then
        printf '%s\n' "$name" >>"$FAILURES_FILE"
    fi
}

run_step() {
    local name=$1
    shift
    printf '\n======== %s ========\n' "$name"
    if "$@"; then
        printf 'OK: %s\n' "$name"
    else
        printf 'FAILED: %s\n' "$name"
        record_failure "$name"
    fi
}

run_prebuild() {
    run_step "architecture (AGENTS.md, standalone header compile)" \
        env CXX="$ARCH_CXX" scripts/check_architecture.sh

    local pre_commit_hooks=(
        trailing-whitespace
        end-of-file-fixer
        mixed-line-ending
        check-yaml
        check-merge-conflict
        check-added-large-files
        check-case-conflict
        clang-format
        cmake-format
        cmake-lint
    )

    for hook in "${pre_commit_hooks[@]}"; do
        run_step "pre-commit:${hook}" pre-commit run "$hook" --all-files
    done
}

run_postbuild() {
    run_step "clang-tidy (compile database targets)" \
        "$RUN_CLANG_TIDY" -quiet -p "$BUILD_DIR" \
        -clang-tidy-binary "$CLANG_TIDY" \
        -warnings-as-errors='*' -j "$(ncpu)" \
        "(^|/)(header_checks|applications|tests)/"

    # Also run misc-include-cleaner directly on each public header as main file
    # (misc-include-cleaner only inspects main files, not included headers).
    run_step "clang-tidy (public headers include-cleaner)" \
        bash -c '
            set -uo pipefail
            clang_tidy=$1
            build_dir=$2
            extra_inc=()
            for cand in \
                "$build_dir/_deps/fmt-src/include" \
                "$build_dir/_deps/highfive-src/include" \
                /opt/homebrew/include \
                /usr/include/hdf5/serial; do
                [[ -d "$cand" ]] && extra_inc+=("-isystem" "$cand")
            done
            if [[ -d "${PSRDADA_INCLUDE_DIR:-}" ]]; then
                extra_inc+=("-isystem" "${PSRDADA_INCLUDE_DIR}")
            elif [[ -d "$build_dir/_deps/psrdada-src/src" ]]; then
                extra_inc+=("-isystem" "$build_dir/_deps/psrdada-src/src")
            fi
            sdk_flags=()
            if command -v xcrun >/dev/null 2>&1; then
                sdk_flags+=("-isysroot" "$(xcrun --show-sdk-path)")
            fi
            rc=0
            pids=()
            for h in $(find include/psrio -name "*.hpp" | sort); do
                if [[ "$h" == "include/psrio/formats/psrdada.hpp" ]]; then
                    has_dada=0
                    for dir in "${PSRDADA_INCLUDE_DIR:-}" "$build_dir/_deps/psrdada-src/src" /usr/local/include /usr/include; do
                        [[ -n "$dir" && -f "$dir/dada_hdu.h" ]] && has_dada=1
                    done
                    [[ $has_dada -eq 0 ]] && continue
                fi
                "$clang_tidy" -quiet -checks="-*,misc-include-cleaner" \
                    -warnings-as-errors="*" "$h" -- \
                    -xc++ -std=c++20 -Wno-pragma-once-outside-header -DFMT_HEADER_ONLY \
                    -Iinclude ${extra_inc[@]+"${extra_inc[@]}"} ${sdk_flags[@]+"${sdk_flags[@]}"} &
                pids+=($!)
            done
            for pid in "${pids[@]}"; do
                if ! wait "$pid"; then
                    rc=1
                fi
            done
            exit $rc
        ' _ "$CLANG_TIDY" "$BUILD_DIR"
}

print_summary() {
    local -a all_failures=()
    if [[ -n $FAILURES_FILE && -f $FAILURES_FILE ]]; then
        while IFS= read -r line; do
            [[ -n "$line" ]] && all_failures+=("$line")
        done <"$FAILURES_FILE"
    elif ((${#failures[@]} > 0)); then
        all_failures=("${failures[@]}")
    fi

    printf '\n======== summary ========\n'
    if ((${#all_failures[@]} == 0)); then
        echo "All static analysis checks passed."
        return 0
    fi

    echo "The following checks failed:"
    local seen=""
    local name
    for name in "${all_failures[@]}"; do
        if [[ "$seen" != *"|$name|"* ]]; then
            echo "  - $name"
            seen+="|$name|"
        fi
    done
    return 1
}

if [[ -n $FAILURES_FILE && $PHASE == prebuild && $CI_AGGREGATE == 1 ]]; then
    : >"$FAILURES_FILE"
fi

case "$PHASE" in
    prebuild)
        run_prebuild
        if [[ $CI_AGGREGATE == 1 ]]; then
            exit 0
        fi
        print_summary
        ;;
    postbuild)
        run_postbuild
        print_summary
        ;;
    all)
        run_prebuild
        run_postbuild
        print_summary
        ;;
    *)
        echo "Unknown PHASE=$PHASE (expected prebuild, postbuild, or all)" >&2
        exit 2
        ;;
esac
