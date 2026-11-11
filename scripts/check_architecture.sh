#!/usr/bin/env bash
# Mechanical checks for the architectural rules in AGENTS.md.
#
# Usage: scripts/check_architecture.sh
#   CXX=<compiler>              compiler for the standalone public-header check (default c++)
#   SKIP_COMPILE=1              skip the standalone public-header check
#   REQUIRE_OPTIONAL_HEADERS=1  fail if optional headers (fbh5.hpp, psrdada.hpp) cannot be compiled
#
# Exits non-zero and prints every violation found.

set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 2

status=0
fail() {
    printf '  %s\n' "$@"
    status=1
}
section() { printf '%s\n' "$1"; }

# All C++ sources under the given directories.
sources() { find "$@" -type f \( -name '*.cpp' -o -name '*.hpp' \) 2>/dev/null; }

# 1. Base headers stay free of heavy optional dependencies (AGENTS.md).
section "base headers do not pull in optional heavy dependencies"
while IFS= read -r line; do fail "$line"; done < <(
    find include/psrio -name '*.hpp' \
        ! -path 'include/psrio/formats/fbh5.hpp' \
        ! -path 'include/psrio/formats/psrdada.hpp' |
        xargs grep -nE '#include [<"](hdf5|H5|highfive|dada_|ipcbuf|ipcio|ascii_header|multilog|fitsio|dmt|CLI/)'
)

section "base headers do not include opt-in format headers (fbh5.hpp, psrdada.hpp)"
while IFS= read -r line; do fail "$line"; done < <(
    find include/psrio -name '*.hpp' \
        ! -path 'include/psrio/formats/fbh5.hpp' \
        ! -path 'include/psrio/formats/psrdada.hpp' |
        xargs grep -nE '#include [<"]psrio/formats/(fbh5|psrdada)\.hpp[>"]'
)

# 2. Every public header uses #pragma once (AGENTS.md Coding Standards rule 4).
section "all public headers use #pragma once"
for h in $(find include/psrio -name '*.hpp' | sort); do
    if ! grep -q '^#pragma once$' "$h"; then
        fail "$h: missing #pragma once"
    fi
done

# 3. Include prefixes: repo headers with a slash must be "psrio/...".
section "include prefixes are psrio/"
while IFS= read -r line; do fail "$line"; done < <(
    sources include applications tests |
        xargs grep -nE '#include "[^"]*/' |
        grep -vE '#include "psrio/'
)

# 4. fmt::format, never std::format: libstdc++ < 13 has no <format> (AGENTS.md).
section "no std::format or <format> (use fmt::format)"
while IFS= read -r line; do fail "$line"; done < <(
    sources include applications tests |
        xargs grep -nE 'std::(v?format|format_to|format_string|formatter)\b|#include <format>'
)

# 5. Public headers open only expected psrio namespaces.
section "public namespaces follow the directory"
common_ns='psrio|psrio::detail|psrio::concepts|psrio::error_check|detail|concepts|[a-z_]+_detail'
for h in $(find include/psrio -name '*.hpp' | sort); do
    rel=${h#include/psrio/}
    dir=$(dirname "$rel")
    case $dir in
        .) ok="${common_ns}|psrio::astro|psrio::digifil" ;;
        common | detail) ok="${common_ns}" ;;
        formats) ok="${common_ns}|psrio::formats::[a-z0-9_]+|psrio::formats::[a-z0-9_]+::detail" ;;
        formats/*)
            sub=${dir#formats/}
            ok="${common_ns}|psrio::formats::${sub}|psrio::formats::${sub}::detail"
            ;;
        *) ok="${common_ns}" ;;
    esac
    while IFS= read -r ns; do fail "$h: namespace $ns"; done < <(
        grep -oE '^namespace [A-Za-z0-9_:]+' "$h" | awk '{print $2}' | grep -vxE "$ok"
    )
done

# 6. Every public header compiles on its own.
if [[ "${SKIP_COMPILE:-0}" != 1 ]]; then
    section "public headers compile standalone"
    cxx=${CXX:-c++}

    extra_includes=()
    # Discover fmt include path if not in default system path
    for cand in \
        "${FMT_INCLUDE_DIR:-}" \
        build/_deps/fmt-src/include \
        build-psrfits/_deps/fmt-src/include \
        build-ci-tidy/_deps/fmt-src/include \
        /opt/homebrew/include; do
        if [[ -n "$cand" && -f "$cand/fmt/format.h" ]]; then
            extra_includes+=("-I$cand")
            break
        fi
    done

    # Discover HDF5 and HighFive include paths
    hdf5_inc=""
    for cand in \
        "${HDF5_INCLUDE_DIR:-}" \
        /opt/homebrew/include \
        /usr/include/hdf5/serial \
        /usr/include; do
        if [[ -n "$cand" && -f "$cand/hdf5.h" ]]; then
            hdf5_inc="$cand"
            break
        fi
    done

    highfive_inc=""
    for cand in \
        "${HIGHFIVE_INCLUDE_DIR:-}" \
        build/_deps/highfive-src/include \
        build-psrfits/_deps/highfive-src/include \
        build-ci-tidy/_deps/highfive-src/include \
        /opt/homebrew/include \
        /usr/local/include \
        /usr/include; do
        if [[ -n "$cand" && -f "$cand/highfive/highfive.hpp" ]]; then
            highfive_inc="$cand"
            break
        fi
    done

    # Discover PSRDADA include path
    psrdada_inc=""
    for cand in \
        "${PSRDADA_INCLUDE_DIR:-}" \
        build/_deps/psrdada-src/src \
        build-psrfits/_deps/psrdada-src/src \
        build-ci-tidy/_deps/psrdada-src/src \
        /usr/local/include \
        /usr/include; do
        if [[ -n "$cand" && -f "$cand/dada_hdu.h" ]]; then
            psrdada_inc="$cand"
            break
        fi
    done

    for h in $(cd include && find psrio -name '*.hpp' | sort); do
        header_flags=()
        if ((${#extra_includes[@]} > 0)); then
            header_flags+=("${extra_includes[@]}")
        fi
        if [[ "$h" == "psrio/formats/fbh5.hpp" ]]; then
            if [[ -n "$hdf5_inc" && -n "$highfive_inc" ]]; then
                header_flags+=("-I$hdf5_inc" "-I$highfive_inc")
            elif [[ "${REQUIRE_OPTIONAL_HEADERS:-0}" == 1 ]]; then
                fail "$h: hdf5.h / highfive.hpp not found and REQUIRE_OPTIONAL_HEADERS=1"
                continue
            else
                continue
            fi
        elif [[ "$h" == "psrio/formats/psrdada.hpp" ]]; then
            if [[ -n "$psrdada_inc" ]]; then
                header_flags+=("-I$psrdada_inc")
            elif [[ "${REQUIRE_OPTIONAL_HEADERS:-0}" == 1 ]]; then
                fail "$h: dada_hdu.h not found and REQUIRE_OPTIONAL_HEADERS=1"
                continue
            else
                continue
            fi
        fi

        if ! out=$(echo "#include \"$h\"" |
            "$cxx" -std=c++20 -Wall -Wextra -Wpedantic -Werror -fsyntax-only \
                -DFMT_HEADER_ONLY -Iinclude ${header_flags[@]+"${header_flags[@]}"} -x c++ - 2>&1); then
            fail "$h: $(echo "$out" | grep -m1 error)"
        fi
    done
fi

if [[ $status -eq 0 ]]; then
    echo "architecture checks passed"
else
    echo "architecture checks FAILED (see AGENTS.md)"
fi
exit $status
