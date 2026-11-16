#include "cli.hpp"

#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <span>
#include <string>
#include <vector>

#if defined(__linux__) && defined(_OPENMP)
#include <stdlib.h> // NOLINT(modernize-deprecated-headers): POSIX setenv
#include <unistd.h>
#endif

namespace {

/// psrio spends much of its time waiting on file reads. OpenMP's default
/// policy keeps idle threads spinning through those waits, which multiplies
/// the CPU time a batch job is charged for no gain. The runtime reads
/// OMP_WAIT_POLICY before main(), so when the user has not chosen a policy,
/// run again with a passive one. A failed exec leaves the defaults in place.
void prefer_passive_openmp_threads([[maybe_unused]] char* const* argv) {
#if defined(__linux__) && defined(_OPENMP)
    // NOLINTBEGIN(concurrency-mt-unsafe): single-threaded at this point
    if (std::getenv("OMP_WAIT_POLICY") != nullptr ||
        std::getenv("GOMP_SPINCOUNT") != nullptr) {
        return;
    }
    if (::setenv("OMP_WAIT_POLICY", "passive", 0) != 0) {
        return;
    }
    // NOLINTEND(concurrency-mt-unsafe)
    // Re-runs this same program with its own, unchanged arguments.
    // NOLINTNEXTLINE(clang-analyzer-optin.taint.GenericTaint)
    (void)::execv("/proc/self/exe", argv);
#endif
}

} // namespace

int main(int argc, char** argv) {
    prefer_passive_openmp_threads(argv);
    try {
        const std::span<char*> raw(argv, static_cast<std::size_t>(argc));
        const std::vector<std::string> args(raw.begin() + 1, raw.end());
        return psrio::cli::run(args, std::cout, std::cerr);
    } catch (const std::exception& error) {
        std::cerr << "psrio: fatal: " << error.what() << '\n';
    } catch (...) {
        std::cerr << "psrio: fatal: unknown error\n";
    }
    return 1;
}
