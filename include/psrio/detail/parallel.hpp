#pragma once

#include <cstdint>

#ifdef _OPENMP
#include <omp.h>
#endif

/// `PSRIO_OMP(omp parallel for ...)` expands to `#pragma omp parallel for ...`
/// when the including translation unit is compiled with OpenMP, and to
/// nothing otherwise. psrio stays header-only with no OpenMP dependency; a
/// consumer that builds with `-fopenmp` (or links `OpenMP::OpenMP_CXX`) gets
/// the parallel loops. Every parallel loop in psrio writes disjoint outputs
/// and keeps each output's arithmetic in serial order, so results do not
/// depend on the thread count.
#ifdef _OPENMP
#define PSRIO_OMP(directive) _Pragma(#directive)
#else
#define PSRIO_OMP(directive)
#endif

namespace psrio::detail {

/// Work, in values touched, below which a loop stays on one thread.
inline constexpr std::uint64_t kParallelMinValues = std::uint64_t{1} << 16U;

/// True when a loop over @p values values is worth splitting across threads.
[[nodiscard]] constexpr bool
parallel_worthwhile(std::uint64_t values) noexcept {
    return values >= kParallelMinValues;
}

} // namespace psrio::detail
