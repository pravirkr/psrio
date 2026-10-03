#pragma once

#include "psrio/detail/exceptions.hpp"

#include <cstdint>
#include <limits>

namespace psrio::detail {

/// Sample index reached by moving @p delta samples from @p current.
///
/// @p limit is the last legal index, inclusive (`nsamples()`). A negative
/// @p delta moves backward.
/// @throws ValidationError if the landing index is outside `[0, limit]`.
[[nodiscard]] inline std::uint64_t
apply_skip(std::uint64_t current, std::int64_t delta, std::uint64_t limit) {
    if (delta > 0) {
        const auto step = static_cast<std::uint64_t>(delta);
        if (step > std::numeric_limits<std::uint64_t>::max() - current ||
            current + step > limit) {
            throw ValidationError(
                "psrio: skip lands past the readable samples");
        }
        return current + step;
    }
    if (delta == 0) {
        if (current > limit) {
            throw ValidationError(
                "psrio: skip lands past the readable samples");
        }
        return current;
    }
    if (delta == std::numeric_limits<std::int64_t>::min()) {
        throw ValidationError("psrio: skip lands before the first sample");
    }
    const auto step = static_cast<std::uint64_t>(-delta);
    if (step > current) {
        throw ValidationError("psrio: skip lands before the first sample");
    }
    return current - step;
}

} // namespace psrio::detail
