#pragma once

#include "psrio/detail/exceptions.hpp"

#include <cmath>
#include <cstdint>
#include <fmt/format.h>
#include <string_view>

namespace psrio::detail {

/// Relative tolerance sigpyproc uses in `numpy.allclose` for the sample gap.
inline constexpr double kSampleGapTolerance = 0.1;

/// Absolute tolerance of `numpy.allclose` (its default `atol`).
inline constexpr double kSampleGapAbsolute = 1.0e-8;

/// Require @p next_tstart to be one sample after the last sample of a file.
///
/// The gap is `(next_tstart - tend) * 86400`, with
/// `tend = tstart + (nsamples - 1) * tsamp / 86400`. It must match @p tsamp
/// within `1e-8 + 0.1 * abs(tsamp)`.
/// @throws FormatError when the times are not finite or the gap is not one
///         sample.
inline void require_sample_gap(double tstart,
                               std::uint64_t nsamples,
                               double tsamp,
                               double next_tstart,
                               std::string_view path) {
    if (!(tsamp > 0.0) || !std::isfinite(tsamp) || !std::isfinite(tstart) ||
        !std::isfinite(next_tstart)) {
        throw FormatError(fmt::format(
            "psrio: {} needs finite tstart and tsamp to continue a stream",
            path));
    }
    const double end =
        tstart + ((static_cast<double>(nsamples) - 1.0) * tsamp / 86400.0);
    const double gap = (next_tstart - end) * 86400.0;
    const double limit =
        kSampleGapAbsolute + (kSampleGapTolerance * std::abs(tsamp));
    if (!(std::abs(gap - tsamp) <= limit)) {
        throw FormatError(fmt::format(
            "psrio: {} is not contiguous; the gap is {} s and one sample is "
            "{} s",
            path, gap, tsamp));
    }
}

} // namespace psrio::detail
