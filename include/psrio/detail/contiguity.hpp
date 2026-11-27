#pragma once

#include "psrio/detail/exceptions.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <string_view>

namespace psrio::detail {

/// Live file-level progress updated by lazy multi-file streams.
struct StreamProgress {
    std::size_t file_index{0};
    std::uint64_t file_first_sample{0};
    std::uint64_t file_samples{0};
    std::uint64_t bytes_before_file{0};
    std::uint64_t file_bytes{0};
};

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

/// Exact split-MJD contiguity check for PSRFITS (`STT_IMJD`, `STT_SMJD`,
/// `STT_OFFS`). Avoids catastrophic cancellation from subtracting two ~60000
/// MJD doubles (~0.63 us ULP), preserving sub-nanosecond precision even for
/// microsecond sampling intervals.
/// @throws FormatError when times are not finite or the gap does not match one
///         sample.
inline void require_split_mjd_gap(std::int64_t imjd1,
                                  double smjd1,
                                  double offs1,
                                  std::uint64_t nsamples,
                                  double tsamp,
                                  std::int64_t imjd2,
                                  double smjd2,
                                  double offs2,
                                  std::string_view path) {
    if (!(tsamp > 0.0) || !std::isfinite(tsamp) || !std::isfinite(smjd1) ||
        !std::isfinite(offs1) || !std::isfinite(smjd2) ||
        !std::isfinite(offs2)) {
        throw FormatError(fmt::format(
            "psrio: {} needs finite start time and tsamp to continue a stream",
            path));
    }
    const double delta = (static_cast<double>(imjd2 - imjd1) * 86400.0) +
                         (smjd2 - smjd1) + (offs2 - offs1);
    const double span  = static_cast<double>(nsamples) * tsamp;
    const double error = delta - span;
    const double limit =
        kSampleGapAbsolute + (kSampleGapTolerance * std::abs(tsamp));
    if (!(std::abs(error) <= limit)) {
        throw FormatError(fmt::format(
            "psrio: {} is not contiguous; the gap is {} s and one sample is "
            "{} s",
            path, error + tsamp, tsamp));
    }
}

} // namespace psrio::detail
