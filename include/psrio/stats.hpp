#pragma once

#include "psrio/block_source.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/float_blocks.hpp"
#include "psrio/detail/parallel.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <limits>
#include <optional>
#include <span>
#include <vector>

namespace psrio::stats {

/// Moments a ChannelStats accumulates.
enum class Moments : std::uint8_t {
    /// Mean, minimum, and maximum (a bandpass).
    kMean,
    /// Up to the fourth central moment (variance, skewness, kurtosis).
    kFull,
};

/// Per-channel moments of a [sample, channel] stream, accumulated block by
/// block in one pass.
///
/// Each block's moments come from double-precision power sums about a
/// per-channel shift (the running mean, or the first sample), and are
/// merged into the running moments with the pairwise update of Pébay
/// (2008). Definitions follow sigpyproc's `ChannelStats`: variance is
/// M2/n, skewness is sqrt(n) M3 / M2^1.5, and kurtosis is the excess
/// n M4 / M2^2 - 3. Columns are split across threads in tiles, so results
/// do not depend on the thread count.
class ChannelStats {
public:
    explicit ChannelStats(std::uint64_t ncols, Moments order = Moments::kFull)
        : m_ncols(ncols),
          m_order(order),
          m_mean(static_cast<std::size_t>(ncols), 0.0),
          m_min(static_cast<std::size_t>(ncols),
                std::numeric_limits<float>::infinity()),
          m_max(static_cast<std::size_t>(ncols),
                -std::numeric_limits<float>::infinity()) {
        if (order == Moments::kFull) {
            m_m2.assign(static_cast<std::size_t>(ncols), 0.0);
            m_m3.assign(static_cast<std::size_t>(ncols), 0.0);
            m_m4.assign(static_cast<std::size_t>(ncols), 0.0);
        }
    }

    [[nodiscard]] std::uint64_t ncols() const noexcept { return m_ncols; }
    [[nodiscard]] Moments order() const noexcept { return m_order; }
    /// Samples accumulated in every column.
    [[nodiscard]] std::uint64_t count() const noexcept { return m_count; }

    /// Add @p nrows samples of ncols() values each, row-major.
    /// @throws ValidationError if @p data does not hold nrows * ncols values.
    void push(std::span<const float> data, std::uint64_t nrows) {
        if (data.size() != static_cast<std::size_t>(nrows * m_ncols)) {
            throw ValidationError(fmt::format(
                "psrio: statistics block has {} values, expected {}",
                data.size(), nrows * m_ncols));
        }
        if (nrows == 0U) {
            return;
        }
        const auto ntiles = (m_ncols + kTile - 1U) / kTile;
        [[maybe_unused]] const bool wide =
            ::psrio::detail::parallel_worthwhile(nrows * m_ncols);
        const bool full = m_order == Moments::kFull;
        PSRIO_OMP(omp parallel for default(none) schedule(static) if (wide)
                      shared(ntiles, full, data, nrows))
        for (std::uint64_t tile = 0; tile < ntiles; ++tile) {
            if (full) {
                push_tile<true>(data, nrows, tile * kTile);
            } else {
                push_tile<false>(data, nrows, tile * kTile);
            }
        }
        m_count += nrows;
    }

    /// Combine with statistics of other samples of the same columns.
    /// @throws ValidationError if the shapes or orders differ.
    void merge(const ChannelStats& other) {
        if (other.m_ncols != m_ncols || other.m_order != m_order) {
            throw ValidationError(
                "psrio: statistics of different shapes do not merge");
        }
        for (std::size_t col = 0; col < static_cast<std::size_t>(m_ncols);
             ++col) {
            Block block{.mean = other.m_mean[col]};
            if (m_order == Moments::kFull) {
                block.m2 = other.m_m2[col];
                block.m3 = other.m_m3[col];
                block.m4 = other.m_m4[col];
            }
            combine(col, other.m_count, block);
            m_min[col] = std::min(m_min[col], other.m_min[col]);
            m_max[col] = std::max(m_max[col], other.m_max[col]);
        }
        m_count += other.m_count;
    }

    [[nodiscard]] std::span<const double> mean() const noexcept {
        return m_mean;
    }
    [[nodiscard]] std::span<const float> minima() const noexcept {
        return m_min;
    }
    [[nodiscard]] std::span<const float> maxima() const noexcept {
        return m_max;
    }
    /// Sum of each column: mean * count.
    [[nodiscard]] std::vector<double> sum() const {
        std::vector<double> out(m_mean.size());
        for (std::size_t col = 0; col < out.size(); ++col) {
            out[col] = m_mean[col] * static_cast<double>(m_count);
        }
        return out;
    }
    /// Population variance, M2 / n (0 with no samples).
    [[nodiscard]] std::vector<double> var() const {
        require_full();
        std::vector<double> out(m_m2.size(), 0.0);
        if (m_count == 0U) {
            return out;
        }
        for (std::size_t col = 0; col < out.size(); ++col) {
            out[col] = m_m2[col] / static_cast<double>(m_count);
        }
        return out;
    }
    [[nodiscard]] std::vector<double> stddev() const {
        auto out = var();
        for (auto& value : out) {
            value = std::sqrt(value);
        }
        return out;
    }
    /// sqrt(n) M3 / M2^1.5 (0 for a constant column).
    [[nodiscard]] std::vector<double> skew() const {
        require_full();
        std::vector<double> out(m_m2.size(), 0.0);
        const double root_n = std::sqrt(static_cast<double>(m_count));
        for (std::size_t col = 0; col < out.size(); ++col) {
            if (m_m2[col] != 0.0) {
                out[col] = root_n * m_m3[col] / std::pow(m_m2[col], 1.5);
            }
        }
        return out;
    }
    /// Excess kurtosis n M4 / M2^2 - 3 (-3 for a constant column, as
    /// sigpyproc).
    [[nodiscard]] std::vector<double> kurtosis() const {
        require_full();
        std::vector<double> out(m_m2.size(), 0.0);
        if (m_count == 0U) {
            return out;
        }
        const auto n = static_cast<double>(m_count);
        for (std::size_t col = 0; col < out.size(); ++col) {
            const double ratio =
                m_m2[col] != 0.0 ? m_m4[col] / (m_m2[col] * m_m2[col]) : 0.0;
            out[col] = (n * ratio) - 3.0;
        }
        return out;
    }

private:
    /// Columns one thread accumulates together.
    static constexpr std::uint64_t kTile = 64;

    /// Central moments of one block of one column.
    struct Block {
        double mean{0.0};
        double m2{0.0};
        double m3{0.0};
        double m4{0.0};
    };

    void require_full() const {
        if (m_order != Moments::kFull) {
            throw ValidationError(
                "psrio: these statistics were accumulated as a mean only");
        }
    }

    template <bool Full>
    void push_tile(std::span<const float> data,
                   std::uint64_t nrows,
                   std::uint64_t first) {
        const auto width = static_cast<std::size_t>(
            std::min<std::uint64_t>(kTile, m_ncols - first));
        const auto col0  = static_cast<std::size_t>(first);
        const auto ncols = static_cast<std::size_t>(m_ncols);
        std::array<double, kTile> shift{};
        std::array<double, kTile> s1{};
        std::array<double, kTile> s2{}; // NOLINT(misc-const-correctness)
        std::array<double, kTile> s3{}; // NOLINT(misc-const-correctness)
        std::array<double, kTile> s4{}; // NOLINT(misc-const-correctness)
        std::array<float, kTile> lo{};
        std::array<float, kTile> hi{};
        for (std::size_t j = 0; j < width; ++j) {
            shift[j] = m_count > 0U ? m_mean[col0 + j]
                                    : static_cast<double>(data[col0 + j]);
            lo[j]    = m_min[col0 + j];
            hi[j]    = m_max[col0 + j];
        }
        const float* base = data.data() + col0;
        for (std::uint64_t row = 0; row < nrows; ++row) {
            const float* values =
                base + (static_cast<std::size_t>(row) * ncols);
            PSRIO_OMP(omp simd)
            for (std::size_t j = 0; j < width; ++j) {
                const float value = values[j];
                const double d    = static_cast<double>(value) - shift[j];
                s1[j] += d;
                if constexpr (Full) {
                    const double d2 = d * d;
                    s2[j] += d2;
                    s3[j] += d2 * d;
                    s4[j] += d2 * d2;
                }
                lo[j] = std::min(lo[j], value);
                hi[j] = std::max(hi[j], value);
            }
        }
        const auto n = static_cast<double>(nrows);
        for (std::size_t j = 0; j < width; ++j) {
            const double m = s1[j] / n;
            // NOLINTNEXTLINE(misc-const-correctness): set when Full
            Block block{.mean = shift[j] + m};
            if constexpr (Full) {
                const double m_sq = m * m;
                block.m2          = s2[j] - (n * m_sq);
                block.m3 = s3[j] - (3.0 * m * s2[j]) + (2.0 * n * m_sq * m);
                block.m4 = s4[j] - (4.0 * m * s3[j]) + (6.0 * m_sq * s2[j]) -
                           (3.0 * n * m_sq * m_sq);
                block.m2 = std::max(block.m2, 0.0);
                block.m4 = std::max(block.m4, 0.0);
            }
            combine(col0 + j, nrows, block);
            m_min[col0 + j] = lo[j];
            m_max[col0 + j] = hi[j];
        }
    }

    /// Merge @p nb samples with moments @p block into column @p col, which
    /// holds m_count samples.
    void combine(std::size_t col, std::uint64_t nb_count, const Block& block) {
        if (nb_count == 0U) {
            return;
        }
        if (m_count == 0U) {
            m_mean[col] = block.mean;
            if (m_order == Moments::kFull) {
                m_m2[col] = block.m2;
                m_m3[col] = block.m3;
                m_m4[col] = block.m4;
            }
            return;
        }
        const auto na    = static_cast<double>(m_count);
        const auto nb    = static_cast<double>(nb_count);
        const double n   = na + nb;
        const double del = block.mean - m_mean[col];
        m_mean[col] += del * nb / n;
        if (m_order != Moments::kFull) {
            return;
        }
        const double m2a  = m_m2[col];
        const double m3a  = m_m3[col];
        const double del2 = del * del;
        const double ab   = na * nb;
        m_m4[col] +=
            block.m4 +
            (del2 * del2 * ab * ((na * na) - ab + (nb * nb)) / (n * n * n)) +
            (6.0 * del2 * ((na * na * block.m2) + (nb * nb * m2a)) / (n * n)) +
            (4.0 * del * ((na * block.m3) - (nb * m3a)) / n);
        m_m3[col] += block.m3 + (del2 * del * ab * (na - nb) / (n * n)) +
                     (3.0 * del * ((na * block.m2) - (nb * m2a)) / n);
        m_m2[col] += block.m2 + (del2 * ab / n);
    }

    std::uint64_t m_ncols;
    Moments m_order;
    std::uint64_t m_count{0};
    std::vector<double> m_mean;
    std::vector<double> m_m2;
    std::vector<double> m_m3;
    std::vector<double> m_m4;
    std::vector<float> m_min;
    std::vector<float> m_max;
};

/// Which samples compute() reads.
struct ComputeOptions {
    std::uint64_t start{0};
    /// Samples to read; all from start when empty.
    std::optional<std::uint64_t> nsamps;
    /// Samples per block.
    std::uint64_t gulp{16384};
    Moments order{Moments::kFull};
};

/// Statistics of every channel (and IF) of @p source, read once with the
/// next block read while the current one is accumulated. Leaves the source
/// positioned after the samples read.
[[nodiscard]] inline ChannelStats compute(BlockSource& source,
                                          const ComputeOptions& options = {}) {
    detail::FloatBlocks blocks(source, options.start, options.nsamps,
                               options.gulp);
    ChannelStats stats(blocks.row(), options.order);
    while (true) {
        const auto block = blocks.next();
        if (block.empty()) {
            break;
        }
        stats.push(block, block.size() / blocks.row());
    }
    return stats;
}

} // namespace psrio::stats
