#pragma once

#include "psrio/detail/exceptions.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace psrio {

/// Sum consecutive groups of @p factor channels, in place.
///
/// @p data holds @p nrows rows of @p nchans values. Row r of the result,
/// nchans / factor values, is written to the front of @p data, so the result
/// is the first nrows * nchans / factor values. Each group is summed in
/// channel order in float.
/// @throws ValidationError if @p factor is 0 or does not divide @p nchans, or
///         @p data is too small.
inline std::span<float> fscrunch_inplace(std::span<float> data,
                                         std::uint64_t nrows,
                                         std::uint64_t nchans,
                                         std::uint64_t factor) {
    if (factor == 0U || nchans % factor != 0U) {
        throw ValidationError(
            "psrio: fscrunch factor must be positive and divide nchans");
    }
    if (data.size() / std::max<std::uint64_t>(nchans, 1U) < nrows) {
        throw ValidationError(
            "psrio: fscrunch buffer is smaller than its rows");
    }
    const auto nout = static_cast<std::size_t>(nchans / factor);
    const auto rows = static_cast<std::size_t>(nrows);
    if (factor == 1U) {
        return data.first(rows * nout);
    }
    const auto width = static_cast<std::size_t>(nchans);
    const auto group = static_cast<std::size_t>(factor);
    for (std::size_t row = 0; row < rows; ++row) {
        const float* in = data.data() + (row * width);
        float* out      = data.data() + (row * nout);
        for (std::size_t ochan = 0; ochan < nout; ++ochan) {
            const float* src = in + (ochan * group);
            float acc        = src[0];
            for (std::size_t k = 1; k < group; ++k) {
                acc += src[k];
            }
            out[ochan] = acc;
        }
    }
    return data.first(rows * nout);
}

/// Streaming sum of consecutive groups of time samples.
///
/// Each push() sums whole groups of @p factor rows of @p row_len values and
/// keeps the remaining rows for the next push, so the output does not depend
/// on how the input is split into blocks. Rows are summed in time order in
/// float.
class TimeScruncher {
public:
    /// @throws ValidationError if @p row_len or @p factor is 0.
    TimeScruncher(std::uint64_t row_len, std::uint64_t factor)
        : m_row_len(static_cast<std::size_t>(row_len)),
          m_factor(static_cast<std::size_t>(factor)) {
        if (row_len == 0U || factor == 0U) {
            throw ValidationError(
                "psrio: time scrunch row length and factor must be positive");
        }
        m_carry.reserve(m_row_len * (m_factor - 1U));
    }

    [[nodiscard]] std::uint64_t factor() const noexcept { return m_factor; }

    /// Rows held over for the next push.
    [[nodiscard]] std::uint64_t pending() const noexcept {
        return m_carry.size() / m_row_len;
    }

    /// Add @p nrows rows from @p data and return the completed output rows.
    /// @p data is not modified. The returned span stays valid until the next
    /// push(); with factor 1 it is the front of @p data itself.
    /// @throws ValidationError if @p data holds fewer than nrows rows.
    std::span<float> push(std::span<float> data, std::uint64_t nrows) {
        if (data.size() / m_row_len < nrows) {
            throw ValidationError(
                "psrio: time scrunch input is smaller than its rows");
        }
        auto rows = static_cast<std::size_t>(nrows);
        data      = data.first(rows * m_row_len);
        if (m_factor == 1U) {
            return data;
        }
        const auto carried = m_carry.size() / m_row_len;
        const auto nout    = (carried + rows) / m_factor;
        m_out.resize(nout * m_row_len);
        std::size_t produced = 0;

        if (carried > 0U) {
            const auto need = m_factor - carried;
            if (rows < need) {
                m_carry.insert(m_carry.end(), data.begin(), data.end());
                return std::span<float>(m_out).first(0);
            }
            const std::span<float> head(m_out.data(), m_row_len);
            std::ranges::copy(std::span<const float>(m_carry).first(m_row_len),
                              head.begin());
            for (std::size_t r = 1; r < carried; ++r) {
                add_row(head, std::span<const float>(m_carry).subspan(
                                  r * m_row_len, m_row_len));
            }
            for (std::size_t r = 0; r < need; ++r) {
                add_row(head, data.subspan(r * m_row_len, m_row_len));
            }
            m_carry.clear();
            data = data.subspan(need * m_row_len);
            rows -= need;
            produced = 1;
        }

        for (; produced < nout; ++produced) {
            const std::span<float> out(m_out.data() + (produced * m_row_len),
                                       m_row_len);
            std::ranges::copy(data.first(m_row_len), out.begin());
            for (std::size_t r = 1; r < m_factor; ++r) {
                add_row(out, data.subspan(r * m_row_len, m_row_len));
            }
            data = data.subspan(m_factor * m_row_len);
            rows -= m_factor;
        }
        m_carry.assign(data.begin(), data.end());
        return std::span<float>(m_out).first(nout * m_row_len);
    }

    /// Drop the rows held over.
    void reset() noexcept { m_carry.clear(); }

private:
    static void add_row(std::span<float> acc, std::span<const float> row) {
        for (std::size_t i = 0; i < acc.size(); ++i) {
            acc[i] += row[i];
        }
    }

    std::size_t m_row_len;
    std::size_t m_factor;
    std::vector<float> m_carry;
    std::vector<float> m_out;
};

} // namespace psrio
