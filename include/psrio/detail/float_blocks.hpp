#pragma once

#include "psrio/block_source.hpp"
#include "psrio/detail/exceptions.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <future>
#include <limits>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace psrio::detail {

/// Float blocks of a BlockSource read one block ahead on another thread, so
/// the caller's work on one block overlaps the read of the next.
///
/// Rows are `nifs * nchans` floats in the source's channel order.
class FloatBlocks {
public:
    /// Read samples [@p start, @p start + @p count) of @p source in blocks of
    /// @p gulp samples. Without @p count, read to the end of the source.
    /// @throws ValidationError if @p gulp is 0 or @p start is past the end.
    FloatBlocks(BlockSource& source,
                std::uint64_t start,
                std::optional<std::uint64_t> count,
                std::uint64_t gulp)
        : m_source(&source),
          m_row(source.nifs() * source.nchans()),
          m_gulp(gulp) {
        if (gulp == 0U) {
            throw ValidationError("psrio: block size must be positive");
        }
        if (source.has_nsamples()) {
            const auto total = source.nsamples();
            if (start > total) {
                throw ValidationError(
                    "psrio: start sample is past the end of the input");
            }
            m_left = std::min(count.value_or(total - start), total - start);
        } else {
            m_left = count.value_or(std::numeric_limits<std::uint64_t>::max());
        }
        source.seek(start);
        launch();
    }

    ~FloatBlocks() { wait(); }

    FloatBlocks(const FloatBlocks&)            = delete;
    FloatBlocks& operator=(const FloatBlocks&) = delete;
    FloatBlocks(FloatBlocks&&)                 = delete;
    FloatBlocks& operator=(FloatBlocks&&)      = delete;

    /// Floats in one sample.
    [[nodiscard]] std::uint64_t row() const noexcept { return m_row; }

    /// The next block, `rows * row()` floats; empty at the end. The span is
    /// valid until the next call.
    std::span<const float> next() {
        if (!m_pending.valid()) {
            return {};
        }
        const auto got = m_pending.get();
        std::swap(m_current, m_next);
        if (got == 0U) {
            m_left = 0;
            return {};
        }
        m_left -= std::min(m_left, got);
        launch();
        return std::span<const float>(m_current).first(
            static_cast<std::size_t>(got * m_row));
    }

private:
    void launch() {
        const auto want = std::min(m_gulp, m_left);
        if (want == 0U) {
            return;
        }
        m_next.resize(static_cast<std::size_t>(want * m_row));
        m_pending =
            std::async(std::launch::async, [source = m_source, want,
                                            block  = std::span<float>(m_next)] {
                return source->read_samples(want, block);
            });
    }

    void wait() noexcept {
        if (m_pending.valid()) {
            m_pending.wait();
        }
    }

    BlockSource* m_source;
    std::uint64_t m_row;
    std::uint64_t m_gulp;
    std::uint64_t m_left{0};
    std::vector<float> m_current;
    std::vector<float> m_next;
    // Last, so it is destroyed (and joined) before the buffers it fills.
    std::future<std::uint64_t> m_pending;
};

} // namespace psrio::detail
