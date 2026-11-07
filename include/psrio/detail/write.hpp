#pragma once

#include "psrio/detail/endian.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <ostream>
#include <span>

namespace psrio::detail {

/// Write @p values to @p out as little-endian IEEE-754 binary32.
/// @return false if the stream failed.
inline bool write_float32_le(std::ostream& out, std::span<const float> values) {
    if (host_is_little_endian()) {
        const auto bytes = std::as_bytes(values);
        out.write(reinterpret_cast<const char*>(bytes.data()),
                  static_cast<std::streamsize>(bytes.size()));
        return static_cast<bool>(out);
    }
    constexpr std::size_t kBatch = 1024;
    std::array<float, kBatch> swapped{};
    for (std::size_t done = 0; done < values.size();) {
        const auto count = std::min(kBatch, values.size() - done);
        for (std::size_t i = 0; i < count; ++i) {
            swapped[i] = to_little_endian(values[done + i]);
        }
        out.write(reinterpret_cast<const char*>(swapped.data()),
                  static_cast<std::streamsize>(count * sizeof(float)));
        done += count;
    }
    return static_cast<bool>(out);
}

} // namespace psrio::detail
