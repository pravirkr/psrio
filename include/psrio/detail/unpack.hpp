#pragma once

#include <cstddef>
#include <cstdint>
#include <span>

namespace psrio::detail {

/// Unpack one 8-bit sample per byte into @p output (identity mapping to float).
inline void unpack_8bit_to_float(std::span<const std::uint8_t> input, std::span<float> output)
{
    if (input.size() != output.size()) {
        return;
    }
    for (std::size_t index = 0; index < input.size(); ++index) {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
        output[index] = static_cast<float>(input[index]);
    }
}

} // namespace psrio::detail
