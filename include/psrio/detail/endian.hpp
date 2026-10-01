#pragma once

#include <bit>
#include <cstdint>
#include <cstring>
#include <type_traits>

namespace psrio::detail {

namespace endian_detail {

constexpr unsigned kBitsPerByte = 8U;

constexpr std::uint64_t kUint64HalfMask         = 0x00000000FFFFFFFFULL;
constexpr std::uint64_t kUint64HighHalfMask     = 0xFFFFFFFF00000000ULL;
constexpr std::uint64_t kUint64WordMask         = 0x0000FFFF0000FFFFULL;
constexpr std::uint64_t kUint64HighWordMask     = 0xFFFF0000FFFF0000ULL;
constexpr std::uint64_t kUint64ByteLaneMask     = 0x00FF00FF00FF00FFULL;
constexpr std::uint64_t kUint64HighByteLaneMask = 0xFF00FF00FF00FF00ULL;

} // namespace endian_detail

constexpr std::uint16_t swap_endian(std::uint16_t value) noexcept {
    const auto wide = static_cast<std::uint32_t>(value);
    return static_cast<std::uint16_t>((wide << endian_detail::kBitsPerByte) |
                                      (wide >> endian_detail::kBitsPerByte));
}

constexpr std::uint32_t swap_endian(std::uint32_t value) noexcept {
    constexpr std::uint32_t kByte0 = 0x000000FFU;
    constexpr std::uint32_t kByte1 = 0x0000FF00U;
    constexpr std::uint32_t kByte2 = 0x00FF0000U;
    constexpr std::uint32_t kByte3 = 0xFF000000U;
    constexpr unsigned kShift24    = endian_detail::kBitsPerByte * 3U;
    return ((value & kByte0) << kShift24) |
           ((value & kByte1) << endian_detail::kBitsPerByte) |
           ((value & kByte2) >> endian_detail::kBitsPerByte) |
           ((value & kByte3) >> kShift24);
}

constexpr std::uint64_t swap_endian(std::uint64_t value) noexcept {
    value = ((value & endian_detail::kUint64HalfMask)
             << (endian_detail::kBitsPerByte * 4U)) |
            ((value & endian_detail::kUint64HighHalfMask) >>
             (endian_detail::kBitsPerByte * 4U));
    value = ((value & endian_detail::kUint64WordMask)
             << (endian_detail::kBitsPerByte * 2U)) |
            ((value & endian_detail::kUint64HighWordMask) >>
             (endian_detail::kBitsPerByte * 2U));
    return ((value & endian_detail::kUint64ByteLaneMask)
            << endian_detail::kBitsPerByte) |
           ((value & endian_detail::kUint64HighByteLaneMask) >>
            endian_detail::kBitsPerByte);
}

/// Swap integer endianness (identity on single-byte types).
template <typename T>
    requires std::is_integral_v<T> && (sizeof(T) == 1 || sizeof(T) == 2 ||
                                       sizeof(T) == 4 || sizeof(T) == 8)
constexpr T swap_endian(T value) noexcept {
    if constexpr (sizeof(T) == 1) {
        return value;
    } else if constexpr (sizeof(T) == 2) {
        return static_cast<T>(swap_endian(static_cast<std::uint16_t>(value)));
    } else if constexpr (sizeof(T) == 4) {
        return static_cast<T>(swap_endian(static_cast<std::uint32_t>(value)));
    } else {
        return static_cast<T>(swap_endian(static_cast<std::uint64_t>(value)));
    }
}

/// @return true when the host uses little-endian byte order.
constexpr bool host_is_little_endian() noexcept {
    return std::endian::native == std::endian::little;
}

/// Interpret @p value as little-endian. Identity on a little-endian host.
template <typename T>
    requires std::is_integral_v<T> && (sizeof(T) == 1 || sizeof(T) == 2 ||
                                       sizeof(T) == 4 || sizeof(T) == 8)
constexpr T from_little_endian(T value) noexcept {
    if (host_is_little_endian()) {
        return value;
    }
    return swap_endian(value);
}

/// Interpret @p value as an IEEE-754 binary32 stored little-endian.
inline float from_little_endian(float value) noexcept {
    if (host_is_little_endian()) {
        return value;
    }
    const auto bits = swap_endian(std::bit_cast<std::uint32_t>(value));
    return std::bit_cast<float>(bits);
}

/// Interpret @p value as an IEEE-754 binary64 stored little-endian.
inline double from_little_endian(double value) noexcept {
    if (host_is_little_endian()) {
        return value;
    }
    const auto bits = swap_endian(std::bit_cast<std::uint64_t>(value));
    return std::bit_cast<double>(bits);
}

/// Convert @p value to little-endian representation.
template <typename T>
    requires std::is_integral_v<T> && (sizeof(T) == 1 || sizeof(T) == 2 ||
                                       sizeof(T) == 4 || sizeof(T) == 8)
constexpr T to_little_endian(T value) noexcept {
    return from_little_endian(value);
}

/// Convert IEEE-754 binary32 to little-endian representation.
inline float to_little_endian(float value) noexcept {
    return from_little_endian(value);
}

/// Convert IEEE-754 binary64 to little-endian representation.
inline double to_little_endian(double value) noexcept {
    return from_little_endian(value);
}

/// Load a little-endian value from an unaligned address.
template <typename T> inline T load_little_endian(const void* source) noexcept {
    T value{};
    std::memcpy(&value, source, sizeof(T));
    return from_little_endian(value);
}

} // namespace psrio::detail
