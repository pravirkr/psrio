#pragma once

#include "psrio/common/types.hpp"
#include "psrio/detail/endian.hpp"
#include "psrio/detail/exceptions.hpp"

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <span>
#include <stdexcept>

namespace psrio::detail {

namespace unpack_detail {

inline void require_count(std::size_t actual, std::size_t expected) {
    if (actual != expected) {
        throw std::invalid_argument(
            "psrio: unpack output size does not match packed input");
    }
}

template <std::size_t NBits> struct LookupTable {
    static constexpr std::size_t kSize     = 256;
    static constexpr std::size_t kElements = 8 / NBits;
    alignas(64) std::array<std::array<std::uint8_t, kElements>, kSize> lsb{};
    alignas(64) std::array<std::array<std::uint8_t, kElements>, kSize> msb{};

    constexpr LookupTable() noexcept {
        const auto mask = static_cast<std::uint8_t>((1U << NBits) - 1U);
        for (std::size_t byte_val = 0; byte_val < kSize; ++byte_val) {
            for (std::size_t elem = 0; elem < kElements; ++elem) {
                const auto lsb_val = static_cast<std::uint8_t>(
                    (byte_val >> (elem * NBits)) & mask);
                lsb[byte_val][elem] = lsb_val;

                const auto msb_val = static_cast<std::uint8_t>(
                    (byte_val >> ((kElements - 1U - elem) * NBits)) & mask);
                msb[byte_val][elem] = msb_val;
            }
        }
    }
};

inline constexpr LookupTable<1> kLookup1Bit{};
inline constexpr LookupTable<2> kLookup2Bit{};
inline constexpr LookupTable<4> kLookup4Bit{};

inline void unpack_sub_byte_u8(std::span<const std::byte> packed,
                               std::span<std::uint8_t> output,
                               int nbits,
                               BitOrder order) {
    if (nbits == 1) {
        const auto& table =
            (order == BitOrder::kLsbFirst) ? kLookup1Bit.lsb : kLookup1Bit.msb;
        for (std::size_t index = 0; index < packed.size(); ++index) {
            const auto byte_val = std::to_integer<std::uint8_t>(packed[index]);
            std::memcpy(output.data() + (index * 8U), table[byte_val].data(),
                        8);
        }
        return;
    }
    if (nbits == 2) {
        const auto& table =
            (order == BitOrder::kLsbFirst) ? kLookup2Bit.lsb : kLookup2Bit.msb;
        for (std::size_t index = 0; index < packed.size(); ++index) {
            const auto byte_val = std::to_integer<std::uint8_t>(packed[index]);
            std::memcpy(output.data() + (index * 4U), table[byte_val].data(),
                        4);
        }
        return;
    }
    const auto& table =
        (order == BitOrder::kLsbFirst) ? kLookup4Bit.lsb : kLookup4Bit.msb;
    for (std::size_t index = 0; index < packed.size(); ++index) {
        const auto byte_val = std::to_integer<std::uint8_t>(packed[index]);
        std::memcpy(output.data() + (index * 2U), table[byte_val].data(), 2);
    }
}

inline void unpack_sub_byte_float(std::span<const std::byte> packed,
                                  std::span<float> output,
                                  int nbits,
                                  BitOrder order) {
    if (nbits == 1) {
        const auto& table =
            (order == BitOrder::kLsbFirst) ? kLookup1Bit.lsb : kLookup1Bit.msb;
        for (std::size_t index = 0; index < packed.size(); ++index) {
            const auto byte_val = std::to_integer<std::uint8_t>(packed[index]);
            const auto& entry   = table[byte_val];
            float* const dest   = output.data() + (index * 8U);
            for (std::size_t j = 0; j < 8; ++j) {
                dest[j] = static_cast<float>(entry[j]);
            }
        }
        return;
    }
    if (nbits == 2) {
        const auto& table =
            (order == BitOrder::kLsbFirst) ? kLookup2Bit.lsb : kLookup2Bit.msb;
        for (std::size_t index = 0; index < packed.size(); ++index) {
            const auto byte_val = std::to_integer<std::uint8_t>(packed[index]);
            const auto& entry   = table[byte_val];
            float* const dest   = output.data() + (index * 4U);
            for (std::size_t j = 0; j < 4; ++j) {
                dest[j] = static_cast<float>(entry[j]);
            }
        }
        return;
    }
    const auto& table =
        (order == BitOrder::kLsbFirst) ? kLookup4Bit.lsb : kLookup4Bit.msb;
    for (std::size_t index = 0; index < packed.size(); ++index) {
        const auto byte_val = std::to_integer<std::uint8_t>(packed[index]);
        const auto& entry   = table[byte_val];
        float* const dest   = output.data() + (index * 2U);
        dest[0]             = static_cast<float>(entry[0]);
        dest[1]             = static_cast<float>(entry[1]);
    }
}

} // namespace unpack_detail

/// Unpack one 8-bit sample per byte into @p output.
/// @throws std::invalid_argument if the spans differ in size.
inline void unpack_8bit_to_float(std::span<const std::uint8_t> input,
                                 std::span<float> output) {
    unpack_detail::require_count(output.size(), input.size());
    for (std::size_t index = 0; index < input.size(); ++index) {
        output[index] = static_cast<float>(input[index]);
    }
}

/// Unpack 1-, 2-, or 4-bit samples according to @p order.
///
/// Under kLsbFirst (DSPSR convention): within each byte, bit 0 is the first
/// 1-bit sample, bits 0–1 are the first 2-bit sample, and the low nibble is the
/// first 4-bit sample. Under kMsbFirst (PRESTO convention): within each byte,
/// bit 7 is the first 1-bit sample, bits 6–7 are the first 2-bit sample, and
/// the high nibble is the first 4-bit sample.
/// @throws std::invalid_argument if @p nbits is not 1, 2, or 4, or if
///         @p output is not exactly the unpacked length of @p packed.
template <typename T>
    requires std::same_as<T, float> || std::same_as<T, std::uint8_t>
inline void unpack_sub_byte(std::span<const std::byte> packed,
                            std::span<T> output,
                            int nbits,
                            BitOrder order = BitOrder::kLsbFirst) {
    if (nbits != 1 && nbits != 2 && nbits != 4) {
        throw std::invalid_argument(
            "psrio: sub-byte unpack nbits must be 1, 2, or 4");
    }
    const auto samples_per_byte = static_cast<std::size_t>(8 / nbits);
    if (packed.size() > (static_cast<std::size_t>(-1) / samples_per_byte)) {
        throw std::invalid_argument("psrio: unpack input is too large");
    }
    unpack_detail::require_count(output.size(),
                                 packed.size() * samples_per_byte);

    if constexpr (std::same_as<T, std::uint8_t>) {
        unpack_detail::unpack_sub_byte_u8(packed, output, nbits, order);
    } else {
        unpack_detail::unpack_sub_byte_float(packed, output, nbits, order);
    }
}

/// Unpack 1-, 2-, or 4-bit samples, least-significant field first (DSPSR
/// convention).
/// @throws std::invalid_argument if @p nbits is not 1, 2, or 4, or if
///         @p output is not exactly the unpacked length of @p packed.
template <typename T>
    requires std::same_as<T, float> || std::same_as<T, std::uint8_t>
inline void
unpack_lsb(std::span<const std::byte> packed, std::span<T> output, int nbits) {
    unpack_sub_byte(packed, output, nbits, BitOrder::kLsbFirst);
}

/// Pack unpacked 1-, 2-, or 4-bit samples into packed bytes according to @p
/// order.
/// @throws std::invalid_argument if @p nbits is not 1, 2, or 4, or lengths do
/// not match.
inline void pack_sub_byte(std::span<const std::uint8_t> unpacked,
                          std::span<std::byte> packed,
                          int nbits,
                          BitOrder order = BitOrder::kLsbFirst) {
    if (nbits != 1 && nbits != 2 && nbits != 4) {
        throw std::invalid_argument(
            "psrio: sub-byte pack nbits must be 1, 2, or 4");
    }
    const auto samples_per_byte = static_cast<std::size_t>(8 / nbits);
    unpack_detail::require_count(unpacked.size(),
                                 packed.size() * samples_per_byte);

    if (nbits == 1) {
        for (std::size_t index = 0; index < packed.size(); ++index) {
            const auto pos = index * 8U;
            std::uint8_t byte_val{0};
            if (order == BitOrder::kMsbFirst) {
                byte_val =
                    static_cast<std::uint8_t>(((unpacked[pos + 0] & 1U) << 7U) |
                                              ((unpacked[pos + 1] & 1U) << 6U) |
                                              ((unpacked[pos + 2] & 1U) << 5U) |
                                              ((unpacked[pos + 3] & 1U) << 4U) |
                                              ((unpacked[pos + 4] & 1U) << 3U) |
                                              ((unpacked[pos + 5] & 1U) << 2U) |
                                              ((unpacked[pos + 6] & 1U) << 1U) |
                                              (unpacked[pos + 7] & 1U));
            } else {
                byte_val =
                    static_cast<std::uint8_t>((unpacked[pos + 0] & 1U) |
                                              ((unpacked[pos + 1] & 1U) << 1U) |
                                              ((unpacked[pos + 2] & 1U) << 2U) |
                                              ((unpacked[pos + 3] & 1U) << 3U) |
                                              ((unpacked[pos + 4] & 1U) << 4U) |
                                              ((unpacked[pos + 5] & 1U) << 5U) |
                                              ((unpacked[pos + 6] & 1U) << 6U) |
                                              ((unpacked[pos + 7] & 1U) << 7U));
            }
            packed[index] = std::byte{byte_val};
        }
        return;
    }
    if (nbits == 2) {
        for (std::size_t index = 0; index < packed.size(); ++index) {
            const auto pos = index * 4U;
            std::uint8_t byte_val{0};
            if (order == BitOrder::kMsbFirst) {
                byte_val =
                    static_cast<std::uint8_t>(((unpacked[pos + 0] & 3U) << 6U) |
                                              ((unpacked[pos + 1] & 3U) << 4U) |
                                              ((unpacked[pos + 2] & 3U) << 2U) |
                                              (unpacked[pos + 3] & 3U));
            } else {
                byte_val =
                    static_cast<std::uint8_t>((unpacked[pos + 0] & 3U) |
                                              ((unpacked[pos + 1] & 3U) << 2U) |
                                              ((unpacked[pos + 2] & 3U) << 4U) |
                                              ((unpacked[pos + 3] & 3U) << 6U));
            }
            packed[index] = std::byte{byte_val};
        }
        return;
    }
    for (std::size_t index = 0; index < packed.size(); ++index) {
        const auto pos = index * 2U;
        std::uint8_t byte_val{0};
        if (order == BitOrder::kMsbFirst) {
            byte_val =
                static_cast<std::uint8_t>(((unpacked[pos + 0] & 0x0FU) << 4U) |
                                          (unpacked[pos + 1] & 0x0FU));
        } else {
            byte_val =
                static_cast<std::uint8_t>((unpacked[pos + 0] & 0x0FU) |
                                          ((unpacked[pos + 1] & 0x0FU) << 4U));
        }
        packed[index] = std::byte{byte_val};
    }
}

/// Pack unpacked 1-, 2-, or 4-bit samples in place within @p buffer.
/// Returns the subspan containing the packed bytes.
template <typename ByteT>
    requires std::same_as<ByteT, std::uint8_t> || std::same_as<ByteT, std::byte>
inline std::span<ByteT> pack_sub_byte_inplace(
    std::span<ByteT> buffer, int nbits, BitOrder order = BitOrder::kLsbFirst) {
    if (nbits != 1 && nbits != 2 && nbits != 4) {
        throw ValidationError("psrio: sub-byte pack nbits must be 1, 2, or 4");
    }
    const auto samples_per_byte = static_cast<std::size_t>(8 / nbits);
    if (buffer.size() % samples_per_byte != 0U) {
        throw ValidationError("psrio: buffer size must be a multiple of "
                              "samples per byte for in-place pack");
    }
    const auto packed_bytes = buffer.size() / samples_per_byte;
    auto* ptr               = reinterpret_cast<std::uint8_t*>(buffer.data());

    if (nbits == 1) {
        for (std::size_t index = 0; index < packed_bytes; ++index) {
            const auto pos = index * 8U;
            std::uint8_t byte_val{0};
            if (order == BitOrder::kMsbFirst) {
                byte_val = static_cast<std::uint8_t>(
                    ((ptr[pos + 0] & 1U) << 7U) | ((ptr[pos + 1] & 1U) << 6U) |
                    ((ptr[pos + 2] & 1U) << 5U) | ((ptr[pos + 3] & 1U) << 4U) |
                    ((ptr[pos + 4] & 1U) << 3U) | ((ptr[pos + 5] & 1U) << 2U) |
                    ((ptr[pos + 6] & 1U) << 1U) | (ptr[pos + 7] & 1U));
            } else {
                byte_val = static_cast<std::uint8_t>(
                    (ptr[pos + 0] & 1U) | ((ptr[pos + 1] & 1U) << 1U) |
                    ((ptr[pos + 2] & 1U) << 2U) | ((ptr[pos + 3] & 1U) << 3U) |
                    ((ptr[pos + 4] & 1U) << 4U) | ((ptr[pos + 5] & 1U) << 5U) |
                    ((ptr[pos + 6] & 1U) << 6U) | ((ptr[pos + 7] & 1U) << 7U));
            }
            ptr[index] = byte_val;
        }
    } else if (nbits == 2) {
        for (std::size_t index = 0; index < packed_bytes; ++index) {
            const auto pos = index * 4U;
            std::uint8_t byte_val{0};
            if (order == BitOrder::kMsbFirst) {
                byte_val = static_cast<std::uint8_t>(
                    ((ptr[pos + 0] & 3U) << 6U) | ((ptr[pos + 1] & 3U) << 4U) |
                    ((ptr[pos + 2] & 3U) << 2U) | (ptr[pos + 3] & 3U));
            } else {
                byte_val = static_cast<std::uint8_t>(
                    (ptr[pos + 0] & 3U) | ((ptr[pos + 1] & 3U) << 2U) |
                    ((ptr[pos + 2] & 3U) << 4U) | ((ptr[pos + 3] & 3U) << 6U));
            }
            ptr[index] = byte_val;
        }
    } else {
        for (std::size_t index = 0; index < packed_bytes; ++index) {
            const auto pos = index * 2U;
            std::uint8_t byte_val{0};
            if (order == BitOrder::kMsbFirst) {
                byte_val = static_cast<std::uint8_t>(
                    ((ptr[pos + 0] & 0x0FU) << 4U) | (ptr[pos + 1] & 0x0FU));
            } else {
                byte_val = static_cast<std::uint8_t>(
                    (ptr[pos + 0] & 0x0FU) | ((ptr[pos + 1] & 0x0FU) << 4U));
            }
            ptr[index] = byte_val;
        }
    }
    return buffer.subspan(0, packed_bytes);
}

/// Unpack 1-, 2-, or 4-bit samples in place backwards within @p buffer.
template <typename ByteT>
    requires std::same_as<ByteT, std::uint8_t> || std::same_as<ByteT, std::byte>
inline void unpack_sub_byte_inplace(std::span<ByteT> buffer,
                                    int nbits,
                                    BitOrder order = BitOrder::kLsbFirst) {
    if (nbits != 1 && nbits != 2 && nbits != 4) {
        throw ValidationError(
            "psrio: sub-byte unpack nbits must be 1, 2, or 4");
    }
    const auto samples_per_byte = static_cast<std::size_t>(8 / nbits);
    if (buffer.size() % samples_per_byte != 0U) {
        throw ValidationError("psrio: buffer size must be a multiple of "
                              "samples per byte for in-place unpack");
    }
    const auto packed_bytes = buffer.size() / samples_per_byte;
    auto* ptr               = reinterpret_cast<std::uint8_t*>(buffer.data());

    if (nbits == 1) {
        const auto& table = (order == BitOrder::kLsbFirst)
                                ? unpack_detail::kLookup1Bit.lsb
                                : unpack_detail::kLookup1Bit.msb;
        for (std::size_t count = packed_bytes; count > 0; --count) {
            const auto index    = count - 1U;
            const auto byte_val = ptr[index];
            const auto pos      = index * 8U;
            std::memcpy(ptr + pos, table[byte_val].data(), 8);
        }
    } else if (nbits == 2) {
        const auto& table = (order == BitOrder::kLsbFirst)
                                ? unpack_detail::kLookup2Bit.lsb
                                : unpack_detail::kLookup2Bit.msb;
        for (std::size_t count = packed_bytes; count > 0; --count) {
            const auto index    = count - 1U;
            const auto byte_val = ptr[index];
            const auto pos      = index * 4U;
            std::memcpy(ptr + pos, table[byte_val].data(), 4);
        }
    } else {
        const auto& table = (order == BitOrder::kLsbFirst)
                                ? unpack_detail::kLookup4Bit.lsb
                                : unpack_detail::kLookup4Bit.msb;
        for (std::size_t count = packed_bytes; count > 0; --count) {
            const auto index    = count - 1U;
            const auto byte_val = ptr[index];
            const auto pos      = index * 2U;
            std::memcpy(ptr + pos, table[byte_val].data(), 2);
        }
    }
}

/// Copy 8-bit samples. When @p signed_samples is true and @p T is float,
/// each byte is interpreted as int8. uint8 output keeps the raw byte.
/// @throws std::invalid_argument if the spans differ in size.
template <typename T>
    requires std::same_as<T, float> || std::same_as<T, std::uint8_t>
inline void unpack_8bit(std::span<const std::byte> packed,
                        std::span<T> output,
                        bool signed_samples) {
    unpack_detail::require_count(output.size(), packed.size());
    if constexpr (std::same_as<T, float>) {
        if (signed_samples) {
            for (std::size_t index = 0; index < packed.size(); ++index) {
                const auto raw = std::to_integer<unsigned char>(packed[index]);
                output[index] =
                    static_cast<float>(static_cast<std::int8_t>(raw));
            }
            return;
        }
        for (std::size_t index = 0; index < packed.size(); ++index) {
            output[index] = static_cast<float>(
                std::to_integer<unsigned char>(packed[index]));
        }
        return;
    }

    (void)signed_samples;
    for (std::size_t index = 0; index < packed.size(); ++index) {
        output[index] = std::to_integer<unsigned char>(packed[index]);
    }
}

/// Unpack little-endian uint16 samples into @p T (uint16 or float).
/// @throws std::invalid_argument if @p packed is not exactly two bytes per
/// output.
template <typename T>
    requires std::same_as<T, float> || std::same_as<T, std::uint16_t>
inline void unpack_16le(std::span<const std::byte> packed,
                        std::span<T> output) {
    if (packed.size() % 2U != 0U) {
        throw std::invalid_argument("psrio: 16-bit payload is not aligned");
    }
    unpack_detail::require_count(output.size(), packed.size() / 2U);
    for (std::size_t index = 0; index < output.size(); ++index) {
        const auto value =
            load_little_endian<std::uint16_t>(packed.data() + (index * 2U));
        output[index] = static_cast<T>(value);
    }
}

/// Unpack little-endian IEEE-754 binary32 samples.
/// @throws std::invalid_argument if @p packed is not exactly four bytes per
/// output.
inline void unpack_32le(std::span<const std::byte> packed,
                        std::span<float> output) {
    if (packed.size() % 4U != 0U) {
        throw std::invalid_argument("psrio: 32-bit payload is not aligned");
    }
    unpack_detail::require_count(output.size(), packed.size() / 4U);
    for (std::size_t index = 0; index < output.size(); ++index) {
        output[index] = load_little_endian<float>(packed.data() + (index * 4U));
    }
}

namespace packed_detail {

/// Reverse 8 individual bits inside a single byte.
[[nodiscard]] constexpr std::uint8_t reverse_bits(std::uint8_t value) noexcept {
    value = static_cast<std::uint8_t>(((value & 0xF0U) >> 4U) |
                                      ((value & 0x0FU) << 4U));
    value = static_cast<std::uint8_t>(((value & 0xCCU) >> 2U) |
                                      ((value & 0x33U) << 2U));
    value = static_cast<std::uint8_t>(((value & 0xAAU) >> 1U) |
                                      ((value & 0x55U) << 1U));
    return value;
}

/// In-byte permutation for 2-bit rows (4 channels/byte). It reverses the four
/// 2-bit fields: ch0<->ch3, ch1<->ch2.
[[nodiscard]] constexpr std::uint8_t
reverse_bit_pairs(std::uint8_t value) noexcept {
    return static_cast<std::uint8_t>(
        ((value & 0x03U) << 6U) | ((value & 0x0CU) << 2U) |
        ((value & 0x30U) >> 2U) | ((value & 0xC0U) >> 6U));
}

/// In-byte permutation for 4-bit rows (2 channels/byte). It swaps the two 4-bit
/// nibbles: ch0<->ch1.
[[nodiscard]] constexpr std::uint8_t swap_nibbles(std::uint8_t value) noexcept {
    return static_cast<std::uint8_t>(((value & 0x0FU) << 4U) |
                                     ((value & 0xF0U) >> 4U));
}

} // namespace packed_detail

inline void reverse_channels_impl(std::span<std::byte> gulp,
                                  std::uint64_t nsamps,
                                  std::uint64_t nchans,
                                  int nbits) {
    if (nbits != 1 && nbits != 2 && nbits != 4 && nbits != 8 && nbits != 16 &&
        nbits != 32) {
        throw ValidationError("psrio: channel reversal supports 1, 2, 4, 8, "
                              "16, and 32-bit samples");
    }
    if (nchans == 0U || nsamps == 0U) {
        throw ValidationError(
            "psrio: channel reversal needs at least one sample and channel");
    }
    const auto width = static_cast<std::uint64_t>(nbits);
    if (nchans > std::numeric_limits<std::uint64_t>::max() / width) {
        throw ValidationError("psrio: channel reversal block is too large");
    }
    if ((nchans * width) % 8U != 0U) {
        throw ValidationError(
            "psrio: nchans * nbits must be a whole number of bytes");
    }
    const auto row_bytes = (nchans * width) / 8U;
    if (row_bytes > std::numeric_limits<std::uint64_t>::max() / nsamps) {
        throw ValidationError("psrio: channel reversal block is too large");
    }
    const auto need = nsamps * row_bytes;
    if (gulp.size() < need) {
        throw ValidationError(
            "psrio: channel reversal buffer is smaller than the block");
    }

    auto* bytes = reinterpret_cast<std::uint8_t*>(gulp.data());

    if (nbits == 16) {
        for (std::uint64_t sample = 0; sample < nsamps; ++sample) {
            auto* row =
                reinterpret_cast<std::uint16_t*>(bytes + (sample * row_bytes));
            std::reverse(row, row + nchans);
        }
        return;
    }

    if (nbits == 32) {
        for (std::uint64_t sample = 0; sample < nsamps; ++sample) {
            auto* row =
                reinterpret_cast<std::uint32_t*>(bytes + (sample * row_bytes));
            std::reverse(row, row + nchans);
        }
        return;
    }

    for (std::uint64_t sample = 0; sample < nsamps; ++sample) {
        std::uint8_t* row = bytes + (sample * row_bytes);
        for (std::uint64_t index = 0; index < row_bytes / 2U; ++index) {
            const auto other        = row_bytes - (1U + index);
            std::uint8_t this_byte  = row[index];
            std::uint8_t other_byte = row[other];
            if (nbits == 1) {
                this_byte  = packed_detail::reverse_bits(this_byte);
                other_byte = packed_detail::reverse_bits(other_byte);
            } else if (nbits == 2) {
                this_byte  = packed_detail::reverse_bit_pairs(this_byte);
                other_byte = packed_detail::reverse_bit_pairs(other_byte);
            } else if (nbits == 4) {
                this_byte  = packed_detail::swap_nibbles(this_byte);
                other_byte = packed_detail::swap_nibbles(other_byte);
            }
            row[other] = this_byte;
            row[index] = other_byte;
        }
        // Handle middle byte on odd-length sub-byte rows.
        if (row_bytes % 2U != 0U) {
            const auto mid = row_bytes / 2U;
            if (nbits == 1) {
                row[mid] = packed_detail::reverse_bits(row[mid]);
            } else if (nbits == 2) {
                row[mid] = packed_detail::reverse_bit_pairs(row[mid]);
            } else if (nbits == 4) {
                row[mid] = packed_detail::swap_nibbles(row[mid]);
            }
        }
    }
}

inline void write_sample_impl(std::span<std::byte> data,
                              std::uint64_t stride,
                              std::uint64_t samps,
                              std::uint64_t chans,
                              SampleType type,
                              float value,
                              BitOrder bit_order) {
    const auto sample_index = (samps * stride) + chans;

    if (type == SampleType::kFloat32) {
        const auto byte_offset = sample_index * sizeof(float);
        if (byte_offset + sizeof(float) > data.size()) {
            throw ValidationError("psrio: write_sample offset out of range");
        }
        std::memcpy(data.data() + byte_offset, &value, sizeof(float));
        return;
    }

    if (type == SampleType::kUInt32) {
        const auto byte_offset = sample_index * sizeof(std::uint32_t);
        if (byte_offset + sizeof(std::uint32_t) > data.size()) {
            throw ValidationError("psrio: write_sample offset out of range");
        }
        const auto raw = static_cast<std::uint32_t>(std::max(0.0F, value));
        std::memcpy(data.data() + byte_offset, &raw, sizeof(std::uint32_t));
        return;
    }

    if (type == SampleType::kUInt16) {
        const auto byte_offset = sample_index * sizeof(std::uint16_t);
        if (byte_offset + sizeof(std::uint16_t) > data.size()) {
            throw ValidationError("psrio: write_sample offset out of range");
        }
        const auto raw = static_cast<std::uint16_t>(std::clamp(
            value, 0.0F,
            static_cast<float>(std::numeric_limits<std::uint16_t>::max())));
        std::memcpy(data.data() + byte_offset, &raw, sizeof(std::uint16_t));
        return;
    }

    if (type == SampleType::kInt16) {
        const auto byte_offset = sample_index * sizeof(std::int16_t);
        if (byte_offset + sizeof(std::int16_t) > data.size()) {
            throw ValidationError("psrio: write_sample offset out of range");
        }
        const auto raw    = static_cast<std::int16_t>(std::clamp(
            value, static_cast<float>(std::numeric_limits<std::int16_t>::min()),
            static_cast<float>(std::numeric_limits<std::int16_t>::max())));
        const auto stored = detail::to_little_endian(raw);
        std::memcpy(data.data() + byte_offset, &stored, sizeof(stored));
        return;
    }

    if (type == SampleType::kUInt8) {
        if (sample_index >= data.size()) {
            throw ValidationError("psrio: write_sample offset out of range");
        }
        const auto raw     = static_cast<std::uint8_t>(std::clamp(
            value, 0.0F,
            static_cast<float>(std::numeric_limits<std::uint8_t>::max())));
        data[sample_index] = static_cast<std::byte>(raw);
        return;
    }

    if (type == SampleType::kInt8) {
        if (sample_index >= data.size()) {
            throw ValidationError("psrio: write_sample offset out of range");
        }
        const auto raw = static_cast<std::int8_t>(std::clamp(
            value, static_cast<float>(std::numeric_limits<std::int8_t>::min()),
            static_cast<float>(std::numeric_limits<std::int8_t>::max())));
        data[sample_index] =
            static_cast<std::byte>(static_cast<std::uint8_t>(raw));
        return;
    }

    // Sub-byte integers: 1, 2, or 4 bits
    const int nbits             = bits_of(type);
    const auto samples_per_byte = 8U / static_cast<std::uint64_t>(nbits);
    const auto byte_offset      = sample_index / samples_per_byte;
    const auto sub_index        = sample_index % samples_per_byte;

    if (byte_offset >= data.size()) {
        throw ValidationError("psrio: write_sample offset out of range");
    }

    const auto shift =
        (bit_order == BitOrder::kLsbFirst)
            ? (sub_index * static_cast<std::uint64_t>(nbits))
            : (8U - ((sub_index + 1U) * static_cast<std::uint64_t>(nbits)));
    const auto mask = static_cast<std::uint8_t>((1U << nbits) - 1U);
    const auto raw  = static_cast<std::uint8_t>(std::max(0.0F, value)) & mask;

    auto current      = std::to_integer<std::uint8_t>(data[byte_offset]);
    current           = static_cast<std::uint8_t>(current & ~(mask << shift));
    current           = static_cast<std::uint8_t>(current | (raw << shift));
    data[byte_offset] = static_cast<std::byte>(current);
}

} // namespace psrio::detail
