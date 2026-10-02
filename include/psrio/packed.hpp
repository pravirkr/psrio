#pragma once

#include "psrio/common/types.hpp"
#include "psrio/detail/endian.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/packed_bits.hpp"

#include <bit>
#include <cstddef>
#include <cstdint>
#include <format>
#include <limits>
#include <span>

/**
 * @file packed.hpp
 * @brief Byte packing and unpacking utilities for packed and floating-point
 * filterbank buffers.
 */

namespace psrio {

// -----------------------------------------------------------------------------
// Bulk Unpacking Routines
// -----------------------------------------------------------------------------

/// Unpack sub-byte samples (1, 2, or 4 bits) into float, uint8_t, or uint16_t.
template <typename T>
inline void unpack_sub_byte(std::span<const std::byte> packed,
                            std::span<T> dest,
                            int nbits,
                            BitOrder order = BitOrder::kLsbFirst) {
    detail::unpack_sub_byte(packed, dest, nbits, order);
}

/// Unpack 8-bit integer samples into float, uint8_t, or uint16_t.
template <typename T>
inline void unpack_8bit(std::span<const std::byte> packed,
                        std::span<T> dest,
                        bool is_signed = false) {
    detail::unpack_8bit(packed, dest, is_signed);
}

/// Unpack 16-bit little-endian samples into float or uint16_t.
template <typename T>
inline void unpack_16le(std::span<const std::byte> packed, std::span<T> dest) {
    detail::unpack_16le(packed, dest);
}

/// Unpack 32-bit little-endian IEEE float samples.
inline void unpack_32le(std::span<const std::byte> packed,
                        std::span<float> dest) {
    detail::unpack_32le(packed, dest);
}

// -----------------------------------------------------------------------------
// Channel Inversion / Reversal
// -----------------------------------------------------------------------------

/// Reverse channel order inside each packed time sample, in place.
///
/// Supports 1, 2, 4, 8, 16, and 32-bit samples. For sub-byte rows, bytes across
/// the row are reversed and the channels within each byte are inverted.
/// An odd-length row also reverses the internal channels of its middle byte.
inline void reverse_channels(std::span<std::byte> gulp,
                             std::uint64_t nsamps,
                             std::uint64_t nchans,
                             int nbits) {
    detail::reverse_channels_impl(gulp, nsamps, nchans, nbits);
}

// -----------------------------------------------------------------------------
// On-the-Fly Packing & Extraction Utilities
// -----------------------------------------------------------------------------

/**
 * @brief Random-access sample extractor for packed or floating-point filterbank
 * buffers.
 *
 * Provides point extraction `extract_sample(t, c)` for row-major buffers
 * without unpacking entire blocks into temporary memory.
 */
class DataExtractor {
public:
    DataExtractor(std::span<const std::byte> data,
                  std::uint64_t nsamps,
                  std::uint64_t nchans,
                  SampleType type,
                  BitOrder bit_order = BitOrder::kLsbFirst)
        : m_data(data),
          m_nsamps(nsamps),
          m_nchans(nchans),
          m_type(type),
          m_bit_order(bit_order) {
        validate_buffer();
    }

    DataExtractor(std::span<const std::byte> data,
                  std::uint64_t nsamps,
                  std::uint64_t nchans,
                  int nbits,
                  bool samples_signed = false,
                  BitOrder bit_order  = BitOrder::kLsbFirst)
        : DataExtractor(data,
                        nsamps,
                        nchans,
                        sample_type_from_nbits(nbits, samples_signed),
                        bit_order) {}

    [[nodiscard]] std::uint64_t nsamps() const noexcept { return m_nsamps; }
    [[nodiscard]] std::uint64_t nchans() const noexcept { return m_nchans; }
    [[nodiscard]] SampleType sample_type() const noexcept { return m_type; }
    [[nodiscard]] BitOrder bit_order() const noexcept { return m_bit_order; }

    /// Extract sample at time index @p t, channel index @p c as float.
    [[nodiscard]] float extract_sample(std::uint64_t t, std::uint64_t c) const {
        if (t >= m_nsamps || c >= m_nchans) {
            throw ValidationError(std::format(
                "psrio: sample indices ({}, {}) out of range ({}, {})", t, c,
                m_nsamps, m_nchans));
        }

        const auto sample_index = (t * m_nchans) + c;

        if (m_type == SampleType::kFloat32) {
            const auto offset = sample_index * sizeof(float);
            const auto raw    = detail::load_little_endian<std::uint32_t>(
                m_data.data() + offset);
            return std::bit_cast<float>(raw);
        }

        if (m_type == SampleType::kUInt32) {
            const auto offset = sample_index * sizeof(std::uint32_t);
            return static_cast<float>(detail::load_little_endian<std::uint32_t>(
                m_data.data() + offset));
        }

        if (m_type == SampleType::kUInt16) {
            const auto offset = sample_index * sizeof(std::uint16_t);
            return static_cast<float>(detail::load_little_endian<std::uint16_t>(
                m_data.data() + offset));
        }

        if (m_type == SampleType::kUInt8) {
            return static_cast<float>(
                std::to_integer<std::uint8_t>(m_data[sample_index]));
        }

        if (m_type == SampleType::kInt8) {
            return static_cast<float>(static_cast<std::int8_t>(
                std::to_integer<std::uint8_t>(m_data[sample_index])));
        }

        // Sub-byte integers: 1, 2, or 4 bits.
        const int nbits             = bits_of(m_type);
        const auto samples_per_byte = 8U / static_cast<std::uint64_t>(nbits);
        const auto byte_offset      = sample_index / samples_per_byte;
        const auto sub_index        = sample_index % samples_per_byte;

        const auto byte_val =
            std::to_integer<std::uint8_t>(m_data[byte_offset]);
        const auto shift =
            (m_bit_order == BitOrder::kLsbFirst)
                ? (sub_index * static_cast<std::uint64_t>(nbits))
                : (8U - ((sub_index + 1U) * static_cast<std::uint64_t>(nbits)));
        const std::uint8_t mask = static_cast<std::uint8_t>((1U << nbits) - 1U);
        return static_cast<float>((byte_val >> shift) & mask);
    }

    /// Extract the complete spectrum at time index @p t into @p dest.
    void extract_spectrum(std::uint64_t t, std::span<float> dest) const {
        if (dest.size() < m_nchans) {
            throw ValidationError(std::format(
                "psrio: destination size {} is smaller than nchans {}",
                dest.size(), m_nchans));
        }
        for (std::uint64_t c = 0; c < m_nchans; ++c) {
            dest[c] = extract_sample(t, c);
        }
    }

private:
    void validate_buffer() const {
        if (m_nchans == 0U || m_nsamps == 0U) {
            return;
        }
        const auto width = static_cast<std::uint64_t>(bits_of(m_type));
        if (m_nchans > std::numeric_limits<std::uint64_t>::max() / width) {
            throw ValidationError("psrio: channel geometry is too large");
        }
        const auto total_bits   = m_nsamps * m_nchans * width;
        const auto bytes_needed = (total_bits + 7U) / 8U;
        if (m_data.size() < bytes_needed) {
            throw ValidationError(std::format(
                "psrio: DataExtractor buffer ({} bytes) too small for "
                "requested geometry (needs {} bytes)",
                m_data.size(), bytes_needed));
        }
    }

    std::span<const std::byte> m_data;
    std::uint64_t m_nsamps{0};
    std::uint64_t m_nchans{0};
    SampleType m_type{SampleType::kUInt8};
    BitOrder m_bit_order{BitOrder::kLsbFirst};
};

/**
 * @brief Write a single sample into a packed row-major buffer at time @p t,
 * channel @p c.
 */
inline void write_sample(std::span<std::byte> data,
                         std::uint64_t stride,
                         std::uint64_t samps,
                         std::uint64_t chans,
                         SampleType type,
                         float value,
                         BitOrder bit_order = BitOrder::kLsbFirst) {
    detail::write_sample_impl(data, stride, samps, chans, type, value,
                              bit_order);
}

} // namespace psrio
