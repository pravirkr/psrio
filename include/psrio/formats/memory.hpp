#pragma once

#include "psrio/astro.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/packed_bits.hpp"
#include "psrio/detail/skip.hpp"
#include "psrio/packed.hpp"

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <format>
#include <limits>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace psrio {

/// Metadata for an in-memory intensity block.
struct MemoryInfo {
    std::string source{"Unknown"};
    std::string telescope{"Unknown"};
    double raj{0.0};
    double dej{0.0};
    std::uint64_t nchans{0};
    std::uint64_t nifs{1};
    int nbits{8};
    bool samples_signed{false};
    BitOrder bit_order{BitOrder::kLsbFirst};
    double tsamp{0.0};
    double tstart{0.0};
    double fch1{0.0};
    double foff{0.0};
    int beam{0};
};

/// Packed intensity samples held in memory (owned or non-owning view).
///
/// Supports both owning `std::vector<std::byte>` storage and non-owning
/// `std::span<const std::byte>` zero-copy viewing.
class MemoryBlock {
public:
    /// Construct with owned byte buffer.
    MemoryBlock(MemoryInfo info, std::vector<std::byte> samples)
        : m_info(std::move(info)),
          m_owned(std::move(samples)) {
        m_bytes = m_owned;
        init();
    }

    /// Construct with non-owning span (zero copy).
    template <typename Span>
        requires(!std::same_as<std::remove_cvref_t<Span>,
                               std::vector<std::byte>>) &&
                    std::convertible_to<Span, std::span<const std::byte>>
    MemoryBlock(MemoryInfo info, Span samples)
        : m_info(std::move(info)),
          m_bytes(samples) {
        init();
    }

    [[nodiscard]] const MemoryInfo& info() const noexcept { return m_info; }

    [[nodiscard]] std::string_view source_name() const noexcept {
        return m_info.source;
    }

    [[nodiscard]] std::string_view telescope() const noexcept {
        return m_info.telescope;
    }

    [[nodiscard]] double raj() const noexcept { return m_info.raj; }
    [[nodiscard]] double dej() const noexcept { return m_info.dej; }

    [[nodiscard]] std::uint64_t nchans() const noexcept {
        return m_info.nchans;
    }
    [[nodiscard]] std::uint64_t nifs() const noexcept { return m_info.nifs; }
    [[nodiscard]] int nbits() const noexcept { return m_info.nbits; }
    [[nodiscard]] SampleType sample_type() const noexcept { return m_type; }
    [[nodiscard]] std::uint64_t bytes_per_sample() const noexcept {
        return m_stride;
    }
    [[nodiscard]] std::uint64_t nsamples() const noexcept { return m_nsamples; }
    static bool has_nsamples() noexcept { return true; }
    [[nodiscard]] double tsamp() const noexcept { return m_info.tsamp; }
    [[nodiscard]] double tstart() const noexcept { return m_info.tstart; }
    [[nodiscard]] double fch1() const noexcept { return m_info.fch1; }
    [[nodiscard]] double foff() const noexcept { return m_info.foff; }
    [[nodiscard]] int beam() const noexcept { return m_info.beam; }
    [[nodiscard]] BitOrder bit_order() const noexcept {
        return m_info.bit_order;
    }

    [[nodiscard]] double bandwidth() const noexcept {
        return std::abs(m_info.foff) * static_cast<double>(m_info.nchans);
    }

    [[nodiscard]] double center_frequency() const noexcept {
        return m_info.fch1 +
               (m_info.foff * (static_cast<double>(m_info.nchans) - 1.0) / 2.0);
    }

    [[nodiscard]] double spectra_rate() const noexcept {
        return m_info.tsamp > 0.0 ? 1.0 / m_info.tsamp : 0.0;
    }

    [[nodiscard]] std::time_t utc_start() const noexcept {
        return astro::mjd_to_time(m_info.tstart);
    }

    [[nodiscard]] std::uint64_t tell() const noexcept { return m_sample; }

    void rewind() noexcept { m_sample = 0; }

    void seek(std::uint64_t sample) {
        if (sample > m_nsamples) {
            throw ValidationError("psrio: seek is past the readable samples");
        }
        m_sample = sample;
    }

    void skip(std::int64_t delta) {
        m_sample = detail::apply_skip(m_sample, delta, m_nsamples);
    }

    void set_fswap(bool enable) noexcept { m_apply_fswap = enable; }
    [[nodiscard]] bool fswap_enabled() const noexcept { return m_apply_fswap; }

    /// Copy up to @p count packed time samples into @p dest.
    std::uint64_t read_block(std::uint64_t count, std::span<std::byte> dest) {
        if (count > 0U &&
            m_stride > (std::numeric_limits<std::uint64_t>::max() / count)) {
            throw ValidationError("psrio: requested sample block is too large");
        }
        const auto bytes_needed = count * m_stride;
        if (dest.size() < bytes_needed) {
            throw ValidationError(std::format(
                "psrio: destination has {} bytes but the block needs {}",
                dest.size(), bytes_needed));
        }
        if (m_sample > m_nsamples) {
            throw ValidationError(
                "psrio: reader cursor is past the readable samples");
        }
        const auto to_read = std::min(count, m_nsamples - m_sample);
        if (to_read == 0U) {
            return 0;
        }
        const auto nbytes = to_read * m_stride;
        std::memcpy(dest.data(), m_bytes.data() + (m_sample * m_stride),
                    static_cast<std::size_t>(nbytes));
        if (m_apply_fswap && to_read > 0U) {
            reverse_channels(dest.first(static_cast<std::size_t>(nbytes)),
                             to_read, m_info.nchans, m_info.nbits);
        }
        m_sample += to_read;
        return to_read;
    }

    /// Zero-copy view of the next @p count time samples.
    [[nodiscard]] std::span<const std::byte> view_block(std::uint64_t count) {
        if (count > 0U &&
            m_stride > (std::numeric_limits<std::uint64_t>::max() / count)) {
            throw ValidationError("psrio: requested sample block is too large");
        }
        const auto bytes_needed = count * m_stride;
        if (m_sample > m_nsamples || count > m_nsamples - m_sample) {
            throw ValidationError(
                "psrio: view_block extends past readable samples");
        }
        const auto view =
            m_bytes.subspan(static_cast<std::size_t>(m_sample * m_stride),
                            static_cast<std::size_t>(bytes_needed));
        m_sample += count;
        return view;
    }

    /// Zero-copy view of the next @p nbytes bytes.
    [[nodiscard]] std::span<const std::byte> view_bytes(std::uint64_t nbytes) {
        if (nbytes == 0U || m_stride == 0U || nbytes % m_stride != 0U) {
            throw ValidationError("psrio: byte request must be a positive "
                                  "multiple of the sample stride");
        }
        const auto count = nbytes / m_stride;
        return view_block(count);
    }

    std::uint64_t read_bytes(std::uint64_t nbytes, std::span<std::byte> dest) {
        if (dest.size() != nbytes) {
            throw ValidationError(
                "psrio: destination size does not match requested bytes");
        }
        if (nbytes == 0U || m_stride == 0U || nbytes % m_stride != 0U) {
            throw ValidationError("psrio: byte request must be a positive "
                                  "multiple of the sample stride");
        }
        const auto samples = nbytes / m_stride;
        if (m_sample > m_nsamples || samples > m_nsamples - m_sample) {
            throw ValidationError(
                "psrio: byte request extends past readable samples");
        }
        const auto copied = read_block(samples, dest);
        if (copied != samples) {
            throw ValidationError(
                "psrio: byte request extends past readable samples");
        }
        return nbytes;
    }

    template <typename T>
        requires std::same_as<T, float> || std::same_as<T, std::uint8_t> ||
                 std::same_as<T, std::uint16_t>
    std::uint64_t read_samples(std::uint64_t count, std::span<T> dest) {
        if constexpr (std::same_as<T, std::uint8_t>) {
            if (m_info.nbits > 8) {
                throw ValidationError(
                    "psrio: uint8_t output requires nbits of 8 or less");
            }
        } else if constexpr (std::same_as<T, std::uint16_t>) {
            if (m_info.nbits != 16) {
                throw ValidationError(
                    "psrio: uint16_t output requires 16-bit samples");
            }
        }

        const auto per_sample = m_info.nifs * m_info.nchans;
        if (count > 0U &&
            per_sample > (std::numeric_limits<std::uint64_t>::max() / count)) {
            throw ValidationError("psrio: requested sample block is too large");
        }
        const auto expected = count * per_sample;
        if (expected != dest.size()) {
            throw ValidationError(std::format(
                "psrio: destination has {} values but the request needs {}",
                dest.size(), expected));
        }
        if (m_sample > m_nsamples) {
            throw ValidationError(
                "psrio: reader cursor is past the readable samples");
        }
        const auto to_read = std::min(count, m_nsamples - m_sample);
        if (to_read == 0U) {
            return 0;
        }

        const auto nbytes = to_read * m_stride;
        const auto raw =
            m_bytes.subspan(static_cast<std::size_t>(m_sample * m_stride),
                            static_cast<std::size_t>(nbytes));
        const auto output =
            dest.first(static_cast<std::size_t>(to_read * per_sample));

        if (m_info.nbits <= 4) {
            if constexpr (std::same_as<T, std::uint16_t>) {
                throw ValidationError(
                    "psrio: uint16_t output requires 16-bit samples");
            } else {
                detail::unpack_sub_byte(raw, output, m_info.nbits,
                                        m_info.bit_order);
            }
        } else if (m_info.nbits == 8) {
            if constexpr (std::same_as<T, std::uint16_t>) {
                throw ValidationError(
                    "psrio: uint16_t output requires 16-bit samples");
            } else {
                detail::unpack_8bit(raw, output, m_info.samples_signed);
            }
        } else if (m_info.nbits == 16) {
            if constexpr (std::same_as<T, std::uint8_t>) {
                throw ValidationError(
                    "psrio: uint8_t output requires nbits of 8 or less");
            } else {
                detail::unpack_16le(raw, output);
            }
        } else if constexpr (std::same_as<T, float>) {
            detail::unpack_32le(raw, output);
        } else {
            throw ValidationError("psrio: 32-bit samples unpack to float");
        }

        m_sample += to_read;
        return to_read;
    }

    [[nodiscard]] std::vector<std::byte> read_bytes(std::uint64_t nbytes) {
        std::vector<std::byte> out(static_cast<std::size_t>(nbytes));
        const auto actual = read_bytes(nbytes, std::span<std::byte>(out));
        out.resize(static_cast<std::size_t>(actual));
        return out;
    }

    [[nodiscard]] std::vector<std::byte> read_block(std::uint64_t count) {
        const auto bytes_needed = count * m_stride;
        std::vector<std::byte> out(static_cast<std::size_t>(bytes_needed));
        const auto actual = read_block(count, std::span<std::byte>(out));
        out.resize(static_cast<std::size_t>(actual * m_stride));
        return out;
    }

    template <typename T = float>
    [[nodiscard]] std::vector<T> read_samples(std::uint64_t count) {
        const auto values_needed = count * m_info.nifs * m_info.nchans;
        std::vector<T> out(static_cast<std::size_t>(values_needed));
        const auto actual = read_samples(count, std::span<T>(out));
        out.resize(
            static_cast<std::size_t>(actual * m_info.nifs * m_info.nchans));
        return out;
    }

private:
    void init() {
        if (m_info.nchans == 0U || m_info.nifs == 0U) {
            throw ValidationError(
                "psrio: memory block needs at least one channel and IF");
        }
        m_type = sample_type_from_nbits(m_info.nbits, m_info.samples_signed);
        const auto width = static_cast<std::uint64_t>(m_info.nbits);
        if (m_info.nchans >
                std::numeric_limits<std::uint64_t>::max() / m_info.nifs ||
            (m_info.nchans * m_info.nifs) >
                std::numeric_limits<std::uint64_t>::max() / width) {
            throw ValidationError("psrio: memory block geometry is too large");
        }
        const auto bits = m_info.nchans * m_info.nifs * width;
        if (bits % 8U != 0U) {
            throw ValidationError(
                "psrio: nifs * nchans * nbits must form a byte stride");
        }
        m_stride = bits / 8U;
        if (m_stride == 0U || m_bytes.size() % m_stride != 0U) {
            throw ValidationError(
                "psrio: memory block byte count must be a multiple of the "
                "sample stride");
        }
        m_nsamples = m_bytes.size() / m_stride;
    }

    MemoryInfo m_info;
    SampleType m_type{SampleType::kUInt8};
    std::vector<std::byte> m_owned;
    std::span<const std::byte> m_bytes;
    std::uint64_t m_stride{0};
    std::uint64_t m_nsamples{0};
    std::uint64_t m_sample{0};
    bool m_apply_fswap{false};
};

} // namespace psrio
