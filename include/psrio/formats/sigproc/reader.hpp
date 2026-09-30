#pragma once

#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/mmap.hpp"
#include "psrio/detail/unpack.hpp"
#include "psrio/formats/sigproc/header.hpp"

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <format>
#include <limits>
#include <span>
#include <string>

namespace psrio::formats::sigproc {

/// Sequential reader for one SIGPROC filterbank file.
///
/// The file is mapped read-only. read(SampleCount) unpacks the next time
/// samples into a caller-owned span. read(ByteCount) and view(ByteCount)
/// expose the raw payload. One time sample is nifs * nchans values in file
/// order, index ipol * nchans + ichan, with no scaling and no Stokes offset.
/// Sub-byte samples are least-significant field first.
///
/// The cursor starts at the first readable sample. seek(nsamples()) is the
/// end. A sample read past the end returns the samples that remain. A byte
/// request that is not a positive multiple of the stride, or that would pass
/// the last readable sample, throws.
class FilterbankReader {
public:
    explicit FilterbankReader(const std::filesystem::path& path,
                              BitOrder bit_order = BitOrder::kLsbFirst);

    FilterbankReader(const FilterbankReader&)                = delete;
    FilterbankReader& operator=(const FilterbankReader&)     = delete;
    FilterbankReader(FilterbankReader&&) noexcept            = default;
    FilterbankReader& operator=(FilterbankReader&&) noexcept = default;
    ~FilterbankReader()                                      = default;

    [[nodiscard]] const FilterbankHeader& header() const noexcept {
        return m_header;
    }

    /// Configured sub-byte bit packing order.
    [[nodiscard]] BitOrder bit_order() const noexcept { return m_bit_order; }

    /// Update the sub-byte bit packing order.
    void set_bit_order(BitOrder order) noexcept { m_bit_order = order; }

    /// Sample index of the next read.
    [[nodiscard]] std::uint64_t tell() const noexcept { return m_sample; }

    void rewind() noexcept { m_sample = 0; }

    /// Move the cursor to @p sample. @p sample may equal nsamples().
    /// @throws ValidationError if @p sample is greater than nsamples().
    void seek(std::uint64_t sample);

    /// Unpack the next @p count time samples into @p dest using @p order for
    /// sub-byte samples.
    ///
    /// @p dest must hold count * nifs * nchans values. T is float for every
    /// supported nbits, uint8_t for nbits <= 8, or uint16_t for nbits == 16.
    /// A short read leaves the unused tail of @p dest untouched and returns
    /// the number of time samples written. Zero means the cursor is at the end.
    /// @throws ValidationError if @p dest has the wrong size or T cannot hold
    ///         this nbits.
    template <typename T>
        requires std::same_as<T, float> || std::same_as<T, std::uint8_t> ||
                 std::same_as<T, std::uint16_t>
    std::uint64_t read(SampleCount count, std::span<T> dest, BitOrder order);

    /// Unpack the next @p count time samples into @p dest using the reader's
    /// bit_order().
    template <typename T>
        requires std::same_as<T, float> || std::same_as<T, std::uint8_t> ||
                 std::same_as<T, std::uint16_t>
    std::uint64_t read(SampleCount count, std::span<T> dest) {
        return read(count, dest, m_bit_order);
    }

    /// Copy the next @p count.value raw payload bytes into @p dest.
    /// @throws ValidationError if the request is not a positive multiple of
    ///         the sample stride, does not fit in the readable samples, or
    ///         @p dest has a different size.
    std::uint64_t read(ByteCount count, std::span<std::byte> dest);

    /// Zero-copy view of the next @p count.value raw payload bytes.
    /// The span refers to the mapped file and stays valid until that mapping
    /// is destroyed. Advances the cursor by the same rule as read(ByteCount).
    [[nodiscard]] std::span<const std::byte> view(ByteCount count);

    /// Convenience allocating read: unpack next @p count time samples into a
    /// new std::vector<T>.
    template <typename T = float>
        requires std::same_as<T, float> || std::same_as<T, std::uint8_t> ||
                 std::same_as<T, std::uint16_t>
    [[nodiscard]] std::vector<T> read_samples(SampleCount count,
                                              BitOrder order) {
        const auto per_sample = values_per_sample();
        std::vector<T> buffer(count.value * per_sample);
        const auto actual = read(count, std::span<T>(buffer), order);
        buffer.resize(static_cast<std::size_t>(actual * per_sample));
        return buffer;
    }

    /// Convenience allocating read using the reader's bit_order().
    template <typename T = float>
        requires std::same_as<T, float> || std::same_as<T, std::uint8_t> ||
                 std::same_as<T, std::uint16_t>
    [[nodiscard]] std::vector<T> read_samples(SampleCount count) {
        return read_samples<T>(count, m_bit_order);
    }

    /// Convenience allocating read for raw bytes.
    [[nodiscard]] std::vector<std::byte> read_bytes(ByteCount count) {
        std::vector<std::byte> buffer(count.value);
        const auto actual = read(count, std::span<std::byte>(buffer));
        buffer.resize(static_cast<std::size_t>(actual));
        return buffer;
    }

private:
    [[nodiscard]] std::uint64_t values_per_sample() const;

    void check_byte_request(std::uint64_t nbytes) const;

    [[nodiscard]] std::span<const std::byte>
    mapped_bytes(std::uint64_t offset, std::uint64_t nbytes) const;

    ::psrio::detail::MappedFile m_file;
    FilterbankHeader m_header;
    std::span<const std::byte> m_payload;
    std::uint64_t m_sample{0};
    BitOrder m_bit_order{BitOrder::kLsbFirst};
};

inline FilterbankReader::FilterbankReader(const std::filesystem::path& path,
                                          BitOrder bit_order)
    : m_bit_order(bit_order) {
    try {
        m_file = ::psrio::detail::MappedFile(path);
        m_file.advise_sequential();
    } catch (const std::exception& ex) {
        throw IoError(ex.what());
    }

    const std::string path_string = path.string();
    m_header = FilterbankHeader::parse(m_file.bytes(), path_string);
    try {
        m_header.validate();
    } catch (const ValidationError& ex) {
        throw ValidationError(path_string + ": " + ex.what());
    }

    const auto mapped = m_file.bytes();
    if (m_header.header_bytes > mapped.size()) {
        throw FormatError(path_string + ": header extends past end of file");
    }
    m_payload = mapped.subspan(static_cast<std::size_t>(m_header.header_bytes));
}

inline void FilterbankReader::seek(std::uint64_t sample) {
    if (sample > m_header.nsamples()) {
        throw ValidationError("psrio: seek is past the readable samples");
    }
    m_sample = sample;
}

inline std::uint64_t FilterbankReader::values_per_sample() const {
    return static_cast<std::uint64_t>(m_header.nifs) *
           static_cast<std::uint64_t>(m_header.nchans);
}

inline std::span<const std::byte>
FilterbankReader::mapped_bytes(std::uint64_t offset,
                               std::uint64_t nbytes) const {
    if (offset > m_payload.size() || nbytes > m_payload.size() - offset) {
        throw FormatError("psrio: sample offset is outside the mapped payload");
    }
    return m_payload.subspan(static_cast<std::size_t>(offset),
                             static_cast<std::size_t>(nbytes));
}

inline void FilterbankReader::check_byte_request(std::uint64_t nbytes) const {
    const auto stride = m_header.bytes_per_sample();
    if (nbytes == 0U || stride == 0U || nbytes % stride != 0U) {
        throw ValidationError("psrio: byte request must be a positive multiple "
                              "of the sample stride");
    }
    if (m_sample > m_header.nsamples()) {
        throw ValidationError(
            "psrio: reader cursor is past the readable samples");
    }
    const auto samples   = nbytes / stride;
    const auto available = m_header.nsamples() - m_sample;
    if (samples > available) {
        throw ValidationError(
            "psrio: byte request extends past the readable samples");
    }
}

template <typename T>
    requires std::same_as<T, float> || std::same_as<T, std::uint8_t> ||
             std::same_as<T, std::uint16_t>
inline std::uint64_t
FilterbankReader::read(SampleCount count, std::span<T> dest, BitOrder order) {
    if constexpr (std::same_as<T, std::uint8_t>) {
        if (m_header.nbits > 8) {
            throw ValidationError(
                "psrio: uint8_t output requires nbits of 8 or less");
        }
    } else if constexpr (std::same_as<T, std::uint16_t>) {
        if (m_header.nbits != 16) {
            throw ValidationError(
                "psrio: uint16_t output requires 16-bit samples");
        }
    }

    const auto per_sample = values_per_sample();
    if (count.value > 0U &&
        per_sample >
            (std::numeric_limits<std::uint64_t>::max() / count.value)) {
        throw ValidationError("psrio: requested sample block is too large");
    }
    const auto expected = count.value * per_sample;
    if (expected != dest.size()) {
        throw ValidationError(std::format(
            "psrio: destination has {} values but the request needs {}",
            dest.size(), expected));
    }
    if (m_sample > m_header.nsamples()) {
        throw ValidationError(
            "psrio: reader cursor is past the readable samples");
    }

    const auto to_read = std::min(count.value, m_header.nsamples() - m_sample);
    if (to_read == 0U) {
        return 0;
    }

    const auto stride = m_header.bytes_per_sample();
    const auto raw    = mapped_bytes(m_sample * stride, to_read * stride);
    const auto output =
        dest.first(static_cast<std::size_t>(to_read * per_sample));

    if (m_header.nbits <= 4) {
        if constexpr (std::same_as<T, std::uint16_t>) {
            throw ValidationError(
                "psrio: uint16_t output requires 16-bit samples");
        } else {
            ::psrio::detail::unpack_sub_byte(raw, output, m_header.nbits,
                                             order);
        }
    } else if (m_header.nbits == 8) {
        if constexpr (std::same_as<T, std::uint16_t>) {
            throw ValidationError(
                "psrio: uint16_t output requires 16-bit samples");
        } else {
            ::psrio::detail::unpack_8bit(raw, output,
                                         m_header.samples_are_signed());
        }
    } else if (m_header.nbits == 16) {
        if constexpr (std::same_as<T, std::uint8_t>) {
            throw ValidationError(
                "psrio: uint8_t output requires nbits of 8 or less");
        } else {
            ::psrio::detail::unpack_16le(raw, output);
        }
    } else if constexpr (std::same_as<T, float>) {
        ::psrio::detail::unpack_32le(raw, output);
    } else {
        throw ValidationError("psrio: 32-bit samples unpack to float");
    }

    m_sample += to_read;
    return to_read;
}

inline std::uint64_t FilterbankReader::read(ByteCount count,
                                            std::span<std::byte> dest) {
    check_byte_request(count.value);
    if (dest.size() != count.value) {
        throw ValidationError(std::format(
            "psrio: destination has {} bytes but the request needs {}",
            dest.size(), count.value));
    }
    const auto stride = m_header.bytes_per_sample();
    const auto raw    = mapped_bytes(m_sample * stride, count.value);
    std::memcpy(dest.data(), raw.data(), static_cast<std::size_t>(count.value));
    m_sample += count.value / stride;
    return count.value;
}

inline std::span<const std::byte> FilterbankReader::view(ByteCount count) {
    check_byte_request(count.value);
    const auto stride = m_header.bytes_per_sample();
    const auto raw    = mapped_bytes(m_sample * stride, count.value);
    m_sample += count.value / stride;
    return raw;
}

} // namespace psrio::formats::sigproc
