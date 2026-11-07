#pragma once

#include "psrio/common/types.hpp"
#include "psrio/detail/endian.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/packed_bits.hpp"
#include "psrio/formats/sigproc/header.hpp"

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fmt/format.h>
#include <fstream>
#include <limits>
#include <span>
#include <utility>
#include <vector>

namespace psrio::formats::sigproc {

/// Sequential writer for one SIGPROC filterbank file.
///
/// The constructor validates the header, truncates @p path, and writes the
/// binary header. The legacy nsamples key is never written: the sample count
/// is the payload length, which a reader derives from the file size. Each
/// write appends whole time samples of nifs * nchans values, index
/// ipol * nchans + ichan, in the on-disk layout that FilterbankReader reads.
///
/// write_samples() quantises float input by rounding to nearest and clamping
/// to the range of the output type; NaN becomes the lower bound. 8-bit output
/// is int8 when the header's signed flag is set. Sub-byte output is packed in
/// the configured BitOrder.
class FilterbankWriter {
public:
    /// @throws ValidationError if @p header does not describe readable
    ///         filterbank data.
    /// @throws IoError if @p path cannot be created or the header write fails.
    explicit FilterbankWriter(std::filesystem::path path,
                              FilterbankHeader header,
                              BitOrder bit_order = BitOrder::kLsbFirst);

    FilterbankWriter(const FilterbankWriter&)                = delete;
    FilterbankWriter& operator=(const FilterbankWriter&)     = delete;
    FilterbankWriter(FilterbankWriter&&) noexcept            = default;
    FilterbankWriter& operator=(FilterbankWriter&&) noexcept = default;

    /// Closes the file. A failure while closing is ignored; call close() to
    /// observe it.
    ~FilterbankWriter() noexcept;

    /// Header as written. header_bytes, data_bytes, and file_bytes track the
    /// bytes written so far.
    [[nodiscard]] const FilterbankHeader& header() const noexcept {
        return m_header;
    }

    [[nodiscard]] const std::filesystem::path& path() const noexcept {
        return m_path;
    }

    /// Time samples written so far.
    [[nodiscard]] std::uint64_t nsamples_written() const noexcept {
        return m_samples_written;
    }

    /// Bytes in one time sample (`nifs * nchans * nbits / 8`).
    [[nodiscard]] std::uint64_t bytes_per_sample() const noexcept {
        return m_stride;
    }

    [[nodiscard]] std::uint64_t nchans() const noexcept {
        return static_cast<std::uint64_t>(m_header.nchans);
    }

    [[nodiscard]] std::uint64_t nifs() const noexcept {
        return static_cast<std::uint64_t>(m_header.nifs);
    }

    [[nodiscard]] int nbits() const noexcept { return m_header.nbits; }

    [[nodiscard]] BitOrder bit_order() const noexcept { return m_bit_order; }

    void set_bit_order(BitOrder order) noexcept { m_bit_order = order; }

    [[nodiscard]] bool is_open() const noexcept { return m_out.is_open(); }

    /// Write @p count packed time samples from the front of @p data.
    /// @throws ValidationError if @p data holds fewer than
    ///         count * bytes_per_sample() bytes.
    /// @throws IoError if the writer is closed or the write fails.
    void write_block(std::uint64_t count, std::span<const std::byte> data);

    /// Write raw payload bytes.
    /// @throws ValidationError if the size is not a multiple of
    ///         bytes_per_sample().
    /// @throws IoError if the writer is closed or the write fails.
    void write_bytes(std::span<const std::byte> data);

    /// Encode and write @p count unpacked time samples from the front of
    /// @p samples, which must hold count * nifs * nchans values.
    ///
    /// float accepts every nbits. uint8_t requires nbits of 8 or less and is
    /// masked to nbits. uint16_t requires 16-bit samples.
    /// @throws ValidationError if the type or size does not fit.
    /// @throws IoError if the writer is closed or the write fails.
    template <typename T>
        requires std::same_as<T, float> || std::same_as<T, std::uint8_t> ||
                 std::same_as<T, std::uint16_t>
    void write_samples(std::uint64_t count,
                       std::span<const T> samples,
                       BitOrder order);

    template <typename T>
        requires std::same_as<T, float> || std::same_as<T, std::uint8_t> ||
                 std::same_as<T, std::uint16_t>
    void write_samples(std::uint64_t count, std::span<const T> samples) {
        write_samples(count, samples, m_bit_order);
    }

    /// Alias for write_samples().
    template <typename T>
        requires std::same_as<T, float> || std::same_as<T, std::uint8_t> ||
                 std::same_as<T, std::uint16_t>
    void write(std::uint64_t count, std::span<const T> samples) {
        write_samples(count, samples, m_bit_order);
    }

    /// Flush buffered bytes to the file.
    /// @throws IoError if the writer is closed or the flush fails.
    void flush();

    /// Flush and close the file. Does nothing when already closed.
    /// @throws IoError if the flush or close fails.
    void close();

private:
    static constexpr std::uint64_t kChunkBytes = std::uint64_t{1} << 22U;

    void require_open() const;
    void append(std::span<const std::byte> bytes, std::uint64_t nsamps);
    [[nodiscard]] std::uint64_t checked_bytes(std::uint64_t count) const;

    template <typename T>
    void encode(std::span<const T> values, BitOrder order);

    std::ofstream m_out;
    std::filesystem::path m_path;
    FilterbankHeader m_header;
    std::uint64_t m_stride{0};
    std::uint64_t m_samples_written{0};
    BitOrder m_bit_order{BitOrder::kLsbFirst};
    std::vector<std::byte> m_bytes;
    std::vector<std::uint8_t> m_levels;
};

namespace writer_detail {

/// Round @p value to the nearest integer in [lo, hi]. NaN maps to lo.
template <typename Out>
[[nodiscard]] inline Out quantise(float value, float lo, float hi) noexcept {
    return static_cast<Out>(std::round(std::fmin(std::fmax(value, lo), hi)));
}

template <typename T>
inline void store_le(std::span<const T> values, std::span<std::byte> out) {
    for (std::size_t index = 0; index < values.size(); ++index) {
        const T stored = ::psrio::detail::to_little_endian(values[index]);
        std::memcpy(out.data() + (index * sizeof(T)), &stored, sizeof(T));
    }
}

} // namespace writer_detail

inline FilterbankWriter::FilterbankWriter(std::filesystem::path path,
                                          FilterbankHeader header,
                                          BitOrder bit_order)
    : m_path(std::move(path)),
      m_header(std::move(header)),
      m_bit_order(bit_order) {
    m_header.declared_nsamples.reset();
    m_header.validate();
    (void)sample_type_from_nbits(m_header.nbits, m_header.samples_are_signed());
    m_stride              = m_header.bytes_per_sample();
    const auto encoded    = m_header.to_bytes();
    m_header.header_bytes = encoded.size();
    m_header.data_bytes   = 0;
    m_header.file_bytes   = encoded.size();

    m_out.open(m_path, std::ios::binary | std::ios::trunc);
    if (!m_out) {
        throw IoError(fmt::format("psrio: cannot create {}", m_path.string()));
    }
    m_out.write(reinterpret_cast<const char*>(encoded.data()),
                static_cast<std::streamsize>(encoded.size()));
    if (!m_out) {
        throw IoError(fmt::format("psrio: failed to write header to {}",
                                  m_path.string()));
    }
}

inline FilterbankWriter::~FilterbankWriter() noexcept {
    try {
        close();
    } catch (...) {
        // Destructors must not throw; close() reports the failure.
    }
}

inline void FilterbankWriter::require_open() const {
    if (!m_out.is_open()) {
        throw IoError(fmt::format("psrio: filterbank writer for {} is closed",
                                  m_path.string()));
    }
}

inline std::uint64_t
FilterbankWriter::checked_bytes(std::uint64_t count) const {
    if (count > std::numeric_limits<std::uint64_t>::max() / m_stride) {
        throw ValidationError("psrio: write sample count is too large");
    }
    return count * m_stride;
}

inline void FilterbankWriter::append(std::span<const std::byte> bytes,
                                     std::uint64_t nsamps) {
    m_out.write(reinterpret_cast<const char*>(bytes.data()),
                static_cast<std::streamsize>(bytes.size()));
    if (!m_out) {
        throw IoError(
            fmt::format("psrio: failed to write data to {}", m_path.string()));
    }
    m_samples_written += nsamps;
    m_header.data_bytes += bytes.size();
    m_header.file_bytes += bytes.size();
}

inline void FilterbankWriter::write_block(std::uint64_t count,
                                          std::span<const std::byte> data) {
    require_open();
    const auto nbytes = checked_bytes(count);
    if (data.size() < nbytes) {
        throw ValidationError(fmt::format(
            "psrio: write_block needs {} bytes for {} samples, got {}", nbytes,
            count, data.size()));
    }
    append(data.first(static_cast<std::size_t>(nbytes)), count);
}

inline void FilterbankWriter::write_bytes(std::span<const std::byte> data) {
    require_open();
    if (data.size() % m_stride != 0U) {
        throw ValidationError(fmt::format(
            "psrio: write_bytes size {} is not a multiple of the {}-byte "
            "sample stride",
            data.size(), m_stride));
    }
    append(data, data.size() / m_stride);
}

template <typename T>
    requires std::same_as<T, float> || std::same_as<T, std::uint8_t> ||
             std::same_as<T, std::uint16_t>
inline void FilterbankWriter::write_samples(std::uint64_t count,
                                            std::span<const T> samples,
                                            BitOrder order) {
    require_open();
    const int bits = m_header.nbits;
    if constexpr (std::same_as<T, std::uint8_t>) {
        if (bits > 8) {
            throw ValidationError(
                "psrio: uint8_t input requires nbits of 8 or less");
        }
    } else if constexpr (std::same_as<T, std::uint16_t>) {
        if (bits != 16) {
            throw ValidationError(
                "psrio: uint16_t input requires 16-bit samples");
        }
    }
    const auto values_per_sample = nifs() * nchans();
    if (samples.size() / values_per_sample < count) {
        throw ValidationError(fmt::format(
            "psrio: write_samples needs {} values per sample for {} samples, "
            "got {}",
            values_per_sample, count, samples.size()));
    }
    (void)checked_bytes(count);

    const auto chunk = std::max<std::uint64_t>(1U, kChunkBytes / m_stride);
    for (std::uint64_t done = 0; done < count;) {
        const auto nsamps = std::min(chunk, count - done);
        const auto values = samples.subspan(
            static_cast<std::size_t>(done * values_per_sample),
            static_cast<std::size_t>(nsamps * values_per_sample));
        const auto nbytes = static_cast<std::size_t>(nsamps * m_stride);

        const bool as_is = (std::same_as<T, std::uint8_t> && bits == 8) ||
                           (::psrio::detail::host_is_little_endian() &&
                            ((std::same_as<T, float> && bits == 32) ||
                             std::same_as<T, std::uint16_t>));
        if (as_is) {
            append(std::as_bytes(values), nsamps);
        } else {
            m_bytes.resize(nbytes);
            encode(values, order);
            append(std::span<const std::byte>(m_bytes).first(nbytes), nsamps);
        }
        done += nsamps;
    }
}

template <typename T>
inline void FilterbankWriter::encode(std::span<const T> values,
                                     BitOrder order) {
    using writer_detail::quantise;
    const int bits = m_header.nbits;
    const std::span<std::byte> out(m_bytes);

    if (bits == 32 || bits == 16) {
        if constexpr (std::same_as<T, float>) {
            if (bits == 32) {
                writer_detail::store_le(values, out);
                return;
            }
            for (std::size_t index = 0; index < values.size(); ++index) {
                const auto level = ::psrio::detail::to_little_endian(
                    quantise<std::uint16_t>(values[index], 0.0F, 65535.0F));
                std::memcpy(out.data() + (index * 2U), &level, 2U);
            }
        } else if constexpr (std::same_as<T, std::uint16_t>) {
            writer_detail::store_le(values, out);
        }
        return;
    }

    if constexpr (std::same_as<T, float>) {
        if (bits == 8) {
            if (m_header.samples_are_signed()) {
                for (std::size_t index = 0; index < values.size(); ++index) {
                    out[index] = static_cast<std::byte>(
                        static_cast<std::uint8_t>(quantise<std::int8_t>(
                            values[index], -128.0F, 127.0F)));
                }
            } else {
                for (std::size_t index = 0; index < values.size(); ++index) {
                    out[index] = static_cast<std::byte>(
                        quantise<std::uint8_t>(values[index], 0.0F, 255.0F));
                }
            }
            return;
        }
        const auto hi =
            static_cast<float>((1U << static_cast<unsigned>(bits)) - 1U);
        m_levels.resize(values.size());
        for (std::size_t index = 0; index < values.size(); ++index) {
            m_levels[index] = quantise<std::uint8_t>(values[index], 0.0F, hi);
        }
        ::psrio::detail::pack_sub_byte(m_levels, out, bits, order);
    } else if constexpr (std::same_as<T, std::uint8_t>) {
        ::psrio::detail::pack_sub_byte(values, out, bits, order);
    }
}

inline void FilterbankWriter::flush() {
    require_open();
    m_out.flush();
    if (!m_out) {
        throw IoError(
            fmt::format("psrio: failed to flush {}", m_path.string()));
    }
}

inline void FilterbankWriter::close() {
    if (!m_out.is_open()) {
        return;
    }
    m_out.close();
    if (!m_out) {
        throw IoError(
            fmt::format("psrio: failed to close {}", m_path.string()));
    }
}

} // namespace psrio::formats::sigproc
