#pragma once

#include "psrio/detail/endian.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/mmap.hpp"
#include "psrio/formats/presto/inf.hpp"

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <format>
#include <span>
#include <string>
#include <vector>

namespace psrio::formats::presto {

/**
 * @brief Streaming reader for PRESTO time-series (`.dat` paired with `.inf`).
 *
 * The raw 32-bit float `.dat` file is memory-mapped read-only. Samples are
 * streamed chunk-by-chunk directly into caller-managed buffers or accessed
 * zero-copy via view_samples().
 */
class TimeSeriesReader {
public:
    explicit TimeSeriesReader(const std::filesystem::path& path);

    TimeSeriesReader(const TimeSeriesReader&)                = delete;
    TimeSeriesReader& operator=(const TimeSeriesReader&)     = delete;
    TimeSeriesReader(TimeSeriesReader&&) noexcept            = default;
    TimeSeriesReader& operator=(TimeSeriesReader&&) noexcept = default;
    ~TimeSeriesReader()                                      = default;

    [[nodiscard]] const Header& header() const noexcept { return m_header; }

    /// Sample index of the next read.
    [[nodiscard]] std::uint64_t tell() const noexcept { return m_sample; }

    void rewind() noexcept { m_sample = 0; }

    /// Move the cursor to @p sample.
    /// @throws ValidationError if @p sample is greater than nsamples().
    void seek(std::uint64_t sample);

    /// Read next @p count time samples into @p dest.
    /// @p dest must hold exactly @p count float values.
    /// Returns the number of time samples read (may be less than count if at
    /// EOF).
    /// @throws ValidationError if @p dest has the wrong size.
    std::uint64_t read_samples(std::uint64_t count, std::span<float> dest);

    /// Alias for read_samples.
    std::uint64_t read(std::uint64_t count, std::span<float> dest) {
        return read_samples(count, dest);
    }

    /// Copy the next @p nbytes raw payload bytes into @p dest.
    /// @throws ValidationError if @p dest has a different size or request
    ///         extends past readable data.
    std::uint64_t read_bytes(std::uint64_t nbytes, std::span<std::byte> dest);

    /// Zero-copy view of the next @p count 32-bit float samples directly in
    /// mapped memory. Advances the cursor by the number of samples viewed.
    [[nodiscard]] std::span<const float> view_samples(std::uint64_t count);

    /// Zero-copy view alias for view_samples.
    [[nodiscard]] std::span<const float> view(std::uint64_t count) {
        return view_samples(count);
    }

    /// Zero-copy view of the next @p nbytes raw payload bytes.
    [[nodiscard]] std::span<const std::byte> view_bytes(std::uint64_t nbytes);

    /// Convenience allocating read: read next @p count samples into
    /// std::vector<float>.
    [[nodiscard]] std::vector<float> read_samples(std::uint64_t count) {
        std::vector<float> buffer(count);
        const auto actual = read_samples(count, std::span<float>(buffer));
        buffer.resize(static_cast<std::size_t>(actual));
        return buffer;
    }

    /// Convenience allocating read for raw bytes.
    [[nodiscard]] std::vector<std::byte> read_bytes(std::uint64_t nbytes) {
        std::vector<std::byte> buffer(nbytes);
        const auto actual = read_bytes(nbytes, std::span<std::byte>(buffer));
        buffer.resize(static_cast<std::size_t>(actual));
        return buffer;
    }

    /// Non-fatal warnings encountered during opening or parsing.
    [[nodiscard]] const std::vector<std::string>& warnings() const noexcept {
        return m_warnings;
    }

private:
    ::psrio::detail::MappedFile m_file;
    Header m_header;
    std::uint64_t m_sample{0};
    std::uint64_t m_available_samples{0};
    std::vector<std::string> m_warnings;
};

inline TimeSeriesReader::TimeSeriesReader(const std::filesystem::path& path) {
    std::filesystem::path dat_path = path;
    std::filesystem::path inf_path = path;

    if (path.extension() == ".inf") {
        dat_path.replace_extension(".dat");
    } else if (path.extension() == ".dat") {
        inf_path.replace_extension(".inf");
    } else {
        dat_path = path.string() + ".dat";
        inf_path = path.string() + ".inf";
    }

    if (std::filesystem::exists(inf_path)) {
        m_header = parse_inf_file(inf_path);
    } else {
        m_header.filename = dat_path.stem().string();
    }

    try {
        m_file = ::psrio::detail::MappedFile(dat_path);
        m_file.advise_sequential();
    } catch (const std::exception& ex) {
        throw IoError(ex.what());
    }

    if (m_file.empty()) {
        throw FormatError(std::format("psrio: PRESTO .dat file '{}' is empty",
                                      dat_path.string()));
    }

    m_available_samples = m_file.size() / sizeof(float);
    if (m_header.nsamples == 0) {
        m_header.nsamples = m_available_samples;
    } else if (m_header.nsamples > m_available_samples) {
        m_warnings.push_back(
            std::format("psrio: .inf declared {} samples but file only has {}",
                        m_header.nsamples, m_available_samples));
        m_header.nsamples = m_available_samples;
    }
}

inline void TimeSeriesReader::seek(std::uint64_t sample) {
    if (sample > m_header.nsamples) {
        throw ValidationError(std::format(
            "psrio: cannot seek to sample {} (total readable samples: {})",
            sample, m_header.nsamples));
    }
    m_sample = sample;
}

inline std::uint64_t TimeSeriesReader::read_samples(std::uint64_t count,
                                                    std::span<float> dest) {
    if (dest.size() != count) {
        throw ValidationError(std::format(
            "psrio: destination has {} values but the request needs {}",
            dest.size(), count));
    }
    if (m_sample > m_header.nsamples) {
        throw ValidationError(
            "psrio: reader cursor is past the readable samples");
    }

    const auto to_read = std::min(count, m_header.nsamples - m_sample);
    if (to_read == 0U) {
        return 0;
    }

    const auto byte_offset = m_sample * sizeof(float);
    const auto byte_count  = to_read * sizeof(float);
    const auto raw         = m_file.bytes().subspan(byte_offset, byte_count);

    if constexpr (std::endian::native == std::endian::little) {
        std::memcpy(dest.data(), raw.data(), byte_count);
    } else {
        for (std::size_t i = 0; i < to_read; ++i) {
            dest[i] = ::psrio::detail::load_little_endian<float>(
                raw.data() + (i * sizeof(float)));
        }
    }

    m_sample += to_read;
    return to_read;
}

inline std::uint64_t TimeSeriesReader::read_bytes(std::uint64_t nbytes,
                                                  std::span<std::byte> dest) {
    if (dest.size() != nbytes) {
        throw ValidationError(std::format(
            "psrio: destination has {} bytes but the request needs {}",
            dest.size(), nbytes));
    }
    const auto sample_stride = sizeof(float);
    if (nbytes % sample_stride != 0) {
        throw ValidationError(
            "psrio: byte request must be a multiple of sample size (4 bytes)");
    }
    const auto samples_requested = nbytes / sample_stride;
    if (m_sample + samples_requested > m_header.nsamples) {
        throw ValidationError(
            "psrio: byte request extends past readable samples");
    }
    const auto byte_offset = m_sample * sample_stride;
    const auto raw         = m_file.bytes().subspan(byte_offset, nbytes);
    std::memcpy(dest.data(), raw.data(), nbytes);
    m_sample += samples_requested;
    return nbytes;
}

inline std::span<const float>
TimeSeriesReader::view_samples(std::uint64_t count) {
    if (m_sample + count > m_header.nsamples) {
        throw ValidationError(
            "psrio: view sample request extends past readable samples");
    }
    const auto byte_offset = m_sample * sizeof(float);
    const auto byte_count  = count * sizeof(float);
    const auto raw         = m_file.bytes().subspan(byte_offset, byte_count);
    m_sample += count;
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    return {reinterpret_cast<const float*>(raw.data()),
            static_cast<std::size_t>(count)};
}

inline std::span<const std::byte>
TimeSeriesReader::view_bytes(std::uint64_t nbytes) {
    const auto sample_stride = sizeof(float);
    if (nbytes % sample_stride != 0) {
        throw ValidationError(
            "psrio: byte request must be a multiple of sample size (4 bytes)");
    }
    const auto samples_requested = nbytes / sample_stride;
    if (m_sample + samples_requested > m_header.nsamples) {
        throw ValidationError(
            "psrio: byte request extends past readable samples");
    }
    const auto byte_offset = m_sample * sample_stride;
    const auto raw         = m_file.bytes().subspan(byte_offset, nbytes);
    m_sample += samples_requested;
    return raw;
}

} // namespace psrio::formats::presto
