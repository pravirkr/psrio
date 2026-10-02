#pragma once

#include "psrio/common/types.hpp"
#include "psrio/detail/endian.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/mmap.hpp"
#include "psrio/detail/packed_bits.hpp"
#include "psrio/detail/skip.hpp"
#include "psrio/formats/sigproc/header.hpp"
#include "psrio/header.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <format>
#include <fstream>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace psrio {

class TimeSeries;

/**
 * @brief Zero-overhead, streaming reader for single-channel time series (.tim,
 * .dat/.inf).
 *
 * Implements the lite time-series interface: reads and parses metadata on
 * construction, while allowing downstream callers to read directly into
 * caller-owned memory buffers (std::span<float>) with zero library-side heap
 * allocations.
 *
 * Models the psrio::concepts::TimeSeriesReader concept.
 */
class TimeSeriesReader {
public:
    explicit TimeSeriesReader(const std::filesystem::path& path,
                              BitOrder bit_order = BitOrder::kLsbFirst)
        : m_bit_order(bit_order) {
        open_file(path);
    }

    TimeSeriesReader(const std::filesystem::path& dat_path,
                     const std::filesystem::path& inf_path,
                     BitOrder bit_order = BitOrder::kLsbFirst)
        : m_bit_order(bit_order) {
        open_presto(dat_path, inf_path);
    }

    TimeSeriesReader(const TimeSeriesReader&)                = delete;
    TimeSeriesReader& operator=(const TimeSeriesReader&)     = delete;
    TimeSeriesReader(TimeSeriesReader&&) noexcept            = default;
    TimeSeriesReader& operator=(TimeSeriesReader&&) noexcept = default;
    ~TimeSeriesReader()                                      = default;

    /// Observational metadata header.
    [[nodiscard]] const Header& header() const noexcept { return m_header; }

    /// Number of time samples available to read.
    [[nodiscard]] std::uint64_t nsamples() const noexcept {
        return m_header.nsamples;
    }

    /// Sampling interval in seconds.
    [[nodiscard]] double dt() const noexcept { return m_header.tsamp; }

    /// Sampling interval in seconds. Same value as `dt()`.
    [[nodiscard]] double tsamp() const noexcept { return m_header.tsamp; }

    /// Total duration in seconds.
    [[nodiscard]] double tobs() const noexcept { return m_header.tobs(); }

    /// Sub-byte bit order.
    [[nodiscard]] BitOrder bit_order() const noexcept { return m_bit_order; }

    void set_bit_order(BitOrder order) noexcept { m_bit_order = order; }

    /// Non-fatal parsing and file-size warnings.
    [[nodiscard]] const std::vector<std::string>& warnings() const noexcept {
        return m_warnings;
    }

    /// Current sample cursor index.
    [[nodiscard]] std::uint64_t tell() const noexcept { return m_sample; }

    /// Reset read cursor to sample 0.
    void rewind() noexcept { m_sample = 0; }

    /// Move cursor to @p sample.
    void seek(std::uint64_t sample) {
        if (sample > m_header.nsamples) {
            throw ValidationError("psrio: seek is past the readable samples");
        }
        m_sample = sample;
    }

    /// Move the cursor by @p delta time samples. Negative moves backward.
    /// A sub-byte series must land on a byte boundary or on `nsamples()`.
    /// @throws ValidationError if the landing index is outside
    ///         `[0, nsamples()]` or is not byte-aligned.
    void skip(std::int64_t delta) {
        const auto landing =
            detail::apply_skip(m_sample, delta, m_header.nsamples);
        if (is_sub_byte() && landing != m_header.nsamples &&
            landing % samples_per_byte() != 0U) {
            throw ValidationError(
                "psrio: sub-byte time series skip must land on a byte "
                "boundary");
        }
        m_sample = landing;
    }

    /// One frequency channel.
    static std::uint64_t nchans() noexcept { return 1; }

    /// IF count from the header, at least 1.
    [[nodiscard]] std::uint64_t nifs() const noexcept {
        return static_cast<std::uint64_t>(m_header.nifs > 0 ? m_header.nifs
                                                            : 1);
    }

    /// Bits per on-disk sample. PRESTO `.dat` samples are 32-bit floats.
    [[nodiscard]] int nbits() const noexcept {
        return m_is_presto ? 32 : m_raw_sigproc.nbits;
    }

    /// On-disk sample type.
    [[nodiscard]] SampleType sample_type() const {
        if (m_is_presto) {
            return SampleType::kFloat32;
        }
        return sample_type_from_nbits(m_raw_sigproc.nbits,
                                      m_raw_sigproc.samples_are_signed());
    }

    /// Bytes charged to one sample by `read_bytes`.
    ///
    /// PRESTO and whole-byte SIGPROC samples return their real width. A
    /// sub-byte series returns 1, which is the `read_bytes` stride.
    /// `read_block` still counts time samples and requires a byte-aligned
    /// count.
    [[nodiscard]] std::uint64_t bytes_per_sample() const {
        if (m_is_presto) {
            return sizeof(float);
        }
        if (m_raw_sigproc.nbits >= 8) {
            return static_cast<std::uint64_t>(m_raw_sigproc.nbits) / 8U;
        }
        return 1;
    }

    /// True. The file length is known after open.
    static bool has_nsamples() noexcept { return true; }

    /// Source name from header.
    [[nodiscard]] std::string_view source_name() const noexcept {
        return m_header.source;
    }

    /// Telescope name from header.
    [[nodiscard]] std::string_view telescope() const noexcept {
        return m_header.telescope;
    }

    /// Right ascension packed sexagesimal (HHMMSS.ss), or 0.0.
    [[nodiscard]] double raj() const noexcept { return m_header.raj; }

    /// Declination packed sexagesimal (+-DDMMSS.ss), or 0.0.
    [[nodiscard]] double dej() const noexcept { return m_header.dej; }

    /// Bandwidth in MHz.
    [[nodiscard]] double bandwidth() const noexcept {
        return m_header.bandwidth();
    }

    /// Centre frequency in MHz.
    [[nodiscard]] double center_frequency() const noexcept {
        return m_header.center_frequency();
    }

    /// MJD of the first sample.
    [[nodiscard]] double tstart() const noexcept { return m_header.tstart; }

    /// Centre frequency of the single channel, in MHz.
    [[nodiscard]] double fch1() const noexcept { return m_header.fch1; }

    /// Signed channel spacing, in MHz.
    [[nodiscard]] double foff() const noexcept { return m_header.foff; }

    /// Beam index from the header.
    [[nodiscard]] int beam() const noexcept { return m_header.ibeam; }

    /// Spectra per second, or 0 when `tsamp` is not positive.
    [[nodiscard]] double spectra_rate() const noexcept {
        return m_header.tsamp > 0.0 ? 1.0 / m_header.tsamp : 0.0;
    }

    /// UTC of `tstart`, as POSIX seconds.
    [[nodiscard]] std::time_t utc_start() const noexcept {
        return astro::mjd_to_time(m_header.tstart);
    }

    /// Read the complete time-series directly into a caller-owned destination
    /// buffer.
    std::uint64_t read_data(std::span<float> dest) {
        if (dest.size() != m_header.nsamples) {
            throw ValidationError(
                std::format("psrio: destination buffer size ({}) does not "
                            "match time-series nsamples ({})",
                            dest.size(), m_header.nsamples));
        }
        rewind();
        return read_samples(m_header.nsamples, dest);
    }

    /// Read next @p count time samples into caller-provided @p dest using @p
    /// order.
    std::uint64_t
    read_samples(std::uint64_t count, std::span<float> dest, BitOrder order) {
        if (count != dest.size()) {
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

        const auto output = dest.first(static_cast<std::size_t>(to_read));

        if (m_is_presto) {
            const auto byte_offset = m_sample * sizeof(float);
            const auto byte_count  = to_read * sizeof(float);
            const auto raw         = mapped_bytes(byte_offset, byte_count);

            if constexpr (std::endian::native == std::endian::little) {
                std::memcpy(output.data(), raw.data(), byte_count);
            } else {
                for (std::size_t i = 0; i < to_read; ++i) {
                    output[i] = detail::load_little_endian<float>(
                        raw.data() + (i * sizeof(float)));
                }
            }
        } else {
            // SIGPROC time-series
            if (m_raw_sigproc.nbits == 32) {
                const auto byte_offset = m_sample * 4U;
                const auto byte_count  = to_read * 4U;
                const auto raw         = mapped_bytes(byte_offset, byte_count);
                detail::unpack_32le(raw, output);
            } else if (m_raw_sigproc.nbits == 16) {
                const auto byte_offset = m_sample * 2U;
                const auto byte_count  = to_read * 2U;
                const auto raw         = mapped_bytes(byte_offset, byte_count);
                detail::unpack_16le(raw, output);
            } else if (m_raw_sigproc.nbits == 8) {
                const auto byte_offset = m_sample;
                const auto byte_count  = to_read;
                const auto raw         = mapped_bytes(byte_offset, byte_count);
                detail::unpack_8bit(raw, output,
                                    m_raw_sigproc.samples_are_signed());
            } else {
                const auto samples_per_byte =
                    8U / static_cast<std::uint32_t>(m_raw_sigproc.nbits);
                if (m_sample % samples_per_byte != 0U) {
                    throw ValidationError("psrio: sub-byte time series seek "
                                          "must be aligned to byte boundary");
                }
                const auto byte_offset = (m_sample * static_cast<std::uint32_t>(
                                                         m_raw_sigproc.nbits)) /
                                         8U;
                const auto byte_count  = ((to_read * static_cast<std::uint32_t>(
                                                         m_raw_sigproc.nbits)) +
                                          7U) /
                                         8U;
                const auto raw         = mapped_bytes(byte_offset, byte_count);
                detail::unpack_sub_byte(raw, output, m_raw_sigproc.nbits,
                                        order);
            }
        }

        m_sample += to_read;
        return to_read;
    }

    /// Read next @p count time samples into caller-provided @p dest using
    /// configured bit order.
    std::uint64_t read_samples(std::uint64_t count, std::span<float> dest) {
        return read_samples(count, dest, m_bit_order);
    }

    /// Alias for read_samples.
    std::uint64_t read(std::uint64_t count, std::span<float> dest) {
        return read_samples(count, dest, m_bit_order);
    }

    /// Convenience allocating read: read all data into an in-memory TimeSeries
    /// object.
    [[nodiscard]] TimeSeries read_data();

    /// Convenience allocating read: read next @p count samples into
    /// std::vector<float>.
    [[nodiscard]] std::vector<float> read_samples(std::uint64_t count) {
        std::vector<float> out(count);
        const auto actual = read_samples(count, std::span<float>(out));
        out.resize(static_cast<std::size_t>(actual));
        return out;
    }

private:
    void open_file(const std::filesystem::path& path) {
        const auto ext = path.extension().string();
        if (ext == ".dat") {
            const auto inf_path =
                path.string().substr(0, path.string().size() - 4) + ".inf";
            open_presto(path, inf_path);
        } else if (ext == ".inf") {
            const auto dat_path =
                path.string().substr(0, path.string().size() - 4) + ".dat";
            open_presto(dat_path, path);
        } else if (ext == ".tim") {
            open_sigproc(path);
        } else {
            // Check companion files or auto-detect
            const auto dat_cand = path.string() + ".dat";
            const auto inf_cand = path.string() + ".inf";
            if (std::filesystem::exists(dat_cand) &&
                std::filesystem::exists(inf_cand)) {
                open_presto(dat_cand, inf_cand);
            } else {
                try {
                    open_sigproc(path);
                } catch (...) {
                    open_presto(path, inf_cand);
                }
            }
        }
    }

    void open_presto(const std::filesystem::path& dat_path,
                     const std::filesystem::path& inf_path) {
        if (!std::filesystem::exists(inf_path)) {
            throw IoError(
                std::format("psrio: companion .inf file '{}' not found",
                            inf_path.string()));
        }
        if (!std::filesystem::exists(dat_path)) {
            throw IoError(
                std::format("psrio: companion .dat file '{}' not found",
                            dat_path.string()));
        }

        m_header          = Header::from_inffile(inf_path);
        m_header.filename = dat_path.string();
        m_is_presto       = true;

        try {
            m_file = detail::MappedFile(dat_path);
            m_file.advise_sequential();
        } catch (const std::exception& ex) {
            throw IoError(ex.what());
        }

        if (m_file.empty()) {
            throw FormatError(std::format(
                "psrio: PRESTO .dat file '{}' is empty", dat_path.string()));
        }

        const auto available_samples = m_file.size() / sizeof(float);
        if (m_header.nsamples == 0U) {
            m_header.nsamples = available_samples;
        } else if (m_header.nsamples > available_samples) {
            m_warnings.push_back(std::format(
                "psrio: .inf declared {} samples but file only has {}",
                m_header.nsamples, available_samples));
            m_header.nsamples = available_samples;
        }

        m_payload = m_file.bytes();
    }

    void open_sigproc(const std::filesystem::path& path) {
        try {
            m_file = detail::MappedFile(path);
            m_file.advise_sequential();
        } catch (const std::exception& ex) {
            throw IoError(ex.what());
        }

        const std::string path_str = path.string();
        m_raw_sigproc =
            formats::sigproc::FilterbankHeader::parse(m_file.bytes(), path_str);
        if (m_raw_sigproc.nchans != 1) {
            throw ValidationError(
                "psrio: SIGPROC time series requires nchans == 1");
        }

        const auto mapped = m_file.bytes();
        if (m_raw_sigproc.header_bytes > mapped.size()) {
            throw FormatError(
                std::format("{}: header extends past end of file", path_str));
        }
        m_payload = mapped.subspan(
            static_cast<std::size_t>(m_raw_sigproc.header_bytes));
        m_header          = Header::from_sigproc(path);
        m_header.filename = path_str;
        m_is_presto       = false;

        if (m_raw_sigproc.nbits <= 4) {
            m_header.nsamples = (m_payload.size() * 8U) /
                                static_cast<std::uint64_t>(m_raw_sigproc.nbits);
        } else {
            m_header.nsamples =
                m_payload.size() /
                (static_cast<std::size_t>(m_raw_sigproc.nbits) / 8U);
        }
    }

    [[nodiscard]] bool is_sub_byte() const noexcept {
        const int nbits = m_raw_sigproc.nbits;
        return !m_is_presto && (nbits == 1 || nbits == 2 || nbits == 4);
    }

    [[nodiscard]] std::uint64_t samples_per_byte() const {
        return 8U / static_cast<std::uint64_t>(m_raw_sigproc.nbits);
    }

    [[nodiscard]] std::uint64_t payload_offset(std::uint64_t sample) const {
        if (m_is_presto) {
            return sample * sizeof(float);
        }
        if (!is_sub_byte()) {
            return sample * bytes_per_sample();
        }
        return (sample * static_cast<std::uint64_t>(m_raw_sigproc.nbits)) / 8U;
    }

    [[nodiscard]] std::uint64_t payload_bytes(std::uint64_t count) const {
        if (!is_sub_byte()) {
            return count * bytes_per_sample();
        }
        return (count * static_cast<std::uint64_t>(m_raw_sigproc.nbits)) / 8U;
    }

    [[nodiscard]] std::uint64_t block_bytes(std::uint64_t count) const {
        if (is_sub_byte()) {
            if (count % samples_per_byte() != 0U) {
                throw ValidationError(
                    "psrio: sub-byte time series block read must cover whole "
                    "bytes");
            }
            return payload_bytes(count);
        }
        const auto stride = bytes_per_sample();
        if (count > 0U &&
            stride > (std::numeric_limits<std::uint64_t>::max() / count)) {
            throw ValidationError("psrio: requested sample block is too large");
        }
        return count * stride;
    }

    [[nodiscard]] std::span<const std::byte>
    mapped_bytes(std::uint64_t offset, std::uint64_t nbytes) const {
        if (offset > m_payload.size() || nbytes > m_payload.size() - offset) {
            throw FormatError("psrio: sample offset is outside mapped payload");
        }
        return m_payload.subspan(static_cast<std::size_t>(offset),
                                 static_cast<std::size_t>(nbytes));
    }

    Header m_header;
    detail::MappedFile m_file;
    std::span<const std::byte> m_payload;
    std::uint64_t m_sample{0};
    bool m_is_presto{false};
    formats::sigproc::FilterbankHeader m_raw_sigproc;
    BitOrder m_bit_order{BitOrder::kLsbFirst};
    std::vector<std::string> m_warnings;
};

/**
 * @brief In-memory 1D time-series container holding float data and
 * observational metadata.
 *
 * Modeled after sigpyproc's TimeSeries. Supports loading from and writing to
 * both SIGPROC (.tim) and PRESTO (.dat / .inf) formats.
 */
class TimeSeries {
public:
    TimeSeries() = default;

    TimeSeries(std::vector<float> data, Header header)
        : m_data(std::move(data)),
          m_header(std::move(header)) {
        m_header.nsamples = m_data.size();
    }

    /// Read-write access to time-series samples.
    [[nodiscard]] std::span<float> data() noexcept { return m_data; }

    /// Read-only access to time-series samples.
    [[nodiscard]] std::span<const float> data() const noexcept {
        return m_data;
    }

    /// Underlying observational metadata.
    [[nodiscard]] const Header& header() const noexcept { return m_header; }
    [[nodiscard]] Header& header() noexcept { return m_header; }

    /// Number of time samples.
    [[nodiscard]] std::uint64_t nsamples() const noexcept {
        return m_data.size();
    }

    /// Sampling interval in seconds.
    [[nodiscard]] double dt() const noexcept { return m_header.tsamp; }

    /// Observation duration in seconds.
    [[nodiscard]] double tobs() const noexcept {
        return static_cast<double>(nsamples()) * dt();
    }

    /// Load time series from a SIGPROC .tim file.
    [[nodiscard]] static TimeSeries
    from_tim(const std::filesystem::path& timfile) {
        TimeSeriesReader reader(timfile);
        return reader.read_data();
    }

    /// Load time series from a PRESTO .dat file (and optional companion .inf
    /// file).
    [[nodiscard]] static TimeSeries from_dat(
        const std::filesystem::path& datfile,
        const std::optional<std::filesystem::path>& inffile = std::nullopt) {
        if (inffile.has_value()) {
            TimeSeriesReader reader(datfile, *inffile);
            return reader.read_data();
        }
        TimeSeriesReader reader(datfile);
        return reader.read_data();
    }

    /// Load time series with auto-format detection.
    [[nodiscard]] static TimeSeries
    from_file(const std::filesystem::path& path) {
        TimeSeriesReader reader(path);
        return reader.read_data();
    }

    /// Open a lite streaming reader without loading data into memory.
    [[nodiscard]] static TimeSeriesReader
    open(const std::filesystem::path& path) {
        return TimeSeriesReader(path);
    }

    /// Write time series in SIGPROC .tim format.
    [[nodiscard]] std::string to_tim(const std::optional<std::filesystem::path>&
                                         filename = std::nullopt) const {
        std::filesystem::path out_path;
        if (filename.has_value()) {
            out_path = *filename;
        } else if (m_header.filename.empty()) {
            out_path = "timeseries.tim";
        } else {
            out_path =
                std::filesystem::path(m_header.filename).stem().string() +
                ".tim";
        }

        std::ofstream out(out_path, std::ios::binary);
        if (!out.is_open()) {
            throw IoError(std::format("psrio: cannot create .tim file '{}'",
                                      out_path.string()));
        }

        const auto write_str = [&](std::string_view s) {
            const auto len    = static_cast<std::uint32_t>(s.size());
            const auto le_len = detail::to_little_endian(len);
            out.write(reinterpret_cast<const char*>(&le_len), sizeof(le_len));
            out.write(s.data(), static_cast<std::streamsize>(s.size()));
        };

        const auto write_i32 = [&](std::string_view key, std::int32_t val) {
            write_str(key);
            const auto le_val = detail::to_little_endian(val);
            out.write(reinterpret_cast<const char*>(&le_val), sizeof(le_val));
        };

        const auto write_f64 = [&](std::string_view key, double val) {
            write_str(key);
            const auto le_val = detail::to_little_endian(val);
            out.write(reinterpret_cast<const char*>(&le_val), sizeof(le_val));
        };

        write_str("HEADER_START");
        write_i32("data_type", 2); // Dedispersed time series
        write_i32("nchans", 1);
        write_i32("nifs", m_header.nifs > 0 ? m_header.nifs : 1);
        write_i32("nbits", 32);
        write_f64("tsamp", m_header.tsamp > 0.0 ? m_header.tsamp : 1.0);
        write_f64("tstart", m_header.tstart);
        write_f64("fch1", m_header.fch1);
        write_f64("foff", m_header.foff);
        write_f64("src_raj", m_header.raj);
        write_f64("src_dej", m_header.dej);
        if (m_header.dm > 0.0) {
            write_f64("refdm", m_header.dm);
        }
        if (!m_header.source.empty()) {
            write_str("source_name");
            write_str(m_header.source);
        }
        write_i32("barycentric", m_header.barycentric);
        write_str("HEADER_END");

        if constexpr (std::endian::native == std::endian::little) {
            out.write(
                reinterpret_cast<const char*>(m_data.data()),
                static_cast<std::streamsize>(m_data.size() * sizeof(float)));
        } else {
            for (const float val : m_data) {
                const auto le_val = detail::to_little_endian(val);
                out.write(reinterpret_cast<const char*>(&le_val),
                          sizeof(le_val));
            }
        }

        return out_path.string();
    }

    /// Write time series in PRESTO .dat format and generate companion .inf
    /// file.
    [[nodiscard]] std::string to_dat(const std::optional<std::filesystem::path>&
                                         basename = std::nullopt) const {
        std::filesystem::path base;
        if (basename.has_value()) {
            base = *basename;
        } else if (m_header.filename.empty()) {
            base = "timeseries";
        } else {
            base = std::filesystem::path(m_header.filename).stem();
        }

        const auto dat_path = base.string() + ".dat";
        const auto inf_path = base.string() + ".inf";

        std::ofstream dat_out(dat_path, std::ios::binary);
        if (!dat_out.is_open()) {
            throw IoError(
                std::format("psrio: cannot create .dat file '{}'", dat_path));
        }

        if constexpr (std::endian::native == std::endian::little) {
            dat_out.write(
                reinterpret_cast<const char*>(m_data.data()),
                static_cast<std::streamsize>(m_data.size() * sizeof(float)));
        } else {
            for (const float val : m_data) {
                const auto le_val = detail::to_little_endian(val);
                dat_out.write(reinterpret_cast<const char*>(&le_val),
                              sizeof(le_val));
            }
        }

        m_header.make_inf(inf_path);
        return dat_path;
    }

private:
    std::vector<float> m_data;
    Header m_header;
};

inline TimeSeries TimeSeriesReader::read_data() {
    std::vector<float> buffer(m_header.nsamples);
    read_data(buffer);
    return {std::move(buffer), m_header};
}

inline Header Header::from_file(const std::filesystem::path& path) {
    const auto ext = path.extension().string();
    if (ext == ".inf") {
        return from_inffile(path);
    }
    if (ext == ".dat") {
        const auto inf_path =
            path.string().substr(0, path.string().size() - 4) + ".inf";
        if (std::filesystem::exists(inf_path)) {
            return from_inffile(inf_path);
        }
    }
    return from_sigproc(path);
}

} // namespace psrio
