#pragma once

#include "psrio/astro.hpp"
#include "psrio/common/types.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/mmap.hpp"
#include "psrio/detail/packed_bits.hpp"
#include "psrio/detail/readahead.hpp"
#include "psrio/detail/rfile.hpp"
#include "psrio/detail/skip.hpp"
#include "psrio/formats/sigproc/header.hpp"
#include "psrio/packed.hpp"

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <exception>
#include <filesystem>
#include <fmt/format.h>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace psrio::formats::sigproc {

/// Sequential reader for one SIGPROC filterbank file.
///
/// The file is read with `pread`: only the header is read at open, and each
/// read copies its samples straight from the file, which streams large files
/// on networked file systems far faster than a memory mapping. view_bytes()
/// and view_block() map the file on first use. read_samples(count, span) or
/// read(count, span)
/// unpacks the next time samples into a caller-owned span. read_bytes(nbytes,
/// span) and view_bytes(nbytes) expose the raw payload. One time sample is nifs
/// * nchans values in file order, index ipol * nchans + ichan, with no scaling
/// and no Stokes offset. Sub-byte samples are unpacked according to the
/// configured BitOrder.
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

    /// Number of frequency channels.
    [[nodiscard]] std::uint64_t nchans() const noexcept {
        return static_cast<std::uint64_t>(m_header.nchans);
    }

    /// Number of IFs stored in each time sample.
    [[nodiscard]] std::uint64_t nifs() const noexcept {
        return static_cast<std::uint64_t>(m_header.nifs);
    }

    /// Bits per on-disk sample.
    [[nodiscard]] int nbits() const noexcept { return m_header.nbits; }

    /// On-disk sample type. 32-bit SIGPROC samples are float.
    [[nodiscard]] SampleType sample_type() const {
        return sample_type_from_nbits(m_header.nbits,
                                      m_header.samples_are_signed());
    }

    /// Bytes in one time sample (`nifs * nchans * nbits / 8`).
    [[nodiscard]] std::uint64_t bytes_per_sample() const {
        return m_header.bytes_per_sample();
    }

    /// Complete time samples a read can return.
    [[nodiscard]] std::uint64_t nsamples() const { return m_header.nsamples(); }

    /// True. A filterbank file has a known length after the header is parsed.
    static bool has_nsamples() noexcept { return true; }

    /// Sample interval in seconds.
    [[nodiscard]] double tsamp() const noexcept {
        return m_header.tsamp.value_or(0.0);
    }

    /// MJD of the first sample.
    [[nodiscard]] double tstart() const noexcept {
        return m_header.tstart.value_or(0.0);
    }

    /// Centre frequency of channel 0, in MHz.
    [[nodiscard]] double fch1() const noexcept {
        return m_header.fch1.value_or(0.0);
    }

    /// Signed channel spacing, in MHz.
    [[nodiscard]] double foff() const noexcept {
        return m_header.foff.value_or(0.0);
    }

    /// Reference DM (`refdm`) the data were dedispersed to, when present.
    [[nodiscard]] std::optional<double> refdm() const noexcept {
        return m_header.refdm;
    }

    /// Beam index from the header, or 0 when `ibeam` is absent.
    [[nodiscard]] int beam() const noexcept {
        return static_cast<int>(m_header.ibeam.value_or(0));
    }

    /// Source name from header, or "Unknown".
    [[nodiscard]] std::string_view source_name() const noexcept {
        return m_header.source_name.has_value()
                   ? std::string_view(*m_header.source_name)
                   : std::string_view("Unknown");
    }

    /// Telescope name derived from telescope_id, or "Unknown".
    [[nodiscard]] std::string_view telescope() const noexcept {
        return m_header.telescope_name();
    }

    /// Right ascension packed sexagesimal (HHMMSS.ss), or 0.0.
    [[nodiscard]] double raj() const noexcept {
        return m_header.src_raj.value_or(0.0);
    }

    /// Declination packed sexagesimal (+-DDMMSS.ss), or 0.0.
    [[nodiscard]] double dej() const noexcept {
        return m_header.src_dej.value_or(0.0);
    }

    /// Total observing bandwidth across all channels, in MHz.
    [[nodiscard]] double bandwidth() const noexcept {
        return std::abs(foff()) * static_cast<double>(nchans());
    }

    /// Centre frequency of the full band, in MHz.
    [[nodiscard]] double center_frequency() const noexcept {
        return fch1() + (foff() * (static_cast<double>(nchans()) - 1.0) / 2.0);
    }

    /// Spectra per second, or 0 when `tsamp` is not positive.
    [[nodiscard]] double spectra_rate() const noexcept {
        const double interval = tsamp();
        return interval > 0.0 ? 1.0 / interval : 0.0;
    }

    /// UTC of `tstart`, as POSIX seconds.
    [[nodiscard]] std::time_t utc_start() const noexcept {
        return astro::mjd_to_time(tstart());
    }

    /// Enable or disable channel reversal (fswap) on read_block and
    /// read_samples. Raw byte reads and views stay in file order.
    void set_fswap(bool enable) noexcept { m_apply_fswap = enable; }

    /// Drop cached file pages once their samples are read. Off by default,
    /// or on when `PSRIO_DROP_CACHE=1` is set at open.
    void set_drop_cache(bool enable) noexcept { m_file.set_drop_cache(enable); }

    /// Once two reads in a row are contiguous, stream the rest of the file on
    /// a background thread, kReadAheadBlocks blocks ahead of the caller. On by
    /// default.
    void set_prefetch(bool enable) noexcept {
        m_prefetch = enable;
        if (!enable) {
            m_ahead.reset();
        }
    }

    /// Stop reading ahead and free its buffers, for example once a stream
    /// has moved past this file.
    void release_buffers() {
        m_ahead.reset();
        m_raw = {};
    }

    /// Blocks of ReadAhead::kBlockBytes read ahead of the caller.
    static constexpr std::size_t kReadAheadBlocks = 4;

    /// Most raw bytes read_samples() stages at once before unpacking.
    static constexpr std::uint64_t kStageBytes = std::uint64_t{16} << 20U;
    [[nodiscard]] bool fswap_enabled() const noexcept { return m_apply_fswap; }

    /// Move the cursor by @p delta time samples. Negative moves backward.
    /// @throws ValidationError if the landing index is outside
    ///         `[0, nsamples()]`.
    void skip(std::int64_t delta);

    /// Copy up to @p count packed time samples into @p dest.
    ///
    /// @p dest must hold at least `count * bytes_per_sample()` bytes. The copy
    /// is raw file order (or channel-reversed if fswap is enabled). Returns
    /// the number of time samples copied, which is short at the end of the
    /// file and zero when the cursor is already there.
    /// @throws ValidationError if @p dest is smaller than the requested block.
    std::uint64_t read_block(std::uint64_t count, std::span<std::byte> dest);

    /// Zero-copy view of the next @p count time samples.
    [[nodiscard]] std::span<const std::byte> view_block(std::uint64_t count);

    /// Convenience allocating read: copy next @p count packed time samples.
    [[nodiscard]] std::vector<std::byte> read_block(std::uint64_t count);

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
    std::uint64_t
    read_samples(std::uint64_t count, std::span<T> dest, BitOrder order);

    /// Unpack the next @p count time samples into @p dest using the reader's
    /// bit_order().
    template <typename T>
        requires std::same_as<T, float> || std::same_as<T, std::uint8_t> ||
                 std::same_as<T, std::uint16_t>
    std::uint64_t read_samples(std::uint64_t count, std::span<T> dest) {
        return read_samples(count, dest, m_bit_order);
    }

    /// Unpack alias: unpack the next @p count time samples into @p dest.
    template <typename T>
        requires std::same_as<T, float> || std::same_as<T, std::uint8_t> ||
                 std::same_as<T, std::uint16_t>
    std::uint64_t read(std::uint64_t count, std::span<T> dest, BitOrder order) {
        return read_samples(count, dest, order);
    }

    /// Unpack alias: unpack the next @p count time samples using reader's
    /// bit_order().
    template <typename T>
        requires std::same_as<T, float> || std::same_as<T, std::uint8_t> ||
                 std::same_as<T, std::uint16_t>
    std::uint64_t read(std::uint64_t count, std::span<T> dest) {
        return read_samples(count, dest, m_bit_order);
    }

    /// Copy the next @p nbytes raw payload bytes into @p dest.
    /// @throws ValidationError if the request is not a positive multiple of
    ///         the sample stride, does not fit in the readable samples, or
    ///         @p dest has a different size.
    std::uint64_t read_bytes(std::uint64_t nbytes, std::span<std::byte> dest);

    /// Zero-copy view of the next @p nbytes raw payload bytes.
    /// The first view maps the file. The span refers to that mapping and
    /// stays valid until the reader is destroyed. Advances the cursor by the
    /// corresponding number of samples.
    [[nodiscard]] std::span<const std::byte> view_bytes(std::uint64_t nbytes);

    /// Zero-copy view alias for view_bytes.
    [[nodiscard]] std::span<const std::byte> view(std::uint64_t nbytes) {
        return view_bytes(nbytes);
    }

    /// Convenience allocating read: unpack next @p count time samples into a
    /// new std::vector<T>.
    template <typename T = float>
        requires std::same_as<T, float> || std::same_as<T, std::uint8_t> ||
                 std::same_as<T, std::uint16_t>
    [[nodiscard]] std::vector<T> read_samples(std::uint64_t count,
                                              BitOrder order) {
        const auto per_sample = values_per_sample();
        std::vector<T> buffer(count * per_sample);
        const auto actual = read_samples(count, std::span<T>(buffer), order);
        buffer.resize(static_cast<std::size_t>(actual * per_sample));
        return buffer;
    }

    /// Convenience allocating read using the reader's bit_order().
    template <typename T = float>
        requires std::same_as<T, float> || std::same_as<T, std::uint8_t> ||
                 std::same_as<T, std::uint16_t>
    [[nodiscard]] std::vector<T> read_samples(std::uint64_t count) {
        return read_samples<T>(count, m_bit_order);
    }

    /// Convenience allocating read for raw bytes.
    [[nodiscard]] std::vector<std::byte> read_bytes(std::uint64_t nbytes) {
        std::vector<std::byte> buffer(nbytes);
        const auto actual = read_bytes(nbytes, std::span<std::byte>(buffer));
        buffer.resize(static_cast<std::size_t>(actual));
        return buffer;
    }

private:
    [[nodiscard]] std::uint64_t values_per_sample() const;

    void check_byte_request(std::uint64_t nbytes) const;

    [[nodiscard]] std::span<const std::byte> mapped_bytes(std::uint64_t offset,
                                                          std::uint64_t nbytes);

    /// Copy payload bytes [offset, offset + dest.size()) into @p dest.
    void read_payload(std::uint64_t offset, std::span<std::byte> dest);

    ::psrio::detail::ReadFile m_file;
    ::psrio::detail::MappedFile m_map;
    FilterbankHeader m_header;
    std::vector<std::byte> m_raw;
    std::unique_ptr<::psrio::detail::ReadAhead> m_ahead;
    std::uint64_t m_read_end{std::numeric_limits<std::uint64_t>::max()};
    bool m_prefetch{true};
    std::uint64_t m_sample{0};
    BitOrder m_bit_order{BitOrder::kLsbFirst};
    bool m_apply_fswap{false};
};

inline FilterbankReader::FilterbankReader(const std::filesystem::path& path,
                                          BitOrder bit_order)
    : m_bit_order(bit_order) {
    try {
        m_file = ::psrio::detail::ReadFile(path);
    } catch (const std::exception& ex) {
        throw IoError(ex.what());
    }

    const std::string path_string = path.string();
    if (m_file.size() == 0U) {
        throw FormatError(fmt::format("{}: file is empty", path_string));
    }
    m_header = FilterbankHeader::read(m_file);
    try {
        m_header.validate();
    } catch (const ValidationError& ex) {
        throw ValidationError(fmt::format("{}: {}", path_string, ex.what()));
    }

    if (m_header.header_bytes > m_file.size()) {
        throw FormatError(
            fmt::format("{}: header extends past end of file", path_string));
    }
}

inline void FilterbankReader::seek(std::uint64_t sample) {
    if (sample > m_header.nsamples()) {
        throw ValidationError("psrio: seek is past the readable samples");
    }
    m_sample = sample;
}

inline void FilterbankReader::skip(std::int64_t delta) {
    m_sample =
        ::psrio::detail::apply_skip(m_sample, delta, m_header.nsamples());
}

inline std::uint64_t FilterbankReader::read_block(std::uint64_t count,
                                                  std::span<std::byte> dest) {
    const auto stride = m_header.bytes_per_sample();
    if (count > 0U &&
        stride > (std::numeric_limits<std::uint64_t>::max() / count)) {
        throw ValidationError("psrio: requested sample block is too large");
    }
    const auto bytes_needed = count * stride;
    if (dest.size() < bytes_needed) {
        throw ValidationError(fmt::format(
            "psrio: destination has {} bytes but the block needs {}",
            dest.size(), bytes_needed));
    }
    if (m_sample > m_header.nsamples()) {
        throw ValidationError(
            "psrio: reader cursor is past the readable samples");
    }
    const auto to_read = std::min(count, m_header.nsamples() - m_sample);
    if (to_read == 0U) {
        return 0;
    }
    const auto nbytes = to_read * stride;
    read_payload(m_sample * stride,
                 dest.first(static_cast<std::size_t>(nbytes)));
    if (m_apply_fswap && to_read > 0U) {
        ::psrio::reverse_channels(dest.first(static_cast<std::size_t>(nbytes)),
                                  to_read * nifs(), nchans(), nbits());
    }
    m_sample += to_read;
    return to_read;
}

inline std::span<const std::byte>
FilterbankReader::view_block(std::uint64_t count) {
    const auto stride = m_header.bytes_per_sample();
    if (count > 0U &&
        stride > (std::numeric_limits<std::uint64_t>::max() / count)) {
        throw ValidationError("psrio: requested sample block is too large");
    }
    if (m_sample > m_header.nsamples() ||
        count > m_header.nsamples() - m_sample) {
        throw ValidationError(
            "psrio: view_block extends past readable samples");
    }
    const auto nbytes = count * stride;
    const auto raw    = mapped_bytes(m_sample * stride, nbytes);
    m_sample += count;
    return raw;
}

inline std::vector<std::byte>
FilterbankReader::read_block(std::uint64_t count) {
    const auto stride = m_header.bytes_per_sample();
    std::vector<std::byte> out(static_cast<std::size_t>(count * stride));
    const auto actual = read_block(count, std::span<std::byte>(out));
    out.resize(static_cast<std::size_t>(actual * stride));
    return out;
}

inline std::uint64_t FilterbankReader::values_per_sample() const {
    return static_cast<std::uint64_t>(m_header.nifs) *
           static_cast<std::uint64_t>(m_header.nchans);
}

inline std::span<const std::byte>
FilterbankReader::mapped_bytes(std::uint64_t offset, std::uint64_t nbytes) {
    if (!m_map.is_mapped()) {
        try {
            m_map = ::psrio::detail::MappedFile(m_file.path());
            m_map.advise_sequential();
        } catch (const std::exception& ex) {
            throw IoError(ex.what());
        }
    }
    const auto mapped = m_map.bytes();
    const auto start  = m_header.header_bytes;
    if (start > mapped.size()) {
        throw FormatError("psrio: header extends past the mapped file");
    }
    const auto payload = mapped.subspan(static_cast<std::size_t>(start));
    if (offset > payload.size() || nbytes > payload.size() - offset) {
        throw FormatError("psrio: sample offset is outside the mapped payload");
    }
    return payload.subspan(static_cast<std::size_t>(offset),
                           static_cast<std::size_t>(nbytes));
}

inline void FilterbankReader::read_payload(std::uint64_t offset,
                                           std::span<std::byte> dest) {
    const auto start = m_header.header_bytes;
    const auto size  = m_file.size() - start;
    if (offset > size || dest.size() > size - offset) {
        throw FormatError("psrio: sample offset is outside the payload");
    }
    const auto at = start + offset;
    if (m_ahead != nullptr && !m_ahead->covers(at, dest.size())) {
        m_ahead.reset();
    }
    if (m_ahead == nullptr && m_prefetch && at == m_read_end && !dest.empty()) {
        m_ahead = std::make_unique<::psrio::detail::ReadAhead>(
            m_file.duplicate(), at, m_file.size(), kReadAheadBlocks);
    }
    if (m_ahead != nullptr) {
        m_ahead->read(at, dest);
    } else {
        m_file.read_at(at, dest);
        m_file.release(at, dest.size());
    }
    m_read_end = at + dest.size();
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
inline std::uint64_t FilterbankReader::read_samples(std::uint64_t count,
                                                    std::span<T> dest,
                                                    BitOrder order) {
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
    if (count > 0U &&
        per_sample > (std::numeric_limits<std::uint64_t>::max() / count)) {
        throw ValidationError("psrio: requested sample block is too large");
    }
    const auto expected = count * per_sample;
    if (expected != dest.size()) {
        throw ValidationError(fmt::format(
            "psrio: destination has {} values but the request needs {}",
            dest.size(), expected));
    }
    if (m_sample > m_header.nsamples()) {
        throw ValidationError(
            "psrio: reader cursor is past the readable samples");
    }

    const auto to_read = std::min(count, m_header.nsamples() - m_sample);
    if (to_read == 0U) {
        return 0;
    }

    const auto stride = m_header.bytes_per_sample();
    const auto output =
        dest.first(static_cast<std::size_t>(to_read * per_sample));
    // Raw bytes are staged at most kStageBytes at a time, so the reader's
    // own buffer stays bounded whatever the caller asks for.
    const auto chunk = std::max<std::uint64_t>(kStageBytes / stride, 1U);
    for (std::uint64_t done = 0; done < to_read;) {
        const auto rows = std::min(chunk, to_read - done);
        m_raw.resize(static_cast<std::size_t>(rows * stride));
        read_payload((m_sample + done) * stride, m_raw);
        const std::span<const std::byte> raw(m_raw);
        const auto out =
            output.subspan(static_cast<std::size_t>(done * per_sample),
                           static_cast<std::size_t>(rows * per_sample));
        if (m_header.nbits <= 4) {
            if constexpr (std::same_as<T, std::uint16_t>) {
                throw ValidationError(
                    "psrio: uint16_t output requires 16-bit samples");
            } else {
                ::psrio::detail::unpack_sub_byte(raw, out, m_header.nbits,
                                                 order);
            }
        } else if (m_header.nbits == 8) {
            if constexpr (std::same_as<T, std::uint16_t>) {
                throw ValidationError(
                    "psrio: uint16_t output requires 16-bit samples");
            } else {
                ::psrio::detail::unpack_8bit(raw, out,
                                             m_header.samples_are_signed());
            }
        } else if (m_header.nbits == 16) {
            if constexpr (std::same_as<T, std::uint8_t>) {
                throw ValidationError(
                    "psrio: uint8_t output requires nbits of 8 or less");
            } else {
                ::psrio::detail::unpack_16le(raw, out);
            }
        } else if constexpr (std::same_as<T, float>) {
            ::psrio::detail::unpack_32le(raw, out);
        } else {
            throw ValidationError("psrio: 32-bit samples unpack to float");
        }
        done += rows;
    }

    if (m_apply_fswap) {
        ::psrio::reverse_unpacked_channels(output, to_read, nchans(), nifs());
    }
    m_sample += to_read;
    return to_read;
}

inline std::uint64_t FilterbankReader::read_bytes(std::uint64_t nbytes,
                                                  std::span<std::byte> dest) {
    check_byte_request(nbytes);
    if (dest.size() != nbytes) {
        throw ValidationError(fmt::format(
            "psrio: destination has {} bytes but the request needs {}",
            dest.size(), nbytes));
    }
    const auto stride = m_header.bytes_per_sample();
    read_payload(m_sample * stride, dest);
    m_sample += nbytes / stride;
    return nbytes;
}

inline std::span<const std::byte>
FilterbankReader::view_bytes(std::uint64_t nbytes) {
    check_byte_request(nbytes);
    const auto stride = m_header.bytes_per_sample();
    const auto raw    = mapped_bytes(m_sample * stride, nbytes);
    m_sample += nbytes / stride;
    return raw;
}

} // namespace psrio::formats::sigproc
