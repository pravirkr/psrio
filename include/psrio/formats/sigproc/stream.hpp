#pragma once

#include "psrio/astro.hpp"
#include "psrio/common/types.hpp"
#include "psrio/detail/contiguity.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/skip.hpp"
#include "psrio/formats/sigproc/header.hpp"
#include "psrio/formats/sigproc/reader.hpp"

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <fmt/format.h>
#include <limits>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

namespace psrio::formats::sigproc {

namespace stream_detail {

template <typename T>
void match_value(const std::optional<T>& left,
                 const std::optional<T>& right,
                 std::string_view key,
                 std::string_view path) {
    if (left != right) {
        throw FormatError(fmt::format(
            "psrio: {} {} does not match the first file", path, key));
    }
}

inline void match_header(const FilterbankHeader& first,
                         const FilterbankHeader& next,
                         std::string_view path) {
    if (next.nchans != first.nchans || next.nifs != first.nifs ||
        next.nbits != first.nbits) {
        throw FormatError(fmt::format(
            "psrio: {} does not match the first file's channel layout", path));
    }
    match_value(next.tsamp, first.tsamp, "tsamp", path);
    match_value(next.fch1, first.fch1, "fch1", path);
    match_value(next.foff, first.foff, "foff", path);
    if (next.frequency_table != first.frequency_table) {
        throw FormatError(fmt::format(
            "psrio: {} frequency table does not match the first file", path));
    }
    match_value(next.data_type, first.data_type, "data_type", path);
    match_value(next.telescope_id, first.telescope_id, "telescope_id", path);
    match_value(next.machine_id, first.machine_id, "machine_id", path);
    match_value(next.source_name, first.source_name, "source_name", path);
    match_value(next.src_raj, first.src_raj, "src_raj", path);
    match_value(next.src_dej, first.src_dej, "src_dej", path);
    match_value(next.az_start, first.az_start, "az_start", path);
    match_value(next.za_start, first.za_start, "za_start", path);
    match_value(next.barycentric, first.barycentric, "barycentric", path);
    match_value(next.pulsarcentric, first.pulsarcentric, "pulsarcentric", path);
    match_value(next.ibeam, first.ibeam, "ibeam", path);
    match_value(next.nbeams, first.nbeams, "nbeams", path);
    match_value(next.refdm, first.refdm, "refdm", path);
    match_value(next.period, first.period, "period", path);
    match_value(next.nbins, first.nbins, "nbins", path);
    match_value(next.npuls, first.npuls, "npuls", path);
    match_value(next.signed_data, first.signed_data, "signed", path);
}

} // namespace stream_detail

/**
 * @brief One time-ordered sequence of SIGPROC filterbank files.
 *
 * Files stay in the order given. Each file after the first must match the
 * first file's layout, and must start one sample after the previous file
 * ends. `nsamples()` is the sum. The header is the first file's header with
 * that combined length.
 */
class FilterbankStream {
public:
    explicit FilterbankStream(std::span<const std::filesystem::path> paths,
                              BitOrder bit_order = BitOrder::kLsbFirst);

    FilterbankStream(const FilterbankStream&)                = delete;
    FilterbankStream& operator=(const FilterbankStream&)     = delete;
    FilterbankStream(FilterbankStream&&) noexcept            = default;
    FilterbankStream& operator=(FilterbankStream&&) noexcept = default;
    ~FilterbankStream()                                      = default;

    [[nodiscard]] const FilterbankHeader& header() const noexcept {
        return m_header;
    }
    [[nodiscard]] BitOrder bit_order() const noexcept { return m_bit_order; }
    void set_bit_order(BitOrder order) noexcept {
        m_bit_order = order;
        for (auto& file : m_files) {
            file.reader.set_bit_order(order);
        }
    }

    [[nodiscard]] std::uint64_t tell() const noexcept { return m_sample; }
    void rewind() { seek(0); }
    [[nodiscard]] std::uint64_t nchans() const noexcept {
        return static_cast<std::uint64_t>(m_header.nchans);
    }
    [[nodiscard]] std::uint64_t nifs() const noexcept {
        return static_cast<std::uint64_t>(m_header.nifs);
    }
    [[nodiscard]] int nbits() const noexcept { return m_header.nbits; }
    [[nodiscard]] SampleType sample_type() const {
        return sample_type_from_nbits(m_header.nbits,
                                      m_header.samples_are_signed());
    }
    [[nodiscard]] std::uint64_t bytes_per_sample() const {
        return m_header.bytes_per_sample();
    }
    [[nodiscard]] std::uint64_t nsamples() const noexcept { return m_nsamples; }
    [[nodiscard]] static bool has_nsamples() noexcept { return true; }
    [[nodiscard]] double tsamp() const noexcept {
        return m_header.tsamp.value_or(0.0);
    }
    [[nodiscard]] double tstart() const noexcept {
        return m_header.tstart.value_or(0.0);
    }
    [[nodiscard]] std::string_view source_name() const noexcept {
        return m_header.source_name.has_value()
                   ? std::string_view(*m_header.source_name)
                   : std::string_view("Unknown");
    }
    [[nodiscard]] std::string_view telescope() const noexcept {
        return m_header.telescope_name();
    }
    [[nodiscard]] double raj() const noexcept {
        return m_header.src_raj.value_or(0.0);
    }
    [[nodiscard]] double dej() const noexcept {
        return m_header.src_dej.value_or(0.0);
    }
    [[nodiscard]] double fch1() const noexcept {
        return m_header.fch1.value_or(0.0);
    }
    [[nodiscard]] double foff() const noexcept {
        return m_header.foff.value_or(0.0);
    }
    [[nodiscard]] double bandwidth() const { return m_header.bandwidth(); }
    [[nodiscard]] double center_frequency() const { return m_header.fcenter(); }
    [[nodiscard]] int beam() const noexcept {
        return static_cast<int>(m_header.ibeam.value_or(0));
    }
    [[nodiscard]] double spectra_rate() const noexcept {
        const double interval = tsamp();
        return interval > 0.0 ? 1.0 / interval : 0.0;
    }
    [[nodiscard]] std::time_t utc_start() const noexcept {
        return astro::mjd_to_time(tstart());
    }

    void set_fswap(bool enable) noexcept {
        m_fswap = enable;
        for (auto& file : m_files) {
            file.reader.set_fswap(enable);
        }
    }
    [[nodiscard]] bool fswap_enabled() const noexcept { return m_fswap; }

    void seek(std::uint64_t sample);
    void skip(std::int64_t delta) {
        m_sample = ::psrio::detail::apply_skip(m_sample, delta, m_nsamples);
    }

    std::uint64_t read_block(std::uint64_t count, std::span<std::byte> dest);
    std::uint64_t read_samples(std::uint64_t count, std::span<float> dest);
    std::uint64_t read_samples(std::uint64_t count,
                               std::span<std::uint8_t> dest);
    std::uint64_t read_samples(std::uint64_t count,
                               std::span<std::uint16_t> dest);
    std::uint64_t read_bytes(std::uint64_t nbytes, std::span<std::byte> dest);

    template <typename T>
        requires std::same_as<T, float> || std::same_as<T, std::uint8_t> ||
                 std::same_as<T, std::uint16_t>
    std::uint64_t read(std::uint64_t count, std::span<T> dest) {
        return read_samples(count, dest);
    }

    [[nodiscard]] std::vector<std::byte> read_block(std::uint64_t count) {
        const auto stride = bytes_per_sample();
        std::vector<std::byte> out(static_cast<std::size_t>(count * stride));
        const auto actual = read_block(count, out);
        out.resize(static_cast<std::size_t>(actual * stride));
        return out;
    }

    [[nodiscard]] std::vector<std::byte> read_bytes(std::uint64_t nbytes) {
        std::vector<std::byte> out(static_cast<std::size_t>(nbytes));
        const auto actual = read_bytes(nbytes, out);
        out.resize(static_cast<std::size_t>(actual));
        return out;
    }

    template <typename T = float>
        requires std::same_as<T, float> || std::same_as<T, std::uint8_t> ||
                 std::same_as<T, std::uint16_t>
    [[nodiscard]] std::vector<T> read_samples(std::uint64_t count) {
        const auto per_sample = nifs() * nchans();
        std::vector<T> buffer(static_cast<std::size_t>(count * per_sample));
        const auto actual = read_samples(count, std::span<T>(buffer));
        buffer.resize(static_cast<std::size_t>(actual * per_sample));
        return buffer;
    }

private:
    struct File {
        FilterbankReader reader;
        std::uint64_t first{0};
        std::uint64_t count{0};
    };

    [[nodiscard]] std::size_t index_at(std::uint64_t sample) const;

    template <typename T, typename Read>
    std::uint64_t read_across(std::uint64_t count,
                              std::span<T> dest,
                              std::uint64_t per_sample,
                              const Read& read);

    std::vector<File> m_files;
    FilterbankHeader m_header{};
    std::uint64_t m_nsamples{0};
    std::uint64_t m_sample{0};
    BitOrder m_bit_order{BitOrder::kLsbFirst};
    bool m_fswap{false};
};

inline FilterbankStream::FilterbankStream(
    std::span<const std::filesystem::path> paths, BitOrder bit_order)
    : m_bit_order(bit_order) {
    if (paths.empty()) {
        throw ValidationError("psrio: a filterbank stream needs a file");
    }
    m_files.reserve(paths.size());
    for (const auto& path : paths) {
        m_files.push_back(File{
            .reader = FilterbankReader(path, bit_order),
            .first  = 0,
            .count  = 0,
        });
    }
    m_header = m_files.front().reader.header();
    if (paths.size() > 1U &&
        (!m_header.tstart.has_value() || !m_header.tsamp.has_value())) {
        throw FormatError(
            "psrio: contiguous filterbanks need tstart and tsamp");
    }
    std::uint64_t cursor = 0;
    for (std::size_t index = 0; index < m_files.size(); ++index) {
        auto& file = m_files[index];
        file.first = cursor;
        file.count = file.reader.nsamples();
        if (paths.size() > 1U && file.count == 0U) {
            throw FormatError(
                fmt::format("psrio: {} has no samples to continue a stream",
                            paths[index].string()));
        }
        if (index > 0U) {
            const auto& header = file.reader.header();
            stream_detail::match_header(m_header, header,
                                        paths[index].string());
            if (!header.tstart.has_value() || !header.tsamp.has_value()) {
                throw FormatError(fmt::format(
                    "psrio: {} needs tstart and tsamp to continue a stream",
                    paths[index].string()));
            }
            ::psrio::detail::require_sample_gap(
                m_files[index - 1U].reader.tstart(), m_files[index - 1U].count,
                m_files[index - 1U].reader.tsamp(), *header.tstart,
                paths[index].string());
        }
        cursor += file.count;
    }
    m_nsamples = cursor;
    m_header.declared_nsamples.reset();
    const auto stride   = m_header.bytes_per_sample();
    m_header.data_bytes = m_nsamples * stride;
    m_header.file_bytes = m_header.header_bytes + m_header.data_bytes;
}

inline void FilterbankStream::seek(std::uint64_t sample) {
    if (sample > m_nsamples) {
        throw ValidationError("psrio: seek is past the readable samples");
    }
    m_sample = sample;
}

inline std::size_t FilterbankStream::index_at(std::uint64_t sample) const {
    for (std::size_t index = 0; index < m_files.size(); ++index) {
        const auto& file = m_files[index];
        if (sample < file.first + file.count) {
            return index;
        }
    }
    return m_files.size() - 1U;
}

template <typename T, typename Read>
inline std::uint64_t FilterbankStream::read_across(std::uint64_t count,
                                                   std::span<T> dest,
                                                   std::uint64_t per_sample,
                                                   const Read& read) {
    if (count > 0U &&
        per_sample > (std::numeric_limits<std::uint64_t>::max() / count)) {
        throw ValidationError("psrio: requested sample block is too large");
    }
    const auto need = count * per_sample;
    if (need > std::numeric_limits<std::size_t>::max() ||
        dest.size() < static_cast<std::size_t>(need)) {
        throw ValidationError(fmt::format(
            "psrio: destination has {} values but the block needs {}",
            dest.size(), need));
    }
    if (m_sample > m_nsamples) {
        throw ValidationError(
            "psrio: reader cursor is past the readable samples");
    }
    const auto to_read = std::min(count, m_nsamples - m_sample);
    std::uint64_t done = 0;
    while (done < to_read) {
        auto& file       = m_files[index_at(m_sample + done)];
        const auto local = (m_sample + done) - file.first;
        const auto room  = file.count - local;
        const auto take  = std::min(room, to_read - done);
        file.reader.seek(local);
        const auto piece =
            dest.subspan(static_cast<std::size_t>(done * per_sample),
                         static_cast<std::size_t>(take * per_sample));
        const auto got = read(file.reader, take, piece);
        if (got == 0U) {
            break;
        }
        done += got;
    }
    m_sample += done;
    return done;
}

inline std::uint64_t FilterbankStream::read_block(std::uint64_t count,
                                                  std::span<std::byte> dest) {
    return read_across(count, dest, bytes_per_sample(),
                       [](FilterbankReader& reader, std::uint64_t take,
                          std::span<std::byte> piece) {
                           return reader.read_block(take, piece);
                       });
}

inline std::uint64_t FilterbankStream::read_samples(std::uint64_t count,
                                                    std::span<float> dest) {
    return read_across(count, dest, nifs() * nchans(),
                       [](FilterbankReader& reader, std::uint64_t take,
                          std::span<float> piece) {
                           return reader.read_samples(take, piece);
                       });
}

inline std::uint64_t
FilterbankStream::read_samples(std::uint64_t count,
                               std::span<std::uint8_t> dest) {
    return read_across(count, dest, nifs() * nchans(),
                       [](FilterbankReader& reader, std::uint64_t take,
                          std::span<std::uint8_t> piece) {
                           return reader.read_samples(take, piece);
                       });
}

inline std::uint64_t
FilterbankStream::read_samples(std::uint64_t count,
                               std::span<std::uint16_t> dest) {
    return read_across(count, dest, nifs() * nchans(),
                       [](FilterbankReader& reader, std::uint64_t take,
                          std::span<std::uint16_t> piece) {
                           return reader.read_samples(take, piece);
                       });
}

inline std::uint64_t FilterbankStream::read_bytes(std::uint64_t nbytes,
                                                  std::span<std::byte> dest) {
    const auto stride = bytes_per_sample();
    if (nbytes == 0U || stride == 0U || nbytes % stride != 0U) {
        throw ValidationError(
            "psrio: byte request must be a positive multiple of the sample "
            "stride");
    }
    if (dest.size() != nbytes) {
        throw ValidationError(fmt::format(
            "psrio: destination has {} bytes but the request needs {}",
            dest.size(), nbytes));
    }
    const auto samples = nbytes / stride;
    if (m_sample > m_nsamples || samples > m_nsamples - m_sample) {
        throw ValidationError(
            "psrio: byte request extends past the readable samples");
    }
    return read_across(samples, dest, stride,
                       [](FilterbankReader& reader, std::uint64_t take,
                          std::span<std::byte> piece) {
                           return reader.read_bytes(
                                      take * reader.bytes_per_sample(), piece) /
                                  reader.bytes_per_sample();
                       }) *
           stride;
}

} // namespace psrio::formats::sigproc
