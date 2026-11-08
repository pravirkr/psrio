#pragma once

#include "psrio/detail/contiguity.hpp"
#include "psrio/detail/skip.hpp"
#include "psrio/formats/psrfits/reader.hpp"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fmt/format.h>
#include <limits>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace psrio::formats::psrfits {

namespace stream_detail {

inline void match_search(const SearchFile& first,
                         const SearchFile& next,
                         std::string_view path) {
    const auto fail = [&](std::string_view what) {
        throw FormatError(fmt::format(
            "psrio: {} {} does not match the first file", path, what));
    };
    if (next.nchans() != first.nchans() ||
        next.input_nbits() != first.input_nbits() ||
        next.input_npol() != first.input_npol() ||
        next.pol_type() != first.pol_type() ||
        next.sums_polarisations() != first.sums_polarisations() ||
        next.subint_samples() != first.subint_samples()) {
        fail("channel layout");
    }
    if (next.tsamp() != first.tsamp() || next.fch1() != first.fch1() ||
        next.foff() != first.foff() ||
        next.source_name() != first.source_name() ||
        next.telescope() != first.telescope() || next.raj() != first.raj() ||
        next.dej() != first.dej() || next.beam() != first.beam()) {
        fail("observing setup");
    }
    const auto left  = first.sigproc_header();
    const auto right = next.sigproc_header();
    if (left.frequency_table != right.frequency_table ||
        left.refdm != right.refdm) {
        fail("frequency table");
    }
}

inline void require_nsuboffs(const SearchFile& previous,
                             const SearchFile& next,
                             std::string_view path) {
    const auto before = previous.nsuboffs();
    const auto after  = next.nsuboffs();
    if (before.has_value() != after.has_value()) {
        throw FormatError(fmt::format(
            "psrio: {} must carry NSUBOFFS when the other files do", path));
    }
    if (!before.has_value()) {
        return;
    }
    if (previous.nrows() > static_cast<std::uint64_t>(
                               std::numeric_limits<std::int64_t>::max()) ||
        *before > std::numeric_limits<std::int64_t>::max() -
                      static_cast<std::int64_t>(previous.nrows())) {
        throw FormatError(fmt::format(
            "psrio: NSUBOFFS of {} does not fit in a signed integer", path));
    }
    const auto expected = *before + static_cast<std::int64_t>(previous.nrows());
    if (*after != expected) {
        throw FormatError(fmt::format("psrio: {} NSUBOFFS is {}, expected {}",
                                      path, *after, expected));
    }
}

} // namespace stream_detail

/**
 * @brief One time-ordered sequence of search-mode PSRFITS files.
 *
 * Files stay in the order given. Each file after the first must match the
 * first file's layout and start one sample after the previous file ends.
 * When `NSUBOFFS` is present it must advance by the previous file's subint
 * count. The calibrated stream is float32 Stokes I, as in `SearchFile`.
 */
class SearchStream {
public:
    explicit SearchStream(std::span<const std::filesystem::path> paths);

    SearchStream(const SearchStream&)                = delete;
    SearchStream& operator=(const SearchStream&)     = delete;
    SearchStream(SearchStream&&) noexcept            = default;
    SearchStream& operator=(SearchStream&&) noexcept = default;
    ~SearchStream()                                  = default;

    [[nodiscard]] const sigproc::FilterbankHeader& sigproc_header() const {
        return m_header;
    }
    [[nodiscard]] int input_nbits() const noexcept { return m_input_nbits; }
    [[nodiscard]] std::uint64_t input_npol() const noexcept { return m_npol; }
    [[nodiscard]] std::uint64_t subint_samples() const noexcept {
        return m_nsblk;
    }
    [[nodiscard]] bool sums_polarisations() const noexcept { return m_sum; }
    [[nodiscard]] std::string_view pol_type() const noexcept {
        return m_pol_type;
    }

    [[nodiscard]] std::uint64_t nchans() const noexcept { return m_nchan; }
    [[nodiscard]] std::uint64_t nifs() const noexcept { return 1; }
    [[nodiscard]] int nbits() const noexcept { return 32; }
    [[nodiscard]] SampleType sample_type() const noexcept {
        return SampleType::kFloat32;
    }
    [[nodiscard]] std::uint64_t bytes_per_sample() const noexcept {
        return m_bytes_per_sample;
    }
    [[nodiscard]] std::uint64_t nsamples() const noexcept { return m_nsamples; }
    [[nodiscard]] bool has_nsamples() const noexcept { return true; }
    [[nodiscard]] double tsamp() const noexcept { return m_tsamp; }
    [[nodiscard]] double tstart() const noexcept { return m_tstart; }
    [[nodiscard]] std::string_view source_name() const noexcept {
        return m_source;
    }
    [[nodiscard]] std::string_view telescope() const noexcept {
        return m_telescope;
    }
    [[nodiscard]] double raj() const noexcept { return m_raj; }
    [[nodiscard]] double dej() const noexcept { return m_dej; }
    [[nodiscard]] double fch1() const noexcept { return m_fch1; }
    [[nodiscard]] double foff() const noexcept { return m_foff; }
    [[nodiscard]] double bandwidth() const noexcept { return m_bandwidth; }
    [[nodiscard]] double center_frequency() const noexcept { return m_centre; }
    [[nodiscard]] int beam() const noexcept { return m_beam; }
    [[nodiscard]] double spectra_rate() const noexcept {
        return m_tsamp > 0.0 ? 1.0 / m_tsamp : 0.0;
    }
    [[nodiscard]] std::time_t utc_start() const noexcept {
        return astro::mjd_to_time(m_tstart);
    }
    [[nodiscard]] std::uint64_t tell() const noexcept { return m_sample; }

    void set_fswap(bool enable) noexcept {
        m_fswap = enable;
        for (auto& file : m_files) {
            file.reader.set_fswap(enable);
        }
    }
    [[nodiscard]] bool fswap_enabled() const noexcept { return m_fswap; }

    void seek(std::uint64_t sample);
    void rewind() { seek(0); }
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

private:
    struct File {
        SearchFile reader;
        std::uint64_t first{0};
        std::uint64_t count{0};
    };

    [[nodiscard]] std::size_t index_at(std::uint64_t sample) const;

    template <typename T, typename Read>
    std::uint64_t read_across(std::uint64_t count,
                              std::span<T> dest,
                              std::uint64_t per_sample,
                              Read&& read);

    std::vector<File> m_files;
    sigproc::FilterbankHeader m_header{};
    std::string m_source;
    std::string m_telescope;
    std::string m_pol_type;
    std::uint64_t m_nchan{0};
    std::uint64_t m_npol{0};
    std::uint64_t m_nsblk{0};
    std::uint64_t m_nsamples{0};
    std::uint64_t m_sample{0};
    std::uint64_t m_bytes_per_sample{0};
    int m_input_nbits{0};
    int m_beam{0};
    bool m_sum{false};
    bool m_fswap{false};
    double m_tsamp{0.0};
    double m_tstart{0.0};
    double m_fch1{0.0};
    double m_foff{0.0};
    double m_bandwidth{0.0};
    double m_centre{0.0};
    double m_raj{0.0};
    double m_dej{0.0};
};

inline SearchStream::SearchStream(
    std::span<const std::filesystem::path> paths) {
    if (paths.empty()) {
        throw ValidationError("psrio: a PSRFITS stream needs a file");
    }
    m_files.reserve(paths.size());
    for (const auto& path : paths) {
        m_files.push_back(File{SearchFile(path), 0, 0});
    }
    const auto& first  = m_files.front().reader;
    m_header           = first.sigproc_header();
    m_source           = std::string(first.source_name());
    m_telescope        = std::string(first.telescope());
    m_pol_type         = std::string(first.pol_type());
    m_nchan            = first.nchans();
    m_npol             = first.input_npol();
    m_nsblk            = first.subint_samples();
    m_bytes_per_sample = first.bytes_per_sample();
    m_input_nbits      = first.input_nbits();
    m_beam             = first.beam();
    m_sum              = first.sums_polarisations();
    m_tsamp            = first.tsamp();
    m_tstart           = first.tstart();
    m_fch1             = first.fch1();
    m_foff             = first.foff();
    m_bandwidth        = first.bandwidth();
    m_centre           = first.center_frequency();
    m_raj              = first.raj();
    m_dej              = first.dej();

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
            stream_detail::match_search(first, file.reader,
                                        paths[index].string());
            stream_detail::require_nsuboffs(m_files[index - 1U].reader,
                                            file.reader, paths[index].string());
            ::psrio::detail::require_sample_gap(
                m_files[index - 1U].reader.tstart(), m_files[index - 1U].count,
                m_files[index - 1U].reader.tsamp(), file.reader.tstart(),
                paths[index].string());
        }
        cursor += file.count;
    }
    m_nsamples = cursor;
    m_header.declared_nsamples.reset();
    m_header.data_bytes = m_nsamples * m_bytes_per_sample;
    m_header.file_bytes = m_header.header_bytes + m_header.data_bytes;
}

inline void SearchStream::seek(std::uint64_t sample) {
    if (sample > m_nsamples) {
        throw ValidationError("psrio: seek is past the readable samples");
    }
    m_sample = sample;
}

inline std::size_t SearchStream::index_at(std::uint64_t sample) const {
    for (std::size_t index = 0; index < m_files.size(); ++index) {
        const auto& file = m_files[index];
        if (sample < file.first + file.count) {
            return index;
        }
    }
    return m_files.size() - 1U;
}

template <typename T, typename Read>
inline std::uint64_t SearchStream::read_across(std::uint64_t count,
                                               std::span<T> dest,
                                               std::uint64_t per_sample,
                                               Read&& read) {
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

inline std::uint64_t SearchStream::read_block(std::uint64_t count,
                                              std::span<std::byte> dest) {
    return read_across(
        count, dest, bytes_per_sample(),
        [](SearchFile& reader, std::uint64_t take, std::span<std::byte> piece) {
            return reader.read_block(take, piece);
        });
}

inline std::uint64_t SearchStream::read_samples(std::uint64_t count,
                                                std::span<float> dest) {
    return read_across(
        count, dest, nchans(),
        [](SearchFile& reader, std::uint64_t take, std::span<float> piece) {
            return reader.read_samples(take, piece);
        });
}

inline std::uint64_t
SearchStream::read_samples(std::uint64_t /*count*/,
                           std::span<std::uint8_t> /*dest*/) {
    throw ValidationError("psrio: calibrated PSRFITS samples unpack to float");
}

inline std::uint64_t
SearchStream::read_samples(std::uint64_t /*count*/,
                           std::span<std::uint16_t> /*dest*/) {
    throw ValidationError("psrio: calibrated PSRFITS samples unpack to float");
}

inline std::uint64_t SearchStream::read_bytes(std::uint64_t nbytes,
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
    return read_block(samples, dest) * stride;
}

} // namespace psrio::formats::psrfits
