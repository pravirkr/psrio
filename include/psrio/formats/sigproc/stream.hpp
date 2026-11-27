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
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <string_view>
#include <utility>
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
 * @brief One contiguous sequence of SIGPROC filterbank files.
 *
 * Files are automatically sorted into chronological order by `tstart` before
 * checking contiguity. Only headers are read at open and file descriptors are
 * closed immediately afterward so streams over hundreds of files keep O(1)
 * open descriptors. `nsamples()` is the sum. The header is the first file's
 * header with that combined length.
 */
class FilterbankStream {
public:
    using OpenProgressCallback =
        std::function<void(std::size_t done, std::size_t total)>;

    explicit FilterbankStream(
        std::span<const std::filesystem::path> paths,
        BitOrder bit_order                  = BitOrder::kLsbFirst,
        const OpenProgressCallback& on_open = {},
        bool check_contiguity               = true,
        bool lazy                           = false,
        std::shared_ptr<::psrio::detail::StreamProgress> progress = {});

    FilterbankStream(const FilterbankStream&)                = delete;
    FilterbankStream& operator=(const FilterbankStream&)     = delete;
    FilterbankStream(FilterbankStream&&) noexcept            = default;
    FilterbankStream& operator=(FilterbankStream&&) noexcept = default;
    ~FilterbankStream()                                      = default;

    [[nodiscard]] const FilterbankHeader& header() const noexcept {
        return m_header;
    }
    [[nodiscard]] std::size_t nfiles() const noexcept { return m_files.size(); }
    [[nodiscard]] std::vector<std::filesystem::path> paths() const {
        std::vector<std::filesystem::path> out;
        out.reserve(m_files.size());
        for (const auto& entry : m_files) {
            out.push_back(entry.path);
        }
        return out;
    }
    [[nodiscard]] std::vector<std::uint64_t> file_ends() const {
        std::vector<std::uint64_t> out;
        out.reserve(m_files.size());
        for (const auto& entry : m_files) {
            out.push_back(entry.first + entry.count);
        }
        return out;
    }
    [[nodiscard]] std::uint64_t total_input_bytes() const noexcept {
        return m_total_bytes;
    }
    [[nodiscard]] bool reordered() const noexcept { return m_reordered; }
    void set_read_ahead_blocks(std::size_t blocks) noexcept {
        m_read_ahead_blocks = blocks;
        for (auto& entry : m_files) {
            if (entry.reader.has_value()) {
                entry.reader->set_read_ahead_blocks(blocks);
            }
        }
    }
    /// Ask the OS to drop consumed pages from the page cache (Linux).
    void set_drop_cache(bool enable) noexcept {
        m_drop_cache = enable;
        for (auto& entry : m_files) {
            if (entry.reader.has_value()) {
                entry.reader->set_drop_cache(enable);
            }
        }
    }
    [[nodiscard]] BitOrder bit_order() const noexcept { return m_bit_order; }
    void set_bit_order(BitOrder order) noexcept {
        m_bit_order = order;
        for (auto& file : m_files) {
            if (file.reader.has_value()) {
                file.reader->set_bit_order(order);
            }
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
    [[nodiscard]] bool has_nsamples() const noexcept {
        return !m_lazy || m_all_opened;
    }
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
            if (file.reader.has_value()) {
                file.reader->set_fswap(enable);
            }
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
        std::optional<FilterbankReader> reader;
        std::filesystem::path path;
        std::size_t input_index{0};
        std::uint64_t first{0};
        std::uint64_t count{0};
        std::uint64_t bytes{0};
        bool opened{false};

        [[nodiscard]] FilterbankReader& open_reader() {
            if (!reader.has_value()) {
                throw ValidationError(
                    "psrio: stream file is not currently open");
            }
            return *reader;
        }
        [[nodiscard]] const FilterbankReader& open_reader() const {
            if (!reader.has_value()) {
                throw ValidationError(
                    "psrio: stream file is not currently open");
            }
            return *reader;
        }
    };

    void ensure_lazy_file(std::size_t idx);
    void publish_progress(std::size_t idx) const noexcept;
    [[nodiscard]] std::size_t index_at(std::uint64_t sample) const;

    template <typename T, typename Read>
    std::uint64_t read_across(std::uint64_t count,
                              std::span<T> dest,
                              std::uint64_t per_sample,
                              const Read& read);

    std::vector<File> m_files;
    FilterbankHeader m_header{};
    std::shared_ptr<::psrio::detail::StreamProgress> m_progress;
    std::uint64_t m_nsamples{0};
    std::uint64_t m_total_bytes{0};
    std::uint64_t m_bytes_before_current{0};
    std::uint64_t m_sample{0};
    std::size_t m_last_index{0};
    std::size_t m_current_idx{0};
    std::size_t m_opened_files{0};
    std::size_t m_read_ahead_blocks{0};
    BitOrder m_bit_order{BitOrder::kLsbFirst};
    bool m_fswap{false};
    bool m_drop_cache{false};
    bool m_reordered{false};
    bool m_lazy{false};
    bool m_all_opened{false};
};

inline void FilterbankStream::publish_progress(std::size_t idx) const noexcept {
    if (!m_progress || idx >= m_files.size()) {
        return;
    }
    const auto& file              = m_files[idx];
    m_progress->file_index        = idx;
    m_progress->file_first_sample = file.first;
    m_progress->file_samples      = file.count;
    m_progress->bytes_before_file = m_bytes_before_current;
    m_progress->file_bytes        = file.bytes;
}

inline void FilterbankStream::ensure_lazy_file(std::size_t idx) {
    auto& entry = m_files[idx];
    if (entry.reader.has_value()) {
        return;
    }
    FilterbankReader reader(entry.path, m_bit_order);
    if (m_read_ahead_blocks > 0U) {
        reader.set_read_ahead_blocks(m_read_ahead_blocks);
    }
    reader.set_drop_cache(m_drop_cache);
    reader.set_fswap(m_fswap);
    if (!entry.opened) {
        if (idx > 0U) {
            stream_detail::match_header(m_header, reader.header(),
                                        entry.path.string());
        }
        const auto count = reader.nsamples();
        if (m_files.size() > 1U && count == 0U) {
            throw FormatError(
                fmt::format("psrio: {} has no samples to continue a stream",
                            entry.path.string()));
        }
        entry.first  = m_nsamples;
        entry.count  = count;
        entry.bytes  = reader.file_size();
        entry.opened = true;
        m_nsamples += count;
        m_total_bytes += entry.bytes;
        ++m_opened_files;
        if (m_opened_files == m_files.size()) {
            m_all_opened = true;
        }
    }
    entry.reader.emplace(std::move(reader));
}

inline FilterbankStream::FilterbankStream(
    std::span<const std::filesystem::path> paths,
    BitOrder bit_order,
    const OpenProgressCallback& on_open,
    bool check_contiguity,
    bool lazy,
    std::shared_ptr<::psrio::detail::StreamProgress> progress)
    : m_progress(std::move(progress)),
      m_bit_order(bit_order) {
    if (paths.empty()) {
        throw ValidationError("psrio: a filterbank stream needs a file");
    }
    if (lazy && !check_contiguity && paths.size() > 1U) {
        m_lazy = true;
        m_files.reserve(paths.size());
        for (std::size_t index = 0; index < paths.size(); ++index) {
            m_files.push_back(File{
                .reader      = std::nullopt,
                .path        = paths[index],
                .input_index = index,
                .first       = 0,
                .count       = 0,
                .bytes       = 0,
                .opened      = false,
            });
        }
        ensure_lazy_file(0U);
        m_header = m_files.front().open_reader().header();
        m_header.declared_nsamples.reset();
        publish_progress(0U);
        return;
    }
    m_files.reserve(paths.size());
    for (std::size_t index = 0; index < paths.size(); ++index) {
        FilterbankReader reader(paths[index], bit_order);
        if (check_contiguity && paths.size() > 1U &&
            (!reader.header().tstart.has_value() ||
             !reader.header().tsamp.has_value())) {
            throw FormatError(fmt::format(
                "psrio: {} needs tstart and tsamp to continue a stream",
                paths[index].string()));
        }
        const auto fsize = reader.file_size();
        reader.close_handle();
        m_files.push_back(File{
            .reader      = std::move(reader),
            .path        = paths[index],
            .input_index = index,
            .first       = 0,
            .count       = 0,
            .bytes       = fsize,
            .opened      = true,
        });
        if (on_open) {
            on_open(index + 1U, paths.size());
        }
    }
    if (check_contiguity && m_files.size() > 1U) {
        std::ranges::stable_sort(m_files, [](const File& left,
                                             const File& right) {
            return left.open_reader().tstart() < right.open_reader().tstart();
        });
        for (std::size_t index = 0; index < m_files.size(); ++index) {
            if (m_files[index].input_index != index) {
                m_reordered = true;
                break;
            }
        }
    }
    m_header             = m_files.front().open_reader().header();
    std::uint64_t cursor = 0;
    for (std::size_t index = 0; index < m_files.size(); ++index) {
        auto& file      = m_files[index];
        const auto& cur = file.open_reader();
        file.first      = cursor;
        file.count      = cur.nsamples();
        m_total_bytes += file.bytes;
        if (m_files.size() > 1U && file.count == 0U) {
            throw FormatError(
                fmt::format("psrio: {} has no samples to continue a stream",
                            file.path.string()));
        }
        if (index > 0U) {
            const auto& prev        = m_files[index - 1U];
            const auto& prev_reader = prev.open_reader();
            const auto& header      = cur.header();
            stream_detail::match_header(m_header, header, file.path.string());
            if (check_contiguity) {
                if (prev_reader.tstart() == cur.tstart()) {
                    throw FormatError(fmt::format(
                        "psrio: {} has the same start time as {} (duplicate "
                        "file?)",
                        file.path.string(), prev.path.string()));
                }
                ::psrio::detail::require_sample_gap(
                    prev_reader.tstart(), prev.count, prev_reader.tsamp(),
                    cur.tstart(), file.path.string());
            }
        }
        cursor += file.count;
    }
    m_nsamples     = cursor;
    m_opened_files = m_files.size();
    m_all_opened   = true;
    m_header.declared_nsamples.reset();
    const auto stride   = m_header.bytes_per_sample();
    m_header.data_bytes = m_nsamples * stride;
    m_header.file_bytes = m_header.header_bytes + m_header.data_bytes;
}

inline void FilterbankStream::seek(std::uint64_t sample) {
    if (m_lazy) {
        if (sample == 0U) {
            m_sample               = 0U;
            m_current_idx          = 0U;
            m_bytes_before_current = 0U;
            publish_progress(0U);
            return;
        }
        while (!m_all_opened && sample > m_nsamples) {
            const auto next = m_opened_files;
            ensure_lazy_file(next);
            m_files[next].reader.reset();
        }
        if (sample > m_nsamples) {
            throw ValidationError("psrio: seek is past the readable samples");
        }
        m_sample               = sample;
        m_current_idx          = index_at(sample);
        m_bytes_before_current = 0U;
        for (std::size_t i = 0; i < m_current_idx; ++i) {
            m_bytes_before_current += m_files[i].bytes;
        }
        return;
    }
    if (sample > m_nsamples) {
        throw ValidationError("psrio: seek is past the readable samples");
    }
    m_sample = sample;
}

inline std::size_t FilterbankStream::index_at(std::uint64_t sample) const {
    if (m_last_index < m_files.size()) {
        const auto& cached = m_files[m_last_index];
        if (sample >= cached.first && sample < cached.first + cached.count) {
            return m_last_index;
        }
        if (m_last_index + 1U < m_files.size()) {
            const auto& next = m_files[m_last_index + 1U];
            if (sample >= next.first && sample < next.first + next.count) {
                return m_last_index + 1U;
            }
        }
    }
    const auto it = std::ranges::upper_bound(
        m_files, sample, std::less<>{},
        [](const File& entry) { return entry.first + entry.count; });
    if (it == m_files.end()) {
        return m_files.size() - 1U;
    }
    return static_cast<std::size_t>(it - m_files.begin());
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
    if (m_lazy) {
        std::uint64_t done = 0;
        while (done < count) {
            while (m_current_idx < m_files.size()) {
                ensure_lazy_file(m_current_idx);
                const auto& cur = m_files[m_current_idx];
                if (m_sample + done < cur.first + cur.count) {
                    break;
                }
                m_files[m_current_idx].reader.reset();
                m_bytes_before_current += cur.bytes;
                ++m_current_idx;
            }
            if (m_current_idx >= m_files.size()) {
                break;
            }
            const auto idx = m_current_idx;
            auto& file     = m_files[idx];
            auto& cur      = file.open_reader();
            publish_progress(idx);
            const auto local = (m_sample + done) - file.first;
            const auto room  = file.count - local;
            const auto take  = std::min(room, count - done);
            cur.seek(local);
            const auto piece =
                dest.subspan(static_cast<std::size_t>(done * per_sample),
                             static_cast<std::size_t>(take * per_sample));
            const auto got = read(cur, take, piece);
            if (got == 0U) {
                break;
            }
            if (local + got == file.count) {
                file.reader.reset();
                m_bytes_before_current += file.bytes;
                m_current_idx = idx + 1U;
                if (m_current_idx < m_files.size()) {
                    ensure_lazy_file(m_current_idx);
                    m_files[m_current_idx].open_reader().prefetch_start();
                    publish_progress(m_current_idx);
                } else if (m_progress) {
                    m_progress->file_index        = idx;
                    m_progress->file_first_sample = file.first;
                    m_progress->file_samples      = file.count;
                    m_progress->bytes_before_file = m_bytes_before_current;
                    m_progress->file_bytes        = 0U;
                }
            }
            done += got;
        }
        m_sample += done;
        return done;
    }
    if (m_sample > m_nsamples) {
        throw ValidationError(
            "psrio: reader cursor is past the readable samples");
    }
    const auto to_read = std::min(count, m_nsamples - m_sample);
    std::uint64_t done = 0;
    while (done < to_read) {
        const auto idx   = index_at(m_sample + done);
        m_last_index     = idx;
        auto& file       = m_files[idx];
        auto& cur        = file.open_reader();
        const auto local = (m_sample + done) - file.first;
        const auto room  = file.count - local;
        const auto take  = std::min(room, to_read - done);
        cur.seek(local);
        const auto piece =
            dest.subspan(static_cast<std::size_t>(done * per_sample),
                         static_cast<std::size_t>(take * per_sample));
        const auto got = read(cur, take, piece);
        if (got == 0U) {
            break;
        }
        if (local + got == file.count) {
            // The stream has moved past this file: close its FD and stop
            // read-ahead.
            cur.close_handle();
            if (idx + 1U < m_files.size()) {
                m_files[idx + 1U].open_reader().prefetch_start();
            }
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
