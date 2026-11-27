#pragma once

#include "psrio/astro.hpp"
#include "psrio/common/types.hpp"
#include "psrio/detail/contiguity.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/skip.hpp"
#include "psrio/formats/psrfits/reader.hpp"
#include "psrio/formats/sigproc/header.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <exception>
#include <filesystem>
#include <fmt/format.h>
#include <functional>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

namespace psrio::formats::psrfits {

namespace stream_detail {

struct FirstMeta {
    std::uint64_t nchans{0};
    int input_nbits{0};
    std::uint64_t input_npol{0};
    std::string pol_type;
    bool sums_polarisations{false};
    std::uint64_t subint_samples{0};
    double tsamp{0.0};
    double fch1{0.0};
    double foff{0.0};
    std::string source_name;
    std::string telescope;
    int beam{0};
    std::vector<double> frequency_table;
    std::optional<double> refdm;

    static FirstMeta from(const SearchFile& first,
                          const sigproc::FilterbankHeader& hdr) {
        return FirstMeta{
            .nchans             = first.nchans(),
            .input_nbits        = first.input_nbits(),
            .input_npol         = first.input_npol(),
            .pol_type           = std::string(first.pol_type()),
            .sums_polarisations = first.sums_polarisations(),
            .subint_samples     = first.subint_samples(),
            .tsamp              = first.tsamp(),
            .fch1               = first.fch1(),
            .foff               = first.foff(),
            .source_name        = std::string(first.source_name()),
            .telescope          = std::string(first.telescope()),
            .beam               = first.beam(),
            .frequency_table    = hdr.frequency_table,
            .refdm              = hdr.refdm,
        };
    }
};

inline void match_search(const FirstMeta& first,
                         const SearchFile& next,
                         std::string_view path) {
    const auto fail = [&](std::string_view detail_msg) {
        throw FormatError(fmt::format(
            "psrio: {} does not match the first file's observing setup ({})",
            path, detail_msg));
    };
    if (next.nchans() != first.nchans) {
        fail(fmt::format("nchans ({} vs {})", next.nchans(), first.nchans));
    }
    if (next.input_nbits() != first.input_nbits) {
        fail(fmt::format("NBITS ({} vs {})", next.input_nbits(),
                         first.input_nbits));
    }
    if (next.input_npol() != first.input_npol) {
        fail(fmt::format("NPOL ({} vs {})", next.input_npol(),
                         first.input_npol));
    }
    if (next.pol_type() != first.pol_type) {
        fail(fmt::format("POL_TYPE ('{}' vs '{}')", next.pol_type(),
                         first.pol_type));
    }
    if (next.sums_polarisations() != first.sums_polarisations) {
        fail("polarisation reduction");
    }
    if (next.subint_samples() != first.subint_samples) {
        fail(fmt::format("NSBLK ({} vs {})", next.subint_samples(),
                         first.subint_samples));
    }
    if (std::abs(next.tsamp() - first.tsamp) >
        (1.0e-9 * std::max(std::abs(first.tsamp), 1.0e-12))) {
        fail(fmt::format("tsamp ({} vs {})", next.tsamp(), first.tsamp));
    }
    if (std::abs(next.fch1() - first.fch1) > detail::kFreqAtolMhz) {
        fail(fmt::format("fch1 ({} vs {})", next.fch1(), first.fch1));
    }
    if (std::abs(next.foff() - first.foff) > detail::kFreqAtolMhz) {
        fail(fmt::format("foff ({} vs {})", next.foff(), first.foff));
    }
    if (next.source_name() != first.source_name) {
        fail(fmt::format("source_name ('{}' vs '{}')", next.source_name(),
                         first.source_name));
    }
    if (next.telescope() != first.telescope) {
        fail(fmt::format("telescope ('{}' vs '{}')", next.telescope(),
                         first.telescope));
    }
    if (next.beam() != first.beam) {
        fail(fmt::format("beam ({} vs {})", next.beam(), first.beam));
    }
    const auto right = next.sigproc_header();
    if (first.frequency_table.size() != right.frequency_table.size()) {
        fail("frequency table length");
    }
    for (std::size_t i = 0; i < first.frequency_table.size(); ++i) {
        if (std::abs(first.frequency_table[i] - right.frequency_table[i]) >
            detail::kFreqAtolMhz) {
            fail(fmt::format("frequency table channel {} ({} vs {})", i,
                             right.frequency_table[i],
                             first.frequency_table[i]));
        }
    }
    if (first.refdm.has_value() != right.refdm.has_value() ||
        (first.refdm.has_value() &&
         std::abs(*first.refdm - *right.refdm) > 1.0e-6)) {
        fail("CHAN_DM");
    }
}

inline void require_nsuboffs(const SearchFile& previous,
                             const SearchFile& next,
                             std::string_view path,
                             bool& all_zero_nsuboffs) {
    const auto before = previous.nsuboffs();
    const auto after  = next.nsuboffs();
    if (before.has_value() != after.has_value()) {
        throw FormatError(fmt::format(
            "psrio: {} must carry NSUBOFFS when the other files do", path));
    }
    if (!before.has_value()) {
        return;
    }
    // Some recorders leave NSUBOFFS=0 in every split file even when STT_*
    // advances contiguously. Allow that only if every file so far has 0.
    if (all_zero_nsuboffs && *before == 0 && *after == 0) {
        return;
    }
    all_zero_nsuboffs = false;
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

[[nodiscard]] inline bool earlier_start(const SearchFile& left,
                                        const SearchFile& right) noexcept {
    if (left.stt_imjd() != right.stt_imjd()) {
        return left.stt_imjd() < right.stt_imjd();
    }
    return (left.stt_smjd() + left.stt_offs()) <
           (right.stt_smjd() + right.stt_offs());
}

} // namespace stream_detail

/**
 * @brief One contiguous sequence of search-mode PSRFITS files.
 *
 * Files are automatically sorted into chronological order by their exact
 * `STT_IMJD`/`STT_SMJD`/`STT_OFFS` start times (matching dspsr's
 * `SerialFiles`), so shell globs over unpadded chunk numbers work out of the
 * box. Only headers are inspected during construction and file descriptors are
 * closed immediately afterward so streams over hundreds or thousands of files
 * keep at most two OS file descriptors open at any time.
 */
class SearchStream {
public:
    using OpenProgressCallback =
        std::function<void(std::size_t done, std::size_t total)>;

    explicit SearchStream(std::span<const std::filesystem::path> paths)
        : SearchStream(paths, DecodeOptions{}, {}, true, false, {}) {}

    /// Every file is decoded with @p options (see SearchFile). When
    /// @p check_contiguity is false (as dspsr's `-cont`), files are kept in
    /// the given order and timestamp / NSUBOFFS gap checks are skipped while
    /// observing-setup compatibility is still verified. When @p lazy is also
    /// true (with @p check_contiguity false), only the first file is opened at
    /// construction and subsequent files are opened and closed one by one as
    /// the stream advances.
    SearchStream(
        std::span<const std::filesystem::path> paths,
        DecodeOptions options,
        OpenProgressCallback on_open                              = {},
        bool check_contiguity                                     = true,
        bool lazy                                                 = false,
        std::shared_ptr<::psrio::detail::StreamProgress> progress = {});

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
    /// Files in the stream, in chronological order.
    [[nodiscard]] std::size_t nfiles() const noexcept { return m_files.size(); }
    /// File @p index of the stream, for its own metadata (rows, NSUBOFFS).
    [[nodiscard]] const SearchFile& file(std::size_t index) const {
        const auto& slot = m_files.at(index).reader;
        if (!slot.has_value()) {
            throw ValidationError("psrio: stream file is not currently open");
        }
        return *slot;
    }
    /// Paths of the stream in chronological order.
    [[nodiscard]] std::vector<std::filesystem::path> paths() const {
        std::vector<std::filesystem::path> out;
        out.reserve(m_files.size());
        for (const auto& entry : m_files) {
            out.push_back(entry.path);
        }
        return out;
    }
    /// Cumulative ending sample index of each file in chronological order.
    [[nodiscard]] std::vector<std::uint64_t> file_ends() const {
        std::vector<std::uint64_t> out;
        out.reserve(m_files.size());
        for (const auto& entry : m_files) {
            out.push_back(entry.first + entry.count);
        }
        return out;
    }
    /// Total on-disk byte size across all files in the stream.
    [[nodiscard]] std::uint64_t total_input_bytes() const noexcept {
        return m_total_bytes;
    }
    /// True if chronological sorting changed the input order.
    [[nodiscard]] bool reordered() const noexcept { return m_reordered; }
    /// Configure read-ahead depth (in ReadAhead::kBlockBytes blocks) per file.
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
    /// `CHAN_DM` of the first file, which every file must match.
    [[nodiscard]] std::optional<double> refdm() const noexcept {
        return m_header.refdm;
    }

    [[nodiscard]] std::uint64_t nchans() const noexcept { return m_nchan; }
    [[nodiscard]] static std::uint64_t nifs() noexcept { return 1; }
    [[nodiscard]] static int nbits() noexcept { return 32; }
    [[nodiscard]] static SampleType sample_type() noexcept {
        return SampleType::kFloat32;
    }
    [[nodiscard]] std::uint64_t bytes_per_sample() const noexcept {
        return m_bytes_per_sample;
    }
    [[nodiscard]] std::uint64_t nsamples() const noexcept { return m_nsamples; }
    [[nodiscard]] bool has_nsamples() const noexcept {
        return !m_lazy || m_all_opened;
    }
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
            if (file.reader.has_value()) {
                file.reader->set_fswap(enable);
            }
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
    static std::uint64_t read_samples(std::uint64_t count,
                                      std::span<std::uint8_t> dest);
    static std::uint64_t read_samples(std::uint64_t count,
                                      std::span<std::uint16_t> dest);
    std::uint64_t read_bytes(std::uint64_t nbytes, std::span<std::byte> dest);

    template <typename T>
        requires std::same_as<T, float> || std::same_as<T, std::uint8_t> ||
                 std::same_as<T, std::uint16_t>
    std::uint64_t read(std::uint64_t count, std::span<T> dest) {
        return read_samples(count, dest);
    }

    [[nodiscard]] std::vector<std::byte> read_block(std::uint64_t count) {
        std::vector<std::byte> out(
            static_cast<std::size_t>(count * m_bytes_per_sample));
        const auto actual = read_block(count, out);
        out.resize(static_cast<std::size_t>(actual * m_bytes_per_sample));
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
        std::vector<T> buffer(static_cast<std::size_t>(count * m_nchan));
        const auto actual = read_samples(count, std::span<T>(buffer));
        buffer.resize(static_cast<std::size_t>(actual * m_nchan));
        return buffer;
    }

private:
    struct File {
        std::optional<SearchFile> reader;
        std::filesystem::path path;
        std::size_t input_index{0};
        std::uint64_t first{0};
        std::uint64_t count{0};
        std::uint64_t bytes{0};
        bool opened{false};

        [[nodiscard]] SearchFile& open_reader() {
            if (!reader.has_value()) {
                throw ValidationError(
                    "psrio: stream file is not currently open");
            }
            return *reader;
        }
        [[nodiscard]] const SearchFile& open_reader() const {
            if (!reader.has_value()) {
                throw ValidationError(
                    "psrio: stream file is not currently open");
            }
            return *reader;
        }
    };

    void init_metadata(const SearchFile& first);
    void ensure_lazy_file(std::size_t idx);
    void publish_progress(std::size_t idx) const noexcept;
    [[nodiscard]] std::size_t index_at(std::uint64_t sample) const;

    template <typename T, typename Read>
    std::uint64_t read_across(std::uint64_t count,
                              std::span<T> dest,
                              std::uint64_t per_sample,
                              const Read& read);

    std::vector<File> m_files;
    sigproc::FilterbankHeader m_header{};
    stream_detail::FirstMeta m_first_meta{};
    DecodeOptions m_decode_options{};
    std::shared_ptr<::psrio::detail::StreamProgress> m_progress;
    std::string m_source;
    std::string m_telescope;
    std::string m_pol_type;
    std::uint64_t m_nchan{0};
    std::uint64_t m_npol{0};
    std::uint64_t m_nsblk{0};
    std::uint64_t m_nsamples{0};
    std::uint64_t m_total_bytes{0};
    std::uint64_t m_bytes_before_current{0};
    std::uint64_t m_sample{0};
    std::uint64_t m_bytes_per_sample{0};
    std::size_t m_last_index{0};
    std::size_t m_current_idx{0};
    std::size_t m_opened_files{0};
    std::size_t m_read_ahead_blocks{0};
    int m_input_nbits{0};
    int m_beam{0};
    bool m_sum{false};
    bool m_fswap{false};
    bool m_drop_cache{false};
    bool m_reordered{false};
    bool m_lazy{false};
    bool m_all_opened{false};
    double m_tsamp{0.0};
    double m_tstart{0.0};
    double m_fch1{0.0};
    double m_foff{0.0};
    double m_bandwidth{0.0};
    double m_centre{0.0};
    double m_raj{0.0};
    double m_dej{0.0};
};

inline void SearchStream::init_metadata(const SearchFile& first) {
    m_header           = first.sigproc_header();
    m_first_meta       = stream_detail::FirstMeta::from(first, m_header);
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
}

inline void SearchStream::publish_progress(std::size_t idx) const noexcept {
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

inline void SearchStream::ensure_lazy_file(std::size_t idx) {
    auto& entry = m_files[idx];
    if (entry.reader.has_value()) {
        return;
    }
    SearchFile reader(entry.path, m_decode_options);
    if (m_read_ahead_blocks > 0U) {
        reader.set_read_ahead_blocks(m_read_ahead_blocks);
    }
    reader.set_drop_cache(m_drop_cache);
    reader.set_fswap(m_fswap);
    if (!entry.opened) {
        if (idx > 0U) {
            stream_detail::match_search(m_first_meta, reader,
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

inline SearchStream::SearchStream(
    std::span<const std::filesystem::path> paths,
    DecodeOptions options,
    OpenProgressCallback on_open,
    bool check_contiguity,
    bool lazy,
    std::shared_ptr<::psrio::detail::StreamProgress> progress)
    : m_decode_options(options),
      m_progress(std::move(progress)) {
    if (paths.empty()) {
        throw ValidationError("psrio: a PSRFITS stream needs a file");
    }
    const auto total                  = paths.size();
    m_decode_options.check_contiguity = check_contiguity;
    options.check_contiguity          = check_contiguity;

    if (lazy && !check_contiguity && total > 1U) {
        // Lazy --cont mode: open only file 0 at startup; stream remaining
        // files one by one as the cursor advances.
        m_lazy = true;
        m_files.reserve(total);
        for (std::size_t index = 0; index < total; ++index) {
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
        init_metadata(m_files.front().open_reader());
        m_header.declared_nsamples.reset();
        publish_progress(0U);
        return;
    }

    std::vector<std::optional<File>> opened(total);

    // Open headers with bounded concurrency (at most 16 open FDs at a time)
    // and close each FD immediately after metadata parsing so 500+ files never
    // exhaust RLIMIT_NOFILE.
    const auto hw =
        std::max<std::size_t>(std::thread::hardware_concurrency(), 1U);
    const std::size_t workers =
        total <= 4U ? 1U : std::min<std::size_t>({total, hw, 16U});

    if (workers == 1U) {
        for (std::size_t index = 0; index < total; ++index) {
            SearchFile reader(paths[index], options);
            const auto fsize = reader.file_size();
            reader.close_handle();
            opened[index].emplace(File{
                .reader      = std::move(reader),
                .path        = paths[index],
                .input_index = index,
                .first       = 0,
                .count       = 0,
                .bytes       = fsize,
                .opened      = true,
            });
            if (on_open) {
                on_open(index + 1U, total);
            }
        }
    } else {
        std::atomic<std::size_t> next_idx{0};
        std::atomic<std::size_t> done_count{0};
        std::mutex cb_mutex;
        std::mutex err_mutex;
        std::exception_ptr first_error;
        std::atomic<bool> failed{false};

        std::vector<std::thread> threads;
        threads.reserve(workers);
        for (std::size_t w = 0; w < workers; ++w) {
            threads.emplace_back([&] {
                while (!failed.load(std::memory_order_relaxed)) {
                    const auto index =
                        next_idx.fetch_add(1U, std::memory_order_relaxed);
                    if (index >= total) {
                        return;
                    }
                    try {
                        SearchFile reader(paths[index], options);
                        const auto fsize = reader.file_size();
                        reader.close_handle();
                        opened[index].emplace(File{
                            .reader      = std::move(reader),
                            .path        = paths[index],
                            .input_index = index,
                            .first       = 0,
                            .count       = 0,
                            .bytes       = fsize,
                            .opened      = true,
                        });
                        const auto done = done_count.fetch_add(
                                              1U, std::memory_order_acq_rel) +
                                          1U;
                        if (on_open) {
                            const std::scoped_lock lock(cb_mutex);
                            on_open(done, total);
                        }
                    } catch (...) {
                        const std::scoped_lock lock(err_mutex);
                        if (!first_error) {
                            first_error = std::current_exception();
                        }
                        failed.store(true, std::memory_order_relaxed);
                        return;
                    }
                }
            });
        }
        for (auto& thread : threads) {
            thread.join();
        }
        if (first_error) {
            std::rethrow_exception(first_error);
        }
    }

    m_files.reserve(total);
    for (auto& slot : opened) {
        if (slot.has_value()) {
            m_files.push_back(std::move(*slot));
        }
    }

    // Sort files chronologically by exact split MJD (matching dspsr
    // SerialFiles) so unpadded shell globs work cleanly. When contiguity
    // checking is disabled (--cont), preserve the caller's order.
    if (check_contiguity) {
        std::ranges::stable_sort(
            m_files, [](const File& left, const File& right) {
                return stream_detail::earlier_start(left.open_reader(),
                                                    right.open_reader());
            });
        for (std::size_t index = 0; index < m_files.size(); ++index) {
            if (m_files[index].input_index != index) {
                m_reordered = true;
                break;
            }
        }
    }

    init_metadata(m_files.front().open_reader());

    std::uint64_t cursor   = 0;
    bool all_zero_nsuboffs = true;
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
            stream_detail::match_search(m_first_meta, cur, file.path.string());
            if (check_contiguity) {
                if (prev_reader.stt_imjd() == cur.stt_imjd() &&
                    prev_reader.stt_smjd() == cur.stt_smjd() &&
                    prev_reader.stt_offs() == cur.stt_offs()) {
                    throw FormatError(fmt::format(
                        "psrio: {} has the same start time as {} (duplicate "
                        "file?)",
                        file.path.string(), prev.path.string()));
                }
                stream_detail::require_nsuboffs(
                    prev_reader, cur, file.path.string(), all_zero_nsuboffs);
                ::psrio::detail::require_split_mjd_gap(
                    prev_reader.stt_imjd(), prev_reader.stt_smjd(),
                    prev_reader.stt_offs(), prev.count, prev_reader.tsamp(),
                    cur.stt_imjd(), cur.stt_smjd(), cur.stt_offs(),
                    file.path.string());
            }
        }
        cursor += file.count;
    }
    m_nsamples     = cursor;
    m_opened_files = m_files.size();
    m_all_opened   = true;
    m_header.declared_nsamples.reset();
    m_header.data_bytes = m_nsamples * m_bytes_per_sample;
    m_header.file_bytes = m_header.header_bytes + m_header.data_bytes;
}

inline void SearchStream::seek(std::uint64_t sample) {
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

inline std::size_t SearchStream::index_at(std::uint64_t sample) const {
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
inline std::uint64_t SearchStream::read_across(std::uint64_t count,
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
            const auto remaining = file.count - (local + got);
            if (remaining > 0U && remaining <= m_nsblk &&
                idx + 1U < m_files.size()) {
                ensure_lazy_file(idx + 1U);
                m_files[idx + 1U].open_reader().prefetch_start();
            }
            if (remaining == 0U) {
                // Close FD and free staging buffers as soon as a file ends.
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
        const auto remaining = file.count - (local + got);
        // Pre-warm the next file's read-ahead near the end of this file so
        // file transitions do not stall on networked/parallel filesystems.
        if (remaining > 0U && remaining <= cur.subint_samples() &&
            idx + 1U < m_files.size()) {
            m_files[idx + 1U].open_reader().prefetch_start();
        }
        if (remaining == 0U) {
            // The stream has moved past this file: close its FD and free its
            // staging so open descriptors stay bounded for 500+ file runs.
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
