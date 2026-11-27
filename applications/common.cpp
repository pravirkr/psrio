#include "common.hpp"

#include <CLI/CLI.hpp>
#include <psrio/common/types.hpp>
#include <psrio/detail/exceptions.hpp>
#include <psrio/formats/psrfits.hpp>
#include <psrio/formats/sigproc/header.hpp>
#include <psrio/formats/sigproc/reader.hpp>
#include <psrio/formats/sigproc/stream.hpp>
#include <psrio/formats/sigproc/writer.hpp>
#include <psrio/header.hpp>
#include <psrio/open.hpp> // IWYU pragma: keep
#include <psrio/timeseries.hpp>

#ifdef PSRIO_CLI_HAS_HDF5
#include <psrio/formats/fbh5.hpp>
#endif

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fmt/format.h>
#include <fstream>
#include <functional>
#include <glob.h>
#include <optional>
#include <ranges>
#include <span>
#include <string>
#include <string_view>
#include <system_error>
#ifndef _OPENMP
#include <thread>
#endif
#include <utility>
#include <variant>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psrio::cli {

namespace {

std::string lower_extension(const std::filesystem::path& path) {
    auto ext = path.extension().string();
    std::ranges::transform(ext, ext.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return ext;
}

#ifdef PSRIO_CLI_HAS_HDF5
sigproc::FilterbankHeader from_fbh5(const formats::fbh5::FilterbankHeader& in) {
    sigproc::FilterbankHeader out;
    out.nchans        = static_cast<std::int32_t>(in.nchans);
    out.nifs          = 1;
    out.nbits         = in.nbits;
    out.tsamp         = in.tsamp;
    out.tstart        = in.tstart;
    out.fch1          = in.fch1;
    out.foff          = in.foff;
    out.data_type     = in.data_type != 0 ? in.data_type : 1;
    out.machine_id    = in.machine_id;
    out.ibeam         = in.ibeam;
    out.barycentric   = in.barycentric;
    out.pulsarcentric = in.pulsarcentric;
    out.az_start      = in.az_start;
    out.za_start      = in.za_start;
    if (in.has_telescope_id) {
        out.telescope_id = in.telescope_id;
    }
    if (!in.source_name.empty()) {
        out.source_name = in.source_name;
    }
    if (in.has_src_raj) {
        out.src_raj = in.src_raj;
    }
    if (in.has_src_dej) {
        out.src_dej = in.src_dej;
    }
    return out;
}
#endif

std::variant<sigproc::FilterbankWriter, TimeSeriesWriter>
make_writer(const std::filesystem::path& path,
            sigproc::FilterbankHeader header,
            bool drop_cache) {
    if (!is_series_path(path)) {
        if (header.data_type == 2) {
            header.data_type = 1;
        }
        sigproc::FilterbankWriter writer(path, std::move(header));
        writer.set_drop_cache(drop_cache);
        return writer;
    }
    if (header.nchans != 1 || header.nifs != 1) {
        throw ValidationError(fmt::format(
            "psrio: a {} time series holds one channel and one IF, not {} "
            "channels and {} IFs",
            path.extension().string(), header.nchans, header.nifs));
    }
    if (header.nbits != 32) {
        throw ValidationError(fmt::format(
            "psrio: a {} time series holds 32-bit floats, not {}-bit samples",
            path.extension().string(), header.nbits));
    }
    auto series      = Header::from_sigproc(header, path.string());
    series.data_type = "time series";
    return TimeSeriesWriter(path, std::move(series));
}

} // namespace

bool is_presto_path(const std::filesystem::path& path) {
    const auto ext = lower_extension(path);
    return ext == ".dat" || ext == ".inf";
}

bool is_tim_path(const std::filesystem::path& path) {
    return lower_extension(path) == ".tim";
}

Input::Input(const std::filesystem::path& path,
             const formats::psrfits::DecodeOptions& decode,
             std::size_t read_ahead_blocks,
             bool drop_cache)
    : m_path(path),
      m_paths{path} {
    std::error_code size_ec;
    m_total_input_bytes = std::filesystem::file_size(path, size_ec);
    if (size_ec) {
        m_total_input_bytes = 0U;
    }
    if (is_tim_path(path)) {
        m_series.emplace(path);
        m_header    = *m_series->sigproc_header();
        m_file_ends = {m_series->nsamples()};
        return;
    }
    if (is_presto_path(path)) {
        m_series.emplace(path);
        m_header           = m_series->header().to_sigproc();
        m_header.data_type = 2;
        m_header.nchans    = 1;
        m_header.nifs      = 1;
        m_header.nbits     = 32;
        m_file_ends        = {m_series->nsamples()};
        return;
    }
    if (formats::psrfits::is_search_path(path)) {
        formats::psrfits::SearchFile file(path, decode);
        if (read_ahead_blocks > 0U) {
            file.set_read_ahead_blocks(read_ahead_blocks);
        }
        file.set_drop_cache(drop_cache);
        m_header            = file.sigproc_header();
        m_psrfits           = true;
        m_stored_nbits      = file.input_nbits();
        m_stored_npol       = file.input_npol();
        m_subint_samples    = file.subint_samples();
        m_sum_feeds         = file.sums_polarisations();
        m_total_input_bytes = file.file_size();
        m_file_ends         = {file.nsamples()};
        m_info              = PsrfitsInfo{
            .frontend = std::string(file.frontend()),
            .backend  = std::string(file.backend()),
            .obs_mode = std::string(file.obs_mode()),
            .pol_type = std::string(file.pol_type()),
            .rows     = {file.nrows()},
            .nsuboffs = {file.nsuboffs()},
        };
        m_block.emplace(std::move(file));
        return;
    }
#ifdef PSRIO_CLI_HAS_HDF5
    if (formats::fbh5::is_hdf5_filterbank(path)) {
        formats::fbh5::FilterbankFile file(path);
        m_header    = from_fbh5(file.header());
        m_file_ends = {file.nsamples()};
        m_block.emplace(std::move(file));
        return;
    }
#else
    if (open_detail::looks_like_hdf5(path)) {
        throw FormatError(fmt::format(
            "psrio: {} is an HDF5 file; this psrio was built without "
            "PSRIO_WITH_HDF5",
            path.string()));
    }
#endif
    sigproc::FilterbankReader reader(path);
    if (read_ahead_blocks > 0U) {
        reader.set_read_ahead_blocks(read_ahead_blocks);
    }
    reader.set_drop_cache(drop_cache);
    m_header            = reader.header();
    m_total_input_bytes = reader.file_size();
    m_file_ends         = {reader.nsamples()};
    m_block.emplace(std::move(reader));
}

Input::Input(std::span<const std::filesystem::path> paths,
             const formats::psrfits::DecodeOptions& decode,
             OpenProgressCallback on_open,
             std::size_t read_ahead_blocks,
             bool drop_cache,
             bool check_contiguity) {
    if (paths.empty()) {
        throw ValidationError("psrio: an input file is required");
    }
    if (paths.size() == 1U) {
        *this = Input(paths.front(), decode, read_ahead_blocks, drop_cache);
        return;
    }
    const bool search = formats::psrfits::is_search_path(paths.front());
    for (const auto& path : paths) {
        if (is_tim_path(path) || is_presto_path(path)) {
            throw ValidationError(
                "psrio: a time series does not join with other files");
        }
        if (formats::psrfits::is_search_path(path) != search) {
            throw ValidationError(
                "psrio: every file in a stream must be the same format");
        }
#ifdef PSRIO_CLI_HAS_HDF5
        if (formats::fbh5::is_hdf5_filterbank(path)) {
            throw ValidationError(
                "psrio: an HDF5 filterbank does not join with other files");
        }
#else
        if (open_detail::looks_like_hdf5(path)) {
            throw FormatError(fmt::format(
                "psrio: {} is an HDF5 file; this psrio was built without "
                "PSRIO_WITH_HDF5",
                path.string()));
        }
#endif
    }
    if (search) {
        formats::psrfits::SearchStream stream(paths, decode, std::move(on_open),
                                              check_contiguity);
        if (read_ahead_blocks > 0U) {
            stream.set_read_ahead_blocks(read_ahead_blocks);
        }
        stream.set_drop_cache(drop_cache);
        m_paths             = stream.paths();
        m_path              = m_paths.front();
        m_file_ends         = stream.file_ends();
        m_total_input_bytes = stream.total_input_bytes();
        m_reordered         = stream.reordered();
        m_header            = stream.sigproc_header();
        m_psrfits           = true;
        m_stored_nbits      = stream.input_nbits();
        m_stored_npol       = stream.input_npol();
        m_subint_samples    = stream.subint_samples();
        m_sum_feeds         = stream.sums_polarisations();
        const auto& first   = stream.file(0);
        m_info              = PsrfitsInfo{
            .frontend = std::string(first.frontend()),
            .backend  = std::string(first.backend()),
            .obs_mode = std::string(first.obs_mode()),
            .pol_type = std::string(stream.pol_type()),
            .rows     = {},
            .nsuboffs = {},
        };
        for (std::size_t index = 0; index < stream.nfiles(); ++index) {
            m_info.rows.push_back(stream.file(index).nrows());
            m_info.nsuboffs.push_back(stream.file(index).nsuboffs());
        }
        m_block.emplace(std::move(stream));
        return;
    }
    sigproc::FilterbankStream stream(paths, BitOrder::kLsbFirst, on_open,
                                     check_contiguity);
    if (read_ahead_blocks > 0U) {
        stream.set_read_ahead_blocks(read_ahead_blocks);
    }
    stream.set_drop_cache(drop_cache);
    m_paths             = stream.paths();
    m_path              = m_paths.front();
    m_file_ends         = stream.file_ends();
    m_total_input_bytes = stream.total_input_bytes();
    m_reordered         = stream.reordered();
    m_header            = stream.header();
    m_block.emplace(std::move(stream));
}

std::size_t Input::file_index_at(std::uint64_t sample) const {
    if (m_file_ends.empty()) {
        return 0U;
    }
    const auto it = std::ranges::upper_bound(m_file_ends, sample);
    if (it == m_file_ends.end()) {
        return m_file_ends.size() - 1U;
    }
    return static_cast<std::size_t>(it - m_file_ends.begin());
}

std::uint64_t Input::nsamples() const {
    if (m_series) {
        return m_series->nsamples();
    }
    return m_block ? m_block->nsamples() : 0U;
}

std::uint64_t Input::nchans() const {
    if (m_series) {
        return 1U;
    }
    return m_block ? m_block->nchans() : 0U;
}

std::uint64_t Input::nifs() const {
    if (m_series) {
        return 1U;
    }
    return m_block ? m_block->nifs() : 0U;
}

int Input::nbits() const {
    if (m_series) {
        return m_series->nbits();
    }
    return m_block ? m_block->nbits() : 0;
}

double Input::tsamp() const {
    if (m_series) {
        return m_series->tsamp();
    }
    return m_block ? m_block->tsamp() : 0.0;
}

BlockSource& Input::block() {
    if (!m_block) {
        throw ValidationError(fmt::format(
            "psrio: {} is a time series; this command needs a filterbank",
            m_path.string()));
    }
    return *m_block;
}

void Input::seek(std::uint64_t sample) {
    if (m_series) {
        m_series->seek(sample);
    } else if (m_block) {
        m_block->seek(sample);
    }
}

std::uint64_t Input::read_floats(std::uint64_t count, std::span<float> dest) {
    const auto values = dest.first(static_cast<std::size_t>(count * row()));
    if (m_series) {
        return m_series->read_samples(count, values);
    }
    return block().read_samples(count, values);
}

Output::Output(const std::filesystem::path& path,
               sigproc::FilterbankHeader header,
               bool drop_cache)
    : m_writer(make_writer(path, std::move(header), drop_cache)) {}

sigproc::FilterbankWriter& Output::sigproc_writer() {
    auto* writer = std::get_if<sigproc::FilterbankWriter>(&m_writer);
    if (writer == nullptr) {
        throw ValidationError("psrio: a time series output takes floats only");
    }
    return *writer;
}

void Output::write_floats(std::uint64_t count, std::span<const float> data) {
    if (auto* writer = std::get_if<sigproc::FilterbankWriter>(&m_writer)) {
        writer->write_samples(count, data);
        return;
    }
    std::get<TimeSeriesWriter>(m_writer).write_samples(
        data.first(static_cast<std::size_t>(count)));
}

std::uint64_t Output::nsamples_written() const {
    return std::visit(
        [](const auto& writer) { return writer.nsamples_written(); }, m_writer);
}

std::filesystem::path Output::path() const {
    return std::visit(
        [](const auto& writer) { return std::filesystem::path(writer.path()); },
        m_writer);
}

void Output::close() {
    std::visit([](auto& writer) { writer.close(); }, m_writer);
}

namespace {

[[nodiscard]] bool has_glob_meta(std::string_view text) noexcept {
    return text.find_first_of("*?[") != std::string_view::npos;
}

void expand_one_spec(std::string_view spec,
                     std::vector<std::filesystem::path>& out);

void read_file_list(const std::filesystem::path& list_path,
                    std::vector<std::filesystem::path>& out) {
    if (!std::filesystem::is_regular_file(list_path)) {
        throw ValidationError(fmt::format(
            "psrio: file list '{}' does not exist", list_path.string()));
    }
    std::ifstream in(list_path);
    if (!in) {
        throw IoError(fmt::format("psrio: cannot open file list '{}'",
                                  list_path.string()));
    }
    const auto base_dir = list_path.parent_path();
    std::string line;
    while (std::getline(in, line)) {
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '#') {
            continue;
        }
        const auto last  = line.find_last_not_of(" \t\r\n");
        const auto entry = line.substr(first, last - first + 1U);
        const std::filesystem::path candidate(entry);
        if (!has_glob_meta(entry) && candidate.is_relative() &&
            !std::filesystem::exists(candidate) && !base_dir.empty() &&
            std::filesystem::exists(base_dir / candidate)) {
            out.push_back(base_dir / candidate);
        } else {
            expand_one_spec(entry, out);
        }
    }
}

void expand_one_spec(std::string_view spec,
                     std::vector<std::filesystem::path>& out) {
    if (spec.empty()) {
        return;
    }
    if (spec.front() == '@') {
        read_file_list(std::filesystem::path(spec.substr(1U)), out);
        return;
    }
    const std::string pattern(spec);
    if (has_glob_meta(spec) && !std::filesystem::exists(pattern)) {
        glob_t matched{};
        const int flags =
#ifdef GLOB_TILDE
            GLOB_TILDE;
#else
            0;
#endif
        // NOLINTNEXTLINE(concurrency-mt-unsafe): single-threaded CLI setup
        const int rc = ::glob(pattern.c_str(), flags, nullptr, &matched);
        struct GlobGuard {
            glob_t* g;
            explicit GlobGuard(glob_t* ptr) noexcept : g(ptr) {}
            GlobGuard(const GlobGuard&)            = delete;
            GlobGuard& operator=(const GlobGuard&) = delete;
            GlobGuard(GlobGuard&&)                 = delete;
            GlobGuard& operator=(GlobGuard&&)      = delete;
            ~GlobGuard() { ::globfree(g); }
        } const guard{&matched};
        if (rc == GLOB_NOMATCH || matched.gl_pathc == 0U) {
            throw ValidationError(
                fmt::format("psrio: pattern '{}' matched no files", spec));
        }
        if (rc != 0) {
            throw IoError(
                fmt::format("psrio: failed to expand pattern '{}'", spec));
        }
        for (std::size_t i = 0; i < matched.gl_pathc; ++i) {
            out.emplace_back(matched.gl_pathv[i]);
        }
        return;
    }
    out.emplace_back(pattern);
}

} // namespace

std::vector<std::filesystem::path>
expand_input_paths(std::span<const std::string> specs,
                   const std::optional<std::filesystem::path>& file_list) {
    std::vector<std::filesystem::path> paths;
    if (file_list.has_value()) {
        read_file_list(*file_list, paths);
    }
    for (const auto& spec : specs) {
        expand_one_spec(spec, paths);
    }
    if (paths.empty()) {
        throw ValidationError(
            "psrio: at least one input file (or --file-list) is required");
    }
    for (const auto& path : paths) {
        std::error_code ec;
        if (!std::filesystem::exists(path, ec)) {
            throw ValidationError(fmt::format(
                "psrio: input file '{}' does not exist", path.string()));
        }
        if (!std::filesystem::is_regular_file(path, ec)) {
            throw ValidationError(fmt::format(
                "psrio: '{}' is not a regular file", path.string()));
        }
    }
    return paths;
}

void require_distinct(const std::filesystem::path& input,
                      const std::filesystem::path& output) {
    std::error_code ec;
    if (std::filesystem::equivalent(input, output, ec)) {
        throw ValidationError(fmt::format(
            "psrio: output {} would overwrite the input", output.string()));
    }
}

int default_threads() {
#ifdef _OPENMP
    return std::max(1, omp_get_max_threads());
#else
    return static_cast<int>(std::max(1U, std::thread::hardware_concurrency()));
#endif
}

void add_threads_option(CLI::App& sub) {
    sub.add_option_function<int>(
           "-j,--threads",
           [](int threads) {
#ifdef _OPENMP
               omp_set_num_threads(threads);
#else
               (void)threads;
#endif
           },
           "Threads for decoding and rescaling (default: OMP_NUM_THREADS, "
           "else all available cores)")
        ->check(CLI::PositiveNumber);
}

void require_distinct(std::span<const std::filesystem::path> inputs,
                      const std::filesystem::path& output) {
    for (const auto& input : inputs) {
        require_distinct(input, output);
    }
}

SampleRange resolve_range(const Context& ctx,
                          std::uint64_t nsamples,
                          double tsamp,
                          std::optional<std::uint64_t> start,
                          std::optional<double> start_time,
                          std::optional<std::uint64_t> count,
                          std::optional<double> duration) {
    const auto to_samples = [tsamp](double seconds, std::string_view what) {
        if (!std::isfinite(seconds) || seconds < 0.0) {
            throw ValidationError(fmt::format(
                "psrio: {} must be a non-negative number of seconds", what));
        }
        return static_cast<std::uint64_t>(std::llround(seconds / tsamp));
    };
    SampleRange range;
    if (start_time) {
        range.start = to_samples(*start_time, "start time");
    } else if (start) {
        range.start = *start;
    }
    if (range.start >= nsamples) {
        throw ValidationError(fmt::format(
            "psrio: start sample {} is outside the input ({} samples)",
            range.start, nsamples));
    }
    const auto available = nsamples - range.start;
    range.count          = available;
    if (duration) {
        range.count = to_samples(*duration, "duration");
    } else if (count) {
        range.count = *count;
    }
    if (range.count == 0U) {
        throw ValidationError("psrio: the selection holds no samples");
    }
    if (range.count > available) {
        ctx.warn("{} samples requested from sample {}, only {} remain",
                 range.count, range.start, available);
        range.count = available;
    }
    return range;
}

std::vector<double>
channel_frequencies(const sigproc::FilterbankHeader& header) {
    std::vector<double> freqs(static_cast<std::size_t>(header.nchans));
    header.channel_frequencies(freqs);
    return freqs;
}

void set_channel_frequencies(sigproc::FilterbankHeader& header,
                             std::span<const double> freqs) {
    header.nchans = static_cast<std::int32_t>(freqs.size());
    header.frequency_table.clear();
    if (freqs.size() == 1U) {
        header.fch1 = freqs.front();
        return;
    }
    const double step = freqs[1] - freqs[0];
    const bool uniform =
        step != 0.0 &&
        std::ranges::all_of(
            std::views::iota(std::size_t{1}, freqs.size()), [&](std::size_t i) {
                const double expected =
                    freqs[0] + (static_cast<double>(i) * step);
                return std::abs(freqs[i] - expected) < 1e-4 * std::abs(step);
            });
    if (uniform) {
        header.fch1 = freqs.front();
        header.foff = step;
    } else {
        header.frequency_table.assign(freqs.begin(), freqs.end());
        header.fch1 = freqs.front();
        header.foff = step;
    }
}

} // namespace psrio::cli
