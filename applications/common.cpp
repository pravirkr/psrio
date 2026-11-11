#include "common.hpp"

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
#include <optional>
#include <ranges>
#include <span>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <variant>
#include <vector>

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
            sigproc::FilterbankHeader header) {
    if (!is_series_path(path)) {
        if (header.data_type == 2) {
            header.data_type = 1;
        }
        return sigproc::FilterbankWriter(path, std::move(header));
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

Input::Input(const std::filesystem::path& path) : m_path(path) {
    if (is_tim_path(path)) {
        m_series.emplace(path);
        m_header = *m_series->sigproc_header();
        return;
    }
    if (is_presto_path(path)) {
        m_series.emplace(path);
        m_header           = m_series->header().to_sigproc();
        m_header.data_type = 2;
        m_header.nchans    = 1;
        m_header.nifs      = 1;
        m_header.nbits     = 32;
        return;
    }
    if (formats::psrfits::is_search_path(path)) {
        formats::psrfits::SearchFile file(path);
        m_header         = file.sigproc_header();
        m_psrfits        = true;
        m_stored_nbits   = file.input_nbits();
        m_stored_npol    = file.input_npol();
        m_subint_samples = file.subint_samples();
        m_sum_feeds      = file.sums_polarisations();
        m_block.emplace(std::move(file));
        return;
    }
#ifdef PSRIO_CLI_HAS_HDF5
    if (formats::fbh5::is_hdf5_filterbank(path)) {
        formats::fbh5::FilterbankFile file(path);
        m_header = from_fbh5(file.header());
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
    m_header = reader.header();
    m_block.emplace(std::move(reader));
}

Input::Input(std::span<const std::filesystem::path> paths) {
    if (paths.empty()) {
        throw ValidationError("psrio: an input file is required");
    }
    if (paths.size() == 1U) {
        *this = Input(paths.front());
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
    m_path = paths.front();
    if (search) {
        formats::psrfits::SearchStream stream(paths);
        m_header         = stream.sigproc_header();
        m_psrfits        = true;
        m_stored_nbits   = stream.input_nbits();
        m_stored_npol    = stream.input_npol();
        m_subint_samples = stream.subint_samples();
        m_sum_feeds      = stream.sums_polarisations();
        m_block.emplace(std::move(stream));
        return;
    }
    sigproc::FilterbankStream stream(paths);
    m_header = stream.header();
    m_block.emplace(std::move(stream));
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
               sigproc::FilterbankHeader header)
    : m_writer(make_writer(path, std::move(header))) {}

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

void require_distinct(const std::filesystem::path& input,
                      const std::filesystem::path& output) {
    std::error_code ec;
    if (std::filesystem::equivalent(input, output, ec)) {
        throw ValidationError(fmt::format(
            "psrio: output {} would overwrite the input", output.string()));
    }
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
