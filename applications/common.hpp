#pragma once

#include <CLI/CLI.hpp>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fmt/format.h>
#include <optional>
#include <ostream>
#include <psrio/block_source.hpp>
#include <psrio/formats/psrfits/reader.hpp>
#include <psrio/formats/sigproc/header.hpp>
#include <psrio/formats/sigproc/writer.hpp>
#include <psrio/timeseries.hpp>
#include <span>
#include <string>
#include <string_view>
#include <variant>

namespace psrio::cli {

namespace sigproc = ::psrio::formats::sigproc;

/// Output streams and verbosity shared by every subcommand.
struct Context {
    std::ostream&
        out; // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)
    std::ostream&
        err; // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)
    int verbosity{0};
    bool quiet{false};

    /// Progress and summary messages, suppressed by --quiet.
    template <typename... Args>
    void info(fmt::format_string<Args...> format, Args&&... args) const {
        if (!quiet) {
            err << fmt::format(format, std::forward<Args>(args)...) << '\n';
        }
    }

    /// Details shown with --verbose.
    template <typename... Args>
    void debug(fmt::format_string<Args...> format, Args&&... args) const {
        if (!quiet && verbosity > 0) {
            err << fmt::format(format, std::forward<Args>(args)...) << '\n';
        }
    }

    template <typename... Args>
    void warn(fmt::format_string<Args...> format, Args&&... args) const {
        err << "psrio: warning: "
            << fmt::format(format, std::forward<Args>(args)...) << '\n';
    }
};

/// Default streaming block, in time samples.
inline constexpr std::uint64_t kDefaultGulp = 65536;

/// True for a PRESTO `.dat` or `.inf` path.
[[nodiscard]] bool is_presto_path(const std::filesystem::path& path);

/// True for a SIGPROC `.tim` path.
[[nodiscard]] bool is_tim_path(const std::filesystem::path& path);

/// True for a `.tim`, `.dat`, or `.inf` time series path.
[[nodiscard]] inline bool is_series_path(const std::filesystem::path& path) {
    return is_tim_path(path) || is_presto_path(path);
}

/// A filterbank or time series opened for streaming.
///
/// SIGPROC filterbanks and FBH5 files stream through a BlockSource; `.tim`
/// and PRESTO `.dat`/`.inf` series through a TimeSeriesReader. header()
/// describes either one in SIGPROC form.
class Input {
public:
    /// @p decode sets how search-mode PSRFITS DATA becomes samples; other
    /// formats ignore it.
    /// @throws psrio::Error if the file cannot be opened.
    explicit Input(const std::filesystem::path& path,
                   const formats::psrfits::DecodeOptions& decode = {});

    /// Open one file, or several filterbank or search-mode PSRFITS files in
    /// the order given. The files must be one contiguous observation.
    /// @throws psrio::Error if a file cannot be opened or the files do not
    ///         join.
    explicit Input(std::span<const std::filesystem::path> paths,
                   const formats::psrfits::DecodeOptions& decode = {});

    [[nodiscard]] const std::filesystem::path& path() const noexcept {
        return m_path;
    }
    [[nodiscard]] const sigproc::FilterbankHeader& header() const noexcept {
        return m_header;
    }
    [[nodiscard]] bool is_series() const noexcept {
        return m_series.has_value();
    }
    [[nodiscard]] std::uint64_t nsamples() const;
    [[nodiscard]] std::uint64_t nchans() const;
    [[nodiscard]] std::uint64_t nifs() const;
    [[nodiscard]] int nbits() const;
    [[nodiscard]] double tsamp() const;
    /// True when the input is a search-mode PSRFITS file.
    [[nodiscard]] bool is_psrfits() const noexcept { return m_psrfits; }
    /// Bits stored in a PSRFITS DATA column. Zero for other inputs.
    [[nodiscard]] int stored_nbits() const noexcept { return m_stored_nbits; }
    /// Polarisation products stored in a PSRFITS DATA column.
    [[nodiscard]] std::uint64_t stored_npol() const noexcept {
        return m_stored_npol;
    }
    /// Samples in one PSRFITS sub-integration. Zero for other inputs.
    [[nodiscard]] std::uint64_t subint_samples() const noexcept {
        return m_subint_samples;
    }
    /// True when the PSRFITS stream sums the first two polarisations.
    [[nodiscard]] bool sums_polarisations() const noexcept {
        return m_sum_feeds;
    }

    /// Search-mode PSRFITS metadata, for inspect.
    struct PsrfitsInfo {
        std::string frontend;
        std::string backend;
        std::string obs_mode;
        std::string pol_type;
        /// SUBINT rows of each file, in time order.
        std::vector<std::uint64_t> rows;
        /// NSUBOFFS of each file, when present.
        std::vector<std::optional<std::int64_t>> nsuboffs;
    };
    /// PSRFITS metadata; empty for other inputs.
    [[nodiscard]] const PsrfitsInfo& psrfits_info() const noexcept {
        return m_info;
    }

    /// Values in one time sample: nifs * nchans.
    [[nodiscard]] std::uint64_t row() const { return nifs() * nchans(); }

    /// The block reader of a filterbank input.
    /// @throws ValidationError for a time series.
    [[nodiscard]] BlockSource& block();

    void seek(std::uint64_t sample);

    /// Read up to @p count samples as float into the front of @p dest.
    /// @return Samples read.
    std::uint64_t read_floats(std::uint64_t count, std::span<float> dest);

private:
    std::filesystem::path m_path;
    sigproc::FilterbankHeader m_header;
    std::optional<BlockSource> m_block;
    std::optional<TimeSeriesReader> m_series;
    bool m_psrfits{false};
    bool m_sum_feeds{false};
    int m_stored_nbits{0};
    std::uint64_t m_stored_npol{0};
    std::uint64_t m_subint_samples{0};
    PsrfitsInfo m_info;
};

/// A SIGPROC filterbank or a `.tim`/`.dat` time series being written.
///
/// The format follows the extension. A `.tim` or `.dat` output holds one
/// channel and one IF as 32-bit floats; anything else is a filterbank with
/// the header's nbits.
class Output {
public:
    /// @throws ValidationError if @p header does not fit the format.
    Output(const std::filesystem::path& path, sigproc::FilterbankHeader header);

    [[nodiscard]] bool is_series() const noexcept {
        return std::holds_alternative<TimeSeriesWriter>(m_writer);
    }

    /// The filterbank writer.
    /// @throws ValidationError for a time series output.
    [[nodiscard]] sigproc::FilterbankWriter& sigproc_writer();

    void write_floats(std::uint64_t count, std::span<const float> data);
    [[nodiscard]] std::uint64_t nsamples_written() const;
    [[nodiscard]] std::filesystem::path path() const;
    void close();

private:
    std::variant<sigproc::FilterbankWriter, TimeSeriesWriter> m_writer;
};

/// @throws ValidationError if @p output names the same file as @p input.
void require_distinct(const std::filesystem::path& input,
                      const std::filesystem::path& output);

/// @throws ValidationError if @p output names the same file as any input.
void require_distinct(std::span<const std::filesystem::path> inputs,
                      const std::filesystem::path& output);

/// Time samples [start, start + count) selected by sample or time options.
struct SampleRange {
    std::uint64_t start{0};
    std::uint64_t count{0};
};

/// Resolve the start and length options against @p nsamples. A length past
/// the end is cut to the end with a warning.
/// @throws ValidationError if the start is outside the data.
[[nodiscard]] SampleRange resolve_range(const Context& ctx,
                                        std::uint64_t nsamples,
                                        double tsamp,
                                        std::optional<std::uint64_t> start,
                                        std::optional<double> start_time,
                                        std::optional<std::uint64_t> count,
                                        std::optional<double> duration);

/// Channel centres of @p header in file order.
[[nodiscard]] std::vector<double>
channel_frequencies(const sigproc::FilterbankHeader& header);

/// Describe one channel as `fch1`/`foff` when @p freqs is a uniform grid, or
/// as a frequency table otherwise.
void set_channel_frequencies(sigproc::FilterbankHeader& header,
                             std::span<const double> freqs);

/// Add `-j,--threads` to @p sub. The count applies to the OpenMP loops in
/// the readers and digifil; without it OpenMP picks its default
/// (`OMP_NUM_THREADS`, else every core this process may run on).
void add_threads_option(CLI::App& sub);

/// Threads a parallel step uses by default: OpenMP's count, which follows
/// `OMP_NUM_THREADS` and `OMP_THREAD_LIMIT` (set by batch systems such as
/// HTCondor), else the hardware thread count.
[[nodiscard]] int default_threads();

/// One key changed by edit_sigproc_header().
struct HeaderChange {
    std::string key;
    std::string before;
    std::string after;
};

/// Set SIGPROC header keys of @p path in place from `key=value`
/// @p assignments, as `spp_header update`: all or nothing, the header size
/// never changes (a shorter string is padded with spaces, a longer one is
/// refused), keys absent from the file cannot be added, and nbits, nchans,
/// and nifs cannot change. With @p dry_run the file is left alone.
/// @throws ValidationError if any assignment cannot be made.
std::vector<HeaderChange>
edit_sigproc_header(const std::filesystem::path& path,
                    std::span<const std::string> assignments,
                    bool dry_run);

/// Register each subcommand on @p app.
void add_inspect(CLI::App& app, Context& ctx);
void add_slice(CLI::App& app, Context& ctx);
void add_splice(CLI::App& app, Context& ctx);
void add_decimate(CLI::App& app, Context& ctx);
void add_convert(CLI::App& app, Context& ctx);
void add_digifil(CLI::App& app, Context& ctx);
void add_cutout(CLI::App& app, Context& ctx);
void add_stats(CLI::App& app, Context& ctx);
void add_clean(CLI::App& app, Context& ctx);
void add_dedisperse(CLI::App& app, Context& ctx);

} // namespace psrio::cli
