#include "common.hpp"

#include <CLI/CLI.hpp>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fmt/format.h>
#include <fstream>
#include <memory>
#include <optional>
#include <psrio/detail/exceptions.hpp>
#include <psrio/stats.hpp>
#include <string>
#include <vector>

namespace psrio::cli {

namespace {

struct StatsOptions {
    std::vector<std::filesystem::path> inputs;
    std::filesystem::path output;
    std::uint64_t start{0};
    std::optional<std::uint64_t> nsamps;
    std::uint64_t gulp{16384};
    std::string method{"mean"};
};

stats::ChannelStats
read_stats(Input& input, const StatsOptions& opts, stats::Moments order) {
    if (input.is_series()) {
        throw ValidationError(
            "psrio: channel statistics need a filterbank, not a time series");
    }
    return stats::compute(input.block(), stats::ComputeOptions{
                                             .start  = opts.start,
                                             .nsamps = opts.nsamps,
                                             .gulp   = opts.gulp,
                                             .order  = order,
                                         });
}

/// Frequency of column @p col (IF-major, channels in file order).
double column_freq(const std::vector<double>& freqs, std::size_t col) {
    return freqs[col % freqs.size()];
}

std::ofstream open_csv(const std::filesystem::path& path) {
    std::ofstream file(path);
    if (!file) {
        throw IoError(fmt::format("psrio: cannot write {}", path.string()));
    }
    return file;
}

void run_stats(const StatsOptions& opts, const Context& ctx) {
    Input input(opts.inputs, {}, {}, 0, ctx.drop_cache, !ctx.force_contiguity);
    const auto result = read_stats(input, opts, stats::Moments::kFull);
    const auto freqs  = channel_frequencies(input.header());
    const auto mean   = result.mean();
    const auto stddev = result.stddev();
    const auto skew   = result.skew();
    const auto kurt   = result.kurtosis();
    const auto lo     = result.minima();
    const auto hi     = result.maxima();
    if (!opts.output.empty()) {
        auto file = open_csv(opts.output);
        file << "chan,freq_mhz,mean,std,skew,kurtosis,min,max\n";
        for (std::size_t col = 0; col < mean.size(); ++col) {
            file << fmt::format("{},{},{},{},{},{},{},{}\n", col,
                                column_freq(freqs, col), mean[col], stddev[col],
                                skew[col], kurt[col], lo[col], hi[col]);
        }
        ctx.info("psrio stats: {} channels over {} samples to {}", mean.size(),
                 result.count(), opts.output.string());
        return;
    }
    ctx.out << "chan\tfreq_mhz\tmean\tstd\tskew\tkurtosis\tmin\tmax\n";
    for (std::size_t col = 0; col < mean.size(); ++col) {
        ctx.out << fmt::format(
            "{}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\n", col,
            column_freq(freqs, col), mean[col], stddev[col], skew[col],
            kurt[col], lo[col], hi[col]);
    }
}

void run_bandpass(const StatsOptions& opts, const Context& ctx) {
    if (opts.method != "mean" && opts.method != "sum") {
        throw ValidationError(fmt::format(
            "psrio: bandpass method must be mean or sum, got {}", opts.method));
    }
    Input input(opts.inputs, {}, {}, 0, ctx.drop_cache, !ctx.force_contiguity);
    const auto result = read_stats(input, opts, stats::Moments::kMean);
    const auto freqs  = channel_frequencies(input.header());
    const auto values =
        opts.method == "sum"
            ? result.sum()
            : std::vector<double>(result.mean().begin(), result.mean().end());
    if (!opts.output.empty()) {
        auto file = open_csv(opts.output);
        file << "freq_mhz," << opts.method << '\n';
        for (std::size_t col = 0; col < values.size(); ++col) {
            file << fmt::format("{:.6f},{:.6f}\n", column_freq(freqs, col),
                                values[col]);
        }
        ctx.info("psrio bandpass: {} channels over {} samples to {}",
                 values.size(), result.count(), opts.output.string());
        return;
    }
    ctx.out << "freq_mhz\t" << opts.method << '\n';
    for (std::size_t col = 0; col < values.size(); ++col) {
        ctx.out << fmt::format("{:.6f}\t{:.6f}\n", column_freq(freqs, col),
                               values[col]);
    }
}

void add_common(CLI::App& sub, StatsOptions& opts, const char* output_help) {
    sub.add_option("inputs", opts.inputs,
                   "Input files in time order: .fil or .sf/.fits, or one .h5")
        ->required()
        ->expected(1, -1)
        ->check(CLI::ExistingFile);
    sub.add_option("-o,--output", opts.output, output_help);
    sub.add_option("-s,--start", opts.start, "First sample")
        ->capture_default_str();
    sub.add_option("-n,--nsamps", opts.nsamps,
                   "Samples to read (default: to the end)");
    sub.add_option("-g,--gulp", opts.gulp, "Samples per block")
        ->capture_default_str()
        ->check(CLI::PositiveNumber);
    add_threads_option(sub);
}

} // namespace

void add_stats(CLI::App& app, Context& ctx) {
    const auto opts = std::make_shared<StatsOptions>();
    auto* sub       = app.add_subcommand(
        "stats",
        "Per-channel mean, std, skewness, kurtosis, min and max (spp_stats)");
    add_common(*sub, *opts,
               "Write CSV (chan, freq_mhz, mean, std, skew, kurtosis, min, "
               "max) instead of a table on stdout");
    sub->callback([opts, &ctx] { run_stats(*opts, ctx); });

    const auto bp_opts = std::make_shared<StatsOptions>();
    auto* bandpass     = app.add_subcommand(
        "bandpass", "Mean or summed bandpass per channel (spp_bandpass)");
    add_common(*bandpass, *bp_opts,
               "Write CSV (freq_mhz, value) instead of a table on stdout");
    bandpass->add_option("-m,--method", bp_opts->method, "mean or sum")
        ->capture_default_str()
        ->check(CLI::IsMember({"mean", "sum"}));
    bandpass->callback([bp_opts, &ctx] { run_bandpass(*bp_opts, ctx); });
}

} // namespace psrio::cli
