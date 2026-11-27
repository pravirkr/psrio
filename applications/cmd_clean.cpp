#include "common.hpp"

#include <CLI/CLI.hpp>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fmt/format.h>
#include <fstream>
#include <functional>
#include <memory>
#include <optional>
#include <psrio/block_source.hpp>
#include <psrio/clean.hpp>
#include <psrio/detail/exceptions.hpp>
#include <psrio/detail/float_blocks.hpp>
#include <psrio/formats/sigproc/header.hpp>
#include <psrio/stats.hpp>
#include <random>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace psrio::cli {

namespace {

struct CleanOptions {
    std::vector<std::filesystem::path> inputs;
    std::filesystem::path output;
    std::uint64_t gulp{16384};
    // rfi
    std::string method{"mad"};
    double threshold{3.0};
    std::filesystem::path mask_out;
    std::filesystem::path stats_out;
    std::optional<std::uint64_t> stats_samples;
    // zap-channels
    std::vector<std::uint64_t> chans;
    std::filesystem::path mask_file;
    double fill{0.0};
    // zap-times
    std::uint64_t zap_start{0};
    std::uint64_t zap_end{0};
    std::string zap_mode{"zero"};
};

/// sigpyproc's default name: the first input's stem plus @p suffix, in the
/// working directory.
std::filesystem::path output_path(const CleanOptions& opts,
                                  std::string_view suffix) {
    if (!opts.output.empty()) {
        return opts.output;
    }
    return {opts.inputs.front().stem().string() + std::string(suffix) + ".fil"};
}

Input open_filterbank(const CleanOptions& opts, const Context& ctx) {
    Input input(opts.inputs, {}, {}, 0, ctx.drop_cache, !ctx.force_contiguity);
    if (input.is_series()) {
        throw ValidationError(
            "psrio: clean needs a filterbank, not a time series");
    }
    return input;
}

/// Write every sample of @p source to @p path with @p header, calling
/// @p edit on each block (first sample index, rows, values) before it is
/// written. Returns the samples written.
std::uint64_t write_stream(
    BlockSource& source,
    const sigproc::FilterbankHeader& header,
    const std::filesystem::path& path,
    std::uint64_t gulp,
    bool drop_cache,
    const std::function<void(std::uint64_t, std::uint64_t, std::span<float>)>&
        edit = {}) {
    Output output(path, header, drop_cache);
    detail::FloatBlocks blocks(source, 0, std::nullopt, gulp);
    std::vector<float> scratch;
    std::uint64_t written = 0;
    while (true) {
        const auto block = blocks.next();
        if (block.empty()) {
            break;
        }
        const auto rows = block.size() / blocks.row();
        if (edit) {
            scratch.assign(block.begin(), block.end());
            edit(written, rows, scratch);
            output.write_floats(rows, scratch);
        } else {
            output.write_floats(rows, block);
        }
        written += rows;
    }
    output.close();
    return written;
}

void write_mask(const std::filesystem::path& path,
                std::span<const std::uint8_t> mask) {
    std::ofstream file(path);
    if (!file) {
        throw IoError(fmt::format("psrio: cannot write {}", path.string()));
    }
    for (std::size_t index = 0; index < mask.size(); ++index) {
        if (mask[index] != 0U) {
            file << index << '\n';
        }
    }
}

std::vector<std::uint64_t> read_mask_file(const std::filesystem::path& path) {
    std::ifstream file(path);
    if (!file) {
        throw IoError(fmt::format("psrio: cannot read {}", path.string()));
    }
    std::vector<std::uint64_t> chans;
    std::string token;
    while (file >> token) {
        if (token.starts_with('#')) {
            std::getline(file, token);
            continue;
        }
        try {
            chans.push_back(std::stoull(token));
        } catch (const std::exception&) {
            throw ValidationError(
                fmt::format("psrio: {}: '{}' is not a channel index",
                            path.string(), token));
        }
    }
    return chans;
}

void write_stats_csv(const std::filesystem::path& path,
                     const stats::ChannelStats& result,
                     std::span<const std::uint8_t> mask,
                     const std::vector<double>& freqs) {
    std::ofstream file(path);
    if (!file) {
        throw IoError(fmt::format("psrio: cannot write {}", path.string()));
    }
    const auto mean   = result.mean();
    const auto stddev = result.stddev();
    const auto skew   = result.skew();
    const auto kurt   = result.kurtosis();
    file << "chan,freq_mhz,mean,std,skew,kurtosis,min,max,masked\n";
    for (std::size_t col = 0; col < mean.size(); ++col) {
        file << fmt::format("{},{},{},{},{},{},{},{},{}\n", col,
                            freqs[col % freqs.size()], mean[col], stddev[col],
                            skew[col], kurt[col], result.minima()[col],
                            result.maxima()[col], mask[col]);
    }
}

std::size_t count_masked(std::span<const std::uint8_t> mask) {
    std::size_t count = 0;
    for (const auto flag : mask) {
        count += flag != 0U ? 1U : 0U;
    }
    return count;
}

void run_rfi(const CleanOptions& opts, const Context& ctx) {
    const auto output = output_path(opts, "_masked");
    require_distinct(opts.inputs, output);
    auto input  = open_filterbank(opts, ctx);
    auto& block = input.block();
    const auto result =
        stats::compute(block, stats::ComputeOptions{
                                  .start  = 0,
                                  .nsamps = opts.stats_samples,
                                  .gulp   = opts.gulp,
                                  .order  = stats::Moments::kFull,
                              });
    const auto method =
        opts.method == "iqrm" ? clean::Method::kIqrm : clean::Method::kMad;
    auto mask         = clean::channel_mask(result, method, opts.threshold);
    const auto fill   = clean::fill_value(result, mask);
    const auto masked = count_masked(mask);
    if (!opts.mask_out.empty()) {
        write_mask(opts.mask_out, mask);
    }
    if (!opts.stats_out.empty()) {
        write_stats_csv(opts.stats_out, result, mask,
                        channel_frequencies(input.header()));
    }
    BlockSource cleaned(
        clean::Source(block, clean::Steps{
                                 .mask            = std::move(mask),
                                 .fill            = static_cast<float>(fill),
                                 .zerodm_bandpass = std::nullopt,
                             }));
    const auto written = write_stream(cleaned, input.header(), output,
                                      opts.gulp, ctx.drop_cache);
    ctx.info("psrio clean rfi: masked {} of {} channels ({:.1f}%) from {} "
             "samples of statistics; wrote {} samples to {}",
             masked, result.ncols(),
             100.0 * static_cast<double>(masked) /
                 static_cast<double>(result.ncols()),
             result.count(), written, output.string());
    ctx.out << output.string() << '\n';
}

void run_zap_channels(const CleanOptions& opts, const Context& ctx) {
    const auto output = output_path(opts, "_masked");
    require_distinct(opts.inputs, output);
    auto chans = opts.chans;
    if (!opts.mask_file.empty()) {
        const auto listed = read_mask_file(opts.mask_file);
        chans.insert(chans.end(), listed.begin(), listed.end());
    }
    if (chans.empty()) {
        throw ValidationError("psrio: give -c/--chans or --mask-file");
    }
    auto input     = open_filterbank(opts, ctx);
    auto& block    = input.block();
    const auto row = input.row();
    std::vector<std::uint8_t> mask(static_cast<std::size_t>(row), 0);
    for (const auto chan : chans) {
        if (chan >= row) {
            throw ValidationError(fmt::format(
                "psrio: channel {} is outside the {} channels", chan, row));
        }
        mask[static_cast<std::size_t>(chan)] = 1;
    }
    const auto masked = count_masked(mask);
    BlockSource cleaned(
        clean::Source(block, clean::Steps{
                                 .mask = std::move(mask),
                                 .fill = static_cast<float>(opts.fill),
                                 .zerodm_bandpass = std::nullopt,
                             }));
    const auto written = write_stream(cleaned, input.header(), output,
                                      opts.gulp, ctx.drop_cache);
    ctx.info("psrio clean zap-channels: zeroed {} channels; wrote {} samples "
             "to {}",
             masked, written, output.string());
    ctx.out << output.string() << '\n';
}

void run_zerodm(const CleanOptions& opts, const Context& ctx) {
    const auto output = output_path(opts, "_noZeroDM");
    require_distinct(opts.inputs, output);
    auto input  = open_filterbank(opts, ctx);
    auto& block = input.block();
    const auto result =
        stats::compute(block, stats::ComputeOptions{
                                  .start  = 0,
                                  .nsamps = std::nullopt,
                                  .gulp   = opts.gulp,
                                  .order  = stats::Moments::kMean,
                              });
    const auto mean = result.mean();
    BlockSource cleaned(clean::Source(
        block,
        clean::Steps{
            .mask            = {},
            .fill            = 0.0F,
            .zerodm_bandpass = std::vector<double>(mean.begin(), mean.end()),
        }));
    const auto written = write_stream(cleaned, input.header(), output,
                                      opts.gulp, ctx.drop_cache);
    ctx.info("psrio clean zerodm: wrote {} samples to {}", written,
             output.string());
    ctx.out << output.string() << '\n';
}

void run_zap_times(const CleanOptions& opts, const Context& ctx) {
    if (opts.zap_end <= opts.zap_start) {
        throw ValidationError(fmt::format("psrio: zap range [{}, {}) is empty",
                                          opts.zap_start, opts.zap_end));
    }
    clean::ZapMode mode = clean::ZapMode::kZero;
    if (opts.zap_mode == "gaussian") {
        mode = clean::ZapMode::kGaussian;
    } else if (opts.zap_mode == "samples") {
        mode = clean::ZapMode::kSamples;
    }
    const auto output = output_path(opts, "_zap");
    require_distinct(opts.inputs, output);
    std::mt19937_64 rng{std::random_device{}()};
    auto input         = open_filterbank(opts, ctx);
    const auto row     = input.row();
    const auto written = write_stream(
        input.block(), input.header(), output, opts.gulp, ctx.drop_cache,
        [&](std::uint64_t first, std::uint64_t rows, std::span<float> values) {
            const auto lo = std::max(opts.zap_start, first);
            const auto hi = std::min(opts.zap_end, first + rows);
            if (lo < hi) {
                clean::zap_rows(values, rows, row, lo - first, hi - first, mode,
                                rng);
            }
        });
    ctx.info("psrio clean zap-times: replaced samples [{}, {}) with {}; "
             "wrote {} samples to {}",
             opts.zap_start, opts.zap_end, opts.zap_mode, written,
             output.string());
    ctx.out << output.string() << '\n';
}

CLI::App* add_sub(CLI::App& clean,
                  const char* name,
                  const char* help,
                  CleanOptions& opts,
                  const char* default_name) {
    auto* sub = clean.add_subcommand(name, help);
    sub->add_option("inputs", opts.inputs,
                    "Input files in time order: .fil or .sf/.fits, or one .h5")
        ->required()
        ->expected(1, -1)
        ->check(CLI::ExistingFile);
    sub->add_option("-o,--outfile", opts.output,
                    fmt::format("Output filterbank (default: <input stem>{}"
                                ".fil in the working directory)",
                                default_name));
    sub->add_option("-g,--gulp", opts.gulp, "Samples per block")
        ->capture_default_str()
        ->check(CLI::PositiveNumber);
    add_threads_option(*sub);
    return sub;
}

} // namespace

void add_clean(CLI::App& app, Context& ctx) {
    auto* clean = app.add_subcommand(
        "clean", "Mask RFI channels, zap samples, or remove zero-DM "
                 "(spp_clean)");
    clean->require_subcommand(1);

    const auto rfi = std::make_shared<CleanOptions>();
    auto* sub = add_sub(*clean, "rfi",
                        "Mask channels whose variance, skewness or kurtosis "
                        "is an outlier, filling them with the median "
                        "channel mean",
                        *rfi, "_masked");
    sub->add_option("-m,--method", rfi->method, "mad or iqrm")
        ->capture_default_str()
        ->check(CLI::IsMember({"mad", "iqrm"}));
    sub->add_option("-t,--threshold", rfi->threshold,
                    "Sigma threshold for outliers")
        ->capture_default_str()
        ->check(CLI::PositiveNumber);
    sub->add_option("-s,--mask-out", rfi->mask_out,
                    "Write the masked channel indices, one per line");
    sub->add_option("--stats-out", rfi->stats_out,
                    "Write the channel statistics and mask as CSV");
    sub->add_option("--stats-samples", rfi->stats_samples,
                    "Estimate statistics from the first N samples only "
                    "(default: all)")
        ->check(CLI::PositiveNumber);
    sub->callback([rfi, &ctx] { run_rfi(*rfi, ctx); });

    const auto zap = std::make_shared<CleanOptions>();
    sub = add_sub(*clean, "zap-channels",
                  "Set chosen channels to a value (zero by default)", *zap,
                  "_masked");
    sub->add_option("-c,--chans", zap->chans, "Channels to zap")
        ->delimiter(',');
    sub->add_option("--mask-file", zap->mask_file,
                    "Text file of channel indices, one per line")
        ->check(CLI::ExistingFile);
    sub->add_option("--fill", zap->fill, "Value for zapped channels")
        ->capture_default_str();
    sub->callback([zap, &ctx] { run_zap_channels(*zap, ctx); });

    const auto times = std::make_shared<CleanOptions>();
    sub = add_sub(*clean, "zap-times", "Replace a range of time samples",
                  *times, "_zap");
    sub->add_option("-s,--start", times->zap_start, "First sample to replace")
        ->required();
    sub->add_option("-e,--end", times->zap_end, "End sample (exclusive)")
        ->required();
    sub->add_option("-m,--mode", times->zap_mode,
                    "zero, gaussian (each channel's mean and std over the "
                    "rest of the block), or samples (random other samples)")
        ->capture_default_str()
        ->check(CLI::IsMember({"zero", "gaussian", "samples"}));
    sub->callback([times, &ctx] { run_zap_times(*times, ctx); });

    const auto zerodm = std::make_shared<CleanOptions>();
    sub               = add_sub(*clean, "zerodm",
                                "Subtract the bandpass-weighted zero-DM series "
                                "(Eatough et al. 2009, as PRESTO)",
                                *zerodm, "_noZeroDM");
    sub->callback([zerodm, &ctx] { run_zerodm(*zerodm, ctx); });
}

} // namespace psrio::cli
