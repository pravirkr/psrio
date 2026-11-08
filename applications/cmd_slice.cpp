#include "common.hpp"

#include <CLI/CLI.hpp>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fmt/format.h>
#include <memory>
#include <optional>
#include <psrio/detail/exceptions.hpp>
#include <span>
#include <vector>

namespace psrio::cli {

namespace {

struct SliceOptions {
    std::filesystem::path input;
    std::filesystem::path output;
    std::optional<std::uint64_t> start;
    std::optional<double> start_time;
    std::optional<std::uint64_t> nsamps;
    std::optional<double> duration;
    std::optional<std::uint64_t> start_chan;
    std::optional<std::uint64_t> nchans;
    std::uint64_t gulp{kDefaultGulp};
};

/// Copy elements [first, first + width) of each @p row_len row of @p in.
template <typename T>
std::span<const T> take_columns(std::span<const T> in,
                                std::vector<T>& out,
                                std::size_t row_len,
                                std::size_t first,
                                std::size_t width) {
    if (width == row_len) {
        return in;
    }
    const auto nrows = in.size() / row_len;
    out.resize(nrows * width);
    for (std::size_t row = 0; row < nrows; ++row) {
        const auto src = in.subspan((row * row_len) + first, width);
        std::ranges::copy(src, out.begin() +
                                   static_cast<std::ptrdiff_t>(row * width));
    }
    return out;
}

void run_slice(const SliceOptions& opts, const Context& ctx) {
    require_distinct(opts.input, opts.output);
    Input input(opts.input);
    const auto range =
        resolve_range(ctx, input.nsamples(), input.tsamp(), opts.start,
                      opts.start_time, opts.nsamps, opts.duration);

    const auto nchans = input.nchans();
    const auto c0     = opts.start_chan.value_or(0U);
    if (c0 >= nchans) {
        throw ValidationError(fmt::format(
            "psrio: start channel {} is outside the {} input channels", c0,
            nchans));
    }
    const auto nc = opts.nchans.value_or(nchans - c0);
    if (nc == 0U || nc > nchans - c0) {
        throw ValidationError(fmt::format(
            "psrio: {} channels from channel {} exceed the {} input channels",
            nc, c0, nchans));
    }

    auto header = input.header();
    if (nc != nchans) {
        const auto freqs = channel_frequencies(header);
        set_channel_frequencies(header, std::span<const double>(freqs).subspan(
                                            static_cast<std::size_t>(c0),
                                            static_cast<std::size_t>(nc)));
    }
    if (header.tstart) {
        header.tstart = input.header().mjd_after_samples(range.start);
    }
    const auto nifs = input.nifs();
    const auto bits = static_cast<std::uint64_t>(input.nbits());
    if (is_series_path(opts.output)) {
        header.nbits = 32;
        header.signed_data.reset();
    } else if ((nifs * nc * bits) % 8U != 0U) {
        throw ValidationError(fmt::format(
            "psrio: {} channels of {} bits do not fill whole bytes", nc, bits));
    }
    Output output(opts.output, header);
    input.seek(range.start);

    const auto gulp          = std::max<std::uint64_t>(opts.gulp, 1U);
    const bool bytes_aligned = (c0 * bits) % 8U == 0U &&
                               (nc * bits) % 8U == 0U &&
                               (nchans * bits) % 8U == 0U;
    const bool lossless      = !input.is_series() && !output.is_series();

    std::uint64_t done = 0;
    if (lossless && bytes_aligned) {
        auto& source         = input.block();
        const auto row_bytes = static_cast<std::size_t>(nchans * bits / 8U);
        std::vector<std::byte> block(
            static_cast<std::size_t>(gulp * source.bytes_per_sample()));
        std::vector<std::byte> picked;
        while (done < range.count) {
            const auto want = std::min(gulp, range.count - done);
            const auto got  = source.read_block(want, block);
            if (got == 0U) {
                break;
            }
            const auto data = take_columns<std::byte>(
                std::span<const std::byte>(block).first(
                    static_cast<std::size_t>(got * source.bytes_per_sample())),
                picked, row_bytes, static_cast<std::size_t>(c0 * bits / 8U),
                static_cast<std::size_t>(nc * bits / 8U));
            output.sigproc_writer().write_block(got, data);
            done += got;
        }
    } else if (lossless) {
        auto& source = input.block();
        std::vector<std::uint8_t> block(
            static_cast<std::size_t>(gulp * nifs * nchans));
        std::vector<std::uint8_t> picked;
        while (done < range.count) {
            const auto want = std::min(gulp, range.count - done);
            const auto got  = source.read_samples(
                want, std::span<std::uint8_t>(block).first(
                          static_cast<std::size_t>(want * nifs * nchans)));
            if (got == 0U) {
                break;
            }
            const auto data = take_columns<std::uint8_t>(
                std::span<const std::uint8_t>(block).first(
                    static_cast<std::size_t>(got * nifs * nchans)),
                picked, static_cast<std::size_t>(nchans),
                static_cast<std::size_t>(c0), static_cast<std::size_t>(nc));
            output.sigproc_writer().write_samples(got, data);
            done += got;
        }
    } else {
        std::vector<float> block(
            static_cast<std::size_t>(gulp * nifs * nchans));
        std::vector<float> picked;
        while (done < range.count) {
            const auto want = std::min(gulp, range.count - done);
            const auto got  = input.read_floats(want, block);
            if (got == 0U) {
                break;
            }
            const auto data = take_columns<float>(
                std::span<const float>(block).first(
                    static_cast<std::size_t>(got * nifs * nchans)),
                picked, static_cast<std::size_t>(nchans),
                static_cast<std::size_t>(c0), static_cast<std::size_t>(nc));
            output.write_floats(got, data);
            done += got;
        }
    }
    output.close();
    ctx.info("psrio slice: wrote {} samples x {} channels from sample {} to {}",
             output.nsamples_written(), nc, range.start,
             output.path().string());
}

} // namespace

void add_slice(CLI::App& app, Context& ctx) {
    const auto opts = std::make_shared<SliceOptions>();
    auto* sub       = app.add_subcommand(
        "slice", "Copy a time and channel range (SIGPROC 'extract')");
    sub->alias("extract");
    sub->add_option("input", opts->input,
                    "Input .fil, .tim, .dat, .sf/.fits, or .h5")
        ->required()
        ->check(CLI::ExistingFile);
    sub->add_option("-o,--output", opts->output,
                    "Output file; .dat writes PRESTO, .tim a time series")
        ->required();
    auto* start = sub->add_option("-s,--start-sample", opts->start,
                                  "First time sample (default 0)");
    sub->add_option("--start-time", opts->start_time, "Start offset in seconds")
        ->excludes(start);
    auto* nsamps =
        sub->add_option("-n,--nsamps", opts->nsamps,
                        "Time samples to copy (default: to the end)");
    sub->add_option("--duration", opts->duration, "Seconds to copy")
        ->excludes(nsamps);
    sub->add_option("--start-chan", opts->start_chan,
                    "First channel, in file order (default 0)");
    sub->add_option("--nchans", opts->nchans,
                    "Channels to copy (default: to the last)")
        ->check(CLI::PositiveNumber);
    sub->add_option("--gulp", opts->gulp, "Time samples per block")
        ->capture_default_str()
        ->check(CLI::PositiveNumber);
    sub->callback([opts, &ctx] { run_slice(*opts, ctx); });
}

} // namespace psrio::cli
