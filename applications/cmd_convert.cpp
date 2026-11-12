#include "common.hpp"

#include <CLI/CLI.hpp>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <psrio/detail/exceptions.hpp>
#include <psrio/digifil.hpp>
#include <span>
#include <vector>

namespace psrio::cli {

namespace {

struct ConvertOptions {
    std::vector<std::filesystem::path> inputs;
    std::filesystem::path output;
    std::optional<int> nbits;
    bool fswap{false};
    std::uint64_t gulp{kDefaultGulp};
};

void copy_bytes(Input& input, Output& output, std::uint64_t gulp) {
    auto& source = input.block();
    std::vector<std::byte> block(
        static_cast<std::size_t>(gulp * source.bytes_per_sample()));
    while (const auto got = source.read_block(gulp, block)) {
        output.sigproc_writer().write_block(got, block);
    }
}

void copy_floats(Input& input, Output& output, std::uint64_t gulp) {
    std::vector<float> block(static_cast<std::size_t>(gulp * input.row()));
    while (const auto got = input.read_floats(gulp, block)) {
        output.write_floats(got,
                            std::span<const float>(block).first(
                                static_cast<std::size_t>(got * input.row())));
    }
}

void run_convert(const ConvertOptions& opts, const Context& ctx) {
    require_distinct(opts.inputs, opts.output);
    Input input(opts.inputs);
    const bool to_series = is_series_path(opts.output);
    const int in_bits    = input.nbits();
    const int bits       = opts.nbits.value_or(to_series ? 32 : in_bits);
    if (to_series && bits != 32) {
        throw ValidationError("psrio: .tim and .dat series hold 32-bit floats");
    }

    auto header = input.header();
    if (opts.fswap) {
        if (input.nchans() < 2U) {
            throw ValidationError("psrio: --fswap needs more than one channel");
        }
        auto freqs = channel_frequencies(header);
        std::ranges::reverse(freqs);
        set_channel_frequencies(header, freqs);
        input.block().set_fswap(true);
    }

    if (bits < in_bits) {
        if (opts.fswap) {
            throw ValidationError(
                "psrio: --fswap cannot be combined with a lower bit depth; "
                "digifil always writes descending frequencies");
        }
        digifil::Config config;
        config.nbits = bits;
        const auto summary =
            digifil::run(input.block(), input.header(), opts.output, config);
        ctx.info("psrio convert: rescaled to {} bits with digifil's defaults "
                 "(see 'psrio digifil'); wrote {} samples to {}",
                 bits, summary.nsamples_out, opts.output.string());
        return;
    }

    const auto gulp = std::max<std::uint64_t>(opts.gulp, 1U);
    if (bits == in_bits && !input.is_series() && !to_series) {
        Output output(opts.output, header);
        copy_bytes(input, output, gulp);
        output.close();
        ctx.info("psrio convert: wrote {} samples to {}",
                 output.nsamples_written(), opts.output.string());
        return;
    }

    if (header.samples_are_signed() && bits == 16) {
        throw ValidationError(
            "psrio: signed 8-bit samples do not fit unsigned 16-bit samples; "
            "use -b 32");
    }
    header.nbits = bits;
    if (header.signed_data && bits != 8) {
        header.signed_data = 0;
    }
    Output output(opts.output, header);
    copy_floats(input, output, gulp);
    output.close();
    ctx.info("psrio convert: wrote {} samples at {} bits to {}",
             output.nsamples_written(), bits, output.path().string());
}

} // namespace

void add_convert(CLI::App& app, Context& ctx) {
    const auto opts = std::make_shared<ConvertOptions>();
    auto* sub       = app.add_subcommand(
        "convert",
        "Change bit depth or channel order, or convert .tim and .dat series");
    sub->add_option("inputs", opts->inputs,
                    "Input files in time order: .fil or .sf/.fits, or one "
                    ".tim, .dat, or .h5")
        ->required()
        ->expected(1, -1)
        ->check(CLI::ExistingFile);
    sub->add_option("-o,--output", opts->output,
                    "Output file; .dat writes PRESTO")
        ->required();
    sub->add_option("-b,--nbits", opts->nbits,
                    "Output bits (default: input bits). Fewer bits rescale "
                    "through digifil")
        ->check(CLI::IsMember({1, 2, 4, 8, 16, 32}));
    sub->add_flag("--fswap", opts->fswap, "Reverse the channel order");
    sub->add_option("--gulp", opts->gulp, "Time samples per block")
        ->capture_default_str()
        ->check(CLI::PositiveNumber);
    sub->callback([opts, &ctx] { run_convert(*opts, ctx); });
}

} // namespace psrio::cli
