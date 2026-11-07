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
#include <psrio/formats/sigproc/header.hpp>
#include <psrio/scrunch.hpp>
#include <span>
#include <vector>

namespace psrio::cli {

namespace {

struct DecimateOptions {
    std::filesystem::path input;
    std::filesystem::path output;
    std::uint64_t time_factor{1};
    std::uint64_t chan_factor{1};
    std::optional<int> nbits;
    bool sum{false};
    std::uint64_t gulp{kDefaultGulp};
};

sigproc::FilterbankHeader decimated_header(const sigproc::FilterbankHeader& in,
                                           std::uint64_t tfac,
                                           std::uint64_t cfac,
                                           int nbits) {
    auto out = in;
    if (out.tsamp) {
        *out.tsamp *= static_cast<double>(tfac);
    }
    if (cfac > 1U) {
        if (!in.frequency_table.empty()) {
            std::vector<double> freqs;
            for (std::size_t i = 0; i < in.frequency_table.size(); i += cfac) {
                double acc = 0.0;
                for (std::size_t k = 0; k < cfac; ++k) {
                    acc += in.frequency_table[i + k];
                }
                freqs.push_back(acc / static_cast<double>(cfac));
            }
            set_channel_frequencies(out, freqs);
        } else if (in.fch1 && in.foff) {
            out.fch1 =
                *in.fch1 + (0.5 * *in.foff * static_cast<double>(cfac - 1U));
            out.foff = *in.foff * static_cast<double>(cfac);
        }
        out.nchans = static_cast<std::int32_t>(
            static_cast<std::uint64_t>(in.nchans) / cfac);
    }
    out.nbits = nbits;
    if (out.signed_data && (nbits != 8 || in.nbits != 8)) {
        out.signed_data = 0;
    }
    return out;
}

void run_decimate(const DecimateOptions& opts, const Context& ctx) {
    require_distinct(opts.input, opts.output);
    Input input(opts.input);
    const auto nchans = input.nchans();
    const auto nifs   = input.nifs();
    const auto tfac   = opts.time_factor;
    const auto cfac   = opts.chan_factor;
    if (nchans % cfac != 0U) {
        throw ValidationError(fmt::format(
            "psrio: channel factor {} does not divide the {} channels", cfac,
            nchans));
    }
    const bool identity = tfac == 1U && cfac == 1U;
    const int nbits     = opts.nbits.value_or(
        identity && !is_series_path(opts.output) ? input.nbits() : 32);

    Output output(opts.output,
                  decimated_header(input.header(), tfac, cfac, nbits));

    const auto gulp       = std::max<std::uint64_t>(opts.gulp, 1U);
    const auto nchans_out = nchans / cfac;
    TimeScruncher scruncher(nifs * nchans_out, tfac);
    const float scale =
        opts.sum ? 1.0F : 1.0F / static_cast<float>(tfac * cfac);
    std::vector<float> block(static_cast<std::size_t>(gulp * nifs * nchans));
    while (true) {
        const auto got = input.read_floats(gulp, block);
        if (got == 0U) {
            break;
        }
        const auto chans = fscrunch_inplace(block, got * nifs, nchans, cfac);
        const auto rows  = scruncher.push(chans, got);
        if (rows.empty()) {
            continue;
        }
        if (scale != 1.0F) {
            std::ranges::transform(rows, rows.begin(),
                                   [scale](float v) { return v * scale; });
        }
        output.write_floats(rows.size() / (nifs * nchans_out), rows);
    }
    if (scruncher.pending() > 0U) {
        ctx.warn("dropped the last {} samples, fewer than the time factor {}",
                 scruncher.pending(), tfac);
    }
    output.close();
    ctx.info("psrio decimate: wrote {} samples x {} channels at {} bits to {}",
             output.nsamples_written(), nchans_out, nbits,
             output.path().string());
}

} // namespace

void add_decimate(CLI::App& app, Context& ctx) {
    const auto opts = std::make_shared<DecimateOptions>();
    auto* sub       = app.add_subcommand(
        "decimate",
        "Average adjacent samples and channels (SIGPROC 'decimate')");
    sub->add_option("input", opts->input, "Input .fil, .tim, .dat, or .h5")
        ->required()
        ->check(CLI::ExistingFile);
    sub->add_option("-o,--output", opts->output,
                    "Output file; .tim or .dat when one channel remains")
        ->required();
    sub->add_option("-t,--time-factor", opts->time_factor,
                    "Time samples combined into one")
        ->capture_default_str()
        ->check(CLI::PositiveNumber);
    sub->add_option("-c,--chan-factor", opts->chan_factor,
                    "Channels combined into one; must divide nchans")
        ->capture_default_str()
        ->check(CLI::PositiveNumber);
    sub->add_option("-b,--nbits", opts->nbits,
                    "Output bits: 8, 16, or 32 (default: input bits without "
                    "decimation, else 32)")
        ->check(CLI::IsMember({8, 16, 32}));
    sub->add_flag("--sum", opts->sum, "Sum instead of averaging");
    sub->add_option("--gulp", opts->gulp, "Time samples per block")
        ->capture_default_str()
        ->check(CLI::PositiveNumber);
    sub->callback([opts, &ctx] { run_decimate(*opts, ctx); });
}

} // namespace psrio::cli
