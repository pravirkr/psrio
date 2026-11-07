#include "common.hpp"

#include <CLI/CLI.hpp>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <psrio/digifil.hpp>

namespace psrio::cli {

namespace {

struct DigifilOptions {
    std::filesystem::path input;
    std::filesystem::path output;
    digifil::Config config;
    std::optional<std::uint64_t> gulp;
};

void run_digifil(const DigifilOptions& opts, const Context& ctx) {
    require_distinct(opts.input, opts.output);
    Input input(opts.input);
    auto config = opts.config;
    config.gulp = opts.gulp;
    const auto summary =
        digifil::run(input.block(), input.header(), opts.output, config);
    ctx.debug("psrio digifil: {} input samples per block", summary.gulp);
    ctx.info("psrio digifil: wrote {} samples at {} bits from {} input samples "
             "to {}",
             summary.nsamples_out, config.nbits, summary.nsamples_in,
             opts.output.string());
}

} // namespace

void add_digifil(CLI::App& app, Context& ctx) {
    const auto opts = std::make_shared<DigifilOptions>();
    auto& config    = opts->config;
    auto* sub       = app.add_subcommand(
        "digifil", "Rescale and redigitise a filterbank, as dspsr's digifil");
    sub->add_option("input", opts->input, "Input .fil or .h5")
        ->required()
        ->check(CLI::ExistingFile);
    sub->add_option("-o,--output", opts->output, "Output filterbank")
        ->required();
    sub->add_option("-b,--nbits", config.nbits,
                    "Output bits: 1, 2, 4, 8, 16, or -32 for float")
        ->capture_default_str()
        ->check(CLI::IsMember({1, 2, 4, 8, 16, -32}));
    sub->add_option("-t,--tscrunch", config.tscrunch,
                    "Time samples summed into one")
        ->capture_default_str();
    sub->add_option("-f,--fscrunch", config.fscrunch,
                    "Channels summed into one; must divide nchans")
        ->capture_default_str();
    sub->add_option("-I,--rescale-seconds", config.rescale_seconds,
                    "Rescale interval in seconds; 0 disables rescaling")
        ->capture_default_str();
    sub->add_flag("-c,--rescale-constant", config.rescale_constant,
                  "Keep the offset and scale from the first interval");
    sub->add_option("-s,--scale", config.scale_fac,
                    "Extra scale applied before digitising")
        ->capture_default_str();
    sub->add_option("-P,--poln-select", config.poln_select,
                    "Keep only this IF (-1 keeps all)")
        ->capture_default_str();
    sub->add_option("--gulp", opts->gulp,
                    "Input samples per block (default: dspsr's block for "
                    "256 MB)")
        ->check(CLI::PositiveNumber);
    sub->callback([opts, &ctx] { run_digifil(*opts, ctx); });
}

} // namespace psrio::cli
